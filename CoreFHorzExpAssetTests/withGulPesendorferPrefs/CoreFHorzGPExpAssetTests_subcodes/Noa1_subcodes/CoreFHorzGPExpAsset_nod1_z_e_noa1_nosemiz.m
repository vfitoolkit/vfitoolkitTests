function output=CoreFHorzGPExpAsset_nod1_z_e_noa1_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Gul-Pesendorfer experienceasset, noa1: no d1, z, e.
% The experience asset a2 is the only endogenous state. n_a is scalar (n_a_justexpasset);
% a_grid is the a2_grid. n_a_big/a_grid_big are unused. n_d=n_d2; d_grid=d2_grid
%
% Methods: BASE only -- with no a1 there is nothing for divide-and-conquer or the grid
% interpolation layer to operate on (same reason the baseline noa1 subcodes have no
% DC/GI/DC+GI blocks). Runs the GP solve over lowmemory {0,1,2} and
% with a ValueFnFromPolicy oracle.
% shocks: {z (markov), e (iid)} -> valid lowmemory {0,1,2}.
%
% Gul-Pesendorfer (temptation and self-control): a temptation fn v alongside the return fn u, with
%   V_j = max_{d,a'} [ u + v + beta*E V_{j+1} ] - max_{d,a'} v
% so Policy maximizes the tempted objective and V nets off the self-control cost.
%
% TEST-FIRST: the toolkit currently has NO Gul-Pesendorfer support for experienceasset
% (the GP dispatcher errors on it), so this errors at the first ValueFnIter call. That is
% expected: this test is written ahead of the toolkit code.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% e
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.pi_e=simoptionsbaseline.pi_e;
simoptions.e_grid=simoptionsbaseline.e_grid;

ReturnFn=@(d2,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_e_noa1_nosemiz(d2,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);

% Experience asset
vfoptions.experienceasset=1;
simoptions.experienceasset=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

%% Gul-Pesendorfer: temptation and self-control
vfoptions.exoticpreferences='GulPesendorfer';
vfoptions.temptationFn=@(d2,a,z,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_nod1_z_e_noa1_nosemiz(d2,a,z,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
% Consumption is tempting: v = lambdaGP*u_c(c) + shiftGP (lambdaGP and shiftGP sit in Params)

%% Solve
vfoptions1=vfoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1 (Policy), this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=2;
[V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=2, this should be zero: %.3e \n',max(abs(V1(:)-V1C(:))))
fprintf('lowmemory=2 (Policy), this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=0;

[V1fromPolicy]=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1fromPolicy(:)-V1(:))))


%% V_Jplus1: use V of period jstar as the terminal value function of a shorter model
% As in the other test banks: solve the model, then solve a shorter model that runs only periods
% 1,...,jstar-1 with Njs=jstar-1 and the age-dependent parameters trimmed to length Njs, and check
% we get the same answer for periods 1,...,jstar-1.
% Gul-Pesendorfer: the continuation value that enters the Bellman equation is V itself (the
% self-control-netted value), so it is V that gets fed back in as vfoptions.V_Jplus1.
% Note: mewj is age-dependent, but is only used for the agent distribution, which is not computed
% here, so it is left alone.
% This tier has no a1 for divide-and-conquer or the grid interp layer to operate on, so it
% defines only vfoptions1. Run at jstar=round(3*N_j/4) and again at jstar=N_j, so the
% retirement periods and the terminal V_Jplus1 branch are covered too.
for jstar=[round(3*N_j/4),N_j]
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    vfoptionsjs=vfoptions1;
    vfoptionsjs.exoticpreferences='GulPesendorfer';
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
    Vbase=Vbase(:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1, this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 (Policy), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=2;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=2, this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=2 (Policy), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;
end

clear Vbase Policybase Vshort Policyshort

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
