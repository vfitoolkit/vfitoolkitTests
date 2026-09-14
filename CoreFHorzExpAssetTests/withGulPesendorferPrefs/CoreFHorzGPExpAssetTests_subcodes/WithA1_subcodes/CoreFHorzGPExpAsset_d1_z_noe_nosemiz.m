function output=CoreFHorzGPExpAsset_d1_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Gul-Pesendorfer experienceasset, single standard endogenous asset a1: with d1 (labour), z, no e.
% a = [a1 (standard), a2 (experience asset)]. n_d=[n_d1,n_d2]; d_grid=[d1_grid; d2_grid]
%
% Gul-Pesendorfer: V_j = max_{choices}[u + v + beta*EV_{j+1}] - max_{choices} v, where v is
% the temptation function. GP only changes the value fn solve (Policy maximizes the tempted
% objective); everything downstream of Policy is standard code.
%
% Methods: base / DC1 / GI1 / DC1_GI1 (single standard asset a1 -> DC1, not DC2A), with a
% ValueFnFromPolicy oracle on every method. GI changes the solution slightly vs base, so GI1
% is not compared to base for equality (only its own ValueFnFromPolicy + DC1_GI1-vs-GI1 are
% checked).
% shocks: {z (markov)} -> valid lowmemory {0,1}.
%
% TEST-FIRST: the toolkit currently has NO Gul-Pesendorfer support for experienceasset at
% all (the GP dispatcher errors on it), so this errors at the first ValueFnIter call. That
% is expected: this test is written ahead of the toolkit code.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();

ReturnFn=@(d1,d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_z_noe_nosemiz(d1,d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

% Experience asset
vfoptions.experienceasset=1;
simoptions.experienceasset=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

%% Gul-Pesendorfer: temptation and self-control
vfoptions.exoticpreferences='GulPesendorfer';
vfoptions.temptationFn=@(d1,d2,a1prime,a1,a2,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_d1_z_noe_nosemiz(d1,d2,a1prime,a1,a2,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
% Consumption is tempting: v = lambdaGP*u_c(c) + shiftGP (lambdaGP and shiftGP sit in Params)

%% Methods: base / DC1 / GI1 / DC1_GI1

% Base
vfoptions1=vfoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% ValueFnFromPolicy oracle on the base method
vfoptionsVFP1=vfoptions1; vfoptionsVFP1.lowmemory=0;
[V1fromPolicy]=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP1);
fprintf('ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1fromPolicy(:)-V1(:))))

% lowmemory on base
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1 (base), this should be zero: %.3e \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1 (base, Policy), this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

% Divide-and-conquer -> DC1, should give the same answer as base
vfoptions2=vfoptions; vfoptions2.divideandconquer=1;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('DC1, this should be zero: %.3e \n',max(abs(V1(:)-V2(:))))
fprintf('DC1 (Policy), this should be zero: %.3e \n',max(abs(Policy1(:)-Policy2(:))))
% ValueFnFromPolicy oracle on DC1
vfoptionsVFP2=vfoptions2; vfoptionsVFP2.lowmemory=0;
[V2fromPolicy]=ValueFnFromPolicy_FHorz(Policy2,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP2);
fprintf('ValueFnFromPolicy (DC1), this should be zero: %.3e \n',max(abs(V2fromPolicy(:)-V2(:))))

% lowmemory on DC1
vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (DC1), this should be zero: %.3e \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (DC1, Policy), this should be zero: %.3e \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=0;

% Grid interpolation -> GI1 (GI changes the solution slightly vs base, so no direct base equality check)
vfoptions3=vfoptions; vfoptions3.gridinterplayer=1; vfoptions3.ngridinterp=5;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
% ValueFnFromPolicy oracle on GI1
vfoptionsVFP3=vfoptions3; vfoptionsVFP3.lowmemory=0;
[V3fromPolicy]=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP3);
fprintf('ValueFnFromPolicy (GI1), this should be zero: %.3e \n',max(abs(V3fromPolicy(:)-V3(:))))

% lowmemory on GI1
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (GI1), this should be zero: %.3e \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (GI1, Policy), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=0;

% DC + GI -> DC1_GI1, should match GI1
vfoptions4=vfoptions; vfoptions4.divideandconquer=1; vfoptions4.gridinterplayer=1; vfoptions4.ngridinterp=5;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('DC1_GI1 vs GI1, this should be zero: %.3e \n',max(abs(V3(:)-V4(:))))
fprintf('DC1_GI1 vs GI1 (Policy), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy4(:))))
% ValueFnFromPolicy oracle on DC1_GI1
vfoptionsVFP4=vfoptions4; vfoptionsVFP4.lowmemory=0;
[V4fromPolicy]=ValueFnFromPolicy_FHorz(Policy4,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP4);
fprintf('ValueFnFromPolicy (DC1_GI1), this should be zero: %.3e \n',max(abs(V4fromPolicy(:)-V4(:))))

% lowmemory on DC1_GI1
vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1 (DC1_GI1), this should be zero: %.3e \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1 (DC1_GI1, Policy), this should be zero: %.3e \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=0;


%% V_Jplus1: use V of period jstar as the terminal value function of a shorter model
% As in the other test banks: solve the model, then solve a shorter model that runs only periods
% 1,...,jstar-1 with Njs=jstar-1 and the age-dependent parameters trimmed to length Njs, and check
% we get the same answer for periods 1,...,jstar-1.
% Note (Gul-Pesendorfer specific): the continuation value that enters the Bellman equation is
% just V itself (the most-tempting term is a within-period object, subtracted after the max),
% so it is V that gets fed back in as vfoptions.V_Jplus1.
% Note: mewj is age-dependent, but is only used for the agent distribution, which is not computed
% here, so it is left alone.

%% V_Jplus1, without divide-and-conquer, without grid interpolation
jstar=round(3*N_j/4);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions1;
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
fprintf('V_Jplus1, lowmemory=1, this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, with divide-and-conquer
jstar=round(2*N_j/3);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions2;
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
Vbase=Vbase(:,:,:,1:Njs);
Policybase=Policybase(:,:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i (with DC)), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i (with DC)), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, with grid interpolation
jstar=round(N_j/2);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions3;
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
Vbase=Vbase(:,:,:,1:Njs);
Policybase=Policybase(:,:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i (with GI)), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i (with GI)), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, with divide-and-conquer and grid interpolation
jstar=N_j;
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions4;
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
Vbase=Vbase(:,:,:,1:Njs);
Policybase=Policybase(:,:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i (with DC+GI)), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i (with DC+GI)), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with DC+GI), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with DC+GI), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

clear Vbase Policybase Vshort Policyshort

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
