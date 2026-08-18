function output=CoreFHorzQHExpAsset_nod1_z_noe_noa1_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Quasi-hyperbolic experienceasset, noa1: no d1, z, no e.
% The experience asset a2 is the only endogenous state. n_a is scalar (n_a_justexpasset);
% a_grid is the a2_grid. n_a_big/a_grid_big are unused. n_d=n_d2; d_grid=d2_grid
%
% Methods: BASE only -- with no a1 there is nothing for divide-and-conquer or the grid
% interpolation layer to operate on (same reason the baseline noa1 subcodes have no
% DC/GI/DC+GI blocks). Runs Naive then Sophisticated, each over lowmemory {0,1} and
% with a ValueFnFromPolicy oracle, then the exponential-discounting cross-tests.
% shocks: {z (markov)} -> valid lowmemory {0,1}.
%
% TEST-FIRST: the toolkit currently has NO quasi-hyperbolic support for experienceasset at
% all (no QH+ExpAsset raws, no dispatcher branch), so this errors at the first ValueFnIter
% call. That is expected: this test is written ahead of the toolkit code.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();

ReturnFn=@(d2,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_noe_noa1_nosemiz(d2,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

% Experience asset
vfoptions.experienceasset=1;
simoptions.experienceasset=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

%% Quasi-Hyperbolic Discounting
vfoptions.exoticpreferences='QuasiHyperbolic';
vfoptions.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;

%% Naive
vfoptions.quasi_hyperbolic='Naive';
vfoptions1=vfoptions;
[V1,Policy1,V1alt,Policy1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

vfoptions1.lowmemory=1;
[V1B,Policy1B,V1Balt,Policy1Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('Naive lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('Naive lowmemory=1 (Valt), this should be zero: %2.8f \n',max(abs(V1alt(:)-V1Balt(:))))
fprintf('Naive lowmemory=1 (Policy), this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

vfoptions1.Policyalt=Policy1alt; % Naive QH: ValueFnFromPolicy reconstructs V from the exponential-discounter argmax (the 4th output of the Naive solve)
[V1fromPolicy,V1altfromPolicy]=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('Naive ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))
fprintf('Naive ValueFnFromPolicy (Valt), this should be zero: %2.8f \n',max(abs(V1altfromPolicy(:)-V1alt(:))))

%% Sophisticated
vfoptions.quasi_hyperbolic='Sophisticated';
vfoptions1=vfoptions;
[V1s,Policy1s,V1salt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

vfoptions1.lowmemory=1;
[V1sB,Policy1sB,V1sBalt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('Sophisticated lowmemory=1, this should be zero: %2.8f \n',max(abs(V1s(:)-V1sB(:))))
fprintf('Sophisticated lowmemory=1 (Valt), this should be zero: %2.8f \n',max(abs(V1salt(:)-V1sBalt(:))))
fprintf('Sophisticated lowmemory=1 (Policy), this should be zero: %2.8f \n',max(abs(Policy1s(:)-Policy1sB(:))))
vfoptions1.lowmemory=0;

vfoptions1.Policyalt=[]; % Sophisticated does not need Policyalt
[V1sfromPolicy,V1saltfromPolicy]=ValueFnFromPolicy_FHorz(Policy1s,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('Sophisticated ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1sfromPolicy(:)-V1s(:))))
fprintf('Sophisticated ValueFnFromPolicy (Valt), this should be zero: %2.8f \n',max(abs(V1saltfromPolicy(:)-V1salt(:))))

%% Versus exponential discounting
% (i) at the actual QH beta0, Naive's continuation value equals the exponential value function
vfoptionsE=vfoptions; vfoptionsE.exoticpreferences='None';
[Vexp,Policyexp]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsE);
vfoptionsN=vfoptions; vfoptionsN.quasi_hyperbolic='Naive';
[VnA,PolicynA,VnAalt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsN);
fprintf('(i) Naive continuation value == exponential (beta0=%g), should be zero: %2.8f \n',Params.beta0,max(abs(VnAalt(:)-Vexp(:))))

% (ii) with beta0=1, Naive main value AND continuation value both equal exponential
% (iii) with beta0=1, Sophisticated main value AND continuation value both equal exponential
beta0_store=Params.beta0; Params.beta0=1;
[Vexp1,~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsE);
[VnA1,~,VnA1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsN);
vfoptionsS=vfoptions; vfoptionsS.quasi_hyperbolic='Sophisticated';
[VsS1,~,VsS1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsS);
fprintf('(ii) Naive V == exponential (beta0=1), should be zero: %2.8f \n',max(abs(VnA1(:)-Vexp1(:))))
fprintf('(ii) Naive Valt == exponential (beta0=1), should be zero: %2.8f \n',max(abs(VnA1alt(:)-Vexp1(:))))
fprintf('(iii) Sophisticated V == exponential (beta0=1), should be zero: %2.8f \n',max(abs(VsS1(:)-Vexp1(:))))
fprintf('(iii) Sophisticated Valt == exponential (beta0=1), should be zero: %2.8f \n',max(abs(VsS1alt(:)-Vexp1(:))))
Params.beta0=beta0_store;

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
