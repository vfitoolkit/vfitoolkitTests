function output=CoreFHorzQHExpAsset_nod1_noz_noe_nosemiz_with2A1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Quasi-hyperbolic experienceasset with TWO standard endogenous assets -> triggers the
% DC2A / GI2A / DC2A_GI2A code paths (used whenever length(n_a1)>1): a1 is divide-conquered,
% the second standard asset a1_2 (return r2) is folded, a2 is the experience asset.
% a = [a1 (standard), a1_2 (standard), a2 (experience asset)]. n_d=n_d2; d_grid=d2_grid
% Note: n_a=[n_a1,n_a1_2,n_a2], a_grid=[a1_grid;a1_2_grid;a2_grid] and Params.r2 are built by
% the calling test script (CoreFHorzExpAsset_setup), exactly as for the baseline counterpart.
%
% Methods: base / DC2A / GI2A / DC2A_GI2A, for Naive then Sophisticated, with a
% ValueFnFromPolicy oracle on every method, then the exponential-discounting cross-tests.
% GI changes the solution slightly vs base, so GI2A is not compared to base for equality
% (only its own ValueFnFromPolicy + DC2A_GI2A-vs-GI2A are checked).
% shocks: {} (none) -> valid lowmemory {0}, so there is no lowmemory sweep on any method
% (same as the baseline counterpart, which has no lowmemory blocks either).
%
% TEST-FIRST: the toolkit currently has NO quasi-hyperbolic support for experienceasset at
% all (no QH+ExpAsset raws, no dispatcher branch), so this errors at the first ValueFnIter
% call. That is expected: this test is written ahead of the toolkit code.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
n_z=0;
z_grid=[];
pi_z=[];

ReturnFn=@(d2,a1prime,a1_2prime,a1,a1_2,a2,r,r2,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz_with2A1(d2,a1prime,a1_2prime,a1,a1_2,a2,r,r2,w,kappa_j,sigma,agej,Jr,pension);

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

%% Naive then Sophisticated: base / DC2A / GI2A / DC2A_GI2A
%% Naive
vfoptions.quasi_hyperbolic='Naive'; qh='Naive';

% 0 shocks (no z, no e, no semiz) -> only lowmemory=0 is valid; no lowmemory>0
% check on any of the four methods.

% Base
vfoptions1=vfoptions;
[V1,Policy1,V1alt,Policy1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% ValueFnFromPolicy oracle on the base method
vfoptionsVFP1=vfoptions1; vfoptionsVFP1.lowmemory=0;
vfoptionsVFP1.Policyalt=Policy1alt; % Naive QH: ValueFnFromPolicy reconstructs V from the exponential-discounter argmax
[V1fromPolicy,V1altfromPolicy]=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP1);
fprintf('%s ValueFnFromPolicy, this should be zero: %2.8f \n',qh,max(abs(V1fromPolicy(:)-V1(:))))
fprintf('%s ValueFnFromPolicy (Valt), this should be zero: %2.8f \n',qh,max(abs(V1altfromPolicy(:)-V1alt(:))))

% Divide-and-conquer -> DC2A, should give the same answer as base
vfoptions2=vfoptions; vfoptions2.divideandconquer=1;
[V2,Policy2,V2alt,Policy2alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('%s DC2A, this should be zero: %2.8f \n',qh,max(abs(V1(:)-V2(:))))
fprintf('%s DC2A (Valt), this should be zero: %2.8f \n',qh,max(abs(V1alt(:)-V2alt(:))))
fprintf('%s DC2A (Policy), this should be zero: %2.8f \n',qh,max(abs(Policy1(:)-Policy2(:))))
fprintf('%s DC2A (Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy1alt(:)-Policy2alt(:))));
% ValueFnFromPolicy oracle on DC2A
vfoptionsVFP2=vfoptions2; vfoptionsVFP2.lowmemory=0;
vfoptionsVFP2.Policyalt=Policy2alt; % Naive QH: ValueFnFromPolicy reconstructs V from the exponential-discounter argmax
[V2fromPolicy,V2altfromPolicy]=ValueFnFromPolicy_FHorz(Policy2,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP2);
fprintf('%s ValueFnFromPolicy (DC2A), this should be zero: %2.8f \n',qh,max(abs(V2fromPolicy(:)-V2(:))))
fprintf('%s ValueFnFromPolicy (DC2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V2altfromPolicy(:)-V2alt(:))))

% Grid interpolation -> GI2A (GI changes the solution slightly vs base, so no direct base equality check)
vfoptions3=vfoptions; vfoptions3.gridinterplayer=1; vfoptions3.ngridinterp=5;
[V3,Policy3,V3alt,Policy3alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
% ValueFnFromPolicy oracle on GI2A
vfoptionsVFP3=vfoptions3; vfoptionsVFP3.lowmemory=0;
vfoptionsVFP3.Policyalt=Policy3alt; % Naive QH: ValueFnFromPolicy reconstructs V from the exponential-discounter argmax
[V3fromPolicy,V3altfromPolicy]=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP3);
fprintf('%s ValueFnFromPolicy (GI2A), this should be zero: %2.8f \n',qh,max(abs(V3fromPolicy(:)-V3(:))))
fprintf('%s ValueFnFromPolicy (GI2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V3altfromPolicy(:)-V3alt(:))))

% DC + GI -> DC2A_GI2A, should match GI2A
vfoptions4=vfoptions; vfoptions4.divideandconquer=1; vfoptions4.gridinterplayer=1; vfoptions4.ngridinterp=5;
[V4,Policy4,V4alt,Policy4alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('%s DC2A_GI2A vs GI2A, this should be zero: %2.8f \n',qh,max(abs(V3(:)-V4(:))))
fprintf('%s DC2A_GI2A vs GI2A (Valt), this should be zero: %2.8f \n',qh,max(abs(V3alt(:)-V4alt(:))))
fprintf('%s DC2A_GI2A vs GI2A (Policy), this should be zero: %2.8f \n',qh,max(abs(Policy3(:)-Policy4(:))))
fprintf('%s DC2A_GI2A vs GI2A (Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy3alt(:)-Policy4alt(:))));
% ValueFnFromPolicy oracle on DC2A_GI2A
vfoptionsVFP4=vfoptions4; vfoptionsVFP4.lowmemory=0;
vfoptionsVFP4.Policyalt=Policy4alt; % Naive QH: ValueFnFromPolicy reconstructs V from the exponential-discounter argmax
[V4fromPolicy,V4altfromPolicy]=ValueFnFromPolicy_FHorz(Policy4,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP4);
fprintf('%s ValueFnFromPolicy (DC2A_GI2A), this should be zero: %2.8f \n',qh,max(abs(V4fromPolicy(:)-V4(:))))
fprintf('%s ValueFnFromPolicy (DC2A_GI2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V4altfromPolicy(:)-V4alt(:))))

%% Sophisticated
vfoptions.quasi_hyperbolic='Sophisticated'; qh='Sophisticated';

% 0 shocks (no z, no e, no semiz) -> only lowmemory=0 is valid; no lowmemory>0
% check on any of the four methods.

% Base
vfoptions1=vfoptions;
[V1,Policy1,V1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% ValueFnFromPolicy oracle on the base method
vfoptionsVFP1=vfoptions1; vfoptionsVFP1.lowmemory=0;
[V1fromPolicy,V1altfromPolicy]=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP1);
fprintf('%s ValueFnFromPolicy, this should be zero: %2.8f \n',qh,max(abs(V1fromPolicy(:)-V1(:))))
fprintf('%s ValueFnFromPolicy (Valt), this should be zero: %2.8f \n',qh,max(abs(V1altfromPolicy(:)-V1alt(:))))

% Divide-and-conquer -> DC2A, should give the same answer as base
vfoptions2=vfoptions; vfoptions2.divideandconquer=1;
[V2,Policy2,V2alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('%s DC2A, this should be zero: %2.8f \n',qh,max(abs(V1(:)-V2(:))))
fprintf('%s DC2A (Valt), this should be zero: %2.8f \n',qh,max(abs(V1alt(:)-V2alt(:))))
fprintf('%s DC2A (Policy), this should be zero: %2.8f \n',qh,max(abs(Policy1(:)-Policy2(:))))
% ValueFnFromPolicy oracle on DC2A
vfoptionsVFP2=vfoptions2; vfoptionsVFP2.lowmemory=0;
[V2fromPolicy,V2altfromPolicy]=ValueFnFromPolicy_FHorz(Policy2,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP2);
fprintf('%s ValueFnFromPolicy (DC2A), this should be zero: %2.8f \n',qh,max(abs(V2fromPolicy(:)-V2(:))))
fprintf('%s ValueFnFromPolicy (DC2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V2altfromPolicy(:)-V2alt(:))))

% Grid interpolation -> GI2A (GI changes the solution slightly vs base, so no direct base equality check)
vfoptions3=vfoptions; vfoptions3.gridinterplayer=1; vfoptions3.ngridinterp=5;
[V3,Policy3,V3alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
% ValueFnFromPolicy oracle on GI2A
vfoptionsVFP3=vfoptions3; vfoptionsVFP3.lowmemory=0;
[V3fromPolicy,V3altfromPolicy]=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP3);
fprintf('%s ValueFnFromPolicy (GI2A), this should be zero: %2.8f \n',qh,max(abs(V3fromPolicy(:)-V3(:))))
fprintf('%s ValueFnFromPolicy (GI2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V3altfromPolicy(:)-V3alt(:))))

% DC + GI -> DC2A_GI2A, should match GI2A
vfoptions4=vfoptions; vfoptions4.divideandconquer=1; vfoptions4.gridinterplayer=1; vfoptions4.ngridinterp=5;
[V4,Policy4,V4alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('%s DC2A_GI2A vs GI2A, this should be zero: %2.8f \n',qh,max(abs(V3(:)-V4(:))))
fprintf('%s DC2A_GI2A vs GI2A (Valt), this should be zero: %2.8f \n',qh,max(abs(V3alt(:)-V4alt(:))))
fprintf('%s DC2A_GI2A vs GI2A (Policy), this should be zero: %2.8f \n',qh,max(abs(Policy3(:)-Policy4(:))))
% ValueFnFromPolicy oracle on DC2A_GI2A
vfoptionsVFP4=vfoptions4; vfoptionsVFP4.lowmemory=0;
[V4fromPolicy,V4altfromPolicy]=ValueFnFromPolicy_FHorz(Policy4,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP4);
fprintf('%s ValueFnFromPolicy (DC2A_GI2A), this should be zero: %2.8f \n',qh,max(abs(V4fromPolicy(:)-V4(:))))
fprintf('%s ValueFnFromPolicy (DC2A_GI2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V4altfromPolicy(:)-V4alt(:))))

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
