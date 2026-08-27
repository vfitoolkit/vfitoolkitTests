function output=CoreFHorzQHExpAssete_nod1_noz_e_with2A1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Quasi-hyperbolic experienceassete + semiz with TWO standard endogenous assets -> triggers
% the DC2A / GI2A / DC2A_GI2A code paths (used whenever length(n_a1)>1): a1_1 is divide-
% conquered, the binary a1_2 is folded, a2 is the experience asset. no d1, no z, e.
% n_d=[n_d2,n_d3]; d_grid=[d2_grid; d3_grid]
%
% Methods: base / DC2A / GI2A / DC2A_GI2A, for Naive then Sophisticated, with a
% ValueFnFromPolicy oracle on every method, then the exponential-discounting cross-tests.
% GI changes the solution slightly vs base, so GI2A is not compared to base for equality
% (only its own ValueFnFromPolicy + DC2A_GI2A-vs-GI2A are checked).
% shocks: {semiz, e (iid)} -> valid lowmemory {0,1,2}.
%
% TEST-FIRST: the toolkit currently has NO quasi-hyperbolic support for experienceassete at
% all (no QH+ExpAssete raws, no dispatcher branch), so this errors at the first ValueFnIter
% call. That is expected: this test is written ahead of the toolkit code.

% n_a=[a1_1 (divide-conquered), a1_2 (multi-point, folded), a2 (experience asset)] and
% a_grid arrive already built from the calling test script; Params.r2 comes from the setup.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
n_z=0;
z_grid=[];
pi_z=[];
% semiz
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
% e
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.pi_e=simoptionsbaseline.pi_e;
simoptions.e_grid=simoptionsbaseline.e_grid;

ReturnFn=@(d2,d3,a1prime,a1_2prime,a1,a1_2,a2,semiz,e,r,r2,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssete_nod1_noz_e_with2A1_semiz(d2,d3,a1prime,a1_2prime,a1,a1_2,a2,semiz,e,r,r2,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% Experience asset (e variant)
vfoptions.experienceassete=1;
simoptions.experienceassete=1;
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

% Base
vfoptions1=vfoptions;
[V1,Policy1,V1alt,Policy1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% ValueFnFromPolicy oracle on the base method
vfoptionsVFP1=vfoptions1; vfoptionsVFP1.lowmemory=0;
vfoptionsVFP1.Policyalt=Policy1alt; % Naive QH: ValueFnFromPolicy reconstructs V from the exponential-discounter argmax
[V1fromPolicy,V1altfromPolicy]=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP1);
fprintf('%s ValueFnFromPolicy, this should be zero: %.3e \n',qh,max(abs(V1fromPolicy(:)-V1(:))))
fprintf('%s ValueFnFromPolicy (Valt), this should be zero: %.3e \n',qh,max(abs(V1altfromPolicy(:)-V1alt(:))))

% lowmemory on base
vfoptions1.lowmemory=1;
[V1B,Policy1B,V1Balt,Policy1Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('%s lowmemory=1 (base), this should be zero: %.3e \n',qh,max(abs(V1(:)-V1B(:))))
fprintf('%s lowmemory=1 (base, Valt), this should be zero: %.3e \n',qh,max(abs(V1alt(:)-V1Balt(:))))
fprintf('%s lowmemory=1 (base, Policy), this should be zero: %.3e \n',qh,max(abs(Policy1(:)-Policy1B(:))))
fprintf('%s lowmemory=1 (base, Policyalt), this should be zero: %.3e \n',qh,max(abs(Policy1alt(:)-Policy1Balt(:))));
vfoptions1.lowmemory=2;
[V1C,Policy1C,V1Calt,Policy1Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('%s lowmemory=2 (base), this should be zero: %.3e \n',qh,max(abs(V1(:)-V1C(:))))
fprintf('%s lowmemory=2 (base, Valt), this should be zero: %.3e \n',qh,max(abs(V1alt(:)-V1Calt(:))))
fprintf('%s lowmemory=2 (base, Policy), this should be zero: %.3e \n',qh,max(abs(Policy1(:)-Policy1C(:))))
fprintf('%s lowmemory=2 (base, Policyalt), this should be zero: %.3e \n',qh,max(abs(Policy1alt(:)-Policy1Calt(:))));
vfoptions1.lowmemory=0;

% Divide-and-conquer -> DC2A, should give the same answer as base
vfoptions2=vfoptions; vfoptions2.divideandconquer=1;
[V2,Policy2,V2alt,Policy2alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('%s DC2A, this should be zero: %.3e \n',qh,max(abs(V1(:)-V2(:))))
fprintf('%s DC2A (Valt), this should be zero: %.3e \n',qh,max(abs(V1alt(:)-V2alt(:))))
fprintf('%s DC2A (Policy), this should be zero: %.3e \n',qh,max(abs(Policy1(:)-Policy2(:))))
fprintf('%s DC2A (Policyalt), this should be zero: %.3e \n',qh,max(abs(Policy1alt(:)-Policy2alt(:))));
% ValueFnFromPolicy oracle on DC2A
vfoptionsVFP2=vfoptions2; vfoptionsVFP2.lowmemory=0;
vfoptionsVFP2.Policyalt=Policy2alt; % Naive QH: ValueFnFromPolicy reconstructs V from the exponential-discounter argmax
[V2fromPolicy,V2altfromPolicy]=ValueFnFromPolicy_FHorz(Policy2,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP2);
fprintf('%s ValueFnFromPolicy (DC2A), this should be zero: %.3e \n',qh,max(abs(V2fromPolicy(:)-V2(:))))
fprintf('%s ValueFnFromPolicy (DC2A, Valt), this should be zero: %.3e \n',qh,max(abs(V2altfromPolicy(:)-V2alt(:))))

% lowmemory on DC2A
vfoptions2.lowmemory=1;
[V2B,Policy2B,V2Balt,Policy2Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('%s lowmemory=1 (DC2A), this should be zero: %.3e \n',qh,max(abs(V2(:)-V2B(:))))
fprintf('%s lowmemory=1 (DC2A, Valt), this should be zero: %.3e \n',qh,max(abs(V2alt(:)-V2Balt(:))))
fprintf('%s lowmemory=1 (DC2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy2(:)-Policy2B(:))))
fprintf('%s lowmemory=1 (DC2A, Policyalt), this should be zero: %.3e \n',qh,max(abs(Policy2alt(:)-Policy2Balt(:))));
vfoptions2.lowmemory=2;
[V2C,Policy2C,V2Calt,Policy2Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('%s lowmemory=2 (DC2A), this should be zero: %.3e \n',qh,max(abs(V2(:)-V2C(:))))
fprintf('%s lowmemory=2 (DC2A, Valt), this should be zero: %.3e \n',qh,max(abs(V2alt(:)-V2Calt(:))))
fprintf('%s lowmemory=2 (DC2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy2(:)-Policy2C(:))))
fprintf('%s lowmemory=2 (DC2A, Policyalt), this should be zero: %.3e \n',qh,max(abs(Policy2alt(:)-Policy2Calt(:))));
vfoptions2.lowmemory=0;

% Grid interpolation -> GI2A (GI changes the solution slightly vs base, so no direct base equality check)
vfoptions3=vfoptions; vfoptions3.gridinterplayer=1; vfoptions3.ngridinterp=5;
[V3,Policy3,V3alt,Policy3alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
% ValueFnFromPolicy oracle on GI2A
vfoptionsVFP3=vfoptions3; vfoptionsVFP3.lowmemory=0;
vfoptionsVFP3.Policyalt=Policy3alt; % Naive QH: ValueFnFromPolicy reconstructs V from the exponential-discounter argmax
[V3fromPolicy,V3altfromPolicy]=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP3);
fprintf('%s ValueFnFromPolicy (GI2A), this should be zero: %.3e \n',qh,max(abs(V3fromPolicy(:)-V3(:))))
fprintf('%s ValueFnFromPolicy (GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V3altfromPolicy(:)-V3alt(:))))

% lowmemory on GI2A
vfoptions3.lowmemory=1;
[V3B,Policy3B,V3Balt,Policy3Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('%s lowmemory=1 (GI2A), this should be zero: %.3e \n',qh,max(abs(V3(:)-V3B(:))))
fprintf('%s lowmemory=1 (GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V3alt(:)-V3Balt(:))))
fprintf('%s lowmemory=1 (GI2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy3(:)-Policy3B(:))))
fprintf('%s lowmemory=1 (GI2A, Policyalt), this should be zero: %.3e \n',qh,max(abs(Policy3alt(:)-Policy3Balt(:))));
vfoptions3.lowmemory=2;
[V3C,Policy3C,V3Calt,Policy3Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('%s lowmemory=2 (GI2A), this should be zero: %.3e \n',qh,max(abs(V3(:)-V3C(:))))
fprintf('%s lowmemory=2 (GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V3alt(:)-V3Calt(:))))
fprintf('%s lowmemory=2 (GI2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy3(:)-Policy3C(:))))
fprintf('%s lowmemory=2 (GI2A, Policyalt), this should be zero: %.3e \n',qh,max(abs(Policy3alt(:)-Policy3Calt(:))));
vfoptions3.lowmemory=0;

% DC + GI -> DC2A_GI2A, should match GI2A
vfoptions4=vfoptions; vfoptions4.divideandconquer=1; vfoptions4.gridinterplayer=1; vfoptions4.ngridinterp=5;
[V4,Policy4,V4alt,Policy4alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('%s DC2A_GI2A vs GI2A, this should be zero: %.3e \n',qh,max(abs(V3(:)-V4(:))))
fprintf('%s DC2A_GI2A vs GI2A (Valt), this should be zero: %.3e \n',qh,max(abs(V3alt(:)-V4alt(:))))
fprintf('%s DC2A_GI2A vs GI2A (Policy), this should be zero: %.3e \n',qh,max(abs(Policy3(:)-Policy4(:))))
fprintf('%s DC2A_GI2A vs GI2A (Policyalt), this should be zero: %.3e \n',qh,max(abs(Policy3alt(:)-Policy4alt(:))));
% ValueFnFromPolicy oracle on DC2A_GI2A
vfoptionsVFP4=vfoptions4; vfoptionsVFP4.lowmemory=0;
vfoptionsVFP4.Policyalt=Policy4alt; % Naive QH: ValueFnFromPolicy reconstructs V from the exponential-discounter argmax
[V4fromPolicy,V4altfromPolicy]=ValueFnFromPolicy_FHorz(Policy4,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP4);
fprintf('%s ValueFnFromPolicy (DC2A_GI2A), this should be zero: %.3e \n',qh,max(abs(V4fromPolicy(:)-V4(:))))
fprintf('%s ValueFnFromPolicy (DC2A_GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V4altfromPolicy(:)-V4alt(:))))

% lowmemory on DC2A_GI2A
vfoptions4.lowmemory=1;
[V4B,Policy4B,V4Balt,Policy4Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('%s lowmemory=1 (DC2A_GI2A), this should be zero: %.3e \n',qh,max(abs(V4(:)-V4B(:))))
fprintf('%s lowmemory=1 (DC2A_GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V4alt(:)-V4Balt(:))))
fprintf('%s lowmemory=1 (DC2A_GI2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy4(:)-Policy4B(:))))
fprintf('%s lowmemory=1 (DC2A_GI2A, Policyalt), this should be zero: %.3e \n',qh,max(abs(Policy4alt(:)-Policy4Balt(:))));
vfoptions4.lowmemory=2;
[V4C,Policy4C,V4Calt,Policy4Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('%s lowmemory=2 (DC2A_GI2A), this should be zero: %.3e \n',qh,max(abs(V4(:)-V4C(:))))
fprintf('%s lowmemory=2 (DC2A_GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V4alt(:)-V4Calt(:))))
fprintf('%s lowmemory=2 (DC2A_GI2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy4(:)-Policy4C(:))))
fprintf('%s lowmemory=2 (DC2A_GI2A, Policyalt), this should be zero: %.3e \n',qh,max(abs(Policy4alt(:)-Policy4Calt(:))));
vfoptions4.lowmemory=0;

%% Sophisticated
vfoptions.quasi_hyperbolic='Sophisticated'; qh='Sophisticated';

% Base
vfoptions1=vfoptions;
[V1,Policy1,V1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% ValueFnFromPolicy oracle on the base method
vfoptionsVFP1=vfoptions1; vfoptionsVFP1.lowmemory=0;
[V1fromPolicy,V1altfromPolicy]=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP1);
fprintf('%s ValueFnFromPolicy, this should be zero: %.3e \n',qh,max(abs(V1fromPolicy(:)-V1(:))))
fprintf('%s ValueFnFromPolicy (Valt), this should be zero: %.3e \n',qh,max(abs(V1altfromPolicy(:)-V1alt(:))))

% lowmemory on base
vfoptions1.lowmemory=1;
[V1B,Policy1B,V1Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('%s lowmemory=1 (base), this should be zero: %.3e \n',qh,max(abs(V1(:)-V1B(:))))
fprintf('%s lowmemory=1 (base, Valt), this should be zero: %.3e \n',qh,max(abs(V1alt(:)-V1Balt(:))))
fprintf('%s lowmemory=1 (base, Policy), this should be zero: %.3e \n',qh,max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=2;
[V1C,Policy1C,V1Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('%s lowmemory=2 (base), this should be zero: %.3e \n',qh,max(abs(V1(:)-V1C(:))))
fprintf('%s lowmemory=2 (base, Valt), this should be zero: %.3e \n',qh,max(abs(V1alt(:)-V1Calt(:))))
fprintf('%s lowmemory=2 (base, Policy), this should be zero: %.3e \n',qh,max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=0;

% Divide-and-conquer -> DC2A, should give the same answer as base
vfoptions2=vfoptions; vfoptions2.divideandconquer=1;
[V2,Policy2,V2alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('%s DC2A, this should be zero: %.3e \n',qh,max(abs(V1(:)-V2(:))))
fprintf('%s DC2A (Valt), this should be zero: %.3e \n',qh,max(abs(V1alt(:)-V2alt(:))))
fprintf('%s DC2A (Policy), this should be zero: %.3e \n',qh,max(abs(Policy1(:)-Policy2(:))))
% ValueFnFromPolicy oracle on DC2A
vfoptionsVFP2=vfoptions2; vfoptionsVFP2.lowmemory=0;
[V2fromPolicy,V2altfromPolicy]=ValueFnFromPolicy_FHorz(Policy2,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP2);
fprintf('%s ValueFnFromPolicy (DC2A), this should be zero: %.3e \n',qh,max(abs(V2fromPolicy(:)-V2(:))))
fprintf('%s ValueFnFromPolicy (DC2A, Valt), this should be zero: %.3e \n',qh,max(abs(V2altfromPolicy(:)-V2alt(:))))

% lowmemory on DC2A
vfoptions2.lowmemory=1;
[V2B,Policy2B,V2Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('%s lowmemory=1 (DC2A), this should be zero: %.3e \n',qh,max(abs(V2(:)-V2B(:))))
fprintf('%s lowmemory=1 (DC2A, Valt), this should be zero: %.3e \n',qh,max(abs(V2alt(:)-V2Balt(:))))
fprintf('%s lowmemory=1 (DC2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=2;
[V2C,Policy2C,V2Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('%s lowmemory=2 (DC2A), this should be zero: %.3e \n',qh,max(abs(V2(:)-V2C(:))))
fprintf('%s lowmemory=2 (DC2A, Valt), this should be zero: %.3e \n',qh,max(abs(V2alt(:)-V2Calt(:))))
fprintf('%s lowmemory=2 (DC2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy2(:)-Policy2C(:))))
vfoptions2.lowmemory=0;

% Grid interpolation -> GI2A (GI changes the solution slightly vs base, so no direct base equality check)
vfoptions3=vfoptions; vfoptions3.gridinterplayer=1; vfoptions3.ngridinterp=5;
[V3,Policy3,V3alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
% ValueFnFromPolicy oracle on GI2A
vfoptionsVFP3=vfoptions3; vfoptionsVFP3.lowmemory=0;
[V3fromPolicy,V3altfromPolicy]=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP3);
fprintf('%s ValueFnFromPolicy (GI2A), this should be zero: %.3e \n',qh,max(abs(V3fromPolicy(:)-V3(:))))
fprintf('%s ValueFnFromPolicy (GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V3altfromPolicy(:)-V3alt(:))))

% lowmemory on GI2A
vfoptions3.lowmemory=1;
[V3B,Policy3B,V3Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('%s lowmemory=1 (GI2A), this should be zero: %.3e \n',qh,max(abs(V3(:)-V3B(:))))
fprintf('%s lowmemory=1 (GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V3alt(:)-V3Balt(:))))
fprintf('%s lowmemory=1 (GI2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=2;
[V3C,Policy3C,V3Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('%s lowmemory=2 (GI2A), this should be zero: %.3e \n',qh,max(abs(V3(:)-V3C(:))))
fprintf('%s lowmemory=2 (GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V3alt(:)-V3Calt(:))))
fprintf('%s lowmemory=2 (GI2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy3(:)-Policy3C(:))))
vfoptions3.lowmemory=0;

% DC + GI -> DC2A_GI2A, should match GI2A
vfoptions4=vfoptions; vfoptions4.divideandconquer=1; vfoptions4.gridinterplayer=1; vfoptions4.ngridinterp=5;
[V4,Policy4,V4alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('%s DC2A_GI2A vs GI2A, this should be zero: %.3e \n',qh,max(abs(V3(:)-V4(:))))
fprintf('%s DC2A_GI2A vs GI2A (Valt), this should be zero: %.3e \n',qh,max(abs(V3alt(:)-V4alt(:))))
fprintf('%s DC2A_GI2A vs GI2A (Policy), this should be zero: %.3e \n',qh,max(abs(Policy3(:)-Policy4(:))))
% ValueFnFromPolicy oracle on DC2A_GI2A
vfoptionsVFP4=vfoptions4; vfoptionsVFP4.lowmemory=0;
[V4fromPolicy,V4altfromPolicy]=ValueFnFromPolicy_FHorz(Policy4,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP4);
fprintf('%s ValueFnFromPolicy (DC2A_GI2A), this should be zero: %.3e \n',qh,max(abs(V4fromPolicy(:)-V4(:))))
fprintf('%s ValueFnFromPolicy (DC2A_GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V4altfromPolicy(:)-V4alt(:))))

% lowmemory on DC2A_GI2A
vfoptions4.lowmemory=1;
[V4B,Policy4B,V4Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('%s lowmemory=1 (DC2A_GI2A), this should be zero: %.3e \n',qh,max(abs(V4(:)-V4B(:))))
fprintf('%s lowmemory=1 (DC2A_GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V4alt(:)-V4Balt(:))))
fprintf('%s lowmemory=1 (DC2A_GI2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=2;
[V4C,Policy4C,V4Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('%s lowmemory=2 (DC2A_GI2A), this should be zero: %.3e \n',qh,max(abs(V4(:)-V4C(:))))
fprintf('%s lowmemory=2 (DC2A_GI2A, Valt), this should be zero: %.3e \n',qh,max(abs(V4alt(:)-V4Calt(:))))
fprintf('%s lowmemory=2 (DC2A_GI2A, Policy), this should be zero: %.3e \n',qh,max(abs(Policy4(:)-Policy4C(:))))
vfoptions4.lowmemory=0;

%% Versus exponential discounting
% (i) at the actual QH beta0, Naive's continuation value equals the exponential value function
vfoptionsE=vfoptions; vfoptionsE.exoticpreferences='None';
[Vexp,Policyexp]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsE);
vfoptionsN=vfoptions; vfoptionsN.quasi_hyperbolic='Naive';
[VnA,PolicynA,VnAalt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsN);
fprintf('(i) Naive continuation value == exponential (beta0=%g), should be zero: %.3e \n',Params.beta0,max(abs(VnAalt(:)-Vexp(:))))

% (ii) with beta0=1, Naive main value AND continuation value both equal exponential
% (iii) with beta0=1, Sophisticated main value AND continuation value both equal exponential
beta0_store=Params.beta0; Params.beta0=1;
[Vexp1,~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsE);
[VnA1,~,VnA1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsN);
vfoptionsS=vfoptions; vfoptionsS.quasi_hyperbolic='Sophisticated';
[VsS1,~,VsS1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsS);
fprintf('(ii) Naive V == exponential (beta0=1), should be zero: %.3e \n',max(abs(VnA1(:)-Vexp1(:))))
fprintf('(ii) Naive Valt == exponential (beta0=1), should be zero: %.3e \n',max(abs(VnA1alt(:)-Vexp1(:))))
fprintf('(iii) Sophisticated V == exponential (beta0=1), should be zero: %.3e \n',max(abs(VsS1(:)-Vexp1(:))))
fprintf('(iii) Sophisticated Valt == exponential (beta0=1), should be zero: %.3e \n',max(abs(VsS1alt(:)-Vexp1(:))))
Params.beta0=beta0_store;

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
