function output=EZRiskyAsset_nod1_noz_noe_semiz_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Epstein-Zin mirror of CoreFHorzRiskyAsset_nod1_noz_noe_semiz_noa1.m. The model is run three times:
%   Case 1: consumption-units (traditional EZ, vfoptions.EZutils=0)
%   Case 2: utility-units, positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%   Case 3: utility-units, negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% followed by the EZ special tests (gamma=1/phi collapse; EZriskaversion=0 collapse;
% EZoneminusbeta=1 vs manual scaling). The u-shock is part of the EZ certainty-equivalent:
% ONE joint CE over (u, semizprime) — the gamma=1/phi collapse only holds exactly
% if the implementation does this, so it is enforced by test.
% Plus special test (v): V_Jplus1 (V of period jstar as the terminal value fn of a shorter model; legs at jstar=round(3*N_j/4) and jstar=N_j).
% Plus special tests (vi)-(ix): (vi) survival probabilities (plumbing, collapses with sj, and
% cross-method with sj); (vii) warm-glow of bequests (De Nardi luxury-good form: terminal-only
% and with declining sj; the sj+warm-glow collapse oracles are DROPPED for riskyasset, see
% (vii).3; the N_j-1 warm-glow identity); (viii) EZmortalityriskaversion; (ix) EZoneminusbeta=2 with sj.
% Note: noa1 riskyasset has no divide-and-conquer and no grid interpolation layer (no a1 to
% refine), so unlike the CoreFHorzTests EZ subcodes there are no DC/GI legs here and no
% GI versions of the collapse tests (those appear in the withA1 subcodes).
% Note: despite noz_noe this is NOT a no-shock model (the u shock and semiz always provide
% genuine risk), so no 'EZ without shocks' warning is expected here.
% n_d=[n_d2,n_d3,n_d4] = [riskyshare, savings, dsemiz];   d_grid=[d2_grid; d3_grid; d4_grid]

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
n_z=0;
z_grid=[];
pi_z=[];
% Semiz
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptions.d_grid=d_grid;
% Riskyasset
vfoptions.riskyasset=1;
vfoptions.refine_d=[0,1,1,1]; % no d1, 1 d2 (riskyshare), 1 d3 (savings), 1 d4 (dsemiz)
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions.n_u=vfoptionsbaseline.n_u;
vfoptions.u_grid=vfoptionsbaseline.u_grid;
vfoptions.pi_u=vfoptionsbaseline.pi_u;
simoptions.riskyasset=vfoptions.riskyasset;
simoptions.refine_d=vfoptions.refine_d;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.n_u=vfoptions.n_u;
simoptions.u_grid=vfoptions.u_grid;
simoptions.pi_u=vfoptions.pi_u;
simoptions.a_grid=a_grid;

jequaloneDist=zeros(n_a_big,vfoptions.n_semiz,'gpuArray');
jequaloneDist(1,ceil(vfoptions.n_semiz/2))=1;

ReturnFn_cons=@(savings,dsemiz,a,semiz,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_cons_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_posU=@(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_positiveUtils_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_negU=@(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_negativeUtils_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% With riskyasset, need to include all d in FnsToEvaluate, even if they are not in the ReturnFn
FnsToEvaluate.assets=@(riskyshare,savings,dsemiz,a,semiz) a;
FnsToEvaluate.savings=@(riskyshare,savings,dsemiz,a,semiz) savings;
FnsToEvaluate.earnings=@(riskyshare,savings,dsemiz,a,semiz,w,kappa_j) w*kappa_j*semiz;


%% Case 1: consumption-units (traditional Epstein-Zin, vfoptions.EZutils=0)
vfoptions.exoticpreferences='EpsteinZin';
vfoptions.EZutils=0;
vfoptions.EZriskaversion='ezgamma';
vfoptions.EZeis='ezphi';
ReturnFn=ReturnFn_cons;

%%
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% V from Policy
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

clear V1 V1B Policy1 Policy1B V1fromPolicy

%% Big a_grid for moments
simoptions1.a_grid=a_grid_big;
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

fprintf('Big a_grid moments (cons-units)\n')
[AllStats1.assets.Mean, AllStats1.savings.Mean, AllStats1.earnings.Mean]

%% Graph
fig=figure(figure_c);
plot(1:1:N_j, AgeConditionalStats1.assets.Mean, 1:1:N_j, AgeConditionalStats1.savings.Mean, 1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('EZ cons-units, age-cond mean: assets / savings / earnings (semiz)')
legend('assets','savings','earnings')

clear V1b Policy1b StationaryDist1 AllStats1 AgeConditionalStats1
simoptions1.a_grid=a_grid;


%% Case 2: utility-units with POSITIVE-valued utility fn (vfoptions.EZutils=1, EZpositiveutility=1)
vfoptions.exoticpreferences='EpsteinZin';
vfoptions.EZutils=1;
vfoptions.EZpositiveutility=1;
vfoptions.EZriskaversion='ezrisk';
ReturnFn=ReturnFn_posU;

%%
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% V from Policy
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

clear V1 V1B Policy1 Policy1B V1fromPolicy

%% Big a_grid for moments
simoptions1.a_grid=a_grid_big;
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

fprintf('Big a_grid moments (positive utils)\n')
[AllStats1.assets.Mean, AllStats1.savings.Mean, AllStats1.earnings.Mean]

%% Graph
fig=figure(100+figure_c); % Case 2, positive-valued utils (Case 1 is figure_c)
plot(1:1:N_j, AgeConditionalStats1.assets.Mean, 1:1:N_j, AgeConditionalStats1.savings.Mean, 1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('EZ positive utils, age-cond mean: assets / savings / earnings (semiz)')
legend('assets','savings','earnings')

clear V1b Policy1b StationaryDist1 AllStats1 AgeConditionalStats1
simoptions1.a_grid=a_grid;


%% Case 3: utility-units with NEGATIVE-valued utility fn (vfoptions.EZutils=1, EZpositiveutility=0)
vfoptions.exoticpreferences='EpsteinZin';
vfoptions.EZutils=1;
vfoptions.EZpositiveutility=0;
vfoptions.EZriskaversion='ezrisk';
ReturnFn=ReturnFn_negU;

%%
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% V from Policy
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

clear V1 V1B Policy1 Policy1B V1fromPolicy

%% Big a_grid for moments
simoptions1.a_grid=a_grid_big;
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

fprintf('Big a_grid moments (negative utils)\n')
[AllStats1.assets.Mean, AllStats1.savings.Mean, AllStats1.earnings.Mean]

%% Graph
fig=figure(200+figure_c); % Case 3, negative-valued utils (Case 1 is figure_c)
plot(1:1:N_j, AgeConditionalStats1.assets.Mean, 1:1:N_j, AgeConditionalStats1.savings.Mean, 1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('EZ negative utils, age-cond mean: assets / savings / earnings (semiz)')
legend('assets','savings','earnings')

clear V1b Policy1b StationaryDist1 AllStats1 AgeConditionalStats1
simoptions1.a_grid=a_grid;


%% Special Epstein-Zin tests
% The u-shock makes these sharper than in CoreFHorzTests: the collapses only hold exactly if
% the u-expectation sits INSIDE the joint certainty-equivalent (one CE over (u,semizprime)),
% so a wrong placement of the u-expectation fails these loudly.

%% (i) Consumption-units: gamma=1/phi collapses Epstein-Zin to standard vNM expected utility
% With vfoptions.EZutils=0 and EZriskaversion=1/EZeis, the EZ value fn is an increasing transform
% of the vNM problem with period utility (x^(1-ezgamma))/(1-ezgamma), which is exactly
% EZRiskyReturnFn_negativeUtils with ezsigma=ezgamma. So: identical Policy (exactly), and
% V_EZ=((1-ezgamma)*V_vNM).^(1/(1-ezgamma)) up to roundoff in the transform (~1e-13 relative).
ezphi_store=Params.ezphi;
ezsigma_store=Params.ezsigma;
Params.ezphi=1/Params.ezgamma;
Params.ezsigma=Params.ezgamma;

vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
[V1a,Policy1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.exoticpreferences='None';
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
V1btransformed=((1-Params.ezgamma)*V1b).^(1/(1-Params.ezgamma));
fprintf('EZ gamma=1/phi collapse, Policy: should give zero: %2.8f \n',max(abs(Policy1a(:)-Policy1b(:))))
fprintf('EZ gamma=1/phi collapse, V after transform (relative): should be roughly 1e-13: %g \n',max(abs(V1a(:)-V1btransformed(:)))/max(abs(V1a(:))))

Params.ezphi=ezphi_store;
Params.ezsigma=ezsigma_store;
clear V1a V1b V1btransformed Policy1a Policy1b

%% (ii) Utility-units: EZriskaversion=0 adds no extra risk aversion, so collapses to standard vNM
% (tests the ezc3/ezc4 sign-handling and the 1+/-crisk branch for both signs of the utility fn)
ezrisk_store=Params.ezrisk;
Params.ezrisk=0;

% positive-valued utility fn
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=1;
vfoptions1.EZriskaversion='ezrisk';
[V2a,Policy2a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.exoticpreferences='None';
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZ (positive utils) with EZriskaversion=0: should give zero: %2.8f \n',max(abs(V2a(:)-V2b(:))))
fprintf('EZ (positive utils) with EZriskaversion=0: should give zero: %2.8f \n',max(abs(Policy2a(:)-Policy2b(:))))

% negative-valued utility fn
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
[V3a,Policy3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.exoticpreferences='None';
[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZ (negative utils) with EZriskaversion=0: should give zero: %2.8f \n',max(abs(V3a(:)-V3b(:))))
fprintf('EZ (negative utils) with EZriskaversion=0: should give zero: %2.8f \n',max(abs(Policy3a(:)-Policy3b(:))))

Params.ezrisk=ezrisk_store;
clear V2a V2b V3a V3b Policy2a Policy2b Policy3a Policy3b

%% (iii) GI versions of the collapse tests: NOT APPLICABLE for noa1 riskyasset (no a1 to
% interpolate). The withA1 subcodes have them.

%% (iv) vfoptions.EZoneminusbeta=1 versus manually scaling the return fn
% - consumption-units: EZ is homogeneous of degree 1 in the return, so EZoneminusbeta=1 with
%   return x should equal EZoneminusbeta=0 with return scaled by (1-beta)^(1/(1-1/ezphi)).
% - utility-units: V=(1-beta)*u+beta*(CE) is the same recursion as default with u scaled by (1-beta).
% (the scale factors enter the scaled return fns as parameters, since GPU arrayfun does
% not support anonymous functions that capture workspace variables)
Params.ezscalefactor=(1-Params.beta)^(1/(1-1/Params.ezphi));
Params.ezscalefactoru=1-Params.beta;

% consumption-units
ReturnFn_cons_scaled=@(savings,dsemiz,a,semiz,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactor) ezscalefactor*EZRiskyReturnFn_cons_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
vfoptions1.EZoneminusbeta=1;
[V4a,Policy4a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.EZoneminusbeta=0;
[V4b,Policy4b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons_scaled,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZoneminusbeta=1 (cons-units) vs manual scaling: should give zero: %2.8f \n',max(abs(V4a(:)-V4b(:))))
fprintf('EZoneminusbeta=1 (cons-units) vs manual scaling: should give zero: %2.8f \n',max(abs(Policy4a(:)-Policy4b(:))))

% utility-units, positive
ReturnFn_posU_scaled=@(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactoru) ezscalefactoru*EZRiskyReturnFn_positiveUtils_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=1;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.EZoneminusbeta=1;
[V5a,Policy5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.EZoneminusbeta=0;
[V5b,Policy5b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU_scaled,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZoneminusbeta=1 (positive utils) vs manual scaling: should give zero: %2.8f \n',max(abs(V5a(:)-V5b(:))))
fprintf('EZoneminusbeta=1 (positive utils) vs manual scaling: should give zero: %2.8f \n',max(abs(Policy5a(:)-Policy5b(:))))

% utility-units, negative
ReturnFn_negU_scaled=@(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactoru) ezscalefactoru*EZRiskyReturnFn_negativeUtils_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.EZoneminusbeta=1;
[V6a,Policy6a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.EZoneminusbeta=0;
[V6b,Policy6b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU_scaled,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZoneminusbeta=1 (negative utils) vs manual scaling: should give zero: %2.8f \n',max(abs(V6a(:)-V6b(:))))
fprintf('EZoneminusbeta=1 (negative utils) vs manual scaling: should give zero: %2.8f \n',max(abs(Policy6a(:)-Policy6b(:))))

clear V4a V4b V5a V5b V6a V6b Policy4a Policy4b Policy5a Policy5b Policy6a Policy6b

%% (v) V_Jplus1: use V of period jstar as the terminal value function of a shorter model
% Solve the model, then solve a shorter model that runs only periods 1,...,jstar-1, giving it
% vfoptions.V_Jplus1=V(:,:,jstar). V_Jplus1 is the value fn of period N_j+1 of the model being
% solved, so the shorter model has Njs=jstar-1 periods, and the age-dependent parameters are
% trimmed to length Njs (only agej and kappa_j: everything else entering the ReturnFn, aprimeFn
% and SemiExoStateFn is a scalar). V and Policy must then be identical to the original model for
% periods 1,...,jstar-1. Run for all three EZ cases. noa1 riskyasset has no DC/GI tiers, so the
% legs are the basic solve plus its lowmemory rung (lowmemory=1, mirroring the basic-solve
% sections above), at jstar=round(3*N_j/4) and again at jstar=N_j (the latter exercises the
% terminal-period V_Jplus1 branches of the EZ riskyasset semiz raws, and covers the retirement
% periods). Both legs use the same vfoptions, so one full-length solve serves both.
% Note: mewj is age-dependent, but is only used for the agent distribution, which is not
% computed here, so it is left alone.
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits ALL the baseline riskyasset settings (aprimeFn, refine_d,
    % n_u/u_grid/pi_u) and semiz settings (n_semiz/semiz_grid/SemiExoStateFn); the EZ fields left
    % on it by Case 3 are all overwritten below (EZeis and EZpositiveutility are only ever read
    % under the matching EZutils setting)
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
    end
    % Full-length solve (shared by both jstar legs)
    [Vfull,Policyfull]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);

    %% V_Jplus1, jstar=round(3*N_j/4)
    jstar=round(3*N_j/4);
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    vfoptionsjs=vfoptionsv;
    vfoptionsjs.V_Jplus1=Vfull(:,:,jstar);
    Vbase=Vfull(:,:,1:Njs);
    Policybase=Policyfull(:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;

    %% V_Jplus1, jstar=N_j (terminal-period V_Jplus1 branches; covers the retirement periods)
    jstar=N_j;
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    vfoptionsjs=vfoptionsv;
    vfoptionsjs.V_Jplus1=Vfull(:,:,jstar);
    Vbase=Vfull(:,:,1:Njs);
    Policybase=Policyfull(:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;
end

clear Vfull Vbase Vshort Policyfull Policybase Policyshort

%% (vi) Survival probabilities: vfoptions.survivalprobability
% sj is the probability of surviving from age j to j+1: it weights the continuation
% certainty-equivalent inside the EZ aggregator (combined with a warm-glow of bequests, 1-sj
% weights the bequest term; that is tested in (vii)). Params.sj is declining with sj(N_j)=0.
% Four legs (noa1 riskyasset has no DC and no GI, so the cross-method legs are the file's own
% lowmemory rungs and ValueFnFromPolicy, mirroring the basic-solve sections above):
%  (vi).1 plumbing: survivalprobability='sjones' with sjones=ones(1,N_j) must be identical to
%         the no-survivalprobability baseline (V and Policy exact; all three EZ cases).
%  (vi).2 cons-units gamma=1/phi collapse with sj: EZ with survivalprobability='sj' equals vNM
%         (negativeUtils return fn, ezsigma=ezgamma) with DiscountFactorParamNames={'beta','sj'}:
%         Policy exact, V via V_EZ=((1-ezgamma)*V_vNM).^(1/(1-ezgamma)).
%  (vi).3 utility-units EZriskaversion=0 collapse with sj (both utility signs): V and Policy
%         exact against vNM with DiscountFactorParamNames={'beta','sj'}.
%  (vi).4 cross-method with sj: basic==lowmemory=1, FromPolicy==V (all three EZ cases; the
%         first exercise of the sj path of the EZ riskyasset ValueFnFromPolicy).

% (vi).1 plumbing: survivalprobability of all ones is the no-survivalprobability baseline
Params2=Params;
Params2.sjones=ones(1,N_j);
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits ALL the baseline riskyasset and semiz settings (see the note in (v));
    % the EZ fields left on it by earlier sections are all overwritten below
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
    end
    [V1a,Policy1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    vfoptionssj=vfoptionsv;
    vfoptionssj.survivalprobability='sjones';
    [V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptionssj);
    fprintf('survivalprobability=ones vs no-survivalprobability [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1a(:)-V1b(:))))
    fprintf('survivalprobability=ones vs no-survivalprobability [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1a(:)-Policy1b(:))))
end
clear V1a V1b Policy1a Policy1b

% (vi).2 cons-units gamma=1/phi collapse with survival probabilities: EZ with sj (and no
% warm-glow, so death is worth zero on both sides) is vNM expected utility with age-dependent
% discounting beta*sj (no GI leg: noa1 riskyasset has no grid interpolation layer)
ezphi_store=Params.ezphi;
ezsigma_store=Params.ezsigma;
Params.ezphi=1/Params.ezgamma;
Params.ezsigma=Params.ezgamma;
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
vfoptions1.survivalprobability='sj';
[V1a,Policy1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions2=vfoptions;
vfoptions2.exoticpreferences='None';
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,{'beta','sj'},[],vfoptions2);
V1btransformed=((1-Params.ezgamma)*V1b).^(1/(1-Params.ezgamma));
fprintf('EZ gamma=1/phi collapse with sj, Policy: should give zero: %2.8f \n',max(abs(Policy1a(:)-Policy1b(:))))
fprintf('EZ gamma=1/phi collapse with sj, V after transform (relative): should be roughly 1e-13: %g \n',max(abs(V1a(:)-V1btransformed(:)))/max(abs(V1a(:))))
Params.ezphi=ezphi_store;
Params.ezsigma=ezsigma_store;
clear V1a V1b V1btransformed Policy1a Policy1b

% (vi).3 utility-units EZriskaversion=0 collapse with survival probabilities (both utility signs)
ezrisk_store=Params.ezrisk;
Params.ezrisk=0;
% positive-valued utility fn
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=1;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.survivalprobability='sj';
[V2a,Policy2a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions2=vfoptions;
vfoptions2.exoticpreferences='None';
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,{'beta','sj'},[],vfoptions2);
fprintf('EZ (positive utils) with EZriskaversion=0 and sj: should give zero: %2.8f \n',max(abs(V2a(:)-V2b(:))))
fprintf('EZ (positive utils) with EZriskaversion=0 and sj: should give zero: %2.8f \n',max(abs(Policy2a(:)-Policy2b(:))))
% negative-valued utility fn
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.survivalprobability='sj';
[V3a,Policy3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions2=vfoptions;
vfoptions2.exoticpreferences='None';
[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,{'beta','sj'},[],vfoptions2);
fprintf('EZ (negative utils) with EZriskaversion=0 and sj: should give zero: %2.8f \n',max(abs(V3a(:)-V3b(:))))
fprintf('EZ (negative utils) with EZriskaversion=0 and sj: should give zero: %2.8f \n',max(abs(Policy3a(:)-Policy3b(:))))
Params.ezrisk=ezrisk_store;
clear V2a V2b V3a V3b Policy2a Policy2b Policy3a Policy3b

% (vi).4 cross-method with survival probabilities: basic==lowmemory=1, FromPolicy==V
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits ALL the baseline riskyasset and semiz settings (see the note in (v));
    % the EZ fields left on it by earlier sections are all overwritten below
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
    end
    vfoptionsv.survivalprobability='sj';
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    vfoptionsv.lowmemory=1;
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('lowmemory=1 with sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('lowmemory=1 with sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    vfoptionsv.lowmemory=0;
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsv);
    fprintf('ValueFnFromPolicy with sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
end
clear V1 V2 Policy1 Policy2 V1fromPolicy

%% (vii) Warm-glow of bequests: vfoptions.WarmGlowBequestsFn (De Nardi luxury-good form)
% The warm-glow fns (in EZRisky_ReturnFns) match each case's units:
%   cons-units:     EZRiskyWarmGlowFn_cons(aprime,wg1,wg2) (a consumption-equivalent, strictly
%                   positive; curvature comes from the EZ preferences)
%   positive utils: EZRiskyWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3) (strictly positive)
%   negative utils: EZRiskyWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3) (strictly negative)
% In the riskyasset case the warm-glow argument is a2prime, the risky-asset holding realized
% via the (d,u) lottery: the toolkit evaluates the warm-glow on the a2 grid and composes it
% through the same lottery as an interior continuation value.
% Five legs:
%  (vii).1 terminal-only default (WarmGlowBequestsFn set, no survivalprobability; the
%          dispatcher prints its assumed-terminal-only warning, which is expected output
%          here): basic==lowmemory=1, FromPolicy==V; and identical to the explicit
%          survivalprobability='sjterm' with sjterm=[ones(1,N_j-1),0].
%  (vii).2 warm-glow with declining sj: basic==lowmemory=1, FromPolicy==V.
%  (vii).3 exact collapse oracles with sj and warm-glow: DROPPED for riskyasset (see below).
%  (vii).4 V_Jplus1 mini-leg with sj and warm-glow (jstar=round(2*N_j/3); basic and lowmemory).
%  (vii).5 N_j-1 warm-glow identity (all three EZ cases; see below).

% (vii).1 terminal-only default
Params2=Params;
Params2.sjterm=[ones(1,N_j-1),0];
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits ALL the baseline riskyasset and semiz settings (see the note in (v));
    % the EZ fields left on it by earlier sections are all overwritten below
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
        vfoptionsv.WarmGlowBequestsFn=@(aprime,wg1,wg2) EZRiskyWarmGlowFn_cons(aprime,wg1,wg2);
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
        vfoptionsv.WarmGlowBequestsFn=@(aprime,wg1,wg2,wg3) EZRiskyWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3);
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
        vfoptionsv.WarmGlowBequestsFn=@(aprime,wg1,wg2,wg3) EZRiskyWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);
    end
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    vfoptionsv.lowmemory=1;
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('lowmemory=1 with terminal warm-glow [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('lowmemory=1 with terminal warm-glow [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    vfoptionsv.lowmemory=0;
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsv);
    fprintf('ValueFnFromPolicy with terminal warm-glow [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    % identical to the explicit terminal-only survival probabilities
    vfoptionssjt=vfoptionsv;
    vfoptionssjt.survivalprobability='sjterm';
    [V5,Policy5]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptionssjt);
    fprintf('terminal-only warm-glow vs explicit sjterm [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V5(:))))
    fprintf('terminal-only warm-glow vs explicit sjterm [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy5(:))))
end
clear V1 V2 V5 Policy1 Policy2 Policy5 V1fromPolicy

% (vii).2 warm-glow with declining survival probabilities
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits ALL the baseline riskyasset and semiz settings (see the note in (v));
    % the EZ fields left on it by earlier sections are all overwritten below
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
        vfoptionsv.WarmGlowBequestsFn=@(aprime,wg1,wg2) EZRiskyWarmGlowFn_cons(aprime,wg1,wg2);
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
        vfoptionsv.WarmGlowBequestsFn=@(aprime,wg1,wg2,wg3) EZRiskyWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3);
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
        vfoptionsv.WarmGlowBequestsFn=@(aprime,wg1,wg2,wg3) EZRiskyWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);
    end
    vfoptionsv.survivalprobability='sj';
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    vfoptionsv.lowmemory=1;
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('lowmemory=1 with warm-glow and sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('lowmemory=1 with warm-glow and sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    vfoptionsv.lowmemory=0;
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsv);
    fprintf('ValueFnFromPolicy with warm-glow and sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
end
clear V1 V2 Policy1 Policy2 V1fromPolicy

% (vii).3 exact collapse oracles with sj and warm-glow: DROPPED for riskyasset.
% The CoreFHorzTests EZ subcodes build the vNM oracle by adding beta*oneminussj*WarmGlowFn(aprime)
% to the period return fn, but a riskyasset ReturnFn has no aprime argument: a2prime is realized
% via the (d,u) lottery AFTER the return fn is evaluated, so the warm-glow term cannot be
% composed into a riskyasset return fn (and the vNM riskyasset raws have no warm-glow support
% to lean on instead). Exactness of the warm-glow composition here is anchored by the
% degenerateu bridge in the cross tests and by the N_j-1 warm-glow identity in (vii).5.

% (vii).4 V_Jplus1 mini-leg with sj and warm-glow: V of period jstar as the terminal value fn
% of a shorter model, with survivalprobability and WarmGlowBequestsFn active (as (v) but only
% one jstar, with the basic solve plus its lowmemory rungs). The age-dependent parameters
% trimmed to Njs now include sj (and oneminussj).
jstar=round(2*N_j/3);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
Paramsjs.sj=Params.sj(1:Njs);
Paramsjs.oneminussj=Params.oneminussj(1:Njs);
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits ALL the baseline riskyasset and semiz settings (see the note in (v));
    % the EZ fields left on it by earlier sections are all overwritten below
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
        vfoptionsv.WarmGlowBequestsFn=@(aprime,wg1,wg2) EZRiskyWarmGlowFn_cons(aprime,wg1,wg2);
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
        vfoptionsv.WarmGlowBequestsFn=@(aprime,wg1,wg2,wg3) EZRiskyWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3);
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
        vfoptionsv.WarmGlowBequestsFn=@(aprime,wg1,wg2,wg3) EZRiskyWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);
    end
    vfoptionsv.survivalprobability='sj';
    [Vfull,Policyfull]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    vfoptionsjs=vfoptionsv;
    vfoptionsjs.V_Jplus1=Vfull(:,:,jstar);
    Vbase=Vfull(:,:,1:Njs);
    Policybase=Policyfull(:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 with sj and warm-glow (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 with sj and warm-glow (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 with sj and warm-glow, lowmemory=1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 with sj and warm-glow, lowmemory=1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;
end
clear Vfull Vbase Vshort Policyfull Policybase Policyshort

% (vii).5 N_j-1 warm-glow identity (ALL THREE EZ CASES)
% A full N_j-period solve with NO warm-glow and NO survival probabilities must equal an
% (N_j-1)-period solve whose WarmGlowBequestsFn is the closed-form terminal value fn of the
% full model (no survivalprobability, so the default terminal warm-glow weight is exactly 1
% and the raws compose it through the (d,u) lottery exactly like an interior continuation
% value; the dispatcher's assumed-terminal-only warning is expected output here too).
% cons-units now included: the terminal warm-glow enters INSIDE the ^ezc7 root (Kraft/Munk/Weiss 2022), so its warm-glow fn is the terminal value fn itself, which for cons-units is just the composite good x ((ezc1*F^ezc2)^ezc7=F for ezc1=1).
% Terminal value derivation for this shape: age N_j is retirement (agej>=Jr), so the budget
% is c=pension+a-savings, independent of the semiz state; V_Nj is attained at
% savings=0 (d3_grid(1)=0, consume everything) and search effort dsemiz=0 (d4_grid(1)=0
% maximizes the (1-searcheffortcost*dsemiz) factor), while riskyshare does not enter the
% return fn, so the terminal composite is x=pension+a, giving
%   cons-units: V_Nj(a)=pension+a (the composite good itself: (ezc1*F^ezc2)^ezc7=F for ezc1=1)
%   positive utils: V_Nj(a)=((1+pension+a)^(1-ezsigma)-1)/(1-ezsigma)
%   negative utils: V_Nj(a)=((pension+a)^(1-ezsigma))/(1-ezsigma)
% all independent of the shocks, so the joint certainty-equivalent over (u,semizprime) reduces
% to the u-lottery over the warm-glow argument a2prime and the identity is exact for V and
% Policy at every age 1,...,N_j-1.
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:N_j-1);
Paramsjs.kappa_j=Params.kappa_j(1:N_j-1);
% consumption-units (traditional Epstein-Zin)
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
[Vfull,Policyfull]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
Vfull=Vfull(:,:,1:N_j-1);
Policyfull=Policyfull(:,:,:,1:N_j-1);
vfoptions1nj=vfoptions1;
vfoptions1nj.WarmGlowBequestsFn=@(aprime,pension) pension+aprime;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j-1,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Paramsjs,DiscountFactorParamNames,[],vfoptions1nj);
fprintf('N_j-1 warm-glow identity [EZ cons-units], this should be zero: %2.8f \n',max(abs(Vfull(:)-Vshort(:))))
fprintf('N_j-1 warm-glow identity [EZ cons-units], this should be zero: %2.8f \n',max(abs(Policyfull(:)-Policyshort(:))))
% positive-valued utility fn
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=1;
vfoptions1.EZriskaversion='ezrisk';
[Vfull,Policyfull]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
Vfull=Vfull(:,:,1:N_j-1);
Policyfull=Policyfull(:,:,:,1:N_j-1);
vfoptions1nj=vfoptions1;
vfoptions1nj.WarmGlowBequestsFn=@(aprime,pension,ezsigma) ((1+pension+aprime)^(1-ezsigma)-1)/(1-ezsigma);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j-1,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Paramsjs,DiscountFactorParamNames,[],vfoptions1nj);
fprintf('N_j-1 warm-glow identity [EZ positive utils], this should be zero: %2.8f \n',max(abs(Vfull(:)-Vshort(:))))
fprintf('N_j-1 warm-glow identity [EZ positive utils], this should be zero: %2.8f \n',max(abs(Policyfull(:)-Policyshort(:))))
% negative-valued utility fn
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
[Vfull,Policyfull]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
Vfull=Vfull(:,:,1:N_j-1);
Policyfull=Policyfull(:,:,:,1:N_j-1);
vfoptions1nj=vfoptions1;
vfoptions1nj.WarmGlowBequestsFn=@(aprime,pension,ezsigma) ((pension+aprime)^(1-ezsigma))/(1-ezsigma);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j-1,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Paramsjs,DiscountFactorParamNames,[],vfoptions1nj);
fprintf('N_j-1 warm-glow identity [EZ negative utils], this should be zero: %2.8f \n',max(abs(Vfull(:)-Vshort(:))))
fprintf('N_j-1 warm-glow identity [EZ negative utils], this should be zero: %2.8f \n',max(abs(Policyfull(:)-Policyshort(:))))
clear Vfull Vshort Policyfull Policyshort

%% (viii) EZmortalityriskaversion (mortality risk aversion distinct from within-period risk aversion)
% All legs use survivalprobability='sj': with sj identically 1 there is no mortality risk for
% the mortality risk aversion to act on, so it would be vacuous (ezc8 only enters the raws
% through the sj-weighted continuation/warm-glow terms). Two legs:
%  (viii).1 identity: EZmortalityriskaversion equal to the within-period risk-aversion
%           coefficient (ezgamma for cons-units, ezrisk for utility-units) gives ezc8=1:
%           identical to leaving EZmortalityriskaversion unset (all three EZ cases).
%  (viii).2 distinct mortality risk aversion (ezmrisk=5): basic==lowmemory=1, FromPolicy==V.

% (viii).1 identity: mortality risk aversion equal to the within-period risk aversion
Params2=Params;
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits ALL the baseline riskyasset and semiz settings (see the note in (v));
    % the EZ fields left on it by earlier sections are all overwritten below
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
    end
    vfoptionsv.survivalprobability='sj';
    if ezcase==1
        Params2.ezmriskeq=Params.ezgamma; % equal to the cons-units risk aversion
    else
        Params2.ezmriskeq=Params.ezrisk; % equal to the utility-units risk aversion
    end
    [V1a,Policy1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    vfoptionsmr=vfoptionsv;
    vfoptionsmr.EZmortalityriskaversion='ezmriskeq';
    [V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptionsmr);
    fprintf('EZmortalityriskaversion=within-period risk aversion vs unset [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1a(:)-V1b(:))))
    fprintf('EZmortalityriskaversion=within-period risk aversion vs unset [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1a(:)-Policy1b(:))))
end
clear V1a V1b Policy1a Policy1b

% (viii).2 distinct mortality risk aversion (Params.ezmrisk=5)
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits ALL the baseline riskyasset and semiz settings (see the note in (v));
    % the EZ fields left on it by earlier sections are all overwritten below
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
    end
    vfoptionsv.survivalprobability='sj';
    vfoptionsv.EZmortalityriskaversion='ezmrisk';
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    vfoptionsv.lowmemory=1;
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('lowmemory=1 with sj and mortality risk aversion [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('lowmemory=1 with sj and mortality risk aversion [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    vfoptionsv.lowmemory=0;
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsv);
    fprintf('ValueFnFromPolicy with sj and mortality risk aversion [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
end
clear V1 V2 Policy1 Policy2 V1fromPolicy

%% (ix) vfoptions.EZoneminusbeta=2 (with survival probabilities)
% EZoneminusbeta=2 puts (1-sj*beta) on the this-period return: the mortality-adjusted version
% of EZoneminusbeta=1. As in (iv):
% - cons-units: equals EZoneminusbeta=0 with the return scaled by the AGE-DEPENDENT factor
%   (1-sj*beta).^(1/(1-1/ezphi)) (the aggregator is homogeneous of degree 1 within each
%   period, so the equivalence holds age by age by backward induction).
% - utility-units: V=(1-sj*beta)*u+beta*(sj-weighted CE) is the same recursion as the
%   default with u scaled by (1-sj*beta).
% (the scale factors enter the scaled return fns as age-dependent parameters, since GPU
% arrayfun does not support anonymous functions that capture workspace variables)
Params2=Params;
Params2.ezscalefactor2=(1-Params.sj*Params.beta).^(1/(1-1/Params.ezphi));
Params2.ezscalefactoru2=1-Params.sj*Params.beta;

% consumption-units
ReturnFn_cons_scaled2=@(savings,dsemiz,a,semiz,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactor2) ezscalefactor2*EZRiskyReturnFn_cons_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
vfoptions1.survivalprobability='sj';
vfoptions1.EZoneminusbeta=2;
[V4a,Policy4a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params2,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.EZoneminusbeta=0;
[V4b,Policy4b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons_scaled2,Params2,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZoneminusbeta=2 with sj (cons-units) vs manual scaling: should give zero: %2.8f \n',max(abs(V4a(:)-V4b(:))))
fprintf('EZoneminusbeta=2 with sj (cons-units) vs manual scaling: should give zero: %2.8f \n',max(abs(Policy4a(:)-Policy4b(:))))

% utility-units, positive
ReturnFn_posU_scaled2=@(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactoru2) ezscalefactoru2*EZRiskyReturnFn_positiveUtils_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=1;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.survivalprobability='sj';
vfoptions1.EZoneminusbeta=2;
[V5a,Policy5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params2,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.EZoneminusbeta=0;
[V5b,Policy5b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU_scaled2,Params2,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZoneminusbeta=2 with sj (positive utils) vs manual scaling: should give zero: %2.8f \n',max(abs(V5a(:)-V5b(:))))
fprintf('EZoneminusbeta=2 with sj (positive utils) vs manual scaling: should give zero: %2.8f \n',max(abs(Policy5a(:)-Policy5b(:))))

% utility-units, negative
ReturnFn_negU_scaled2=@(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactoru2) ezscalefactoru2*EZRiskyReturnFn_negativeUtils_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.survivalprobability='sj';
vfoptions1.EZoneminusbeta=2;
[V6a,Policy6a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params2,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.EZoneminusbeta=0;
[V6b,Policy6b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU_scaled2,Params2,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZoneminusbeta=2 with sj (negative utils) vs manual scaling: should give zero: %2.8f \n',max(abs(V6a(:)-V6b(:))))
fprintf('EZoneminusbeta=2 with sj (negative utils) vs manual scaling: should give zero: %2.8f \n',max(abs(Policy6a(:)-Policy6b(:))))

clear V4a V4b V5a V5b V6a V6b Policy4a Policy4b Policy5a Policy5b Policy6a Policy6b

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
