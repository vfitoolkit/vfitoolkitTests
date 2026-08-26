function output=EZRiskyAsset_d1_z_e_nosemiz_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Epstein-Zin mirror of CoreFHorzRiskyAsset_d1_z_e_nosemiz_noa1.m. The model is run three times:
%   Case 1: consumption-units (traditional EZ, vfoptions.EZutils=0)
%   Case 2: utility-units, positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%   Case 3: utility-units, negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% followed by the EZ special tests (gamma=1/phi collapse; EZriskaversion=0 collapse;
% EZoneminusbeta=1 vs manual scaling). The u-shock is part of the EZ certainty-equivalent:
% ONE joint CE over (u, zprime, eprime) — the gamma=1/phi collapse only holds exactly if the
% implementation does this, so it is enforced by test.
% Section (v) then tests vfoptions.V_Jplus1: a shorter model handed V of period jstar as its
% terminal value fn must reproduce the full solve (the jstar=N_j leg runs the raws' V_Jplus1 terminal branches).
% Sections (vi)-(ix) then run the sj/warm-glow special tests:
%   (vi)   survival probabilities: sjones plumbing; the (i)/(ii) collapse oracles with the declining
%          sj (the vNM oracle discounts by DiscountFactorParamNames={'beta','sj'}); FromPolicy
%          and lowmemory agreement with sj
%   (vii)  warm-glow of bequests (De Nardi luxury-good form): terminal-only default (vs explicit
%          sjterm) and with the declining sj; a V_Jplus1 mini-leg; the N_j-1 terminal warm-glow
%          identity (all three EZ cases). The vNM collapse oracles with warm-glow are DROPPED for
%          riskyasset (no vNM warm-glow support; WG(aprime) cannot enter a riskyasset ReturnFn).
%   (viii) EZmortalityriskaversion: identity when set equal to the case's own risk aversion; a
%          distinct mortality risk aversion (ezmrisk=5) with FromPolicy and lowmemory agreement
%   (ix)   EZoneminusbeta=2 versus manually scaling the return fn by the age-dependent (1-sj*beta) factor
% Note: noa1 riskyasset has no divide-and-conquer and no grid interpolation layer (no a1 to
% refine), so unlike the CoreFHorzTests EZ subcodes there are no DC/GI legs here and no
% GI versions of the collapse tests (those appear in the withA1 subcodes).
% n_d=[n_d1,n_d2,n_d3] = [labour, riskyshare, savings];   d_grid=[d1_grid; d2_grid; d3_grid]

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=vfoptions.n_e;
simoptions.e_grid=vfoptions.e_grid;
simoptions.pi_e=vfoptions.pi_e;
% Riskyasset
vfoptions.riskyasset=1;
vfoptions.refine_d=[1,1,1]; % 1 d1 (labour), 1 d2 (riskyshare), 1 d3 (savings)
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
simoptions.d_grid=d_grid;

jequaloneDist=zeros([n_a_big,n_z,vfoptions.n_e],'gpuArray');
jequaloneDist(1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;

ReturnFn_cons=@(h,savings,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension) EZRiskyReturnFn_cons_d1_z_e_nosemiz(h,savings,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension);
ReturnFn_posU=@(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZRiskyReturnFn_positiveUtils_d1_z_e_nosemiz(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
ReturnFn_negU=@(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZRiskyReturnFn_negativeUtils_d1_z_e_nosemiz(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);

% With riskyasset, need to include all d in FnsToEvaluate, even if they are not in the ReturnFn
FnsToEvaluate.assets=@(h,riskyshare,savings,a,z,e) a;
FnsToEvaluate.savings=@(h,riskyshare,savings,a,z,e) savings;
FnsToEvaluate.hours=@(h,riskyshare,savings,a,z,e) h;
FnsToEvaluate.earnings=@(h,riskyshare,savings,a,z,e,w,kappa_j) w*kappa_j*h*z*e;


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

vfoptions1.lowmemory=2;
[V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(V1(:)-V1C(:))))
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=0;

clear V1 V1B V1C Policy1 Policy1B Policy1C V1fromPolicy

%% Big a_grid for moments
simoptions1.a_grid=a_grid_big;
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

fprintf('Big a_grid moments (cons-units)\n')
[AllStats1.assets.Mean, AllStats1.savings.Mean, AllStats1.hours.Mean, AllStats1.earnings.Mean]

%% Graph
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:N_j, AgeConditionalStats1.assets.Mean, 1:1:N_j, AgeConditionalStats1.savings.Mean)
title('EZ cons-units, age-cond mean: assets / savings (d1+z+e)')
legend('assets','savings')
subplot(2,1,2); plot(1:1:N_j, AgeConditionalStats1.hours.Mean, 1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('EZ cons-units, age-cond mean: hours / earnings')
legend('hours','earnings')

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

vfoptions1.lowmemory=2;
[V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(V1(:)-V1C(:))))
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=0;

clear V1 V1B V1C Policy1 Policy1B Policy1C V1fromPolicy

%% Big a_grid for moments
simoptions1.a_grid=a_grid_big;
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

fprintf('Big a_grid moments (positive utils)\n')
[AllStats1.assets.Mean, AllStats1.savings.Mean, AllStats1.hours.Mean, AllStats1.earnings.Mean]

%% Graph
fig=figure(100+figure_c); % Case 2, positive-valued utils (Case 1 is figure_c)
subplot(2,1,1); plot(1:1:N_j, AgeConditionalStats1.assets.Mean, 1:1:N_j, AgeConditionalStats1.savings.Mean)
title('EZ positive utils, age-cond mean: assets / savings (d1+z+e)')
legend('assets','savings')
subplot(2,1,2); plot(1:1:N_j, AgeConditionalStats1.hours.Mean, 1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('EZ positive utils, age-cond mean: hours / earnings')
legend('hours','earnings')

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

vfoptions1.lowmemory=2;
[V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(V1(:)-V1C(:))))
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=0;

clear V1 V1B V1C Policy1 Policy1B Policy1C V1fromPolicy

%% Big a_grid for moments
simoptions1.a_grid=a_grid_big;
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

fprintf('Big a_grid moments (negative utils)\n')
[AllStats1.assets.Mean, AllStats1.savings.Mean, AllStats1.hours.Mean, AllStats1.earnings.Mean]

%% Graph
fig=figure(200+figure_c); % Case 3, negative-valued utils (Case 1 is figure_c)
subplot(2,1,1); plot(1:1:N_j, AgeConditionalStats1.assets.Mean, 1:1:N_j, AgeConditionalStats1.savings.Mean)
title('EZ negative utils, age-cond mean: assets / savings (d1+z+e)')
legend('assets','savings')
subplot(2,1,2); plot(1:1:N_j, AgeConditionalStats1.hours.Mean, 1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('EZ negative utils, age-cond mean: hours / earnings')
legend('hours','earnings')

clear V1b Policy1b StationaryDist1 AllStats1 AgeConditionalStats1
simoptions1.a_grid=a_grid;


%% Special Epstein-Zin tests
% The u-shock makes these sharper than in CoreFHorzTests: the collapses only hold exactly if
% the u-expectation sits INSIDE the joint certainty-equivalent (one CE over (u,zprime,eprime)), so a
% wrong placement of the u-expectation fails these loudly.

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
ReturnFn_cons_scaled=@(h,savings,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension,ezscalefactor) ezscalefactor*EZRiskyReturnFn_cons_d1_z_e_nosemiz(h,savings,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension);
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
ReturnFn_posU_scaled=@(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension,ezscalefactoru) ezscalefactoru*EZRiskyReturnFn_positiveUtils_d1_z_e_nosemiz(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
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
ReturnFn_negU_scaled=@(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension,ezscalefactoru) ezscalefactoru*EZRiskyReturnFn_negativeUtils_d1_z_e_nosemiz(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
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
% vfoptions.V_Jplus1=V(:,:,:,jstar). V_Jplus1 is the value fn of period N_j+1 of the model being
% solved, so the shorter model has Njs=jstar-1 periods, and the age-dependent parameters are
% trimmed to length Njs. V and Policy must then be identical to the original model for periods
% 1,...,jstar-1. Two jstar values: jstar=round(3*N_j/4), and jstar=N_j so that the shorter model
% also covers the retirement periods (and is the first-ever run of the EZ riskyasset raws'
% V_Jplus1 terminal branches). Each jstar is run with the same lowmemory rungs as the standard
% sections above (lowmemory=0, 1 and 2). Run for all three EZ cases.
% Note: mewj is age-dependent, but is only used for the agent distribution, which is not
% computed here, so it is left alone. The aprimeFn parameter (r) is not age-dependent, so
% only agej and kappa_j get trimmed.
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits all the riskyasset settings (refine_d, aprimeFn, n_u, u_grid, pi_u) and the e-shock (n_e, e_grid, pi_e)
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
    % One full solve serves both jstar values (the vfoptions are identical for both)
    [Vfull,Policyfull]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    for jstar=[round(3*N_j/4),N_j]
        Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
        Paramsjs=Params;
        Paramsjs.agej=Params.agej(1:Njs);
        Paramsjs.kappa_j=Params.kappa_j(1:Njs);
        vfoptionsjs=vfoptionsv;
        vfoptionsjs.V_Jplus1=Vfull(:,:,:,jstar);
        Vbase=Vfull(:,:,:,1:Njs);
        Policybase=Policyfull(:,:,:,:,1:Njs);
        [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
        fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
        fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
        % lowmemory (the V_Jplus1 branch of the raws has its own lowmemory sub-branches)
        vfoptionsjs.lowmemory=1;
        [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
        fprintf('V_Jplus1 (jstar=%i), lowmemory=1 [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
        fprintf('V_Jplus1 (jstar=%i), lowmemory=1 [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
        vfoptionsjs.lowmemory=0;
        vfoptionsjs.lowmemory=2;
        [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
        fprintf('V_Jplus1 (jstar=%i), lowmemory=2 [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
        fprintf('V_Jplus1 (jstar=%i), lowmemory=2 [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
        vfoptionsjs.lowmemory=0;
    end
end

clear Vfull Policyfull Vbase Vshort Policybase Policyshort

%% (vi) Survival probabilities (vfoptions.survivalprobability)
% Params.sj is declining with sj(N_j)=0 (set in CoreFHorzRiskyAssetEZTests.m). When
% vfoptions.survivalprobability is not set the EZ codes use sj=ones(N_j,1) internally, so
% survivalprobability='sjones' (all ones) must reproduce the baseline solve exactly (pure
% plumbing). Then the (i)/(ii) collapse oracles are repeated WITH the declining sj: once EZ
% collapses to vNM, the survival probability is just an age-dependent discount factor, so the
% vNM oracle uses DiscountFactorParamNames={'beta','sj'} (multiplicative discount factors -- no
% vNM warm-glow is involved).
Params2=Params;
Params2.sjones=ones(1,N_j);
DiscountFactorParamNames2={'beta','sj'}; % beta*sj: age-dependent discounting for the vNM oracles

% (vi).1 plumbing: survivalprobability='sjones' vs not setting survivalprobability, all three EZ cases
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits all the riskyasset settings (refine_d, aprimeFn, n_u, u_grid, pi_u) and the e-shock (n_e, e_grid, pi_e)
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
    [V1a,Policy1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptionsv);
    vfoptionsv.survivalprobability='sjones';
    [V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('survivalprobability all-ones vs unset [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1a(:)-V1b(:))))
    fprintf('survivalprobability all-ones vs unset [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1a(:)-Policy1b(:))))
end
clear V1a V1b Policy1a Policy1b

% (vi).2 cons-units gamma=1/phi collapse WITH the declining sj (as (i), plus survivalprobability='sj';
% the vNM oracle discounts by beta*sj). Policy exact; V via the usual transform.
% the grid interpolation layer repeat of the main-EZ-bank version is replaced by the file's
% own lowmemory rungs (noa1 riskyasset has no DC and no GI).
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
vfoptions1=rmfield(vfoptions1,'survivalprobability');
vfoptions1.exoticpreferences='None';
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames2,[],vfoptions1);
V1btransformed=((1-Params.ezgamma)*V1b).^(1/(1-Params.ezgamma));
fprintf('EZ gamma=1/phi collapse with sj, Policy [EZ cons-units]: should give zero: %2.8f \n',max(abs(Policy1a(:)-Policy1b(:))))
fprintf('EZ gamma=1/phi collapse with sj, V after transform (relative) [EZ cons-units]: should be roughly 1e-13: %g \n',max(abs(V1a(:)-V1btransformed(:)))/max(abs(V1a(:))))
% the same on the lowmemory=1 rung (both sides of the collapse solved with lowmemory=1)
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
vfoptions1.survivalprobability='sj';
vfoptions1.lowmemory=1;
[V1c,Policy1c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1=rmfield(vfoptions1,'survivalprobability');
vfoptions1.exoticpreferences='None';
[V1d,Policy1d]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames2,[],vfoptions1);
V1dtransformed=((1-Params.ezgamma)*V1d).^(1/(1-Params.ezgamma));
fprintf('EZ gamma=1/phi collapse with sj (lowmemory=1), Policy [EZ cons-units]: should give zero: %2.8f \n',max(abs(Policy1c(:)-Policy1d(:))))
fprintf('EZ gamma=1/phi collapse with sj (lowmemory=1), V after transform (relative) [EZ cons-units]: should be roughly 1e-13: %g \n',max(abs(V1c(:)-V1dtransformed(:)))/max(abs(V1c(:))))
% the same on the lowmemory=2 rung (both sides of the collapse solved with lowmemory=2)
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
vfoptions1.survivalprobability='sj';
vfoptions1.lowmemory=2;
[V1c,Policy1c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1=rmfield(vfoptions1,'survivalprobability');
vfoptions1.exoticpreferences='None';
[V1d,Policy1d]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames2,[],vfoptions1);
V1dtransformed=((1-Params.ezgamma)*V1d).^(1/(1-Params.ezgamma));
fprintf('EZ gamma=1/phi collapse with sj (lowmemory=2), Policy [EZ cons-units]: should give zero: %2.8f \n',max(abs(Policy1c(:)-Policy1d(:))))
fprintf('EZ gamma=1/phi collapse with sj (lowmemory=2), V after transform (relative) [EZ cons-units]: should be roughly 1e-13: %g \n',max(abs(V1c(:)-V1dtransformed(:)))/max(abs(V1c(:))))

Params.ezphi=ezphi_store;
Params.ezsigma=ezsigma_store;
clear V1a V1b V1btransformed Policy1a Policy1b V1c V1d V1dtransformed Policy1c Policy1d

% (vi).3 utility-units EZriskaversion=0 WITH the declining sj (as (ii), plus survivalprobability='sj';
% the vNM oracle discounts by beta*sj). V and Policy exact, both signs of the utility fn.
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
vfoptions1=rmfield(vfoptions1,'survivalprobability');
vfoptions1.exoticpreferences='None';
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames2,[],vfoptions1);
fprintf('EZ with EZriskaversion=0 and sj [EZ positive utils]: should give zero: %2.8f \n',max(abs(V2a(:)-V2b(:))))
fprintf('EZ with EZriskaversion=0 and sj [EZ positive utils]: should give zero: %2.8f \n',max(abs(Policy2a(:)-Policy2b(:))))

% negative-valued utility fn
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.survivalprobability='sj';
[V3a,Policy3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1=rmfield(vfoptions1,'survivalprobability');
vfoptions1.exoticpreferences='None';
[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames2,[],vfoptions1);
fprintf('EZ with EZriskaversion=0 and sj [EZ negative utils]: should give zero: %2.8f \n',max(abs(V3a(:)-V3b(:))))
fprintf('EZ with EZriskaversion=0 and sj [EZ negative utils]: should give zero: %2.8f \n',max(abs(Policy3a(:)-Policy3b(:))))

Params.ezrisk=ezrisk_store;
clear V2a V2b V3a V3b Policy2a Policy2b Policy3a Policy3b

% (vi).4 cross-method agreement with the declining sj, all three EZ cases: ValueFnFromPolicy
% on the basic solve (the first exercise of the FromPolicy sj path), plus the file's lowmemory
% rungs (lowmemory=0, 1 and 2); noa1 riskyasset has no DC/GI methods to cross.
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits all the riskyasset settings (refine_d, aprimeFn, n_u, u_grid, pi_u) and the e-shock (n_e, e_grid, pi_e)
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
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsv);
    fprintf('sj, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    vfoptionsv.lowmemory=1;
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('sj, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('sj, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    vfoptionsv.lowmemory=0;
    vfoptionsv.lowmemory=2;
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('sj, lowmemory=2 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('sj, lowmemory=2 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    vfoptionsv.lowmemory=0;
end
clear V1 V2 Policy1 Policy2 V1fromPolicy

%% (vii) Warm-glow of bequests (De Nardi luxury-good form), vfoptions.WarmGlowBequestsFn
% Each EZ case gets a warm-glow fn matching its sign/units convention (see EZRisky_ReturnFns):
% cons-units gets EZRiskyWarmGlowFn_cons (a consumption-equivalent, strictly positive; curvature
% comes from the EZ preferences), positive/negative utils get EZRiskyWarmGlowFn_positiveUtils/
% EZRiskyWarmGlowFn_negativeUtils (strictly positive/strictly negative, nonzero at a2prime=0 so
% the WG==0 mask conventions are avoided). In riskyasset the toolkit evaluates the warm-glow at
% a2prime, the risky-asset holding realized via the (d,u) lottery, so the warm-glow goes through
% the same lottery as the continuation-value expectation.
WGFn_cons=@(aprime,wg1,wg2) EZRiskyWarmGlowFn_cons(aprime,wg1,wg2);
WGFn_posU=@(aprime,wg1,wg2,wg3) EZRiskyWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3);
WGFn_negU=@(aprime,wg1,wg2,wg3) EZRiskyWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);
Params2.sjterm=[ones(1,N_j-1),0]; % the sj the toolkit defaults to when WarmGlowBequestsFn is set but survivalprobability is not

% (vii).1 terminal-only warm-glow: WarmGlowBequestsFn set but NO survivalprobability, so the
% toolkit defaults to warm-glow only at the end of the final period (the V1 solve and the
% FromPolicy call each print the expected warning). All three EZ cases: FromPolicy, and the
% identity with explicitly passing survivalprobability='sjterm'.
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits all the riskyasset settings (refine_d, aprimeFn, n_u, u_grid, pi_u) and the e-shock (n_e, e_grid, pi_e)
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        vfoptionsv.WarmGlowBequestsFn=WGFn_cons;
        ReturnFn=ReturnFn_cons;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        vfoptionsv.WarmGlowBequestsFn=WGFn_posU;
        ReturnFn=ReturnFn_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        vfoptionsv.WarmGlowBequestsFn=WGFn_negU;
        ReturnFn=ReturnFn_negU;
    end
    vfoptions1=vfoptionsv;
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptions1);
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,vfoptions1);
    fprintf('warm-glow terminal-only, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    vfoptions5=vfoptions1;
    vfoptions5.survivalprobability='sjterm';
    [V5,Policy5]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptions5);
    fprintf('warm-glow terminal-only vs explicit sjterm [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V5(:))))
    fprintf('warm-glow terminal-only vs explicit sjterm [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy5(:))))
end
clear V1 V5 Policy1 Policy5 V1fromPolicy

% (vii).2 warm-glow at every age: WarmGlowBequestsFn plus the declining sj (warm-glow weight
% 1-sj(j) at every age, weight exactly one in the final period as sj(N_j)=0). All three EZ
% cases: FromPolicy, plus the file's lowmemory rungs (lowmemory=0, 1 and 2).
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits all the riskyasset settings (refine_d, aprimeFn, n_u, u_grid, pi_u) and the e-shock (n_e, e_grid, pi_e)
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        vfoptionsv.WarmGlowBequestsFn=WGFn_cons;
        ReturnFn=ReturnFn_cons;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        vfoptionsv.WarmGlowBequestsFn=WGFn_posU;
        ReturnFn=ReturnFn_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        vfoptionsv.WarmGlowBequestsFn=WGFn_negU;
        ReturnFn=ReturnFn_negU;
    end
    vfoptionsv.survivalprobability='sj';
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsv);
    fprintf('warm-glow with sj, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    vfoptionsv.lowmemory=1;
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('warm-glow with sj, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('warm-glow with sj, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    vfoptionsv.lowmemory=0;
    vfoptionsv.lowmemory=2;
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('warm-glow with sj, lowmemory=2 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('warm-glow with sj, lowmemory=2 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    vfoptionsv.lowmemory=0;
end
clear V1 V2 Policy1 Policy2 V1fromPolicy

% (vii).3 exact collapse oracles with sj AND warm-glow: DROPPED for riskyasset. The vNM
% riskyasset solvers have no warm-glow support, and a WG(aprime) term cannot be folded into a
% riskyasset ReturnFn either (aprime is stochastic via the (d,u) lottery, so it is not available
% to the return fn). Warm-glow exactness is instead anchored by the degenerateu cross-test
% bridge and the N_j-1 terminal warm-glow identity in (vii).5 below.

% (vii).4 V_Jplus1 mini-leg with sj and warm-glow on: as (v) but only one jstar. The
% age-dependent sj and oneminussj are trimmed to the shorter horizon too (the toolkit checks
% length(sj) against the N_j of the model being solved).
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits all the riskyasset settings (refine_d, aprimeFn, n_u, u_grid, pi_u) and the e-shock (n_e, e_grid, pi_e)
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        vfoptionsv.WarmGlowBequestsFn=WGFn_cons;
        ReturnFn=ReturnFn_cons;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        vfoptionsv.WarmGlowBequestsFn=WGFn_posU;
        ReturnFn=ReturnFn_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        vfoptionsv.WarmGlowBequestsFn=WGFn_negU;
        ReturnFn=ReturnFn_negU;
    end
    vfoptionsv.survivalprobability='sj';
    jstar=round(2*N_j/3);
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    Paramsjs.sj=Params.sj(1:Njs);
    Paramsjs.oneminussj=Params.oneminussj(1:Njs);
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    vfoptionsjs=vfoptionsv;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
    Vbase=Vbase(:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of the raws has its own lowmemory sub-branches)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i), lowmemory=1 [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i), lowmemory=1 [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;
    % lowmemory (the V_Jplus1 branch of the raws has its own lowmemory sub-branches)
    vfoptionsjs.lowmemory=2;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i), lowmemory=2 [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i), lowmemory=2 [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;
end
clear Vbase Vshort Policybase Policyshort

% (vii).5 the N_j-1 terminal warm-glow identity, ALL THREE EZ CASES. Solve the full N_j model
% with no warm-glow and no survival probability; then solve an N_j-1 model whose terminal-only
% warm-glow (weight exactly 1, since with no survivalprobability the default sj=[ones,0]
% applies; the expected warning prints) is the CLOSED FORM of the full model's terminal value
% fn: age N_j is retired, so the terminal budget is c=pension+a-savings (no (1+r)*a term:
% riskyasset returns are realised via aprimeFn, not in the budget), utility is increasing in c
% so the optimal terminal choice is savings=0 (d3_grid(1)=0) and h=0
% (retirement pays no wage, so full leisure; the composite good is then x=c^varphi), giving
% terminal V=u((pension+a)^varphi), which is shock-independent.
% V and Policy at ages 1,...,N_j-1 must then reproduce the full model's: the terminal V depends
% only on the risky asset, and the shorter model's warm-glow goes through the same (d,u)
% lottery on a2prime as the full model's certainty-equivalent over V_Nj(aprime).
% cons-units now included: the terminal warm-glow enters INSIDE the ^ezc7 root (Kraft/Munk/Weiss 2022), so its warm-glow fn is the terminal value fn itself, which for cons-units is just the composite good x ((ezc1*F^ezc2)^ezc7=F for ezc1=1).
WGterm_cons=@(aprime,varphi,pension) (aprime+pension)^varphi;
WGterm_posU=@(aprime,ezsigma,varphi,pension) ((1+(aprime+pension)^varphi)^(1-ezsigma)-1)/(1-ezsigma);
WGterm_negU=@(aprime,ezsigma,varphi,pension) (((aprime+pension)^varphi)^(1-ezsigma))/(1-ezsigma);
Njs2=N_j-1;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs2);
Paramsjs.kappa_j=Params.kappa_j(1:Njs2);
% consumption-units (traditional Epstein-Zin)
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
Vbase=Vbase(:,:,:,1:Njs2);
Policybase=Policybase(:,:,:,:,1:Njs2);
vfoptions1.WarmGlowBequestsFn=WGterm_cons;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs2,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Paramsjs,DiscountFactorParamNames,[],vfoptions1);
fprintf('N_j-1 terminal warm-glow identity [EZ cons-units], this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('N_j-1 terminal warm-glow identity [EZ cons-units], this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
% positive-valued utility fn
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=1;
vfoptions1.EZriskaversion='ezrisk';
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
Vbase=Vbase(:,:,:,1:Njs2);
Policybase=Policybase(:,:,:,:,1:Njs2);
vfoptions1.WarmGlowBequestsFn=WGterm_posU;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs2,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Paramsjs,DiscountFactorParamNames,[],vfoptions1);
fprintf('N_j-1 terminal warm-glow identity [EZ positive utils], this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('N_j-1 terminal warm-glow identity [EZ positive utils], this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
% negative-valued utility fn
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
Vbase=Vbase(:,:,:,1:Njs2);
Policybase=Policybase(:,:,:,:,1:Njs2);
vfoptions1.WarmGlowBequestsFn=WGterm_negU;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs2,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Paramsjs,DiscountFactorParamNames,[],vfoptions1);
fprintf('N_j-1 terminal warm-glow identity [EZ negative utils], this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('N_j-1 terminal warm-glow identity [EZ negative utils], this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
clear Vbase Vshort Policybase Policyshort

%% (viii) vfoptions.EZmortalityriskaversion (double Epstein-Zin: separate mortality risk aversion)
% Only meaningful together with survival probabilities, so survivalprobability='sj' throughout
% this section (with sj==1 the ^ezc8 exactly cancels the ^ezc6 in the survival aggregator and
% the option is vacuous).
% (viii).1 identity: setting EZmortalityriskaversion EQUAL to the case's own risk aversion
% (ezgamma for cons-units, ezrisk for utility-units) gives ezc8=1, which must reproduce the
% unset default exactly.
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits all the riskyasset settings (refine_d, aprimeFn, n_u, u_grid, pi_u) and the e-shock (n_e, e_grid, pi_e)
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
        Params2.ezmriskeq=Params.ezgamma; % the cons-units risk aversion
    else
        Params2.ezmriskeq=Params.ezrisk; % the utility-units risk aversion
    end
    [V1a,Policy1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptionsv);
    vfoptionsv.EZmortalityriskaversion='ezmriskeq';
    [V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('EZmortalityriskaversion equal to own risk aversion vs unset [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1a(:)-V1b(:))))
    fprintf('EZmortalityriskaversion equal to own risk aversion vs unset [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1a(:)-Policy1b(:))))
end
clear V1a V1b Policy1a Policy1b

% (viii).2 a DISTINCT mortality risk aversion (ezmrisk=5): the first nontrivial exercise of the
% ^ezc8 sites. All three EZ cases: FromPolicy, plus the file's lowmemory rungs (lowmemory=0, 1 and 2).
for ezcase=1:3
    vfoptionsv=vfoptions; % inherits all the riskyasset settings (refine_d, aprimeFn, n_u, u_grid, pi_u) and the e-shock (n_e, e_grid, pi_e)
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
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsv);
    fprintf('EZmortalityriskaversion distinct, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    vfoptionsv.lowmemory=1;
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('EZmortalityriskaversion distinct, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('EZmortalityriskaversion distinct, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    vfoptionsv.lowmemory=0;
    vfoptionsv.lowmemory=2;
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    fprintf('EZmortalityriskaversion distinct, lowmemory=2 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('EZmortalityriskaversion distinct, lowmemory=2 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    vfoptionsv.lowmemory=0;
end
clear V1 V2 Policy1 Policy2 V1fromPolicy

%% (ix) vfoptions.EZoneminusbeta=2
% The (1-sj*beta) generalization of (iv): EZoneminusbeta=2 puts a (1-sj*beta)* term on the this
% period return, an age-dependent factor since sj varies with age (survivalprobability='sj'
% throughout; without it =2 is the same as =1). As in (iv):
% - consumption-units: EZ is homogeneous of degree 1 in the return, so EZoneminusbeta=2 with return
%   x should equal EZoneminusbeta=0 with return scaled by (1-sj*beta)^(1/(1-1/ezphi));
% - utility-units: V=(1-sj*beta)*u+beta*(CE) is the same recursion as default with u scaled by (1-sj*beta).
% The age-dependent scale factors enter the scaled return fns as trailing (age-dependent) parameters,
% since GPU arrayfun does not support anonymous functions that capture workspace variables.
Params2.ezscalefactor2=(1-Params.sj*Params.beta).^(1/(1-1/Params.ezphi));
Params2.ezscalefactoru2=1-Params.sj*Params.beta;

% consumption-units
ReturnFn_cons_scaled2=@(h,savings,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension,ezscalefactor2) ezscalefactor2*EZRiskyReturnFn_cons_d1_z_e_nosemiz(h,savings,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension);
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
fprintf('EZoneminusbeta=2 vs manual (1-sj*beta) scaling [EZ cons-units]: should give zero: %2.8f \n',max(abs(V4a(:)-V4b(:))))
fprintf('EZoneminusbeta=2 vs manual (1-sj*beta) scaling [EZ cons-units]: should give zero: %2.8f \n',max(abs(Policy4a(:)-Policy4b(:))))

% utility-units, positive
ReturnFn_posU_scaled2=@(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension,ezscalefactoru2) ezscalefactoru2*EZRiskyReturnFn_positiveUtils_d1_z_e_nosemiz(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
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
fprintf('EZoneminusbeta=2 vs manual (1-sj*beta) scaling [EZ positive utils]: should give zero: %2.8f \n',max(abs(V5a(:)-V5b(:))))
fprintf('EZoneminusbeta=2 vs manual (1-sj*beta) scaling [EZ positive utils]: should give zero: %2.8f \n',max(abs(Policy5a(:)-Policy5b(:))))

% utility-units, negative
ReturnFn_negU_scaled2=@(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension,ezscalefactoru2) ezscalefactoru2*EZRiskyReturnFn_negativeUtils_d1_z_e_nosemiz(h,savings,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
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
fprintf('EZoneminusbeta=2 vs manual (1-sj*beta) scaling [EZ negative utils]: should give zero: %2.8f \n',max(abs(V6a(:)-V6b(:))))
fprintf('EZoneminusbeta=2 vs manual (1-sj*beta) scaling [EZ negative utils]: should give zero: %2.8f \n',max(abs(Policy6a(:)-Policy6b(:))))

clear V4a V4b V5a V5b V6a V6b Policy4a Policy4b Policy5a Policy5b Policy6a Policy6b

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
