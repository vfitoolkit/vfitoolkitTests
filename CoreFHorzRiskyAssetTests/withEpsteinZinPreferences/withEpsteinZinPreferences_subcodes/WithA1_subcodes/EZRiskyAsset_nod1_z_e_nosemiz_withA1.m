function output=EZRiskyAsset_nod1_z_e_nosemiz_withA1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Epstein-Zin mirror of CoreFHorzRiskyAsset_nod1_z_e_nosemiz_withA1.m. The model is run three times:
%   Case 1: consumption-units (traditional EZ, vfoptions.EZutils=0)
%   Case 2: utility-units, positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%   Case 3: utility-units, negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% followed by the EZ special tests (gamma=1/phi collapse; EZriskaversion=0 collapse; the same
% collapses under the grid interpolation layer; EZoneminusbeta=1 vs manual scaling). The
% u-shock is part of the EZ certainty-equivalent: ONE joint CE over (u, zprime, eprime) — the
% gamma=1/phi collapse only holds exactly if the implementation does this, so it is enforced
% by test.
% Plus special test (v): V_Jplus1 (V of period jstar as the terminal value fn of a shorter model).
% Plus special tests (vi)-(ix):
%   (vi)   survival probabilities: sjones plumbing; the (i)/(ii) collapse oracles with declining sj
%          (the vNM riskyasset oracle discounts by DiscountFactorParamNames={'beta','sj'}), each
%          also under the grid interpolation layer; cross-method agreement (basic/DC/GI/DC+GI,
%          FromPolicy) with sj
%   (vii)  warm-glow of bequests (De Nardi luxury-good form, evaluated at a2prime): terminal-only
%          default and with declining sj — the main exactness content is method agreement across
%          the four independent implementations plus FromPolicy; a V_Jplus1 mini-leg. The vNM
%          warm-glow oracle is dropped (no vNM riskyasset warm-glow support; see (vii).3) and the
%          N_j-1 terminal warm-glow identity is noa1-only (see (vii).5)
%   (viii) EZmortalityriskaversion: identity when set equal to the case's own risk aversion;
%          cross-method agreement with a distinct mortality risk aversion (ezmrisk=5)
%   (ix)   EZoneminusbeta=2 versus manually scaling the return fn by the age-dependent (1-sj*beta) factor
% With a1 (a=[a1,a2]: a1=safe asset, a2=risky asset) the riskyasset supports divide-and-conquer
% and the grid interpolation layer (unlike the noa1 case), so each EZ case runs the full
% basic/DC/GI/DC+GI solver ladder with lowmemory and ValueFnFromPolicy checks, then compares
% big-a_grid moments with/without GI.
% TEST-FIRST: the EZ riskyasset DC/GI/DC+GI tiers and the riskyasset branch of
% ValueFnFromPolicy_FHorz_EpsteinZin do not exist yet, so those legs are EXPECTED to error
% until implemented.
% n_d=[n_d2,n_d3] = [riskyshare, savings];   d_grid=[d2_grid; d3_grid]

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.pi_e=simoptionsbaseline.pi_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
% Riskyasset
vfoptions.riskyasset=1;
vfoptions.refine_d=[0,1,1]; % no d1, 1 d2 (riskyshare), 1 d3 (savings)
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
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

jequaloneDist=zeros([n_a_big,n_z,vfoptions.n_e],'gpuArray');
jequaloneDist(1,1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;

ReturnFn_cons=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,r_a1,agej,Jr,pension) EZRiskyReturnFn_cons_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,r_a1,agej,Jr,pension);
ReturnFn_posU=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension) EZRiskyReturnFn_positiveUtils_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension);
ReturnFn_negU=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension) EZRiskyReturnFn_negativeUtils_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension);

% With riskyasset, need to include all d in FnsToEvaluate, even if they are not in the ReturnFn
FnsToEvaluate.a1=@(riskyshare,savings,a1prime,a1,a2,z,e) a1;
FnsToEvaluate.a2=@(riskyshare,savings,a1prime,a1,a2,z,e) a2;
FnsToEvaluate.savings=@(riskyshare,savings,a1prime,a1,a2,z,e) savings;
FnsToEvaluate.earnings=@(riskyshare,savings,a1prime,a1,a2,z,e,w,kappa_j) w*kappa_j*z*e;


%% Case 1: consumption-units (traditional Epstein-Zin, vfoptions.EZutils=0)
vfoptions.exoticpreferences='EpsteinZin';
vfoptions.EZutils=0;
vfoptions.EZriskaversion='ezgamma';
vfoptions.EZeis='ezphi';
ReturnFn=ReturnFn_cons;

%% Basic VFI
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% V from Policy
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

%% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))

% lowmemory tests
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=2;
[V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(V1(:)-V1C(:))))
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=0;

% lowmemory tests
vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=2;
[V2C,Policy2C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=2 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2C(:))))
fprintf('lowmemory=2 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2C(:))))
vfoptions2.lowmemory=0;

clear V1 V2 V1B V2B V1C V2C Policy1 Policy2 Policy1B Policy2B Policy1C Policy2C V1fromPolicy

%% Solve with grid-interpolation (a valid, more-accurate solution)
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

% V from Policy under the grid interpolation layer
V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy with grid interp, this should be zero: %2.8f \n',max(abs(V3fromPolicy(:)-V3(:))))

%% Solve with divide-and-conquer + grid interp (should match grid-interp)
vfoptions4=vfoptions;
vfoptions4.divideandconquer=1;
vfoptions4.gridinterplayer=1;
vfoptions4.ngridinterp=5;
simoptions4=simoptions;
simoptions4.gridinterplayer=vfoptions4.gridinterplayer;
simoptions4.ngridinterp=vfoptions4.ngridinterp;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);

fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %2.8f \n',max(abs(V3(:)-V4(:))))
fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy4(:))))

% lowmemory tests
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=2;
[V3C,Policy3C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=2 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3C(:))))
fprintf('lowmemory=2 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3C(:))))
vfoptions3.lowmemory=0;

% lowmemory tests
vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=2;
[V4C,Policy4C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=2 (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4C(:))))
fprintf('lowmemory=2 (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4C(:))))
vfoptions4.lowmemory=0;

clear V3 V4 V3B V4B V3C V4C Policy3 Policy4 Policy3B Policy4B Policy3C Policy4C V3fromPolicy

%% Big a_grid: the moments should be essentially the same with/without grid interpolation
simoptions1.a_grid=a_grid_big;
[~,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

simoptions3.a_grid=a_grid_big;
[~,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy3b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions3);
AllStats3=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);
AgeConditionalStats3=LifeCycleProfiles_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('With/without grid interp, should get much the same moments (for big a_grid; cons-units) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.a1.Mean,AllStats3.a1.Mean] %#ok<NOPRT>
[AllStats1.a2.Mean,AllStats3.a2.Mean] %#ok<NOPRT>
[AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean] %#ok<NOPRT>
[AgeConditionalStats1.a1.StdDeviation; AgeConditionalStats3.a1.StdDeviation] %#ok<NOPRT>

clear Policy1b Policy3b StationaryDist1 StationaryDist3

% Also with divide-and-conquer
simoptions2.a_grid=a_grid_big;
[~,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy2b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AllStats2=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);
AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);

simoptions4.a_grid=a_grid_big;
[~,Policy4b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy4b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions4);
AllStats4=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);
AgeConditionalStats4=LifeCycleProfiles_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);

fprintf('With/without grid interp (with DC), should get much the same moments (for big a_grid; cons-units) \n')
[AllStats2.a1.Mean,AllStats4.a1.Mean] %#ok<NOPRT>
[AllStats2.a2.Mean,AllStats4.a2.Mean] %#ok<NOPRT>
[AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean] %#ok<NOPRT>
[AgeConditionalStats2.a1.StdDeviation; AgeConditionalStats4.a1.StdDeviation] %#ok<NOPRT>

clear Policy2b Policy4b StationaryDist2 StationaryDist4

%% Graph
figure(figure_c);
subplot(3,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('EZ cons-units: Earnings Mean'); legend('1','2','3','4')
subplot(3,1,2); plot(1:1:N_j,AgeConditionalStats1.a1.StdDeviation, 1:1:N_j,AgeConditionalStats2.a1.StdDeviation, 1:1:N_j,AgeConditionalStats3.a1.StdDeviation, 1:1:N_j,AgeConditionalStats4.a1.StdDeviation)
title('EZ cons-units: a1 (safe) Std Dev'); legend('1','2','3','4')
subplot(3,1,3); plot(1:1:N_j,AgeConditionalStats1.a2.Mean, 1:1:N_j,AgeConditionalStats2.a2.Mean, 1:1:N_j,AgeConditionalStats3.a2.Mean, 1:1:N_j,AgeConditionalStats4.a2.Mean)
title('EZ cons-units: a2 (risky) Mean'); legend('1','2','3','4')

clear AllStats1 AllStats2 AllStats3 AllStats4 AgeConditionalStats1 AgeConditionalStats2 AgeConditionalStats3 AgeConditionalStats4


%% Case 2: utility-units with POSITIVE-valued utility fn (vfoptions.EZutils=1, EZpositiveutility=1)
vfoptions.exoticpreferences='EpsteinZin';
vfoptions.EZutils=1;
vfoptions.EZpositiveutility=1;
vfoptions.EZriskaversion='ezrisk';
ReturnFn=ReturnFn_posU;

%% Basic VFI
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% V from Policy
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

%% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))

% lowmemory tests
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=2;
[V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(V1(:)-V1C(:))))
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=0;

% lowmemory tests
vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=2;
[V2C,Policy2C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=2 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2C(:))))
fprintf('lowmemory=2 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2C(:))))
vfoptions2.lowmemory=0;

clear V1 V2 V1B V2B V1C V2C Policy1 Policy2 Policy1B Policy2B Policy1C Policy2C V1fromPolicy

%% Solve with grid-interpolation (a valid, more-accurate solution)
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

% V from Policy under the grid interpolation layer
V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy with grid interp, this should be zero: %2.8f \n',max(abs(V3fromPolicy(:)-V3(:))))

%% Solve with divide-and-conquer + grid interp (should match grid-interp)
vfoptions4=vfoptions;
vfoptions4.divideandconquer=1;
vfoptions4.gridinterplayer=1;
vfoptions4.ngridinterp=5;
simoptions4=simoptions;
simoptions4.gridinterplayer=vfoptions4.gridinterplayer;
simoptions4.ngridinterp=vfoptions4.ngridinterp;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);

fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %2.8f \n',max(abs(V3(:)-V4(:))))
fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy4(:))))

% lowmemory tests
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=2;
[V3C,Policy3C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=2 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3C(:))))
fprintf('lowmemory=2 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3C(:))))
vfoptions3.lowmemory=0;

% lowmemory tests
vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=2;
[V4C,Policy4C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=2 (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4C(:))))
fprintf('lowmemory=2 (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4C(:))))
vfoptions4.lowmemory=0;

clear V3 V4 V3B V4B V3C V4C Policy3 Policy4 Policy3B Policy4B Policy3C Policy4C V3fromPolicy

%% Big a_grid: the moments should be essentially the same with/without grid interpolation
simoptions1.a_grid=a_grid_big;
[~,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

simoptions3.a_grid=a_grid_big;
[~,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy3b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions3);
AllStats3=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);
AgeConditionalStats3=LifeCycleProfiles_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('With/without grid interp, should get much the same moments (for big a_grid; positive utils) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.a1.Mean,AllStats3.a1.Mean] %#ok<NOPRT>
[AllStats1.a2.Mean,AllStats3.a2.Mean] %#ok<NOPRT>
[AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean] %#ok<NOPRT>
[AgeConditionalStats1.a1.StdDeviation; AgeConditionalStats3.a1.StdDeviation] %#ok<NOPRT>

clear Policy1b Policy3b StationaryDist1 StationaryDist3

% Also with divide-and-conquer
simoptions2.a_grid=a_grid_big;
[~,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy2b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AllStats2=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);
AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);

simoptions4.a_grid=a_grid_big;
[~,Policy4b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy4b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions4);
AllStats4=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);
AgeConditionalStats4=LifeCycleProfiles_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);

fprintf('With/without grid interp (with DC), should get much the same moments (for big a_grid; positive utils) \n')
[AllStats2.a1.Mean,AllStats4.a1.Mean] %#ok<NOPRT>
[AllStats2.a2.Mean,AllStats4.a2.Mean] %#ok<NOPRT>
[AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean] %#ok<NOPRT>
[AgeConditionalStats2.a1.StdDeviation; AgeConditionalStats4.a1.StdDeviation] %#ok<NOPRT>

clear Policy2b Policy4b StationaryDist2 StationaryDist4

%% Graph
figure(100+figure_c); % Case 2, positive-valued utils (Case 1 is figure_c)
subplot(3,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('EZ positive utils: Earnings Mean'); legend('1','2','3','4')
subplot(3,1,2); plot(1:1:N_j,AgeConditionalStats1.a1.StdDeviation, 1:1:N_j,AgeConditionalStats2.a1.StdDeviation, 1:1:N_j,AgeConditionalStats3.a1.StdDeviation, 1:1:N_j,AgeConditionalStats4.a1.StdDeviation)
title('EZ positive utils: a1 (safe) Std Dev'); legend('1','2','3','4')
subplot(3,1,3); plot(1:1:N_j,AgeConditionalStats1.a2.Mean, 1:1:N_j,AgeConditionalStats2.a2.Mean, 1:1:N_j,AgeConditionalStats3.a2.Mean, 1:1:N_j,AgeConditionalStats4.a2.Mean)
title('EZ positive utils: a2 (risky) Mean'); legend('1','2','3','4')

clear AllStats1 AllStats2 AllStats3 AllStats4 AgeConditionalStats1 AgeConditionalStats2 AgeConditionalStats3 AgeConditionalStats4


%% Case 3: utility-units with NEGATIVE-valued utility fn (vfoptions.EZutils=1, EZpositiveutility=0)
vfoptions.exoticpreferences='EpsteinZin';
vfoptions.EZutils=1;
vfoptions.EZpositiveutility=0;
vfoptions.EZriskaversion='ezrisk';
ReturnFn=ReturnFn_negU;

%% Basic VFI
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% V from Policy
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

%% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))

% lowmemory tests
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=2;
[V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(V1(:)-V1C(:))))
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=0;

% lowmemory tests
vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=2;
[V2C,Policy2C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=2 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2C(:))))
fprintf('lowmemory=2 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2C(:))))
vfoptions2.lowmemory=0;

clear V1 V2 V1B V2B V1C V2C Policy1 Policy2 Policy1B Policy2B Policy1C Policy2C V1fromPolicy

%% Solve with grid-interpolation (a valid, more-accurate solution)
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

% V from Policy under the grid interpolation layer
V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy with grid interp, this should be zero: %2.8f \n',max(abs(V3fromPolicy(:)-V3(:))))

%% Solve with divide-and-conquer + grid interp (should match grid-interp)
vfoptions4=vfoptions;
vfoptions4.divideandconquer=1;
vfoptions4.gridinterplayer=1;
vfoptions4.ngridinterp=5;
simoptions4=simoptions;
simoptions4.gridinterplayer=vfoptions4.gridinterplayer;
simoptions4.ngridinterp=vfoptions4.ngridinterp;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);

fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %2.8f \n',max(abs(V3(:)-V4(:))))
fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy4(:))))

% lowmemory tests
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=2;
[V3C,Policy3C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=2 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3C(:))))
fprintf('lowmemory=2 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3C(:))))
vfoptions3.lowmemory=0;

% lowmemory tests
vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=2;
[V4C,Policy4C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=2 (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4C(:))))
fprintf('lowmemory=2 (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4C(:))))
vfoptions4.lowmemory=0;

clear V3 V4 V3B V4B V3C V4C Policy3 Policy4 Policy3B Policy4B Policy3C Policy4C V3fromPolicy

%% Big a_grid: the moments should be essentially the same with/without grid interpolation
simoptions1.a_grid=a_grid_big;
[~,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

simoptions3.a_grid=a_grid_big;
[~,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy3b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions3);
AllStats3=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);
AgeConditionalStats3=LifeCycleProfiles_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('With/without grid interp, should get much the same moments (for big a_grid; negative utils) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.a1.Mean,AllStats3.a1.Mean] %#ok<NOPRT>
[AllStats1.a2.Mean,AllStats3.a2.Mean] %#ok<NOPRT>
[AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean] %#ok<NOPRT>
[AgeConditionalStats1.a1.StdDeviation; AgeConditionalStats3.a1.StdDeviation] %#ok<NOPRT>

clear Policy1b Policy3b StationaryDist1 StationaryDist3

% Also with divide-and-conquer
simoptions2.a_grid=a_grid_big;
[~,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy2b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AllStats2=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);
AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);

simoptions4.a_grid=a_grid_big;
[~,Policy4b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy4b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions4);
AllStats4=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);
AgeConditionalStats4=LifeCycleProfiles_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);

fprintf('With/without grid interp (with DC), should get much the same moments (for big a_grid; negative utils) \n')
[AllStats2.a1.Mean,AllStats4.a1.Mean] %#ok<NOPRT>
[AllStats2.a2.Mean,AllStats4.a2.Mean] %#ok<NOPRT>
[AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean] %#ok<NOPRT>
[AgeConditionalStats2.a1.StdDeviation; AgeConditionalStats4.a1.StdDeviation] %#ok<NOPRT>

clear Policy2b Policy4b StationaryDist2 StationaryDist4

%% Graph
figure(200+figure_c); % Case 3, negative-valued utils (Case 1 is figure_c)
subplot(3,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('EZ negative utils: Earnings Mean'); legend('1','2','3','4')
subplot(3,1,2); plot(1:1:N_j,AgeConditionalStats1.a1.StdDeviation, 1:1:N_j,AgeConditionalStats2.a1.StdDeviation, 1:1:N_j,AgeConditionalStats3.a1.StdDeviation, 1:1:N_j,AgeConditionalStats4.a1.StdDeviation)
title('EZ negative utils: a1 (safe) Std Dev'); legend('1','2','3','4')
subplot(3,1,3); plot(1:1:N_j,AgeConditionalStats1.a2.Mean, 1:1:N_j,AgeConditionalStats2.a2.Mean, 1:1:N_j,AgeConditionalStats3.a2.Mean, 1:1:N_j,AgeConditionalStats4.a2.Mean)
title('EZ negative utils: a2 (risky) Mean'); legend('1','2','3','4')

clear AllStats1 AllStats2 AllStats3 AllStats4 AgeConditionalStats1 AgeConditionalStats2 AgeConditionalStats3 AgeConditionalStats4


%% Special Epstein-Zin tests
% The u-shock makes these sharper than in CoreFHorzTests: the collapses only hold exactly if
% the u-expectation sits INSIDE the joint certainty-equivalent (one CE over (u, zprime, eprime)), so a
% wrong placement of the u-expectation fails these loudly.

%% (i) Consumption-units: gamma=1/phi collapses Epstein-Zin to standard vNM expected utility
% With vfoptions.EZutils=0 and EZriskaversion=1/EZeis, the EZ value fn is an increasing transform
% of the vNM problem with period utility (x^(1-ezgamma))/(1-ezgamma), which is exactly
% EZRiskyReturnFn_negativeUtils_nod1_z_e_nosemiz_withA1 with ezsigma=ezgamma. So: identical Policy
% (exactly), and V_EZ=((1-ezgamma)*V_vNM).^(1/(1-ezgamma)) up to roundoff in the transform
% (~1e-13 relative).
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

%% (iii) The same three collapse tests, under the grid interpolation layer (withA1 subcodes
% only: the noa1 riskyasset has no GI). The EZ GI layer linearly interpolates the transformed
% continuation EV — here one joint CE over (u, zprime, eprime) — and linear interpolation commutes with
% the affine collapse transform, so these should be exact just like the no-GI versions.
vfoptionsGI=vfoptions;
vfoptionsGI.gridinterplayer=1;
vfoptionsGI.ngridinterp=5;

ezphi_store=Params.ezphi;
ezsigma_store=Params.ezsigma;
Params.ezphi=1/Params.ezgamma;
Params.ezsigma=Params.ezgamma;
vfoptions1=vfoptionsGI;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
[V1a,Policy1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.exoticpreferences='None';
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
V1btransformed=((1-Params.ezgamma)*V1b).^(1/(1-Params.ezgamma));
fprintf('EZ gamma=1/phi collapse (GI), Policy: should give zero: %2.8f \n',max(abs(Policy1a(:)-Policy1b(:))))
fprintf('EZ gamma=1/phi collapse (GI), V after transform (relative): should be roughly 1e-13: %g \n',max(abs(V1a(:)-V1btransformed(:)))/max(abs(V1a(:))))
Params.ezphi=ezphi_store;
Params.ezsigma=ezsigma_store;

ezrisk_store=Params.ezrisk;
Params.ezrisk=0;
vfoptions1=vfoptionsGI;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=1;
vfoptions1.EZriskaversion='ezrisk';
[V2a,Policy2a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.exoticpreferences='None';
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZ (positive utils) with EZriskaversion=0 (GI): should give zero: %2.8f \n',max(abs(V2a(:)-V2b(:))))
fprintf('EZ (positive utils) with EZriskaversion=0 (GI): should give zero: %2.8f \n',max(abs(Policy2a(:)-Policy2b(:))))
vfoptions1=vfoptionsGI;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
[V3a,Policy3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.exoticpreferences='None';
[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZ (negative utils) with EZriskaversion=0 (GI): should give zero: %2.8f \n',max(abs(V3a(:)-V3b(:))))
fprintf('EZ (negative utils) with EZriskaversion=0 (GI): should give zero: %2.8f \n',max(abs(Policy3a(:)-Policy3b(:))))
Params.ezrisk=ezrisk_store;
clear V1a V1b V1btransformed V2a V2b V3a V3b Policy1a Policy1b Policy2a Policy2b Policy3a Policy3b

%% (iv) vfoptions.EZoneminusbeta=1 versus manually scaling the return fn
% - consumption-units: EZ is homogeneous of degree 1 in the return, so EZoneminusbeta=1 with
%   return x should equal EZoneminusbeta=0 with return scaled by (1-beta)^(1/(1-1/ezphi)).
% - utility-units: V=(1-beta)*u+beta*(CE) is the same recursion as default with u scaled by (1-beta).
% (the scale factors enter the scaled return fns as parameters, since GPU arrayfun does
% not support anonymous functions that capture workspace variables)
Params.ezscalefactor=(1-Params.beta)^(1/(1-1/Params.ezphi));
Params.ezscalefactoru=1-Params.beta;

% consumption-units
ReturnFn_cons_scaled=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,r_a1,agej,Jr,pension,ezscalefactor) ezscalefactor*EZRiskyReturnFn_cons_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,r_a1,agej,Jr,pension);
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
ReturnFn_posU_scaled=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,ezscalefactoru) ezscalefactoru*EZRiskyReturnFn_positiveUtils_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension);
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
ReturnFn_negU_scaled=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,ezscalefactoru) ezscalefactoru*EZRiskyReturnFn_negativeUtils_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension);
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
% vfoptions.V_Jplus1=V(:,:,:,:,jstar). V_Jplus1 is the value fn of period N_j+1 of the model being
% solved, so the shorter model has Njs=jstar-1 periods, and the age-dependent parameters are
% trimmed to length Njs. V and Policy must then be identical to the original model for periods
% 1,...,jstar-1. Each of the four solution methods gets a different jstar (the last of them uses
% jstar=N_j, so that one of them covers the retirement periods). Run for all three EZ cases.
% Note: mewj is age-dependent, but is only used for the agent distribution, which is not
% computed here, so it is left alone.
for ezcase=1:3
    vfoptionsv=vfoptions; % baseline vfoptions: the riskyasset settings and n_e/e_grid/pi_e
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
    % The four solution methods, composed with this case's EZ vfoptions
    vfoptions1=vfoptionsv;
    vfoptions2=vfoptionsv;
    vfoptions2.divideandconquer=1;
    vfoptions3=vfoptionsv;
    vfoptions3.gridinterplayer=1;
    vfoptions3.ngridinterp=5;
    vfoptions4=vfoptionsv;
    vfoptions4.divideandconquer=1;
    vfoptions4.gridinterplayer=1;
    vfoptions4.ngridinterp=5;

    %% V_Jplus1, without divide-and-conquer, without grid interpolation
    jstar=round(3*N_j/4);
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    vfoptionsjs=vfoptions1;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,:,jstar);
    Vbase=Vbase(:,:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=2;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=2 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=2 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;

    %% V_Jplus1, with divide-and-conquer
    jstar=round(2*N_j/3);
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    vfoptionsjs=vfoptions2;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,:,jstar);
    Vbase=Vbase(:,:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1 (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=2;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=2 (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=2 (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;

    %% V_Jplus1, with grid interpolation
    jstar=round(N_j/2);
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    vfoptionsjs=vfoptions3;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,:,jstar);
    Vbase=Vbase(:,:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i, with GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i, with GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1 (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=2;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=2 (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=2 (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;

    %% V_Jplus1, with divide-and-conquer and grid interpolation
    jstar=N_j;
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    vfoptionsjs=vfoptions4;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,:,jstar);
    Vbase=Vbase(:,:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i, with DC+GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i, with DC+GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=2;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=2 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=2 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;
end

clear Vbase Vshort Policybase Policyshort

%% (vi) Survival probabilities (vfoptions.survivalprobability)
% Params.sj is declining with sj(N_j)=0 (set in CoreFHorzRiskyAssetEZTests.m). When
% vfoptions.survivalprobability is not set the EZ codes use sj=ones(N_j,1) internally, so
% survivalprobability='sjones' (all ones) must reproduce the baseline solve exactly (pure
% plumbing). Then the (i)/(ii) collapse oracles are repeated WITH the declining sj: once EZ
% collapses to vNM, the survival probability is just an age-dependent discount factor, so the
% vNM riskyasset oracle uses DiscountFactorParamNames={'beta','sj'}. The withA1 riskyasset has
% GI, so both collapses are also run under the grid interpolation layer (mirroring (iii)).
Params2=Params;
Params2.sjones=ones(1,N_j);
DiscountFactorParamNames2={'beta','sj'}; % beta*sj: age-dependent discounting for the vNM oracles

% (vi).1 plumbing: survivalprobability='sjones' vs not setting survivalprobability, all three EZ cases
for ezcase=1:3
    vfoptionsv=vfoptions; % baseline vfoptions: the riskyasset settings and n_e/e_grid/pi_e
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
% the vNM riskyasset oracle discounts by beta*sj). Policy exact; V via the usual transform. Also
% repeated under the grid interpolation layer (the EZ GI layer linearly interpolates the transformed
% continuation EV, and linear interpolation commutes with the affine collapse transform, so it
% should stay exact).
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
% the same under the grid interpolation layer
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
vfoptions1.survivalprobability='sj';
vfoptions1.gridinterplayer=1;
vfoptions1.ngridinterp=5;
[V1c,Policy1c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1=rmfield(vfoptions1,'survivalprobability');
vfoptions1.exoticpreferences='None';
[V1d,Policy1d]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames2,[],vfoptions1);
V1dtransformed=((1-Params.ezgamma)*V1d).^(1/(1-Params.ezgamma));
fprintf('EZ gamma=1/phi collapse with sj (GI), Policy [EZ cons-units]: should give zero: %2.8f \n',max(abs(Policy1c(:)-Policy1d(:))))
fprintf('EZ gamma=1/phi collapse with sj (GI), V after transform (relative) [EZ cons-units]: should be roughly 1e-13: %g \n',max(abs(V1c(:)-V1dtransformed(:)))/max(abs(V1c(:))))

Params.ezphi=ezphi_store;
Params.ezsigma=ezsigma_store;
clear V1a V1b V1btransformed V1c V1d V1dtransformed Policy1a Policy1b Policy1c Policy1d

% (vi).3 utility-units EZriskaversion=0 WITH the declining sj (as (ii), plus survivalprobability='sj';
% the vNM riskyasset oracle discounts by beta*sj). V and Policy exact, both signs of the utility fn;
% each also under the grid interpolation layer.
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
% positive-valued utility fn, under the grid interpolation layer
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=1;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.survivalprobability='sj';
vfoptions1.gridinterplayer=1;
vfoptions1.ngridinterp=5;
[V2c,Policy2c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1=rmfield(vfoptions1,'survivalprobability');
vfoptions1.exoticpreferences='None';
[V2d,Policy2d]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames2,[],vfoptions1);
fprintf('EZ with EZriskaversion=0 and sj (GI) [EZ positive utils]: should give zero: %2.8f \n',max(abs(V2c(:)-V2d(:))))
fprintf('EZ with EZriskaversion=0 and sj (GI) [EZ positive utils]: should give zero: %2.8f \n',max(abs(Policy2c(:)-Policy2d(:))))

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
% negative-valued utility fn, under the grid interpolation layer
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.survivalprobability='sj';
vfoptions1.gridinterplayer=1;
vfoptions1.ngridinterp=5;
[V3c,Policy3c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1=rmfield(vfoptions1,'survivalprobability');
vfoptions1.exoticpreferences='None';
[V3d,Policy3d]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames2,[],vfoptions1);
fprintf('EZ with EZriskaversion=0 and sj (GI) [EZ negative utils]: should give zero: %2.8f \n',max(abs(V3c(:)-V3d(:))))
fprintf('EZ with EZriskaversion=0 and sj (GI) [EZ negative utils]: should give zero: %2.8f \n',max(abs(Policy3c(:)-Policy3d(:))))

Params.ezrisk=ezrisk_store;
clear V2a V2b V2c V2d V3a V3b V3c V3d Policy2a Policy2b Policy2c Policy2d Policy3a Policy3b Policy3c Policy3d

% (vi).4 cross-method agreement with the declining sj, all three EZ cases: basic==DC, GI==DC+GI,
% and ValueFnFromPolicy on the basic and GI solves (the first exercise of the FromPolicy sj path).
for ezcase=1:3
    vfoptionsv=vfoptions; % baseline vfoptions: the riskyasset settings and n_e/e_grid/pi_e
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
    vfoptions1=vfoptionsv;
    vfoptions2=vfoptionsv;
    vfoptions2.divideandconquer=1;
    vfoptions3=vfoptionsv;
    vfoptions3.gridinterplayer=1;
    vfoptions3.ngridinterp=5;
    vfoptions4=vfoptions3;
    vfoptions4.divideandconquer=1;
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    [V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    [V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
    V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
    fprintf('sj, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    fprintf('sj, ValueFnFromPolicy (GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3fromPolicy(:)-V3(:))))
    fprintf('sj, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('sj, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    fprintf('sj, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V4(:))))
    fprintf('sj, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy4(:))))
end
clear V1 V2 V3 V4 Policy1 Policy2 Policy3 Policy4 V1fromPolicy V3fromPolicy

%% (vii) Warm-glow of bequests (De Nardi luxury-good form), vfoptions.WarmGlowBequestsFn
% Each EZ case gets a warm-glow fn matching its sign/units convention (see EZRisky_ReturnFns):
% cons-units gets EZRiskyWarmGlowFn_cons (a consumption-equivalent, strictly positive; curvature
% comes from the EZ preferences), positive/negative utils get
% EZRiskyWarmGlowFn_positiveUtils/EZRiskyWarmGlowFn_negativeUtils (strictly positive/strictly
% negative, nonzero at zero bequest so the WG==0 mask conventions are avoided). In the riskyasset
% banks the warm-glow argument is a2prime (the risky-asset holding realized via the (d,u) lottery;
% the toolkit evaluates the warm-glow on the a2 grid, a1prime-independent).
WGFn_cons=@(aprime,wg1,wg2) EZRiskyWarmGlowFn_cons(aprime,wg1,wg2);
WGFn_posU=@(aprime,wg1,wg2,wg3) EZRiskyWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3);
WGFn_negU=@(aprime,wg1,wg2,wg3) EZRiskyWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);
Params2.sjterm=[ones(1,N_j-1),0]; % the sj the toolkit defaults to when WarmGlowBequestsFn is set but survivalprobability is not

% (vii).1 terminal-only warm-glow: WarmGlowBequestsFn set but NO survivalprobability, so the toolkit
% defaults to warm-glow only at the end of the final period (each solve prints the expected warning).
% This is the main exactness content for withA1: four independent implementations (basic/DC/GI/DC+GI)
% must agree with the warm-glow active. All three EZ cases: FromPolicy (basic and GI), basic==DC,
% GI==DC+GI, and the identity with explicitly passing survivalprobability='sjterm'.
for ezcase=1:3
    vfoptionsv=vfoptions; % baseline vfoptions: the riskyasset settings and n_e/e_grid/pi_e
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
    vfoptions2=vfoptionsv;
    vfoptions2.divideandconquer=1;
    vfoptions3=vfoptionsv;
    vfoptions3.gridinterplayer=1;
    vfoptions3.ngridinterp=5;
    vfoptions4=vfoptions3;
    vfoptions4.divideandconquer=1;
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptions1);
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptions2);
    [V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptions3);
    [V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptions4);
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,vfoptions1);
    V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,vfoptions3);
    fprintf('warm-glow terminal-only, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    fprintf('warm-glow terminal-only, ValueFnFromPolicy (GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3fromPolicy(:)-V3(:))))
    fprintf('warm-glow terminal-only, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('warm-glow terminal-only, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    fprintf('warm-glow terminal-only, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V4(:))))
    fprintf('warm-glow terminal-only, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy4(:))))
    vfoptions5=vfoptions1;
    vfoptions5.survivalprobability='sjterm';
    [V5,Policy5]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptions5);
    fprintf('warm-glow terminal-only vs explicit sjterm [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V5(:))))
    fprintf('warm-glow terminal-only vs explicit sjterm [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy5(:))))
end
clear V1 V2 V3 V4 V5 Policy1 Policy2 Policy3 Policy4 Policy5 V1fromPolicy V3fromPolicy

% (vii).2 warm-glow at every age: WarmGlowBequestsFn plus the declining sj (warm-glow weight 1-sj(j)
% at every age, weight exactly one in the final period as sj(N_j)=0). All three EZ cases: FromPolicy
% (basic and GI), basic==DC, GI==DC+GI, and the file's lowmemory rungs.
for ezcase=1:3
    vfoptionsv=vfoptions; % baseline vfoptions: the riskyasset settings and n_e/e_grid/pi_e
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
    vfoptions1=vfoptionsv;
    vfoptions2=vfoptionsv;
    vfoptions2.divideandconquer=1;
    vfoptions3=vfoptionsv;
    vfoptions3.gridinterplayer=1;
    vfoptions3.ngridinterp=5;
    vfoptions4=vfoptions3;
    vfoptions4.divideandconquer=1;
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    [V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    [V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
    V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
    fprintf('warm-glow with sj, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    fprintf('warm-glow with sj, ValueFnFromPolicy (GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3fromPolicy(:)-V3(:))))
    fprintf('warm-glow with sj, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('warm-glow with sj, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    fprintf('warm-glow with sj, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V4(:))))
    fprintf('warm-glow with sj, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy4(:))))
    % lowmemory rungs, mirroring the file's tier legs
    vfoptions1.lowmemory=1;
    [V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    fprintf('warm-glow with sj, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V1B(:))))
    fprintf('warm-glow with sj, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy1B(:))))
    vfoptions1.lowmemory=2;
    [V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    fprintf('warm-glow with sj, lowmemory=2 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V1C(:))))
    fprintf('warm-glow with sj, lowmemory=2 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy1C(:))))
    vfoptions1.lowmemory=0;
    vfoptions2.lowmemory=1;
    [V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    fprintf('warm-glow with sj, lowmemory=1 (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V2(:)-V2B(:))))
    fprintf('warm-glow with sj, lowmemory=1 (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy2(:)-Policy2B(:))))
    vfoptions2.lowmemory=2;
    [V2C,Policy2C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    fprintf('warm-glow with sj, lowmemory=2 (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V2(:)-V2C(:))))
    fprintf('warm-glow with sj, lowmemory=2 (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy2(:)-Policy2C(:))))
    vfoptions2.lowmemory=0;
    vfoptions3.lowmemory=1;
    [V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    fprintf('warm-glow with sj, lowmemory=1 (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V3B(:))))
    fprintf('warm-glow with sj, lowmemory=1 (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy3B(:))))
    vfoptions3.lowmemory=2;
    [V3C,Policy3C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    fprintf('warm-glow with sj, lowmemory=2 (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V3C(:))))
    fprintf('warm-glow with sj, lowmemory=2 (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy3C(:))))
    vfoptions3.lowmemory=0;
    vfoptions4.lowmemory=1;
    [V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    fprintf('warm-glow with sj, lowmemory=1 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V4(:)-V4B(:))))
    fprintf('warm-glow with sj, lowmemory=1 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy4(:)-Policy4B(:))))
    vfoptions4.lowmemory=2;
    [V4C,Policy4C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    fprintf('warm-glow with sj, lowmemory=2 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V4(:)-V4C(:))))
    fprintf('warm-glow with sj, lowmemory=2 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy4(:)-Policy4C(:))))
    vfoptions4.lowmemory=0;
end
clear V1 V2 V3 V4 Policy1 Policy2 Policy3 Policy4 V1fromPolicy V3fromPolicy V1B V2B V3B V4B V1C V2C V3C V4C Policy1B Policy2B Policy3B Policy4B Policy1C Policy2C Policy3C Policy4C

% (vii).3 exact vNM collapse oracles with sj AND warm-glow: DROPPED for the riskyasset banks.
% In the main EZ bank the utility-units oracle adds beta*oneminussj*WG(aprime) to the return fn,
% but with a riskyasset aprime (here a2prime) is stochastic — realized via the (d,u) lottery after
% the choice is made — so WG(a2prime) cannot enter a riskyasset ReturnFn (which only sees the
% choices, not the realized a2prime). There is no vNM riskyasset warm-glow support to collapse to.
% Exactness of the warm-glow riskyasset code is instead anchored by the degenerateu cross-test
% bridge (EZRiskyAsset_CrossTests_degenerateu_*: riskyasset warm-glow vs plain-EZ warm-glow when
% the u-lottery is degenerate), plus the method-agreement checks in (vii).1/(vii).2 above.

% (vii).4 V_Jplus1 mini-leg with sj and warm-glow on: as (v) but only one jstar and only the basic
% and DC methods. The age-dependent sj and oneminussj are trimmed to the shorter horizon too (the
% toolkit checks length(sj) against the N_j of the model being solved).
for ezcase=1:3
    vfoptionsv=vfoptions; % baseline vfoptions: the riskyasset settings and n_e/e_grid/pi_e
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
    vfoptions1=vfoptionsv;
    vfoptions2=vfoptionsv;
    vfoptions2.divideandconquer=1;
    jstar=round(2*N_j/3);
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    Paramsjs.sj=Params.sj(1:Njs);
    Paramsjs.oneminussj=Params.oneminussj(1:Njs);
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    vfoptionsjs=vfoptions1;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,:,jstar);
    Vbase=Vbase(:,:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % with divide-and-conquer
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    vfoptionsjs=vfoptions2;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,:,jstar);
    Vbase=Vbase(:,:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
end
clear Vbase Vshort Policybase Policyshort

% (vii).5 the N_j-1 terminal warm-glow identity: EXCLUDED for withA1 — the terminal value fn depends
% on a1 as well as a2, which the a2prime-only warm-glow cannot represent, so it is a noa1-only test.

%% (viii) vfoptions.EZmortalityriskaversion (double Epstein-Zin: separate mortality risk aversion)
% Only meaningful together with survival probabilities, so survivalprobability='sj' throughout this
% section (with sj==1 the ^ezc8 exactly cancels the ^ezc6 in the survival aggregator and the option
% is vacuous).
% (viii).1 identity: setting EZmortalityriskaversion EQUAL to the case's own risk aversion (ezgamma
% for cons-units, ezrisk for utility-units) gives ezc8=1, which must reproduce the unset default
% exactly.
for ezcase=1:3
    vfoptionsv=vfoptions; % baseline vfoptions: the riskyasset settings and n_e/e_grid/pi_e
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
% ^ezc8 sites. All three EZ cases: FromPolicy (basic and GI), basic==DC, GI==DC+GI.
for ezcase=1:3
    vfoptionsv=vfoptions; % baseline vfoptions: the riskyasset settings and n_e/e_grid/pi_e
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
    vfoptions1=vfoptionsv;
    vfoptions2=vfoptionsv;
    vfoptions2.divideandconquer=1;
    vfoptions3=vfoptionsv;
    vfoptions3.gridinterplayer=1;
    vfoptions3.ngridinterp=5;
    vfoptions4=vfoptions3;
    vfoptions4.divideandconquer=1;
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    [V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    [V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
    V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
    fprintf('EZmortalityriskaversion distinct, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    fprintf('EZmortalityriskaversion distinct, ValueFnFromPolicy (GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3fromPolicy(:)-V3(:))))
    fprintf('EZmortalityriskaversion distinct, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('EZmortalityriskaversion distinct, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    fprintf('EZmortalityriskaversion distinct, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V4(:))))
    fprintf('EZmortalityriskaversion distinct, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy4(:))))
end
clear V1 V2 V3 V4 Policy1 Policy2 Policy3 Policy4 V1fromPolicy V3fromPolicy

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
ReturnFn_cons_scaled2=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,r_a1,agej,Jr,pension,ezscalefactor2) ezscalefactor2*EZRiskyReturnFn_cons_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,r_a1,agej,Jr,pension);
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
ReturnFn_posU_scaled2=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,ezscalefactoru2) ezscalefactoru2*EZRiskyReturnFn_positiveUtils_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension);
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
ReturnFn_negU_scaled2=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,ezscalefactoru2) ezscalefactoru2*EZRiskyReturnFn_negativeUtils_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension);
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
