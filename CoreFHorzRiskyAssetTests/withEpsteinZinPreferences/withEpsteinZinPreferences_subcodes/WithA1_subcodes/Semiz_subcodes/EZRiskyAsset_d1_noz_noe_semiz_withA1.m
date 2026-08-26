function output=EZRiskyAsset_d1_noz_noe_semiz_withA1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c) %#ok<INUSD>
% Epstein-Zin mirror of CoreFHorzRiskyAsset_d1_noz_noe_semiz_withA1.m (riskyasset with a1+a2:
% a1=standard safe asset earning r_a1, a2=risky asset). The model is run three times:
%   Case 1: consumption-units (traditional EZ, vfoptions.EZutils=0)
%   Case 2: utility-units, positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%   Case 3: utility-units, negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% Each case runs the full withA1 solver ladder: basic + DC + GI + DC+GI, each with its
% lowmemory legs, plus ValueFnFromPolicy checks (basic and under GI), then big-a_grid
% moments and a figure. Then the EZ special tests:
%   (i)   consumption-units: gamma=1/phi collapses EZ to standard vNM
%   (ii)  utility-units: EZriskaversion=0 collapses EZ to standard vNM
%   (iii) the same collapse tests, under the grid interpolation layer
%   (iv)  vfoptions.EZoneminusbeta=1 versus manually scaling the return fn
%   (v)   V_Jplus1: V of period jstar as the terminal value fn of a shorter (jstar-1)-period model
%   (vi)  survival probabilities (vfoptions.survivalprobability): the collapse tests with sj + FromPolicy
%   (vii) warm-glow of bequests (vfoptions.WarmGlowBequestsFn; the riskyasset warm-glow fns are
%         a2prime-only): terminal-only and with declining sj, plus a V_Jplus1 mini-leg
%   (viii) vfoptions.EZmortalityriskaversion (separate mortality risk aversion)
%   (ix)  vfoptions.EZoneminusbeta=2 versus manually scaling the return fn by (1-sj*beta)
% The u-shock is part of the EZ certainty-equivalent: ONE joint CE over
% (u, semizprime) — the collapse tests only hold exactly if the implementation
% does this, so it is enforced by test.
% Note: this noz_noe shape is NOT a no-shock model (the u shock and semiz provide genuine
% risk), so no 'EZ without shocks' warning is expected.
% TEST-FIRST: the EZ riskyasset withA1+semiz raws and the EZ riskyasset DC/GI tiers do not
% exist yet, so this file is EXPECTED TO ERROR until they are implemented (see
% EZRiskyAsset_coverage_proposal.md).
% n_d=[n_d1,n_d2,n_d3,n_d4] = [labor h, riskyshare, savings, semiz decision];   n_a=[n_a1,n_a2]

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
n_z=0;
z_grid=[];
pi_z=[];
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;

% Riskyasset
vfoptions.riskyasset=1;
vfoptions.refine_d=[1,1,1,1]; % 1 d1 (labor h), 1 d2 (riskyshare), 1 d3 (savings), 1 d4 (semiz decision)
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

jequaloneDist=zeros([n_a_big,vfoptions.n_semiz],'gpuArray');
jequaloneDist(1,1,ceil(vfoptions.n_semiz/2))=1;

ReturnFn_cons=@(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_cons_d1_noz_noe_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_posU=@(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_positiveUtils_d1_noz_noe_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_negU=@(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_negativeUtils_d1_noz_noe_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);

% With riskyasset, need to include all d in FnsToEvaluate, even if they are not in the ReturnFn
FnsToEvaluate.a1=@(h,riskyshare,savings,dsemiz,a1prime,a1,a2,semiz) a1;
FnsToEvaluate.a2=@(h,riskyshare,savings,dsemiz,a1prime,a1,a2,semiz) a2;
FnsToEvaluate.savings=@(h,riskyshare,savings,dsemiz,a1prime,a1,a2,semiz) savings;
FnsToEvaluate.earnings=@(h,riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,w,kappa_j) w*kappa_j*h*semiz;


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
vfoptions1.lowmemory=0;

% lowmemory tests
vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=0;

%%
clear V1 V2 V1B V2B Policy1 Policy2 Policy1B Policy2B V1fromPolicy
%% Solve with grid-interpolation (a valid, more-accurate solution)
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

% V from Policy
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
vfoptions3.lowmemory=0;

% lowmemory tests
vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=0;

%%
clear V3 V4 V3B V4B Policy3 Policy4 Policy3B Policy4B V3fromPolicy
%% Big a_grid for accurate StationaryDist moments
fprintf('Big a_grid moments (EZ cons-units)\n')
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

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
disp([AllStats1.a1.Mean,AllStats3.a1.Mean])
disp([AllStats1.a2.Mean,AllStats3.a2.Mean])
disp([AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean])
disp([AgeConditionalStats1.a1.StdDeviation; AgeConditionalStats3.a1.StdDeviation])

clear Policy1b Policy3b StationaryDist1 StationaryDist3

% Also with DC
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

fprintf('With/without grid interp (with DC), should get much the same moments \n')
disp([AllStats2.a1.Mean,AllStats4.a1.Mean])
disp([AllStats2.a2.Mean,AllStats4.a2.Mean])
disp([AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean])
disp([AgeConditionalStats2.a1.StdDeviation; AgeConditionalStats4.a1.StdDeviation])

clear Policy2b Policy4b StationaryDist2 StationaryDist4

%% Age-conditional plots
figure(figure_c); % Case 1, consumption-units
subplot(3,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('Earnings Mean (EZ cons-units)'); legend('1','2','3','4')
subplot(3,1,2); plot(1:1:N_j,AgeConditionalStats1.a1.StdDeviation, 1:1:N_j,AgeConditionalStats2.a1.StdDeviation, 1:1:N_j,AgeConditionalStats3.a1.StdDeviation, 1:1:N_j,AgeConditionalStats4.a1.StdDeviation)
title('a1 (safe) Std Dev'); legend('1','2','3','4')
subplot(3,1,3); plot(1:1:N_j,AgeConditionalStats1.a2.Mean, 1:1:N_j,AgeConditionalStats2.a2.Mean, 1:1:N_j,AgeConditionalStats3.a2.Mean, 1:1:N_j,AgeConditionalStats4.a2.Mean)
title('a2 (risky) Mean'); legend('1','2','3','4')

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
vfoptions1.lowmemory=0;

% lowmemory tests
vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=0;

%%
clear V1 V2 V1B V2B Policy1 Policy2 Policy1B Policy2B V1fromPolicy
%% Solve with grid-interpolation (a valid, more-accurate solution)
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

% V from Policy
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
vfoptions3.lowmemory=0;

% lowmemory tests
vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=0;

%%
clear V3 V4 V3B V4B Policy3 Policy4 Policy3B Policy4B V3fromPolicy
%% Big a_grid for accurate StationaryDist moments
fprintf('Big a_grid moments (EZ positive utils)\n')
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

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
disp([AllStats1.a1.Mean,AllStats3.a1.Mean])
disp([AllStats1.a2.Mean,AllStats3.a2.Mean])
disp([AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean])
disp([AgeConditionalStats1.a1.StdDeviation; AgeConditionalStats3.a1.StdDeviation])

clear Policy1b Policy3b StationaryDist1 StationaryDist3

% Also with DC
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

fprintf('With/without grid interp (with DC), should get much the same moments \n')
disp([AllStats2.a1.Mean,AllStats4.a1.Mean])
disp([AllStats2.a2.Mean,AllStats4.a2.Mean])
disp([AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean])
disp([AgeConditionalStats2.a1.StdDeviation; AgeConditionalStats4.a1.StdDeviation])

clear Policy2b Policy4b StationaryDist2 StationaryDist4

%% Age-conditional plots
figure(100+figure_c); % Case 2, positive-valued utils (Case 1 is figure_c)
subplot(3,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('Earnings Mean (EZ positive utils)'); legend('1','2','3','4')
subplot(3,1,2); plot(1:1:N_j,AgeConditionalStats1.a1.StdDeviation, 1:1:N_j,AgeConditionalStats2.a1.StdDeviation, 1:1:N_j,AgeConditionalStats3.a1.StdDeviation, 1:1:N_j,AgeConditionalStats4.a1.StdDeviation)
title('a1 (safe) Std Dev'); legend('1','2','3','4')
subplot(3,1,3); plot(1:1:N_j,AgeConditionalStats1.a2.Mean, 1:1:N_j,AgeConditionalStats2.a2.Mean, 1:1:N_j,AgeConditionalStats3.a2.Mean, 1:1:N_j,AgeConditionalStats4.a2.Mean)
title('a2 (risky) Mean'); legend('1','2','3','4')

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
vfoptions1.lowmemory=0;

% lowmemory tests
vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=0;

%%
clear V1 V2 V1B V2B Policy1 Policy2 Policy1B Policy2B V1fromPolicy
%% Solve with grid-interpolation (a valid, more-accurate solution)
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

% V from Policy
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
vfoptions3.lowmemory=0;

% lowmemory tests
vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=0;

%%
clear V3 V4 V3B V4B Policy3 Policy4 Policy3B Policy4B V3fromPolicy
%% Big a_grid for accurate StationaryDist moments
fprintf('Big a_grid moments (EZ negative utils)\n')
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

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
disp([AllStats1.a1.Mean,AllStats3.a1.Mean])
disp([AllStats1.a2.Mean,AllStats3.a2.Mean])
disp([AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean])
disp([AgeConditionalStats1.a1.StdDeviation; AgeConditionalStats3.a1.StdDeviation])

clear Policy1b Policy3b StationaryDist1 StationaryDist3

% Also with DC
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

fprintf('With/without grid interp (with DC), should get much the same moments \n')
disp([AllStats2.a1.Mean,AllStats4.a1.Mean])
disp([AllStats2.a2.Mean,AllStats4.a2.Mean])
disp([AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean])
disp([AgeConditionalStats2.a1.StdDeviation; AgeConditionalStats4.a1.StdDeviation])

clear Policy2b Policy4b StationaryDist2 StationaryDist4

%% Age-conditional plots
figure(200+figure_c); % Case 3, negative-valued utils (Case 1 is figure_c)
subplot(3,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('Earnings Mean (EZ negative utils)'); legend('1','2','3','4')
subplot(3,1,2); plot(1:1:N_j,AgeConditionalStats1.a1.StdDeviation, 1:1:N_j,AgeConditionalStats2.a1.StdDeviation, 1:1:N_j,AgeConditionalStats3.a1.StdDeviation, 1:1:N_j,AgeConditionalStats4.a1.StdDeviation)
title('a1 (safe) Std Dev'); legend('1','2','3','4')
subplot(3,1,3); plot(1:1:N_j,AgeConditionalStats1.a2.Mean, 1:1:N_j,AgeConditionalStats2.a2.Mean, 1:1:N_j,AgeConditionalStats3.a2.Mean, 1:1:N_j,AgeConditionalStats4.a2.Mean)
title('a2 (risky) Mean'); legend('1','2','3','4')

clear AllStats1 AllStats2 AllStats3 AllStats4 AgeConditionalStats1 AgeConditionalStats2 AgeConditionalStats3 AgeConditionalStats4


%% Special Epstein-Zin tests
% The u-shock makes these sharper than in CoreFHorzTests: the collapses only hold exactly if
% the u-expectation sits INSIDE the joint certainty-equivalent (one CE over (u, semizprime)),
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

%% (iii) The same three collapse tests, under the grid interpolation layer
% The EZ GI layer linearly interpolates the transformed continuation EV, and linear
% interpolation commutes with the affine collapse transform, so these should be exact just
% like the no-GI versions. (For riskyasset withA1 the GI layer interpolates over a1prime.)
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
ReturnFn_cons_scaled=@(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactor) ezscalefactor*EZRiskyReturnFn_cons_d1_noz_noe_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
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
ReturnFn_posU_scaled=@(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactoru) ezscalefactoru*EZRiskyReturnFn_positiveUtils_d1_noz_noe_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
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
ReturnFn_negU_scaled=@(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactoru) ezscalefactoru*EZRiskyReturnFn_negativeUtils_d1_noz_noe_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
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
% trimmed to length Njs (only agej and kappa_j: everything else entering the ReturnFn, aprimeFn
% and SemiExoStateFn is age-independent). V and Policy must then be identical to the original
% model for periods 1,...,jstar-1. Each of the four solution methods gets a different jstar (the
% last of them uses jstar=N_j, so that one of them covers the retirement periods). Run for all
% three EZ cases. vfoptionsjs inherits all the riskyasset and semiz settings (aprimeFn/n_u/
% u_grid/pi_u, n_semiz/semiz_grid/SemiExoStateFn) via this file's per-method vfoptions.
% Note: mewj is age-dependent, but is only used for the agent distribution, which is not
% computed here, so it is left alone.
for ezcase=1:3
    vfoptionsv=vfoptions; % carries the riskyasset and semiz settings
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
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
    Vbase=Vbase(:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (mirroring the lowmemory values this file uses for this method in the standard sections)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;

    %% V_Jplus1, with divide-and-conquer
    jstar=round(2*N_j/3);
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    vfoptionsjs=vfoptions2;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
    Vbase=Vbase(:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (mirroring the lowmemory values this file uses for this method in the standard sections)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1 (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;

    %% V_Jplus1, with grid interpolation
    jstar=round(N_j/2);
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    vfoptionsjs=vfoptions3;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
    Vbase=Vbase(:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i, with GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i, with GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (mirroring the lowmemory values this file uses for this method in the standard sections)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1 (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;

    %% V_Jplus1, with divide-and-conquer and grid interpolation
    jstar=N_j;
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    vfoptionsjs=vfoptions4;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
    Vbase=Vbase(:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i, with DC+GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i, with DC+GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (mirroring the lowmemory values this file uses for this method in the standard sections)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;
end

clear Vbase Vshort Policybase Policyshort

%% (vi) Survival probabilities: vfoptions.survivalprobability
% sj is the probability of surviving from age j to j+1: it weights the continuation
% certainty-equivalent inside the EZ aggregator (combined with a warm-glow of bequests, 1-sj
% weights the bequest term; that is tested in (vii)). Params.sj is declining with sj(N_j)=0.
% Three legs (this file's collapse tests (i)-(iii), now with sj):
%  (vi).1 cons-units gamma=1/phi collapse with sj: EZ with survivalprobability='sj' equals vNM
%         (negativeUtils return fn, ezsigma=ezgamma) with DiscountFactorParamNames={'beta','sj'}:
%         Policy exact, V via V_EZ=((1-ezgamma)*V_vNM).^(1/(1-ezgamma)); also under the grid
%         interpolation layer.
%  (vi).2 utility-units EZriskaversion=0 collapse with sj (both utility signs): V and Policy
%         exact against vNM with DiscountFactorParamNames={'beta','sj'}; also under the grid
%         interpolation layer.
%  (vi).3 ValueFnFromPolicy with sj: FromPolicy==V (all three EZ cases; basic and under GI —
%         the first exercise of the sj path of the riskyasset EZ ValueFnFromPolicy).
% The u-shock is part of the joint certainty-equivalent here too, so these collapses also
% guard the placement of the u-expectation under sj.

% (vi).1 cons-units gamma=1/phi collapse with survival probabilities: EZ with sj (and no
% warm-glow, so death is worth zero on both sides) is vNM expected utility with age-dependent
% discounting beta*sj
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
% and the same under the grid interpolation layer
vfoptions1.gridinterplayer=1;
vfoptions1.ngridinterp=5;
vfoptions2.gridinterplayer=1;
vfoptions2.ngridinterp=5;
[V1a,Policy1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,{'beta','sj'},[],vfoptions2);
V1btransformed=((1-Params.ezgamma)*V1b).^(1/(1-Params.ezgamma));
fprintf('EZ gamma=1/phi collapse with sj (GI), Policy: should give zero: %2.8f \n',max(abs(Policy1a(:)-Policy1b(:))))
fprintf('EZ gamma=1/phi collapse with sj (GI), V after transform (relative): should be roughly 1e-13: %g \n',max(abs(V1a(:)-V1btransformed(:)))/max(abs(V1a(:))))
Params.ezphi=ezphi_store;
Params.ezsigma=ezsigma_store;
clear V1a V1b V1btransformed Policy1a Policy1b

% (vi).2 utility-units EZriskaversion=0 collapse with survival probabilities (both utility signs)
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
% and the same under the grid interpolation layer
vfoptions1.gridinterplayer=1;
vfoptions1.ngridinterp=5;
vfoptions2.gridinterplayer=1;
vfoptions2.ngridinterp=5;
[V2a,Policy2a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,{'beta','sj'},[],vfoptions2);
fprintf('EZ (positive utils) with EZriskaversion=0 and sj (GI): should give zero: %2.8f \n',max(abs(V2a(:)-V2b(:))))
fprintf('EZ (positive utils) with EZriskaversion=0 and sj (GI): should give zero: %2.8f \n',max(abs(Policy2a(:)-Policy2b(:))))
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
% and the same under the grid interpolation layer
vfoptions1.gridinterplayer=1;
vfoptions1.ngridinterp=5;
vfoptions2.gridinterplayer=1;
vfoptions2.ngridinterp=5;
[V3a,Policy3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,{'beta','sj'},[],vfoptions2);
fprintf('EZ (negative utils) with EZriskaversion=0 and sj (GI): should give zero: %2.8f \n',max(abs(V3a(:)-V3b(:))))
fprintf('EZ (negative utils) with EZriskaversion=0 and sj (GI): should give zero: %2.8f \n',max(abs(Policy3a(:)-Policy3b(:))))
Params.ezrisk=ezrisk_store;
clear V2a V2b V3a V3b Policy2a Policy2b Policy3a Policy3b

% (vi).3 ValueFnFromPolicy with survival probabilities: FromPolicy==V (basic and under GI)
for ezcase=1:3
    vfoptionsv=vfoptions; % carries the riskyasset and semiz settings
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
    fprintf('ValueFnFromPolicy with sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    vfoptionsGI2=vfoptionsv;
    vfoptionsGI2.gridinterplayer=1;
    vfoptionsGI2.ngridinterp=5;
    [V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsGI2);
    V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsGI2);
    fprintf('ValueFnFromPolicy with sj (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3fromPolicy(:)-V3(:))))
end
clear V1 V3 Policy1 Policy3 V1fromPolicy V3fromPolicy

%% (vii) Warm-glow of bequests: vfoptions.WarmGlowBequestsFn (De Nardi luxury-good form)
% The warm-glow fns (in EZRisky_ReturnFns) match each case's units, and in the riskyasset
% banks take a2prime (the risky-asset holding realized via the (d,u) lottery) as their only
% asset argument — the toolkit evaluates them on the a2 grid:
%   cons-units:     EZRiskyWarmGlowFn_cons(a2prime,wg1,wg2) (a consumption-equivalent, strictly
%                   positive; curvature comes from the EZ preferences)
%   positive utils: EZRiskyWarmGlowFn_positiveUtils(a2prime,wg1,wg2,wg3) (strictly positive)
%   negative utils: EZRiskyWarmGlowFn_negativeUtils(a2prime,wg1,wg2,wg3) (strictly negative)
% Five legs (two of them dropped for riskyasset withA1, see the comments at .3 and .5):
%  (vii).1 terminal-only default (WarmGlowBequestsFn set, no survivalprobability; the
%          dispatcher prints its assumed-terminal-only warning, which is expected output
%          here): basic==DC, GI==DC+GI, FromPolicy==V; and identical to the explicit
%          survivalprobability='sjterm' with sjterm=[ones(1,N_j-1),0].
%  (vii).2 warm-glow with declining sj: basic==DC, GI==DC+GI, FromPolicy==V, plus the
%          lowmemory rungs of each method.
%  (vii).3 exact collapse oracles with sj and warm-glow: DROPPED (see below).
%  (vii).4 V_Jplus1 mini-leg with sj and warm-glow (jstar=round(2*N_j/3); basic and DC).
%  (vii).5 N_j-1 warm-glow identity: EXCLUDED (see below).

% (vii).1 terminal-only default
Params2=Params;
Params2.sjterm=[ones(1,N_j-1),0];
for ezcase=1:3
    vfoptionsv=vfoptions; % carries the riskyasset and semiz settings
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
        vfoptionsv.WarmGlowBequestsFn=@(a2prime,wg1,wg2) EZRiskyWarmGlowFn_cons(a2prime,wg1,wg2);
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
        vfoptionsv.WarmGlowBequestsFn=@(a2prime,wg1,wg2,wg3) EZRiskyWarmGlowFn_positiveUtils(a2prime,wg1,wg2,wg3);
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
        vfoptionsv.WarmGlowBequestsFn=@(a2prime,wg1,wg2,wg3) EZRiskyWarmGlowFn_negativeUtils(a2prime,wg1,wg2,wg3);
    end
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
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    fprintf('Divide-and-conquer with terminal warm-glow [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('Divide-and-conquer with terminal warm-glow [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    [V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    [V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    fprintf('Divide-and-conquer with terminal warm-glow (with Grid Interp Layer) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V4(:))))
    fprintf('Divide-and-conquer with terminal warm-glow (with Grid Interp Layer) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy4(:))))
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
    fprintf('ValueFnFromPolicy with terminal warm-glow [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    % identical to the explicit terminal-only survival probabilities
    vfoptionssjt=vfoptionsv;
    vfoptionssjt.survivalprobability='sjterm';
    [V5,Policy5]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params2,DiscountFactorParamNames,[],vfoptionssjt);
    fprintf('terminal-only warm-glow vs explicit sjterm [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V5(:))))
    fprintf('terminal-only warm-glow vs explicit sjterm [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy5(:))))
end
clear V1 V2 V3 V4 V5 Policy1 Policy2 Policy3 Policy4 Policy5 V1fromPolicy

% (vii).2 warm-glow with declining survival probabilities
for ezcase=1:3
    vfoptionsv=vfoptions; % carries the riskyasset and semiz settings
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
        vfoptionsv.WarmGlowBequestsFn=@(a2prime,wg1,wg2) EZRiskyWarmGlowFn_cons(a2prime,wg1,wg2);
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
        vfoptionsv.WarmGlowBequestsFn=@(a2prime,wg1,wg2,wg3) EZRiskyWarmGlowFn_positiveUtils(a2prime,wg1,wg2,wg3);
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
        vfoptionsv.WarmGlowBequestsFn=@(a2prime,wg1,wg2,wg3) EZRiskyWarmGlowFn_negativeUtils(a2prime,wg1,wg2,wg3);
    end
    vfoptionsv.survivalprobability='sj';
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
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    fprintf('Divide-and-conquer with warm-glow and sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('Divide-and-conquer with warm-glow and sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    [V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    [V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    fprintf('Divide-and-conquer with warm-glow and sj (with Grid Interp Layer) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V4(:))))
    fprintf('Divide-and-conquer with warm-glow and sj (with Grid Interp Layer) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy4(:))))
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
    fprintf('ValueFnFromPolicy with warm-glow and sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    % lowmemory rungs (mirroring the lowmemory values this file uses in the standard sections)
    vfoptions1.lowmemory=1;
    [V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    fprintf('lowmemory=1 with warm-glow and sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V1B(:))))
    fprintf('lowmemory=1 with warm-glow and sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy1B(:))))
    vfoptions1.lowmemory=2;
    [V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    fprintf('lowmemory=2 with warm-glow and sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V1C(:))))
    fprintf('lowmemory=2 with warm-glow and sj [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy1C(:))))
    vfoptions1.lowmemory=0;
    vfoptions2.lowmemory=1;
    [V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    fprintf('lowmemory=1 with warm-glow and sj (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V2(:)-V2B(:))))
    fprintf('lowmemory=1 with warm-glow and sj (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy2(:)-Policy2B(:))))
    vfoptions2.lowmemory=2;
    [V2C,Policy2C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    fprintf('lowmemory=2 with warm-glow and sj (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V2(:)-V2C(:))))
    fprintf('lowmemory=2 with warm-glow and sj (with DC) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy2(:)-Policy2C(:))))
    vfoptions2.lowmemory=0;
    vfoptions3.lowmemory=1;
    [V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    fprintf('lowmemory=1 with warm-glow and sj (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V3B(:))))
    fprintf('lowmemory=1 with warm-glow and sj (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy3B(:))))
    vfoptions3.lowmemory=2;
    [V3C,Policy3C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    fprintf('lowmemory=2 with warm-glow and sj (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V3C(:))))
    fprintf('lowmemory=2 with warm-glow and sj (with GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy3C(:))))
    vfoptions3.lowmemory=0;
    vfoptions4.lowmemory=1;
    [V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    fprintf('lowmemory=1 with warm-glow and sj (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V4(:)-V4B(:))))
    fprintf('lowmemory=1 with warm-glow and sj (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy4(:)-Policy4B(:))))
    vfoptions4.lowmemory=2;
    [V4C,Policy4C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    fprintf('lowmemory=2 with warm-glow and sj (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V4(:)-V4C(:))))
    fprintf('lowmemory=2 with warm-glow and sj (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy4(:)-Policy4C(:))))
    vfoptions4.lowmemory=0;
end
clear V1 V2 V3 V4 V1B V1C V2B V2C V3B V3C V4B V4C Policy1 Policy2 Policy3 Policy4 Policy1B Policy1C Policy2B Policy2C Policy3B Policy3C Policy4B Policy4C V1fromPolicy

% (vii).3 exact collapse oracles with sj and warm-glow: DROPPED for riskyasset. The main-bank
% version builds a vNM oracle by adding beta*oneminussj*WarmGlowFn(aprime) to the period
% return fn, but there is no vNM riskyasset warm-glow support, and the warm-glow CANNOT enter
% a riskyasset ReturnFn: its argument is a2prime, which is realized via the (d,u) lottery
% AFTER the period return is evaluated, so no ReturnFn(...,a1prime,a1,a2,...) can see it.
% Exactness of the sj+warm-glow composition is instead anchored by the degenerateu bridge
% (the plain-vs-riskyasset cross-tests, where the lottery is degenerate).

% (vii).4 V_Jplus1 mini-leg with sj and warm-glow: V of period jstar as the terminal value fn
% of a shorter model, with survivalprobability and WarmGlowBequestsFn active (basic and DC).
% The age-dependent parameters trimmed to Njs now include sj (and oneminussj).
jstar=round(2*N_j/3);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
Paramsjs.sj=Params.sj(1:Njs);
Paramsjs.oneminussj=Params.oneminussj(1:Njs);
for ezcase=1:3
    vfoptionsv=vfoptions; % carries the riskyasset and semiz settings
    vfoptionsv.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptionsv.EZutils=0;
        vfoptionsv.EZriskaversion='ezgamma';
        vfoptionsv.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
        vfoptionsv.WarmGlowBequestsFn=@(a2prime,wg1,wg2) EZRiskyWarmGlowFn_cons(a2prime,wg1,wg2);
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=1;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
        vfoptionsv.WarmGlowBequestsFn=@(a2prime,wg1,wg2,wg3) EZRiskyWarmGlowFn_positiveUtils(a2prime,wg1,wg2,wg3);
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptionsv.EZutils=1;
        vfoptionsv.EZpositiveutility=0;
        vfoptionsv.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
        vfoptionsv.WarmGlowBequestsFn=@(a2prime,wg1,wg2,wg3) EZRiskyWarmGlowFn_negativeUtils(a2prime,wg1,wg2,wg3);
    end
    vfoptionsv.survivalprobability='sj';
    % basic
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsv);
    vfoptionsjs=vfoptionsv;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
    Vbase=Vbase(:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 with sj and warm-glow (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 with sj and warm-glow (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % with divide-and-conquer
    vfoptionsvDC=vfoptionsv;
    vfoptionsvDC.divideandconquer=1;
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsvDC);
    vfoptionsjs=vfoptionsvDC;
    vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
    Vbase=Vbase(:,:,:,1:Njs);
    Policybase=Policybase(:,:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 with sj and warm-glow (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 with sj and warm-glow (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
end
clear Vbase Vshort Policybase Policyshort

% (vii).5 N_j-1 warm-glow identity: EXCLUDED for withA1 — the terminal value fn depends on a1
% as well as a2, but WarmGlowBequestsFn only takes the single (a2prime) asset argument, so no
% closed-form warm-glow can reproduce V_Nj.

%% (viii) EZmortalityriskaversion (mortality risk aversion distinct from within-period risk aversion)
% All legs use survivalprobability='sj': with sj identically 1 there is no mortality risk for
% the mortality risk aversion to act on, so it would be vacuous (ezc8 only enters the raws
% through the sj-weighted continuation/warm-glow terms). Two legs:
%  (viii).1 identity: EZmortalityriskaversion equal to the within-period risk-aversion
%           coefficient (ezgamma for cons-units, ezrisk for utility-units) gives ezc8=1:
%           identical to leaving EZmortalityriskaversion unset (all three EZ cases).
%  (viii).2 distinct mortality risk aversion (ezmrisk=5): basic==DC, GI==DC+GI, FromPolicy==V.

% (viii).1 identity: mortality risk aversion equal to the within-period risk aversion
Params2=Params;
for ezcase=1:3
    vfoptionsv=vfoptions; % carries the riskyasset and semiz settings
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
    vfoptionsv=vfoptions; % carries the riskyasset and semiz settings
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
    vfoptions4=vfoptionsv;
    vfoptions4.divideandconquer=1;
    vfoptions4.gridinterplayer=1;
    vfoptions4.ngridinterp=5;
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    fprintf('Divide-and-conquer with sj and mortality risk aversion [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('Divide-and-conquer with sj and mortality risk aversion [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    [V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    [V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    fprintf('Divide-and-conquer with sj and mortality risk aversion (with Grid Interp Layer) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V4(:))))
    fprintf('Divide-and-conquer with sj and mortality risk aversion (with Grid Interp Layer) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy4(:))))
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
    fprintf('ValueFnFromPolicy with sj and mortality risk aversion [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
end
clear V1 V2 V3 V4 Policy1 Policy2 Policy3 Policy4 V1fromPolicy

%% (ix) vfoptions.EZoneminusbeta=2 (with survival probabilities)
% EZoneminusbeta=2 puts (1-sj*beta) on the this-period return: the mortality-adjusted version
% of the EZoneminusbeta=1 tested in (iv). Oracles:
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
ReturnFn_cons_scaled2=@(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactor2) ezscalefactor2*EZRiskyReturnFn_cons_d1_noz_noe_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
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
ReturnFn_posU_scaled2=@(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactoru2) ezscalefactoru2*EZRiskyReturnFn_positiveUtils_d1_noz_noe_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
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
ReturnFn_negU_scaled2=@(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost,ezscalefactoru2) ezscalefactoru2*EZRiskyReturnFn_negativeUtils_d1_noz_noe_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
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
