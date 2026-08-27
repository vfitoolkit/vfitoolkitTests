function output=EZFHorz_nod_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Epstein-Zin mirror of QHDFHorz_nod_noz_e_nosemiz.m (which mirrors CoreFHorz_nod_noz_e_nosemiz.m).
% Where QH runs the model twice (Naive/Sophisticated), EZ runs it three times:
%   Case 1: consumption-units (traditional EZ, vfoptions.EZutils=0)
%   Case 2: utility-units, positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%   Case 3: utility-units, negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% followed by the EZ special tests (gamma=1/phi collapse; EZriskaversion=0 collapse; EZoneminusbeta=1).
% plus special test (v): V_Jplus1 (V of period jstar as the terminal value fn of a shorter model).
% plus special tests (vi)-(ix):
%   (vi)   survival probabilities: sjones plumbing; the (i)/(ii) collapse oracles with declining sj
%          (the vNM oracle discounts by DiscountFactorParamNames={'beta','sj'}); cross-method agreement with sj
%   (vii)  warm-glow of bequests (De Nardi luxury-good form): terminal-only default and with declining sj;
%          exact vNM collapse oracles with sj+warm-glow; a V_Jplus1 mini-leg; the N_j-1 terminal
%          warm-glow identity (all three EZ cases)
%   (viii) EZmortalityriskaversion: identity when set equal to the case's own risk aversion; cross-method
%          agreement with a distinct mortality risk aversion
%   (ix)   EZoneminusbeta=2 versus manually scaling the return fn by the age-dependent (1-sj*beta) factor
% Two-endogenous-state (with2A) versions of these tests will be added later.

% Setup vfoptions and simoptions
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;
% Do the current setup
n_d=0;
d_grid=[];
n_z=0;
z_grid=[];
pi_z=[];
% zeros assets, mid points for any shocks
jequaloneDist=zeros(n_a_big,vfoptions.n_e,'gpuArray'); % Note: based on n_a_big, not n_a
jequaloneDist(1,ceil(vfoptions.n_e/2))=1; % no assets, midpoint shock

ReturnFn_cons=@(aprime,a,e,r,w,kappa_j,agej,Jr,pension) EZReturnFn_cons_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,agej,Jr,pension);
ReturnFn_posU=@(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension) EZReturnFn_positiveUtils_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_negU=@(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension) EZReturnFn_negativeUtils_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(aprime,a,e) a;
FnsToEvaluate.earnings=@(aprime,a,e,w,kappa_j) w*kappa_j*e;


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

% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=0;

%%
% V from Policy
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

clear V1 V1B V2 V2B Policy1 Policy1B Policy2 Policy2B V1fromPolicy

%% Solve with grid-interpolation
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

% Solve with divide-and-conquer, should give same answer
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

% lowmemory
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=0;

%%
% V from Policy
V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy (GI), this should be zero: %2.8f \n',max(abs(V3fromPolicy(:)-V3(:))))

clear V3 V3B V4 V4B Policy3 Policy3B Policy4 Policy4B V3fromPolicy

%% Use a really big a_grid, then the moments should be essentially the same with/without grid interpolation

[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy3b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions3);
AllStats3=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);
AgeConditionalStats3=LifeCycleProfiles_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.assets.Mean,AllStats3.assets.Mean]
[AllStats1.earnings.Gini,AllStats3.earnings.Gini]
[AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean]
[AgeConditionalStats1.assets.StdDeviation; AgeConditionalStats3.assets.StdDeviation]

clear V1b V3b Policy1b Policy3b StationaryDist1 StationaryDist3

% This is also true if using divide-and-conquer
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy2b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AllStats2=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);
AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);

[V4b,Policy4b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy4b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions4);
AllStats4=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);
AgeConditionalStats4=LifeCycleProfiles_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);

fprintf('With/without grid interp, should get much the same moments (for big a_grid; with divide-and-conquer) \n')
[AllStats2.assets.Mean,AllStats4.assets.Mean]
[AllStats2.earnings.Gini,AllStats4.earnings.Gini]
[AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean]
[AgeConditionalStats2.assets.StdDeviation; AgeConditionalStats4.assets.StdDeviation]

clear V2b V4b StationaryDist2 StationaryDist4 % Policy2b Policy4b

%% Do some graphs of the age-conditional to see them
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('Earnings Mean')
legend('1','2','3','4')
subplot(2,1,2); plot(1:1:N_j,AgeConditionalStats1.assets.StdDeviation, 1:1:N_j,AgeConditionalStats2.assets.StdDeviation, 1:1:N_j,AgeConditionalStats3.assets.StdDeviation, 1:1:N_j,AgeConditionalStats4.assets.StdDeviation)
title('Assets Std Dev')
legend('1','2','3','4')
% ylim([0,0.01]) % If you want to make graph look nicer


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

% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=0;

%%
% V from Policy
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

clear V1 V1B V2 V2B Policy1 Policy1B Policy2 Policy2B V1fromPolicy

%% Solve with grid-interpolation
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

% Solve with divide-and-conquer, should give same answer
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

% lowmemory
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=0;

%%
% V from Policy
V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy (GI), this should be zero: %2.8f \n',max(abs(V3fromPolicy(:)-V3(:))))

clear V3 V3B V4 V4B Policy3 Policy3B Policy4 Policy4B V3fromPolicy

%% Use a really big a_grid, then the moments should be essentially the same with/without grid interpolation

[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy3b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions3);
AllStats3=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);
AgeConditionalStats3=LifeCycleProfiles_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.assets.Mean,AllStats3.assets.Mean]
[AllStats1.earnings.Gini,AllStats3.earnings.Gini]
[AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean]
[AgeConditionalStats1.assets.StdDeviation; AgeConditionalStats3.assets.StdDeviation]

clear V1b V3b Policy1b Policy3b StationaryDist1 StationaryDist3

% This is also true if using divide-and-conquer
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy2b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AllStats2=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);
AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);

[V4b,Policy4b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy4b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions4);
AllStats4=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);
AgeConditionalStats4=LifeCycleProfiles_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);

fprintf('With/without grid interp, should get much the same moments (for big a_grid; with divide-and-conquer) \n')
[AllStats2.assets.Mean,AllStats4.assets.Mean]
[AllStats2.earnings.Gini,AllStats4.earnings.Gini]
[AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean]
[AgeConditionalStats2.assets.StdDeviation; AgeConditionalStats4.assets.StdDeviation]

clear V2b V4b StationaryDist2 StationaryDist4 % Policy2b Policy4b

%% Do some graphs of the age-conditional to see them
fig=figure(100+figure_c); % Case 2, positive-valued utils (Case 1 is figure_c)
subplot(2,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('Earnings Mean')
legend('1','2','3','4')
subplot(2,1,2); plot(1:1:N_j,AgeConditionalStats1.assets.StdDeviation, 1:1:N_j,AgeConditionalStats2.assets.StdDeviation, 1:1:N_j,AgeConditionalStats3.assets.StdDeviation, 1:1:N_j,AgeConditionalStats4.assets.StdDeviation)
title('Assets Std Dev')
legend('1','2','3','4')
% ylim([0,0.01]) % If you want to make graph look nicer


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

% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=0;

%%
% V from Policy
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

clear V1 V1B V2 V2B Policy1 Policy1B Policy2 Policy2B V1fromPolicy

%% Solve with grid-interpolation
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

% Solve with divide-and-conquer, should give same answer
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

% lowmemory
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=0;

%%
% V from Policy
V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy (GI), this should be zero: %2.8f \n',max(abs(V3fromPolicy(:)-V3(:))))

clear V3 V3B V4 V4B Policy3 Policy3B Policy4 Policy4B V3fromPolicy

%% Use a really big a_grid, then the moments should be essentially the same with/without grid interpolation

[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy3b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions3);
AllStats3=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);
AgeConditionalStats3=LifeCycleProfiles_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.assets.Mean,AllStats3.assets.Mean]
[AllStats1.earnings.Gini,AllStats3.earnings.Gini]
[AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean]
[AgeConditionalStats1.assets.StdDeviation; AgeConditionalStats3.assets.StdDeviation]

clear V1b V3b Policy1b Policy3b StationaryDist1 StationaryDist3

% This is also true if using divide-and-conquer
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy2b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AllStats2=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);
AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);

[V4b,Policy4b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy4b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions4);
AllStats4=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);
AgeConditionalStats4=LifeCycleProfiles_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);

fprintf('With/without grid interp, should get much the same moments (for big a_grid; with divide-and-conquer) \n')
[AllStats2.assets.Mean,AllStats4.assets.Mean]
[AllStats2.earnings.Gini,AllStats4.earnings.Gini]
[AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean]
[AgeConditionalStats2.assets.StdDeviation; AgeConditionalStats4.assets.StdDeviation]

clear V2b V4b StationaryDist2 StationaryDist4 % Policy2b Policy4b

%% Do some graphs of the age-conditional to see them
fig=figure(200+figure_c); % Case 3, negative-valued utils (Case 1 is figure_c)
subplot(2,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('Earnings Mean')
legend('1','2','3','4')
subplot(2,1,2); plot(1:1:N_j,AgeConditionalStats1.assets.StdDeviation, 1:1:N_j,AgeConditionalStats2.assets.StdDeviation, 1:1:N_j,AgeConditionalStats3.assets.StdDeviation, 1:1:N_j,AgeConditionalStats4.assets.StdDeviation)
title('Assets Std Dev')
legend('1','2','3','4')
% ylim([0,0.01]) % If you want to make graph look nicer


%% Special Epstein-Zin tests (these play the role the beta0=1 tests play in the QH test suite)

%% (i) Consumption-units: gamma=1/phi collapses Epstein-Zin to standard vNM expected utility
% With vfoptions.EZutils=0 and EZriskaversion=1/EZeis, the EZ value fn is an increasing transform
% of the vNM problem with period utility (x^(1-ezgamma))/(1-ezgamma), which is exactly
% EZReturnFn_negativeUtils with ezsigma=ezgamma. So: identical Policy (exactly), and
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
% The EZ GI layer linearly interpolates the transformed continuation EV=E[(ezc4*V)^ezc5], and
% linear interpolation commutes with the affine collapse transform, so these should be exact
% just like the no-GI versions.
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

%% (iv) vfoptions.EZoneminusbeta=1
% (Flag retired 2026-08-25: these tests were originally flagged may-fail against an ezc1
% double-application in the raws; that was resolved by removing the outer ezc1 from the raws
% [2026-08-11], and the tests have passed exactly on GPU ever since.) Oracles:
% - consumption-units: EZ is homogeneous of degree 1 in the return, so EZoneminusbeta=1 with
%   return x should equal EZoneminusbeta=0 with return scaled by (1-beta)^(1/(1-1/ezphi)).
% - utility-units: V=(1-beta)*u+beta*(CE) is the same recursion as default with u scaled by (1-beta).
% (the scale factors enter the scaled return fns as parameters, since GPU arrayfun does
% not support anonymous functions that capture workspace variables)
Params.ezscalefactor=(1-Params.beta)^(1/(1-1/Params.ezphi));
Params.ezscalefactoru=1-Params.beta;

% consumption-units
ReturnFn_cons_scaled=@(aprime,a,e,r,w,kappa_j,agej,Jr,pension,ezscalefactor) ezscalefactor*EZReturnFn_cons_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,agej,Jr,pension);
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
ReturnFn_posU_scaled=@(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension,ezscalefactoru) ezscalefactoru*EZReturnFn_positiveUtils_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension);
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
ReturnFn_negU_scaled=@(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension,ezscalefactoru) ezscalefactoru*EZReturnFn_negativeUtils_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension);
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
% trimmed to length Njs. V and Policy must then be identical to the original model for periods
% 1,...,jstar-1. Each of the four solution methods gets a different jstar (the last of them uses
% jstar=N_j, so that one of them covers the retirement periods). Run for all three EZ cases.
% Note: mewj is age-dependent, but is only used for the agent distribution, which is not
% computed here, so it is left alone.
for ezcase=1:3
    vfoptionsv=struct();
    vfoptionsv.n_e=vfoptionsbaseline.n_e;
    vfoptionsv.e_grid=vfoptionsbaseline.e_grid;
    vfoptionsv.pi_e=vfoptionsbaseline.pi_e;
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
    vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
    Vbase=Vbase(:,:,1:Njs);
    Policybase=Policybase(:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
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
    vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
    Vbase=Vbase(:,:,1:Njs);
    Policybase=Policybase(:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
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
    vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
    Vbase=Vbase(:,:,1:Njs);
    Policybase=Policybase(:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i, with GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i, with GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
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
    vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
    Vbase=Vbase(:,:,1:Njs);
    Policybase=Policybase(:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i, with DC+GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i, with DC+GI) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
    vfoptionsjs.lowmemory=1;
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1, lowmemory=1 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1, lowmemory=1 (with DC+GI) [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policybase(:)-Policyshort(:))))
    vfoptionsjs.lowmemory=0;
end

clear Vbase Vshort Policybase Policyshort

%% (vi) Survival probabilities (vfoptions.survivalprobability)
% Params.sj is declining with sj(N_j)=0 (set in CoreFHorzEZTests.m). When vfoptions.survivalprobability
% is not set the EZ codes use sj=ones(N_j,1) internally, so survivalprobability='sjones' (all ones)
% must reproduce the baseline solve exactly (pure plumbing). Then the (i)/(ii) collapse oracles are
% repeated WITH the declining sj: once EZ collapses to vNM, the survival probability is just an
% age-dependent discount factor, so the vNM oracle uses DiscountFactorParamNames={'beta','sj'}.
Params2=Params;
Params2.sjones=ones(1,N_j);
DiscountFactorParamNames2={'beta','sj'}; % beta*sj: age-dependent discounting for the vNM oracles

% (vi).1 plumbing: survivalprobability='sjones' vs not setting survivalprobability, all three EZ cases
for ezcase=1:3
    vfoptionsv=struct();
    vfoptionsv.n_e=vfoptionsbaseline.n_e;
    vfoptionsv.e_grid=vfoptionsbaseline.e_grid;
    vfoptionsv.pi_e=vfoptionsbaseline.pi_e;
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
% the vNM oracle discounts by beta*sj). Policy exact; V via the usual transform. Also repeated under
% the grid interpolation layer (linear interpolation commutes with the affine collapse transform,
% so it should stay exact).
ezphi_store=Params.ezphi;
ezsigma_store=Params.ezsigma;
Params.ezphi=1/Params.ezgamma;
Params.ezsigma=Params.ezgamma;

vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
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
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
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
% the vNM oracle discounts by beta*sj). V and Policy exact, both signs of the utility fn.
ezrisk_store=Params.ezrisk;
Params.ezrisk=0;

% positive-valued utility fn
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
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
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
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

% (vi).4 cross-method agreement with the declining sj, all three EZ cases: basic==DC, GI==DC+GI,
% and ValueFnFromPolicy on the basic solve (the first exercise of the FromPolicy sj path).
for ezcase=1:3
    vfoptionsv=struct();
    vfoptionsv.n_e=vfoptionsbaseline.n_e;
    vfoptionsv.e_grid=vfoptionsbaseline.e_grid;
    vfoptionsv.pi_e=vfoptionsbaseline.pi_e;
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
    fprintf('sj, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    fprintf('sj, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('sj, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    fprintf('sj, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V4(:))))
    fprintf('sj, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy4(:))))
end
clear V1 V2 V3 V4 Policy1 Policy2 Policy3 Policy4 V1fromPolicy

%% (vii) Warm-glow of bequests (De Nardi luxury-good form), vfoptions.WarmGlowBequestsFn
% Each EZ case gets a warm-glow fn matching its sign/units convention (see EZ_ReturnFns): cons-units
% gets EZWarmGlowFn_cons (a consumption-equivalent, strictly positive; curvature comes from the EZ
% preferences), positive/negative utils get EZWarmGlowFn_positiveUtils/EZWarmGlowFn_negativeUtils
% (strictly positive/strictly negative, nonzero at aprime=0 so the WG==0 mask conventions are avoided).
WGFn_cons=@(aprime,wg1,wg2) EZWarmGlowFn_cons(aprime,wg1,wg2);
WGFn_posU=@(aprime,wg1,wg2,wg3) EZWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3);
WGFn_negU=@(aprime,wg1,wg2,wg3) EZWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);
Params2.sjterm=[ones(1,N_j-1),0]; % the sj the toolkit defaults to when WarmGlowBequestsFn is set but survivalprobability is not

% (vii).1 terminal-only warm-glow: WarmGlowBequestsFn set but NO survivalprobability, so the toolkit
% defaults to warm-glow only at the end of the final period (each solve prints the expected warning).
% All three EZ cases: FromPolicy, basic==DC, GI==DC+GI, and the identity with explicitly passing
% survivalprobability='sjterm'.
for ezcase=1:3
    vfoptionsv=struct();
    vfoptionsv.n_e=vfoptionsbaseline.n_e;
    vfoptionsv.e_grid=vfoptionsbaseline.e_grid;
    vfoptionsv.pi_e=vfoptionsbaseline.pi_e;
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
    fprintf('warm-glow terminal-only, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
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
clear V1 V2 V3 V4 V5 Policy1 Policy2 Policy3 Policy4 Policy5 V1fromPolicy

% (vii).2 warm-glow at every age: WarmGlowBequestsFn plus the declining sj (warm-glow weight 1-sj(j)
% at every age, weight exactly one in the final period as sj(N_j)=0). All three EZ cases: FromPolicy,
% basic==DC, GI==DC+GI.
for ezcase=1:3
    vfoptionsv=struct();
    vfoptionsv.n_e=vfoptionsbaseline.n_e;
    vfoptionsv.e_grid=vfoptionsbaseline.e_grid;
    vfoptionsv.pi_e=vfoptionsbaseline.pi_e;
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
    fprintf('warm-glow with sj, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    fprintf('warm-glow with sj, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('warm-glow with sj, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    fprintf('warm-glow with sj, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V4(:))))
    fprintf('warm-glow with sj, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy4(:))))
    % lowmemory rung on the sj+warm-glow solve (mirroring the lowmemory values the standard
    % sections exercise): covers the main-loop lowmemory warm-glow sites of the e-raws
    vfoptions1.lowmemory=1;
    [V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    fprintf('warm-glow with sj, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V1B(:))))
    fprintf('warm-glow with sj, lowmemory=1 [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy1B(:))))
    vfoptions1.lowmemory=0;
end
clear V1 V2 V3 V4 V1B Policy1 Policy2 Policy3 Policy4 Policy1B V1fromPolicy

% (vii).3 exact collapse oracles with sj AND warm-glow.
% (vii).3a cons-units gamma=1/phi with sj and warm-glow: cons-units became valid under the
% inside-the-root terminal warm-glow convention (Kraft-Munk-Weiss 2022): the terminal age is now the
% interior recursion with the warm-glow as its continuation, so the combined-wrapper vNM oracle's
% W-recursion (W=F^x+beta*sj*E[W']+beta*(1-sj)*WG^x, x=1-ezgamma) holds at every age INCLUDING N_j.
% As in (vi).2 the vNM oracle is the negativeUtils return fn with ezsigma=ezgamma and
% DiscountFactorParamNames={'beta','sj'}, here plus beta*oneminussj*(WG_cons^(1-ezsigma))/(1-ezsigma)
% (the cons-units warm-glow through the same crosswalk as the return fn; the scale factors enter as
% trailing Params args). Policy exact; V via V_EZ=((1-ezgamma)*V_vNM).^(1/(1-ezgamma)).
ezphi_store=Params.ezphi;
ezsigma_store=Params.ezsigma;
Params.ezphi=1/Params.ezgamma;
Params.ezsigma=Params.ezgamma;
ReturnFn_consWGoracle=@(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension,beta,oneminussj,wg1,wg2) EZReturnFn_negativeUtils_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension)+beta*oneminussj*((EZWarmGlowFn_cons(aprime,wg1,wg2)^(1-ezsigma))/(1-ezsigma));
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
vfoptions1.survivalprobability='sj';
vfoptions1.WarmGlowBequestsFn=WGFn_cons;
[V1a,Policy1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
vfoptions1.exoticpreferences='None';
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_consWGoracle,Params,DiscountFactorParamNames2,[],vfoptions1);
V1btransformed=((1-Params.ezgamma)*V1b).^(1/(1-Params.ezgamma));
fprintf('EZ gamma=1/phi collapse, sj+warm-glow, Policy [EZ cons-units]: should give zero: %2.8f \n',max(abs(Policy1a(:)-Policy1b(:))))
fprintf('EZ gamma=1/phi collapse, sj+warm-glow, V after transform (relative) [EZ cons-units]: should be roughly 1e-13: %g \n',max(abs(V1a(:)-V1btransformed(:)))/max(abs(V1a(:))))
Params.ezphi=ezphi_store;
Params.ezsigma=ezsigma_store;
clear V1a V1b V1btransformed Policy1a Policy1b

% (vii).3b utility-units EZriskaversion=0 with sj and warm-glow: the vNM oracle adds
% beta*oneminussj*WG(aprime) to the return fn (the warm-glow is already in utility units).
% V and Policy exact directly (no transform), both signs of the utility fn.
ezrisk_store=Params.ezrisk;
Params.ezrisk=0;
ReturnFn_posUWGoracle=@(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension,beta,oneminussj,wg1,wg2,wg3) EZReturnFn_positiveUtils_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension)+beta*oneminussj*EZWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3);
ReturnFn_negUWGoracle=@(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension,beta,oneminussj,wg1,wg2,wg3) EZReturnFn_negativeUtils_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension)+beta*oneminussj*EZWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);
% positive-valued utility fn
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=1;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.survivalprobability='sj';
vfoptions1.WarmGlowBequestsFn=WGFn_posU;
[V2a,Policy2a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
vfoptions1.exoticpreferences='None';
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posUWGoracle,Params,DiscountFactorParamNames2,[],vfoptions1);
fprintf('EZ with EZriskaversion=0, sj+warm-glow [EZ positive utils]: should give zero: %2.8f \n',max(abs(V2a(:)-V2b(:))))
fprintf('EZ with EZriskaversion=0, sj+warm-glow [EZ positive utils]: should give zero: %2.8f \n',max(abs(Policy2a(:)-Policy2b(:))))
% negative-valued utility fn
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.survivalprobability='sj';
vfoptions1.WarmGlowBequestsFn=WGFn_negU;
[V3a,Policy3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
vfoptions1.exoticpreferences='None';
[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negUWGoracle,Params,DiscountFactorParamNames2,[],vfoptions1);
fprintf('EZ with EZriskaversion=0, sj+warm-glow [EZ negative utils]: should give zero: %2.8f \n',max(abs(V3a(:)-V3b(:))))
fprintf('EZ with EZriskaversion=0, sj+warm-glow [EZ negative utils]: should give zero: %2.8f \n',max(abs(Policy3a(:)-Policy3b(:))))
Params.ezrisk=ezrisk_store;
clear V2a V2b V3a V3b Policy2a Policy2b Policy3a Policy3b

% (vii).4 V_Jplus1 mini-leg with sj and warm-glow on: as (v) but only one jstar and only the basic
% and DC methods. The age-dependent sj and oneminussj are trimmed to the shorter horizon too (the
% toolkit checks length(sj) against the N_j of the model being solved).
for ezcase=1:3
    vfoptionsv=struct();
    vfoptionsv.n_e=vfoptionsbaseline.n_e;
    vfoptionsv.e_grid=vfoptionsbaseline.e_grid;
    vfoptionsv.pi_e=vfoptionsbaseline.pi_e;
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
    vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
    Vbase=Vbase(:,:,1:Njs);
    Policybase=Policybase(:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
    % with divide-and-conquer
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    vfoptionsjs=vfoptions2;
    vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
    Vbase=Vbase(:,:,1:Njs);
    Policybase=Policybase(:,:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 with sj+warm-glow (jstar=%i, with DC) [EZ %s], this should be zero: %2.8f \n',jstar,casestr,max(abs(Policybase(:)-Policyshort(:))))
end
clear Vbase Vshort Policybase Policyshort

% (vii).5 the N_j-1 terminal warm-glow identity, ALL THREE EZ CASES. Solve the full N_j model with
% no warm-glow and no survival probability; then solve an N_j-1 model whose terminal-only warm-glow
% (weight exactly 1, since with no survivalprobability the default sj=[ones,0] applies) is the
% CLOSED FORM of the full model's terminal value fn: age N_j is retired, so the terminal budget
% c=(1+r)*a+pension-aprime is z/e-independent, and the optimal terminal choices are aprime=0
% (utility is increasing in c, and a_grid(1)=0). V and Policy at ages 1,...,N_j-1 must then
% reproduce the full model's (the terminal V is shock-independent, so the certainty-equivalent the
% full model takes over it is over a degenerate distribution).
% cons-units became valid under the inside-the-root terminal warm-glow convention (Kraft-Munk-Weiss
% 2022): the cons-units terminal V is (ezc1*F^ezc2)^ezc7=F itself (ezc1=1, ezc2*ezc7=1), so its
% warm-glow fn is the same closed-form terminal composite with NO utility wrap.
WGterm_cons=@(aprime,r,pension) (1+r)*aprime+pension; % the terminal composite itself (cons-units terminal V=(ezc1*F^ezc2)^ezc7=F)
WGterm_posU=@(aprime,ezsigma,r,pension) ((1+((1+r)*aprime+pension))^(1-ezsigma)-1)/(1-ezsigma);
WGterm_negU=@(aprime,ezsigma,r,pension) (((1+r)*aprime+pension)^(1-ezsigma))/(1-ezsigma);
Njs2=N_j-1;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs2);
Paramsjs.kappa_j=Params.kappa_j(1:Njs2);
% consumption-units
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
Vbase=Vbase(:,:,1:Njs2);
Policybase=Policybase(:,:,:,1:Njs2);
vfoptions1.WarmGlowBequestsFn=WGterm_cons;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs2,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Paramsjs,DiscountFactorParamNames,[],vfoptions1);
fprintf('N_j-1 terminal warm-glow identity [EZ cons-units], this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('N_j-1 terminal warm-glow identity [EZ cons-units], this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
% positive-valued utility fn
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=1;
vfoptions1.EZriskaversion='ezrisk';
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
Vbase=Vbase(:,:,1:Njs2);
Policybase=Policybase(:,:,:,1:Njs2);
vfoptions1.WarmGlowBequestsFn=WGterm_posU;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs2,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Paramsjs,DiscountFactorParamNames,[],vfoptions1);
fprintf('N_j-1 terminal warm-glow identity [EZ positive utils], this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('N_j-1 terminal warm-glow identity [EZ positive utils], this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
% negative-valued utility fn
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
Vbase=Vbase(:,:,1:Njs2);
Policybase=Policybase(:,:,:,1:Njs2);
vfoptions1.WarmGlowBequestsFn=WGterm_negU;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs2,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Paramsjs,DiscountFactorParamNames,[],vfoptions1);
fprintf('N_j-1 terminal warm-glow identity [EZ negative utils], this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('N_j-1 terminal warm-glow identity [EZ negative utils], this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
clear Vbase Vshort Policybase Policyshort

%% (viii) vfoptions.EZmortalityriskaversion (double Epstein-Zin: separate mortality risk aversion)
% Only meaningful together with survival probabilities, so survivalprobability='sj' throughout this
% section (with sj==1 the ^ezc8 exactly cancels the ^ezc6 in the survival aggregator and the option
% is vacuous).
% (viii).1 identity: setting EZmortalityriskaversion EQUAL to the case's own risk aversion (ezgamma
% for cons-units, ezrisk for utility-units) gives ezc8=1, which must reproduce the unset default
% exactly.
for ezcase=1:3
    vfoptionsv=struct();
    vfoptionsv.n_e=vfoptionsbaseline.n_e;
    vfoptionsv.e_grid=vfoptionsbaseline.e_grid;
    vfoptionsv.pi_e=vfoptionsbaseline.pi_e;
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
% ^ezc8 sites. All three EZ cases: FromPolicy, basic==DC, GI==DC+GI.
for ezcase=1:3
    vfoptionsv=struct();
    vfoptionsv.n_e=vfoptionsbaseline.n_e;
    vfoptionsv.e_grid=vfoptionsbaseline.e_grid;
    vfoptionsv.pi_e=vfoptionsbaseline.pi_e;
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
    fprintf('EZmortalityriskaversion distinct, ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))
    fprintf('EZmortalityriskaversion distinct, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1(:)-V2(:))))
    fprintf('EZmortalityriskaversion distinct, basic vs DC [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy1(:)-Policy2(:))))
    fprintf('EZmortalityriskaversion distinct, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V3(:)-V4(:))))
    fprintf('EZmortalityriskaversion distinct, GI vs DC+GI [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(Policy3(:)-Policy4(:))))
end
clear V1 V2 V3 V4 Policy1 Policy2 Policy3 Policy4 V1fromPolicy

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
ReturnFn_cons_scaled2=@(aprime,a,e,r,w,kappa_j,agej,Jr,pension,ezscalefactor2) ezscalefactor2*EZReturnFn_cons_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,agej,Jr,pension);
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
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
ReturnFn_posU_scaled2=@(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension,ezscalefactoru2) ezscalefactoru2*EZReturnFn_positiveUtils_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension);
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
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
ReturnFn_negU_scaled2=@(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension,ezscalefactoru2) ezscalefactoru2*EZReturnFn_negativeUtils_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension);
vfoptions1=struct();
vfoptions1.n_e=vfoptionsbaseline.n_e;
vfoptions1.e_grid=vfoptionsbaseline.e_grid;
vfoptions1.pi_e=vfoptionsbaseline.pi_e;
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
