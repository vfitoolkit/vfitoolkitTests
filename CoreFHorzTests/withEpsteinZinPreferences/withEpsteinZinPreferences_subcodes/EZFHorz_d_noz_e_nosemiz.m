function output=EZFHorz_d_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Epstein-Zin mirror of QHDFHorz_d_noz_e_nosemiz.m (which mirrors CoreFHorz_d_noz_e_nosemiz.m).
% Where QH runs the model twice (Naive/Sophisticated), EZ runs it three times:
%   Case 1: consumption-units (traditional EZ, vfoptions.EZutils=0)
%   Case 2: utility-units, positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%   Case 3: utility-units, negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% followed by the EZ special tests (gamma=1/phi collapse; EZriskaversion=0 collapse; EZoneminusbeta=1).
% Two-endogenous-state (with2A) versions of these tests will be added later.

% Setup vfoptions and simoptions
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;
% Do the current setup
n_z=0;
z_grid=[];
pi_z=[];
% zeros assets, mid points for any shocks
jequaloneDist=zeros(n_a_big,vfoptions.n_e,'gpuArray'); % Note: based on n_a_big, not n_a
jequaloneDist(1,ceil(vfoptions.n_e/2))=1; % no assets, midpoint shock

ReturnFn_cons=@(d,aprime,a,e,r,w,kappa_j,varphi,agej,Jr,pension) EZReturnFn_cons_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,varphi,agej,Jr,pension);
ReturnFn_posU=@(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_positiveUtils_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
ReturnFn_negU=@(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_negativeUtils_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(d,aprime,a,e) a;
FnsToEvaluate.earnings=@(d,aprime,a,e,w,kappa_j) w*kappa_j*e*d;


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

%% (iv) vfoptions.EZoneminusbeta=1, tested against the DOCUMENTED behaviour
% FLAGGED, POSSIBLY-FAILING: reading the EZ raws, ezc1 is applied both to the this-period return
% AND again to the whole aggregated RHS before the final root (e.g. lines 114 and 117 of
% ValueFnIter_FHorz_EpsteinZin_nod_raw.m), i.e. an extra (1-beta) scaling beyond the documented
% '(1-beta)* term on the this period return'. These tests are written against the DOCUMENTED
% behaviour so a failure here drives a decision (fix the code, or document the normalization).
% See 'Open issues' in EZtests_full_coverage_proposal.md. Details:
% - consumption-units: EZ is homogeneous of degree 1 in the return, so EZoneminusbeta=1 with
%   return x should equal EZoneminusbeta=0 with return scaled by (1-beta)^(1/(1-1/ezphi)).
% - utility-units: V=(1-beta)*u+beta*(CE) is the same recursion as default with u scaled by (1-beta).
% (the scale factors enter the scaled return fns as parameters, since GPU arrayfun does
% not support anonymous functions that capture workspace variables)
Params.ezscalefactor=(1-Params.beta)^(1/(1-1/Params.ezphi));
Params.ezscalefactoru=1-Params.beta;

% consumption-units
ReturnFn_cons_scaled=@(d,aprime,a,e,r,w,kappa_j,varphi,agej,Jr,pension,ezscalefactor) ezscalefactor*EZReturnFn_cons_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,varphi,agej,Jr,pension);
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
vfoptions1.EZoneminusbeta=1;
[V4a,Policy4a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.EZoneminusbeta=0;
[V4b,Policy4b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons_scaled,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZoneminusbeta=1 (cons-units) vs manual scaling: should give zero (FLAGGED, may fail): %2.8f \n',max(abs(V4a(:)-V4b(:))))
fprintf('EZoneminusbeta=1 (cons-units) vs manual scaling: should give zero (FLAGGED, may fail): %2.8f \n',max(abs(Policy4a(:)-Policy4b(:))))

% utility-units, positive
ReturnFn_posU_scaled=@(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension,ezscalefactoru) ezscalefactoru*EZReturnFn_positiveUtils_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=1;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.EZoneminusbeta=1;
[V5a,Policy5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.EZoneminusbeta=0;
[V5b,Policy5b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU_scaled,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZoneminusbeta=1 (positive utils) vs manual scaling: should give zero (FLAGGED, may fail): %2.8f \n',max(abs(V5a(:)-V5b(:))))
fprintf('EZoneminusbeta=1 (positive utils) vs manual scaling: should give zero (FLAGGED, may fail): %2.8f \n',max(abs(Policy5a(:)-Policy5b(:))))

% utility-units, negative
ReturnFn_negU_scaled=@(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension,ezscalefactoru) ezscalefactoru*EZReturnFn_negativeUtils_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
vfoptions1=vfoptions;
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=1;
vfoptions1.EZpositiveutility=0;
vfoptions1.EZriskaversion='ezrisk';
vfoptions1.EZoneminusbeta=1;
[V6a,Policy6a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.EZoneminusbeta=0;
[V6b,Policy6b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU_scaled,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('EZoneminusbeta=1 (negative utils) vs manual scaling: should give zero (FLAGGED, may fail): %2.8f \n',max(abs(V6a(:)-V6b(:))))
fprintf('EZoneminusbeta=1 (negative utils) vs manual scaling: should give zero (FLAGGED, may fail): %2.8f \n',max(abs(Policy6a(:)-Policy6b(:))))

clear V4a V4b V5a V5b V6a V6b Policy4a Policy4b Policy5a Policy5b Policy6a Policy6b

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
