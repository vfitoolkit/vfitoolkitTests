function output=CoreInfHorz_nod_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% Do the current setup
n_d=0; d_grid=[];

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;

%% Baseline VFI
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

% Note: divide-and-conquer is not tested for InfHorz (too slow to be usable);
% it is covered in the InfHorz-TPath test bank.

%% lowmemory: lowmemory=1 should give the same V and Policy as lowmemory=0
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

%% Grid-interpolation layer
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

PolicyVals3=PolicyInd2Val_InfHorz(Policy3,n_d,n_a,n_z,d_grid,a_grid,vfoptions3);

V3fromPolicy=ValueFnFromPolicy_InfHorz(Policy3,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy with grid interp, this should be zero: %2.8f \n',max(abs(V3fromPolicy(:)-V3(:))))

%% lowmemory with grid interp
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=0;

%%
clear V1 V3 V1B V3B Policy1 Policy3 Policy1B Policy3B PolicyVals1 PolicyVals3 V1fromPolicy V3fromPolicy

%% Big a_grid: the stationary dist and moments should be essentially the same with/without grid interpolation
[~,Policy1b]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_InfHorz(Policy1b,n_d,n_a_big,n_z,pi_z,simoptions1,Params,[]);
AllStats1=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
AggVars1=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
AutoCorr1=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,simoptions1);
CrossSectionCorr1=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);

[~,Policy3b]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist3=StationaryDist_InfHorz(Policy3b,n_d,n_a_big,n_z,pi_z,simoptions3,Params,[]);
AllStats3=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions3);
AggVars3=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.assets.Mean,AllStats3.assets.Mean]
[AllStats1.earnings.Gini,AllStats3.earnings.Gini]
[AllStats1.assets.StdDeviation,AllStats3.assets.StdDeviation]
[AggVars1.assets.Mean,AggVars3.assets.Mean]

%% Check the remaining dist commands run without issue
ValuesOnGrid1=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
ProbDensityFns1=EvalFnOnAgentDist_ProbDensityFn_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
fprintf('ValuesOnGrid, ProbDensityFn all ran \n')

%% SimTimeSeries: cross-validate the AgentDist statistics with a long simulated time series
% (Monte Carlo, so these are roughly-equal checks, not machine precision.)
simoptionsTS=simoptions1;
simoptionsTS.simperiods=10^5;
simoptionsTS.burnin=10^3;
TimeSeries1=SimTimeSeriesValues_InfHorz(Policy1b,FnsToEvaluate,Params,n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,simoptionsTS);
tsA=TimeSeries1.assets(:); tsE=TimeSeries1.earnings(:);

fprintf('SimTimeSeries mean should roughly match AllStats/AggVars mean (Monte Carlo) \n')
[AllStats1.assets.Mean, AggVars1.assets.Mean, mean(tsA)]
[AllStats1.earnings.Mean, AggVars1.earnings.Mean, mean(tsE)]
fprintf('SimTimeSeries std should roughly match AllStats std \n')
[AllStats1.assets.StdDeviation, std(tsA,1)]
[AllStats1.earnings.StdDeviation, std(tsE,1)]
fprintf('SimTimeSeries lag-1 autocorrelation should roughly match AutoCorrTransProbs.AutoCorrelation \n')
tempA=corrcoef(tsA(1:end-1),tsA(2:end)); tempE=corrcoef(tsE(1:end-1),tsE(2:end));
[AutoCorr1.assets.AutoCorrelation, tempA(1,2)]
[AutoCorr1.earnings.AutoCorrelation, tempE(1,2)]
fprintf('SimTimeSeries cross-section corr(earnings,assets) should roughly match CrossSectionCovarCorr \n')
tempC=corrcoef(tsE,tsA);
[CrossSectionCorr1.earnings.CorrelationWith.assets, tempC(1,2)]

%% Plot the stationary dist over assets (with/without grid interp)
fig=figure(figure_c);
plot(a_grid_big,cumsum(sum(StationaryDist1,2)), a_grid_big,cumsum(sum(StationaryDist3,2)))
title('CDF of assets: without vs with grid interp'); legend('1','3')

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
