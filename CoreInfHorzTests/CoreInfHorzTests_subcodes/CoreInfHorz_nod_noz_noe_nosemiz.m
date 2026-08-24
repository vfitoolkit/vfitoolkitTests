function output=CoreInfHorz_nod_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% Do the current setup
n_d=0; d_grid=[];
n_z=0; z_grid=[]; pi_z=[];

ReturnFn=@(aprime,a,r,w,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(aprime,a) a;
FnsToEvaluate.earnings=@(aprime,a,w) w;

%% Baseline VFI
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

% Note: divide-and-conquer is not tested for InfHorz (too slow to be usable);
% it is covered in the InfHorz-TPath test bank.

% Note: lowmemory is not tested in the noz cases (lowmemory=1 loops over z,
% of which there is none here); it is tested in the z cases.

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

%%
clear V1 V3 Policy1 Policy3 PolicyVals1 PolicyVals3 V1fromPolicy V3fromPolicy

%% Big a_grid: the stationary dist and moments should be essentially the same with/without grid interpolation
[~,Policy1b]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_InfHorz(Policy1b,n_d,n_a_big,n_z,pi_z,simoptions1,Params,[]);
AllStats1=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
AggVars1=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);

[~,Policy3b]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist3=StationaryDist_InfHorz(Policy3b,n_d,n_a_big,n_z,pi_z,simoptions3,Params,[]);
AllStats3=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions3);
AggVars3=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.assets.Mean,AllStats3.assets.Mean]
[AllStats1.earnings.Mean,AllStats3.earnings.Mean]
[AggVars1.assets.Mean,AggVars3.assets.Mean]

%% Check the various other dist commands run without issue
% (No exogenous shock => the stationary dist is degenerate, so these are 'runs
%  without error' checks here; the z cases cross-validate the actual numbers.)
AutoCorrTransProbs1=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,simoptions1);
CrossSectionCorr1=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
ValuesOnGrid1=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
ProbDensityFns1=EvalFnOnAgentDist_ProbDensityFn_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
fprintf('AutoCorrTransProbs, CrossSectionCovarCorr, ValuesOnGrid, ProbDensityFn all ran \n')

%% SimTimeSeries: the time-series mean should reproduce the AggVars mean (Monte Carlo, so approximate)
simoptionsTS=simoptions1;
simoptionsTS.simperiods=10^4;
TimeSeries1=SimTimeSeriesValues_InfHorz(Policy1b,FnsToEvaluate,Params,n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,simoptionsTS);
fprintf('SimTimeSeries mean should roughly match AggVars mean (Monte Carlo) \n')
[AggVars1.assets.Mean, mean(TimeSeries1.assets(:))]
[AggVars1.earnings.Mean, mean(TimeSeries1.earnings(:))]

%% SimPanelValues: cross-validate the AgentDist statistics with a simulated panel
% SimPanelValues simulates many agents drawn from an initial distribution, where SimTimeSeries above
% follows one agent for a long time, so it is a different code path and worth checking separately.
% Starting the panel from the stationary distribution means it stays there, so pooling over agents
% and periods should reproduce the same moments.
% Note: a panel is not reproducible run to run (the parfor workers do not take the client's rng), so
% these are roughly-equal Monte Carlo checks and must never be compared as an exact zero.
simoptionsPanel=simoptions1;
simoptionsPanel.numbersims=10^4; % number of agents
simoptionsPanel.simperiods=100; % periods per agent
simoptionsPanel.burnin=0; % InitialDist is the stationary dist, so there is nothing to burn in
SimPanel1=SimPanelValues_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,[],Params,n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,simoptionsPanel);
spA=SimPanel1.assets(:); spE=SimPanel1.earnings(:); % [simperiods,numbersims] each, pooled

fprintf('SimPanelValues mean should roughly match AllStats/AggVars mean (Monte Carlo) \n')
[AllStats1.assets.Mean, AggVars1.assets.Mean, mean(spA)]
[AllStats1.earnings.Mean, AggVars1.earnings.Mean, mean(spE)]
fprintf('SimPanelValues std should roughly match AllStats std \n')
[AllStats1.assets.StdDeviation, std(spA,1)]
[AllStats1.earnings.StdDeviation, std(spE,1)]
fprintf('SimPanelValues period-1 cross-section is a straight draw from the stationary dist, so its mean should match too \n')
[AllStats1.assets.Mean, mean(SimPanel1.assets(1,:))]
[AllStats1.earnings.Mean, mean(SimPanel1.earnings(1,:))]
% Note: this model has no z and no e, so it is deterministic and the stationary distribution is a
% point mass. Every correlation below is therefore 0/0 and prints as NaN on both sides. NaN paired
% with NaN is still the two agreeing; it is the shock-bearing figs that make these checks bite.
fprintf('SimPanelValues lag-1 autocorrelation (within agent, pooled over agents) should roughly match AutoCorrTransProbs.AutoCorrelation \n')
spAx=SimPanel1.assets(1:end-1,:); spAy=SimPanel1.assets(2:end,:);
spEx=SimPanel1.earnings(1:end-1,:); spEy=SimPanel1.earnings(2:end,:);
tempPA=corrcoef(spAx(:),spAy(:)); tempPE=corrcoef(spEx(:),spEy(:));
[AutoCorrTransProbs1.assets.AutoCorrelation, tempPA(1,2)]
[AutoCorrTransProbs1.earnings.AutoCorrelation, tempPE(1,2)]
fprintf('SimPanelValues cross-section corr(earnings,assets) should roughly match CrossSectionCovarCorr \n')
tempPC=corrcoef(spE,spA);
[CrossSectionCorr1.earnings.CorrelationWith.assets, tempPC(1,2)]


%% SimPanelValues WITH the grid interpolation layer
% The panel above uses the non-GI solve, so simoptions.gridinterplayer is 0 and SimPanelIndexes_InfHorz
% only ever takes its gridinterplayer==0 branches. Running the panel on the GI solve as well is what
% reaches the GI branches, where Policy carries the extra L2index and L2flag rows. Those branches are
% much less travelled: a two-endogenous-state model with the grid interpolation layer is exactly the
% case that had a reshape using 1:end-1 where the policy has l_a+2 rows, and that then read the L2flag
% as if it were the L2 index.
AutoCorr3=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,simoptions3);
simoptionsPanel3=simoptions3;
simoptionsPanel3.numbersims=10^4;
simoptionsPanel3.simperiods=100;
simoptionsPanel3.burnin=0;
SimPanel3=SimPanelValues_InfHorz(StationaryDist3,Policy3b,FnsToEvaluate,[],Params,n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,simoptionsPanel3);
sp3A=SimPanel3.assets(:); sp3E=SimPanel3.earnings(:);

fprintf('SimPanelValues with GI, mean should roughly match AllStats/AggVars mean (Monte Carlo) \n')
[AllStats3.assets.Mean, AggVars3.assets.Mean, mean(sp3A)]
[AllStats3.earnings.Mean, AggVars3.earnings.Mean, mean(sp3E)]
fprintf('SimPanelValues with GI, std should roughly match AllStats std \n')
[AllStats3.assets.StdDeviation, std(sp3A,1)]
[AllStats3.earnings.StdDeviation, std(sp3E,1)]
fprintf('SimPanelValues with GI, period-1 cross-section mean should match the stationary dist \n')
[AllStats3.assets.Mean, mean(SimPanel3.assets(1,:))]
[AllStats3.earnings.Mean, mean(SimPanel3.earnings(1,:))]
% Note: this model has no z and no e, so it is deterministic and the stationary distribution is a
% point mass. Every correlation below is therefore 0/0 and prints as NaN on both sides. NaN paired
% with NaN is still the two agreeing; it is the shock-bearing figs that make these checks bite.
fprintf('SimPanelValues with GI, lag-1 autocorrelation should roughly match AutoCorrTransProbs (of these, the one most sensitive to the interpolation weights) \n')
sp3Ax=SimPanel3.assets(1:end-1,:); sp3Ay=SimPanel3.assets(2:end,:);
sp3Ex=SimPanel3.earnings(1:end-1,:); sp3Ey=SimPanel3.earnings(2:end,:);
tempP3A=corrcoef(sp3Ax(:),sp3Ay(:)); tempP3E=corrcoef(sp3Ex(:),sp3Ey(:));
[AutoCorr3.assets.AutoCorrelation, tempP3A(1,2)]
[AutoCorr3.earnings.AutoCorrelation, tempP3E(1,2)]


%% Plot the stationary dist over assets (with/without grid interp)
fig=figure(figure_c);
plot(a_grid_big,cumsum(sum(StationaryDist1,2)), a_grid_big,cumsum(sum(StationaryDist3,2)))
title('CDF of assets: without vs with grid interp'); legend('1','3')

%% Howards iteration
% Howards improvement iterations are just an accelerator for the value function iteration, so
% turning them off (vfoptions.howards=0, which is then pure value function iteration) must give
% the same V and Policy.
vfoptions1_noH=vfoptions1;
vfoptions1_noH.howards=0;
[V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
[V1noH,Policy1noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_noH);
fprintf('howards=0 (pure VFI), this should be zero: %2.8f \n',max(abs(V1(:)-V1noH(:))))
fprintf('howards=0 (pure VFI), this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1noH(:))))

% Same again, with the grid interpolation layer
vfoptions3_noH=vfoptions3;
vfoptions3_noH.howards=0;
[V3,Policy3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
[V3noH,Policy3noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3_noH);
fprintf('howards=0 (pure VFI, with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3noH(:))))
fprintf('howards=0 (pure VFI, with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3noH(:))))

clear V1 V3 V1noH V3noH Policy1 Policy3 Policy1noH Policy3noH

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
