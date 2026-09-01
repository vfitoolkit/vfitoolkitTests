function output=CoreInfHorz_nod_z_noe_nosemiz_with2A(Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Two-endogenous-state InfHorz test (fig 11), mirroring the fig-3 (nod, z) battery.
% Model: simplified Kitao (2008), no taxes, partial equilibrium (fixed r,w).
%   Endogenous states: a (assets), e (occupation: 0=worker, 1=entrepreneur)
%   Exogenous markov z: eta (labor productivity), theta (entrepreneurial ability)
% This exercises the InfHorz two-endogenous-state (GI2A) code paths.
%
% Note: like the published Kitao replication, an entrepreneur with (near) zero assets
% has no feasible positive consumption, so V takes -Inf at those states. The toolkit's
% iterated VFI and GI handle this (agents never optimally transition there).

DF=DiscountFactorParamNames;

%% Grids and model parameters (self-contained)
n_asset=101;
n_asset_big=501; % for the with/without GI moment check
n_e=2;           % occupation: worker/entrepreneur
n_eta=3;         % labor productivity
n_theta=2;       % entrepreneurial ability

n_a=[n_asset,n_e];
n_a_big=[n_asset_big,n_e];
n_z=[n_eta,n_theta];

% Assets: cubic-spaced, puts points near zero
assetmaxfactor=60;
asset_grid=assetmaxfactor*(linspace(0,1,n_asset).^3)';
asset_grid_big=assetmaxfactor*(linspace(0,1,n_asset_big).^3)';
e_grid=[0;1]; % 0=worker, 1=entrepreneur
a_grid=[asset_grid; e_grid];
a_grid_big=[asset_grid_big; e_grid];

% Exogenous shocks
[eta_grid,pi_eta]=discretizeAR1_FarmerToda(0,0.9,0.2,n_eta);
eta_grid=exp(eta_grid);
theta_grid=[0; 1.5];               % entrepreneurial ability (0 => no production)
pi_theta=[0.9,0.1; 0.2,0.8];       % persistent
z_grid=[eta_grid; theta_grid];
pi_z=kron(pi_theta,pi_eta);        % eta varies fastest (matches n_z=[n_eta,n_theta])

n_d=0; d_grid=[];

% Kitao parameters (no taxes)
Params.r=0.04; Params.w=1.4;
Params.alpha=0.36; Params.delta=0.06; Params.upsilon=0.88;
Params.upsilon1=Params.alpha*Params.upsilon;
Params.upsilon2=(1-Params.alpha)*Params.upsilon;
Params.leverage=0.5; % max borrowing leverage (was 'd' in the Kitao code)
Params.phi=0.05;     % extra borrowing cost
% Params.sigma, Params.beta come from the shared setup

ReturnFn=@(aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi) ...
    ReturnFn_nod_z_noe_nosemiz_with2A(aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi);

% FnsToEvaluate (functions of (aprime,eprime,a,e,eta,theta))
FnsToEvaluate.assets=@(aprime,eprime,a,e,eta,theta) a;
FnsToEvaluate.entrepreneur=@(aprime,eprime,a,e,eta,theta) e; % fraction who are entrepreneurs

%% Baseline VFI
vfoptions=struct(); simoptions=struct();
vfoptions.verbose_advice=0; % 2-endo GI would otherwise sound the postGI advice
vfoptions1=vfoptions; simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions1);

PolicyVals1=PolicyInd2Val_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions1);
fprintf('with2A: ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1fromPolicy(:)-V1(:))))

%% Grid-interpolation layer (GI2A: interpolates the asset dimension)
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions3);

PolicyVals3=PolicyInd2Val_InfHorz(Policy3,n_d,n_a,n_z,d_grid,a_grid,vfoptions3);

V3fromPolicy=ValueFnFromPolicy_InfHorz(Policy3,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions3);
fprintf('with2A: ValueFnFromPolicy with grid interp, this should be zero: %.3e \n',max(abs(V3fromPolicy(:)-V3(:))))

clear V1 V3 Policy1 Policy3 PolicyVals1 PolicyVals3 V1fromPolicy V3fromPolicy

%% Big asset grid: stationary dist and moments should be essentially the same with/without grid interp
[~,Policy1b]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions1);
StationaryDist1=StationaryDist_InfHorz(Policy1b,n_d,n_a_big,n_z,pi_z,simoptions1,Params,[]);
AllStats1=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
AggVars1=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);

[~,Policy3b]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions3);
StationaryDist3=StationaryDist_InfHorz(Policy3b,n_d,n_a_big,n_z,pi_z,simoptions3,Params,[]);
AllStats3=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions3);
AggVars3=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('with2A: with/without grid interp, should get much the same moments (for big asset grid) \n')
fprintf('with2A: StationaryDist with/without grid interp, close to zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.assets.Mean,AllStats3.assets.Mean]
[AllStats1.entrepreneur.Mean,AllStats3.entrepreneur.Mean]  % fraction entrepreneurs
[AggVars1.assets.Mean,AggVars3.assets.Mean]

%% Check the remaining InfHorz dist commands run for two endogenous states
AutoCorr1=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,simoptions1);
CrossSectionCorr1=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
ValuesOnGrid1=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
ProbDensityFns1=EvalFnOnAgentDist_ProbDensityFn_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
fprintf('with2A: AutoCorrTransProbs, CrossSectionCovarCorr, ValuesOnGrid, ProbDensityFn all ran \n')

%% SimTimeSeries: a long simulated series should reproduce the AgentDist mean (Monte Carlo)
simoptionsTS=simoptions1;
simoptionsTS.simperiods=10^5;
simoptionsTS.burnin=10^3;
TimeSeries1=SimTimeSeriesValues_InfHorz(Policy1b,FnsToEvaluate,Params,n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,simoptionsTS);
fprintf('with2A: SimTimeSeries mean should roughly match AllStats/AggVars mean (Monte Carlo) \n')
[AllStats1.assets.Mean, AggVars1.assets.Mean, mean(TimeSeries1.assets(:))]
[AllStats1.entrepreneur.Mean, AggVars1.entrepreneur.Mean, mean(TimeSeries1.entrepreneur(:))]

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
spA=SimPanel1.assets(:); spE=SimPanel1.entrepreneur(:); % [simperiods,numbersims] each, pooled

fprintf('with2A: SimPanelValues mean should roughly match AllStats/AggVars mean (Monte Carlo) \n')
[AllStats1.assets.Mean, AggVars1.assets.Mean, mean(spA)]
[AllStats1.entrepreneur.Mean, AggVars1.entrepreneur.Mean, mean(spE)]
fprintf('with2A: SimPanelValues std should roughly match AllStats std \n')
[AllStats1.assets.StdDeviation, std(spA,1)]
[AllStats1.entrepreneur.StdDeviation, std(spE,1)]
fprintf('with2A: SimPanelValues period-1 cross-section is a straight draw from the stationary dist, so its mean should match too \n')
[AllStats1.assets.Mean, mean(SimPanel1.assets(1,:))]
[AllStats1.entrepreneur.Mean, mean(SimPanel1.entrepreneur(1,:))]
fprintf('with2A: SimPanelValues lag-1 autocorrelation (within agent, pooled over agents) should roughly match AutoCorrTransProbs.AutoCorrelation \n')
spAx=SimPanel1.assets(1:end-1,:); spAy=SimPanel1.assets(2:end,:);
spEx=SimPanel1.entrepreneur(1:end-1,:); spEy=SimPanel1.entrepreneur(2:end,:);
tempPA=corrcoef(spAx(:),spAy(:)); tempPE=corrcoef(spEx(:),spEy(:));
[AutoCorr1.assets.AutoCorrelation, tempPA(1,2)]
[AutoCorr1.entrepreneur.AutoCorrelation, tempPE(1,2)]
fprintf('with2A: SimPanelValues cross-section corr(entrepreneur,assets) should roughly match CrossSectionCovarCorr \n')
tempPC=corrcoef(spE,spA);
[CrossSectionCorr1.entrepreneur.CorrelationWith.assets, tempPC(1,2)]


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
sp3A=SimPanel3.assets(:); sp3E=SimPanel3.entrepreneur(:);

fprintf('with2A: SimPanelValues with GI, mean should roughly match AllStats/AggVars mean (Monte Carlo) \n')
[AllStats3.assets.Mean, AggVars3.assets.Mean, mean(sp3A)]
[AllStats3.entrepreneur.Mean, AggVars3.entrepreneur.Mean, mean(sp3E)]
fprintf('with2A: SimPanelValues with GI, std should roughly match AllStats std \n')
[AllStats3.assets.StdDeviation, std(sp3A,1)]
[AllStats3.entrepreneur.StdDeviation, std(sp3E,1)]
fprintf('with2A: SimPanelValues with GI, period-1 cross-section mean should match the stationary dist \n')
[AllStats3.assets.Mean, mean(SimPanel3.assets(1,:))]
[AllStats3.entrepreneur.Mean, mean(SimPanel3.entrepreneur(1,:))]
fprintf('with2A: SimPanelValues with GI, lag-1 autocorrelation should roughly match AutoCorrTransProbs (of these, the one most sensitive to the interpolation weights) \n')
sp3Ax=SimPanel3.assets(1:end-1,:); sp3Ay=SimPanel3.assets(2:end,:);
sp3Ex=SimPanel3.entrepreneur(1:end-1,:); sp3Ey=SimPanel3.entrepreneur(2:end,:);
tempP3A=corrcoef(sp3Ax(:),sp3Ay(:)); tempP3E=corrcoef(sp3Ex(:),sp3Ey(:));
[AutoCorr3.assets.AutoCorrelation, tempP3A(1,2)]
[AutoCorr3.entrepreneur.AutoCorrelation, tempP3E(1,2)]


%% Plot: asset distribution and entrepreneur share, with/without grid interp
fig=figure(figure_c);
% marginal over assets: StationaryDist is [n_asset,n_e,n_eta,n_theta], so sum over dimensions 2,3,4
assetdist1=sum(sum(sum(StationaryDist1,4),3),2);
assetdist3=sum(sum(sum(StationaryDist3,4),3),2);
subplot(2,1,1); plot(asset_grid_big,cumsum(assetdist1), asset_grid_big,cumsum(assetdist3))
title('with2A: CDF of assets (without vs with grid interp)'); legend('1','3')
subplot(2,1,2); plot(1:2,[AllStats1.entrepreneur.Mean,AllStats3.entrepreneur.Mean],'o')
title('with2A: entrepreneur share')

%% Howards iteration
% Howards improvement iterations are just an accelerator for the value function iteration, so
% turning them off (vfoptions.howards=0, which is then pure value function iteration) must give
% the same V and Policy.
vfoptions1_noH=vfoptions1;
vfoptions1_noH.howards=0;
[V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions1);
[V1noH,Policy1noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions1_noH);
fprintf('howards=0 (pure VFI), this should be zero: %.3e \n',max(abs(V1(:)-V1noH(:))))
fprintf('howards=0 (pure VFI), this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1noH(:))))

% Same again, with the grid interpolation layer
vfoptions3_noH=vfoptions3;
vfoptions3_noH.howards=0;
[V3,Policy3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions3);
[V3noH,Policy3noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions3_noH);
fprintf('howards=0 (pure VFI, with GI), this should be zero: %.3e \n',max(abs(V3(:)-V3noH(:))))
fprintf('howards=0 (pure VFI, with GI), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy3noH(:))))

clear V1 V3 V1noH V3noH Policy1 Policy3 Policy1noH Policy3noH

%%
output=struct();

end
