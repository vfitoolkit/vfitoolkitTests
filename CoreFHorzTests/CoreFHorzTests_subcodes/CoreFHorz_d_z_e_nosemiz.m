function output=CoreFHorz_d_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)

% Setup vfoptions and simoptions
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;
% zeros assets, mid points for any shocks
jequaloneDist=zeros(n_a_big,n_z,vfoptions.n_e,'gpuArray'); % Note: based on n_a_big, not n_a
jequaloneDist(1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1; % no assets, midpoint shock

ReturnFn=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(d,aprime,a,z,e) a;
FnsToEvaluate.earnings=@(d,aprime,a,z,e,w,kappa_j) w*kappa_j*z*e*d;


%%
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1fromPolicy(:)-V1(:))))

% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %.3e \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy2(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=0;

% lowmemory2
vfoptions1.lowmemory=2;
[V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=2, this should be zero: %.3e \n',max(abs(V1(:)-V1C(:))))
fprintf('lowmemory=2, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=2;
[V2C,Policy2C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=2 (with DC), this should be zero: %.3e \n',max(abs(V2(:)-V2C(:))))
fprintf('lowmemory=2 (with DC), this should be zero: %.3e \n',max(abs(Policy2(:)-Policy2C(:))))
vfoptions2.lowmemory=0;

%%
clear V1 V2 V1B V2B V1C V2C Policy1 Policy2 Policy1B Policy2B Policy1C Policy2C PolicyVals1 V1fromPolicy

%% Solve with grid-interpolation
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

PolicyVals3=PolicyInd2Val_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions3);

V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy with grid interp, this should be zero: %.3e \n',max(abs(V3fromPolicy(:)-V3(:))))

% Solve with divide-and-conquer, should give same answer
vfoptions4=vfoptions;
vfoptions4.divideandconquer=1;
vfoptions4.gridinterplayer=1;
vfoptions4.ngridinterp=5;
simoptions4=simoptions;
simoptions4.gridinterplayer=vfoptions4.gridinterplayer;
simoptions4.ngridinterp=vfoptions4.ngridinterp;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);

fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %.3e \n',max(abs(V3(:)-V4(:))))
fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy4(:))))

% lowmemory
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1  (with DC+GI), this should be zero: %.3e \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1  (with DC+GI), this should be zero: %.3e \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=0;

% lowmemory2
vfoptions3.lowmemory=2;
[V3C,Policy3C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=2 (with GI), this should be zero: %.3e \n',max(abs(V3(:)-V3C(:))))
fprintf('lowmemory=2 (with GI), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy3C(:))))
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=2;
[V4C,Policy4C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=2  (with DC+GI), this should be zero: %.3e \n',max(abs(V4(:)-V4C(:))))
fprintf('lowmemory=2  (with DC+GI), this should be zero: %.3e \n',max(abs(Policy4(:)-Policy4C(:))))
vfoptions4.lowmemory=0;


%%
clear V3 V4 V3B V4B V3C V4C Policy3 Policy4 Policy3B Policy4B Policy3C Policy4C PolicyVals3 V3fromPolicy

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
fprintf('StationaryDist with/without grid interp, this should be close to zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.assets.Mean,AllStats3.assets.Mean]
[AllStats1.earnings.Gini,AllStats3.earnings.Gini]
[AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean]
[AgeConditionalStats1.assets.StdDeviation; AgeConditionalStats3.assets.StdDeviation]

clear V1b V3b Policy1b Policy3b StationaryDist1 StationaryDist3

% This is also true if using divideand-conquer
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy2b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AllStats2=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);
AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);
AutoCorrTransProbs2=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,simoptions2);
CrossSectionCorr2=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);

[V4b,Policy4b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy4b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions4);
AllStats4=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);
AgeConditionalStats4=LifeCycleProfiles_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);
AutoCorrTransProbs4=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,simoptions4);
CrossSectionCorr4=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);

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


%% Sim panel and check it gives the same age conditional stats
% With and without grid interpolation layer
SimPanelValues2=SimPanelValues_FHorz_Case1(jequaloneDist,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z, simoptions2);
SimPanelValues4=SimPanelValues_FHorz_Case1(jequaloneDist,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z, simoptions4);

% Do two comparisons to age conditional stats
fprintf('Without grid interp, sim panel data should give roughly the same age conditional stats \n')
[AgeConditionalStats2.earnings.Mean; mean(SimPanelValues2.earnings,2)']
[AgeConditionalStats2.assets.Mean; mean(SimPanelValues2.assets,2)']
fprintf('With grid interp, sim panel data should give roughly the same age conditional stats \n')
[AgeConditionalStats4.earnings.Mean; mean(SimPanelValues4.earnings,2)']
[AgeConditionalStats4.assets.Mean; mean(SimPanelValues4.assets,2)']

%% Compare AllStats, AutoCorrTransProbs, CrossSectionCovarCorr to the panel data simulation
% Note: mewj is uniform, so pooling all panel observations uses the same age-weighting as the stationary dist

% AllStats vs pooled panel data
fprintf('Without grid interp, sim panel data should give roughly the same AllStats \n')
[AllStats2.earnings.Mean, mean(SimPanelValues2.earnings(:))]
[AllStats2.assets.Mean, mean(SimPanelValues2.assets(:))]
[AllStats2.earnings.StdDeviation, std(SimPanelValues2.earnings(:),1)]
[AllStats2.assets.StdDeviation, std(SimPanelValues2.assets(:),1)]
fprintf('With grid interp, sim panel data should give roughly the same AllStats \n')
[AllStats4.earnings.Mean, mean(SimPanelValues4.earnings(:))]
[AllStats4.assets.Mean, mean(SimPanelValues4.assets(:))]
[AllStats4.earnings.StdDeviation, std(SimPanelValues4.earnings(:),1)]
[AllStats4.assets.StdDeviation, std(SimPanelValues4.assets(:),1)]

% AutoCorrTransProbs vs panel data autocovariance/autocorrelation (index jj is the transition from age jj to jj+1)
PanelAutoCovar2_earnings=zeros(1,N_j-1); PanelAutoCorr2_earnings=zeros(1,N_j-1);
PanelAutoCovar2_assets=zeros(1,N_j-1); PanelAutoCorr2_assets=zeros(1,N_j-1);
PanelAutoCovar4_earnings=zeros(1,N_j-1); PanelAutoCorr4_earnings=zeros(1,N_j-1);
PanelAutoCovar4_assets=zeros(1,N_j-1); PanelAutoCorr4_assets=zeros(1,N_j-1);
for jj=1:N_j-1
    temp=cov(SimPanelValues2.earnings(jj,:)',SimPanelValues2.earnings(jj+1,:)',1);
    PanelAutoCovar2_earnings(jj)=temp(1,2);
    temp=corrcoef(SimPanelValues2.earnings(jj,:)',SimPanelValues2.earnings(jj+1,:)');
    PanelAutoCorr2_earnings(jj)=temp(1,2);
    temp=cov(SimPanelValues2.assets(jj,:)',SimPanelValues2.assets(jj+1,:)',1);
    PanelAutoCovar2_assets(jj)=temp(1,2);
    temp=corrcoef(SimPanelValues2.assets(jj,:)',SimPanelValues2.assets(jj+1,:)');
    PanelAutoCorr2_assets(jj)=temp(1,2);
    temp=cov(SimPanelValues4.earnings(jj,:)',SimPanelValues4.earnings(jj+1,:)',1);
    PanelAutoCovar4_earnings(jj)=temp(1,2);
    temp=corrcoef(SimPanelValues4.earnings(jj,:)',SimPanelValues4.earnings(jj+1,:)');
    PanelAutoCorr4_earnings(jj)=temp(1,2);
    temp=cov(SimPanelValues4.assets(jj,:)',SimPanelValues4.assets(jj+1,:)',1);
    PanelAutoCovar4_assets(jj)=temp(1,2);
    temp=corrcoef(SimPanelValues4.assets(jj,:)',SimPanelValues4.assets(jj+1,:)');
    PanelAutoCorr4_assets(jj)=temp(1,2);
end
fprintf('Without grid interp, sim panel data should give roughly the same autocovariance/autocorrelation \n')
[AutoCorrTransProbs2.earnings.AutoCovariance; PanelAutoCovar2_earnings]
[AutoCorrTransProbs2.assets.AutoCovariance; PanelAutoCovar2_assets]
[AutoCorrTransProbs2.earnings.AutoCorrelation; PanelAutoCorr2_earnings]
[AutoCorrTransProbs2.assets.AutoCorrelation; PanelAutoCorr2_assets]
fprintf('With grid interp, sim panel data should give roughly the same autocovariance/autocorrelation \n')
[AutoCorrTransProbs4.earnings.AutoCovariance; PanelAutoCovar4_earnings]
[AutoCorrTransProbs4.assets.AutoCovariance; PanelAutoCovar4_assets]
[AutoCorrTransProbs4.earnings.AutoCorrelation; PanelAutoCorr4_earnings]
[AutoCorrTransProbs4.assets.AutoCorrelation; PanelAutoCorr4_assets]

% CrossSectionCovarCorr vs pooled panel data
fprintf('Without grid interp, sim panel data should give roughly the same cross-section covariance/correlation \n')
temp=cov(SimPanelValues2.earnings(:),SimPanelValues2.assets(:),1);
[CrossSectionCorr2.earnings.CovarianceWith.assets, temp(1,2)]
temp=corrcoef(SimPanelValues2.earnings(:),SimPanelValues2.assets(:));
[CrossSectionCorr2.earnings.CorrelationWith.assets, temp(1,2)]
fprintf('With grid interp, sim panel data should give roughly the same cross-section covariance/correlation \n')
temp=cov(SimPanelValues4.earnings(:),SimPanelValues4.assets(:),1);
[CrossSectionCorr4.earnings.CovarianceWith.assets, temp(1,2)]
temp=corrcoef(SimPanelValues4.earnings(:),SimPanelValues4.assets(:));
[CrossSectionCorr4.earnings.CorrelationWith.assets, temp(1,2)]


%% Check the various other commands run without issue
% with grid options is likely a touch trickier
vfoptions5=vfoptions;
vfoptions5.gridinterplayer=1;
vfoptions5.ngridinterp=5;
simoptions5=simoptions;
simoptions5.gridinterplayer=vfoptions5.gridinterplayer;
simoptions5.ngridinterp=vfoptions5.ngridinterp;
[V5,Policy5]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions5);
jequaloneDist5=zeros(n_a,n_z,vfoptions.n_e,'gpuArray'); % small-grid init for Policy5
jequaloneDist5(1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1; % no assets, midpoint shock
StationaryDist5=StationaryDist_FHorz_Case1(jequaloneDist5,AgeWeightParamNames,Policy5,n_d,n_a,n_z,N_j,pi_z,Params,simoptions5);
% AllStats and LifeCycleProfiles were already used
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist5,Policy5,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions5);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy5, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions5);


%% V_Jplus1: use V of period jstar as the terminal value function of a shorter model
% Solve the model, then solve a shorter model that runs only periods 1,...,jstar-1, giving it
% vfoptions.V_Jplus1=V(:,:,:,jstar). V_Jplus1 is the value fn of period N_j+1 of the model being
% solved, so the shorter model has Njs=jstar-1 periods, and the age-dependent parameters are
% trimmed to length Njs. V and Policy must then be identical to the original model for periods
% 1,...,jstar-1. Each of the four solution methods gets a different jstar (the last of them uses
% jstar=N_j, so that one of them covers the retirement periods).
% Note: mewj is age-dependent, but is only used for the agent distribution, which is not
% computed here, so it is left alone.

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
fprintf('V_Jplus1 (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1, this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1, this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=2;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=2, this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=2, this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
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
fprintf('V_Jplus1 (jstar=%i, with DC), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, with DC), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=2;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=2 (with DC), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=2 (with DC), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
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
fprintf('V_Jplus1 (jstar=%i, with GI), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, with GI), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=2;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=2 (with GI), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=2 (with GI), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
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
fprintf('V_Jplus1 (jstar=%i, with DC+GI), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, with DC+GI), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with DC+GI), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with DC+GI), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=2;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=2 (with DC+GI), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=2 (with DC+GI), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, with age-dependent shocks
% pi_z_J slice j is the transition from period j to period j+1, so the shorter model is given
% slices 1:Njs (the last of these is the transition into the V_Jplus1 period).
% pi_e_J column j is the distribution of the e realized in period j, so the shorter model is
% given columns 1:jstar (the last of these is the distribution of e in the V_Jplus1 period).
jstar=round(N_j/3);
Njs=jstar-1;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
pi_z_J=pi_z.*ones(1,1,N_j);
pi_z_J(:,:,1:2:N_j)=0.5*pi_z_J(:,:,1:2:N_j)+0.5*eye(n_z); % make it genuinely age-dependent
pi_e_J=vfoptions.pi_e.*ones(1,N_j);
pi_e_J(:,1:2:N_j)=0.5*pi_e_J(:,1:2:N_j)+0.5/vfoptions.n_e; % make it genuinely age-dependent
vfoptionsjs=vfoptions1;
vfoptionsjs.pi_e=pi_e_J;
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Vbase(:,:,:,jstar);
vfoptionsjs.pi_e=pi_e_J(:,1:jstar);
Vbase=Vbase(:,:,:,1:Njs);
Policybase=Policybase(:,:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z_J(:,:,1:Njs),ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 with age-dependent shocks (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 with age-dependent shocks (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))

clear Vbase Policybase Vshort Policyshort


%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end