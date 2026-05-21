function output=CoreFHorzExpAssetU_d1_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% Do the current setup
n_z=0;
z_grid=[];
pi_z=[];
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.pi_e=simoptionsbaseline.pi_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
% zeros assets, mid points for any shocks
jequaloneDist=zeros([n_a_big,vfoptions.n_e],'gpuArray'); % Note: based on n_a_big, not n_a
jequaloneDist(1,1,ceil(vfoptions.n_e/2))=1;

ReturnFn=@(d1,d2,a1prime,a1,a2,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_noz_e_nosemiz(d1,d2,a1prime,a1,a2,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(d1,d2,a1prime,a1,a2,e) a1;
FnsToEvaluate.humancapital=@(d1,d2,a1prime,a1,a2,e) a2;
FnsToEvaluate.earnings=@(d1,d2,a1prime,a1,a2,e,w,kappa_j) w*kappa_j*d1*d2*a2*e;

% Experience asset u
vfoptions.experienceassetu=1;
simoptions.experienceassetu=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;
vfoptions.n_u=vfoptionsbaseline.n_u;
vfoptions.u_grid=vfoptionsbaseline.u_grid;
vfoptions.pi_u=vfoptionsbaseline.pi_u;
simoptions.n_u=vfoptions.n_u;
simoptions.u_grid=vfoptions.u_grid;
simoptions.pi_u=vfoptions.pi_u;


%%
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))

% 
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=0;

%%
clear V1 V2 V1B V2B Policy1 Policy2 Policy1B Policy2B PolicyVals1 V1fromPolicy
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
fprintf('ValueFnFromPolicy with grid interp, this should be zero: %2.8f \n',max(abs(V3fromPolicy(:)-V3(:))))

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

% 
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=0;

%%
clear V3 V4 V3B V4B Policy3 Policy4 Policy3B Policy4B PolicyVals3 V3fromPolicy
%% Use a really big a_grid, then the moments should be essentially the same with/without grid interpolation

simoptions1.a_grid=a_grid_big;
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

simoptions3.a_grid=a_grid_big;
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

% This is also true if using divideand-conquer
simoptions2.a_grid=a_grid_big;
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy2b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AllStats2=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);
AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);

simoptions4.a_grid=a_grid_big;
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
subplot(3,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('Earnings Mean')
legend('1','2','3','4')
subplot(3,1,2); plot(1:1:N_j,AgeConditionalStats1.assets.StdDeviation, 1:1:N_j,AgeConditionalStats2.assets.StdDeviation, 1:1:N_j,AgeConditionalStats3.assets.StdDeviation, 1:1:N_j,AgeConditionalStats4.assets.StdDeviation)
title('Assets Std Dev')
legend('1','2','3','4')
subplot(3,1,3); plot(1:1:N_j,AgeConditionalStats1.humancapital.Mean, 1:1:N_j,AgeConditionalStats2.humancapital.Mean, 1:1:N_j,AgeConditionalStats3.humancapital.Mean, 1:1:N_j,AgeConditionalStats4.humancapital.Mean)
title('Human capital mean')
legend('1','2','3','4')
% ylim([0,0.01]) % If you want to make graph look nicer


%% Sim panel and check it gives the same age conditional stats
% With and without grid interpolation layer
SimPanelValues2=SimPanelValues_FHorz_Case1(jequaloneDist,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z, simoptions2);
SimPanelValues4=SimPanelValues_FHorz_Case1(jequaloneDist,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z, simoptions4);

% Do two comparisons to age conditional stats
fprintf('Without grid interp, sim panel data should give roughly the same age conditional stats \n')
[AgeConditionalStats2.earnings.Mean; mean(SimPanelValues2.earnings,2)']
[AgeConditionalStats2.assets.Mean; mean(SimPanelValues2.assets,2)']
[AgeConditionalStats2.humancapital.Mean; mean(SimPanelValues2.humancapital,2)']
fprintf('With grid interp, sim panel data should give roughly the same age conditional stats \n')
[AgeConditionalStats4.earnings.Mean; mean(SimPanelValues4.earnings,2)']
[AgeConditionalStats4.assets.Mean; mean(SimPanelValues4.assets,2)']
[AgeConditionalStats4.humancapital.Mean; mean(SimPanelValues4.humancapital,2)']

%% Check the various other commands run without issue
% with grid options is likely a touch trickier
vfoptions5=vfoptions;
vfoptions5.gridinterplayer=1;
vfoptions5.ngridinterp=5;
simoptions5=simoptions;
simoptions5.gridinterplayer=vfoptions5.gridinterplayer;
simoptions5.ngridinterp=vfoptions5.ngridinterp;
[V5,Policy5]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions5);
jequaloneDist5=zeros([n_a,vfoptions.n_e],'gpuArray'); % small-grid init for Policy5
jequaloneDist5(1,1,ceil(vfoptions.n_e/2))=1;
StationaryDist5=StationaryDist_FHorz_Case1(jequaloneDist5,AgeWeightParamNames,Policy5,n_d,n_a,n_z,N_j,pi_z,Params,simoptions5);
% AllStats and LifeCycleProfiles were already used
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist5,Policy5,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions5);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy5, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions5);

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end