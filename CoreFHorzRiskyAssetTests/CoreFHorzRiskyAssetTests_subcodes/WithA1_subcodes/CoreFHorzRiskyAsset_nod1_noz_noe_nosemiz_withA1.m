function output=CoreFHorzRiskyAsset_nod1_noz_noe_nosemiz_withA1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% RiskyAsset, with a1 + a2 (a1=safe, a2=risky). Variant: nod1_noz_noe_nosemiz.
% NOTE: GI and DC+GI blocks are PENDING toolkit support for RiskyAsset and will currently error at runtime.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
n_z=0;
z_grid=[];
pi_z=[];

% Riskyasset
vfoptions.riskyasset=1;
vfoptions.refine_d=[0,1,1];
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

jequaloneDist=zeros([n_a_big],'gpuArray');
jequaloneDist(1,1)=1;

ReturnFn=@(savings,a1prime,a1,a2,r,w,kappa_j,sigma,r_a1,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz_withA1(savings,a1prime,a1,a2,r,w,kappa_j,sigma,r_a1,agej,Jr,pension);

FnsToEvaluate.a1=@(savings,a1prime,a1,a2) a1;
FnsToEvaluate.a2=@(savings,a1prime,a1,a2) a2;
FnsToEvaluate.savings=@(savings,a1prime,a1,a2) savings;
FnsToEvaluate.earnings=@(savings,a1prime,a1,a2,w,kappa_j) w*kappa_j;

%% Basic VFI
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

%% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))

%%
clear V1 V2 Policy1 Policy2 PolicyVals1 V1fromPolicy
%% Solve with grid-interpolation (PENDING toolkit support for RiskyAsset)
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

%% Solve with divide-and-conquer + grid interp (PENDING toolkit support for RiskyAsset)
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

%%
clear V3 V4 Policy3 Policy4 PolicyVals3 V3fromPolicy
%% Big a_grid for accurate StationaryDist moments
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
[AllStats1.a1.Mean,AllStats3.a1.Mean]
[AllStats1.a2.Mean,AllStats3.a2.Mean]
[AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean]
[AgeConditionalStats1.a1.StdDeviation; AgeConditionalStats3.a1.StdDeviation]

clear V1b V3b Policy1b Policy3b StationaryDist1 StationaryDist3

% Also with DC
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

fprintf('With/without grid interp (with DC), should get much the same moments \n')
[AllStats2.a1.Mean,AllStats4.a1.Mean]
[AllStats2.a2.Mean,AllStats4.a2.Mean]
[AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean]
[AgeConditionalStats2.a1.StdDeviation; AgeConditionalStats4.a1.StdDeviation]

clear V2b V4b StationaryDist2 StationaryDist4

%% Age-conditional plots
fig=figure(figure_c);
subplot(3,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('Earnings Mean'); legend('1','2','3','4')
subplot(3,1,2); plot(1:1:N_j,AgeConditionalStats1.a1.StdDeviation, 1:1:N_j,AgeConditionalStats2.a1.StdDeviation, 1:1:N_j,AgeConditionalStats3.a1.StdDeviation, 1:1:N_j,AgeConditionalStats4.a1.StdDeviation)
title('a1 (safe) Std Dev'); legend('1','2','3','4')
subplot(3,1,3); plot(1:1:N_j,AgeConditionalStats1.a2.Mean, 1:1:N_j,AgeConditionalStats2.a2.Mean, 1:1:N_j,AgeConditionalStats3.a2.Mean, 1:1:N_j,AgeConditionalStats4.a2.Mean)
title('a2 (risky) Mean'); legend('1','2','3','4')

%% Sim panel cross-check
SimPanelValues2=SimPanelValues_FHorz_Case1(jequaloneDist,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,simoptions2);
SimPanelValues4=SimPanelValues_FHorz_Case1(jequaloneDist,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,simoptions4);

fprintf('Without grid interp, sim panel data should give roughly the same age conditional stats \n')
[AgeConditionalStats2.earnings.Mean; mean(SimPanelValues2.earnings,2)']
[AgeConditionalStats2.a1.Mean; mean(SimPanelValues2.a1,2)']
[AgeConditionalStats2.a2.Mean; mean(SimPanelValues2.a2,2)']
fprintf('With grid interp, sim panel data should give roughly the same age conditional stats \n')
[AgeConditionalStats4.earnings.Mean; mean(SimPanelValues4.earnings,2)']
[AgeConditionalStats4.a1.Mean; mean(SimPanelValues4.a1,2)']
[AgeConditionalStats4.a2.Mean; mean(SimPanelValues4.a2,2)']

%% Check the various other commands run without issue
vfoptions5=vfoptions;
vfoptions5.gridinterplayer=1;
vfoptions5.ngridinterp=5;
simoptions5=simoptions;
simoptions5.gridinterplayer=vfoptions5.gridinterplayer;
simoptions5.ngridinterp=vfoptions5.ngridinterp;
[V5,Policy5]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions5);
jequaloneDist5=zeros([n_a],'gpuArray');
jequaloneDist5(1,1)=1;
StationaryDist5=StationaryDist_FHorz_Case1(jequaloneDist5,AgeWeightParamNames,Policy5,n_d,n_a,n_z,N_j,pi_z,Params,simoptions5);
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist5,Policy5,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions5);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy5,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions5);

output=struct();

end
