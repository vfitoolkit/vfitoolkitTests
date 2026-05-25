function output=CoreFHorzRiskyAsset_d1_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% RiskyAsset, with d1 (labour), z, no e.
% n_d=[n_d1,n_d2,n_d3] = [labour, riskyshare, savings];   d_grid=[d1_grid; d2_grid; d3_grid]

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% Riskyasset
vfoptions.riskyasset=1;
vfoptions.refine_d=[1,1,1]; % 1 d1 (labour), 1 d2 (riskyshare), 1 d3 (savings)
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
simoptions.a_grid=a_grid;
simoptions.d_grid=d_grid;

jequaloneDist=zeros([n_a_big,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

ReturnFn=@(h,savings,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d1_z_noe_nosemiz(h,savings,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);

FnsToEvaluate.assets=@(h,savings,a,z) a;
FnsToEvaluate.savings=@(h,savings,a,z) savings;
FnsToEvaluate.hours=@(h,savings,a,z) h;
FnsToEvaluate.earnings=@(h,savings,a,z,w,kappa_j) w*kappa_j*h*z;

%% Basic VFI
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

%%
clear V1 V1B Policy1 Policy1B PolicyVals1 V1fromPolicy
%% Big a_grid for moments
simoptions1.a_grid=a_grid_big;
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

fprintf('Big a_grid moments\n')
[AllStats1.assets.Mean, AllStats1.savings.Mean, AllStats1.hours.Mean, AllStats1.earnings.Mean]

%% Graph
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:N_j, AgeConditionalStats1.assets.Mean, 1:1:N_j, AgeConditionalStats1.savings.Mean)
title('Age-conditional mean: assets / savings (d1+z)')
legend('assets','savings')
subplot(2,1,2); plot(1:1:N_j, AgeConditionalStats1.hours.Mean, 1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('Age-conditional mean: hours / earnings')
legend('hours','earnings')

%% Sim panel cross-check
SimPanelValues1=SimPanelValues_FHorz_Case1(jequaloneDist,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,simoptions1);
fprintf('SimPanel mean assets by age (should roughly match AgeConditionalStats1):\n')
[AgeConditionalStats1.assets.Mean; mean(SimPanelValues1.assets,2)']

%% Check the various other commands run without issue
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

%%
output=struct();

end
