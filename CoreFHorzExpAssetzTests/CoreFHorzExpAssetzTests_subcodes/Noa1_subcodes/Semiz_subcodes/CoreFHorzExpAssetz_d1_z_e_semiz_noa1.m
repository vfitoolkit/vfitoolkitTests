function output=CoreFHorzExpAssetz_d1_z_e_semiz_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% ExpAssetz noa1 + semiz: with d1 (labour), z, e. Experience asset is the only endogenous state.
% n_a is scalar (n_a_justexpasset); a_grid is the a2_grid. n_a_big/a_grid_big are unused.
% n_d=[n_d1,n_d2,n_d3]; d_grid=[d1_grid; d2_grid; d3_grid]
%
% Differs from the withA1 version:
%   - no DC/GI/DC+GI blocks (irrelevant without a1)
%   - jequaloneDist built on n_a directly (no big-grid distinction)
% PENDING TOOLKIT SUPPORT (test-first): the ExpAssetzSemiExo noa1 raws do not exist
% yet, so this errors at the ValueFnIter call.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% semiz
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
% e
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.pi_e=simoptionsbaseline.pi_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
% Experience asset (z variant)
vfoptions.experienceassetz=1;
simoptions.experienceassetz=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;
simoptions.z_grid=z_grid;

jequaloneDist=zeros([n_a,vfoptions.n_semiz,n_z,vfoptions.n_e],'gpuArray');
jequaloneDist(1,ceil(vfoptions.n_semiz/2),ceil(n_z/2),ceil(vfoptions.n_e/2))=1;

ReturnFn=@(d1,d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetz_d1_z_e_semiz_noa1(d1,d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);

FnsToEvaluate.assets=@(d1,d2,d3,a,semiz,z,e) a;
FnsToEvaluate.earnings=@(d1,d2,d3,a,semiz,z,e,w,kappa_j) w*kappa_j*d1*d2*a*semiz*z*e;

%% Basic VFI
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

% lowmemory variants
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=2;
[V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(V1(:)-V1C(:))))
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=3;
[V1D,Policy1D]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=3, this should be zero: %2.8f \n',max(abs(V1(:)-V1D(:))))
fprintf('lowmemory=3, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1D(:))))
vfoptions1.lowmemory=0;

%%
clear V1 V1B V1C V1D Policy1B Policy1C Policy1D PolicyVals1 V1fromPolicy

%% StationaryDist + EvalFn (use n_a directly -- no big-grid distinction in noa1)
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);

fprintf('noa1 moments (d1+semiz+z+e)\n')
[AllStats1.assets.Mean, AllStats1.earnings.Mean]

%% Sim panel cross-check
SimPanelValues1=SimPanelValues_FHorz_Case1(jequaloneDist,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions1);
fprintf('SimPanel mean assets (should roughly match AgeConditionalStats1):\n')
[AgeConditionalStats1.assets.Mean; mean(SimPanelValues1.assets,2)']

%% Graph
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:N_j, AgeConditionalStats1.assets.Mean)
title('Age-conditional experience asset (d1+semiz+z+e, noa1)')
subplot(2,1,2); plot(1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('Age-conditional earnings')

%% Other commands
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);

%%
output=struct();

end
