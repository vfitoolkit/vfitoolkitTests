function output=CoreFHorzExpAssetsemiz_nod1_noz_noe_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% experienceassetsemiz noa1: aprime depends on (d2,a2,semiz). semiz always present. No ordinary z, no e.
% The experience asset a2 is the ONLY endogenous state.
% n_a is scalar (n_a_justexpasset); a_grid is the a2_grid. n_a_big/a_grid_big are unused.
% n_d input = [n_d2, n_d3]. d_grid = [d2_grid; d3_grid].
%
% Differs from the withA1 version:
%   - no DC/GI/DC+GI blocks (irrelevant without a1)
%   - jequaloneDist built on n_a directly (no big-grid distinction)

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
% Experience asset (semiz variant)
vfoptions.experienceassetsemiz=1;
simoptions.experienceassetsemiz=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

% zero experience asset, mid point for semiz
jequaloneDist=zeros([n_a,vfoptions.n_semiz],'gpuArray');
jequaloneDist(1,ceil(vfoptions.n_semiz/2))=1;

ReturnFn=@(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetsemiz_nod1_noz_noe_noa1(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% Setup some FnsToEvaluate
FnsToEvaluate.humancapital=@(d2,d3,a,semiz) a;
FnsToEvaluate.earnings=@(d2,d3,a,semiz,w,kappa_j) w*kappa_j*d2*a*semiz;

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
vfoptions1.lowmemory=0;

%%
clear V1 V1B Policy1B PolicyVals1 V1fromPolicy

%% StationaryDist + EvalFn (use n_a directly -- no big-grid distinction in noa1)
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);

fprintf('noa1 moments (semiz, noz, noe)\n')
[AllStats1.humancapital.Mean, AllStats1.earnings.Mean]

%% Sim panel cross-check
SimPanelValues1=SimPanelValues_FHorz_Case1(jequaloneDist,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions1);
fprintf('SimPanel means should roughly match the age-conditional stats:\n')
[AgeConditionalStats1.humancapital.Mean; mean(SimPanelValues1.humancapital,2)']

%% Graph
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:N_j, AgeConditionalStats1.humancapital.Mean)
title('Age-conditional human capital (experienceassetsemiz, noa1, nod1, noz, noe)')
subplot(2,1,2); plot(1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('Age-conditional earnings')

%% Other commands
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);

%%
output=struct();

end
