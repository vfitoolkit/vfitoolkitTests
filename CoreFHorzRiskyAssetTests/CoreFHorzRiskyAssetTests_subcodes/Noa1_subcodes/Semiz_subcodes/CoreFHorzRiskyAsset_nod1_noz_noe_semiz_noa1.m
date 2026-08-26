function output=CoreFHorzRiskyAsset_nod1_noz_noe_semiz_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% RiskyAsset + semiz, no d1, no z, no e.
% n_d=[n_d2,n_d3,n_d4] = [riskyshare, savings, dsemiz];   d_grid=[d2_grid; d3_grid; d4_grid]

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
n_z=0;
z_grid=[];
pi_z=[];

% Semiz
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptions.d_grid=d_grid;

% Riskyasset
vfoptions.riskyasset=1;
vfoptions.refine_d=[0,1,1,1]; % no d1, 1 d2 (riskyshare), 1 d3 (savings), 1 d4 (dsemiz)
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

% Initial dist (mass at a=0, midpoint of semiz)
jequaloneDist=zeros(n_a_big,vfoptions.n_semiz,'gpuArray');
jequaloneDist(1,ceil(vfoptions.n_semiz/2))=1;

ReturnFn=@(savings,dsemiz,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% With riskyasset, need to include all d in FnsToEvaluate, even if they are not in the ReturnFn
FnsToEvaluate.assets=@(riskyshare,savings,dsemiz,a,semiz) a;
FnsToEvaluate.savings=@(riskyshare,savings,dsemiz,a,semiz) savings;
FnsToEvaluate.earnings=@(riskyshare,savings,dsemiz,a,semiz,w,kappa_j) w*kappa_j*semiz;

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
%% V_Jplus1: use V of period jstar as the terminal value function of a shorter model
% Solve the model, then solve a shorter model that runs only periods 1,...,jstar-1, giving it
% vfoptions.V_Jplus1=Vbase(:,:,jstar). V_Jplus1 is the value fn of period N_j+1 of the model
% being solved, so the shorter model has Njs=jstar-1 periods, and the age-dependent parameters
% (agej, kappa_j) are trimmed to length Njs. V and Policy must then be identical to the original
% model for periods 1,...,Njs. Two legs: jstar=round(3*N_j/4), then jstar=N_j (so the V_Jplus1
% terminal branches cover the retirement periods); each leg repeats the lowmemory rungs above.
% Note: mewj is age-dependent, but is only used for the agent distribution, so is left alone.

% First leg: jstar=round(3*N_j/4)
jstar=round(3*N_j/4);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptionsjs=vfoptions1;
vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i), lowmemory=1, this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i), lowmemory=1, this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

% Second leg: jstar=N_j (the V_Jplus1 terminal branches now cover the retirement periods)
jstar=N_j;
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptionsjs=vfoptions1;
vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i), lowmemory=1, this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i), lowmemory=1, this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

clear Vbase Vshort Policybase Policyshort

%% Big a_grid for moments
simoptions1.a_grid=a_grid_big;
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

fprintf('Big a_grid moments\n')
[AllStats1.assets.Mean, AllStats1.savings.Mean, AllStats1.earnings.Mean]
fprintf('Age-conditional assets.Mean:\n')
AgeConditionalStats1.assets.Mean

%% Graph
fig=figure(figure_c);
plot(1:1:N_j, AgeConditionalStats1.assets.Mean, 1:1:N_j, AgeConditionalStats1.savings.Mean, 1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('Age-conditional mean: assets / savings / earnings (with semiz)')
legend('assets','savings','earnings')

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
