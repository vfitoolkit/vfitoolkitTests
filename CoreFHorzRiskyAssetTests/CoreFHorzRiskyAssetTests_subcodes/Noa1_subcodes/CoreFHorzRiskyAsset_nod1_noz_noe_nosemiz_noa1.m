function output=CoreFHorzRiskyAsset_nod1_noz_noe_nosemiz_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% RiskyAsset, no d1, no z, no e, no semiz.
% n_d=[n_d2,n_d3] = [riskyshare, savings];   d_grid=[d2_grid; d3_grid]

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
n_z=0;
z_grid=[];
pi_z=[];

% Riskyasset
vfoptions.riskyasset=1;
vfoptions.refine_d=[0,1,1]; % no d1, 1 d2, 1 d3
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

% Initial dist (concentrate mass at a=0)
jequaloneDist=zeros(n_a_big,1,'gpuArray');
jequaloneDist(1)=1;

ReturnFn=@(savings,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(savings,a,r,w,kappa_j,sigma,agej,Jr,pension);

% With riskyasset, need to include all d in FnsToEvaluate, even if they are not in the ReturnFn
FnsToEvaluate.assets=@(riskyshare,savings,a) a;
FnsToEvaluate.savings=@(riskyshare,savings,a) savings;
FnsToEvaluate.earnings=@(riskyshare,savings,a,w,kappa_j) w*kappa_j;

%% Basic VFI
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1fromPolicy(:)-V1(:))))
% no lowmemory here: with none of z/e/semiz, only lowmemory=0 exists

%%
clear V1 Policy1 PolicyVals1 V1fromPolicy

%% Use a really big a_grid for more accurate moments
% Riskyasset doesn't support DC or grid interpolation layer, so there's just the one path here.
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
plot(1:1:N_j, AgeConditionalStats1.assets.Mean, 1:1:N_j, AgeConditionalStats1.savings.Mean)
title('Age-conditional mean: assets / savings')
legend('assets','savings')

%% Sim panel cross-check
SimPanelValues1=SimPanelValues_FHorz_Case1(jequaloneDist,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,simoptions1);
fprintf('SimPanel mean assets by age (should roughly match AgeConditionalStats1):\n')
[AgeConditionalStats1.assets.Mean; mean(SimPanelValues1.assets,2)']

%% Check the various other commands run without issue
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

%% V_Jplus1: use V of period jstar as the terminal value function of a shorter model
% Solve the model, then solve a shorter model that runs only periods 1,...,jstar-1, giving it
% vfoptions.V_Jplus1=Vbase(:,jstar). V_Jplus1 is the value fn of period N_j+1 of the model being
% solved, so the shorter model has Njs=jstar-1 periods, and the age-dependent parameters are
% trimmed to length Njs. V and Policy must then be identical to the original model for periods
% 1,...,jstar-1. Run at jstar=round(3*N_j/4) and again at jstar=N_j (so the retirement periods,
% and the terminal V_Jplus1 branch, are also covered).
% Note: mewj is age-dependent, but is only used for the agent distribution, so it is left alone.
for jstar=[round(3*N_j/4),N_j]
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    [Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    vfoptionsjs=vfoptions1; % inherit all the riskyasset settings
    vfoptionsjs.V_Jplus1=Vbase(:,jstar);
    Vbase=Vbase(:,1:Njs);
    Policybase=Policybase(:,:,1:Njs);
    [Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
    % no lowmemory here: with none of z/e/semiz, only lowmemory=0 exists
end
clear Vbase Vshort Policybase Policyshort Paramsjs vfoptionsjs

%%
output=struct(); % Not currently used for anything.

end
