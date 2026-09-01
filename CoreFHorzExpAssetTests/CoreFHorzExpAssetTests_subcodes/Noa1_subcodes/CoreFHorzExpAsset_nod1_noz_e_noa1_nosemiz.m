function output=CoreFHorzExpAsset_nod1_noz_e_noa1_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% ExpAsset noa1: no d1, no z, e. Experience asset is the only endogenous state.
% n_a is scalar (n_a_justexpasset); a_grid is the a2_grid. n_a_big/a_grid_big are unused.
% n_d=n_d2; d_grid=d2_grid
%
% Differs from the a1 version:
%   - no DC/GI/DC+GI blocks (irrelevant without a1)
%   - jequaloneDist built on n_a directly (no big-grid distinction)

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
n_z=0;
z_grid=[];
pi_z=[];
% e
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.pi_e=simoptionsbaseline.pi_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
% Experience asset
vfoptions.experienceasset=1;
simoptions.experienceasset=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

jequaloneDist=zeros([n_a,vfoptions.n_e],'gpuArray');
jequaloneDist(1,ceil(vfoptions.n_e/2))=1;

ReturnFn=@(d2,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_e_noa1_nosemiz(d2,a,e,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate.assets=@(d2,a,e) a;
FnsToEvaluate.earnings=@(d2,a,e,w,kappa_j) w*kappa_j*d2*a*e;

%% Basic VFI
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1fromPolicy(:)-V1(:))))

% lowmemory variants
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

%%
clear V1 V1B Policy1B PolicyVals1 V1fromPolicy

%% StationaryDist + EvalFn (use n_a directly -- no big-grid distinction in noa1)
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);

fprintf('noa1 moments (e)\n')
[AllStats1.assets.Mean, AllStats1.earnings.Mean]

%% Sim panel cross-check
SimPanelValues1=SimPanelValues_FHorz_Case1(jequaloneDist,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions1);
fprintf('SimPanel mean assets (should roughly match AgeConditionalStats1):\n')
[AgeConditionalStats1.assets.Mean; mean(SimPanelValues1.assets,2)']

%% Graph
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:N_j, AgeConditionalStats1.assets.Mean)
title('Age-conditional assets (e, noa1)')
subplot(2,1,2); plot(1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('Age-conditional earnings')

%% Other commands
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);


%% V_Jplus1: use V of period jstar as the terminal value function of a shorter model
% Solve the model, then solve a shorter model that runs only periods 1,...,jstar-1, giving it
% vfoptions.V_Jplus1=V(:,:,jstar). V_Jplus1 is the value fn of period N_j+1 of the model being
% solved, so the shorter model has Njs=jstar-1 periods, and the age-dependent parameters are
% trimmed to length Njs (agej and kappa_j; the aprimeFn parameters phi1 and phi2 are scalars).
% V and Policy must then be identical to the original model for periods 1,...,jstar-1. Each of
% the four solution methods gets a different jstar (the last of them uses jstar=N_j, so that
% one of them covers the retirement periods).
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
vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1, this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1, this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, with divide-and-conquer
jstar=round(2*N_j/3);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
vfoptionsjs=vfoptions2;
vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, with DC), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, with DC), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, with grid interpolation
jstar=round(N_j/2);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
vfoptionsjs=vfoptions3;
vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, with GI), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, with GI), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, with divide-and-conquer and grid interpolation
jstar=N_j;
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
vfoptionsjs=vfoptions4;
vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, with DC+GI), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, with DC+GI), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with DC+GI), this should be zero: %.3e \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with DC+GI), this should be zero: %.3e \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, with age-dependent shocks
% pi_e_J column j is the distribution of the e realized in period j, so the shorter model is
% given columns 1:jstar (the last of these is the distribution of e in the V_Jplus1 period).
jstar=round(N_j/3);
Njs=jstar-1;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
pi_e_J=vfoptions.pi_e.*ones(1,N_j);
pi_e_J(:,1:2:N_j)=0.5*pi_e_J(:,1:2:N_j)+0.5/vfoptions.n_e; % make it genuinely age-dependent
vfoptionsjs=vfoptions1;
vfoptionsjs.pi_e=pi_e_J;
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
vfoptionsjs.pi_e=pi_e_J(:,1:jstar);
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 with age-dependent shocks (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 with age-dependent shocks (jstar=%i), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))

clear Vbase Policybase Vshort Policyshort

%%
output=struct();

end
