function output=CoreFHorzExpAsset_d1_noz_noe_noa1_semiz_with2A2(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% with2A2: TWO experience assets (vfoptions.experienceasset=2), and no a1, so the two
% experience assets are the only endogenous states.
%   n_a    = [n_a2_1, n_a2_2]           (n_a_2A2_justexpasset)
%   a_grid = [a2_1_grid; a2_2_grid]
% aprimeFn takes the 'whicha' selector between the a2 inputs and the parameters:
%   aprimeFn(d2,a2_1,a2_2,whicha,phi1,phi2,phi3,phi4)
%   whicha=1 -> a2_1' = phi1*(1-d2)+(1-phi2)*a2_1        (as in the one-experience-asset tier)
%   whicha=2 -> a2_2' = phi3*d2*a2_1+(1-phi4)*a2_2       (deliberately coupled to a2_1)
% a2_2 is 'cumulated lifetime earnings' and enters the ReturnFn through the retirement pension
% (pension+pensionrate*a2_2), so V genuinely depends on it.
% Otherwise identical to the one-experience-asset version of this test.
% ExpAsset noa1 + semiz: with d1 (labour), no z, no e. The two experience assets are the only endogenous states.
% n_a is [n_a2_1,n_a2_2] (n_a_2A2_justexpasset); a_grid is [a2_1_grid;a2_2_grid]. n_a_big/a_grid_big are unused.
% n_d=[n_d1,n_d2,n_d3]; d_grid=[d1_grid; d2_grid; d3_grid]
%
% Differs from the a1 version:
%   - no DC/GI/DC+GI blocks (irrelevant without a1)
%   - jequaloneDist built on n_a directly (no big-grid distinction)
% PENDING TOOLKIT SUPPORT: ExpAsset+SemiExo+noa1 errors immediately at ValueFnIter call
% ('Have not implemented experience assets with semi-exogenous shocks, without also having a standard asset')


% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
n_z=0;
z_grid=[];
pi_z=[];
% semiz
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
% Experience asset
vfoptions.experienceasset=2;
simoptions.experienceasset=2;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn_2A2;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

jequaloneDist=zeros([n_a,vfoptions.n_semiz],'gpuArray');
jequaloneDist(1,1,ceil(vfoptions.n_semiz/2))=1;

ReturnFn=@(d1,d2,d3,a2_1,a2_2,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,pensionrate,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_noa1_semiz_with2A2(d1,d2,d3,a2_1,a2_2,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,pensionrate,uempbenefit,searcheffortcost);

FnsToEvaluate.assets=@(d1,d2,d3,a2_1,a2_2,semiz) a2_1;
FnsToEvaluate.lifetimeearnings=@(d1,d2,d3,a2_1,a2_2,semiz) a2_2;
FnsToEvaluate.earnings=@(d1,d2,d3,a2_1,a2_2,semiz,w,kappa_j) w*kappa_j*d1*d2*a2_1*semiz;

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

fprintf('noa1 moments (d1+semiz)\n')
[AllStats1.assets.Mean, AllStats1.earnings.Mean]

%% Sim panel cross-check
SimPanelValues1=SimPanelValues_FHorz_Case1(jequaloneDist,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions1);
fprintf('SimPanel mean assets (should roughly match AgeConditionalStats1):\n')
[AgeConditionalStats1.assets.Mean; mean(SimPanelValues1.assets,2)']
[AgeConditionalStats1.lifetimeearnings.Mean; mean(SimPanelValues1.lifetimeearnings,2)']

%% Graph
fig=figure(figure_c);
subplot(3,1,1); plot(1:1:N_j, AgeConditionalStats1.assets.Mean)
title('Age-conditional assets (d1+semiz, noa1)')
subplot(3,1,2); plot(1:1:N_j, AgeConditionalStats1.earnings.Mean)
title('Age-conditional earnings')
subplot(3,1,3); plot(1:1:N_j, AgeConditionalStats1.lifetimeearnings.Mean)
title('Second experience asset (lifetime earnings) mean')

%% Other commands
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions1);


%% V_Jplus1: use V of period jstar as the terminal value function of a shorter model
% Solve the model, then solve a shorter model that runs only periods 1,...,jstar-1, giving it
% vfoptions.V_Jplus1=V(:,:,:,jstar). V_Jplus1 is the value fn of period N_j+1 of the model being
% solved, so the shorter model has Njs=jstar-1 periods, and the age-dependent parameters are
% trimmed to length Njs (agej and kappa_j; the aprimeFn parameters phi1..phi4 are scalars).
% V and Policy must then be identical to the original model for periods 1,...,jstar-1.
% This tier has no a1 for divide-and-conquer or the grid interp layer to operate on, so it
% defines only vfoptions1. Run at jstar=round(3*N_j/4) and again at jstar=N_j (so the
% retirement periods, and the terminal V_Jplus1 branch, are also covered), each with the
% same lowmemory rungs as above.
% Note: mewj is age-dependent, but is only used for the agent distribution, which is not
% computed here, so it is left alone.
for jstar=[round(3*N_j/4),N_j]
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
    fprintf('V_Jplus1 (jstar=%i), lowmemory=1, this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i), lowmemory=1, this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
end

clear Vbase Policybase Vshort Policyshort

%%
output=struct();

end
