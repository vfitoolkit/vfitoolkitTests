function output=CoreFHorzRiskyAsset_nod1_z_noe_semiz_withA1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% RiskyAsset, with a1 + a2 (a1=safe, a2=risky). Variant: nod1_z_noe_semiz.
% GI and DC+GI are supported for RiskyAsset (grid-interp VFI + StationaryDist); confirm the printed diffs with a MATLAB run.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;

% Riskyasset
vfoptions.riskyasset=1;
vfoptions.refine_d=[0,1,1,1];
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

jequaloneDist=zeros([n_a_big,vfoptions.n_semiz,n_z],'gpuArray');
jequaloneDist(1,1,ceil(vfoptions.n_semiz/2),ceil(n_z/2))=1;

ReturnFn=@(savings,dsemiz,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);

% With riskyasset, need to include all d in FnsToEvaluate, even if they are not in the ReturnFn
FnsToEvaluate.a1=@(riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,z) a1;
FnsToEvaluate.a2=@(riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,z) a2;
FnsToEvaluate.savings=@(riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,z) savings;
FnsToEvaluate.earnings=@(riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,z,w,kappa_j) w*kappa_j;

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

% lowmemory tests
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=2;
[V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(V1(:)-V1C(:))))
fprintf('lowmemory=2, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=0;

% lowmemory tests
vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=2;
[V2C,Policy2C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=2 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2C(:))))
fprintf('lowmemory=2 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2C(:))))
vfoptions2.lowmemory=0;

%%
clear V1 V2 V1B V2B V1C V2C Policy1 Policy2 Policy1B Policy2B Policy1C Policy2C PolicyVals1 V1fromPolicy
%% Solve with grid-interpolation (a valid, more-accurate solution)
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

%% Solve with divide-and-conquer + grid interp (should match grid-interp)
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

% lowmemory tests
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=2;
[V3C,Policy3C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=2 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3C(:))))
fprintf('lowmemory=2 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3C(:))))
vfoptions3.lowmemory=0;

% lowmemory tests
vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=2;
[V4C,Policy4C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=2  (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4C(:))))
fprintf('lowmemory=2  (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4C(:))))
vfoptions4.lowmemory=0;

%%
clear V3 V4 V3B V4B V3C V4C Policy3 Policy4 Policy3B Policy4B Policy3C Policy4C PolicyVals3 V3fromPolicy
%% V_Jplus1: use V of period jstar as the terminal value function of a shorter model
% Solve the model, then solve a shorter model that runs only periods 1,...,jstar-1, giving it
% vfoptions.V_Jplus1=Vbase(:,:,:,:,jstar). V_Jplus1 is the value fn of period N_j+1 of the model
% being solved, so the shorter model has Njs=jstar-1 periods, and the age-dependent parameters
% (agej, kappa_j) are trimmed to length Njs. V and Policy must then be identical to the original
% model for periods 1,...,Njs. Each of the four solution methods gets a different jstar (the last
% of them uses jstar=N_j, so the V_Jplus1 terminal branches cover the retirement periods); each
% method repeats the lowmemory rungs above. vfoptionsjs inherits all the riskyasset and semiz
% settings (aprimeFn/n_u/u_grid/pi_u, n_semiz/semiz_grid/SemiExoStateFn); the semiz setup depends
% only on the constant probfindjob/problosejob, so nothing else J-dependent needs trimming.
% Note: mewj is age-dependent, but is only used for the agent distribution, so is left alone.

% V_Jplus1, without divide-and-conquer, without grid interpolation
jstar=round(3*N_j/4);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptionsjs=vfoptions1;
vfoptionsjs.V_Jplus1=Vbase(:,:,:,:,jstar);
Vbase=Vbase(:,:,:,:,1:Njs);
Policybase=Policybase(:,:,:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1, this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1, this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=2;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=2, this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=2, this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

% V_Jplus1, with divide-and-conquer
jstar=round(2*N_j/3);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
vfoptionsjs=vfoptions2;
vfoptionsjs.V_Jplus1=Vbase(:,:,:,:,jstar);
Vbase=Vbase(:,:,:,:,1:Njs);
Policybase=Policybase(:,:,:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, with DC), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, with DC), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=2;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=2 (with DC), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=2 (with DC), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

% V_Jplus1, with grid interpolation
jstar=round(N_j/2);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
vfoptionsjs=vfoptions3;
vfoptionsjs.V_Jplus1=Vbase(:,:,:,:,jstar);
Vbase=Vbase(:,:,:,:,1:Njs);
Policybase=Policybase(:,:,:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, with GI), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, with GI), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=2;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=2 (with GI), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=2 (with GI), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

% V_Jplus1, with divide-and-conquer and grid interpolation
jstar=N_j;
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
vfoptionsjs=vfoptions4;
vfoptionsjs.V_Jplus1=Vbase(:,:,:,:,jstar);
Vbase=Vbase(:,:,:,:,1:Njs);
Policybase=Policybase(:,:,:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, with DC+GI), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, with DC+GI), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (with DC+GI), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=2;
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=2 (with DC+GI), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=2 (with DC+GI), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
vfoptionsjs.lowmemory=0;

clear Vbase Vshort Policybase Policyshort Paramsjs vfoptionsjs
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
jequaloneDist5=zeros([n_a,vfoptions.n_semiz,n_z],'gpuArray');
jequaloneDist5(1,1,ceil(vfoptions.n_semiz/2),ceil(n_z/2))=1;
StationaryDist5=StationaryDist_FHorz_Case1(jequaloneDist5,AgeWeightParamNames,Policy5,n_d,n_a,n_z,N_j,pi_z,Params,simoptions5);
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist5,Policy5,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions5);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy5,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions5);

output=struct();

end
