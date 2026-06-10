function output=CoreFHorzTPathExpAsset_d1_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c)

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% ExperienceAsset (a2)
vfoptions.experienceasset=vfoptionsbaseline.experienceasset;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.experienceasset=simoptionsbaseline.experienceasset;
simoptions.aprimeFn=simoptionsbaseline.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

% zero assets, mid points for any shocks
jequaloneDist_big=zeros([n_a_big,n_z],'gpuArray');
jequaloneDist_big(1,1,ceil(n_z/2))=1;
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,1,ceil(n_z/2))=1;

ReturnFn=@(d1,d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_z_noe_nosemiz(d1,d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(d1,d2,a1prime,a1,a2,z) a1+a2;
FnsToEvaluate.earnings=@(d1,d2,a1prime,a1,a2,z,w,kappa_j) w*kappa_j*d1*d2*a2*z;


% Need period T for V and Policy
V_final=zeros([n_a,n_z,N_j],'gpuArray');
Policy_final=ones([3,n_a,n_z,N_j],'gpuArray');
Policy_final_GI=ones([5,n_a,n_z,N_j],'gpuArray');
% big versions of them
V_final_big=zeros([n_a_big,n_z,N_j],'gpuArray');
Policy_final_big=ones([3,n_a_big,n_z,N_j],'gpuArray');
Policy_final_big_GI=ones([5,n_a_big,n_z,N_j],'gpuArray');

%% With vs Without fastOLG
transpathoptionsslow.fastOLG=0;
vfoptions1=vfoptions;
simoptions1=simoptions;
[VPath1slow,PolicyPath1slow]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions1);

[VPath1fast,PolicyPath1fast]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

fprintf('fastOLG, this should be zero: %2.8f \n',max(abs(VPath1slow(:)-VPath1fast(:))))
fprintf('fastOLG, this should be zero: %2.8f \n',max(abs(PolicyPath1slow(:)-PolicyPath1fast(:))))

clear VPath1fast VPath1slow PolicyPath1fast PolicyPath1slow

%% With and without divide-and-conquer (both with slowOLG)
[VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions1);

% PolicyVals1=PolicyInd2Val_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions1);

% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[VPath2,PolicyPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions2);

fprintf('Divide-and-conquer (slowOLG), this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('Divide-and-conquer (slowOLG), this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))

% lowmemory should give same answer
vfoptions1.lowmemory=1;
[VPath1B,PolicyPath1B]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions1);
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[VPath2B,PolicyPath2B]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions2);
vfoptions2.lowmemory=0;

fprintf('low memory (slowOLG), this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath1B(:))))
fprintf('low memory (slowOLG), this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath1B(:))))
fprintf('low memory (slowOLG), this should be zero: %2.8f \n',max(abs(VPath2(:)-VPath2B(:))))
fprintf('low memory (slowOLG), this should be zero: %2.8f \n',max(abs(PolicyPath2(:)-PolicyPath2B(:))))

clear VPath1 VPath2 PolicyPath1 PolicyPath2 VPath1B VPath2B PolicyPath1B PolicyPath2B % PolicyVals1


%% With and without divide-and-conquer (both with fastOLG)
[VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

% Solve with divide-and-conquer, should give same answer
[VPath2,PolicyPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);

fprintf('Divide-and-conquer (fastOLG), this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('Divide-and-conquer (fastOLG), this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))

% lowmemory should give same answer
vfoptions1.lowmemory=1;
[VPath1B,PolicyPath1B]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[VPath2B,PolicyPath2B]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
vfoptions2.lowmemory=0;

fprintf('low memory (fastOLG), this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath1B(:))))
fprintf('low memory (fastOLG), this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath1B(:))))
fprintf('low memory (fastOLG), this should be zero: %2.8f \n',max(abs(VPath2(:)-VPath2B(:))))
fprintf('low memory (fastOLG), this should be zero: %2.8f \n',max(abs(PolicyPath2(:)-PolicyPath2B(:))))

%%
clear VPath1 VPath2 PolicyPath1 PolicyPath2 VPath1B VPath2B PolicyPath1B PolicyPath2B

%% Solve with grid-interpolation. With and without divide-and-conquer (both with slowOLG)
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[VPath3,PolicyPath3]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions3);

% PolicyVals3=PolicyInd2Val_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions3);

% Solve with divide-and-conquer, should give same answer
vfoptions4=vfoptions;
vfoptions4.divideandconquer=1;
vfoptions4.gridinterplayer=1;
vfoptions4.ngridinterp=5;
simoptions4=simoptions;
simoptions4.gridinterplayer=vfoptions4.gridinterplayer;
simoptions4.ngridinterp=vfoptions4.ngridinterp;
[VPath4,PolicyPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions4);

fprintf('Divide-and-conquer (with GI, slowOLG), this should be zero: %2.8f \n',max(abs(VPath3(:)-VPath4(:))))
fprintf('Divide-and-conquer (with GI, slowOLG), this should be zero: %2.8f \n',max(abs(PolicyPath3(:)-PolicyPath4(:))))

% lowmemory should give same answer
vfoptions3.lowmemory=1;
[VPath3B,PolicyPath3B]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions3);
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[VPath4B,PolicyPath4B]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions4);
vfoptions4.lowmemory=0;

fprintf('low memory (with GI, slowOLG), this should be zero: %2.8f \n',max(abs(VPath3(:)-VPath3B(:))))
fprintf('low memory (with GI, slowOLG), this should be zero: %2.8f \n',max(abs(PolicyPath3(:)-PolicyPath3B(:))))
fprintf('low memory (with GI, slowOLG), this should be zero: %2.8f \n',max(abs(VPath4(:)-VPath4B(:))))
fprintf('low memory (with GI, slowOLG), this should be zero: %2.8f \n',max(abs(PolicyPath4(:)-PolicyPath4B(:))))

clear VPath3 VPath4 PolicyPath3 PolicyPath4 VPath3B VPath4B PolicyPath3B PolicyPath4B % PolicyVals3


%% Solve with grid-interpolation. With and without divide-and-conquer (both with fastOLG)
[VPath3,PolicyPath3]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions3);

% Solve with divide-and-conquer, should give same answer
[VPath4,PolicyPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);

fprintf('Divide-and-conquer (with GI, fastOLG), this should be zero: %2.8f \n',max(abs(VPath3(:)-VPath4(:))))
fprintf('Divide-and-conquer (with GI, fastOLG), this should be zero: %2.8f \n',max(abs(PolicyPath3(:)-PolicyPath4(:))))

% lowmemory should give same answer
vfoptions3.lowmemory=1;
[VPath3B,PolicyPath3B]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions3);
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[VPath4B,PolicyPath4B]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
vfoptions4.lowmemory=0;

fprintf('low memory (with GI, fastOLG), this should be zero: %2.8f \n',max(abs(VPath3(:)-VPath3B(:))))
fprintf('low memory (with GI, fastOLG), this should be zero: %2.8f \n',max(abs(PolicyPath3(:)-PolicyPath3B(:))))
fprintf('low memory (with GI, fastOLG), this should be zero: %2.8f \n',max(abs(VPath4(:)-VPath4B(:))))
fprintf('low memory (with GI, fastOLG), this should be zero: %2.8f \n',max(abs(PolicyPath4(:)-PolicyPath4B(:))))

%%
clear VPath3 VPath4 PolicyPath3 PolicyPath4 VPath3B VPath4B PolicyPath3B PolicyPath4B

%% Use a really big a_grid, then the moments should be essentially the same with/without grid interpolation
simoptions2.a_grid=a_grid_big; % ExpAsset codes slice a2_grid out of simoptions.a_grid using n_a1, so must match n_a_big
simoptions4.a_grid=a_grid_big;
[VPath2b,PolicyPath2b]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_big, Policy_final_big, Params, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDist_initial_big=StationaryDist_FHorz_Case1(jequaloneDist_big,AgeWeightParamNames,PolicyPath2b(:,:,:,:,:,1),n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AgentDistPath2=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_big, jequaloneDist_big, PricePath, ParamPath, PolicyPath2b, AgeWeightParamNames,n_d,n_a_big,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions2);
AggVarsPath2=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath2, PolicyPath2b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, transpathoptionsbaseline, simoptions2);

[VPath4b,PolicyPath4b]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_big, Policy_final_big_GI, Params, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
AgentDistPath4=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_big, jequaloneDist_big, PricePath, ParamPath, PolicyPath4b, AgeWeightParamNames,n_d,n_a_big,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions4);
AggVarsPath4=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath4, PolicyPath4b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, transpathoptionsbaseline, simoptions4);
simoptions2.a_grid=a_grid; % restore for the small-grid section below
simoptions4.a_grid=a_grid;

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(AgentDistPath2(:)-AgentDistPath4(:))))
[AggVarsPath2.earnings.Mean; AggVarsPath4.earnings.Mean]
[AggVarsPath2.assets.Mean; AggVarsPath4.assets.Mean]

%% Do some graphs of the AggVars path to see them
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:T,AggVarsPath2.earnings.Mean, 1:1:T,AggVarsPath4.earnings.Mean)
title('Earnings Mean')
legend('1','2')
subplot(2,1,2); plot(1:1:T,AggVarsPath2.assets.Mean, 1:1:T,AggVarsPath4.assets.Mean)
title('Assets Mean')
legend('1','2')

clear VPath2b VPath4b AggVarsPath4 PolicyPath4b AgentDistPath2 AgentDistPath4

%% If the path is all constant, should just get same answer as when we dont have a TPath
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
AgentDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions1);
PricePathConstant.r=Params.r*ones(1,T);
ParamPathConstant.sigma=Params.sigma*ones(1,T);
[VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePathConstant, ParamPathConstant, T, V1, Policy1, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
AgentDistPath1=AgentDistOnTransPath_Case1_FHorz(AgentDist1, jequaloneDist, PricePathConstant, ParamPathConstant, PolicyPath1, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions1);
V1_rep=repmat(V1,1,1,1,1,1,T);
fprintf('Do nothing TPath, this should be zero, V: %2.8f \n',max(abs(VPath1(:)-V1_rep(:))))
Policy1_rep=repmat(Policy1,1,1,1,1,1,1,T);
fprintf('Do nothing TPath, this should be zero, Policy: %2.8f \n',max(abs(PolicyPath1(:)-Policy1_rep(:))))
AgentDist1_rep=repmat(AgentDist1,1,1,1,1,1,T);
fprintf('Do nothing TPath, this should be zero, AgentDist: %2.8f \n',max(abs(AgentDistPath1(:)-AgentDist1_rep(:))))

clear V1 Policy1 VPath1 PolicyPath1

%% Run the GE transition path, but with transpathoptions.maxit=1 -- just shape-check
[~,PolicyPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyPath2(:,:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptions2);
clear PolicyPath2

transpathoptions.maxiter=1;
GeneralEqmEqns.dummy=@(earnings) 0;

transpathoptions.fastOLG=1;
PricePath2=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

transpathoptions.fastOLG=0;
PricePath2B=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

fprintf('One iter of TPath, with/without fastOLG, this should be zero: %2.8f \n',max(abs(PricePath2.r-PricePath2B.r)))

% Big grid, uses vfoptions2 with divide-and-conquer
transpathoptions.fastOLG=1;
PricePath3A=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final_big, AgentDist_initial_big, jequaloneDist_big, n_d, n_a_big, n_z, N_j, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions2, vfoptions2);

% vfoptions4 has divide-and-conquer and grid interpolation layer
PricePath3B=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final_big, AgentDist_initial_big, jequaloneDist_big, n_d, n_a_big, n_z, N_j, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions4, vfoptions4);

fprintf('One iter of TPath, with/without GI, this should be close to zero: %2.8f \n',max(abs(PricePath3A.r-PricePath3B.r)))

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
