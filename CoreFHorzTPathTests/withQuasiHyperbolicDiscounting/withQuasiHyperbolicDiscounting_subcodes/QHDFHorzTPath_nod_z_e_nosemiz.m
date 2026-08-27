function output=QHDFHorzTPath_nod_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c)

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% Do the current setup
n_d=0;
d_grid=[];

% iid shocks
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;

% zeros assets, mid points for any shocks
jequaloneDist_big=zeros([n_a_big,n_z,vfoptions.n_e],'gpuArray'); % Note: based on n_a_big
jequaloneDist_big(1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;
jequaloneDist=zeros([n_a,n_z,vfoptions.n_e],'gpuArray'); % Note: based on n_a
jequaloneDist(1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;

ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(aprime,a,z,e) a;
FnsToEvaluate.earnings=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;

%% Quasi-Hyperbolic Discounting common setup
vfoptions.exoticpreferences='QuasiHyperbolic';
vfoptions.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;

%% ============================================================
%% Naive
%% ============================================================
vfoptions.quasi_hyperbolic='Naive';

% Steady-state solve to obtain V_final (=Valt_final) and Policy_final for the TPath
vfoptions1=vfoptions;
simoptions1=simoptions;
[~,Policy_finalA,Valt_finalA,Policyalt_finalA]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

vfoptions2=vfoptions; vfoptions2.divideandconquer=1;
simoptions2=simoptions;

vfoptions3=vfoptions; vfoptions3.gridinterplayer=1; vfoptions3.ngridinterp=5;
simoptions3=simoptions; simoptions3.gridinterplayer=1; simoptions3.ngridinterp=5;
[~,Policy_finalA_GI,Valt_finalA_GI,Policyalt_finalA_GI]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

vfoptions4=vfoptions; vfoptions4.divideandconquer=1; vfoptions4.gridinterplayer=1; vfoptions4.ngridinterp=5;
simoptions4=simoptions; simoptions4.gridinterplayer=1; simoptions4.ngridinterp=5;

% big-grid versions (used downstream for big-grid moments)
[~,Policy_finalA_big,Valt_finalA_big,Policyalt_finalA_big]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
[~,Policy_finalA_big_GI,Valt_finalA_big_GI,Policyalt_finalA_big_GI]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

%% With vs Without fastOLG (Naive)
transpathoptionsslow=transpathoptionsbaseline; transpathoptionsslow.fastOLG=0;
[VPath1slow,PolicyPath1slow,ValtPath1slow,PolicyaltPath1slow]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA, Policy_finalA, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions1);
[VPath1fast,PolicyPath1fast,ValtPath1fast,PolicyaltPath1fast]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA, Policy_finalA, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

fprintf('fastOLG (Naive), this should be zero: %.3e \n',max(abs(VPath1slow(:)-VPath1fast(:))))
fprintf('fastOLG (Naive), this should be zero: %.3e \n',max(abs(PolicyPath1slow(:)-PolicyPath1fast(:))))
fprintf('fastOLG (Naive, Valt), this should be zero: %.3e \n',max(abs(ValtPath1slow(:)-ValtPath1fast(:))))
fprintf('fastOLG (Naive, Policyalt), this should be zero: %.3e \n',max(abs(PolicyaltPath1slow(:)-PolicyaltPath1fast(:))))

clear VPath1fast VPath1slow PolicyPath1fast PolicyPath1slow ValtPath1fast ValtPath1slow PolicyaltPath1fast PolicyaltPath1slow

%% Divide-and-conquer (slowOLG, Naive)
[VPath1,PolicyPath1,ValtPath1,PolicyaltPath1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA, Policy_finalA, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions1);
[VPath2,PolicyPath2,ValtPath2,PolicyaltPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA, Policy_finalA, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions2);

fprintf('Divide-and-conquer (slowOLG, Naive), this should be zero: %.3e \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('Divide-and-conquer (slowOLG, Naive), this should be zero: %.3e \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))
fprintf('Divide-and-conquer (slowOLG, Naive, Valt), this should be zero: %.3e \n',max(abs(ValtPath1(:)-ValtPath2(:))))
fprintf('Divide-and-conquer (slowOLG, Naive, Policyalt), this should be zero: %.3e \n',max(abs(PolicyaltPath1(:)-PolicyaltPath2(:))))

clear VPath1 VPath2 PolicyPath1 PolicyPath2 ValtPath1 ValtPath2 PolicyaltPath1 PolicyaltPath2

%% Divide-and-conquer (fastOLG, Naive)
[VPath1,PolicyPath1,ValtPath1,PolicyaltPath1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA, Policy_finalA, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
[VPath2,PolicyPath2,ValtPath2,PolicyaltPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA, Policy_finalA, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);

fprintf('Divide-and-conquer (fastOLG, Naive), this should be zero: %.3e \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('Divide-and-conquer (fastOLG, Naive), this should be zero: %.3e \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))
fprintf('Divide-and-conquer (fastOLG, Naive, Valt), this should be zero: %.3e \n',max(abs(ValtPath1(:)-ValtPath2(:))))
fprintf('Divide-and-conquer (fastOLG, Naive, Policyalt), this should be zero: %.3e \n',max(abs(PolicyaltPath1(:)-PolicyaltPath2(:))))

clear VPath1 VPath2 PolicyPath1 PolicyPath2 ValtPath1 ValtPath2 PolicyaltPath1 PolicyaltPath2

%% Grid-interpolation +/- divide-and-conquer (slowOLG, Naive)
[VPath3,PolicyPath3,ValtPath3,PolicyaltPath3]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA_GI, Policy_finalA_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions3);
[VPath4,PolicyPath4,ValtPath4,PolicyaltPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA_GI, Policy_finalA_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions4);

fprintf('Divide-and-conquer (with GI, slowOLG, Naive), this should be zero: %.3e \n',max(abs(VPath3(:)-VPath4(:))))
fprintf('Divide-and-conquer (with GI, slowOLG, Naive), this should be zero: %.3e \n',max(abs(PolicyPath3(:)-PolicyPath4(:))))
fprintf('Divide-and-conquer (with GI, slowOLG, Naive, Valt), this should be zero: %.3e \n',max(abs(ValtPath3(:)-ValtPath4(:))))

%% Grid-interpolation +/- divide-and-conquer (fastOLG, Naive)
[VPath3,PolicyPath3,ValtPath3,PolicyaltPath3]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA_GI, Policy_finalA_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions3);
[VPath4,PolicyPath4,ValtPath4,PolicyaltPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA_GI, Policy_finalA_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);

fprintf('Divide-and-conquer (with GI, fastOLG, Naive), this should be zero: %.3e \n',max(abs(VPath3(:)-VPath4(:))))
fprintf('Divide-and-conquer (with GI, fastOLG, Naive), this should be zero: %.3e \n',max(abs(PolicyPath3(:)-PolicyPath4(:))))
fprintf('Divide-and-conquer (with GI, fastOLG, Naive, Valt), this should be zero: %.3e \n',max(abs(ValtPath3(:)-ValtPath4(:))))

clear VPath3 VPath4 PolicyPath3 PolicyPath4 ValtPath3 ValtPath4 PolicyaltPath3 PolicyaltPath4

%% Big a_grid: moments should be ~the same with/without grid interpolation (Naive)
[VPath2b,PolicyPath2b,~,~]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA_big, Policy_finalA_big, Params, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDist_initial_big=StationaryDist_FHorz_Case1(jequaloneDist_big,AgeWeightParamNames,PolicyPath2b(:,:,:,:,:,1),n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AgentDistPath2=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_big, jequaloneDist_big, PricePath, ParamPath, PolicyPath2b, AgeWeightParamNames,n_d,n_a_big,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions2);
AggVarsPath2=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath2, PolicyPath2b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, transpathoptionsbaseline, simoptions2);

[VPath4b,PolicyPath4b,~,~]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA_big_GI, Policy_finalA_big_GI, Params, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
AgentDistPath4=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_big, jequaloneDist_big, PricePath, ParamPath, PolicyPath4b, AgeWeightParamNames,n_d,n_a_big,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions4);
AggVarsPath4=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath4, PolicyPath4b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, transpathoptionsbaseline, simoptions4);

fprintf('With/without grid interp, should get much the same moments (for big a_grid, Naive) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %.3e \n',max(abs(AgentDistPath2(:)-AgentDistPath4(:))))
[AggVarsPath2.earnings.Mean; AggVarsPath4.earnings.Mean]
[AggVarsPath2.assets.Mean; AggVarsPath4.assets.Mean]

%% Naive figure
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:T,AggVarsPath2.earnings.Mean, 1:1:T,AggVarsPath4.earnings.Mean)
title('Earnings Mean (Naive)')
legend('1','2')
subplot(2,1,2); plot(1:1:T,AggVarsPath2.assets.Mean, 1:1:T,AggVarsPath4.assets.Mean)
title('Assets Mean (Naive)')
legend('1','2')

clear VPath2b VPath4b AggVarsPath2 AggVarsPath4 PolicyPath4b AgentDistPath2 AgentDistPath4

%% Constant-path identity (Naive): VPath should be repmat of steady-state V
[V1,Policy1,V1alt,Policy1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
PricePathConstant.r=Params.r*ones(1,T);
ParamPathConstant.sigma=Params.sigma*ones(1,T);
[VPath1,PolicyPath1,ValtPath1,PolicyaltPath1]=ValueFnOnTransPath_Case1_FHorz(PricePathConstant, ParamPathConstant, T, V1alt, Policy1, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

fprintf('Do nothing TPath (Naive), this should be zero, V: %.3e \n',max(abs(VPath1(:)-reshape(repmat(V1,1,1,1,1,T),[],1))))
fprintf('Do nothing TPath (Naive), this should be zero, Policy: %.3e \n',max(abs(PolicyPath1(:)-reshape(repmat(Policy1,1,1,1,1,1,T),[],1))))
fprintf('Do nothing TPath (Naive), this should be zero, Valt: %.3e \n',max(abs(ValtPath1(:)-reshape(repmat(V1alt,1,1,1,1,T),[],1))))
fprintf('Do nothing TPath (Naive), this should be zero, Policyalt: %.3e \n',max(abs(PolicyaltPath1(:)-reshape(repmat(Policy1alt,1,1,1,1,1,T),[],1))))

clear V1 Policy1 V1alt Policy1alt VPath1 PolicyPath1 ValtPath1 PolicyaltPath1

%% One-iter TransitionPath_Case1_FHorz shape test (Naive)
[~,PolicyPath2,~,~]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA, Policy_finalA, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyPath2(:,:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptions2);
clear PolicyPath2

transpathoptions=transpathoptionsbaseline;
transpathoptions.maxiter=1;
GeneralEqmEqns.dummy=@(earnings) 0;

transpathoptions.fastOLG=1;
PricePath2=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

transpathoptions.fastOLG=0;
PricePath2B=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

fprintf('One iter of TPath (Naive), with/without fastOLG, this should be zero: %.3e \n',max(abs(PricePath2.r-PricePath2B.r)))

% Big grid via DC + GI
transpathoptions.fastOLG=1;
PricePath3A=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA_big, AgentDist_initial_big, jequaloneDist_big, n_d, n_a_big, n_z, N_j, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions2, vfoptions2);
PricePath3B=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, Valt_finalA_big_GI, AgentDist_initial_big, jequaloneDist_big, n_d, n_a_big, n_z, N_j, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions4, vfoptions4);

fprintf('One iter of TPath (Naive), with/without GI, this should be close to zero: %.3e \n',max(abs(PricePath3A.r-PricePath3B.r)))


%% ============================================================
%% Sophisticated
%% ============================================================
vfoptions.quasi_hyperbolic='Sophisticated';

vfoptions1=vfoptions;
simoptions1=simoptions;
[~,Policy_finalS,Vunderbar_finalS]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

vfoptions2=vfoptions; vfoptions2.divideandconquer=1;
simoptions2=simoptions;

vfoptions3=vfoptions; vfoptions3.gridinterplayer=1; vfoptions3.ngridinterp=5;
simoptions3=simoptions; simoptions3.gridinterplayer=1; simoptions3.ngridinterp=5;
[~,Policy_finalS_GI,Vunderbar_finalS_GI]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

vfoptions4=vfoptions; vfoptions4.divideandconquer=1; vfoptions4.gridinterplayer=1; vfoptions4.ngridinterp=5;
simoptions4=simoptions; simoptions4.gridinterplayer=1; simoptions4.ngridinterp=5;

[~,Policy_finalS_big,Vunderbar_finalS_big]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
[~,Policy_finalS_big_GI,Vunderbar_finalS_big_GI]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

%% With vs Without fastOLG (Sophisticated)
[VPath1slow,PolicyPath1slow,ValtPath1slow]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS, Policy_finalS, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions1);
[VPath1fast,PolicyPath1fast,ValtPath1fast]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS, Policy_finalS, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

fprintf('fastOLG (Sophisticated), this should be zero: %.3e \n',max(abs(VPath1slow(:)-VPath1fast(:))))
fprintf('fastOLG (Sophisticated), this should be zero: %.3e \n',max(abs(PolicyPath1slow(:)-PolicyPath1fast(:))))
fprintf('fastOLG (Sophisticated, Vunderbar), this should be zero: %.3e \n',max(abs(ValtPath1slow(:)-ValtPath1fast(:))))

clear VPath1fast VPath1slow PolicyPath1fast PolicyPath1slow ValtPath1fast ValtPath1slow

%% Divide-and-conquer (slowOLG, Sophisticated)
[VPath1,PolicyPath1,ValtPath1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS, Policy_finalS, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions1);
[VPath2,PolicyPath2,ValtPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS, Policy_finalS, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions2);

fprintf('Divide-and-conquer (slowOLG, Sophisticated), this should be zero: %.3e \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('Divide-and-conquer (slowOLG, Sophisticated), this should be zero: %.3e \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))
fprintf('Divide-and-conquer (slowOLG, Sophisticated, Vunderbar), this should be zero: %.3e \n',max(abs(ValtPath1(:)-ValtPath2(:))))

clear VPath1 VPath2 PolicyPath1 PolicyPath2 ValtPath1 ValtPath2

%% Divide-and-conquer (fastOLG, Sophisticated)
[VPath1,PolicyPath1,ValtPath1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS, Policy_finalS, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
[VPath2,PolicyPath2,ValtPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS, Policy_finalS, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);

fprintf('Divide-and-conquer (fastOLG, Sophisticated), this should be zero: %.3e \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('Divide-and-conquer (fastOLG, Sophisticated), this should be zero: %.3e \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))
fprintf('Divide-and-conquer (fastOLG, Sophisticated, Vunderbar), this should be zero: %.3e \n',max(abs(ValtPath1(:)-ValtPath2(:))))

clear VPath1 VPath2 PolicyPath1 PolicyPath2 ValtPath1 ValtPath2

%% Grid-interpolation +/- divide-and-conquer (slowOLG, Sophisticated)
[VPath3,PolicyPath3,ValtPath3]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS_GI, Policy_finalS_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions3);
[VPath4,PolicyPath4,ValtPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS_GI, Policy_finalS_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions4);

fprintf('Divide-and-conquer (with GI, slowOLG, Sophisticated), this should be zero: %.3e \n',max(abs(VPath3(:)-VPath4(:))))
fprintf('Divide-and-conquer (with GI, slowOLG, Sophisticated), this should be zero: %.3e \n',max(abs(PolicyPath3(:)-PolicyPath4(:))))
fprintf('Divide-and-conquer (with GI, slowOLG, Sophisticated, Vunderbar), this should be zero: %.3e \n',max(abs(ValtPath3(:)-ValtPath4(:))))

%% Grid-interpolation +/- divide-and-conquer (fastOLG, Sophisticated)
[VPath3,PolicyPath3,ValtPath3]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS_GI, Policy_finalS_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions3);
[VPath4,PolicyPath4,ValtPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS_GI, Policy_finalS_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);

fprintf('Divide-and-conquer (with GI, fastOLG, Sophisticated), this should be zero: %.3e \n',max(abs(VPath3(:)-VPath4(:))))
fprintf('Divide-and-conquer (with GI, fastOLG, Sophisticated), this should be zero: %.3e \n',max(abs(PolicyPath3(:)-PolicyPath4(:))))
fprintf('Divide-and-conquer (with GI, fastOLG, Sophisticated, Vunderbar), this should be zero: %.3e \n',max(abs(ValtPath3(:)-ValtPath4(:))))

clear VPath3 VPath4 PolicyPath3 PolicyPath4 ValtPath3 ValtPath4

%% Big a_grid: moments should be ~the same with/without grid interpolation (Sophisticated)
[VPath2b,PolicyPath2b,~]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS_big, Policy_finalS_big, Params, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDist_initial_big=StationaryDist_FHorz_Case1(jequaloneDist_big,AgeWeightParamNames,PolicyPath2b(:,:,:,:,:,1),n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AgentDistPath2=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_big, jequaloneDist_big, PricePath, ParamPath, PolicyPath2b, AgeWeightParamNames,n_d,n_a_big,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions2);
AggVarsPath2=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath2, PolicyPath2b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, transpathoptionsbaseline, simoptions2);

[VPath4b,PolicyPath4b,~]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS_big_GI, Policy_finalS_big_GI, Params, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
AgentDistPath4=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_big, jequaloneDist_big, PricePath, ParamPath, PolicyPath4b, AgeWeightParamNames,n_d,n_a_big,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions4);
AggVarsPath4=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath4, PolicyPath4b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, transpathoptionsbaseline, simoptions4);

fprintf('With/without grid interp, should get much the same moments (for big a_grid, Sophisticated) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %.3e \n',max(abs(AgentDistPath2(:)-AgentDistPath4(:))))
[AggVarsPath2.earnings.Mean; AggVarsPath4.earnings.Mean]
[AggVarsPath2.assets.Mean; AggVarsPath4.assets.Mean]

clear VPath2b VPath4b AggVarsPath2 AggVarsPath4 PolicyPath4b AgentDistPath2 AgentDistPath4

%% Constant-path identity (Sophisticated)
[V1,Policy1,V1under]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
[VPath1,PolicyPath1,ValtPath1]=ValueFnOnTransPath_Case1_FHorz(PricePathConstant, ParamPathConstant, T, V1under, Policy1, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

fprintf('Do nothing TPath (Sophisticated), this should be zero, V: %.3e \n',max(abs(VPath1(:)-reshape(repmat(V1,1,1,1,1,T),[],1))))
fprintf('Do nothing TPath (Sophisticated), this should be zero, Policy: %.3e \n',max(abs(PolicyPath1(:)-reshape(repmat(Policy1,1,1,1,1,1,T),[],1))))
fprintf('Do nothing TPath (Sophisticated), this should be zero, Vunderbar: %.3e \n',max(abs(ValtPath1(:)-reshape(repmat(V1under,1,1,1,1,T),[],1))))

clear V1 Policy1 V1under VPath1 PolicyPath1 ValtPath1

%% One-iter TransitionPath_Case1_FHorz shape test (Sophisticated)
transpathoptions=transpathoptionsbaseline;
transpathoptions.maxiter=1;

transpathoptions.fastOLG=1;
PricePath2=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

transpathoptions.fastOLG=0;
PricePath2B=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

fprintf('One iter of TPath (Sophisticated), with/without fastOLG, this should be zero: %.3e \n',max(abs(PricePath2.r-PricePath2B.r)))

transpathoptions.fastOLG=1;
PricePath3A=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS_big, AgentDist_initial_big, jequaloneDist_big, n_d, n_a_big, n_z, N_j, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions2, vfoptions2);
PricePath3B=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, Vunderbar_finalS_big_GI, AgentDist_initial_big, jequaloneDist_big, n_d, n_a_big, n_z, N_j, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions4, vfoptions4);

fprintf('One iter of TPath (Sophisticated), with/without GI, this should be close to zero: %.3e \n',max(abs(PricePath3A.r-PricePath3B.r)))


%% ============================================================
%% Versus exponential discounting (beta0=1)
%% ============================================================
Params.beta0=1;

% Basic (no need to test divide-and-conquer again; established above)
vfoptions1.exoticpreferences='None';
[Va_final,Pa_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
[VPathA,PolicyPathA]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Va_final, Pa_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

vfoptions1.exoticpreferences='QuasiHyperbolic';
vfoptions1.quasi_hyperbolic='Naive';
[~,Pb_final,Vbalt_final,~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
[VPathB,PolicyPathB,ValtPathB,PolicyaltPathB]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vbalt_final, Pb_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

vfoptions1.quasi_hyperbolic='Sophisticated';
[~,Pc_final,Vcunder_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
[VPathC,PolicyPathC,ValtPathC]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vcunder_final, Pc_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

fprintf('QH-TPath with beta0=1 (Naive vs Exp), this should be zero: V    %.3e \n',max(abs(VPathA(:)-VPathB(:))))
fprintf('QH-TPath with beta0=1 (Naive vs Exp), this should be zero: P    %.3e \n',max(abs(PolicyPathA(:)-PolicyPathB(:))))
fprintf('QH-TPath with beta0=1 (Naive vs Exp, Valt), this should be zero: %.3e \n',max(abs(VPathA(:)-ValtPathB(:))))
fprintf('QH-TPath with beta0=1 (Soph vs Exp), this should be zero:  V    %.3e \n',max(abs(VPathA(:)-VPathC(:))))
fprintf('QH-TPath with beta0=1 (Soph vs Exp), this should be zero:  P    %.3e \n',max(abs(PolicyPathA(:)-PolicyPathC(:))))
fprintf('QH-TPath with beta0=1 (Soph vs Exp, Vunderbar), this should be zero: %.3e \n',max(abs(VPathA(:)-ValtPathC(:))))


%% Since QH preferences have no impact beyond Policy, there is no point testing the AgentDist/EvalFn pipelines beyond what's done above.

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
