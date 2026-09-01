function output=CoreFHorzTPathExpAssetz_CrossTest2_d1(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)

% Cross-test 2 for TPath with experienceassetz: 'fake' experienceassetz whose aprimeFn ignores z, vs plain experienceasset.
% Both models have z as an ordinary Markov shock (used in the ReturnFn) -- only the asset-type machinery differs.
% Should give same VPath, PolicyPath, AgentDistPath, AggVarsPath.

ReturnFn=@(d1,d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_TPathExpAssetz_d1_z_noe_nosemiz(d1,d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

FnsToEvaluate.assets=@(d1,d2,a1prime,a1,a2,z) a1+a2;
FnsToEvaluate.earnings=@(d1,d2,a1prime,a1,a2,z,w,kappa_j) w*kappa_j*d1*d2*a2*z;

% 'fake' experienceassetz aprimeFn: takes z but ignores it
aprimeFn_fakez=@(d2,a2,z,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;
% Plain experienceasset aprimeFn
aprimeFn_plain=@(d2,a2,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;

% Model A: 'fake' experienceassetz that ignores z
vfoptionsA=struct();
vfoptionsA.experienceassetz=1;
vfoptionsA.aprimeFn=aprimeFn_fakez;
simoptionsA=struct();
simoptionsA.experienceassetz=1;
simoptionsA.aprimeFn=aprimeFn_fakez;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;
simoptionsA.z_grid=z_grid;

% Model B: plain experienceasset (with same z shock present, but aprimeFn doesn't depend on z)
vfoptionsB=struct();
vfoptionsB.experienceasset=1;
vfoptionsB.aprimeFn=aprimeFn_plain;
simoptionsB=struct();
simoptionsB.experienceasset=1;
simoptionsB.aprimeFn=aprimeFn_plain;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;

% zero assets, mid points for any shocks
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,1,ceil(n_z/2))=1;

% Need period T for V and Policy
V_final=zeros([n_a,n_z,N_j],'gpuArray');
Policy_final=ones([3,n_a,n_z,N_j],'gpuArray');
Policy_final_GI=ones([5,n_a,n_z,N_j],'gpuArray');

%% Plain (no GI), slowOLG and fastOLG
transpathoptionsslow.fastOLG=0;
[VPathAslow,PolicyPathAslow]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptionsA);
[VPathBslow,PolicyPathBslow]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptionsB);

fprintf('CrossTest2 (fake-z expassetz vs plain expasset, slowOLG), this should be zero: V %.3e, Policy %.3e \n', max(abs(VPathAslow(:)-VPathBslow(:))), max(abs(PolicyPathAslow(:)-PolicyPathBslow(:))))

clear VPathAslow VPathBslow PolicyPathAslow PolicyPathBslow

[VPathA,PolicyPathA]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsA);
[VPathB,PolicyPathB]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsB);

fprintf('CrossTest2 (fake-z expassetz vs plain expasset, fastOLG), this should be zero: V %.3e, Policy %.3e \n', max(abs(VPathA(:)-VPathB(:))), max(abs(PolicyPathA(:)-PolicyPathB(:))))

%% AgentDist and AggVars on the path (using the fastOLG PolicyPaths)
AgentDist_initialA=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyPathA(:,:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);
AgentDistPathA=AgentDistOnTransPath_Case1_FHorz(AgentDist_initialA, jequaloneDist, PricePath, ParamPath, PolicyPathA, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptionsA);
AggVarsPathA=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPathA, PolicyPathA, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptionsA);

AgentDist_initialB=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyPathB(:,:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptionsB);
AgentDistPathB=AgentDistOnTransPath_Case1_FHorz(AgentDist_initialB, jequaloneDist, PricePath, ParamPath, PolicyPathB, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptionsB);
AggVarsPathB=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPathB, PolicyPathB, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptionsB);

fprintf('CrossTest2 (fake-z expassetz vs plain expasset), this should be zero: AgentDistPath %.3e, AggVars earnings %.3e, AggVars assets %.3e \n', max(abs(AgentDistPathA(:)-AgentDistPathB(:))), max(abs(AggVarsPathA.earnings.Mean(:)-AggVarsPathB.earnings.Mean(:))), max(abs(AggVarsPathA.assets.Mean(:)-AggVarsPathB.assets.Mean(:))))

clear VPathA VPathB PolicyPathA PolicyPathB AgentDistPathA AgentDistPathB AggVarsPathA AggVarsPathB

%% With grid interpolation layer (fastOLG)
vfoptionsA.gridinterplayer=1;
vfoptionsA.ngridinterp=5;
vfoptionsB.gridinterplayer=1;
vfoptionsB.ngridinterp=5;
[VPathA,PolicyPathA]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsA);
[VPathB,PolicyPathB]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsB);

fprintf('CrossTest2 (fake-z expassetz vs plain expasset, with GI, fastOLG), this should be zero: V %.3e, Policy %.3e \n', max(abs(VPathA(:)-VPathB(:))), max(abs(PolicyPathA(:)-PolicyPathB(:))))

output=struct();

end
