function output=CoreFHorzTPathPType_NivsNames(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i,transpathoptionsbaseline)
% Test 1: TPath with N_i (numeric) vs Names_i (cell) should give identical
% output. Both types are genuinely different (different kappa_j via kappa_j_pt).
% Only z is used (no d, no e, no semiz).
% The N_i solve passes jequaloneDist broadcast (one [n_a,n_z] array for all
% types), the Names_i solve passes it as a Names_i-keyed structure, so this
% also tests both accepted jequaloneDist forms.

n_d=0;
d_grid=[];

% Per-PType parameter: kappa_j_pt is [N_i x N_j], so toolkit picks row ii for type ii.
% ReturnFn declares kappa_j_pt as a parameter (rather than kappa_j), so type-specific
% earnings are read from Params.kappa_j_pt.
ReturnFn=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

% Trivial period-T V and Policy (as in CoreFHorzTPathTests), but as Names-keyed
% structs since the PType TPath commands require V_final.(typename).
names_A={'ptype001','ptype002'}; % auto-names used when N_i is numeric
V_final_A=struct(); Policy_final_A=struct();
V_final_B=struct(); Policy_final_B=struct();
jequaloneDist_B=struct();
for ii=1:N_i
    V_final_A.(names_A{ii})=zeros([n_a,n_z,N_j],'gpuArray');
    Policy_final_A.(names_A{ii})=ones([1,n_a,n_z,N_j],'gpuArray');
    V_final_B.(Names_i{ii})=zeros([n_a,n_z,N_j],'gpuArray');
    Policy_final_B.(Names_i{ii})=ones([1,n_a,n_z,N_j],'gpuArray');
    jequaloneDist_B.(Names_i{ii})=jequaloneDist;
end

%% Solve with N_i (numeric) — auto-names ptype001, ptype002; jequaloneDist broadcast
[VPath_A,PolicyPath_A]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_A, Policy_final_A, Params, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
Policy_init_A=struct();
for ii=1:N_i
    Policy_init_A.(names_A{ii})=PolicyPath_A.(names_A{ii})(:,:,:,:,1);
end
AgentDist_initial_A=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_init_A,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
AgentDistPath_A=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_A, jequaloneDist, PricePath, ParamPath, PolicyPath_A, AgeWeightParamNames,n_d,n_a,n_z,N_j,N_i,pi_z, T,Params, transpathoptionsbaseline, simoptions);
AggVarsPath_A=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_A, PolicyPath_A, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);

%% Solve with Names_i (cell); jequaloneDist as Names_i-keyed struct
[VPath_B,PolicyPath_B]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_B, Policy_final_B, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
Policy_init_B=struct();
for ii=1:N_i
    Policy_init_B.(Names_i{ii})=PolicyPath_B.(Names_i{ii})(:,:,:,:,1);
end
AgentDist_initial_B=StationaryDist_Case1_FHorz_PType(jequaloneDist_B,AgeWeightParamNames,PTypeDistParamNames,Policy_init_B,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AgentDistPath_B=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_B, jequaloneDist_B, PricePath, ParamPath, PolicyPath_B, AgeWeightParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params, transpathoptionsbaseline, simoptions);
AggVarsPath_B=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_B, PolicyPath_B, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);

%% Compare per-type (struct fields are ordered by ii, regardless of name)
for ii=1:N_i
    nA=names_A{ii};
    nB=Names_i{ii};
    fprintf('N_i vs Names_i, VPath        (type %d), this should be zero: %2.8f \n',ii,max(abs(VPath_A.(nA)(:)-VPath_B.(nB)(:))))
    fprintf('N_i vs Names_i, PolicyPath   (type %d), this should be zero: %2.8f \n',ii,max(abs(PolicyPath_A.(nA)(:)-PolicyPath_B.(nB)(:))))
    fprintf('N_i vs Names_i, AgentDistPath(type %d), this should be zero: %2.8f \n',ii,max(abs(AgentDistPath_A.(nA)(:)-AgentDistPath_B.(nB)(:))))
    fprintf('N_i vs Names_i, AggVarsPath assets   breakdown (type %d), this should be zero: %2.8f \n',ii,max(abs(AggVarsPath_A.assets.(nA).Mean(:)-AggVarsPath_B.assets.(nB).Mean(:))))
    fprintf('N_i vs Names_i, AggVarsPath earnings breakdown (type %d), this should be zero: %2.8f \n',ii,max(abs(AggVarsPath_A.earnings.(nA).Mean(:)-AggVarsPath_B.earnings.(nB).Mean(:))))
end
fprintf('N_i vs Names_i, ptweights, this should be zero: %2.8f \n',max(abs(AgentDistPath_A.ptweights(:)-AgentDistPath_B.ptweights(:))))
fprintf('N_i vs Names_i, AggVarsPath assets.Mean,   this should be zero: %2.8f \n',max(abs(AggVarsPath_A.assets.Mean(:)  -AggVarsPath_B.assets.Mean(:))))
fprintf('N_i vs Names_i, AggVarsPath earnings.Mean, this should be zero: %2.8f \n',max(abs(AggVarsPath_A.earnings.Mean(:)-AggVarsPath_B.earnings.Mean(:))))

output=struct();

end
