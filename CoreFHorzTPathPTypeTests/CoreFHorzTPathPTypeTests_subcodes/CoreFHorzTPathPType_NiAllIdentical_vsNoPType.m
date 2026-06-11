function output=CoreFHorzTPathPType_NiAllIdentical_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,transpathoptionsbaseline)
% Test 2: PType TPath with N_i=2 identical types should equal a single
% no-PType TPath solve: each type's VPath/PolicyPath/AgentDistPath matches
% the no-PType paths, and since ptypeweights sum to one the aggregated
% AggVarsPath does too. Only z is used (no d, no e, no semiz).

n_d=0;
d_grid=[];

% ReturnFn declares kappa_j (not kappa_j_pt), so both types are identical.
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

% Trivial period-T V and Policy (as in CoreFHorzTPathTests)
V_final=zeros([n_a,n_z,N_j],'gpuArray');
Policy_final=ones([1,n_a,n_z,N_j],'gpuArray');
% Names-keyed structs for the PType commands
names_PT={'ptype001','ptype002'}; % auto-names used when N_i is numeric
V_final_PT=struct(); Policy_final_PT=struct();
for ii=1:N_i
    V_final_PT.(names_PT{ii})=V_final;
    Policy_final_PT.(names_PT{ii})=Policy_final;
end

%% PType solve
[VPath_PT,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
Policy_init_PT=struct();
for ii=1:N_i
    Policy_init_PT.(names_PT{ii})=PolicyPath_PT.(names_PT{ii})(:,:,:,:,1);
end
AgentDist_initial_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_init_PT,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
AgentDistPath_PT=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath, ParamPath, PolicyPath_PT, AgeWeightParamNames,n_d,n_a,n_z,N_j,N_i,pi_z, T,Params, transpathoptionsbaseline, simoptions);
AggVarsPath_PT=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);

%% Single no-PType solve
[VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyPath(:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
AgentDistPath=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath, PolicyPath, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions);
AggVarsPath=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath, PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);

%% Compare per-type slices to the single no-PType solve
for ii=1:N_i
    nm=names_PT{ii};
    fprintf('All identical vs no-PType, VPath        (type %d), this should be zero: %2.8f \n',ii,max(abs(VPath_PT.(nm)(:)-VPath(:))))
    fprintf('All identical vs no-PType, PolicyPath   (type %d), this should be zero: %2.8f \n',ii,max(abs(PolicyPath_PT.(nm)(:)-PolicyPath(:))))
    fprintf('All identical vs no-PType, AgentDistPath(type %d), this should be zero: %2.8f \n',ii,max(abs(AgentDistPath_PT.(nm)(:)-AgentDistPath(:))))
end

%% Aggregated AggVarsPath equals the no-PType AggVarsPath (ptypeweights sum to one)
fprintf('All identical vs no-PType, AggVarsPath assets.Mean,   this should be zero: %2.8f \n',max(abs(AggVarsPath_PT.assets.Mean(:)  -AggVarsPath.assets.Mean(:))))
fprintf('All identical vs no-PType, AggVarsPath earnings.Mean, this should be zero: %2.8f \n',max(abs(AggVarsPath_PT.earnings.Mean(:)-AggVarsPath.earnings.Mean(:))))

output=struct();

end
