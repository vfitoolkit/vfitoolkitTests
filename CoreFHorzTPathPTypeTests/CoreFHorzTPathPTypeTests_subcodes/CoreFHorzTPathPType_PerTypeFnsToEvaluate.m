function output=CoreFHorzTPathPType_PerTypeFnsToEvaluate(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i,transpathoptionsbaseline)
% Test 5: two PTypes that are actually identical. Solve once, then evaluate
% FnsToEvaluate two ways: (A) a single function returning the endogenous
% state a, broadcast to both types; (B) per-type functions (struct-of-structs,
% one sub-entry per type under the same field name, together covering all
% agents) each returning a. Both must give the same answer.
% Checked for EvalFnOnTransPath_AggVars_Case1_FHorz_PType (aggregate and
% per-type breakdowns) and for one iteration of TransitionPath_Case1_FHorz_PType.
% Only z is used (no d, no e, no semiz).
%
% AllStats and AgeConditionalStats on the transition path are NOT YET
% IMPLEMENTED for PType (see commented-out blocks at the end).

n_d=0;
d_grid=[];

% ReturnFn declares kappa_j (not kappa_j_pt), so both types are identical.
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

% (A) single FnsToEvaluate, broadcast to both types
FnsToEvaluate_single.assets=@(aprime,a,z) a;
% (B) per-type FnsToEvaluate: same field name, one sub-entry per type
FnsToEvaluate_pertype.assets.(Names_i{1})=@(aprime,a,z) a;
FnsToEvaluate_pertype.assets.(Names_i{2})=@(aprime,a,z) a;

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

% Trivial period-T V and Policy (as in CoreFHorzTPathTests)
V_final_PT=struct(); Policy_final_PT=struct();
for ii=1:N_i
    V_final_PT.(Names_i{ii})=zeros([n_a,n_z,N_j],'gpuArray');
    Policy_final_PT.(Names_i{ii})=ones([1,n_a,n_z,N_j],'gpuArray');
end

%% Solve the path machinery once (FnsToEvaluate does not enter these, so both cases share it)
[~,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
Policy_init_PT=struct();
for ii=1:N_i
    Policy_init_PT.(Names_i{ii})=PolicyPath_PT.(Names_i{ii})(:,:,:,:,1);
end
AgentDist_initial_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_init_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AgentDistPath_PT=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath, ParamPath, PolicyPath_PT, AgeWeightParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params, transpathoptionsbaseline, simoptions);

%% AggVars on the transition path: single vs per-type FnsToEvaluate
AggVarsPath_A=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate_single, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
AggVarsPath_B=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate_pertype, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);

fprintf('per-type FnsToEvaluate, AggVarsPath assets.Mean, this should be zero: %.3e \n',max(abs(AggVarsPath_A.assets.Mean(:)-AggVarsPath_B.assets.Mean(:))))
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('per-type FnsToEvaluate, AggVarsPath assets breakdown (type %s), this should be zero: %.3e \n',nm,max(abs(AggVarsPath_A.assets.(nm).Mean(:)-AggVarsPath_B.assets.(nm).Mean(:))))
end

%% One iteration of the GE transition path: single vs per-type FnsToEvaluate
transpathoptions=transpathoptionsbaseline;
transpathoptions.maxiter=1;

GeneralEqmEqns.dummy=@(assets) 0;

PricePath_A=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_single, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
PricePath_B=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_pertype, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);

fprintf('per-type FnsToEvaluate, one iter of GE TPath, this should be zero: %.3e \n',max(abs(PricePath_A.r(:)-PricePath_B.r(:))))

%% AllStats on the transition path: single vs per-type FnsToEvaluate
% NOT YET IMPLEMENTED: there is no EvalFnOnTransPath_AllStats_Case1_FHorz_PType
% (nor a non-PType FHorz EvalFnOnTransPath_AllStats; only InfHorz has one).
% Uncomment once the command exists.
% AllStatsPath_A=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate_single, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
% AllStatsPath_B=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate_pertype, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
% fprintf('per-type FnsToEvaluate, AllStatsPath assets.Mean, this should be zero: %.3e \n',max(abs(AllStatsPath_A.assets.Mean(:)-AllStatsPath_B.assets.Mean(:))))
% fprintf('per-type FnsToEvaluate, AllStatsPath assets.Gini, this should be zero: %.3e \n',max(abs(AllStatsPath_A.assets.Gini(:)-AllStatsPath_B.assets.Gini(:))))

%% AgeConditionalStats on the transition path: single vs per-type FnsToEvaluate
% NOT YET IMPLEMENTED: LifeCycleProfiles_TransPath_FHorz_Case1 exists but has
% no PType version. Uncomment once LifeCycleProfiles_TransPath_FHorz_Case1_PType
% exists.
% AgeCondStatsPath_A=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate_single, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, simoptions);
% AgeCondStatsPath_B=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate_pertype, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, simoptions);
% fprintf('per-type FnsToEvaluate, AgeCondStatsPath assets.Mean, this should be zero: %.3e \n',max(abs(AgeCondStatsPath_A.assets.Mean(:)-AgeCondStatsPath_B.assets.Mean(:))))

output=struct();

end
