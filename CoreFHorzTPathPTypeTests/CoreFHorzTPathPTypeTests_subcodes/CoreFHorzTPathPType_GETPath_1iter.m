function output=CoreFHorzTPathPType_GETPath_1iter(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i)
% Test 4: Run the GE transition path (TransitionPath_Case1_FHorz_PType) with
% transpathoptions.maxiter=1, so it ends after one iteration. With N_i=2
% identical types (and ptypeweights summing to one) the aggregates equal the
% no-PType aggregates, so the one-iteration price path update must equal the
% no-PType TransitionPath_Case1_FHorz one. Checked with fastOLG=1 and
% fastOLG=0. Only z is used (no d, no e, no semiz).

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
names_PT={'ptype001','ptype002'}; % auto-names used when N_i is numeric
V_final_PT=struct(); Policy_final_PT=struct();
for ii=1:N_i
    V_final_PT.(names_PT{ii})=V_final;
    Policy_final_PT.(names_PT{ii})=Policy_final;
end

%% Create the initial agent distributions (from the first period of the policy path)
transpathoptions.fastOLG=1;
[~,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyPath(:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
clear PolicyPath

[~,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
Policy_init_PT=struct();
for ii=1:N_i
    Policy_init_PT.(names_PT{ii})=PolicyPath_PT.(names_PT{ii})(:,:,:,:,1);
end
AgentDist_initial_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_init_PT,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
clear PolicyPath_PT

%% One iteration of the GE transition path, with fastOLG=1
transpathoptions.maxiter=1;

GeneralEqmEqns.dummy=@(earnings) 0;

transpathoptions.fastOLG=1;
PricePath_PT1=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, N_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
PricePath_NoPT1=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

fprintf('One iter of GE TPath (fastOLG=1), PType vs no-PType, this should be zero: %2.8f \n',max(abs(PricePath_PT1.r(:)-PricePath_NoPT1.r(:))))

%% And again with fastOLG=0
transpathoptions.fastOLG=0;
PricePath_PT0=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, N_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
PricePath_NoPT0=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

fprintf('One iter of GE TPath (fastOLG=0), PType vs no-PType, this should be zero: %2.8f \n',max(abs(PricePath_PT0.r(:)-PricePath_NoPT0.r(:))))

fprintf('One iter of GE TPath, PType with/without fastOLG, this should be zero: %2.8f \n',max(abs(PricePath_PT1.r(:)-PricePath_PT0.r(:))))

output=struct();

end
