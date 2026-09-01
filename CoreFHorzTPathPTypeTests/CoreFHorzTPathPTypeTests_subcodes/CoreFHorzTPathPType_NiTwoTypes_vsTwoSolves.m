function output=CoreFHorzTPathPType_NiTwoTypes_vsTwoSolves(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,transpathoptionsbaseline,figure_c)
% Test 3: PType TPath with N_i=2 genuinely different types should equal
% solving the same transition twice without PType, once per type, and
% stacking. Per-type VPath/PolicyPath/AgentDistPath must match the
% corresponding single-type solves; the aggregated AggVarsPath.Mean must
% equal the ptypeweight-weighted combination of the two per-type means.
% Only z is used (no d, no e, no semiz).

N_i=2;

n_d=0;
d_grid=[];

% Per-type kappa_j (already set up as Params.kappa_j_pt = [2 x N_j])
% PType solve uses kappa_j_pt; single-type solves overwrite Params.kappa_j per type.
ReturnFn_PT=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);
ReturnFn_NoPT=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate_PT.assets=@(aprime,a,z) a;
FnsToEvaluate_PT.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate_NoPT.assets=@(aprime,a,z) a;
FnsToEvaluate_NoPT.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;

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

%% PType solve
[VPath_PT,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn_PT, transpathoptionsbaseline, vfoptions);
Policy_init_PT=struct();
for ii=1:N_i
    Policy_init_PT.(names_PT{ii})=PolicyPath_PT.(names_PT{ii})(:,:,:,:,1);
end
AgentDist_initial_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_init_PT,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
AgentDistPath_PT=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath, ParamPath, PolicyPath_PT, AgeWeightParamNames,n_d,n_a,n_z,N_j,N_i,pi_z, T,Params, transpathoptionsbaseline, simoptions);
AggVarsPath_PT=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate_PT, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);

%% Two separate single-type solves
ptw=Params.(PTypeDistParamNames{:});
VPath_solo=cell(N_i,1); PolicyPath_solo=cell(N_i,1); AgentDistPath_solo=cell(N_i,1); AggVarsPath_solo=cell(N_i,1);
for ii=1:N_i
    Params_ii=Params;
    Params_ii.kappa_j=Params.kappa_j_pt(ii,:);
    [VPath_solo{ii},PolicyPath_solo{ii}]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params_ii, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn_NoPT, transpathoptionsbaseline, vfoptions);
    AgentDist_initial_ii=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyPath_solo{ii}(:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    AgentDistPath_solo{ii}=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_ii, jequaloneDist, PricePath, ParamPath, PolicyPath_solo{ii}, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params_ii, transpathoptionsbaseline, simoptions);
    AggVarsPath_solo{ii}=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate_NoPT, AgentDistPath_solo{ii}, PolicyPath_solo{ii}, PricePath, ParamPath, Params_ii, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
end

%% Compare per-type slices
for ii=1:N_i
    nm=names_PT{ii};
    fprintf('Two PTypes vs two solves, VPath        (type %d), this should be zero: %.3e \n',ii,max(abs(VPath_PT.(nm)(:)-VPath_solo{ii}(:))))
    fprintf('Two PTypes vs two solves, PolicyPath   (type %d), this should be zero: %.3e \n',ii,max(abs(PolicyPath_PT.(nm)(:)-PolicyPath_solo{ii}(:))))
    fprintf('Two PTypes vs two solves, AgentDistPath(type %d), this should be zero: %.3e \n',ii,max(abs(AgentDistPath_PT.(nm)(:)-AgentDistPath_solo{ii}(:))))
end

%% Aggregated AggVarsPath.Mean equals ptw-weighted combination of solo mean paths
agg_assets_mean  =ptw(1)*AggVarsPath_solo{1}.assets.Mean(:)  +ptw(2)*AggVarsPath_solo{2}.assets.Mean(:);
agg_earnings_mean=ptw(1)*AggVarsPath_solo{1}.earnings.Mean(:)+ptw(2)*AggVarsPath_solo{2}.earnings.Mean(:);
fprintf('Two PTypes vs two solves, AggVarsPath assets.Mean,   this should be zero: %.3e \n',max(abs(AggVarsPath_PT.assets.Mean(:)  -agg_assets_mean)))
fprintf('Two PTypes vs two solves, AggVarsPath earnings.Mean, this should be zero: %.3e \n',max(abs(AggVarsPath_PT.earnings.Mean(:)-agg_earnings_mean)))

%% Do some graphs of the AggVars path to see them
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:T,AggVarsPath_PT.earnings.Mean, 1:1:T,agg_earnings_mean)
title('Earnings Mean (PType vs weighted two solves)')
legend('PType','weighted solo')
subplot(2,1,2); plot(1:1:T,AggVarsPath_PT.assets.Mean, 1:1:T,agg_assets_mean)
title('Assets Mean (PType vs weighted two solves)')
legend('PType','weighted solo')

output=struct();

end
