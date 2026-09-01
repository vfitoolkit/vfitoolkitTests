function output=CoreFHorzPType_NiTwoTypes_vsTwoSolves(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test 4: PType solve with N_i=2 genuinely different types should equal solving
% the same model twice without PType, once per type, and stacking. Per-type
% V/Policy/Dist must match the corresponding single-type solves; aggregated
% AllStats.Mean must equal the ptypeweight-weighted combination of the two
% per-type means.

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

%% PType solve
[V_PT,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
V_PT_vfp=ValueFnFromPolicy_FHorz_PType(Policy_PT,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
AllStats_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);

%% Two separate single-type solves
names_PT=fieldnames(V_PT);
ptw=Params.(PTypeDistParamNames{:});
V_solo=cell(N_i,1); V_solo_vfp=cell(N_i,1); Policy_solo=cell(N_i,1); Dist_solo=cell(N_i,1); AllStats_solo=cell(N_i,1);
for ii=1:N_i
    Params_ii=Params;
    Params_ii.kappa_j=Params.kappa_j_pt(ii,:);
    [V_solo{ii},Policy_solo{ii}]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,[],vfoptions);
    V_solo_vfp{ii}=ValueFnFromPolicy_FHorz(Policy_solo{ii},n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,vfoptions);
    Dist_solo{ii}=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo{ii},n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    AllStats_solo{ii}=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist_solo{ii},Policy_solo{ii},FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
end

%% Compare per-type slices
for ii=1:N_i
    nA=names_PT{ii};
    fprintf('Two PTypes vs two solves, V    (type %d), this should be zero: %.3e \n',ii,max(abs(V_PT.(nA)(:)-V_solo{ii}(:))))
    fprintf('Two PTypes vs two solves, Pol  (type %d), this should be zero: %.3e \n',ii,max(abs(Policy_PT.(nA)(:)-Policy_solo{ii}(:))))
    fprintf('Two PTypes vs two solves, Dist (type %d), this should be zero: %.3e \n',ii,max(abs(StationaryDist_PT.(nA)(:)-Dist_solo{ii}(:))))
    fprintf('Two PTypes vs two solves, VFP  vs V (PType,   type %d), this should be zero: %.3e \n',ii,max(abs(V_PT_vfp.(nA)(:)-V_PT.(nA)(:))))
    fprintf('Two PTypes vs two solves, VFP  vs V (noPType, type %d), this should be zero: %.3e \n',ii,max(abs(V_solo_vfp{ii}(:)-V_solo{ii}(:))))
end

%% Aggregated AllStats.Mean equals ptw-weighted combination of solo means
agg_assets_mean   = ptw(1)*AllStats_solo{1}.assets.Mean   + ptw(2)*AllStats_solo{2}.assets.Mean;
agg_earnings_mean = ptw(1)*AllStats_solo{1}.earnings.Mean + ptw(2)*AllStats_solo{2}.earnings.Mean;
fprintf('Two PTypes vs two solves, AllStats assets.Mean,   this should be zero: %.3e \n',abs(AllStats_PT.assets.Mean   -agg_assets_mean))
fprintf('Two PTypes vs two solves, AllStats earnings.Mean, this should be zero: %.3e \n',abs(AllStats_PT.earnings.Mean -agg_earnings_mean))

output=struct();

end
