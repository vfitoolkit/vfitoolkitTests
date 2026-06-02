function output=CoreFHorzPType_NivsNames(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i)
% Test 1: solving with N_i (numeric) vs Names_i (cell) should give identical
% output. Both types are genuinely different (different kappa_j via kappa_j_pt).
% Only z is used (no d, no e, no semiz).

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

%% Solve with N_i (numeric) — auto-names ptype001, ptype002
[V_A,Policy_A]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
V_A_vfp=ValueFnFromPolicy_FHorz_PType(Policy_A,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_A=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_A,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);

%% Solve with Names_i (cell)
[V_B,Policy_B]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
V_B_vfp=ValueFnFromPolicy_FHorz_PType(Policy_B,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_B=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_B,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_B,Policy_B,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

%% Compare per-type (struct fields are ordered by ii, regardless of name)
names_A=fieldnames(V_A);
for ii=1:N_i
    nA=names_A{ii};
    nB=Names_i{ii};
    fprintf('N_i vs Names_i, V    (type %d), this should be zero: %2.8f \n',ii,max(abs(V_A.(nA)(:)-V_B.(nB)(:))))
    fprintf('N_i vs Names_i, Pol  (type %d), this should be zero: %2.8f \n',ii,max(abs(Policy_A.(nA)(:)-Policy_B.(nB)(:))))
    fprintf('N_i vs Names_i, Dist (type %d), this should be zero: %2.8f \n',ii,max(abs(StationaryDist_A.(nA)(:)-StationaryDist_B.(nB)(:))))
    fprintf('N_i vs Names_i, VFP  vs V (N_i,    type %d), this should be zero: %2.8f \n',ii,max(abs(V_A_vfp.(nA)(:)-V_A.(nA)(:))))
    fprintf('N_i vs Names_i, VFP  vs V (Names_i,type %d), this should be zero: %2.8f \n',ii,max(abs(V_B_vfp.(nB)(:)-V_B.(nB)(:))))
end
fprintf('N_i vs Names_i, ptweights, this should be zero: %2.8f \n',max(abs(StationaryDist_A.ptweights-StationaryDist_B.ptweights)))
fprintf('N_i vs Names_i, AllStats assets.Mean,   this should be zero: %2.8f \n',abs(AllStats_A.assets.Mean   -AllStats_B.assets.Mean))
fprintf('N_i vs Names_i, AllStats earnings.Gini, this should be zero: %2.8f \n',abs(AllStats_A.earnings.Gini -AllStats_B.earnings.Gini))

output=struct();

end
