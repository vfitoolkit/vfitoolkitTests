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
    fprintf('N_i vs Names_i, V    (type %d), this should be zero: %.3e \n',ii,max(abs(V_A.(nA)(:)-V_B.(nB)(:))))
    fprintf('N_i vs Names_i, Pol  (type %d), this should be zero: %.3e \n',ii,max(abs(Policy_A.(nA)(:)-Policy_B.(nB)(:))))
    fprintf('N_i vs Names_i, Dist (type %d), this should be zero: %.3e \n',ii,max(abs(StationaryDist_A.(nA)(:)-StationaryDist_B.(nB)(:))))
    fprintf('N_i vs Names_i, VFP  vs V (N_i,    type %d), this should be zero: %.3e \n',ii,max(abs(V_A_vfp.(nA)(:)-V_A.(nA)(:))))
    fprintf('N_i vs Names_i, VFP  vs V (Names_i,type %d), this should be zero: %.3e \n',ii,max(abs(V_B_vfp.(nB)(:)-V_B.(nB)(:))))
end
fprintf('N_i vs Names_i, ptweights, this should be zero: %.3e \n',max(abs(StationaryDist_A.ptweights-StationaryDist_B.ptweights)))
fprintf('N_i vs Names_i, AllStats assets.Mean,   this should be zero: %.3e \n',abs(AllStats_A.assets.Mean   -AllStats_B.assets.Mean))
fprintf('N_i vs Names_i, AllStats earnings.Gini, this should be zero: %.3e \n',abs(AllStats_A.earnings.Gini -AllStats_B.earnings.Gini))

%% Extension: per-type DiscountFactor (vector form) + per-type ReturnFn param (struct form)
% Above tested kappa_j_pt as an [N_i x N_j] matrix. This section adds two more
% per-type parameter paths not hit by the kappa_j_pt test:
%   - beta as a length-N_i vector (DiscountFactorParamNames)
%   - sigma as a struct keyed by type names (ReturnFn parameter)
% N_i (auto-names ptype001/2) and Names_i must give the same answer when the
% sigma struct keys are renamed to match each naming scheme.

Params.beta_pt=[0.95; 0.92];
DiscountFactorParamNames2={'beta_pt'};

%% Solve with N_i — sigma struct keyed by auto-names
Params.sigma=struct();
Params.sigma.ptype001=2;
Params.sigma.ptype002=3;
[V_A2,Policy_A2]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames2,vfoptions);
V_A2_vfp=ValueFnFromPolicy_FHorz_PType(Policy_A2,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames2,vfoptions);
StationaryDist_A2=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_A2,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
AllStats_A2=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_A2,Policy_A2,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);

%% Solve with Names_i — sigma struct keyed by Names_i
Params.sigma=struct();
Params.sigma.(Names_i{1})=2;
Params.sigma.(Names_i{2})=3;
[V_B2,Policy_B2]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames2,vfoptions);
V_B2_vfp=ValueFnFromPolicy_FHorz_PType(Policy_B2,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames2,vfoptions);
StationaryDist_B2=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_B2,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AllStats_B2=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_B2,Policy_B2,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

%% Compare per-type
for ii=1:N_i
    nA=names_A{ii};
    nB=Names_i{ii};
    fprintf('N_i vs Names_i (vec beta, struct sigma), V    (type %d), this should be zero: %.3e \n',ii,max(abs(V_A2.(nA)(:)-V_B2.(nB)(:))))
    fprintf('N_i vs Names_i (vec beta, struct sigma), Pol  (type %d), this should be zero: %.3e \n',ii,max(abs(Policy_A2.(nA)(:)-Policy_B2.(nB)(:))))
    fprintf('N_i vs Names_i (vec beta, struct sigma), Dist (type %d), this should be zero: %.3e \n',ii,max(abs(StationaryDist_A2.(nA)(:)-StationaryDist_B2.(nB)(:))))
    fprintf('N_i vs Names_i (vec beta, struct sigma), VFP vs V (N_i,    type %d), this should be zero: %.3e \n',ii,max(abs(V_A2_vfp.(nA)(:)-V_A2.(nA)(:))))
    fprintf('N_i vs Names_i (vec beta, struct sigma), VFP vs V (Names_i,type %d), this should be zero: %.3e \n',ii,max(abs(V_B2_vfp.(nB)(:)-V_B2.(nB)(:))))
end
fprintf('N_i vs Names_i (vec beta, struct sigma), AllStats assets.Mean,   this should be zero: %.3e \n',abs(AllStats_A2.assets.Mean   -AllStats_B2.assets.Mean))
fprintf('N_i vs Names_i (vec beta, struct sigma), AllStats earnings.Gini, this should be zero: %.3e \n',abs(AllStats_A2.earnings.Gini -AllStats_B2.earnings.Gini))

%% Per-type FnsToEvaluate (struct-of-structs, different fieldnames per type) ≡ single FnsToEvaluate
% Split each of the single 'assets'/'earnings' functions into two per-type fields
% (different field names per type) where each sub-function evaluates the same
% expression. AllStats normalizes Mean by the mass of *relevant* types only, so
% grouped means differ across the two setups, but the per-type breakdown
% (AllStats.<fn>.<typename>) for the relevant type must match the single-fn case.

FnsToEvaluate_PT=struct();
FnsToEvaluate_PT.assets_low.(Names_i{1})    = @(aprime,a,z) a;
FnsToEvaluate_PT.assets_high.(Names_i{2})   = @(aprime,a,z) a;
FnsToEvaluate_PT.earnings_low.(Names_i{1})  = @(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate_PT.earnings_high.(Names_i{2}) = @(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;

AllStats_B2_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_B2,Policy_B2,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

fprintf('per-type FnsToEvaluate, assets   breakdown (type %s, single==per-type), this should be zero: %.3e \n',Names_i{1},abs(AllStats_B2.assets.(Names_i{1}).Mean   -AllStats_B2_PT.assets_low.(Names_i{1}).Mean))
fprintf('per-type FnsToEvaluate, assets   breakdown (type %s, single==per-type), this should be zero: %.3e \n',Names_i{2},abs(AllStats_B2.assets.(Names_i{2}).Mean   -AllStats_B2_PT.assets_high.(Names_i{2}).Mean))
fprintf('per-type FnsToEvaluate, earnings breakdown (type %s, single==per-type), this should be zero: %.3e \n',Names_i{1},abs(AllStats_B2.earnings.(Names_i{1}).Mean -AllStats_B2_PT.earnings_low.(Names_i{1}).Mean))
fprintf('per-type FnsToEvaluate, earnings breakdown (type %s, single==per-type), this should be zero: %.3e \n',Names_i{2},abs(AllStats_B2.earnings.(Names_i{2}).Mean -AllStats_B2_PT.earnings_high.(Names_i{2}).Mean))

output=struct();

end
