function output=CoreFHorzPType_NiAllIdentical_vsNoPType(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i)
% Test 3: PType solve with N_i types that are all identical (no per-type
% parameter dependence) should give the same per-type V, Policy, and (conditional)
% StationaryDist as a single no-PType solve. The aggregated AllStats should also
% match. Holds for arbitrary ptypeweights.

n_d=0;
d_grid=[];

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

%% PType solve, N_i types, all identical
[V_PT,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
V_PT_vfp=ValueFnFromPolicy_FHorz_PType(Policy_PT,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
AllStats_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);

%% Single no-PType solve
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
V0_vfp=ValueFnFromPolicy_FHorz(Policy0,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy0,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
AllStats0=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist0,Policy0,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

%% Compare per type
names=fieldnames(V_PT);
for ii=1:N_i
    nA=names{ii};
    fprintf('Identical PType vs no-PType, V    (type %d), this should be zero: %2.8f \n',ii,max(abs(V_PT.(nA)(:)-V0(:))))
    fprintf('Identical PType vs no-PType, Pol  (type %d), this should be zero: %2.8f \n',ii,max(abs(Policy_PT.(nA)(:)-Policy0(:))))
    fprintf('Identical PType vs no-PType, Dist (type %d), this should be zero: %2.8f \n',ii,max(abs(StationaryDist_PT.(nA)(:)-StationaryDist0(:))))
    fprintf('Identical PType vs no-PType, VFP  vs V (PType,   type %d), this should be zero: %2.8f \n',ii,max(abs(V_PT_vfp.(nA)(:)-V_PT.(nA)(:))))
end
fprintf('Identical PType vs no-PType, VFP  vs V (noPType),         this should be zero: %2.8f \n',max(abs(V0_vfp(:)-V0(:))))

%% Aggregated AllStats should equal the no-PType AllStats
fprintf('Identical PType vs no-PType, AllStats assets.Mean,    this should be zero: %2.8f \n',abs(AllStats_PT.assets.Mean   -AllStats0.assets.Mean))
fprintf('Identical PType vs no-PType, AllStats earnings.Mean,  this should be zero: %2.8f \n',abs(AllStats_PT.earnings.Mean -AllStats0.earnings.Mean))
fprintf('Identical PType vs no-PType, AllStats earnings.Gini,  this should be zero: %2.8f \n',abs(AllStats_PT.earnings.Gini -AllStats0.earnings.Gini))

output=struct();

end
