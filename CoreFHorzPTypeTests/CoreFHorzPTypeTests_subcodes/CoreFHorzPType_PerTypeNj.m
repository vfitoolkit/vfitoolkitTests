function output=CoreFHorzPType_PerTypeNj(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames)
% Test: PType with per-type N_j (different lifespans) should equal stacking
% individual non-PType solves, each run with that type's N_j and the matching
% age-length-dependent parameters. Exercises the struct-form of N_j, the
% per-type slicing of age-dependent params (kappa_j, mewj, agej), and the
% struct-form of jequaloneDist (required when N_j differs across types).
%
% Types differ only in lifespan:
%   long : N_j      (baseline)
%   short: N_j - 5

Names_i={'long','short'};
N_i=length(Names_i);

n_d=0;
d_grid=[];

% Per-type lifespans
N_j_PT.long  = N_j;
N_j_PT.short = N_j-5;

% Age-length-dependent params must be per-type structs because lengths differ
kappa_j_base=Params.kappa_j;
Params.kappa_j=struct();
Params.kappa_j.long  = kappa_j_base;
Params.kappa_j.short = kappa_j_base(1:N_j_PT.short);
Params.mewj=struct();
Params.mewj.long  = ones(1,N_j_PT.long )/N_j_PT.long;
Params.mewj.short = ones(1,N_j_PT.short)/N_j_PT.short;
Params.agej=struct();
Params.agej.long  = 1:N_j_PT.long;
Params.agej.short = 1:N_j_PT.short;
AgeWeightParamNames={'mewj'};

% Per-type weights (overrides what was set in setup)
Params.ptypeweights=[0.6; 0.4];

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;

% jequaloneDist must be a struct because state-space-time shapes differ across types
jequaloneDist_mat=zeros(n_a,n_z,'gpuArray');
jequaloneDist_mat(1,ceil(n_z/2))=1;
jequaloneDist.long  = jequaloneDist_mat;
jequaloneDist.short = jequaloneDist_mat;

vfoptions=struct();
simoptions=struct();

%% PType solve with per-type N_j
[V_PT,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j_PT,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
V_PT_vfp=ValueFnFromPolicy_FHorz_PType(Policy_PT,n_d,n_a,n_z,N_j_PT,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j_PT,Names_i,pi_z,Params,simoptions);
AllStats_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,N_j_PT,Names_i,d_grid,a_grid,z_grid,simoptions);

%% N_i individual non-PType solves, each with its own N_j and age-dependent params
V_solo=struct(); Policy_solo=struct(); Dist_solo=struct(); AllStats_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Params_ii=Params;
    Params_ii.kappa_j=Params.kappa_j.(nm);
    Params_ii.mewj   =Params.mewj.(nm);
    Params_ii.agej   =Params.agej.(nm);
    [V_solo.(nm),Policy_solo.(nm)]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j_PT.(nm),d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist.(nm),AgeWeightParamNames,Policy_solo.(nm),n_d,n_a,n_z,N_j_PT.(nm),pi_z,Params_ii,simoptions);
    AllStats_solo.(nm)=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,N_j_PT.(nm),d_grid,a_grid,z_grid,simoptions);
end

%% Compare per-type slices
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('PerTypeNj, V    (type %s, N_j=%d), this should be zero: %2.8f \n',nm,N_j_PT.(nm),max(abs(V_PT.(nm)(:)-V_solo.(nm)(:))))
    fprintf('PerTypeNj, Pol  (type %s, N_j=%d), this should be zero: %2.8f \n',nm,N_j_PT.(nm),max(abs(Policy_PT.(nm)(:)-Policy_solo.(nm)(:))))
    fprintf('PerTypeNj, Dist (type %s, N_j=%d), this should be zero: %2.8f \n',nm,N_j_PT.(nm),max(abs(StationaryDist_PT.(nm)(:)-Dist_solo.(nm)(:))))
    fprintf('PerTypeNj, VFP  vs V (type %s, N_j=%d), this should be zero: %2.8f \n',nm,N_j_PT.(nm),max(abs(V_PT_vfp.(nm)(:)-V_PT.(nm)(:))))
end

%% Aggregated AllStats.Mean equals ptw-weighted combination of solo means
ptw=Params.(PTypeDistParamNames{:});
agg_assets_mean   = ptw(1)*AllStats_solo.long.assets.Mean   + ptw(2)*AllStats_solo.short.assets.Mean;
agg_earnings_mean = ptw(1)*AllStats_solo.long.earnings.Mean + ptw(2)*AllStats_solo.short.earnings.Mean;
fprintf('PerTypeNj, AllStats assets.Mean,   this should be zero: %2.8f \n',abs(AllStats_PT.assets.Mean   -agg_assets_mean))
fprintf('PerTypeNj, AllStats earnings.Mean, this should be zero: %2.8f \n',abs(AllStats_PT.earnings.Mean -agg_earnings_mean))

output=struct();

end
