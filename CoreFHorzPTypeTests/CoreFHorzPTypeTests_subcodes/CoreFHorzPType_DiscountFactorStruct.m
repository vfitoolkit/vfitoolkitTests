function output=CoreFHorzPType_DiscountFactorStruct(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,AgeWeightParamNames,PTypeDistParamNames)
% Test: DiscountFactorParamNames given as a struct keyed by Names_i, so the two types discount
% with different LISTS of parameters: 'low' with {'beta'}, 'high' with {'beta','sj'}, where sj is
% an age-dependent survival probability. (Part 1 covers a per-type VALUE of beta, by a vector;
% this is the other route, a per-type set of names.)
% Per-type V, Policy, StationaryDist and ValueFnFromPolicy must equal solo solves of each type
% with its own DiscountFactorParamNames.

Names_i={'low','high'};
N_i=length(Names_i);

n_d=0;
d_grid=[];

Params.sj=[0.99*ones(1,N_j-5), 0.9*ones(1,5)];

DiscountFactorParamNames_PT=struct();
DiscountFactorParamNames_PT.low={'beta'};
DiscountFactorParamNames_PT.high={'beta','sj'};

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

%% PType
[V_PT,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames_PT,vfoptions);
V_PT_vfp=ValueFnFromPolicy_FHorz_PType(Policy_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames_PT,vfoptions);
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

%% Solo solves
V_solo=struct(); Policy_solo=struct(); Dist_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};
    [V_solo.(nm),Policy_solo.(nm)]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames_PT.(nm),[],vfoptions);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo.(nm),n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
end

%% Compare
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('DiscountFactorParamNames struct, V    (type %s), this should be zero: %.3e \n',nm,max(abs(V_PT.(nm)(:)-V_solo.(nm)(:))))
    fprintf('DiscountFactorParamNames struct, Pol  (type %s), this should be zero: %.3e \n',nm,max(abs(Policy_PT.(nm)(:)-Policy_solo.(nm)(:))))
    fprintf('DiscountFactorParamNames struct, Dist (type %s), this should be zero: %.3e \n',nm,max(abs(StationaryDist_PT.(nm)(:)-Dist_solo.(nm)(:))))
    fprintf('DiscountFactorParamNames struct, VFP  vs V (type %s), this should be zero: %.3e \n',nm,max(abs(V_PT_vfp.(nm)(:)-V_PT.(nm)(:))))
end
% sj must actually matter, or the two types are the same model
fprintf('DiscountFactorParamNames struct, V of the two types differ, this should NOT be zero: %.3e \n',max(abs(V_solo.low(:)-V_solo.high(:))))

output=struct();

end
