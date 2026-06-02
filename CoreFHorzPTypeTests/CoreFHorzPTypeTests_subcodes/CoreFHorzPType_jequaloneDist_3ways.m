function output=CoreFHorzPType_jequaloneDist_3ways(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i)
% Test: jequaloneDist accepts three input shapes for PType StationaryDist:
%   (i)   [n_a,n_z]            — same dist for every PType
%   (ii)  [n_a,n_z,N_i]        — last dim indexes PTypes; slice masses override PTypeDistParamNames
%   (iii) struct with Names_i  — one [n_a,n_z] field (mass 1) per type; uses PTypeDistParamNames
% V and Policy do not depend on jequaloneDist, so only StationaryDist and
% downstream AllStats are compared. Per-type kappa_j_pt makes the two types
% genuinely different (so the comparison actually exercises per-type setup).

n_d=0;
d_grid=[];

ReturnFn_PT=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;

vfoptions=struct();
simoptions=struct();

ptw=Params.(PTypeDistParamNames{:}); % [0.6; 0.4] from setup

%% Solve V/Policy once (independent of jequaloneDist)
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);

%% Part A: three semantically identical jequaloneDists (same dist for every type)
jd_base=zeros(n_a,n_z,'gpuArray');
jd_base(1,ceil(n_z/2))=1;

% (i) 2D matrix, mass 1; uses ptypeweights from Params
jd_i=jd_base;

% (ii) 3D matrix; slice masses become the ptypeweights (override). Set them
% to match Params.ptypeweights so the result agrees with (i) and (iii).
jd_ii=zeros(n_a,n_z,N_i,'gpuArray');
for ii=1:N_i
    jd_ii(:,:,ii)=ptw(ii)*jd_base;
end

% (iii) struct, each field mass 1; uses ptypeweights from Params
jd_iii=struct();
for ii=1:N_i
    jd_iii.(Names_i{ii})=jd_base;
end

Dist_i  =StationaryDist_Case1_FHorz_PType(jd_i,  AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
Dist_ii =StationaryDist_Case1_FHorz_PType(jd_ii, AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
Dist_iii=StationaryDist_Case1_FHorz_PType(jd_iii,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);

AllStats_i  =EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(Dist_i,  Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
AllStats_ii =EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(Dist_ii, Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
AllStats_iii=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(Dist_iii,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);

for ii=1:N_i
    nA=Names_i{ii};
    fprintf('jequaloneDist 3 ways, Dist (type %s) (i)  vs (ii),  this should be zero: %2.8f \n',nA,max(abs(Dist_i.(nA)(:) -Dist_ii.(nA)(:))))
    fprintf('jequaloneDist 3 ways, Dist (type %s) (i)  vs (iii), this should be zero: %2.8f \n',nA,max(abs(Dist_i.(nA)(:) -Dist_iii.(nA)(:))))
    fprintf('jequaloneDist 3 ways, Dist (type %s) (ii) vs (iii), this should be zero: %2.8f \n',nA,max(abs(Dist_ii.(nA)(:)-Dist_iii.(nA)(:))))
end
fprintf('jequaloneDist 3 ways, AllStats assets.Mean   (i) vs (ii),  this should be zero: %2.8f \n',abs(AllStats_i.assets.Mean   -AllStats_ii.assets.Mean))
fprintf('jequaloneDist 3 ways, AllStats assets.Mean   (i) vs (iii), this should be zero: %2.8f \n',abs(AllStats_i.assets.Mean   -AllStats_iii.assets.Mean))
fprintf('jequaloneDist 3 ways, AllStats earnings.Mean (i) vs (ii),  this should be zero: %2.8f \n',abs(AllStats_i.earnings.Mean -AllStats_ii.earnings.Mean))
fprintf('jequaloneDist 3 ways, AllStats earnings.Mean (i) vs (iii), this should be zero: %2.8f \n',abs(AllStats_i.earnings.Mean -AllStats_iii.earnings.Mean))

%% Part B: per-type-different bases — approach (i) cannot represent this; check (ii) vs (iii)
jdA=zeros(n_a,n_z,'gpuArray'); jdA(1,ceil(n_z/2))=1;
jdB=zeros(n_a,n_z,'gpuArray'); jdB(1,1)=1;

jd_ii_diff=zeros(n_a,n_z,N_i,'gpuArray');
jd_ii_diff(:,:,1)=ptw(1)*jdA;
jd_ii_diff(:,:,2)=ptw(2)*jdB;

jd_iii_diff=struct();
jd_iii_diff.(Names_i{1})=jdA;
jd_iii_diff.(Names_i{2})=jdB;

Dist_ii_d =StationaryDist_Case1_FHorz_PType(jd_ii_diff, AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
Dist_iii_d=StationaryDist_Case1_FHorz_PType(jd_iii_diff,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);

AllStats_ii_d =EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(Dist_ii_d, Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
AllStats_iii_d=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(Dist_iii_d,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);

for ii=1:N_i
    nA=Names_i{ii};
    fprintf('jequaloneDist 3 ways (per-type diff), Dist (type %s) (ii) vs (iii), this should be zero: %2.8f \n',nA,max(abs(Dist_ii_d.(nA)(:)-Dist_iii_d.(nA)(:))))
end
fprintf('jequaloneDist 3 ways (per-type diff), AllStats assets.Mean   (ii) vs (iii), this should be zero: %2.8f \n',abs(AllStats_ii_d.assets.Mean   -AllStats_iii_d.assets.Mean))
fprintf('jequaloneDist 3 ways (per-type diff), AllStats earnings.Mean (ii) vs (iii), this should be zero: %2.8f \n',abs(AllStats_ii_d.earnings.Mean -AllStats_iii_d.earnings.Mean))

output=struct();

end
