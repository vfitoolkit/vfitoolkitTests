function output=CoreFHorzPType_SharedPiZNzEqNi(n_a,N_j,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: n_z equal to N_i (both 2), where the last dimension of pi_z has length N_i whether or not pi_z depends on ptype.
%   shared        pi_z [n_z,n_z]       a 2D pi_z is shared by the types; it used to be sliced as if its last dimension were ptype
%   struct        pi_z.(name)          the same matrix for both types, as a struct
%   trailing dim  pi_z [n_z,n_z,N_i]   a different matrix per type: still sliced by ptype
% Per-type V, Policy, StationaryDist and AllStats Mean must equal solo solves with that type's pi_z.
% The types differ by kappa_j_pt; in the trailing-dim case they also differ by pi_z (checked to change V,
% so that reading the wrong type's pi_z would show).

N_i=2;
n_z=2;
n_d=0;
d_grid=[];

z_grid=gpuArray([0.8;1.2]);
pi_z_shared=gpuArray([0.9,0.1;0.3,0.7]); % not symmetric, so a transposed or sliced read changes the answer
pi_z_pt=cat(3,gpuArray([0.9,0.1;0.3,0.7]),gpuArray([0.6,0.4;0.2,0.8])); % [n_z,n_z,N_i]

ReturnFn_PT=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);
ReturnFn_NoPT=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate_PT.assets=@(aprime,a,z) a;
FnsToEvaluate_PT.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate_NoPT.assets=@(aprime,a,z) a;
FnsToEvaluate_NoPT.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnNames=fieldnames(FnsToEvaluate_NoPT);

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,1)=1;

vfoptions=struct();
simoptions=struct();

%% Solo solves: each type with the shared pi_z, and with its own pi_z
V_solo=cell(N_i,2); Policy_solo=cell(N_i,2); Dist_solo=cell(N_i,2); AllStats_solo=cell(N_i,2);
for ii=1:N_i
    Params_ii=Params;
    Params_ii.kappa_j=Params.kappa_j_pt(ii,:);
    for pc=1:2
        if pc==1
            pi_z_ii=pi_z_shared;
        else
            pi_z_ii=pi_z_pt(:,:,ii);
        end
        [V_solo{ii,pc},Policy_solo{ii,pc}]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_ii,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,[],vfoptions);
        Dist_solo{ii,pc}=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo{ii,pc},n_d,n_a,n_z,N_j,pi_z_ii,Params_ii,simoptions);
        AllStats_solo{ii,pc}=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist_solo{ii,pc},Policy_solo{ii,pc},FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    end
end

%% The three PType forms
casenames={'shared 2D pi_z','shared pi_z as struct','trailing-dim pi_z'};
casepc=[1,1,2]; % which solo solves each case is compared to
for cc=1:3
    if cc==1
        pi_z_PT=pi_z_shared;
    elseif cc==2
        pi_z_PT=struct();
        for ii=1:N_i
            pi_z_PT.(['ptype',sprintf('%03i',ii)])=pi_z_shared;
        end
    elseif cc==3
        pi_z_PT=pi_z_pt;
    end
    [V_PT,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z_PT,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
    StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,N_i,pi_z_PT,Params,simoptions);
    AllStats_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
    names_PT=fieldnames(V_PT);
    pc=casepc(cc);
    for ii=1:N_i
        nA=names_PT{ii};
        fprintf('n_z==N_i (%s), V    (type %d), this should be zero: %.3e \n',casenames{cc},ii,max(abs(V_PT.(nA)(:)-V_solo{ii,pc}(:))))
        fprintf('n_z==N_i (%s), Pol  (type %d), this should be zero: %.3e \n',casenames{cc},ii,max(abs(Policy_PT.(nA)(:)-Policy_solo{ii,pc}(:))))
        fprintf('n_z==N_i (%s), Dist (type %d), this should be zero: %.3e \n',casenames{cc},ii,max(abs(StationaryDist_PT.(nA)(:)-Dist_solo{ii,pc}(:))))
        for ff=1:length(FnNames)
            fprintf('n_z==N_i (%s), AllStats %s Mean (type %d), this should be zero: %.3e \n',casenames{cc},FnNames{ff},ii,abs(AllStats_PT.(FnNames{ff}).(nA).Mean-AllStats_solo{ii,pc}.(FnNames{ff}).Mean))
        end
    end
end

%% The per-type pi_z really differ (else the trailing-dim case cannot tell a wrong slice), and so do the types
fprintf('n_z==N_i, V with type 2 own pi_z differs from shared pi_z, this should NOT be zero: %.3e \n',max(abs(V_solo{2,1}(:)-V_solo{2,2}(:))))
fprintf('n_z==N_i, V of the two types differ (shared pi_z), this should NOT be zero: %.3e \n',max(abs(V_solo{1,1}(:)-V_solo{2,1}(:))))

output=struct();

end
