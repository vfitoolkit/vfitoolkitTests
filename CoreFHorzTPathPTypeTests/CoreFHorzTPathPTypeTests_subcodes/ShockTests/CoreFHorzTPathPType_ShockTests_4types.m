function output=CoreFHorzTPathPType_ShockTests_4types(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline)
% Shock-combination cross test for PType TPath. One PType solve with 4 types
% covering every (z,e) combination should give the same per-type VPath,
% PolicyPath and AgentDistPath as 4 separate non-PType TPath solves.
%
% Per-type n_z, z_grid, pi_z, ReturnFn, jequaloneDist, V_final, Policy_final
% and the e pieces of vfoptions/simoptions are all passed as Names_i-keyed
% structures (required when these differ across types).
%
% semiz is not yet part of the TPath tests, so (unlike the core FHorz PType
% ShockTests with 8 types) it is not included here. All types are without d.

Names_i={'none','z','e','ze'};
N_i=length(Names_i);

% Arbitrary unequal type weights
Params.ptypeweights=(1:N_i)'/sum(1:N_i);

n_d=0;
d_grid=[];

% Pull the e pieces out of vfoptionsbaseline for convenience
n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

%% Per-type structs for the single PType call
n_z_PT=struct();      z_grid_PT=struct();      pi_z_PT=struct();
ReturnFn_PT=struct(); jequaloneDist_PT=struct();
V_final_PT=struct();  Policy_final_PT=struct();

vfopt_PT=struct();      simopt_PT=struct();
vfopt_PT.n_e=struct();  vfopt_PT.e_grid=struct();  vfopt_PT.pi_e=struct();
simopt_PT.n_e=struct(); simopt_PT.e_grid=struct(); simopt_PT.pi_e=struct();

% type 1: none
nm='none';
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist_PT.(nm)=zeros(n_a,1,'gpuArray');
jequaloneDist_PT.(nm)(1)=1;
V_final_PT.(nm)=zeros([n_a,N_j],'gpuArray');
Policy_final_PT.(nm)=ones([1,n_a,N_j],'gpuArray');

% type 2: z only
nm='z';
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist_PT.(nm)=zeros(n_a,n_z,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_z/2))=1;
V_final_PT.(nm)=zeros([n_a,n_z,N_j],'gpuArray');
Policy_final_PT.(nm)=ones([1,n_a,n_z,N_j],'gpuArray');

% type 3: e only
nm='e';
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e;simopt_PT.e_grid.(nm)=e_grid;simopt_PT.pi_e.(nm)=pi_e;
jequaloneDist_PT.(nm)=zeros(n_a,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_e/2))=1;
V_final_PT.(nm)=zeros([n_a,n_e,N_j],'gpuArray');
Policy_final_PT.(nm)=ones([1,n_a,n_e,N_j],'gpuArray');

% type 4: z + e
nm='ze';
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e;simopt_PT.e_grid.(nm)=e_grid;simopt_PT.pi_e.(nm)=pi_e;
jequaloneDist_PT.(nm)=zeros(n_a,n_z,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_z/2),ceil(n_e/2))=1;
V_final_PT.(nm)=zeros([n_a,n_z,n_e,N_j],'gpuArray');
Policy_final_PT.(nm)=ones([1,n_a,n_z,n_e,N_j],'gpuArray');

%% PType solve (one call, four types)
[VPath_PT,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z_PT, N_j, Names_i, d_grid, a_grid,z_grid_PT, pi_z_PT, DiscountFactorParamNames, ReturnFn_PT, transpathoptionsbaseline, vfopt_PT);

% First-period policy of each type (the number of dimensions differs by type)
Policy_init_PT=struct();
Policy_init_PT.none=PolicyPath_PT.none(:,:,:,1);
Policy_init_PT.z=PolicyPath_PT.z(:,:,:,:,1);
Policy_init_PT.e=PolicyPath_PT.e(:,:,:,:,1);
Policy_init_PT.ze=PolicyPath_PT.ze(:,:,:,:,:,1);

AgentDist_initial_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist_PT,AgeWeightParamNames,PTypeDistParamNames,Policy_init_PT,n_d,n_a,n_z_PT,N_j,Names_i,pi_z_PT,Params,simopt_PT);
AgentDistPath_PT=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist_PT, PricePath, ParamPath, PolicyPath_PT, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z_PT,N_j,Names_i,pi_z_PT, T,Params, transpathoptionsbaseline, simopt_PT);

%% (PA8) Statistics along the path, with (z,e) different for each ptype
% FnsToEvaluate by ptype (their inputs differ with the shocks). Each ptype's AllStats and life-cycle
% profiles must equal the commands without permanent types on its solo paths (checked in the loop below).
FnsToEvaluate_PT=struct();
FnsToEvaluate_PT.assets.none=@(aprime,a) a;
FnsToEvaluate_PT.assets.z=@(aprime,a,z) a;
FnsToEvaluate_PT.assets.e=@(aprime,a,e) a;
FnsToEvaluate_PT.assets.ze=@(aprime,a,z,e) a;
FnsToEvaluate_PT.earnings.none=@(aprime,a,w,kappa_j) w*kappa_j;
FnsToEvaluate_PT.earnings.z=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate_PT.earnings.e=@(aprime,a,e,w,kappa_j) w*kappa_j*e;
FnsToEvaluate_PT.earnings.ze=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
AllStats_PT=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate_PT, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z_PT, N_j, Names_i, d_grid, a_grid,z_grid_PT, transpathoptionsbaseline, simopt_PT);
LCP_PT=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate_PT, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z_PT, N_j, Names_i, d_grid, a_grid,z_grid_PT, transpathoptionsbaseline, simopt_PT);
OutputTypes={'Mean','Median','Variance','LorenzCurve','Gini','QuantileCutoffs','QuantileMeans'};
fnnames=fieldnames(FnsToEvaluate_PT);

%% Four independent single-type solves
for ii=1:N_i
    nm=Names_i{ii};

    % Per-type vfoptions/simoptions for the non-PType call
    vfopt_solo=struct();
    simopt_solo=struct();
    if isfield(vfopt_PT.n_e,nm)
        vfopt_solo.n_e=vfopt_PT.n_e.(nm);
        vfopt_solo.e_grid=vfopt_PT.e_grid.(nm);
        vfopt_solo.pi_e=vfopt_PT.pi_e.(nm);
        simopt_solo.n_e=vfopt_solo.n_e;
        simopt_solo.e_grid=vfopt_solo.e_grid;
        simopt_solo.pi_e=vfopt_solo.pi_e;
    end

    [VPath_solo,PolicyPath_solo]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_PT.(nm), Policy_final_PT.(nm), Params, n_d, n_a, n_z_PT.(nm), N_j, d_grid, a_grid,z_grid_PT.(nm), pi_z_PT.(nm), DiscountFactorParamNames, ReturnFn_PT.(nm), transpathoptionsbaseline, vfopt_solo);

    fprintf('ShockTests 4 types, VPath        (type %s), this should be zero: %.3e \n',nm,max(abs(VPath_PT.(nm)(:)-VPath_solo(:))))
    fprintf('ShockTests 4 types, PolicyPath   (type %s), this should be zero: %.3e \n',nm,max(abs(PolicyPath_PT.(nm)(:)-PolicyPath_solo(:))))

    % PolicyPath equality has just been checked, so seeding the solo initial
    % distribution from Policy_init_PT.(nm) is the same as using the first
    % period of PolicyPath_solo (and avoids per-type slicing again).
    AgentDist_initial_solo=StationaryDist_FHorz_Case1(jequaloneDist_PT.(nm),AgeWeightParamNames,Policy_init_PT.(nm),n_d,n_a,n_z_PT.(nm),N_j,pi_z_PT.(nm),Params,simopt_solo);
    AgentDistPath_solo=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_solo, jequaloneDist_PT.(nm), PricePath, ParamPath, PolicyPath_solo, AgeWeightParamNames,n_d,n_a,n_z_PT.(nm),N_j,pi_z_PT.(nm), T,Params, transpathoptionsbaseline, simopt_solo);

    fprintf('ShockTests 4 types, AgentDistInit(type %s), this should be zero: %.3e \n',nm,max(abs(AgentDist_initial_PT.(nm)(:)-AgentDist_initial_solo(:))))
    fprintf('ShockTests 4 types, AgentDistPath(type %s), this should be zero: %.3e \n',nm,max(abs(AgentDistPath_PT.(nm)(:)-AgentDistPath_solo(:))))

    % (PA8) this ptype's stats vs the commands without permanent types on its solo paths
    FnsToEvaluate_solo=struct();
    FnsToEvaluate_solo.assets=FnsToEvaluate_PT.assets.(nm);
    FnsToEvaluate_solo.earnings=FnsToEvaluate_PT.earnings.(nm);
    REF=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluate_solo, AgentDistPath_solo, PolicyPath_solo, PricePath, ParamPath, Params, T, n_d, n_a, n_z_PT.(nm), N_j, d_grid, a_grid,z_grid_PT.(nm), transpathoptionsbaseline, simopt_solo);
    REFL=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate_solo, AgentDistPath_solo, PolicyPath_solo, PricePath, ParamPath, Params, T, n_d, n_a, n_z_PT.(nm), N_j, d_grid, a_grid, z_grid_PT.(nm), transpathoptionsbaseline, simopt_solo);
    % every field of the reference, to any depth, found in the PType output with the ptype name after the FnsToEvaluate name
    maxdiff=0; nanmismatch=0;
    queue={};
    REFnames=fieldnames(REF);
    for kk=1:length(REFnames)
        queue{end+1}=REFnames(kk);
    end
    while ~isempty(queue)
        fieldpath=queue{1};
        queue(1)=[];
        O=getfield(REF,fieldpath{:});
        if isstruct(O)
            Onames=fieldnames(O);
            for kk=1:length(Onames)
                queue{end+1}=[fieldpath,Onames(kk)];
            end
        else
            Ppath=[fieldpath(1),{nm},fieldpath(2:end)];
            P=AllStats_PT;
            found=1;
            for kk=1:length(Ppath)
                if isstruct(P) && isfield(P,Ppath{kk})
                    P=P.(Ppath{kk});
                else
                    found=0;
                    break
                end
            end
            if found==0
                nanmismatch=nanmismatch+1;
            elseif iscell(O) % text (the comments): must be the same text
                nanmismatch=nanmismatch+(~isequal(P,O));
            else
                maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
                nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
            end
        end
    end
    fprintf('ShockTests 4 types, (PA8) EvalFnOnTransPath_AllStats_PType vs without permanent types (type %s), every stat, this should be zero: %.3e \n',nm,maxdiff)
    fprintf('ShockTests 4 types, (PA8) EvalFnOnTransPath_AllStats_PType vs without permanent types (type %s), positions NaN (or text) in one but not the other, or missing, this should be zero: %i \n',nm,nanmismatch)
    maxdiff=0; nanmismatch=0;
    for ff=1:length(fnnames)
        for oo=1:length(OutputTypes)
            PI=gather(LCP_PT.(fnnames{ff}).(nm).initialages.(OutputTypes{oo}));
            PB=gather(LCP_PT.(fnnames{ff}).(nm).bornduringtranstion.(OutputTypes{oo}));
            OI=gather(REFL.(fnnames{ff}).initialages.(OutputTypes{oo}));
            OB=gather(REFL.(fnnames{ff}).bornduringtranstion.(OutputTypes{oo}));
            maxdiff=max([maxdiff; abs(PI(:)-OI(:)); abs(PB(:)-OB(:))],[],'omitnan');
            nanmismatch=nanmismatch+sum(isnan(PI(:))~=isnan(OI(:)))+sum(isnan(PB(:))~=isnan(OB(:)));
        end
    end
    fprintf('ShockTests 4 types, (PA8) LifeCycleProfiles_TransPath_PType vs without permanent types (type %s), all seven outputs, this should be zero: %.3e \n',nm,maxdiff)
    fprintf('ShockTests 4 types, (PA8) LifeCycleProfiles_TransPath_PType vs without permanent types (type %s), entries NaN in one but not the other, this should be zero: %i \n',nm,nanmismatch)
end

output=struct();

end
