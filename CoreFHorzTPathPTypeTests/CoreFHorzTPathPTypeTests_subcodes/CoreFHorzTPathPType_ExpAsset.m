function output=CoreFHorzTPathPType_ExpAsset(T,PricePath,ParamPath,n_z,N_j,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline)
% Non-standard endogenous states and the grid interpolation layer, by permanent type, on the transition path.
% Three types, each its own model (inputs by Names_i structures throughout):
%   typeA: an experience asset (vfoptions.experienceasset), whose aprimeFn parameter phi1 is on the ParamPath
%          (a structure by ptype; the other types have it constant on the path)
%   typeB: an experience asset that depends on z (vfoptions.experienceassetz), with z_grid and pi_z varying
%          over the path (z_grid.typeB [n_z,T], pi_z.typeB [n_z,n_z,T])
%   typeC: the model with just a (no d), solved with the grid interpolation layer (vfoptions.gridinterplayer)
% The models are those of CoreFHorzTPathExpAssetTests, CoreFHorzTPathExpAssetzTests and CoreFHorzTPathTests,
% on smaller grids. Each type must be exactly what the commands without permanent types give with that type's
% inputs.
%   (X1) ValueFnOnTransPath and AgentDistOnTransPath (PType) vs solo, by ptype, fastOLG=0 and fastOLG=1; and the
%        PType paths with fastOLG=0 vs fastOLG=1 (an independent check: the fastOLG=0 agent dist with the grid
%        interpolation layer is code of its own, which until 2026-10-27 did not exist)
%   (X2) EvalFnOnTransPath_AggVars (PType) vs solo, by ptype
%   (X3) EvalFnOnTransPath_AllStats and LifeCycleProfiles_TransPath (PType) vs solo, by ptype
%   (X4) general eqm: an interest rate and a general eqm condition by ptype (GEptype), which depends on
%        earnings: shooting at maxiter=2 vs solo, by ptype, fastOLG=1 and 0; Anderson with memory=0 vs
%        shooting at maxiter=3 (and with the default memory it must differ)
%   (X5) AgentDistOnTransPath_Case1_FHorz (no permanent types) with the grid interpolation layer, simoptions.fastOLG=0
%        vs fastOLG=1 (value fn and agent dist paths, transpathoptions.fastOLG the same), for each of no z/no e, z, e,
%        z and e, without and with an experience asset: every branch of the fastOLG=0 grid interpolation code
%   and, for scale, that the phi1 path moves typeA's agent dist, and the z path typeB's (must not be zero)

OutputTypes={'Mean','Median','Variance','LorenzCurve','Gini','QuantileCutoffs','QuantileMeans'};

Names_i={'typeA','typeB','typeC'};
N_i=3;
PTypeDistParamNames={'ptweights_ea'};
Params.ptweights_ea=[0.3;0.3;0.4];
Params.phi1=0.3; % ratio at which (1-d2) (times z for typeB) is converted into a2
Params.phi2=0.03; % depreciation rate of a2

%% Grids
n_a1=51;
n_a2=9;
a1_grid=5*linspace(0,1,n_a1)'.^3;
a2_grid=linspace(0,10,n_a2)';
n_d2=3;
d2_grid=linspace(0,1,n_d2)';

driftB=1+0.1*(1:1:T)/T; % [1,T]
pimixB=0.3*(1:1:T)/T;
z_gridT_B=z_grid.*driftB; % [n_z,T]
pi_zT_B=(1-reshape(pimixB,[1,1,T])).*pi_z+reshape(pimixB,[1,1,T]).*ones(n_z,n_z)/n_z; % [n_z,n_z,T]
phi1Path=0.3+0.1*min((0:T-1)/(T/2),1); % rises from 0.3 to 0.4 over the first half of the path

%% Each type's inputs
n_d_PT=struct('typeA',n_d2,'typeB',n_d2,'typeC',0);
d_grid_PT=struct('typeA',d2_grid,'typeB',d2_grid,'typeC',[]);
n_a_PT=struct('typeA',[n_a1,n_a2],'typeB',[n_a1,n_a2],'typeC',n_a1);
a_grid_PT=struct('typeA',[a1_grid;a2_grid],'typeB',[a1_grid;a2_grid],'typeC',a1_grid);
z_grid_PT=struct('typeA',z_grid,'typeB',z_gridT_B,'typeC',z_grid);
pi_z_PT=struct('typeA',pi_z,'typeB',pi_zT_B,'typeC',pi_z);
ReturnFn_PT=struct();
ReturnFn_PT.typeA=@(d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_noe_nosemiz(d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_PT.typeB=@(d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_TPathExpAssetz_nod1_z_noe_nosemiz(d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_PT.typeC=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate=struct();
FnsToEvaluate.assets.typeA=@(d2,a1prime,a1,a2,z) a1+a2;
FnsToEvaluate.assets.typeB=@(d2,a1prime,a1,a2,z) a1+a2;
FnsToEvaluate.assets.typeC=@(aprime,a,z) a;
FnsToEvaluate.earnings.typeA=@(d2,a1prime,a1,a2,z,w,kappa_j) w*kappa_j*d2*a2*z;
FnsToEvaluate.earnings.typeB=@(d2,a1prime,a1,a2,z,w,kappa_j) w*kappa_j*d2*a2*z;
FnsToEvaluate.earnings.typeC=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
fnnames=fieldnames(FnsToEvaluate);
jequaloneDist_PT=struct();
jequaloneDist_PT.typeA=zeros([n_a1,n_a2,n_z],'gpuArray');
jequaloneDist_PT.typeA(1,1,ceil(n_z/2))=1;
jequaloneDist_PT.typeB=jequaloneDist_PT.typeA;
jequaloneDist_PT.typeC=zeros([n_a1,n_z],'gpuArray');
jequaloneDist_PT.typeC(1,ceil(n_z/2))=1;

% vfoptions and simoptions, each option a structure by ptype
vfoptions=struct();
vfoptions.experienceasset=struct('typeA',1,'typeB',0,'typeC',0);
vfoptions.experienceassetz=struct('typeA',0,'typeB',1,'typeC',0);
vfoptions.aprimeFn=struct();
vfoptions.aprimeFn.typeA=@(d2,a2,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;
vfoptions.aprimeFn.typeB=@(d2,a2,z,phi1,phi2) phi1*(1-d2)*z+(1-phi2)*a2;
vfoptions.gridinterplayer=struct('typeA',0,'typeB',0,'typeC',1);
vfoptions.ngridinterp=struct('typeC',5);
simoptions=vfoptions;
simoptions.d_grid=struct('typeA',d2_grid,'typeB',d2_grid);
simoptions.a_grid=struct('typeA',[a1_grid;a2_grid],'typeB',[a1_grid;a2_grid]);
simoptions.z_grid=struct('typeB',z_gridT_B); % the aprimeFn of typeB depends on z, so the sim codes need z_grid
ParamPathE=ParamPath;
ParamPathE.phi1=struct('typeA',phi1Path,'typeB',Params.phi1*ones(1,T),'typeC',Params.phi1*ones(1,T));

% The same inputs for each type on its own, without permanent types
vfopt_solo=cell(N_i,1); simopt_solo=cell(N_i,1); ParamPath_solo=cell(N_i,1);
for ii=1:N_i
    nm=Names_i{ii};
    vfopt_solo{ii}=struct();
    vfopt_solo{ii}.experienceasset=vfoptions.experienceasset.(nm);
    vfopt_solo{ii}.experienceassetz=vfoptions.experienceassetz.(nm);
    vfopt_solo{ii}.gridinterplayer=vfoptions.gridinterplayer.(nm);
    if isfield(vfoptions.aprimeFn,nm)
        vfopt_solo{ii}.aprimeFn=vfoptions.aprimeFn.(nm);
    end
    if isfield(vfoptions.ngridinterp,nm)
        vfopt_solo{ii}.ngridinterp=vfoptions.ngridinterp.(nm);
    end
    simopt_solo{ii}=vfopt_solo{ii};
    if isfield(simoptions.d_grid,nm)
        simopt_solo{ii}.d_grid=simoptions.d_grid.(nm);
        simopt_solo{ii}.a_grid=simoptions.a_grid.(nm);
    end
    if isfield(simoptions.z_grid,nm)
        simopt_solo{ii}.z_grid=simoptions.z_grid.(nm);
    end
    ParamPath_solo{ii}=ParamPath;
    ParamPath_solo{ii}.phi1=ParamPathE.phi1.(nm);
end

%% Initial agent dist and V_final of each type: the stationary ones with the plain z_grid, at Params
V_final_PT=struct(); Policy_final_PT=struct(); AgentDist_initial_PT=struct();
for ii=1:N_i
    nm=Names_i{ii};
    simopt0=simopt_solo{ii};
    if isfield(simopt0,'z_grid')
        simopt0.z_grid=z_grid;
    end
    [V_final_PT.(nm),Policy_final_PT.(nm)]=ValueFnIter_Case1_FHorz(n_d_PT.(nm),n_a_PT.(nm),n_z,N_j,d_grid_PT.(nm),a_grid_PT.(nm),z_grid,pi_z,ReturnFn_PT.(nm),Params,DiscountFactorParamNames,[],vfopt_solo{ii});
    AgentDist_initial_PT.(nm)=StationaryDist_FHorz_Case1(jequaloneDist_PT.(nm),AgeWeightParamNames,Policy_final_PT.(nm),n_d_PT.(nm),n_a_PT.(nm),n_z,N_j,pi_z,Params,simopt0);
end
AgentDist_initial_PT.ptweights=Params.ptweights_ea;

%% (X1)/(X2) the value fn and agent dist paths, and AggVars
PolicyPath_PT0=[]; AgentDistPath_PT0=[];
for fastOLG=[0,1]
    transpathoptions=struct();
    transpathoptions.fastOLG=fastOLG;
    simoptionsF=simoptions;
    simoptionsF.fastOLG=fastOLG;
    [VPath_PT,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPathE, T, V_final_PT, Policy_final_PT, Params, n_d_PT, n_a_PT, n_z, N_j, Names_i, d_grid_PT, a_grid_PT,z_grid_PT, pi_z_PT, DiscountFactorParamNames, ReturnFn_PT, transpathoptions, vfoptions);
    AgentDistPath_PT=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist_PT, PricePath, ParamPathE, PolicyPath_PT, AgeWeightParamNames, PTypeDistParamNames,n_d_PT,n_a_PT,n_z,N_j,Names_i,pi_z_PT, T,Params, transpathoptions, simoptionsF);
    AggVarsPath_PT=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPathE, Params, T, n_d_PT, n_a_PT, n_z, N_j, Names_i, d_grid_PT, a_grid_PT,z_grid_PT, transpathoptions, simoptionsF);
    for ii=1:N_i
        nm=Names_i{ii};
        simopt_soloF=simopt_solo{ii};
        simopt_soloF.fastOLG=fastOLG;
        FnsToEvaluate_solo=struct('assets',FnsToEvaluate.assets.(nm),'earnings',FnsToEvaluate.earnings.(nm));
        [VPath_solo,PolicyPath_solo]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath_solo{ii}, T, V_final_PT.(nm), Policy_final_PT.(nm), Params, n_d_PT.(nm), n_a_PT.(nm), n_z, N_j, d_grid_PT.(nm), a_grid_PT.(nm),z_grid_PT.(nm), pi_z_PT.(nm), DiscountFactorParamNames, ReturnFn_PT.(nm), transpathoptions, vfopt_solo{ii});
        AgentDistPath_solo=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_PT.(nm), jequaloneDist_PT.(nm), PricePath, ParamPath_solo{ii}, PolicyPath_solo, AgeWeightParamNames,n_d_PT.(nm),n_a_PT.(nm),n_z,N_j,pi_z_PT.(nm), T,Params, transpathoptions, simopt_soloF);
        AggVarsPath_solo=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate_solo, AgentDistPath_solo, PolicyPath_solo, PricePath, ParamPath_solo{ii}, Params, T, n_d_PT.(nm), n_a_PT.(nm), n_z, N_j, d_grid_PT.(nm), a_grid_PT.(nm),z_grid_PT.(nm), transpathoptions, simopt_soloF);
        fprintf('ExpAsset/GI by ptype (fastOLG=%i): (X1) PType vs solo (type %s), VPath, this should be zero: %.3e \n',fastOLG,nm,gather(max(abs(VPath_PT.(nm)(:)-VPath_solo(:)))))
        fprintf('ExpAsset/GI by ptype (fastOLG=%i): (X1) PType vs solo (type %s), PolicyPath, this should be zero: %.3e \n',fastOLG,nm,gather(max(abs(PolicyPath_PT.(nm)(:)-PolicyPath_solo(:)))))
        fprintf('ExpAsset/GI by ptype (fastOLG=%i): (X1) PType vs solo (type %s), AgentDistPath, this should be zero: %.3e \n',fastOLG,nm,gather(max(abs(AgentDistPath_PT.(nm)(:)-AgentDistPath_solo(:)))))
        fprintf('ExpAsset/GI by ptype (fastOLG=%i): (X2) PType vs solo (type %s), AggVars Mean of assets and earnings, this should be zero: %.3e \n',fastOLG,nm,max([max(abs(gather(AggVarsPath_PT.assets.(nm).Mean(:)-AggVarsPath_solo.assets.Mean(:)))),max(abs(gather(AggVarsPath_PT.earnings.(nm).Mean(:)-AggVarsPath_solo.earnings.Mean(:))))]))
        if fastOLG==0 && ii==1
            % For scale: the phi1 path must move typeA (the plain run is typeA's solo with phi1 constant)
            [~,PolicyPath_c]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_PT.(nm), Policy_final_PT.(nm), Params, n_d_PT.(nm), n_a_PT.(nm), n_z, N_j, d_grid_PT.(nm), a_grid_PT.(nm),z_grid_PT.(nm), pi_z_PT.(nm), DiscountFactorParamNames, ReturnFn_PT.(nm), transpathoptions, vfopt_solo{ii});
            AgentDistPath_c=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_PT.(nm), jequaloneDist_PT.(nm), PricePath, ParamPath, PolicyPath_c, AgeWeightParamNames,n_d_PT.(nm),n_a_PT.(nm),n_z,N_j,pi_z_PT.(nm), T,Params, transpathoptions, simopt_soloF);
            fprintf('ExpAsset/GI by ptype: for scale, the phi1 path moves the typeA agent dist path, this should NOT be zero: %.3e \n',gather(max(abs(AgentDistPath_solo(:)-AgentDistPath_c(:)))))
            clear PolicyPath_c AgentDistPath_c
        end
        if fastOLG==0 && ii==2
            % For scale: the z path must move typeB (the plain run is typeB's solo with the plain z_grid and pi_z)
            simopt_c=simopt_soloF;
            simopt_c.z_grid=z_grid;
            [~,PolicyPath_c]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath_solo{ii}, T, V_final_PT.(nm), Policy_final_PT.(nm), Params, n_d_PT.(nm), n_a_PT.(nm), n_z, N_j, d_grid_PT.(nm), a_grid_PT.(nm),z_grid, pi_z, DiscountFactorParamNames, ReturnFn_PT.(nm), transpathoptions, vfopt_solo{ii});
            AgentDistPath_c=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_PT.(nm), jequaloneDist_PT.(nm), PricePath, ParamPath_solo{ii}, PolicyPath_c, AgeWeightParamNames,n_d_PT.(nm),n_a_PT.(nm),n_z,N_j,pi_z, T,Params, transpathoptions, simopt_c);
            fprintf('ExpAsset/GI by ptype: for scale, the z path moves the typeB agent dist path, this should NOT be zero: %.3e \n',gather(max(abs(AgentDistPath_solo(:)-AgentDistPath_c(:)))))
            clear PolicyPath_c AgentDistPath_c
        end
    end
    if fastOLG==0
        PolicyPath_PT0=PolicyPath_PT;
        AgentDistPath_PT0=AgentDistPath_PT;
    else
        for ii=1:N_i
            nm=Names_i{ii};
            fprintf('ExpAsset/GI by ptype: (X1) fastOLG=1 vs fastOLG=0 (type %s), PolicyPath, this should be zero: %.3e \n',nm,gather(max(abs(PolicyPath_PT.(nm)(:)-PolicyPath_PT0.(nm)(:)))))
            fprintf('ExpAsset/GI by ptype: (X1) fastOLG=1 vs fastOLG=0 (type %s), AgentDistPath, this should be zero: %.3e \n',nm,gather(max(abs(AgentDistPath_PT.(nm)(:)-AgentDistPath_PT0.(nm)(:)))))
        end
    end
    clear VPath_PT VPath_solo PolicyPath_solo AgentDistPath_solo
end

%% (X3) the stats, by ptype, vs solo
transpathoptions=struct();
transpathoptions.fastOLG=0;
simoptionsF=simoptions;
simoptionsF.fastOLG=0;
AllStats_PT=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_PT0, PolicyPath_PT0, PricePath, ParamPathE, Params, T, n_d_PT, n_a_PT, n_z, N_j, Names_i, d_grid_PT, a_grid_PT,z_grid_PT, transpathoptions, simoptionsF);
LCP_PT=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate, AgentDistPath_PT0, PolicyPath_PT0, PricePath, ParamPathE, Params, T, n_d_PT, n_a_PT, n_z, N_j, Names_i, d_grid_PT, a_grid_PT,z_grid_PT, transpathoptions, simoptionsF);
for ii=1:N_i
    nm=Names_i{ii};
    simopt_soloF=simopt_solo{ii};
    simopt_soloF.fastOLG=0;
    FnsToEvaluate_solo=struct('assets',FnsToEvaluate.assets.(nm),'earnings',FnsToEvaluate.earnings.(nm));
    REF=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluate_solo, AgentDistPath_PT0.(nm), PolicyPath_PT0.(nm), PricePath, ParamPath_solo{ii}, Params, T, n_d_PT.(nm), n_a_PT.(nm), n_z, N_j, d_grid_PT.(nm), a_grid_PT.(nm),z_grid_PT.(nm), transpathoptions, simopt_soloF);
    REFL=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate_solo, AgentDistPath_PT0.(nm), PolicyPath_PT0.(nm), PricePath, ParamPath_solo{ii}, Params, T, n_d_PT.(nm), n_a_PT.(nm), n_z, N_j, d_grid_PT.(nm), a_grid_PT.(nm), z_grid_PT.(nm), transpathoptions, simopt_soloF);
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
    fprintf('ExpAsset/GI by ptype: (X3) EvalFnOnTransPath_AllStats_PType vs without permanent types (type %s), every stat, this should be zero: %.3e \n',nm,maxdiff)
    fprintf('ExpAsset/GI by ptype: (X3) EvalFnOnTransPath_AllStats_PType vs without permanent types (type %s), positions NaN (or text) in one but not the other, or missing, this should be zero: %i \n',nm,nanmismatch)
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
    fprintf('ExpAsset/GI by ptype: (X3) LifeCycleProfiles_TransPath_PType vs without permanent types (type %s), all seven outputs, this should be zero: %.3e \n',nm,maxdiff)
    fprintf('ExpAsset/GI by ptype: (X3) LifeCycleProfiles_TransPath_PType vs without permanent types (type %s), entries NaN in one but not the other, this should be zero: %i \n',nm,nanmismatch)
end
clear AllStats_PT LCP_PT PolicyPath_PT0 AgentDistPath_PT0

%% (X4) general eqm: an interest rate and a general eqm condition by ptype, which depends on earnings
% The agent dist in the general eqm path is fastOLG=1 (its default): the experience asset's a2prime
% probabilities are a grid interpolation layer, which the general eqm path does not allow with simoptions.fastOLG=0
GeneralEqmEqns.CapitalMarket=@(r,assets,earnings) r-0.02-0.01*assets/(1+earnings);
PricePath_PT=struct();
PricePath_PT.r=repmat(PricePath.r,N_i,1); % N_i-by-T, a price that depends on ptype
simoptionsG=simoptions;
simoptionsG.fastOLG=1;
for fastOLG=[1,0]
    transpathoptionsGE=struct();
    transpathoptionsGE.fastOLG=fastOLG;
    transpathoptionsGE.GEnewprice=3; % the shooting algorithm
    transpathoptionsGE.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
    transpathoptionsGE.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
    transpathoptionsGE.verbose=0;
    transpathoptionsGE_PT=transpathoptionsGE;
    transpathoptionsGE_PT.GEptype={'CapitalMarket'};
    [PricePath_GE,GEcondnPath_GE]=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPathE, T, V_final_PT, AgentDist_initial_PT, jequaloneDist_PT, n_d_PT, n_a_PT, n_z, N_j, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsGE_PT, simoptionsG, vfoptions);
    for ii=1:N_i
        nm=Names_i{ii};
        simopt_soloG=simopt_solo{ii};
        simopt_soloG.fastOLG=1;
        FnsToEvaluate_solo=struct('assets',FnsToEvaluate.assets.(nm),'earnings',FnsToEvaluate.earnings.(nm));
        [PricePath_solo,GEcondnPath_solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath_solo{ii}, T, V_final_PT.(nm), AgentDist_initial_PT.(nm), jequaloneDist_PT.(nm), n_d_PT.(nm), n_a_PT.(nm), n_z, N_j, d_grid_PT.(nm),a_grid_PT.(nm),z_grid_PT.(nm), pi_z_PT.(nm), ReturnFn_PT.(nm), FnsToEvaluate_solo, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsGE, simopt_soloG, vfopt_solo{ii});
        fprintf('ExpAsset/GI by ptype (fastOLG=%i): (X4) GE by ptype, shooting at maxiter=2, PType vs solo, PricePath r (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(PricePath_GE.r.(nm)(:)-PricePath_solo.r(:))))
        fprintf('ExpAsset/GI by ptype (fastOLG=%i): (X4) GE by ptype, shooting at maxiter=2, PType vs solo, GEcondnPath (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(GEcondnPath_GE.CapitalMarket.(nm)(:)-GEcondnPath_solo.CapitalMarket(:))))
    end
    if fastOLG==1
        % Anderson: not separable by ptype, so against shooting. memory=0 is shooting exactly; the default memory is not.
        transpathoptionsS3=transpathoptionsGE_PT;
        transpathoptionsS3.maxiter=3;
        PricePathS3=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPathE, T, V_final_PT, AgentDist_initial_PT, jequaloneDist_PT, n_d_PT, n_a_PT, n_z, N_j, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsS3, simoptionsG, vfoptions);
        transpathoptionsA=rmfield(transpathoptionsS3,'GEnewprice3');
        transpathoptionsA.GEnewprice=2;
        transpathoptionsA.GEnewprice2.howtoupdate=transpathoptionsS3.GEnewprice3.howtoupdate;
        transpathoptionsA.anderson.warmup=0;
        transpathoptionsA.anderson.safeguard=0;
        transpathoptionsA.anderson.memory=0;
        PricePathA0=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPathE, T, V_final_PT, AgentDist_initial_PT, jequaloneDist_PT, n_d_PT, n_a_PT, n_z, N_j, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptionsG, vfoptions);
        transpathoptionsA=rmfield(transpathoptionsA,'anderson');
        transpathoptionsA.anderson.warmup=0;
        transpathoptionsA.anderson.safeguard=0;
        PricePathA5=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPathE, T, V_final_PT, AgentDist_initial_PT, jequaloneDist_PT, n_d_PT, n_a_PT, n_z, N_j, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptionsG, vfoptions);
        maxdiff0=0; maxdiff5=0;
        for ii=1:N_i
            nm=Names_i{ii};
            maxdiff0=max(maxdiff0,max(abs(PricePathA0.r.(nm)(:)-PricePathS3.r.(nm)(:))));
            maxdiff5=max(maxdiff5,max(abs(PricePathA5.r.(nm)(:)-PricePathS3.r.(nm)(:))));
        end
        fprintf('ExpAsset/GI by ptype: (X4) Anderson with memory=0 vs shooting, both at maxiter=3, all three types, this should be zero: %.3e \n',maxdiff0)
        fprintf('ExpAsset/GI by ptype: (X4) Anderson with the default memory vs shooting, both at maxiter=3, all three types, this should NOT be zero: %.3e \n',maxdiff5)
    end
end

%% (X5) AgentDistOnTransPath_Case1_FHorz with the grid interpolation layer: simoptions.fastOLG=0 vs fastOLG=1
% The fastOLG=0 agent dist with the grid interpolation layer (and with it and an experience asset together) has
% one branch for each of no z/no e, z, e, z and e; the fastOLG=1 one is checked elsewhere (CoreFHorzTPathTests,
% CoreFHorzTPathExpAssetTests). Each model: the stationary dist at Params is the initial dist, then the value fn
% and agent dist paths with transpathoptions.fastOLG=simoptions.fastOLG=0 and =1, which must agree.
n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;
modelnames={'no z, no e','z','e','z and e','no z, no e, experience asset','z, experience asset','e, experience asset','z and e, experience asset'};
for mm=1:8
    usez=any(mm==[2,4,6,8]);
    usee=any(mm==[3,4,7,8]);
    useexpasset=(mm>=5);
    vfoptions5=struct();
    vfoptions5.gridinterplayer=1;
    vfoptions5.ngridinterp=5;
    if usee
        vfoptions5.n_e=n_e;
        vfoptions5.e_grid=e_grid;
        vfoptions5.pi_e=pi_e;
    end
    if usez
        n_z5=n_z; z_grid5=z_grid; pi_z5=pi_z;
    else
        n_z5=0; z_grid5=[]; pi_z5=[];
    end
    if useexpasset
        n_d5=n_d2; d_grid5=d2_grid;
        n_a5=[n_a1,n_a2]; a_grid5=[a1_grid;a2_grid];
        vfoptions5.experienceasset=1;
        vfoptions5.aprimeFn=vfoptions.aprimeFn.typeA;
    else
        n_d5=0; d_grid5=[];
        n_a5=n_a1; a_grid5=a1_grid;
    end
    simoptions5=vfoptions5;
    if useexpasset
        simoptions5.d_grid=d_grid5;
        simoptions5.a_grid=a_grid5;
    end
    switch mm
        case 1
            ReturnFn5=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
        case 2
            ReturnFn5=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
        case 3
            ReturnFn5=@(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
        case 4
            ReturnFn5=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
        case 5
            ReturnFn5=@(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension);
        case 6
            ReturnFn5=@(d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_noe_nosemiz(d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension);
        case 7
            ReturnFn5=@(d2,a1prime,a1,a2,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_e_nosemiz(d2,a1prime,a1,a2,e,r,w,kappa_j,sigma,agej,Jr,pension);
        case 8
            ReturnFn5=@(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_e_nosemiz(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
    end
    % zero assets, mid points for any shocks
    jsize=n_a5;
    jidx=num2cell(ones(1,length(n_a5)));
    if usez
        jsize=[jsize,n_z];
        jidx{end+1}=ceil(n_z/2);
    end
    if usee
        jsize=[jsize,n_e];
        jidx{end+1}=ceil(n_e/2);
    end
    if isscalar(jsize)
        jsize=[jsize,1];
    end
    jequaloneDist5=zeros(jsize,'gpuArray');
    jequaloneDist5(jidx{:})=1;

    [V_final5,Policy_final5]=ValueFnIter_Case1_FHorz(n_d5,n_a5,n_z5,N_j,d_grid5,a_grid5,z_grid5,pi_z5,ReturnFn5,Params,DiscountFactorParamNames,[],vfoptions5);
    AgentDist_initial5=StationaryDist_FHorz_Case1(jequaloneDist5,AgeWeightParamNames,Policy_final5,n_d5,n_a5,n_z5,N_j,pi_z5,Params,simoptions5);
    PolicyPath5=cell(2,1); AgentDistPath5=cell(2,1);
    for fastOLG=[0,1]
        transpathoptions5=struct();
        transpathoptions5.fastOLG=fastOLG;
        simoptions5.fastOLG=fastOLG;
        [~,PolicyPath5{fastOLG+1}]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final5, Policy_final5, Params, n_d5, n_a5, n_z5, N_j, d_grid5, a_grid5,z_grid5, pi_z5, DiscountFactorParamNames, ReturnFn5, transpathoptions5, vfoptions5);
        AgentDistPath5{fastOLG+1}=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial5, jequaloneDist5, PricePath, ParamPath, PolicyPath5{fastOLG+1}, AgeWeightParamNames,n_d5,n_a5,n_z5,N_j,pi_z5, T,Params, transpathoptions5, simoptions5);
    end
    % The L2 index is the second last row of the policy (the last is the L2 flag)
    L2row=PolicyPath5{1}(end-1,:);
    fprintf('ExpAsset/GI by ptype: (X5) grid interpolation, %s, fastOLG=0 vs fastOLG=1, PolicyPath, this should be zero: %.3e \n',modelnames{mm},gather(max(abs(PolicyPath5{1}(:)-PolicyPath5{2}(:)))))
    fprintf('ExpAsset/GI by ptype: (X5) grid interpolation, %s, fastOLG=0 vs fastOLG=1, AgentDistPath, this should be zero: %.3e \n',modelnames{mm},gather(max(abs(AgentDistPath5{1}(:)-AgentDistPath5{2}(:)))))
    fprintf('ExpAsset/GI by ptype: (X5) grid interpolation, %s, for scale, the share of policies strictly between two grid points, this should NOT be zero: %.3e \n',modelnames{mm},gather(mean(L2row>1 & L2row<vfoptions5.ngridinterp+2)))
    clear PolicyPath5 AgentDistPath5 V_final5 Policy_final5 AgentDist_initial5
end


output=struct();

end
