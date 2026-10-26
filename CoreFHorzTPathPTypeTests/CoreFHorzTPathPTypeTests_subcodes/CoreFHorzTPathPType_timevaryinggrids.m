function output=CoreFHorzTPathPType_timevaryinggrids(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline)
% z and e grids (and pi_z, pi_e) that vary over the transition path, differently for each permanent type.
% Two types, both on the model with z and e of CoreFHorzTPath_nod_z_e_nosemiz:
%   typeA: input grids that vary over the path: z_grid.typeA [n_z,T] and pi_z.typeA [n_z,n_z,T], and
%          e_grid.typeA, pi_e.typeA [n_e,T] in vfoptions/simoptions; the grids drift up, by 1+0.1*t/T
%   typeB: ExogShockFn and EiidShockFn (by ptype) whose scale and mix are on the ParamPath (a structure by
%          ptype); the grids drift down, by 1-0.1*t/T
% and in both pi_z and pi_e move part of the way towards uniform over the path.
% Each type must be exactly what the commands without permanent types give with that type's inputs
% (those commands are checked against the stationary commands period by period in CoreFHorzTPathTests
% doPart(8)).
%   (G1) ValueFnOnTransPath and AgentDistOnTransPath (PType) vs solo, by ptype, fastOLG=0 and fastOLG=1
%   (G2) the PType value fn and agent dist paths, fastOLG=1 vs fastOLG=0
%   (G3) EvalFnOnTransPath_AllStats and LifeCycleProfiles_TransPath (PType) vs solo, by ptype
%   (G4) the same PType stats, pooled and by ptype, vs the stationary PType commands period by period, given
%        that period's grids of each ptype directly (not via the shock fns or the transition-path shock setup)
%   (G5) general eqm path with an interest rate and a general eqm condition by ptype (GEptype), which depends
%        on earnings: shooting at maxiter=2 vs solo, by ptype, fastOLG=1 and 0; Anderson with memory=0 vs
%        shooting at maxiter=3 (and with the default memory it must differ)
%   and, for scale, that the two types' earnings differ (must not be zero, or the checks cannot tell the
%   types' grids apart)

OutputTypes={'Mean','Median','Variance','LorenzCurve','Gini','QuantileCutoffs','QuantileMeans'};
scalarstat=[1,1,1,0,1,0,0]; % Mean, Median, Variance, Gini are one number per age

Names_i={'typeA','typeB'};
N_i=2;
PTypeDistParamNames={'ptweights_tvg'};
Params.ptweights_tvg=[0.55;0.45];
n_d=0;
d_grid=[];

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z,e) a;
FnsToEvaluate.earnings=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
fnnames=fieldnames(FnsToEvaluate);
jequaloneDist=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist(1,ceil(n_z/2),ceil(n_e/2))=1;

%% The grids of each type over the path
driftA=1+0.1*(1:1:T)/T; % [1,T]
driftB=1-0.1*(1:1:T)/T;
pimixA=0.3*(1:1:T)/T;
pimixB=0.2*(1:1:T)/T;
z_gridT_A=z_grid.*driftA; % [n_z,T]
pi_zT_A=(1-reshape(pimixA,[1,1,T])).*pi_z+reshape(pimixA,[1,1,T]).*ones(n_z,n_z)/n_z; % [n_z,n_z,T]
e_gridT_A=e_grid.*driftA; % [n_e,T]
pi_eT_A=(1-pimixA).*pi_e+pimixA.*ones(n_e,1)/n_e; % [n_e,T]

% PType inputs: typeA by input grids, typeB by shock fns
z_grid_PT=struct('typeA',z_gridT_A,'typeB',z_grid);
pi_z_PT=struct('typeA',pi_zT_A,'typeB',pi_z);
vfoptions=struct();
vfoptions.n_e=n_e;
vfoptions.e_grid=struct('typeA',e_gridT_A,'typeB',e_grid);
vfoptions.pi_e=struct('typeA',pi_eT_A,'typeB',pi_e);
vfoptions.ExogShockFn.typeB=@(zscale,zpimix) deal(z_grid*zscale,(1-zpimix)*pi_z+zpimix*ones(n_z,n_z)/n_z);
vfoptions.EiidShockFn.typeB=@(escale,epimix) deal(e_grid*escale,(1-epimix)*pi_e+epimix*ones(n_e,1)/n_e);
simoptions=vfoptions;
Params.zscale=1; Params.zpimix=0; Params.escale=1; Params.epimix=0;
ParamPathG=ParamPath;
ParamPathG.zscale=struct('typeA',ones(1,T),'typeB',driftB); % typeA does not use them
ParamPathG.zpimix=struct('typeA',zeros(1,T),'typeB',pimixB);
ParamPathG.escale=struct('typeA',ones(1,T),'typeB',driftB);
ParamPathG.epimix=struct('typeA',zeros(1,T),'typeB',pimixB);

% Solo inputs, without permanent types
vfopt_solo=cell(N_i,1); ParamPath_solo=cell(N_i,1); z_grid_solo=cell(N_i,1); pi_z_solo=cell(N_i,1);
vfopt_solo{1}=struct('n_e',n_e);
vfopt_solo{1}.e_grid=e_gridT_A;
vfopt_solo{1}.pi_e=pi_eT_A;
ParamPath_solo{1}=ParamPath;
z_grid_solo{1}=z_gridT_A;
pi_z_solo{1}=pi_zT_A;
vfopt_solo{2}=struct('n_e',n_e);
vfopt_solo{2}.e_grid=e_grid;
vfopt_solo{2}.pi_e=pi_e;
vfopt_solo{2}.ExogShockFn=vfoptions.ExogShockFn.typeB;
vfopt_solo{2}.EiidShockFn=vfoptions.EiidShockFn.typeB;
ParamPath_solo{2}=ParamPath;
ParamPath_solo{2}.zscale=driftB;
ParamPath_solo{2}.zpimix=pimixB;
ParamPath_solo{2}.escale=driftB;
ParamPath_solo{2}.epimix=pimixB;
z_grid_solo{2}=z_grid;
pi_z_solo{2}=pi_z;

%% Initial agent dist and V_final: the stationary ones with the plain grids (the same for both types)
vfoptions0=struct('n_e',n_e);
vfoptions0.e_grid=e_grid;
vfoptions0.pi_e=pi_e;
[V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions0);
AgentDist_init=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params,vfoptions0);
V_final_PT=struct('typeA',V_final,'typeB',V_final);
Policy_final_PT=struct('typeA',Policy_final,'typeB',Policy_final);
AgentDist_initial_PT=struct('typeA',AgentDist_init,'typeB',AgentDist_init);
AgentDist_initial_PT.ptweights=Params.ptweights_tvg;

%% (G1)/(G2) the value fn and agent dist paths
VPath_PT=cell(2,1); PolicyPath_PT=cell(2,1); AgentDistPath_PT=cell(2,1);
for fastOLG=[0,1]
    transpathoptions=struct();
    transpathoptions.fastOLG=fastOLG;
    [VPath_PT{fastOLG+1},PolicyPath_PT{fastOLG+1}]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPathG, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid_PT, pi_z_PT, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    AgentDistPath_PT{fastOLG+1}=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath, ParamPathG, PolicyPath_PT{fastOLG+1}, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z_PT, T,Params, transpathoptions, simoptions);
    for ii=1:N_i
        nm=Names_i{ii};
        [VPath_solo,PolicyPath_solo]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath_solo{ii}, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid_solo{ii}, pi_z_solo{ii}, DiscountFactorParamNames, ReturnFn, transpathoptions, vfopt_solo{ii});
        AgentDistPath_solo=AgentDistOnTransPath_Case1_FHorz(AgentDist_init, jequaloneDist, PricePath, ParamPath_solo{ii}, PolicyPath_solo, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z_solo{ii}, T,Params, transpathoptions, vfopt_solo{ii});
        fprintf('Time-varying grids by ptype (fastOLG=%i): (G1) PType vs solo (type %s), VPath, this should be zero: %.3e \n',fastOLG,nm,max(abs(VPath_PT{fastOLG+1}.(nm)(:)-VPath_solo(:))))
        fprintf('Time-varying grids by ptype (fastOLG=%i): (G1) PType vs solo (type %s), PolicyPath, this should be zero: %.3e \n',fastOLG,nm,max(abs(PolicyPath_PT{fastOLG+1}.(nm)(:)-PolicyPath_solo(:))))
        fprintf('Time-varying grids by ptype (fastOLG=%i): (G1) PType vs solo (type %s), AgentDistPath, this should be zero: %.3e \n',fastOLG,nm,max(abs(AgentDistPath_PT{fastOLG+1}.(nm)(:)-AgentDistPath_solo(:))))
    end
end
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('Time-varying grids by ptype: (G2) fastOLG=1 vs fastOLG=0 (type %s), VPath, this should be zero: %.3e \n',nm,max(abs(VPath_PT{2}.(nm)(:)-VPath_PT{1}.(nm)(:))))
    fprintf('Time-varying grids by ptype: (G2) fastOLG=1 vs fastOLG=0 (type %s), AgentDistPath, this should be zero: %.3e \n',nm,max(abs(AgentDistPath_PT{2}.(nm)(:)-AgentDistPath_PT{1}.(nm)(:))))
end

%% (G3) the stats, by ptype, vs solo
transpathoptions=struct();
transpathoptions.fastOLG=0;
AllStats_PT=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_PT{1}, PolicyPath_PT{1}, PricePath, ParamPathG, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid_PT, transpathoptions, simoptions);
LCP_PT=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate, AgentDistPath_PT{1}, PolicyPath_PT{1}, PricePath, ParamPathG, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid_PT, transpathoptions, simoptions);
for ii=1:N_i
    nm=Names_i{ii};
    REF=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluate, AgentDistPath_PT{1}.(nm), PolicyPath_PT{1}.(nm), PricePath, ParamPath_solo{ii}, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid_solo{ii}, transpathoptions, vfopt_solo{ii});
    REFL=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate, AgentDistPath_PT{1}.(nm), PolicyPath_PT{1}.(nm), PricePath, ParamPath_solo{ii}, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid_solo{ii}, transpathoptions, vfopt_solo{ii});
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
    fprintf('Time-varying grids by ptype: (G3) EvalFnOnTransPath_AllStats_PType vs without permanent types (type %s), every stat, this should be zero: %.3e \n',nm,maxdiff)
    fprintf('Time-varying grids by ptype: (G3) EvalFnOnTransPath_AllStats_PType vs without permanent types (type %s), positions NaN (or text) in one but not the other, or missing, this should be zero: %i \n',nm,nanmismatch)
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
    fprintf('Time-varying grids by ptype: (G3) LifeCycleProfiles_TransPath_PType vs without permanent types (type %s), all seven outputs, this should be zero: %.3e \n',nm,maxdiff)
    fprintf('Time-varying grids by ptype: (G3) LifeCycleProfiles_TransPath_PType vs without permanent types (type %s), entries NaN in one but not the other, this should be zero: %i \n',nm,nanmismatch)
end
fprintf('Time-varying grids by ptype: for scale, mean earnings of typeA minus typeB, max over the path, this should NOT be zero: %.3e \n',max(abs(AllStats_PT.earnings.typeA.Mean-AllStats_PT.earnings.typeB.Mean)))

%% (G4) the stats vs the stationary PType commands period by period, given each period's grids directly
PricePathNamesS=fieldnames(PricePath);
ParamPathNamesS=fieldnames(ParamPath);
pi_eT_B=(1-pimixB).*pi_e+pimixB.*ones(n_e,1)/n_e; % [n_e,T]
AScell=cell(1,T);
LCPcell=cell(1,T);
for tt=1:T
    Params_tt=Params;
    for kk=1:length(PricePathNamesS)
        Params_tt.(PricePathNamesS{kk})=PricePath.(PricePathNamesS{kk})(tt);
    end
    for kk=1:length(ParamPathNamesS)
        Params_tt.(ParamPathNamesS{kk})=ParamPath.(ParamPathNamesS{kk})(tt);
    end
    StationaryDist_tt=struct(); Policy_tt=struct();
    for ii=1:N_i
        StationaryDist_tt.(Names_i{ii})=AgentDistPath_PT{1}.(Names_i{ii})(:,:,:,:,tt);
        Policy_tt.(Names_i{ii})=PolicyPath_PT{1}.(Names_i{ii})(:,:,:,:,:,tt);
    end
    StationaryDist_tt.ptweights=Params.ptweights_tvg;
    z_grid_tt=struct('typeA',z_gridT_A(:,tt),'typeB',z_grid*driftB(tt));
    simoptions_tt=struct('n_e',n_e);
    simoptions_tt.e_grid=struct('typeA',e_gridT_A(:,tt),'typeB',e_grid*driftB(tt));
    simoptions_tt.pi_e=struct('typeA',pi_eT_A(:,tt),'typeB',pi_eT_B(:,tt));
    AScell{tt}=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_tt,Policy_tt,FnsToEvaluate,Params_tt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid_tt,simoptions_tt);
    LCPcell{tt}=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_tt,Policy_tt,FnsToEvaluate,Params_tt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid_tt,simoptions_tt);
end
maxdiff=0; nanmismatch=0;
for tt=1:T
    AS=AScell{tt};
    queue={{}};
    while ~isempty(queue)
        fieldpath=queue{1};
        queue(1)=[];
        if isempty(fieldpath)
            O=AS;
        else
            O=getfield(AS,fieldpath{:});
        end
        if isstruct(O)
            Onames=fieldnames(O);
            for kk=1:length(Onames)
                queue{end+1}=[fieldpath,Onames(kk)];
            end
        elseif ischar(O) || isstring(O) || iscell(O) % text (the comments): must be the same text
            P=getfield(AllStats_PT,fieldpath{:});
            nanmismatch=nanmismatch+(~isequal(P{tt},O));
        else
            O=gather(O);
            P=getfield(AllStats_PT,fieldpath{:});
            P=P(:,tt);
            maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
            nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
        end
    end
end
fprintf('Time-varying grids by ptype: (G4) EvalFnOnTransPath_AllStats_PType vs the stationary PType AllStats with that period''s grids, pooled and by ptype, every stat, this should be zero: %.3e \n',maxdiff)
fprintf('Time-varying grids by ptype: (G4) EvalFnOnTransPath_AllStats_PType vs the stationary PType AllStats with that period''s grids, positions NaN (or text) in one but not the other, this should be zero: %i \n',nanmismatch)
groupnames=[{''},Names_i];
maxdiff=0; nanmismatch=0;
for gg=1:length(groupnames)
    for ff=1:length(fnnames)
        if isempty(groupnames{gg})
            OUTff=LCP_PT.(fnnames{ff});
        else
            OUTff=LCP_PT.(fnnames{ff}).(groupnames{gg});
        end
        for oo=1:length(OutputTypes)
            PI=gather(OUTff.initialages.(OutputTypes{oo}));
            PB=gather(OUTff.bornduringtranstion.(OutputTypes{oo}));
            OI=nan(size(PI));
            OB=nan(size(PB));
            % Initial cohorts: aged jj in period 1, so aged jj+tt-1 in period tt, alive while that is at most N_j
            for tt=1:min(T,N_j)
                if isempty(groupnames{gg})
                    temp=gather(LCPcell{tt}.(fnnames{ff}).(OutputTypes{oo}));
                else
                    temp=gather(LCPcell{tt}.(fnnames{ff}).(groupnames{gg}).(OutputTypes{oo}));
                end
                for jj=1:N_j-tt+1
                    if scalarstat(oo)==1
                        OI(jj,jj+tt-1)=temp(jj+tt-1);
                    else
                        OI(:,jj,jj+tt-1)=temp(:,jj+tt-1);
                    end
                end
            end
            % Cohorts born during the transition: aged 1 in period bb, so aged age in period bb+age-1
            for bb=1:T
                for age=1:N_j
                    tt=bb+age-1;
                    if tt<=T
                        if isempty(groupnames{gg})
                            temp=gather(LCPcell{tt}.(fnnames{ff}).(OutputTypes{oo}));
                        else
                            temp=gather(LCPcell{tt}.(fnnames{ff}).(groupnames{gg}).(OutputTypes{oo}));
                        end
                        if scalarstat(oo)==1
                            OB(bb,age)=temp(age);
                        else
                            OB(:,bb,age)=temp(:,age);
                        end
                    end
                end
            end
            maxdiff=max([maxdiff; abs(PI(:)-OI(:)); abs(PB(:)-OB(:))],[],'omitnan');
            nanmismatch=nanmismatch+sum(isnan(PI(:))~=isnan(OI(:)))+sum(isnan(PB(:))~=isnan(OB(:)));
        end
    end
end
fprintf('Time-varying grids by ptype: (G4) LifeCycleProfiles_TransPath_PType vs the per-period stationary PType profiles with that period''s grids, pooled and by ptype, all seven outputs, this should be zero: %.3e \n',maxdiff)
fprintf('Time-varying grids by ptype: (G4) LifeCycleProfiles_TransPath_PType vs the per-period stationary PType profiles, entries NaN in one but not the other (the alive pattern), this should be zero: %i \n',nanmismatch)
clear AScell LCPcell

%% (G5) general eqm: an interest rate and a general eqm condition by ptype, which depends on earnings
GeneralEqmEqns.CapitalMarket=@(r,assets,earnings) r-0.02-0.01*assets/(1+earnings);
PricePath_PT=struct();
PricePath_PT.r=repmat(PricePath.r,N_i,1); % N_i-by-T, a price that depends on ptype
for fastOLG=[1,0]
    transpathoptionsGE=struct();
    transpathoptionsGE.fastOLG=fastOLG;
    transpathoptionsGE.GEnewprice=3; % the shooting algorithm
    transpathoptionsGE.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
    transpathoptionsGE.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
    transpathoptionsGE.verbose=0;
    transpathoptionsGE_PT=transpathoptionsGE;
    transpathoptionsGE_PT.GEptype={'CapitalMarket'};
    [PricePath_GE,GEcondnPath_GE]=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPathG, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid_PT, pi_z_PT, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsGE_PT, simoptions, vfoptions);
    for ii=1:N_i
        nm=Names_i{ii};
        [PricePath_solo,GEcondnPath_solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath_solo{ii}, T, V_final, AgentDist_init, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid_solo{ii}, pi_z_solo{ii}, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsGE, vfopt_solo{ii}, vfopt_solo{ii});
        fprintf('Time-varying grids by ptype (fastOLG=%i): (G5) GE by ptype, shooting at maxiter=2, PType vs solo, PricePath r (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(PricePath_GE.r.(nm)(:)-PricePath_solo.r(:))))
        fprintf('Time-varying grids by ptype (fastOLG=%i): (G5) GE by ptype, shooting at maxiter=2, PType vs solo, GEcondnPath (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(GEcondnPath_GE.CapitalMarket.(nm)(:)-GEcondnPath_solo.CapitalMarket(:))))
    end
    if fastOLG==1
        % Anderson: not separable by ptype, so against shooting. memory=0 is shooting exactly; the default memory is not.
        transpathoptionsS3=transpathoptionsGE_PT;
        transpathoptionsS3.maxiter=3;
        PricePathS3=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPathG, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid_PT, pi_z_PT, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsS3, simoptions, vfoptions);
        transpathoptionsA=rmfield(transpathoptionsS3,'GEnewprice3');
        transpathoptionsA.GEnewprice=2;
        transpathoptionsA.GEnewprice2.howtoupdate=transpathoptionsS3.GEnewprice3.howtoupdate;
        transpathoptionsA.anderson.warmup=0;
        transpathoptionsA.anderson.safeguard=0;
        transpathoptionsA.anderson.memory=0;
        PricePathA0=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPathG, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid_PT, pi_z_PT, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
        transpathoptionsA=rmfield(transpathoptionsA,'anderson');
        transpathoptionsA.anderson.warmup=0;
        transpathoptionsA.anderson.safeguard=0;
        PricePathA5=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPathG, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid_PT, pi_z_PT, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
        fprintf('Time-varying grids by ptype: (G5) Anderson with memory=0 vs shooting, both at maxiter=3, both types, this should be zero: %.3e \n',max([max(abs(PricePathA0.r.typeA(:)-PricePathS3.r.typeA(:))),max(abs(PricePathA0.r.typeB(:)-PricePathS3.r.typeB(:)))]))
        fprintf('Time-varying grids by ptype: (G5) Anderson with the default memory vs shooting, both at maxiter=3, both types, this should NOT be zero: %.3e \n',max([max(abs(PricePathA5.r.typeA(:)-PricePathS3.r.typeA(:))),max(abs(PricePathA5.r.typeB(:)-PricePathS3.r.typeB(:)))]))
    end
end

output=struct();

end
