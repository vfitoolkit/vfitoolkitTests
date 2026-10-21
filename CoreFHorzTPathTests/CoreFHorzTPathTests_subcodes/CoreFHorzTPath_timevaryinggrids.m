function output=CoreFHorzTPath_timevaryinggrids(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,n_d2_semiz,d2_grid_semiz,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline)
% z and e grids that vary over the transition path: z_grid as [n_z,T] and simoptions.e_grid as [n_e,T]
% (the 'time-varying, age-independent, stacked' input shapes of ExogShockSetup_FHorz_TPath).
%
% THE METHOD. EvalFnOnTransPath_AllStats_Case1_FHorz and LifeCycleProfiles_TransPath_FHorz_Case1 are,
% period by period, the stationary command applied to that period's slice of the agent dist and policy
% paths, at that period's prices, parameters AND grids. So the oracle is the stationary command called on
% each period by hand, given that period's column of the grid (z_gridT(:,tt), e_gridT(:,tt)) directly,
% not via the transition-path shock setup the commands use.
%
% Models: (1) z, (2) e, (3) z and e (no d), (4) semiz and z, (5) semiz, z and e (binary d2 for semiz). Each
% grid drifts up over the path, by a factor 1+0.1*t/T, and in the drift run pi_z and pi_e also move over the
% path (a third of the way towards uniform by period T), so the per-period transition matrices are exercised
% with real variation over t.
% Three runs per model, each the full chain ValueFnOnTransPath -> AgentDistOnTransPath -> stats:
%   drift:    the drifting grids (and moving pi_z, pi_e)
%   constant: grids that are time-varying in form but constant over t (the plain grid repeated T times)
%   plain:    the ordinary time-invariant grids
% The value fn and agent dist paths in V1-V4 use transpathoptions.fastOLG=0; V5 checks fastOLG=1 against them.
%   (V1) drift: AllStats along the path vs the stationary AllStats period by period, with that period's grids
%   (V2) drift: LifeCycleProfiles along the path vs the stationary profiles period by period (the cohort
%        bookkeeping, as in figs 1-16), and NaN exactly where no cohort is
%   (V3) drift: EvalFnOnTransPath_AggVars Mean vs the AllStats Mean
%   (V4) constant vs plain: PolicyPath, AgentDistPath, AllStats and LifeCycleProfiles all identical (the
%        time-varying route through value fn, agent dist and stats agrees with the ordinary one)
%   and, for scale, how far the drift moves earnings from the plain run (must not be zero, or V1-V3 are empty)
%   (V5) drift and constant: the value fn and agent dist paths with fastOLG=1 equal those with fastOLG=0
%   (V6) one iteration of TransitionPath_Case1_FHorz with the drifting grids, fastOLG=1 vs fastOLG=0, with a
%        general eqm condition that depends on earnings; for z and e also with only one of them drifting
%   (V7) pi_e that depends on age (not on the path): value fn and agent dist paths, fastOLG=1 vs fastOLG=0
%   (V8) models 1-3: ExogShockFn/EiidShockFn returning the scaled grid and mixed pi, with the scale and mix in
%        ParamPath, reproduce the drift run (value fn, agent dist, AllStats, LifeCycleProfiles; fastOLG=0 and 1;
%        one GE iteration), and with them held fixed (not in ParamPath) reproduce the plain run
%   (V10) transpathoptions.fastOLG=0 with simoptions.fastOLG=1 (the default in TransitionPath_Case1_FHorz):
%        the agent dist path, and one GE iteration, against both fastOLG=0

OutputTypes={'Mean','Median','Variance','LorenzCurve','Gini','QuantileCutoffs','QuantileMeans'};
scalarstat=[1,1,1,0,1,0,0]; % Mean, Median, Variance, Gini are one number per age

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;
driftfactor=1+0.1*(1:1:T)/T; % [1,T]
z_gridT=z_grid.*driftfactor; % [n_z,T]
e_gridT=e_grid.*driftfactor; % [n_e,T]
pimixT=0.3*(1:1:T)/T; % [1,T]
pi_zT=(1-reshape(pimixT,[1,1,T])).*pi_z+reshape(pimixT,[1,1,T]).*ones(n_z,n_z)/n_z; % [n_z,n_z,T]
pi_eT=(1-pimixT).*pi_e+pimixT.*ones(n_e,1)/n_e; % [n_e,T]

PricePathNamesV=fieldnames(PricePath);
ParamPathNamesV=fieldnames(ParamPath);

for mm=1:5
    %% The model
    % vfoptionsM, simoptionsM: the model's own options (the semi-exogenous state), on which every run's options are built
    vfoptionsM=struct();
    simoptionsM=struct();
    n_d=0;
    d_grid=[];
    if mm>=4
        n_d=n_d2_semiz;
        d_grid=d2_grid_semiz;
        vfoptionsM.n_semiz=vfoptionsbaseline.n_semiz;
        vfoptionsM.semiz_grid=vfoptionsbaseline.semiz_grid;
        vfoptionsM.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
        simoptionsM.n_semiz=simoptionsbaseline.n_semiz;
        simoptionsM.semiz_grid=simoptionsbaseline.semiz_grid;
        simoptionsM.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
        simoptionsM.d_grid=d_grid;
        n_semiz=vfoptionsbaseline.n_semiz;
    end
    if mm==1
        modelname='z';
        n_zm=n_z; pi_zm=pi_z;
        usee=0;
        jequaloneDist=zeros([n_a,n_z],'gpuArray');
        jequaloneDist(1,ceil(n_z/2))=1;
        ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
    elseif mm==2
        modelname='e';
        n_zm=0; pi_zm=[];
        usee=1;
        jequaloneDist=zeros([n_a,n_e],'gpuArray');
        jequaloneDist(1,ceil(n_e/2))=1;
        ReturnFn=@(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(aprime,a,e) a;
        FnsToEvaluate.earnings=@(aprime,a,e,w,kappa_j) w*kappa_j*e;
    elseif mm==3
        modelname='z and e';
        n_zm=n_z; pi_zm=pi_z;
        usee=1;
        jequaloneDist=zeros([n_a,n_z,n_e],'gpuArray');
        jequaloneDist(1,ceil(n_z/2),ceil(n_e/2))=1;
        ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(aprime,a,z,e) a;
        FnsToEvaluate.earnings=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
    elseif mm==4
        modelname='semiz and z';
        n_zm=n_z; pi_zm=pi_z;
        usee=0;
        jequaloneDist=zeros([n_a,n_semiz,n_z],'gpuArray');
        jequaloneDist(1,ceil(n_semiz/2),ceil(n_z/2))=1;
        ReturnFn=@(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,z) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,w,kappa_j) w*kappa_j*z*semiz;
    else
        modelname='semiz, z and e';
        n_zm=n_z; pi_zm=pi_z;
        usee=1;
        jequaloneDist=zeros([n_a,n_semiz,n_z,n_e],'gpuArray');
        jequaloneDist(1,ceil(n_semiz/2),ceil(n_z/2),ceil(n_e/2))=1;
        ReturnFn=@(d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,z,e) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,e,w,kappa_j) w*kappa_j*z*e*semiz;
    end
    fnnames=fieldnames(FnsToEvaluate);

    % The initial dist: the stationary dist with the plain grids (the same for all three runs)
    vfoptions0=vfoptionsM;
    simoptions0=simoptionsM;
    if usee==1
        vfoptions0.n_e=n_e; vfoptions0.e_grid=e_grid; vfoptions0.pi_e=pi_e;
        simoptions0.n_e=n_e; simoptions0.e_grid=e_grid; simoptions0.pi_e=pi_e;
    end
    if n_zm>0
        z_grid0=z_grid;
    else
        z_grid0=[];
    end
    [~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zm,N_j,d_grid,a_grid,z_grid0,pi_zm,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions0);
    AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_zm,N_j,pi_zm,Params,simoptions0);

    %% The three runs: drift, constant, plain
    VPathR=cell(1,3);
    PolicyPathR=cell(1,3);
    AgentDistPathR=cell(1,3);
    AllStatsPathR=cell(1,3);
    LCPPathR=cell(1,3);
    for cc=1:3
        % The grids handed to the path commands, and the final-period grids (for V_final)
        if cc==1
            z_gridX=z_gridT; e_gridX=e_gridT;
            z_gridF=z_gridT(:,T); e_gridF=e_gridT(:,T);
        elseif cc==2
            z_gridX=repmat(z_grid,1,T); e_gridX=repmat(e_grid,1,T);
            z_gridF=z_grid; e_gridF=e_grid;
        else
            z_gridX=z_grid; e_gridX=e_grid;
            z_gridF=z_grid; e_gridF=e_grid;
        end
        if n_zm==0
            z_gridX=[]; z_gridF=[];
        end
        % The transition probabilities handed to the path commands, and the final-period ones (for V_final): they
        % move over the path in the drift run only
        if cc==1
            pi_zX=pi_zT; pi_eX=pi_eT;
            pi_zF=pi_zT(:,:,T); pi_eF=pi_eT(:,T);
        else
            pi_zX=pi_z; pi_eX=pi_e;
            pi_zF=pi_z; pi_eF=pi_e;
        end
        if n_zm==0
            pi_zX=[]; pi_zF=[];
        end
        vfoptionsF=vfoptionsM;
        vfoptionsX=vfoptionsM;
        simoptionsX=simoptionsM;
        simoptionsX.fastOLG=0;
        if usee==1
            vfoptionsF.n_e=n_e; vfoptionsF.e_grid=e_gridF; vfoptionsF.pi_e=pi_eF;
            vfoptionsX.n_e=n_e; vfoptionsX.e_grid=e_gridX; vfoptionsX.pi_e=pi_eX;
            simoptionsX.n_e=n_e; simoptionsX.e_grid=e_gridX; simoptionsX.pi_e=pi_eX;
        end
        transpathoptions=transpathoptionsbaseline;
        transpathoptions.fastOLG=0;

        [V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zm,N_j,d_grid,a_grid,z_gridF,pi_zF,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsF);
        if cc==1 % kept for (V8)
            V_finalD=V_final; Policy_finalD=Policy_final;
        elseif cc==3
            V_finalP=V_final; Policy_finalP=Policy_final;
        end
        [VPathR{cc},PolicyPathR{cc}]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_gridX, pi_zX, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionsX);
        AgentDistPathR{cc}=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath, PolicyPathR{cc}, AgeWeightParamNames,n_d,n_a,n_zm,N_j,pi_zX, T,Params, transpathoptions, simoptionsX);
        % (V10) the agent dist path with transpathoptions.fastOLG=0 but simoptions.fastOLG=1
        if cc==1
            simoptionsS1=simoptionsX;
            simoptionsS1.fastOLG=1;
            AgentDistPathS1=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath, PolicyPathR{cc}, AgeWeightParamNames,n_d,n_a,n_zm,N_j,pi_zX, T,Params, transpathoptions, simoptionsS1);
            fprintf('Shocks over the path, %s: (V10) drifting grids, AgentDistPath with transpathoptions.fastOLG=0 and simoptions.fastOLG=1 vs both 0, this should be zero: %.3e \n',modelname,gather(max(abs(AgentDistPathS1(:)-AgentDistPathR{cc}(:)))))
            clear AgentDistPathS1
        end
        % (V5) the same value fn and agent dist paths with fastOLG=1, for the two runs whose grids are time-varying in form
        if cc<=2
            transpathoptionsF1=transpathoptions;
            transpathoptionsF1.fastOLG=1;
            simoptionsF1=simoptionsX;
            simoptionsF1.fastOLG=1;
            [VPathF1,PolicyPathF1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_gridX, pi_zX, DiscountFactorParamNames, ReturnFn, transpathoptionsF1, vfoptionsX);
            AgentDistPathF1=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath, PolicyPathF1, AgeWeightParamNames,n_d,n_a,n_zm,N_j,pi_zX, T,Params, transpathoptionsF1, simoptionsF1);
            if cc==1
                runname='drifting grids';
            else
                runname='grids time-varying in form but constant';
            end
            fprintf('Shocks over the path, %s: (V5) %s, fastOLG=1 vs fastOLG=0, VPath, this should be zero: %.3e \n',modelname,runname,gather(max(abs(VPathF1(:)-VPathR{cc}(:)))))
            fprintf('Shocks over the path, %s: (V5) %s, fastOLG=1 vs fastOLG=0, PolicyPath, this should be zero: %.3e \n',modelname,runname,gather(max(abs(PolicyPathF1(:)-PolicyPathR{cc}(:)))))
            fprintf('Shocks over the path, %s: (V5) %s, fastOLG=1 vs fastOLG=0, AgentDistPath, this should be zero: %.3e \n',modelname,runname,gather(max(abs(AgentDistPathF1(:)-AgentDistPathR{cc}(:)))))
            clear VPathF1 PolicyPathF1 AgentDistPathF1
        end
        AllStatsPathR{cc}=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluate, AgentDistPathR{cc}, PolicyPathR{cc}, PricePath, ParamPath, Params, T, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_gridX, transpathoptions, simoptionsX);
        LCPPathR{cc}=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate, AgentDistPathR{cc}, PolicyPathR{cc}, PricePath, ParamPath, Params, T, n_d, n_a, n_zm, N_j, d_grid, a_grid, z_gridX, transpathoptions, simoptionsX);
        if cc==1
            AggVarsPathD=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPathR{cc}, PolicyPathR{cc}, PricePath, ParamPath, Params, T, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_gridX, transpathoptions, simoptionsX);
        end
    end

    %% The oracle for the drift run: the stationary commands on each period's slice, at that period's prices, parameters and grids
    AgentDistPath=AgentDistPathR{1};
    PolicyPath=PolicyPathR{1};
    AllStatsPath=AllStatsPathR{1};
    LCPPath=LCPPathR{1};
    sliceAgentDist=repmat({':'},1,ndims(AgentDistPath)-1);
    slicePolicy=repmat({':'},1,ndims(PolicyPath)-1);
    AScell=cell(1,T);
    LCPcell=cell(1,T);
    for tt=1:T
        Params_tt=Params;
        for kk=1:length(PricePathNamesV)
            Params_tt.(PricePathNamesV{kk})=PricePath.(PricePathNamesV{kk})(tt);
        end
        for kk=1:length(ParamPathNamesV)
            Params_tt.(ParamPathNamesV{kk})=ParamPath.(ParamPathNamesV{kk})(tt);
        end
        simoptions_tt=simoptionsM;
        if usee==1
            simoptions_tt.n_e=n_e; simoptions_tt.e_grid=e_gridT(:,tt); simoptions_tt.pi_e=pi_e;
        end
        if n_zm>0
            z_grid_tt=z_gridT(:,tt);
        else
            z_grid_tt=[];
        end
        AScell{tt}=EvalFnOnAgentDist_AllStats_FHorz_Case1(AgentDistPath(sliceAgentDist{:},tt),PolicyPath(slicePolicy{:},tt),FnsToEvaluate,Params_tt,[],n_d,n_a,n_zm,N_j,d_grid,a_grid,z_grid_tt,simoptions_tt);
        LCPcell{tt}=LifeCycleProfiles_FHorz_Case1(AgentDistPath(sliceAgentDist{:},tt),PolicyPath(slicePolicy{:},tt),FnsToEvaluate,Params_tt,[],n_d,n_a,n_zm,N_j,d_grid,a_grid,z_grid_tt,simoptions_tt);
    end

    % (V1) AllStats: every stat of every FnsToEvaluate, every period
    maxdiff=0; nanmismatch=0;
    for tt=1:T
        AS=AScell{tt};
        % every field of AllStats, to any depth (MoreInequality is a struct; the comments are text)
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
                P=getfield(AllStatsPath,fieldpath{:});
                nanmismatch=nanmismatch+(~isequal(P{tt},O));
            else
                O=gather(O);
                P=getfield(AllStatsPath,fieldpath{:});
                P=P(:,tt);
                maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
                nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
            end
        end
    end
    fprintf('Shocks over the path, %s: (V1) EvalFnOnTransPath_AllStats vs the stationary AllStats period by period with that period''s grids, every stat, this should be zero: %.3e \n',modelname,maxdiff)
    fprintf('Shocks over the path, %s: (V1) EvalFnOnTransPath_AllStats vs the stationary AllStats period by period with that period''s grids, positions NaN (or text) in one but not the other, this should be zero: %i \n',modelname,nanmismatch)

    % (V2) Life-cycle profiles: the cohort bookkeeping
    maxdiffI=0; maxdiffB=0; nanmisI=0; nanmisB=0;
    for ff=1:length(fnnames)
        for oo=1:length(OutputTypes)
            PI=gather(LCPPath.(fnnames{ff}).initialages.(OutputTypes{oo}));
            PB=gather(LCPPath.(fnnames{ff}).bornduringtranstion.(OutputTypes{oo}));
            OI=nan(size(PI));
            OB=nan(size(PB));
            % Initial cohorts: aged jj in period 1, so aged jj+tt-1 in period tt, alive while that is at most N_j
            for tt=1:min(T,N_j)
                temp=gather(LCPcell{tt}.(fnnames{ff}).(OutputTypes{oo}));
                for jj=1:N_j-tt+1
                    if scalarstat(oo)==1
                        OI(jj,jj+tt-1)=temp(jj+tt-1);
                    else
                        OI(:,jj,jj+tt-1)=temp(:,jj+tt-1);
                    end
                end
            end
            % Cohorts born during the transition: aged 1 in period ii, so aged age in period ii+age-1
            for ii=1:T
                for age=1:N_j
                    tt=ii+age-1;
                    if tt<=T
                        temp=gather(LCPcell{tt}.(fnnames{ff}).(OutputTypes{oo}));
                        if scalarstat(oo)==1
                            OB(ii,age)=temp(age);
                        else
                            OB(:,ii,age)=temp(:,age);
                        end
                    end
                end
            end
            maxdiffI=max([maxdiffI; abs(PI(:)-OI(:))],[],'omitnan');
            maxdiffB=max([maxdiffB; abs(PB(:)-OB(:))],[],'omitnan');
            nanmisI=nanmisI+sum(isnan(PI(:))~=isnan(OI(:)));
            nanmisB=nanmisB+sum(isnan(PB(:))~=isnan(OB(:)));
        end
    end
    fprintf('Shocks over the path, %s: (V2) LifeCycleProfiles_TransPath initialages vs the per-period stationary profiles with that period''s grids, all seven outputs, this should be zero: %.3e \n',modelname,maxdiffI)
    fprintf('Shocks over the path, %s: (V2) LifeCycleProfiles_TransPath bornduringtranstion vs the per-period stationary profiles with that period''s grids, all seven outputs, this should be zero: %.3e \n',modelname,maxdiffB)
    fprintf('Shocks over the path, %s: (V2) LifeCycleProfiles_TransPath entries NaN in one but not the other (the alive pattern), initialages %i and bornduringtranstion %i, this should be zero: %i \n',modelname,nanmisI,nanmisB,nanmisI+nanmisB)

    % (V3) the AggVars command against the AllStats Mean
    maxdiff=0;
    for ff=1:length(fnnames)
        maxdiff=max(maxdiff,max(abs(AllStatsPath.(fnnames{ff}).Mean-gather(AggVarsPathD.(fnnames{ff}).Mean))));
    end
    fprintf('Shocks over the path, %s: (V3) EvalFnOnTransPath_AggVars Mean vs EvalFnOnTransPath_AllStats Mean, with the drifting grids, this should be zero: %.3e \n',modelname,maxdiff)

    % (V4) constant-over-t time-varying grids vs the plain grids: the whole chain must agree
    fprintf('Shocks over the path, %s: (V4) grids time-varying in form but constant vs plain grids, PolicyPath, this should be zero: %.3e \n',modelname,gather(max(abs(PolicyPathR{2}(:)-PolicyPathR{3}(:)))))
    fprintf('Shocks over the path, %s: (V4) grids time-varying in form but constant vs plain grids, AgentDistPath, this should be zero: %.3e \n',modelname,gather(max(abs(AgentDistPathR{2}(:)-AgentDistPathR{3}(:)))))
    for ss=1:2
        if ss==1
            SC=AllStatsPathR{2}; SP=AllStatsPathR{3}; statname='EvalFnOnTransPath_AllStats';
        else
            SC=LCPPathR{2}; SP=LCPPathR{3}; statname='LifeCycleProfiles_TransPath';
        end
        maxdiff=0; nanmismatch=0;
        % every field, to any depth, all periods at once
        queue={{}};
        while ~isempty(queue)
            fieldpath=queue{1};
            queue(1)=[];
            if isempty(fieldpath)
                O=SP; P=SC;
            else
                O=getfield(SP,fieldpath{:}); P=getfield(SC,fieldpath{:});
            end
            if isstruct(O)
                Onames=fieldnames(O);
                for kk=1:length(Onames)
                    queue{end+1}=[fieldpath,Onames(kk)];
                end
            elseif iscell(O)
                nanmismatch=nanmismatch+(~isequal(P,O));
            else
                O=gather(O); P=gather(P);
                maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
                nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
            end
        end
        fprintf('Shocks over the path, %s: (V4) grids time-varying in form but constant vs plain grids, %s, every stat, this should be zero: %.3e \n',modelname,statname,maxdiff)
        fprintf('Shocks over the path, %s: (V4) grids time-varying in form but constant vs plain grids, %s, positions NaN (or text) in one but not the other, this should be zero: %i \n',modelname,statname,nanmismatch)
    end

    % For scale: the drift must actually move things
    fprintf('Shocks over the path, %s: for scale, the drifting grids move mean earnings from the plain run by up to %.3e (must not be zero, or V1-V3 are empty) \n',modelname,max(abs(AllStatsPathR{1}.earnings.Mean-AllStatsPathR{3}.earnings.Mean)))

    % (V6) One iteration of the GE transition path with the drifting grids, fastOLG=1 vs fastOLG=0 (as in the figs,
    % a shape check). The general eqm condition depends on earnings, so the price update after one iteration
    % depends on the grids. With z and e: both drifting, only z drifting, only e drifting (when only one of them
    % varies the other has no _T version, which is its own branch of the shock setup).
    if n_zm>0 && usee==1
        nGEcases=3;
    else
        nGEcases=1;
    end
    GeneralEqmEqns=struct();
    GeneralEqmEqns.earningsgap=@(earnings) earnings-0.5;
    for gg=1:nGEcases
        if gg==1
            z_gridG=z_gridT; e_gridG=e_gridT; GEcasename='drifting grids';
        elseif gg==2
            z_gridG=z_gridT; e_gridG=e_grid; GEcasename='only the z grid drifting';
        else
            z_gridG=z_grid; e_gridG=e_gridT; GEcasename='only the e grid drifting';
        end
        if n_zm==0
            z_gridG=[]; z_gridGF=[];
        else
            z_gridGF=z_gridG(:,end);
        end
        vfoptionsG=vfoptionsM;
        vfoptionsGF=vfoptionsM;
        simoptionsG=simoptionsM;
        if usee==1
            vfoptionsG.n_e=n_e; vfoptionsG.e_grid=e_gridG; vfoptionsG.pi_e=pi_e;
            vfoptionsGF.n_e=n_e; vfoptionsGF.e_grid=e_gridG(:,end); vfoptionsGF.pi_e=pi_e;
            simoptionsG.n_e=n_e; simoptionsG.e_grid=e_gridG; simoptionsG.pi_e=pi_e;
        end
        V_finalG=ValueFnIter_Case1_FHorz(n_d,n_a,n_zm,N_j,d_grid,a_grid,z_gridGF,pi_zm,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsGF);
        transpathoptionsG=transpathoptionsbaseline;
        transpathoptionsG.maxiter=1;
        transpathoptionsG.GEnewprice=2; % Anderson acceleration
        transpathoptionsG.GEnewprice2.howtoupdate={'earningsgap','r',0,0.1};
        transpathoptionsG.fastOLG=1;
        simoptionsG.fastOLG=1;
        PricePathG1=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_finalG, AgentDist_initial, jequaloneDist, n_d, n_a, n_zm, N_j, d_grid,a_grid,z_gridG, pi_zm, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsG, simoptionsG, vfoptionsG);
        transpathoptionsG.fastOLG=0;
        simoptionsG.fastOLG=0;
        PricePathG0=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_finalG, AgentDist_initial, jequaloneDist, n_d, n_a, n_zm, N_j, d_grid,a_grid,z_gridG, pi_zm, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsG, simoptionsG, vfoptionsG);
        fprintf('Shocks over the path, %s: (V6) one iteration of the GE path, %s, fastOLG=1 vs fastOLG=0, this should be zero: %.3e \n',modelname,GEcasename,max(abs(PricePathG1.r-PricePathG0.r)))
        fprintf('Shocks over the path, %s: for scale, that one iteration moves r by up to %.3e (must not be zero, or (V6) is empty) \n',modelname,max(abs(PricePathG1.r-PricePath.r)))
        % (V10) transpathoptions.fastOLG=0 with simoptions.fastOLG left to its default (=1 in TransitionPath_Case1_FHorz)
        simoptionsG=rmfield(simoptionsG,'fastOLG');
        PricePathGS=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_finalG, AgentDist_initial, jequaloneDist, n_d, n_a, n_zm, N_j, d_grid,a_grid,z_gridG, pi_zm, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsG, simoptionsG, vfoptionsG);
        fprintf('Shocks over the path, %s: (V10) one iteration of the GE path, %s, transpathoptions.fastOLG=0 with the default simoptions.fastOLG vs both 0, this should be zero: %.3e \n',modelname,GEcasename,max(abs(PricePathGS.r-PricePathG0.r)))
        if gg==1 % kept for (V8)
            V_finalGD=V_finalG;
            PricePathGD1=PricePathG1;
        end
    end

    % (V7) pi_e that depends on age (but not on the path), value fn and agent dist paths, fastOLG=1 vs fastOLG=0.
    % The fastOLG value fn wants pi_e as (a-j,e) with a fastest; the setup built it with j fastest, which only
    % matters when pi_e depends on age.
    if usee==1 && mm<=3
        pi_eJ=zeros(n_e,N_j);
        for jj=1:N_j
            pi_eJ(:,jj)=(1-0.5*(jj-1)/(N_j-1))*pi_e+0.5*(jj-1)/(N_j-1)*ones(n_e,1)/n_e; % moves from pi_e at age 1 to halfway to uniform at age N_j
        end
        vfoptionsA=struct();
        vfoptionsA.n_e=n_e; vfoptionsA.e_grid=e_grid; vfoptionsA.pi_e=pi_eJ;
        simoptionsA=struct();
        simoptionsA.n_e=n_e; simoptionsA.e_grid=e_grid; simoptionsA.pi_e=pi_eJ;
        [V_finalA,Policy_finalA]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zm,N_j,d_grid,a_grid,z_grid0,pi_zm,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsA);
        transpathoptionsA=transpathoptionsbaseline;
        transpathoptionsA.fastOLG=0;
        simoptionsA.fastOLG=0;
        [VPathA0,PolicyPathA0]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_finalA, Policy_finalA, Params, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_grid0, pi_zm, DiscountFactorParamNames, ReturnFn, transpathoptionsA, vfoptionsA);
        AgentDistPathA0=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath, PolicyPathA0, AgeWeightParamNames,n_d,n_a,n_zm,N_j,pi_zm, T,Params, transpathoptionsA, simoptionsA);
        transpathoptionsA.fastOLG=1;
        simoptionsA.fastOLG=1;
        [VPathA1,PolicyPathA1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_finalA, Policy_finalA, Params, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_grid0, pi_zm, DiscountFactorParamNames, ReturnFn, transpathoptionsA, vfoptionsA);
        AgentDistPathA1=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath, PolicyPathA1, AgeWeightParamNames,n_d,n_a,n_zm,N_j,pi_zm, T,Params, transpathoptionsA, simoptionsA);
        fprintf('Shocks over the path, %s: (V7) pi_e depending on age, fastOLG=1 vs fastOLG=0, VPath, this should be zero: %.3e \n',modelname,gather(max(abs(VPathA1(:)-VPathA0(:)))))
        fprintf('Shocks over the path, %s: (V7) pi_e depending on age, fastOLG=1 vs fastOLG=0, PolicyPath, this should be zero: %.3e \n',modelname,gather(max(abs(PolicyPathA1(:)-PolicyPathA0(:)))))
        fprintf('Shocks over the path, %s: (V7) pi_e depending on age, fastOLG=1 vs fastOLG=0, AgentDistPath, this should be zero: %.3e \n',modelname,gather(max(abs(AgentDistPathA1(:)-AgentDistPathA0(:)))))
        fprintf('Shocks over the path, %s: for scale, making pi_e depend on age moves VPath by up to %.3e (must not be zero, or (V7) is empty) \n',modelname,gather(max(abs(VPathA0(:)-VPathR{3}(:)))))
        clear VPathA0 VPathA1 PolicyPathA0 PolicyPathA1 AgentDistPathA0 AgentDistPathA1
    end

    % (V8) ExogShockFn/EiidShockFn instead of the input grids. The shock fns return the plain grid times a scale,
    % and pi mixed towards uniform by a weight; with the scale and weight in ParamPath (driftfactor, pimixT) they
    % are exactly the drift run, and held fixed at 1 and 0 (not in ParamPath) exactly the plain run.
    if mm<=3
        vfoptionsS=vfoptionsM;
        simoptionsS=simoptionsM;
        ParamPathS=ParamPath;
        ParamsS=Params;
        if n_zm>0
            vfoptionsS.ExogShockFn=@(zscale,zpimix) deal(z_grid*zscale,(1-zpimix)*pi_z+zpimix*ones(n_z,n_z)/n_z);
            simoptionsS.ExogShockFn=vfoptionsS.ExogShockFn;
            ParamPathS.zscale=driftfactor;
            ParamPathS.zpimix=pimixT;
            ParamsS.zscale=1;
            ParamsS.zpimix=0;
        end
        if usee==1
            vfoptionsS.n_e=n_e; vfoptionsS.e_grid=e_grid; vfoptionsS.pi_e=pi_e;
            vfoptionsS.EiidShockFn=@(escale,epimix) deal(e_grid*escale,(1-epimix)*pi_e+epimix*ones(n_e,1)/n_e);
            simoptionsS.n_e=n_e; simoptionsS.e_grid=e_grid; simoptionsS.pi_e=pi_e;
            simoptionsS.EiidShockFn=vfoptionsS.EiidShockFn;
            ParamPathS.escale=driftfactor;
            ParamPathS.epimix=pimixT;
            ParamsS.escale=1;
            ParamsS.epimix=0;
        end
        transpathoptionsS=transpathoptionsbaseline;
        transpathoptionsS.fastOLG=0;
        simoptionsS.fastOLG=0;
        % The shock fns depending on ParamPath, vs the drift run
        [VPathS,PolicyPathS]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPathS, T, V_finalD, Policy_finalD, ParamsS, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_grid0, pi_zm, DiscountFactorParamNames, ReturnFn, transpathoptionsS, vfoptionsS);
        AgentDistPathS=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPathS, PolicyPathS, AgeWeightParamNames,n_d,n_a,n_zm,N_j,pi_zm, T,ParamsS, transpathoptionsS, simoptionsS);
        AllStatsPathS=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluate, AgentDistPathS, PolicyPathS, PricePath, ParamPathS, ParamsS, T, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_grid0, transpathoptionsS, simoptionsS);
        LCPPathS=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate, AgentDistPathS, PolicyPathS, PricePath, ParamPathS, ParamsS, T, n_d, n_a, n_zm, N_j, d_grid, a_grid, z_grid0, transpathoptionsS, simoptionsS);
        fprintf('Shocks over the path, %s: (V8) shock fns depending on ParamPath vs the drifting input grids, VPath, this should be zero: %.3e \n',modelname,gather(max(abs(VPathS(:)-VPathR{1}(:)))))
        fprintf('Shocks over the path, %s: (V8) shock fns depending on ParamPath vs the drifting input grids, PolicyPath, this should be zero: %.3e \n',modelname,gather(max(abs(PolicyPathS(:)-PolicyPathR{1}(:)))))
        fprintf('Shocks over the path, %s: (V8) shock fns depending on ParamPath vs the drifting input grids, AgentDistPath, this should be zero: %.3e \n',modelname,gather(max(abs(AgentDistPathS(:)-AgentDistPathR{1}(:)))))
        for ss=1:2
            if ss==1
                SC=AllStatsPathS; SP=AllStatsPathR{1}; statname='EvalFnOnTransPath_AllStats';
            else
                SC=LCPPathS; SP=LCPPathR{1}; statname='LifeCycleProfiles_TransPath';
            end
            maxdiff=0; nanmismatch=0;
            % every field, to any depth, all periods at once
            queue={{}};
            while ~isempty(queue)
                fieldpath=queue{1};
                queue(1)=[];
                if isempty(fieldpath)
                    O=SP; P=SC;
                else
                    O=getfield(SP,fieldpath{:}); P=getfield(SC,fieldpath{:});
                end
                if isstruct(O)
                    Onames=fieldnames(O);
                    for kk=1:length(Onames)
                        queue{end+1}=[fieldpath,Onames(kk)];
                    end
                elseif iscell(O)
                    nanmismatch=nanmismatch+(~isequal(P,O));
                else
                    O=gather(O); P=gather(P);
                    maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
                    nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
                end
            end
            fprintf('Shocks over the path, %s: (V8) shock fns depending on ParamPath vs the drifting input grids, %s, every stat, this should be zero: %.3e \n',modelname,statname,maxdiff)
            fprintf('Shocks over the path, %s: (V8) shock fns depending on ParamPath vs the drifting input grids, %s, positions NaN (or text) in one but not the other, this should be zero: %i \n',modelname,statname,nanmismatch)
        end
        % fastOLG=1, against the same (fastOLG=0) drift run
        transpathoptionsS.fastOLG=1;
        simoptionsS.fastOLG=1;
        [VPathS1,PolicyPathS1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPathS, T, V_finalD, Policy_finalD, ParamsS, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_grid0, pi_zm, DiscountFactorParamNames, ReturnFn, transpathoptionsS, vfoptionsS);
        AgentDistPathS1=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPathS, PolicyPathS1, AgeWeightParamNames,n_d,n_a,n_zm,N_j,pi_zm, T,ParamsS, transpathoptionsS, simoptionsS);
        fprintf('Shocks over the path, %s: (V8) shock fns depending on ParamPath, fastOLG=1, vs the drifting input grids, VPath, this should be zero: %.3e \n',modelname,gather(max(abs(VPathS1(:)-VPathR{1}(:)))))
        fprintf('Shocks over the path, %s: (V8) shock fns depending on ParamPath, fastOLG=1, vs the drifting input grids, PolicyPath, this should be zero: %.3e \n',modelname,gather(max(abs(PolicyPathS1(:)-PolicyPathR{1}(:)))))
        fprintf('Shocks over the path, %s: (V8) shock fns depending on ParamPath, fastOLG=1, vs the drifting input grids, AgentDistPath, this should be zero: %.3e \n',modelname,gather(max(abs(AgentDistPathS1(:)-AgentDistPathR{1}(:)))))
        clear VPathS PolicyPathS AgentDistPathS AllStatsPathS LCPPathS VPathS1 PolicyPathS1 AgentDistPathS1
        % One GE iteration, against the (V6) drifting-grid one; that used the plain pi, so the mixing weight is zero here
        ParamPathSG=ParamPathS;
        if n_zm>0
            ParamPathSG.zpimix=zeros(1,T);
        end
        if usee==1
            ParamPathSG.epimix=zeros(1,T);
        end
        transpathoptionsG1=transpathoptionsG; % maxiter=1, Anderson, as in (V6)
        transpathoptionsG1.fastOLG=1;
        PricePathSG1=TransitionPath_Case1_FHorz(PricePath, ParamPathSG, T, V_finalGD, AgentDist_initial, jequaloneDist, n_d, n_a, n_zm, N_j, d_grid,a_grid,z_grid0, pi_zm, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsS, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsG1, simoptionsS, vfoptionsS);
        fprintf('Shocks over the path, %s: (V8) one iteration of the GE path, fastOLG=1, shock fns depending on ParamPath vs the drifting input grids, this should be zero: %.3e \n',modelname,max(abs(PricePathSG1.r-PricePathGD1.r)))
        % The shock fns held fixed (scale 1, weight 0, not in ParamPath), vs the plain run
        transpathoptionsS.fastOLG=0;
        simoptionsS.fastOLG=0;
        [VPathS,PolicyPathS]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_finalP, Policy_finalP, ParamsS, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_grid0, pi_zm, DiscountFactorParamNames, ReturnFn, transpathoptionsS, vfoptionsS);
        AgentDistPathS=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath, PolicyPathS, AgeWeightParamNames,n_d,n_a,n_zm,N_j,pi_zm, T,ParamsS, transpathoptionsS, simoptionsS);
        fprintf('Shocks over the path, %s: (V8) shock fns not depending on the path vs the plain grids, VPath, this should be zero: %.3e \n',modelname,gather(max(abs(VPathS(:)-VPathR{3}(:)))))
        fprintf('Shocks over the path, %s: (V8) shock fns not depending on the path vs the plain grids, PolicyPath, this should be zero: %.3e \n',modelname,gather(max(abs(PolicyPathS(:)-PolicyPathR{3}(:)))))
        fprintf('Shocks over the path, %s: (V8) shock fns not depending on the path vs the plain grids, AgentDistPath, this should be zero: %.3e \n',modelname,gather(max(abs(AgentDistPathS(:)-AgentDistPathR{3}(:)))))
        clear VPathS PolicyPathS AgentDistPathS
    end

    clear VPathR PolicyPathR AgentDistPathR AllStatsPathR LCPPathR AgentDistPath PolicyPath AllStatsPath LCPPath AScell LCPcell AggVarsPathD
end

output=struct();

end
