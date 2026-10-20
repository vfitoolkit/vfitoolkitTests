function output=CoreFHorzTPath_stats(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,n_d2_semiz,d2_grid_semiz,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline)
% Statistics along an FHorz transition path: LifeCycleProfiles_TransPath_FHorz_Case1 and
% EvalFnOnTransPath_AllStats_Case1_FHorz.
%
% THE METHOD. Both commands are, period by period, the stationary command applied to that period's
% slice of the agent dist and policy paths, with the Parameters at that period's prices and parameters.
% So the oracle is exactly that: LifeCycleProfiles_FHorz_Case1 and EvalFnOnAgentDist_AllStats_FHorz_Case1
% called on each period by hand. Every number the path commands return must match, and must sit where
% it belongs: for the life-cycle profiles that is the cohort bookkeeping,
%   initialages(jj,age):        the cohort aged jj in period 1, now aged age, so period age-jj+1
%   bornduringtranstion(ii,age): the cohort aged 1 in period ii, now aged age, so period ii+age-1
% and an entry that does not exist (cohort not alive) must be NaN while every one that does must not
% be, which is checked as a count of positions that are NaN in one and not the other.
%
% Models: (1) the life-cycle model without d (all checks); (2) without z, (3) with z and e, and
% (4) with a semi-exogenous state (AllStats only: LifeCycleProfiles_TransPath_FHorz_Case1 does not
% handle semiz), which repeat (L1),(L2),(A1) since each has its own shape handling.
%   (L1) life-cycle profiles, all seven outputs, vs the per-period stationary profiles
%   (L2) the alive pattern: NaN exactly where no cohort is
%   (L3) on a constant (do-nothing) path from the stationary dist, every period's profile is the stationary one
%   (A1) AllStats, every stat, every period, vs the per-period stationary AllStats
%   (A2) the AllStats Mean vs EvalFnOnTransPath_AggVars_Case1_FHorz
%   (A3) on the do-nothing path, every period's AllStats is the stationary one
%   (A4) A1 with simoptions.conditionalrestrictions (restricted sample mass and restricted stats)
%   (A5) _tminus1 and _tplus1 in the FnsToEvaluate vs the same functions on hand-shifted paths

OutputTypes={'Mean','Median','Variance','LorenzCurve','Gini','QuantileCutoffs','QuantileMeans'};
scalarstat=[1,1,1,0,1,0,0]; % Mean, Median, Variance, Gini are one number per age

for mm=1:4
    %% The model
    vfoptions=struct();
    simoptions=struct();
    if mm==1
        modelname='no d';
        n_d=0; d_grid=[]; n_zm=n_z; z_gridm=z_grid; pi_zm=pi_z;
        jequaloneDist=zeros([n_a,n_z],'gpuArray');
        jequaloneDist(1,ceil(n_z/2))=1;
        ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
    elseif mm==2
        modelname='no z';
        n_d=0; d_grid=[]; n_zm=0; z_gridm=[]; pi_zm=[];
        jequaloneDist=zeros(n_a,1,'gpuArray');
        jequaloneDist(1)=1;
        ReturnFn=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(aprime,a) a;
        FnsToEvaluate.earnings=@(aprime,a,w,kappa_j) w*kappa_j;
    elseif mm==3
        modelname='z and e';
        n_d=0; d_grid=[]; n_zm=n_z; z_gridm=z_grid; pi_zm=pi_z;
        vfoptions.n_e=vfoptionsbaseline.n_e;
        vfoptions.e_grid=vfoptionsbaseline.e_grid;
        vfoptions.pi_e=vfoptionsbaseline.pi_e;
        simoptions.n_e=simoptionsbaseline.n_e;
        simoptions.e_grid=simoptionsbaseline.e_grid;
        simoptions.pi_e=simoptionsbaseline.pi_e;
        jequaloneDist=zeros([n_a,n_z,vfoptions.n_e],'gpuArray');
        jequaloneDist(1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;
        ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(aprime,a,z,e) a;
        FnsToEvaluate.earnings=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
    else
        modelname='semiz';
        n_d=n_d2_semiz; d_grid=d2_grid_semiz; n_zm=n_z; z_gridm=z_grid; pi_zm=pi_z;
        vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
        vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
        vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
        simoptions.n_semiz=simoptionsbaseline.n_semiz;
        simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
        simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
        simoptions.d_grid=d_grid;
        jequaloneDist=zeros([n_a,vfoptions.n_semiz,n_z],'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_semiz/2),ceil(n_z/2))=1;
        ReturnFn=@(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,z) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,w,kappa_j) w*kappa_j*z*semiz;
    end
    doLCP=(mm<=3);

    [V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zm,N_j,d_grid,a_grid,z_gridm,pi_zm,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_final,n_d,n_a,n_zm,N_j,pi_zm,Params,simoptions);

    %% Two paths: the transition (PricePath, ParamPath), and a do-nothing one from the stationary dist (model 1 only)
    if mm==1
        ncases=2;
    else
        ncases=1;
    end
    for cc=1:ncases
        if cc==1
            PricePathX=PricePath;
            ParamPathX=ParamPath;
            casename='the path';
        else
            PricePathX.r=Params.r*ones(1,T);
            ParamPathX.w=Params.w*ones(1,T);
            casename='do-nothing path';
        end
        transpathoptions=transpathoptionsbaseline;
        [~,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePathX, ParamPathX, T, V_final, Policy_final, Params, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_gridm, pi_zm, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
        AgentDistPath=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePathX, ParamPathX, PolicyPath, AgeWeightParamNames,n_d,n_a,n_zm,N_j,pi_zm, T,Params, transpathoptions, simoptions);

        AllStatsPath=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluate, AgentDistPath, PolicyPath, PricePathX, ParamPathX, Params, T, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_gridm, transpathoptions, simoptions);
        if doLCP
            LCPPath=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate, AgentDistPath, PolicyPath, PricePathX, ParamPathX, Params, T, n_d, n_a, n_zm, N_j, d_grid, a_grid, z_gridm, transpathoptions, simoptions);
        end

        %% The oracle, period by period
        % On the transition: the stationary commands on each period's slice, at that period's prices and
        % parameters. On the do-nothing path: the stationary commands on the stationary dist itself.
        sliceAgentDist=repmat({':'},1,ndims(AgentDistPath)-1);
        slicePolicy=repmat({':'},1,ndims(PolicyPath)-1);
        AScell=cell(1,T);
        LCPcell=cell(1,T);
        if cc==2
            AS_stat=EvalFnOnAgentDist_AllStats_FHorz_Case1(AgentDist_initial,Policy_final,FnsToEvaluate,Params,[],n_d,n_a,n_zm,N_j,d_grid,a_grid,z_gridm,simoptions);
            LCP_stat=LifeCycleProfiles_FHorz_Case1(AgentDist_initial,Policy_final,FnsToEvaluate,Params,[],n_d,n_a,n_zm,N_j,d_grid,a_grid,z_gridm,simoptions);
        end
        for tt=1:T
            if cc==1
                Params_tt=Params;
                Params_tt.r=PricePathX.r(tt);
                Params_tt.w=ParamPathX.w(tt);
                AScell{tt}=EvalFnOnAgentDist_AllStats_FHorz_Case1(AgentDistPath(sliceAgentDist{:},tt),PolicyPath(slicePolicy{:},tt),FnsToEvaluate,Params_tt,[],n_d,n_a,n_zm,N_j,d_grid,a_grid,z_gridm,simoptions);
                if doLCP
                    LCPcell{tt}=LifeCycleProfiles_FHorz_Case1(AgentDistPath(sliceAgentDist{:},tt),PolicyPath(slicePolicy{:},tt),FnsToEvaluate,Params_tt,[],n_d,n_a,n_zm,N_j,d_grid,a_grid,z_gridm,simoptions);
                end
            else
                AScell{tt}=AS_stat;
                LCPcell{tt}=LCP_stat;
            end
        end

        %% (A1)/(A3) AllStats: every stat of every FnsToEvaluate, every period
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
        if cc==1
            fprintf('Stats, %s: (A1) EvalFnOnTransPath_AllStats vs the stationary AllStats period by period, every stat, this should be zero: %.3e \n',modelname,maxdiff)
            fprintf('Stats, %s: (A1) EvalFnOnTransPath_AllStats vs the stationary AllStats period by period, every stat, positions NaN (or text) in one but not the other, this should be zero: %i \n',modelname,nanmismatch)
        else
            fprintf('Stats, %s: (A3) on the %s, EvalFnOnTransPath_AllStats vs the stationary AllStats, every period, this should be zero: %.3e \n',modelname,casename,maxdiff)
            fprintf('Stats, %s: (A3) on the %s, EvalFnOnTransPath_AllStats vs the stationary AllStats, every period, positions NaN (or text) in one but not the other, this should be zero: %i \n',modelname,casename,nanmismatch)
        end

        %% (L1)/(L2)/(L3) Life-cycle profiles: the cohort bookkeeping
        if doLCP
            maxdiffI=0; maxdiffB=0; nanmisI=0; nanmisB=0;
            fnnames=fieldnames(FnsToEvaluate);
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
            if cc==1
                labelL='(L1)'; labelP='(L2)';
            else
                labelL='(L3)'; labelP='(L3)';
            end
            fprintf('Stats, %s: %s LifeCycleProfiles_TransPath initialages vs the %s profiles, all seven outputs, this should be zero: %.3e \n',modelname,labelL,casename,maxdiffI)
            fprintf('Stats, %s: %s LifeCycleProfiles_TransPath bornduringtranstion vs the %s profiles, all seven outputs, this should be zero: %.3e \n',modelname,labelL,casename,maxdiffB)
            fprintf('Stats, %s: %s LifeCycleProfiles_TransPath entries NaN in one but not the other (the alive pattern), initialages %i and bornduringtranstion %i, this should be zero: %i \n',modelname,labelP,nanmisI,nanmisB,nanmisI+nanmisB)
        end

        %% Model 1, the transition only: (A2), (A4), (A5)
        if mm==1 && cc==1
            % (A2) the AllStats Mean against the AggVars command
            AggVarsPath=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath, PolicyPath, PricePathX, ParamPathX, Params, T, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_gridm, transpathoptions, simoptions);
            fprintf('Stats, %s: (A2) EvalFnOnTransPath_AllStats Mean vs EvalFnOnTransPath_AggVars, this should be zero: %.3e \n',modelname,max([max(abs(AllStatsPath.assets.Mean-gather(AggVarsPath.assets.Mean))),max(abs(AllStatsPath.earnings.Mean-gather(AggVarsPath.earnings.Mean)))]))

            % (A4) with conditional restrictions
            simoptionsR=simoptions;
            simoptionsR.conditionalrestrictions.working=@(aprime,a,z,kappa_j) (kappa_j>0);
            simoptionsR.conditionalrestrictions.wealthy=@(aprime,a,z) (a>1);
            AllStatsPathR=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluate, AgentDistPath, PolicyPath, PricePathX, ParamPathX, Params, T, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_gridm, transpathoptions, simoptionsR);
            maxdiff=0; nanmismatch=0; nfields=0;
            for tt=1:T
                Params_tt=Params;
                Params_tt.r=PricePathX.r(tt);
                Params_tt.w=ParamPathX.w(tt);
                AS=EvalFnOnAgentDist_AllStats_FHorz_Case1(AgentDistPath(sliceAgentDist{:},tt),PolicyPath(slicePolicy{:},tt),FnsToEvaluate,Params_tt,[],n_d,n_a,n_zm,N_j,d_grid,a_grid,z_gridm,simoptionsR);
                % every field, to any depth, including the restricted ones
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
                        P=getfield(AllStatsPathR,fieldpath{:});
                        nanmismatch=nanmismatch+(~isequal(P{tt},O)); nfields=nfields+1;
                    else
                        O=gather(O);
                        P=getfield(AllStatsPathR,fieldpath{:});
                        P=P(:,tt);
                        maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
                        nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:))); nfields=nfields+1;
                    end
                end
            end
            fprintf('Stats, %s: (A4) with conditionalrestrictions, %i stat-periods compared (the restricted ones included), this should be zero: %.3e \n',modelname,nfields,maxdiff)
            fprintf('Stats, %s: (A4) with conditionalrestrictions, %i stat-periods compared (the restricted ones included), positions NaN (or text) in one but not the other, this should be zero: %i \n',modelname,nfields,nanmismatch)
            fprintf('Stats, %s: for scale, the restricted sample masses move between %.3f and %.3f (working) and %.3f and %.3f (wealthy) over the path \n',modelname,min(AllStatsPathR.working.RestrictedSampleMass),max(AllStatsPathR.working.RestrictedSampleMass),min(AllStatsPathR.wealthy.RestrictedSampleMass),max(AllStatsPathR.wealthy.RestrictedSampleMass))

            % (A5) _tminus1 and _tplus1 inside the FnsToEvaluate, against the same functions on hand-shifted paths
            initialvalues.r=0.04; % PricePath.r(1) is 0.03
            initialvalues.w=1.2;  % ParamPath.w(1) is 1.1
            initialvalues.assets=0.5;
            transpathoptionsL=transpathoptions;
            transpathoptionsL.initialvalues=initialvalues;
            ParamPathO=ParamPathX;
            ParamPathO.rlag=[initialvalues.r, PricePathX.r(1:T-1)];
            ParamPathO.rlead=[PricePathX.r(2:T), PricePathX.r(T)]; % _tplus1 in period T is period T
            ParamPathO.wlag=[initialvalues.w, ParamPathX.w(1:T-1)];
            ParamPathO.Alag=[initialvalues.assets, AllStatsPath.assets.Mean(1:T-1)];
            ParamsO=Params;
            ParamsO.rlag=ParamPathO.rlag(1); ParamsO.rlead=ParamPathO.rlead(1); ParamsO.wlag=ParamPathO.wlag(1); ParamsO.Alag=ParamPathO.Alag(1);
            FnsToEvaluateL=FnsToEvaluate;
            FnsToEvaluateL.x=@(aprime,a,z,r_tminus1,r_tplus1,w_tminus1) a*r_tminus1+aprime*r_tplus1+w_tminus1*z;
            FnsToEvaluateL.y=@(aprime,a,z,assets_tminus1) a*assets_tminus1;
            FnsToEvaluateO=FnsToEvaluate;
            FnsToEvaluateO.x=@(aprime,a,z,rlag,rlead,wlag) a*rlag+aprime*rlead+wlag*z;
            FnsToEvaluateO.y=@(aprime,a,z,Alag) a*Alag;
            AllStatsPathL=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluateL, AgentDistPath, PolicyPath, PricePathX, ParamPathO, Params, T, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_gridm, transpathoptionsL, simoptions);
            AllStatsPathO=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluateO, AgentDistPath, PolicyPath, PricePathX, ParamPathO, ParamsO, T, n_d, n_a, n_zm, N_j, d_grid, a_grid,z_gridm, transpathoptions, simoptions);
            maxdiff=0; nanmismatch=0;
            lagnames={'x','y'};
            % every field of x and y, to any depth, all periods at once (numeric fields are numel-by-T, text ones 1-by-T cells)
            queue={};
            for ll=1:2
                queue{end+1}={lagnames{ll}};
            end
            while ~isempty(queue)
                fieldpath=queue{1};
                queue(1)=[];
                O=getfield(AllStatsPathO,fieldpath{:});
                P=getfield(AllStatsPathL,fieldpath{:});
                if isstruct(O)
                    Onames=fieldnames(O);
                    for kk=1:length(Onames)
                        queue{end+1}=[fieldpath,Onames(kk)];
                    end
                elseif iscell(O)
                    nanmismatch=nanmismatch+(~isequal(P,O));
                else
                    maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
                    nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
                end
            end
            fprintf('Stats, %s: (A5) _tminus1/_tplus1 (prices, a parameter, an aggregate) in the FnsToEvaluate vs the shifted paths, every stat, this should be zero: %.3e \n',modelname,maxdiff)
            fprintf('Stats, %s: (A5) _tminus1/_tplus1 (prices, a parameter, an aggregate) in the FnsToEvaluate vs the shifted paths, every stat, positions NaN (or text) in one but not the other, this should be zero: %i \n',modelname,nanmismatch)
            fprintf('Stats, %s: for scale, r_tminus1 differs from r by up to %.3e and w_tminus1 from w by up to %.3e over the path (must not be zero, or (A5) is empty) \n',modelname,max(abs(ParamPathO.rlag-PricePathX.r)),max(abs(ParamPathO.wlag-ParamPathX.w)))
        end
    end
end

output=struct();

end
