function output=CoreFHorzPType_CondlRest(n_d,n_d2_semiz,n_a,n_z,N_j,d_grid,d2_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,vfoptionsbaseline,simoptionsbaseline)
% Test: simoptions.conditionalrestrictions in every PType stats command, per type and grouped, against a brute force:
%   EvalFnOnAgentDist_AllStats_FHorz_Case1_PType                              every statistic
%   LifeCycleProfiles_FHorz_Case1_PType                                       every statistic, per age (the default), in age groups of five, in one age group
%   EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType                       Mean, StdDeviation, covariances, correlations (restrictions new 2026-10-07)
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType   the same per age group (restrictions new 2026-10-07)
%   EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType                          Mean and StdDeviation at each age; pass-through and pooled pairs
% and the RestrictedSampleMass outputs of all five, on three configs:
%   d_z_e            two types, DIFFERENT age weights per type (mewj a structure), ptweights 0.7/0.3
%   nod_noz_noe      three types, the middle one of zero weight (ptweights 0.6/0/0.4): it is not in any grouped pool, but its
%                    per-type (within-type) stats are still reported
%   nod1_z_noe_semiz two types, ptweights 0.4/0.6, with a restriction on the semi-exogenous state
%
% Restrictions (pages are the unrestricted output first, then the restrictions in fieldnames order):
%   positiveassets (a>0), working (agej<Jr; with parameters), atJ1 (agej==J1; mass at one age only), always (1),
%   employed (semiz==1; semiz config only), and for the CrossSection and AutoCorr commands also typespecific: a different
%   restriction for each type (a structure with a field per type), which AllStats and LifeCycleProfiles PType do not take.
%
% Brute force: every type's values come from the single-type ValuesOnGrid on that type's PType Policy and parameters, and its
% weights from its StationaryDist (mass one within the type, including that type's age weights). A type's restricted weights
% are its weights times the restriction indicator; the GROUPED population of a page is every type's (restricted) weights times
% its ptweight, concatenated. The statistics are then computed from their definitions, exactly as in the single-type test
% (CoreFHorzTests part 12, CoreFHorz_CrossSecCondlRest): quantiles as the smallest value with cumulative mass >= p (an exact tie
% to the lower value), Lorenz curve as the integral of the quantile function, Gini as half the mean absolute difference over
% the mean, means/variances/covariances as plain weighted sums. So the grouped stats are those of the pooled population, which is
% the definition, and a type of zero ptweight drops out of it automatically.
%
% Also: the 'always' page equals the unrestricted output, every leaf, for every command; a zero-mass restriction in the two
% CrossSection commands (warning, RestrictedSampleMass zero, NaN); whichcombos with the page and the type dimension in the two
% CrossSection commands; AutoCorr per-type output equals the single-type command on the type (every leaf, restricted pages too,
% and with the grid interpolation layer on the d_z_e config), and its grouped restricted pair outputs equal the pooling of the
% per-type pair outputs with weights ptweights times the pair mass.

Params.J1=floor(N_j/2);
agegroupings=1:5:N_j;
ngroups=length(agegroupings);
npts=100; nq=20; tol=10^(-12); % the defaults of AllStats and LifeCycleProfiles

n_d_withd=n_d; d_grid_withd=d_grid;
n_z_withz=n_z; z_grid_withz=z_grid; pi_z_withz=pi_z;

statnames={'Mean','RatioMeanToMedian','Median','MoreInequality.Percentile50th','MoreInequality.Percentile90th','MoreInequality.Percentile95th','MoreInequality.Percentile99th', ...
    'Variance','StdDeviation','Gini','LorenzCurve','Minimum','Maximum','QuantileCutoffs','QuantileMeans', ...
    'MoreInequality.Top1share','MoreInequality.Top5share','MoreInequality.Top10share','MoreInequality.Bottom50share'};
statgroup=[1,1,2,2,2,2,2,3,3,4,4,5,5,6,6,7,7,7,7];
groupnames={'Mean/RatioMeanToMedian','Median/percentiles','Variance/StdDeviation','Gini/LorenzCurve','Minimum/Maximum','QuantileCutoffs/QuantileMeans','Top/Bottom shares'};

output=struct();

for cc=1:3
    %% Setup of this config
    Params_PT=Params;
    FnsToEvaluate=struct();
    CondlRest=struct();
    TypeSpecific=struct();
    NeverRest=struct();
    if cc==1
        configname='d_z_e';
        Names_i={'low','high'};
        ptw=[0.7;0.3];
        n_d=n_d_withd; d_grid=d_grid_withd;
        n_z=n_z_withz; z_grid=z_grid_withz; pi_z=pi_z_withz;
        vfoptions=struct();
        vfoptions.n_e=vfoptionsbaseline.n_e; vfoptions.e_grid=vfoptionsbaseline.e_grid; vfoptions.pi_e=vfoptionsbaseline.pi_e;
        simoptions=struct();
        simoptions.n_e=simoptionsbaseline.n_e; simoptions.e_grid=simoptionsbaseline.e_grid; simoptions.pi_e=simoptionsbaseline.pi_e;
        jequaloneDist=zeros(n_a,n_z,vfoptions.n_e,'gpuArray');
        jequaloneDist(1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;
        Params_PT.kappa_j_pt=[Params.kappa_j; 1.5*Params.kappa_j];
        Params_PT.mewj=struct();
        Params_PT.mewj.low=Params.mewj;
        Params_PT.mewj.high=0.9.^(0:N_j-1); Params_PT.mewj.high=Params_PT.mewj.high/sum(Params_PT.mewj.high);
        ReturnFn=@(d,aprime,a,z,e,r,w,kappa_j_pt,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j_pt,sigma,eta,varphi,agej,Jr,pension);
        FnsToEvaluate.assets=@(d,aprime,a,z,e) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,e,w,kappa_j_pt) w*kappa_j_pt*z*e*d;
        FnsToEvaluate.consumption=@(d,aprime,a,z,e,r,w,kappa_j_pt) (1+r)*a+w*kappa_j_pt*z*e*d-aprime;
        CondlRest.positiveassets=@(d,aprime,a,z,e) (a>0);
        CondlRest.working=@(d,aprime,a,z,e,agej,Jr) (agej<Jr);
        CondlRest.atJ1=@(d,aprime,a,z,e,agej,J1) (agej==J1);
        CondlRest.always=@(d,aprime,a,z,e) 1;
        TypeSpecific.low=@(d,aprime,a,z,e) (a>0.5);
        TypeSpecific.high=@(d,aprime,a,z,e,agej,Jr) (agej>=Jr);
        NeverRest.never=@(d,aprime,a,z,e) 0;
    elseif cc==2
        configname='nod_noz_noe';
        Names_i={'low','ghost','high'};
        ptw=[0.6;0;0.4];
        n_d=0; d_grid=[];
        n_z=0; z_grid=[]; pi_z=[];
        vfoptions=struct();
        simoptions=struct();
        jequaloneDist=zeros(n_a,1,'gpuArray');
        jequaloneDist(1)=1;
        Params_PT.kappa_j_pt=[Params.kappa_j; 1.2*Params.kappa_j; 1.5*Params.kappa_j];
        ReturnFn=@(aprime,a,r,w,kappa_j_pt,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j_pt,sigma,agej,Jr,pension);
        FnsToEvaluate.assets=@(aprime,a) a;
        FnsToEvaluate.earnings=@(aprime,a,w,kappa_j_pt) w*kappa_j_pt;
        FnsToEvaluate.savings=@(aprime,a) aprime-a;
        CondlRest.positiveassets=@(aprime,a) (a>0);
        CondlRest.working=@(aprime,a,agej,Jr) (agej<Jr);
        CondlRest.atJ1=@(aprime,a,agej,J1) (agej==J1);
        CondlRest.always=@(aprime,a) 1;
        TypeSpecific.low=@(aprime,a) (a>0.5);
        TypeSpecific.ghost=@(aprime,a) (a>0.2);
        TypeSpecific.high=@(aprime,a,agej,Jr) (agej>=Jr);
        NeverRest.never=@(aprime,a) 0;
    else
        configname='nod1_z_noe_semiz';
        Names_i={'low','high'};
        ptw=[0.4;0.6];
        n_d=n_d2_semiz; d_grid=d2_grid_semiz;
        n_z=n_z_withz; z_grid=z_grid_withz; pi_z=pi_z_withz;
        vfoptions=struct();
        vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
        vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
        vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
        simoptions=struct();
        simoptions.n_semiz=simoptionsbaseline.n_semiz;
        simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
        simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
        simoptions.d_grid=d_grid;
        jequaloneDist=zeros(n_a,vfoptions.n_semiz,n_z,'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_semiz/2),ceil(n_z/2))=1;
        Params_PT.kappa_j_pt=[Params.kappa_j; 1.5*Params.kappa_j];
        ReturnFn=@(d2,aprime,a,semiz,z,r,w,kappa_j_pt,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j_pt,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,z) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,w,kappa_j_pt) w*kappa_j_pt*z*semiz;
        FnsToEvaluate.consumption=@(d2,aprime,a,semiz,z,r,w,kappa_j_pt,uempbenefit) (1+r)*a+w*kappa_j_pt*z*semiz+uempbenefit*(1-semiz)-aprime;
        CondlRest.positiveassets=@(d2,aprime,a,semiz,z) (a>0);
        CondlRest.working=@(d2,aprime,a,semiz,z,agej,Jr) (agej<Jr);
        CondlRest.atJ1=@(d2,aprime,a,semiz,z,agej,J1) (agej==J1);
        CondlRest.always=@(d2,aprime,a,semiz,z) 1;
        CondlRest.employed=@(d2,aprime,a,semiz,z) (semiz==1);
        TypeSpecific.low=@(d2,aprime,a,semiz,z) (semiz==1);
        TypeSpecific.high=@(d2,aprime,a,semiz,z) (a>0.5);
        NeverRest.never=@(d2,aprime,a,semiz,z) 0;
    end
    fprintf('\n-- CondlRest PType, config %s --\n',configname)
    N_i=length(Names_i);
    Params_PT.(PTypeDistParamNames{1})=ptw;
    FnNames=fieldnames(FnsToEvaluate);
    nFns=length(FnNames);
    RestrNames=fieldnames(CondlRest);
    nRestr=length(RestrNames);
    npages=1+nRestr; % AllStats and LifeCycleProfiles
    CondlRestX=CondlRest;
    CondlRestX.typespecific=TypeSpecific;
    RestrNamesX=fieldnames(CondlRestX);
    npagesX=npages+1; % the CrossSection and AutoCorr commands (typespecific is the last page)
    PageNamesX=[{'unrestricted'};RestrNamesX];

    %% Solve
    [~,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_PT,DiscountFactorParamNames,vfoptions);
    StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_PT,simoptions);

    %% Brute-force inputs, type by type: SD{ii} (N_ss,N_j), V{ii} (N_ss,N_j,nFns), R{ii} (N_ss,N_j,number of restrictions in CondlRestX)
    SD=cell(N_i,1); V=cell(N_i,1); R=cell(N_i,1); Params_type=cell(N_i,1); CondlRestX_type=cell(N_i,1);
    for ii=1:N_i
        nm=Names_i{ii};
        Params_ii=Params_PT;
        Params_ii.kappa_j_pt=Params_PT.kappa_j_pt(ii,:);
        if isstruct(Params_PT.mewj)
            Params_ii.mewj=Params_PT.mewj.(nm);
        end
        Params_type{ii}=Params_ii;
        CondlRestX_ii=struct(); % the restrictions of this type (typespecific is its own function)
        for rr=1:length(RestrNamesX)
            if isstruct(CondlRestX.(RestrNamesX{rr}))
                CondlRestX_ii.(RestrNamesX{rr})=CondlRestX.(RestrNamesX{rr}).(nm);
            else
                CondlRestX_ii.(RestrNamesX{rr})=CondlRestX.(RestrNamesX{rr});
            end
        end
        CondlRestX_type{ii}=CondlRestX_ii;
        VOG=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_PT.(nm),FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
        ROG=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_PT.(nm),CondlRestX_ii,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
        SD{ii}=reshape(gather(double(StationaryDist_PT.(nm))),[],N_j);
        N_ss=size(SD{ii},1);
        V{ii}=zeros(N_ss,N_j,nFns);
        for ff=1:nFns
            V{ii}(:,:,ff)=reshape(gather(double(VOG.(FnNames{ff}))),[N_ss,N_j]);
        end
        R{ii}=zeros(N_ss,N_j,length(RestrNamesX));
        for rr=1:length(RestrNamesX)
            R{ii}(:,:,rr)=double(reshape(gather(ROG.(RestrNamesX{rr})),[N_ss,N_j])~=0);
        end
    end
    groups5=cell(1,ngroups);
    for kk=1:ngroups
        if kk<ngroups
            groups5{kk}=agegroupings(kk):agegroupings(kk+1)-1;
        else
            groups5{kk}=agegroupings(kk):N_j;
        end
    end

    simoptions_R=simoptions; simoptions_R.conditionalrestrictions=CondlRest;
    simoptions_RX=simoptions; simoptions_RX.conditionalrestrictions=CondlRestX;

    %% ===================== AllStats and LifeCycleProfiles PType: every statistic, per type and grouped, vs brute force =====================
    AS=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_R);
    LCP1=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_R);
    simoptions_R5=simoptions_R; simoptions_R5.agegroupings=agegroupings;
    LCP5=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_R5);
    simoptions_Rone=simoptions_R; simoptions_Rone.agegroupings=1;
    LCPone=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_Rone);
    for cmd=1:4
        if cmd==1
            cmdname='AllStats PType'; Out=AS; groups={1:N_j};
        elseif cmd==2
            cmdname='LifeCycleProfiles PType per age'; Out=LCP1; groups=num2cell(1:N_j);
        elseif cmd==3
            cmdname='LifeCycleProfiles PType age groups of five'; Out=LCP5; groups=groups5;
        else
            cmdname='LifeCycleProfiles PType one age group'; Out=LCPone; groups={1:N_j};
        end
        for tt=1:N_i+1 % the types, then grouped
            if tt<=N_i
                unitname=['type ',Names_i{tt}];
            else
                unitname='grouped';
            end
            maxd=zeros(1,7); mism=zeros(1,7); worst=repmat({''},1,7);
            for pp=1:npages
                if pp==1
                    OutP=Out;
                else
                    OutP=Out.(RestrNames{pp-1});
                end
                for ff=1:nFns
                    for kk=1:numel(groups)
                        % the weights and values of this unit (one type, or every type times its ptweight) on this page and age group
                        v=zeros(0,1); w=zeros(0,1);
                        for ii=1:N_i
                            if tt<=N_i && ii~=tt
                                continue
                            end
                            if pp==1
                                Wii=SD{ii}(:,groups{kk});
                            else
                                Wii=SD{ii}(:,groups{kk}).*R{ii}(:,groups{kk},pp-1);
                            end
                            if tt>N_i
                                Wii=ptw(ii)*Wii;
                            end
                            v=[v; reshape(V{ii}(:,groups{kk},ff),[],1)]; %#ok<AGROW>
                            w=[w; Wii(:)]; %#ok<AGROW>
                        end
                        keep=(w>0); w=w(keep); v=v(keep);
                        B=cell(1,numel(statnames));
                        if isempty(w)
                            for ss=1:numel(statnames)
                                if strcmp(statnames{ss},'LorenzCurve')
                                    B{ss}=nan(npts,1);
                                elseif strcmp(statnames{ss},'QuantileCutoffs')
                                    B{ss}=nan(nq+1,1);
                                elseif strcmp(statnames{ss},'QuantileMeans')
                                    B{ss}=nan(nq,1);
                                else
                                    B{ss}=NaN;
                                end
                            end
                        else
                            w=w/sum(w);
                            [u,~,ic]=unique(v);
                            m=accumarray(ic,w);
                            F=cumsum(m); cS=cumsum(u.*m); S=cS(end); Ftot=F(end);
                            pq=[0.5; tol; 1-tol; 0.99; 0.95; 0.90; (1:nq-1)'/nq];
                            dq=10^(-12)*ones(size(pq)); dq(2:3)=0;
                            Qv=zeros(numel(pq),1);
                            for ii=1:numel(pq)
                                k=find(F>=pq(ii)-dq(ii),1);
                                if isempty(k)
                                    k=numel(u);
                                end
                                Qv(ii)=u(k);
                            end
                            pl=[(1:npts-1)'/npts; 0.99; 0.95; 0.90; 0.5; (1:nq-1)'/nq];
                            Lv=zeros(numel(pl),1);
                            for ii=1:numel(pl)
                                k=find(F>=pl(ii),1);
                                if isempty(k)
                                    k=numel(u);
                                end
                                Lv(ii)=(k>1)*cS(max(k-1,1))+(pl(ii)-(k>1)*F(max(k-1,1)))*u(k);
                            end
                            Mean=sum(u.*m);
                            Var=sum(m.*(u-Mean).^2);
                            if numel(u)==1
                                Gini=0; Lorenz=(1:npts)'/npts; shares=[0.01;0.05;0.1;0.5];
                            elseif u(1)<0
                                Gini=NaN; Lorenz=nan(npts,1); shares=nan(4,1);
                            else
                                Gini=sum(m.*u.*([0;F(1:end-1)]-(Ftot-F)))/(S*Ftot);
                                Lorenz=[Lv(1:npts-1);S]/S;
                                shares=[1-Lv(npts)/S; 1-Lv(npts+1)/S; 1-Lv(npts+2)/S; Lv(npts+3)/S];
                            end
                            B={Mean, Mean/Qv(1), Qv(1), Qv(1), Qv(6), Qv(5), Qv(4), ...
                                Var, sqrt(Var), Gini, Lorenz, Qv(2), Qv(3), [Qv(2);Qv(7:end);Qv(3)], diff([0;Lv(npts+4:end);S])*nq, ...
                                shares(1), shares(2), shares(3), shares(4)};
                        end
                        % compare with the command's output for this unit
                        for ss=1:numel(statnames)
                            gg=statgroup(ss);
                            if tt<=N_i
                                if isfield(OutP.(FnNames{ff}),Names_i{tt})
                                    X=OutP.(FnNames{ff}).(Names_i{tt});
                                else
                                    mism(gg)=mism(gg)+1;
                                    continue
                                end
                            else
                                X=OutP.(FnNames{ff});
                            end
                            parts=strsplit(statnames{ss},'.');
                            present=true;
                            for ii=1:numel(parts)
                                if isstruct(X) && isfield(X,parts{ii})
                                    X=X.(parts{ii});
                                else
                                    present=false;
                                end
                            end
                            if ~present
                                mism(gg)=mism(gg)+1;
                                continue
                            end
                            X=gather(X);
                            a=X(:,kk); b=B{ss}(:);
                            if numel(a)~=numel(b)
                                mism(gg)=mism(gg)+1;
                                continue
                            end
                            mism(gg)=mism(gg)+sum(isnan(a)~=isnan(b))+sum(xor(isinf(a),isinf(b)));
                            d=abs(a-b); d=d(~isnan(d));
                            if ~isempty(d) && max(d)>maxd(gg)
                                maxd(gg)=max(d);
                                worst{gg}=sprintf('%s, page %s, %s, group %i',statnames{ss},PageNamesX{pp},FnNames{ff},kk);
                            end
                        end
                    end
                end
            end
            for gg=1:7
                fprintf('CondlRest PType %s %s, %s, %s vs brute force on every page, this should be zero: %.3e, NaN/Inf mismatches + missing fields %i (worst: %s) \n',cmdname,configname,unitname,groupnames{gg},maxd(gg),mism(gg),worst{gg})
            end
        end
    end
    % Restricted sample masses: AllStats .(type) a within-type share and .TotalAllPTypes a population share; LifeCycleProfiles
    % .(type) the within-type mass at each age (it includes the type's age weights), .ByAge sum_ii ptweights(ii)*that, .ByPType
    % its sum over ages, .Total the population share
    dA=0; dL=0;
    for rr=1:nRestr
        RSMtype=zeros(N_i,N_j);
        for ii=1:N_i
            RSMtype(ii,:)=sum(SD{ii}.*R{ii}(:,:,rr),1);
            dA=max(dA,abs(gather(AS.(RestrNames{rr}).RestrictedSampleMass.(Names_i{ii}))-sum(RSMtype(ii,:))));
            dL=max(dL,max(abs(reshape(gather(LCP1.(RestrNames{rr}).RestrictedSampleMass.(Names_i{ii})),1,[])-RSMtype(ii,:))));
        end
        dA=max(dA,abs(gather(AS.(RestrNames{rr}).RestrictedSampleMass.TotalAllPTypes)-sum(ptw.*sum(RSMtype,2))));
        dL=max([dL,max(abs(reshape(gather(LCP1.(RestrNames{rr}).RestrictedSampleMass.ByAge),1,[])-sum(ptw.*RSMtype,1))), ...
            max(abs(reshape(gather(LCP1.(RestrNames{rr}).RestrictedSampleMass.ByPType),[],1)-sum(RSMtype,2))), ...
            abs(gather(LCP1.(RestrNames{rr}).RestrictedSampleMass.Total)-sum(ptw.*sum(RSMtype,2)))]);
    end
    fprintf('CondlRest PType AllStats %s: RestrictedSampleMass per type and TotalAllPTypes vs brute force, this should be zero: %.3e \n',configname,dA)
    fprintf('CondlRest PType LifeCycleProfiles %s: RestrictedSampleMass per type, ByAge, ByPType, Total vs brute force, this should be zero: %.3e \n',configname,dL)

    %% ===================== CrossSection PType and AgeConditional CrossSection PType: per type and grouped vs brute force =====================
    CS=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_RX);
    simoptions_RX5=simoptions_RX; simoptions_RX5.agegroupings=agegroupings;
    ACS=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_RX5);
    for cmd=1:2
        if cmd==1
            cmdname='CrossSection PType'; Out=CS; groups={1:N_j};
        else
            cmdname='AgeConditional CrossSection PType'; Out=ACS; groups=groups5;
        end
        ng=numel(groups);
        for tt=1:N_i+1
            if tt<=N_i
                unitname=['type ',Names_i{tt}];
            else
                unitname='grouped';
            end
            dm=0; dc=0; dnamed=0; nanmis=0; nok=0; nall=0;
            for pp=1:npagesX
                if pp==1
                    OutP=Out;
                else
                    OutP=Out.(RestrNamesX{pp-1});
                end
                % the command's matrices for this unit (the per-type matrices cover the functions relevant to the type, here all of them)
                if tt<=N_i
                    Mt=gather(OutP.CovarianceMatrix_ptype.(Names_i{tt})); Rt=gather(OutP.CorrelationMatrix_ptype.(Names_i{tt}));
                else
                    Mt=gather(OutP.CovarianceMatrix); Rt=gather(OutP.CorrelationMatrix);
                end
                Mt=reshape(Mt,[nFns,nFns,ng]); Rt=reshape(Rt,[nFns,nFns,ng]);
                for kk=1:ng
                    v=zeros(0,nFns); w=zeros(0,1);
                    for ii=1:N_i
                        if tt<=N_i && ii~=tt
                            continue
                        end
                        if pp==1
                            Wii=SD{ii}(:,groups{kk});
                        else
                            Wii=SD{ii}(:,groups{kk}).*R{ii}(:,groups{kk},pp-1);
                        end
                        if tt>N_i
                            Wii=ptw(ii)*Wii;
                        end
                        v=[v; reshape(V{ii}(:,groups{kk},:),[],nFns)]; %#ok<AGROW>
                        w=[w; Wii(:)]; %#ok<AGROW>
                    end
                    mass=sum(w);
                    if mass==0
                        Bm=nan(nFns,1); Bs=nan(nFns,1); BC=nan(nFns,nFns); BR=nan(nFns,nFns); ok=false(nFns,nFns);
                    else
                        w=w/mass;
                        Bm=(v'*w);
                        dev=v-Bm';
                        Bs=sqrt((dev.^2)'*w);
                        BC=dev'*(dev.*w);
                        BR=BC./(Bs*Bs'); BR(logical(eye(nFns)))=1;
                        sok=(Bs>1e-8*max(1,abs(Bm)));
                        ok=(sok*sok')>0; ok(logical(eye(nFns)))=true;
                    end
                    % Mean/StdDeviation (named), covariance matrix, correlation matrix where both std devs are not rounding noise
                    for ff=1:nFns
                        if tt<=N_i
                            if isfield(OutP,FnNames{ff}) && isfield(OutP.(FnNames{ff}),Names_i{tt})
                                P=OutP.(FnNames{ff}).(Names_i{tt});
                            else
                                nanmis=nanmis+1;
                                continue
                            end
                        else
                            P=OutP.(FnNames{ff});
                        end
                        a=[gather(P.Mean(kk));gather(P.StdDeviation(kk))]; b=[Bm(ff);Bs(ff)];
                        nanmis=nanmis+sum(isnan(a)~=isnan(b)); d=abs(a-b); dm=max([dm;d(~isnan(d))]);
                        % named fields vs the matrices
                        dnamed=dnamed+(isnan(gather(P.(FnNames{ff})(kk)))~=isnan(Rt(ff,ff,kk)));
                        for ff2=1:nFns
                            if ff2~=ff
                                a=gather(P.CovarianceWith.(FnNames{ff2})(kk)); b=Mt(ff,ff2,kk);
                                dnamed=dnamed+(isnan(a)~=isnan(b)); if ~isnan(a) && ~isnan(b), dnamed=max(dnamed,abs(a-b)); end
                                a=gather(P.CorrelationWith.(FnNames{ff2})(kk)); b=Rt(ff,ff2,kk);
                                dnamed=dnamed+(isnan(a)~=isnan(b)); if ~isnan(a) && ~isnan(b), dnamed=max(dnamed,abs(a-b)); end
                            end
                        end
                    end
                    a=reshape(Mt(:,:,kk),[],1); b=BC(:);
                    nanmis=nanmis+sum(isnan(a)~=isnan(b)); d=abs(a-b); dm=max([dm;d(~isnan(d))]);
                    a=reshape(Rt(:,:,kk),[],1); b=BR(:);
                    nanmis=nanmis+sum(isnan(a(ok(:))));
                    d=abs(a(ok(:))-b(ok(:))); dc=max([dc;d(~isnan(d))]);
                    nok=nok+sum(ok(:)); nall=nall+numel(ok);
                end
            end
            fprintf('CondlRest PType %s %s, %s, every page: Mean/StdDeviation/CovarianceMatrix vs brute force, this should be zero: %.3e, NaN mismatches %i \n',cmdname,configname,unitname,dm,nanmis)
            fprintf('CondlRest PType %s %s, %s, every page: CorrelationMatrix vs brute force (where both std devs are nonzero, %i of %i entries), this should be zero: %.3e \n',cmdname,configname,unitname,nok,nall,dc)
            fprintf('CondlRest PType %s %s, %s, every page: named fields vs matrices, this should be zero: %.3e \n',cmdname,configname,unitname,dnamed)
        end
    end
    % Restricted sample masses of the two CrossSection commands
    dC=0; dA=0;
    for rr=1:length(RestrNamesX)
        RSMtype=zeros(N_i,N_j);
        for ii=1:N_i
            RSMtype(ii,:)=sum(SD{ii}.*R{ii}(:,:,rr),1);
            dC=max(dC,abs(gather(CS.(RestrNamesX{rr}).RestrictedSampleMass.(Names_i{ii}))-sum(RSMtype(ii,:))));
            dA=max(dA,max(abs(reshape(gather(ACS.(RestrNamesX{rr}).RestrictedSampleMass.(Names_i{ii})),1,[])-RSMtype(ii,:))));
        end
        dC=max([dC,max(abs(reshape(gather(CS.(RestrNamesX{rr}).RestrictedSampleMass.ByPType),[],1)-sum(RSMtype,2))),abs(gather(CS.(RestrNamesX{rr}).RestrictedSampleMass.Total)-sum(ptw.*sum(RSMtype,2)))]);
        dA=max([dA,max(abs(reshape(gather(ACS.(RestrNamesX{rr}).RestrictedSampleMass.ByAge),1,[])-sum(ptw.*RSMtype,1))), ...
            max(abs(reshape(gather(ACS.(RestrNamesX{rr}).RestrictedSampleMass.ByPType),[],1)-sum(RSMtype,2))), ...
            abs(gather(ACS.(RestrNamesX{rr}).RestrictedSampleMass.Total)-sum(ptw.*sum(RSMtype,2)))]);
    end
    fprintf('CondlRest PType CrossSection %s: RestrictedSampleMass per type, ByPType, Total vs brute force, this should be zero: %.3e \n',configname,dC)
    fprintf('CondlRest PType AgeConditional CrossSection %s: RestrictedSampleMass per type, ByAge, ByPType, Total vs brute force, this should be zero: %.3e \n',configname,dA)

    %% ===================== AutoCorr PType: Mean and StdDeviation at each age, per type and grouped, vs brute force =====================
    simoptions_acp=simoptions_RX; simoptions_acp.timehorizons=[1,2]; % horizons 1 and 2 (since 2026-10-14 exactly the horizons listed are computed)
    ACP=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_acp);
    for tt=1:N_i+1
        if tt<=N_i
            unitname=['type ',Names_i{tt}];
        else
            unitname='grouped';
        end
        dm=0; nanmis=0;
        for pp=1:npagesX
            if pp==1
                OutP=ACP;
            else
                OutP=ACP.(RestrNamesX{pp-1});
            end
            for ff=1:nFns
                Bm=nan(1,N_j); Bs=nan(1,N_j);
                for jj=1:N_j
                    v=zeros(0,1); w=zeros(0,1);
                    for ii=1:N_i
                        if tt<=N_i && ii~=tt
                            continue
                        end
                        if pp==1
                            Wii=SD{ii}(:,jj);
                        else
                            Wii=SD{ii}(:,jj).*R{ii}(:,jj,pp-1);
                        end
                        if tt>N_i
                            Wii=ptw(ii)*Wii;
                        end
                        v=[v; V{ii}(:,jj,ff)]; %#ok<AGROW>
                        w=[w; Wii]; %#ok<AGROW>
                    end
                    if sum(w)>0
                        w=w/sum(w);
                        Bm(jj)=sum(w.*v);
                        Bs(jj)=sqrt(sum(w.*(v-Bm(jj)).^2));
                    end
                end
                if tt<=N_i
                    P=OutP.(FnNames{ff}).(Names_i{tt});
                else
                    P=OutP.(FnNames{ff});
                end
                a=[reshape(gather(P.Mean),1,[]),reshape(gather(P.StdDeviation),1,[])]; b=[Bm,Bs];
                if numel(a)~=numel(b)
                    nanmis=nanmis+1;
                else
                    nanmis=nanmis+sum(isnan(a)~=isnan(b)); d=abs(a-b); dm=max([dm,d(~isnan(d))]);
                end
            end
        end
        fprintf('CondlRest PType AutoCorr %s, %s, every page: Mean/StdDeviation at each age vs brute force, this should be zero: %.3e, NaN mismatches %i \n',configname,unitname,dm,nanmis)
    end
    dR=0;
    for rr=1:length(RestrNamesX)
        RSMtype=zeros(N_i,N_j);
        for ii=1:N_i
            RSMtype(ii,:)=sum(SD{ii}.*R{ii}(:,:,rr),1);
            dR=max(dR,max(abs(reshape(gather(ACP.(RestrNamesX{rr}).RestrictedSampleMass.(Names_i{ii})),1,[])-RSMtype(ii,:))));
        end
        dR=max([dR,max(abs(reshape(gather(ACP.(RestrNamesX{rr}).RestrictedSampleMass.ByAge),1,[])-sum(ptw.*RSMtype,1))), ...
            max(abs(reshape(gather(ACP.(RestrNamesX{rr}).RestrictedSampleMass.ByPType),[],1)-sum(RSMtype,2))), ...
            abs(gather(ACP.(RestrNamesX{rr}).RestrictedSampleMass.Total)-sum(ptw.*sum(RSMtype,2)))]);
    end
    fprintf('CondlRest PType AutoCorr %s: RestrictedSampleMass per type, ByAge, ByPType, Total vs brute force, this should be zero: %.3e \n',configname,dR)

    % Per type: the single-type command on the type, every leaf (restricted pages included); with GI as well on the d_z_e config
    for gi=0:(cc==1)
        if gi==0
            Policy_c=Policy_PT; Dist_c=StationaryDist_PT; simoptions_c=simoptions_acp; ACP_c=ACP; gistr='';
        else
            vfoptions_gi=vfoptions; vfoptions_gi.gridinterplayer=1; vfoptions_gi.ngridinterp=5;
            simoptions_c=simoptions_acp; simoptions_c.gridinterplayer=1; simoptions_c.ngridinterp=5;
            [~,Policy_c]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_PT,DiscountFactorParamNames,vfoptions_gi);
            Dist_c=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_c,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_PT,simoptions_c);
            ACP_c=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(Dist_c,Policy_c,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_c);
            gistr=' GI';
        end
        ACP_solo=cell(N_i,1);
        maxd=0; mis=0;
        for ii=1:N_i
            nm=Names_i{ii};
            simoptions_ii=simoptions_c; simoptions_ii.conditionalrestrictions=CondlRestX_type{ii};
            ACP_solo{ii}=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_c.(nm),Policy_c.(nm),FnsToEvaluate,Params_type{ii},[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_ii);
            for pp=1:npagesX
                for ff=1:nFns
                    if pp==1
                        [p0,v0]=CoreStructLeaves(ACP_solo{ii}.(FnNames{ff})); [p1,v1]=CoreStructLeaves(ACP_c.(FnNames{ff}).(nm));
                    else
                        [p0,v0]=CoreStructLeaves(ACP_solo{ii}.(RestrNamesX{pp-1}).(FnNames{ff})); [p1,v1]=CoreStructLeaves(ACP_c.(RestrNamesX{pp-1}).(FnNames{ff}).(nm));
                    end
                    for k=1:numel(p0)
                        [present,idx]=ismember(p0{k},p1);
                        if present && isequal(size(v0{k}),size(v1{idx}))
                            a=v0{k}; b=v1{idx}; mis=mis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                        else
                            mis=mis+1;
                        end
                    end
                end
            end
        end
        fprintf('CondlRest PType AutoCorr%s %s: per type vs the single-type command on that type, every page, every leaf, this should be zero: %.3e, mismatches %i \n',gistr,configname,maxd,mis)
    end

    % Grouped restricted pair outputs (horizons 1 and 2) vs the pooling of the per-type pair outputs, weights ptweights times the pair mass (ACP_solo is now
    % the last of the loop above: without GI on configs 2-3, with GI on config 1, so pool against the matching ACP_c)
    dmx=0; dcov=0; dcor=0; dpm=0; nanmis=0;
    hstrs={'','_horizon2'}; horizons=[1,2]; % (timehorizons=[1,2] above)
    for pp=2:npagesX
        rname=RestrNamesX{pp-1};
        for ff=1:nFns
            fn=FnNames{ff};
            for hh=1:numel(horizons)
                T=N_j-horizons(hh);
                PM=zeros(T,N_i); MX=zeros(T,N_i); MY=zeros(T,N_i); SX=zeros(T,N_i); SY=zeros(T,N_i); C=zeros(T,N_i);
                for ii=1:N_i
                    PM(:,ii)=ptw(ii)*reshape(gather(ACP_solo{ii}.(rname).(fn).(['PairMass',hstrs{hh}])),[],1);
                    MX(:,ii)=reshape(gather(ACP_solo{ii}.(rname).(fn).(['PairMean_j',hstrs{hh}])),[],1);
                    MY(:,ii)=reshape(gather(ACP_solo{ii}.(rname).(fn).(['PairMean_jplusk',hstrs{hh}])),[],1);
                    SX(:,ii)=reshape(gather(ACP_solo{ii}.(rname).(fn).(['PairStdDeviation_j',hstrs{hh}])),[],1);
                    SY(:,ii)=reshape(gather(ACP_solo{ii}.(rname).(fn).(['PairStdDeviation_jplusk',hstrs{hh}])),[],1);
                    C(:,ii)=reshape(gather(ACP_solo{ii}.(rname).(fn).(['AutoCovariance',hstrs{hh}])),[],1);
                end
                PM(isnan(PM))=0;
                none=(PM==0); MX(none)=0; MY(none)=0; SX(none)=0; SY(none)=0; C(none)=0; % a type with no pairs contributes nothing (its moments are NaN)
                Wt=PM./sum(PM,2);
                mx=sum(Wt.*MX,2); my=sum(Wt.*MY,2);
                vx=sum(Wt.*(SX.^2+(MX-mx).^2),2); vy=sum(Wt.*(SY.^2+(MY-my).^2),2);
                cG=sum(Wt.*(C+(MX-mx).*(MY-my)),2);
                rG=cG./sqrt(vx.*vy);
                nopairs=(sum(PM,2)==0); mx(nopairs)=NaN; my(nopairs)=NaN; cG(nopairs)=NaN; rG(nopairs)=NaN;
                corrok=~nopairs & sqrt(vx.*vy)>1e-10; % correlations compared where the std devs are not rounding noise
                G=ACP_c.(rname).(fn);
                x=reshape(gather(G.(['PairMass',hstrs{hh}])),[],1); y=sum(PM,2); d=abs(x-y); dpm=max([dpm;d(~isnan(d))]); nanmis=nanmis+sum(isnan(x)~=isnan(y));
                x=[reshape(gather(G.(['PairMean_j',hstrs{hh}])),[],1);reshape(gather(G.(['PairMean_jplusk',hstrs{hh}])),[],1)]; y=[mx;my]; d=abs(x-y); dmx=max([dmx;d(~isnan(d))]); nanmis=nanmis+sum(isnan(x)~=isnan(y));
                x=reshape(gather(G.(['AutoCovariance',hstrs{hh}])),[],1); d=abs(x-cG); dcov=max([dcov;d(~isnan(d))]); nanmis=nanmis+sum(isnan(x)~=isnan(cG));
                x=reshape(gather(G.(['AutoCorrelation',hstrs{hh}])),[],1); d=abs(x(corrok)-rG(corrok)); dcor=max([dcor;d(~isnan(d))]); nanmis=nanmis+sum(isnan(x(corrok)));
            end
        end
    end
    fprintf('CondlRest PType AutoCorr%s %s: grouped restricted PairMass vs pooled per-type outputs, this should be zero: %.3e \n',gistr,configname,dpm)
    fprintf('CondlRest PType AutoCorr%s %s: grouped restricted PairMean_j/_jplusk vs pooled per-type outputs, this should be zero: %.3e \n',gistr,configname,dmx)
    fprintf('CondlRest PType AutoCorr%s %s: grouped restricted AutoCovariance vs pooled per-type outputs, this should be zero: %.3e \n',gistr,configname,dcov)
    fprintf('CondlRest PType AutoCorr%s %s: grouped restricted AutoCorrelation vs pooled per-type outputs (where the std devs are nonzero), this should be zero: %.3e, NaN mismatches %i \n',gistr,configname,dcor,nanmis)

    %% ===================== The always page equals the unrestricted output, every leaf, every command =====================
    for cmd=1:5
        if cmd==1
            cmdname='AllStats PType'; Out=AS;
        elseif cmd==2
            cmdname='LifeCycleProfiles PType'; Out=LCP1;
        elseif cmd==3
            cmdname='CrossSection PType'; Out=CS;
        elseif cmd==4
            cmdname='AgeConditional CrossSection PType'; Out=ACS;
        else
            cmdname='AutoCorr PType'; Out=ACP;
        end
        maxd=0; mis=0;
        names0=FnNames;
        if cmd==3 || cmd==4
            names0=[FnNames;{'CovarianceMatrix';'CorrelationMatrix';'CovarianceMatrix_ptype';'CorrelationMatrix_ptype'}];
        end
        for ff=1:numel(names0)
            if ~isfield(Out,names0{ff})
                continue % (e.g. the AutoCorr command has no matrices)
            end
            if isstruct(Out.(names0{ff}))
                [p0,v0]=CoreStructLeaves(Out.(names0{ff})); [pA,vA]=CoreStructLeaves(Out.always.(names0{ff}));
            else
                p0={''}; v0={gather(Out.(names0{ff}))}; pA={''}; vA={gather(Out.always.(names0{ff}))};
            end
            for k=1:numel(p0)
                [present,idx]=ismember(p0{k},pA);
                if present && isequal(size(v0{k}),size(vA{idx}))
                    a=v0{k}; b=vA{idx}; mis=mis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                else
                    mis=mis+1;
                end
            end
            if cmd~=5 % (the restricted AutoCorr output also has the pair fields, PairMass, PairMean_j, ..., which the unrestricted one does not)
                mis=mis+sum(~ismember(pA,p0));
            end
        end
        fprintf('CondlRest PType %s %s: page always vs the unrestricted output, every leaf, this should be zero: %.3e, mismatches %i \n',cmdname,configname,maxd,mis)
    end

    %% ===================== A zero-mass restriction in the two CrossSection commands =====================
    simoptions_N=simoptions; simoptions_N.conditionalrestrictions=NeverRest;
    lastwarn('');
    CSN=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_N);
    warned=~isempty(lastwarn);
    notnan=sum(~isnan(gather(CSN.never.CovarianceMatrix(:))))+sum(~isnan(gather(CSN.never.CorrelationMatrix(:))));
    for ff=1:nFns
        notnan=notnan+~isnan(gather(CSN.never.(FnNames{ff}).Mean))+~isnan(gather(CSN.never.(FnNames{ff}).StdDeviation))+~isnan(gather(CSN.never.(FnNames{ff}).(FnNames{ff})));
    end
    fprintf('CondlRest PType CrossSection %s: zero-mass restriction: no warning + RestrictedSampleMass.Total + grouped entries not NaN, this should be zero: %.3e \n',configname,(1-warned)+abs(gather(CSN.never.RestrictedSampleMass.Total))+notnan)
    simoptions_N5=simoptions_N; simoptions_N5.agegroupings=agegroupings;
    lastwarn('');
    ACSN=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_N5);
    warned=~isempty(lastwarn);
    notnan=sum(~isnan(gather(ACSN.never.CovarianceMatrix(:))))+sum(~isnan(gather(ACSN.never.CorrelationMatrix(:))));
    for ff=1:nFns
        notnan=notnan+sum(~isnan(gather(ACSN.never.(FnNames{ff}).Mean)))+sum(~isnan(gather(ACSN.never.(FnNames{ff}).StdDeviation)))+sum(~isnan(gather(ACSN.never.(FnNames{ff}).(FnNames{ff}))));
    end
    fprintf('CondlRest PType AgeConditional CrossSection %s: zero-mass restriction: no warning + RestrictedSampleMass.Total + grouped entries not NaN, this should be zero: %.3e \n',configname,(1-warned)+abs(gather(ACSN.never.RestrictedSampleMass.Total))+notnan)

    %% ===================== whichcombos with the page and the type dimension in the two CrossSection commands =====================
    % (a) a checkerboard over (pair, page, type slot): the grouped matrices equal the full run where the grouped slot selects, NaN
    % elsewhere; each type's matrices equal the full run where the type's own slot or the grouped slot selects, NaN elsewhere (the
    % diagonal wherever the function is evaluated); (b) inputs without the type or the page dimension equal their repmats; (c) a wrong
    % number of pages errors. The same for the age-conditional command with an age-group dimension.
    for cmd=1:2
        if cmd==1
            cmdname='CrossSection PType'; Full=CS; ng=1; simoptions_W=simoptions_RX;
        else
            cmdname='AgeConditional CrossSection PType'; Full=ACS; ng=ngroups; simoptions_W=simoptions_RX5;
        end
        wc=zeros(nFns,nFns,ng,npagesX,N_i+1);
        for ss=1:N_i+1
            for pp=1:npagesX
                for kk=1:ng
                    for ff1=1:nFns
                        for ff2=1:nFns
                            if ff1~=ff2
                                wc(ff1,ff2,kk,pp,ss)=mod(ff1+ff2+kk+pp+ss,2);
                            end
                        end
                    end
                end
            end
        end
        wc(:,:,:,:,N_i+1)=wc(:,:,:,:,N_i+1).*reshape(mod((1:npagesX),2),[1,1,1,npagesX]); % the grouped slot only on the odd pages
        if cmd==1
            simoptions_W.whichcombos=reshape(wc,[nFns,nFns,npagesX,N_i+1]);
            Sub=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_W);
        else
            simoptions_W.whichcombos=wc;
            Sub=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_W);
        end
        maxd=0; bad=0;
        for pp=1:npagesX
            if pp==1
                FP=Full; SP=Sub;
            else
                FP=Full.(RestrNamesX{pp-1}); SP=Sub.(RestrNamesX{pp-1});
            end
            for tt=1:N_i+1
                if tt<=N_i
                    sel=max(wc(:,:,:,pp,tt),wc(:,:,:,pp,N_i+1)); % the type's own slot, or forced by the grouped slot
                    Mf=reshape(gather(FP.CovarianceMatrix_ptype.(Names_i{tt})),[nFns,nFns,ng]); Rf=reshape(gather(FP.CorrelationMatrix_ptype.(Names_i{tt})),[nFns,nFns,ng]);
                    if ~isfield(SP,'CovarianceMatrix_ptype') || ~isfield(SP.CovarianceMatrix_ptype,Names_i{tt})
                        bad=bad+1;
                        continue
                    end
                    Ms=reshape(gather(SP.CovarianceMatrix_ptype.(Names_i{tt})),[nFns,nFns,ng]); Rs=reshape(gather(SP.CorrelationMatrix_ptype.(Names_i{tt})),[nFns,nFns,ng]);
                else
                    sel=wc(:,:,:,pp,N_i+1);
                    Mf=reshape(gather(FP.CovarianceMatrix),[nFns,nFns,ng]); Rf=reshape(gather(FP.CorrelationMatrix),[nFns,nFns,ng]);
                    Ms=reshape(gather(SP.CovarianceMatrix),[nFns,nFns,ng]); Rs=reshape(gather(SP.CorrelationMatrix),[nFns,nFns,ng]);
                end
                on=false(nFns,nFns,ng);
                for kk=1:ng
                    on(:,:,kk)=(max(sel(:,:,kk),diag(any(sel(:,:,kk),2)))==1);
                end
                d=[abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]; d=d(~isnan(d)); maxd=max([maxd;d]);
                bad=bad+sum(isnan(Mf(on))~=isnan(Ms(on)))+sum(isnan(Rf(on))~=isnan(Rs(on)))+sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on)));
            end
        end
        fprintf('CondlRest PType %s %s, whichcombos with the page and type dimensions: selected entries vs full run, this should be zero: %.3e, unselected not NaN + NaN mismatches %i \n',cmdname,configname,maxd,bad)
        % (b) broadcasts: without the type dimension = repmat over the type slots; without the page dimension = repmat over the pages
        wcb=wc(:,:,:,:,1);
        for bb=1:2
            if bb==1 % no type dimension
                if cmd==1
                    inA=reshape(wcb,[nFns,nFns,npagesX]); inB=repmat(reshape(wcb,[nFns,nFns,npagesX]),[1,1,1,N_i+1]);
                else
                    inA=wcb; inB=repmat(wcb,[1,1,1,1,N_i+1]);
                end
                bstr='without the type dimension vs its repmat over the type slots';
            else % no page and no type dimension
                if cmd==1
                    inA=wcb(:,:,1,1); inB=repmat(wcb(:,:,1,1),[1,1,npagesX,N_i+1]);
                else
                    inA=wcb(:,:,:,1); inB=repmat(wcb(:,:,:,1),[1,1,1,npagesX,N_i+1]);
                end
                bstr='without the page and type dimensions vs its repmat';
            end
            simoptions_W.whichcombos=inA;
            if cmd==1
                A=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_W);
            else
                A=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_W);
            end
            simoptions_W.whichcombos=inB;
            if cmd==1
                B=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_W);
            else
                B=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_W);
            end
            [pa,va]=CoreStructLeaves(A); [pb,vb]=CoreStructLeaves(B);
            maxd=0; mis=0;
            for k=1:numel(pa)
                [present,idx]=ismember(pa{k},pb);
                if present && isequal(size(va{k}),size(vb{idx}))
                    a=va{k}; b=vb{idx}; mis=mis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                else
                    mis=mis+1;
                end
            end
            fprintf('CondlRest PType %s %s, a whichcombos %s, every leaf, this should be zero: %.3e, mismatches %i \n',cmdname,configname,bstr,maxd,mis+sum(~ismember(pb,pa)))
        end
        % (c) the wrong number of pages
        errored=0;
        if cmd==1
            simoptions_W.whichcombos=ones(nFns,nFns,npagesX-1,N_i+1);
            try
                EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_W);
            catch
                errored=1;
            end
        else
            simoptions_W.whichcombos=ones(nFns,nFns,ngroups,npagesX-1,N_i+1);
            try
                EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_W);
            catch
                errored=1;
            end
        end
        fprintf('CondlRest PType %s %s, a whichcombos with the wrong number of pages must error, this should be zero: %.3e \n',cmdname,configname,1-errored)
    end

    output.(configname).AllStats=AS;
    output.(configname).LifeCycleProfiles=LCP1;
    output.(configname).CrossSection=CS;
    output.(configname).AgeConditionalCrossSection=ACS;
    output.(configname).AutoCorr=ACP;
end

end
