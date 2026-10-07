function output=CoreInfHorz_CondlRest(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% Test: simoptions.conditionalrestrictions in the InfHorz stats commands
%   EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz   output .(restriction).(fn).Mean etc, .CovarianceMatrix, .RestrictedSampleMass (new 2026-10-07)
%   EvalFnOnAgentDist_AllStats_InfHorz                every statistic
% on two models: with d and z (one endogenous state), and without d, with z and two endogenous states (assets and occupation,
% the with2A model of fig 11: the cross-section command reshapes the two endogenous states into one joint grid).
% (EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz does not take conditional restrictions.)
%
% Checks:
%   (1) every restricted (and the unrestricted) Mean, StdDeviation, CovarianceMatrix, CorrelationMatrix and RestrictedSampleMass of
%       the cross-section command vs a brute force from ValuesOnGrid and the agent distribution (the restricted distribution
%       normalised to mass one); correlations only where both brute-force std devs are not rounding noise
%   (2) its named fields equal its matrices; its restricted Mean/StdDeviation equal AllStats under the same restrictions
%   (3) every statistic of AllStats (Mean, RatioMeanToMedian, Median, percentiles, Variance, StdDeviation, Gini, LorenzCurve,
%       Minimum, Maximum, QuantileCutoffs, QuantileMeans, top/bottom shares) on every page vs a brute force written from the
%       definitions, as in CoreFHorzTests part 12 (quantile = smallest value with cumulative mass >= p, an exact tie to the lower
%       value; Lorenz curve = integral of the quantile function; Gini = half the mean absolute difference over the mean)
%   (4) the always page equals the unrestricted output, every leaf (both commands); the unrestricted cross-section output with
%       restrictions equals a run without them
%   (5) a restriction of zero mass in the cross-section command: warning, RestrictedSampleMass zero, NaN

npts=100; nq=20; tol=10^(-12); % the defaults of AllStats

statnames={'Mean','RatioMeanToMedian','Median','MoreInequality.Percentile50th','MoreInequality.Percentile90th','MoreInequality.Percentile95th','MoreInequality.Percentile99th', ...
    'Variance','StdDeviation','Gini','LorenzCurve','Minimum','Maximum','QuantileCutoffs','QuantileMeans', ...
    'MoreInequality.Top1share','MoreInequality.Top5share','MoreInequality.Top10share','MoreInequality.Bottom50share'};
statgroup=[1,1,2,2,2,2,2,3,3,4,4,5,5,6,6,7,7,7,7];
groupnames={'Mean/RatioMeanToMedian','Median/percentiles','Variance/StdDeviation','Gini/LorenzCurve','Minimum/Maximum','QuantileCutoffs/QuantileMeans','Top/Bottom shares'};

n_d_withd=n_d; d_grid_withd=d_grid; n_a_1A=n_a; a_grid_1A=a_grid; n_z_1A=n_z; z_grid_1A=z_grid; pi_z_1A=pi_z;

output=struct();

for cc=1:2
    %% Setup of this model
    Params_c=Params;
    vfoptions=struct();
    simoptions=struct();
    FnsToEvaluate=struct();
    CondlRest=struct();
    NeverRest=struct();
    if cc==1
        configname='d_z';
        n_d=n_d_withd; d_grid=d_grid_withd;
        n_a=n_a_1A; a_grid=a_grid_1A;
        n_z=n_z_1A; z_grid=z_grid_1A; pi_z=pi_z_1A;
        ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a,z) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,w) w*d*z;
        FnsToEvaluate.savings=@(d,aprime,a,z) aprime-a; % negative values: Gini, Lorenz curve and shares are NaN
        CondlRest.positiveassets=@(d,aprime,a,z) (a>0);
        CondlRest.highz=@(d,aprime,a,z) (z>1);
        CondlRest.highincome=@(d,aprime,a,z,r,w) (r*a+w*d*z>0.5); % a restriction with parameters
        CondlRest.always=@(d,aprime,a,z) 1;
        NeverRest.never=@(d,aprime,a,z) 0;
    else
        configname='nod_z_with2A';
        % The with2A model of fig 11 (simplified Kitao 2008, no taxes, partial equilibrium), on its own grids
        n_asset=101; n_occ=2; n_eta=3; n_theta=2;
        n_a=[n_asset,n_occ];
        n_z=[n_eta,n_theta];
        asset_grid=60*(linspace(0,1,n_asset).^3)';
        a_grid=[asset_grid; [0;1]];
        [eta_grid,pi_eta]=discretizeAR1_FarmerToda(0,0.9,0.2,n_eta);
        eta_grid=exp(eta_grid);
        theta_grid=[0; 1.5];
        pi_theta=[0.9,0.1; 0.2,0.8];
        z_grid=[eta_grid; theta_grid];
        pi_z=kron(pi_theta,pi_eta);
        n_d=0; d_grid=[];
        Params_c.r=0.04; Params_c.w=1.4;
        Params_c.alpha=0.36; Params_c.delta=0.06; Params_c.upsilon=0.88;
        Params_c.upsilon1=Params_c.alpha*Params_c.upsilon;
        Params_c.upsilon2=(1-Params_c.alpha)*Params_c.upsilon;
        Params_c.leverage=0.5;
        Params_c.phi=0.05;
        vfoptions.verbose_advice=0;
        ReturnFn=@(aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi) ...
            ReturnFn_nod_z_noe_nosemiz_with2A(aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi);
        FnsToEvaluate.assets=@(aprime,eprime,a,e,eta,theta) a;
        FnsToEvaluate.entrepreneur=@(aprime,eprime,a,e,eta,theta) e;
        FnsToEvaluate.savings=@(aprime,eprime,a,e,eta,theta) aprime-a;
        CondlRest.entrepreneurs=@(aprime,eprime,a,e,eta,theta) (e==1);
        CondlRest.positiveassets=@(aprime,eprime,a,e,eta,theta) (a>0);
        CondlRest.able=@(aprime,eprime,a,e,eta,theta) (theta>0);
        CondlRest.richworkers=@(aprime,eprime,a,e,eta,theta,w) (e==0)*(a>w); % a restriction with parameters
        CondlRest.always=@(aprime,eprime,a,e,eta,theta) 1;
        NeverRest.never=@(aprime,eprime,a,e,eta,theta) 0;
    end
    fprintf('\n-- CondlRest InfHorz, model %s --\n',configname)
    FnNames=fieldnames(FnsToEvaluate);
    nFns=length(FnNames);
    RestrNames=fieldnames(CondlRest);
    nRestr=length(RestrNames);
    npages=1+nRestr;
    PageNames=[{'unrestricted'};RestrNames];

    %% Solve
    [~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_c,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params_c,[]);

    %% Brute-force inputs: SD (N_ss,1), V (N_ss,nFns), R (N_ss,nRestr)
    VOG=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
    ROG=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy,CondlRest,Params_c,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
    SD=reshape(gather(double(StationaryDist)),[],1);
    N_ss=numel(SD);
    V=zeros(N_ss,nFns);
    for ff=1:nFns
        V(:,ff)=reshape(gather(double(VOG.(FnNames{ff}))),[N_ss,1]);
    end
    R=zeros(N_ss,nRestr);
    for rr=1:nRestr
        R(:,rr)=double(reshape(gather(ROG.(RestrNames{rr})),[N_ss,1])~=0);
    end

    simoptions_R=simoptions; simoptions_R.conditionalrestrictions=CondlRest;
    CS=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_R);
    AS=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_R);

    %% (1)-(2) the cross-section command, page by page
    for pp=1:npages
        if pp==1
            P=CS; Q=AS; w=SD;
        else
            P=CS.(RestrNames{pp-1}); Q=AS.(RestrNames{pp-1}); w=SD.*R(:,pp-1);
        end
        mass=sum(w);
        if mass==0 % a restriction with no mass: no function fields, NaN matrices, RestrictedSampleMass zero
            bad=sum(isfield(P,FnNames))+sum(~isnan(gather(P.CovarianceMatrix(:))))+sum(~isnan(gather(P.CorrelationMatrix(:))))+abs(gather(P.RestrictedSampleMass));
            fprintf('CondlRest InfHorz CrossSection %s, page %s has no mass: fn fields present + matrix entries not NaN + RestrictedSampleMass, this should be zero: %.3e \n',configname,PageNames{pp},bad)
            continue
        end
        w=w/mass;
        Bm=V'*w;
        dev=V-Bm';
        Bs=sqrt((dev.^2)'*w);
        BC=dev'*(dev.*w);
        BR=BC./(Bs*Bs'); BR(logical(eye(nFns)))=1;
        sok=(Bs>1e-8*max(1,abs(Bm)));
        ok=(sok*sok')>0; ok(logical(eye(nFns)))=true;
        dm=0; nanmis=0; das=0; dnamed=0;
        for ff=1:nFns
            a=[gather(P.(FnNames{ff}).Mean);gather(P.(FnNames{ff}).StdDeviation)]; b=[Bm(ff);Bs(ff)];
            nanmis=nanmis+sum(isnan(a)~=isnan(b)); d=abs(a-b); dm=max([dm;d(~isnan(d))]);
            das=max([das,abs(gather(P.(FnNames{ff}).Mean)-gather(Q.(FnNames{ff}).Mean)),abs(gather(P.(FnNames{ff}).StdDeviation)-gather(Q.(FnNames{ff}).StdDeviation))]);
        end
        Mt=gather(P.CovarianceMatrix); Rt=gather(P.CorrelationMatrix);
        d=abs(Mt(:)-BC(:)); nanmis=nanmis+sum(isnan(Mt(:))~=isnan(BC(:))); dm=max([dm;d(~isnan(d))]);
        d=abs(Rt(ok)-BR(ok)); nanmis=nanmis+sum(isnan(Rt(ok))); dc=max([0;d(~isnan(d))]);
        for ff1=1:nFns
            dnamed=max(dnamed,abs(gather(P.(FnNames{ff1}).(FnNames{ff1}))-1));
            for ff2=1:nFns
                if ff1~=ff2
                    a=gather(P.(FnNames{ff1}).CovarianceWith.(FnNames{ff2})); b=Mt(ff1,ff2);
                    dnamed=dnamed+(isnan(a)~=isnan(b)); if ~isnan(a) && ~isnan(b), dnamed=max(dnamed,abs(a-b)); end
                    a=gather(P.(FnNames{ff1}).CorrelationWith.(FnNames{ff2})); b=Rt(ff1,ff2);
                    dnamed=dnamed+(isnan(a)~=isnan(b)); if ~isnan(a) && ~isnan(b), dnamed=max(dnamed,abs(a-b)); end
                end
            end
        end
        fprintf('CondlRest InfHorz CrossSection %s, page %s: Mean/StdDeviation/CovarianceMatrix vs brute force, this should be zero: %.3e, NaN mismatches %i \n',configname,PageNames{pp},dm,nanmis)
        fprintf('CondlRest InfHorz CrossSection %s, page %s: CorrelationMatrix vs brute force (where both std devs are nonzero, %i of %i entries), this should be zero: %.3e \n',configname,PageNames{pp},sum(ok(:)),nFns^2,dc)
        fprintf('CondlRest InfHorz CrossSection %s, page %s: named fields vs matrices, this should be zero: %.3e \n',configname,PageNames{pp},dnamed)
        fprintf('CondlRest InfHorz CrossSection %s, page %s: Mean/StdDeviation vs AllStats, this should be zero: %.3e \n',configname,PageNames{pp},das)
        if pp>1
            fprintf('CondlRest InfHorz CrossSection %s, page %s: RestrictedSampleMass vs brute force, this should be zero: %.3e; AllStats RestrictedSampleMass vs brute force, this should be zero: %.3e \n',configname,PageNames{pp},abs(gather(P.RestrictedSampleMass)-mass),abs(gather(Q.RestrictedSampleMass)-mass))
        end
    end

    %% (3) every statistic of AllStats, every page, vs brute force
    maxd=zeros(1,7); mism=zeros(1,7); worst=repmat({''},1,7);
    for pp=1:npages
        if pp==1
            OutP=AS; w=SD;
        else
            OutP=AS.(RestrNames{pp-1}); w=SD.*R(:,pp-1);
        end
        if sum(w)==0 % a restriction with no mass: AllStats reports only its RestrictedSampleMass
            mism(1)=mism(1)+sum(isfield(OutP,FnNames));
            continue
        end
        for ff=1:nFns
            v=V(:,ff);
            keep=(w>0); wk=w(keep); vk=v(keep);
            wk=wk/sum(wk);
            [u,~,ic]=unique(vk);
            m=accumarray(ic,wk);
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
            for ss=1:numel(statnames)
                gg=statgroup(ss);
                X=OutP.(FnNames{ff});
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
                a=gather(X(:)); b=B{ss}(:);
                if numel(a)~=numel(b)
                    mism(gg)=mism(gg)+1;
                    continue
                end
                mism(gg)=mism(gg)+sum(isnan(a)~=isnan(b))+sum(xor(isinf(a),isinf(b)));
                d=abs(a-b); d=d(~isnan(d));
                if ~isempty(d) && max(d)>maxd(gg)
                    maxd(gg)=max(d);
                    worst{gg}=sprintf('%s, page %s, %s',statnames{ss},PageNames{pp},FnNames{ff});
                end
            end
        end
    end
    for gg=1:7
        fprintf('CondlRest InfHorz AllStats %s, %s vs brute force on every page, this should be zero: %.3e, NaN/Inf mismatches + missing fields %i (worst: %s) \n',configname,groupnames{gg},maxd(gg),mism(gg),worst{gg})
    end

    %% (4) the always page equals the unrestricted output, every leaf; the unrestricted cross-section output equals a run without restrictions
    CS0=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
    for cmd=1:3
        if cmd==1
            cmdname='CrossSection: page always vs a run without restrictions'; A=CS0; B=CS.always;
        elseif cmd==2
            cmdname='CrossSection: unrestricted output with restrictions vs a run without restrictions'; A=CS0; B=CS;
        else
            cmdname='AllStats: page always vs the unrestricted output'; A=struct(); B=struct();
            for ff=1:nFns
                A.(FnNames{ff})=AS.(FnNames{ff}); B.(FnNames{ff})=AS.always.(FnNames{ff});
            end
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
        if cmd~=2 % (with restrictions the output also has the restriction subtrees)
            mis=mis+sum(~ismember(pb,[pa,{'.RestrictedSampleMass'}]));
        end
        fprintf('CondlRest InfHorz %s %s, every leaf, this should be zero: %.3e, mismatches %i \n',cmdname,configname,maxd,mis)
    end

    %% (5) a restriction of zero mass in the cross-section command
    simoptions_N=simoptions; simoptions_N.conditionalrestrictions=NeverRest;
    lastwarn('');
    CSN=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_N);
    warned=~isempty(lastwarn);
    notnan=sum(~isnan(gather(CSN.never.CovarianceMatrix(:))))+sum(~isnan(gather(CSN.never.CorrelationMatrix(:))))+sum(isfield(CSN.never,FnNames));
    fprintf('CondlRest InfHorz CrossSection %s: zero-mass restriction: no warning + RestrictedSampleMass + entries not NaN + fn fields present, this should be zero: %.3e \n',configname,(1-warned)+abs(gather(CSN.never.RestrictedSampleMass))+notnan)

    output.(configname).CrossSection=CS;
    output.(configname).AllStats=AS;
end

end
