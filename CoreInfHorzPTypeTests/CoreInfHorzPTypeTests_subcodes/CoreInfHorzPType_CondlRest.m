function output=CoreInfHorzPType_CondlRest(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% Test (CoreInfHorzPTypeTests part 9): simoptions.conditionalrestrictions in EvalFnOnAgentDist_AllStats_InfHorz_PType, per type and grouped, against a brute force
% (the InfHorz PType stats commands are AllStats and AggVars; AggVars takes no restrictions, and there is no InfHorz PType
% CrossSection or AutoCorr command). On two models:
%   d_z    with d and z, two types (ptweights 0.6/0.4)
%   nod_z  without d, with z, three types, the middle one of zero weight (ptweights 0.5/0/0.5): it is not in the grouped
%          population, but its per-type (within-type) stats are still reported
% The types differ by a per-type wage w_pt, which enters the return function, a FnsToEvaluate (earnings) and a restriction
% (highincome), so the per-type parameters must reach the functions the command evaluates for each type.
%
% Checks:
%   (1) every statistic (Mean, RatioMeanToMedian, Median, percentiles, Variance, StdDeviation, Gini, LorenzCurve, Minimum, Maximum,
%       QuantileCutoffs, QuantileMeans, top/bottom shares), per type and grouped, on every page (unrestricted and each restriction)
%       vs a brute force written from the definitions, as in CoreFHorzTests part 12 and CoreFHorzPTypeTests part 17: every
%       type's values from the single-type ValuesOnGrid on that type's Policy and parameters, its weights from its
%       StationaryDist (mass one within the type), the grouped population every type's (restricted) weights times its ptweight
%   (2) the RestrictedSampleMass, per type (a within-type share) and TotalAllPTypes (a population share)
%   (3) each type's output equals the single-type AllStats on that type, under the same restrictions, every leaf
%   (4) the always page equals the unrestricted output, every leaf, per type and grouped

npts=100; nq=20; tol=10^(-12); % the defaults of AllStats

statnames={'Mean','RatioMeanToMedian','Median','MoreInequality.Percentile50th','MoreInequality.Percentile90th','MoreInequality.Percentile95th','MoreInequality.Percentile99th', ...
    'Variance','StdDeviation','Gini','LorenzCurve','Minimum','Maximum','QuantileCutoffs','QuantileMeans', ...
    'MoreInequality.Top1share','MoreInequality.Top5share','MoreInequality.Top10share','MoreInequality.Bottom50share'};
statgroup=[1,1,2,2,2,2,2,3,3,4,4,5,5,6,6,7,7,7,7];
groupnames={'Mean/RatioMeanToMedian','Median/percentiles','Variance/StdDeviation','Gini/LorenzCurve','Minimum/Maximum','QuantileCutoffs/QuantileMeans','Top/Bottom shares'};

PTypeDistParamNames={'ptypemass'};
n_d_withd=n_d; d_grid_withd=d_grid;

output=struct();

for cc=1:2
    Params_PT=Params;
    FnsToEvaluate=struct();
    CondlRest=struct();
    vfoptions=struct();
    simoptions=struct();
    if cc==1
        configname='d_z';
        Names_i={'low','high'};
        ptw=[0.6;0.4];
        Params_PT.w_pt=[1;1.4];
        n_d=n_d_withd; d_grid=d_grid_withd;
        ReturnFn=@(d,aprime,a,z,r,w_pt,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w_pt,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a,z) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,w_pt) w_pt*d*z;
        FnsToEvaluate.savings=@(d,aprime,a,z) aprime-a; % negative values: Gini, Lorenz curve and shares are NaN
        CondlRest.positiveassets=@(d,aprime,a,z) (a>0);
        CondlRest.highz=@(d,aprime,a,z) (z>1);
        CondlRest.highincome=@(d,aprime,a,z,r,w_pt) (r*a+w_pt*d*z>0.5); % a per-type parameter
        CondlRest.always=@(d,aprime,a,z) 1;
    else
        configname='nod_z';
        Names_i={'low','ghost','high'};
        ptw=[0.5;0;0.5];
        Params_PT.w_pt=[1;1.2;1.4];
        n_d=0; d_grid=[];
        ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
        FnsToEvaluate.consumption=@(aprime,a,z,r,w_pt) (1+r)*a+w_pt*z-aprime;
        CondlRest.positiveassets=@(aprime,a,z) (a>0);
        CondlRest.highz=@(aprime,a,z) (z>1);
        CondlRest.highincome=@(aprime,a,z,w_pt) (w_pt*z>1.1); % a per-type parameter
        CondlRest.always=@(aprime,a,z) 1;
    end
    fprintf('\n-- CondlRest InfHorz PType, model %s --\n',configname)
    N_i=length(Names_i);
    Params_PT.(PTypeDistParamNames{1})=ptw;
    FnNames=fieldnames(FnsToEvaluate);
    nFns=length(FnNames);
    RestrNames=fieldnames(CondlRest);
    nRestr=length(RestrNames);
    npages=1+nRestr;
    PageNames=[{'unrestricted'};RestrNames];

    %% Solve
    [~,Policy_PT]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_PT,DiscountFactorParamNames,vfoptions);
    StationaryDist_PT=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,Names_i,pi_z,Params_PT,simoptions);
    simoptions_R=simoptions; simoptions_R.conditionalrestrictions=CondlRest;
    AS=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params_PT,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_R);

    %% Brute-force inputs, type by type
    SD=cell(N_i,1); V=cell(N_i,1); R=cell(N_i,1); AS_solo=cell(N_i,1);
    for ii=1:N_i
        nm=Names_i{ii};
        Params_ii=Params_PT;
        Params_ii.w_pt=Params_PT.w_pt(ii);
        VOG=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy_PT.(nm),FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
        ROG=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy_PT.(nm),CondlRest,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
        SD{ii}=reshape(gather(double(StationaryDist_PT.(nm))),[],1);
        N_ss=numel(SD{ii});
        V{ii}=zeros(N_ss,nFns);
        for ff=1:nFns
            V{ii}(:,ff)=reshape(gather(double(VOG.(FnNames{ff}))),[N_ss,1]);
        end
        R{ii}=zeros(N_ss,nRestr);
        for rr=1:nRestr
            R{ii}(:,rr)=double(reshape(gather(ROG.(RestrNames{rr})),[N_ss,1])~=0);
        end
        % (3) the single-type command on this type
        AS_solo{ii}=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist_PT.(nm),Policy_PT.(nm),FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_R);
    end

    %% (1) every statistic, per type and grouped, every page, vs brute force
    for tt=1:N_i+1 % the types, then grouped
        if tt<=N_i
            unitname=['type ',Names_i{tt}];
        else
            unitname='grouped';
        end
        maxd=zeros(1,7); mism=zeros(1,7); worst=repmat({''},1,7);
        for pp=1:npages
            if pp==1
                OutP=AS;
            else
                OutP=AS.(RestrNames{pp-1});
            end
            for ff=1:nFns
                v=zeros(0,1); w=zeros(0,1);
                for ii=1:N_i
                    if tt<=N_i && ii~=tt
                        continue
                    end
                    if pp==1
                        Wii=SD{ii};
                    else
                        Wii=SD{ii}.*R{ii}(:,pp-1);
                    end
                    if tt>N_i
                        Wii=ptw(ii)*Wii;
                    end
                    v=[v; V{ii}(:,ff)]; %#ok<AGROW>
                    w=[w; Wii]; %#ok<AGROW>
                end
                keep=(w>0); w=w(keep); v=v(keep);
                if isempty(w)
                    continue % (no restricted mass: nothing to compare; the restrictions here all have mass)
                end
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
            fprintf('CondlRest InfHorz PType AllStats %s, %s, %s vs brute force on every page, this should be zero: %.3e, NaN/Inf mismatches + missing fields %i (worst: %s) \n',configname,unitname,groupnames{gg},maxd(gg),mism(gg),worst{gg})
        end
    end

    %% (2) the restricted sample masses
    dR=0;
    for rr=1:nRestr
        RSMtype=zeros(N_i,1);
        for ii=1:N_i
            RSMtype(ii)=sum(SD{ii}.*R{ii}(:,rr));
            dR=max(dR,abs(gather(AS.(RestrNames{rr}).RestrictedSampleMass.(Names_i{ii}))-RSMtype(ii)));
        end
        dR=max(dR,abs(gather(AS.(RestrNames{rr}).RestrictedSampleMass.TotalAllPTypes)-sum(ptw.*RSMtype)));
    end
    fprintf('CondlRest InfHorz PType AllStats %s: RestrictedSampleMass per type and TotalAllPTypes vs brute force, this should be zero: %.3e \n',configname,dR)

    %% (3) each type's output equals the single-type command on that type, every page, every leaf
    maxd=0; mis=0;
    for ii=1:N_i
        for pp=1:npages
            for ff=1:nFns
                if pp==1
                    [p0,v0]=CoreStructLeaves(AS_solo{ii}.(FnNames{ff})); [p1,v1]=CoreStructLeaves(AS.(FnNames{ff}).(Names_i{ii}));
                else
                    [p0,v0]=CoreStructLeaves(AS_solo{ii}.(RestrNames{pp-1}).(FnNames{ff})); [p1,v1]=CoreStructLeaves(AS.(RestrNames{pp-1}).(FnNames{ff}).(Names_i{ii}));
                end
                for k=1:numel(p0)
                    [present,idx]=ismember(p0{k},p1);
                    if present && isequal(size(v0{k}),size(v1{idx}))
                        a=v0{k}; b=v1{idx}; mis=mis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                    else
                        mis=mis+1;
                    end
                end
                mis=mis+sum(~ismember(p1,p0));
            end
        end
    end
    fprintf('CondlRest InfHorz PType AllStats %s: per type vs the single-type command on that type, every page, every leaf, this should be zero: %.3e, mismatches %i \n',configname,maxd,mis)

    %% (4) the always page equals the unrestricted output, every leaf (per type and grouped)
    maxd=0; mis=0;
    for ff=1:nFns
        [p0,v0]=CoreStructLeaves(AS.(FnNames{ff})); [pA,vA]=CoreStructLeaves(AS.always.(FnNames{ff}));
        for k=1:numel(p0)
            [present,idx]=ismember(p0{k},pA);
            if present && isequal(size(v0{k}),size(vA{idx}))
                a=v0{k}; b=vA{idx}; mis=mis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
            else
                mis=mis+1;
            end
        end
        mis=mis+sum(~ismember(pA,p0));
    end
    fprintf('CondlRest InfHorz PType AllStats %s: page always vs the unrestricted output, every leaf, this should be zero: %.3e, mismatches %i \n',configname,maxd,mis)

    output.(configname).AllStats=AS;
end

end
