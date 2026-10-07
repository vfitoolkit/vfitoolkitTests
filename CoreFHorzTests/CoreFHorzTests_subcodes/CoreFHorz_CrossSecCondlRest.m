function output=CoreFHorz_CrossSecCondlRest(n_d,n_d2_semiz,n_a,n_z,N_j,d_grid,d2_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Test: simoptions.conditionalrestrictions in the FHorz stats commands (doPart(12) holds all the conditional-restriction tests of
% CoreFHorzTests: this, and the AutoCorrTransProbs brute forces TestFnsToEvaluate_AutoCorrBruteForce(_semiz); the restriction pages
% of whichcombos are in doPart(13)):
%   EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz                       output .(restriction).(fn).Mean etc, .CovarianceMatrix, .RestrictedSampleMass
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz   the same per age group, RestrictedSampleMass 1-by-N_j
%   EvalFnOnAgentDist_AllStats_FHorz_Case1                              every statistic
%   LifeCycleProfiles_FHorz_Case1                                       every statistic, per age (the default), in age groups of five, in one age group
% on three configs: with d, z and e; without d, z or e (the N_z==0 branch); without d1, with z and semiz.
%
% Pages are the unrestricted output first, then the restrictions in the fieldnames order of simoptions.conditionalrestrictions:
%   positiveassets (a>0)       restricted mass differs by age, and is zero at age 1 (everyone starts with zero assets)
%   working        (agej<Jr)   a restriction with parameters; the last age group is all retired, so has zero restricted mass
%   atJ1           (agej==J1)  mass at one age only: in the age-conditional versions every age group but one has zero restricted mass
%   always         (1)         must reproduce the unrestricted output exactly
%
% Checks, CrossSectionCovarCorr and its AgeConditional version:
%   (1) every restricted (and the unrestricted) Mean, StdDeviation, CovarianceMatrix, CorrelationMatrix and RestrictedSampleMass vs a
%       brute-force calculation from ValuesOnGrid and the agent distribution (the restricted distribution normalised to mass one,
%       within each age group for the age-conditional version). Correlations are compared only where both brute-force std devs are
%       nonzero: a function that is constant on the restricted sample has a std dev that is rounding noise, and so a correlation
%       that is noise divided by noise.
%   (2) the named fields (.CovarianceWith/.CorrelationWith/self-correlation) equal the matrices.
%   (3) restricted Mean/StdDeviation vs AllStats (pooled) and LifeCycleProfiles (age groups) under the same restrictions.
%   (4) the 'always' page equals the unrestricted page, every leaf; the unrestricted output with restrictions equals a run without.
%   (5) a restriction with zero mass: warning, RestrictedSampleMass zero, everything NaN.
%   (6) simoptions.whichcombos with a page dimension: selected entries equal the full run, unselected are NaN; a whichcombos without
%       the page dimension equals its repmat over the pages; a wrong number of pages errors.
% Checks, AllStats and LifeCycleProfiles:
%   (7) every statistic (Mean, RatioMeanToMedian, Median, percentiles, Variance, StdDeviation, Gini, LorenzCurve, Minimum, Maximum,
%       QuantileCutoffs, QuantileMeans, top/bottom shares) on every page vs a brute force written from the definitions (see the
%       block), for AllStats and for LifeCycleProfiles per age, in age groups of five and in one age group of every age.
%   (8) the 'always' page equals the unrestricted output, every leaf.
%   (9) the atJ1 page of AllStats equals the unrestricted LifeCycleProfiles at age J1, every leaf.
%   (10) LifeCycleProfiles in one age group of every age equals AllStats, every page, every leaf (the age group weights each age by
%       its restricted mass).

Params.J1=floor(N_j/2);
agegroupings=1:5:N_j;
ngroups=length(agegroupings);

n_d_withd=n_d; d_grid_withd=d_grid;
n_z_withz=n_z; z_grid_withz=z_grid; pi_z_withz=pi_z;

for cc=1:3
    %% Setup of this config
    FnsToEvaluate=struct();
    CondlRest=struct();
    NeverRest=struct();
    if cc==1
        configname='d_z_e';
        n_d=n_d_withd; d_grid=d_grid_withd;
        n_z=n_z_withz; z_grid=z_grid_withz; pi_z=pi_z_withz;
        vfoptions=struct();
        vfoptions.n_e=vfoptionsbaseline.n_e;
        vfoptions.e_grid=vfoptionsbaseline.e_grid;
        vfoptions.pi_e=vfoptionsbaseline.pi_e;
        simoptions=struct();
        simoptions.n_e=simoptionsbaseline.n_e;
        simoptions.e_grid=simoptionsbaseline.e_grid;
        simoptions.pi_e=simoptionsbaseline.pi_e;
        jequaloneDist=zeros(n_a,n_z,vfoptions.n_e,'gpuArray');
        jequaloneDist(1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;
        ReturnFn=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
        FnsToEvaluate.assets=@(d,aprime,a,z,e) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,e,w,kappa_j) w*kappa_j*z*e*d;
        FnsToEvaluate.consumption=@(d,aprime,a,z,e,r,w,kappa_j) (1+r)*a+w*kappa_j*z*e*d-aprime;
        CondlRest.positiveassets=@(d,aprime,a,z,e) (a>0);
        CondlRest.working=@(d,aprime,a,z,e,agej,Jr) (agej<Jr);
        CondlRest.atJ1=@(d,aprime,a,z,e,agej,J1) (agej==J1);
        CondlRest.always=@(d,aprime,a,z,e) 1;
        NeverRest.never=@(d,aprime,a,z,e) 0;
    elseif cc==2
        configname='nod_noz_noe';
        n_d=0; d_grid=[];
        n_z=0; z_grid=[]; pi_z=[];
        vfoptions=struct();
        simoptions=struct();
        jequaloneDist=zeros(n_a,1,'gpuArray');
        jequaloneDist(1)=1;
        ReturnFn=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate.assets=@(aprime,a) a;
        FnsToEvaluate.earnings=@(aprime,a,w,kappa_j) w*kappa_j;
        FnsToEvaluate.savings=@(aprime,a) aprime-a;
        CondlRest.positiveassets=@(aprime,a) (a>0);
        CondlRest.working=@(aprime,a,agej,Jr) (agej<Jr);
        CondlRest.atJ1=@(aprime,a,agej,J1) (agej==J1);
        CondlRest.always=@(aprime,a) 1;
        NeverRest.never=@(aprime,a) 0;
    else
        configname='nod1_z_noe_semiz';
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
        ReturnFn=@(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,z) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,w,kappa_j) w*kappa_j*z*semiz;
        FnsToEvaluate.consumption=@(d2,aprime,a,semiz,z,r,w,kappa_j) (1+r)*a+w*kappa_j*z*semiz-aprime;
        CondlRest.positiveassets=@(d2,aprime,a,semiz,z) (a>0);
        CondlRest.working=@(d2,aprime,a,semiz,z,agej,Jr) (agej<Jr);
        CondlRest.atJ1=@(d2,aprime,a,semiz,z,agej,J1) (agej==J1);
        CondlRest.always=@(d2,aprime,a,semiz,z) 1;
        NeverRest.never=@(d2,aprime,a,semiz,z) 0;
    end
    fprintf('\n-- CondlRest CrossSection, config %s --\n',configname)
    FnNames=fieldnames(FnsToEvaluate);
    nFns=length(FnNames);
    RestrNames=fieldnames(CondlRest);
    nRestr=length(RestrNames);
    npages=1+nRestr;
    PageNames=[{'unrestricted'};RestrNames];

    %% Solve
    [~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

    %% Brute force, from ValuesOnGrid and the agent distribution
    VOG=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    ROG=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy,CondlRest,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    SD=reshape(gather(double(StationaryDist)),[],N_j);
    N_ss=size(SD,1);
    V=zeros(N_ss,N_j,nFns);
    for ff=1:nFns
        V(:,:,ff)=reshape(gather(double(VOG.(FnNames{ff}))),[N_ss,N_j]);
    end
    R=zeros(N_ss,N_j,nRestr);
    for rr=1:nRestr
        R(:,:,rr)=reshape(gather(double(ROG.(RestrNames{rr}))),[N_ss,N_j]);
    end
    % pooled (kk=0) and per age group (kk=1:ngroups); a page/group with no mass stays NaN
    BMean=nan(nFns,npages); BStd=nan(nFns,npages); BCov=nan(nFns,nFns,npages); BCorr=nan(nFns,nFns,npages); BCorrOK=false(nFns,nFns,npages);
    AMean=nan(nFns,ngroups,npages); AStd=nan(nFns,ngroups,npages); ACov=nan(nFns,nFns,ngroups,npages); ACorr=nan(nFns,nFns,ngroups,npages); ACorrOK=false(nFns,nFns,ngroups,npages);
    BMass=nan(1,npages); AMassJ=nan(npages,N_j);
    for pp=1:npages
        if pp==1
            W=SD;
        else
            W=SD.*R(:,:,pp-1);
        end
        BMass(pp)=sum(W(:));
        AMassJ(pp,:)=sum(W,1);
        for kk=0:ngroups
            if kk==0
                jj=1:N_j;
            else
                j1=agegroupings(kk);
                if kk<ngroups
                    jend=agegroupings(kk+1)-1;
                else
                    jend=N_j;
                end
                jj=j1:jend;
            end
            w=reshape(W(:,jj),[],1);
            mass=sum(w);
            if mass==0
                continue
            end
            w=w/mass;
            m=zeros(nFns,1); s=zeros(nFns,1); dev=zeros(numel(w),nFns);
            for ff=1:nFns
                v=reshape(V(:,jj,ff),[],1);
                m(ff)=sum(w.*v);
                dev(:,ff)=v-m(ff);
                s(ff)=sqrt(sum(w.*dev(:,ff).^2));
            end
            C=dev'*(dev.*w);
            Cr=C./(s*s');
            Cr(logical(eye(nFns)))=1;
            sok=(s>1e-8*max(1,abs(m))); % std dev is not rounding noise
            ok=(sok*sok')>0;
            ok(logical(eye(nFns)))=true;
            if kk==0
                BMean(:,pp)=m; BStd(:,pp)=s; BCov(:,:,pp)=C; BCorr(:,:,pp)=Cr; BCorrOK(:,:,pp)=ok;
            else
                AMean(:,kk,pp)=m; AStd(:,kk,pp)=s; ACov(:,:,kk,pp)=C; ACorr(:,:,kk,pp)=Cr; ACorrOK(:,:,kk,pp)=ok;
            end
        end
    end

    simoptions_R=simoptions;
    simoptions_R.conditionalrestrictions=CondlRest;
    simoptions_RA=simoptions_R;
    simoptions_RA.agegroupings=agegroupings;

    %% ===================== EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz =====================
    CS=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_R);
    AS=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_R);

    % (1)-(3) per page
    for pp=1:npages
        if pp==1
            P=CS; Q=AS;
        else
            P=CS.(RestrNames{pp-1}); Q=AS.(RestrNames{pp-1});
        end
        % (1) vs brute force
        dm=0; nanmis=0; dnamed=0; das=0;
        for ff=1:nFns
            dm=max([dm,abs(gather(P.(FnNames{ff}).Mean)-BMean(ff,pp)),abs(gather(P.(FnNames{ff}).StdDeviation)-BStd(ff,pp))]);
            nanmis=nanmis+isnan(gather(P.(FnNames{ff}).Mean))+isnan(gather(P.(FnNames{ff}).StdDeviation));
            % (3) vs AllStats
            das=max([das,abs(gather(P.(FnNames{ff}).Mean)-gather(Q.(FnNames{ff}).Mean)),abs(gather(P.(FnNames{ff}).StdDeviation)-gather(Q.(FnNames{ff}).StdDeviation))]);
        end
        Mt=gather(P.CovarianceMatrix); Rt=gather(P.CorrelationMatrix);
        d=abs(Mt(:)-reshape(BCov(:,:,pp),[],1));
        nanmis=nanmis+sum(isnan(Mt(:))~=isnan(reshape(BCov(:,:,pp),[],1)));
        dm=max([dm;d(~isnan(d))]);
        okpp=BCorrOK(:,:,pp);
        Bc=BCorr(:,:,pp);
        dc=abs(Rt(okpp)-Bc(okpp));
        nanmis=nanmis+sum(isnan(Rt(okpp)));
        dc=max([0;dc(~isnan(dc))]);
        % (2) named fields vs the matrices
        for ff1=1:nFns
            dnamed=max(dnamed,abs(gather(P.(FnNames{ff1}).(FnNames{ff1}))-1));
            for ff2=1:nFns
                if ff1~=ff2
                    a=gather(P.(FnNames{ff1}).CovarianceWith.(FnNames{ff2})); b=Mt(ff1,ff2);
                    if isnan(a)~=isnan(b)
                        dnamed=dnamed+1;
                    elseif ~isnan(a)
                        dnamed=max(dnamed,abs(a-b));
                    end
                    a=gather(P.(FnNames{ff1}).CorrelationWith.(FnNames{ff2})); b=Rt(ff1,ff2);
                    if isnan(a)~=isnan(b)
                        dnamed=dnamed+1;
                    elseif ~isnan(a)
                        dnamed=max(dnamed,abs(a-b));
                    end
                end
            end
        end
        fprintf('CondlRest CrossSection %s, page %s: Mean/StdDeviation/CovarianceMatrix vs brute force, this should be zero: %.3e, NaN mismatches %i \n',configname,PageNames{pp},dm,nanmis)
        fprintf('CondlRest CrossSection %s, page %s: CorrelationMatrix vs brute force (where both std devs are nonzero, %i of %i entries), this should be zero: %.3e \n',configname,PageNames{pp},sum(okpp(:)),nFns^2,dc)
        fprintf('CondlRest CrossSection %s, page %s: named fields vs matrices, this should be zero: %.3e \n',configname,PageNames{pp},dnamed)
        fprintf('CondlRest CrossSection %s, page %s: Mean/StdDeviation vs AllStats, this should be zero: %.3e \n',configname,PageNames{pp},das)
        if pp>1
            fprintf('CondlRest CrossSection %s, page %s: RestrictedSampleMass vs brute force, this should be zero: %.3e \n',configname,PageNames{pp},abs(gather(P.RestrictedSampleMass)-BMass(pp)))
        end
    end

    % (4) the 'always' page equals the unrestricted page; the unrestricted output equals a run without restrictions
    CS0=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    [p0,v0]=CoreStructLeaves(CS0);
    [pA,vA]=CoreStructLeaves(CS.always);
    [pU,vU]=CoreStructLeaves(CS);
    maxd=0; mis=0; maxd2=0; mis2=0;
    for k=1:numel(p0)
        [present,idx]=ismember(p0{k},pA);
        if present && isequal(size(v0{k}),size(vA{idx}))
            a=v0{k}; b=vA{idx}; mis=mis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        else
            mis=mis+1;
        end
        [present,idx]=ismember(p0{k},pU);
        if present && isequal(size(v0{k}),size(vU{idx}))
            a=v0{k}; b=vU{idx}; mis2=mis2+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd2=max([maxd2;d]);
        else
            mis2=mis2+1;
        end
    end
    mis=mis+sum(~ismember(pA,[p0,{'.RestrictedSampleMass'}]));
    fprintf('CondlRest CrossSection %s: page always vs a run without restrictions, every leaf, this should be zero: %.3e, mismatches %i \n',configname,maxd,mis)
    fprintf('CondlRest CrossSection %s: unrestricted output with restrictions vs a run without restrictions, every leaf, this should be zero: %.3e, mismatches %i \n',configname,maxd2,mis2)

    % (5) a restriction with zero mass
    simoptions_N=simoptions; simoptions_N.conditionalrestrictions=NeverRest;
    lastwarn('');
    CSN=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_N);
    warned=~isempty(lastwarn);
    notnan=sum(~isnan(gather(CSN.never.CovarianceMatrix(:))))+sum(~isnan(gather(CSN.never.CorrelationMatrix(:))))+sum(isfield(CSN.never,FnNames));
    fprintf('CondlRest CrossSection %s: zero-mass restriction: no warning + RestrictedSampleMass + entries not NaN + fn fields present, this should be zero: %.3e \n',configname,(1-warned)+abs(gather(CSN.never.RestrictedSampleMass))+notnan)

    % (6) whichcombos with a page dimension
    wc=zeros(nFns,nFns,npages);
    for pp=1:npages
        for ff1=1:nFns
            for ff2=1:nFns
                if ff1~=ff2
                    wc(ff1,ff2,pp)=mod(ff1+ff2+pp,2);
                end
            end
        end
    end
    wc(1,1,1)=1; % the unrestricted page: assets' own stats by the diagonal alone
    simoptions_W=simoptions_R; simoptions_W.whichcombos=wc;
    CSW=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_W);
    maxd=0; bad=0;
    for pp=1:npages
        if pp==1
            P=CS; S=CSW;
        else
            P=CS.(RestrNames{pp-1}); S=CSW.(RestrNames{pp-1});
        end
        fnon=any(wc(:,:,pp),2);
        on=(max(wc(:,:,pp),diag(fnon))==1);
        Mf=gather(P.CovarianceMatrix); Ms=gather(S.CovarianceMatrix); Rf=gather(P.CorrelationMatrix); Rs=gather(S.CorrelationMatrix);
        d=[abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]; d=d(~isnan(d)); maxd=max([maxd;d]);
        bad=bad+sum(isnan(Mf(on))~=isnan(Ms(on)))+sum(isnan(Rf(on))~=isnan(Rs(on)))+sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on)));
        for ff=1:nFns
            bad=bad+(isfield(S,FnNames{ff})~=fnon(ff)); % a function with nothing selected on this page has no fields
            if fnon(ff) && isfield(S,FnNames{ff})
                maxd=max([maxd,abs(gather(S.(FnNames{ff}).Mean)-gather(P.(FnNames{ff}).Mean)),abs(gather(S.(FnNames{ff}).StdDeviation)-gather(P.(FnNames{ff}).StdDeviation))]);
            end
        end
    end
    fprintf('CondlRest CrossSection %s, whichcombos with a page dimension: selected entries vs full run, this should be zero: %.3e, unselected not NaN + NaN mismatches + field presence mismatches %i \n',configname,maxd,bad)
    wc2=[1,1,0; 0,0,1; 0,0,1];
    simoptions_W.whichcombos=wc2;
    A=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_W);
    simoptions_W.whichcombos=repmat(wc2,[1,1,npages]);
    B=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_W);
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
    fprintf('CondlRest CrossSection %s, a whichcombos without the page dimension vs its repmat over the pages, every leaf, this should be zero: %.3e, mismatches %i \n',configname,maxd,mis+sum(~ismember(pb,pa)))
    errored=0;
    simoptions_bad=simoptions_R; simoptions_bad.whichcombos=ones(nFns,nFns,npages-1);
    try
        EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_bad);
    catch
        errored=1;
    end
    fprintf('CondlRest CrossSection %s, a whichcombos with the wrong number of pages must error, this should be zero: %.3e \n',configname,1-errored)

    %% ===================== EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz =====================
    ACS=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_RA);
    LCP=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_RA);

    % (1)-(3) per page
    for pp=1:npages
        if pp==1
            P=ACS; Q=LCP;
        else
            P=ACS.(RestrNames{pp-1}); Q=LCP.(RestrNames{pp-1});
        end
        dm=0; nanmis=0; dnamed=0; dlcp=0; nanmislcp=0;
        for ff=1:nFns
            a=[gather(P.(FnNames{ff}).Mean);gather(P.(FnNames{ff}).StdDeviation)]; b=[AMean(ff,:,pp);AStd(ff,:,pp)];
            nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); dm=max([dm;d(~isnan(d))]);
            % (3) vs LifeCycleProfiles
            c=[gather(Q.(FnNames{ff}).Mean);gather(Q.(FnNames{ff}).StdDeviation)];
            nanmislcp=nanmislcp+sum(isnan(a(:))~=isnan(c(:))); d=abs(a(:)-c(:)); dlcp=max([dlcp;d(~isnan(d))]);
        end
        Mt=gather(P.CovarianceMatrix); Rt=gather(P.CorrelationMatrix);
        Bm=ACov(:,:,:,pp); Bc=ACorr(:,:,:,pp); okpp=ACorrOK(:,:,:,pp);
        nanmis=nanmis+sum(isnan(Mt(:))~=isnan(Bm(:)));
        d=abs(Mt(:)-Bm(:)); dm=max([dm;d(~isnan(d))]);
        nanmis=nanmis+sum(isnan(Rt(okpp)));
        dc=abs(Rt(okpp)-Bc(okpp)); dc=max([0;dc(~isnan(dc))]);
        % (2) named fields vs the matrices
        for ff1=1:nFns
            a=gather(P.(FnNames{ff1}).(FnNames{ff1})); b=reshape(Rt(ff1,ff1,:),[1,ngroups]);
            dnamed=dnamed+sum(isnan(a)~=isnan(b)); d=abs(a-b); dnamed=max([dnamed,d(~isnan(d))]);
            for ff2=1:nFns
                if ff1~=ff2
                    a=gather(P.(FnNames{ff1}).CovarianceWith.(FnNames{ff2})); b=reshape(Mt(ff1,ff2,:),[1,ngroups]);
                    dnamed=dnamed+sum(isnan(a)~=isnan(b)); d=abs(a-b); dnamed=max([dnamed,d(~isnan(d))]);
                    a=gather(P.(FnNames{ff1}).CorrelationWith.(FnNames{ff2})); b=reshape(Rt(ff1,ff2,:),[1,ngroups]);
                    dnamed=dnamed+sum(isnan(a)~=isnan(b)); d=abs(a-b); dnamed=max([dnamed,d(~isnan(d))]);
                end
            end
        end
        fprintf('CondlRest AgeConditional CrossSection %s, page %s: Mean/StdDeviation/CovarianceMatrix vs brute force, this should be zero: %.3e, NaN mismatches %i \n',configname,PageNames{pp},dm,nanmis)
        fprintf('CondlRest AgeConditional CrossSection %s, page %s: CorrelationMatrix vs brute force (where both std devs are nonzero, %i of %i entries), this should be zero: %.3e \n',configname,PageNames{pp},sum(okpp(:)),numel(okpp),dc)
        fprintf('CondlRest AgeConditional CrossSection %s, page %s: named fields vs matrices, this should be zero: %.3e \n',configname,PageNames{pp},dnamed)
        fprintf('CondlRest AgeConditional CrossSection %s, page %s: Mean/StdDeviation vs LifeCycleProfiles, this should be zero: %.3e, NaN mismatches %i \n',configname,PageNames{pp},dlcp,nanmislcp)
        if pp>1
            fprintf('CondlRest AgeConditional CrossSection %s, page %s: RestrictedSampleMass (1-by-N_j) vs brute force, this should be zero: %.3e, size mismatch %i \n',configname,PageNames{pp},max(abs(reshape(gather(P.RestrictedSampleMass),1,[])-AMassJ(pp,1:numel(P.RestrictedSampleMass)))),~isequal(size(P.RestrictedSampleMass),[1,N_j]))
        end
    end
    fprintf('CondlRest AgeConditional CrossSection %s: age groups with no restricted mass (atJ1 should have %i, working 1), information only: positiveassets %i, working %i, atJ1 %i \n',configname,ngroups-1,sum(isnan(AMean(1,:,2))),sum(isnan(AMean(1,:,3))),sum(isnan(AMean(1,:,4))))

    % (4) the 'always' page equals the unrestricted page; the unrestricted output equals a run without restrictions
    simoptions_A0=simoptions; simoptions_A0.agegroupings=agegroupings;
    ACS0=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_A0);
    [p0,v0]=CoreStructLeaves(ACS0);
    [pA,vA]=CoreStructLeaves(ACS.always);
    [pU,vU]=CoreStructLeaves(ACS);
    maxd=0; mis=0; maxd2=0; mis2=0;
    for k=1:numel(p0)
        [present,idx]=ismember(p0{k},pA);
        if present && isequal(size(v0{k}),size(vA{idx}))
            a=v0{k}; b=vA{idx}; mis=mis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        else
            mis=mis+1;
        end
        [present,idx]=ismember(p0{k},pU);
        if present && isequal(size(v0{k}),size(vU{idx}))
            a=v0{k}; b=vU{idx}; mis2=mis2+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd2=max([maxd2;d]);
        else
            mis2=mis2+1;
        end
    end
    mis=mis+sum(~ismember(pA,[p0,{'.RestrictedSampleMass'}]));
    fprintf('CondlRest AgeConditional CrossSection %s: page always vs a run without restrictions, every leaf, this should be zero: %.3e, mismatches %i \n',configname,maxd,mis)
    fprintf('CondlRest AgeConditional CrossSection %s: unrestricted output with restrictions vs a run without restrictions, every leaf, this should be zero: %.3e, mismatches %i \n',configname,maxd2,mis2)

    % (5) a restriction with zero mass
    simoptions_N=simoptions_A0; simoptions_N.conditionalrestrictions=NeverRest;
    lastwarn('');
    ACSN=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_N);
    warned=~isempty(lastwarn);
    notnan=sum(~isnan(gather(ACSN.never.CovarianceMatrix(:))))+sum(~isnan(gather(ACSN.never.CorrelationMatrix(:))));
    for ff=1:nFns
        notnan=notnan+sum(~isnan(gather(ACSN.never.(FnNames{ff}).Mean)))+sum(~isnan(gather(ACSN.never.(FnNames{ff}).StdDeviation)));
    end
    fprintf('CondlRest AgeConditional CrossSection %s: zero-mass restriction: no warning + RestrictedSampleMass + Mean/StdDev/matrix entries not NaN, this should be zero: %.3e \n',configname,(1-warned)+sum(abs(gather(ACSN.never.RestrictedSampleMass)))+notnan)

    % (6) whichcombos with a page dimension
    wca=zeros(nFns,nFns,ngroups,npages);
    for pp=1:npages
        for kk=1:ngroups
            for ff1=1:nFns
                for ff2=1:nFns
                    if ff1~=ff2
                        wca(ff1,ff2,kk,pp)=mod(ff1+ff2+kk+pp,2);
                    end
                end
            end
        end
    end
    wca(1,1,:,1)=1; % the unrestricted page: assets' own stats by the diagonal alone
    simoptions_W=simoptions_RA; simoptions_W.whichcombos=wca;
    ACSW=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_W);
    maxd=0; bad=0;
    for pp=1:npages
        if pp==1
            P=ACS; S=ACSW;
        else
            P=ACS.(RestrNames{pp-1}); S=ACSW.(RestrNames{pp-1});
        end
        Mf=gather(P.CovarianceMatrix); Ms=gather(S.CovarianceMatrix); Rf=gather(P.CorrelationMatrix); Rs=gather(S.CorrelationMatrix);
        on=false(nFns,nFns,ngroups);
        for kk=1:ngroups
            on(:,:,kk)=(max(wca(:,:,kk,pp),diag(any(wca(:,:,kk,pp),2)))==1);
        end
        d=[abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]; d=d(~isnan(d)); maxd=max([maxd;d]);
        bad=bad+sum(isnan(Mf(on))~=isnan(Ms(on)))+sum(isnan(Rf(on))~=isnan(Rs(on)))+sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on)));
        for ff=1:nFns
            fnon=reshape(any(wca(ff,:,:,pp),2),[1,ngroups]);
            a=[gather(S.(FnNames{ff}).Mean);gather(S.(FnNames{ff}).StdDeviation)]; b=[gather(P.(FnNames{ff}).Mean);gather(P.(FnNames{ff}).StdDeviation)];
            d=abs(a(:,fnon)-b(:,fnon)); d=d(~isnan(d)); maxd=max([maxd;d]);
            bad=bad+sum(sum(isnan(a(:,fnon))~=isnan(b(:,fnon))))+sum(sum(~isnan(a(:,~fnon))));
        end
    end
    fprintf('CondlRest AgeConditional CrossSection %s, whichcombos with a page dimension: selected entries vs full run, this should be zero: %.3e, unselected not NaN + NaN mismatches %i \n',configname,maxd,bad)
    wca3=wca(:,:,:,1);
    simoptions_W.whichcombos=wca3;
    A=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_W);
    simoptions_W.whichcombos=repmat(wca3,[1,1,1,npages]);
    B=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_W);
    simoptions_W.whichcombos=wc2;
    C2=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_W);
    simoptions_W.whichcombos=repmat(wc2,[1,1,ngroups,npages]);
    D2=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_W);
    [pa,va]=CoreStructLeaves(A); [pb,vb]=CoreStructLeaves(B);
    [pc,vc]=CoreStructLeaves(C2); [pd,vd]=CoreStructLeaves(D2);
    maxd=0; mis=0; maxd2=0; mis2=0;
    for k=1:numel(pa)
        [present,idx]=ismember(pa{k},pb);
        if present && isequal(size(va{k}),size(vb{idx}))
            a=va{k}; b=vb{idx}; mis=mis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        else
            mis=mis+1;
        end
    end
    for k=1:numel(pc)
        [present,idx]=ismember(pc{k},pd);
        if present && isequal(size(vc{k}),size(vd{idx}))
            a=vc{k}; b=vd{idx}; mis2=mis2+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd2=max([maxd2;d]);
        else
            mis2=mis2+1;
        end
    end
    fprintf('CondlRest AgeConditional CrossSection %s, a [nFns,nFns,ngroups] whichcombos vs its repmat over the pages, every leaf, this should be zero: %.3e, mismatches %i \n',configname,maxd,mis+sum(~ismember(pb,pa)))
    fprintf('CondlRest AgeConditional CrossSection %s, a [nFns,nFns] whichcombos vs its repmat over age groups and pages, every leaf, this should be zero: %.3e, mismatches %i \n',configname,maxd2,mis2+sum(~ismember(pd,pc)))
    errored=0;
    simoptions_bad=simoptions_RA; simoptions_bad.whichcombos=ones(nFns,nFns,ngroups,npages-1);
    try
        EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_bad);
    catch
        errored=1;
    end
    fprintf('CondlRest AgeConditional CrossSection %s, a whichcombos with the wrong number of pages must error, this should be zero: %.3e \n',configname,1-errored)

    %% ===================== AllStats and LifeCycleProfiles: every statistic vs brute force =====================
    % The brute force is written from the definitions of the statistics, not from the code of StatsFromWeightedGrid: drop the points
    % with zero mass, aggregate the rest onto their unique values u with masses m (normalised to one), cumulative mass F.
    %   quantile Q(p) = smallest u with F>=p-1e-12 (the 1e-12 resolves an exact tie to the lower value); Minimum/Maximum are Q at
    %   tolerance/1-tolerance without the 1e-12; Median and the percentiles are Q(0.5), Q(0.9), ...; QuantileCutoffs are
    %   [Minimum; Q(1/nq),...,Q(1-1/nq); Maximum]
    %   Lorenz curve (not normalised) L(p) = integral of Q from 0 to p (the mass below the value that straddles p, plus the part of
    %   that value up to p); LorenzCurve is L at 1/npoints,...,1 over L(1); QuantileMeans are the differences of L at the quantile
    %   cutoffs times nq; Top X share is 1-L(1-X)/L(1), Bottom50share is L(0.5)/L(1)
    %   Gini = sum_u m_u u (mass below u - mass above u)/(L(1)*total mass), which is half the mean absolute difference over the mean
    %   one value only: Gini 0, Lorenz curve the 45 degree line; a negative value (and more than one value): Gini, Lorenz curve and
    %   the shares are NaN
    % Checked against AllStats (all ages), LifeCycleProfiles per age (the default), in age groups of five, and in one age group of
    % every age, on every page.
    LCP1=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_R);
    simoptions_R1=simoptions_R; simoptions_R1.agegroupings=1;
    LCPone=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_R1);
    npts=100; nq=20; tol=10^(-12); % the defaults of both commands
    statnames={'Mean','RatioMeanToMedian','Median','MoreInequality.Percentile50th','MoreInequality.Percentile90th','MoreInequality.Percentile95th','MoreInequality.Percentile99th', ...
        'Variance','StdDeviation','Gini','LorenzCurve','Minimum','Maximum','QuantileCutoffs','QuantileMeans', ...
        'MoreInequality.Top1share','MoreInequality.Top5share','MoreInequality.Top10share','MoreInequality.Bottom50share'};
    statgroup=[1,1,2,2,2,2,2,3,3,4,4,5,5,6,6,7,7,7,7];
    groupnames={'Mean/RatioMeanToMedian','Median/percentiles','Variance/StdDeviation','Gini/LorenzCurve','Minimum/Maximum','QuantileCutoffs/QuantileMeans','Top/Bottom shares'};
    for cmd=1:4
        if cmd==1
            cmdname='AllStats'; Out=AS; groups={1:N_j};
        elseif cmd==2
            cmdname='LifeCycleProfiles per age'; Out=LCP1; groups=num2cell(1:N_j);
        elseif cmd==3
            cmdname='LifeCycleProfiles age groups of five'; Out=LCP;
            groups=cell(1,ngroups);
            for kk=1:ngroups
                if kk<ngroups
                    groups{kk}=agegroupings(kk):agegroupings(kk+1)-1;
                else
                    groups{kk}=agegroupings(kk):N_j;
                end
            end
        else
            cmdname='LifeCycleProfiles one age group'; Out=LCPone; groups={1:N_j};
        end
        maxd=zeros(1,7); mism=zeros(1,7); worst=repmat({''},1,7);
        for pp=1:npages
            if pp==1
                OutP=Out; W=SD;
            else
                OutP=Out.(RestrNames{pp-1}); W=SD.*R(:,:,pp-1);
            end
            for ff=1:nFns
                for kk=1:numel(groups)
                    w=reshape(W(:,groups{kk}),[],1);
                    v=reshape(V(:,groups{kk},ff),[],1);
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
                        % quantiles: median, min, max, p99, p95, p90, then the quantile cutoffs
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
                        % Lorenz curve (not normalised) at 1/npts,...,(npts-1)/npts, then 0.99, 0.95, 0.90, 0.5, then the quantile cutoffs
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
                    % compare with the command's output
                    for ss=1:numel(statnames)
                        parts=strsplit(statnames{ss},'.');
                        X=OutP.(FnNames{ff});
                        present=true;
                        for ii=1:numel(parts)
                            if isstruct(X) && isfield(X,parts{ii})
                                X=X.(parts{ii});
                            else
                                present=false;
                            end
                        end
                        gg=statgroup(ss);
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
                            worst{gg}=sprintf('%s, page %s, %s, group %i',statnames{ss},PageNames{pp},FnNames{ff},kk);
                        end
                    end
                end
            end
        end
        for gg=1:7
            fprintf('CondlRest %s %s, %s vs brute force on every page, this should be zero: %.3e, NaN/Inf mismatches + missing fields %i (worst: %s) \n',cmdname,configname,groupnames{gg},maxd(gg),mism(gg),worst{gg})
        end
    end

    % The always page equals the unrestricted output, every leaf (AllStats; LifeCycleProfiles per age)
    for cmd=1:2
        if cmd==1
            cmdname='AllStats'; Out=AS;
        else
            cmdname='LifeCycleProfiles per age'; Out=LCP1;
        end
        maxd=0; mis=0;
        for ff=1:nFns
            [p0,v0]=CoreStructLeaves(Out.(FnNames{ff}));
            [pA,vA]=CoreStructLeaves(Out.always.(FnNames{ff}));
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
        fprintf('CondlRest %s %s: page always vs the unrestricted output, every leaf, this should be zero: %.3e, mismatches %i \n',cmdname,configname,maxd,mis)
    end

    % The atJ1 page of AllStats equals LifeCycleProfiles (unrestricted, per age) at age J1, every leaf
    maxd=0; mis=0;
    for ff=1:nFns
        [pA,vA]=CoreStructLeaves(AS.atJ1.(FnNames{ff}));
        [pL,vL]=CoreStructLeaves(LCP1.(FnNames{ff}));
        for k=1:numel(pA)
            [present,idx]=ismember(pA{k},pL);
            if present && size(vL{idx},2)==N_j && size(vL{idx},1)==numel(vA{k})
                a=vA{k}(:); b=vL{idx}(:,Params.J1); mis=mis+sum(isnan(a)~=isnan(b)); d=abs(a-b); d=d(~isnan(d)); maxd=max([maxd;d]);
            else
                mis=mis+1;
            end
        end
    end
    fprintf('CondlRest AllStats %s: page atJ1 vs LifeCycleProfiles at age J1, every leaf, this should be zero: %.3e, mismatches %i \n',configname,maxd,mis)

    % LifeCycleProfiles with one age group of every age equals AllStats, on every page, every leaf (the restricted distribution of
    % the age group weights each age by its restricted mass)
    maxd=0; mis=0;
    for pp=1:npages
        if pp==1
            OutA=AS; OutL=LCPone;
        else
            OutA=AS.(RestrNames{pp-1}); OutL=LCPone.(RestrNames{pp-1});
        end
        for ff=1:nFns
            [pA,vA]=CoreStructLeaves(OutA.(FnNames{ff}));
            [pL,vL]=CoreStructLeaves(OutL.(FnNames{ff}));
            for k=1:numel(pA)
                [present,idx]=ismember(pA{k},pL);
                if present && numel(vL{idx})==numel(vA{k})
                    a=vA{k}(:); b=vL{idx}(:); mis=mis+sum(isnan(a)~=isnan(b)); d=abs(a-b); d=d(~isnan(d)); maxd=max([maxd;d]);
                else
                    mis=mis+1;
                end
            end
        end
    end
    fprintf('CondlRest LifeCycleProfiles one age group %s vs AllStats, every page, every leaf, this should be zero: %.3e, mismatches %i \n',configname,maxd,mis)

    output.(configname).CrossSection=CS;
    output.(configname).AgeConditionalCrossSection=ACS;
    output.(configname).AllStats=AS;
    output.(configname).LifeCycleProfiles=LCP1;
end

end
