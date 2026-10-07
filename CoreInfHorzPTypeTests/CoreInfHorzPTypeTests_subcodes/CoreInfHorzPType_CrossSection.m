function output=CoreInfHorzPType_CrossSection(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test: EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType (new), with conditional restrictions (one of them different for
% each type), on a model with d and z (two types) and one without d (three types, the middle one of zero weight):
%   per type: CrossSectionCorr.(fn).(typename) and .CovarianceMatrix_ptype.(typename) equal the single-type command on that
%     type, every leaf, and the same under each restriction (.(restriction).(fn).(typename))
%   grouped: Mean, StdDeviation, covariances and correlations (named and matrices) equal the pooling of the single-type outputs,
%     with weights ptweights (unrestricted) and ptweights times the type's restricted sample mass (restricted); the grouped
%     covariance is the within-type covariance plus the between-type part about the pooled means
%   RestrictedSampleMass.(typename), .ByPType and .Total; groupptypesforstats=0 drops the grouped fields
% (The exact reference for the grouped outputs, via z-identity, is part 14.) (InfHorz counterpart of CoreFHorzPTypeTests parts 18-19.)

n_d_withd=n_d; d_grid_withd=d_grid;
for cc=1:2
    FnsToEvaluate=struct(); CondlRest=struct(); vfoptions=struct();
    if cc==1
        mstr='with d, two types';
        Names_i={'low','high'}; ptw=[0.6;0.4]; Params.w_pt=[1;1.4];
        n_d=n_d_withd; d_grid=d_grid_withd;
        ReturnFn=@(d,aprime,a,z,r,w_pt,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w_pt,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a,z) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,w_pt) w_pt*d*z;
        FnsToEvaluate.hours=@(d,aprime,a,z) d;
        CondlRest.positiveassets=@(d,aprime,a,z) (a>0);
        CondlRest.highincome=@(d,aprime,a,z,r,w_pt) (r*a+w_pt*d*z>0.5);
        CondlRest.typespecific.low=@(d,aprime,a,z) (z>1);
        CondlRest.typespecific.high=@(d,aprime,a,z) (a>0.5);
    else
        mstr='without d, three types (a zero-weight type)';
        Names_i={'low','ghost','high'}; ptw=[0.5;0;0.5]; Params.w_pt=[1;1.2;1.4];
        n_d=0; d_grid=[];
        ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
        FnsToEvaluate.consumption=@(aprime,a,z,r,w_pt) (1+r)*a+w_pt*z-aprime;
        CondlRest.positiveassets=@(aprime,a,z) (a>0);
        CondlRest.highz=@(aprime,a,z) (z>1);
        CondlRest.typespecific.low=@(aprime,a,z) (a>0.2);
        CondlRest.typespecific.ghost=@(aprime,a,z) (a>0.2);
        CondlRest.typespecific.high=@(aprime,a,z) (z<1);
    end
    N_i=numel(Names_i);
    Params.(PTypeDistParamNames{1})=ptw;
    FnNames=fieldnames(FnsToEvaluate); nF=numel(FnNames);
    RestrNames=fieldnames(CondlRest); nR=numel(RestrNames);
    [~,Policy_PT]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    Dist_PT=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,Names_i,pi_z,Params,simoptionsPT);
    simoptions_R=struct(); simoptions_R.conditionalrestrictions=CondlRest;
    CS=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_R);
    simoptions_ng=simoptions_R; simoptions_ng.groupptypesforstats=0;
    CSng=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_ng);

    % single-type on each type, with that type's restrictions
    CS_solo=cell(N_i,1);
    for ii=1:N_i
        Params_ii=Params; Params_ii.w_pt=Params.w_pt(ii);
        simoptions_ii=struct();
        for rr=1:nR
            if isstruct(CondlRest.(RestrNames{rr}))
                simoptions_ii.conditionalrestrictions.(RestrNames{rr})=CondlRest.(RestrNames{rr}).(Names_i{ii});
            else
                simoptions_ii.conditionalrestrictions.(RestrNames{rr})=CondlRest.(RestrNames{rr});
            end
        end
        CS_solo{ii}=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(Dist_PT.(Names_i{ii}),Policy_PT.(Names_i{ii}),FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_ii);
    end
    % per type, every leaf, every page
    maxd=0; mis=0;
    for ii=1:N_i
        nm=Names_i{ii};
        for pp=1:1+nR
            if pp==1
                Sp=CS_solo{ii}; Pp=CS;
            else
                Sp=CS_solo{ii}.(RestrNames{pp-1}); Pp=CS.(RestrNames{pp-1});
            end
            for ff=0:nF
                if ff==0
                    a0={gather(Sp.CovarianceMatrix),gather(Sp.CorrelationMatrix)}; b0={gather(Pp.CovarianceMatrix_ptype.(nm)),gather(Pp.CorrelationMatrix_ptype.(nm))};
                    for k=1:2
                        if isequal(size(a0{k}),size(b0{k}))
                            mis=mis+sum(isnan(a0{k}(:))~=isnan(b0{k}(:))); d=abs(a0{k}(:)-b0{k}(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                        else
                            mis=mis+1;
                        end
                    end
                    continue
                end
                [p0,v0]=CoreStructLeaves(Sp.(FnNames{ff})); [p1,v1]=CoreStructLeaves(Pp.(FnNames{ff}).(nm));
                for k=1:numel(p0)
                    [present,idx]=ismember(p0{k},p1);
                    if present && isequal(size(v0{k}),size(v1{idx}))
                        x=v0{k}; y=v1{idx}; mis=mis+sum(isnan(x(:))~=isnan(y(:))); d=abs(x(:)-y(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                    else
                        mis=mis+1;
                    end
                end
            end
        end
    end
    fprintf('CrossSection PType %s, per type vs the single-type command on that type, every page, every leaf, this should be zero: %.3e, mismatches %i \n',mstr,maxd,mis)
    % grouped vs the pooling of the single-type outputs; restricted sample masses
    dm=0; dc=0; nanmis=0; dR=0;
    for pp=1:1+nR
        if pp==1
            Pp=CS;
        else
            Pp=CS.(RestrNames{pp-1});
        end
        wt=ptw;
        if pp>1
            rsm=zeros(N_i,1);
            for ii=1:N_i
                rsm(ii)=gather(CS_solo{ii}.(RestrNames{pp-1}).RestrictedSampleMass);
                dR=max(dR,abs(gather(Pp.RestrictedSampleMass.(Names_i{ii}))-rsm(ii)));
            end
            dR=max([dR,max(abs(reshape(gather(Pp.RestrictedSampleMass.ByPType),[],1)-rsm)),abs(gather(Pp.RestrictedSampleMass.Total)-sum(ptw.*rsm))]);
            wt=ptw.*rsm;
        end
        in=(wt>0); pw=wt(in)/sum(wt(in)); idx=find(in);
        M=zeros(nF,numel(idx)); S=zeros(nF,numel(idx)); C=zeros(nF,nF,numel(idx));
        for kk=1:numel(idx)
            if pp==1
                Sp=CS_solo{idx(kk)};
            else
                Sp=CS_solo{idx(kk)}.(RestrNames{pp-1});
            end
            for ff=1:nF
                M(ff,kk)=gather(Sp.(FnNames{ff}).Mean); S(ff,kk)=gather(Sp.(FnNames{ff}).StdDeviation);
            end
            C(:,:,kk)=gather(Sp.CovarianceMatrix);
        end
        mG=M*pw;
        CG=zeros(nF,nF);
        for kk=1:numel(idx)
            CG=CG+pw(kk)*(C(:,:,kk)+(M(:,kk)-mG)*(M(:,kk)-mG)');
        end
        sG=sqrt(diag(CG));
        RG=CG./(sG*sG');
        for ff=1:nF
            a=[gather(Pp.(FnNames{ff}).Mean);gather(Pp.(FnNames{ff}).StdDeviation)]; b=[mG(ff);sG(ff)];
            nanmis=nanmis+sum(isnan(a)~=isnan(b)); d=abs(a-b); dm=max([dm;d(~isnan(d))]);
            for ff2=1:nF
                if ff2~=ff
                    a=gather(Pp.(FnNames{ff}).CovarianceWith.(FnNames{ff2})); dm=max(dm,abs(a-CG(ff,ff2))); nanmis=nanmis+isnan(a);
                    if sG(ff)*sG(ff2)>1e-10
                        dc=max(dc,abs(gather(Pp.(FnNames{ff}).CorrelationWith.(FnNames{ff2}))-RG(ff,ff2)));
                    end
                end
            end
        end
        a=gather(Pp.CovarianceMatrix); d=abs(a(:)-CG(:)); dm=max([dm;d(~isnan(d))]); nanmis=nanmis+sum(isnan(a(:)));
    end
    fprintf('CrossSection PType %s, grouped Mean/StdDeviation/covariances vs the pooled single-type outputs, every page, this should be zero: %.3e, NaN mismatches %i \n',mstr,dm,nanmis)
    fprintf('CrossSection PType %s, grouped correlations vs the pooled single-type outputs, every page, this should be zero: %.3e \n',mstr,dc)
    fprintf('CrossSection PType %s, RestrictedSampleMass per type, ByPType, Total, this should be zero: %.3e \n',mstr,dR)
    bad=0;
    for ff=1:nF
        bad=bad+isfield(CSng.(FnNames{ff}),'Mean')+isfield(CSng,'CovarianceMatrix');
    end
    fprintf('CrossSection PType %s, groupptypesforstats=0 has no grouped fields, this should be zero: %i \n',mstr,bad)
end
output=struct();
end
