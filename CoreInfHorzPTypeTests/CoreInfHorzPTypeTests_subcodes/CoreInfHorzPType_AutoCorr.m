function output=CoreInfHorzPType_AutoCorr(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test: EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType (new), with simoptions.timehorizons=[2,3], transition probabilities for
% assets, and conditional restrictions (one of them different for each type), on a model with d and z (two types) and one without
% d (three types, the middle one of zero weight):
%   per type: CorrTransProbs.(fn).(typename) equals the single-type command on that type, every leaf (TransitionProbs included),
%     and the same under each restriction
%   grouped, unrestricted: Mean and StdDeviation are the pooling with weights ptweights; the auto-covariance at horizon k is
%     sum_i ptweights(i)*(AutoCovariance_i + (Mean_i-Mean)^2) (each type is stationary, so its mean is the same now and k periods on),
%     the auto-correlation that over the pooled variance
%   grouped, restricted: Mean and StdDeviation pool with weights ptweights times the restricted sample mass; the pair outputs pool
%     the types' pair populations, weights ptweights times PairMass (so the grouped PairMass is sum_i ptweights(i)*PairMass_i, a
%     population share), the pair means and std devs and the covariance about the pooled pair means
%   RestrictedSampleMass.(typename), .ByPType and .Total; groupptypesforstats=0 drops the grouped fields
% (The exact reference for the grouped outputs, via z-identity, is part 14.)

timehorizons=[2,3]; horizons=[1,2,3];
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
        CondlRest.positiveassets=@(d,aprime,a,z) (a>0);
        CondlRest.typespecific.low=@(d,aprime,a,z) (z>1);
        CondlRest.typespecific.high=@(d,aprime,a,z) (a>0.5);
    else
        mstr='without d, three types (a zero-weight type)';
        Names_i={'low','ghost','high'}; ptw=[0.5;0;0.5]; Params.w_pt=[1;1.2;1.4];
        n_d=0; d_grid=[];
        ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.consumption=@(aprime,a,z,r,w_pt) (1+r)*a+w_pt*z-aprime;
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
    simoptions_A=struct(); simoptions_A.timehorizons=timehorizons; simoptions_A.transprobs={'assets'}; simoptions_A.conditionalrestrictions=CondlRest;
    ACP=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_A);
    simoptions_ng=simoptions_A; simoptions_ng.groupptypesforstats=0;
    ACPng=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ng);

    ACP_solo=cell(N_i,1);
    for ii=1:N_i
        Params_ii=Params; Params_ii.w_pt=Params.w_pt(ii);
        simoptions_ii=simoptions_A; simoptions_ii.conditionalrestrictions=struct();
        for rr=1:nR
            if isstruct(CondlRest.(RestrNames{rr}))
                simoptions_ii.conditionalrestrictions.(RestrNames{rr})=CondlRest.(RestrNames{rr}).(Names_i{ii});
            else
                simoptions_ii.conditionalrestrictions.(RestrNames{rr})=CondlRest.(RestrNames{rr});
            end
        end
        ACP_solo{ii}=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(Dist_PT.(Names_i{ii}),Policy_PT.(Names_i{ii}),FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptions_ii);
    end
    % per type, every leaf, every page
    maxd=0; mis=0;
    for ii=1:N_i
        for pp=1:1+nR
            for ff=1:nF
                if pp==1
                    Sp=ACP_solo{ii}.(FnNames{ff}); Pp=ACP.(FnNames{ff}).(Names_i{ii});
                else
                    Sp=ACP_solo{ii}.(RestrNames{pp-1}).(FnNames{ff}); Pp=ACP.(RestrNames{pp-1}).(FnNames{ff}).(Names_i{ii});
                end
                [p0,v0]=CoreStructLeaves(Sp); [p1,v1]=CoreStructLeaves(Pp);
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
    fprintf('AutoCorr PType %s, per type vs the single-type command on that type, every page, every leaf, this should be zero: %.3e, mismatches %i \n',mstr,maxd,mis)
    % grouped, unrestricted
    dm=0; dc=0; nanmis=0;
    for ff=1:nF
        fn=FnNames{ff};
        M=zeros(N_i,1); S=zeros(N_i,1);
        for ii=1:N_i
            M(ii)=gather(ACP_solo{ii}.(fn).Mean); S(ii)=gather(ACP_solo{ii}.(fn).StdDeviation);
        end
        in=(ptw>0); pw=ptw(in)/sum(ptw(in));
        mG=sum(pw.*M(in)); vG=sum(pw.*(S(in).^2+(M(in)-mG).^2));
        a=[gather(ACP.(fn).Mean);gather(ACP.(fn).StdDeviation)]; b=[mG;sqrt(vG)];
        nanmis=nanmis+sum(isnan(a)~=isnan(b)); d=abs(a-b); dm=max([dm;d(~isnan(d))]);
        for hh=1:numel(horizons)
            Cv=zeros(N_i,1);
            for ii=1:N_i
                if horizons(hh)==1
                    Cv(ii)=gather(ACP_solo{ii}.(fn).AutoCovariance);
                else
                    Cv(ii)=gather(ACP_solo{ii}.(fn).(['tperiods',num2str(horizons(hh))]).AutoCovariance);
                end
            end
            cG=sum(pw.*(Cv(in)+(M(in)-mG).^2));
            if horizons(hh)==1
                G=ACP.(fn);
            else
                G=ACP.(fn).(['tperiods',num2str(horizons(hh))]);
            end
            dm=max(dm,abs(gather(G.AutoCovariance)-cG)); nanmis=nanmis+isnan(gather(G.AutoCovariance));
            if vG>1e-10
                dc=max(dc,abs(gather(G.AutoCorrelation)-cG/vG));
            end
        end
    end
    fprintf('AutoCorr PType %s, grouped unrestricted Mean/StdDeviation/AutoCovariance (horizons 1,2,3) vs the pooled single-type outputs, this should be zero: %.3e, NaN mismatches %i \n',mstr,dm,nanmis)
    fprintf('AutoCorr PType %s, grouped unrestricted AutoCorrelation vs the pooled single-type outputs, this should be zero: %.3e \n',mstr,dc)
    % grouped, restricted
    dm=0; dc=0; nanmis=0; dR=0;
    for rr=1:nR
        rname=RestrNames{rr};
        rsm=zeros(N_i,1);
        for ii=1:N_i
            rsm(ii)=gather(ACP_solo{ii}.(rname).RestrictedSampleMass);
            dR=max(dR,abs(gather(ACP.(rname).RestrictedSampleMass.(Names_i{ii}))-rsm(ii)));
        end
        dR=max([dR,max(abs(reshape(gather(ACP.(rname).RestrictedSampleMass.ByPType),[],1)-rsm)),abs(gather(ACP.(rname).RestrictedSampleMass.Total)-sum(ptw.*rsm))]);
        for ff=1:nF
            fn=FnNames{ff};
            M=zeros(N_i,1); S=zeros(N_i,1);
            for ii=1:N_i
                M(ii)=gather(ACP_solo{ii}.(rname).(fn).Mean); S(ii)=gather(ACP_solo{ii}.(rname).(fn).StdDeviation);
            end
            w=ptw.*rsm; in=(w>0); pw=w(in)/sum(w(in));
            mG=sum(pw.*M(in)); sG=sqrt(sum(pw.*(S(in).^2+(M(in)-mG).^2)));
            a=[gather(ACP.(rname).(fn).Mean);gather(ACP.(rname).(fn).StdDeviation)]; b=[mG;sG];
            nanmis=nanmis+sum(isnan(a)~=isnan(b)); d=abs(a-b); dm=max([dm;d(~isnan(d))]);
            for hh=1:numel(horizons)
                PM=zeros(N_i,1); MX=PM; MY=PM; SX=PM; SY=PM; C=PM;
                for ii=1:N_i
                    if horizons(hh)==1
                        Q=ACP_solo{ii}.(rname).(fn);
                    else
                        Q=ACP_solo{ii}.(rname).(fn).(['tperiods',num2str(horizons(hh))]);
                    end
                    PM(ii)=ptw(ii)*gather(Q.PairMass); MX(ii)=gather(Q.PairMean_t); MY(ii)=gather(Q.PairMean_tplusk);
                    SX(ii)=gather(Q.PairStdDeviation_t); SY(ii)=gather(Q.PairStdDeviation_tplusk); C(ii)=gather(Q.AutoCovariance);
                end
                in=(PM>0); pw=PM(in)/sum(PM(in));
                mx=sum(pw.*MX(in)); my=sum(pw.*MY(in));
                vx=sum(pw.*(SX(in).^2+(MX(in)-mx).^2)); vy=sum(pw.*(SY(in).^2+(MY(in)-my).^2));
                cG=sum(pw.*(C(in)+(MX(in)-mx).*(MY(in)-my)));
                if horizons(hh)==1
                    G=ACP.(rname).(fn);
                else
                    G=ACP.(rname).(fn).(['tperiods',num2str(horizons(hh))]);
                end
                a=gather([G.PairMass;G.PairMean_t;G.PairMean_tplusk;G.PairStdDeviation_t;G.PairStdDeviation_tplusk;G.AutoCovariance]);
                b=[sum(PM);mx;my;sqrt(vx);sqrt(vy);cG];
                nanmis=nanmis+sum(isnan(a)~=isnan(b)); d=abs(a-b); dm=max([dm;d(~isnan(d))]);
                if sqrt(vx*vy)>1e-10
                    dc=max(dc,abs(gather(G.AutoCorrelation)-cG/sqrt(vx*vy)));
                end
            end
        end
    end
    fprintf('AutoCorr PType %s, grouped restricted Mean/StdDeviation/pair outputs (horizons 1,2,3) vs the pooled single-type outputs, this should be zero: %.3e, NaN mismatches %i \n',mstr,dm,nanmis)
    fprintf('AutoCorr PType %s, grouped restricted AutoCorrelation vs the pooled single-type outputs, this should be zero: %.3e \n',mstr,dc)
    fprintf('AutoCorr PType %s, RestrictedSampleMass per type, ByPType, Total, this should be zero: %.3e \n',mstr,dR)
    bad=0;
    for ff=1:nF
        bad=bad+isfield(ACPng.(FnNames{ff}),'Mean')+isfield(ACPng.(FnNames{ff}),'AutoCovariance');
    end
    fprintf('AutoCorr PType %s, groupptypesforstats=0 has no grouped fields, this should be zero: %i \n',mstr,bad)
end
output=struct();
end
