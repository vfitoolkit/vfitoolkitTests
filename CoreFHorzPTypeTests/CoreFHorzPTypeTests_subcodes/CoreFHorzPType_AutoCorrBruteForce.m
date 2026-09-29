function output=CoreFHorzPType_AutoCorrBruteForce(n_a,n_z,n_e,N_j,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: EvalFnOnAgentDist_AutoCorrTransProbs_FHorz with simoptions.timehorizons and simoptions.conditionalrestrictions
% (added 2026-09-29), and its PType wrapper, against a brute-force calculation written here independently of the
% toolkit's transition-matrix code.
%
% Brute force: for each age jj the one-period transition kernel T_jj (N_states x N_states, dense) is built by loops
% from the Policy indexes and pi_z/pi_e (with the two-point grid-interpolation weights when gridinterplayer=1), the
% k-period kernel is the product T_j0*...*T_{j0+k-1}, and the joint distribution of (state at j0, state at j0+k) is
% dist_j0(s)*Pk(s,s'). Every moment is then a plain sum over that joint distribution: means, std devs and the
% covariance of x_j0 and x_{j0+k}; under a restriction the joint distribution is masked to the pairs that satisfy
% it at both ages and renormalised (the pair population), which is the definition the command documents.
%
% Model: no d, markov z, iid e (state (a,z,e), N_states=n_a*n_z*n_e). Checks, all exact identities:
%  A1. horizon-1 outputs (Mean, StdDeviation, AutoCovariance, AutoCorrelation, TransitionProbs) unchanged vs the
%      pre-change command (a renamed copy of the command at toolkit commit d6ebef06 kept in this bank), on a model
%      without e and on the model with e, lowmemory=0 and =1.
%  A2. horizons 1,2,3,5, unrestricted and restricted (rich: a>0.5; highz: z>1), vs the brute force: AutoCovariance,
%      AutoCorrelation and, restricted, PairMass, PairMean_j, PairMean_jplusk, PairStdDeviation_j,
%      PairStdDeviation_jplusk, Mean, StdDeviation.
%  A3. lowmemory=1 equals lowmemory=0 (every field).
%  A4. gridinterplayer=1: the same brute-force checks with the two-point kernel.
%  A5. a cohort distribution started at age 10 (simoptions.jequaloneDistAge=10, from the age-10 slice of the full
%      run): NaN at every age/pair before 10, equal to the full run from 10 on, no error.
%  B.  PType (two types with different kappa_j and DIFFERENT age weights, ptweights 0.6/0.4): per-type outputs equal
%      the solo command; grouped Mean/StdDeviation/AutoCovariance/AutoCorrelation (unrestricted and restricted, all
%      horizons) equal a by-hand pooling of the two types' brute-force joint distributions with weights
%      ptweights*(type's mass at j0); RestrictedSampleMass.ByAge by hand; a type of zero weight leaves the grouped
%      output equal to the other type alone; groupptypesforstats=0 drops the grouped fields.
% Zero checks are printed as %.3e; NaN is compared as a pattern (count of entries NaN in one but not the other).

fprintf('\n========== CoreFHorzPType_AutoCorrBruteForce ==========\n')

n_d=0;
d_grid=[];
N_s=n_a*n_z*n_e;

vfoptions=struct();
vfoptions.n_e=n_e; vfoptions.e_grid=e_grid; vfoptions.pi_e=pi_e;
simoptions=struct();
simoptions.n_e=n_e; simoptions.e_grid=e_grid; simoptions.pi_e=pi_e;

ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate.assets=@(aprime,a,z,e) a;
FnsToEvaluate.earnings=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
FnsToEvaluate.consumption=@(aprime,a,z,e,r,w,kappa_j) (1+r)*a+w*kappa_j*z*e-aprime;
FnsToEvaluate.evar=@(aprime,a,z,e) e;
FnNames=fieldnames(FnsToEvaluate);

Restrictions.rich=@(aprime,a,z,e) (a>0.5);
Restrictions.highz=@(aprime,a,z,e) (z>1);
RestNames=fieldnames(Restrictions);

jequaloneDist=zeros(n_a,n_z,n_e,'gpuArray');
jequaloneDist(1,ceil(n_z/2),ceil(n_e/2))=1;

timehorizons=[2,3,5];
horizons=[1,timehorizons];
hstr={'','_k2','_k3','_k5'};
Kmax=max(horizons);

pi_z_cpu=gather(pi_z);
pi_e_cpu=gather(pi_e(:));

%% Solve, with and without the grid interpolation layer
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

vfoptions_gi=vfoptions; vfoptions_gi.gridinterplayer=1; vfoptions_gi.ngridinterp=5;
simoptions_gi=simoptions; simoptions_gi.gridinterplayer=1; simoptions_gi.ngridinterp=5;
[~,Policy_gi]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_gi);
StationaryDist_gi=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_gi,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_gi);

%% A1. Horizon-1 outputs unchanged vs the pre-change command
% (i) a model without e
vfoptions_noe=struct();
simoptions_noe=struct();
ReturnFn_noe=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate_noe.assets=@(aprime,a,z) a;
FnsToEvaluate_noe.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate_noe.consumption=@(aprime,a,z,r,w,kappa_j) (1+r)*a+w*kappa_j*z-aprime;
jequaloneDist_noe=zeros(n_a,n_z,'gpuArray');
jequaloneDist_noe(1,ceil(n_z/2))=1;
[~,Policy_noe]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_noe,Params,DiscountFactorParamNames,[],vfoptions_noe);
StationaryDist_noe=StationaryDist_FHorz_Case1(jequaloneDist_noe,AgeWeightParamNames,Policy_noe,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_noe);
for lowmem=0:1
    simoptions_a1=simoptions_noe;
    simoptions_a1.lowmemory=lowmem;
    simoptions_a1.transprobs={'assets','earnings'};
    ACP_new=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist_noe,Policy_noe,FnsToEvaluate_noe,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_a1);
    ACP_old=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_prechange(StationaryDist_noe,Policy_noe,FnsToEvaluate_noe,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_a1);
    fnames_noe=fieldnames(FnsToEvaluate_noe);
    for ff=1:length(fnames_noe)
        fn=fnames_noe{ff};
        statnames={'Mean','StdDeviation','AutoCovariance','AutoCorrelation'};
        for ss=1:length(statnames)
            x=gather(ACP_new.(fn).(statnames{ss})(:)); y=gather(ACP_old.(fn).(statnames{ss})(:));
            dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force A1 (no e, lowmemory=%i), %s %s new vs pre-change command, this should be zero: %.3e, NaN-pattern mismatches %i \n',lowmem,fn,statnames{ss},max([dd;0]),sum(isnan(x)~=isnan(y)))
        end
        if isfield(ACP_new.(fn),'TransitionProbs')
            dtp=0;
            for jj=1:N_j-1
                if ~isempty(ACP_old.(fn).TransitionProbs{jj})
                    dtp=max(dtp,max(abs(ACP_new.(fn).TransitionProbs{jj}(:)-ACP_old.(fn).TransitionProbs{jj}(:))));
                end
            end
            fprintf('AutoCorr brute force A1 (no e, lowmemory=%i), %s TransitionProbs new vs pre-change command, this should be zero: %.3e \n',lowmem,fn,dtp)
        end
    end
end
% (ii) the model with e
for lowmem=0:1
    simoptions_a1=simoptions;
    simoptions_a1.lowmemory=lowmem;
    ACP_new=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_a1);
    ACP_old=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_prechange(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_a1);
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        statnames={'Mean','StdDeviation','AutoCovariance','AutoCorrelation'};
        for ss=1:length(statnames)
            x=gather(ACP_new.(fn).(statnames{ss})(:)); y=gather(ACP_old.(fn).(statnames{ss})(:));
            dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force A1 (with e, lowmemory=%i), %s %s new vs pre-change command, this should be zero: %.3e, NaN-pattern mismatches %i \n',lowmem,fn,statnames{ss},max([dd;0]),sum(isnan(x)~=isnan(y)))
        end
    end
end

%% A2-A4. Brute force, without and with the grid interpolation layer
for cfg=1:2
    if cfg==1
        cfgstr='no GI';
        Policy_c=Policy; Dist_c=StationaryDist; simoptions_c=simoptions; ngridinterp=0;
    else
        cfgstr='GI';
        Policy_c=Policy_gi; Dist_c=StationaryDist_gi; simoptions_c=simoptions_gi; ngridinterp=simoptions_gi.ngridinterp;
    end
    simoptions_c.timehorizons=timehorizons;
    simoptions_c.conditionalrestrictions=Restrictions;
    ACP=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_c,Policy_c,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_c);
    simoptions_c1=simoptions_c; simoptions_c1.lowmemory=1;
    ACP1=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_c,Policy_c,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_c1);

    % Values of every function and restriction on the grid, (N_s, N_j)
    ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_c,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_c);
    Vals=cell(length(FnNames),1);
    for ff=1:length(FnNames)
        Vals{ff}=gpuArray(reshape(gather(ValuesOnGrid.(FnNames{ff})),[N_s,N_j]));
    end
    RestOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_c,Restrictions,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_c);
    Rvals=cell(length(RestNames),1);
    for rr=1:length(RestNames)
        Rvals{rr}=gpuArray(double(reshape(gather(RestOnGrid.(RestNames{rr})),[N_s,N_j])~=0));
    end
    Dist_u=gpuArray(reshape(gather(Dist_c),[N_s,N_j])); % unnormalised (includes the age weights)
    agemass=gather(sum(Dist_u,1));

    % One-period kernels T{jj}: state s=(a,z,e) in column-major order, s' likewise; T(s,s')=1{a'=aprime(s,jj)} pi_z(z,z') pi_e(e')
    Policy_cpu=gather(Policy_c);
    T=cell(N_j-1,1);
    for jj=1:N_j-1
        Tjj=zeros(N_s,N_s);
        for zz=1:n_z
            for ee=1:n_e
                rows=(1:n_a)'+n_a*(zz-1)+n_a*n_z*(ee-1);
                aplow=reshape(Policy_cpu(1,:,zz,ee,jj),[n_a,1]);
                if ngridinterp==0
                    pup=zeros(n_a,1);
                else
                    L2index=reshape(Policy_cpu(2,:,zz,ee,jj),[n_a,1]);
                    L2flag=reshape(Policy_cpu(3,:,zz,ee,jj),[n_a,1]);
                    pup=(L2index-1)/(ngridinterp+1);
                    pup(L2flag==1)=0; % all weight on the lower grid point
                    pup(L2flag==3)=1; % all weight on the upper grid point
                end
                for zz2=1:n_z
                    for ee2=1:n_e
                        prob=pi_z_cpu(zz,zz2)*pi_e_cpu(ee2);
                        cols=aplow+n_a*(zz2-1)+n_a*n_z*(ee2-1);
                        Tjj(rows+(cols-1)*N_s)=Tjj(rows+(cols-1)*N_s)+prob*(1-pup);
                        if ngridinterp>0
                            Tjj(rows+cols*N_s)=Tjj(rows+cols*N_s)+prob*pup; % upper point is aplow+1, so one column further
                        end
                    end
                end
            end
        end
        T{jj}=gpuArray(Tjj);
    end
    fprintf('AutoCorr brute force (%s), kernels row sums minus one, this should be zero: %.3e \n',cfgstr,max(cellfun(@(M) gather(max(abs(sum(M,2)-1))),T)))

    % Restricted per-age Mean and StdDeviation
    for rr=1:length(RestNames)
        for ff=1:length(FnNames)
            MeanR_bf=nan(1,N_j); StdR_bf=nan(1,N_j);
            for jj=1:N_j
                mr=Dist_u(:,jj).*Rvals{rr}(:,jj);
                if sum(mr)>0
                    mr=mr/sum(mr);
                    MeanR_bf(jj)=gather(sum(mr.*Vals{ff}(:,jj)));
                    StdR_bf(jj)=gather(sqrt(sum(mr.*(Vals{ff}(:,jj)-MeanR_bf(jj)).^2)));
                end
            end
            x=gather(ACP.(RestNames{rr}).(FnNames{ff}).Mean(:)); y=MeanR_bf(:); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force A2 (%s), %s %s Mean vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,RestNames{rr},FnNames{ff},max([dd;0]),sum(isnan(x)~=isnan(y)))
            x=gather(ACP.(RestNames{rr}).(FnNames{ff}).StdDeviation(:)); y=StdR_bf(:); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force A2 (%s), %s %s StdDeviation vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,RestNames{rr},FnNames{ff},max([dd;0]),sum(isnan(x)~=isnan(y)))
        end
    end

    % Pairs (j0, j0+k) for every horizon
    nR=length(RestNames);
    AutoCov_bf=nan(length(FnNames),length(horizons),N_j-1);
    AutoCorr_bf=nan(length(FnNames),length(horizons),N_j-1);
    AutoCovR_bf=nan(length(FnNames),length(horizons),N_j-1,nR);
    AutoCorrR_bf=nan(length(FnNames),length(horizons),N_j-1,nR);
    PairMass_bf=nan(length(FnNames),length(horizons),N_j-1,nR);
    PairMeanx_bf=nan(length(FnNames),length(horizons),N_j-1,nR);
    PairMeany_bf=nan(length(FnNames),length(horizons),N_j-1,nR);
    PairStdx_bf=nan(length(FnNames),length(horizons),N_j-1,nR);
    PairStdy_bf=nan(length(FnNames),length(horizons),N_j-1,nR);
    for j0=1:N_j-1
        Pk=eye(N_s,'gpuArray');
        for kk=1:min(Kmax,N_j-j0)
            Pk=Pk*T{j0+kk-1};
            hh=find(horizons==kk);
            if isempty(hh)
                continue
            end
            J=(Dist_u(:,j0)/agemass(j0)).*Pk; % joint distribution of (state at j0, state at j0+kk), mass one
            for ff=1:length(FnNames)
                x=Vals{ff}(:,j0); y=Vals{ff}(:,j0+kk);
                mj=sum(J,2); mk=sum(J,1)';
                mux=sum(mj.*x); muy=sum(mk.*y);
                sdx=sqrt(sum(mj.*(x-mux).^2)); sdy=sqrt(sum(mk.*(y-muy).^2));
                cov=sum(sum(J.*((x-mux)*(y-muy)')));
                AutoCov_bf(ff,hh,j0)=gather(cov);
                if sdx*sdy>1e-15
                    AutoCorr_bf(ff,hh,j0)=gather(cov/(sdx*sdy));
                end
                for rr=1:nR
                    Jr=J.*(Rvals{rr}(:,j0)*Rvals{rr}(:,j0+kk)');
                    pm=sum(Jr(:));
                    PairMass_bf(ff,hh,j0,rr)=gather(pm*agemass(j0));
                    if pm>0
                        Jr=Jr/pm;
                        mj=sum(Jr,2); mk=sum(Jr,1)';
                        mux=sum(mj.*x); muy=sum(mk.*y);
                        sdx=sqrt(sum(mj.*(x-mux).^2)); sdy=sqrt(sum(mk.*(y-muy).^2));
                        cov=sum(sum(Jr.*((x-mux)*(y-muy)')));
                        PairMeanx_bf(ff,hh,j0,rr)=gather(mux); PairMeany_bf(ff,hh,j0,rr)=gather(muy);
                        PairStdx_bf(ff,hh,j0,rr)=gather(sdx); PairStdy_bf(ff,hh,j0,rr)=gather(sdy);
                        AutoCovR_bf(ff,hh,j0,rr)=gather(cov);
                        if sdx*sdy>1e-15
                            AutoCorrR_bf(ff,hh,j0,rr)=gather(cov/(sdx*sdy));
                        end
                    end
                end
            end
        end
    end

    % Compare
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        for hh=1:length(horizons)
            kk=horizons(hh);
            x=gather(ACP.(fn).(['AutoCovariance',hstr{hh}])(:)); y=reshape(AutoCov_bf(ff,hh,1:N_j-kk),[],1); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force A2 (%s), %s AutoCovariance%s vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,fn,hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(y)))
            x=gather(ACP.(fn).(['AutoCorrelation',hstr{hh}])(:)); y=reshape(AutoCorr_bf(ff,hh,1:N_j-kk),[],1); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force A2 (%s), %s AutoCorrelation%s vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,fn,hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(y)))
            for rr=1:nR
                rn=RestNames{rr};
                bfnames={'AutoCovariance','AutoCorrelation','PairMass','PairMean_j','PairMean_jplusk','PairStdDeviation_j','PairStdDeviation_jplusk'};
                bfvals={AutoCovR_bf,AutoCorrR_bf,PairMass_bf,PairMeanx_bf,PairMeany_bf,PairStdx_bf,PairStdy_bf};
                for bb=1:length(bfnames)
                    x=gather(ACP.(rn).(fn).([bfnames{bb},hstr{hh}])(:)); y=reshape(bfvals{bb}(ff,hh,1:N_j-kk,rr),[],1); dd=abs(x-y); dd=dd(~isnan(dd));
                    fprintf('AutoCorr brute force A2 (%s), %s %s %s%s vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,rn,fn,bfnames{bb},hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(y)))
                end
            end
        end
    end

    % A3. lowmemory=1 vs lowmemory=0, every field
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        statnames=fieldnames(ACP.(fn));
        for ss=1:length(statnames)
            x=gather(ACP1.(fn).(statnames{ss})(:)); y=gather(ACP.(fn).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force A3 (%s), %s %s lowmemory=1 vs =0, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,fn,statnames{ss},max([dd;0]),sum(isnan(x)~=isnan(y)))
        end
        for rr=1:nR
            rn=RestNames{rr};
            statnames=fieldnames(ACP.(rn).(fn));
            for ss=1:length(statnames)
                x=gather(ACP1.(rn).(fn).(statnames{ss})(:)); y=gather(ACP.(rn).(fn).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
                fprintf('AutoCorr brute force A3 (%s), %s %s %s lowmemory=1 vs =0, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,rn,fn,statnames{ss},max([dd;0]),sum(isnan(x)~=isnan(y)))
            end
        end
    end
    for rr=1:nR
        rn=RestNames{rr};
        fprintf('AutoCorr brute force A3 (%s), %s RestrictedSampleMass lowmemory=1 vs =0, this should be zero: %.3e \n',cfgstr,rn,max(abs(gather(ACP1.(rn).RestrictedSampleMass(:))-gather(ACP.(rn).RestrictedSampleMass(:)))))
        fprintf('AutoCorr brute force A2 (%s), %s RestrictedSampleMass vs brute force, this should be zero: %.3e \n',cfgstr,rn,max(abs(gather(ACP.(rn).RestrictedSampleMass(:))-gather(sum(Dist_u.*Rvals{rr},1))')))
    end

    if cfg==1
        ACP_full=ACP; % kept for A5
    end
end

%% A5. Cohort distribution started at age 10: NaN before, equal to the full run from age 10 on
j0c=10;
jequaloneDist10=reshape(gather(StationaryDist),[N_s,N_j]);
jequaloneDist10=gpuArray(reshape(jequaloneDist10(:,j0c)/sum(jequaloneDist10(:,j0c)),[n_a,n_z,n_e]));
simoptions_c10=simoptions;
simoptions_c10.jequaloneDistAge=j0c;
Dist10=StationaryDist_FHorz_Case1(jequaloneDist10,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_c10);
simoptions_c10.timehorizons=timehorizons;
simoptions_c10.conditionalrestrictions=Restrictions;
ACP10=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist10,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_c10);
fprintf('AutoCorr brute force A5, cohort from age %i, Dist mass at ages 1..%i, this should be zero: %.3e \n',j0c,j0c-1,sum(sum(reshape(gather(Dist10),[N_s,N_j]),1)*[ones(j0c-1,1);zeros(N_j-j0c+1,1)]))
for ff=1:length(FnNames)
    fn=FnNames{ff};
    fprintf('AutoCorr brute force A5, cohort from age %i, %s Mean non-NaN entries before age %i, this should be zero: %i \n',j0c,fn,j0c,sum(~isnan(gather(ACP10.(fn).Mean(1:j0c-1)))))
    fprintf('AutoCorr brute force A5, cohort from age %i, %s Mean from age %i on vs full run, this should be zero: %.3e \n',j0c,fn,j0c,max(abs(gather(ACP10.(fn).Mean(j0c:end))-gather(ACP_full.(fn).Mean(j0c:end)))))
    for hh=1:length(horizons)
        kk=horizons(hh);
        fprintf('AutoCorr brute force A5, cohort from age %i, %s AutoCovariance%s non-NaN entries before age %i, this should be zero: %i \n',j0c,fn,hstr{hh},j0c,sum(~isnan(gather(ACP10.(fn).(['AutoCovariance',hstr{hh}])(1:j0c-1)))))
        x=gather(ACP10.(fn).(['AutoCovariance',hstr{hh}])(j0c:N_j-kk)); y=gather(ACP_full.(fn).(['AutoCovariance',hstr{hh}])(j0c:N_j-kk)); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force A5, cohort from age %i, %s AutoCovariance%s from age %i on vs full run, this should be zero: %.3e, NaN-pattern mismatches %i \n',j0c,fn,hstr{hh},j0c,max([dd,0]),sum(isnan(x)~=isnan(y)))
        x=gather(ACP10.rich.(fn).(['AutoCovariance',hstr{hh}])(j0c:N_j-kk)); y=gather(ACP_full.rich.(fn).(['AutoCovariance',hstr{hh}])(j0c:N_j-kk)); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force A5, cohort from age %i, rich %s AutoCovariance%s from age %i on vs full run, this should be zero: %.3e, NaN-pattern mismatches %i \n',j0c,fn,hstr{hh},j0c,max([dd,0]),sum(isnan(x)~=isnan(y)))
        fprintf('AutoCorr brute force A5, cohort from age %i, rich %s PairMass%s non-NaN entries before age %i, this should be zero: %i \n',j0c,fn,hstr{hh},j0c,sum(~isnan(gather(ACP10.rich.(fn).(['PairMass',hstr{hh}])(1:j0c-1)))))
    end
end

%% B. PType: two types with different kappa_j and different age weights
Names_i={'low','high'};
N_i=2;
Params_PT=Params;
Params_PT.kappa_j_pt=[Params.kappa_j; 1.3*Params.kappa_j];
Params_PT.mewj=struct();
Params_PT.mewj.low=Params.mewj;
Params_PT.mewj.high=0.9.^(0:N_j-1); Params_PT.mewj.high=Params_PT.mewj.high/sum(Params_PT.mewj.high);
Params_PT.(PTypeDistParamNames{1})=[0.6;0.4];
ptw=[0.6;0.4];
ReturnFn_PT=@(aprime,a,z,e,r,w,kappa_j_pt,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j_pt,sigma,agej,Jr,pension);
FnsToEvaluate_PT.assets=@(aprime,a,z,e) a;
FnsToEvaluate_PT.earnings=@(aprime,a,z,e,w,kappa_j_pt) w*kappa_j_pt*z*e;
FnsToEvaluate_PT.consumption=@(aprime,a,z,e,r,w,kappa_j_pt) (1+r)*a+w*kappa_j_pt*z*e-aprime;
FnsToEvaluate_PT.evar=@(aprime,a,z,e) e;
timehorizons_PT=[2,3];
horizons_PT=[1,2,3];
hstr_PT={'','_k2','_k3'};

simoptions_PT=simoptions;
simoptions_PT.timehorizons=timehorizons_PT;
simoptions_PT.conditionalrestrictions=Restrictions;
[~,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params_PT,DiscountFactorParamNames,vfoptions);
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_PT,simoptions);
ACP_PT=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_PT);
simoptions_PT_ng=simoptions_PT; simoptions_PT_ng.groupptypesforstats=0;
ACP_PT_ng=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_PT_ng);

% Solo solves, and their brute-force joint distributions
ACP_solo=struct(); Dist_solo=struct(); Vals_solo=struct(); Rvals_solo=struct(); T_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Params_ii=Params;
    Params_ii.kappa_j=Params_PT.kappa_j_pt(ii,:);
    Params_ii.mewj=Params_PT.mewj.(nm);
    [~,Policy_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_ii,n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    ACP_solo.(nm)=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_solo.(nm),Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_PT);
    fprintf('AutoCorr brute force B, Dist (type %s) PType vs solo, this should be zero: %.3e \n',nm,max(abs(gather(StationaryDist_PT.(nm)(:))-gather(Dist_solo.(nm)(:)))))
    ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    Vals_solo.(nm)=cell(length(FnNames),1);
    for ff=1:length(FnNames)
        Vals_solo.(nm){ff}=gpuArray(reshape(gather(ValuesOnGrid.(FnNames{ff})),[N_s,N_j]));
    end
    RestOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_ii,Restrictions,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    Rvals_solo.(nm)=cell(length(RestNames),1);
    for rr=1:length(RestNames)
        Rvals_solo.(nm){rr}=gpuArray(double(reshape(gather(RestOnGrid.(RestNames{rr})),[N_s,N_j])~=0));
    end
    Policy_cpu=gather(Policy_ii);
    T_solo.(nm)=cell(N_j-1,1);
    for jj=1:N_j-1
        Tjj=zeros(N_s,N_s);
        for zz=1:n_z
            for ee=1:n_e
                rows=(1:n_a)'+n_a*(zz-1)+n_a*n_z*(ee-1);
                aplow=reshape(Policy_cpu(1,:,zz,ee,jj),[n_a,1]);
                for zz2=1:n_z
                    for ee2=1:n_e
                        cols=aplow+n_a*(zz2-1)+n_a*n_z*(ee2-1);
                        Tjj(rows+(cols-1)*N_s)=Tjj(rows+(cols-1)*N_s)+pi_z_cpu(zz,zz2)*pi_e_cpu(ee2);
                    end
                end
            end
        end
        T_solo.(nm){jj}=gpuArray(Tjj);
    end
end

% Per type vs solo, every field (unrestricted and restricted)
for ii=1:N_i
    nm=Names_i{ii};
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        statnames=fieldnames(ACP_solo.(nm).(fn));
        for ss=1:length(statnames)
            x=gather(ACP_PT.(fn).(nm).(statnames{ss})(:)); y=gather(ACP_solo.(nm).(fn).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force B, %s %s (type %s) PType vs solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,statnames{ss},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
            x=gather(ACP_PT_ng.(fn).(nm).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force B, %s %s (type %s) groupptypesforstats=0 vs solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,statnames{ss},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
        end
        for rr=1:length(RestNames)
            rn=RestNames{rr};
            statnames=fieldnames(ACP_solo.(nm).(rn).(fn));
            for ss=1:length(statnames)
                x=gather(ACP_PT.(rn).(fn).(nm).(statnames{ss})(:)); y=gather(ACP_solo.(nm).(rn).(fn).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
                fprintf('AutoCorr brute force B, %s %s %s (type %s) PType vs solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',rn,fn,statnames{ss},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
            end
        end
    end
    for rr=1:length(RestNames)
        rn=RestNames{rr};
        fprintf('AutoCorr brute force B, %s RestrictedSampleMass (type %s) PType vs solo, this should be zero: %.3e \n',rn,nm,max(abs(gather(ACP_PT.(rn).RestrictedSampleMass.(nm)(:))-gather(ACP_solo.(nm).(rn).RestrictedSampleMass(:)))))
    end
end
for ff=1:length(FnNames)
    fprintf('AutoCorr brute force B, %s groupptypesforstats=0 has no grouped Mean, this should be zero: %i \n',FnNames{ff},isfield(ACP_PT_ng.(FnNames{ff}),'Mean'))
end

% Grouped vs by-hand pooling of the two types' brute-force joint distributions
Dist_u_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Dist_u_solo.(nm)=gpuArray(reshape(gather(Dist_solo.(nm)),[N_s,N_j]));
end
for ff=1:length(FnNames)
    fn=FnNames{ff};
    % Grouped Mean and StdDeviation at each age, by hand: weight ptw(ii)*(mass of type ii at age jj)
    MeanG_bf=nan(1,N_j); StdG_bf=nan(1,N_j);
    MeanRG_bf=nan(length(RestNames),N_j); StdRG_bf=nan(length(RestNames),N_j);
    for jj=1:N_j
        wmass=zeros(N_i,1); mx=zeros(N_i,1); mxx=zeros(N_i,1);
        wmassR=zeros(N_i,length(RestNames)); mxR=zeros(N_i,length(RestNames)); mxxR=zeros(N_i,length(RestNames));
        for ii=1:N_i
            nm=Names_i{ii};
            x=Vals_solo.(nm){ff}(:,jj);
            wmass(ii)=gather(ptw(ii)*sum(Dist_u_solo.(nm)(:,jj)));
            mx(ii)=gather(ptw(ii)*sum(Dist_u_solo.(nm)(:,jj).*x));
            mxx(ii)=gather(ptw(ii)*sum(Dist_u_solo.(nm)(:,jj).*x.^2));
            for rr=1:length(RestNames)
                mr=Dist_u_solo.(nm)(:,jj).*Rvals_solo.(nm){rr}(:,jj);
                wmassR(ii,rr)=gather(ptw(ii)*sum(mr));
                mxR(ii,rr)=gather(ptw(ii)*sum(mr.*x));
                mxxR(ii,rr)=gather(ptw(ii)*sum(mr.*x.^2));
            end
        end
        if sum(wmass)>0
            MeanG_bf(jj)=sum(mx)/sum(wmass);
            StdG_bf(jj)=sqrt(max(sum(mxx)/sum(wmass)-MeanG_bf(jj)^2,0));
        end
        for rr=1:length(RestNames)
            if sum(wmassR(:,rr))>0
                MeanRG_bf(rr,jj)=sum(mxR(:,rr))/sum(wmassR(:,rr));
                StdRG_bf(rr,jj)=sqrt(max(sum(mxxR(:,rr))/sum(wmassR(:,rr))-MeanRG_bf(rr,jj)^2,0));
            end
        end
    end
    x=gather(ACP_PT.(fn).Mean(:)); y=MeanG_bf(:); dd=abs(x-y); dd=dd(~isnan(dd));
    fprintf('AutoCorr brute force B, %s grouped Mean vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,max([dd;0]),sum(isnan(x)~=isnan(y)))
    x=gather(ACP_PT.(fn).StdDeviation(:)); y=StdG_bf(:); dd=abs(x-y); dd=dd(~isnan(dd));
    fprintf('AutoCorr brute force B, %s grouped StdDeviation vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,max([dd;0]),sum(isnan(x)~=isnan(y)))
    for rr=1:length(RestNames)
        rn=RestNames{rr};
        x=gather(ACP_PT.(rn).(fn).Mean(:)); y=MeanRG_bf(rr,:)'; dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force B, %s %s grouped Mean vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',rn,fn,max([dd;0]),sum(isnan(x)~=isnan(y)))
        x=gather(ACP_PT.(rn).(fn).StdDeviation(:)); y=StdRG_bf(rr,:)'; dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force B, %s %s grouped StdDeviation vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',rn,fn,max([dd;0]),sum(isnan(x)~=isnan(y)))
    end

    % Grouped pairs at each horizon, by hand: pool the types' unnormalised joint distributions with weight ptw(ii)
    for hh=1:length(horizons_PT)
        kk=horizons_PT(hh);
        AutoCovG_bf=nan(1,N_j-kk); AutoCorrG_bf=nan(1,N_j-kk);
        AutoCovRG_bf=nan(length(RestNames),N_j-kk); AutoCorrRG_bf=nan(length(RestNames),N_j-kk); PairMassG_bf=nan(length(RestNames),N_j-kk);
        PairMeanxG_bf=nan(length(RestNames),N_j-kk); PairMeanyG_bf=nan(length(RestNames),N_j-kk); PairStdxG_bf=nan(length(RestNames),N_j-kk); PairStdyG_bf=nan(length(RestNames),N_j-kk);
        for j0=1:N_j-kk
            Ju=cell(N_i,1); xs=cell(N_i,1); ys=cell(N_i,1);
            for ii=1:N_i
                nm=Names_i{ii};
                Pk=eye(N_s,'gpuArray');
                for s2=0:kk-1
                    Pk=Pk*T_solo.(nm){j0+s2};
                end
                Ju{ii}=ptw(ii)*(Dist_u_solo.(nm)(:,j0).*Pk); % unnormalised joint of type ii, weighted by its ptweight
                xs{ii}=Vals_solo.(nm){ff}(:,j0);
                ys{ii}=Vals_solo.(nm){ff}(:,j0+kk);
            end
            % unrestricted
            total=0; sx=0; sy=0;
            for ii=1:N_i
                total=total+sum(Ju{ii}(:)); sx=sx+sum(sum(Ju{ii},2).*xs{ii}); sy=sy+sum(sum(Ju{ii},1)'.*ys{ii});
            end
            if total>0
                mux=sx/total; muy=sy/total;
                vx=0; vy=0; cv=0;
                for ii=1:N_i
                    vx=vx+sum(sum(Ju{ii},2).*(xs{ii}-mux).^2); vy=vy+sum(sum(Ju{ii},1)'.*(ys{ii}-muy).^2);
                    cv=cv+sum(sum(Ju{ii}.*((xs{ii}-mux)*(ys{ii}-muy)')));
                end
                AutoCovG_bf(j0)=gather(cv/total);
                if sqrt(vx*vy)/total>1e-15
                    AutoCorrG_bf(j0)=gather(cv/sqrt(vx*vy));
                end
            end
            % restricted
            for rr=1:length(RestNames)
                total=0; sx=0; sy=0; Jr=cell(N_i,1);
                for ii=1:N_i
                    nm=Names_i{ii};
                    Jr{ii}=Ju{ii}.*(Rvals_solo.(nm){rr}(:,j0)*Rvals_solo.(nm){rr}(:,j0+kk)');
                    total=total+sum(Jr{ii}(:)); sx=sx+sum(sum(Jr{ii},2).*xs{ii}); sy=sy+sum(sum(Jr{ii},1)'.*ys{ii});
                end
                PairMassG_bf(rr,j0)=gather(total);
                if total>0
                    mux=sx/total; muy=sy/total;
                    vx=0; vy=0; cv=0;
                    for ii=1:N_i
                        vx=vx+sum(sum(Jr{ii},2).*(xs{ii}-mux).^2); vy=vy+sum(sum(Jr{ii},1)'.*(ys{ii}-muy).^2);
                        cv=cv+sum(sum(Jr{ii}.*((xs{ii}-mux)*(ys{ii}-muy)')));
                    end
                    PairMeanxG_bf(rr,j0)=gather(mux); PairMeanyG_bf(rr,j0)=gather(muy);
                    PairStdxG_bf(rr,j0)=gather(sqrt(vx/total)); PairStdyG_bf(rr,j0)=gather(sqrt(vy/total));
                    AutoCovRG_bf(rr,j0)=gather(cv/total);
                    if sqrt(vx*vy)/total>1e-15
                        AutoCorrRG_bf(rr,j0)=gather(cv/sqrt(vx*vy));
                    end
                end
            end
        end
        x=gather(ACP_PT.(fn).(['AutoCovariance',hstr_PT{hh}])(:)); y=AutoCovG_bf(:); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force B, %s grouped AutoCovariance%s vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,hstr_PT{hh},max([dd;0]),sum(isnan(x)~=isnan(y)))
        x=gather(ACP_PT.(fn).(['AutoCorrelation',hstr_PT{hh}])(:)); y=AutoCorrG_bf(:); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force B, %s grouped AutoCorrelation%s vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,hstr_PT{hh},max([dd;0]),sum(isnan(x)~=isnan(y)))
        for rr=1:length(RestNames)
            rn=RestNames{rr};
            bfnames={'AutoCovariance','AutoCorrelation','PairMass','PairMean_j','PairMean_jplusk','PairStdDeviation_j','PairStdDeviation_jplusk'};
            bfvals={AutoCovRG_bf,AutoCorrRG_bf,PairMassG_bf,PairMeanxG_bf,PairMeanyG_bf,PairStdxG_bf,PairStdyG_bf};
            for bb=1:length(bfnames)
                x=gather(ACP_PT.(rn).(fn).([bfnames{bb},hstr_PT{hh}])(:)); y=bfvals{bb}(rr,:)'; dd=abs(x-y); dd=dd(~isnan(dd));
                fprintf('AutoCorr brute force B, %s %s grouped %s%s vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',rn,fn,bfnames{bb},hstr_PT{hh},max([dd;0]),sum(isnan(x)~=isnan(y)))
            end
        end
    end
end
for rr=1:length(RestNames)
    rn=RestNames{rr};
    ByAge_bf=zeros(1,N_j);
    for ii=1:N_i
        nm=Names_i{ii};
        ByAge_bf=ByAge_bf+ptw(ii)*gather(sum(Dist_u_solo.(nm).*Rvals_solo.(nm){rr},1));
    end
    fprintf('AutoCorr brute force B, %s RestrictedSampleMass.ByAge vs by hand, this should be zero: %.3e \n',rn,max(abs(gather(ACP_PT.(rn).RestrictedSampleMass.ByAge(:))-ByAge_bf(:))))
    fprintf('AutoCorr brute force B, %s RestrictedSampleMass.Total vs by hand, this should be zero: %.3e \n',rn,abs(gather(ACP_PT.(rn).RestrictedSampleMass.Total)-sum(ByAge_bf)))
end

% A type of zero weight: grouped equals the other type alone
Params_PT0=Params_PT;
Params_PT0.(PTypeDistParamNames{1})=[0;1];
StationaryDist_PT0=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_PT0,simoptions);
ACP_PT0=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_PT0,Policy_PT,FnsToEvaluate_PT,Params_PT0,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_PT);
for ff=1:length(FnNames)
    fn=FnNames{ff};
    statnames={'Mean','StdDeviation'};
    for hh=1:length(horizons_PT)
        statnames=[statnames,{['AutoCovariance',hstr_PT{hh}],['AutoCorrelation',hstr_PT{hh}]}]; %#ok<AGROW>
    end
    for ss=1:length(statnames)
        x=gather(ACP_PT0.(fn).(statnames{ss})(:)); y=gather(ACP_solo.high.(fn).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force B, zero-weight type, %s grouped %s vs the other type alone, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,statnames{ss},max([dd;0]),sum(isnan(x)~=isnan(y)))
    end
    for rr=1:length(RestNames)
        rn=RestNames{rr};
        statnames=fieldnames(ACP_solo.high.(rn).(fn));
        for ss=1:length(statnames)
            x=gather(ACP_PT0.(rn).(fn).(statnames{ss})(:)); y=gather(ACP_solo.high.(rn).(fn).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force B, zero-weight type, %s %s grouped %s vs the other type alone, this should be zero: %.3e, NaN-pattern mismatches %i \n',rn,fn,statnames{ss},max([dd;0]),sum(isnan(x)~=isnan(y)))
        end
    end
end

output=struct();

end
