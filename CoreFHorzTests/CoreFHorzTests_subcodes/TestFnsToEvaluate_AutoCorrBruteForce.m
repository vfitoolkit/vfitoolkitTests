function output=TestFnsToEvaluate_AutoCorrBruteForce(n_a,n_z,n_e,N_j,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,AgeWeightParamNames)
% Fourth config of TestFnsToEvaluate: no d, with z, with e. EvalFnOnAgentDist_AutoCorrTransProbs_FHorz with
% simoptions.timehorizons and simoptions.conditionalrestrictions (added 2026-09-29) against a brute-force calculation
% written here independently of the toolkit's transition-matrix code. (The PType wrapper is tested the same way in
% CoreFHorzPTypeTests part 21.) Every check also feeds fail_count, as in the other TestFnsToEvaluate subcodes.
%
% Brute force: for each age jj the one-period transition kernel T_jj (N_states x N_states, dense) is built by loops
% from the Policy indexes and pi_z/pi_e (with the two-point grid-interpolation weights when gridinterplayer=1), the
% k-period kernel is the product T_j0*...*T_{j0+k-1}, and the joint distribution of (state at j0, state at j0+k) is
% dist_j0(s)*Pk(s,s'). Every moment is then a plain sum over that joint distribution: means, std devs and the
% covariance of x_j0 and x_{j0+k}; under a restriction the joint distribution is masked to the pairs that satisfy
% it at both ages and renormalised (the pair population), which is the definition the command documents.
%
% Model: no d, markov z, iid e (state (a,z,e), N_states=n_a*n_z*n_e). Checks, all exact identities:
%  A1. horizons 1,2,3,5, unrestricted and restricted (rich: a>0.5; highz: z>1), vs the brute force: AutoCovariance,
%      AutoCorrelation and, restricted, PairMass, PairMean_j, PairMean_jplusk, PairStdDeviation_j,
%      PairStdDeviation_jplusk, Mean, StdDeviation (PairMass is zero when nobody at the start age satisfies the
%      restriction, NaN only when the age has no mass).
%  A2. lowmemory=1 equals lowmemory=0 (every field).
%  A3. gridinterplayer=1: the same brute-force checks with the two-point kernel.
%  A4. a cohort distribution started at age 10 (simoptions.jequaloneDistAge=10, from the age-10 slice of the full
%      run): NaN at every age/pair before 10, equal to the full run from 10 on, no error.
% Zero checks are printed as %.3e; NaN is compared as a pattern (count of entries NaN in one but not the other).

fprintf('\n========== TestFnsToEvaluate_AutoCorrBruteForce ==========\n')

% Counter incremented for every test that exceeds its tolerance (reported at end of subcode)
fail_count=0;
TOL_EXACT=1e-10;
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

%% A1-A3. Brute force, without and with the grid interpolation layer
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
    fprintf('AutoCorr brute force (%s), kernels row sums minus one, this should be zero: %.3e \n',cfgstr,max(cellfun(@(M) gather(max(abs(sum(M,2)-1))),T))); fail_count=fail_count+(max(cellfun(@(M) gather(max(abs(sum(M,2)-1))),T))>TOL_EXACT);

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
            fprintf('AutoCorr brute force A1 (%s), %s %s Mean vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,RestNames{rr},FnNames{ff},max([dd;0]),sum(isnan(x)~=isnan(y))); fail_count=fail_count+(max([dd;0])>TOL_EXACT)+(sum(isnan(x)~=isnan(y))>0);
            x=gather(ACP.(RestNames{rr}).(FnNames{ff}).StdDeviation(:)); y=StdR_bf(:); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force A1 (%s), %s %s StdDeviation vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,RestNames{rr},FnNames{ff},max([dd;0]),sum(isnan(x)~=isnan(y))); fail_count=fail_count+(max([dd;0])>TOL_EXACT)+(sum(isnan(x)~=isnan(y))>0);
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
            fprintf('AutoCorr brute force A1 (%s), %s AutoCovariance%s vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,fn,hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(y))); fail_count=fail_count+(max([dd;0])>TOL_EXACT)+(sum(isnan(x)~=isnan(y))>0);
            x=gather(ACP.(fn).(['AutoCorrelation',hstr{hh}])(:)); y=reshape(AutoCorr_bf(ff,hh,1:N_j-kk),[],1); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force A1 (%s), %s AutoCorrelation%s vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,fn,hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(y))); fail_count=fail_count+(max([dd;0])>TOL_EXACT)+(sum(isnan(x)~=isnan(y))>0);
            for rr=1:nR
                rn=RestNames{rr};
                bfnames={'AutoCovariance','AutoCorrelation','PairMass','PairMean_j','PairMean_jplusk','PairStdDeviation_j','PairStdDeviation_jplusk'};
                bfvals={AutoCovR_bf,AutoCorrR_bf,PairMass_bf,PairMeanx_bf,PairMeany_bf,PairStdx_bf,PairStdy_bf};
                for bb=1:length(bfnames)
                    x=gather(ACP.(rn).(fn).([bfnames{bb},hstr{hh}])(:)); y=reshape(bfvals{bb}(ff,hh,1:N_j-kk,rr),[],1); dd=abs(x-y); dd=dd(~isnan(dd));
                    fprintf('AutoCorr brute force A1 (%s), %s %s %s%s vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,rn,fn,bfnames{bb},hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(y))); fail_count=fail_count+(max([dd;0])>TOL_EXACT)+(sum(isnan(x)~=isnan(y))>0);
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
            fprintf('AutoCorr brute force A2 (%s), %s %s lowmemory=1 vs =0, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,fn,statnames{ss},max([dd;0]),sum(isnan(x)~=isnan(y))); fail_count=fail_count+(max([dd;0])>TOL_EXACT)+(sum(isnan(x)~=isnan(y))>0);
        end
        for rr=1:nR
            rn=RestNames{rr};
            statnames=fieldnames(ACP.(rn).(fn));
            for ss=1:length(statnames)
                x=gather(ACP1.(rn).(fn).(statnames{ss})(:)); y=gather(ACP.(rn).(fn).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
                fprintf('AutoCorr brute force A2 (%s), %s %s %s lowmemory=1 vs =0, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,rn,fn,statnames{ss},max([dd;0]),sum(isnan(x)~=isnan(y))); fail_count=fail_count+(max([dd;0])>TOL_EXACT)+(sum(isnan(x)~=isnan(y))>0);
            end
        end
    end
    for rr=1:nR
        rn=RestNames{rr};
        fprintf('AutoCorr brute force A2 (%s), %s RestrictedSampleMass lowmemory=1 vs =0, this should be zero: %.3e \n',cfgstr,rn,max(abs(gather(ACP1.(rn).RestrictedSampleMass(:))-gather(ACP.(rn).RestrictedSampleMass(:))))); fail_count=fail_count+(max(abs(gather(ACP1.(rn).RestrictedSampleMass(:))-gather(ACP.(rn).RestrictedSampleMass(:))))>TOL_EXACT);
        fprintf('AutoCorr brute force A1 (%s), %s RestrictedSampleMass vs brute force, this should be zero: %.3e \n',cfgstr,rn,max(abs(gather(ACP.(rn).RestrictedSampleMass(:))-gather(sum(Dist_u.*Rvals{rr},1))'))); fail_count=fail_count+(max(abs(gather(ACP.(rn).RestrictedSampleMass(:))-gather(sum(Dist_u.*Rvals{rr},1))'))>TOL_EXACT);
    end

    if cfg==1
        ACP_full=ACP; % kept for A4
    end
end

%% A4. Cohort distribution started at age 10: NaN before, equal to the full run from age 10 on
j0c=10;
jequaloneDist10=reshape(gather(StationaryDist),[N_s,N_j]);
jequaloneDist10=gpuArray(reshape(jequaloneDist10(:,j0c)/sum(jequaloneDist10(:,j0c)),[n_a,n_z,n_e]));
simoptions_c10=simoptions;
simoptions_c10.jequaloneDistAge=j0c;
Dist10=StationaryDist_FHorz_Case1(jequaloneDist10,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_c10);
simoptions_c10.timehorizons=timehorizons;
simoptions_c10.conditionalrestrictions=Restrictions;
ACP10=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist10,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_c10);
fprintf('AutoCorr brute force A4, cohort from age %i, Dist mass at ages 1..%i, this should be zero: %.3e \n',j0c,j0c-1,sum(sum(reshape(gather(Dist10),[N_s,N_j]),1)*[ones(j0c-1,1);zeros(N_j-j0c+1,1)])); fail_count=fail_count+(sum(sum(reshape(gather(Dist10),[N_s,N_j]),1)*[ones(j0c-1,1);zeros(N_j-j0c+1,1)])>TOL_EXACT);
for ff=1:length(FnNames)
    fn=FnNames{ff};
    fprintf('AutoCorr brute force A4, cohort from age %i, %s Mean non-NaN entries before age %i, this should be zero: %i \n',j0c,fn,j0c,sum(~isnan(gather(ACP10.(fn).Mean(1:j0c-1))))); fail_count=fail_count+(sum(~isnan(gather(ACP10.(fn).Mean(1:j0c-1))))>0);
    fprintf('AutoCorr brute force A4, cohort from age %i, %s Mean from age %i on vs full run, this should be zero: %.3e \n',j0c,fn,j0c,max(abs(gather(ACP10.(fn).Mean(j0c:end))-gather(ACP_full.(fn).Mean(j0c:end))))); fail_count=fail_count+(max(abs(gather(ACP10.(fn).Mean(j0c:end))-gather(ACP_full.(fn).Mean(j0c:end))))>TOL_EXACT);
    for hh=1:length(horizons)
        kk=horizons(hh);
        fprintf('AutoCorr brute force A4, cohort from age %i, %s AutoCovariance%s non-NaN entries before age %i, this should be zero: %i \n',j0c,fn,hstr{hh},j0c,sum(~isnan(gather(ACP10.(fn).(['AutoCovariance',hstr{hh}])(1:j0c-1))))); fail_count=fail_count+(sum(~isnan(gather(ACP10.(fn).(['AutoCovariance',hstr{hh}])(1:j0c-1))))>0);
        x=gather(ACP10.(fn).(['AutoCovariance',hstr{hh}])(j0c:N_j-kk)); y=gather(ACP_full.(fn).(['AutoCovariance',hstr{hh}])(j0c:N_j-kk)); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force A4, cohort from age %i, %s AutoCovariance%s from age %i on vs full run, this should be zero: %.3e, NaN-pattern mismatches %i \n',j0c,fn,hstr{hh},j0c,max([dd,0]),sum(isnan(x)~=isnan(y))); fail_count=fail_count+(max([dd,0])>TOL_EXACT)+(sum(isnan(x)~=isnan(y))>0);
        x=gather(ACP10.rich.(fn).(['AutoCovariance',hstr{hh}])(j0c:N_j-kk)); y=gather(ACP_full.rich.(fn).(['AutoCovariance',hstr{hh}])(j0c:N_j-kk)); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force A4, cohort from age %i, rich %s AutoCovariance%s from age %i on vs full run, this should be zero: %.3e, NaN-pattern mismatches %i \n',j0c,fn,hstr{hh},j0c,max([dd,0]),sum(isnan(x)~=isnan(y))); fail_count=fail_count+(max([dd,0])>TOL_EXACT)+(sum(isnan(x)~=isnan(y))>0);
        fprintf('AutoCorr brute force A4, cohort from age %i, rich %s PairMass%s non-NaN entries before age %i, this should be zero: %i \n',j0c,fn,hstr{hh},j0c,sum(~isnan(gather(ACP10.rich.(fn).(['PairMass',hstr{hh}])(1:j0c-1))))); fail_count=fail_count+(sum(~isnan(gather(ACP10.rich.(fn).(['PairMass',hstr{hh}])(1:j0c-1))))>0);
    end
end


%% Summary for this subcode
fprintf('\nTestFnsToEvaluate_AutoCorrBruteForce: %d tests outside tolerance.\n', fail_count)
output=struct();
output.fail_count=fail_count;
end
