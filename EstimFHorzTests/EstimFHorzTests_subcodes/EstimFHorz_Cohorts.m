function output=EstimFHorz_Cohorts(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 6: cohorts entering at different ages (estimoptions.cohortagejshifter) in EstimateLifeCycleModel_MethodOfMoments.
% jequaloneDist carries an extra last dimension of N_j: slice j is the cohort entering at age j. Targets are by cohort
% (TargetMoments.AgeConditionalStats.cohortC.fn.stat, .AllStats.cohortC.fn.stat, .AutoCorr.cohortC.fn.stat), each cohort's moments
% being the statistics of its own agent distribution (StationaryDist_FHorz_Case1 with simoptions.jequaloneDistAge at the entry age,
% age weights mewj(j) at age j, the population's, unless simoptions.agemass_withCohort says otherwise; the drivers warn when that default is
% used. Until 2026-10-07 the default was mewj(entry age) flat at every age, which this bank's uniform mewj could not distinguish; the PType
% bank's geometric mewj did, through its check (c), and the default was changed so that one cohort at age 1 IS the standard estimator).
%   (a) two cohorts (entry ages 1 and 5, masses 0.6 and 0.4): at the truth (skipestimation=1) the objective is zero and the current
%       moments equal targets the test builds by that rule from the agent-distribution and stats commands (exact)
%   (b) cohortagejshifter given as the vector of entry ages [1,5] and simoptions.agemass_withCohort given explicitly (the same weights as the default): identical (exact)
%   (c) one cohort entering at age 1 with all the mass: J and current moments equal the standard estimator's on the same targets (exact)
%   (d) started away (beta 0.90, pension 0.30), lsqnonlin: recovered (close)
%   (e) must error: the new target kinds with cohorts; logmoments by name with cohorts; a jequaloneDist without the cohort dimension
%   (f) the cutoff of each cohort's agent distribution at its last targeted age (simoptions.lastDistAge, derived by the estimator; cohort 2's
%       targets stop four ages before N_j): identical J and moments to forcing N_j with estimoptions.cohortlastage (exact); a cutoff before the entry age must error
%   (g) simoptions.lastDistAge in StationaryDist_FHorz_Case1 itself (entry age 5, and the standard entry age 1): equal to the full
%       distribution up to the cutoff, zero after it, mass = the age weights up to the cutoff (exact); larger than N_j must error
%   (h) a matrix target, QuantileCutoffs (90th-percentile row), for the cohort entering at age 5: at the truth the objective is zero and the
%       current moments equal the hand-built ones; the derived cutoff equals forcing N_j (exact); before entry or without N_j columns must error
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
Jr=Params.Jr;
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);
estimoptions=struct(); estimoptions.verbose=0;
entryages=[1,5]; masses=[0.6,0.4];

%% The cohort initial distribution and the hand-built targets at the truth
jequaloneDist_c=zeros(n_a,n_z,N_j,'gpuArray');
jequaloneDist_c(1,ceil(n_z/2),1)=masses(1); % the cohort entering at age 1: no assets, the middle shock
jequaloneDist_c(1,ceil(n_z/2),5)=masses(2); % the cohort entering at age 5
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
mewjc=zeros(2,N_j);
for cc=1:2
    mewjc(cc,:)=Params.mewj; % the default rule (2026-10-07): the population age weights from the entry age on (ages before entry carry no mass)
end
TM=struct(); hand=[];
for cc=1:2
    Params_c=Params; Params_c.mewj=mewjc(cc,:);
    simoptions_c=simoptions; simoptions_c.jequaloneDistAge=entryages(cc); simoptions_c.agemass_withCohort=mewjc;
    simoptions_cfull=simoptions_c; simoptions_cfull.whichstats=ones(1,7);
    slice=reshape(jequaloneDist_c(:,:,entryages(cc))/masses(cc),[n_a,n_z]);
    Dist_c=StationaryDist_FHorz_Case1(slice,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params_c,simoptions_c);
    LCP_c=LifeCycleProfiles_FHorz_Case1(Dist_c,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_cfull);
    t=gather(LCP_c.assets.Mean); t(1:entryages(cc))=NaN; t(t<1e-6)=NaN; % (no mass before entry, no assets at entry)
    if cc==1
        TM.AgeConditionalStats.cohort1.assets.Mean=t; t1=t; % t1 kept for check (h)
        All_c=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist_c,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_cfull);
        TM.AllStats.cohort1.assets.Gini=gather(All_c.assets.Gini);
        hand=[hand; TM.AllStats.cohort1.assets.Gini; t(:)]; % within a cohort: AllStats, then AgeConditionalStats, then AutoCorr
    else
        t(end-3:end)=NaN; % the last four ages are not targeted, so the cohort's distribution need not be computed there (check (f))
        TM.AgeConditionalStats.cohort2.assets.Mean=t;
        AC_c=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_c,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_c);
        ta=gather(AC_c.earnings.AutoCorrelation); ta(1:entryages(cc)-1)=NaN; ta(Jr-1:end)=NaN; % (pairs starting before entry have no mass; earnings are zero in retirement)
        cutoff2=max(find(~isnan(t),1,'last'),find(~isnan(ta),1,'last')+1); % the last age cohort 2's targets need (a horizon-1 pair starting at jj needs age jj+1)
        Dist_c2=Dist_c; LCP_c2=LCP_c; % kept for checks (g) and (h)
        TM.AutoCorr.cohort2.earnings.AutoCorrelation=ta;
        hand=[hand; t(:); ta(:)];
    end
end
hand=hand(~isnan(hand));
nT=length(hand); W=eye(nT); Omega=diag(0.01*(1+(1:nT)'/nT));
EstimParamNames={'beta','pension'};
truth=[Params.beta;Params.pension];

%% (a) auto-detected cohorts, (b) the explicit vector
res=cell(2,1); tform=zeros(2,1);
for form=1:2
    estimoptions_t=estimoptions; estimoptions_t.skipestimation=1; simoptions_t=simoptions;
    if form==1
        estimoptions_t.cohortagejshifter=1; lab='cohortagejshifter=1 (entry ages from the mass), default age weights (with a warning)';
    else
        estimoptions_t.cohortagejshifter=entryages; simoptions_t.agemass_withCohort=mewjc; lab='cohortagejshifter=[1,5], simoptions.agemass_withCohort given (the same weights)';
    end
    ran=0; objv=NaN; dmom=NaN;
    try
        t0=tic;
        [~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_t,vfoptions,simoptions_t);
        tform(form)=toc(t0);
        ran=1; objv=estsummary.objectivefnval; res{form}=estsummary;
        if length(estsummary.currentmomentvec)==nT
            dmom=CoreMaxAbs(estsummary.currentmomentvec(:)-hand);
        else
            dmom=abs(length(estsummary.currentmomentvec)-nT);
        end
    catch err
        fprintf('Cohorts %s: the estimation command errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('Cohorts (%s) two cohorts (entry ages 1 and 5), %s, %i target entries, at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current moments vs the hand-built cohort moments, this should be zero: %.3e \n',char('a'+form-1),lab,nT,1-ran,objv,dmom)
end
if ~isempty(res{1}) && ~isempty(res{2})
    fprintf('Cohorts (b) vector vs auto-detected entry ages: J, this should be zero: %.3e; current moments, this should be zero: %.3e \n',CoreMaxAbs(res{1}.variousmatrices.J(:)-res{2}.variousmatrices.J(:)),CoreMaxAbs(res{1}.currentmomentvec(:)-res{2}.currentmomentvec(:)))
else
    fprintf('Cohorts (b) vector vs auto-detected entry ages: one of the two did not run, this should be zero: %.3e \n',1)
end

%% (f) the cutoff of each cohort's distribution at its last targeted age (simoptions.lastDistAge, derived by SetupTargetMoments_FHorz_withCohorts)
% Cohort 1 has an AllStats target so it runs to N_j; cohort 2's targets end four ages before N_j. Forcing both to N_j (estimoptions.cohortlastage) must give identical J and moments.
estimoptions_f=estimoptions; estimoptions_f.skipestimation=1; estimoptions_f.cohortagejshifter=1; estimoptions_f.cohortlastage=[N_j,N_j];
ran=0; dJ=NaN; dmom=NaN; tf=NaN;
try
    t0=tic;
    [~,~,es_full]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_f,vfoptions,simoptions);
    tf=toc(t0); ran=1;
    if ~isempty(res{1})
        dJ=CoreMaxAbs(res{1}.variousmatrices.J(:)-es_full.variousmatrices.J(:)); dmom=CoreMaxAbs(res{1}.currentmomentvec(:)-es_full.currentmomentvec(:));
    end
catch err
    fprintf('Cohorts (f): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Cohorts (f) cohort 2 cutoff at age %i of %i (derived, %.1f s) vs both cohorts forced to N_j (estimoptions.cohortlastage, %.1f s): did not run, this should be zero: %.3e; cutoff is before N_j, this should be zero: %.3e; J, this should be zero: %.3e; current moments, this should be zero: %.3e \n',cutoff2,N_j,tform(1),tf,1-ran,double(cutoff2>=N_j),dJ,dmom)
estimoptions_f.cohortlastage=[N_j,entryages(2)-1]; % a cutoff before the entry age must error
errored=0;
try
    EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_f,vfoptions,simoptions);
catch
    errored=1;
end
fprintf('Cohorts (f) estimoptions.cohortlastage before a cohort''s entry age must error, this should be zero: %.3e \n',1-errored)

%% (h) a matrix target, QuantileCutoffs ((nquantiles+1) x N_j), for the cohort entering at age 5
% Only the 90th-percentile row at ages 6 to N_j-4. Until 2026-10-15 SetupTargetMoments_FHorz_withCohorts reshaped the stacked age-conditional
% block as (N_j, moments), which misreads a matrix target: it refused this one (entries at later ages fell on ages before entry) and its last age
% came out before N_j-4, which would have cut the cohort's distribution too early. The ages of a target now come from its own shape (age is the column).
Qall=gather(LCP_c2.assets.QuantileCutoffs);
r90=round(0.9*(size(Qall,1)-1))+1; % row 1 is the minimum, row k+1 the k/nquantiles quantile
Q=nan(size(Qall)); Q(r90,entryages(2)+1:N_j-4)=Qall(r90,entryages(2)+1:N_j-4);
temp=reshape(Q(:),N_j,[]); % the old reshape-based reading of this target
oldrefused=any(~isnan(temp(1:entryages(2)-1,:)),'all');
TMh=struct(); TMh.AgeConditionalStats.cohort1.assets.Mean=t1; TMh.AgeConditionalStats.cohort2.assets.QuantileCutoffs=Q;
handh=[t1(:); Q(:)]; handh=handh(~isnan(handh));
nTh=length(handh); Wh=eye(nTh); Omegah=diag(0.01*(1+(1:nTh)'/nTh));
estimoptions_h=estimoptions; estimoptions_h.skipestimation=1; estimoptions_h.cohortagejshifter=entryages;
simoptions_h=simoptions; simoptions_h.agemass_withCohort=mewjc;
ran=0; objv=NaN; dmom=NaN; es_h=[];
try
    [~,~,es_h]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMh,Wh,Omegah,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_h,vfoptions,simoptions_h);
    ran=1; objv=es_h.objectivefnval;
    if length(es_h.currentmomentvec)==nTh
        dmom=CoreMaxAbs(es_h.currentmomentvec(:)-handh);
    else
        dmom=abs(length(es_h.currentmomentvec)-nTh);
    end
catch err
    fprintf('Cohorts (h): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Cohorts (h) a QuantileCutoffs target (row %i, the 90th percentile, ages %i to %i) for the cohort entering at age %i, %i target entries, at the truth: the old reshape-based check refused it, this should be zero: %.3e; did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current moments vs the hand-built cohort moments, this should be zero: %.3e \n',r90,entryages(2)+1,N_j-4,entryages(2),nTh,double(~oldrefused),1-ran,objv,dmom)
estimoptions_h.cohortlastage=[N_j,N_j]; % the derived cutoff of cohort 2 (age N_j-4) vs forcing N_j
ran=0; dJ=NaN; dmom=NaN;
try
    [~,~,es_hfull]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMh,Wh,Omegah,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_h,vfoptions,simoptions_h);
    ran=1;
    if ~isempty(es_h)
        dJ=CoreMaxAbs(es_h.variousmatrices.J(:)-es_hfull.variousmatrices.J(:)); dmom=CoreMaxAbs(es_h.currentmomentvec(:)-es_hfull.currentmomentvec(:));
    end
catch err
    fprintf('Cohorts (h) forced to N_j: the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Cohorts (h) the QuantileCutoffs target with the derived cutoff vs both cohorts forced to N_j (estimoptions.cohortlastage): did not run, this should be zero: %.3e; J, this should be zero: %.3e; current moments, this should be zero: %.3e \n',1-ran,dJ,dmom)
estimoptions_h=rmfield(estimoptions_h,'cohortlastage');
errored=zeros(1,2);
TMbad=TMh; TMbad.AgeConditionalStats.cohort2.assets.QuantileCutoffs(r90,entryages(2)-1)=1; % a matrix target at an age before entry must error
Wb=eye(nTh+1); Omegab=diag(0.01*(1+(1:nTh+1)'/(nTh+1))); % (sized for the extra entry, so only the entry age can be the reason it errors)
try
    EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMbad,Wb,Omegab,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_h,vfoptions,simoptions_h);
catch
    errored(1)=1;
end
TMbad=TMh; TMbad.AgeConditionalStats.cohort2.assets.QuantileCutoffs=Q(:,1:N_j-1); % a matrix target without N_j columns must error
try
    EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMbad,Wh,Omegah,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_h,vfoptions,simoptions_h);
catch
    errored(2)=1;
end
fprintf('Cohorts (h) a QuantileCutoffs target before the entry age, and one without N_j columns, must each error, this should be zero: %.3e \n',2-sum(errored))

%% (c) one cohort at age 1 with all the mass equals the standard estimator
jequaloneDist_1=zeros(n_a,n_z,N_j,'gpuArray'); jequaloneDist_1(:,:,1)=jequaloneDist;
Dist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
LCPFull=LifeCycleProfiles_FHorz_Case1(Dist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
t=gather(LCPFull.assets.Mean); t(1)=NaN; t(t<1e-6)=NaN;
TMs=struct(); TMs.AllStats.assets.Gini=gather(AllFull.assets.Gini); TMs.AgeConditionalStats.assets.Mean=t;
TM1=struct(); TM1.AllStats.cohort1.assets.Gini=TMs.AllStats.assets.Gini; TM1.AgeConditionalStats.cohort1.assets.Mean=t;
nT1=sum(~isnan(t))+1; W1=eye(nT1); Omega1=diag(0.01*(1+(1:nT1)'/nT1));
estimoptions_s=estimoptions; estimoptions_s.skipestimation=1;
estimoptions_1=estimoptions_s; estimoptions_1.cohortagejshifter=1;
ran=0; dJ=NaN; dmom=NaN;
try
    [~,~,es_std]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMs,W1,Omega1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_s,vfoptions,simoptions);
    [~,~,es_coh]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM1,W1,Omega1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_1,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_1,vfoptions,simoptions);
    ran=1; dJ=CoreMaxAbs(es_std.variousmatrices.J(:)-es_coh.variousmatrices.J(:)); dmom=CoreMaxAbs(es_std.currentmomentvec(:)-es_coh.currentmomentvec(:));
catch err
    fprintf('Cohorts (c): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Cohorts (c) one cohort entering at age 1 with all the mass vs the standard estimator: did not run, this should be zero: %.3e; J, this should be zero: %.3e; current moments, this should be zero: %.3e \n',1-ran,dJ,dmom)

%% (d) started away
Params_away=Params; Params_away.beta=0.90; Params_away.pension=0.30;
estimoptions_d=estimoptions; estimoptions_d.fminalgo=8; estimoptions_d.cohortagejshifter=1;
ran=0; perr=NaN; tt=NaN; objv=NaN; rec=[NaN;NaN];
try
    t0=tic;
    [EstimParams,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_d,vfoptions,simoptions);
    tt=toc(t0); ran=1; rec=[EstimParams.beta;EstimParams.pension]; perr=CoreMaxAbs(rec-truth); objv=estsummary.objectivefnval;
catch err
    fprintf('Cohorts (d): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Cohorts (d) two cohorts, started at beta=0.90, pension=0.30 (%.1f s), recovered %.6f, %.6f, objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',tt,rec(1),rec(2),objv,1-ran,perr)

%% (e) must error
errored=0;
TMbad=TM; TMbad.CrossSectionCovarCorr.assets.CorrelationWith.earnings=0.3; % the new target kinds are not implemented with cohorts
try
    EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMbad,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_1,vfoptions,simoptions);
catch
    errored=errored+1;
end
estimoptions_bad=estimoptions_1; estimoptions_bad.logmoments=struct(); estimoptions_bad.logmoments.AllStats.cohort1.assets.Gini=1; % logmoments by name is not implemented with cohorts
try
    EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_bad,vfoptions,simoptions);
catch
    errored=errored+1;
end
try
    EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_1,vfoptions,simoptions); % no cohort dimension
catch
    errored=errored+1;
end
fprintf('Cohorts (e) the new target kinds, logmoments by name, and a jequaloneDist without the cohort dimension must all error with cohorts, this should be zero: %.3e \n',3-errored)

%% (g) simoptions.lastDistAge in the agent distribution itself: equal to the full distribution up to the cutoff, zero after it
% (Params_c, simoptions_c and slice are cohort 2's from the loop above: entry age 5)
L=N_j-4;
simoptions_g=simoptions_c; simoptions_g.lastDistAge=L;
t0=tic; Dist_g=StationaryDist_FHorz_Case1(slice,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params_c,simoptions_g); tg=toc(t0);
t0=tic; StationaryDist_FHorz_Case1(slice,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params_c,simoptions_c); tfull=toc(t0);
Dg=reshape(Dist_g,[n_a*n_z,N_j]); Df=reshape(Dist_c2,[n_a*n_z,N_j]);
fprintf('Cohorts (g) cohort 2 (entry age 5) distribution with simoptions.lastDistAge=%i (%.1f s) vs without (%.1f s): up to the cutoff, this should be zero: %.3e; after the cutoff, this should be zero: %.3e; total mass vs the age weights up to the cutoff, this should be zero: %.3e \n',L,tg,tfull,CoreMaxAbs(Dg(:,1:L)-Df(:,1:L)),CoreMaxAbs(Dg(:,L+1:end)),abs(gather(sum(Dg(:)))-sum(mewjc(2,entryages(2):L))))
% the same with the standard first-period distribution (jequaloneDistAge=1): the cutoff alone routes through the sliced solver
simoptions_g=simoptions; simoptions_g.lastDistAge=L;
t0=tic; Dist_g1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_g); tg1=toc(t0);
Dg=reshape(Dist_g1,[n_a*n_z,N_j]); Df=reshape(Dist,[n_a*n_z,N_j]);
fprintf('Cohorts (g) the standard distribution (entry age 1) with simoptions.lastDistAge=%i (%.1f s): up to the cutoff, this should be zero: %.3e; after the cutoff, this should be zero: %.3e \n',L,tg1,CoreMaxAbs(Dg(:,1:L)-Df(:,1:L)),CoreMaxAbs(Dg(:,L+1:end)))
simoptions_g.lastDistAge=N_j+1; errored=0; % must error
try
    StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_g);
catch
    errored=1;
end
fprintf('Cohorts (g) simoptions.lastDistAge larger than N_j must error, this should be zero: %.3e \n',1-errored)

output=struct();
output.TM=TM;

end
