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
        TM.AgeConditionalStats.cohort1.assets.Mean=t;
        All_c=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist_c,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_cfull);
        TM.AllStats.cohort1.assets.Gini=gather(All_c.assets.Gini);
        hand=[hand; TM.AllStats.cohort1.assets.Gini; t(:)]; % within a cohort: AllStats, then AgeConditionalStats, then AutoCorr
    else
        TM.AgeConditionalStats.cohort2.assets.Mean=t;
        AC_c=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_c,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_c);
        ta=gather(AC_c.earnings.AutoCorrelation); ta(1:entryages(cc)-1)=NaN; ta(Jr-1:end)=NaN; % (pairs starting before entry have no mass; earnings are zero in retirement)
        TM.AutoCorr.cohort2.earnings.AutoCorrelation=ta;
        hand=[hand; t(:); ta(:)];
    end
end
hand=hand(~isnan(hand));
nT=length(hand); W=eye(nT); Omega=diag(0.01*(1+(1:nT)'/nT));
EstimParamNames={'beta','pension'};
truth=[Params.beta;Params.pension];

%% (a) auto-detected cohorts, (b) the explicit vector
res=cell(2,1);
for form=1:2
    estimoptions_t=estimoptions; estimoptions_t.skipestimation=1; simoptions_t=simoptions;
    if form==1
        estimoptions_t.cohortagejshifter=1; lab='cohortagejshifter=1 (entry ages from the mass), default age weights (with a warning)';
    else
        estimoptions_t.cohortagejshifter=entryages; simoptions_t.agemass_withCohort=mewjc; lab='cohortagejshifter=[1,5], simoptions.agemass_withCohort given (the same weights)';
    end
    ran=0; objv=NaN; dmom=NaN;
    try
        [~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_t,vfoptions,simoptions_t);
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

output=struct();
output.TM=TM;

end
