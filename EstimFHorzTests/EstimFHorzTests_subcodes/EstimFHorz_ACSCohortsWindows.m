function output=EstimFHorz_ACSCohortsWindows(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 11: cohort AgeConditionalStats targets (TargetMoments.AgeConditionalStats.cohortC...) of the window and log rows
% (added 2026-10-14), fnname_log, fnname_MperiodAverage and fnname_MperiodAverage_log, in EstimateLifeCycleModel_MethodOfMoments with
% estimoptions.cohortagejshifter (the route of the De Nardi, French and Jones (2010) medical-expense profiles: the 90th percentile of
% two-period averages, by cohort). SetupTargetMoments_FHorz_withCohorts records each cohort's window lengths and log rows, and the cohort
% objective sets simoptions.ntimeperiods, .ntimeperiodstype, .logwhichcombos and .pi_z for the cohort's LifeCycleProfiles call.
%   (a) two cohorts (entry ages 1 and 5): at the truth (skipestimation=1) the objective is zero and the current moments equal targets
%       built from each cohort's agent distribution with LifeCycleProfiles at the level and at the log (exact)
%   (b) the cutoff of cohort 2's agent distribution at its last targeted age, which for an M-period window is M-1 ages after its
%       start: identical J and moments to forcing N_j with estimoptions.cohortlastage (exact)
%   (c) must error: MperiodAverage and MperiodSum rows together in a cohort; a statistic and its log at the same age
%   (d) estimoptions.whichcombos=0 (every statistic any target names, of every row, at every age, on every restriction page) vs the default =1,
%       under which the cohort's LifeCycleProfiles call computes only the targeted (row, age, restriction) combinations and their statistics
%       (simoptions.whichcombos and a per-combination whichstats from the selectors SetupTargetMoments_FHorz builds per cohort): identical
%       objective, J and moments (exact), both timed (until 2026-10-23 the cohort objective never set simoptions.whichcombos, so =1 silently
%       did the work of =0: every window at every start age on every page, the DFJ2010 slowdown)
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.inc=@(aprime,a,z,w,kappa_j) 1+a+w*kappa_j*z; % positive at every state
FnsInc.inc=FnsToEvaluate.inc;
Jr=Params.Jr;
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
estimoptions=struct(); estimoptions.verbose=0;
estimoptions.eedefault=3; % the runs compared below share one step (the default 'auto' chooses per run)
entryages=[1,5]; masses=[0.6,0.4];

%% The cohort initial distribution and the hand-built targets at the truth
jequaloneDist_c=zeros(n_a,n_z,N_j,'gpuArray');
jequaloneDist_c(1,ceil(n_z/2),1)=masses(1);
jequaloneDist_c(1,ceil(n_z/2),5)=masses(2);
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
mewjc=repmat(Params.mewj,2,1); % the default rule: the population age weights from the entry age on
TM=struct(); hand=[];
for cc=1:2
    Params_c=Params; Params_c.mewj=mewjc(cc,:);
    simoptions_c=simoptions; simoptions_c.jequaloneDistAge=entryages(cc); simoptions_c.agemass_withCohort=mewjc;
    slice=reshape(jequaloneDist_c(:,:,entryages(cc))/masses(cc),[n_a,n_z]);
    Dist_c=StationaryDist_FHorz_Case1(slice,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params_c,simoptions_c);
    simoptions_L=simoptions_c; simoptions_L.ntimeperiods=[1,2]; simoptions_L.pi_z=pi_z; simoptions_L.whichstats=ones(1,7); simoptions_L.warnzerorestrictedmass=0;
    ACL=LifeCycleProfiles_FHorz_Case1(Dist_c,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_L); % the level rows
    simoptions_G=simoptions_L; simoptions_G.logwhichcombos=1;
    ACG=LifeCycleProfiles_FHorz_Case1(Dist_c,Policy,FnsInc,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_G); % the log rows of inc
    if cc==1
        t=gather(ACL.earnings_2periodAverage.MoreInequality.Percentile90th); t(Jr-2:end)=NaN; % the 90th percentile of 2-period averages while working
        TM.AgeConditionalStats.cohort1.earnings_2periodAverage.MoreInequality.Percentile90th=t;
        t2=gather(ACL.earnings.Mean); t2(Jr:end)=NaN;
        TM.AgeConditionalStats.cohort1.earnings.Mean=t2;
        hand=[hand; t(:); t2(:)]; %#ok<AGROW>
    else
        t1=gather(ACG.inc_2periodAverage_log.Mean); t1(1:entryages(cc)-1)=NaN; t1(N_j-4:end)=NaN; % the last targeted window starts at N_j-5 (so ends at N_j-4)
        t2=gather(ACG.rich.inc_log.StdDeviation); t2(1:entryages(cc)-1)=NaN; t2(13:end)=NaN;
        TM.AgeConditionalStats.cohort2.inc_2periodAverage_log.Mean=t1;
        TM.AgeConditionalStats.cohort2.rich.inc_log.StdDeviation=t2;
        hand=[hand; t1(:); t2(:)]; %#ok<AGROW>
    end
end
hand=hand(~isnan(hand));
nT=length(hand); W=eye(nT); Omega=diag(0.01*(1+(1:nT)'/nT));
EstimParamNames={'beta','pension'};

%% (a) at the truth
estimoptions_a=estimoptions; estimoptions_a.skipestimation=1; estimoptions_a.cohortagejshifter=entryages;
simoptions_a=simoptions; simoptions_a.agemass_withCohort=mewjc;
ran=0; objv=NaN; dmom=NaN; es_a=[]; ta=NaN;
try
    t0=tic;
    [~,~,es_a]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_a,vfoptions,simoptions_a);
    ta=toc(t0); ran=1; objv=es_a.objectivefnval;
    if length(es_a.currentmomentvec)==nT
        dmom=CoreMaxAbs(es_a.currentmomentvec(:)-hand);
    else
        dmom=abs(length(es_a.currentmomentvec)-nT);
    end
catch err
    fprintf('ACSCohortsWindows (a): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ACSCohortsWindows (a) two cohorts, AgeConditionalStats targets of a 2-period window (the 90th percentile), of the log of one, and of a log under a restriction, %i target entries, at the truth: did not run, this should be zero: %i; objective, this should be zero: %.3e; current moments vs the hand-built ones, this should be zero: %.3e \n',nT,1-ran,objv,dmom)

%% (b) cohort 2's cutoff includes the M-1 ages of its last window: identical to forcing N_j
estimoptions_b=estimoptions_a; estimoptions_b.cohortlastage=[N_j,N_j];
ran=0; dJ=NaN; dmom=NaN;
try
    [~,~,es_b]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_b,vfoptions,simoptions_a);
    ran=1;
    if ~isempty(es_a)
        dJ=CoreMaxAbs(es_a.variousmatrices.J(:)-es_b.variousmatrices.J(:)); dmom=CoreMaxAbs(es_a.currentmomentvec(:)-es_b.currentmomentvec(:));
    end
catch err
    fprintf('ACSCohortsWindows (b): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ACSCohortsWindows (b) derived cutoffs vs both cohorts forced to N_j: did not run, this should be zero: %i; J, this should be zero: %.3e; current moments, this should be zero: %.3e \n',1-ran,dJ,dmom)

%% (d) estimoptions.whichcombos=0 (every statistic of every row at every age on every page) vs the default =1 (only the targeted combinations: the per-cohort selectors)
estimoptions_d=estimoptions_a; estimoptions_d.whichcombos=0;
ran=0; objv=NaN; dJ=NaN; dmom=NaN; td=NaN;
try
    t0=tic;
    [~,~,es_d]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_d,vfoptions,simoptions_a);
    td=toc(t0); ran=1; objv=es_d.objectivefnval;
    if ~isempty(es_a)
        dJ=CoreMaxAbs(es_a.variousmatrices.J(:)-es_d.variousmatrices.J(:)); dmom=CoreMaxAbs(es_a.currentmomentvec(:)-es_d.currentmomentvec(:));
    end
catch err
    fprintf('ACSCohortsWindows (d): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ACSCohortsWindows (d) estimoptions.whichcombos=0 (%.1f s) vs the default =1 (%.1f s; only the targeted combinations): did not run, this should be zero: %i; objective, this should be zero: %.3e; J, this should be zero: %.3e; current moments, this should be zero: %.3e \n',td,ta,1-ran,objv,dJ,dmom)

%% (c) must error
errcases={'MperiodAverage and MperiodSum rows together in a cohort','a statistic and its log at the same age'};
errmsgs={'both MperiodAverage and MperiodSum','both a statistic'};
for ee=1:2
    TMbad=TM;
    if ee==1
        TMbad.AgeConditionalStats.cohort2.inc_2periodSum.Mean=TM.AgeConditionalStats.cohort2.inc_2periodAverage_log.Mean;
    else
        TMbad.AgeConditionalStats.cohort2.inc_2periodAverage.Mean=TM.AgeConditionalStats.cohort2.inc_2periodAverage_log.Mean; % the level where the log is targeted
    end
    errok=0;
    try
        EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMbad,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_a,vfoptions,simoptions_a);
    catch ME
        errok=contains(ME.message,errmsgs{ee});
        if ~errok
            fprintf('   (the error was: %s) \n',ME.message);
        end
    end
    fprintf('ACSCohortsWindows (c) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok);
end

output.hand=hand;
end
