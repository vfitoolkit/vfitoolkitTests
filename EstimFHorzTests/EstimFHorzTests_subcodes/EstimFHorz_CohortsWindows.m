function output=EstimFHorz_CohortsWindows(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 10: cohort AutoCorr targets (TargetMoments.AutoCorr.cohortC...) of the window and log rows (added 2026-10-14),
% fnname_log, fnname_MperiodAverage and fnname_MperiodAverage_log, in EstimateLifeCycleModel_MethodOfMoments with
% estimoptions.cohortagejshifter. SetupTargetMoments_FHorz_withCohorts reads the rows off the names, and the cohort objective sets
% simoptions.ntimeperiods, .ntimeperiodstype and .logwhichcombos for each cohort's AutoCorr call.
%   (a) two cohorts (entry ages 1 and 5): at the truth (skipestimation=1) the objective is zero and the current moments equal targets
%       built from each cohort's agent distribution with the command at the level and at the log (exact)
%   (b) the cutoff of cohort 2's agent distribution at its last targeted age, which for an M-period window is M-1 ages after its
%       start: identical J and moments to forcing N_j with estimoptions.cohortlastage (exact)
%   (c) must error: MperiodAverage and MperiodSum rows together in a cohort; a statistic and its log at the same start age
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.inc=@(aprime,a,z,w,kappa_j) 1+a+w*kappa_j*z; % positive at every state
FnsInc.inc=FnsToEvaluate.inc;
Jr=Params.Jr;
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
estimoptions=struct(); estimoptions.verbose=0;
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
    simoptions_L=simoptions_c; simoptions_L.ntimeperiods=[1,2]; simoptions_L.timehorizons=2;
    ACL=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_c,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_L); % the level rows
    simoptions_G=simoptions_L; simoptions_G.logwhichcombos=1;
    ACG=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_c,Policy,FnsInc,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_G); % the log rows of inc
    if cc==1
        t=gather(ACL.earnings_2periodAverage.AutoCorrelation_horizon2); t(Jr-3:end)=NaN; % consecutive 2-period windows while working
        TM.AutoCorr.cohort1.earnings_2periodAverage.AutoCorrelation_horizon2=t;
        hand=[hand; t(:)]; %#ok<AGROW>
    else
        t1=gather(ACG.inc_2periodAverage_log.Mean); t1(1:entryages(cc)-1)=NaN; t1(N_j-4:end)=NaN; % the last targeted window starts at N_j-5 (so ends at N_j-4)
        t2=gather(ACG.inc_log.AutoCorrelation_horizon2); t2(1:entryages(cc)-1)=NaN; t2(13:end)=NaN;
        TM.AutoCorr.cohort2.inc_2periodAverage_log.Mean=t1;
        TM.AutoCorr.cohort2.inc_log.AutoCorrelation_horizon2=t2;
        hand=[hand; t1(:); t2(:)]; %#ok<AGROW>
    end
end
hand=hand(~isnan(hand));
nT=length(hand); W=eye(nT); Omega=diag(0.01*(1+(1:nT)'/nT));
EstimParamNames={'beta','pension'};

%% (a) at the truth
estimoptions_a=estimoptions; estimoptions_a.skipestimation=1; estimoptions_a.cohortagejshifter=entryages;
simoptions_a=simoptions; simoptions_a.agemass_withCohort=mewjc;
ran=0; objv=NaN; dmom=NaN; es_a=[];
try
    [~,~,es_a]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_a,vfoptions,simoptions_a);
    ran=1; objv=es_a.objectivefnval;
    if length(es_a.currentmomentvec)==nT
        dmom=CoreMaxAbs(es_a.currentmomentvec(:)-hand);
    else
        dmom=abs(length(es_a.currentmomentvec)-nT);
    end
catch err
    fprintf('CohortsWindows (a): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('CohortsWindows (a) two cohorts, targets of 2-period windows (level and log) and of a log, %i target entries, at the truth: did not run, this should be zero: %i; objective, this should be zero: %.3e; current moments vs the hand-built ones, this should be zero: %.3e \n',nT,1-ran,objv,dmom)

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
    fprintf('CohortsWindows (b): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('CohortsWindows (b) derived cutoffs vs both cohorts forced to N_j: did not run, this should be zero: %i; J, this should be zero: %.3e; current moments, this should be zero: %.3e \n',1-ran,dJ,dmom)

%% (c) must error
errcases={'MperiodAverage and MperiodSum rows together in a cohort','a statistic and its log at the same start age'};
errmsgs={'both MperiodAverage and MperiodSum','both a statistic of'};
for ee=1:2
    TMbad=TM;
    if ee==1
        TMbad.AutoCorr.cohort2.inc_2periodSum.Mean=TM.AutoCorr.cohort2.inc_2periodAverage_log.Mean;
    else
        TMbad.AutoCorr.cohort2.inc_2periodAverage.Mean=TM.AutoCorr.cohort2.inc_2periodAverage_log.Mean; % the level where the log is targeted
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
    fprintf('CohortsWindows (c) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok);
end

output.hand=hand;
end
