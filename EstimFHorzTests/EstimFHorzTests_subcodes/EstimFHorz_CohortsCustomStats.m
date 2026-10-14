function output=EstimFHorz_CohortsCustomStats(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 12: TargetMoments.CustomModelStats.cohortC.(name) in EstimateLifeCycleModel_MethodOfMoments with
% estimoptions.cohortagejshifter (added 2026-10-14). The same estimoptions.CustomModelStats function as without cohorts
% (EstimFHorz_CohortsCustomStats_Fn), called once per cohort on that cohort's agent distribution, with simoptions.jequaloneDistAge
% its entry age; its targets are stacked last within the cohort (after AgeConditionalStats and AutoCorr), and the cohort's
% distribution runs to N_j. Its stats include the De Nardi, French and Jones (2010) autocorrelation (the full-cell normalisation).
%   (a) two cohorts (entry ages 1 and 5), CustomModelStats targets in both, with AgeConditionalStats targets in cohort 1 and AutoCorr
%       targets in cohort 2: at the truth (skipestimation=1) the objective is zero and the current moments equal targets built with
%       the same function on each cohort's agent distribution, in the stacking order (exact); entryage checks the cohort is told
%   (b) identical J and moments to forcing both cohorts to N_j with estimoptions.cohortlastage (exact)
%   (c) must error: CustomModelStats targets not under cohortC; without estimoptions.CustomModelStats; a targeted name the
%       function does not return; a target of the wrong length
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.inc=@(aprime,a,z,w,kappa_j) 1+a+w*kappa_j*z; % positive at every state
Jr=Params.Jr;
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
estimoptions=struct(); estimoptions.verbose=0;
estimoptions.CustomModelStats=@EstimFHorz_CohortsCustomStats_Fn;
entryages=[1,5]; masses=[0.6,0.4];

%% The cohort initial distribution and the hand-built targets at the truth
jequaloneDist_c=zeros(n_a,n_z,N_j,'gpuArray');
jequaloneDist_c(1,ceil(n_z/2),1)=masses(1);
jequaloneDist_c(1,ceil(n_z/2),5)=masses(2);
[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
mewjc=repmat(Params.mewj,2,1); % the default rule: the population age weights from the entry age on
TM=struct(); hand=[];
for cc=1:2
    Params_c=Params; Params_c.mewj=mewjc(cc,:);
    simoptions_c=simoptions; simoptions_c.jequaloneDistAge=entryages(cc); simoptions_c.agemass_withCohort=mewjc; simoptions_c.warnzerorestrictedmass=0;
    slice=reshape(jequaloneDist_c(:,:,entryages(cc))/masses(cc),[n_a,n_z]);
    Dist_c=StationaryDist_FHorz_Case1(slice,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params_c,simoptions_c);
    CS=EstimFHorz_CohortsCustomStats_Fn(V,Policy,Dist_c,Params_c,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,estimoptions,vfoptions,simoptions_c);
    cohortstr=['cohort',num2str(cc)];
    if cc==1
        simoptions_L=simoptions_c; simoptions_L.whichstats=[1,0,0,0,0,0,0];
        ACS=LifeCycleProfiles_FHorz_Case1(Dist_c,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_L);
        t=gather(ACS.earnings.Mean); t(Jr:end)=NaN;
        TM.AgeConditionalStats.cohort1.earnings.Mean=t;
        hand=[hand; t(:)]; %#ok<AGROW>
    else
        simoptions_A=simoptions_c; simoptions_A.timehorizons=2;
        ACP=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_c,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_A);
        t=gather(ACP.assets.AutoCorrelation_horizon2); t(1:entryages(cc)-1)=NaN; t(12:end)=NaN;
        TM.AutoCorr.cohort2.assets.AutoCorrelation_horizon2=t;
        hand=[hand; t(:)]; %#ok<AGROW>
    end
    % the CustomModelStats targets, stacked last within the cohort, in the order given
    t1=CS.entryage; TM.CustomModelStats.(cohortstr).entryage=t1;
    t2=CS.hasV; TM.CustomModelStats.(cohortstr).hasV=t2;
    t3=CS.incmean; t3(1:entryages(cc)-1)=NaN; TM.CustomModelStats.(cohortstr).incmean=t3;
    t4=CS.dfjac; t4(1:entryages(cc)-1)=NaN; t4(14:end)=NaN; TM.CustomModelStats.(cohortstr).dfjac=t4;
    hand=[hand; t1(:); t2(:); t3(:); t4(:)]; %#ok<AGROW>
end
fprintf('CohortsCustomStats hand-built: the entry ages read by the function, %i and %i (should be 1 and 5); the DFJ autocorrelation finite at %i entries of cohort 2 (should be positive) \n',TM.CustomModelStats.cohort1.entryage,TM.CustomModelStats.cohort2.entryage,sum(~isnan(TM.CustomModelStats.cohort2.dfjac)));
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
    fprintf('CohortsCustomStats (a): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('CohortsCustomStats (a) two cohorts, CustomModelStats targets (entry age, hasV, mean of inc, the DFJ autocorrelation) after AgeConditionalStats (cohort 1) and AutoCorr (cohort 2), %i target entries, at the truth: did not run, this should be zero: %i; objective, this should be zero: %.3e; current moments vs the hand-built ones, this should be zero: %.3e \n',nT,1-ran,objv,dmom)

%% (b) the cohorts with CustomModelStats targets run to N_j: identical to forcing it
estimoptions_b=estimoptions_a; estimoptions_b.cohortlastage=[N_j,N_j];
ran=0; dJ=NaN; dmom=NaN;
try
    [~,~,es_b]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_b,vfoptions,simoptions_a);
    ran=1;
    if ~isempty(es_a)
        dJ=CoreMaxAbs(es_a.variousmatrices.J(:)-es_b.variousmatrices.J(:)); dmom=CoreMaxAbs(es_a.currentmomentvec(:)-es_b.currentmomentvec(:));
    end
catch err
    fprintf('CohortsCustomStats (b): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('CohortsCustomStats (b) derived cutoffs vs both cohorts forced to N_j: did not run, this should be zero: %i; J, this should be zero: %.3e; current moments, this should be zero: %.3e \n',1-ran,dJ,dmom)

%% (c) must error
errcases={'CustomModelStats targets not under cohortC','no estimoptions.CustomModelStats','a targeted name the function does not return','a target of the wrong length'};
errmsgs={'must only contain cohort1','needs estimoptions.CustomModelStats','does not return a field of that name','entries, but its target'};
for ee=1:4
    TMbad=TM; estimoptions_c=estimoptions_a;
    nTb=nT+(ee>=3); % cases 3 and 4 add one targeted entry
    Wb=eye(nTb); Omegab=diag(0.01*(1+(1:nTb)'/nTb));
    switch ee
        case 1
            TMbad.CustomModelStats.incmean=TM.CustomModelStats.cohort1.incmean;
        case 2
            estimoptions_c=rmfield(estimoptions_c,'CustomModelStats');
        case 3
            TMbad.CustomModelStats.cohort2.nosuchstat=1;
        case 4
            TMbad.CustomModelStats.cohort1.hasV=[1;1];
    end
    errok=0;
    try
        EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMbad,Wb,Omegab,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_c,vfoptions,simoptions_a);
    catch ME
        errok=contains(ME.message,errmsgs{ee});
        if ~errok
            fprintf('   (the error was: %s) \n',ME.message);
        end
    end
    fprintf('CohortsCustomStats (c) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok);
end

output.hand=hand;
end
