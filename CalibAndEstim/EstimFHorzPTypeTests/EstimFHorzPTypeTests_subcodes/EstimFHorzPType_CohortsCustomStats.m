function output=EstimFHorzPType_CohortsCustomStats(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% EstimFHorzPTypeTests part 10: TargetMoments.CustomModelStats.cohortC.(name) in EstimateLifeCycleModel_PType_MethodOfMoments with
% estimoptions.cohortagejshifter (added 2026-10-14): the same estimoptions.CustomModelStats function as without cohorts
% (EstimFHorzPType_CohortsCustomStats_Fn, the PType signature with Names_i), called once per cohort on that cohort's PType agent
% distribution, with simoptions.jequaloneDistAge its entry age; stacked last within the cohort; the cohort's distribution runs to N_j.
%   (a) two cohorts (entry ages 1 and 5), CustomModelStats targets in both (with the grouped DFJ autocorrelation), after grouped
%       AgeConditionalStats (cohort 1) and AutoCorr (cohort 2) targets: at the truth the objective is zero and the current moments
%       equal the hand-built ones (exact)
%   (b) identical J and moments to forcing both cohorts to N_j with estimoptions.cohortlastage
%   (c) must error: a targeted name the function does not return
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate.inc=@(aprime,a,z,w,kappa_j_pt) 1+a+w*kappa_j_pt*z; % positive at every state
Jr=Params.Jr;
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
estimoptions=struct(); estimoptions.verbose=0;
estimoptions.eedefault=3; % the runs compared below share one step (the default 'auto' chooses per run)
estimoptions.CustomModelStats=@EstimFHorzPType_CohortsCustomStats_Fn;
entryages=[1,5];
typemass=[0.4,0.3; 0.2,0.1];
typemass=typemass./sum(typemass,2); % (a type's jequaloneDist has mass one across its slices)

%% The cohort initial distributions and the hand-built targets at the truth
jequaloneDist_c=struct();
for ii=1:N_i
    jequaloneDist_c.(Names_i{ii})=zeros(n_a,n_z,N_j,'gpuArray');
    for cc=1:2
        jequaloneDist_c.(Names_i{ii})(1,ceil(n_z/2),entryages(cc))=typemass(ii,cc);
    end
end
[V,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
mewjc=repmat(Params.mewj,2,1); % the default rule: the population age weights from the entry age on
ptweights=Params.(PTypeDistParamNames{1}); ptweights=ptweights(:);
TM=struct(); hand=[];
for cc=1:2
    Params_c=Params; Params_c.mewj=mewjc(cc,:); Params_c.(PTypeDistParamNames{1})=ptweights.*typemass(:,cc);
    simoptions_c=simoptions; simoptions_c.jequaloneDistAge=entryages(cc); simoptions_c.agemass_withCohort=mewjc; simoptions_c.warnzerorestrictedmass=0;
    jeq_c=struct();
    for ii=1:N_i
        jeq_c.(Names_i{ii})=reshape(jequaloneDist_c.(Names_i{ii})(:,:,entryages(cc))/typemass(ii,cc),[n_a,n_z]);
    end
    Dist_c=StationaryDist_Case1_FHorz_PType(jeq_c,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_c,simoptions_c);
    CS=EstimFHorzPType_CohortsCustomStats_Fn(V,Policy,Dist_c,Params_c,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,estimoptions,vfoptions,simoptions_c);
    cohortstr=['cohort',num2str(cc)];
    if cc==1
        simoptions_L=simoptions_c; simoptions_L.whichstats=[1,0,0,0,0,0,0];
        ACS=LifeCycleProfiles_FHorz_Case1_PType(Dist_c,Policy,FnsToEvaluate,Params_c,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_L);
        t=gather(ACS.earnings.Mean); t(Jr:end)=NaN;
        TM.AgeConditionalStats.cohort1.earnings.Mean=t;
        hand=[hand; t(:)]; %#ok<AGROW>
    else
        simoptions_A=simoptions_c; simoptions_A.timehorizons=2;
        ACP=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(Dist_c,Policy,FnsToEvaluate,Params_c,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_A);
        t=gather(ACP.assets.AutoCorrelation_horizon2); t(1:entryages(cc)-1)=NaN; t(12:end)=NaN;
        TM.AutoCorr.cohort2.assets.AutoCorrelation_horizon2=t;
        hand=[hand; t(:)]; %#ok<AGROW>
    end
    t1=CS.entryage; TM.CustomModelStats.(cohortstr).entryage=t1;
    t2=CS.hasV; TM.CustomModelStats.(cohortstr).hasV=t2;
    t3=CS.incmean; t3(1:entryages(cc)-1)=NaN; TM.CustomModelStats.(cohortstr).incmean=t3;
    t4=CS.dfjac; t4(1:entryages(cc)-1)=NaN; t4(14:end)=NaN; TM.CustomModelStats.(cohortstr).dfjac=t4;
    hand=[hand; t1(:); t2(:); t3(:); t4(:)]; %#ok<AGROW>
end
fprintf('CohortsCustomStats PType hand-built: the entry ages read by the function, %i and %i (should be 1 and 5); the grouped DFJ autocorrelation finite at %i entries of cohort 2 (should be positive) \n',TM.CustomModelStats.cohort1.entryage,TM.CustomModelStats.cohort2.entryage,sum(~isnan(TM.CustomModelStats.cohort2.dfjac)));
hand=hand(~isnan(hand));
nT=length(hand); W=eye(nT); Omega=diag(0.01*(1+(1:nT)'/nT));
EstimParamNames={'beta','pension'};

%% (a) at the truth
estimoptions_a=estimoptions; estimoptions_a.skipestimation=1; estimoptions_a.cohortagejshifter=entryages;
simoptions_a=simoptions; simoptions_a.agemass_withCohort=mewjc;
ran=0; objv=NaN; dmom=NaN; es_a=[];
try
    [~,~,es_a]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_a,vfoptions,simoptions_a);
    ran=1; objv=es_a.objectivefnval;
    if length(es_a.currentmomentvec)==nT
        dmom=CoreMaxAbs(es_a.currentmomentvec(:)-hand);
    else
        dmom=abs(length(es_a.currentmomentvec)-nT);
    end
catch err
    fprintf('CohortsCustomStats PType (a): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('CohortsCustomStats PType (a) two cohorts, CustomModelStats targets after grouped AgeConditionalStats (cohort 1) and AutoCorr (cohort 2), %i target entries, at the truth: did not run, this should be zero: %i; objective, this should be zero: %.3e; current moments vs the hand-built ones, this should be zero: %.3e \n',nT,1-ran,objv,dmom)

%% (b) identical to forcing N_j
estimoptions_b=estimoptions_a; estimoptions_b.cohortlastage=[N_j,N_j];
ran=0; dJ=NaN; dmom=NaN;
try
    [~,~,es_b]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_b,vfoptions,simoptions_a);
    ran=1;
    if ~isempty(es_a)
        dJ=CoreMaxAbs(es_a.variousmatrices.J(:)-es_b.variousmatrices.J(:)); dmom=CoreMaxAbs(es_a.currentmomentvec(:)-es_b.currentmomentvec(:));
    end
catch err
    fprintf('CohortsCustomStats PType (b): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('CohortsCustomStats PType (b) derived cutoffs vs both cohorts forced to N_j: did not run, this should be zero: %i; J, this should be zero: %.3e; current moments, this should be zero: %.3e \n',1-ran,dJ,dmom)

%% (c) must error: a targeted name the function does not return
TMbad=TM;
TMbad.CustomModelStats.cohort2.nosuchstat=1;
Wb=eye(nT+1); Omegab=diag(0.01*(1+(1:nT+1)'/(nT+1)));
errok=0;
try
    EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TMbad,Wb,Omegab,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_a,vfoptions,simoptions_a);
catch ME
    errok=contains(ME.message,'does not return a field of that name');
    if ~errok
        fprintf('   (the error was: %s) \n',ME.message);
    end
end
fprintf('CohortsCustomStats PType (c) a targeted name the function does not return errors with the intended message, this should be zero: %i \n',1-errok);

output.hand=hand;
end
