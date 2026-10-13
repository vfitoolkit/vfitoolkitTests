function output=EstimFHorzPType_CohortsWindows(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% EstimFHorzPTypeTests part 8: cohort AutoCorr targets (TargetMoments.AutoCorr.cohortC...) of the window and log rows (added 2026-10-14)
% with permanent types, grouped and per type, in EstimateLifeCycleModel_PType_MethodOfMoments with estimoptions.cohortagejshifter (the
% route of the De Nardi, French and Jones (2010) medical-expense autocorrelations: two-period averages, their log, by cohort and type).
%   (a) two cohorts (entry ages 1 and 5): at the truth (skipestimation=1) the objective is zero and the current moments equal targets
%       built from each cohort's PType agent distribution with the PType command at the level and at the log (exact)
%   (b) the cutoff of cohort 2's agent distribution (its last window runs M-1 ages after its start): identical J and moments to forcing
%       N_j with estimoptions.cohortlastage (exact)
%   (c) must error: the level for one type where the log is targeted (grouped) at the same start age (simoptions.logwhichcombos has
%       no type dimension)
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate.inc=@(aprime,a,z,w,kappa_j_pt) 1+a+w*kappa_j_pt*z; % positive at every state
FnsInc.inc=FnsToEvaluate.inc;
Jr=Params.Jr;
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
estimoptions=struct(); estimoptions.verbose=0;
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
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
mewjc=repmat(Params.mewj,2,1); % the default rule: the population age weights from the entry age on
ptweights=Params.(PTypeDistParamNames{1}); ptweights=ptweights(:);
TM=struct(); hand=[];
for cc=1:2
    Params_c=Params; Params_c.mewj=mewjc(cc,:); Params_c.(PTypeDistParamNames{1})=ptweights.*typemass(:,cc);
    simoptions_c=simoptions; simoptions_c.jequaloneDistAge=entryages(cc); simoptions_c.agemass_withCohort=mewjc;
    jeq_c=struct();
    for ii=1:N_i
        jeq_c.(Names_i{ii})=reshape(jequaloneDist_c.(Names_i{ii})(:,:,entryages(cc))/typemass(ii,cc),[n_a,n_z]);
    end
    Dist_c=StationaryDist_Case1_FHorz_PType(jeq_c,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_c,simoptions_c);
    simoptions_L=simoptions_c; simoptions_L.ntimeperiods=[1,2]; simoptions_L.timehorizons=2;
    ACL=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(Dist_c,Policy,FnsToEvaluate,Params_c,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_L);
    simoptions_G=simoptions_L; simoptions_G.logwhichcombos=1;
    ACG=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(Dist_c,Policy,FnsInc,Params_c,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_G);
    if cc==1
        t=gather(ACL.earnings_2periodAverage.AutoCorrelation_horizon2); t(Jr-3:end)=NaN; % grouped, consecutive 2-period windows while working
        tl=gather(ACL.earnings_2periodAverage.low.AutoCorrelation_horizon2); tl(Jr-3:end)=NaN; % and type low's own
        TM.AutoCorr.cohort1.earnings_2periodAverage.AutoCorrelation_horizon2=t;
        TM.AutoCorr.cohort1.earnings_2periodAverage.low.AutoCorrelation_horizon2=tl;
        hand=[hand; t(:); tl(:)]; %#ok<AGROW>
    else
        t1=gather(ACG.inc_2periodAverage_log.Mean); t1(1:entryages(cc)-1)=NaN; t1(N_j-4:end)=NaN; % grouped; the last targeted window starts at N_j-5
        t2=gather(ACG.inc_log.high.AutoCorrelation_horizon2); t2(1:entryages(cc)-1)=NaN; t2(13:end)=NaN; % type high, the log of inc
        TM.AutoCorr.cohort2.inc_2periodAverage_log.Mean=t1;
        TM.AutoCorr.cohort2.inc_log.high.AutoCorrelation_horizon2=t2;
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
    [~,~,es_a]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_a,vfoptions,simoptions_a);
    ran=1; objv=es_a.objectivefnval;
    if length(es_a.currentmomentvec)==nT
        dmom=CoreMaxAbs(es_a.currentmomentvec(:)-hand);
    else
        dmom=abs(length(es_a.currentmomentvec)-nT);
    end
catch err
    fprintf('CohortsWindows PType (a): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('CohortsWindows PType (a) two cohorts, grouped and per-type targets of 2-period windows (level and log) and of a log, %i target entries, at the truth: did not run, this should be zero: %i; objective, this should be zero: %.3e; current moments vs the hand-built ones, this should be zero: %.3e \n',nT,1-ran,objv,dmom)

%% (b) cohort 2's cutoff includes the M-1 ages of its last window: identical to forcing N_j
estimoptions_b=estimoptions_a; estimoptions_b.cohortlastage=[N_j,N_j];
ran=0; dJ=NaN; dmom=NaN;
try
    [~,~,es_b]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_b,vfoptions,simoptions_a);
    ran=1;
    if ~isempty(es_a)
        dJ=CoreMaxAbs(es_a.variousmatrices.J(:)-es_b.variousmatrices.J(:)); dmom=CoreMaxAbs(es_a.currentmomentvec(:)-es_b.currentmomentvec(:));
    end
catch err
    fprintf('CohortsWindows PType (b): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('CohortsWindows PType (b) derived cutoffs vs both cohorts forced to N_j: did not run, this should be zero: %i; J, this should be zero: %.3e; current moments, this should be zero: %.3e \n',1-ran,dJ,dmom)

%% (c) must error: the level for type low where the grouped log is targeted, at the same start ages
TMbad=TM;
TMbad.AutoCorr.cohort2.inc_2periodAverage.low.Mean=TM.AutoCorr.cohort2.inc_2periodAverage_log.Mean;
errok=0;
try
    EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TMbad,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist_c,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_a,vfoptions,simoptions_a);
catch ME
    errok=contains(ME.message,'both a statistic of');
    if ~errok
        fprintf('   (the error was: %s) \n',ME.message);
    end
end
fprintf('CohortsWindows PType (c) the level for one type where the grouped log is targeted errors with the intended message, this should be zero: %i \n',1-errok);

output.hand=hand;
end
