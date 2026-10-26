function output=CalibInfHorzPType_Runtime(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn)
% CalibInfHorzPTypeTests part 6: the runtime of one objective evaluation with caliboptions.whichcombos=1 against =0 (informational).
%
% The objective is called directly (as the driver calls it), timed over several evaluations, for two target sets:
%   sparse: one restricted per-type AllStats mean, one per-type autocorrelation and one grouped cross-section covariance; with four
%           FnsToEvaluate and two restrictions in simoptions
%   dense:  Mean, Median, Variance and Gini of every function, grouped and per type, under every restriction (AllStats), the
%           autocorrelation of every function grouped and per type under every restriction, and both grouped cross-section matrices
%           under every restriction
% Prints the time per evaluation and the ratio =1/=0 for each. No identity is checked: the moments were shown identical in part 2.
% The model solve (two value functions and distributions) is inside the timing, as it is in a calibration.

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w_pt) (1+r)*a+w_pt*z-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);
FnNames=fieldnames(FnsToEvaluate);
nFns=length(FnNames);
nreps=5;

vfoptions=struct();
simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);
RestrNames=fieldnames(simoptions.conditionalrestrictions);

[z_gridvals,pi_zg,vfoptions]=ExogShockSetup_InfHorz_PType(n_z,z_grid,pi_z,Names_i,Params,vfoptions,3,0);
simoptions.e_gridvals=vfoptions.e_gridvals;
simoptions.pi_e=vfoptions.pi_e;
vfoptions.alreadygridvals=1;
simoptions.alreadygridvals=1;

% the two target sets (values do not matter for timing; sizes do)
TMsparse=struct();
TMsparse.AllStats.rich.earnings.low.Mean=1;
TMsparse.AutoCorrTransProbs.assets.high.AutoCorrelation=0.9;
TMsparse.CrossSectionCovarCorr.assets.CovarianceWith.earnings=0.1;
TMdense=struct();
for ff=1:nFns
    TMdense.AllStats.(FnNames{ff}).Mean=1;
    TMdense.AllStats.(FnNames{ff}).Median=1;
    TMdense.AllStats.(FnNames{ff}).Variance=1;
    TMdense.AllStats.(FnNames{ff}).Gini=1;
    TMdense.AutoCorrTransProbs.(FnNames{ff}).AutoCorrelation=0.9;
    for ii=1:N_i
        TMdense.AllStats.(FnNames{ff}).(Names_i{ii}).Mean=1;
        TMdense.AllStats.(FnNames{ff}).(Names_i{ii}).Median=1;
        TMdense.AllStats.(FnNames{ff}).(Names_i{ii}).Variance=1;
        TMdense.AllStats.(FnNames{ff}).(Names_i{ii}).Gini=1;
        TMdense.AutoCorrTransProbs.(FnNames{ff}).(Names_i{ii}).AutoCorrelation=0.9;
    end
    for rr=1:length(RestrNames)
        TMdense.AllStats.(RestrNames{rr}).(FnNames{ff}).Mean=1;
        TMdense.AllStats.(RestrNames{rr}).(FnNames{ff}).Median=1;
        TMdense.AllStats.(RestrNames{rr}).(FnNames{ff}).Variance=1;
        TMdense.AllStats.(RestrNames{rr}).(FnNames{ff}).Gini=1;
        TMdense.AutoCorrTransProbs.(RestrNames{rr}).(FnNames{ff}).AutoCorrelation=0.9;
        for ii=1:N_i
            TMdense.AllStats.(RestrNames{rr}).(FnNames{ff}).(Names_i{ii}).Mean=1;
            TMdense.AllStats.(RestrNames{rr}).(FnNames{ff}).(Names_i{ii}).Gini=1;
            TMdense.AutoCorrTransProbs.(RestrNames{rr}).(FnNames{ff}).(Names_i{ii}).AutoCorrelation=0.9;
        end
    end
end
Mfull=ones(nFns,nFns); Mfull(logical(eye(nFns)))=NaN;
TMdense.CrossSectionCovarCorr.CovarianceMatrix=Mfull;
TMdense.CrossSectionCovarCorr.CorrelationMatrix=Mfull;
for rr=1:length(RestrNames)
    TMdense.CrossSectionCovarCorr.(RestrNames{rr}).CovarianceMatrix=Mfull;
    TMdense.CrossSectionCovarCorr.(RestrNames{rr}).CorrelationMatrix=Mfull;
end
CalibParamNames={'beta'};
calibparamsvec0=Params.beta;
nCalibParams=1;
nCalibParamsFinder=[1,0];
calibparamsvecindex=[0;1];
calibparamssizes=[1,1];
calibomitparams_counter=0;
calibomitparamsmatrix={};
caliboptions=struct();
caliboptions.verbose=0;
caliboptions.constrainpositive={};
caliboptions.constrainpositivemethod='softplus';
caliboptions.constrain0to1={};
caliboptions.constrainAtoB={};
caliboptions.logmoments=0;
caliboptions.metric='sum_squared';
caliboptions.simulatemoments=0;
caliboptions.vectoroutput=0;
caliboptions.useCustomModelStats=0;
caliboptions.calibrateshocks=0;
[calibparamsvec0,caliboptions]=ParameterConstraints_PType_TransformParamsToUnconstrained(calibparamsvec0,calibparamsvecindex,CalibParamNames,nCalibParamsFinder,caliboptions,1);

output=struct();
for ss=1:2
    if ss==1
        TM=TMsparse; lab='sparse';
    else
        TM=TMdense; lab='dense';
    end
    [targetmomentvec,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,autocorrmomentnames,autocorrcummomentsizes,AutoCorrStats_whichstats,FnsToEvaluate_AutoCorr,crosssecmomentnames,crossseccummomentsizes,CrossSecStats_whichstats,FnsToEvaluate_CrossSec,cmsmomentnames,cmscummomentsizes,selectors,autocorrtimehorizons]=SetupTargetMoments_InfHorz(TM,FnsToEvaluate,1,simoptions,Names_i);
    caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
    caliboptions.selectors=selectors;
    caliboptions.autocorrtimehorizons=autocorrtimehorizons;
    argcell={CalibParamNames,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,pi_zg,ReturnFn,Params,DiscountFactorParamNames,PTypeDistParamNames,[],FnsToEvaluate,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,targetmomentvec,allstatmomentnames,autocorrmomentnames,crosssecmomentnames,cmsmomentnames,allstatcummomentsizes,autocorrcummomentsizes,crossseccummomentsizes,cmscummomentsizes,AllStats_whichstats,AutoCorrStats_whichstats,CrossSecStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_AutoCorr,FnsToEvaluate_CrossSec,nCalibParams,nCalibParamsFinder,calibparamsvecindex,calibparamssizes,calibomitparams_counter,calibomitparamsmatrix};
    times=zeros(2,1);
    for wc=[1,0]
        caliboptions_t=caliboptions; caliboptions_t.whichcombos=wc;
        CalibrateInfHorzAgentModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions); % warm-up (first-call costs)
        tic;
        for rep=1:nreps
            CalibrateInfHorzAgentModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions);
        end
        times(2-wc)=toc/nreps;
    end
    fprintf('Runtime PType, %s targets (%i targeted entries): one objective evaluation takes %.3f s with whichcombos=1 and %.3f s with whichcombos=0, ratio %.3f (model solve included) \n',lab,sum(~isnan(targetmomentvec)),times(1),times(2),times(1)/times(2))
    output.(lab)=times;
end

end
