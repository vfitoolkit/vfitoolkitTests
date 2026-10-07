function output=CalibFHorzPType_Runtime(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% CalibFHorzPTypeTests part 6: the runtime of one objective evaluation with caliboptions.whichcombos=1 against =0 (informational).
%
% The objective is called directly (as the driver calls it), timed over several evaluations, for two target sets:
%   sparse: the age-conditional mean of earnings of the rich, type low only, at ages 6-10, plus the grouped AllStats mean of assets;
%           with four FnsToEvaluate and two restrictions in simoptions
%   dense:  Mean, Median, Variance and Gini of every function at every age, grouped and per type, unrestricted and under both
%           restrictions, plus the AllStats Mean and Gini of every function, grouped and per type, under every restriction
% Prints the time per evaluation and the ratio =1/=0 for each. No identity is checked: the moments were shown identical in part 2.
% The model solve (two value functions and distributions) is inside the timing, as it is in a calibration.

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w,kappa_j_pt,agej,Jr,pension) (1+r)*a+w*kappa_j_pt*z*(agej<Jr)+pension*(agej>=Jr)-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);
FnNames=fieldnames(FnsToEvaluate);
nreps=5;

vfoptions=struct();
simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);
RestrNames=fieldnames(simoptions.conditionalrestrictions);

[z_gridvals_J,pi_z_J,vfoptions]=ExogShockSetup_FHorz_PType(n_z,z_grid,pi_z,N_j,Names_i,Params,vfoptions,3,0);
simoptions.e_gridvals_J=vfoptions.e_gridvals_J;
simoptions.pi_e_J=vfoptions.pi_e_J;
vfoptions.alreadygridvals=1;
simoptions.alreadygridvals=1;

% the two target sets (values do not matter for timing; sizes do)
TMsparse=struct();
t=nan(1,N_j); t(6:10)=1; TMsparse.AgeConditionalStats.rich.earnings.low.Mean=t;
TMsparse.AllStats.assets.Mean=1;
TMdense=struct();
slots=[Names_i,{''}]; % per type, then the grouped stats (no type level)
for ff=1:length(FnNames)
    for ss=1:length(slots)
        if isempty(slots{ss})
            TMdense.AgeConditionalStats.(FnNames{ff}).Mean=ones(1,N_j);
            TMdense.AgeConditionalStats.(FnNames{ff}).Median=ones(1,N_j);
            TMdense.AgeConditionalStats.(FnNames{ff}).Variance=ones(1,N_j);
            TMdense.AgeConditionalStats.(FnNames{ff}).Gini=ones(1,N_j);
            TMdense.AllStats.(FnNames{ff}).Mean=1;
            TMdense.AllStats.(FnNames{ff}).Gini=1;
            for rr=1:length(RestrNames)
                TMdense.AgeConditionalStats.(RestrNames{rr}).(FnNames{ff}).Mean=ones(1,N_j);
                TMdense.AgeConditionalStats.(RestrNames{rr}).(FnNames{ff}).Median=ones(1,N_j);
                TMdense.AgeConditionalStats.(RestrNames{rr}).(FnNames{ff}).Variance=ones(1,N_j);
                TMdense.AgeConditionalStats.(RestrNames{rr}).(FnNames{ff}).Gini=ones(1,N_j);
                TMdense.AllStats.(RestrNames{rr}).(FnNames{ff}).Mean=1;
                TMdense.AllStats.(RestrNames{rr}).(FnNames{ff}).Gini=1;
            end
        else
            TMdense.AgeConditionalStats.(FnNames{ff}).(slots{ss}).Mean=ones(1,N_j);
            TMdense.AgeConditionalStats.(FnNames{ff}).(slots{ss}).Median=ones(1,N_j);
            TMdense.AgeConditionalStats.(FnNames{ff}).(slots{ss}).Variance=ones(1,N_j);
            TMdense.AgeConditionalStats.(FnNames{ff}).(slots{ss}).Gini=ones(1,N_j);
            TMdense.AllStats.(FnNames{ff}).(slots{ss}).Mean=1;
            TMdense.AllStats.(FnNames{ff}).(slots{ss}).Gini=1;
            for rr=1:length(RestrNames)
                TMdense.AgeConditionalStats.(RestrNames{rr}).(FnNames{ff}).(slots{ss}).Mean=ones(1,N_j);
                TMdense.AgeConditionalStats.(RestrNames{rr}).(FnNames{ff}).(slots{ss}).Median=ones(1,N_j);
                TMdense.AgeConditionalStats.(RestrNames{rr}).(FnNames{ff}).(slots{ss}).Variance=ones(1,N_j);
                TMdense.AgeConditionalStats.(RestrNames{rr}).(FnNames{ff}).(slots{ss}).Gini=ones(1,N_j);
                TMdense.AllStats.(RestrNames{rr}).(FnNames{ff}).(slots{ss}).Mean=1;
                TMdense.AllStats.(RestrNames{rr}).(FnNames{ff}).(slots{ss}).Gini=1;
            end
        end
    end
end

CalibParamNames={'beta'};
calibparamsvec0=Params.beta;
nCalibParams=1;
nCalibParamsFinder=[1,0];
calibparamsvecindex=[0;1];
calibparamssizes=[1,1];
calibomitparams_counter=0;
calibomitparamsmatrix=zeros(N_j,1);
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
caliboptions.calibsemiexo=0;
[calibparamsvec0,caliboptions]=ParameterConstraints_PType_TransformParamsToUnconstrained(calibparamsvec0,calibparamsvecindex,CalibParamNames,nCalibParamsFinder,caliboptions,1);

output=struct();
for ss=1:2
    if ss==1
        TM=TMsparse; lab='sparse';
    else
        TM=TMdense; lab='dense';
    end
    [targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,1,N_j,simoptions,Names_i);
    caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
    caliboptions.selectors=selectors;
    argcell={CalibParamNames,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,nCalibParams,nCalibParamsFinder,calibparamsvecindex,calibparamssizes,calibomitparams_counter,calibomitparamsmatrix};
    times=zeros(2,1);
    for wc=[1,0]
        caliboptions_t=caliboptions; caliboptions_t.whichcombos=wc;
        CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions); % warm-up (first-call costs)
        tic;
        for rep=1:nreps
            CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions);
        end
        times(2-wc)=toc/nreps;
    end
    fprintf('Runtime PType, %s targets (%i targeted entries): one objective evaluation takes %.3f s with whichcombos=1 and %.3f s with whichcombos=0, ratio %.3f (model solve included) \n',lab,sum(~isnan(targetmomentvec)),times(1),times(2),times(1)/times(2))
    output.(lab)=times;
end

end
