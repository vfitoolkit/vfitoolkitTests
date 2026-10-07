function output=CalibFHorz_MomentVector(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% CalibFHorzTests part 2: the moment vector that CalibrateLifeCycleModel_objectivefn builds, against a hand-built one.
%
% One solve at the truth, then full LifeCycleProfiles and AllStats runs (whichstats all ones). The targets are taken from those
% runs, with NaN where an entry is omitted: every statistic kind, four functions, two conditional restrictions, matrix-valued
% age-conditional targets (LorenzCurve, QuantileCutoffs). The objective is called directly, as CalibrateLifeCycleModel calls it,
% with vectoroutput=1 (the vector of model moments at the targeted entries) under caliboptions.whichcombos=1 and =0.
% Checks (all exact):
%   (a) each of the three vectors equals the hand-built vector (walk TargetMoments, read the same field of the full run,
%       vectorise column-major, drop the NaN targets) and has no NaN
%   (b) the three vectors equal SetupTargetMoments_FHorz's own targetmomentvec at the targeted entries (targets are the model)
%   (c) the selectors from SetupTargetMoments_FHorz handed straight to LifeCycleProfiles and AllStats: omitted ages are NaN,
%       untargeted statistics of a targeted function are NaN or absent, a function targeted only under a restriction has no
%       unrestricted struct
% (d_grid is empty, n_d is zero; DiscountFactorParamNames etc. are the bank's standard inputs.)

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) (1+r)*a+w*kappa_j*z*(agej<Jr)+pension*(agej>=Jr)-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);
Jr=Params.Jr;

vfoptions=struct();
simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);

%% Grids as CalibrateLifeCycleModel sets them up, then one solve and the full stats
[z_gridvals_J,pi_z_J,vfoptions]=ExogShockSetup_FHorz(n_z,z_grid,pi_z,N_j,Params,vfoptions,3,0);
simoptions.e_gridvals_J=vfoptions.e_gridvals_J;
simoptions.pi_e_J=vfoptions.pi_e_J;
vfoptions.alreadygridvals=1;
simoptions.alreadygridvals=1;
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_full);

%% Targets from the full runs, NaN where omitted
TM=struct();
TM.AllStats.assets.Mean=gather(AllFull.assets.Mean);
TM.AllStats.assets.Variance=gather(AllFull.assets.Variance);
TM.AllStats.assets.RatioMeanToMedian=gather(AllFull.assets.RatioMeanToMedian);
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);
TM.AllStats.assets.Minimum=gather(AllFull.assets.Minimum);
TM.AllStats.assets.Maximum=gather(AllFull.assets.Maximum);
TM.AllStats.assets.MoreInequality.Top1share=gather(AllFull.assets.MoreInequality.Top1share);
TM.AllStats.earnings.LorenzCurve=gather(AllFull.earnings.LorenzCurve);
TM.AllStats.earnings.QuantileCutoffs=gather(AllFull.earnings.QuantileCutoffs);
TM.AllStats.earnings.QuantileMeans=gather(AllFull.earnings.QuantileMeans);
TM.AllStats.richind.Mean=gather(AllFull.richind.Mean);
TM.AllStats.rich.earnings.Mean=gather(AllFull.rich.earnings.Mean);
TM.AllStats.rich.assets.StdDeviation=gather(AllFull.rich.assets.StdDeviation);
TM.AllStats.poor.earnings.Median=gather(AllFull.poor.earnings.Median);
TM.AllStats.poor.consumption.Mean=gather(AllFull.poor.consumption.Mean);
t=gather(LCPFull.assets.Mean); t(1)=NaN; TM.AgeConditionalStats.assets.Mean=t; % (everyone has zero assets at age 1)
t=gather(LCPFull.assets.Median); t(1)=NaN; TM.AgeConditionalStats.assets.Median=t;
t=gather(LCPFull.assets.Gini); t(1:2)=NaN; TM.AgeConditionalStats.assets.Gini=t;
t=gather(LCPFull.assets.QuantileCutoffs); t(:,1)=NaN; TM.AgeConditionalStats.assets.QuantileCutoffs=t; % matrix: age groups as columns
t=gather(LCPFull.earnings.Mean); t(Jr:end)=NaN; TM.AgeConditionalStats.earnings.Mean=t; % (earnings are zero in retirement)
t=gather(LCPFull.earnings.StdDeviation); t(Jr:end)=NaN; TM.AgeConditionalStats.earnings.StdDeviation=t;
t=gather(LCPFull.earnings.LorenzCurve); t(:,Jr:end)=NaN; TM.AgeConditionalStats.earnings.LorenzCurve=t; % matrix
TM.AgeConditionalStats.consumption.Mean=gather(LCPFull.consumption.Mean);
TM.AgeConditionalStats.consumption.RatioMeanToMedian=gather(LCPFull.consumption.RatioMeanToMedian);
t=gather(LCPFull.rich.earnings.Mean); t(1:3)=NaN; TM.AgeConditionalStats.rich.earnings.Mean=t;
TM.AgeConditionalStats.poor.assets.Mean=gather(LCPFull.poor.assets.Mean);
TM.AgeConditionalStats.rich.consumption.Minimum=gather(LCPFull.rich.consumption.Minimum);

% The hand-built vector: walk TargetMoments in field order and read the same field of the full run (the full-run structs have the
% same nesting as the targets, for restrictions and for MoreInequality alike); then the target vector by the same walk
hand=[]; tvec=[];
a1vec=fieldnames(TM.AllStats);
for a1=1:length(a1vec)
    a2vec=fieldnames(TM.AllStats.(a1vec{a1}));
    for a2=1:length(a2vec)
        tval=TM.AllStats.(a1vec{a1}).(a2vec{a2});
        if isstruct(tval)
            a3vec=fieldnames(tval);
            for a3=1:length(a3vec)
                v=gather(AllFull.(a1vec{a1}).(a2vec{a2}).(a3vec{a3})); hand=[hand; v(:)];
                v=tval.(a3vec{a3}); tvec=[tvec; v(:)];
            end
        else
            v=gather(AllFull.(a1vec{a1}).(a2vec{a2})); hand=[hand; v(:)];
            tvec=[tvec; tval(:)];
        end
    end
end
a1vec=fieldnames(TM.AgeConditionalStats);
for a1=1:length(a1vec)
    a2vec=fieldnames(TM.AgeConditionalStats.(a1vec{a1}));
    for a2=1:length(a2vec)
        tval=TM.AgeConditionalStats.(a1vec{a1}).(a2vec{a2});
        if isstruct(tval)
            a3vec=fieldnames(tval);
            for a3=1:length(a3vec)
                v=gather(LCPFull.(a1vec{a1}).(a2vec{a2}).(a3vec{a3})); hand=[hand; v(:)];
                v=tval.(a3vec{a3}); tvec=[tvec; v(:)];
            end
        else
            v=gather(LCPFull.(a1vec{a1}).(a2vec{a2})); hand=[hand; v(:)];
            tvec=[tvec; tval(:)];
        end
    end
end
actual=~isnan(tvec);
fprintf('MomentVector: %i target entries, %i of them targeted (the rest NaN); the model is NaN at a targeted entry, this should be zero: %.3e \n',length(tvec),sum(actual),sum(isnan(hand(actual))))

%% The objective, called as CalibrateLifeCycleModel calls it
CalibParamNames={'beta'};
calibparamsvec0=Params.beta;
calibparamsvecindex=[0;1];
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
caliboptions.vectoroutput=1;
caliboptions.useCustomModelStats=0;
caliboptions.calibrateshocks=0;
[calibparamsvec0,caliboptions]=ParameterConstraints_TransformParamsToUnconstrained(calibparamsvec0,calibparamsvecindex,CalibParamNames,caliboptions,1);
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,0,N_j,simoptions);
caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
argcell={CalibParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,[],Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,[],usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,calibparamsvecindex,calibomitparams_counter,calibomitparamsmatrix};

if length(targetmomentvec)==length(tvec)
    d=abs(targetmomentvec-tvec); d=d(:); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(tvec));
else
    d=1; patternbad=abs(length(targetmomentvec)-length(tvec));
end
fprintf('MomentVector: SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',max([d;0]),patternbad)

caliboptions_1=caliboptions; caliboptions_1.whichcombos=1; caliboptions_1.selectors=selectors;
m1=gather(CalibrateLifeCycleModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_1,vfoptions,simoptions));
caliboptions_0=caliboptions; caliboptions_0.whichcombos=0; caliboptions_0.selectors=selectors;
m0=gather(CalibrateLifeCycleModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_0,vfoptions,simoptions));
handA=hand(actual);
tA=targetmomentvec(actual);
for mm=1:2 % (caliboptions.whichcombos is always set by the driver, 1 or 0; since 2026-10-07 the objective takes no other value)
    if mm==1
        m=m1; lab='whichcombos=1';
    else
        m=m0; lab='whichcombos=0';
    end
    if length(m)==length(handA)
        d=max([abs(m(:)-handA(:));0]); dt=max([abs(m(:)-tA(:));0]);
    else
        d=1; dt=1;
    end
    fprintf('MomentVector (a) %s: moment vector vs hand-built, this should be zero: %.3e; NaN entries + length mismatch, this should be zero: %.3e \n',lab,d,sum(isnan(m))+abs(length(m)-length(handA)))
    fprintf('MomentVector (b) %s: moment vector vs the targets (targets are the model), this should be zero: %.3e \n',lab,dt)
end

%% (c) the selectors handed straight to the two stats commands
simsel=simoptions; simsel.whichcombos=selectors.ACStats.whichcombos; simsel.whichstats=selectors.ACStats.whichstats;
LCPsel=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_ACStats,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simsel);
bad=0;
bad=bad+sum(~isnan(gather(LCPsel.assets.Mean(1)))); % omitted age
bad=bad+sum(~isnan(gather(LCPsel.earnings.Mean(Jr:end)))); % omitted ages
bad=bad+sum(~isnan(gather(LCPsel.rich.earnings.Mean(1:3)))); % omitted ages under a restriction
bad=bad+sum(~isnan(gather(LCPsel.assets.QuantileCutoffs(:,1)))); % omitted column of a matrix target
if isfield(LCPsel.assets,'Variance')
    bad=bad+sum(~isnan(gather(LCPsel.assets.Variance))); % assets Variance is not an age-conditional target
end
if isfield(LCPsel.consumption,'Gini')
    bad=bad+sum(~isnan(gather(LCPsel.consumption.Gini))); % consumption Gini is not a target
end
bad=bad+isfield(LCPsel,'richind'); % richind has no age-conditional target
d=abs(gather(LCPsel.consumption.Mean)-gather(LCPFull.consumption.Mean)); d=d(:); d=d(~isnan(d));
fprintf('MomentVector (c) LifeCycleProfiles with the selectors: targeted consumption Mean vs full run, this should be zero: %.3e; untargeted entries not NaN or present, this should be zero: %.3e \n',max([d;0]),bad)
simselA=simoptions; simselA.whichcombos=selectors.AllStats.whichcombos; simselA.whichstats=selectors.AllStats.whichstats;
AllSel=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_AllStats,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simselA);
bad=0;
bad=bad+isfield(AllSel,'consumption'); % consumption is targeted only under poor
bad=bad+~isfield(AllSel.poor,'consumption');
if isfield(AllSel.poor,'consumption')
    bad=bad+isfield(AllSel.poor.consumption,'Median'); % only its Mean is targeted
end
bad=bad+isfield(AllSel.richind,'Gini'); % only richind's Mean is targeted
bad=bad+isfield(AllSel.rich,'consumption'); % nothing of consumption under rich
fprintf('MomentVector (c) AllStats with the selectors: targeted assets Mean vs full run, this should be zero: %.3e; untargeted structs present, this should be zero: %.3e \n',abs(gather(AllSel.assets.Mean)-gather(AllFull.assets.Mean)),bad)

output=struct();
output.TargetMoments=TM;
output.m1=m1;
output.hand=handA;

end
