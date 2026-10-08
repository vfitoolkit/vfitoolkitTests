function output=CalibFHorzPType_MomentVector(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% CalibFHorzPTypeTests part 2: the moment vector that CalibrateLifeCycleModel_PType_objectivefn builds, against a hand-built one.
%
% One solve at the truth, then full PType LifeCycleProfiles and AllStats runs (whichstats all ones). The targets are taken from those
% runs in every layout (grouped, per type, restricted, restricted per type, MoreInequality sub-statistics, matrix-valued age profiles),
% with NaN where an entry is omitted. The objective is called directly, as the driver calls it, with vectoroutput=1 under
% caliboptions.whichcombos=1 and =0. Checks (all exact):
%   (a) each of the three vectors equals the hand-built vector (walk TargetMoments, read the same field of the full run, vectorise
%       column-major, drop the NaN targets) and has no NaN
%   (b) the three vectors equal SetupTargetMoments_FHorz's own targetmomentvec at the targeted entries (targets are the model)
%   (c) the selectors handed straight to the two PType commands: omitted ages NaN, an untargeted type's statistic absent or NaN unless
%       a grouped target forces it (then present and equal), an untargeted function absent

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w,kappa_j_pt,agej,Jr,pension) (1+r)*a+w*kappa_j_pt*z*(agej<Jr)+pension*(agej>=Jr)-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);
Jr=Params.Jr;

vfoptions=struct();
simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);

%% Grids as the driver sets them up, then one solve and the full stats
[z_gridvals_J,pi_z_J,vfoptions]=ExogShockSetup_FHorz_PType(n_z,z_grid,pi_z,N_j,Names_i,Params,vfoptions,3,0);
simoptions.e_gridvals_J=vfoptions.e_gridvals_J;
simoptions.pi_e_J=vfoptions.pi_e_J;
vfoptions.alreadygridvals=1;
simoptions.alreadygridvals=1;
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params,simoptions);
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simoptions_full);

%% Targets from the full runs, NaN where omitted
TM=struct();
TM.AllStats.assets.Mean=gather(AllFull.assets.Mean);
TM.AllStats.assets.low.Mean=gather(AllFull.assets.low.Mean);
TM.AllStats.assets.high.Variance=gather(AllFull.assets.high.Variance);
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);
TM.AllStats.assets.MoreInequality.Top1share=gather(AllFull.assets.MoreInequality.Top1share);
TM.AllStats.assets.low.MoreInequality.Top10share=gather(AllFull.assets.low.MoreInequality.Top10share);
TM.AllStats.earnings.low.LorenzCurve=gather(AllFull.earnings.low.LorenzCurve);
TM.AllStats.earnings.QuantileCutoffs=gather(AllFull.earnings.QuantileCutoffs);
TM.AllStats.richind.Mean=gather(AllFull.richind.Mean);
TM.AllStats.rich.earnings.Mean=gather(AllFull.rich.earnings.Mean);
TM.AllStats.rich.earnings.high.Mean=gather(AllFull.rich.earnings.high.Mean);
TM.AllStats.poor.assets.low.StdDeviation=gather(AllFull.poor.assets.low.StdDeviation);
TM.AllStats.poor.consumption.Mean=gather(AllFull.poor.consumption.Mean);
t=gather(LCPFull.assets.Mean); t(1)=NaN; TM.AgeConditionalStats.assets.Mean=t;
t=gather(LCPFull.assets.low.Mean); t(1)=NaN; TM.AgeConditionalStats.assets.low.Mean=t;
t=gather(LCPFull.assets.high.Median); t(1)=NaN; TM.AgeConditionalStats.assets.high.Median=t;
t=gather(LCPFull.assets.QuantileCutoffs); t(:,1)=NaN; TM.AgeConditionalStats.assets.QuantileCutoffs=t; % matrix: ages as columns
t=gather(LCPFull.earnings.Mean); t(Jr:end)=NaN; TM.AgeConditionalStats.earnings.Mean=t;
t=gather(LCPFull.earnings.low.StdDeviation); t(Jr:end)=NaN; TM.AgeConditionalStats.earnings.low.StdDeviation=t;
t=gather(LCPFull.earnings.high.LorenzCurve); t(:,Jr:end)=NaN; TM.AgeConditionalStats.earnings.high.LorenzCurve=t; % matrix
TM.AgeConditionalStats.consumption.Mean=gather(LCPFull.consumption.Mean);
TM.AgeConditionalStats.consumption.low.RatioMeanToMedian=gather(LCPFull.consumption.low.RatioMeanToMedian);
t=gather(LCPFull.rich.earnings.Mean); t(1:3)=NaN; TM.AgeConditionalStats.rich.earnings.Mean=t;
t=gather(LCPFull.rich.earnings.low.Mean); t(1:3)=NaN; TM.AgeConditionalStats.rich.earnings.low.Mean=t;
TM.AgeConditionalStats.poor.assets.high.Mean=gather(LCPFull.poor.assets.high.Mean);
TM.AgeConditionalStats.poor.earnings.high.Mean=nan(1,N_j); % a per-type target that is NaN in full: it must warn and be dropped, with whichcombos=1 nothing is computed for it (2026-10-09)

% The hand-built vector: walk TargetMoments (two to four levels) and read the same field of the full run (same nesting); the target
% vector by the same walk
hand=[]; tvec=[];
for ss=1:2
    if ss==1
        TMs=TM.AllStats; Fs=AllFull;
    else
        TMs=TM.AgeConditionalStats; Fs=LCPFull;
    end
    a1vec=fieldnames(TMs);
    for a1=1:length(a1vec)
        a2vec=fieldnames(TMs.(a1vec{a1}));
        for a2=1:length(a2vec)
            t2=TMs.(a1vec{a1}).(a2vec{a2});
            if isstruct(t2)
                a3vec=fieldnames(t2);
                for a3=1:length(a3vec)
                    t3=t2.(a3vec{a3});
                    if isstruct(t3)
                        a4vec=fieldnames(t3);
                        for a4=1:length(a4vec)
                            v=gather(Fs.(a1vec{a1}).(a2vec{a2}).(a3vec{a3}).(a4vec{a4})); hand=[hand; v(:)];
                            v=t3.(a4vec{a4}); tvec=[tvec; v(:)];
                        end
                    else
                        v=gather(Fs.(a1vec{a1}).(a2vec{a2}).(a3vec{a3})); hand=[hand; v(:)];
                        tvec=[tvec; t3(:)];
                    end
                end
            else
                v=gather(Fs.(a1vec{a1}).(a2vec{a2})); hand=[hand; v(:)];
                tvec=[tvec; t2(:)];
            end
        end
    end
end
actual=~isnan(tvec);
fprintf('MomentVector PType: %i target entries, %i of them targeted (the rest NaN); the model is NaN at a targeted entry, this should be zero: %.3e \n',length(tvec),sum(actual),sum(isnan(hand(actual))))

%% The objective, called as the driver calls it (one common parameter)
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
caliboptions.vectoroutput=1;
caliboptions.useCustomModelStats=0;
caliboptions.calibrateshocks=0;
caliboptions.calibsemiexo=0;
[calibparamsvec0,caliboptions]=ParameterConstraints_PType_TransformParamsToUnconstrained(calibparamsvec0,calibparamsvecindex,CalibParamNames,nCalibParamsFinder,caliboptions,1);
lastwarn('');
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,1,N_j,simoptions,Names_i);
wmsg=lastwarn;
fprintf('MomentVector: the target that is NaN in full makes SetupTargetMoments_FHorz warn (the warning names it), did not, this should be zero: %.3e \n',~contains(wmsg,'poor.earnings.high.Mean is NaN in full'))
caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
argcell={CalibParamNames,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,nCalibParams,nCalibParamsFinder,calibparamsvecindex,calibparamssizes,calibomitparams_counter,calibomitparamsmatrix};

if length(targetmomentvec)==length(tvec)
    d=abs(targetmomentvec-tvec); d=d(:); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(tvec));
else
    d=1; patternbad=abs(length(targetmomentvec)-length(tvec));
end
fprintf('MomentVector PType: SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',max([d;0]),patternbad)

caliboptions_1=caliboptions; caliboptions_1.whichcombos=1; caliboptions_1.selectors=selectors;
m1=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_1,vfoptions,simoptions));
caliboptions_0=caliboptions; caliboptions_0.whichcombos=0; caliboptions_0.selectors=selectors;
m0=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_0,vfoptions,simoptions));
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
    fprintf('MomentVector PType (a) %s: moment vector vs hand-built, this should be zero: %.3e; NaN entries + length mismatch, this should be zero: %.3e \n',lab,d,sum(isnan(m))+abs(length(m)-length(handA)))
    fprintf('MomentVector PType (b) %s: moment vector vs the targets (targets are the model), this should be zero: %.3e \n',lab,dt)
end

%% (c) the selectors handed straight to the two PType commands
simsel=simoptions; simsel.whichcombos=selectors.ACStats.whichcombos; simsel.whichstats=selectors.ACStats.whichstats;
LCPsel=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate_ACStats,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simsel);
bad=0; maxd=0;
bad=bad+sum(~isnan(gather(LCPsel.assets.Mean(1)))); % omitted age of the grouped target
bad=bad+sum(~isnan(gather(LCPsel.assets.low.Mean(1)))); % omitted age of a per-type target
bad=bad+sum(~isnan(gather(LCPsel.earnings.Mean(Jr:end))));
bad=bad+sum(~isnan(gather(LCPsel.rich.earnings.low.Mean(1:3))));
% the grouped assets Mean target forces the high type's Mean (byproduct): present and equal to the full run at the targeted ages
if isfield(LCPsel.assets,'high') && isfield(LCPsel.assets.high,'Mean')
    a=gather(LCPsel.assets.high.Mean(2:end)); b=gather(LCPFull.assets.high.Mean(2:end)); maxd=max([maxd;abs(a(:)-b(:))]); bad=bad+sum(isnan(a(:)));
else
    bad=bad+1;
end
% the high type's assets Variance is not targeted anywhere, nor the grouped one: absent or NaN
if isfield(LCPsel.assets,'high') && isfield(LCPsel.assets.high,'Variance')
    bad=bad+sum(~isnan(gather(LCPsel.assets.high.Variance)));
end
bad=bad+isfield(LCPsel,'richind'); % richind has no age-conditional target
a=gather(LCPsel.consumption.Mean); b=gather(LCPFull.consumption.Mean); maxd=max([maxd;abs(a(:)-b(:))]);
fprintf('MomentVector PType (c) LifeCycleProfiles PType with the selectors: targeted and forced entries vs full run, this should be zero: %.3e; untargeted entries not NaN or present, forced entries missing, this should be zero: %.3e \n',maxd,bad)
simselA=simoptions; simselA.whichcombos=selectors.AllStats.whichcombos; simselA.whichstats=selectors.AllStats.whichstats;
AllSel=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate_AllStats,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simselA);
bad=0; maxd=0;
bad=bad+~isfield(AllSel.poor,'consumption'); % targeted under poor (grouped), so present
bad=bad+isfield(AllSel,'consumption'); % consumption is not targeted unrestricted
if isfield(AllSel.poor,'consumption')
    bad=bad+~isfield(AllSel.poor.consumption,'low'); % the grouped poor consumption Mean forces both types' (byproducts)
    bad=bad+~isfield(AllSel.poor.consumption,'high');
    if isfield(AllSel.poor.consumption,'high')
        maxd=max(maxd,abs(gather(AllSel.poor.consumption.high.Mean)-gather(AllFull.poor.consumption.high.Mean)));
    end
end
bad=bad+isfield(AllSel.richind,'Gini'); % only richind's grouped Mean is targeted
bad=bad+isfield(AllSel.rich,'assets'); % nothing of assets under rich
maxd=max(maxd,abs(gather(AllSel.assets.low.Mean)-gather(AllFull.assets.low.Mean)));
fprintf('MomentVector PType (c) AllStats PType with the selectors: targeted and forced entries vs full run, this should be zero: %.3e; untargeted structs present, forced structs missing, this should be zero: %.3e \n',maxd,bad)

output=struct();
output.TargetMoments=TM;
output.m1=m1;

end
