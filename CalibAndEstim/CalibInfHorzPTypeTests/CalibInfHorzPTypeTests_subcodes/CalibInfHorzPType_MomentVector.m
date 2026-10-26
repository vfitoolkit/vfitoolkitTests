function output=CalibInfHorzPType_MomentVector(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn)
% CalibInfHorzPTypeTests part 2: the moment vector that CalibrateInfHorzAgentModel_PType_objectivefn builds, against a hand-built one.
%
% One solve at the truth, then full runs of the three InfHorz PType stats commands (AllStats with whichstats all ones,
% AutoCorrTransProbs with horizons 2 and 3, CrossSectionCovarCorr), all under two conditional restrictions. The targets are taken
% from those runs, grouped and per type, with NaN where an entry is omitted, in every layout (two to five names), including a
% per-type AutoCorr target that is NaN in full (it must warn and be dropped). The objective is called directly, as the driver calls
% it, with vectoroutput=1 under caliboptions.whichcombos=1 and =0.
% Checks (all exact):
%   (a) each of the two vectors equals the hand-built vector (walk TargetMoments, read the same field of the full run, a _kK suffix
%       read from the command's .tperiodsK nesting, grouped or per type, vectorise column-major, drop the NaN targets) and has no NaN
%   (b) the two vectors equal SetupTargetMoments_InfHorz's own targetmomentvec at the targeted entries (targets are the model)
%   (c) the selectors from SetupTargetMoments_InfHorz handed straight to the three PType commands: a grouped target forces the
%       per-type statistics (present and equal to the full run), an untargeted function or type is absent, an unselected pair or
%       omitted matrix entry is NaN, the targeted entries equal the full runs
% (d_grid is empty, n_d is zero; the other inputs are the bank's standard ones.)

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w_pt) (1+r)*a+w_pt*z-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);

vfoptions=struct();
simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);

%% Grids as the driver sets them up, then one solve and the full stats
[z_gridvals,pi_zg,vfoptions]=ExogShockSetup_InfHorz_PType(n_z,z_grid,pi_z,Names_i,Params,vfoptions,3,0);
simoptions.e_gridvals=vfoptions.e_gridvals;
simoptions.pi_e=vfoptions.pi_e;
vfoptions.alreadygridvals=1;
simoptions.alreadygridvals=1;
[~,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,pi_zg,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_zg,Params,simoptions);
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,simoptions_full);
simoptions_ac=simoptions;
simoptions_ac.timehorizons=[2,3];
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,pi_zg,simoptions_ac);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,simoptions);
Full={AllFull,ACFull,CSFull};
kinds={'AllStats','AutoCorrTransProbs','CrossSectionCovarCorr'};

%% Targets from the full runs, grouped and per type, NaN where omitted
TM=struct();
TM.AllStats.assets.Mean=gather(AllFull.assets.Mean);
TM.AllStats.assets.low.Variance=gather(AllFull.assets.low.Variance);
TM.AllStats.assets.RatioMeanToMedian=gather(AllFull.assets.RatioMeanToMedian);
TM.AllStats.assets.high.Gini=gather(AllFull.assets.high.Gini);
TM.AllStats.assets.Minimum=gather(AllFull.assets.Minimum);
TM.AllStats.assets.MoreInequality.Top1share=gather(AllFull.assets.MoreInequality.Top1share);
TM.AllStats.earnings.LorenzCurve=gather(AllFull.earnings.LorenzCurve);
TM.AllStats.earnings.low.QuantileCutoffs=gather(AllFull.earnings.low.QuantileCutoffs);
TM.AllStats.earnings.QuantileMeans=gather(AllFull.earnings.QuantileMeans);
TM.AllStats.richind.Mean=gather(AllFull.richind.Mean);
TM.AllStats.rich.earnings.Mean=gather(AllFull.rich.earnings.Mean);
TM.AllStats.rich.assets.low.StdDeviation=gather(AllFull.rich.assets.low.StdDeviation);
TM.AllStats.poor.earnings.high.Median=gather(AllFull.poor.earnings.high.Median);
TM.AllStats.poor.consumption.Mean=gather(AllFull.poor.consumption.Mean);
TM.AllStats.rich.consumption.StdDeviation=gather(AllFull.rich.consumption.StdDeviation); % a grouped StdDeviation as the ONLY target of its (function, restriction): the grouped value is rebuilt from the per-type Means and StdDeviations, so the selector must force the per-type Means too (the FHorz PType commands had this defect, fixed 2026-10-07)
TM.AutoCorrTransProbs.assets.AutoCorrelation=gather(ACFull.assets.AutoCorrelation);
TM.AutoCorrTransProbs.assets.low.AutoCorrelation_k3=gather(ACFull.assets.low.tperiods3.AutoCorrelation);
TM.AutoCorrTransProbs.earnings.AutoCovariance=gather(ACFull.earnings.AutoCovariance);
TM.AutoCorrTransProbs.earnings.high.Mean=gather(ACFull.earnings.high.Mean);
TM.AutoCorrTransProbs.rich.consumption.AutoCorrelation_k2=gather(ACFull.rich.consumption.tperiods2.AutoCorrelation);
TM.AutoCorrTransProbs.poor.assets.low.StdDeviation=gather(ACFull.poor.assets.low.StdDeviation);
TM.AutoCorrTransProbs.rich.assets.PairMass=gather(ACFull.rich.assets.PairMass);
TM.AutoCorrTransProbs.poor.earnings.high.AutoCorrelation=NaN; % omitted in full (warns, dropped)
TM.CrossSectionCovarCorr.assets.CovarianceWith.earnings=gather(CSFull.assets.CovarianceWith.earnings);
TM.CrossSectionCovarCorr.earnings.high.CorrelationWith.consumption=gather(CSFull.earnings.high.CorrelationWith.consumption);
TM.CrossSectionCovarCorr.consumption.low.StdDeviation=gather(CSFull.consumption.low.StdDeviation);
TM.CrossSectionCovarCorr.consumption.Mean=gather(CSFull.consumption.Mean);
TM.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption=gather(CSFull.rich.assets.CorrelationWith.consumption);
TM.CrossSectionCovarCorr.rich.earnings.low.CovarianceWith.consumption=gather(CSFull.rich.earnings.low.CovarianceWith.consumption); % five names
t=gather(CSFull.poor.CovarianceMatrix); t(1,4)=NaN; t(4,1)=NaN; TM.CrossSectionCovarCorr.poor.CovarianceMatrix=t; % grouped matrix under poor, one pair omitted
t=gather(CSFull.CorrelationMatrix); t(logical(eye(4)))=NaN; TM.CrossSectionCovarCorr.CorrelationMatrix=t; % grouped matrix, the diagonal (always one) omitted

% The hand-built vector: a depth-first walk of TargetMoments in field order, reading the same field of the full run (the full-run
% structs have the same nesting as the targets, grouped or per type; an AutoCorr statistic with a _kK suffix is .tperiodsK.(statistic))
hand=[]; tvec=[];
for kk=1:3
    stack={{}};
    while ~isempty(stack)
        p=stack{end}; stack(end)=[];
        node=TM.(kinds{kk});
        for k=1:numel(p)
            node=node.(p{k});
        end
        if isstruct(node)
            fn=fieldnames(node);
            for k=numel(fn):-1:1
                stack{end+1}=[p,fn(k)]; %#ok<AGROW>
            end
        else
            tvec=[tvec; node(:)]; %#ok<AGROW>
            q=p;
            if kk==2
                tk=regexp(q{end},'_k(\d+)$','tokens','once');
                if ~isempty(tk)
                    q=[q(1:end-1),{['tperiods',tk{1}]},{q{end}(1:end-length(tk{1})-2)}];
                end
            end
            v=Full{kk};
            for k=1:numel(q)
                v=v.(q{k});
            end
            v=gather(v);
            hand=[hand; v(:)]; %#ok<AGROW>
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
caliboptions.vectoroutput=1;
caliboptions.useCustomModelStats=0;
caliboptions.calibrateshocks=0;
[calibparamsvec0,caliboptions]=ParameterConstraints_PType_TransformParamsToUnconstrained(calibparamsvec0,calibparamsvecindex,CalibParamNames,nCalibParamsFinder,caliboptions,1);
ran=0; wmsg='';
lastwarn('');
try
    [targetmomentvec,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,autocorrmomentnames,autocorrcummomentsizes,AutoCorrStats_whichstats,FnsToEvaluate_AutoCorr,crosssecmomentnames,crossseccummomentsizes,CrossSecStats_whichstats,FnsToEvaluate_CrossSec,cmsmomentnames,cmscummomentsizes,selectors,autocorrtimehorizons]=SetupTargetMoments_InfHorz(TM,FnsToEvaluate,1,simoptions,Names_i);
    ran=1; wmsg=lastwarn;
catch err
    fprintf('MomentVector PType: SetupTargetMoments_InfHorz errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('MomentVector PType: SetupTargetMoments_InfHorz with simoptions and Names_i, did not run, this should be zero: %.3e; the target that is NaN in full warns (the warning names it), did not, this should be zero: %.3e \n',1-ran,~contains(wmsg,'poor.earnings.high.AutoCorrelation is NaN in full'))
output=struct();
output.TargetMoments=TM;
output.hand=hand(actual);
if ran==0
    return
end
caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
caliboptions.selectors=selectors;
caliboptions.autocorrtimehorizons=autocorrtimehorizons;
argcell={CalibParamNames,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,pi_zg,ReturnFn,Params,DiscountFactorParamNames,PTypeDistParamNames,[],FnsToEvaluate,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,targetmomentvec,allstatmomentnames,autocorrmomentnames,crosssecmomentnames,cmsmomentnames,allstatcummomentsizes,autocorrcummomentsizes,crossseccummomentsizes,cmscummomentsizes,AllStats_whichstats,AutoCorrStats_whichstats,CrossSecStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_AutoCorr,FnsToEvaluate_CrossSec,nCalibParams,nCalibParamsFinder,calibparamsvecindex,calibparamssizes,calibomitparams_counter,calibomitparamsmatrix};

if length(targetmomentvec)==length(tvec)
    d=abs(targetmomentvec-tvec); d=d(:); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(tvec));
else
    d=1; patternbad=abs(length(targetmomentvec)-length(tvec));
end
fprintf('MomentVector PType: the horizons read off the _kK suffixes are [%s] (hand [2 3]) \n',num2str(autocorrtimehorizons(:)'))
fprintf('MomentVector PType: SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e; horizons vs hand, this should be zero: %.3e \n',max([d;0]),patternbad,~isequal(autocorrtimehorizons(:)',[2,3]))

handA=hand(actual);
tA=targetmomentvec(actual);
for mm=1:2
    caliboptions_m=caliboptions; caliboptions_m.whichcombos=2-mm;
    if mm==1
        lab='whichcombos=1';
    else
        lab='whichcombos=0';
    end
    ranm=0; m=[];
    try
        m=gather(CalibrateInfHorzAgentModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_m,vfoptions,simoptions));
        ranm=1;
    catch err
        fprintf('MomentVector PType (a) %s: the objective errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ranm==1 && length(m)==length(handA)
        d=CoreMaxAbs(m(:)-handA(:)); dt=CoreMaxAbs(m(:)-tA(:)); lenbad=0;
    else
        d=1; dt=1; lenbad=1;
    end
    fprintf('MomentVector PType (a) %s: did not run, this should be zero: %.3e; moment vector vs hand-built, this should be zero: %.3e; length mismatch, this should be zero: %.3e \n',lab,1-ranm,d,lenbad)
    fprintf('MomentVector PType (b) %s: moment vector vs the targets (targets are the model), this should be zero: %.3e \n',lab,dt)
    if mm==1
        output.m1=m;
    end
end

%% (c) the selectors handed straight to the three PType commands
% AllStats
simselA=simoptions; simselA.whichcombos=selectors.AllStats.whichcombos; simselA.whichstats=selectors.AllStats.whichstats;
AllSel=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate_AllStats,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,simselA);
bad=0; d=[];
d=[d; gather(AllSel.assets.Mean)-gather(AllFull.assets.Mean)]; % grouped target
if isfield(AllSel.assets,'low') && isfield(AllSel.assets.low,'Mean') % the grouped Mean forces both types' Means (byproducts)
    d=[d; gather(AllSel.assets.low.Mean)-gather(AllFull.assets.low.Mean)];
else
    bad=bad+1;
end
if isfield(AllSel.assets,'high') && isfield(AllSel.assets.high,'Gini') % a per-type target
    d=[d; gather(AllSel.assets.high.Gini)-gather(AllFull.assets.high.Gini)];
else
    bad=bad+1;
end
bad=bad+isfield(AllSel,'consumption'); % consumption is targeted only under poor (grouped)
bad=bad+~isfield(AllSel.poor,'consumption');
if isfield(AllSel.rich,'consumption') && isfield(AllSel.rich.consumption,'StdDeviation') % the grouped StdDeviation as the only target of consumption under rich
    d=[d; gather(AllSel.rich.consumption.StdDeviation)-gather(AllFull.rich.consumption.StdDeviation)];
else
    bad=bad+1;
end
bad=bad+isfield(AllSel.richind,'Gini'); % only richind's grouped Mean is targeted
if isfield(AllSel.poor,'earnings') && isfield(AllSel.poor.earnings,'high') % a per-type target under a restriction
    d=[d; gather(AllSel.poor.earnings.high.Median)-gather(AllFull.poor.earnings.high.Median)];
else
    bad=bad+1;
end
fprintf('MomentVector PType (c) AllStats PType with the selectors: targeted and forced entries vs full run, this should be zero: %.3e; untargeted structs present, forced structs missing, this should be zero: %.3e \n',CoreMaxAbs(d),bad)
% AutoCorrTransProbs
simselC=simoptions; simselC.whichcombos=selectors.AutoCorr.whichcombos; simselC.timehorizons=autocorrtimehorizons;
ACSel=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate_AutoCorr,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,pi_zg,simselC);
bad=0; d=[];
d=[d; gather(ACSel.assets.AutoCorrelation)-gather(ACFull.assets.AutoCorrelation)]; % grouped target
if isfield(ACSel.assets,'low') && isfield(ACSel.assets.low,'tperiods3') % per-type target at horizon 3
    d=[d; gather(ACSel.assets.low.tperiods3.AutoCorrelation)-gather(ACFull.assets.low.tperiods3.AutoCorrelation)];
else
    bad=bad+1;
end
if isfield(ACSel.rich,'consumption') && isfield(ACSel.rich.consumption,'tperiods2') % grouped target under rich at horizon 2
    d=[d; gather(ACSel.rich.consumption.tperiods2.AutoCorrelation)-gather(ACFull.rich.consumption.tperiods2.AutoCorrelation)];
else
    bad=bad+1;
end
bad=bad+isfield(ACSel,'richind'); % richind has no AutoCorr target (not in the reduced set)
bad=bad+isfield(ACSel,'consumption'); % consumption is targeted only under rich
bad=bad+isfield(ACSel.poor,'earnings'); % the only target of earnings under poor is NaN in full, so nothing is computed for it
fprintf('MomentVector PType (c) AutoCorrTransProbs PType with the selector and horizons: targeted entries vs full run, this should be zero: %.3e; untargeted structs present, targeted ones missing, this should be zero: %.3e \n',CoreMaxAbs(d),bad)
% CrossSectionCovarCorr
simselX=simoptions; simselX.whichcombos=selectors.CrossSec.whichcombos;
CSSel=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate_CrossSec,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,simselX);
bad=0; d=[];
d=[d; gather(CSSel.assets.CovarianceWith.earnings)-gather(CSFull.assets.CovarianceWith.earnings)]; % grouped pair
if isfield(CSSel.earnings,'high') && isfield(CSSel.earnings.high,'CorrelationWith') % per-type pair
    d=[d; gather(CSSel.earnings.high.CorrelationWith.consumption)-gather(CSFull.earnings.high.CorrelationWith.consumption)];
else
    bad=bad+1;
end
if isfield(CSSel.rich,'earnings') && isfield(CSSel.rich.earnings,'low') && isfield(CSSel.rich.earnings.low,'CovarianceWith') % five-name target
    d=[d; gather(CSSel.rich.earnings.low.CovarianceWith.consumption)-gather(CSFull.rich.earnings.low.CovarianceWith.consumption)];
else
    bad=bad+1;
end
if isfield(CSSel.rich,'earnings')
    bad=bad+isfield(CSSel.rich.earnings,'high'); % the high type evaluates nothing of earnings on the rich page
end
Mp=gather(CSSel.poor.CovarianceMatrix);
bad=bad+~isnan(Mp(1,4))+~isnan(Mp(4,1)); % the omitted off-diagonal entries stay NaN
d=[d; Mp(1,2)-gather(CSFull.poor.CovarianceMatrix(1,2))];
Mc=gather(CSSel.CorrelationMatrix); Mcf=gather(CSFull.CorrelationMatrix);
d=[d; Mc(~eye(4))-Mcf(~eye(4))]; % every off-diagonal pair is targeted unrestricted (grouped)
fprintf('MomentVector PType (c) CrossSectionCovarCorr PType with the pair selector: targeted pairs and matrix entries vs full run, this should be zero: %.3e; an unevaluated type present, a targeted one missing, or omitted matrix entries not NaN, this should be zero: %.3e \n',CoreMaxAbs(d),bad)

end
