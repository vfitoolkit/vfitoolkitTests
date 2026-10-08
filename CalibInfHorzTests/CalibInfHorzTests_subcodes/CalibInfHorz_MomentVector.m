function output=CalibInfHorz_MomentVector(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn)
% CalibInfHorzTests part 2: the moment vector that CalibrateInfHorzAgentModel_objectivefn builds, against a hand-built one.
%
% One solve at the truth, then full AllStats (whichstats all ones), AutoCorrTransProbs (horizons 2 and 3) and CrossSectionCovarCorr
% runs, all under two conditional restrictions. The targets are taken from those runs, with NaN where an entry is omitted: every
% statistic kind, four functions, restricted targets of all three kinds, a horizon suffix, two matrix targets. The objective is
% called directly, as CalibrateInfHorzAgentModel calls it, with vectoroutput=1 (the vector of model moments at the targeted entries)
% under caliboptions.whichcombos=1 and =0 (the selectors and the horizons travel in caliboptions, as caliboptions.selectors and
% caliboptions.autocorrtimehorizons).
% Checks (all exact):
%   (a) each of the two vectors equals the hand-built vector (walk TargetMoments, read the same field of the full run, a
%       _kK suffix read from the command's .tperiodsK nesting, vectorise column-major, drop the NaN targets) and has no NaN
%   (b) the two vectors equal SetupTargetMoments_InfHorz's own targetmomentvec at the targeted entries (targets are the model)
%   (c) the selectors from SetupTargetMoments_InfHorz handed straight to the three commands: untargeted statistics of a targeted
%       function are NaN or absent, a function targeted only under a restriction has no unrestricted struct, an unselected pair
%       is NaN, and the targeted entries equal the full runs
% (d_grid is empty, n_d is zero; DiscountFactorParamNames etc. are the bank's standard inputs.)

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w) (1+r)*a+w*z-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);

vfoptions=struct();
simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);

%% Grids as CalibrateInfHorzAgentModel sets them up, then one solve and the full stats
[z_gridvals,pi_zg,vfoptions]=ExogShockSetup_InfHorz(n_z,z_grid,pi_z,Params,vfoptions,3,0);
simoptions.e_gridvals=vfoptions.e_gridvals;
simoptions.pi_e=vfoptions.pi_e;
vfoptions.alreadygridvals=1;
simoptions.alreadygridvals=1;
[~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zg,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_zg,simoptions,Params);
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions_full);
simoptions_ac=simoptions;
simoptions_ac.timehorizons=[2,3];
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zg,simoptions_ac);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions);
Full={AllFull,ACFull,CSFull};
kinds={'AllStats','AutoCorrTransProbs','CrossSectionCovarCorr'};

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
TM.AutoCorrTransProbs.assets.AutoCorrelation=gather(ACFull.assets.AutoCorrelation);
TM.AutoCorrTransProbs.assets.AutoCorrelation_k3=gather(ACFull.assets.tperiods3.AutoCorrelation);
TM.AutoCorrTransProbs.earnings.AutoCovariance=gather(ACFull.earnings.AutoCovariance);
TM.AutoCorrTransProbs.earnings.Mean=gather(ACFull.earnings.Mean);
TM.AutoCorrTransProbs.rich.consumption.AutoCorrelation_k2=gather(ACFull.rich.consumption.tperiods2.AutoCorrelation);
TM.AutoCorrTransProbs.poor.assets.StdDeviation=gather(ACFull.poor.assets.StdDeviation);
TM.AutoCorrTransProbs.rich.assets.PairMass=gather(ACFull.rich.assets.PairMass);
TM.AutoCorrTransProbs.poor.earnings.AutoCorrelation=NaN; % omitted
TM.CrossSectionCovarCorr.assets.CovarianceWith.earnings=gather(CSFull.assets.CovarianceWith.earnings);
TM.CrossSectionCovarCorr.earnings.CorrelationWith.consumption=gather(CSFull.earnings.CorrelationWith.consumption);
TM.CrossSectionCovarCorr.consumption.StdDeviation=gather(CSFull.consumption.StdDeviation);
TM.CrossSectionCovarCorr.consumption.Mean=gather(CSFull.consumption.Mean);
TM.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption=gather(CSFull.rich.assets.CorrelationWith.consumption);
t=gather(CSFull.poor.CovarianceMatrix); t(1,4)=NaN; t(4,1)=NaN; t(2,2)=NaN; TM.CrossSectionCovarCorr.poor.CovarianceMatrix=t; % matrix under poor, three entries omitted
t=gather(CSFull.CorrelationMatrix); t(logical(eye(4)))=NaN; TM.CrossSectionCovarCorr.CorrelationMatrix=t; % unrestricted matrix, the diagonal (always one) omitted

% The hand-built vector: a depth-first walk of TargetMoments in field order, reading the same field of the full run (the full-run
% structs have the same nesting as the targets, for restrictions, matrices and MoreInequality alike; an AutoCorr statistic with a _kK
% suffix is the command's .tperiodsK.(statistic)); then the target vector by the same walk
hand=[]; tvec=[];
for kk=1:3
    stack={{}}; % paths still to visit (last in, first out, children pushed in reverse so the walk is in field order)
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
            if kk==2 % AutoCorrTransProbs: stat_kK is under .tperiodsK
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
fprintf('MomentVector: %i target entries, %i of them targeted (the rest NaN); the model is NaN at a targeted entry, this should be zero: %.3e \n',length(tvec),sum(actual),sum(isnan(hand(actual))))

%% The objective, called as CalibrateInfHorzAgentModel calls it
CalibParamNames={'beta'};
calibparamsvec0=Params.beta;
calibparamsvecindex=[0;1];
calibomitparams_counter=0;
calibomitparamsmatrix={}; % (a cell since 2026-10-08, one omitted-values vector per masked parameter; empty here)
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
ran=0; wmsg='';
lastwarn('');
try
    [targetmomentvec,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,autocorrmomentnames,autocorrcummomentsizes,AutoCorrStats_whichstats,FnsToEvaluate_AutoCorr,crosssecmomentnames,crossseccummomentsizes,CrossSecStats_whichstats,FnsToEvaluate_CrossSec,cmsmomentnames,cmscummomentsizes,selectors,autocorrtimehorizons]=SetupTargetMoments_InfHorz(TM,FnsToEvaluate,0,simoptions);
    ran=1; wmsg=lastwarn;
catch err
    fprintf('MomentVector: SetupTargetMoments_InfHorz errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('MomentVector: SetupTargetMoments_InfHorz with simoptions, did not run, this should be zero: %.3e \n',1-ran)
fprintf('MomentVector: the target that is NaN in full (poor.earnings.AutoCorrelation) makes SetupTargetMoments_InfHorz warn (the warning names it), did not, this should be zero: %.3e \n',~contains(wmsg,'poor.earnings.AutoCorrelation is NaN in full'))
output=struct();
output.TargetMoments=TM;
output.hand=hand(actual);
if ran==0
    return
end
caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
caliboptions.selectors=selectors;
caliboptions.autocorrtimehorizons=autocorrtimehorizons;
argcell={CalibParamNames,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zg,ReturnFn,[],Params,DiscountFactorParamNames,[],FnsToEvaluate,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,targetmomentvec,allstatmomentnames,autocorrmomentnames,crosssecmomentnames,cmsmomentnames,allstatcummomentsizes,autocorrcummomentsizes,crossseccummomentsizes,cmscummomentsizes,AllStats_whichstats,AutoCorrStats_whichstats,CrossSecStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_AutoCorr,FnsToEvaluate_CrossSec,calibparamsvecindex,calibomitparams_counter,calibomitparamsmatrix};

if length(targetmomentvec)==length(tvec)
    d=abs(targetmomentvec-tvec); d=d(:); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(tvec));
else
    d=1; patternbad=abs(length(targetmomentvec)-length(tvec));
end
fprintf('MomentVector: the horizons read off the _kK suffixes are [%s] (hand [2 3]) \n',num2str(autocorrtimehorizons(:)'))
fprintf('MomentVector: SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e; horizons vs hand, this should be zero: %.3e \n',max([d;0]),patternbad,~isequal(autocorrtimehorizons(:)',[2,3]))

handA=hand(actual);
tA=targetmomentvec(actual);
for mm=1:2 % (caliboptions.whichcombos is always set by the driver, 1 or 0)
    caliboptions_m=caliboptions; caliboptions_m.whichcombos=2-mm;
    if mm==1
        lab='whichcombos=1';
    else
        lab='whichcombos=0';
    end
    ranm=0; m=[];
    try
        m=gather(CalibrateInfHorzAgentModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_m,vfoptions,simoptions));
        ranm=1;
    catch err
        fprintf('MomentVector (a) %s: the objective errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ranm==1 && length(m)==length(handA)
        d=CoreMaxAbs(m(:)-handA(:)); dt=CoreMaxAbs(m(:)-tA(:)); lenbad=0;
    else
        d=1; dt=1; lenbad=1;
    end
    fprintf('MomentVector (a) %s: did not run, this should be zero: %.3e; moment vector vs hand-built, this should be zero: %.3e; length mismatch, this should be zero: %.3e \n',lab,1-ranm,d,lenbad)
    fprintf('MomentVector (b) %s: moment vector vs the targets (targets are the model), this should be zero: %.3e \n',lab,dt)
    if mm==1
        output.m1=m;
    end
end

%% (c) the selectors handed straight to the three commands
% AllStats
simselA=simoptions; simselA.whichcombos=selectors.AllStats.whichcombos; simselA.whichstats=selectors.AllStats.whichstats;
AllSel=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate_AllStats,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simselA);
bad=0;
bad=bad+isfield(AllSel,'consumption'); % consumption is targeted only under poor
bad=bad+~isfield(AllSel.poor,'consumption');
if isfield(AllSel.poor,'consumption')
    bad=bad+isfield(AllSel.poor.consumption,'Median'); % only its Mean is targeted
end
bad=bad+isfield(AllSel.richind,'Gini'); % only richind's Mean is targeted
bad=bad+isfield(AllSel.rich,'consumption'); % nothing of consumption under rich
if isfield(AllSel.earnings,'Mean')
    bad=bad+~isnan(gather(AllSel.earnings.Mean)); % earnings unrestricted: only Lorenz/quantile targets (StatsFromWeightedGrid gives Minimum/Maximum with the quantiles, not the Mean)
end
fprintf('MomentVector (c) AllStats with the selectors: targeted assets Mean and rich earnings Mean vs full run, this should be zero: %.3e; untargeted structs present, this should be zero: %.3e \n',CoreMaxAbs([gather(AllSel.assets.Mean)-gather(AllFull.assets.Mean);gather(AllSel.rich.earnings.Mean)-gather(AllFull.rich.earnings.Mean)]),bad)
% AutoCorrTransProbs
simselC=simoptions; simselC.whichcombos=selectors.AutoCorr.whichcombos; simselC.timehorizons=autocorrtimehorizons;
ACSel=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate_AutoCorr,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zg,simselC);
bad=0;
bad=bad+isfield(ACSel,'richind'); % richind has no AutoCorr target (and is not in the reduced set)
bad=bad+isfield(ACSel.poor,'earnings'); % earnings under poor: the only target is NaN
bad=bad+isfield(ACSel,'consumption'); % consumption is targeted only under rich
bad=bad+~isfield(ACSel.rich,'consumption');
d=[gather(ACSel.assets.AutoCorrelation)-gather(ACFull.assets.AutoCorrelation); gather(ACSel.assets.tperiods3.AutoCorrelation)-gather(ACFull.assets.tperiods3.AutoCorrelation)];
if isfield(ACSel.rich,'consumption')
    d=[d; gather(ACSel.rich.consumption.tperiods2.AutoCorrelation)-gather(ACFull.rich.consumption.tperiods2.AutoCorrelation)];
end
fprintf('MomentVector (c) AutoCorrTransProbs with the selector and horizons: targeted assets AutoCorrelation (horizons 1 and 3) and rich consumption horizon 2 vs full run, this should be zero: %.3e; untargeted structs present, this should be zero: %.3e \n',CoreMaxAbs(d),bad)
% CrossSectionCovarCorr
simselX=simoptions; simselX.whichcombos=selectors.CrossSec.whichcombos;
CSSel=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate_CrossSec,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simselX);
bad=0;
% under rich only the assets-consumption pair is targeted: earnings is not evaluated on that page, so its struct is absent (the command
% creates a function's struct on a page only when one of its pairs or its own stats are wanted there; a pair of two evaluated functions
% that is not wanted is NaN); the run of 2026-10-09 expected a NaN pair here and crashed on the absent struct
bad=bad+isfield(CSSel.rich,'earnings');
if isfield(CSSel.rich.assets,'CovarianceWith')
    bad=bad+isfield(CSSel.rich.assets.CovarianceWith,'earnings');
end
Mp=gather(CSSel.poor.CovarianceMatrix);
bad=bad+~isnan(Mp(1,4))+~isnan(Mp(4,1)); % the omitted off-diagonal entries stay NaN; the omitted diagonal entry (2,2) is NOT NaN: earnings is evaluated on the poor page for its other pairs, and the command reports an evaluated function's own variance on the diagonal regardless of the selector (its documented byproduct rule; the run of 2026-10-09 expected NaN there)
d=[gather(CSSel.assets.CovarianceWith.earnings)-gather(CSFull.assets.CovarianceWith.earnings); gather(CSSel.rich.assets.CorrelationWith.consumption)-gather(CSFull.rich.assets.CorrelationWith.consumption); Mp(1,2)-gather(CSFull.poor.CovarianceMatrix(1,2)); Mp(3,3)-gather(CSFull.poor.CovarianceMatrix(3,3))];
Mc=gather(CSSel.CorrelationMatrix); Mcf=gather(CSFull.CorrelationMatrix);
d=[d; Mc(~eye(4))-Mcf(~eye(4))]; % every off-diagonal pair is targeted unrestricted
fprintf('MomentVector (c) CrossSectionCovarCorr with the pair selector: targeted pairs and matrix entries vs full run, this should be zero: %.3e; an unevaluated function present or omitted matrix entries not NaN, this should be zero: %.3e \n',CoreMaxAbs(d),bad)

end
