function output=CalibFHorz_AutoCorrWindows(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% CalibFHorzTests part 9: AutoCorrTransProbs targets of the window and log rows (added 2026-10-14), fnname_log,
% fnname_MperiodAverage, fnname_MperiodSum and fnname_MperiodAverage_log/_MperiodSum_log, the rows of the command's
% simoptions.ntimeperiods and simoptions.logwhichcombos. SetupTargetMoments_FHorz reads them off the target names and puts
% simoptions.ntimeperiods, .ntimeperiodstype and .logwhichcombos in selectors.AutoCorr, and the objective passes them on.
%   (a) parsing, no model solved: selectors.AutoCorr (whichcombos with rows in blocks of the window lengths, logwhichcombos,
%       ntimeperiods, ntimeperiodstype) and the horizons against hand-built ones; four invalid target sets error (MperiodAverage
%       and MperiodSum together, a statistic and its log at the same start age, a 1-period window, a name that is no row)
%   (b) one solve at the truth, full runs of the command at the level and at the log; targets from them with NaN omissions;
%       the objective's moment vector (vectoroutput=1) under caliboptions.whichcombos=1 and =0 equals the hand walk of the
%       targets through the full runs and the targets themselves (the targets are the model)
% The log of a 2-period window needs horizons of at least 2 (the command refuses overlapping windows), so no target here is at
% horizon 1. Each objective run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.inc=@(aprime,a,z,w,kappa_j) 1+a+w*kappa_j*z; % positive at every state, so its log and the log of its windows exist
Jr=Params.Jr;

simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);

%% (a) parsing, no model solved
TM=struct();
t=ones(1,N_j); t(Jr:end)=NaN; TM.AutoCorrTransProbs.earnings.Mean=t;                                 % earnings (M=1), ages 1..Jr-1
t=nan(1,N_j-2); t(1:Jr-4)=1; TM.AutoCorrTransProbs.earnings_2periodAverage.AutoCorrelation_horizon2=t; % consecutive 2-period windows, start ages 1..Jr-4
t=nan(1,N_j-4); t(2:10)=1; TM.AutoCorrTransProbs.inc_2periodAverage_log.AutoCorrelation_horizon4=t;   % the log of the 2-period average, start ages 2..10
t=ones(1,N_j); t(N_j)=NaN; TM.AutoCorrTransProbs.inc_2periodAverage_log.Mean=t;                        % its mean, ages 1..N_j-1
TM.AutoCorrTransProbs.rich.inc_log.StdDeviation=ones(1,N_j);                                           % rich, the log of inc (M=1), every age
TM.AutoCorrTransProbs.rich.assets_2periodAverage.PairMass_horizon2=ones(1,N_j-2);                     % rich, 2-period windows of assets
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors,~,~,~,FnsToEvaluate_AutoCorr,autocorrtimehorizons]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,0,N_j,simoptions);
% hand-built: the functions assets, earnings, inc (FnsToEvaluate order), rows in blocks of ntimeperiods=[1,2], pages unrestricted and rich
wc=zeros(6,N_j,2); lw=zeros(6,N_j,2);
wc(2,1:Jr-1,1)=1;            % earnings, M=1: Mean
wc(5,1:Jr-4,1)=1;            % earnings, M=2: AutoCorrelation_horizon2
wc(6,1:N_j-1,1)=1;           % inc, M=2: the log (Mean ages 1..N_j-1, which include the AutoCorrelation_horizon4 start ages 2..10)
lw(6,1:N_j-1,1)=1;
wc(3,:,2)=1; lw(3,:,2)=1;    % rich, inc, M=1: the log
wc(4,1:N_j-2,2)=1;           % rich, assets, M=2: PairMass_horizon2
s=selectors.AutoCorr;
bad=~isequal(fieldnames(FnsToEvaluate_AutoCorr),{'assets';'earnings';'inc'})+~isequal(s.ntimeperiods,[1,2])+~strcmp(s.ntimeperiodstype,'average')+~isequal(autocorrtimehorizons(:)',[2,4]);
if isequal(size(s.whichcombos),size(wc)) && isequal(size(s.logwhichcombos),size(lw))
    d=max([abs(s.whichcombos(:)-wc(:));abs(s.logwhichcombos(:)-lw(:))]);
else
    d=1;
end
fprintf('AutoCorrWindows (a) reduced FnsToEvaluate, ntimeperiods [%s] (hand [1 2]), ntimeperiodstype %s (hand average), horizons [%s] (hand [2 4]): mismatches, this should be zero: %i; selectors whichcombos and logwhichcombos vs hand, this should be zero: %.3e \n',num2str(s.ntimeperiods),s.ntimeperiodstype,num2str(autocorrtimehorizons(:)'),bad,d)
% invalid target sets
errcases={'MperiodAverage and MperiodSum together','a statistic and its log at the same start age','a 1-period window','a name that is no row'};
errmsgs={'both MperiodAverage and MperiodSum','both a statistic and its log','single period','is not one of the FnsToEvaluate'};
errok=zeros(1,4);
for ee=1:4
    TMbad=TM;
    switch ee
        case 1
            TMbad.AutoCorrTransProbs.earnings_2periodSum.Mean=ones(1,N_j);
        case 2
            t=nan(1,N_j); t(3:5)=1; TMbad.AutoCorrTransProbs.inc_2periodAverage.Mean=t; % the level at ages where the log is targeted
        case 3
            TMbad.AutoCorrTransProbs.earnings_1periodAverage.Mean=ones(1,N_j);
        case 4
            TMbad.AutoCorrTransProbs.rich.nosuchfn_2periodAverage.Mean=ones(1,N_j);
    end
    try
        SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
    catch ME
        errok(ee)=contains(ME.message,errmsgs{ee});
        if ~errok(ee)
            fprintf('   (the error was: %s) \n',ME.message);
        end
    end
    fprintf('AutoCorrWindows (a) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok(ee));
end

%% (b) one solve at the truth, full runs at the level and at the log, the objective's moment vector
vfoptions=struct();
[z_gridvals_J,pi_z_J,vfoptions]=ExogShockSetup_FHorz(n_z,z_grid,pi_z,N_j,Params,vfoptions,3,0);
simoptions.e_gridvals_J=vfoptions.e_gridvals_J;
simoptions.pi_e_J=vfoptions.pi_e_J;
vfoptions.alreadygridvals=1;
simoptions.alreadygridvals=1;
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
simoptions_L=simoptions; simoptions_L.ntimeperiods=[1,2]; simoptions_L.timehorizons=[2,4];
ACL=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,simoptions_L); % the level rows
FnsInc.inc=FnsToEvaluate.inc;
simoptions_G=simoptions_L; simoptions_G.logwhichcombos=1;
ACG=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsInc,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,simoptions_G); % the log rows of inc
ACFull=ACL; % the level rows and the log rows together
ACFull.inc_log=ACG.inc_log;
ACFull.inc_2periodAverage_log=ACG.inc_2periodAverage_log;
ACFull.rich.inc_log=ACG.rich.inc_log;
% the targets: the structure of (a), the values from the full runs, NaN where omitted
TMb=struct();
t=gather(ACFull.earnings.Mean); t(Jr:end)=NaN; TMb.AutoCorrTransProbs.earnings.Mean=t;
t=gather(ACFull.earnings_2periodAverage.AutoCorrelation_horizon2); t(Jr-3:end)=NaN; TMb.AutoCorrTransProbs.earnings_2periodAverage.AutoCorrelation_horizon2=t;
t=gather(ACFull.inc_2periodAverage_log.AutoCorrelation_horizon4); t([1,11:end])=NaN; TMb.AutoCorrTransProbs.inc_2periodAverage_log.AutoCorrelation_horizon4=t;
TMb.AutoCorrTransProbs.inc_2periodAverage_log.Mean=gather(ACFull.inc_2periodAverage_log.Mean);
TMb.AutoCorrTransProbs.rich.inc_log.StdDeviation=gather(ACFull.rich.inc_log.StdDeviation);
TMb.AutoCorrTransProbs.rich.assets_2periodAverage.PairMass_horizon2=gather(ACFull.rich.assets_2periodAverage.PairMass_horizon2);
% the hand walk: the same fields of the full runs, in the target order
hand=[]; tvec=[];
a1vec=fieldnames(TMb.AutoCorrTransProbs);
for a1=1:length(a1vec)
    a2vec=fieldnames(TMb.AutoCorrTransProbs.(a1vec{a1}));
    for a2=1:length(a2vec)
        t2=TMb.AutoCorrTransProbs.(a1vec{a1}).(a2vec{a2});
        if ~isstruct(t2)
            v=gather(ACFull.(a1vec{a1}).(a2vec{a2})); hand=[hand; v(:)]; tvec=[tvec; t2(:)]; %#ok<AGROW>
        else
            a3vec=fieldnames(t2);
            for a3=1:length(a3vec)
                v=gather(ACFull.(a1vec{a1}).(a2vec{a2}).(a3vec{a3})); t3=t2.(a3vec{a3}); hand=[hand; v(:)]; tvec=[tvec; t3(:)]; %#ok<AGROW>
            end
        end
    end
end
actual=~isnan(tvec);
fprintf('AutoCorrWindows (b) %i target entries, %i of them targeted; the model is NaN at a targeted entry, this should be zero: %i \n',length(tvec),sum(actual),sum(isnan(hand(actual))))

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
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TMb,FnsToEvaluate,0,N_j,simoptions);
caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
caliboptions.selectors=selectors;
argcell={CalibParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,[],Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,[],usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,calibparamsvecindex,calibomitparams_counter,calibomitparamsmatrix};
if length(targetmomentvec)==length(tvec)
    d=abs(targetmomentvec-tvec); d=d(:); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(tvec));
else
    d=1; patternbad=abs(length(targetmomentvec)-length(tvec));
end
fprintf('AutoCorrWindows (b) SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %i \n',CoreMaxAbs(d),patternbad)
handA=hand(actual); tA=targetmomentvec(actual);
for wcc=[1,0]
    caliboptions_t=caliboptions; caliboptions_t.whichcombos=wcc;
    ran=0; m=[];
    try
        m=gather(CalibrateLifeCycleModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions));
        ran=1;
    catch err
        fprintf('AutoCorrWindows (b) whichcombos=%i: the objective errored with: %s \n',wcc,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ran==1 && length(m)==length(handA)
        d=CoreMaxAbs(m(:)-handA(:)); dt=CoreMaxAbs(m(:)-tA(:)); nanbad=sum(isnan(m));
    else
        d=1; dt=1; nanbad=1;
    end
    fprintf('AutoCorrWindows (b) whichcombos=%i: did not run, this should be zero: %i; moment vector vs hand-built, this should be zero: %.3e; vs the targets, this should be zero: %.3e; NaN entries or length mismatch, this should be zero: %i \n',wcc,1-ran,d,dt,nanbad)
end

output.selectors=selectors;
end
