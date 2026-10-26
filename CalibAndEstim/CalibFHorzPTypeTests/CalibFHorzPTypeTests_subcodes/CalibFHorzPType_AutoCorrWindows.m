function output=CalibFHorzPType_AutoCorrWindows(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% CalibFHorzPTypeTests part 9: AutoCorrTransProbs targets of the window and log rows (added 2026-10-14) with permanent types,
% grouped (.(row).(stat)) and per type (.(row).(type).(stat)): fnname_log, fnname_MperiodAverage, fnname_MperiodAverage_log, ...,
% the rows of simoptions.ntimeperiods and simoptions.logwhichcombos. SetupTargetMoments_FHorz puts simoptions.ntimeperiods,
% .ntimeperiodstype and .logwhichcombos in selectors.AutoCorr (logwhichcombos without the type dimension: the PType command applies
% it to every type and to the grouped stats), and the objective passes them on.
%   (a) parsing, no model solved: selectors.AutoCorr (whichcombos [6,N_j,2,N_i+1], logwhichcombos [6,N_j,2], ntimeperiods,
%       ntimeperiodstype) and the horizons against hand-built ones; four invalid target sets error (MperiodAverage and MperiodSum
%       together; the level for one type where the log is targeted for another, at the same start age; a 1-period window; a name
%       that is no row)
%   (b) one solve at the truth, full runs of the PType command at the level and at the log; targets from them with NaN omissions;
%       the objective's moment vector (vectoroutput=1) under caliboptions.whichcombos=1 and =0 equals the hand walk and the targets
% No target is at horizon 1 (the log of a 2-period window refuses overlapping windows). Each objective run is wrapped.

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate.inc=@(aprime,a,z,w,kappa_j_pt) 1+a+w*kappa_j_pt*z; % positive at every state
Jr=Params.Jr;
ilow=find(strcmp(Names_i,'low')); ihigh=find(strcmp(Names_i,'high')); igrp=N_i+1;

simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);

%% (a) parsing, no model solved
TM=struct();
t=ones(1,N_j); t(Jr:end)=NaN; TM.AutoCorrTransProbs.earnings.Mean=t;                                      % grouped, M=1
t=nan(1,N_j-2); t(1:Jr-4)=1; TM.AutoCorrTransProbs.earnings_2periodAverage.low.AutoCorrelation_horizon2=t; % low, 2-period windows
t=nan(1,N_j-4); t(2:10)=1; TM.AutoCorrTransProbs.inc_2periodAverage_log.AutoCorrelation_horizon4=t;        % grouped, the log of the 2-period average
t=ones(1,N_j); t(N_j)=NaN; TM.AutoCorrTransProbs.inc_2periodAverage_log.high.Mean=t;                        % high, its mean
TM.AutoCorrTransProbs.rich.inc_log.StdDeviation=ones(1,N_j);                                                % grouped, rich, the log of inc
TM.AutoCorrTransProbs.rich.assets_2periodAverage.low.PairMass_horizon2=ones(1,N_j-2);                      % low, rich, 2-period windows of assets
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors,~,~,~,FnsToEvaluate_AutoCorr,autocorrtimehorizons]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,1,N_j,simoptions,Names_i);
% hand-built: assets, earnings, inc (FnsToEvaluate order), rows in blocks of ntimeperiods=[1,2], pages unrestricted and rich, then the type slots
wc=zeros(6,N_j,2,N_i+1); lw=zeros(6,N_j,2);
wc(2,1:Jr-1,1,igrp)=1;          % earnings, M=1, grouped Mean
wc(5,1:Jr-4,1,ilow)=1;          % earnings, M=2, low AutoCorrelation_horizon2
wc(6,2:10,1,igrp)=1;            % inc, M=2, the log: grouped AutoCorrelation_horizon4
wc(6,1:N_j-1,1,ihigh)=1;        % and high Mean
lw(6,1:N_j-1,1)=1;              % the log combinations (no type dimension)
wc(3,:,2,igrp)=1; lw(3,:,2)=1;  % rich, inc, M=1, the log
wc(4,1:N_j-2,2,ilow)=1;         % rich, assets, M=2, low PairMass_horizon2
s=selectors.AutoCorr;
bad=~isequal(fieldnames(FnsToEvaluate_AutoCorr),{'assets';'earnings';'inc'})+~isequal(s.ntimeperiods,[1,2])+~strcmp(s.ntimeperiodstype,'average')+~isequal(autocorrtimehorizons(:)',[2,4]);
if isequal(size(s.whichcombos),size(wc)) && isequal(size(s.logwhichcombos),size(lw))
    d=max([abs(s.whichcombos(:)-wc(:));abs(s.logwhichcombos(:)-lw(:))]);
else
    d=1;
end
fprintf('AutoCorrWindows PType (a) reduced FnsToEvaluate, ntimeperiods [%s] (hand [1 2]), ntimeperiodstype %s (hand average), horizons [%s] (hand [2 4]): mismatches, this should be zero: %i; selectors whichcombos and logwhichcombos vs hand, this should be zero: %.3e \n',num2str(s.ntimeperiods),s.ntimeperiodstype,num2str(autocorrtimehorizons(:)'),bad,d)
errcases={'MperiodAverage and MperiodSum together','the level for one type where another has the log, at the same start age','a 1-period window','a name that is no row'};
errmsgs={'both MperiodAverage and MperiodSum','both a statistic and its log','single period','is not one of the FnsToEvaluate'};
for ee=1:4
    TMbad=TM;
    switch ee
        case 1
            TMbad.AutoCorrTransProbs.earnings_2periodSum.low.Mean=ones(1,N_j);
        case 2
            t=nan(1,N_j); t(3:5)=1; TMbad.AutoCorrTransProbs.inc_2periodAverage.low.Mean=t; % the level for low where high has the log
        case 3
            TMbad.AutoCorrTransProbs.earnings_1periodAverage.Mean=ones(1,N_j);
        case 4
            TMbad.AutoCorrTransProbs.rich.nosuchfn_2periodAverage.Mean=ones(1,N_j);
    end
    errok=0;
    try
        SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
    catch ME
        errok=contains(ME.message,errmsgs{ee});
        if ~errok
            fprintf('   (the error was: %s) \n',ME.message);
        end
    end
    fprintf('AutoCorrWindows PType (a) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok);
end

%% (b) one solve at the truth, full runs at the level and the log, the objective's moment vector
vfoptions=struct();
[z_gridvals_J,pi_z_J,vfoptions]=ExogShockSetup_FHorz_PType(n_z,z_grid,pi_z,N_j,Names_i,Params,vfoptions,3,0);
simoptions.e_gridvals_J=vfoptions.e_gridvals_J;
simoptions.pi_e_J=vfoptions.pi_e_J;
vfoptions.alreadygridvals=1;
simoptions.alreadygridvals=1;
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params,simoptions);
simoptions_L=simoptions; simoptions_L.ntimeperiods=[1,2]; simoptions_L.timehorizons=[2,4];
ACL=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,simoptions_L);
FnsInc.inc=FnsToEvaluate.inc;
simoptions_G=simoptions_L; simoptions_G.logwhichcombos=1;
ACG=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsInc,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,simoptions_G);
ACFull=ACL; % the level rows and the log rows together (each with its grouped and per-type outputs)
ACFull.inc_log=ACG.inc_log;
ACFull.inc_2periodAverage_log=ACG.inc_2periodAverage_log;
ACFull.rich.inc_log=ACG.rich.inc_log;
TMb=struct();
t=gather(ACFull.earnings.Mean); t(Jr:end)=NaN; TMb.AutoCorrTransProbs.earnings.Mean=t;
t=gather(ACFull.earnings_2periodAverage.low.AutoCorrelation_horizon2); t(Jr-3:end)=NaN; TMb.AutoCorrTransProbs.earnings_2periodAverage.low.AutoCorrelation_horizon2=t;
t=gather(ACFull.inc_2periodAverage_log.AutoCorrelation_horizon4); t([1,11:end])=NaN; TMb.AutoCorrTransProbs.inc_2periodAverage_log.AutoCorrelation_horizon4=t;
TMb.AutoCorrTransProbs.inc_2periodAverage_log.high.Mean=gather(ACFull.inc_2periodAverage_log.high.Mean);
TMb.AutoCorrTransProbs.rich.inc_log.StdDeviation=gather(ACFull.rich.inc_log.StdDeviation);
TMb.AutoCorrTransProbs.rich.assets_2periodAverage.low.PairMass_horizon2=gather(ACFull.rich.assets_2periodAverage.low.PairMass_horizon2);
% the hand walk (two to four names), in the target order
hand=[]; tvec=[];
a1vec=fieldnames(TMb.AutoCorrTransProbs);
for a1=1:length(a1vec)
    t1=TMb.AutoCorrTransProbs.(a1vec{a1});
    a2vec=fieldnames(t1);
    for a2=1:length(a2vec)
        t2=t1.(a2vec{a2});
        if ~isstruct(t2)
            v=gather(ACFull.(a1vec{a1}).(a2vec{a2})); hand=[hand; v(:)]; tvec=[tvec; t2(:)]; %#ok<AGROW>
        else
            a3vec=fieldnames(t2);
            for a3=1:length(a3vec)
                t3=t2.(a3vec{a3});
                if ~isstruct(t3)
                    v=gather(ACFull.(a1vec{a1}).(a2vec{a2}).(a3vec{a3})); hand=[hand; v(:)]; tvec=[tvec; t3(:)]; %#ok<AGROW>
                else
                    a4vec=fieldnames(t3);
                    for a4=1:length(a4vec)
                        v=gather(ACFull.(a1vec{a1}).(a2vec{a2}).(a3vec{a3}).(a4vec{a4})); t4=t3.(a4vec{a4}); hand=[hand; v(:)]; tvec=[tvec; t4(:)]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
actual=~isnan(tvec);
fprintf('AutoCorrWindows PType (b) %i target entries, %i of them targeted; the model is NaN at a targeted entry, this should be zero: %i \n',length(tvec),sum(actual),sum(isnan(hand(actual))))

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
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TMb,FnsToEvaluate,1,N_j,simoptions,Names_i);
caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
caliboptions.selectors=selectors;
argcell={CalibParamNames,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,nCalibParams,nCalibParamsFinder,calibparamsvecindex,calibparamssizes,calibomitparams_counter,calibomitparamsmatrix};
if length(targetmomentvec)==length(tvec)
    d=abs(targetmomentvec-tvec); d=d(:); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(tvec));
else
    d=1; patternbad=abs(length(targetmomentvec)-length(tvec));
end
fprintf('AutoCorrWindows PType (b) SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %i \n',CoreMaxAbs(d),patternbad)
handA=hand(actual); tA=targetmomentvec(actual);
for wcc=[1,0]
    caliboptions_t=caliboptions; caliboptions_t.whichcombos=wcc;
    ran=0; m=[];
    try
        m=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions));
        ran=1;
    catch err
        fprintf('AutoCorrWindows PType (b) whichcombos=%i: the objective errored with: %s \n',wcc,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ran==1 && length(m)==length(handA)
        d=CoreMaxAbs(m(:)-handA(:)); dt=CoreMaxAbs(m(:)-tA(:)); nanbad=sum(isnan(m));
    else
        d=1; dt=1; nanbad=1;
    end
    fprintf('AutoCorrWindows PType (b) whichcombos=%i: did not run, this should be zero: %i; moment vector vs hand-built, this should be zero: %.3e; vs the targets, this should be zero: %.3e; NaN entries or length mismatch, this should be zero: %i \n',wcc,1-ran,d,dt,nanbad)
end

output.selectors=selectors;
end
