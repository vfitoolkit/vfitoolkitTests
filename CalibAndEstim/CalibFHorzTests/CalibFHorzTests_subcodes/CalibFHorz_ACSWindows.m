function output=CalibFHorz_ACSWindows(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% CalibFHorzTests part 10: AgeConditionalStats targets of the window and log rows (added 2026-10-14), fnname_log,
% fnname_MperiodAverage, fnname_MperiodSum and fnname_MperiodAverage_log/_MperiodSum_log, the rows of LifeCycleProfiles'
% simoptions.ntimeperiods and simoptions.logwhichcombos (a window target is by start age, 1 x N_j). SetupTargetMoments_FHorz reads them
% off the target names into selectors.ACStats (rows in blocks of ntimeperiods, logwhichcombos, ntimeperiodstype, agegroupings) and,
% when the single-period targets are by age groups, into a second call selectors.ACStatsW; the objective makes the call(s).
%   (a) parsing, no model solved: the selectors against hand-built ones; the split with agegroupings=1:5:N_j, and windows only with
%       agegroupings=1:5:N_j (forced to 1:N_j); five invalid target sets error (MperiodAverage and MperiodSum together, a statistic and
%       its log at the same age, a 1-period window, a name that is no row, a window targeted at a start age where it runs past N_j)
%   (b) one solve at the truth, full LifeCycleProfiles runs at the level and at the log; targets from them with NaN omissions; the
%       objective's moment vector (vectoroutput=1) under caliboptions.whichcombos=1 and =0 equals the hand walk of the targets through
%       the full runs and the targets themselves; with the default agegroupings (one call) and with agegroupings=1:5:N_j (two calls)
%   (c) LifeCycleProfiles itself: a restricted 3-period window that every path leaves (mid-window and at its last age) gives
%       RestrictedSampleMass 0 and NaN stats, with the unrestricted page in the same list and without it
% Each objective run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.inc=@(aprime,a,z,w,kappa_j) 1+a+w*kappa_j*z; % positive at every state, so its log and the log of its windows exist
Jr=Params.Jr;

simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);

%% (a) parsing, no model solved
TM=struct();
t=ones(1,N_j); t(Jr:end)=NaN; TM.AgeConditionalStats.earnings.Mean=t;                                        % earnings (M=1), ages 1..Jr-1
t=nan(1,N_j); t(1:Jr-3)=1; TM.AgeConditionalStats.earnings_2periodAverage.MoreInequality.Percentile90th=t;    % the 90th percentile of the 2-period average, start ages 1..Jr-3
t=nan(1,N_j); t(2:10)=1; TM.AgeConditionalStats.inc_2periodAverage_log.Mean=t;                               % the log of the 2-period average, start ages 2..10
TM.AgeConditionalStats.rich.inc_log.StdDeviation=ones(1,N_j);                                                 % rich, the log of inc (M=1), every age
t=ones(1,N_j); t(N_j-1:N_j)=NaN; TM.AgeConditionalStats.rich.assets_3periodAverage.Median=t;                 % rich, 3-period averages of assets, start ages 1..N_j-2
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,0,N_j,simoptions);
[~,~,~,~,~,~,~,~,~,~,~,FnsToEvaluate_ACStats]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,0,N_j,simoptions);
% hand-built: the functions assets, earnings, inc, rows in blocks of ntimeperiods=[1,2,3], pages unrestricted and rich
ws=zeros(9,N_j,2,7); lw=zeros(9,N_j,2);
ws(2,1:Jr-1,1,1)=1;           % earnings, M=1: Mean
ws(5,1:Jr-3,1,7)=1;           % earnings, M=2: MoreInequality
ws(6,2:10,1,1)=1; lw(6,2:10,1)=1;   % inc, M=2: the log, Mean
ws(3,:,2,3)=1; lw(3,:,2)=1;   % rich, inc, M=1: the log, StdDeviation
ws(7,1:N_j-2,2,2)=1;          % rich, assets, M=3: Median
s=selectors.ACStats;
bad=~isequal(fieldnames(FnsToEvaluate_ACStats),{'assets';'earnings';'inc'})+~(isfield(s,'haswindowlog') && s.haswindowlog==1);
if bad==0
    bad=bad+~isequal(s.ntimeperiods,[1,2,3])+~strcmp(s.ntimeperiodstype,'average')+~isempty(s.agegroupings)+(s.split~=0)+isfield(selectors,'ACStatsW');
end
if bad==0 && isequal(size(s.whichstats),size(ws)) && isequal(size(s.logwhichcombos),size(lw)) && isequal(size(s.whichcombos),size(lw))
    d=max([abs(s.whichstats(:)-ws(:));abs(s.logwhichcombos(:)-lw(:));abs(s.whichcombos(:)-reshape(any(ws,4),[],1))]);
else
    d=1;
end
fprintf('ACSWindows (a) reduced FnsToEvaluate, haswindowlog, ntimeperiods [1 2 3], average, no agegroupings, no split: mismatches, this should be zero: %i; selectors whichstats, whichcombos and logwhichcombos vs hand, this should be zero: %.3e \n',bad,d)
% the split: the single-period target by the age groups 1:5:N_j (four groups), the windows by start age in a second call
simoptions_g=simoptions; simoptions_g.agegroupings=1:5:N_j;
TMg=TM;
TMg.AgeConditionalStats.earnings.Mean=[1,1,1,NaN];
TMg.AgeConditionalStats.rich.inc_log.StdDeviation=ones(1,4);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,selg]=SetupTargetMoments_FHorz(TMg,FnsToEvaluate,0,N_j,simoptions_g);
bad=0;
if isfield(selg.ACStats,'split') && selg.ACStats.split==1 && isfield(selg,'ACStatsW')
    sA=selg.ACStats; sW=selg.ACStatsW;
    bad=bad+~isequal(sA.ntimeperiods,1)+~isempty(sA.agegroupings)+~isequal(size(sA.whichstats),[3,4,2,7])+~isequal(sW.ntimeperiods,[2,3])+~isequal(sW.agegroupings,1:N_j)+~isequal(size(sW.whichstats),[6,N_j,2,7]);
    if bad==0
        bad=bad+~isequal(reshape(sA.whichstats(2,:,1,1),1,[]),[1,1,1,0])+~isequal(reshape(sA.logwhichcombos(3,:,2),1,[]),[1,1,1,1])+~isequal(reshape(sW.whichstats(2,:,1,7),1,[]),reshape(ws(5,:,1,7),1,[]))+~isequal(reshape(sW.logwhichcombos(3,:,1),1,[]),reshape(lw(6,:,1),1,[]))+~isequal(reshape(sW.whichstats(4,:,2,2),1,[]),reshape(ws(7,:,2,2),1,[]));
    end
else
    bad=1;
end
fprintf('ACSWindows (a) agegroupings 1:5:N_j with single-period and window targets: the split into ACStats (age groups) and ACStatsW (start ages), mismatches vs hand, this should be zero: %i \n',bad)
% windows only, with agegroupings 1:5:N_j: one call, by start age
TMw=struct();
TMw.AgeConditionalStats.earnings_2periodAverage=TM.AgeConditionalStats.earnings_2periodAverage;
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,selw]=SetupTargetMoments_FHorz(TMw,FnsToEvaluate,0,N_j,simoptions_g);
bad=~(isfield(selw.ACStats,'split') && selw.ACStats.split==0 && isequal(selw.ACStats.agegroupings,1:N_j) && isequal(selw.ACStats.ntimeperiods,2) && isequal(size(selw.ACStats.whichstats),[1,N_j,2,7]));
fprintf('ACSWindows (a) windows only with agegroupings 1:5:N_j: one call with agegroupings forced to 1:N_j, mismatches, this should be zero: %i \n',bad)
% invalid target sets
errcases={'MperiodAverage and MperiodSum together','a statistic and its log at the same age','a 1-period window','a name that is no row','a window targeted where it runs past N_j'};
errmsgs={'both MperiodAverage and MperiodSum','both a statistic and its log','single period','is not one of the FnsToEvaluate','runs past age N_j'};
for ee=1:5
    TMbad=TM;
    switch ee
        case 1
            TMbad.AgeConditionalStats.earnings_2periodSum.Mean=ones(1,N_j-1);
        case 2
            t=nan(1,N_j); t(3:5)=1; TMbad.AgeConditionalStats.inc_2periodAverage.Mean=t; % the level at ages where the log is targeted
        case 3
            TMbad.AgeConditionalStats.earnings_1periodAverage.Mean=ones(1,N_j);
        case 4
            TMbad.AgeConditionalStats.rich.nosuchfn_2periodAverage.Mean=ones(1,N_j);
        case 5
            TMbad.AgeConditionalStats.assets_3periodAverage.Mean=ones(1,N_j); % start ages N_j-1 and N_j run past N_j
    end
    errok=0;
    try
        SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
    catch ME
        errok=contains(ME.message,errmsgs{ee});
        if ~errok
            fprintf('   (the error was: %s) \n',ME.message);
        end
    end
    fprintf('ACSWindows (a) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok);
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
simoptions_L=simoptions; simoptions_L.ntimeperiods=[1,2,3]; simoptions_L.pi_z=pi_z_J; simoptions_L.whichstats=ones(1,7);
ACL=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_L); % the level rows
FnsInc.inc=FnsToEvaluate.inc;
simoptions_G=simoptions_L; simoptions_G.logwhichcombos=1;
ACG=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsInc,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_G); % the log rows of inc
simoptions_A=simoptions; simoptions_A.agegroupings=1:5:N_j; simoptions_A.whichstats=ones(1,7);
ACA=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_A); % by age groups
FnsIncA.inc=FnsToEvaluate.inc;
simoptions_AG=simoptions_A; simoptions_AG.logwhichcombos=1;
ACAG=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsIncA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_AG); % the log of inc, by age groups
for scen=1:2
    ACFull=ACL; % the level rows and the log rows together
    ACFull.inc_log=ACG.inc_log;
    ACFull.inc_2periodAverage_log=ACG.inc_2periodAverage_log;
    ACFull.rich.inc_log=ACG.rich.inc_log;
    simoptions_s=simoptions;
    if scen==1
        scenstr='one call (default agegroupings)';
    else
        scenstr='two calls (agegroupings 1:5:N_j)';
        ACFull.earnings=ACA.earnings; % the single-period rows by age groups
        ACFull.rich.inc_log=ACAG.rich.inc_log;
        simoptions_s.agegroupings=1:5:N_j;
    end
    % the targets: the structure of (a), the values from the full runs, NaN where omitted
    TMb=struct();
    t=gather(ACFull.earnings.Mean);
    if scen==1
        t(Jr:end)=NaN;
    else
        t(4)=NaN;
    end
    TMb.AgeConditionalStats.earnings.Mean=t;
    t=gather(ACFull.earnings_2periodAverage.MoreInequality.Percentile90th); t(Jr-2:end)=NaN; TMb.AgeConditionalStats.earnings_2periodAverage.MoreInequality.Percentile90th=t;
    t=gather(ACFull.inc_2periodAverage_log.Mean); t([1,11:end])=NaN; TMb.AgeConditionalStats.inc_2periodAverage_log.Mean=t;
    TMb.AgeConditionalStats.rich.inc_log.StdDeviation=gather(ACFull.rich.inc_log.StdDeviation);
    t=gather(ACFull.rich.assets_3periodAverage.Median); t(N_j-1:N_j)=NaN; TMb.AgeConditionalStats.rich.assets_3periodAverage.Median=t;
    % the hand walk: the same fields of the full runs, in the target order
    hand=[]; tvec=[];
    a1vec=fieldnames(TMb.AgeConditionalStats);
    for a1=1:length(a1vec)
        a2vec=fieldnames(TMb.AgeConditionalStats.(a1vec{a1}));
        for a2=1:length(a2vec)
            t2=TMb.AgeConditionalStats.(a1vec{a1}).(a2vec{a2});
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
    fprintf('ACSWindows (b) %s: %i target entries, %i of them targeted; the model is NaN at a targeted entry, this should be zero: %i \n',scenstr,length(tvec),sum(actual),sum(isnan(hand(actual))))

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
    [targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TMb,FnsToEvaluate,0,N_j,simoptions_s);
    caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
    caliboptions.selectors=selectors;
    argcell={CalibParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,[],Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,[],usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,calibparamsvecindex,calibomitparams_counter,calibomitparamsmatrix};
    if length(targetmomentvec)==length(tvec)
        d=abs(targetmomentvec-tvec); d=d(:); d=d(~isnan(d));
        patternbad=sum(isnan(targetmomentvec)~=isnan(tvec));
    else
        d=1; patternbad=abs(length(targetmomentvec)-length(tvec));
    end
    fprintf('ACSWindows (b) %s: SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %i; split as expected, this should be zero: %i \n',scenstr,CoreMaxAbs(d),patternbad,~(selectors.ACStats.split==(scen==2)))
    handA=hand(actual); tA=targetmomentvec(actual);
    for wcc=[1,0]
        caliboptions_t=caliboptions; caliboptions_t.whichcombos=wcc;
        ran=0; m=[];
        try
            m=gather(CalibrateLifeCycleModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions_s));
            ran=1;
        catch err
            fprintf('ACSWindows (b) %s, whichcombos=%i: the objective errored with: %s \n',scenstr,wcc,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        end
        if ran==1 && length(m)==length(handA)
            d=CoreMaxAbs(m(:)-handA(:)); dt=CoreMaxAbs(m(:)-tA(:)); nanbad=sum(isnan(m));
        else
            d=1; dt=1; nanbad=1;
        end
        fprintf('ACSWindows (b) %s, whichcombos=%i: did not run, this should be zero: %i; moment vector vs hand-built, this should be zero: %.3e; vs the targets, this should be zero: %.3e; NaN entries or length mismatch, this should be zero: %i \n',scenstr,wcc,1-ran,d,dt,nanbad)
    end
end

%% (c) a restricted window that every path leaves (added 2026-11-07)
% young (agej<Jr) holds at every state before Jr and at none from Jr on. A 3-period window of assets starting at Jr-1 loses every path at
% its middle age Jr; one starting at Jr-2 loses them all at its last age. Both must give RestrictedSampleMass 0 and NaN stats, with no
% error. Run (i): both pages, so the unrestricted page keeps the list alive while the young page empties. Run (ii): the unrestricted page
% off (whichcombos), so the whole list empties mid-window (the merge of an empty list used to error) and is then extended while empty.
FnsC.assets=FnsToEvaluate.assets;
simoptions_C=simoptions;
simoptions_C.conditionalrestrictions=struct();
simoptions_C.conditionalrestrictions.young=@(aprime,a,z,agej,Jr) (agej<Jr);
simoptions_C.ntimeperiods=3;
simoptions_C.pi_z=pi_z_J;
simoptions_C.whichstats=[1,1,0,0,0,0,0];
simoptions_C0=simoptions_C;
simoptions_C0=rmfield(simoptions_C0,'conditionalrestrictions'); % the same run without the restriction
nrun=0; % the runs that completed
try
    ACC1=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsC,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_C);
    simoptions_C2=simoptions_C;
    simoptions_C2.whichcombos=cat(3,zeros(1,N_j),ones(1,N_j)); % [nFns, N_j, 1+nRestr]: only the young page
    ACC2=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsC,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_C2);
    ACC0=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsC,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_C0);
    nrun=3;
catch ME
    fprintf('   (the error was: %s) \n',ME.message);
end
fprintf('ACSWindows (c) restricted window that every path leaves: runs errored, this should be zero: %i \n',3-nrun)
if nrun==3
    Y1=ACC1.young.assets_3periodAverage;
    Y2=ACC2.young.assets_3periodAverage;
    U1=ACC1.assets_3periodAverage;
    U0=ACC0.assets_3periodAverage;
    StationaryDist=reshape(StationaryDist,[numel(StationaryDist)/N_j,N_j]); % for the start-age masses below
    dead=[Jr-2,Jr-1]; % start ages whose window every path leaves (mid-window for Jr-1, at the last age for Jr-2)
    for rr=1:2
        if rr==1
            Y=Y1; runstr='(i) both pages';
        else
            Y=Y2; runstr='(ii) young page only';
        end
        m=gather(Y.RestrictedSampleMass(dead)); mn=gather(Y.Mean(dead)); md=gather(Y.Median(dead));
        fprintf('ACSWindows (c) %s, start ages Jr-2, Jr-1: RestrictedSampleMass, this should be zero: %.3e; Mean and Median not NaN, this should be zero: %i \n',runstr,max(abs(m)),sum(~isnan([mn,md])))
        % non-vacuous: at Jr-3 the window ends at Jr-1, so every path stays young: the young stats are the unrestricted ones, and the mass is the start-age mass
        j0=Jr-3;
        fprintf('ACSWindows (c) %s, start age Jr-3 (all paths stay young): Mean vs unrestricted, this should be zero: %.3e; Median vs unrestricted, this should be zero: %.3e; RestrictedSampleMass vs the start-age mass, this should be zero: %.3e \n',runstr,abs(gather(Y.Mean(j0)-U0.Mean(j0))),abs(gather(Y.Median(j0)-U0.Median(j0))),abs(gather(Y.RestrictedSampleMass(j0))-gather(sum(StationaryDist(:,j0)))))
        fprintf('ACSWindows (c) %s, start age Jr-3: RestrictedSampleMass is positive (else the checks above are vacuous), this should NOT be zero: %.3e \n',runstr,gather(Y.RestrictedSampleMass(j0)))
    end
    % the two runs agree on the young page, and the unrestricted page of run (i) is unaffected by the dead young page in its list
    % (the Means are NaN by design at the dead start ages and past N_j-2, so the NaN pattern is compared on its own and the values where both are finite)
    dY=abs(gather(Y1.Mean(:))-gather(Y2.Mean(:))); dU=abs(gather(U1.Mean(:))-gather(U0.Mean(:))); dUm=abs(gather(U1.Median(:))-gather(U0.Median(:)));
    fprintf('ACSWindows (c) young page, run (i) vs (ii): Mean, this should be zero: %.3e, NaN-pattern mismatches %i; RestrictedSampleMass, this should be zero: %.3e \n',max([0;dY(~isnan(dY))]),sum(isnan(gather(Y1.Mean(:)))~=isnan(gather(Y2.Mean(:)))),CoreMaxAbs(gather(Y1.RestrictedSampleMass(1:N_j-2))-gather(Y2.RestrictedSampleMass(1:N_j-2))))
    fprintf('ACSWindows (c) unrestricted page of run (i) vs the run without the restriction: Mean, this should be zero: %.3e, NaN-pattern mismatches %i; Median, this should be zero: %.3e, NaN-pattern mismatches %i \n',max([0;dU(~isnan(dU))]),sum(isnan(gather(U1.Mean(:)))~=isnan(gather(U0.Mean(:)))),max([0;dUm(~isnan(dUm))]),sum(isnan(gather(U1.Median(:)))~=isnan(gather(U0.Median(:)))))
end

output.selectors=selectors;
end
