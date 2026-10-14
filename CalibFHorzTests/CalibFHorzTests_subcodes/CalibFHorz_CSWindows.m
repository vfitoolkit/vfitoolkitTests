function output=CalibFHorz_CSWindows(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% CalibFHorzTests part 11: CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets of the windows of
% simoptions.ntimeperiods (added 2026-10-14): fnname_MperiodAverage (or _MperiodSum) with Mean, StdDeviation,
% CovarianceWith.(fn2_MperiodAverage), CorrelationWith.(fn2_MperiodAverage), and the matrices CovarianceMatrix_MperiodAverage/
% CorrelationMatrix_MperiodAverage, under a restriction too; an age-conditional window target is by start age (1 x N_j).
% SetupTargetMoments_FHorz reads them off the target names into selectors.CrossSec/.AgeCrossSec (whichcombos block diagonal in blocks of
% ntimeperiods, ntimeperiodstype; agegroupings and split for the age-conditional one) and, when the age-conditional single-period
% targets are by age groups, into a second call selectors.AgeCrossSecW; the objective makes the call(s).
%   (a) parsing, no model solved: the selectors against hand-built ones (both commands); the split with agegroupings=1:5:N_j; windows
%       only with agegroupings=1:5:N_j (forced to 1:N_j); six invalid target sets error (MperiodAverage and MperiodSum together; a pair
%       of windows of different lengths; a pair of a window and a single-period function; a 1-period window; a window longer than N_j;
%       an age-conditional window target by age groups instead of start ages)
%   (b) one solve at the truth, full runs of the two commands with the windows; targets from them with NaN omissions; the objective's
%       moment vector (vectoroutput=1) under caliboptions.whichcombos=1 and =0 equals the hand walk of the targets through the full
%       runs and the targets themselves; with the default agegroupings (one call each) and with agegroupings=1:5:N_j (the
%       age-conditional windows in a second call)
% Each objective run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.inc=@(aprime,a,z,w,kappa_j) 1+a+w*kappa_j*z;
Jr=Params.Jr;

simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);

%% (a) parsing, no model solved
TM=struct();
TM.CrossSectionCovarCorr.assets.StdDeviation=1;                                                    % M=1, page 1
TM.CrossSectionCovarCorr.earnings_2periodAverage.CorrelationWith.inc_2periodAverage=1;             % M=2, page 1
TM.CrossSectionCovarCorr.rich.assets_3periodAverage.Mean=1;                                        % M=3, rich
t=ones(1,N_j); t(Jr:end)=NaN; TM.AgeConditionalCrossSectionCovarCorr.earnings.Mean=t;             % M=1, ages 1..Jr-1
t=nan(1,N_j); t(2:10)=1; TM.AgeConditionalCrossSectionCovarCorr.inc_2periodAverage.CovarianceWith.assets_2periodAverage=t; % M=2, start ages 2..10
t=nan(3,3,N_j); t(2,3,1:Jr-2)=1; TM.AgeConditionalCrossSectionCovarCorr.rich.CorrelationMatrix_2periodAverage=t;            % M=2, rich, the (earnings,inc) entry at start ages 1..Jr-2
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors,~,~,~,~,~,~,~,~,FnsX,~,~,~,FnsXA]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,0,N_j,simoptions);
% hand-built: the functions assets, earnings, inc; CrossSec in blocks of [1,2,3] (9 x 9 x 2 pages), AgeCrossSec of [1,2] (6 x 6 x N_j x 2)
wc=zeros(9,9,2);
wc(1,1,1)=1;                 % assets, M=1
wc(5,6,1)=1; wc(6,5,1)=1;    % (earnings, inc), M=2
wc(7,7,2)=1;                 % rich, assets, M=3
wca=zeros(6,6,N_j,2);
wca(2,2,1:Jr-1,1)=1;                          % earnings, M=1
wca(6,4,2:10,1)=1; wca(4,6,2:10,1)=1;         % (inc, assets), M=2
wca(5,6,1:Jr-2,2)=1; wca(6,5,1:Jr-2,2)=1;     % rich, (earnings, inc), M=2
s=selectors.CrossSec; sa=selectors.AgeCrossSec;
bad=~isequal(fieldnames(FnsX),{'assets';'earnings';'inc'})+~isequal(fieldnames(FnsXA),{'assets';'earnings';'inc'});
bad=bad+~(isfield(s,'ntimeperiods') && isequal(s.ntimeperiods,[1,2,3]) && strcmp(s.ntimeperiodstype,'average'));
bad=bad+~(isfield(sa,'ntimeperiods') && isequal(sa.ntimeperiods,[1,2]) && strcmp(sa.ntimeperiodstype,'average') && isempty(sa.agegroupings) && sa.split==0 && ~isfield(selectors,'AgeCrossSecW'));
if bad==0 && isequal(size(s.whichcombos),size(wc)) && isequal(size(sa.whichcombos),size(wca))
    d=max([abs(s.whichcombos(:)-wc(:));abs(sa.whichcombos(:)-wca(:))]);
else
    d=1;
end
fprintf('CSWindows (a) reduced FnsToEvaluate, ntimeperiods [1 2 3] and [1 2], average, no agegroupings, no split: mismatches, this should be zero: %i; selectors whichcombos vs hand (both commands), this should be zero: %.3e \n',bad,d)
% the split: the single-period age-conditional target by the age groups 1:5:N_j (four groups), the windows by start age in a second call
simoptions_g=simoptions; simoptions_g.agegroupings=1:5:N_j;
TMg=TM;
TMg.AgeConditionalCrossSectionCovarCorr.earnings.Mean=[1,1,1,NaN];
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,selg]=SetupTargetMoments_FHorz(TMg,FnsToEvaluate,0,N_j,simoptions_g);
bad=0;
if isfield(selg.AgeCrossSec,'split') && selg.AgeCrossSec.split==1 && isfield(selg,'AgeCrossSecW')
    sA=selg.AgeCrossSec; sW=selg.AgeCrossSecW;
    bad=bad+~isequal(sA.ntimeperiods,1)+~isempty(sA.agegroupings)+~isequal(size(sA.whichcombos),[3,3,4,2])+~isequal(sW.ntimeperiods,2)+~isequal(sW.agegroupings,1:N_j)+~isequal(size(sW.whichcombos),[3,3,N_j,2]);
    if bad==0
        bad=bad+~isequal(reshape(sA.whichcombos(2,2,:,1),1,[]),[1,1,1,0])+~isequal(sW.whichcombos,wca(4:6,4:6,:,:));
    end
    bad=bad+~isequal(selg.CrossSec.ntimeperiods,[1,2,3]); % (the pooled command has no age groups)
else
    bad=1;
end
fprintf('CSWindows (a) agegroupings 1:5:N_j with single-period and window age-conditional targets: the split into AgeCrossSec (age groups) and AgeCrossSecW (start ages), mismatches vs hand, this should be zero: %i \n',bad)
% windows only, with agegroupings 1:5:N_j: one call, by start age
TMw=struct();
TMw.AgeConditionalCrossSectionCovarCorr.inc_2periodAverage=TM.AgeConditionalCrossSectionCovarCorr.inc_2periodAverage;
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,selw]=SetupTargetMoments_FHorz(TMw,FnsToEvaluate,0,N_j,simoptions_g);
bad=~(isfield(selw.AgeCrossSec,'split') && selw.AgeCrossSec.split==0 && isequal(selw.AgeCrossSec.agegroupings,1:N_j) && isequal(selw.AgeCrossSec.ntimeperiods,2) && isequal(size(selw.AgeCrossSec.whichcombos),[2,2,N_j,2]));
fprintf('CSWindows (a) age-conditional windows only with agegroupings 1:5:N_j: one call with agegroupings forced to 1:N_j, mismatches, this should be zero: %i \n',bad)
% invalid target sets
errcases={'MperiodAverage and MperiodSum together','a pair of windows of different lengths','a pair of a window and a single-period function','a 1-period window','a window longer than N_j','an age-conditional window target by age groups'};
errmsgs={'both MperiodAverage and MperiodSum','of the same window as the first','of the same window as the first','at least 2 periods','cannot be longer than N_j','start ages'};
for ee=1:6
    TMbad=TM;
    so=simoptions;
    switch ee
        case 1
            TMbad.CrossSectionCovarCorr.earnings_2periodSum.Mean=1;
        case 2
            TMbad.CrossSectionCovarCorr.earnings_2periodAverage.CovarianceWith.inc_3periodAverage=1;
        case 3
            TMbad.AgeConditionalCrossSectionCovarCorr.earnings_2periodAverage.CorrelationWith.inc=ones(1,N_j);
        case 4
            TMbad.CrossSectionCovarCorr.earnings_1periodAverage.Mean=1;
        case 5
            TMbad.CrossSectionCovarCorr.(['earnings_',num2str(N_j+1),'periodAverage']).Mean=1;
        case 6
            TMbad=TMg;
            TMbad.AgeConditionalCrossSectionCovarCorr.inc_2periodAverage.CovarianceWith.assets_2periodAverage=ones(1,4);
            so=simoptions_g;
    end
    errok=0;
    try
        SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,so);
    catch ME
        errok=contains(ME.message,errmsgs{ee});
        if ~errok
            fprintf('   (the error was: %s) \n',ME.message);
        end
    end
    fprintf('CSWindows (a) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok);
end

%% (b) one solve at the truth, full runs of the two commands, the objective's moment vector
vfoptions=struct();
[z_gridvals_J,pi_z_J,vfoptions]=ExogShockSetup_FHorz(n_z,z_grid,pi_z,N_j,Params,vfoptions,3,0);
simoptions.e_gridvals_J=vfoptions.e_gridvals_J;
simoptions.pi_e_J=vfoptions.pi_e_J;
vfoptions.alreadygridvals=1;
simoptions.alreadygridvals=1;
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
simoptions_L=simoptions; simoptions_L.ntimeperiods=[1,2,3]; simoptions_L.pi_z=pi_z_J;
CSF=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_L); % pooled, every window
simoptions_L.ntimeperiods=[1,2];
ACF=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_L); % by age, every window
simoptions_A=simoptions; simoptions_A.agegroupings=1:5:N_j;
ACA=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_A); % the single period by age groups
kinds={'CrossSectionCovarCorr','AgeConditionalCrossSectionCovarCorr'};
for scen=1:2
    ACFull=ACF;
    simoptions_s=simoptions;
    if scen==1
        scenstr='one call each (default agegroupings)';
    else
        scenstr='two age-conditional calls (agegroupings 1:5:N_j)';
        ACFull.earnings=ACA.earnings; % the single-period age-conditional fields by age groups
        ACFull.assets=ACA.assets;
        simoptions_s.agegroupings=1:5:N_j;
    end
    % the targets: the structure of (a) (plus a sum-free matrix target of the pooled command), the values from the full runs, NaN where omitted
    TMb=struct();
    TMb.CrossSectionCovarCorr.assets.StdDeviation=gather(CSF.assets.StdDeviation);
    TMb.CrossSectionCovarCorr.earnings_2periodAverage.CorrelationWith.inc_2periodAverage=gather(CSF.earnings_2periodAverage.CorrelationWith.inc_2periodAverage);
    TMb.CrossSectionCovarCorr.rich.assets_3periodAverage.Mean=gather(CSF.rich.assets_3periodAverage.Mean);
    t=gather(CSF.CovarianceMatrix_3periodAverage); t(1,:)=NaN; TMb.CrossSectionCovarCorr.CovarianceMatrix_3periodAverage=t;
    t=gather(ACFull.earnings.Mean);
    if scen==1
        t(Jr:end)=NaN;
    else
        t(4)=NaN;
    end
    TMb.AgeConditionalCrossSectionCovarCorr.earnings.Mean=t;
    t=gather(ACFull.inc_2periodAverage.CovarianceWith.assets_2periodAverage); t([1,11:end])=NaN; TMb.AgeConditionalCrossSectionCovarCorr.inc_2periodAverage.CovarianceWith.assets_2periodAverage=t;
    t=nan(3,3,N_j); temp=gather(ACFull.rich.CorrelationMatrix_2periodAverage); t(2,3,1:Jr-2)=temp(2,3,1:Jr-2); TMb.AgeConditionalCrossSectionCovarCorr.rich.CorrelationMatrix_2periodAverage=t;
    t=gather(ACFull.rich.assets_2periodAverage.StdDeviation); t(N_j)=NaN; TMb.AgeConditionalCrossSectionCovarCorr.rich.assets_2periodAverage.StdDeviation=t;
    % the hand walk: every leaf of the targets in depth-first fieldnames order (as SetupTargetMoments reads them), the same path in the full runs
    hand=[]; tvec=[];
    for kk=1:2
        if kk==1
            Full=CSF;
        else
            Full=ACFull;
        end
        paths={{}};
        while ~isempty(paths)
            p=paths{1}; paths(1)=[];
            node=TMb.(kinds{kk}); v=Full;
            for ll=1:length(p)
                node=node.(p{ll}); v=v.(p{ll});
            end
            if isstruct(node)
                fn=fieldnames(node);
                children=cell(1,length(fn));
                for ll=1:length(fn)
                    children{ll}=[p,fn(ll)];
                end
                paths=[children,paths]; %#ok<AGROW>
            else
                v=gather(v); hand=[hand; v(:)]; tvec=[tvec; node(:)]; %#ok<AGROW>
            end
        end
    end
    actual=~isnan(tvec);
    fprintf('CSWindows (b) %s: %i target entries, %i of them targeted; the model is NaN at a targeted entry, this should be zero: %i \n',scenstr,length(tvec),sum(actual),sum(isnan(hand(actual))))

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
    fprintf('CSWindows (b) %s: SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %i; split as expected, this should be zero: %i \n',scenstr,CoreMaxAbs(d),patternbad,~(selectors.AgeCrossSec.split==(scen==2)))
    handA=hand(actual); tA=targetmomentvec(actual);
    for wcc=[1,0]
        caliboptions_t=caliboptions; caliboptions_t.whichcombos=wcc;
        ran=0; m=[];
        try
            m=gather(CalibrateLifeCycleModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions_s));
            ran=1;
        catch err
            fprintf('CSWindows (b) %s, whichcombos=%i: the objective errored with: %s \n',scenstr,wcc,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        end
        if ran==1 && length(m)==length(handA)
            d=CoreMaxAbs(m(:)-handA(:)); dt=CoreMaxAbs(m(:)-tA(:)); nanbad=sum(isnan(m));
        else
            d=1; dt=1; nanbad=1;
        end
        fprintf('CSWindows (b) %s, whichcombos=%i: did not run, this should be zero: %i; moment vector vs hand-built, this should be zero: %.3e; vs the targets, this should be zero: %.3e; NaN entries or length mismatch, this should be zero: %i \n',scenstr,wcc,1-ran,d,dt,nanbad)
    end
end

output.selectors=selectors;
end
