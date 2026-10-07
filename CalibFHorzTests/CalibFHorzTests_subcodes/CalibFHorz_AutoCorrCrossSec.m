function output=CalibFHorz_AutoCorrCrossSec(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% CalibFHorzTests part 8: the three target kinds added 2026-10-07, TargetMoments.AutoCorrTransProbs (EvalFnOnAgentDist_AutoCorrTransProbs_FHorz),
% .CrossSectionCovarCorr (EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz) and .AgeConditionalCrossSectionCovarCorr (the age-conditional
% version), standing alongside AllStats and AgeConditionalStats in the target vector (in that order, then CustomModelStats).
%   (a) parsing, no model solved: the target vector against a hand-built one; the flags, the horizons read off the _kK suffixes
%       (autocorrtimehorizons), the cumulative sizes; the selectors against hand-built ones (AutoCorr [nFns,N_j,1+nRestr] by start age,
%       CrossSec [nFns,nFns] pair-shaped and symmetric, AgeCrossSec [nFns,nFns,ngroups]; a matrix target makes the reduced set the full
%       one); the same targets in another field order give the same vector and selectors (the kinds have a fixed order); a call without
%       the new fields gives empty new outputs; the three-input legacy call with the new fields errors; ten invalid inputs error
%   (b) one solve at the truth, full runs of the three commands; targets from them with NaN omissions (start ages, pairs, age groups,
%       matrix entries) and simoptions.agegroupings for the age-conditional kinds (the AutoCorr targets are per age regardless) and a
%       user simoptions.timehorizons that the targets' horizons are joined to; the objective's moment vector (vectoroutput=1) under
%       caliboptions.whichcombos=1 and =0 equals the hand-built one and the targets (the targets are the model)
%   (c) the selectors handed straight to the three commands: untargeted start ages NaN, a function targeted only under a restriction
%       has no unrestricted output, an unwanted pair NaN, an unwanted age group NaN, targeted entries equal to the full run
%   (d) round trips through CalibrateLifeCycleModel (grid interpolation layer): rho_z (ExogShockFn) from the earnings autocorrelation
%       profile, at the truth (exact) and from 0.7 (close); beta from the assets-earnings correlation and the age profile of the
%       assets-consumption covariance, at the truth (exact) and from 0.90 (close), whichcombos=0 identical to =1 (exact);
%       caliboptions.logmoments by name on an AutoCorrTransProbs target (the driver's by-name walk for the new kinds)
%   (e) runtime of one objective evaluation with a sparse AutoCorr target, whichcombos=1 vs =0 (informational)
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) (1+r)*a+w*kappa_j*z*(agej<Jr)+pension*(agej>=Jr)-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);
FnNames=fieldnames(FnsToEvaluate);
Jr=Params.Jr;

simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);
simoptions.agegroupings=1:5:N_j;
ngroups=length(simoptions.agegroupings);
npages=3; % unrestricted, rich, poor

%% (a) parsing, no model solved
TM=struct();
TM.AllStats.assets.Mean=1.1;
TM.AgeConditionalStats.assets.Mean=[NaN,2,3,4];
t=ones(1,N_j); t(Jr:end)=NaN; TM.AutoCorrTransProbs.earnings.Mean=t;                    % ages 1..Jr-1
t=ones(1,N_j-1); t(Jr-1:end)=NaN; TM.AutoCorrTransProbs.earnings.AutoCorrelation=t;    % start ages 1..Jr-2, horizon 1
t=nan(1,N_j-3); t(2:4)=1; TM.AutoCorrTransProbs.earnings.AutoCovariance_k3=t;          % start ages 2..4, horizon 3
t=ones(1,N_j-1); t(1:5)=NaN; TM.AutoCorrTransProbs.rich.assets.AutoCorrelation=t;      % rich, start ages 6..N_j-1
TM.AutoCorrTransProbs.rich.assets.PairMass_k2=ones(1,N_j-2);                           % rich, start ages 1..N_j-2, horizon 2
TM.AutoCorrTransProbs.poor.consumption.StdDeviation=ones(1,N_j);                       % poor, every age
TM.CrossSectionCovarCorr.assets.CorrelationWith.earnings=0.3;
TM.CrossSectionCovarCorr.consumption.Mean=1;
TM.AgeConditionalCrossSectionCovarCorr.assets.CovarianceWith.consumption=[NaN,1,1,1];
t=nan(4,4,ngroups); t(1,2,2:3)=1; t(4,4,1)=1; TM.AgeConditionalCrossSectionCovarCorr.CorrelationMatrix=t; % a matrix: over every FnsToEvaluate

[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,0,N_j,simoptions);

% the hand-built target vector: AllStats, AgeConditionalStats, AutoCorrTransProbs, CrossSectionCovarCorr, AgeConditionalCrossSectionCovarCorr, each in its field order, vectorised column-major
hand=[1.1; [NaN;2;3;4]; TM.AutoCorrTransProbs.earnings.Mean(:); TM.AutoCorrTransProbs.earnings.AutoCorrelation(:); TM.AutoCorrTransProbs.earnings.AutoCovariance_k3(:); ...
    TM.AutoCorrTransProbs.rich.assets.AutoCorrelation(:); TM.AutoCorrTransProbs.rich.assets.PairMass_k2(:); TM.AutoCorrTransProbs.poor.consumption.StdDeviation(:); ...
    0.3; 1; [NaN;1;1;1]; TM.AgeConditionalCrossSectionCovarCorr.CorrelationMatrix(:)];
if length(targetmomentvec)==length(hand)
    d=abs(targetmomentvec-hand); d=d(:); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(hand));
else
    d=1; patternbad=abs(length(targetmomentvec)-length(hand));
end
fprintf('AutoCorrCrossSec (a) targetmomentvec vs hand-built (%i entries), this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',length(hand),max([d;0]),patternbad)
handAC=N_j+(N_j-1)+(N_j-3)+(N_j-1)+(N_j-2)+N_j; handX=2; handAX=ngroups+16*ngroups;
fprintf('AutoCorrCrossSec (a) flags usingautocorr, usingcrosssec, usingagecrosssec = %i %i %i (hand 1 1 1), usingallstats, usinglcp, usingcustomstats = %i %i %i (hand 1 1 0); cumulative sizes end AutoCorr %i (hand %i), CrossSec %i (hand %i), AgeCrossSec %i (hand %i), this should be zero: %.3e \n',usingautocorr,usingcrosssec,usingagecrosssec,usingallstats,usinglcp,usingcustomstats,autocorrcummomentsizes(end),handAC,crossseccummomentsizes(end),handX,agecrossseccummomentsizes(end),handAX,abs(usingautocorr-1)+abs(usingcrosssec-1)+abs(usingagecrosssec-1)+abs(usingallstats-1)+abs(usinglcp-1)+abs(usingcustomstats)+abs(autocorrcummomentsizes(end)-handAC)+abs(crossseccummomentsizes(end)-handX)+abs(agecrossseccummomentsizes(end)-handAX)+abs(allstatcummomentsizes(end)-1)+abs(acscummomentsizes(end)-4)+abs(cmscummomentsizes))
fprintf('AutoCorrCrossSec (a) autocorrtimehorizons [%s] (hand [2 3]), names tables of %i, %i, %i rows (hand 6, 2, 2): horizons, this should be zero: %.3e; rows, this should be zero: %.3e \n',num2str(autocorrtimehorizons),size(autocorrmomentnames,1),size(crosssecmomentnames,1),size(agecrosssecmomentnames,1),~isequal(autocorrtimehorizons(:)',[2,3]),abs(size(autocorrmomentnames,1)-6)+abs(size(crosssecmomentnames,1)-2)+abs(size(agecrosssecmomentnames,1)-2))

% the selectors: AutoCorr over assets, earnings, consumption (FnsToEvaluate order; richind has no AutoCorr target)
FnAC=fieldnames(FnsToEvaluate_AutoCorr);
wcAC=zeros(3,N_j,npages);
wcAC(2,1:Jr-1,1)=1;        % earnings unrestricted: Mean ages 1..Jr-1, AutoCorrelation start ages 1..Jr-2, AutoCovariance_k3 start ages 2..4
wcAC(1,1:N_j-1,2)=1;       % assets under rich: AutoCorrelation start ages 6..N_j-1 and PairMass_k2 start ages 1..N_j-2
wcAC(3,:,3)=1;             % consumption under poor: StdDeviation every age
fnok=isequal(FnAC,{'assets';'earnings';'consumption'});
fprintf('AutoCorrCrossSec (a) AutoCorr reduced FnsToEvaluate are assets, earnings, consumption in that order, this should be zero: %.3e \n',1-fnok)
if fnok
    fprintf('AutoCorrCrossSec (a) AutoCorr whichcombos [%s] vs hand, this should be zero: %.3e; size, this should be zero: %.3e \n',num2str(size(selectors.AutoCorr.whichcombos)),max(abs(selectors.AutoCorr.whichcombos(:)-wcAC(:))),~isequal(size(selectors.AutoCorr.whichcombos),[3,N_j,npages]))
end
% CrossSec over assets, earnings, consumption (no matrix target, richind untargeted); AgeCrossSec over all four (the matrix target)
FnX=fieldnames(FnsToEvaluate_CrossSec);
wcX=zeros(3,3); wcX(1,2)=1; wcX(2,1)=1; wcX(3,3)=1; % (earnings is in the set through its pair with assets; its own stats are byproducts, not targets)
FnAX=fieldnames(FnsToEvaluate_AgeCrossSec);
wcAX=zeros(4,4,ngroups); wcAX(1,3,2:4)=1; wcAX(3,1,2:4)=1; wcAX(1,2,2:3)=1; wcAX(2,1,2:3)=1; wcAX(4,4,1)=1;
fnok=isequal(FnX,{'assets';'earnings';'consumption'}); fnok2=isequal(FnAX,FnNames);
fprintf('AutoCorrCrossSec (a) CrossSec reduced FnsToEvaluate are assets, earnings, consumption, this should be zero: %.3e; AgeCrossSec reduced FnsToEvaluate are all four (a matrix target), this should be zero: %.3e \n',1-fnok,1-fnok2)
if fnok
    fprintf('AutoCorrCrossSec (a) CrossSec whichcombos [%s] vs hand, this should be zero: %.3e; size, this should be zero: %.3e \n',num2str(size(selectors.CrossSec.whichcombos)),max(abs(selectors.CrossSec.whichcombos(:)-wcX(:))),~isequal(size(selectors.CrossSec.whichcombos),[3,3]))
end
if fnok2
    fprintf('AutoCorrCrossSec (a) AgeCrossSec whichcombos [%s] vs hand, this should be zero: %.3e; size, this should be zero: %.3e \n',num2str(size(selectors.AgeCrossSec.whichcombos)),max(abs(selectors.AgeCrossSec.whichcombos(:)-wcAX(:))),~isequal(size(selectors.AgeCrossSec.whichcombos,1:3),[4,4,ngroups]))
end

% another field order: the new kinds first, the restricted AutoCorr targets first; the kinds have a fixed order in the vector, so it must be identical
TM2=struct();
TM2.AgeConditionalCrossSectionCovarCorr=TM.AgeConditionalCrossSectionCovarCorr;
TM2.CrossSectionCovarCorr=TM.CrossSectionCovarCorr;
TM2.AutoCorrTransProbs.poor=TM.AutoCorrTransProbs.poor;
TM2.AutoCorrTransProbs.rich=TM.AutoCorrTransProbs.rich;
TM2.AutoCorrTransProbs.earnings=TM.AutoCorrTransProbs.earnings;
TM2.AgeConditionalStats=TM.AgeConditionalStats;
TM2.AllStats=TM.AllStats;
[tv2,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors2,~,~,acsz2,FnAC2,ath2,~,~,xsz2,FnX2,~,~,axsz2,FnAX2]=SetupTargetMoments_FHorz(TM2,FnsToEvaluate,0,N_j,simoptions);
bad=~isequal(size(tv2),size(targetmomentvec))+~isequal(fieldnames(FnAC2),FnAC)+~isequal(fieldnames(FnX2),FnX)+~isequal(fieldnames(FnAX2),FnAX)+~isequal(ath2(:),autocorrtimehorizons(:))+~isequal(acsz2(end),autocorrcummomentsizes(end))+~isequal(xsz2(end),crossseccummomentsizes(end))+~isequal(axsz2(end),agecrossseccummomentsizes(end));
maxd=0;
if isequal(size(tv2),size(targetmomentvec))
    s1=sort(tv2(~isnan(tv2))); s2=sort(targetmomentvec(~isnan(targetmomentvec)));
    if length(s1)==length(s2)
        maxd=max([maxd;abs(s1-s2)]);
    else
        bad=bad+1;
    end
    % the AutoCorr block is a permutation (the restricted targets come first within the kind); the other kinds are in the same order
    maxd=max([maxd;abs(tv2(1:5)-targetmomentvec(1:5))]); % AllStats and AgeConditionalStats (the NaN entry compares NaN, dropped below)
end
maxd=max([maxd;abs(selectors2.AutoCorr.whichcombos(:)-selectors.AutoCorr.whichcombos(:));abs(selectors2.CrossSec.whichcombos(:)-selectors.CrossSec.whichcombos(:));abs(selectors2.AgeCrossSec.whichcombos(:)-selectors.AgeCrossSec.whichcombos(:))]);
maxd=maxd(~isnan(maxd));
fprintf('AutoCorrCrossSec (a) another field order (new kinds first): selectors, horizons and sorted target vector vs the first order, this should be zero: %.3e; sizes or function sets differ, this should be zero: %.3e \n',max([maxd;0]),bad)

% without the new fields: empty new outputs (and the AllStats/AgeConditionalStats outputs as before)
TM0=struct(); TM0.AllStats=TM.AllStats; TM0.AgeConditionalStats=TM.AgeConditionalStats;
[tv0,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors0,ua0,acn0,acs0,FnAC0,ath0,ux0,xn0,xs0,FnX0,uax0,axn0,axs0,FnAX0]=SetupTargetMoments_FHorz(TM0,FnsToEvaluate,0,N_j,simoptions);
bad=abs(ua0)+abs(ux0)+abs(uax0)+~isequal(acs0,0)+~isequal(xs0,0)+~isequal(axs0,0)+~isempty(fieldnames(FnAC0))+~isempty(fieldnames(FnX0))+~isempty(fieldnames(FnAX0))+~isempty(ath0);
bad=bad+~(iscell(acn0) && size(acn0,2)==3 && isempty(acn0{1,1}))+~(iscell(xn0) && isempty(xn0{1,1}))+~(iscell(axn0) && isempty(axn0{1,1}));
bad=bad+~isempty(selectors0.AutoCorr.whichcombos)+~isempty(selectors0.CrossSec.whichcombos)+~isempty(selectors0.AgeCrossSec.whichcombos);
d=abs(tv0-targetmomentvec(1:5)); d=d(~isnan(d));
fprintf('AutoCorrCrossSec (a) targets without the new fields: new outputs empty and flags zero, this should be zero: %.3e; target vector vs the first five entries, this should be zero: %.3e \n',bad,max([d;0]))

% invalid inputs must error: the legacy call with the new fields, and ten bad targets
errored=0;
try
    SetupTargetMoments_FHorz(TM,FnsToEvaluate,0); % three-input call (no N_j, no simoptions) with the new fields
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AutoCorrTransProbs.earnings.Skewness=ones(1,N_j); % unknown statistic
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AutoCorrTransProbs.earnings.AutoCorrelation=ones(1,N_j); % wrong length (1 x N_j-1)
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AutoCorrTransProbs.earnings.TransitionProbs=1; % TransitionProbs cannot be targeted
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AutoCorrTransProbs.wealth.Mean=ones(1,N_j); % unknown name
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.SkewnessMatrix=ones(4); % a one-level target that is not CovarianceMatrix/CorrelationMatrix
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.assets.CovarianceWith=1; % no second function
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.CovarianceMatrix=ones(3,3); % a matrix target is over every FnsToEvaluate (four)
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.rich.assets.Mean=1; % the cross-section commands have no conditional restrictions
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AgeConditionalCrossSectionCovarCorr.assets.CovarianceWith.earnings=ones(1,N_j); % N_j entries but ngroups age groups
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.earnings.earnings=1; % the self-correlation is always one, not a target (dropped 2026-10-08)
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
fprintf('AutoCorrCrossSec (a) the legacy call with the new fields and ten invalid targets must all error, this should be zero: %.3e \n',11-errored)

%% (b) one solve at the truth, full runs of the three commands, the objective's moment vector
vfoptions=struct();
[z_gridvals_J,pi_z_J,vfoptions]=ExogShockSetup_FHorz(n_z,z_grid,pi_z,N_j,Params,vfoptions,3,0);
simoptions.e_gridvals_J=vfoptions.e_gridvals_J;
simoptions.pi_e_J=vfoptions.pi_e_J;
vfoptions.alreadygridvals=1;
simoptions.alreadygridvals=1;
simoptions.timehorizons=2; % the user's horizon; the targets name horizon 3 as well, the objective must do both
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_full);
simoptions_fullAC=rmfield(simoptions,'agegroupings'); simoptions_fullAC.timehorizons=[2,3]; % (the command has no age bins)
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,simoptions_fullAC);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions);
AXFull=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions);

% the targets: the same structure as in (a), values from the full runs, NaN where omitted
TMb=struct();
TMb.AllStats.assets.Mean=gather(AllFull.assets.Mean);
t=gather(LCPFull.assets.Mean); t(1)=NaN; TMb.AgeConditionalStats.assets.Mean=t;
t=gather(ACFull.earnings.Mean); t(Jr:end)=NaN; TMb.AutoCorrTransProbs.earnings.Mean=t;
t=gather(ACFull.earnings.AutoCorrelation); t(Jr-1:end)=NaN; TMb.AutoCorrTransProbs.earnings.AutoCorrelation=t;
t=gather(ACFull.earnings.AutoCovariance_k3); t([1,5:end])=NaN; TMb.AutoCorrTransProbs.earnings.AutoCovariance_k3=t;
t=gather(ACFull.rich.assets.AutoCorrelation); t(1:5)=NaN; TMb.AutoCorrTransProbs.rich.assets.AutoCorrelation=t;
TMb.AutoCorrTransProbs.rich.assets.PairMass_k2=gather(ACFull.rich.assets.PairMass_k2);
TMb.AutoCorrTransProbs.poor.consumption.StdDeviation=gather(ACFull.poor.consumption.StdDeviation);
TMb.CrossSectionCovarCorr.assets.CorrelationWith.earnings=gather(CSFull.assets.CorrelationWith.earnings);
TMb.CrossSectionCovarCorr.consumption.Mean=gather(CSFull.consumption.Mean);
t=gather(AXFull.assets.CovarianceWith.consumption); t(1)=NaN; TMb.AgeConditionalCrossSectionCovarCorr.assets.CovarianceWith.consumption=t;
t=gather(AXFull.CorrelationMatrix); m=nan(size(t)); m(1,2,2:3)=t(1,2,2:3); m(4,4,1)=t(4,4,1); TMb.AgeConditionalCrossSectionCovarCorr.CorrelationMatrix=m;

% the hand-built vector: walk TargetMoments kind by kind (the order of the vector) and read the same field of the full run (same nesting)
kinds={'AllStats','AgeConditionalStats','AutoCorrTransProbs','CrossSectionCovarCorr','AgeConditionalCrossSectionCovarCorr'};
fulls={AllFull,LCPFull,ACFull,CSFull,AXFull};
hand=[]; tvec=[];
for kk=1:5
    TMk=TMb.(kinds{kk}); Fk=fulls{kk};
    a1vec=fieldnames(TMk);
    for a1=1:length(a1vec)
        t1=TMk.(a1vec{a1});
        if ~isstruct(t1) % one name: a matrix
            v=gather(Fk.(a1vec{a1})); hand=[hand; v(:)]; tvec=[tvec; t1(:)];
        else
            a2vec=fieldnames(t1);
            for a2=1:length(a2vec)
                t2=t1.(a2vec{a2});
                if ~isstruct(t2)
                    v=gather(Fk.(a1vec{a1}).(a2vec{a2})); hand=[hand; v(:)]; tvec=[tvec; t2(:)];
                else
                    a3vec=fieldnames(t2);
                    for a3=1:length(a3vec)
                        v=gather(Fk.(a1vec{a1}).(a2vec{a2}).(a3vec{a3})); hand=[hand; v(:)]; t3=t2.(a3vec{a3}); tvec=[tvec; t3(:)];
                    end
                end
            end
        end
    end
end
actual=~isnan(tvec);
fprintf('AutoCorrCrossSec (b) %i target entries, %i of them targeted (the rest NaN); the model is NaN at a targeted entry, this should be zero: %.3e \n',length(tvec),sum(actual),sum(isnan(hand(actual))))

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
fprintf('AutoCorrCrossSec (b) SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',max([d;0]),patternbad)
handA=hand(actual); tA=targetmomentvec(actual);
for wc=[1,0]
    caliboptions_t=caliboptions; caliboptions_t.whichcombos=wc;
    ran=0; m=[]; 
    try
        m=gather(CalibrateLifeCycleModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions));
        ran=1;
    catch err
        fprintf('AutoCorrCrossSec (b) whichcombos=%i: the objective errored with: %s \n',wc,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ran==1 && length(m)==length(handA)
        d=max([abs(m(:)-handA(:));0]); dt=max([abs(m(:)-tA(:));0]); nanbad=sum(isnan(m));
    else
        d=1; dt=1; nanbad=1;
    end
    fprintf('AutoCorrCrossSec (b) whichcombos=%i: did not run, this should be zero: %.3e; moment vector vs hand-built, this should be zero: %.3e; vs the targets (targets are the model), this should be zero: %.3e; NaN entries or length mismatch, this should be zero: %.3e \n',wc,1-ran,d,dt,nanbad)
end

%% (c) the selectors handed straight to the three commands
simsel=simoptions_fullAC; simsel.whichcombos=selectors.AutoCorr.whichcombos;
ACsel=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate_AutoCorr,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,simsel);
bad=0; d=[];
bad=bad+sum(~isnan(gather(ACsel.earnings.AutoCorrelation(Jr:end))));            % unselected start ages, horizon 1 (start ages 1..Jr-1 are selected by the Mean target; at Jr-1 the pair is computed, whatever its value)
bad=bad+sum(~isnan(gather(ACsel.earnings.AutoCovariance_k3(Jr:end))));           % unselected start ages, horizon 3 (every horizon follows the start ages)
bad=bad+isfield(ACsel,'assets')+isfield(ACsel,'consumption')+isfield(ACsel,'richind'); % assets only under rich, consumption only under poor, richind nowhere
bad=bad+isfield(ACsel.rich,'earnings')+isfield(ACsel.rich,'consumption')+isfield(ACsel.poor,'assets')+isfield(ACsel.poor,'earnings');
bad=bad+sum(~isnan(gather(ACsel.rich.assets.AutoCorrelation(1:5))));             % unselected start ages under rich
d=[d; abs(gather(ACsel.earnings.Mean(1:Jr-1))-gather(ACFull.earnings.Mean(1:Jr-1)))'];
d=[d; abs(gather(ACsel.earnings.AutoCorrelation(1:Jr-2))-gather(ACFull.earnings.AutoCorrelation(1:Jr-2)))'];
d=[d; abs(gather(ACsel.earnings.AutoCovariance_k3(2:4))-gather(ACFull.earnings.AutoCovariance_k3(2:4)))'];
d=[d; abs(gather(ACsel.rich.assets.AutoCorrelation(6:end))-gather(ACFull.rich.assets.AutoCorrelation(6:end)))'];
d=[d; abs(gather(ACsel.rich.assets.PairMass_k2)-gather(ACFull.rich.assets.PairMass_k2))'];
d=[d; abs(gather(ACsel.poor.consumption.StdDeviation)-gather(ACFull.poor.consumption.StdDeviation))'];
d=d(~isnan(d));
fprintf('AutoCorrCrossSec (c) AutoCorrTransProbs with the selectors: targeted entries vs full run, this should be zero: %.3e; untargeted entries not NaN or present, this should be zero: %.3e; horizon fields k2 and k3 both present, this should be zero: %.3e \n',max([d;0]),bad,~isfield(ACsel.earnings,'AutoCovariance_k2')+~isfield(ACsel.earnings,'AutoCovariance_k3'))
simselX=simoptions; simselX.whichcombos=selectors.CrossSec.whichcombos;
CSsel=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate_CrossSec,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simselX);
bad=0;
bad=bad+~isnan(gather(CSsel.assets.CorrelationWith.consumption));                 % an unwanted pair of two evaluated functions
bad=bad+~isnan(gather(CSsel.CovarianceMatrix(1,3)))+~isnan(gather(CSsel.CovarianceMatrix(2,3)));
bad=bad+isfield(CSsel,'richind');                                                 % not in the reduced set at all
d=[abs(gather(CSsel.assets.CorrelationWith.earnings)-gather(CSFull.assets.CorrelationWith.earnings)); abs(gather(CSsel.consumption.Mean)-gather(CSFull.consumption.Mean)); abs(gather(CSsel.CorrelationMatrix(1,2))-gather(CSFull.CorrelationMatrix(1,2)))];
fprintf('AutoCorrCrossSec (c) CrossSectionCovarCorr with the selectors: targeted entries vs full run, this should be zero: %.3e; unwanted pairs not NaN or function present, this should be zero: %.3e \n',max(d),bad)
simselAX=simoptions; simselAX.whichcombos=selectors.AgeCrossSec.whichcombos;
AXsel=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate_AgeCrossSec,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simselAX);
bad=0;
bad=bad+~isnan(gather(AXsel.assets.CovarianceWith.consumption(1)));               % the omitted age group
bad=bad+~isnan(gather(AXsel.CorrelationMatrix(1,2,1)))+~isnan(gather(AXsel.CorrelationMatrix(1,2,4)));
bad=bad+sum(~isnan(gather(AXsel.CorrelationMatrix(2,3,:))));                      % a pair wanted in no age group
bad=bad+sum(~isnan(gather(AXsel.richind.Mean(2:end))));                           % richind wanted only in age group 1
d=[abs(gather(AXsel.assets.CovarianceWith.consumption(2:4))-gather(AXFull.assets.CovarianceWith.consumption(2:4)))'; reshape(abs(gather(AXsel.CorrelationMatrix(1,2,2:3))-gather(AXFull.CorrelationMatrix(1,2,2:3))),[],1); abs(gather(AXsel.CorrelationMatrix(4,4,1))-1); abs(gather(AXsel.richind.Mean(1))-gather(AXFull.richind.Mean(1)))];
fprintf('AutoCorrCrossSec (c) AgeConditionalCrossSectionCovarCorr with the selectors: targeted entries vs full run, this should be zero: %.3e; unwanted entries not NaN, this should be zero: %.3e \n',max(d),bad)

%% (d) round trips through CalibrateLifeCycleModel (grid interpolation layer, as in part 4)
vfoptions_r=struct(); vfoptions_r.gridinterplayer=1; vfoptions_r.ngridinterp=15;
simoptions_r=struct(); simoptions_r.gridinterplayer=1; simoptions_r.ngridinterp=15;
caliboptions_r=struct(); caliboptions_r.verbose=0; caliboptions_r.fminalgo=8;
FnsToEvaluate_r.assets=FnsToEvaluate.assets;
FnsToEvaluate_r.earnings=FnsToEvaluate.earnings;
FnsToEvaluate_r.consumption=FnsToEvaluate.consumption;
% (d1) rho_z from the autocorrelation profile of earnings (earnings are w*kappa_j*z, so this is the autocorrelation of z, informative about rho_z)
Params_r=Params; Params_r.rho_z=0.9; Params_r.sigma_z=0.03; % (the bank's z grid is the Farmer-Toda discretisation at these values)
[z_grid_r,pi_z_r]=CalibFHorz_ExogShockFn(Params_r.rho_z,Params_r.sigma_z);
vfoptions_s=vfoptions_r; vfoptions_s.ExogShockFn=@(rho_z,sigma_z) CalibFHorz_ExogShockFn(rho_z,sigma_z);
simoptions_s=simoptions_r; simoptions_s.ExogShockFn=vfoptions_s.ExogShockFn;
[~,Policy_r]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_r,pi_z_r,ReturnFn,Params_r,DiscountFactorParamNames,[],vfoptions_r);
StationaryDist_r=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_r,n_d,n_a,n_z,N_j,pi_z_r,Params_r,simoptions_r);
ACr=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist_r,Policy_r,FnsToEvaluate_r,Params_r,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_r,pi_z_r,simoptions_r);
TMd1=struct();
t=gather(ACr.earnings.AutoCorrelation); t(Jr-1:end)=NaN; TMd1.AutoCorrTransProbs.earnings.AutoCorrelation=t; % the working ages (earnings are zero in retirement)
caliboptions_d1=caliboptions_r; caliboptions_d1.constrain0to1={'rho_z'};
for ss=1:2
    Params_start=Params_r; lab='started at the truth';
    if ss==2
        Params_start.rho_z=0.7; lab='started at rho_z=0.7 (truth 0.9)';
    end
    ran=0; objv=NaN; perr=NaN; tt=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateLifeCycleModel({'rho_z'},TMd1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_r,pi_z_r,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate_r,caliboptions_d1,vfoptions_s,simoptions_s);
        tt=toc(t0);
        ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.rho_z-0.9);
    catch err
        fprintf('AutoCorrCrossSec (d1) %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ss==1
        fprintf('AutoCorrCrossSec (d1) rho_z from the earnings autocorrelation profile, %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; rho_z moved, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr)
    else
        fprintf('AutoCorrCrossSec (d1) rho_z from the earnings autocorrelation profile, %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; rho_z error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
    end
end
% (d2) beta from the assets-earnings correlation and the age profile of the assets-consumption covariance
[~,Policy_r]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_r);
StationaryDist_r=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_r,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_r);
CSr=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist_r,Policy_r,FnsToEvaluate_r,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_r);
AXr=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist_r,Policy_r,FnsToEvaluate_r,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_r);
TMd2=struct();
TMd2.CrossSectionCovarCorr.assets.CorrelationWith.earnings=gather(CSr.assets.CorrelationWith.earnings);
t=gather(AXr.assets.CovarianceWith.consumption); t(1)=NaN; TMd2.AgeConditionalCrossSectionCovarCorr.assets.CovarianceWith.consumption=t; % (no assets at age 1)
CalibParams_d2=[]; calibsummary_d2=[];
for ss=1:3
    Params_start=Params; lab='started at the truth'; caliboptions_t=caliboptions_r;
    if ss>=2
        Params_start.beta=0.90; lab='started at beta=0.90 (truth 0.95)';
    end
    if ss==3
        caliboptions_t.whichcombos=0; lab=[lab,', whichcombos=0'];
    end
    ran=0; objv=NaN; perr=NaN; tt=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateLifeCycleModel({'beta'},TMd2,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate_r,caliboptions_t,vfoptions_r,simoptions_r);
        tt=toc(t0);
        ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.beta-Params.beta);
    catch err
        fprintf('AutoCorrCrossSec (d2) %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ss==1
        fprintf('AutoCorrCrossSec (d2) beta from the assets-earnings correlation and the assets-consumption covariance profile, %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr)
    elseif ss==2
        fprintf('AutoCorrCrossSec (d2) beta from the assets-earnings correlation and the assets-consumption covariance profile, %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; beta error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
        if ran==1
            CalibParams_d2=CalibParams; calibsummary_d2=calibsummary;
        end
    elseif ran==1 && ~isempty(CalibParams_d2)
        fprintf('AutoCorrCrossSec (d2) whichcombos=0 vs =1 from the same start (%.1f s): beta, this should be zero: %.3e; objective, this should be zero: %.3e; current moments, this should be zero: %.3e \n',tt,abs(CalibParams.beta-CalibParams_d2.beta),abs(calibsummary.objvalue-calibsummary_d2.objvalue),max(abs(calibsummary.currentmomentvec(:)-calibsummary_d2.currentmomentvec(:))))
    else
        fprintf('AutoCorrCrossSec (d2) whichcombos=0 vs =1 from the same start: one of the two did not run, this should be zero: %.3e \n',1)
    end
end
% (d3) logmoments by name on an AutoCorrTransProbs target: the age profile of mean earnings (positive while working) as log targets, the assets-earnings correlation plain
TMd3=struct();
t=gather(ACr.earnings.Mean); t(Jr:end)=NaN; TMd3.AutoCorrTransProbs.earnings.Mean=log(t); % (ACr is at the truth with the bank's grids: rho_z=0.9, sigma_z=0.03)
TMd3.CrossSectionCovarCorr.assets.CorrelationWith.earnings=gather(CSr.assets.CorrelationWith.earnings);
caliboptions_d3=caliboptions_r; caliboptions_d3.logmoments=struct(); caliboptions_d3.logmoments.AutoCorrTransProbs.earnings.Mean=1;
expected=[ones(N_j,1);0]; % one entry per target entry, the NaN (omitted) entries included, as the other logmoments forms (part 3)
ran=0; objv=NaN; perr=NaN; dl=NaN; dm=NaN;
try
    [CalibParams,calibsummary]=CalibrateLifeCycleModel({'beta'},TMd3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate_r,caliboptions_d3,vfoptions_r,simoptions_r);
    ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.beta-Params.beta);
    if length(calibsummary.logmoments)==length(expected)
        dl=max(abs(calibsummary.logmoments(:)-expected));
    else
        dl=abs(length(calibsummary.logmoments)-length(expected));
    end
    expectedm=[log(gather(ACr.earnings.Mean(1:Jr-1)))';TMd3.CrossSectionCovarCorr.assets.CorrelationWith.earnings]; % the current moments are over the targeted entries only
    if length(calibsummary.currentmomentvec)==length(expectedm)
        dm=max(abs(calibsummary.currentmomentvec(:)-expectedm));
    else
        dm=abs(length(calibsummary.currentmomentvec)-length(expectedm));
    end
catch err
    fprintf('AutoCorrCrossSec (d3): the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('AutoCorrCrossSec (d3) logmoments by name on AutoCorrTransProbs.earnings.Mean, started at the truth: did not run, this should be zero: %.3e; logmoments vector (one per target entry, NaN entries included) vs expected, this should be zero: %.3e; current moments (targeted entries) vs the log targets, this should be zero: %.3e; objective, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',1-ran,dl,dm,objv,perr)

%% (e) runtime (informational): a sparse AutoCorr target, whichcombos=1 vs =0
TMe=struct();
t=nan(1,N_j-1); t(3:7)=1; TMe.AutoCorrTransProbs.earnings.AutoCorrelation=t;
simoptions_e=rmfield(simoptions,'agegroupings'); simoptions_e=rmfield(simoptions_e,'timehorizons');
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TMe,FnsToEvaluate,0,N_j,simoptions_e);
caliboptions_e=caliboptions; caliboptions_e.vectoroutput=0; caliboptions_e.weights=ones(sum(~isnan(targetmomentvec)),1); caliboptions_e.selectors=selectors;
argcell={CalibParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,[],Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,[],usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,calibparamsvecindex,calibomitparams_counter,calibomitparamsmatrix};
nreps=3; times=zeros(2,1);
for wc=[1,0]
    caliboptions_t=caliboptions_e; caliboptions_t.whichcombos=wc;
    CalibrateLifeCycleModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions_e); % warm-up
    t0=tic;
    for rep=1:nreps
        CalibrateLifeCycleModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions_e);
    end
    times(2-wc)=toc(t0)/nreps;
end
fprintf('AutoCorrCrossSec (e) runtime per objective evaluation, one AutoCorr target at five start ages (four FnsToEvaluate, two restrictions): whichcombos=1 %.3f s, whichcombos=0 %.3f s, ratio %.2f (informational, the model solve is inside the timing) \n',times(1),times(2),times(1)/times(2))

output=struct();
output.TargetMoments=TMb;
output.selectors=selectors;

end
