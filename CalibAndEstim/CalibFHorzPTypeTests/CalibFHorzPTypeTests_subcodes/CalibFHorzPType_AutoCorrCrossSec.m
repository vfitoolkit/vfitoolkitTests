function output=CalibFHorzPType_AutoCorrCrossSec(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% CalibFHorzPTypeTests part 8: the three target kinds added 2026-10-07/08 with permanent types, TargetMoments.AutoCorrTransProbs
% (EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType), .CrossSectionCovarCorr (EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType) and
% .AgeConditionalCrossSectionCovarCorr (the age-conditional version), alongside AllStats and AgeConditionalStats in the target vector
% (in that order, then CustomModelStats). A target is grouped (fn.stat, restriction.fn.stat, fn.CovarianceWith.fn2, the matrices) or a
% permanent type's own (fn.type.stat, restriction.fn.type.stat, fn.type.CovarianceWith.fn2); the selectors carry a trailing type dimension
% of N_i+1 (the types in Names_i order, then the grouped slot, which makes the commands compute every type at that combination).
%   (a) parsing, no model solved: the target vector against a hand-built one; flags, horizons, sizes; the selectors against hand-built
%       ones with the type slots (AutoCorr [nFns,N_j,1+nRestr,N_i+1], CrossSec [nFns,nFns,1+nRestr,N_i+1], AgeCrossSec
%       [nFns,nFns,ngroups,1+nRestr,N_i+1]; the cross-section targets take conditional restrictions as pages, up to five names);
%       another field order gives the same vector and selectors; a call without the new fields gives empty new outputs; the legacy call
%       with the new fields and fourteen invalid inputs error (incl. the self-correlation, the per-type matrices, per-type N_j, agejshifter)
%   (b) one solve at the truth, full runs of the three PType commands; targets from them in every layout with NaN omissions, with
%       simoptions.agegroupings for the age-conditional kinds and a user simoptions.timehorizons; the objective's moment vector under
%       caliboptions.whichcombos=1 and =0 equals the hand-built one and the targets
%   (c) the selectors handed straight to the three PType commands: a type's targeted entries and the entries a grouped target forces are
%       present and equal to the full run, untargeted start ages / pairs / age groups NaN, untargeted (fn, page, type) absent
%   (d) round trips through CalibrateLifeCycleModel_PType (grid interpolation layer): common rho_z (ExogShockFn) from the grouped earnings
%       autocorrelation profile, at the truth (exact) and from 0.7 (close); beta and a per-type pension from per-type cross-section targets (one among the rich),
%       at the truth (exact) and from a wrong start (close), whichcombos=0 identical (exact); logmoments by name on a four-level name
%   (e) runtime of one objective evaluation with a sparse per-type AutoCorr target, whichcombos=1 vs =0 (informational)
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w,kappa_j_pt,agej,Jr,pension) (1+r)*a+w*kappa_j_pt*z*(agej<Jr)+pension*(agej>=Jr)-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);
FnNames=fieldnames(FnsToEvaluate);
Jr=Params.Jr;
ilow=find(strcmp(Names_i,'low')); ihigh=find(strcmp(Names_i,'high')); igrp=N_i+1;

simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);
simoptions.agegroupings=1:5:N_j;
ngroups=length(simoptions.agegroupings);
npages=3; % unrestricted, rich, poor

%% (a) parsing, no model solved
TM=struct();
TM.AllStats.assets.Mean=1.1;
TM.AgeConditionalStats.assets.low.Mean=[NaN,2,3,4];
t=ones(1,N_j); t(Jr:end)=NaN; TM.AutoCorrTransProbs.earnings.Mean=t;                        % grouped, ages 1..Jr-1
t=ones(1,N_j-1); t(Jr-1:end)=NaN; TM.AutoCorrTransProbs.earnings.low.AutoCorrelation=t;    % low, start ages 1..Jr-2
t=nan(1,N_j-3); t(2:4)=1; TM.AutoCorrTransProbs.earnings.high.AutoCovariance_horizon3=t;         % high, start ages 2..4, horizon 3
t=ones(1,N_j-1); t(1:5)=NaN; TM.AutoCorrTransProbs.rich.assets.AutoCorrelation=t;          % grouped, rich, start ages 6..N_j-1
TM.AutoCorrTransProbs.rich.assets.low.PairMass_horizon2=ones(1,N_j-2);                           % low, rich, start ages 1..N_j-2, horizon 2
TM.AutoCorrTransProbs.poor.consumption.StdDeviation=ones(1,N_j);                           % grouped, poor, every age
TM.CrossSectionCovarCorr.assets.CorrelationWith.earnings=0.3;                               % grouped pair
TM.CrossSectionCovarCorr.assets.high.CovarianceWith.consumption=0.2;                        % high pair
TM.CrossSectionCovarCorr.consumption.low.Mean=1;                                            % low own
TM.CrossSectionCovarCorr.rich.assets.CovarianceWith.consumption=0.1;                        % grouped pair under rich (page 2)
TM.CrossSectionCovarCorr.rich.consumption.low.Mean=1;                                       % low own under rich (page 2)
TM.AgeConditionalCrossSectionCovarCorr.assets.low.CovarianceWith.consumption=[NaN,1,1,1];   % low pair by age group
t=nan(4,4,ngroups); t(1,2,2:3)=1; t(4,4,1)=1; TM.AgeConditionalCrossSectionCovarCorr.CorrelationMatrix=t; % grouped matrix: over every FnsToEvaluate
TM.AgeConditionalCrossSectionCovarCorr.poor.earnings.high.CovarianceWith.consumption=[1,NaN,1,1]; % high pair under poor (page 3): five names
t=nan(4,4,ngroups); t(2,3,1:2)=1; TM.AgeConditionalCrossSectionCovarCorr.rich.CovarianceMatrix=t; % grouped matrix under rich (page 2)

[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,1,N_j,simoptions,Names_i);

hand=[1.1; [NaN;2;3;4]; TM.AutoCorrTransProbs.earnings.Mean(:); TM.AutoCorrTransProbs.earnings.low.AutoCorrelation(:); TM.AutoCorrTransProbs.earnings.high.AutoCovariance_horizon3(:); ...
    TM.AutoCorrTransProbs.rich.assets.AutoCorrelation(:); TM.AutoCorrTransProbs.rich.assets.low.PairMass_horizon2(:); TM.AutoCorrTransProbs.poor.consumption.StdDeviation(:); ...
    0.3; 0.2; 1; 0.1; 1; [NaN;1;1;1]; TM.AgeConditionalCrossSectionCovarCorr.CorrelationMatrix(:); [1;NaN;1;1]; TM.AgeConditionalCrossSectionCovarCorr.rich.CovarianceMatrix(:)];
if length(targetmomentvec)==length(hand)
    d=abs(targetmomentvec-hand); d=d(:); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(hand));
else
    d=1; patternbad=abs(length(targetmomentvec)-length(hand));
end
fprintf('AutoCorrCrossSec PType (a) targetmomentvec vs hand-built (%i entries), this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',length(hand),max([d;0]),patternbad)
handAC=N_j+(N_j-1)+(N_j-3)+(N_j-1)+(N_j-2)+N_j; handX=5; handAX=2*ngroups+32*ngroups;
fprintf('AutoCorrCrossSec PType (a) flags usingautocorr, usingcrosssec, usingagecrosssec = %i %i %i (hand 1 1 1); cumulative sizes end AutoCorr %i (hand %i), CrossSec %i (hand %i), AgeCrossSec %i (hand %i), this should be zero: %.3e \n',usingautocorr,usingcrosssec,usingagecrosssec,autocorrcummomentsizes(end),handAC,crossseccummomentsizes(end),handX,agecrossseccummomentsizes(end),handAX,abs(usingautocorr-1)+abs(usingcrosssec-1)+abs(usingagecrosssec-1)+abs(usingallstats-1)+abs(usinglcp-1)+abs(usingcustomstats)+abs(autocorrcummomentsizes(end)-handAC)+abs(crossseccummomentsizes(end)-handX)+abs(agecrossseccummomentsizes(end)-handAX)+abs(allstatcummomentsizes(end)-1)+abs(acscummomentsizes(end)-4)+abs(cmscummomentsizes))
fprintf('AutoCorrCrossSec PType (a) autocorrtimehorizons [%s] (hand [1 2 3]: the unsuffixed targets are horizon 1), names tables of %i, %i, %i rows and %i, %i, %i columns (hand 6, 5, 4 rows; 4, 5, 5 columns): horizons, this should be zero: %.3e; rows and columns, this should be zero: %.3e \n',num2str(autocorrtimehorizons),size(autocorrmomentnames,1),size(crosssecmomentnames,1),size(agecrosssecmomentnames,1),size(autocorrmomentnames,2),size(crosssecmomentnames,2),size(agecrosssecmomentnames,2),~isequal(autocorrtimehorizons(:)',[1,2,3]),abs(size(autocorrmomentnames,1)-6)+abs(size(crosssecmomentnames,1)-5)+abs(size(agecrosssecmomentnames,1)-4)+abs(size(autocorrmomentnames,2)-4)+abs(size(crosssecmomentnames,2)-5)+abs(size(agecrosssecmomentnames,2)-5))

% the selectors, with the type slots
FnAC=fieldnames(FnsToEvaluate_AutoCorr);
wcAC=zeros(3,N_j,npages,N_i+1);
wcAC(2,1:Jr-1,1,igrp)=1;       % earnings grouped: Mean ages 1..Jr-1
wcAC(2,1:Jr-2,1,ilow)=1;       % earnings low: AutoCorrelation start ages 1..Jr-2
wcAC(2,2:4,1,ihigh)=1;         % earnings high: AutoCovariance_horizon3 start ages 2..4
wcAC(1,6:N_j-1,2,igrp)=1;      % assets grouped under rich: start ages 6..N_j-1
wcAC(1,1:N_j-2,2,ilow)=1;      % assets low under rich: PairMass_horizon2 start ages 1..N_j-2
wcAC(3,:,3,igrp)=1;            % consumption grouped under poor: every age
fnok=isequal(FnAC,{'assets';'earnings';'consumption'});
fprintf('AutoCorrCrossSec PType (a) AutoCorr reduced FnsToEvaluate are assets, earnings, consumption in that order, this should be zero: %.3e \n',1-fnok)
if fnok
    fprintf('AutoCorrCrossSec PType (a) AutoCorr whichcombos [%s] vs hand, this should be zero: %.3e; size, this should be zero: %.3e \n',num2str(size(selectors.AutoCorr.whichcombos)),max(abs(selectors.AutoCorr.whichcombos(:)-wcAC(:))),~isequal(size(selectors.AutoCorr.whichcombos,1:4),[3,N_j,npages,N_i+1]))
end
FnX=fieldnames(FnsToEvaluate_CrossSec);
wcX=zeros(3,3,npages,N_i+1); wcX(1,2,1,igrp)=1; wcX(2,1,1,igrp)=1; wcX(1,3,1,ihigh)=1; wcX(3,1,1,ihigh)=1; wcX(3,3,1,ilow)=1; % page 1
wcX(1,3,2,igrp)=1; wcX(3,1,2,igrp)=1; wcX(3,3,2,ilow)=1; % page 2 (rich)
FnAX=fieldnames(FnsToEvaluate_AgeCrossSec);
wcAX=zeros(4,4,ngroups,npages,N_i+1); wcAX(1,3,2:4,1,ilow)=1; wcAX(3,1,2:4,1,ilow)=1; wcAX(1,2,2:3,1,igrp)=1; wcAX(2,1,2:3,1,igrp)=1; wcAX(4,4,1,1,igrp)=1; % page 1
wcAX(2,3,[1,3,4],3,ihigh)=1; wcAX(3,2,[1,3,4],3,ihigh)=1; % page 3 (poor)
wcAX(2,3,1:2,2,igrp)=1; wcAX(3,2,1:2,2,igrp)=1; % page 2 (rich)
fnok=isequal(FnX,{'assets';'earnings';'consumption'}); fnok2=isequal(FnAX,FnNames);
fprintf('AutoCorrCrossSec PType (a) CrossSec reduced FnsToEvaluate are assets, earnings, consumption, this should be zero: %.3e; AgeCrossSec reduced FnsToEvaluate are all four (a matrix target), this should be zero: %.3e \n',1-fnok,1-fnok2)
if fnok
    fprintf('AutoCorrCrossSec PType (a) CrossSec whichcombos [%s] vs hand, this should be zero: %.3e; size, this should be zero: %.3e \n',num2str(size(selectors.CrossSec.whichcombos)),max(abs(selectors.CrossSec.whichcombos(:)-wcX(:))),~isequal(size(selectors.CrossSec.whichcombos,1:4),[3,3,npages,N_i+1]))
end
if fnok2
    fprintf('AutoCorrCrossSec PType (a) AgeCrossSec whichcombos [%s] vs hand, this should be zero: %.3e; size, this should be zero: %.3e \n',num2str(size(selectors.AgeCrossSec.whichcombos)),max(abs(selectors.AgeCrossSec.whichcombos(:)-wcAX(:))),~isequal(size(selectors.AgeCrossSec.whichcombos,1:5),[4,4,ngroups,npages,N_i+1]))
end

% another field order (the new kinds first, per-type targets before grouped): same vector up to a permutation within a kind, same selectors
TM2=struct();
TM2.AgeConditionalCrossSectionCovarCorr=TM.AgeConditionalCrossSectionCovarCorr;
TM2.CrossSectionCovarCorr.rich=TM.CrossSectionCovarCorr.rich;
TM2.CrossSectionCovarCorr.consumption=TM.CrossSectionCovarCorr.consumption;
TM2.CrossSectionCovarCorr.assets.high=TM.CrossSectionCovarCorr.assets.high;
TM2.CrossSectionCovarCorr.assets.CorrelationWith=TM.CrossSectionCovarCorr.assets.CorrelationWith;
TM2.AutoCorrTransProbs.poor=TM.AutoCorrTransProbs.poor;
TM2.AutoCorrTransProbs.rich.assets.low=TM.AutoCorrTransProbs.rich.assets.low;
TM2.AutoCorrTransProbs.rich.assets.AutoCorrelation=TM.AutoCorrTransProbs.rich.assets.AutoCorrelation;
TM2.AutoCorrTransProbs.earnings.high=TM.AutoCorrTransProbs.earnings.high;
TM2.AutoCorrTransProbs.earnings.low=TM.AutoCorrTransProbs.earnings.low;
TM2.AutoCorrTransProbs.earnings.Mean=TM.AutoCorrTransProbs.earnings.Mean;
TM2.AgeConditionalStats=TM.AgeConditionalStats;
TM2.AllStats=TM.AllStats;
[tv2,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors2,~,~,acsz2,FnAC2,ath2,~,~,xsz2,FnX2,~,~,axsz2,FnAX2]=SetupTargetMoments_FHorz(TM2,FnsToEvaluate,1,N_j,simoptions,Names_i);
bad=~isequal(size(tv2),size(targetmomentvec))+~isequal(fieldnames(FnAC2),FnAC)+~isequal(fieldnames(FnX2),FnX)+~isequal(fieldnames(FnAX2),FnAX)+~isequal(ath2(:),autocorrtimehorizons(:))+~isequal(acsz2(end),autocorrcummomentsizes(end))+~isequal(xsz2(end),crossseccummomentsizes(end))+~isequal(axsz2(end),agecrossseccummomentsizes(end));
maxd=0;
if isequal(size(tv2),size(targetmomentvec))
    s1=sort(tv2(~isnan(tv2))); s2=sort(targetmomentvec(~isnan(targetmomentvec)));
    if length(s1)==length(s2)
        maxd=max([maxd;abs(s1-s2)]);
    else
        bad=bad+1;
    end
end
maxd=max([maxd;abs(selectors2.AutoCorr.whichcombos(:)-selectors.AutoCorr.whichcombos(:));abs(selectors2.CrossSec.whichcombos(:)-selectors.CrossSec.whichcombos(:));abs(selectors2.AgeCrossSec.whichcombos(:)-selectors.AgeCrossSec.whichcombos(:))]);
fprintf('AutoCorrCrossSec PType (a) another field order (new kinds first, per-type before grouped): selectors, horizons and sorted target vector vs the first order, this should be zero: %.3e; sizes or function sets differ, this should be zero: %.3e \n',maxd,bad)

% without the new fields: empty new outputs (four-column AutoCorr and five-column cross-section name tables)
TM0=struct(); TM0.AllStats=TM.AllStats; TM0.AgeConditionalStats=TM.AgeConditionalStats;
[tv0,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors0,ua0,acn0,acs0,FnAC0,ath0,ux0,xn0,xs0,FnX0,uax0,axn0,axs0,FnAX0]=SetupTargetMoments_FHorz(TM0,FnsToEvaluate,1,N_j,simoptions,Names_i);
bad=abs(ua0)+abs(ux0)+abs(uax0)+~isequal(acs0,0)+~isequal(xs0,0)+~isequal(axs0,0)+~isempty(fieldnames(FnAC0))+~isempty(fieldnames(FnX0))+~isempty(fieldnames(FnAX0))+~isempty(ath0);
bad=bad+~(iscell(acn0) && size(acn0,2)==4 && isempty(acn0{1,1}))+~(iscell(xn0) && size(xn0,2)==5 && isempty(xn0{1,1}))+~(iscell(axn0) && size(axn0,2)==5 && isempty(axn0{1,1}));
bad=bad+~isempty(selectors0.AutoCorr.whichcombos)+~isempty(selectors0.CrossSec.whichcombos)+~isempty(selectors0.AgeCrossSec.whichcombos);
d=abs(tv0-targetmomentvec(1:5)); d=d(~isnan(d));
fprintf('AutoCorrCrossSec PType (a) targets without the new fields: new outputs empty and flags zero, this should be zero: %.3e; target vector vs the first five entries, this should be zero: %.3e \n',bad,max([d;0]))

% invalid inputs must error: the legacy call with the new fields, and fourteen bad inputs
errored=0;
try
    SetupTargetMoments_FHorz(TM,FnsToEvaluate,1); % three-input call with the new fields
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AutoCorrTransProbs.earnings.low.Skewness=ones(1,N_j); % unknown statistic
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AutoCorrTransProbs.earnings.AutoCorrelation=ones(1,N_j); % wrong length (1 x N_j-1)
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AutoCorrTransProbs.earnings.low.TransitionProbs=1; % TransitionProbs cannot be targeted
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AutoCorrTransProbs.earnings.low=ones(1,N_j); % a type with no statistic
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AutoCorrTransProbs.wealth.Mean=ones(1,N_j); % unknown name
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.earnings.earnings=1; % the self-correlation is always one, not a target
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.CovarianceMatrix_ptype.low=ones(4); % the per-type matrices cannot be targeted (they shrink to the type's evaluated functions)
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.assets.low.CovarianceWith=1; % no second function
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.CovarianceMatrix=ones(3,3); % a matrix target is over every FnsToEvaluate (four)
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AgeConditionalCrossSectionCovarCorr.assets.CovarianceWith.earnings=ones(1,N_j); % N_j entries but ngroups age groups
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
N_j_struct=struct(); N_j_struct.low=N_j; N_j_struct.high=N_j; % per-type N_j: the AutoCorr PType command has no whichcombos then
try
    SetupTargetMoments_FHorz(TM,FnsToEvaluate,1,N_j_struct,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.rich.assets=1; % a restriction, then a function with no statistic
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.rich=ones(3,3); % a restriction with nothing under it
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.agejshifter=[0,2]; % agejshifter: neither the AutoCorr nor the age-conditional cross-section PType command implements it
try
    SetupTargetMoments_FHorz(TM,FnsToEvaluate,1,N_j,simoptions_bad,Names_i);
catch
    errored=errored+1;
end
fprintf('AutoCorrCrossSec PType (a) the legacy call with the new fields and fourteen invalid inputs must all error, this should be zero: %.3e \n',15-errored)

%% (b) one solve at the truth, full runs of the three PType commands, the objective's moment vector
vfoptions=struct();
[z_gridvals_J,pi_z_J,vfoptions]=ExogShockSetup_FHorz_PType(n_z,z_grid,pi_z,N_j,Names_i,Params,vfoptions,3,0);
simoptions.e_gridvals_J=vfoptions.e_gridvals_J;
simoptions.pi_e_J=vfoptions.pi_e_J;
vfoptions.alreadygridvals=1;
simoptions.alreadygridvals=1;
simoptions.timehorizons=2; % the user's horizon; the targets name horizon 3 as well
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params,simoptions);
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simoptions_full);
simoptions_fullAC=rmfield(simoptions,'agegroupings'); simoptions_fullAC.timehorizons=[1,2,3]; % (the command has no age bins; exactly the horizons listed are computed)
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,simoptions_fullAC);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simoptions);
AXFull=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simoptions);

TMb=struct();
TMb.AllStats.assets.Mean=gather(AllFull.assets.Mean);
t=gather(LCPFull.assets.low.Mean); t(1)=NaN; TMb.AgeConditionalStats.assets.low.Mean=t;
t=gather(ACFull.earnings.Mean); t(Jr:end)=NaN; TMb.AutoCorrTransProbs.earnings.Mean=t;
t=gather(ACFull.earnings.low.AutoCorrelation); t(Jr-1:end)=NaN; TMb.AutoCorrTransProbs.earnings.low.AutoCorrelation=t;
t=gather(ACFull.earnings.high.AutoCovariance_horizon3); t([1,5:end])=NaN; TMb.AutoCorrTransProbs.earnings.high.AutoCovariance_horizon3=t;
t=gather(ACFull.rich.assets.AutoCorrelation); t(1:5)=NaN; TMb.AutoCorrTransProbs.rich.assets.AutoCorrelation=t;
TMb.AutoCorrTransProbs.rich.assets.low.PairMass_horizon2=gather(ACFull.rich.assets.low.PairMass_horizon2);
TMb.AutoCorrTransProbs.poor.consumption.StdDeviation=gather(ACFull.poor.consumption.StdDeviation);
TMb.CrossSectionCovarCorr.assets.CorrelationWith.earnings=gather(CSFull.assets.CorrelationWith.earnings);
TMb.CrossSectionCovarCorr.assets.high.CovarianceWith.consumption=gather(CSFull.assets.high.CovarianceWith.consumption);
TMb.CrossSectionCovarCorr.consumption.low.Mean=gather(CSFull.consumption.low.Mean);
TMb.CrossSectionCovarCorr.rich.assets.CovarianceWith.consumption=gather(CSFull.rich.assets.CovarianceWith.consumption);
TMb.CrossSectionCovarCorr.rich.consumption.low.Mean=gather(CSFull.rich.consumption.low.Mean);
t=gather(AXFull.assets.low.CovarianceWith.consumption); t(1)=NaN; TMb.AgeConditionalCrossSectionCovarCorr.assets.low.CovarianceWith.consumption=t;
t=gather(AXFull.CorrelationMatrix); m=nan(size(t)); m(1,2,2:3)=t(1,2,2:3); m(4,4,1)=t(4,4,1); TMb.AgeConditionalCrossSectionCovarCorr.CorrelationMatrix=m;
t=gather(AXFull.poor.earnings.high.CovarianceWith.consumption); t(2)=NaN; TMb.AgeConditionalCrossSectionCovarCorr.poor.earnings.high.CovarianceWith.consumption=t;
t=gather(AXFull.rich.CovarianceMatrix); m=nan(size(t)); m(2,3,1:2)=t(2,3,1:2); TMb.AgeConditionalCrossSectionCovarCorr.rich.CovarianceMatrix=m;

% the hand-built vector: walk TargetMoments kind by kind (one to five names) and read the same field of the full run (same nesting)
kinds={'AllStats','AgeConditionalStats','AutoCorrTransProbs','CrossSectionCovarCorr','AgeConditionalCrossSectionCovarCorr'};
fulls={AllFull,LCPFull,ACFull,CSFull,AXFull};
hand=[]; tvec=[];
for kk=1:5
    TMk=TMb.(kinds{kk}); Fk=fulls{kk};
    a1vec=fieldnames(TMk);
    for a1=1:length(a1vec)
        t1=TMk.(a1vec{a1});
        if ~isstruct(t1)
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
                        t3=t2.(a3vec{a3});
                        if ~isstruct(t3)
                            v=gather(Fk.(a1vec{a1}).(a2vec{a2}).(a3vec{a3})); hand=[hand; v(:)]; tvec=[tvec; t3(:)];
                        else
                            a4vec=fieldnames(t3);
                            for a4=1:length(a4vec)
                                t4=t3.(a4vec{a4});
                                if ~isstruct(t4)
                                    v=gather(Fk.(a1vec{a1}).(a2vec{a2}).(a3vec{a3}).(a4vec{a4})); hand=[hand; v(:)]; tvec=[tvec; t4(:)];
                                else
                                    a5vec=fieldnames(t4);
                                    for a5=1:length(a5vec)
                                        v=gather(Fk.(a1vec{a1}).(a2vec{a2}).(a3vec{a3}).(a4vec{a4}).(a5vec{a5})); hand=[hand; v(:)]; t5=t4.(a5vec{a5}); tvec=[tvec; t5(:)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
actual=~isnan(tvec);
fprintf('AutoCorrCrossSec PType (b) %i target entries, %i of them targeted (the rest NaN); the model is NaN at a targeted entry, this should be zero: %.3e \n',length(tvec),sum(actual),sum(isnan(hand(actual))))

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
fprintf('AutoCorrCrossSec PType (b) SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',max([d;0]),patternbad)
handA=hand(actual); tA=targetmomentvec(actual);
for wc=[1,0]
    caliboptions_t=caliboptions; caliboptions_t.whichcombos=wc;
    ran=0; m=[];
    try
        m=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions));
        ran=1;
    catch err
        fprintf('AutoCorrCrossSec PType (b) whichcombos=%i: the objective errored with: %s \n',wc,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ran==1 && length(m)==length(handA)
        d=max([abs(m(:)-handA(:));0]); dt=max([abs(m(:)-tA(:));0]); nanbad=sum(isnan(m));
    else
        d=1; dt=1; nanbad=1;
    end
    fprintf('AutoCorrCrossSec PType (b) whichcombos=%i: did not run, this should be zero: %.3e; moment vector vs hand-built, this should be zero: %.3e; vs the targets (targets are the model), this should be zero: %.3e; NaN entries or length mismatch, this should be zero: %.3e \n',wc,1-ran,d,dt,nanbad)
end

%% (c) the selectors handed straight to the three PType commands
simsel=simoptions_fullAC; simsel.whichcombos=selectors.AutoCorr.whichcombos;
ACsel=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate_AutoCorr,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,simsel);
bad=0; d=[];
% grouped earnings (start ages 1..Jr-1 from the Mean target): AutoCorrelation NaN from Jr; the per-type earnings are forced there, with low's own start ages inside that range
bad=bad+sum(~isnan(gather(ACsel.earnings.AutoCorrelation(Jr:end))));
bad=bad+sum(~isnan(gather(ACsel.earnings.low.AutoCorrelation(Jr:end))));
bad=bad+sum(~isnan(gather(ACsel.earnings.high.AutoCovariance_horizon3(Jr:end))));
d=[d; abs(gather(ACsel.earnings.Mean(1:Jr-1))-gather(ACFull.earnings.Mean(1:Jr-1)))'];
d=[d; abs(gather(ACsel.earnings.low.AutoCorrelation(1:Jr-2))-gather(ACFull.earnings.low.AutoCorrelation(1:Jr-2)))'];
d=[d; abs(gather(ACsel.earnings.high.AutoCovariance_horizon3(2:4))-gather(ACFull.earnings.high.AutoCovariance_horizon3(2:4)))'];
% assets: unrestricted nothing (no struct at all), under rich grouped start ages 6..N_j-1 (NaN at 1..5, the high type forced there), low's PairMass_horizon2 at every start age
bad=bad+isfield(ACsel,'assets')+isfield(ACsel,'consumption')+isfield(ACsel,'richind');
bad=bad+sum(~isnan(gather(ACsel.rich.assets.AutoCorrelation(1:5))))+sum(~isnan(gather(ACsel.rich.assets.high.AutoCorrelation(1:5))));
d=[d; abs(gather(ACsel.rich.assets.AutoCorrelation(6:end))-gather(ACFull.rich.assets.AutoCorrelation(6:end)))'];
d=[d; abs(gather(ACsel.rich.assets.high.AutoCorrelation(6:end))-gather(ACFull.rich.assets.high.AutoCorrelation(6:end)))'];
d=[d; abs(gather(ACsel.rich.assets.low.PairMass_horizon2)-gather(ACFull.rich.assets.low.PairMass_horizon2))'];
% consumption under poor: grouped at every age, both types forced
d=[d; abs(gather(ACsel.poor.consumption.StdDeviation)-gather(ACFull.poor.consumption.StdDeviation))'];
bad=bad+~isfield(ACsel.poor.consumption,'low')+~isfield(ACsel.poor.consumption,'high');
% nothing else present
bad=bad+isfield(ACsel.rich,'earnings')+isfield(ACsel.rich,'consumption')+isfield(ACsel.poor,'assets')+isfield(ACsel.poor,'earnings');
d=d(~isnan(d));
fprintf('AutoCorrCrossSec PType (c) AutoCorrTransProbs PType with the selectors: targeted and forced entries vs full run, this should be zero: %.3e; untargeted entries not NaN or present, forced entries missing, this should be zero: %.3e; horizon fields k2 and k3 both present, this should be zero: %.3e \n',max([d;0]),bad,~isfield(ACsel.earnings,'AutoCovariance_horizon2')+~isfield(ACsel.earnings,'AutoCovariance_horizon3'))
simselX=simoptions; simselX.whichcombos=selectors.CrossSec.whichcombos;
CSsel=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate_CrossSec,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simselX);
bad=0; d=[];
d=[d; abs(gather(CSsel.assets.CorrelationWith.earnings)-gather(CSFull.assets.CorrelationWith.earnings))]; % grouped pair
bad=bad+~isnan(gather(CSsel.assets.CorrelationWith.consumption));   % grouped pair not asked
bad=bad+~isnan(gather(CSsel.consumption.Mean));                      % grouped own stats of consumption not asked
d=[d; abs(gather(CSsel.consumption.low.Mean)-gather(CSFull.consumption.low.Mean))];                                        % low's own target
d=[d; abs(gather(CSsel.assets.low.CorrelationWith.earnings)-gather(CSFull.assets.low.CorrelationWith.earnings))];           % forced on low by the grouped pair
bad=bad+~isnan(gather(CSsel.assets.low.CovarianceWith.consumption)); % not asked for low (both evaluated for low)
d=[d; abs(gather(CSsel.assets.high.CovarianceWith.consumption)-gather(CSFull.assets.high.CovarianceWith.consumption))];     % high's own target
d=[d; abs(gather(CSsel.assets.high.CorrelationWith.earnings)-gather(CSFull.assets.high.CorrelationWith.earnings))];         % forced on high
bad=bad+~isnan(gather(CSsel.earnings.high.CovarianceWith.consumption)); % not asked for high (both evaluated for high)
bad=bad+isfield(CSsel,'richind');
% the rich page: the grouped pair (assets, consumption) forces both types there, low's own consumption Mean; earnings has nothing on it; the poor page has nothing at all (the grouped leaves of every page are always present, NaN where not asked)
d=[d; abs(gather(CSsel.rich.assets.CovarianceWith.consumption)-gather(CSFull.rich.assets.CovarianceWith.consumption)); abs(gather(CSsel.rich.consumption.low.Mean)-gather(CSFull.rich.consumption.low.Mean))];
bad=bad+~isfield(CSsel.rich.assets,'low')+~isfield(CSsel.rich.assets,'high'); % forced on both types by the grouped pair
bad=bad+isfield(CSsel.rich.earnings,'low')+isfield(CSsel.rich.earnings,'high')+~isnan(gather(CSsel.rich.earnings.Mean)); % earnings: no type evaluates it on the rich page, the grouped leaves are there but NaN
bad=bad+isfield(CSsel.poor.assets,'low')+isfield(CSsel.poor.assets,'high')+~isnan(gather(CSsel.poor.assets.Mean))+~isnan(gather(CSsel.poor.assets.CovarianceWith.consumption)); % the poor page: nothing selected, so no per-type output and NaN grouped leaves
fprintf('AutoCorrCrossSec PType (c) CrossSectionCovarCorr PType with the selectors (incl. the rich page): targeted and forced entries vs full run, this should be zero: %.3e; unwanted pairs or own stats not NaN, function present, this should be zero: %.3e \n',max(d),bad)
simselAX=simoptions; simselAX.whichcombos=selectors.AgeCrossSec.whichcombos;
AXsel=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate_AgeCrossSec,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simselAX);
bad=0; d=[];
d=[d; reshape(abs(gather(AXsel.CorrelationMatrix(1,2,2:3))-gather(AXFull.CorrelationMatrix(1,2,2:3))),[],1); abs(gather(AXsel.CorrelationMatrix(4,4,1))-1)];
bad=bad+~isnan(gather(AXsel.CorrelationMatrix(1,2,1)))+~isnan(gather(AXsel.CorrelationMatrix(1,2,4)))+sum(~isnan(gather(AXsel.CorrelationMatrix(2,3,:))));
d=[d; abs(gather(AXsel.richind.Mean(1))-gather(AXFull.richind.Mean(1)))]; bad=bad+sum(~isnan(gather(AXsel.richind.Mean(2:end)))); % grouped own stats of richind in age group 1 only
d=[d; abs(gather(AXsel.assets.low.CovarianceWith.consumption(2:4))-gather(AXFull.assets.low.CovarianceWith.consumption(2:4)))']; % low's own target
bad=bad+~isnan(gather(AXsel.assets.low.CovarianceWith.consumption(1)));
d=[d; abs(gather(AXsel.assets.low.CovarianceWith.earnings(2:3))-gather(AXFull.assets.low.CovarianceWith.earnings(2:3)))'];     % forced on low by the grouped matrix
bad=bad+sum(~isnan(gather(AXsel.assets.low.CovarianceWith.earnings([1,4]))));
d=[d; abs(gather(AXsel.assets.high.CovarianceWith.earnings(2:3))-gather(AXFull.assets.high.CovarianceWith.earnings(2:3)))'];   % forced on high
bad=bad+sum(~isnan(gather(AXsel.assets.high.CovarianceWith.earnings([1,4]))));
bad=bad+(isfield(AXsel.consumption,'high') && any(~isnan(gather(AXsel.consumption.high.Mean)))); % consumption is not evaluated for high on page 1 (the age-conditional command preallocates every function on every page, so the struct may exist, all NaN)
% the poor page: high's pair (earnings, consumption) at age groups 1, 3, 4 (NaN at 2), low not evaluated there (its struct, if present, all NaN); the rich page: the grouped matrix entries (2,3) at groups 1:2 only
d=[d; abs(gather(AXsel.poor.earnings.high.CovarianceWith.consumption([1,3,4]))-gather(AXFull.poor.earnings.high.CovarianceWith.consumption([1,3,4])))'; reshape(abs(gather(AXsel.rich.CovarianceMatrix(2,3,1:2))-gather(AXFull.rich.CovarianceMatrix(2,3,1:2))),[],1)];
bad=bad+~isnan(gather(AXsel.poor.earnings.high.CovarianceWith.consumption(2)))+(isfield(AXsel.poor.earnings,'low') && any(~isnan(gather(AXsel.poor.earnings.low.Mean))))+sum(~isnan(gather(AXsel.rich.CovarianceMatrix(2,3,3:4))))+sum(~isnan(gather(AXsel.rich.CovarianceMatrix(1,2,:))));
d=d(~isnan(d));
fprintf('AutoCorrCrossSec PType (c) AgeConditionalCrossSectionCovarCorr PType with the selectors (incl. the rich and poor pages): targeted and forced entries vs full run, this should be zero: %.3e; unwanted entries not NaN, unevaluated type present, this should be zero: %.3e \n',max([d;0]),bad)

%% (d) round trips through CalibrateLifeCycleModel_PType (grid interpolation layer)
vfoptions_r=struct(); vfoptions_r.gridinterplayer=1; vfoptions_r.ngridinterp=15;
simoptions_r=struct(); simoptions_r.gridinterplayer=1; simoptions_r.ngridinterp=15;
caliboptions_r=struct(); caliboptions_r.verbose=0; caliboptions_r.fminalgo=8;
FnsToEvaluate_r.assets=FnsToEvaluate.assets;
FnsToEvaluate_r.earnings=FnsToEvaluate.earnings;
FnsToEvaluate_r.consumption=FnsToEvaluate.consumption;
% (d1) common rho_z from the grouped autocorrelation profile of earnings (earnings are w*kappa_j_pt*z, so this is the pooled autocorrelation of z)
Params_r=Params; Params_r.rho_z=0.9; Params_r.sigma_z=0.03;
[z_grid_r,pi_z_r]=CalibFHorzPType_ExogShockFn(Params_r.rho_z,Params_r.sigma_z);
vfoptions_s=vfoptions_r; vfoptions_s.ExogShockFn=@(rho_z,sigma_z) CalibFHorzPType_ExogShockFn(rho_z,sigma_z);
simoptions_s=simoptions_r; simoptions_s.ExogShockFn=vfoptions_s.ExogShockFn;
[~,Policy_r]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid_r,pi_z_r,ReturnFn,Params_r,DiscountFactorParamNames,vfoptions_r);
StationaryDist_r=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_r,n_d,n_a,n_z,N_j,Names_i,pi_z_r,Params_r,simoptions_r);
ACr=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_r,Policy_r,FnsToEvaluate_r,Params_r,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid_r,pi_z_r,simoptions_r);
TMd1=struct();
t=gather(ACr.earnings.AutoCorrelation); t(Jr-1:end)=NaN; TMd1.AutoCorrTransProbs.earnings.AutoCorrelation=t; % the working ages
caliboptions_d1=caliboptions_r; caliboptions_d1.constrain0to1={'rho_z'};
for ss=1:2
    Params_start=Params_r; lab='started at the truth';
    if ss==2
        Params_start.rho_z=0.7; lab='started at rho_z=0.7 (truth 0.9)';
    end
    ran=0; objv=NaN; perr=NaN; tt=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'rho_z'},TMd1,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid_r,pi_z_r,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate_r,caliboptions_d1,vfoptions_s,simoptions_s);
        tt=toc(t0);
        ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.rho_z-0.9);
    catch err
        fprintf('AutoCorrCrossSec PType (d1) %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ss==1
        fprintf('AutoCorrCrossSec PType (d1) common rho_z from the grouped earnings autocorrelation profile, %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; rho_z moved, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr)
    else
        fprintf('AutoCorrCrossSec PType (d1) common rho_z from the grouped earnings autocorrelation profile, %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; rho_z error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
    end
end
% (d2) beta (common) and a per-type pension (truth low 0.5, high 0.4) from per-type cross-section targets: the assets-earnings correlation of
% each type and the age profile of the assets-consumption covariance of each type
Params_pt=Params; Params_pt.pension=struct(); Params_pt.pension.low=0.5; Params_pt.pension.high=0.4;
[~,Policy_pt]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_pt,DiscountFactorParamNames,vfoptions_r);
StationaryDist_pt=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_pt,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_pt,simoptions_r);
simoptions_r2=simoptions_r; simoptions_r2.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5); % one restriction: a restricted per-type cross-section target too
CSpt=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_pt,Policy_pt,FnsToEvaluate_r,Params_pt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_r2);
AXpt=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_pt,Policy_pt,FnsToEvaluate_r,Params_pt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_r);
TMd2=struct();
TMd2.CrossSectionCovarCorr.assets.low.CorrelationWith.earnings=gather(CSpt.assets.low.CorrelationWith.earnings);
TMd2.CrossSectionCovarCorr.assets.high.CorrelationWith.earnings=gather(CSpt.assets.high.CorrelationWith.earnings);
TMd2.CrossSectionCovarCorr.rich.assets.low.CorrelationWith.earnings=gather(CSpt.rich.assets.low.CorrelationWith.earnings);
t=gather(AXpt.assets.low.CovarianceWith.consumption); t(1)=NaN; TMd2.AgeConditionalCrossSectionCovarCorr.assets.low.CovarianceWith.consumption=t; % (no assets at age 1)
t=gather(AXpt.assets.high.CovarianceWith.consumption); t(1)=NaN; TMd2.AgeConditionalCrossSectionCovarCorr.assets.high.CovarianceWith.consumption=t;
truth_pt=[0.95;0.5;0.4];
CalibParams_d2=[]; calibsummary_d2=[];
for ss=1:3
    Params_start=Params_pt; lab='started at the truth'; caliboptions_t=caliboptions_r;
    if ss>=2
        Params_start.beta=0.90; Params_start.pension.low=0.30; Params_start.pension.high=0.30; lab='started at beta=0.90, pensions 0.30 (truth 0.95, 0.5, 0.4)';
    end
    if ss==3
        caliboptions_t.whichcombos=0; lab=[lab,', whichcombos=0'];
    end
    ran=0; objv=NaN; perr=NaN; tt=NaN; rec=nan(3,1);
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'beta','pension'},TMd2,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate_r,caliboptions_t,vfoptions_r,simoptions_r2);
        tt=toc(t0);
        ran=1; objv=calibsummary.objvalue;
        if isstruct(CalibParams.pension)
            rec=[CalibParams.beta;CalibParams.pension.low;CalibParams.pension.high];
        end
        perr=max(abs(rec-truth_pt));
    catch err
        fprintf('AutoCorrCrossSec PType (d2) %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ss==1
        fprintf('AutoCorrCrossSec PType (d2) beta and a per-type pension from per-type cross-section targets (one among the rich), %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; max parameter move, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr)
    elseif ss==2
        fprintf('AutoCorrCrossSec PType (d2) beta and a per-type pension from per-type cross-section targets (one among the rich), %s (%.1f s), recovered %.6f, %.6f, %.6f, objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',lab,tt,rec(1),rec(2),rec(3),objv,1-ran,perr)
        if ran==1
            CalibParams_d2=rec; calibsummary_d2=calibsummary;
        end
    elseif ran==1 && ~isempty(CalibParams_d2)
        fprintf('AutoCorrCrossSec PType (d2) whichcombos=0 vs =1 from the same start (%.1f s): parameters, this should be zero: %.3e; objective, this should be zero: %.3e; current moments, this should be zero: %.3e \n',tt,max(abs(rec-CalibParams_d2)),abs(calibsummary.objvalue-calibsummary_d2.objvalue),max(abs(calibsummary.currentmomentvec(:)-calibsummary_d2.currentmomentvec(:))))
    else
        fprintf('AutoCorrCrossSec PType (d2) whichcombos=0 vs =1 from the same start: one of the two did not run, this should be zero: %.3e \n',1)
    end
end
% (d3) logmoments by name on a four-level name: the low type's age profile of mean earnings (positive while working) as log targets, a grouped
% correlation plain; at the truth (common parameters, the bank's grids)
TMd3=struct();
t=gather(ACr.earnings.low.Mean); t(Jr:end)=NaN; TMd3.AutoCorrTransProbs.earnings.low.Mean=log(t); % (ACr is at the truth: rho_z=0.9, sigma_z=0.03 are the bank's grid)
CSr=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_r,Policy_r,FnsToEvaluate_r,Params_r,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid_r,simoptions_r);
TMd3.CrossSectionCovarCorr.assets.CorrelationWith.earnings=gather(CSr.assets.CorrelationWith.earnings);
caliboptions_d3=caliboptions_r; caliboptions_d3.logmoments=struct(); caliboptions_d3.logmoments.AutoCorrTransProbs.earnings.low.Mean=1;
expected=[ones(N_j,1);0]; % one entry per target entry, the NaN (omitted) entries included (part 3's convention)
expectedm=[log(gather(ACr.earnings.low.Mean(1:Jr-1)))';TMd3.CrossSectionCovarCorr.assets.CorrelationWith.earnings]; % the current moments are over the targeted entries only
ran=0; objv=NaN; perr=NaN; dl=NaN; dm=NaN;
try
    [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'beta'},TMd3,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate_r,caliboptions_d3,vfoptions_r,simoptions_r);
    ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.beta-Params.beta);
    if length(calibsummary.logmoments)==length(expected)
        dl=max(abs(calibsummary.logmoments(:)-expected));
    else
        dl=abs(length(calibsummary.logmoments)-length(expected));
    end
    if length(calibsummary.currentmomentvec)==length(expectedm)
        dm=max(abs(calibsummary.currentmomentvec(:)-expectedm));
    else
        dm=abs(length(calibsummary.currentmomentvec)-length(expectedm));
    end
catch err
    fprintf('AutoCorrCrossSec PType (d3): the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('AutoCorrCrossSec PType (d3) logmoments by name on AutoCorrTransProbs.earnings.low.Mean (four levels), started at the truth: did not run, this should be zero: %.3e; logmoments vector (one per target entry, NaN entries included) vs expected, this should be zero: %.3e; current moments (targeted entries) vs the log targets, this should be zero: %.3e; objective, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',1-ran,dl,dm,objv,perr)

%% (e) runtime (informational): a sparse per-type AutoCorr target, whichcombos=1 vs =0
TMe=struct();
t=nan(1,N_j-1); t(3:7)=1; TMe.AutoCorrTransProbs.earnings.low.AutoCorrelation=t;
simoptions_e=rmfield(simoptions,'agegroupings'); simoptions_e=rmfield(simoptions_e,'timehorizons');
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TMe,FnsToEvaluate,1,N_j,simoptions_e,Names_i);
caliboptions_e=caliboptions; caliboptions_e.vectoroutput=0; caliboptions_e.weights=ones(sum(~isnan(targetmomentvec)),1); caliboptions_e.selectors=selectors;
argcell={CalibParamNames,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,nCalibParams,nCalibParamsFinder,calibparamsvecindex,calibparamssizes,calibomitparams_counter,calibomitparamsmatrix};
nreps=3; times=zeros(2,1);
for wc=[1,0]
    caliboptions_t=caliboptions_e; caliboptions_t.whichcombos=wc;
    CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions_e); % warm-up
    t0=tic;
    for rep=1:nreps
        CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_t,vfoptions,simoptions_e);
    end
    times(2-wc)=toc(t0)/nreps;
end
fprintf('AutoCorrCrossSec PType (e) runtime per objective evaluation, one per-type AutoCorr target at five start ages (four FnsToEvaluate, two restrictions, two types): whichcombos=1 %.3f s, whichcombos=0 %.3f s, ratio %.2f (informational, the model solve is inside the timing) \n',times(1),times(2),times(1)/times(2))

output=struct();
output.TargetMoments=TMb;
output.selectors=selectors;

end
