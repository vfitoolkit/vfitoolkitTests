function output=CalibFHorzPType_TargetParsing(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% CalibFHorzPTypeTests part 1: SetupTargetMoments_FHorz with permanent types (useptype=1). No model is solved.
%
% TargetMoments with arbitrary values of the right sizes, in every layout: fn.stat (grouped), fn.type.stat, restriction.fn.stat,
% restriction.fn.type.stat, fn.MoreInequality.substat, NaN-omitted age groups, matrix-valued age-conditional targets, agegroupings.
% The selectors carry a trailing type dimension of N_i+1 (types in Names_i order, then the grouped slot).
% Checks (all exact):
%   (a) targetmomentvec equals the hand-built concatenation; sizes
%   (b) AllStats selectors [nFns,1+nRestr,N_i+1] and [..,7] equal the hand-built arrays
%   (c) AgeConditionalStats selectors [nFns,ngroups,1+nRestr,N_i+1] and [..,7] equal the hand-built arrays (ages off where NaN;
%       LorenzCurve overriding a Gini-only 3 within a slot)
%   (d) per-type agegroupings of different lengths: the age axis is the longest, a shorter type fills its own groups;
%       agejshifter: a type's groups sit shifted on the common axis
%   (e) the union 7-vectors
%   (f) two other field orders give the same selectors, unions and function sets, and a permutation of the target vector
%   (g) the three-input legacy call (no N_j, simoptions, Names_i) gives the same unions and function sets and empty selectors
%   (h) five invalid inputs must error: an unknown second-level name, a type with no statistic, an unknown function, a per-type
%       age-conditional target of the wrong length, a name that is both a function and a restriction
% (the model inputs are the bank's standard ones and unused here)

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w,kappa_j_pt,agej,Jr,pension) (1+r)*a+w*kappa_j_pt*z*(agej<Jr)+pension*(agej>=Jr)-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);
simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);
simoptions.agegroupings=1:5:N_j;
ngroups=length(simoptions.agegroupings);
npages=3; % unrestricted, rich, poor
ilow=1; ihigh=2; igrp=N_i+1; % the slots (Names_i={'low','high'})

%% The targets (arbitrary values of the right sizes)
TM=struct();
TM.AllStats.earnings.Gini=0.4;                                   % grouped
TM.AllStats.earnings.low.LorenzCurve=linspace(0,1,100)';         % type low
TM.AllStats.poor.assets.StdDeviation=0.2;                        % restricted, grouped
TM.AllStats.assets.Mean=1.1;                                     % grouped
TM.AllStats.assets.high.QuantileCutoffs=(1:21)';                 % type high
TM.AllStats.assets.MoreInequality.Top1share=0.1;                 % grouped, three levels
TM.AllStats.rich.consumption.low.Minimum=0.3;                    % restricted, type low (four levels)
TM.AllStats.rich.consumption.RatioMeanToMedian=1.05;             % restricted, grouped
TM.AgeConditionalStats.assets.low.Mean=[NaN,2,3,4];
TM.AgeConditionalStats.assets.Gini=[NaN,0.3,0.3,NaN];
TM.AgeConditionalStats.assets.high.LorenzCurve=[nan(100,1),ones(100,3)];
TM.AgeConditionalStats.rich.earnings.Mean=[NaN,NaN,1,1];
TM.AgeConditionalStats.earnings.low.QuantileMeans=[ones(20,2),nan(20,2)];
TM.AgeConditionalStats.poor.consumption.high.Variance=[1,1,1,1];

[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,1,N_j,simoptions,Names_i);

%% (a) the target vector
hand=[0.4; linspace(0,1,100)'; 0.2; 1.1; (1:21)'; 0.1; 0.3; 1.05; ...
    [NaN;2;3;4]; [NaN;0.3;0.3;NaN]; reshape([nan(100,1),ones(100,3)],[],1); [NaN;NaN;1;1]; reshape([ones(20,2),nan(20,2)],[],1); [1;1;1;1]];
if length(targetmomentvec)==length(hand)
    d=abs(targetmomentvec-hand); d=d(:); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(hand));
else
    d=1; patternbad=abs(length(targetmomentvec)-length(hand));
end
fprintf('TargetParsing PType (a) targetmomentvec vs hand-built, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',max([d;0]),patternbad)
fprintf('TargetParsing PType (a) flags usingallstats, usinglcp, usingcustomstats = %i %i %i; cumulative sizes AllStats end %i (hand 127), AgeConditionalStats end %i (hand 496), name tables have four columns: this should be zero: %.3e \n',usingallstats,usinglcp,usingcustomstats,allstatcummomentsizes(end),acscummomentsizes(end),abs(usingallstats-1)+abs(usinglcp-1)+abs(usingcustomstats)+abs(allstatcummomentsizes(end)-127)+abs(acscummomentsizes(end)-496)+abs(cmscummomentsizes)+(~isempty(cmsmomentnames{1}))+abs(size(allstatmomentnames,2)-4)+abs(size(acsmomentnames,2)-4))

%% (b) AllStats selectors: functions assets, earnings, consumption (FnsToEvaluate order; richind is not targeted)
FnA=fieldnames(FnsToEvaluate_AllStats);
wcA=zeros(3,npages,N_i+1); wsA=zeros(3,npages,N_i+1,7);
wsA(1,1,igrp,[1,7])=1;   % assets grouped unrestricted: Mean, MoreInequality
wsA(1,1,ihigh,6)=1;      % assets type high unrestricted: QuantileCutoffs
wsA(1,3,igrp,3)=1;       % assets grouped under poor: StdDeviation
wsA(2,1,igrp,4)=3;       % earnings grouped unrestricted: Gini only
wsA(2,1,ilow,4)=1;       % earnings type low unrestricted: LorenzCurve
wsA(3,2,ilow,5)=1;       % consumption type low under rich: Minimum
wsA(3,2,igrp,[1,2])=1;   % consumption grouped under rich: RatioMeanToMedian
wcA(1,1,igrp)=1; wcA(1,1,ihigh)=1; wcA(1,3,igrp)=1; wcA(2,1,igrp)=1; wcA(2,1,ilow)=1; wcA(3,2,ilow)=1; wcA(3,2,igrp)=1;
fnok=isequal(FnA,{'assets';'earnings';'consumption'});
fprintf('TargetParsing PType (b) AllStats reduced FnsToEvaluate are assets, earnings, consumption in that order, this should be zero: %.3e \n',1-fnok)
if fnok
    fprintf('TargetParsing PType (b) AllStats selectors with the type dimension, sizes [%s] and [%s]: whichcombos vs hand, this should be zero: %.3e; whichstats vs hand, this should be zero: %.3e; wrong sizes, this should be zero: %.3e \n',num2str(size(selectors.AllStats.whichcombos)),num2str(size(selectors.AllStats.whichstats)),max(abs(selectors.AllStats.whichcombos(:)-wcA(:))),max(abs(selectors.AllStats.whichstats(:)-wsA(:))),~isequal(size(selectors.AllStats.whichcombos),[3,npages,N_i+1])+~isequal(size(selectors.AllStats.whichstats),[3,npages,N_i+1,7]))
end

%% (c) AgeConditionalStats selectors
FnL=fieldnames(FnsToEvaluate_ACStats);
wcL=zeros(3,ngroups,npages,N_i+1); wsL=zeros(3,ngroups,npages,N_i+1,7);
wsL(1,2:4,1,ilow,1)=1;     % assets type low: Mean at ages 2-4
wsL(1,2:3,1,igrp,4)=3;     % assets grouped: Gini at ages 2,3
wsL(1,2:4,1,ihigh,4)=2;    % assets type high: LorenzCurve at ages 2-4
wsL(2,3:4,2,igrp,1)=1;     % earnings grouped under rich: Mean at ages 3,4
wsL(2,1:2,1,ilow,6)=1;     % earnings type low: QuantileMeans at ages 1,2
wsL(3,:,3,ihigh,3)=1;      % consumption type high under poor: Variance at every age
wcL(1,2:4,1,ilow)=1; wcL(1,2:3,1,igrp)=1; wcL(1,2:4,1,ihigh)=1; wcL(2,3:4,2,igrp)=1; wcL(2,1:2,1,ilow)=1; wcL(3,:,3,ihigh)=1;
fnok=isequal(FnL,{'assets';'earnings';'consumption'});
fprintf('TargetParsing PType (c) AgeConditionalStats reduced FnsToEvaluate are assets, earnings, consumption in that order, this should be zero: %.3e \n',1-fnok)
if fnok
    fprintf('TargetParsing PType (c) AgeConditionalStats selectors with the type dimension, sizes [%s] and [%s]: whichcombos vs hand, this should be zero: %.3e; whichstats vs hand, this should be zero: %.3e; wrong sizes, this should be zero: %.3e \n',num2str(size(selectors.ACStats.whichcombos)),num2str(size(selectors.ACStats.whichstats)),max(abs(selectors.ACStats.whichcombos(:)-wcL(:))),max(abs(selectors.ACStats.whichstats(:)-wsL(:))),~isequal(size(selectors.ACStats.whichcombos),[3,ngroups,npages,N_i+1])+~isequal(size(selectors.ACStats.whichstats),[3,ngroups,npages,N_i+1,7]))
end

%% (d) per-type agegroupings of different lengths, and agejshifter
simoptions_d=simoptions;
simoptions_d.agegroupings=struct();
simoptions_d.agegroupings.low=1:5:N_j;   % 4 groups
simoptions_d.agegroupings.high=1:10:N_j; % 2 groups: the age axis has 4, the high type fills the first 2
TMd=struct();
TMd.AgeConditionalStats.assets.low.Mean=[1,NaN,3,4];
TMd.AgeConditionalStats.assets.high.Mean=[1,2];
TMd.AgeConditionalStats.assets.Mean=[1,2,3,4];
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors_d]=SetupTargetMoments_FHorz(TMd,FnsToEvaluate,1,N_j,simoptions_d,Names_i);
wcd=zeros(1,4,npages,N_i+1); wcd(1,[1,3,4],1,ilow)=1; wcd(1,1:2,1,ihigh)=1; wcd(1,:,1,igrp)=1;
if isequal(size(selectors_d.ACStats.whichcombos),size(wcd))
    dd=max(abs(selectors_d.ACStats.whichcombos(:)-wcd(:)));
else
    dd=1;
end
fprintf('TargetParsing PType (d) per-type agegroupings (low 4 groups, high 2): selectors sized by the longest, size [%s], whichcombos vs hand, this should be zero: %.3e \n',num2str(size(selectors_d.ACStats.whichcombos)),dd)
simoptions_s=simoptions; simoptions_s=rmfield(simoptions_s,'agegroupings'); % default: each age separately, N_j groups
simoptions_s.agejshifter=struct(); simoptions_s.agejshifter.low=0; simoptions_s.agejshifter.high=2; % the high type starts two ages later on the common axis
TMs=struct();
TMs.AgeConditionalStats.assets.high.Mean=ones(1,N_j); % its N_j ages sit at positions 3..N_j+2, the last two beyond the axis
TMs.AgeConditionalStats.assets.low.Mean=ones(1,N_j);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors_s]=SetupTargetMoments_FHorz(TMs,FnsToEvaluate,1,N_j,simoptions_s,Names_i);
wcs=zeros(1,N_j,npages,N_i+1); wcs(1,:,1,ilow)=1; wcs(1,3:N_j,1,ihigh)=1;
if isequal(size(selectors_s.ACStats.whichcombos),size(wcs))
    dd=max(abs(selectors_s.ACStats.whichcombos(:)-wcs(:)));
else
    dd=1;
end
fprintf('TargetParsing PType (d) agejshifter (high shifted by two ages): the high slot is on from age three on the common axis, whichcombos vs hand, this should be zero: %.3e \n',dd)

%% (e) the union 7-vectors
fprintf('TargetParsing PType (e) union whichstats vs hand ([1 1 1 1 1 1 1] for AllStats; [1 0 1 2 0 1 0] for AgeConditionalStats, quantiles being entry six): AllStats, this should be zero: %.3e; AgeConditionalStats, this should be zero: %.3e \n',max(abs(AllStats_whichstats(:)-[1;1;1;1;1;1;1])),max(abs(ACStats_whichstats(:)-[1;0;1;2;0;1;0])))

%% (f) two other field orders of the same targets
TM2=struct(); % restricted and per-type targets first
TM2.AllStats.rich.consumption.RatioMeanToMedian=1.05;
TM2.AllStats.rich.consumption.low.Minimum=0.3;
TM2.AllStats.poor.assets.StdDeviation=0.2;
TM2.AllStats.assets.high.QuantileCutoffs=(1:21)';
TM2.AllStats.assets.MoreInequality.Top1share=0.1;
TM2.AllStats.assets.Mean=1.1;
TM2.AllStats.earnings.low.LorenzCurve=linspace(0,1,100)';
TM2.AllStats.earnings.Gini=0.4;
TM2.AgeConditionalStats.poor.consumption.high.Variance=[1,1,1,1];
TM2.AgeConditionalStats.earnings.low.QuantileMeans=[ones(20,2),nan(20,2)];
TM2.AgeConditionalStats.rich.earnings.Mean=[NaN,NaN,1,1];
TM2.AgeConditionalStats.assets.high.LorenzCurve=[nan(100,1),ones(100,3)];
TM2.AgeConditionalStats.assets.Gini=[NaN,0.3,0.3,NaN];
TM2.AgeConditionalStats.assets.low.Mean=[NaN,2,3,4];
TM3=struct(); % restricted in the middle
TM3.AllStats.assets.Mean=1.1;
TM3.AllStats.assets.MoreInequality.Top1share=0.1;
TM3.AllStats.assets.high.QuantileCutoffs=(1:21)';
TM3.AllStats.rich.consumption.low.Minimum=0.3;
TM3.AllStats.rich.consumption.RatioMeanToMedian=1.05;
TM3.AllStats.earnings.Gini=0.4;
TM3.AllStats.earnings.low.LorenzCurve=linspace(0,1,100)';
TM3.AllStats.poor.assets.StdDeviation=0.2;
TM3.AgeConditionalStats.earnings.low.QuantileMeans=[ones(20,2),nan(20,2)];
TM3.AgeConditionalStats.assets.Gini=[NaN,0.3,0.3,NaN];
TM3.AgeConditionalStats.assets.low.Mean=[NaN,2,3,4];
TM3.AgeConditionalStats.assets.high.LorenzCurve=[nan(100,1),ones(100,3)];
TM3.AgeConditionalStats.poor.consumption.high.Variance=[1,1,1,1];
TM3.AgeConditionalStats.rich.earnings.Mean=[NaN,NaN,1,1];
orderbad=0; maxd=0;
for oo=1:2
    if oo==1
        TMo=TM2;
    else
        TMo=TM3;
    end
    [tv_o,~,~,~,~,~,AllStats_whichstats_o,FnsA_o,~,~,ACStats_whichstats_o,FnsL_o,~,~,selectors_o]=SetupTargetMoments_FHorz(TMo,FnsToEvaluate,1,N_j,simoptions,Names_i);
    orderbad=orderbad+~isequal(fieldnames(FnsA_o),FnA)+~isequal(fieldnames(FnsL_o),FnL)+~isequal(size(tv_o),size(targetmomentvec));
    if isequal(size(tv_o),size(targetmomentvec))
        s1=sort(tv_o(~isnan(tv_o))); s2=sort(targetmomentvec(~isnan(targetmomentvec)));
        if length(s1)==length(s2)
            maxd=max([maxd;abs(s1(:)-s2(:))]);
        else
            orderbad=orderbad+1;
        end
    end
    maxd=max([maxd;abs(AllStats_whichstats_o(:)-AllStats_whichstats(:));abs(ACStats_whichstats_o(:)-ACStats_whichstats(:))]);
    maxd=max([maxd;abs(selectors_o.AllStats.whichcombos(:)-selectors.AllStats.whichcombos(:));abs(selectors_o.AllStats.whichstats(:)-selectors.AllStats.whichstats(:))]);
    maxd=max([maxd;abs(selectors_o.ACStats.whichcombos(:)-selectors.ACStats.whichcombos(:));abs(selectors_o.ACStats.whichstats(:)-selectors.ACStats.whichstats(:))]);
end
fprintf('TargetParsing PType (f) two other field orders of the same targets: selectors, unions and sorted target vectors vs the first order, this should be zero: %.3e; function sets or sizes differ, this should be zero: %.3e \n',maxd,orderbad)

%% (g) the three-input legacy call (the estimation commands): unions and function sets, selectors empty
[tv_legacy,~,~,~,~,~,AllStats_whichstats_legacy,FnsA_legacy,~,~,ACStats_whichstats_legacy,FnsL_legacy,~,~,selectors_legacy]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,1);
if length(tv_legacy)==length(targetmomentvec)
    d=abs(tv_legacy-targetmomentvec); d=d(:); d=d(~isnan(d));
else
    d=1;
end
fprintf('TargetParsing PType (g) three-input call: target vector, unions and function sets vs the six-input call, this should be zero: %.3e; selectors not empty, this should be zero: %.3e \n',max([d;abs(AllStats_whichstats_legacy(:)-AllStats_whichstats(:));abs(ACStats_whichstats_legacy(:)-ACStats_whichstats(:))])+~isequal(fieldnames(FnsA_legacy),FnA)+~isequal(fieldnames(FnsL_legacy),FnL),~isempty(fieldnames(selectors_legacy)))

%% (h) invalid inputs must error
errored=0;
TMbad=TM; TMbad.AllStats.assets.mid.Mean=1; % unknown second-level name (neither a type nor a statistic)
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AllStats.assets.low=1; % a type with no statistic
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AllStats.wealth.low.Mean=1; % unknown function
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AgeConditionalStats.assets.high.Median=ones(1,N_j); % N_j entries but ngroups age groups for this type
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,1,N_j,simoptions,Names_i);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.conditionalrestrictions.assets=@(aprime,a,z) (a>1); % a name that is both a function and a restriction
try
    SetupTargetMoments_FHorz(TM,FnsToEvaluate,1,N_j,simoptions_bad,Names_i);
catch
    errored=errored+1;
end
fprintf('TargetParsing PType (h) five invalid inputs must all error, this should be zero: %.3e \n',5-errored)

output=struct();
output.targetmomentvec=targetmomentvec;
output.selectors=selectors;

end
