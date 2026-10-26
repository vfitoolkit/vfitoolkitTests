function output=CalibInfHorzPType_TargetParsing(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn)
% CalibInfHorzPTypeTests part 1: SetupTargetMoments_InfHorz with permanent types (useptype=1). No model is solved.
%
% TargetMoments with arbitrary values of the right sizes, in every layout: grouped fn.stat, per-type fn.type.stat, restricted
% restriction.fn.stat and restriction.fn.type.stat, fn.MoreInequality.substat, AutoCorr names with a _kK horizon suffix, cross-section
% pairs grouped and per type (fn1.type.CovarianceWith.fn2), a grouped CovarianceMatrix under a restriction with NaN entries, and a
% wholly NaN per-type AutoCorr target (which must warn and select nothing).
% The call is SetupTargetMoments_InfHorz(TargetMoments,FnsToEvaluate,1,simoptions,Names_i): the selectors (20th output) carry a
% trailing type dimension of N_i+1 (types in Names_i order, then the grouped slot), the horizons are the 21st output. The names tables
% have four columns (AllStats, AutoCorr) and five (CrossSec). Checks (all exact):
%   (a) targetmomentvec equals the hand-built concatenation (field order, kinds in the order AllStats, AutoCorrTransProbs,
%       CrossSectionCovarCorr), NaN pattern included; the flags, cumulative sizes and name-table shapes
%   (b) AllStats selectors [nFns,1+nRestr,N_i+1] and [..,7] equal the hand-built arrays; the reduced FnsToEvaluate
%   (c) the AutoCorr selector [nFns,1+nRestr,N_i+1], the wholly NaN target leaving its slot off; the horizons; the reduced set
%   (d) the cross-section selector [nFns,nFns,1+nRestr,N_i+1], symmetric, a matrix target selecting its non-NaN entries in the
%       grouped slot and making the reduced set the full one
%   (e) the union 7-vector for AllStats (quantile targets switch on entry 6)
%   (f) two other field orders of the same targets give the same selectors, unions and function sets, and a permutation of the
%       target vector
%   (g) the three-input legacy call (no simoptions, no Names_i) on the unrestricted targets (grouped and per type, a type being any
%       second-level name that is not a statistic) gives the same target vector, unions and function sets, and an empty selectors
%   (h) eight invalid inputs must error: an unknown second-level name (neither a type nor a statistic), a type with no statistic,
%       an unknown function, a TransitionProbs AutoCorr target, the per-type self-correlation, a per-type matrix
%       (CovarianceMatrix_ptype), a matrix target of the wrong size, a name that is both a function and a restriction
% (the model inputs are the bank's standard ones and unused here)

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w_pt) (1+r)*a+w_pt*z-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);
simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);
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
TM.AutoCorrTransProbs.assets.AutoCorrelation=0.9;                % grouped
TM.AutoCorrTransProbs.assets.low.AutoCorrelation_k3=0.7;         % type low, horizon 3
TM.AutoCorrTransProbs.earnings.AutoCovariance=0.01;              % grouped
TM.AutoCorrTransProbs.rich.consumption.high.AutoCorrelation_k2=0.5; % restricted, type high, horizon 2
TM.AutoCorrTransProbs.poor.assets.Mean=0.3;                      % restricted, grouped
TM.AutoCorrTransProbs.poor.earnings.low.AutoCorrelation=NaN;     % omitted in full: warns, selects nothing
Mpoor=nan(4,4); Mpoor(1,1)=0.5; Mpoor(2,4)=0.1;
TM.CrossSectionCovarCorr.assets.CovarianceWith.earnings=0.2;     % grouped pair
TM.CrossSectionCovarCorr.earnings.high.CorrelationWith.consumption=0.5; % type high pair
TM.CrossSectionCovarCorr.consumption.low.StdDeviation=0.3;       % type low own stat
TM.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption=0.4;  % restricted, grouped pair
TM.CrossSectionCovarCorr.poor.CovarianceMatrix=Mpoor;            % restricted, grouped matrix
TM.CrossSectionCovarCorr.richind.Mean=0.3;                       % grouped own stat

hand=[0.4; linspace(0,1,100)'; 0.2; 1.1; (1:21)'; 0.1; 0.3; 1.05; ...   % AllStats, 127 entries
    0.9; 0.7; 0.01; 0.5; 0.3; NaN; ...                                  % AutoCorrTransProbs, 6 entries
    0.2; 0.5; 0.3; 0.4; Mpoor(:); 0.3];                                 % CrossSectionCovarCorr, 21 entries

ran=0; wmsg='';
lastwarn('');
try
    [targetmomentvec,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,autocorrmomentnames,autocorrcummomentsizes,~,FnsToEvaluate_AutoCorr,crosssecmomentnames,crossseccummomentsizes,~,FnsToEvaluate_CrossSec,cmsmomentnames,cmscummomentsizes,selectors,autocorrtimehorizons]=SetupTargetMoments_InfHorz(TM,FnsToEvaluate,1,simoptions,Names_i);
    ran=1; wmsg=lastwarn;
catch err
    fprintf('TargetParsing PType: SetupTargetMoments_InfHorz with simoptions and Names_i errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('TargetParsing PType: SetupTargetMoments_InfHorz(TargetMoments,FnsToEvaluate,1,simoptions,Names_i) with 21 outputs, did not run, this should be zero: %.3e; the wholly NaN target warns (the warning names it), did not, this should be zero: %.3e \n',1-ran,~contains(wmsg,'poor.earnings.low.AutoCorrelation is NaN in full'))

if ran==1
    %% (a) the target vector
    if length(targetmomentvec)==length(hand)
        d=abs(targetmomentvec-hand); d=d(:); d=d(~isnan(d));
        patternbad=sum(isnan(targetmomentvec)~=isnan(hand));
    else
        d=1; patternbad=abs(length(targetmomentvec)-length(hand));
    end
    fprintf('TargetParsing PType (a) targetmomentvec vs hand-built, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',max([d;0]),patternbad)
    fprintf('TargetParsing PType (a) flags usingallstats, usingautocorr, usingcrosssec, usingcustomstats = %i %i %i %i; cumulative sizes AllStats end %i (hand 127), AutoCorr end %i (hand 6), CrossSec end %i (hand 21), this should be zero: %.3e \n',usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatcummomentsizes(end),autocorrcummomentsizes(end),crossseccummomentsizes(end),abs(usingallstats-1)+abs(usingautocorr-1)+abs(usingcrosssec-1)+abs(usingcustomstats)+abs(allstatcummomentsizes(end)-127)+abs(autocorrcummomentsizes(end)-6)+abs(crossseccummomentsizes(end)-21)+abs(cmscummomentsizes)+(~isempty(cmsmomentnames{1})))
    fprintf('TargetParsing PType (a) names tables: AllStats %i rows by %i columns (hand 8 by 4), AutoCorr %i by %i (hand 6 by 4), CrossSec %i by %i (hand 6 by 5), this should be zero: %.3e \n',size(allstatmomentnames,1),size(allstatmomentnames,2),size(autocorrmomentnames,1),size(autocorrmomentnames,2),size(crosssecmomentnames,1),size(crosssecmomentnames,2),abs(size(allstatmomentnames,1)-8)+abs(size(allstatmomentnames,2)-4)+abs(size(autocorrmomentnames,1)-6)+abs(size(autocorrmomentnames,2)-4)+abs(size(crosssecmomentnames,1)-6)+abs(size(crosssecmomentnames,2)-5))

    %% (b) AllStats selectors: functions assets, earnings, consumption (FnsToEvaluate order; richind has no AllStats target)
    FnA=fieldnames(FnsToEvaluate_AllStats);
    wcA=zeros(3,npages,N_i+1); wsA=zeros(3,npages,N_i+1,7);
    wsA(1,1,igrp,[1,7])=1;   % assets grouped unrestricted: Mean, MoreInequality
    wsA(1,1,ihigh,6)=1;      % assets type high unrestricted: QuantileCutoffs (entry 6)
    wsA(1,3,igrp,3)=1;       % assets grouped under poor: StdDeviation
    wsA(2,1,igrp,4)=3;       % earnings grouped unrestricted: Gini only
    wsA(2,1,ilow,4)=1;       % earnings type low unrestricted: LorenzCurve
    wsA(3,2,ilow,5)=1;       % consumption type low under rich: Minimum
    wsA(3,2,igrp,[1,2])=1;   % consumption grouped under rich: RatioMeanToMedian
    wcA(1,1,igrp)=1; wcA(1,1,ihigh)=1; wcA(1,3,igrp)=1; wcA(2,1,igrp)=1; wcA(2,1,ilow)=1; wcA(3,2,ilow)=1; wcA(3,2,igrp)=1;
    fnok=isequal(FnA,{'assets';'earnings';'consumption'});
    fprintf('TargetParsing PType (b) AllStats reduced FnsToEvaluate are assets, earnings, consumption in that order, this should be zero: %.3e \n',1-fnok)
    if fnok && isfield(selectors,'AllStats')
        sizebad=~isequal(size(selectors.AllStats.whichcombos),[3,npages,N_i+1])+~isequal(size(selectors.AllStats.whichstats),[3,npages,N_i+1,7]);
        fprintf('TargetParsing PType (b) AllStats selector sizes [%s] and [%s], this should be zero: %.3e \n',num2str(size(selectors.AllStats.whichcombos)),num2str(size(selectors.AllStats.whichstats)),sizebad)
        if sizebad==0
            fprintf('TargetParsing PType (b) AllStats whichcombos with the type dimension vs hand, this should be zero: %.3e; whichstats vs hand, this should be zero: %.3e \n',CoreMaxAbs(selectors.AllStats.whichcombos(:)-wcA(:)),CoreMaxAbs(selectors.AllStats.whichstats(:)-wsA(:)))
        end
    else
        fprintf('TargetParsing PType (b) selectors.AllStats missing or the function set wrong, this should be zero: %.3e \n',1)
    end

    %% (c) AutoCorr selector
    FnC=fieldnames(FnsToEvaluate_AutoCorr);
    wcC=zeros(3,npages,N_i+1);
    wcC(1,1,igrp)=1;         % assets grouped: AutoCorrelation
    wcC(1,1,ilow)=1;         % assets low: AutoCorrelation_k3
    wcC(2,1,igrp)=1;         % earnings grouped: AutoCovariance
    wcC(3,2,ihigh)=1;        % consumption high under rich: AutoCorrelation_k2
    wcC(1,3,igrp)=1;         % assets grouped under poor: Mean
                             % earnings low under poor: NaN in full, so off
    fnok=isequal(FnC,{'assets';'earnings';'consumption'});
    fprintf('TargetParsing PType (c) the horizons read off the _kK suffixes are [%s] (hand [2 3]) \n',num2str(autocorrtimehorizons(:)'))
    fprintf('TargetParsing PType (c) AutoCorr reduced FnsToEvaluate are assets, earnings, consumption in that order, this should be zero: %.3e; horizons vs hand, this should be zero: %.3e \n',1-fnok,~isequal(autocorrtimehorizons(:)',[2,3]))
    if fnok && isfield(selectors,'AutoCorr')
        sizebad=~isequal(size(selectors.AutoCorr.whichcombos),[3,npages,N_i+1]);
        fprintf('TargetParsing PType (c) AutoCorr selector size [%s], this should be zero: %.3e \n',num2str(size(selectors.AutoCorr.whichcombos)),sizebad)
        if sizebad==0
            fprintf('TargetParsing PType (c) AutoCorr whichcombos with the type dimension vs hand (the wholly NaN target leaves its slot off), this should be zero: %.3e \n',CoreMaxAbs(selectors.AutoCorr.whichcombos(:)-wcC(:)))
        end
    else
        fprintf('TargetParsing PType (c) selectors.AutoCorr missing or the function set wrong, this should be zero: %.3e \n',1)
    end

    %% (d) CrossSec selector: the matrix target makes the reduced set all four functions
    FnX=fieldnames(FnsToEvaluate_CrossSec);
    wcX=zeros(4,4,npages,N_i+1);
    wcX(1,2,1,igrp)=1; wcX(2,1,1,igrp)=1;      % assets-earnings covariance, grouped
    wcX(2,3,1,ihigh)=1; wcX(3,2,1,ihigh)=1;    % earnings-consumption correlation, type high
    wcX(3,3,1,ilow)=1;                         % consumption StdDeviation, type low (own stats on the diagonal)
    wcX(4,4,1,igrp)=1;                         % richind Mean, grouped
    wcX(1,3,2,igrp)=1; wcX(3,1,2,igrp)=1;      % under rich: assets-consumption correlation, grouped
    wcX(1,1,3,igrp)=1; wcX(2,4,3,igrp)=1; wcX(4,2,3,igrp)=1; % under poor: the non-NaN entries of the grouped matrix, symmetrised
    fnok=isequal(FnX,{'assets';'earnings';'consumption';'richind'});
    fprintf('TargetParsing PType (d) CrossSec reduced FnsToEvaluate are all four (a matrix target), in order, this should be zero: %.3e \n',1-fnok)
    if fnok && isfield(selectors,'CrossSec')
        sizebad=~isequal(size(selectors.CrossSec.whichcombos),[4,4,npages,N_i+1]);
        fprintf('TargetParsing PType (d) CrossSec selector size [%s], this should be zero: %.3e \n',num2str(size(selectors.CrossSec.whichcombos)),sizebad)
        if sizebad==0
            fprintf('TargetParsing PType (d) CrossSec whichcombos with the type dimension vs hand (symmetric, matrix entries selected where not NaN), this should be zero: %.3e \n',CoreMaxAbs(selectors.CrossSec.whichcombos(:)-wcX(:)))
        end
    else
        fprintf('TargetParsing PType (d) selectors.CrossSec missing or the function set wrong, this should be zero: %.3e \n',1)
    end

    %% (e) the union 7-vector for AllStats
    fprintf('TargetParsing PType (e) AllStats union whichstats [%s] vs hand [1 1 1 1 1 1 1] (quantiles being entry six), this should be zero: %.3e \n',num2str(AllStats_whichstats(:)'),CoreMaxAbs(AllStats_whichstats(:)-[1;1;1;1;1;1;1]))

    %% (f) two other field orders of the same targets
    TM2=struct(); % kinds in another order, restricted and per-type targets first
    TM2.CrossSectionCovarCorr.poor.CovarianceMatrix=Mpoor;
    TM2.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption=0.4;
    TM2.CrossSectionCovarCorr.richind.Mean=0.3;
    TM2.CrossSectionCovarCorr.consumption.low.StdDeviation=0.3;
    TM2.CrossSectionCovarCorr.earnings.high.CorrelationWith.consumption=0.5;
    TM2.CrossSectionCovarCorr.assets.CovarianceWith.earnings=0.2;
    TM2.AutoCorrTransProbs.poor.earnings.low.AutoCorrelation=NaN;
    TM2.AutoCorrTransProbs.poor.assets.Mean=0.3;
    TM2.AutoCorrTransProbs.rich.consumption.high.AutoCorrelation_k2=0.5;
    TM2.AutoCorrTransProbs.earnings.AutoCovariance=0.01;
    TM2.AutoCorrTransProbs.assets.low.AutoCorrelation_k3=0.7;
    TM2.AutoCorrTransProbs.assets.AutoCorrelation=0.9;
    TM2.AllStats.rich.consumption.RatioMeanToMedian=1.05;
    TM2.AllStats.rich.consumption.low.Minimum=0.3;
    TM2.AllStats.poor.assets.StdDeviation=0.2;
    TM2.AllStats.assets.high.QuantileCutoffs=(1:21)';
    TM2.AllStats.assets.MoreInequality.Top1share=0.1;
    TM2.AllStats.assets.Mean=1.1;
    TM2.AllStats.earnings.low.LorenzCurve=linspace(0,1,100)';
    TM2.AllStats.earnings.Gini=0.4;
    TM3=struct(); % restricted in the middle
    TM3.AllStats.assets.Mean=1.1;
    TM3.AllStats.assets.MoreInequality.Top1share=0.1;
    TM3.AllStats.assets.high.QuantileCutoffs=(1:21)';
    TM3.AllStats.rich.consumption.low.Minimum=0.3;
    TM3.AllStats.rich.consumption.RatioMeanToMedian=1.05;
    TM3.AllStats.earnings.Gini=0.4;
    TM3.AllStats.earnings.low.LorenzCurve=linspace(0,1,100)';
    TM3.AllStats.poor.assets.StdDeviation=0.2;
    TM3.AutoCorrTransProbs.earnings.AutoCovariance=0.01;
    TM3.AutoCorrTransProbs.rich.consumption.high.AutoCorrelation_k2=0.5;
    TM3.AutoCorrTransProbs.assets.AutoCorrelation=0.9;
    TM3.AutoCorrTransProbs.assets.low.AutoCorrelation_k3=0.7;
    TM3.AutoCorrTransProbs.poor.assets.Mean=0.3;
    TM3.AutoCorrTransProbs.poor.earnings.low.AutoCorrelation=NaN;
    TM3.CrossSectionCovarCorr.richind.Mean=0.3;
    TM3.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption=0.4;
    TM3.CrossSectionCovarCorr.assets.CovarianceWith.earnings=0.2;
    TM3.CrossSectionCovarCorr.earnings.high.CorrelationWith.consumption=0.5;
    TM3.CrossSectionCovarCorr.consumption.low.StdDeviation=0.3;
    TM3.CrossSectionCovarCorr.poor.CovarianceMatrix=Mpoor;
    orderbad=0; maxd=0;
    for oo=1:2
        if oo==1
            TMo=TM2;
        else
            TMo=TM3;
        end
        try
            [tv_o,~,~,~,~,~,~,AllStats_whichstats_o,FnsA_o,~,~,~,FnsC_o,~,~,~,FnsX_o,~,~,selectors_o,horizons_o]=SetupTargetMoments_InfHorz(TMo,FnsToEvaluate,1,simoptions,Names_i);
            orderbad=orderbad+~isequal(fieldnames(FnsA_o),FnA)+~isequal(fieldnames(FnsC_o),FnC)+~isequal(fieldnames(FnsX_o),FnX);
            orderbad=orderbad+~isequal(size(tv_o),size(targetmomentvec))+~isequal(horizons_o(:)',autocorrtimehorizons(:)');
            if isequal(size(tv_o),size(targetmomentvec))
                s1=sort(tv_o(~isnan(tv_o))); s2=sort(targetmomentvec(~isnan(targetmomentvec)));
                if length(s1)==length(s2)
                    maxd=max([maxd;abs(s1(:)-s2(:))]);
                else
                    orderbad=orderbad+1;
                end
            end
            maxd=max([maxd;abs(AllStats_whichstats_o(:)-AllStats_whichstats(:))]);
            if isequal(size(selectors_o.AllStats.whichcombos),size(selectors.AllStats.whichcombos)) && isequal(size(selectors_o.AutoCorr.whichcombos),size(selectors.AutoCorr.whichcombos)) && isequal(size(selectors_o.CrossSec.whichcombos),size(selectors.CrossSec.whichcombos))
                maxd=max([maxd;abs(selectors_o.AllStats.whichcombos(:)-selectors.AllStats.whichcombos(:));abs(selectors_o.AllStats.whichstats(:)-selectors.AllStats.whichstats(:))]);
                maxd=max([maxd;abs(selectors_o.AutoCorr.whichcombos(:)-selectors.AutoCorr.whichcombos(:));abs(selectors_o.CrossSec.whichcombos(:)-selectors.CrossSec.whichcombos(:))]);
            else
                orderbad=orderbad+1;
            end
        catch err
            fprintf('TargetParsing PType (f) field order %i: SetupTargetMoments_InfHorz errored with: %s \n',oo,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
            orderbad=orderbad+1;
        end
    end
    fprintf('TargetParsing PType (f) two other field orders of the same targets: selectors, unions and sorted target vectors vs the first order, this should be zero: %.3e; function sets, sizes or horizons differ, this should be zero: %.3e \n',maxd,orderbad)
end % ran

%% (g) the three-input legacy call on the unrestricted targets (grouped and per type) equals the five-input call on them, with an empty selectors
TMu=struct();
TMu.AllStats.earnings.Gini=0.4;
TMu.AllStats.earnings.low.LorenzCurve=linspace(0,1,100)';
TMu.AllStats.assets.Mean=1.1;
TMu.AllStats.assets.high.QuantileCutoffs=(1:21)';
TMu.AllStats.assets.MoreInequality.Top1share=0.1;
TMu.AutoCorrTransProbs.assets.AutoCorrelation=0.9;
TMu.AutoCorrTransProbs.assets.low.AutoCorrelation_k3=0.7;
TMu.AutoCorrTransProbs.earnings.AutoCovariance=0.01;
TMu.CrossSectionCovarCorr.assets.CovarianceWith.earnings=0.2;
TMu.CrossSectionCovarCorr.earnings.high.CorrelationWith.consumption=0.5;
TMu.CrossSectionCovarCorr.consumption.low.StdDeviation=0.3;
TMu.CrossSectionCovarCorr.richind.Mean=0.3;
handu=[0.4; linspace(0,1,100)'; 1.1; (1:21)'; 0.1; 0.9; 0.7; 0.01; 0.2; 0.5; 0.3; 0.3];
ranl=0; d=1; setbad=1; selbad=1;
try
    [tv_l,~,~,~,~,~,~,AllStats_whichstats_l,FnsA_l,~,~,~,FnsC_l,~,~,~,FnsX_l]=SetupTargetMoments_InfHorz(TMu,FnsToEvaluate,1); % the legacy call: at most 19 outputs
    ranl=1;
    d=CoreMaxAbs(tv_l(:)-handu(:))*(length(tv_l)==length(handu))+abs(length(tv_l)-length(handu));
    setbad=~isequal(fieldnames(FnsA_l),{'assets';'earnings'})+~isequal(fieldnames(FnsC_l),{'assets';'earnings'})+~isequal(fieldnames(FnsX_l),{'assets';'earnings';'consumption';'richind'})+CoreMaxAbs(AllStats_whichstats_l(:)-[1;0;0;1;0;1;1]);
catch err
    fprintf('TargetParsing PType (g): the three-input call errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
try
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors_l,horizons_l]=SetupTargetMoments_InfHorz(TMu,FnsToEvaluate,1);
    selbad=~isempty(fieldnames(selectors_l))+~isequal(horizons_l(:)',3);
catch err
    fprintf('TargetParsing PType (g): the three-input call with 21 outputs errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('TargetParsing PType (g) three-input legacy call on the unrestricted targets, grouped and per type: did not run, this should be zero: %.3e; target vector vs hand, this should be zero: %.3e; function sets or AllStats union wrong, this should be zero: %.3e; selectors not empty or horizons wrong, this should be zero: %.3e \n',1-ranl,d,setbad,selbad)

%% (h) invalid inputs must error
errored=0;
TMbad=TM; TMbad.AllStats.assets.mid.Mean=1; % unknown second-level name (neither a type nor a statistic)
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,1,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AllStats.assets.low=1; % a type with no statistic
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,1,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AllStats.wealth.low.Mean=1; % unknown function
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,1,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AutoCorrTransProbs.assets.low.TransitionProbs=eye(3); % not targetable
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,1,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.earnings.low.CorrelationWith.earnings=1; % the per-type self-correlation
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,1,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.CovarianceMatrix_ptype.low=ones(4); % the per-type matrices cannot be targeted
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,1,simoptions,Names_i);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.poor.CovarianceMatrix=nan(3,3); % four functions, so the matrix must be 4 by 4
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,1,simoptions,Names_i);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.conditionalrestrictions.assets=@(aprime,a,z) (a>1); % a name that is both a function and a restriction
try
    SetupTargetMoments_InfHorz(TM,FnsToEvaluate,1,simoptions_bad,Names_i);
catch
    errored=errored+1;
end
fprintf('TargetParsing PType (h) eight invalid inputs must all error, this should be zero: %.3e \n',8-errored)

output=struct();
output.ran=ran;
if ran==1
    output.targetmomentvec=targetmomentvec;
    output.selectors=selectors;
end

end
