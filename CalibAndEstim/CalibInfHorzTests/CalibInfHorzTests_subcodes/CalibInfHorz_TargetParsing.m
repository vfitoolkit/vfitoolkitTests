function output=CalibInfHorz_TargetParsing(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn)
% CalibInfHorzTests part 1: SetupTargetMoments_InfHorz, the target plumbing of CalibrateInfHorzAgentModel. No model is solved.
%
% TargetMoments with arbitrary values of the right sizes: AllStats, AutoCorrTransProbs and CrossSectionCovarCorr, unrestricted and
% under the two conditional restrictions rich and poor, every statistic kind (AllStats: Mean, RatioMeanToMedian, StdDeviation, Gini,
% LorenzCurve, Minimum, QuantileCutoffs, MoreInequality sub-statistic; AutoCorr: Mean, AutoCovariance, AutoCorrelation with and
% without a _kK horizon suffix; cross-section: CovarianceWith, CorrelationWith, Mean, StdDeviation, a CovarianceMatrix under a
% restriction with NaN entries), a NaN-omitted AutoCorr target.
% The call is SetupTargetMoments_InfHorz(TargetMoments,FnsToEvaluate,useptype,simoptions): with simoptions (its conditional
% restrictions) the selectors for caliboptions.whichcombos=1 are built, and returned as the 20th output with the AutoCorr horizons
% as the 21st. Checks (all exact):
%   (a) targetmomentvec equals the hand-built concatenation (TargetMoments field order, AllStats then AutoCorrTransProbs then
%       CrossSectionCovarCorr, matrices vectorised column-major), NaN pattern included; the flags and cumulative sizes
%   (b) the AllStats selectors: whichcombos [nFns,1+nRestr] and whichstats [nFns,1+nRestr,7] equal the hand-built arrays (a Gini
%       target overridden by a LorenzCurve one with a 1, the AllStats code for the Lorenz curve); the reduced FnsToEvaluate
%   (c) the AutoCorr selector: whichcombos [nFns,1+nRestr] (no age dimension in InfHorz), a NaN target leaving its combination off;
%       the horizons [2,3] read off the _kK suffixes; the reduced FnsToEvaluate
%   (d) the cross-section selector: pair-shaped whichcombos [nFns,nFns,1+nRestr], symmetric, the diagonal being a function's own
%       Mean/StdDeviation, a matrix target selecting its non-NaN entries; a matrix target makes the reduced set the full one
%   (e) the union 7-vector for AllStats (quantile targets switch on entry 6, not 5)
%   (f) two other field orders of the same targets give the same selectors, unions and function sets, and target vectors that are
%       permutations of each other
%   (g) the three-input legacy call (no simoptions) on the unrestricted targets gives the same target vector, unions and function
%       sets as the four-input call on them, and an empty selectors
%   (h) seven invalid inputs must error: an unknown kind name, an unknown AllStats statistic, a TransitionProbs AutoCorr target,
%       the self-correlation .(fn).CorrelationWith.(fn), a restricted AllStats target with no statistic level, a name that is both
%       a function and a restriction, a matrix target of the wrong size
% (n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, Params, DiscountFactorParamNames, ReturnFn are the bank's standard inputs, unused here.)

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w) (1+r)*a+w*z-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);

simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);
npages=3; % unrestricted, rich, poor

%% The targets (arbitrary values of the right sizes)
TM=struct();
TM.AllStats.earnings.Gini=0.4;
TM.AllStats.earnings.LorenzCurve=linspace(0,1,100)';
TM.AllStats.poor.assets.StdDeviation=0.2;
TM.AllStats.assets.Mean=1.1;
TM.AllStats.assets.QuantileCutoffs=(1:21)';
TM.AllStats.assets.MoreInequality.Top1share=0.1;
TM.AllStats.rich.consumption.Minimum=0.3;
TM.AllStats.rich.consumption.RatioMeanToMedian=1.05;
TM.AutoCorrTransProbs.assets.AutoCorrelation=0.9;
TM.AutoCorrTransProbs.assets.AutoCorrelation_k3=0.7;
TM.AutoCorrTransProbs.earnings.AutoCovariance=0.01;
TM.AutoCorrTransProbs.rich.consumption.AutoCorrelation_k2=0.5;
TM.AutoCorrTransProbs.poor.assets.Mean=0.3;
TM.AutoCorrTransProbs.poor.earnings.AutoCorrelation=NaN; % omitted: the only poor target of earnings, so that combination stays off
Mpoor=nan(4,4); Mpoor(1,1)=0.5; Mpoor(2,4)=0.1; % a matrix target under poor: the variance of assets and the earnings-richind covariance
TM.CrossSectionCovarCorr.assets.CovarianceWith.earnings=0.2;
TM.CrossSectionCovarCorr.earnings.CorrelationWith.consumption=0.5;
TM.CrossSectionCovarCorr.consumption.StdDeviation=0.3;
TM.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption=0.4;
TM.CrossSectionCovarCorr.poor.CovarianceMatrix=Mpoor;
TM.CrossSectionCovarCorr.richind.Mean=0.3;

hand=[0.4; linspace(0,1,100)'; 0.2; 1.1; (1:21)'; 0.1; 0.3; 1.05; ...   % AllStats, 127 entries
    0.9; 0.7; 0.01; 0.5; 0.3; NaN; ...                                  % AutoCorrTransProbs, 6 entries
    0.2; 0.5; 0.3; 0.4; Mpoor(:); 0.3];                                 % CrossSectionCovarCorr, 21 entries

ran=0;
try
    [targetmomentvec,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,autocorrmomentnames,autocorrcummomentsizes,~,FnsToEvaluate_AutoCorr,crosssecmomentnames,crossseccummomentsizes,~,FnsToEvaluate_CrossSec,cmsmomentnames,cmscummomentsizes,selectors,autocorrtimehorizons]=SetupTargetMoments_InfHorz(TM,FnsToEvaluate,0,simoptions);
    ran=1;
catch err
    fprintf('TargetParsing: SetupTargetMoments_InfHorz with simoptions errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('TargetParsing: SetupTargetMoments_InfHorz(TargetMoments,FnsToEvaluate,0,simoptions) with 21 outputs, did not run, this should be zero: %.3e \n',1-ran)

if ran==1
    %% (a) the target vector
    if length(targetmomentvec)==length(hand)
        d=abs(targetmomentvec-hand); d=d(:); d=d(~isnan(d));
        patternbad=sum(isnan(targetmomentvec)~=isnan(hand));
    else
        d=1; % (a length mismatch is a failure in itself)
        patternbad=abs(length(targetmomentvec)-length(hand));
    end
    fprintf('TargetParsing (a) targetmomentvec vs hand-built, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',max([d;0]),patternbad)
    fprintf('TargetParsing (a) flags usingallstats, usingautocorr, usingcrosssec, usingcustomstats = %i %i %i %i; cumulative sizes AllStats end %i (hand 127), AutoCorr end %i (hand 6), CrossSec end %i (hand 21), this should be zero: %.3e \n',usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatcummomentsizes(end),autocorrcummomentsizes(end),crossseccummomentsizes(end),abs(usingallstats-1)+abs(usingautocorr-1)+abs(usingcrosssec-1)+abs(usingcustomstats)+abs(allstatcummomentsizes(end)-127)+abs(autocorrcummomentsizes(end)-6)+abs(crossseccummomentsizes(end)-21)+abs(cmscummomentsizes)+(~isempty(cmsmomentnames{1})))
    fprintf('TargetParsing (a) names tables: AllStats %i rows (hand 8), AutoCorr %i rows (hand 6), CrossSec %i rows (hand 6), this should be zero: %.3e \n',size(allstatmomentnames,1),size(autocorrmomentnames,1),size(crosssecmomentnames,1),abs(size(allstatmomentnames,1)-8)+abs(size(autocorrmomentnames,1)-6)+abs(size(crosssecmomentnames,1)-6))

    %% (b) AllStats selectors: functions assets, earnings, consumption (FnsToEvaluate order; richind has no AllStats target)
    FnA=fieldnames(FnsToEvaluate_AllStats);
    wcA=zeros(3,npages); wsA=zeros(3,npages,7);
    wcA(1,1)=1; wsA(1,1,[1,6,7])=1;      % assets unrestricted: Mean, QuantileCutoffs (entry 6), MoreInequality (entry 7)
    wcA(1,3)=1; wsA(1,3,3)=1;            % assets under poor: StdDeviation
    wcA(2,1)=1; wsA(2,1,4)=1;            % earnings unrestricted: Gini (3) overridden by LorenzCurve (1)
    wcA(3,2)=1; wsA(3,2,[1,2,5])=1;      % consumption under rich: Minimum (5), RatioMeanToMedian (Mean and Median)
    fnok=isequal(FnA,{'assets';'earnings';'consumption'});
    fprintf('TargetParsing (b) AllStats reduced FnsToEvaluate are assets, earnings, consumption in that order, this should be zero: %.3e \n',1-fnok)
    if fnok && isfield(selectors,'AllStats')
        fprintf('TargetParsing (b) AllStats selector sizes [%s] and [%s], this should be zero: %.3e \n',num2str(size(selectors.AllStats.whichcombos)),num2str(size(selectors.AllStats.whichstats)),~isequal(size(selectors.AllStats.whichcombos),[3,npages])+~isequal(size(selectors.AllStats.whichstats),[3,npages,7]))
        if isequal(size(selectors.AllStats.whichcombos),[3,npages]) && isequal(size(selectors.AllStats.whichstats),[3,npages,7])
            fprintf('TargetParsing (b) AllStats whichcombos vs hand, this should be zero: %.3e; whichstats vs hand, this should be zero: %.3e \n',CoreMaxAbs(selectors.AllStats.whichcombos(:)-wcA(:)),CoreMaxAbs(selectors.AllStats.whichstats(:)-wsA(:)))
        end
    else
        fprintf('TargetParsing (b) selectors.AllStats missing or the function set wrong, this should be zero: %.3e \n',1)
    end

    %% (c) AutoCorr selector: functions assets, earnings, consumption; horizons 2 and 3
    FnC=fieldnames(FnsToEvaluate_AutoCorr);
    wcC=zeros(3,npages);
    wcC(1,1)=1;                          % assets unrestricted: AutoCorrelation and AutoCorrelation_k3
    wcC(1,3)=1;                          % assets under poor: Mean
    wcC(2,1)=1;                          % earnings unrestricted: AutoCovariance
    wcC(3,2)=1;                          % consumption under rich: AutoCorrelation_k2
                                         % earnings under poor: the only target is NaN, so off
    fnok=isequal(FnC,{'assets';'earnings';'consumption'});
    fprintf('TargetParsing (c) the horizons read off the _kK suffixes are [%s] (hand [2 3]) \n',num2str(autocorrtimehorizons(:)')) % (digits stay out of the check lines: CoreSummary reads the numbers before a wording as that check's)
    fprintf('TargetParsing (c) AutoCorr reduced FnsToEvaluate are assets, earnings, consumption in that order, this should be zero: %.3e; horizons vs hand, this should be zero: %.3e \n',1-fnok,~isequal(autocorrtimehorizons(:)',[2,3]))
    if fnok && isfield(selectors,'AutoCorr')
        fprintf('TargetParsing (c) AutoCorr selector size [%s] (hand [3 3]), this should be zero: %.3e \n',num2str(size(selectors.AutoCorr.whichcombos)),~isequal(size(selectors.AutoCorr.whichcombos),[3,npages]))
        if isequal(size(selectors.AutoCorr.whichcombos),[3,npages])
            fprintf('TargetParsing (c) AutoCorr whichcombos vs hand (the NaN target leaves earnings under poor off), this should be zero: %.3e \n',CoreMaxAbs(selectors.AutoCorr.whichcombos(:)-wcC(:)))
        end
    else
        fprintf('TargetParsing (c) selectors.AutoCorr missing or the function set wrong, this should be zero: %.3e \n',1)
    end

    %% (d) CrossSec selector: the matrix target makes the reduced set all four functions
    FnX=fieldnames(FnsToEvaluate_CrossSec);
    wcX=zeros(4,4,npages);
    wcX(1,2,1)=1; wcX(2,1,1)=1;          % assets-earnings covariance
    wcX(2,3,1)=1; wcX(3,2,1)=1;          % earnings-consumption correlation
    wcX(3,3,1)=1;                        % consumption StdDeviation (own stats on the diagonal)
    wcX(4,4,1)=1;                        % richind Mean
    wcX(1,3,2)=1; wcX(3,1,2)=1;          % under rich: assets-consumption correlation
    wcX(1,1,3)=1; wcX(2,4,3)=1; wcX(4,2,3)=1; % under poor: the non-NaN entries of the matrix, symmetrised
    fnok=isequal(FnX,{'assets';'earnings';'consumption';'richind'});
    fprintf('TargetParsing (d) CrossSec reduced FnsToEvaluate are all four (a matrix target), in order, this should be zero: %.3e \n',1-fnok)
    if fnok && isfield(selectors,'CrossSec')
        fprintf('TargetParsing (d) CrossSec selector size [%s] (hand [4 4 3]), this should be zero: %.3e \n',num2str(size(selectors.CrossSec.whichcombos)),~isequal(size(selectors.CrossSec.whichcombos),[4,4,npages]))
        if isequal(size(selectors.CrossSec.whichcombos),[4,4,npages])
            fprintf('TargetParsing (d) CrossSec whichcombos vs hand (symmetric, matrix entries selected where not NaN), this should be zero: %.3e \n',CoreMaxAbs(selectors.CrossSec.whichcombos(:)-wcX(:)))
        end
    else
        fprintf('TargetParsing (d) selectors.CrossSec missing or the function set wrong, this should be zero: %.3e \n',1)
    end

    %% (e) the union 7-vector for AllStats
    fprintf('TargetParsing (e) AllStats union whichstats [%s] vs hand [1 1 1 1 1 1 1] (quantiles being entry six), this should be zero: %.3e \n',num2str(AllStats_whichstats(:)'),CoreMaxAbs(AllStats_whichstats(:)-[1;1;1;1;1;1;1]))

    %% (f) two other field orders of the same targets
    TM2=struct(); % restricted targets first, kinds in another order
    TM2.CrossSectionCovarCorr.poor.CovarianceMatrix=Mpoor;
    TM2.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption=0.4;
    TM2.CrossSectionCovarCorr.richind.Mean=0.3;
    TM2.CrossSectionCovarCorr.consumption.StdDeviation=0.3;
    TM2.CrossSectionCovarCorr.earnings.CorrelationWith.consumption=0.5;
    TM2.CrossSectionCovarCorr.assets.CovarianceWith.earnings=0.2;
    TM2.AutoCorrTransProbs.poor.earnings.AutoCorrelation=NaN;
    TM2.AutoCorrTransProbs.poor.assets.Mean=0.3;
    TM2.AutoCorrTransProbs.rich.consumption.AutoCorrelation_k2=0.5;
    TM2.AutoCorrTransProbs.earnings.AutoCovariance=0.01;
    TM2.AutoCorrTransProbs.assets.AutoCorrelation_k3=0.7;
    TM2.AutoCorrTransProbs.assets.AutoCorrelation=0.9;
    TM2.AllStats.rich.consumption.RatioMeanToMedian=1.05;
    TM2.AllStats.rich.consumption.Minimum=0.3;
    TM2.AllStats.poor.assets.StdDeviation=0.2;
    TM2.AllStats.assets.MoreInequality.Top1share=0.1;
    TM2.AllStats.assets.QuantileCutoffs=(1:21)';
    TM2.AllStats.assets.Mean=1.1;
    TM2.AllStats.earnings.LorenzCurve=linspace(0,1,100)';
    TM2.AllStats.earnings.Gini=0.4;
    TM3=struct(); % a restricted target between two unrestricted ones, the last first-level name with a different second level (the leak the FHorz round found)
    TM3.AllStats.assets.Mean=1.1;
    TM3.AllStats.assets.QuantileCutoffs=(1:21)';
    TM3.AllStats.assets.MoreInequality.Top1share=0.1;
    TM3.AllStats.rich.consumption.Minimum=0.3;
    TM3.AllStats.rich.consumption.RatioMeanToMedian=1.05;
    TM3.AllStats.earnings.Gini=0.4;
    TM3.AllStats.earnings.LorenzCurve=linspace(0,1,100)';
    TM3.AllStats.poor.assets.StdDeviation=0.2;
    TM3.AutoCorrTransProbs.earnings.AutoCovariance=0.01;
    TM3.AutoCorrTransProbs.rich.consumption.AutoCorrelation_k2=0.5;
    TM3.AutoCorrTransProbs.assets.AutoCorrelation=0.9;
    TM3.AutoCorrTransProbs.assets.AutoCorrelation_k3=0.7;
    TM3.AutoCorrTransProbs.poor.assets.Mean=0.3;
    TM3.AutoCorrTransProbs.poor.earnings.AutoCorrelation=NaN;
    TM3.CrossSectionCovarCorr.richind.Mean=0.3;
    TM3.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption=0.4;
    TM3.CrossSectionCovarCorr.assets.CovarianceWith.earnings=0.2;
    TM3.CrossSectionCovarCorr.earnings.CorrelationWith.consumption=0.5;
    TM3.CrossSectionCovarCorr.consumption.StdDeviation=0.3;
    TM3.CrossSectionCovarCorr.poor.CovarianceMatrix=Mpoor;
    orderbad=0; maxd=0;
    for oo=1:2
        if oo==1
            TMo=TM2;
        else
            TMo=TM3;
        end
        try
            [tv_o,~,~,~,~,~,~,AllStats_whichstats_o,FnsA_o,~,~,~,FnsC_o,~,~,~,FnsX_o,~,~,selectors_o,horizons_o]=SetupTargetMoments_InfHorz(TMo,FnsToEvaluate,0,simoptions);
            orderbad=orderbad+~isequal(fieldnames(FnsA_o),FnA)+~isequal(fieldnames(FnsC_o),FnC)+~isequal(fieldnames(FnsX_o),FnX);
            orderbad=orderbad+~isequal(size(tv_o),size(targetmomentvec))+~isequal(horizons_o(:)',autocorrtimehorizons(:)');
            if isequal(size(tv_o),size(targetmomentvec))
                s1=sort(tv_o(~isnan(tv_o))); s2=sort(targetmomentvec(~isnan(targetmomentvec)));
                if length(s1)==length(s2)
                    maxd=max([maxd;abs(s1-s2)]);
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
            fprintf('TargetParsing (f) field order %i: SetupTargetMoments_InfHorz errored with: %s \n',oo,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
            orderbad=orderbad+1;
        end
    end
    fprintf('TargetParsing (f) two other field orders of the same targets: selectors, unions and sorted target vectors vs the first order, this should be zero: %.3e; function sets, sizes or horizons differ, this should be zero: %.3e \n',maxd,orderbad)
end % ran

%% (g) the three-input legacy call on the unrestricted targets equals the four-input call on them, with an empty selectors
TMu=struct();
TMu.AllStats.earnings.Gini=0.4;
TMu.AllStats.earnings.LorenzCurve=linspace(0,1,100)';
TMu.AllStats.assets.Mean=1.1;
TMu.AllStats.assets.QuantileCutoffs=(1:21)';
TMu.AllStats.assets.MoreInequality.Top1share=0.1;
TMu.AutoCorrTransProbs.assets.AutoCorrelation=0.9;
TMu.AutoCorrTransProbs.assets.AutoCorrelation_k3=0.7;
TMu.AutoCorrTransProbs.earnings.AutoCovariance=0.01;
TMu.CrossSectionCovarCorr.assets.CovarianceWith.earnings=0.2;
TMu.CrossSectionCovarCorr.earnings.CorrelationWith.consumption=0.5;
TMu.CrossSectionCovarCorr.consumption.StdDeviation=0.3;
TMu.CrossSectionCovarCorr.richind.Mean=0.3;
handu=[0.4; linspace(0,1,100)'; 1.1; (1:21)'; 0.1; 0.9; 0.7; 0.01; 0.2; 0.5; 0.3; 0.3];
ranl=0; d=1; setbad=1; selbad=1;
try
    [tv_l,~,~,~,~,~,~,AllStats_whichstats_l,FnsA_l,~,~,~,FnsC_l,~,~,~,FnsX_l]=SetupTargetMoments_InfHorz(TMu,FnsToEvaluate,0); % the legacy call: 19 outputs at most
    ranl=1;
    d=CoreMaxAbs(tv_l(:)-handu(:))*(length(tv_l)==length(handu))+abs(length(tv_l)-length(handu));
    setbad=~isequal(fieldnames(FnsA_l),{'assets';'earnings'})+~isequal(fieldnames(FnsC_l),{'assets';'earnings'})+~isequal(fieldnames(FnsX_l),{'assets';'earnings';'consumption';'richind'})+CoreMaxAbs(AllStats_whichstats_l(:)-[1;0;0;1;0;1;1]);
catch err
    fprintf('TargetParsing (g): the three-input call errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
try
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,selectors_l,horizons_l]=SetupTargetMoments_InfHorz(TMu,FnsToEvaluate,0);
    selbad=~isempty(fieldnames(selectors_l))+~isequal(horizons_l(:)',3);
catch err
    fprintf('TargetParsing (g): the three-input call with 21 outputs errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('TargetParsing (g) three-input legacy call on the unrestricted targets: did not run, this should be zero: %.3e; target vector vs hand, this should be zero: %.3e; function sets or AllStats union wrong, this should be zero: %.3e; selectors not empty or horizons wrong, this should be zero: %.3e \n',1-ranl,d,setbad,selbad)

%% (h) invalid inputs must error
errored=0;
TMbad=TM; TMbad.LifeCycleProfiles.assets.Mean=1; % unknown kind name
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,0,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AllStats.assets.Skewness=1; % unknown statistic
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,0,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AutoCorrTransProbs.assets.TransitionProbs=eye(3); % not targetable
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,0,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.assets.CorrelationWith.assets=1; % the self-correlation
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,0,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AllStats.rich.assets=0.3; % a restricted target with no statistic level
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,0,simoptions);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.conditionalrestrictions.assets=@(aprime,a,z) (a>1); % a name that is both a function and a restriction
try
    SetupTargetMoments_InfHorz(TM,FnsToEvaluate,0,simoptions_bad);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.CrossSectionCovarCorr.poor.CovarianceMatrix=nan(3,3); % four functions, so the matrix must be 4 by 4
try
    SetupTargetMoments_InfHorz(TMbad,FnsToEvaluate,0,simoptions);
catch
    errored=errored+1;
end
fprintf('TargetParsing (h) seven invalid inputs must all error, this should be zero: %.3e \n',7-errored)

output=struct();
output.ran=ran;
if ran==1
    output.targetmomentvec=targetmomentvec;
    output.selectors=selectors;
end

end
