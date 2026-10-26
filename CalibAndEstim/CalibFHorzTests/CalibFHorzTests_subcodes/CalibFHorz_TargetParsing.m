function output=CalibFHorz_TargetParsing(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% CalibFHorzTests part 1: SetupTargetMoments_FHorz, the target plumbing of CalibrateLifeCycleModel. No model is solved.
%
% TargetMoments with arbitrary values of the right sizes: AllStats and AgeConditionalStats, unrestricted and under the two
% conditional restrictions rich and poor, every statistic kind (Mean, Median, RatioMeanToMedian, Variance/StdDeviation, Gini,
% LorenzCurve, Minimum/Maximum, QuantileCutoffs/QuantileMeans, MoreInequality sub-statistic), NaN-omitted ages, matrix-valued
% age-conditional targets (LorenzCurve, QuantileMeans: age groups as columns), and simoptions.agegroupings.
% Checks (all exact):
%   (a) targetmomentvec equals the hand-built concatenation (TargetMoments field order, values vectorised column-major)
%   (b) the AllStats selectors: whichcombos [nFns,1+nRestr] and whichstats [nFns,1+nRestr,7] equal the hand-built arrays
%   (c) the AgeConditionalStats selectors: whichcombos [nFns,ngroups,1+nRestr] and whichstats [..,7] equal the hand-built arrays,
%       with an age off where the target is NaN and a LorenzCurve target overriding a Gini-only 3 with a 2
%   (d) the union 7-vectors (quantile targets switch on entry 6, not 5)
%   (e) the reduced FnsToEvaluate hold every targeted function, in FnsToEvaluate order, for three field orders of the same targets;
%       the selectors and unions are identical across the orders, and the target vectors are permutations of each other
%   (f) the three-input legacy call (no N_j, no simoptions) gives the same unions and function sets and an empty selectors
%   (g) five invalid inputs must error: an unknown first-level name, an unknown statistic, an age-conditional target of the wrong
%       length, a restricted target with no statistic level, a name that is both a function and a restriction
% (n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, Params, DiscountFactorParamNames, AgeWeightParamNames, ReturnFn, jequaloneDist are
% the bank's standard inputs and unused here.)

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) (1+r)*a+w*kappa_j*z*(agej<Jr)+pension*(agej>=Jr)-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);

simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);
simoptions.agegroupings=1:5:N_j;
ngroups=length(simoptions.agegroupings);
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
TM.AgeConditionalStats.assets.Mean=[NaN,2,3,4];
TM.AgeConditionalStats.assets.Gini=[NaN,0.3,0.3,NaN];
TM.AgeConditionalStats.assets.LorenzCurve=[nan(100,1),ones(100,3)];
TM.AgeConditionalStats.rich.earnings.Mean=[NaN,NaN,1,1];
TM.AgeConditionalStats.earnings.QuantileMeans=[ones(20,2),nan(20,2)];
TM.AgeConditionalStats.poor.consumption.Variance=[1,1,1,1];

[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,0,N_j,simoptions);

%% (a) the target vector
hand=[0.4; linspace(0,1,100)'; 0.2; 1.1; (1:21)'; 0.1; 0.3; 1.05; ...
    [NaN;2;3;4]; [NaN;0.3;0.3;NaN]; reshape([nan(100,1),ones(100,3)],[],1); [NaN;NaN;1;1]; reshape([ones(20,2),nan(20,2)],[],1); [1;1;1;1]];
if length(targetmomentvec)==length(hand)
    d=abs(targetmomentvec-hand); d=d(:); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(hand));
else
    d=1; % (a length mismatch is a failure in itself)
    patternbad=abs(length(targetmomentvec)-length(hand));
end
fprintf('TargetParsing (a) targetmomentvec vs hand-built, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',max([d;0]),patternbad)
fprintf('TargetParsing (a) flags usingallstats, usinglcp, usingcustomstats = %i %i %i; cumulative sizes AllStats end %i (hand 127), AgeConditionalStats end %i (hand 496), this should be zero: %.3e \n',usingallstats,usinglcp,usingcustomstats,allstatcummomentsizes(end),acscummomentsizes(end),abs(usingallstats-1)+abs(usinglcp-1)+abs(usingcustomstats)+abs(allstatcummomentsizes(end)-127)+abs(acscummomentsizes(end)-496)+abs(cmscummomentsizes)+(~isempty(cmsmomentnames{1})))

%% (b) AllStats selectors: functions assets, earnings, consumption (FnsToEvaluate order; richind is not targeted)
FnA=fieldnames(FnsToEvaluate_AllStats);
wcA=zeros(3,npages); wsA=zeros(3,npages,7);
wcA(1,1)=1; wsA(1,1,[1,6,7])=1;      % assets unrestricted: Mean, QuantileCutoffs, MoreInequality
wcA(1,3)=1; wsA(1,3,3)=1;            % assets under poor: StdDeviation
wcA(2,1)=1; wsA(2,1,4)=1;            % earnings unrestricted: Gini (3) overridden by LorenzCurve (1)
wcA(3,2)=1; wsA(3,2,[1,2,5])=1;      % consumption under rich: Minimum, RatioMeanToMedian (Mean and Median)
fnok=isequal(FnA,{'assets';'earnings';'consumption'});
fprintf('TargetParsing (b) AllStats reduced FnsToEvaluate are assets, earnings, consumption in that order, this should be zero: %.3e \n',1-fnok)
if fnok
    fprintf('TargetParsing (b) AllStats whichcombos vs hand, this should be zero: %.3e; whichstats vs hand, this should be zero: %.3e \n',max(abs(selectors.AllStats.whichcombos(:)-wcA(:))),max(abs(selectors.AllStats.whichstats(:)-wsA(:))))
    fprintf('TargetParsing (b) AllStats selector sizes [%s] and [%s], this should be zero: %.3e \n',num2str(size(selectors.AllStats.whichcombos)),num2str(size(selectors.AllStats.whichstats)),~isequal(size(selectors.AllStats.whichcombos),[3,npages])+~isequal(size(selectors.AllStats.whichstats),[3,npages,7]))
end

%% (c) AgeConditionalStats selectors
FnL=fieldnames(FnsToEvaluate_ACStats);
wcL=zeros(3,ngroups,npages); wsL=zeros(3,ngroups,npages,7);
wsL(1,2:4,1,1)=1;                    % assets Mean at ages 2-4
wsL(1,2:4,1,4)=2;                    % assets Gini (3 at ages 2,3) overridden by the LorenzCurve matrix (2 at ages 2-4)
wsL(2,3:4,2,1)=1;                    % earnings under rich: Mean at ages 3,4
wsL(2,1:2,1,6)=1;                    % earnings QuantileMeans at ages 1,2
wsL(3,:,3,3)=1;                      % consumption under poor: Variance at every age
wcL(1,2:4,1)=1; wcL(2,3:4,2)=1; wcL(2,1:2,1)=1; wcL(3,:,3)=1;
fnok=isequal(FnL,{'assets';'earnings';'consumption'});
fprintf('TargetParsing (c) AgeConditionalStats reduced FnsToEvaluate are assets, earnings, consumption in that order, this should be zero: %.3e \n',1-fnok)
if fnok
    fprintf('TargetParsing (c) AgeConditionalStats whichcombos vs hand, this should be zero: %.3e; whichstats vs hand, this should be zero: %.3e \n',max(abs(selectors.ACStats.whichcombos(:)-wcL(:))),max(abs(selectors.ACStats.whichstats(:)-wsL(:))))
    fprintf('TargetParsing (c) AgeConditionalStats selector sizes [%s] and [%s], this should be zero: %.3e \n',num2str(size(selectors.ACStats.whichcombos)),num2str(size(selectors.ACStats.whichstats)),~isequal(size(selectors.ACStats.whichcombos),[3,ngroups,npages])+~isequal(size(selectors.ACStats.whichstats),[3,ngroups,npages,7]))
end

%% (d) the union 7-vectors
fprintf('TargetParsing (d) union whichstats vs hand ([1 1 1 1 1 1 1] for AllStats; [1 0 1 2 0 1 0] for AgeConditionalStats, quantiles being entry six): AllStats, this should be zero: %.3e; AgeConditionalStats, this should be zero: %.3e \n',max(abs(AllStats_whichstats(:)-[1;1;1;1;1;1;1])),max(abs(ACStats_whichstats(:)-[1;0;1;2;0;1;0])))

%% (e) three field orders of the same targets
TM2=struct(); % restricted targets first
TM2.AllStats.rich.consumption.RatioMeanToMedian=1.05;
TM2.AllStats.rich.consumption.Minimum=0.3;
TM2.AllStats.poor.assets.StdDeviation=0.2;
TM2.AllStats.assets.MoreInequality.Top1share=0.1;
TM2.AllStats.assets.QuantileCutoffs=(1:21)';
TM2.AllStats.assets.Mean=1.1;
TM2.AllStats.earnings.LorenzCurve=linspace(0,1,100)';
TM2.AllStats.earnings.Gini=0.4;
TM2.AgeConditionalStats.poor.consumption.Variance=[1,1,1,1];
TM2.AgeConditionalStats.rich.earnings.Mean=[NaN,NaN,1,1];
TM2.AgeConditionalStats.earnings.QuantileMeans=[ones(20,2),nan(20,2)];
TM2.AgeConditionalStats.assets.LorenzCurve=[nan(100,1),ones(100,3)];
TM2.AgeConditionalStats.assets.Gini=[NaN,0.3,0.3,NaN];
TM2.AgeConditionalStats.assets.Mean=[NaN,2,3,4];
TM3=struct(); % a restricted target between two unrestricted ones (the last first-level name has a different second level)
TM3.AllStats.assets.Mean=1.1;
TM3.AllStats.assets.QuantileCutoffs=(1:21)';
TM3.AllStats.assets.MoreInequality.Top1share=0.1;
TM3.AllStats.rich.consumption.Minimum=0.3;
TM3.AllStats.rich.consumption.RatioMeanToMedian=1.05;
TM3.AllStats.earnings.Gini=0.4;
TM3.AllStats.earnings.LorenzCurve=linspace(0,1,100)';
TM3.AllStats.poor.assets.StdDeviation=0.2;
TM3.AgeConditionalStats.earnings.QuantileMeans=[ones(20,2),nan(20,2)];
TM3.AgeConditionalStats.assets.Mean=[NaN,2,3,4];
TM3.AgeConditionalStats.assets.Gini=[NaN,0.3,0.3,NaN];
TM3.AgeConditionalStats.assets.LorenzCurve=[nan(100,1),ones(100,3)];
TM3.AgeConditionalStats.poor.consumption.Variance=[1,1,1,1];
TM3.AgeConditionalStats.rich.earnings.Mean=[NaN,NaN,1,1];
orderbad=0; maxd=0;
for oo=1:2
    if oo==1
        TMo=TM2;
    else
        TMo=TM3;
    end
    [tv_o,~,~,~,~,~,AllStats_whichstats_o,FnsA_o,~,~,ACStats_whichstats_o,FnsL_o,~,~,selectors_o]=SetupTargetMoments_FHorz(TMo,FnsToEvaluate,0,N_j,simoptions);
    orderbad=orderbad+~isequal(fieldnames(FnsA_o),FnA)+~isequal(fieldnames(FnsL_o),FnL);
    orderbad=orderbad+~isequal(size(tv_o),size(targetmomentvec));
    if isequal(size(tv_o),size(targetmomentvec))
        s1=sort(tv_o(~isnan(tv_o))); s2=sort(targetmomentvec(~isnan(targetmomentvec)));
        if length(s1)==length(s2)
            maxd=max([maxd;abs(s1-s2)]);
        else
            orderbad=orderbad+1;
        end
    end
    maxd=max([maxd;abs(AllStats_whichstats_o(:)-AllStats_whichstats(:));abs(ACStats_whichstats_o(:)-ACStats_whichstats(:))]);
    maxd=max([maxd;abs(selectors_o.AllStats.whichcombos(:)-selectors.AllStats.whichcombos(:));abs(selectors_o.AllStats.whichstats(:)-selectors.AllStats.whichstats(:))]);
    maxd=max([maxd;abs(selectors_o.ACStats.whichcombos(:)-selectors.ACStats.whichcombos(:));abs(selectors_o.ACStats.whichstats(:)-selectors.ACStats.whichstats(:))]);
end
fprintf('TargetParsing (e) two other field orders of the same targets: selectors, unions and sorted target vectors vs the first order, this should be zero: %.3e; function sets or sizes differ, this should be zero: %.3e \n',maxd,orderbad)

%% (f) the three-input legacy call (the estimation commands): unions and function sets, selectors empty
% (without simoptions a first-level name that is not a function is taken to be a restriction, as before)
[tv_legacy,~,~,~,~,~,AllStats_whichstats_legacy,FnsA_legacy,~,~,ACStats_whichstats_legacy,FnsL_legacy,~,~,selectors_legacy]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,0);
d=abs(tv_legacy-targetmomentvec); d=d(:); d=d(~isnan(d));
fprintf('TargetParsing (f) three-input call: target vector, unions and function sets vs the five-input call, this should be zero: %.3e; selectors not empty, this should be zero: %.3e \n',max([d;abs(AllStats_whichstats_legacy(:)-AllStats_whichstats(:));abs(ACStats_whichstats_legacy(:)-ACStats_whichstats(:))])+~isequal(fieldnames(FnsA_legacy),FnA)+~isequal(fieldnames(FnsL_legacy),FnL),~isempty(fieldnames(selectors_legacy)))

%% (g) invalid inputs must error
errored=0;
TMbad=TM; TMbad.AllStats.wealth.Mean=1; % unknown first-level name
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AllStats.assets.Skewness=1; % unknown statistic
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AgeConditionalStats.assets.Median=ones(1,N_j); % N_j entries but ngroups age groups
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
TMbad=TM; TMbad.AllStats.rich.assets=0.3; % a restricted target with no statistic level
try
    SetupTargetMoments_FHorz(TMbad,FnsToEvaluate,0,N_j,simoptions);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.conditionalrestrictions.assets=@(aprime,a,z) (a>1); % a name that is both a function and a restriction
try
    SetupTargetMoments_FHorz(TM,FnsToEvaluate,0,N_j,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('TargetParsing (g) five invalid inputs must all error, this should be zero: %.3e \n',5-errored)

output=struct();
output.targetmomentvec=targetmomentvec;
output.selectors=selectors;
output.allstatmomentnames=allstatmomentnames;
output.acsmomentnames=acsmomentnames;

end
