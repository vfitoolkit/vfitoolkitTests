function fail_count=TestStatsFromWeightedGrid()
% Direct unit tests of StatsFromWeightedGrid(), which computes every AllStats and LifeCycleProfiles
% statistic in the toolkit. Everything here is hand-built distributions with known answers, no model.
% Runs on the GPU because that is how the toolkit always calls it, and prints runtimes at realistic
% sizes (the timing lines are not checks).
%
% Median: the median is the smallest value with cumulative mass >=0.5 (the quantile function at 0.5).
%   Formally m is a median if P(X<=m)>=0.5 and P(X>=m)>=0.5; when the cumulative mass hits exactly 0.5
%   at some value, everything from that value to the next is a median and the toolkit takes the lower.
%   U1 and U2 are distributions where the old rule (the point whose cumulative mass is NEAREST 0.5)
%   returned a value that is not a median.
% Gini: checked against the pairwise definition sum_i sum_j w_i*w_j*|y_i-y_j|/(2*mean), against the
%   unweighted closed form 2*sum(i*y_i)/(n*sum(y))-(n+1)/n (QuantEcon.py PR 937), and against the old
%   toolkit formula; plus analytical cases, invariance to rescaling the weights, and npoints=0.
% Top/bottom shares (whichstats(7)): analytical, and independent of npoints (they used to be read off the
%   Lorenz curve at round(0.99*npoints) etc, which crashed for npoints=0 and was wrong unless npoints was a multiple of 100).

fail_count=0;
TOL_EXACT=1e-10; % same bar as TestFnsToEvaluate (analytical-exact moments)

npoints=100;
nquantiles=20;
tolerance=10^(-12);
whichstats=ones(1,7); % AllStats default

%% U1. Median: cumulative masses 0.1, 0.49, 1.0. The old rule picked 2 (cumulative mass 0.49 is nearest 0.5), but P(X<=2)=0.49<0.5
y=gpuArray([1;2;3]); w=gpuArray([0.1;0.39;0.51]);
AllStats=StatsFromWeightedGrid(y,w,npoints,nquantiles,tolerance,0,whichstats);
err=abs(gather(AllStats.Median)-3);
fprintf('U1 Median of {1,2,3} w={0.1,0.39,0.51} minus 3, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
err=abs(gather(AllStats.MoreInequality.Percentile50th)-3);
fprintf('U1 Percentile50th minus 3, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U2. Median: cumulative masses 0.1, 0.4, 1.0. The median is 3 (the 'largest value with mass <=0.5' would be 2, which is not a median)
y=gpuArray([1;2;3]); w=gpuArray([0.1;0.3;0.6]);
AllStats=StatsFromWeightedGrid(y,w,npoints,nquantiles,tolerance,0,whichstats);
err=abs(gather(AllStats.Median)-3);
fprintf('U2 Median of {1,2,3} w={0.1,0.3,0.6} minus 3, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U3. Median at an exact tie: four values of mass 0.25 each (dyadic, so the cumulative mass is exactly 0.5 at 2). Both 2 and 3 are medians, toolkit takes the lower
y=gpuArray([1;2;3;4]); w=gpuArray(0.25*ones(4,1));
AllStats=StatsFromWeightedGrid(y,w,npoints,nquantiles,tolerance,0,whichstats);
err=abs(gather(AllStats.Median)-2);
fprintf('U3 Median at exact tie {1,2,3,4} equal mass minus 2, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U4. Median with unsorted input and a zero-weight point (presorted=0). Sorted: 2,3,4,5 with mass 0.25,0.25,0.3,0.2 (value 1 has zero mass), exact tie at 3
y=gpuArray([5;1;4;2;3]); w=gpuArray([0.2;0;0.3;0.25;0.25]);
AllStats=StatsFromWeightedGrid(y,w,npoints,nquantiles,tolerance,0,whichstats);
err=abs(gather(AllStats.Median)-3);
fprintf('U4 Median unsorted with zero-weight point minus 3, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
err=abs(gather(AllStats.Minimum)-2);
fprintf('U4 Minimum ignores the zero-weight point (minus 2), should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% Random distributions used by U5-U9 (fixed seed; lognormal values, random weights normalised to mass one)
rng(1,'twister')
n=2000;
y_cpu=exp(randn(n,1));
w_cpu=rand(n,1); w_cpu=w_cpu/sum(w_cpu);
y=gpuArray(y_cpu); w=gpuArray(w_cpu);
AllStats=StatsFromWeightedGrid(y,w,npoints,nquantiles,tolerance,0,whichstats);

%% U5. Median of the random distribution satisfies the definition, and is the SMALLEST value that does (no ties, so it is unique)
m=gather(AllStats.Median);
err=max(0,0.5-sum(w_cpu(y_cpu<=m)))+max(0,0.5-sum(w_cpu(y_cpu>=m)))+(sum(w_cpu(y_cpu<m))>=0.5)+(~any(y_cpu==m));
fprintf('U5 Median of random distribution satisfies the median definition, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U6. Gini vs the pairwise definition sum_i sum_j w_i*w_j*|y_i-y_j|/(2*mean)
Gini_pairwise=gather(sum(sum((w*w').*abs(y-y')))/(2*sum(w.*y)));
err=abs(gather(AllStats.Gini)-Gini_pairwise);
fprintf('U6 Gini vs pairwise definition (random weights), should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U7. Gini vs the old toolkit formula (trapezoid via sum of L_i*F_{i-1}-L_{i-1}*F_i), which it replaces
[SortedValues,sortindex]=sort(y); SortedWeights=w(sortindex);
CumSumSortedWeights=cumsum(SortedWeights);
CumSumWeightedSortedValues=cumsum(SortedValues.*SortedWeights);
CumSumWeightedSortedValues=CumSumWeightedSortedValues/CumSumWeightedSortedValues(end);
Gini_old=gather(sum(CumSumWeightedSortedValues(2:end).*CumSumSortedWeights(1:end-1)- CumSumWeightedSortedValues(1:end-1).*CumSumSortedWeights(2:end)));
err=abs(gather(AllStats.Gini)-Gini_old);
fprintf('U7 Gini vs old toolkit formula, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U8. Gini with equal weights vs the unweighted closed form 2*sum(i*y_i)/(n*sum(y))-(n+1)/n (QuantEcon.py PR 937)
AllStats_eq=StatsFromWeightedGrid(y,gpuArray(ones(n,1)/n),npoints,nquantiles,tolerance,0,whichstats);
ys=sort(y_cpu);
Gini_closedform=2*sum((1:n)'.*ys)/(n*sum(ys))-(n+1)/n;
err=abs(gather(AllStats_eq.Gini)-Gini_closedform);
fprintf('U8 Gini equal weights vs unweighted closed form, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U9. Gini is unchanged when the weights do not sum to one (rescaled by 0.7)
% Gini only (whichstats(4)=3): the percentiles and shares are defined relative to mass one, so are not meaningful for these weights
AllStats_scaled=StatsFromWeightedGrid(y,0.7*w,npoints,nquantiles,tolerance,0,[0,0,0,3,0,0,0]);
err=abs(gather(AllStats_scaled.Gini)-gather(AllStats.Gini));
fprintf('U9 Gini with weights rescaled by 0.7 minus Gini, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U10. presorted=1 (sorted, no zero weights) and presorted=2 (sorted, may contain zero weights) give the same Median and Gini as presorted=0
AllStats_ps1=StatsFromWeightedGrid(SortedValues,SortedWeights,npoints,nquantiles,tolerance,1,whichstats);
AllStats_ps2=StatsFromWeightedGrid([0;SortedValues],[0;SortedWeights],npoints,nquantiles,tolerance,2,whichstats);
err=abs(gather(AllStats_ps1.Gini)-gather(AllStats.Gini))+abs(gather(AllStats_ps1.Median)-gather(AllStats.Median));
fprintf('U10 presorted=1 vs presorted=0 (Gini+Median), should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
err=abs(gather(AllStats_ps2.Gini)-gather(AllStats.Gini))+abs(gather(AllStats_ps2.Median)-gather(AllStats.Median));
fprintf('U10 presorted=2 vs presorted=0 (Gini+Median), should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U11. npoints=0 still gives the Gini (it used to leave AllStats.Gini unset), and an empty Lorenz curve; for whichstats(4)=1, 2 and 3
% (whichstats(7) is off here, U16 tests it with npoints=0)
for w4=1:3
    whichstats_np0=[1,1,1,w4,1,1,0];
    AllStats_np0=StatsFromWeightedGrid(y,w,0,nquantiles,tolerance,0,whichstats_np0);
    if isfield(AllStats_np0,'Gini')
        err=abs(gather(AllStats_np0.Gini)-gather(AllStats.Gini));
    else
        err=1;
    end
    if w4<3
        err=err+~isfield(AllStats_np0,'LorenzCurve');
        if isfield(AllStats_np0,'LorenzCurve')
            err=err+~isempty(AllStats_np0.LorenzCurve);
        end
    end
    fprintf('U11 npoints=0, whichstats(4)=%i: Gini present and equal to npoints=100 Gini, should be zero: %.3e\n',w4,err); fail_count=fail_count+(err>TOL_EXACT);
end

%% U12. whichstats(4)=2 (vectorised Lorenz) and whichstats(4)=3 (Gini only) give the same Gini as whichstats(4)=1
AllStats_w42=StatsFromWeightedGrid(y,w,npoints,nquantiles,tolerance,0,[1,1,1,2,1,2,1]);
AllStats_w43=StatsFromWeightedGrid(y,w,npoints,nquantiles,tolerance,0,[1,1,1,3,1,1,0]);
err=abs(gather(AllStats_w42.Gini)-gather(AllStats.Gini))+abs(gather(AllStats_w43.Gini)-gather(AllStats.Gini));
fprintf('U12 Gini with whichstats(4)=2 and =3 vs =1, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U13. Analytical Gini: two-point {0,1} with P(1)=0.3 has Gini 1-0.3=0.7; discrete uniform on 1..20 has Gini (n-1)/(3n)=19/60
AllStats_2pt=StatsFromWeightedGrid(gpuArray([0;1]),gpuArray([0.7;0.3]),npoints,nquantiles,tolerance,0,whichstats);
err=abs(gather(AllStats_2pt.Gini)-0.7);
fprintf('U13 Gini of two-point {0,1} with P(1)=0.3 minus 0.7, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
AllStats_unif=StatsFromWeightedGrid(gpuArray((1:20)'),gpuArray(ones(20,1)/20),npoints,nquantiles,tolerance,0,whichstats);
err=abs(gather(AllStats_unif.Gini)-19/60);
fprintf('U13 Gini of discrete uniform 1..20 minus 19/60, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U14. Degenerate cases: a point mass has Gini 0; a negative value makes the Gini NaN (printed as 0/1 as the summary skips NaN)
AllStats_pm=StatsFromWeightedGrid(gpuArray(3*ones(5,1)),gpuArray(0.2*ones(5,1)),npoints,nquantiles,tolerance,0,whichstats);
err=abs(gather(AllStats_pm.Gini))+abs(gather(AllStats_pm.Median)-3);
fprintf('U14 Point mass: Gini and Median-3, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
AllStats_neg=StatsFromWeightedGrid(gpuArray([-1;2;3]),gpuArray([0.2;0.2;0.6]),npoints,nquantiles,tolerance,0,whichstats);
err=~isnan(gather(AllStats_neg.Gini));
fprintf('U14 Negative value: Gini is not NaN (0 means it is NaN), should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
err=abs(gather(AllStats_neg.Median)-3);
fprintf('U14 Negative value: Median minus 3, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);


%% U15. Analytical top/bottom shares: discrete uniform on 1..20 (total 210, mean 10.5).
% Bottom 50% is values 1..10 (55/210); top 10% is values 19,20 (39/210); top 5% is all value 20 (1/10.5); top 1% is 1% of mass at value 20 (0.2/10.5)
MI=AllStats_unif.MoreInequality;
err=abs(gather(MI.Bottom50share)-55/210)+abs(gather(MI.Top10share)-39/210)+abs(gather(MI.Top5share)-1/10.5)+abs(gather(MI.Top1share)-0.2/10.5);
fprintf('U15 Top/bottom shares of discrete uniform 1..20 vs analytical, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U16. Shares agree with the npoints=100 Lorenz curve (L(0.99), L(0.95), L(0.90), L(0.5) are its points 99, 95, 90, 50), for whichstats(4)=1 and 2
MI=AllStats.MoreInequality;
LC=gather(AllStats.LorenzCurve);
err=abs(gather(MI.Top1share)-(1-LC(99)))+abs(gather(MI.Top5share)-(1-LC(95)))+abs(gather(MI.Top10share)-(1-LC(90)))+abs(gather(MI.Bottom50share)-LC(50));
fprintf('U16 Shares vs npoints=100 Lorenz curve (whichstats(4)=1), should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
LC2=gather(AllStats_w42.LorenzCurve);
err=abs(gather(MI.Top1share)-(1-LC2(99)))+abs(gather(MI.Top5share)-(1-LC2(95)))+abs(gather(MI.Top10share)-(1-LC2(90)))+abs(gather(MI.Bottom50share)-LC2(50));
fprintf('U16 Shares vs npoints=100 Lorenz curve (whichstats(4)=2), should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% U17. Shares do not depend on npoints: npoints=0 (used to crash) and npoints=30 (used to give Top1share=0) match npoints=100
% Also whichstats(2)=0 and whichstats(4)=0, which whichstats(7) no longer needs
for np=[0,30]
    AllStats_np=StatsFromWeightedGrid(y,w,np,nquantiles,tolerance,0,[1,0,1,0,1,1,1]);
    MInp=AllStats_np.MoreInequality;
    err=abs(gather(MInp.Top1share)-gather(MI.Top1share))+abs(gather(MInp.Top5share)-gather(MI.Top5share))+abs(gather(MInp.Top10share)-gather(MI.Top10share))+abs(gather(MInp.Bottom50share)-gather(MI.Bottom50share))+abs(gather(MInp.Percentile50th)-gather(AllStats.Median));
    fprintf('U17 npoints=%i with whichstats(2)=0 and (4)=0: shares and Percentile50th vs npoints=100, should be zero: %.3e\n',np,err); fail_count=fail_count+(err>TOL_EXACT);
end


%% Runtimes on the GPU at realistic sizes (these are timings, not checks)
% For each size: the whole StatsFromWeightedGrid call with the AllStats default whichstats, and with the
% LifeCycleProfiles default [1,1,1,2,1,2,1]; then the median and Gini formulas alone, new vs old, on the
% already-sorted arrays (so this isolates what the rewrite changed).
fprintf('Runtimes of StatsFromWeightedGrid on the GPU (gputimeit, seconds): \n')
for N=[10^5, 10^6, 10^7]
    y=exp(randn(N,1,'gpuArray'));
    w=rand(N,1,'gpuArray'); w=w/sum(w);
    t_all=gputimeit(@() StatsFromWeightedGrid(y,w,npoints,nquantiles,tolerance,0,ones(1,7)));
    t_lc=gputimeit(@() StatsFromWeightedGrid(y,w,npoints,nquantiles,tolerance,0,[1,1,1,2,1,2,1]));
    t_gini=gputimeit(@() StatsFromWeightedGrid(y,w,npoints,nquantiles,tolerance,0,[0,0,0,3,0,0,0]));

    [SortedValues,sortindex]=sort(y); SortedWeights=w(sortindex);
    WeightedSortedValues=SortedValues.*SortedWeights;
    CumSumSortedWeights=cumsum(SortedWeights);
    CumSumSortedWeightedValues=cumsum(WeightedSortedValues);
    t_mediannew=gputimeit(@() SortedValues(find(CumSumSortedWeights>=0.5,1,'first')));
    t_medianold=gputimeit(@() min(abs(CumSumSortedWeights-0.5)),2);
    % The two Gini formulas are several statements each, so time them with tic/toc (wait(gpuDevice) makes the GPU finish before toc)
    nreps=20;
    wait(gpuDevice); tic;
    for rr=1:nreps
        Gini_new=sum(WeightedSortedValues.*(2*CumSumSortedWeights-SortedWeights-CumSumSortedWeights(end)))/(CumSumSortedWeightedValues(end)*CumSumSortedWeights(end));
    end
    wait(gpuDevice); t_gininew=toc/nreps;
    wait(gpuDevice); tic;
    for rr=1:nreps
        CumSumWeightedSortedValues=cumsum(WeightedSortedValues);
        CumSumWeightedSortedValues=CumSumWeightedSortedValues/CumSumWeightedSortedValues(end);
        Gini_old=sum(CumSumWeightedSortedValues(2:end).*CumSumSortedWeights(1:end-1)- CumSumWeightedSortedValues(1:end-1).*CumSumSortedWeights(2:end));
    end
    wait(gpuDevice); t_giniold=toc/nreps;
    fprintf('  N=%1.0e: whole call AllStats default %.4f, LifeCycle default %.4f, Gini only %.4f | median new %.5f old %.5f | Gini formula new %.5f old %.5f \n',N,t_all,t_lc,t_gini,t_mediannew,t_medianold,t_gininew,t_giniold)
end

end
