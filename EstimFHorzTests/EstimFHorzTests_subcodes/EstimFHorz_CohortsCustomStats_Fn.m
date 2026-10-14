function CustomStats=EstimFHorz_CohortsCustomStats_Fn(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,estimoptions,vfoptions,simoptions)
% EstimFHorzTests part 12: the CustomModelStats function of the cohort estimation (the same signature as without cohorts).
%   entryage   the cohort's entry age, read off simoptions.jequaloneDistAge (checks the estimator tells the function which cohort)
%   hasV       1 when V is given (the cohort objective keeps V when CustomModelStats is used)
%   incmean    the mean of inc by age (LifeCycleProfiles), 1 x N_j
%   dfjac      the autocorrelation of the log of the 2-period average of inc at horizon 2 under the restriction rich, normalised as
%              De Nardi, French and Jones (2010): [AutoCovariance + (PairMean_j-Mean_j)(PairMean_jplusk-Mean_j+2)]/(StdDeviation_j*StdDeviation_j+2),
%              the pair cross-product about the full-cell means over the full-cell standard deviations, 1 x N_j-2
CustomStats.entryage=simoptions.jequaloneDistAge;
CustomStats.hasV=double(~isempty(V));
FnsInc.inc=FnsToEvaluate.inc;
simoptions_l=simoptions; simoptions_l.whichstats=[1,0,0,0,0,0,0];
ACS=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsInc,Parameters,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_l);
CustomStats.incmean=gather(ACS.inc.Mean);
simoptions_a=simoptions; simoptions_a.ntimeperiods=2; simoptions_a.timehorizons=2; simoptions_a.logwhichcombos=1;
ACP=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsInc,Parameters,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,simoptions_a);
R=ACP.rich.inc_2periodAverage_log;
m=gather(R.Mean); s=gather(R.StdDeviation);
num=gather(R.AutoCovariance_horizon2)+(gather(R.PairMean_j_horizon2)-m(1:N_j-2)).*(gather(R.PairMean_jplusk_horizon2)-m(3:N_j));
CustomStats.dfjac=num./(s(1:N_j-2).*s(3:N_j));
end
