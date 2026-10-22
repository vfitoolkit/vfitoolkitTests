% Tests for EstimateLifeCycleModel_PType_MethodOfMoments (the FHorz method-of-moments estimator with permanent types). (EstimateLifeCycleModel_MomentDerivatives
% has no permanent-type version; it is tested in EstimFHorzTests part 7.)
%
% The model is the CoreFHorzPTypeTests part-22 model (two permanent types, low and high, with a per-type earnings profile kappa_j_pt,
% type weights 0.7 and 0.3 and geometric age weights) on small grids, as in CalibFHorzPTypeTests, with the
% grid interpolation layer on so that the moments vary smoothly with the parameters (finite differences at a relative step of 1.5e-4
% are then meaningful). Targets are generated from the model at known parameter values, so the truth is known exactly: an estimation
% started at the truth returns the truth with a zero objective, and the standard-error machinery can be checked at the truth with
% estimoptions.skipestimation=1 without running an optimiser.
%
% Checks are identities ("this should be zero") wherever two routes compute the same object: the Jacobian against finite differences
% the test takes itself from direct objective calls at the reported up/down parameter vectors, Sigma, the standard deviations, the
% confidence intervals and the sensitivity matrix recomputed from the reported J, W and Omega, the objective against
% (M_d-M_m)'W(M_d-M_m) at the estimate. Only the optimiser recoveries from a wrong start are "close to zero".
%
% Built 2026-10-07 as the fourth bank of the whichcombos programme (after CalibFHorzTests, CalibFHorzPTypeTests and EstimFHorzPTypeTests):
% the estimation commands had no tests. Parts 1 and 2 were run against the drivers as they stood first, to show their defects failing before the
% fixes (the omit-mask transpose, the Jacobian steps of a vector-valued parameter, the logmoments by-name and per-name forms, the
% sensitivity to calibrated parameters).
%
% Bootstrap standard errors and iterated GMM are not implemented in the toolkit (they error), and the banks check only that they do.
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): the standard-error machinery at the truth (skipestimation=1): J vs finite differences at the four step sizes, Sigma, std devs, confidence intervals, sensitivity matrix, identification, efficientW, logmoments forms, inputs that must error
% doPart(2): point estimates: start at the truth (exact) and away (close), lsqnonlin and fminsearch, a diagonal W (objective identity), constraints (J unaffected), omitestimparam on the per-type kappa_j_pt (J of a vector-valued per-type parameter), a per-type pension as a structure
% doPart(3): ExogShockFn common and per type, EiidShockFn, with estimated shock parameters (J with the grids rebuilt), jequaloneDist as a function (both usergrids), CustomModelStats (both usergrids)
% doPart(4): the AutoCorrTransProbs, CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets, grouped and per type, one restricted: truth, J, whichcombos 1 vs 0 identical, a recovery
% doPart(5): sensitivity to calibrated parameters (Jorgensen) vs finite differences, Lambda*Jcalib identity; an unidentified pair (rank(J) short)
% doPart(6): cohorts entering at different ages, by type: hand-built cohort moments (grouped and per type), the vector and the auto-detected entry ages, one cohort equals the standard estimator, a recovery, must-error inputs, the last-age cutoff of a cohort's distribution (simoptions.lastDistAge / estimoptions.cohortlastage), a QuantileCutoffs (matrix) target, whichcombos=0 vs the default =1 identical (both timed)
% doPart(7): runtime of skipestimation with whichcombos=1 vs =0 (informational)
% doPart(8): cohort AutoCorr targets of the window and log rows, grouped and per type (fnname_MperiodAverage, fnname_log, fnname_MperiodAverage_log): truth, the cutoff with M-1, a refusal, whichcombos=0 vs the default =1 identical (both timed)
% doPart(9): cohort AgeConditionalStats targets of the window and log rows, grouped and per type (the 90th percentile of fnname_2periodAverage, fnname_MperiodAverage_log, fnname_log): truth, the cutoff with M-1, a refusal, whichcombos=0 vs the default =1 identical (both timed)
% doPart(10): cohort CustomModelStats targets with permanent types (the same function as without cohorts, with the grouped DFJ full-cell autocorrelation): truth, the N_j cutoff, a refusal
% doPart(11): CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets of the windows with permanent types (grouped and per type, one restricted, an age-conditional window by start age): truth, J, whichcombos 1 vs 0 identical, a recovery
%
% Parts are independent: the setup, the addpaths and the model preamble all sit OUTSIDE the if-blocks and so always run, and no
% part reads another part's output. Any subset can be run, in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/EstimFHorzPTypeTestsdiary.txt','file')
    delete('./TestOutput/EstimFHorzPTypeTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/EstimFHorzPTypeTestsdiary.txt
fprintf('EstimFHorzPTypeTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, CoreStructLeaves, shared by the Core banks
addpath('./EstimFHorzPTypeTests_subcodes/')
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/') % the return functions are shared with CoreFHorzTests
addpath('../CalibFHorzPTypeTests/CalibFHorzPTypeTests_subcodes/') % the shock, initial-distribution and custom-statistics functions are shared with CalibFHorzPTypeTests
addpath('../CoreFHorzTests/CoreFHorzTests_Setup/') % CoreFHorzSetup_SemiExoStateFn, for a later semi-exogenous part

%% The model: the CoreFHorzPTypeTests part-22 model on small grids (shared with withMaxLikelihood/)
EstimFHorzPType_setup

%% ===== doPart(1): the standard-error machinery at the truth =====
if doPart(1)==1
    fprintf('\n===== doPart(1): standard errors at the truth (skipestimation=1): Jacobian, Sigma, confidence intervals, sensitivity, identification, efficientW, logmoments, errors =====\n')
    output=EstimFHorzPType_StdErrors(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(1)

%% ===== doPart(2): point estimates =====
if doPart(2)==1
    fprintf('\n===== doPart(2): point estimates: truth and away, lsqnonlin and fminsearch, a diagonal W, constraints, omitestimparam =====\n')
    output=EstimFHorzPType_PointEstimates(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(2)

%% ===== doPart(3): shocks with estimated parameters (common and per type), jequaloneDist as a function, CustomModelStats =====
if doPart(3)==1
    fprintf('\n===== doPart(3): shocks with estimated parameters (common and per type), jequaloneDist as a function, CustomModelStats =====\n')
    output=EstimFHorzPType_ShocksCustomStats(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(3)

%% ===== doPart(4): AutoCorr and cross-section targets, grouped and per type: truth, J, whichcombos 1 vs 0, a recovery =====
if doPart(4)==1
    fprintf('\n===== doPart(4): AutoCorr and cross-section targets, grouped and per type: truth, J, whichcombos 1 vs 0, a recovery =====\n')
    output=EstimFHorzPType_NewTargets(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(4)

%% ===== doPart(5): sensitivity to calibrated parameters and local identification =====
if doPart(5)==1
    fprintf('\n===== doPart(5): sensitivity to calibrated parameters and local identification =====\n')
    output=EstimFHorzPType_Sensitivity(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(5)

%% ===== doPart(6): cohorts entering at different ages, by type =====
if doPart(6)==1
    fprintf('\n===== doPart(6): cohorts entering at different ages, by type =====\n')
    output=EstimFHorzPType_Cohorts(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(6)

%% ===== doPart(7): runtime, whichcombos=1 vs =0 (informational) =====
if doPart(7)==1
    fprintf('\n===== doPart(7): runtime, whichcombos=1 vs =0 (informational) =====\n')
    output=EstimFHorzPType_Runtime(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(7)

%% ===== doPart(8): cohort AutoCorr targets of the window and log rows =====
if doPart(8)==1
    fprintf('\n===== doPart(8): cohort AutoCorr targets of the window and log rows, grouped and per type =====\n')
    output=EstimFHorzPType_CohortsWindows(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(8)

%% ===== doPart(9): cohort AgeConditionalStats targets of the window and log rows =====
if doPart(9)==1
    fprintf('\n===== doPart(9): cohort AgeConditionalStats targets of the window and log rows, grouped and per type =====\n')
    output=EstimFHorzPType_ACSCohortsWindows(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(9)

%% ===== doPart(10): cohort CustomModelStats targets =====
if doPart(10)==1
    fprintf('\n===== doPart(10): cohort CustomModelStats targets, with permanent types =====\n')
    output=EstimFHorzPType_CohortsCustomStats(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(10)

%% ===== doPart(11): CrossSec window targets, grouped and per type: truth, J, whichcombos 1 vs 0, a recovery =====
if doPart(11)==1
    fprintf('\n===== doPart(11): CrossSec window targets, grouped and per type: truth, J, whichcombos 1 vs 0, a recovery =====\n')
    output=EstimFHorzPType_CSWindows(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(11)

% All looks good!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/EstimFHorzPTypeTestsdiary.txt')

diary off
