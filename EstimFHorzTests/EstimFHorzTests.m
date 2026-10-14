% Tests for EstimateLifeCycleModel_MethodOfMoments (the single-type FHorz method-of-moments estimator) and EstimateLifeCycleModel_MomentDerivatives.
%
% The model is the CoreFHorzTests nod_z_noe life-cycle model (ReturnFn_nod_z_noe_nosemiz) on small grids, as in CalibFHorzTests, with the
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
% Built 2026-10-07 as the third bank of the whichcombos programme (after CalibFHorzTests and CalibFHorzPTypeTests): the estimation
% commands had no tests. Parts 1 and 2 were run against the drivers as they stood first, to show their defects failing before the
% fixes (the omit-mask transpose, the Jacobian steps of a vector-valued parameter, the logmoments by-name and per-name forms, the
% sensitivity to calibrated parameters).
%
% Bootstrap standard errors and iterated GMM are not implemented in the toolkit (they error), and the banks check only that they do.
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): the standard-error machinery at the truth (skipestimation=1): J vs finite differences at the four step sizes, Sigma, std devs, confidence intervals, sensitivity matrix, identification, efficientW, logmoments forms, inputs that must error
% doPart(2): point estimates: start at the truth (exact) and away (close), lsqnonlin and fminsearch, a diagonal W (objective identity), constraints (J unaffected), omitestimparam (J of a vector-valued parameter)
% doPart(3): ExogShockFn and EiidShockFn with estimated shock parameters (J with the grids rebuilt), jequaloneDist as a function (both usergrids), CustomModelStats (both usergrids)
% doPart(4): the AutoCorrTransProbs, CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets (one restricted each): truth, J, whichcombos 1 vs 0 identical, a recovery
% doPart(5): sensitivity to calibrated parameters (Jorgensen) vs finite differences, Lambda*Jcalib identity; an unidentified pair (rank(J) short)
% doPart(6): cohorts entering at different ages: hand-built cohort moments, the vector and the auto-detected entry ages, one cohort equals the standard estimator, a recovery, must-error inputs, the last-age cutoff of a cohort's distribution (simoptions.lastDistAge / estimoptions.cohortlastage)
% doPart(7): EstimateLifeCycleModel_MomentDerivatives vs finite differences from full AllStats/LifeCycleProfiles runs; omit mask; ExogShockFn; the sorted output
% doPart(8): runtime of skipestimation with whichcombos=1 vs =0 (informational)
% doPart(9): parameters DERIVED by ParametrizeParamsFn feeding the shocks and the jequaloneDist function (built 2026-10-12 with the
%            fix it pins: calibrateshocks used to name-match only EstimParamNames, silently freezing shocks whose dependence runs
%            through the derived parameters; now ParseParametrizeParamsFn / estimoptions.ParametrizeParamsFn_OutputNames feed the
%            match, and a user-set estimoptions.calibrateshocks always wins)
%
% Parts are independent: the setup, the addpaths and the model preamble all sit OUTSIDE the if-blocks and so always run, and no
% part reads another part's output. Any subset can be run, in any combination. Anything added to this bank later must keep that true.
% doPart(10): cohort AutoCorr targets of the window and log rows (fnname_MperiodAverage, fnname_log, fnname_MperiodAverage_log): truth, the cutoff with M-1, refusals
% doPart(11): cohort AgeConditionalStats targets of the window and log rows (the 90th percentile of fnname_2periodAverage, fnname_MperiodAverage_log, fnname_log under a restriction): truth, the cutoff with M-1, refusals
doPart=[1,1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/EstimFHorzTestsdiary.txt','file')
    delete('./TestOutput/EstimFHorzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/EstimFHorzTestsdiary.txt
fprintf('EstimFHorzTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, CoreStructLeaves, shared by the Core banks
addpath('./EstimFHorzTests_subcodes/')
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/') % the return functions are shared with CoreFHorzTests
addpath('../CalibFHorzTests/CalibFHorzTests_subcodes/') % the shock, initial-distribution and custom-statistics functions are shared with CalibFHorzTests

%% The model: the CoreFHorzTests nod_z_noe life-cycle model on small grids
n_d=0;
d_grid=[];
n_a=101;
a_grid=5*linspace(0,1,n_a)'.^3;
n_z=5;
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
z_grid=exp(z_grid);
N_j=20;

Params.beta=0.95; % discount factor
DiscountFactorParamNames={'beta'};
Params.mewj=ones(1,N_j)/N_j;
AgeWeightParamNames={'mewj'};
Params.sigma=2; % CES utility param for consumption
Params.w=1;
Params.r=0.05;
Params.Jr=16;
Params.pension=0.5;
Params.agej=1:1:N_j;
Params.kappa_j=[0.5:0.1:1,ones(1,9),zeros(1,5)]; % deterministic earnings profile, zero in retirement

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

jequaloneDist=zeros(n_a,n_z,'gpuArray'); % everyone starts with no assets at the middle shock
jequaloneDist(1,ceil(n_z/2))=1;

%% ===== doPart(1): the standard-error machinery at the truth =====
if doPart(1)==1
    fprintf('\n===== doPart(1): standard errors at the truth (skipestimation=1): Jacobian, Sigma, confidence intervals, sensitivity, identification, efficientW, logmoments, errors =====\n')
    output=EstimFHorz_StdErrors(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(1)

%% ===== doPart(2): point estimates =====
if doPart(2)==1
    fprintf('\n===== doPart(2): point estimates: truth and away, lsqnonlin and fminsearch, a diagonal W, constraints, omitestimparam =====\n')
    output=EstimFHorz_PointEstimates(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(2)

%% ===== doPart(3): shocks with estimated parameters, jequaloneDist as a function, CustomModelStats =====
if doPart(3)==1
    fprintf('\n===== doPart(3): shocks with estimated parameters, jequaloneDist as a function, CustomModelStats =====\n')
    output=EstimFHorz_ShocksCustomStats(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(3)

%% ===== doPart(4): AutoCorr and cross-section targets: truth, J, whichcombos 1 vs 0, a recovery =====
if doPart(4)==1
    fprintf('\n===== doPart(4): AutoCorr and cross-section targets: truth, J, whichcombos 1 vs 0, a recovery =====\n')
    output=EstimFHorz_NewTargets(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(4)

%% ===== doPart(5): sensitivity to calibrated parameters and local identification =====
if doPart(5)==1
    fprintf('\n===== doPart(5): sensitivity to calibrated parameters and local identification =====\n')
    output=EstimFHorz_Sensitivity(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(5)

%% ===== doPart(6): cohorts entering at different ages =====
if doPart(6)==1
    fprintf('\n===== doPart(6): cohorts entering at different ages =====\n')
    output=EstimFHorz_Cohorts(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(6)

%% ===== doPart(7): EstimateLifeCycleModel_MomentDerivatives =====
if doPart(7)==1
    fprintf('\n===== doPart(7): EstimateLifeCycleModel_MomentDerivatives =====\n')
    output=EstimFHorz_MomentDerivatives(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(7)

%% ===== doPart(8): runtime, whichcombos=1 vs =0 (informational) =====
if doPart(8)==1
    fprintf('\n===== doPart(8): runtime, whichcombos=1 vs =0 (informational) =====\n')
    output=EstimFHorz_Runtime(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(8)

%% ===== doPart(9): ParametrizeParamsFn-derived parameters feeding the shocks and the jequaloneDist function =====
if doPart(9)==1
    fprintf('\n===== doPart(9): ParametrizeParamsFn-derived parameters feeding the shocks and the jequaloneDist function =====\n')
    output=EstimFHorz_ParametrizeParams(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(9)

%% ===== doPart(10): cohort AutoCorr targets of the window and log rows =====
if doPart(10)==1
    fprintf('\n===== doPart(10): cohort AutoCorr targets of the window and log rows =====\n')
    output=EstimFHorz_CohortsWindows(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(10)

%% ===== doPart(11): cohort AgeConditionalStats targets of the window and log rows =====
if doPart(11)==1
    fprintf('\n===== doPart(11): cohort AgeConditionalStats targets of the window and log rows =====\n')
    output=EstimFHorz_ACSCohortsWindows(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(11)

% All looks good!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/EstimFHorzTestsdiary.txt')

diary off
