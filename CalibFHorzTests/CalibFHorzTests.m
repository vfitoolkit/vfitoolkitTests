% Tests for CalibrateLifeCycleModel (the single-type FHorz calibration command) and its target plumbing.
%
% The model is the CoreFHorzTests nod_z_noe life-cycle model (ReturnFn_nod_z_noe_nosemiz, from CoreFHorzTests/CoreFHorz_ReturnFns)
% on small grids, so that one calibration takes seconds. Targets are generated from the model at known parameter values, so the
% truth is known exactly: a calibration started at the truth must return the truth with an objective of zero, and one started
% away from it must get back within the tolerances.
%
% Checks are identities ("this should be zero") wherever two routes compute the same object: the target vector against a
% hand-built one, the objective's moment vector against a hand-built one from full LifeCycleProfiles/AllStats runs, the objective
% forms against their formulas, a calibration started at the truth. Only the optimiser recoveries from a wrong start are "close
% to zero" (held to CoreSummary's looser bar), because they are only as exact as the tolerances.
%
% What the bank is for, 2026-10-06: the wiring of simoptions.whichcombos (and a per-combination simoptions.whichstats) into the
% calibration, through caliboptions.whichcombos (=1 default: the stats commands compute only the targeted combinations; =0: every
% statistic of every targeted function). Every moment must be identical either way. The bank also pins four pre-existing defects
% fixed in the same programme: quantile targets switched on Minimum/Maximum instead of the quantiles; the reduced FnsToEvaluate
% depended on the field order of TargetMoments; caliboptions.logmoments by name read the AgeConditionalStats names from the AllStats
% branch; and the by-target-name vector form of logmoments had an off-by-one.
%
% Added 2026-10-07 (part 8): three more kinds of target standing alongside AllStats and AgeConditionalStats, TargetMoments.AutoCorrTransProbs
% (the autocorrelation command, horizons read off the _horizonK suffixes of the target names, 1 for none), .CrossSectionCovarCorr and
% .AgeConditionalCrossSectionCovarCorr (the two cross-sectional covariance/correlation commands, matrix targets with NaN accepted). The
% target vector is AllStats, AgeConditionalStats, AutoCorrTransProbs, CrossSectionCovarCorr, AgeConditionalCrossSectionCovarCorr,
% CustomModelStats. caliboptions.whichcombos is now always 1 or 0 (the objective takes no other value).
%
% Not covered (later programmes): CalibrateLifeCycleModel_PType, the estimation commands (MethodOfMoments, MomentDerivatives,
% cohorts), the OLG/GE calibration commands, the panel-simulation moments, CMA-ES (fminalgo=4).
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): SetupTargetMoments_FHorz: target vector, selectors (whichcombos/whichstats), reduced FnsToEvaluate, field-order independence, invalid inputs
% doPart(2): the objective's moment vector vs a hand-built one from full LifeCycleProfiles/AllStats runs; whichcombos 1 and 0
% doPart(3): objective forms (sum_squared, sum_logratiosquared, MethodOfMoments, residual vector, logmoments, weights) vs their formulas; logmoments plumbing through CalibrateLifeCycleModel
% doPart(4): round trips through CalibrateLifeCycleModel: start at the truth (exact) and away from it (close), lsqnonlin and fminsearch, parameter constraints, omitcalibparam, ParametrizeParamsFn
% doPart(5): shock grids that depend on a calibrated parameter (ExogShockFn), jequaloneDist as a function, CustomModelStats under both usergrids settings
% doPart(6): runtime of one objective evaluation, whichcombos=1 vs =0 (informational)
% doPart(7): an iid shock through EiidShockFn with its dispersion calibrated, alone and together with ExogShockFn
% doPart(8): the AutoCorrTransProbs, CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets: parsing and selectors, the moment vector vs full runs, the selectors handed to the three commands, round trips (rho_z, beta), logmoments by name, runtime
% doPart(9): the AutoCorrTransProbs targets of the window and log rows (fnname_log, fnname_MperiodAverage, fnname_MperiodAverage_log, ...): parsing and selectors (ntimeperiods, logwhichcombos), the moment vector vs full runs at the level and the log
% doPart(10): the AgeConditionalStats targets of the window and log rows (fnname_log, fnname_MperiodAverage, ...): parsing and selectors (blocks of ntimeperiods, logwhichcombos, the split into a second call when the single-period targets are by age groups), the moment vector vs full LifeCycleProfiles runs
% doPart(11): the CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets of the windows (fnname_MperiodAverage, CovarianceMatrix_MperiodAverage, ...): parsing and selectors (block-diagonal whichcombos in blocks of ntimeperiods, the age-conditional split into a second call when the single-period targets are by age groups), the moment vector vs full runs
%
% Parts are independent: the setup, the addpaths and the model preamble all sit OUTSIDE the if-blocks and so always run, and no
% part reads another part's output. Any subset can be run, in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CalibFHorzTestsdiary.txt','file')
    delete('./TestOutput/CalibFHorzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CalibFHorzTestsdiary.txt
fprintf('CalibFHorzTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, CoreStructLeaves, shared by the Core banks
addpath('./CalibFHorzTests_subcodes/')
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/') % the return functions are shared with CoreFHorzTests

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

%% ===== doPart(1): SetupTargetMoments_FHorz =====
if doPart(1)==1
    fprintf('\n===== doPart(1): SetupTargetMoments_FHorz: target vector, selectors, reduced FnsToEvaluate, field order, invalid inputs =====\n')
    %% 1. No model is solved: TargetMoments with arbitrary values of the right sizes, across AllStats and AgeConditionalStats, unrestricted
    % and under two conditional restrictions, every statistic kind, NaN-omitted ages, matrix-valued age-conditional targets (LorenzCurve,
    % QuantileMeans), and agegroupings. The target vector, the per-combination selectors, the union whichstats and the reduced
    % FnsToEvaluate are compared with hand-built versions; three field orders of the same targets must give the same selectors and
    % function sets; the three-input legacy call still works; five invalid inputs must error.
    output=CalibFHorz_TargetParsing(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(1): SetupTargetMoments_FHorz

%% ===== doPart(2): the objective's moment vector =====
if doPart(2)==1
    fprintf('\n===== doPart(2): the objective moment vector vs hand-built from full LifeCycleProfiles/AllStats runs =====\n')
    %% 2. One solve at the truth. Targets are taken from full LifeCycleProfiles and AllStats runs (every statistic kind, two restrictions,
    % NaN-omitted ages, matrix-valued age-conditional targets). CalibrateLifeCycleModel_objectivefn with vectoroutput=1 is called with
    % caliboptions.whichcombos=1 and =0; each vector must equal the hand-built one (the NaN targets dropped). The selectors from
    % SetupTargetMoments_FHorz are then handed straight to the two stats commands: untargeted combinations must come back NaN or absent.
    output=CalibFHorz_MomentVector(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(2): the objective's moment vector

%% ===== doPart(3): objective forms =====
if doPart(3)==1
    fprintf('\n===== doPart(3): objective forms vs their formulas; logmoments and weights through CalibrateLifeCycleModel =====\n')
    %% 3. With perturbed targets, the objective (two calibrated parameters, so the division by their number shows) in each form:
    % sum_squared, sum_logratiosquared, MethodOfMoments with a weight matrix, the lsqnonlin residual vector (its squared norm is the
    % sum_squared objective, its last entry the zero penalty), logmoments, against the formulas on the hand-built moment vector. Then,
    % through CalibrateLifeCycleModel started at the truth with model-generated targets: the objective is zero, calibsummary's moments
    % equal the targets, and caliboptions.logmoments in each of its four forms (scalar, per target entry, per target name, by name)
    % becomes the expected vector; a weights vector of the wrong length errors; whichcombos=0 gives the same moments.
    output=CalibFHorz_ObjectiveForms(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(3): objective forms

%% ===== doPart(4): round trips =====
if doPart(4)==1
    fprintf('\n===== doPart(4): round trips through CalibrateLifeCycleModel =====\n')
    %% 4. Targets generated at the truth. Started at the truth, lsqnonlin and fminsearch return the truth and an objective of zero
    % (exact). Started away, both recover the truth within the tolerances (close), and whichcombos=0 gives exactly the same answer as
    % =1. Constraints (constrainpositive softplus and log, constrain0to1, constrainAtoB), omitcalibparam on the age-dependent kappa_j,
    % and ParametrizeParamsFn each recover from a wrong start.
    output=CalibFHorz_RoundTrips(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(4): round trips

%% ===== doPart(5): shocks that depend on calibrated parameters, jequaloneDist as a function, CustomModelStats =====
if doPart(5)==1
    fprintf('\n===== doPart(5): ExogShockFn with a calibrated parameter, jequaloneDist as a function, CustomModelStats =====\n')
    %% 5. vfoptions.ExogShockFn with the autocorrelation of z calibrated (the shock grids are rebuilt every evaluation): at the truth
    % the objective is zero, from a wrong start rho_z is recovered. jequaloneDist given as a function (both jequaloneDist_usergrids
    % settings) reproduces the matrix version exactly. CustomModelStats targets: a by-hand mean of assets (equal to the AllStats one) and
    % a check that the z_grid handed to the custom function is the user's grid (usergrids=1) or the internal joint grid (usergrids=0).
    output=CalibFHorz_ShocksCustomStats(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(5)

%% ===== doPart(6): runtime =====
if doPart(6)==1
    fprintf('\n===== doPart(6): runtime of one objective evaluation, whichcombos=1 vs =0 (informational) =====\n')
    %% 6. One objective evaluation, timed, for a sparse target set (one statistic of one function under one restriction at five ages)
    % and for a dense one, with caliboptions.whichcombos=1 and =0. Prints the ratios; no identity is checked here.
    output=CalibFHorz_Runtime(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(6): runtime

%% ===== doPart(7): EiidShockFn =====
if doPart(7)==1
    fprintf('\n===== doPart(7): an iid shock through EiidShockFn with its dispersion calibrated, alone and with ExogShockFn =====\n')
    output=CalibFHorz_EiidShock(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(7): EiidShockFn

%% ===== doPart(8): AutoCorrTransProbs, CrossSectionCovarCorr, AgeConditionalCrossSectionCovarCorr targets =====
if doPart(8)==1
    fprintf('\n===== doPart(8): AutoCorrTransProbs, CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets: parsing, moment vector, selectors, round trips, logmoments, runtime =====\n')
    output=CalibFHorz_AutoCorrCrossSec(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(8): the three new target kinds

%% ===== doPart(9): AutoCorrTransProbs targets of the window and log rows =====
if doPart(9)==1
    fprintf('\n===== doPart(9): AutoCorrTransProbs targets of the window and log rows: parsing and selectors, the moment vector vs full runs at the level and the log =====\n')
    output=CalibFHorz_AutoCorrWindows(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(9): the window and log rows

%% ===== doPart(10): AgeConditionalStats targets of the window and log rows =====
if doPart(10)==1
    fprintf('\n===== doPart(10): AgeConditionalStats targets of the window and log rows: parsing and selectors, the split, the moment vector vs full LifeCycleProfiles runs =====\n')
    output=CalibFHorz_ACSWindows(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(10): the AgeConditionalStats window and log rows

%% ===== doPart(11): CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets of the windows =====
if doPart(11)==1
    fprintf('\n===== doPart(11): CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets of the windows: parsing and selectors, the split, the moment vector vs full runs =====\n')
    output=CalibFHorz_CSWindows(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist);
end % doPart(11): the CrossSec window targets

% All looks good!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CalibFHorzTestsdiary.txt')

diary off
