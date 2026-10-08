% Tests for CalibrateInfHorzAgentModel (the single-type InfHorz calibration command) and its target plumbing.
%
% The model is the CoreInfHorzTests nod_z model (ReturnFn_nod_z_noe_nosemiz, from CoreInfHorzTests/CoreInfHorz_ReturnFns) on small
% grids (with more income risk than CoreInfHorz_setup, see the preamble), so that one calibration takes seconds. Targets are generated from the model at known parameter values,
% so the truth is known exactly: a calibration started at the truth must return the truth with an objective of zero, and one
% started away from it must get back within the tolerances.
%
% Checks are identities ("this should be zero") wherever two routes compute the same object: the target vector against a hand-built
% one, the objective's moment vector against a hand-built one from full AllStats/AutoCorrTransProbs/CrossSectionCovarCorr runs, the
% objective forms against their formulas, a calibration started at the truth. Only the optimiser recoveries from a wrong start are
% "close to zero" (held to CoreSummary's looser bar), because they are only as exact as the tolerances.
%
% What the bank is for, 2026-10-07 (built test-first, against the unextended command; the InfHorz counterpart of CalibFHorzTests):
%   - the wiring of simoptions.whichcombos (and a per-combination simoptions.whichstats) into the InfHorz calibration, through
%     caliboptions.whichcombos (=1 default: the three stats commands compute only the targeted combinations; =0: every statistic of
%     every targeted function), with the selectors built by SetupTargetMoments_InfHorz from TargetMoments and simoptions
%   - conditional restrictions on the AutoCorrTransProbs and CrossSectionCovarCorr targets (.(restriction).(fn).(stat), the commands
%     support them since 2026-10-07; the objective used to strip them), horizons read off the _kK suffix of an AutoCorr target name,
%     and the CovarianceMatrix/CorrelationMatrix targets (NaN omits an entry)
%   - calibsummary.currentmomentvec, .targetmomentvec and .logmoments, as CalibrateLifeCycleModel reports them
%   - pre-existing defects pinned here: quantile targets switch on Minimum/Maximum (entry 5) instead of the quantiles (entry 6); the
%     reduced FnsToEvaluate read only the LAST first-level name's second level, so they depend on the field order of TargetMoments; a
%     two-level CrossSectionCovarCorr name as the FIRST cross-section target reads an empty third name; the per-target-name vector form
%     of caliboptions.logmoments reads acsmomentnames, which only exists in the FHorz driver; logmoments by name reads a restricted
%     (three-level) name as two levels, and the cross-section kind under the wrong field name (CrossSecCovarCorr)
%
% Not covered (later programmes): CalibrateInfHorzAgentModel_PType (its own bank, CalibInfHorzPTypeTests), the OLG/GE calibration
% commands, CMA-ES (fminalgo=4). The InfHorz value function has no iid e shock, so EiidShockFn is not exercised.
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): SetupTargetMoments_InfHorz: target vector, selectors (whichcombos/whichstats) for the three commands, horizons, reduced FnsToEvaluate, field-order independence, legacy call, invalid inputs
% doPart(2): the objective's moment vector vs a hand-built one from full AllStats/AutoCorrTransProbs/CrossSectionCovarCorr runs; whichcombos 1 and 0; the selectors handed straight to the three commands
% doPart(3): objective forms (sum_squared, sum_logratiosquared, MethodOfMoments, residual vector, logmoments, weights) vs their formulas; logmoments and weights plumbing through CalibrateInfHorzAgentModel
% doPart(4): round trips through CalibrateInfHorzAgentModel: start at the truth (exact) and away from it (close), lsqnonlin and fminsearch, parameter constraints, omitcalibparam, ParametrizeParamsFn
% doPart(5): shock grids that depend on a calibrated parameter (ExogShockFn, rho_z), CustomModelStats under both usergrids settings and together with ExogShockFn
% doPart(6): round trips on the AutoCorrTransProbs and CrossSectionCovarCorr targets: restricted, with a horizon suffix (rho_z), a four-level cross-section name (beta), matrix targets, logmoments by their names
% doPart(7): runtime of one objective evaluation, whichcombos=1 vs =0 (informational)
%
% Parts are independent: the setup, the addpaths and the model preamble all sit OUTSIDE the if-blocks and so always run, and no
% part reads another part's output. Any subset can be run, in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CalibInfHorzTestsdiary.txt','file')
    delete('./TestOutput/CalibInfHorzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CalibInfHorzTestsdiary.txt
fprintf('CalibInfHorzTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, CoreStructLeaves, CoreMaxAbs, shared by the Core banks
addpath('./CalibInfHorzTests_subcodes/')
addpath('../CoreInfHorzTests/CoreInfHorz_ReturnFns/') % the return functions are shared with CoreInfHorzTests

%% The model: the CoreInfHorzTests nod_z model (CoreInfHorz_setup parameters), with more income risk and a longer asset grid
% CoreInfHorz_setup has an innovation standard deviation of 0.03 and a_grid up to 5. With so little risk the borrowing constraint
% binds for everyone once beta is a little below the truth (run of 2026-10-08: at beta=0.90 the mean of assets was flat in beta, so
% lsqnonlin's finite-difference gradient was exactly zero and it declared the wrong start a local minimum; fminsearch, whose first
% simplex steps are larger, recovered). The risk here is 0.1, so assets stay positive and the moments respond to beta over the whole
% range the round trips start from; the grid is longer for the fatter right tail. Targets are generated from this same model, so
% the truth is still exact.
n_d=0;
d_grid=[];
n_a=101;
a_grid=8*linspace(0,1,n_a)'.^3;
n_z=5;
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.1,n_z);
z_grid=exp(z_grid);

Params.beta=0.95; % discount factor
DiscountFactorParamNames={'beta'};
Params.sigma=2; % CES utility param for consumption
Params.w=1;
Params.r=0.05;

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);

%% ===== doPart(1): SetupTargetMoments_InfHorz =====
if doPart(1)==1
    fprintf('\n===== doPart(1): SetupTargetMoments_InfHorz: target vector, selectors, horizons, reduced FnsToEvaluate, field order, legacy call, invalid inputs =====\n')
    %% 1. No model is solved: TargetMoments with arbitrary values of the right sizes across AllStats, AutoCorrTransProbs and
    % CrossSectionCovarCorr, unrestricted and under two conditional restrictions, every statistic kind, horizon suffixes, a matrix
    % target with NaN entries. The target vector, the per-command selectors, the horizons, the union whichstats and the reduced
    % FnsToEvaluate are compared with hand-built versions; three field orders of the same targets must give the same selectors and
    % function sets; the three-input legacy call still works for unrestricted targets; seven invalid inputs must error.
    output=CalibInfHorz_TargetParsing(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn);
end % doPart(1): SetupTargetMoments_InfHorz

%% ===== doPart(2): the objective's moment vector =====
if doPart(2)==1
    fprintf('\n===== doPart(2): the objective moment vector vs hand-built from full AllStats/AutoCorrTransProbs/CrossSectionCovarCorr runs =====\n')
    %% 2. One solve at the truth. Targets are taken from full runs of the three commands (every statistic kind, two restrictions,
    % two horizons, matrix targets with NaN). CalibrateInfHorzAgentModel_objectivefn with vectoroutput=1 is called with
    % caliboptions.whichcombos=1 and =0; each vector must equal the hand-built one (the NaN targets dropped). The selectors from
    % SetupTargetMoments_InfHorz are then handed straight to the three commands: untargeted combinations must come back NaN or absent.
    output=CalibInfHorz_MomentVector(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn);
end % doPart(2): the objective's moment vector

%% ===== doPart(3): objective forms =====
if doPart(3)==1
    fprintf('\n===== doPart(3): objective forms vs their formulas; logmoments and weights through CalibrateInfHorzAgentModel =====\n')
    %% 3. With perturbed targets, the objective (two calibrated parameters, so the division by their number shows) in each form:
    % sum_squared, sum_logratiosquared, MethodOfMoments with a weight matrix, the lsqnonlin residual vector (its squared norm is the
    % sum_squared objective, its last entry the zero penalty), logmoments, against the formulas on the hand-built moment vector. Then,
    % through CalibrateInfHorzAgentModel started at the truth with model-generated targets: the objective is zero, calibsummary's
    % moments equal the targets, and caliboptions.logmoments in each of its four forms (scalar, per target entry, per target name, by
    % name) becomes the expected vector; a weights vector of the wrong length errors; whichcombos=0 gives the same moments.
    output=CalibInfHorz_ObjectiveForms(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn);
end % doPart(3): objective forms

%% ===== doPart(4): round trips =====
if doPart(4)==1
    fprintf('\n===== doPart(4): round trips through CalibrateInfHorzAgentModel =====\n')
    %% 4. Targets generated at the truth. Started at the truth, lsqnonlin and fminsearch return the truth and an objective of zero
    % (exact). Started away, both recover the truth within the tolerances (close), and whichcombos=0 gives exactly the same answer as
    % =1. Constraints (constrainpositive softplus and log, constrain0to1, constrainAtoB), omitcalibparam on a vector parameter, and
    % ParametrizeParamsFn each recover from a wrong start.
    output=CalibInfHorz_RoundTrips(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn);
end % doPart(4): round trips

%% ===== doPart(5): shocks that depend on calibrated parameters, CustomModelStats =====
if doPart(5)==1
    fprintf('\n===== doPart(5): ExogShockFn with a calibrated parameter, CustomModelStats =====\n')
    %% 5. vfoptions.ExogShockFn with the autocorrelation of z calibrated (the shock grids are rebuilt every evaluation): at the truth
    % the objective is zero, from a wrong start rho_z is recovered. CustomModelStats targets: a by-hand mean of assets (equal to the
    % AllStats one) and a check that the z_grid handed to the custom function is the user's grid (usergrids=1) or the internal joint
    % grid (usergrids=0); and usergrids=1 together with ExogShockFn (the user's grid as rebuilt from the current parameters).
    output=CalibInfHorz_ShocksCustomStats(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn);
end % doPart(5)

%% ===== doPart(6): restricted AutoCorrTransProbs and CrossSectionCovarCorr targets, round trips and logmoments =====
if doPart(6)==1
    fprintf('\n===== doPart(6): round trips on the AutoCorrTransProbs and CrossSectionCovarCorr targets (restricted, horizons, matrices), logmoments by name =====\n')
    %% 6. rho_z (ExogShockFn) from the autocorrelation of assets among the rich at horizons one and two and the autocorrelation of
    % earnings; beta from the assets-consumption correlation among the rich and the assets-earnings covariance, with whichcombos=0
    % identical; the two matrix targets at the truth; caliboptions.logmoments by a three-level and a four-level name.
    output=CalibInfHorz_RestrictedTargets(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn);
end % doPart(6)

%% ===== doPart(7): runtime =====
if doPart(7)==1
    fprintf('\n===== doPart(7): runtime of one objective evaluation, whichcombos=1 vs =0 (informational) =====\n')
    %% 7. One objective evaluation, timed, for a sparse target set and for a dense one, with caliboptions.whichcombos=1 and =0.
    % Prints the ratios; no identity is checked here.
    output=CalibInfHorz_Runtime(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn);
end % doPart(7)

% All looks good!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CalibInfHorzTestsdiary.txt')

diary off
