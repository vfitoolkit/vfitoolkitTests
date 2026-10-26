% Tests for CalibrateInfHorzAgentModel_PType (the InfHorz calibration command with permanent types) and its target plumbing.
%
% The PType counterpart of CalibInfHorzTests, on the CoreInfHorzPTypeTests model: two permanent types, low and high, that differ by
% a per-type wage w_pt (a vector with one entry per type, which PType reads as a type-dependent parameter), type masses 0.6 and 0.4,
% on the single-type bank's grids and income risk (innovation sd 0.1, a_grid to 8; see CalibInfHorzTests), so that one calibration
% takes seconds. Targets are generated from the model at known parameter values, so the truth is known exactly.
%
% What is specific to permanent types: a target has two to five names (fn.stat, fn.type.stat, restriction.fn.stat,
% restriction.fn.type.stat, restriction.fn1.type.CovarianceWith.fn2), and the selectors that caliboptions.whichcombos=1 hands to the
% three InfHorz PType stats commands carry a trailing type dimension of N_i+1: one slot per type in the order of Names_i, then the
% grouped statistics. A grouped slot forces the per-type statistics of that combination for every type (the grouped Mean, StdDeviation,
% Minimum and Maximum are built from them), and those are reported. Calibrated parameters can be common or per type (as a structure
% over the types, or as a vector/matrix with N_i as one dimension), and the shock functions can differ by type.
%
% Checks are identities ("this should be zero") wherever two routes compute the same object; only the optimiser recoveries from a
% wrong start are "close to zero".
%
% Built test-first 2026-10-09 against the unextended command, which has no whichcombos wiring, errors on AutoCorrTransProbs and
% CrossSectionCovarCorr targets ("not yet implemented"), reports only calibsummary.objvalue, cannot hold a vector parameter under
% omitcalibparam, and reads a structure of per-type ExogShockFn as one function. Not targetable: the self-correlation fn.fn, the
% per-type matrices (CovarianceMatrix_ptype), TransitionProbs. The InfHorz value function has no iid e shock, so EiidShockFn is not
% exercised; there is no jequaloneDist.
%
% Not covered: the OLG/GE calibration commands, CMA-ES (fminalgo=4).
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): SetupTargetMoments_InfHorz with permanent types: target vector, type-slot selectors for the three commands, horizons, field order, legacy call, invalid inputs
% doPart(2): the objective's moment vector vs hand-built from full PType AllStats/AutoCorrTransProbs/CrossSectionCovarCorr runs; whichcombos 1 and 0; the selectors handed to the PType commands
% doPart(3): objective forms vs their formulas; logmoments (incl. four- and five-level names) and weights through CalibrateInfHorzAgentModel_PType
% doPart(4): round trips: common and per-type parameters (vector, matrix and structure forms), constraints, omitcalibparam, ParametrizeParamsFn, the type weights calibrated
% doPart(5): ExogShockFn common and per type with a calibrated parameter, CustomModelStats under both grid settings
% doPart(6): runtime of one objective evaluation, whichcombos=1 vs =0 (informational)
%
% Parts are independent: the setup, the addpaths and the model preamble all sit OUTSIDE the if-blocks and so always run, and no
% part reads another part's output. Any subset can be run, in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CalibInfHorzPTypeTestsdiary.txt','file')
    delete('./TestOutput/CalibInfHorzPTypeTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CalibInfHorzPTypeTestsdiary.txt
fprintf('CalibInfHorzPTypeTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, CoreStructLeaves, CoreMaxAbs, shared by the Core banks
addpath('./CalibInfHorzPTypeTests_subcodes/')
addpath('../../CoreInfHorzTests/CoreInfHorz_ReturnFns/') % the return functions are shared with CoreInfHorzTests

%% The model: the CoreInfHorzPTypeTests model (two types by wage) on the CalibInfHorzTests grids and income risk
n_d=0;
d_grid=[];
n_a=101;
a_grid=8*linspace(0,1,n_a)'.^3;
n_z=5;
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.1,n_z);
z_grid=exp(z_grid);
Names_i={'low','high'};

Params.beta=0.95; % discount factor
DiscountFactorParamNames={'beta'};
Params.ptypemass=[0.6;0.4];
PTypeDistParamNames={'ptypemass'};
Params.sigma=2; % CES utility param for consumption
Params.w_pt=[1;1.4]; % the wage, per type (N_i entries: PType reads it as type-dependent)
Params.r=0.05;

ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);

%% ===== doPart(1): SetupTargetMoments_InfHorz with permanent types =====
if doPart(1)==1
    fprintf('\n===== doPart(1): SetupTargetMoments_InfHorz with permanent types: target vector, type-slot selectors, horizons, field order, legacy call, invalid inputs =====\n')
    output=CalibInfHorzPType_TargetParsing(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn);
end % doPart(1)

%% ===== doPart(2): the objective's moment vector =====
if doPart(2)==1
    fprintf('\n===== doPart(2): the objective moment vector vs hand-built from full PType AllStats/AutoCorrTransProbs/CrossSectionCovarCorr runs =====\n')
    output=CalibInfHorzPType_MomentVector(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn);
end % doPart(2)

%% ===== doPart(3): objective forms =====
if doPart(3)==1
    fprintf('\n===== doPart(3): objective forms vs their formulas; logmoments and weights through CalibrateInfHorzAgentModel_PType =====\n')
    output=CalibInfHorzPType_ObjectiveForms(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn);
end % doPart(3)

%% ===== doPart(4): round trips =====
if doPart(4)==1
    fprintf('\n===== doPart(4): round trips through CalibrateInfHorzAgentModel_PType =====\n')
    output=CalibInfHorzPType_RoundTrips(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn);
end % doPart(4)

%% ===== doPart(5): shocks, CustomModelStats =====
if doPart(5)==1
    fprintf('\n===== doPart(5): ExogShockFn common and per type with a calibrated parameter, CustomModelStats =====\n')
    output=CalibInfHorzPType_ShocksCustomStats(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn);
end % doPart(5)

%% ===== doPart(6): runtime =====
if doPart(6)==1
    fprintf('\n===== doPart(6): runtime of one objective evaluation, whichcombos=1 vs =0 (informational) =====\n')
    output=CalibInfHorzPType_Runtime(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn);
end % doPart(6)

% All looks good!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CalibInfHorzPTypeTestsdiary.txt')

diary off
