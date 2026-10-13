% Tests for CalibrateLifeCycleModel_PType (the FHorz calibration command with permanent types) and its target plumbing.
%
% The PType counterpart of CalibFHorzTests, on the CoreFHorzPTypeTests part-21 model: two permanent types, low and high, with a
% per-type deterministic earnings profile (kappa_j_pt), type weights 0.7 and 0.3 and geometric age weights, on small grids so that
% one calibration takes seconds. Targets are generated from the model at known parameter values, so the truth is known exactly.
%
% What is specific to permanent types: a target has two to four names (fn.stat, fn.type.stat, restriction.fn.stat,
% restriction.fn.type.stat), and the selectors that caliboptions.whichcombos=1 hands to the two stats commands carry a trailing type
% dimension of N_i+1: one slot per type in the order of Names_i, then the grouped statistics. A grouped slot forces the per-type
% statistics of that combination for every type (the grouped Mean, StdDeviation, Minimum and Maximum are built from them), and
% those are reported. Calibrated parameters can be common or per type (as a structure over the types, or as a matrix with N_i as
% one dimension), and the shock functions can differ by type.
%
% Checks are identities ("this should be zero") wherever two routes compute the same object; only the optimiser recoveries from a
% wrong start are "close to zero". Parts 4 and 5 use the grid interpolation layer so the moments vary smoothly with the parameters.
%
% Pre-existing defects pinned here (fixed 2026-10-07 with the wiring): the omit-mask transpose test on the wrong variable, the
% by-name and per-target-name forms of logmoments, a per-type (structure) ExogShockFn erroring in the calibrated-shocks detection.
%
% Added 2026-10-08 (part 8): the AutoCorrTransProbs, CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets with permanent
% types, grouped (fn.stat, fn.CovarianceWith.fn2, the grouped matrices) or per type (fn.type.stat, fn.type.CovarianceWith.fn2); the
% selectors of the three PType commands carry the same trailing type dimension. Not targetable: the self-correlation fn.fn, the per-type
% matrices (they shrink to the type's evaluated functions), TransitionProbs; not supported: per-type N_j or agejshifter with these kinds.
%
% Not covered (later programmes): the estimation commands, the OLG/GE calibration commands, CMA-ES.
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): SetupTargetMoments_FHorz with permanent types: target vector, type-slot selectors, per-type age groups and agejshifter, field order, legacy call, invalid inputs
% doPart(2): the objective's moment vector vs hand-built from full PType LifeCycleProfiles/AllStats runs; whichcombos 1 and 0; the selectors handed to the PType commands
% doPart(3): objective forms vs their formulas; logmoments (incl. a four-level name) and weights through CalibrateLifeCycleModel_PType
% doPart(4): round trips: common and per-type parameters (structure and matrix forms), constraints, omitcalibparam, ParametrizeParamsFn
% doPart(5): ExogShockFn common and per type with a calibrated parameter, jequaloneDist as a function, CustomModelStats under both grid settings
% doPart(6): runtime of one objective evaluation, whichcombos=1 vs =0 (informational)
% doPart(7): extensions: the type weights calibrated (via ParametrizeParamsFn, and the vector itself), per-type N_j with and without agejshifter, an iid shock through EiidShockFn (alone and with ExogShockFn), a semi-exogenous state with SemiExoStateFn common and per type
% doPart(8): the AutoCorrTransProbs, CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets with permanent types (grouped and per-type names, type-slot selectors): parsing, the moment vector vs full PType runs, the selectors handed to the three PType commands, round trips (common rho_z, beta + a per-type pension), logmoments by a four-level name, runtime
% doPart(9): the AutoCorrTransProbs targets of the window and log rows (fnname_log, fnname_MperiodAverage, fnname_MperiodAverage_log, ...) with permanent types, grouped and per type: parsing and selectors (ntimeperiods, logwhichcombos), the moment vector vs full runs at the level and the log
%
% Parts are independent: the setup, the addpaths and the model preamble all sit OUTSIDE the if-blocks and so always run, and no
% part reads another part's output. Any subset can be run, in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CalibFHorzPTypeTestsdiary.txt','file')
    delete('./TestOutput/CalibFHorzPTypeTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CalibFHorzPTypeTestsdiary.txt
fprintf('CalibFHorzPTypeTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, CoreStructLeaves, shared by the Core banks
addpath('./CalibFHorzPTypeTests_subcodes/')
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/') % the return functions are shared with CoreFHorzTests
addpath('../CoreFHorzTests/CoreFHorzTests_Setup/') % CoreFHorzSetup_SemiExoStateFn, for the semi-exogenous model of part 7

%% The model: the CoreFHorzPTypeTests part-21 model on small grids
n_d=0;
d_grid=[];
n_a=101;
a_grid=5*linspace(0,1,n_a)'.^3;
n_z=5;
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
z_grid=exp(z_grid);
N_j=20;
Names_i={'low','high'};

Params.beta=0.95; % discount factor
DiscountFactorParamNames={'beta'};
Params.mewj=0.9.^(0:N_j-1); Params.mewj=Params.mewj/sum(Params.mewj);
AgeWeightParamNames={'mewj'};
Params.ptypeweights=[0.7;0.3];
PTypeDistParamNames={'ptypeweights'};
Params.sigma=2; % CES utility param for consumption
Params.w=1;
Params.r=0.05;
Params.Jr=16;
Params.pension=0.5;
Params.agej=1:1:N_j;
Params.kappa_j=[0.5:0.1:1,ones(1,9),zeros(1,5)]; % deterministic earnings profile, zero in retirement
Params.kappa_j_pt=[Params.kappa_j; 1.5*Params.kappa_j]; % per type: N_i is the first dimension

ReturnFn=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);

jequaloneDist=zeros(n_a,n_z,'gpuArray'); % everyone starts with no assets at the middle shock (both types)
jequaloneDist(1,ceil(n_z/2))=1;

%% ===== doPart(1): SetupTargetMoments_FHorz with permanent types =====
if doPart(1)==1
    fprintf('\n===== doPart(1): SetupTargetMoments_FHorz with permanent types: target vector, type-slot selectors, age groups, field order, invalid inputs =====\n')
    output=CalibFHorzPType_TargetParsing(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(1)

%% ===== doPart(2): the objective's moment vector =====
if doPart(2)==1
    fprintf('\n===== doPart(2): the objective moment vector vs hand-built from full PType LifeCycleProfiles/AllStats runs =====\n')
    output=CalibFHorzPType_MomentVector(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(2)

%% ===== doPart(3): objective forms =====
if doPart(3)==1
    fprintf('\n===== doPart(3): objective forms vs their formulas; logmoments and weights through CalibrateLifeCycleModel_PType =====\n')
    output=CalibFHorzPType_ObjectiveForms(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(3)

%% ===== doPart(4): round trips =====
if doPart(4)==1
    fprintf('\n===== doPart(4): round trips through CalibrateLifeCycleModel_PType =====\n')
    output=CalibFHorzPType_RoundTrips(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(4)

%% ===== doPart(5): shocks, jequaloneDist as a function, CustomModelStats =====
if doPart(5)==1
    fprintf('\n===== doPart(5): ExogShockFn common and per type with a calibrated parameter, jequaloneDist as a function, CustomModelStats =====\n')
    output=CalibFHorzPType_ShocksCustomStats(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(5)

%% ===== doPart(6): runtime =====
if doPart(6)==1
    fprintf('\n===== doPart(6): runtime of one objective evaluation, whichcombos=1 vs =0 (informational) =====\n')
    output=CalibFHorzPType_Runtime(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(6)

%% ===== doPart(7): extensions =====
if doPart(7)==1
    fprintf('\n===== doPart(7): type weights calibrated, per-type N_j, EiidShockFn, semi-exogenous state =====\n')
    output=CalibFHorzPType_Extensions(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(7)

%% ===== doPart(8): AutoCorrTransProbs, CrossSectionCovarCorr, AgeConditionalCrossSectionCovarCorr targets =====
if doPart(8)==1
    fprintf('\n===== doPart(8): AutoCorrTransProbs, CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets with permanent types: parsing, moment vector, selectors, round trips, logmoments, runtime =====\n')
    output=CalibFHorzPType_AutoCorrCrossSec(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(8)

%% ===== doPart(9): AutoCorrTransProbs targets of the window and log rows =====
if doPart(9)==1
    fprintf('\n===== doPart(9): AutoCorrTransProbs targets of the window and log rows with permanent types: parsing and selectors, the moment vector vs full runs at the level and the log =====\n')
    output=CalibFHorzPType_AutoCorrWindows(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist);
end % doPart(9): the window and log rows

% All looks good!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CalibFHorzPTypeTestsdiary.txt')

diary off
