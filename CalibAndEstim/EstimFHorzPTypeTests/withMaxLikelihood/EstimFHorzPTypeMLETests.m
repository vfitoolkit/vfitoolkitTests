% Tests for maximum likelihood estimation of life-cycle models with permanent types: LogLikelihood_FHorz_PType (the finite mixture over
% the types of the single-type LogLikelihood_FHorz) and EstimateLifeCycleModel_PType_MaximumLikelihood.
%
% A sub-bank of EstimFHorzPTypeTests: it uses the same model (EstimFHorzPTypeTests_subcodes/EstimFHorzPType_setup.m, two permanent types
% that differ in their earnings profile kappa_j_pt, with masses 0.7 and 0.3), with the grid interpolation layer, shares its TestOutput
% folder, and points to the EstimFHorzPTypeTests, EstimFHorzTests/withMaxLikelihood and CoreFHorzTests subcodes rather than copying them.
%
%   doPart(1): LogLikelihood_FHorz_PType vs the single-type command called by hand type by type, with the mixture, the posterior, observed
%              types and the conditional likelihood by hand; jequaloneDist forms; refusals
%   doPart(2): EstimateLifeCycleModel_PType_MaximumLikelihood at the truth (skipestimation=1), with a common parameter, a per-type
%              parameter and the type masses (via ParametrizeParamsFn) estimated: the scores, OPG, Hessian, sandwich, standard errors and
%              confidence intervals in the per-type layout, recomputed from direct calls; the per-type parameter as a matrix; refusals
%   doPart(3): EstimateLifeCycleModel_PType_MaximumLikelihood point estimates (CMA-ES from the truth and from away): "close to zero"
% No Monte Carlo. The data are panels simulated by hand (EstimFHorzPTypeMLE_SimulatePTypePanel, a plain loop with a fixed seed), which
% is allowed in a test bank but never in the toolkit.
%
% Built 2026-10-16 as Phase 3 of the maximum likelihood project (test-first, with LogLikelihood_FHorz_PType and
% EstimateLifeCycleModel_PType_MaximumLikelihood).
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): LogLikelihood_FHorz_PType: mixture, posterior, zero-mass type, observed types, conditional, classification, jequaloneDist forms, refusals
% doPart(2): EstimateLifeCycleModel_PType_MaximumLikelihood standard errors at the truth (skipestimation=1), matrix form, refusals
% doPart(3): EstimateLifeCycleModel_PType_MaximumLikelihood point estimates: CMA-ES from the truth and from away
%
% Parts are independent: the setup, the addpaths and the model preamble all sit OUTSIDE the if-blocks and so always run, and no
% part reads another part's output. Any subset can be run, in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1];

%% Diary of the command window output (in the parent bank's TestOutput folder)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/EstimFHorzPTypeMLETestsdiary.txt','file')
    delete('../TestOutput/EstimFHorzPTypeMLETestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/EstimFHorzPTypeMLETestsdiary.txt
fprintf('EstimFHorzPTypeMLETests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../../SharedSubcodes/') % CoreSummary, CoreMaxAbs, shared by the Core banks
addpath('../EstimFHorzPTypeTests_subcodes/') % EstimFHorzPType_setup, the model shared with EstimFHorzPTypeTests
addpath('../../../CoreFHorzTests/CoreFHorz_ReturnFns/') % the return functions are shared with CoreFHorzTests
addpath('../../EstimFHorzTests/withMaxLikelihood/withMaxLikelihood_subcodes/') % EstimFHorzMLE_SimulatePanel, the single-type panel simulator
addpath('./withMaxLikelihood_subcodes/')

%% The model: the EstimFHorzPTypeTests model, with the grid interpolation layer
EstimFHorzPType_setup
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;

%% ===== doPart(1): LogLikelihood_FHorz_PType =====
if doPart(1)==1
    fprintf('\n===== doPart(1): LogLikelihood_FHorz_PType: the mixture over the types vs single-type calls by hand =====\n')
    output=EstimFHorzPTypeMLE_LogLik(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions);
end % doPart(1)

%% ===== doPart(2): standard errors at the truth =====
if doPart(2)==1
    fprintf('\n===== doPart(2): EstimateLifeCycleModel_PType_MaximumLikelihood standard errors at the truth (skipestimation=1) =====\n')
    output=EstimFHorzPTypeMLE_StdErrors(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions);
end % doPart(2)

%% ===== doPart(3): point estimates =====
if doPart(3)==1
    fprintf('\n===== doPart(3): EstimateLifeCycleModel_PType_MaximumLikelihood point estimates: CMA-ES from the truth and from away =====\n')
    output=EstimFHorzPTypeMLE_PointEstimates(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions);
end % doPart(3)

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/EstimFHorzPTypeMLETestsdiary.txt')

diary off
