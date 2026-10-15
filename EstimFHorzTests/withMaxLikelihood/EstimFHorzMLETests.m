% Tests for maximum likelihood estimation of life-cycle models: LogLikelihood_FHorz (the log-likelihood of panel data, evaluated
% exactly on the model grid by treating the model as a hidden Markov model and running the forward algorithm), and later
% EstimateLifeCycleModel_MaximumLikelihood.
%
% A sub-bank of EstimFHorzTests: it uses the same model (EstimFHorzTests_subcodes/EstimFHorz_setup.m, the CoreFHorzTests nod_z_noe
% life-cycle model on small grids, with the grid interpolation layer as in the EstimFHorzTests subcodes), shares its TestOutput
% folder, and points to the EstimFHorzTests and CoreFHorzTests subcodes rather than copying them.
%
% Checks are identities ("this should be zero") between two routes to the same number:
%   doPart(1): the forward algorithm against an enumeration of every latent path (never merging paths), on tiny models with a
%              hand-built Policy, with the model values, lotteries and emission densities all computed by hand
%   doPart(2): the EstimFHorzTests model solved (as is, and with an iid e added): the population distribution against
%              StationaryDist_FHorz_Case1, a cross-section against the closed-form mixture over StationaryDist, and the data forms /
%              batching / conditional initial condition against each other
%   doPart(3): inputs that must be refused, each with its intended message
%   doPart(4): EstimateLifeCycleModel_MaximumLikelihood at the truth (skipestimation=1): the scores, OPG, Hessian, sandwich, standard
%              errors, confidence intervals, information criteria and diagnostics recomputed from direct LogLikelihood_FHorz calls
%   doPart(5): EstimateLifeCycleModel_MaximumLikelihood point estimates (CMA-ES, constrained, fminsearch): "close to zero" recoveries
%   doPart(6): AR(1) measurement error components vs an enumeration of every latent path (the AR state folded into the Markov chain)
%   doPart(7): AR(1) measurement error on the solved model: population unchanged, cross-section closed form, the estimator's scores
%   doPart(8): AR(1) measurement error point estimates: rho, the AR innovation s.d. and the iid s.d. recovered (CMA-ES)
% No Monte Carlo. Parts 1-3 use deterministic observations placed near the model values; parts 4-5 use ONE panel simulated by hand
% (EstimFHorzMLE_SimulatePanel, a plain loop with a fixed seed), which is allowed in a test bank but never in the toolkit.
%
% Built 2026-10-15 as Phase 1 of the maximum likelihood project (test-first, with LogLikelihood_FHorz); parts 4-5 are Phase 2
% (EstimateLifeCycleModel_MaximumLikelihood).
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): brute force: six tiny models (no shocks; z; e with grid interpolation; d,z,e with grid interpolation; two endogenous states; conditional + agej0 + batches), every emission distribution
% doPart(2): the EstimFHorzTests model (z) and the same with an iid e: PopDist vs StationaryDist, cross-section vs closed form, all-missing, Data.age vs Data.agej0, batches, conditional, runtime
% doPart(3): refusals
% doPart(4): EstimateLifeCycleModel_MaximumLikelihood standard errors at the truth (skipestimation=1), stderrors options, omitestimparam, refinecheck, refusals
% doPart(5): EstimateLifeCycleModel_MaximumLikelihood point estimates: CMA-ES from the truth, CMA-ES from away with constraints, fminsearch from away
% doPart(6): AR(1) measurement error: brute force, two tiny models (one AR variable with z; two AR variables, one log, conditional)
% doPart(7): AR(1) measurement error on the solved model: PopDist, cross-section closed form, all missing, estimator scores at the truth, refusals
% doPart(8): AR(1) measurement error point estimates (CMA-ES from away), and the fit without the AR(1) component (informational)
%
% Parts are independent: the setup, the addpaths and the model preamble all sit OUTSIDE the if-blocks and so always run, and no
% part reads another part's output. Any subset can be run, in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1];

%% Diary of the command window output (in the parent bank's TestOutput folder)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/EstimFHorzMLETestsdiary.txt','file')
    delete('../TestOutput/EstimFHorzMLETestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/EstimFHorzMLETestsdiary.txt
fprintf('EstimFHorzMLETests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, CoreMaxAbs, shared by the Core banks
addpath('../EstimFHorzTests_subcodes/') % EstimFHorz_setup, the model shared with EstimFHorzTests
addpath('../../CoreFHorzTests/CoreFHorz_ReturnFns/') % the return functions are shared with CoreFHorzTests
addpath('./withMaxLikelihood_subcodes/')

%% The model: the EstimFHorzTests model, with the grid interpolation layer as in the EstimFHorzTests subcodes
EstimFHorz_setup
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
% measurement error parameters
Params.sig_c=0.2;
Params.sig_y=0.1;
Params.nu=6;

%% ===== doPart(1): brute force =====
if doPart(1)==1
    fprintf('\n===== doPart(1): brute force: the forward algorithm vs an enumeration of every latent path, six tiny models =====\n')
    output=EstimFHorzMLE_BruteForce();
end % doPart(1)

%% ===== doPart(2): solved model =====
if doPart(2)==1
    fprintf('\n===== doPart(2): solved model: PopDist vs StationaryDist, cross-section vs closed form, missing data, data forms, batches, conditional =====\n')
    output=EstimFHorzMLE_SolvedModel(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions);
end % doPart(2)

%% ===== doPart(3): refusals =====
if doPart(3)==1
    fprintf('\n===== doPart(3): refusals =====\n')
    output=EstimFHorzMLE_Refusals();
end % doPart(3)

%% ===== doPart(4): standard errors at the truth =====
if doPart(4)==1
    fprintf('\n===== doPart(4): EstimateLifeCycleModel_MaximumLikelihood standard errors at the truth (skipestimation=1): scores, OPG, Hessian, sandwich, options, refusals =====\n')
    output=EstimFHorzMLE_StdErrors(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions);
end % doPart(4)

%% ===== doPart(5): point estimates =====
if doPart(5)==1
    fprintf('\n===== doPart(5): EstimateLifeCycleModel_MaximumLikelihood point estimates: CMA-ES, constrained, fminsearch =====\n')
    output=EstimFHorzMLE_PointEstimates(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions);
end % doPart(5)

%% ===== doPart(6): AR(1) measurement error, brute force =====
if doPart(6)==1
    fprintf('\n===== doPart(6): AR(1) measurement error: the forward algorithm vs an enumeration of every latent path =====\n')
    output=EstimFHorzMLE_BruteForceAR();
end % doPart(6)

%% ===== doPart(7): AR(1) measurement error on the solved model =====
if doPart(7)==1
    fprintf('\n===== doPart(7): AR(1) measurement error on the solved model: PopDist, cross-section closed form, all missing, estimator scores, refusals =====\n')
    output=EstimFHorzMLE_AR(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions);
end % doPart(7)

%% ===== doPart(8): AR(1) measurement error point estimates =====
if doPart(8)==1
    fprintf('\n===== doPart(8): AR(1) measurement error point estimates =====\n')
    output=EstimFHorzMLE_ARPointEstimates(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions);
end % doPart(8)

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/EstimFHorzMLETestsdiary.txt')

diary off
