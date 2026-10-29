% Tests for the general-eqm calibration commands: CalibrateOLGModel (FHorz) and CalibrateBIHAModel (InfHorz), each in its two modes,
% caliboptions.jointoptimization=1 (the GE prices calibrated alongside the parameters, the GE conditions part of the objective) and =0
% (nested: every objective evaluation solves the general eqm with HeteroAgentStationaryEqm and evaluates the targets there).
%
% The model is the CoreStationaryGeneralEqm model (Aiyagari-style, endogenous labour, a government; GE prices r, Tr, tau_c solving
% CapitalMarket, GovBudget, ConsTax; ReturnFns from CoreStationaryGeneralEqm/CoreStationaryGeneralEqm_ReturnFns) on smaller grids. The
% calibrated parameters are beta (identified by aggregate capital) and varphi (by mean hours). Targets are generated from the model at
% the truth and its general eqm, so the truth is known. Where a calibration has to solve many general eqms (nested round trips), a
% one-price version of the model is used: tau=G=Tr=tau_c=0, so only r and CapitalMarket remain.
%
% Checks are identities ("this should be zero") wherever two routes compute the same object: the objective's moment vector against one
% built by hand from full AllStats/LifeCycleProfiles/AutoCorrTransProbs/CrossSectionCovarCorr runs; its GE conditions against
% HeteroAgentStationaryEqm with maxiter=0 at the same prices; the objective forms against their formulas; joint against nested at the
% same prices. Comparisons with an equilibrium, which is only as exact as its tolerance, are "should be near zero", and optimiser
% recoveries from a wrong start are "close to zero".
%
% What the bank pins, 2026-10-07 (all found by reading, fixed in the same change, none seen failing on a GPU run first):
%   both  - target setup was the legacy 3-input SetupTargetMoments call: no whichcombos, and (OLG) none of the AutoCorrTransProbs,
%           CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets that CalibrateLifeCycleModel has
%         - caliboptions.logmoments: by name read AgeConditionalStats from the AllStats names; the per-name vector form used undefined
%           sizes with an off-by-one (BIHA: referred to acsmomentnames, which InfHorz does not have)
%         - nested mode handed HeteroAgentStationaryEqm a heteroagentoptions whose GE-price constraints had already been turned from
%           names into 0/1 vectors (and the internal grids/options); it now gets the user's own
%         - joint mode: V was only kept for caliboptions.CustomModelStats, so heteroagentoptions.CustomModelStats alone crashed
%         - nested mode did not return the GE prices; both modes now put them in CalibParams, and calibsummary has currentmomentvec,
%           targetmomentvec, logmoments and GeneralEqmConditions
%         - the scalar objective divided the targets part by the number of parameters INCLUDING the GE prices in joint mode; both modes
%           now divide by the number of calibrated (non-GE) parameters
%         - the verbose print of the joint objective had three placeholders and two arguments, and printed the residual VECTORS under
%           lsqnonlin (dozens of lines); the penalty message tested a vector. Print-only (also fixed in the two PType joint and the two
%           PType nested objectives)
%   OLG   - omitcalibparam: a row omit mask was never transposed; heteroagentoptions.multiGEweights as a column made a matrix
%   BIHA  - the nested objective used undefined Params and z_grid and the cell of GE eqns where the structure is needed: it could
%           never have run; CustomModelStats was handed the internal vfoptions/simoptions (alreadygridvals=1); a two-name
%           CrossSectionCovarCorr target read three names
%   both  - joint mode joined the GE prices' constraint flags onto the calibrated parameters' but not their per-parameter softplus
%           flags, so a constrainpositive GE price indexed past the end; joined now, so each GE price also keeps
%           heteroagentoptions.constrainpositivemethod
%   Added the same day (user decisions):
%   both  - under lsqnonlin the residual vector weighted the GE conditions by relativeGEweight and multiGEweights, so its sum of squares
%           weighted them by their SQUARES, and the targets part was not divided as in the scalar form; now the vector carries the
%           square roots, and its sum of squares is exactly the scalar objective (also in the two PType joint objectives)
%   all calibration drivers (LifeCycle, InfHorzAgent, OLG, BIHA, each with PType) - MethodOfMoments under lsqnonlin handed the
%           objective the ELEMENTWISE sqrt of the weighting matrix, right only for a diagonal one; now chol(W,'upper'), as the
%           estimation commands do (parts 2 and 8 check the reported lsqnonlin objective against the formula at the solution)
%   HeteroAgentStationaryEqm_InfHorz handed heteroagentoptions.CustomModelStats its internal vfoptions/simoptions (alreadygridvals=1);
%           now the user's own, as FHorz does (parts 5 and 11, nested)
%   the PType drivers had the same missing concatenation of the softplus flags (not exercised here until the PType parts are added)
%   2026-10-13: CalibrateBIHAModel ported to the rewritten SetupTargetMoments_InfHorz (CalibInfHorzTests, 461873b2): caliboptions.whichcombos
%           (default 1) with its selectors, restricted AutoCorr/CrossSec targets, _horizonK horizons, and an omit mask of any length (was 1-by-1);
%           its objectives take CalibrateInfHorzAgentModel_objectivefn's stats and moment-vector blocks. Part 7 checks whichcombos 1 vs 0.
%
% Not covered yet: the PType versions (CalibrateOLGModel_PType, CalibrateBIHAModel_PType) - they go into this same bank later, with a
% round trip for a per-type parameter given as a matrix; the full InfHorz target battery (quantiles, restricted AutoCorr/CrossSec,
% horizons) is CalibInfHorzTests' job and is not repeated here; CMA-ES (fminalgo=4).
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1):  FHorz (OLG): the objective's moment vector and GE conditions vs hand-built (joint), whichcombos 1 and 0
% doPart(2):  FHorz (OLG): objective forms vs their formulas, joint and nested scaling, multiGEweights as row/column, logmoments plumbing, verbose prints
% doPart(3):  FHorz (OLG): joint vs nested at the same prices; the inner GE solve; constraints reach the inner solver; GE prices in CalibParams
% doPart(4):  FHorz (OLG): round trips, joint (three prices) and nested (one price), lsqnonlin and fminsearch, constraints, omitcalibparam
% doPart(5):  FHorz (OLG): ExogShockFn with a calibrated parameter, CustomModelStats (caliboptions, both usergrids; heteroagentoptions), intermediateEqns, jequaloneDist as a function
% doPart(6):  FHorz (OLG): runtime of one objective evaluation: joint whichcombos=1 vs =0, nested (informational)
% doPart(7):  InfHorz (BIHA): the objective's moment vector and GE conditions vs hand-built (joint)
% doPart(8):  InfHorz (BIHA): objective forms, joint and nested scaling, multiGEweights as row/column, logmoments plumbing, verbose prints
% doPart(9):  InfHorz (BIHA): joint vs nested at the same prices; the inner GE solve; constraints reach the inner solver; GE prices in CalibParams
% doPart(10): InfHorz (BIHA): round trips, joint (three prices) and nested (one price), lsqnonlin and fminsearch, constraints
% doPart(11): InfHorz (BIHA): ExogShockFn with a calibrated parameter, CustomModelStats (caliboptions, both usergrids; heteroagentoptions), intermediateEqns
%
% Parts are independent: the setup, the addpaths and the model preamble all sit OUTSIDE the if-blocks and so always run, and no
% part reads another part's output. Any subset can be run, in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CalibStationaryGeneralEqmTestsdiary.txt','file')
    delete('./TestOutput/CalibStationaryGeneralEqmTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CalibStationaryGeneralEqmTestsdiary.txt
fprintf('CalibStationaryGeneralEqmTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CalibStationaryGeneralEqmTests_subcodes/')
addpath('../../CoreStationaryGeneralEqm/CoreStationaryGeneralEqm_ReturnFns/') % the return functions are shared with CoreStationaryGeneralEqm

%% The model: the CoreStationaryGeneralEqm model on smaller grids
n_d=11;    % hours of labor
n_a=201;   % assets
n_z=5;     % productivity shock
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
z_grid=exp(z_grid);
pi_z_stat=pi_z^1000; pi_z_stat=pi_z_stat(1,:)';
z_grid=z_grid./sum(z_grid.*pi_z_stat); % (the ExogShockFn of parts 5 and 11 does the same normalisation)
d_grid=linspace(0,1,n_d)';
a_grid=30*(linspace(0,1,n_a).^3)';

Params.beta=0.96; % discount factor (calibrated)
DiscountFactorParamNames={'beta'};
Params.sigma=2;   % CRRA on consumption
Params.eta=2;     % curvature of disutility of labor (Frisch=1/eta)
Params.varphi=1;  % weight on disutility of labor (calibrated)
Params.alpha=0.36; % capital share
Params.delta=0.05; % depreciation
Params.A=1;        % TFP
Params.tau=0.2; % labor income tax rate (fixed; the transfer Tr adjusts to balance this budget)
Params.G=0.05;  % exogenous government spending (fixed; the consumption tax tau_c adjusts to fund it)
Params.r=0.04;    % GE price (initial guess; each part solves the general eqm at the truth first)
Params.Tr=0.2;    % GE price
Params.tau_c=0.1; % GE price
Params.rho_z=0.9;    % (only used by the ExogShockFn of parts 5 and 11)
Params.sigma_z=0.03;
GEPriceParamNames={'r','Tr','tau_c'};

% InfHorz (BIHA): 101 hours points, not 11. Hours are a discrete choice (the grid interpolation layer is on assets only), so the
% general eqm conditions move in steps; on 11 points fminsearch stopped the InfHorz truth equilibrium at GE conditions 7.5e-3
% (2026-10-13; the FHorz one reached 1.4e-6), which every near-zero check of parts 8-11 is measured against. On 51
% points it reached 4.8e-4, too close to the 1e-3 bar for the nested inner solves of part 9 (which start 10% off).
n_d_InfHorz=101;
d_grid_InfHorz=linspace(0,1,n_d_InfHorz)';
% and a standard Aiyagari income process, sigma_z=0.2 (the FHorz parts keep 0.03). With sigma_z=0.03 there is almost no precautionary
% saving, so the InfHorz equilibrium r sat 3e-4 below 1/beta-1, on the near-vertical stretch of asset supply: aggregate K went from
% 1.5 to 9.3 as r moved by 5e-4, and to the grid top for r+0.01. Joint mode probes r and beta off the equilibrium (lsqnonlin's
% 0.01 finite-difference step, fminsearch's 5% simplex), so every probe hit that wall (K at the grid top, or beta>1 and no value
% function): all the joint round trips failed while nested worked, and the nested ones were noisy (2026-10-15 verbose diagnostic run).
Params_InfHorz=Params;
Params_InfHorz.sigma_z=0.2; % (the ExogShockFn of part 11 reads it)
[z_grid_InfHorz,pi_z_InfHorz]=discretizeAR1_FarmerToda(0,Params_InfHorz.rho_z,Params_InfHorz.sigma_z,n_z);
z_grid_InfHorz=exp(z_grid_InfHorz);
pi_z_stat_InfHorz=pi_z_InfHorz^1000; pi_z_stat_InfHorz=pi_z_stat_InfHorz(1,:)';
z_grid_InfHorz=z_grid_InfHorz./sum(z_grid_InfHorz.*pi_z_stat_InfHorz); % (the ExogShockFn of part 11 does the same normalisation)

% FHorz (OLG)
N_j=21;
Params.agej=1:1:N_j;
Params.kappa_j=1+0.3*(Params.agej-1)/(N_j-1)-0.4*((Params.agej-1)/(N_j-1)).^2;
Params.kappa_j=Params.kappa_j/mean(Params.kappa_j);
Params.mewj=ones(1,N_j)/N_j;
AgeWeightParamNames={'mewj'};
jequaloneDist=zeros([n_a,n_z],'gpuArray'); % everyone starts with zero assets, z from its stationary distribution
jequaloneDist(1,:)=shiftdim(pi_z_stat,-1);

vfoptions=struct();
vfoptions.gridinterplayer=1; % finite differences (lsqnonlin) need moments that move smoothly with the parameters
vfoptions.ngridinterp=15;
simoptions=struct();
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;

%% ===== doPart(1): FHorz moment vector =====
if doPart(1)==1
    fprintf('\n===== doPart(1): FHorz (OLG) objective moment vector and GE conditions vs hand-built =====\n')
    %% 1. Solve the general eqm at the truth. Targets from full AllStats, LifeCycleProfiles, AutoCorrTransProbs and CrossSectionCovarCorr
    % runs at that equilibrium (every statistic kind, a conditional restriction, NaN-omitted ages). CalibrateOLGModel_Joint_objectivefn,
    % called as CalibrateOLGModel calls it, with vectoroutput=1 at (truth, equilibrium prices): the moment block must equal the hand-built
    % vector and the targets, the GE block must equal HeteroAgentStationaryEqm_Case1_FHorz with maxiter=0 at the same prices; under
    % caliboptions.whichcombos=1 and =0.
    output=CalibSGE_FHorz_MomentVector(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions);
end % doPart(1)

%% ===== doPart(2): FHorz objective forms =====
if doPart(2)==1
    fprintf('\n===== doPart(2): FHorz (OLG) objective forms, scaling, multiGEweights, logmoments, verbose =====\n')
    %% 2. With perturbed targets and non-trivial multiGEweights and relativeGEweight, the joint objective in each form (sum_squared,
    % sum_logratiosquared, MethodOfMoments with a weight matrix; scalar and lsqnonlin residual vector) against its formula on the
    % objective's own moment and GE vectors; the targets part divided by the number of calibrated (non-GE) parameters, and equal to the
    % nested objective at the same prices. Through CalibrateOLGModel: multiGEweights as a column gives exactly what the row gives; the
    % four forms of caliboptions.logmoments become the expected vector; one evaluation with verbose=1 in each of the scalar and vector
    % forms (the print was broken).
    output=CalibSGE_FHorz_ObjectiveForms(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions);
end % doPart(2)

%% ===== doPart(3): FHorz joint vs nested =====
if doPart(3)==1
    fprintf('\n===== doPart(3): FHorz (OLG) joint vs nested at the same prices; the inner GE solve; GE prices in CalibParams =====\n')
    %% 3. The nested objective at the truth: its inner GE prices against the equilibrium (near), its moments against hand-built at its
    % own prices (exact), and the joint objective at those prices gives the same moments and the same GE conditions (exact). Constraints
    % on the GE prices reach the inner solver: 'log' and 'softplus' on Tr both find the equilibrium, but not along the same path (should
    % NOT be zero). CalibrateOLGModel in nested mode returns the GE prices in CalibParams. Joint mode with a 'log' constrainpositive GE
    % price alongside 'softplus' calibrated parameters runs (it indexed past the end of the softplus flags).
    output=CalibSGE_FHorz_JointNested(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions);
end % doPart(3)

%% ===== doPart(4): FHorz round trips =====
if doPart(4)==1
    fprintf('\n===== doPart(4): FHorz (OLG) round trips through CalibrateOLGModel =====\n')
    %% 4. Targets at the truth. Joint (three prices, 51 hours points: on 11 the stepped objective stalls both optimisers) from a wrong
    % start with lsqnonlin; joint with constraints on both the calibrated parameters and the GE prices (constrain0to1, constrainpositive,
    % constrainAtoB); omitcalibparam on part of kappa_j with a row mask. One-price model: nested with lsqnonlin and fminsearch, joint with
    % lsqnonlin and fminsearch, and joint vs nested. Each must recover the truth
    % and its equilibrium (close), and a fresh GE solve at the calibrated parameters must reproduce the reported GE prices (near).
    output=CalibSGE_FHorz_RoundTrips(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions);
end % doPart(4)

%% ===== doPart(5): FHorz shocks, CustomModelStats, intermediateEqns, jequaloneDist as a function =====
if doPart(5)==1
    fprintf('\n===== doPart(5): FHorz (OLG) ExogShockFn, CustomModelStats, intermediateEqns, jequaloneDist as a function =====\n')
    output=CalibSGE_FHorz_Extras(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions);
end % doPart(5)

%% ===== doPart(6): FHorz runtime =====
if doPart(6)==1
    fprintf('\n===== doPart(6): FHorz (OLG) runtime of one objective evaluation (informational) =====\n')
    output=CalibSGE_FHorz_Runtime(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions);
end % doPart(6)

%% ===== doPart(7): InfHorz moment vector =====
if doPart(7)==1
    fprintf('\n===== doPart(7): InfHorz (BIHA) objective moment vector and GE conditions vs hand-built =====\n')
    output=CalibSGE_InfHorz_MomentVector(n_d_InfHorz,n_a,n_z,d_grid_InfHorz,a_grid,z_grid_InfHorz,pi_z_InfHorz,Params_InfHorz,DiscountFactorParamNames,GEPriceParamNames,vfoptions,simoptions);
end % doPart(7)

%% ===== doPart(8): InfHorz objective forms =====
if doPart(8)==1
    fprintf('\n===== doPart(8): InfHorz (BIHA) objective forms, scaling, multiGEweights, logmoments, verbose =====\n')
    output=CalibSGE_InfHorz_ObjectiveForms(n_d_InfHorz,n_a,n_z,d_grid_InfHorz,a_grid,z_grid_InfHorz,pi_z_InfHorz,Params_InfHorz,DiscountFactorParamNames,GEPriceParamNames,vfoptions,simoptions);
end % doPart(8)

%% ===== doPart(9): InfHorz joint vs nested =====
if doPart(9)==1
    fprintf('\n===== doPart(9): InfHorz (BIHA) joint vs nested at the same prices; the inner GE solve; GE prices in CalibParams =====\n')
    output=CalibSGE_InfHorz_JointNested(n_d_InfHorz,n_a,n_z,d_grid_InfHorz,a_grid,z_grid_InfHorz,pi_z_InfHorz,Params_InfHorz,DiscountFactorParamNames,GEPriceParamNames,vfoptions,simoptions);
end % doPart(9)

%% ===== doPart(10): InfHorz round trips =====
if doPart(10)==1
    fprintf('\n===== doPart(10): InfHorz (BIHA) round trips through CalibrateBIHAModel =====\n')
    output=CalibSGE_InfHorz_RoundTrips(n_d_InfHorz,n_a,n_z,d_grid_InfHorz,a_grid,z_grid_InfHorz,pi_z_InfHorz,Params_InfHorz,DiscountFactorParamNames,GEPriceParamNames,vfoptions,simoptions);
end % doPart(10)

%% ===== doPart(11): InfHorz shocks, CustomModelStats, intermediateEqns =====
if doPart(11)==1
    fprintf('\n===== doPart(11): InfHorz (BIHA) ExogShockFn, CustomModelStats, intermediateEqns =====\n')
    output=CalibSGE_InfHorz_Extras(n_d_InfHorz,n_a,n_z,d_grid_InfHorz,a_grid,z_grid_InfHorz,pi_z_InfHorz,Params_InfHorz,DiscountFactorParamNames,GEPriceParamNames,vfoptions,simoptions);
end % doPart(11)

% All looks good!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CalibStationaryGeneralEqmTestsdiary.txt')

diary off
