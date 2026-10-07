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
%
% Parts are independent: the setup, the addpaths and the model preamble all sit OUTSIDE the if-blocks and so always run, and no
% part reads another part's output. Any subset can be run, in any combination. Anything added to this bank later must keep that true.
doPart=[1,1];

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

% All looks good!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/EstimFHorzTestsdiary.txt')

diary off
