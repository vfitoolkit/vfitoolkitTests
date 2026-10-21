% Setup for the Core InfHorz tests
%
% Uses the same d, a, z in every model that needs them, so results are comparable.
%
% Differences from the FHorz setup:
%   - No N_j, no AgeWeightParamNames, no agej/Jr/pension, no kappa_j age profile
%   - Discount: just beta (no age-dependent factor)
%   - semiz (the semi-exogenous state, whose transitions depend on the last decision variable d2) is as in the FHorz
%     setup: n_semiz=2 employed/not-employed, a binary d2 (search effort), and an optional d1 (labor supply)
%
% Note on the noz cases: with no exogenous shock the InfHorz stationary
% distribution collapses to a single mass point, so the moment/autocorrelation/
% cross-section/time-series statistics are degenerate. For those cases the stats
% commands are run as 'does it execute' checks; the z cases provide the
% meaningful cross-validation of the statistics.

n_d=9;
n_a=101;
n_a_big=1001; % to test Grid Interpolation
n_z=5;
n_e=3;
n_semiz=2; % hardcoded into CoreInfHorzSetup_SemiExoStateFn
n_d_semiz=[9,2]; % n_d for semiz models with d1
n_d2_semiz=2; % n_d for semiz models without d1

d_grid=linspace(0,1,n_d)';
d_grid_semiz=[linspace(0,1,n_d)'; 0;1]; % d grid for semiz models, binary d2, with d1
d2_grid_semiz=[0;1]; % d grid for semiz models, binary d2, without d1
a_grid=5*linspace(0,1,n_a)'.^3;
a_grid_big=5*linspace(0,1,n_a_big)'.^3; % to test Grid Interpolation (same grid, just more points)

% setup z
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
z_grid=exp(z_grid);

% setup e (iid). Passed to the value fn iteration via vfoptions.n_e/e_grid/pi_e (and to the
% distribution/evaluation commands via the same three simoptions fields)
[e_grid,pi_e]=discretizeAR1_FarmerToda(0,0,0.1,n_e);
pi_e=pi_e(1,:)'; % iid, so all rows are the same; keep one, as a column
e_grid=exp(e_grid);

%% Parameters
Params.beta=0.95; % discount factor
DiscountFactorParamNames={'beta'};

% Preferences
Params.sigma=2; % CES utility param for consumption
Params.eta=1.5; % curvature of leisure
Params.varphi=0.8; % relative weight of leisure in utility

% Prices
Params.w=1;
Params.r=0.05;

% semiz (as in the FHorz setup)
Params.uempbenefit=0.2;
Params.searcheffortcost=0.6;
Params.probfindjob=0.7;
Params.problosejob=0.3;

% vfoptions/simoptions baselines. The e and the semiz are carried on the baselines, and the
% subcodes that use them copy them across; the subcodes that do not use them simply ignore them.
vfoptionsbaseline=struct();
vfoptionsbaseline.n_e=n_e;
vfoptionsbaseline.e_grid=e_grid;
vfoptionsbaseline.pi_e=pi_e;
simoptionsbaseline=struct();
simoptionsbaseline.n_e=n_e;
simoptionsbaseline.e_grid=e_grid;
simoptionsbaseline.pi_e=pi_e;
% semiz
vfoptionsbaseline.n_semiz=n_semiz;
vfoptionsbaseline.semiz_grid=[0; 1]; % interpretation: 1 is employed, 0 is not-employed
vfoptionsbaseline.SemiExoStateFn=@(n,nprime,dsemiz,probfindjob,problosejob) CoreInfHorzSetup_SemiExoStateFn(n,nprime,dsemiz,probfindjob,problosejob);
simoptionsbaseline.n_semiz=vfoptionsbaseline.n_semiz;
simoptionsbaseline.semiz_grid=vfoptionsbaseline.semiz_grid;
simoptionsbaseline.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
