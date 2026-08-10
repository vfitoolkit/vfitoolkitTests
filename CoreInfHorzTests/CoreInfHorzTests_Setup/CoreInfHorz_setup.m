% Setup for the Core InfHorz tests
%
% Uses the same d, a, z in every model that needs them, so results are comparable.
%
% Differences from the FHorz setup:
%   - No N_j, no AgeWeightParamNames, no agej/Jr/pension, no kappa_j age profile
%   - Discount: just beta (no age-dependent factor)
%   - No e/semiz (this bank only sweeps with/without d and with/without z)
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

d_grid=linspace(0,1,n_d)';
a_grid=5*linspace(0,1,n_a)'.^3;
a_grid_big=5*linspace(0,1,n_a_big)'.^3; % to test Grid Interpolation (same grid, just more points)

% setup z
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
z_grid=exp(z_grid);

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

% vfoptions/simoptions baselines (nothing e/semiz specific here, but keep the
% same calling convention as the other test banks)
vfoptionsbaseline=struct();
simoptionsbaseline=struct();
