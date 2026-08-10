% Setup for the Core InfHorz VFI-algorithm tests
%
% This bank tests ONLY ValueFnIter_InfHorz, checking that the different
% value-function-iteration algorithm options give the same V and Policy.
% Uses the same d, a, z, e in every model that needs them, so results are comparable.
%
% Differences from the main CoreInfHorz setup:
%   - adds an iid shock e (for the with-e models)
%   - n_a_big=1500 (used for the big-grid with/without-GI comparison)

n_d=9;
n_a=101;
n_a_big=1500; % for the big-grid 'with vs without GI give very similar V' test
n_z=5;
n_e=3;

d_grid=linspace(0,1,n_d)';
a_grid=5*linspace(0,1,n_a)'.^3;
a_grid_big=5*linspace(0,1,n_a_big)'.^3; % same grid shape, just many more points

% setup z (markov)
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
z_grid=exp(z_grid);

% setup e (iid). Passed to the value fn iteration via vfoptions.n_e/e_grid/pi_e
[e_grid,pi_e]=discretizeAR1_FarmerToda(0,0,0.1,n_e);
pi_e=pi_e(1,:)';
e_grid=exp(e_grid);
vfoptionsbaseline.n_e=n_e;
vfoptionsbaseline.e_grid=e_grid;
vfoptionsbaseline.pi_e=pi_e;
% (simoptions not really needed here since this bank does not simulate, but keep for consistency)
simoptionsbaseline.n_e=n_e;
simoptionsbaseline.e_grid=e_grid;
simoptionsbaseline.pi_e=pi_e;

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
