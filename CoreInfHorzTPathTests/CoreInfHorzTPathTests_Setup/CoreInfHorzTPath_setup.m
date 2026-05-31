% Setup for InfHorz TPath tests
%
% Differences from the FHorz setup:
%   - No N_j, no AgeWeightParamNames, no agej/Jr/pension, no kappa_j age profile
%   - Discount: just beta (no age-dependent factor)
%   - No fastOLG (that's an OLG-specific optimization)
%   - Initial agent dist for TPath = stationary distribution at period-0 prices
%     (computed inside each subcode from the period-0 VFI solution)

n_d=9;
n_a=101;
n_a_big=1001; % to test Grid Interpolation
n_z=5;
n_e=3;

d_grid=linspace(0,1,n_d)';
a_grid=5*linspace(0,1,n_a)'.^3;
a_grid_big=5*linspace(0,1,n_a_big)'.^3; % to test Grid Interpolation (same grid, just more points)

% setup z
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
z_grid=exp(z_grid);

% setup e
[e_grid,pi_e]=discretizeAR1_FarmerToda(0,0,0.1,n_e);
pi_e=pi_e(1,:)';
e_grid=exp(e_grid);
vfoptionsbaseline.n_e=n_e;
vfoptionsbaseline.e_grid=e_grid;
vfoptionsbaseline.pi_e=pi_e;

simoptionsbaseline.n_e=vfoptionsbaseline.n_e;
simoptionsbaseline.e_grid=vfoptionsbaseline.e_grid;
simoptionsbaseline.pi_e=vfoptionsbaseline.pi_e;


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


%% Transition path aspects

T=100;

% PricePath, ParamPath
PricePath.r=[linspace(0.03,0.05,floor(T/2)),0.05*ones(1,T-floor(T/2))];
ParamPath.w=[linspace(1.1,1,floor(T/2)),1*ones(1,T-floor(T/2))];

transpathoptionsbaseline=struct();
