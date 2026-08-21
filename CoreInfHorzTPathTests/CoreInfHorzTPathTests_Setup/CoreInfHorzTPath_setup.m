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

% Finer grids used ONLY by the general eqm block of the 1A subcodes. Aggregate labour
% N=E[d*z] and aggregate assets K are step functions of the prices when d and a are coarse,
% which makes the GE objective piecewise-constant and stalls fminsearch. n_d=9 in particular
% was giving a CapitalMarket residual ~170x worse than the nod variant of the same model.
% Everything outside the GE block keeps the grids above, so those checks stay comparable.
n_d_GE=51;
n_a_GE=201;
d_grid_GE=linspace(0,1,n_d_GE)';
a_grid_GE=5*linspace(0,1,n_a_GE)'.^3;
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


%% Two-endogenous-state (2A) model, used by the _with2A subcodes
% Simplified Kitao (2008), no taxes: the same model CoreInfHorzTests uses for its 2A subcodes.
%   Endogenous states: asset, occupation (0=worker, 1=entrepreneur)
%   Exogenous markov z: eta (labor productivity), theta (entrepreneurial ability)
% n_a_2A(2)=2 is the occupation state, which is what triggers the DC2A/GI2A code paths.
n_asset_2A=101;
n_asset_2A_big=501; % for the with/without grid interp moment check
n_occ_2A=2;
n_a_2A=[n_asset_2A,n_occ_2A];
n_a_2A_big=[n_asset_2A_big,n_occ_2A];

assetmax_2A=60;
asset_grid_2A=assetmax_2A*(linspace(0,1,n_asset_2A).^3)'; % cubic-spaced, puts points near zero
% Finer asset grid used ONLY by the general eqm block of the 2A subcodes (same reason as the
% 1A ones above). n_l_2A=51 is already fine enough, so the labour grid is left alone.
n_asset_2A_GE=201;
asset_grid_2A_GE=assetmax_2A*(linspace(0,1,n_asset_2A_GE).^3)';
asset_grid_2A_big=assetmax_2A*(linspace(0,1,n_asset_2A_big).^3)';
occ_grid_2A=[0;1]; % 0=worker, 1=entrepreneur
a_grid_2A=[asset_grid_2A; occ_grid_2A];
n_a_2A_GE=[n_asset_2A_GE,n_occ_2A];      % GE-block grid for the 2A subcodes
a_grid_2A_GE=[asset_grid_2A_GE; occ_grid_2A];
a_grid_2A_big=[asset_grid_2A_big; occ_grid_2A];

n_eta_2A=3;   % labor productivity
n_theta_2A=2; % entrepreneurial ability
n_z_2A=[n_eta_2A,n_theta_2A];
[eta_grid_2A,pi_eta_2A]=discretizeAR1_FarmerToda(0,0.9,0.2,n_eta_2A);
eta_grid_2A=exp(eta_grid_2A);
theta_grid_2A=[0; 1.5];             % theta=0 => no production
pi_theta_2A=[0.9,0.1; 0.2,0.8];     % persistent
z_grid_2A=[eta_grid_2A; theta_grid_2A];
pi_z_2A=kron(pi_theta_2A,pi_eta_2A); % eta varies fastest (matches n_z_2A=[n_eta_2A,n_theta_2A])

% Kitao parameters (no taxes). Params.sigma and Params.beta are shared with the 1A subcodes.
Params.alpha=0.36;
Params.delta=0.06;
Params.upsilon=0.88;
Params.upsilon1=Params.alpha*Params.upsilon;
Params.upsilon2=(1-Params.alpha)*Params.upsilon;
Params.leverage=0.5; % max borrowing leverage
Params.phi=0.05;     % extra borrowing cost

% Decision variable for the _with2A subcode that has one: labor supply (Bruggemann 2021).
% Only the d variant uses n_d_2A/d_grid_2A; the nod variant sets n_d=0 itself.
n_l_2A=51;
n_d_2A=n_l_2A;
l_grid_2A=linspace(0,1.5,n_l_2A)';
d_grid_2A=l_grid_2A;
Params.xi=0.716;         % weight on disutility of labor
Params.sigma2=1.7;       % inverse Frisch elasticity
Params.lbar=l_grid_2A(21); % fixed labor supply for entrepreneurs (=0.6, and exactly a grid point,
                           % which matters: the ReturnFn gives -Inf unless abs(l-lbar)<1e-9)

% Same shape and length as the 1A path, but centred on the entrepreneur model's calibration
% (r=0.04, w=1.4) so the model stays well behaved along the path. As with the 1A path, the
% terminal values are what the subcode uses for the period-0 (terminal) VFI.
PricePath_2A.r=[linspace(0.03,0.04,floor(T/2)),0.04*ones(1,T-floor(T/2))];
ParamPath_2A.w=[linspace(1.5,1.4,floor(T/2)),1.4*ones(1,T-floor(T/2))];
