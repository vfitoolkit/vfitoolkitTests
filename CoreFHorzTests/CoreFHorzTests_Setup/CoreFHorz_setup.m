% Setup so that use the same d,a,z,e,semiz in all the models that use them

% I put everything in vfoptionsbaseline and simoptionsbaseline, so can then
% copy out of these later for each case.

n_d=9;
n_d_semiz=[9,2]; % n_d for semiz models with d1
n_d2_semiz=2; % n_d for semiz models without d1
n_a=101;
n_a_big=1001; % to test Grid Interpolation
n_z=5;
n_semiz=2; % hardcoded into SemiExoStateFn
n_e=3;

N_j=20;

d_grid=linspace(0,1,n_d)';
d_grid_semiz=[linspace(0,1,n_d)'; 0;1]; % n_d for semiz models, binary d2, with d1
d2_grid_semiz=[0;1]; % n_d for semiz models, binary d2, without d1

a_grid=5*linspace(0,1,n_a)'.^3;
a_grid_big=5*linspace(0,1,n_a_big)'.^3; % to test Grid Interpolation (same grid, just more points)

% with2A: a genuine SECOND standard endogenous asset (4 grid points) is set up in the
% with2A section of CoreFHorzTests.m (n_a_2A=[n_a,4], phi1/phi2 preferences). Having
% TWO standard endogenous states triggers the DC2A/GI2A/DC2A_GI2A code paths.

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

% setup semiz
vfoptionsbaseline.n_semiz=n_semiz;
vfoptionsbaseline.semiz_grid=[0; 1]; % interpretation: 1 is employed, 0 is not-employed
vfoptionsbaseline.SemiExoStateFn=@(n,nprime,dsemiz,probfindjob,problosejob) CoreFHorzSetup_SemiExoStateFn(n,nprime,dsemiz,probfindjob,problosejob);

% We also need to tell simoptions about the semi-exogenous states
simoptionsbaseline.n_semiz=vfoptionsbaseline.n_semiz;
simoptionsbaseline.semiz_grid=vfoptionsbaseline.semiz_grid;
simoptionsbaseline.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptionsbaseline.n_semiz=vfoptionsbaseline.n_semiz;

%% Now some parameters that models use

Params.beta=0.95; % discount factor
DiscountFactorParamNames={'beta'};

Params.mewj=ones(1,N_j)/N_j;
AgeWeightParamNames={'mewj'};

% Preferences
Params.sigma=2; % CES utility param for consumption
Params.eta=1.5; % curvature of leisure
Params.varphi=0.8; % relative weight of leisure in utility

% Prices
Params.w=1;
Params.r=0.05;

% Retirement
Params.Jr=16;
Params.pension=0.5;
Params.agej=1:1:N_j;

% Earings
Params.kappa_j=[0.5:0.1:1,ones(1,9),zeros(1,5)];

% When using semiz
Params.uempbenefit=0.2;
Params.searcheffortcost=0.6;
Params.probfindjob=0.7;
Params.problosejob=0.3;
