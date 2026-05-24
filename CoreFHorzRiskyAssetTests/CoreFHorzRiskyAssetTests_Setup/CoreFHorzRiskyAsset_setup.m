% Setup: shared params, grids, and options for all RiskyAsset test variants.
%
% Decision-variable categories under vfoptions.refine_d (riskyasset convention):
%   d1: in ReturnFn but not aprimeFn (here: labour supply, h)
%   d2: in aprimeFn but not ReturnFn (here: riskyshare)
%   d3: in both ReturnFn and aprimeFn (here: savings)
%   d4: in ReturnFn but not aprimeFn, AND determines semiz transitions (only with semiz)

n_d1=3;  % labour supply
n_d2=3;  % riskyshare
n_d3=11; % savings
n_d4=2;  % semiz decision (only used in semiz variants)

% Composite n_d for the four variant patterns
n_d_withoutd1=[n_d2,n_d3];
n_d_withd1=[n_d1,n_d2,n_d3];
n_d_withoutd1semiz=[n_d2,n_d3,n_d4];
n_d_withd1semiz=[n_d1,n_d2,n_d3,n_d4];

n_a=51;     % single endogenous state (the riskyasset itself)
n_a_big=201; % for more accurate StationaryDist moments

n_z=5;
n_e=3;
n_u=3;     % iid between-period shock to risky-asset return
n_semiz=2; % hardcoded into SemiExoStateFn

N_j=20;

% d grids
d1_grid=linspace(0,1,n_d1)';
d2_grid=linspace(0,1,n_d2)';
d3_grid=10*linspace(0,1,n_d3)'.^3; % savings, cubic spacing for finer near zero
d4_grid=[0;1];                     % semiz decision, binary

d_grid_withoutd1=[d2_grid; d3_grid];
d_grid_withd1=[d1_grid; d2_grid; d3_grid];
d_grid_withoutd1semiz=[d2_grid; d3_grid; d4_grid];
d_grid_withd1semiz=[d1_grid; d2_grid; d3_grid; d4_grid];

% Asset grid (single endogenous state, the riskyasset)
a_grid=10*linspace(0,1,n_a)'.^3;
a_grid_big=10*linspace(0,1,n_a_big)'.^3;

% setup z
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
z_grid=exp(z_grid);

% setup e (within-period iid)
[e_grid,pi_e]=discretizeAR1_FarmerToda(0,0,0.1,n_e);
pi_e=pi_e(1,:)';
e_grid=exp(e_grid);
vfoptionsbaseline.n_e=n_e;
vfoptionsbaseline.e_grid=e_grid;
vfoptionsbaseline.pi_e=pi_e;
simoptionsbaseline.n_e=vfoptionsbaseline.n_e;
simoptionsbaseline.e_grid=vfoptionsbaseline.e_grid;
simoptionsbaseline.pi_e=vfoptionsbaseline.pi_e;

% setup u (between-period iid shock to risky asset return)
Params.rp=0.03;       % mean excess return
Params.sigma_u=0.025; % stdev of innovations
[u_grid,pi_u]=discretizeAR1_FarmerToda(Params.rp,0,Params.sigma_u,n_u);
pi_u=pi_u(1,:)'; % iid
vfoptionsbaseline.n_u=n_u;
vfoptionsbaseline.u_grid=u_grid;
vfoptionsbaseline.pi_u=pi_u;
simoptionsbaseline.n_u=vfoptionsbaseline.n_u;
simoptionsbaseline.u_grid=vfoptionsbaseline.u_grid;
simoptionsbaseline.pi_u=vfoptionsbaseline.pi_u;

% setup semiz
vfoptionsbaseline.n_semiz=n_semiz;
vfoptionsbaseline.semiz_grid=[0;1]; % interpretation: 1 is employed, 0 is not-employed
vfoptionsbaseline.SemiExoStateFn=@(n,nprime,dsemiz,probfindjob,problosejob) CoreFHorzRiskyAssetSetup_SemiExoStateFn(n,nprime,dsemiz,probfindjob,problosejob);
simoptionsbaseline.n_semiz=vfoptionsbaseline.n_semiz;
simoptionsbaseline.semiz_grid=vfoptionsbaseline.semiz_grid;
simoptionsbaseline.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

%% Risky asset
vfoptionsbaseline.riskyasset=1;
vfoptionsbaseline.aprimeFn=@(riskyshare,savings,u,r) aprimeFn_CoreTestRiskyAsset(riskyshare,savings,u,r);
% aprime = (1+r)*(1-riskyshare)*savings + (1+r+u)*riskyshare*savings   (when savings>0)
simoptionsbaseline.riskyasset=vfoptionsbaseline.riskyasset;
simoptionsbaseline.aprimeFn=vfoptionsbaseline.aprimeFn;

%% Parameters
Params.beta=0.95;
DiscountFactorParamNames={'beta'};

Params.mewj=ones(1,N_j)/N_j;
AgeWeightParamNames={'mewj'};

% Preferences
Params.sigma=2;
Params.eta=1.5;
Params.varphi=0.8;

% Prices
Params.w=1;
Params.r=0.05;

% Retirement
Params.Jr=16;
Params.pension=0.5;
Params.agej=1:1:N_j;

% Earnings
Params.kappa_j=[0.5:0.1:1,ones(1,9),zeros(1,5)];

% For semiz models
Params.uempbenefit=0.2;
Params.searcheffortcost=0.6;
Params.probfindjob=0.7;
Params.problosejob=0.3;
