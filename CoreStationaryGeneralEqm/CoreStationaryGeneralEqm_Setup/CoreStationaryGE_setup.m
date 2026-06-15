% Setup for the Core Stationary General Eqm tests.
%
% Baseline model (shared by InfHorz and FHorz halves):
%   Aiyagari-style incomplete markets with endogenous labor and a government.
%   Household state (a,z), decision d=hours, choice aprime.
%     Budget: (1+tau_c)*c+aprime=(1+r)*a+(1-tau)*w*z*d+Tr  (FHorz multiplies earnings by kappa_j)
%   Firm Cobb-Douglas Y=A*K^alpha*N^(1-alpha):
%     r = alpha*A*(K/N)^(alpha-1)-delta ,  w = (1-alpha)*A*(K/N)^alpha
%   The wage w is hardcoded from r via the firm FOC (computed inside the
%   ReturnFn and wherever else it is needed), so it is NOT a GE price.
%   Government: fixed labor tax tau funds a lump-sum transfer Tr (Tr adjusts);
%   a consumption tax tau_c funds exogenous government spending G (tau_c adjusts).
%
% Three GE prices: r, Tr, tau_c
% Three GE eqns:   CapitalMarket (pins r), GovBudget tau*w*N=Tr (pins Tr),
%                  ConsTax tau_c*C=G (pins tau_c)
%
% The equilibrium is interior so that the parameter-constraint tests are
% non-binding (r in (0,1), Tr>0, tau_c in [0,1]).

%% Grids
n_d=101;    % hours of labor
n_a=501;   % assets
n_z=5;     % productivity shock

% setup z (AR(1) labor productivity)
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
z_grid=exp(z_grid);
% normalize so E[z]=1 under the stationary distribution
pi_z_stat=pi_z^1000; pi_z_stat=pi_z_stat(1,:)';
z_grid=z_grid./sum(z_grid.*pi_z_stat);

% hours grid
d_grid=linspace(0,1,n_d)';

%% Parameters
Params.beta=0.96; % discount factor
DiscountFactorParamNames={'beta'};

% Preferences
Params.sigma=2;   % CRRA on consumption
Params.eta=2;     % curvature of disutility of labor (Frisch=1/eta)
Params.varphi=1;  % weight on disutility of labor

% Production
Params.alpha=0.36; % capital share
Params.delta=0.05; % depreciation
Params.A=1;        % TFP

% Government
Params.tau=0.2; % labor income tax rate (fixed; the transfer Tr adjusts to balance this budget)
Params.G=0.05;  % exogenous government spending (fixed; the consumption tax tau_c adjusts to fund it)

% Prices (GE variables) -- initial guesses
Params.r=0.04;
Params.Tr=0.2;    % lump-sum transfer level (GE price: set by Tr=tau*w*N)
Params.tau_c=0.1; % consumption tax rate (GE price: set by tau_c*C=G)

GEPriceParamNames={'r','Tr','tau_c'};

%% Asset grid (scaled off the no-risk steady-state capital)
a_grid=30*(linspace(0,1,n_a).^3)';

%% FHorz-specific setup (ignored by the InfHorz half)
N_j=41; % number of periods (model age)
Params.agej=1:1:N_j;

% age-profile of labor productivity (hump shaped), normalized to mean 1
Params.kappa_j=1+0.3*(Params.agej-1)/(N_j-1)-0.4*((Params.agej-1)/(N_j-1)).^2;
Params.kappa_j=Params.kappa_j/mean(Params.kappa_j);

% age weights: stationary OLG with no population growth and no death
Params.mewj=ones(1,N_j)/N_j;
AgeWeightParamNames={'mewj'};

% newborns: zero assets, productivity drawn from stationary dist of pi_z
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,:)=shiftdim(pi_z_stat,-1);

%% Baseline options
vfoptionsbaseline.gridinterplayer=1;
vfoptionsbaseline.ngridinterp=20;
simoptionsbaseline.gridinterplayer=vfoptionsbaseline.gridinterplayer;
simoptionsbaseline.ngridinterp=vfoptionsbaseline.ngridinterp;

heteroagentoptionsbaseline.verbose=1;
heteroagentoptionsbaseline.toleranceGEprices=1e-5;
heteroagentoptionsbaseline.toleranceGEcondns=1e-5;
