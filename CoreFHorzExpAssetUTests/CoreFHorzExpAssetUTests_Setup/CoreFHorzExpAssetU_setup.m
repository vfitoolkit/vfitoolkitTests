% Setup so that use the same d,a,z,e,semiz in all the models that use them

% I put everything in vfoptionsbaseline and simoptionsbaseline, so can then
% copy out of these later for each case.

n_d1=7; % labour supply
n_d2=3; % d2 decision for experience asset
n_d3=2; % d3 for semiz
n_d_withoutd1=n_d2;
n_d_withd1=[n_d1,n_d2];
n_d_withoutd1semiz=[n_d2; n_d3];
n_d_withd1semiz=[n_d1,n_d2,n_d3];
n_a_justexpasset=13;
n_a=[101,n_a_justexpasset];
n_a_big=[1001,n_a_justexpasset]; % to test Grid Interpolation
n_z=5;
n_semiz=2; % hardcoded into SemiExoStateFn
n_e=3;

N_j=20;

d1_grid=linspace(0,1,n_d1)'; % d1, labour supply
d2_grid=linspace(0,1,n_d2)'; % d2 for the 
d3_grid=[0;1]; % d3 for semiz models
d_grid_withoutd1=d2_grid;
d_grid_withd1=[d1_grid; d2_grid];
d_grid_withoutd1semiz=[d2_grid; d3_grid];
d_grid_withd1semiz=[d1_grid; d2_grid; d3_grid];

a1_grid=5*linspace(0,1,n_a(1))'.^3;
a1_grid_big=5*linspace(0,1,n_a_big(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a2_grid=linspace(0,10,n_a(2))'; % the experience asset
a_grid=[a1_grid;a2_grid];
a_grid_big=[a1_grid_big;a2_grid];
a_grid_justexpasset=a2_grid;

% with2A1: a genuine SECOND standard endogenous asset a1_2 (multi-point), inserted between the
% liquid asset a1 and the experience asset a2. Two standard assets -> length(n_a1)>1 -> the
% dispatcher routes to DC2A / GI2A / DC2A_GI2A. Grid layout: a = [a1, a1_2, a2].
n_a1_2=4; % multi-point second standard asset (a genuine grid, not binary)
a1_2_grid=2*linspace(0,1,n_a1_2)'.^2; % capped at 2; return r2>r so it gets used up to the cap
n_a_2A1=[n_a(1),n_a1_2,n_a_justexpasset];                       % [a1 (divide-conquered), a1_2 (folded), a2 (experience)]
n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];                 % 'big' a1 grid for the grid-interp moment comparison (smaller than ExpAsset's 201: ExpAssetU carries an extra u dimension in the VFI, so trim to avoid OOM; tune per-variant if needed)
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1=[a1_grid;a1_2_grid;a2_grid];
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

% with2A2: a SECOND EXPERIENCE asset a2_2, so vfoptions.experienceassetu is the integer COUNT of
% a2 dimensions (2 here) rather than a flag. Grid layout: a = [a1, a2_1, a2_2], and the noa1 tier
% is just a = [a2_1, a2_2].
n_a2_1=7; % first experience asset  (human capital, the one the existing tier has)
n_a2_2=5; % second experience asset (cumulated lifetime earnings)
a2_1_grid=linspace(0,10,n_a2_1)';
a2_2_grid=linspace(0,10,n_a2_2)';
n_a_2A2=[51,n_a2_1,n_a2_2]; % base a1 grid; the driver's 'notsobig' grid is FINER in a1
a_grid_2A2=[5*linspace(0,1,n_a_2A2(1))'.^3;a2_1_grid;a2_2_grid];
n_a_2A2_justexpasset=[n_a2_1,n_a2_2];          % noa1 tier: the two experience assets are all there is
a_grid_2A2_justexpasset=[a2_1_grid;a2_2_grid];
% NOTE on sizing: this family carries an extra u dimension through the VFI, so the with2A2 arrays
% are N_u times those of the plain-ExpAsset with2A2 tier at the same a1 grid. 51/151/101 matches
% the ExpAsset(e) with2A2 tier and fits here, but if a variant out-of-memories, trim its a1 grid
% rather than the a2 grids -- the a2 dims are what the tier exists to exercise.


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
vfoptionsbaseline.SemiExoStateFn=@(n,nprime,dsemiz,probfindjob,problosejob) CoreFHorzExpAssetSetup_SemiExoStateFn(n,nprime,dsemiz,probfindjob,problosejob);

% We also need to tell simoptions about the semi-exogenous states
simoptionsbaseline.n_semiz=vfoptionsbaseline.n_semiz;
simoptionsbaseline.semiz_grid=vfoptionsbaseline.semiz_grid;
simoptionsbaseline.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

%% Experience Asset U
vfoptionsbaseline.experienceassetu=1;
vfoptionsbaseline.aprimeFn=@(d2,a2,u,phi1,phi2) aprimeFn_CoreTestExpAssetU(d2,a2,u,phi1,phi2);
% a2prime=u*(phi1*(1-d2)+(1-phi2)*a2);

n_u=7;
Params.stddev_u=0.1;
[u_grid,pi_u]=discretizeAR1_FarmerToda(0,0,Params.stddev_u,n_u);
pi_u=pi_u(1,:)'; % iid
u_grid=exp(u_grid);

vfoptionsbaseline.n_u=n_u;
vfoptionsbaseline.u_grid=u_grid;
vfoptionsbaseline.pi_u=pi_u;
simoptionsbaseline.n_u=vfoptionsbaseline.n_u;
simoptionsbaseline.u_grid=vfoptionsbaseline.u_grid;
simoptionsbaseline.pi_u=vfoptionsbaseline.pi_u;



simoptionsbaseline.experienceassetu=vfoptionsbaseline.experienceassetu;
simoptionsbaseline.aprimeFn=vfoptionsbaseline.aprimeFn;

% with2A2 aprimeFn. GPU arrayfun gives scalar output only, so the builders call aprimeFn once per
% a2 dimension with a 'whicha' integer selector. NOTE the argument order: for experienceassetU the
% selector sits AFTER u and before the parameters, i.e. aprimeFn(d2,a2_1,a2_2,u,whicha,params...).
% Verified against the 32 arrayfun call sites in the l_a2==2 branch of
% CreateExperienceAssetuFnMatrix and its arity check nargin==l_d+l_a2+l_u+(l_a2>=2)+nparams.
vfoptionsbaseline.aprimeFn_2A2=@(d2,a2_1,a2_2,u,whicha,phi1,phi2,phi3,phi4) ...
    (whicha==1)*(u*(phi1*(1-d2)+(1-phi2)*a2_1)) + ...
    (whicha==2)*(u*(phi3*d2*a2_1+(1-phi4)*a2_2));
% a2_1' = u*(phi1*(1-d2)+(1-phi2)*a2_1)  -- EXACTLY the law of motion of the one-experience-asset
%         tier (aprimeFn_CoreTestExpAssetU), so the inert-second-asset cross-test reduces to that
%         tier bit-for-bit.
% a2_2' = u*(phi3*d2*a2_1+(1-phi4)*a2_2) -- deliberately COUPLED to a2_1, and carrying the same u
%         shock, so it is not a payoff-irrelevant passenger and a stride bug cannot cancel out.
simoptionsbaseline.aprimeFn_2A2=vfoptionsbaseline.aprimeFn_2A2;
% simoptionsbaseline.a_grid=
% simoptionsbaseline.d_grid=

Params.phi1=0.3; % ratio at which (1-d2) is converted into a2
Params.phi2=0.03; % depreciation rate
Params.phi3=0.2; % rate at which d2*a2_1 accumulates into a2_2 (with2A2 only)
Params.phi4=0.2; % depreciation of a2_2 (with2A2 only)
Params.pensionrate=0.1; % a2_2 raises the retirement pension -- this is what makes V genuinely
                        % depend on a2_2, rather than it being a payoff-irrelevant passenger.

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
Params.r2=0.06; % return on the second standard asset a1_2 (with2A1 tests); slightly higher than r

% Retirement
Params.Jr=16;
Params.pension=0.5;
Params.agej=1:1:N_j;

% Earings
Params.kappa_j=[0.5:0.1:1,ones(1,9),zeros(1,5)];

% When using semiz
Params.uempbenefit=0.2;
Params.transfer=0.1; % experienceassetu noa1 tiers only: keeps working-age c>0 at a=0, where
                     % earnings are exactly 0 because a is the unconsumable experience asset.
                     % Small relative to pension (0.5) and peak earnings, so it does not swamp
                     % the variation the checks rely on.
Params.searcheffortcost=0.6;
Params.probfindjob=0.7;
Params.problosejob=0.3;
