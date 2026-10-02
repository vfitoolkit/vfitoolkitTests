% Setup so that use the same d,a,z,e in all the models that use them
% experienceassetze: aprime depends on (d2,a2,z,e)

n_d1=7; % labour supply
n_d2=3; % d2 decision for experience asset
n_d3=2; % d3 decision for semiz
n_d_withoutd1=n_d2;
n_d_withd1=[n_d1,n_d2];
n_d_withoutd1semiz=[n_d2,n_d3];
n_d_withd1semiz=[n_d1,n_d2,n_d3];
n_a_justexpasset=13;
n_a=[101,n_a_justexpasset];
n_a_big=[1001,n_a_justexpasset]; % to test Grid Interpolation
n_z=5;
n_e=3;
n_semiz=2; % hardcoded into SemiExoStateFn

N_j=20;

d1_grid=linspace(0,1,n_d1)'; % d1, labour supply
d2_grid=linspace(0,1,n_d2)'; % d2 for the experience asset
d3_grid=[0;1]; % d3 for semiz
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

% with2A1: a SECOND standard endogenous asset a1_2, inserted between the liquid asset a1 and the
% experience asset a2. Two standard assets -> length(n_a1)>1 -> the dispatcher routes to
% DC2A / GI2A / DC2A_GI2A. Grid layout: a = [a1, a1_2, a2]. Built here rather than inside the
% subcodes, so that n_a_2A1 holds the real 2A1 shape and not a withA1 shape.
n_a1_2=2;
a1_2_grid=[0;1]; % binary second asset (capped high-return asset)
n_a_2A1=[51,n_a1_2,n_a_justexpasset];
a_grid_2A1=[5*linspace(0,1,n_a_2A1(1))'.^3;a1_2_grid;a2_grid];

% with2A2: a SECOND EXPERIENCE asset a2_2, so vfoptions.experienceassetze is the integer COUNT of
% a2 dimensions (2 here) rather than a flag. Grid layout: a = [a1, a2_1, a2_2], and the noa1 tier
% is just a = [a2_1, a2_2].
n_a2_1=7; % first experience asset  (human capital, the one the existing tier has)
n_a2_2=5; % second experience asset (cumulated lifetime earnings)
a2_1_grid=linspace(0,10,n_a2_1)';
a2_2_grid=linspace(0,10,n_a2_2)';
n_a_2A2=[51,n_a2_1,n_a2_2]; % base a1 grid, same 51 as n_a_2A1; the driver's 'notsobig' grid is FINER in a1
a_grid_2A2=[5*linspace(0,1,n_a_2A2(1))'.^3;a2_1_grid;a2_2_grid];
n_a_2A2_justexpasset=[n_a2_1,n_a2_2];          % noa1 tier: the two experience assets are all there is
a_grid_2A2_justexpasset=[a2_1_grid;a2_2_grid];
% NOTE the argument order. For experienceassetZE the 'whicha' selector sits AFTER z and e, i.e.
% aprimeFn(d2,a2_1,a2_2,z,e,whicha,params...). That is NOT where it sits for the plain experience
% asset (there it follows the a2 inputs directly). Verified against the 128 arrayfun call sites in
% the l_a2==2 branch of CreateExperienceAssetzeFnMatrix, and against its arity check
% nargin==l_d+l_a2+l_z+l_e+(l_a2>=2)+nparams. The selector exists because GPU arrayfun is
% scalar-output only, so the builder calls aprimeFn once per a2 dimension.
vfoptionsbaseline.aprimeFn_2A2=@(d2,a2_1,a2_2,z,e,whicha,phi1,phi2,phi3,phi4) ...
    (whicha==1)*(phi1*(1-d2)*z*e+(1-phi2)*a2_1) + ...
    (whicha==2)*(phi3*d2*a2_1+(1-phi4)*a2_2);
% a2_1' = phi1*(1-d2)*z*e+(1-phi2)*a2_1  -- EXACTLY the law of motion of the one-experience-asset
%         tier, so the inert-second-asset cross-test reduces to that tier bit-for-bit.
% a2_2' = phi3*d2*a2_1+(1-phi4)*a2_2     -- deliberately COUPLED to a2_1. With an uncoupled second
%         asset a stride bug can cancel out and the checks go vacuous.
Params.phi3=0.2; % rate at which d2*a2_1 accumulates into a2_2
Params.phi4=0.2; % depreciation of a2_2
Params.pensionrate=0.1; % a2_2 raises the retirement pension -- this is what makes V genuinely
                        % depend on a2_2, rather than it being a payoff-irrelevant passenger.



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

%% Semiz
vfoptionsbaseline.n_semiz=n_semiz;
vfoptionsbaseline.semiz_grid=[0; 1]; % interpretation: 1 is employed, 0 is not-employed
vfoptionsbaseline.SemiExoStateFn=@(n,nprime,dsemiz,probfindjob,problosejob) CoreFHorzExpAssetzeSetup_SemiExoStateFn(n,nprime,dsemiz,probfindjob,problosejob);

simoptionsbaseline.n_semiz=vfoptionsbaseline.n_semiz;
simoptionsbaseline.semiz_grid=vfoptionsbaseline.semiz_grid;
simoptionsbaseline.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

%% Experience Asset (ze variant)
vfoptionsbaseline.experienceassetze=1;
vfoptionsbaseline.aprimeFn=@(d2,a2,z,e,phi1,phi2) aprimeFn_CoreTestExpAssetze(d2,a2,z,e,phi1,phi2);
% a2prime=phi1*(1-d2)*z*e+(1-phi2)*a2;


simoptionsbaseline.experienceassetze=vfoptionsbaseline.experienceassetze;
simoptionsbaseline.aprimeFn=vfoptionsbaseline.aprimeFn;

Params.phi1=0.3; % ratio at which (1-d2)*z*e is converted into a2
Params.phi2=0.03; % depreciation rate

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
Params.r2=0.08; % return on the binary asset a1_2 (higher than r, so it is used up to the cap)

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
