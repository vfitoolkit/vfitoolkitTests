% QUASI-HYPERBOLIC discounting tests of the core VFI Toolkit FHorz ExpAssetsemiz commands.
% Mirrors CoreFHorzExpAssetsemizTests.m one-for-one: the same 24 combinations in the same figure
% order -- noa1 1-8, withA1 9-16, with2A1 17-24 -- over
%   with/without d1, with/without z, with/without e (semiz is inherent to this family),
%   noa1 / withA1 / with2A1, and every solution method the baseline runs.
%
% Each subcode runs Naive AND Sophisticated, the full lowmemory ladder on each method, a
% ValueFnFromPolicy oracle, and the exponential cross-tests: (i) the Naive continuation value
% equals the exponential value fn at the actual beta0; (ii)/(iii) with beta0=1 both the main
% value fn and the continuation value equal exponential, for Naive and Sophisticated.
%
% TEST-FIRST: the toolkit currently has NO quasi-hyperbolic support for experienceassetsemiz
% (ValueFnIter_Case1_FHorz has no QH branch for it), so every figure errors at its first
% ValueFnIter call. These were written ahead of the toolkit code, on purpose.
%
% The subcodes draw no figures, so only the diary is saved.
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzQHExpAssetsemizTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzQHExpAssetsemizTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzQHExpAssetsemizTestsdiary.txt

addpath('./CoreFHorzQHExpAssetsemizTests_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzQHExpAssetsemizTests_subcodes/WithA1_subcodes/Semiz_subcodes/')
addpath('../CoreFHorzExpAssetsemizTests_Setup/')
addpath('../CoreFHorzExpAssetsemiz_ReturnFns/')


%% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzExpAssetsemiz_setup

Params.beta0=0.9; % additional today-tomorrow (present-bias) discount factor
vfoptionsbaseline.QHadditionaldiscount={'beta0'};


%% ================= WITHOUT a1 (figs 1-8): experience asset is the only endogenous state =================
% No DC/GI/DC+GI blocks (irrelevant without a1).
% Pass n_a_justexpasset as n_a, a_grid_justexpasset as a_grid. n_a_big/a_grid_big slots unused.

addpath('./CoreFHorzQHExpAssetsemizTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzQHExpAssetsemizTests_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('../CoreFHorzExpAssetsemiz_ReturnFns/Noa1_ReturnFns/')

%% noa1 nosemiz (8 variants)









%% noa1 nosemiz cross-tests
% Markov-as-iid equivalence cross-test
% ExpAssetsemiz noa1 with degenerate aprimeFn(d2,a2)=d2 should match a standard 1-endo model with d2
% no a1 vs model with a1 but where it is ignored


%% noa1 (8 variants)

%% without d1, without z, without e, noa1, semiz
figure_c=1;
output=CoreFHorzQHExpAssetsemiz_nod1_noz_noe_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, without e, noa1, semiz
figure_c=2;
output=CoreFHorzQHExpAssetsemiz_d1_noz_noe_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, without e, noa1, semiz
figure_c=3;
output=CoreFHorzQHExpAssetsemiz_nod1_z_noe_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, without e, noa1, semiz
figure_c=4;
output=CoreFHorzQHExpAssetsemiz_d1_z_noe_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, without z, with e, noa1, semiz
figure_c=5;
output=CoreFHorzQHExpAssetsemiz_nod1_noz_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, with e, noa1, semiz
figure_c=6;
output=CoreFHorzQHExpAssetsemiz_d1_noz_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e, noa1, semiz
figure_c=7;
output=CoreFHorzQHExpAssetsemiz_nod1_z_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, noa1, semiz
figure_c=8;
output=CoreFHorzQHExpAssetsemiz_d1_z_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% noa1 semiz cross-tests
% Markov-as-iid equivalence cross-test
% semiz state as a plain Markov should give same answer
% ExpAssetsemiz noa1 with degenerate aprimeFn(d2,a2)=d2 should match a standard 1-endo model with d2
% no a1 vs model with a1 but where it is ignored




%% ================= WITH a1 (figs 9-16) =================
% Reset the setup (the without-a1 half above left the workspace alone, but re-run for safety/independence)
CoreFHorzExpAssetsemiz_setup

%% a1 nosemiz (8 variants)

% looks good

% looks good

% looks good

%% with d1, with z, without e, without semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

% looks good

% looks good

%% with d1, without z, with e, without semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

% looks good

%% without d1, with z, with e, without semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

% looks good

%% with d1, with z, with e, without semiz
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

% looks good

%% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid

% all looking good :)

%% Do a test with a 'fake experience asset' and compare to a standard endogneous asset

% all looking good :)


%% Worth doing a 'clear all' here, but not necessary.
% Mainly is so you can run second half independent of first half

%% That is all the without semiz, now with semiz
% From here on, it is the eight with semiz
% From here on, use n_d_semiz and d_grid_semiz as the inputs (instead of n_d and d_grid)


% Uses the same setup, which already had a semi-exogenous state, just that it wasn't used.
CoreFHorzExpAssetsemiz_setup

% For models without d1, use:
% n_d2_semiz and d2_grid_semiz (as n_d and d_grid)
% For models with d1, use:
% n_d_semiz and d_grid_semiz (as n_d and d_grid)

%% a1 semiz (8 variants)

%% without d1, without z, without e, with semiz
figure_c=9;
output=CoreFHorzQHExpAssetsemiz_nod1_noz_noe(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, without z, without e, with semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=10;
output=CoreFHorzQHExpAssetsemiz_d1_noz_noe(n_d_withd1,n_a,n_a_notsobig,0,N_j,d_grid_withd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

% looks good (as good as it can be expected to given the n_a_notsobig)

%% without d1, with z, without e, with semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=11;
output=CoreFHorzQHExpAssetsemiz_nod1_z_noe(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

% looks good (as good as it can be expected to given the n_a_notsobig)

%% with d1, with z, without e, with semiz
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=12;
output=CoreFHorzQHExpAssetsemiz_d1_z_noe(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

% looks good (as good as it can be expected to given the n_a_notsobig)

%% without d1, without z, with e, with semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=13;
output=CoreFHorzQHExpAssetsemiz_nod1_noz_e(n_d_withoutd1,n_a,n_a_notsobig,0,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

% looks good (as good as it can be expected to given the n_a_notsobig)

%% with d1, without z, with e, with semiz
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=14;
output=CoreFHorzQHExpAssetsemiz_d1_noz_e(n_d_withd1,n_a,n_a_notsobig,0,N_j,d_grid_withd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

% looks good (as good as it can be expected to given the n_a_notsobig)

%% without d1, with z, with e, with semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=15;
output=CoreFHorzQHExpAssetsemiz_nod1_z_e(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

% looks good

%% with d1, with z, with e, with semiz
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=16;
output=CoreFHorzQHExpAssetsemiz_d1_z_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

% looks good (as good as it can be expected to given the n_a_notsobig)

%% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid

% all looking good :)

%% Now some further cross-tests, using a semi-exo that is really just a markov

% all looking good :)

%% Do a test with a 'fake experience asset' and compare to a standard endogneous asset

% all looking good :)


%% ================= WITH 2a1 (figs 33-48): two standard endogenous assets (triggers DC2A/GI2A/DC2A_GI2A) =================
% A genuine multi-point second standard asset a1_2 (return r2) is spliced between the liquid asset
% a1 and the experience asset a2: n_a_2A1=[n_a1, n_a1_2, n_a2] (built in CoreFHorzExpAssetsemiz_setup).
% Two standard assets -> length(n_a1)>1 -> DC2A / GI2A / DC2A_GI2A.
% The 8 semiz variants (figs 41-48) exercise the ExpAssetSemiExo DC2A/GI2A/DC2A_GI2A family
% (written test-first; the toolkit raws now exist).
CoreFHorzExpAssetsemiz_setup
addpath('./CoreFHorzQHExpAssetsemizTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzQHExpAssetsemizTests_subcodes/With2A1_subcodes/Semiz_subcodes/')
addpath('../CoreFHorzExpAssetsemiz_ReturnFns/With2A1_ReturnFns/')

%% 2a1 nosemiz (8 variants)
n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];


n_a_2A1_notsobig=[101,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];


n_a_2A1_notsobig=[75,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS

%% 2a1 nosemiz cross-tests: a degenerate second asset a1_2 (single point {0}) reduces to the with-a1 model

%% 2a1 semiz (8 variants)
n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

figure_c=17;
output=CoreFHorzQHExpAssetsemiz_nod1_noz_noe_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
figure_c=18;
output=CoreFHorzQHExpAssetsemiz_d1_noz_noe_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

n_a_2A1_notsobig=[101,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

figure_c=19;
output=CoreFHorzQHExpAssetsemiz_nod1_z_noe_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

n_a_2A1_notsobig=[75,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

figure_c=20;
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
output=CoreFHorzQHExpAssetsemiz_d1_z_noe_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
figure_c=21;
output=CoreFHorzQHExpAssetsemiz_nod1_noz_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
figure_c=22;
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
output=CoreFHorzQHExpAssetsemiz_d1_noz_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
figure_c=23;
output=CoreFHorzQHExpAssetsemiz_nod1_z_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
figure_c=24;
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
output=CoreFHorzQHExpAssetsemiz_d1_z_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS


diary off

%% THINGS NOT CHECKED
% Check using two decision variables in any of d1 or d3 (the decision variables that are not in experience asset)
