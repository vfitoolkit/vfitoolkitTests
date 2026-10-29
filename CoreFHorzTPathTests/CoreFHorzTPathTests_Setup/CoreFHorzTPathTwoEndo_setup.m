% Setup for the TWO ENDOGENOUS STATE parts, doPart(12)-(13). Run AFTER CoreFHorzTPath_setup.
%
% These parts used to be a bank of their own, CoreFHorzTPathTwoEndoTests, whose setup was a copy of
% CoreFHorzTPath_setup differing ONLY in what is defined here: the asset grids (two endogenous states, a1
% and a2) and three parameters (the second asset's preference parameters phi1 and phi2, and Jr=11 rather
% than 16). Everything else (d, z, e, semiz, N_j, T, PricePath, ParamPath, the baseline options) is shared.
% So this defines only the differences, under _2A names, and leaves every one-asset variable untouched:
% the parts are independent and the one-asset parts read n_a, a_grid and Params.

n_a_2A=[101,4]; % two endogenous states
n_a_big_2A=[1001,4]; % to test Grid Interpolation

% Two endogenous states
a1_grid_2A=5*linspace(0,1,n_a_2A(1))'.^3;
a1_grid_big_2A=5*linspace(0,1,n_a_big_2A(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a2_grid_2A=[0;1;2;3]; % n_a_2A(2)
a_grid_2A=[a1_grid_2A; a2_grid_2A];
a_grid_big_2A=[a1_grid_big_2A; a2_grid_2A];

% A middle-sized a1 grid for the with d, z and e model, whose big-grid version does not fit in memory
n_a_notsobig_2A=[251,4];
a_grid_notsobig_2A=[a_grid_big_2A(round(linspace(1,n_a_big_2A(1),n_a_notsobig_2A(1)))'); a_grid_big_2A(n_a_big_2A(1)+1:end)];

% Parameters: the one-asset ones, plus preferences for the second endogenous state, and an earlier retirement
Params_2A=Params;
Params_2A.phi1=3;
Params_2A.phi2=0.1;
Params_2A.Jr=11;
