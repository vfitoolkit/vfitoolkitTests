% Setup for CoreFHorzPTypeTests.
% Builds on CoreFHorz_setup (same d, a, z, e, semiz grids and Params), then
% adds PType-specific pieces.

run('../CoreFHorzTests/CoreFHorzTests_Setup/CoreFHorz_setup.m')

%% PType-specific
N_i=2;
Names_i={'low','high'};

PTypeDistParamNames={'ptypeweights'};
Params.ptypeweights=[0.6; 0.4];

% A parameter that genuinely differs across PTypes (used in tests 2 and 4)
% Row ii is the kappa_j profile for PType ii.
Params.kappa_j_pt=[Params.kappa_j; 1.2*Params.kappa_j];
