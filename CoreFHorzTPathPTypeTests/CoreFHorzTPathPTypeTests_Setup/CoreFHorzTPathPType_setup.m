% Setup for CoreFHorzTPathPTypeTests.
% Builds on CoreFHorzTPath_setup (same d, a, z, e grids and Params, plus T,
% PricePath, ParamPath and transpathoptionsbaseline), then adds PType-specific
% pieces.

run('../CoreFHorzTPathTests/CoreFHorzTPathTests_Setup/CoreFHorzTPath_setup.m')

%% PType-specific
N_i=2;
Names_i={'low','high'};

PTypeDistParamNames={'ptypeweights'};
Params.ptypeweights=[0.6; 0.4];

% A parameter that genuinely differs across PTypes (used in tests 1 and 3)
% Row ii is the kappa_j profile for PType ii.
Params.kappa_j_pt=[Params.kappa_j; 1.2*Params.kappa_j];
