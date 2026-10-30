% Setup for CoreInfHorzTPathPTypeTests.
% Builds on CoreInfHorzTPath_setup (same d, a, z, e grids and Params, plus T, PricePath, ParamPath and
% transpathoptionsbaseline), then adds PType-specific pieces.

run('../CoreInfHorzTPathTests/CoreInfHorzTPathTests_Setup/CoreInfHorzTPath_setup.m')

%% PType-specific
N_i=2;
Names_i={'low','high'};

PTypeDistParamNames={'ptypeweights'};
Params.ptypeweights=[0.6; 0.4];
