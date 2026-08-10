% Tests of the core VFI Toolkit InfHorz VALUE FUNCTION ITERATION algorithms.
% with/without d, with/without z, with/without e
%
% Only ValueFnIter_InfHorz is being tested here (not StationaryDist/stats).
% InfHorz has several value-function-iteration algorithm options; they all
% solve the same Bellman equation so should give the same V and Policy.
%
% For each of the 8 models we check:
%   Part 1: WITHOUT grid interpolation (GI), all give the same V and Policy:
%             howardsgreedy=0,1 and howardssparse=0,1, plus howards=0 (pure VFI,
%             no Howards acceleration) and lowmemory=1
%   Part 2: WITH GI, all (implemented) give the same V and Policy:
%             howardsgreedy=0,1,2,3 (2,3 are HowardMix), and sparse=1/lowmemory=1
%   Part 3: WITH GI, preGI and postGI give the same V and Policy
%   Part 4: with a big n_a (=1500), with and without GI give very similar V
%
% NOT tested here:
%   - divide-and-conquer: not usable for InfHorz (too slow); tested in the InfHorz-TPath bank
%   - semiz: not implemented for InfHorz value function iteration
% Some (option x GI) combinations are not yet implemented in the toolkit; those
% lines are included but commented out (see CoreInfHorzVFIAlgo_algocompare.m).

addpath('./CoreInfHorzVFIAlgoTests_subcodes/')
addpath('./CoreInfHorzVFIAlgoTests_Setup/')
addpath('./CoreInfHorzVFIAlgo_ReturnFns/')
% Setup so that use the same d,a,z,e in all the models that use them
CoreInfHorzVFIAlgo_setup

%% without d
% output=CoreInfHorzVFIAlgo_nod_noz_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
% output=CoreInfHorzVFIAlgo_nod_noz_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreInfHorzVFIAlgo_nod_z_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
% output=CoreInfHorzVFIAlgo_nod_z_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);

%% with d
% output=CoreInfHorzVFIAlgo_d_noz_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
% output=CoreInfHorzVFIAlgo_d_noz_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreInfHorzVFIAlgo_d_z_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
% output=CoreInfHorzVFIAlgo_d_z_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);

%% Cross-test: a d variable that does nothing (all d grid points equal) vs no d (parameter=d value)
% Should give the same V (and aprime policy). Run with and without GI; a GI-only difference
% isolates the with-d (Refine) grid-interpolation code path.
output=CoreInfHorzVFIAlgo_CrossTest_dummyd(n_a,n_a_big,n_z,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);

%% DIAGNOSTIC scans --- Used to figure out defaults for postGI vfoptions settings
% Characterisation of the postGI window/basin behaviour -- these print NON-zero by design
% Based on these, I decided postGIrepeat is fairly useless
% And set a conservatively large maxaprimediff, will be a bit slower than
% things could be, but trying to set something large enough tn ensure that it always converges.
output=CoreInfHorzVFIAlgo_ScanMaxaprimediff(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
output=CoreInfHorzVFIAlgo_ScanPostGIrepeat(n_z,z_grid,pi_z,Params,DiscountFactorParamNames);
