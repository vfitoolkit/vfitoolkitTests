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
% The Howards accelerators (same answer as no Howards, and are they actually faster) are checked
% separately, in the diagnostic scan CoreInfHorzVFIAlgo_ScanHowardsSettings at the end.
%
% NOT tested here:
%   - divide-and-conquer: not usable for InfHorz (too slow); tested in the InfHorz-TPath bank
%   - semiz: not implemented for InfHorz value function iteration
% Some (option x GI) combinations are not yet implemented in the toolkit; those
% lines are included but commented out (see CoreInfHorzVFIAlgo_algocompare.m).
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Diary of the command window output (written to the parent bank's TestOutput folder)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreInfHorzVFIAlgoTestsdiary.txt','file')
    delete('../TestOutput/CoreInfHorzVFIAlgoTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreInfHorzVFIAlgoTestsdiary.txt

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

% NOTE: the diagnostic scans below are commented out. They were used to set the postGI and
% Howards defaults, those defaults are now in place, and the scans are slow. Two of them
% (ScanPostGIrepeat and ScanPostGIrepeat_nod) also need a one line change before they will
% run again: they sweep postGIrepeat>0 without setting vfoptions.howardssparse, so they now
% pick up the howardssparse=1 default and hit the sparse raws, which do not implement
% postGIrepeat>0. Adding vfo.howardssparse=0 where they build vfo fixes that.
% ScanHowardsSettingsWith2A is left running: it is the only coverage of the two endogenous
% state code paths in this bank.

%% DIAGNOSTIC scans --- Used to figure out defaults for postGI vfoptions settings
% Characterisation of the postGI window/basin behaviour -- these print NON-zero by design
% Based on these, I decided postGIrepeat is fairly useless
% And set a conservatively large maxaprimediff, will be a bit slower than
% things could be, but trying to set something large enough tn ensure that it always converges.
output=CoreInfHorzVFIAlgo_ScanMaxaprimediff(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
output=CoreInfHorzVFIAlgo_ScanPostGIrepeat(n_z,z_grid,pi_z,Params,DiscountFactorParamNames);
% Same again but for a nod model (no n_d loop, and maxaprimediff centred on 5, which is the nod default)
output=CoreInfHorzVFIAlgo_ScanPostGIrepeat_nod(n_z,z_grid,pi_z,Params,DiscountFactorParamNames);

%% Howards settings: do the Howards accelerators give the same answer as no Howards, and are they faster?
% z (no e) models, with and without d, at n_z=5,15,25,75 and n_a=100,200,500.
% Iterated Howards is tried at howards=40,80,120 so it can be compared against greedy Howards.
% (n_z is swept because the vfoptions.howardssparse default triggers on N_z>100, but the timings
%  had only ever been done at n_z=5. This scan builds its own z, so takes no z inputs.)
output=CoreInfHorzVFIAlgo_ScanHowardsSettings(Params,DiscountFactorParamNames);

%% Howards defaults on big grids: is the vfoptions.howardssparse default set on the right rule?
% Without the grid interpolation layer that default is 'N_a>1200 && N_z>100', and the scan above
% never fires it (it stops at n_a=500, n_z=75). Five (n_z,n_a) cells straddling the two conditions,
% with a cut-down config set (no howardsgreedy>0, no howards=120) since those are already settled.
% These grids are big enough to run out of GPU memory, so the solves are wrapped in try/catch.
output=CoreInfHorzVFIAlgo_ScanHowardsBigGrids(Params,DiscountFactorParamNames);

%% Howards settings for models with TWO endogenous states
% Part A: without GI, sweep the howardssparse/lowmemory configs (same raws as one endogenous state,
%         just a bigger N_a, but two endogenous states is how models realistically get a big N_a).
% Part B: with GI, sweep only the number of Howards iterations. Two endogenous states only implements
%         howardsgreedy=0 with howardssparse=0 there, so there is no config comparison to make.
% Part C: check the defaults dispatch, ie solve with nothing set and confirm it runs and matches.
output=CoreInfHorzVFIAlgo_ScanHowardsSettingsWith2A(Params,DiscountFactorParamNames);

diary off
