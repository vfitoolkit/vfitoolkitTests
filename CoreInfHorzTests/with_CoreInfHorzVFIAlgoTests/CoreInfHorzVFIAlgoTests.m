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

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): without d
% doPart(2): with d
% doPart(3): cross-test: a dummy d vs no d
% doPart(4): DIAGNOSTIC scans (slow; set to 0 to skip)
% doPart(5): Howards settings scan
% doPart(6): Howards defaults on big grids
% doPart(7): Howards settings with two endogenous states
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1];

%% Diary of the command window output (written to the parent bank's TestOutput folder)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreInfHorzVFIAlgoTestsdiary.txt','file')
    delete('../TestOutput/CoreInfHorzVFIAlgoTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreInfHorzVFIAlgoTestsdiary.txt
fprintf('CoreInfHorzVFIAlgoTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreInfHorzVFIAlgoTests_subcodes/')
addpath('./CoreInfHorzVFIAlgoTests_Setup/')
addpath('./CoreInfHorzVFIAlgo_ReturnFns/')
% Setup so that use the same d,a,z,e in all the models that use them
CoreInfHorzVFIAlgo_setup

%% ===== doPart(1): without d =====
if doPart(1)==1
    fprintf('\n===== doPart(1): without d =====\n')
    %% without d
    % output=CoreInfHorzVFIAlgo_nod_noz_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    % output=CoreInfHorzVFIAlgo_nod_noz_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_nod_z_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    % output=CoreInfHorzVFIAlgo_nod_z_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(1): without d

%% ===== doPart(2): with d =====
if doPart(2)==1
    fprintf('\n===== doPart(2): with d =====\n')
    %% with d
    % output=CoreInfHorzVFIAlgo_d_noz_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    % output=CoreInfHorzVFIAlgo_d_noz_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_d_z_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    % output=CoreInfHorzVFIAlgo_d_z_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(2): with d

%% ===== doPart(3): cross-test: a dummy d vs no d =====
if doPart(3)==1
    fprintf('\n===== doPart(3): cross-test: a dummy d vs no d =====\n')
    %% Cross-test: a d variable that does nothing (all d grid points equal) vs no d (parameter=d value)
    % Should give the same V (and aprime policy). Run with and without GI; a GI-only difference
    % isolates the with-d (Refine) grid-interpolation code path.
    output=CoreInfHorzVFIAlgo_CrossTest_dummyd(n_a,n_a_big,n_z,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(3): cross-test: a dummy d vs no d

% NOTE: the diagnostic scans below are LIVE. This note used to say they were commented out; they
% are not, and the run of the bank in ./TestOutput shows SCAN 1 and SCAN 2 producing their tables.
% They were used to set the postGI and Howards defaults, those defaults are now in place, and the
% scans are slow - doPart(4) is how you skip them now.
% Two of them (ScanPostGIrepeat and ScanPostGIrepeat_nod) still need a one line change: they sweep
% postGIrepeat>0 without setting vfoptions.howardssparse, so they pick up the howardssparse=1
% default and hit the sparse raws, which do not implement postGIrepeat>0. That is what the 'n/a'
% cells at maxaprimediff=30,50 in the SCAN 2 tables are. Adding vfo.howardssparse=0 where they
% build vfo fixes it.
% ScanHowardsSettingsWith2A is left running: it is the only coverage of the two endogenous
% state code paths in this bank.

%% ===== doPart(4): DIAGNOSTIC scans (slow; set to 0 to skip) =====
if doPart(4)==1
    fprintf('\n===== doPart(4): DIAGNOSTIC scans (slow; set to 0 to skip) =====\n')
    %% DIAGNOSTIC scans --- Used to figure out defaults for postGI vfoptions settings
    % Characterisation of the postGI window/basin behaviour -- these print NON-zero by design
    % Based on these, I decided postGIrepeat is fairly useless
    % And set a conservatively large maxaprimediff, will be a bit slower than
    % things could be, but trying to set something large enough tn ensure that it always converges.
    output=CoreInfHorzVFIAlgo_ScanMaxaprimediff(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
    output=CoreInfHorzVFIAlgo_ScanPostGIrepeat(n_z,z_grid,pi_z,Params,DiscountFactorParamNames);
    % Same again but for a nod model (no n_d loop, and maxaprimediff centred on 5, which is the nod default)
    output=CoreInfHorzVFIAlgo_ScanPostGIrepeat_nod(n_z,z_grid,pi_z,Params,DiscountFactorParamNames);
end % doPart(4): DIAGNOSTIC scans (slow; set to 0 to skip)

%% ===== doPart(5): Howards settings scan =====
if doPart(5)==1
    fprintf('\n===== doPart(5): Howards settings scan =====\n')
    %% Howards settings: do the Howards accelerators give the same answer as no Howards, and are they faster?
    % z (no e) models, with and without d, at n_z=5,15,25,75 and n_a=100,200,500.
    % Iterated Howards is tried at howards=40,80,120 so it can be compared against greedy Howards.
    % (n_z is swept because the vfoptions.howardssparse default triggers on N_z>100, but the timings
    %  had only ever been done at n_z=5. This scan builds its own z, so takes no z inputs.)
    output=CoreInfHorzVFIAlgo_ScanHowardsSettings(Params,DiscountFactorParamNames);
end % doPart(5): Howards settings scan

%% ===== doPart(6): Howards defaults on big grids =====
if doPart(6)==1
    fprintf('\n===== doPart(6): Howards defaults on big grids =====\n')
    %% Howards defaults on big grids: is the vfoptions.howardssparse default set on the right rule?
    % Without the grid interpolation layer that default is 'N_a>1200 && N_z>100', and the scan above
    % never fires it (it stops at n_a=500, n_z=75). Five (n_z,n_a) cells straddling the two conditions,
    % with a cut-down config set (no howardsgreedy>0, no howards=120) since those are already settled.
    % These grids are big enough to run out of GPU memory, so the solves are wrapped in try/catch.
    output=CoreInfHorzVFIAlgo_ScanHowardsBigGrids(Params,DiscountFactorParamNames);
end % doPart(6): Howards defaults on big grids

%% ===== doPart(7): Howards settings with two endogenous states =====
if doPart(7)==1
    fprintf('\n===== doPart(7): Howards settings with two endogenous states =====\n')
    %% Howards settings for models with TWO endogenous states
    % Part A: without GI, sweep the howardssparse/lowmemory configs (same raws as one endogenous state,
    %         just a bigger N_a, but two endogenous states is how models realistically get a big N_a).
    % Part B: with GI, sweep only the number of Howards iterations. Two endogenous states only implements
    %         howardsgreedy=0 with howardssparse=0 there, so there is no config comparison to make.
    % Part C: check the defaults dispatch, ie solve with nothing set and confirm it runs and matches.
    output=CoreInfHorzVFIAlgo_ScanHowardsSettingsWith2A(Params,DiscountFactorParamNames);
end % doPart(7): Howards settings with two endogenous states

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/CoreInfHorzVFIAlgoTestsdiary.txt')

diary off
