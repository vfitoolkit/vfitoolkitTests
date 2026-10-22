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
%             howardsgreedy=0,1,2,3 (2,3 are HowardMix), sparse=1/lowmemory=1, and postGIrepeat=1
%   Part 3: WITH GI, preGI and postGI give the same V and Policy
%   Part 4: with a big n_a (=1500), with and without GI give very similar V
% With a d, Part 1 also checks solnmethod='purediscretization' (the with-d raws) against the default refinement.
% Without semiz, Part 1 also checks vfoptions.endogenousvfi=1 (endogenous VFI) for each solver, plus a discriminator that a
% plain-VFI twin would fail (see CoreInfHorzVFIAlgo_algocompare); Parts 2, 3 and 3b check it with GI (postGI, preGI, two
% endogenous states); doPart(10) times it, with and without GI.
% Without e and semiz, Part 1 also checks outputkron=1 against the default.
% The models without z and without e stop after Part 1 (their grid interpolation layer is a single postGI solve).
% The e models run Part 1, and of Part 2 all but the HowardMix (howardsgreedy=2,3) checks (with e the grid interpolation
% layer is so far only postGI, howardsgreedy=0 or 1, one endogenous state); they skip Parts 3-5.
% The Howards accelerators (same answer as no Howards, and are they actually faster) are checked
% separately, in the diagnostic scan CoreInfHorzVFIAlgo_ScanHowardsSettings at the end.
%
% NOT tested here:
%   - divide-and-conquer: not usable for InfHorz (too slow); tested in the InfHorz-TPath bank
%   - semiz with the grid interpolation layer beyond postGI (doPart(8) and doPart(9) run a short GI part for the semiz
%     models: reference, howards=0, and lowmemory=1 with a d1 or sparse and greedy without)
% Some (option x GI) combinations are not yet implemented in the toolkit; those
% lines are included but commented out (see CoreInfHorzVFIAlgo_algocompare.m).
%
% doPart(11) is a DELIBERATE TRIGGER, not an algorithm comparison: a model whose V contains -Inf (an asset grid below the
% natural borrowing limit; and a -Inf cap on the state at the top z), to fire the grid interpolation layer's Howards step
% putting a zero interpolation weight on a -Inf node (0*(-Inf)=NaN, then priced at 0). See InfHorzGI_InfNaN_proposal.md at
% the toolkit root and CoreInfHorzVFIAlgo_InfNaNcompare. It runs the base postGI raws (nod/d, z/noz, two endogenous states)
% and the _e raws as a positive control (Figs 1-8, 11-16); the semiz raws and their endogenous twins (Figs 21-28, 31-34); and the sophisticated
% quasi-hyperbolic postGI raws at beta0=1 (Figs 41-44). These are the only figures in this bank.
%
% Only doPart(11) draws figures (saved to ../TestOutput as they are drawn); the other parts save only the diary.

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
% doPart(8): semiz (with/without d1, with/without z; Part 1 and a short postGI part)
% doPart(9): semiz with e (with/without d1, with/without z; Part 1 and a short postGI part)
% doPart(10): DIAGNOSTIC endogenous VFI runtimes (vfoptions.endogenousvfi=1 against plain VFI; slow)
% doPart(11): -Inf in V with the grid interpolation layer (deliberate trigger of the Howards zero-weight-on-a--Inf-node defect)
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1,1,1];

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
    output=CoreInfHorzVFIAlgo_nod_noz_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_nod_noz_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_nod_z_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_nod_z_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(1): without d

%% ===== doPart(2): with d =====
if doPart(2)==1
    fprintf('\n===== doPart(2): with d =====\n')
    %% with d
    output=CoreInfHorzVFIAlgo_d_noz_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_d_noz_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_d_z_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_d_z_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
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
    % And how the number of Howards steps should depend on beta (the defaults were chosen at beta=0.95)
    output=CoreInfHorzVFIAlgo_ScanHowardsBeta(Params,DiscountFactorParamNames);
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

%% ===== doPart(8): semiz =====
if doPart(8)==1
    fprintf('\n===== doPart(8): semiz =====\n')
    %% With a semi-exogenous state semiz (its transitions depend on the decision variable d2), with/without d1, with/without
    % z. Part 1, and a short GI part (postGI: reference, howards=0, and lowmemory=1 with a d1 or sparse and greedy without).
    output=CoreInfHorzVFIAlgo_nod1_noz_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_nod1_z_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_d1_noz_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_d1_z_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(8): semiz

if doPart(9)==1
    fprintf('\n===== doPart(9): semiz with e =====\n')
    %% As doPart(8), with an iid e as well.
    output=CoreInfHorzVFIAlgo_nod1_noz_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_nod1_z_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_d1_noz_e_semiz(n_d_semiz,n_a,n_a_big,n_z,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzVFIAlgo_d1_z_e_semiz(n_d_semiz,n_a,n_a_big,n_z,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(9): semiz with e

%% ===== doPart(10): DIAGNOSTIC endogenous VFI runtimes =====
if doPart(10)==1
    fprintf('\n===== doPart(10): DIAGNOSTIC endogenous VFI runtimes =====\n')
    %% Does vfoptions.endogenousvfi=1 earn its keep? Plain vs endogenous VFI runtimes across beta, rho_z and n_z
    % (the V/Policy checks on endogenous VFI are in Part 1 of every model, in CoreInfHorzVFIAlgo_algocompare)
    output=CoreInfHorzVFIAlgo_ScanEndogenousVFI(Params,DiscountFactorParamNames);
end % doPart(10): DIAGNOSTIC endogenous VFI runtimes

%% ===== doPart(11): -Inf in V with the grid interpolation layer (deliberate trigger) =====
if doPart(11)==1
    fprintf('\n===== doPart(11): -Inf in V with the grid interpolation layer (deliberate trigger) =====\n')
    %% GPU probes first (InfHorzGI_InfNaN_proposal.md, section 5 step 0). Two things the proposal could only verify on CPU
    % MATLAB decide whether the sparse Howards raws (class B) and the preGI interp1 Howards (class C) share the defect:
    %  - interp1 on a gpuArray with a -Inf node: CPU gives -Inf on the open intervals touching that node, exact values at the
    %    nodes, and no NaN. A NaN here would make the main-loop interp1 a second entry point for the defect.
    %  - sparse() on gpuArray inputs with an explicit zero weight, times a vector holding -Inf: CPU drops the zero at
    %    construction, so the product is 5 where a dense [1 0]*[5;-Inf] is NaN. If the GPU keeps the zero, the sparse
    %    Howards raws do 0*(-Inf) too and join the fix.
    probeinterp=gather(interp1(gpuArray((1:3)'),gpuArray([-5;-Inf;-3]),gpuArray(linspace(1,3,9)')));
    fprintf('GPU probe: interp1 with a -Inf middle node, 9 query points from node 1 to node 3: %s \n',mat2str(probeinterp',4))
    fprintf('GPU probe: interp1 with a -Inf node gives no NaN,                      should be one: %i \n',~any(isnan(probeinterp)))
    probeS=sparse(gpuArray([1;1]),gpuArray([1;2]),gpuArray([1;0]),1,2);
    probeSv=gather(full(probeS*gpuArray([5;-Inf])));
    fprintf('GPU probe: sparse() with an explicit zero weight: nnz %i (CPU drops it, 1); times [5;-Inf]: %s (CPU gives 5) \n',nnz(probeS),num2str(probeSv))
    fprintf('GPU probe (informational, a property of the GPU not of the toolkit; the sparse raws no longer rely on it): explicit zero weight on a -Inf node is dropped: %i \n',probeSv==5)
    fprintf('GPU probe: dense [1 0]*[5;-Inf] on the GPU: %s (CPU gives NaN); 0*(-Inf) on the GPU: %s (CPU gives NaN) \n',num2str(gather(gpuArray([1,0])*gpuArray([5;-Inf]))),num2str(gather(gpuArray(0)*gpuArray(-Inf))))
    probeS2=sparse(gpuArray([1;2]),gpuArray([1;2]),gpuArray([1;1]),2,3); % no explicit zero at all; column 3 is touched by no row
    fprintf('GPU probe (informational): sparse with NO explicit zero, times [5;6;-Inf] (CPU gives 5 6): %s \n',num2str(gather(full(probeS2*gpuArray([5;6;-Inf])))'))

    %% Device 1: the asset grid extends below the natural borrowing limit (fires the zero weight on a dead LOWER node)
    % nod/d, z/noz, e (positive control), two endogenous states. Figs 1-8.
    output=CoreInfHorzVFIAlgo_InfNaN_belowNBL(n_d,d_grid,n_z,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,1);
    for figure_c=1:8
        exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzVFIAlgoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    end
    %% Device 2: a -Inf cap on the state at the top z (fires the zero weight on a dead UPPER node)
    % nod/d with z, e (positive control), two endogenous states. Figs 11-16.
    output=CoreInfHorzVFIAlgo_InfNaN_capAtTopZ(n_d,d_grid,n_z,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,11);
    for figure_c=11:16
        exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzVFIAlgoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    end
    %% The semiz models (SemiExo postGI raws, their EndogenousVFI and RelativeVFI twins): device 1 on all eight variants (Figs 21-28),
    % device 2 on the four with z (Figs 31-34)
    output=CoreInfHorzVFIAlgo_InfNaN_belowNBL_semiz(n_d2_semiz,n_d_semiz,d2_grid_semiz,d_grid_semiz,n_z,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,21);
    for figure_c=21:28
        exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzVFIAlgoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    end
    output=CoreInfHorzVFIAlgo_InfNaN_capAtTopZ_semiz(n_d2_semiz,n_d_semiz,d2_grid_semiz,d_grid_semiz,n_z,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,31);
    for figure_c=31:34
        exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzVFIAlgoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    end
    %% Sophisticated quasi-hyperbolic at beta0=1 (the two QuasiHyperbolicS postGI raws), both devices, nod and d (Figs 41-44)
    output=CoreInfHorzVFIAlgo_InfNaN_QH(n_d,d_grid,n_z,z_grid,pi_z,Params,DiscountFactorParamNames,41);
    for figure_c=41:44
        exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzVFIAlgoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    end
end % doPart(11): -Inf in V with the grid interpolation layer (deliberate trigger)

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/CoreInfHorzVFIAlgoTestsdiary.txt')

diary off
