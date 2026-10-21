% Implement tests of the core VFI Toolkit InfHorz commands, under QUASI-HYPERBOLIC discounting.
% with/without d
% with/without grid interpolation
% Naive and Sophisticated
%
% This is the QH mirror of CoreInfHorzTests.m. It reuses that bank's setup and ReturnFns.
%
% Per the exotic-preferences convention, this bank tests V and Policy ONLY. There is no agent
% distribution, no panel, no statistics: everything downstream of Policy is computed by code
% that never sees vfoptions.exoticpreferences, and the exponential bank already covers it.
%
% INTERFACE
% vfoptions.QHadditionaldiscount is CHAR ONLY ('beta0'). The cell form {'beta0'} errors, and so
% does a non-scalar beta0 (age-dependent beta0 was deliberately removed). Both are checked by the
% guards cross-test, which is the only place in any bank that still writes the cell form, and
% writes it deliberately so that the refusal is exercised.
% This is no longer an InfHorz-only interface: the char-only guard is now in five dispatchers
% across ValueFnIter/ValueFnFromPolicy/TransitionPaths, every bank declares the char form, and no
% live code brace-indexes the field anywhere in the toolkit. An earlier version of this comment
% warned that the FHorz QH banks still used the cell form and not to copy them; that migration has
% since happened and the warning no longer applies.
%
% NOT tested here:
%   - noz: InfHorz QH has no noz code path at all (every raw indexes pi_z and loops over z).
%   - e / semiz: no InfHorz QH support.
%   - two endogenous states (2A): not supported on any InfHorz QH path, and since 2026-10-21 the
%     interface guard sits ahead of every path rather than only inside the GI block, so both entry
%     points refuse it with the same message. Guard 5 and test 7 of the guards cross-test assert
%     that, including that the message is the interface one and not a failure inside a raw.
% CONFIRMED ONCE, NOT KEPT AS CHECKS (2026-10-21 run): the dispatcher used to force both
% vfoptions.maxiter=1000 and vfoptions.verbose=1 for Sophisticated, and the verbose forcing printed
% nothing because no quasi-hyperbolic raw read the flag. After making maxiter a default and
% threading verbose into the four Sophisticated raws, a one-off pair of checks confirmed that a
% user-set maxiter changes the answer (max|V(maxiter=3)-V(maxiter=10)|=1.708, where the override
% would have forced both to 1000 and made them identical), that verbose=1 prints the progress line
% and verbose=0 does not, and that verbose does not change the solution (exactly 0). Those two
% checks were then removed: they pin option plumbing no model result depends on, and one of them
% tested evalc as much as the toolkit. Do not re-add them; re-run them by hand if the forcing ever
% comes back.
%
% Tested as of 2026-10-21, having previously been an unsupported shape:
%   - a two-dimensional exogenous state. The exotic-preference branches were being handed z_grid
%     where every other solver branch gets z_gridvals, so multi-dim z errored; test 6 of the guards
%     cross-test now solves it and checks it against the equivalent single z.
%   - divide-and-conquer: not usable for InfHorz at all.
%
% KNOWN FAILING CHECKS: none. The bank went fully green on 2026-10-21 at 213 checks, 0 failed.
% The history is worth keeping, because each round of failures was a different KIND of mistake:
% 16 failures on 2026-09-29, of which a Howards accelerator bug in the Sophisticated raws was real
% and fixed, while the rest were Sophisticated convergence DIAGNOSTICS that carried check wording
% and so were being graded as checks (they now say DIAGNOSTIC and are not graded). Then 3 failures
% from the Naive postGI maxaprimediff window being short of the width actually needed -- resolved
% not by widening the default but by removing the window from the quasi-hyperbolic step entirely
% (it now goes through ValueFnIter_Case1_FHorz at N_j=1). Then 2 failures that were not
% disagreement at all but the FHorz and InfHorz families encoding a coarse-grid-point optimum
% differently, fixed by comparing decoded policy values instead of raw indices.
% Two of those four rounds were the bank measuring a representation rather than a quantity.
% Nothing here asserts tight Sophisticated convergence: the iteration cycles (about half the state
% space changes policy every iteration, measured 2026-09-07) and the postGI Sophisticated answer is
% separately limited by the kinked-Vunderbar collapse at beta0<1, which no window width fixes.
%
% ValueFnFromPolicy_InfHorz
% Every subcode calls ValueFnFromPolicy_InfHorz for both Naive and Sophisticated: eight checks,
% a second reconstruction oracle independent of the hand-written Bellman identities in figs 1-4.
% These were written test-first and the toolkit support landed 2026-08-30 (varargout signature,
% quasi-hyperbolic dispatch ahead of the grid-interpolation hand-off, and the new
% ValueFnFromPolicy_InfHorz_QuasiHyperbolic), so they are ordinary checks now, no try/catch.
%   Naive:         [V,Valt] = [Vtilde, exponential-discounter value at Policyalt];
%                  requires vfoptions.Policyalt (the 4th ValueFnIter output for Naive).
%                  Expect near zero, ~1e-8.
%   Sophisticated: [V,Valt] = [Vhat, Vunderbar]; Policy alone suffices.
%                  Expect NOT near zero: ValueFnFromPolicy solves the Vunderbar fixed point
%                  properly while the solver stops at the forced maxiter=1000, so the gap
%                  measures that non-convergence rather than any fault. Labelled
%                  CONVERGENCE-DEPENDENT in the diary.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): solver variants (figs 1-4)
% doPart(2): cross-test: beta0=1 vs exponential
% doPart(3): cross-test: z relabel/duplicate
% doPart(4): cross-test: Howards on vs off, Sophisticated
% doPart(5): cross-test: interface guards
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1];

%% Diary of the command window output (figures are saved into the parent bank's TestOutput)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreInfHorzQHTestsdiary.txt','file')
    delete('../TestOutput/CoreInfHorzQHTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreInfHorzQHTestsdiary.txt
fprintf('CoreInfHorzQHTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('../CoreInfHorzTests_Setup/')
addpath('../CoreInfHorz_ReturnFns/')

addpath('./withQuasiHyperbolicDiscounting_subcodes/')
addpath('./withQuasiHyperbolicDiscounting_subcodes/CrossTests/')

% Setup so that use the same d,a,z in all the models that use them
CoreInfHorz_setup

%% Finer grids than the parent bank uses
% CoreInfHorz_setup gives n_a=101, n_d=9. That is too coarse here: the quasi-hyperbolic policy
% sits a long way from the exponential one, so on a coarse grid the postGI search window is
% measured in a handful of very unevenly spaced points and the Sophisticated fixed point chatters
% between adjacent grid points instead of converging. Override AFTER the setup call, so the
% shared CoreInfHorz_setup.m is untouched and the exponential bank keeps its own grids.
n_a=251;
n_d=51;
a_grid=5*linspace(0,1,n_a)'.^3;  % same cubic spacing as the parent bank, more points
d_grid=linspace(0,1,n_d)';
% n_a_big stays at 1001: (1001-1)/(251-1)=4, so a_grid_big(1:4:end)==a_grid exactly, which is what
% the grid-interpolation accuracy comparison in fig 3 relies on.
% n_z stays at 5, and z is deliberately tight (log-z stationary std 0.0688, so z in [0.87,1.15]).

Params.beta0=0.9; % additional today-tomorrow discount factor
vfoptionsbaseline.QHadditionaldiscount='beta0'; % CHAR, not cell: the cell form now errors

% Note: CoreInfHorz_setup puts n_e/e_grid/pi_e onto vfoptionsbaseline. InfHorz QH has no e
% support, so the subcodes deliberately start from vfoptions=struct() and copy across only
% QHadditionaldiscount. Do not pass vfoptionsbaseline itself into ValueFnIter_InfHorz here.

% vfoptions.exoticpreferences='QuasiHyperbolic';
% vfoptions.quasi_hyperbolic='Naive';
% vfoptions.quasi_hyperbolic='Sophisticated';

%% The only functions worth testing are the value fn ones: after you have Policy, everything
%% downstream anyway ignores the quasi-hyperbolic discounting.

%% ===== doPart(1): solver variants (figs 1-4) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): solver variants (figs 1-4) =====\n')
    %% without d, with z (pure discretization)
    figure_c=1;
    output=QHDInfHorz_nod_z(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreInfHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% with d, with z (refinement, which is the default, AND plain pure discretization)
    figure_c=2;
    output=QHDInfHorz_d_z(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreInfHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% without d, with z, with the grid interpolation layer
    figure_c=3;
    output=QHDInfHorz_nod_z_GI(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreInfHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% with d, with z, with the grid interpolation layer
    figure_c=4;
    output=QHDInfHorz_d_z_GI(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreInfHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
end % doPart(1): solver variants (figs 1-4)

%% ===== doPart(2): cross-test: beta0=1 vs exponential =====
if doPart(2)==1
    fprintf('\n===== doPart(2): cross-test: beta0=1 vs exponential =====\n')
    %% Cross-test: beta0=1 must reproduce the standard exponential-discounting solution
    % This is the single strongest check in the bank. It runs every one of the eight raws against
    % the exponential solve of the identical model, and it is what catches a wrong discount factor
    % (bug 1) directly and loudly.
    output=QHDInfHorz_CrossTests_beta0one(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(2): cross-test: beta0=1 vs exponential

%% ===== doPart(3): cross-test: z relabel/duplicate =====
if doPart(3)==1
    fprintf('\n===== doPart(3): cross-test: z relabel/duplicate =====\n')
    %% Cross-test: relabel/duplicate the z states, the solution must follow
    output=QHDInfHorz_CrossTests_zdisguise(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(3): cross-test: z relabel/duplicate

%% ===== doPart(4): cross-test: Howards on vs off, Sophisticated =====
if doPart(4)==1
    fprintf('\n===== doPart(4): cross-test: Howards on vs off, Sophisticated =====\n')
    %% Cross-test: Howards on vs off, Sophisticated
    % Regression test for the Howards fix of 2026-09-01. Must use beta0<1: at beta0=1 the bug it targets is invisible, so the beta0=1
    % cross-test above cannot see it. Not a plain on==off check -- see the file header for why.
    output=QHDInfHorz_CrossTests_howards(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(4): cross-test: Howards on vs off, Sophisticated

%% ===== doPart(5): cross-test: interface guards =====
if doPart(5)==1
    fprintf('\n===== doPart(5): cross-test: interface guards =====\n')
    %% Cross-test: the interface guards, and the paths that are known not to be supported
    % Everything here is a try/catch check, so a guard that stops firing shows up as a failed
    % check rather than as a crash that aborts the rest of the bank.
    output=QHDInfHorz_CrossTests_guards(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(5): cross-test: interface guards

%% Done!
% Terminating banner. The bank's last output before this is test 7 of the guards cross-test,
% which PRINTS A CAUGHT MATLAB ERROR STRING as data -- that is its expected result, not a crash.
% Without this banner a completed run and a run that died in the guards cross-test look identical
% at the tail of the diary. Per the 2026-09-20 convention, a stable file size is not evidence that
% a Dropbox-synced diary is complete; this line is.
fprintf('\n===== CoreInfHorzQHTests: RUN COMPLETE ===== \n')
fprintf('If you do not see this line at the end of the diary, the run did not finish. \n')
fprintf('Note: test 7 of the guards cross-test prints a CAUGHT error message as data. That is the \n')
fprintf('interface guard firing, which is what the check beside it asserts, not a failure. Every \n')
fprintf('check states its own expected value inline. \n')

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/CoreInfHorzQHTestsdiary.txt')

diary off
