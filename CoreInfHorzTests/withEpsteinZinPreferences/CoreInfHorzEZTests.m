% Implement tests of the core VFI Toolkit InfHorz commands, under EPSTEIN-ZIN preferences.
% with/without d
% with/without low memory
% three preference cases (see below)
%
% This is the EZ mirror of CoreInfHorzTests.m, and the EZ sibling of CoreInfHorzQHTests.m. It
% reuses the parent bank's setup (CoreInfHorz_setup) and has its own ReturnFns (EZ_ReturnFns/),
% because the parent's (c^(1-sigma)-1)/(1-sigma) is mixed-sign and so unusable under Epstein-Zin.
%
% Per the exotic-preferences convention, this bank tests V and Policy ONLY. There is no agent
% distribution, no panel, no statistics: everything downstream of Policy is computed by code
% that never sees vfoptions.exoticpreferences, and the exponential bank already covers it.
%
% Where QH runs each model twice (Naive, Sophisticated), EZ runs each model three times, exactly
% as CoreFHorzEZTests does:
%   Case 1: consumption-units (traditional Epstein-Zin, vfoptions.EZutils=0)
%   Case 2: utility-units with positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%   Case 3: utility-units with negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% Each subcode draws Case 1 as figure_c, Case 2 as 100+figure_c, and Case 3 as 200+figure_c.
%
% WHAT THE INFHORZ EPSTEIN-ZIN SOLVER IS (2026-10-22)
% One dispatcher, ValueFnIter_InfHorz_EpsteinZin, and four raws: with/without d, lowmemory 0/1.
% That is the whole of it. In particular:
%   - vfoptions.solnmethod is NOT a second code path: the with-d raw takes a flat max over
%     (d,aprime), there is no Epstein-Zin Refine raw, so 'purediscretization_refinement' (the
%     default with d) and 'purediscretization' reach the same raw. A check between them would be
%     two names for one solve, so there is none.
%   - vfoptions.howardsgreedy and vfoptions.howardssparse are not read by any Epstein-Zin raw
%     (only howards and maxhowards are). Not checked, for the same reason.
%   - the grid interpolation layer is NOT implemented. Until 2026-10-22 vfoptions.gridinterplayer=1
%     was silently ignored (coarse solve, one-row Policy); it now errors, and the guards cross-test
%     asserts that. A GI implementation is the planned Tier 2 of this bank (see
%     ../../CoreInfHorzEZTests_proposal.md); when it lands, figs 3-4 (nod_z_GI, d_z_GI) go in as a sixth part
%     and guard 10 is retired.
%   - a feasible return of exactly zero is treated as infeasible by every Epstein-Zin raw, in both
%     horizons (temp2(ReturnMatrix==0)=-Inf). The three ReturnFn families of parts 1-5 are strictly
%     signed so those parts never meet it. Part 6 (added 2026-10-27) is the deliberate trigger.
%   - V=-Inf entries (a state with no feasible choice) are not handled consistently by the
%     negative-utility branch (the -Inf passes through the (ezc4*V)^ezc5 transform untouched and
%     comes out of ezc3*beta*(.) as +Inf). Same in the FHorz raws. The models of parts 1-5 have
%     every state feasible (a_grid starts at 0, income w*z>0), so V is finite everywhere there.
%     Part 6 is the deliberate trigger for this too.
%
% NOT tested here (unsupported, each refused with an interface message that the guards
% cross-test asserts):
%   - noz: every Epstein-Zin InfHorz raw loops over z, and with n_z=0 that loop never ran, so
%     the solve returned V=0 and Policy=0 silently. Now an error (2026-10-22).
%   - e / semiz / relativevfi: refused by ValueFnIter_InfHorz itself.
%   - survivalprobability / WarmGlowBequestsFn / EZmortalityriskaversion / EZoneminusbeta=2:
%     all need ages; EZmortalityriskaversion and EZoneminusbeta=2 used to be ignored, now errors.
%   - EZutils=2: mentioned in an old comment as 'traditional EZ with endogenous labor', never
%     implemented anywhere; it used to die on an undefined variable inside a raw, now an error.
%   - divide-and-conquer: not usable for InfHorz at all.
%
% TESTED AS WORKING, having never been tried (2026-10-22):
%   - two endogenous states (2A), on the plain path. The Epstein-Zin raws only ever use
%     prod(n_a), unlike the quasi-hyperbolic raws, which fail inside repelem on a vector n_a
%     and are guarded. Test 11 of the guards cross-test solves an a2-ignored model and checks it
%     against the one-endogenous-state solve. If that test fails on the GPU, the decision is to
%     turn it into a guard (as QH has), not to extend the raws.
%   - a two-dimensional exogenous state: the z_gridvals defect that the QH bank found in
%     87bca184 had the same shape in the Epstein-Zin branch and was fixed with it, with no
%     Epstein-Zin coverage. The z-disguise cross-test now covers it, on both lowmemory paths.
%
% CONFIRMED ONCE, NOT KEPT AS CHECKS (2026-10-27 run): vfoptions.maxiter and vfoptions.verbose were
% threaded into the four Epstein-Zin raws on 2026-10-22 (no raw read either before; the iteration
% had no cap at all). A one-off block in the guards cross-test confirmed that a user-set maxiter
% changes the answer (max|V(maxiter=3)-V(maxiter=10)|=8.03e-3; against the converged V, 8.03e-3 and
% 1.94e-6; both solves printed the not-converged warning), that verbose=1 prints the
% 'ValueFnIter (Epstein-Zin)' progress line, and that verbose does not change the solution
% (exactly 0). The block was then deleted, as the QH equivalents were on 2026-10-21: they pin
% option plumbing no model result depends on. Do not re-add them; re-run by hand if needed.
%
% Policy comparisons between two DIFFERENT solves of the same problem (the collapse cross-tests)
% are not identities in infinite horizon: two fixed-point iterations can legitimately pick
% different grid points where the objective is within tolerance. So every such comparison prints
% the raw index difference as a DIAGNOSTIC, and grades instead the Epstein-Zin VALUE of the other
% solve's policy against V (via ValueFnFromPolicy_InfHorz): if the policies differ only at ties,
% the two values agree to convergence tolerance, and if either policy is wrong they do not. The
% QH bank paid two GPU runs to learn that a raw index comparison measures a representation.
% Same-model disguises (split/permuted/degenerate z, lowmemory, 2A-ignored) are the identical
% computation up to floating-point order and are graded on the raw indices, as QH does.
%
% ValueFnFromPolicy_InfHorz
% Supports Epstein-Zin as of 2026-10-22 (ValueFnFromPolicy_InfHorz_EpsteinZin, dispatched ahead
% of the grid-interpolation hand-off exactly as the quasi-hyperbolic one is). Every subcode calls
% it at the solver's own Policy: V should be reproduced to convergence tolerance. It is also the
% discriminator described above.
%
% Convergence tolerance: every solve in this bank sets vfoptions.tolerance=1e-10 (default 1e-9).
% Two independent solves of one problem each sit within about tolerance/(1-beta) of the fixed
% point, so with beta=0.95 a 1e-9 tolerance leaves a worst-case gap of 4e-8 against CoreSummary's
% 1e-7 identity bar: too close. At 1e-10 the worst case is 4e-9.
%
% z IS OVERRIDDEN after CoreInfHorz_setup (see below): the parent's z is too tight for
% risk aversion to move the policy, and the bank's anti-vacuity lines need it to.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): solver variants (figs 1-2)
% doPart(2): brute-force oracle (no figure)
% doPart(3): cross-test: vNM collapses, EZoneminusbeta scaling, homogeneity
% doPart(4): cross-test: z disguises (split, permute, degenerate, two-dimensional)
% doPart(5): cross-test: interface guards and the 2A probe
% doPart(6): feasibility conventions: -Inf states and zero returns (deliberate trigger, no figures)
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1];

%% Diary of the command window output (figures are saved into the parent bank's TestOutput)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreInfHorzEZTestsdiary.txt','file')
    delete('../TestOutput/CoreInfHorzEZTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreInfHorzEZTestsdiary.txt
fprintf('CoreInfHorzEZTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('../CoreInfHorzTests_Setup/')
addpath('./EZ_ReturnFns/')

addpath('./withEpsteinZinPreferences_subcodes/')
addpath('./withEpsteinZinPreferences_subcodes/CrossTests/')

% Setup so that use the same d,a,z in all the models that use them
CoreInfHorz_setup

%% A riskier z than the parent bank uses
% CoreInfHorz_setup gives z from discretizeAR1_FarmerToda(0,0.9,0.03,5): log-z stationary std
% 0.069, z in [0.87,1.15]. With that little risk the Epstein-Zin policy may not move a single
% grid point away from the vNM one, and the bank's anti-vacuity lines (max|Policy-Policyref|
% should be well above zero) would be reporting the calibration rather than the solver. Override
% AFTER the setup call, so the shared CoreInfHorz_setup.m is untouched and the exponential bank
% keeps its own z. The innovation std goes from 0.03 to 0.2 (stationary std 0.46). Nothing else
% changes: n_a=101, n_d=9, the cubic a_grid, beta=0.95, r=0.05, w=1, varphi=0.8.
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.2,n_z);
z_grid=exp(z_grid);

%% Epstein-Zin preference parameters (the EZ vfoptions themselves are set inside each subcode)
% Same values as CoreFHorzEZTests.
Params.ezgamma=3; % consumption-units risk aversion (vfoptions.EZriskaversion when EZutils=0; must be >1)
Params.ezphi=0.5; % consumption-units elasticity of intertemporal substitution (vfoptions.EZeis; cannot be 1)
Params.ezrisk=3; % utility-units additional risk aversion (vfoptions.EZriskaversion when EZutils=1)
Params.ezsigma=2; % curvature of the positiveUtils/negativeUtils utility fns (>1; the positiveUtils
% family uses a (1+x) shift so it stays strictly positive despite ezsigma>1)

% vfoptions.exoticpreferences='EpsteinZin';
% Case 1: vfoptions.EZutils=0; vfoptions.EZriskaversion='ezgamma'; vfoptions.EZeis='ezphi';
% Case 2: vfoptions.EZutils=1; vfoptions.EZpositiveutility=1; vfoptions.EZriskaversion='ezrisk';
% Case 3: vfoptions.EZutils=1; vfoptions.EZpositiveutility=0; vfoptions.EZriskaversion='ezrisk';

% Note: CoreInfHorz_setup puts n_e/e_grid/pi_e and the semiz fields onto vfoptionsbaseline. InfHorz
% Epstein-Zin supports neither, so the subcodes start from vfoptions=struct() and never see
% vfoptionsbaseline (it is not even passed in).

%% The only functions worth testing are the value fn ones: after you have Policy, everything
%% downstream anyway ignores the Epstein-Zin preferences.

%% ===== doPart(1): solver variants (figs 1-2) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): solver variants (figs 1-2) =====\n')
    %% without d, with z
    figure_c=1;
    output=EZInfHorz_nod_z(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreInfHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreInfHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% with d, with z
    figure_c=2;
    output=EZInfHorz_d_z(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreInfHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreInfHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)
end % doPart(1): solver variants (figs 1-2)

%% ===== doPart(2): brute-force oracle =====
if doPart(2)==1
    fprintf('\n===== doPart(2): brute-force oracle =====\n')
    %% Brute-force oracle: an independent implementation of the Epstein-Zin recursion
    % Solves a tiny model (n_a=5, n_z=2, no d) twice: by explicit loops written directly from the
    % Epstein-Zin FORMULAS, iterated to a fixed point (gamma/phi/risk enter directly; no ezc
    % coefficients, no toolkit subcodes, nothing shared with the toolkit implementation), and
    % through ValueFnIter_InfHorz with the bank's standard EZ vfoptions. All three cases. This is
    % the only test class in the bank that can catch 'toolkit and reference both wrong in the same
    % way', because every other check compares two toolkit code paths, or the toolkit against an
    % identity evaluated at the toolkit's own answer.
    output=EZInfHorz_bruteforce_oracle(Params,DiscountFactorParamNames);
end % doPart(2): brute-force oracle

%% ===== doPart(3): cross-test: vNM collapses, EZoneminusbeta scaling, homogeneity =====
if doPart(3)==1
    fprintf('\n===== doPart(3): cross-test: vNM collapses, EZoneminusbeta scaling, homogeneity =====\n')
    %% Cross-test: the special Epstein-Zin identities (these play the role the beta0=1 test plays in the QH bank)
    % (i)   consumption-units: gamma=1/phi collapses EZ to standard vNM (V related by
    %       V_EZ=((1-gamma)*V_vNM)^(1/(1-gamma)); same policy up to ties)
    % (ii)  utility-units: EZriskaversion=0 collapses EZ to standard vNM (positive and negative cases)
    % (iv)  vfoptions.EZoneminusbeta=1 versus manually scaling the return fn
    % (xii) homogeneity of degree one in consumption units, at gamma=3 and at gamma=20 (overflow stress)
    % (i) and (ii) again with d.
    output=EZInfHorz_CrossTests_collapse(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
end % doPart(3): cross-test: vNM collapses, EZoneminusbeta scaling, homogeneity

%% ===== doPart(4): cross-test: z disguises =====
if doPart(4)==1
    fprintf('\n===== doPart(4): cross-test: z disguises =====\n')
    %% Cross-test: the same model disguised as a different z process, the solution must follow
    % split (each z state duplicated), permuted, degenerate (all z equal, vs a 1-point z), and a
    % two-dimensional z vs the equivalent single z (on both lowmemory paths: the lowmemory raws
    % convert z to gridvals themselves).
    output=EZInfHorz_CrossTests_zdisguise(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
end % doPart(4): cross-test: z disguises

%% ===== doPart(5): cross-test: interface guards and the 2A probe =====
if doPart(5)==1
    fprintf('\n===== doPart(5): cross-test: interface guards and the 2A probe =====\n')
    %% Cross-test: the interface guards, the paths that are known not to be supported, and the 2A probe
    % Everything here is a try/catch check, so a guard that stops firing shows up as a failed
    % check rather than as a crash that aborts the rest of the bank.
    output=EZInfHorz_CrossTests_guards(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
end % doPart(5): cross-test: interface guards and the 2A probe

%% ===== doPart(6): feasibility conventions: -Inf states and zero returns =====
if doPart(6)==1
    fprintf('\n===== doPart(6): feasibility conventions: -Inf states and zero returns =====\n')
    %% DELIBERATE TRIGGER (2026-10-27) for the two conventions the header above documents as 'not tested':
    % (A) a -Inf continuation passes through the negative-utility transform untouched and comes out of
    %     ezc3*beta*(.) as +Inf; (B) a return of exactly zero at a feasible choice is rewritten to -Inf.
    % Model A: a_grid extended below the natural borrowing limit (the vNM solve gives the reference dead set).
    % Model B: a subsistence floor yielding a zero return, floorval=0 vs 1e-10 for continuity.
    % Written to run RED on the raws as they stood; the FHorz twin is part 8 of CoreFHorzEZTests.
    output=EZInfHorz_CrossTests_feasibility(n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
end % doPart(6): feasibility conventions: -Inf states and zero returns

%% Done!
% Terminating banner. The bank's last output before this is the guards cross-test, which PRINTS
% CAUGHT MATLAB ERROR STRINGS as data -- those are its expected results, not crashes. Without this
% banner a completed run and a run that died in the guards cross-test look identical at the tail
% of the diary. Per the 2026-09-20 convention, a stable file size is not evidence that a
% Dropbox-synced diary is complete; this line is.
fprintf('\n===== CoreInfHorzEZTests: RUN COMPLETE ===== \n')
fprintf('If you do not see this line at the end of the diary, the run did not finish. \n')
fprintf('Note: the guards cross-test prints CAUGHT error messages as data. That is the interface \n')
fprintf('guard firing, which is what the check beside each asserts, not a failure. Every check \n')
fprintf('states its own expected value inline. \n')

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/CoreInfHorzEZTestsdiary.txt')

diary off
