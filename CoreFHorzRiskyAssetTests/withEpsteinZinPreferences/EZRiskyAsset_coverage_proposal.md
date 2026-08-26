# EZ-RiskyAsset test bank: coverage, decisions, and test-first state

Approved 2026-08-12. Mirrors `CoreFHorzTests/withEpsteinZinPreferences/` grafted onto
`CoreFHorzRiskyAssetTests/`. Scope: noa1 (figs 1-16) + withA1 (figs 17-32). 2a1 (two standard
endogenous assets alongside the risky asset) is OUT OF SCOPE for now — add later.

## V_Jplus1 legs — BUILT 2026-08-25 (unrun), BOTH riskyasset banks

64 subcodes (vNM bank: Noa1/Noa1-Semiz/WithA1/WithA1-Semiz x8; EZ bank: same four folders x8;
With2A1 excluded — other session's lane). Convention = the GPU-validated main-EZ-bank (v)
template: full solve, shorter model of Njs=jstar-1 periods with vfoptions.V_Jplus1=Vbase(...,jstar),
agej/kappa_j trimmed (mewj alone), exact V/Policy reproduction. noa1 shapes (no DC/GI): basic legs
at jstar=round(3*N_j/4) AND jstar=N_j (terminal-branch coverage), each with the file's own
lowmemory rungs. withA1 shapes: four methods (basic/DC/GI/DC+GI) at jstars 15/13/10/20, lowmemory
rungs per method. vNM sections are flat; EZ sections loop ezcase=1:3. All lint-clean vs baseline
(only the template's own trailing lowmemory=0 mlint note). FIRST-EVER execution of the V_Jplus1
branches in both the vNM riskyasset raws (incl. the batch items 5-7 ladder/reshape fixes) and the
EZ riskyasset raws — expect the branch-drift class here on run 1. After green: sj/warm-glow round
for the EZ riskyasset bank (approved next step).

Run 1 (vNM bank, 2026-08-26): several figs' legs green (both jstars, batch-fixed ladders
confirmed), crashed d1_z_e_noa1 jstar=15 lm1: Raw/noa1_e_raw:169 used ReturnMatrix where the
V_Jplus1 lm1 e-loop defines ReturnMatrix_e — FIXED. Full 132-file normalize-and-diff sweep of
every vNM riskyasset V_Jplus1 block vs its own main loop found 3 more latent V_Jplus1-lm0 bugs
(missing aprimeIndex(:) flatten + aprimeProbs repmat/reshape in withA1 blocks — the same
u-lottery z-offset class fixed in the EZ withA1 raws during the EZ build; error when N_u~=N_z,
SILENT when N_u==N_z): Raw/RiskyAsset_raw.m:101-111, RiskyAssetSemiExo/noz_raw.m:114-119,
RiskyAssetSemiExo/nod1_noz_raw.m:104-109 — all FIXED to mirror each file's main loop. Remaining
128 files verified line-identical to main loops after normalization. Awaiting rerun.

## sj/warm-glow round — BUILT 2026-08-26 (unrun; user-approved design)

Sections (vi)-(ix) added to all 32 subcodes + a 3-variant degenerateu bridge; driver gained the
main-bank sj/WG calibration (sj=[linspace(1,0.6,N_j-1),0], oneminussj, wg1=2, wg2=1, wg3=ezsigma,
ezmrisk=5); 3 new EZRiskyWarmGlowFn_{cons,positiveUtils,negativeUtils}.m (De Nardi forms, a2prime
argument). Design deltas vs the main bank: (vii).3 wrapper collapse oracles DROPPED everywhere
(vNM riskyasset has no WG/sj-solver support, and WG(aprime) cannot enter a riskyasset ReturnFn —
a2prime is realized via the (d,u) lottery after the return fn); exactness instead anchored by
(a) the degenerateu bridge (CrossTests: +sj / +warm-glow / +sj+warm-glow vs the GPU-green
main-family EZ solvers; cons-units INCLUDED — same terminal convention both sides) and (b) the
(vii).5 N_j-1 identity, INCLUDED in all 16 noa1 subcodes (closed forms use the riskyasset budget:
x=pension+a, or (pension+a)^varphi for d1; semiz/z/e verified absent from the retirement branch of
all 24 semiz ReturnFns), EXCLUDED in withA1 (terminal V depends on a1). (vi) collapses use vNM
riskyasset references with DiscountFactorParamNames={'beta','sj'}. withA1 exactness content =
basic/DC/GI/DC+GI four-method agreement with WG active. (vii).4 V_Jplus1+sj+WG mini-legs trim
agej/kappa_j/sj/oneminussj. ~2100 new checks. All lint-clean vs baseline (only the house-pattern
lowmemory=0 mlint note). FIRST-EVER execution of the EZ riskyasset raws'/FromPolicy's sj+WG code
paths (incl. WG-through-the-(d,u)-lottery) — expect debug iterations on run 1.

Run 1 (2026-08-26): fig 1 fully green (229 checks); crashed fig 2 (vii).2 lm1. Fix round: crash =
nod1_noa1_raw:140 terminal-lm1 WG refine row-vs-column (WGmatrix_onlyd3 row [1,N_d3] added
unshifted to [N_d3,N_a2] ReturnMatrix_z) — FIXED, + full 10-agent audit of every WG block in all
64 EZ riskyasset raws (base+3 tiers, nosemiz+SemiExo) vs their GPU-green lm0 blocks. 6 more bugs
FIXED: nod1_noa1_e_raw:182 (undefined ReturnMatrix_onlyd3 + same composition class, terminal lm2),
nod1_raw:165 (same class, terminal lm1; the with-d1 lookalikes are valid — both operands refined),
and 4 SemiExo e-variant base raws (noa1_e:101, noa1_noz_e:94, nod1_noa1_e:95, nod1_noa1_noz_e:89)
whose terminal-V_Jplus1 WG repmat guard 'lm0||lm1' was copied from non-e siblings but their
lm>=1 V_Jplus1 branch consumes full-(both)z temp4 at every rung → OOB mask crash at lm>=2; made
unconditional. All 6 tier folders + dispatchers audited CLEAN (WGcol column/pre-collapse-add
patterns structurally immune). Non-bug observations: noa1_noz/noa1 (and SemiExo noa1_noz_semiz)
omit the entireRHS==0→-Inf guard in VJplus1/interior blocks uniformly incl. green lm0 (cross-file
inconsistency, not drift). Awaiting rerun.

Run 2 (2026-08-26): END-TO-END, 6850 checks, 0 errors; degenerateu bridge FULLY GREEN (incl.
cons-units WG legs — the WG-through-lottery machinery validated); only 16 nonzero: basic-vs-DC,
cons-units only, the 4 with-d1 nosemiz withA1 shapes, identical V-diff 1.16573502. Fix round 2
INVERTED the suspect: the DC tier was CORRECT — the BASIC raws' terminal warmglow==1 block
misapplied the ezc9 min-trick to the POST-ezc7 (final-units) return object, so with ezc9=-1
(cons-units) basic picked the welfare-MINIMIZING d1 at terminal (h=0.5 not h=0; V diff =
c^0.8*(1-0.5^0.2), shape-independent in retirement). The WG(d2) refine is pre-ezc7-units and
correct everywhere — untouched. 49 sites/20 files fixed (bitwise no-op for ezc9=+1): 8 nosemiz
basic + 8 semiz basic + 4 semiz DC1 raws — the semiz basic AND semiz DC shared the defect
(cancelled in basic==DC checks; the bank compares basic==DC and GI==DC+GI, never basic==GI, which
is why semiz stayed 'green'); both sides moved to the joint-max convention in lockstep, matching
the semiz GI/DC_GI tiers. Nosemiz DC1 + all GI/DC_GI verified correct, unedited. V_Jplus1-branch
WG refines act on pre-ezc7 objects — correct on both sides, no change. Residual notes: (1) the
terminal convention adds the pre-ezc7-units WG term to the post-ezc7 return without an outer wrap
(inherited from the plain-EZ family; uniform across methods — convention quirk, not a bug);
(2) test-bank blind spot: with WG the bank never directly compares basic to GI (only basic==DC,
GI==DC+DC+GI pairs + the nod1-only degenerateu bridge) — a shared basic/DC defect can cancel, as
it did for semiz. Awaiting run 3.

## Terminal warm-glow convention fix — BUILT 2026-08-26 (unrun; user-approved)

Per Kraft/Munk/Weiss 2022 (JBF 138, 106428) + the Intro's EZ appendix: the warm-glow enters
INSIDE the ^ezc7 root at the terminal age (the bequest is the terminal condition, composed like
any continuation value). The old no-V_Jplus1 terminal blocks added it AFTER the root (wrong for
cons-units; coincided for utility-units; the V_Jplus1-branch terminal was already correct).
TOOLKIT SWEEP (~151 files, 10+2 agents, all lint-identical-to-baseline, git hunks confined to
preamble+terminal, V_Jplus1/interior untouched): main EZ 8 base + 24 tiers + 32 SemiExo;
riskyasset EZ 16 base + 24 tiers + 16 SemiExo base + 24 SemiExo tiers; 7 FromPolicy files (14
terminal blocks). Recipe: preamble WG transform → temp4 form (drop ezc3*beta); terminal rung =
the file's own interior/V_Jplus1 composition with EV→WG (pre-root temp2; riskyasset refines
d1/d2 with the ezc9 trick pre-root, ezc3 riding the WG refine — supersedes the run-2 'plain max
on post-root' terminal edits); warmglow==0 bitwise-preserved (branch split, or shared path with
DiscountedWG=zeros where the family's interior does that). Utility-units bitwise unchanged.
TESTS: main bank (vii).3a cons-units collapse legs RESTORED + (vii).5 cons-units legs added
(closed form (1+r)aprime+pension[, ^varphi]) = 64 new checks; riskyasset bank (vii).5 cons-units
legs (pension+aprime[, ^varphi]) = 32 new checks. Expect on GPU: previously-green utility-units
unchanged; new cons-units oracles zero; degenerateu bridge stays green (both sides moved).

## Decisions (user-approved)

1. **u-shock CE convention: JOINT.** The u-expectation is taken inside ONE joint
   certainty-equivalent over (u, semiz', z', e'). The gamma=1/phi and EZriskaversion=0 collapse
   tests only hold exactly under this convention, so it is enforced by test, not just documented.
2. Three EZ cases per subcode (cons-units / positiveUtils / negativeUtils), same parameters as the
   main EZ bank (ezgamma=3, ezphi=0.5, ezrisk=3, ezsigma=2).
3. Build order: full test bank first (test-first), then a toolkit-changes proposal (the Changes
   A-F analogue) once the tests define the target.

## Structure

- Driver `CoreFHorzRiskyAssetEZTests.m`: figs 1-32 (16 shapes noa1, then the same 16 withA1),
  three figures per fig-number (figure_c / 100+figure_c / 200+figure_c), then cross tests.
- Subcodes `EZRiskyAsset_<shape>[_withA1].m` in `withEpsteinZinPreferences_subcodes/`
  (Semiz_subcodes/, WithA1_subcodes/, WithA1_subcodes/Semiz_subcodes/ mirroring the vNM bank).
- ReturnFns `EZRiskyReturnFn_{cons,positiveUtils,negativeUtils}_<shape>[_withA1].m` in
  `EZRisky_ReturnFns/` (mirroring the vNM ReturnFn tree; 3 case-variants of each vNM ReturnFn).
- Per case: basic solve + ValueFnFromPolicy check + lowmemory legs (+ for withA1: DC, GI, DC+GI
  ladders with their lowmemory legs) + big-a_grid dist/AllStats/LifeCycleProfiles + figure.
- Special tests per subcode: (i) gamma=1/phi collapse (V_EZ=((1-ezgamma)V_vNM)^(1/(1-ezgamma)),
  Policy exact); (ii) EZriskaversion=0 collapse (both utility signs); (iii) the collapses under GI
  [withA1 subcodes ONLY — noa1 riskyasset has no DC/GI since there is no a1 to refine];
  (iv) EZoneminusbeta=1 vs manual scaling (scale factors as Params.ezscalefactor/ezscalefactoru —
  GPU arrayfun rejects workspace-captured anonymous functions).
- NOTE: the noz_noe shapes are NOT no-shock models here (u always provides risk), so no
  'EZ without shocks' warning is expected anywhere in this bank.

## Cross tests

Included (all EZ on both sides unless stated): zase (8 variants — the riskyasset analogue of the
main bank's CrossTests3 nested-vs-joint CE), semizasz (4), plainvswithA1 (4), d2recon regression
guard (1), and NEW degenerateu (1): riskyasset with riskyshare grid {0}, Params.r=0
(stored/restored) and d3 savings grid = a_grid (so aprime=savings exactly, on-grid) == the
plain (non-riskyasset) EZ savings model from the CoreFHorzTests EZ solvers — ties the two banks
together through a completely different code path. u is kept at baseline (n_u=3): it is genuinely
irrelevant under riskyshare=0, which additionally checks that a degenerate-in-u lottery drops out
of the joint CE. Runs cons-units and negativeUtils cases; compares V exactly and savings-policy
VALUES (PolicyInd2Val on both sides).
Cross tests run all three EZ cases via a for ezcase loop (matching the main EZ bank's convention),
except degenerateu (cases 1 and 3).

OMITTED: the ExpAssetu cross family (riskyasset == degenerate experienceassetu). Needs
EZ-ExpAssetu solvers, which do not exist (separate project). RESTORE when they exist.

Also skipped relative to the vNM bank: SimPanelValues / AggVars / ValuesOnGrid smoke checks —
those commands are Policy-driven and preference-free, already covered by the vNM bank.

## Test-first state (updated as toolkit work lands)

Existing toolkit pieces (NEVER GPU-run; carry the Change B outer-ezc1 removal and two parse-error
fixes from eb44124b — even the 'existing' shapes are unverified):
- `ValueFnIter/FHorz/RiskyAsset/EpsteinZin/`: 8 base raws, {d1,nod1} x {noa1,a1} x {noe,e}, all
  with-z. `ValueFnIter_FHorz_EpsteinZin_RiskyAsset(.m/_semiz.m)` dispatchers + 1 semiz raw
  (`EpsteinZinSemiExo/..._nod1_semiz_raw.m`).

IMPLEMENTATION STATUS (2026-08-18): the NOSEMIZ half is CODE-COMPLETE, awaiting GPU:
1. [DONE] noz shapes — 8 new EZ base raws in EpsteinZin/ + dispatcher stub replacement (figs 1,2,5,6
   + withA1 counterparts; semiz noz still pending).
2. [DONE] DC / GI / DC+GI tiers for withA1 EZ riskyasset — EpsteinZin/DivideConquer/ (dispatcher+8),
   GridInterpLayer/ (dispatcher+8), DivideConquerGridInterpLayer/ (dispatcher+8); tier routing
   added to ValueFnIter_FHorz_EpsteinZin_RiskyAsset.m; GI/DC_GI interpolate the TRANSFORMED EV
   before ^ezc6 (exact collapses under GI).
3. [NOA1 SEMIZ DONE 2026-08-18, awaiting GPU (figs 9-16)] semiz noa1 half implemented: the 7
   on-disk noa1 semiz raws audited against the full inherited-bug checklist (all verified CLEAN —
   the killed agents' grafts were post-fix-consistent) + noa1_e_semiz_raw created (~700 lines,
   full 0/1/2/3 ladders); `_semiz` dispatcher rewritten: 8 noa1 shapes routed (verified
   arg-for-arg vs on-disk signatures; UnKron3/4 noa1 tail mirroring vNM, fixing the old tail's
   missing-n_semiz e-reshape), all 8 withA1 semiz shapes now error cleanly ('unverified, carry
   known bugs, do not route'), DC/GI guarded; FromPolicy
   `ValueFnFromPolicy_FHorz_EpsteinZin_RiskyAsset_SemiExo.m` created (joint CE over
   (u,semiz',z',e'); e->z->per-d4-semiz->lookup->u ordering; a2-only WG; isnan clears incl. two
   the vNM lacks) and dispatched from the EZ FromPolicy parent (plain riskyasset unchanged).
   Post-GPU fix (run 4, 2026-08-18): fig 9 solved but crashed in the dispatcher UnKron — the 8
   noa1 semiz raws ended with a collapse-to-joint-index line (old exemplar convention) while the
   new dispatcher tail expects the vNM multi-row form; collapse lines removed, raws now return
   Policy3/Policy4 rows directly (the GPU-proven vNM pairing with UnKron3/4).
   GPU run 5 (2026-08-18): 385 checks ALL PASS — FIGS 9-16 (noa1 semiz) COMPLETE, incl. the semiz
   FromPolicy and all lowmemory rungs; run then reached fig 17 (first-ever withA1 execution) and
   crashed in the DC leg: the top EZ dispatcher built n_d1=[] (empty) for nod1 shapes instead of
   vNM's guarded n_d1=0, so the tier dispatchers' prod([])=1 misrouted nod1 into the with-d1 raw.
   Fixed by adopting the vNM guard block for n_d1/n_d2/n_d3.
   GPU run 6 (2026-08-18): 714 checks ALL PASS — figs 1-16 + the withA1 figs through d1_z_noe
   (incl. first successful DC/GI/DC+GI tier executions); crashed at the withA1 e-shape in
   nod1_e_raw:547: the EV u-lottery z-offset built as a row vector `aprimeIndex+N_a*((1:N_z)-1)`
   collides with the N_u dim — the GPU-validated pattern flattens with aprimeIndex(:), tiles
   aprimeProbs by N_z, and reshapes after masking. 10 blocks fixed across the 4 withA1 base raws
   (4 each in the e-raws incl. their main loops; 1 each in raw/nod1_raw where it lurked in the
   V_Jplus1 terminal branch).
   GPU run 7 (2026-08-18): 888 checks ALL PASS — FIGS 1-24 COMPLETE (noa1 nosemiz+semiz, withA1
   nosemiz incl. all DC/GI/DC+GI tier legs, FromPolicy(-GI), collapses exact under GI). Run stops
   at fig 25's deliberate withA1-semiz pause error. EVERYTHING IMPLEMENTED IS NOW GPU-VALIDATED.
   Cross tests not yet reached (they sit after fig 32 in the driver).
   *** GPU RUN 8 (2026-08-19): FULL BANK PASS — 1750 checks, zero errors, figs 1-32 + ALL cross
   tests (zase x8, semizasz x4, plainvswithA1 x4, degenerateu, d2recon), first-ever execution of
   the withA1-semiz half incl. its DC/GI/DC+GI tiers and FromPolicy semiz(-GI). Collapses at
   ~1e-16 throughout. THE EZ RISKYASSET PROJECT'S TEST BANK IS FULLY GREEN. ***
   Test-side note from the run: the d2recon guard's own diagnostic reports dsemiz is CONSTANT
   under the EZ calibration (all 3 cases) — the regression guard is vacuous here; strengthening
   the semiz/portfolio coupling in that cross test would make it bite (same applies to checking
   the vNM d2recon calibration).
   (run 7 aftermath, for the record: one intermediate run falsely appeared to regress figs 1-2 —
   cause was a stale .claude/worktrees checkout on the GPU machine shadowing the toolkit via
   genpath, not a code change.)

   FINAL BATCH COMPLETE 2026-08-18 — THE EZ RISKYASSET IMPLEMENTATION IS DONE (awaiting GPU of
   figs 25-32 + cross tests): withA1 semiz raws repaired (nod1_semiz needed 10 fix classes; the
   rest mostly WG/lottery; noz_e created) and all 8 now return multi-row Policy; the _semiz
   dispatcher routes all 16 shapes + the 3 semiz tiers (UnKron5/4 withA1 tail); the 3 semiz tiers
   built (EpsteinZinSemiExo/{DivideConquer,GridInterpLayer,DivideConquerGridInterpLayer}/:
   dispatcher + 8 raws each, 40-arg signatures token-exact vs the dispatcher calls); FromPolicy
   semiz-GI created + dispatched. Also fixed (EZ-side): the nosemiz raws' terminal-V_Jplus1
   warm-glow shaping drift. Final sweep: 93 files, zero non-baseline lint. New vNM finds recorded: RiskyAssetSemiExo V_Jplus1 ladders shorter than
   main-loop ladders in 5 noa1 raws (V(:,..,N_j) silently left 0 under high lowmemory+V_Jplus1);
   FromPolicy RiskyAsset_SemiExo missing pi_e isnan clear.
4. [DONE, nosemiz] ValueFnFromPolicy: ValueFnFromPolicy_FHorz_EpsteinZin_RiskyAsset(.m/_GI.m) +
   dispatch from ValueFnFromPolicy_FHorz_EpsteinZin (semiz+riskyasset errors cleanly).
5. [DONE, nosemiz] full lowmemory ladders mirrored from the vNM sources throughout.
Verification: 61-file lint sweep clean vs baseline; dispatcher raw-calls verified token-exact
against on-disk signatures in all three tiers; matrix-creator calls byte-identical to vNM (one
deliberate systematic exception: EZ terminal layer-2 interp flag 3 vs vNM 2, needed for the
refined warm-glow — all siblings consistent).

Shared commands expected to need NO changes: StationaryDist / AllStats / LifeCycleProfiles /
PolicyInd2Val (Policy-driven, preference-free).

## Policies (carried over from the EZ-SemiExo round)

- House style: no local helper functions; inline everything; full lowmemory ladders.
- **vNM riskyasset bug list**: any NEW bugs found in the vNM riskyasset files while implementing
  are LISTED below (file, line, symptom, correct pattern), NOT fixed in place; batch corrections
  proposed after the EZ-riskyasset implementation is finished.

### NEW vNM findings post-batch (MAIN family, not riskyasset; found during the main-family
### EZ-SemiExo step-2/3 build 2026-08-19; note the main bank's own md files
### appear to have been moved/removed, so recorded here)
### MINI-BATCH APPLIED 2026-08-25 (user approved): all five finds A-E fixed in the vNM files,
### concept order D (offset), E (flag 3->6 + noz_e shiftdim removals), A/B (isnan clears:
### pi_e in both FromPolicy files + 4 EVnext_atpolicy clears in _GI), C (lowmem1 4-D prealloc).
### All lint-clean vs baseline. Awaiting vNM GPU validation run; note D/E sites are in blocks
### the current vNM bank never reaches (validated indirectly via the GPU-green EZ mirrors).
A. `ValueFnFromPolicy/SemiExo/ValueFnFromPolicy_FHorz_SemiExo_GI.m` — 0*(-Inf) NaN class (same
   as applied batch items 1/8/10): no isnan clear after the pi_e e'-integration (lines 173-174;
   own z'/semiz' sums have clears) and none after the L2-interpolated weighted sums (217/227/
   241/254; weight exactly 0 whenever L2=1 times a -Inf endpoint → NaN into V unmasked).
B. `ValueFnFromPolicy/SemiExo/ValueFnFromPolicy_FHorz_SemiExo.m:156` — same missing pi_e clear
   (own clears at 169-170/177-178/185-186). The new EZ SemiExo(-GI) FromPolicy files include all
   these clears.
C. `ValueFnIter/FHorz/SemiExo/DivideConquer/ValueFnIter_FHorz_SemiExo_DC1_e_raw.m:47-49` —
   lowmemory==1 preallocation of V_ford2_jj/Policy_ford2_jj omits the N_e dimension but is
   indexed (...,e_c,d2_c): auto-grow hides it when N_d2<=N_e; errors when N_d2>N_e. Same class
   as the base-SemiExo e-raw preallocation bug fixed (with approval) in the Class A round —
   evidently also present in the DC tier. (The new EZ DC1_e raw uses a correct 4-D prealloc.)

D. `ValueFnIter/FHorz/SemiExo/DivideConquerGridInterpLayer/ValueFnIter_FHorz_SemiExo_DC1_GI1_nod1_noz_raw.m:72,127`
   — terminal no-V_Jplus1 blocks use the no-decision COLUMN offset `loweredge+(0:1:maxgap(ii))'`
   where loweredge is [N_d2,1,1,N_semiz] (d2 in dim1, full-n_d2 flag-3 creator) → dim1 clash
   whenever maxgap(ii)>0; stale "1-by-1-by-n_z" comment at :71. Found 2026-08-25 while fixing the
   identical (latent, warm-glow-triggered) drift in the EZ mirror
   `EpsteinZinSemiExo/DivideConquerGridInterpLayer/..._DC1_GI1_nod1_noz_raw.m:116,185` (FIXED, EZ
   only). In vNM the terminal maxindex1 comes straight from the return matrix (no warm-glow gate),
   so maxgap>0 is reachable — implying the vNM nod1_noz DC_GI no-V_Jplus1 leg has never been
   executed by any test. Rest of the vNM SemiExo DC/GI/DC_GI tier verified row-form correct both
   directions. Fix: row offset `(0:1:maxgap(ii))` as in its own siblings (e.g. DC1_GI1_nod1_raw).

E. vNM SemiExo DCGI e-raws pass row-form aprime arrays to CreateReturnFnMatrix_Disc_DC1_e with
   flag 3 — but the _e creator's flag 3 is the COLUMN-form flag (N_aprime=size(...,1)); the
   row-form analog is flag 6. Errors at the creator reshape when maxgap+1~=N_d2, and reshapes
   SILENTLY WRONG when maxgap+1==N_d2. All in `ValueFnIter/FHorz/SemiExo/DivideConquerGridInterpLayer/`:
   DC1_GI1_e_raw.m:110,170,234,296; DC1_GI1_nod1_e_raw.m:99,155,215,275;
   DC1_GI1_noz_e_raw.m:93,151,211; DC1_GI1_nod1_noz_e_raw.m:86,142,202. (vNM DivideConquer/ and
   GridInterpLayer/ e-raws have no flag-3 sites — clean.) Found 2026-08-25 fixing the identical
   drift in the EZ mirrors (14 sites flag 3→6 across the 4 EpsteinZinSemiExo DCGI e-raws, plus 3
   stray shiftdim(...,1) midpoint-write wrappers removed in the noz_e file; FIXED, EZ only).
   Fix: flag 3→6 at the listed sites (non-e creator's flag 3 IS row-form — do not touch non-e raws).

### vNM RiskyAsset bug list — BATCH APPLIED 2026-08-19 (all 10 items; see
vNM_RiskyAsset_batchcorrections_proposal.md for the per-item before/afters). Applied in concept
order: crashes (6+7, 3, 4), silent-wrongs (2, then 1/8/10 isnan clears), silent-zeros (5: five
terminal-ladder guard widenings ==1 -> >=1; 9: eight GI-noz guard widenings), cosmetics bundled
per touched file. 22 touched files lint clean vs baseline. Item 2 verified behaviour-preserving
for the current bank (terminal corner solution); everything else lives in unexercised branches.
Leftover-A ladder sweep 2026-08-19: family-wide fall-through audit (818 ladders parsed, 224
files): 290 sites widened in 120 files (272 ladder finals ==k -> >=k, verified most-looped rung
each; 18 support gates incl. 12 special_n_semiz gates the item-9 batch had left stale — those
would have crashed at lm3). Includes the base SemiExo noz raws + all EZ mirrors + the DC2A/GI2A
with2A-session files (one-token edits, flagged to user). 149 at-family-max ladders + 2 local
midpoint ladders + 47 sub-max shape gates verified not-a-bug. Lint byte-identical on all 120.
Bank behaviour unchanged by construction (only previously-dead lowmemory values gain a branch).

vNM bank rerun 2026-08-19: CONFIRMED behaviour-preserving — 600 checks all exact zeros through
figs + all cross tests, loose-check max byte-identical to the pre-fix baseline (0.0239...). (The
run later crashed in a with2A1 fig's SimPanel section — unrelated concurrent project, per user.)

1. `ValueFnFromPolicy/RiskyAsset/ValueFnFromPolicy_FHorz_RiskyAsset_GI.m` — per_u sums at lines
   226, 241, 257, 274 (and the e-collapse at line 204) have no isnan clear: zero interpolation
   weights (w_a1_upper=0 whenever L2=1) times -Inf corner values give 0*(-Inf)=NaN which
   propagates into V for all earlier ages. Correct pattern: the plain file's
   `EVnext_atpolicy(isnan(...))=0` after the pi_u sum (ValueFnFromPolicy_FHorz_RiskyAsset.m
   lines 237/254). Also dead code: GI lines 122-126 preallocate a1_lower that line 132 overwrites.
   (The new EZ FromPolicy GI file deliberately includes the isnan clears.)

2. **Terminal-period Policy decode uses wrong N_d (wrong whenever N_d2>1)** — the no-V_Jplus1
   branch maxes a ReturnMatrix whose first dim excludes d2 but decodes with N_d including d2:
   `Raw/ValueFnIter_FHorz_RiskyAsset_nod1_raw.m` lines 52-55, 65-68 (first dim N_d3*N_a1, decode
   uses N_d=N_d2*N_d3; correct: rem(maxindex-1,N_d3)+1 and ceil(maxindex/N_d3));
   `Raw/..._nod1_noz_raw.m` lines 45-48 (same, plus cosmetic ceil(dindex) no-op);
   `Raw/..._raw.m` lines 53-57, 67-71 (first dim N_d1*N_d3*N_a1, decode uses N_d1*N_d2*N_d3;
   correct pattern is noz_e_raw lines 52-56 which uses N_d1*N_d3);
   `Raw/..._noz_raw.m` lines 46-50 (same as raw.m).
3. **Stray extra n_e argument (runtime 'too many input arguments')** in the terminal-V_Jplus1
   lowmemory==1 e-loop of all four noz_e raws: `Raw/..._nod1_noz_e_raw.m:119`,
   `Raw/..._nod1_noa1_noz_e_raw.m:106`, `Raw/..._noa1_noz_e_raw.m:110`, `Raw/..._noz_e_raw.m:132`
   — each passes `..., n_e, special_n_e, ...` (9 args) to the 8-arg CreateReturnFnMatrix_Case2_Disc.
   Correct pattern: the main-loop call in the same file (drop the n_e).
4. **V_Jplus1 reshaped to [N_a2,1] instead of [N_a,1]** — `Raw/..._nod1_noz_raw.m:53` and
   `Raw/..._noz_raw.m:53`; errors whenever N_a1>1 (all other 14 raws use [N_a,...] correctly,
   and aprimeIndex ranges to N_a).

5. `ValueFnIter/FHorz/RiskyAsset/RiskyAssetSemiExo/` — V_Jplus1 terminal ladders SHORTER than the
   main-loop lowmemory ladders in 5 noa1 raws (nod1_noa1: V_Jplus1 does 0/1 vs loop 0/1/2;
   noa1_e + nod1_noa1_e: 0/1 vs 0/1/2/3; noa1_noz_e + nod1_noa1_noz_e: 0/1 vs 0/1/2): under
   V_Jplus1 + a high lowmemory the terminal age is silently never solved (V(:,...,N_j) stays 0,
   no error). noa1_raw covers 0/1/2 fully; the EZ semiz raws mirror this vNM control flow
   faithfully per the list-don't-fix policy.
6. `RiskyAssetSemiExo/ValueFnIter_FHorz_RiskyAssetSemiExo_nod1_raw.m:97` — V_Jplus1 reshaped
   [N_a,N_z], should be [N_a,N_bothz].
7. Same file (terminal V_Jplus1 lm0/lm1, ~132-138/~183-188) and `..._RiskyAssetSemiExo_raw.m`
   (~127-137/~178-188): unflattened EV u-lottery `EV(aprimeIndex+N_a*((1:N_bothz)-1))` with matrix
   aprimeIndex (dim clash with N_u) + aprimeProbs missing the N_bothz expansion — the same files'
   main loops use the correct flattened `aprimeIndex(:)` form.
8. `ValueFnFromPolicy/RiskyAsset/ValueFnFromPolicy_FHorz_RiskyAsset_SemiExo_GI.m` — same missing
   pi_e isnan clear as item 6-below (lines 224-227); dead identical if/else at 118-122; unused
   l_a2/l_aprime/a1_grid.
9. `RiskyAssetSemiExo/GridInterpLayer/` noz-family lowmemory fall-through (silent zeros):
   GI1_noz_raw/GI1_nod1_noz_raw ladders end at `elseif lowmemory==1` and the noz_e terminal
   ladders at `elseif lowmemory==2`, with no final else and no upstream range validation — an
   out-of-range lowmemory (values the with-z siblings accept) executes no branch and silently
   leaves V at zeros. (EZ tier mirrors verbatim per policy, so inherits the behaviour.) Also
   lint hygiene: unused a2ind/special_n_d4/zind assignments and stale %#ok pragmas across 6 files.
10. `ValueFnFromPolicy/RiskyAsset/ValueFnFromPolicy_FHorz_RiskyAsset_SemiExo.m:215` — missing isnan
   clear after the pi_e sum (`V_next=sum(V_next.*shiftdim(pi_e_J(:,jj+1),-2),3);` with no
   `V_next(isnan(...))=0`): 0*(-Inf) NaNs propagate into the z/semiz sums where downstream clears
   then zero entries that should be -Inf/finite. Every other probability-weighted sum in the file
   has its clear. (The new EZ FromPolicy SemiExo file includes the clear.) Cosmetic in same file:
   dead identical if/else at 119-123; unused l_a/N_a2/a1_grid/l_a2/l_aprime.

EZ-side latent find — FIXED 2026-08-18: the terminal-V_Jplus1 warm-glow shaping in raw/nod1_raw/
nod1_e_raw used `repmat(WGmatrix,1,N_a[,N_z[,N_e]])` (dim-2 N_a mismatches temp4's
[N_d23*N_a1,1,N_z] and breaks the WGmatrix(becareful) logical indexing); rewritten to each file's
GPU-validated main-loop shaping (`WGmatrix.*ones(1,1,N_z)` at the vectorized rungs, column at the
per-z rungs; no e dim since the e-expectation precedes). e_raw already had the correct form.
V_Jplus1+warmglow branch only (untested by the bank).

### EZ-side bugs found in the PRE-EXISTING EZ riskyasset raws — ALL FIXED 2026-08-18
(items 1-6 below plus two more found while fixing: the WG composite-index bug (item 2) was in ALL
FOUR withA1 raws (raw, nod1, e, nod1_e; both the top and in-loop WG blocks — rewritten to index
by a2primeIndex only, lottery at [N_d23,N_u], then repmat to the (d,a1prime) rows, with the isnan
clears the e-raws were also missing); the missing ^ezc8 (items 3/5) was at all 6 main-loop temp4
sites in EACH e-raw (12 lines fixed); the V_Jplus1 [N_a2,...] reshape (item 4) was in all four
withA1 raws incl. both e-raws; and nod1_e_raw:44 had the same eind N_z-vs-N_e bug as e_raw:46
(dead code there, fixed anyway). All four files checkcode-clean apart from benign
now-dead composite-index assignments left in place.)

Second fix round 2026-08-18 (found by the tier-build agents' sweeps): the WG composite-index bug
had also been inherited by the 4 NEW noz withA1 base raws (grafted from the pre-fix exemplars) —
all 8 blocks fixed with the same a2-only pattern; and the missing-^ezc8 class was ALSO present in
the two PRE-EXISTING noa1 e-raws (noa1_e, nod1_noa1_e; 6 sites each) — fixed. Sweep confirms zero
remaining missing-ezc8 sites anywhere in EpsteinZin/ or EpsteinZinSemiExo/.

Third fix round 2026-08-18 (triggered by the first GPU run of the EZ bank, which crashed at fig 3's
lowmemory=1 leg with 'Unrecognized function or variable d2index' in nod1_noa1_raw:149): full
structural audit of all 16 nosemiz base raws' terminal/lowmemory branches; 17 fixes in 5 files
(11 files clean):
- nod1_noa1_raw + noa1_raw: terminal no-V_Jplus1 lm1 z_c-loop shared Policy assignments referenced
  d2index (and in noa1_raw also d1index) which only exist under warmglow==1 — restructured per the
  lm0 pattern (the GPU crash); plus out-of-bounds `+N_d3*zind`-style offsets on per-z_c refines in
  the V_Jplus1/jj lm1 branches (would be the NEXT crash); plus undefined n_d/d_grid in jj-lm1
  CreateReturnFnMatrix calls.
- raw + nod1_raw: `Policy(row,:,:,...)` all-z writes inside z_c loops -> `(row,:,z_c,...)`
  (size-mismatch crash class).
- raw + e_raw: terminal warmglow d1index lookups missing the N_a1 stride factor (SILENT
  wrong-answer bugs, in bounds).
All 16 files lint-identical to baseline after edits. GPU note: many of these live in V_Jplus1 and
lowmemory branches the bank does not currently exercise (no V_Jplus1 leg, warmglow untested) —
same systematic gap as recorded for the main bank.

Fourth fix round 2026-08-18 (second GPU run: 141 checks pass — figs through d1_z_e's basic legs —
then 'Index exceeds matrix dimension' in CreateReturnFnMatrix_Case2_Disc_e from noa1_e_raw:94):
four pre-existing EZ raws passed raw stacked GRIDS to the ReturnFn matrix creators where the vNM
counterparts pass joint GRIDVALS (the creators index d_gridvals(:,2) whenever there are 2+ d/a
variables): noa1_e_raw (9 sites; missing d13_gridvals), e_raw (9; missing d13a1_gridvals +
a12_gridvals), nod1_e_raw (9; missing d3a1_gridvals + a1a2_gridvals), raw (6; missing
d13a1_gridvals + a12_gridvals). Fixed: CreateGridvals added to each preamble, all 33 call sites
swapped to the vNM convention. nod1_raw/noa1_raw and all NEW noz/tier raws already used gridvals.
All four files lint clean.

Composite-WG audit 2026-08-18 (final): base raws fixed (a2-only+repmat); DC tier raws are CORRECT
via a different route (they pre-expand WGmatrix with repelem(WGmatrix,N_a1,1) before the composite
indexing — in-bounds and equivalent); GI/DC_GI tiers clean by construction; noa1 files correct
as-is (aprimeIndex IS the a2 index there). REMAINING (for the paused semiz batch): the 7 withA1
semiz raws in EpsteinZinSemiExo/ — including the PRE-EXISTING nod1_semiz_raw that served as the
exemplar — have composite WG indexing with NO expansion guard (OOB when N_a1>1 with warmglow).
The semiz batch must fix these 7 (plus whatever the 2 missing raws need) when it resumes.

In `ValueFnIter/FHorz/RiskyAsset/EpsteinZin/ValueFnIter_FHorz_RiskyAsset_EpsteinZin_e_raw.m`:
1. Line 46: `eind=shiftdim(0:1:N_z-1,-2)` should use N_e — wrong d1index lookups whenever
   N_e~=N_z (used at lines 131, 313, 586).
2. Warm-glow blocks (lines 72-90, 478-507): WGmatrix built over n_a2 (length N_a2) but indexed
   with the composite aprimeIndex (values up to N_a1*N_a2) — out of bounds whenever N_a1>1 with
   warmglow. (The noa1 raw indexes with a2-only indices and is fine.)
3. Line 561 vs noa1 raw line 370: the jj<N_j warmglow temp4 omits `.^ezc8`
   (`(sj*temp4+(1-sj)*WGmatrix).^ezc6`), inconsistent with its own N_j block (line 287).
4. `EpsteinZin_nod1_raw.m` / `EpsteinZin_raw.m`: V_Jplus1 reshaped to [N_a2,N_z] instead of
   [N_a,N_z] (same class as vNM bug 4; errors whenever N_a1>1 with V_Jplus1).
5. `EpsteinZin_nod1_e_raw.m`: main-loop temp4 drops the `.^ezc8(jj)` factor that its own terminal
   block and all no-e exemplars apply (same class as EZ bug 3 in _e_raw).
6. `EpsteinZin_nod1_e_raw.m` lines 131-141: terminal no-V_Jplus1 lowmemory==0 branch writes
   Policy rows 3/4 (silently growing the 3-row Policy) when warmglow==1 and misassigns rows when
   warmglow==0.

## Still to add (later)

- 2a1 EZ riskyasset tests (two standard endogenous assets + risky asset).
- ExpAssetu cross family once EZ-ExpAssetu exists.
- No coverage for sj/warm-glow/EZmortalityriskaversion/EZutils=2 (same gap as the main EZ bank).
