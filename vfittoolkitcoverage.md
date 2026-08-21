# VFI Toolkit — test coverage

Coverage of the FHorz core test banks. Tier is assigned from each subcode's directory
(`Noa1_subcodes` / `With2A1_subcodes` / else), not its filename — filenames misattribute
cross-tests. Format in the tier columns is `variants + cross-testsx`.

Last updated: 2026-08-19

## CoreFHorzTests + the ExpAsset family + RiskyAsset

| bank | noa1 | withA1 | with2A | var | x-tests | panel vs dist | V_Jplus1 | core raws | QH | EZ |
|---|---|---|---|---|---|---|---|---|---|---|
| CoreFHorzTests | — | 16+6x | 16+4x | 32 | 10 | 32/32 | 32 + 32 QH | 142 | 42 | 24 |
| ExpAsset | 16+10x | 16+14x | 16+2x | 48 | 26 | 48/48 | 0 | 128 | 48 ⚠ | — |
| ExpAssetU | 16+10x | 16+14x | 16+2x | 48 | 26 | 48/48 | 0 | 128 | none | — |
| ExpAssete | 8+16x | 8+8x | 8+2x | 24 | 26 | 24/24 | 0 | 64 | 24 | — |
| ExpAssetz | 8+16x | 8+8x | 8+4x | 24 | 28 | 24/24 | 0 | 64 | 24 | — |
| ExpAssetze | 4+4x | 4+16x | 4+2x | 12 | 22 | 12/12 | 0 | 32 | 12 | — |
| ExpAssetsemiz | 8+4x | 8+6x | 8+2x | 24 | 12 | 24/24 | 0 | 64 | none | — |
| RiskyAsset | 16+0x | 16+25x | 16+4x | 48 | 29 | 48/48 | 0 | 128 | none | 50 |
| **total** | | | | **260** | **179** | **260/260** | **64** | | **150** | **74** |

439 subtests in the main banks, plus 150 QH and 74 EZ = **663**.

**The RiskyAsset with2A tier is complete and GPU-green** (2026-08-19).

All 48 2A raws exist — `{DC2A, GI2A, DC2A_GI2A}` × `{nosemiz, SemiExo}` × 8 shock/decision
combinations — plus 2A routing in all six tier dispatchers. Figs 33–48 ran **754 checks with every
exact check zero**: 64 DC2A, 32 GI2A, 64 DC2A_GI2A, 432 lowmemory-ladder, and the 4
degenerate-`a1_2` cross-tests. The only non-zeros are the 32 `StationaryDist with/without grid
interp … close to zero` convergence lines, which are not exact checks.

This closes the gap that used to be the entire difference between ExperienceAsset's core raws and
RiskyAsset's: both families now stand at **128 core raws**.

Supporting toolkit changes made along the way, all of them general rather than RiskyAsset-specific:

- `CreateReturnFnMatrix_Case2_Disc{,_e,_noz}` extended from `l_d<=4` to `l_d<=6`. The *base* raws
  pack `a1` into the d slot, so each extra `a1` dimension costs a d slot and
  `[d1,d3,d4,a1_1,a1_2]`=5 was blocked. (`l_d==6` is still unreachable until
  `SubCodes/CreateGridvals.m` gains an `l_x==6` case.)
- `EvalFnOnSimPanelIndex` gained an `l_daprime==6` block plus an `else` guard. `l_daprime` is
  `length(n_d)+length(n_aprime)`, so any model with four decision variables and two chosen assets
  reaches 6; previously it fell through with `Values` unassigned, which surfaced as a confusing
  error from the caller rather than an unsupported-case message.
- `UnKronPolicyIndexes6_FHorz_{z,z_e}` added, and all six 2A tiers now store `a1_1prime` and
  `a1_2prime` as **separate** Policy channels. Two SemiExo tiers briefly used a joint channel;
  both conventions produce the identical unkronned Policy, so this was consistency, not a fix.
- `level1n=min(level1n,n_a1)` → `n_a1(1)` in the 4 DC/DC_GI dispatchers (a vector `level1n` was a
  real bug), and `if n_a1>0` → `if prod(n_a1)>0` in `ValueFnIter_FHorz_RiskyAsset.m`.

`ValueFnFromPolicy` needed no 2A work — it was already generic.

Directory-name trap when recounting: the baseline uses `With2A_subcodes`, the ExpAsset family
uses `With2A1_subcodes`. A tier rule that matches only one spelling silently collapses the
baseline's 16 `with2A` variants into `withA1`.

### Column meanings

- **noa1 / withA1 / with2A** — the three asset tiers. `noa1` does not apply to the baseline.
  Variant counts are uniform within a bank by construction: 16 per tier where three shock
  axes are free (`{nod,d} × {noz,z} × {noe,e}`, doubled for semiz), 8 where the family's
  driving shock fixes one axis, 4 for ExpAssetze where both `z` and `e` are mandatory.
- **panel vs dist** — variants that simulate a panel *and* compare it against the agent
  distribution (`AgeConditionalStats` / `AllStats`). Cross-tests are excluded: they assert
  bit-exact V/Policy/Dist agreement, which a Monte-Carlo panel cannot improve on.
  These are printed side-by-side rows, eyeballed, not thresholded assertions.
- **V_Jplus1** — subcodes that exercise the terminal-condition branch.
- **core raws** — solver raws in the family, excluding exotic-preference variants
  (which the QH/EZ mirrors cover). This is the surface the main bank dispatches into.
- **QH / EZ** — subtests in the bank's `withQuasiHyperbolicDiscounting` /
  `withEpsteinZinPreferences` mirror. `—` means no mirror; `⚠` see below.

### Scope of the exotic-preference mirrors

**The QH and EZ mirrors only need to test `V` and `Policy`.** That is the design principle,
and it is why the blank `panel vs dist` entries for those mirrors are not a coverage hole.

Alternative preferences change how the value function is formed, and so they change `Policy`.
They change nothing downstream of it. Every other model object — the agent distribution, the
panel simulation, and all model statistics — is computed *from* `Policy` by code that never
sees `vfoptions.exoticpreferences`. So once a QH or EZ subtest has established that `V` and
`Policy` are right, the whole downstream chain is already covered by the baseline
exponential-discounting / von-Neumann-Morgenstern banks, which exercise that same code on the
same shapes. Adding dist/panel/stats blocks to a QH or EZ subtest re-tests baseline code at
exotic-preference prices rather than extending coverage.

The exception is `ValueFnFromPolicy`, which *does* have to know the preference type — hence
its heavy use in the QH mirror (141 of 155 files), including the `Valt` reconstruction oracle
for Naive.

#### What the mirrors actually contain today

Panel simulation is absent everywhere — `SimPanelValues` appears in 0 of 155 QH and 0 of 220
EZ subcodes. Dist-derived output is not absent, but it is confined to two banks (TPath banks
excluded throughout):

| mirror / bank | files | asserted | displayed |
|---|---|---|---|
| QH — CoreFHorzTests | 43 | 120 | 512 |
| QH — ExpAsset, ExpAssete, ExpAssetz, ExpAssetze | 112 | **0** | **0** |
| EZ — CoreFHorzTests | 73 | 120 | 384 |
| EZ — RiskyAsset | 147 | 100 | 240 |

*asserted* = `fprintf` checks whose compared value is dist-derived. *displayed* = unsuppressed
`[...]` rows that dump moments to screen with no assertion.

**The four ExpAsset-family QH mirrors already follow the principle exactly** — no dist, no
stats, no display, only `V` / `Policy` / `ValueFnFromPolicy`. The older baseline mirrors and
the RiskyAsset EZ bank predate it.

What the asserted checks compare:

- **`StationaryDist` across solver tiers** (64 QH / 96 EZ) — `StationaryDist with/without grid
  interp, this should be close to zero`, on a deliberately big `a_grid`. Base vs GI and DC vs
  DC+GI, i.e. a solver-tier comparison, not a validation of the distribution code.
- **`StationaryDist` in cross-tests** (56 QH / 112 EZ) — `z as e`, `z and e 1/2`, `semiz as z`,
  `single-point z vs no z`, `(z,e) vs merged joint-markov`, `plainvswithA1`, `semizasz`. Two
  formulations that should coincide, compared on the dist as well as on `V`/`Policy`.
- **`AllStats.Mean`** — only 12 lines, all EZ `CrossTest plainvswithA1`.
- **`LifeCycleProfiles` / `AgeConditionalStats` are never asserted anywhere.** They feed
  figures (160 QH / 312 EZ `plot`/`subplot` calls) and the displayed rows.

The displayed rows are four shapes in QH, repeated 128× each — `[AllStatsN.assets.Mean,
AllStatsM.assets.Mean]`, `[AllStatsN.earnings.Gini, …]`, `[AgeConditionalStatsN.earnings.Mean;
…]`, `[AgeConditionalStatsN.assets.StdDeviation; …]` — under the heading `With/without grid
interp, should get much the same moments (for big a_grid)`. EZ has ten shapes, adding `a1`/`a2`
variants for the two-asset tier. They reach the diaries: `CoreFHorzQHTestsdiary.txt` carries
512 `ans =` blocks against 64 `close to zero` lines.

**These are not covered by the "Policy determines everything downstream" argument, and should
stay.** They are a grid-interpolation convergence check on the exotic-preference solver: with a
big enough `a_grid`, the GI tier should reproduce the brute-force answer. `V` and `Policy`
cannot carry that check — GI returns policies off the coarse grid, so the two tiers do not
agree elementwise by construction — which is exactly why the comparison is made at the moment
level. The moments are the readout, not the object under test. That makes this a genuine
exotic-preference-specific test: it asks whether *this* solver's GI tier converges.

The caveat is that the moment rows are displayed, not asserted — no `fprintf`, no tolerance,
so a regression is caught by eye or not at all. Only the accompanying `StationaryDist …
close to zero` line is a real check. Worth converting the moment rows to tolerance checks if
the eyeballing ever becomes a problem.

The cross-test `StationaryDist` comparisons (56 QH / 112 EZ) are the group that *is* largely
redundant under the principle above: the two formulations agree exactly on `V` and `Policy`
on the same line, so the dist adds little. Harmless, but new subtests need not repeat it.

### Raw counts behind the "core raws" column

Recounted 2026-08-18 against the **working tree**, which includes 124 uncommitted raws.
ExperienceAssetz re-counted 2026-08-19 after the QH semiz work landed (commit `2f291d33`):
168 total, 104 QH (40 nosemiz + 64 semiz), 64 core.

| family | total raws | QH | EZ | core |
|---|---|---|---|---|
| baseline (FHorz, excl. asset families) | 406 | 224 | 40 | 142 |
| ExperienceAsset | 128 | 0 | 0 | 128 |
| ExperienceAssetu | 128 | 0 | 0 | 128 |
| ExperienceAssete | 192 | 128 | 0 | 64 |
| ExperienceAssetz | 168 | 104 | 0 | 64 |
| ExperienceAssetze | 84 | 52 | 0 | 32 |
| ExperienceAssetsemiz | 64 | 0 | 0 | 64 |
| RiskyAsset | 208 | 0 | 80 | 128 |

Note the toolkit is inconsistent about where exotic-preference raws live: ExpAssete/z/ze
keep them in a `QuasiHyperbolic/` subdir of their own family, while the *baseline* QH raws
sit under `ExoticPrefs/QuasiHyperbolic/`.

### Known open items

- **V_Jplus1**: now 64 subcodes — all 32 in `CoreFHorzTests_subcodes` and all 32 in its QH
  mirror (committed `daa90d5`). Still zero coverage across the whole ExpAsset family and
  RiskyAsset, against 660
  ExpAsset-family raws that all carry a `V_Jplus1` branch (ExperienceAsset 128,
  ExperienceAssetu 128, ExperienceAssete 160, ExperienceAssetz 96, ExperienceAssetze 84,
  ExperienceAssetsemiz 64). Those branches are age-shifted copies of the in-loop code — the
  shape that produced the `jj`/`N_j` bug. Largest remaining gap on any axis.
- ⚠ **One QH bank is still test-first against an unwritten solver.** Status as of
  2026-08-19, read off the working tree rather than the banners:
  - `CoreFHorzQHExpAssetTests` (48 subtests) — genuinely unsupported, and the only remaining
    test-first bank. `ExperienceAsset` has **0 QH raws** and `ValueFnIter_Case1_FHorz` has no QH
    branch under `vfoptions.experienceasset>=1`. Every subtest errors at its first `ValueFnIter`
    call. Same for `ExperienceAssetu` and `ExperienceAssetsemiz` (no mirrors written). Closing
    this is solver work, not test work.
  - `CoreFHorzQHExpAsseteTests` (24 subtests) — **done and GPU-green.** The toolkit has 128 QH
    raws under `ExperienceAssete/QuasiHyperbolic/`, four `QuasiHyperbolicExpAssete*` dispatchers
    plus four SemiExo ones, and live routing at `ValueFnIter_Case1_FHorz.m:401`/`:411` (committed
    `b6bc5acf`). The 2026-08-18 run is green: **1624 checks, 24/24 figures, no errors** — Naive and
    Sophisticated on V / Valt / Policy / Policyalt across each lowmemory ladder, plus 24
    `ValueFnFromPolicy` oracle checks. Two checks print `0.00000001` rather than `0.00000000`
    (`Sophisticated ValueFnFromPolicy (Valt)` and `... (DC1, Valt)`); that is display-precision
    rounding on the reconstruction oracle, not a defect — but it does mean a plain
    grep-for-nonzero flags this diary, so read those two lines before concluding anything.
    One toolkit limit remains but is not reachable from this bank:
    `QuasiHyperbolicExpAsseteSemiExo_{DC,GI,DC_GI}` error on `N_a1==0`, and the noa1+semiz
    block (figs 5-8) is base-method only.
  - `CoreFHorzQHExpAssetzTests` (24 subtests) — **done and GPU-green as of 2026-08-19** (commit
    `2f291d33`). The stale banner on this bank should be removed. `ExperienceAssetz` now carries
    104 QH raws (40 nosemiz + 64 semiz), the `QuasiHyperbolicExpAssetzSemiExo` dispatcher and its
    three `{DC,GI,DC_GI}` sub-dispatchers all exist, and no `not yet implemented` error survives
    anywhere in the family. The run is **1268 checks, 24/24 figures, every one zero, no errors**,
    covering all four solver tiers x nod1/with-d1 x no-e/with-e x Naive/Sophisticated, including
    the with-e 2A paths at lowmemory 3 and the beta0=1 degeneracy checks that collapse QH onto
    the exponential. `ValueFnFromPolicy_FHorz_QuasiHyperbolic_ExpAssetz_SemiExo_GI` was added to
    close the last gap (ExpAssetz had been the only family missing a SemiExo_GI value-from-policy).
  - `CoreFHorzQHExpAssetzeTests` (12 subtests) — passes. No banner; GPU-validated, 776 checks,
    12/12 figures.
  - `CoreFHorzQHTests` (baseline, 32 subtests) — 4016 checks, none nonzero. Note this diary uses
    a different closing-marker style from the ExpAsset-family banks, so the per-figure counting
    used above does not apply to it.
  - `CoreFHorzTPathQHTests` — bank exists but has **never produced a diary**; status unknown.
- **QH**: no mirror at all for ExpAssetU, ExpAssetsemiz or RiskyAsset. For ExpAssetsemiz
  and ExpAssetU this is a *toolkit* gap — neither family has any QH raws — so closing it
  means solver code, not tests.
- **EZ**: mirrors exist only for the baseline and RiskyAsset. The entire ExpAsset family
  has no Epstein-Zin coverage, and the toolkit has no EZ raws for those families either.
- **Panel means are not reproducible across calls.** `SimPanelIndexes_FHorz_*` simulates each
  agent's lifecycle inside an unconditional `parfor`, so the shock draws happen on the workers
  and a client-side `rng(1)` does not reset their streams. Two calls on identical models give
  different panels. (The seedpoint draw *is* on the client and so is reproducible.) The noa1
  cross-tests originally printed a panel-mean difference under a `this should be zero` label;
  as of 2026-08-18 it is split out into its own `sim panel means should roughly match` line
  across 32 files / 56 sites, so the diaries are once again clean under a grep-for-nonzero.
  Anything comparing panels must be written as a tolerance check, never an exact zero.
- **QH banks run no panel simulation** at all — by design, see "Scope of the exotic-preference
  mirrors" above. Not a gap.
- **RiskyAsset with2A: DONE and GPU-green** (2026-08-19). All 48 raws
  (`{DC2A, GI2A, DC2A_GI2A}` × `{nosemiz, SemiExo}` × 8) plus 2A routing in all six dispatchers.
  Figs 33-48, 754 checks, every exact check zero. RiskyAsset core raws 80 -> 128, matching
  ExperienceAsset exactly; the two families now cover the same 8 shock/decision combinations
  across {base, DC, GI, DC+GI} × {nosemiz, SemiExo} with a 2A variant in each of the 6 DC/GI
  directories. Base dirs are unaffected: brute force Krons the a1 dimensions, so no separate 2A
  raw is needed there. See the status section near the top for the supporting toolkit changes.
- **RiskyAsset calls its cross-tests as bare `fn(...)`, not `output=fn(...)`.** Every other
  bank uses the `output=` form. Any cross-bank tally that greps `^output=` will silently
  report 0 cross-tests for this bank instead of 25.
- **Missing Epstein-Zin raw**: `RiskyAsset/EpsteinZinSemiExo` has 7 of the 8 noa1
  combinations — `RiskyAsset_EpsteinZin_noa1_e_semiz` (d1, z, e, semiz) is absent while all
  seven siblings exist. Looks like an oversight rather than a deliberate omission.
- RiskyAsset is the one bank that recently ran clean to completion.
- **OOM ceiling**: every ExpAsset bank's largest `d1_z_e` + 2A/semiz variant exceeds GPU
  memory. Where it lands matters: in ExpAsset/ExpAssetU it hits a *leading* brute-force
  baseline and takes the rest of the script with it; in ExpAssete/ExpAssetz/ExpAssetsemiz
  it hits the *trailing* big-grid block, after that subtest's exact checks have printed.
