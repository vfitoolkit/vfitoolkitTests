# VFI Toolkit — test coverage

Coverage of the FHorz core test banks. Tier is assigned from each subcode's directory
(`Noa1_subcodes` / `With2A1_subcodes` / else), not its filename — filenames misattribute
cross-tests. Format in the tier columns is `variants + cross-testsx`.

Last updated: 2026-08-18

## CoreFHorzTests + the ExpAsset family + RiskyAsset

| bank | noa1 | withA1 | with2A | var | x-tests | panel vs dist | V_Jplus1 | core raws | QH | EZ |
|---|---|---|---|---|---|---|---|---|---|---|
| CoreFHorzTests | — | 16+6x | 16+4x | 32 | 10 | 32/32 | 32 + 32 QH | 142 | 42 | 24 |
| ExpAsset | 16+10x | 16+14x | 16+2x | 48 | 26 | 48/48 | 0 | 128 | 48 ⚠ | — |
| ExpAssetU | 16+10x | 16+14x | 16+2x | 48 | 26 | 48/48 | 0 | 128 | none | — |
| ExpAssete | 8+16x | 8+8x | 8+2x | 24 | 26 | 24/24 | 0 | 64 | 24 ⚠ | — |
| ExpAssetz | 8+16x | 8+8x | 8+4x | 24 | 28 | 24/24 | 0 | 64 | 24 ⚠ | — |
| ExpAssetze | 4+4x | 4+16x | 4+2x | 12 | 22 | 12/12 | 0 | 32 | 12 | — |
| ExpAssetsemiz | 8+4x | 8+6x | 8+2x | 24 | 12 | 24/24 | 0 | 64 | none | — |
| RiskyAsset | 16+0x | 16+25x | 16+4x ⚠ | 48 | 29 | 48/48 | 0 | 80 | none | 50 |
| **total** | | | | **260** | **179** | **260/260** | **64** | | **150** | **74** |

439 subtests in the main banks, plus 150 QH and 74 EZ = **663**.

⚠ **The RiskyAsset with2A tier is partly test-first** (added 2026-08-18). The 16 subtests +
4 cross-tests exist. Of the 48 2A raws the toolkit needs, the 24 nosemiz ones (`_DC2A_`,
`_GI2A_`, `_DC2A_GI2A_` under `RiskyAsset/{DivideConquer,GridInterpLayer,
DivideConquerGridInterpLayer}`) landed the same day, so figs 33-40 should run — though only
`DC2A_nod1_noz` has been GPU-proven so far. The 24 SemiExo twins do not exist:
`RiskyAsset/RiskyAssetSemiExo/*` still holds only the `_DC1_`/`_GI1_`/`_DC1_GI1_` raws, so
every semiz DC / GI / DC+GI block (figs 41-48) errors until stage D lands. Base + lowmemory
pass throughout (the base raw is dimension-generic in `n_a1`).

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

| family | total raws | QH | EZ | core |
|---|---|---|---|---|
| baseline (FHorz, excl. asset families) | 406 | 224 | 40 | 142 |
| ExperienceAsset | 128 | 0 | 0 | 128 |
| ExperienceAssetu | 128 | 0 | 0 | 128 |
| ExperienceAssete | 192 | 128 | 0 | 64 |
| ExperienceAssetz | 96 | 32 | 0 | 64 |
| ExperienceAssetze | 84 | 52 | 0 | 32 |
| ExperienceAssetsemiz | 64 | 0 | 0 | 64 |
| RiskyAsset | 132 | 0 | 52 | 80 |

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
- ⚠ **Three QH banks are test-first against solvers in varying states.** Status as of
  2026-08-18, read off the working tree rather than the banners:
  - `CoreFHorzQHExpAssetTests` (48 subtests) — genuinely unsupported. `ExperienceAsset` has
    **0 QH raws** and `ValueFnIter_Case1_FHorz` has no QH branch under
    `vfoptions.experienceasset>=1`. Every subtest errors at its first `ValueFnIter` call.
    Same for `ExperienceAssetu` and `ExperienceAssetsemiz` (no mirrors written).
  - `CoreFHorzQHExpAsseteTests` (24 subtests) — **its banner is now stale.** It claims "NO
    quasi-hyperbolic support for experienceassete at all", but the working tree has 128 QH
    raws, four `QuasiHyperbolicExpAssete*` dispatchers plus four SemiExo ones, and live
    routing at `ValueFnIter_Case1_FHorz.m:401` / `:411`. The whole
    `ExperienceAssete/QuasiHyperbolic/` directory is untracked, so the banner was accurate
    when written and the solvers landed after it. Needs a run to find out where it stands.
  - `CoreFHorzQHExpAssetzTests` (24 subtests) — partial, and the banner is accurate: figs 1–8
    (noa1) and 13–16, 21–24 (semiz) error. `ValueFnIter_FHorz_QuasiHyperbolicExpAssetz.m`
    still raises `noa1 variant not yet implemented` at lines 281/289 and the `_e` version at
    136/144, and there is no `QuasiHyperbolicExpAssetzSemiExo` dispatcher at all.
  - `CoreFHorzQHExpAssetzeTests` (12 subtests) — the one that passes. No banner; GPU-validated.
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
- **RiskyAsset with2A is now covered on the test side only.** As of 2026-08-18 the bank has a
  full with2A1 tier (16 variants, figs 33–48, plus 4 degenerate-`a1_2` cross-tests), but the
  toolkit still has zero `*2A*` raws for the family, so this is a *toolkit* gap. Phase 2 is 48
  raws — `{DC2A, GI2A, DC2A_GI2A}` × `{nosemiz, SemiExo}` × 8 shock/decision combos — plus a 2A
  branch in the 6 DC/GI dispatchers, a `level1n=min(level1n,n_a1)` → `n_a1(1)` fix in each of
  them, `if n_a1>0` → `if prod(n_a1)>0` in `ValueFnIter_FHorz_RiskyAsset.m:174`, the UnKron
  level bump, and `ValueFnFromPolicy` 2A support. That 48 is also the entire difference between ExperienceAsset's 128 core raws and
  RiskyAsset's 80 — both cover the same 8 shock/decision combinations across
  {base, DC, GI, DC+GI} × {nosemiz, SemiExo}, but ExperienceAsset has a DC2A/GI2A/DC2A_GI2A
  variant in each of the 6 DC/GI directories (8 × 6 = 48). Base dirs are unaffected: brute
  force Krons the a1 dimensions, so no separate 2A raw is needed there.
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
