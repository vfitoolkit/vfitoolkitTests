# VFI Toolkit — test coverage

Coverage of the FHorz core test banks. The transition-path banks are covered separately, in
the FHorz transition-path section at the end. Tier is assigned from each subcode's directory
(`Noa1_subcodes` / `With2A1_subcodes` / else), not its filename — filenames misattribute
cross-tests. Format in the tier columns is `variants + cross-testsx`.

`ResidAsset` is the one bank with a single tier, so it has no tier subdirectories at all: its 16
withA1 subcodes sit at the top level of `..._subcodes/` and in `Semiz_subcodes/`. A tier rule
that keys off directory names must fall through to withA1 for it, not skip it.

Last updated: 2026-09-05

## CoreFHorzTests + the ExpAsset family + RiskyAsset + ResidAsset

| bank | noa1 | withA1 | with2A | var | x-tests | panel vs dist | V_Jplus1 | core raws | QH | EZ | AA |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CoreFHorzTests | — | 16+6x | 16+4x | 32 | 10 | 32/32 | 32 + 32 QH + 16 EZ | 142 | 42 | 24 | 16 |
| ExpAsset | 16+10x | 16+14x | 16+2x | 48 | 26 | 48/48 | 48 + 48 QH | 128 | 48 | — | — |
| ExpAssetU | 16+10x | 16+14x | 16+2x | 48 | 26 | 48/48 | 48 + 48 QH | 128 | 48 | — | — |
| ExpAssete | 8+16x | 8+8x | 8+2x | 24 | 26 | 24/24 | 24 + 24 QH | 64 | 24 | — | — |
| ExpAssetz | 8+16x | 8+8x | 8+4x | 24 | 28 | 24/24 | 24 + 24 QH | 64 | 24 | — | — |
| ExpAssetze | 4+4x | 4+16x | 4+2x | 12 | 22 | 12/12 | 12 + 12 QH | 32 | 12 | — | — |
| ExpAssetsemiz | 8+4x | 8+6x | 8+2x | 24 | 12 | 24/24 | 24 + 24 QH | 64 | 24 | — | — |
| RiskyAsset | 16+0x | 16+25x | 16+4x | 48 | 29 | 48/48 | 32 + 32 EZ | 128 | none | 50 | — |
| ResidAsset | n/a | 16+22x | — | 16 | 22 | 16/16 ‡ | 16 | 4 ‡‡ | none | — | — |
| **total** | | | | **276** | **201** | **276/276** | **520** | | **222** | **74** | **16** |

477 subtests in the main banks, plus 222 QH, 74 EZ and 16 AA = **789**.

**Every experience-asset family now has a QH mirror with a solver behind it, and every mirror is
complete at raw level.** ExpAssetU (48) and ExpAssetsemiz (24) were the last two showing `none`;
both are implemented and have run. ExpAssetz and ExpAssetze were the last two whose mirrors were
short of a full raw set — the missing nosemiz 1A tiers landed in `6031762b` (2026-09-01) and both
banks are green. All six families now sit at exactly 2× their core raw count.

‡ written but not runnable — `SimulateTimeSeries/` has no `residualasset` support at all.
‡‡ the bank is written test-first against 64 raws; the toolkit has 4 (this said 2 until
2026-09-01 — the raw-count table below had already been updated to 4 and this footnote had not).
See below.

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

## ResidAsset — new bank, written test-first (2026-08-25)

`CoreFHorzResidAssetTests` exists: **70 files, 16 variants, 22 cross-tests, 661 exact checks**,
figures 1–16. It was briefly the only test-first bank in the suite; `CoreFHorzAmbiguityTests`
(9 files) joined it on 2026-08-28 — and left again on 2026-09-01, when its solver wave landed and
the bank ran GPU-green first try. **ResidAsset is once more the only bank whose solver does not
yet exist** — 4 of a needed 64 raws.

| | |
|---|---|
| exact checks in the 16 variants | 528 |
| exact checks in the 22 cross-tests | 133 |
| of which `V_Jplus1` | 240 |
| of which the `lowmemory` ladder | 304 |
| `close to zero` convergence lines | 16 |
| anti-vacuity guards | 16 |

**Nothing here has been run, and nothing can be.** `CreateResidualAssetFnMatrix_Case1.m:11`
references `d_gridvals` before it is defined on line 14 (the parameter is named `d_grid`), and
both raws call it at every `jj`, so *every* residual-asset solve currently errors. That defect
has been live in the working tree for months; the family had no test bank, which is exactly why
nothing caught it. One-line fix: rename the parameter to `d_gridvals` and delete line 14 — the
callers already pass gridvals, so re-converting would double-convert.

### Why there is no noa1 tier

`rprimeFn` is a function of `(d, a1prime, a1, ...)`, so the residual asset is the residual of the
standard asset's *choice* and there has to be a standard asset for it to be a residual of. With
`isscalar(n_a)` the dispatcher sets `n_a1=0`, `N_a=prod(0)=0`, and the return-fn matrix
collapses. Structural, like `CoreFHorzTests` — not a missing raw. `with2A1` is deferred.

### What the bank assumes and the toolkit does not have

| | today | assumed |
|---|---|---|
| VFI raws | 2 (`_raw`, `_nod_raw`) | 64 |
| solver tiers | base only | base, DC1, GI1, DC1_GI1 |
| `z` | mandatory (explicit `error`) | optional |
| `e` | raw calls commented out | supported |
| `semiz` | nothing | supported |
| `lowmemory` | 0 and 1, no `else` | full ladder |
| SimPanel | none | supported |
| `ValueFnFromPolicy` | no `residualasset` branch | supported |
| StationaryDist | `z`/`e`; errors on `N_z==0` and on `gridinterplayer==1` | all combos + GI |

64 = 8 shock/decision combos × {base, DC1, GI1, DC1_GI1} × {nosemiz, SemiExo}, which is exact
parity with ExperienceAsset's own withA1 tier (64 of its 128).

### The cross-tests, and why this bank leans on them

With one solver tier and a one-step `lowmemory` ladder there is very little *internal* agreement
to check, so the oracles do the work here more than in any other bank:

| # | test | files | checks |
|---|---|---|---|
| 1 | z as e (and z/e halves) | 4 | 36 |
| 2 | semiz as z | 2 | 6 |
| 3 | single-point z vs genuine noz | 4 | 12 |
| 4 | **ResidAsset vs ExperienceAsset** | 4 | 48 |
| 5 | plain vs frozen residual asset | 4 | 16 |
| 6 | rprime-on-grid vs two-endogenous-state | 2 | 9 |
| 7 | degenerate semiz vs nosemiz | 2 | 6 |

Cross-test 4 is the load-bearing one. `rprimeFn` cannot see `r` and `aprimeFn(d2,a2)` can see
`a2`, so the two families do not nest in general — **but they coincide exactly when the
transition is memoryless**, i.e. `phi1*(1-d2)` on both sides. Both then map `d2` to the same
value by the same lower-point-plus-probability scheme, share the same ReturnFn (`r` and `a2`
occupy the same argument slot), and produce the same Policy channels. It runs at all four solver
tiers, so it is also the acceptance test for each tier of raws as they land. If it ever returns
~1e-16 rather than 0, diff `CreateResidualAssetFnMatrix_Case1` against
`CreateExperienceAssetFnMatrix_Case1` before suspecting the residual asset: they are separate
code computing the same weights.

### Conventions this bank pins down

- **`ReturnFn(d..., a1prime, a1, r, semiz, z, e, params...)`** and **`rprimeFn` = the same list
  minus `r`.** That absence is the family's defining feature. `semiz` and `e` must both be
  inputs, because `r` is a budget residual and both shift the budget. The
  `rprimeFnParamNames` split currently assumes `(l_d+l_a1+l_a1+l_z)` and needs
  `(l_d+l_a1+l_a1+l_semiz+l_z+l_e)` — in the VFI *and* in `StationaryDist_FHorz_ResidAsset`,
  which keeps its own copy (with a comment claiming the leading inputs are `(d2,a2)`).
- **The residual asset has no decision of its own.** No `d2` analog: `n_d` is `0` / `d1` /
  `dsemiz` / `[d1,dsemiz]`, semiz decision last. `nod1_noz_noe_nosemiz` has no `d` at all.
- **`n_r>=2` is required** — `r_griddiff=r_grid(2:end)-r_grid(1:end-1)` is empty at `n_r=1`. The
  single-point-degenerate trick that works for `z` does not work for `r`, which is why
  cross-test 5 freezes `rprimeFn` rather than shrinking the grid.
- **`rprime` depends on `a1prime`**, so the arrays scale with `n_a(1)` **squared** — ExpAsset's
  `aprimeFn` arrays scale with `n_d2*n_a2`. `n_a_big` is 151, not the 1001 the ExpAsset bank
  uses, and `d1_z_e_semiz` is the expected OOM point.
- **The anti-vacuity guard is new here and worth copying.** `habit=0.2` is what makes the
  residual asset payoff-relevant; without it `V` would not vary in `r` and all 661 checks would
  pass on a degenerate model. Every variant prints `Residual asset matters, this should be well
  above zero`, filtered to finite entries because `c > habit*clag` creates real `-Inf` regions.
  The RiskyAsset bank's vacuous `d2recon` calibration is the cautionary case this guards against.


### Column meanings

- **noa1 / withA1 / with2A** — the three asset tiers. `noa1` does not apply to the baseline,
  nor to ResidAsset (`n/a` in both cases, rather than `0`).
  Variant counts are uniform within a bank by construction: 16 per tier where three shock
  axes are free (`{nod,d} × {noz,z} × {noe,e}`, doubled for semiz), 8 where the family's
  driving shock fixes one axis, 4 for ExpAssetze where both `z` and `e` are mandatory.
- **panel vs dist** — variants that simulate a panel *and* compare it against the agent
  distribution (`AgeConditionalStats` / `AllStats`). Cross-tests are excluded: they assert
  bit-exact V/Policy/Dist agreement, which a Monte-Carlo panel cannot improve on.
  These are printed side-by-side rows, eyeballed, not thresholded assertions.
- **V_Jplus1** — subcodes that exercise the terminal-condition branch.
- **core raws** — solver raws in the family, excluding exotic-preference variants
  (which the QH/EZ/AA mirrors cover). This is the surface the main bank dispatches into.
- **QH / EZ / AA** — subtests in the bank's `withQuasiHyperbolicDiscounting` /
  `withEpsteinZinPreferences` / `withAmbiguityAversion` mirror. `—` means no mirror. (No `⚠`
  remains anywhere: the last test-first QH bank, ExpAsset, got its solver on 2026-08-22, and the
  AA mirror got its solver on 2026-09-01 and its with2A tier on 2026-09-02 — see the
  AmbiguityAversion section. AA is 16 rather than a QH-style 42 because ambiguity has no noz+noe
  variant (deliberately an error) and no semiz tier (deferred by design, not test-first); its
  with2A tier is 6+2x rather than 8+2x for the same noz+noe reason.)

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

Panel simulation is absent everywhere — `SimPanelValues` appears in 0 of 155 QH, 0 of 220
EZ and 0 of 17 AA subcodes. Dist-derived output is not absent, but it is confined to two banks (TPath banks
excluded throughout):

| mirror / bank | files | asserted | displayed |
|---|---|---|---|
| QH — CoreFHorzTests | 43 | 120 | 512 |
| QH — ExpAsset, ExpAssete, ExpAssetz, ExpAssetze | 112 | **0** | **0** |
| EZ — CoreFHorzTests | 73 | 120 | 384 |
| EZ — RiskyAsset | 147 | 100 | 240 |
| AA — CoreFHorzTests | 17 | 12 | 96 |

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
168 total, 104 QH (40 nosemiz + 64 semiz), 64 core. ExperienceAsset re-counted 2026-08-25 after
the QH work landed (`5cada730`): 400 total, 256 QH, 128 core — it had 0 QH raws when this table
was first written. ResidualAsset added 2026-08-25: 2 raws, and they do not currently run.

Re-counted 2026-08-31. ExperienceAssetu went 128 → 384 (`71856d75` nosemiz, `2aaecbc6` SemiExo);
ExperienceAssetsemiz 64 → 192 (`d5cef0af`); ResidualAsset 2 → 4.

Re-counted 2026-09-01: AmbiguityAversion's 6 old-style raws were deleted and replaced by 24
modern ones ({plain, DC1, GI1, DC1_GI1} × {nod,d} × {z, noz_e, z_e}), all mechanical transforms
of the modern exponential raws; the AA column (and the baseline total, 406 → 430) is new. The
old table's baseline total never counted the ambiguity (or GulPesendorfer) raws.

Re-counted again 2026-09-01 after `6031762b`: ExperienceAssetz 168 → 192 (QH 104 → 128) and
ExperienceAssetze 84 → 96 (QH 52 → 64), closing the last two raw-level QH gaps. The
ExperienceAsset row also drops 400 → 384: that figure was counting all `.m` files in the family,
including its 16 dispatchers, while every other row counted `*_raw.m` only. **This table is
`*_raw.m` only** — dispatchers are not raws and are counted in the prose, not here.

Re-counted 2026-09-02: GulPesendorfer's 7 old-style plain raws (never counted, and never
exercised — the run-1 diary caught a pre-split-Policy2-vs-UnKron mismatch and a 3-D kron in
them) grew to 32 modern ones ({plain, DC1, GI1, DC1_GI1} × {nod,d} × {noz_noe, z, noz_e, z_e};
unlike AA, GP works with no shocks, so 4 shock combos not 3), all GPU-green against the new
withGulPesendorferPrefs bank (266 zero-checks, figs 1-8). The GP column (and the baseline
total, 430 → 462) is new. No 2A/semiz GP raws yet (tier dispatchers error on those).

Re-counted 2026-09-03: GulPesendorfer 32 → 88 (baseline total 462 → 518): the semiz tier
(32 raws: {plain, DC1, GI1, DC1_GI1} × {d1,nod1} × {z,noz} × {e,noe}, per-d2 most-tempting
collection) and the 2A tiers (24 raws: {DC2A, GI2A, DC2A_GI2A} × 8 combos) landed together,
plus 4 SemiExo dispatchers, 2A branches in the GP DC/GI/DC_GI dispatchers, and
ValueFnFromPolicy GP semiz + GI2A support. GPU-green against the extended bank (746
zero-checks, figs 1-24, incl. a just-a-markov semiz-vs-markov equivalence at baseline
temptation that is exact-zero at all four tiers). GP is now the only exotic preference at
full {plain,DC,GI,DC+GI} × {1A,2A} × {nosemiz,semiz(1A)} coverage. Note: GP holds the
return matrix and its temptation twin simultaneously (~2x core memory), so the biggest
bank cases run their moments blocks on reduced grids.

Re-counted 2026-09-03 (later the same day): GulPesendorfer 88 → 112 (baseline total 518 → 542):
the semiz × 2A cross tier landed (24 raws: {DC2A, GI2A, DC2A_GI2A} × {d1,nod1} × {z,noz} ×
{e,noe} under SemiExo/, per-d2 most-tempting collection with the joint (d1,a1prime,a2prime)
max, interpolation on a1prime only), plus 2A branches in the GP SemiExo tier dispatchers and
2A support in ValueFnFromPolicy GP SemiExo. GP is the first family with the full
{plain,DC,GI,DC+GI} × {1A,2A} × {nosemiz,semiz} cube. Bank extension: 8 semiz-with2A cases +
2 cross tests (figs 25-32). GPU-GREEN 2026-09-07 (266 checks all at/below the ULP floor;
just-a-markov pins the semiz-2A machinery to the validated nosemiz-2A family, 12 exact
zeros + 4 at 8e-15). The one run-1 find: the GP SemiExo_DC dispatcher had wrongly bumped
its 2A UnKron level — DC raws return the joint aprime kron index, so DC keeps the 1A
UnKron calls (only GI/DC_GI bump); fixed to match the core dispatcher.

Re-counted 2026-09-14: GulPesendorfer 112 → 152, with the 40 new raws sitting in the
ExperienceAsset family row (384 → 424): the GulPesendorferExpAsset family ({plain incl
noa1, DC1, GI1, DC1_GI1} × the 1A1 nosemiz combos) landed with its own sub-bank
(CoreFHorzExpAssetTests/withGulPesendorferPrefs/, 508 checks, GPU-green 2026-09-14; 24
checks sit at 3-4 ULP of the poor-corner V~1.1e7, accepted as the scale-adjusted ULP
floor). Run-1 find: experience assets create states where EVERY choice is infeasible
(earnings scale with a2, a2 grid starts at 0) — GP's V=max(-Inf)-max(-Inf)=NaN then
corrupted continuations via the EV NaN-scrub; fixed by guarding every MostTempting
(-Inf -> 0, so V=-Inf exactly as the standard solver). The GP-core 112 raws carry the
same latent pattern but their models cannot trigger it (earnings never vanish); a
guarded port there is a recorded follow-up. With2A1/SemiExo GP-ExpAsset tiers and the
u/e/z/ze/semiz siblings remain gaps (dispatchers error). Unblocks KLM2021.

| family | total raws | QH | EZ | AA | GP | core |
|---|---|---|---|---|---|---|
| baseline (FHorz, excl. asset families) | 542 | 224 | 40 | 24 | 112 | 142 |
| ExperienceAsset | 424 | 256 | 0 | 0 | 40 | 128 |
| ExperienceAssetu | 384 | 256 | 0 | 0 | 0 | 128 |
| ExperienceAssete | 192 | 128 | 0 | 0 | 0 | 64 |
| ExperienceAssetz | 192 | 128 | 0 | 0 | 0 | 64 |
| ExperienceAssetze | 96 | 64 | 0 | 0 | 0 | 32 |
| ExperienceAssetsemiz | 192 | 128 | 0 | 0 | 0 | 64 |
| RiskyAsset | 208 | 0 | 80 | 0 | 0 | 128 |
| ResidualAsset | 4 | 0 | 0 | 0 | 0 | 4 |

**Every experience-asset family is now at exactly 2× its core count** — one Naive and one
Sophisticated raw per exponential raw, which is the expected ratio. The two families that were
short (ExperienceAssetz at 104 of 128, ExperienceAssetze at 52 of 64, both missing only nosemiz
1A tiers) were closed in `6031762b`; see below.

Note the toolkit is inconsistent about where exotic-preference raws live: ExpAssete/z/ze
keep them in a `QuasiHyperbolic/` subdir of their own family, while the *baseline* QH raws
sit under `ExoticPrefs/QuasiHyperbolic/`.

### QH ExperienceAssetu closed (2026-08-31)

`ExperienceAssetu` was the last family with **zero** QH raws. It now has 256 — the nosemiz half in
`71856d75`, the SemiExo half in `2aaecbc6` — plus 8 dispatchers and 4 `ValueFnFromPolicy` subfns.
The 48-config bank runs end to end: 1933 checks, no value-function failure.

Each of the 16 tiers was generated from its own exponential source and held to four instruments:
byte-exact reproduction of the committed QH ExpAsset references, an inverse check reducing back to
the ExpAssetu source, the structural gate, and a semantic checker mutation-tested against seeded
faults. **Measured across those tiers the structural gate caught 0–2 of every 15–34 seeded faults,
and the inverse check roughly half.** Treat gate-clean as a syntax precondition, never as evidence.

The dispatchers get a different instrument: a call-graph check calibrated on the committed QH
ExpAsset tree (277 calls, 256 raws, 0 findings) that verifies every callee exists, argument counts
match the callee's signature, Naive/Sophisticated branches reach the matching raw, and every raw is
reached exactly once. It does not catch a wrong *value* passed through a correct signature: the
SemiExo dispatchers derived `aprimeFnParamNames` with the nosemiz spelling, leaving `l_u` computed
but unused so that `u` was taken as a parameter name. Every call site was well-formed. What caught
it was an **assigned-but-never-used** check over the wired dispatchers, now part of that pass.

Two defects found in shared code along the way, both now fixed:

- **`aprimeProbs` accumulated its `skipinterp` zeroing across the `d3` loop** — built once per age
  outside the loop, mutated in place inside it, so a cell zeroed at one `d3` stayed zero at the
  next. 1116 sites over 192 SemiExo raws in both the ExpAsset and ExpAssetu families, the
  GPU-validated reference family included. Latent in ExpAsset (whose base tier already rebuilt per
  iteration) and live in ExpAssetu (whose base tier did not).
- **The GI `ValueFnFromPolicy` files never cleared NaN** after combining 2×2 corner weights with
  continuation values. `w_a1_upper` is exactly zero whenever the policy sits on a lower grid point,
  so `0*(-Inf)` at an infeasible corner produced NaN in `V`. 58 guards across 22 files, in the
  idiom the RiskyAsset and SemiExo GI files already used.

Both fixes are in `d7dfcd7c` and were validated by four bank runs: QH ExpAsset identical to its
pre-fix baseline over all 1933 checks, ExpAsset 514 exact-zero checks all zero, QH ExpAssetU and
ExpAssetU value-exact throughout.

**Every remaining non-zero exact check in the FHorz suite is a Policy tie** (plus, since
2026-09-01, the AA bank's 18 ULP-floor V lines — see the AmbiguityAversion section), not a wrong
value —
demonstrated rather than assumed: one was measured at 30 differing entries of 1 260 480, confined
to a single policy channel, every difference exactly ±1, with the achieved value bitwise identical
on both sides.

### QH ExperienceAssetz/ze 1A tiers closed (2026-09-01)

The last raw-level QH gap. Both families had full SemiExo mirrors but only the base method on the
nosemiz side, so `DC1`, `GI1` and `DC1_GI1` were missing: ExpAssetz 24 raws short, ExpAssetze 12.
`6031762b` adds all 36, and both families now sit at 2× core.

The leaf sets differ because the families differ. ExpAssetz has `z` structural in
`aprimeFn(d2,a2,z,…)`, so there are no `noz` variants and the leaves are `{d1,nod1}×{e,noe}` — 24
raws. ExpAssetze has **exactly one `e`**, serving both the structural role and the i.i.d. one, so
the leaves are `{d1,nod1}` alone — 12 raws.

Wiring took two shapes, and which one applies is forced by where the 2A routing already lives:

- **ExpAssetz** already had `{DC, GI, DC_GI}` sub-dispatchers carrying 2A. They gained 1A routing
  with their 2A paths left byte-identical (5250/5250/5290 bytes unchanged), and the base
  dispatcher's not-yet-implemented guards came out.
- **ExpAssetze** had no sub-dispatchers; its 2A routing sat directly in the family dispatcher.
  Three sub-dispatchers were created and the 2A block moved into them verbatim — 53 substantive
  lines, +13/−82 on the family dispatcher. This was the riskiest edit in the batch, because it
  relocates GPU-validated code. One ordering decision was deliberate: the 2A block stays ahead of
  the 1A `level1n` setup, since `min(vfoptions.level1n, n_a1)` against a vector would change what
  2A's own `min(…, n_a1DC)` produces.

The six `withA1` bank subcodes were extended from base-method-only to the full tier ladder
(`3c64211`, +244 checks in ExpAssetz, +144 in ExpAssetze). Both banks then ran green:

| bank | configs | checks | exactly 0 | at ULP floor | errors |
|---|---|---|---|---|---|
| QH ExpAssetz | 24/24 | 1512 | 1248 | 264 | 0 |
| QH ExpAssetze | 12/12 | 920 | 772 | 148 | 0 |

Three discriminators separate "at the floor" from "wrong", and all three are clean in both runs:
every `lowmemory` check exactly zero, every β₀=1 exponential cross-check exactly zero, and **zero
Policy differences anywhere** — not one argmax tie to explain away. All 412 non-zeros sit on
`ValueFnFromPolicy` (which re-evaluates `V` from `Policy`, so a different summation order) or on
the `DC1_GI1 vs GI1` / `DC2A_GI2A vs GI2A` method comparisons, and every distinct magnitude in
both diaries is dyadic (`2^-33` … `2^-28`), which is what representation-level difference looks
like and is not what a real numerical fault looks like. That covers the relocated 2A too: it
contributed 234 checks in ze with every method and β₀ comparison bitwise zero.

### V_Jplus1 project, bank 1: CoreFHorzExpAssetTests (2026-09-02)

First bank of the twelve-bank effort to put runtime coverage behind the `V_Jplus1` branches of the
experience-asset families. All 48 subcodes gained a `V_Jplus1` section (`d246b10`, `164c0a3`),
**700 checks across 42 completed variants, all exactly zero** after one real bug was fixed.

There was no toolkit gate to open: the dispatchers pass `vfoptions` straight through and the raws
read `vfoptions.V_Jplus1` directly, so every branch was already reachable. This is a test-writing
project; the toolkit work is only whatever the tests break.

**The donor is `CoreFHorzTests_subcodes`**, which already had a complete 32-subcode precedent over
the identical d×z×e×semiz structure at 1-asset and 2-asset tiers. Generated blocks were validated
by exact reproduction — 10/10 variants matched the donor's code lines byte-for-byte — so only the
3-asset `with2A1` tier had no precedent, and it differs by one colon. Three rules parameterise the
block, with `a` = asset dims and `s` = shock dims:

- V slice takes `a+s` colons, Policy `a+s+1`
- the lowmemory ladder runs `1..s`, and there is no ladder at all when `s=0`
- the age-dependent-shocks block appears iff z or e; semiz adds nothing to it

`pi_z_J` trims to `1:Njs` but `pi_e_J` trims to `1:jstar` — a pi_e column is the distribution of e
*realised* in period j, so the shorter model needs the column for the V_Jplus1 period too.

**The noa1 tier takes a different shape.** With no standard asset a1 there is nothing for DC or the
grid interp layer to operate on, so those 16 subcodes define only `vfoptions1`; giving them the
four-method block killed the first run at `Unrecognized function or variable 'vfoptions2'`. They now
use the `CoreFHorzRiskyAsset_*_noa1.m` shape — one method looped over `jstar=[round(3*N_j/4),N_j]`.
Generalises: **check that every `vfoptionsN` a generated block references is defined in the target
file**, not merely that the block's own shape is right.

#### The bug it caught

Twelve ExperienceAsset raws built `aprimeProbs` in their `V_Jplus1` branch as
`repmat(a2primeProbs,N_a1,1,1,N_z)` where their own in-loop code uses `repmat(...,N_a1,1,N_z)`.
In this family `a2primeProbs` is 2-D `[N_d2,N_a2]`, so the extra factor inserts an interior
singleton — `[D,A2,1,Z]` instead of `[D,A2,Z]` — which **broadcasts** against `Vlower [D,A2,Z]` to
give `[D,A2,Z,Z]`. The later `sum(...,3)` restores the correct shape, so it never errored; it
returned wrong values wherever `skipinterp` was non-empty, since that zeroing is what makes the
probabilities z-dependent. The line's own comment still documented the correct shape.

It surfaced as 6 of 642 checks non-zero — `8.631e-04`, `1.302e-03`, and a Policy difference of
exactly `1.000e+00` — in `d1_z_noe_nosemiz` and its `with2A1` sibling, both at the base method,
`jstar=15`. Every other non-zero in that diary was at the ULP floor. Fixed in `070afc87`.

**Method note, because the naive check gives a flood of false positives.** The 4-factor form is
*correct* in the z-bearing families, where `a2primeProbs` is 3-D, and 252 sites use it
legitimately. The valid test is per-file: **diff each file's `V_Jplus1` branch against its own
in-loop line**, since the in-loop code is exercised by every existing test. That sweep over 1082
raws returned exactly 12 — the `d1+z+noe` leaves of {base, DC1, GI1, DC1_GI1} × {exponential, QH-N,
QH-S}. This is the opposite of the "fewer `repelem` factors is legal" rule recorded elsewhere:
*fewer* is harmless because trailing singletons drop, *more* is not.

8 of the 12 are QH raws whose bank has no `V_Jplus1` tests yet, so those fixes ride on this
evidence rather than their own — which makes QH ExpAsset the natural bank 2.

#### Bank 2: CoreFHorzQHExpAssetTests (2026-09-02)

All 48 subcodes gained a V_Jplus1 section, **+2560 checks** (bank 1933 → ~4493), committed
`f684de0`. The donor is the **QH mirror of CoreFHorzTests**, and the generated code reproduces it
exactly across all 12 (a,s) combinations it covers. First run: **2056 checks reached over 42 of 48
subcodes, all exactly zero** — including the eight QH raws fixed in `070afc87`, which until then
rested on the exponential bank's evidence alone.

The quasi-hyperbolic specifics, all inherited from that donor:

- **`Valt` is what gets fed back as `vfoptions.V_Jplus1`, never `V`.** The continuation value in
  the Bellman equation is the standard-discounted one — `V_std` when Naive, `Vunderbar` when
  Sophisticated — which is the *third* output. `V` is Vtilde/Vhat, the QH-discounted object.
- Eight blocks per subcode: {Naive, Sophisticated} × {base, DC, GI, DC+GI}. Naive returns four
  outputs and gets four checks at each jstar (V, Policy, Valt, **Policyalt**); Sophisticated
  returns three. Both get three per lowmemory rung. So **28+24s** checks per subcode, or
  **14+12s** for noa1.
- The section sits **before** `%% Versus exponential discounting`, which sets `beta0=1` and would
  collapse QH onto exponential, disabling exactly what these tests exercise.
- `vfoptions1..4` are redefined in the Sophisticated section, so every block sets
  `quasi_hyperbolic` explicitly rather than inheriting whatever they currently hold.

**Deliberate gap:** the QH donor has no age-dependent-shocks block in any of its 32 subcodes, so
bank 2 has none either. QH + age-dependent `pi_z`/`pi_e` + `V_Jplus1` is untested across every QH
bank, and closing it would need a new block with no validated precedent.

The three non-zero lines in this bank's diary are the long-standing DC2A argmax ties — Naive
Policy 2, Naive Policyalt 1, Sophisticated Policy 1, unchanged from the 2026-08-22 run, with V and
Valt agreeing to `9.313e-10` and `1.863e-09`. Note this diary no longer carries the per-tie
`[diag]` lines that the 2026-08-22 run had; the V/Valt magnitudes beside the Policy difference are
what identify them as ties.

**Two caveats on bank 1's reach.** The run OOMs at config 40 (`d1_z_e_nosemiz_with2A1`), a
pre-existing ceiling; because the call is not commented out, the OOM *aborts the script*, so
configs 41–48 only run when the post-40 part is invoked separately. And one non-zero remains in the
diary — `Divide-and-conquer (DC2A)` Policy `1.000e+00` with V agreeing to `1.863e-09` — which is an
argmax tie, not a defect. The contrast is the argument: the real bug showed `1.302e-03` alongside
the same Policy difference, four orders of magnitude away.

#### Bank 3: CoreFHorzExpAssetUTests (2026-09-03) — five defects

All 48 subcodes gained a V_Jplus1 section (**+896 checks**, `900ad5b`), generated from the same
CoreFHorzTests donor as bank 1. `u` needs no special handling: it is iid, `vfoptions.pi_u` carries
no age dimension, and it is integrated out inside the a2prime transition, so it adds no dimension
to V or Policy.

**This bank took five runs, each stopped by a different toolkit defect — all five in V_Jplus1
branches, none reachable by any bank before this project.** Fixed in `ca9cdc37` (13 files):

| # | defect | sites |
|---|---|---|
| 1 | `DiscountedEVinterp` never expanded over `N_d1` in the `lowmemory==0` sub-branch | 1 |
| 2 | `n_a1,n_a2` passed twice — 19 args to a 17-parameter builder | 1 |
| 3 | `_Disc_noz` (13 params) called with the 15-argument `_Disc` list | 8 |
| 4 | reshape target contradicting its own builder's `n_e`/`special_n_e` | 17 |
| 5 | `zind` used in a `lowmemory==3` block that never defines it | 2 |

Defects 3 and 5 are **wave-1 ExpAssetu defects fixed in August and never swept into the sibling
families** — the recurring lesson that V_Jplus1 drift replicates across siblings.

Defect 4 is the subtlest and gives the most reusable rule: **a reshape target must carry `N_e`
exactly when the nearest preceding `CreateReturnFnMatrix_*` call is passed `n_e` rather than
`special_n_e`.** That is a purely local check, needing no comparison to the in-loop code. One of
the two wrong forms was `[...,N_bothz,1]`, where the trailing `,1` is a `reshape` no-op — so it
silently dropped the e dimension instead of erroring.

Each fix was generalised by a whole-toolkit sweep, and each now returns zero: 390 raws for repelem
counts, 45 854 calls for argument arity, 614 indexed reshapes against their builders, 7010
`for z_c` blocks for an undefined `zind`.

**Three of those sweeps first indicted GPU-validated reference code, and every time the checker
was wrong** — an unanchored regex matching `entireEV_d2=repelem(EV_d2,...)` as a substring, a
comma counter that ignored `[...]` brackets (20 380 false hits of 45 854), and a `for`/`end`
balance check whose region ran to end-of-file and counted the function's own terminating `end`.
When a sweep accuses the reference implementation, suspect the sweep first.

**Validation standing:** the bank exercises defects 1, 2 and 5. The 8 sites in defect 3 sit in a
conditionally-dead DC arm (`maxgap(ii)==0`), and the 15 ExpAssetz sites in defect 4 are in a
family whose V_Jplus1 bank does not exist yet — so 23 of the 29 fixed sites are unvalidated by any
run, the same standing as the `loweredge` fixes.

Final run: 42 of 48 subcodes, **700 V_Jplus1 checks, 694 exactly zero and 3 at the ULP floor**. The
3 non-zero are Policy checks in `d1_z_noe_semiz` against V agreeing to `2.220e-16`; that subcode's
*own* pre-existing Divide-and-conquer and lowmemory checks show the identical difference of 2, so
the tie predates these tests. Two OOMs in the with2A1 tier, the known ceiling.

The bank's two zero-shock noa1 subcodes also lost their `lowmemory=1` checks: with no shocks there
is nothing to loop over, all **192** zero-shock FHorz raws never mention `lowmemory`, and the check
compared a solve to itself. Eight such vacuous checks remain in other banks (AmbRiskyAsset 4,
GPFHorz 2, and 2 more), flagged to their owners.

#### Bank 4: CoreFHorzQHExpAssetUTests (2026-09-03) — clean first run

All 48 subcodes gained a V_Jplus1 section (**+2560 checks**, `ee6bef2`), from bank 2's QH template
plus bank 3's ExpAssetU facts. **First run green with no toolkit changes at all**: 39 of 48
subcodes, 1876 checks, 1857 exactly zero and 12 at the ULP floor.

The 7 above the floor are Policy/Policyalt differences of 1 or 2 in `d1_z_noe_semiz`, each paired
with a V check at `1.110e-16` or `2.220e-16`. That subcode's *own* pre-existing lowmemory and DC1
checks show the identical difference of 2 against V at `7.105e-15`, so the tie predates these
tests — and it is the same model that tied in the exponential ExpAssetU bank. One OOM, an
allocation failure on `entireRHS_tilde` in the QH twin.

**Why this bank was clean when bank 3 found five defects, and what it implies.** Four of bank 3's
five defects were in *exponential* ExpAssetu SemiExo raws, so the obvious expectation was that
their QH mirrors carried the same drift. Re-running all three class sweeps over the 256 QH
ExpAssetu raws (128 SemiExo) gave zero hits before the run, and the run confirmed it. The reason
is provenance: those QH raws were **generated** in August from the committed QH ExpAsset references
plus the `u` delta, whereas the exponential raws are older hand-maintained code. **Generated and
hand-maintained code carry different defect distributions — a defect found in one should not be
assumed to cross to the other, in either direction.**

#### Bank 5: CoreFHorzExpAsseteTests (2026-09-04) — clean, and a naming fix

All 24 subcodes gained a V_Jplus1 section (**+544 checks**, `642a14a`). **First run green with no
toolkit changes**: all 24 subcodes ran, 494 checks reached, every one exactly zero, and **not one
check anywhere in the diary above the ULP floor** — no ties either, which no other bank has
managed. One OOM, an `arrayfun` allocation failure in the largest configuration
(`z_e_with2A1_semiz`), which is why 494 of the 544 were reached.

**This family has 24 variants, not 48**, because `e` is structural in `aprimeFn(d2,a2,e)`: there is
no `noe` variant, and the shock count runs 1–3 rather than 0–3.

**A naming inconsistency fixed first (`0b23d38`).** ExpAssete was the only bank with a
semiz/nosemiz distinction that left the nosemiz case unmarked, writing
`CoreFHorzExpAssete_d1_z_e.m` where every other bank writes `..._d1_z_e_nosemiz.m`.
(`CoreFHorzExpAssetsemizTests` is also unmarked but correctly so — that family is always semiz.)
The fix renamed **36 files** — 12 subcodes, 12 QH-mirror subcodes and 12 ReturnFns — and rewrote
108 references across 56 files, including four `CoreFHorzExpAssetz` cross-tests that call ExpAssete
ReturnFns, so **the rename crosses bank boundaries**. The reference rewrite needs a negative
lookahead: `CoreFHorzExpAssete_d1_noz_e` is a prefix of `..._noa1`, `..._with2A1` and `..._semiz`.

Two benefits beyond consistency: the ExpAsset generator now works on this family **unchanged**, so
the two banks share one generator; and the run validated the rename by execution, with zero
unresolved function names across all 24 subcodes.

**This bank weakens a heuristic worth not over-trusting.** Banks 2 and 4 (generated QH code) were
clean and banks 1 and 3 (hand-maintained exponential code) carried six defects, which suggested
defects track provenance. Bank 5 is hand-maintained exponential code and is completely clean. The
provenance split is real for the QH mirrors but is not a reliable predictor for the exponential
families.

#### Bank 6: CoreFHorzQHExpAsseteTests (2026-09-04) — the first complete bank

All 24 subcodes gained a V_Jplus1 section (**+1520 checks**). **The cleanest run in the project**:
all 24 subcodes ran, **all 1520 checks were reached and every one is exactly zero**, nothing
anywhere in the diary above the ULP floor, and **no errors at all** — not even the with2A1 OOM
that stopped every previous bank short.

Bank 2's QH generator applied **unchanged**, which the `_nosemiz` rename made possible: without it
this family would have needed a token override. The two QH banks now share one generator, as the
two exponential banks share the other.

One method difference worth keeping: QH banks compute no agent distribution, so there is no
`jequaloneDist` to read asset dimensions from. Take them from the setup instead —
`n_a_justexpasset=13` (a=1), `n_a=[101,13]` (a=2), `n_a_2A1=[51,2,13]` (a=3) — after confirming the
QH main script passes the same variables as its exponential twin.

#### Bank 7: CoreFHorzExpAssetzTests (2026-09-04) — validates 3 of the 15 dormant fixes

All 24 subcodes gained a V_Jplus1 section (**+544 checks**). Green first run, no toolkit changes:
all 24 subcodes ran, 494 checks reached, every one exactly zero, nothing above the ULP floor. One
OOM in the largest configuration (`d1_z_e_semiz_with2A1`).

**24 variants, not 48**, because `z` is structural in `aprimeFn(d2,a2,z)`: no `noz` variant, shock
count 1–3. Bank 1's generator applied unchanged.

**What it validated, and what it did not.** Of the 15 reshape sites fixed in `ca9cdc37` that no run
had ever exercised, this bank reaches **3** — those in `ExpAssetzSemiExo_GI1_e_raw` (1) and
`ExpAssetzSemiExo_DC1_GI1_e_raw` (2). Both `_z_e_semiz_withA1` subcodes completed their GI and
DC+GI blocks with **34 of 34 checks exactly zero**, which is direct evidence for those three. **The
other 12 are in QH ExpAssetz raws and remain unvalidated until the QH bank runs** — a green bank 7
does not clear them.

#### Bank 8: CoreFHorzQHExpAssetzTests (2026-09-04) — two more defects, in generated code

Written with bank 2's QH generator unchanged (+1520 checks, 8/8 exact reproduction). The first run
reached 16 of 24 subcodes, 903 checks with 900 exactly zero, and surfaced **two defects — the
seventh and eighth of the project, and the first found in generated QH code.** After the fix
(`01e115aa`) the bank **completes entirely**: all 24 subcodes, all 1520 checks reached, every one
exactly zero, no errors and no OOM. **It also retires the last 12 of the 15 `ca9cdc37` reshape
sites**, verified by checking that the four `_z_e_semiz_with*A1` subcodes completed their GI and
DC+GI blocks under both preferences with 50 of 50 checks zero each — so with bank 7's 3, that
backlog is now empty.

Both are the same class: **the `e`-stride term in `maxindexfull` must match whether `e` is looped
at that lowmemory level**, and both are in the *Sophisticated* QH ExpAssetzSemiExo GI raws — the
`Vunderbar` machinery that only the Sophisticated dual exercises (the Naive raws have no
`V_ford3_under` at all and are structurally immune).

| site | lm | defect | symptom |
|---|---|---|---|
| `...S_GI1_e_raw:263` | 0 | `(0:1:(1)-1)` where `(0:1:(N_e)-1)` belongs | **silent**: V 5.083e-01, Valt 1.943e+00, Policy 44 |
| `...S_DC1_GI1_e_raw:385` | 0 | same | latent in that run |
| `...S_DC1_GI1_e_raw:469` | 1 | an extra `N_e` stride where e *is* looped | `Index exceeds matrix dimension` |

`(0:1:(1)-1)` evaluates to `0`, so at `lowmemory==0` — where e is not looped — **every e slice read
from e=1's block**, which is why it returned plausible wrong numbers rather than erroring. All
three now match their own in-loop counterparts byte-for-byte.

**Sweep discipline, again.** The first pass reported **297** mismatches including hits in
GPU-validated in-loop code, because one of its three rules ("e not looped but stride lacks `N_e`")
is unsound: DC level-2 code indexes through `curraindex` and does not follow it. Restricted to the
two rules that are wrong *on their face* — a literal `(1)` where a dimension size belongs, and an
`N_e` stride in an e-looped branch — it gives **7 of 3014 sites**, of which **4 are harmless**:
`(0:1:(1)-1)` in branches where e *is* looped is a no-op equivalent to omitting the term, and
their in-loop twins carry identical text. That leaves the 3 above.

#### Bank 9: CoreFHorzExpAssetzeTests (2026-09-04) — clean

All 12 subcodes gained a V_Jplus1 section (**+312 checks**, `d4de909`). Green first run, no toolkit
changes: all 12 subcodes ran, 270 checks reached, every one exactly zero, nothing above the ULP
floor. One OOM in the largest configuration.

**12 variants — the smallest bank in the project** — because both `z` and `e` are structural in
`aprimeFn(d2,a2,z,e)`: no `noz` or `noe` leaves, so the grid is only
`{d1,nod1} × {semiz,nosemiz} × 3 tiers` and the shock count is always 2 or 3.

#### Bank 10: CoreFHorzQHExpAssetzeTests (2026-09-04) — clean

All 12 subcodes gained a V_Jplus1 section. Green first run, no toolkit changes: **880 V_Jplus1
checks, every one exactly zero**, largest residual anywhere in the bank `5.329e-15`, no OOM. (It was
recorded in the summary table at the time — `a58bebe` — but never written up here.)

#### Bank 11: CoreFHorzExpAssetsemizTests (2026-09-07) — the largest defect of the project

All 24 subcodes gained a V_Jplus1 section (**+528 checks**). Every V_Jplus1 check that ran was
exactly zero, first time and after the fix below. But the bank also reported
`ValueFnFromPolicy, this should be zero: Inf` in **all 8 noa1 subcodes** — not a V_Jplus1 failure,
and worth chasing because it was confined to one tier of one family while the same family's withA1
tier read `1.776e-15`.

**The cause was not in ValueFnFromPolicy — it was in the solver, and it is toolkit-wide.** Every
asset-interpolation site combines neighbouring nodes with weights summing to one. When a weight is
exactly `0` and its node is `-Inf`, `0*(-Inf)` is `NaN`. Exact zero weights are routine: a2prime
landing on a grid point, off either end of the grid, and the `skipinterp` guard, which zeroes the
weight *precisely when the two nodes are equal* — so the guard written to handle the degenerate case
manufactured a NaN in exactly that case. One age later `EV(isnan(EV))=0` turned the NaN into a
continuation value of **zero**, pricing an infeasible dead end at 0. With CRRA utility negative
throughout, zero beats almost any real continuation, so the policy was drawn *towards* infeasible
states. The fingerprint was a spurious exact `-4.0000` = `F` at `c=uempbenefit=0.2` plus `beta*0`.

Only this tier fired it because it needs both interpolation nodes infeasible, which needs an
absorbing zero-consumption region: `a2_grid` starts at 0 and `aprimeFn` maps `(a2=0,semiz=0)` back to
0, so `a2=0` is absorbing and at `semiz=1` gives `c=0` for every `d2`. The withA1 model escapes via
the standard asset.

Fixed at all 684 ExpAssetsemiz solver sites and the 8 ValueFnFromPolicy ones (`ac60a3ae`): **zero
each product term BEFORE summing, never after**. On the 13x2x20 reproducer, `V` went from 70 NaN to
0 and the 44 states where one side was infinite and the other finite went to 0. Remaining scope —
2578 solver sites across 398 files in the other five families, 60 ValueFnFromPolicy sites, RiskyAsset
unsurveyed — is in `InterpZeroWeightNaN_proposal.md` in the toolkit repo.

Final run: **482 V_Jplus1 checks, all exactly zero**, all 8 `Inf` cleared to `4.441e-16`–`4.441e-15`.
Fig 24 (`d1_z_e_with2A1`) OOMs in its big-`a_grid` section, which sits *before* its V_Jplus1 block, so
34 of its checks never run — GPU capacity, not correctness. Moving that section after the block would
recover them.

#### Bank 12: CoreFHorzQHExpAssetsemizTests (2026-09-07) — closes the project

All 24 subcodes gained a V_Jplus1 section (**+1520 checks**). **1520/1520 exactly zero, full run, no
errors, no OOM.** Largest residual anywhere in the bank is `1.863e-09` = `2^-29`, one ULP at the
known |V| ≈ 8.4e6 poor corner. Its 32 noa1 `ValueFnFromPolicy ... Inf` — the QH counterpart of bank
11's 8 — were cleared by the same solver fix.

### V_Jplus1 project closed (2026-09-07)

Twelve banks, **8 toolkit defects, every one of them in a SemiExo path**, two of which returned
silently wrong answers rather than erroring. Suite V_Jplus1 coverage went **160 → 520 subcodes**.
Commits `337ae6b` (tests) and `ac60a3ae` (toolkit).

Two lessons outlast the project:

- **`max(abs(A-B))` ignores NaN.** Wherever both sides are NaN the check passes silently, so a
  `0.000e+00` line can be comparing nothing. This is what hid the interpolation defect for so long,
  and it means the noa1 checks in these banks were partly vacuous before the fix. When a result looks
  too clean, count the finite / `-Inf` / `NaN` census on each side separately rather than trusting a
  single max.
- **A sweep needs a loose match and a reconciled count.** A strict end-of-line regex found 672 of the
  684 ExpAssetsemiz sites; the 12 it missed carried a trailing `% [N_d2, N_a2, N_bothz]` shape
  comment — and all 12 were in the `noa1_e` raws, the tier that actually fires. A partial sweep would
  have left the most exposed files unfixed while reporting success.

Deliberate gap, unchanged: **QH + age-dependent `pi_z`/`pi_e` + `V_Jplus1` is untested across every
QH bank**, because the QH donor has no age-dependent block.

### AmbiguityAversion closed (2026-09-01)

`CoreFHorzAmbiguityTests` was written test-first on 2026-08-28 (9 files: 6 variants + 2
cross-test files, 156 exact checks + 8 warning checks + 1 error assert) and got its full solver
wave on 2026-09-01. First GPU run was green: 156/156, all asserts. The spec is
`AmbiguityAversion_testbank_proposal.md` in the toolkit repo.

What the bank covers: `{nod,d} × {z, e, z&e}` (no noz+noe — ambiguity with no shocks is
deliberately an error, and the bank asserts that it errors), each variant at all four solver
tiers (plain/DC/GI/DC+GI) with the full lowmemory ladder (=2 in the z&e variants) and
`ValueFnFromPolicy` at the plain and GI tiers, plus the big-`a_grid` GI-convergence
moments/dist section in the baseline-mirror style. Cross-tests, each in a markov-z and an iid-e
flavour: identical priors = vNM (with a flat-vs-`_J` input-form twin), 3-vs-9 duplicated priors
(min ignores multiplicity — catches count-weighting), an unambiguously-worse pi binding from
either prior slot (prior ordering irrelevant), age-varying `n_ambiguity` via `V_Jplus1`
(single-prior second half = vNM there; first half = short solve seeded with the vNM `V`), and
the pi-consistency warning firing/staying silent. Deferred by design (not test-first): semiz,
with2A, and a combined z&e cross-test flavour.

The solver wave: all 24 raws are mechanical transforms of the modern exponential raws (donor +
per-prior loop + min at the EV sites, nothing else), which retired the 6 old-style raws — one of
which, the z&e variant, had a real bug: its z-stage prior loop clobbered its own EV base, so
priors 2+ multiplied garbage. Also new: the `ValueFnIter_FHorz_AmbiguityAversion` dispatcher +
three level-2 dispatchers, `ExogShockSetup_FHorz_AmbiguityAversion` (each prior's pi runs
through the standard `ExogShockSetup_FHorz` pipeline via one recursive `gridpiboth=2` call per
prior, so every accepted pi shape and timing/trim convention applies to the priors for free; it
also throws the pi-consistency `warning()`), and `ValueFnFromPolicy_FHorz_AmbiguityAversion`.

Two conventions worth knowing when reading its diary or extending it:

- **GI is conditional on the prior** (like d/z/e): each prior's EV is interpolated over aprime
  and the min over priors is taken afterwards
  (`EVinterp=min(interp1(a_grid,ambEV,aprime_grid),[],4)`); same in `ValueFnFromPolicy` (per-prior
  L2-weighted lookup, then min). At grid nodes the orders coincide, which is why GI level 1 and
  all of DC (where EV is only ever read pointwise at grid nodes) just use the pointwise min.
- **The 18 non-zero exact checks are a new benign class**, distinct from Policy ties: 10
  `ValueFnFromPolicy` lines at 1e-15..1e-14 (FromPolicy does the z-expectation by matrix
  multiply where the raws broadcast-and-sum), and 8 cross-test-3 `V` lines (retirement flattens
  `V` in z, making the worse and original priors' EVs mathematically equal there, so `min` is
  free to return either one's rounding — the accompanying Policy lines are exactly zero).

**The with2A tier followed on 2026-09-02** (spec: `AmbiguityAversion_with2A_proposal.md`,
test-first, one debug iteration): 6 with2A variants (figs 7-12, the QH bank's 2A setup reused) +
2 cross-test files (identical-priors-=-vNM and V_Jplus1 both AT ALL FOUR TIERS, worse-pi-binds at
plain), plus a companion fix repeating the main-tier cross-test 4 at the DC/GI/DC+GI tiers. The
bank now totals 416 checks and is fully green: 374 exact zeros, 42 ULP-floor lines in three
benign classes — the two from the main tier plus GI2A-vs-DC2A_GI2A at ~1e-14 in the six 2A
subcodes, which is the donors' own discount-before-vs-after-interpolation asymmetry surfacing
under `%.3e` (the exponential banks hide it under `%2.8f`). Toolkit side: 18
`AmbAverse_{DC2A,GI2A,DC2A_GI2A}` raws, 2A branches in the three level-2 dispatchers, and the
GI2A block in `ValueFnFromPolicy_FHorz_AmbiguityAversion` (a2prime folded into the linear index,
QH-style).

**The riskyasset tier followed on 2026-09-02** (spec: `AmbiguityAversion_RiskyAsset_proposal.md`
in the toolkit repo; bank: `CoreFHorzRiskyAssetTests/withAmbiguityAversion/`, 220 checks, fully
green after three debug iterations; originally 240 — the 20 lowmemory checks in the four noz+noe
subcodes were removed 2026-09-02 as vacuous: with zero shocks lowmemory is silently discarded, so
each check re-solved the identical problem and compared it to itself). u is treated as
AMBIGUITY, not risk: `ambiguity_pi_u` ([N_u, max(n_ambiguity)]) is mandatory in this
combination — the agent does not know the risky
return distribution — and noz+noe becomes a valid ambiguity model (pure return ambiguity), unlike
the standard-asset family. Sequential mins innermost-first: e on the grid, z at the a2 lottery
(the lottery — the riskyasset analog of the standard family's aprime interpolation — is
conditional on the prior), u last; the d2 riskyshare refinement max comes after all the mins.
Tiers: noa1 (base only, structural) and withA1 (all four methods); 40 raws + 4 dispatchers +
FromPolicy base/GI + the setup subfn's riskyasset block. T1 (identical priors = vNM) is bit-exact
at all four methods in both d-variants via running-argmin component tracking. The worse-pi_u-binds
cross test is the marquee economics: the ambiguity-averse investor behaves as if facing the worst
return distribution. Three iteration finds worth remembering: dead `aprimeProbsK=aprimeProbs;`
initializers (the wrap assumed a variable the riskyasset donors never define — they build the
probs fresh inside the segment); the e-min stage emitting `EV` where the z-stack reads `EVpre`;
and the squeeze-form u-min collapsing `EV` in place (priors 2+ summed garbage — same bug class
the main tier caught in its e-only raws). One documented DESIGN DEVIATION: the withA1 a1
grid-interp layer is min-then-interp (all mins and the d2 max at the coarse a1prime nodes, then
`interp1` of the worst-case EV) — moving the mins to the fine grid would drag the d2 max there
too and break T1 against the vNM donors' max-at-coarse-then-interp; FromPolicy GI mirrors the
raws with per-corner min pipelines. Deferred: semiz (as everywhere in AA) and with2A1 for
riskyasset (explicit errors).

The former honest gap (DC/GI-tier `V_Jplus1` branches unexercised) is CLOSED — and closing it
paid immediately: the companion GI leg caught a real bug on its first run. `interp1` returns the
**query's** shape when its value input is a vector, and the e-only GI1 raws' per-prior `ambEV`
`[N_a,n_ambiguity(jj)]` collapses to a vector exactly when age-varying `n_ambiguity` hits 1 —
row query, row output, `min(...,[],2)` crushed it to a scalar, index error. Fixed by transposing
the query to a column in all 24 GI1/DC1_GI1 interp sites (bit-identical where the input stays a
matrix/N-D; the 2A GI families were immune, their `a1prime_grid` being column-shaped already).
Every `V_Jplus1` check in the family is now an exact zero at every tier, 1A and 2A. The trap
generalizes: any raw that interpolates a per-prior/per-something EV array whose trailing
dimension can reach 1 has this failure mode.

### Reading the diaries: the ULP floor

`%2.8f` cannot distinguish a true zero from one unit in the last place, and `V` legitimately
reaches ~1e7 at the poor corner of these grids (`sigma=2` makes the ReturnFn `1-1/c`, and the
cubic `a1_grid` puts its second point at 5e-6). One ULP there is `2^-28 = 3.7e-09`. A diary
reporting `0.00000001` may be showing 1.5 ULP.

All eight **FHorz QH banks** were converted to `%.3e` in `bf158c7` (15 700 checks) and now surface
that floor honestly — expect small non-zeros as normal. The exponential banks still print `%2.8f`,
so their "exactly zero" counts include an unknown number of sub-5e-9 values. Their true state is at
least as good as reported, but not verified to be exactly zero.

**How much the formatter was hiding, measured.** QH ExpAssetze was recorded in this doc as
"776 checks, 12/12 figures", all zero, on 2026-08-18 under `%2.8f`. Its first `%.3e` run
(2026-09-01) reports 148 non-zeros — every one of them at the floor, worst case `5.588e-09`,
i.e. 1.5 ULP. Nothing changed in the numbers; the earlier all-zero reading was an artifact of
the format string. Do not read a pre-`bf158c7` "all zero" as stronger than "all below 1e-8".

That run also shows what the floor looks like when you have enough of it to see the shape: the
non-zeros take only ten distinct values across 920 checks, all dyadic. A fault does not produce
a dyadic histogram. Three checks stay exactly zero even at the floor and are the ones worth
grepping — `lowmemory`, the β₀=1 exponential cross-checks, and anything ending `(Policy)`.

Any relative-error diagnostic in these banks must divide by `max(abs(V(isfinite(V))))`: the
ReturnFn returns `-Inf` wherever `c<=0`, so the plain max is `Inf` and silently reports 0.

### Known open items

- **V_Jplus1**: now 160 subcodes (recounted 2026-09-01 off the working tree — the earlier "80"
  missed the EZ side entirely): all 32 in `CoreFHorzTests_subcodes`, all 32 in its QH mirror
  (committed `daa90d5`), all 16 in its EZ mirror, 32 in the RiskyAsset main bank plus 32 in the
  RiskyAsset EZ mirror (both banks' noa1 and withA1 tiers incl. semiz; the with2A1 tier has
  none), and all 16 in the ResidAsset bank (240 checks, unrun).
  ResidAsset was the first family to get V_Jplus1 coverage from day one rather than retrofitted;
  do the same for any new bank (AmbiguityAversion did: cross-test 4 covers its V_Jplus1
  branches, since 2026-09-02 at ALL solver tiers in both its 1A and 2A models — and that
  coverage caught a real interp1 shape bug on its first run; see its section). So all three
  baseline preference mirrors now have V_Jplus1 coverage: QH and EZ per-variant (retrofitted),
  AA via cross-test 4 at every tier. **The ExpAsset family is no longer a gap at all: the
  twelve-bank V_Jplus1 project closed on 2026-09-07** (see the section below), giving every
  variant subcode in all six families — and each family's QH mirror — a V_Jplus1 block. The
  suite total went 160 → 520 subcodes.

  **The exposed-raw count was wrong and is now measured.** This entry previously said 660
  ExpAsset-family raws carry a `V_Jplus1` branch, itemised per family; those per-family figures
  were internally inconsistent (some were core-only, some were not) and the total was low. A
  direct `grep -l V_Jplus1` over every `*_raw.m` in the six families returns **1440 of 1440** —
  i.e. *every* raw carries one: ExperienceAsset 384, ExperienceAssetu
  384, ExperienceAssete 192, ExperienceAssetz 192, ExperienceAssetze 96, ExperienceAssetsemiz
  192. Those branches are age-shifted copies of the in-loop code — the shape that produced the
  `jj`/`N_j` bug. Still the largest gap on any axis, by more than double what this doc claimed,
  though bank 1 has now put runtime coverage behind part of ExperienceAsset's 384.

  The QH ExperienceAssetu port put a number on the cost. Eleven real defects were found in
  shipped exponential code by diffing each `V_Jplus1` branch against the in-loop code of the same
  file, per lowmemory branch — a missing `squeeze` after the `pi_u` contraction, `EVpre`
  contracting dim 2 where `shiftdim(...,-2)` requires dim 3, a dead no-`u` lottery block calling
  the wrong builder one argument short, a 13-parameter function handed 15 arguments, a `loweredge`
  index list one short of a 6-D `maxindex1` (fusing `z` with `e`), an unbound `zind` on the
  `lowmemory==3` path, and a missing `repelem` pair. **Two of them returned silently wrong answers
  rather than erroring.** None was reachable by any bank, because no ExpAsset-family bank sets
  `V_Jplus1`. That is what this gap costs.

  **The `loweredge` defect then recurred, and the recurrence was also undetectable.** The same
  five-subscript index against a 6-D `maxindex1` was found in four more `DC1_e` raws on
  2026-09-01 and fixed in `6031762b` — two of them in *shipped exponential* code
  (`ExpAsset_DC1_e_raw`, `ExpAssetz_DC1_e_raw`), one in each of the QH ExpAsset N/S pair. All
  five sites are in `V_Jplus1` branches, so the two green bank runs that same day say nothing
  about them: **this fix shipped untested, necessarily.** It is confident by inspection only —
  the terminal and backward arms at the same position in the same files already use six
  subscripts. The operational lesson: `V_Jplus1` drift replicates across sibling families, so
  when one of these turns up, sweep the siblings immediately rather than fixing the file in hand.
- ~~**Staleness.**~~ **Closed 2026-09-01.** All three exposed banks — `QHExpAssetz`,
  `QHExpAssetze`, `QHExpAssete` — were re-run that day and are green. They had last run 18–19
  August, before the beta0 refactor, the QH `ValueFnFromPolicy` hierarchy change and the two
  fixes in `d7dfcd7c`; the hierarchy change had **turned on** their
  `experienceasset{e,z,ze}`+semiz paths, which previously routed through the generic SemiExo
  router. Same destination and same arguments, and the runs confirm it: `QHExpAssete` returns
  **exactly the same 1624 checks across 24/24 configs** as its August run, with no errors. No
  FHorz bank now predates the refactor.
- **No QH bank is test-first any more.** `CoreFHorzQHExpAssetTests` was the last one, and it
  closed on 2026-08-22. Status re-read off the working tree 2026-08-25, not off the banners:
  - `CoreFHorzQHExpAssetTests` (48 subtests) — **implemented and substantially validated.**
    `ExperienceAsset` now carries **256 QH raws** under
    `ValueFnIter/FHorz/ExperienceAsset/QuasiHyperbolic/` (8 directories x 32:
    `{base, DivideConquer, GridInterpLayer, DivideConquerGridInterpLayer}` x
    `{nosemiz, ExpAssetSemiExo}`), plus 8 dispatchers, with live routing at
    `ValueFnIter_Case1_FHorz.m:388`/`:398`. Landed across `d634db8b`, `b6bc5acf` and `5cada730`
    (the last being the 2A1 tiers, 96 raws across DC2A/GI2A/DC2A_GI2A). The earlier claim in this
    doc that the family had "0 QH raws" and that "every subtest errors at its first ValueFnIter
    call" was true on 2026-08-19 and is now wrong — do not trust it.
    The 2026-08-22 run reached **figs 1-43 of 48, 217 checks**, then aborted in
    `CoreFHorzQHExpAsset_d1_z_noe_semiz_with2A1` (fig 44) on a genuine
    `Out of memory on device` from `gpuArray/arrayfun` in `CreateReturnFnMatrix_ExpAsset_Disc`.
    That is GPU capacity, not correctness — the same OOM ceiling every ExpAsset bank hits on its
    largest 2A1+semiz variant. Figs 44-48 remain unrun for that reason.
    **Three lines in that diary are nonzero and all three are fine**: `Naive DC2A (Policy)` = 2,
    `Naive DC2A (Policyalt)` = 1, `Sophisticated DC2A (Policy)` = 1. Each is followed by its own
    diagnostic — `differs at 13 of 210080 state-age points; max|dV| there = 0 (0 => tied optima,
    not a wrong argmax)` — so `V` agrees exactly and only the argmax tie-break differs, at 13 and
    3 points respectively, all at age 14 in dim-1 index range 87..100 of 101. A plain
    grep-for-nonzero flags this diary; read the nine `[diag]` lines before concluding anything.
  - `CoreFHorzQHExpAsseteTests` (24 subtests) — **done and GPU-green.** The toolkit has 128 QH
    raws under `ExperienceAssete/QuasiHyperbolic/`, four `QuasiHyperbolicExpAssete*` dispatchers
    plus four SemiExo ones, and live routing at `ValueFnIter_Case1_FHorz.m:401`/`:411` (committed
    `b6bc5acf`). Re-run 2026-09-01 and green: **1624 checks, 24/24 figures, no errors** — Naive and
    Sophisticated on V / Valt / Policy / Policyalt across each lowmemory ladder, plus 24
    `ValueFnFromPolicy` oracle checks. Identical check count to the 2026-08-18 run, which
    confirms the beta0 refactor and the `ValueFnFromPolicy` hierarchy change left this family
    alone. 1296 exactly zero and 328 at the ULP floor (worst `5.588e-09`), every one of them on
    `ValueFnFromPolicy` or a cross-method comparison; the August run reported the same values
    under `%2.8f`, where all but two of them rounded to `0.00000000`.
    One toolkit limit remains but is not reachable from this bank:
    `QuasiHyperbolicExpAsseteSemiExo_{DC,GI,DC_GI}` error on `N_a1==0`, and the noa1+semiz
    block (figs 5-8) is base-method only.
  - `CoreFHorzQHExpAssetzTests` (24 subtests) — **done and GPU-green, re-run 2026-09-01.**
    `ExperienceAssetz` carries **128 QH raws** (64 nosemiz + 64 semiz) after `6031762b` added the
    24 nosemiz 1A raws; the `QuasiHyperbolicExpAssetzSemiExo` dispatcher and the three
    `{DC,GI,DC_GI}` sub-dispatchers all exist, and no `not yet implemented` error survives
    anywhere in the family. The current run is **1512 checks, 24/24 figures, no errors** — 1248
    exactly zero and 264 at the ULP floor (worst `3.725e-09` = 1 ULP), all of them on
    `ValueFnFromPolicy` or cross-method comparisons. Covers all four solver tiers x nod1/with-d1
    x no-e/with-e x Naive/Sophisticated, including the with-e 2A paths at lowmemory 3 and the
    beta0=1 degeneracy checks that collapse QH onto the exponential.
    `ValueFnFromPolicy_FHorz_QuasiHyperbolic_ExpAssetz_SemiExo_GI` was added earlier
    (`2f291d33`) to close the last gap — ExpAssetz had been the only family missing a SemiExo_GI
    value-from-policy. The earlier 1268-check figure predates the 1A tiers and the `%.3e`
    conversion; do not compare the two runs' zero counts directly.
  - `CoreFHorzQHExpAssetzeTests` (12 subtests) — **done and GPU-green, re-run 2026-09-01.**
    `ExperienceAssetze` carries **64 QH raws** after `6031762b` added the 12 nosemiz 1A raws and
    three sub-dispatchers, the latter taking over the 2A routing that had lived in the family
    dispatcher. **920 checks, 12/12 figures, no errors**, 772 exactly zero and 148 at the floor.
    Its previously recorded "776 checks, all zero" was measured under `%2.8f` — see the ULP-floor
    section; the numbers did not change, the format string did.
  - `CoreFHorzQHTests` (baseline, 32 subtests) — 4016 checks, none nonzero. Note this diary uses
    a different closing-marker style from the ExpAsset-family banks, so the per-figure counting
    used above does not apply to it.
  - `CoreFHorzTPathQHTests` — bank exists (440 check sites) but has **never produced a
    diary**. See the FHorz transition-path section at the end of this doc: it is one of five
    TPath banks in that state.
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
- **ResidualAsset blocker: FIXED in the working tree 2026-08-25, not yet run.**
  `CreateResidualAssetFnMatrix_Case1` used `d_gridvals` before defining it; the parameter is now
  named `d_gridvals` and the redundant `CreateGridvals` line is gone (the callers already pass
  gridvals). `LifeCycleModel38` is the smoke test and has not been run yet.
  Also fixed alongside it: the `rprimeFnParamNames` split now counts `semiz` and `e`
  (`l_d+l_a1+l_a1+l_semiz+l_z+l_e`) in the VFI *and* the StationaryDist copy, and
  `divideandconquer=1`/`gridinterplayer=1` now error for `residualasset` instead of being
  silently ignored.
- **ResidualAsset: 60 raws to write.** 4 exist as of 2026-08-25 (all unrun): the two z-only
  nosemiz raws (`_raw`, `_nod1_raw` — renamed from `_nod_raw` for ExperienceAsset parity, both
  now closing their `lowmemory` ladder with an `else error`), plus two new SemiExo z-only raws
  and a `ValueFnIter_FHorz_ResidAssetSemiExo` dispatcher. Remaining order:
  the 4 `noz`/`e` base raws + 4 SemiExo → DC1 (16) → GI1 (16) → DC1_GI1 (16), with cross-test 4
  (vs ExperienceAsset) as the acceptance test at each tier.
- **`ValueFnFromPolicy_FHorz` is the real gate on the ResidAsset bank, not the raws.** Every
  variant subcode calls it four lines after its first solve, so with no `residualasset` branch
  there *no* variant subtest can print anything, however many raws exist. Until it is written,
  the only runnable things are `LifeCycleModel38` and the base-solver tier of cross-test 4.
  It should be done before, not after, the remaining raws.
- **The `noz`/`e` base raws need two builder files first.** ExperienceAsset has the full
  `{_Disc, _Disc_e, _Disc_noz}` ReturnFnMatrix family; ResidualAsset has only `_Disc`, so
  `CreateReturnFnMatrix_ResidAsset_Disc_noz` and `..._Disc_e` have to be written (~1,650 lines
  of arrayfun dispatch between them), and `CreateResidualAssetFnMatrix_Case1` has to dispatch to
  them rather than hardcoding `_Disc`.
- **ResidualAsset downstream gaps**, in the order the subcodes need them: `StationaryDist`
  `N_z==0` paths and a `..._ResidAssetSemiExo` sibling and the `gridinterplayer==1` path (all
  three currently `error`); semiz and GI shapes for `EvalFnOnAgentDist`/`LifeCycleProfiles`; a
  `residualasset` branch in `ValueFnFromPolicy_FHorz`; and SimPanel, which is the only one that
  is genuinely new code rather than a shape fix (a residual-asset panel has to evolve `r` by
  drawing the lower/upper `r_grid` point with `rprimeProbs`).
- **`divideandconquer=1` and `gridinterplayer=1` are silently ignored with `residualasset`** —
  `ValueFnIter_Case1_FHorz` returns at line 544 before either is consulted, so the user gets a
  base solve and no warning. `StationaryDist` *does* error on `gridinterplayer==1`, so the two
  halves of the toolkit disagree about whether GI is an error or a no-op.
- **`lowmemory=2` silently returns `V=zeros` in both residual-asset raws** — the
  `if lowmemory==0 ... elseif lowmemory==1` has no `else`. Every new raw should close the ladder
  with an `else error(...)`.
- **`CreaterprimePolicyResidualAsset_Case1` converted to `z_gridvals` 2026-08-25 (unrun).** It was
  the only member of the `CreateXprimePolicy*` family taking a stacked `z_grid` — every sibling
  (`...ExperienceAssetz`, `...ze`, `...InheritanceAsset`, and the `_J` variants) takes gridvals.
  It sliced marginals onto *separate* dimensions, which silently tensor-products them and so
  cannot represent a joint grid: a discretised VAR1 would have been read as garbage with no error.
  Now every z component sits on one dimension, `shiftdim(z_gridvals(:,k),-l_a-l_r)`. Invisible at
  `l_z==1`, which is every residual-asset model that exists, which is why it never bit.
  Fixed in the same pass: `z5vals` was referenced by 20 `arrayfun` branches and never assigned;
  and `StationaryDist_FHorz_ResidAsset` passed a z-only grid alongside `n_ze`, so the helper
  indexed past the end of it as soon as the model had an `e`. The subfn now takes `z_gridvals_J`
  (matching `StationaryDist_FHorz_ExpAssetzSemiExo`) and builds the joint (z,e) gridvals itself.
- **`d1vals(1,1,1,1)=d_grid(1)` removed 2026-08-25.** It sat in the `l_d==1 && l_a==1 && l_z==1`
  branch of `CreaterprimePolicyResidualAsset_Case1` — one of ~100 — and overwrote a single
  policy-derived `d` value with the first grid point, commented only `% Requires special
  treatment`. No reading of the surrounding code makes it correct; it looks like a debugging
  leftover. If it was in fact papering over something, cross-test 4 (ResidAsset vs ExperienceAsset)
  will catch it: a one-element `rprime` discrepancy shows up as a non-zero `StationaryDist`
  comparison against a family that does not have the line.
- **One residual-asset defect still open.** `CreaterprimePolicyResidualAsset_Case1` guards its
  `a4` case with `if l_a>=1` (should be `l_a>=4`), so `l_a==3` reaches an undefined `a4grid`.
  Unreachable from the bank as scoped (needs three standard endogenous assets).
  Also still cosmetic-but-misleading: `CreateResidualAssetFnMatrix_Case1` names its parameter
  `z_grid` while every caller passes `z_gridvals_J(:,:,jj)`.
- **OOM ceiling**: every ExpAsset bank's largest `d1_z_e` + 2A/semiz variant exceeds GPU
  memory. Where it lands matters: in ExpAsset/ExpAssetU it hits a *leading* brute-force
  baseline and takes the rest of the script with it; in ExpAssete/ExpAssetz/ExpAssetsemiz
  it hits the *trailing* big-grid block, after that subtest's exact checks have printed.

## FHorz transition-path (TPath) banks

Everything above is the *stationary* FHorz suite. The transition-path banks are a separate set,
excluded from the main table and from every count in it (the `panel vs dist` note already says
"TPath banks excluded throughout"). Surveyed 2026-09-05.

| bank | variants | cross-tests | subcodes | check sites | of which `close to zero` | diary |
|---|---|---|---|---|---|---|
| `CoreFHorzTPathTests` | 8 nosemiz + 8 semiz | 6 | 22 | 554 | 32 | **2026-08-27, complete** |
| `CoreFHorzTPathTests/withQuasiHyperbolicDiscounting` | 8 | 2 | 10 | 440 | 32 | never |
| `CoreFHorzTPathExpAssetTests` | 8 | — | 8 | 256 | 16 | never |
| `CoreFHorzTPathExpAssetzTests` | 4 | 2 | 6 | 168 | 8 | never |
| `CoreFHorzTPathTwoEndoTests` | 8 | — | 8 | 256 | 16 | never |
| `CoreFHorzTPathPTypeTests` | 7 tests | — | 7 | 44 ‡ | 9 | never |
| **total** | | | **61** | **1718** | **113** | |

‡ static `fprintf` sites. The PType checks sit inside `for ii=1:N_i` loops, so at runtime it is
about **73** — `N_i=2` everywhere except `ShockTests_4types`, which uses 4.

**1718 check sites are written; 554 of them — 32% — have ever executed.** That is the headline.
Five of the six banks have never produced a diary, and `CoreInfHorzTPathTests` (out of scope here)
is the only other TPath bank in the suite that has.

### What a TPath variant tests that a stationary variant does not

The extra axis is **`transpathoptions.fastOLG` ∈ {0,1}**, crossed with the same
`{base, DC, GI, DC+GI}` ladder and the same `lowmemory` ladder as the stationary banks — hence
the `Divide-and-conquer (slowOLG)` / `(fastOLG)` / `(with GI, slowOLG)` / `(with GI, fastOLG)`
quartet on every label. On top of that each variant runs four oracles:

- **`Do nothing TPath`** — a constant `PricePath`/`ParamPath` must reproduce the stationary solve
  exactly, checked on `V`, `Policy` and `AgentDist`. This is the load-bearing oracle of the whole
  TPath suite: it is the only one that pins the path code to an independently-computed answer
  rather than to another tier of itself.
- **one iteration of the GE transition path** (`transpathoptions.maxiter=1`), with/without
  fastOLG and with/without GI. Explicitly a shape check — the subcodes say so — the GE core
  being tested in `CoreStationaryGeneralEqm`.
- **big-`a_grid` GI convergence** — `StationaryDist with/without grid interp … close to zero`
  plus the eyeballed moment rows, the same construction as the stationary banks' QH/EZ mirrors.
- **cross-tests** — `z as e`, `z and e 1/2`, `semiz as z`, and (ExpAssetz) a fake-`z`
  experienceassetz against plain experienceasset.

The QH mirror adds the Naive/Sophisticated split with `Valt`/`Vunderbar`/`Policyalt` channels and
a `beta0=1` collapse onto exponential, exactly mirroring the stationary QH banks.

### The one run there is

`CoreFHorzTPathTestsdiary.txt`, 2026-08-27, is **clean and complete**: all 16 figures, all 554
check sites reached (the static and diary counts agree exactly), **534 printing `0.00000000`**,
16 `close to zero` GI-convergence lines, and 4 lines at `0.00000001`. No errors, no OOM, no
aborted script — this is the only TPath bank known to run end to end.

Two caveats on reading it:

- **It predates the `%.3e` conversion.** All five TPath banks were moved onto `%.3e` in `917bb7c`
  (2026-09-01), but this diary is from 08-27, so it is a `%2.8f` reading and the ULP-floor
  section above applies in full: its "534 exactly zero" means "534 below 5e-9", and the four
  `0.00000001` lines are floor-level, not faults. The next run will report honest small non-zeros
  and that is expected, not a regression.
- **It is still representative.** The only substantive change to the bank since is
  `QHadditionaldiscount` going from `{'beta0'}` to `'beta0'` in the *QH* driver (`b119a89`);
  everything else in the 33 touched files is format strings.

### The header comment about semiz was stale (fixed 2026-09-05)

`CoreFHorzTPathTests.m` said `with/without semiz [tests built; toolkit does not yet implement
FHorz TPath with semiz]` at the top and `NOTE: these require the toolkit to implement FHorz TPath
with semiz (not yet done)` above the second half. **Both were wrong.**
`ValueFnIter/TransPathFHorz/subcodes/ValueFnOnTransPath_FHorz_SemiExo.m` exists and is dispatched
from `ValueFnOnTransPath_Case1_FHorz.m:178`, and the 08-27 run executed figures 9–16 and all 231
semiz check sites with no errors — the comments would have led someone to skip the working half
of the script. Both replaced with the dispatch fact.

The `NOT YET IMPLEMENTED` comments in `CoreFHorzTPathPTypeTests` were checked at the same time and
are all still accurate: neither `EvalFnOnTransPath_AllStats_Case1_FHorz_PType` nor
`LifeCycleProfiles_TransPath_FHorz_Case1_PType` exists, and nor does a non-PType
`EvalFnOnTransPath_AllStats_Case1_FHorz` (only InfHorz has one).

### Six silent-wrong-answer holes on the TPath side — all guarded 2026-09-06

The TPath entry points dispatch on `experienceasset` / `experienceassetz` / `n_semiz` and nothing
else, so every other non-standard combination fell through to the standard-endogenous-state code
and returned a different model without complaint. Four such combinations existed; all now error.

| combination | was | now |
|---|---|---|
| `riskyasset` on a TPath | zero mentions anywhere in the TPath tree | `error` |
| `residualasset` on a TPath | zero mentions anywhere in the TPath tree | `error` |
| `experienceasset(z)` + `semiz` | semiz tested *before* asset type, so it reached `ValueFnOnTransPath_FHorz_SemiExo`, which has zero knowledge of experience assets | `error` |
| `QuasiHyperbolic` + `semiz` (compute path) | QH dispatched before the semiz branch; the QH subfn never mentions semiz | `error` |
| `QuasiHyperbolic` + `semiz` (GE path) | the GE path has its *own* QH implementation, which never reads `n_semiz` | `error` |
| `QuasiHyperbolic` + `experienceasset(z)` (GE path) | the GE substep dispatcher returns at the expasset branch *before* the QH branch, and the ExpAsset substep has zero mentions of `exoticpreferences` — so QH was silently dropped, while the compute-only path errors on the same model | `error` |

Guards added at all five entry points: `ValueFnOnTransPath_Case1_FHorz`,
`ValueFnOnTransPath_FHorz_QuasiHyperbolic`, `AgentDistOnTransPath_Case1_FHorz`,
`TransitionPath_Case1_FHorz`, and `TransitionPath_Case1_FHorz_PType` — the last needs its own copy
because it builds per-type option structures and drives the substeps directly instead of routing
through `TransitionPath_Case1_FHorz`. The PType value-fn and agent-dist commands do route through
their guarded non-PType versions, so they inherit.

### riskyasset and residualasset: the detail

`riskyasset` and `residualasset` appear in **zero files** under `ValueFnIter/TransPathFHorz/`,
`StationaryDist/TransPathFHorz/`, `EvaluateFnOnAgentDist/TransPathFHorz/` and `TransitionPaths/`.
Not "errors with a message" — absent. `ValueFnOnTransPath_Case1_FHorz` never reads
`vfoptions.riskyasset`, so such a model is dispatched down the standard-endogenous-state path with
`n_u`/`pi_u`/`aprimeFn`/`refine_d` simply unread, and a different model comes back without
complaint. Every neighbouring unsupported combination guards — `experienceasset`+noa1
(`ValueFnOnTransPath_FHorz_ExpAsset.m:49`), `experienceassetz`+z-varying
(`TransitionPath_Case1_FHorz.m:366`), `experienceasset`+QH
(`ValueFnOnTransPath_FHorz_QuasiHyperbolic.m:21`) — these two do not.

`TPath_RiskyAsset_proposal.md` (toolkit repo, 2026-09-05) scopes both the guards and a full
RiskyAsset TPath family: 80 SingleStep raws, a 16-variant test bank, and the reason `withA1` has
to be built before `noa1`.

### QH + semiz is silently wrong, not an error

`ValueFnOnTransPath_Case1_FHorz.m:108-111` dispatches to
`ValueFnOnTransPath_FHorz_QuasiHyperbolic` and **returns**, which is *before* the semiz branch at
`:178`. The QH subfunction never mentions `semiz` — zero occurrences in the file. So
`exoticpreferences='QuasiHyperbolic'` together with `vfoptions.n_semiz` does not error; it
computes the no-semiz answer and hands it back.

The QH driver's own banner already calls this combination `[NOT SUPPORTED by
ValueFnOnTransPath_FHorz_QuasiHyperbolic]`, so the intent is documented — it is the *enforcement*
that is missing. The idiom to copy is three lines away in the same file:

```matlab
if vfoptions.experienceasset>=1
    error('ValueFnOnTransPath_FHorz_QuasiHyperbolic: experienceasset not yet supported')
end
```

`experienceasset` gets a guard; `semiz` does not. One `prod(vfoptions.n_semiz)>0` error beside it
closes it.

### Open items recorded in the drivers themselves

These are annotations the drivers already carry, not new findings — but they are the only record
of them, and they are the reason four of these banks have no diary worth keeping:

- **`CoreFHorzTPathExpAssetTests`, figures 2 and 4** — `RUNS BUT: Policy differs by 2, Claude
  claims it is just about how DC handles indifferent policies different from without DC`. That
  is the standard argmax-tie story and is probably right, but it is asserted, not demonstrated:
  the stationary banks settle these by showing `max|dV|` at the differing points is 0. The
  subcode has the diagnostic block for it (`(d2_nonDC, d2_DC) -> count`, diff counts by `t` and
  by `j`); nobody has run it and written the answer down.
- **`CoreFHorzTPathExpAssetzTests`, figure 4** — `Some of the lowmemory are not quite right,
  seems to be just the Policy (V seems fine); at first I thought L2flag, but appears to impact
  the DC (without GI) so that is not the reason.` **This one is unexplained and is the most
  likely real defect in the TPath suite.** A `lowmemory` disagreement is the check the stationary
  banks treat as never tied — "every `lowmemory` check exactly zero" is one of the three
  discriminators that separate the ULP floor from a fault.
- **OOM is routine here.** Four of the largest configurations replace `a_grid_big` with a
  hand-built `n_a_notsobig` (501, 501, 301, 251 points) to fit. Same ceiling as the stationary
  ExpAsset banks, hit earlier because a path holds `T` periods at once.
- **`CoreFHorzTPathPTypeTests`** is blocked on two toolkit functions that do not exist:
  `EvalFnOnTransPath_AllStats_Case1_FHorz_PType` and
  `LifeCycleProfiles_TransPath_FHorz_Case1_PType`. The corresponding blocks in
  `..._PerTypeFnsToEvaluate.m` are commented out awaiting them. Also noted there: **per-type
  `N_j` is not allowed on a TPath**, unlike the stationary PType commands, so there is
  deliberately no `PerTypeNj` test.

### Coverage gaps relative to the stationary suite

Whole axes that the stationary banks cover and the TPath banks do not:

| axis | stationary | TPath |
|---|---|---|
| asset tiers | noa1 / withA1 / with2A1 | withA1 only (plus a separate TwoEndo bank) |
| families | ExpAsset, U, e, z, ze, semiz, RiskyAsset, ResidAsset | ExpAsset, ExpAssetz only |
| exotic prefs | QH, EZ, AA, GP | QH only, and unrun |
| semiz | every family | main bank only (ExpAsset/ExpAssetz drivers say `NOT YET IMPLEMENTED`) |
| `V_Jplus1` | 520 subcodes | n/a — the path's terminal condition is `V_final`, exercised by every variant |

The `V_Jplus1` row is the one place the TPath banks are structurally *better* off: `V_final` is a
mandatory input to every path solve, so there is no separate rarely-taken branch of the kind that
produced eleven defects in the stationary families.

**The cheapest thing to do next is run the five banks that have never run.** They are written,
committed, and on `%.3e`; the main bank's clean 554 says the shared setup and the shared oracle
shapes work. Order by expected yield: ExpAssetz (a recorded unexplained `lowmemory` Policy
difference), then ExpAsset (an asserted-but-unverified tie), then the QH mirror (440 sites, never
executed once), then TwoEndo and PType.
