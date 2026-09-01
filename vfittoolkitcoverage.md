# VFI Toolkit — test coverage

Coverage of the FHorz core test banks. Tier is assigned from each subcode's directory
(`Noa1_subcodes` / `With2A1_subcodes` / else), not its filename — filenames misattribute
cross-tests. Format in the tier columns is `variants + cross-testsx`.

`ResidAsset` is the one bank with a single tier, so it has no tier subdirectories at all: its 16
withA1 subcodes sit at the top level of `..._subcodes/` and in `Semiz_subcodes/`. A tier rule
that keys off directory names must fall through to withA1 for it, not skip it.

Last updated: 2026-09-01

## CoreFHorzTests + the ExpAsset family + RiskyAsset + ResidAsset

| bank | noa1 | withA1 | with2A | var | x-tests | panel vs dist | V_Jplus1 | core raws | QH | EZ | AA |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CoreFHorzTests | — | 16+6x | 16+4x | 32 | 10 | 32/32 | 32 + 32 QH + 16 EZ | 142 | 42 | 24 | 8 |
| ExpAsset | 16+10x | 16+14x | 16+2x | 48 | 26 | 48/48 | 0 | 128 | 48 | — | — |
| ExpAssetU | 16+10x | 16+14x | 16+2x | 48 | 26 | 48/48 | 0 | 128 | 48 | — | — |
| ExpAssete | 8+16x | 8+8x | 8+2x | 24 | 26 | 24/24 | 0 | 64 | 24 | — | — |
| ExpAssetz | 8+16x | 8+8x | 8+4x | 24 | 28 | 24/24 | 0 | 64 | 24 | — | — |
| ExpAssetze | 4+4x | 4+16x | 4+2x | 12 | 22 | 12/12 | 0 | 32 | 12 | — | — |
| ExpAssetsemiz | 8+4x | 8+6x | 8+2x | 24 | 12 | 24/24 | 0 | 64 | 24 | — | — |
| RiskyAsset | 16+0x | 16+25x | 16+4x | 48 | 29 | 48/48 | 32 + 32 EZ | 128 | none | 50 | — |
| ResidAsset | n/a | 16+22x | — | 16 | 22 | 16/16 ‡ | 16 | 4 ‡‡ | none | — | — |
| **total** | | | | **276** | **201** | **276/276** | **160** | | **222** | **74** | **8** |

477 subtests in the main banks, plus 222 QH, 74 EZ and 8 AA = **781**.

**Every experience-asset family now has a QH mirror with a solver behind it, and every mirror is
complete at raw level.** ExpAssetU (48) and ExpAssetsemiz (24) were the last two showing `none`;
both are implemented and have run. ExpAssetz and ExpAssetze were the last two whose mirrors were
short of a full raw set — the missing nosemiz 1A tiers landed in `6031762b` (2026-09-01) and both
banks are green. All six families now sit at exactly 2× their core raw count.

‡ written but not runnable — `SimulateTimeSeries/` has no `residualasset` support at all.
‡‡ the bank is written test-first against 64 raws; the toolkit has 2. See below.

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
  AA mirror got its solver on 2026-09-01 — see the AmbiguityAversion section. AA is 8 rather
  than a QH-style 42 because ambiguity has no noz+noe variant (deliberately an error), and its
  semiz and with2A tiers are deferred by design, not test-first.)

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
EZ and 0 of 9 AA subcodes. Dist-derived output is not absent, but it is confined to two banks (TPath banks
excluded throughout):

| mirror / bank | files | asserted | displayed |
|---|---|---|---|
| QH — CoreFHorzTests | 43 | 120 | 512 |
| QH — ExpAsset, ExpAssete, ExpAssetz, ExpAssetze | 112 | **0** | **0** |
| EZ — CoreFHorzTests | 73 | 120 | 384 |
| EZ — RiskyAsset | 147 | 100 | 240 |
| AA — CoreFHorzTests | 9 | 6 | 48 |

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

| family | total raws | QH | EZ | AA | core |
|---|---|---|---|---|---|
| baseline (FHorz, excl. asset families) | 430 | 224 | 40 | 24 | 142 |
| ExperienceAsset | 384 | 256 | 0 | 0 | 128 |
| ExperienceAssetu | 384 | 256 | 0 | 0 | 128 |
| ExperienceAssete | 192 | 128 | 0 | 0 | 64 |
| ExperienceAssetz | 192 | 128 | 0 | 0 | 64 |
| ExperienceAssetze | 96 | 64 | 0 | 0 | 32 |
| ExperienceAssetsemiz | 192 | 128 | 0 | 0 | 64 |
| RiskyAsset | 208 | 0 | 80 | 0 | 128 |
| ResidualAsset | 4 | 0 | 0 | 0 | 4 |

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

One honest gap: cross-test 4 exercises the `V_Jplus1` branches at the **plain** tier only (z and
e flavours, nod and d); the 18 DC/GI/DC+GI raws' `V_Jplus1` branches have no runtime coverage
yet. They are the same mechanical transform of donors whose `V_Jplus1` branches the baseline
bank does cover, but given what the `V_Jplus1` audits keep finding, a DC/GI variant of
cross-test 4 is the natural next check if this family ever misbehaves.

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
  do the same for any new bank (AmbiguityAversion did: cross-test 4 covers its plain-tier
  V_Jplus1 branches from day one, though not yet its DC/GI tiers — see its section). So among
  the baseline preference mirrors, QH and EZ have full V_Jplus1 variant coverage (retrofitted)
  while AA has cross-test-4 coverage of the plain tier only. Still zero coverage across the
  whole ExpAsset family (re-verified 2026-08-31: the ExpAsset and ExpAssetU banks mention
  `V_Jplus1` in 0 subcodes; the ExpAssetz and ExpAssetze banks likewise, re-checked 2026-09-01).

  **The exposed-raw count was wrong and is now measured.** This entry previously said 660
  ExpAsset-family raws carry a `V_Jplus1` branch, itemised per family; those per-family figures
  were internally inconsistent (some were core-only, some were not) and the total was low. A
  direct `grep -l V_Jplus1` over every `*_raw.m` in the six families returns **1440 of 1440** —
  i.e. *every* raw carries one, with zero test coverage: ExperienceAsset 384, ExperienceAssetu
  384, ExperienceAssete 192, ExperienceAssetz 192, ExperienceAssetze 96, ExperienceAssetsemiz
  192. Those branches are age-shifted copies of the in-loop code — the shape that produced the
  `jj`/`N_j` bug. Largest remaining gap on any axis, by more than double what this doc claimed.

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
