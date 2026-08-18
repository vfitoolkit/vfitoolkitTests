# VFI Toolkit — test coverage

Coverage of the FHorz core test banks. Tier is assigned from each subcode's directory
(`Noa1_subcodes` / `With2A1_subcodes` / else), not its filename — filenames misattribute
cross-tests. Format in the tier columns is `variants + cross-testsx`.

Last updated: 2026-08-18

## CoreFHorzTests + the ExpAsset family + RiskyAsset

| bank | noa1 | withA1 | with2A | var | x-tests | panel vs dist | V_Jplus1 | core raws | QH | EZ |
|---|---|---|---|---|---|---|---|---|---|---|
| CoreFHorzTests | — | 16+6x | 16+4x | 32 | 10 | 32/32 | 32 | 128 | 42 | 24 |
| ExpAsset | 16+10x | 16+14x | 16+2x | 48 | 26 | 48/48 | 0 | 128 | 48 ⚠ | — |
| ExpAssetU | 16+10x | 16+14x | 16+2x | 48 | 26 | 48/48 | 0 | 128 | none | — |
| ExpAssete | 8+16x | 8+8x | 8+2x | 24 | 26 | 24/24 | 0 | 64 | 24 | — |
| ExpAssetz | 8+16x | 8+8x | 8+4x | 24 | 28 | 24/24 | 0 | 64 | 24 | — |
| ExpAssetze | 4+4x | 4+16x | 4+2x | 12 | 22 | 12/12 | 0 | 32 | 12 | — |
| ExpAssetsemiz | 8+4x | 8+6x | 8+2x | 24 | 12 | 24/24 | 0 | 64 | none | — |
| RiskyAsset | 16+0x | 16+25x | — | 32 | 25 | 32/32 | 0 | 80 | none | 32 |
| **total** | | | | **244** | **175** | **244/244** | 32 | | | |

419 subtests.

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

### Raw counts behind the "core raws" column

| family | total raws | QH | EZ | core |
|---|---|---|---|---|
| baseline (FHorz, excl. asset families) | 392 | 224 | 40 | 128 |
| ExperienceAsset | 128 | 0 | 0 | 128 |
| ExperienceAssetu | 128 | 0 | 0 | 128 |
| ExperienceAssete | 180 | 116 | 0 | 64 |
| ExperienceAssetz | 96 | 32 | 0 | 64 |
| ExperienceAssetze | 84 | 52 | 0 | 32 |
| ExperienceAssetsemiz | 64 | 0 | 0 | 64 |
| RiskyAsset | 122 | 0 | 42 | 80 |

Note the toolkit is inconsistent about where exotic-preference raws live: ExpAssete/z/ze
keep them in a `QuasiHyperbolic/` subdir of their own family, while the *baseline* QH raws
sit under `ExoticPrefs/QuasiHyperbolic/`.

### Known open items

- **V_Jplus1**: zero coverage across the whole ExpAsset family and RiskyAsset, against 660
  ExpAsset-family raws that all carry a `V_Jplus1` branch (ExperienceAsset 128,
  ExperienceAssetu 128, ExperienceAssete 160, ExperienceAssetz 96, ExperienceAssetze 84,
  ExperienceAssetsemiz 64). Those branches are age-shifted copies of the in-loop code — the
  shape that produced the `jj`/`N_j` bug. Largest remaining gap on any axis.
- ⚠ **ExpAsset QH bank tests code that does not exist.** `CoreFHorzQHExpAssetTests` has 48
  subtests, but there are **no QuasiHyperbolic raws for the plain `ExperienceAsset` family**
  (nor for `ExperienceAssetu` or `ExperienceAssetsemiz`). The bank has never been run — it
  has no `TestOutput`. Either it is test-first against unwritten solvers, or it is stale.
- **QH**: no mirror at all for ExpAssetU, ExpAssetsemiz or RiskyAsset. For ExpAssetsemiz
  this is a *toolkit* gap — the family has no QH raws — so closing it means solver code.
- **EZ**: mirrors exist only for the baseline and RiskyAsset. The entire ExpAsset family
  has no Epstein-Zin coverage, and the toolkit has no EZ raws for those families either.
- **QH banks run no panel simulation** at all.
- **RiskyAsset has no with2A tier**, and this is a *toolkit* gap: the family has zero `*2A*`
  raws. That is also the entire difference between ExperienceAsset's 128 core raws and
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
