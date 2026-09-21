# doPart for the test banks

> **Status 2026-09-22. ROLLED OUT to EVERY test bank: 39 banks, 227 parts.** Nothing is left out.
> The first three (`CoreFHorzTests` 10 parts, `CoreFHorzExpAssetTests` 19, `CoreFHorzExpAssetUTests`
> 14) were mapped by hand and GPU-run on 2026-09-21; the remaining 25 were mapped with a boundary
> proposer validated against those three hand maps, and all 28 pass the static checks. `CoreFHorzTPathExpAssetTests` and `CoreFHorzTPathTwoEndoTests` were
> initially left out as single-section banks and then converted anyway, at the user's request, so
> that the tiers they will gain later slot straight in.
>
> The seven CoreInfHorz banks followed on 2026-09-22 (27 parts). `CoreSummary` is wired into all 35
> and tested against all 31 existing diaries; the InfHorz banks needed four further check formats
> added to it, one of which had been making a whole bank read as zero checks.
>
> Decisions taken against the first draft: no try/catch (a test bank should die on an error); plain
> diary filename always, no doPart stamp.

## The short version

`DiscretizationMethodTests` already had a `doPart`, written the same way the replications write one.
Nothing needed inventing: the pattern was generalised from it, with one bank-specific rule (which
lines stay outside the `if`).

Granularity, as suggested: **one `doPart` entry per group of 8 variants, plus one per cross-test
block.** For `CoreFHorzTests` that is 10 entries; for `CoreFHorzExpAssetTests`, 19.

## Why

Four things the banks currently have no answer for:

1. **A bank that dies partway takes the rest of the run with it.** There is no try/catch anywhere,
   so an out-of-memory in one tier costs every tier after it, and `diary off` never fires.
   `CoreFHorzExpAssetTests` halts at fig 40 on a pre-existing OOM; figs 41-48 have never run, not
   because they are wrong but because they sit downstream of it. With `doPart` the remainder is one
   edit and one job.
2. **Splitting a bank across GPU jobs.** The big banks (ExpAsset 80 figures, ExpAssetU / QHExpAsset
   / RiskyAsset 48 each) are single all-or-nothing runs. Several tiers are currently listed as
   "blocked on GPU capacity" when what is actually needed is to run them separately.
3. **Iterating on one tier.** While a new tier is being written there is no reason to re-run the
   ones that already pass. This is exactly the reason `DiscretizationMethodTests` grew its `doPart`.
4. **Commenting out is the current mechanism, and it is worse.** Six main scripts carry
   commented-out `output=...` calls today:

   | file | commented | live |
   |---|---|---|
   | `CoreFHorzExpAssetTests.m` | 4 | 110 |
   | `CoreFHorzExpAssetUTests.m` | 1 | 73 |
   | `CoreInfHorzTests.m` | 2 | 12 |
   | `CoreInfHorzVFIAlgoTests.m` | 6 | 9 |
   | `CoreInfHorzTPathTests.m` | 6 | 4 |
   | `CoreInfHorzTPathAlgoTests.m` | 3 | 3 |

   Two different things are mixed in there: temporarily skipped (should be a `doPart`) and
   permanently disabled because it OOMs or the toolkit support does not exist (should be a `doPart`
   entry defaulting to `0` with the reason written next to it). Either way the state is invisible
   unless you read the whole file, it does not appear in the diary, and it gets committed.
   `CoreFHorzExpAssetTests` currently has an addpath whose comment reads "only needed while figs
   1-48 are commented out" — that is the failure mode.

## The convention

Copied from `DiscretizationMethodTests.m` and the replications. The diary block is untouched; the
only additions are the `doPart` block above it and one `fprintf` below it.

`doPart` is declared before the diary block, where it is the first thing in the file after the
header comment and so the easiest thing to edit. The `fprintf` goes immediately after `diary` opens,
so the run is labelled inside its own log.

```matlab
%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1):  nosemiz (figs 1-8)
% doPart(2):  nosemiz cross-tests
% doPart(3):  semiz (figs 9-16)
% ...
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so Fig 20
% is the same test whatever doPart says, and a png from a previous run is never overwritten by a
% different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzTestsdiary.txt
fprintf('CoreFHorzTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))
```

and each part:

```matlab
%% ===== doPart(1): nosemiz (figs 1-8) =====
if doPart(1)==1
    %% without d, without z, without e, without semiz
    figure_c=1;
    output=CoreFHorz_nod_noz_noe_nosemiz(n_d,n_a,...,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    ...
end % doPart(1): nosemiz (figs 1-8)
```

Three details worth fixing now rather than discovering later:

- **The banner line above each `if`.** `%% ===== doPart(1): nosemiz (figs 1-8) =====`. It means the
  `if` starts its own editor cell rather than hanging off the end of the previous section, and the
  MATLAB section navigator then lists the parts. The per-variant `%%` headings stay where they are,
  inside the block — which also means Run-Section on one of them still runs that single variant,
  ignoring `doPart` entirely, which is convenient rather than a problem.
- **The closing comment on every `end`.** `end % doPart(1): nosemiz (figs 1-8)`. With 19 parts and
  40-line bodies this is the difference between a readable file and an unreadable one.
  `DiscretizationMethodTests` already does it (`end % doPart(1), which is P0`).
- **`doPart` is hard-assigned in the file, not `if ~exist('doPart','var')`.** The `~exist` form lets
  you set it at the prompt without editing the file, but it also means a stale `doPart` left in the
  base workspace silently changes a later run — and these banks are run repeatedly in one MATLAB
  session. The replications all hard-assign; so should these.

## Granularity: the part map

One entry per 8-variant group, one per cross-test block. The `%%` section headings in the banks
already mark exactly these boundaries, so the map is read off the file rather than invented.

`CoreFHorzTests` (32 figures, 5 cross-test blocks) becomes 10 parts:

| part | contents |
|---|---|
| 1 | figs 1-8, nosemiz (8 variants) |
| 2 | nosemiz cross-tests (`CrossTests_nod`, `CrossTests_d`) |
| 3 | figs 9-16, semiz |
| 4 | semiz cross-tests (`CrossTests_nod1_semiz`, `_d1_semiz`) |
| 5 | semiz cross-tests 2 (semi-exo that is really a markov) |
| 6 | figs 17-24, with2A nosemiz |
| 7 | with2A nosemiz cross-tests |
| 8 | figs 25-32, with2A semiz |
| 9 | with2A semiz cross-tests |
| 10 | `TestFnsToEvaluate` |

Cross-tests get their own entry rather than being folded into the 8 above them. They are the
machine-precision checks, they are often the ones you want to re-run alone after a toolkit change,
and in at least one tier (ExpAsset 2A2) they are the whole tier. The cost is a slightly longer
vector.

`CoreFHorzExpAssetTests` (80 figures, 9 cross-test blocks) becomes 19, following its own headings:
noa1 nosemiz / noa1 nosemiz cross / noa1 semiz / noa1 semiz cross / a1 nosemiz / a1 nosemiz cross /
a1 nosemiz fake-expasset cross / a1 semiz / a1 semiz cross / a1 semiz cross2 / 2a1 nosemiz /
2a1 nosemiz cross / 2a1 semiz / 2A2 cross / 2A2 noa1 nosemiz / ... .

## The invariant that makes this safe

**Guarded blocks contain solves and nothing else. Setup, `addpath`, and the grid/parameter
preambles that sit under a `%% ===== SECTION =====` heading stay unguarded.**

This differs from `DiscretizationMethodTests`, where each part calls its own `DiscSetup_Pn` inside
the `if`. The Core banks share one setup script across every tier, so per-part setup would mean
repeating the same `CoreFHorz_setup` call 19 times. Leaving the preambles unguarded costs
microseconds and achieves the same independence — provided the preambles are genuinely additive.

That gives one mechanical audit per bank, and it is the only thing that can go wrong: **for each
proposed part, every variable it reads must be assigned either in an unguarded preamble or inside
that part.** Two checks done already:

- `CoreFHorzTests`: the with2A preamble (`n_a_2A`, `n_a_2A_big`, `a2_grid_2A`, `a_grid_2A`,
  `a_grid_2A_big`, `Params.phi1`, `Params.phi2`) defines new names only — it overwrites nothing
  that parts 1-5 use — so running it unconditionally is safe. Every `n_a_notsobig` /
  `a1_grid_notsobig` / `a_grid_notsobig` is redefined immediately above the call that uses it
  (lines 234, 249, 264, 309, 322, 330, 338), so no figure depends on a definition in another part.
- `CoreFHorzExpAssetTests`: the three section preambles each re-call `CoreFHorzExpAsset_setup` plus
  their own `addpath`s. Left unguarded they execute in the same order as now, so behaviour is
  unchanged; the `2a1`/`2A2` grid definitions likewise define new names.

Where a bank fails this audit the fix is to move the definition up into the unguarded preamble, not
to add a guard to it.

## Figure numbers

No change. Every Core bank already writes `figure_c` as a literal, which is the right thing and
better than the DMT scheme (a running counter with ten numbers reserved per part). Keep it:
**never convert `figure_c` to a running counter, and never renumber to close a gap.** A partial run
must produce pngs with the same numbers as a full run, or the saved figures stop being comparable
across runs. Gaps are already normal and intended (`CoreInfHorzTPathTests` has only figs 3 and 4
live by design).

## Diary: unchanged

The diary block keeps the plain filename it already has, on a partial run as much as a full one
(`CoreFHorzTestsdiary.txt`, etc.), exactly as `DiscretizationMethodTests` does. A doPart-stamped
filename for partial runs was considered and rejected: one predictable path per bank is worth more
than never overwriting a partial run's log.

The consequence to be aware of: **a partial run overwrites the full run's diary.** Copy a diary you
want to keep before running a subset. The new header line says which parts a diary came from:

```
CoreFHorzTests, run started 20-Sep-2026 14:02:11, doPart=[1111111111]
```

Because the filename does not change, `diary` keeps its command syntax and the block is untouched
apart from that one added `fprintf`.

## Migrating the currently commented-out blocks

The committed default vector must reproduce **current** behaviour, not "everything on":

- a section that currently runs → `1`
- a section that is currently commented out in full → `0`, with the reason written into the legend
  at the top (`% doPart(12): 2a1 nosemiz — OFF: out-of-memory on the 151-point a1 grid`), and the
  code uncommented so it stays live and greppable
- a call commented out inside an otherwise-live section → leave it commented; `doPart` is for
  sections, not individual calls

`DiscretizationMethodTests` already ships `doPart(9)=0` with a paragraph explaining why, which is
the model to follow.

## Rollout and validation

Pure control flow, so almost all of the validation needs no GPU. Every bank passed all of this:

1. **Insertion-only diffs.** Wrapping a group re-indents its body, so compare with `diff -w` and
   require *zero* deleted or changed lines. 28/28 clean: every change is an added banner, `if`,
   `end`, `fprintf`, legend, addpath or CoreSummary line.
2. **checkcode message classes unchanged.** Not merely "clean" - the set of distinct messages must
   not grow against the pre-edit file. 28/28 unchanged (the only message anywhere is "Value assigned
   to variable might be unused"; its *count* falls, because MATLAB can no longer prove the last
   `output=` in a part is dead once it sits inside a conditional).
3. **Every part self-contained.** Per part, collect identifiers read before assignment and require
   each to come from the unguarded preamble, the setup script, a function file on the addpath, or a
   builtin; and separately require that no `addpath` or `*_setup` call sits inside a part. 205/205
   parts pass. The checker was mutation-tested: moving the `with2A` grid preamble inside part 6 of
   `CoreFHorzTests` makes it report exactly `doPart(7) reads but nothing supplies: a_grid_2A,
   a_grid_2A_big, n_a_2A, n_a_2A_big`, and the same for parts 8 and 9.
4. **Nothing outside a part.** Every `figure_c`, `output=` and `exportgraphics` line must sit inside
   some part. 28/28 exact.
5. **Labels consistent.** The label must be identical in all four places it appears (the `%%`
   banner, the printed `fprintf`, the closing `end` comment, the legend). 28/28, the only difference
   being the deliberate "WILL ERROR" note on `CoreFHorzExpAssetUTests` part 14.
6. **GPU: the default-vector run.** Done for the first three on 2026-09-21 - see below. Outstanding
   for the other 25.

### The boundary proposer, and why it was trusted

Mapping 25 banks by hand was not realistic, so the boundaries were proposed mechanically and the
proposer was validated by **reproducing the three hand maps exactly** (10, 19 and 14 parts, matching
starts, full coverage). Building it that way paid for itself twice:

- Its first version ended a part at the last `output=` line, leaving the following `exportgraphics`
  **outside** the part - which would have saved a png for a figure that was never drawn. The
  comparison against the hand map caught it immediately.
- It initially merged the 32-figure 2A2 region into one part and collapsed `CoreFHorzPTypeTests` to
  a single part. Both are now handled by the rule that a figure run breaks at a heading once it
  holds 8 figures, or if it holds none.

The rules it encodes are the ones in this document: a part is a maximal run of test calls; `addpath`
and `*_setup` break a part and stay unguarded; a figure group and a cross-test block are separate
parts; a part starts at the `%%` heading above its first call and ends at its last call.

### The default vector is not always all-ones

`CoreFHorzExpAssetTests` now ships `doPart=[0,...,0,1,1,1]` - parts 1-16 off, so a run goes straight
to the 2A2 tiers being worked on. That was set by the user on 2026-09-22, and it is the mechanism
being used as intended rather than a deviation: the committed vector records what a run of this bank
should currently do, and the legend carries the reason and a `TO RESTORE:` line with the all-ones
vector. `CoreInfHorzEntryExitTests` ships `[0,1,1]` for the same kind of reason. Every other bank is
all-ones.

### What the GPU runs of 2026-09-21 showed

| bank | parts run | verdict |
|---|---|---|
| `CoreFHorzTests` | 10 of 10 | **1602 checks, 0 failed, ALL CHECKS PASSED** |
| `CoreFHorzExpAssetTests` | 18 of 19 | died in part 18; no verdict line |
| `CoreFHorzExpAssetUTests` | 14 of 14 banners | died inside part 14; no verdict line |

`CoreFHorzTests` is the clean result: 1602 checks, exactly the pre-edit baseline count, so the wrap
changed nothing.

The other two died in tiers already documented as unimplemented, not in anything the doPart edit
touched. `CoreFHorzExpAssetTests` stopped at the first fig-65 solve with `Arrays have incompatible
sizes of 202 and 303` in `ValueFnIter_FHorz_ExpAsset_DC1_nod1_noz_raw` - the 2A2 with-a1 tier, which
that bank's own header says is written test-first and waiting on the toolkit port; parts 16 and 17
(the noa1 2A2 tiers) ran clean. `CoreFHorzExpAssetUTests` died inside part 14, the tier whose legend
entry already said it would.

**The consequence worth knowing: a bank that dies never reaches `CoreSummary`, so it produces no
verdict line at all.** With no try/catch that is the intended behaviour, and `doPart` is the remedy:
set the known-broken tier to zero and the rest of the bank runs to the end and reports. For these
two that is `doPart(18:19)=0` and `doPart(14)=0` respectively.

Order for anything added later: map with the proposer, review its boundaries, apply, then run checks
1-5 before committing. Banks to skip remain those with a single natural section.

## Part counts, as built

28 CoreFHorz banks, 178 parts (plus 27 across the 7 CoreInfHorz banks below: 205 in all). Baseline check counts are from `CoreSummary` on each bank's current
diary; compare them after any future run, because a falling count means checks stopped executing and
nothing else in the output would say so.

| bank | parts | figs | checks in its current diary |
|---|---|---|---|
| CoreFHorzExpAssetTests | 19 | 80 | 1799 (run died in part 18) |

| CoreFHorzExpAssetUTests | 14 | 48 | 1394 (run died in part 14) |
| CoreFHorzExpAssetzTests | 12 | 24 | 985 |
| CoreFHorzExpAsseteTests | 11 | 24 | 985 |
| CoreFHorzExpAssetzeTests | 11 | 12 | 577 |
| CoreFHorzTests | 10 | 32 | **1602, 0 failed** |
| CoreFHorzQHTests | 9 | 32 | 4464 |
| CoreFHorzGulPesendorferTests | 8 | 32 | 1034 |
| CoreFHorzPTypeTests | 7 | 0 | (no diary yet) |
| CoreFHorzEZTests | 7 | 16 | 4306 |
| CoreFHorzExpAssetsemizTests | 6 | 24 | 919 |
| CoreFHorzQHExpAssetTests | 6 | 48 | 3809 |
| CoreFHorzQHExpAssetUTests | 6 | 48 | 3809 |
| CoreFHorzRiskyAssetTests | 6 | 48 | 1411 |
| CoreFHorzTPathRiskyAssetTests | 5 | 16 | (no diary yet) |
| CoreFHorzTPathTests | 5 | 16 | 492 |
| CoreFHorzResidAssetTests | 4 | 16 | (no diary yet) |
| CoreFHorzRiskyAssetEZTests | 4 | 32 | 7233 |
| CoreFHorzAmbiguityTests | 4 | 12 | 428 |
| CoreFHorzQHExpAsseteTests | 4 | 24 | 3144 |
| CoreFHorzGPExpAssetTests | 3 | 16 | 506 |
| CoreFHorzQHExpAssetsemizTests | 3 | 24 | 3144 |
| CoreFHorzQHExpAssetzTests | 3 | 24 | 3032 |
| CoreFHorzTPathPTypeTests | 3 | 2 | (no diary yet) |
| CoreFHorzQHExpAssetzeTests | 2 | 12 | 1800 |
| CoreFHorzRiskyAssetAmbiguityTests | 2 | 1 | 220 |
| CoreFHorzTPathExpAssetzTests | 2 | 4 | 166 |
| CoreFHorzTPathQHTests | 2 | 8 | (no diary yet) |

Not converted, single section each: `CoreFHorzTPathExpAssetTests`, `CoreFHorzTPathTwoEndoTests`.

### CoreInfHorz banks (added 2026-09-22)

| bank | parts | figs | checks in its current diary |
|---|---|---|---|
| CoreInfHorzVFIAlgoTests | 7 | 0 | 70 (was reading as **0** before CoreSummary learned its format) |
| CoreInfHorzTests | 5 | 8 | 95 |
| CoreInfHorzQHTests | 5 | 4 | 190, 8 flagged |
| CoreInfHorzInheritAssetTests | 3 | 2 | 30 |
| CoreInfHorzEntryExitTests | 3 | 5 | 3 |
| CoreInfHorzTPathTests | 2 | 4 | 189, 8 flagged |
| CoreInfHorzTPathAlgoTests | 2 | 2 | 11 |

Two of these needed a hand override rather than the proposer's suggestion:

- **`CoreInfHorzTPathTests`** was proposed as one part. It is two tiers - one endogenous state
  (figs 3-4 live, 1-2 and 5-8 commented out) and two endogenous states (figs 9-10) - but the tier
  holds fewer than 8 figures, so the 8-figure break never fired. Split by hand.
- **`CoreInfHorzEntryExitTests` already had a hand-rolled doPart**: `dofigs1to4=0; if dofigs1to4==1
  ... end`, whose comment gives the same reason this document does ("A flag rather than
  commented-out lines, so restoring is one character and cannot be half-done"). That flag was folded
  into `doPart(1)`, which **defaults to 0** to preserve the existing behaviour, with the original
  explanation moved verbatim into the legend. This is the one file in the whole rollout whose diff
  is not insertion-only; the nine removed lines are exactly the flag block and nothing else.

### CoreStationaryGeneralEqm (added 2026-09-22)

Nine parts: the `constrainpositivemethod` round trips, then each of the four families (InfHorz,
FHorz, InfHorz PType, FHorz PType) split into its fminalgo-agreement and its
parameter-constraint-invariance exercise.

| part | runs | cost |
|---|---|---|
| 1 | constrainpositivemethod round trips | no model, no GE solve - seconds |
| 2 / 3 | InfHorz fminalgo / constraint invariance | 6 / 7 GE solves |
| 4 / 5 | FHorz fminalgo / constraint invariance | 6 / 7 GE solves |
| 6 / 7 | InfHorz PType fminalgo / constraints | 6 / 7 GE solves |
| 8 / 9 | FHorz PType fminalgo / constraints | 6 / 7 GE solves |

52 GE solves in all. The proposer is blind to this bank - it names its outputs `output1`,
`output2`, `output1ptype` rather than `output` - so the map was made by hand. Independence was
checked directly: no section takes another's output as an argument, and the transform block's
`cpm_*` variables are never read after it.

Two things make it the best case for `doPart` in the whole repo. The cost is wildly concentrated:
the diary is **620,430 lines** and roughly 518,000 of them come from the InfHorz fminalgo section
alone (the stochastic CMA-ES `fminalgo=4` run), so `doPart(2)=0` removes about 83% of the output.
And its last run **did not finish** - the final banner is InfHorz PType fminalgo, then 16,639 lines
ending mid-solve with no terminating banner, so parts 8 and 9 never reported. That matches the
FHorz PType crash that has since been fixed; `doPart=[0,0,0,0,0,0,0,1,1]` now runs just those two.

Nothing is now without a `doPart`. The only two `.m` files at bank level that have none are
`CoreFHorzTests/TestFnsToEvaluate.m`, which is a subcode called by `CoreFHorzTests` and is already
its `doPart(10)`, and `DiscretizationMethodTests/DiscSummary.m`, which is that bank's summariser
function rather than a test script.

## CoreSummary: the verdict line (BUILT 2026-09-20)

`SharedSubcodes/CoreSummary.m`, the Core-bank counterpart of `DiscSummary`. Each of the three banks
now ends with

```matlab
%% One verdict for the whole run
CoreSummary('./TestOutput/CoreFHorzTestsdiary.txt')

diary off
```

and carries `addpath('../SharedSubcodes/')` with its other addpaths. It closes the diary, reads it
back, and appends one verdict to it.

### Why the check forms cannot share one bar

The survey that came first mattered more than the code. Grepping every `fprintf` turns up a dominant
form and six others, several of which must not be held to the same bar:

| form | what it is | bar |
|---|---|---|
| `this should be zero: 2.665e-15` | identity — two ways of computing the same object | `restol` |
| `should give zero: 0.000e+00` | identity too (QH `beta0=1`, EZ `riskaversion=0`) | `restol` |
| `this should be close to zero: 4.839e-03` | **not** an identity | `closetol` |
| `V should be ~0: 0.00000001` | identity, InfHorz VFI-algo spelling | `restol` |
| `Pol should be  0: 0.00000000` | identity (exact Policy match) | `restol` |
| `this should be one: 1` | a flag | fails at 0 |
| `should be very similar (small): 0.00008` | **not** an identity | `closetol` |
| `should be WELL ABOVE zero: 0.4286` | inverted | fails NEAR zero |
| `should NOT be zero: 1.119e-06` | inverted | fails NEAR zero |
| `should be near zero: 3.49e-06` | its own population | `neartol` = 1e-3 |
| `should be near zero (loose: <reason>)` | the bank saying so itself | `closetol` |
| `should be small, r: 0.0727, w: 0.3281` | loose, and **two** values | `closetol` |

The third one is the catch. In the ExpAssetU diary every one of the 23 values above 1e-6 is
`StationaryDist with/without grid interp`, running up to 2.5e-2 — the agent distribution genuinely
differs with the interpolation layer, because interpolation puts mass between grid points. A single
bar would have reported 23 failures on a bank that is fine, which is precisely the kind of
false-confidence-destroying noise that would make the summariser get ignored. The bank's own wording
("close to") already drew the distinction; the summariser just has to respect it.

Two more forms the parser has to survive: a line may carry **several** checks
(`Cross test (noa1): this should be zero: V 0.000e+00, Policy 0.000e+00, Dist 6.939e-18` is three,
and is counted as three), and anything in parentheses after the value is diagnostic detail, not a
check (`(rel 1.3e-16, max|V|=1.1e+07, worst at [3] of [5])`, `(2 of 100 entries differ)`), so the
line is cut at the first `(` after the colon.

### Where restol comes from

Not chosen by feel. Pooling the identity residuals in the three current diaries — 3031 checks — the
distribution is sharply bimodal with an empty decade in it:

| band | count |
|---|---|
| exactly zero | 2524 |
| below 1e-12 (ordinary double noise) | 468 |
| 1e-12 to 1e-10 | **0** |
| 1e-10 to 1e-8 (the V~1e7 ULP band: `ValueFnFromPolicy` with grid interp) | 39 |
| 1e-8 and above | **0** |

The largest residual anywhere is 5.588e-09. So `restol=1e-7` sits in empty space, a factor of 18
above anything these banks have printed and five orders below the `close to zero` scale. A first
draft used 1e-8, which the data then showed was only a factor of 1.8 above that largest residual —
too tight to survive ordinary drift. Re-derive the bar if a bank starts printing into the gap; do
not quietly raise it to make a run pass.

Passing identities are reported in three bands (exactly zero / below 1e-10 / in the ULP band) so
that drift *towards* the bar is visible long before anything fails. A failing identity is counted as
a failure and in none of the bands.

### Per part

Each part now prints its own banner as its first statement:

```matlab
if doPart(1)==1
    fprintf('\n===== doPart(1): nosemiz (figs 1-8) =====\n')
```

so the diary says what is running, and `CoreSummary` attributes counts and failures to the part they
sit in. A partial run's summary is therefore self-describing — it lists only the parts that ran. A
bank without banners still works; everything is attributed to `(whole run)`.

### What it prints

```
1602 checks, 0 failed.
  identities ("should be zero", bar 1e-07): 1586
    passing: 1166 exactly zero, 420 below 1e-10 (double noise), 0 in the V~1e7 ULP band above it
    largest: 2.842e-14  (diary line 9606)
      ValueFnFromPolicy, this should be zero: 2.842e-14
  "close to zero" checks (bar 0.1, NOT identities): 16
    largest: 5.000e-02  (diary line 8)

Checks by part (only the parts that actually ran appear here):
  doPart(1): nosemiz (figs 1-8)                         ...
```

Failures are listed with diary line number and part.

### Three parser bugs, all found by running it on diaries it was not written against

It was built against 3 diaries and then run against all 31. Each round found a defect, and all three
are the same mistake in different clothes - **reading a number that is not the check's value**:

1. **Digits inside identifiers.** `CoreFHorzRiskyAssetTests` reported 12 failures on lines like
   `... should be zero: V 0.000e+00, Policy 0.000e+00, Dist 0.000e+00, AllStats.a2.Mean 0.000e+00` -
   every value an exact zero. The culprit was the **`2` in `AllStats.a2.Mean`**. Values are now
   matched only where a number stands on its own. That bank went 12 failures to 0, and its count
   corrected from 1423 to 1411.
2. **A tolerance quoted in an annotation.** `CoreInfHorzTPathTests` prints
   `this should be zero, AgentDist (note: tolerance=1e-6): 4.907e-06`. The first colon after the
   wording is the one inside the annotation, so the `1e-6` was read as a second value and the line
   was counted - and failed - twice. When the text after that colon still contains a colon, the
   value is now taken from after the last one. That bank went from 197 checks/14 failures to
   189/8.
3. **Whole formats unseen.** `CoreInfHorzVFIAlgoTests` reported **0 checks across 1929 lines**,
   because it writes `V should be ~0: 0.00000001` and `Pol should be  0: 0.00000000` rather than
   "should be zero". Four further forms were added: those two, `this should be one:` (a flag,
   which fails at 0), `should be very similar (small):` (loose), and `should be WELL ABOVE zero:`
   (inverted - it fails by being NEAR zero). Also, values are now allowed to be labelled after the
   wording (`this should be zero, Policy: 6.000e+00`), which is what made a Policy difference of 6
   on a constant transition path visible for the first time.

A summariser that silently reports "0 checks" is worse than no summariser, so form coverage is worth
checking whenever it meets a bank it has not seen. The additions are strictly additive: every
CoreFHorz baseline is unchanged (`CoreFHorzTests` still 1602, ExpAsset 1799, RiskyAsset 1411), no
bank lost checks, and no new failures appeared in the CoreFHorz set.

### neartol, and a bar that must not be loosened

`should be near zero` is CoreStationaryGeneralEqm's wording for comparing two GE solves, and its
values sit on neither the identity nor the distribution scale, so it gets its own bar. That bank's
plain `near zero` values run 1.1e-16 to 7.6e-5 and then jump to two outliers at 6.5e-2 and 6.9e-2 -
a gap of a factor of 854. `neartol=1e-3` sits in it, 13x above the highest ordinary value and 65x
below the outliers.

The outliers matter: they are `fminalgo 1 vs 8`, two optimizers disagreeing about the same
equilibrium by about 7%, which is the known "lsqnonlin wrong on InfHorz" problem. A bar loose enough
to pass them would hide exactly the thing worth seeing. Note the bank marks its own known-loose
checks `(loose: <reason>)` and those are taken at `closetol` instead - it is saying so itself.

### What it flags in the CoreInfHorz banks

`CoreInfHorzQHTests` (8) and `CoreInfHorzTPathTests` (8). On inspection these split three ways:

- **Characterisations written as identities.** Five of the QH ones are
  `maxaprimediff=125 (the WHOLE grid) vs the default, Policy, this should be zero: 4.900e+01`, and
  three more are `howards on vs off ... should be zero to tolerance IF it converged`. Comparing a
  deliberately different setting, and a conditional, are not identities - but they are worded
  exactly like one, so nothing downstream can tell them apart. If they should stop being flagged,
  the fix is in the bank's wording (as `close to zero`, or an explicit "non-zero by design" as
  CoreInfHorzVFIAlgoTests already uses for its diagnostic scans).
- **Failures against the bank's own stated tolerance.** `Constant TPath ... AgentDist (note:
  tolerance=1e-6): 4.907e-06` exceeds the 1e-6 the line itself names.
- **A Policy difference of 6 on a constant transition path**, twice
  (`Constant TPath (with DC and GI), this should be zero, Policy: 6.000e+00`). This one had never
  been visible to any reader who scanned for `should be zero:`, and it sits squarely in the area
  already flagged as "TPath banks are silently wrong on a path, not errors".

None of these was caused by this work; they are pre-existing diary contents that the summariser has
made visible.

### How it was tested

- **On all three real diaries** (`CoreFHorzTestsdiary.txt`, and the ExpAsset / ExpAssetU ones):
  1602, 183 and 1285 checks, 0 failed. These are real data, so the parser is exercised against
  every format the banks actually emit, not against a guess at them.
- **On a synthetic diary** carrying part banners, both bar classes, a multi-value line, both
  parenthetical forms, an eyeball-only line and a stray data line, plus three deliberate failures
  (an identity at 1.234e-03, a Policy mismatch at 1.0, a loose check at 0.5). All three are caught
  and attributed to the right part; the eyeball line and the data line are correctly ignored; the
  multi-value line counts as 3.
- **Called twice on the same diary**: identical output, so the "stop at a summary already there"
  guard works and re-running it on a finished diary does not double-count the echoed failures.
- **`checkcode` clean.**

The check counts are the thing to watch across runs. 1602 / 183 / 1285 are the numbers to compare
the first GPU runs against — a smaller count means checks stopped executing, which nothing else in
the output would tell you.
