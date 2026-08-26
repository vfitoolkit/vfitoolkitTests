# vNM RiskyAsset batch-corrections proposal

Proposal for fixing the 10 bugs recorded (list-don't-fix policy) in the vNM (non-Epstein-Zin)
riskyasset files during the EZ-riskyasset build. Source list: section
'### vNM RiskyAsset bug list' of `EZRiskyAsset_coverage_proposal.md` (same directory).

Every item was re-verified against the toolkit on disk on 2026-08-19
(repo `/home/kmarshallbanana/Dropbox/sandpit/Claude/VFIToolkit-matlab`). All 10 bugs still exist
exactly as recorded; every recorded line number is still accurate (no drift). Severity, trigger,
and bank-impact analysis below use the CURRENT vNM CoreFHorzRiskyAssetTests bank
(`.../vfitoolkitTests/CoreFHorzRiskyAssetTests/`), which:

- never passes `vfoptions.V_Jplus1` (0 mentions in the whole bank),
- tests lowmemory only at the values each shape's ladder covers (nosemiz: 0/1; semiz z legs:
  0/1/2(/3 for z+e); semiz noz legs: 0/1 (noe) and 0/1/2 (e)),
- has no infeasible states: the ReturnFns give `c = pension (or w*kappa_j*z) + (1+r_a1)*a1 + a2
  - a1prime - savings`, and `pension>0`, so the zero-choice corner is always feasible and V is
  finite everywhere (no -Inf anywhere in any V the bank produces),
- at terminal age the optimum is always the zero-savings/zero-a1prime corner (utility strictly
  increasing in c, c strictly decreasing in savings and a1prime, no bequest motive).

Those four facts are what keep the current bank green despite these bugs. The last GPU run's
diary (`TestOutput/CoreFHorzRiskyAssetTestsdiary.txt`) shows every "should be zero" check at
0.00000000, consistent with the analysis below.

Severity classes: CRASH (errors when the branch runs) / SILENT-WRONG (wrong numbers, no error) /
SILENT-ZERO (V left at zeros, no error) / COSMETIC.

---

## Summary table

| # | File(s) | Class | Bank exercises the branch? | Fix changes live bank results? |
|---|---------|-------|---------------------------|-------------------------------|
| 1 | `ValueFnFromPolicy/RiskyAsset/ValueFnFromPolicy_FHorz_RiskyAsset_GI.m` (4 sites + dead code) | SILENT-WRONG (latent NaN) | Path yes (every withA1 GI FromPolicy leg); trigger no (no -Inf V) | No |
| 2 | `ValueFnIter/FHorz/RiskyAsset/Raw/`: `_nod1_raw` (2 blocks), `_nod1_noz_raw` (1), `_raw` (2), `_noz_raw` (1) | SILENT-WRONG (escalates to downstream out-of-range crash when it bites) | YES — default terminal branch of every basic + lowmemory withA1 leg | **No** (masked by the terminal corner solution — see item 2) |
| 3 | Same dir: `_nod1_noz_e_raw:119`, `_nod1_noa1_noz_e_raw:106`, `_noa1_noz_e_raw:110`, `_noz_e_raw:132` | CRASH | No (needs V_Jplus1 + lowmemory=1) | No |
| 4 | Same dir: `_nod1_noz_raw:53`, `_noz_raw:53` | CRASH | No (needs V_Jplus1) | No |
| 5 | `.../RiskyAssetSemiExo/`: `_nod1_noa1_raw`, `_noa1_e_raw`, `_nod1_noa1_e_raw`, `_noa1_noz_e_raw`, `_nod1_noa1_noz_e_raw` | SILENT-ZERO | No (needs V_Jplus1 + lowmemory≥2) | No |
| 6 | `.../RiskyAssetSemiExo/ValueFnIter_FHorz_RiskyAssetSemiExo_nod1_raw.m:97` | CRASH | No (needs V_Jplus1) | No |
| 7 | Same file (2 blocks) + `..._RiskyAssetSemiExo_raw.m` (2 blocks) | CRASH | No (needs V_Jplus1) | No |
| 8 | `ValueFnFromPolicy/RiskyAsset/ValueFnFromPolicy_FHorz_RiskyAsset_SemiExo_GI.m` (1 site + cosmetics) | SILENT-WRONG (latent NaN) | Path yes (semiz+e GI FromPolicy legs); trigger no | No |
| 9 | `.../RiskyAssetSemiExo/GridInterpLayer/`: `GI1_noz_raw` (3 ladders), `GI1_nod1_noz_raw` (3), `GI1_noz_e_raw` (1), `GI1_nod1_noz_e_raw` (1) + lint | SILENT-ZERO | No (bank only passes ladder-covered lowmemory values) | No |
| 10 | `ValueFnFromPolicy/RiskyAsset/ValueFnFromPolicy_FHorz_RiskyAsset_SemiExo.m:215` (+ cosmetics) | SILENT-WRONG (latent NaN) | Path yes (every semiz+e FromPolicy leg); trigger no | No |

Bottom line: **no item changes any result the current vNM bank produces.** Items 3/4/6/7 are
crashes in V_Jplus1-only branches; 5/9 are silent zeros in lowmemory values the bank never
passes; 1/8/10 need -Inf next-period values the bank's calibration never generates; item 2 is
live code but exactly masked by the bank's terminal corner solution. A clean re-run of the bank
after applying the whole batch should reproduce the current diary bit-for-bit (up to the known
SimPanel parfor non-reproducibility).

---

## Suggested application order

Crashes and silent-wrongs first; cosmetics only bundled into files the batch already touches.

1. **Items 6 + 7 together** — both live in the terminal-V_Jplus1 branch of the same two
   RiskyAssetSemiExo files (`nod1_raw`, `raw`); item 6 (line 97) crashes before item 7's lines
   are reached in `nod1_raw`, so fixing one without the other leaves the branch broken.
2. **Item 3** — four one-line argument fixes (CRASH class).
3. **Item 4** — two one-line reshape fixes (CRASH class).
4. **Item 2** — the only live-code item; 6 blocks in 4 files, two-token change per block.
   After this, re-run the withA1 nosemiz figs and confirm the diary zeros are unchanged
   (they should be — see the masking analysis).
5. **Items 1, 8, 10** — the three FromPolicy isnan clears (SILENT-WRONG class, one insertion
   pattern each).
6. **Item 5** — five one-token ladder-guard fixes (`==1` → `>=1`) (SILENT-ZERO class).
7. **Item 9** — eight one-token ladder-guard fixes across 4 GI files (SILENT-ZERO class).
8. **Cosmetics, bundled per touched file only**: item 1's dead prealloc (same file as its fix),
   item 2's `ceil(dindex)` no-op (same file as its fix), items 8/10's dead if/else + unused
   locals (same files as their fixes), item 9's unused `a2ind` + stale `%#ok<NASGU>` on
   `special_n_d4` (only in the 4 GI noz files the batch touches; the 4 with-z GI1 siblings have
   the same lint but are NOT otherwise touched, so leave them).

---

## Item 1 — GI FromPolicy: missing isnan clear after the per-u corner-interpolation sums

**File**: `ValueFnFromPolicy/RiskyAsset/ValueFnFromPolicy_FHorz_RiskyAsset_GI.m`

**Verified locations** (all as recorded, no drift): the four `per_u` weighted sums at lines
226 (noz-noe), 241 (noz-e), 257 (z-noe), 274 (z-e), whose `EVnext_atpolicy` results (lines
227/242/258/275) feed straight into `V(...,jj)` with no isnan clear. The e-collapse at line 204
also has no clear (matching the plain file's line 209, which likewise has none — the plain
file's convention is that the final per-u clear catches everything).

**Severity**: SILENT-WRONG. Under GI, `w_a1_upper=(L2-1)/(n2short+1)` is exactly 0 whenever
L2==1 (policy sits on a coarse grid point — extremely common). If the zero-weight corner
(`a1_upper`, or an `a2` corner with zero `a2primeProbs` weight) holds `-Inf` (infeasible next
state), `0*(-Inf)=NaN` enters `per_u`, survives the pi_u sum, and propagates into V for ALL
earlier ages.

**Trigger**: gridinterplayer FromPolicy (withA1, any shock shape) on a model whose V' contains
-Inf (any calibration with genuinely infeasible states). **Bank**: exercises the code path in
every withA1 nosemiz leg ("ValueFnFromPolicy with grid interp" checks), but its calibration has
no -Inf states, so the NaN never fires. Diary shows zeros.

**Fix** (insert one line after each of the four per_u sums; pattern = the plain
`ValueFnFromPolicy_FHorz_RiskyAsset.m` lines 237/254, and the EZ GI file's
`EVnextOfPolicy(isnan(...))=0` at its lines 264/286/310/336). Representative site (lines
226-228, the noz-noe block):

BEFORE
```matlab
            per_u=wa1l.*wa2l.*EV_LL + wa1l.*wa2u.*EV_LU + wa1u.*wa2l.*EV_UL + wa1u.*wa2u.*EV_UU;
            EVnext_atpolicy=sum(per_u .* shiftdim(pi_u,-1), 2);
            V(:,jj)=F_jj+beta*EVnext_atpolicy;
```

AFTER
```matlab
            per_u=wa1l.*wa2l.*EV_LL + wa1l.*wa2u.*EV_LU + wa1u.*wa2l.*EV_UL + wa1u.*wa2u.*EV_UU;
            EVnext_atpolicy=sum(per_u .* shiftdim(pi_u,-1), 2);
            EVnext_atpolicy(isnan(EVnext_atpolicy))=0; % zero corner weights times -Inf next-states
            V(:,jj)=F_jj+beta*EVnext_atpolicy;
```

Sites (textually identical apart from the pi_u shiftdim offset and sum dim, and the V
assignment indexing): after line 227 (`shiftdim(pi_u,-1)`, sum dim 2), after line 242
(`-2`, dim 3), after line 258 (`-2`, dim 3), after line 275 (`-3`, dim 4). The line-204
e-collapse needs no separate clear once these are in (any NaN it produces is caught by the
per-u clear, same as in the plain file).

**Bundled cosmetic** (same file): delete the dead preallocation at lines 122-126
(`if N_z==0 && N_e==0, a1_lower=ones(N_a,N_j,'gpuArray'); else ... end`) — line 132
(`a1_lower=a1_mid;`) unconditionally overwrites it.

**Behaviour impact**: none on the current bank (no -Inf V, so the inserted clears never touch a
value); changes results only for calibrations with infeasible states, where current output is
NaN-poisoned.

---

## Item 2 — Base raws: terminal-period Policy decode uses the wrong N_d (LIVE default branch)

**Files/blocks** (all verified at the recorded lines):

| File | Blocks | First dim of ReturnMatrix | Wrong divisor used | Correct divisor |
|------|--------|--------------------------|--------------------|-----------------|
| `Raw/ValueFnIter_FHorz_RiskyAsset_nod1_raw.m` | lines 52-55 (lm0), 65-68 (lm1) | `N_d3*N_a1` | `N_d` (=`N_d2*N_d3`, file line 12) | `N_d3` |
| `Raw/ValueFnIter_FHorz_RiskyAsset_nod1_noz_raw.m` | lines 45-48 | `N_d3*N_a1` | `N_d` (=`N_d2*N_d3`, line 11) | `N_d3` |
| `Raw/ValueFnIter_FHorz_RiskyAsset_raw.m` | lines 53-57 (lm0), 67-71 (lm1) | `N_d1*N_d3*N_a1` | `N_d` (=`N_d1*N_d2*N_d3`, line 14) | `N_d1*N_d3` |
| `Raw/ValueFnIter_FHorz_RiskyAsset_noz_raw.m` | lines 46-50 | `N_d1*N_d3*N_a1` | `N_d` (=`N_d1*N_d2*N_d3`, line 13) | `N_d1*N_d3` |

The no-V_Jplus1 terminal branch maxes a ReturnMatrix whose first dimension EXCLUDES d2 (d2 is
aprimeFn-only and there is no continuation value at terminal), but decodes `maxindex` with the
full `N_d` that INCLUDES d2. Known-good patterns, both verified: `noz_e_raw` lines 52-56 (uses
`N_d1*N_d3` throughout), and each nod1 file's own V_Jplus1/main-loop decode
(`rem(maxindex-1,N_d3)+1` / `ceil(maxindex/N_d3)`, e.g. `nod1_raw` lines 123-125, 213-215).

**Severity**: SILENT-WRONG when it bites — and worse: since `rem(maxindex-1,N_d2*N_d3)+1` can
reach `N_d2*N_d3 > N_d3`, the decoded d3 row can be OUT OF RANGE for `d3_grid`, so downstream
commands (PolicyInd2Val, StationaryDist, FromPolicy) would then index garbage or error.

**Trigger**: any withA1 riskyasset solve with `N_d2>1` whose terminal-age optimum is NOT in the
first `N_d3` (nod1) / `N_d1*N_d3` (d1) rows — i.e. whose optimal terminal a1prime is above the
bottom grid point.

**Bank exercised?** YES — this is the default terminal branch, hit by every basic and
lowmemory=1 withA1 leg (8 nosemiz subcodes plus the noz variants), with `N_d2=3 > 1`.

**Does the fix change live bank results? NO — analysed explicitly:** the two decodes coincide
whenever `maxindex <= N_d3` (nod1 files) / `maxindex <= N_d1*N_d3` (d1 files), because then
`rem(maxindex-1,N_dwrong)+1 = maxindex = rem(maxindex-1,N_dright)+1`-composed values and both
a1prime decodes give 1. The bank's terminal age is retirement (`agej>=Jr`), consumption
`c = pension + (1+r_a1)*a1 + a2 - a1prime - savings` is strictly decreasing in a1prime with
strictly increasing utility and no bequest motive, so the terminal argmax always has a1prime at
its bottom grid point → `maxindex` lies in the first `N_d3` (resp. `N_d1*N_d3`) rows → the
buggy and correct decodes produce IDENTICAL Policy at N_j (ties are impossible to matter: both
decodes read the same maxindex from the same ReturnMatrix). This is why the diary's
FromPolicy/DC/lowmemory checks are all zero today. The masking is knife-edge, though: any
terminal-age reason to hold a1 (warm glow, operative terminal constraint, a ReturnFn where
a1prime relaxes something) breaks it immediately, with out-of-range d3 as the failure mode.

**Fix** — minimal two-token change per block (`N_d` → `N_d3`, twice). Representative
(`nod1_raw` lines 52-55, lm0 block):

BEFORE
```matlab
        dindex=rem(maxindex-1,N_d)+1;
        Policy(1,:,:,N_j)=1; % is meaningless anyway
        Policy(2,:,:,N_j)=shiftdim(dindex,-1);
        Policy(3,:,:,N_j)=shiftdim(ceil(maxindex/N_d),-1);
```

AFTER
```matlab
        dindex=rem(maxindex-1,N_d3)+1;
        Policy(1,:,:,N_j)=1; % is meaningless anyway
        Policy(2,:,:,N_j)=shiftdim(dindex,-1);
        Policy(3,:,:,N_j)=shiftdim(ceil(maxindex/N_d3),-1);
```

Site list:
- `nod1_raw` 52-55 and 65-68: `N_d`→`N_d3` at lines 52/55 and 65/68 (the two blocks are
  textually identical up to `(:,:,N_j)` vs `(:,z_c,N_j)` indexing).
- `nod1_noz_raw` 45-48: same substitution at 45/48; **bundled cosmetic**: line 47's
  `shiftdim(ceil(dindex),-1)` — `ceil` of an integer index is a no-op, drop it to
  `shiftdim(dindex,-1)`.
- `raw` 53-57 and 67-71: `N_d`→`N_d1*N_d3` at lines 53/57 and 67/71 (dindex sub-decode lines
  54-56 are already correct given a correct dindex). Identical pair up to z indexing.
- `noz_raw` 46-50: `N_d`→`N_d1*N_d3` at lines 46/50.

---

## Item 3 — noz_e raws: stray extra `n_e` argument (runtime crash) in the terminal-V_Jplus1 lm1 e-loop

**Verified locations** (exactly as recorded):
- `Raw/ValueFnIter_FHorz_RiskyAsset_nod1_noz_e_raw.m:119`
- `Raw/ValueFnIter_FHorz_RiskyAsset_nod1_noa1_noz_e_raw.m:106`
- `Raw/ValueFnIter_FHorz_RiskyAsset_noa1_noz_e_raw.m:110`
- `Raw/ValueFnIter_FHorz_RiskyAsset_noz_e_raw.m:132`

Each passes 9 arguments (`..., n_e, special_n_e, ...`) to
`CreateReturnFnMatrix_Case2_Disc(ReturnFn, n_d, n_a, n_z, d_gridvals, a_gridvals, z_gridvals,
ReturnFnParamsVec)` — an 8-argument function (verified signature,
`ReturnFnMatrix/CreateReturnFnMatrix_Case2_Disc.m:1`).

**Severity**: CRASH — "Too many input arguments" the moment the branch runs.

**Trigger**: `vfoptions.V_Jplus1` + `vfoptions.lowmemory=1` in a noz+e riskyasset model (each
of the four shapes). **Bank**: never passes V_Jplus1 → unexercised.

**Fix**: drop the `n_e, ` (keep `special_n_e`), matching the main-loop lm1 call in the SAME
file (verified: `nod1_noz_e_raw:203`, `nod1_noa1_noz_e_raw:179`, `noa1_noz_e_raw:185`,
`noz_e_raw:222`). Representative (`nod1_noz_e_raw:119`):

BEFORE
```matlab
           ReturnMatrix_e=CreateReturnFnMatrix_Case2_Disc(ReturnFn, [n_d3,n_a1], [n_a1,n_a2], n_e, special_n_e, d3a1_gridvals, a1a2_gridvals, e_val, ReturnFnParamsVec);
```

AFTER
```matlab
           ReturnMatrix_e=CreateReturnFnMatrix_Case2_Disc(ReturnFn, [n_d3,n_a1], [n_a1,n_a2], special_n_e, d3a1_gridvals, a1a2_gridvals, e_val, ReturnFnParamsVec);
```

The other three sites are the same one-token deletion with their own gridvals arguments
(`n_d3, n_a` + `d3_gridvals, a_gridvals` in nod1_noa1; `n_d13, n_a` + `d13_gridvals,
a_gridvals` in noa1; `[n_d13,n_a1], [n_a1,n_a2]` + `d13a1_gridvals, a12_gridvals` in noz_e).

**Behaviour impact**: none on the bank (branch unreachable without V_Jplus1); converts a crash
into the intended lowmemory=1 terminal solve.

---

## Item 4 — noz raws: V_Jplus1 reshaped to [N_a2,1] instead of [N_a,1]

**Verified locations** (as recorded):
- `Raw/ValueFnIter_FHorz_RiskyAsset_nod1_noz_raw.m:53`:
  `EV=reshape(vfoptions.V_Jplus1,[N_a2,1]); % Using V_Jplus1`
- `Raw/ValueFnIter_FHorz_RiskyAsset_noz_raw.m:53`:
  `V_Jplus1=reshape(vfoptions.V_Jplus1,[N_a2,1]);    % First, switch V_Jplus1 into Kron form`

These are withA1 raws: V_Jplus1 has `N_a=N_a1*N_a2` elements, and the subsequent u-lottery
indexes it up to `N_a` (`aprimeIndex=repelem((1:1:N_a1)',N_d23,N_u)+N_a1*repmat(a2primeIndex-1,
N_a1,1)`, lines 59-62 in both files). All other 14 raws reshape with `N_a` (verified e.g.
`nod1_raw:75` `[N_a,N_z]`, `noz_e_raw:74` `[N_a,N_e]`, `e_raw:102` `[N_a,N_z,N_e]`).

**Severity**: CRASH — reshape element-count mismatch whenever `N_a1>1` (always, in these
withA1 shapes; if a user somehow passed `n_a1=1` the reshape succeeds and the indexing is then
coincidentally consistent).

**Trigger**: `vfoptions.V_Jplus1` in a noz withA1 riskyasset model. **Bank**: no V_Jplus1 →
unexercised.

**Fix** (one token per file, `N_a2`→`N_a`). Representative (`noz_raw:53`):

BEFORE
```matlab
    V_Jplus1=reshape(vfoptions.V_Jplus1,[N_a2,1]);    % First, switch V_Jplus1 into Kron form
```

AFTER
```matlab
    V_Jplus1=reshape(vfoptions.V_Jplus1,[N_a,1]);    % First, switch V_Jplus1 into Kron form
```

Same one-token change at `nod1_noz_raw:53` (variable is named `EV` there).

**Behaviour impact**: none on the bank.

---

## Item 5 — SemiExo noa1 raws: terminal V_Jplus1 lowmemory ladders shorter than the main-loop ladders

**Verified locations** (all in `ValueFnIter/FHorz/RiskyAsset/RiskyAssetSemiExo/`):

| File | Terminal V_Jplus1 ladder (lines) | Main-loop ladder (lines) | Missing terminal values |
|------|-----------------------------------|--------------------------|------------------------|
| `..._nod1_noa1_raw.m` | 0/1 (101, 137) | 0/1/2 (200, 236, 276) | 2 |
| `..._noa1_e_raw.m` | 0/1 (135, 176) | 0/1/2/3 (248, 289, 334, 382) | 2, 3 |
| `..._nod1_noa1_e_raw.m` | 0/1 (122, 158) | 0/1/2/3 (224, 260, 299, 341) | 2, 3 |
| `..._noa1_noz_e_raw.m` | 0/1 (113, 154) | 0/1/2 (226, 267, 312) | 2 |
| `..._nod1_noa1_noz_e_raw.m` | 0/1 (102, 138) | 0/1/2 (204, 240, 279) | 2 |

Exemplar that does it right: `..._noa1_raw.m`, whose terminal V_Jplus1 ladder covers 0/1/2
fully (lines 98, 139, 185).

**Severity**: SILENT-ZERO. With V_Jplus1 + a lowmemory value the main loop accepts but the
terminal ladder doesn't, the terminal age executes NO branch: `V(:,:,N_j)` stays 0 (and its
Policy stays 1s), then the main loop runs normally off that zero terminal — every age silently
wrong, no error.

**Trigger**: `vfoptions.V_Jplus1` + lowmemory 2 (or 3 for the z+e files) in the matching
SemiExo noa1 shape. **Bank**: no V_Jplus1 → unexercised (its semiz legs do use lowmemory up to
2/3, but only without V_Jplus1).

**Fix** — minimal: widen the final terminal-ladder guard from `==1` to `>=1`. This is valid
because each file's terminal lm1 block is already the most-looped variant for its shape
(verified: `nod1_noa1_raw`'s terminal lm1 at 137-172 loops over ALL bothz — the semantics of
the main loop's lm2; the e-files' terminal lm1 loops e inside d4 with the EV work per-d4 —
correct, if less memory-lean than a dedicated lm2/lm3 block, for any lowmemory>=1).
Representative (`nod1_noa1_raw:137`):

BEFORE
```matlab
    elseif vfoptions.lowmemory==1
```

AFTER
```matlab
    elseif vfoptions.lowmemory>=1 % terminal lm1 already loops over bothz, so it also serves lowmemory==2
```

Site list (one token each): `nod1_noa1_raw:137`, `noa1_e_raw:176`, `nod1_noa1_e_raw:158`,
`noa1_noz_e_raw:154`, `nod1_noa1_noz_e_raw:138`.

Alternative (bigger diff, memory parity): port each file's main-loop lm2/lm3 blocks into the
terminal ladder with `jj`→`N_j` and `V(:,:,jj+1)`→`V_Jplus1` (`EVpre` in the e-files), the way
`noa1_raw` does for lm2. Same numerical results; only worth it if terminal-age memory at
lowmemory 2/3 ever matters in practice.

Note: the EZ semiz raws mirrored this control flow verbatim under the list-don't-fix policy, so
whichever form is chosen here should be replicated there when the paused EZ-semiz batch resumes.

**Behaviour impact**: none on the bank.

---

## Item 6 — SemiExo nod1_raw: V_Jplus1 reshaped [N_a,N_z] instead of [N_a,N_bothz]

**Verified location**:
`ValueFnIter/FHorz/RiskyAsset/RiskyAssetSemiExo/ValueFnIter_FHorz_RiskyAssetSemiExo_nod1_raw.m:97`:
`V_Jplus1=reshape(vfoptions.V_Jplus1,[N_a,N_z]);` — a SemiExo file's V_Jplus1 has
`N_a*N_bothz` (=`N_a*N_semiz*N_z`) elements. Sibling `..._raw.m:98` uses `[N_a,N_bothz]`
correctly (as do `noa1_raw:83` and `nod1_noa1_raw:87`).

**Severity**: CRASH — reshape element mismatch whenever `N_semiz>1` (always, in a semiz model).
It fires at the top of the terminal-V_Jplus1 branch, before item 7's lines in this file.

**Trigger**: `vfoptions.V_Jplus1` in a nod1 withA1 riskyasset+semiz model. **Bank**: no
V_Jplus1 → unexercised.

**Fix** (one token):

BEFORE
```matlab
    V_Jplus1=reshape(vfoptions.V_Jplus1,[N_a,N_z]);    % First, switch V_Jplus1 into Kron form
```

AFTER
```matlab
    V_Jplus1=reshape(vfoptions.V_Jplus1,[N_a,N_bothz]);    % First, switch V_Jplus1 into Kron form
```

**Behaviour impact**: none on the bank. Must land together with item 7 (same branch).

---

## Item 7 — SemiExo nod1_raw + raw: unflattened EV u-lottery in the terminal-V_Jplus1 lm0/lm1 blocks

**Verified locations** (both files in `ValueFnIter/FHorz/RiskyAsset/RiskyAssetSemiExo/`):
- `..._nod1_raw.m`: lm0 block lines 132-142 (skipinterp 132, aprimeProbs 133-134, EV1/EV2
  137-138, reshapes 141-142); lm1 block lines 183-191 (same statements, fewer comments).
- `..._raw.m`: lm0 block lines 127-135; lm1 block lines 178-186. (Textually the same statements
  as nod1_raw's minus its extra comments.)

The bug: `EV(aprimeIndex+N_a*((1:1:N_bothz)-1))` with MATRIX `aprimeIndex` (`[N_d23*N_a1,N_u]`)
— implicit expansion clashes dim 2 (`N_u` vs `N_bothz`) → error whenever `N_u~=N_bothz` (and a
reshape-size error even in the freak `N_u==N_bothz` case). Additionally
`aprimeProbs=repmat(a2primeProbs,N_a1,1)` lacks the `N_bothz` expansion the flattened form
needs. The same files' MAIN LOOPS use the correct flattened form (verified: `nod1_raw` 324-327,
330-331 and 376-377, 381-382; `raw` 307-310, 312-313 and 354-357, 359-360).

**Severity**: CRASH. **Trigger**: `vfoptions.V_Jplus1` + lowmemory 0 or 1 in a withA1
riskyasset+semiz model (nod1 or d1 shape). (The lm2 blocks of both files use a per-z_c `EV_z`
and are fine.) **Bank**: no V_Jplus1 → unexercised.

**Fix**: copy the same file's main-loop flattened pattern verbatim. Representative
(`nod1_raw` lm0 terminal, lines 130-142):

BEFORE
```matlab
            % Seems like interpolation has trouble due to numerical precision rounding errors when the two points being interpolated are equal
            % So I will add a check for when this happens, and then overwrite those (by setting aprimeProbs to zero)
            skipinterp=logical(EV(aprimeIndex+N_a*((1:1:N_bothz)-1))==EV(aprimeplus1Index+N_a*((1:1:N_bothz)-1))); % Note, probably just do this off of a2prime values
            aprimeProbs=repmat(a2primeProbs,N_a1,1);  % [N_d*N_a1,N_u]
            aprimeProbs(skipinterp)=0;

            % Switch EV from being in terms of aprime to being in terms of d (in expectation because of the u shocks)
            EV1=EV(aprimeIndex+N_a*((1:1:N_bothz)-1)); % (d,a1prime,u,z), the lower aprime
            EV2=EV((aprimeplus1Index)+N_a*((1:1:N_bothz)-1)); % (d,a1prime,u,z), the upper aprime

            % Apply the aprimeProbs
            EV1=reshape(EV1,[N_d23*N_a1,N_u,N_bothz]).*aprimeProbs; % probability of lower grid point
            EV2=reshape(EV2,[N_d23*N_a1,N_u,N_bothz]).*(1-aprimeProbs); % probability of upper grid point
```

AFTER (the main-loop form, `nod1_raw` 324-335)
```matlab
            % Seems like interpolation has trouble due to numerical precision rounding errors when the two points being interpolated are equal
            % So I will add a check for when this happens, and then overwrite those (by setting aprimeProbs to zero)
            skipinterp=logical(EV(aprimeIndex(:)+N_a*((1:1:N_bothz)-1))==EV(aprimeplus1Index(:)+N_a*((1:1:N_bothz)-1))); % Note, probably just do this off of a2prime values
            aprimeProbs=repmat(a2primeProbs,N_a1,N_bothz);
            aprimeProbs(skipinterp)=0;
            aprimeProbs=reshape(aprimeProbs,[N_d23*N_a1,N_u,N_bothz]);

            % Switch EV from being in terms of aprime to being in terms of d (in expectation because of the u shocks)
            EV1=EV(aprimeIndex(:)+N_a*((1:1:N_bothz)-1)); % (d,u,z), the lower aprime
            EV2=EV(aprimeplus1Index(:)+N_a*((1:1:N_bothz)-1)); % (d,u,z), the upper aprime

            % Apply the aprimeProbs
            EV1=reshape(EV1,[N_d23*N_a1,N_u,N_bothz]).*aprimeProbs; % probability of lower grid point
            EV2=reshape(EV2,[N_d23*N_a1,N_u,N_bothz]).*(1-aprimeProbs); % probability of upper grid point
```

Site list (4 blocks, statements textually identical across sites): `nod1_raw` 132-142 (lm0) and
183-191 (lm1); `raw` 127-135 (lm0) and 178-186 (lm1). In each, three changes: `aprimeIndex`→
`aprimeIndex(:)` and `aprimeplus1Index`→`aprimeplus1Index(:)` in the three lookup lines,
`repmat(a2primeProbs,N_a1,1)`→`repmat(a2primeProbs,N_a1,N_bothz)`, and add the
`aprimeProbs=reshape(aprimeProbs,[N_d23*N_a1,N_u,N_bothz]);` line after the skipinterp zeroing.

**Behaviour impact**: none on the bank.

---

## Item 8 — SemiExo GI FromPolicy: missing isnan clear after the pi_e collapse (+ cosmetics)

**File**: `ValueFnFromPolicy/RiskyAsset/ValueFnFromPolicy_FHorz_RiskyAsset_SemiExo_GI.m`

**Verified location**: Step 3a, lines 223-227 (record said 224-227 — the sum is at line 225;
counting includes the surrounding assignment lines; effectively no drift):

**Severity**: SILENT-WRONG (latent). `0*(-Inf)=NaN` from a zero pi_e weight against an -Inf V'
poisons the whole e-sum; the downstream z/semiz-sum clears (lines 235, 245, 254) then convert
that NaN to 0 — turning states that should carry -Inf (or a finite value) into 0.

**Trigger**: riskyasset+semiz+e GI FromPolicy on a model with -Inf V' AND a zero entry in
pi_e_J. **Bank**: exercises the path (withA1 semiz e legs, GI FromPolicy check) but has neither
-Inf V nor zero pi_e entries → never fires.

**Fix** (pattern = the EZ SemiExo FromPolicy file, verified at its lines 267-269):

BEFORE (lines 223-227)
```matlab
            V_next=V(:,:,:,jj+1);
            V_next=sum(V_next .* shiftdim(vfoptions.pi_e_J(:,jj+1), -2), 3);
            V_next=reshape(V_next, [N_a, N_shocks]);
```

AFTER
```matlab
            V_next=V(:,:,:,jj+1);
            V_next=sum(V_next .* shiftdim(vfoptions.pi_e_J(:,jj+1), -2), 3);
            V_next(isnan(V_next))=0; % -Inf times zero e'-probability
            V_next=reshape(V_next, [N_a, N_shocks]);
```

**Bundled cosmetics** (verified): dead identical if/else at lines 118-122 (both branches are
`PolicyValuesPermute=permute(PolicyValues,[2,3,1,4]);` — collapse to the single statement);
assignment-only locals `a1_grid` (line 71), `l_a2` (73), `l_aprime` (74) — delete (each has
exactly 1 mention in the file). (`l_a` at line 54 is also assignment-only but was not on the
recorded list; delete or keep at the implementer's discretion.)

**Behaviour impact**: none on the bank.

---

## Item 9 — SemiExo GridInterpLayer noz family: lowmemory fall-through leaves V at zeros (+ lint)

**Directory**: `ValueFnIter/FHorz/RiskyAsset/RiskyAssetSemiExo/GridInterpLayer/`

**Verified ladder map** (no final `else` anywhere; the GI dispatcher
`ValueFnIter_FHorz_RiskyAssetSemiExo_GI.m` does no lowmemory range validation):

| File | Terminal no-VJp1 | Terminal VJp1 | Main loop | Fall-through values (with-z sibling accepts them) |
|------|------------------|---------------|-----------|--------------------------------------------------|
| `GI1_noz_raw.m` | 0/1 (68, 93) | 0/1 (165, 228) | 0/1 (351, 414) | 2 (sibling `GI1_raw`: 0/1/2) |
| `GI1_nod1_noz_raw.m` | 0/1 (63, 88) | 0/1 (154, 215) | 0/1 (328, 389) | 2 (sibling `GI1_nod1_raw`: 0/1/2) |
| `GI1_noz_e_raw.m` | 0/1/2 (71, 97, 126) | 0/>=1 (200, 263 — total) | 0/>=1 (380, 443 — total) | 3, terminal-no-VJp1 only (sibling `GI1_e_raw`: 0/1/2/3) |
| `GI1_nod1_noz_e_raw.m` | 0/1/2 (68, 94, 123) | 0/>=1 (199, 260 — total) | 0/>=1 (376, 437 — total) | 3, terminal-no-VJp1 only |

**Severity**: SILENT-ZERO. In the two noe files an out-of-range lowmemory executes NO branch in
any of the three ladders → the entire V returns as zeros, no error. In the two e files only the
terminal no-V_Jplus1 ladder falls through → terminal age stays zeros and every earlier age is
then silently wrong (their V_Jplus1/main ladders already use the `>=1` catch-all form).

**Trigger**: riskyasset+semiz noz GI solve with lowmemory=2 (noe shapes) or lowmemory=3 (e
shapes) — values a user carrying settings over from the with-z siblings would pass. **Bank**:
its GI noz legs pass only lowmemory 1 (noe) / 1-2 (e), all ladder-covered → unexercised.

**Fix** — minimal, using the family's own in-file good pattern (the `elseif
vfoptions.lowmemory>=1` catch-alls already present in the two e files): widen the final guard
of every fall-through ladder. Representative (`GI1_noz_raw:414`, main loop):

BEFORE
```matlab
    elseif vfoptions.lowmemory==1
```

AFTER
```matlab
    elseif vfoptions.lowmemory>=1
```

Site list (one token each, 8 sites): `GI1_noz_raw` 93, 228, 414; `GI1_nod1_noz_raw` 88, 215,
389 (all `==1`→`>=1`); `GI1_noz_e_raw` 126 and `GI1_nod1_noz_e_raw` 123 (`==2`→`>=2`).

**Bundled lint** (verified, and refined vs the record): `a2ind` is assignment-only in ALL 8 GI1
files (e.g. `GI1_noz_raw:56`, `GI1_nod1_noz_raw:51`) — delete in the 4 noz files this batch
touches; the `%#ok<NASGU>` pragmas on `special_n_d4` are STALE (the variable IS used, 13-25
mentions per file) — drop the pragma in the touched files. Two recorded claims did not
reproduce: `zind` IS used (3-5 mentions, e.g. `GI1_nod1_noz_raw` 191/365), and the
`n_d13`/`d13_grid` NASGU pragmas in `GI1_noz_raw` are accurate (those really are unused), so
leave both alone. Leave the 4 with-z GI1 files' identical lint untouched (not in this batch).

**Adjacent observation (NOT part of the recorded 10, listed for completeness)**: the base
(non-GI) SemiExo noz raws (`..._RiskyAssetSemiExo_noz_raw.m`, `..._nod1_noz_raw.m`) have the
same ladder shape (all three ladders end at `elseif vfoptions.lowmemory==1`, no else) while
their with-z siblings accept lowmemory=2 — the same silent-zero fall-through class. If desired,
apply the same `==1`→`>=1` widening there in the same pass; strictly optional and out of scope
for this batch.

**Behaviour impact**: none on the bank.

---

## Item 10 — SemiExo FromPolicy (non-GI): missing isnan clear after the pi_e sum (+ cosmetics)

**File**: `ValueFnFromPolicy/RiskyAsset/ValueFnFromPolicy_FHorz_RiskyAsset_SemiExo.m`

**Verified location** (as recorded): line 215,
`V_next=sum(V_next .* shiftdim(vfoptions.pi_e_J(:,jj+1), -2), 3);` with no clear before the
reshape at 216. Every other probability-weighted sum in the file has its clear (z' at 225,
semiz' at 235/244, per-u in Step 4). Same mechanism as item 8: the NaN survives to the
z'/semiz' sums whose clears then zero entries that should be -Inf/finite.

**Severity**: SILENT-WRONG (latent). **Trigger**: riskyasset+semiz+e FromPolicy with -Inf V'
and a zero pi_e entry. **Bank**: path exercised in every semiz e leg's FromPolicy check;
trigger absent → never fires.

**Fix** (identical insertion to item 8; pattern = EZ SemiExo FromPolicy lines 267-269):

BEFORE (lines 213-217)
```matlab
            V_next=V(:,:,:,jj+1);
            V_next=sum(V_next .* shiftdim(vfoptions.pi_e_J(:,jj+1), -2), 3);
            V_next=reshape(V_next, [N_a, N_shocks]);
```

AFTER
```matlab
            V_next=V(:,:,:,jj+1);
            V_next=sum(V_next .* shiftdim(vfoptions.pi_e_J(:,jj+1), -2), 3);
            V_next(isnan(V_next))=0; % -Inf times zero e'-probability
            V_next=reshape(V_next, [N_a, N_shocks]);
```

**Bundled cosmetics** (verified): dead identical if/else at lines 119-123 (both branches
`permute(PolicyValues,[2,3,1,4])`); assignment-only locals `l_a` (line 58), `N_a2` (73),
`a1_grid` (74), `l_a2` (76), `l_aprime` (77) — each has exactly 1 mention; delete.

**Behaviour impact**: none on the bank.

---

## Closing note — test-bank work (NOT part of this fix batch)

Two test-side suggestions already on record accompany this batch; they are vfitoolkitTests
work, to be done separately from (ideally after) the toolkit fixes:

1. **A V_Jplus1 test leg.** Items 3, 4, 5, 6, and 7 all live in terminal-V_Jplus1 branches the
   bank never enters, and item 2's masking means the default terminal branch is only tested at
   a corner. A leg that solves the last J2 ages with `vfoptions.V_Jplus1` taken from a full-run
   V (the jstar convention already used by the CoreFHorz V_Jplus1 subcodes) and checks equality
   with the full solve, swept over each shape's full lowmemory ladder, would have caught all
   five items — and is the only way to GPU-validate their fixes. Follow the existing
   CoreFHorzTests V_Jplus1-block pattern (32 subcodes, tiers, age-dependent-shock variant).

2. **Strengthen the d2recon cross-test calibration whose dsemiz is constant.** The GPU run of
   the EZ bank reported that the d2recon regression guard's own diagnostic shows dsemiz is
   CONSTANT under the calibration (all cases), making the guard vacuous; the same applies to
   the vNM d2recon calibration (`CoreFHorzRiskyAssetTests_subcodes/CrossTests/
   CoreFHorzRiskyAsset_CrossTests_d2recon_nod1_z_semiz_withA1.m`). Strengthening the
   semiz/portfolio coupling so the optimal d4 actually varies across states would make the
   reconstruction check bite.

Additionally (new, from item 2's analysis): the bank's terminal-age corner solution is what
masks item 2 — a variant whose terminal age has a reason to choose a1prime above the bottom
grid point (e.g. terminal-age income large enough relative to pension, or a warm-glow-style
term) would make the terminal Policy decode observable in the FromPolicy zero-checks. Worth
folding into whichever of the two suggestions above is done first.
