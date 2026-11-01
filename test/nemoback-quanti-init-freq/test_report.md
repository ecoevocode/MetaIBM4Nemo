# quanti_init_freq (NemoBackend, Python redraw of founder alleles) -- test report

Generated 2026-10-31 01:22 (re-run after continuous/pleiotropy ValueError guard; prior runs 2026-10-30 20:07, 22:53) by `test_quanti_init_freq.py` + `test_ecology.py` (this folder). Scope: the uncommitted change in `metaibm/indfactory.py` on top of `05f050b` -- `indfactory(..., metaibm_nemo_params={'quanti_init_freq': ...})`, `_use_quanti_init_freq()`, `_set_quanti_init_alleles()` (np.random draw `U >= p` -> allele 1, `TTrait.set_allele_bit`, `set_value()`), the `_quanti_seq_pos` cache and the one-line hook in `make_founder()`. Test-only: nothing outside this folder was edited, nothing committed, metaibm/ imported read-only (`sys.dont_write_bytecode`, `-B`; the existing `metaibm/__pycache__` dates from Oct 9-29 and was not touched).

## Environment

- HEAD `05f050b` (`FIX a bug: in habitat.py: hab_initialize(reproduce_mode='mix') raised UnboundLocalError ...`), branch `features/developing`; `metaibm/indfactory.py` at HEAD is the same as at `7fa64f6` (the first run), so the HEAD == tree regressions compare against the same baseline
- conda env `metaibm-nemo-test3` (pynemo with `TTrait.set_allele_bit`; `metaibm-nemo-test2` lacks it), python `3.11.18`, numpy `2.4.6`, `macOS-15.7.3-arm64`, Apple M2, 8 CPUs
- every check group in its own child process (Nemo's RAND, `Individual::currentID` and prototypes are process-global); Nemo output in a temp dir, removed afterwards

| file | md5 |
|---|---|
| `metaibm/indfactory.py` (working tree, under test) | `039174c4c2ad2baa435e9a0f0d1d9d87` |
| `test/nemoback-quanti-init-freq/test_quanti_init_freq.py` | `4eae956329eee5f1f5fcc43cc8372bd9` |
| `test/nemoback-quanti-init-freq/test_ecology.py` | `0edf45d9bf7944037ed4d47eda7e2a15` |

Run:

```
cd test/nemoback-quanti-init-freq
PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test3/bin/python -B test_quanti_init_freq.py   # ~25 s, exit 1 on any FAIL
PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test3/bin/python -B test_ecology.py           # ~12 s
```

(one group: append its name, e.g. `test_quanti_init_freq.py stats rng`).

## Contract under test

- The user passes `quanti_init_freq` only as `0.8` | `[0.8]` | `[[0.8]]` | `[p_1, ..., p_L]` (L = all quanti loci, trait by trait). Malformed shapes (two rows, wrong length) are out of scope and not tested.
- `p` = frequency of allele 0 (`+a` in `diallelic` / `diallelic_HC`, the ancestral allele 0 in `diallelic_asymmetrical`), as in Nemo's `LCE_QuantiInit` (`U < p` -> allele 0).
- Continuous allele models and active pleiotropy (`full`/`f`/`var`/`variable`/`v`) are explicitly rejected with `ValueError` in `_use_quanti_init_freq()`. Combinations left to the user (D3 `quanti_init_model != 0`, D4 asymmetrical + scalar `quanti_allele_value`) are recorded as NOTE, not judged.

## Results: 116 PASS, 0 FAIL, 7 NOTE

| script / group | covers | PASS | FAIL | NOTE |
|---|---|---|---|---|
| `test_quanti_init_freq.py` | | **103** | **0** | **6** |
| - forms | the four forms, 3 di-allelic models x L = 5, 3+2, 20+20, 100+100, 4 seeds: identical founders, and identical to the reference `np.random.random((2, L)) >= p` (strand-major); 200 consecutive founders with distinct per-locus p | 14 | 0 | 0 |
| - mapping | p in {0, 1}: column k -> trait-by-trait locus k -> sequence position; 5 loci layouts (1 to 3 traits), 4 patterns, scalar and locus-specific allele values, both genders: bits, allele values, phenotype exact; `_quanti_seq_pos == range(L)` | 16 | 0 | 0 |
| - stats | 6000 founders, p = 0.05..0.95: allele-0 frequency, per strand, Hardy-Weinberg, independence of strands and loci, genotypic mean and variance (analytic), Ve residual N(0, Ve) | 11 | 0 | 0 |
| - rng | exactly 2L np uniforms per founder, none without the parameter (4 ways to leave it out), Python `random` untouched, reproducible in-process and across processes, alleles independent of indfactory seed / gender / home_habitat | 11 | 0 | 0 |
| - activation | `diallelic`, `diallelic_HC`, `diallelic_asymmetrical` (and `quanti_pleiotropy='no'/'n'`) redraw; the 3 continuous models and pleiotropy (`full`/`f`/`var`/`variable`/`v`) raise `ValueError` without touching `np.random` | 8 | 0 | 0 |
| - lifecycle | `make_clone` / `make_copy` keep the redrawn alleles; offspring Mendelian (random-p parents and homozygous parents), offspring phenotype recomputed, offspring never redrawn | 11 | 0 | 0 |
| - genetics | dominance (k 0.3, 1.5) and dense epistasis: the phenotype `set_value()` leaves equals the one Nemo computes at birth for the same alleles (zero-mutation clone), with Ve too | 4 | 0 | 0 |
| - recycle | founders built on recycled Nemo objects are fully redrawn, phenotype fresh, counters | 3 | 0 | 0 |
| - multi_species | 4 species (3 redrawn, 1 continuous) interleaved: each consumes only its slice of one np stream; per-species caches | 2 | 0 | 0 |
| - set_allele_bit (absorbed) | the `TTrait.set_allele_bit` wrapper on a running Nemo sim: symmetric, asymmetrical, continuous no-op, non-quanti no-op (from `Claude outputs/test_set_allele_bit.py`, one sim per process) | 4 | 0 | 0 |
| - nemo_lce | Nemo's own LCE `quanti_init` (patch of 2 = Bernoulli branch, of 8 = exact-count branch) vs MetaIBM founders, 5 configs incl. locus-specific allele values: bits and phenotypes identical | 10 | 0 | 1 |
| - regression | no `quanti_init_freq`: working tree == HEAD bit for bit (founders, offspring, clones, copies, IDs, np and random state), 3 configs incl. example-like continuous | 9 | 0 | 0 |
| - user_notes / init_value / perf | D3, D4, D7, `quanti_init_value` on a di-allelic species, cost | 0 | 0 | 5 |
| `test_ecology.py` | | **13** | **0** | **1** |
| - init_stream | `meta_initialize` -> `patch_initialize` -> `hab_initialize` -> `make_founder`: 400 founders equal the reference np stream walked in patch / habitat / row / col order; invariants; both genders; control species untouched; minted == recycled | 5 | 0 | 0 |
| - cycle_asexual / sexual / mix | example2's life cycle on a mainland + 4 islands, 30 steps: every step every individual phenotype == genotype, the p = 0 species stays all allele 1, species <-> factory; minted == recycled at the end | 6 | 0 | 0 |
| - repro | sexual life cycle, same seeds, two processes: identical landscape md5 | 1 | 0 | 0 |
| - regression_cycle | no species with the parameter: 20-step sexual life cycle identical on HEAD and on the working tree | 1 | 0 | 0 |
| - mix_init_note | `meta_initialize(reproduce_mode='mix')` (outside this change; fixed in `05f050b`, now accepted) | 0 | 0 | 1 |

## Findings

No defect found in the change under the stated contract. The redraw matches Nemo's `LCE_QuantiInit` allele for allele where both are deterministic (p in {0, 1}), follows `U < p -> allele 0` exactly, and leaves every run without the parameter bit-identical to HEAD.

NOTEs (behaviour recorded, decisions are the user's):

1. **Continuous models and pleiotropy now raise `ValueError`**: `_use_quanti_init_freq()` rejects `continuous`, `continuous_in_place`, `continuous_HC` and `quanti_pleiotropy` in `('full', 'f', 'var', 'variable', 'v')` (`np.random` untouched; verified as PASS in `activation`).
2. **D3** `quanti_init_model` 0/1/2/3 + `quanti_init_freq`: founders bit-identical for one np seed, the redraw overwrites every allele.
3. **D4** `diallelic_asymmetrical` with scalar `quanti_allele_value 0.1`: p = 0 gives phenotype (-0.6, -0.4), i.e. symmetric +-a; with the matrix `{{0.1}}` it gives (0.6, 0.4).
4. **D7** full pleiotropy (2 traits x 3 loci): all input forms (`0.8`, `[0.8]`, `[[0.8]]`, `[0.8]*3`, `[0.8]*6`) now uniformly raise `ValueError` in `_use_quanti_init_freq()`.
5. `quanti_init_value` set on a di-allelic species: Nemo prints at boot `***ERROR*** "quanti_allele_model diallelic" needs to be initialized with "quanti_init_freq" provided by LCE "quanti_init" instead of "quanti_init_value".` (non-fatal); the founder is still redrawn.
6. Nemo LCE vs MetaIBM at p = 0.3, N = 2000: Nemo allocates exact counts (0.3 at every locus), MetaIBM draws Bernoulli (0.295-0.307) -- the documented difference (D7, expected frequency vs exact allocation).
7. Cost: `make_founder` 1.7 us without the parameter, ~27-29 us with it at L = 40, ~129 us at L = 200 (2L `set_allele_bit` SWIG calls + `set_value`); continuous example-like founder 2.2 us.
8. Outside this change: `meta_initialize(reproduce_mode='mix')` raised `UnboundLocalError: cannot access local variable 'gender'` in the first run (`habitat.hab_initialize` set gender only for `'sexual'` / `'asexual'`). Fixed in `05f050b` (`if reproduce_mode in ('sexual', 'mix')`): this run records `reproduce_mode="mix" accepted (400 founders)`. `cycle_mix` still initializes with `'sexual'`, so its results are unchanged.

Observations:

- Nemo's offspring layout in this pynemo build: strand 0 carries the father's gamete, strand 1 the mother's (homozygous mother all allele 0 x father all allele 1 -> every offspring strand 0 all 1, strand 1 all 0; random-p parents 400/400 consistent with "strand 0 from father").
- Test pitfall met on the way (test code, fixed): `patch.add_habitat` draws the microsite environments from `np.random`, so `np.random.seed` must come before building the landscape, otherwise two runs diverge at the first `meta_dead_selection`.
- `Claude outputs/test_set_allele_bit.py` (pytest, pynemo high-level API) is absorbed as the four `sab_*` groups; the original file has since been deleted by the user.

## All checks

Copied verbatim from the two runs above.

### test_quanti_init_freq.py

| status | id | detail |
|---|---|---|
| PASS | `forms/diallelic/loci(5,)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/diallelic/loci(3, 2)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/diallelic/loci(20, 20)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/diallelic/loci(100, 100)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/diallelic_HC/loci(5,)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/diallelic_HC/loci(3, 2)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/diallelic_HC/loci(20, 20)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/diallelic_HC/loci(100, 100)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/diallelic_asymmetrical/loci(5,)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/diallelic_asymmetrical/loci(3, 2)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/diallelic_asymmetrical/loci(20, 20)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/diallelic_asymmetrical/loci(100, 100)` | four forms x 4 seeds x 3 founders == reference (U >= p, strand-major) |
| PASS | `forms/per-locus/diallelic` | 200 consecutive founders, p=[0.05, 0.3, 0.5, 0.7, 0.95]: bit-identical to the reference stream |
| PASS | `forms/per-locus/diallelic_asymmetrical` | 200 consecutive founders, p=[0.05, 0.3, 0.5, 0.7, 0.95]: bit-identical to the reference stream |
| PASS | `mapping/diallelic/loci(4,)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic/loci(3, 2)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic/loci(2, 3)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic/loci(2, 2, 2)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic/loci(1, 4)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic_HC/loci(4,)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic_HC/loci(3, 2)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic_HC/loci(2, 3)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic_HC/loci(2, 2, 2)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic_HC/loci(1, 4)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic_asymmetrical/loci(4,)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic_asymmetrical/loci(3, 2)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic_asymmetrical/loci(2, 3)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic_asymmetrical/loci(2, 2, 2)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/diallelic_asymmetrical/loci(1, 4)` | 4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact |
| PASS | `mapping/seq_pos-identity` | _quanti_seq_pos=[0, 1, 2, 3, 4] |
| PASS | `stats/diallelic/allele0-freq` | p=[0.05, 0.3, 0.5, 0.7, 0.95] observed=[0.046, 0.3025, 0.5056, 0.7043, 0.9467] max\|z\|=2.01 |
| PASS | `stats/diallelic/per-strand` | strand0=[0.0483, 0.3038, 0.5035, 0.7088, 0.9433] strand1=[0.0437, 0.3012, 0.5077, 0.6998, 0.95] max\|z\|=2.37 |
| PASS | `stats/diallelic/hardy-weinberg` | chi2 (df 2, known p) per locus=[5.47, 0.36, 2.85, 1.48, 2.81] (limit 18.4 ~ P 1e-4) |
| PASS | `stats/diallelic/independence` | max \|r\| between allele columns (strands x loci)=0.0322 (limit 0.0581) |
| PASS | `stats/diallelic/genotypic-mean-var` | mean=[-0.2584, 0.2604] (expect [-0.26, 0.26]) var/expected=[0.989, 0.996] |
| PASS | `stats/diallelic_asymmetrical/allele0-freq` | p=[0.05, 0.3, 0.5, 0.7, 0.95] observed=[0.046, 0.3025, 0.5056, 0.7043, 0.9467] max\|z\|=2.01 |
| PASS | `stats/diallelic_asymmetrical/per-strand` | strand0=[0.0483, 0.3038, 0.5035, 0.7088, 0.9433] strand1=[0.0437, 0.3012, 0.5077, 0.6998, 0.95] max\|z\|=2.37 |
| PASS | `stats/diallelic_asymmetrical/hardy-weinberg` | chi2 (df 2, known p) per locus=[5.47, 0.36, 2.85, 1.48, 2.81] (limit 18.4 ~ P 1e-4) |
| PASS | `stats/diallelic_asymmetrical/independence` | max \|r\| between allele columns (strands x loci)=0.0322 (limit 0.0581) |
| PASS | `stats/diallelic_asymmetrical/genotypic-mean-var` | mean=[0.4292, 0.0698] (expect [0.43, 0.07]) var/expected=[0.989, 0.996] |
| PASS | `stats/Ve-residual` | phenotype - genotype: mean=[0.0004, 0.0005] sd=[0.1012, 0.0978] (expect 0, 0.1); corr(E, G)=-0.0263 |
| PASS | `rng/np-consumption/L5` | 7 founders consume exactly 7 x 2L = 70 uniforms |
| PASS | `rng/python-random-untouched/L5` | random.getstate() unchanged by make_founder |
| PASS | `rng/np-consumption/L40` | 7 founders consume exactly 7 x 2L = 560 uniforms |
| PASS | `rng/python-random-untouched/L40` | random.getstate() unchanged by make_founder |
| PASS | `rng/no-freq/no-argument` | np.random untouched, _quanti_seq_pos stays None, _use_quanti_init_freq() False |
| PASS | `rng/no-freq/None` | np.random untouched, _quanti_seq_pos stays None, _use_quanti_init_freq() False |
| PASS | `rng/no-freq/empty-dict` | np.random untouched, _quanti_seq_pos stays None, _use_quanti_init_freq() False |
| PASS | `rng/no-freq/freq-None` | np.random untouched, _quanti_seq_pos stays None, _use_quanti_init_freq() False |
| PASS | `rng/in-process-repro` | two factories, same seeds: bits+phenotypes md5 f3afd17dfa81 |
| PASS | `rng/independent-of-nemo-seed-gender-home` | indfactory seed 1/999/123456 x gender x home_habitat 0/7/65535: identical alleles for one np seed |
| PASS | `activation/diallelic` | redrawn: p=0 -> every allele is allele 1 |
| PASS | `activation/diallelic_HC` | redrawn: p=0 -> every allele is allele 1 |
| PASS | `activation/diallelic_asymmetrical` | redrawn: p=0 -> every allele is allele 1 |
| PASS | `activation/continuous` | raises ValueError (species 'sp4': quanti_init_freq requires a di-allelic model without pleiotropy), np.random untouched |
| PASS | `activation/continuous_in_place` | raises ValueError (species 'sp5': quanti_init_freq requires a di-allelic model without pleiotropy), np.random untouched |
| PASS | `activation/continuous_HC` | raises ValueError (species 'sp6': quanti_init_freq requires a di-allelic model without pleiotropy), np.random untouched |
| PASS | `activation/pleiotropy-rejected` | full/f/var/variable/v all raise ValueError, np.random untouched |
| PASS | `activation/pleiotropy-no-accepted` | quanti_pleiotropy='no'/'n' accepted and redrawn |
| PASS | `lifecycle/diallelic/clone-copy` | make_clone and make_copy keep all 2L alleles (mutation rate 0) |
| PASS | `lifecycle/diallelic/mendelian` | 400 offspring of random-p parents, consistent with: {'s0 from mother': 34, 's0 from father': 400, 'no consistent origin': 0} |
| PASS | `lifecycle/diallelic/mendelian-homozygous-parents` | mother all allele 0 x father all allele 1 -> 200/200 offspring heterozygous at every locus; strand 0 patterns [(1, 1, 1, 1, 1)], strand 1 patterns [(0, 0, 0, 0, 0)] |
| PASS | `lifecycle/diallelic/offspring-phenotype` | offspring phenotype == sum of its allele values (Ve 0) |
| PASS | `lifecycle/diallelic/offspring-not-redrawn` | make_offspring/make_clone draw no np.random (only founders are redrawn) |
| PASS | `lifecycle/diallelic_asymmetrical/clone-copy` | make_clone and make_copy keep all 2L alleles (mutation rate 0) |
| PASS | `lifecycle/diallelic_asymmetrical/mendelian` | 400 offspring of random-p parents, consistent with: {'s0 from mother': 29, 's0 from father': 400, 'no consistent origin': 0} |
| PASS | `lifecycle/diallelic_asymmetrical/mendelian-homozygous-parents` | mother all allele 0 x father all allele 1 -> 200/200 offspring heterozygous at every locus; strand 0 patterns [(1, 1, 1, 1, 1)], strand 1 patterns [(0, 0, 0, 0, 0)] |
| PASS | `lifecycle/diallelic_asymmetrical/offspring-phenotype` | offspring phenotype == sum of its allele values (Ve 0) |
| PASS | `lifecycle/diallelic_asymmetrical/offspring-not-redrawn` | make_offspring/make_clone draw no np.random (only founders are redrawn) |
| PASS | `lifecycle/phenotype-views` | phenotype_set == get_indi_phenotype_ls == (-0.6, -0.4) for p=0 |
| PASS | `genetics/dominance k=0.3, diallelic (3, 2)` | 60 founders: \|founder - Nemo clone\| max 0.00e+00; 60/60 differ from the additive sum (the map is active) |
| PASS | `genetics/dominance k=1.5, asymmetrical (2, 2)` | 60 founders: \|founder - Nemo clone\| max 0.00e+00; 53/60 differ from the additive sum (the map is active) |
| PASS | `genetics/epistasis dense eps=0.5, diallelic (5,)` | 60 founders: \|founder - Nemo clone\| max 0.00e+00; 39/60 differ from the additive sum (the map is active) |
| PASS | `genetics/dominance + Ve 0.01 (clone redraws Ve: compare genotype part only)` | founder - clone = difference of two independent Ve draws: mean 0.0027 sd 0.1343 (expect 0, 0.1414) |
| PASS | `recycle/redraw-on-reused-objects` | 300 founders after 300 recycles == reference stream (no stale alleles) |
| PASS | `recycle/counters` | minted=600 recycled=300 |
| PASS | `recycle/phenotype-fresh` | reused objects: phenotype recomputed (-0.6, -0.4) |
| PASS | `multi_species/interleaved-streams` | 4 species (3 redrawn, 1 continuous) interleaved 100 founders: each matches its slice of one np stream |
| PASS | `multi_species/caches` | per-species _quanti_seq_pos lengths 5/4/40, continuous stays None |
| NOTE | `user/D3-init-model` | quanti_init_model 0/1/2/3 + quanti_init_freq: founders bit-identical for one np seed -> the redraw overwrites every allele |
| NOTE | `user/D4-asym-scalar` | p=0 (all allele 1): matrix {{0.1}} -> [0.6, 0.4], scalar 0.1 -> [-0.6, -0.4] (scalar is read as symmetric +-a) |
| NOTE | `user/D7-full-pleiotropy` | 2 traits x 3 loci: 0.8: ValueError; [0.8]: ValueError; [[0.8]]: ValueError; [0.8]*num_locus(3): ValueError; [0.8]*seq_len(6): ValueError |
| NOTE | `perf/make_founder` | diallelic L=40, no freq: 1.7 us/founder; diallelic L=40, freq [p]*40: 28.7 us/founder; diallelic L=40, freq 0.5: 27.3 us/founder; diallelic L=200, freq 0.5: 128.9 us/founder; continuous L=40 (example-like): 2.2 us/founder |
| PASS | `set_allele_bit/symmetric-diallelic` | strand 0 -> allele 1: bits read back, +0.5 -> -0.5 on strand 0 only, phenotype 10 -> 0 (expect 10 -> 0) |
| PASS | `set_allele_bit/asymmetrical-diallelic` | allele 0 = 0, allele 1 = a: strand 0 -> allele 1 gives phenotype 5 (expect 5) |
| PASS | `set_allele_bit/continuous-noop` | continuous model: TTQuanti_continuous::set_allele_bit is empty -> values and phenotype unchanged |
| PASS | `set_allele_bit/non-quanti-noop` | ntrl trait (slot 0): every allele written with its opposite -> no crash, nothing changes; quant phenotype 10 (expect 10) |
| PASS | `nemo_lce/diallelic/loci(3, 2)/[1, 0, 0, 1, 0]/cap2` | Nemo LCE (2 inds) == MetaIBM founders: bits [[0, 1, 1, 0, 1], [0, 1, 1, 0, 1]], phenotype [-0.2, 0.0] |
| PASS | `nemo_lce/diallelic/loci(3, 2)/[1, 0, 0, 1, 0]/cap8` | Nemo LCE (8 inds) == MetaIBM founders: bits [[0, 1, 1, 0, 1], [0, 1, 1, 0, 1]], phenotype [-0.2, 0.0] |
| PASS | `nemo_lce/diallelic/loci(2, 3)/[0, 1, 1, 0, 1]/cap2` | Nemo LCE (2 inds) == MetaIBM founders: bits [[1, 0, 0, 1, 0], [1, 0, 0, 1, 0]], phenotype [0.2, 0.8] |
| PASS | `nemo_lce/diallelic/loci(2, 3)/[0, 1, 1, 0, 1]/cap8` | Nemo LCE (8 inds) == MetaIBM founders: bits [[1, 0, 0, 1, 0], [1, 0, 0, 1, 0]], phenotype [0.2, 0.8] |
| PASS | `nemo_lce/diallelic_asymmetrical/loci(3, 2)/[1, 0, 0, 1, 0]/cap2` | Nemo LCE (2 inds) == MetaIBM founders: bits [[0, 1, 1, 0, 1], [0, 1, 1, 0, 1]], phenotype [0.4, 0.2] |
| PASS | `nemo_lce/diallelic_asymmetrical/loci(3, 2)/[1, 0, 0, 1, 0]/cap8` | Nemo LCE (8 inds) == MetaIBM founders: bits [[0, 1, 1, 0, 1], [0, 1, 1, 0, 1]], phenotype [0.4, 0.2] |
| PASS | `nemo_lce/diallelic_asymmetrical/loci(4,)/[0, 0, 1, 1]/cap2` | Nemo LCE (2 inds) == MetaIBM founders: bits [[1, 1, 0, 0], [1, 1, 0, 0]], phenotype [0.3] |
| PASS | `nemo_lce/diallelic_asymmetrical/loci(4,)/[0, 0, 1, 1]/cap8` | Nemo LCE (8 inds) == MetaIBM founders: bits [[1, 1, 0, 0], [1, 1, 0, 0]], phenotype [0.3] |
| PASS | `nemo_lce/diallelic_HC/loci(5,)/[0, 0, 0, 0, 0]/cap2` | Nemo LCE (2 inds) == MetaIBM founders: bits [[1, 1, 1, 1, 1], [1, 1, 1, 1, 1]], phenotype [-1.0] |
| PASS | `nemo_lce/diallelic_HC/loci(5,)/[0, 0, 0, 0, 0]/cap8` | Nemo LCE (8 inds) == MetaIBM founders: bits [[1, 1, 1, 1, 1], [1, 1, 1, 1, 1]], phenotype [-1.0] |
| NOTE | `nemo_lce/p0.3-exact-vs-bernoulli` | Nemo LCE (N=2000, exact counts) allele-0 freq=[0.3, 0.3, 0.3, 0.3, 0.3]; MetaIBM (N=2000, Bernoulli) =[0.3048, 0.2948, 0.307, 0.3033, 0.3048] |
| PASS | `regression/HEAD==tree/asym-default` | HEAD d8ba51254f25 \| tree d8ba51254f25 |
| PASS | `regression/HEAD==tree/asym-default/np` | HEAD e13d2ece11ee \| tree e13d2ece11ee |
| PASS | `regression/HEAD==tree/asym-default/random` | HEAD 44dca204ba92 \| tree 44dca204ba92 |
| PASS | `regression/HEAD==tree/continuous-example` | HEAD 015666ef15a9 \| tree 015666ef15a9 |
| PASS | `regression/HEAD==tree/continuous-example/np` | HEAD e13d2ece11ee \| tree e13d2ece11ee |
| PASS | `regression/HEAD==tree/continuous-example/random` | HEAD 44dca204ba92 \| tree 44dca204ba92 |
| PASS | `regression/HEAD==tree/diallelic-default` | HEAD 1223031e392d \| tree 1223031e392d |
| PASS | `regression/HEAD==tree/diallelic-default/np` | HEAD e13d2ece11ee \| tree e13d2ece11ee |
| PASS | `regression/HEAD==tree/diallelic-default/random` | HEAD 44dca204ba92 \| tree 44dca204ba92 |
| NOTE | `user/init-value-diallelic` | founder still redrawn (all allele 1); Nemo boot messages: see parent detail \| Nemo printed: ['***ERROR*** "quanti_allele_model diallelic" needs to be initialized with "quanti_init_freq" provided by LCE "quanti_init" instead of "quanti_init_value".'] |
| PASS | `rng/cross-process-repro` | two child processes: md5 f3afd17dfa81 / f3afd17dfa81 |

PASS 103  FAIL 0  NOTE 6  (24s)  indfactory.py md5 039174c4c2ad2baa435e9a0f0d1d9d87

### test_ecology.py

| status | id | detail |
|---|---|---|
| PASS | `ecology/init/stream` | 400 founders via meta_initialize: species by habitat env, sp1-sp3 bits == reference np stream in patch/habitat/row/col order |
| PASS | `ecology/init/invariants` | 400 founders: phenotype == genotype, sp3 all allele 1; allele-0 freq sp1=[0.18, 0.345, 0.61, 0.825, 0.475] sp2=[0.535, 0.515, 0.45, 0.465] |
| PASS | `ecology/init/genders` | sexual init: both genders redrawn alike {'sp1': ['female', 'male'], 'sp2': ['female', 'male'], 'sp3': ['female', 'male'], 'sp4': ['female', 'male']} |
| PASS | `ecology/init/sp4-untouched` | control species never enters the redraw |
| PASS | `ecology/init/recycled` | after meta_clear_up: minted == recycled per species {'sp1': (100, 100), 'sp2': (100, 100), 'sp3': (100, 100), 'sp4': (100, 100)} |
| PASS | `ecology/asexual/invariants` | 30 steps, mainland + islands, every step: phenotype == genotype, sp3 all allele 1, species <-> factory; (t, alive, allele-0 freq) [(0, 400, {'sp1': [0.15, 0.44, 0.59, 0.77, 0.55], 'sp2': [0.48, 0.535, 0.48, 0.46]}), (10, 500, {'sp1': [0.121, 0.338, 0.596, 0.833, 0.517], 'sp2': [0.37, 0.602, 0.531, 0.492]}), (20, 790, {'sp1': [0.245, 0.383, 0.461, 0.721, 0.518], 'sp2': [0.415, 0.635, 0.442, 0.514]}), (30, 790, {'sp1': [0.208, 0.326, 0.547, 0.758, 0.5], 'sp2': [0.415, 0.657, 0.458, 0.54]})] |
| PASS | `ecology/asexual/minted==recycled` | {'sp1': (1317, 1317), 'sp2': (1408, 1408), 'sp3': (1361, 1361), 'sp4': (1192, 1192)} |
| PASS | `ecology/sexual/invariants` | 30 steps, mainland + islands, every step: phenotype == genotype, sp3 all allele 1, species <-> factory; (t, alive, allele-0 freq) [(0, 400, {'sp1': [0.15, 0.44, 0.59, 0.77, 0.55], 'sp2': [0.48, 0.535, 0.48, 0.46]}), (10, 447, {'sp1': [0.134, 0.361, 0.69, 0.87, 0.426], 'sp2': [0.576, 0.602, 0.555, 0.475]}), (20, 594, {'sp1': [0.137, 0.282, 0.786, 0.878, 0.477], 'sp2': [0.534, 0.568, 0.601, 0.382]}), (30, 790, {'sp1': [0.085, 0.32, 0.797, 0.882, 0.392], 'sp2': [0.592, 0.52, 0.721, 0.379]})] |
| PASS | `ecology/sexual/minted==recycled` | {'sp1': (1105, 1105), 'sp2': (1102, 1102), 'sp3': (1066, 1066), 'sp4': (739, 739)} |
| PASS | `ecology/mix/invariants` | 30 steps, mainland + islands, every step: phenotype == genotype, sp3 all allele 1, species <-> factory; (t, alive, allele-0 freq) [(0, 400, {'sp1': [0.15, 0.44, 0.59, 0.77, 0.55], 'sp2': [0.48, 0.535, 0.48, 0.46]}), (10, 427, {'sp1': [0.13, 0.491, 0.63, 0.801, 0.514], 'sp2': [0.59, 0.557, 0.457, 0.481]}), (20, 625, {'sp1': [0.1, 0.476, 0.664, 0.896, 0.336], 'sp2': [0.563, 0.453, 0.503, 0.514]}), (30, 736, {'sp1': [0.124, 0.456, 0.712, 0.952, 0.308], 'sp2': [0.587, 0.437, 0.561, 0.469]})] |
| PASS | `ecology/mix/minted==recycled` | {'sp1': (980, 980), 'sp2': (1113, 1113), 'sp3': (1045, 1045), 'sp4': (879, 879)} |
| NOTE | `ecology/note/meta_initialize-mix` | reproduce_mode="mix" accepted (400 founders) |
| PASS | `ecology/repro/cross-process` | sexual life cycle, 15 steps, same seeds, two processes: md5 01eb54e9764a / 01eb54e9764a |
| PASS | `ecology/regression/HEAD==tree` | no species with quanti_init_freq, sexual life cycle 20 steps (mainland + islands): HEAD md5 20f1feaa7b3c \| tree md5 20f1feaa7b3c |

PASS 13  FAIL 0  NOTE 1  (12s)  indfactory.py md5 039174c4c2ad2baa435e9a0f0d1d9d87
