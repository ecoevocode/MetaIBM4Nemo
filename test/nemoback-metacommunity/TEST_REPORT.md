# metacommunity.py test report (NemoBackend, sexual)

Generated 2026-10-19 14:11 by `test/nemoback-metacommunity/test_metacommunity.py all`. Scope: `metaibm/metacommunity.py` and its call chain metacommunity -> patch -> habitat -> indfactory -> individual (asex / mix and the extension offspring network only probed).

## Environment

- python `3.11.17`, pynemo `2.4.1`, numpy `2.4.6`, `macOS-15.7.3-arm64-arm-64bit`
- HEAD `7a1adbb` (`MIGRATE patch: NemoBackend (sexual) + test/nemoback-patch`), branch `features/developing`

| file | md5 | state |
|---|---|---|
| `metaibm/metacommunity.py` | `5a3b0ba84fb03bdd0e7b156af6897582` | modified (working tree) |
| `metaibm/patch.py` | `29e70117ed847b7b66eacadc77bf84cf` | committed |
| `metaibm/habitat.py` | `cc2c2001c092d149325b6d1f42e4d2a7` | committed |
| `metaibm/indfactory.py` | `eacaf91df4dbf3ab7dc72100323355be` | committed |
| `metaibm/individual.py` | `2914e40044bb48e0b10364227282d915` | committed |
| `extension/global_habitat_network.py` | `e9a76892a27f0e2ef95e10d9ecd6afd4` | committed |

- original reference (auto-located: newest commit whose metacommunity.py still calls `patch_clear_up_offspring_and_immigrant_pool`): `git show 7a1adbb:metaibm/metacommunity.py` (blob `c42246d`). The last commit that changed metacommunity.py is `ef779f8`, so `7a1adbb:metacommunity.py` == `ef779f8:metacommunity.py`. Original bodies are compiled into a subclass of the CURRENT metacommunity (they run on the current NemoBackend patch / habitat); `meta_initialize` and the two combined clear-ups (removed patch APIs) are never compiled.

## Step 0 -- baseline (existing suites, run unchanged)

| suite | script | result | time |
|---|---|---|---|
| test/nemo-indfactory | accept_genetics.py | PASS | 0.6 s |
| test/nemo-indfactory | accept_guards.py | PASS | 0.6 s |
| test/nemo-indfactory | accept_lifecycle.py | PASS | 0.6 s |
| test/nemo-indfactory | smoke.py | PASS | 0.6 s |
| test/nemoback-habitat (run_all.py) | smoke.py | PASS |  |
| test/nemoback-habitat (run_all.py) | test_bookkeeping.py | PASS | 6.2 s |
| test/nemoback-habitat (run_all.py) | test_nemo.py | PASS | 15.5 s |
| test/nemoback-habitat | **run_all.py total** | PASS | 22.5 s |
| test/nemoback-patch (run_all.py --no-baseline) | smoke.py | PASS |  |
| test/nemoback-patch (run_all.py --no-baseline) | test_bookkeeping.py | PASS | 4.0 s |
| test/nemoback-patch (run_all.py --no-baseline) | test_nemo.py | PASS | 14.2 s |
| test/nemoback-patch (run_all.py --no-baseline) | test_chain.py | **FAIL** | 0.7 s |
| test/nemoback-patch | **run_all.py total** | **FAIL** | 20.9 s |

- both folders were backed up before and restored after their run (their run_all.py rewrites TEST_REPORT.md + results/): restored byte-identical (md5 of every file + git status) = habitat `True`, patch `True`; their fresh reports are kept as `results/baseline_habitat_TEST_REPORT.md` / `results/baseline_patch_TEST_REPORT.md`.

`test/nemoback-patch/test_chain.py` failed; last lines:
```
4 metaibm/simulator.py:95            build_empty_metacommunity_from_patch_habitat_csv           add_habitat                                                              no global_index passed -> every habitat gets homepatch 0
Traceback (most recent call last):
  File "/Users/jianhaolin/Library/CloudStorage/OneDrive-UniversityofHelsinki/MetaIBM4Nemo_workspace/MetaIBM4Nemo-features-developing/test/nemoback-patch/test_chain.py", line 279, in <module>
    assert removed == {'patch_clear_up_offspring_and_immigrant_pool', 'patch_clear_up_offspring_marker_and_immigrant_marker_pool'}, removed
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError: set()
```

## Smoke: PASS (3.6 s)

| step | return | invariants |
|---|---|---|
| meta_initialize (mainland + islands) | `mainland initialization done!  islands initialization done!` | True |
| pool: mainland.meta_dead_selection | `[Dead selection] in mainland: there are 253 individuals dead in selection; there are 1547 individuals in the m` | True |
| pool: mainland.meta_mainland_sexual_birth_mutate_germinate | `[Birth process] in mainland: there are 253 individuals germinating from local habitat; there are 1800 individu` | True |
| pool: islands.meta_colonize_from_propagules_rains | `[Colonization process]: there are 0 individuals colonizing the metacommunity from mainland; there are 150 indi` | True |
| pool: islands.meta_dead_selection | `[Dead selection] in islands: there are 18 individuals dead in selection; there are 132 individuals in the isla` | True |
| pool: islands.meta_sex_reproduce_mutate_into_offspring_pool | `islands: there are 60 individuals born into the offspring_pool; there are 60 individuals in the offspring_pool` | True |
| pool: islands.meta_dispersal_within_patch_from_offspring_to_immigrant_pool | `[Dispersal within process] islands: there are 17 individuals disperse into habs_immigrant_pool within patches;` | True |
| pool: islands.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool | `[Dispersal among patches] in islands: there are 9 individuals disperse into habs_immigrant_pool among patches;` | True |
| pool: islands.meta_local_germinate_from_offspring_and_immigrant_pool | `[Germination process] in islands: there are 18 individuals germinating from local offspring_pool and immigrant` | True |
| pool: islands.meta_dormancy_process_from_offspring_pool_and_immigrant_pool | `[Dormancy process] in islands: there are 0 dormancy survived in the dormancy pool; there are 0 dormancy elimin` | True |
| pool: islands.meta_disturbance_process_in_habitat | `[Disturbance process] occurred in []` | True |
| pool: islands.meta_clear_up_offspring_and_immigrant_pool | `0` | True |
| marker: mainland.meta_dead_selection | `[Dead selection] in mainland: there are 236 individuals dead in selection; there are 1564 individuals in the m` | True |
| marker: mainland.meta_mainland_sexual_birth_mutate_germinate | `[Birth process] in mainland: there are 236 individuals germinating from local habitat; there are 1800 individu` | True |
| marker: islands.meta_dead_selection | `[Dead selection] in islands: there are 25 individuals dead in selection; there are 125 individuals in the isla` | True |
| marker: islands.meta_sex_reproduce_calculation_with_offspring_marker_pool | `[Reproduction into offspring_marker_pool] in islands: there are 40 offspring_marker born into the offspring_ma` | True |
| marker: islands.meta_colonize_from_propagules_rains | `[Colonization process]: there are 2 individuals colonizing the metacommunity from mainland; there are 127 indi` | True |
| marker: islands.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool | `[Dispersal within process] in islands: there are 11 offspring marker disperse into habs_immigrant_marker_pool ` | True |
| marker: islands.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool | `[Dispersal among patches] in islands: there are 6 individuals disperse into habs_immigrant_marker_pool among p` | True |
| marker: islands.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool | `[Germination & Birth process] islands: there are 23 individuals germinating_and_birth from local offspring_mar` | True |
| marker: islands.meta_disturbance_process_in_patches | `[Disturbance process] occurred in ['p2']` | True |
| marker: islands.meta_clear_up_offspring_marker_and_immigrant_marker_pool | `12` | True |
| teardown (meta_clear_up_* on both, minted == recycled) | `1843` | True |
| two fresh processes, seed 7 | `state 9e2c9cff090cf7af, rng e6daecc36f3d2f69, 12 per-step digests identical` | True |

## Tests: 3 / 3 PASS

| script | covers | result | time |
|---|---|---|---|
| bookkeeping | R1 is_remove defaults, R3 else branches, R2 is_remove=False guards (state + both RNGs + counters + currentID), R4 marker modes, clear-ups, meta_initialize, R6 tag audit + docstring rule, patch_num < 2 logs, unchanged helpers | PASS | 2.1 s |
| nemo | one fresh process per section: equivalence with the original bodies (one worker process per function and variant) / colonization statistics / Nemo level / edges / probes | PASS | 29.5 s |
| chain | long runs (3 schedules x 3 seeds x 200 steps, invariants after every call, RSS, C++ addresses, teardown) + fuzz (>= 3000 random operations) | PASS | 20.9 s |

## bookkeeping: PASS

| check | result | numbers |
|---|---|---|
| R1 every is_remove default True | PASS | 14 functions (original defaults: False) |
| R3 guard first; unreachable else == `pass  # ... TODO(copy_from)` | PASS | 5 guarded meta functions, 6 else branches |
| R2 is_remove=False -> NotImplementedError before any RNG use / state change | PASS | 352 calls over 11 paths (3- and 1-patch islands); raised in: habitat layer, meta layer, patch layer |
| R4 marker paths with is_remove False and True | PASS | {'within True': 15, 'within False': 15, 'among True': 15, 'among False': 15, 'birth True': 15, 'birth False': 15}, 764 births |
| clear-ups: 6 split + 2 wrappers, int returns, dormancy kept by O+I wrapper | PASS | 64 cases |
| meta_initialize fills every site with the right species, homepatch, re-init | PASS | 90 founders; bad inputs (info): {'mean (0.9, 0.9) in p2.h1': 'ValueError', "indfactory 'sp2' missing (first sp2 habitat p0.h1)": 'KeyError'} |
| R6 NemoBackend tags == functions reaching recycle / make_founder / make_offspring / tagged methods | PASS | 16 tagged: dispersal_among_patches_from_offspring_pool_and_dormancy_pool, meta_clear_up_dormancy_pool, meta_clear_up_immigrant_pool, meta_clear_up_microsite_individuals, meta_clear_up_offspring_and_immigrant_pool, meta_clear_up_offspring_pool, meta_dead_selection, meta_dispersal_within_patch_from_offspring_and_dormancy_pool, meta_disturbance_process_in_habitat, meta_disturbance_process_in_patches, meta_dormancy_process_from_offspring_pool_and_immigrant_pool, meta_dormancy_process_from_offspring_pool_to_dormancy_pool, meta_initialize, meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool, meta_mainland_sexual_birth_mutate_germinate, meta_sex_reproduce_mutate_into_offspring_pool |
| R6 untagged functions keep the original docstring; new ones one-line | PASS | 59 kept, 0 new one-line () |
| patch_num < 2 log strings (info) | PASS | 4 functions now format the metacommunity name |
| unchanged getters / matrix helpers / CSV writer == original bodies (same seed) | PASS | 34 functions x 6 seeds |

<details><summary>R6 tag audit: why each tagged function is NemoBackend</summary>

| function | reason (calls) |
|---|---|
| dispersal_among_patches_from_offspring_pool_and_dormancy_pool | recycle |
| meta_clear_up_dormancy_pool | patch_clear_up_dormancy_pool |
| meta_clear_up_immigrant_pool | patch_clear_up_immigrant_pool |
| meta_clear_up_microsite_individuals | patch_clear_up_microsite_individuals |
| meta_clear_up_offspring_and_immigrant_pool | self.meta_clear_up_immigrant_pool, self.meta_clear_up_offspring_pool |
| meta_clear_up_offspring_pool | patch_clear_up_offspring_pool |
| meta_dead_selection | patch_dead_selection |
| meta_dispersal_within_patch_from_offspring_and_dormancy_pool | patch_dipersal_within_from_offspring_and_dormancy_pool |
| meta_disturbance_process_in_habitat | habitat_disturbance_process |
| meta_disturbance_process_in_patches | patch_disturbance_process |
| meta_dormancy_process_from_offspring_pool_and_immigrant_pool | patch_dormancy_process_from_offspring_pool_and_immigrant_pool |
| meta_dormancy_process_from_offspring_pool_to_dormancy_pool | patch_dormancy_process_from_offspring_pool_to_dormancy_pool |
| meta_initialize | patch_initialize |
| meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool | hab_sex_reproduce_mutate_with_num |
| meta_mainland_sexual_birth_mutate_germinate | patch_sexual_birth_germinate |
| meta_sex_reproduce_mutate_into_offspring_pool | patch_sex_reproduce_mutate_into_offspring_pool |

</details>

## nemo: PASS

### Equivalence with the original bodies (one worker process per function and variant, same seed + same world)

| key | seeds | state | random + np.random | return | left all owners: same IDs | new: all left recycled | extra |
|---|---|---|---|---|---|---|---|
| among OD->sites | 120 | 120 | 120 | 120 | 120 | 120 | recycled == dropped by v3 on 120 seeds, 36 extras (2 natural) |
| among O->I | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| among marker True | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| among marker False | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| among OD->I | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| within marker True | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| within marker False | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| within O->I | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| within OD->I | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| within OD->sites | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| germ O+D | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| germ O+I | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| germ O+I+D | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| germ+birth markers True | 30 | 30 | 30 | 30 | 30 | 30 | births 289, recycled 0 |
| germ+birth markers False | 30 | 30 | 30 | 30 | 30 | 30 | births 289, recycled 0 |
| dead selection | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 357 |
| dormancy O->D | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 416 |
| dormancy O+I | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 279 |
| sex reproduce | 30 | 30 | 30 | 30 | 30 | 30 | births 229, recycled 364 |
| mainland sexual birth | 30 | 30 | 30 | 30 | 30 | 30 | births 779, recycled 0 |
| marker calc sex | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| marker calc asex | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| marker calc mix | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 0 |
| disturbance patches | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 1039 |
| disturbance habitat | 30 | 30 | 30 | 30 | 30 | 30 | births 0, recycled 1092 |

### Colonization (statistical)

| check | result | numbers |
|---|---|---|
| count = floor(n) + Bernoulli(frac), n=1.8 | PASS | 1500 seeds, counts {'2': 1185, '1': 315}, mean 1.79, P(extra) 0.79 vs 0.8, binomial p = 0.3329 |
| count = floor(n) + Bernoulli(frac), n=6.4 | PASS | 1500 seeds, counts {'6': 893, '7': 607}, mean 6.4047, P(extra) 0.4047 vs 0.4, binomial p = 0.7122 |
| count = floor(n) + Bernoulli(frac), n=0.3 | PASS | 1500 seeds, counts {'0': 1047, '1': 453}, mean 0.302, P(extra) 0.302 vs 0.3, binomial p = 0.8658 |
| each mainland individual equally likely | PASS | 2000 seeds, 12788 picks / 72 individuals, chi2 = 44.2 (df 71), p = 0.9948 |
| each empty island site equally likely | PASS | 2000 seeds, 12774 placements / 24 sites, chi2 = 19.0 (df 23), p = 0.6997 |
| pairwise: (female, male), one species, one mainland habitat -> two empty sites of one island habitat | PASS | 1500 seeds, pairs {'4': 324, '3': 1176}, P(4 pairs) = 0.216 vs 0.2, p = 0.1213; pairs by mainland habitat {'h1': 2419, 'h0': 2405} |
| mainland loses exactly what islands gain; one np.random.uniform per call | PASS | checked on every call above |
| capacity-limited | PASS | {'single': [2], 'pairwise': [0]} |

### Nemo level

| check | result | numbers |
|---|---|---|
| home_habitat == birth habitat global_index | PASS | offspring pool 190, mainland birth 2714, marker birth 155 (117 in another habitat, 41 in another patch than the germination site) |
| marker birth: newborns per (germination, birth habitat) == markers consumed | PASS | every marker-birth call |
| recycle() only on an individual owned by nobody; recycled never read again | PASS | 3205 recycle() calls |
| minted - recycled == held per factory; teardown minted == recycled | PASS | {'sp1': [1715, 1715], 'sp2': [1782, 1782]} |

### Edges (pre-existing, recorded)

| case | result |
|---|---|
| rain 30 > 18 mainland individuals | ValueError: Sample larger than population or is negative, state unchanged |
| pairwise rain 30 > 2 x mainland pairs | ValueError: Sample larger than population or is negative, state unchanged |
| empty mainland, rain 0.4 (200 seeds) | ValueError in 80 / 200 calls (expected ~80), state unchanged |
| stale markers after a disturbance (30 seeds) | {'ValueError': 30, 'ok': 0} |

### Probes

| probe | result |
|---|---|
| meta_mainland_asexual_birth_mutate_germinate | TypeError: cannot pickle 'swig_runtime_data5.SwigPyObject' object |
| meta_mainland_mixed_birth_mutate_germinate | TypeError: cannot pickle 'swig_runtime_data5.SwigPyObject' object |
| meta_asex_reproduce_mutate_into_offspring_pool | TypeError: cannot pickle 'swig_runtime_data5.SwigPyObject' object |
| meta_mix_reproduce_mutate_into_offspring_pool | TypeError: cannot pickle 'swig_runtime_data5.SwigPyObject' object |
| meta_asex_reproduce_calculation_into_offspring_marker_pool | returned '[Reproduction into offspring_marker_pool] in islands: there ' |
| meta_mix_reproduce_calculation_with_offspring_marker_pool | returned '[Reproduction into offspring_marker_pool] in islands: there ' |
| germinate-and-birth with asexual marker | {'raised': 19, 'born_before_raise': 24, 'seeds_with_partial_births': 12} |
| germinate-and-birth with mix_asexual marker | {'raised': 19, 'born_before_raise': 39, 'seeds_with_partial_births': 16} |
| germinate-and-birth with mix_sexual marker | {'raised': 19, 'born_before_raise': 33, 'seeds_with_partial_births': 15} |
| extension: offspring_network_two_owner_individuals | 14 |
| extension: check_world | invariant broken: extension (a) an object is held twice |
| extension: marker_network (source before, source after, immigrant markers) | [36, 36, 20] |

## chain: PASS

| run | steps | calls (invariants after each) | time | births | islands alive mean / capacity | RSS MB (t=20 -> 199) | distinct C++ addresses / alive high-water mark | teardown minted = recycled |
|---|---|---|---|---|---|---|---|---|
| pool+dormancy seed=1 | 200 | 2200 | 11 s | 65288 | 266.2 / 288 (min 186) | 131.7 -> 131.8 -> 131.8 -> 131.8 | sp1 1135 / 1135, sp2 1136 / 1136 | sp1 32620 = 32620, sp2 32668 = 32668 |
| pool+dormancy seed=2 | 200 | 2200 | 11 s | 65711 | 272.4 / 288 (min 144) | 132.6 -> 132.7 -> 132.8 -> 132.8 | sp1 1135 / 1135, sp2 1143 / 1143 | sp1 32652 = 32652, sp2 33059 = 33059 |
| pool+dormancy seed=3 | 200 | 2200 | 11 s | 65213 | 263.9 / 288 (min 178) | 132.4 -> 132.5 -> 132.6 -> 132.6 | sp1 1127 / 1127, sp2 1148 / 1148 | sp1 32590 = 32590, sp2 32623 = 32623 |
| pool+clear-ups seed=1 | 200 | 2000 | 10 s | 64327 | 244.6 / 288 (min 163) | 132.2 -> 132.3 -> 132.3 -> 132.3 | sp1 1093 / 1093, sp2 1086 / 1086 | sp1 32221 = 32221, sp2 32106 = 32106 |
| pool+clear-ups seed=2 | 200 | 2000 | 10 s | 61643 | 206.8 / 288 (min 67) | 131.8 -> 131.8 -> 131.8 -> 131.9 | sp1 1089 / 1089, sp2 1091 / 1091 | sp1 30440 = 30440, sp2 31203 = 31203 |
| pool+clear-ups seed=3 | 200 | 2000 | 10 s | 60346 | 184.4 / 288 (min 68) | 131.7 -> 131.7 -> 131.8 -> 131.8 | sp1 1090 / 1090, sp2 1079 / 1079 | sp1 31437 = 31437, sp2 28909 = 28909 |
| marker seed=1 | 200 | 2000 | 8 s | 57174 | 244.9 / 288 (min 151) | 130.6 -> 130.8 -> 130.8 -> 130.8 | sp1 1070 / 1070, sp2 1057 / 1057 | sp1 28908 = 28908, sp2 28266 = 28266 |
| marker seed=2 | 200 | 2000 | 8 s | 57477 | 245.8 / 288 (min 148) | 131.0 -> 131.1 -> 131.2 -> 131.2 | sp1 1069 / 1069, sp2 1060 / 1060 | sp1 28842 = 28842, sp2 28635 = 28635 |
| marker seed=3 | 200 | 2000 | 8 s | 56748 | 217.3 / 288 (min 88) | 130.9 -> 131.0 -> 131.0 -> 131.0 | sp1 1077 / 1077, sp2 1053 / 1053 | sp1 28392 = 28392, sp2 28356 = 28356 |

- fuzz: 3300 random operations (36 kinds), invariants after each, minted == recycled after every run; allowed pre-existing exceptions {'colonize: rain > mainland population': 35, 'germinate-and-birth: stale marker': 1}
- fuzz operations: among_O->I 79, among_OD->I 87, among_OD->sites 89, among_marker 88, clear_D 83, clear_I 96, clear_IM 90, clear_O 106, clear_O+I 92, clear_OM 98, clear_OM+IM 98, clear_micro 84, colonize 104, colonize_pairwise 89, dead 106, dist_habitat 102, dist_patches 85, dorm_O+I 92, dorm_O->D 92, germ+birth 80, germ_O+D 89, germ_O+I 85, germ_O+I+D 101, getters 82, guard is_remove=False 82, initialize 106, marker_calc 93, ml_birth 96, ml_dead 80, ml_initialize 100, network_marker 83, reproduce 104, within_O->I 77, within_OD->I 90, within_OD->sites 101, within_marker 91

<details><summary>bookkeeping output</summary>

```
1 R1: 14 metacommunity functions take is_remove, every default is True (original: ['False'])
2 R3: the 5 meta-layer guards are the first statement (after the docstring); their 6 `else` branches are `pass  # unreachable ... TODO(copy_from)`; marker paths keep a live is_remove=False branch
3 R2: 352 guarded calls (11 paths x positional + keyword, 3- and 1-patch islands); NotImplementedError before any change: state digest, random + np.random state, minted / recycled, Nemo currentID all unchanged
4 R4: marker within / among / germinate-and-birth run with is_remove=True (markers move; births consume them) and False (source kept, copies sent; markers kept after birth): {'within True': 15, 'within False': 15, 'among True': 15, 'among False': 15, 'birth True': 15, 'birth False': 15}, 764 births
5 clear-ups: 64 cases over the 6 split functions + 2 wrappers; return == number cleared (int), recycled == former contents, only those pools emptied, meta_clear_up_offspring_and_immigrant_pool keeps the dormancy pool (identity), marker clear-ups recycle nothing, second call returns 0
6 meta_initialize(indfactories, species_2_phenotype_ls, reproduce_mode="sexual"): every site filled with the species of its habitat mean, getHome == global_index, re-initialize recycles residents and keeps pools, asexual -> all female; bad inputs (info): {'mean (0.9, 0.9) in p2.h1': ('ValueError', 45), "indfactory 'sp2' missing (first sp2 habitat p0.h1)": ('KeyError', 9)}
7 R6: 16 NemoBackend-tagged functions == exactly those calling recycle / make_founder / make_offspring or a tagged patch / habitat method (or a tagged metacommunity wrapper); 59 untagged functions keep the original docstring, 0 new untagged ones are one-line
8 patch_num < 2: 4 log strings now contain the metacommunity name (v3: literal "%s"), e.g. '[Dispersal among patches] in solo: patch_num < 2, there are 0 individuals disperse among patches'
9 unchanged getters / matrix helpers / data IO (34 functions incl. meta_distribution_data_all_time_to_csv_gz to a temp dir): same seed -> identical values to the original bodies on 6 seeds
PASSED
```
</details>
<details><summary>nemo output</summary>

```
==== section equivalence
1 among OD->sites            dispersal_among_patches_from_offspring_pool_and_dormancy_pool                    seeds 120  state=120 rng=120 return=120 left=120 OK
1 among O->I                 dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool                    seeds  30  state=30 rng=30 return=30 left=30 OK
1 among marker True          dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool      seeds  30  state=30 rng=30 return=30 left=30 OK
1 among marker False         dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool      seeds  30  state=30 rng=30 return=30 left=30 OK
1 among OD->I                dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool  seeds  30  state=30 rng=30 return=30 left=30 OK
1 within marker True         meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool       seeds  30  state=30 rng=30 return=30 left=30 OK
1 within marker False        meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool       seeds  30  state=30 rng=30 return=30 left=30 OK
1 within O->I                meta_dispersal_within_patch_from_offspring_to_immigrant_pool                     seeds  30  state=30 rng=30 return=30 left=30 OK
1 within OD->I               meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool        seeds  30  state=30 rng=30 return=30 left=30 OK
1 within OD->sites           meta_dispersal_within_patch_from_offspring_and_dormancy_pool                     seeds  30  state=30 rng=30 return=30 left=30 OK
1 germ O+D                   meta_local_germinate_from_offspring_and_dormancy_pool                            seeds  30  state=30 rng=30 return=30 left=30 OK
1 germ O+I                   meta_local_germinate_from_offspring_and_immigrant_pool                           seeds  30  state=30 rng=30 return=30 left=30 OK
1 germ O+I+D                 meta_local_germinate_from_offspring_immigrant_and_dormancy_pool                  seeds  30  state=30 rng=30 return=30 left=30 OK
1 germ+birth markers True    meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool   seeds  30  state=30 rng=30 return=30 left=30 OK
1 germ+birth markers False   meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool   seeds  30  state=30 rng=30 return=30 left=30 OK
1 dead selection             meta_dead_selection                                                              seeds  30  state=30 rng=30 return=30 left=30 OK
1 dormancy O->D              meta_dormancy_process_from_offspring_pool_to_dormancy_pool                       seeds  30  state=30 rng=30 return=30 left=30 OK
1 dormancy O+I               meta_dormancy_process_from_offspring_pool_and_immigrant_pool                     seeds  30  state=30 rng=30 return=30 left=30 OK
1 sex reproduce              meta_sex_reproduce_mutate_into_offspring_pool                                    seeds  30  state=30 rng=30 return=30 left=30 OK
1 mainland sexual birth      meta_mainland_sexual_birth_mutate_germinate                                      seeds  30  state=30 rng=30 return=30 left=30 OK
1 marker calc sex            meta_sex_reproduce_calculation_with_offspring_marker_pool                        seeds  30  state=30 rng=30 return=30 left=30 OK
1 marker calc asex           meta_asex_reproduce_calculation_into_offspring_marker_pool                       seeds  30  state=30 rng=30 return=30 left=30 OK
1 marker calc mix            meta_mix_reproduce_calculation_with_offspring_marker_pool                        seeds  30  state=30 rng=30 return=30 left=30 OK
1 disturbance patches        meta_disturbance_process_in_patches                                              seeds  30  state=30 rng=30 return=30 left=30 OK
1 disturbance habitat        meta_disturbance_process_in_habitat                                              seeds  30  state=30 rng=30 return=30 left=30 OK
1 equivalence: 25 functions / modes, every seed identical to 7a1adbb (state digest, random + np.random state, return value, who left every owner); among OD->sites: recycled == dropped by v3 (36 extras)
---- section equivalence PASSED
==== section colonization
2a count: floor(n) + Bernoulli(frac) over 1500 seeds each: {'n=1.8': (0.79, 0.8, 0.3329), 'n=6.4': (0.4047, 0.4, 0.7122), 'n=0.3': (0.302, 0.3, 0.8658)}
2b mainland pick uniform: 12788 picks over 72 individuals, chi2 = 44.2 (df 71), p = 0.995
2c empty island site uniform: 12774 placements over 24 empty sites, chi2 = 19.0 (df 23), p = 0.700
2d pairwise over 1500 seeds: pairs {4: 324, 3: 1176} (P(4) = 0.216, expected 0.2, p = 0.121), every pair (female, male) of one species from one mainland habitat into two empty sites of one island habitat; colonists keep their mainland homepatch
2e capacity: 2 empty sites in two habitats -> single colonization places exactly 2, pairwise places 0 (no habitat with 2 empty sites); mainland loses only what lands
---- section colonization PASSED
==== section nemo_level
3 nemo level: homepatch == birth habitat global_index for 190 offspring-pool births, 2714 mainland births, 155 marker births (117 born in another habitat than they germinate in, 41 in another patch): home = marker habitat, newborns per (germination, birth habitat) == markers consumed
3 recycle: 3205 recycle() calls, each on an individual owned by nobody; recycled individuals raise on every read; minted - recycled == held after every step; teardown {'sp1': (1715, 1715), 'sp2': (1782, 1782)}
---- section nemo_level PASSED
==== section edges
4 edges: colonization {'rain 30 > 18 mainland individuals': ('ValueError: Sample larger than population or is negative', 'state unchanged'), 'pairwise rain 30 > 2 x mainland pairs': ('ValueError: Sample larger than population or is negative', 'state unchanged'), 'empty mainland, rain 0.4 (200 seeds)': ('ValueError in 80 / 200 calls (expected ~80)', 'state unchanged')}; stale markers after a disturbance -> {'ValueError': 30, 'ok': 0}
---- section edges PASSED
==== section probes
5.1 asex / mix probes: ['meta_mainland_asexual_birth_mutate_germinate: TypeError', 'meta_mainland_mixed_birth_mutate_germinate: TypeError', 'meta_asex_reproduce_mutate_into_offspring_pool: TypeError', 'meta_mix_reproduce_mutate_into_offspring_pool: TypeError', 'meta_asex_reproduce_calculation_into_offspring_marker_pool: returned', 'meta_mix_reproduce_calculation_with_offspring_marker_pool: returned']
5.2 asexual / mix markers in germinate-and-birth: {'asexual': {'raised': 19, 'born_before_raise': 24, 'seeds_with_partial_births': 12}, 'mix_asexual': {'raised': 19, 'born_before_raise': 39, 'seeds_with_partial_births': 16}, 'mix_sexual': {'raised': 19, 'born_before_raise': 33, 'seeds_with_partial_births': 15}}
5.3 extension offspring network: 14 individuals end up with two owners (invariant broken: extension (a) an object is held twice); marker network copies markers (36, 36, 20)
---- section probes PASSED
PASSED
```
</details>
<details><summary>chain output</summary>

```
1 pool+dormancy   seed=1: 2200 calls (invariants after each) in 11s, births 65288, islands alive mean 266 / 288 (min 186), RSS {'20': 131.7, '100': 131.8, '150': 131.8, '199': 131.8} MB, distinct addresses vs alive high-water mark {'sp1': [1135, 1135], 'sp2': [1136, 1136]}, teardown {'sp1': [32620, 32620], 'sp2': [32668, 32668]}
1 pool+dormancy   seed=2: 2200 calls (invariants after each) in 11s, births 65711, islands alive mean 272 / 288 (min 144), RSS {'20': 132.6, '100': 132.7, '150': 132.8, '199': 132.8} MB, distinct addresses vs alive high-water mark {'sp1': [1135, 1135], 'sp2': [1143, 1143]}, teardown {'sp1': [32652, 32652], 'sp2': [33059, 33059]}
1 pool+dormancy   seed=3: 2200 calls (invariants after each) in 11s, births 65213, islands alive mean 264 / 288 (min 178), RSS {'20': 132.4, '100': 132.5, '150': 132.6, '199': 132.6} MB, distinct addresses vs alive high-water mark {'sp1': [1127, 1127], 'sp2': [1148, 1148]}, teardown {'sp1': [32590, 32590], 'sp2': [32623, 32623]}
1 pool+clear-ups  seed=1: 2000 calls (invariants after each) in 10s, births 64327, islands alive mean 245 / 288 (min 163), RSS {'20': 132.2, '100': 132.3, '150': 132.3, '199': 132.3} MB, distinct addresses vs alive high-water mark {'sp1': [1093, 1093], 'sp2': [1086, 1086]}, teardown {'sp1': [32221, 32221], 'sp2': [32106, 32106]}
1 pool+clear-ups  seed=2: 2000 calls (invariants after each) in 10s, births 61643, islands alive mean 207 / 288 (min 67), RSS {'20': 131.8, '100': 131.8, '150': 131.8, '199': 131.9} MB, distinct addresses vs alive high-water mark {'sp1': [1089, 1089], 'sp2': [1091, 1091]}, teardown {'sp1': [30440, 30440], 'sp2': [31203, 31203]}
1 pool+clear-ups  seed=3: 2000 calls (invariants after each) in 10s, births 60346, islands alive mean 184 / 288 (min 68), RSS {'20': 131.7, '100': 131.7, '150': 131.8, '199': 131.8} MB, distinct addresses vs alive high-water mark {'sp1': [1090, 1090], 'sp2': [1079, 1079]}, teardown {'sp1': [31437, 31437], 'sp2': [28909, 28909]}
1 marker          seed=1: 2000 calls (invariants after each) in 8s, births 57174, islands alive mean 245 / 288 (min 151), RSS {'20': 130.6, '100': 130.8, '150': 130.8, '199': 130.8} MB, distinct addresses vs alive high-water mark {'sp1': [1070, 1070], 'sp2': [1057, 1057]}, teardown {'sp1': [28908, 28908], 'sp2': [28266, 28266]}
1 marker          seed=2: 2000 calls (invariants after each) in 8s, births 57477, islands alive mean 246 / 288 (min 148), RSS {'20': 131.0, '100': 131.1, '150': 131.2, '199': 131.2} MB, distinct addresses vs alive high-water mark {'sp1': [1069, 1069], 'sp2': [1060, 1060]}, teardown {'sp1': [28842, 28842], 'sp2': [28635, 28635]}
1 marker          seed=3: 2000 calls (invariants after each) in 8s, births 56748, islands alive mean 217 / 288 (min 88), RSS {'20': 130.9, '100': 131.0, '150': 131.0, '199': 131.0} MB, distinct addresses vs alive high-water mark {'sp1': [1077, 1077], 'sp2': [1053, 1053]}, teardown {'sp1': [28392, 28392], 'sp2': [28356, 28356]}
2 fuzz: 3300 random operations (36 kinds) over every migrated metacommunity method + guarded is_remove=False calls + getters, invariants after each, minted == recycled after every run; allowed pre-existing exceptions {'colonize: rain > mainland population': 35, 'germinate-and-birth: stale marker': 1}
PASSED
```
</details>

## Findings (13: bug 0, risk 1, info 8, known 4)

Severity: **bug** = wrong behaviour in the code under test; **risk** = correct today, fragile / surprising; **info** = observation or intended behaviour change vs v3; **known** = expected / next migration round / pre-existing. Suggested fixes are suggestions only -- nothing was applied.

| severity | function | finding | evidence | suggested fix (not applied) |
|---|---|---|---|---|
| risk | meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool | a marker whose birth habitat has no female-male pair any more raises ValueError (random.sample of an empty pair list in hab_sex_reproduce_mutate_with_num) part-way through the loop: earlier markers are already born. Not reachable in the experiments/model-sloss*.py order (disturbance after germination), reachable if a disturbance / clear-up runs between marker dispersal and germination; pre-existing in v3. Invariants hold afterwards. | marker calc -> among-patch marker dispersal (m=0.8) -> disturbance of p0 -> germinate-and-birth: ValueError in 30 / 30 seeds | skip (or drop) a marker whose birth habitat has no parent pair, or clear the markers of a disturbed habitat everywhere |
| info | test/nemoback-patch/test_chain.py (Step 0) | stale since metacommunity.py was migrated (suite not edited): its section 4 static scan asserts that metacommunity.py still calls the REMOVED patch_clear_up_offspring_and_immigrant_pool / patch_clear_up_offspring_marker_and_immigrant_marker_pool and patch_initialize with the old 6-argument signature (the "next migration round" findings of that report); the migration removed exactly those calls, so `removed == set()` and the assertion fails. Sections 1-3 of that script (call chains, ordering rule, homepatch) passed before the scan. Not a code defect. |     assert removed == {'patch_clear_up_offspring_and_immigrant_pool', 'patch_clear_up_offspring_marker_and_immigrant_marker_pool'}, removed /            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ / AssertionError: set() | in test/nemoback-patch/test_chain.py scan(): expect removed == set() and changed == set() once metacommunity.py is migrated (or compare against the git version of metacommunity.py the patch migration was tested with) |
| info | meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool, dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool | marker dispersal default is_remove False -> True: with the defaults, markers now MOVE from the offspring- to the immigrant-marker pool (total conserved); v3 defaults kept every marker in the source pool AND sent a copy (more births). Both modes stay seed-identical to v3 when is_remove is passed explicitly (test_nemo equivalence). | within, 20 seeds x 60 markers (offspring, immigrant) after the call with each version's default: now [600, 600], v3 [1200, 600]; among, 20 seeds x 60 markers (offspring, immigrant) after the call with each version's default: now [600, 600], v3 [1200, 600] |  |
| info | meta_initialize | no up-front validation (non-atomic by design in this version): a bad habitat raises inside patch_initialize after every earlier habitat was already filled; invariants and factory balance still hold. | mean (0.9, 0.9) in p2.h1 -> ValueError after 45 founders, filled ['p0.h0', 'p0.h1', 'p1.h0', 'p1.h1', 'p2.h0']; indfactory 'sp2' missing (first sp2 habitat p0.h1) -> KeyError after 9 founders, filled ['p0.h0'] | validate every habitat (species_2_phenotype_ls.index(mean_env_ls), indfactory key, 0 <= global_index <= 65535, reproduce_mode) before the first patch_initialize |
| info | meta_clear_up_immigrant_marker_pool, meta_clear_up_offspring_marker_pool | docstring rule (untagged new functions: one-line docstring) not met; documentation only, no runtime effect. | {'meta_clear_up_offspring_marker_pool': 'NOT one-line', 'meta_clear_up_immigrant_marker_pool': 'NOT one-line'} | add a one-line docstring, e.g. ''' empty <pool> in every patch (markers are tuples: nothing to recycle) ''' |
| info | dispersal_among_patches_from_offspring_pool_and_dormancy_pool, dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool, dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool, dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool | the 4 `patch_num < 2` early-return log strings now format the metacommunity name (v3 returned the literal "%s"); no RNG used, nothing changed. | [Dispersal among patches] in solo: patch_num < 2, there are 0 individuals disperse among patches |  |
| info | dispersal_among_patches_from_offspring_pool_and_dormancy_pool | migrants left without an empty site are now recycled (v3 removed them from their source pool and dropped them, alive and owned by nobody); otherwise seed-identical to v3. | 120 seeds: state / random / np.random / return identical; recycled set == the migrants v3 dropped on every seed (36 extras, 2 of them without forcing the number matrix) |  |
| info | meta_colonize_from_propagules_rains, pairwise_sexual_colonization_from_prpagules_rains | colonization now MOVES the propagules (is_remove=True default: del_individual on the mainland, no recycle) -- v3 default kept them on the mainland too (two owners); the stochastic rounding is int(n + U) and the mainland shuffle is gone, so runs are not seed-identical to v3; tested statistically instead. | P(extra): {'n=1.8': (0.79, 0.8), 'n=6.4': (0.4047, 0.4), 'n=0.3': (0.302, 0.3)}; mainland pick chi2 p = 0.995; island site chi2 p = 0.700; pairwise P(4 pairs) = 0.216 (p = 0.121); np.random use = exactly one uniform per call |  |
| info | meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool | asexual / mix_asexual / mix_sexual markers raise NotImplementedError inside the reproduce_mode branch (TODO(makeClone) / TODO(mixed)), after earlier markers of the same call were already born (partial state by design); bookkeeping stays consistent. | asexual: raised in 19/20 seeds (the marker is not always reached: fewer empty sites than markers), 24 births before the raise; mix_asexual: raised in 19/20 seeds (the marker is not always reached: fewer empty sites than markers), 39 births before the raise; mix_sexual: raised in 19/20 seeds (the marker is not always reached: fewer empty sites than markers), 33 births before the raise |  |
| known | meta_colonize_from_propagules_rains, pairwise_sexual_colonization_from_prpagules_rains | random.sample raises ValueError when int(n + U) exceeds the mainland population (or its female-male pairs), e.g. an empty mainland with rain 0.4 fails in ~40% of calls; nothing is moved (pre-existing, v3 identical). | rain 30 > 18 mainland individuals -> ValueError: Sample larger than population or is negative, state unchanged; pairwise rain 30 > 2 x mainland pairs -> ValueError: Sample larger than population or is negative, state unchanged; empty mainland, rain 0.4 (200 seeds) -> ValueError in 80 / 200 calls (expected ~80), state unchanged | sample min(k, len(population)) like sample_offspring_without_replacement does |
| known | meta_mainland_asexual_birth_mutate_germinate, meta_mainland_mixed_birth_mutate_germinate, meta_asex_reproduce_mutate_into_offspring_pool, meta_mix_reproduce_mutate_into_offspring_pool, meta_asex_reproduce_calculation_into_offspring_marker_pool, meta_mix_reproduce_calculation_with_offspring_marker_pool | asex / mix meta functions are not migrated: meta_mainland_asexual_birth_mutate_germinate -> TypeError; meta_mainland_mixed_birth_mutate_germinate -> TypeError; meta_asex_reproduce_mutate_into_offspring_pool -> TypeError; meta_mix_reproduce_mutate_into_offspring_pool -> TypeError; meta_asex_reproduce_calculation_into_offspring_marker_pool -> returns (markers only, no individual); meta_mix_reproduce_calculation_with_offspring_marker_pool -> returns (markers only, no individual) (TypeError = copy.deepcopy of a SwigPyObject in the habitat layer). |  |  |
| known | extension/global_habitat_network.py: dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool | not migrated: random.sample(h_i.offspring_pool, n) without removal + immigrant_pool += -> the SAME individual is in the source offspring_pool and a target immigrant_pool (two owners); a later clear-up would recycle it twice. | 14 individuals with two owners after one call; check_world: invariant broken: extension (a) an object is held twice | use patch / habitat sample_offspring_without_replacement (h_i_object.sample_offspring_without_replacement(n)) as the migrated metacommunity functions do |
| known | extension/global_habitat_network.py: dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool | markers are copied, not moved (v3 is_remove=False semantics; tuples, no ownership problem): offspring markers 36 -> 36, immigrant markers +20. |  |  |

## Method coverage

`vs original` compares each method with `7a1adbb:metaibm/metacommunity.py` (AST + source text). `scripts` = test files of this suite that name the method (common = the shared schedules).

| method | vs original | tag | scope | scripts |
|---|---|---|---|---|
| `__init__` | changed |  | in scope | common, test_metacommunity |
| `get_data` | unchanged |  | in scope | test_metacommunity |
| `__str__` | unchanged |  | in scope | test_metacommunity |
| `update_disp_current_matrix` | unchanged |  | in scope | test_metacommunity |
| `add_patch` | unchanged |  | in scope | common |
| `reshape_habitat_data_in_patch` | unchanged |  | in scope | test_metacommunity |
| `meta_show_species_distribution` | unchanged |  | plot (unchanged, not run) | - |
| `meta_show_species_phenotype_distribution` | unchanged |  | plot (unchanged, not run) | - |
| `meta_show_environment_distribution` | unchanged |  | plot (unchanged, not run) | - |
| `meta_show_two_environment_distribution` | unchanged |  | plot (unchanged, not run) | - |
| `meta_initialize` | changed | NemoBackend | in scope | common, test_metacommunity |
| `get_meta_individual_num` | unchanged |  | in scope | test_metacommunity |
| `get_meta_empty_sites_ls` | unchanged |  | in scope | test_metacommunity |
| `show_meta_empty_sites_num` | unchanged |  | in scope | test_metacommunity |
| `get_meta_pairwise_empty_sites_ls` | unchanged |  | in scope | test_metacommunity |
| `meta_get_occupied_location_ls` | unchanged |  | in scope | test_metacommunity |
| `get_meta_pairwise_occupied_sites_ls` | unchanged |  | in scope | test_metacommunity |
| `meta_offspring_pool_individual_num` | unchanged |  | in scope | test_metacommunity |
| `meta_immigrant_pool_individual_num` | unchanged |  | in scope | test_metacommunity |
| `meta_dormancy_pool_individual_num` | unchanged |  | in scope | test_metacommunity |
| `meta_offspring_marker_pool_marker_num` | unchanged |  | in scope | test_metacommunity |
| `meta_immigrant_marker_pool_marker_num` | unchanged |  | in scope | test_metacommunity |
| `meta_offset_environmental_values` | unchanged |  | in scope | test_metacommunity |
| `meta_dead_selection` | docstring / comment only | NemoBackend | in scope | common, test_metacommunity |
| `meta_mainland_asexual_birth_mutate_germinate` | unchanged |  | asex/mix (probe) | test_metacommunity |
| `meta_mainland_sexual_birth_mutate_germinate` | docstring / comment only | NemoBackend | in scope | common, test_metacommunity |
| `meta_mainland_mixed_birth_mutate_germinate` | unchanged |  | asex/mix (probe) | test_metacommunity |
| `meta_colonize_from_propagules_rains` | changed |  | in scope | common, test_metacommunity |
| `pairwise_sexual_colonization_from_prpagules_rains` | changed |  | in scope | common, test_metacommunity |
| `meta_asex_reproduce_mutate_into_offspring_pool` | unchanged |  | asex/mix (probe) | test_metacommunity |
| `meta_sex_reproduce_mutate_into_offspring_pool` | docstring / comment only | NemoBackend | in scope | common, test_metacommunity |
| `meta_mix_reproduce_mutate_into_offspring_pool` | unchanged |  | asex/mix (probe) | test_metacommunity |
| `meta_asex_reproduce_calculation_into_offspring_marker_pool` | unchanged |  | asex/mix (probe) | test_metacommunity |
| `meta_sex_reproduce_calculation_with_offspring_marker_pool` | unchanged |  | in scope | common, test_metacommunity |
| `meta_mix_reproduce_calculation_with_offspring_marker_pool` | unchanged |  | asex/mix (probe) | test_metacommunity |
| `matrix_around` | unchanged |  | in scope | test_metacommunity |
| `dispersal_kernel_strength` | unchanged |  | in scope | test_metacommunity |
| `calculate_dispersal_kernel_strength_matrix` | unchanged |  | in scope | test_metacommunity |
| `normalized_calculate_dispersal_kernel_strength_matrix` | unchanged |  | in scope | test_metacommunity |
| `get_disp_among_rate_matrix` | unchanged |  | in scope | test_metacommunity |
| `get_offspring_marker_pool_num_matrix` | unchanged |  | in scope | test_metacommunity |
| `get_offspring_pool_num_matrix` | unchanged |  | in scope | test_metacommunity |
| `get_dormance_pool_num_matrix` | unchanged |  | in scope | test_metacommunity |
| `get_patch_empty_sites_num_matrix` | unchanged |  | in scope | test_metacommunity |
| `get_emigrants_matrix` | unchanged |  | in scope | test_metacommunity |
| `get_immigrants_matrix` | unchanged |  | in scope | test_metacommunity |
| `get_dispersal_among_num_matrix` | unchanged |  | in scope | test_metacommunity |
| `dispersal_among_patches_from_offspring_pool_and_dormancy_pool` | changed | NemoBackend | in scope | common, test_metacommunity |
| `dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool` | changed |  | in scope | common, test_metacommunity |
| `dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool` | changed |  | in scope | common, test_metacommunity |
| `dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool` | changed |  | in scope | common, test_metacommunity |
| `meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool` | changed |  | in scope | common, test_metacommunity |
| `meta_dispersal_within_patch_from_offspring_to_immigrant_pool` | changed |  | in scope | common, test_metacommunity |
| `meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool` | changed |  | in scope | common, test_metacommunity |
| `meta_dispersal_within_patch_from_offspring_and_dormancy_pool` | changed | NemoBackend | in scope | common, test_metacommunity |
| `meta_local_germinate_from_offspring_and_dormancy_pool` | changed |  | in scope | common, test_metacommunity |
| `meta_local_germinate_from_offspring_and_immigrant_pool` | changed |  | in scope | common, test_metacommunity |
| `meta_local_germinate_from_offspring_immigrant_and_dormancy_pool` | changed |  | in scope | common, test_metacommunity |
| `meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool` | changed | NemoBackend | in scope | common, test_metacommunity |
| `meta_dormancy_process_from_offspring_pool_to_dormancy_pool` | docstring / comment only | NemoBackend | in scope | common, test_metacommunity |
| `meta_dormancy_process_from_offspring_pool_and_immigrant_pool` | docstring / comment only | NemoBackend | in scope | common, test_metacommunity |
| `meta_clear_up_microsite_individuals` | new | NemoBackend | in scope | common, test_metacommunity |
| `meta_clear_up_offspring_pool` | new | NemoBackend | in scope | common, test_metacommunity |
| `meta_clear_up_immigrant_pool` | new | NemoBackend | in scope | common, test_metacommunity |
| `meta_clear_up_dormancy_pool` | new | NemoBackend | in scope | common, test_metacommunity |
| `meta_clear_up_offspring_marker_pool` | new |  | in scope | test_metacommunity |
| `meta_clear_up_immigrant_marker_pool` | new |  | in scope | test_metacommunity |
| `meta_clear_up_offspring_and_immigrant_pool` | changed | NemoBackend | in scope | common, test_metacommunity |
| `meta_clear_up_offspring_marker_and_immigrant_marker_pool` | changed |  | in scope | common, test_metacommunity |
| `meta_disturbance_process_in_patches` | docstring / comment only | NemoBackend | in scope | common, test_metacommunity |
| `meta_disturbance_process_in_habitat` | docstring / comment only | NemoBackend | in scope | common, test_metacommunity |
| `get_meta_microsites_optimum_sp_id_val` | unchanged |  | in scope | test_metacommunity |
| `get_meta_microsite_environment_values` | unchanged |  | in scope | test_metacommunity |
| `get_meta_microsites_individuals_sp_id_values` | unchanged |  | in scope | test_metacommunity |
| `get_meta_microsites_individuals_phenotype_values` | unchanged |  | in scope | test_metacommunity |
| `columns_patch_habitat_microsites_id` | unchanged |  | in scope | test_metacommunity |
| `meta_distribution_data_all_time_to_csv_gz` | unchanged |  | in scope | test_metacommunity |

Installed by `extension/global_habitat_network.py` (not part of the class body, not migrated): the marker network dispersal runs in test_chain (marker schedule, odd steps) and the fuzz; the offspring network dispersal is probed only (two owners, KNOWN).

## Repository state

```
 M metaibm/metacommunity.py
?? test/metaibm-nemo-bench/
?? test/metaibm-nemo-run/
?? test/nemoback-metacommunity/
 metaibm/metacommunity.py | 181 +++++++++++++++++++++++++++++------------------
 1 file changed, 112 insertions(+), 69 deletions(-)
```
