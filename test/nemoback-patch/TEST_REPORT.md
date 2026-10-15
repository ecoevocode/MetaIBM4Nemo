# patch.py test report (NemoBackend, sexual)

Generated 2026-10-15 14:27 by `test/nemoback-patch/run_all.py`. Scope: `metaibm/patch.py` and its call chain patch -> habitat -> indfactory -> individual (metacommunity layer only scanned statically; asex / mix only probed).

## Environment

- python `3.11.17`, pynemo `2.4.1`, numpy `2.4.6`, `macOS-15.7.3-arm64-arm-64bit`
- HEAD `ff1a7bb` (`FIX habitat: own copies of mean_env_ls / var_env_ls`), branch `features/developing`

| file | md5 | state |
|---|---|---|
| `metaibm/patch.py` | `29e70117ed847b7b66eacadc77bf84cf` | modified (working tree) |
| `metaibm/habitat.py` | `cc2c2001c092d149325b6d1f42e4d2a7` | committed |
| `metaibm/indfactory.py` | `eacaf91df4dbf3ab7dc72100323355be` | committed |
| `metaibm/individual.py` | `2914e40044bb48e0b10364227282d915` | committed |

- original reference: `git show baace34:metaibm/patch.py` (blob `0e90181`; verified to still contain `hab_clear_up_offspring_and_immigrant_pool`). `baace34` did not touch patch.py itself: the last commit that changed it is `ef779f8`, so `baace34:patch.py` == `ef779f8:patch.py` == `HEAD:patch.py`. Full original stack = patch.py + habitat.py from `baace34`.

## Step 0 -- baseline (existing suites, run unchanged)

| suite | script | result | time |
|---|---|---|---|
| test/nemo-indfactory | accept_genetics.py | PASS | 0.6 s |
| test/nemo-indfactory | accept_guards.py | PASS | 0.6 s |
| test/nemo-indfactory | accept_lifecycle.py | PASS | 0.6 s |
| test/nemo-indfactory | smoke.py | PASS | 0.6 s |
| test/nemoback-habitat (run_all.py) | smoke.py | PASS | 0.7 s |
| test/nemoback-habitat (run_all.py) | test_bookkeeping.py | PASS | 6.2 s |
| test/nemoback-habitat (run_all.py) | test_nemo.py | PASS | 14.2 s |

- habitat suite folder backed up before and restored after its run (its run_all.py rewrites TEST_REPORT.md + results/): restored clean = `True`; its fresh report is kept as `results/baseline_habitat_TEST_REPORT.md`.
- habitat suite original reference: `ORIGINAL_COMMIT = 'baace34'`; `metaibm/habitat.py` is committed (last commit `ff1a7bb FIX habitat: own copies of mean_env_ls / var_env_ls`).

## Smoke: PASS

| step | return | invariants |
|---|---|---|
| patch_initialize | `None` | True |
| patch_dead_selection | `9` | True |
| patch_sex_reproduce_mutate_into_offspring_pool | `28` | True |
| patch_dispersal_within_from_offspring_pool_to_immigrant_pool | `11` | True |
| patch_local_germinate_from_offspring_and_immigrant_pool | `9` | True |
| patch_dormancy_process_from_offspring_pool_and_immigrant_pool | `(0, 0, 19, 19)` | True |
| patch_dead_selection | `21` | True |
| patch_sex_reproduce_mutate_into_offspring_pool | `21` | True |
| patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool | `20` | True |
| patch_local_germinate_from_offspring_immigrant_and_dormancy_pool | `21` | True |
| patch_clear_up_offspring_pool | `8` | True |
| patch_clear_up_immigrant_pool | `7` | True |
| patch_dead_selection | `28` | True |
| patch_sex_reproduce_mutate_into_offspring_pool | `19` | True |
| patch_dipersal_within_from_offspring_and_dormancy_pool | `9` | True |
| patch_local_germinate_from_offspring_and_dormancy_pool | `13` | True |
| patch_dormancy_process_from_offspring_pool_to_dormancy_pool | `(0, 0, 1, 1)` | True |
| patch_sex_reproduce_calculation_into_offspring_marker_pool | `26` | True |
| patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool | `13` | True |
| patch_clear_up_offspring_marker_pool | `13` | True |
| patch_clear_up_immigrant_marker_pool | `13` | True |
| patch_dead_selection | `22` | True |
| patch_sexual_birth_germinate | `19` | True |
| patch_disturbance_process | `66` | True |
| patch_initialize (re-initialize, dormancy kept) | `None` | True |
| patch_clear_up_dormancy_pool | `1` | True |
| patch_clear_up_microsite_individuals | `75` | True |
| teardown (recycle_everything) | `0` | True |

## Tests: 3 / 3 PASS

| script | covers | result | time |
|---|---|---|---|
| test_bookkeeping.py | six clear-ups, disturbance, is_remove=False guards, equivalence with the original (bodies on the current habitat + full original stack), forced / natural extras, sampling uniformity, edge shapes (fake individuals, no Nemo) | PASS | 4.2 s |
| test_nemo.py | real Nemo individuals, one fresh process per section: initialize / reproduction / two_species / lifecycle / getters / memory / reproducibility / fuzz / edges_known | PASS | 13.9 s |
| test_chain.py | call chain patch -> habitat -> indfactory -> individual (method wrapping), ordering rule, offspring homepatch, static scan of callers | PASS | 0.7 s |

## Key numbers

- equivalence with `baace34` (fake individuals): {'O->I': 400, 'OD->I': 400, 'OD->sites': 400, 'marker True': 400, 'marker False': 400, 'germ': 450, 'dorm': 300, 'dead': 300, 'marker_birth': 150, 'unchanged_helpers': 150, 'dorm_K0_dormants_dropped_vs_original': 746}
- OD->sites extras recycled == dropped by the original: 138 over the sweep; forced extras {'cases': 200, 'extras_recycled': 6601}
- natural extras (P = 0.25): {'extras1': 149, 'arrive0': 150, 'arrive1': 301, 'arrive2': 149, 'seeds': 600, 'p_extras': 0.248, 'expected': 0.25}
- uniformity chi2 (loose bound df + 6 sqrt(2 df)): {'OD->sites migrants (df14)': 13.1, 'OD->sites destination sites (df5)': 2.2, 'O->I migrants (df11)': 10.5}
- memory: RSS {'50': 118.9, '160': 119.1, '319': 119.2} MB, address reuse {'late_births': 18898, 'late_births_on_recycled_address': 0.9981, 'distinct_addresses_vs_alive_high_water_mark': {'50': {'sp1': [457, 457], 'sp2': [213, 213]}, '100': {'sp1': [475, 475], 'sp2': [213, 213]}, '160': {'sp1': [495, 495], 'sp2': [213, 213]}, '240': {'sp1': [529, 529], 'sp2': [213, 213]}, '319': {'sp1': [531, 531], 'sp2': [213, 213]}}, 'growth_after_t160': 0.0508}, teardown (minted, recycled) {'sp1': [34267, 34267], 'sp2': [4378, 4378]}
- fuzz: 3300 operations
- reproducibility: identical nemo_digest over 50 steps in two fresh processes
- chain: 640 recycle() calls checked for the ordering rule, 377 make_offspring calls checked for homepatch

<details><summary>test_bookkeeping.py output</summary>

```
1 six clear-ups: 180 cases; return == pool size, only that pool emptied, recycled == former contents (in order), markers recycle nothing
2 disturbance: 40 cases; return == sum of habitat counts, dormancy untouched, marker pools emptied, == original apart from markers (original returns None)
3 is_remove=False: 120 cases; the three individual-pool dispersals raise NotImplementedError before any RNG use / state change
4 equivalence with baace34 (original bodies on the current habitat + full original stack): {'O->I': 400, 'OD->I': 400, 'OD->sites': 400, 'marker True': 400, 'marker False': 400, 'germ': 450, 'dorm': 300, 'dead': 300, 'marker_birth': 150, 'unchanged_helpers': 150, 'dorm_K0_dormants_dropped_vs_original': 746}; OD->sites recycled == dropped by the original (138 extras over all seeds)
5 forced extras: 200 cases, 6601 extras recycled; return == empties filled (== all empty sites when migrants >= empties), invariants hold
6 natural extras over 600 seeds: P(extra) = 0.248 (expected 0.25, 4-sigma), arrivals 0/1/2 = 150/301/149; invariants every time
7 uniformity chi2: {'OD->sites migrants (df14)': 13.1, 'OD->sites destination sites (df5)': 2.2, 'O->I migrants (df11)': 10.5}
8 1 habitat 3x3                errors none | warnings ['O->I: divide by zero encountered in divide', 'O->I: invalid value encountered in divide', 'OD->I: divide by zero encountered in divide', 'OD->I: invalid value encountered in divide', 'OD->sites: divide by zero encountered in divide', 'OD->sites: divide by zero encountered in divide; invalid value encountered in divide', 'OD->sites: invalid value encountered in divide', 'marker: divide by zero encountered in divide', 'marker: invalid value encountered in divide', 'rate_matrix: divide by zero encountered in divide']
8 2 habitats 3x3               errors none | warnings ['OD->sites: invalid value encountered in divide']
8 3 habitats 2x5 (non-square)  errors none | warnings ['OD->sites: invalid value encountered in divide']
8 4 habitats 1x1               errors none | warnings ['OD->sites: invalid value encountered in divide']
8 2 habitats full              errors none | warnings ['OD->sites: invalid value encountered in divide']
8 3 habitats empty             errors none | warnings ['OD->sites: invalid value encountered in divide']
8 empty patch (0 habitats)     errors {'migrants_to_patch': 'ZeroDivisionError: integer division or modulo by zero', 'migrants_marker_to_patch': 'ZeroDivisionError: integer division or modulo by zero', 'sp_id_values': "UnboundLocalError: cannot access local variable 'h_object' where it is not asso", 'env_values': "UnboundLocalError: cannot access local variable 'h_object' where it is not asso", 'optimum_sp_id': "UnboundLocalError: cannot access local variable 'h_object' where it is not asso"} | warnings none
8 unequal sizes 2x2 + 3x3      errors {'sp_id_values': 'ValueError: cannot reshape array of size 13 into shape (2,9)', 'env_values': 'ValueError: cannot reshape array of size 13 into shape (2,9)', 'optimum_sp_id': 'ValueError: cannot reshape array of size 13 into shape (2,9)'} | warnings ['OD->sites: invalid value encountered in divide']
PASSED
```
</details>
<details><summary>test_nemo.py output</summary>

```
==== section initialize
1 patch_initialize: every microsite filled, species per habitat from the table, getHome == global_index, original gender draw order, no np.random; re-initialize recycled 62 residents, kept 60 pool members; asexual -> all female; bad mean (first habitat) -> ValueError, nothing minted
---- section initialize PASSED
==== section reproduction
2 reproduction via patch: 514 offspring-pool births + 323 direct births; parents same species + same habitat, mother's factory, getHome == mating habitat (20/21/22), per-locus inheritance OK, |pheno - sum(alleles)| = 0e+00, sex ratio z = -1.42; fractional rate per habitat lo..lo+1
---- section reproduction PASSED
==== section two_species
3 two species in different habitats + dispersal: habitats that held both species ['h0', 'h1', 'h2'], 2947 offspring checked, 0 cross-species / cross-habitat pairs, both factories balanced
---- section two_species PASSED
==== section lifecycle
4 lifecycle: 12 configurations x 40 steps, invariants after every phase, factory balance every step; OD->sites extras recycled {'K=0 OD->sites disturbance=False': 2, 'K=0 OD->sites disturbance=True': 1, 'K=30 OD->sites disturbance=False': 1, 'K=30 OD->sites disturbance=True': 0}
---- section lifecycle PASSED
==== section getters
5 getters with real individuals: sp_id, phenotype (phenotype_set == get_indi_phenotype_ls), environment, optimum sp_id, get_data identity, matrices, pools, pairwise lists OK; offset -> re-initialize ValueError; shared var list -> own copies, var [0.035, 0.035, 0.035]
---- section getters PASSED
==== section memory
6 memory: 320 steps, 38645 births; RSS {50: 118.9, 160: 119.1, 319: 119.2} MB; late births on a recycled address 99.8%, per factory (distinct addresses, alive high-water mark) {50: {'sp1': (457, 457), 'sp2': (213, 213)}, 100: {'sp1': (475, 475), 'sp2': (213, 213)}, 160: {'sp1': (495, 495), 'sp2': (213, 213)}, 240: {'sp1': (529, 529), 'sp2': (213, 213)}, 319: {'sp1': (531, 531), 'sp2': (213, 213)}} (growth after t=160: 5.1%); reference individuals unchanged; recycled individuals inert; teardown minted == recycled {'sp1': (34267, 34267), 'sp2': (4378, 4378)}
---- section memory PASSED
==== section reproducibility
7 reproducibility: identical nemo_digest over 50 steps in two fresh processes ; seed 8 differs
---- section reproducibility PASSED
==== section fuzz
8 fuzz: 3300 random patch-level operations (25 kinds) on two patches, invariants after each, minted == recycled after every run
---- section fuzz PASSED
==== section edges_known
9.1 mean not in table (h2)               {'error': 'ValueError', 'minted': 18, 'recycled': 18, 'per_habitat (indi_num, old residents kept)': {'h0': (9, 0), 'h1': (9, 0), 'h2': (0, 0)}}
9.1 missing indfactory key (h1 -> sp2)   {'error': 'KeyError', 'minted': 9, 'recycled': 9, 'per_habitat (indi_num, old residents kept)': {'h0': (9, 0), 'h1': (9, 9), 'h2': (0, 0)}}
9.1 global_index 65536 (h1)              {'error': 'TypeError', 'minted': 9, 'recycled': 18, 'per_habitat (indi_num, old residents kept)': {'h0': (9, 0), 'h1': (0, 0), 'h2': (0, 0)}}
9.2 global_index {0: 'ok, getHome()=0', 65535: 'ok, getHome()=65535', 65536: 'TypeError, minted 0', -1: 'TypeError, minted 0'} | inputs {'tuple mean_env_ls in table': 'ValueError', "reproduce_mode='Sexual'": 'UnboundLocalError', 'sexual_birth_rate 1.5': 'ValueError'}
9.4 asex / mix probes: ['asexual_birth_germinate: TypeError', 'mixed_birth_germinate: TypeError', 'asex_reproduce_calculation_into_offspring_marker_pool: returned 25', 'mix_reproduce_calculation_into_offspring_marker_pool: returned 24', 'asex_reproduce_mutate_into_offspring_pool: TypeError', 'mix_reproduce_mutate_into_offspring_pool: TypeError']
---- section edges_known PASSED
PASSED
```
</details>
<details><summary>test_chain.py output</summary>

```
1 call chains: 51 public patch methods called (asex/mix probed in test_nemo.py); 28 expected edge sets observed; make_founder only from patch_initialize, make_offspring only from the two sexual reproduction methods, indfactory.recycle only via individual.recycle; pool-to-pool dispersal + marker functions never reach recycle / indfactory
2 ordering rule: 640 recycle() calls, each on an individual already out of every microsite / pool, none touched afterwards
3 make_offspring: 377 calls, home_habitat == global_index of the habitat owning both parents, mother's factory, same species (roots: ['patch_sex_reproduce_mutate_into_offspring_pool', 'patch_sexual_birth_germinate'])
4 metaibm/metacommunity.py:193       meta_initialize                                            patch_initialize                                                         SIGNATURE CHANGED: call with 6 positional args does not bind to (self, indfactories, species_2_phenotype_ls, reproduce_mode='sexual') (too many positional arguments)
4 metaibm/metacommunity.py:311       meta_mainland_asexual_birth_mutate_germinate               patch_asexual_birth_germinate                                            asex / mix: not migrated (TypeError on Nemo individuals)
4 metaibm/metacommunity.py:333       meta_mainland_mixed_birth_mutate_germinate                 patch_mixed_birth_germinate                                              asex / mix: not migrated (TypeError on Nemo individuals)
4 metaibm/metacommunity.py:413       meta_asex_reproduce_mutate_into_offspring_pool             patch_asex_reproduce_mutate_into_offspring_pool                          asex / mix: not migrated (TypeError on Nemo individuals)
4 metaibm/metacommunity.py:431       meta_mix_reproduce_mutate_into_offspring_pool              patch_mix_reproduce_mutate_into_offspring_pool                           asex / mix: not migrated (TypeError on Nemo individuals)
4 metaibm/metacommunity.py:441       meta_asex_reproduce_calculation_into_offspring_marker_pool patch_asex_reproduce_calculation_into_offspring_marker_pool              asex / mix: not migrated (TypeError on Nemo individuals)
4 metaibm/metacommunity.py:461       meta_mix_reproduce_calculation_with_offspring_marker_pool  patch_mix_reproduce_calculation_into_offspring_marker_pool               asex / mix: not migrated (TypeError on Nemo individuals)
4 metaibm/metacommunity.py:730       meta_dispersal_within_patch_from_offspring_to_immigrant_pool patch_dispersal_within_from_offspring_pool_to_immigrant_pool             caller default is_remove=False -> NotImplementedError with the caller's defaults
4 metaibm/metacommunity.py:739       meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool caller default is_remove=False -> NotImplementedError with the caller's defaults
4 metaibm/metacommunity.py:748       meta_dispersal_within_patch_from_offspring_and_dormancy_pool patch_dipersal_within_from_offspring_and_dormancy_pool                   caller default is_remove=False -> NotImplementedError with the caller's defaults
4 metaibm/metacommunity.py:760       meta_local_germinate_from_offspring_and_dormancy_pool      patch_local_germinate_from_offspring_and_dormancy_pool                   caller default is_remove=False -> NotImplementedError with the caller's defaults
4 metaibm/metacommunity.py:771       meta_local_germinate_from_offspring_and_immigrant_pool     patch_local_germinate_from_offspring_and_immigrant_pool                  caller default is_remove=False -> NotImplementedError with the caller's defaults
4 metaibm/metacommunity.py:782       meta_local_germinate_from_offspring_immigrant_and_dormancy_pool patch_local_germinate_from_offspring_immigrant_and_dormancy_pool         caller default is_remove=False -> NotImplementedError with the caller's defaults
4 metaibm/metacommunity.py:864       meta_clear_up_offspring_and_immigrant_pool                 patch_clear_up_offspring_and_immigrant_pool                              REMOVED from patch (split into the six patch_clear_up_* functions)
4 metaibm/metacommunity.py:868       meta_clear_up_offspring_marker_and_immigrant_marker_pool   patch_clear_up_offspring_marker_and_immigrant_marker_pool                REMOVED from patch (split into the six patch_clear_up_* functions)
4 metaibm/simulator.py:58            build_empty_mainland_from_species_csv                      add_habitat                                                              no global_index passed -> every habitat gets homepatch 0
4 metaibm/simulator.py:95            build_empty_metacommunity_from_patch_habitat_csv           add_habitat                                                              no global_index passed -> every habitat gets homepatch 0
PASSED
```
</details>

## Performance (record only, `--perf`): seconds per call, mean of 3

| patch function | 4 habitats x 10x10 | 4 habitats x 20x20 | 4 habitats x 40x40 |
|---|---|---|---|
| `patch_initialize` | 0.0016 | 0.0069 | 0.0430 |
| `patch_dead_selection` | 0.0033 | 0.0148 | 0.0856 |
| `patch_sex_reproduce_mutate_into_offspring_pool` | 0.0010 | 0.0040 | 0.0177 |
| `patch_dispersal_within_from_offspring_pool_to_immigrant_pool` | 0.0002 | 0.0003 | 0.0020 |
| `patch_local_germinate_from_offspring_and_immigrant_pool` | 0.0001 | 0.0006 | 0.0036 |
| `patch_dormancy_process_from_offspring_pool_and_immigrant_pool` | 0.0000 | 0.0000 | 0.0001 |
| `patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool` | 0.0001 | 0.0003 | 0.0023 |
| `patch_dipersal_within_from_offspring_and_dormancy_pool` | 0.0002 | 0.0004 | 0.0016 |
| `patch_local_germinate_from_offspring_immigrant_and_dormancy_pool` | 0.0001 | 0.0004 | 0.0027 |
| `patch_dormancy_process_from_offspring_pool_to_dormancy_pool` | 0.0000 | 0.0001 | 0.0005 |
| `patch_sexual_birth_germinate` | 0.0001 | 0.0004 | 0.0017 |
| `patch_sex_reproduce_calculation_into_offspring_marker_pool` | 0.0001 | 0.0002 | 0.0010 |
| `patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool` | 0.0001 | 0.0002 | 0.0007 |
| `get_patch_microsites_individals_phenotype_values` | 0.0005 | 0.0022 | 0.0091 |
| `get_patch_microsites_individals_sp_id_values` | 0.0005 | 0.0022 | 0.0122 |
| `patch_disturbance_process` | 0.0003 | 0.0011 | 0.0044 |
| `patch_clear_up_dormancy_pool` | 0.0000 | 0.0002 | 0.0007 |

## Findings (18: bug 0, risk 1, info 5, known 12)

Severity: **bug** = wrong behaviour in the code under test; **risk** = correct today, fragile / surprising; **info** = observation; **known** = expected / next migration round / pre-existing.

| severity | finding | minimal reproducer | suggested fix (not applied) |
|---|---|---|---|
| risk | patch_initialize is not atomic: it calls hab_initialize habitat by habitat, so a bad habitat late in the order raises after earlier habitats were already re-initialized (residents recycled, new founders minted); with an out-of-range global_index the failing habitat itself is left EMPTY (hab_initialize recycles residents before the first make_founder raises). Observed: {'mean not in table (h2)': ('ValueError', 'minted 18', 'recycled 18'), 'missing indfactory key (h1 -> sp2)': ('KeyError', 'minted 9', 'recycled 9'), 'global_index 65536 (h1)': ('TypeError', 'minted 9', 'recycled 18')}. Invariants and factory balance still hold. | `p with h0 (0.2,0.2), h1 (0.5,0.5), h2 (0.9,0.9) all occupied; p.patch_initialize(F, [[0.2,0.2],[0.5,0.5]]) -> ValueError after h0/h1 were refilled (18 minted, 18 recycled)` | validate first: for every habitat resolve species_2_phenotype_ls.index(mean_env_ls), check the indfactory key and 0 <= global_index <= 65535 (and reproduce_mode), then initialize; or validate global_index in habitat.__init__ / add_habitat |
| info | patch_dispersal_within_*: when a column of the emigrant matrix sums to 0 (no offspring / dormants anywhere else) get_habs_immigrants_matrix divides 0/0 -> numpy RuntimeWarning, NaN entries are rounded to 0 by patch_matrix_around (pre-existing, harmless; seen in 7 of 8 edge shapes). |  |  |
| info | np.matrix is used by every dispersal matrix helper: numpy emits PendingDeprecationWarning "the matrix subclass is not the recommended way" (560 times in the edge section; pre-existing, will break when numpy removes np.matrix). |  |  |
| info | habitat layer (not patch): hab_dormancy_process_from_offspring_pool_to_dormancy_pool with dormancy_pool_max_size == 0 now empties (recycles) a non-empty dormancy_pool; the original kept it and reported eliminate=0 (746 dormants over the equivalence sweep). Only reachable when K is set to 0 with dormants present; intended by the habitat migration (comment: K changed mid-run). |  |  |
| info | patch_initialize finds species by exact equality species_2_phenotype_ls.index(mean_env_ls): a tuple row vs list mean -> ValueError; reproduce_mode="Sexual" -> UnboundLocalError (gender unbound in hab_initialize, after residents were recycled); no validation (habitat-level, as in the habitat report). |  |  |
| info | patch_initialize now returns None (the original returned 0); metacommunity.meta_initialize ignores the return value. |  |  |
| known | get_patch_microsites_individals_sp_id_values / _phenotype_values / get_patch_microsites_environment_values / get_patch_microsites_optimum_sp_id_value_array reshape with the LAST habitat's size: unequal habitat sizes -> ValueError; an empty patch -> UnboundLocalError (pre-existing, unchanged from the original). |  |  |
| known | get_patch_dispersal_within_rate_matrix divides by hab_num-1: with 1 habitat numpy only warns (divide by zero encountered in divide) and the 1x1 matrix is still right because the diagonal is overwritten; every dispersal returns 0 (pre-existing). |  |  |
| known | split() / migrants_to_patch_into_habs_immigrant_pool index pieces by h_object.index: hab_index 1,2 (not 0..n-1) -> IndexError; an empty patch -> ZeroDivisionError (pre-existing). |  |  |
| known | Marker dispersal with is_remove=False keeps the marker at the source and sends a copy (model semantics; can inflate births). Both modes are identical to the original. |  |  |
| known | add_habitat(global_index=...) is not validated; out of the unsigned-short range Nemo refuses it at the first make_founder: {0: 'ok, getHome()=0', 65535: 'ok, getHome()=65535', 65536: 'TypeError, minted 0', -1: 'TypeError, minted 0'}. |  |  |
| known | sexual_birth_rate > 1 -> ValueError from random.sample at habitat level (one offspring per pair per call); invariants hold afterwards. |  |  |
| known | asex / mix patch functions are not migrated: patch_asexual_birth_germinate -> TypeError; patch_mixed_birth_germinate -> TypeError; patch_asex_reproduce_calculation_into_offspring_marker_pool -> returned 25; patch_mix_reproduce_calculation_into_offspring_marker_pool -> returned 24; patch_asex_reproduce_mutate_into_offspring_pool -> TypeError; patch_mix_reproduce_mutate_into_offspring_pool -> TypeError (TypeError = cannot pickle SwigPyObject in copy.deepcopy; the two marker calculations create no individual). |  |  |
| known | metaibm/metacommunity.py (next migration round) still calls the removed patch_clear_up_offspring_and_immigrant_pool / patch_clear_up_offspring_marker_and_immigrant_marker_pool and patch_initialize with the old 6-argument signature: metaibm/metacommunity.py:193 patch_initialize; metaibm/metacommunity.py:864 patch_clear_up_offspring_and_immigrant_pool; metaibm/metacommunity.py:868 patch_clear_up_offspring_marker_and_immigrant_marker_pool. |  |  |
| known | metacommunity meta_dispersal_within_patch_* / meta_local_germinate_* default is_remove=False and pass it through, so with their defaults they now raise NotImplementedError (patch / habitat default is True): metaibm/metacommunity.py:730 meta_dispersal_within_patch_from_offspring_to_immigrant_pool, metaibm/metacommunity.py:739 meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool, metaibm/metacommunity.py:748 meta_dispersal_within_patch_from_offspring_and_dormancy_pool, metaibm/metacommunity.py:760 meta_local_germinate_from_offspring_and_dormancy_pool, metaibm/metacommunity.py:771 meta_local_germinate_from_offspring_and_immigrant_pool, metaibm/metacommunity.py:782 meta_local_germinate_from_offspring_immigrant_and_dormancy_pool. |  |  |
| known | metaibm/simulator.py builds habitats with patch.add_habitat(...) without global_index, so every habitat would get Nemo homepatch 0 (next round): metaibm/simulator.py:58, metaibm/simulator.py:95. |  |  |
| known | metacommunity still calls the unmigrated asex / mix patch functions: metaibm/metacommunity.py:311 patch_asexual_birth_germinate, metaibm/metacommunity.py:333 patch_mixed_birth_germinate, metaibm/metacommunity.py:413 patch_asex_reproduce_mutate_into_offspring_pool, metaibm/metacommunity.py:431 patch_mix_reproduce_mutate_into_offspring_pool, metaibm/metacommunity.py:441 patch_asex_reproduce_calculation_into_offspring_marker_pool, metaibm/metacommunity.py:461 patch_mix_reproduce_calculation_into_offspring_marker_pool. |  |  |
| known | metacommunity among-patch dispersal (next round) with is_remove=False samples individuals with replacement from get_patch_offspring_pool() / get_patch_offspring_and_dormancy_pool() and appends the SAME objects to another patch's immigrant pools -> two owners (breaks the one-owner rule once individuals are Nemo-backed): metaibm/metacommunity.py:612, metaibm/metacommunity.py:646, metaibm/metacommunity.py:707. |  |  |

## Method coverage

`vs original` compares each method with `baace34:metaibm/patch.py` (AST + source text). `scripts` = test scripts that name the method; `reached via` = patch methods that call it internally (from the chain test).

| method | vs original | scope | scripts | reached via |
|---|---|---|---|---|
| `__init__` | changed | in scope | test_chain, test_nemo | constructor |
| `get_data` | unchanged | in scope | test_chain, test_nemo |  |
| `add_habitat` | changed | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `patch_initialize` | changed | in scope | smoke, test_chain, test_nemo |  |
| `patch_dead_selection` | docstring / comment only | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_individual_num` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_empty_sites_ls` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo | patch_empty_sites_num |
| `patch_empty_sites_num` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_pairwise_empty_sites_ls` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_pairwise_occupied_sites_ls` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_occupied_sites_ls` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_offspring_marker_pool` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_offspring_pool` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_dormancy_pool` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_offspring_and_dormancy_pool` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `sample_offspring_without_replacement` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_microsites_individals_sp_id_values` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_microsites_individals_phenotype_values` | unchanged | in scope | test_chain, test_nemo |  |
| `get_patch_microsites_environment_values` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `patch_offset_environmental_values` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `patch_asexual_birth_germinate` | unchanged | asex/mix (probe only) | test_nemo |  |
| `patch_sexual_birth_germinate` | docstring / comment only | in scope | smoke, test_chain, test_nemo |  |
| `patch_mixed_birth_germinate` | unchanged | asex/mix (probe only) | test_nemo |  |
| `patch_asex_reproduce_calculation_into_offspring_marker_pool` | unchanged | asex/mix (probe only) | test_nemo |  |
| `patch_sex_reproduce_calculation_into_offspring_marker_pool` | unchanged | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_mix_reproduce_calculation_into_offspring_marker_pool` | unchanged | asex/mix (probe only) | test_nemo |  |
| `patch_asex_reproduce_mutate_into_offspring_pool` | unchanged | asex/mix (probe only) | test_nemo |  |
| `patch_sex_reproduce_mutate_into_offspring_pool` | docstring / comment only | in scope | smoke, test_chain, test_nemo |  |
| `patch_mix_reproduce_mutate_into_offspring_pool` | unchanged | asex/mix (probe only) | test_nemo |  |
| `split` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo | migrants_marker_to_patch_into_habs_immigrant_marker_pool, migrants_to_patch_into_habs_immigrant_pool |
| `migrants_to_patch_into_habs_immigrant_pool` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `migrants_marker_to_patch_into_habs_immigrant_marker_pool` | unchanged | in scope | test_bookkeeping, test_chain |  |
| `patch_matrix_around` | unchanged | in scope | test_bookkeeping, test_chain | 5 patch methods, e.g. get_dispersal_within_num_matrix, patch_dipersal_within_from_offspring_and_dormancy_pool |
| `get_patch_dispersal_within_rate_matrix` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo | 7 patch methods, e.g. get_dispersal_within_num_matrix, get_habs_emigrants_matrix |
| `get_patch_habs_offspring_num_matrix` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo | 6 patch methods, e.g. get_dispersal_within_num_matrix, get_habs_emigrants_matrix |
| `get_patch_habs_offspring_marker_num_matrix` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo | patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool |
| `get_patch_habs_dormancy_num_matrix` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo | 5 patch methods, e.g. get_dispersal_within_num_matrix, get_habs_emigrants_matrix |
| `get_patch_habs_empty_sites_num_matrix` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo | 3 patch methods, e.g. get_dispersal_within_num_matrix, get_habs_immigrants_matrix |
| `get_habs_emigrants_matrix` | unchanged | in scope | test_chain, test_nemo | 4 patch methods, e.g. get_dispersal_within_num_matrix, get_habs_immigrants_matrix |
| `get_habs_immigrants_matrix` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo | get_dispersal_within_num_matrix, patch_dipersal_within_from_offspring_and_dormancy_pool |
| `get_dispersal_within_num_matrix` | unchanged | in scope | test_chain, test_nemo | patch_dipersal_within_from_offspring_and_dormancy_pool |
| `patch_dipersal_within_from_offspring_and_dormancy_pool` | changed | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool` | changed | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_dispersal_within_from_offspring_pool_to_immigrant_pool` | changed | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool` | changed | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_local_germinate_from_offspring_and_dormancy_pool` | changed | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_local_germinate_from_offspring_and_immigrant_pool` | changed | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_local_germinate_from_offspring_immigrant_and_dormancy_pool` | changed | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_dormancy_process_from_offspring_pool_to_dormancy_pool` | docstring / comment only | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_dormancy_process_from_offspring_pool_and_immigrant_pool` | docstring / comment only | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_clear_up_microsite_individuals` | new | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_clear_up_offspring_pool` | new | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_clear_up_immigrant_pool` | new | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_clear_up_dormancy_pool` | new | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_clear_up_offspring_marker_pool` | new | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_clear_up_immigrant_marker_pool` | new | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `patch_disturbance_process` | changed | in scope | smoke, test_bookkeeping, test_chain, test_nemo |  |
| `get_patch_microsites_optimum_sp_id_value_array` | unchanged | in scope | test_bookkeeping, test_chain, test_nemo |  |
| `patch_clear_up_offspring_and_immigrant_pool` | removed | in scope | - | replaced by the six patch_clear_up_* (callers: metacommunity.py, next round) |
| `patch_clear_up_offspring_marker_and_immigrant_marker_pool` | removed | in scope | - | replaced by the six patch_clear_up_* (callers: metacommunity.py, next round) |

## Call-chain matrix (observed)

What each public patch method reached while it was the outermost call (union over all calls in `test_chain.py`).

| patch method | habitat | indfactory | individual | patch (internal) |
|---|---|---|---|---|
| `add_habitat` | __init__, _empty_microsite_grid | - | - | - |
| `patch_initialize` | _empty_microsite_grid, add_individual, hab_clear_up_microsite_individuals, hab_initialize | make_founder, recycle | gender, recycle | - |
| `patch_dead_selection` | del_individual, get_microsite_env_val_ls, hab_dead_selection, survival_rate | recycle | gender, get_indi_phenotype_ls, recycle | - |
| `patch_empty_sites_num` | - | - | - | get_patch_empty_sites_ls |
| `get_patch_pairwise_empty_sites_ls` | get_hab_pairwise_empty_site_pos_ls | - | - | - |
| `get_patch_pairwise_occupied_sites_ls` | get_hab_pairwise_occupied_site_pos_ls | - | - | - |
| `get_patch_microsites_individals_phenotype_values` | - | - | phenotype_set | - |
| `patch_offset_environmental_values` | hab_offset_environment_values, hab_reset_environment_values | - | - | - |
| `patch_sexual_birth_germinate` | add_individual, hab_sex_reproduce_mutate_with_num, hab_sexual_pairwise_parents_ls, hab_sexual_pairwise_parents_num, hab_sexual_reprodece_germinate | make_offspring | gender | - |
| `patch_sex_reproduce_calculation_into_offspring_marker_pool` | hab_sex_reproduce_calculation_into_offspring_marker_pool, hab_sexual_pairwise_parents_ls, hab_sexual_pairwise_parents_num | - | - | - |
| `patch_sex_reproduce_mutate_into_offspring_pool` | hab_clear_up_offspring_pool, hab_sex_reproduce_mutate_into_offspring_pool, hab_sex_reproduce_mutate_with_num, hab_sexual_pairwise_parents_ls, hab_sexual_pairwise_parents_num | make_offspring, recycle | recycle | - |
| `migrants_to_patch_into_habs_immigrant_pool` | - | - | - | split |
| `migrants_marker_to_patch_into_habs_immigrant_marker_pool` | - | - | - | split |
| `get_habs_emigrants_matrix` | - | - | - | get_patch_dispersal_within_rate_matrix, get_patch_habs_dormancy_num_matrix, get_patch_habs_offspring_num_matrix |
| `get_habs_immigrants_matrix` | - | - | - | get_habs_emigrants_matrix, get_patch_dispersal_within_rate_matrix, get_patch_habs_dormancy_num_matrix, get_patch_habs_empty_sites_num_matrix, get_patch_habs_offspring_num_matrix |
| `get_dispersal_within_num_matrix` | - | - | - | get_habs_emigrants_matrix, get_habs_immigrants_matrix, get_patch_dispersal_within_rate_matrix, get_patch_habs_dormancy_num_matrix, get_patch_habs_empty_sites_num_matrix, get_patch_habs_offspring_num_matrix, patch_matrix_around |
| `patch_dipersal_within_from_offspring_and_dormancy_pool` | add_individual, sample_offspring_without_replacement | recycle | gender, recycle | get_dispersal_within_num_matrix, get_habs_emigrants_matrix, get_habs_immigrants_matrix, get_patch_dispersal_within_rate_matrix, get_patch_habs_dormancy_num_matrix, get_patch_habs_empty_sites_num_matrix, get_patch_habs_offspring_num_matrix, patch_matrix_around |
| `patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool` | sample_offspring_without_replacement | - | - | get_patch_dispersal_within_rate_matrix, get_patch_habs_offspring_marker_num_matrix, patch_matrix_around |
| `patch_dispersal_within_from_offspring_pool_to_immigrant_pool` | sample_offspring_without_replacement | - | - | get_patch_dispersal_within_rate_matrix, get_patch_habs_offspring_num_matrix, patch_matrix_around |
| `patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool` | sample_offspring_without_replacement | - | - | get_habs_emigrants_matrix, get_patch_dispersal_within_rate_matrix, get_patch_habs_dormancy_num_matrix, get_patch_habs_offspring_num_matrix, patch_matrix_around |
| `patch_local_germinate_from_offspring_and_dormancy_pool` | _germinate_from_pools, add_individual, hab_local_germinate_from_offspring_and_dormancy_pool | - | gender | - |
| `patch_local_germinate_from_offspring_and_immigrant_pool` | _germinate_from_pools, add_individual, hab_local_germinate_from_offspring_and_immigrant_pool | - | gender | - |
| `patch_local_germinate_from_offspring_immigrant_and_dormancy_pool` | _germinate_from_pools, add_individual, hab_local_germinate_from_offspring_immigrant_and_dormancy_pool | - | gender | - |
| `patch_dormancy_process_from_offspring_pool_to_dormancy_pool` | _recycle_dropped, hab_dormancy_process_from_offspring_pool_to_dormancy_pool | recycle | recycle | - |
| `patch_dormancy_process_from_offspring_pool_and_immigrant_pool` | _recycle_dropped, hab_dormancy_process_from_offspring_pool_and_immigrant_pool | recycle | recycle | - |
| `patch_clear_up_microsite_individuals` | _empty_microsite_grid, hab_clear_up_microsite_individuals | recycle | recycle | - |
| `patch_clear_up_offspring_pool` | hab_clear_up_offspring_pool | recycle | recycle | - |
| `patch_clear_up_immigrant_pool` | hab_clear_up_immigrant_pool | recycle | recycle | - |
| `patch_clear_up_dormancy_pool` | hab_clear_up_dormancy_pool | recycle | recycle | - |
| `patch_clear_up_offspring_marker_pool` | hab_clear_up_offspring_marker_pool | - | - | - |
| `patch_clear_up_immigrant_marker_pool` | hab_clear_up_immigrant_marker_pool | - | - | - |
| `patch_disturbance_process` | _empty_microsite_grid, hab_clear_up_immigrant_marker_pool, hab_clear_up_immigrant_pool, hab_clear_up_microsite_individuals, hab_clear_up_offspring_marker_pool, hab_clear_up_offspring_pool, habitat_disturbance_process | recycle | recycle | - |
| `get_patch_microsites_optimum_sp_id_value_array` | survival_rate | - | - | - |

Pure patch-level (no habitat / indfactory / individual call): `get_data`, `get_patch_individual_num`, `get_patch_empty_sites_ls`, `get_patch_occupied_sites_ls`, `get_patch_offspring_marker_pool`, `get_patch_offspring_pool`, `get_patch_dormancy_pool`, `get_patch_offspring_and_dormancy_pool`, `sample_offspring_without_replacement`, `get_patch_microsites_individals_sp_id_values`, `get_patch_microsites_environment_values`, `split`, `patch_matrix_around`, `get_patch_dispersal_within_rate_matrix`, `get_patch_habs_offspring_num_matrix`, `get_patch_habs_offspring_marker_num_matrix`, `get_patch_habs_dormancy_num_matrix`, `get_patch_habs_empty_sites_num_matrix`

Not called in the chain process (asex / mix, probed in separate processes by `test_nemo.py` edges_known): `patch_asexual_birth_germinate`, `patch_mixed_birth_germinate`, `patch_asex_reproduce_calculation_into_offspring_marker_pool`, `patch_mix_reproduce_calculation_into_offspring_marker_pool`, `patch_asex_reproduce_mutate_into_offspring_pool`, `patch_mix_reproduce_mutate_into_offspring_pool`

### Static scan: patch calls in metacommunity.py / simulator.py / extension/*.py

| location | caller | patch method | issue |
|---|---|---|---|
| metaibm/metacommunity.py:193 | meta_initialize | `patch_initialize` | SIGNATURE CHANGED: call with 6 positional args does not bind to (self, indfactories, species_2_phenotype_ls, reproduce_mode='sexual') (too many positional arguments) |
| metaibm/metacommunity.py:311 | meta_mainland_asexual_birth_mutate_germinate | `patch_asexual_birth_germinate` | asex / mix: not migrated (TypeError on Nemo individuals) |
| metaibm/metacommunity.py:333 | meta_mainland_mixed_birth_mutate_germinate | `patch_mixed_birth_germinate` | asex / mix: not migrated (TypeError on Nemo individuals) |
| metaibm/metacommunity.py:413 | meta_asex_reproduce_mutate_into_offspring_pool | `patch_asex_reproduce_mutate_into_offspring_pool` | asex / mix: not migrated (TypeError on Nemo individuals) |
| metaibm/metacommunity.py:431 | meta_mix_reproduce_mutate_into_offspring_pool | `patch_mix_reproduce_mutate_into_offspring_pool` | asex / mix: not migrated (TypeError on Nemo individuals) |
| metaibm/metacommunity.py:441 | meta_asex_reproduce_calculation_into_offspring_marker_pool | `patch_asex_reproduce_calculation_into_offspring_marker_pool` | asex / mix: not migrated (TypeError on Nemo individuals) |
| metaibm/metacommunity.py:461 | meta_mix_reproduce_calculation_with_offspring_marker_pool | `patch_mix_reproduce_calculation_into_offspring_marker_pool` | asex / mix: not migrated (TypeError on Nemo individuals) |
| metaibm/metacommunity.py:730 | meta_dispersal_within_patch_from_offspring_to_immigrant_pool | `patch_dispersal_within_from_offspring_pool_to_immigrant_pool` | caller default is_remove=False -> NotImplementedError with the caller's defaults |
| metaibm/metacommunity.py:739 | meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool | `patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool` | caller default is_remove=False -> NotImplementedError with the caller's defaults |
| metaibm/metacommunity.py:748 | meta_dispersal_within_patch_from_offspring_and_dormancy_pool | `patch_dipersal_within_from_offspring_and_dormancy_pool` | caller default is_remove=False -> NotImplementedError with the caller's defaults |
| metaibm/metacommunity.py:760 | meta_local_germinate_from_offspring_and_dormancy_pool | `patch_local_germinate_from_offspring_and_dormancy_pool` | caller default is_remove=False -> NotImplementedError with the caller's defaults |
| metaibm/metacommunity.py:771 | meta_local_germinate_from_offspring_and_immigrant_pool | `patch_local_germinate_from_offspring_and_immigrant_pool` | caller default is_remove=False -> NotImplementedError with the caller's defaults |
| metaibm/metacommunity.py:782 | meta_local_germinate_from_offspring_immigrant_and_dormancy_pool | `patch_local_germinate_from_offspring_immigrant_and_dormancy_pool` | caller default is_remove=False -> NotImplementedError with the caller's defaults |
| metaibm/metacommunity.py:864 | meta_clear_up_offspring_and_immigrant_pool | `patch_clear_up_offspring_and_immigrant_pool` | REMOVED from patch (split into the six patch_clear_up_* functions) |
| metaibm/metacommunity.py:868 | meta_clear_up_offspring_marker_and_immigrant_marker_pool | `patch_clear_up_offspring_marker_and_immigrant_marker_pool` | REMOVED from patch (split into the six patch_clear_up_* functions) |
| metaibm/simulator.py:58 | build_empty_mainland_from_species_csv | `add_habitat` | no global_index passed -> every habitat gets homepatch 0 |
| metaibm/simulator.py:95 | build_empty_metacommunity_from_patch_habitat_csv | `add_habitat` | no global_index passed -> every habitat gets homepatch 0 |
| metaibm/metacommunity.py:612 |  | `random.sample(get_patch_..._pool())` | among-patch is_remove=False: same individual in two patches |
| metaibm/metacommunity.py:646 |  | `random.sample(get_patch_..._pool())` | among-patch is_remove=False: same individual in two patches |
| metaibm/metacommunity.py:707 |  | `random.sample(get_patch_..._pool())` | among-patch is_remove=False: same individual in two patches |

`extension/global_habitat_network.py` reads `patch_object.set` / habitat pools directly and calls no patch method: nothing to flag.

## Repository state

```
 M metaibm/patch.py
?? test/nemoback-patch/
 metaibm/patch.py | 136 ++++++++++++++++++++++++++++++++++++++++++-------------
 1 file changed, 105 insertions(+), 31 deletions(-)
```
