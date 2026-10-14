# habitat.py test report (NemoBackend, sexual)

Generated 2026-10-14 15:44 by `test/nemoback-habitat/run_all.py`. Scope: `metaibm/habitat.py` (patch / metacommunity layers not tested; asex / mix only probed).

- python `3.11.17`, pynemo `2.4.1`, `macOS-15.7.3-arm64-arm-64bit`
- habitat.py md5 `cc2c2001c092d149325b6d1f42e4d2a7`; original reference = `git show baace34:metaibm/habitat.py` (HEAD `c3711cd`)

## Smoke: PASS

| step | return |
|---|---|
| hab_initialize | `0` |
| hab_dead_selection | `5` |
| hab_sex_reproduce_mutate_into_offspring_pool | `13` |
| hab_local_germinate_from_offspring_and_immigrant_pool | `5` |
| hab_clear_up_offspring_pool + immigrant_pool | `8` |
| reproduce | `16` |
| hab_local_germinate_from_offspring_immigrant_and_dormancy_pool | `4` |
| hab_dormancy_process_from_offspring_pool_to_dormancy_pool | `(0, 0, 10, 10)` |
| reproduce | `17` |
| hab_dormancy_process_from_offspring_pool_and_immigrant_pool | `(0, 10, 10, 10)` |
| hab_sex_reproduce_calculation_into_offspring_marker_pool | `17` |
| hab_sexual_reprodece_germinate | `3` |
| habitat_disturbance_process | `53` |
| teardown (recycle_everything) | `10` |

## Tests: 2 / 2 PASS

| script | covers | result | time |
|---|---|---|---|
| test_bookkeeping.py | pool / microsite bookkeeping vs original, dormancy branches, clear-up, sampling uniformity (fake individuals, no Nemo) | PASS | 6.0 s |
| test_nemo.py | real Nemo individuals: reproduction | lifecycle | memory | edges_known (one section each) | PASS | 13.7 s |

## Key numbers

<details><summary>test_bookkeeping.py</summary>

```
1 germination == original: 550 cases; is_remove=False touches nothing
2 dormancy: table 11/11, sweep 2x2500, == original 2x2500 (eliminate diff {'A2': 1652, 'B1': 1580, 'B3': 114, 'B2': 650}), runtime K change consistent
3 clear-up x6, disturbance == original x10 (only marker pools differ), two-owner duplicates raise
4 dead selection == original: 480 cases; recycled == dead
5 uniformity chi2: {'germination_candidates(df49)': 33.8, 'dormancy_A2(df24)': 13.2, 'sample_without_replacement(df9)': 5.9}
PASSED
```
</details>
<details><summary>test_nemo.py</summary>

```
2 death probability vs survival_rate: z = [-0.5, 0.13, -1.29]
3 parents + per-locus inheritance OK, |pheno - sum(alleles)| = 0e+00, sex ratio z=1.67, parent chi2 54.8/70.7, mutation 24/2480
4 count rule (fractional birth {0.37: [0.76, 0.745], 0.81: [0.88, 0.88]}), markers mint nothing, direct-birth germination places every newborn
5 two species: same-species parents, offspring minted by the mother's factory, both factories balanced
6 single sex / one pair / rate 0 / leftover offspring recycled / rate 1.5 -> ValueError
---- section reproduction PASSED
==== section lifecycle
1 v3 loops: 18 configurations, invariants after every phase, recovery from dormancy in every dormancy run
2 reproducibility in two fresh processes: {'individual L=12 seed=7': 'identical over 50 steps', 'dormB L=12 seed=7': 'identical over 50 steps'}
3 selection response: {'selection_0.35': [0.198, 0.299], 'neutral': [0.198, 0.205], 'offset_+0.15_at_gen10': [0.198, 0.299]}
4 three habitats, two species: {'moved': 237, 'habitats_with_two_species': []}
5 fuzz: 4500 random operations, invariants after each
---- section lifecycle PASSED
==== section memory
reference individuals unchanged after 4445 births / 4405 recycles; pointers + IDs unique; recycled individuals inert
long run L=30 x 400 steps, 136319 births: RSS {50: 132.7, 200: 132.7, 399: 132.7} MB
---- section memory PASSED
==== section edges_known
1 shapes 1x1 / 1x9 / 9x1 / 3x7 / 7x3 and an empty habitat: every function, invariants hold, empty -> all zero
2 constructor defaults / positional K, global_index bounds {0: 'ok, getHome()=0', 65535: 'ok, getHome()=65535', 65536: 'TypeError', -1: 'TypeError'}, bad species inputs {'mean not in species table': 'ValueError', 'missing indfactory key': 'KeyError', 'tuple mean_env_ls': 'no error'}
3 __str__, env helpers, survival_rate formula, pairwise empty / occupied helpers OK with Nemo individuals
4 asex / mix probes: ['hab_asex_reproduce_mutate_with_num: TypeError', 'hab_asex_reproduce_mutate_into_offspring_pool: TypeError', 'hab_asexual_reprodece_germinate: TypeError', 'hab_mix_asex_reproduce_mutate_with_num: TypeError', 'hab_mix_sex_reproduce_mutate_with_num: TypeError', 'hab_mix_reproduce_mutate_into_offspring_pool: TypeError', 'hab_mixed_reproduce_germinate: TypeError', 'hab_asex_reproduce_calculation_into_offspring_marker_pool: returned 13', 'hab_mix_reproduce_calculation_into_offspring_marker_pool: returned 11']
5 caller-contract violations recorded
---- section edges_known PASSED
PASSED
```
</details>

## Findings / known behaviour (12)

| severity | finding |
|---|---|
| info | sexual_birth_rate > 1: int(pairs*rate) > pairs -> ValueError from random.sample (one offspring per pair per call, as in the original); the old offspring_pool is already recycled, invariants hold. |
| info | global_index is not validated by habitat.__init__; out of the unsigned-short range Nemo refuses it at the first make_founder: {65536: 'TypeError', -1: 'TypeError'} (nothing minted). |
| info | hab_initialize finds the species by exact equality species_2_phenotype_ls.index(mean_env_ls) and names it "sp<k+1>": {'mean not in species table': 'ValueError', 'missing indfactory key': 'KeyError', 'tuple mean_env_ls': 'no error'}; no founder is minted on these errors. A tuple mean_env_ls now works: habitat.__init__ copies it into a list. hab_reset/offset_environment_values change mean_env_ls in place, so hab_initialize after an environment change also raises (as in the original). |
| info | unknown dead-selection method string ("gaussian") -> no error from survival_rate (no validation; unchanged from the original). |
| info | hab_sexual_pairwise_parents_ls() / _num() shuffle the position lists in self.species_category in place, so even a count consumes Python RNG (same as the original; bookkeeping stays consistent). |
| info | asex / mix functions are not migrated: hab_asex_reproduce_mutate_with_num -> TypeError: cannot pickle 'swig_runtime_data5.SwigPyObject' object; hab_asex_reproduce_mutate_into_offspring_pool -> TypeError: cannot pickle 'swig_runtime_data5.SwigPyObject' object; hab_asexual_reprodece_germinate -> TypeError: cannot pickle 'swig_runtime_data5.SwigPyObject' object; hab_mix_asex_reproduce_mutate_with_num -> TypeError: cannot pickle 'swig_runtime_data5.SwigPyObject' object; hab_mix_sex_reproduce_mutate_with_num -> TypeError: cannot pickle 'swig_runtime_data5.SwigPyObject' object; hab_mix_reproduce_mutate_into_offspring_pool -> TypeError: cannot pickle 'swig_runtime_data5.SwigPyObject' object; hab_mixed_reproduce_germinate -> TypeError: cannot pickle 'swig_runtime_data5.SwigPyObject' object; hab_asex_reproduce_calculation_into_offspring_marker_pool -> returned 13; hab_mix_reproduce_calculation_into_offspring_marker_pool -> returned 11 |
| info | germinating an already-recycled individual -> RuntimeError, invariants broken: add_individual writes grid / empty / occupied / indi_num before it reads gender. |
| info | add_individual on an occupied microsite only prints and returns; the individual stays with the caller (hab_initialize resets the microsites first, so it never hits this). |
| info | _recycle_dropped detects a two-owner individual only when it is dropped; if kept it may end up twice in dormancy_pool (copies=2 here). Callers must never create two owners. |
| info | clear-up functions empty a pool before its recycle loop: if a recycle raises (two owners), the rest of that pool is neither recycled nor held. The error stops the run anyway. |
| info | sample_offspring_without_replacement() hands ownership to the caller (patch dispersal): the returned individuals must go into a pool or be recycled -- to verify when patch.py is migrated. |
| info | habitat_disturbance_process() returns recycled individuals + cleared markers; the original returned None and kept the marker pools. |

## Method coverage

| method | vs original | scope | scripts | reached via |
|---|---|---|---|---|
| `__init__` | changed | in scope | 1 | constructor |
| `__str__` | unchanged | in scope | 0 | str(hab) |
| `hab_reset_environment_values` | unchanged | in scope | 1 |  |
| `hab_offset_environment_values` | unchanged | in scope | 1 |  |
| `_empty_microsite_grid` | unchanged | in scope | 0 | __init__ |
| `add_individual` | unchanged | in scope | 2 |  |
| `del_individual` | changed | in scope | 1 |  |
| `get_hab_pairwise_empty_site_pos_ls` | unchanged | in scope | 1 |  |
| `get_hab_pairwise_occupied_site_pos_ls` | unchanged | in scope | 1 |  |
| `sample_offspring_without_replacement` | unchanged | in scope | 2 |  |
| `hab_initialize` | changed | in scope | 2 |  |
| `get_microsite_env_val_ls` | unchanged | in scope | 1 | dead selection |
| `survival_rate` | unchanged | in scope | 1 | dead selection |
| `hab_dead_selection` | changed | in scope | 3 |  |
| `hab_asex_reproduce_mutate_with_num` | unchanged | asex/mix (probe only) | 1 |  |
| `hab_sexual_pairwise_parents_ls` | unchanged | in scope | 0 | ..._num() |
| `hab_sexual_pairwise_parents_num` | unchanged | in scope | 1 |  |
| `hab_sex_reproduce_mutate_with_num` | changed | in scope | 1 |  |
| `hab_mixed_sexual_pairwise_parents_ls` | unchanged | asex/mix (probe only) | 0 |  |
| `hab_mixed_sexual_pairwse_parents_num` | unchanged | asex/mix (probe only) | 0 |  |
| `hab_mixed_asexual_parent_num` | unchanged | asex/mix (probe only) | 0 |  |
| `hab_mix_asex_reproduce_mutate_with_num` | unchanged | asex/mix (probe only) | 1 |  |
| `hab_mix_sex_reproduce_mutate_with_num` | unchanged | asex/mix (probe only) | 1 |  |
| `hab_asexual_reprodece_germinate` | unchanged | asex/mix (probe only) | 1 |  |
| `hab_sexual_reprodece_germinate` | changed | in scope | 2 |  |
| `hab_mixed_reproduce_germinate` | unchanged | asex/mix (probe only) | 1 |  |
| `hab_asex_reproduce_calculation_into_offspring_marker_pool` | unchanged | asex/mix (probe only) | 1 |  |
| `hab_sex_reproduce_calculation_into_offspring_marker_pool` | unchanged | in scope | 2 |  |
| `hab_mix_reproduce_calculation_into_offspring_marker_pool` | unchanged | asex/mix (probe only) | 1 |  |
| `hab_asex_reproduce_mutate_into_offspring_pool` | unchanged | asex/mix (probe only) | 1 |  |
| `hab_sex_reproduce_mutate_into_offspring_pool` | changed | in scope | 2 |  |
| `hab_mix_reproduce_mutate_into_offspring_pool` | unchanged | asex/mix (probe only) | 1 |  |
| `_germinate_from_pools` | new | in scope | 0 | the 3 germination wrappers |
| `hab_local_germinate_from_offspring_and_dormancy_pool` | changed | in scope | 1 |  |
| `hab_local_germinate_from_offspring_and_immigrant_pool` | changed | in scope | 2 |  |
| `hab_local_germinate_from_offspring_immigrant_and_dormancy_pool` | changed | in scope | 3 |  |
| `_recycle_dropped` | new | in scope | 0 | dormancy processes |
| `hab_dormancy_process_from_offspring_pool_to_dormancy_pool` | changed | in scope | 3 |  |
| `hab_dormancy_process_from_offspring_pool_and_immigrant_pool` | changed | in scope | 2 |  |
| `hab_clear_up_microsite_individuals` | new | in scope | 3 |  |
| `hab_clear_up_offspring_pool` | new | in scope | 4 |  |
| `hab_clear_up_immigrant_pool` | new | in scope | 4 |  |
| `hab_clear_up_dormancy_pool` | new | in scope | 2 |  |
| `hab_clear_up_offspring_marker_pool` | new | in scope | 2 |  |
| `hab_clear_up_immigrant_marker_pool` | new | in scope | 2 |  |
| `habitat_disturbance_process` | changed | in scope | 3 |  |

## Repository state

```
 M metaibm/habitat.py
 M metaibm/patch.py
 M test/nemoback-habitat/common.py
 M test/nemoback-habitat/run_all.py
 M test/nemoback-habitat/test_nemo.py
?? test/nemoback-patch/
 metaibm/habitat.py |   4 +-
 metaibm/patch.py   | 132 +++++++++++++++++++++++++++++++++++++++++------------
 2 files changed, 105 insertions(+), 31 deletions(-)
```
