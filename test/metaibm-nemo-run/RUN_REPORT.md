# MetaIBM-Nemo full simulations (sexual, experiments/*.py interfaces)

Generated 2026-10-18 22:44 by `test/metaibm-nemo-run/run_all.py`.

- python `3.11.17`, numpy `2.4.6`, `macOS-15.7.3-arm64-arm-64bit`; HEAD `7a1adbb` (`MIGRATE patch: NemoBackend (sexual) + test/nemoback-patch`), branch `features/developing`

| file | md5 | state |
|---|---|---|
| `metaibm/metacommunity.py` | `b262a4394d5de41a9a902543f3b1a6e3` | modified (working tree) |
| `metaibm/patch.py` | `29e70117ed847b7b66eacadc77bf84cf` | committed |
| `metaibm/habitat.py` | `cc2c2001c092d149325b6d1f42e4d2a7` | committed |
| `metaibm/indfactory.py` | `eacaf91df4dbf3ab7dc72100323355be` | committed |
| `metaibm/individual.py` | `2914e40044bb48e0b10364227282d915` | committed |
| `extension/global_habitat_network.py` | `e9a76892a27f0e2ef95e10d9ecd6afd4` | committed |

Nothing in metaibm/, extension/ or experiments/ was edited: the run scripts copy each experiment's landscape generator and time loop, with only these interface changes: sexual functions instead of the (not migrated) mix / asex ones, `meta_initialize(indfactories, species_2_phenotype_ls, reproduce_mode)` with one Nemo indfactory per species, a distinct `global_index` per habitat, CSV rows every `--save-every` steps (the experiments write every step: ~1 GB per run at 102400 sites), invariant checks every 100 steps, teardown at the end. The global-habitat-network dispersal is not tested for now (GRFE uses the experiment's commented patch-level alternative).

## Runs

| run | command | status | steps | s / step (after t=100) | alive at last check | species on islands | RSS MB first -> last simulation check (max) | teardown minted == recycled (individuals) |
|---|---|---|---|---|---|---|---|---|
| model | run_model.py  | done | 2000 / 2000 | 5.42 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 102400} | 2 | 399.2 -> 1307.3 (max 1373.3); after plots 2026.7 | OK 29667226 |
| sloss-p16 | run_sloss.py --patch-num 16 | done | 2000 / 2000 | 5.60 | {'mainland': 6394, 'metacommunity': 102400} | 13 | 375.1 -> 1460.4 (max 1465.6); after plots 2365.7 | OK 26763403 |
| sloss-p64-hetero | run_sloss.py --patch-num 64 --hetero | done | 2000 / 2000 | 3.87 | {'mainland': 6392, 'metacommunity': 102400} | 11 | 718.0 -> 1370.9 (max 1435.8); after plots 2164.8 | OK 26772117 |
| sloss-p1 | run_sloss.py --patch-num 1 | stopped by user (SIGTERM, rc=-15) | 100 / 2000 |  | {'mainland': 6392, 'metacommunity': 45259} | 12 | 338.9 -> 594.4 (max 594.4) | - |
| sloss-pool-p16 | run_sloss.py --patch-num 16 --schedule pool --dormancy-total 10240 | done | 2000 / 2000 | 2.42 | {'mainland': 6394, 'metacommunity': 102400} | 11 | 612.0 -> 1553.5 (max 1803.0); after plots 2600.3 | OK 85196430 |
| grfe-p16 | run_grfe.py --patch-num 16 | done | 2000 / 2000 | 2.82 | {'mainland': 6394, 'metacommunity': 102400} | 11 | 517.2 -> 1218.9 (max 1427.0); after plots 1980.1 | OK 25785388 |

## High-disturbance short run (`run_disturbance.py`: sloss-p16 landscape, patch_dist_rate 0.01, 300 steps)

- 40 disturbances (expected 48), every one checked: the disturbed patches' individuals (microsites + offspring / immigrant pools) all recycled with the exact per-factory count, their marker pools emptied, dormancy kept, every other patch untouched, full invariants after the call
- recycled by disturbance 177662 individuals (mean 4442 per event), 40681 markers cleared; teardown minted == recycled: True; 2042.5 s
- recolonization of a disturbed patch (6400 sites), mean individuals after 1 / 5 / 10 / 25 / 50 / 100 steps: 159 / 772 / 1941 / 4304 / 4792 / 5642

## Details

### model

- interface: experiments/model.py: 1 trait, 2 mainlands (100x100), 100 islands x 1 habitat 32x32, two colonizations / step, environment +0.1 every 100 steps
- invariant checks (one owner, alive, minted - recycled == held, microsite bookkeeping): 22, all passed
- files: {'rep=0-meta_species_distribution_all_time.csv.gz': 0.1, 'rep=0-meta_phenotype_all_time.csv.gz': 27.7}; total 10488.4 s
- phase seconds: meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool 5093.6, (plots) 1844.4, meta_dead_selection 1586.3, dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool 504.9, mainland1.meta_dead_selection 440.0, mainland2.meta_dead_selection 424.4, (data: phenotype csv.gz) 267.2, (data: sp_id csv.gz) 125.5, meta_colonize_from_propagules_rains 65.9, meta_sex_reproduce_calculation_with_offspring_marker_pool 56.7
- environment changes: 19, final island mean [1.9000000000000006]

| t | elapsed s | RSS MB | alive | minted | recycled |
|---|---|---|---|---|---|
| 1 | 2.5 | 399.2 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 20} | 22168 | 2168 |
| 200 | 678.2 | 960.3 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 102400} | 2154158 | 2031778 |
| 400 | 1976.1 | 1285.5 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 102400} | 5353253 | 5230873 |
| 600 | 3139.4 | 1298.2 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 101378} | 8415365 | 8294007 |
| 800 | 4059.3 | 1371.2 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 102400} | 11432320 | 11309940 |
| 1000 | 5055.2 | 1045.9 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 102400} | 14481655 | 14359275 |
| 1200 | 5913.1 | 1286.3 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 102400} | 17516924 | 17394544 |
| 1400 | 6713.4 | 1295.7 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 102400} | 20558342 | 20435962 |
| 1600 | 7446.8 | 1302.7 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 102400} | 23611042 | 23488662 |
| 1800 | 8041.5 | 1307.1 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 102400} | 26643668 | 26521288 |
| 2000 | 8643.2 | 1307.3 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 102400} | 29667226 | 29544846 |
| 2000 | 10488.1 | 2026.7 | {'mainland1': 9990, 'mainland2': 9990, 'metacommunity': 102400} | 29667226 | 29544846 |

### sloss-p16

- interface: experiments/model-sloss.py: 16 species, mainland 16 x 20x20, 16 patches x 4 habitats 40x40 (discrete environment), pairwise colonization
- invariant checks (one owner, alive, minted - recycled == held, microsite bookkeeping): 22, all passed
- files: {'rep=0-meta_species_distribution_all_time.csv.gz': 1.7, 'rep=0-meta_x_axis_phenotype_all_time.csv.gz': 26.4, 'rep=0-meta_y_axis_phenotype_all_time.csv.gz': 27.1}; total 10889.4 s
- phase seconds: meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool 7859.6, meta_dead_selection 1817.1, (data: phenotype csv.gz) 508.5, (plots) 226.7, (data: sp_id csv.gz) 136.1, pairwise_sexual_colonization_from_prpagules_rains 97.9, mainland.meta_dead_selection 89.5, meta_sex_reproduce_calculation_with_offspring_marker_pool 54.9, dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool 46.5, meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool 18.4

| t | elapsed s | RSS MB | alive | minted | recycled |
|---|---|---|---|---|---|
| 1 | 2.9 | 375.1 | {'mainland': 6394, 'metacommunity': 6} | 7168 | 768 |
| 200 | 1026.6 | 903.7 | {'mainland': 6394, 'metacommunity': 102400} | 2259528 | 2150734 |
| 400 | 2737.8 | 1341.3 | {'mainland': 6394, 'metacommunity': 102400} | 5120432 | 5011638 |
| 600 | 4058.7 | 1457.6 | {'mainland': 6394, 'metacommunity': 102400} | 7883578 | 7774784 |
| 800 | 5382.6 | 1387.1 | {'mainland': 6392, 'metacommunity': 102400} | 10619114 | 10510322 |
| 1000 | 6480.1 | 1444.3 | {'mainland': 6392, 'metacommunity': 102400} | 13339876 | 13231084 |
| 1200 | 7477.0 | 1451.0 | {'mainland': 6392, 'metacommunity': 102400} | 16042117 | 15933325 |
| 1400 | 8277.3 | 1454.9 | {'mainland': 6394, 'metacommunity': 102400} | 18733143 | 18624349 |
| 1600 | 9081.3 | 1457.3 | {'mainland': 6394, 'metacommunity': 102400} | 21417928 | 21309134 |
| 1800 | 9881.9 | 1458.4 | {'mainland': 6394, 'metacommunity': 102400} | 24093968 | 23985174 |
| 2000 | 10662.1 | 1460.4 | {'mainland': 6394, 'metacommunity': 102400} | 26763403 | 26654609 |
| 2000 | 10889.1 | 2365.7 | {'mainland': 6394, 'metacommunity': 102400} | 26763403 | 26654609 |

### sloss-p64-hetero

- interface: experiments/model-sloss.py: 64 patches x 4 habitats 20x20, continuous environment (is_heterogeneity=True)
- invariant checks (one owner, alive, minted - recycled == held, microsite bookkeeping): 22, all passed
- files: {'rep=0-meta_species_distribution_all_time.csv.gz': 1.9, 'rep=0-meta_x_axis_phenotype_all_time.csv.gz': 26.4, 'rep=0-meta_y_axis_phenotype_all_time.csv.gz': 28.4}; total 7425.9 s
- phase seconds: meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool 2708.8, meta_dead_selection 1825.9, (plots) 1481.5, (data: phenotype csv.gz) 604.0, dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool 270.2, (data: sp_id csv.gz) 153.8, pairwise_sexual_colonization_from_prpagules_rains 124.0, mainland.meta_dead_selection 112.5, meta_sex_reproduce_calculation_with_offspring_marker_pool 69.3, meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool 33.0

| t | elapsed s | RSS MB | alive | minted | recycled |
|---|---|---|---|---|---|
| 1 | 3.6 | 718.0 | {'mainland': 6394, 'metacommunity': 6} | 7168 | 768 |
| 200 | 343.6 | 915.8 | {'mainland': 6394, 'metacommunity': 102400} | 1806542 | 1697748 |
| 400 | 998.7 | 956.7 | {'mainland': 6394, 'metacommunity': 102400} | 4681443 | 4572649 |
| 600 | 1729.8 | 1373.5 | {'mainland': 6394, 'metacommunity': 102400} | 7482174 | 7373380 |
| 800 | 2474.8 | 1413.9 | {'mainland': 6394, 'metacommunity': 102400} | 10259339 | 10150545 |
| 1000 | 3152.5 | 1332.4 | {'mainland': 6394, 'metacommunity': 102400} | 13025080 | 12916286 |
| 1200 | 3708.1 | 1378.5 | {'mainland': 6394, 'metacommunity': 102400} | 15780676 | 15671882 |
| 1400 | 4268.4 | 1423.3 | {'mainland': 6394, 'metacommunity': 102400} | 18535137 | 18426343 |
| 1600 | 4849.9 | 1387.5 | {'mainland': 6394, 'metacommunity': 102400} | 21283710 | 21174916 |
| 1800 | 5454.4 | 1361.5 | {'mainland': 6394, 'metacommunity': 102400} | 24029300 | 23920506 |
| 2000 | 5943.6 | 1370.9 | {'mainland': 6392, 'metacommunity': 102400} | 26772117 | 26663325 |
| 2000 | 7425.6 | 2164.8 | {'mainland': 6392, 'metacommunity': 102400} | 26772117 | 26663325 |

### sloss-p1

- interface: experiments/model-sloss.py: 1 patch x 16 habitats 80x80 (among-patch dispersal takes the patch_num < 2 branch)
- invariant checks (one owner, alive, minted - recycled == held, microsite bookkeeping): 2, all passed
- files: None; total None s
- phase seconds: meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool 274.0, meta_dead_selection 27.0, pairwise_sexual_colonization_from_prpagules_rains 22.0, (data: sp_id csv.gz) 21.6, (data: phenotype csv.gz) 8.9, mainland.meta_dead_selection 4.7, mainland.meta_mainland_sexual_birth_mutate_germinate 1.0, meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool 0.4, (invariant checks) 0.4, meta_sex_reproduce_calculation_with_offspring_marker_pool 0.3

| t | elapsed s | RSS MB | alive | minted | recycled |
|---|---|---|---|---|---|
| 1 | 5.0 | 338.9 | {'mainland': 6394, 'metacommunity': 6} | 7168 | 768 |
| 100 | 365.1 | 594.4 | {'mainland': 6392, 'metacommunity': 45259} | 335072 | 283421 |

### sloss-pool-p16

- interface: experiments/model-sloss.py with its commented individual-pool lines (sexual): offspring pool -> pairwise colonization -> within / among dispersal to immigrant pools -> germinate O+I+D -> dormancy (10240 total, 160 / habitat) -> disturbance -> clear-ups (queued by queue_runs.py pool)
- invariant checks (one owner, alive, minted - recycled == held, microsite bookkeeping): 22, all passed
- files: {'rep=0-meta_species_distribution_all_time.csv.gz': 1.4, 'rep=0-meta_x_axis_phenotype_all_time.csv.gz': 26.4, 'rep=0-meta_y_axis_phenotype_all_time.csv.gz': 26.9}; total 4720.7 s
- phase seconds: meta_dead_selection 2100.4, meta_sex_reproduce_mutate_into_offspring_pool 831.9, (data: phenotype csv.gz) 595.2, (plots) 337.3, pairwise_sexual_colonization_from_prpagules_rains 164.8, meta_local_germinate_from_offspring_immigrant_and_dormancy_pool 164.0, (data: sp_id csv.gz) 162.8, mainland.meta_dead_selection 100.8, meta_dormancy_process_from_offspring_pool_and_immigrant_pool 86.6, dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool 85.2

| t | elapsed s | RSS MB | alive | minted | recycled |
|---|---|---|---|---|---|
| 1 | 4.2 | 612.0 | {'mainland': 6394, 'metacommunity': 6} | 7168 | 768 |
| 200 | 411.8 | 1703.2 | {'mainland': 6394, 'metacommunity': 102400} | 5332979 | 5213945 |
| 400 | 918.6 | 1772.6 | {'mainland': 6394, 'metacommunity': 102400} | 14125897 | 14006863 |
| 600 | 1367.4 | 1789.8 | {'mainland': 6394, 'metacommunity': 102400} | 22970349 | 22851315 |
| 800 | 1815.6 | 1792.9 | {'mainland': 6394, 'metacommunity': 102400} | 31838806 | 31719772 |
| 1000 | 2262.8 | 1803.0 | {'mainland': 6394, 'metacommunity': 102400} | 40718512 | 40599478 |
| 1200 | 2772.9 | 1469.3 | {'mainland': 6394, 'metacommunity': 102400} | 49600811 | 49481777 |
| 1400 | 3212.5 | 1545.4 | {'mainland': 6394, 'metacommunity': 102400} | 58494269 | 58375235 |
| 1600 | 3602.5 | 1550.5 | {'mainland': 6394, 'metacommunity': 102400} | 67391871 | 67272837 |
| 1800 | 3993.5 | 1552.5 | {'mainland': 6394, 'metacommunity': 102400} | 76286925 | 76167891 |
| 2000 | 4382.5 | 1553.5 | {'mainland': 6394, 'metacommunity': 102400} | 85196430 | 85077396 |
| 2000 | 4720.2 | 2600.3 | {'mainland': 6394, 'metacommunity': 102400} | 85196430 | 85077396 |

### grfe-p16

- interface: experiments/model-sloss-GRFE.py: landscape from the CSV layouts, 16 patches x 64 habitats 10x10, among-patch marker dispersal (rho=0.2)
- invariant checks (one owner, alive, minted - recycled == held, microsite bookkeeping): 22, all passed
- files: {'rep=0-meta_species_distribution_all_time.csv.gz': 1.0, 'rep=0-meta_x_axis_phenotype_all_time.csv.gz': 26.8, 'rep=0-meta_y_axis_phenotype_all_time.csv.gz': 26.7}; total 5439.7 s
- phase seconds: meta_dead_selection 1731.7, meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool 1098.5, meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool 955.8, (data: phenotype csv.gz) 617.5, (plots) 412.2, (data: sp_id csv.gz) 170.8, pairwise_sexual_colonization_from_prpagules_rains 139.1, mainland.meta_dead_selection 114.0, meta_sex_reproduce_calculation_with_offspring_marker_pool 86.3, dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool 66.1

| t | elapsed s | RSS MB | alive | minted | recycled |
|---|---|---|---|---|---|
| 1 | 6.7 | 517.2 | {'mainland': 6394, 'metacommunity': 6} | 7168 | 768 |
| 200 | 252.9 | 794.5 | {'mainland': 6394, 'metacommunity': 93622} | 1034886 | 934870 |
| 400 | 744.2 | 945.1 | {'mainland': 6394, 'metacommunity': 102400} | 3930435 | 3821641 |
| 600 | 1346.8 | 1257.3 | {'mainland': 6394, 'metacommunity': 102400} | 6691072 | 6582278 |
| 800 | 1966.1 | 1327.1 | {'mainland': 6394, 'metacommunity': 102400} | 9431346 | 9322552 |
| 1000 | 2540.9 | 1342.6 | {'mainland': 6392, 'metacommunity': 102400} | 12169401 | 12060609 |
| 1200 | 3108.3 | 1343.0 | {'mainland': 6394, 'metacommunity': 102400} | 14896067 | 14787273 |
| 1400 | 3567.8 | 1366.6 | {'mainland': 6394, 'metacommunity': 102400} | 17621638 | 17512844 |
| 1600 | 4033.2 | 1398.0 | {'mainland': 6394, 'metacommunity': 102400} | 20342242 | 20233448 |
| 1800 | 4498.3 | 1426.5 | {'mainland': 6394, 'metacommunity': 102400} | 23065338 | 22956544 |
| 2000 | 5026.6 | 1218.9 | {'mainland': 6394, 'metacommunity': 102400} | 25785388 | 25676594 |
| 2000 | 5439.3 | 1980.1 | {'mainland': 6394, 'metacommunity': 102400} | 25785388 | 25676594 |

## Repository state

```
 M metaibm/metacommunity.py
?? test/metaibm-nemo-run/
?? test/nemoback-metacommunity/
```
