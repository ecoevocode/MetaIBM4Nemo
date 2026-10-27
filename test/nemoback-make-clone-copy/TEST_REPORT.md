# make_copy / make_clone migration -- test report

Generated 2026-10-27 09:20 by `test/nemoback-make-clone-copy/run_all.py`. Scope: the uncommitted migration on top of `140acc4` -- `individual.make_copy()` (Nemo `Individual::operator=`), `indfactory.make_clone()` (Nemo `IndFactory::makeClonedOffsprg`) and the core changes that use them (A: 10 `is_remove=False` copy sites, B: asexual / mixed births, C: marker births, D: NemoBackend tags). Test-only: nothing outside this folder was edited, nothing committed.

## Environment

- date `2026-10-27 09:20`, HEAD `140acc4` (`UPDATE matrix_around / patch_matrix_around: vectorized per-row systematic rounding`), branch `features/developing`
- python `3.11.17`, pynemo `2.4.1`, numpy `2.4.6`, machine `macOS-15.7.3-arm64-arm-64bit, Apple M2, 8 CPUs`
- conda env `metaibm-nemo-test2`; every script in its own process with `PYTHONDONTWRITEBYTECODE=1`

| changed file (working tree) | md5 |
|---|---|
| `extension/global_habitat_network.py` | `658b74471a45a9a6ee63257b0cf3d6b0` |
| `metaibm/habitat.py` | `2c026d6b5315b4a2dd42923e741e6bbf` |
| `metaibm/indfactory.py` | `e31efcb199af3adf6a8f3c470576f1c1` |
| `metaibm/individual.py` | `4b269122e36a3bf43e97a58982cc9846` |
| `metaibm/metacommunity.py` | `89fc3f547e063968ff5205d84338851c` |
| `metaibm/patch.py` | `7170e7ca8666b530ba20dc9da877fbab` |

- decisions under test: D1 `is_remove` defaults stay True; D2 births into `offspring_pool` round with `int(x + U)`; D3 `hab_mixed_reproduce_germinate` makes exactly min(empty, a + s) offspring (asex_num = min(a, np.around(E a / (a + s))), sex_num = min(s, E - asex_num)), all placed; D4 `mutation_rate` / `pheno_var_ls` kept in the interface, unused; D5 a birth habitat without parents raises `ValueError` (v3 behaviour, not fixed here)

## Results: 454 checks, 0 failed

| script / part | covers | checks | failed | findings | time | result |
|---|---|---|---|---|---|---|
| `test_core.py unit` | static checks S1-S25; make_copy (operator= fields + alleles, ID, Ve); make_clone (pedigree, mutation law); guards; lifecycle (recycling, RSS); rng (two processes, Python imitation); bench | 86 | 0 | 1 | 9.4 s | PASS |
|  | - static | 25 | 0 | 0 |  |  |
|  | - copy | 18 | 0 | 1 |  |  |
|  | - clone | 19 | 0 | 0 |  |  |
|  | - guards | 13 | 0 | 0 |  |  |
|  | - lifecycle | 6 | 0 | 0 |  |  |
|  | - rng | 5 | 0 | 0 |  |  |
| `test_core.py core` | A: the 10 is_remove=False sites (fake individuals vs v3 seed by seed, colonization False vs True, real Nemo individuals); B + C: births with_num, into sites / pools, mixed split, D wrappers, marker births, error-path probes | 354 | 0 | 4 | 6.2 s | PASS |
|  | - copy a1: fake individuals vs v3 | 128 | 0 | 0 |  |  |
|  | - copy a2: colonization False vs True | 32 | 0 | 0 |  |  |
|  | - copy b: real Nemo individuals | 140 | 0 | 0 |  |  |
|  | - births 1 with_num | 15 | 0 | 0 |  |  |
|  | - births 2 asexual birth into sites | 2 | 0 | 1 |  |  |
|  | - births 3 mixed birth into sites | 4 | 0 | 1 |  |  |
|  | - births 4 into offspring_pool + D wrappers | 10 | 0 | 0 |  |  |
|  | - births 5 marker births (C) | 20 | 0 | 0 |  |  |
|  | - births 6 probes | 3 | 0 | 2 |  |  |
| `test_simulation.py` | experiments/model.py, both branches, scaled: invariants, teardown, digests, NemoBackend vs pure v3 + genetics controls | 14 | 0 | 1 | 36.6 s | PASS |

`test_core.py unit` part bench: microseconds per call incl. recycle of the result (20000 reps, Ve > 0, mu = 1e-4)

| loci | make_copy | make_clone | make_offspring | make_founder | choice + recycle only |
|---|---|---|---|---|---|
| 20 | 1.95 | 1.77 | 2.07 | 1.91 | 0.14 |
| 80 | 1.99 | 1.99 | 2.45 | 2.73 | 0.16 |
| 320 | 2.31 | 2.49 | 4.02 | 5.75 | 0.15 |

## Findings (6 from the tests + 0 regression failures of category 3)

| severity | function(s) | finding | evidence | from |
|---|---|---|---|---|
| info | `individual.make_copy (all 10 is_remove=False sites)` | a copy keeps everything of its source but the ID: home (birth habitat), mother / father IDs, pedigree class, mating counters, age -- a dispersed or colonizing copy is NOT re-homed to its destination (as the v3 alias was the same object); Nemo computes the pedigree class of its offspring from these homes / parent IDs. MetaIBM never reads the copied _mother / _father raw pointers (IDs only). | {"source home / copy-of-copy home": [7, 7], "source mother_id / copy-of-copy mother_id": [1, 1], "counters copied": [3, 0, 0, 0, 0], "pedigree class of copy (home 7) x father (home 4) offspring": 0} | test_core.py unit / copy |
| info | `hab_mixed_reproduce_germinate` | split changed vs v3 (intended): when the empty sites cannot take both expectations, v3 rounded E*a/(a+s) and E*s/(a+s) each with np.around (round half to even), so an x.5 tie made E+1 offspring (one without a site) or only E-1 (one site left empty). Now asex_num = min(a, np.around(E*a/(a+s))) and sex_num = min(s, E - asex_num): exactly min(E, a + s) offspring, all placed, nothing recycled; on a tie the asexual part keeps the round-half-even result and the sexual part takes the rest. | {"v3 overshoot -> now (E, A, P, v3 (asex, sex), now (asex, sex)) e.g.": [[3, 2, 2, [2, 2], [2, 1]], [3, 3, 3, [2, 2], [2, 1]], [3, 4, 4, [2, 2], [2, 1]], [3, 5, 5, [2, 2], [2, 1]]], "n_overshoot": 8, "v3 undershoot -> now e.g.": [[1, 1, 1, [0, 0], [0, 1]], [1, 2, 2, [0, 0], [0, 1]], [1, 3, 3, [0, 0], [0, 1]], [1, 4, 4, [0, 0], [0, 1]]], "n_undershoot": 14} | test_core.py core / births 3 mixed birth into sites |
| info | `hab_mix_asex_reproduce_mutate_with_num, hab_mixed_reproduce_germinate` | asexual_parent_pos_ls is only rebuilt by hab_dead_selection: if a parent leaves its site in between (del_individual), the next mixed birth may pick the empty site -> AttributeError on None (v3: deepcopy(None) -> AttributeError too). Error path only: the clones already made in that call are lost (minted, never held or recycled), which matters only if a caller catches the error and continues. Not reachable in the experiments/ order (dead selection right before every birth). | {"hab_mix_asex_reproduce_mutate_with_num(all asexual parents)": {"error": "AttributeError: 'NoneType' object has no attribute 'indfactory'", "clones minted before the error": 8, "minted but neither held nor recycled": 8}, "hab_mixed_reproduce_germinate(1.0, 0.0)": {"error": "no error", "clones minted before the error": 11, "minted but neither held nor recycled": 0}} | test_core.py core / births 6 probes |
| info | `test_simulation asexual branch: occupancy (last 30 steps)` | NemoBackend (v3 semantics, 5 seeds) vs pure v3 (5 seeds) differ by 3.0 SE (statistical comparison, not a FAIL); controls vs v3: init_variance 0 -> 3.2 SE, init_variance 0 + v3-sized mutation effects -> 2.5 SE | {"nemo mean, SE": [116.02, 33.8253], "v3 mean, SE": [13.9533, 0.9788], "nemo init_variance 0 mean, SE": [58.5333, 13.829], "nemo init_variance 0 + mutation_variance 0.025^2 mean, SE": [8.6867, 1.877]} | test_simulation.py |
| known | `hab_asexual_reprodece_germinate, hab_asex_ / hab_mix_reproduce_mutate_into_offspring_pool, hab_mixed_reproduce_germinate` | asexual (or mixed-asexual) birth rate > 1: the number of births can exceed the number of (asexual) parents -> random.sample(..., num) raises ValueError (one offspring per parent per call, as for sexual_birth_rate > 1 / TODO(offspring_num)); the v3 bodies sample the same way. experiments/ use 0.5. | {"hab_asexual_reprodece_germinate(2.0)": {"error": "ValueError: Sample larger than population or is negative", "minted": 0, "recycled": 0, "residents unchanged": true, "offspring_pool": "3 -> 3", "check_held": true}, "hab_asex_reproduce_mutate_into_offspring_pool(2.0)": {"error": "ValueError: Sample larger than population or is negative", "minted": 0, "recycled": 3, "residents unchanged": true, "offspring_pool": "3 -> 0", "check_held": true}, "hab_mixed_reproduce_germinate(2.0, 0.0)": {"error": "ValueError: Sample larger than population or is negative", "minted": 0, "recycled": 0, "residents unchanged": true, "offspring_pool": "3 -> 3", "check_held": true}, "hab_mix_reproduce_mutate_into_off | test_core.py core / births 2 asexual birth into sites |
| known | `meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool (D5)` | a marker whose birth habitat has no (asexual / paired) parent any more raises ValueError from random.sample when the shuffled loop reaches it (markers handled before it are already born, see the per-seed counts); now also for the restored asexual / mix_asexual / mix_sexual markers. Not reachable in the experiments/model*.py order; v3 behaviour, not fixed in this round. | {"asexual": {"errors": {"ValueError: Sample larger than population or is negative": 10}, "born before the error (per seed)": [1, 1, 1, 3, 0, 2, 2, 2, 2, 1], "minted == born": true, "check_world": true}, "mix_asexual": {"errors": {"ValueError: Sample larger than population or is negative": 10}, "born before the error (per seed)": [2, 1, 1, 3, 0, 2, 2, 2, 2, 1], "minted == born": true, "check_world": true}, "mix_sexual": {"errors": {"ValueError: Sample larger than population or is negative": 10}, "born before the error (per seed)": [1, 1, 1, 3, 2, 2, 2, 0, 2, 3], "minted == born": true, "check_world": true}, "sexual": {"errors": {"ValueError: Sample larger than population or is negative": 10}, | test_core.py core / births 6 probes |


## Simulation (`test_simulation.py`)

- experiments/model.py loop, scaled: islands 4x4 patches x 1 habitat 8x8 (1024 sites), mainlands 20x20, 150 steps, environment +0.1 at t = 29 / 59 / 89 / 119, rain 0.10 per mainland per step (area-scaled), per-capita rates as in model.py; 52 runs, 36 s wall with 6 processes
- modes: `v3sem` = NemoBackend with `is_remove=False` passed to the three calls whose v3 default was False (2 colonizations, marker dispersal, marker germinate-and-birth); `default` = NemoBackend, defaults (True); `novar` = v3sem with `quanti_init_variance = 0`; `v3gen` = novar + `quanti_mutation_variance = 0.025^2` (size of a v3 bit flip); `v3pure` = pure v3 (git archive baace34 metaibm extension; metaibm/individual.py = baace34^ (67a742f); indfactory.py removed), v3 API, `is_remove=False`

| branch / kind / mode | runs | calls per run (copy / clone / sexual) | mean step | max step | run | RSS MB start / t=30 / end (max) | final island occupancy |
|---|---|---|---|---|---|---|---|
| asexual nemo v3sem | 8 | 27 / 17372 / 0 | 9.6 ms | 45 ms | 2.6 s | 126.6 / 127.2 / 127.9 | 91.2 |
| asexual nemo default | 3 | 0 / 18849 / 0 | 11.2 ms | 62 ms | 3.1 s | 124.9 / 125.7 / 126.9 | 69.0 |
| asexual nemo novar | 5 | 29 / 17095 / 0 | 8.1 ms | 35 ms | 2.3 s | 124.5 / 125.2 / 126.2 | 51.2 |
| asexual nemo v3gen | 5 | 30 / 17543 / 0 | 8.4 ms | 96 ms | 2.5 s | 126.0 / 126.6 / 127.7 | 0.8 |
| asexual v3pure v3sem | 5 | - | 17.4 ms | 80 ms | 3.9 s | 112.4 / 113.6 / 114.2 | 0.0 |
| sexual nemo v3sem | 8 | 34 / 6448 / 6657 | 7.9 ms | 82 ms | 2.3 s | 126.0 / 126.6 / 126.8 | 0.1 |
| sexual nemo default | 3 | 0 / 6322 / 7008 | 8.2 ms | 35 ms | 2.5 s | 125.2 / 125.8 / 126.2 | 0.3 |
| sexual nemo novar | 5 | 28 / 7856 / 5006 | 7.9 ms | 70 ms | 2.3 s | 126.4 / 127.0 / 127.3 | 0.0 |
| sexual nemo v3gen | 5 | 31 / 7891 / 4799 | 6.7 ms | 75 ms | 1.8 s | 128.2 / 128.9 / 129.2 | 0.2 |
| sexual v3pure v3sem | 5 | - | 13.5 ms | 74 ms | 3.1 s | 111.8 / 112.6 / 112.8 | 0.2 |

### asexual branch: NemoBackend vs pure v3 (5 seeds each, mean +- SE)

| metric | NemoBackend v3sem | pure v3 | z | novar | z | v3gen | z |
|---|---|---|---|---|---|---|---|
| occupancy (time mean, of 1024 sites) | 116.349 +- 20.120 (n 5) | 154.043 +- 25.507 (n 5) | 1.2 | 121.845 +- 24.468 (n 5) | 0.9 | 132.468 +- 26.151 (n 5) | 0.6 |
| occupancy (last 30 steps) | 116.020 +- 33.825 (n 5) | 13.953 +- 0.979 (n 5) | **3.0** | 58.533 +- 13.829 (n 5) | **3.2** | 8.687 +- 1.877 (n 5) | 2.5 |
| sp2 share on islands (time mean) | 0.010 +- 0.007 (n 5) | 0.022 +- 0.003 (n 5) | 1.6 | 0.015 +- 0.007 (n 5) | 0.8 | 0.035 +- 0.012 (n 5) | 1.1 |
| final phenotype mean sp1 (islands) | 0.197 +- 0.014 (n 4) | - | - | 0.172 +- 0.019 (n 4) | - | 0.027 +- - (n 1) | - |
| final phenotype mean sp2 (islands) | - | - | - | - | - | 1.026 +- - (n 1) | - |
| final phenotype mean sp1 (mainland1) | 0.005 +- 0.003 (n 5) | 0.001 +- 0.001 (n 5) | 1.3 | 0.001 +- 0.001 (n 5) | 0.6 | -0.000 +- 0.001 (n 5) | 2.2 |
| final phenotype mean sp2 (mainland2) | 1.000 +- 0.001 (n 5) | 0.999 +- 0.001 (n 5) | 1.0 | 1.001 +- 0.001 (n 5) | 2.3 | 1.000 +- 0.000 (n 5) | 2.1 |
| mainland1 size (time mean) | 400.000 +- 0.000 (n 5) | 400.000 +- 0.000 (n 5) | - | 400.000 +- 0.000 (n 5) | - | 400.000 +- 0.000 (n 5) | - |
| mainland2 size (time mean) | 400.000 +- 0.000 (n 5) | 400.000 +- 0.000 (n 5) | - | 400.000 +- 0.000 (n 5) | - | 400.000 +- 0.000 (n 5) | - |

| island occupancy / sp1 phenotype | t=0 | t=29 | t=59 | t=89 | t=119 | t=149 | founder phenotype SD (ml1 sp1, ml2 sp2) |
|---|---|---|---|---|---|---|---|
| nemo v3sem | 0.2 / -0.042 | 61.2 / -0.012 | 141.2 / 0.006 | 220.0 / 0.019 | 179.4 / 0.104 | 98.8 / 0.197 | 0.0511 / 0.0533 |
| nemo novar | 0.2 / 0.022 | 73.8 / -0.001 | 214.2 / 0.002 | 283.4 / 0.008 | 119.4 / 0.061 | 51.2 / 0.172 | 0.0248 / 0.0246 |
| nemo v3gen | 0.2 / 0.022 | 82.0 / -0.001 | 235.6 / 0.002 | 367.4 / 0.004 | 89.2 / 0.005 | 0.8 / 0.027 | 0.0248 / 0.0246 |
| v3pure v3sem | 0.2 / -0.008 | 100.6 / -0.002 | 270.0 / 0.005 | 389.6 / 0.007 | 120.0 / 0.018 | 0.0 / - | 0.025 / 0.0243 |
| environment mean | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.4 |  |

### sexual branch: NemoBackend vs pure v3 (5 seeds each, mean +- SE)

| metric | NemoBackend v3sem | pure v3 | z | novar | z | v3gen | z |
|---|---|---|---|---|---|---|---|
| occupancy (time mean, of 1024 sites) | 5.275 +- 2.217 (n 5) | 10.944 +- 8.160 (n 5) | 0.7 | 8.701 +- 6.319 (n 5) | 0.2 | 4.668 +- 2.787 (n 5) | 0.7 |
| occupancy (last 30 steps) | 0.227 +- 0.041 (n 5) | 3.733 +- 3.342 (n 5) | 1.0 | 0.180 +- 0.043 (n 5) | 1.1 | 0.320 +- 0.065 (n 5) | 1.0 |
| sp2 share on islands (time mean) | 0.113 +- 0.028 (n 5) | 0.163 +- 0.061 (n 5) | 0.7 | 0.108 +- 0.025 (n 5) | 0.8 | 0.142 +- 0.029 (n 5) | 0.3 |
| final phenotype mean sp1 (islands) | -0.004 +- - (n 1) | - | - | - | - | -0.009 +- - (n 1) | - |
| final phenotype mean sp2 (islands) | - | 1.000 +- - (n 1) | - | - | - | - | - |
| final phenotype mean sp1 (mainland1) | -0.002 +- 0.004 (n 5) | 0.002 +- 0.000 (n 5) | 1.0 | -0.000 +- 0.000 (n 5) | **3.7** | 0.001 +- 0.000 (n 5) | 1.5 |
| final phenotype mean sp2 (mainland2) | 0.999 +- 0.002 (n 5) | 0.997 +- 0.001 (n 5) | 0.7 | 1.000 +- 0.001 (n 5) | 2.8 | 0.999 +- 0.001 (n 5) | 1.7 |
| mainland1 size (time mean) | 400.000 +- 0.000 (n 5) | 400.000 +- 0.000 (n 5) | - | 400.000 +- 0.000 (n 5) | - | 400.000 +- 0.000 (n 5) | - |
| mainland2 size (time mean) | 400.000 +- 0.000 (n 5) | 399.997 +- 0.002 (n 5) | 1.6 | 400.000 +- 0.000 (n 5) | 1.6 | 400.000 +- 0.000 (n 5) | 1.6 |

| island occupancy / sp1 phenotype | t=0 | t=29 | t=59 | t=89 | t=119 | t=149 | founder phenotype SD (ml1 sp1, ml2 sp2) |
|---|---|---|---|---|---|---|---|
| nemo v3sem | 0.2 / 0.020 | 34.4 / 0.013 | 9.0 / 0.038 | 0.2 / 0.034 | 0.2 / -0.013 | 0.2 / -0.004 | 0.0511 / 0.0533 |
| nemo novar | 0.2 / -0.005 | 33.6 / -0.001 | 14.0 / 0.011 | 13.0 / 0.015 | 0.2 / -0.008 | 0.0 / - | 0.0248 / 0.0246 |
| nemo v3gen | 0.2 / -0.005 | 30.2 / 0.013 | 3.8 / 0.028 | 0.8 / 0.016 | 0.0 / - | 0.2 / -0.009 | 0.0248 / 0.0246 |
| v3pure v3sem | 0.2 / 0.029 | 23.2 / 0.002 | 15.8 / 0.013 | 14.0 / 0.027 | 13.2 / 0.065 | 0.2 / - | 0.0248 / 0.0255 |
| environment mean | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.4 |  |

### is_remove=False (v3 semantics) vs default True, NemoBackend, 3 seeds (mean +- SE)

| branch | metric | is_remove=False | default True |
|---|---|---|---|
| asexual | occupancy (time mean, of 1024 sites) | 126.447 +- 32.978 | 166.558 +- 54.192 |
| asexual | occupancy (last 30 steps) | 103.100 +- 52.054 | 81.411 +- 63.289 |
| asexual | sp2 share on islands (time mean) | 0.004 +- 0.002 | 0.011 +- 0.005 |
| asexual | final phenotype mean sp1 (islands) | 0.174 +- 0.001 | 0.133 +- 0.060 |
| asexual | final phenotype mean sp2 (islands) | - +- - | - +- - |
| asexual | final phenotype mean sp1 (mainland1) | 0.008 +- 0.004 | -0.006 +- 0.004 |
| asexual | final phenotype mean sp2 (mainland2) | 1.001 +- 0.001 | 1.009 +- 0.005 |
| asexual | mainland1 size (time mean) | 400.000 +- 0.000 | 399.884 +- 0.008 |
| asexual | mainland2 size (time mean) | 400.000 +- 0.000 | 399.896 +- 0.012 |
| sexual | occupancy (time mean, of 1024 sites) | 3.931 +- 3.025 | 8.769 +- 6.839 |
| sexual | occupancy (last 30 steps) | 0.256 +- 0.068 | 0.456 +- 0.149 |
| sexual | sp2 share on islands (time mean) | 0.155 +- 0.019 | 0.148 +- 0.070 |
| sexual | final phenotype mean sp1 (islands) | - +- - | - +- - |
| sexual | final phenotype mean sp2 (islands) | - +- - | 0.923 +- - |
| sexual | final phenotype mean sp1 (mainland1) | 0.002 +- 0.005 | -0.005 +- 0.005 |
| sexual | final phenotype mean sp2 (mainland2) | 1.002 +- 0.002 | 0.998 +- 0.003 |
| sexual | mainland1 size (time mean) | 400.000 +- 0.000 | 399.907 +- 0.029 |
| sexual | mainland2 size (time mean) | 400.000 +- 0.000 | 399.880 +- 0.015 |

**Reading the simulation tables**

- asexual branch, occupancy (last 30 steps): NemoBackend 116 +- 34 vs pure v3 14 +- 1 (z 3.0). The trajectories show where it comes from: NemoBackend island sp1 follows the environmental shift (sp1 phenotype 0.006 at t=59 -> 0.197 at t=149, environment 0.2 -> 0.4) and keeps islands occupied; v3 sp1 stays near 0 (0.018 at t=119) and the islands are empty at the end. Occupancy NemoBackend 0.2 -> 61.2 -> 141.2 -> 220.0 -> 179.4 -> 98.8; v3 0.2 -> 100.6 -> 270.0 -> 389.6 -> 120.0 -> 0.0.
- the two controls separate genetics from the copy / clone code: without founder genetic variance (`novar`: founder phenotype SD 0.0248, as v3 0.025) NemoBackend still adapts (sp1 0.172 at the end, z 3.2); with v3-sized mutation effects as well (`v3gen`: quanti_mutation_variance 0.025^2 instead of 0.01, i.e. the per-birth mutational variance of v3's bit flips) NemoBackend reproduces v3 -- occupancy 0.2 -> 82.0 -> 235.6 -> 367.4 -> 89.2 -> 0.8 vs v3 0.2 -> 100.6 -> 270.0 -> 389.6 -> 120.0 -> 0.0, every metric within 2.5 SE. So the difference is the default Nemo genetic parameters of the NemoBackend runs (quanti_init_variance 0.002025, quanti_mutation_variance 0.01, as in test/metaibm-nemo-run/simulate.py), not the migration; the clone / copy semantics themselves are pinned by the unit tests.
- sexual branch, mainland1 size (time mean): NemoBackend 400.000 (z -) vs v3 400.000. Mixed births on the mainlands, measured in the 8 NemoBackend v3sem runs: of 2400 calls 2400 filled every empty site, 0 left a site empty although asex + sex >= empty (the v3 split's np.around undershoot; 0 expected with the current split, see the hab_mixed_reproduce_germinate finding) and 0 had too few parents; v3gen control 400.000.
- the scaled sexual (mixed) islands stay almost empty in every variant (time-mean occupancy NemoBackend 5.3, v3 10.9 of 1024): single colonists rarely find a mate or fit well enough to clone, so this branch's island statistics are weak; its births (6448 clones + 6657 sexual offspring per run) are still exercised, and the mixed births are tested directly in test_core.py (births parts).
- is_remove=False vs default True: with False the colonists are copies (27 make_copy per run) and the mainlands stay full; with True (0 copies) colonists leave the mainland, which then averages 399.88 sites instead of 400.
- memory: RSS is flat within each run (start / t=30 / end in the performance table); teardown minted == recycled in every NemoBackend run.

## Regression: existing suites in throw-away clones

- `/tmp/mibg` (migrated): clone of `/Users/jianhaolin/Library/CloudStorage/OneDrive-UniversityofHelsinki/MetaIBM4Nemo_workspace/MetaIBM4Nemo`, HEAD `140acc4` (worktree HEAD `140acc4`); 6 changed files copied, md5 == working tree: **True**; `git status`: `M extension/global_habitat_network.py`, ` M metaibm/habitat.py`, ` M metaibm/indfactory.py`, ` M metaibm/individual.py`, ` M metaibm/metacommunity.py`, ` M metaibm/patch.py`
- `/tmp/mibg-base` (base): clone of `/Users/jianhaolin/Library/CloudStorage/OneDrive-UniversityofHelsinki/MetaIBM4Nemo_workspace/MetaIBM4Nemo`, HEAD `140acc4` (worktree HEAD `140acc4`); no copy (pristine 140acc4; md5 differ from the working tree as expected: True); `git status`: clean

| suite / script | pristine 140acc4 | migrated | (1) known at 140acc4 | (2) expected from the migration | (3) other | only at 140acc4 (masked / not reached) |
|---|---|---|---|---|---|---|
| nemo-indfactory/smoke.py | rc 0, 0.6 s | rc 0, 0.5 s | 0 | 0 | 0 | 0 |
| nemo-indfactory/accept_genetics.py | rc 0, 0.6 s | rc 0, 0.5 s | 0 | 0 | 0 | 0 |
| nemo-indfactory/accept_guards.py | rc 0, 0.6 s | rc 0, 0.6 s | 0 | 0 | 0 | 0 |
| nemo-indfactory/accept_lifecycle.py | rc 0, 0.6 s | rc 0, 264.2 s | 0 | 0 | 0 | 0 |
| nemoback-habitat/run_all.py | rc 0, 19.9 s | rc 1, 1.1 s | 0 | 1 | 0 | 0 |
| nemoback-habitat/test_bookkeeping.py (direct) | rc 0, 5.8 s | rc 1, 0.8 s | 0 | 1 | 0 | 0 |
| nemoback-habitat/test_nemo.py (direct) | rc 0, 13.8 s | rc 0, 14.0 s | 0 | 0 | 0 | 0 |
| nemoback-patch/run_all.py | rc 1, 37.0 s | rc 1, 4.6 s | 0 | 1 | 0 | 1 |
| nemoback-patch/test_bookkeeping.py (direct) | rc 1, 0.7 s | rc 1, 0.8 s | 0 | 1 | 0 | 1 |
| nemoback-patch/test_nemo.py (direct) | rc 0, 12.4 s | rc 1, 13.4 s | 0 | 3 | 0 | 0 |
| nemoback-patch/test_chain.py (direct) | rc 0, 0.6 s | rc 0, 0.7 s | 0 | 0 | 0 | 0 |
| nemoback-metacommunity/test_metacommunity.py all | rc 1, 92.2 s | rc 1, 57.6 s | 1 | 5 | 0 | 2 |
| nemoback-extension/test_extension.py all --quick | rc 0, 291.9 s | rc 0, 210.4 s | 22 | 13 | 0 | 4 |

- nemo-indfactory smoke digest: `2b9b330388a3bfcdc7239fb45f7903076ebfaded` (expected `2b9b330388a3bfcdc7239fb45f7903076ebfaded`): **OK**

### Every failure of the migrated clone

| suite | cat | why | kind | failure (test frame / code / error, or check name, or summary line) |
|---|---|---|---|---|
| nemoback-habitat/run_all.py | 2 | summary line of a script / part whose own failures are listed separately | line | smoke.py FAIL 1.0s |
| nemoback-habitat/test_bookkeeping.py (direct) | 2 | asserts the removed is_remove=False NotImplementedError guard | traceback | test_bookkeeping.py:<module> / try: getattr(h, meth)(is_remove=False); raise AssertionError('is_remove=False accepted') / AttributeError: 'fake_indi' object has no attribute 'make_copy' |
| nemoback-patch/run_all.py | 2 | summary line of a script / part whose own failures are listed separately | line | smoke.py FAIL 0.6s |
| nemoback-patch/test_bookkeeping.py (direct) | 2 | asserts the removed is_remove=False NotImplementedError guard | traceback | test_bookkeeping.py:<module> / try: getattr(p, fn)(0.5, False); raise AssertionError('%s accepted is_remove=False' % fn) / AttributeError: 'fake_indi' object has no attribute 'make_copy' |
| nemoback-patch/test_nemo.py (direct) | 2 | asserts that asex / mix still raise TypeError (now migrated) | traceback | test_nemo.py:section_edges_known / if 'calculation' not in n: assert r.startswith('TypeError'), (n, r) / AssertionError: ('patch_asexual_birth_germinate', 'returned 21') |
| nemoback-patch/test_nemo.py (direct) | 2 | summary line of a script / part whose own failures are listed separately | line | ---- section edges_known FAILED |
| nemoback-patch/test_nemo.py (direct) | 2 | summary line of a script / part whose own failures are listed separately | line | FAILED sections: ['edges_known'] |
| nemoback-metacommunity/test_metacommunity.py all | 2 | asserts that asex / mix still raise TypeError (now migrated) | traceback | test_metacommunity.py:section_probes / if 'calculation' not in n: assert r.startswith(('TypeError', 'NotImplementedError')), (n, r) / AssertionError: ('meta_mainland_asexual_birth_mutate_germinate', "returned '[Birth process] in mainland: there are 4 individuals germina'") |
| nemoback-metacommunity/test_metacommunity.py all | 2 | asserts the guard / the reserved else branch of the previous round | traceback | test_metacommunity.py:bookkeeping_main / assert sorted(n for n, v in r3.items() if v['guard_first']) == sorted(GUARDED_META), sorted(n for n, v in r3.items() if v['guard_first']) / AssertionError: [] |
| nemoback-metacommunity/test_metacommunity.py all | 2 | asserts the removed is_remove=False NotImplementedError guard | traceback | test_metacommunity.py:guard / try: rng.choice(GUARDED)(ml, isl); raise AssertionError('is_remove=False accepted') / AssertionError: is_remove=False accepted |
| nemoback-metacommunity/test_metacommunity.py all | 1 | also fails in the pristine 140acc4 clone | line | bookkeeping            FAIL 0.7s |
| nemoback-metacommunity/test_metacommunity.py all | 2 | summary line of a script / part whose own failures are listed separately | line | nemo                   FAIL 28.0s |
| nemoback-metacommunity/test_metacommunity.py all | 2 | summary line of a script / part whose own failures are listed separately | line | chain                  FAIL 19.8s |
| nemoback-extension/test_extension.py all --quick | 2 | Step 0 of the extension re-runs a lower suite: its failures are classified under that suite above | check | baseline: Step 0 test/nemoback-habitat run_all.py (smoke.py FAIL) / {'smoke.py': ['smoke.py FAIL 0.6s']} |
| nemoback-extension/test_extension.py all --quick | 2 | Step 0 of the extension re-runs a lower suite: its failures are classified under that suite above | check | baseline: Step 0 test/nemoback-patch run_all.py --no-baseline (smoke.py FAIL) / {'smoke.py': ['smoke.py FAIL 0.6s']} |
| nemoback-extension/test_extension.py all --quick | 2 | Step 0 of the extension re-runs a lower suite: its failures are classified under that suite above | check | baseline: Step 0 test/nemoback-metacommunity test_metacommunity.py all (smoke PASS, bookkeeping FAIL, nemo FAIL, chain FAIL) / {'bookkeeping': ["1 R1: 14 metacommunity functions take is_remove, every default is True (original: ['False'])", 'Traceback (most recent call last):', '  File "/private/tmp/mibg/test/nemoback-metacommunity/test_metacommunity.py", line 1504, in <module>', '    PARTS[part]()', '  File "/private |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | traceback | dispersal: common.py:original_metacommunity / assert body.count(old) == 1, 'patch for %s matches %d times: %r' % (k, body.count(old), old) / AssertionError: patch for dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool matches 0 times: 'migrants_to_j_indi_object_ls += random.sample(h_i_object.offspring_pool, migrants_i_j_num)' |
| nemoback-extension/test_extension.py all --quick | 2 | asserts the removed is_remove=False NotImplementedError guard | check | fuzz: G 4 random operations (3 kinds, 3 processes x 10 runs x 110 ops): invariants after each, teardown minted == recycled after every run, no unexpected exception / ['run 0 op 1 NETWORK offspring is_remove=False (guard): AssertionError: no NotImplementedError', 'run 0 op 1 NETWORK offspring is_remove=False (guard): AssertionError: no NotImplementedError', 'run 0 op 2 NETWORK offspring is_remove=False (guard): Assert |
| nemoback-extension/test_extension.py all --quick | 2 | asserts the removed is_remove=False NotImplementedError guard | check | fuzz: G guarded is_remove=False calls (offspring network): NotImplementedError + full snapshot unchanged / 0 |
| nemoback-extension/test_extension.py all --quick | 1 | also fails in the pristine 140acc4 clone | check | marker: D1 is_remove=False == v3 (copy) on every one of 60 seeds (state, random + np.random, log, which markers reach which pool in which order) / {'digest': 11, 'rng': 40, 'ret': 43, 'new': 41} |
| nemoback-extension/test_extension.py all --quick | 1 | also fails in the pristine 140acc4 clone | check | marker: D1 reference sanity: v3 copy removes nothing /  |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | traceback | perf: common.py:original_metacommunity / assert body.count(old) == 1, 'patch for %s matches %d times: %r' % (k, body.count(old), old) / AssertionError: patch for dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool matches 0 times: 'migrants_to_j_indi_object_ls += random.sample(h_i_object.offspring_pool, migrants_i_j_num)' |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | traceback | registry: common.py:original_metacommunity / body = _cut(srcs, fns[k], '    ') / KeyError: 'global_habitat_id_idx_registry' |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | check | sim: E run equiv_net_v3 produced its result file /  = C.original_metacommunity(ext_names=[C.ORIGINAL_REGISTRY, MARK], meta_names=['add_patch'])           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   File "/private/tmp/mibg/test/nemoback-extension/common.py", line 91, in original_metacommunity     body = |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | check | sim: E run equiv_grfe_v3 produced its result file /  = C.original_metacommunity(ext_names=[C.ORIGINAL_REGISTRY, MARK], meta_names=['add_patch'])           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   File "/private/tmp/mibg/test/nemoback-extension/common.py", line 91, in original_metacommunity     body = |
| nemoback-extension/test_extension.py all --quick | 2 | asserts the previous round's diff / working tree / tag set | check | static: A1 every changed line of extension / add_patch / docs is one of the agreed changes (or blank / comment / layout) / [['extension/global_habitat_network.py', '-', "''' dispersal from habitat_i offspring_pool to habitat_j immigrant_pool in the global habitat network '''", 'UNEXPECTED'], ['extension/global_habitat_network.py', '+', "''' NemoBackend: dispersal from habitat_i offspring_pool to habitat_j immigrant_p |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | check | static: A1 agreed change present: registry rename / 0 |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | check | static: A1 agreed change present: record stores habitat.global_index / 0 |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | check | static: A1 agreed change present: distance matrix recomputed inside the registry / 0 |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | check | static: A1 agreed change present: is_remove in the signature / 0 |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | check | static: A1 agreed change present: sampling / 0 |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | check | static: A1 agreed change present: install name / 0 |
| nemoback-extension/test_extension.py all --quick | 1 | extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location | check | static: A1 agreed change present: TYPE_CHECKING import / 0 |
| nemoback-extension/test_extension.py all --quick | 1 | also fails in the pristine 140acc4 clone | check | static: A1 working tree: modified == extension, metacommunity.py, docs (+ the two Step 0 test files); untracked == test/nemoback-extension/ only / [' M extension/global_habitat_network.py', ' M metaibm/habitat.py', ' M metaibm/indfactory.py', ' M metaibm/individual.py', ' M metaibm/metacommunity.py', ' M metaibm/patch.py', ' M test/nemoback-habitat/TEST_REPORT.md', ' M test/nemoback-metacommunity/TEST_REPORT.md', ' M |
| nemoback-extension/test_extension.py all --quick | 2 | asserts the guard / the reserved else branch of the previous round | check | static: A3 offspring: guard is the first statement (before patch_num < 2 and any RNG / state change) / if self.patch_num < 2: |
| nemoback-extension/test_extension.py all --quick | 2 | asserts the guard / the reserved else branch of the previous round | check | static: A3 offspring: guard text == metacommunity.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool / None |
| nemoback-extension/test_extension.py all --quick | 2 | asserts the guard / the reserved else branch of the previous round | check | static: A3 offspring: else branch is `pass  # unreachable (guard above); TODO(copy_from) ...` /  |
| nemoback-extension/test_extension.py all --quick | 2 | asserts the guard / the reserved else branch of the previous round | check | static: A3 offspring: no random.sample left in code (# comments stripped first) /  |
| nemoback-extension/test_extension.py all --quick | 2 | asserts the previous round's diff / working tree / tag set | check | static: A7 no NemoBackend tag on the registry / the two network dispersals (they touch no Nemo object) / {'global_habitat_network_registry': False, 'dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool': True, 'dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool': False} |
| nemoback-extension/test_extension.py all --quick | 2 | asserts the previous round's diff / working tree / tag set | check | static: A7 no function of the extension is tagged NemoBackend / ['dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool'] |
| nemoback-extension/test_extension.py all --quick | 1 | also fails in the pristine 140acc4 clone | line | [all] baseline rc=1 63.3s ['BASELINE FAILED: 8/11 checks, 1 findings (62.67 s)'] |
| nemoback-extension/test_extension.py all --quick | 1 | also fails in the pristine 140acc4 clone | line | [all] static rc=1 2.4s ['STATIC FAILED: 14/29 checks, 7 findings (1.65 s)'] |
| nemoback-extension/test_extension.py all --quick | 1 | also fails in the pristine 140acc4 clone | line | [all] perf rc=1 0.7s ['PERF FAILED: 0/0 checks, 0 findings (0.07 s)'] |
| nemoback-extension/test_extension.py all --quick | 1 | also fails in the pristine 140acc4 clone | line | [all] dispersal rc=1 1.1s ['DISPERSAL FAILED: 0/0 checks, 0 findings (0.09 s)'] |
| nemoback-extension/test_extension.py all --quick | 1 | also fails in the pristine 140acc4 clone | line | [all] registry rc=1 1.4s ['REGISTRY FAILED: 5/5 checks, 0 findings (0.34 s)'] |
| nemoback-extension/test_extension.py all --quick | 2 | summary line of a script / part whose own failures are listed separately | line | [all] fuzz rc=1 3.5s ['FUZZ FAILED: 0/2 checks, 0 findings (2.51 s)'] |
| nemoback-extension/test_extension.py all --quick | 1 | also fails in the pristine 140acc4 clone | line | [all] marker rc=1 36.1s ['MARKER FAILED: 15/17 checks, 2 findings (35.02 s)'] |
| nemoback-extension/test_extension.py all --quick | 1 | also fails in the pristine 140acc4 clone | line | [all] sim rc=1 141.6s ['SIM FAILED: 32/34 checks, 1 findings (140.48 s)'] |

### Failing at 140acc4 but not in the migrated run (the migrated script stopped earlier, or the failure point moved)

| suite | failure at 140acc4 |
|---|---|
| nemoback-patch/run_all.py | test_bookkeeping.py    FAIL 0.7s |
| nemoback-patch/test_bookkeeping.py (direct) | test_bookkeeping.py:<module> / assert same, tag + ' RNG streams differ' / AssertionError: O->I seed=0 rate=0.5 vs patch RNG streams differ |
| nemoback-metacommunity/test_metacommunity.py all | test_bookkeeping.py:<module> / assert same, tag + ' RNG streams differ' / AssertionError: O->I seed=0 rate=0.5 vs patch RNG streams differ |
| nemoback-metacommunity/test_metacommunity.py all | test_metacommunity.py:bookkeeping_main / assert json.dumps(a, default=str) == json.dumps(b, default=str), 'helpers differ from the original seed=%d' % seed / AssertionError: helpers differ from the original seed=0 |
| nemoback-extension/test_extension.py all --quick | baseline: Step 0 test/nemoback-patch run_all.py --no-baseline (smoke.py PASS, test_bookkeeping.py FAIL, test_nemo.py PASS, test_chain.py PASS) / {'test_bookkeeping.py': ['1 six clear-ups: 180 cases; return == pool size, only that pool emptied, recycled == former contents (in order), markers recycle  |
| nemoback-extension/test_extension.py all --quick | baseline: Step 0 test/nemoback-metacommunity test_metacommunity.py all (smoke PASS, bookkeeping FAIL, nemo PASS, chain PASS) / {'bookkeeping': ['8 patch_num < 2: 4 log strings now contain the metacommunity name (v3: literal "%s"), e.g. \'[Dispersal among patches] in solo: patch_num < 2, there are 0  |
| nemoback-extension/test_extension.py all --quick | static: A1 agreed change present: guard / 0 |
| nemoback-extension/test_extension.py all --quick | static: A1 agreed change present: else branch / 0 |

**Reading the regression tables**

- 48 failures of the migrated clone: 23 already fail at 140acc4 (1), 25 assert what this migration changed on purpose (2), **0 other (3)**.
- coverage the old suites lose until their tests are updated: these scripts / sections stop at an expected (2) failure, so their later checks did not run against the migrated code -- `nemoback-habitat/test_bookkeeping.py (direct)` at `test_bookkeeping.py:<module>` (`try: getattr(h, meth)(is_remove=False); raise AssertionError('is_remove=False accepted')`); `nemoback-patch/test_bookkeeping.py (direct)` at `test_bookkeeping.py:<module>` (`try: getattr(p, fn)(0.5, False); raise AssertionError('%s accepted is_remove=False' % fn)`); `nemoback-patch/test_nemo.py (direct)` at `test_nemo.py:section_edges_known` (`if 'calculation' not in n: assert r.startswith('TypeError'), (n, r)`); `nemoback-metacommunity/test_metacommunity.py all` at `test_metacommunity.py:section_probes` (`if 'calculation' not in n: assert r.startswith(('TypeError', 'NotImplementedError')), (n, r)`); `nemoback-metacommunity/test_metacommunity.py all` at `test_metacommunity.py:bookkeeping_main` (`assert sorted(n for n, v in r3.items() if v['guard_first']) == sorted(GUARDED_META), sorted(n for n, v in r3.items() if `); `nemoback-metacommunity/test_metacommunity.py all` at `test_metacommunity.py:guard` (`try: rng.choice(GUARDED)(ml, isl); raise AssertionError('is_remove=False accepted')`). The behaviour they would have checked is covered by the new tests above.
- the known 140acc4 failures of patch test_bookkeeping (v3 class comparison: "RNG streams differ") and metacommunity bookkeeping ("helpers differ from the original") are not reached in the migrated run (the scripts stop earlier, at the guard assertions); they are listed under "Failing at 140acc4 but not in the migrated run".
- extension suite: `ORIGINAL_EXT_COMMIT = git log -1 -- extension/global_habitat_network.py`, which is `53152d5` (the migration commit) at 140acc4 instead of `1356387`; its perf / dispersal / registry parts stop inside `original_metacommunity` in BOTH trees (same failure point, category 1). `test_extension.py all` exits 0 although parts fail.
- the expected (2) failures need test updates, not source changes: drop the NotImplementedError / TypeError expectations, give hab_common.fake_indi a make_copy(), and refresh the diff / tag / working-tree assertions of the previous round.

## Notes on the method

- **v3 comparison (A)**: the v3 bodies are compiled verbatim from git into a subclass of the CURRENT class (habitat / patch `baace34`, metacommunity `7a1adbb`, extension `1356387`, module-level function indented into the class), so both twins run on the same lower layers; the twins are built from one seed and called with one seed. Fake individuals carry `.origin`, so a v3 alias and a current copy compare equal in the layout. The two direct-germination dispersals rarely overflow the empty sites on random worlds, so they are also run with the number matrix pinned (`[forced overflow]`, no RNG, same matrix for both twins) -- there the unplaced copies are recycled on ~185 / 200 seeds.
- **make_copy semantics**: `test_core.py unit` (part rng) rebuilds copies and clones in Python without `operator=` / `makeClonedOffsprg` (setters, `set_allele_value`, `set_value`, then the clone steps of `breed_cloning`) and gets the same IDs, fields, alleles and phenotypes record by record -- so `operator=` copies exactly what is listed (all fields except the ID; `_mother` / `_father` raw pointers too, never read by MetaIBM) and RAND is consumed in the same order (one Ve draw in `operator=`, one more for a clone after mutation).
- **test the tests**: broken variants (alias germination, remove instead of copy, a make_copy that drops the age) were monkeypatched in a scratch process; the copy-path harness flagged all three.
- **pure v3**: `git archive baace34 metaibm` alone is not v3 -- `baace34` itself introduced the Nemo-backed `individual.py` (+ `indfactory.py`); `/tmp/v3pure` therefore takes `metaibm/individual.py` from `baace34^` (`67a742f`) and `extension/` from `baace34` (== `1356387`), as the v3 metacommunity imports it.
- **regression classification**: category 1 = the same failure point (test frame + code line, or extension check name) also fails in the pristine 140acc4 clone; category 2 = the failure asserts what the migration intentionally changed (patterns listed in `run_all.py:EXPECTED`); category 3 = everything else. A script that aborts at its first failure hides its later checks: those parts are not exercised against the migrated code by the old suites (the new tests here cover the migrated behaviour).

