# disp-trait sampling of emigrants -- test report

Generated on 2026-11-11 at 13:41 by `test_disp_sampling.py` and `test_disp_lifecycle.py` in this folder.

**Scope.** Dispersal that samples emigrants by their disp value v (`fdisp`; `mdisp` for males if `sex_specific`). It is switched on by the new kwargs `use_disp_trait=False, sex_specific=False, is_rescale=True, disp_cost=None` of the dispersal functions.
A species can fix one sex's v with Nemo's two rows in `indfactory(..., metaibm_nemo_params={'dispersal_fixed_trait': 'female' | 'male', 'dispersal_fixed_rate': [0, 1]})`: that sex's v is the fixed rate, the other sex reads its trait as above. Unset (default), nothing changes.

| scheme | kwargs | who goes from source i to sink j |
|---|---|---|
| 1 (old MetaIBM, default) | `use_disp_trait=False` | e_ij individuals, uniform |
| 2 | `use_disp_trait=True` (`is_rescale=True`) | e_ij individuals, each draw P = v_k / sum(v) of the remaining |
| 3 | `is_rescale=False` | each individual leaves source_i (v_k) and goes to sink_j (m_ij): P = v_k · m_ij |
| 4 | `is_rescale=False, disp_cost=c` or `(c_fem, c_mal)` | scheme 3, then dies on the way with P = c (removed and recycled) |

Code under test is the working tree on top of HEAD `b3d71bc`, uncommitted:

- `metaibm/habitat.py`:
  - `hab_sample_emigrants` (line 181), whose v line reads `dispersal_fixed_trait` / `dispersal_fixed_rate`
  - `_hab_sample_by_number` (195)
  - `_hab_sample_by_rate` (201)
- `metaibm/patch.py`:
  - `patch_sample_emigrants` (148), whose v line reads `dispersal_fixed_trait` / `dispersal_fixed_rate`
  - `_patch_sample_by_number` (162)
  - `_patch_sample_by_rate` (168)
  - the 3 within functions (384, 442, 470)
- `metaibm/metacommunity.py`:
  - the 3 among functions (576, 617, 680)
  - the 3 `meta_dispersal_within_patch_...` wrappers (722, 731, 740)
- `extension/global_habitat_network.py`:
  - `dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool` (116)
- `metaibm/indfactory.py`:
  - `_check_dispersal_fixed_params` (76), called in `_boot` (86): ValueError for a `dispersal_fixed_trait` other than 'female' / 'male', or a `dispersal_fixed_rate` outside [0, 1]

**Side effects: none.**

- The runs used a copy of the working tree; its `metaibm/`, `extension/` and test files have the md5s listed below, the same as the working tree. Nothing in the repository was touched, and nothing was committed.
- `metaibm/` was imported read-only (`sys.dont_write_bytecode`, `-B`). This folder holds only these 3 files, with no `__pycache__`.
- Nemo output went to temp directories, which were removed afterwards.

## Environment

- HEAD is `b3d71bc`, on branch `features/developing`. `habitat.py`, `patch.py`, `metacommunity.py`, `global_habitat_network.py` and `indfactory.py` have the uncommitted changes listed above.
- Conda env `metaibm-nemo-test6`:
  - Python 3.11.18, numpy 2.4.6, `macOS-15.7.3-arm64`.
  - pynemo `_nemo.cpython-311-darwin.so` md5 `201e8c8968f2e8b2e323b4827e0a0112`. It has `TTrait.get_value()`, which `individual.get_trait_value` needs.
- Each check group, life cycle and simulation runs in its own child process with fixed seeds (`random.seed` + `np.random.seed` + `indfactory(seed=...)`). Nemo's RAND and `Individual::currentID` are process-global.
- Species: `trait_keys=('fdisp', 'mdisp')`, with these disp init rows:
  - spU: `disp_init_distribution uniform {{0, 0}}`
  - spFM: fdisp 0.2 / mdisp 0.8
  - spFixF / spFixM: spFM with `dispersal_fixed_trait` 'female' / 'male' and `dispersal_fixed_rate` 0.6
  - sp00: v = 0
  - spN: a species without any trait, for the KeyError check

| file | md5 |
|---|---|
| `metaibm/habitat.py` (working tree, under test) | `29b13e1845a5d2ea87fe39b3432030d5` |
| `metaibm/patch.py` (working tree, under test) | `713a94897c2025ce8ef44c4f2fe1f46e` |
| `metaibm/metacommunity.py` (working tree, under test) | `1bb5e8981365fa79c800195075291a7c` |
| `extension/global_habitat_network.py` (working tree, under test) | `08c8ee10384bebcabde44e43f0e0eec2` |
| `metaibm/indfactory.py` (working tree, under test) | `ecf3a9110c0841f7849f0bf205b845e4` |
| `test/nemoback-traits/disp-trait/test_disp_sampling.py` | `d0a4bb8e4bf15d23b9593751b0f02c36` |
| `test/nemoback-traits/disp-trait/test_disp_lifecycle.py` | `7cc1230dbb9cbad16e5acdb43baa1c46` |

Run:

```
cd test/nemoback-traits/disp-trait
PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test6/bin/python -B test_disp_sampling.py    # ~14 s, exit 1 on any FAIL
PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test6/bin/python -B test_disp_lifecycle.py   # ~78 s
```

To run one group, append its name, e.g. `test_disp_sampling.py fixed` or `test_disp_lifecycle.py regression_head`.

## Oracles

All oracles live in the test code only.

- **frequency**: repeated single calls (up to 30 000) compared with v_k / sum(v), v_k · m_ij and the share of each sex, within 4 standard errors.
- **per-individual expectation**: the real among function, with one source and 3 destinations in loop order (patch2 → patch4), m_1j = (0.40276, 0.08987, 0.00738) at m = 0.5.
  - With a_kj = v_k · m_1j, `is_remove=True`:
    - P(→ j) = a_kj (1 − c) Π_{j' before j} (1 − a_kj')
    - P(dead) = c Σ_j a_kj Π_{j' before j} (1 − a_kj')
    - P(stay) = Π_j (1 − a_kj)
  - `is_remove=False`: P(copy at j) = a_kj (1 − c), independently for each j.
  - Each count is summed over 40 000 individuals and must have |z| < 4.
- **digest**: md5 of the full state after every step:
  - residents with id, species, fdisp and mdisp;
  - the ids in every offspring, immigrant and dormancy pool.
- **HEAD copy**: the working tree's `metaibm/` + `extension/` in a temp dir, with the 4 changed files taken from `git show HEAD:<path>`.
- **bookkeeping**: every individual minted and not recycled sits in exactly one place, a microsite or one pool. Ids are unique, and their count equals Σ(minted − recycled).
- **fixed rate**: 3000 females + 3000 males per species in one offspring pool, one call at m_ij = 1; each count of arrived / dead individuals per species and sex is Binomial(3000, p) and must have |z| < 4.

## Results

`test_disp_sampling.py`: PASS 26, FAIL 0, NOTE 6. `test_disp_lifecycle.py`: PASS 31, FAIL 0, NOTE 5.
Every number of the groups `units`, `exact`, `order` and of `test_disp_lifecycle.py` is the same as before `dispersal_fixed_trait` / `dispersal_fixed_rate` were added (2026-11-09 report).

### Per-pair API (`units`)

| status | check | measured |
|---|---|---|
| PASS | U1 scheme 2, e_ij = 1: P = v_k / sum(v); the emigrant leaves the pool | 10 individuals, 30 000 draws, all within 4 SE (e.g. 0.1823 / 0.183, 0.0175 / 0.0166) |
| PASS | U2 scheme 3: P(k → j) = v_k · m_ij (m 0.4) | e.g. 0.2826 / 0.2828, 0.2703 / 0.2665 |
| PASS | U3 scheme 4, (c_fem 0.2, c_mal 0.6): arrive / dead (recycled) / pool size | female arrive 614 (exp 603), male arrive 292 (exp 300), dead 575 (exp 601), pool left 4519 = 6000 − 906 − 575 |
| PASS | U4 `is_remove=False`: pool untouched, copies only for survivors, nothing recycled | scheme 4 copies 526 (exp 499), minted = copies, recycled 0; scheme 2 copies 5 |
| PASS | U5 `sex_specific=True`: share of males 0.8 | 0.8014 |
| PASS | U5 `sex_specific=False`: share of males 0.5 | 0.5060 |
| PASS | U6 v = 0 never drawn; e = 8 > 5 with v > 0 → 5; empty pool; e = 0 | drawn 5, all v > 0; v = 0 drawn by scheme 3: 0 |
| PASS | U7 `patch_sample_emigrants`: all habitats, offspring + dormancy pools | drawn 25, left 15, no overlap |
| PASS | U8 species without fdisp / mdisp | KeyError, pool untouched |

### Exact probabilities through the among function (`exact`, plan items 5 / 6)

| status | check | largest \|z\| |
|---|---|---|
| PASS | scheme 3, `is_remove=True`: P(→ j), P(stay) | 1.49 (stay, female) |
| PASS | scheme 4, c 0.3, `is_remove=True`: P(→ j), P(dead), P(stay) | 1.86 (patch2, male) |
| PASS | scheme 4, (c_fem 0.1, c_mal 0.5), `is_remove=True` | 2.73 (stay, male) |
| PASS | scheme 3, `is_remove=False`: copies at each j; copies at both patch2 and patch3 (independence) | 1.60; both: obs 427, exp 425.1 |
| PASS | scheme 4, c 0.3, `is_remove=False`: copies only for survivors | 1.38; both: obs 227, exp 207.2 |

### Destination order (`order`, plan item 8)

The among loop serves destinations in order. With `is_remove=True`, a later destination draws from a pool that the earlier ones have already thinned. Mean v of the emigrants per destination, at m = 0.5 (pool mean v ≈ 0.50):

| status | scheme | patch2 (first) | patch3 | patch4 |
|---|---|---|---|---|
| PASS | 1 (no difference expected) | 0.503 ± 0.002 | 0.500 ± 0.003 | 0.494 ± 0.012 |
| NOTE | 2 | 0.567 ± 0.001 | 0.533 ± 0.003 | 0.523 ± 0.011 |
| NOTE | 3 | 0.579 ± 0.002 | 0.562 ± 0.005 | 0.539 ± 0.020 |
| NOTE | 4, c 0.3 | 0.590 ± 0.002 | 0.563 ± 0.006 | 0.573 ± 0.021 |
| NOTE | 3, `is_remove=False` (each j draws from the full pool) | 0.580 ± 0.002 | 0.583 ± 0.004 | 0.599 ± 0.015 |

The leave rate at m = 1 is checked as follows:

- **PASS**: with a single scale, scheme 3, P(stay) = Π_j (1 − v_k m_1j). Observed 21 852 stays, expected 21 833.8, z +0.20.
- **NOTE**: the overall P(leave) is 0.454, against a target of Σv/N = 0.500. For v ≥ 0.8 it is 0.748, against a target of 0.867.

### Fixed rate per sex (`fixed`)

spFM founders: fdisp 0.2, mdisp 0.8, no mutation; the fixed rate is 0.6. Scheme 3 at m_ij = 1 unless the check says otherwise.

| status | check | measured (obs / exp, z) |
|---|---|---|
| PASS | F1 fixed 'female': the female leaves with 0.6 (not fdisp 0.2), the male with mdisp 0.8 | female 1835 / 1800, +1.30; male 2385 / 2400, −0.68 |
| PASS | F2 fixed 'male': the male leaves with 0.6 (not mdisp 0.8), the female with fdisp 0.2 | female 609 / 600, +0.41; male 1804 / 1800, +0.15 |
| PASS | F3 fixed 'female', `sex_specific=False`: the male reads fdisp 0.2 | female 1780 / 1800, −0.75; male 596 / 600, −0.18 |
| PASS | F4 fixed 'female', `disp_cost` 0.5: the fixed sex pays the cost too | female arrived 883 / 900, dead 906 / 900; male arrived 1194 / 1200, dead 1213 / 1200; all \|z\| ≤ 0.68 |
| PASS | F5 scheme 2, e_ij = 1: the fixed rate is the draw weight, share of males 0.8 / (0.8 + 0.6) | 0.5733 / 0.5714 (20 000 draws) |
| PASS | F6 one pool, spFixF + spFM: each species by its own setting | spFixF female 1775 / 1800, male 2378 / 2400; spFM female 624 / 600, male 2359 / 2400; largest \|z\| 1.87 |
| PASS | F7 `patch_sample_emigrants` | female 1825 / 1800, +0.93; male 2373 / 2400, −1.23 |
| PASS | F8 `dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool`, m 1 (two patches) | female 1790 / 1800, −0.37; male 2409 / 2400, +0.41 |
| PASS | F9 `use_disp_trait=False` (m 0.5, same seeds): spFixF and spFM lose the same individuals | 500 / 500, same positions |
| PASS | F10 ValueError for trait 'both', rate 1.5, rate −0.1; 'male' + rate 0 and None / None boot | as expected |

### Life cycles (`regression_head`, `cost_none`, `bookkeeping`, plan items 1–3)

The life cycle has 4 patches × 2 habitats of 8 × 8, dormancy 10, 2 species, mutation on, seed 1 and 16 steps. Each step runs:

- dead selection 0.2 (neutral);
- reproduction at 0.5;
- one of 4 dispersal variants, rotating, so every one of the 7 dispersal functions is called, with `is_remove` True for steps 0–3 and 8–11 and False for steps 4–7 and 12–15;
- germination and dormancy.

| status | check | digest (first 12) |
|---|---|---|
| PASS | regression_head, sexual: default path, HEAD == working tree | `2020a9596db6` / `2020a9596db6` |
| PASS | regression_head, asexual: default path, HEAD == working tree | `401b16f18c01` / `401b16f18c01` |
| PASS | cost_none, sexual: scheme 3, `disp_cost` None == 0 == (0, 0) | `adb64e801b49` ×3 |
| PASS | cost_none, sexual: scheme 2 ignores `disp_cost` (none == 0.3 == (0.1, 0.5)) | `c7db4d73ecf0` ×3 |
| PASS | cost_none, asexual: scheme 3, `disp_cost` None == 0 == (0, 0) | `85ca1541a65e` ×3 |
| PASS | cost_none, asexual: scheme 2 ignores `disp_cost` | `f1d7f746632f` ×3 |

Bookkeeping passed in all 22 life cycles: 11 configurations × sexual / asexual. The configurations are:

- scheme 1;
- scheme 2: plain, `sex_specific`, with c 0.3, and with (0.1, 0.5);
- scheme 3: plain, `sex_specific`, with c 0 and with (0, 0);
- scheme 4: c 0.3, and `sex_specific` with (0.1, 0.5).

Individuals placed at the end, sexual:

| scheme | placed |
|---|---|
| 1 | 582 |
| 2 | 569 |
| 2, `sex_specific` | 563 |
| 3 | 359 |
| 3, `sex_specific` | 409 |
| 4, c 0.3 | 333 |
| 4, (0.1, 0.5) | 299 |

Asexual: 592 placed in every configuration. NOTE: in the asexual runs `sex_specific` has no effect (digest equal with and without), because asexual founders are all female and mdisp is never read.

### Range expansion (`full_sim`, plan item 10)

10 patches in a line, each with 2 habitats of 12 × 12. Only patch1 is seeded. Sexual, birth rate 1.0, death 0.1, among m 0.05 (gaussian, sigma 1), within m 0.2, disp traits mutate. 150 steps, 5 seeds per scheme, means over the seeds:

| scheme | migrants / (m × offspring) | migrant v − pool v | step patch 10 is half full | mean fdisp t = 0 / 50 / 100 / 150 |
|---|---|---|---|---|
| 1 | 1.001 | +0.001 | 81.2 | 0.497 / 0.467 / 0.455 / 0.454 |
| 2 | 1.000 | +0.057 | 79.2 | 0.497 / 0.588 / 0.620 / 0.627 |
| 3 | 0.631 | +0.057 | 93.0 | 0.497 / 0.632 / 0.675 / 0.677 |
| 4, c 0.3 | 0.405 | +0.070 | 125.6 | 0.497 / 0.575 / 0.629 / 0.640 |

| status | check |
|---|---|
| PASS | scheme 1: \|ratio − 1\| < 0.03 and \|migrant v − pool v\| < 0.01 |
| PASS | scheme 2: \|ratio − 1\| < 0.03 and migrant v − pool v > 0.02 |
| PASS | schemes 3 and 4: migrant v − pool v > 0.02; scheme 4 ratio < scheme 3 ratio |

## Known behaviour (accepted this round, recorded only)

- **Destination-order bias, schemes 2–4 with `is_remove=True`.** Destinations served earlier get emigrants with higher v (table above). With `is_remove=False` each j draws from the full pool, and the bias is gone.
- **Leave-rate gap, schemes 3 / 4.** Each destination is a separate draw, so P(leave) = 1 − Π_j (1 − v_k m_ij) instead of v_k · m. The gap grows with v² and m². At m = 0.5 the stays are only ~2 % above the ideal; at m = 1 an individual with v ≥ 0.8 leaves with 0.748 instead of 0.867. Serving the destinations in random order would not change this gap.
- **`disp_cost` with scheme 2** (`is_rescale=True`) is ignored silently. No one dies on the way, and the digest is equal with and without it.
