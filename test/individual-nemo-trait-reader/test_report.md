# individual.get_trait_value / get_trait_alleles -- test report

Generated on 2026-11-03 at 20:49 by `test_trait_readers.py` and `test_ecology_readers.py` in this folder.

**Scope.** The two public readers of `metaibm/individual.py`:

- `get_trait_value(trait_key, **kwargs)` (lines 68-75). This is the regular interface for the value Nemo stores, for all 7 traits.
- `get_trait_alleles(trait_key, method, **kwargs)` (lines 77-116, with its five `_get_*` helpers).

Both are tested in the working tree on top of `9b626dc`. The tree also contains the uncommitted pynemo additions `TTrait.get_value()` and `TTrait.get_wolb_value()` in metaibm-nemo `swig/nemo.i`, built into conda env `metaibm-nemo-test6`.

**Side effects: none.**

- No file outside this folder was edited, and nothing was committed.
- `metaibm/` was imported read-only (`sys.dont_write_bytecode`, `-B`). The `metaibm/__pycache__` timestamps are the same before and after the run.
- Nemo output went to temp directories, which were removed afterwards.

## Environment

- HEAD is `9b626dc`, on branch `features/developing`. `metaibm/individual.py` has uncommitted changes:
  - `get_trait_value` added;
  - strict `**kwargs` in both readers;
  - the legacy-section header and two docstrings changed.
- Conda env `metaibm-nemo-test6`:
  - Python 3.11.18, numpy 2.4.6, `macOS-15.7.3-arm64`.
  - pynemo `_nemo.cpython-311-darwin.so` md5 `201e8c8968f2e8b2e323b4827e0a0112`, built 2026-11-03 19:32 from metaibm-nemo `640b9da` plus the uncommitted `swig/nemo.i`.
- `metaibm-nemo-test5` is used only for `regression_env`. Its pynemo is the same, but without `get_value` / `get_wolb_value`.
- Each check group runs in its own child process. Nemo's RAND, `Individual::currentID`, the trait prototypes and the static genetic map are process-global. Inside one process, every species keeps the same loci counts and chromosome layout.

| file | md5 |
|---|---|
| `metaibm/individual.py` (working tree, under test) | `abde21ac769cb6c139b6579acdaedc0a` |
| `test/individual-nemo-trait-reader/test_trait_readers.py` | `02613c9d3ee6e011a84aa2fda76f16c8` |
| `test/individual-nemo-trait-reader/test_ecology_readers.py` | `2355e6179d17f74b79abdb5dacf87a01` |

Run:

```
cd test/individual-nemo-trait-reader
PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test6/bin/python -B test_trait_readers.py     # ~16 s, exit 1 on any FAIL
PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test6/bin/python -B test_ecology_readers.py   # ~28 s
```

To run one group, append its name, e.g. `test_trait_readers.py delet quant_layout:mapped`.

## Oracles

All oracles live in the test code only.

- **raw**: direct pynemo calls for each sequence position:
  - `TTrait.get_allele(pos, strand)` and `get_allele_value(pos, strand)`;
  - quant positions come from Nemo's own `get_locus_seq_pos`;
  - the locus counts of ntrl, delet and dmi come from `nemo_params`, not from the reader's `get_locus_number()`.
- **stored**: the value at the `TTrait::getValue()` address, read with `ctypes` (`c_double`, or `c_bool` for wolb). This is read-only.
- **formula**: values recomputed from the allele tables, following the Nemo source:
  - **delet**:
    - model 1: `W = (1-s)^hmz (1-hs)^htz * scale`;
    - model 2: `W = 1 - s (hmz + h htz) * scale`;
    - see `ttdeletmutations_bitstring.cc:907-921`.
  - **dmi**:
    - `W = prod_pairs (1 + T[pair][4 cA + cB])`, where `c = 2 strand0 + strand1`;
    - table rows are recycled over the pairs, and a 9-column table is expanded to 16 columns;
    - see `ttbdmi.cc:113-154, 646-648`.
  - **disp**: the mean of the two alleles.
  - **quant dominance**: `G = sum(a0 + a1 + k|a1 - a0|)`, see `ttquanti.cc:2571`.
- **legacy**: `get_indi_phenotype_ls()`, `phenotype_set`, and `TTrait.get_genotypic_value(d)` (G).

## Results: 117 PASS, 0 FAIL, 10 NOTE

| script / group | covers | PASS | FAIL | NOTE |
|---|---|---|---|---|
| `test_trait_readers.py` | | **109** | **0** | **10** |
| - routing | Type, dtype and shape of all 8 valid `(trait_key, method)` pairs, and of `get_trait_value` for all 7 traits: float, float, float, float, float, bool, None. The 6 invalid pairs and 5 unknown key/method values raise `ValueError`, and so do 5 unknown `get_trait_value` keys (including `'disp'` and `''`). Both readers work on a species without quant. After `recycle()`, all 15 calls raise `RuntimeError`. | 21 | 0 | 2 |
| - kwargs | 34 calls cover every branch of both readers: stray `trait_dim`, misspelled `trait_dims`, `trait=`, positional `trait_dim`, and extra `x=`. Each raises `TypeError` from the function signature. An invalid pair with kwargs still raises `ValueError`. Omitting `trait_dim` means axis 0. The caller's kwargs dict is not changed by `pop`. | 5 | 0 | 0 |
| - quant_layout (7 cases, one process each) | Cases: continuous, diallelic, diallelic_asymmetrical (`{{a}}`), full pleiotropy (40 loci per axis), `quanti_loci_per_trait {{30,10}}`, three axes `{{14,13,13}}`, and mapped (`quanti_recombination_rate 0.05`). Checks: alleles equal raw on every axis, with shape `(loci_d, 2)`; the value equals `get_indi_phenotype_ls()[d]`, `phenotype_set[t_d]` and the stored `_phenotypes[d]`; sum equals G equals P (Ve 0); di-allelic value by state (`+a/-a`, or `0/a`); continuous states are all 0; strands are exact (column 1 maternal). The maternal-strand switch rate between adjacent loci is 0.505 for free recombination and 0.057 for the mapped case. | 44 | 0 | 1 |
| - quant_genetics | Dominance `k = 0.3`: value equals G equals `sum(a0+a1+k|a1-a0|)`, and differs from the additive sum (min gap 0.028). Ve 0.5 over 4000 founders: value equals stored P; `P - G` has mean `[-0.0085, 0.0076]` and variance `[0.510, 0.498]`, both within 4 se. | 2 | 0 | 0 |
| - ntrl | Byte version (`ntrl_all 4`) and bitstring version (`ntrl_all 2`): int `(12, 2)` equals raw, all states are seen, strands are exact, and `get_trait_value('ntrl') is None`. | 6 | 0 | 0 |
| - delet | Fitness models 1 and 2, each with scaling factor 1 and 0.9, over 140 individuals per case (founders with `delet_init_freq 0.3`, plus 3 generations). Value equals formula equals stored, with max gap 1.1e-16; hmz ranges over 0-12 and htz over 6-22. Strands are exact (0 violations as read, 2077 if swapped). | 6 | 0 | 0 |
| - dmi | 16-column table with 3 rows (recycled over 4 pairs) and a 9-column table: value equals formula equals stored for 260 individuals each, and all 16 genotype codes are seen at every pair. Strand check on opposite-homozygote parents: 0.831 as read, 0.003 swapped. | 3 | 0 | 0 |
| - disp | Founder init: fem 0.37, mal 0.81, on both genders. `set_value()` writes and `get_trait_value` reads: after `set_allele_value`, the value stays stale until `set_value()`. Strands: column 0 is maternal (reversed, as documented). Under mutation 0.5 with mean 0.3, the value equals the allele mean equals stored, and the alleles stay in [0, 1] (the 640b9da fix). | 6 | 0 | 0 |
| - wolb | The value is bool and equals the alleles reader and stored (154 of 320 infected). With transmission 1, a kid is infected exactly when its mother is; of 82 kids with an infected father and an uninfected mother, none is infected. A clone keeps its mother's state. `set_allele_value(0,0,1/0)` gives True/False. With transmission 0.5, 0.509 of 4000 kids are infected. | 5 | 0 | 0 |
| - stored | One 7-trait species with every mutation on and Ve 0.5, 240 individuals: `get_trait_value` equals stored for quant (2 axes), delet, dmi, fdisp, mdisp and wolb; ntrl is None; all 8 allele pairs equal raw. | 3 | 0 | 0 |
| - lifecycle | `make_copy`: everything equal. `make_clone` without mutation: equal to the mother. `make_clone` with mutation: stored values are current (equal to a fresh copy). 20 of 20 newborns built on recycled Nemo objects read their own data. Old wrappers raise `RuntimeError`. Minted equals recycled. | 6 | 0 | 0 |
| - out_of_range | `trait_dim` = `traits_num`, and `trait_dim = -1`, in both readers. Recorded only; see Findings. | 0 | 0 | 5 |
| - regression_head | `get_trait_alleles` (8 pairs) plus `get_indi_phenotype_ls`, `phenotype_set` and gender, 240 individuals, same seeds: HEAD `individual.py` and the working tree are bit-identical (md5 `b35737683cc4`). | 1 | 0 | 0 |
| - regression_env | The same dump under test5's pynemo and test6's pynemo is bit-identical. On test5, `get_trait_value('delet')` raises `AttributeError` (dependency, NOTE). | 1 | 0 | 1 |
| - perf | µs per call; see below. | 0 | 0 | 1 |
| `test_ecology_readers.py` | | **8** | **0** | **0** |
| - cycle_asexual / sexual / mix | example2's life cycle (dead selection, births, dispersal, colonization, germination) on a mainland of 4 × 100 sites plus 4 islands of 100 sites, 2 species with all 7 traits and every mutation on, 30 steps. At every step, for every living individual (400-790), both readers equal stored and the oracles: quant legacy, the delet/dmi formulas, the disp mean in [0, 1], and wolb. Minted equals recycled per species (about 2400-2800 each). | 6 | 0 | 0 |
| - recycled_refs | 1002 wrappers of individuals that died in their step raise `RuntimeError` from both readers, while their Nemo objects now hold newborns that read correctly. | 1 | 0 | 0 |
| - repro | Sexual life cycle, same seeds, two processes: the landscape digest (values plus allele tables) is identical (`1623ce1df727`). | 1 | 0 | 0 |

Per-call cost on test6:

| call | µs |
|---|---|
| `get_trait_value('quant', trait_dim=1)` | 0.70 |
| `get_trait_value('delet')` | 0.62 |
| `get_trait_value('wolb')` | 0.63 |
| `get_trait_value('ntrl')` | 0.63 |
| `get_indi_phenotype_ls()` (legacy) | 0.96 |
| `get_trait_alleles('quant', 'get_allele_value', ...)`, 20 loci | 17.2 |
| `get_trait_alleles('delet', ...)`, 30 loci | 18.7 |
| `get_trait_alleles('fdisp', ...)` | 1.30 |
| `get_trait_alleles('wolb', ...)` | 0.72 |

## Findings

No defect was found in either reader under the stated contract. Every value equals what Nemo stores, and every allele table equals the raw pynemo calls. Both agree with the formulas recomputed from the Nemo source. They also hold across mutation, cloning, copying and memory reuse through recycling, and through MetaIBM's full life cycle. The uncommitted changes leave `get_trait_alleles` and the legacy readers bit-identical to HEAD.

Behaviour recorded, not judged (NOTE):

1. **Out-of-range `trait_dim` fails silently.** This follows the agreed "no defensive raises" rule; Nemo's `assert`s are compiled out by `-DNDEBUG`. On a species with 2 quanti axes:
   - `get_trait_value('quant', trait_dim=2)` reads past `_phenotypes` and returns a different denormal on each run (e.g. `2.16e-314`, `5e-324`).
   - `get_trait_alleles('quant', 'get_allele_value', trait_dim=2)` returns a `(20, 2)` table equal to **axis 0**. `get_num_loci(2)` returns 20 and the positions read are 0-19.
   - `trait_dim = -1` raises `OverflowError` from SWIG (`argument 2 of type 'unsigned int'`) in both readers.

   If you ever want a guard, the right place is the routing line of each reader. That is your decision; nothing was changed.
2. **A trait absent from the species raises `KeyError`, not `ValueError`.** Examples are `'delet'` on a quant-only species, or `'quant'` on a species without quant. The error comes from `get_trait()` → `indfactory.traits_index`, before routing reaches the `else`.
3. **`diallelic_asymmetrical` with a scalar `quanti_allele_value`.** Nemo stores `+a/-a`, not `0/a`; only the matrix form `{{a}}` gives `0/a`. This is the "D4" combination of `test/nemoback-quanti-init-freq`. The readers return exactly what Nemo stores (sum equals G), so this is a parameter-form matter, not a reader issue.
4. **Dependency.** `get_trait_value` for delet, dmi, fdisp, mdisp and wolb needs a pynemo with `TTrait.get_value` / `get_wolb_value`, i.e. test6 or later. Under test5 it raises `AttributeError`. quant and ntrl work on both.
5. **wolb in the ecology runs.** Only founder *females* were infected (half of them). The infected fraction therefore rises from about 0.24 towards 0.5 within a few steps, because every kid takes its mother's state, and then drifts. Trajectory in the sexual run: 0.24 → 0.50 → 0.63 → 0.70. In the asexual run, all founders are female and the fraction stays near 0.5-0.6. This is the expected maternal transmission, not a reader effect.
