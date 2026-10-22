# extension/global_habitat_network.py test report (NemoBackend, sexual)

Generated 2026-10-21 22:06 by `test/nemoback-extension/test_extension.py all --full` (mode `full`). Scope: `extension/global_habitat_network.py` (the registry and the two network dispersals migrated to NemoBackend) and `metacommunity.add_patch`, on the NemoBackend call chain metacommunity -> patch -> habitat -> indfactory -> individual; sexual only.

## Environment

- python `3.11.17`, pynemo `2.4.1`, numpy `2.4.6`, `macOS-15.7.3-arm64-arm-64bit`
- HEAD `fe5f2b0` (`MIGRATE metacommunity: NemoBackend (sexual) + test/nemoback-metacommunity`), branch `features/developing`

| file | md5 | state |
|---|---|---|
| `extension/global_habitat_network.py` | `a301a867ba6953f570bf84ac6800123c` | modified (working tree) |
| `metaibm/metacommunity.py` | `5ad90e7b68f037158cbe6a34077987eb` | modified (working tree) |
| `metaibm/patch.py` | `29e70117ed847b7b66eacadc77bf84cf` | committed |
| `metaibm/habitat.py` | `cc2c2001c092d149325b6d1f42e4d2a7` | committed |
| `metaibm/indfactory.py` | `eacaf91df4dbf3ab7dc72100323355be` | committed |
| `metaibm/individual.py` | `2914e40044bb48e0b10364227282d915` | committed |
| `docs-developer/extension-global-habitat-network.md` | `fbba7309632ae4367a5e7492a13b117a` | modified (working tree) |

- originals: extension = `git show 1356387:extension/global_habitat_network.py` (last commit that touched it), add_patch = `git show fe5f2b0:metaibm/metacommunity.py` (HEAD). Their bodies are compiled into a subclass of the CURRENT metacommunity (they run on the NemoBackend patch / habitat). The record attribute is read from `_init_global_habitat_network_fields`: `global_habitat_id_2_index_dir` (both `global_habitat_id_2_index_dir` and `global_habitat_id_2_gidx_dir` accepted).
- mode `full`: experiment-size simulations (2000 steps, net habitats 40x40, GRFE 10x10), more seeds and repetitions

## Step 0 -- test-script fixes and baseline

Four test-script edits (no source code touched):

| file | edit |
|---|---|
| `test/nemoback-extension/test_extension.py` | S3 (no random.sample left): strip # comments (tokenize) before scanning; the agreed else comment mentions random.sample |
| `test/nemoback-extension/test_extension.py` | marker_bookkeeping: world m1 shares the factories with the later worlds -> recycle m1 (teardown clear-ups; C.release: without teardown's minted == recycled assert, which cannot hold while m is alive) before check_world([m]) |
| `test/nemoback-metacommunity/test_metacommunity.py` | R1: skip methods with fn.__module__ == 'extension.global_habitat_network' when collecting is_remove defaults |
| `test/nemoback-patch/test_chain.py` | lines 279-280: removed == set() and changed == set() (metacommunity.py is migrated); the stale "next migration round" finding deleted |

| suite | script | result | time |
|---|---|---|---|
| test/nemo-indfactory | accept_genetics.py | PASS | 0.5 s |
| test/nemo-indfactory | accept_guards.py | PASS | 0.5 s |
| test/nemo-indfactory | accept_lifecycle.py | PASS | 0.6 s |
| test/nemo-indfactory | smoke.py | PASS | 0.5 s |
| test/nemoback-habitat (run_all.py) | smoke.py | PASS |  |
| test/nemoback-habitat (run_all.py) | test_bookkeeping.py | PASS | 5.9 s |
| test/nemoback-habitat (run_all.py) | test_nemo.py | PASS | 13.8 s |
| test/nemoback-habitat | **run_all.py total** | PASS | 20.4 s |
| test/nemoback-patch (run_all.py --no-baseline) | smoke.py | PASS |  |
| test/nemoback-patch (run_all.py --no-baseline) | test_bookkeeping.py | PASS | 4.1 s |
| test/nemoback-patch (run_all.py --no-baseline) | test_nemo.py | PASS | 13.0 s |
| test/nemoback-patch (run_all.py --no-baseline) | test_chain.py | PASS | 0.7 s |
| test/nemoback-patch | **run_all.py --no-baseline total** | PASS | 19.5 s |
| test/nemoback-metacommunity (test_metacommunity.py all) | smoke | PASS |  |
| test/nemoback-metacommunity (test_metacommunity.py all) | bookkeeping | PASS | 2.1 s |
| test/nemoback-metacommunity (test_metacommunity.py all) | nemo | PASS | 39.8 s |
| test/nemoback-metacommunity (test_metacommunity.py all) | chain | PASS | 34.2 s |
| test/nemoback-metacommunity | **test_metacommunity.py all total** | PASS | 125.2 s |

- every folder was backed up before and restored after its run (their runners rewrite TEST_REPORT.md + results/): restored byte-identical (md5 of every file + git status vs before) = nemoback-habitat `True`, nemoback-patch `True`, nemoback-metacommunity `True`; test/nemo-indfactory untouched = True. Fresh reports kept as `results/baseline_{habitat,patch,metacommunity}_TEST_REPORT.md`.

## Smoke: PASS (0.02 s)

| flow | step | return | invariants |
|---|---|---|---|
| marker, network is_remove=True (default) | mainland.meta_dead_selection | [Dead selection] in mainland: there are 21 individuals dead in selection; there are 107 individuals in the mainland; there are 21  | True |
| marker, network is_remove=True (default) | mainland.meta_mainland_sexual_birth_mutate_germinate | [Birth process] in mainland: there are 21 individuals germinating from local habitat; there are 128 individuals in the mainland; t | True |
| marker, network is_remove=True (default) | islands.meta_dead_selection | [Dead selection] in islands: there are 20 individuals dead in selection; there are 76 individuals in the islands; there are 20 emp | True |
| marker, network is_remove=True (default) | islands.meta_sex_reproduce_calculation_with_offspring_marker_pool | [Reproduction into offspring_marker_pool] in islands: there are 26 offspring_marker born into the offspring_marker_pool; there are | True |
| marker, network is_remove=True (default) | islands.meta_colonize_from_propagules_rains | [Colonization process]: there are 1 individuals colonizing the metacommunity from mainland; there are 77 individuals in the metaco | True |
| marker, network is_remove=True (default) | islands.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool | [Dispersal within process] in islands: there are 6 offspring marker disperse into habs_immigrant_marker_pool within patches; there | True |
| marker, network is_remove=True (default) | dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool | [Dispersal among patches in global habitat network] in islands: there are 11 offspring markers disperse into immigrant_marker_pool | True |
| marker, network is_remove=True (default) | islands.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool | [Germination & Birth process] islands: there are 17 individuals germinating_and_birth from local offspring_marker_pool and immigra | True |
| marker, network is_remove=True (default) | islands.meta_disturbance_process_in_patches | [Disturbance process] occurred in ['p2'] | True |
| marker, network is_remove=True (default) | islands.meta_clear_up_offspring_marker_and_immigrant_marker_pool | 9 | True |
| marker, network is_remove=False | mainland.meta_dead_selection | [Dead selection] in mainland: there are 19 individuals dead in selection; there are 108 individuals in the mainland; there are 20  | True |
| marker, network is_remove=False | mainland.meta_mainland_sexual_birth_mutate_germinate | [Birth process] in mainland: there are 20 individuals germinating from local habitat; there are 128 individuals in the mainland; t | True |
| marker, network is_remove=False | islands.meta_dead_selection | [Dead selection] in islands: there are 9 individuals dead in selection; there are 55 individuals in the islands; there are 41 empt | True |
| marker, network is_remove=False | islands.meta_sex_reproduce_calculation_with_offspring_marker_pool | [Reproduction into offspring_marker_pool] in islands: there are 19 offspring_marker born into the offspring_marker_pool; there are | True |
| marker, network is_remove=False | islands.meta_colonize_from_propagules_rains | [Colonization process]: there are 2 individuals colonizing the metacommunity from mainland; there are 57 individuals in the metaco | True |
| marker, network is_remove=False | islands.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool | [Dispersal within process] in islands: there are 6 offspring marker disperse into habs_immigrant_marker_pool within patches; there | True |
| marker, network is_remove=False | dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool | [Dispersal among patches in global habitat network] in islands: there are 7 offspring markers disperse into immigrant_marker_pool  | True |
| marker, network is_remove=False | islands.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool | [Germination & Birth process] islands: there are 12 individuals germinating_and_birth from local offspring_marker_pool and immigra | True |
| marker, network is_remove=False | islands.meta_disturbance_process_in_patches | [Disturbance process] occurred in [] | True |
| marker, network is_remove=False | islands.meta_clear_up_offspring_marker_and_immigrant_marker_pool | 14 | True |
| pool, network offspring (default) | mainland.meta_dead_selection | [Dead selection] in mainland: there are 19 individuals dead in selection; there are 107 individuals in the mainland; there are 21  | True |
| pool, network offspring (default) | mainland.meta_mainland_sexual_birth_mutate_germinate | [Birth process] in mainland: there are 21 individuals germinating from local habitat; there are 128 individuals in the mainland; t | True |
| pool, network offspring (default) | islands.meta_colonize_from_propagules_rains | [Colonization process]: there are 1 individuals colonizing the metacommunity from mainland; there are 70 individuals in the metaco | True |
| pool, network offspring (default) | islands.meta_dead_selection | [Dead selection] in islands: there are 12 individuals dead in selection; there are 58 individuals in the islands; there are 38 emp | True |
| pool, network offspring (default) | islands.meta_sex_reproduce_mutate_into_offspring_pool | islands: there are 18 individuals born into the offspring_pool; there are 18 individuals in the offspring_pool | True |
| pool, network offspring (default) | islands.meta_dispersal_within_patch_from_offspring_to_immigrant_pool | [Dispersal within process] islands: there are 5 individuals disperse into habs_immigrant_pool within patches; there are 5 individu | True |
| pool, network offspring (default) | dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool | [Dispersal among patches in global habitat network] in islands: there are 9 individuals disperse into immigrant_pool among patches | True |
| pool, network offspring (default) | islands.meta_local_germinate_from_offspring_and_immigrant_pool | [Germination process] in islands: there are 11 individuals germinating from local offspring_pool and immigrant_pool in the local h | True |
| pool, network offspring (default) | islands.meta_dormancy_process_from_offspring_pool_and_immigrant_pool | [Dormancy process] in islands: there are 0 dormancy survived in the dormancy pool; there are 0 dormancy eliminated from the dorman | True |
| pool, network offspring (default) | islands.meta_disturbance_process_in_habitat | [Disturbance process] occurred in [] | True |
| pool, network offspring (default) | islands.meta_clear_up_offspring_and_immigrant_pool | 0 | True |

## Tests

| part | covers | result | checks | time |
|---|---|---|---|---|
| smoke | one step of each flow (marker True / False, pool) with invariants after every call | PASS | 1 / 1 | 0.7 s |
| static | A: diff whitelist, signatures, pool_name tuple, guard / else, TYPE_CHECKING, old name, callers, tags | PASS | 29 / 29 | 2.6 s |
| registry | B: order, global_index record, distance matrix (brute force + v3 + one recompute), edges, rate-matrix properties x 5 kernels, underflow | PASS | 22 / 22 | 118.5 s |
| dispersal | C: offspring network dispersal -- guard, patch_num < 2, 200-seed equivalence, conservation, v3 statistics, capping, distance law, Nemo level, chains | PASS | 42 / 42 | 220.6 s |
| marker | D: marker network dispersal -- both modes vs v3 (200 seeds each), germinate-and-birth homes, ValueError risk, default effect on GRFE births, chains | PASS | 17 / 17 | 235.3 s |
| sim | E: full simulations (2 scenarios x 3 groups + patch-level comparator), reproducibility, whole-run equivalence vs v3, ecology, disturbance | PASS | 32 / 32 | 22295.4 s |
| perf | F: one call at H = 16 / 64 / 256 / 1024 vs v3, cost breakdown, zero-skip variant | PASS | 5 / 5 | 120.2 s |
| fuzz | G: random mixes of every meta operation + both network dispersals + guarded calls | PASS | 2 / 2 | 8.4 s |

- `all --full` wall clock 22586.9 s (parts run in parallel processes; perf ran alone)

- last `all --quick` run (2026-10-21 15:49, 591.6 s): baseline PASS (11/11), smoke PASS (1/1), static PASS (29/29), perf PASS (5/5), sim PASS (36/36), registry PASS (22/22), dispersal PASS (42/42), marker PASS (17/17), fuzz PASS (2/2)

## A static: PASS

| check | result | detail |
|---|---|---|
| A1 every changed line of extension / add_patch / docs is one of the agreed changes (or blank / comment / layout) | PASS | [] |
| A1 agreed change present: registry rename | PASS | 2 |
| A1 agreed change present: record stores habitat.global_index | PASS | 2 |
| A1 agreed change present: distance matrix recomputed inside the registry | PASS | 1 |
| A1 agreed change present: is_remove in the signature | PASS | 4 |
| A1 agreed change present: guard | PASS | 1 |
| A1 agreed change present: else branch | PASS | 1 |
| A1 agreed change present: sampling | PASS | 5 |
| A1 agreed change present: install name | PASS | 2 |
| A1 agreed change present: TYPE_CHECKING import | PASS | 2 |
| A1 working tree: modified == extension, metacommunity.py, docs (+ the two Step 0 test files); untracked == test/nemoback-extension/ only | PASS | [' M docs-developer/extension-global-habitat-network.md', ' M extension/global_habitat_network.py', ' M metaibm/metacommunity.py', ' M test/nemoback-metacommunity/test_metacommunity.py', ' M test/nemoback-patch/test_chain.py', '?? test/nemoback-extension/'] |
| A2 rk_from_offspring_pool_to_immigrant_pool: is_remove=True right after method, before **kwargs | PASS | (self: 'metacommunity', total_disp_among_rate, method='uniform', is_remove=True, **kwargs) |
| A2 ing_marker_pool_to_immigrant_marker_pool: is_remove=True right after method, before **kwargs | PASS | (self: 'metacommunity', total_disp_among_rate, method='uniform', is_remove=True, **kwargs) |
| A3 marker: pool_name=('offspring_marker_pool',) is a TUPLE (a bare string would iterate its characters) | PASS | h_i_object.sample_offspring_without_replacement(migrants_i_j_num, pool_name=('offspring_marker_pool',)) |
| A3 marker: no guard (is_remove=False stays valid for markers) | PASS |  |
| A3 offspring: sample_offspring_without_replacement(migrants_i_j_num) with the default pool (no pool_name) | PASS | h_i_object.sample_offspring_without_replacement(migrants_i_j_num) |
| A3 offspring: guard is the first statement (before patch_num < 2 and any RNG / state change) | PASS | if not is_remove: |
| A3 offspring: guard text == metacommunity.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool | PASS | 'is_remove=False (sampling with replacement: the original stays in its source pool, a copy disperses) waits for pynemo Individual::operator=' |
| A3 offspring: else branch is `pass  # unreachable (guard above); TODO(copy_from) ...` | PASS |  |
| A3 offspring: no random.sample left in code (# comments stripped first) | PASS |  |
| A4 `from metaibm import patch, metacommunity` sits under `if TYPE_CHECKING:` (module level) | PASS |  |
| A4 runtime: the names are not bound in the extension module (TYPE_CHECKING is False) | PASS |  |
| A4 fresh process, import metaibm first: imports work, nothing bound at runtime | PASS | rc=0 'OK True False False'  |
| A4 fresh process, import extension.global_habitat_network first: imports work, nothing bound at runtime | PASS | rc=0 'OK True False False'  |
| A4 annotations are strings (never evaluated) | PASS | {'global_habitat_network_registry': {'patch_object': 'patch'}, 'dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool': {'self': 'metacommunity'}, 'dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_t |
| A5 old name global_habitat_id_idx_registry: no residue outside this suite (repo-wide incl. experiments, examples, .ipynb, docs, test) | PASS | [] |
| A6 no live caller passes is_remove=False to the offspring version (it raises) | PASS | [] |
| A7 no NemoBackend tag on the registry / the two network dispersals (they touch no Nemo object) | PASS | {'global_habitat_network_registry': False, 'dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool': False, 'dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool': False} |
| A7 no function of the extension is tagged NemoBackend | PASS | [] |

Diff lines by label: TYPE_CHECKING import x2, comment only x4, registry rename x2, record stores habitat.global_index x2, distance matrix recomputed inside the registry x1, is_remove in the signature x4, layout: patch_num < 2 early return folded into one line x8, blank line x2, guard x1, sampling x5, else branch x1, install name x2, add_patch: registry call renamed x2, add_patch: distance-matrix call moved into the registry x2, docs heading renamed x2

<details><summary>every +/- line</summary>

| file | +/- | line | label |
|---|---|---|---|
| extension/global_habitat_network.py | + | from typing import TYPE_CHECKING | TYPE_CHECKING import |
| extension/global_habitat_network.py | + | if TYPE_CHECKING: from metaibm import patch, metacommunity | TYPE_CHECKING import |
| extension/global_habitat_network.py | - | self.global_habitat_id_ls = []                           # element: (patch_id, habitat_id) | comment only |
| extension/global_habitat_network.py | - | self.global_habitat_id_2_index_dir = {}                  # {(patch_id, habitat_id): global_habitat_index}; row/col index for global_habitat_distance_matrix | comment only |
| extension/global_habitat_network.py | + | self.global_habitat_id_ls = []                           # list of (patch_id, habitat_id) in add_patch x add_habitat order; row/col order of the network matrice | comment only |
| extension/global_habitat_network.py | + | self.global_habitat_id_2_index_dir = {}                  # {(patch_id, habitat_id): global_habitat_index}; cache of translating id to habitat_object.global_inde | comment only |
| extension/global_habitat_network.py | - | def global_habitat_id_idx_registry(self, patch_object): | registry rename |
| extension/global_habitat_network.py | + | def global_habitat_network_registry(self, patch_object:'patch'): | registry rename |
| extension/global_habitat_network.py | - | self.global_habitat_id_2_index_dir[(patch_id, h_id)] = len(self.global_habitat_id_ls) - 1 | record stores habitat.global_index |
| extension/global_habitat_network.py | + | self.global_habitat_id_2_index_dir[(patch_id, h_id)] = h_object.global_index | record stores habitat.global_index |
| extension/global_habitat_network.py | + | self.update_global_habitat_distance_matrix() | distance matrix recomputed inside the registry |
| extension/global_habitat_network.py | - | def dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool(self, total_disp_among_rate, method='uniform', **kwargs): | is_remove in the signature |
| extension/global_habitat_network.py | + | def dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool(self:'metacommunity', total_disp_among_rate, method='uniform', is_re | is_remove in the signature |
| extension/global_habitat_network.py | - | if self.patch_num < 2: | layout: patch_num < 2 early return folded into one line |
| extension/global_habitat_network.py | - | log_info = '[Dispersal among patches in global habitat network] in %s: patch_num < 2, there are 0 individuals disperse into immigrant_pool among habitats \n' %  | layout: patch_num < 2 early return folded into one line |
| extension/global_habitat_network.py | - | return log_info | layout: patch_num < 2 early return folded into one line |
| extension/global_habitat_network.py | + |  | blank line |
| extension/global_habitat_network.py | + | if not is_remove: raise NotImplementedError('is_remove=False (sampling with replacement: the original stays in its source pool, a copy disperses) waits for pyne | guard |
| extension/global_habitat_network.py | + | if self.patch_num < 2: return '[Dispersal among patches in global habitat network] in %s: patch_num < 2, there are 0 individuals disperse into immigrant_pool am | layout: patch_num < 2 early return folded into one line |
| extension/global_habitat_network.py | - | migrants_to_j_indi_object_ls += random.sample(h_i_object.offspring_pool, migrants_i_j_num) | sampling |
| extension/global_habitat_network.py | + | if is_remove: migrants_to_j_indi_object_ls += h_i_object.sample_offspring_without_replacement(migrants_i_j_num) | sampling |
| extension/global_habitat_network.py | + | else: pass  # unreachable (guard above); TODO(copy_from): migrants += copies of random.sample(offspring_pool, n), the originals stay in their source pools; need | else branch |
| extension/global_habitat_network.py | - | def dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool(self, total_disp_among_rate, method='uniform', **kwarg | is_remove in the signature |
| extension/global_habitat_network.py | + | def dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool(self:'metacommunity', total_disp_among_rate, method='u | is_remove in the signature |
| extension/global_habitat_network.py | - | if self.patch_num < 2: | layout: patch_num < 2 early return folded into one line |
| extension/global_habitat_network.py | - | log_info = '[Dispersal among patches in global habitat network] in %s: patch_num < 2, there are 0 offspring markers disperse into immigrant_marker_pool among ha | layout: patch_num < 2 early return folded into one line |
| extension/global_habitat_network.py | - | return log_info | layout: patch_num < 2 early return folded into one line |
| extension/global_habitat_network.py | + |  | blank line |
| extension/global_habitat_network.py | + | if self.patch_num < 2: return '[Dispersal among patches in global habitat network] in %s: patch_num < 2, there are 0 offspring markers disperse into immigrant_m | layout: patch_num < 2 early return folded into one line |
| extension/global_habitat_network.py | - | migrants_to_j_offspring_marker_ls += random.sample(h_i_object.offspring_marker_pool, migrants_i_j_num) | sampling |
| extension/global_habitat_network.py | + | if is_remove: migrants_to_j_offspring_marker_ls += h_i_object.sample_offspring_without_replacement(migrants_i_j_num, pool_name=('offspring_marker_pool',)) | sampling |
| extension/global_habitat_network.py | + | else: migrants_to_j_offspring_marker_ls += random.sample(h_i_object.offspring_marker_pool, migrants_i_j_num) | sampling |
| extension/global_habitat_network.py | - | cls.global_habitat_id_idx_registry = global_habitat_id_idx_registry | install name |
| extension/global_habitat_network.py | + | cls.global_habitat_network_registry = global_habitat_network_registry | install name |
| metaibm/metacommunity.py | - | self.global_habitat_id_idx_registry(patch_object)                # global_habitat index in a ls, index for habitat-network matrix | add_patch: registry call renamed |
| metaibm/metacommunity.py | - | self.update_global_habitat_distance_matrix()                     # calculate global habitats distance matrix | add_patch: distance-matrix call moved into the registry |
| metaibm/metacommunity.py | - | #self.incremental_update_global_habitat_distance_matrix() | add_patch: distance-matrix call moved into the registry |
| metaibm/metacommunity.py | + | self.global_habitat_network_registry(patch_object)               # register attributes for global-habitat-network | add_patch: registry call renamed |
| docs-developer/extension-global-habitat-network.md | - | ### `global_habitat_id_idx_registry` | docs heading renamed |
| docs-developer/extension-global-habitat-network.md | + | ### `global_habitat_network_registry` | docs heading renamed |

</details>

Callers of the two network dispersals (repo-wide incl. untracked, excluding this suite):

| file:line | function | kind | passes is_remove | text |
|---|---|---|---|---|
| docs-developer/extension-global-habitat-network.md:266 | offspring | docs / report text | False | ### `dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool` |
| docs-developer/extension-global-habitat-network.md:287 | marker | docs / report text | False | ### `dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool` |
| experiments/model-sloss-GRFE.py:277 | marker | live call | False | write_in_logger_info += meta_obj.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_ |
| experiments/model-sloss-GRFE.py:304 | marker | live call | False | write_in_logger_info += meta_obj.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_ |
| experiments/model-sloss-global-habitat-network.py:261 | marker | live call | False | write_in_logger_info += meta_obj.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_ |
| experiments/model-sloss-global-habitat-network.py:288 | marker | commented out | False | #write_in_logger_info += meta_obj.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool |
| extension/global_habitat_network.py:142 | offspring | definition / install | True | def dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool(self:'metacommunit |
| extension/global_habitat_network.py:170 | marker | definition / install | True | def dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool(self |
| extension/global_habitat_network.py:221 | offspring | definition / install | False | cls.dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool = dispersal_among_ |
| extension/global_habitat_network.py:222 | marker | definition / install | False | cls.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool = di |
| playgrounds/model-simulator-GRFE.py:81 | marker | simulator config (method name) | False | {'target': 'islands', 'method': 'dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_ |
| playgrounds/model-simulator-GRFE.py:193 | marker | simulator config (method name) | False | {'target': 'islands', 'method': 'dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_ |
| test/lecacy_v3.1.0-v3.3.0/model-sloss-global-habitat-network.py:261 | marker | legacy v3 test | False | write_in_logger_info += meta_obj.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_ |
| test/lecacy_v3.1.0-v3.3.0/model-sloss-global-habitat-network.py:288 | marker | commented out | False | #write_in_logger_info += meta_obj.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool |
| test/lecacy_v3.1.0-v3.3.0/test_global_habitat_network.py:341 | offspring | legacy v3 test | False | log_info = meta.dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool( |
| test/lecacy_v3.1.0-v3.3.0/test_global_habitat_network.py:367 | marker | legacy v3 test | False | log_info = meta.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_mark |
| test/lecacy_v3.1.0-v3.3.0/test_global_habitat_network.py:389 | offspring | legacy v3 test | False | log_info = meta.dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool( |
| test/lecacy_v3.3.1/model-sloss-GRFE.py:277 | marker | legacy v3 test | False | write_in_logger_info += meta_obj.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_ |
| test/lecacy_v3.3.1/model-sloss-GRFE.py:304 | marker | legacy v3 test | False | write_in_logger_info += meta_obj.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_ |
| test/metaibm-nemo-run/simulate.py:310 | marker | test suite | False | write_in_logger_info += mon('dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_i |
| test/nemoback-metacommunity/TEST_REPORT.md:325 | offspring | docs / report text | False | / known / extension/global_habitat_network.py: dispersal_among_patches_in_global_habitat_network_from_offsprin |
| test/nemoback-metacommunity/TEST_REPORT.md:326 | marker | docs / report text | True | / known / extension/global_habitat_network.py: dispersal_among_patches_in_global_habitat_network_from_offsprin |
| test/nemoback-metacommunity/common.py:262 | marker | test suite | False | (('islands.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_po |
| test/nemoback-metacommunity/test_metacommunity.py:915 | offspring | test suite | False | isl.dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool(0.6, method='expon |
| test/nemoback-metacommunity/test_metacommunity.py:925 | marker | test suite | False | isl2.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool(0.6 |
| test/nemoback-metacommunity/test_metacommunity.py:929 | offspring | test suite | False | if not twice: finding('info', 'extension/global_habitat_network.py: dispersal_among_patches_in_global_habitat_ |
| test/nemoback-metacommunity/test_metacommunity.py:930 | offspring | test suite | False | else: finding('known','extension/global_habitat_network.py: dispersal_among_patches_in_global_habitat_network_ |
| test/nemoback-metacommunity/test_metacommunity.py:932 | marker | test suite | True | finding('known', 'extension/global_habitat_network.py: dispersal_among_patches_in_global_habitat_network_from_ |
| test/nemoback-metacommunity/test_metacommunity.py:1110 | marker | test suite | False | 'network_marker': among('dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immig |

Old name `global_habitat_id_idx_registry`: test/nemoback-extension/TEST_REPORT.md:137:| A5 old name global_habitat_id_idx_registry: n (this suite (v3 reference name)); test/nemoback-extension/TEST_REPORT.md:154:| extension/global_habitat_network.py | - | def (this suite (v3 reference name)); test/nemoback-extension/TEST_REPORT.md:180:| extension/global_habitat_network.py | - | cls (this suite (v3 reference name)); test/nemoback-extension/TEST_REPORT.md:182:| metaibm/metacommunity.py | - | self.global_ha (this suite (v3 reference name)); test/nemoback-extension/TEST_REPORT.md:186:| docs-developer/extension-global-habitat-netwo (this suite (v3 reference name)); test/nemoback-extension/TEST_REPORT.md:225:Old name `global_habitat_id_idx_registry`: test (this suite (v3 reference name)); test/nemoback-extension/common.py:46:ORIGINAL_REGISTRY = 'global_habitat_id_idx_registry' (this suite (v3 reference name)); test/nemoback-extension/test_extension.py:146:          ('registry rename', r'^def global_ (this suite (v3 reference name)); test/nemoback-extension/test_extension.py:160:          ('install name', r'^cls\.global_ha (this suite (v3 reference name)); test/nemoback-extension/test_extension.py:162:    META: [('add_patch: registry call rename (this suite (v3 reference name)); test/nemoback-extension/test_extension.py:281:    check('A5 old name global_habitat_id_idx (this suite (v3 reference name))

- `typing.get_type_hints(global_habitat_network_registry)` at runtime: NameError: name 'patch' is not defined (string annotations under TYPE_CHECKING are for static tools only)

## B registry and matrices: PASS

| check | result | detail |
|---|---|---|
| B0 new name installed | PASS |  |
| B0 old name gone | PASS |  |
| B0 install list == original with the one rename | PASS | diff [] |
| B0 add_patch calls only the registry (matrix recomputed inside it) | PASS | ['self.global_habitat_network_registry(patch_object)               # register attributes for global-habitat-network'] |
| B0 record attribute name read from _init_global_habitat_network_fields: global_habitat_id_2_index_dir | PASS | global_habitat_id_2_index_dir |
| B1 fixed world (3 patches, 2/3/4 habitats, unsorted names, global_index incl. 0 and 65535) == original and brute force | PASS |  |
| B1 fixed world: rate matrix x 5 kernels x m in (0.0, 0.01, 0.1, 0.5, 0.9, 1.0): diag == 1 - m, same patch == 0, rows sum to 1, no NaN | PASS | [] |
| B2 1000 random worlds (2-5 patches, 1-4 habitats, unsorted names, random global_index): order, record, distance (v3 + brute force), one recompute per add_patch | PASS | [] |
| B2 1000 random worlds: rate matrix == v3 (md5) and diag / same patch / row sum properties (every 4th world: 6 rates x 5 kernels) | PASS | [] |
| B3 single patch == original | PASS |  |
| B3 patch without habitats (first) == original, contributes no row | PASS | [] |
| B3 patch without habitats (middle) == original, contributes no row | PASS | [] |
| B3 patch without habitats (last) == original, contributes no row | PASS | [] |
| B3 duplicate global_index: recorded as passed, no validation, no error (agreed: record only) | PASS | {('pB', 'h0'): 7, ('pB', 'h1'): 7, ('pA', 'h0'): 7} [] |
| B3 habitat added after add_patch is not registered and the matrices keep their size (same as v3) | PASS |  |
| B4 the record (global_habitat_id_2_index_dir/global_habitat_id_2_gidx_dir) is only touched in the extension (init + registry), never read as a matrix index | PASS | ['extension/global_habitat_network.py:11', 'extension/global_habitat_network.py:19'] |
| B6 exponential(rho=2) underflows to exactly 0.0 from d = 372.567 (expected 745.13 / rho = 372.57: exp(-rho d) underflows first) | PASS | 372.5666095510334 |
| B6 not reachable in the experiments (max habitat distance << underflow distance) | PASS | {'grfe: max habitat distance (patch_num 16)': 43.84, 'grfe: exponential rho=0.2 zero from': 3717.618905948377, 'net: max habitat distance (16 patches x 4, coordinates 0..7)': 9.9, 'net: exponential rho=2 zero from': 372.5666095510334} |
| B6 control (A at d=300, no underflow): no NaN, A emits and B <-> C exchange | PASS | {'distance_A_to_B': 300.0, 'nan_rows_rate_matrix': [], 'nan_cells_O_times_M': 0, 'cells': 16, 'warnings': [], 'A_sent': 24035, 'A_received': 0, 'B_C_exchanged': 72033, 'logs': ['[Dispersal among patches in global habitat network] in islands: there are N indivi |
| B6 far habitat (A at d=400): its rate-matrix row is NaN | PASS | [['pA', 'h0']] |
| B6 far habitat: O * M is a MATRIX PRODUCT, 0 * NaN = NaN spreads the NaN row into every cell outside A's own column (12 / 16 NaN; A's diagonal is reset to 1 - m) -> after matrix_around nobody emits: A 0, B <-> C 0 (3 dispersal modes x 400 seeds), the log just says 0 | PASS | {'distance_A_to_B': 400.0, 'nan_rows_rate_matrix': [['pA', 'h0']], 'nan_cells_O_times_M': 12, 'cells': 16, 'warnings': ['RuntimeWarning: invalid value encountered in divide'], 'A_sent': 0, 'A_received': 0, 'B_C_exchanged': 0, 'logs': ['[Dispersal among patches |
| B6 patch_num == 2 but the other patch has no habitat: all rows NaN, 0 markers disperse, no error (no warning in the log) | PASS | {'nan_rows': 2, 'log': '[Dispersal among patches in global habitat network] in islands: there are 0 offspring markers disperse into immigrant_marker_pool among patches in global-habitat-network; there are 0 offspring markers in the immigrant_marker_pools in th |

Kernel underflow -- smallest distance where `dispersal_kernel_strength` returns exactly 0.0 (bisection):

| kernel | params | zero from distance |
|---|---|---|
| exponential | {'rho': 0.2} | 3718 |
| exponential | {'rho': 1} | 745.1 |
| exponential | {'rho': 2} | 372.6 |
| exponential | {'rho': 5} | 149 |
| exponential | {'rho': 10} | 74.51 |
| gaussian | {'sigma': 0.5} | 19.3 |
| gaussian | {'sigma': 1.5} | 57.91 |
| gaussian | {'sigma': 5} | 193 |
| cauchy | {'gamma': 1.0} | 1.341e+154 |
| power_law | {'alpha': 2.0, 'r0': 1.0} | 6.362e+161 |
| power_law | {'alpha': 8.0, 'r0': 1.0} | 2.824e+40 |
| uniform | {} | never (< 1e300) |

- experiments: {'grfe: max habitat distance (patch_num 16)': 43.84, 'grfe: exponential rho=0.2 zero from': 3717.618905948377, 'net: max habitat distance (16 patches x 4, coordinates 0..7)': 9.9, 'net: exponential rho=2 zero from': 372.5666095510334}
- pipeline (A at d=400 from B/C, exponential rho=2): {'far (d=400)': {'distance_A_to_B': 400.0, 'nan_rows_rate_matrix': [['pA', 'h0']], 'nan_cells_O_times_M': 12, 'cells': 16, 'warnings': ['RuntimeWarning: invalid value encountered in divide'], 'A_sent': 0, 'A_received': 0, 'B_C_exchanged': 0, 'logs': ['[Dispersal among patches in global habitat network] in islands: there are N individuals disperse into immigrant_pool among patches in global-habitat-network; there are N individuals in the immigrant pools in the metacommunity', '[Dispersal among patches in global habitat network] in islands: there are N offspring markers disperse into immigrant_marker_pool among patches in global-habitat-network; there are N offspring markers in the immigrant_marker_pools in the metacommunity']}, 'control (d=300)': {'distance_A_to_B': 300.0, 'nan_rows_rate_matrix': [], 'nan_cells_O_times_M': 0, 'cells': 16, 'warnings': [], 'A_sent': 24035, 'A_received': 0, 'B_C_exchanged': 72033, 'logs': ['[Dispersal among patches in global habitat network] in islands: there are N individuals disperse into immigrant_pool among patches in global-habitat-network; there are N individuals in the immigrant pools in the metacommunity', '[Dispersal among patches in global habitat network] in islands: there are N offspring markers disperse into immigrant_marker_pool among patches in global-habitat-network; there are N offspring markers in the immigrant_marker_pools in the metacommunity']}}
- other patch without habitats: {'nan_rows': 2, 'log': '[Dispersal among patches in global habitat network] in islands: there are 0 offspring markers disperse into immigrant_marker_pool among patches in global-habitat-network; there are 0 offspring markers in the immigrant_marker_pools in the metacommunity'}

## C offspring network dispersal: PASS

| check | result | detail |
|---|---|---|
| C1 is_remove=False (keyword, patch_num=3): NotImplementedError, text as metacommunity, state / random / np.random / minted / recycled / currentID unchanged | PASS | 'is_remove=False (sampling with replacement: the original stays in its source pool, a copy disperses) waits for pynemo I |
| C1 is_remove=False (positional, patch_num=3): NotImplementedError, text as metacommunity, state / random / np.random / minted / recycled / currentID unchanged | PASS | 'is_remove=False (sampling with replacement: the original stays in its source pool, a copy disperses) waits for pynemo I |
| C1 is_remove=False (kernel kwargs, patch_num=3): NotImplementedError, text as metacommunity, state / random / np.random / minted / recycled / currentID unchanged | PASS | 'is_remove=False (sampling with replacement: the original stays in its source pool, a copy disperses) waits for pynemo I |
| C1 is_remove=False (keyword, patch_num=1): NotImplementedError, text as metacommunity, state / random / np.random / minted / recycled / currentID unchanged | PASS | 'is_remove=False (sampling with replacement: the original stays in its source pool, a copy disperses) waits for pynemo I |
| C1 is_remove=False (positional, patch_num=1): NotImplementedError, text as metacommunity, state / random / np.random / minted / recycled / currentID unchanged | PASS | 'is_remove=False (sampling with replacement: the original stays in its source pool, a copy disperses) waits for pynemo I |
| C1 is_remove=False (kernel kwargs, patch_num=1): NotImplementedError, text as metacommunity, state / random / np.random / minted / recycled / currentID unchanged | PASS | 'is_remove=False (sampling with replacement: the original stays in its source pool, a copy disperses) waits for pynemo I |
| C2 patch_num < 2 (default is_remove): returns the log at once (== v3 text), nothing changed | PASS | [Dispersal among patches in global habitat network] in isl1: patch_num < 2, there are 0 individuals  |
| C4 every new immigrant came from an offspring_pool of ANOTHER patch | PASS | moved 14 |
| C4 migrants left their source pool (one owner per individual) | PASS |  |
| C4 offspring -n, immigrant +n, n == logged counter | PASS | off 17->3 imm 12->26 log 14 |
| C8 Nemo level: the dispersal call itself mints nothing, recycles nothing, currentID unchanged | PASS | (228, 66, 229) |
| C8 migrants keep their birth habitat: getHome() unchanged == global_index of the SOURCE habitat (!= destination) | PASS |  |
| C8 germination after the network dispersal: 14 immigrants germinate, each in its destination habitat, the rest stay in that habitat's immigrant_pool | PASS | [Germination process] in isl3: there are 29 individuals germinating from local offspring_pool and immigrant_pool in the  |
| C8 germination moves, it does not mint / recycle | PASS |  |
| C8 world invariants after dispersal + germination (one owner, alive, minted - recycled == held, microsites) | PASS |  |
| C8 clean-up: every held individual recycled exactly once (162), none twice, minted == recycled | PASS |  |
| C3 reference (v3 + removal): 200 seeds, invariants hold, never two pools | PASS |  |
| C3 current: 200 seeds, invariants hold, never two pools | PASS | [] |
| C3 current == v3 + removal on every one of 200 seeds (state digest, random + np.random, log, who moved where in which order) | PASS | {'digest': 200, 'rng': 200, 'ret': 200, 'moved': 200} |
| C4 current: no same-patch migrant, counter == moved, no Nemo birth / recycle / currentID change (every seed) | PASS |  |
| C6 forced over-request in the equivalence worlds: capping exercised | PASS | {'natural_seeds': 134, 'natural_seeds_capped': 35, 'natural_capped_total': 72, 'natural_requested_total': 1063, 'forced_seeds': 66, 'forced_capped_total': 3227} |
| C5 m=0.1: expected immigrants per destination, current (move) == v3 (copy) (Welch, Bonferroni alpha 0.0001667; 5000 calls each) | PASS | {'mean_current': [0.773, 0.757, 2.173, 2.109, 0.759, 0.724], 'mean_v3': [0.78, 0.762, 2.188, 2.144, 0.764, 0.722], 'p_current_vs_v3': [0.6013, 0.7312, 0.4436, 0.0736, 0.7465, 0.9135], 'p_current_vs_current': [0.4891, 0.2456, 0.2226, 0.7492, 0.2203, 0.8172], 'p |
| C5 m=0.5: expected immigrants per destination, current (move) == v3 (copy) (Welch, Bonferroni alpha 0.0001667; 5000 calls each) | PASS | {'mean_current': [3.87, 3.74, 10.812, 10.637, 3.783, 3.661], 'mean_v3': [3.884, 3.763, 10.833, 10.66, 3.782, 3.647], 'p_current_vs_v3': [0.3667, 0.2405, 0.1919, 0.1635, 0.9794, 0.4684], 'p_current_vs_current': [0.9698, 0.7515, 0.1928, 0.628, 0.3277, 0.6239], ' |
| C5 m=0.9: current total <= v3 total (removal caps what v3 over-copied) | PASS | {'mean_current': [6.972, 6.742, 19.449, 19.134, 6.822, 6.552], 'mean_v3': [6.978, 6.747, 19.459, 19.173, 6.804, 6.564], 'p_current_vs_v3': [0.7357, 0.8149, 0.5434, 0.0106, 0.3002, 0.5339], 'p_current_vs_current': [0.7909, 0.9031, 0.8797, 0.9475, 0.4684, 0.4599 |
| C5 m=1.0: current total <= v3 total (removal caps what v3 over-copied) | PASS | {'mean_current': [7.743, 7.494, 21.616, 20.705, 7.507, 6.301], 'mean_v3': [7.754, 7.481, 21.629, 21.306, 7.551, 7.317], 'p_current_vs_v3': [0.4462, 0.3577, 0.4203, 0.0, 0.0055, 0.0], 'p_current_vs_current': [0.2826, 0.5197, 0.0734, 0.1395, 0.7248, 0.3859], 'pi |
| C5 sanity: the test detects nothing between two independent current samples (m=0.5) | PASS |  |
| C5 which source individual migrates is uniform (chi2, habitat 0, m=0.5) | PASS | 0.7549 |
| C6 truncation quantified (m x pool): requests beyond the pool are cut; realized <= pool always | PASS | [(1, 1.0, 0.2548), (3, 1.0, 0.3403), (10, 1.0, 0.3257), (50, 1.0, 0.3257)] |
| C6 shortfall falls on the LATER destinations of the j loop (order bias) | PASS | [(1, 0.1, [0.0, 0.0002, 0.0007, 0.0005, 0.0005, 0.0005, 0.001, 0.0008]), (1, 0.5, [0.0, 0.0042, 0.0058, 0.0105, 0.0143, 0.0187, 0.0215, 0.018]), (1, 0.9, [0.0, 0.0143, 0.022, 0.033, 0.0425, 0.0527, 0.0562, 0.0603])] |
| C7 uniform {}: migrants per distance ~ rate matrix (chi2 p=1, 99954 migrants) | PASS | {'observed': [16666, 16665, 16664, 16637, 16661, 16661], 'expected': [16666.7, 16666.7, 16666.7, 16666.7, 16666.7, 16666.7]} |
| C7 exponential {'rho': 2}: migrants per distance ~ rate matrix (chi2 p=0.635, 99967 migrants) | PASS | {'observed': [86461, 11701, 1582, 198, 23, 2], 'expected': [86467.0, 11702.0, 1583.7, 214.3, 29.0, 3.9]} |
| C7 gaussian {'sigma': 1.5}: migrants per distance ~ rate matrix (chi2 p=0.904, 99974 migrants) | PASS | {'observed': [58024, 29777, 9818, 2069, 260, 26], 'expected': [58026.5, 29791.8, 9807.2, 2070.0, 280.1, 24.3]} |
| C7 cauchy {'gamma': 1.0}: migrants per distance ~ rate matrix (chi2 p=0.998, 99949 migrants) | PASS | {'observed': [54083, 21635, 10840, 6336, 4141, 2914], 'expected': [54094.3, 21637.7, 10818.9, 6364.0, 4161.1, 2924.0]} |
| C7 power_law {'alpha': 2.0, 'r0': 1.0}: migrants per distance ~ rate matrix (chi2 p=1, 99998 migrants) | PASS | {'observed': [48860, 21713, 12232, 7805, 5409, 3979], 'expected': [48847.5, 21710.0, 12211.9, 7815.6, 5427.5, 3987.5]} |
| C7 exponential {'rho': 0.5}: migrants per distance ~ rate matrix (chi2 p=1, 100004 migrants) | PASS | {'observed': [41412, 25119, 15247, 9219, 5597, 3410], 'expected': [41408.5, 25115.6, 15233.4, 9239.5, 5604.0, 3399.0]} |
| C7 exponential {'rho': 1.0}: migrants per distance ~ rate matrix (chi2 p=0.989, 99981 migrants) | PASS | {'observed': [63368, 23337, 8568, 3136, 1140, 432], 'expected': [63369.1, 23312.2, 8576.1, 3155.0, 1160.6, 427.0]} |
| C7 exponential rho=2: slope of log(migrants) vs distance = -2.1184 ~ -rho | PASS | -2.1184 |
| C7 exponential rho=0.5: slope of log(migrants) vs distance = -0.4998 ~ -rho | PASS | -0.4998 |
| C7 exponential rho=1.0: slope of log(migrants) vs distance = -1.0001 ~ -rho | PASS | -1.0001 |
| C9 chain seed=1 x 500 steps (pool schedule, network offspring dispersal every step): invariants after every call, teardown minted == recycled | PASS | {"ok": true, "ops": 5500, "network_migrants": 898, "seconds": 7.7} |
| C9 chain seed=2 x 500 steps (pool schedule, network offspring dispersal every step): invariants after every call, teardown minted == recycled | PASS | {"ok": true, "ops": 5500, "network_migrants": 617, "seconds": 7.5} |
| C9 chain seed=3 x 500 steps (pool schedule, network offspring dispersal every step): invariants after every call, teardown minted == recycled | PASS | {"ok": true, "ops": 5500, "network_migrants": 1283, "seconds": 8.2} |

Expected immigrants per destination habitat, current (move) vs v3 (copy), 5000 calls each, offspring pools [12, 11, 14, 10, 14, 12], exponential rho=1:

| m | mean current | mean v3 | p (Welch) per destination | total current | total v3 | pick uniform p |
|---|---|---|---|---|---|---|
| 0.1 | [0.773, 0.757, 2.173, 2.109, 0.759, 0.724] | [0.78, 0.762, 2.188, 2.144, 0.764, 0.722] | [0.6013, 0.7312, 0.4436, 0.0736, 0.7465, 0.9135] | 7.294 | 7.361 | 0.213 |
| 0.5 | [3.87, 3.74, 10.812, 10.637, 3.783, 3.661] | [3.884, 3.763, 10.833, 10.66, 3.782, 3.647] | [0.3667, 0.2405, 0.1919, 0.1635, 0.9794, 0.4684] | 36.504 | 36.57 | 0.7549 |
| 0.9 | [6.972, 6.742, 19.449, 19.134, 6.822, 6.552] | [6.978, 6.747, 19.459, 19.173, 6.804, 6.564] | [0.7357, 0.8149, 0.5434, 0.0106, 0.3002, 0.5339] | 65.672 | 65.724 | 0.9994 |
| 1.0 | [7.743, 7.494, 21.616, 20.705, 7.507, 6.301] | [7.754, 7.481, 21.629, 21.306, 7.551, 7.317] | [0.4462, 0.3577, 0.4203, 0.0, 0.0055, 0.0] | 71.367 | 73.039 | 1.0 |

Capping (C6): one source habitat with `pool` offspring, 8 destinations at equal distance (uniform kernel), 6000 calls per cell. `P(excess)` = P(rounded request > pool); `shortfall by destination order` = mean requested - realized for j = 1..8 in loop order:

| pool | m | P(excess) | mean truncated | requested | realized | share truncated | shortfall by destination order |
|---|---|---|---|---|---|---|---|
| 1 | 0.01 | 0.0 | 0.0 | 0.011 | 0.011 | 0.0 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] |
| 1 | 0.1 | 0.004 | 0.0042 | 0.104 | 0.1 | 0.0402 | [0.0, 0.0002, 0.0007, 0.0005, 0.0005, 0.0005, 0.001, 0.0008] |
| 1 | 0.5 | 0.081 | 0.093 | 0.495 | 0.402 | 0.188 | [0.0, 0.0042, 0.0058, 0.0105, 0.0143, 0.0187, 0.0215, 0.018] |
| 1 | 0.9 | 0.2222 | 0.281 | 0.885 | 0.604 | 0.3175 | [0.0, 0.0143, 0.022, 0.033, 0.0425, 0.0527, 0.0562, 0.0603] |
| 1 | 1.0 | 0.2548 | 0.3327 | 0.977 | 0.645 | 0.3404 | [0.0, 0.0172, 0.0278, 0.0413, 0.0495, 0.061, 0.066, 0.0698] |
| 3 | 0.01 | 0.0 | 0.0 | 0.033 | 0.033 | 0.0 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] |
| 3 | 0.1 | 0.0002 | 0.0002 | 0.304 | 0.304 | 0.0005 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0002] |
| 3 | 0.5 | 0.0443 | 0.051 | 1.48 | 1.429 | 0.0345 | [0.0, 0.0, 0.0, 0.0003, 0.0042, 0.0077, 0.0155, 0.0233] |
| 3 | 0.9 | 0.2625 | 0.3802 | 2.682 | 2.301 | 0.1418 | [0.0, 0.0, 0.0, 0.0117, 0.0398, 0.0783, 0.104, 0.1463] |
| 3 | 1.0 | 0.3403 | 0.5227 | 2.986 | 2.463 | 0.175 | [0.0, 0.0, 0.0, 0.018, 0.0587, 0.1097, 0.1427, 0.1937] |
| 10 | 0.01 | 0.0 | 0.0 | 0.104 | 0.104 | 0.0 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] |
| 10 | 0.1 | 0.0 | 0.0 | 0.977 | 0.977 | 0.0 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] |
| 10 | 0.5 | 0.0 | 0.0 | 4.956 | 4.956 | 0.0 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] |
| 10 | 0.9 | 0.067 | 0.0778 | 8.977 | 8.899 | 0.0087 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0002, 0.007, 0.0707] |
| 10 | 1.0 | 0.3257 | 0.4683 | 9.982 | 9.514 | 0.0469 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0033, 0.0767, 0.3883] |
| 50 | 0.01 | 0.0 | 0.0 | 0.495 | 0.495 | 0.0 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] |
| 50 | 0.1 | 0.0 | 0.0 | 4.956 | 4.956 | 0.0 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] |
| 50 | 0.5 | 0.0 | 0.0 | 24.977 | 24.977 | 0.0 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] |
| 50 | 0.9 | 0.0 | 0.0 | 44.956 | 44.956 | 0.0 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] |
| 50 | 1.0 | 0.3257 | 0.4683 | 49.982 | 49.514 | 0.0094 | [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4683] |

Distance law (C7): source pool 200, m=0.5, destinations at d = 1..6, 1000 calls per kernel:

| kernel | params | observed | expected | chi2 / df | p | log-slope |
|---|---|---|---|---|---|---|
| uniform | {} | [16666, 16665, 16664, 16637, 16661, 16661] | [16666.7, 16666.7, 16666.7, 16666.7, 16666.7, 16666.7] | 0.04 / 5 | 1.0 |  |
| exponential | {'rho': 2} | [86461, 11701, 1582, 198, 23, 2] | [86467.0, 11702.0, 1583.7, 214.3, 29.0, 3.9] | 3.43 / 5 | 0.6347 | -2.1184 |
| gaussian | {'sigma': 1.5} | [58024, 29777, 9818, 2069, 260, 26] | [58026.5, 29791.8, 9807.2, 2070.0, 280.1, 24.3] | 1.58 / 5 | 0.9036 |  |
| cauchy | {'gamma': 1.0} | [54083, 21635, 10840, 6336, 4141, 2914] | [54094.3, 21637.7, 10818.9, 6364.0, 4161.1, 2924.0] | 0.27 / 5 | 0.9981 |  |
| power_law | {'alpha': 2.0, 'r0': 1.0} | [48860, 21713, 12232, 7805, 5409, 3979] | [48847.5, 21710.0, 12211.9, 7815.6, 5427.5, 3987.5] | 0.13 / 5 | 0.9997 |  |
| exponential | {'rho': 0.5} | [41412, 25119, 15247, 9219, 5597, 3410] | [41408.5, 25115.6, 15233.4, 9239.5, 5604.0, 3399.0] | 0.1 / 5 | 0.9998 | -0.4998 |
| exponential | {'rho': 1.0} | [63368, 23337, 8568, 3136, 1140, 432] | [63369.1, 23312.2, 8576.1, 3155.0, 1160.6, 427.0] | 0.57 / 5 | 0.9893 | -1.0001 |

- capping in the equivalence worlds: {'natural_seeds': 134, 'natural_seeds_capped': 35, 'natural_capped_total': 72, 'natural_requested_total': 1063, 'forced_seeds': 66, 'forced_capped_total': 3227}

## D marker network dispersal: PASS

| check | result | detail |
|---|---|---|
| D1 patch_num < 2, is_remove=True: same log as v3, nothing changed | PASS | [Dispersal among patches in global habitat network] in islands: patch_num < 2, there are 0 |
| D1 patch_num < 2, is_remove=False: same log as v3, nothing changed | PASS | [Dispersal among patches in global habitat network] in islands: patch_num < 2, there are 0 |
| D1 is_remove=True: markers move (offspring -n, immigrant +n, n == log, multiset of markers conserved, old immigrants kept in order) | PASS | n=3 new=3 |
| D1 is_remove=True: no Nemo birth / recycle / currentID change | PASS |  |
| D1 is_remove=False: markers copy (offspring_marker pools untouched, immigrant +n, n == log) | PASS | n=3 new=3 |
| D1 is_remove=False: no Nemo birth / recycle / currentID change | PASS |  |
| D2 germinate-and-birth after the network marker dispersal: each newborn sits in the habitat whose pool held its marker; consumed markers == newborns per habitat; getHome() == global_index of the PARENTS' (birth) habitat (354 newborns, 132 from another patch) | PASS | [] |
| D3 known risk reproduced: network marker dispersal -> disturbance of the birth patch -> germinate-and-birth raises ValueError part-way (18 / 40 seeds), invariants hold afterwards | PASS |  |
| D1 is_remove=True == v3 + removal on every one of 200 seeds (state, random + np.random, log, which markers reach which pool in which order) | PASS | {'digest': 200, 'rng': 200, 'ret': 200, 'new': 200} |
| D1 is_remove=False == v3 (copy) on every one of 200 seeds (state, random + np.random, log, which markers reach which pool in which order) | PASS | {'digest': 200, 'rng': 200, 'ret': 200, 'new': 200} |
| D1 is_remove=True: total markers conserved and counter == removed from the offspring pools (every seed) | PASS |  |
| D1 reference sanity: v3 copy removes nothing | PASS |  |
| D5 chain seed=1 x 300 steps, is_remove=true: invariants after every call, teardown minted == recycled | PASS | {"ok": true, "ops": 3000, "network_markers": 1046, "seconds": 4.5} |
| D5 chain seed=2 x 300 steps, is_remove=true: invariants after every call, teardown minted == recycled | PASS | {"ok": true, "ops": 3000, "network_markers": 233, "seconds": 3.8} |
| D5 chain seed=3 x 300 steps, is_remove=false: invariants after every call, teardown minted == recycled | PASS | {"ok": true, "ops": 3000, "network_markers": 3562, "seconds": 5.2} |
| D4 GRFE (full scale, 1024 habitats, filled with nearest-species founders): identical state before the network marker dispersal in both modes (8 seeds) | PASS |  |
| D4 births with is_remove=True <= births with False (fewer markers available) on every seed | PASS | 1.0 |

Default False -> True, GRFE step 0 from the same state (full scale): births True 117622 vs False 117623 (ratio 1.0)

| seed | births_true | births_false | network_markers | offspring_markers_after_true | offspring_markers_after_false | immigrant_markers_after_true | immigrant_markers_after_false | empty_sites |
|---|---|---|---|---|---|---|---|---|
| 1 | 14702 | 14702 | 346 | 35785 | 36131 | 4246 | 4246 | 14702 |
| 2 | 14680 | 14680 | 337 | 35786 | 36123 | 4388 | 4388 | 14680 |
| 3 | 14768 | 14769 | 352 | 35778 | 36130 | 4354 | 4354 | 14769 |
| 4 | 14730 | 14730 | 333 | 35804 | 36137 | 4288 | 4288 | 14730 |
| 5 | 14754 | 14754 | 370 | 35643 | 36013 | 4438 | 4438 | 14754 |
| 6 | 14667 | 14667 | 340 | 35740 | 36080 | 4449 | 4449 | 14669 |
| 7 | 14710 | 14710 | 360 | 35720 | 36080 | 4338 | 4338 | 14710 |
| 8 | 14611 | 14611 | 354 | 35714 | 36068 | 4488 | 4488 | 14611 |

## E full simulations: PASS

Scenarios copied from `experiments/model-sloss-global-habitat-network.py` (`net`: 16 patches x 4 habitats, rho 2, m_among 0.1) and `experiments/model-sloss-GRFE.py` (`grfe`: 16 patches x 64 habitats from the layout CSV, rho 0.2, m_among 0.01); sexual functions; 16 Nemo indfactories; distinct global_index per habitat (mainland 0.., metacommunity 100..). The metacommunity starts empty and is colonized from the mainland, as in the experiments. `net`'s sexual branch calls the patch-level marker dispersal in the experiment (network line commented out): here the network line is used. Groups: mT = marker + network is_remove default (True), mF = marker + is_remove=False, pool = individual pools + network offspring dispersal, mPatch = marker + the patch-level dispersal (the comparator). Scale full.

| check | result | detail |
|---|---|---|
| E1 repro_net_a: 500 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 repro_net_b: 500 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 repro_grfe_a: 30 steps, distinct global_index per habitat, invariants every 10 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 repro_grfe_b: 30 steps, distinct global_index per habitat, invariants every 10 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 equiv_net_cur: 1000 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 equiv_net_v3: 1000 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 equiv_grfe_cur: 50 steps, distinct global_index per habitat, invariants every 10 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 equiv_grfe_v3: 50 steps, distinct global_index per habitat, invariants every 10 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 disturb_mT: 300 steps, distinct global_index per habitat, invariants every 50 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 disturb_pool: 300 steps, distinct global_index per habitat, invariants every 50 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 net_mT_s1: 2000 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 grfe_mT_s1: 2000 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 net_mF_s1: 2000 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 grfe_mF_s1: 2000 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 net_pool_s1: 2000 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 grfe_pool_s1: 2000 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 net_mPatch_s1: 2000 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E1 grfe_mPatch_s1: 2000 steps, distinct global_index per habitat, invariants every 100 steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled | PASS |  |
| E2 reproducible: net mT seed 7, two fresh processes, digest of the whole state + RNG after each of 500 steps identical | PASS | first mismatch None |
| E2 reproducible: grfe mT seed 7, two fresh processes, digest of the whole state + RNG after each of 30 steps identical | PASS | first mismatch None |
| E3 whole run, net marker flow, is_remove=False: current == v3 extension bodies (registry + add_patch + marker dispersal) after every one of 1000 steps | PASS | first mismatch t=None |
| E3 whole run, grfe marker flow, is_remove=False: current == v3 extension bodies (registry + add_patch + marker dispersal) after every one of 50 steps | PASS | first mismatch t=None |
| E4 net: network marker dispersal vs the commented patch-level alternative -- alive / richness within 15% (2nd half means) | PASS | {'alive_mean_2nd_half': 0.001, 'richness_mean_2nd_half': 0.076, 'among_patch_migrants_mean_2nd_half': 0.002} |
| E5 net/mT: immigrants per habitat pair decay with distance (slope of log flow [-1.877], kernel rho 2) | PASS | [-1.8769459390949073] |
| E5 net/mF: immigrants per habitat pair decay with distance (slope of log flow [-1.897], kernel rho 2) | PASS | [-1.897405388163783] |
| E5 net/pool: immigrants per habitat pair decay with distance (slope of log flow [-1.921], kernel rho 2) | PASS | [-1.921106115053819] |
| E4 grfe: network marker dispersal vs the commented patch-level alternative -- alive / richness within 15% (2nd half means) | PASS | {'alive_mean_2nd_half': 0.0, 'richness_mean_2nd_half': 0.112, 'among_patch_migrants_mean_2nd_half': 0.003} |
| E5 grfe/mT: immigrants per habitat pair decay with distance (slope of log flow [-0.174], kernel rho 0.2) | PASS | [-0.17445544472718724] |
| E5 grfe/mF: immigrants per habitat pair decay with distance (slope of log flow [-0.167], kernel rho 0.2) | PASS | [-0.1674714917592595] |
| E5 grfe/pool: immigrants per habitat pair decay with distance (slope of log flow [-0.173], kernel rho 0.2) | PASS | [-0.17251609930467413] |
| E7 high disturbance (mT, patch_dist_rate 0.05): 218 events; each recycles exactly what the patch held (microsites + pools), leaves it empty, keeps dormancy, other patches untouched | PASS | {'events': 218, 'recycled_ok': 218, 'empty_after': 218, 'others_untouched': 218, 'dormancy_kept': 218, 'recycled_total': 225858, 'steps_to_50pct_recolonized': [0, 2, 0, 11, 15, 1, 5, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 3, 0, 0, 35, 0, 0, 0, 8, 0, 9, 0, 8, 48,  |
| E7 high disturbance (pool, patch_dist_rate 0.05): 250 events; each recycles exactly what the patch held (microsites + pools), leaves it empty, keeps dormancy, other patches untouched | PASS | {'events': 250, 'recycled_ok': 250, 'empty_after': 250, 'others_untouched': 250, 'dormancy_kept': 250, 'recycled_total': 369595, 'steps_to_50pct_recolonized': [0, 2, 0, 3, 0, 0, 0, 0, 2, 0, 0, 8, 2, 0, 0, 16, 0, 3, 7, 2, 5, 10, 2, 19, 0, 6, 1, 11, 19, 1, 9, 2, |

| run | status | H | sites | steps | total s | s/step | among-dispersal s (total) | RSS MB first/last | alive end | occupancy end | species end | teardown ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| repro_net_a | done | 64 | 6400 | 500 | 248.1 | 0.336 | 42.78 | [125.5, 174.8] | 6400 | 1.0 | 12 | True |
| repro_net_b | done | 64 | 6400 | 500 | 247.5 | 0.334 | 42.19 | [126.2, 179.7] | 6400 | 1.0 | 12 | True |
| repro_grfe_a | done | 1024 | 9216 | 30 | 446.7 | 14.184 | 401.11 | [131.6, 288.5] | 78 | 0.0085 | 12 | True |
| repro_grfe_b | done | 1024 | 9216 | 30 | 446.6 | 14.195 | 400.95 | [129.0, 286.4] | 78 | 0.0085 | 12 | True |
| equiv_net_cur | done | 64 | 6400 | 1000 | 392.0 | 0.258 | 52.2 | [125.5, 182.4] | 6400 | 1.0 | 13 | True |
| equiv_net_v3 | done | 64 | 6400 | 1000 | 392.9 | 0.259 | 52.08 | [128.0, 178.8] | 6400 | 1.0 | 13 | True |
| equiv_grfe_cur | done | 1024 | 9216 | 50 | 610.8 | 11.780 | 552.62 | [134.2, 309.9] | 71 | 0.0077 | 11 | True |
| equiv_grfe_v3 | done | 1024 | 9216 | 50 | 616.7 | 11.896 | 557.6 | [111.2, 318.7] | 71 | 0.0077 | 11 | True |
| disturb_mT | done | 64 | 102400 | 300 | 673.3 | 2.238 | 17.23 | [165.9, 519.4] | 45867 | 0.4479 | 10 | True |
| disturb_pool | done | 64 | 102400 | 300 | 354.7 | 1.176 | 20.98 | [165.6, 669.1] | 42860 | 0.4186 | 10 | True |
| net_mT_s1 | done | 64 | 102400 | 2000 | 15755.3 | 7.870 | 311.52 | [168.7, 880.7] | 102400 | 1.0 | 12 | True |
| grfe_mT_s1 | done | 1024 | 102400 | 2000 | 21891.0 | 10.930 | 17815.64 | [285.3, 902.0] | 102400 | 1.0 | 14 | True |
| net_mF_s1 | done | 64 | 102400 | 2000 | 15913.6 | 7.950 | 68.38 | [170.3, 792.8] | 102400 | 1.0 | 13 | True |
| grfe_mF_s1 | done | 1024 | 102400 | 2000 | 19288.2 | 9.628 | 14863.65 | [287.7, 909.3] | 102400 | 1.0 | 14 | True |
| net_pool_s1 | done | 64 | 102400 | 2000 | 5544.5 | 2.763 | 418.33 | [171.7, 990.4] | 102400 | 1.0 | 15 | True |
| grfe_pool_s1 | done | 1024 | 102400 | 2000 | 21584.2 | 10.777 | 17635.8 | [287.7, 1018.4] | 102400 | 1.0 | 13 | True |
| net_mPatch_s1 | done | 64 | 102400 | 2000 | 15920.3 | 7.953 | 76.28 | [168.9, 848.9] | 102400 | 1.0 | 11 | True |
| grfe_mPatch_s1 | done | 1024 | 102400 | 2000 | 5487.6 | 2.726 | 89.35 | [268.5, 728.5] | 102400 | 1.0 | 12 | True |

Ecology (2nd-half means over steps, averaged over seeds):

| scenario/group | runs | alive | occupancy end | species | births / step | among-patch migrants / step | decay slope | s/step |
|---|---|---|---|---|---|---|---|---|
| net/mT | 1 | 102400.0 | 1.0 | 12.85 | 12969.99 | 3911.777 | [-1.877] | 7.8704 |
| net/mF | 1 | 102400.0 | 1.0 | 11.34 | 13170.26 | 3914.103 | [-1.897] | 7.9496 |
| net/pool | 1 | 102400.0 | 1.0 | 13.26 | 12760.75 | 3916.817 | [-1.921] | 2.7634 |
| net/mPatch | 1 | 102341.0 | 1.0 | 11.94 | 13254.81 | 3903.254 | [-0.787] | 7.9527 |
| grfe/mT | 1 | 102399.9 | 1.0 | 13.75 | 13033.69 | 365.783 | [-0.174] | 10.9301 |
| grfe/mF | 1 | 102400.0 | 1.0 | 13.82 | 12961.1 | 363.429 | [-0.167] | 9.6275 |
| grfe/pool | 1 | 102400.0 | 1.0 | 12.59 | 12914.75 | 364.151 | [-0.173] | 10.7766 |
| grfe/mPatch | 1 | 102400.0 | 1.0 | 12.37 | 13063.09 | 364.728 | [-0.137] | 2.7257 |

- network vs patch-level (relative differences): {'net': {'alive_mean_2nd_half': 0.001, 'richness_mean_2nd_half': 0.076, 'among_patch_migrants_mean_2nd_half': 0.002}, 'grfe': {'alive_mean_2nd_half': 0.0, 'richness_mean_2nd_half': 0.112, 'among_patch_migrants_mean_2nd_half': 0.003}}
- is_remove True vs False (ratios mT / mF): {'net': {'births_ratio': 0.9848, 'alive_ratio': 1.0, 'migrants_ratio': 0.9994}, 'grfe': {'births_ratio': 1.0056, 'alive_ratio': 1.0, 'migrants_ratio': 1.0065}}

Distance decay (immigrants per habitat pair by distance bin):

- `disturb_mT`: slope -1.906; [[1, 465.52193], [2, 68.38393], [3, 9.42353], [4, 1.3913], [5, 0.20612], [6, 0.01339], [7, 0.00926]]
- `disturb_pool`: slope -1.96; [[1, 710.34211], [2, 101.45089], [3, 13.73824], [4, 1.97981], [5, 0.28059], [6, 0.03348], [7, 0.00617]]
- `net_mT_s1`: slope -1.877; [[1, 21055.78509], [2, 3336.66667], [3, 467.14706], [4, 72.33075], [5, 10.28324], [6, 1.51786], [7, 0.29938], [8, 0.04167]]
- `grfe_mT_s1`: slope -0.174; [[1, 6.3196], [2, 5.17587], [3, 4.34808], [4, 3.68316], [5, 2.98356], [6, 2.50373], [7, 2.08409], [8, 1.72761], [9, 1.42681], [10, 1.16407], [11, 0.97122], [12, 0.7896]]
- `net_mF_s1`: slope -1.897; [[1, 21444.0], [2, 3400.40774], [3, 476.76618], [4, 73.38975], [5, 10.39096], [6, 1.65179], [7, 0.25309]]
- `grfe_mF_s1`: slope -0.167; [[1, 6.36447], [2, 5.1531], [3, 4.2994], [4, 3.70402], [5, 3.00293], [6, 2.50119], [7, 2.08623], [8, 1.71951], [9, 1.42098], [10, 1.16671], [11, 0.97674], [12, 0.79752]]
- `net_pool_s1`: slope -1.921; [[1, 21336.87719], [2, 3382.4494], [3, 474.98382], [4, 71.81366], [5, 10.33644], [6, 1.60491], [7, 0.25309], [8, 0.02778]]
- `grfe_pool_s1`: slope -0.173; [[1, 6.38004], [2, 5.15116], [3, 4.34808], [4, 3.65988], [5, 2.97024], [6, 2.49857], [7, 2.06651], [8, 1.72965], [9, 1.40951], [10, 1.15493], [11, 0.97037], [12, 0.78934]]
- `net_mPatch_s1`: slope -0.787; [[1, 5792.67105], [2, 4775.2128], [3, 2924.45294], [4, 799.85093], [5, 427.84309], [6, 117.42411], [7, 86.51235], [8, 34.70833], [9, 16.2]]
- `grfe_mPatch_s1`: slope -0.137; [[1, 1.84341], [2, 1.88663], [3, 1.86154], [4, 1.86147], [5, 1.79503], [6, 1.7684], [7, 1.71191], [8, 1.66541], [9, 1.54207], [10, 1.37361], [11, 1.24979], [12, 1.07198]]

Disturbance (patch_dist_rate 0.05): {'mT': {'events': 218, 'recycled_ok': 218, 'empty_after': 218, 'others_untouched': 218, 'dormancy_kept': 218, 'recycled_total': 225858, 'steps_to_50pct_recolonized': [0, 2, 0, 11, 15, 1, 5, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 3, 0, 0, 35, 0, 0, 0, 8, 0, 9, 0, 8, 48, 0, 0, 36, 6, 10, 24, 3, 7, 0, 1, 0, 1, 0, 21, 0, 7, 0, 0, 11, 0, 6, 0, 7, 1, 2, 1, 1, 1, 2, 24, 8, 5, 3, 2, 2, 2, 1, 1, 3, 1, 2, 1, 4, 3, 1, 6, 1, 1, 1, 8, 10, 29, 41, 20, 7, 11, 1, 5, 1, 10, 4, 23, 1, 2, 5, 2, 15, 3, 46, 4, 3, 2, 7, 76, 23, 2, 80, 16, 2, 1, 2, 5, 2, 33, 12, 7, 1, 8, 2, 12, 19, 3, 2, 16, 7, 4, 12, 11, 0, 4, 2, 16, 4, 67, 6, 1, 2, 2, 1, 53, 32, 7, 43, 2, 1, 25, 5, 1, 24, 8, 13, 2, 2, 1, 1, 3, 2, 12, 8, 8, 1, 2, 11, 3, 7, 9, 1, 9, 11, 4, 1, 3, 3, 9, 38, 6, 2, 2, 6, 11, 3, 2, 1, 2, 14, 2, 5, 9, 24, 2, 25, 9, 26, 24, 20, 7, 5, 10, 1, 18, 11, 5, 12, 17, None, 11, None, 3, 2, None, None, None, None, 3, None, None], 'not_recolonized_by_end': 8}, 'pool': {'events': 250, 'recycled_ok': 250, 'empty_after': 250, 'others_untouched': 250, 'dormancy_kept': 250, 'recycled_total': 369595, 'steps_to_50pct_recolonized': [0, 2, 0, 3, 0, 0, 0, 0, 2, 0, 0, 8, 2, 0, 0, 16, 0, 3, 7, 2, 5, 10, 2, 19, 0, 6, 1, 11, 19, 1, 9, 2, 0, 0, 0, 0, 0, 0, 8, 2, 0, 0, 1, 1, 2, 0, 0, 4, 0, 1, 1, 3, 1, 3, 3, 15, 36, 13, 6, 81, 4, 2, 8, 0, 0, 1, 2, 2, 0, 2, 19, 1, 8, 14, 35, 5, 0, 12, 2, 16, 5, 2, 2, 0, 5, 16, 1, 17, 9, 1, 1, 1, 2, 4, 3, 2, 13, 1, 1, 1, 1, 1, 2, 3, 4, 14, 10, 3, 9, 3, 1, 1, 2, 2, 2, 2, 2, 16, 2, 2, 26, 15, 2, 17, 15, 7, 1, 18, 1, 3, 2, 9, 4, 3, 2, 1, 2, 17, 10, 2, 6, 23, 15, 50, 11, 14, 10, 15, 2, 3, 13, 18, 9, 5, None, 15, 2, 16, 5, 13, 18, 12, 20, 8, 1, 28, 16, 4, 15, 14, 46, 2, 44, 8, 25, 10, 6, 4, 2, 35, 37, 2, 2, 12, 27, 5, 23, 5, 41, 13, 2, None, None, 20, 3, 2, 16, 13, 18, 20, 12, 11, None, 3, 4, 1, 2, 40, 4, 9, 2, 11, 5, None, 5, 2, 2, 2, 3, 16, 4, 23, None, 1, 9, None, 4, 9, 16, 2, 13, None, None, 11, 10, 13, 1, 3, None, None, None, 1, None, 1, 2, None, 2, None, None, None], 'not_recolonized_by_end': 17}}

<details><summary>phase seconds per run</summary>

- `repro_net_a`: {'(digest)': 78.58, 'meta_dead_selection': 55.12, 'network_marker': 42.78, 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 39.65, 'pairwise_sexual_colonization_from_prpagules_rains': 11.61, 'mainland.meta_dead_selection': 10.44, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 3.49, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 2.65}
- `repro_net_b`: {'(digest)': 79.48, 'meta_dead_selection': 54.46, 'network_marker': 42.19, 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 39.17, 'pairwise_sexual_colonization_from_prpagules_rains': 11.85, 'mainland.meta_dead_selection': 10.69, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 3.36, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 2.78}
- `repro_grfe_a`: {'network_marker': 401.11, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 21.0, '(digest)': 1.75, 'pairwise_sexual_colonization_from_prpagules_rains': 1.3, 'meta_dead_selection': 0.99, 'mainland.meta_dead_selection': 0.45, 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 0.45, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 0.13}
- `repro_grfe_b`: {'network_marker': 400.95, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 21.46, '(digest)': 1.53, 'pairwise_sexual_colonization_from_prpagules_rains': 1.25, 'meta_dead_selection': 0.79, 'mainland.meta_dead_selection': 0.61, 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 0.46, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 0.21}
- `equiv_net_cur`: {'(digest)': 132.87, 'meta_dead_selection': 92.85, 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 68.31, 'network_marker': 52.2, 'mainland.meta_dead_selection': 16.34, 'pairwise_sexual_colonization_from_prpagules_rains': 14.11, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 5.04, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 4.77}
- `equiv_net_v3`: {'(digest)': 132.23, 'meta_dead_selection': 94.29, 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 68.08, 'network_marker': 52.08, 'mainland.meta_dead_selection': 16.11, 'pairwise_sexual_colonization_from_prpagules_rains': 14.11, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 5.54, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 4.61}
- `equiv_grfe_cur`: {'network_marker': 552.62, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 31.51, '(digest)': 2.49, 'pairwise_sexual_colonization_from_prpagules_rains': 2.0, 'meta_dead_selection': 1.09, 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 0.72, 'mainland.meta_dead_selection': 0.64, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 0.24}
- `equiv_grfe_v3`: {'network_marker': 557.6, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 32.26, '(digest)': 2.43, 'pairwise_sexual_colonization_from_prpagules_rains': 2.31, 'meta_dead_selection': 1.03, 'mainland.meta_dead_selection': 0.65, 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 0.57, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 0.22}
- `disturb_mT`: {'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 417.23, 'meta_dead_selection': 101.59, 'pairwise_sexual_colonization_from_prpagules_rains': 101.28, 'mainland.meta_dead_selection': 22.9, 'network_marker': 17.23, 'mainland.meta_mainland_sexual_birth_mutate_germinate': 4.22, 'meta_disturbance_process_in_patches': 3.61, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 1.82}
- `disturb_pool`: {'meta_dead_selection': 136.49, 'pairwise_sexual_colonization_from_prpagules_rains': 94.89, 'meta_sex_reproduce_mutate_into_offspring_pool': 35.2, 'meta_local_germinate_from_offspring_immigrant_and_dormancy_pool': 26.65, 'mainland.meta_dead_selection': 22.41, 'network_offspring': 20.98, 'meta_disturbance_process_in_patches': 6.25, 'mainland.meta_mainland_sexual_birth_mutate_germinate': 4.28}
- `net_mT_s1`: {'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 12206.59, 'meta_dead_selection': 2790.64, 'network_marker': 311.52, 'pairwise_sexual_colonization_from_prpagules_rains': 152.61, 'mainland.meta_dead_selection': 136.91, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 77.91, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 27.96, 'mainland.meta_mainland_sexual_birth_mutate_germinate': 27.49}
- `grfe_mT_s1`: {'network_marker': 17815.64, 'meta_dead_selection': 1666.8, 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 1048.15, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 958.21, 'pairwise_sexual_colonization_from_prpagules_rains': 149.27, 'mainland.meta_dead_selection': 112.22, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 79.89, 'mainland.meta_mainland_sexual_birth_mutate_germinate': 24.77}
- `net_mF_s1`: {'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 12558.93, 'meta_dead_selection': 2849.47, 'pairwise_sexual_colonization_from_prpagules_rains': 142.15, 'mainland.meta_dead_selection': 135.59, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 78.16, 'network_marker': 68.38, 'mainland.meta_mainland_sexual_birth_mutate_germinate': 28.38, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 28.21}
- `grfe_mF_s1`: {'network_marker': 14863.65, 'meta_dead_selection': 1816.07, 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 1147.22, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 1028.68, 'pairwise_sexual_colonization_from_prpagules_rains': 154.33, 'mainland.meta_dead_selection': 121.79, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 87.98, 'mainland.meta_mainland_sexual_birth_mutate_germinate': 29.2}
- `net_pool_s1`: {'meta_dead_selection': 3129.58, 'meta_sex_reproduce_mutate_into_offspring_pool': 1035.97, 'network_offspring': 418.33, 'pairwise_sexual_colonization_from_prpagules_rains': 328.5, 'meta_local_germinate_from_offspring_immigrant_and_dormancy_pool': 229.45, 'mainland.meta_dead_selection': 148.66, 'meta_dormancy_process_from_offspring_pool_and_immigrant_pool': 106.27, 'meta_dispersal_within_patch_from_offspring_to_immigrant_pool': 77.69}
- `grfe_pool_s1`: {'network_offspring': 17635.8, 'meta_dead_selection': 1917.32, 'meta_dispersal_within_patch_from_offspring_to_immigrant_pool': 949.63, 'meta_sex_reproduce_mutate_into_offspring_pool': 584.62, 'pairwise_sexual_colonization_from_prpagules_rains': 144.3, 'mainland.meta_dead_selection': 110.77, 'meta_local_germinate_from_offspring_immigrant_and_dormancy_pool': 98.11, 'meta_dormancy_process_from_offspring_pool_and_immigrant_pool': 70.93}
- `net_mPatch_s1`: {'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 12604.78, 'meta_dead_selection': 2799.55, 'pairwise_sexual_colonization_from_prpagules_rains': 149.35, 'mainland.meta_dead_selection': 133.65, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 76.65, 'patch_marker': 76.28, 'mainland.meta_mainland_sexual_birth_mutate_germinate': 28.83, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 27.45}
- `grfe_mPatch_s1`: {'meta_dead_selection': 2226.91, 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool': 1377.51, 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool': 1264.86, 'pairwise_sexual_colonization_from_prpagules_rains': 195.11, 'mainland.meta_dead_selection': 148.85, 'meta_sex_reproduce_calculation_with_offspring_marker_pool': 105.59, 'patch_marker': 89.35, 'mainland.meta_mainland_sexual_birth_mutate_germinate': 35.8}

</details>

## F performance: PASS

| check | result | detail |
|---|---|---|
| F H=16: zero-skip variant gives the identical result (state, RNG, log) -- the fix suggestion is behaviour-neutral | PASS | {'_pool_to_immigrant_marker_pool': True, 'fspring_pool_to_immigrant_pool': True} |
| F H=64: zero-skip variant gives the identical result (state, RNG, log) -- the fix suggestion is behaviour-neutral | PASS | {'_pool_to_immigrant_marker_pool': True, 'fspring_pool_to_immigrant_pool': True} |
| F H=256: zero-skip variant gives the identical result (state, RNG, log) -- the fix suggestion is behaviour-neutral | PASS | {'_pool_to_immigrant_marker_pool': True, 'fspring_pool_to_immigrant_pool': True} |
| F H=1024: zero-skip variant gives the identical result (state, RNG, log) -- the fix suggestion is behaviour-neutral | PASS | {'_pool_to_immigrant_marker_pool': True, 'fspring_pool_to_immigrant_pool': True} |
| F one call recorded at H = 16 / 64 / 256 / 1024 for current, v3, reference and the test-only zero-skip variant | PASS |  |

Median seconds per call (pools of 20 per habitat, m=0.1, exponential rho=2; perf ran alone):

| H (patches x habitats) | marker current is_remove=True | marker current is_remove=False | marker v3 (copy) | marker zero-skip variant (test-only) | offspring current | offspring v3 + removal | offspring zero-skip variant (test-only) |
|---|---|---|---|---|---|---|---|
| 16 (4x4) | 0.00112 | 0.00099 | 0.00098 | 0.00093 | 0.00112 | 0.00102 | 0.00094 |
| 64 (16x4) | 0.01869 | 0.01765 | 0.01862 | 0.01721 | 0.02197 | 0.02114 | 0.01939 |
| 256 (64x4) | 0.41232 | 0.36465 | 0.36855 | 0.32849 | 0.38887 | 0.35819 | 0.3106 |
| 1024 (16x64) | 4.41084 | 3.80742 | 3.87187 | 3.50292 | 4.40733 | 3.89913 | 3.38888 |

Breakdown of one current marker call (share of the total):

| H | total s | kernel strength matrix | normalize + rate (rest) | matrix_around | pool-size matrix | loop + sampling (rest) | sample calls | with n=0 |
|---|---|---|---|---|---|---|---|---|
| 16 | 0.0012 | 10% | 1% | 61% | 0% | 27% | 192 | 157 |
| 64 | 0.026 | 12% | 0% | 58% | 0% | 30% | 3840 | 3708 |
| 256 | 0.4026 | 10% | 0% | 62% | 0% | 29% | 64512 | 63993 |
| 1024 | 4.4249 | 11% | 0% | 59% | 0% | 29% | 983040 | 981003 |

Pool-size sweep (H=256, marker call, median s):

| pool per habitat | marker current is_remove=True | marker v3 (copy) | marker zero-skip variant (test-only) | current / v3 |
|---|---|---|---|---|
| 20 | 0.3635 | 0.3244 | 0.294 | 1.12 |
| 100 | 0.4477 | 0.3176 | 0.3103 | 1.41 |
| 400 | 0.7402 | 0.3261 | 0.3297 | 2.27 |

Build (all add_patch calls, distance matrix recomputed each time): [{'H': 16, 'patches': 4, 'current_s': 0.001, 'v3_s': 0.001}, {'H': 64, 'patches': 16, 'current_s': 0.018, 'v3_s': 0.018}, {'H': 256, 'patches': 64, 'current_s': 1.562, 'v3_s': 1.258}, {'H': 1024, 'patches': 16, 'current_s': 3.784, 'v3_s': 3.737}]

## G fuzz: PASS

| check | result | detail |
|---|---|---|
| G 9900 random operations (32 kinds, 3 processes x 30 runs x 110 ops): invariants after each, teardown minted == recycled after every run, no unexpected exception | PASS | [] |
| G guarded is_remove=False calls (offspring network): NotImplementedError + full snapshot unchanged | PASS | 892 |

- operations by kind: {'NETWORK markers is_remove=True': 884, 'NETWORK offspring is_remove=False (guard)': 892, 'patch-level among markers': 242, 'clear dormancy': 209, 'clear offspring': 232, 'dead selection': 246, 'network getters': 192, 'NETWORK markers is_remove=False': 899, 'clear immigrant': 265, 'within offspring+dormancy (sites)': 230, 'disturbance patches': 242, 'pairwise colonize': 232, 'clear offspring+immigrant': 215, 'dormancy offspring+immigrant': 224, 'NETWORK offspring (is_remove=True)': 917, 'patch-level among offspring->immigrant': 246, 'disturbance habitats': 198, 'germinate-and-birth from markers': 209, 'germinate offspring+dormancy': 216, 'offspring markers': 205, 'within markers': 224, 'dormancy offspring->dormancy': 243, 'within offspring+dormancy->immigrant': 239, 'germinate offspring+immigrant': 209, 'mainland dead + birth': 199, 'within offspring->immigrant': 212, 'patch-level among offspring+dormancy (sites)': 233, 'clear markers': 234, 'colonize': 233, 'germinate offspring+immigrant+dormancy': 227, 'patch-level among offspring+dormancy->immigrant': 223, 'reproduce into offspring_pool': 229}
- allowed pre-existing exceptions: {'germinate-and-birth: marker birth habitat without a parent pair (known risk D3)': 19, 'colonize: rain > mainland population': 21}

## Findings (19: bug 0, risk 3, info 10, known 2, todo(accepted) 4)

Severity: **bug** = wrong behaviour of the code under test; **risk** = correct today, fragile / surprising / costly; **info** = observation or intended change vs v3; **known** = agreed / pre-existing; **todo(accepted)** = style items the author decided to keep. Every suggested fix is a suggestion only -- **not applied**.

| severity | function | finding | evidence | suggested fix |
|---|---|---|---|---|
| risk | normalized_calculate_global_habitat_network_dispersal_kernel_strength_matrix -> both network dispersals | kernel underflow switches the WHOLE network dispersal off: when every other-patch kernel value of ONE habitat underflows to 0.0, its row of D / D.sum() is 0/0 = NaN; the dispersal number matrix is the matrix PRODUCT O * M (O diagonal), so 0 * NaN = NaN reaches every cell outside that habitat's own column; matrix_around turns NaN into 0 -> no habitat emits anything, silently (a numpy RuntimeWarning only, the log reports 0). Not reachable in the experiments (GRFE max habitat distance 43.8 vs underflow from 3718 at rho=0.2; net 9.9 vs 372.6 at rho=2); reachable with large coordinates, a large rho or a small sigma (exponential: d > 745.13 / rho; gaussian: d > 38.6 sigma; cauchy / power_law: no practical limit), and whenever a habitat has no habitat in another patch (e.g. patches without habitats, B6 last check). Pre-existing: the matrix code is unchanged from v3 (B2: rate matrices identical). | {"thresholds": [["exponential", {"rho": 0.2}, 3717.62], ["exponential", {"rho": 1}, 745.133], ["exponential", {"rho": 2}, 372.567], ["exponential", {"rho": 5}, 149.027], ["exponential", {"rho": 10}, 74.5133], ["gaussian", {"sigma": 0.5}, 19.302], ["gaussian", {"sigma": 1.5}, 57.906], ["gaussian", {"sigma": 5}, 193.02], ["cauchy", {"gamma": 1.0}, 1.34078e+154], ["power_law", {"alpha": 2.0, "r0": 1.0}, 6.36242e+161], ["power_law", {"alpha": 8.0, "r0": 1.0}, 2.82427e+40], ["uniform", {}, null]], "far": {"nan_rows_rate_matrix": [["pA", "h0"]], "nan_cells_O_times_M": 12, "cells": 16, "A_sent": 0, "B_C_exchanged": 0, "warnings": ["RuntimeWarning: invalid value encountered in divide"]}, "control":  | scale rows with np.multiply(pool_sizes[:, None], M) (element-wise, not O * M) and treat rows with D.sum() == 0 as "no dispersal" (log it); optionally compute the exponential kernel relative to the row minimum distance so it cannot underflow (not applied) |
| risk | meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool (+ network marker dispersal) | a marker whose birth habitat has no female-male pair any more raises ValueError (random.sample on an empty pair list) part-way through the loop -- earlier markers already born. Reachable through the network dispersal when a disturbance / clear-up of the birth patch runs between the network marker dispersal and germination; not in the experiments' order (disturbance after germination). Pre-existing (v3). | "network marker dispersal (m=0.9) -> p0.patch_disturbance_process() -> germinate-and-birth: ValueError in 18 / 40 seeds, births before the raise [0, 2, 4, 3, 2, 3, 7, 0]" | skip (or drop) a marker whose birth habitat has no parent pair; or clear the markers of a disturbed habitat in every immigrant_marker_pool (not applied) |
| risk | dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool, dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool | performance: the migrated loop calls habitat.sample_offspring_without_replacement for EVERY (i, j) pair of different patches, also when migrants_i_j_num == 0, and each call rebuilds the whole (pool_name, obj) list of the source pool: O(H^2 x pool) per call (v3: random.sample(pool, 0) is O(1)). H=1024, pools of 20: 4.41 s vs 3.87 s for v3 (x1.14); H=256, pools of 20 / 100 / 400: x1.12 / x1.41 / x2.27 vs v3. Skipping zero counts is bit-identical (random.sample(x, 0) draws nothing) and is faster than v3 (H=1024: 3.50 s). | {"H=16": {"marker current is_remove=True": 0.00112, "marker v3 (copy)": 0.00098, "marker zero-skip variant (test-only)": 0.00093, "offspring current": 0.00112, "offspring zero-skip variant (test-only)": 0.00094}, "H=64": {"marker current is_remove=True": 0.01869, "marker v3 (copy)": 0.01862, "marker zero-skip variant (test-only)": 0.01721, "offspring current": 0.02197, "offspring zero-skip variant (test-only)": 0.01939}, "H=256": {"marker current is_remove=True": 0.41232, "marker v3 (copy)": 0.36855, "marker zero-skip variant (test-only)": 0.32849, "offspring current": 0.38887, "offspring zero-skip variant (test-only)": 0.3106}, "H=1024": {"marker current is_remove=True": 4.41084, "marker v3 | in both network dispersals: `if migrants_i_j_num == 0: continue` before sampling (identical results), or sample once per source row (not applied) |
| info | test/nemoback-metacommunity (Step 0) | its own nested Step 0 reports a restore flag False (habitat `True`, patch `False`): its check requires `git status --short -- <suite>` to be EMPTY, and test/nemoback-patch/test_chain.py now carries the Step 0 edit (modified). The md5 manifests match (this suite's own restore check of the same folder is True). Not a defect. |  | compare git status before / after instead of requiring it to be empty (test-side only) (not applied) |
| info | dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool, dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool | diff lines outside the agreed list, non-functional: comment only x4, layout: patch_num < 2 early return folded into one line x8 (behaviour pinned by C2 / MB1 patch_num < 2 and the equivalence tests) | ["- self.global_habitat_id_ls = []                           # element: (patch_id, habitat_id)", "- self.global_habitat_id_2_index_dir = {}                  # {(patch_id, habitat_id): global_habitat_index}; row", "+ self.global_habitat_id_ls = []                           # list of (patch_id, habitat_id) in add_patch x add_h", "+ self.global_habitat_id_2_index_dir = {}                  # {(patch_id, habitat_id): global_habitat_index}; cac", "- if self.patch_num < 2:", "- log_info = '[Dispersal among patches in global habitat network] in %s: patch_num < 2, there are 0 individuals ", "- return log_info", "+ if self.patch_num < 2: return '[Dispersal among patches in global habitat network] in % |  |
| info | dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool, dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool | legacy v3 test folders still call the network dispersal functions (not part of any suite run here; the offspring calls now MOVE instead of copy) | ["test/lecacy_v3.1.0-v3.3.0/model-sloss-global-habitat-network.py:261", "test/lecacy_v3.1.0-v3.3.0/test_global_habitat_network.py:341", "test/lecacy_v3.1.0-v3.3.0/test_global_habitat_network.py:367", "test/lecacy_v3.1.0-v3.3.0/test_global_habitat_network.py:389", "test/lecacy_v3.3.1/model-sloss-GRFE.py:277", "test/lecacy_v3.3.1/model-sloss-GRFE.py:304"] |  |
| info | dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool | live callers in experiments / playgrounds use the marker version with the DEFAULT is_remove -> they now MOVE markers (v3 copied); quantified in D4 / E | ["experiments/model-sloss-GRFE.py:277 (live call)", "experiments/model-sloss-GRFE.py:304 (live call)", "experiments/model-sloss-global-habitat-network.py:261 (live call)", "playgrounds/model-simulator-GRFE.py:81 (simulator config (method name))", "playgrounds/model-simulator-GRFE.py:193 (simulator config (method name))"] |  |
| info | global_habitat_network_registry | duplicate global_index values are recorded without validation (agreed: the registry records, it does not assign or check); several habitats then share one Nemo homepatch | [[["pB", "h0"], 7], [["pB", "h1"], 7], [["pA", "h0"], 7]] |  |
| info | dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool | default is_remove False -> True: with the default the network marker dispersal now MOVES markers (v3 copied: the source pool kept every marker and a copy travelled). GRFE step 0 from a metacommunity filled with founders (full scale, 8 seeds): births 117622 (True) vs 117623 (False), ratio 1.0; network markers per call [346, 337, 352, 333, 370, 340, 360, 354]. Long-run effect: see E (mT vs mF). | [{"seed": 1, "births_true": 14702, "births_false": 14702, "network_markers": 346, "offspring_markers_after_true": 35785, "offspring_markers_after_false": 36131, "immigrant_markers_after_true": 4246, "immigrant_markers_after_false": 4246, "empty_sites": 14702}, {"seed": 2, "births_true": 14680, "births_false": 14680, "network_markers": 337, "offspring_markers_after_true": 35786, "offspring_markers_after_false": 36123, "immigrant_markers_after_true": 4388, "immigrant_markers_after_false": 4388, "empty_sites": 14680}, {"seed": 3, "births_true": 14768, "births_false": 14769, "network_markers": 352, "offspring_markers_after_true": 35778, "offspring_markers_after_false": 36130, "immigrant_markers_ |  |
| info | dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool vs dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool | net: same totals, different spatial pattern -- with the same rho the network kernel acts on habitat-to-habitat distances, the patch-level one on patch-centre distances: log-slope of immigrants per habitat pair vs distance [-1.88] (network) vs [-0.79] (patch-level); alive / richness / among-patch migrants differ by {'alive_mean_2nd_half': 0.001, 'richness_mean_2nd_half': 0.076, 'among_patch_migrants_mean_2nd_half': 0.002} |  |  |
| info | dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool vs dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool | grfe: same totals, different spatial pattern -- with the same rho the network kernel acts on habitat-to-habitat distances, the patch-level one on patch-centre distances: log-slope of immigrants per habitat pair vs distance [-0.17] (network) vs [-0.14] (patch-level); alive / richness / among-patch migrants differ by {'alive_mean_2nd_half': 0.0, 'richness_mean_2nd_half': 0.112, 'among_patch_migrants_mean_2nd_half': 0.003} |  |  |
| info | global habitat network (H=1024) | cost of one marker call at H=1024: kernel strength matrix 11%, normalize + rate (rest) 0%, matrix_around 59%, pool-size matrix 0%, loop + sampling (rest) 29%; the kernel matrix is a Python double loop of H^2 dispersal_kernel_strength calls and matrix_around draws H^2 single uniforms; both are recomputed on every call although the distance matrix is fixed. | {"kernel strength matrix": 0.4988, "normalize + rate (rest)": 0.0095, "matrix_around": 2.6288, "pool-size matrix": 0.0011, "loop + sampling (rest)": 1.2867} | cache the normalized kernel per (method, kwargs); vectorize the kernel (numpy on the distance matrix) and matrix_around (one uniform array) (not applied) |
| info | global_habitat_network_registry -> update_global_habitat_distance_matrix | every add_patch recomputes the full H x H distance matrix in a Python double loop (same as v3): building GRFE (16 patches, 1024 habitats) costs 3.8 s; incremental_update_global_habitat_distance_matrix exists but is unused. | [{"H": 16, "patches": 4, "current_s": 0.001, "v3_s": 0.001}, {"H": 64, "patches": 16, "current_s": 0.018, "v3_s": 0.018}, {"H": 256, "patches": 64, "current_s": 1.562, "v3_s": 1.258}, {"H": 1024, "patches": 16, "current_s": 3.784, "v3_s": 3.737}] | use the incremental update in the registry (or one numpy computation after the last add_patch) (not applied) |
| known | global_habitat_network_registry | a habitat added to a patch after add_patch is not registered (no row / column in the network, distance matrix not updated); unchanged from v3 | "pA.add_habitat(\"late\") after add_patch -> not in global_habitat_id_ls, distance matrix stays 9 x 9" | validation module: forbid add_habitat after add_patch, or re-register the patch (not applied) |
| known | dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool | matrix_around rounds each (i, j) separately: with high m the rounded row can exceed the offspring pool; sample_offspring_without_replacement caps with min() and the excess is cut from the destinations visited LAST in the j loop (v3 copied, never short). Agreed: not fixed, quantified. | {"P(excess) m=1.0 by pool 1/3/10/50": [0.2548, 0.3403, 0.3257, 0.3257], "P(excess) m=0.9": [0.2222, 0.2625, 0.067, 0.0], "P(excess) m=0.5": [0.081, 0.0443, 0.0, 0.0], "mean truncated m=1.0": [0.3327, 0.5227, 0.4683, 0.4683], "share of request truncated m=1.0": [0.3404, 0.175, 0.0469, 0.0094], "shortfall by destination order (pool 3, m=1.0)": [0.0, 0.0, 0.0, 0.018, 0.0587, 0.1097, 0.1427, 0.1937]} | (not requested) sample the whole row once (e.g. multinomial of min(pool, sum) over destinations, or shuffle the j order) (not applied) |
| todo(accepted) | global_habitat_network_registry | docstring kept as is (no three-line INPUT / RETURN format) -- accepted by the author | [] |  |
| todo(accepted) | dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool | docstring kept as is (no three-line INPUT / RETURN format) -- accepted by the author | ["dispersal from habitat_i offspring_pool to habitat_j immigrant_pool in the global habitat network"] |  |
| todo(accepted) | dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool | docstring kept as is (no three-line INPUT / RETURN format) -- accepted by the author | ["dispersal from habitat_i offspring_marker_pool to habitat_j immigrant_marker_pool in the global habitat network"] |  |
| todo(accepted) | global_habitat_network_registry, dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool, dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool | annotation spacing kept as is (no space after the colon) -- accepted by the author | ["patch_object:'patch'", "self:'metacommunity'"] |  |

## Method coverage

`vs v3` compares each function with `1356387:extension/global_habitat_network.py` (AST without docstrings / annotations, and source text). `tests` = parts of this suite that name the function.

| function | vs v3 | tag | tests |
|---|---|---|---|
| _init_global_habitat_network_fields | docstring / comment / annotation only |  | registry, common (schedules) |
| global_habitat_network_registry | changed (renamed) |  | static, registry, perf |
| update_global_habitat_distance_matrix | unchanged |  | static, perf, common (schedules) |
| incremental_update_global_habitat_distance_matrix | unchanged |  | static, perf |
| get_global_habitat_network_dormancy_pool_num_matrix | unchanged |  | fuzz |
| get_global_habitat_network_offspring_pool_num_matrix | unchanged |  | fuzz |
| get_global_habitat_network_offspring_marker_pool_num_matrix | unchanged |  | registry, perf, fuzz |
| get_global_habitat_network_empty_sites_num_matrix | unchanged |  | fuzz |
| calculate_global_habitat_network_dispersal_kernel_strength_matrix | unchanged |  | perf |
| normalized_calculate_global_habitat_network_dispersal_kernel_strength_matrix | unchanged |  | registry |
| global_habitat_dispersal_among_rate_matrix | unchanged |  | registry, dispersal, marker, perf |
| get_global_habitat_network_emigrants_matrix | unchanged |  | fuzz |
| get_global_habitat_network_immigrants_matrix | unchanged |  | fuzz |
| get_global_habitat_network_disp_among_num_matrix | unchanged |  | fuzz |
| dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool | changed |  | static, registry, dispersal, perf, fuzz, simulate_ext (E), common (schedules) |
| dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool | changed |  | smoke, static, registry, marker, perf, fuzz, sim, simulate_ext (E), common (schedules) |
| install_global_habitat_network_methods | changed |  | registry |

- `metacommunity.add_patch`: changed (one registry call; the distance matrix is recomputed inside the registry) -- static A1 / B0, every world build in B-G.

## Repository state

```
 M docs-developer/extension-global-habitat-network.md
 M extension/global_habitat_network.py
 M metaibm/metacommunity.py
 M test/nemoback-metacommunity/test_metacommunity.py
 M test/nemoback-patch/test_chain.py
?? test/nemoback-extension/
 docs-developer/extension-global-habitat-network.md |  2 +-
 extension/global_habitat_network.py                | 34 ++++++++++++----------
 metaibm/metacommunity.py                           |  4 +--
 test/nemoback-metacommunity/test_metacommunity.py  |  2 +-
 test/nemoback-patch/test_chain.py                  |  5 ++--
 5 files changed, 24 insertions(+), 23 deletions(-)
```
