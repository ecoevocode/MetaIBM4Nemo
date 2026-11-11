# -*- coding: utf-8 -*-
''' test_disp_sampling.py -- disp-trait sampling of emigrants: hab_sample_emigrants / patch_sample_emigrants and their scheme methods

    Run (pynemo with TTrait.get_value -> conda env metaibm-nemo-test6):
        PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test6/bin/python -B test_disp_sampling.py [group ...]

    Schemes (use_disp_trait=True): scheme 2 = is_rescale=True, draw e_ij weighted by v_k; scheme 3 = is_rescale=False, leaves source_i (v_k),
    goes to sink_j (m_ij); scheme 4 = scheme 3 + dies on the way (disp_cost c or (c_fem, c_mal)). Scheme 1 = use_disp_trait=False (old MetaIBM).

    Groups:
        units   one call of hab_sample_emigrants / patch_sample_emigrants: P = v_k / sum(v), v_k * m_ij, arrive / dead / stay, copies, sex_specific, edges
        exact   through the real among function, one source -> 3 destinations in loop order; per-individual expectation, |z| < 4
        order   mean v of the emigrants per destination (destination-order bias, accepted this round) and the leave-rate gap at m = 1 (NOTE)
        fixed   metaibm_nemo_params dispersal_fixed_trait / dispersal_fixed_rate: the named sex leaves with the fixed rate; cost, scheme 2, mixed species, 3 entry points, validation

    Every group runs in its own child process (Nemo's RAND and Individual::currentID are process-global).
    A child prints one '@@ {json}' line per check: status PASS | FAIL | NOTE (NOTE = behaviour recorded, not judged).
    Nemo writes only into temporary directories; metaibm/ is imported read-only (no __pycache__).
    The harness below (Ctx, run_child, parent_main) is shared with test_disp_lifecycle.py.
'''
import sys
sys.dont_write_bytecode = True                  # never drop __pycache__ into metaibm/ or this folder
import os, json, time, shutil, hashlib, tempfile, subprocess, random
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[2]
GROUPS = ('units', 'exact', 'order', 'fixed')

NON_GENETIC = (('replicates', 1), ('generations', 1), ('patch_number', 1), ('patch_capacity', 1), ('aging', 1), ('run_mode', 'silent'))
MUTATION = (('disp_mutation_rate', 0.02), ('disp_mutation_mean', 0.05))
NO_MUTATION = (('disp_mutation_rate', 0.0), ('disp_mutation_mean', 0.05))
DISP_INIT = {
    'spU':  (('disp_init_distribution', 'uniform'), ('disp_init_dist_params', '{{0, 0}}')),   # every allele ~ U(0, 1)
    'spFM': (('disp_init_rate_fem', 0.2), ('disp_init_rate_mal', 0.8)),                        # fdisp 0.2, mdisp 0.8
    'sp00': (('disp_init_rate', 0.0),),                                                        # v = 0
}

##****************************** child side: shared harness ******************************
class Ctx:
    ''' one child process: seeds, imports from root, a temporary Nemo cache, check() / data() reporting, species and landscape builders '''
    def __init__(self, group, root, seed):
        sys.path.insert(0, str(root))
        random.seed(seed)
        import numpy as np
        np.random.seed(seed)
        from pynemo import _nemo
        from pynemo.params import SimParameters
        import metaibm
        from metaibm.indfactory import indfactory
        _nemo.cvar.SILENT_RUN = True
        self.np, self.SimParameters, self.indfactory, self.metaibm = np, SimParameters, indfactory, metaibm
        self.group = group
        self.cache = tempfile.mkdtemp(prefix='disp_%s_' % group.split(':')[0])

    def check(self, check_id, ok, detail='', note=False):
        status = 'NOTE' if note else ('PASS' if ok else 'FAIL')
        print('@@ ' + json.dumps({'group': self.group, 'id': check_id, 'status': status, 'detail': str(detail)}), flush=True)

    def data(self, key, value):
        print('@@DATA ' + json.dumps({'key': key, 'value': value}), flush=True)

    def species(self, name, init_key=None, seed=1, mutation=True, metaibm_nemo_params=None):
        ''' an indfactory with fdisp + mdisp (init_key in DISP_INIT), or without any trait (init_key None); metaibm_nemo_params goes to indfactory as is '''
        params = self.SimParameters()
        rows = NON_GENETIC + (('root_dir', self.cache), ('logfile', self.cache + '/nemo.log'), ('filename', name))
        if init_key is not None: rows += (MUTATION if mutation else NO_MUTATION) + DISP_INIT[init_key]
        for key, value in rows: params.add(key, value)
        return self.indfactory(name, (), params, trait_keys=('fdisp', 'mdisp') if init_key is not None else (), seed=seed, metaibm_nemo_params=metaibm_nemo_params)

    def landscape(self, locations, habs_per_patch, length, width, dormancy=0):
        ''' patches at the given locations, each with habs_per_patch habitats of length x width '''
        meta_obj = self.metaibm.metacommunity(metacommunity_name='meta')
        global_index = 0
        for patch_index, location in enumerate(locations):
            patch_obj = self.metaibm.patch(patch_name='patch%d' % (patch_index + 1), patch_index=patch_index, location=location)
            for hab_index in range(habs_per_patch):
                patch_obj.add_habitat(hab_name='h%d' % (hab_index + 1), hab_index=hab_index, hab_location=(location[0] + 0.1 * hab_index, location[1]), num_env_types=1,
                                      env_types_name=('env1',), mean_env_ls=[0.5], var_env_ls=[0.0], length=length, width=width, dormancy_pool_max_size=dormancy, global_index=global_index)
                global_index += 1
            meta_obj.add_patch(patch_name='patch%d' % (patch_index + 1), patch_object=patch_obj)
        return meta_obj

    def fill(self, meta_obj, factories, patch_names=None, sexual=True):
        ''' every microsite of the named patches gets a founder of a random species (random gender if sexual, else female) '''
        np = self.np
        names = list(factories)
        for patch_name, patch_obj in meta_obj.set.items():
            if patch_names is not None and patch_name not in patch_names: continue
            for h_object in patch_obj.set.values():
                for row in range(h_object.length):
                    for col in range(h_object.width):
                        factory = factories[names[np.random.randint(len(names))]]
                        gender = np.random.choice(('male', 'female')) if sexual else 'female'
                        h_object.add_individual(factory.make_founder(gender, home_habitat=h_object.global_index), row, col)

    def close(self):
        shutil.rmtree(self.cache, ignore_errors=True)

def z_check(observed, expected, variance):
    ''' RETURN: (|z| < 4, 'obs .. exp .. z ..') '''
    z_value = (observed - expected) / variance ** 0.5 if variance > 0 else 0.0
    return abs(z_value) < 4, 'obs %d exp %.1f z %+.2f' % (observed, expected, z_value)

##****************************** groups ******************************
def group_units(c):
    ''' direct calls of the per-pair API on one habitat / one patch '''
    np = c.np
    factories = {key: c.species(key, key, seed=1, mutation=False) for key in ('spU', 'spFM', 'sp00')}
    meta = c.landscape([(0, 0), (0, 1)], habs_per_patch=2, length=4, width=4, dormancy=100)
    h1, h2, patch1 = meta.set['patch1'].set['h1'], meta.set['patch1'].set['h2'], meta.set['patch1']
    reps = 30000
    def within_4se(simulated, expected, standard_error): return bool(np.all(np.abs(simulated - expected) < 4 * standard_error + 1e-12))

    # U1 scheme 2, e_ij = 1
    pool_ls = [factories['spU'].make_founder('female') for _ in range(10)]
    v_array = np.array([indi.get_trait_value('fdisp') for indi in pool_ls]); hit_array = np.zeros(10); left_ok = True
    for _ in range(reps):
        h1.offspring_pool[:] = pool_ls
        emigrant_ls = h1.hab_sample_emigrants(1)
        hit_array[pool_ls.index(emigrant_ls[0])] += 1
        left_ok &= len(h1.offspring_pool) == 9 and emigrant_ls[0] not in h1.offspring_pool
    expected = v_array / v_array.sum()
    c.check('U1 scheme 2, e_ij=1: P = v_k / sum(v); the emigrant leaves the pool', left_ok and within_4se(hit_array / reps, expected, np.sqrt(expected * (1 - expected) / reps)),
            'sim %s / exp %s' % (np.round(hit_array / reps, 4).tolist(), np.round(expected, 4).tolist()))

    # U2 scheme 3, one call
    rate = 0.4; hit_array = np.zeros(10)
    for _ in range(reps):
        h1.offspring_pool[:] = pool_ls
        for indi in h1.hab_sample_emigrants(rate, is_rescale=False): hit_array[pool_ls.index(indi)] += 1
    expected = v_array * rate
    c.check('U2 scheme 3: P(k -> j) = v_k * m_ij', within_4se(hit_array / reps, expected, np.sqrt(expected * (1 - expected) / reps)),
            'sim %s / exp %s' % (np.round(hit_array / reps, 4).tolist(), np.round(expected, 4).tolist()))

    # U3 scheme 4 with per-sex cost
    pool_fm = [factories['spU'].make_founder(g) for g in ['female'] * 3000 + ['male'] * 3000]
    v_fm = np.array([indi.get_trait_value('fdisp') for indi in pool_fm]); rate = 0.5
    h1.offspring_pool[:] = pool_fm
    recycled_before = factories['spU'].recycled_num
    emigrant_ls = h1.hab_sample_emigrants(rate, is_rescale=False, disp_cost=(0.2, 0.6))
    dead_num = factories['spU'].recycled_num - recycled_before
    female_arrive, male_arrive = sum(i.gender == 'female' for i in emigrant_ls), sum(i.gender == 'male' for i in emigrant_ls)
    leave_female, leave_male = (v_fm[:3000] * rate).sum(), (v_fm[3000:] * rate).sum()
    ok = abs(female_arrive - leave_female * 0.8) < 4 * np.sqrt(leave_female) and abs(male_arrive - leave_male * 0.4) < 4 * np.sqrt(leave_male) \
         and abs(dead_num - (leave_female * 0.2 + leave_male * 0.6)) < 4 * np.sqrt(leave_female + leave_male) and len(h1.offspring_pool) == 6000 - len(emigrant_ls) - dead_num
    c.check('U3 scheme 4, (c_fem 0.2, c_mal 0.6): arrive / dead (recycled) / pool size', ok,
            'female arrive %d (exp %.0f), male arrive %d (exp %.0f), dead %d (exp %.0f), pool left %d' % (female_arrive, leave_female * 0.8, male_arrive, leave_male * 0.4, dead_num, leave_female * 0.2 + leave_male * 0.6, len(h1.offspring_pool)))
    for indi in emigrant_ls + h1.offspring_pool: indi.recycle()
    h1.offspring_pool[:] = []

    # U4 is_remove=False
    pool_big = [factories['spU'].make_founder('female') for _ in range(4000)]
    v_big = np.array([indi.get_trait_value('fdisp') for indi in pool_big])
    h1.offspring_pool[:] = pool_big
    minted_before, recycled_before = factories['spU'].minted_num, factories['spU'].recycled_num
    copy_ls = h1.hab_sample_emigrants(0.5, is_rescale=False, disp_cost=0.5, is_remove=False)
    expected_copies = (v_big * 0.5 * 0.5).sum()
    ok = h1.offspring_pool == pool_big and not set(map(id, copy_ls)) & set(map(id, pool_big)) and factories['spU'].minted_num - minted_before == len(copy_ls) \
         and factories['spU'].recycled_num == recycled_before and abs(len(copy_ls) - expected_copies) < 4 * np.sqrt(expected_copies)
    copy_ls2 = h1.hab_sample_emigrants(5, is_remove=False)
    ok = ok and len(copy_ls2) == 5 and h1.offspring_pool == pool_big and len({indi.individual_id for indi in copy_ls2}) == 5
    c.check('U4 is_remove=False: pool untouched; copies only for survivors (m 0.5, c 0.5); nothing recycled', ok,
            'scheme 4 copies %d (exp %.0f), recycled %d; scheme 2 copies %d' % (len(copy_ls), expected_copies, factories['spU'].recycled_num - recycled_before, len(copy_ls2)))
    for indi in copy_ls + copy_ls2 + pool_big: indi.recycle()
    h1.offspring_pool[:] = []

    # U5 sex_specific: fdisp 0.2 / mdisp 0.8
    pool_sex = [factories['spFM'].make_founder(g) for g in ['female'] * 50 + ['male'] * 50]
    for sex_specific, expected_share in ((True, 0.8), (False, 0.5)):
        male_num = 0
        for _ in range(20000):
            h1.offspring_pool[:] = pool_sex
            male_num += h1.hab_sample_emigrants(1, sex_specific=sex_specific)[0].gender == 'male'
        c.check('U5 sex_specific=%s: share of males ~ %.1f' % (sex_specific, expected_share), abs(male_num / 20000 - expected_share) < 4 * np.sqrt(expected_share * (1 - expected_share) / 20000), 'sim %.4f' % (male_num / 20000))

    # U6 edges
    pool_zero = [factories['sp00'].make_founder('female') for _ in range(5)] + [factories['spU'].make_founder('female') for _ in range(5)]
    h1.offspring_pool[:] = pool_zero
    emigrant_ls = h1.hab_sample_emigrants(8)
    zero_picked = 0
    for _ in range(3000):
        h1.offspring_pool[:] = pool_zero
        zero_picked += sum(indi.species_id == 'sp00' for indi in h1.hab_sample_emigrants(1.0, is_rescale=False))
    h1.offspring_pool[:] = []
    ok = len(emigrant_ls) == 5 and all(indi.species_id == 'spU' for indi in emigrant_ls) and zero_picked == 0 and h1.hab_sample_emigrants(3) == [] and h1.hab_sample_emigrants(0.5, is_rescale=False) == []
    h1.offspring_pool[:] = pool_zero
    ok = ok and h1.hab_sample_emigrants(0) == [] and len(h1.offspring_pool) == 10
    c.check('U6 v=0 never drawn; e=8 > 5 with v>0 -> 5; empty pool; e=0', ok, 'drawn %d (all v>0: %s); v=0 drawn by scheme 3: %d' % (len(emigrant_ls), all(i.species_id == 'spU' for i in emigrant_ls), zero_picked))
    h1.offspring_pool[:] = []

    # U7 patch version: all habitats, offspring + dormancy
    pool_a = [factories['spU'].make_founder('female') for _ in range(20)]; pool_b = [factories['spU'].make_founder('female') for _ in range(20)]
    h1.offspring_pool[:] = pool_a[:10]; h1.dormancy_pool[:] = pool_a[10:]; h2.offspring_pool[:] = pool_b[:10]; h2.dormancy_pool[:] = pool_b[10:]
    emigrant_ls = patch1.patch_sample_emigrants(25, pool_name=('offspring_pool', 'dormancy_pool'))
    left_ls = h1.offspring_pool + h1.dormancy_pool + h2.offspring_pool + h2.dormancy_pool
    ok = len(emigrant_ls) == 25 and len(left_ls) == 15 and not set(map(id, emigrant_ls)) & set(map(id, left_ls)) and set(map(id, emigrant_ls + left_ls)) == set(map(id, pool_a + pool_b))
    c.check('U7 patch_sample_emigrants: draws over all habitats and both pools, removes from each', ok, 'drawn %d, left %d' % (len(emigrant_ls), len(left_ls)))
    for pool in (h1.offspring_pool, h1.dormancy_pool, h2.offspring_pool, h2.dormancy_pool): pool[:] = []

    # U8 species without disp traits
    no_disp = c.species('spN', None, seed=1)
    pool_mix = [factories['spU'].make_founder('female') for _ in range(5)] + [no_disp.make_founder('female')]
    h1.offspring_pool[:] = pool_mix
    try:
        h1.hab_sample_emigrants(2); raised = False
    except KeyError:
        raised = True
    c.check('U8 species without fdisp/mdisp -> KeyError, pool untouched', raised and h1.offspring_pool == pool_mix, 'raised KeyError: %s' % raised)
    h1.offspring_pool[:] = []

def among_setup(c, seed):
    ''' 4 patches in a line, 1 habitat each; the source is patch1, destinations patch2 -> patch4 in loop order '''
    factory = c.species('spU', 'spU', seed=seed, mutation=False)
    meta = c.landscape([(0, k) for k in range(4)], habs_per_patch=1, length=25, width=25)
    return factory, meta, meta.set['patch1'].set['h1'], [meta.set['patch%d' % (j + 1)].set['h1'] for j in range(1, 4)]

KERNEL = dict(method='gaussian', sigma=1.0)

def group_exact(c):
    ''' plan items 5 / 6: per-individual expectations through dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool '''
    np = c.np
    factory, meta, source_hab, dest_habs = among_setup(c, 5)
    m_among = 0.5
    rate_row = np.asarray(meta.get_disp_among_rate_matrix(m_among, **KERNEL))[0, 1:4]
    reps, pool_size = 100, 400
    c.check('exact/m_1j', True, 'm_1j (patch2, patch3, patch4) = %s' % np.round(rate_row, 5).tolist(), note=True)

    def run_case(label, cost_female, cost_male, is_remove, disp_cost):
        v_ls, gender_ls, outcome_ls, copy_sets_ls, contract_ok = [], [], [], [], True
        for _ in range(reps):
            pool_ls = [factory.make_founder('female' if k % 2 == 0 else 'male') for k in range(pool_size)]
            source_hab.offspring_pool[:] = list(pool_ls)
            v_array = np.array([indi.get_trait_value('fdisp') for indi in pool_ls])       # read before: the dead are recycled
            gender_ls += [indi.gender for indi in pool_ls]
            recycled_before = factory.recycled_num
            meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(m_among, is_remove=is_remove, use_disp_trait=True, is_rescale=False, disp_cost=disp_cost, **KERNEL)
            if is_remove:
                dest_of_id = {id(indi): j for j, h_object in enumerate(dest_habs) for indi in h_object.immigrant_pool}
                stay_id_set = {id(indi) for indi in source_hab.offspring_pool}
                outcome_ls += [str(dest_of_id[id(indi)]) if id(indi) in dest_of_id else ('stay' if id(indi) in stay_id_set else 'dead') for indi in pool_ls]
                contract_ok &= factory.recycled_num - recycled_before == sum(1 for indi in pool_ls if id(indi) not in dest_of_id and id(indi) not in stay_id_set)   # the dead are recycled
            else:
                contract_ok &= source_hab.offspring_pool == pool_ls and factory.recycled_num == recycled_before                                                 # pool untouched, nothing recycled
                index_of_v = {v: k for k, v in enumerate(v_array)}
                contract_ok &= len(index_of_v) == pool_size                                                                                                        # copies are mapped back by their disp value
                copy_sets = [set() for _ in range(pool_size)]
                for j, h_object in enumerate(dest_habs):
                    for copy in h_object.immigrant_pool: copy_sets[index_of_v[copy.get_trait_value('fdisp')]].add(j)
                copy_sets_ls += copy_sets
            v_ls += v_array.tolist()
            meta.meta_clear_up_offspring_and_immigrant_pool()
        v_all, gender_array = np.array(v_ls), np.array(gender_ls)
        cost_array = np.where(gender_array == 'male', cost_male, cost_female)
        leave_prob = v_all[:, None] * rate_row[None, :]                                                   # a_kj = v_k * m_1j
        lines, ok = [], contract_ok
        if is_remove:
            still_here = np.cumprod(np.hstack([np.ones((len(v_all), 1)), 1 - leave_prob[:, :-1]]), axis=1)   # prod over earlier destinations of (1 - a_kj')
            reach_prob = leave_prob * still_here
            prob_by_outcome = {str(j): reach_prob[:, j] * (1 - cost_array) for j in range(3)}
            prob_by_outcome['dead'] = (reach_prob * cost_array[:, None]).sum(axis=1)
            prob_by_outcome['stay'] = np.prod(1 - leave_prob, axis=1)
            outcome_array = np.array(outcome_ls)
            for outcome, prob in prob_by_outcome.items():
                for sex in ('female', 'male'):
                    sex_mask = gender_array == sex
                    line_ok, text = z_check(np.sum((outcome_array == outcome) & sex_mask), prob[sex_mask].sum(), (prob[sex_mask] * (1 - prob[sex_mask])).sum())
                    ok &= line_ok; lines.append('%s %s: %s' % ('patch%d' % (int(outcome) + 2) if outcome.isdigit() else outcome, sex, text))
        else:
            for j in range(3):
                prob = leave_prob[:, j] * (1 - cost_array)
                line_ok, text = z_check(sum(j in s for s in copy_sets_ls), prob.sum(), (prob * (1 - prob)).sum())
                ok &= line_ok; lines.append('copies at patch%d: %s' % (j + 2, text))
            both_prob = leave_prob[:, 0] * leave_prob[:, 1] * (1 - cost_array) ** 2
            line_ok, text = z_check(sum({0, 1} <= s for s in copy_sets_ls), both_prob.sum(), (both_prob * (1 - both_prob)).sum())
            ok &= line_ok; lines.append('copies at both patch2 and patch3: %s' % text)
        c.check(label, ok, '; '.join(lines))

    run_case('exact/scheme 3, is_remove=True: P(-> j), P(stay) = prod_j (1 - v_k m_1j)', 0.0, 0.0, True, None)
    run_case('exact/scheme 4, c 0.3, is_remove=True: P(-> j), P(dead), P(stay)', 0.3, 0.3, True, 0.3)
    run_case('exact/scheme 4, (c_fem 0.1, c_mal 0.5), is_remove=True', 0.1, 0.5, True, (0.1, 0.5))
    run_case('exact/scheme 3, is_remove=False: copies independent per j', 0.0, 0.0, False, None)
    run_case('exact/scheme 4, c 0.3, is_remove=False: copies only for survivors', 0.3, 0.3, False, 0.3)

def group_order(c):
    ''' plan item 8: mean v of the emigrants per destination, in loop order; the leave-rate gap at m = 1 '''
    np = c.np
    factory, meta, source_hab, dest_habs = among_setup(c, 5)
    reps, pool_size = 100, 400
    for label, kwargs in (('scheme 1', dict()), ('scheme 2', dict(use_disp_trait=True)), ('scheme 3', dict(use_disp_trait=True, is_rescale=False)),
                          ('scheme 4 c 0.3', dict(use_disp_trait=True, is_rescale=False, disp_cost=0.3)), ('scheme 3 is_remove=False', dict(use_disp_trait=True, is_rescale=False, is_remove=False))):
        v_at_dest = [[] for _ in range(3)]
        for _ in range(reps):
            source_hab.offspring_pool[:] = [factory.make_founder('female') for _ in range(pool_size)]
            meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(0.5, **KERNEL, **kwargs)
            for j, h_object in enumerate(dest_habs): v_at_dest[j] += [indi.get_trait_value('fdisp') for indi in h_object.immigrant_pool]
            meta.meta_clear_up_offspring_and_immigrant_pool()
        summary = '; '.join('patch%d %.3f +- %.3f (n %d)' % (j + 2, np.mean(v_ls), np.std(v_ls) / np.sqrt(len(v_ls)), len(v_ls)) for j, v_ls in enumerate(v_at_dest))
        if label == 'scheme 1':
            mean_2, mean_3 = np.mean(v_at_dest[0]), np.mean(v_at_dest[1])
            pooled_se = np.sqrt(np.var(v_at_dest[0]) / len(v_at_dest[0]) + np.var(v_at_dest[1]) / len(v_at_dest[1]))
            c.check('order/scheme 1: no destination-order difference in mean v (patch2 vs patch3)', abs(mean_2 - mean_3) < 4 * pooled_se, summary)
        else:
            c.check('order/%s: mean v per destination in loop order' % label, True, summary, note=True)

    # leave-rate gap at m = 1 (single scale, Nemo semantics): P(leave) = 1 - prod_j (1 - v_k m_1j) instead of v_k
    rate_row = np.asarray(meta.get_disp_among_rate_matrix(1.0, **KERNEL))[0, 1:4]
    v_ls, stay_num = [], 0
    for _ in range(reps):
        pool_ls = [factory.make_founder('female') for _ in range(pool_size)]
        source_hab.offspring_pool[:] = list(pool_ls)
        v_ls += [indi.get_trait_value('fdisp') for indi in pool_ls]
        meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(1.0, use_disp_trait=True, is_rescale=False, **KERNEL)
        stay_num += len(source_hab.offspring_pool)
        meta.meta_clear_up_offspring_and_immigrant_pool()
    v_all = np.array(v_ls)
    stay_prob = np.prod(1 - v_all[:, None] * rate_row[None, :], axis=1)
    line_ok, text = z_check(stay_num, stay_prob.sum(), (stay_prob * (1 - stay_prob)).sum())
    c.check('order/m=1, scheme 3: P(stay) = prod_j (1 - v_k m_1j)', line_ok, 'stay %s' % text)
    high_mask = v_all >= 0.8
    c.check('order/m=1 leave-rate gap (sequential destinations)', True,
            'P(leave) overall: target sum(v)/N %.3f, expected 1 - prod %.3f; v >= 0.8: target %.3f, expected %.3f'
            % (v_all.mean(), 1 - stay_prob.mean(), v_all[high_mask].mean(), 1 - stay_prob[high_mask].mean()), note=True)

def group_fixed(c):
    ''' metaibm_nemo_params dispersal_fixed_trait / dispersal_fixed_rate: the named sex leaves with the fixed rate instead of reading its fdisp / mdisp '''
    np = c.np
    fixed_rate, sex_num = 0.6, 3000                                                         # spFM founders: fdisp 0.2, mdisp 0.8, no mutation
    fixed_female = c.species('spFixF', 'spFM', seed=7, mutation=False, metaibm_nemo_params={'dispersal_fixed_trait': 'female', 'dispersal_fixed_rate': fixed_rate})
    fixed_male = c.species('spFixM', 'spFM', seed=7, mutation=False, metaibm_nemo_params={'dispersal_fixed_trait': 'male', 'dispersal_fixed_rate': fixed_rate})
    unfixed = c.species('spFM', 'spFM', seed=7, mutation=False)
    meta = c.landscape([(0, 0), (0, 1)], habs_per_patch=1, length=4, width=4)
    patch1, source_hab, sink_hab = meta.set['patch1'], meta.set['patch1'].set['h1'], meta.set['patch2'].set['h1']

    def outcome_counts(factories, run_dispersal):
        ''' sex_num females + sex_num males of each factory in the source's offspring pool, run_dispersal() moves emigrants into the sink's immigrant pool;
            RETURN: {(species_id, gender, 'arrived' | 'stayed' | 'dead'): count}; every individual goes back to Nemo '''
        pool_ls = [factory.make_founder(gender) for factory in factories for gender in ['female'] * sex_num + ['male'] * sex_num]
        key_ls = [(indi.species_id, indi.gender) for indi in pool_ls]                       # read before: the dead are recycled
        source_hab.offspring_pool[:] = list(pool_ls)
        run_dispersal()
        arrived_id_set, stayed_id_set = {id(indi) for indi in sink_hab.immigrant_pool}, {id(indi) for indi in source_hab.offspring_pool}
        count_dict = {}
        for indi, key in zip(pool_ls, key_ls):
            outcome_key = key + ('arrived' if id(indi) in arrived_id_set else ('stayed' if id(indi) in stayed_id_set else 'dead'),)
            count_dict[outcome_key] = count_dict.get(outcome_key, 0) + 1
        meta.meta_clear_up_offspring_and_immigrant_pool()
        return count_dict

    def binomial_checks(check_id, count_dict, expected_ls):
        ''' expected_ls: (species_id, gender, outcome, probability per individual) -> one z check per line, all |z| < 4 '''
        ok, line_ls = True, []
        for species_id, gender, outcome, probability in expected_ls:
            line_ok, text = z_check(count_dict.get((species_id, gender, outcome), 0), sex_num * probability, sex_num * probability * (1 - probability))
            ok &= line_ok; line_ls.append('%s %s %s (p %.2f): %s' % (species_id, gender, outcome, probability, text))
        c.check(check_id, ok, '; '.join(line_ls))

    # F1-F4 one call of hab_sample_emigrants, scheme 3 / 4, m_ij = 1
    def hab_call(**kwargs): return lambda: sink_hab.immigrant_pool.extend(source_hab.hab_sample_emigrants(1.0, is_rescale=False, **kwargs))
    binomial_checks('F1 fixed female 0.6, sex_specific: female leaves with 0.6 (not fdisp 0.2), male with mdisp 0.8',
                    outcome_counts([fixed_female], hab_call(sex_specific=True)), [('spFixF', 'female', 'arrived', 0.6), ('spFixF', 'male', 'arrived', 0.8)])
    binomial_checks('F2 fixed male 0.6, sex_specific: male leaves with 0.6 (not mdisp 0.8), female with fdisp 0.2',
                    outcome_counts([fixed_male], hab_call(sex_specific=True)), [('spFixM', 'female', 'arrived', 0.2), ('spFixM', 'male', 'arrived', 0.6)])
    binomial_checks('F3 fixed female 0.6, sex_specific=False: male reads fdisp 0.2',
                    outcome_counts([fixed_female], hab_call(sex_specific=False)), [('spFixF', 'female', 'arrived', 0.6), ('spFixF', 'male', 'arrived', 0.2)])
    binomial_checks('F4 fixed female 0.6, disp_cost 0.5: the fixed sex pays the cost too',
                    outcome_counts([fixed_female], hab_call(sex_specific=True, disp_cost=0.5)),
                    [('spFixF', 'female', 'arrived', 0.3), ('spFixF', 'female', 'dead', 0.3), ('spFixF', 'male', 'arrived', 0.4), ('spFixF', 'male', 'dead', 0.4)])

    # F5 scheme 2, e_ij = 1: the fixed rate is the draw weight
    pool_sex = [fixed_female.make_founder(g) for g in ['female'] * 50 + ['male'] * 50]
    reps, male_num = 20000, 0
    for _ in range(reps):
        source_hab.offspring_pool[:] = pool_sex
        male_num += source_hab.hab_sample_emigrants(1, sex_specific=True)[0].gender == 'male'
    expected_share = 0.8 / (0.8 + fixed_rate)
    c.check('F5 scheme 2, fixed female 0.6: share of males ~ 0.8 / (0.8 + 0.6)', abs(male_num / reps - expected_share) < 4 * np.sqrt(expected_share * (1 - expected_share) / reps),
            'sim %.4f / exp %.4f' % (male_num / reps, expected_share))
    source_hab.offspring_pool[:] = []
    for indi in pool_sex: indi.recycle()

    # F6 two species in one pool: each by its own metaibm_nemo_params
    binomial_checks('F6 one pool, fixed + unfixed species: each species by its own setting',
                    outcome_counts([fixed_female, unfixed], hab_call(sex_specific=True)),
                    [('spFixF', 'female', 'arrived', 0.6), ('spFixF', 'male', 'arrived', 0.8), ('spFM', 'female', 'arrived', 0.2), ('spFM', 'male', 'arrived', 0.8)])

    # F7 patch_sample_emigrants; F8 the among function at m = 1 (two patches: m_12 = 1)
    binomial_checks('F7 patch_sample_emigrants, fixed female 0.6',
                    outcome_counts([fixed_female], lambda: sink_hab.immigrant_pool.extend(patch1.patch_sample_emigrants(1.0, sex_specific=True, is_rescale=False))),
                    [('spFixF', 'female', 'arrived', 0.6), ('spFixF', 'male', 'arrived', 0.8)])
    binomial_checks('F8 dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool, m 1, fixed female 0.6',
                    outcome_counts([fixed_female], lambda: meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(1.0, method='uniform', use_disp_trait=True, is_rescale=False, sex_specific=True)),
                    [('spFixF', 'female', 'arrived', 0.6), ('spFixF', 'male', 'arrived', 0.8)])

    # F9 use_disp_trait=False: the same individuals leave as for the unfixed species (same seeds, same pool order)
    picked_index_ls = []
    for factory in (fixed_female, unfixed):
        pool_ls = [factory.make_founder(gender) for gender in ['female'] * 500 + ['male'] * 500]
        source_hab.offspring_pool[:] = list(pool_ls)
        random.seed(11); np.random.seed(11)
        meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(0.5, method='uniform')
        arrived_id_set = {id(indi) for indi in sink_hab.immigrant_pool}
        picked_index_ls.append([k for k, indi in enumerate(pool_ls) if id(indi) in arrived_id_set])
        meta.meta_clear_up_offspring_and_immigrant_pool()
    c.check('F9 use_disp_trait=False: fixed and unfixed species lose the same individuals', picked_index_ls[0] == picked_index_ls[1], 'picked %d / %d' % (len(picked_index_ls[0]), len(picked_index_ls[1])))

    # F10 the checks Nemo makes on each parameter
    outcome_ls = []
    for label, params in (("trait 'both'", {'dispersal_fixed_trait': 'both', 'dispersal_fixed_rate': 0.1}), ('rate 1.5', {'dispersal_fixed_trait': 'female', 'dispersal_fixed_rate': 1.5}),
                          ('rate -0.1', {'dispersal_fixed_trait': 'female', 'dispersal_fixed_rate': -0.1})):
        try:
            c.species('spBad', 'spFM', seed=7, mutation=False, metaibm_nemo_params=params); outcome_ls.append((label, 'booted'))
        except ValueError:
            outcome_ls.append((label, 'ValueError'))
    for label, params in (("'male', rate 0", {'dispersal_fixed_trait': 'male', 'dispersal_fixed_rate': 0}), ('None / None', {'dispersal_fixed_trait': None, 'dispersal_fixed_rate': None})):
        try:
            c.species('spOk', 'spFM', seed=7, mutation=False, metaibm_nemo_params=params); outcome_ls.append((label, 'booted'))
        except ValueError:
            outcome_ls.append((label, 'ValueError'))
    c.check('F10 ValueError for a trait other than female / male and a rate outside [0, 1]; valid values boot',
            [outcome for _, outcome in outcome_ls] == ['ValueError'] * 3 + ['booted'] * 2, '; '.join('%s: %s' % pair for pair in outcome_ls))

GROUP_FUNCS = {'units': (group_units, 1), 'exact': (group_exact, 5), 'order': (group_order, 5), 'fixed': (group_fixed, 7)}

##****************************** parent side: shared runner ******************************
def run_child(script, group, root=PROJECT_ROOT, args=(), timeout=1200):
    ''' run one group of `script` in a fresh child process; RETURN: (checks, data, raw output, exit code, seconds) '''
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    t0 = time.time()
    try:
        result = subprocess.run([sys.executable, '-I', '-B', str(script), '--child', group, str(root)] + [str(a) for a in args],
                                capture_output=True, text=True, env=env, cwd=tempfile.gettempdir(), timeout=timeout)
        raw, code = result.stdout + result.stderr, result.returncode
    except subprocess.TimeoutExpired as error:
        raw, code = (error.stdout.decode() if isinstance(error.stdout, bytes) else (error.stdout or '')), 'timeout %ds' % timeout
    checks = [json.loads(line[3:]) for line in raw.splitlines() if line.startswith('@@ ')]
    data = {d['key']: d['value'] for d in (json.loads(line[7:]) for line in raw.splitlines() if line.startswith('@@DATA '))}
    return checks, data, raw, code, time.time() - t0

def md5(path): return hashlib.md5(Path(path).read_bytes()).hexdigest()

def print_header():
    import pynemo
    so_name = [f for f in os.listdir(os.path.dirname(pynemo.__file__)) if f.startswith('_nemo')][0]
    print('python %s | pynemo %s (%s md5 %s)' % (sys.executable, os.path.dirname(pynemo.__file__), so_name, md5(Path(pynemo.__file__).parent / so_name)))
    for rel in ('metaibm/habitat.py', 'metaibm/patch.py', 'metaibm/metacommunity.py', 'extension/global_habitat_network.py'):
        print('  %-38s md5 %s' % (rel, md5(PROJECT_ROOT / rel)))

def parent_report(results, t_all):
    ''' print the markdown table; RETURN: exit code (1 on any FAIL) '''
    counts = {s: sum(1 for x in results if x['status'] == s) for s in ('PASS', 'FAIL', 'NOTE')}
    print('\n| status | id | detail |\n|---|---|---|')
    for x in results: print('| %s | `%s` | %s |' % (x['status'], x['id'], x['detail'].replace('|', '\\|')))
    print('\nPASS %d  FAIL %d  NOTE %d  (%.0fs)' % (counts['PASS'], counts['FAIL'], counts['NOTE'], time.time() - t_all))
    return 1 if counts['FAIL'] else 0

def add_child(results, group, checks, code, raw):
    results.extend(checks)
    if code != 0: results.append({'group': group, 'id': '%s/child-exit' % group, 'status': 'FAIL', 'detail': 'exit %s: %s' % (code, raw[-600:].replace('\n', ' '))})

def main(groups):
    print_header()
    results, t_all = [], time.time()
    for group in groups:
        checks, _, raw, code, seconds = run_child(Path(__file__).resolve(), group)
        add_child(results, group, checks, code, raw)
        print('  %-10s %5.1fs  %s' % (group, seconds, ' '.join(x['status'][0] for x in checks)), flush=True)
    return parent_report(results, t_all)

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--child':
        group_name, root = sys.argv[2], sys.argv[3]
        group_func, seed = GROUP_FUNCS[group_name]
        ctx = Ctx(group_name, root, seed)
        try: group_func(ctx)
        finally: ctx.close()
    else:
        sys.exit(main(sys.argv[1:] or list(GROUPS)))
