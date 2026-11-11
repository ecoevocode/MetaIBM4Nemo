# -*- coding: utf-8 -*-
''' test_disp_lifecycle.py -- disp-trait sampling inside whole life cycles: the 7 dispersal functions (3 among, 3 within, 1 network)

    Run (pynemo with TTrait.get_value -> conda env metaibm-nemo-test6):
        PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test6/bin/python -B test_disp_lifecycle.py [group ...]

    Groups:
        regression_head   default path (no new kwargs): HEAD's habitat / patch / metacommunity / global_habitat_network vs the working tree,
                          a seeded life cycle through all 7 dispersal functions (is_remove True / False), full-state md5 must be identical
        cost_none         disp_cost=0 / (0, 0) == None (scheme 3); with is_rescale=True (scheme 2) disp_cost has no effect -- digest equalities
        bookkeeping       schemes 2 / 3 / 4 x sex_specific x cost x sexual / asexual: every individual is recycled or sits in exactly one place
        full_sim          range expansion over 10 patches, 150 steps, 5 seeds x 4 schemes: migrant numbers, migrant v, spread, fdisp evolution

    Each life cycle / simulation runs in its own child process with the same seeds (Nemo's RAND and Individual::currentID are process-global).
    Shares the harness of test_disp_sampling.py (Ctx, run_child, parent_report); Nemo writes only into temporary directories.
'''
import sys
sys.dont_write_bytecode = True                  # never drop __pycache__ into metaibm/ or this folder
import os, json, time, shutil, hashlib, tempfile, subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from test_disp_sampling import Ctx, PROJECT_ROOT, run_child, add_child, parent_report, print_header

GROUPS = ('regression_head', 'cost_none', 'bookkeeping', 'full_sim')
CHANGED_FILES = ('metaibm/habitat.py', 'metaibm/patch.py', 'metaibm/metacommunity.py', 'extension/global_habitat_network.py')
LIFECYCLE_SEED, LIFECYCLE_STEPS = 1, 16
CONFIGS = (                                                                       # (label, dispersal kwargs) for every life cycle
    ('scheme 1 (default)', {}),
    ('scheme 2', {'use_disp_trait': True}),
    ('scheme 2 sex_specific', {'use_disp_trait': True, 'sex_specific': True}),
    ('scheme 2 + disp_cost 0.3', {'use_disp_trait': True, 'disp_cost': 0.3}),
    ('scheme 2 + disp_cost (0.1, 0.5)', {'use_disp_trait': True, 'disp_cost': [0.1, 0.5]}),
    ('scheme 3', {'use_disp_trait': True, 'is_rescale': False}),
    ('scheme 3 sex_specific', {'use_disp_trait': True, 'is_rescale': False, 'sex_specific': True}),
    ('scheme 3 disp_cost 0', {'use_disp_trait': True, 'is_rescale': False, 'disp_cost': 0}),
    ('scheme 3 disp_cost (0, 0)', {'use_disp_trait': True, 'is_rescale': False, 'disp_cost': [0, 0]}),
    ('scheme 4 c 0.3', {'use_disp_trait': True, 'is_rescale': False, 'disp_cost': 0.3}),
    ('scheme 4 sex_specific (0.1, 0.5)', {'use_disp_trait': True, 'is_rescale': False, 'sex_specific': True, 'disp_cost': [0.1, 0.5]}),
)
FULL_SCHEMES = (('scheme 1', {}), ('scheme 2', {'use_disp_trait': True}), ('scheme 3', {'use_disp_trait': True, 'is_rescale': False}),
                ('scheme 4 c 0.3', {'use_disp_trait': True, 'is_rescale': False, 'disp_cost': 0.3}))
FULL_SEEDS, FULL_STEPS, FULL_PATCHES, HALF_FULL = (1, 2, 3, 4, 5), 150, 10, 144
M_AMONG, M_WITHIN = 0.05, 0.2

##****************************** child side ******************************
def state_digest(meta_obj, md5):
    ''' feed the full state (residents with id / species / fdisp / mdisp, and every pool) into md5 '''
    for patch_name, patch_obj in meta_obj.set.items():
        for h_id, h_object in patch_obj.set.items():
            for row_id, row in enumerate(h_object.set['microsite_individuals']):
                for col_id, indi in enumerate(row):
                    if indi is not None:
                        md5.update(('%s|%s|%d|%d|%d|%s|%.15g|%.15g;' % (patch_name, h_id, row_id, col_id, indi.individual_id, indi.species_id,
                                    indi.get_trait_value('fdisp'), indi.get_trait_value('mdisp'))).encode())
            for pool_name in ('offspring_pool', 'immigrant_pool', 'dormancy_pool'):
                md5.update(('%s|%s|%s:%s;' % (patch_name, h_id, pool_name, ','.join(str(indi.individual_id) for indi in getattr(h_object, pool_name)))).encode())

def child_lifecycle(c, mode, extra):
    ''' a seeded life cycle that rotates through all 7 dispersal functions (is_remove True / False); RETURN via @@DATA: digest + bookkeeping '''
    factories = {'spU': c.species('spU', 'spU', seed=LIFECYCLE_SEED), 'spFM': c.species('spFM', 'spFM', seed=LIFECYCLE_SEED)}
    meta = c.landscape([(0, 0), (0, 1), (1, 0), (1, 1)], habs_per_patch=2, length=8, width=8, dormancy=10)
    c.fill(meta, factories, sexual=(mode == 'sexual'))
    md5 = hashlib.md5(); state_digest(meta, md5)
    kernel = dict(method='gaussian', sigma=1.5)
    for step in range(LIFECYCLE_STEPS):
        meta.meta_dead_selection(0.2, None, method='neutral_uniform')
        if mode == 'sexual': meta.meta_sex_reproduce_mutate_into_offspring_pool(0.5)
        else:                meta.meta_asex_reproduce_mutate_into_offspring_pool(0.5)
        variant, is_remove = step % 4, (step // 4) % 2 == 0
        if variant == 0:
            meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(0.3, is_remove=is_remove, **kernel, **extra)
            meta.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(0.3, is_remove=is_remove, **extra)
            meta.meta_local_germinate_from_offspring_and_immigrant_pool()
        elif variant == 1:
            meta.dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool(0.3, is_remove=is_remove, **kernel, **extra)
            meta.meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool(0.3, is_remove=is_remove, **extra)
            meta.meta_local_germinate_from_offspring_immigrant_and_dormancy_pool()
        elif variant == 2:
            meta.dispersal_among_patches_from_offspring_pool_and_dormancy_pool(0.3, is_remove=is_remove, **kernel, **extra)
            meta.meta_dispersal_within_patch_from_offspring_and_dormancy_pool(0.3, is_remove=is_remove, **extra)
            meta.meta_local_germinate_from_offspring_and_dormancy_pool()
        else:
            meta.dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool(0.3, is_remove=is_remove, **kernel, **extra)
            meta.meta_local_germinate_from_offspring_and_immigrant_pool()
        meta.meta_dormancy_process_from_offspring_pool_to_dormancy_pool()
        state_digest(meta, md5)
    id_ls = [indi.individual_id for p in meta.set.values() for h in p.set.values() for row in h.set['microsite_individuals'] for indi in row if indi is not None]
    id_ls += [indi.individual_id for p in meta.set.values() for h in p.set.values() for pool in ('offspring_pool', 'immigrant_pool', 'dormancy_pool') for indi in getattr(h, pool)]
    alive_num = sum(f.minted_num - f.recycled_num for f in factories.values())
    c.data('digest', md5.hexdigest())
    c.data('placed', len(id_ls)); c.data('unique', len(id_ls) == len(set(id_ls))); c.data('alive', alive_num)

def child_full_sim(c, extra):
    ''' range expansion along 10 patches (only patch1 seeded), sexual, disp traits mutate and are inherited; RETURN via @@DATA: the time series '''
    np = c.np
    factories = {'spU': c.species('spU', 'spU', seed=c.seed)}
    meta = c.landscape([(0, k) for k in range(FULL_PATCHES)], habs_per_patch=2, length=12, width=12)
    c.fill(meta, factories, patch_names=['patch1'])
    def residents(patch_index):
        return [indi for h in meta.set['patch%d' % (patch_index + 1)].set.values() for row in h.set['microsite_individuals'] for indi in row if indi is not None]
    def pool_v(pool_name):
        return [indi.get_trait_value('fdisp') for p in meta.set.values() for h in p.set.values() for indi in getattr(h, pool_name)]
    record = {'N': [], 'v': [], 'mig_ratio': [], 'mig_v': [], 'pool_v': []}
    for step in range(FULL_STEPS + 1):
        record['N'].append([len(residents(p)) for p in range(FULL_PATCHES)])
        record['v'].append([float(np.mean([indi.get_trait_value('fdisp') for indi in residents(p)])) if residents(p) else None for p in range(FULL_PATCHES)])
        if step == FULL_STEPS: break
        meta.meta_dead_selection(0.1, None, method='neutral_uniform')
        meta.meta_sex_reproduce_mutate_into_offspring_pool(1.0)
        offspring_v = pool_v('offspring_pool')
        meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(M_AMONG, method='gaussian', sigma=1.0, **extra)
        migrant_v = pool_v('immigrant_pool')                                   # the among-patch migrants of this step
        if offspring_v:
            record['mig_ratio'].append(len(migrant_v) / (M_AMONG * len(offspring_v)))    # migrants / MetaIBM's expected number
            record['pool_v'].append(float(np.mean(offspring_v)))
            record['mig_v'].append(float(np.mean(migrant_v)) if migrant_v else None)
        meta.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(M_WITHIN, **extra)
        meta.meta_local_germinate_from_offspring_and_immigrant_pool()
        meta.meta_clear_up_offspring_and_immigrant_pool()
    c.data('record', record)

##****************************** parent side ******************************
def head_copy():
    ''' metaibm/ + extension/ of the working tree with HEAD's 4 changed files, in a temp dir '''
    copy_dir = Path(tempfile.mkdtemp(prefix='disp_head_'))
    for pkg in ('metaibm', 'extension'): shutil.copytree(PROJECT_ROOT / pkg, copy_dir / pkg, ignore=shutil.ignore_patterns('__pycache__'))
    for rel in CHANGED_FILES:
        head_text = subprocess.run(['git', 'show', 'HEAD:' + rel], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True).stdout
        (copy_dir / rel).write_text(head_text)
    return copy_dir

class Runs:
    ''' memoised child runs: one life cycle per (root, mode, kwargs) '''
    def __init__(self, results):
        self.results, self.cache = results, {}
    def lifecycle(self, root, mode, extra):
        key = (str(root), mode, json.dumps(extra, sort_keys=True))
        if key not in self.cache:
            checks, data, raw, code, _ = run_child(Path(__file__).resolve(), 'lifecycle', root, args=(mode, json.dumps(extra)))
            add_child(self.results, 'lifecycle', checks, code, raw)
            data['bookkeeping_ok'] = bool(data.get('unique')) and data.get('placed') == data.get('alive')
            self.cache[key] = data
        return self.cache[key]

def check(results, group, check_id, ok, detail, note=False):
    results.append({'group': group, 'id': check_id, 'status': 'NOTE' if note else ('PASS' if ok else 'FAIL'), 'detail': detail})

def short(data): return '%s (%s placed, bookkeeping %s)' % (str(data.get('digest'))[:12], data.get('placed'), 'OK' if data.get('bookkeeping_ok') else 'MISMATCH')

def group_regression_head(results, runs):
    head_dir = head_copy()
    try:
        for mode in ('sexual', 'asexual'):
            head_data, tree_data = runs.lifecycle(head_dir, mode, {}), runs.lifecycle(PROJECT_ROOT, mode, {})
            check(results, 'regression_head', 'regression_head/%s: default path, HEAD == working tree' % mode,
                  head_data.get('digest') and head_data.get('digest') == tree_data.get('digest') and head_data['bookkeeping_ok'] and tree_data['bookkeeping_ok'],
                  'HEAD %s / working tree %s' % (short(head_data), short(tree_data)))
    finally:
        shutil.rmtree(head_dir, ignore_errors=True)

def group_cost_none(results, runs):
    config = dict(CONFIGS)
    for mode in ('sexual', 'asexual'):
        scheme3 = [runs.lifecycle(PROJECT_ROOT, mode, config[label]) for label in ('scheme 3', 'scheme 3 disp_cost 0', 'scheme 3 disp_cost (0, 0)')]
        check(results, 'cost_none', 'cost_none/%s: scheme 3, disp_cost None == 0 == (0, 0)' % mode,
              scheme3[0].get('digest') and len({d.get('digest') for d in scheme3}) == 1, ' / '.join(short(d) for d in scheme3))
        scheme2 = [runs.lifecycle(PROJECT_ROOT, mode, config[label]) for label in ('scheme 2', 'scheme 2 + disp_cost 0.3', 'scheme 2 + disp_cost (0.1, 0.5)')]
        check(results, 'cost_none', 'cost_none/%s: scheme 2 ignores disp_cost (none == 0.3 == (0.1, 0.5))' % mode,
              scheme2[0].get('digest') and len({d.get('digest') for d in scheme2}) == 1, ' / '.join(short(d) for d in scheme2))

def group_bookkeeping(results, runs):
    for mode in ('sexual', 'asexual'):
        for label, extra in CONFIGS:
            data = runs.lifecycle(PROJECT_ROOT, mode, extra)
            check(results, 'bookkeeping', 'bookkeeping/%s: %s' % (mode, label), data['bookkeeping_ok'],
                  '%s placed, ids unique %s, minted - recycled %s' % (data.get('placed'), data.get('unique'), data.get('alive')))
    config = dict(CONFIGS)
    same = all(runs.lifecycle(PROJECT_ROOT, 'asexual', config[a]).get('digest') == runs.lifecycle(PROJECT_ROOT, 'asexual', config[b]).get('digest')
               for a, b in (('scheme 2', 'scheme 2 sex_specific'), ('scheme 3', 'scheme 3 sex_specific')))
    check(results, 'bookkeeping', 'bookkeeping/asexual: sex_specific has no effect (asexual founders are all female, mdisp is never read)', True,
          'digest with sex_specific == without, schemes 2 and 3: %s' % same, note=True)

def group_full_sim(results):
    import numpy as np
    summary = {}
    for label, extra in FULL_SCHEMES:
        record_ls = []
        for seed in FULL_SEEDS:
            checks, data, raw, code, _ = run_child(Path(__file__).resolve(), 'full_sim', PROJECT_ROOT, args=(seed, json.dumps(extra)))
            add_child(results, 'full_sim', checks, code, raw)
            if 'record' in data: record_ls.append(data['record'])
        if len(record_ls) != len(FULL_SEEDS): continue
        ratio = float(np.mean([np.mean(r['mig_ratio']) for r in record_ls]))
        delta_v = float(np.mean([np.nanmean([m - p for m, p in zip(r['mig_v'], r['pool_v']) if m is not None]) for r in record_ls]))
        half_full_step = np.array([[next((t for t, n in enumerate(r['N']) if n[p] >= HALF_FULL), np.nan) for p in range(FULL_PATCHES)] for r in record_ls], dtype=float)
        def meta_v(r, t): return float(np.sum([n * v for n, v in zip(r['N'][t], r['v'][t]) if v is not None]) / np.sum(r['N'][t]))
        fdisp_time = [round(float(np.mean([meta_v(r, t) for r in record_ls])), 3) for t in (0, 50, 100, 150)]
        summary[label] = dict(ratio=ratio, delta_v=delta_v)
        check(results, 'full_sim', 'full_sim/%s: summary over %d seeds' % (label, len(FULL_SEEDS)), True,
              'migrants / (m x offspring) %.3f; migrant v - pool v %+.3f; step patch 2..10 is half full %s; mean fdisp t=0/50/100/150 %s'
              % (ratio, delta_v, np.round(np.nanmean(half_full_step, axis=0)[1:], 1).tolist(), fdisp_time), note=True)
    if len(summary) != len(FULL_SCHEMES): return
    s1, s2, s3, s4 = (summary[label] for label, _ in FULL_SCHEMES)
    check(results, 'full_sim', 'full_sim/scheme 1: migrant number = m x offspring, migrants are a random sample of v',
          abs(s1['ratio'] - 1) < 0.03 and abs(s1['delta_v']) < 0.01, 'ratio %.3f, migrant v - pool v %+.3f' % (s1['ratio'], s1['delta_v']))
    check(results, 'full_sim', 'full_sim/scheme 2: migrant number unchanged, migrants have higher v',
          abs(s2['ratio'] - 1) < 0.03 and s2['delta_v'] > 0.02, 'ratio %.3f, migrant v - pool v %+.3f' % (s2['ratio'], s2['delta_v']))
    check(results, 'full_sim', 'full_sim/schemes 3 and 4: migrants have higher v; the cost lowers the migrant number',
          s3['delta_v'] > 0.02 and s4['delta_v'] > 0.02 and s4['ratio'] < s3['ratio'],
          'scheme 3 ratio %.3f, v diff %+.3f; scheme 4 ratio %.3f, v diff %+.3f' % (s3['ratio'], s3['delta_v'], s4['ratio'], s4['delta_v']))

def main(groups):
    print_header()
    results, t_all = [], time.time()
    runs = Runs(results)
    for group in groups:
        t0, before = time.time(), len(results)
        if group == 'regression_head': group_regression_head(results, runs)
        elif group == 'cost_none':     group_cost_none(results, runs)
        elif group == 'bookkeeping':   group_bookkeeping(results, runs)
        elif group == 'full_sim':      group_full_sim(results)
        print('  %-16s %5.1fs  %s' % (group, time.time() - t0, ' '.join(x['status'][0] for x in results[before:])), flush=True)
    return parent_report(results, t_all)

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--child':
        group_name, root = sys.argv[2], sys.argv[3]
        if group_name == 'lifecycle':
            ctx = Ctx('lifecycle', root, LIFECYCLE_SEED)
            try: child_lifecycle(ctx, sys.argv[4], json.loads(sys.argv[5]))
            finally: ctx.close()
        else:
            seed = int(sys.argv[4])
            ctx = Ctx('full_sim', root, seed); ctx.seed = seed
            try: child_full_sim(ctx, json.loads(sys.argv[5]))
            finally: ctx.close()
    else:
        sys.exit(main(sys.argv[1:] or list(GROUPS)))
