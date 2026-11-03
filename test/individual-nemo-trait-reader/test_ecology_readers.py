# -*- coding: utf-8 -*-
''' test_ecology_readers.py -- get_trait_value / get_trait_alleles through MetaIBM's own life cycle (metacommunity -> patch -> habitat -> indfactory),
    where individuals are born, die, are recycled and their Nemo objects reused every step.

    Run (conda env metaibm-nemo-test6):
        PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test6/bin/python -B test_ecology_readers.py [group ...]

    Landscape (example2-like, smaller): a mainland of one patch with four 10 x 10 habitats and four one-habitat 10 x 10 islands;
    species sp1 / sp2 alternate over the habitats (matched by the habitat mean environment = the species' quanti_init_value).
    Both species carry all 7 Nemo traits with the same layout (quant continuous 2 x 20 loci, ntrl 12, delet 30, dmi 8 = 4 pairs, disp, wolb),
    every mutation on; quant Ve 0; wolb transmission 1, half of the founder females infected after meta_initialize.
    Every step, every living individual: get_trait_value == stored == oracle, get_trait_alleles == raw (oracles from test_trait_readers.py).
'''
import sys
sys.dont_write_bytecode = True
import os, json, time, hashlib, random
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from test_trait_readers import (Ctx, PROJECT_ROOT, DBL, VALID_PAIRS, DMI16, quant_rows, ntrl_rows, delet_rows, dmi_rows, disp_rows, wolb_rows,
                                raw_table, stored, delet_formula, dmi_table16, dmi_formula, raises)

GROUPS = ('cycle_asexual', 'cycle_sexual', 'cycle_mix', 'recycled_refs', 'repro')
ENVS = [[0.3, -0.2], [-0.1, 0.2]]               # founder mean phenotype of sp1 / sp2 = habitat mean environment
S_, H_ = 0.05, 0.3                              # delet effects / dominance (delet_rows defaults)
T16 = dmi_table16(DMI16)

##****************************** child side ******************************
def species_traits(k, mut=0.25):
    q = quant_rows('continuous', mut=0.01 * mut, extra=[('quanti_init_value', '{{%r, %r}}' % tuple(ENVS[k]))])
    return {'quant': q, 'ntrl': ntrl_rows(4, 0.05 * mut), 'delet': delet_rows(0.2 * mut, model=1, init=0.3, s=S_, h=H_), 'dmi': dmi_rows(0.2 * mut), 'disp': disp_rows(0.2 * mut, 0.3), 'wolb': wolb_rows(1.0)}

def build(c, seed=12):
    import metaibm
    facs = {}
    for k in range(2):
        f = c.species(species_traits(k), seed=seed + k)
        f.species_id = 'sp%d' % (k + 1)                                    # habitat.hab_initialize maps species by 'sp<k+1>'
        facs[f.species_id] = f
    mainland = metaibm.metacommunity(metacommunity_name='mainland')
    patch = metaibm.patch(patch_name='mainland', patch_index=0, location=(0, 0))
    for h in range(4):
        patch.add_habitat(hab_name='h%d' % (h + 1), hab_index=h, hab_location=(h // 2, h % 2), num_env_types=2, env_types_name=('env1', 'env2'),
                          mean_env_ls=list(ENVS[h % 2]), var_env_ls=[0.025 ** 2] * 2, length=10, width=10, global_index=4 + h)
    mainland.add_patch(patch_name='mainland', patch_object=patch)
    islands = metaibm.metacommunity(metacommunity_name='islands')
    for h in range(4):
        p = metaibm.patch(patch_name='patch%d' % (h + 1), patch_index=h, location=(h // 2, h % 2))
        p.add_habitat(hab_name='h%d' % (h + 1), hab_index=0, hab_location=(h // 2, h % 2), num_env_types=2, env_types_name=('env1', 'env2'),
                      mean_env_ls=list(ENVS[h % 2]), var_env_ls=[0.025 ** 2] * 2, length=10, width=10, global_index=h)
        islands.add_patch(patch_name='patch%d' % (h + 1), patch_object=p)
    return facs, mainland, islands

def individuals(meta):
    for pid, p in meta.set.items():
        for hid, h in p.set.items():
            for (r, cc) in sorted(h.occupied_site_pos_ls):
                yield (meta.metacommunity_name, pid, hid, r, cc), h.set['microsite_individuals'][r][cc]

def invariants(c, metas):
    ''' every living individual: both readers against stored + oracles; RETURN (n, bad, stats) '''
    np = c.np
    bad, n, infected, W = [], 0, 0, []
    for meta in metas:
        for where, z in individuals(meta):
            n += 1
            if z.indfactory.species_id != z.species_id: bad.append('%s species/factory' % (where,))
            legacy = z.get_indi_phenotype_ls()
            for d in range(2):
                v = z.get_trait_value('quant', trait_dim=d)
                if not (v == legacy[d] == stored(z, 'quant', d)): bad.append('%s quant%d value' % (where, d))
            st = z.get_trait_alleles('delet', 'get_allele'); v = z.get_trait_value('delet')
            if not (abs(v - delet_formula(st, S_, H_, 1)[0]) < 1e-12 and v == stored(z, 'delet')): bad.append('%s delet' % (where,))
            W.append(v)
            v = z.get_trait_value('dmi')
            if not (abs(v - dmi_formula(z.get_trait_alleles('dmi', 'get_allele'), T16)[0]) < 1e-12 and v == stored(z, 'dmi')): bad.append('%s dmi' % (where,))
            for k in ('fdisp', 'mdisp'):
                a = z.get_trait_alleles(k, 'get_allele_value'); v = z.get_trait_value(k)
                if not (abs(v - a.mean()) < 1e-15 and v == stored(z, k) and a.min() >= 0 and a.max() <= 1): bad.append('%s %s' % (where, k))
            w = z.get_trait_value('wolb'); infected += w
            if not (type(w) is bool and w == z.get_trait_alleles('wolb', 'get_allele') == stored(z, 'wolb')): bad.append('%s wolb' % (where,))
            if z.get_trait_value('ntrl') is not None: bad.append('%s ntrl' % (where,))
            for k, m in VALID_PAIRS:
                if k == 'quant':
                    if not all(np.array_equal(z.get_trait_alleles(k, m, trait_dim=d), raw_table(np, z, k, m, d)) for d in range(2)): bad.append('%s alleles quant %s' % (where, m))
                elif k != 'wolb' and not np.array_equal(z.get_trait_alleles(k, m), raw_table(np, z, k, m)): bad.append('%s alleles %s' % (where, k))
    return n, bad, {'alive': n, 'infected': round(infected / n, 3) if n else None, 'delet_W_mean': round(float(np.mean(W)), 4) if W else None}

def infect_half(meta):
    k = 0
    for where, z in individuals(meta):
        if z.gender == 'female':
            if k % 2 == 0: z.get_trait('wolb').set_allele_value(0, 0, 1.0)
            k += 1

def digest(c, metas):
    np = c.np
    rows = []
    for meta in metas:
        for where, z in individuals(meta):
            r = [list(map(str, where)), z.species_id, z.individual_id, [repr(z.get_trait_value('quant', trait_dim=d)) for d in range(2)]] + [repr(z.get_trait_value(k)) for k in DBL + ('wolb',)]
            r += [z.get_trait_alleles(k, m, trait_dim=0).tolist() if k == 'quant' else (z.get_trait_alleles(k, m) if k == 'wolb' else z.get_trait_alleles(k, m).tolist()) for k, m in VALID_PAIRS]
            rows.append(r)
    return hashlib.md5(json.dumps(rows).encode()).hexdigest()

def cycle(c, mode, steps=30, on_step=None):
    ''' example2's life cycle; invariants every step; RETURN (digest, bad, trajectory, (minted, recycled) per species) '''
    np = c.np
    np.random.seed(7); random.seed(7)                                       # before build(): add_habitat draws the microsite environments from np.random
    facs, mainland, islands = build(c)
    mainland.meta_initialize(facs, ENVS, reproduce_mode='asexual' if mode == 'asexual' else 'sexual')
    infect_half(mainland)
    n0, bad0, st0 = invariants(c, [mainland])
    all_bad, traj = list(bad0), [(0, st0)]
    for t in range(1, steps + 1):
        before = [z for meta in (mainland, islands) for _, z in individuals(meta)]
        mainland.meta_dead_selection(0.1, 0.5, method='niche_gaussian')
        if mode == 'asexual': mainland.meta_mainland_asexual_birth_mutate_germinate(0.5)
        elif mode == 'sexual': mainland.meta_mainland_sexual_birth_mutate_germinate(1)
        else: mainland.meta_mainland_mixed_birth_mutate_germinate(0.5, 1)
        islands.meta_dead_selection(0.1, 0.5, method='niche_gaussian')
        if mode == 'asexual': islands.meta_asex_reproduce_mutate_into_offspring_pool(0.5)
        elif mode == 'sexual': islands.meta_sex_reproduce_mutate_into_offspring_pool(1)
        else: islands.meta_mix_reproduce_mutate_into_offspring_pool(0.5, 1)
        islands.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(0.1, method='uniform', is_remove=True)
        if mode == 'asexual': islands.meta_colonize_from_propagules_rains(mainland, 10, is_remove=True)
        else: islands.pairwise_sexual_colonization_from_prpagules_rains(mainland, 10, is_remove=True)
        islands.meta_local_germinate_from_offspring_and_immigrant_pool(is_remove=True)
        islands.meta_clear_up_offspring_and_immigrant_pool()
        if on_step: on_step(t, before)
        n, bad, st = invariants(c, [mainland, islands])
        all_bad += ['t=%d %s' % (t, b) for b in bad]
        if t % 10 == 0: traj.append((t, st))
    d = digest(c, [mainland, islands])
    for meta in (mainland, islands):
        meta.meta_clear_up_microsite_individuals(); meta.meta_clear_up_offspring_and_immigrant_pool(); meta.meta_clear_up_dormancy_pool()
    balanced = {s: (f.minted_num, f.recycled_num) for s, f in facs.items()}
    for f in facs.values(): f.purge_recycling_pool()
    return d, all_bad, traj, balanced

def _g_cycle(c, mode):
    d, bad, traj, balanced = cycle(c, mode)
    c.check('ecology/%s/readers-every-step' % mode, not bad, '30 steps, mainland + islands, every living individual every step: get_trait_value == stored == oracle (quant legacy, delet / dmi formula, disp mean in [0, 1], wolb), get_trait_alleles == raw for 8 pairs; (t, stats) %s' % traj if not bad else '; '.join(bad[:5]))
    c.check('ecology/%s/minted==recycled' % mode, all(m == r for m, r in balanced.values()), str(balanced))
    c.data('digest/%s' % mode, d)

def g_cycle_asexual(c): _g_cycle(c, 'asexual')
def g_cycle_sexual(c): _g_cycle(c, 'sexual')
def g_cycle_mix(c): _g_cycle(c, 'mix')

def g_recycled_refs(c):
    ''' wrappers of individuals that died during a step: both readers raise RuntimeError on them, while their Nemo objects carry newborns '''
    stats = {'dead': 0, 'bad': []}
    def on_step(t, before):
        for z in before:
            if z.is_alive: continue
            stats['dead'] += 1
            e = [raises(z.get_trait_value, k) for k in ('quant', 'delet', 'wolb', 'ntrl')] + [raises(z.get_trait_alleles, 'quant', 'get_allele_value'), raises(z.get_trait_alleles, 'dmi', 'get_allele')]
            if not all(x and x.startswith('RuntimeError') for x in e): stats['bad'].append('t=%d %r' % (t, e))
    d, bad, traj, balanced = cycle(c, 'sexual', steps=15, on_step=on_step)
    c.check('ecology/recycled_refs', stats['dead'] > 0 and not stats['bad'] and not bad, '%d wrappers of individuals that died in their step: get_trait_value x4 and get_trait_alleles x2 -> RuntimeError; living ones stay correct'
            % stats['dead'] if not stats['bad'] else '; '.join(stats['bad'][:3]))

def g_repro(c):
    d, bad, _, _ = cycle(c, 'sexual', steps=15)
    c.data('digest/repro', d)
    c.check('ecology/repro/run', not bad, 'one run of 15 steps', note=True)

##****************************** parent side ******************************
def main(groups):
    results, t_all = [], time.time()
    for g in groups:
        if g == 'repro':
            _, d1, r1, k1, t1 = run_child_eco('repro'); _, d2, r2, k2, t2 = run_child_eco('repro')
            ok = k1 == 0 and k2 == 0 and d1.get('digest/repro') and d1.get('digest/repro') == d2.get('digest/repro')
            results.append({'group': 'repro', 'id': 'ecology/repro/cross-process', 'status': 'PASS' if ok else 'FAIL',
                            'detail': 'sexual life cycle, same seeds, two processes: landscape digest (values + allele tables) %s / %s' % (str(d1.get('digest/repro'))[:12], str(d2.get('digest/repro'))[:12])})
            print('  %-16s %5.1fs  %s' % (g, t1 + t2, results[-1]['status'][0]), flush=True)
            continue
        ch, _, raw, code, dt = run_child_eco(g)
        results += ch
        if code != 0: results.append({'group': g, 'id': '%s/child-exit' % g, 'status': 'FAIL', 'detail': 'exit %s: %s' % (code, raw[-600:].replace('\n', ' '))})
        print('  %-16s %5.1fs  %s' % (g, dt, ' '.join(x['status'][0] for x in ch)), flush=True)
    n = {s: sum(1 for x in results if x['status'] == s) for s in ('PASS', 'FAIL', 'NOTE')}
    print('\n| status | id | detail |\n|---|---|---|')
    for x in results: print('| %s | `%s` | %s |' % (x['status'], x['id'], x['detail'].replace('|', '\\|')))
    print('\nPASS %d  FAIL %d  NOTE %d  (%.0fs)' % (n['PASS'], n['FAIL'], n['NOTE'], time.time() - t_all))
    return 1 if n['FAIL'] else 0

def run_child_eco(group):
    ''' run_child of test_trait_readers.py, pointed at this file '''
    import subprocess, tempfile
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    t0 = time.time()
    r = subprocess.run([sys.executable, '-I', '-B', str(Path(__file__).resolve()), '--child', group, str(PROJECT_ROOT)], capture_output=True, text=True, env=env, cwd=tempfile.gettempdir(), timeout=1800)
    raw = r.stdout + r.stderr
    checks = [json.loads(l[3:]) for l in raw.splitlines() if l.startswith('@@ ')]
    data = {d['key']: d['value'] for d in (json.loads(l[7:]) for l in raw.splitlines() if l.startswith('@@DATA '))}
    return checks, data, raw, r.returncode, time.time() - t0

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--child':
        ctx = Ctx(sys.argv[2], Path(sys.argv[3]))
        try: globals()['g_' + sys.argv[2]](ctx)
        finally: ctx.close()
    else:
        sys.exit(main([g for g in sys.argv[1:]] or list(GROUPS)))
