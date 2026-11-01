# -*- coding: utf-8 -*-
''' test_ecology.py -- quanti_init_freq through MetaIBM's own life cycle (metacommunity -> patch -> habitat -> indfactory.make_founder)

    Run (conda env metaibm-nemo-test3):
        PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test3/bin/python test_ecology.py [group ...]

    Landscape (example2-like, smaller): a mainland of one patch with four 10 x 10 habitats (one species each, global_index 4-7)
    and four one-habitat 10 x 10 islands (global_index 0-3). Species, all two traits, mutation 0, Ve 0:
        sp1 diallelic              loci (3, 2)  quanti_init_freq [0.2, 0.4, 0.6, 0.8, 0.5]   ([p_1..p_L] form)
        sp2 diallelic_asymmetrical loci (2, 2)  quanti_init_freq 0.5                         (scalar form)
        sp3 diallelic              loci (3, 2)  quanti_init_freq [[0.0]]                     (monomorphic: every allele is allele 1)
        sp4 continuous (example-like, quanti_init_model 2), no quanti_init_freq              (control)
    Helpers (Ctx, head_copy) come from test_quanti_init_freq.py in this folder.
'''
import sys
sys.dont_write_bytecode = True
import os, json, time, shutil, hashlib, tempfile, subprocess, random
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from test_quanti_init_freq import Ctx, head_copy

GROUPS = ('init_stream', 'cycle_asexual', 'cycle_sexual', 'cycle_mix', 'mix_init_note', 'repro', 'regression_cycle')
SPECIES = {   # species_id: (factory kwargs, quanti_init_freq or None, founder mean phenotype = the habitat environment)
    'sp1': (dict(model='diallelic', loci=(3, 2)), [0.2, 0.4, 0.6, 0.8, 0.5], [-0.12, 0.12]),
    'sp2': (dict(model='diallelic_asymmetrical', loci=(2, 2)), 0.5, [0.2, 0.2]),
    'sp3': (dict(model='diallelic', loci=(3, 2)), [[0.0]], [-0.6, -0.4]),
    'sp4': (dict(model='continuous', loci=(20, 20), mut=0.0, extra=[('quanti_init_model', 2), ('quanti_init_value', '{{0.5, 0.5}}'), ('quanti_init_variance', '{{0.01, 0.01}}')]), None, [0.5, 0.5]),
}

##****************************** child side ******************************
def build(c, use_freq=True, seed=12):
    ''' factories + mainland + islands; use_freq=False builds every species without metaibm_nemo_params (HEAD-compatible) '''
    import metaibm
    facs = {}
    for sp, (kw, freq, _) in SPECIES.items():
        f = c.factory(freq=freq, seed=seed, **kw) if (use_freq and freq is not None) else c.factory(seed=seed, **kw)
        f.species_id = sp                                                    # Ctx names factories sp<n>; the landscape maps species by 'sp<k+1>'
        facs[sp] = f
    envs = [SPECIES[sp][2] for sp in SPECIES]
    mainland = metaibm.metacommunity(metacommunity_name='mainland')
    patch = metaibm.patch(patch_name='mainland', patch_index=0, location=(0, 0))
    for k, env in enumerate(envs):
        patch.add_habitat(hab_name='h%d' % (k + 1), hab_index=k, hab_location=(k // 2, k % 2), num_env_types=2, env_types_name=('env1', 'env2'),
                          mean_env_ls=list(env), var_env_ls=[0.025 ** 2] * 2, length=10, width=10, global_index=4 + k)
    mainland.add_patch(patch_name='mainland', patch_object=patch)
    islands = metaibm.metacommunity(metacommunity_name='islands')
    for k, env in enumerate(envs):
        p = metaibm.patch(patch_name='patch%d' % (k + 1), patch_index=k, location=(k // 2, k % 2))
        p.add_habitat(hab_name='h%d' % (k + 1), hab_index=0, hab_location=(k // 2, k % 2), num_env_types=2, env_types_name=('env1', 'env2'),
                      mean_env_ls=list(env), var_env_ls=[0.025 ** 2] * 2, length=10, width=10, global_index=k)
        islands.add_patch(patch_name='patch%d' % (k + 1), patch_object=p)
    return facs, envs, mainland, islands

def individuals(meta):
    for pid, p in meta.set.items():
        for hid, h in p.set.items():
            for (r, cc) in sorted(h.occupied_site_pos_ls):
                yield (pid, hid, r, cc), h.set['microsite_individuals'][r][cc]

def invariants(c, metas, sp3_all1=True):
    ''' every diallelic individual: phenotype == sum of allele values (Ve 0); sp3 stays all allele 1 (when it was redrawn with p=0); species <-> factory '''
    np = c.np
    bad, n, freq = [], 0, {'sp1': [], 'sp2': []}
    for meta in metas:
        for where, ind in individuals(meta):
            n += 1
            if ind.indfactory.species_id != ind.species_id: bad.append('%s species/factory mismatch' % (where,))
            if ind.species_id == 'sp4': continue
            if not np.allclose(ind.get_indi_phenotype_ls(), c.genotype(ind), atol=1e-12): bad.append('%s phenotype != genotype' % (where,))
            b = c.bits(ind)
            if sp3_all1 and ind.species_id == 'sp3' and not (b == 1).all(): bad.append('%s sp3 bits %s' % (where, b.tolist()))
            if ind.species_id in freq: freq[ind.species_id].append(1 - b)
    return n, bad, {sp: (np.array(v).mean(axis=(0, 1)).round(3).tolist() if v else None) for sp, v in freq.items()}

def digest(c, metas):
    rows = []
    for meta in metas:
        for where, ind in individuals(meta):
            rows.append([list(map(str, where)), ind.species_id, ind.individual_id, c.values(ind).tolist(), [repr(x) for x in ind.get_indi_phenotype_ls()]])
    return hashlib.md5(json.dumps(rows).encode()).hexdigest()

def g_init_stream(c):
    ''' meta_initialize -> patch_initialize -> hab_initialize -> make_founder: every founder equals the reference np stream walked in the same order '''
    np = c.np
    facs, envs, mainland, _ = build(c)
    np.random.seed(2024); random.seed(2024)
    mainland.meta_initialize(facs, envs, reproduce_mode='sexual')
    np.random.seed(2024)
    ok, n, msgs = True, 0, []
    for pid, p in mainland.set.items():
        for hid, h in p.set.items():
            sp = 'sp%d' % (envs.index(h.mean_env_ls) + 1)
            kw, freq, _ = SPECIES[sp]
            for r in range(h.length):
                for cc in range(h.width):
                    ind = h.set['microsite_individuals'][r][cc]; n += 1
                    if ind.species_id != sp: ok = False; msgs.append('%s %s species %s' % (hid, (r, cc), ind.species_id))
                    if freq is None: continue
                    L = sum(kw['loci'])
                    want = (np.random.random((2, L)) >= np.broadcast_to(np.asarray(freq, dtype=float).reshape(-1), (L,))).astype(int)
                    if not (c.bits(ind) == want).all(): ok = False; msgs.append('%s %s bits' % (hid, (r, cc)))
    c.check('ecology/init/stream', ok and n == 400, '400 founders via meta_initialize: species by habitat env, sp1-sp3 bits == reference np stream in patch/habitat/row/col order' if ok else '; '.join(msgs[:5]))
    n_ind, bad, fr = invariants(c, [mainland])
    c.check('ecology/init/invariants', not bad, '%d founders: phenotype == genotype, sp3 all allele 1; allele-0 freq sp1=%s sp2=%s' % (n_ind, fr['sp1'], fr['sp2']) if not bad else '; '.join(bad[:5]))
    genders = {sp: sorted({ind.gender for _, ind in individuals(mainland) if ind.species_id == sp}) for sp in SPECIES}
    c.check('ecology/init/genders', all(g == ['female', 'male'] for g in genders.values()), 'sexual init: both genders redrawn alike %s' % genders)
    c.check('ecology/init/sp4-untouched', facs['sp4']._quanti_seq_pos is None and all(facs[s]._quanti_seq_pos is not None for s in ('sp1', 'sp2', 'sp3')), 'control species never enters the redraw')
    mainland.meta_clear_up_microsite_individuals()
    c.check('ecology/init/recycled', all(f.minted_num == f.recycled_num for f in facs.values()), 'after meta_clear_up: minted == recycled per species %s' % {s: (f.minted_num, f.recycled_num) for s, f in facs.items()})

def cycle(c, mode, steps=30, use_freq=True):
    ''' example2's life cycle on a smaller landscape; checks invariants every step; RETURN: digest of the final landscape '''
    np = c.np
    np.random.seed(7); random.seed(7)                                       # before build(): add_habitat draws the microsite environments from np.random
    facs, envs, mainland, islands = build(c, use_freq=use_freq)
    mainland.meta_initialize(facs, envs, reproduce_mode='asexual' if mode == 'asexual' else 'sexual')
    n0, bad0, fr0 = invariants(c, [mainland], use_freq)
    all_bad, traj = list(bad0), [(0, n0, fr0)]
    for t in range(1, steps + 1):
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
        n, bad, fr = invariants(c, [mainland, islands], use_freq)
        all_bad += ['t=%d %s' % (t, b) for b in bad]
        if t % 10 == 0: traj.append((t, n, fr))
    d = digest(c, [mainland, islands])
    for meta in (mainland, islands):
        meta.meta_clear_up_microsite_individuals(); meta.meta_clear_up_offspring_and_immigrant_pool(); meta.meta_clear_up_dormancy_pool()
    balanced = {s: (f.minted_num, f.recycled_num) for s, f in facs.items()}
    for f in facs.values(): f.purge_recycling_pool()
    return d, all_bad, traj, balanced

def _g_cycle(c, mode):
    d, bad, traj, balanced = cycle(c, mode)
    c.check('ecology/%s/invariants' % mode, not bad, '30 steps, mainland + islands, every step: phenotype == genotype, sp3 all allele 1, species <-> factory; (t, alive, allele-0 freq) %s' % traj if not bad else '; '.join(bad[:5]))
    c.check('ecology/%s/minted==recycled' % mode, all(m == r for m, r in balanced.values()), str(balanced))
    c.data('digest/%s' % mode, d)

def g_cycle_asexual(c): _g_cycle(c, 'asexual')
def g_cycle_sexual(c): _g_cycle(c, 'sexual')
def g_cycle_mix(c): _g_cycle(c, 'mix')

def g_repro(c):
    d, bad, _, _ = cycle(c, 'sexual', steps=15)
    c.data('digest/repro', d)

def g_mix_init_note(c):
    ''' outside this change: hab_initialize has no 'mix' branch (founder gender undefined) -- recorded once '''
    facs, envs, mainland, _ = build(c)
    try:
        mainland.meta_initialize(facs, envs, reproduce_mode='mix')
        c.check('ecology/note/meta_initialize-mix', True, 'reproduce_mode="mix" accepted (%d founders)' % mainland.get_meta_individual_num(), note=True)
    except Exception as e:
        c.check('ecology/note/meta_initialize-mix', True, 'pre-existing, not part of quanti_init_freq: meta_initialize(reproduce_mode="mix") raises %s: %s (habitat.hab_initialize sets gender only for sexual/asexual); the cycle_mix group initializes with "sexual"' % (type(e).__name__, e), note=True)

def g_regression_cycle(c):
    ''' no species uses quanti_init_freq: the whole sexual life cycle is identical on HEAD and on the working tree (parent compares) '''
    d, bad, _, balanced = cycle(c, 'sexual', steps=20, use_freq=False)
    c.data('digest/regression', d)
    c.data('regression/bad', len(bad))

##****************************** parent side ******************************
def run_child(group, root=PROJECT_ROOT, py=sys.executable):
    t0 = time.time()
    r = subprocess.run([py, '-I', '-B', str(Path(__file__).resolve()), '--child', group, str(root)], capture_output=True, text=True, cwd=tempfile.gettempdir())
    raw = r.stdout + r.stderr
    checks = [json.loads(l[3:]) for l in raw.splitlines() if l.startswith('@@ ')]
    data = {d['key']: d['value'] for d in (json.loads(l[7:]) for l in raw.splitlines() if l.startswith('@@DATA '))}
    if r.returncode != 0: checks.append({'group': group, 'id': '%s/child-exit' % group, 'status': 'FAIL', 'detail': 'exit %d: %s' % (r.returncode, raw[-800:])})
    return checks, data, time.time() - t0

def main(groups):
    md5 = hashlib.md5((PROJECT_ROOT / 'metaibm' / 'indfactory.py').read_bytes()).hexdigest()
    print('python %s | indfactory.py md5 %s' % (sys.executable, md5))
    results, t_all = [], time.time()
    for g in groups:
        if g == 'repro':
            _, d1, t1 = run_child('repro'); _, d2, t2 = run_child('repro')
            a, b = d1.get('digest/repro'), d2.get('digest/repro')
            results.append({'group': g, 'id': 'ecology/repro/cross-process', 'status': 'PASS' if a and a == b else 'FAIL', 'detail': 'sexual life cycle, 15 steps, same seeds, two processes: md5 %s / %s' % (str(a)[:12], str(b)[:12])})
            dt = t1 + t2
        elif g == 'regression_cycle':
            hd = head_copy()
            try:
                ch, dh, t1 = run_child(g, root=hd); cw, dw, t2 = run_child(g); results += ch + cw
                a, b = dh.get('digest/regression'), dw.get('digest/regression')
                results.append({'group': g, 'id': 'ecology/regression/HEAD==tree', 'status': 'PASS' if a and a == b and dw.get('regression/bad') == 0 else 'FAIL',
                                'detail': 'no species with quanti_init_freq, sexual life cycle 20 steps (mainland + islands): HEAD md5 %s | tree md5 %s' % (str(a)[:12], str(b)[:12])})
            finally: shutil.rmtree(hd, ignore_errors=True)
            dt = t1 + t2
        else:
            checks, data, dt = run_child(g); results += checks
        print('  %-18s %6.1fs' % (g, dt), flush=True)
    n = {s: sum(1 for x in results if x['status'] == s) for s in ('PASS', 'FAIL', 'NOTE')}
    print('\n| status | id | detail |\n|---|---|---|')
    for x in results: print('| %s | `%s` | %s |' % (x['status'], x['id'], x['detail'].replace('|', '\\|')))
    print('\nPASS %d  FAIL %d  NOTE %d  (%.0fs)  indfactory.py md5 %s' % (n['PASS'], n['FAIL'], n['NOTE'], time.time() - t_all, md5))
    return 1 if n['FAIL'] else 0

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--child':
        ctx = Ctx(sys.argv[2], Path(sys.argv[3]))
        try: globals()['g_' + sys.argv[2]](ctx)
        finally: ctx.close()
    else:
        sys.exit(main([g for g in sys.argv[1:]] or list(GROUPS)))
