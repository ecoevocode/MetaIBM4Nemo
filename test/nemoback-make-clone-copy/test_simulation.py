# Scaled-down full simulations of experiments/model.py, both branches (reproduce_mode 'asexual'; 'sexual' = mixed births on the mainlands + mix markers on the islands).
# The time loop is copied from experiments/model.py; landscape + species from test/metaibm-nemo-run/simulate.py (importlib): model_generating_empty_metacommunity,
# GlobalIndex, make_indfactories. Scale: islands 4x4 patches x 1 habitat 8x8 (1024 sites), two mainlands 20x20, 150 steps, environment +0.1 every 30 steps (t = 29, 59, 89, 119);
# rain scaled by island area (10 x 1024 / 102400 = 0.1 per mainland per step); per-capita rates unchanged.
#   modes  v3sem    NemoBackend, the three calls whose v3 default was is_remove=False (2 colonizations, marker dispersal, marker germinate-and-birth) get is_remove=False
#          default  NemoBackend, no is_remove passed (all True) -- control
#          novar    NemoBackend as v3sem but quanti_init_variance = 0 (founders without standing genetic variance, as v3's all-0 / all-1 binary genotypes) -- control for the v3 comparison
#          v3gen    as novar + quanti_mutation_variance = 0.025^2 (a v3 bit flip moves mean(genotype) by 1/40 = 0.025) -- second control: Nemo genetics tuned to v3's
#          v3pure   pure v3 (git archive baace34 metaibm extension -> /tmp/v3pure, individual.py from baace34^ = 67a742f: baace34 itself already holds the Nemo-backed one),
#                   v3 API meta_initialize(traits_num, pheno_names_ls, pheno_var_ls, geno_len_ls, reproduce_mode, species_2_phenotype_ls), same loop, is_remove=False
#   checks  mc.check_world every 10 steps, teardown minted == recycled, digest identical for the same seed in two processes (3 seeds per branch);
#           v3sem (5 seeds) vs v3pure (5 seeds): island occupancy (time mean), sp2 share, final phenotype mean per species -> |diff| > 3 SE is a FINDING, not a FAIL;
#           per-step wall time + RSS recorded
# usage: python test_simulation.py [--steps N]                    (main: all workers, at most 6 at once)
#        python test_simulation.py --worker {nemo|v3pure} BRANCH MODE SEED REP OUT.json
import os, sys, json, time, math, random, hashlib, subprocess, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
V3PURE = '/tmp/v3pure'
SIM_DIR = '/tmp/nemoback-make-clone-copy/sim'                    # raw per-run records (outside the repo); the summary goes to C.RESULTS/test_simulation.json
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
STEPS = int(sys.argv[sys.argv.index('--steps') + 1]) if '--steps' in sys.argv else 150
SEEDS_STATS, SEEDS_REPRO, SEEDS_DEFAULT = (0, 1, 2, 3, 4), (0, 1, 2), (0, 1, 2)
MAX_PROCS = 6

# ---------------------------------------------------------------- experiments/model.py parameters, scaled
PATCH_NUM, HAB_NUM, META_LEN, ML_LEN = 16, 1, 32, 20
CHANGE_EVERY = 30
P = dict(environment_types_num=1, environment_types_name=['environment'], environment_variation_ls=[0.025], env_name_ls=['environment'], delta_var_ls=[0],
         species_2_phenotype_ls=[[0.0], [1.0]], traits_num=1, pheno_names_ls=['phenotype'], pheno_var_ls=[0.025], geno_len_ls=[20],
         dispersal_amomg_rate=0.0001, propagules_rain_num=10 * (META_LEN * META_LEN) / (320 * 320), base_dead_rate=0.1, fitness_wid=0.5, asexual_birth_rate=0.5, sexual_birth_rate=1,
         patch_dist_rate=0.00001, mutation_rate=0.0001, environment_mean_value=0.0, delta_mean_ls=[0.1])


def rss_mb():
    return int(subprocess.run(['ps', '-o', 'rss=', '-p', str(os.getpid())], capture_output=True, text=True).stdout) / 1024


def island_stats(meta):
    sp = {}; ph = {}
    for p in meta.set.values():
        for h in p.set.values():
            for row in h.set['microsite_individuals']:
                for x in row:
                    if x is None: continue
                    sp[x.species_id] = sp.get(x.species_id, 0) + 1
                    ph.setdefault(x.species_id, []).append(x.get_indi_phenotype_ls()[0])
    return sp, ph


def time_loop(reproduce_mode, mainland1, mainland2, meta_object, rm, steps, on_step, bargs=()):
    ''' the loop of experiments/model.py main() (logging / csv / plots dropped); rm = {'is_remove': False} or {} for the three calls whose v3 default was False;
        bargs = (mutation_rate, pheno_var_ls) for pure v3, () for NemoBackend (removed from its birth signatures) '''
    time_step = -1
    all_time_step = steps
    start_change_time, end_change_time = time_step + CHANGE_EVERY, all_time_step - CHANGE_EVERY
    change_time_step = CHANGE_EVERY
    g = P
    for time_step in range(all_time_step):
        t0 = time.perf_counter()
        if reproduce_mode == 'asexual':
            mainland1.meta_dead_selection(g['base_dead_rate'], g['fitness_wid'], method='niche_gaussian')
            mainland1.meta_mainland_asexual_birth_mutate_germinate(g['asexual_birth_rate'], *bargs)
            mainland2.meta_dead_selection(g['base_dead_rate'], g['fitness_wid'], method='niche_gaussian')
            mainland2.meta_mainland_asexual_birth_mutate_germinate(g['asexual_birth_rate'], *bargs)
            meta_object.meta_dead_selection(g['base_dead_rate'], g['fitness_wid'], method='niche_gaussian')
            meta_object.meta_asex_reproduce_calculation_into_offspring_marker_pool(g['asexual_birth_rate'])
            if (time_step + 1) % change_time_step == 0 and start_change_time <= time_step <= end_change_time: meta_object.meta_offset_environmental_values(g['env_name_ls'], g['delta_mean_ls'], g['delta_var_ls'])
            meta_object.meta_colonize_from_propagules_rains(mainland1, g['propagules_rain_num'], **rm)
            meta_object.meta_colonize_from_propagules_rains(mainland2, g['propagules_rain_num'], **rm)
            meta_object.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool(g['dispersal_amomg_rate'], **rm)
            meta_object.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(*bargs, **rm)
            meta_object.meta_disturbance_process_in_patches(g['patch_dist_rate'])
            meta_object.meta_clear_up_offspring_marker_and_immigrant_marker_pool()
        elif reproduce_mode == 'sexual':
            mainland1.meta_dead_selection(g['base_dead_rate'], g['fitness_wid'], method='niche_gaussian')
            mainland1.meta_mainland_mixed_birth_mutate_germinate(g['asexual_birth_rate'], g['sexual_birth_rate'], *bargs)
            mainland2.meta_dead_selection(g['base_dead_rate'], g['fitness_wid'], method='niche_gaussian')
            mainland2.meta_mainland_mixed_birth_mutate_germinate(g['asexual_birth_rate'], g['sexual_birth_rate'], *bargs)
            meta_object.meta_dead_selection(g['base_dead_rate'], g['fitness_wid'], method='niche_gaussian')
            meta_object.meta_mix_reproduce_calculation_with_offspring_marker_pool(g['asexual_birth_rate'], g['sexual_birth_rate'])
            if (time_step + 1) % change_time_step == 0 and start_change_time <= time_step <= end_change_time: meta_object.meta_offset_environmental_values(g['env_name_ls'], g['delta_mean_ls'], g['delta_var_ls'])
            meta_object.meta_colonize_from_propagules_rains(mainland1, g['propagules_rain_num'], **rm)
            meta_object.meta_colonize_from_propagules_rains(mainland2, g['propagules_rain_num'], **rm)
            meta_object.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool(g['dispersal_amomg_rate'], **rm)
            meta_object.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(*bargs, **rm)
            meta_object.meta_disturbance_process_in_patches(g['patch_dist_rate'])
            meta_object.meta_clear_up_offspring_marker_and_immigrant_marker_pool()
        on_step(time_step, time.perf_counter() - t0)


def landscape(gen, gidx_factory):
    patch_num_x_axis = int(math.sqrt(PATCH_NUM))
    patch_location_ls = [(i, j) for i in range(patch_num_x_axis) for j in range(patch_num_x_axis)]
    hab_len = int(META_LEN / patch_num_x_axis / int(math.sqrt(HAB_NUM)))
    meta = gen('metacommunity', PATCH_NUM, patch_location_ls, HAB_NUM, hab_len, hab_len, 0, P['environment_mean_value'], P['environment_types_num'], P['environment_types_name'], P['environment_variation_ls'], gidx_factory(100))
    ml1 = gen('mainland1', 1, [(0, 0)], 1, ML_LEN, ML_LEN, 0, 0.0, 1, ['environment'], [0.025], gidx_factory(0))
    ml2 = gen('mainland2', 1, [(0, 0)], 1, ML_LEN, ML_LEN, 0, 1.0, 1, ['environment'], [0.025], gidx_factory(1))
    return meta, ml1, ml2


def run(kind, branch, mode, seed, rep, out_path):
    rm = {} if mode == 'default' else {'is_remove': False}
    rec = {'kind': kind, 'branch': branch, 'mode': mode, 'seed': seed, 'rep': rep, 'steps': STEPS, 'step_s': [], 'rss_mb': [], 'occupancy': [], 'sp': [], 'mainland_n': [], 'checks': [], 'error': None}
    import numpy as np
    if kind == 'nemo':
        sys.path.insert(0, HERE)
        import test_core as C
        spec = importlib.util.spec_from_file_location('simulate_run', os.path.join(REPO, 'test', 'metaibm-nemo-run', 'simulate.py'))
        sim = importlib.util.module_from_spec(spec); spec.loader.exec_module(sim)
        gen = lambda *a: sim.model_generating_empty_metacommunity(*a[:-1], gidx=a[-1])[0]
        random.seed(seed); np.random.seed(seed)
        meta, ml1, ml2 = landscape(gen, sim.GlobalIndex)
        if mode in ('novar', 'v3gen'):
            F = {'sp%d' % (k + 1): C.indfactory('sp%d' % (k + 1), list(P['pheno_names_ls']), sim.nemo_params('sp%d' % (k + 1), ph, loci=P['geno_len_ls'][0], mutation_rate=P['mutation_rate'], init_variance=0.0,
                                                                                    mutation_variance=(1 / (2 * P['geno_len_ls'][0])) ** 2 if mode == 'v3gen' else 0.01,
                                                                                    env_variance=P['pheno_var_ls'][0] ** 2), seed=1000 * seed + 1 + k) for k, ph in enumerate(P['species_2_phenotype_ls'])}
        else:
            F = sim.make_indfactories(P['species_2_phenotype_ls'], P['pheno_names_ls'], seed=1000 * seed + 1, loci=P['geno_len_ls'][0], mutation_rate=P['mutation_rate'], pheno_var=P['pheno_var_ls'][0])
        fs = list(F.values())
        calls = rec['calls'] = {'make_copy': 0, 'make_clone': 0, 'make_offspring': 0}
        def counting(cls, name):
            orig = getattr(cls, name)
            def wrapped(*a, **k):
                calls[name] += 1; return orig(*a, **k)
            setattr(cls, name, wrapped)
        counting(C.individual, 'make_copy'); counting(C.indfactory, 'make_clone'); counting(C.indfactory, 'make_offspring')
        short = rec['mixed_birth_short'] = {'calls': 0, 'filled': 0, 'site left empty although asex + sex >= empty': 0, 'too few parents (asex + sex < empty)': 0}
        orig_mixed = C.habitat.hab_mixed_reproduce_germinate
        def mixed(h, *a, **k):
            empty = len(h.empty_site_pos_ls)
            pairs = sum(min(len(d['female']), len(d['male'])) for d in h.species_category_for_sexual_parents_pos.values() if 'female' in d and 'male' in d)
            want = int(np.around(len(h.asexual_parent_pos_ls) * a[0])) + int(np.around(pairs * a[1]))
            born = orig_mixed(h, *a, **k)
            short['calls'] += 1
            if born == empty: short['filled'] += 1
            elif want >= empty: short['site left empty although asex + sex >= empty'] += 1
            else: short['too few parents (asex + sex < empty)'] += 1
            return born
        C.habitat.hab_mixed_reproduce_germinate = mixed
        ml1.meta_initialize(F, P['species_2_phenotype_ls'], branch); ml2.meta_initialize(F, P['species_2_phenotype_ls'], branch)
    else:
        sys.path.insert(0, V3PURE)
        import metaibm
        assert os.path.dirname(metaibm.__file__) == os.path.join(V3PURE, 'metaibm'), metaibm.__file__
        from metaibm.metacommunity import metacommunity
        from metaibm.patch import patch
        def gen(meta_name, patch_num, patch_location_ls, hab_num, hab_length, hab_width, dormancy_pool_max_size, environment_mean_value, environment_types_num, environment_types_name, environment_variation_ls, _gidx):
            ''' experiments/model.py generating_empty_metacommunity (v3 add_habitat: no global_index) '''
            meta_object = metacommunity(metacommunity_name=meta_name)
            hab_num_x_axis, hab_num_y_axis = int(np.sqrt(hab_num)), int(np.sqrt(hab_num))
            for i in range(patch_num):
                p = patch('patch%d' % (i + 1), i, patch_location_ls[i])
                patch_x_loc, patch_y_loc = patch_location_ls[i]
                for j in range(hab_num):
                    p.add_habitat(hab_name='h%s' % str(j + 1), hab_index=j, hab_location=(patch_x_loc * hab_num_x_axis + j // hab_num_y_axis, patch_y_loc * hab_num_y_axis + j % hab_num_y_axis),
                                  num_env_types=environment_types_num, env_types_name=environment_types_name, mean_env_ls=[environment_mean_value], var_env_ls=environment_variation_ls,
                                  length=hab_length, width=hab_width, dormancy_pool_max_size=dormancy_pool_max_size)
                meta_object.add_patch(patch_name=p.name, patch_object=p)
            return meta_object
        random.seed(seed); np.random.seed(seed)
        meta, ml1, ml2 = landscape(gen, lambda start: None)
        for ml in (ml1, ml2): ml.meta_initialize(P['traits_num'], P['pheno_names_ls'], P['pheno_var_ls'], P['geno_len_ls'], branch, P['species_2_phenotype_ls'])
    rec['rss_mb'].append(round(rss_mb(), 1))
    sd = lambda v: math.sqrt(sum((x - sum(v) / len(v)) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else None
    rec['founder_pheno_sd'] = {'mainland1 sp1': sd(island_stats(ml1)[1].get('sp1', [])), 'mainland2 sp2': sd(island_stats(ml2)[1].get('sp2', []))}
    rec['island_pheno'], rec['island_env'] = [], []

    def on_step(t, dt):
        rec['step_s'].append(round(dt, 5)); rec['rss_mb'].append(round(rss_mb(), 1))
        sp, ph = island_stats(meta)
        rec['occupancy'].append(meta.get_meta_individual_num()); rec['sp'].append([sp.get('sp1', 0), sp.get('sp2', 0)])
        rec['island_pheno'].append({k: sum(v) / len(v) for k, v in ph.items()})
        rec['island_env'].append(sum(h.mean_env_ls[0] for p in meta.set.values() for h in p.set.values()) / PATCH_NUM / HAB_NUM)
        rec['mainland_n'].append([ml1.get_meta_individual_num(), ml2.get_meta_individual_num()])
        if kind == 'nemo' and ((t + 1) % 10 == 0 or t == 0):
            try: C.mc.check_world([ml1, ml2, meta], fs, 't=%d' % t); rec['checks'].append([t, True, None])
            except AssertionError as e: rec['checks'].append([t, False, str(e)])
    t0 = time.time()
    try: time_loop(branch, ml1, ml2, meta, rm, STEPS, on_step, bargs=(P['mutation_rate'], P['pheno_var_ls']) if kind == 'v3pure' else ())
    except Exception as e:
        import traceback
        rec['error'] = '%s: %s | %s' % (type(e).__name__, e, traceback.format_exc()[-1500:])
    rec['seconds'] = round(time.time() - t0, 2)
    sp, ph = island_stats(meta)
    rec['final_sp'] = sp; rec['final_pheno_mean'] = {k: sum(v) / len(v) for k, v in ph.items()}
    rec['final_mainland_pheno_mean'] = {m.metacommunity_name: (lambda d: {k: sum(v) / len(v) for k, v in d.items()})(island_stats(m)[1]) for m in (ml1, ml2)}
    if kind == 'nemo':
        rec['digest'] = C.mc.digest([ml1, ml2, meta])
        rec['minted_before_teardown'] = {f.species_id: [f.minted_num, f.recycled_num] for f in fs}
        try: C.mc.teardown([ml1, ml2, meta], fs); rec['teardown'] = {f.species_id: [f.minted_num, f.recycled_num] for f in fs}; rec['teardown_ok'] = all(f.minted_num == f.recycled_num for f in fs)
        except AssertionError as e: rec['teardown'] = str(e); rec['teardown_ok'] = False
    with open(out_path, 'w') as fh: json.dump(rec, fh)
    print('JSON ' + json.dumps({'out': out_path, 'error': rec['error']}), flush=True)


# ---------------------------------------------------------------- main
def ensure_v3pure():
    ''' /tmp/v3pure = git archive baace34 metaibm extension + metaibm/individual.py from baace34^ (pure-Python individual) '''
    marker = os.path.join(V3PURE, '.made_by_nemoback_make_clone_copy')
    if os.path.exists(marker): return open(marker).read()
    os.makedirs(V3PURE, exist_ok=True)
    arc = subprocess.run(['git', 'archive', 'baace34', 'metaibm', 'extension'], capture_output=True, cwd=REPO, check=True).stdout
    subprocess.run(['tar', '-x', '-C', V3PURE], input=arc, check=True)
    ind = subprocess.run(['git', 'show', 'baace34^:metaibm/individual.py'], capture_output=True, cwd=REPO, check=True).stdout
    open(os.path.join(V3PURE, 'metaibm', 'individual.py'), 'wb').write(ind)
    os.remove(os.path.join(V3PURE, 'metaibm', 'indfactory.py'))
    note = 'git archive baace34 metaibm extension; metaibm/individual.py = baace34^ (%s); indfactory.py removed' % subprocess.run(['git', 'rev-parse', '--short', 'baace34^'], capture_output=True, text=True, cwd=REPO).stdout.strip()
    open(marker, 'w').write(note)
    return note


def launch(jobs):
    from concurrent.futures import ThreadPoolExecutor
    os.makedirs(SIM_DIR, exist_ok=True)
    def one(job):
        kind, branch, mode, seed, rep = job
        out = os.path.join(SIM_DIR, '%s_%s_%s_s%d_r%d.json' % (kind, branch, mode, seed, rep))
        r = subprocess.run([sys.executable, os.path.abspath(__file__), '--worker', kind, branch, mode, str(seed), str(rep), out] + (['--steps', str(STEPS)] if '--steps' in sys.argv else []),
                           capture_output=True, text=True, env=ENV, cwd=HERE)
        if r.returncode != 0 or not os.path.exists(out): return job, {'error': 'worker rc=%d: %s' % (r.returncode, (r.stderr or r.stdout)[-1500:])}
        return job, json.load(open(out))
    with ThreadPoolExecutor(MAX_PROCS) as ex: return dict(ex.map(one, jobs))


def mean_se(xs):
    xs = [x for x in xs if x is not None]
    if len(xs) < 2: return (xs[0] if xs else None), None, len(xs)
    m = sum(xs) / len(xs); sd = math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))
    return m, sd / math.sqrt(len(xs)), len(xs)


def metrics(r):
    occ = r['occupancy']; sp = r['sp']
    share = [b / (a + b) for a, b in sp if a + b > 0]
    return {'occupancy (time mean, of %d sites)' % (META_LEN * META_LEN): sum(occ) / len(occ) if occ else None,
            'occupancy (last 30 steps)': sum(occ[-30:]) / len(occ[-30:]) if occ else None,
            'sp2 share on islands (time mean)': sum(share) / len(share) if share else None,
            'final phenotype mean sp1 (islands)': r['final_pheno_mean'].get('sp1'), 'final phenotype mean sp2 (islands)': r['final_pheno_mean'].get('sp2'),
            'final phenotype mean sp1 (mainland1)': r['final_mainland_pheno_mean']['mainland1'].get('sp1'), 'final phenotype mean sp2 (mainland2)': r['final_mainland_pheno_mean']['mainland2'].get('sp2'),
            'mainland1 size (time mean)': sum(m[0] for m in r['mainland_n']) / len(r['mainland_n']), 'mainland2 size (time mean)': sum(m[1] for m in r['mainland_n']) / len(r['mainland_n'])}


def main():
    sys.path.insert(0, HERE)
    import test_core as C
    note = ensure_v3pure()
    C.section('setup'); print('v3pure:', note)
    jobs = []
    for b in ('asexual', 'sexual'):
        jobs += [('nemo', b, 'v3sem', s, 0) for s in SEEDS_STATS] + [('nemo', b, 'v3sem', s, 1) for s in SEEDS_REPRO]
        jobs += [('nemo', b, 'default', s, 0) for s in SEEDS_DEFAULT] + [('v3pure', b, 'v3sem', s, 0) for s in SEEDS_STATS] + [('nemo', b, m, s, 0) for m in ('novar', 'v3gen') for s in SEEDS_STATS]
    t0 = time.time(); R = launch(jobs); wall = time.time() - t0
    print('%d runs in %.0f s (%d at once)' % (len(jobs), wall, MAX_PROCS))
    perf, comparison, tvf = {}, {}, {}
    for b in ('asexual', 'sexual'):
        C.section('branch %s' % b)
        nemo = [R[j] for j in jobs if j[0] == 'nemo' and j[1] == b]
        errs = [(j, R[j]['error']) for j in jobs if j[1] == b and R[j].get('error')]
        C.ok('%s: all %d runs (NemoBackend v3sem / default, pure v3) finished %d steps without an exception' % (b, sum(1 for j in jobs if j[1] == b), STEPS), not errs, errs[:3] or None)
        bad_chk = [(r['mode'], r['seed'], r['rep'], c) for r in nemo if 'checks' in r for c in r['checks'] if not c[1]]
        n_chk = sum(len(r.get('checks', [])) for r in nemo)
        C.ok('%s: mc.check_world every 10 steps in every NemoBackend run (%d checks)' % (b, n_chk), not bad_chk and n_chk > 0, bad_chk[:3] or None)
        C.ok('%s: teardown minted == recycled in every NemoBackend run' % b, all(r.get('teardown_ok') for r in nemo), {('%s s%d r%d' % (r['mode'], r['seed'], r['rep'])): r.get('teardown') for r in nemo})
        for s in SEEDS_REPRO:
            a, c = R[('nemo', b, 'v3sem', s, 0)], R[('nemo', b, 'v3sem', s, 1)]
            C.ok('%s: seed %d in two processes -> identical mc.digest at the end' % (b, s), a.get('digest') and a.get('digest') == c.get('digest'), (a.get('digest'), c.get('digest')))
        C.ok('%s: different seeds -> different digests' % b, len({R[('nemo', b, 'v3sem', s, 0)].get('digest') for s in SEEDS_REPRO}) == len(SEEDS_REPRO))
        # performance
        for kind, mode in (('nemo', 'v3sem'), ('nemo', 'default'), ('nemo', 'novar'), ('nemo', 'v3gen'), ('v3pure', 'v3sem')):
            rs = [R[j] for j in jobs if j[:3] == (kind, b, mode) and 'step_s' in R[j]]
            if not rs: continue
            steps = [x for r in rs for x in r['step_s']]
            perf['%s %s %s' % (b, kind, mode)] = {'mixed birth calls not filling every empty site (sum over runs)': (lambda d: {k: sum(r['mixed_birth_short'][k] for r in rs) for k in d})(rs[0]['mixed_birth_short']) if 'mixed_birth_short' in rs[0] else None, 'runs': len(rs), 'calls per run (mean)': {k: round(sum(r['calls'][k] for r in rs) / len(rs)) for k in rs[0]['calls']} if 'calls' in rs[0] else None, 'mean step ms': round(1000 * sum(steps) / len(steps), 2), 'max step ms': round(1000 * max(steps), 1),
                                                  'run s (mean)': round(sum(r['seconds'] for r in rs) / len(rs), 1),
                                                  'RSS MB start -> step 30 -> end (max over runs)': [max(r['rss_mb'][0] for r in rs), max(r['rss_mb'][min(30, len(r['rss_mb']) - 1)] for r in rs), max(r['rss_mb'][-1] for r in rs)],
                                                  'final islands occupancy (mean)': round(sum(r['occupancy'][-1] for r in rs) / len(rs), 1)}
        # NemoBackend (v3 semantics) vs pure v3
        mN = [metrics(R[('nemo', b, 'v3sem', s, 0)]) for s in SEEDS_STATS if not R[('nemo', b, 'v3sem', s, 0)].get('error')]
        mV = [metrics(R[('v3pure', b, 'v3sem', s, 0)]) for s in SEEDS_STATS if not R[('v3pure', b, 'v3sem', s, 0)].get('error')]
        m0 = [metrics(R[('nemo', b, 'novar', s, 0)]) for s in SEEDS_STATS if not R[('nemo', b, 'novar', s, 0)].get('error')]
        mG = [metrics(R[('nemo', b, 'v3gen', s, 0)]) for s in SEEDS_STATS if not R[('nemo', b, 'v3gen', s, 0)].get('error')]
        zf = lambda a, sa, v, sv: abs(a - v) / math.sqrt(sa ** 2 + sv ** 2) if None not in (a, v, sa, sv) and (sa or sv) else None
        rows = {}
        for k in (mN[0] if mN else {}):
            (a, sa, na), (v, sv, nv), (c, sc, nc), (g, sg, ng) = mean_se([m[k] for m in mN]), mean_se([m[k] for m in mV]), mean_se([m[k] for m in m0]), mean_se([m[k] for m in mG])
            rows[k] = {'nemo': [a, sa, na], 'v3': [v, sv, nv], 'z': zf(a, sa, v, sv), 'nemo novar': [c, sc, nc], 'z novar vs v3': zf(c, sc, v, sv), 'nemo v3gen': [g, sg, ng], 'z v3gen vs v3': zf(g, sg, v, sv)}
        comparison[b] = rows
        traj = {}
        for kind, mode in (('nemo', 'v3sem'), ('nemo', 'novar'), ('nemo', 'v3gen'), ('v3pure', 'v3sem')):
            rs = [R[(kind, b, mode, s_, 0)] for s_ in SEEDS_STATS if 'island_pheno' in R[(kind, b, mode, s_, 0)]]
            row = {}
            for t in (0, 29, 59, 89, 119, STEPS - 1):
                if t >= STEPS: continue
                v = [r['island_pheno'][t].get('sp1') for r in rs if r['island_pheno'][t].get('sp1') is not None]
                row['t=%d' % t] = {'env': round(rs[0]['island_env'][t], 2) if rs else None, 'occupancy': round(sum(r['occupancy'][t] for r in rs) / len(rs), 1) if rs else None,
                                   'sp1 phenotype': round(sum(v) / len(v), 3) if v else None, 'runs with sp1': len(v)}
            row['founder phenotype SD'] = {k: round(sum(r['founder_pheno_sd'][k] for r in rs) / len(rs), 4) for k in ('mainland1 sp1', 'mainland2 sp2')} if rs else None
            traj['%s %s' % (kind, mode)] = row
        comparison[b + ' trajectory'] = traj
        # v3 semantics vs default (True) -- informational
        mD = [metrics(R[('nemo', b, 'default', s, 0)]) for s in SEEDS_DEFAULT if not R[('nemo', b, 'default', s, 0)].get('error')]
        mN3 = [metrics(R[('nemo', b, 'v3sem', s, 0)]) for s in SEEDS_DEFAULT if not R[('nemo', b, 'v3sem', s, 0)].get('error')]
        tvf[b] = {k: {'is_remove=False (v3 semantics)': mean_se([m[k] for m in mN3])[:2], 'default True': mean_se([m[k] for m in mD])[:2]} for k in (mN3[0] if mN3 else {})}
    for b in ('asexual', 'sexual'):
        for k, r in comparison[b].items():
            if r['z'] is not None and r['z'] > 3:
                C.finding('info', 'test_simulation %s branch: %s' % (b, k), 'NemoBackend (v3 semantics, %d seeds) vs pure v3 (%d seeds) differ by %.1f SE (statistical comparison, not a FAIL); controls vs v3: init_variance 0 -> %s SE, init_variance 0 + v3-sized mutation effects -> %s SE' % (r['nemo'][2], r['v3'][2], r['z'], None if r['z novar vs v3'] is None else round(r['z novar vs v3'], 1), None if r['z v3gen vs v3'] is None else round(r['z v3gen vs v3'], 1)),
                          {'nemo mean, SE': [round(x, 4) if isinstance(x, float) else x for x in r['nemo'][:2]], 'v3 mean, SE': [round(x, 4) if isinstance(x, float) else x for x in r['v3'][:2]],
                           'nemo init_variance 0 mean, SE': [round(x, 4) if isinstance(x, float) else x for x in r['nemo novar'][:2]],
                           'nemo init_variance 0 + mutation_variance 0.025^2 mean, SE': [round(x, 4) if isinstance(x, float) else x for x in r['nemo v3gen'][:2]],
                           'trajectory (mean over seeds)': comparison[b + ' trajectory']})
    return C.summary('test_simulation', {'steps': STEPS, 'wall_s': round(wall, 1), 'v3pure': note, 'params': dict(P, PATCH_NUM=PATCH_NUM, HAB_NUM=HAB_NUM, META_LEN=META_LEN, ML_LEN=ML_LEN, CHANGE_EVERY=CHANGE_EVERY),
                                             'perf': perf, 'nemo_vs_v3': comparison, 'v3sem_vs_default': tvf, 'runs': {'%s %s %s s%d r%d' % j: {k: R[j].get(k) for k in ('seconds', 'digest', 'error', 'teardown_ok')} for j in jobs}})


if __name__ == '__main__':
    if len(sys.argv) > 7 and sys.argv[1] == '--worker':
        run(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6]), sys.argv[7])
    else:
        raise SystemExit(main())
