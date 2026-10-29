#!/usr/bin/env python3
# -*- coding: utf-8 -*-
''' Benchmark MetaIBM-Nemo (this worktree, NemoBackend individuals) vs pure MetaIBM (MetaIBM-bug-fixing checkout: imported read only, bytecode writing disabled).
Same landscape, parameters, seed and call sequence on both sides; is_remove=True passed explicitly wherever the function takes it (same model semantics on both sides:
propagules / migrants / markers move). Per case: setup (build + initialize), s / step, seconds per phase, births / deaths, RSS.
  interfaces (as the experiments: islands start empty, marker schedule, CSV rows every 10 steps into a temp dir, no plots):
      model             experiments/model.py: 1 trait, 2 mainlands 100x100, 100 islands x 1 habitat 32x32, two colonizations / step, environment offset every 100 steps
      sloss-p16         experiments/model-sloss.py: 16 patches x 4 habitats 40x40, discrete environment
      sloss-p64-hetero  experiments/model-sloss.py: 64 patches x 4 habitats 20x20, continuous environment
      sloss-p1          experiments/model-sloss.py: 1 patch x 16 habitats 80x80 (among-patch dispersal: patch_num < 2 branch)
      grfe-p16          experiments/model-sloss-GRFE.py: CSV landscape, 16 patches x 64 habitats 10x10 (patch-level among dispersal, rho=0.2)
  sweeps (islands 16 x 16 habitats 20x20 fully occupied, both schedules): one ecological (dead, wid, birth, rain, dispersal, disturbance) or genetic (loci, traits,
      mutation) parameter changed at a time from BASE -- see SWEEP
  boundaries (islands fully occupied from step 0, total 102400 sites in 16 patches, one habitat size each; marker AND individual-pool schedule):
      L10 (16 x 64 habitats 10x10), L20 (16 x 16 x 20x20), L40 (16 x 4 x 40x40), L80 (16 x 1 x 80x80); pool schedule with dormancy (10% of the habitat size)
    python bench.py --impl nemo|pure --case CASE [--steps N]     one case -> bench/<case>__<impl>.json
    python bench.py --all                                          every case, pure then nemo, one process at a time (run on an idle machine), then BENCH_REPORT.md
    python bench.py --report                                       BENCH_REPORT.md from bench/*.json '''
import os, sys, json, time, random, argparse, inspect, tempfile, shutil, subprocess, platform, datetime
sys.dont_write_bytecode = True
os.environ.setdefault('MPLBACKEND', 'Agg')
HERE = os.path.dirname(os.path.abspath(__file__))
NEMO_REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
PURE_REPO = '/Users/jianhaolin/Library/CloudStorage/OneDrive-UniversityofHelsinki/MetaIBM_workspaces/Macbook_Air/MetaIBM-bug-fixing'
EXPERIMENTS = os.path.join(NEMO_REPO, 'experiments')
BENCH = os.path.join(HERE, 'bench')
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
INTERFACES = ['model', 'sloss-p16', 'sloss-p64-hetero', 'sloss-p1', 'grfe-p16']
BOUNDARIES = ['%s-%s' % (L, s) for L in ('L10', 'L20', 'L40', 'L80') for s in ('marker', 'pool')]
BASE = dict(dead=0.1, wid=0.5, birth=1.0, rain=6.4, m_within=0.1, m_among=0.1, dist=0.00001, loci=20, traits=2, mut=0.0001)
SWEEP = {   # one parameter changed at a time from BASE, islands 16 patches x 16 habitats 20x20 fully occupied (eco = ecological, gen = genetic)
    'base': {}, 'dead0.05': dict(dead=0.05), 'dead0.3': dict(dead=0.3), 'wid0.2': dict(wid=0.2), 'wid2.0': dict(wid=2.0), 'birth0.5': dict(birth=0.5),
    'rain64': dict(rain=64.0), 'disp0.5': dict(m_within=0.5, m_among=0.5), 'dist0.01': dict(dist=0.01),
    'loci5': dict(loci=5), 'loci80': dict(loci=80), 'loci320': dict(loci=320), 'traits1': dict(traits=1), 'traits4': dict(traits=4), 'mut0.01': dict(mut=0.01)}
SWEEPS = ['sweep-%s-%s' % (s, v) for v in SWEEP for s in ('marker', 'pool')]
STEPS = {'interface': 100, 'boundary': 30, 'sweep': 20}
TABLE2 = [[i / 10, j / 10] for i in range(2, 10, 2) for j in range(2, 10, 2)]
VALS = [0.2, 0.4, 0.6, 0.8]


def species_table(n):
    ''' n quanti traits: 16 species on the first two axes (0.2..0.8), 0.5 on the others; 1 trait: 4 species '''
    return [[v] for v in VALS] if n == 1 else [row + [0.5] * (n - 2) for row in TABLE2]


def env_names(n):
    return ('x_axis_environment', 'y_axis_environment') if n == 2 else tuple('env%d' % k for k in range(n))


def rss_mb():
    return int(subprocess.run(['ps', '-o', 'rss=', '-p', str(os.getpid())], capture_output=True, text=True).stdout) / 1024


# ------------------------------------------------------------------ landscapes (impl-independent; global_index only where add_habitat takes it)
class Builder():
    def __init__(self, metacommunity, patch):
        self.M, self.P = metacommunity, patch
        self.takes_gi = 'global_index' in inspect.signature(patch.add_habitat).parameters
        self.gi = 0

    def hab(self, p, name, idx, loc, mean, L, W=None, K=0, n_env=2, env_names=('x_axis_environment', 'y_axis_environment'), var=(0.025, 0.025)):
        kw = dict(hab_name=name, hab_index=idx, hab_location=loc, num_env_types=n_env, env_types_name=env_names, mean_env_ls=list(mean), var_env_ls=list(var),
                  length=L, width=L if W is None else W, dormancy_pool_max_size=K)
        if self.takes_gi: kw['global_index'] = self.gi
        self.gi += 1
        p.add_habitat(**kw)

    def sloss_mainland(self, n=2):
        m = self.M(metacommunity_name='mainland'); p = self.P('patch0', 0, (0, 0)); tab = species_table(n)
        for j in range(16):
            mean = [VALS[j // 4], VALS[j % 4]] if n == 2 else tab[j % len(tab)]
            self.hab(p, 'h%d' % j, j, (j // 4, j % 4), mean, 20, n_env=n, env_names=env_names(n), var=[0.025] * n)
        m.add_patch(patch_name='patch0', patch_object=p); return m

    def sloss_islands(self, patch_num, hetero):
        ''' experiments/model-sloss.py generating_empty_metacommunity '''
        self.gi = max(self.gi, 100)
        m = self.M(metacommunity_name='metacommunity')
        n = int(patch_num ** 0.5); hn = 16 if patch_num == 1 else 4; h1 = int(hn ** 0.5); L = int(320 / n / h1)
        lo, hi = 0.1, 0.9
        for i, (px, py) in enumerate([(a, b) for a in range(n) for b in range(n)]):
            p = self.P('patch%d' % (i + 1), i, (px, py))
            for j in range(hn):
                hx, hy = px * h1 + j // h1, py * h1 + j % h1
                if not hetero: mean = [VALS[int(hx * 4 / (n * h1)) % 4], VALS[int(hy * 4 / (n * h1)) % 4]]
                else:
                    st = (hi - lo) / (n * h1); mean = [lo + st / 2 + st * hx, lo + st / 2 + st * hy]
                self.hab(p, 'h%d' % j, j, (hx, hy), mean, L)
            m.add_patch(patch_name=p.name, patch_object=p)
        return m

    def boundary_islands(self, L, dorm, n=2):
        ''' 16 patches, 102400 sites, habitats of LxL; habitat means cycle through the species table (so the islands can be initialized: full occupancy) '''
        self.gi = max(self.gi, 100)
        m = self.M(metacommunity_name='metacommunity'); hn = 102400 // (16 * L * L); h1 = int(round(hn ** 0.5)); k = 0
        for i, (px, py) in enumerate([(a, b) for a in range(4) for b in range(4)]):
            p = self.P('patch%d' % (i + 1), i, (px, py))
            for j in range(hn):
                tab = species_table(n)
                self.hab(p, 'h%d' % j, j, (px * h1 + j // h1, py * h1 + j % h1), tab[k % len(tab)], L, K=int(0.1 * L * L) if dorm else 0, n_env=n, env_names=env_names(n), var=[0.025] * n); k += 1
            m.add_patch(patch_name=p.name, patch_object=p)
        return m

    def model_world(self):
        ''' experiments/model.py: 100 islands x 1 habitat 32x32 (env 0.0), mainland1 (env 0.0) and mainland2 (env 1.0) 100x100 '''
        one = dict(n_env=1, env_names=['environment'], var=[0.025])
        m1 = self.M(metacommunity_name='mainland1'); p = self.P('patch1', 0, (0, 0)); self.gi = 0; self.hab(p, 'h1', 0, (0, 0), [0.0], 100, **one); m1.add_patch(patch_name='patch1', patch_object=p)
        m2 = self.M(metacommunity_name='mainland2'); p = self.P('patch1', 0, (0, 0)); self.gi = 1; self.hab(p, 'h1', 0, (0, 0), [1.0], 100, **one); m2.add_patch(patch_name='patch1', patch_object=p)
        self.gi = 100
        isl = self.M(metacommunity_name='metacommunity')
        for i, (px, py) in enumerate([(a, b) for a in range(10) for b in range(10)]):
            p = self.P('patch%d' % (i + 1), i, (px, py)); self.hab(p, 'h1', 0, (px, py), [0.0], 32, **one); isl.add_patch(patch_name=p.name, patch_object=p)
        return m1, m2, isl

    def grfe_islands(self, patch_num=16, same=True):
        ''' experiments/model-sloss-GRFE.py generate_empty_metacommunity (CSV layout, read only) '''
        import pandas as pd
        lay = pd.read_csv(os.path.join(EXPERIMENTS, 'patch_habitat_layouts.csv')); lay = lay[lay['patch_num'] == patch_num].copy()

        def env(f):
            df = pd.read_csv(os.path.join(EXPERIMENTS, f), header=[0, 1]); df.columns = pd.MultiIndex.from_tuples(df.columns)
            sub = df.loc[:, [('meta', 'x_index'), ('meta', 'y_index'), ('patch_num=%d' % patch_num, str(bool(same)))]].copy(); sub.columns = ['x_index', 'y_index', 'value']
            sub['x_index'], sub['y_index'] = sub['x_index'].astype(int), sub['y_index'].astype(int); return sub
        e = env('32x32_habitats_env1.csv').merge(env('32x32_habitats_env2.csv'), on=['x_index', 'y_index']); e.columns = ['x_index', 'y_index', 'env_axis1', 'env_axis2']
        lay = lay.merge(e, left_on=['habitat_x_location', 'habitat_y_location'], right_on=['x_index', 'y_index'], how='left')
        self.gi = max(self.gi, 100)
        m = self.M(metacommunity_name='metacommunity')
        for _, pr in lay[['patch_id', 'patch_index', 'patch_location_x', 'patch_location_y']].drop_duplicates().sort_values(['patch_index', 'patch_id']).iterrows():
            p = self.P(pr['patch_id'], int(pr['patch_index']), (float(pr['patch_location_x']), float(pr['patch_location_y'])))
            for _, hr in lay[lay['patch_id'] == pr['patch_id']].sort_values(['habitat_index', 'habitat_id']).iterrows():
                self.hab(p, str(hr['habitat_id']), int(hr['habitat_index']), (int(hr['habitat_x_location']), int(hr['habitat_y_location'])), [float(hr['env_axis1']), float(hr['env_axis2'])], 10,
                         env_names=('env_axis1', 'env_axis2'))
            m.add_patch(patch_name=p.name, patch_object=p)
        return m


# ------------------------------------------------------------------ one case
def run(impl, case, steps, save_every=10, seed=0):
    repo = NEMO_REPO if impl == 'nemo' else PURE_REPO
    sys.path.insert(0, repo)
    import numpy as np
    from metaibm.patch import patch
    from metaibm.metacommunity import metacommunity
    mc_file = sys.modules['metaibm.metacommunity'].__file__
    assert mc_file.startswith(repo), mc_file
    random.seed(seed); np.random.seed(seed)
    B = Builder(metacommunity, patch)
    kind = 'interface' if case in INTERFACES else 'boundary' if case in BOUNDARIES else 'sweep'
    steps = steps or STEPS[kind]
    P = dict(BASE, **(SWEEP[case.split('-', 2)[2]] if kind == 'sweep' else {}))
    ntraits = 1 if case == 'model' else P['traits']
    table = [[0.0], [1.0]] if case == 'model' else species_table(ntraits)
    names = ('phenotype',) if case == 'model' else ('x_axis_phenotype', 'y_axis_phenotype') if ntraits == 2 else tuple('pheno%d' % k for k in range(ntraits))
    var, geno = (0.025,) * ntraits, (P['loci'],) * ntraits
    t_setup = time.perf_counter()
    if impl == 'nemo':
        sys.path.insert(0, os.path.join(HERE, '..', 'metaibm-nemo-run'))
        import simulate as R                                   # Nemo species (nemo_params / make_indfactories) shared with the full simulations
        F = R.make_indfactories(table, names, seed=1, loci=P['loci'], mutation_rate=P['mut'])
        probe = next(iter(F.values())).make_founder('female', 0)
        loci_per_trait = [probe.get_trait('quant').get_num_loci(k) for k in range(ntraits)]
        probe.recycle()
        assert loci_per_trait == [P['loci']] * ntraits, ('Nemo loci per trait %s != geno_len %d' % (loci_per_trait, P['loci']))   # same genome size as the pure side
        init = lambda m: m.meta_initialize(F, table, 'sexual')
    else:
        F = None
        init = lambda m: m.meta_initialize(ntraits, names, var, geno, 'sexual', table)
    t_build = time.perf_counter()
    if case == 'model': ml1, ml2, isl = B.model_world(); mainlands = [ml1, ml2]
    else:
        mainlands = [B.sloss_mainland(ntraits)]
        if case.startswith('sloss'): isl = B.sloss_islands(int(case.split('-p')[1].split('-')[0]), 'hetero' in case)
        elif case.startswith('grfe'): isl = B.grfe_islands()
        elif kind == 'sweep': isl = B.boundary_islands(20, case.split('-')[1] == 'pool', ntraits)
        else: isl = B.boundary_islands(int(case[1:3]), case.endswith('pool'))
    build_s = time.perf_counter() - t_build
    t_init = time.perf_counter()
    for m in mainlands: init(m)
    if kind != 'interface': init(isl)
    init_s = time.perf_counter() - t_init
    setup_s = time.perf_counter() - t_setup
    rm = {'is_remove': True}
    phase = {}
    tmp = tempfile.mkdtemp(prefix='bench_csv_')
    cols = list(isl.columns_patch_habitat_microsites_id())
    trait_names = list(names)

    def T(name, fn, *a, **kw):
        t = time.perf_counter(); r = fn(*a, **kw); phase[name] = phase.get(name, 0.0) + time.perf_counter() - t; return r
    pool = case.endswith('pool') or (kind == 'sweep' and case.split('-')[1] == 'pool')
    among = dict(method='exponential', rho=0.2) if case.startswith('grfe') else dict(method='uniform') if case == 'model' else dict(method='exponential', rho=2)
    m_among, m_within = (0.0001, None) if case == 'model' else (0.01, 0.1) if case.startswith('grfe') else (P['m_among'], P['m_within'])
    rain = 10 if case == 'model' else P['rain']
    dead, wid, birth, mut, dist = P['dead'], P['wid'], P['birth'], P['mut'], P['dist']
    bargs = () if impl == 'nemo' else (mut, var)              # mutation_rate / pheno_var_ls: pure MetaIBM only (removed from the NemoBackend birth signatures)
    rows, t0 = [], time.perf_counter()
    for t in range(steps):
        ts = time.perf_counter()
        for k, ml in enumerate(mainlands):
            T('mainland dead selection', ml.meta_dead_selection, dead, wid, method='niche_gaussian')
            T('mainland sexual birth', ml.meta_mainland_sexual_birth_mutate_germinate, birth, *bargs)
        T('islands dead selection', isl.meta_dead_selection, dead, wid, method='niche_gaussian')
        if not pool:
            T('marker calculation', isl.meta_sex_reproduce_calculation_with_offspring_marker_pool, birth)
            if case == 'model' and (t + 1) % 100 == 0: T('environment offset', isl.meta_offset_environmental_values, ['environment'], [0.1], [0])
            if case == 'model':
                for ml in mainlands: T('colonization', isl.meta_colonize_from_propagules_rains, ml, rain, **rm)
            else:
                T('colonization', isl.pairwise_sexual_colonization_from_prpagules_rains, mainlands[0], rain, **rm)
            if m_within is not None: T('marker dispersal within', isl.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool, m_within, **rm)
            T('marker dispersal among', isl.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool, m_among, **among, **rm)
            T('germinate + birth from markers', isl.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool, *bargs, **rm)
            T('disturbance', isl.meta_disturbance_process_in_patches, dist)
            T('marker clear-up', isl.meta_clear_up_offspring_marker_and_immigrant_marker_pool)
        else:
            T('reproduce into offspring pool', isl.meta_sex_reproduce_mutate_into_offspring_pool, birth, *bargs)
            T('colonization', isl.pairwise_sexual_colonization_from_prpagules_rains, mainlands[0], rain, **rm)
            T('dispersal within -> immigrant pool', isl.meta_dispersal_within_patch_from_offspring_to_immigrant_pool, m_within, **rm)
            T('dispersal among -> immigrant pool', isl.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool, m_among, method='exponential', rho=2, **rm)
            T('germinate O+I+D', isl.meta_local_germinate_from_offspring_immigrant_and_dormancy_pool, **rm)
            T('dormancy', isl.meta_dormancy_process_from_offspring_pool_and_immigrant_pool)
            T('disturbance', isl.meta_disturbance_process_in_patches, dist)
            T('clear-up O+I', isl.meta_clear_up_offspring_and_immigrant_pool)
        if kind == 'interface' and (t + 1) % save_every == 0:
            T('data: sp_id csv.gz', lambda: isl.meta_distribution_data_all_time_to_csv_gz(isl.get_meta_microsites_individuals_sp_id_values(), os.path.join(tmp, 'sp.csv.gz'), ['t%d' % t], cols, mode='a'))
            for tr in trait_names:
                T('data: phenotype csv.gz', lambda: isl.meta_distribution_data_all_time_to_csv_gz(isl.get_meta_microsites_individuals_phenotype_values(tr), os.path.join(tmp, tr + '.csv.gz'), ['t%d' % t], cols, mode='a'))
        rows.append((t + 1, round(time.perf_counter() - ts, 4), isl.get_meta_individual_num(), sum(m.get_meta_individual_num() for m in mainlands)))
        if (t + 1) % 10 == 0: print('[%s %s] t=%d/%d %.1fs islands=%d rss=%.0fMB' % (impl, case, t + 1, steps, time.perf_counter() - t0, rows[-1][2], rss_mb()), flush=True)
    total = time.perf_counter() - t0
    shutil.rmtree(tmp)
    out = {'impl': impl, 'case': case, 'kind': kind, 'steps': steps, 'params': P if case != 'model' else None, 'repo': repo, 'metacommunity_py': mc_file, 'build_s': round(build_s, 2), 'init_s': round(init_s, 2), 'setup_s': round(setup_s, 2),
           'total_s': round(total, 2), 's_per_step': round(total / steps, 4), 'phase_s': {k: round(v, 3) for k, v in sorted(phase.items(), key=lambda kv: -kv[1])},
           'rows': rows, 'rss_mb_end': round(rss_mb(), 1), 'islands_alive_mean': round(sum(r[2] for r in rows) / len(rows), 1), 'load_at_end': os.getloadavg()}
    out['loci_per_trait'] = P['loci'] if case != 'model' else 20
    if F: out.update(nemo_minted=sum(f.minted_num for f in F.values()), nemo_recycled=sum(f.recycled_num for f in F.values()), nemo_loci_per_trait=loci_per_trait)
    os.makedirs(BENCH, exist_ok=True)
    json.dump(out, open(os.path.join(BENCH, '%s__%s.json' % (case, impl)), 'w'), indent=1)
    print('DONE %s %s: %.3f s/step over %d steps (build %.1f s, init %.1f s), islands alive mean %.0f' % (impl, case, out['s_per_step'], steps, build_s, init_s, out['islands_alive_mean']))


# ------------------------------------------------------------------ report
def report():
    R = {}
    for f in sorted(os.listdir(BENCH)):
        if f.endswith('.json'): d = json.load(open(os.path.join(BENCH, f))); R[(d['case'], d['impl'])] = d
    sp = lambda a, b: '%.2fx' % (a / b) if b else '-'
    L = ['# Benchmark: MetaIBM-Nemo vs pure MetaIBM (sexual)', '', 'Generated %s by `test/metaibm-nemo-bench/bench.py` on `%s`, python `%s`.' % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), platform.platform(), sys.version.split()[0]), '',
         '- nemo = `%s` (NemoBackend); pure = `%s` (pure-Python MetaIBM, read only)' % (NEMO_REPO, PURE_REPO),
         '- same landscape / parameters / seed / call sequence on both sides; `is_remove=True` passed wherever the function takes it (same semantics: things move); one case at a time on an idle machine (load at the end of each case is listed)',
         '- genome size is identical on both sides: Nemo quanti_loci = loci per trait x traits (asserted per case via get_num_loci; an earlier run that passed loci per trait as the Nemo total was discarded)',
         '- speed-up = pure / nemo (> 1: MetaIBM-Nemo faster). Population sizes evolve stochastically and differ between the two runs: compare also "s per 1000 island individuals per step".', '']
    L += ['## Interfaces (experiment setups: islands start empty, marker schedule, %d steps, CSV every 10 steps)' % STEPS['interface'], '',
          '| case | pure s/step | nemo s/step | speed-up | islands alive mean pure / nemo | per 1000 individuals pure / nemo | setup pure / nemo s | RSS end pure / nemo MB |', '|---|---|---|---|---|---|---|---|']
    for c in INTERFACES:
        a, b = R.get((c, 'pure')), R.get((c, 'nemo'))
        if not (a and b): L.append('| %s | %s | %s | | | | | |' % (c, a and a['s_per_step'], b and b['s_per_step'])); continue
        L.append('| %s | %.3f | %.3f | **%s** | %.0f / %.0f | %.4f / %.4f (%s) | %.1f / %.1f | %.0f / %.0f |' % (c, a['s_per_step'], b['s_per_step'], sp(a['s_per_step'], b['s_per_step']), a['islands_alive_mean'], b['islands_alive_mean'],
                 1000 * a['s_per_step'] / max(1, a['islands_alive_mean']), 1000 * b['s_per_step'] / max(1, b['islands_alive_mean']), sp(a['s_per_step'] / max(1, a['islands_alive_mean']), b['s_per_step'] / max(1, b['islands_alive_mean'])),
                 a['setup_s'], b['setup_s'], a['rss_mb_end'], b['rss_mb_end']))
    L += ['', '## Boundaries (islands fully occupied from step 0, 102400 sites in 16 patches, %d steps): habitat size x schedule' % STEPS['boundary'], '',
          '| habitat | schedule | pure s/step | nemo s/step | speed-up | init 102400 + 6400 founders pure / nemo s | RSS end pure / nemo MB |', '|---|---|---|---|---|---|---|']
    for c in BOUNDARIES:
        a, b = R.get((c, 'pure')), R.get((c, 'nemo'))
        Lh, s = c.split('-')
        if not (a and b): L.append('| %s | %s | %s | %s | | | |' % (Lh, s, a and a['s_per_step'], b and b['s_per_step'])); continue
        L.append('| %sx%s | %s | %.3f | %.3f | **%s** | %.1f / %.1f | %.0f / %.0f |' % (Lh[1:], Lh[1:], s, a['s_per_step'], b['s_per_step'], sp(a['s_per_step'], b['s_per_step']), a['init_s'], b['init_s'], a['rss_mb_end'], b['rss_mb_end']))
    L += ['', '## Parameter sweeps (islands 16 x 16 habitats 20x20 fully occupied, %d steps; one parameter changed from %s)' % (STEPS['sweep'], BASE), '',
          '| variant | change | schedule | pure s/step | nemo s/step | speed-up | init pure / nemo s | RSS end pure / nemo MB |', '|---|---|---|---|---|---|---|---|']
    for v in SWEEP:
        for s_ in ('marker', 'pool'):
            c = 'sweep-%s-%s' % (s_, v); a, b = R.get((c, 'pure')), R.get((c, 'nemo'))
            ch = ', '.join('%s=%s' % kv for kv in SWEEP[v].items()) or 'baseline'
            if not (a and b): L.append('| %s | %s | %s | %s | %s | | | |' % (v, ch, s_, a and a['s_per_step'], b and b['s_per_step'])); continue
            L.append('| %s | %s | %s | %.3f | %.3f | **%s** | %.1f / %.1f | %.0f / %.0f |' % (v, ch, s_, a['s_per_step'], b['s_per_step'], sp(a['s_per_step'], b['s_per_step']), a['init_s'], b['init_s'], a['rss_mb_end'], b['rss_mb_end']))
    L += ['', '## Per phase (seconds per step)', '']
    for c in INTERFACES + BOUNDARIES + SWEEPS:
        a, b = R.get((c, 'pure')), R.get((c, 'nemo'))
        if not (a and b): continue
        L += ['<details><summary>%s</summary>' % c, '', '| phase | pure | nemo | speed-up |', '|---|---|---|---|']
        for k in sorted(set(a['phase_s']) | set(b['phase_s']), key=lambda k: -a['phase_s'].get(k, 0)):
            x, y = a['phase_s'].get(k, 0) / a['steps'], b['phase_s'].get(k, 0) / b['steps']
            L.append('| %s | %.4f | %.4f | %s |' % (k, x, y, sp(x, y)))
        L += ['', '- load at end: pure %s, nemo %s' % (tuple(round(v, 1) for v in a['load_at_end']), tuple(round(v, 1) for v in b['load_at_end'])), '</details>', '']
    open(os.path.join(HERE, 'BENCH_REPORT.md'), 'w').write('\n'.join(L) + '\n')
    print('report:', os.path.join(HERE, 'BENCH_REPORT.md'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--impl', choices=('nemo', 'pure')); ap.add_argument('--case', choices=INTERFACES + BOUNDARIES + SWEEPS); ap.add_argument('--steps', type=int, default=0)
    ap.add_argument('--all', action='store_true'); ap.add_argument('--only', nargs='*', default=[]); ap.add_argument('--report', action='store_true')
    a = ap.parse_args()
    if a.report: report(); sys.exit(0)
    if a.all:
        for c in (a.only or INTERFACES + BOUNDARIES + SWEEPS):
            for impl in ('pure', 'nemo'):                     # one process at a time: no competition for cores
                t = time.time()
                p = subprocess.run([sys.executable, os.path.abspath(__file__), '--impl', impl, '--case', c] + (['--steps', str(a.steps)] if a.steps else []), cwd=HERE, env=ENV)
                print('%s %s rc=%d %.0fs' % (c, impl, p.returncode, time.time() - t), flush=True)
        report(); sys.exit(0)
    run(a.impl, a.case, a.steps)
