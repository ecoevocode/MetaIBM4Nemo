#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" E: full sexual simulations with the global-habitat-network extension (NemoBackend). Nothing in metaibm/, extension/ or experiments/ is edited.
The landscape generators and time loops of experiments/model-sloss-global-habitat-network.py ('net') and experiments/model-sloss-GRFE.py ('grfe') are copied
here (via test/metaibm-nemo-run/simulate.py) with only these interface changes:
  * sexual only: meta_mainland_sexual_birth_mutate_germinate / meta_sex_reproduce_calculation_with_offspring_marker_pool replace the experiments' mixed functions
    (asex / mix are not migrated); meta_initialize(indfactories, species_2_phenotype_ls, 'sexual') with one Nemo indfactory per species
  * every habitat gets a distinct global_index (mainland 0.., metacommunity 100..)
  * 'net': the experiment's sexual branch calls the PATCH-level marker dispersal (network line commented out); here the network line is used (as in its asexual branch)
Groups:  mT   marker flow + network marker dispersal, is_remove default (True)       mF   ... is_remove=False (v3 behaviour; variant v3 = v3 bodies compiled in)
         pool individual-pool flow + network offspring dispersal (is_remove=True)     mPatch  marker flow + the experiments' commented patch-level marker dispersal
Scales:  quick = net 80x80 metacommunity (habitats 10x10), grfe habitats 3x3, mainland habitats 8x8;  full = the experiments' sizes (net 40x40, grfe 10x10, mainland 20x20)
    python simulate_ext.py run --scenario net|grfe --group mT|mF|pool|mPatch --steps N --seed S --scale quick|full [--variant cur|v3] [--digest-every K]
                               [--check-every 100] [--dist RATE] [--watch-dist] --out results/sim_<name>.json
    python simulate_ext.py onestep --scenario grfe --scale quick --seed S --mode true|false --out FILE      (D4: births of step 0 from the same state)
"""
import os, sys, re, json, time, random, argparse, tempfile, traceback
from collections import defaultdict

sys.dont_write_bytecode = True
os.environ.setdefault('MPLBACKEND', 'Agg')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import warnings                                                  # noqa: E402
warnings.filterwarnings('ignore', category=RuntimeWarning); warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
import common as C                                               # noqa: E402  (imports metaibm via the metacommunity suite)
import numpy as np                                               # noqa: E402
import pandas as pd                                              # noqa: E402
from metaibm.indfactory import indfactory                        # noqa: E402
from metaibm.patch import patch                                  # noqa: E402
from metaibm.metacommunity import metacommunity                  # noqa: E402
from pynemo.params import SimParameters                          # noqa: E402

EXPERIMENTS = os.path.join(C.REPO, 'experiments')
NEMO_CACHE = tempfile.mkdtemp(prefix='nemo_ext_sim_')
OFF, MARK, PATCH_MARK = C.OFF, C.MARK, 'dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool'
SCEN = {'net': dict(src='experiments/model-sloss-global-habitat-network.py', patch_num=16, rho=2, among=0.1, within=0.1, dist=0.00001),
        'grfe': dict(src='experiments/model-sloss-GRFE.py', patch_num=16, rho=0.2, among=0.01, within=0.1, dist=0.00001)}
SCALE = {'quick': dict(meta_length=80, grfe_hab=3, mainland_hab=8), 'full': dict(meta_length=320, grfe_hab=10, mainland_hab=20)}


# ====================================================================== species = Nemo indfactories (as test/metaibm-nemo-run/simulate.py)
def nemo_params(name, init_ls, loci=20, mutation_rate=0.0001, mutation_variance=0.01, init_variance=0.002025, env_variance=0.025 ** 2):
    ''' len(init_ls) quanti traits, `loci` loci PER TRAIT (Nemo quanti_loci is the total), environmental variance = pheno_var ** 2 '''
    fmt = lambda v: '{{%s}}' % ', '.join('%g' % x for x in v)
    rows = [('replicates', 1), ('generations', 1), ('patch_number', 1), ('patch_capacity', 1), ('aging', 1),
            ('run_mode', 'overwrite'), ('root_dir', NEMO_CACHE), ('logfile', NEMO_CACHE + '/nemo.log'), ('filename', name),
            ('quanti_traits', len(init_ls)), ('quanti_loci', loci * len(init_ls)), ('quanti_allele_model', 'continuous'),
            ('quanti_init_model', 2), ('quanti_init_value', fmt(init_ls)), ('quanti_init_variance', fmt([init_variance] * len(init_ls))),
            ('quanti_mutation_rate', mutation_rate), ('quanti_mutation_variance', mutation_variance),
            ('quanti_environmental_variance', fmt([env_variance] * len(init_ls)))]
    p = SimParameters()
    for k, v in rows: p.add(k, v)
    return p


def make_indfactories(species_2_phenotype_ls, pheno_names_ls, seed, loci=20, mutation_rate=0.0001, pheno_var=0.025):
    ''' species 'sp<k+1>' <-> species_2_phenotype_ls[k] (the convention of habitat.hab_initialize) '''
    return {'sp%d' % (k + 1): indfactory('sp%d' % (k + 1), list(pheno_names_ls), nemo_params('sp%d' % (k + 1), ph, loci=loci, mutation_rate=mutation_rate, env_variance=pheno_var ** 2), seed=seed + k)
            for k, ph in enumerate(species_2_phenotype_ls)}


class GlobalIndex():
    ''' distinct habitat.global_index values (Nemo homepatch, 0..65535) across mainland + metacommunity '''
    def __init__(self, start=0): self.next = start
    def __call__(self):
        self.next += 1
        assert self.next - 1 <= 65535
        return self.next - 1


# ====================================================================== landscapes (copied from the experiments via simulate.py; + global_index, + cls)
def new_round(_float, _len=0):
    if isinstance(_float, float):
        if str(_float)[::-1].find('.') <= _len: return (_float)
        if str(_float)[-1] == '5': return (round(float(str(_float)[:-1] + '6'), _len))
        else: return (round(_float, _len))
    else: return (round(_float, _len))


def sloss_generating_empty_metacommunity(meta_name, patch_num, patch_location_ls, hab_num, hab_length, hab_width, dormancy_pool_max_size,
                                         x_axis_environment_values_ls, y_axis_environment_values_ls, environment_types_num, environment_types_name,
                                         environment_variation_ls, is_heterogeneity=False, gidx=None, cls=metacommunity):
    log_info = 'generating empty metacommunity ... \n'
    meta_object = cls(metacommunity_name=meta_name)
    patch_num_x_axis, patch_num_y_axis = int(np.sqrt(patch_num)), int(np.sqrt(patch_num))
    hab_num_x_axis, hab_num_y_axis = int(np.sqrt(hab_num)), int(np.sqrt(hab_num))
    for i in range(0, patch_num):
        patch_name = 'patch%d' % (i + 1)
        location = patch_location_ls[i]
        p = patch(patch_name, i, location)
        patch_x_loc, patch_y_loc = location[0], location[1]
        for j in range(hab_num):
            habitat_name = 'h%s' % str(j)
            hab_x_loc, hab_y_loc = patch_x_loc * hab_num_x_axis + j // hab_num_y_axis, patch_y_loc * hab_num_y_axis + j % hab_num_y_axis
            hab_location = (hab_x_loc, hab_y_loc)
            if is_heterogeneity == False:
                x_axis_environment_means_value = x_axis_environment_values_ls[int(hab_x_loc * 4 / (patch_num_x_axis * hab_num_x_axis)) % 4]
                y_axis_environment_means_value = y_axis_environment_values_ls[int(hab_y_loc * 4 / (patch_num_y_axis * hab_num_y_axis)) % 4]
            p.add_habitat(hab_name=habitat_name, hab_index=j, hab_location=hab_location, num_env_types=environment_types_num, env_types_name=environment_types_name,
                          mean_env_ls=[x_axis_environment_means_value, y_axis_environment_means_value], var_env_ls=environment_variation_ls, length=hab_length, width=hab_width,
                          dormancy_pool_max_size=dormancy_pool_max_size, global_index=gidx())
        meta_object.add_patch(patch_name=patch_name, patch_object=p)
    return meta_object, log_info


def sloss_generate_empty_mainland(meta_name, patch_num, patch_location_ls, hab_num, hab_length, hab_width, dormancy_pool_max_size,
                                  micro_environment_values_ls, macro_environment_values_ls, environment_types_num, environment_types_name, environment_variation_ls, gidx=None):
    log_info = 'generating empty mainland ... \n'
    meta_object = metacommunity(metacommunity_name=meta_name)
    for i in range(0, patch_num):
        patch_name = 'patch%d' % (i)
        location = patch_location_ls[i]
        patch_x_loc, patch_y_loc = location[0], location[1]
        p = patch(patch_name, i, location)
        hab_num_length, hab_num_width = int(np.sqrt(hab_num)), int(np.sqrt(hab_num))
        for j in range(hab_num):
            hab_location = (patch_x_loc * hab_num_length + j // hab_num_width, patch_y_loc * hab_num_width + j % hab_num_width)
            p.add_habitat(hab_name='h%s' % str(j), hab_index=j, hab_location=hab_location, num_env_types=environment_types_num, env_types_name=environment_types_name,
                          mean_env_ls=[micro_environment_values_ls[j // hab_num_width], macro_environment_values_ls[j % hab_num_width]], var_env_ls=environment_variation_ls,
                          length=hab_length, width=hab_width, dormancy_pool_max_size=dormancy_pool_max_size, global_index=gidx())
        meta_object.add_patch(patch_name=patch_name, patch_object=p)
    return meta_object, log_info


def read_landscape_configuration(layout_csv, env1_csv, env2_csv, patch_num, is_same_heterogeneity):
    def env_grid(env_csv):
        df = pd.read_csv(env_csv, header=[0, 1]); df.columns = pd.MultiIndex.from_tuples(df.columns)
        key = (f'patch_num={patch_num}', str(bool(is_same_heterogeneity)))
        sub = df.loc[:, [('meta', 'x_index'), ('meta', 'y_index'), key]].copy()
        sub.columns = ['x_index', 'y_index', 'value']
        sub['x_index'], sub['y_index'] = sub['x_index'].astype(int), sub['y_index'].astype(int)
        return sub
    df_layout = pd.read_csv(layout_csv); df_layout = df_layout[df_layout['patch_num'] == patch_num].copy()
    env = env_grid(env1_csv).merge(env_grid(env2_csv), on=['x_index', 'y_index'], how='inner')
    env.columns = ['x_index', 'y_index', 'env_axis1', 'env_axis2']
    return df_layout.merge(env, left_on=['habitat_x_location', 'habitat_y_location'], right_on=['x_index', 'y_index'], how='left')


def grfe_generate_empty_metacommunity(meta_name, patch_num, is_same_heterogeneity, num_env_types, env_types_name, var_env_ls, hab_length, hab_width, dormancy_pool_max_size, gidx=None, cls=metacommunity):
    layout = read_landscape_configuration(os.path.join(EXPERIMENTS, 'patch_habitat_layouts.csv'), os.path.join(EXPERIMENTS, '32x32_habitats_env1.csv'),
                                          os.path.join(EXPERIMENTS, '32x32_habitats_env2.csv'), patch_num, is_same_heterogeneity)
    meta_obj = cls(metacommunity_name=meta_name)
    reg = layout[['patch_id', 'patch_index', 'patch_location_x', 'patch_location_y']].drop_duplicates().sort_values(['patch_index', 'patch_id'])
    for _, pr in reg.iterrows():
        patch_name = pr['patch_id']
        patch_object = patch(patch_name, int(pr['patch_index']), (float(pr['patch_location_x']), float(pr['patch_location_y'])))
        for _, hr in layout[layout['patch_id'] == patch_name].copy().sort_values(['habitat_index', 'habitat_id']).iterrows():
            patch_object.add_habitat(hab_name=str(hr['habitat_id']), hab_index=int(hr['habitat_index']), hab_location=(int(hr['habitat_x_location']), int(hr['habitat_y_location'])),
                                     num_env_types=num_env_types, env_types_name=env_types_name, mean_env_ls=[float(hr['env_axis1']), float(hr['env_axis2'])], var_env_ls=var_env_ls,
                                     length=hab_length, width=hab_width, dormancy_pool_max_size=dormancy_pool_max_size, global_index=gidx())
        meta_obj.add_patch(patch_name=patch_name, patch_object=patch_object)
    return meta_obj, ''


SPECIES_2_PHENOTYPE = [[i / 10, j / 10] for i in range(2, 10, 2) for j in range(2, 10, 2)]


def build(scenario, scale, seed, cls=metacommunity):
    ''' as the experiments: only the mainland is initialized, the metacommunity starts empty and is colonized by propagule rain '''
    np.random.seed(seed); random.seed(seed)
    S = SCALE[scale]
    etn = ('x_axis_environment', 'y_axis_environment') if scenario == 'net' else ('env_axis1', 'env_axis2')
    gidx = GlobalIndex(0)
    ml, _ = sloss_generate_empty_mainland('mainland', 1, [(0, 0)], 16, S['mainland_hab'], S['mainland_hab'], 0, [0.2, 0.4, 0.6, 0.8], [0.2, 0.4, 0.6, 0.8], 2, etn, [0.025, 0.025], gidx=gidx)
    gidx.next = 100
    t = time.time()
    if scenario == 'net':
        patch_num = SCEN['net']['patch_num']
        px = int(np.sqrt(patch_num))
        hab_num_in_patch = 16 if patch_num == 1 else 4
        hab_len = int(S['meta_length'] / px / int(np.sqrt(hab_num_in_patch)))
        meta, _ = sloss_generating_empty_metacommunity('metacommunity', patch_num, [(i, j) for i in range(px) for j in range(px)], hab_num_in_patch, hab_len, hab_len,
                                                       int(new_round(0 / (patch_num * hab_num_in_patch))), [0.2, 0.4, 0.6, 0.8], [0.2, 0.4, 0.6, 0.8], 2, etn, [0.025, 0.025], False, gidx=gidx, cls=cls)
        pheno_names = ('x_axis_phenotype', 'y_axis_phenotype')
    else:
        meta, _ = grfe_generate_empty_metacommunity('metacommunity', SCEN['grfe']['patch_num'], True, 2, etn, [0.025, 0.025], S['grfe_hab'], S['grfe_hab'], 0, gidx=gidx, cls=cls)
        pheno_names = ('phenotype_axis1', 'phenotype_axis2')
    build_s = time.time() - t
    F = make_indfactories(SPECIES_2_PHENOTYPE, pheno_names, seed=1000 * seed + 1, loci=20, mutation_rate=0.0001, pheno_var=0.025)
    ml.meta_initialize(F, SPECIES_2_PHENOTYPE, 'sexual')
    return ml, meta, F, build_s


# ====================================================================== one time step (the experiments' sexual branch, network line for the among-patch dispersal)
RX_BIRTH = re.compile(r'there are (\d+) individuals germinating')
RX_NET = re.compile(r'there are (\d+) (?:individuals|offspring markers) disperse into immigrant(?:_marker)?_pool among patches in global-habitat-network')
RX_PATCH = re.compile(r'there are (\d+) individuals disperse into habs_immigrant_marker_pool among patches')
RX_COL = re.compile(r'there are (\d+) individuals colonizing')


def step_ops(group, ml, meta, sc, variant, dist):
    rho, among, within = sc['rho'], sc['among'], sc['within']
    if group == 'mT': net = ('network_marker', lambda: getattr(meta, MARK)(among, method='exponential', rho=rho))
    elif group == 'mF':
        net = ('network_marker', (lambda: getattr(meta, MARK)(among, method='exponential', rho=rho)) if variant == 'v3' else (lambda: getattr(meta, MARK)(among, method='exponential', is_remove=False, rho=rho)))
    elif group == 'mPatch': net = ('patch_marker', lambda: getattr(meta, PATCH_MARK)(among, method='exponential', rho=rho))
    if group == 'pool':
        return [('mainland.meta_dead_selection', lambda: ml.meta_dead_selection(0.1, 0.5, method='niche_gaussian')),
                ('mainland.meta_mainland_sexual_birth_mutate_germinate', lambda: ml.meta_mainland_sexual_birth_mutate_germinate(1)),
                ('meta_dead_selection', lambda: meta.meta_dead_selection(0.1, 0.5, method='niche_gaussian')),
                ('meta_sex_reproduce_mutate_into_offspring_pool', lambda: meta.meta_sex_reproduce_mutate_into_offspring_pool(1)),
                ('pairwise_sexual_colonization_from_prpagules_rains', lambda: meta.pairwise_sexual_colonization_from_prpagules_rains(ml, 6400 * 0.001)),
                ('meta_dispersal_within_patch_from_offspring_to_immigrant_pool', lambda: meta.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(within)),
                ('network_offspring', lambda: getattr(meta, OFF)(among, method='exponential', rho=rho)),
                ('meta_local_germinate_from_offspring_immigrant_and_dormancy_pool', lambda: meta.meta_local_germinate_from_offspring_immigrant_and_dormancy_pool()),
                ('meta_dormancy_process_from_offspring_pool_and_immigrant_pool', lambda: meta.meta_dormancy_process_from_offspring_pool_and_immigrant_pool()),
                ('meta_disturbance_process_in_patches', lambda: meta.meta_disturbance_process_in_patches(dist)),
                ('meta_clear_up_offspring_and_immigrant_pool', lambda: meta.meta_clear_up_offspring_and_immigrant_pool())]
    return [('mainland.meta_dead_selection', lambda: ml.meta_dead_selection(0.1, 0.5, method='niche_gaussian')),
            ('mainland.meta_mainland_sexual_birth_mutate_germinate', lambda: ml.meta_mainland_sexual_birth_mutate_germinate(1)),
            ('meta_dead_selection', lambda: meta.meta_dead_selection(0.1, 0.5, method='niche_gaussian')),
            ('meta_sex_reproduce_calculation_with_offspring_marker_pool', lambda: meta.meta_sex_reproduce_calculation_with_offspring_marker_pool(1)),
            ('pairwise_sexual_colonization_from_prpagules_rains', lambda: meta.pairwise_sexual_colonization_from_prpagules_rains(ml, 6400 * 0.001)),
            ('meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool', lambda: meta.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(within)),
            net,
            ('meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool', lambda: meta.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool()),
            ('meta_disturbance_process_in_patches', lambda: meta.meta_disturbance_process_in_patches(dist)),
            ('meta_clear_up_offspring_marker_and_immigrant_marker_pool', lambda: meta.meta_clear_up_offspring_marker_and_immigrant_marker_pool())]


# ====================================================================== observation helpers
class Network():
    ''' row index / distance / patch of every registered habitat (from the registry), global_index -> row '''
    def __init__(self, meta):
        self.ls = list(meta.global_habitat_id_ls); self.row = {k: i for i, k in enumerate(self.ls)}
        self.D = np.asarray(meta.global_habitat_distance_matrix, dtype=float)
        self.g2row = {meta.set[p].set[h].global_index: i for i, (p, h) in enumerate(self.ls)}
        self.habs = [meta.set[p].set[h] for p, h in self.ls]
        self.pairs = defaultdict(int)                               # distance bin -> number of (i, j) pairs in different patches
        for i, (pi, _) in enumerate(self.ls):
            for j, (pj, _) in enumerate(self.ls):
                if pi != pj: self.pairs[int(self.D[i, j])] += 1

    def immigrants(self, flow):
        ''' (source row, dest row) of every among-patch immigrant now in the immigrant pools: markers carry (birth patch, birth habitat), individuals their Nemo home (= birth habitat global_index) '''
        out = []
        for j, h in enumerate(self.habs):
            pj = self.ls[j][0]
            if flow == 'pool':
                for x in h.immigrant_pool:
                    i = self.g2row.get(x.nemo_indi_obj.getHome())
                    if i is not None and self.ls[i][0] != pj: out.append((i, j))
            else:
                for mk in h.immigrant_marker_pool:
                    if mk[0] != pj: out.append((self.row[(mk[0], mk[1])], j))
        return out


def species_counts(m):
    c = defaultdict(int)
    for p in m.set.values():
        for h in p.set.values():
            for s, d in h.species_category.items(): c[s] += sum(len(v) for v in d.values())
    return {k: v for k, v in c.items() if v}


def sites(m): return sum(h.size for p in m.set.values() for h in p.set.values())


def patch_alive(p): return sum(h.indi_num for h in p.set.values())


def patch_signature(p):
    ''' cheap fingerprint of a patch: occupant ids per habitat + pool lengths (an untouched patch keeps it) '''
    return [(hn, tuple(sorted(id(x) for x in C.meta_common.micro(h))), len(h.offspring_pool), len(h.immigrant_pool), len(h.dormancy_pool), len(h.offspring_marker_pool), len(h.immigrant_marker_pool))
            for hn, h in p.set.items()]


def watch_disturbance(meta, F, events, cur_t):
    ''' wrap every patch's patch_disturbance_process (instance attribute): record what it held, what was recycled, that it is empty afterwards and the other patches untouched '''
    fs = list(F.values())
    for pn, p in meta.set.items():
        real = p.patch_disturbance_process
        def wrapped(pn=pn, p=p, real=real):
            others = {qn: patch_signature(q) for qn, q in meta.set.items() if qn != pn}
            micro = sum(h.indi_num for h in p.set.values()); pools = sum(len(h.offspring_pool) + len(h.immigrant_pool) for h in p.set.values())
            dorm = sum(len(h.dormancy_pool) for h in p.set.values()); markers = sum(len(h.offspring_marker_pool) + len(h.immigrant_marker_pool) for h in p.set.values())
            rec0 = sum(f.recycled_num for f in fs)
            ret = real()
            rec = sum(f.recycled_num for f in fs) - rec0
            events.append({'t': cur_t[0], 'patch': pn, 'held_micro': micro, 'held_pools': pools, 'held_dormancy': dorm, 'markers': markers, 'recycled': rec, 'returned': ret,
                           'recycled_ok': rec == micro + pools and ret == micro + pools + markers, 'empty_after': patch_alive(p) == 0 and all(not h.offspring_pool and not h.immigrant_pool for h in p.set.values()),
                           'dormancy_kept': sum(len(h.dormancy_pool) for h in p.set.values()) == dorm,
                           'others_untouched': all(patch_signature(meta.set[qn]) == sig for qn, sig in others.items()), 'pre_alive': micro, 'recolonize': []})
            return ret
        p.patch_disturbance_process = wrapped


# ====================================================================== run
def run(a):
    sc = dict(SCEN[a.scenario]); dist = a.dist if a.dist is not None else sc['dist']
    flow = 'pool' if a.group == 'pool' else 'marker'
    cls = metacommunity
    if a.variant == 'v3':
        assert a.group == 'mF', 'v3 bodies only for the copy-mode marker run (the v3 offspring version gives two owners)'
        cls = C.original_metacommunity(ext_names=[C.ORIGINAL_REGISTRY, MARK], meta_names=['add_patch'])
    out = {'args': vars(a), 'status': 'running', 'scenario_src': sc['src'], 'rows': [], 'digests': [], 'phase_seconds': {}, 'error': None}
    path = a.out
    def write():
        tmp = path + '.tmp'
        with open(tmp, 'w') as fh: json.dump(out, fh, default=str)
        os.replace(tmp, path)
    t0 = time.time()
    ml, meta, F, build_s = build(a.scenario, a.scale, a.seed, cls)
    fs = list(F.values())
    net = Network(meta)
    out.update(H=len(net.ls), patches=meta.patch_num, sites=sites(meta), mainland_sites=sites(ml), build_s=round(build_s, 2), init_s=round(time.time() - t0, 2),
               pairs_by_distance=dict(net.pairs), gidx_distinct=len(net.g2row) == len(net.ls))
    phase = defaultdict(float)
    series = {k: [] for k in ('alive', 'births', 'migrants', 'colonizers', 'step_s', 'richness')}
    flow_bins = defaultdict(int); in_by_hab = [0] * len(net.ls); dist_events = []; cur_t = [0]
    if a.watch_dist: watch_disturbance(meta, F, dist_events, cur_t)
    name, t = '', -1
    try:
        for t in range(a.steps):
            cur_t[0] = t; ts = time.perf_counter()
            births = migrants = col = 0
            for name, op in step_ops(a.group, ml, meta, sc, a.variant, dist):
                tp = time.perf_counter(); ret = op(); phase[name] += time.perf_counter() - tp
                ret = ret if isinstance(ret, str) else ''
                if name in ('network_marker', 'network_offspring', 'patch_marker'):
                    mt = (RX_PATCH if name == 'patch_marker' else RX_NET).search(ret); migrants = int(mt.group(1)) if mt else -1
                    tp = time.perf_counter()
                    for i, j in net.immigrants(flow): flow_bins[int(net.D[i, j])] += 1; in_by_hab[j] += 1
                    phase['(observe immigrants)'] += time.perf_counter() - tp
                mt = RX_BIRTH.search(ret)
                if mt and not name.startswith('mainland'): births = int(mt.group(1))
                mt = RX_COL.search(ret)
                if mt: col = int(mt.group(1))
            series['step_s'].append(round(time.perf_counter() - ts, 4))
            series['alive'].append(meta.get_meta_individual_num()); series['births'].append(births); series['migrants'].append(migrants); series['colonizers'].append(col)
            series['richness'].append(len(species_counts(meta)))
            for e in dist_events:
                if e['t'] < t and len(e['recolonize']) < 400: e['recolonize'].append(patch_alive(meta.set[e['patch']]))
            if a.digest_every and (t % a.digest_every == 0 or t == a.steps - 1):
                tp = time.perf_counter(); out['digests'].append([t, C.digest([ml, meta]), C.rng_digest()]); phase['(digest)'] += time.perf_counter() - tp
            if t == 0 or (t + 1) % a.check_every == 0 or t == a.steps - 1:
                tp = time.perf_counter(); C.check_world([ml, meta], fs, '%s t=%d' % (a.out, t)); phase['(invariant checks)'] += time.perf_counter() - tp
                out['rows'].append({'t': t, 'elapsed_s': round(time.time() - t0, 1), 'rss_mb': round(C.meta_common.rss_mb(), 1), 'alive': meta.get_meta_individual_num(),
                                    'mainland_alive': ml.get_meta_individual_num(), 'occupancy': round(meta.get_meta_individual_num() / out['sites'], 4),
                                    'species': species_counts(meta), 'minted': C.minted(F), 'recycled': C.recycled(F), 'held': len(C.held([ml, meta])), 'invariants': 'ok'})
                out.update(series=series, phase_seconds={k: round(v, 2) for k, v in sorted(phase.items(), key=lambda kv: -kv[1])}, flow_bins=dict(flow_bins),
                           immigrants_by_habitat=in_by_hab, disturbances=dist_events, steps_done=t + 1)
                write()
        tp = time.perf_counter()
        n = C.teardown([ml, meta], fs)
        bal = {f.species_id: [f.minted_num, f.recycled_num] for f in fs}
        for f in fs: f.purge_recycling_pool()
        phase['(teardown)'] += time.perf_counter() - tp
        out.update(status='done', teardown_recycled=n, teardown_balance=bal, teardown_ok=all(x == y for x, y in bal.values()))
    except BaseException as e:
        out.update(status='error', error='%s: %s' % (type(e).__name__, str(e)[:600]), at='t=%d %s' % (t, name), traceback=traceback.format_exc()[-3000:])
    out.update(series=series, phase_seconds={k: round(v, 2) for k, v in sorted(phase.items(), key=lambda kv: -kv[1])}, flow_bins=dict(flow_bins), immigrants_by_habitat=in_by_hab,
               disturbances=dist_events, steps_done=t + 1, total_s=round(time.time() - t0, 1),
               habitat_distance_to_mainland=None, habitats=[[p, h] for p, h in net.ls])
    write()
    sys.stdout.flush(); os._exit(0 if out['status'] == 'done' else 1)          # skip interpreter teardown (a broken run may hold individuals owned twice)


def onestep(a):
    ''' D4: the metacommunity filled with founders (test-side: GRFE habitat means are continuous, meta_initialize needs an exact species match -- each microsite gets a
        founder of the species nearest to the habitat mean), then step 0 of the marker flow up to germination with the network marker dispersal is_remove = a.mode;
        the state before the network call is identical for both modes '''
    sc = SCEN[a.scenario]
    ml, meta, F, _ = build(a.scenario, a.scale, a.seed)
    for p in meta.set.values():
        for h in p.set.values():
            k = int(np.argmin([(x - h.mean_env_ls[0]) ** 2 + (y - h.mean_env_ls[1]) ** 2 for x, y in SPECIES_2_PHENOTYPE]))
            for row in range(h.length):
                for col in range(h.width):
                    h.add_individual(F['sp%d' % (k + 1)].make_founder(random.sample(('male', 'female'), 1)[0], h.global_index), row, col)
    ops = step_ops('mT' if a.mode == 'true' else 'mF', ml, meta, sc, 'cur', sc['dist'])
    rec = {}
    for name, op in ops:
        if name == 'network_marker':
            rec['before'] = [C.digest([ml, meta]), C.rng_digest()]
            rec['offspring_markers_before'] = meta.meta_offspring_marker_pool_marker_num(); rec['immigrant_markers_before'] = meta.meta_immigrant_marker_pool_marker_num()
        ret = op()
        if name == 'network_marker':
            mt = RX_NET.search(ret); rec['network_markers'] = int(mt.group(1)) if mt else None
            rec['offspring_markers_after'] = meta.meta_offspring_marker_pool_marker_num(); rec['immigrant_markers_after'] = meta.meta_immigrant_marker_pool_marker_num()
            rec['empty_sites'] = meta.show_meta_empty_sites_num()
        if name.startswith('meta_local_germinate'):
            rec['births'] = int(RX_BIRTH.search(ret).group(1)); break
    C.check_world([ml, meta], list(F.values()), 'onestep')
    rec.update(seed=a.seed, mode=a.mode)
    json.dump(rec, open(a.out, 'w'))
    sys.stdout.flush(); os._exit(0)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('--scenario', choices=list(SCEN), required=True); r.add_argument('--group', choices=['mT', 'mF', 'pool', 'mPatch'], required=True)
    r.add_argument('--steps', type=int, default=300); r.add_argument('--seed', type=int, default=1); r.add_argument('--scale', choices=list(SCALE), default='quick')
    r.add_argument('--variant', choices=['cur', 'v3'], default='cur'); r.add_argument('--digest-every', type=int, default=0); r.add_argument('--check-every', type=int, default=100)
    r.add_argument('--dist', type=float, default=None); r.add_argument('--watch-dist', action='store_true'); r.add_argument('--out', required=True)
    o = sub.add_parser('onestep')
    o.add_argument('--scenario', default='grfe'); o.add_argument('--scale', default='quick'); o.add_argument('--seed', type=int, default=1)
    o.add_argument('--mode', choices=['true', 'false'], required=True); o.add_argument('--out', required=True)
    a = ap.parse_args()
    if a.cmd == 'run': run(a)
    else: onestep(a)


if __name__ == '__main__':
    main()
