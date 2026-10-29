# Shared helpers for test/nemoback-metacommunity. Every script runs in its own process (Nemo keeps process-wide static state: RAND, Individual::currentID).
# Reuses test/nemoback-patch/common.py (loaded via importlib as `patch_common`) and, through it, test/nemoback-habitat/common.py (`hab_common`, also importlib):
# both files are called common.py. Nothing in metaibm/ or extension/ is edited; original (pre-migration) method bodies are compiled from git into a subclass.
import os, sys, ast, json, random, hashlib, inspect, subprocess, importlib.util
import numpy as np

sys.dont_write_bytecode = True                                  # never write __pycache__ into the other suites' folders
os.environ.setdefault('MPLBACKEND', 'Agg')                      # metacommunity.py imports matplotlib.pyplot
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
PATCH_SUITE = os.path.join(REPO, 'test', 'nemoback-patch')
HAB_SUITE = os.path.join(REPO, 'test', 'nemoback-habitat')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


patch_common = _load('patch_common', os.path.join(PATCH_SUITE, 'common.py'))      # loads test/nemoback-habitat/common.py itself (importlib, as `hab_common`)
hab_common = patch_common.hab_common                                               # the same module object patch_common uses (one RECYCLED list, one NEMO_CACHE)
assert os.path.samefile(hab_common.__file__, os.path.join(HAB_SUITE, 'common.py'))

from metaibm.metacommunity import metacommunity                 # noqa: E402
from metaibm.patch import patch                                 # noqa: E402
from metaibm.habitat import habitat                             # noqa: E402
from metaibm.indfactory import indfactory                       # noqa: E402
from metaibm.individual import individual                       # noqa: E402

make_species, seed_all, chi2_ok, rss_mb, fake_indi = hab_common.make_species, hab_common.seed_all, hab_common.chi2_ok, hab_common.rss_mb, hab_common.fake_indi
binom_ok, rng_state, same_rng_state = patch_common.binom_ok, patch_common.rng_state, patch_common.same_rng_state
POOLS, MARKER_POOLS = patch_common.POOLS, patch_common.MARKER_POOLS
TABLE = [[0.2, 0.2], [0.5, 0.5]]                                # species_2_phenotype_ls: 'sp1' <-> [0.2, 0.2], 'sp2' <-> [0.5, 0.5]
SPECIES_MEANS = [(0.2, 0.2), (0.5, 0.5)]
FINDINGS = []


def git(*a):
    return subprocess.run(['git'] + list(a), capture_output=True, text=True, cwd=REPO).stdout


# ------------------------------------------------------------------ original (pre-migration) metacommunity.py
ORIGINAL_MARK = 'patch_clear_up_offspring_and_immigrant_pool'


def _find_original_commit():
    ''' newest commit (walking back from HEAD) whose metaibm/metacommunity.py still calls patch_clear_up_offspring_and_immigrant_pool (expected 7a1adbb) '''
    for c in git('log', '--format=%h', 'HEAD').split():
        if ORIGINAL_MARK in git('show', c + ':metaibm/metacommunity.py'): return c
    raise RuntimeError('no commit with the pre-migration metacommunity.py found')


ORIGINAL_COMMIT = _find_original_commit()
NEVER_COMPILE = {'meta_initialize', 'meta_clear_up_offspring_and_immigrant_pool', 'meta_clear_up_offspring_marker_and_immigrant_marker_pool'}   # call removed APIs


def original_source():
    return git('show', '%s:metaibm/metacommunity.py' % ORIGINAL_COMMIT)


def class_functions(src, cls_name='metacommunity'):
    cls = [n for n in ast.parse(src).body if isinstance(n, ast.ClassDef) and n.name == cls_name][0]
    return {n.name: n for n in cls.body if isinstance(n, ast.FunctionDef)}


def original_metacommunity(names):
    ''' subclass of the CURRENT metacommunity whose `names` methods are the original bodies from git (compiled verbatim); they run on the current (NemoBackend) patch / habitat '''
    bad = set(names) & NEVER_COMPILE
    if bad: raise ValueError('original %s call removed patch APIs: never compile them' % sorted(bad))
    src = original_source(); lines = src.splitlines(); fn = class_functions(src)
    body = 'class original_metacommunity(metacommunity):\n' + '\n'.join('\n'.join(lines[fn[k].lineno - 1:fn[k].end_lineno]) for k in names) + '\n'
    ns = {'metacommunity': metacommunity, 'patch': patch, 'random': random, 'np': np}
    exec(compile(body, '<git %s:metaibm/metacommunity.py>' % ORIGINAL_COMMIT, 'exec'), ns)
    return ns['original_metacommunity']


# ------------------------------------------------------------------ building worlds
def species(seed=1, **kw):
    ''' the two species of every test: sp1 (founders at 0.2) and sp2 (founders at 0.5), one boot sim each '''
    return {'sp1': make_species('sp1', 0.2, seed=seed, **kw), 'sp2': make_species('sp2', 0.5, seed=seed + 1, **kw)}


def make_islands(n_patches=3, n_habs=2, L=4, W=None, K=10, gidx0=100, means=None, name='islands', cls=metacommunity, prefix='p'):
    ''' metacommunity of n_patches patches (patch.index k, location (2k, 0)) x n_habs habitats; habitat h of patch k: hab_index h, location (2k, h),
        global_index gidx0 + k*n_habs + h (distinct per habitat), mean means[k][h] (default alternates the two species), dormancy_pool_max_size K.
        Seeds nothing: callers seed np.random first (habitat.__init__ draws the environment). '''
    m = cls(name)
    for k in range(n_patches):
        p = patch('%s%d' % (prefix, k), k, (2.0 * k, 0.0))
        for h in range(n_habs):
            mean = list(means[k][h]) if means is not None else list(SPECIES_MEANS[(k + h) % 2])
            p.add_habitat('h%d' % h, h, (2.0 * k, float(h)), 2, ['env1', 'env2'], mean, [0.025, 0.025], L, L if W is None else W, K, global_index=gidx0 + k * n_habs + h)
        m.add_patch(p.name, p)
    return m


def make_mainland(n_habs=2, L=30, W=None, gidx0=0, means=None, cls=metacommunity):
    ''' realistic mainland: one patch 'm0' with n_habs habitats of L x W (default 2 x 30x30 = 1800 sites), no dormancy, one species per habitat '''
    return make_islands(1, n_habs, L=L, W=W, K=0, gidx0=gidx0, means=[means] if means is not None else [[SPECIES_MEANS[h % 2] for h in range(n_habs)]], name='mainland', cls=cls, prefix='m')


def rain_for(mainland, rate=0.001):
    ''' propagules_rain_num as experiments/model-sloss*.py (6400 sites x 0.001 = 6.4), scaled to the mainland size '''
    return sum(h.size for p in mainland.set.values() for h in p.set.values()) * rate


def fill_pools(m, K_birth=1.0):
    ''' drive an initialized metacommunity into a state with empty sites + non-empty offspring / immigrant / dormancy / marker pools (current code only) '''
    m.meta_dead_selection(0.3, 0.5, 'niche_gaussian')
    m.meta_sex_reproduce_mutate_into_offspring_pool(K_birth)
    m.meta_dormancy_process_from_offspring_pool_and_immigrant_pool()
    m.meta_sex_reproduce_mutate_into_offspring_pool(K_birth)
    m.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(0.4)
    m.meta_sex_reproduce_calculation_with_offspring_marker_pool(0.8)
    m.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(0.5)
    m.meta_dead_selection(0.2, 0.5, 'niche_gaussian')
    return m


# ------------------------------------------------------------------ holdings, invariants
def _as_list(x):
    if x is None: return []
    if isinstance(x, dict): return list(x.values())
    return list(x) if isinstance(x, (list, tuple)) else [x]


def all_patches(metas):
    return [m.set[n] for m in _as_list(metas) for n in m.patch_id_ls]


def all_habs(metas):
    return [p.set[n] for p in all_patches(metas) for n in p.hab_id_ls]


def held(metas):
    return patch_common.held(all_patches(metas))


def held_ids(metas):
    return {id(x) for x in held(metas)}


def micro(h):
    return hab_common.microsite_individuals(h)


def check_world(metas, factories, tag=''):
    ''' R5 over a list of metacommunities (mainland + islands): one owner (microsite or pool) for every individual, everything held is alive,
        minted - recycled == held per factory, microsite bookkeeping consistent (via patch_common / hab_common); plus metacommunity bookkeeping and marker sanity '''
    for m in _as_list(metas):
        assert m.patch_num == len(m.set) == len(m.patch_id_ls) and list(m.set) == m.patch_id_ls, '%s %s patch bookkeeping' % (tag, m.metacommunity_name)
        assert all(m.patch_id_2_index_dir[n] == m.set[n].index == i for i, n in enumerate(m.patch_id_ls)), '%s %s patch index' % (tag, m.metacommunity_name)
        for p in all_patches(m):
            for h in p.set.values():
                for q in MARKER_POOLS:
                    for mk in getattr(h, q):
                        assert isinstance(mk, tuple) and len(mk) == 3 and mk[0] in m.set and mk[1] in m.set[mk[0]].set, '%s marker %r in %s.%s.%s' % (tag, mk, m.metacommunity_name, p.name, q)
    return patch_common.check_invariants(all_patches(metas), factories, tag)


def minted(F): return sum(f.minted_num for f in _as_list(F))
def recycled(F): return sum(f.recycled_num for f in _as_list(F))
def indi_n(metas): return sum(h.indi_num for h in all_habs(metas))
def pool_n(metas, names=POOLS): return sum(len(getattr(h, q)) for h in all_habs(metas) for q in names)


def teardown(metas, factories):
    ''' meta_clear_up_microsite_individuals / offspring / immigrant / dormancy (+ markers) on every metacommunity, then minted == recycled per factory '''
    n = 0
    for m in _as_list(metas):
        n += m.meta_clear_up_microsite_individuals() + m.meta_clear_up_offspring_pool() + m.meta_clear_up_immigrant_pool() + m.meta_clear_up_dormancy_pool()
        m.meta_clear_up_offspring_marker_and_immigrant_marker_pool()
    for f in _as_list(factories):
        assert f.minted_num == f.recycled_num, 'teardown: %s minted %d != recycled %d' % (f.species_id, f.minted_num, f.recycled_num)
    return n


# ------------------------------------------------------------------ digests (cross-process comparable)
def _ident(x):
    if hasattr(x, 'label'): return x.label
    return [x.individual_id, x.gender, x.nemo_indi_obj.getHome(), [round(v, 12) for v in x.get_indi_phenotype_ls()]]


def state_rows(metas):
    ''' full state: per habitat the grid (ID, gender, home, phenotype), every pool (in order), marker pools, empty / occupied lists (order matters for later RNG use), species_category '''
    rows = []
    for m in _as_list(metas):
        for p in all_patches(m):
            for hn in p.hab_id_ls:
                h = p.set[hn]
                rows.append([m.metacommunity_name, p.name, hn, [[None if x is None else _ident(x) for x in row] for row in h.set['microsite_individuals']],
                             [[_ident(x) for x in getattr(h, q)] for q in POOLS], [[list(mk) for mk in getattr(h, q)] for q in MARKER_POOLS],
                             [list(s) for s in h.empty_site_pos_ls], [list(s) for s in h.occupied_site_pos_ls], h.indi_num,
                             [[s, [[g, [list(v) for v in vs]] for g, vs in d.items()]] for s, d in h.species_category.items()]])
    return rows


def digest(metas):
    return hashlib.sha1(json.dumps(state_rows(metas)).encode()).hexdigest()[:16]


def rng_digest():
    s = np.random.get_state()
    return hashlib.sha1(repr(random.getstate()).encode() + s[1].tobytes() + repr((s[0], s[2], s[3], s[4])).encode()).hexdigest()[:16]


def full_snapshot(metas, factories):
    ''' everything a guarded call must leave untouched: state digest, both RNG states, factory counters, Nemo's static currentID '''
    return (digest(metas), rng_state(), tuple((f.minted_num, f.recycled_num) for f in _as_list(factories)), indfactory.get_current_individual_id())


def same_snapshot(a, b):
    return a[0] == b[0] and same_rng_state(a[1], b[1]) and a[2:] == b[2:]


# ------------------------------------------------------------------ v3-style schedules (experiments/model-sloss*.py order, sexual)
POOL_VARIANTS = {   # within, among, germinate, dormancy
    'I': ('meta_dispersal_within_patch_from_offspring_to_immigrant_pool', 'dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool',
          'meta_local_germinate_from_offspring_and_immigrant_pool', 'meta_dormancy_process_from_offspring_pool_and_immigrant_pool'),
    'DI': ('meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool', 'dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool',
           'meta_local_germinate_from_offspring_immigrant_and_dormancy_pool', 'meta_dormancy_process_from_offspring_pool_and_immigrant_pool'),
    'sites': ('meta_dispersal_within_patch_from_offspring_and_dormancy_pool', 'dispersal_among_patches_from_offspring_pool_and_dormancy_pool',
              'meta_local_germinate_from_offspring_and_dormancy_pool', 'meta_dormancy_process_from_offspring_pool_to_dormancy_pool')}


def pool_schedule(mainland, islands, t, rain, dorm=True, pairwise=None, m_within=0.3, m_among=0.2, dist=0.05, kernel=('exponential', {'rho': 2})):
    ''' one step of the individual-pool schedule as (name, callable) pairs: mainland dead selection + sexual birth -> colonization -> islands dead selection
        -> reproduce into offspring_pool -> dispersal within / among -> germinate -> dormancy (or clear-ups) -> disturbance -> clear-ups '''
    v = list(POOL_VARIANTS)[t % 3]
    within, among, germ, dormancy = POOL_VARIANTS[v]
    pw = (t % 2 == 1) if pairwise is None else pairwise
    col = 'pairwise_sexual_colonization_from_prpagules_rains' if pw else 'meta_colonize_from_propagules_rains'
    ops = [('mainland.meta_dead_selection', lambda: mainland.meta_dead_selection(0.1, 0.5, 'niche_gaussian')),
           ('mainland.meta_mainland_sexual_birth_mutate_germinate', lambda: mainland.meta_mainland_sexual_birth_mutate_germinate(1.0)),
           ('islands.' + col, lambda: getattr(islands, col)(mainland, rain)),
           ('islands.meta_dead_selection', lambda: islands.meta_dead_selection(0.1, 0.5, 'niche_gaussian')),
           ('islands.meta_sex_reproduce_mutate_into_offspring_pool', lambda: islands.meta_sex_reproduce_mutate_into_offspring_pool(1.0)),
           ('islands.' + within, lambda: getattr(islands, within)(m_within)),
           ('islands.' + among, lambda: getattr(islands, among)(m_among, method=kernel[0], **kernel[1])),
           ('islands.' + germ, lambda: getattr(islands, germ)())]
    if dorm: ops.append(('islands.' + dormancy, lambda: getattr(islands, dormancy)()))
    ops.append(('islands.meta_disturbance_process_in_' + ('patches' if t % 2 else 'habitat'),
                (lambda: islands.meta_disturbance_process_in_patches(dist)) if t % 2 else (lambda: islands.meta_disturbance_process_in_habitat(dist))))
    ops.append(('islands.meta_clear_up_offspring_and_immigrant_pool', islands.meta_clear_up_offspring_and_immigrant_pool))
    return ops


def marker_schedule(mainland, islands, t, rain, pairwise=None, m_within=0.3, m_among=0.2, dist=0.05, kernel=('exponential', {'rho': 2}), birth=0.8, network=False):
    ''' one step of the marker schedule: mainland dead selection + sexual birth -> colonization -> islands dead selection -> offspring markers
        -> marker dispersal within / among (network=True: on odd t the extension's global-habitat-network marker dispersal, as model-sloss-global-habitat-network.py)
        -> germinate + birth from markers -> disturbance -> marker clear-ups '''
    pw = (t % 2 == 1) if pairwise is None else pairwise
    col = 'pairwise_sexual_colonization_from_prpagules_rains' if pw else 'meta_colonize_from_propagules_rains'
    return [('mainland.meta_dead_selection', lambda: mainland.meta_dead_selection(0.1, 0.5, 'niche_gaussian')),
            ('mainland.meta_mainland_sexual_birth_mutate_germinate', lambda: mainland.meta_mainland_sexual_birth_mutate_germinate(1.0)),
            ('islands.meta_dead_selection', lambda: islands.meta_dead_selection(0.1, 0.5, 'niche_gaussian')),
            ('islands.meta_sex_reproduce_calculation_with_offspring_marker_pool', lambda: islands.meta_sex_reproduce_calculation_with_offspring_marker_pool(birth)),
            ('islands.' + col, lambda: getattr(islands, col)(mainland, rain)),
            ('islands.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool', lambda: islands.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(m_within)),
            (('islands.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool', lambda: islands.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool(m_among, method=kernel[0], **kernel[1]))
             if network and t % 2 else
             ('islands.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool', lambda: islands.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool(m_among, method=kernel[0], **kernel[1]))),
            ('islands.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool', lambda: islands.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool()),
            ('islands.meta_disturbance_process_in_patches', lambda: islands.meta_disturbance_process_in_patches(dist)),
            ('islands.meta_clear_up_offspring_marker_and_immigrant_marker_pool', islands.meta_clear_up_offspring_marker_and_immigrant_marker_pool)]


# ------------------------------------------------------------------ reporting
def finding(severity, func, text, evidence=None, fix=None):
    ''' severity: bug | risk | info | known; func = function(s) concerned; evidence = numbers / reproducer; fix = suggestion only (never applied) '''
    FINDINGS.append({'severity': severity, 'func': func, 'text': text, 'evidence': evidence, 'fix': fix})
    print('[FINDING/%s] %s: %s' % (severity, func, text))


def save(name, data):
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    data = dict(data); data['findings'] = FINDINGS
    with open(os.path.join(HERE, 'results', name + '.json'), 'w') as fh:
        json.dump(data, fh, indent=1, default=str)
