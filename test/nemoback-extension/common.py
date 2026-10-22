# Shared helpers for test/nemoback-extension (extension/global_habitat_network.py, sexual). Every script runs in its own process (Nemo keeps process-wide static state).
# Reuses test/nemoback-metacommunity/common.py (loaded via importlib as `meta_common`; through it the patch and habitat suites' helpers).
# Nothing in metaibm/ or extension/ is edited; original (pre-migration) bodies are compiled from git into a subclass of the CURRENT metacommunity.
import os, sys, ast, json, random, hashlib, subprocess, importlib.util
import numpy as np

sys.dont_write_bytecode = True                                  # never write __pycache__ into the other suites' folders
os.environ.setdefault('MPLBACKEND', 'Agg')                      # metacommunity.py imports matplotlib.pyplot
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
META_SUITE = os.path.join(REPO, 'test', 'nemoback-metacommunity')
EXT_FILE = 'extension/global_habitat_network.py'
META_FILE = 'metaibm/metacommunity.py'


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


meta_common = _load('meta_common', os.path.join(META_SUITE, 'common.py'))         # inserts REPO into sys.path, imports metaibm, loads patch_common / hab_common
patch_common, hab_common = meta_common.patch_common, meta_common.hab_common

from metaibm.metacommunity import metacommunity                 # noqa: E402
from metaibm.patch import patch                                 # noqa: E402
import extension.global_habitat_network as ext                  # noqa: E402

species, seed_all, check_world, teardown, digest, rng_digest = meta_common.species, meta_common.seed_all, meta_common.check_world, meta_common.teardown, meta_common.digest, meta_common.rng_digest
rng_state, same_rng_state, full_snapshot, same_snapshot = meta_common.rng_state, meta_common.same_rng_state, meta_common.full_snapshot, meta_common.same_snapshot
held, held_ids, all_habs, minted, recycled = meta_common.held, meta_common.held_ids, meta_common.all_habs, meta_common.minted, meta_common.recycled
POOLS, MARKER_POOLS, TABLE = meta_common.POOLS, meta_common.MARKER_POOLS, meta_common.TABLE
make_islands, make_mainland, rain_for, fill_pools, indfactory = meta_common.make_islands, meta_common.make_mainland, meta_common.rain_for, meta_common.fill_pools, meta_common.indfactory
KERNELS = [('uniform', {}), ('exponential', {'rho': 2}), ('gaussian', {'sigma': 1.5}), ('cauchy', {'gamma': 1.0}), ('power_law', {'alpha': 2.0, 'r0': 1.0})]
FINDINGS = []


def git(*a):
    return subprocess.run(['git'] + list(a), capture_output=True, text=True, cwd=REPO).stdout


# ------------------------------------------------------------------ originals from git
ORIGINAL_EXT_COMMIT = git('log', '-1', '--format=%h', '--', EXT_FILE).strip()       # last commit that touched the extension (v3 version; 1356387 when written)
ORIGINAL_META_COMMIT = git('rev-parse', '--short', 'HEAD').strip()                  # add_patch as committed (registry call + separate distance-matrix call)
ORIGINAL_REGISTRY = 'global_habitat_id_idx_registry'


def _functions(src, cls_name=None):
    body = ast.parse(src).body
    if cls_name is not None: body = [n for n in body if isinstance(n, ast.ClassDef) and n.name == cls_name][0].body
    return {n.name: n for n in body if isinstance(n, ast.FunctionDef)}


def _cut(src, node, indent=''):
    ''' source of one def; dedented to module level, re-indented by `indent` '''
    lines = src.splitlines()[node.lineno - 1:node.end_lineno]
    pad = len(lines[0]) - len(lines[0].lstrip())
    return '\n'.join(indent + l[pad:] if l.strip() else '' for l in lines)


def current_ext_source(): return open(os.path.join(REPO, EXT_FILE)).read()
def original_ext_source(): return git('show', '%s:%s' % (ORIGINAL_EXT_COMMIT, EXT_FILE))


# ------------------------------------------------------------------ the {(patch_id, habitat_id): global_index} record: name read from the code (may be renamed)
GIDX_DICT_NAMES = ('global_habitat_id_2_index_dir', 'global_habitat_id_2_gidx_dir')


def gidx_dict_name():
    ''' the one dict attribute that _init_global_habitat_network_fields creates (global_habitat_id_2_index_dir today, global_habitat_id_2_gidx_dir after a rename) '''
    fn = _functions(current_ext_source())['_init_global_habitat_network_fields']
    names = [t.attr for n in ast.walk(fn) if isinstance(n, ast.Assign) and isinstance(n.value, ast.Dict) for t in n.targets if isinstance(t, ast.Attribute)]
    assert len(names) == 1, 'expected one dict field in _init_global_habitat_network_fields, found %s' % names
    return names[0]


GIDX_DICT = gidx_dict_name()
def gidx_dict(m): return getattr(m, GIDX_DICT)
def original_meta_source(): return git('show', '%s:%s' % (ORIGINAL_META_COMMIT, META_FILE))


def original_metacommunity(ext_names=(), meta_names=(), patches=None):
    ''' subclass of the CURRENT metacommunity whose ext_names (module-level functions of the extension) and meta_names (metacommunity methods) are the original
        bodies from git, compiled verbatim. patches = {name: [(old, new), ...]} text substitutions applied to a body before compiling (each must match exactly once). '''
    ext_src, meta_src = original_ext_source(), original_meta_source()
    ext_fn, meta_fn = _functions(ext_src), _functions(meta_src, 'metacommunity')
    parts = []
    for srcs, fns, names in ((ext_src, ext_fn, ext_names), (meta_src, meta_fn, meta_names)):
        for k in names:
            body = _cut(srcs, fns[k], '    ')
            for old, new in (patches or {}).get(k, []):
                assert body.count(old) == 1, 'patch for %s matches %d times: %r' % (k, body.count(old), old)
                body = body.replace(old, new)
            parts.append(body)
    code = 'class original_metacommunity(metacommunity):\n' + '\n\n'.join(parts) + '\n'
    ns = {'metacommunity': metacommunity, 'patch': patch, 'random': random, 'np': np}
    exec(compile(code, '<git %s:%s + %s:%s>' % (ORIGINAL_EXT_COMMIT, EXT_FILE, ORIGINAL_META_COMMIT, META_FILE), 'exec'), ns)
    cls = ns['original_metacommunity']
    if GIDX_DICT != 'global_habitat_id_2_index_dir':            # the v3 registry writes the v3 name; after a rename the current __init__ no longer creates it
        cls.global_habitat_id_2_index_dir = property(lambda self: self.__dict__.setdefault('_v3_global_habitat_id_2_index_dir', {}))
    return cls


# ------------------------------------------------------------------ worlds
def build_world(spec, cls=metacommunity, seed=0, name='islands'):
    ''' spec = [(patch_name, patch_index, (x, y), [(hab_name, hab_index, (x, y), global_index), ...]), ...] added in list order.
        Seeds both RNGs first (habitat.__init__ draws the environment); returns (metacommunity, per-add_patch counts of update_global_habitat_distance_matrix calls) '''
    seed_all(seed)
    m = cls(name)
    calls, per_add = [0], []
    real = m.update_global_habitat_distance_matrix
    def counted(*a, **k):
        calls[0] += 1; return real(*a, **k)
    m.update_global_habitat_distance_matrix = counted                       # instance attribute: counts calls from registry or add_patch alike
    for pn, pidx, ploc, habs in spec:
        p = patch(pn, pidx, ploc)
        for hn, hidx, hloc, gidx in habs:
            mean = list(meta_common.SPECIES_MEANS[(pidx + hidx) % 2])
            p.add_habitat(hn, hidx, hloc, 2, ['env1', 'env2'], mean, [0.025, 0.025], 3, 3, 5, global_index=gidx)
        before = calls[0]; m.add_patch(pn, p); per_add.append(calls[0] - before)
    del m.update_global_habitat_distance_matrix
    return m, per_add


def random_spec(rng, n_patches=None):
    ''' random world: unsorted patch / habitat names, 1..4 habitats per patch, random locations, shuffled non-contiguous global_index in 0..65535 '''
    n_patches = n_patches or rng.randint(2, 5)
    sizes = [rng.randint(1, 4) for _ in range(n_patches)]
    gidx = rng.sample(range(65536), sum(sizes)); k = 0
    pnames = rng.sample(['pZ', 'pA', 'p10', 'p2', 'pm', 'pB', 'q'], n_patches)
    spec = []
    for pi, (pn, n) in enumerate(zip(pnames, sizes)):
        px, py = rng.uniform(-10, 10), rng.uniform(-10, 10)
        habs = []
        for hi, hn in enumerate(rng.sample(['hz', 'h0', 'ha', 'h10', 'h2', 'x'], n)):
            habs.append((hn, hi, (px + rng.uniform(-1, 1), py + rng.uniform(-1, 1)), gidx[k])); k += 1
        spec.append((pn, pi, (px, py), habs))
    return spec


def release(metas):
    ''' teardown without the factory-balance assertion: recycle everything one metacommunity holds while others still hold individuals of the same factories '''
    n = 0
    for m in (metas if isinstance(metas, (list, tuple)) else [metas]):
        n += m.meta_clear_up_microsite_individuals() + m.meta_clear_up_offspring_pool() + m.meta_clear_up_immigrant_pool() + m.meta_clear_up_dormancy_pool()
        m.meta_clear_up_offspring_marker_and_immigrant_marker_pool()
    return n


def brute_distance(m):
    locs = [m.set[a].set[b].location for a, b in m.global_habitat_id_ls]
    return np.array([[np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2) for (x2, y2) in locs] for (x1, y1) in locs])


def md5(a): return hashlib.md5(np.ascontiguousarray(np.asarray(a)).tobytes()).hexdigest()[:12]


# ------------------------------------------------------------------ schedules (one step, as (name, callable) pairs; shared by smoke / chains / fuzz)
OFF = 'dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool'
MARK = 'dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool'


def pool_ops(ml, isl, rain, m_among, kern, kw, dist=0.05, net_call=None):
    ''' individual-pool schedule (v3 order) with the network offspring dispersal; net_call(isl) overrides the network call (reference classes) '''
    net = net_call or (lambda: getattr(isl, OFF)(m_among, method=kern, **kw))
    return [('mainland.meta_dead_selection', lambda: ml.meta_dead_selection(0.1, 0.5, 'niche_gaussian')),
            ('mainland.meta_mainland_sexual_birth_mutate_germinate', lambda: ml.meta_mainland_sexual_birth_mutate_germinate(1.0, None, None)),
            ('islands.meta_colonize_from_propagules_rains', lambda: isl.meta_colonize_from_propagules_rains(ml, rain)),
            ('islands.meta_dead_selection', lambda: isl.meta_dead_selection(0.1, 0.5, 'niche_gaussian')),
            ('islands.meta_sex_reproduce_mutate_into_offspring_pool', lambda: isl.meta_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)),
            ('islands.meta_dispersal_within_patch_from_offspring_to_immigrant_pool', lambda: isl.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(0.3)),
            (OFF, net),
            ('islands.meta_local_germinate_from_offspring_and_immigrant_pool', lambda: isl.meta_local_germinate_from_offspring_and_immigrant_pool()),
            ('islands.meta_dormancy_process_from_offspring_pool_and_immigrant_pool', lambda: isl.meta_dormancy_process_from_offspring_pool_and_immigrant_pool()),
            ('islands.meta_disturbance_process_in_habitat', lambda: isl.meta_disturbance_process_in_habitat(dist)),
            ('islands.meta_clear_up_offspring_and_immigrant_pool', lambda: isl.meta_clear_up_offspring_and_immigrant_pool())]


def marker_ops(ml, isl, rain, m_among, kern, kw, is_remove=True, dist=0.02, net_call=None):
    ''' marker schedule (experiments/model-sloss-GRFE.py order, sexual) with the network marker dispersal '''
    net = net_call or (lambda: getattr(isl, MARK)(m_among, method=kern, is_remove=is_remove, **kw))
    return [('mainland.meta_dead_selection', lambda: ml.meta_dead_selection(0.1, 0.5, 'niche_gaussian')),
            ('mainland.meta_mainland_sexual_birth_mutate_germinate', lambda: ml.meta_mainland_sexual_birth_mutate_germinate(1.0, None, None)),
            ('islands.meta_dead_selection', lambda: isl.meta_dead_selection(0.1, 0.5, 'niche_gaussian')),
            ('islands.meta_sex_reproduce_calculation_with_offspring_marker_pool', lambda: isl.meta_sex_reproduce_calculation_with_offspring_marker_pool(0.8)),
            ('islands.meta_colonize_from_propagules_rains', lambda: isl.meta_colonize_from_propagules_rains(ml, rain)),
            ('islands.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool', lambda: isl.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(0.3)),
            (MARK, net),
            ('islands.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool', lambda: isl.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None)),
            ('islands.meta_disturbance_process_in_patches', lambda: isl.meta_disturbance_process_in_patches(dist)),
            ('islands.meta_clear_up_offspring_marker_and_immigrant_marker_pool', lambda: isl.meta_clear_up_offspring_marker_and_immigrant_marker_pool())]


# ------------------------------------------------------------------ statistics (no scipy in the test environment)
import math


def _gammaincc(a, x):
    ''' regularized upper incomplete gamma Q(a, x) (Numerical Recipes: series for x < a + 1, continued fraction otherwise) '''
    if x <= 0: return 1.0
    if x < a + 1:
        ap, s, d = a, 1.0 / a, 1.0 / a
        for _ in range(10000):
            ap += 1; d *= x / ap; s += d
            if abs(d) < abs(s) * 1e-15: break
        return max(0.0, 1.0 - s * math.exp(-x + a * math.log(x) - math.lgamma(a)))
    b, c, d = x + 1 - a, 1e300, 1.0 / (x + 1 - a); h = d
    for i in range(1, 10000):
        an = -i * (i - a); b += 2
        d = an * d + b; d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c; c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d; de = d * c; h *= de
        if abs(de - 1) < 1e-15: break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def chi2_sf(x, df): return _gammaincc(df / 2.0, x / 2.0) if df > 0 else 1.0
def norm_sf(z): return 0.5 * math.erfc(z / math.sqrt(2))


def chi2_gof(obs, exp, min_exp=5.0):
    ''' goodness of fit; bins with expectation < min_exp are pooled into one; returns (stat, df, p) '''
    o, e, ro, re_ = [], [], 0.0, 0.0
    for a, b in zip(obs, exp):
        if b >= min_exp: o.append(a); e.append(b)
        else: ro += a; re_ += b
    if re_ > 0: o.append(ro); e.append(re_)
    tot_o, tot_e = sum(o), sum(e)
    e = [x * tot_o / tot_e for x in e] if tot_e else e                # compare shapes (totals may differ by the capping)
    stat = sum((a - b) ** 2 / b for a, b in zip(o, e) if b > 0)
    df = len(o) - 1
    return stat, df, chi2_sf(stat, df)


def welch_p(xs, ys):
    ''' two-sided p of the difference of means (normal approximation, large samples) '''
    nx, ny = len(xs), len(ys)
    mx, my = sum(xs) / nx, sum(ys) / ny
    vx = sum((v - mx) ** 2 for v in xs) / max(nx - 1, 1); vy = sum((v - my) ** 2 for v in ys) / max(ny - 1, 1)
    se = math.sqrt(vx / nx + vy / ny)
    return 1.0 if se == 0 else 2 * norm_sf(abs(mx - my) / se)


# ------------------------------------------------------------------ reporting
def finding(severity, func, text, evidence=None, fix=None):
    ''' severity: bug | risk | info | known | todo (accepted); fix = suggestion only (never applied) '''
    FINDINGS.append({'severity': severity, 'func': func, 'text': text, 'evidence': evidence, 'fix': fix})
    print('[FINDING/%s] %s: %s' % (severity, func, text))


def save(name, data):
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    data = dict(data); data['findings'] = FINDINGS
    with open(os.path.join(HERE, 'results', name + '.json'), 'w') as fh:
        json.dump(data, fh, indent=1, default=str)
