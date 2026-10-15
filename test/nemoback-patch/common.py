# Shared helpers for test/nemoback-patch. Every script runs in its own process (Nemo keeps process-wide static state: RAND, Individual::currentID).
# Reuses test/nemoback-habitat/common.py (loaded as `hab_common`, both files are called common.py) for fake individuals, Nemo params and invariants.
import os, sys, ast, json, random, hashlib, inspect, subprocess, importlib.util
import numpy as np

sys.dont_write_bytecode = True                                  # never write __pycache__ into the habitat suite's folder
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
HAB_SUITE = os.path.join(REPO, 'test', 'nemoback-habitat')

_spec = importlib.util.spec_from_file_location('hab_common', os.path.join(HAB_SUITE, 'common.py'))
hab_common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hab_common)                            # inserts REPO into sys.path, imports metaibm + indfactory

from metaibm.patch import patch                                 # noqa: E402
from metaibm.habitat import habitat                             # noqa: E402
from metaibm.indfactory import indfactory                       # noqa: E402
from metaibm.individual import individual                       # noqa: E402

fake_indi, make_species, seed_all, chi2_ok, rss_mb = hab_common.fake_indi, hab_common.make_species, hab_common.seed_all, hab_common.chi2_ok, hab_common.rss_mb
RECYCLED = hab_common.RECYCLED                                  # labels of recycled fake individuals (shared list object)
POOLS = ('offspring_pool', 'immigrant_pool', 'dormancy_pool')
MARKER_POOLS = ('offspring_marker_pool', 'immigrant_marker_pool')
FINDINGS = []


def git(*a):
    return subprocess.run(['git'] + list(a), capture_output=True, text=True, cwd=REPO).stdout


# ------------------------------------------------------------------ original (pre-NemoBackend) commit
def _find_original_commit():
    ''' baace34 is expected to hold the pre-NemoBackend patch.py (still has the combined clear-up); otherwise walk git log -- metaibm/patch.py '''
    if 'hab_clear_up_offspring_and_immigrant_pool' in git('show', 'baace34:metaibm/patch.py'): return 'baace34'
    for c in git('log', '--format=%h', '--', 'metaibm/patch.py').split():
        if 'hab_clear_up_offspring_and_immigrant_pool' in git('show', c + ':metaibm/patch.py'): return c
    raise RuntimeError('no commit with the pre-NemoBackend patch.py found')


ORIGINAL_COMMIT = _find_original_commit()


def original_source(module='patch'):
    return git('show', '%s:metaibm/%s.py' % (ORIGINAL_COMMIT, module))


def original_patch(names):
    ''' subclass of the current patch whose `names` methods are the original bodies from git (compiled verbatim); runs on the current (NemoBackend) habitat '''
    import re, copy
    src = original_source('patch'); lines = src.splitlines()
    fn = {n.name: n for n in ast.walk(ast.parse(src)) if isinstance(n, ast.FunctionDef)}
    body = 'class original_patch(patch):\n' + '\n'.join('\n'.join(lines[fn[k].lineno - 1:fn[k].end_lineno]) for k in names) + '\n'
    ns = {'patch': patch, 'habitat': habitat, 'random': random, 'np': np, 're': re, 'copy': copy}
    exec(compile(body, '<git %s:metaibm/patch.py>' % ORIGINAL_COMMIT, 'exec'), ns)
    return ns['original_patch']


_STACK = {}


def original_stack():
    ''' full original stack: (patch, habitat) classes compiled from the original commit's patch.py + habitat.py (pre-NemoBackend on both layers) '''
    if not _STACK:
        hsrc = original_source('habitat').replace('from .individual import individual', 'individual = _INDIVIDUAL')
        nh = {'__name__': 'orig_habitat', '_INDIVIDUAL': individual}
        exec(compile(hsrc, '<git %s:metaibm/habitat.py>' % ORIGINAL_COMMIT, 'exec'), nh)
        psrc = original_source('patch').replace('from .habitat import habitat', 'habitat = _ORIG_HABITAT')
        np_ = {'__name__': 'orig_patch', '_ORIG_HABITAT': nh['habitat']}
        exec(compile(psrc, '<git %s:metaibm/patch.py>' % ORIGINAL_COMMIT, 'exec'), np_)
        _STACK.update(patch=np_['patch'], habitat=nh['habitat'])
    return _STACK['patch'], _STACK['habitat']


# ------------------------------------------------------------------ building patches
def make_patch(n_habs=3, L=4, K=0, gidx0=0, means=None, W=None, cls=patch, name='p1', sizes=None):
    ''' patch with n_habs habitats h0..h{n-1} (hab_index 0..n-1, global_index gidx0+k); means[k] = mean_env_ls of habitat k (default (0.2, 0.2));
        sizes[k] = (L, W) overrides L / W per habitat. Seeds nothing: callers seed np.random first (habitat.__init__ draws the environment). '''
    p = cls(name, 0, (0, 0))
    takes_gi = 'global_index' in inspect.signature(cls.add_habitat).parameters
    for k in range(n_habs):
        mean = list(means[k]) if means is not None else [0.2, 0.2]
        l, w = sizes[k] if sizes is not None else (L, L if W is None else W)
        kw = dict(global_index=gidx0 + k) if takes_gi else {}
        p.add_habitat('h%d' % k, k, (0, k), 2, ['env1', 'env2'], mean, [0.025, 0.025], l, w, K, **kw)
    return p


def fill_fakes(h, tag=None, occ=0, O=0, I=0, D=0, markers=(0, 0), species=lambda i: 'sp%d' % (1 + i % 2), positions=None):
    ''' fill habitat h with fake individuals: occ microsites (first occ sites row-major, or `positions`), O/I/D pool members, marker pools '''
    tag = h.name if tag is None else tag
    pos = positions if positions is not None else [(k // h.width, k % h.width) for k in range(occ)]
    for k, (i, j) in enumerate(pos):
        h.add_individual(fake_indi('%s.m%d' % (tag, k), species(k), 'female' if k % 3 else 'male'), i, j)
    h.offspring_pool = [fake_indi('%s.o%d' % (tag, i), species(i)) for i in range(O)]
    h.immigrant_pool = [fake_indi('%s.i%d' % (tag, i), species(i)) for i in range(I)]
    h.dormancy_pool = [fake_indi('%s.d%d' % (tag, i), species(i)) for i in range(D)]
    h.offspring_marker_pool = [('p1', h.name, 'sexual')] * markers[0]
    h.immigrant_marker_pool = [('p0', 'hX', 'sexual')] * markers[1]
    return h


def clear_recycled():
    del RECYCLED[:]


# ------------------------------------------------------------------ holdings, snapshots, invariants
def habs(p):
    return [p.set[n] for n in p.hab_id_ls]


def _as_patches(ps):
    return ps if isinstance(ps, (list, tuple)) else [ps]


def held(ps):
    out = []
    for p in _as_patches(ps):
        for h in habs(p): out += hab_common.held_individuals(h)
    return out


def _lab(x):
    return x.label if hasattr(x, 'label') else repr(x)


def snapshot(p):
    ''' full bookkeeping state of a patch of fake individuals (labels, list orders included) '''
    out = []
    for h in habs(p):
        cat = {s: {g: list(v) for g, v in d.items()} for s, d in h.species_category.items()}
        out.append((h.name, hab_common.labels(h), tuple(tuple(_lab(x) for x in getattr(h, q)) for q in POOLS), tuple(tuple(getattr(h, q)) for q in MARKER_POOLS),
                    tuple(h.empty_site_pos_ls), tuple(h.occupied_site_pos_ls), h.indi_num, repr(cat)))
    return out


def held_labels(p):
    return {_lab(x) for x in held(p)}


def nemo_digest(ps):
    ''' sha1 of the state of real Nemo individuals: per habitat the grid (ID, gender, home, phenotype) and each pool (IDs) '''
    rows = []
    for p in _as_patches(ps):
        for h in habs(p):
            grid = [None if x is None else (x.individual_id, x.gender, x.nemo_indi_obj.getHome(), [round(v, 12) for v in x.get_indi_phenotype_ls()])
                    for row in h.set['microsite_individuals'] for x in row]
            rows.append((p.name, h.name, grid, [[x.individual_id for x in getattr(h, q)] for q in POOLS], [len(getattr(h, q)) for q in MARKER_POOLS]))
    return hashlib.sha1(json.dumps(rows).encode()).hexdigest()[:16]


def check_invariants(ps, factories, tag=''):
    ''' habitat invariants (via hab_common) over every habitat of every patch + patch bookkeeping (hab_num / hab_id_ls / set) '''
    all_h = []
    for p in _as_patches(ps):
        assert p.hab_num == len(p.set) == len(p.hab_id_ls) and list(p.set) == p.hab_id_ls, '%s patch bookkeeping' % tag
        all_h += habs(p)
    return hab_common.check_invariants(all_h, factories, tag)


def factories_balanced(ps, factories, tag=''):
    hs = held(ps)
    for f in factories:
        n = sum(1 for x in hs if getattr(x, 'indfactory', None) is f)
        assert f.minted_num - f.recycled_num == n, '%s %s minted-recycled=%d held=%d' % (tag, f.species_id, f.minted_num - f.recycled_num, n)
    return True


def recycle_everything(ps):
    ''' teardown with the patch's own clear-up functions; returns the number recycled '''
    return sum(p.patch_clear_up_microsite_individuals() + p.patch_clear_up_offspring_pool() + p.patch_clear_up_immigrant_pool() + p.patch_clear_up_dormancy_pool()
               for p in _as_patches(ps))


# ------------------------------------------------------------------ RNG + matrix forcing
def rng_state():
    s = np.random.get_state()
    return random.getstate(), (s[0], s[1].copy(), s[2], s[3], s[4])


def same_rng_state(a, b):
    return a[0] == b[0] and a[1][0] == b[1][0] and np.array_equal(a[1][1], b[1][1]) and a[1][2:] == b[1][2:]


def force_within_matrix(p, M):
    ''' monkeypatch (instance only) the within-dispersal number matrix to M: get_dispersal_within_num_matrix (direct-germination variant)
        and patch_matrix_around (immigrant-pool / marker variants); no np.random is used by the forced versions '''
    M = np.matrix(np.array(M, dtype=float))
    p.get_dispersal_within_num_matrix = lambda disp_within_rate: M.copy()
    p.patch_matrix_around = lambda matrix: M.copy()
    return p


# ------------------------------------------------------------------ Nemo species
def nemo_species(name='sp1', init=0.2, seed=1, **kw):
    return make_species(name, init, seed=seed, **kw)


# ------------------------------------------------------------------ reporting
def finding(severity, text, repro=None, fix=None):
    ''' severity: bug | risk | info | known '''
    FINDINGS.append({'severity': severity, 'text': text, 'repro': repro, 'fix': fix})
    print('[FINDING/%s] %s' % (severity, text))


def save(name, data):
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    data = dict(data); data['findings'] = FINDINGS
    with open(os.path.join(HERE, 'results', name + '.json'), 'w') as fh:
        json.dump(data, fh, indent=1, default=str)


def binom_ok(k, n, p, z=4.0):
    ''' |k/n - p| within z standard errors '''
    import math
    return abs(k / n - p) <= z * math.sqrt(p * (1 - p) / n) + 1e-12
