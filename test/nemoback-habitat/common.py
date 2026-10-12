# Shared helpers for test/nemoback-habitat. Every script runs in its own process (Nemo keeps process-wide static state: RAND, Individual::currentID).
import os, sys, ast, json, random, tempfile, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, REPO)
import metaibm                                                  # noqa: F401,E402  normal package import (metaibm/__init__ imports fine)
from metaibm.habitat import habitat                             # noqa: E402
from metaibm.indfactory import indfactory                       # noqa: E402

NEMO_CACHE = tempfile.mkdtemp(prefix='nemo_cache_')
POOLS = ('offspring_pool', 'immigrant_pool', 'dormancy_pool')
ORIGINAL_COMMIT = 'HEAD'                                        # pre-NemoBackend habitat.py (the migration is the working-tree change)
FINDINGS = []


def seed_all(seed):
    random.seed(seed); np.random.seed(seed)


def nemo_params(sp_log_name, init_value, mutation_rate=0.0001, env_variance=0.002025):
    ''' as test/nemo-indfactory/common.py: 2 quanti traits, 40 loci, continuous alleles '''
    from pynemo.params import SimParameters
    params = SimParameters()
    rows = [('replicates', 1), ('generations', 1), ('patch_number', 1), ('patch_capacity', 1), ('aging', 1),
            ('run_mode', 'overwrite'), ('root_dir', NEMO_CACHE), ('logfile', NEMO_CACHE + '/nemo.log'), ('filename', sp_log_name),
            ('quanti_traits', 2), ('quanti_loci', 40), ('quanti_allele_model', 'continuous'),
            ('quanti_init_model', 2), ('quanti_init_value', '{{%g, %g}}' % (init_value, init_value)), ('quanti_init_variance', '{{0.002025, 0.002025}}'),
            ('quanti_mutation_rate', mutation_rate), ('quanti_mutation_variance', 0.01)]
    if env_variance is not None: rows.append(('quanti_environmental_variance', '{{%g, %g}}' % (env_variance, env_variance)))
    for k, v in rows: params.add(k, v)
    return params


def make_species(name='sp1', init=0.2, seed=1, **kw):
    return indfactory(name, ['pheno1', 'pheno2'], nemo_params(name, init, **kw), seed=seed)


def make_hab(L=12, K=0, global_index=3, mean=(0.2, 0.2), cls=habitat, name='h1', W=None):
    return cls(name, 0, (0, 0), 2, ['env1', 'env2'], list(mean), [0.025, 0.025], L, L if W is None else W, dormancy_pool_max_size=K, global_index=global_index)


# ------------------------------------------------------------------ fake individuals (bookkeeping without Nemo)
RECYCLED = []


class fake_indi():
    ''' identity-compared stand-in for individual (no __eq__ / __hash__); recycle() twice raises like the real one '''
    def __init__(self, label, species_id='sp1', gender='female', pheno=(0.2, 0.2)):
        self.label, self.species_id, self.gender, self.pheno, self.is_alive = label, species_id, gender, list(pheno), True

    def get_indi_phenotype_ls(self):
        if not self.is_alive: raise RuntimeError('read after recycle')
        return list(self.pheno)

    def recycle(self):
        if not self.is_alive: raise RuntimeError('fake_indi %r recycled twice' % (self.label,))
        self.is_alive = False
        RECYCLED.append(self.label)

    def __repr__(self):
        return '<fake %s>' % (self.label,)


def fakes(prefix, n, species_id='sp1'):
    return [fake_indi('%s%d' % (prefix, i), species_id, 'female' if i % 2 == 0 else 'male') for i in range(n)]


def labels(hab):
    return [[None if x is None else x.label for x in row] for row in hab.set['microsite_individuals']]


# ------------------------------------------------------------------ holdings + invariants
def microsite_individuals(hab):
    return [x for row in hab.set['microsite_individuals'] for x in row if x is not None]


def held_individuals(hab):
    held = microsite_individuals(hab)
    for p in POOLS: held += list(getattr(hab, p))
    return held


def _as_list(x):
    if x is None: return []
    if isinstance(x, dict): return list(x.values())
    return list(x) if isinstance(x, (list, tuple)) else [x]


def check_invariants(hab, factories, tag=''):
    ''' (a) nothing held twice; (b) everything held is alive; (c) minted - recycled == held per factory (summed over habitats);
        (d) indi_num / occupied / empty consistent; (e) species_category matches the microsites '''
    all_held = []
    for h in _as_list(hab):
        all_held += held_individuals(h)
        occ = [(i, j) for i in range(h.length) for j in range(h.width) if h.set['microsite_individuals'][i][j] is not None]
        assert h.indi_num == len(occ) == len(h.occupied_site_pos_ls) and sorted(h.occupied_site_pos_ls) == occ, '%s (d) indi_num/occupied' % tag
        assert sorted(h.empty_site_pos_ls + h.occupied_site_pos_ls) == [(i, j) for i in range(h.length) for j in range(h.width)], '%s (d) empty+occupied != grid' % tag
        cat = []
        for sp_id, by_g in h.species_category.items():
            for g, pos_ls in by_g.items():
                for (i, j) in pos_ls:
                    x = h.set['microsite_individuals'][i][j]
                    assert x is not None and x.species_id == sp_id and x.gender == g, '%s (e) species_category at %r' % (tag, (i, j))
                    cat.append((i, j))
        assert sorted(cat) == occ, '%s (e) species_category != microsites' % tag
    assert len(set(map(id, all_held))) == len(all_held), '%s (a) an object is held twice' % tag
    assert all(x.is_alive for x in all_held), '%s (b) a held individual is recycled' % tag
    for f in _as_list(factories):
        n = sum(1 for x in all_held if getattr(x, 'indfactory', None) is f)
        assert f.minted_num - f.recycled_num == n, '%s (c) %s minted-recycled=%d held=%d' % (tag, f.species_id, f.minted_num - f.recycled_num, n)
    return True


def recycle_everything(habs):
    ''' end-of-run teardown with habitat's own clear-up functions '''
    return sum(h.hab_clear_up_microsite_individuals() + h.hab_clear_up_offspring_pool() + h.hab_clear_up_immigrant_pool() + h.hab_clear_up_dormancy_pool()
               for h in _as_list(habs))


# ------------------------------------------------------------------ original (pre-NemoBackend) code, read verbatim from git
def original_source():
    return subprocess.run(['git', 'show', ORIGINAL_COMMIT + ':metaibm/habitat.py'], capture_output=True, text=True, cwd=REPO).stdout


def original_habitat(names):
    ''' subclass of habitat whose `names` methods are the original bodies from git (compiled verbatim) '''
    import math, copy
    src = original_source(); lines = src.splitlines()
    fn = {n.name: n for n in ast.walk(ast.parse(src)) if isinstance(n, ast.FunctionDef)}
    body = 'class original_habitat(habitat):\n' + '\n'.join('\n'.join(lines[fn[k].lineno - 1:fn[k].end_lineno]) for k in names) + '\n'
    ns = {'habitat': habitat, 'random': random, 'np': np, 'math': math, 'copy': copy}
    exec(compile(body, '<git %s:metaibm/habitat.py>' % ORIGINAL_COMMIT, 'exec'), ns)
    return ns['original_habitat']


# ------------------------------------------------------------------ reporting
def finding(severity, text):
    FINDINGS.append({'severity': severity, 'text': text})
    print('[FINDING/%s] %s' % (severity, text))


def save(name, data):
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    data = dict(data); data['findings'] = FINDINGS
    with open(os.path.join(HERE, 'results', name + '.json'), 'w') as fh:
        json.dump(data, fh, indent=1, default=str)


def chi2_ok(counts, expected):
    ''' chi-square goodness of fit with a loose bound (df + 6*sqrt(2df)): catches structural bias, not 1-in-20 noise '''
    import math
    df = len(counts) - 1
    stat = sum((c - e) ** 2 / e for c, e in zip(counts, expected))
    return stat <= df + 6 * math.sqrt(2 * df), round(stat, 1)


def rss_mb():
    return int(subprocess.run(['ps', '-o', 'rss=', '-p', str(os.getpid())], capture_output=True, text=True).stdout) / 1024
