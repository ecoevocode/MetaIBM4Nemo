# Tests of the make_copy / make_clone migration below the simulation level. One file, three sections; every test run is its own process (Nemo keeps
# process-wide static state: RAND, Individual::currentID):
#   helpers  shared by the two runs below, by test_simulation.py and by run_all.py (all use them as C.<name>)
#   unit     static checks, make_copy / make_clone semantics, guards, lifecycle, rng, bench                        -> RESULTS/test_unit.json
#   core     the 10 is_remove=False copy sites (A), births (B + C) through habitat / patch / metacommunity         -> RESULTS/test_core.json
# RESULTS = /tmp/nemoback-make-clone-copy/results (outside the repo); run_all.py turns them into TEST_REPORT.md.
# usage: python test_core.py unit [static] [copy] [clone] [guards] [lifecycle] [rng] [bench]      (default: all parts)
#        python test_core.py core [copy] [births]                                                  (default: both)
#        python test_core.py --worker {native|imitation} SEED                                      (rng worker of unit / rng)
import os, re, io, sys, ast, copy, json, math, time, random, hashlib, platform, textwrap, tokenize, warnings, traceback, subprocess, importlib.util
from collections import Counter
import numpy as np


# ###################################################################### helpers
# Shared helpers for test/nemoback-make-clone-copy: individual.make_copy (Nemo Individual::operator=) / indfactory.make_clone (IndFactory::makeClonedOffsprg)
# and the core migration that uses them (is_remove=False copies, asexual / mixed births, marker births).
# Every script runs in its own process (Nemo keeps process-wide static state: RAND, Individual::currentID).
# Reuses test/nemoback-metacommunity/common.py (`mc`) and, through it, test/nemoback-patch/common.py (`pc`) and test/nemoback-habitat/common.py (`hc`), all loaded
# with importlib (every one is called common.py). Nothing outside this folder is edited: v3 method bodies are compiled from git into subclasses of the current classes.

sys.dont_write_bytecode = True                                  # never write __pycache__ (here or in the other suites' folders)
os.environ.setdefault('MPLBACKEND', 'Agg')                      # metacommunity.py imports matplotlib.pyplot
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
RESULTS = '/tmp/nemoback-make-clone-copy/results'                # outside the repo
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mc = _load('mc_common', os.path.join(REPO, 'test', 'nemoback-metacommunity', 'common.py'))   # loads patch_common + hab_common itself
pc = mc.patch_common
hc = mc.hab_common
assert os.path.samefile(hc.__file__, os.path.join(REPO, 'test', 'nemoback-habitat', 'common.py'))

from metaibm.habitat import habitat                             # noqa: E402
from metaibm.patch import patch                                 # noqa: E402
from metaibm.metacommunity import metacommunity                 # noqa: E402
from metaibm.indfactory import indfactory                       # noqa: E402
from metaibm.individual import individual                       # noqa: E402

POOLS, MARKER_POOLS = mc.POOLS, mc.MARKER_POOLS
seed_all, rng_state, same_rng_state = hc.seed_all, pc.rng_state, pc.same_rng_state
RECYCLED = hc.RECYCLED                                          # labels of recycled fake individuals (one list shared by hc / pc / mc)


def git(*a):
    return subprocess.run(['git'] + list(a), capture_output=True, text=True, cwd=REPO).stdout


def md5(rel):
    return hashlib.md5(open(os.path.join(REPO, rel), 'rb').read()).hexdigest()


CHANGED = git('diff', '--name-only').split()                    # the uncommitted migration (6 files expected)
HEAD = git('rev-parse', '--short', 'HEAD').strip()


# ------------------------------------------------------------------ v3 bodies compiled from git
V3_COMMIT = {'metaibm/habitat.py': 'baace34', 'metaibm/patch.py': 'baace34', 'metaibm/metacommunity.py': '7a1adbb', 'extension/global_habitat_network.py': '1356387'}


def v3_source(path):
    src = git('show', '%s:%s' % (V3_COMMIT[path], path))
    if not src: raise RuntimeError('git show %s:%s is empty' % (V3_COMMIT[path], path))
    return src


def v3_subclass(base, path, names):
    ''' subclass of the CURRENT class `base` whose `names` methods are the v3 bodies of `path` (verbatim from git, commit V3_COMMIT[path]);
        extension/*.py holds module-level functions: they are indented into the class body. The v3 bodies run on the current (NemoBackend) lower layers. '''
    src = v3_source(path); lines = src.splitlines(); tree = ast.parse(src)
    if path.startswith('extension/'):
        fn = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
        block = lambda n: textwrap.indent('\n'.join(lines[n.lineno - 1:n.end_lineno]), '    ')
    else:
        cls = [n for n in tree.body if isinstance(n, ast.ClassDef)][0]
        fn = {n.name: n for n in cls.body if isinstance(n, ast.FunctionDef)}
        block = lambda n: '\n'.join(lines[n.lineno - 1:n.end_lineno])
    missing = [k for k in names if k not in fn]
    if missing: raise KeyError('%s:%s has no %s' % (V3_COMMIT[path], path, missing))
    name = 'v3_' + base.__name__
    body = 'class %s(_BASE):\n' % name + '\n'.join(block(fn[k]) for k in names) + '\n'
    ns = {'_BASE': base, base.__name__: base, 'random': random, 'np': np, 'math': math, 'copy': copy}
    exec(compile(body, '<git %s:%s>' % (V3_COMMIT[path], path), 'exec'), ns)
    return ns[name]


# ------------------------------------------------------------------ fake individuals: make_copy()
FAKE_COPIES = []                                                # every fake copy made (for counting copies per call)


def _fake_make_copy(self):
    ''' stand-in for individual.make_copy(): a NEW fake_indi labelled label + '*', .origin = the label of the root original '''
    if not self.is_alive: raise RuntimeError('make_copy of recycled fake_indi %r' % (self.label,))
    c = hc.fake_indi(self.label + '*', self.species_id, self.gender, self.pheno)
    c.origin = getattr(self, 'origin', self.label)
    FAKE_COPIES.append(c)
    return c


hc.fake_indi.make_copy = _fake_make_copy


def origin(x):
    return getattr(x, 'origin', x.label)


def is_copy(x):
    return x.label.endswith('*')


# ------------------------------------------------------------------ worlds = habitat | patch | metacommunity | list of them
def habs_of(world):
    ''' [(key, habitat)] in a fixed order; key = (metacommunity name, patch name, habitat name) as far as known '''
    out = []
    for w in (world if isinstance(world, (list, tuple)) else [world]):
        if isinstance(w, habitat): out.append(((w.name,), w))
        elif isinstance(w, patch): out += [((w.name, hn), w.set[hn]) for hn in w.hab_id_ls]
        elif isinstance(w, metacommunity): out += [((w.metacommunity_name, pn, hn), w.set[pn].set[hn]) for pn in w.patch_id_ls for hn in w.set[pn].hab_id_ls]
        else: raise TypeError(type(w))
    return out


def layout(world):
    ''' bookkeeping of a world of fake individuals with every individual replaced by its origin label: grid, the three pools (in order), marker pools,
        empty / occupied lists (order matters for later RNG use), indi_num, species_category '''
    o = lambda x: None if x is None else origin(x)
    return [(key, [[o(x) for x in row] for row in h.set['microsite_individuals']], [[o(x) for x in getattr(h, q)] for q in POOLS],
             [list(getattr(h, q)) for q in MARKER_POOLS], list(h.empty_site_pos_ls), list(h.occupied_site_pos_ls), h.indi_num,
             repr({s: {g: list(v) for g, v in d.items()} for s, d in h.species_category.items()}))
            for key, h in habs_of(world)]


def pool_ids(world, pools=POOLS):
    ''' {(key, pool): [id(x), ...]} -- object identity of every pool member, in order '''
    return {(key, q): [id(x) for x in getattr(h, q)] for key, h in habs_of(world) for q in pools}


def grid_ids(world):
    return {(key, i, j): id(x) for key, h in habs_of(world) for i, row in enumerate(h.set['microsite_individuals']) for j, x in enumerate(row) if x is not None}


def individuals(world):
    ''' every held individual: microsites (row-major) + pools, per habitat '''
    out = []
    for _, h in habs_of(world):
        out += [x for row in h.set['microsite_individuals'] for x in row if x is not None]
        for q in POOLS: out += list(getattr(h, q))
    return out


def locations(world):
    ''' {id(x): where} for every held individual: ('grid', key, i, j) or ('pool', key, pool, index) '''
    loc = {}
    for key, h in habs_of(world):
        for i, row in enumerate(h.set['microsite_individuals']):
            for j, x in enumerate(row):
                if x is not None: loc[id(x)] = ('grid', key, i, j)
        for q in POOLS:
            for k, x in enumerate(getattr(h, q)): loc[id(x)] = ('pool', key, q, k)
    return loc


def copied_state(x):
    ''' everything Individual::operator= copies (all but the ID) + MetaIBM species_id / age: sex, Nemo age, mother / father IDs, home, pedigree class,
        fecundity, mating + realized-fecundity counters (5 classes), and every allele of every quanti axis (both strands). Phenotype excluded: redrawn when Ve > 0. '''
    n = x.nemo_indi_obj
    return (x.species_id, x.gender, x.age, n.getAge(), n.getMotherID(), n.getFatherID(), n.getHome(), n.getPedigreeClass(), n.getFecundity(),
            tuple(n.getMatings(c) for c in range(5)), tuple(n.getRealizedFecundity(c) for c in range(5)),
            tuple(tuple(x.genotype_table(d)) for d in range(x.traits_num)))


def species_params(name, init=0.2, traits=2, loci=40, mutation_rate=0.0001, env_variance=0.002025, mutation_variance=0.01):
    ''' as hc.nemo_params (continuous quanti alleles, founders at init) with configurable traits / loci (loci = TOTAL over the traits, Nemo's quanti_loci);
        env_variance None = no quanti_environmental_variance parameter at all '''
    from pynemo.params import SimParameters
    fmt = lambda v: '{{%s}}' % ', '.join('%g' % x for x in v)
    rows = [('replicates', 1), ('generations', 1), ('patch_number', 1), ('patch_capacity', 1), ('aging', 1),
            ('run_mode', 'overwrite'), ('root_dir', hc.NEMO_CACHE), ('logfile', hc.NEMO_CACHE + '/nemo.log'), ('filename', name),
            ('quanti_traits', traits), ('quanti_loci', loci), ('quanti_allele_model', 'continuous'),
            ('quanti_init_model', 2), ('quanti_init_value', fmt([init] * traits)), ('quanti_init_variance', fmt([0.002025] * traits)),
            ('quanti_mutation_rate', mutation_rate), ('quanti_mutation_variance', mutation_variance)]
    if env_variance is not None: rows.append(('quanti_environmental_variance', fmt([env_variance] * traits)))
    p = SimParameters()
    for k, v in rows: p.add(k, v)
    return p


def make_species(name='sp1', init=0.2, seed=1, traits=2, **kw):
    return indfactory(name, ['pheno%d' % (k + 1) for k in range(traits)], species_params(name, init, traits=traits, **kw), seed=seed)


def genotypic_values(x):
    ''' additive genotypic value per quanti axis = sum of both allele values over the loci of that axis (Ve = 0 -> phenotype) '''
    return [sum(a + b for a, b in x.genotype_table(d)) for d in range(x.traits_num)]


def alleles(x):
    return tuple(tuple(x.genotype_table(d)) for d in range(x.traits_num))


def check_held(world, factories, tag=''):
    ''' minted - recycled == individuals held by the world (per factory), every held one alive, nothing held twice; raises AssertionError '''
    held = individuals(world)
    assert len({id(x) for x in held}) == len(held), '%s an individual is held twice' % tag
    assert all(x.is_alive for x in held), '%s a held individual is recycled' % tag
    for f in (factories.values() if isinstance(factories, dict) else factories):
        n = sum(1 for x in held if getattr(x, 'indfactory', None) is f)
        assert f.minted_num - f.recycled_num == n, '%s %s minted - recycled = %d, held %d' % (tag, f.species_id, f.minted_num - f.recycled_num, n)
    return True


class CopyTracker():
    ''' records every individual.make_copy() / indfactory.make_clone() result while active (wraps the class methods, calls the originals unchanged);
        states=True also records (id(source), copied_state(source), copied_state(copy)) at the moment of each copy '''
    def __init__(self, states=False):
        self.copies, self.clones, self.states, self._with_states = [], [], [], states

    def __enter__(self):
        self._mc, self._cl = individual.make_copy, indfactory.make_clone
        tr = self
        def make_copy(indi_self):
            c = tr._mc(indi_self); tr.copies.append((indi_self, c))
            if tr._with_states: tr.states.append((id(indi_self), copied_state(indi_self), copied_state(c)))
            return c
        def make_clone(fac_self, mother_indi_obj, home_habitat=0):
            c = tr._cl(fac_self, mother_indi_obj, home_habitat); tr.clones.append((mother_indi_obj, c)); return c
        individual.make_copy, indfactory.make_clone = make_copy, make_clone
        return self

    def __exit__(self, *exc):
        individual.make_copy, indfactory.make_clone = self._mc, self._cl
        return False


# ------------------------------------------------------------------ checks, findings, results
CHECKS, FINDINGS = [], []
_T0 = time.time()
PART = ['']                                                     # current part of a run (unit: static / copy / clone / ...), recorded with every check


def part(name):
    PART[0] = name
    print('\n#### %s' % name, flush=True)


def ok(name, cond, evidence=None):
    ''' record one check; evidence = numbers / first failing seed (kept for PASS too) '''
    CHECKS.append({'name': name, 'ok': bool(cond), 'evidence': evidence, 'part': PART[0]})
    print('%s  %s%s' % ('PASS' if cond else 'FAIL', name, '' if evidence in (None, '') else '  -- %s' % (evidence,)), flush=True)
    return bool(cond)


def guarded(name, fn, *a, **kw):
    ''' run fn; an exception is a FAIL of check `name` (with the traceback tail as evidence) instead of aborting the script '''
    try: return fn(*a, **kw)
    except Exception as e:
        ok(name, False, '%s: %s | %s' % (type(e).__name__, e, traceback.format_exc().strip().splitlines()[-3:]))
        return None


def finding(severity, func, text, evidence=None):
    ''' severity: bug | risk | info | known '''
    FINDINGS.append({'severity': severity, 'func': func, 'text': text, 'evidence': evidence, 'part': PART[0]})
    print('[FINDING/%s] %s: %s%s' % (severity, func, text, '' if evidence is None else '  {%s}' % (evidence,)), flush=True)


def section(title):
    print('\n== %s' % title, flush=True)


def summary(script, extra=None):
    ''' print the totals, write results/<script>.json, return the exit code (1 if any check failed) '''
    n_fail = sum(1 for c in CHECKS if not c['ok'])
    os.makedirs(RESULTS, exist_ok=True)
    parts = {}
    for c in CHECKS:
        p = parts.setdefault(c['part'], {'checks': 0, 'failed': 0, 'findings': 0}); p['checks'] += 1; p['failed'] += not c['ok']
    for f in FINDINGS: parts.setdefault(f['part'], {'checks': 0, 'failed': 0, 'findings': 0})['findings'] += 1
    data = {'script': script, 'checks': len(CHECKS), 'failed': n_fail, 'seconds': round(time.time() - _T0, 1), 'parts': parts, 'results': CHECKS, 'findings': FINDINGS}
    if extra: data.update(extra)
    with open(os.path.join(RESULTS, script + '.json'), 'w') as fh: json.dump(data, fh, indent=1, default=str)
    print('\nSUMMARY %s: %d checks, %d failed, %d findings, %.1f s' % (script, len(CHECKS), n_fail, len(FINDINGS), time.time() - _T0), flush=True)
    return 1 if n_fail else 0


# ------------------------------------------------------------------ one fresh process per scenario
def run_worker(script, *args, timeout=3600):
    ''' run `python script --worker args...` in a fresh process; the worker prints one line 'JSON <payload>' '''
    r = subprocess.run([sys.executable, script, '--worker'] + [str(a) for a in args], capture_output=True, text=True, env=ENV, cwd=HERE, timeout=timeout)
    for line in r.stdout.splitlines()[::-1]:
        if line.startswith('JSON '): return json.loads(line[5:])
    raise RuntimeError('worker %s %s failed (rc=%d): %s' % (os.path.basename(script), args, r.returncode, (r.stderr or r.stdout)[-2000:]))


def worker_reply(payload):
    print('JSON ' + json.dumps(payload, default=str), flush=True)


def rss_mb():
    return hc.rss_mb()


def env_info():
    import pynemo  # noqa: F401
    from importlib.metadata import version
    return {'date': time.strftime('%Y-%m-%d %H:%M'), 'head': HEAD, 'head_subject': git('log', '-1', '--format=%s').strip(), 'branch': git('rev-parse', '--abbrev-ref', 'HEAD').strip(),
            'changed': {f: md5(f) for f in CHANGED}, 'python': platform.python_version(), 'pynemo': version('pynemo'), 'numpy': np.__version__,
            'machine': '%s, %s, %s CPUs' % (platform.platform(), subprocess.run(['sysctl', '-n', 'machdep.cpu.brand_string'], capture_output=True, text=True).stdout.strip() or platform.processor(), os.cpu_count())}


C = sys.modules[__name__]                                       # the helpers above, under the C. prefix the tests below use


# ###################################################################### unit
# Unit-level tests of the make_copy / make_clone migration (one process; part rng starts its own worker processes). Parts:
#   static     S1-S25 on the source only: no NotImplementedError / TODO(makeClone|copy_from|mixed) / .mutation( / genotype_set / deepcopy(individual) left;
#              signatures of the 6 changed files == HEAD (individual + make_copy, indfactory + make_clone); every is_remove default True (D1); the migration list
#              carries the NemoBackend tag and reaches Nemo in the call graph; make_copy / make_clone bodies, the 10 copy sites, the clone / offspring sites, D4
#   copy       individual.make_copy() = Nemo Individual::operator=: new ID, every other field + every allele copied, original untouched, age; phenotype identical
#              without Ve / with Ve = 0, a fresh environmental draw with Ve > 0 (variance == Ve)
#   clone      indfactory.make_clone() = IndFactory::makeClonedOffsprg: mu = 0 -> identical alleles; father 0, pedigree 4, age 0, home, gender kept, mother's
#              class-4 counters +1; number of changed alleles per clone == the law of Binomial(2L, mu) hits (mean, variance, distribution)
#   guards     bad inputs fail before Nemo is called and change nothing (factory counters, static ID, the mother)
#   lifecycle  generations of copy / clone / sexual birth + recycle: IDs unique, recycling pool reused, minted - recycled == alive, RSS flat, minted == recycled
#   rng        fresh processes: native twice -> identical; native == Python imitation without operator= / makeClonedOffsprg, record by record
#   bench      time per call for 20 / 80 / 320 loci (table only, no checks)

EXTRA = {}


# ====================================================================== static
LAYERS = ['metaibm/habitat.py', 'metaibm/patch.py', 'metaibm/metacommunity.py', 'extension/global_habitat_network.py']
SCANNED = sorted([os.path.join('metaibm', f) for f in os.listdir(os.path.join(C.REPO, 'metaibm')) if f.endswith('.py')] +
                 [os.path.join('extension', f) for f in os.listdir(os.path.join(C.REPO, 'extension')) if f.endswith('.py')])
NEMO_SEEDS = {'make_copy', 'make_clone', 'make_offspring', 'make_founder', 'recycle'}

COPY_SITES = {   # A: the 10 is_remove=False branches -> number of make_copy() calls in each
    'metaibm/habitat.py': {'_germinate_from_pools': 1},
    'metaibm/patch.py': {'patch_dipersal_within_from_offspring_and_dormancy_pool': 1, 'patch_dispersal_within_from_offspring_pool_to_immigrant_pool': 1,
                         'patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool': 1},
    'metaibm/metacommunity.py': {'meta_colonize_from_propagules_rains': 1, 'pairwise_sexual_colonization_from_prpagules_rains': 2, 'dispersal_among_patches_from_offspring_pool_and_dormancy_pool': 1,
                                 'dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool': 1, 'dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool': 1},
    'extension/global_habitat_network.py': {'dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool': 1}}
BIRTH_SITES = {'hab_asex_reproduce_mutate_with_num': 'make_clone', 'hab_mix_asex_reproduce_mutate_with_num': 'make_clone',
               'hab_sex_reproduce_mutate_with_num': 'make_offspring', 'hab_mix_sex_reproduce_mutate_with_num': 'make_offspring'}
B_FUNCS = ['hab_asex_reproduce_mutate_with_num', 'hab_mix_asex_reproduce_mutate_with_num', 'hab_mix_sex_reproduce_mutate_with_num', 'hab_asexual_reprodece_germinate',
           'hab_mixed_reproduce_germinate', 'hab_asex_reproduce_mutate_into_offspring_pool', 'hab_mix_reproduce_mutate_into_offspring_pool']
MIGRATION = {   # A + B + C + D (outer functions that got the tag)
    'metaibm/habitat.py': ['_germinate_from_pools', 'hab_local_germinate_from_offspring_and_dormancy_pool', 'hab_local_germinate_from_offspring_and_immigrant_pool',
                           'hab_local_germinate_from_offspring_immigrant_and_dormancy_pool'] + B_FUNCS,
    'metaibm/patch.py': list(COPY_SITES['metaibm/patch.py']) + ['patch_local_germinate_from_offspring_and_dormancy_pool', 'patch_local_germinate_from_offspring_and_immigrant_pool',
                         'patch_local_germinate_from_offspring_immigrant_and_dormancy_pool', 'patch_asexual_birth_germinate', 'patch_mixed_birth_germinate',
                         'patch_asex_reproduce_mutate_into_offspring_pool', 'patch_mix_reproduce_mutate_into_offspring_pool'],
    'metaibm/metacommunity.py': list(COPY_SITES['metaibm/metacommunity.py']) + ['meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool',
                         'meta_mainland_asexual_birth_mutate_germinate', 'meta_mainland_mixed_birth_mutate_germinate', 'meta_asex_reproduce_mutate_into_offspring_pool',
                         'meta_mix_reproduce_mutate_into_offspring_pool', 'meta_dispersal_within_patch_from_offspring_to_immigrant_pool',
                         'meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool', 'meta_local_germinate_from_offspring_and_dormancy_pool',
                         'meta_local_germinate_from_offspring_and_immigrant_pool', 'meta_local_germinate_from_offspring_immigrant_and_dormancy_pool'],
    'extension/global_habitat_network.py': list(COPY_SITES['extension/global_habitat_network.py'])}


def read(rel, rev=None):
    return C.git('show', '%s:%s' % (rev, rel)) if rev else open(os.path.join(C.REPO, rel)).read()


def functions(src):
    ''' {qualified name: FunctionDef} for module-level functions and class methods '''
    out = {}
    for n in ast.parse(src).body:
        if isinstance(n, ast.FunctionDef): out[n.name] = n
        elif isinstance(n, ast.ClassDef):
            for m in n.body:
                if isinstance(m, ast.FunctionDef): out['%s.%s' % (n.name, m.name)] = m
    return out


def short(q): return q.split('.')[-1]


def signature(fn):
    return ast.unparse(fn.args)


def calls(fn):
    out = []
    for n in ast.walk(fn):
        if isinstance(n, ast.Call):
            if isinstance(n.func, ast.Attribute): out.append(n.func.attr)
            elif isinstance(n.func, ast.Name): out.append(n.func.id)
    return out


def sec_static():
    C.part('static')
    SRC = {f: read(f) for f in SCANNED}

    # ============================================================ S1-S5 leftovers of the reserved / pre-Nemo code
    C.section('S1-S5 leftovers in metaibm/ + extension/')
    hits = [(f, t.start[0]) for f in SCANNED for t in tokenize.generate_tokens(io.StringIO(SRC[f]).readline) if t.type == tokenize.NAME and t.string == 'NotImplementedError']
    C.ok('S1 no NotImplementedError (code) in metaibm/ + extension/', not hits, hits or '%d files scanned' % len(SCANNED))
    hits = [(f, i + 1, m.group(0)) for f in SCANNED for i, line in enumerate(SRC[f].splitlines()) for m in re.finditer(r'TODO\((makeClone|copy_from|mixed)\)', line)]
    other_todo = sorted({m.group(0) for f in SCANNED for m in re.finditer(r'TODO\(\w+\)', SRC[f])})
    C.ok('S2 no TODO(makeClone|copy_from|mixed) in metaibm/ + extension/', not hits, hits or 'other TODO tags left: %s' % other_todo)
    hits = [(f, n.lineno) for f in SCANNED for n in ast.walk(ast.parse(SRC[f])) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == 'mutation']
    C.ok('S3 no .mutation( call in metaibm/ + extension/', not hits, hits or None)
    hits = [(f, n.lineno) for f in SCANNED for n in ast.walk(ast.parse(SRC[f])) if (isinstance(n, ast.Attribute) and n.attr == 'genotype_set') or (isinstance(n, ast.Name) and n.id == 'genotype_set')]
    C.ok('S4 no genotype_set in the code of metaibm/ + extension/', not hits, hits or None)
    ALLOWED_DEEPCOPY = {'self.empty_site_pos_ls', 'self.species_category'}      # position lists / dicts, no individuals
    dc = [(f, n.lineno, ast.unparse(n.args[0]) if n.args else '') for f in SCANNED for n in ast.walk(ast.parse(SRC[f]))
          if isinstance(n, ast.Call) and ((isinstance(n.func, ast.Attribute) and n.func.attr == 'deepcopy') or (isinstance(n.func, ast.Name) and n.func.id == 'deepcopy'))]
    bad = [d for d in dc if d[2] not in ALLOWED_DEEPCOPY]
    C.ok('S5 no deepcopy of an individual in metaibm/ + extension/', not bad, bad or 'remaining deepcopy args: %s' % sorted({d[2] for d in dc}))

    # ============================================================ S6-S11 signatures equal to HEAD
    C.section('S6-S11 signatures vs HEAD')
    ADDED = {'metaibm/individual.py': {'individual.make_copy'}, 'metaibm/indfactory.py': {'indfactory.make_clone'}}
    for k, f in enumerate(sorted(C.CHANGED)):
        now, head = functions(read(f)), functions(read(f, 'HEAD'))
        added, removed = sorted(set(now) - set(head)), sorted(set(head) - set(now))
        changed = sorted(q for q in set(now) & set(head) if signature(now[q]) != signature(head[q]))
        good = set(added) == ADDED.get(f, set()) and not removed and not changed
        C.ok('S%d signatures of %s == HEAD (+ %s)' % (6 + k, f, sorted(ADDED.get(f, [])) or 'nothing'), good,
             {'added': added, 'removed': removed, 'changed': [(q, signature(head[q]), signature(now[q])) for q in changed]} if not good else '%d functions' % len(now))
    if len(C.CHANGED) != 6: C.ok('changed files == the 6 of the migration', False, C.CHANGED)

    # ============================================================ S12 is_remove defaults
    C.section('S12 is_remove defaults (D1)')
    def param_defaults(fn):
        a = fn.args; pos = a.posonlyargs + a.args
        d = dict(zip([x.arg for x in pos][len(pos) - len(a.defaults):], a.defaults))
        d.update({x.arg: v for x, v in zip(a.kwonlyargs, a.kw_defaults)})
        return {x.arg: (ast.unparse(d[x.arg]) if d.get(x.arg) is not None else None) for x in pos + a.kwonlyargs}


    defaults = [(f, q, param_defaults(fn)['is_remove']) for f in SCANNED for q, fn in functions(SRC[f]).items() if 'is_remove' in param_defaults(fn)]
    not_true = [d for d in defaults if d[2] != 'True']
    C.ok('S12 every is_remove parameter defaults to True', not not_true and defaults, not_true or '%d functions with is_remove, all default True' % len(defaults))

    # ============================================================ S13-S20 NemoBackend tag + call graph
    C.section('S13-S20 migration list: NemoBackend tag + reaches Nemo')
    ALL_FN = {}
    for f in LAYERS + ['metaibm/individual.py', 'metaibm/indfactory.py']:
        for q, fn in functions(SRC[f]).items(): ALL_FN[(f, q)] = fn
    reach = set(NEMO_SEEDS)
    while True:
        new = {short(q) for (f, q), fn in ALL_FN.items() if short(q) not in reach and set(calls(fn)) & reach}
        if not new: break
        reach |= new
    for k, f in enumerate(LAYERS):
        fns = {short(q): fn for q, fn in functions(SRC[f]).items()}
        missing = [n for n in MIGRATION[f] if n not in fns]
        untagged = [n for n in MIGRATION[f] if n in fns and 'NemoBackend' not in (ast.get_docstring(fns[n]) or '')]
        C.ok('S%d %s: %d migration-list functions tagged NemoBackend' % (13 + k, f, len(MIGRATION[f])), not missing and not untagged, {'missing': missing, 'untagged': untagged} if missing or untagged else None)
    for k, f in enumerate(LAYERS):
        fns = {short(q): fn for q, fn in functions(SRC[f]).items()}
        cold = [n for n in MIGRATION[f] if n not in reach]
        C.ok('S%d %s: every migration-list function reaches Nemo in the call graph' % (17 + k, f), not cold, cold or None)
    # tag audit over the whole layer files (info only)
    tag_cold, hot_untagged = [], []
    for f in LAYERS:
        for q, fn in functions(SRC[f]).items():
            tagged = 'NemoBackend' in (ast.get_docstring(fn) or '')
            if tagged and short(q) not in reach: tag_cold.append((f, short(q)))
            if not tagged and short(q) in reach and not short(q).startswith('_'): hot_untagged.append((f, short(q)))
    if tag_cold: C.finding('info', 'tag audit', 'functions tagged NemoBackend that never reach Nemo in the call graph (tag describes the layer, not a Nemo call)', tag_cold)
    if hot_untagged: C.finding('info', 'tag audit', 'functions that reach Nemo (via make_copy / make_clone / make_offspring / make_founder / recycle) but carry no NemoBackend tag', hot_untagged)

    # ============================================================ S21-S25 bodies
    C.section('S21-S25 bodies')
    fi = functions(SRC['metaibm/individual.py'])['individual.make_copy']
    cl = calls(fi)
    C.ok('S21 individual.make_copy: _require_alive -> getNewIndividual().make_copy(nemo_indi_obj) -> _wrap_nemo_individual, age copied',
         cl[:1] == ['_require_alive'] and 'getNewIndividual' in cl and 'make_copy' in cl and '_wrap_nemo_individual' in cl and 'copy_indi_obj.age = self.age' in ast.unparse(fi), cl)
    fc = functions(SRC['metaibm/indfactory.py'])['indfactory.make_clone']
    cl = calls(fc)
    C.ok('S22 indfactory.make_clone: _require_own(mother) before makeClonedOffsprg(mother.nemo_indi_obj, home_habitat), result wrapped',
         cl and cl[0] == '_require_own' and 'makeClonedOffsprg' in cl and '_wrap_nemo_individual' in cl and signature(fc) == 'self, mother_indi_obj, home_habitat=0', (cl, signature(fc)))
    sites, want = {}, {}
    for f in LAYERS:
        for q, fn in functions(SRC[f]).items():
            n = calls(fn).count('make_copy')
            if n: sites[(f, short(q))] = n
        for name, n in COPY_SITES[f].items(): want[(f, name)] = n
    C.ok('S23 make_copy() is called exactly in the 10 is_remove=False branches (11 call sites)', sites == want,
         {'unexpected': {k: v for k, v in sites.items() if want.get(k) != v}, 'missing': {k: v for k, v in want.items() if sites.get(k) != v}} if sites != want else '%d functions, %d call sites' % (len(sites), sum(sites.values())))
    hab = {short(q): fn for q, fn in functions(SRC['metaibm/habitat.py']).items()}
    bad = []
    for name, method in BIRTH_SITES.items():
        body = ast.unparse(hab[name])
        if method not in calls(hab[name]) or 'home_habitat=self.global_index' not in body: bad.append(name)
    clone_sites = sorted(short(q) for (f, q), fn in ALL_FN.items() if f in LAYERS and 'make_clone' in calls(fn))
    C.ok('S24 asexual births -> make_clone, sexual births -> make_offspring, all with home_habitat=self.global_index; make_clone called only in the two asexual with_num',
         not bad and clone_sites == ['hab_asex_reproduce_mutate_with_num', 'hab_mix_asex_reproduce_mutate_with_num'], {'bad': bad, 'make_clone sites': clone_sites})
    used = []
    for name in B_FUNCS:
        fn = hab[name]
        params = [a.arg for a in fn.args.args]
        if 'mutation_rate' not in params or 'pheno_var_ls' not in params: used.append((name, 'parameter missing'))
        reads = [n for n in ast.walk(fn) if isinstance(n, ast.Name) and n.id in ('mutation_rate', 'pheno_var_ls') and isinstance(n.ctx, ast.Load)]
        passes = [n for n in ast.walk(fn) if isinstance(n, ast.Call) for a in n.args if isinstance(a, ast.Name) and a.id in ('mutation_rate', 'pheno_var_ls')]
        if len(reads) != len(passes): used.append((name, 'read outside a pass-through call'))
    C.ok('S25 D4: mutation_rate / pheno_var_ls stay in the B signatures and are only passed through (never used)', not used, used or '%d functions' % len(B_FUNCS))


# ====================================================================== copy
def sec_copy():
    C.part('copy')
    random.seed(3)
    F = C.make_species('sp_copy', 0.2, seed=5)                     # Ve = 0.002025, mu = 1e-4
    G = C.make_species('sp_other', 0.5, seed=6)

    # ============================================================ 1. fields, alleles, ID, counters
    C.section('1. fields / alleles / ID / counters')
    mother, father = F.make_founder('female', 3), F.make_founder('male', 4)
    x = F.make_offspring(mother, father, 'female', home_habitat=7)        # real parent IDs, pedigree class, home 7
    for _ in range(3): F.make_offspring(x, father, 'male').recycle()        # x's mating / realized-fecundity counters move
    x.age = 5; x.nemo_indi_obj.setAge(5)
    x.nemo_indi_obj.setFecundity(2.5)
    before, ph, xid = C.copied_state(x), x.get_indi_phenotype_ls(), x.individual_id
    next_id, m0, r0 = indfactory.get_current_individual_id(), F.minted_num, F.recycled_num
    c = x.make_copy()
    C.ok('copy: every operator= field + alleles == the original (sex, age, mother/father ID, home, pedigree, fecundity, 5+5 counters, %d alleles)' % (2 * sum(len(t) for t in C.alleles(x))),
         C.copied_state(c) == before, {'counters': before[9], 'pedigree': before[7], 'home': before[6]})
    C.ok('copy: new ID = the static counter before the call, counter +1', c.individual_id == next_id != xid and indfactory.get_current_individual_id() == next_id + 1, (xid, c.individual_id, next_id))
    C.ok('copy: minted +1, recycled unchanged', (F.minted_num - m0, F.recycled_num - r0) == (1, 0), (F.minted_num - m0, F.recycled_num - r0))
    C.ok('copy: original unchanged (state, phenotype, ID, alive)', C.copied_state(x) == before and x.get_indi_phenotype_ls() == ph and x.individual_id == xid and x.is_alive)
    C.ok('copy: same species / factory / pheno names, MetaIBM age and Nemo age copied (5)',
         c.indfactory is F and c.species_id == x.species_id and c.pheno_names_ls == x.pheno_names_ls and c.age == 5 and c.nemo_indi_obj.getAge() == 5)
    C.ok('copy: its own C++ object', c is not x and int(c.nemo_indi_obj.this) != int(x.nemo_indi_obj.this))
    gx = C.genotypic_values(x)
    C.ok('copy: genotypic values equal (phenotype - G is the Ve draw, |e| < 6 sd)', C.genotypic_values(c) == gx and all(abs(p - g) < 6 * math.sqrt(0.002025) for p, g in zip(c.get_indi_phenotype_ls(), gx)))
    x.age = 9
    C.ok('copy: MetaIBM age is a snapshot (later changes to the original do not propagate)', c.age == 5)
    c2 = c.make_copy()
    C.ok('copy of a copy: same copied state, third ID', C.copied_state(c2)[1:] == C.copied_state(c)[1:] and len({x.individual_id, c.individual_id, c2.individual_id}) == 3)
    c.recycle()
    C.ok('recycling the copy leaves the original alive and unchanged', x.is_alive and C.copied_state(x)[3:] == before[3:] and x.get_indi_phenotype_ls() == ph)
    kid = F.make_offspring(c2, father, 'male', home_habitat=9)
    C.ok('a copy breeds like any individual (mother_id = the copy\'s ID)', kid.mother_id == c2.individual_id and kid.father_id == father.individual_id)
    male = F.make_founder('male', 11); cm = male.make_copy()
    C.ok('copy keeps gender male, founder sentinel parent IDs and home', cm.gender == 'male' and (cm.mother_id, cm.father_id, cm.nemo_indi_obj.getHome()) == (male.mother_id, male.father_id, 11) and cm.mother_id >= 2 ** 63)
    C.finding('info', 'individual.make_copy (all 10 is_remove=False sites)', 'a copy keeps everything of its source but the ID: home (birth habitat), mother / father IDs, pedigree class, '
              'mating counters, age -- a dispersed or colonizing copy is NOT re-homed to its destination (as the v3 alias was the same object); Nemo computes the pedigree class '
              'of its offspring from these homes / parent IDs. MetaIBM never reads the copied _mother / _father raw pointers (IDs only).',
              {'source home / copy-of-copy home': (7, c2.nemo_indi_obj.getHome()), 'source mother_id / copy-of-copy mother_id': (C.copied_state(x)[4], C.copied_state(c2)[4]), 'counters copied': before[9],
               'pedigree class of copy (home 7) x father (home 4) offspring': kid.nemo_indi_obj.getPedigreeClass()})
    for y in (kid, c2, cm, male, x, mother, father): y.recycle()

    # ============================================================ 2. phenotype: no Ve / Ve = 0 identical, Ve > 0 redrawn
    C.section('2. phenotype and Ve')
    for label, ve in (('no Ve parameter', None), ('Ve = 0', 0.0)):
        S = C.make_species('sp_ve_%s' % ('none' if ve is None else '0'), 0.3, seed=8, env_variance=ve)
        xs = [S.make_founder('female') for _ in range(20)]
        same = all(y.make_copy().get_indi_phenotype_ls() == y.get_indi_phenotype_ls() for y in xs for _ in range(5))
        C.ok('%s: copy phenotype == original phenotype exactly (20 individuals x 5 copies)' % label, same)
    VE, N = 0.01, 3000
    S = C.make_species('sp_ve', 0.3, seed=9, env_variance=VE)
    y = S.make_founder('female'); gy = C.genotypic_values(y); py = y.get_indi_phenotype_ls()
    e = [[], []]
    for _ in range(N):
        cy = y.make_copy()
        for d, (p, g) in enumerate(zip(cy.get_indi_phenotype_ls(), gy)): e[d].append(p - g)
        cy.recycle()
    for d in range(2):
        mean = sum(e[d]) / N; var = sum((v - mean) ** 2 for v in e[d]) / (N - 1)
        C.ok('Ve = %g, axis %d: copy phenotype - G ~ N(0, Ve) over %d copies (mean within 4 SE, variance within 4 SE)' % (VE, d, N),
             abs(mean) < 4 * math.sqrt(VE / N) and abs(var - VE) < 4 * VE * math.sqrt(2 / (N - 1)), {'mean': round(mean, 5), 'var': round(var, 5), 'Ve': VE})
    C.ok('Ve > 0: the original\'s phenotype is never redrawn by copying', y.get_indi_phenotype_ls() == py)
    C.ok('balance after the Ve runs: minted - recycled == 1 (the original)', S.minted_num - S.recycled_num == 1, (S.minted_num, S.recycled_num))


# ====================================================================== clone
def sec_clone():
    C.part('clone')
    random.seed(4)

    # ============================================================ 1. mu = 0: an exact genetic copy with clone pedigree
    C.section('1. mu = 0')
    F = C.make_species('sp_clone0', 0.2, seed=11, mutation_rate=0.0)
    m_f, m_m = F.make_founder('female', 2), F.make_founder('male', 3)
    for mother in (m_f, m_m):
        mother.age = 4; mother.nemo_indi_obj.setAge(4)
        st = C.copied_state(mother); ph = mother.get_indi_phenotype_ls(); mid = mother.individual_id
        next_id, m0 = indfactory.get_current_individual_id(), F.minted_num
        k = F.make_clone(mother, home_habitat=17)
        n = k.nemo_indi_obj
        C.ok('%s mother: alleles identical (mu = 0)' % mother.gender, C.alleles(k) == C.alleles(mother))
        C.ok('%s mother: mother_id = mother, father_id 0, pedigree class 4, age 0 (MetaIBM + Nemo), home 17, gender kept' % mother.gender,
             (k.mother_id, k.father_id, n.getPedigreeClass(), k.age, n.getAge(), n.getHome(), k.gender) == (mid, 0, 4, 0, 0, 17, mother.gender),
             (k.mother_id, k.father_id, n.getPedigreeClass(), k.age, n.getAge(), n.getHome(), k.gender))
        C.ok('%s mother: clone counters reset to 0, fecundity copied' % mother.gender,
             all(n.getMatings(c) == 0 and n.getRealizedFecundity(c) == 0 for c in range(5)) and n.getFecundity() == mother.nemo_indi_obj.getFecundity())
        C.ok('%s mother: new ID = static counter, minted +1' % mother.gender, k.individual_id == next_id and F.minted_num == m0 + 1, (k.individual_id, next_id))
        nm = mother.nemo_indi_obj
        C.ok('%s mother: mother\'s class-4 mating + realized-fecundity counters +1, nothing else of the mother changed' % mother.gender,
             nm.getMatings(4) == st[9][4] + 1 and nm.getRealizedFecundity(4) == st[10][4] + 1 and C.copied_state(mother)[:9] == st[:9]
             and C.copied_state(mother)[11] == st[11] and mother.get_indi_phenotype_ls() == ph and mother.individual_id == mid)
        C.ok('%s mother: clone genotypic value == mother\'s, phenotype recomputed (|pheno - G| < 6 sd of Ve)' % mother.gender,
             C.genotypic_values(k) == C.genotypic_values(mother) and all(abs(p - g) < 6 * math.sqrt(0.002025) for p, g in zip(k.get_indi_phenotype_ls(), C.genotypic_values(k))))
    kk = F.make_clone(k, home_habitat=18)
    C.ok('clone of a clone: mother_id = the first clone, alleles still the founder\'s (mu = 0)', kk.mother_id == k.individual_id and C.alleles(kk) == C.alleles(m_m))
    kid = F.make_offspring(F.make_clone(m_f, 5), m_m, 'female', home_habitat=5)
    C.ok('a clone breeds sexually; its sexual offspring has pedigree class != 4', kid.nemo_indi_obj.getPedigreeClass() != 4 and kid.father_id == m_m.individual_id, kid.nemo_indi_obj.getPedigreeClass())

    # ============================================================ 2. mutation count ~ Binomial(2L, mu) (Nemo TProtoQuanti::get_num_mutations)
    C.section('2. mutation count')
    MU, L, N = 0.01, 40, 4000
    M = C.make_species('sp_clone_mu', 0.2, seed=12, mutation_rate=MU, loci=L)
    mom = M.make_founder('female', 1); ref = [a for t in C.alleles(mom) for pair in t for a in pair]
    counts = []
    for _ in range(N):
        k = M.make_clone(mom, 1)
        counts.append(sum(1 for a, b in zip(ref, [a for t in C.alleles(k) for pair in t for a in pair]) if a != b))
        k.recycle()
    TWO_L = 2 * L
    # exact law of the number of DISTINCT mutated allele slots: N_mut ~ Binomial(2L, mu), each hit lands on a uniform slot (locus, strand) of 2L
    pn = [math.comb(TWO_L, n) * MU ** n * (1 - MU) ** (TWO_L - n) for n in range(TWO_L + 1)]
    pd = [0.0] * (TWO_L + 1)
    occ = [1.0] + [0.0] * TWO_L                                     # P(d distinct | n hits), updated hit by hit
    for n in range(TWO_L + 1):
        for d in range(TWO_L + 1): pd[d] += pn[n] * occ[d]
        occ = [occ[d] * d / TWO_L + (occ[d - 1] * (TWO_L - d + 1) / TWO_L if d else 0.0) for d in range(TWO_L + 1)]
    e_mean = sum(d * p for d, p in enumerate(pd)); e_var = sum((d - e_mean) ** 2 * p for d, p in enumerate(pd)); e_m4 = sum((d - e_mean) ** 4 * p for d, p in enumerate(pd))
    s_mean = sum(counts) / N; s_var = sum((c - s_mean) ** 2 for c in counts) / (N - 1)
    C.ok('changed alleles per clone: mean == E[distinct slots of Binomial(2L=%d, mu=%g) hits] within 4 SE (%d clones)' % (TWO_L, MU, N),
         abs(s_mean - e_mean) < 4 * math.sqrt(e_var / N), {'sample': round(s_mean, 4), 'expected': round(e_mean, 4), 'binomial 2L*mu': TWO_L * MU})
    C.ok('changed alleles per clone: variance within 4 SE', abs(s_var - e_var) < 4 * math.sqrt((e_m4 - e_var ** 2) / N), {'sample': round(s_var, 4), 'expected': round(e_var, 4), 'binomial': round(TWO_L * MU * (1 - MU), 4)})
    obs = [counts.count(d) for d in range(3)] + [sum(1 for c in counts if c >= 3)]
    exp = [N * pd[d] for d in range(3)] + [N * sum(pd[3:])]
    okc, stat = C.hc.chi2_ok(obs, exp)
    C.ok('changed alleles per clone: distribution {0, 1, 2, >=3} fits (chi2, loose bound)', okc, {'observed': obs, 'expected': [round(e, 1) for e in exp], 'chi2': stat})
    C.ok('the mother is never mutated by cloning (%d clones)' % N, [a for t in C.alleles(mom) for pair in t for a in pair] == ref and mom.nemo_indi_obj.getMatings(4) == N)
    C.ok('balance: minted - recycled == 1 (the mother)', M.minted_num - M.recycled_num == 1, (M.minted_num, M.recycled_num))


# ====================================================================== guards
def sec_guards():
    C.part('guards')
    F = C.make_species('sp_g1', 0.2, seed=21)
    G = C.make_species('sp_g2', 0.5, seed=22)
    mom = F.make_founder('female', 1)
    dead = F.make_founder('male', 1); dead.recycle()
    other = G.make_founder('female', 1)


    def snap():
        return (F.minted_num, F.recycled_num, G.minted_num, G.recycled_num, indfactory.get_current_individual_id(), C.copied_state(mom), mom.get_indi_phenotype_ls(),
                C.copied_state(other), other.get_indi_phenotype_ls())


    CASES = [('make_clone(recycled mother)', RuntimeError, lambda: F.make_clone(dead, 1)),
             ('make_clone(mother of another species)', RuntimeError, lambda: F.make_clone(other, 1)),
             ('make_clone(None)', AttributeError, lambda: F.make_clone(None, 1)),
             ('make_clone(fake individual)', AttributeError, lambda: F.make_clone(C.hc.fake_indi('x'), 1)),
             ('make_clone(home_habitat=65536)', (OverflowError, TypeError), lambda: F.make_clone(mom, 65536)),
             ('make_clone(home_habitat=-1)', (OverflowError, TypeError), lambda: F.make_clone(mom, -1)),
             ('make_clone(home_habitat="3")', TypeError, lambda: F.make_clone(mom, '3')),
             ('make_clone(home_habitat=2.5)', TypeError, lambda: F.make_clone(mom, 2.5)),
             ('recycled.make_copy()', RuntimeError, lambda: dead.make_copy()),
             ('make_copy() of an individual whose nemo_indi_obj is gone', (RuntimeError, AttributeError), lambda: dead.make_copy())]
    for label, exc, fn in CASES:
        before = snap(); err = None
        try: fn()
        except Exception as e: err = e
        C.ok('%s -> %s, nothing changed' % (label, exc.__name__ if isinstance(exc, type) else '/'.join(x.__name__ for x in exc)),
             isinstance(err, exc) and snap() == before, '%s: %s' % (type(err).__name__, err) if err is not None else 'no error')

    nid = indfactory.get_current_individual_id()
    k = F.make_clone(mom, 65535)
    C.ok('after the refused calls: next clone gets the expected ID, home 65535 accepted (upper bound)', k.individual_id == nid and k.nemo_indi_obj.getHome() == 65535, (k.individual_id, nid))
    c = k.make_copy(); c.recycle()
    try: c.recycle(); C.ok('copy recycled twice -> RuntimeError', False)
    except RuntimeError: C.ok('copy recycled twice -> RuntimeError', True)
    C.ok('balance: minted - recycled == alive (mom, other, clone)', (F.minted_num - F.recycled_num, G.minted_num - G.recycled_num) == (2, 1), (F.minted_num, F.recycled_num, G.minted_num, G.recycled_num))


# ====================================================================== lifecycle
def sec_lifecycle():
    C.part('lifecycle')
    random.seed(31)
    F = C.make_species('sp_life', 0.2, seed=31, mutation_rate=0.001)
    N, GENS = 300, 60
    pop = [F.make_founder('female' if i % 2 else 'male', 1) for i in range(N)]
    seen = {x.individual_id for x in pop}
    freed, reused, born, dup_id, unbalanced = set(), 0, 0, 0, []
    rss, t0 = {}, time.time()
    for gen in range(1, GENS + 1):
        fem = [x for x in pop if x.gender == 'female']; mal = [x for x in pop if x.gender == 'male']
        new = []
        for i in range(N):
            r = i % 3
            if r == 0: y = random.choice(pop).make_copy()
            elif r == 1: y = F.make_clone(random.choice(pop), gen % 7)
            else: y = F.make_offspring(random.choice(fem), random.choice(mal), random.choice(('female', 'male')), gen % 7)
            reused += int(y.nemo_indi_obj.this) in freed; born += 1
            dup_id += y.individual_id in seen; seen.add(y.individual_id)
            new.append(y)
        for x in pop:
            freed.add(int(x.nemo_indi_obj.this)); x.recycle()
        pop = new
        if F.minted_num - F.recycled_num != N: unbalanced.append((gen, F.minted_num, F.recycled_num))
    el = time.time() - t0
    C.ok('%d generations x %d births (copy / clone / sexual): IDs unique' % (GENS, N), dup_id == 0, dup_id)
    C.ok('minted - recycled == %d alive after every generation' % N, not unbalanced, unbalanced[:3] or None)
    C.ok('births placed in recycled memory >= births - N (recycling pool reused by copies and clones too)', reused >= born - N, '%d / %d' % (reused, born))
    # phase 2: the same mix without the test's own ID / address sets (they grow by themselves), RSS sampled along the way; the allocator warms up for a few hundred
    # generations (RSS rises < 1 MB until ~gen 200-800, then stays flat to gen 2000 -- measured) and may still take a single ~0.4 MB step later, so the check is
    # growth per birth over a long window: gen 600 -> 3000 = 720k births, < 1 MB (a leak of > ~1.5 bytes per birth would show; operator= reallocates 640 bytes per copy)
    G2 = 3000
    for gen in range(1, G2 + 1):
        fem = [x for x in pop if x.gender == 'female']; mal = [x for x in pop if x.gender == 'male']
        new = [pop[i].make_copy() if i % 3 == 0 else F.make_clone(pop[i], 1) if i % 3 == 1 else F.make_offspring(random.choice(fem), random.choice(mal), 'female' if i % 2 else 'male', 1) for i in range(N)]
        for x in pop: x.recycle()
        pop = new
        if gen in (50, 200, 600, 1000, 1500, 2000, 2500, G2): rss[gen] = round(C.rss_mb(), 2)
    C.ok('RSS does not grow with births after the allocator warm-up (gen 600 -> %d, %d births: growth < 1 MB)' % (G2, (G2 - 600) * N), rss[G2] - rss[600] < 1.0,
         {'rss_mb': rss, 'bytes per birth': round((rss[G2] - rss[600]) * 2 ** 20 / ((G2 - 600) * N), 3)})
    C.ok('minted - recycled == %d after phase 2' % N, F.minted_num - F.recycled_num == N, (F.minted_num, F.recycled_num))
    for x in pop: x.recycle()
    C.ok('after the last recycle: minted == recycled', F.minted_num == F.recycled_num, (F.minted_num, F.recycled_num))
    F.purge_recycling_pool()
    print('%d births in %.2f s (%.1f us each incl. recycle)' % (born, el, 1e6 * el / born))
    EXTRA['lifecycle'] = {'births': born, 'seconds_births': round(el, 2), 'rss_mb': rss}


# ====================================================================== rng (fresh processes)
STEPS, POP = 300, 30


def imit_assign(F, x):
    ''' Individual::operator= rebuilt from setters: every field but the ID, every allele, then set_value() (the phenotype draw operator= makes) '''
    s = x.nemo_indi_obj
    n = F._nemo_indfactory.getNewIndividual()
    n.setSex(s.getSex()); n.setAge(s.getAge()); n.setMotherID(s.getMotherID()); n.setFatherID(s.getFatherID()); n.setHome(s.getHome())
    n.setPedigreeClass(s.getPedigreeClass()); n.setFecundity(s.getFecundity())
    for c in range(5):
        assert s.getMatings(c) == s.getRealizedFecundity(c)          # only DidHaveABaby moves them here
        for _ in range(s.getMatings(c)): n.DidHaveABaby(c)
    idx = F.traits_index['quant']; tn, ts = n.getTrait(idx), s.getTrait(idx)
    seq = sorted({ts.get_locus_seq_pos(l, d) for d in range(F.traits_num) for l in range(ts.get_num_loci(d))})
    for pos in seq:
        for a in (0, 1): tn.set_allele_value(pos, a, ts.get_allele_value(pos, a))
    tn.set_value()
    return n


def imit_copy(F, x):
    c = F._wrap_nemo_individual(imit_assign(F, x)); c.age = x.age
    return c


def imit_clone(F, m, home):
    n = imit_assign(F, m)
    n.reset_counters(); n.setFather(None); n.setFatherID(0); n.setMother(m.nemo_indi_obj); n.setMotherID(m.nemo_indi_obj.getID())
    n.setIsSelfed(True); n.setHome(home); n.setAge(0)
    m.nemo_indi_obj.DidHaveABaby(n.getPedigreeClass())
    n.create(False, True)
    return F._wrap_nemo_individual(n)


def record(op, y, F):
    n = y.nemo_indi_obj
    return [op, y.individual_id, y.gender, y.age, n.getAge(), str(n.getMotherID()), str(n.getFatherID()), n.getHome(), n.getPedigreeClass(), n.getFecundity(),
            [n.getMatings(c) for c in range(5)], [repr(v) for t in C.alleles(y) for pair in t for v in pair], [repr(v) for v in y.get_indi_phenotype_ls()],
            F.minted_num, F.recycled_num, indfactory.get_current_individual_id()]


def scenario(mode, seed):
    random.seed(seed)
    F = C.make_species('sp_rng', 0.2, seed=seed, mutation_rate=0.01, env_variance=0.01)
    pop = [F.make_founder(random.choice(('female', 'male')), random.randrange(100)) for _ in range(POP)]
    out = []
    for step in range(STEPS):
        r, x = random.random(), random.choice(pop)
        x.age = random.randrange(10); x.nemo_indi_obj.setAge(x.age)
        if r < 0.4: op, y = 'copy', (x.make_copy() if mode == 'native' else imit_copy(F, x))
        elif r < 0.8: op, y = 'clone', (F.make_clone(x, step % 50) if mode == 'native' else imit_clone(F, x, step % 50))
        else:
            fem = [z for z in pop if z.gender == 'female'] or [x]; mal = [z for z in pop if z.gender == 'male'] or [x]
            op, y = 'sex', F.make_offspring(random.choice(fem), random.choice(mal), random.choice(('female', 'male')), step % 50)
        out.append(record(op, y, F))
        k = random.randrange(POP); pop[k].recycle(); pop[k] = y
    return out


def sec_rng():
    C.part('rng')
    me = os.path.abspath(__file__)
    for seed in (1, 2):
        a, b, im = (C.run_worker(me, m, seed) for m in ('native', 'native', 'imitation'))
        ops = [r[0] for r in a['records']]
        C.ok('seed %d: native twice in fresh processes -> identical (%d records: %d copies, %d clones, %d sexual)' % (seed, STEPS, ops.count('copy'), ops.count('clone'), ops.count('sex')),
             a['digest'] == b['digest'], (a['digest'][:12], b['digest'][:12]))
        first = next((i for i, (r1, r2) in enumerate(zip(a['records'], im['records'])) if r1 != r2), None)
        C.ok('seed %d: native == Python imitation without operator= / makeClonedOffsprg, record by record (IDs, fields, %d alleles, phenotypes, counters)' % (seed, len(a['records'][0][11])),
             first is None and len(a['records']) == len(im['records']),
             None if first is None else {'first_diff': first, 'native': a['records'][first][:11], 'imitation': im['records'][first][:11]})
    c1 = C.run_worker(me, 'native', 3)
    C.ok('different seed -> different record', c1['digest'] != a['digest'])


# ====================================================================== bench
def sec_bench():
    C.part('bench')
    REPS = 20000
    random.seed(1)
    rows = []
    for loci in (20, 80, 320):
        F = C.make_species('sp_bench_%d' % loci, 0.2, seed=loci, loci=loci)
        pop = [F.make_founder('female' if i % 2 else 'male', 1) for i in range(100)]
        fem, mal = pop[1::2], pop[0::2]
        ops = {'make_copy': lambda: random.choice(pop).make_copy(),
               'make_clone': lambda: F.make_clone(random.choice(pop), 1),
               'make_offspring': lambda: F.make_offspring(random.choice(fem), random.choice(mal), 'female', 1),
               'make_founder': lambda: F.make_founder('female', 1),
               'choice + recycle only': lambda: random.choice(pop) and None}
        row = {'loci': loci}
        for name, fn in ops.items():
            for _ in range(200):                                            # warm the recycling pool
                y = fn()
                if y is not None: y.recycle()
            t = time.perf_counter()
            for _ in range(REPS):
                y = fn()
                if y is not None: y.recycle()
            row[name] = round(1e6 * (time.perf_counter() - t) / REPS, 2)
        rows.append(row)
        for x in pop: x.recycle()
        F.purge_recycling_pool()

    cols = ['make_copy', 'make_clone', 'make_offspring', 'make_founder', 'choice + recycle only']
    print('\nus per call (incl. recycle of the result), %d reps, Ve > 0, mu = 1e-4' % REPS)
    print('| loci | ' + ' | '.join(cols) + ' |'); print('|---' * (len(cols) + 1) + '|')
    for r in rows: print('| %d | ' % r['loci'] + ' | '.join('%.2f' % r[c] for c in cols) + ' |')
    EXTRA['bench_us_per_call'] = rows; EXTRA['bench_reps'] = REPS


UNIT_PARTS = {'static': sec_static, 'copy': sec_copy, 'clone': sec_clone, 'guards': sec_guards, 'lifecycle': sec_lifecycle, 'rng': sec_rng, 'bench': sec_bench}


# ###################################################################### core
# Core-layer tests with fake and real Nemo individuals, in two parts:
# copy (A): the 10 is_remove=False sites -- the original stays where it is, a make_copy() germinates / disperses / colonizes.
#   part a (fake individuals, no Nemo), seed by seed:
#     a1  8 functions (_germinate_from_pools via its 3 wrappers, 3 patch, 3 metacommunity x uniform / exponential rho=2, the extension x both kernels) against the v3
#         bodies compiled into a subclass (habitat / patch = baace34, metacommunity = 7a1adbb, extension = 1356387): two identical worlds from one seed, one call each
#         with is_remove=False (v3 extension: no is_remove), same seed -> layout (origin labels) == v3, random + np.random == v3, source pools untouched,
#         other pools only appended to, residents stay, everything new is a copy (placed + recycled == copies made), RECYCLED holds copies only
#     a2  the 2 colonizations (not against v3: v3 shuffles the mainland list first): same world, same seed, is_remove=False vs True -> same picks (by origin) on the
#         same sites, same RNG; False: mainland untouched, copies placed; True: the picked mainland sites emptied, the originals placed; integer and fractional rain
#   part b (real Nemo individuals): all 10 functions with is_remove=False on mc.make_mainland + mc.make_islands + mc.fill_pools worlds -> originals in place with
#     unchanged copied_state + phenotype, each new individual's copied_state == an original of the source pool(s), minted += copies, recycled += copies not placed,
#     mc.check_world after every call, mc.teardown (minted == recycled) after every world
# births (B + C) with real Nemo individuals: asexual / mixed births (make_clone, make_offspring) and marker births.
#   1  with_num: hab_asex_ / hab_mix_asex_ / hab_mix_sex_reproduce_mutate_with_num -- the mothers (pairs) are exactly random.sample(list, num) of the same RNG state
#      (asex: occupied_site_pos_ls, mix asex: asexual_parent_pos_ls, mix sex: pairs rebuilt from species_category_for_sexual_parents_pos as the code does, + gender draws);
#      clones: father_id 0, pedigree 4, home = global_index, gender + species = mother's; sexual: same-species female x male of the mixed-sexual parents, pedigree != 4
#   2  hab_asexual_reprodece_germinate: births == min(empty, int(indi_num x rate)) on every seed, every newborn a clone of a resident
#   3  hab_mixed_reproduce_germinate: asex_num = min(a, np.around(E a / (a + s))), sex_num = min(s, E - asex_num) -> births == min(E, a + s), every offspring made is
#      placed (nothing recycled), clones placed == asex_num; the grid reaches the x.5 ties where the v3 split gave E + 1 / E - 1;
#      the special case 3 empty sites, asexual expectation 3, sexual expectation 3 -> 2 clones + 1 sexual, 3 made, 3 born, 0 recycled
#   4  hab_asex_ / hab_mix_reproduce_mutate_into_offspring_pool: pool size floor or ceil of the expectation, mean within 4 SE over 400 seeds, last round's offspring recycled
#      (check_held); + the patch / metacommunity wrappers of the asex / mix births (D) on a whole metacommunity
#   5  marker births (C): 'asexual' / 'mix_asexual' / 'mix_sexual' (+ 'sexual' control) markers placed by hand -> mother from the birth habitat (mix_asexual: from its
#      asexual_parent_pos_ls), home = birth habitat's global_index, births == min(empty, markers), is_remove=True removes the born markers / False keeps all, no leak
#   6  probes (recorded as findings, not FAIL): D5 for every marker kind; a stale asexual_parent_pos_ls before a mixed birth

F = FS = None                                                   # the two species of the births part (births_species())


def births_species():
    global F, FS
    F = mc.species(41); FS = list(F.values())


# ====================================================================== copy (A): the 10 is_remove=False sites
warnings.filterwarnings('ignore', category=RuntimeWarning)       # 0/0 in the dispersal matrices of empty habitats (NaN -> 0 by the rounding, as designed)
SEEDS = range(200)
GERM = {'hab_local_germinate_from_offspring_and_dormancy_pool': ('offspring_pool', 'dormancy_pool'),
        'hab_local_germinate_from_offspring_and_immigrant_pool': ('offspring_pool', 'immigrant_pool'),
        'hab_local_germinate_from_offspring_immigrant_and_dormancy_pool': ('offspring_pool', 'immigrant_pool', 'dormancy_pool')}
PATCH = {'patch_dipersal_within_from_offspring_and_dormancy_pool': ('offspring_pool', 'dormancy_pool'),
         'patch_dispersal_within_from_offspring_pool_to_immigrant_pool': ('offspring_pool',),
         'patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool': ('offspring_pool', 'dormancy_pool')}
META = {'dispersal_among_patches_from_offspring_pool_and_dormancy_pool': ('offspring_pool', 'dormancy_pool'),
        'dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool': ('offspring_pool',),
        'dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool': ('offspring_pool', 'dormancy_pool')}
EXT = 'dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool'
COLON = ('meta_colonize_from_propagules_rains', 'pairwise_sexual_colonization_from_prpagules_rains')
KERNELS = (('uniform', {}), ('exponential', {'rho': 2}))
RES_KEYS = ['layout == v3', 'random + np.random == v3', 'source pools untouched (same objects, same order)', 'other pools only appended to', 'residents stay',
            'everything new is a copy (placed + recycled == copies made)', 'RECYCLED holds copies only']


def rate(s):
    return round(0.1 + 0.8 * ((s * 7919) % 101) / 100, 4)


# ============================================================ part a: worlds of fake individuals
def fill(h, rng, tag, occ=(0.0, 1.0), pools=True):
    sites = [(i, j) for i in range(h.length) for j in range(h.width)]
    pos = rng.sample(sites, rng.randint(int(len(sites) * occ[0]), int(len(sites) * occ[1])))
    n = len(sites)
    O, I, D = (rng.randint(0, 2 * n), rng.randint(0, n), rng.randint(0, n)) if pools else (0, 0, 0)
    pc.fill_fakes(h, tag, positions=pos, O=O, I=I, D=D)


def build_hab(s, cls):
    rng = random.Random(s); np.random.seed(s)
    h = hc.make_hab(L=rng.randint(1, 6), W=rng.randint(1, 6), K=10, cls=cls)
    fill(h, rng, 'h')
    return h


def build_patch(s, cls):
    rng = random.Random(s); np.random.seed(s)
    n = rng.randint(2, 4)
    p = pc.make_patch(n, K=10, cls=cls, sizes=[(rng.randint(1, 5), rng.randint(1, 5)) for _ in range(n)])
    for h in pc.habs(p): fill(h, rng, h.name)
    return p


def build_meta(s, cls):
    rng = random.Random(s); np.random.seed(s)
    m = mc.make_islands(rng.randint(2, 4), rng.randint(1, 3), L=rng.randint(1, 4), W=rng.randint(1, 4), K=10, cls=cls)
    for key, h in C.habs_of(m): fill(h, rng, '%s.%s' % key[1:])
    return m


def compare(label, build, v3_cls, call_cur, call_v3, sources, cur_cls, prep=None):
    bad, first, st = Counter(), {}, Counter()
    for s in SEEDS:
        cur, old = build(s, cur_cls), build(s, v3_cls)
        if prep: prep(cur, s); prep(old, s)
        if C.layout(cur) != C.layout(old): bad['twins'] += 1; first.setdefault('twins', s); continue
        pools0, grid0 = C.pool_ids(cur), C.grid_ids(cur)
        del C.RECYCLED[:]; del C.FAKE_COPIES[:]
        C.seed_all(10 ** 6 + s); call_cur(cur, s); r_cur = C.rng_state()
        rec, made = list(C.RECYCLED), list(C.FAKE_COPIES)
        C.seed_all(10 ** 6 + s); call_v3(old, s); r_old = C.rng_state()
        pools1, grid1 = C.pool_ids(cur), C.grid_ids(cur)
        by_id = {id(x): x for x in C.individuals(cur)}
        new_ids = [grid1[k] for k in grid1 if k not in grid0] + [i for k in pools1 for i in pools1[k][len(pools0[k]):]]
        res = dict(zip(RES_KEYS, [
            C.layout(cur) == C.layout(old),
            C.same_rng_state(r_cur, r_old),
            all(pools1[k] == pools0[k] for k in pools0 if k[1] in sources),
            all(pools1[k][:len(pools0[k])] == pools0[k] for k in pools0),
            all(grid1.get(k) == v for k, v in grid0.items()),
            all(C.is_copy(by_id[i]) for i in new_ids) and set(new_ids) == {id(c) for c in made if c.is_alive} and len(new_ids) + len(rec) == len(made),
            all(lab.endswith('*') for lab in rec)]))
        for k, v in res.items():
            if not v: bad[k] += 1; first.setdefault(k, s)
        st['copies'] += len(made); st['placed'] += len(new_ids); st['recycled'] += len(rec); st['seeds with recycling'] += bool(rec)
    for k in RES_KEYS:
        C.ok('a1 %s: %s (%d seeds)' % (label, k, len(SEEDS)), bad[k] == 0 and bad['twins'] == 0,
             {'failing seeds': bad[k], 'first': first.get(k), 'twins differ': bad['twins']} if bad[k] or bad['twins'] else (dict(st) if k == RES_KEYS[0] else None))
    C.ok('a1 %s: not vacuous (copies made on %d seeds)' % (label, len(SEEDS)), st['copies'] > 0 and st['placed'] > 0, dict(st))
    return dict(st)


def force_patch(p, s):
    ''' pin the within-patch number matrix (no np.random) to random counts up to each source pool size: migrants often exceed the empty sites -> copies recycled '''
    rng = random.Random(-1 - s); hs = pc.habs(p)
    pc.force_within_matrix(p, [[0 if i == j else rng.randint(0, len(hs[i].offspring_pool) + len(hs[i].dormancy_pool)) for j in range(len(hs))] for i in range(len(hs))])


def force_meta(m, s):
    ''' same for the among-patch number matrix (instance attribute, read by both the current and the v3 body) '''
    rng = random.Random(-1 - s); ps = mc.all_patches(m)
    M = np.matrix([[0 if i == j else rng.randint(0, len(ps[i].get_patch_offspring_and_dormancy_pool())) for j in range(len(ps))] for i in range(len(ps))], dtype=float)
    m.get_dispersal_among_num_matrix = lambda *a, **k: M.copy()


def part_a1():
    C.part('copy a1: fake individuals vs v3')
    C.section('a1: fake individuals vs v3, seed by seed')
    out = {}
    V3H = C.v3_subclass(habitat, 'metaibm/habitat.py', list(GERM))
    for f, src in GERM.items():
        out[f] = compare(f, build_hab, V3H, lambda h, s, f=f: getattr(h, f)(is_remove=False), lambda h, s, f=f: getattr(h, f)(False), src, habitat)
    V3P = C.v3_subclass(patch, 'metaibm/patch.py', list(PATCH))
    for f, src in PATCH.items():
        out[f] = compare(f, build_patch, V3P, lambda p, s, f=f: getattr(p, f)(rate(s), is_remove=False), lambda p, s, f=f: getattr(p, f)(rate(s), False), src, patch)
    f = 'patch_dipersal_within_from_offspring_and_dormancy_pool'
    out[f + ' forced'] = compare(f + ' [forced overflow]', build_patch, V3P, lambda p, s: getattr(p, f)(0.5, is_remove=False), lambda p, s: getattr(p, f)(0.5, False), PATCH[f], patch, prep=force_patch)
    V3M = C.v3_subclass(metacommunity, 'metaibm/metacommunity.py', list(META))
    f = 'dispersal_among_patches_from_offspring_pool_and_dormancy_pool'
    out[f + ' forced'] = compare(f + ' [forced overflow]', build_meta, V3M, lambda m, s: getattr(m, f)(0.5, is_remove=False), lambda m, s: getattr(m, f)(0.5, is_remove=False), META[f], metacommunity, prep=force_meta)
    V3E = C.v3_subclass(metacommunity, 'extension/global_habitat_network.py', [EXT])
    for meth, kw in KERNELS:
        for f, src in META.items():
            out['%s %s' % (f, meth)] = compare('%s [%s]' % (f, meth), build_meta, V3M, lambda m, s, f=f, meth=meth, kw=kw: getattr(m, f)(rate(s), method=meth, is_remove=False, **kw),
                                               lambda m, s, f=f, meth=meth, kw=kw: getattr(m, f)(rate(s), method=meth, is_remove=False, **kw), src, metacommunity)
        out['%s %s' % (EXT, meth)] = compare('%s [%s]' % (EXT, meth), build_meta, V3E, lambda m, s, meth=meth, kw=kw: getattr(m, EXT)(rate(s), method=meth, is_remove=False, **kw),
                                             lambda m, s, meth=meth, kw=kw: getattr(m, EXT)(rate(s), method=meth, **kw), ('offspring_pool',), metacommunity)
    return out


def build_colon(s):
    rng = random.Random(s); np.random.seed(s)
    ml = mc.make_mainland(n_habs=2, L=rng.randint(4, 6))
    for key, h in C.habs_of(ml): fill(h, rng, 'ML.%s' % key[2], occ=(0.7, 1.0), pools=False)
    isl = mc.make_islands(rng.randint(1, 3), rng.randint(1, 2), L=rng.randint(2, 4), K=0)
    for key, h in C.habs_of(isl): fill(h, rng, '%s.%s' % key[1:], occ=(0.0, 0.7), pools=False)
    return ml, isl


def part_a2():
    C.part('copy a2: colonization False vs True')
    C.section('a2: colonization, is_remove=False vs True in twin worlds')
    keys = ['same picks (by origin) on the same island sites', 'random + np.random identical', 'False: mainland untouched (layout + objects)', 'False: every colonist is a copy',
            'True: exactly the picked mainland sites emptied, the rest untouched', 'True: every colonist is the mainland object itself', 'nothing recycled']
    out = {}
    for f in COLON:
        for kind in ('integer', 'fractional'):
            bad, first, st = Counter(), {}, Counter()
            for s in SEEDS:
                rain = random.Random(s).randint(1, 7) + (0.0 if kind == 'integer' else 0.37)
                (mlA, isA), (mlB, isB) = build_colon(s), build_colon(s)
                ml_layout0, ml_grid0, isl_grid0 = C.layout(mlA), C.grid_ids(mlA), C.grid_ids(isA)
                ml_lab = {(key, i, j): x.label for key, h in C.habs_of(mlB) for i, row in enumerate(h.set['microsite_individuals']) for j, x in enumerate(row) if x is not None}
                ml_objB = {x.label: x for x in C.individuals(mlB)}
                del C.RECYCLED[:]
                C.seed_all(2 * 10 ** 6 + s); getattr(isA, f)(mlA, rain, is_remove=False); rA = C.rng_state()
                C.seed_all(2 * 10 ** 6 + s); getattr(isB, f)(mlB, rain, is_remove=True); rB = C.rng_state()
                gA, gB = C.grid_ids(isA), C.grid_ids(isB)
                objA = {id(x): x for x in C.individuals(isA)}; objB = {id(x): x for x in C.individuals(isB)}
                colA = [objA[gA[k]] for k in gA if k not in isl_grid0]; colB = [objB[gB[k]] for k in gB if k not in isl_grid0]
                picked = {C.origin(x) for x in colB}
                mlB_now = {(key, i, j): x.label for key, h in C.habs_of(mlB) for i, row in enumerate(h.set['microsite_individuals']) for j, x in enumerate(row) if x is not None}
                res = dict(zip(keys, [
                    C.layout(isA) == C.layout(isB) and len(colA) == len(colB),
                    C.same_rng_state(rA, rB),
                    C.layout(mlA) == ml_layout0 and C.grid_ids(mlA) == ml_grid0,
                    all(C.is_copy(x) for x in colA),
                    mlB_now == {k: v for k, v in ml_lab.items() if v not in picked} and len(picked) == len(colB),
                    all(x is ml_objB[x.label] for x in colB),
                    not C.RECYCLED]))
                for k, v in res.items():
                    if not v: bad[k] += 1; first.setdefault(k, s)
                st['colonists'] += len(colA); st['seeds'] += 1
            for k in keys:
                C.ok('a2 %s, %s rain: %s (%d seeds)' % (f, kind, k, len(SEEDS)), bad[k] == 0, {'failing seeds': bad[k], 'first': first.get(k)} if bad[k] else (dict(st) if k == keys[0] else None))
            C.ok('a2 %s, %s rain: not vacuous' % (f, kind), st['colonists'] > 0, dict(st))
            out['%s %s' % (f, kind)] = dict(st)
    return out


# ============================================================ part b: real Nemo individuals
def nemo_world(seed, F):
    C.seed_all(seed)
    ml = mc.make_mainland(L=6); isl = mc.make_islands(3, 2, L=5, K=10, gidx0=100)
    ml.meta_initialize(F, mc.TABLE); isl.meta_initialize(F, mc.TABLE)
    mc.fill_pools(isl)
    return ml, isl


def pool_members(world, pools):
    return [x for _, h in C.habs_of(world) for q in pools for x in getattr(h, q)]


def micro_members(world):
    return [x for _, h in C.habs_of(world) for row in h.set['microsite_individuals'] for x in row if x is not None]


def B_CALLS():
    hab_calls = [('habitat._germinate_from_pools via %s' % f, (lambda ml, isl, f=f: [getattr(h, f)(is_remove=False) for _, h in C.habs_of(isl)]), (lambda ml, isl, src=src: pool_members(isl, src)))
                 for f, src in GERM.items()]
    patch_calls = [('patch.%s' % f, (lambda ml, isl, f=f: [getattr(p, f)(0.5, is_remove=False) for p in mc.all_patches(isl)]), (lambda ml, isl, src=src: pool_members(isl, src)))
                   for f, src in PATCH.items()]
    meta_calls = [('metacommunity.%s [%s]' % (f, meth), (lambda ml, isl, f=f, meth=meth, kw=kw: getattr(isl, f)(0.5, method=meth, is_remove=False, **kw)), (lambda ml, isl, src=src: pool_members(isl, src)))
                  for f, src in META.items() for meth, kw in KERNELS]
    col_calls = [('metacommunity.%s (rain %s)' % (f, r), (lambda ml, isl, f=f, r=r: getattr(isl, f)(ml, r, is_remove=False)), (lambda ml, isl: micro_members(ml))) for f in COLON for r in (6, 12.4)]
    ext_calls = [('extension.%s [%s]' % (EXT, meth), (lambda ml, isl, meth=meth, kw=kw: getattr(isl, EXT)(0.5, method=meth, is_remove=False, **kw)), (lambda ml, isl: pool_members(isl, ('offspring_pool',))))
                 for meth, kw in KERNELS]
    fp, fm = 'patch_dipersal_within_from_offspring_and_dormancy_pool', 'dispersal_among_patches_from_offspring_pool_and_dormancy_pool'
    forced = [('patch.%s [forced overflow]' % fp, (lambda ml, isl: [force_patch(p, k) or getattr(p, fp)(0.5, is_remove=False) for k, p in enumerate(mc.all_patches(isl))]),
               (lambda ml, isl: pool_members(isl, PATCH[fp]))),
              ('metacommunity.%s [forced overflow]' % fm, (lambda ml, isl: force_meta(isl, 0) or getattr(isl, fm)(0.5, is_remove=False)), (lambda ml, isl: pool_members(isl, META[fm])))]
    return hab_calls + patch_calls + meta_calls + col_calls + ext_calls + forced


def part_b(n_seeds=5):
    C.part('copy b: real Nemo individuals')
    C.section('b: real Nemo individuals, is_remove=False')
    F = mc.species(7); fs = list(F.values())
    keys = ['originals stay in place with unchanged copied_state + phenotype', 'each new individual has the copied_state of an original in the source pool(s)',
            'each copy == its own source at creation (tracker), source taken from the source pool(s)', 'minted += copies, recycled += copies not placed (and they are dead)',
            'mc.check_world after the call', 'mc.teardown: minted == recycled']
    out = {}
    for label, call, sources in B_CALLS():
        bad, first, st = Counter(), {}, Counter()
        for s in range(n_seeds):
            ml, isl = nemo_world(500 + s, F); world = [ml, isl]
            held0 = C.individuals(world); loc0 = C.locations(world)
            state0 = {id(x): (C.copied_state(x), x.get_indi_phenotype_ls()) for x in held0}
            src = sources(ml, isl); src_ids = {id(x) for x in src}; src_states = {state0[id(x)][0] for x in src}
            m0, r0 = mc.minted(fs), mc.recycled(fs)
            with C.CopyTracker(states=True) as tr:
                C.seed_all(9000 + s); call(ml, isl)
            created = tr.states
            held1 = C.individuals(world); loc1 = C.locations(world)
            new = [x for x in held1 if id(x) not in state0]
            copies = [c for _, c in tr.copies]; copy_ids = {id(c) for c in copies}
            res = dict(zip(keys, [
                all(loc1.get(i) == loc0[i] for i in loc0) and all(x.is_alive and (C.copied_state(x), x.get_indi_phenotype_ls()) == state0[id(x)] for x in held0),
                all(C.copied_state(x) in src_states for x in new),
                len(created) == len(copies) and all(a == b for _, a, b in created) and all(i in src_ids for i, _, _ in created),
                mc.minted(fs) - m0 == len(copies) and mc.recycled(fs) - r0 == len(copies) - len(new) and {id(x) for x in new} <= copy_ids
                and all(not c.is_alive for c in copies if id(c) not in {id(x) for x in new}),
                C.guarded('b %s seed %d: check_world' % (label, s), mc.check_world, world, fs, label) is True,
                C.guarded('b %s seed %d: teardown' % (label, s), mc.teardown, world, fs) is not None]))
            for k, v in res.items():
                if not v: bad[k] += 1; first.setdefault(k, s)
            st['copies'] += len(copies); st['placed'] += len(new); st['recycled (not placed)'] += len(copies) - len(new)
        for k in keys:
            C.ok('b %s: %s (%d worlds)' % (label, k, n_seeds), bad[k] == 0, {'failing worlds': bad[k], 'first': first.get(k)} if bad[k] else (dict(st) if k == keys[0] else None))
        C.ok('b %s: not vacuous' % label, st['copies'] > 0 and st['placed'] > 0, dict(st))
        out[label] = dict(st)
    return out


# ====================================================================== births (B + C)
def set_rng(state):
    random.setstate(state[0]); np.random.set_state(state[1])


def residents(h):
    return {(i, j): x for i, row in enumerate(h.set['microsite_individuals']) for j, x in enumerate(row) if x is not None}


def mixed_hab(s, L=8, W=8, occ=1.0, d=0.05, gidx=7):
    ''' habitat (mean 0.2) with founders of sp1 (~70 %, fits: mostly asexual parents) and sp2 (misfits: sexual parents) on a random `occ` share of the sites, after one dead selection '''
    rng = random.Random(s); C.seed_all(s)
    h = hc.make_hab(L=L, W=W, K=0, global_index=gidx)
    for (i, j) in rng.sample([(i, j) for i in range(L) for j in range(W)], int(L * W * occ)):
        sp = 'sp1' if rng.random() < 0.7 else 'sp2'
        h.add_individual(F[sp].make_founder(rng.choice(('female', 'male')), gidx), i, j)
    h.hab_dead_selection(d, 0.5, 'niche_gaussian')
    return h


def pairs_like_code(cat):
    ''' hab_mixed_sexual_pairwise_parents_ls() replayed on a copy: lookups first, then shuffle females, shuffle males, zip '''
    out = []
    for sp, d in cat.items():
        if 'female' not in d or 'male' not in d: continue
        random.shuffle(d['female']); random.shuffle(d['male'])
        out += list(zip(d['female'], d['male']))
    return out


def is_clone_of(k, mother, gidx):
    n = k.nemo_indi_obj
    return (k.mother_id == mother.individual_id and k.father_id == 0 and n.getPedigreeClass() == 4 and n.getHome() == gidx and k.gender == mother.gender
            and k.species_id == mother.species_id and k.indfactory is mother.indfactory)


# ============================================================ 1. with_num
def sec_with_num(n_seeds=150):
    C.part('births 1 with_num')
    C.section('1. with_num: parents == random.sample of the same RNG state')
    bad, first, st = Counter(), {}, Counter()
    def fail(k, s):
        bad[k] += 1; first.setdefault(k, s)
    for s in range(n_seeds):
        h = mixed_hab(s); res = residents(h); g = h.global_index
        for kind, lst, fn in (('asex', h.occupied_site_pos_ls, h.hab_asex_reproduce_mutate_with_num), ('mix asex', h.asexual_parent_pos_ls, h.hab_mix_asex_reproduce_mutate_with_num)):
            num = random.Random(s).randint(0, len(lst))
            s0 = C.rng_state(); exp = random.sample(list(lst), num); s_exp = C.rng_state(); set_rng(s0)
            m0 = F['sp1'].minted_num + F['sp2'].minted_num
            kids = fn(None, None, num)
            if [k.mother_id for k in kids] != [res[p].individual_id for p in exp]: fail('%s: mothers == random.sample(%s, num)' % (kind, 'occupied_site_pos_ls' if kind == 'asex' else 'asexual_parent_pos_ls'), s)
            if not C.same_rng_state(C.rng_state(), s_exp): fail('%s: Python RNG consumed exactly as random.sample' % kind, s)
            if not all(is_clone_of(k, res[p], g) for k, p in zip(kids, exp)): fail('%s: clone = father 0, pedigree 4, home = global_index, gender + species of the mother' % kind, s)
            if F['sp1'].minted_num + F['sp2'].minted_num - m0 != num: fail('%s: minted += num' % kind, s)
            st[kind] += len(kids); st[kind + ' male mothers'] += sum(k.gender == 'male' for k in kids)
            for k in kids: k.recycle()
        cat0 = copy.deepcopy(h.species_category_for_sexual_parents_pos)
        n_pairs = sum(min(len(d.get('female', [])), len(d.get('male', []))) for d in cat0.values() if 'female' in d and 'male' in d)
        num = random.Random(-s).randint(0, n_pairs)
        s0 = C.rng_state(); cat = copy.deepcopy(cat0); exp = random.sample(pairs_like_code(cat), num); genders = [random.sample(('male', 'female'), 1)[0] for _ in exp]
        s_exp = C.rng_state(); set_rng(s0)
        kids = h.hab_mix_sex_reproduce_mutate_with_num(None, None, num)
        if [(k.mother_id, k.father_id) for k in kids] != [(res[f].individual_id, res[m].individual_id) for f, m in exp]: fail('mix sex: (mother, father) == random.sample(pairs rebuilt like the code, num)', s)
        if [k.gender for k in kids] != genders: fail('mix sex: offspring genders == the random.sample(("male", "female"), 1) draws', s)
        if not C.same_rng_state(C.rng_state(), s_exp) or h.species_category_for_sexual_parents_pos != cat: fail('mix sex: RNG + in-place shuffle of species_category_for_sexual_parents_pos as replayed', s)
        legal = {(res[f].individual_id, res[m].individual_id) for sp, d in cat0.items() for f in d.get('female', []) for m in d.get('male', [])}
        if not all((k.mother_id, k.father_id) in legal and res_sp_ok(k, res, cat0) for k in kids): fail('mix sex: every (mother, father) a same-species female x male pair of species_category_for_sexual_parents_pos', s)
        if not all(k.nemo_indi_obj.getPedigreeClass() != 4 and k.nemo_indi_obj.getHome() == g for k in kids): fail('mix sex: pedigree class != 4, home = global_index', s)
        st['mix sex'] += len(kids)
        for k in kids: k.recycle()
        if not C.guarded('with_num seed %d: check_held' % s, C.check_held, h, FS): fail('check_held', s)
        h.hab_clear_up_microsite_individuals()
    keys = ['asex: mothers == random.sample(occupied_site_pos_ls, num)', 'asex: Python RNG consumed exactly as random.sample', 'asex: clone = father 0, pedigree 4, home = global_index, gender + species of the mother',
            'asex: minted += num', 'mix asex: mothers == random.sample(asexual_parent_pos_ls, num)', 'mix asex: Python RNG consumed exactly as random.sample',
            'mix asex: clone = father 0, pedigree 4, home = global_index, gender + species of the mother', 'mix asex: minted += num',
            'mix sex: (mother, father) == random.sample(pairs rebuilt like the code, num)', 'mix sex: offspring genders == the random.sample(("male", "female"), 1) draws',
            'mix sex: RNG + in-place shuffle of species_category_for_sexual_parents_pos as replayed', 'mix sex: every (mother, father) a same-species female x male pair of species_category_for_sexual_parents_pos',
            'mix sex: pedigree class != 4, home = global_index', 'check_held']
    for k in keys: C.ok('1 %s (%d seeds)' % (k, n_seeds), bad[k] == 0, {'failing seeds': bad[k], 'first': first.get(k)} if bad[k] else (dict(st) if k == keys[0] else None))
    C.ok('1 not vacuous: clones from male mothers too, sexual offspring born', st['asex'] > 0 and st['mix asex'] > 0 and st['mix sex'] > 0 and st['asex male mothers'] > 0, dict(st))
    extra = set(bad) - set(keys)
    if extra: C.ok('1 unexpected failure keys', False, sorted(extra))


def res_sp_ok(k, res, cat):
    f = [p for d in cat.values() for p in d.get('female', []) if res[p].individual_id == k.mother_id]
    m = [p for d in cat.values() for p in d.get('male', []) if res[p].individual_id == k.father_id]
    return len(f) == 1 and len(m) == 1 and res[f[0]].species_id == res[m[0]].species_id == k.species_id and res[f[0]].gender == 'female' and res[m[0]].gender == 'male'


# ============================================================ 2. asexual birth into empty sites
def sec_asexual_germinate(n_seeds=200):
    C.part('births 2 asexual birth into sites')
    C.section('2. hab_asexual_reprodece_germinate')
    bad, first, st = Counter(), {}, Counter()
    for s in range(n_seeds):
        rng = random.Random(s)
        h = mixed_hab(s, L=rng.randint(2, 7), W=rng.randint(2, 7), occ=rng.uniform(0.05, 1.0))
        rate = round(rng.uniform(0, 1.0), 3); res = residents(h); empty0, n0 = len(h.empty_site_pos_ls), h.indi_num
        m0, r0 = C.mc.minted(FS), C.mc.recycled(FS)
        born = h.hab_asexual_reprodece_germinate(rate, None, None)
        new = [x for p, x in residents(h).items() if p not in res]
        want = min(empty0, int(n0 * rate))
        ok = (born == len(new) == want == h.indi_num - n0 and C.mc.minted(FS) - m0 == want and C.mc.recycled(FS) == r0
              and all(any(is_clone_of(x, m, h.global_index) for m in res.values() if m.individual_id == x.mother_id) for x in new))
        held = C.guarded('2 seed %d: check_held' % s, C.check_held, h, FS)
        if not ok or not held: bad['x'] += 1; first.setdefault('x', (s, born, want, empty0, n0, rate))
        st['births'] += born; st['capped by empty sites'] += int(n0 * rate) > empty0
        h.hab_clear_up_microsite_individuals()
    over_rate_probe()
    C.ok('2 births == min(empty, int(indi_num x rate)) on every seed; each newborn a clone of a resident (father 0, pedigree 4, home, gender); minted += births, nothing recycled; check_held (%d seeds)' % n_seeds,
         bad['x'] == 0, {'failing': bad['x'], 'first (seed, born, want, empty, indi_num, rate)': first.get('x')} if bad['x'] else dict(st))


def over_rate_probe():
    ''' rate > 1: int(indi_num x rate) may exceed the residents -> random.sample(occupied, num) raises; is the habitat left untouched? (v3 body: the same random.sample) '''
    out = {}
    v3 = C.v3_source('metaibm/habitat.py')
    same_in_v3 = 'parent_pos_ls = random.sample(self.occupied_site_pos_ls, num)' in v3 and 'parent_pos_ls = random.sample(self.asexual_parent_pos_ls, num)' in v3
    for label, fn in (('hab_asexual_reprodece_germinate(2.0)', lambda h: h.hab_asexual_reprodece_germinate(2.0, None, None)),
                      ('hab_asex_reproduce_mutate_into_offspring_pool(2.0)', lambda h: h.hab_asex_reproduce_mutate_into_offspring_pool(2.0, None, None)),
                      ('hab_mixed_reproduce_germinate(2.0, 0.0)', lambda h: h.hab_mixed_reproduce_germinate(2.0, 0.0, None, None)),
                      ('hab_mix_reproduce_mutate_into_offspring_pool(2.0, 0.0)', lambda h: h.hab_mix_reproduce_mutate_into_offspring_pool(2.0, 0.0, None, None))):
        h = mixed_hab(4242, L=6, W=6, occ=0.3)
        h.offspring_pool = [F['sp1'].make_founder('female', 7) for _ in range(3)]                        # leftovers of a previous round
        lay0, m0, r0 = (dict((p, x.individual_id) for p, x in residents(h).items()), len(h.offspring_pool)), mc.minted(FS), mc.recycled(FS)
        try: fn(h); err = 'no error'
        except ValueError as e: err = 'ValueError: %s' % e
        held_ok = C.guarded('2 probe %s: check_held' % label, C.check_held, h, FS) is True
        out[label] = {'error': err, 'minted': mc.minted(FS) - m0, 'recycled': mc.recycled(FS) - r0, 'residents unchanged': dict((p, x.individual_id) for p, x in residents(h).items()) == lay0[0],
                      'offspring_pool': '%d -> %d' % (lay0[1], len(h.offspring_pool)), 'check_held': held_ok}
        h.hab_clear_up_offspring_pool(); h.hab_clear_up_microsite_individuals()
    C.ok('2 probe: birth rate > 1 with fewer parents than births -> ValueError from random.sample (identical in the v3 bodies); recorded, invariants hold',
         all(v['error'].startswith('ValueError') and v['check_held'] for v in out.values()) and same_in_v3, out)
    C.finding('known', 'hab_asexual_reprodece_germinate, hab_asex_ / hab_mix_reproduce_mutate_into_offspring_pool, hab_mixed_reproduce_germinate',
              'asexual (or mixed-asexual) birth rate > 1: the number of births can exceed the number of (asexual) parents -> random.sample(..., num) raises ValueError '
              '(one offspring per parent per call, as for sexual_birth_rate > 1 / TODO(offspring_num)); the v3 bodies sample the same way. experiments/ use 0.5.', out)


# ============================================================ 3. mixed birth into empty sites
def mixed_expectations(h, ra, rs):
    pairs = sum(min(len(d['female']), len(d['male'])) for d in h.species_category_for_sexual_parents_pos.values() if 'female' in d and 'male' in d)
    return int(np.around(len(h.asexual_parent_pos_ls) * ra)), int(np.around(pairs * rs)), len(h.empty_site_pos_ls)


def mixed_numbers(h, ra, rs):
    ''' the current split: asex_num = min(a, np.around(E a / (a + s))), sex_num = min(s, E - asex_num) '''
    a_exp, s_exp, empty = mixed_expectations(h, ra, rs)
    a = min(a_exp, int(np.around(empty * a_exp / (a_exp + s_exp)))) if a_exp + s_exp > 0 else 0
    return a, min(s_exp, empty - a), empty


def v3_mixed_numbers(h, ra, rs):
    ''' the v3 split (both parts rounded on their own), for the comparison only '''
    a_exp, s_exp, empty = mixed_expectations(h, ra, rs)
    if empty < a_exp + s_exp: return int(np.around(empty * a_exp / (a_exp + s_exp))), int(np.around(empty * s_exp / (a_exp + s_exp))), empty
    return a_exp, s_exp, empty


def sec_mixed_germinate(n_seeds=200):
    C.part('births 3 mixed birth into sites')
    C.section('3. hab_mixed_reproduce_germinate')
    bad, first, st = Counter(), {}, Counter()
    for s in range(n_seeds):
        rng = random.Random(s)
        h = mixed_hab(s, L=rng.randint(2, 7), W=rng.randint(2, 7), occ=rng.uniform(0.3, 1.0), d=rng.choice((0.02, 0.05, 0.2)))
        ra, rs = round(rng.uniform(0, 1.0), 3), round(rng.uniform(0, 1.0), 3)
        a_exp, s_exp, _ = mixed_expectations(h, ra, rs); a, x, empty = mixed_numbers(h, ra, rs); res = residents(h)
        m0, r0 = C.mc.minted(FS), C.mc.recycled(FS)
        born = h.hab_mixed_reproduce_germinate(ra, rs, None, None)
        new = [y for p, y in residents(h).items() if p not in res]
        clones = sum(1 for y in new if y.father_id == 0)
        want = min(empty, a_exp + s_exp)
        ok = (born == len(new) == want == a + x and C.mc.minted(FS) - m0 == born and C.mc.recycled(FS) - r0 == 0 and clones == a
              and all(y.nemo_indi_obj.getHome() == h.global_index for y in new))
        held = C.guarded('3 seed %d: check_held' % s, C.check_held, h, FS)
        if not ok or not held: bad['x'] += 1; first.setdefault('x', (s, born, want, a, x, empty, clones, C.mc.recycled(FS) - r0))
        st['births'] += born; st['seeds with a + s > empty'] += a_exp + s_exp > empty
        h.hab_clear_up_microsite_individuals()
    C.ok('3 births == min(empty, a + s), minted == births, nothing recycled, clones placed == asex_num, home = global_index, check_held (%d seeds)' % n_seeds,
         bad['x'] == 0, {'failing': bad['x'], 'first (seed, born, want, asex, sex, empty, clones, recycled)': first.get('x')} if bad['x'] else dict(st))
    # deterministic grid: E empty sites x A asexual parents x P sexual pairs (rates 1.0), every A + P > E -- reaches the np.around ties (x.5) both ways
    bad, first, st, over, under = Counter(), {}, Counter(), [], []
    for E in range(1, 10):
        for A in range(1, 7):
            for P in range(1, 7):
                if A + P <= E: continue
                C.seed_all(E * 100 + A * 10 + P)
                h = hc.make_hab(L=6, W=5, K=0, global_index=11); sites = [(i, j) for i in range(6) for j in range(5)]
                for k, (i, j) in enumerate(sites[:30 - E]): h.add_individual(F['sp1'].make_founder('female' if k % 2 else 'male', 11), i, j)
                res = residents(h); fem = [p for p in res if res[p].gender == 'female']; mal = [p for p in res if res[p].gender == 'male']
                h.asexual_parent_pos_ls = list(res)[:A]
                h.species_category_for_sexual_parents_pos = {'sp1': {'female': fem[-P:], 'male': mal[-P:]}}
                a, x, empty = mixed_numbers(h, 1.0, 1.0); va, vx, _ = v3_mixed_numbers(h, 1.0, 1.0)
                m0, r0 = C.mc.minted(FS), C.mc.recycled(FS)
                born = h.hab_mixed_reproduce_germinate(1.0, 1.0, None, None)
                new = [y for p, y in residents(h).items() if p not in res]
                clones = sum(1 for y in new if y.father_id == 0)
                ok = born == len(new) == E == a + x and C.mc.minted(FS) - m0 == born and C.mc.recycled(FS) - r0 == 0 and clones == a and h.empty_site_pos_ls == []
                if not ok or not C.guarded('3 grid E=%d A=%d P=%d: check_held' % (E, A, P), C.check_held, h, FS): bad['x'] += 1; first.setdefault('x', (E, A, P, born, a, x, clones))
                if va + vx > E: over.append((E, A, P, (va, vx), (a, x)))
                if va + vx < E: under.append((E, A, P, (va, vx), (a, x)))
                st['configs'] += 1; st['births'] += born
                h.hab_clear_up_microsite_individuals()
    C.ok('3 grid E 1..9 x A 1..6 x P 1..6 (A + P > E, rates 1.0): births == E (every empty site filled), minted == births, nothing recycled, clones placed == asex_num',
         bad['x'] == 0 and over and under, {'failing': bad['x'], 'first (E, A, P, born, asex, sex, clones)': first.get('x')} if bad['x'] else dict(st, v3_overshoot_configs=len(over), v3_undershoot_configs=len(under)))
    C.finding('info', 'hab_mixed_reproduce_germinate', 'split changed vs v3 (intended): when the empty sites cannot take both expectations, v3 rounded E*a/(a+s) and E*s/(a+s) '
              'each with np.around (round half to even), so an x.5 tie made E+1 offspring (one without a site) or only E-1 (one site left empty). Now asex_num = '
              'min(a, np.around(E*a/(a+s))) and sex_num = min(s, E - asex_num): exactly min(E, a + s) offspring, all placed, nothing recycled; on a tie the asexual '
              'part keeps the round-half-even result and the sexual part takes the rest.',
              {'v3 overshoot -> now (E, A, P, v3 (asex, sex), now (asex, sex)) e.g.': over[:4], 'n_overshoot': len(over),
               'v3 undershoot -> now e.g.': under[:4], 'n_undershoot': len(under)})
    # the former D3 case: 3 empty sites, asexual expectation 3, sexual expectation 3 (v3: np.around(1.5) = 2 + 2 = 4 for 3 sites)
    C.seed_all(77)
    h = hc.make_hab(L=3, W=4, K=0, global_index=9)
    sites = [(i, j) for i in range(3) for j in range(4)]
    for k, (i, j) in enumerate(sites[:9]): h.add_individual(F['sp1'].make_founder('female' if k < 6 else 'male', 9), i, j)
    h.asexual_parent_pos_ls = sites[:3]
    h.species_category_for_sexual_parents_pos = {'sp1': {'female': sites[3:6], 'male': sites[6:9]}}
    a, x, empty = mixed_numbers(h, 1.0, 1.0)
    m0, r0 = C.mc.minted(FS), C.mc.recycled(FS); res = residents(h)
    born = h.hab_mixed_reproduce_germinate(1.0, 1.0, None, None)
    new = [y for p, y in residents(h).items() if p not in res]
    C.ok('3 tie case: 3 empty, asexual 3 x 1.0, sexual 3 pairs x 1.0 -> 2 clones + 1 sexual: 3 made, 3 born, 0 recycled (v3: 2 + 2 = 4 made)',
         (a, x, empty) == (2, 1, 3) and born == 3 and C.mc.minted(FS) - m0 == 3 and C.mc.recycled(FS) - r0 == 0 and sum(y.father_id == 0 for y in new) == 2 and h.empty_site_pos_ls == [],
         {'asex, sex, empty': (a, x, empty), 'born': born, 'minted': C.mc.minted(FS) - m0, 'recycled': C.mc.recycled(FS) - r0, 'clones placed': sum(y.father_id == 0 for y in new)})
    C.guarded('3 tie case: check_held', C.check_held, h, FS) and C.ok('3 tie case: check_held', True)
    h.hab_clear_up_microsite_individuals()


# ============================================================ 4. into offspring_pool
def sec_into_pool(n_seeds=400):
    C.part('births 4 into offspring_pool + D wrappers')
    C.section('4. hab_asex_ / hab_mix_reproduce_mutate_into_offspring_pool')
    h = mixed_hab(1234, L=8, W=8, d=0.05)
    n, A = h.indi_num, len(h.asexual_parent_pos_ls)
    P = sum(min(len(d['female']), len(d['male'])) for d in h.species_category_for_sexual_parents_pos.values() if 'female' in d and 'male' in d)
    ra = round((math.floor(n * 0.4) + 0.7) / n, 6); xa = n * ra                     # asex expectation ~ k + 0.7
    ma, ms = round((math.floor(A * 0.5) + 0.35) / A, 6), round((math.floor(P * 0.8) + 0.6) / max(P, 1), 6)
    xma, xms = A * ma, P * ms
    res = residents(h); res_ids = {x.individual_id for x in res.values()}
    for label, call, parts in (('asex', lambda: h.hab_asex_reproduce_mutate_into_offspring_pool(ra, None, None), [('clones', xa)]),
                               ('mix', lambda: h.hab_mix_reproduce_mutate_into_offspring_pool(ma, ms, None, None), [('clones', xma), ('sexual', xms)])):
        sizes = {p: [] for p, _ in parts}; bad, first = Counter(), {}
        prev = list(h.offspring_pool)
        for s in range(n_seeds):
            C.seed_all(s)
            r = call()
            got = {'clones': [y for y in h.offspring_pool if y.father_id == 0], 'sexual': [y for y in h.offspring_pool if y.father_id != 0]}
            for p, x in parts:
                k = len(got[p]); sizes[p].append(k)
                if k not in (math.floor(x), math.ceil(x)): bad['size'] += 1; first.setdefault('size', (s, p, k, x))
            if r != len(h.offspring_pool) or (label == 'asex' and got['sexual']): bad['return'] += 1; first.setdefault('return', s)
            if any(y.is_alive for y in prev): bad['leftover alive'] += 1; first.setdefault('leftover alive', s)
            if not all(y.mother_id in res_ids for y in h.offspring_pool): bad['parents'] += 1; first.setdefault('parents', s)
            if label == 'mix' and not all(y.mother_id in {res[p].individual_id for p in h.asexual_parent_pos_ls} for y in got['clones']): bad['parents'] += 1; first.setdefault('parents', s)
            if not C.guarded('4 %s seed %d: check_held' % (label, s), C.check_held, h, FS): bad['held'] += 1; first.setdefault('held', s)
            prev = list(h.offspring_pool)
        C.ok('4 %s: pool size per part only floor / ceil of the expectation, return == len(offspring_pool), parents are residents (mix clones: asexual parents) (%d seeds)' % (label, n_seeds),
             bad['size'] == 0 and bad['return'] == 0 and bad['parents'] == 0, {k: (bad[k], first.get(k)) for k in bad} or None)
        C.ok('4 %s: the previous round\'s offspring are all recycled; check_held after every round' % label, bad['leftover alive'] == 0 and bad['held'] == 0, {k: (bad[k], first.get(k)) for k in bad} or None)
        for p, x in parts:
            f = x - math.floor(x); mean = sum(sizes[p]) / n_seeds; se = math.sqrt(f * (1 - f) / n_seeds)
            C.ok('4 %s %s: mean %.4f vs expectation %.4f within 4 SE (SE %.4f)' % (label, p, mean, x, se), abs(mean - x) <= 4 * se + 1e-12, {'floor/ceil counts': dict(Counter(sizes[p]))})
    h.hab_clear_up_offspring_pool(); h.hab_clear_up_microsite_individuals()
    C.ok('4 balance after the habitat is cleared: minted == recycled', C.mc.minted(FS) == C.mc.recycled(FS), (C.mc.minted(FS), C.mc.recycled(FS)))
    # D: the patch / metacommunity wrappers of the asexual / mixed births on a whole metacommunity (mainland-style: no dormancy)
    C.seed_all(55)
    ml = mc.make_mainland(n_habs=2, L=10, gidx0=0); ml.meta_initialize(F, mc.TABLE, 'sexual')
    rows = []
    for name, fn in (('meta_dead_selection', lambda: ml.meta_dead_selection(0.05, 0.5, 'niche_gaussian')),
                     ('meta_mainland_asexual_birth_mutate_germinate', lambda: ml.meta_mainland_asexual_birth_mutate_germinate(0.5, None, None)),
                     ('meta_dead_selection', lambda: ml.meta_dead_selection(0.05, 0.5, 'niche_gaussian')),
                     ('meta_mainland_mixed_birth_mutate_germinate', lambda: ml.meta_mainland_mixed_birth_mutate_germinate(0.5, 1.0, None, None)),
                     ('meta_dead_selection', lambda: ml.meta_dead_selection(0.05, 0.5, 'niche_gaussian')),
                     ('meta_asex_reproduce_mutate_into_offspring_pool', lambda: ml.meta_asex_reproduce_mutate_into_offspring_pool(0.5, None, None)),
                     ('meta_mix_reproduce_mutate_into_offspring_pool', lambda: ml.meta_mix_reproduce_mutate_into_offspring_pool(0.5, 1.0, None, None)),
                     ('meta_local_germinate_from_offspring_and_immigrant_pool', lambda: ml.meta_local_germinate_from_offspring_and_immigrant_pool())):
        n0, o0 = mc.indi_n(ml), mc.pool_n(ml, ('offspring_pool',))
        log = C.guarded('4 D wrapper %s' % name, fn)
        rows.append((name, n0, mc.indi_n(ml), o0, mc.pool_n(ml, ('offspring_pool',)), C.guarded('4 D wrapper %s: check_world' % name, mc.check_world, ml, FS, name) is True))
    C.ok('4 D: patch / metacommunity asex + mix wrappers run on Nemo individuals (births / pool sizes > 0), invariants after each', all(r[-1] for r in rows) and rows[1][2] > rows[1][1] and rows[3][2] > rows[3][1] and rows[5][4] > 0 and rows[6][4] > 0,
         [r[:5] for r in rows])
    C.guarded('4 D: teardown', mc.teardown, ml, FS) is not None and C.ok('4 D: teardown minted == recycled', True)


# ============================================================ 5. marker births
def sec_markers(n_seeds=20):
    C.part('births 5 marker births (C)')
    C.section('5. marker births (asexual / mix_asexual / mix_sexual + sexual control)')
    BIRTH = {'asexual': ('p0', 'h0'), 'mix_asexual': ('p0', 'h1'), 'mix_sexual': ('p1', 'h0'), 'sexual': ('p1', 'h1')}
    for is_remove in (True, False):
        bad, first, st = Counter(), {}, Counter()
        for s in range(n_seeds):
            C.seed_all(300 + s); rng = random.Random(s)
            isl = mc.make_islands(3, 2, L=5, K=0, gidx0=100); isl.meta_initialize(F, mc.TABLE)
            isl.meta_dead_selection(0.05, 0.5, 'niche_gaussian')
            H = lambda p, hn: isl.set[p].set[hn]
            hb = {k: H(*v) for k, v in BIRTH.items()}
            hb['mix_asexual'].asexual_parent_pos_ls = list(hb['mix_asexual'].occupied_site_pos_ls)[:rng.randint(1, 8)]
            hb['mix_sexual'].species_category_for_sexual_parents_pos = {sp: {g: list(v)[:rng.randint(1, 4)] for g, v in d.items() if v} for sp, d in hb['mix_sexual'].species_category.items()}
            dest = [H('p2', 'h0'), H('p2', 'h1')]
            dest[0].hab_clear_up_microsite_individuals()                                                   # 25 empty
            for (i, j) in rng.sample(list(dest[1].occupied_site_pos_ls), 10):                              # 10 empty
                y = dest[1].set['microsite_individuals'][i][j]; dest[1].del_individual(i, j); y.recycle()
            for h in dest:
                h.offspring_marker_pool = [('p%d' % int(BIRTH[k][0][1]), BIRTH[k][1], k) for k in BIRTH for _ in range(rng.randint(1, 4))]
                h.immigrant_marker_pool = [('p%d' % int(BIRTH[k][0][1]), BIRTH[k][1], k) for k in BIRTH for _ in range(rng.randint(1, 4))]
            ids = {k: {x.individual_id: x for x in residents(h).values()} for k, h in hb.items()}
            asex_par = {residents(hb['mix_asexual'])[p].individual_id for p in hb['mix_asexual'].asexual_parent_pos_ls}
            msx = hb['mix_sexual'].species_category_for_sexual_parents_pos; rmsx = residents(hb['mix_sexual'])
            msx_pairs = {(rmsx[f].individual_id, rmsx[m].individual_id) for d in msx.values() for f in d.get('female', []) for m in d.get('male', []) if rmsx[f].species_id == rmsx[m].species_id}
            gidx2kind = {h.global_index: k for k, h in hb.items()}
            before = [(residents(h), list(h.offspring_marker_pool), list(h.immigrant_marker_pool), len(h.empty_site_pos_ls)) for h in dest]
            m0, r0 = mc.minted(FS), mc.recycled(FS)
            C.seed_all(900 + s)
            isl.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None, is_remove=is_remove)
            born_total = 0
            for h, (res0, om0, im0, e0) in zip(dest, before):
                new = [y for p, y in residents(h).items() if p not in res0]
                born_total += len(new)
                kinds = Counter()
                for y in new:
                    k = gidx2kind.get(y.nemo_indi_obj.getHome()); kinds[k] += 1
                    if k is None: bad['home = birth habitat global_index'] += 1; first.setdefault('home = birth habitat global_index', s); continue
                    if k == 'asexual': good = y.father_id == 0 and y.mother_id in ids[k] and y.nemo_indi_obj.getPedigreeClass() == 4 and y.gender == ids[k][y.mother_id].gender
                    elif k == 'mix_asexual': good = y.father_id == 0 and y.mother_id in asex_par and y.nemo_indi_obj.getPedigreeClass() == 4
                    elif k == 'mix_sexual': good = (y.mother_id, y.father_id) in msx_pairs and y.nemo_indi_obj.getPedigreeClass() != 4
                    else: good = y.mother_id in ids[k] and y.father_id in ids[k] and ids[k][y.mother_id].gender == 'female' and ids[k][y.father_id].gender == 'male' and ids[k][y.mother_id].species_id == ids[k][y.father_id].species_id
                    if not good: bad['parents (%s)' % k] += 1; first.setdefault('parents (%s)' % k, s)
                    st[k] += 1
                if len(new) != min(e0, len(om0) + len(im0)): bad['births == min(empty, markers)'] += 1; first.setdefault('births == min(empty, markers)', s)
                left = Counter(m[2] for m in h.offspring_marker_pool + h.immigrant_marker_pool); was = Counter(m[2] for m in om0 + im0)
                if is_remove and (left != was - kinds or len(h.offspring_marker_pool) + len(h.immigrant_marker_pool) != len(om0) + len(im0) - len(new)):
                    bad['is_remove=True: exactly the born markers removed'] += 1; first.setdefault('is_remove=True: exactly the born markers removed', s)
                if not is_remove and (h.offspring_marker_pool != om0 or h.immigrant_marker_pool != im0):
                    bad['is_remove=False: marker pools unchanged'] += 1; first.setdefault('is_remove=False: marker pools unchanged', s)
                st['capped'] += e0 < len(om0) + len(im0)
            if mc.minted(FS) - m0 != born_total or mc.recycled(FS) != r0: bad['minted += births, nothing recycled'] += 1; first.setdefault('minted += births, nothing recycled', s)
            if not C.guarded('5 seed %d: check_world' % s, mc.check_world, isl, FS, 'markers') or C.guarded('5 seed %d: teardown' % s, mc.teardown, isl, FS) is None:
                bad['check_world + teardown'] += 1; first.setdefault('check_world + teardown', s)
        keys = ['home = birth habitat global_index', 'parents (asexual)', 'parents (mix_asexual)', 'parents (mix_sexual)', 'parents (sexual)', 'births == min(empty, markers)',
                'is_remove=True: exactly the born markers removed' if is_remove else 'is_remove=False: marker pools unchanged', 'minted += births, nothing recycled', 'check_world + teardown']
        for k in keys:
            C.ok('5 is_remove=%s: %s (%d worlds)' % (is_remove, k, n_seeds), bad[k] == 0, {'failing': bad[k], 'first': first.get(k)} if bad[k] else (dict(st) if k == keys[0] else None))
        C.ok('5 is_remove=%s: all four marker kinds born, destination capped at least once' % is_remove, all(st[k] > 0 for k in BIRTH) and st['capped'] > 0, dict(st))
        extra = set(bad) - set(keys)
        if extra: C.ok('5 is_remove=%s: unexpected failure keys' % is_remove, False, sorted(extra))


# ============================================================ 6. probes of the error paths (recorded, not FAIL)
def sec_probes():
    C.part('births 6 probes')
    C.section('6. probes: D5 (birth habitat without parents) and stale parent lists')
    BIRTH = {'asexual': ('p0', 'h0'), 'mix_asexual': ('p0', 'h1'), 'mix_sexual': ('p1', 'h0'), 'sexual': ('p1', 'h1')}
    d5 = {}
    for kind, (pn, hn) in BIRTH.items():
        rows = []
        for sd in range(10):
            C.seed_all(700 + sd); isl = mc.make_islands(3, 2, L=4, K=0, gidx0=100); isl.meta_initialize(F, mc.TABLE); isl.meta_dead_selection(0.05, 0.5, 'niche_gaussian')
            dest = isl.set['p2'].set['h0']; dest.hab_clear_up_microsite_individuals()
            good = [('p0', 'h0', 'asexual')] * 3 if kind == 'sexual' else [('p1', 'h1', 'sexual')] * 3   # markers of a habitat that still has parents
            dest.offspring_marker_pool = good + [(pn, hn, kind)]
            isl.set[pn].set[hn].hab_clear_up_microsite_individuals()                               # the birth habitat loses every parent (e.g. a disturbance)
            m0, n0 = mc.minted(FS), dest.indi_num
            try: isl.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None); err = 'no error'
            except Exception as e: err = '%s: %s' % (type(e).__name__, e)
            rows.append((err, dest.indi_num - n0, mc.minted(FS) - m0, C.guarded('6 D5 %s seed %d: check_world' % (kind, sd), mc.check_world, isl, FS, kind) is True))
            C.guarded('6 D5 %s seed %d: teardown' % (kind, sd), mc.teardown, isl, FS)
        d5[kind] = {'errors': dict(Counter(r[0] for r in rows)), 'born before the error (per seed)': [r[1] for r in rows], 'minted == born': all(r[1] == r[2] for r in rows), 'check_world': all(r[3] for r in rows)}
    C.ok('6 D5 probe: a marker whose birth habitat has no parent -> ValueError for every kind, nothing minted for it, invariants hold',
         all(list(v['errors']) == ['ValueError: Sample larger than population or is negative'] and v['minted == born'] and v['check_world'] for v in d5.values()), d5)
    C.finding('known', 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool (D5)',
              'a marker whose birth habitat has no (asexual / paired) parent any more raises ValueError from random.sample when the shuffled loop reaches it (markers handled before it are already born, see the per-seed counts); '
              'now also for the restored asexual / mix_asexual / mix_sexual markers. Not reachable in the experiments/model*.py order; v3 behaviour, not fixed in this round.', d5)
    # stale parent list: a parent leaves after dead selection (del_individual, as a colonization with is_remove=True does on a mainland), then a mixed birth runs
    out, lost_total = {}, 0
    for label, call in (('hab_mix_asex_reproduce_mutate_with_num(all asexual parents)', lambda h: h.hab_mix_asex_reproduce_mutate_with_num(None, None, len(h.asexual_parent_pos_ls))),
                        ('hab_mixed_reproduce_germinate(1.0, 0.0)', lambda h: h.hab_mixed_reproduce_germinate(1.0, 0.0, None, None))):
        h = mixed_hab(808, L=8, W=8, d=0.05)
        stale = h.asexual_parent_pos_ls[len(h.asexual_parent_pos_ls) // 2]
        gone = h.set['microsite_individuals'][stale[0]][stale[1]]; h.del_individual(*stale)
        m0, r0 = mc.minted(FS), mc.recycled(FS)
        unheld0 = (m0 - r0) - len(C.individuals(h)) - 1                                           # 1 = the parent that left
        random.seed(1)
        try: call(h); err = 'no error'
        except Exception as e: err = '%s: %s' % (type(e).__name__, e)
        lost = (mc.minted(FS) - mc.recycled(FS)) - len(C.individuals(h)) - 1 - unheld0; lost_total += lost
        out[label] = {'error': err, 'clones minted before the error': mc.minted(FS) - m0, 'minted but neither held nor recycled': lost}
        gone.recycle(); h.hab_clear_up_microsite_individuals()
    C.ok('6 stale asexual_parent_pos_ls probe recorded', True, out)
    if any(v['error'] != 'no error' for v in out.values()):
        C.finding('info', 'hab_mix_asex_reproduce_mutate_with_num, hab_mixed_reproduce_germinate',
                  'asexual_parent_pos_ls is only rebuilt by hab_dead_selection: if a parent leaves its site in between (del_individual), the next mixed birth may pick the empty '
                  'site -> AttributeError on None (v3: deepcopy(None) -> AttributeError too). Error path only: the clones already made in that call are lost (minted, never held or '
                  'recycled), which matters only if a caller catches the error and continues. Not reachable in the experiments/ order (dead selection right before every birth).', out)
    C.ok('6 after the probes (everything cleared): minted - recycled == the clones lost on the error path above', mc.minted(FS) - mc.recycled(FS) == lost_total,
         {'minted - recycled': mc.minted(FS) - mc.recycled(FS), 'lost': lost_total})


CORE_PARTS = ('copy', 'births')


def run_unit(wanted):
    for name in wanted: UNIT_PARTS[name]()
    return C.summary('test_unit' + ('' if len(wanted) == len(UNIT_PARTS) else '_' + '_'.join(wanted)), EXTRA)


def run_core(wanted):
    extra = {}
    if 'copy' in wanted: extra['a1'] = part_a1(); extra['a2'] = part_a2(); extra['b'] = part_b()
    if 'births' in wanted:
        births_species()
        sec_with_num(); sec_asexual_germinate(); sec_mixed_germinate(); sec_into_pool(); sec_markers(); sec_probes()
    return C.summary('test_core' + ('' if len(wanted) == len(CORE_PARTS) else '_' + wanted[0]), extra)


if __name__ == '__main__':
    if len(sys.argv) > 3 and sys.argv[1] == '--worker':
        out = scenario(sys.argv[2], int(sys.argv[3]))
        C.worker_reply({'digest': hashlib.sha1(json.dumps(out).encode()).hexdigest(), 'records': out})
        raise SystemExit(0)
    if sys.argv[1:2] == ['unit']: raise SystemExit(run_unit([a for a in sys.argv[2:] if a in UNIT_PARTS] or list(UNIT_PARTS)))
    if sys.argv[1:2] == ['core']: raise SystemExit(run_core([a for a in sys.argv[2:] if a in CORE_PARTS] or list(CORE_PARTS)))
    raise SystemExit('usage: python test_core.py {unit|core} [part ...]  |  python test_core.py --worker {native|imitation} SEED')
