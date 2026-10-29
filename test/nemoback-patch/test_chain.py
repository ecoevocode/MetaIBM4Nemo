# Call chain patch -> habitat -> indfactory -> individual, instrumented by wrapping methods in this test process only (nothing in metaibm/ is edited):
#   1 per-function call chain (observed matrix + expected edges)   2 ordering rule: bookkeeping before recycle(), nothing touches a recycled individual
#   3 make_offspring(home_habitat) == global_index of the habitat that owns both parents   4 static scan of patch calls in metacommunity.py / simulator.py / extension/*.py
import os, re, ast, sys, glob, json, random, inspect, functools, warnings
import numpy as np
import common as C
from common import patch, habitat, indfactory, individual, make_patch, habs, held, check_invariants, recycle_everything, force_within_matrix, nemo_species, seed_all, finding, save, REPO, POOLS

warnings.filterwarnings('ignore', category=RuntimeWarning); warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
TABLE = [[0.2, 0.2], [0.5, 0.5]]
STACK, EDGES, CALLS, VIOLATIONS, HOME_CHECKS = [], {}, {}, [], []
WATCH, RECYCLED_REFS = [], []          # patches whose holdings are scanned at every recycle(); strong refs to recycled individuals (no id() reuse)
ROOT = [None]                          # outermost patch method of the current call


def holds(x):
    for p in WATCH:
        for h in habs(p):
            if any(y is x for row in h.set['microsite_individuals'] for y in row): return '%s.%s microsite' % (p.name, h.name)
            for q in POOLS:
                if any(y is x for y in getattr(h, q)): return '%s.%s %s' % (p.name, h.name, q)
    return None


def owner_of(x):
    for p in WATCH:
        for h in habs(p):
            if any(y is x for row in h.set['microsite_individuals'] for y in row): return h
    return None


def pre_hook(node, args):
    if node == 'individual.recycle':
        x = args[0]
        where = holds(x)
        if where: VIOLATIONS.append('recycle() while still held in %s (root %s)' % (where, ROOT[0]))
        if not x.is_alive: VIOLATIONS.append('recycle() of an already recycled individual (root %s)' % ROOT[0])
    elif node.startswith('individual.') and not args[0].is_alive:
        VIOLATIONS.append('%s touched a recycled individual (root %s)' % (node, ROOT[0]))
    elif node == 'indfactory.make_offspring':
        f, mother, father = args[0], args[1], args[2]
        home = args[4] if len(args) > 4 else ARGS_KW[0].get('home_habitat', 0)
        hm, hf = owner_of(mother), owner_of(father)
        HOME_CHECKS.append((ROOT[0], home, None if hm is None else hm.global_index, hm is hf and hm is not None, mother.indfactory is f, mother.species_id == father.species_id))


ARGS_KW = [{}]


def post_hook(node, args):
    if node == 'individual.recycle': RECYCLED_REFS.append(args[0])


def wrap(cls, name, layer):
    raw = cls.__dict__[name]
    node = '%s.%s' % (layer, name)
    if isinstance(raw, property):
        fget = raw.fget

        @functools.wraps(fget)
        def g(self):
            return _call(node, fget, (self,), {})
        setattr(cls, name, property(g, raw.fset, raw.fdel, raw.__doc__))
        return
    if isinstance(raw, (staticmethod, classmethod)): return
    fn = raw

    @functools.wraps(fn)
    def w(*a, **kw):
        return _call(node, fn, a, kw)
    setattr(cls, name, w)


def _call(node, fn, a, kw):
    parent = STACK[-1] if STACK else None
    if parent is None and node.startswith('patch.'): ROOT[0] = node
    root = ROOT[0]
    if root is not None:
        EDGES.setdefault(root, set()).add((parent, node))
        CALLS.setdefault(root, {}).setdefault(node, 0); CALLS[root][node] += 1
    ARGS_KW[0] = kw
    pre_hook(node, a)
    STACK.append(node)
    try: return fn(*a, **kw)
    finally:
        STACK.pop()
        post_hook(node, a)
        if not STACK: ROOT[0] = None


def instrument():
    for name, v in list(vars(habitat).items()):
        if inspect.isfunction(v): wrap(habitat, name, 'habitat')
    for name in ('make_founder', 'make_offspring', 'recycle'): wrap(indfactory, name, 'indfactory')
    for name in ('recycle', 'gender', 'get_indi_phenotype_ls', 'phenotype_set'): wrap(individual, name, 'individual')
    for name, v in list(vars(patch).items()):
        if inspect.isfunction(v) and not name.startswith('__'): wrap(patch, name, 'patch')


PUBLIC = [n for n, v in vars(patch).items() if inspect.isfunction(v) and not n.startswith('_')]
ASEX_MIX = [n for n in PUBLIC if re.search(r'asex|mix', n)]
instrument()

# ================================================================ drive every in-scope public patch method (several states, so every branch is hit)
seed_all(5); F = {'sp1': nemo_species('sp1', 0.2, seed=5), 'sp2': nemo_species('sp2', 0.5, seed=6)}; FS = list(F.values())
np.random.seed(5); p = make_patch(3, L=5, K=15, gidx0=30, means=[(0.2, 0.2), (0.5, 0.5), (0.2, 0.2)], name='p1')
np.random.seed(6); q = make_patch(2, L=4, K=5, gidx0=40, means=[(0.2, 0.2), (0.5, 0.5)], name='p2')
WATCH[:] = [p, q]
p.patch_initialize(F, TABLE); q.patch_initialize(F, TABLE)
p.patch_initialize(F, TABLE)                                  # occupied -> residents recycled first
D = ('patch_dispersal_within_from_offspring_pool_to_immigrant_pool', 'patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool', 'patch_dipersal_within_from_offspring_and_dormancy_pool')
for t in range(30):
    p.patch_dead_selection(0.2, 0.5, 'niche_gaussian')
    p.patch_sex_reproduce_mutate_into_offspring_pool(1.0)
    p.patch_sex_reproduce_mutate_into_offspring_pool(1.0)        # leftovers recycled
    p.patch_sex_reproduce_calculation_into_offspring_marker_pool(0.7)
    p.patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool(0.5, t % 2 == 0)
    getattr(p, D[t % 3])(0.6)
    [p.patch_local_germinate_from_offspring_and_dormancy_pool, p.patch_local_germinate_from_offspring_and_immigrant_pool, p.patch_local_germinate_from_offspring_immigrant_and_dormancy_pool][t % 3]()
    (p.patch_dormancy_process_from_offspring_pool_to_dormancy_pool if t % 2 else p.patch_dormancy_process_from_offspring_pool_and_immigrant_pool)()
    p.patch_clear_up_offspring_pool(); p.patch_clear_up_immigrant_pool(); p.patch_clear_up_offspring_marker_pool(); p.patch_clear_up_immigrant_marker_pool()
    if t % 7 == 6: p.patch_disturbance_process()
    p.patch_dead_selection(0.3, 0.5, 'niche_gaussian'); p.patch_sexual_birth_germinate(1.0)
    q.migrants_to_patch_into_habs_immigrant_pool(p.sample_offspring_without_replacement(3, ('offspring_pool', 'dormancy_pool')))
    q.migrants_marker_to_patch_into_habs_immigrant_marker_pool([('p1', 'h0', 'sexual')] * 3)
    check_invariants([p, q], FS, 'chain t=%d' % t)
# forced extras: OD->sites with more migrants than empty sites (individual.recycle reached from the patch method itself)
p.patch_initialize(F, TABLE); p.patch_dead_selection(0.05, 5.0, 'niche_gaussian'); p.patch_sex_reproduce_mutate_into_offspring_pool(1.0)   # few empty sites, many offspring
force_within_matrix(p, [[0, 50, 50], [50, 0, 50], [50, 50, 0]])
r0 = sum(f.recycled_num for f in FS)
p.patch_dipersal_within_from_offspring_and_dormancy_pool(0.5)
extras = sum(f.recycled_num for f in FS) - r0
del p.get_dispersal_within_num_matrix, p.patch_matrix_around
assert extras > 0
check_invariants([p, q], FS, 'forced extras')
# explicit states for the recycle branches of dormancy B / clear-ups
for _ in range(3):
    p.patch_sex_reproduce_mutate_into_offspring_pool(1.0); p.patch_dormancy_process_from_offspring_pool_and_immigrant_pool()   # O + D > K: old dormants dropped
assert sum(len(h.dormancy_pool) for h in habs(p)) > 0
p.patch_sex_reproduce_mutate_into_offspring_pool(1.0); assert p.patch_clear_up_offspring_pool() > 0
p.migrants_to_patch_into_habs_immigrant_pool([F['sp1'].make_founder('male', 30) for _ in range(6)]); assert p.patch_clear_up_immigrant_pool() == 6
assert p.patch_clear_up_dormancy_pool() > 0
check_invariants([p, q], FS, 'explicit recycle branches')
# getters / helpers
p.get_data(); p.get_patch_individual_num(); p.get_patch_empty_sites_ls(); p.patch_empty_sites_num(); p.get_patch_pairwise_empty_sites_ls(); p.get_patch_pairwise_occupied_sites_ls()
p.get_patch_occupied_sites_ls(); p.get_patch_offspring_marker_pool(); p.get_patch_offspring_pool(); p.get_patch_dormancy_pool(); p.get_patch_offspring_and_dormancy_pool()
p.get_patch_microsites_individals_sp_id_values(); p.get_patch_microsites_individals_phenotype_values('pheno1'); p.get_patch_microsites_environment_values('env1')
p.get_patch_microsites_optimum_sp_id_value_array(0.1, 0.5, TABLE); p.get_patch_dispersal_within_rate_matrix(0.3); p.get_patch_habs_offspring_num_matrix()
p.get_patch_habs_offspring_marker_num_matrix(); p.get_patch_habs_dormancy_num_matrix(); p.get_patch_habs_empty_sites_num_matrix(); p.get_habs_emigrants_matrix(0.3)
p.get_habs_immigrants_matrix(0.3); p.get_dispersal_within_num_matrix(0.3); p.patch_matrix_around(np.matrix([[0.5, 1.2]])); p.split(list(range(5)), 3)
np.random.seed(0); p.patch_offset_environmental_values(['env1'], [0.0])
np.random.seed(7); r = patch('p3', 0, (0, 0)); r.add_habitat('h0', 0, (0, 0), 2, ['env1', 'env2'], [0.2, 0.2], [0.025, 0.025], 2, 2, 0, global_index=9)
recycle_everything([p, q])
assert all(f.minted_num == f.recycled_num for f in FS)

# ================================================================ 1. chain matrix + expected edges
def reached(root, layer):
    return sorted({n.split('.', 1)[1] for (_, n) in EDGES.get(root, ()) if n.startswith(layer + '.')})


matrix = {}
for m in PUBLIC:
    root = 'patch.' + m
    matrix[m] = {'habitat': reached(root, 'habitat'), 'indfactory': reached(root, 'indfactory'), 'individual': reached(root, 'individual'),
                 'patch (internal)': [n for n in reached(root, 'patch') if n != m]}
not_called = [m for m in PUBLIC if 'patch.' + m not in EDGES and m not in ASEX_MIX]
assert not not_called, 'public in-scope patch methods never called: %s' % not_called

E = lambda root: EDGES['patch.' + root]
EXPECT = {
    'patch_dead_selection': [('patch.patch_dead_selection', 'habitat.hab_dead_selection'), ('habitat.hab_dead_selection', 'individual.get_indi_phenotype_ls'),
                             ('habitat.hab_dead_selection', 'habitat.del_individual'), ('habitat.hab_dead_selection', 'individual.recycle'), ('individual.recycle', 'indfactory.recycle')],
    'patch_initialize': [('patch.patch_initialize', 'habitat.hab_initialize'), ('habitat.hab_initialize', 'habitat.hab_clear_up_microsite_individuals'),
                         ('habitat.hab_initialize', 'indfactory.make_founder'), ('habitat.hab_initialize', 'habitat.add_individual'),
                         ('habitat.hab_clear_up_microsite_individuals', 'individual.recycle'), ('individual.recycle', 'indfactory.recycle')],
    'patch_dipersal_within_from_offspring_and_dormancy_pool': [('patch.patch_dipersal_within_from_offspring_and_dormancy_pool', 'habitat.sample_offspring_without_replacement'),
                         ('patch.patch_dipersal_within_from_offspring_and_dormancy_pool', 'habitat.add_individual'), ('patch.patch_dipersal_within_from_offspring_and_dormancy_pool', 'individual.recycle'),
                         ('individual.recycle', 'indfactory.recycle')],
    'patch_dispersal_within_from_offspring_pool_to_immigrant_pool': [('patch.patch_dispersal_within_from_offspring_pool_to_immigrant_pool', 'habitat.sample_offspring_without_replacement')],
    'patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool': [('patch.patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool', 'habitat.sample_offspring_without_replacement')],
    'patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool': [('patch.patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool', 'habitat.sample_offspring_without_replacement')],
    'patch_sex_reproduce_mutate_into_offspring_pool': [('patch.patch_sex_reproduce_mutate_into_offspring_pool', 'habitat.hab_sex_reproduce_mutate_into_offspring_pool'),
                         ('habitat.hab_sex_reproduce_mutate_into_offspring_pool', 'habitat.hab_clear_up_offspring_pool'), ('habitat.hab_clear_up_offspring_pool', 'individual.recycle'),
                         ('habitat.hab_sex_reproduce_mutate_into_offspring_pool', 'habitat.hab_sex_reproduce_mutate_with_num'), ('habitat.hab_sex_reproduce_mutate_with_num', 'indfactory.make_offspring')],
    'patch_sexual_birth_germinate': [('patch.patch_sexual_birth_germinate', 'habitat.hab_sexual_reprodece_germinate'), ('habitat.hab_sexual_reprodece_germinate', 'habitat.hab_sex_reproduce_mutate_with_num'),
                         ('habitat.hab_sex_reproduce_mutate_with_num', 'indfactory.make_offspring'), ('habitat.hab_sexual_reprodece_germinate', 'habitat.add_individual')],
    'patch_local_germinate_from_offspring_and_dormancy_pool': [('habitat.hab_local_germinate_from_offspring_and_dormancy_pool', 'habitat._germinate_from_pools'), ('habitat._germinate_from_pools', 'habitat.add_individual')],
    'patch_local_germinate_from_offspring_and_immigrant_pool': [('habitat.hab_local_germinate_from_offspring_and_immigrant_pool', 'habitat._germinate_from_pools'), ('habitat._germinate_from_pools', 'habitat.add_individual')],
    'patch_local_germinate_from_offspring_immigrant_and_dormancy_pool': [('habitat.hab_local_germinate_from_offspring_immigrant_and_dormancy_pool', 'habitat._germinate_from_pools')],
    'patch_dormancy_process_from_offspring_pool_to_dormancy_pool': [('habitat.hab_dormancy_process_from_offspring_pool_to_dormancy_pool', 'habitat._recycle_dropped'), ('habitat._recycle_dropped', 'individual.recycle')],
    'patch_dormancy_process_from_offspring_pool_and_immigrant_pool': [('habitat.hab_dormancy_process_from_offspring_pool_and_immigrant_pool', 'habitat._recycle_dropped'), ('habitat._recycle_dropped', 'individual.recycle')],
    'patch_disturbance_process': [('patch.patch_disturbance_process', 'habitat.habitat_disturbance_process')] + [('habitat.habitat_disturbance_process', 'habitat.hab_clear_up_' + s) for s in
                         ('microsite_individuals', 'offspring_pool', 'immigrant_pool', 'offspring_marker_pool', 'immigrant_marker_pool')],
    'patch_clear_up_microsite_individuals': [('patch.patch_clear_up_microsite_individuals', 'habitat.hab_clear_up_microsite_individuals'), ('habitat.hab_clear_up_microsite_individuals', 'individual.recycle')],
    'patch_clear_up_offspring_pool': [('patch.patch_clear_up_offspring_pool', 'habitat.hab_clear_up_offspring_pool'), ('habitat.hab_clear_up_offspring_pool', 'individual.recycle')],
    'patch_clear_up_immigrant_pool': [('patch.patch_clear_up_immigrant_pool', 'habitat.hab_clear_up_immigrant_pool'), ('habitat.hab_clear_up_immigrant_pool', 'individual.recycle')],
    'patch_clear_up_dormancy_pool': [('patch.patch_clear_up_dormancy_pool', 'habitat.hab_clear_up_dormancy_pool'), ('habitat.hab_clear_up_dormancy_pool', 'individual.recycle')],
    'patch_clear_up_offspring_marker_pool': [('patch.patch_clear_up_offspring_marker_pool', 'habitat.hab_clear_up_offspring_marker_pool')],
    'patch_clear_up_immigrant_marker_pool': [('patch.patch_clear_up_immigrant_marker_pool', 'habitat.hab_clear_up_immigrant_marker_pool')],
    'get_patch_microsites_individals_phenotype_values': [('patch.get_patch_microsites_individals_phenotype_values', 'individual.phenotype_set')],
    'get_patch_microsites_optimum_sp_id_value_array': [('patch.get_patch_microsites_optimum_sp_id_value_array', 'habitat.survival_rate')],
    'patch_offset_environmental_values': [('patch.patch_offset_environmental_values', 'habitat.hab_offset_environment_values'), ('habitat.hab_offset_environment_values', 'habitat.hab_reset_environment_values')],
    'get_patch_pairwise_empty_sites_ls': [('patch.get_patch_pairwise_empty_sites_ls', 'habitat.get_hab_pairwise_empty_site_pos_ls')],
    'get_patch_pairwise_occupied_sites_ls': [('patch.get_patch_pairwise_occupied_sites_ls', 'habitat.get_hab_pairwise_occupied_site_pos_ls')],
    'migrants_to_patch_into_habs_immigrant_pool': [('patch.migrants_to_patch_into_habs_immigrant_pool', 'patch.split')],
    'patch_matrix_around': [], 'add_habitat': [('patch.add_habitat', 'habitat.__init__')],
}
missing = {m: [e for e in es if e not in E(m)] for m, es in EXPECT.items()}
missing = {m: v for m, v in missing.items() if v}
assert not missing, 'expected call edges not observed: %s' % missing
# forbidden: who may mint / recycle
mint_roots = {m for m in PUBLIC if matrix[m]['indfactory'] and set(matrix[m]['indfactory']) & {'make_founder', 'make_offspring'}}
founder_roots = {m for m in PUBLIC if 'make_founder' in matrix[m]['indfactory']}
assert founder_roots == {'patch_initialize'}, founder_roots
assert mint_roots == {'patch_initialize', 'patch_sex_reproduce_mutate_into_offspring_pool', 'patch_sexual_birth_germinate'}, mint_roots
assert all(parent == 'individual.recycle' for es in EDGES.values() for parent, n in es if n == 'indfactory.recycle'), 'indfactory.recycle reached other than via individual.recycle'
no_indiv_layer = ['patch_dispersal_within_from_offspring_pool_to_immigrant_pool', 'patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool',
                  'patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool', 'patch_clear_up_offspring_marker_pool', 'patch_clear_up_immigrant_marker_pool',
                  'patch_sex_reproduce_calculation_into_offspring_marker_pool', 'sample_offspring_without_replacement', 'migrants_to_patch_into_habs_immigrant_pool']
assert all(not matrix[m]['indfactory'] and 'recycle' not in matrix[m]['individual'] for m in no_indiv_layer), {m: matrix[m] for m in no_indiv_layer}
print('1 call chains: %d public patch methods called (asex/mix probed in test_nemo.py); %d expected edge sets observed; make_founder only from patch_initialize, '
      'make_offspring only from the two sexual reproduction methods, indfactory.recycle only via individual.recycle; pool-to-pool dispersal + marker functions never reach recycle / indfactory' % (len(PUBLIC) - len(ASEX_MIX), len(EXPECT)))

# ================================================================ 2. ordering rule
assert not VIOLATIONS, VIOLATIONS[:10]
assert len(RECYCLED_REFS) == sum(f.recycled_num for f in FS) and all(not x.is_alive and x.nemo_indi_obj is None for x in RECYCLED_REFS)
print('2 ordering rule: %d recycle() calls, each on an individual already out of every microsite / pool, none touched afterwards' % len(RECYCLED_REFS))

# ================================================================ 3. homepatch of every offspring
bad = [c for c in HOME_CHECKS if not (c[3] and c[1] == c[2] and c[4] and c[5])]
assert HOME_CHECKS and not bad, bad[:5]
print('3 make_offspring: %d calls, home_habitat == global_index of the habitat owning both parents, mother\'s factory, same species (roots: %s)' %
      (len(HOME_CHECKS), sorted({c[0][6:] for c in HOME_CHECKS})))

# ================================================================ 4. static scan of callers
def scan():
    import importlib
    OP = C.original_stack()[0]
    cur = {n: inspect.signature(getattr(patch, n)) for n in PUBLIC}
    orig = {n for n, v in vars(OP).items() if inspect.isfunction(v)}
    raises_false = set(D) | {'patch_local_germinate_from_offspring_and_dormancy_pool', 'patch_local_germinate_from_offspring_and_immigrant_pool', 'patch_local_germinate_from_offspring_immigrant_and_dormancy_pool'}
    rows = []
    files = [os.path.join(REPO, 'metaibm', 'metacommunity.py'), os.path.join(REPO, 'metaibm', 'simulator.py')] + sorted(glob.glob(os.path.join(REPO, 'extension', '*.py')))
    for path in files:
        src = open(path).read(); tree = ast.parse(src)
        parents = {}
        for node in ast.walk(tree):
            for ch in ast.iter_child_nodes(node): parents[ch] = node
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)): continue
            recv, m = ast.unparse(node.func.value), node.func.attr
            if 'patch' not in recv or not (m in cur or m in orig): continue
            loc = '%s:%d' % (os.path.relpath(path, REPO), node.lineno)
            fdef = node
            while fdef in parents and not isinstance(fdef, ast.FunctionDef): fdef = parents[fdef]
            caller = fdef.name if isinstance(fdef, ast.FunctionDef) else '<module>'
            if m not in cur:
                rows.append((loc, caller, m, 'REMOVED from patch (split into the six patch_clear_up_* functions)')); continue
            try: cur[m].bind(None, *[0] * len(node.args), **{k.arg: 0 for k in node.keywords if k.arg})
            except TypeError as e:
                rows.append((loc, caller, m, 'SIGNATURE CHANGED: call with %d positional args does not bind to %s (%s)' % (len(node.args), str(cur[m]), e))); continue
            passes = [a for a in node.args if isinstance(a, ast.Name) and a.id == 'is_remove'] + [k for k in node.keywords if k.arg == 'is_remove']
            if m in raises_false and passes and isinstance(fdef, ast.FunctionDef):
                names = [a.arg for a in fdef.args.args]; defaults = dict(zip(names[len(names) - len(fdef.args.defaults):], fdef.args.defaults))
                d = defaults.get('is_remove')
                if d is not None and isinstance(d, ast.Constant) and d.value is False:
                    rows.append((loc, caller, m, 'caller default is_remove=False -> NotImplementedError with the caller\'s defaults')); continue
            if re.search(r'asex|mix', m):
                rows.append((loc, caller, m, 'asex / mix: not migrated (TypeError on Nemo individuals)')); continue
            if m == 'add_habitat' and not any(k.arg == 'global_index' for k in node.keywords) and len(node.args) < 11:
                rows.append((loc, caller, m, 'no global_index passed -> every habitat gets homepatch 0')); continue
    return rows


rows = scan()
for r in rows: print('4 %-34s %-58s %-72s %s' % r)
removed = {r[2] for r in rows if r[3].startswith('REMOVED')}
changed = {r[2] for r in rows if r[3].startswith('SIGNATURE')}
assert removed == set(), removed                                   # metacommunity.py is migrated: no call to a removed patch API left
assert changed == set(), changed                                   # ... and no call with an old signature
f_rm = [r for r in rows if 'NotImplementedError' in r[3]]
if f_rm: finding('known', 'metacommunity meta_dispersal_within_patch_* / meta_local_germinate_* default is_remove=False and pass it through, so with their defaults they now raise NotImplementedError (patch / habitat default is True): %s.' % ', '.join('%s %s' % (r[0], r[1]) for r in f_rm))
f_gi = [r for r in rows if 'global_index' in r[3]]
if f_gi: finding('known', 'metaibm/simulator.py builds habitats with patch.add_habitat(...) without global_index, so every habitat would get Nemo homepatch 0 (next round): %s.' % ', '.join(r[0] for r in f_gi))
f_am = [r for r in rows if r[3].startswith('asex')]
if f_am: finding('known', 'metacommunity still calls the unmigrated asex / mix patch functions: %s.' % ', '.join('%s %s' % (r[0], r[2]) for r in f_am))
src_mc = open(os.path.join(REPO, 'metaibm', 'metacommunity.py')).read().splitlines()
among = ['metaibm/metacommunity.py:%d' % (i + 1) for i, l in enumerate(src_mc) if re.search(r'random\.sample\(patch_i_object\.get_patch_(offspring_and_dormancy_pool|offspring_pool)\(\)', l)]
if among: finding('known', 'metacommunity among-patch dispersal (next round) with is_remove=False samples individuals with replacement from get_patch_offspring_pool() / get_patch_offspring_and_dormancy_pool() and appends the SAME objects to another patch\'s immigrant pools -> two owners (breaks the one-owner rule once individuals are Nemo-backed): %s.' % ', '.join(among))
rows += [(loc, '', 'random.sample(get_patch_..._pool())', 'among-patch is_remove=False: same individual in two patches') for loc in among]
save('chain', {'matrix': matrix, 'expected_edges': {m: [list(e) for e in es] for m, es in EXPECT.items()}, 'call_counts': {k[6:]: v for k, v in CALLS.items()},
               'recycle_calls': len(RECYCLED_REFS), 'make_offspring_calls': len(HOME_CHECKS), 'forced_extras_recycled': extras, 'static_scan': rows})
print('PASSED')
