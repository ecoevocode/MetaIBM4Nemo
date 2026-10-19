# test/nemoback-metacommunity: tests of the migrated metaibm/metacommunity.py (NemoBackend, sexual). Every part runs in its own process (Nemo keeps process-wide
# static state: RAND, Individual::currentID); helpers live in common.py. Nothing in metaibm/ or extension/ is edited.
#   python test_metacommunity.py smoke         one short run of each schedule (pool + marker, mainland colonization), invariants after every call; two fresh processes -> same digest
#   python test_metacommunity.py bookkeeping   R1 is_remove defaults, R3 else branches, R2 is_remove=False guards, R4 marker modes, clear-ups, meta_initialize, R6 tag audit, logs, helpers
#   python test_metacommunity.py nemo [SECTION ...]   equivalence with the original bodies (git) / colonization statistics / Nemo level / edges / probes (one process per section)
#   python test_metacommunity.py chain         long runs (3 schedules x 3 seeds x 200 steps, invariants after every call, RSS, C++ addresses) + fuzz (>= 3000 random ops)
#   python test_metacommunity.py all [--perf] [--no-baseline]   Step 0 baseline (existing suites, unchanged), smoke, bookkeeping, nemo, chain, coverage -> TEST_REPORT.md
import os, re, ast, sys, json, glob, math, time, random, shutil, hashlib, inspect, platform, tempfile, functools, subprocess, traceback, warnings, datetime
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import common as C
from common import (metacommunity, patch, individual, make_islands, make_mainland, rain_for, fill_pools, check_world, teardown, digest, rng_digest, full_snapshot, same_snapshot,
                    all_habs, all_patches, held, held_ids, micro, minted, recycled, indi_n, pool_n, pool_schedule, marker_schedule, species, seed_all, chi2_ok, binom_ok, rss_mb,
                    finding, save, HERE, REPO, POOLS, MARKER_POOLS, TABLE)

sys.dont_write_bytecode = True
warnings.filterwarnings('ignore', category=RuntimeWarning); warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
SELF = os.path.basename(__file__)
NOISE = ('setting random seed', '    outputs: ')
KERNELS = [('uniform', {}), ('exponential', {'rho': 2}), ('gaussian', {'sigma': 1.5}), ('cauchy', {'gamma': 1.0}), ('power_law', {'alpha': 2.0, 'r0': 1.0})]


def _spawn(args):
    ''' run this file with args in a fresh process; the worker prints 'JSON <payload>' '''
    pr = subprocess.run([sys.executable, os.path.abspath(__file__)] + list(args), capture_output=True, text=True, cwd=HERE)
    line = [l for l in pr.stdout.splitlines() if l.startswith('JSON ')]
    if not line: return {'error': 'rc=%d %s' % (pr.returncode, (pr.stdout + pr.stderr)[-2500:])}
    return json.loads(line[-1][5:])


# ############################################################################
# smoke
# ############################################################################
T0 = time.time()
SMOKE_STEPS, CURRENT = [], ['import']
N_STEPS = 6


def smoke_run(seed, record=True):
    seed_all(seed); F = species(seed + 1); fs = list(F.values())
    np.random.seed(seed); mainland = make_mainland(); islands = make_islands(3, 2, L=5, K=10, gidx0=100)
    world = [mainland, islands]
    CURRENT[0] = 'meta_initialize'
    r = mainland.meta_initialize(F, C.TABLE) + islands.meta_initialize(F, C.TABLE)
    check_world(world, fs, 'init')
    if record: SMOKE_STEPS.append({'step': 'meta_initialize (mainland + islands)', 'return': r.strip().replace('\n', ' '), 'invariants': True})
    rain = rain_for(mainland)
    trace = []
    for schedule in ('pool', 'marker'):
        for t in range(N_STEPS):
            ops = pool_schedule(mainland, islands, t, rain, dorm=t % 2 == 0) if schedule == 'pool' else marker_schedule(mainland, islands, t, rain)
            for name, fn in ops:
                CURRENT[0] = '%s t=%d %s' % (schedule, t, name)
                r = fn()
                ok = check_world(world, fs, CURRENT[0])
                if record and t == 0: SMOKE_STEPS.append({'step': '%s: %s' % (schedule, name), 'return': (r.strip().replace('\n', ' ')[:110] if isinstance(r, str) else repr(r)), 'invariants': ok})
            trace.append(digest(world))
    d = {'state': digest(world), 'rng': rng_digest(), 'trace': trace, 'minted': {f.species_id: f.minted_num for f in fs}, 'islands_alive': C.indi_n(islands)}
    CURRENT[0] = 'teardown'
    n = teardown(world, fs)
    if record: SMOKE_STEPS.append({'step': 'teardown (meta_clear_up_* on both, minted == recycled)', 'return': repr(n), 'invariants': True})
    for f in fs: f.purge_recycling_pool()
    return d


def smoke_main():
    d = smoke_run(7)
    print('[smoke] in-process run: state digest %s, rng %s, minted %s, islands alive %d' % (d['state'], d['rng'], d['minted'], d['islands_alive']))
    CURRENT[0] = 'two fresh processes, same seed'
    outs = []
    for _ in range(2):
        pr = subprocess.run([sys.executable, os.path.abspath(__file__), 'smoke', '--digest', '7'], capture_output=True, text=True, cwd=C.HERE)
        line = [l for l in pr.stdout.splitlines() if l.startswith('DIGEST ')]
        assert line, pr.stdout[-500:] + pr.stderr[-2000:]
        outs.append(json.loads(line[-1][7:]))
    assert outs[0] == outs[1], 'same seed, two processes: digests differ'
    assert outs[0]['state'] == d['state'] and outs[0]['trace'] == d['trace'], 'subprocess digest != in-process digest'
    SMOKE_STEPS.append({'step': 'two fresh processes, seed 7', 'return': 'state %s, rng %s, %d per-step digests identical' % (outs[0]['state'], outs[0]['rng'], len(outs[0]['trace'])), 'invariants': True})
    print('[smoke] two fresh processes with seed 7: identical digests (%s, %d steps)' % (outs[0]['state'], len(outs[0]['trace'])))
    return {'digest': d, 'processes': outs}


def smoke_cli():
        if len(sys.argv) > 2 and sys.argv[1] == '--digest':
            print('DIGEST ' + json.dumps(smoke_run(int(sys.argv[2]), record=False))); sys.exit(0)
        out = {'steps': SMOKE_STEPS}
        try:
            out['result'] = smoke_main(); out['ok'] = True
        except BaseException:
            out.update(ok=False, failed_step=CURRENT[0], traceback=traceback.format_exc())
            print('SMOKE FAILED at %s\n%s' % (CURRENT[0], out['traceback']))
        out['seconds'] = round(time.time() - T0, 2)
        save('smoke', out)
        if not out['ok']: sys.exit(1)
        print('SMOKE PASSED (%.2f s)' % out['seconds'])


# ############################################################################
# bookkeeping
# ############################################################################
def bookkeeping_main():
    out = {}
    SRC = open(os.path.join(REPO, 'metaibm', 'metacommunity.py')).read()
    FN = C.class_functions(SRC)
    OFN = C.class_functions(C.original_source())
    seed_all(1); F = species(1); FS = list(F.values())


    def world(seed, n_patches=3, L=4, K=10, pools=True, mainland_L=8):
        np.random.seed(seed); random.seed(seed)
        ml = make_mainland(L=mainland_L); isl = make_islands(n_patches, 2, L=L, K=K, gidx0=100)
        ml.meta_initialize(F, TABLE); isl.meta_initialize(F, TABLE)
        if pools: fill_pools(isl)
        ml.meta_dead_selection(0.2, 0.5, 'niche_gaussian')
        return ml, isl


    # ================================================================ 1. R1: every is_remove default is True
    defaults = {}
    for name, fn in vars(metacommunity).items():
        if inspect.isfunction(fn) and 'is_remove' in inspect.signature(fn).parameters:
            defaults[name] = inspect.signature(fn).parameters['is_remove'].default
    odefaults = {n: next((ast.literal_eval(d) for a, d in zip(f.args.args[len(f.args.args) - len(f.args.defaults):], f.args.defaults) if a.arg == 'is_remove'), None) for n, f in OFN.items() if any(a.arg == 'is_remove' for a in f.args.args)}
    assert defaults and all(v is True for v in defaults.values()), {k: v for k, v in defaults.items() if v is not True}
    assert set(defaults) == set(odefaults), (set(defaults) ^ set(odefaults))
    out['R1'] = {'functions': len(defaults), 'now': {k: v for k, v in defaults.items()}, 'original': odefaults}
    print('1 R1: %d metacommunity functions take is_remove, every default is True (original: %s)' % (len(defaults), sorted(set(map(str, odefaults.values())))))

    # ================================================================ 2. R3 (static): unreachable else branches are `pass  # ... TODO(copy_from)`; guards come first
    GUARDED_META = ['meta_colonize_from_propagules_rains', 'pairwise_sexual_colonization_from_prpagules_rains', 'dispersal_among_patches_from_offspring_pool_and_dormancy_pool',
                    'dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool', 'dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool']
    LINES = SRC.splitlines()
    r3 = {}
    for name, f in FN.items():
        body = f.body[1:] if ast.get_docstring(f) is not None else f.body
        guard = bool(body) and isinstance(body[0], ast.If) and ast.unparse(body[0].test) == 'not is_remove' and isinstance(body[0].body[0], ast.Raise) and 'NotImplementedError' in ast.unparse(body[0].body[0])
        ifs = [n for n in ast.walk(f) if isinstance(n, ast.If) and ast.unparse(n.test) == 'is_remove']
        if not (guard or ifs): continue
        branches = []
        for n in ifs:
            if guard:
                ok = len(n.orelse) == 1 and isinstance(n.orelse[0], ast.Pass) and 'TODO(copy_from)' in LINES[n.orelse[0].lineno - 1] and 'unreachable' in LINES[n.orelse[0].lineno - 1]
            else:
                ok = not n.orelse or not (len(n.orelse) == 1 and isinstance(n.orelse[0], ast.Pass))    # unguarded (marker) paths keep a real is_remove=False branch (or none)
            branches.append((n.lineno, ok))
        r3[name] = {'guard_first': guard, 'is_remove_branches': branches}
    assert sorted(n for n, v in r3.items() if v['guard_first']) == sorted(GUARDED_META), sorted(n for n, v in r3.items() if v['guard_first'])
    assert all(ok for v in r3.values() for _, ok in v['is_remove_branches']), {n: v for n, v in r3.items() if not all(ok for _, ok in v['is_remove_branches'])}
    n_else = sum(len(v['is_remove_branches']) for n, v in r3.items() if v['guard_first'])
    out['R3'] = r3
    print('2 R3: the 5 meta-layer guards are the first statement (after the docstring); their %d `else` branches are `pass  # unreachable ... TODO(copy_from)`; marker paths keep a live is_remove=False branch' % n_else)

    # ================================================================ 3. R2: is_remove=False raises before any RNG use / state change
    GUARDS = [('meta_colonize_from_propagules_rains', 'meta', lambda ml, isl, kw: isl.meta_colonize_from_propagules_rains(ml, rain_for(ml), **kw), lambda ml, isl: isl.meta_colonize_from_propagules_rains(ml, rain_for(ml), False)),
              ('pairwise_sexual_colonization_from_prpagules_rains', 'meta', lambda ml, isl, kw: isl.pairwise_sexual_colonization_from_prpagules_rains(ml, 6.4, **kw), lambda ml, isl: isl.pairwise_sexual_colonization_from_prpagules_rains(ml, 6.4, False)),
              ('dispersal_among_patches_from_offspring_pool_and_dormancy_pool', 'meta', lambda ml, isl, kw: isl.dispersal_among_patches_from_offspring_pool_and_dormancy_pool(0.5, method='exponential', rho=2, **kw), lambda ml, isl: isl.dispersal_among_patches_from_offspring_pool_and_dormancy_pool(0.5, 'uniform', False)),
              ('dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool', 'meta', lambda ml, isl, kw: isl.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(0.5, **kw), lambda ml, isl: isl.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(0.5, 'uniform', False)),
              ('dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool', 'meta', lambda ml, isl, kw: isl.dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool(0.5, **kw), lambda ml, isl: isl.dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool(0.5, 'uniform', False)),
              ('meta_dispersal_within_patch_from_offspring_to_immigrant_pool', 'patch', lambda ml, isl, kw: isl.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(0.5, **kw), lambda ml, isl: isl.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(0.5, False)),
              ('meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool', 'patch', lambda ml, isl, kw: isl.meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool(0.5, **kw), lambda ml, isl: isl.meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool(0.5, False)),
              ('meta_dispersal_within_patch_from_offspring_and_dormancy_pool', 'patch', lambda ml, isl, kw: isl.meta_dispersal_within_patch_from_offspring_and_dormancy_pool(0.5, **kw), lambda ml, isl: isl.meta_dispersal_within_patch_from_offspring_and_dormancy_pool(0.5, False)),
              ('meta_local_germinate_from_offspring_and_dormancy_pool', 'habitat', lambda ml, isl, kw: isl.meta_local_germinate_from_offspring_and_dormancy_pool(**kw), lambda ml, isl: isl.meta_local_germinate_from_offspring_and_dormancy_pool(False)),
              ('meta_local_germinate_from_offspring_and_immigrant_pool', 'habitat', lambda ml, isl, kw: isl.meta_local_germinate_from_offspring_and_immigrant_pool(**kw), lambda ml, isl: isl.meta_local_germinate_from_offspring_and_immigrant_pool(False)),
              ('meta_local_germinate_from_offspring_immigrant_and_dormancy_pool', 'habitat', lambda ml, isl, kw: isl.meta_local_germinate_from_offspring_immigrant_and_dormancy_pool(**kw), lambda ml, isl: isl.meta_local_germinate_from_offspring_immigrant_and_dormancy_pool(False))]
    n3, layer_seen = 0, {}
    for seed in range(12):
        for n_p in ((3, 1) if seed < 4 else (3,)):                                       # 1-patch islands: the guard must come before the patch_num < 2 early return
            ml, isl = world(seed, n_patches=n_p)
            for name, layer, kw_call, pos_call in GUARDS:
                for call in (lambda: kw_call(ml, isl, {'is_remove': False}), lambda: pos_call(ml, isl)):
                    snap = full_snapshot([ml, isl], FS)
                    try: call(); raise AssertionError('%s accepted is_remove=False' % name)
                    except NotImplementedError as e: msg = str(e)
                    assert same_snapshot(snap, full_snapshot([ml, isl], FS)), '%s seed=%d patches=%d: guard changed state / RNG / counters' % (name, seed, n_p)
                    layer_seen.setdefault(name, set()).add(('copy_from' if 'operator=' in msg else 'makeClone' if 'makeClonedOffsprg' in msg else '?', layer))
                    n3 += 1
            check_world([ml, isl], FS, 'R2 seed=%d' % seed)
            teardown([ml, isl], FS)
    assert all(len(v) == 1 for v in layer_seen.values()), layer_seen
    assert all(next(iter(v))[0] == ('copy_from' if next(iter(v))[1] == 'meta' else 'makeClone') for v in layer_seen.values()), layer_seen
    out['R2'] = {'calls': n3, 'paths': {k: sorted(v)[0] for k, v in layer_seen.items()}}
    print('3 R2: %d guarded calls (11 paths x positional + keyword, 3- and 1-patch islands); NotImplementedError before any change: state digest, random + np.random state, minted / recycled, Nemo currentID all unchanged' % n3)

    # ================================================================ 4. R4: marker paths accept is_remove False and True
    r4 = {'within True': 0, 'within False': 0, 'among True': 0, 'among False': 0, 'birth True': 0, 'birth False': 0}
    born_tot = 0
    for seed in range(15):
        for mode in (True, False):
            ml, isl = world(100 + seed, pools=False)
            isl.meta_dead_selection(0.3, 0.5, 'niche_gaussian'); isl.meta_sex_reproduce_calculation_with_offspring_marker_pool(1.0)
            o0, i0 = pool_n(isl, ('offspring_marker_pool',)), pool_n(isl, ('immigrant_marker_pool',))
            assert o0 > 0, 'no markers'
            isl.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(0.6, mode)
            o1, i1 = pool_n(isl, ('offspring_marker_pool',)), pool_n(isl, ('immigrant_marker_pool',))
            assert (o1 + i1 == o0 + i0 and o1 <= o0) if mode else (o1 == o0 and i1 >= i0), ('within', mode, o0, i0, o1, i1)
            r4['within %s' % mode] += 1
            isl.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool(0.6, method='exponential', rho=2, is_remove=mode)
            o2, i2 = pool_n(isl, ('offspring_marker_pool',)), pool_n(isl, ('immigrant_marker_pool',))
            assert (o2 + i2 == o1 + i1 and o2 <= o1) if mode else (o2 == o1 and i2 >= i1), ('among', mode, o1, i1, o2, i2)
            r4['among %s' % mode] += 1
            check_world([ml, isl], FS, 'R4 dispersal %s' % mode)
            before = {(h, mk) for h in all_habs(isl) for q in MARKER_POOLS for mk in getattr(h, q)}
            mk_n, n0, m0, empt = pool_n(isl, MARKER_POOLS), indi_n(isl), minted(F), [len(h.empty_site_pos_ls) for h in all_habs(isl)]
            exp = sum(min(e, len(h.offspring_marker_pool) + len(h.immigrant_marker_pool)) for e, h in zip(empt, all_habs(isl)))
            isl.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None, mode)
            born = indi_n(isl) - n0
            assert born == minted(F) - m0 == exp, ('birth', mode, born, exp)
            assert pool_n(isl, MARKER_POOLS) == (mk_n - born if mode else mk_n), ('birth markers', mode)
            r4['birth %s' % mode] += 1; born_tot += born
            check_world([ml, isl], FS, 'R4 birth %s' % mode)
            teardown([ml, isl], FS)
    out['R4'] = dict(r4, births=born_tot)
    # defaults: marker dispersal default False -> True. Twins without individuals (markers are plain tuples), same seed, called with their own defaults
    OM = C.original_metacommunity(['meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool', 'dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool'])
    dflt = {}
    for label, call in (('within', lambda m: m.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(0.5)),
                        ('among', lambda m: m.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool(0.5))):
        tot = {'now': [0, 0], 'v3': [0, 0]}
        for seed in range(20):
            for key, cls in (('now', metacommunity), ('v3', OM)):
                np.random.seed(seed); m = make_islands(3, 2, L=2, cls=cls)
                for h in all_habs(m): h.offspring_marker_pool = [('p0', 'h0', 'sexual')] * 10
                seed_all(700 + seed); call(m)
                tot[key][0] += pool_n(m, ('offspring_marker_pool',)); tot[key][1] += pool_n(m, ('immigrant_marker_pool',))
        assert sum(tot['now']) == 20 * 60 and tot['v3'][0] == 20 * 60 and tot['v3'][1] > 0, (label, tot)
        dflt[label] = tot
    out['R4_defaults'] = dflt
    finding('info', 'meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool, dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool',
            'marker dispersal default is_remove False -> True: with the defaults, markers now MOVE from the offspring- to the immigrant-marker pool (total conserved); v3 defaults kept every marker in the source pool AND sent a copy (more births). Both modes stay seed-identical to v3 when is_remove is passed explicitly (test_nemo equivalence).',
            evidence='; '.join('%s, 20 seeds x 60 markers (offspring, immigrant) after the call with each version\'s default: now %s, v3 %s' % (k, v['now'], v['v3']) for k, v in dflt.items()))
    print('4 R4: marker within / among / germinate-and-birth run with is_remove=True (markers move; births consume them) and False (source kept, copies sent; markers kept after birth): %s, %d births' % (r4, born_tot))

    # ================================================================ 5. clear-ups
    CLEAR = {'meta_clear_up_microsite_individuals': ('microsites',), 'meta_clear_up_offspring_pool': ('offspring_pool',), 'meta_clear_up_immigrant_pool': ('immigrant_pool',),
             'meta_clear_up_dormancy_pool': ('dormancy_pool',), 'meta_clear_up_offspring_marker_pool': ('offspring_marker_pool',), 'meta_clear_up_immigrant_marker_pool': ('immigrant_marker_pool',),
             'meta_clear_up_offspring_and_immigrant_pool': ('offspring_pool', 'immigrant_pool'), 'meta_clear_up_offspring_marker_and_immigrant_marker_pool': ('offspring_marker_pool', 'immigrant_marker_pool')}
    n5, ret_types = 0, {}
    for seed in range(8):
        for fn, what in CLEAR.items():
            ml, isl = world(200 + seed, K=(0, 5, 20)[seed % 3])
            isl.meta_sex_reproduce_calculation_with_offspring_marker_pool(0.9); isl.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(0.5)
            hs = all_habs(isl)
            content = {q: [x for h in hs for x in (micro(h) if q == 'microsites' else getattr(h, q))] for q in ('microsites',) + POOLS + MARKER_POOLS}
            kept = {q: [list(getattr(h, q)) for h in hs] for q in POOLS + MARKER_POOLS if q not in what}
            kept_micro = None if 'microsites' in what else [list(micro(h)) for h in hs]
            r0, ml_d = recycled(F), digest([ml])
            r = getattr(isl, fn)()
            exp = sum(len(content[q]) for q in what)
            ret_types[fn] = type(r).__name__
            assert r == exp, '%s seed=%d return %r != %d' % (fn, seed, r, exp)
            ind = sum(len(content[q]) for q in what if q not in MARKER_POOLS)
            assert recycled(F) - r0 == ind and all(not x.is_alive for q in what if q not in MARKER_POOLS for x in content[q]), '%s recycled count' % fn
            for q in what:
                if q == 'microsites': assert all(h.indi_num == 0 and h.occupied_site_pos_ls == [] and h.species_category == {} and len(h.empty_site_pos_ls) == h.size for h in hs)
                else: assert all(getattr(h, q) == [] for h in hs), fn + ' pool not emptied'
            for q, v in kept.items(): assert all(len(getattr(h, q)) == len(b) and all(x is y for x, y in zip(getattr(h, q), b)) for h, b in zip(hs, v)), '%s touched %s' % (fn, q)
            if kept_micro is not None: assert all(len(micro(h)) == len(b) and all(x is y for x, y in zip(micro(h), b)) for h, b in zip(hs, kept_micro)), fn + ' touched microsites'
            assert digest([ml]) == ml_d and getattr(isl, fn)() == 0, fn + ' mainland touched / second call != 0'
            check_world([ml, isl], FS, fn)
            teardown([ml, isl], FS); n5 += 1
    assert all(t == 'int' for t in ret_types.values()), ret_types
    out['clear_ups'] = {'cases': n5, 'returns': ret_types}
    print('5 clear-ups: %d cases over the 6 split functions + 2 wrappers; return == number cleared (int), recycled == former contents, only those pools emptied, '
          'meta_clear_up_offspring_and_immigrant_pool keeps the dormancy pool (identity), marker clear-ups recycle nothing, second call returns 0' % n5)

    # ================================================================ 6. meta_initialize
    np.random.seed(300); random.seed(300)
    means = [[(0.2, 0.2), (0.5, 0.5)], [(0.5, 0.5), (0.5, 0.5)], [(0.2, 0.2), (0.2, 0.2)]]
    isl = make_islands(3, 2, L=5, W=3, K=10, gidx0=400, means=means)
    m0 = minted(F)
    log = isl.meta_initialize(F, TABLE)
    assert log == 'islands initialization done! \n' and minted(F) - m0 == 90
    for h, mean in zip(all_habs(isl), [m for row in means for m in row]):
        sp = 'sp%d' % (TABLE.index(list(mean)) + 1)
        assert h.indi_num == h.size and h.empty_site_pos_ls == [] and all(x.species_id == sp and x.indfactory is F[sp] and x.nemo_indi_obj.getHome() == h.global_index for x in micro(h)), h.name
    genders = [x.gender for x in held(isl)]
    check_world(isl, FS, 'init')
    fill_pools(isl)
    residents, pools = held(isl), {(i, q): list(getattr(h, q)) for i, h in enumerate(all_habs(isl)) for q in POOLS}
    residents = [x for h in all_habs(isl) for x in micro(h)]
    m0, r0 = minted(F), recycled(F)
    isl.meta_initialize(F, TABLE)
    assert minted(F) - m0 == 90 and recycled(F) - r0 == len(residents) and not any(x.is_alive for x in residents)
    assert all(all(a is b for a, b in zip(getattr(h, q), pools[(i, q)])) and len(getattr(h, q)) == len(pools[(i, q)]) for i, h in enumerate(all_habs(isl)) for q in POOLS), 're-initialize touched pools'
    isl.meta_initialize(F, TABLE, 'asexual'); assert all(x.gender == 'female' for h in all_habs(isl) for x in micro(h))
    check_world(isl, FS, 're-init'); teardown(isl, FS)
    init_bad = {}
    for label, mns, fac in (('mean (0.9, 0.9) in p2.h1', [[(0.2, 0.2), (0.5, 0.5)], [(0.5, 0.5), (0.2, 0.2)], [(0.2, 0.2), (0.9, 0.9)]], F),
                            ("indfactory 'sp2' missing (first sp2 habitat p0.h1)", None, {'sp1': F['sp1']})):
        np.random.seed(301); isl = make_islands(3, 2, L=3, K=0, gidx0=500, means=mns)
        m0 = minted(F)
        try: isl.meta_initialize(fac, TABLE); err = 'no error'
        except Exception as e: err = type(e).__name__
        filled = ['%s.%s' % (p.name, h.name) for p in C.all_patches(isl) for h in p.set.values() if h.indi_num == h.size]
        init_bad[label] = {'error': err, 'minted_before_error': minted(F) - m0, 'filled_habitats': filled}
        check_world(isl, FS, label); teardown(isl, FS)
    assert init_bad['mean (0.9, 0.9) in p2.h1']['error'] == 'ValueError' and len(init_bad['mean (0.9, 0.9) in p2.h1']['filled_habitats']) == 5
    assert init_bad["indfactory 'sp2' missing (first sp2 habitat p0.h1)"]['error'] == 'KeyError'
    out['meta_initialize'] = {'filled': 90, 'female_of_90': genders.count('female'), 'bad_inputs': init_bad}
    finding('info', 'meta_initialize', 'no up-front validation (non-atomic by design in this version): a bad habitat raises inside patch_initialize after every earlier habitat was already filled; invariants and factory balance still hold.',
            evidence='; '.join('%s -> %s after %d founders, filled %s' % (k, v['error'], v['minted_before_error'], v['filled_habitats']) for k, v in init_bad.items()),
            fix='validate every habitat (species_2_phenotype_ls.index(mean_env_ls), indfactory key, 0 <= global_index <= 65535, reproduce_mode) before the first patch_initialize')
    print('6 meta_initialize(indfactories, species_2_phenotype_ls, reproduce_mode="sexual"): every site filled with the species of its habitat mean, getHome == global_index, '
          're-initialize recycles residents and keeps pools, asexual -> all female; bad inputs (info): %s' % {k: (v['error'], v['minted_before_error']) for k, v in init_bad.items()})

    # ================================================================ 7. R6: NemoBackend docstring tags (static) + docstring rule
    def tagged(f): return (ast.get_docstring(f) or '').startswith('NemoBackend')
    P_T = {n for n, f in C.class_functions(open(os.path.join(REPO, 'metaibm', 'patch.py')).read(), 'patch').items() if tagged(f)}
    H_T = {n for n, f in C.class_functions(open(os.path.join(REPO, 'metaibm', 'habitat.py')).read(), 'habitat').items() if tagged(f)}
    NEMO_CALLS = {'recycle', 'make_founder', 'make_offspring'}
    calls = {n: {c.func.attr for c in ast.walk(f) if isinstance(c, ast.Call) and isinstance(c.func, ast.Attribute)} for n, f in FN.items()}
    self_calls = {n: {c.func.attr for c in ast.walk(f) if isinstance(c, ast.Call) and isinstance(c.func, ast.Attribute) and ast.unparse(c.func.value) == 'self'} for n, f in FN.items()}
    expect = {n for n in FN if calls[n] & (NEMO_CALLS | P_T | H_T)}
    while True:                                                                                   # a wrapper of a tagged metacommunity method is NemoBackend too
        more = {n for n in FN if self_calls[n] & expect} - expect
        if not more: break
        expect |= more
    actual = {n for n, f in FN.items() if tagged(f)}
    why = {n: sorted(calls[n] & (NEMO_CALLS | P_T | H_T)) or ['self.' + x for x in sorted(self_calls[n] & expect)] for n in expect}
    doc_rule = {}
    for n, f in FN.items():
        if n in actual: continue
        d = ast.get_docstring(f, clean=False)
        if n in OFN: doc_rule[n] = 'kept' if d == ast.get_docstring(OFN[n], clean=False) else 'CHANGED'
        else: doc_rule[n] = 'one-line' if d is not None and '\n' not in d.strip() else 'NOT one-line'
    out['R6'] = {'expected': sorted(expect), 'actual': sorted(actual), 'why': why, 'docstring_rule': doc_rule, 'patch_tagged': sorted(P_T), 'habitat_tagged': sorted(H_T)}
    assert actual == expect, {'missing tag': sorted(expect - actual), 'extra tag': sorted(actual - expect)}
    bad_doc = {n: v for n, v in doc_rule.items() if v in ('CHANGED', 'NOT one-line')}
    if bad_doc:                                                                                   # documentation only, no runtime effect: reported, not failed
        finding('info', ', '.join(sorted(bad_doc)), 'docstring rule (untagged new functions: one-line docstring) not met; documentation only, no runtime effect.',
                evidence=bad_doc, fix="add a one-line docstring, e.g. ''' empty <pool> in every patch (markers are tuples: nothing to recycle) '''")
    print('7 R6: %d NemoBackend-tagged functions == exactly those calling recycle / make_founder / make_offspring or a tagged patch / habitat method (or a tagged metacommunity wrapper); '
          '%d untagged functions keep the original docstring, %d new untagged ones are one-line' % (len(actual), sum(v == 'kept' for v in doc_rule.values()), sum(v == 'one-line' for v in doc_rule.values())))

    # ================================================================ 8. patch_num < 2 log strings
    ONE = ['dispersal_among_patches_from_offspring_pool_and_dormancy_pool', 'dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool',
           'dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool', 'dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool']
    O1 = C.original_metacommunity(ONE)
    logs = {}
    for name in ONE:
        np.random.seed(5); a = make_islands(1, 2, L=2, name='solo'); np.random.seed(5); b = make_islands(1, 2, L=2, name='solo', cls=O1)
        st = C.rng_state()
        la = getattr(a, name)(0.5); lb = getattr(b, name)(0.5, is_remove=True)
        assert C.same_rng_state(st, C.rng_state()) and 'solo' in la and '%s' not in la and '%s' in lb and la == lb.replace('%s', 'solo', 1), (name, la, lb)
        logs[name] = la.strip()
    out['patch_num_lt_2_logs'] = logs
    finding('info', ', '.join(ONE), 'the 4 `patch_num < 2` early-return log strings now format the metacommunity name (v3 returned the literal "%s"); no RNG used, nothing changed.', evidence=list(logs.values())[0])
    print('8 patch_num < 2: 4 log strings now contain the metacommunity name (v3: literal "%%s"), e.g. %r' % list(logs.values())[0])

    # ================================================================ 9. unchanged getters / matrix helpers on Nemo individuals (sanity, same seed -> same result as the original bodies)
    import tempfile
    HELPERS = ['get_data', '__str__', 'get_meta_individual_num', 'get_meta_empty_sites_ls', 'show_meta_empty_sites_num', 'get_meta_pairwise_empty_sites_ls', 'meta_get_occupied_location_ls',
               'get_meta_pairwise_occupied_sites_ls', 'meta_offspring_pool_individual_num', 'meta_immigrant_pool_individual_num', 'meta_dormancy_pool_individual_num',
               'meta_offspring_marker_pool_marker_num', 'meta_immigrant_marker_pool_marker_num', 'matrix_around', 'dispersal_kernel_strength', 'calculate_dispersal_kernel_strength_matrix',
               'normalized_calculate_dispersal_kernel_strength_matrix', 'get_disp_among_rate_matrix', 'get_offspring_marker_pool_num_matrix', 'get_offspring_pool_num_matrix',
               'get_dormance_pool_num_matrix', 'get_patch_empty_sites_num_matrix', 'get_emigrants_matrix', 'get_immigrants_matrix', 'get_dispersal_among_num_matrix',
               'get_meta_microsites_optimum_sp_id_val', 'get_meta_microsite_environment_values', 'get_meta_microsites_individuals_sp_id_values', 'get_meta_microsites_individuals_phenotype_values',
               'columns_patch_habitat_microsites_id', 'meta_offset_environmental_values', 'update_disp_current_matrix', 'reshape_habitat_data_in_patch']
    OH = C.original_metacommunity(HELPERS)


    def helper_values(m):
        with tempfile.TemporaryDirectory() as td:
            v = [sorted(m.get_data()), len(str(m)), m.get_meta_individual_num(), m.get_meta_empty_sites_ls(), m.show_meta_empty_sites_num(), m.get_meta_pairwise_empty_sites_ls(), m.meta_get_occupied_location_ls(),
                 m.get_meta_pairwise_occupied_sites_ls(), m.meta_offspring_pool_individual_num(), m.meta_immigrant_pool_individual_num(), m.meta_dormancy_pool_individual_num(),
                 m.meta_offspring_marker_pool_marker_num(), m.meta_immigrant_marker_pool_marker_num(), m.matrix_around(np.matrix([[0.5, 1.2], [np.nan, 2.7]])).tolist(),
                 [m.dispersal_kernel_strength(1.5, k, **kw) for k, kw in (('uniform', {}), ('gaussian', {'sigma': 1.0}), ('exponential', {'rho': 2}), ('cauchy', {'gamma': 1}), ('power_law', {'alpha': 2, 'r0': 1}))],
                 m.calculate_dispersal_kernel_strength_matrix('exponential', rho=2).tolist(), m.normalized_calculate_dispersal_kernel_strength_matrix('cauchy', gamma=1).tolist(),
                 m.get_disp_among_rate_matrix(0.3, 'gaussian', sigma=1.0).tolist(), m.get_offspring_marker_pool_num_matrix().tolist(), m.get_offspring_pool_num_matrix().tolist(),
                 m.get_dormance_pool_num_matrix().tolist(), m.get_patch_empty_sites_num_matrix().tolist(), m.get_emigrants_matrix(0.3).tolist(), np.nan_to_num(m.get_immigrants_matrix(0.3)).tolist(),
                 m.get_dispersal_among_num_matrix(0.3, 'exponential', rho=2).tolist(), m.get_meta_microsites_optimum_sp_id_val(0.1, 0.5, TABLE).tolist(),
                 m.get_meta_microsite_environment_values('env1').tolist(), np.nan_to_num(m.get_meta_microsites_individuals_sp_id_values(), nan=-1).tolist(),
                 np.nan_to_num(m.get_meta_microsites_individuals_phenotype_values('pheno2'), nan=-1).tolist(), [c.tolist() for c in m.columns_patch_habitat_microsites_id()],
                 m.update_disp_current_matrix().tolist(), [a.tolist() for a in m.reshape_habitat_data_in_patch(__import__('pandas').DataFrame(m.set['p0'].get_patch_microsites_individals_sp_id_values()), 1, 2, 4, 4, 'lower')]]
            cols = m.columns_patch_habitat_microsites_id()
            df = m.meta_distribution_data_all_time_to_csv_gz(m.get_meta_microsites_individuals_sp_id_values(), os.path.join(td, 'x.csv.gz'), ['t0'], cols[2]); v.append(df.shape)
        try: m.dispersal_kernel_strength(1.0, 'bogus'); v.append('no error')
        except ValueError: v.append('ValueError')
        return v


    n9 = 0
    for seed in range(6):
        ml, isl = world(600 + seed, L=4)
        st0 = digest([ml, isl])
        seed_all(9000 + seed); a = helper_values(isl)                                 # read-only helpers: run them on the same world with the current class ...
        isl.__class__ = OH
        seed_all(9000 + seed); b = helper_values(isl)                                 # ... and with the original bodies
        isl.__class__ = metacommunity
        assert digest([ml, isl]) == st0, 'a getter changed the state'
        assert json.dumps(a, default=str) == json.dumps(b, default=str), 'helpers differ from the original seed=%d' % seed
        assert a[2] == indi_n(isl) and a[-1] == 'ValueError' and np.allclose(np.array(a[17]).sum(axis=1), 1)
        np.random.seed(seed); e1 = make_islands(3, 2, L=3); np.random.seed(seed); e2 = make_islands(3, 2, L=3, cls=OH)      # meta_offset_environmental_values writes the environment: two empty twins
        seed_all(seed); e1.meta_offset_environmental_values(['env1', 'env2'], [0.05, -0.1], [0.01, 0.0]); seed_all(seed); e2.meta_offset_environmental_values(['env1', 'env2'], [0.05, -0.1], [0.01, 0.0])
        assert all(np.array_equal(x.set[k], y.set[k]) and x.mean_env_ls == y.mean_env_ls and x.var_env_ls == y.var_env_ls for x, y in zip(all_habs(e1), all_habs(e2)) for k in ('env1', 'env2'))
        check_world([ml, isl], FS, 'helpers %d' % seed)
        teardown([ml, isl], FS); n9 += 1
    out['helpers'] = {'seeds': n9, 'functions': len(HELPERS) + 1}
    print('9 unchanged getters / matrix helpers / data IO (%d functions incl. meta_distribution_data_all_time_to_csv_gz to a temp dir): same seed -> identical values to the original bodies on %d seeds' % (len(HELPERS) + 1, n9))

    for f in FS: assert f.minted_num == f.recycled_num
    save('bookkeeping', out)
    print('PASSED')


# ############################################################################
# nemo
# ############################################################################
SECTIONS = ('equivalence', 'colonization', 'nemo_level', 'edges', 'probes')


# ------------------------------------------------------------------ statistics (no scipy in the env)
def chi2_sf(x, k):
    ''' P(chi2_k >= x): regularized upper incomplete gamma Q(k/2, x/2) (series / continued fraction, Numerical Recipes) '''
    a, x = k / 2.0, x / 2.0
    if x <= 0: return 1.0
    gl = math.lgamma(a)
    if x < a + 1:
        s, term, n = 1.0 / a, 1.0 / a, a
        for _ in range(10000):
            n += 1; term *= x / n; s += term
            if abs(term) < abs(s) * 1e-14: break
        return max(0.0, 1.0 - s * math.exp(-x + a * math.log(x) - gl))
    b, c, d = x + 1 - a, 1e300, 1 / (x + 1 - a); h = d
    for i in range(1, 10000):
        an = -i * (i - a); b += 2; d = an * d + b; d = 1 / (d if abs(d) > 1e-300 else 1e-300); c = b + an / c; c = c if abs(c) > 1e-300 else 1e-300
        delta = d * c; h *= delta
        if abs(delta - 1) < 1e-14: break
    return math.exp(-x + a * math.log(x) - gl) * h


def chi2_test(counts, expected):
    ok, stat = chi2_ok(counts, expected)
    return ok, stat, round(chi2_sf(sum((c - e) ** 2 / e for c, e in zip(counts, expected)), len(counts) - 1), 4)


def binom_p(k, n, p):
    z = (k - n * p) / math.sqrt(n * p * (1 - p))
    return round(math.erfc(abs(z) / math.sqrt(2)), 4)


# ================================================================ 1. equivalence with the original bodies
def kern(s): return KERNELS[s % len(KERNELS)]
def rate(s): return (0.1, 0.37, 0.5, 0.9, 1.0)[(s // 5) % 5]
def T(o): return {'is_remove': True} if o else {}           # originals defaulted to is_remove=False: call them with True explicitly


def among(fn, mode=None):
    def call(m, ml, s, o):
        k, kw = kern(s)
        extra = {'is_remove': mode} if mode is not None else T(o)
        return getattr(m, fn)(rate(s), method=k, **extra, **kw)
    return call


def within(fn, mode=None):
    return lambda m, ml, s, o: getattr(m, fn)(rate(s), *((mode,) if mode is not None else ((True,) if o else ())))


def germ(fn):
    return lambda m, ml, s, o: getattr(m, fn)(*((True,) if o else ()))


EQUIV = {   # key: (function, world, call(meta, mainland, seed, is_original))
    'among OD->sites': ('dispersal_among_patches_from_offspring_pool_and_dormancy_pool', 'pools', among('dispersal_among_patches_from_offspring_pool_and_dormancy_pool')),
    'among O->I': ('dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool', 'pools', among('dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool')),
    'among marker True': ('dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool', 'pools', among('dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool', True)),
    'among marker False': ('dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool', 'pools', among('dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool', False)),
    'among OD->I': ('dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool', 'pools', among('dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool')),
    'within marker True': ('meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool', 'pools', within('meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool', True)),
    'within marker False': ('meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool', 'pools', within('meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool', False)),
    'within O->I': ('meta_dispersal_within_patch_from_offspring_to_immigrant_pool', 'pools', within('meta_dispersal_within_patch_from_offspring_to_immigrant_pool')),
    'within OD->I': ('meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool', 'pools', within('meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool')),
    'within OD->sites': ('meta_dispersal_within_patch_from_offspring_and_dormancy_pool', 'pools', within('meta_dispersal_within_patch_from_offspring_and_dormancy_pool')),
    'germ O+D': ('meta_local_germinate_from_offspring_and_dormancy_pool', 'pools', germ('meta_local_germinate_from_offspring_and_dormancy_pool')),
    'germ O+I': ('meta_local_germinate_from_offspring_and_immigrant_pool', 'pools', germ('meta_local_germinate_from_offspring_and_immigrant_pool')),
    'germ O+I+D': ('meta_local_germinate_from_offspring_immigrant_and_dormancy_pool', 'pools', germ('meta_local_germinate_from_offspring_immigrant_and_dormancy_pool')),
    'germ+birth markers True': ('meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool', 'markers', lambda m, ml, s, o: m.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None, True)),
    'germ+birth markers False': ('meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool', 'markers', lambda m, ml, s, o: m.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None, False)),
    'dead selection': ('meta_dead_selection', 'pools', lambda m, ml, s, o: m.meta_dead_selection((0.1, 0.4)[s % 2], (0.5, 0.3)[s % 2], ('niche_gaussian', 'neutral_uniform')[s % 4 // 3])),
    'dormancy O->D': ('meta_dormancy_process_from_offspring_pool_to_dormancy_pool', 'pools', lambda m, ml, s, o: m.meta_dormancy_process_from_offspring_pool_to_dormancy_pool()),
    'dormancy O+I': ('meta_dormancy_process_from_offspring_pool_and_immigrant_pool', 'pools', lambda m, ml, s, o: m.meta_dormancy_process_from_offspring_pool_and_immigrant_pool()),
    'sex reproduce': ('meta_sex_reproduce_mutate_into_offspring_pool', 'pools', lambda m, ml, s, o: m.meta_sex_reproduce_mutate_into_offspring_pool(rate(s), None, None)),
    'mainland sexual birth': ('meta_mainland_sexual_birth_mutate_germinate', 'mainland', lambda m, ml, s, o: ml.meta_mainland_sexual_birth_mutate_germinate((0.5, 1.0)[s % 2], None, None)),
    'marker calc sex': ('meta_sex_reproduce_calculation_with_offspring_marker_pool', 'calc', lambda m, ml, s, o: m.meta_sex_reproduce_calculation_with_offspring_marker_pool(rate(s))),
    'marker calc asex': ('meta_asex_reproduce_calculation_into_offspring_marker_pool', 'calc', lambda m, ml, s, o: m.meta_asex_reproduce_calculation_into_offspring_marker_pool(rate(s))),
    'marker calc mix': ('meta_mix_reproduce_calculation_with_offspring_marker_pool', 'calc', lambda m, ml, s, o: m.meta_mix_reproduce_calculation_with_offspring_marker_pool(rate(s), 1 - rate(s) / 2)),
    'disturbance patches': ('meta_disturbance_process_in_patches', 'pools', lambda m, ml, s, o: m.meta_disturbance_process_in_patches((0.3, 0.6)[s % 2])),
    'disturbance habitat': ('meta_disturbance_process_in_habitat', 'pools', lambda m, ml, s, o: m.meta_disturbance_process_in_habitat((0.3, 0.6)[s % 2])),
}
N_SEEDS = {'among OD->sites': 120}
FORCED = 'among OD->sites'                                  # every 3rd seed: the among-patch number matrix is forced high (instance override) -> extras guaranteed


def build(kind, s, cls):
    ''' one world, built and prepared with the CURRENT class (identical in both variants), then switched to `cls` (the original-bodies subclass in the orig worker) '''
    rng = random.Random(s)
    seed_all(s)
    ml = make_mainland(L=8) if kind == 'mainland' else None
    isl = make_islands(3, 2, L=rng.choice((3, 4, 5)), W=rng.choice((3, 4)), K=rng.choice((0, 5, 20)), gidx0=100)
    isl.meta_initialize(EQ_F, TABLE)
    if ml is not None:
        ml.meta_initialize(EQ_F, TABLE); ml.meta_dead_selection(0.3, 0.5, 'niche_gaussian')
    if kind == 'pools': fill_pools(isl)
    elif kind == 'calc': isl.meta_dead_selection(0.3, 0.5, 'niche_gaussian')
    elif kind == 'markers':
        isl.meta_dead_selection(0.3, 0.5, 'niche_gaussian'); isl.meta_sex_reproduce_calculation_with_offspring_marker_pool(rate(s))
        isl.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(0.5); isl.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool(0.6, method='exponential', rho=2)
    isl.__class__ = cls
    if ml is not None: ml.__class__ = cls
    return ml, isl


EQ_F = None


def equiv_worker(key, variant):
    ''' one process: run `key` on N seeds with the current code (variant 'new') or the original bodies (variant 'orig'); returns per-seed records '''
    global EQ_F
    seed_all(1); EQ_F = species(1); fs = list(EQ_F.values())
    fn, kind, call = EQUIV[key]
    cls = C.original_metacommunity([fn]) if variant == 'orig' else metacommunity
    recs = []
    for s in range(N_SEEDS.get(key, 30)):
        ml, isl = build(kind, s, cls)
        world = [isl] + ([ml] if ml is not None else [])
        if key == FORCED and s % 3 == 2:
            M = np.matrix([[0 if i == j else 6 for j in range(isl.patch_num)] for i in range(isl.patch_num)], dtype=float)
            isl.get_dispersal_among_num_matrix = lambda *a, **k: M.copy()
        ids = {id(x): x.individual_id for x in held(world)}
        objs = {id(x): x for x in held(world)}
        m0, r0 = minted(EQ_F), recycled(EQ_F)
        seed_all(10000 + s)
        ret = call(isl, ml, s, variant == 'orig')
        after = held_ids(world)
        left = sorted(ids[i] for i in objs if i not in after)
        left_alive = sorted(ids[i] for i in objs if i not in after and objs[i].is_alive)
        rec = {'seed': s, 'digest': digest(world), 'rng': rng_digest(), 'ret': ret, 'left': left, 'left_alive': left_alive, 'minted': minted(EQ_F) - m0, 'recycled': recycled(EQ_F) - r0}
        for i in objs:                                                         # the original drops extra migrants without recycling them: release them here (harness cleanup)
            if i not in after and objs[i].is_alive: objs[i].recycle()
        check_world(world, fs, '%s %s seed=%d' % (key, variant, s))
        teardown(world, fs)
        recs.append(rec)
    return recs


def section_equivalence():
    out = {}
    jobs = [(k, v) for k in EQUIV for v in ('new', 'orig')]
    with ThreadPoolExecutor(max_workers=6) as ex:
        res = dict(zip(jobs, ex.map(lambda kv: _spawn(['nemo', '--equiv', kv[0], kv[1]]), jobs)))
    bad = []
    for key in EQUIV:
        a, b = res[(key, 'new')], res[(key, 'orig')]
        if isinstance(a, dict) or isinstance(b, dict):
            bad.append((key, 'worker failed', a if isinstance(a, dict) else b)); continue
        n = len(a)
        st = {'seeds': n, 'state': sum(x['digest'] == y['digest'] for x, y in zip(a, b)), 'rng': sum(x['rng'] == y['rng'] for x, y in zip(a, b)),
              'return': sum(x['ret'] == y['ret'] for x, y in zip(a, b)), 'left_same': sum(x['left'] == y['left'] for x, y in zip(a, b)),
              'new_left_all_recycled': sum(not x['left_alive'] for x in a), 'orig_left_alive': sum(len(y['left_alive']) for y in b),
              'births': sum(x['minted'] for x in a), 'recycled_new': sum(x['recycled'] for x in a), 'recycled_orig_in_call': sum(y['recycled'] for y in b)}
        if key == FORCED:
            st['recycled == dropped by v3'] = sum(x['left'] == y['left_alive'] and len(x['left']) == len(y['left_alive']) for x, y in zip(a, b))
            st['extras'] = sum(len(y['left_alive']) for y in b); st['extras_natural'] = sum(len(y['left_alive']) for y in b if y['seed'] % 3 != 2)
            ok = all(st[k] == n for k in ('state', 'rng', 'return', 'left_same', 'new_left_all_recycled', 'recycled == dropped by v3')) and st['extras'] > 0
        else:
            ok = all(st[k] == n for k in ('state', 'rng', 'return', 'left_same', 'new_left_all_recycled')) and st['orig_left_alive'] == 0
        st['ok'] = ok
        out[key] = st
        if not ok:
            first = next((x['seed'] for x, y in zip(a, b) if (x['digest'], x['rng'], x['ret'], x['left']) != (y['digest'], y['rng'], y['ret'], y['left'])), None)
            bad.append((key, st, 'first differing seed %s' % first))
        print('1 %-26s %-80s seeds %3d  state=%d rng=%d return=%d left=%d %s' % (key, EQUIV[key][0], n, st['state'], st['rng'], st['return'], st['left_same'], 'OK' if ok else 'DIFF'))
    assert not bad, bad
    e = out[FORCED]
    finding('info', 'dispersal_among_patches_from_offspring_pool_and_dormancy_pool', 'migrants left without an empty site are now recycled (v3 removed them from their source pool and dropped them, alive and owned by nobody); otherwise seed-identical to v3.',
            evidence='%d seeds: state / random / np.random / return identical; recycled set == the migrants v3 dropped on every seed (%d extras, %d of them without forcing the number matrix)' % (e['seeds'], e['extras'], e['extras_natural']))
    print('1 equivalence: %d functions / modes, every seed identical to %s (state digest, random + np.random state, return value, who left every owner); %s: recycled == dropped by v3 (%d extras)' % (len(EQUIV), C.ORIGINAL_COMMIT, FORCED, e['extras']))
    return out


# ================================================================ 2. colonization (statistical)
def posmap(ml):
    return {id(x): (p.name, h.name, i, j) for p in all_patches(ml) for h in p.set.values() for (i, j) in h.occupied_site_pos_ls for x in [h.set['microsite_individuals'][i][j]]}


class AddLog():
    ''' records island add_individual calls in order (instance attribute shadows the class method; nothing in metaibm/ is edited) '''
    def __init__(self, isl):
        self.log = []
        for p in all_patches(isl):
            for h in p.set.values():
                h.add_individual = self._wrap(p, h)

    def _wrap(self, p, h):
        orig = h.add_individual
        def add(indi_object, len_id, wid_id):
            self.log.append((p.name, h.name, (len_id, wid_id), (len_id, wid_id) in h.empty_site_pos_ls, indi_object))
            return orig(indi_object, len_id, wid_id)
        return add


def colonize_once(ml, isl, n, pairwise, pm, addlog):
    ''' one call; checks mainland loses exactly what islands gain and the np.random use (one uniform); then puts every colonist back to its mainland site '''
    before_ml = {id(x): x for x in held(ml)}
    del addlog.log[:]
    st = np.random.get_state()
    (isl.pairwise_sexual_colonization_from_prpagules_rains if pairwise else isl.meta_colonize_from_propagules_rains)(ml, n)
    after = np.random.get_state(); np.random.set_state(st); np.random.uniform(0, 1, 1)
    assert np.array_equal(np.random.get_state()[1], after[1]) and np.random.get_state()[2] == after[2], 'colonization used np.random other than one uniform draw'
    np.random.set_state(after)
    lost = set(before_ml) - held_ids(ml)
    gained = [e for e in addlog.log]
    assert {id(e[4]) for e in gained} == lost and len(gained) == len(lost), 'mainland loses != islands gain'
    assert all(e[3] for e in gained), 'colonist placed on an occupied site'
    for pn, hn, (i, j), _, x in gained:                                           # undo: back to the mainland site it came from
        isl.set[pn].set[hn].del_individual(i, j)
        a, b, r, c = pm[id(x)]; ml.set[a].set[b].add_individual(x, r, c)
    return gained


def section_colonization():
    out = {}
    seed_all(11); F = species(11); fs = list(F.values())
    np.random.seed(11); ml = make_mainland(); ml.meta_initialize(F, TABLE)                       # 2 habitats x 30x30 = 1800
    np.random.seed(12); isl = make_islands(3, 2, L=6, K=0, gidx0=100)                            # 216 empty sites
    pm, al = posmap(ml), AddLog(isl)
    rain = rain_for(ml)
    # (a) number colonizing = floor(n) + Bernoulli(frac)
    cnt = {}
    for n in (rain, 6.4, 0.3):
        N, extra, counts = 1500, 0, {}
        for s in range(N):
            seed_all(20000 + s)
            g = colonize_once(ml, isl, n, False, pm, al)
            assert len(g) in (int(n), int(n) + 1), (n, len(g))
            extra += len(g) == int(n) + 1; counts[len(g)] = counts.get(len(g), 0) + 1
        frac = n - int(n)
        ok = binom_ok(extra, N, frac)
        cnt['n=%g' % n] = {'seeds': N, 'counts': counts, 'P(extra)': round(extra / N, 4), 'frac': round(frac, 4), 'mean': round(sum(k * v for k, v in counts.items()) / N, 4), 'p_binom': binom_p(extra, N, frac)}
        assert ok, cnt['n=%g' % n]
    check_world([ml, isl], fs, 'count')
    out['count'] = cnt
    print('2a count: floor(n) + Bernoulli(frac) over 1500 seeds each: %s' % {k: (v['P(extra)'], v['frac'], v['p_binom']) for k, v in cnt.items()})
    # (b) each mainland individual equally likely (small mainland so every individual is hit often)
    np.random.seed(13); sm = make_mainland(L=6); sm.meta_initialize(F, TABLE)                  # 72 individuals
    pm2 = posmap(sm); hits = {i: 0 for i in pm2}
    for s in range(2000):
        seed_all(30000 + s)
        for e in colonize_once(sm, isl, 6.4, False, pm2, al): hits[id(e[4])] += 1
    tot = sum(hits.values())
    ok_b, stat_b, p_b = chi2_test(list(hits.values()), [tot / len(hits)] * len(hits))
    out['mainland_uniform'] = {'seeds': 2000, 'individuals': len(hits), 'picks': tot, 'chi2': stat_b, 'df': len(hits) - 1, 'p': p_b}
    assert ok_b, out['mainland_uniform']
    print('2b mainland pick uniform: %d picks over %d individuals, chi2 = %.1f (df %d), p = %.3f' % (tot, len(hits), stat_b, len(hits) - 1, p_b))
    # (c) each empty island site equally likely (islands partly occupied)
    isl.meta_initialize(F, TABLE)
    rng = random.Random(14)
    for h in all_habs(isl):
        for (i, j) in rng.sample(list(h.occupied_site_pos_ls), 4):
            x = h.set['microsite_individuals'][i][j]; h.del_individual(i, j); x.recycle()
    empt = [(p.name, h.name, s) for p in all_patches(isl) for h in p.set.values() for s in h.empty_site_pos_ls]
    site_hits = {e: 0 for e in empt}
    for s in range(2000):
        seed_all(40000 + s)
        for pn, hn, ij, _, _ in colonize_once(ml, isl, 6.4, False, pm, al): site_hits[(pn, hn, ij)] += 1
    tot = sum(site_hits.values())
    ok_c, stat_c, p_c = chi2_test(list(site_hits.values()), [tot / len(empt)] * len(empt))
    out['island_site_uniform'] = {'seeds': 2000, 'empty_sites': len(empt), 'placements': tot, 'chi2': stat_c, 'df': len(empt) - 1, 'p': p_c}
    assert ok_c, out['island_site_uniform']
    print('2c empty island site uniform: %d placements over %d empty sites, chi2 = %.1f (df %d), p = %.3f' % (tot, len(empt), stat_c, len(empt) - 1, p_c))
    # (d) pairwise
    isl.meta_clear_up_microsite_individuals()
    pairs_n, extra, homes, by_mh = {}, 0, set(), {}
    N = 1500
    for s in range(N):
        seed_all(50000 + s)
        g = colonize_once(ml, isl, 6.4, True, pm, al)
        assert len(g) % 2 == 0 and len(g) // 2 in (3, 4), len(g)
        for (pa, ha, sa, _, f), (pb, hb, sb, _, m) in zip(g[0::2], g[1::2]):
            assert f.gender == 'female' and m.gender == 'male' and f.species_id == m.species_id, 'pair not (female, male) of one species'
            assert pm[id(f)][:2] == pm[id(m)][:2], 'pair from two mainland habitats'
            assert (pa, ha) == (pb, hb) and sa != sb, 'pair not in two empty sites of one island habitat'
            by_mh[pm[id(f)][1]] = by_mh.get(pm[id(f)][1], 0) + 1
            homes.add(f.nemo_indi_obj.getHome() == ml.set[pm[id(f)][0]].set[pm[id(f)][1]].global_index)
        extra += len(g) == 8; pairs_n[len(g) // 2] = pairs_n.get(len(g) // 2, 0) + 1
    ok_d = binom_ok(extra, N, 0.2)
    out['pairwise'] = {'seeds': N, 'pairs': pairs_n, 'P(4 pairs)': round(extra / N, 4), 'expected': 0.2, 'p_binom': binom_p(extra, N, 0.2), 'pairs_by_mainland_habitat': by_mh, 'colonists_keep_mainland_home': homes == {True}}
    assert ok_d and homes == {True}, out['pairwise']
    print('2d pairwise over %d seeds: pairs %s (P(4) = %.3f, expected 0.2, p = %.3f), every pair (female, male) of one species from one mainland habitat into two empty sites of one island habitat; '
          'colonists keep their mainland homepatch' % (N, pairs_n, extra / N, binom_p(extra, N, 0.2)))
    # (e) capacity-limited: fewer empty sites than propagules
    isl.meta_initialize(F, TABLE)
    hs = all_habs(isl)
    for h in hs[:2]:
        x = h.set['microsite_individuals'][0][0]; h.del_individual(0, 0); x.recycle()
    lim = {}
    for pw in (False, True):
        sizes = set()
        for s in range(200):
            seed_all(60000 + s); sizes.add(len(colonize_once(ml, isl, 6.4, pw, pm, al)))
        lim['pairwise' if pw else 'single'] = sorted(sizes)
    assert lim == {'single': [2], 'pairwise': [0]}, lim
    out['capacity_limited'] = lim
    print('2e capacity: 2 empty sites in two habitats -> single colonization places exactly 2, pairwise places 0 (no habitat with 2 empty sites); mainland loses only what lands')
    check_world([ml, isl, sm], fs, 'colonization end')
    teardown([ml, isl, sm], fs)
    finding('info', 'meta_colonize_from_propagules_rains, pairwise_sexual_colonization_from_prpagules_rains', 'colonization now MOVES the propagules (is_remove=True default: del_individual on the mainland, no recycle) -- v3 default kept them on the mainland too (two owners); '
            'the stochastic rounding is int(n + U) and the mainland shuffle is gone, so runs are not seed-identical to v3; tested statistically instead.',
            evidence='P(extra): %s; mainland pick chi2 p = %.3f; island site chi2 p = %.3f; pairwise P(4 pairs) = %.3f (p = %.3f); np.random use = exactly one uniform per call' % (
                {k: (v['P(extra)'], v['frac']) for k, v in cnt.items()}, p_b, p_c, extra / N, binom_p(extra, N, 0.2)))
    return out


# ================================================================ 3. Nemo level
WATCH, VIOL, RECYCLED_REFS = [], [], []


def instrument_recycle():
    raw = individual.recycle

    @functools.wraps(raw)
    def rec(self):
        for m in WATCH:
            for h in all_habs(m):
                if any(y is self for row in h.set['microsite_individuals'] for y in row) or any(y is self for q in POOLS for y in getattr(h, q)):
                    VIOL.append('recycle() while still owned by %s.%s' % (m.metacommunity_name, h.name))
        r = raw(self)
        RECYCLED_REFS.append(self)
        return r
    individual.recycle = rec


def section_nemo_level():
    out = {}
    instrument_recycle()
    seed_all(71); F = species(71); fs = list(F.values())
    np.random.seed(71); ml = make_mainland(L=12); isl = make_islands(3, 2, L=5, K=15, gidx0=200)
    WATCH[:] = [ml, isl]
    ml.meta_initialize(F, TABLE); isl.meta_initialize(F, TABLE)
    rain = rain_for(ml, 0.01)
    home = {'offspring_pool': 0, 'mainland_birth': 0, 'marker_birth': 0, 'marker_birth_cross_habitat': 0, 'marker_birth_cross_patch': 0}
    patch_of = {id(h): p.name for p in all_patches(isl) for h in p.set.values()}
    for t in range(40):
        r0 = recycled(F)
        ml.meta_dead_selection(0.2, 0.5, 'niche_gaussian'); assert all(not x.is_alive for x in RECYCLED_REFS)
        res = {x.individual_id: h for h in all_habs(ml) for x in micro(h)}; ids0 = set(map(id, held(ml))); m0 = minted(F)
        ml.meta_mainland_sexual_birth_mutate_germinate(1.0, None, None)
        new = [(h, x) for h in all_habs(ml) for x in micro(h) if id(x) not in ids0]
        assert len(new) == minted(F) - m0 and all(res[x.mother_id] is h and res[x.father_id] is h and x.nemo_indi_obj.getHome() == h.global_index for h, x in new), 'mainland birth homepatch / parents'
        home['mainland_birth'] += len(new)
        getattr(isl, ('meta_colonize_from_propagules_rains', 'pairwise_sexual_colonization_from_prpagules_rains')[t % 2])(ml, rain)
        isl.meta_dead_selection(0.2, 0.5, 'niche_gaussian')
        res = {x.individual_id: h for h in all_habs(isl) for x in micro(h)}
        if t % 2 == 0:
            m0 = minted(F)
            isl.meta_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
            kids = [(h, x) for h in all_habs(isl) for x in h.offspring_pool]
            assert len(kids) == minted(F) - m0 and all(res[x.mother_id] is h and res[x.father_id] is h and x.nemo_indi_obj.getHome() == h.global_index for h, x in kids), 'offspring pool homepatch'
            home['offspring_pool'] += len(kids)
            v = C.POOL_VARIANTS[list(C.POOL_VARIANTS)[t // 2 % 3]]
            getattr(isl, v[0])(0.4); getattr(isl, v[1])(0.4, method='exponential', rho=2); getattr(isl, v[2])(); getattr(isl, v[3])()
        else:
            isl.meta_sex_reproduce_calculation_with_offspring_marker_pool(1.0)
            isl.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(0.5)
            isl.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool(0.5, method='uniform')
            hab_of = {(p.name, h.name): h for p in all_patches(isl) for h in p.set.values()}
            before_mk = {id(h): [mk for q in MARKER_POOLS for mk in getattr(h, q)] for h in all_habs(isl)}
            ids0, m0 = set(map(id, held(isl))), minted(F)
            isl.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None)
            consumed = {}
            for h in all_habs(isl):
                left = [mk for q in MARKER_POOLS for mk in getattr(h, q)]
                for mk in set(before_mk[id(h)]):
                    k = before_mk[id(h)].count(mk) - left.count(mk)
                    if k: consumed[(id(h), id(hab_of[mk[:2]]))] = consumed.get((id(h), id(hab_of[mk[:2]])), 0) + k
            born, res_after = {}, {x.individual_id: (h, x) for h in all_habs(isl) for x in micro(h)}   # parents never move during the call; a newborn may itself become a parent
            for h in all_habs(isl):
                for x in micro(h):
                    if id(x) in ids0: continue
                    (B, mother), (B2, _) = res_after[x.mother_id], res_after[x.father_id]
                    assert B2 is B and x.nemo_indi_obj.getHome() == B.global_index and x.indfactory is mother.indfactory, 'marker birth: home != birth habitat'
                    born[(id(h), id(B))] = born.get((id(h), id(B)), 0) + 1
                    home['marker_birth'] += 1; home['marker_birth_cross_habitat'] += B is not h
                    home['marker_birth_cross_patch'] += patch_of[id(B)] != patch_of[id(h)]
            assert born == consumed and sum(born.values()) == minted(F) - m0, ('newborns per (germination habitat, birth habitat) != markers consumed', born, consumed)
            isl.meta_clear_up_offspring_marker_and_immigrant_marker_pool()
        isl.meta_disturbance_process_in_habitat(0.03)
        isl.meta_clear_up_offspring_and_immigrant_pool()
        check_world([ml, isl], fs, 'nemo_level t=%d' % t)
        assert minted(F) - recycled(F) == len(held([ml, isl]))
    assert home['marker_birth_cross_habitat'] > 0 and home['marker_birth_cross_patch'] > 0, home
    assert not VIOL, VIOL[:5]
    for x in RECYCLED_REFS[:300]:
        assert not x.is_alive and x.nemo_indi_obj is None
        for f in (lambda: x.gender, lambda: x.individual_id, lambda: x.get_indi_phenotype_ls(), lambda: x.recycle()):
            try: f(); raise AssertionError('recycled individual still readable')
            except RuntimeError: pass
    n_rec = len(RECYCLED_REFS)
    out['home'] = home; out['recycle_calls'] = n_rec
    out['factories'] = {f.species_id: (f.minted_num, f.recycled_num) for f in fs}
    teardown([ml, isl], fs)
    assert not VIOL
    out['teardown'] = {f.species_id: (f.minted_num, f.recycled_num) for f in fs}
    print('3 nemo level: homepatch == birth habitat global_index for %(offspring_pool)d offspring-pool births, %(mainland_birth)d mainland births, %(marker_birth)d marker births '
          '(%(marker_birth_cross_habitat)d born in another habitat than they germinate in, %(marker_birth_cross_patch)d in another patch): home = marker habitat, newborns per (germination, birth habitat) == markers consumed' % home)
    print('3 recycle: %d recycle() calls, each on an individual owned by nobody; recycled individuals raise on every read; minted - recycled == held after every step; teardown %s' % (n_rec, out['teardown']))
    return out


# ================================================================ 4. edges (pre-existing behaviour, recorded)
def section_edges():
    out = {}
    seed_all(81); F = species(81); fs = list(F.values())
    np.random.seed(81); ml = make_mainland(L=3); isl = make_islands(3, 2, L=4, K=5, gidx0=100)
    ml.meta_initialize(F, TABLE); isl.meta_initialize(F, TABLE)
    isl.meta_dead_selection(0.5, 0.5, 'niche_gaussian')
    res = {}
    for label, call in (('rain 30 > 18 mainland individuals', lambda: isl.meta_colonize_from_propagules_rains(ml, 30)),
                        ('pairwise rain 30 > 2 x mainland pairs', lambda: isl.pairwise_sexual_colonization_from_prpagules_rains(ml, 30))):
        d0, m0 = digest([ml, isl]), minted(F)
        try: call(); r = 'no error'
        except Exception as e: r = '%s: %s' % (type(e).__name__, str(e)[:50])
        res[label] = (r, 'state unchanged' if digest([ml, isl]) == d0 else 'STATE CHANGED')
        check_world([ml, isl], fs, label)
    ml.meta_clear_up_microsite_individuals()
    errs = 0
    for s in range(200):
        seed_all(s)
        try: isl.meta_colonize_from_propagules_rains(ml, 0.4)
        except ValueError: errs += 1
    res['empty mainland, rain 0.4 (200 seeds)'] = ('ValueError in %d / 200 calls (expected ~80)' % errs, 'state unchanged')
    assert 50 < errs < 110
    out['colonization'] = res
    finding('known', 'meta_colonize_from_propagules_rains, pairwise_sexual_colonization_from_prpagules_rains', 'random.sample raises ValueError when int(n + U) exceeds the mainland population (or its female-male pairs), e.g. an empty mainland with rain 0.4 fails in ~40% of calls; nothing is moved (pre-existing, v3 identical).',
            evidence='; '.join('%s -> %s, %s' % (k, v[0], v[1]) for k, v in res.items()), fix='sample min(k, len(population)) like sample_offspring_without_replacement does')
    # stale markers: a marker whose birth habitat lost its parents (disturbance between marker dispersal and germinate-and-birth)
    ml.meta_initialize(F, TABLE); isl.meta_initialize(F, TABLE)
    stale = {'ValueError': 0, 'ok': 0}
    for s in range(30):
        seed_all(500 + s); isl.meta_initialize(F, TABLE)
        isl.meta_dead_selection(0.3, 0.5, 'niche_gaussian'); isl.meta_sex_reproduce_calculation_with_offspring_marker_pool(1.0)
        isl.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool(0.8, method='uniform')
        isl.set['p0'].patch_disturbance_process()                                 # p0 loses its parents; p0's markers sitting in p1 / p2 survive
        try: isl.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None); stale['ok'] += 1
        except ValueError: stale['ValueError'] += 1
        check_world([ml, isl], fs, 'stale markers %d' % s)
        isl.meta_clear_up_offspring_marker_and_immigrant_marker_pool()
    assert stale['ValueError'] > 0, stale
    out['stale_markers'] = stale
    finding('risk', 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool', 'a marker whose birth habitat has no female-male pair any more raises ValueError (random.sample of an empty pair list in hab_sex_reproduce_mutate_with_num) part-way through the loop: earlier markers are already born. '
            'Not reachable in the experiments/model-sloss*.py order (disturbance after germination), reachable if a disturbance / clear-up runs between marker dispersal and germination; pre-existing in v3. Invariants hold afterwards.',
            evidence='marker calc -> among-patch marker dispersal (m=0.8) -> disturbance of p0 -> germinate-and-birth: ValueError in %d / 30 seeds' % stale['ValueError'],
            fix='skip (or drop) a marker whose birth habitat has no parent pair, or clear the markers of a disturbed habitat everywhere')
    teardown([ml, isl], fs)
    print('4 edges: colonization %s; stale markers after a disturbance -> %s' % (res, stale))
    return out


# ================================================================ 5. probes
ASEX_MIX = [('meta_mainland_asexual_birth_mutate_germinate', lambda ml, isl: ml.meta_mainland_asexual_birth_mutate_germinate(0.5, 0.01, [0.01, 0.01])),
            ('meta_mainland_mixed_birth_mutate_germinate', lambda ml, isl: ml.meta_mainland_mixed_birth_mutate_germinate(0.5, 0.5, 0.01, [0.01, 0.01])),
            ('meta_asex_reproduce_mutate_into_offspring_pool', lambda ml, isl: isl.meta_asex_reproduce_mutate_into_offspring_pool(0.5, 0.01, [0.01, 0.01])),
            ('meta_mix_reproduce_mutate_into_offspring_pool', lambda ml, isl: isl.meta_mix_reproduce_mutate_into_offspring_pool(0.5, 0.5, 0.01, [0.01, 0.01])),
            ('meta_asex_reproduce_calculation_into_offspring_marker_pool', lambda ml, isl: isl.meta_asex_reproduce_calculation_into_offspring_marker_pool(0.5)),
            ('meta_mix_reproduce_calculation_with_offspring_marker_pool', lambda ml, isl: isl.meta_mix_reproduce_calculation_with_offspring_marker_pool(0.5, 0.5))]


def probe(i):
    seed_all(3); F = species(3)
    np.random.seed(3); ml = make_mainland(L=6); isl = make_islands(2, 2, L=5, K=5, gidx0=100)
    ml.meta_initialize(F, TABLE); isl.meta_initialize(F, TABLE)
    for m in (ml, isl):
        for _ in range(2): m.meta_dead_selection(0.0, 5.0, 'niche_gaussian')
        m.meta_dead_selection(0.0, 0.3, 'niche_gaussian')                    # empties + both asexual parents (fit) and sexual parents (unfit) for the mix paths
    name, fn = ASEX_MIX[i]
    try: r = fn(ml, isl); res = 'returned %r' % (r.strip()[:60],)
    except Exception as e: res = '%s: %s' % (type(e).__name__, str(e)[:70])
    return [name, res]


def section_probes():
    out = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        probes = list(ex.map(lambda i: _spawn(['nemo', '--probe', str(i)]), range(len(ASEX_MIX))))
    probes = [p if isinstance(p, list) else [ASEX_MIX[i][0], 'process died: %s' % p.get('error', '')[:80]] for i, p in enumerate(probes)]
    out['asex_mix'] = probes
    for n, r in probes:
        if 'calculation' not in n: assert r.startswith(('TypeError', 'NotImplementedError')), (n, r)
    finding('known', ', '.join(n for n, _ in probes), 'asex / mix meta functions are not migrated: ' + '; '.join('%s -> %s' % (n, r.split(':')[0] if not r.startswith('returned') else 'returns (markers only, no individual)') for n, r in probes) + ' (TypeError = copy.deepcopy of a SwigPyObject in the habitat layer).')
    print('5.1 asex / mix probes:', ['%s: %s' % (n, r.split(':')[0] if not r.startswith('returned') else 'returned') for n, r in probes])
    # asexual / mix markers in germinate-and-birth: NotImplementedError inside the branch, earlier markers may already be born
    seed_all(4); F = species(4); fs = list(F.values())
    part = {}
    for mode in ('asexual', 'mix_asexual', 'mix_sexual'):
        part[mode] = {'raised': 0, 'born_before_raise': 0, 'seeds_with_partial_births': 0}
        for s in range(20):
            seed_all(700 + s); isl = make_islands(2, 2, L=4, K=0, gidx0=100); isl.meta_initialize(F, TABLE)
            isl.meta_dead_selection(0.4, 0.5, 'niche_gaussian'); isl.meta_sex_reproduce_calculation_with_offspring_marker_pool(1.0)
            h = all_habs(isl)[0]; h.immigrant_marker_pool.append(('p1', 'h0', mode))
            n0 = indi_n(isl)
            try: isl.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None); raised = False
            except NotImplementedError: raised = True
            b = indi_n(isl) - n0
            if raised: part[mode]['raised'] += 1; part[mode]['born_before_raise'] += b; part[mode]['seeds_with_partial_births'] += b > 0
            check_world(isl, fs, 'partial %s' % mode); teardown(isl, fs)
        assert part[mode]['raised'] > 0 and part[mode]['seeds_with_partial_births'] > 0, part
    out['marker_birth_partial'] = part
    finding('info', 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool', 'asexual / mix_asexual / mix_sexual markers raise NotImplementedError inside the reproduce_mode branch (TODO(makeClone) / TODO(mixed)), after earlier markers of the same call were already born (partial state by design); bookkeeping stays consistent.',
            evidence='; '.join('%s: raised in %d/20 seeds (the marker is not always reached: fewer empty sites than markers), %d births before the raise' % (k, v['raised'], v['born_before_raise']) for k, v in part.items()))
    print('5.2 asexual / mix markers in germinate-and-birth:', part)
    # extension: offspring network dispersal samples without removal -> two owners (KNOWN, demonstrated, does not fail the suite)
    seed_all(5); np.random.seed(5); isl = make_islands(3, 2, L=4, K=0, gidx0=100); isl.meta_initialize(F, TABLE)
    isl.meta_dead_selection(0.3, 0.5, 'niche_gaussian'); isl.meta_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
    isl.dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool(0.6, method='exponential', rho=2)
    hs = all_habs(isl)
    off = {id(x) for h in hs for x in h.offspring_pool}
    twice = [x for h in hs for x in h.immigrant_pool if id(x) in off]
    try: check_world(isl, fs, 'extension'); inv = 'invariants hold'
    except AssertionError as e: inv = 'invariant broken: %s' % str(e)[:60]
    for h in hs: h.immigrant_pool = [x for x in h.immigrant_pool if id(x) not in off]   # harness cleanup only: drop the second owner before teardown (double recycle otherwise)
    teardown(isl, fs)
    isl2 = make_islands(3, 2, L=4, K=0, gidx0=100); isl2.meta_initialize(F, TABLE); isl2.meta_sex_reproduce_calculation_with_offspring_marker_pool(1.0)
    o0 = pool_n(isl2, ('offspring_marker_pool',))
    isl2.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool(0.6, method='exponential', rho=2)
    mk = (o0, pool_n(isl2, ('offspring_marker_pool',)), pool_n(isl2, ('immigrant_marker_pool',)))
    teardown(isl2, fs)
    out['extension'] = {'offspring_network_two_owner_individuals': len(twice), 'check_world': inv, 'marker_network (source before, source after, immigrant markers)': mk}
    if not twice: finding('info', 'extension/global_habitat_network.py: dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool', 'no two-owner individual observed any more (was KNOWN): re-check the probe.', evidence=inv)
    else: finding('known','extension/global_habitat_network.py: dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool', 'not migrated: random.sample(h_i.offspring_pool, n) without removal + immigrant_pool += -> the SAME individual is in the source offspring_pool and a target immigrant_pool (two owners); a later clear-up would recycle it twice.',
            evidence='%d individuals with two owners after one call; check_world: %s' % (len(twice), inv), fix='use patch / habitat sample_offspring_without_replacement (h_i_object.sample_offspring_without_replacement(n)) as the migrated metacommunity functions do')
    finding('known', 'extension/global_habitat_network.py: dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool', 'markers are copied, not moved (v3 is_remove=False semantics; tuples, no ownership problem): offspring markers %d -> %d, immigrant markers +%d.' % mk)
    print('5.3 extension offspring network: %d individuals end up with two owners (%s); marker network copies markers %s' % (len(twice), inv, mk))
    return out


# ================================================================ entry points
def run_section(name):
    try:
        res = globals()['section_' + name](); ok = True
    except Exception:
        res, ok = {'error': traceback.format_exc()}, False
        print('---- section %s FAILED\n%s' % (name, res['error']))
    save('nemo_' + name, {'ok': ok, 'result': res})
    return ok


def nemo_cli():
        if len(sys.argv) > 3 and sys.argv[1] == '--equiv':
            print('JSON ' + json.dumps(equiv_worker(sys.argv[2], sys.argv[3]))); sys.exit(0)
        if len(sys.argv) > 2 and sys.argv[1] == '--probe':
            print('JSON ' + json.dumps(probe(int(sys.argv[2])))); sys.exit(0)
        if len(sys.argv) > 2 and sys.argv[1] == '--section':
            sys.exit(0 if run_section(sys.argv[2]) else 1)
        wanted = sys.argv[1:] or list(SECTIONS)
        failed, results, findings = [], {}, []
        for name in wanted:
            print('==== section %s' % name, flush=True)
            pr = subprocess.run([sys.executable, os.path.abspath(__file__), 'nemo', '--section', name], capture_output=True, text=True, cwd=HERE)
            for l in (pr.stdout + pr.stderr).splitlines():
                if l.strip() and not l.startswith(NOISE): print(l)
            try: data = json.load(open(os.path.join(HERE, 'results', 'nemo_%s.json' % name)))
            except Exception: data = {'ok': False, 'result': {'error': 'no result file (rc=%d)' % pr.returncode}}
            results[name] = data.get('result'); findings += data.get('findings', [])
            ok = pr.returncode == 0 and data.get('ok')
            if not ok: failed.append(name)
            print('---- section %s %s' % (name, 'PASSED' if ok else 'FAILED'), flush=True)
        C.FINDINGS[:] = findings
        results['failed_sections'] = failed
        save('nemo', results)
        if failed: sys.exit('FAILED sections: %s' % failed)
        print('PASSED')


# ############################################################################
# chain
# ############################################################################
CONFIGS = ('pool+dormancy', 'pool+clear-ups', 'marker')
SEEDS = (1, 2, 3)
CHAIN_STEPS = 200


# ================================================================ long runs
def long_run(cfg, seed, steps=CHAIN_STEPS):
    seed_all(seed); F = species(10 * seed); fs = list(F.values())
    np.random.seed(seed); ml = make_mainland(); isl = make_islands(4, 2, L=6, K=10 if cfg == 'pool+dormancy' else 0, gidx0=100)
    world = [ml, isl]
    ml.meta_initialize(F, TABLE); isl.meta_initialize(F, TABLE)
    check_world(world, fs, 'init')
    rain = rain_for(ml)                                                                     # 1800 x 0.001 = 1.8
    seen = {f.species_id: set() for f in fs}
    hw = {f.species_id: 0 for f in fs}

    def scan():
        for x in held(world): seen[x.species_id].add(int(x.nemo_indi_obj.this))
        for f in fs: hw[f.species_id] = max(hw[f.species_id], f.minted_num - f.recycled_num)
    scan()
    calls, per_op, rss, alive, distinct_at = 0, {}, {}, [], {}
    t0 = time.time()
    for t in range(steps):
        if cfg == 'marker': ops = marker_schedule(ml, isl, t, rain, birth=1.0, dist=0.02, network=True)
        else: ops = pool_schedule(ml, isl, t, rain, dorm=cfg == 'pool+dormancy', dist=0.02)
        for name, fn in ops:
            fn(); calls += 1; per_op[name] = per_op.get(name, 0) + 1
            check_world(world, fs, '%s seed=%d t=%d %s' % (cfg, seed, t, name))
            assert all(len(h.dormancy_pool) <= h.dormancy_pool_max_size for h in all_habs(isl))
            scan()
        alive.append(indi_n(isl))
        if t in (20, 100, 150, steps - 1): rss[t] = round(rss_mb(), 1); distinct_at[t] = {s: (len(seen[s]), hw[s]) for s in seen}
    out = {'cfg': cfg, 'seed': seed, 'steps': steps, 'calls': calls, 'per_op': per_op, 'rss_mb': rss, 'distinct_vs_hw': distinct_at, 'seconds': round(time.time() - t0, 1),
           'islands_alive': {'min': min(alive), 'mean': round(sum(alive) / len(alive), 1), 'last': alive[-1], 'capacity': sum(h.size for h in all_habs(isl))},
           'births': minted(F), 'mainland_alive_last': indi_n(ml)}
    assert all(n <= h for n, h in distinct_at[steps - 1].values()), ('distinct C++ addresses > alive high-water mark', distinct_at)
    assert sum(distinct_at[steps - 1][s][0] - distinct_at[100][s][0] for s in seen) < 0.01 * out['births'], ('address set keeps growing', distinct_at)
    assert rss[steps - 1] - rss[20] < 40, ('RSS grows', rss)
    assert min(alive[steps // 2:]) > 0, 'islands went extinct: the run does not exercise the island methods'
    out['teardown_recycled'] = teardown(world, fs)
    out['teardown'] = {f.species_id: (f.minted_num, f.recycled_num) for f in fs}
    for f in fs: f.purge_recycling_pool()
    return out


# ================================================================ fuzz
GUARDED = [lambda ml, isl: isl.meta_colonize_from_propagules_rains(ml, 1.8, False), lambda ml, isl: isl.pairwise_sexual_colonization_from_prpagules_rains(ml, 1.8, is_remove=False),
           lambda ml, isl: isl.dispersal_among_patches_from_offspring_pool_and_dormancy_pool(0.5, is_remove=False), lambda ml, isl: isl.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(0.5, is_remove=False),
           lambda ml, isl: isl.dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool(0.5, is_remove=False),
           lambda ml, isl: isl.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(0.5, False), lambda ml, isl: isl.meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool(0.5, False),
           lambda ml, isl: isl.meta_dispersal_within_patch_from_offspring_and_dormancy_pool(0.5, False), lambda ml, isl: isl.meta_local_germinate_from_offspring_and_dormancy_pool(False),
           lambda ml, isl: isl.meta_local_germinate_from_offspring_and_immigrant_pool(False), lambda ml, isl: isl.meta_local_germinate_from_offspring_immigrant_and_dormancy_pool(False)]


def _has_pair(h):
    return any(d.get('female') and d.get('male') for d in h.species_category.values())


def _pairs(m):
    return sum(min(len(d.get('female', [])), len(d.get('male', []))) for h in all_habs(m) for d in h.species_category.values())


def fuzz(n_runs=30, n_ops=110):
    seed_all(81); F = species(81); fs = list(F.values())
    counts, allowed, n = {}, {}, 0
    for run in range(n_runs):
        rng = random.Random(run); seed_all(run)
        ml = make_mainland(n_habs=rng.choice((1, 2, 3)), L=rng.choice((3, 6)))
        isl = make_islands(rng.choice((1, 2, 3, 4)), rng.choice((1, 2, 3)), L=rng.choice((1, 2, 4)), W=rng.choice((1, 3)), K=rng.choice((0, 1, 10)), gidx0=100)
        world = [ml, isl]
        ml.meta_initialize(F, TABLE); isl.meta_initialize(F, TABLE)
        r = lambda: rng.choice((0.0, 0.2, 0.5, 0.9, 1.0))
        k = lambda: rng.choice(KERNELS)

        def among(fn, **kw):
            def call():
                kn, kk = k(); return getattr(isl, fn)(r(), method=kn, **kw, **kk)
            return call

        def colonize(pairwise):
            def call():
                rain = rng.choice((0.0, 0.4, 1.8, 6.4))
                pop = _pairs(ml) if pairwise else len(ml.meta_get_occupied_location_ls())
                need = int((rain / 2 if pairwise else rain)) + 1
                try: return (isl.pairwise_sexual_colonization_from_prpagules_rains if pairwise else isl.meta_colonize_from_propagules_rains)(ml, rain)
                except ValueError:
                    if need > pop: allowed['colonize: rain > mainland population'] = allowed.get('colonize: rain > mainland population', 0) + 1; return 'allowed ValueError'
                    raise
            return call

        def germ_birth():
            hab = {(p.name, h.name): h for p in C.all_patches(isl) for h in p.set.values()}
            stale = any(not _has_pair(hab[mk[:2]]) for h in all_habs(isl) for q in MARKER_POOLS for mk in getattr(h, q))
            try: return isl.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None, rng.random() < 0.7)
            except ValueError:
                if stale: allowed['germinate-and-birth: stale marker'] = allowed.get('germinate-and-birth: stale marker', 0) + 1; return 'allowed ValueError'
                raise

        def guard():
            snap = full_snapshot(world, fs)
            try: rng.choice(GUARDED)(ml, isl); raise AssertionError('is_remove=False accepted')
            except NotImplementedError: pass
            assert same_snapshot(snap, full_snapshot(world, fs)), 'guard changed state'

        def getters():
            assert isl.get_meta_individual_num() == indi_n(isl) and isl.show_meta_empty_sites_num() == sum(len(h.empty_site_pos_ls) for h in all_habs(isl))
            assert sorted(isl.get_meta_empty_sites_ls() + isl.meta_get_occupied_location_ls()) == sorted((p.name, h.name, i, j) for p in C.all_patches(isl) for h in p.set.values() for i in range(h.length) for j in range(h.width))
            assert isl.meta_offspring_pool_individual_num() == pool_n(isl, ('offspring_pool',)) and isl.meta_immigrant_pool_individual_num() == pool_n(isl, ('immigrant_pool',))
            assert isl.meta_dormancy_pool_individual_num() == pool_n(isl, ('dormancy_pool',)) and isl.meta_offspring_marker_pool_marker_num() + isl.meta_immigrant_marker_pool_marker_num() == pool_n(isl, MARKER_POOLS)
            for a, b in isl.get_meta_pairwise_occupied_sites_ls():
                x, y = isl.set[a[0]].set[a[1]].set['microsite_individuals'][a[2]][a[3]], isl.set[b[0]].set[b[1]].set['microsite_individuals'][b[2]][b[3]]
                assert a[:2] == b[:2] and x.gender == 'female' and y.gender == 'male' and x.species_id == y.species_id
            for a, b in isl.get_meta_pairwise_empty_sites_ls(): assert a[:2] == b[:2] and a != b
            sp = isl.get_meta_microsites_individuals_sp_id_values(); ph = isl.get_meta_microsites_individuals_phenotype_values('pheno1')
            assert sp.shape == ph.shape == (1, sum(h.size for h in all_habs(isl))) and int(np.sum(~np.isnan(sp))) == indi_n(isl)
            if isl.patch_num > 1:
                N = isl.get_dispersal_among_num_matrix(0.3); assert (N >= 0).all() and (N == np.round(N)).all()

        OPS = {'ml_dead': lambda: ml.meta_dead_selection(rng.choice((0.0, 0.1, 0.5, 1.0)), rng.choice((0.2, 0.5, 2.0)), rng.choice(('niche_gaussian', 'neutral_uniform'))),
               'ml_birth': lambda: ml.meta_mainland_sexual_birth_mutate_germinate(rng.choice((0.3, 1.0)), None, None),
               'colonize': colonize(False), 'colonize_pairwise': colonize(True),
               'dead': lambda: isl.meta_dead_selection(rng.choice((0.0, 0.1, 0.5, 1.0)), rng.choice((0.2, 0.5, 2.0)), rng.choice(('niche_gaussian', 'neutral_uniform'))),
               'reproduce': lambda: isl.meta_sex_reproduce_mutate_into_offspring_pool(rng.choice((0.0, 0.3, 0.77, 1.0)), None, None),
               'marker_calc': lambda: isl.meta_sex_reproduce_calculation_with_offspring_marker_pool(rng.random()),
               'within_marker': lambda: isl.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(r(), rng.random() < 0.5),
               'within_O->I': lambda: isl.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(r()),
               'within_OD->I': lambda: isl.meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool(r()),
               'within_OD->sites': lambda: isl.meta_dispersal_within_patch_from_offspring_and_dormancy_pool(r()),
               'among_O->I': among('dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool'),
               'among_OD->I': among('dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool'),
               'among_OD->sites': among('dispersal_among_patches_from_offspring_pool_and_dormancy_pool'),
               'among_marker': lambda: among('dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool', is_remove=rng.random() < 0.5)(),
               'network_marker': among('dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool'),
               'germ_O+D': isl.meta_local_germinate_from_offspring_and_dormancy_pool, 'germ_O+I': isl.meta_local_germinate_from_offspring_and_immigrant_pool,
               'germ_O+I+D': isl.meta_local_germinate_from_offspring_immigrant_and_dormancy_pool, 'germ+birth': germ_birth,
               'dorm_O->D': isl.meta_dormancy_process_from_offspring_pool_to_dormancy_pool, 'dorm_O+I': isl.meta_dormancy_process_from_offspring_pool_and_immigrant_pool,
               'dist_patches': lambda: isl.meta_disturbance_process_in_patches(rng.choice((0.0, 0.3, 1.0))), 'dist_habitat': lambda: isl.meta_disturbance_process_in_habitat(rng.choice((0.0, 0.3, 1.0))),
               'clear_micro': isl.meta_clear_up_microsite_individuals, 'clear_O': isl.meta_clear_up_offspring_pool, 'clear_I': isl.meta_clear_up_immigrant_pool,
               'clear_D': isl.meta_clear_up_dormancy_pool, 'clear_OM': isl.meta_clear_up_offspring_marker_pool, 'clear_IM': isl.meta_clear_up_immigrant_marker_pool,
               'clear_O+I': isl.meta_clear_up_offspring_and_immigrant_pool, 'clear_OM+IM': isl.meta_clear_up_offspring_marker_and_immigrant_marker_pool,
               'initialize': lambda: isl.meta_initialize(F, TABLE, rng.choice(('sexual', 'asexual'))),
               'ml_initialize': lambda: ml.meta_initialize(F, TABLE, rng.choice(('sexual', 'sexual', 'asexual'))),
               'guard is_remove=False': guard, 'getters': getters}
        names = list(OPS)
        for step in range(n_ops):
            op = rng.choice(names)
            m0, r0 = minted(F), recycled(F)
            ret = OPS[op]()
            if op.startswith('clear_') and isinstance(ret, int) and 'M' not in op[6:]: assert recycled(F) - r0 == ret, (op, ret)
            check_world(world, fs, 'fuzz run=%d step=%d op=%s' % (run, step, op))
            assert all(len(h.dormancy_pool) <= h.dormancy_pool_max_size for h in all_habs(isl))
            counts[op] = counts.get(op, 0) + 1; n += 1
        teardown(world, fs)
    assert n >= 3000 and len(counts) == len(names), (n, sorted(set(names) - set(counts)))
    return {'operations': n, 'kinds': len(counts), 'per_op': counts, 'allowed_exceptions': allowed}


# ================================================================ entry points
def chain_main():
    out = {'long': {}, 'fuzz': None}
    t0 = time.time()
    jobs = [(c, s) for c in CONFIGS for s in SEEDS] + [('fuzz', 0)]
    with ThreadPoolExecutor(max_workers=6) as ex:
        res = dict(zip(jobs, ex.map(lambda j: _spawn(['chain', '--fuzz'] if j[0] == 'fuzz' else ['chain', '--run', j[0], str(j[1])]), jobs)))
    failed = []
    for (c, s), r in res.items():
        if 'error' in r: failed.append((c, s, r['error'])); continue
        if c == 'fuzz':
            out['fuzz'] = r
            print('2 fuzz: %d random operations (%d kinds) over every migrated metacommunity method + guarded is_remove=False calls + getters, invariants after each, minted == recycled after every run; '
                  'allowed pre-existing exceptions %s' % (r['operations'], r['kinds'], r['allowed_exceptions']))
            continue
        out['long']['%s seed=%d' % (c, s)] = r
        d = r['distinct_vs_hw'][str(CHAIN_STEPS - 1)]
        print('1 %-15s seed=%d: %d calls (invariants after each) in %.0fs, births %d, islands alive mean %.0f / %d (min %d), RSS %s MB, distinct addresses vs alive high-water mark %s, teardown %s' %
              (c, s, r['calls'], r['seconds'], r['births'], r['islands_alive']['mean'], r['islands_alive']['capacity'], r['islands_alive']['min'], r['rss_mb'], d, r['teardown']))
    for c, s, e in failed: print('FAILED %s seed=%s\n%s' % (c, s, e))
    # the fuzz's allowed exceptions are the pre-existing ones recorded as findings by test_nemo.py (edges): counted here, not re-reported
    out['seconds'] = round(time.time() - t0, 1)
    out['failed'] = [(c, s) for c, s, _ in failed]
    save('chain', out)
    if failed: sys.exit('FAILED: %s' % out['failed'])
    print('PASSED')


def chain_cli():
        if len(sys.argv) > 3 and sys.argv[1] == '--run':
            print('JSON ' + json.dumps(long_run(sys.argv[2], int(sys.argv[3])))); sys.exit(0)
        if len(sys.argv) > 1 and sys.argv[1] == '--fuzz':
            print('JSON ' + json.dumps(fuzz())); sys.exit(0)
        chain_main()


# ############################################################################
# all: Step 0 baseline, every part, coverage, TEST_REPORT.md
# ############################################################################
TESTS = [('bookkeeping', 'R1 is_remove defaults, R3 else branches, R2 is_remove=False guards (state + both RNGs + counters + currentID), R4 marker modes, clear-ups, meta_initialize, R6 tag audit + docstring rule, patch_num < 2 logs, unchanged helpers'),
         ('nemo', 'one fresh process per section: equivalence with the original bodies (one worker process per function and variant) / colonization statistics / Nemo level / edges / probes'),
         ('chain', 'long runs (3 schedules x 3 seeds x 200 steps, invariants after every call, RSS, C++ addresses, teardown) + fuzz (>= 3000 random operations)')]
FILES = ['metaibm/metacommunity.py', 'metaibm/patch.py', 'metaibm/habitat.py', 'metaibm/indfactory.py', 'metaibm/individual.py', 'extension/global_habitat_network.py']
RA_NOISE = ('setting random seed from input value', 'outputs: ', 'the microsite in the habitat is occupied.', 'the microsite in the habitat is empty.')
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
BASE = [None]


def ra_run(script, cwd=HERE, keep=40, args=()):
    t = time.time()
    p = subprocess.run([sys.executable, os.path.join(cwd, script)] + list(args), capture_output=True, text=True, cwd=cwd, env=ENV)
    lines = [l for l in (p.stdout + '\n' + p.stderr).splitlines() if l.strip() and not any(n in l for n in RA_NOISE) and not l.startswith('[FINDING')]
    return {'rc': p.returncode, 's': time.time() - t, 'tail': lines[-keep:]}


def load(name):
    try: return json.load(open(os.path.join(HERE, 'results', name + '.json')))
    except Exception: return {}


def git(*a):
    return subprocess.run(['git'] + list(a), capture_output=True, text=True, cwd=REPO).stdout.rstrip()


def md5(rel):
    return hashlib.md5(open(os.path.join(REPO, rel), 'rb').read()).hexdigest()


# ------------------------------------------------------------------ Step 0 baseline
def _manifest(d):
    ''' {relative path: md5} of every file in the folder (git ignores results/ and __pycache__, so git status alone cannot prove the restore) '''
    return {os.path.relpath(os.path.join(r, f), d): hashlib.md5(open(os.path.join(r, f), 'rb').read()).hexdigest() for r, _, fs in os.walk(d) for f in fs}


def _run_suite_restored(rel, args):
    ''' run <rel>/run_all.py with a backup of the whole folder, restore it afterwards (its run_all rewrites TEST_REPORT.md + results/); keep its fresh report in our results/ '''
    d = os.path.join(REPO, rel)
    before = _manifest(d)
    tmp = tempfile.mkdtemp(prefix='suite_backup_')
    shutil.copytree(d, os.path.join(tmp, 'suite'), copy_function=shutil.copy2)
    try:
        r = ra_run('run_all.py', cwd=d, keep=40, args=args)
        rep = open(os.path.join(d, 'TEST_REPORT.md')).read() if os.path.exists(os.path.join(d, 'TEST_REPORT.md')) else ''
    finally:
        shutil.rmtree(d); shutil.copytree(os.path.join(tmp, 'suite'), d, copy_function=shutil.copy2); shutil.rmtree(tmp)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    open(os.path.join(HERE, 'results', 'baseline_%s_TEST_REPORT.md' % os.path.basename(rel).replace('nemoback-', '')), 'w').write(rep)
    rows = re.findall(r'^\| (\S+\.py) \| .*? \| (PASS|\*\*FAIL\*\*) \| ([\d.]+) s \|$', rep, re.M)
    sm = re.search(r'^## Smoke: (PASS|\*\*FAIL\*\*)', rep, re.M)
    if sm: rows = [('smoke.py', sm.group(1), '')] + rows
    if not rows:                                                                          # the habitat run_all prints '<script> PASS|FAIL <s>s'
        rows = [(m.group(1), m.group(2), m.group(3)) for l in r['tail'] for m in [re.match(r'(\S+\.py)\s+(PASS|FAIL)\s+([\d.]+)s', l)] if m]
    fails = {}
    for script, res, _ in rows:
        if 'FAIL' in res:
            blk = re.search(r'%s output</summary>\n\n```\n(.*?)```' % re.escape(script), rep, re.S)
            fails[script] = (blk.group(1).strip().splitlines()[-6:] if blk else r['tail'][-6:])
    return {'rc': r['rc'], 's': round(r['s'], 1), 'rows': [(a, 'FAIL' if 'FAIL' in b else 'PASS', c) for a, b, c in rows], 'fails': fails,
            'restored_clean': git('status', '--short', '--', rel) == '' and _manifest(d) == before, 'files': len(before)}


def baseline():
    out = {'indfactory': {}}
    d = os.path.join(REPO, 'test', 'nemo-indfactory')
    for s in sorted(os.path.basename(x) for x in glob.glob(os.path.join(d, '*.py')) if not x.endswith('common.py')):
        r = ra_run(s, cwd=d, keep=3); out['indfactory'][s] = (r['rc'], round(r['s'], 1), r['tail'])
    out['habitat'] = _run_suite_restored('test/nemoback-habitat', [])
    out['patch'] = _run_suite_restored('test/nemoback-patch', ['--no-baseline'])
    return out


# ------------------------------------------------------------------ method coverage
def _strip_doc(fn):
    fn = ast.parse(ast.unparse(fn)).body[0]
    if fn.body and isinstance(fn.body[0], ast.Expr) and isinstance(getattr(fn.body[0], 'value', None), ast.Constant) and isinstance(fn.body[0].value.value, str): fn.body = fn.body[1:] or [ast.Pass()]
    return ast.dump(fn)


def coverage(orig_commit):
    src = open(os.path.join(REPO, 'metaibm', 'metacommunity.py')).read(); osrc = git('show', orig_commit + ':metaibm/metacommunity.py')
    def fns(s): return [n for n in [c for c in ast.parse(s).body if isinstance(c, ast.ClassDef) and c.name == 'metacommunity'][0].body if isinstance(n, ast.FunctionDef)]
    cur, orig = fns(src), {n.name: n for n in fns(osrc)}
    lines, olines = src.splitlines(), osrc.splitlines()
    texts = {os.path.basename(p): open(p).read() for p in glob.glob(os.path.join(HERE, '*.py')) if not p.endswith('run_all.py')}
    rows, missing = [], []
    for n in cur:
        m = n.name
        if m not in orig: status = 'new'
        elif ast.dump(orig[m]) == ast.dump(n) and '\n'.join(lines[n.lineno - 1:n.end_lineno]) == '\n'.join(olines[orig[m].lineno - 1:orig[m].end_lineno]): status = 'unchanged'
        elif _strip_doc(orig[m]) == _strip_doc(n): status = 'docstring / comment only'
        else: status = 'changed'
        tag = 'NemoBackend' if (ast.get_docstring(n) or '').startswith('NemoBackend') else ''
        scope = 'plot (unchanged, not run)' if m.startswith('meta_show_') else 'asex/mix (probe)' if re.search(r'asex|mix', m) else 'in scope'
        hits = sorted(f.replace('.py', '') for f, t in texts.items() if re.search(r'\b%s\b' % re.escape(m), t))
        if m in ('__init__', '__str__'): hits = hits or ['constructor / str()']
        rows.append((m, status, tag, scope, ', '.join(hits) or '-'))
        if scope == 'in scope' and not hits: missing.append(m)
    for m in orig:
        if m not in {n.name for n in cur}: rows.append((m, 'removed', '', 'in scope', '-'))
    return rows, missing


# ------------------------------------------------------------------ perf (record only)
def perf():
    sys.path.insert(0, HERE)
    import random, warnings
    import numpy as np
    warnings.filterwarnings('ignore')
    import common as C
    C.seed_all(91); F = C.species(91)
    out = {}
    for L in (10, 20):
        np.random.seed(L); ml = C.make_mainland(); isl = C.make_islands(4, 2, L=L, K=L * L // 2, gidx0=100)
        ml.meta_initialize(F, C.TABLE); isl.meta_initialize(F, C.TABLE)
        calls = [('meta_initialize (islands)', lambda: isl.meta_initialize(F, C.TABLE)), ('meta_dead_selection (mainland 1800)', lambda: ml.meta_dead_selection(0.1, 0.5, 'niche_gaussian')),
                 ('meta_mainland_sexual_birth_mutate_germinate', lambda: ml.meta_mainland_sexual_birth_mutate_germinate(1.0, None, None)),
                 ('meta_colonize_from_propagules_rains (6.4)', lambda: isl.meta_colonize_from_propagules_rains(ml, 6.4)),
                 ('pairwise_sexual_colonization_from_prpagules_rains (6.4)', lambda: isl.pairwise_sexual_colonization_from_prpagules_rains(ml, 6.4)),
                 ('meta_dead_selection (islands)', lambda: isl.meta_dead_selection(0.2, 0.5, 'niche_gaussian')),
                 ('meta_sex_reproduce_mutate_into_offspring_pool', lambda: isl.meta_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)),
                 ('meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool', lambda: isl.meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool(0.3)),
                 ('dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool', lambda: isl.dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool(0.2, method='exponential', rho=2)),
                 ('meta_local_germinate_from_offspring_immigrant_and_dormancy_pool', isl.meta_local_germinate_from_offspring_immigrant_and_dormancy_pool),
                 ('meta_dormancy_process_from_offspring_pool_and_immigrant_pool', isl.meta_dormancy_process_from_offspring_pool_and_immigrant_pool),
                 ('meta_dead_selection (islands) ', lambda: isl.meta_dead_selection(0.2, 0.5, 'niche_gaussian')),
                 ('meta_sex_reproduce_mutate_into_offspring_pool ', lambda: isl.meta_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)),
                 ('dispersal_among_patches_from_offspring_pool_and_dormancy_pool', lambda: isl.dispersal_among_patches_from_offspring_pool_and_dormancy_pool(0.2, method='exponential', rho=2)),
                 ('meta_sex_reproduce_calculation_with_offspring_marker_pool', lambda: isl.meta_sex_reproduce_calculation_with_offspring_marker_pool(1.0)),
                 ('dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool', lambda: isl.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool(0.2, method='exponential', rho=2)),
                 ('meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool', lambda: isl.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None)),
                 ('meta_disturbance_process_in_patches', lambda: isl.meta_disturbance_process_in_patches(0.5)),
                 ('meta_clear_up_offspring_and_immigrant_pool', isl.meta_clear_up_offspring_and_immigrant_pool),
                 ('meta_clear_up_offspring_marker_and_immigrant_marker_pool', isl.meta_clear_up_offspring_marker_and_immigrant_marker_pool)]
        tm = {}
        for rep in range(3):
            for name, fn in calls:
                t = time.perf_counter(); fn(); tm[name.strip()] = tm.get(name.strip(), 0) + (time.perf_counter() - t) / 3
            isl.meta_initialize(F, C.TABLE)
        out['mainland 2x30x30 + islands 4x2x%dx%d' % (L, L)] = {k: round(v, 4) for k, v in tm.items()}
        C.teardown([ml, isl], list(F.values()))
    return out


# ------------------------------------------------------------------ report helpers
def esc(s): return str(s if s is not None else '').replace('|', '/').replace('\n', ' ')
def table(head, rows): return ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)] + ['| ' + ' | '.join(esc(c) for c in r) + ' |' for r in rows]
def pf(ok): return 'PASS' if ok else '**FAIL**'


def section_bookkeeping(bk):
    if not bk: return ['(no results)']
    r2, r4, cu, mi, r6 = bk.get('R2', {}), bk.get('R4', {}), bk.get('clear_ups', {}), bk.get('meta_initialize', {}), bk.get('R6', {})
    rows = [('R1 every is_remove default True', 'PASS', '%d functions (original defaults: False)' % bk.get('R1', {}).get('functions', 0)),
            ('R3 guard first; unreachable else == `pass  # ... TODO(copy_from)`', 'PASS', '5 guarded meta functions, %d else branches' % sum(len(v['is_remove_branches']) for v in bk.get('R3', {}).values() if v['guard_first'])),
            ('R2 is_remove=False -> NotImplementedError before any RNG use / state change', 'PASS', '%d calls over 11 paths (3- and 1-patch islands); raised in: %s' % (r2.get('calls', 0), ', '.join(sorted({'%s layer' % v[1] for v in r2.get('paths', {}).values()})))),
            ('R4 marker paths with is_remove False and True', 'PASS', '%s, %d births' % ({k: v for k, v in r4.items() if k != 'births'}, r4.get('births', 0))),
            ('clear-ups: 6 split + 2 wrappers, int returns, dormancy kept by O+I wrapper', 'PASS', '%d cases' % cu.get('cases', 0)),
            ('meta_initialize fills every site with the right species, homepatch, re-init', 'PASS', '90 founders; bad inputs (info): %s' % {k: v['error'] for k, v in mi.get('bad_inputs', {}).items()}),
            ('R6 NemoBackend tags == functions reaching recycle / make_founder / make_offspring / tagged methods', 'PASS', '%d tagged: %s' % (len(r6.get('actual', [])), ', '.join(r6.get('actual', [])))),
            ('R6 untagged functions keep the original docstring; new ones one-line', 'PASS', '%d kept, %d new one-line (%s)' % (sum(v == 'kept' for v in r6.get('docstring_rule', {}).values()), sum(v == 'one-line' for v in r6.get('docstring_rule', {}).values()), ', '.join(k for k, v in r6.get('docstring_rule', {}).items() if v == 'one-line'))),
            ('patch_num < 2 log strings (info)', 'PASS', '4 functions now format the metacommunity name'),
            ('unchanged getters / matrix helpers / CSV writer == original bodies (same seed)', 'PASS', '%s functions x %s seeds' % (bk.get('helpers', {}).get('functions'), bk.get('helpers', {}).get('seeds')))]
    L = table(('check', 'result', 'numbers'), rows)
    why = r6.get('why', {})
    if why: L += ['', '<details><summary>R6 tag audit: why each tagged function is NemoBackend</summary>', ''] + table(('function', 'reason (calls)'), [(k, ', '.join(v)) for k, v in sorted(why.items())]) + ['', '</details>']
    return L


def section_nemo(nm):
    if not nm: return ['(no results)']
    L = []
    eq = nm.get('equivalence') or {}
    if 'error' in eq: return ['```', eq['error'][-3000:], '```']
    L += ['### Equivalence with the original bodies (one worker process per function and variant, same seed + same world)', '']
    L += table(('key', 'seeds', 'state', 'random + np.random', 'return', 'left all owners: same IDs', 'new: all left recycled', 'extra'), [
        (k, v['seeds'], v['state'], v['rng'], v['return'], v['left_same'], v['new_left_all_recycled'],
         ('recycled == dropped by v3 on %d seeds, %d extras (%d natural)' % (v.get('recycled == dropped by v3'), v.get('extras'), v.get('extras_natural')) if 'extras' in v else 'births %d, recycled %d' % (v['births'], v['recycled_new'])))
        for k, v in eq.items()])
    co = nm.get('colonization') or {}
    if co and 'error' not in co:
        L += ['', '### Colonization (statistical)', '']
        rows = [('count = floor(n) + Bernoulli(frac), %s' % k, 'PASS', '%d seeds, counts %s, mean %s, P(extra) %s vs %s, binomial p = %s' % (v['seeds'], v['counts'], v['mean'], v['P(extra)'], v['frac'], v['p_binom'])) for k, v in co['count'].items()]
        mu, iu, pw = co['mainland_uniform'], co['island_site_uniform'], co['pairwise']
        rows += [('each mainland individual equally likely', 'PASS', '%d seeds, %d picks / %d individuals, chi2 = %s (df %d), p = %s' % (mu['seeds'], mu['picks'], mu['individuals'], mu['chi2'], mu['df'], mu['p'])),
                 ('each empty island site equally likely', 'PASS', '%d seeds, %d placements / %d sites, chi2 = %s (df %d), p = %s' % (iu['seeds'], iu['placements'], iu['empty_sites'], iu['chi2'], iu['df'], iu['p'])),
                 ('pairwise: (female, male), one species, one mainland habitat -> two empty sites of one island habitat', 'PASS', '%d seeds, pairs %s, P(4 pairs) = %s vs 0.2, p = %s; pairs by mainland habitat %s' % (pw['seeds'], pw['pairs'], pw['P(4 pairs)'], pw['p_binom'], pw['pairs_by_mainland_habitat'])),
                 ('mainland loses exactly what islands gain; one np.random.uniform per call', 'PASS', 'checked on every call above'),
                 ('capacity-limited', 'PASS', str(co.get('capacity_limited')))]
        L += table(('check', 'result', 'numbers'), rows)
    nl = nm.get('nemo_level') or {}
    if nl and 'error' not in nl:
        h = nl['home']
        L += ['', '### Nemo level', ''] + table(('check', 'result', 'numbers'), [
            ('home_habitat == birth habitat global_index', 'PASS', 'offspring pool %d, mainland birth %d, marker birth %d (%d in another habitat, %d in another patch than the germination site)' % (h['offspring_pool'], h['mainland_birth'], h['marker_birth'], h['marker_birth_cross_habitat'], h['marker_birth_cross_patch'])),
            ('marker birth: newborns per (germination, birth habitat) == markers consumed', 'PASS', 'every marker-birth call'),
            ('recycle() only on an individual owned by nobody; recycled never read again', 'PASS', '%d recycle() calls' % nl['recycle_calls']),
            ('minted - recycled == held per factory; teardown minted == recycled', 'PASS', str(nl['teardown']))])
    ed = nm.get('edges') or {}
    if ed and 'error' not in ed:
        L += ['', '### Edges (pre-existing, recorded)', ''] + table(('case', 'result'), [(k, '%s, %s' % tuple(v)) for k, v in ed['colonization'].items()] + [('stale markers after a disturbance (30 seeds)', str(ed['stale_markers']))])
    pr = nm.get('probes') or {}
    if pr and 'error' not in pr:
        L += ['', '### Probes', ''] + table(('probe', 'result'), [(n, r) for n, r in pr['asex_mix']] + [('germinate-and-birth with %s marker' % k, str(v)) for k, v in pr['marker_birth_partial'].items()] +
                                             [('extension: ' + k, str(v)) for k, v in pr['extension'].items()])
    for s in ('equivalence', 'colonization', 'nemo_level', 'edges', 'probes'):
        if isinstance(nm.get(s), dict) and 'error' in nm[s]: L += ['', '**section %s failed:**' % s, '```', nm[s]['error'][-3000:], '```']
    return L


def section_chain(ch):
    if not ch or not ch.get('long'): return ['(no results)']
    rows = []
    for k, r in ch['long'].items():
        d = r['distinct_vs_hw'][str(max(map(int, r['distinct_vs_hw'])))]
        rows.append((k, r['steps'], r['calls'], '%.0f s' % r['seconds'], r['births'], '%s / %s (min %s)' % (r['islands_alive']['mean'], r['islands_alive']['capacity'], r['islands_alive']['min']),
                     ' -> '.join('%s' % v for v in r['rss_mb'].values()), ', '.join('%s %d / %d' % (s, a, b) for s, (a, b) in d.items()), ', '.join('%s %d = %d' % (s, a, b) for s, (a, b) in r['teardown'].items())))
    L = table(('run', 'steps', 'calls (invariants after each)', 'time', 'births', 'islands alive mean / capacity', 'RSS MB (t=20 -> 199)', 'distinct C++ addresses / alive high-water mark', 'teardown minted = recycled'), rows)
    fz = ch.get('fuzz') or {}
    if fz: L += ['', '- fuzz: %d random operations (%d kinds), invariants after each, minted == recycled after every run; allowed pre-existing exceptions %s' % (fz['operations'], fz['kinds'], fz['allowed_exceptions']),
                 '- fuzz operations: %s' % ', '.join('%s %d' % kv for kv in sorted(fz['per_op'].items()))]
    return L


def ra_main():
    T0 = time.time()
    for f in glob.glob(os.path.join(HERE, 'results', '*')): os.remove(f)
    sys.path.insert(0, HERE)
    try:
        import importlib.metadata as md; pv = md.version('pynemo')
    except Exception: pv = '?'
    orig_commit = subprocess.run([sys.executable, '-c', 'import common; print(common.ORIGINAL_COMMIT)'], capture_output=True, text=True, cwd=HERE, env=ENV).stdout.split()[-1]
    last_mc = git('log', '-1', '--format=%h', '--', 'metaibm/metacommunity.py')
    L = ['# metacommunity.py test report (NemoBackend, sexual)', '',
         'Generated %s by `test/nemoback-metacommunity/test_metacommunity.py all`. Scope: `metaibm/metacommunity.py` and its call chain metacommunity -> patch -> habitat -> indfactory -> individual (asex / mix and the extension offspring network only probed).' % datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), '',
         '## Environment', '', '- python `%s`, pynemo `%s`, numpy `%s`, `%s`' % (sys.version.split()[0], pv, __import__('numpy').__version__, platform.platform()),
         '- HEAD `%s` (`%s`), branch `%s`' % (git('rev-parse', '--short', 'HEAD'), git('log', '-1', '--format=%s'), git('rev-parse', '--abbrev-ref', 'HEAD')), '']
    L += table(('file', 'md5', 'state'), [('`%s`' % f, '`%s`' % md5(f), 'modified (working tree)' if git('status', '--short', '--', f) else 'committed') for f in FILES])
    L += ['', '- original reference (auto-located: newest commit whose metacommunity.py still calls `patch_clear_up_offspring_and_immigrant_pool`): `git show %s:metaibm/metacommunity.py` (blob `%s`). '
          'The last commit that changed metacommunity.py is `%s`, so `%s:metacommunity.py` == `%s:metacommunity.py`. Original bodies are compiled into a subclass of the CURRENT metacommunity '
          '(they run on the current NemoBackend patch / habitat); `meta_initialize` and the two combined clear-ups (removed patch APIs) are never compiled.' % (orig_commit, git('rev-parse', '--short', orig_commit + ':metaibm/metacommunity.py'), last_mc, orig_commit, last_mc), '']
    findings_extra = []
    # ---- Step 0
    if '--no-baseline' not in sys.argv:
        t = time.time(); b = baseline(); BASE[0] = b; print('baseline %.1fs' % (time.time() - t))
        rows = [('test/nemo-indfactory', s, pf(rc == 0), '%.1f s' % sec) for s, (rc, sec, _) in b['indfactory'].items()]
        for key, rel in (('habitat', 'test/nemoback-habitat'), ('patch', 'test/nemoback-patch')):
            rows += [('%s (run_all.py%s)' % (rel, ' --no-baseline' if key == 'patch' else ''), s, pf(res == 'PASS'), ('%s s' % sec) if sec else '') for s, res, sec in b[key]['rows']]
            rows.append(('%s' % rel, '**run_all.py total**', pf(b[key]['rc'] == 0), '%.1f s' % b[key]['s']))
        L += ['## Step 0 -- baseline (existing suites, run unchanged)', ''] + table(('suite', 'script', 'result', 'time'), rows)
        L += ['', '- both folders were backed up before and restored after their run (their run_all.py rewrites TEST_REPORT.md + results/): restored byte-identical (md5 of every file + git status) = habitat `%s`, patch `%s`; '
              'their fresh reports are kept as `results/baseline_habitat_TEST_REPORT.md` / `results/baseline_patch_TEST_REPORT.md`.' % (b['habitat']['restored_clean'], b['patch']['restored_clean'])]
        for key in ('habitat', 'patch'):
            for script, tail in b[key]['fails'].items():
                L += ['', '`test/nemoback-%s/%s` failed; last lines:' % (key, script), '```'] + tail + ['```']
        L.append('')
    # ---- smoke
    smoke = ra_run(SELF, args=['smoke']); sj = load('smoke')
    print('smoke %s %.1fs' % ('PASS' if smoke['rc'] == 0 else 'FAIL', smoke['s']))
    L += ['## Smoke: %s (%.1f s)' % (pf(smoke['rc'] == 0), smoke['s']), ''] + table(('step', 'return', 'invariants'), [(s['step'], '`%s`' % s['return'], s['invariants']) for s in sj.get('steps', [])]) + ['']
    if smoke['rc'] != 0:
        L += ['Failed step `%s`:' % sj.get('failed_step'), '```', sj.get('traceback', '\n'.join(smoke['tail'])), '```', 'Other tests not run.']
        open(os.path.join(HERE, 'TEST_REPORT.md'), 'w').write('\n'.join(L) + '\n'); print('SMOKE FAILED -- stop'); sys.exit(1)
    res = {}
    for f, _ in TESTS:
        res[f] = ra_run(SELF, keep=80, args=[f]); print('%-22s %s %.1fs' % (f, 'PASS' if res[f]['rc'] == 0 else 'FAIL', res[f]['s']))
    rows, missing = coverage(orig_commit)
    perf_res = None
    if '--perf' in sys.argv:
        t = time.time(); perf_res = perf(); print('perf (record only) %.1fs' % (time.time() - t))
    n_pass = sum(r['rc'] == 0 for r in res.values())
    L += ['## Tests: %d / %d PASS' % (n_pass, len(res)), ''] + table(('script', 'covers', 'result', 'time'), [(f, d, pf(res[f]['rc'] == 0), '%.1f s' % res[f]['s']) for f, d in TESTS]) + ['']
    bk, nm, ch = load('bookkeeping'), load('nemo'), load('chain')
    L += ['## bookkeeping: %s' % pf(res['bookkeeping']['rc'] == 0), ''] + section_bookkeeping(bk) + ['']
    L += ['## nemo: %s' % pf(res['nemo']['rc'] == 0), ''] + section_nemo(nm) + ['']
    L += ['## chain: %s' % pf(res['chain']['rc'] == 0), ''] + section_chain(ch) + ['']
    for f, _ in TESTS:
        L.append('<details><summary>%s output</summary>' % f); L += ['', '```'] + res[f]['tail'] + ['```', '</details>']
    if perf_res:
        names = list(next(iter(perf_res.values())).keys())
        L += ['', '## Performance (record only, `--perf`): seconds per call, mean of 3', ''] + table(['metacommunity call'] + list(perf_res), [['`%s`' % n] + ['%.4f' % perf_res[k][n] for k in perf_res] for n in names])
    # ---- findings
    fl = []
    if BASE[0]:
        for key in ('habitat', 'patch'):
            for script, tail in BASE[0][key]['fails'].items():
                fl.append({'severity': 'info', 'func': 'test/nemoback-%s/%s (Step 0)' % (key, script), 'text': BASE_NOTES.get((key, script), 'existing suite fails on the current tree (suite not edited).'),
                           'evidence': ' / '.join(tail[-3:]), 'fix': BASE_FIX.get((key, script))})
    for n in ('smoke', 'bookkeeping', 'nemo', 'chain'):
        fl += load(n).get('findings', [])
    seen, uniq = set(), []
    for x in fl:
        if x['text'] not in seen: seen.add(x['text']); uniq.append(x)
    order = {'bug': 0, 'risk': 1, 'info': 2, 'known': 3}
    uniq.sort(key=lambda x: order.get(x['severity'], 9))
    cnt = {s: sum(x['severity'] == s for x in uniq) for s in order}
    L += ['', '## Findings (%d: %s)' % (len(uniq), ', '.join('%s %d' % kv for kv in cnt.items())), '',
          'Severity: **bug** = wrong behaviour in the code under test; **risk** = correct today, fragile / surprising; **info** = observation or intended behaviour change vs v3; **known** = expected / next migration round / pre-existing. Suggested fixes are suggestions only -- nothing was applied.', '']
    L += table(('severity', 'function', 'finding', 'evidence', 'suggested fix (not applied)'), [(x['severity'], x.get('func'), x['text'], x.get('evidence'), x.get('fix')) for x in uniq])
    L += ['', '## Method coverage%s' % ('' if not missing else ' -- **not exercised: %s**' % missing), '',
          '`vs original` compares each method with `%s:metaibm/metacommunity.py` (AST + source text). `scripts` = test files of this suite that name the method (common = the shared schedules).' % orig_commit, '']
    L += table(('method', 'vs original', 'tag', 'scope', 'scripts'), [('`%s`' % r[0],) + r[1:] for r in rows])
    L += ['', 'Installed by `extension/global_habitat_network.py` (not part of the class body, not migrated): the marker network dispersal runs in test_chain (marker schedule, odd steps) and the fuzz; the offspring network dispersal is probed only (two owners, KNOWN).']
    L += ['', '## Repository state', '', '```', git('status', '--short'), git('diff', '--stat'), '```']
    open(os.path.join(HERE, 'TEST_REPORT.md'), 'w').write('\n'.join(L) + '\n')
    total = time.time() - T0
    print('coverage: %s' % ('all in-scope methods exercised' if not missing else 'MISSING %s' % missing))
    print('SUMMARY: tests %d/%d PASS, smoke PASS, findings %s, total %.1fs' % (n_pass, len(res), cnt, total))
    print('report:', os.path.join(HERE, 'TEST_REPORT.md'))
    sys.exit(1 if (n_pass < len(res) or missing) else 0)


BASE_NOTES = {('patch', 'test_chain.py'): 'stale since metacommunity.py was migrated (suite not edited): its section 4 static scan asserts that metacommunity.py still calls the REMOVED patch_clear_up_offspring_and_immigrant_pool / '
                                          'patch_clear_up_offspring_marker_and_immigrant_marker_pool and patch_initialize with the old 6-argument signature (the "next migration round" findings of that report); '
                                          'the migration removed exactly those calls, so `removed == set()` and the assertion fails. Sections 1-3 of that script (call chains, ordering rule, homepatch) passed before the scan. Not a code defect.'}
BASE_FIX = {('patch', 'test_chain.py'): 'in test/nemoback-patch/test_chain.py scan(): expect removed == set() and changed == set() once metacommunity.py is migrated (or compare against the git version of metacommunity.py the patch migration was tested with)'}


def ra_cli():
        ra_main()


# ############################################################################
# command line
# ############################################################################
PARTS = {'smoke': smoke_cli, 'bookkeeping': bookkeeping_main, 'nemo': nemo_cli, 'chain': chain_cli, 'all': ra_cli}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in PARTS: sys.exit('usage: python test_metacommunity.py {%s} [args]  (see the header of this file)' % '|'.join(PARTS))
    part = sys.argv.pop(1)                       # the remaining argv is what the part's own command line expects
    PARTS[part]()
