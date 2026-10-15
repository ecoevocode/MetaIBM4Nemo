# Patch-level pool / microsite bookkeeping with fake individuals (no Nemo individuals are created):
#   1 six clear-ups   2 disturbance   3 is_remove=False guards   4 equivalence with the original (git) on many seeds / shapes
#   5 forced extras   6 natural extras (binomial)   7 sampling uniformity (chi-square)   8 edge shapes (warnings / errors recorded)
import sys, math, random, warnings, traceback
import numpy as np
import common as C
warnings.filterwarnings('ignore', category=RuntimeWarning); warnings.filterwarnings('ignore', category=PendingDeprecationWarning)   # recorded explicitly in section 8
from common import (patch, make_patch, fill_fakes, habs, snapshot, held, held_labels, check_invariants, rng_state, same_rng_state, force_within_matrix,
                    original_patch, original_stack, clear_recycled, RECYCLED, POOLS, MARKER_POOLS, finding, save, chi2_ok, binom_ok)

out = {}
OP, OH = original_stack()                                  # full original stack (patch + habitat from git)
DISP = {'O->I': 'patch_dispersal_within_from_offspring_pool_to_immigrant_pool',
        'OD->I': 'patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool',
        'OD->sites': 'patch_dipersal_within_from_offspring_and_dormancy_pool'}
MARK = 'patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool'
GERM = ['patch_local_germinate_from_offspring_and_dormancy_pool', 'patch_local_germinate_from_offspring_and_immigrant_pool', 'patch_local_germinate_from_offspring_immigrant_and_dormancy_pool']
DORM = ['patch_dormancy_process_from_offspring_pool_to_dormancy_pool', 'patch_dormancy_process_from_offspring_pool_and_immigrant_pool']
ORIG_BODY = original_patch(list(DISP.values()) + [MARK] + GERM + DORM + ['patch_dead_selection', 'patch_disturbance_process'])   # original bodies on the current habitat


def random_patch(seed, cls=patch, n=None, L=None, W=None, K=None, fill=None):
    ''' reproducible random patch of fakes: shape, occupancy (scattered sites), pools and markers all drawn from `seed` '''
    rng = random.Random(seed)
    n = n if n is not None else rng.choice((2, 3, 4, 5))
    L = L if L is not None else rng.choice((1, 2, 3, 5)); W = W if W is not None else rng.choice((1, 3, 4))
    K = K if K is not None else rng.choice((0, 3, 20))
    np.random.seed(seed)
    p = make_patch(n, L=L, W=W, K=K, cls=cls)
    for h in habs(p):
        sites = [(i, j) for i in range(L) for j in range(W)]
        f = fill if fill is not None else rng.random()
        pos = [s for s in sites if rng.random() < f]
        rng.shuffle(pos)
        fill_fakes(h, positions=pos, O=rng.choice((0, 1, 4, 9, 25)), I=rng.choice((0, 2, 7)), D=rng.choice((0, 3, 11)), markers=(rng.randint(0, 12), rng.randint(0, 3)))
    return p


# ================================================================ 1. six clear-ups
CLEAR = {'patch_clear_up_microsite_individuals': 'microsites', 'patch_clear_up_offspring_pool': 'offspring_pool', 'patch_clear_up_immigrant_pool': 'immigrant_pool',
         'patch_clear_up_dormancy_pool': 'dormancy_pool', 'patch_clear_up_offspring_marker_pool': 'offspring_marker_pool', 'patch_clear_up_immigrant_marker_pool': 'immigrant_marker_pool'}
n1 = 0
for seed in range(30):
    for fn, what in CLEAR.items():
        p = random_patch(seed); clear_recycled()
        before = snapshot(p)
        if what == 'microsites': content = [x.label for h in habs(p) for row in h.set['microsite_individuals'] for x in row if x is not None]
        elif what in POOLS: content = [x.label for h in habs(p) for x in getattr(h, what)]
        else: content = [m for h in habs(p) for m in getattr(h, what)]
        r = getattr(p, fn)()
        after = snapshot(p)
        tag = '%s seed=%d' % (fn, seed)
        assert r == len(content), tag + ' return != pool size'
        for b, a, h in zip(before, after, habs(p)):
            if what == 'microsites':
                assert a[1] == [[None] * h.width for _ in range(h.length)] and a[5] == () and a[6] == 0 and a[7] == '{}' and sorted(a[4]) == sorted(b[4] + b[5]), tag
                assert a[2:4] == b[2:4] and h.asexual_parent_pos_ls == [] and h.species_category_for_sexual_parents_pos == {}, tag + ' other parts touched'
            else:
                k = list(POOLS + MARKER_POOLS).index(what)
                pools_b, pools_a = list(b[2]) + list(b[3]), list(a[2]) + list(a[3])
                assert pools_a[k] == (), tag + ' pool not emptied'
                assert [x for i, x in enumerate(pools_a) if i != k] == [x for i, x in enumerate(pools_b) if i != k] and a[:2] == b[:2] and a[4:] == b[4:], tag + ' other parts touched'
        if what in MARKER_POOLS: assert RECYCLED == [], tag + ' markers recycled something'
        else: assert RECYCLED == content, tag + ' recycled != former contents'
        assert getattr(p, fn)() == 0, tag + ' second call != 0'
        check_invariants(p, None, tag)
        n1 += 1
out['clear_up_cases'] = n1
print('1 six clear-ups: %d cases; return == pool size, only that pool emptied, recycled == former contents (in order), markers recycle nothing' % n1)

# ================================================================ 2. disturbance
n2 = 0
for seed in range(40):
    p, o = random_patch(seed), random_patch(seed, cls=OP)
    clear_recycled()
    exp = sum(h.indi_num + len(h.offspring_pool) + len(h.immigrant_pool) + len(h.offspring_marker_pool) + len(h.immigrant_marker_pool) for h in habs(p))
    dorm = [list(h.dormancy_pool) for h in habs(p)]
    cut = [x.label for h in habs(p) for x in [y for row in h.set['microsite_individuals'] for y in row if y is not None] + h.offspring_pool + h.immigrant_pool]
    st = rng_state()
    r = p.patch_disturbance_process()
    assert same_rng_state(st, rng_state()), 'disturbance used RNG'
    ro = o.patch_disturbance_process()
    assert r == exp and ro is None, (seed, r, exp)
    assert all(h.dormancy_pool == d and all(x.is_alive for x in d) for h, d in zip(habs(p), dorm)), 'dormancy touched'
    assert all(h.offspring_marker_pool == [] and h.immigrant_marker_pool == [] for h in habs(p)), 'marker pools not emptied'
    sp, so = snapshot(p), snapshot(o)
    assert [s[:3] + s[4:] for s in sp] == [s[:3] + s[4:] for s in so], 'disturbance != original apart from markers'
    assert sorted(RECYCLED) == sorted(cut), 'recycled != microsites + O + I'
    assert p.patch_disturbance_process() == 0
    check_invariants(p, None, 'disturbance %d' % seed)
    n2 += 1
out['disturbance_cases'] = n2
print('2 disturbance: %d cases; return == sum of habitat counts, dormancy untouched, marker pools emptied, == original apart from markers (original returns None)' % n2)

# ================================================================ 3. is_remove=False guards
n3 = 0
for seed in range(20):
    for key, fn in DISP.items():
        p = random_patch(seed, n=3); clear_recycled()
        b, st = snapshot(p), rng_state()
        try: getattr(p, fn)(0.5, False); raise AssertionError('%s accepted is_remove=False' % fn)
        except NotImplementedError: pass
        try: getattr(p, fn)(0.5, is_remove=False); raise AssertionError('%s accepted is_remove=False (kw)' % fn)
        except NotImplementedError: pass
        assert snapshot(p) == b and same_rng_state(st, rng_state()) and RECYCLED == [], key + ' guard changed state / RNG'
        n3 += 1
    for fn in GERM:
        p = random_patch(seed, n=3); b, st = snapshot(p), rng_state()
        try: getattr(p, fn)(False); raise AssertionError('%s accepted is_remove=False' % fn)
        except NotImplementedError: pass
        assert snapshot(p) == b and same_rng_state(st, rng_state()), fn + ' guard changed state / RNG'
        n3 += 1
out['is_remove_false_cases'] = n3
print('3 is_remove=False: %d cases; the three individual-pool dispersals raise NotImplementedError before any RNG use / state change' % n3)
# germinate(False) raising inside the first habitat: did an earlier habitat germinate? (habitat raises per call, patch loops habitats)
p = random_patch(3, n=3, fill=0.3); b = snapshot(p)
try: p.patch_local_germinate_from_offspring_and_immigrant_pool(False)
except NotImplementedError: pass
out['germinate_false_state_unchanged'] = snapshot(p) == b
assert out['germinate_false_state_unchanged'], 'patch germinate(is_remove=False) changed state before raising'

# ================================================================ 4. equivalence with the original
def run_pair(seed, fn, args_new, args_old, cls_old, rate_seed=None, mk=None):
    mk = mk or (lambda c: random_patch(seed, cls=c))
    pn, po = mk(patch), mk(cls_old)
    hb = held_labels(pn)
    clear_recycled()
    random.seed(1000 + seed); np.random.seed(2000 + seed)
    rn = getattr(pn, fn)(*args_new); sn = rng_state(); rec = list(RECYCLED)
    random.seed(1000 + seed); np.random.seed(2000 + seed)
    ro = getattr(po, fn)(*args_old); so = rng_state()
    return pn, po, rn, ro, same_rng_state(sn, so), rec, hb


eq = {k: 0 for k in list(DISP) + ['marker True', 'marker False'] + ['germ', 'dorm', 'dead', 'marker_birth', 'unchanged_helpers']}
extras_seen = 0
for seed in range(400):
    rate = random.Random(seed).choice((0.0, 0.1, 0.37, 0.5, 0.9, 1.0))
    for key, fn in DISP.items():
        for cls_old in (ORIG_BODY, OP):
            pn, po, rn, ro, same, rec, hb = run_pair(seed, fn, (rate,), (rate, True), cls_old)
            tag = '%s seed=%d rate=%g vs %s' % (key, seed, rate, cls_old.__name__)
            assert same, tag + ' RNG streams differ'
            assert snapshot(pn) == snapshot(po) and rn == ro, tag + ' state / return differ'
            if key != 'OD->sites':
                assert rec == [], tag + ' recycled something'
            else:
                dropped = hb - held_labels(po)
                assert set(rec) == dropped and len(rec) == len(dropped), tag + ' recycled != dropped by the original'
                assert not (set(rec) & held_labels(pn)), tag + ' recycled individual still held'
                extras_seen += len(rec)
            check_invariants(pn, None, tag)
        eq[key] += 1
    for mode in (True, False):
        for cls_old in (ORIG_BODY, OP):
            pn, po, rn, ro, same, rec, _ = run_pair(seed, MARK, (rate, mode), (rate, mode), cls_old)
            assert same and snapshot(pn) == snapshot(po) and rn == ro and rec == [], 'marker %s seed=%d' % (mode, seed)
        eq['marker %s' % mode] += 1
    if seed < 150:
        for fn in GERM:
            for cls_old in (ORIG_BODY, OP):
                pn, po, rn, ro, same, rec, _ = run_pair(seed, fn, (), (True,), cls_old)
                assert same and snapshot(pn) == snapshot(po) and rn == ro and rec == [], '%s seed=%d vs %s' % (fn, seed, cls_old.__name__)
            eq['germ'] += 1
        for fn in DORM:
            pn, po, rn, ro, same, rec, hb = run_pair(seed, fn, (), (), ORIG_BODY)          # same habitat layer -> identical incl. return + recycled
            assert same and snapshot(pn) == snapshot(po) and rn == ro, '%s seed=%d vs body' % (fn, seed)
            pn, po, rn, ro, same, rec, hb = run_pair(seed, fn, (), (), OP)                 # original habitat layer -> same pools, eliminate count may differ
            k0 = 0
            for a, b, h in zip(snapshot(pn), snapshot(po), habs(po)):
                if fn == DORM[0] and h.dormancy_pool_max_size == 0 and b[2][2]:                # habitat layer: K == 0 now empties a non-empty dormancy pool (original kept it)
                    assert a[2] == ((), (), ()) and b[2][:2] == ((), ()), fn
                    a, b, k0 = a[:2] + a[3:], b[:2] + b[3:], k0 + len(b[2][2])
                assert a == b and same, '%s seed=%d vs full original' % (fn, seed)
            assert rn[0] == ro[0] and rn[2] == ro[2] and rn[3] == ro[3] - k0, '%s seed=%d returns vs full original' % (fn, seed)
            eq['dorm_K0_dormants_dropped_vs_original'] = eq.get('dorm_K0_dormants_dropped_vs_original', 0) + k0
            assert set(rec) == hb - held_labels(pn), '%s seed=%d recycled != dropped' % (fn, seed)
            eq['dorm'] += 1
        for d, w, m in ((0.1, 0.5, 'niche_gaussian'), (0.3, 0.5, 'neutral_uniform')):
            pn, po, rn, ro, same, rec, hb = run_pair(seed, 'patch_dead_selection', (d, w, m), (d, w, m), OP)
            assert same and snapshot(pn) == snapshot(po) and rn == ro == len(rec) and set(rec) == hb - held_labels(pn), 'dead seed=%d' % seed
            eq['dead'] += 1
        pn, po, rn, ro, same, rec, _ = run_pair(seed, 'patch_sex_reproduce_calculation_into_offspring_marker_pool', (rate,), (rate,), OP)
        assert same and snapshot(pn) == snapshot(po) and rn == ro, 'marker birth seed=%d' % seed
        eq['marker_birth'] += 1
        # unchanged helpers / getters on identical states
        pn, po = random_patch(seed), random_patch(seed, cls=OP)
        random.seed(seed); np.random.seed(seed)
        vn = [pn.get_patch_individual_num(), pn.get_patch_empty_sites_ls(), pn.patch_empty_sites_num(), pn.get_patch_pairwise_empty_sites_ls(), pn.get_patch_pairwise_occupied_sites_ls(),
              pn.get_patch_occupied_sites_ls(), len(pn.get_patch_offspring_marker_pool()), [x.label for x in pn.get_patch_offspring_pool()], [x.label for x in pn.get_patch_dormancy_pool()],
              [x.label for x in pn.get_patch_offspring_and_dormancy_pool()], pn.get_patch_habs_offspring_num_matrix().tolist(), pn.get_patch_habs_offspring_marker_num_matrix().tolist(),
              pn.get_patch_habs_dormancy_num_matrix().tolist(), pn.get_patch_habs_empty_sites_num_matrix().tolist(), [x.label for x in pn.sample_offspring_without_replacement(3, ('offspring_pool', 'dormancy_pool'))]]
        sn_ = rng_state(); random.seed(seed); np.random.seed(seed)
        vo = [po.get_patch_individual_num(), po.get_patch_empty_sites_ls(), po.patch_empty_sites_num(), po.get_patch_pairwise_empty_sites_ls(), po.get_patch_pairwise_occupied_sites_ls(),
              po.get_patch_occupied_sites_ls(), len(po.get_patch_offspring_marker_pool()), [x.label for x in po.get_patch_offspring_pool()], [x.label for x in po.get_patch_dormancy_pool()],
              [x.label for x in po.get_patch_offspring_and_dormancy_pool()], po.get_patch_habs_offspring_num_matrix().tolist(), po.get_patch_habs_offspring_marker_num_matrix().tolist(),
              po.get_patch_habs_dormancy_num_matrix().tolist(), po.get_patch_habs_empty_sites_num_matrix().tolist(), [x.label for x in po.sample_offspring_without_replacement(3, ('offspring_pool', 'dormancy_pool'))]]
        assert vn == vo and same_rng_state(sn_, rng_state()) and snapshot(pn) == snapshot(po), 'helpers seed=%d' % seed
        if all(h.size == habs(pn)[0].size for h in habs(pn)):
            assert np.array_equal(pn.get_patch_microsites_individals_sp_id_values(), po.get_patch_microsites_individals_sp_id_values(), equal_nan=True)
            assert np.array_equal(pn.get_patch_microsites_environment_values('env1'), po.get_patch_microsites_environment_values('env1'))
            assert np.array_equal(pn.get_patch_microsites_optimum_sp_id_value_array(0.1, 0.5, [[0.2, 0.2], [0.5, 0.5]]), po.get_patch_microsites_optimum_sp_id_value_array(0.1, 0.5, [[0.2, 0.2], [0.5, 0.5]]))
        mig = [C.fake_indi('mig%d' % i) for i in range(7)]
        random.seed(seed); pn.migrants_to_patch_into_habs_immigrant_pool(list(mig)); pn.migrants_marker_to_patch_into_habs_immigrant_marker_pool([('q', 'h', 's')] * 5)
        random.seed(seed); po.migrants_to_patch_into_habs_immigrant_pool(list(mig)); po.migrants_marker_to_patch_into_habs_immigrant_marker_pool([('q', 'h', 's')] * 5)
        assert snapshot(pn) == snapshot(po), 'migrants_to_patch seed=%d' % seed
        eq['unchanged_helpers'] += 1
out['equivalence_cases'] = eq
out['OD->sites_extras_recycled_total'] = extras_seen
print('4 equivalence with %s (original bodies on the current habitat + full original stack): %s; OD->sites recycled == dropped by the original (%d extras over all seeds)' % (C.ORIGINAL_COMMIT, eq, extras_seen))
assert extras_seen > 0, 'equivalence sweep never produced an extra migrant'

# ================================================================ 5. forced extras
def expected_counts(p, M):
    avail = [len(h.offspring_pool) + len(h.dormancy_pool) for h in habs(p)]
    empties = [len(h.empty_site_pos_ls) for h in habs(p)]
    placed, extras = [], []
    for j in range(p.hab_num):
        mig = 0
        for i in range(p.hab_num):
            if i == j: continue
            t = min(int(M[i][j]), avail[i]); avail[i] -= t; mig += t
        placed.append(min(empties[j], mig)); extras.append(mig - placed[-1])
    return placed, extras, empties


n5, tot_extra = 0, 0
for seed in range(200):
    rng = random.Random(seed)
    p = random_patch(seed, n=rng.choice((2, 3, 4)), L=3, W=3, fill=rng.choice((0.5, 0.8, 1.0)))
    M = [[0 if i == j else rng.randint(0, 9) for j in range(p.hab_num)] for i in range(p.hab_num)]
    if seed % 2 == 0:
        for h in habs(p): h.offspring_pool += [C.fake_indi('%s.x%d' % (h.name, k)) for k in range(30)]
        M = [[0 if i == j else 9 for j in range(p.hab_num)] for i in range(p.hab_num)]
    placed, extras, empties = expected_counts(p, M)
    force_within_matrix(p, M); clear_recycled()
    pool_b = {x.label for h in habs(p) for q in ('offspring_pool', 'dormancy_pool') for x in getattr(h, q)}
    sites_b = {x.label for h in habs(p) for row in h.set['microsite_individuals'] for x in row if x is not None}
    r = p.patch_dipersal_within_from_offspring_and_dormancy_pool(0.5)
    sites_a = {x.label for h in habs(p) for row in h.set['microsite_individuals'] for x in row if x is not None}
    pool_a = {x.label for h in habs(p) for q in ('offspring_pool', 'dormancy_pool') for x in getattr(h, q)}
    new_in_sites = sites_a - sites_b
    assert r == sum(placed) == len(new_in_sites) and len(RECYCLED) == sum(extras), (seed, r, placed, extras)
    assert [e - len(h.empty_site_pos_ls) for e, h in zip(empties, habs(p))] == placed
    assert set(RECYCLED) | new_in_sites == pool_b - pool_a and not (set(RECYCLED) & held_labels(p)), 'taken != placed + recycled'
    if seed % 2 == 0: assert r == sum(empties), 'migrants >= empties everywhere but return != empty sites'
    check_invariants(p, None, 'forced %d' % seed)
    n5 += 1; tot_extra += len(RECYCLED)
out['forced_extras'] = {'cases': n5, 'extras_recycled': tot_extra}
print('5 forced extras: %d cases, %d extras recycled; return == empties filled (== all empty sites when migrants >= empties), invariants hold' % (n5, tot_extra))

# ================================================================ 6. natural extras
N6, freq = 600, {'extras1': 0, 'arrive0': 0, 'arrive1': 0, 'arrive2': 0}
for seed in range(N6):
    np.random.seed(7)
    p = make_patch(3, L=3, K=0)
    h0, h1, h2 = habs(p)
    fill_fakes(h0, occ=8); fill_fakes(h1, occ=9, O=100); fill_fakes(h2, occ=9, O=100)
    clear_recycled()
    random.seed(seed); np.random.seed(10000 + seed)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r = p.patch_dipersal_within_from_offspring_and_dormancy_pool(0.5)
    arrivals = r + len(RECYCLED)
    assert r <= 1 and len(RECYCLED) <= 1 and h0.indi_num == 8 + r and len(h1.offspring_pool) + len(h2.offspring_pool) == 200 - arrivals
    freq['extras1'] += len(RECYCLED) == 1; freq['arrive%d' % arrivals] += 1
    check_invariants(p, None, 'natural %d' % seed)
ok6 = binom_ok(freq['extras1'], N6, 0.25) and binom_ok(freq['arrive0'], N6, 0.25) and binom_ok(freq['arrive1'], N6, 0.5)
out['natural_extras'] = dict(freq, seeds=N6, p_extras=round(freq['extras1'] / N6, 3), expected=0.25)
assert ok6, out['natural_extras']
print('6 natural extras over %d seeds: P(extra) = %.3f (expected 0.25, 4-sigma), arrivals 0/1/2 = %d/%d/%d; invariants every time' % (N6, freq['extras1'] / N6, freq['arrive0'], freq['arrive1'], freq['arrive2']))

# ================================================================ 7. sampling uniformity
T = 3000
pick, site_hits, pick_i = {}, {}, {}
for t in range(T):
    np.random.seed(t % 5)
    p = make_patch(2, L=3, K=0)
    fill_fakes(habs(p)[0], occ=9, O=10, D=5); fill_fakes(habs(p)[1], occ=3)
    force_within_matrix(p, [[0, 4], [0, 0]]); clear_recycled(); random.seed(20000 + t)
    p.patch_dipersal_within_from_offspring_and_dormancy_pool(0.5)
    for (i, j) in habs(p)[1].occupied_site_pos_ls:
        x = habs(p)[1].set['microsite_individuals'][i][j]
        if not x.label.startswith('h1.m'): pick[x.label] = pick.get(x.label, 0) + 1; site_hits[(i, j)] = site_hits.get((i, j), 0) + 1
    q = make_patch(3, L=2, K=0)
    fill_fakes(habs(q)[0], O=12); force_within_matrix(q, [[0, 2, 1], [0, 0, 0], [0, 0, 0]]); random.seed(40000 + t)
    q.patch_dispersal_within_from_offspring_pool_to_immigrant_pool(0.5)
    for h in habs(q)[1:]:
        for x in h.immigrant_pool: pick_i[x.label] = pick_i.get(x.label, 0) + 1
cand = ['h0.o%d' % i for i in range(10)] + ['h0.d%d' % i for i in range(5)]
ok_a, s_a = chi2_ok([pick.get(l, 0) for l in cand], [T * 4 / 15] * 15)
empt = [(i, j) for i in range(3) for j in range(3)][3:]
ok_b, s_b = chi2_ok([site_hits.get(s, 0) for s in empt], [T * 4 / 6] * 6)
ok_c, s_c = chi2_ok([pick_i.get('h0.o%d' % i, 0) for i in range(12)], [T * 3 / 12] * 12)
out['uniformity_chi2'] = {'OD->sites migrants (df14)': s_a, 'OD->sites destination sites (df5)': s_b, 'O->I migrants (df11)': s_c}
assert ok_a and ok_b and ok_c, out['uniformity_chi2']
print('7 uniformity chi2:', out['uniformity_chi2'])

# ================================================================ 8. edge shapes
def every_function(p):
    calls = [('dead_selection', lambda: p.patch_dead_selection(0.1, 0.5, 'niche_gaussian'))]
    for rate in (0.0, 1.0, 0.4):
        calls += [('%s m=%g' % (k, rate), (lambda fn, r: lambda: getattr(p, fn)(r))(fn, rate)) for k, fn in DISP.items()]
        calls += [('marker m=%g %s' % (rate, mode), (lambda r, mo: lambda: p.patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool(r, mo))(rate, mode)) for mode in (True, False)]
    calls += [(fn[6:], getattr(p, fn)) for fn in GERM + DORM]
    calls += [('marker_birth', lambda: p.patch_sex_reproduce_calculation_into_offspring_marker_pool(0.7)),
              ('migrants_to_patch', lambda: p.migrants_to_patch_into_habs_immigrant_pool([C.fake_indi('in%d' % k) for k in range(3)])),
              ('migrants_marker_to_patch', lambda: p.migrants_marker_to_patch_into_habs_immigrant_marker_pool([('q', 'h', 's')] * 3)),
              ('sp_id_values', p.get_patch_microsites_individals_sp_id_values), ('env_values', lambda: p.get_patch_microsites_environment_values('env1')),
              ('optimum_sp_id', lambda: p.get_patch_microsites_optimum_sp_id_value_array(0.1, 0.5, [[0.2, 0.2]])),
              ('rate_matrix m=0.3', lambda: p.get_patch_dispersal_within_rate_matrix(0.3)), ('offset_env', lambda: p.patch_offset_environmental_values(['env1'], [0.0])),
              ('disturbance', p.patch_disturbance_process)] + [(fn[6:], getattr(p, fn)) for fn in CLEAR]
    res = {}
    for name, fn in calls:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            try: r = fn(); res[name] = 'ok'
            except Exception as e: res[name] = '%s: %s' % (type(e).__name__, str(e)[:60])
        msgs = sorted({str(x.message)[:40] for x in w if 'matrix subclass' not in str(x.message)})
        MATRIX_WARN[0] += sum('matrix subclass' in str(x.message) for x in w)
        if msgs: res[name] += ' [warn: %s]' % '; '.join(msgs)
        check_invariants(p, None, name)
    return res


edges, MATRIX_WARN = {}, [0]
for label, kw in (('1 habitat 3x3', dict(n=1, L=3, W=3)), ('2 habitats 3x3', dict(n=2, L=3, W=3)), ('3 habitats 2x5 (non-square)', dict(n=3, L=2, W=5)),
                  ('4 habitats 1x1', dict(n=4, L=1, W=1)), ('2 habitats full', dict(n=2, L=2, W=2, fill=1.0)), ('3 habitats empty', dict(n=3, L=3, W=3, fill=0.0))):
    edges[label] = every_function(random_patch(5, **kw))
np.random.seed(0); p0 = make_patch(0)
edges['empty patch (0 habitats)'] = every_function(p0)
np.random.seed(0); pu = make_patch(2, sizes=[(2, 2), (3, 3)]); [fill_fakes(h, occ=2, O=3) for h in habs(pu)]
edges['unequal sizes 2x2 + 3x3'] = every_function(pu)
np.random.seed(0); pi = patch('p', 0, (0, 0))
for k in (1, 2): pi.add_habitat('h%d' % k, k, (0, k), 2, ['env1', 'env2'], [0.2, 0.2], [0.025, 0.025], 2, 2)
try: pi.migrants_to_patch_into_habs_immigrant_pool([C.fake_indi('a'), C.fake_indi('b')]); split_res = 'ok'
except Exception as e: split_res = type(e).__name__
out['edges'] = edges; out['split_hab_index_1_2'] = split_res
for label, res in edges.items():
    bad = {k: v for k, v in res.items() if not v.startswith('ok')}
    warn = sorted({k.split(' ')[0] + ': ' + v.split('[warn: ')[1].rstrip(']') for k, v in res.items() if '[warn: ' in v})
    print('8 %-28s errors %s | warnings %s' % (label, bad if bad else 'none', warn if warn else 'none'))
# pre-existing behaviour (verify + record, not regressions)
u = edges['unequal sizes 2x2 + 3x3']
assert all(u[k].startswith('ValueError') for k in ('sp_id_values', 'env_values', 'optimum_sp_id')), u
one = edges['1 habitat 3x3']
assert 'divide by zero' in one['rate_matrix m=0.3'] and one['OD->sites m=0.4'].startswith('ok'), one
assert split_res == 'IndexError'
e0 = edges['empty patch (0 habitats)']
finding('known', 'get_patch_microsites_individals_sp_id_values / _phenotype_values / get_patch_microsites_environment_values / get_patch_microsites_optimum_sp_id_value_array reshape with the LAST habitat\'s size: unequal habitat sizes -> %s; an empty patch -> %s (pre-existing, unchanged from the original).' % (u['sp_id_values'].split(':')[0], e0['sp_id_values'].split(':')[0]))
finding('known', 'get_patch_dispersal_within_rate_matrix divides by hab_num-1: with 1 habitat numpy only warns (%s) and the 1x1 matrix is still right because the diagonal is overwritten; every dispersal returns 0 (pre-existing).' % one['rate_matrix m=0.3'].split('[warn: ')[-1].rstrip(']'))
finding('known', 'split() / migrants_to_patch_into_habs_immigrant_pool index pieces by h_object.index: hab_index 1,2 (not 0..n-1) -> %s; an empty patch -> %s (pre-existing).' % (split_res, e0['migrants_to_patch'].split(':')[0]))
finding('info', 'patch_dispersal_within_*: when a column of the emigrant matrix sums to 0 (no offspring / dormants anywhere else) get_habs_immigrants_matrix divides 0/0 -> numpy RuntimeWarning, NaN entries are rounded to 0 by patch_matrix_around (pre-existing, harmless; seen in %d of %d edge shapes).' % (sum(any('invalid value' in v for v in r.values()) for r in edges.values()), len(edges)))
finding('info', 'np.matrix is used by every dispersal matrix helper: numpy emits PendingDeprecationWarning "the matrix subclass is not the recommended way" (%d times in the edge section; pre-existing, will break when numpy removes np.matrix).' % MATRIX_WARN[0])
finding('info', 'habitat layer (not patch): hab_dormancy_process_from_offspring_pool_to_dormancy_pool with dormancy_pool_max_size == 0 now empties (recycles) a non-empty dormancy_pool; the original kept it and reported eliminate=0 (%d dormants over the equivalence sweep). Only reachable when K is set to 0 with dormants present; intended by the habitat migration (comment: K changed mid-run).' % eq['dorm_K0_dormants_dropped_vs_original'])
import inspect as _inspect
_sig = {fn: _inspect.signature(getattr(patch, fn)).parameters['is_remove'].default for fn in list(DISP.values()) + [MARK]}
_defaults = {k: ('required' if v is _inspect.Parameter.empty else v) for k, v in zip(list(DISP) + ['marker'], _sig.values())}
assert _defaults == {'O->I': True, 'OD->I': True, 'OD->sites': True, 'marker': True}, 'is_remove defaults of the four within-patch dispersals: %s' % _defaults
out['is_remove_defaults'] = _defaults
finding('known', 'Marker dispersal with is_remove=False keeps the marker at the source and sends a copy (model semantics; can inflate births). Both modes are identical to the original.')
save('bookkeeping', out)
print('PASSED')
