# Pool / microsite bookkeeping with fake individuals (no Nemo):
#   1 germination == original   2 dormancy: exact table, random sweep, == original   3 clear-up + disturbance (== original), two-owner detection
#   4 dead selection == original   5 sampling uniformity (germination, dormancy, sample_offspring_without_replacement)
import math, random
import numpy as np
import common as C
from common import habitat, fake_indi, make_hab, labels, check_invariants, microsite_individuals, original_habitat, save, chi2_ok, POOLS

out = {}
ORIG = original_habitat(['hab_local_germinate_from_offspring_and_dormancy_pool', 'hab_local_germinate_from_offspring_and_immigrant_pool',
                         'hab_local_germinate_from_offspring_immigrant_and_dormancy_pool', 'hab_dormancy_process_from_offspring_pool_to_dormancy_pool',
                         'hab_dormancy_process_from_offspring_pool_and_immigrant_pool', 'habitat_disturbance_process', 'hab_dead_selection'])


# ================================================================ 1. germination == original
GERM = {'O+D': ('hab_local_germinate_from_offspring_and_dormancy_pool', ('offspring_pool', 'dormancy_pool')),
        'O+I': ('hab_local_germinate_from_offspring_and_immigrant_pool', ('offspring_pool', 'immigrant_pool')),
        'O+I+D': ('hab_local_germinate_from_offspring_immigrant_and_dormancy_pool', POOLS)}


def germ_hab(cls, seed, occupied, sizes):
    np.random.seed(seed)
    h = make_hab(L=6, K=1000, cls=cls)
    for k in range(occupied): h.add_individual(fake_indi('m%d' % k, 'sp1', 'female' if k % 2 else 'male'), k // 6, k % 6)
    for p, n in zip(POOLS, sizes): setattr(h, p, [fake_indi('%s%d' % (p[0], i), 'sp%d' % (1 + i % 2), 'female' if i % 3 else 'male') for i in range(n)])
    return h


n_cases = 0
for key, (meth, used) in GERM.items():
    sz = lambda per, drop=None: tuple((0 if p == drop else per) if p in used else 3 for p in POOLS)
    per = 12 // len(used)
    scen = [(20, sz(2)), (36 - per * len(used), sz(per)), (30, sz(10)), (36, sz(5)), (26, sz(0))] + [(26, sz(8, p)) for p in used]
    for occupied, sizes in scen:
        for seed in range(100, 125):
            new, old = germ_hab(habitat, seed, occupied, sizes), germ_hab(ORIG, seed, occupied, sizes)
            random.seed(seed); rn = getattr(new, meth)(); sn = random.getstate()
            random.seed(seed); ro = getattr(old, meth)(True); so = random.getstate()
            tag = '%s occ=%d %s seed=%d' % (key, occupied, sizes, seed)
            assert rn == ro and sn == so and labels(new) == labels(old), tag
            assert all([x.label for x in getattr(new, p)] == [x.label for x in getattr(old, p)] for p in POOLS), tag + ' pools'
            sites = set(map(id, microsite_individuals(new)))
            assert not any(id(x) in sites for p in POOLS for x in getattr(new, p)), tag + ' held twice'
            check_invariants(new, None, tag)
            n_cases += 1
    h = germ_hab(habitat, 1, 20, (4, 4, 4))
    rs, ns, empty, pools = random.getstate(), np.random.get_state(), list(h.empty_site_pos_ls), {p: list(getattr(h, p)) for p in POOLS}
    try: getattr(h, meth)(is_remove=False); raise AssertionError('is_remove=False accepted')
    except NotImplementedError: pass
    assert random.getstate() == rs and np.array_equal(np.random.get_state()[1], ns[1]) and h.empty_site_pos_ls == empty
    assert all(getattr(h, p) == pools[p] for p in POOLS)
out['germination_equivalence_cases'] = n_cases
print('1 germination == original: %d cases; is_remove=False touches nothing' % n_cases)

# ================================================================ 2. dormancy
A, B = 'hab_dormancy_process_from_offspring_pool_to_dormancy_pool', 'hab_dormancy_process_from_offspring_pool_and_immigrant_pool'


def dorm(cls, meth, K, D, O, I, seed=0):
    np.random.seed(seed)
    h = make_hab(L=2, K=K, cls=cls)
    h.dormancy_pool, h.offspring_pool, h.immigrant_pool = ([fake_indi('%s%d' % (c, i)) for i in range(n)] for c, n in (('d', D), ('o', O), ('i', I)))
    snap = (h.offspring_pool + h.immigrant_pool + h.dormancy_pool, list(h.dormancy_pool), list(h.offspring_pool), list(h.immigrant_pool))
    n0 = len(C.RECYCLED); random.seed(seed)
    r = getattr(h, meth)()
    return h, r, len(C.RECYCLED) - n0, snap


TABLE = [('A1', A, 0, 0, 10, 5, (0, 0, 0, 0), 15), ('A2', A, 100, 30, 150, 7, (0, 30, 100, 100), 87), ('A3', A, 100, 30, 40, 5, (30, 0, 40, 70), 5),
         ('A4', A, 100, 80, 40, 5, (60, 20, 40, 100), 25), ('B1', B, 100, 30, 150, 7, (0, 30, 100, 100), 87), ('B2', B, 100, 30, 60, 80, (0, 30, 100, 100), 70),
         ('B3', B, 100, 60, 30, 20, (50, 10, 50, 100), 10), ('B4', B, 100, 20, 30, 20, (20, 0, 50, 70), 0), ('B5', B, 0, 0, 10, 0, (0, 0, 0, 0), 10),
         ('B6', B, 0, 0, 0, 5, (0, 0, 0, 0), 5), ('B7', B, 0, 0, 0, 0, (0, 0, 0, 0), 0)]
out['dormancy_table'] = []
for name, meth, K, D, O, I, exp, exp_rec in TABLE:
    _, r, rec, _ = dorm(habitat, meth, K, D, O, I, 1)
    out['dormancy_table'].append([name, K, D, O, I, list(r), rec])
    assert tuple(r) == exp and rec == exp_rec, (name, r, rec)

rng = random.Random(2024)
for meth in (A, B):
    for t in range(2500):
        K = rng.randint(0, 200); D = rng.randint(0, K); O = rng.randint(0, 300); I = rng.randint(0, 300)
        h, r, rec, (before, old, offs, imms) = dorm(habitat, meth, K, D, O, I, t)
        after, ids = h.dormancy_pool, set(map(id, h.dormancy_pool))
        tag = '%s K=%d D=%d O=%d I=%d' % (meth[-14:], K, D, O, I)
        assert len(after) <= K and len(after) == r[3] and r[0] + r[1] == D and rec == len(before) - len(after), tag
        assert h.offspring_pool == [] and h.immigrant_pool == [] and len(ids) == len(after) and all(x.is_alive for x in after), tag
        assert all(x.is_alive == (id(x) in ids) for x in before), tag + ' alive != kept'
        if K and O <= K: assert all(id(x) in ids for x in offs), tag + ' offspring priority'
        if meth == B and O + I <= K: assert all(id(x) in ids for x in imms), tag + ' immigrant priority'
        if meth == A: assert not any(id(x) in ids for x in imms), tag + ' immigrant entered in A'
        assert sum(id(x) in set(map(id, old)) for x in after) == r[0], tag + ' survival'

diff = {}
rng = random.Random(77)
for meth in (A, B):
    for t in range(2500):
        K = rng.randint(0, 200); D = rng.randint(0, K); O = rng.randint(0, 300); I = rng.randint(0, 300)
        hn, rn, _, _ = dorm(habitat, meth, K, D, O, I, t)
        ho, ro, _, _ = dorm(ORIG, meth, K, D, O, I, t)
        assert [x.label for x in hn.dormancy_pool] == [x.label for x in ho.dormancy_pool] and rn[0] == ro[0] and rn[2:] == ro[2:], (meth, K, D, O, I)
        if rn[1] != ro[1]:
            br = ('A2' if O > K else 'A4') if meth == A else ('B1' if O > K else 'B2' if O + I > K else 'B3' if O + I + D > K else 'B4')
            assert D < K and br in ('A2', 'B1', 'B2', 'B3'), (meth, K, D, O, I, br)
            diff[br] = diff.get(br, 0) + 1
out['dormancy_eliminate_diff_vs_original'] = diff

out['runtime_K'] = {}
for K_new in (5, 0):
    for meth in (A, B):
        h = make_hab(L=2, K=20)
        h.dormancy_pool, h.offspring_pool = [fake_indi('d%d' % i) for i in range(15)], [fake_indi('o%d' % i) for i in range(2)]
        h.dormancy_pool_max_size = K_new
        n0 = len(C.RECYCLED); random.seed(0); r = getattr(h, meth)()
        assert len(h.dormancy_pool) <= K_new and r[0] + r[1] == 15 and len(C.RECYCLED) - n0 == 17 - len(h.dormancy_pool), (K_new, meth, r)
        out['runtime_K']['%s K20->%d' % (meth[-14:], K_new)] = list(r)
print('2 dormancy: table 11/11, sweep 2x2500, == original 2x2500 (eliminate diff %s), runtime K change consistent' % diff)

# ================================================================ 3. clear-up, disturbance, two owners
def filled(cls=habitat, seed=0, occ=20, O=7, I=5, D=4, marks=(3, 2)):
    np.random.seed(seed)
    h = make_hab(L=6, K=50, cls=cls)
    for k in range(occ): h.add_individual(fake_indi('m%d' % k, 'sp%d' % (1 + k % 2), 'female' if k % 3 else 'male'), k // 6, k % 6)
    h.offspring_pool, h.immigrant_pool, h.dormancy_pool = [fake_indi('o%d' % i) for i in range(O)], [fake_indi('i%d' % i) for i in range(I)], [fake_indi('d%d' % i) for i in range(D)]
    h.offspring_marker_pool, h.immigrant_marker_pool = [('p', 'h1', 'sexual')] * marks[0], [('p', 'h0', 'sexual')] * marks[1]
    h.asexual_parent_pos_ls, h.species_category_for_sexual_parents_pos = [(0, 0)], {'sp1': {'female': [(0, 2)]}}
    return h


for fn, pool, n in (('hab_clear_up_offspring_pool', 'offspring_pool', 7), ('hab_clear_up_immigrant_pool', 'immigrant_pool', 5), ('hab_clear_up_dormancy_pool', 'dormancy_pool', 4)):
    h = filled(); xs = list(getattr(h, pool)); others = {p: list(getattr(h, p)) for p in POOLS if p != pool}
    n0 = len(C.RECYCLED)
    assert getattr(h, fn)() == n and len(C.RECYCLED) - n0 == n and getattr(h, pool) == [] and not any(x.is_alive for x in xs)
    assert all(getattr(h, p) == v for p, v in others.items()) and h.indi_num == 20 and getattr(h, fn)() == 0
    check_invariants(h, None, fn)
h = filled(); ms = microsite_individuals(h)
assert h.hab_clear_up_microsite_individuals() == 20 and not any(x.is_alive for x in ms) and h.indi_num == 0 and h.occupied_site_pos_ls == []
assert h.species_category == {} and h.asexual_parent_pos_ls == [] and h.species_category_for_sexual_parents_pos == {} and len(h.offspring_pool) == 7
check_invariants(h, None, 'clear micro')
h = filled(); n0 = len(C.RECYCLED)
assert h.hab_clear_up_offspring_marker_pool() == 3 and h.hab_clear_up_immigrant_marker_pool() == 2 and len(C.RECYCLED) == n0 and h.offspring_marker_pool == []

for seed in range(10):
    kw = dict(seed=seed, occ=5 + 3 * seed, O=seed, I=2 * seed % 7, D=seed % 4, marks=(seed % 3, seed % 2))
    h, o = filled(**kw), filled(cls=ORIG, **kw)
    dorm_before = list(h.dormancy_pool)
    r = h.habitat_disturbance_process(); assert o.habitat_disturbance_process() is None
    assert r == kw['occ'] + kw['O'] + kw['I'] + sum(kw['marks']) and h.habitat_disturbance_process() == 0
    assert h.dormancy_pool == dorm_before and all(x.is_alive for x in dorm_before)
    assert labels(h) == labels(o) and h.species_category == o.species_category and all([x.label for x in getattr(h, p)] == [x.label for x in getattr(o, p)] for p in POOLS)
    assert h.offspring_marker_pool == [] and len(o.offspring_marker_pool) == kw['marks'][0]       # only difference: new version empties the marker pools
    check_invariants(h, None, 'disturbance %d' % seed)

x = fake_indi('dup'); h = filled(O=0, I=0, D=0); h.offspring_pool, h.immigrant_pool = [x], [x]
try: h.hab_clear_up_offspring_pool(); h.hab_clear_up_immigrant_pool(); raise AssertionError('two owners not detected')
except RuntimeError: pass
x = fake_indi('dup'); h = make_hab(L=2, K=5); h.offspring_pool = [fake_indi('o%d' % i) for i in range(8)]; h.immigrant_pool = [x]; h.dormancy_pool = [x]
try: h.hab_dormancy_process_from_offspring_pool_to_dormancy_pool(); raise AssertionError('two owners not detected')
except RuntimeError: pass
print('3 clear-up x6, disturbance == original x10 (only marker pools differ), two-owner duplicates raise')

# ================================================================ 4. dead selection == original
def dgrid(cls, seed, fill):
    np.random.seed(seed); rng = random.Random(seed)
    h = make_hab(L=8, cls=cls)
    for i in range(8):
        for j in range(8):
            if rng.random() < fill:
                h.add_individual(fake_indi('m%d_%d' % (i, j), 'sp%d' % rng.randint(1, 2), rng.choice(('male', 'female')), (0.2 + rng.gauss(0, .3), 0.2 + rng.gauss(0, .3))), i, j)
    return h


n_dead = 0
for seed in range(40):
    for d, w, m in ((0.1, 0.5, 'niche_gaussian'), (0.0, 0.2, 'niche_gaussian'), (0.3, 0.5, 'neutral_uniform'), (0.0, 5.0, 'niche_gaussian')):
        for fill in (1.0, 0.6, 0.0):
            hn, ho = dgrid(habitat, seed, fill), dgrid(ORIG, seed, fill)
            before = {x.label: x for x in microsite_individuals(hn)}
            np.random.seed(seed + 1000); n0 = len(C.RECYCLED); rn = hn.hab_dead_selection(d, w, m); sn = np.random.get_state()
            np.random.seed(seed + 1000); ro = ho.hab_dead_selection(d, w, m); so = np.random.get_state()
            assert rn == ro and labels(hn) == labels(ho) and np.array_equal(sn[1], so[1]) and sn[2] == so[2]
            assert hn.asexual_parent_pos_ls == ho.asexual_parent_pos_ls and hn.species_category_for_sexual_parents_pos == ho.species_category_for_sexual_parents_pos
            dead = set(before) - {x.label for x in microsite_individuals(hn)}
            assert set(C.RECYCLED[n0:]) == dead and len(C.RECYCLED) - n0 == rn
            check_invariants(hn, None, 'dead')
            n_dead += 1
out['dead_selection_equivalence_cases'] = n_dead
print('4 dead selection == original: %d cases; recycled == dead' % n_dead)

# ================================================================ 5. sampling uniformity
T, hits = 3000, {}
for t in range(T):
    np.random.seed(t % 7)
    g = make_hab(L=5)
    for k in range(15): g.add_individual(fake_indi('m%d' % k), k // 5, k % 5)
    g.offspring_pool, g.immigrant_pool, g.dormancy_pool = ([fake_indi('%s%d' % (c, i)) for i in range(n)] for c, n in (('o', 30), ('i', 10), ('d', 10)))
    random.seed(10000 + t); g.hab_local_germinate_from_offspring_immigrant_and_dormancy_pool()
    for x in microsite_individuals(g):
        if x.label[0] != 'm': hits[x.label] = hits.get(x.label, 0) + 1
cand = ['%s%d' % (c, i) for c, n in (('o', 30), ('i', 10), ('d', 10)) for i in range(n)]
ok1, s1 = chi2_ok([hits.get(l, 0) for l in cand], [T * 10 / 50] * 50)
keep = {}
for t in range(T):
    g = make_hab(L=2, K=10); g.offspring_pool = [fake_indi('o%d' % i) for i in range(25)]
    random.seed(30000 + t); g.hab_dormancy_process_from_offspring_pool_to_dormancy_pool()
    for x in g.dormancy_pool: keep[x.label] = keep.get(x.label, 0) + 1
ok2, s2 = chi2_ok([keep.get('o%d' % i, 0) for i in range(25)], [T * 10 / 25] * 25)
pick = {}
for t in range(T):
    g = make_hab(L=2); g.offspring_pool, g.dormancy_pool = [fake_indi('o%d' % i) for i in range(6)], [fake_indi('d%d' % i) for i in range(4)]
    random.seed(50000 + t)
    for p in g.sample_offspring_without_replacement(3, ('offspring_pool', 'dormancy_pool')): pick[p.label] = pick.get(p.label, 0) + 1
    assert len(g.offspring_pool) + len(g.dormancy_pool) == 7
ok3, s3 = chi2_ok([pick.get(l, 0) for l in ['o%d' % i for i in range(6)] + ['d%d' % i for i in range(4)]], [T * 3 / 10] * 10)
out['uniformity_chi2'] = {'germination_candidates(df49)': s1, 'dormancy_A2(df24)': s2, 'sample_without_replacement(df9)': s3}
assert ok1 and ok2 and ok3, out['uniformity_chi2']
print('5 uniformity chi2:', out['uniformity_chi2'])
save('bookkeeping', out)
print('PASSED')
