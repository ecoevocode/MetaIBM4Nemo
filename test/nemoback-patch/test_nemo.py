# Patch-level tests with real Nemo individuals. Every section runs in its own fresh process (Nemo keeps process-wide static state):
#   initialize       fill, species table, homepatch, re-initialize (residents recycled, pools kept), bad inputs (atomicity across habitats)
#   reproduction     patch_sex_reproduce_mutate_into_offspring_pool + patch_sexual_birth_germinate: parents, factory, homepatch, per-locus inheritance, sex ratio
#   two_species      two species in different habitats + within-patch dispersal: mixing, never cross-species pairs, both factories balanced
#   lifecycle        v3-style loops over 12 configurations (K, dispersal variant, disturbance), invariants after every phase, return values vs counters
#   getters          data getters with real individuals (sp_id, phenotype, environment, optimum sp_id, get_data, matrices, pools, pairwise lists)
#   memory           long-run RSS, recycling-pool address reuse, reference individuals unchanged, recycled individuals inert, teardown balance
#   reproducibility  two fresh processes with the same seed -> identical nemo_digest over 50 steps
#   fuzz             >= 3000 random patch-level operations on two patches, invariants after each
#   edges_known      bad inputs, asex / mix probes (one subprocess each), known behaviour (recorded)
import os, sys, json, math, random, subprocess, traceback, warnings
import numpy as np
import common as C
from common import (patch, make_patch, habs, held, check_invariants, factories_balanced, recycle_everything, nemo_digest, nemo_species, seed_all,
                    chi2_ok, binom_ok, rss_mb, finding, save, HERE, POOLS)

warnings.filterwarnings('ignore', category=RuntimeWarning); warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
TABLE = [[0.2, 0.2], [0.5, 0.5]]
SECTIONS = ('initialize', 'reproduction', 'two_species', 'lifecycle', 'getters', 'memory', 'reproducibility', 'fuzz', 'edges_known')
micro = C.hab_common.microsite_individuals


def species(seed, **kw):
    return {'sp1': nemo_species('sp1', 0.2, seed=seed, **kw), 'sp2': nemo_species('sp2', 0.5, seed=seed + 1, **kw)}


def minted(F): return sum(f.minted_num for f in F.values())
def recycled(F): return sum(f.recycled_num for f in F.values())
def pool_n(p, names=POOLS): return sum(len(getattr(h, q)) for h in habs(p) for q in names)
def indi_n(p): return sum(h.indi_num for h in habs(p))


# ------------------------------------------------------------------ one v3-style step (shared by lifecycle / memory / reproducibility)
DISP = {'O->I': 'patch_dispersal_within_from_offspring_pool_to_immigrant_pool', 'OD->I': 'patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool',
        'OD->sites': 'patch_dipersal_within_from_offspring_and_dormancy_pool'}


def v3_step(p, F, variant, dorm, t, rng, disturb_p=0.0, rate=0.5, check=True, tag='', stats=None):
    fs = list(F.values()); stats = stats if stats is not None else {}
    chk = (lambda ph: check_invariants(p, fs, '%s t=%d %s' % (tag, t, ph))) if check else (lambda ph: None)
    r0 = recycled(F)
    d = p.patch_dead_selection(0.1, 0.5, 'niche_gaussian'); assert d == recycled(F) - r0, tag + ' dead != recycled'; chk('dead selection')
    b = p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None); assert b == pool_n(p, ('offspring_pool',)), tag + ' reproduce return'; chk('reproduce')
    m0 = minted(F)
    p.patch_sex_reproduce_calculation_into_offspring_marker_pool(0.6)
    p.patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool(rate, t % 2 == 0)
    assert minted(F) == m0, tag + ' markers minted individuals'
    p.patch_clear_up_offspring_marker_pool(); p.patch_clear_up_immigrant_marker_pool()
    r0, i0, n0, pool0 = recycled(F), pool_n(p, ('immigrant_pool',)), indi_n(p), pool_n(p)
    r = getattr(p, DISP[variant])(rate)
    if variant == 'OD->sites':
        extras = pool0 - pool_n(p) - r
        assert indi_n(p) - n0 == r and recycled(F) - r0 == extras >= 0, tag + ' OD->sites counts'
        stats['extras'] = stats.get('extras', 0) + extras
    else:
        assert pool_n(p, ('immigrant_pool',)) - i0 == r and recycled(F) == r0 and pool_n(p) == pool0, tag + ' dispersal counts'
    stats['migrants'] = stats.get('migrants', 0) + r
    chk('dispersal')
    n0, r0 = indi_n(p), recycled(F)
    if variant == 'OD->sites': g = p.patch_local_germinate_from_offspring_and_dormancy_pool()
    elif dorm: g = p.patch_local_germinate_from_offspring_immigrant_and_dormancy_pool()
    else: g = p.patch_local_germinate_from_offspring_and_immigrant_pool()
    assert g == indi_n(p) - n0 and recycled(F) == r0, tag + ' germinate counts'
    chk('germinate')
    if dorm:
        res = (p.patch_dormancy_process_from_offspring_pool_to_dormancy_pool if variant == 'OD->sites' else p.patch_dormancy_process_from_offspring_pool_and_immigrant_pool)()
        assert res[3] == pool_n(p, ('dormancy_pool',)) and all(len(h.dormancy_pool) <= h.dormancy_pool_max_size for h in habs(p)), tag + ' dormancy'
        chk('dormancy')
    else:
        r0, o, i = recycled(F), pool_n(p, ('offspring_pool',)), pool_n(p, ('immigrant_pool',))
        assert p.patch_clear_up_offspring_pool() == o and p.patch_clear_up_immigrant_pool() == i and recycled(F) - r0 == o + i, tag + ' clear-up counts'
        chk('clear-up')
    if disturb_p and rng.random() < disturb_p:
        dorm_ids = [list(map(id, h.dormancy_pool)) for h in habs(p)]
        r0, exp = recycled(F), indi_n(p) + pool_n(p, ('offspring_pool', 'immigrant_pool'))
        assert p.patch_disturbance_process() == exp == recycled(F) - r0, tag + ' disturbance count'
        assert [list(map(id, h.dormancy_pool)) for h in habs(p)] == dorm_ids, tag + ' disturbance touched dormancy'
        stats['disturbances'] = stats.get('disturbances', 0) + 1
        chk('disturbance')
    if check: factories_balanced(p, fs, tag)
    return stats


# ================================================================ sections
def section_initialize():
    out = {}
    seed_all(1); F = species(1)
    means = [(0.2, 0.2), (0.5, 0.5), (0.2, 0.2)]
    np.random.seed(1); p = make_patch(3, L=6, W=5, K=20, gidx0=10, means=means)
    random.seed(5); np_state = np.random.get_state()[1].copy()
    r = p.patch_initialize(F, TABLE)
    assert r is None
    assert np.array_equal(np.random.get_state()[1], np_state), 'patch_initialize used np.random'
    random.seed(5)
    exp_genders = [random.sample(('male', 'female'), 1)[0] for _ in range(90)]
    assert [x.gender for h in habs(p) for x in micro(h)] == exp_genders, 'gender draw order != habitat-by-habitat original order'
    for h, m in zip(habs(p), means):
        sp = 'sp%d' % (TABLE.index(list(m)) + 1)
        assert h.indi_num == h.size == 30 and h.empty_site_pos_ls == []
        assert all(x.species_id == sp and x.indfactory is F[sp] and x.nemo_indi_obj.getHome() == h.global_index for x in micro(h)), h.name
    nf = sum(g == 'female' for g in exp_genders)
    check_invariants(p, list(F.values()), 'init')
    out['init_female_of_90'] = nf
    # fill the pools, then re-initialize: residents recycled, pools kept, factories balanced
    p.patch_dead_selection(0.3, 0.5, 'niche_gaussian'); p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
    p.patch_dormancy_process_from_offspring_pool_and_immigrant_pool()
    p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
    p.migrants_to_patch_into_habs_immigrant_pool([F['sp2'].make_founder('female', 99) for _ in range(6)])
    residents = [x for h in habs(p) for x in micro(h)]
    pools = {(h.name, q): list(getattr(h, q)) for h in habs(p) for q in POOLS}
    m0, r0 = minted(F), recycled(F)
    p.patch_initialize(F, TABLE)
    assert minted(F) - m0 == 90 and recycled(F) - r0 == len(residents), 're-initialize counts'
    assert not any(x.is_alive for x in residents) and all(x.nemo_indi_obj is None for x in residents)
    assert all(len(getattr(h, q)) == len(pools[(h.name, q)]) and all(a is b for a, b in zip(getattr(h, q), pools[(h.name, q)])) for h in habs(p) for q in POOLS), 'pools touched'
    assert sum(len(v) for v in pools.values()) > 0
    check_invariants(p, list(F.values()), 're-init')
    out['reinit'] = {'residents_recycled': len(residents), 'pool_members_kept': sum(len(v) for v in pools.values())}
    p.patch_initialize(F, TABLE, 'asexual')
    assert all(x.gender == 'female' for h in habs(p) for x in micro(h))
    check_invariants(p, list(F.values()), 'asexual init')
    recycle_everything(p)
    # mean not in the table: first habitat bad -> nothing minted / recycled
    np.random.seed(2); q = make_patch(2, L=3, means=[(0.9, 0.9), (0.2, 0.2)])
    m0, r0 = minted(F), recycled(F)
    try: q.patch_initialize(F, TABLE); res = 'no error'
    except ValueError: res = 'ValueError'
    assert res == 'ValueError' and minted(F) == m0 and recycled(F) == r0 and indi_n(q) == 0
    out['bad_mean_first_habitat'] = 'ValueError, nothing minted'
    for f in F.values(): assert f.minted_num == f.recycled_num
    print('1 patch_initialize: every microsite filled, species per habitat from the table, getHome == global_index, original gender draw order, no np.random; '
          're-initialize recycled %(residents_recycled)d residents, kept %(pool_members_kept)d pool members; asexual -> all female; bad mean (first habitat) -> ValueError, nothing minted' % out['reinit'])
    return out


def section_reproduction():
    out = {}
    seed_all(21)
    sg = nemo_species('sp1', 0.3, seed=21, mutation_rate=0.0, env_variance=0.0)
    F = {'sp1': sg}
    np.random.seed(21); p = make_patch(3, L=8, gidx0=20, means=[(0.3, 0.3)] * 3)
    p.patch_initialize(F, [[0.3, 0.3]])
    err, n_off, nf, homes = 0.0, 0, 0, set()
    for rnd in range(6):
        by_h = {h.name: {x.individual_id: x for x in micro(h)} for h in habs(p)}
        tables = {i: (x.gender, x.species_id, [x.genotype_table(a) for a in range(2)]) for d in by_h.values() for i, x in d.items()}
        m0 = sg.minted_num
        r = p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
        assert r == sum(len(h.offspring_pool) for h in habs(p)) == sg.minted_num - m0, 'reproduce return != offspring minted'
        for h in habs(p):
            for k in h.offspring_pool:
                assert k.mother_id in by_h[h.name] and k.father_id in by_h[h.name], 'parents not from the mating habitat'
                m, f = by_h[h.name][k.mother_id], by_h[h.name][k.father_id]
                assert m.gender == 'female' and f.gender == 'male' and m.species_id == f.species_id == k.species_id and k.indfactory is m.indfactory
                assert k.nemo_indi_obj.getHome() == h.global_index; homes.add(h.global_index)
                for axis in range(2):
                    tab = k.genotype_table(axis)
                    err = max(err, abs(k.get_indi_phenotype_ls()[axis] - sum(a + b for a, b in tab)))
                    for (c0, c1), (m0_, m1), (f0, f1) in zip(tab, tables[k.mother_id][2][axis], tables[k.father_id][2][axis]):
                        assert (c1 in (m0_, m1) and c0 in (f0, f1)) or (c0 in (m0_, m1) and c1 in (f0, f1)), 'allele not inherited'
                n_off += 1; nf += k.gender == 'female'
        check_invariants(p, sg, 'reproduce round %d' % rnd)
        p.patch_dead_selection(0.2, 0.5, 'niche_gaussian')
        p.patch_local_germinate_from_offspring_and_immigrant_pool(); p.patch_clear_up_offspring_pool()
    assert err < 1e-9 and homes == {20, 21, 22}
    # direct birth into microsites
    n_dir, nf_dir = 0, 0
    for rnd in range(6):
        np.random.seed(300 + rnd); p.patch_dead_selection(0.3, 0.5, 'niche_gaussian')
        before = {h.name: {x.individual_id: x for x in micro(h)} for h in habs(p)}
        ids_before = {id(x) for h in habs(p) for x in micro(h)}
        empty, m0 = sum(len(h.empty_site_pos_ls) for h in habs(p)), sg.minted_num
        r = p.patch_sexual_birth_germinate(1.0, None, None)
        assert r == sg.minted_num - m0 and r <= empty, 'direct birth count'
        new = [(h, x) for h in habs(p) for x in micro(h) if id(x) not in ids_before]
        assert len(new) == r
        for h, k in new:
            assert k.mother_id in before[h.name] and k.father_id in before[h.name] and k.nemo_indi_obj.getHome() == h.global_index and k.indfactory is sg
            for axis in range(2): err = max(err, abs(k.get_indi_phenotype_ls()[axis] - sum(a + b for a, b in k.genotype_table(axis))))
            nf_dir += k.gender == 'female'
        n_dir += r
        check_invariants(p, sg, 'direct birth %d' % rnd)
    assert err < 1e-9
    tot, fem = n_off + n_dir, nf + nf_dir
    z = (fem - tot / 2) / math.sqrt(tot / 4)
    assert abs(z) < 4, z
    # fractional rate: per habitat lo or lo+1
    pairs = [h.hab_sexual_pairwise_parents_num() for h in habs(p)]
    for t in range(50):
        np.random.seed(900 + t); r = p.patch_sex_reproduce_mutate_into_offspring_pool(0.37, None, None)
        assert all(int(n * 0.37) <= len(h.offspring_pool) <= int(n * 0.37) + 1 for n, h in zip(pairs, habs(p))) and r == pool_n(p, ('offspring_pool',))
    check_invariants(p, sg, 'fractional')
    recycle_everything(p)
    assert sg.minted_num == sg.recycled_num
    out.update(offspring=n_off, direct_births=n_dir, max_abs_pheno_minus_allele_sum=err, female_z=round(z, 2))
    print('2 reproduction via patch: %d offspring-pool births + %d direct births; parents same species + same habitat, mother\'s factory, getHome == mating habitat (20/21/22), '
          'per-locus inheritance OK, |pheno - sum(alleles)| = %.0e, sex ratio z = %.2f; fractional rate per habitat lo..lo+1' % (n_off, n_dir, err, z))
    return out


def section_two_species():
    seed_all(31); F = species(31)
    np.random.seed(31); p = make_patch(3, L=7, K=25, gidx0=40, means=[(0.2, 0.2), (0.5, 0.5), (0.2, 0.2)])
    p.patch_initialize(F, TABLE)
    mixed, cross, checked, variants = set(), 0, 0, list(DISP)
    for t in range(60):
        p.patch_dead_selection(0.15, 0.5, 'niche_gaussian')
        by_h = {h.name: {x.individual_id: x for x in micro(h)} for h in habs(p)}
        p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
        for h in habs(p):
            for k in h.offspring_pool:
                m, f = by_h[h.name][k.mother_id], by_h[h.name][k.father_id]
                cross += not (m.species_id == f.species_id == k.species_id and k.indfactory is F[k.species_id] and k.nemo_indi_obj.getHome() == h.global_index)
                checked += 1
        getattr(p, DISP[variants[t % 3]])(0.6)
        check_invariants(p, list(F.values()), 'two species t=%d dispersal' % t)
        p.patch_local_germinate_from_offspring_immigrant_and_dormancy_pool()
        p.patch_dormancy_process_from_offspring_pool_and_immigrant_pool()
        for h in habs(p):
            if len({x.species_id for x in micro(h)}) > 1: mixed.add(h.name)
        check_invariants(p, list(F.values()), 'two species t=%d' % t)
    assert cross == 0 and mixed, (cross, mixed)
    recycle_everything(p)
    assert all(f.minted_num == f.recycled_num for f in F.values())
    print('3 two species in different habitats + dispersal: habitats that held both species %s, %d offspring checked, 0 cross-species / cross-habitat pairs, both factories balanced' % (sorted(mixed), checked))
    return {'mixed_habitats': sorted(mixed), 'offspring_checked': checked}


def section_lifecycle():
    out = {}
    n = 0
    for K in (0, 30):
        for variant in DISP:
            for disturb in (False, True):
                seed = 100 + n; seed_all(seed); F = species(seed); rng = random.Random(seed)
                np.random.seed(seed); p = make_patch(3, L=6, K=K, gidx0=50, means=[(0.2, 0.2), (0.5, 0.5), (0.2, 0.2)])
                p.patch_initialize(F, TABLE)
                st = {}
                for t in range(40):
                    v3_step(p, F, variant, K > 0, t, rng, disturb_p=0.15 if disturb else 0.0, rate=0.9, tag='K=%d %s dist=%s' % (K, variant, disturb), stats=st)
                    if disturb and t == 20 and not st.get('disturbances'):
                        p.patch_disturbance_process(); st['disturbances'] = 1
                recycle_everything(p)
                assert all(f.minted_num == f.recycled_num for f in F.values())
                st['minted'] = minted(F)
                out['K=%d %s disturbance=%s' % (K, variant, disturb)] = st
                n += 1
    assert any(v.get('extras') for v in out.values()), 'no OD->sites run produced extra migrants'
    print('4 lifecycle: %d configurations x 40 steps, invariants after every phase, factory balance every step; OD->sites extras recycled %s' %
          (n, {k: v.get('extras') for k, v in out.items() if 'sites' in k}))
    return out


def section_getters():
    seed_all(41); F = species(41)
    np.random.seed(41); p = make_patch(3, L=4, W=5, K=10, gidx0=60, means=[(0.2, 0.2), (0.5, 0.5), (0.2, 0.2)])
    p.patch_initialize(F, TABLE)
    p.patch_dead_selection(0.3, 0.5, 'niche_gaussian')
    sp = p.get_patch_microsites_individals_sp_id_values()
    assert sp.shape == (3, 20)
    for k, h in enumerate(habs(p)):
        for i in range(h.length):
            for j in range(h.width):
                x, v = h.set['microsite_individuals'][i][j], sp[k, i * h.width + j]
                assert (np.isnan(v) if x is None else v == int(x.species_id[2:])), 'sp_id value'
    for a, name in enumerate(('pheno1', 'pheno2')):
        ph = p.get_patch_microsites_individals_phenotype_values(name)
        for k, h in enumerate(habs(p)):
            for i in range(h.length):
                for j in range(h.width):
                    x, v = h.set['microsite_individuals'][i][j], ph[k, i * h.width + j]
                    assert np.isnan(v) if x is None else v == x.phenotype_set[name] == x.get_indi_phenotype_ls()[a], 'phenotype value'
    env = p.get_patch_microsites_environment_values('env2')
    assert np.array_equal(env, np.array([h.set['env2'].flatten() for h in habs(p)]))
    opt = p.get_patch_microsites_optimum_sp_id_value_array(0.1, 0.5, TABLE)
    assert opt.tolist() == [[1] * 20, [2] * 20, [1] * 20]
    d = p.get_data(); assert list(d) == p.hab_id_ls and all(d[h.name] is h.set for h in habs(p))
    assert p.get_patch_individual_num() == indi_n(p) and p.patch_empty_sites_num() == sum(len(h.empty_site_pos_ls) for h in habs(p))
    assert sorted(p.get_patch_empty_sites_ls() + p.get_patch_occupied_sites_ls()) == sorted((h.name, i, j) for h in habs(p) for i in range(4) for j in range(5))
    for (a, b) in p.get_patch_pairwise_occupied_sites_ls():
        x, y = p.set[a[0]].set['microsite_individuals'][a[1]][a[2]], p.set[b[0]].set['microsite_individuals'][b[1]][b[2]]
        assert a[0] == b[0] and x.gender == 'female' and y.gender == 'male' and x.species_id == y.species_id
    for (a, b) in p.get_patch_pairwise_empty_sites_ls():
        assert a[0] == b[0] and a != b and p.set[a[0]].set['microsite_individuals'][a[1]][a[2]] is None
    p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None); p.patch_dormancy_process_from_offspring_pool_to_dormancy_pool(); p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
    op, dp, odp = p.get_patch_offspring_pool(), p.get_patch_dormancy_pool(), p.get_patch_offspring_and_dormancy_pool()
    assert [id(x) for x in op] == [id(x) for h in habs(p) for x in h.offspring_pool] and [id(x) for x in dp] == [id(x) for h in habs(p) for x in h.dormancy_pool]
    assert sorted(map(id, odp)) == sorted(map(id, op + dp))
    assert np.diag(p.get_patch_habs_offspring_num_matrix()).tolist() == [len(h.offspring_pool) for h in habs(p)]
    assert np.diag(p.get_patch_habs_dormancy_num_matrix()).tolist() == [len(h.dormancy_pool) for h in habs(p)]
    assert np.diag(p.get_patch_habs_empty_sites_num_matrix()).tolist() == [len(h.empty_site_pos_ls) for h in habs(p)]
    R = p.get_patch_dispersal_within_rate_matrix(0.3); assert np.allclose(R.sum(axis=1), 1) and np.allclose(np.diag(R), 0.7)
    E = p.get_habs_emigrants_matrix(0.3); assert np.allclose(E.sum(axis=1).A1, [len(h.offspring_pool) + len(h.dormancy_pool) for h in habs(p)])
    I = p.get_habs_immigrants_matrix(0.3); assert np.allclose(np.nansum(I, axis=0).A1, [len(h.empty_site_pos_ls) for h in habs(p)])
    N = p.get_dispersal_within_num_matrix(0.3); assert (N == np.round(N)).all() and (N >= 0).all()
    p.patch_sex_reproduce_calculation_into_offspring_marker_pool(0.8); assert len(p.get_patch_offspring_marker_pool()) == sum(len(h.offspring_marker_pool) for h in habs(p)) > 0
    assert np.diag(p.get_patch_habs_offspring_marker_num_matrix()).tolist() == [len(h.offspring_marker_pool) for h in habs(p)]
    picked = p.sample_offspring_without_replacement(5, ('offspring_pool', 'dormancy_pool'))
    assert len(picked) == 5 and not any(x is y for x in picked for y in held(p))
    p.migrants_to_patch_into_habs_immigrant_pool(picked)
    check_invariants(p, list(F.values()), 'getters')
    # environment offset: env values follow, patch_initialize after an offset finds no species (habitat behaviour, as the original)
    np.random.seed(0); p.patch_offset_environmental_values(['env1', 'env2'], [0.1, 0.1])
    assert np.allclose(p.get_patch_microsites_environment_values('env1').mean(axis=1), [0.3, 0.6, 0.3], atol=0.05)
    try: p.patch_initialize(F, TABLE); r = 'no error'
    except ValueError: r = 'ValueError'
    assert r == 'ValueError'
    recycle_everything(p)
    # regression: one var_env_ls list shared by every habitat (as metaibm/simulator.py passes environment_variation_ls) + delta_var_ls
    shared_var = [0.025, 0.025]
    q = patch('q', 0, (0, 0))
    for k in range(3): q.add_habitat('h%d' % k, k, (0, k), 2, ['env1', 'env2'], [0.2, 0.2], shared_var, 2, 2)
    np.random.seed(0); q.patch_offset_environmental_values(['env1'], [0.0], [0.01])
    vars_after = [round(h.var_env_ls[0], 4) for h in habs(q)]
    assert all(abs(h.var_env_ls[0] - 0.035) < 1e-12 for h in habs(q)), 'shared var list: delta applied more than once %s' % vars_after
    assert len({id(h.var_env_ls) for h in habs(q)} | {id(shared_var)}) == 4 and shared_var == [0.025, 0.025], 'habitats share var_env_ls (or wrote the caller\'s list)'
    assert all(f.minted_num == f.recycled_num for f in F.values())
    print('5 getters with real individuals: sp_id, phenotype (phenotype_set == get_indi_phenotype_ls), environment, optimum sp_id, get_data identity, matrices, pools, pairwise lists OK; offset -> re-initialize ValueError; shared var list -> own copies, var %s' % vars_after)
    return {'shared_var_after_offset': vars_after}


def section_memory():
    out = {}
    # (a) reference individuals in an untouched patch stay unchanged while another patch churns
    seed_all(51); F = species(51)
    np.random.seed(51); ref = make_patch(1, L=4, gidx0=1, name='ref'); ref.patch_initialize(F, TABLE)
    state = lambda x: (x.gender, x.get_indi_phenotype_ls(), [x.genotype_table(a) for a in range(2)], x.nemo_indi_obj.getHome(), x.mother_id, x.father_id)
    snap = {x.individual_id: state(x) for x in held(ref)}
    np.random.seed(52); p = make_patch(3, L=12, K=50, gidx0=70, means=[(0.2, 0.2), (0.5, 0.5), (0.2, 0.2)])
    p.patch_initialize(F, TABLE)
    seen = {f.species_id: {int(x.nemo_indi_obj.this) for x in held([p, ref]) if x.indfactory is f} for f in F.values()}   # per factory: each IndFactory has its own recycling pool
    reuse_late, births_late, distinct_at = 0, 0, {}
    rss, rng, dead_sample = {}, random.Random(52), []
    hw = {f.species_id: f.minted_num - f.recycled_num for f in F.values()}
    STEPS = 320
    for t in range(STEPS):
        live_before = set(map(id, held(p)))
        p.patch_dead_selection(0.2, 0.5, 'niche_gaussian')
        p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
        born = [k for h in habs(p) for k in h.offspring_pool]
        addrs = [(k.species_id, int(k.nemo_indi_obj.this)) for k in born]
        live = held([p, ref])
        for f in F.values(): hw[f.species_id] = max(hw[f.species_id], f.minted_num - f.recycled_num)   # alive high-water mark per factory (peak is right after reproduction)
        assert len({int(x.nemo_indi_obj.this) for x in live}) == len(live), 'two live individuals share one C++ Individual at t=%d' % t
        if t >= STEPS // 2:
            reuse_late += sum(a in seen[s] for s, a in addrs); births_late += len(addrs)
        for s, a in addrs: seen[s].add(a)
        if t in (50, 100, 160, 240, STEPS - 1): distinct_at[t] = {s: (len(seen[s]), hw[s]) for s in seen}
        p.patch_dipersal_within_from_offspring_and_dormancy_pool(0.4)
        p.patch_local_germinate_from_offspring_and_dormancy_pool()
        p.patch_dormancy_process_from_offspring_pool_to_dormancy_pool()
        if t == 100: p.patch_disturbance_process()                  # one disturbance only: frequent ones consume the dormancy bank and drive the patch extinct (model dynamics)
        if t % 10 == 0: check_invariants([p, ref], list(F.values()), 'memory t=%d' % t)
        dead_sample += [k for k in born if not k.is_alive][:2]
        if t in (50, 160, STEPS - 1): rss[t] = round(rss_mb(), 1)
    assert all(state(x) == snap[x.individual_id] for x in held(ref)), 'a living individual changed while memory was recycled / reused'
    for x in dead_sample[:40]:
        assert x.nemo_indi_obj is None and not x.is_alive
        for fn in (lambda: x.gender, lambda: x.get_indi_phenotype_ls(), lambda: x.individual_id, lambda: x.phenotype_set, lambda: x.recycle()):
            try: fn(); raise AssertionError('recycled individual still usable')
            except RuntimeError: pass
    reuse = reuse_late / max(births_late, 1)
    tot = lambda t: sum(v[0] for v in distinct_at[t].values())
    growth = (tot(STEPS - 1) - tot(160)) / tot(160)
    out['address_reuse'] = {'late_births': births_late, 'late_births_on_recycled_address': round(reuse, 4), 'distinct_addresses_vs_alive_high_water_mark': distinct_at, 'growth_after_t160': round(growth, 4)}
    out['rss_mb'] = rss; out['births'] = minted(F)
    assert births_late > 5000 and indi_n(p) > 0, 'no steady state in the second half: %d births, %d alive' % (births_late, indi_n(p))
    assert reuse > 0.95, out['address_reuse']
    assert all(n <= h for d in distinct_at.values() for n, h in d.values()), out['address_reuse']   # an IndFactory allocates only when its recycling pool is empty: distinct addresses <= its alive high-water mark
    assert tot(STEPS - 1) - tot(160) < 0.01 * births_late, out['address_reuse']                       # growth (if any) tracks new population peaks, not births
    assert rss[STEPS - 1] - rss[50] < 40, rss
    recycle_everything([p, ref])
    assert all(f.minted_num == f.recycled_num for f in F.values()), 'teardown: minted != recycled'
    out['teardown'] = {f.species_id: (f.minted_num, f.recycled_num) for f in F.values()}
    for f in F.values(): f.purge_recycling_pool()
    print('6 memory: %d steps, %d births; RSS %s MB; late births on a recycled address %.1f%%, per factory (distinct addresses, alive high-water mark) %s (growth after t=160: %.1f%%); reference individuals unchanged; '
          'recycled individuals inert; teardown minted == recycled %s' % (STEPS, out['births'], rss, 100 * reuse, distinct_at, 100 * growth, out['teardown']))
    return out


def trace_run(seed, steps=50):
    seed_all(seed); F = species(seed); rng = random.Random(seed)
    np.random.seed(seed); p = make_patch(3, L=8, K=20, gidx0=80, means=[(0.2, 0.2), (0.5, 0.5), (0.2, 0.2)])
    p.patch_initialize(F, TABLE)
    trace = []
    for t in range(steps):
        v3_step(p, F, list(DISP)[t % 3], True, t, rng, disturb_p=0.1, check=False)
        trace.append(nemo_digest(p))
    return trace


def section_reproducibility():
    out = {}
    runs = {}
    for seed in (7, 7, 8):
        pr = subprocess.run([sys.executable, os.path.join(HERE, 'test_nemo.py'), '--trace', str(seed)], capture_output=True, text=True, cwd=HERE)
        line = [l for l in pr.stdout.splitlines() if l.startswith('TRACE ')]
        assert line, pr.stdout[-500:] + pr.stderr[-2000:]
        runs.setdefault(seed, []).append(json.loads(line[-1][6:]))
    assert runs[7][0] == runs[7][1] and len(runs[7][0]) == 50, 'not reproducible across fresh processes'
    assert runs[8][0] != runs[7][0], 'different seed gave the same trace'
    out['seed7'] = 'identical nemo_digest over %d steps in two fresh processes' % len(runs[7][0]); out['seed8_differs'] = True
    print('7 reproducibility:', out['seed7'], '; seed 8 differs')
    return out


def section_fuzz():
    out = {}
    seed_all(81); F = species(81); fs = list(F.values())
    n_ops, counts = 0, {}
    for seed in range(30):
        rng = random.Random(seed); np.random.seed(seed); random.seed(seed)
        means = lambda n: [rng.choice(((0.2, 0.2), (0.5, 0.5))) for _ in range(n)]
        n1, n2 = rng.choice((1, 2, 3, 4)), rng.choice((2, 3))
        p = make_patch(n1, L=rng.choice((1, 2, 4)), W=rng.choice((1, 3)), K=rng.choice((0, 1, 10)), gidx0=seed * 10, means=means(n1), name='p')
        q = make_patch(n2, L=3, K=rng.choice((0, 5)), gidx0=500 + seed * 10, means=means(n2), name='q')
        p.patch_initialize(F, TABLE); q.patch_initialize(F, TABLE)

        def kill_one():
            occ = p.get_patch_occupied_sites_ls()
            if occ:
                h_id, i, j = rng.choice(occ); h = p.set[h_id]; x = h.set['microsite_individuals'][i][j]; h.del_individual(i, j); x.recycle()
        OPS = {'dead': lambda: p.patch_dead_selection(rng.choice((0.0, 0.1, 0.5, 1.0)), rng.choice((0.2, 0.5, 2.0)), rng.choice(('niche_gaussian', 'neutral_uniform'))),
               'reproduce': lambda: p.patch_sex_reproduce_mutate_into_offspring_pool(rng.choice((0.0, 0.3, 0.77, 1.0)), None, None),
               'birth_germinate': lambda: p.patch_sexual_birth_germinate(rng.choice((0.5, 1.0)), None, None),
               'marker_birth': lambda: p.patch_sex_reproduce_calculation_into_offspring_marker_pool(rng.random()),
               'marker_disp': lambda: p.patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool(rng.random(), rng.random() < 0.5),
               'O->I': lambda: p.patch_dispersal_within_from_offspring_pool_to_immigrant_pool(rng.choice((0.0, 0.3, 1.0))),
               'OD->I': lambda: p.patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool(rng.choice((0.0, 0.3, 1.0))),
               'OD->sites': lambda: p.patch_dipersal_within_from_offspring_and_dormancy_pool(rng.choice((0.0, 0.3, 0.8, 1.0))),
               'germ_OD': p.patch_local_germinate_from_offspring_and_dormancy_pool, 'germ_OI': p.patch_local_germinate_from_offspring_and_immigrant_pool,
               'germ_OID': p.patch_local_germinate_from_offspring_immigrant_and_dormancy_pool,
               'dormA': p.patch_dormancy_process_from_offspring_pool_to_dormancy_pool, 'dormB': p.patch_dormancy_process_from_offspring_pool_and_immigrant_pool,
               'clear_micro': p.patch_clear_up_microsite_individuals, 'clear_O': p.patch_clear_up_offspring_pool, 'clear_I': p.patch_clear_up_immigrant_pool,
               'clear_D': p.patch_clear_up_dormancy_pool, 'clear_OM': p.patch_clear_up_offspring_marker_pool, 'clear_IM': p.patch_clear_up_immigrant_marker_pool,
               'disturbance': p.patch_disturbance_process, 'initialize': lambda: p.patch_initialize(F, TABLE, rng.choice(('sexual', 'asexual'))),
               'immigrants_in': lambda: p.migrants_to_patch_into_habs_immigrant_pool([rng.choice(fs).make_founder(rng.choice(('male', 'female')), rng.randint(0, 9)) for _ in range(rng.randint(0, 6))]),
               'emigrate_to_q': lambda: q.migrants_to_patch_into_habs_immigrant_pool(p.sample_offspring_without_replacement(rng.randint(0, 5), rng.choice((('offspring_pool',), ('offspring_pool', 'dormancy_pool'))))),
               'q_step': lambda: (q.patch_dead_selection(0.2, 0.5, 'niche_gaussian'), q.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None),
                                  q.patch_local_germinate_from_offspring_immigrant_and_dormancy_pool(), q.patch_dormancy_process_from_offspring_pool_and_immigrant_pool()),
               'kill_one': kill_one}
        names = list(OPS)
        for step in range(110):
            k = rng.choice(names); OPS[k]()
            check_invariants([p, q], fs, 'fuzz seed=%d step=%d op=%s' % (seed, step, k))
            assert all(len(h.dormancy_pool) <= h.dormancy_pool_max_size for h in habs(p) + habs(q))
            counts[k] = counts.get(k, 0) + 1; n_ops += 1
        recycle_everything([p, q])
        assert all(f.minted_num == f.recycled_num for f in fs)
    out['operations'] = n_ops; out['per_op'] = counts
    assert n_ops >= 3000 and len(counts) == len(names)
    print('8 fuzz: %d random patch-level operations (%d kinds) on two patches, invariants after each, minted == recycled after every run' % (n_ops, len(counts)))
    return out


ASEX_MIX = [('patch_asexual_birth_germinate', lambda p: p.patch_asexual_birth_germinate(0.5, 0.01, [0.01, 0.01])),
            ('patch_mixed_birth_germinate', lambda p: p.patch_mixed_birth_germinate(0.5, 0.5, 0.01, [0.01, 0.01])),
            ('patch_asex_reproduce_calculation_into_offspring_marker_pool', lambda p: p.patch_asex_reproduce_calculation_into_offspring_marker_pool(0.5)),
            ('patch_mix_reproduce_calculation_into_offspring_marker_pool', lambda p: p.patch_mix_reproduce_calculation_into_offspring_marker_pool(0.5, 0.5)),
            ('patch_asex_reproduce_mutate_into_offspring_pool', lambda p: p.patch_asex_reproduce_mutate_into_offspring_pool(0.5, 0.01, [0.01, 0.01])),
            ('patch_mix_reproduce_mutate_into_offspring_pool', lambda p: p.patch_mix_reproduce_mutate_into_offspring_pool(0.5, 0.5, 0.01, [0.01, 0.01]))]


def probe(i):
    seed_all(3); F = species(3)
    np.random.seed(3); p = make_patch(2, L=6, means=[(0.2, 0.2), (0.5, 0.5)]); p.patch_initialize(F, TABLE)
    for h in habs(p):
        for k in range(10): x = h.set['microsite_individuals'][k // 6][k % 6]; h.del_individual(k // 6, k % 6); x.recycle()
    for _ in range(3): p.patch_dead_selection(0.0, 5.0, 'niche_gaussian')
    p.patch_dead_selection(0.0, 0.5, 'niche_gaussian')
    name, fn = ASEX_MIX[i]
    try: r = fn(p); res = 'returned %r' % (r,)
    except Exception as e: res = '%s: %s' % (type(e).__name__, str(e)[:70])
    return [name, res]


def section_edges_known():
    out = {}
    seed_all(61); F = species(61); fs = list(F.values())
    # (1) patch_initialize is not atomic across habitats: a bad later habitat raises after earlier habitats were re-initialized
    res = {}
    for label, means, kw, factories in (('mean not in table (h2)', [(0.2, 0.2), (0.5, 0.5), (0.9, 0.9)], {}, F),
                                        ('missing indfactory key (h1 -> sp2)', [(0.2, 0.2), (0.5, 0.5), (0.2, 0.2)], {}, {'sp1': F['sp1']}),
                                        ('global_index 65536 (h1)', [(0.2, 0.2)] * 3, {'gi': [0, 65536, 2]}, F)):
        np.random.seed(1); p = patch('p', 0, (0, 0))
        gis = kw.get('gi', [0, 1, 2])
        for k in range(3): p.add_habitat('h%d' % k, k, (0, k), 2, ['env1', 'env2'], list(means[k]), [0.025, 0.025], 3, 3, 0, global_index=gis[k])
        p.set['h0'].hab_initialize(F, TABLE)                       # residents from an earlier initialization
        if label.startswith('mean'): p.set['h1'].hab_initialize(F, TABLE)
        else:
            p.set['h1'].global_index = 1; p.set['h1'].hab_initialize(F, TABLE); p.set['h1'].global_index = gis[1]
        before = {h.name: list(micro(h)) for h in habs(p)}                # keep references: id() of collected wrappers is reused
        m0, r0 = minted(F), recycled(F)
        try: p.patch_initialize(factories, TABLE); err = 'no error'
        except Exception as e: err = type(e).__name__
        after = {h.name: (h.indi_num, sum(any(x is y for y in before[h.name]) for x in micro(h))) for h in habs(p)}
        check_invariants(p, fs, label)
        res[label] = {'error': err, 'minted': minted(F) - m0, 'recycled': recycled(F) - r0, 'per_habitat (indi_num, old residents kept)': after}
        recycle_everything(p)
    out['initialize_atomicity'] = res
    for label, r in res.items(): print('9.1 %-36s %s' % (label, r))
    assert res['mean not in table (h2)']['error'] == 'ValueError' and res['mean not in table (h2)']['minted'] == 18
    assert res['global_index 65536 (h1)']['error'] == 'TypeError' and res['global_index 65536 (h1)']['per_habitat (indi_num, old residents kept)']['h1'] == (0, 0)
    finding('risk', 'patch_initialize is not atomic: it calls hab_initialize habitat by habitat, so a bad habitat late in the order raises after earlier habitats were already '
            're-initialized (residents recycled, new founders minted); with an out-of-range global_index the failing habitat itself is left EMPTY (hab_initialize recycles residents before the first make_founder raises). Observed: %s. Invariants and factory balance still hold.' % {k: (v['error'], 'minted %d' % v['minted'], 'recycled %d' % v['recycled']) for k, v in res.items()},
            repro="p with h0 (0.2,0.2), h1 (0.5,0.5), h2 (0.9,0.9) all occupied; p.patch_initialize(F, [[0.2,0.2],[0.5,0.5]]) -> ValueError after h0/h1 were refilled (18 minted, 18 recycled)",
            fix='validate first: for every habitat resolve species_2_phenotype_ls.index(mean_env_ls), check the indfactory key and 0 <= global_index <= 65535 (and reproduce_mode), then initialize; or validate global_index in habitat.__init__ / add_habitat')
    # (2) global_index bounds through add_habitat
    gi_res = {}
    for gi in (0, 65535, 65536, -1):
        np.random.seed(2); p = make_patch(1, L=2, gidx0=gi)
        m0 = minted(F)
        try: p.patch_initialize(F, TABLE); gi_res[gi] = 'ok, getHome()=%d' % micro(p.set['h0'])[0].nemo_indi_obj.getHome()
        except Exception as e: gi_res[gi] = '%s, minted %d' % (type(e).__name__, minted(F) - m0)
        check_invariants(p, fs, 'gi'); recycle_everything(p)
    assert gi_res[0].startswith('ok') and gi_res[65535] == 'ok, getHome()=65535' and gi_res[65536].startswith('TypeError, minted 0') and gi_res[-1].startswith('TypeError, minted 0')
    finding('known', 'add_habitat(global_index=...) is not validated; out of the unsigned-short range Nemo refuses it at the first make_founder: %s.' % gi_res)
    # (3) other inputs
    np.random.seed(3); p = make_patch(2, L=2)
    errs = {}
    for label, call in (('tuple mean_env_ls in table', lambda: make_patch(1, L=2, means=[(0.2, 0.2)]).patch_initialize(F, [(0.2, 0.2)])),
                        ("reproduce_mode='Sexual'", lambda: p.patch_initialize(F, TABLE, 'Sexual')),
                        ('sexual_birth_rate 1.5', lambda: (p.patch_initialize(F, TABLE), p.patch_sex_reproduce_mutate_into_offspring_pool(1.5, None, None)))):
        m0 = minted(F)
        try: call(); errs[label] = 'no error'
        except Exception as e: errs[label] = type(e).__name__
        check_invariants(p, fs, label)
    recycle_everything(p)
    out['inputs'] = errs
    finding('known', 'sexual_birth_rate > 1 -> %s from random.sample at habitat level (one offspring per pair per call); invariants hold afterwards.' % errs['sexual_birth_rate 1.5'])
    finding('info', 'patch_initialize finds species by exact equality species_2_phenotype_ls.index(mean_env_ls): a tuple row vs list mean -> %s; reproduce_mode="Sexual" -> %s (gender unbound in hab_initialize, after residents were recycled); no validation (habitat-level, as in the habitat report).' % (errs['tuple mean_env_ls in table'], errs["reproduce_mode='Sexual'"]))
    finding('info', 'patch_initialize now returns None (the original returned 0); metacommunity.meta_initialize ignores the return value.')
    # (4) asex / mix probes (not migrated), each in its own process (a SWIG deepcopy could crash it)
    probes = []
    for i in range(len(ASEX_MIX)):
        pr = subprocess.run([sys.executable, os.path.join(HERE, 'test_nemo.py'), '--probe', str(i)], capture_output=True, text=True, cwd=HERE)
        line = [l for l in pr.stdout.splitlines() if l.startswith('PROBE ')]
        probes.append(json.loads(line[-1][6:]) if line else [ASEX_MIX[i][0], 'process died rc=%d' % pr.returncode])
    out['asex_mix_probes'] = probes
    for n, r in probes:
        if 'calculation' not in n: assert r.startswith('TypeError'), (n, r)
    finding('known', 'asex / mix patch functions are not migrated: ' + '; '.join('%s -> %s' % (n, r.split(':')[0] if r.startswith('TypeError') else r) for n, r in probes) + ' (TypeError = cannot pickle SwigPyObject in copy.deepcopy; the two marker calculations create no individual).')
    print('9.2 global_index', gi_res, '| inputs', errs)
    print('9.4 asex / mix probes:', ['%s: %s' % (n[6:], r.split(':')[0]) for n, r in probes])
    assert all(f.minted_num == f.recycled_num for f in fs)
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


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == '--trace':
        print('TRACE ' + json.dumps(trace_run(int(sys.argv[2])))); sys.exit(0)
    if len(sys.argv) > 2 and sys.argv[1] == '--probe':
        print('PROBE ' + json.dumps(probe(int(sys.argv[2])))); sys.exit(0)
    if len(sys.argv) > 2 and sys.argv[1] == '--section':
        sys.exit(0 if run_section(sys.argv[2]) else 1)
    wanted = sys.argv[1:] or list(SECTIONS)
    failed, results, findings = [], {}, []
    for name in wanted:
        print('==== section %s' % name, flush=True)
        pr = subprocess.run([sys.executable, os.path.abspath(__file__), '--section', name], capture_output=True, text=True, cwd=HERE)
        for l in (pr.stdout + pr.stderr).splitlines():
            if l.strip() and not l.startswith(('setting random seed', '    outputs: ')): print(l)
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
