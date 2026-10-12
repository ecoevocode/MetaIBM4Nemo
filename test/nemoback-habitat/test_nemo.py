# All tests that use real Nemo individuals, one section each (sections run in this process and report independently):
#   reproduction  founders, death probability == survival_rate, parents + per-locus inheritance, mutation, sex ratio, count rule, markers, two species, degenerate cases
#   lifecycle     v3 loops + dormancy recovery, cross-process reproducibility (fresh subprocesses), selection response, multi-habitat, random-order fuzz
#   memory        living individuals unchanged under recycle / reuse, unique pointers + IDs, recycled individuals inert, long-run RSS
#   edges_known   degenerate shapes, constructor / input guards, legacy helpers, asex / mix probes (subprocesses), known behaviour (recorded only)
import os, sys, json, math, copy, random, hashlib, subprocess, traceback
import numpy as np
from common import (habitat, fake_indi, make_species, make_hab, check_invariants, held_individuals, microsite_individuals, recycle_everything,
                    save, seed_all, finding, chi2_ok, rss_mb, HERE)

ASEX_MIX = [('hab_asex_reproduce_mutate_with_num', lambda h: h.hab_asex_reproduce_mutate_with_num(0.01, [0.01, 0.01], 2)),
            ('hab_asex_reproduce_mutate_into_offspring_pool', lambda h: h.hab_asex_reproduce_mutate_into_offspring_pool(0.5, 0.01, [0.01, 0.01])),
            ('hab_asexual_reprodece_germinate', lambda h: h.hab_asexual_reprodece_germinate(0.5, 0.01, [0.01, 0.01])),
            ('hab_mix_asex_reproduce_mutate_with_num', lambda h: h.hab_mix_asex_reproduce_mutate_with_num(0.01, [0.01, 0.01], 1)),
            ('hab_mix_sex_reproduce_mutate_with_num', lambda h: h.hab_mix_sex_reproduce_mutate_with_num(0.01, [0.01, 0.01], 1)),
            ('hab_mix_reproduce_mutate_into_offspring_pool', lambda h: h.hab_mix_reproduce_mutate_into_offspring_pool(0.5, 0.5, 0.01, [0.01, 0.01])),
            ('hab_mixed_reproduce_germinate', lambda h: h.hab_mixed_reproduce_germinate(0.5, 0.5, 0.01, [0.01, 0.01])),
            ('hab_asex_reproduce_calculation_into_offspring_marker_pool', lambda h: h.hab_asex_reproduce_calculation_into_offspring_marker_pool('p', 0.5)),
            ('hab_mix_reproduce_calculation_into_offspring_marker_pool', lambda h: h.hab_mix_reproduce_calculation_into_offspring_marker_pool('p', 0.5, 0.5))]
# ================================================================ 1. v3 loops
def v3_run(L, seed, mode, steps=50, check=True):
    seed_all(seed)
    sp = make_species('sp1', 0.2, seed=seed)
    K = 0 if mode == 'individual' else 40
    h = make_hab(L=L, K=K, global_index=4)
    h.hab_initialize({'sp1': sp}, [[0.2, 0.2]])
    trace, st, pending = [], {'disturbances': 0, 'recoveries': 0}, False
    for t in range(steps):
        tag = '%s L=%d seed=%d t=%d' % (mode, L, seed, t)
        r0 = sp.recycled_num
        assert h.hab_dead_selection(0.1, 0.5, 'niche_gaussian') == sp.recycled_num - r0, tag + ' dead != recycled'
        h.hab_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
        if mode != 'individual':
            for _ in range(random.randint(0, 15)): h.immigrant_pool.append(sp.make_founder(random.choice(('male', 'female')), 9))
        pools = ('offspring_pool', 'immigrant_pool') + (('dormancy_pool',) if mode != 'individual' else ())
        empty, cand = len(h.empty_site_pos_ls), sum(len(getattr(h, p)) for p in pools)
        dorm_ids = set(map(id, h.dormancy_pool))
        g = (h.hab_local_germinate_from_offspring_and_immigrant_pool() if mode == 'individual' else h.hab_local_germinate_from_offspring_immigrant_and_dormancy_pool())
        assert g == min(empty, cand), tag + ' germinate'
        if pending:
            assert any(id(x) in dorm_ids for x in microsite_individuals(h)), tag + ' no recovery from dormancy'
            st['recoveries'] += 1; pending = False
        if check: check_invariants(h, sp, tag + ' germinate')
        if mode == 'individual':
            h.hab_clear_up_offspring_pool(); h.hab_clear_up_immigrant_pool()
        else:
            r = (h.hab_dormancy_process_from_offspring_pool_to_dormancy_pool if mode == 'dormA' else h.hab_dormancy_process_from_offspring_pool_and_immigrant_pool)()
            assert len(h.dormancy_pool) <= K and r[3] == len(h.dormancy_pool)
        if random.random() < (0.02 if mode == 'individual' else 0.05) or (mode != 'individual' and t == steps // 2):
            h.habitat_disturbance_process(); st['disturbances'] += 1
            pending = mode != 'individual' and len(h.dormancy_pool) > 0
        if check: check_invariants(h, sp, tag + ' end of step')
        ids = sorted(x.individual_id for x in held_individuals(h))
        trace.append([h.indi_num, hashlib.sha1(repr(ids).encode()).hexdigest()[:16], round(sum(sum(x.get_indi_phenotype_ls()) for x in microsite_individuals(h)), 10)])
    recycle_everything(h)
    assert sp.minted_num == sp.recycled_num
    st['minted'] = sp.minted_num
    sp.purge_recycling_pool()
    return trace, st


# ---- child-process entry points: --trace (reproducibility) and --probe (asex / mix, a SWIG deepcopy could crash the process)
if len(sys.argv) > 1 and sys.argv[1] == '--trace':
    print('TRACE ' + json.dumps(v3_run(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], check=False)[0])); sys.exit(0)
if len(sys.argv) > 2 and sys.argv[1] == '--probe':
    seed_all(3); s = make_species('sp1', 0.2, seed=3); h = make_hab(L=6); h.hab_initialize({'sp1': s}, [[0.2, 0.2]])
    for k in range(10): x = h.set['microsite_individuals'][k // 6][k % 6]; h.del_individual(k // 6, k % 6); x.recycle()
    for _ in range(3): h.hab_dead_selection(0.0, 5.0, 'niche_gaussian')
    h.hab_dead_selection(0.0, 0.5, 'niche_gaussian')
    name, fn = ASEX_MIX[int(sys.argv[2])]
    try: r = fn(h); res = 'returned %r' % (r,)
    except Exception as e: res = '%s: %s' % (type(e).__name__, str(e)[:70])
    print('PROBE ' + json.dumps([name, res])); sys.exit(0)


def section_reproduction():
    out = {}
    seed_all(11)
    sp1, sp2 = make_species('sp1', 0.5, seed=11), make_species('sp2', 0.2, seed=12)
    sg = make_species('sp_g', 0.3, seed=13, mutation_rate=0.0, env_variance=0.0)       # genetics: no mutation, Ve = 0
    sm = make_species('sp_m', 0.3, seed=14, mutation_rate=0.01, env_variance=0.0)
    FACT = [sp1, sp2, sg, sm]
    TABLE = [[0.5, 0.5], [0.2, 0.2]]


    # ================================================================ 1. hab_initialize
    for gi in (0, 7):
        h = make_hab(L=12, global_index=gi, mean=(0.2, 0.2))
        random.seed(gi); np_state = np.random.get_state()[1].copy()
        h.hab_initialize({'sp1': sp1, 'sp2': sp2}, TABLE, 'sexual')
        ms = microsite_individuals(h)
        nf = sum(x.gender == 'female' for x in ms)
        assert len(ms) == 144 and abs(nf - 72) <= 0.3 * 72, nf
        assert all(x.species_id == 'sp2' and x.indfactory is sp2 and x.nemo_indi_obj.getHome() == gi for x in ms)
        random.seed(gi)
        assert [x.gender for x in ms] == [random.sample(('male', 'female'), 1)[0] for _ in range(144)], 'gender draw order != original'
        assert np.array_equal(np.random.get_state()[1], np_state), 'hab_initialize used np.random'
        check_invariants(h, FACT, 'init gi=%d' % gi)
        out['init_female_of_144_gi%d' % gi] = nf
        recycle_everything(h)
        h = make_hab(L=12, global_index=gi, mean=(0.5, 0.5))
        h.hab_initialize({'sp1': sp1, 'sp2': sp2}, TABLE, 'asexual')
        assert all(x.gender == 'female' and x.species_id == 'sp1' for x in microsite_individuals(h))
        recycle_everything(h)
    h = make_hab(L=3)
    h.hab_initialize({'sp2': sp2}, TABLE)
    x = h.set['microsite_individuals'][1][1]; h.del_individual(1, 1); x.recycle()
    m0, r0, residents = sp2.minted_num, sp2.recycled_num, [x for x in microsite_individuals(h)]
    h.hab_initialize({'sp2': sp2}, TABLE)                                   # reset: residents recycled, every site gets a new founder
    assert sp2.minted_num - m0 == 9 and sp2.recycled_num - r0 == 8 and h.indi_num == 9
    assert not any(r.is_alive for r in residents) and not any(r is y for r in residents for y in microsite_individuals(h))
    check_invariants(h, FACT, 'init on occupied')
    m0 = sp2.minted_num; recycle_everything(h)
    try: make_hab(L=3).hab_initialize({'sp2': sp2}, TABLE, 'Sexual'); r = 'no error'
    except Exception as e: r = type(e).__name__
    assert r != 'no error' and sp2.minted_num == m0, 'invalid reproduce_mode accepted'
    out['invalid_reproduce_mode'] = r
    print('1 hab_initialize: sex ratio %s, home == global_index, species from table, original draw order, occupied sites reset, invalid mode -> %s' % (out['init_female_of_144_gi0'], r))

    # ================================================================ 2. death probability == survival_rate
    for d, w, method, env in ((0.1, 0.5, 'niche_gaussian', (0.2, 0.2)), (0.0, 0.15, 'niche_gaussian', (0.35, 0.3)), (0.4, 0.5, 'neutral_uniform', (0.2, 0.2))):
        exp_s = var_s = obs_s = 0.0
        for t in range(60):
            np.random.seed(500 + t)
            h = make_hab(L=10, mean=env)
            h.hab_initialize({'sp1': sp2}, [list(env)])
            ps = [h.survival_rate(d, h.set['microsite_individuals'][i][j].get_indi_phenotype_ls(), h.get_microsite_env_val_ls(i, j), w, method) for i in range(10) for j in range(10)]
            r0 = sp2.recycled_num
            dead = h.hab_dead_selection(d, w, method)
            assert dead == sp2.recycled_num - r0
            exp_s += sum(ps); var_s += sum(p * (1 - p) for p in ps); obs_s += 100 - dead
            recycle_everything(h)
        z = (obs_s - exp_s) / math.sqrt(var_s)
        out['death_%s_d%g_w%g_z' % (method, d, w)] = round(z, 2)
        assert abs(z) < 4, z
    print('2 death probability vs survival_rate: z =', [v for k, v in out.items() if k.startswith('death_')])

    # ================================================================ 3. reproduction: parents, genetics, mutation, sex ratio, parent sampling
    h = make_hab(L=12, global_index=5, mean=(0.3, 0.3))
    h.hab_initialize({'sp1': sg}, [[0.3, 0.3]])
    by_id = {x.individual_id: x for x in microsite_individuals(h)}
    err = 0.0
    for num in (0, 1, 10, h.hab_sexual_pairwise_parents_num()):
        kids = h.hab_sex_reproduce_mutate_with_num(None, None, num)
        assert len(kids) == num
        for k in kids:
            m, f = by_id[k.mother_id], by_id[k.father_id]
            assert k.nemo_indi_obj.getHome() == 5 and m.gender == 'female' and f.gender == 'male' and k.indfactory is m.indfactory
            for axis in range(2):
                tab = k.genotype_table(axis)
                err = max(err, abs(k.get_indi_phenotype_ls()[axis] - sum(a + b for a, b in tab)))
                for (c0, c1), (m0_, m1), (f0, f1) in zip(tab, m.genotype_table(axis), f.genotype_table(axis)):
                    assert (c1 in (m0_, m1) and c0 in (f0, f1)) or (c0 in (m0_, m1) and c1 in (f0, f1)), 'allele not inherited'
            k.recycle()
    assert err < 1e-9
    out['max_abs_phenotype_minus_allele_sum'] = err
    females = [i for i, x in by_id.items() if x.gender == 'female']; males = [i for i, x in by_id.items() if x.gender == 'male']
    nf = tot = 0; mu, fu = dict.fromkeys(females, 0), dict.fromkeys(males, 0)
    T, NUM, pairs = 300, 8, h.hab_sexual_pairwise_parents_num()
    for t in range(T):
        kids = h.hab_sex_reproduce_mutate_with_num(None, None, NUM)
        assert len({k.mother_id for k in kids}) == NUM and len({k.father_id for k in kids}) == NUM, 'a parent bred twice in one batch'
        for k in kids:
            nf += k.gender == 'female'; tot += 1; mu[k.mother_id] += 1; fu[k.father_id] += 1; k.recycle()
    z = (nf - tot / 2) / math.sqrt(tot / 4)
    okm, sm_ = chi2_ok(list(mu.values()), [T * NUM / max(len(females), pairs)] * len(females))
    okf, sf_ = chi2_ok(list(fu.values()), [T * NUM / max(len(males), pairs)] * len(males))
    out.update(offspring_female_z=round(z, 2), mother_chi2=sm_, father_chi2=sf_)
    assert abs(z) < 4 and okm and okf, (z, sm_, sf_)
    recycle_everything(h)
    h = make_hab(L=8, mean=(0.3, 0.3))
    h.hab_initialize({'sp1': sm}, [[0.3, 0.3]])
    by_id = {x.individual_id: x for x in microsite_individuals(h)}
    h.hab_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
    new = copies = 0
    for k in h.offspring_pool:
        for axis in range(2):
            par = {a for row in by_id[k.mother_id].genotype_table(axis) + by_id[k.father_id].genotype_table(axis) for a in row}
            for a, b in k.genotype_table(axis): copies += 2; new += (a not in par) + (b not in par)
    out['mutation_0.01_new_alleles'] = [new, copies]
    assert new > 0
    recycle_everything(h)
    print('3 parents + per-locus inheritance OK, |pheno - sum(alleles)| = %.0e, sex ratio z=%.2f, parent chi2 %s/%s, mutation %d/%d' % (err, z, sm_, sf_, new, copies))

    # ================================================================ 4. count rule, markers, direct-birth germination, RNG contract
    h = make_hab(L=10, global_index=1)
    h.hab_initialize({'sp2': sp2}, TABLE)
    pairs = h.hab_sexual_pairwise_parents_num()
    out['birth_fraction'] = {}
    for rate in (0.37, 0.81):
        lo = int(pairs * rate); frac = pairs * rate - lo; extra = 0
        for t in range(400):
            np.random.seed(1000 + t)
            r = h.hab_sex_reproduce_mutate_into_offspring_pool(rate, None, None)
            assert r in (lo, lo + 1) and len(h.offspring_pool) == r and all(k.nemo_indi_obj.getHome() == 1 for k in h.offspring_pool)
            extra += r - lo
        assert abs(extra / 400 - frac) < 4 * math.sqrt(frac * (1 - frac) / 400) + 1e-9
        out['birth_fraction'][rate] = [round(frac, 3), extra / 400]
        m0 = sp2.minted_num
        for t in range(100):
            np.random.seed(5000 + t)
            r = h.hab_sex_reproduce_calculation_into_offspring_marker_pool('patchX', rate)
            assert r in (lo, lo + 1) and h.offspring_marker_pool == [('patchX', 'h1', 'sexual')] * r
        assert sp2.minted_num == m0
    np.random.seed(3); h.hab_sex_reproduce_mutate_into_offspring_pool(0.37, None, None); st = np.random.get_state()
    np.random.seed(3); np.random.uniform(0, 1, 1)
    assert np.array_equal(np.random.get_state()[1], st[1]) and np.random.get_state()[2] == st[2], 'np.random use != one uniform per call'
    h.hab_clear_up_offspring_pool()
    for t in range(5):
        np.random.seed(t); h.hab_dead_selection(0.3, 0.5, 'niche_gaussian')
        empty, m0 = len(h.empty_site_pos_ls), sp2.minted_num
        st = np.random.get_state()[1].copy()
        r = h.hab_sexual_reprodece_germinate(1.0, None, None)
        assert r == sp2.minted_num - m0 and r <= empty and np.array_equal(np.random.get_state()[1], st), 'direct birth leaked or used np.random'
        check_invariants(h, FACT, 'reprodece_germinate %d' % t)
    recycle_everything(h)
    print('4 count rule (fractional birth %s), markers mint nothing, direct-birth germination places every newborn' % out['birth_fraction'])

    # ================================================================ 5. two species in one habitat
    h = make_hab(L=10, global_index=2, mean=(0.5, 0.5))
    h.hab_initialize({'sp1': sp1, 'sp2': sp2}, TABLE)
    for k, (i, j) in enumerate([(i, j) for i in range(10) for j in range(10) if (i + j) % 2 == 0]):
        x = h.set['microsite_individuals'][i][j]; h.del_individual(i, j); x.recycle()
        h.add_individual(sp2.make_founder('female' if k % 2 else 'male', 2), i, j)
    by_id = {x.individual_id: x for x in microsite_individuals(h)}
    h.hab_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
    assert {k.species_id for k in h.offspring_pool} == {'sp1', 'sp2'}
    assert all(by_id[k.mother_id].species_id == by_id[k.father_id].species_id == k.species_id and k.indfactory is by_id[k.mother_id].indfactory for k in h.offspring_pool)
    h.hab_local_germinate_from_offspring_and_immigrant_pool(); h.hab_clear_up_offspring_pool()
    check_invariants(h, FACT, 'two species')
    recycle_everything(h)
    print('5 two species: same-species parents, offspring minted by the mother\'s factory, both factories balanced')

    # ================================================================ 6. degenerate parent sets and rates
    h = make_hab(L=4); h.hab_initialize({'sp2': sp2}, TABLE, 'asexual')
    m0 = sp2.minted_num
    assert h.hab_sexual_pairwise_parents_num() == 0 and h.hab_sex_reproduce_mutate_into_offspring_pool(0.99, None, None) == 0 and h.hab_sexual_reprodece_germinate(1.0, None, None) == 0 and sp2.minted_num == m0
    recycle_everything(h)
    h = make_hab(L=3); h.add_individual(sp2.make_founder('female', 3), 0, 0); h.add_individual(sp2.make_founder('male', 3), 2, 2)
    assert h.hab_sex_reproduce_mutate_into_offspring_pool(1.0, None, None) == 1
    recycle_everything(h)
    h = make_hab(L=6); h.hab_initialize({'sp2': sp2}, TABLE); pairs = h.hab_sexual_pairwise_parents_num()
    assert h.hab_sex_reproduce_mutate_into_offspring_pool(0.0, None, None) == 0
    h.hab_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
    first = list(h.offspring_pool)
    h.hab_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)                       # leftovers recycled, not dropped
    assert first and not any(k.is_alive for k in first)
    try: h.hab_sex_reproduce_mutate_into_offspring_pool(1.5, None, None); r15 = 'no error'
    except ValueError: r15 = 'ValueError'
    check_invariants(h, FACT, 'rate 1.5 aftermath')
    out['rate_1.5'] = r15
    finding('info', 'sexual_birth_rate > 1: int(pairs*rate) > pairs -> %s from random.sample (one offspring per pair per call, as in the original); the old offspring_pool is already recycled, invariants hold.' % r15)
    recycle_everything(h)
    for f in FACT: assert f.minted_num == f.recycled_num, f
    print('6 single sex / one pair / rate 0 / leftover offspring recycled / rate 1.5 -> %s' % r15)
    return out


def section_lifecycle():
    out = {'v3': {}}
    for mode in ('individual', 'dormA', 'dormB'):
        for L in (6, 12, 30):
            for seed in (1, 2):
                _, st = v3_run(L, seed, mode)
                out['v3']['%s L=%d seed=%d' % (mode, L, seed)] = st
                if mode != 'individual': assert st['recoveries'] >= 1
    print('1 v3 loops: %d configurations, invariants after every phase, recovery from dormancy in every dormancy run' % len(out['v3']))

    # ================================================================ 2. cross-process reproducibility
    out['reproducibility'] = {}
    for L, seed, mode in ((12, 7, 'individual'), (12, 7, 'dormB')):
        tr = [subprocess.run([sys.executable, os.path.join(HERE, 'test_nemo.py'), '--trace', str(L), str(seed), mode], capture_output=True, text=True, cwd=HERE) for _ in range(2)]
        tr = [json.loads([l for l in p.stdout.splitlines() if l.startswith('TRACE ')][-1][6:]) for p in tr]
        assert tr[0] == tr[1], 'not reproducible: %s L=%d seed=%d' % (mode, L, seed)
        out['reproducibility']['%s L=%d seed=%d' % (mode, L, seed)] = 'identical over %d steps' % len(tr[0])
    print('2 reproducibility in two fresh processes:', out['reproducibility'])


    # ================================================================ 3. selection response
    def evolve(sp, env, method, gens=60, seed=0, offset_at=None):
        seed_all(seed)
        h = make_hab(L=20, mean=(0.2, 0.2))
        h.hab_initialize({'sp1': sp}, [[0.2, 0.2]])
        if tuple(env) != (0.2, 0.2): h.hab_reset_environment_values(['env1', 'env2'], list(env), [0.025, 0.025])
        traj = []
        for g in range(gens):
            if g == offset_at: h.hab_offset_environment_values(['env1', 'env2'], [0.15, 0.15])
            h.hab_dead_selection(0.1, 0.5, method)
            h.hab_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
            h.hab_local_germinate_from_offspring_and_immigrant_pool(); h.hab_clear_up_offspring_pool()
            traj.append(float(np.mean([x.get_indi_phenotype_ls()[0] for x in microsite_individuals(h)])))
        check_invariants(h, sp, 'evolve'); recycle_everything(h)
        return traj


    ss = make_species('ss', 0.2, seed=43, mutation_rate=0.001)
    sel, neu, off = evolve(ss, (0.35, 0.35), 'niche_gaussian', seed=1), evolve(ss, (0.35, 0.35), 'neutral_uniform', seed=1), evolve(ss, (0.2, 0.2), 'niche_gaussian', seed=2, offset_at=10)
    out['selection'] = {'selection_0.35': [round(sel[0], 3), round(sel[-1], 3)], 'neutral': [round(neu[0], 3), round(neu[-1], 3)], 'offset_+0.15_at_gen10': [round(off[9], 3), round(off[-1], 3)]}
    assert sel[-1] - sel[0] > 0.06 and abs(neu[-1] - neu[0]) < 0.06 and off[-1] - off[9] > 0.05, out['selection']
    print('3 selection response:', out['selection'])

    # ================================================================ 4. several habitats, manual dispersal (no patch layer)
    seed_all(71)
    s1, s2 = make_species('sp1', 0.2, seed=71), make_species('sp2', 0.5, seed=72)
    habs = [make_hab(L=8, K=30, global_index=g, mean=m, name='h%d' % g) for g, m in ((10, (0.2, 0.2)), (11, (0.2, 0.2)), (12, (0.5, 0.5)))]
    for h in habs: h.hab_initialize({'sp1': s1, 'sp2': s2}, [[0.2, 0.2], [0.5, 0.5]])
    moved = []
    for t in range(40):
        for h in habs:
            h.hab_dead_selection(0.1, 0.5, 'niche_gaussian'); h.hab_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
            assert all(k.nemo_indi_obj.getHome() == h.global_index for k in h.offspring_pool)
        for i, h in enumerate(habs):
            mig = h.sample_offspring_without_replacement(random.randint(0, 4), ('offspring_pool', 'dormancy_pool'))
            habs[(i + 1) % 3].immigrant_pool += mig
            moved += [(m, m.nemo_indi_obj.getHome()) for m in mig]
        check_invariants(habs, [s1, s2], 'dispersal t=%d' % t)
        for h in habs: h.hab_local_germinate_from_offspring_immigrant_and_dormancy_pool(); h.hab_dormancy_process_from_offspring_pool_and_immigrant_pool()
        if t == 20:
            snap = [(list(h.dormancy_pool), h.indi_num) for h in habs[1:]]
            habs[0].habitat_disturbance_process()
            assert all(h.dormancy_pool == s[0] and h.indi_num == s[1] for h, s in zip(habs[1:], snap)), 'disturbance leaked into another habitat'
        check_invariants(habs, [s1, s2], 't=%d' % t)
    assert all(m.nemo_indi_obj.getHome() == home for m, home in moved if m.is_alive), 'migrant lost its birth home'
    out['multi_habitat'] = {'moved': len(moved), 'habitats_with_two_species': [h.name for h in habs if len({x.species_id for x in microsite_individuals(h)}) > 1]}
    recycle_everything(habs)
    assert s1.minted_num == s1.recycled_num and s2.minted_num == s2.recycled_num
    print('4 three habitats, two species:', out['multi_habitat'])

    # ================================================================ 5. fuzz: random operation orders
    seed_all(111)
    f1, f2 = make_species('sp1', 0.2, seed=111), make_species('sp2', 0.3, seed=112)
    F = {'sp1': f1, 'sp2': f2}; TAB = [[0.2, 0.2], [0.3, 0.3]]
    n_ops = 0
    for seed in range(30):
        rng = random.Random(seed); seed_all(seed)
        L, K = rng.choice((1, 2, 4, 7, 10)), rng.choice((0, 1, 10, 40))
        h, other = make_hab(L=L, K=K, global_index=seed, mean=rng.choice(((0.2, 0.2), (0.3, 0.3)))), make_hab(L=3, K=5, global_index=100 + seed)
        h.hab_initialize(F, TAB)

        def kill_one():
            if h.occupied_site_pos_ls:
                i, j = rng.choice(h.occupied_site_pos_ls); x = h.set['microsite_individuals'][i][j]; h.del_individual(i, j); x.recycle()
        OPS = [lambda: h.hab_dead_selection(rng.choice((0.0, 0.1, 0.5, 1.0)), rng.choice((0.2, 0.5, 2.0)), rng.choice(('niche_gaussian', 'neutral_uniform'))),
               lambda: h.hab_sex_reproduce_mutate_into_offspring_pool(rng.choice((0.0, 0.3, 0.77, 1.0)), None, None),
               lambda: h.hab_sex_reproduce_calculation_into_offspring_marker_pool('p', rng.random()),
               lambda: h.hab_sexual_reprodece_germinate(rng.choice((0.5, 1.0)), None, None),
               h.hab_local_germinate_from_offspring_and_immigrant_pool, h.hab_local_germinate_from_offspring_and_dormancy_pool,
               h.hab_local_germinate_from_offspring_immigrant_and_dormancy_pool, h.hab_dormancy_process_from_offspring_pool_to_dormancy_pool,
               h.hab_dormancy_process_from_offspring_pool_and_immigrant_pool, h.habitat_disturbance_process, h.hab_clear_up_offspring_pool,
               h.hab_clear_up_immigrant_pool, h.hab_clear_up_dormancy_pool, h.hab_clear_up_microsite_individuals, h.hab_clear_up_offspring_marker_pool,
               lambda: h.immigrant_pool.extend(rng.choice((f1, f2)).make_founder(rng.choice(('male', 'female')), rng.randint(0, 9)) for _ in range(rng.randint(0, 6))),
               lambda: other.immigrant_pool.extend(h.sample_offspring_without_replacement(rng.randint(0, 5), rng.choice((('offspring_pool',), ('offspring_pool', 'dormancy_pool'))))),
               kill_one, lambda: h.hab_initialize(F, TAB, rng.choice(('sexual', 'asexual')))]
        for step in range(150):
            k = rng.randrange(len(OPS)); OPS[k]()
            check_invariants([h, other], [f1, f2], 'fuzz seed=%d step=%d op=%d' % (seed, step, k))
            assert len(h.dormancy_pool) <= K
            n_ops += 1
        recycle_everything([h, other])
        assert f1.minted_num == f1.recycled_num and f2.minted_num == f2.recycled_num
    out['fuzz_operations'] = n_ops
    print('5 fuzz: %d random operations, invariants after each' % n_ops)
    return out


def section_memory():
    out = {}
    seed_all(51)
    sp = make_species('sp1', 0.2, seed=51)
    ref = make_hab(L=5, global_index=9)                       # 25 individuals that stay untouched for the whole test
    ref.hab_initialize({'sp1': sp}, [[0.2, 0.2]])
    state = lambda x: (x.gender, x.get_indi_phenotype_ls(), [x.genotype_table(a) for a in range(2)], x.nemo_indi_obj.getHome(), x.mother_id, x.father_id)
    snap = {x.individual_id: state(x) for x in microsite_individuals(ref)}

    churn = make_hab(L=20, K=60, global_index=4)
    churn.hab_initialize({'sp1': sp}, [[0.2, 0.2]])
    max_id = max(x.individual_id for x in held_individuals(churn) + held_individuals(ref))
    dead = []
    for t in range(150):
        churn.hab_dead_selection(0.3, 0.4, 'niche_gaussian')
        churn.hab_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
        born = list(churn.offspring_pool)
        assert all(k.individual_id > max_id for k in born), 'ID reused at t=%d' % t
        max_id = max([max_id] + [k.individual_id for k in born])
        for _ in range(random.randint(0, 10)): churn.immigrant_pool.append(sp.make_founder('female', 7))
        churn.hab_local_germinate_from_offspring_immigrant_and_dormancy_pool()
        (churn.hab_dormancy_process_from_offspring_pool_to_dormancy_pool if t % 2 else churn.hab_dormancy_process_from_offspring_pool_and_immigrant_pool)()
        if t % 37 == 36: churn.habitat_disturbance_process()
        held = held_individuals(churn) + held_individuals(ref)
        assert len({int(x.nemo_indi_obj.this) for x in held}) == len(held), 'two held individuals share one C++ Individual at t=%d' % t
        assert len({x.individual_id for x in held}) == len(held), 'duplicate live ID at t=%d' % t
        check_invariants([churn, ref], sp, 't=%d' % t)
        dead += [k for k in born if not k.is_alive][:3]
    assert all(state(x) == snap[x.individual_id] for x in microsite_individuals(ref)), 'a living individual changed while memory was recycled / reused'
    for x in dead[:30]:
        assert x.nemo_indi_obj is None
        for fn in (lambda: x.gender, lambda: x.get_indi_phenotype_ls(), lambda: x.individual_id, lambda: x.recycle()):
            try: fn(); raise AssertionError('recycled individual still usable')
            except RuntimeError: pass
    out['churn'] = {'births': sp.minted_num, 'recycles': sp.recycled_num}
    recycle_everything([churn, ref])
    assert sp.minted_num == sp.recycled_num
    sp.purge_recycling_pool()
    print('reference individuals unchanged after %(births)d births / %(recycles)d recycles; pointers + IDs unique; recycled individuals inert' % out['churn'])

    seed_all(52)
    sp2 = make_species('sp2', 0.2, seed=52)
    h = make_hab(L=30, K=200, global_index=1)
    h.hab_initialize({'sp1': sp2}, [[0.2, 0.2]])
    rss = {}
    for t in range(400):
        h.hab_dead_selection(0.2, 0.5, 'niche_gaussian')
        h.hab_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
        h.hab_local_germinate_from_offspring_immigrant_and_dormancy_pool()
        h.hab_dormancy_process_from_offspring_pool_to_dormancy_pool()
        if t in (50, 200, 399): rss[t] = round(rss_mb(), 1)
    check_invariants(h, sp2, 'long run')
    out['long_run'] = {'births': sp2.minted_num, 'rss_mb': rss}
    assert rss[399] - rss[50] < 40, 'RSS keeps growing in a steady-state run'
    recycle_everything(h)
    assert sp2.minted_num == sp2.recycled_num
    print('long run L=30 x 400 steps, %d births: RSS %s MB' % (sp2.minted_num, rss))
    return out


def section_edges_known():
    out = {}
    seed_all(61)
    sp = make_species('sp1', 0.2, seed=61)
    F = {'sp1': sp}

    # ================================================================ 1. shapes + empty habitat: every NemoBackend function, invariants after each
    def every_function(h, tag):
        calls = [lambda: h.hab_dead_selection(0.1, 0.5, 'niche_gaussian'), lambda: h.hab_sex_reproduce_mutate_into_offspring_pool(1.0, None, None),
                 h.hab_local_germinate_from_offspring_and_immigrant_pool, h.hab_local_germinate_from_offspring_and_dormancy_pool,
                 h.hab_local_germinate_from_offspring_immigrant_and_dormancy_pool, h.hab_dormancy_process_from_offspring_pool_to_dormancy_pool,
                 h.hab_dormancy_process_from_offspring_pool_and_immigrant_pool, lambda: h.hab_sex_reproduce_calculation_into_offspring_marker_pool('p', 0.7),
                 lambda: h.hab_sexual_reprodece_germinate(1.0, None, None), h.habitat_disturbance_process, h.hab_clear_up_microsite_individuals,
                 h.hab_clear_up_offspring_pool, h.hab_clear_up_immigrant_pool, h.hab_clear_up_dormancy_pool, h.hab_clear_up_offspring_marker_pool,
                 h.hab_clear_up_immigrant_marker_pool]
        res = []
        for c in calls:
            res.append(c()); check_invariants(h, sp, tag)
        return res


    for L, W in ((1, 1), (1, 9), (9, 1), (3, 7), (7, 3)):
        np.random.seed(L * 10 + W)
        h = make_hab(L=L, W=W, K=5, global_index=1)
        assert len(h.set['microsite_individuals']) == L and all(len(r) == W for r in h.set['microsite_individuals']) and h.set['env1'].shape == (L, W)
        assert h.empty_site_pos_ls == [(i, j) for i in range(L) for j in range(W)]
        for rep in range(3):
            h.hab_initialize(F, [[0.2, 0.2]]); assert h.indi_num == L * W
            every_function(h, '%dx%d' % (L, W))
        recycle_everything(h)
    res = every_function(make_hab(L=4, K=5), 'empty')
    assert all(v in (0, (0, 0, 0, 0)) for v in res), res
    print('1 shapes 1x1 / 1x9 / 9x1 / 3x7 / 7x3 and an empty habitat: every function, invariants hold, empty -> all zero')

    # ================================================================ 2. constructor + inputs
    h = habitat('h', 0, (0, 0), 2, ['env1', 'env2'], [0.2, 0.2], [0.025, 0.025], 4, 4)
    assert h.global_index == 0 and h.dormancy_pool_max_size == 0
    h = habitat('h', 0, (0, 0), 2, ['env1', 'env2'], [0.2, 0.2], [0.025, 0.025], 4, 5, 25)        # 10th positional = dormancy_pool_max_size (as patch.add_habitat passes it)
    assert h.dormancy_pool_max_size == 25 and h.global_index == 0 and h.size == 20
    gi_res = {}
    for gi in (0, 65535, 65536, -1):
        h = make_hab(L=2, global_index=gi)
        try: h.hab_initialize(F, [[0.2, 0.2]]); gi_res[gi] = 'ok, getHome()=%d' % microsite_individuals(h)[0].nemo_indi_obj.getHome()
        except Exception as e: gi_res[gi] = type(e).__name__
        check_invariants(h, sp, 'gi'); recycle_everything(h)
    assert gi_res[0] == 'ok, getHome()=0' and gi_res[65535] == 'ok, getHome()=65535'
    finding('info', 'global_index is not validated by habitat.__init__; out of the unsigned-short range Nemo refuses it at the first make_founder: %s (nothing minted).' % {k: v for k, v in gi_res.items() if k in (65536, -1)})
    errs = {}
    for label, hab, args in (('mean not in species table', make_hab(L=2), (F, [[0.5, 0.5]])), ('missing indfactory key', make_hab(L=2), ({'spX': sp}, [[0.2, 0.2]])),
                             ('tuple mean_env_ls', habitat('h', 0, (0, 0), 2, ['env1', 'env2'], (0.2, 0.2), [0.025, 0.025], 2, 2), (F, [[0.2, 0.2]]))):
        try: hab.hab_initialize(*args); errs[label] = 'no error'
        except Exception as e: errs[label] = type(e).__name__
    assert sp.minted_num == sp.recycled_num
    finding('info', 'hab_initialize finds the species by exact equality species_2_phenotype_ls.index(mean_env_ls) and names it "sp<k+1>": %s; no founder is minted on these errors. hab_reset/offset_environment_values change mean_env_ls in place, so hab_initialize after an environment change also raises (as in the original).' % errs)
    h = make_hab(L=3)
    try: h.hab_dead_selection(0.1, 0.5, 'gaussian'); r = 'no error'
    except Exception as e: r = type(e).__name__
    finding('info', 'unknown dead-selection method string ("gaussian") -> %s from survival_rate (no validation; unchanged from the original).' % r)
    print('2 constructor defaults / positional K, global_index bounds %s, bad species inputs %s' % (gi_res, errs))

    # ================================================================ 3. legacy helpers with Nemo individuals
    h = make_hab(L=6, K=5)
    h.hab_initialize(F, [[0.2, 0.2]])
    assert '<individual species_id=' in str(h)
    assert h.get_microsite_env_val_ls(2, 3) == [h.set['env1'][2][3], h.set['env2'][2][3]]
    for d, w, ph, env in ((0.1, 0.5, [0.2, 0.3], [0.25, 0.1]), (0.5, 0.2, [1.0, 1.0], [0.0, 0.0])):
        assert abs(h.survival_rate(d, ph, env, w) - (1 - d) * math.exp(-((ph[0] - env[0]) / w) ** 2) * math.exp(-((ph[1] - env[1]) / w) ** 2)) < 1e-15
        assert h.survival_rate(d, ph, env, w, 'neutral_uniform') == 1 - d
    np.random.seed(0)
    h.hab_reset_environment_values(['env2'], [0.6], [0.0]); assert np.allclose(h.set['env2'], 0.6) and h.mean_env_ls == [0.2, 0.6]
    h.hab_offset_environment_values(['env1', 'env2'], [0.1, -0.1], [-0.025, 0.0]); assert np.allclose(h.set['env1'], 0.3) and np.allclose(h.set['env2'], 0.5)
    check_invariants(h, sp, 'env change')
    for k in range(13): x = h.set['microsite_individuals'][k // 6][k % 6]; h.del_individual(k // 6, k % 6); x.recycle()
    empty = list(h.empty_site_pos_ls)
    pe = h.get_hab_pairwise_empty_site_pos_ls(); flat = [p for pr in pe for p in pr]
    assert len(pe) == 6 and len(set(flat)) == 12 and set(flat) <= set(empty) and h.empty_site_pos_ls == empty
    cat = copy.deepcopy(h.species_category)
    po = h.get_hab_pairwise_occupied_site_pos_ls(); flat = [p for pr in po for p in pr]
    assert all(h.set['microsite_individuals'][f[0]][f[1]].gender == 'female' and h.set['microsite_individuals'][m[0]][m[1]].gender == 'male' for f, m in po)
    assert len(set(flat)) == len(flat) and h.species_category == cat
    random.seed(1); h.hab_sexual_pairwise_parents_num()
    assert all(sorted(h.species_category[s][g]) == sorted(cat[s][g]) for s in cat for g in cat[s])
    finding('info', 'hab_sexual_pairwise_parents_ls() / _num() shuffle the position lists in self.species_category in place, so even a count consumes Python RNG (same as the original; bookkeeping stays consistent).')
    check_invariants(h, sp, 'helpers'); recycle_everything(h)
    print('3 __str__, env helpers, survival_rate formula, pairwise empty / occupied helpers OK with Nemo individuals')

    # ================================================================ 4. asex / mix (not migrated): one probe per function, each in its own process
    probes = []
    for i in range(len(ASEX_MIX)):
        p = subprocess.run([sys.executable, os.path.join(HERE, 'test_nemo.py'), '--probe', str(i)], capture_output=True, text=True, cwd=HERE)
        line = [l for l in p.stdout.splitlines() if l.startswith('PROBE ')]
        probes.append(json.loads(line[-1][6:]) if line else [ASEX_MIX[i][0], 'process died rc=%d' % p.returncode])
    out['asex_mix_probes'] = probes
    finding('info', 'asex / mix functions are not migrated: ' + '; '.join('%s -> %s' % (n, r) for n, r in probes))
    print('4 asex / mix probes:', ['%s: %s' % (n, r.split(':')[0]) for n, r in probes])

    # ================================================================ 5. caller-contract violations (what the state looks like; these errors stop a run)
    def state(h):
        try: check_invariants(h, sp, 's'); return 'invariants hold'
        except (AssertionError, RuntimeError) as e: return 'invariants broken'


    h = make_hab(L=3); x = sp.make_founder('female', 3); x.recycle(); h.offspring_pool = [x]
    try: h.hab_local_germinate_from_offspring_and_immigrant_pool(); r = 'no error'
    except RuntimeError: r = 'RuntimeError'
    finding('info', 'germinating an already-recycled individual -> %s, %s: add_individual writes grid / empty / occupied / indi_num before it reads gender.' % (r, state(h)))
    h = make_hab(L=2); a, b = fake_indi('a'), fake_indi('b'); h.add_individual(a, 0, 0); h.add_individual(b, 0, 0)
    assert h.indi_num == 1 and b.is_alive
    finding('info', 'add_individual on an occupied microsite only prints and returns; the individual stays with the caller (hab_initialize resets the microsites first, so it never hits this).')
    x = fake_indi('dup'); h = make_hab(L=2, K=10); h.offspring_pool = [x]; h.dormancy_pool = [x]
    h.hab_dormancy_process_from_offspring_pool_to_dormancy_pool()
    finding('info', '_recycle_dropped detects a two-owner individual only when it is dropped; if kept it may end up twice in dormancy_pool (copies=%d here). Callers must never create two owners.' % sum(y is x for y in h.dormancy_pool))
    finding('info', 'clear-up functions empty a pool before its recycle loop: if a recycle raises (two owners), the rest of that pool is neither recycled nor held. The error stops the run anyway.')
    finding('info', 'sample_offspring_without_replacement() hands ownership to the caller (patch dispersal): the returned individuals must go into a pool or be recycled -- to verify when patch.py is migrated.')
    finding('info', 'habitat_disturbance_process() returns recycled individuals + cleared markers; the original returned None and kept the marker pools.')
    assert sp.minted_num - sp.recycled_num == 0
    print('5 caller-contract violations recorded')
    return out


if __name__ == '__main__':
    results, failed = {}, []
    for name in ('reproduction', 'lifecycle', 'memory', 'edges_known'):
        print('==== section %s' % name)
        try:
            results[name] = globals()['section_' + name]()
            print('---- section %s PASSED' % name)
        except Exception:
            failed.append(name); results[name] = {'error': traceback.format_exc()}
            print('---- section %s FAILED\n%s' % (name, results[name]['error']))
    results['failed_sections'] = failed
    save('nemo', results)
    if failed: sys.exit('FAILED sections: %s' % failed)
    print('PASSED')

