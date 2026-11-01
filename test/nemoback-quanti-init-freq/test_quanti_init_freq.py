# -*- coding: utf-8 -*-
''' test_quanti_init_freq.py -- indfactory(metaibm_nemo_params={'quanti_init_freq': ...}): founder quanti alleles redrawn in Python (NemoBackend)

    Run (pynemo with TTrait.set_allele_bit -> conda env metaibm-nemo-test3):
        PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test3/bin/python test_quanti_init_freq.py [group ...]

    Contract under test (the user follows it): quanti_init_freq is one of 0.8 | [0.8] | [[0.8]] | [p_1, ..., p_L]
        (L = all quanti loci, trait by trait; p = frequency of allele 0 = +a in 'diallelic', the ancestral allele in 'diallelic_asymmetrical').
    Every group runs in its own child process (Nemo's RAND, Individual::currentID and prototypes are process-global).
    A child prints one '@@ {json}' line per check: status PASS | FAIL | NOTE (NOTE = behaviour documented, not judged).
    Nothing is written outside a temporary directory; metaibm/ is imported read-only (no __pycache__).
'''
import sys
sys.dont_write_bytecode = True                  # never drop __pycache__ into the user's metaibm/
import os, json, time, shutil, hashlib, tempfile, subprocess, random
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
GROUPS = ('forms', 'mapping', 'stats', 'rng', 'activation', 'lifecycle', 'genetics', 'recycle', 'multi_species',
          'nemo_lce', 'metaibm_lce', 'regression', 'user_notes', 'init_value_note', 'perf',
          'sab_symmetric', 'sab_asymmetrical', 'sab_continuous', 'sab_ntrl')
UNSET = object()

##****************************** child side: helpers ******************************
def _import_stack(root):
    sys.path.insert(0, str(root))
    import numpy as np
    from pynemo import _nemo, nemo
    from pynemo.params import SimParameters
    from pynemo.sim import NemoSim
    from metaibm.indfactory import indfactory
    _nemo.cvar.SILENT_RUN = True
    return np, nemo, SimParameters, NemoSim, indfactory

class Ctx:
    ''' one child process: imports, a temp Nemo cache, check() reporting '''
    def __init__(self, group, root):
        self.group = group
        self.np, self.nemo, self.SimParameters, self.NemoSim, self.indfactory = _import_stack(root)
        self.cache = tempfile.mkdtemp(prefix='qif_%s_' % group)
        self.n = 0

    def check(self, cid, ok, detail='', note=False):
        status = 'NOTE' if note else ('PASS' if ok else 'FAIL')
        print('@@ ' + json.dumps({'group': self.group, 'id': cid, 'status': status, 'detail': detail}), flush=True)
        return ok

    def data(self, key, value):
        print('@@DATA ' + json.dumps({'key': key, 'value': value}), flush=True)

    def params(self, model='diallelic', loci=(3, 2), allele_value='{{0.1}}', ve=None, mut=0.0, pleio=None, extra=(), lce=None):
        ''' one species' SimParameters; loci = loci per trait (pleio: (L,)*T); lce = quanti_init_freq string -> Nemo's own LCE quanti_init '''
        self.n += 1
        p = self.SimParameters()
        cap = 1 if lce is None else lce[1]
        rows = [('replicates', 1), ('generations', 1), ('patch_number', 1), ('patch_capacity', cap), ('aging', 1 if lce is None else 2),
                ('run_mode', 'silent'), ('root_dir', self.cache), ('logfile', self.cache + '/nemo.log'), ('filename', 'sp%d' % self.n),
                ('quanti_traits', len(loci)), ('quanti_loci', loci[0] if pleio else sum(loci)), ('quanti_allele_model', model),
                ('quanti_mutation_rate', mut), ('quanti_recombination_rate', 0.5)]
        if len(loci) > 1 and not pleio: rows.append(('quanti_loci_per_trait', '{{%s}}' % ', '.join(str(x) for x in loci)))
        if pleio: rows.append(('quanti_pleiotropy', pleio))
        if model.startswith('diallelic'):
            if allele_value is not None: rows.append(('quanti_allele_value', allele_value))
        else: rows.append(('quanti_mutation_variance', 0.001))
        if ve is not None: rows.append(('quanti_environmental_variance', ve))
        if lce is not None: rows += [('quanti_init', 1), ('quanti_init_freq', lce[0])]
        for k, v in rows + list(extra): p.add(k, v)
        return p

    def factory(self, freq=UNSET, seed=1, mnp=UNSET, loci=(3, 2), **kw):
        ''' one indfactory; freq -> metaibm_nemo_params={'quanti_init_freq': freq}; mnp passes metaibm_nemo_params verbatim (None/{} allowed) '''
        if mnp is UNSET: mnp = UNSET if freq is UNSET else {'quanti_init_freq': freq}
        p = self.params(loci=loci, **kw)
        extra = {} if mnp is UNSET else {'metaibm_nemo_params': mnp}
        return self.indfactory('sp%d' % self.n, pheno_name_ls=tuple('t%d' % i for i in range(len(loci))), nemo_params=p, trait_keys=('quant',), seed=seed, **extra)

    @staticmethod
    def positions(trait, traits_num):
        ''' sequence positions trait by trait, read off Nemo's tables (same order as quanti_init_freq) '''
        return [trait.get_locus_seq_pos(l, d) for d in range(traits_num) for l in range(trait.get_num_loci(d))]

    def bits(self, ind):
        t = ind.get_trait('quant')
        P = self.positions(t, ind.traits_num)
        return self.np.array([[t.get_allele(pos, s) for pos in P] for s in (0, 1)], dtype=int)

    def values(self, ind):
        t = ind.get_trait('quant')
        P = self.positions(t, ind.traits_num)
        return self.np.array([[t.get_allele_value(pos, s) for pos in P] for s in (0, 1)], dtype=float)

    def genotype(self, ind):
        ''' additive genotypic value per trait = sum of both strands' allele values over that trait's loci '''
        t = ind.get_trait('quant')
        return [sum(t.get_allele_value(t.get_locus_seq_pos(l, d), s) for l in range(t.get_num_loci(d)) for s in (0, 1)) for d in range(ind.traits_num)]

    def close(self):
        shutil.rmtree(self.cache, ignore_errors=True)

def _expected_values(np, model, a, bits):
    ''' allele values from allele identities: diallelic 0 -> +a, 1 -> -a; diallelic_asymmetrical 0 -> 0, 1 -> a (a per locus) '''
    a = np.broadcast_to(np.asarray(a, dtype=float), bits.shape[1:])
    if model == 'diallelic_asymmetrical': return np.where(bits == 0, 0.0, a)
    return np.where(bits == 0, a, -a)

FORMS = lambda L: (('0.8', 0.8), ('[0.8]', [0.8]), ('[[0.8]]', [[0.8]]), ('[0.8]*L', [0.8] * L))

##****************************** groups ******************************
def g_forms(c):
    ''' the four forms draw the same founders, and every draw equals the reference U >= p computed from np.random directly '''
    np = c.np
    for model in ('diallelic', 'diallelic_HC', 'diallelic_asymmetrical'):
        for loci in ((5,), (3, 2), (20, 20), (100, 100)):
            L = sum(loci)
            facs = {name: c.factory(freq=f, model=model, loci=loci) for name, f in FORMS(L)}
            ok, bad = True, []
            for seed in (1, 7, 2026, 99991):
                ref = None
                for name, fac in facs.items():
                    np.random.seed(seed)
                    got = np.array([c.bits(fac.make_founder(g)) for g in ('female', 'male', 'female')])
                    np.random.seed(seed)
                    want = np.array([(np.random.random((2, L)) >= 0.8).astype(int) for _ in range(3)])
                    if not (got == want).all(): ok = False; bad.append('%s seed %d' % (name, seed))
                    if ref is None: ref = got
                    elif not (got == ref).all(): ok = False; bad.append('%s != 0.8 seed %d' % (name, seed))
            c.check('forms/%s/loci%s' % (model, loci), ok, 'four forms x 4 seeds x 3 founders == reference (U >= p, strand-major)' if ok else '; '.join(bad))
    # [p_1..p_L] with distinct p: 200 consecutive founders == consecutive reference draws
    p = [0.05, 0.3, 0.5, 0.7, 0.95]
    for model in ('diallelic', 'diallelic_asymmetrical'):
        fac = c.factory(freq=p, model=model, loci=(3, 2))
        np.random.seed(5)
        got = np.array([c.bits(fac.make_founder('female' if i % 2 else 'male')) for i in range(200)])
        np.random.seed(5)
        want = np.array([(np.random.random((2, 5)) >= np.array(p)).astype(int) for _ in range(200)])
        c.check('forms/per-locus/%s' % model, (got == want).all(), '200 consecutive founders, p=%s: bit-identical to the reference stream' % p)

def g_mapping(c):
    ''' deterministic p in {0, 1}: column k -> trait-by-trait locus k -> its sequence position; phenotype == analytic genotype (Ve unset = 0) '''
    np = c.np
    for model in ('diallelic', 'diallelic_HC', 'diallelic_asymmetrical'):
        for loci in ((4,), (3, 2), (2, 3), (2, 2, 2), (1, 4)):
            L = sum(loci)
            pats = {'all0(p=1)': [1.0] * L, 'all1(p=0)': [0.0] * L, 'alternating': [float(k % 2 == 0) for k in range(L)],
                    'by-trait': [float(d % 2 == 0) for d in range(len(loci)) for _ in range(loci[d])]}
            avals = [0.1 * (k + 1) for k in range(L)]                                       # locus-specific allele values: a wrong column->locus map shows
            ok, msgs = True, []
            for av_name, av in (('scalar-matrix', '{{0.1}}'), ('locus-specific', '{{%s}}' % ', '.join('%g' % x for x in avals))):
                a = 0.1 if av_name == 'scalar-matrix' else avals
                for pname, p in pats.items():
                    fac = c.factory(freq=p, model=model, loci=loci, allele_value=av)
                    for g in ('female', 'male'):
                        ind = fac.make_founder(g)
                        b = c.bits(ind)
                        want_b = np.tile((1 - np.array(p)).astype(int), (2, 1))
                        vals = _expected_values(np, model, a, want_b)
                        want_ph = [vals[:, sum(loci[:d]):sum(loci[:d + 1])].sum() for d in range(len(loci))]
                        ph = ind.get_indi_phenotype_ls()
                        if not (b == want_b).all(): ok = False; msgs.append('%s %s %s bits %s' % (av_name, pname, g, b.tolist()))
                        if not np.allclose(ph, want_ph, atol=1e-12): ok = False; msgs.append('%s %s %s pheno %s != %s' % (av_name, pname, g, ph, want_ph))
                        if not np.allclose(c.values(ind), vals, atol=1e-12): ok = False; msgs.append('%s %s values' % (av_name, pname))
            c.check('mapping/%s/loci%s' % (model, loci), ok, '4 patterns x 2 allele-value layouts x 2 genders: bits, allele values and phenotype exact' if ok else '; '.join(msgs[:4]))
    # the cache equals Nemo's no-pleiotropy layout: position k == k
    fac = c.factory(freq=0.5, loci=(3, 2)); fac.make_founder('female')
    c.check('mapping/seq_pos-identity', fac._quanti_seq_pos == list(range(5)), '_quanti_seq_pos=%s' % fac._quanti_seq_pos)

def g_stats(c):
    ''' many founders: allele frequencies, Hardy-Weinberg, independence, genotypic mean/variance, Ve '''
    np = c.np
    N = 6000
    p = np.array([0.05, 0.3, 0.5, 0.7, 0.95])
    for model in ('diallelic', 'diallelic_asymmetrical'):
        fac = c.factory(freq=list(p), model=model, loci=(3, 2), allele_value='{{0.1}}')
        np.random.seed(17)
        inds = [fac.make_founder('female' if i % 2 else 'male') for i in range(N)]
        B = np.array([c.bits(i) for i in inds])                                             # (N, strand, locus)
        a0 = 1 - B
        f = a0.mean(axis=(0, 1)); z = (f - p) / np.sqrt(p * (1 - p) / (2 * N))
        c.check('stats/%s/allele0-freq' % model, (abs(z) < 4).all(), 'p=%s observed=%s max|z|=%.2f' % (p.tolist(), f.round(4).tolist(), abs(z).max()))
        fs = a0.mean(axis=0); zs = (fs - p) / np.sqrt(p * (1 - p) / N)
        c.check('stats/%s/per-strand' % model, (abs(zs) < 4).all(), 'strand0=%s strand1=%s max|z|=%.2f' % (fs[0].round(4).tolist(), fs[1].round(4).tolist(), abs(zs).max()))
        n00 = (a0.sum(axis=1) == 2).mean(axis=0); n01 = (a0.sum(axis=1) == 1).mean(axis=0); n11 = (a0.sum(axis=1) == 0).mean(axis=0)
        exp = np.array([p ** 2, 2 * p * (1 - p), (1 - p) ** 2])
        chi = (N * (np.array([n00, n01, n11]) - exp) ** 2 / exp).sum(axis=0)
        c.check('stats/%s/hardy-weinberg' % model, (chi < 18.4).all(), 'chi2 (df 2, known p) per locus=%s (limit 18.4 ~ P 1e-4)' % chi.round(2).tolist())
        X = a0.reshape(N, 10).astype(float)                                                 # 10 allele columns
        R = np.corrcoef(X, rowvar=False); off = R[~np.eye(10, dtype=bool)]
        c.check('stats/%s/independence' % model, abs(off).max() < 4.5 / np.sqrt(N), 'max |r| between allele columns (strands x loci)=%.4f (limit %.4f)' % (abs(off).max(), 4.5 / np.sqrt(N)))
        G = np.array([ind.get_indi_phenotype_ls() for ind in inds])
        a = 0.1
        if model == 'diallelic': mu1, v1 = a * (2 * p - 1), 4 * a * a * p * (1 - p)
        else: mu1, v1 = a * (1 - p), a * a * p * (1 - p)
        mu = [2 * mu1[:3].sum(), 2 * mu1[3:].sum()]; var = [2 * v1[:3].sum(), 2 * v1[3:].sum()]
        zm = (G.mean(axis=0) - mu) / np.sqrt(np.array(var) / N); rv = G.var(axis=0) / var
        c.check('stats/%s/genotypic-mean-var' % model, (abs(zm) < 4).all() and (abs(rv - 1) < 0.1).all(), 'mean=%s (expect %s) var/expected=%s' % (G.mean(axis=0).round(4).tolist(), np.round(mu, 4).tolist(), rv.round(3).tolist()))
    # Ve: phenotype = new genotype + a fresh N(0, Ve) draw
    fac = c.factory(freq=[0.5] * 5, loci=(3, 2), ve=0.01)
    np.random.seed(3)
    inds = [fac.make_founder('female') for _ in range(N)]
    P = np.array([i.get_indi_phenotype_ls() for i in inds]); G = np.array([c.genotype(i) for i in inds]); E = P - G
    c.check('stats/Ve-residual', (abs(E.mean(axis=0)) < 4 * 0.1 / np.sqrt(N)).all() and (abs(E.std(axis=0) / 0.1 - 1) < 0.05).all() and abs(np.corrcoef(E[:, 0], G[:, 0])[0, 1]) < 4 / np.sqrt(N),
            'phenotype - genotype: mean=%s sd=%s (expect 0, 0.1); corr(E, G)=%.4f' % (E.mean(axis=0).round(4).tolist(), E.std(axis=0).round(4).tolist(), np.corrcoef(E[:, 0], G[:, 0])[0, 1]))

def g_rng(c):
    ''' np.random: exactly 2L uniforms per redrawn founder, none otherwise; Python random untouched; seeds reproduce '''
    np = c.np
    for loci in ((3, 2), (20, 20)):
        L = sum(loci)
        fac = c.factory(freq=0.4, loci=loci)
        np.random.seed(11); random.seed(11); r0 = random.getstate()
        for _ in range(7): fac.make_founder('male')
        s1 = np.random.get_state()
        np.random.seed(11); np.random.random(2 * L * 7); s2 = np.random.get_state()
        c.check('rng/np-consumption/L%d' % L, (s1[1] == s2[1]).all() and s1[2] == s2[2], '7 founders consume exactly 7 x 2L = %d uniforms' % (14 * L))
        c.check('rng/python-random-untouched/L%d' % L, random.getstate() == r0, 'random.getstate() unchanged by make_founder')
    for name, kw in (('no-argument', {}), ('None', {'mnp': None}), ('empty-dict', {'mnp': {}}), ('freq-None', {'mnp': {'quanti_init_freq': None}})):
        fac = c.factory(**kw)
        np.random.seed(4); before = np.random.get_state()
        for _ in range(5): fac.make_founder('female')
        after = np.random.get_state()
        c.check('rng/no-freq/%s' % name, (before[1] == after[1]).all() and before[2] == after[2] and fac._quanti_seq_pos is None and not fac._use_quanti_init_freq(),
                'np.random untouched, _quanti_seq_pos stays None, _use_quanti_init_freq() False')
    # same np seed -> same founders (Ve > 0 too: Nemo RAND follows the indfactory seed)
    digs = []
    for _ in range(2):
        fac = c.factory(freq=[0.2, 0.4, 0.6, 0.8, 0.5], loci=(3, 2), ve=0.01, seed=42)
        np.random.seed(8)
        inds = [fac.make_founder('female' if i % 3 else 'male') for i in range(50)]
        digs.append(hashlib.md5(json.dumps([[c.bits(i).tolist(), [repr(x) for x in i.get_indi_phenotype_ls()]] for i in inds]).encode()).hexdigest())
    c.data('rng/digest', digs[0])
    c.check('rng/in-process-repro', digs[0] == digs[1], 'two factories, same seeds: bits+phenotypes md5 %s' % digs[0][:12])
    # the redraw reads np.random only: indfactory seed (Nemo RAND), gender and home_habitat do not change the alleles
    got = []
    for seed, g, home in ((1, 'female', 0), (999, 'male', 7), (123456, 'female', 65535)):
        fac = c.factory(freq=[0.3, 0.5, 0.7, 0.2, 0.9], loci=(3, 2), seed=seed)
        np.random.seed(55)
        got.append([c.bits(fac.make_founder(g, home_habitat=home)).tolist() for _ in range(10)])
    c.check('rng/independent-of-nemo-seed-gender-home', got[0] == got[1] == got[2], 'indfactory seed 1/999/123456 x gender x home_habitat 0/7/65535: identical alleles for one np seed')

def g_activation(c):
    ''' which allele models / pleiotropy settings redraw vs raise ValueError '''
    np = c.np
    for model in ('diallelic', 'diallelic_HC', 'diallelic_asymmetrical'):
        fac = c.factory(freq=0.0, model=model)
        b = c.bits(fac.make_founder('female'))
        c.check('activation/%s' % model, fac._use_quanti_init_freq() and (b == 1).all(), 'redrawn: p=0 -> every allele is allele 1')
    for model in ('continuous', 'continuous_in_place', 'continuous_HC'):
        fa = c.factory(freq=0.0, model=model, seed=5, extra=[('quanti_init_model', 2), ('quanti_init_variance', 0.01)])
        np.random.seed(2); s0 = np.random.get_state()
        err = None
        try: fa.make_founder('female')
        except ValueError as e: err = str(e)
        s1 = np.random.get_state()
        c.check('activation/%s' % model, err is not None and 'quanti_init_freq' in err and (s0[1] == s1[1]).all(), 'raises ValueError (%s), np.random untouched' % err)
    pleio_errs = {}
    for pl in ('full', 'f', 'var', 'variable', 'v'):
        ex = [('quanti_pleio_matrix', '{{1, 1}{1, 0}{0, 1}}')] if pl.startswith('v') else []
        fa = c.factory(freq=0.5, loci=(3, 3), pleio=pl, extra=ex)
        np.random.seed(2); s0 = np.random.get_state()
        try:
            fa.make_founder('female')
            pleio_errs[pl] = None
        except ValueError as e:
            pleio_errs[pl] = str(e) if (s0[1] == np.random.get_state()[1]).all() else 'np touched'
    c.check('activation/pleiotropy-rejected', all(v and v != 'np touched' for v in pleio_errs.values()), 'full/f/var/variable/v all raise ValueError, np.random untouched')
    ok_no = []
    for pl in ('no', 'n'):
        fa = c.factory(freq=0.0, loci=(3, 2), extra=[('quanti_pleiotropy', pl)])
        ok_no.append(fa._use_quanti_init_freq() and (c.bits(fa.make_founder('female')) == 1).all())
    c.check('activation/pleiotropy-no-accepted', all(ok_no), "quanti_pleiotropy='no'/'n' accepted and redrawn")

def g_lifecycle(c):
    ''' Nemo inheritance after a redrawn founder (mutation 0): clone / copy identical, offspring Mendelian, phenotypes recomputed, offspring not redrawn '''
    np = c.np
    for model in ('diallelic', 'diallelic_asymmetrical'):
        fac = c.factory(freq=[0.5] * 5, model=model, loci=(3, 2))
        np.random.seed(21)
        moms = [fac.make_founder('female') for _ in range(20)]; dads = [fac.make_founder('male') for _ in range(20)]
        ok_clone = all((c.bits(fac.make_clone(m)) == c.bits(m)).all() for m in moms)
        ok_copy = all((c.bits(m.make_copy()) == c.bits(m)).all() for m in moms)
        c.check('lifecycle/%s/clone-copy' % model, ok_clone and ok_copy, 'make_clone and make_copy keep all 2L alleles (mutation rate 0)')
        conv = {'s0 from mother': 0, 's0 from father': 0, 'no consistent origin': 0}
        pheno_ok, n_off = True, 0
        state = np.random.get_state()
        for i in range(400):
            m, f = moms[i % 20], dads[(i * 7) % 20]
            o = fac.make_offspring(m, f, 'female' if i % 2 else 'male'); n_off += 1
            ob, mb, fb = c.bits(o), c.bits(m), c.bits(f)
            in_m = lambda s: ((ob[s] == mb[0]) | (ob[s] == mb[1])).all()
            in_f = lambda s: ((ob[s] == fb[0]) | (ob[s] == fb[1])).all()
            a, b = in_m(0) and in_f(1), in_f(0) and in_m(1)
            conv['s0 from mother'] += int(a); conv['s0 from father'] += int(b); conv['no consistent origin'] += int(not (a or b))
            if not np.allclose(o.get_indi_phenotype_ls(), c.genotype(o), atol=1e-12): pheno_ok = False
        after = np.random.get_state()
        c.check('lifecycle/%s/mendelian' % model, conv['no consistent origin'] == 0, '%d offspring of random-p parents, consistent with: %s' % (n_off, conv))
        # homozygous parents (test-only: bits set by hand on two redrawn founders): every offspring must be heterozygous at every locus
        mom, dad = fac.make_founder('female'), fac.make_founder('male')
        for ind, bit in ((mom, False), (dad, True)):
            t = ind.get_trait('quant')
            for pos in c.positions(t, 2):
                for s in (0, 1): t.set_allele_bit(pos, s, bit)
            t.set_value()
        kids = [c.bits(fac.make_offspring(mom, dad, 'female' if i % 2 else 'male')) for i in range(200)]
        s_mom = {tuple(int(x) for x in k[0]) for k in kids}, {tuple(int(x) for x in k[1]) for k in kids}
        het = all(((k[0] + k[1]) == 1).all() for k in kids)
        c.check('lifecycle/%s/mendelian-homozygous-parents' % model, het, 'mother all allele 0 x father all allele 1 -> 200/200 offspring heterozygous at every locus; strand 0 patterns %s, strand 1 patterns %s' % (sorted(s_mom[0]), sorted(s_mom[1])))
        c.check('lifecycle/%s/offspring-phenotype' % model, pheno_ok, 'offspring phenotype == sum of its allele values (Ve 0)')
        c.check('lifecycle/%s/offspring-not-redrawn' % model, (state[1] == after[1]).all() and state[2] == after[2], 'make_offspring/make_clone draw no np.random (only founders are redrawn)')
    # founders' phenotype reaches selection: get_indi_phenotype_ls / phenotype_set agree
    fac = c.factory(freq=0.0, loci=(3, 2))
    ind = fac.make_founder('male')
    c.check('lifecycle/phenotype-views', np.allclose(list(ind.phenotype_set.values()), ind.get_indi_phenotype_ls()) and np.allclose(ind.get_indi_phenotype_ls(), [-0.6, -0.4]),
            'phenotype_set == get_indi_phenotype_ls == (-0.6, -0.4) for p=0')

def g_genetics(c):
    ''' non-additive genotype-phenotype maps: the phenotype set_value() leaves on a redrawn founder equals the one Nemo computes at birth for the same alleles
        (a zero-mutation clone: makeClonedOffsprg recomputes it), so the redraw is right whatever Nemo's G-P map '''
    np = c.np
    cases = (('dominance k=0.3, diallelic (3, 2)', dict(model='diallelic', loci=(3, 2), extra=[('quanti_dominance_effects', 0.3)])),
             ('dominance k=1.5, asymmetrical (2, 2)', dict(model='diallelic_asymmetrical', loci=(2, 2), extra=[('quanti_dominance_effects', 1.5)])),
             ('epistasis dense eps=0.5, diallelic (5,)', dict(model='diallelic', loci=(5,), extra=[('quanti_epistatic_effects', '{{0.5}}')])),
             ('dominance + Ve 0.01 (clone redraws Ve: compare genotype part only)', dict(model='diallelic', loci=(3, 2), ve=0.01, extra=[('quanti_dominance_effects', 0.3)])))
    for name, kw in cases:
        fac = c.factory(freq=0.5, **kw)
        np.random.seed(61)
        F = [fac.make_founder('female') for _ in range(60)]
        if kw.get('ve') is None:
            diffs = [max(abs(a - b) for a, b in zip(f.get_indi_phenotype_ls(), fac.make_clone(f).get_indi_phenotype_ls())) for f in F]
            nonadd = sum(1 for f in F if not np.allclose(f.get_indi_phenotype_ls(), c.genotype(f), atol=1e-9))
            c.check('genetics/%s' % name, max(diffs) < 1e-12 and nonadd > 0, '60 founders: |founder - Nemo clone| max %.2e; %d/60 differ from the additive sum (the map is active)' % (max(diffs), nonadd))
        else:
            P = np.array([f.get_indi_phenotype_ls() for f in F]); Q = np.array([fac.make_clone(f).get_indi_phenotype_ls() for f in F])
            d = (P - Q).ravel()
            c.check('genetics/%s' % name, abs(d.mean()) < 4 * np.sqrt(0.02 / d.size) and abs(d.std() / np.sqrt(0.02) - 1) < 0.25,
                    'founder - clone = difference of two independent Ve draws: mean %.4f sd %.4f (expect 0, %.4f)' % (d.mean(), d.std(), np.sqrt(0.02)))

def g_recycle(c):
    ''' founders built on recycled Nemo objects are fully redrawn; counters balance '''
    np = c.np
    fac = c.factory(freq=[0.3, 0.7, 0.5, 0.1, 0.9], loci=(3, 2), ve=0.01)
    np.random.seed(1)
    first = [fac.make_founder('female') for _ in range(300)]
    for ind in first: fac.recycle(ind)
    np.random.seed(77)
    got = np.array([c.bits(fac.make_founder('male')) for _ in range(300)])
    np.random.seed(77)
    want = np.array([(np.random.random((2, 5)) >= np.array([0.3, 0.7, 0.5, 0.1, 0.9])).astype(int) for _ in range(300)])
    c.check('recycle/redraw-on-reused-objects', (got == want).all(), '300 founders after 300 recycles == reference stream (no stale alleles)')
    c.check('recycle/counters', fac.minted_num == 600 and fac.recycled_num == 300, 'minted=%d recycled=%d' % (fac.minted_num, fac.recycled_num))
    fac2 = c.factory(freq=0.0, loci=(3, 2))
    olds = [fac2.make_founder('female') for _ in range(50)]
    for o in olds: fac2.recycle(o)
    news = [fac2.make_founder('female') for _ in range(50)]
    c.check('recycle/phenotype-fresh', all(np.allclose(n.get_indi_phenotype_ls(), [-0.6, -0.4]) for n in news), 'reused objects: phenotype recomputed (-0.6, -0.4)')
    fac2.purge_recycling_pool()

def g_multi_species(c):
    ''' species with and without quanti_init_freq interleaved: each draws only its own stream, caches are per species '''
    np = c.np
    A = c.factory(freq=[0.2, 0.8, 0.5, 0.5, 0.1], loci=(3, 2))
    B = c.factory(freq=0.6, loci=(4,), model='diallelic_asymmetrical')
    C = c.factory(loci=(3, 2), model='continuous', extra=[('quanti_init_model', 2), ('quanti_init_variance', 0.01)])
    D = c.factory(freq=[0.5] * 40, loci=(20, 20))
    order = ['A', 'C', 'B', 'D', 'A', 'B', 'C', 'C', 'D', 'A'] * 10
    np.random.seed(123)
    got = {'A': [], 'B': [], 'D': []}
    for s in order:
        ind = {'A': A, 'B': B, 'C': C, 'D': D}[s].make_founder('female')
        if s != 'C': got[s].append(c.bits(ind))
    np.random.seed(123)
    want = {'A': [], 'B': [], 'D': []}
    P = {'A': (np.array([0.2, 0.8, 0.5, 0.5, 0.1]), 5), 'B': (0.6, 4), 'D': (np.full(40, 0.5), 40)}
    for s in order:
        if s != 'C': want[s].append((np.random.random((2, P[s][1])) >= P[s][0]).astype(int))
    ok = all((np.array(got[s]) == np.array(want[s])).all() for s in got)
    c.check('multi_species/interleaved-streams', ok, '4 species (3 redrawn, 1 continuous) interleaved 100 founders: each matches its slice of one np stream')
    c.check('multi_species/caches', A._quanti_seq_pos == list(range(5)) and B._quanti_seq_pos == list(range(4)) and D._quanti_seq_pos == list(range(40)) and C._quanti_seq_pos is None,
            'per-species _quanti_seq_pos lengths 5/4/40, continuous stays None')

LCE_CASES = [  # (model, loci, allele_value, freq pattern) -- p in {0, 1}: Nemo's both branches are deterministic
    ('diallelic', (3, 2), '{{0.1}}', [1, 0, 0, 1, 0]),
    ('diallelic', (2, 3), '{{0.1, 0.2, 0.3, 0.4, 0.5}}', [0, 1, 1, 0, 1]),
    ('diallelic_asymmetrical', (3, 2), '{{0.1}}', [1, 0, 0, 1, 0]),
    ('diallelic_asymmetrical', (4,), '{{0.05, 0.1, 0.2, 0.4}}', [0, 0, 1, 1]),
    ('diallelic_HC', (5,), '{{0.1}}', [0, 0, 0, 0, 0]),
]

def g_nemo_lce(c):
    ''' Nemo's own LCE quanti_init (patch of 2 and of 8: Bernoulli and exact-count branches) -> bits and phenotypes, as data '''
    np, nemo = c.np, c.nemo
    for i, (model, loci, av, pat) in enumerate(LCE_CASES):
        for cap in (2, 8):
            p = c.params(model=model, loci=loci, allele_value=av, lce=('{{%s}}' % ', '.join(str(x) for x in pat), cap))
            sim = c.NemoSim(p, seed=1)
            for _ in sim.step_through(1, keep_open=True):
                pop = sim._pop; ti = pop.getTraitIndex('quant')
                inds = [pop.get(sx, ax, j, 0) for sx in (nemo.FEM, nemo.MAL) for ax in (nemo.ADLTx, nemo.OFFSx) for j in range(pop.size(sx, ax))]
                out = []
                for ind in inds:
                    t = ind.getTrait(ti); P = c.positions(t, len(loci))
                    out.append({'bits': [[t.get_allele(pos, s) for pos in P] for s in (0, 1)], 'pheno': [t.get_phenotype(d) for d in range(len(loci))]})
                c.data('lce/%d/cap%d' % (i, cap), out)
    # Nemo exact allocation vs MetaIBM Bernoulli at p = 0.3 (documented difference, D7)
    p = c.params(model='diallelic', loci=(5,), lce=('{{0.3}}', 2000))
    sim = c.NemoSim(p, seed=1)
    for _ in sim.step_through(1, keep_open=True):
        pop = sim._pop; ti = pop.getTraitIndex('quant')
        inds = [pop.get(sx, ax, j, 0) for sx in (nemo.FEM, nemo.MAL) for ax in (nemo.ADLTx, nemo.OFFSx) for j in range(pop.size(sx, ax))]
        B = np.array([[[ind.getTrait(ti).get_allele(pos, s) for pos in range(5)] for s in (0, 1)] for ind in inds])
        c.data('lce/p0.3/freq', {'n': len(inds), 'freq_allele0': (1 - B).mean(axis=(0, 1)).tolist()})

def g_metaibm_lce(c):
    ''' MetaIBM founders for the same LCE_CASES (compared with g_nemo_lce by the parent) + the p=0.3 Bernoulli frequency '''
    np = c.np
    for i, (model, loci, av, pat) in enumerate(LCE_CASES):
        fac = c.factory(freq=[float(x) for x in pat], model=model, loci=loci, allele_value=av)
        out = []
        for g in ('female', 'male'):
            ind = fac.make_founder(g)
            out.append({'bits': c.bits(ind).tolist(), 'pheno': ind.get_indi_phenotype_ls()})
        c.data('mib/%d' % i, out)
    fac = c.factory(freq=0.3, model='diallelic', loci=(5,))
    np.random.seed(9)
    B = np.array([c.bits(fac.make_founder('female' if i % 2 else 'male')) for i in range(2000)])
    c.data('mib/p0.3/freq', {'n': 2000, 'freq_allele0': (1 - B).mean(axis=(0, 1)).tolist()})

def g_regression(c):
    ''' without quanti_init_freq the working tree must equal HEAD bit for bit (parent runs this group on both and compares) '''
    np, nemo = c.np, c.nemo
    out = {}
    cases = (('continuous-example', dict(model='continuous', loci=(20, 20), ve='{{0.0025, 0.0025}}', mut=0.0001,
                                         extra=[('quanti_init_model', 2), ('quanti_init_value', '{{0.5, 0.5}}'), ('quanti_init_variance', '{{0.01, 0.01}}')])),
             ('diallelic-default', dict(model='diallelic', loci=(3, 2), ve=0.01, mut=0.001)),
             ('asym-default', dict(model='diallelic_asymmetrical', loci=(4,), mut=0.001)))
    for name, kw in cases:
        fac = c.factory(seed=31, **kw)
        np.random.seed(31); random.seed(31)
        F = [fac.make_founder('female' if i % 2 else 'male', home_habitat=i % 3) for i in range(40)]
        moms = [x for x in F if x.gender == 'female']; dads = [x for x in F if x.gender == 'male']
        O = [fac.make_offspring(moms[i % len(moms)], dads[(i * 3) % len(dads)], 'female' if i % 2 else 'male') for i in range(30)]
        K = [fac.make_clone(moms[i % len(moms)]) for i in range(10)]
        Y = [F[i].make_copy() for i in range(5)]
        rows = [[i.individual_id, i.gender, c.values(i).tolist(), [repr(x) for x in i.get_indi_phenotype_ls()]] for i in F + O + K + Y]
        out[name] = hashlib.md5(json.dumps(rows).encode()).hexdigest()
        out[name + '/np'] = hashlib.md5(np.random.get_state()[1].tobytes()).hexdigest()
        out[name + '/random'] = hashlib.md5(repr(random.getstate()).encode()).hexdigest()
    c.data('regression', out)

def g_user_notes(c):
    ''' combinations left to the user (D3, D4, D7): behaviour recorded as NOTE, not judged '''
    np = c.np
    # D3: quanti_init_model != 0 with quanti_init_freq -- the redraw overwrites every allele
    res = {}
    for im in (0, 1, 2, 3):
        try:
            fac = c.factory(freq=[0.3, 0.6, 0.5, 0.8, 0.2], loci=(3, 2), extra=[('quanti_init_model', im)])
            np.random.seed(14)
            res[im] = [c.bits(fac.make_founder('female')).tolist() for _ in range(20)]
        except Exception as e: res[im] = '%s: %s' % (type(e).__name__, e)
    same = all(res[im] == res[0] for im in res)
    c.check('user/D3-init-model', same, 'quanti_init_model 0/1/2/3 + quanti_init_freq: founders bit-identical for one np seed -> the redraw overwrites every allele' if same else str({k: (v if isinstance(v, str) else 'ok') for k, v in res.items()}), note=True)
    # D4: diallelic_asymmetrical with a scalar quanti_allele_value is silently symmetric
    ph = {}
    for av in ('{{0.1}}', '0.1'):
        fac = c.factory(freq=0.0, loci=(3, 2), model='diallelic_asymmetrical', allele_value=av)
        ph[av] = fac.make_founder('female').get_indi_phenotype_ls()
    c.check('user/D4-asym-scalar', True, 'p=0 (all allele 1): matrix {{0.1}} -> %s, scalar 0.1 -> %s (scalar is read as symmetric +-a)' % (np.round(ph['{{0.1}}'], 6).tolist(), np.round(ph['0.1'], 6).tolist()), note=True)
    # D7: full pleiotropy (not supported this round) -- what the four forms do
    msgs = []
    for name, f in (('0.8', 0.8), ('[0.8]', [0.8]), ('[[0.8]]', [[0.8]]), ('[0.8]*num_locus(3)', [0.8] * 3), ('[0.8]*seq_len(6)', [0.8] * 6)):
        try:
            fac = c.factory(freq=f, loci=(3, 3), pleio='full')
            np.random.seed(2)
            ind = fac.make_founder('female')
            msgs.append('%s: ok, %d positions written' % (name, len(fac._quanti_seq_pos)))
        except Exception as e: msgs.append('%s: %s' % (name, type(e).__name__))
    c.check('user/D7-full-pleiotropy', True, '2 traits x 3 loci: ' + '; '.join(msgs), note=True)

def g_init_value_note(c):
    ''' quanti_init_value set on a di-allelic species (Nemo prints ***ERROR*** at boot, non-fatal); the redraw still applies '''
    fac = c.factory(freq=0.0, loci=(3, 2), extra=[('quanti_init_value', '{{0.5, 0.5}}')])
    b = c.bits(fac.make_founder('female'))
    c.check('user/init-value-diallelic', (b == 1).all(), 'founder still redrawn (all allele 1); Nemo boot messages: see parent detail', note=True)

def g_perf(c):
    ''' cost of the redraw: 2L set_allele_bit SWIG calls + set_value per founder '''
    np = c.np
    rows = []
    for name, kw in (('diallelic L=40, no freq', dict(loci=(20, 20))), ('diallelic L=40, freq [p]*40', dict(loci=(20, 20), freq=[0.5] * 40)),
                     ('diallelic L=40, freq 0.5', dict(loci=(20, 20), freq=0.5)), ('diallelic L=200, freq 0.5', dict(loci=(100, 100), freq=0.5)),
                     ('continuous L=40 (example-like)', dict(loci=(20, 20), model='continuous', extra=[('quanti_init_model', 2), ('quanti_init_variance', 0.01)]))):
        fac = c.factory(**kw)
        fac.make_founder('female')
        n = 2000; t0 = time.perf_counter()
        for _ in range(n): fac.recycle(fac.make_founder('female'))
        rows.append('%s: %.1f us/founder' % (name, (time.perf_counter() - t0) / n * 1e6))
    c.check('perf/make_founder', True, '; '.join(rows), note=True)

##*********** TTrait.set_allele_bit on its own (absorbed from Claude outputs/test_set_allele_bit.py) ***********
# A running 1-patch, 1-trait, 10-locus Nemo simulation (breed, no mutation); at generation 2 the first adult's alleles are edited.
# One simulation per child process: each case is its own group.
SAB_LOCI = 10

def _sab_run(c, model, allele_value=None, ntrl=False, body=None):
    p = c.SimParameters()
    for k, v in [('run_mode', 'silent'), ('random_seed', 42), ('root_dir', c.cache), ('logfile', c.cache + '/nemo.log'), ('filename', 'bit_test'),
                 ('replicates', 1), ('generations', 3), ('patch_number', 1), ('patch_capacity', 20),
                 ('breed', 1), ('mating_system', 1), ('mating_isWrightFisher', None), ('quanti_traits', 1), ('quanti_loci', SAB_LOCI),
                 ('quanti_mutation_rate', 0), ('quanti_allele_model', model), ('quanti_init_model', 0)]:
        p.add(k, v)
    if allele_value is not None: p.add('quanti_allele_value', allele_value)
    else: p.add('quanti_init_value', '{{2.4}}'); p.add('quanti_mutation_variance', 0.1)
    if ntrl:
        for k, v in [('ntrl_loci', 5), ('ntrl_all', 2), ('ntrl_mutation_rate', 0), ('ntrl_mutation_model', 1)]: p.add(k, v)
    sim, out = c.NemoSim(p), None
    for gen in sim.step_through(3):
        if gen == 2:
            ind = sim.population[0].adults[0]
            t = ind.trait(0)
            out = body(ind, t, t._t)
    return out

def _sab_edit(ind, t, tt):
    ''' strand 0 set to allele 1 at every locus, then set_value() '''
    pos = [tt.get_locus_seq_pos(l, 0) for l in range(SAB_LOCI)]
    r = {'before_bits': [tt.get_allele(p, s) for s in (0, 1) for p in pos], 'before_vals': [tt.get_allele_value(p, s) for s in (0, 1) for p in pos], 'before_pheno': t.phenotype(0)}
    for p in pos: tt.set_allele_bit(p, 0, True)
    t.set_value()
    r.update({'after_bits': [tt.get_allele(p, s) for s in (0, 1) for p in pos], 'after_vals': [tt.get_allele_value(p, s) for s in (0, 1) for p in pos], 'after_pheno': t.phenotype(0)})
    return r

def g_sab_symmetric(c):
    r = _sab_run(c, 'diallelic', 0.5, body=_sab_edit)
    L = SAB_LOCI
    ok = r['before_bits'] == [0] * (2 * L) and r['after_bits'] == [1] * L + [0] * L and r['after_vals'] == [-0.5] * L + [0.5] * L and abs(r['before_pheno'] - 2 * L * 0.5) < 1e-9 and abs(r['after_pheno']) < 1e-9
    c.check('set_allele_bit/symmetric-diallelic', ok, 'strand 0 -> allele 1: bits read back, +0.5 -> -0.5 on strand 0 only, phenotype %g -> %g (expect 10 -> 0)' % (r['before_pheno'], r['after_pheno']))

def g_sab_asymmetrical(c):
    r = _sab_run(c, 'diallelic_asymmetrical', '{{0.5}}', body=_sab_edit)
    L = SAB_LOCI
    ok = r['after_vals'] == [0.5] * L + [0.0] * L and abs(r['before_pheno']) < 1e-9 and abs(r['after_pheno'] - L * 0.5) < 1e-9
    c.check('set_allele_bit/asymmetrical-diallelic', ok, 'allele 0 = 0, allele 1 = a: strand 0 -> allele 1 gives phenotype %g (expect %g)' % (r['after_pheno'], L * 0.5))

def g_sab_continuous(c):
    r = _sab_run(c, 'continuous', body=_sab_edit)
    ok = r['after_vals'] == r['before_vals'] and abs(r['after_pheno'] - r['before_pheno']) < 1e-9
    c.check('set_allele_bit/continuous-noop', ok, 'continuous model: TTQuanti_continuous::set_allele_bit is empty -> values and phenotype unchanged')

def g_sab_ntrl(c):
    def body(ind, t, tt):
        before = [tt.get_allele(l, s) for s in (0, 1) for l in range(5)]
        for l in range(5):
            for s in (0, 1): tt.set_allele_bit(l, s, not tt.get_allele(l, s))
        return {'before': before, 'after': [tt.get_allele(l, s) for s in (0, 1) for l in range(5)], 'quant': ind.trait(1).phenotype(0)}
    r = _sab_run(c, 'diallelic', 0.5, ntrl=True, body=body)
    ok = r['after'] == r['before'] and abs(r['quant'] - 2 * SAB_LOCI * 0.5) < 1e-9
    c.check('set_allele_bit/non-quanti-noop', ok, 'ntrl trait (slot 0): every allele written with its opposite -> no crash, nothing changes; quant phenotype %g (expect 10)' % r['quant'])

##****************************** parent side ******************************
def run_child(group, root=PROJECT_ROOT, py=sys.executable):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    t0 = time.time()
    r = subprocess.run([py, '-I', '-B', str(Path(__file__).resolve()), '--child', group, str(root)], capture_output=True, text=True, env=env, cwd=tempfile.gettempdir())
    raw = r.stdout + r.stderr
    checks = [json.loads(l[3:]) for l in raw.splitlines() if l.startswith('@@ ')]
    data = {d['key']: d['value'] for d in (json.loads(l[7:]) for l in raw.splitlines() if l.startswith('@@DATA '))}
    if r.returncode != 0: checks.append({'group': group, 'id': '%s/child-exit' % group, 'status': 'FAIL', 'detail': 'exit %d: %s' % (r.returncode, raw[-800:])})
    return checks, data, raw, time.time() - t0

def head_copy():
    ''' metaibm/ + extension/ of the working tree with HEAD's indfactory.py, in a temp dir '''
    d = Path(tempfile.mkdtemp(prefix='qif_head_'))
    for pkg in ('metaibm', 'extension'): shutil.copytree(PROJECT_ROOT / pkg, d / pkg, ignore=shutil.ignore_patterns('__pycache__'))
    head = subprocess.run(['git', 'show', 'HEAD:metaibm/indfactory.py'], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True).stdout
    (d / 'metaibm' / 'indfactory.py').write_text(head)
    return d

def main(groups):
    md5 = hashlib.md5((PROJECT_ROOT / 'metaibm' / 'indfactory.py').read_bytes()).hexdigest()
    import pynemo
    print('python %s | pynemo %s | indfactory.py md5 %s' % (sys.executable, os.path.dirname(pynemo.__file__), md5))
    results, datas, t_all = [], {}, time.time()
    for g in groups:
        if g in ('nemo_lce', 'metaibm_lce', 'regression', 'init_value_note'): continue
        checks, data, raw, dt = run_child(g); results += checks; datas.update(data)
        print('  %-14s %5.1fs  %s' % (g, dt, ' '.join(x['status'][0] for x in checks)), flush=True)
    if 'nemo_lce' in groups or 'metaibm_lce' in groups:
        cn, dn, _, t1 = run_child('nemo_lce'); cm, dm, _, t2 = run_child('metaibm_lce'); results += cn + cm
        for i, (model, loci, av, pat) in enumerate(LCE_CASES):
            for cap in (2, 8):
                nemo_inds = dn.get('lce/%d/cap%d' % (i, cap), []); mib = dm.get('mib/%d' % i, [])
                ok = bool(nemo_inds) and bool(mib) and all(x['bits'] == mib[0]['bits'] and all(abs(a - b) < 1e-12 for a, b in zip(x['pheno'], mib[0]['pheno'])) for x in nemo_inds + mib)
                results.append({'group': 'nemo_lce', 'id': 'nemo_lce/%s/loci%s/%s/cap%d' % (model, loci, pat, cap), 'status': 'PASS' if ok else 'FAIL',
                                'detail': 'Nemo LCE (%d inds) == MetaIBM founders: bits %s, phenotype %s' % (len(nemo_inds), mib[0]['bits'] if mib else None, [round(v, 6) for v in mib[0]['pheno']] if mib else None)})
        fn, fm = dn.get('lce/p0.3/freq'), dm.get('mib/p0.3/freq')
        results.append({'group': 'nemo_lce', 'id': 'nemo_lce/p0.3-exact-vs-bernoulli', 'status': 'NOTE',
                        'detail': 'Nemo LCE (N=%d, exact counts) allele-0 freq=%s; MetaIBM (N=%d, Bernoulli) =%s' % (fn['n'], [round(x, 4) for x in fn['freq_allele0']], fm['n'], [round(x, 4) for x in fm['freq_allele0']]) if fn and fm else 'missing data'})
        print('  %-14s %5.1fs' % ('nemo/metaibm lce', t1 + t2), flush=True)
    if 'regression' in groups:
        hd = head_copy()
        try:
            ch, dh, _, t1 = run_child('regression', root=hd); cw, dw, _, t2 = run_child('regression'); results += ch + cw
            h, w = dh.get('regression', {}), dw.get('regression', {})
            for k in sorted(h):
                results.append({'group': 'regression', 'id': 'regression/HEAD==tree/%s' % k, 'status': 'PASS' if h.get(k) == w.get(k) else 'FAIL', 'detail': 'HEAD %s | tree %s' % (h.get(k, '')[:12], w.get(k, '')[:12])})
            if not h: results.append({'group': 'regression', 'id': 'regression/data', 'status': 'FAIL', 'detail': 'no data'})
        finally: shutil.rmtree(hd, ignore_errors=True)
        print('  %-14s %5.1fs' % ('regression', t1 + t2), flush=True)
    if 'init_value_note' in groups:
        ci, _, raw, dt = run_child('init_value_note')
        errs = sorted(set(l.strip() for l in raw.splitlines() if 'ERROR' in l))
        for x in ci: x['detail'] += ' | Nemo printed: %s' % (errs or 'nothing')
        results += ci
        print('  %-14s %5.1fs' % ('init_value', dt), flush=True)
    if 'rng' in groups:     # cross-process reproducibility
        _, d1, _, _ = run_child('rng'); _, d2, _, _ = run_child('rng')
        results.append({'group': 'rng', 'id': 'rng/cross-process-repro', 'status': 'PASS' if d1.get('rng/digest') and d1.get('rng/digest') == d2.get('rng/digest') else 'FAIL',
                        'detail': 'two child processes: md5 %s / %s' % (str(d1.get('rng/digest'))[:12], str(d2.get('rng/digest'))[:12])})
    n = {s: sum(1 for x in results if x['status'] == s) for s in ('PASS', 'FAIL', 'NOTE')}
    print('\n| status | id | detail |\n|---|---|---|')
    for x in results: print('| %s | `%s` | %s |' % (x['status'], x['id'], x['detail'].replace('|', '\\|')))
    print('\nPASS %d  FAIL %d  NOTE %d  (%.0fs)  indfactory.py md5 %s' % (n['PASS'], n['FAIL'], n['NOTE'], time.time() - t_all, md5))
    return 1 if n['FAIL'] else 0

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--child':
        ctx = Ctx(sys.argv[2], Path(sys.argv[3]))
        try: globals()['g_' + sys.argv[2]](ctx)
        finally: ctx.close()
    else:
        sys.exit(main([g for g in sys.argv[1:]] or list(GROUPS)))
