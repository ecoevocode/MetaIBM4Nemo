# -*- coding: utf-8 -*-
''' test_trait_readers.py -- individual.get_trait_value(trait_key, **kwargs) and individual.get_trait_alleles(trait_key, method, **kwargs)

    Run (pynemo with TTrait.get_value / get_wolb_value -> conda env metaibm-nemo-test6):
        PYTHONDONTWRITEBYTECODE=1 ~/miniconda3/envs/metaibm-nemo-test6/bin/python -B test_trait_readers.py [group ...]

    Every group runs in its own child process (Nemo's RAND, Individual::currentID, the trait prototypes and the
    static genetic map are process-global); inside one process every species keeps the same loci counts and layout.
    A child prints one '@@ {json}' line per check: status PASS | FAIL | NOTE (NOTE = behaviour recorded, not judged).
    Nemo writes only into temporary directories; metaibm/ is imported read-only (no __pycache__).

    Oracles (test-only):
        raw      direct pynemo calls per sequence position (TTrait.get_allele / get_allele_value), locus counts from nemo_params
        stored   the value at the TTrait::getValue() address, read with ctypes (c_double; c_bool for wolb)
        formula  delet / dmi fitness and disp rate recomputed from the allele tables (Nemo source: ttdeletmutations_bitstring.cc
                 viability_multi / viability_epist, ttbdmi.cc viability_diplo + the 9 -> 16 table expansion, ttdispersal.h set_value)
        legacy   get_indi_phenotype_ls(), phenotype_set, TTrait.get_genotypic_value(d)
'''
import sys
sys.dont_write_bytecode = True                  # never drop __pycache__ into the user's metaibm/
import os, json, time, shutil, hashlib, tempfile, subprocess, random, ctypes, math
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
PY_ALT = Path.home() / 'miniconda3' / 'envs' / 'metaibm-nemo-test5' / 'bin' / 'python'     # pynemo without get_value (regression_env)

KEYS = ('quant', 'ntrl', 'delet', 'dmi', 'fdisp', 'mdisp', 'wolb')
DBL = ('delet', 'dmi', 'fdisp', 'mdisp')
VALID_PAIRS = (('quant', 'get_allele_value'), ('quant', 'get_allele'), ('ntrl', 'get_allele'), ('delet', 'get_allele'), ('dmi', 'get_allele'),
               ('fdisp', 'get_allele_value'), ('mdisp', 'get_allele_value'), ('wolb', 'get_allele'))
QUANT_CASES = ('continuous', 'diallelic', 'asymmetrical', 'pleiotropy', 'per_trait', 'three_axes', 'mapped')
OOR_CASES = ('value_dim2', 'value_dim-1', 'alleles_values_dim2', 'alleles_states_dim2', 'alleles_dim-1')
GROUPS = (('routing', 'kwargs') + tuple('quant_layout:' + q for q in QUANT_CASES) +
          ('quant_genetics', 'ntrl', 'delet', 'dmi', 'disp', 'wolb', 'stored', 'lifecycle', 'out_of_range', 'regression_head', 'regression_env', 'perf'))

DMI16 = [[round(-0.02 * (j + 1) - 0.005 * r, 6) for j in range(16)] for r in range(3)]            # 3 distinct rows of 16 (4 pairs -> row 3 recycles row 0)
DMI9 = [[0.0, -0.05, -0.1, -0.15, -0.2, -0.25, -0.3, -0.35, -0.4]]                                  # 9-column form, one row for every pair
DMI9_INDEX = (0, 1, 1, 2, 3, 4, 4, 5, 3, 4, 4, 5, 6, 7, 7, 8)                                       # ttbdmi.cc: 9 -> 16 expansion

##****************************** child side: harness ******************************
def _import_stack(root):
    sys.path.insert(0, str(root))
    import numpy as np
    from pynemo import _nemo, nemo
    from pynemo.params import SimParameters
    from metaibm.indfactory import indfactory
    _nemo.cvar.SILENT_RUN = True
    return np, nemo, SimParameters, indfactory

class Ctx:
    ''' one child process: imports, a temporary Nemo cache, check() / data() reporting, species builders '''
    def __init__(self, group, root):
        self.group = group
        self.np, self.nemo, self.SimParameters, self.indfactory = _import_stack(root)
        self.cache = tempfile.mkdtemp(prefix='ntr_%s_' % group.replace(':', '_'))
        self.n = 0

    def check(self, cid, ok, detail='', note=False):
        status = 'NOTE' if note else ('PASS' if ok else 'FAIL')
        print('@@ ' + json.dumps({'group': self.group, 'id': cid, 'status': status, 'detail': str(detail)}), flush=True)
        return ok

    def data(self, key, value):
        print('@@DATA ' + json.dumps({'key': key, 'value': value}), flush=True)

    def close(self):
        shutil.rmtree(self.cache, ignore_errors=True)

    def species(self, traits, pheno=('t0', 't1'), seed=1, mnp=None):
        ''' one indfactory; traits = {'quant' | 'ntrl' | 'delet' | 'dmi' | 'disp' | 'wolb': rows}; 'disp' builds fdisp + mdisp '''
        self.n += 1
        p = self.SimParameters()
        rows = [('replicates', 1), ('generations', 1), ('patch_number', 1), ('patch_capacity', 1), ('aging', 1), ('run_mode', 'silent'),
                ('root_dir', self.cache), ('logfile', self.cache + '/nemo.log'), ('filename', 'sp%d' % self.n)]
        for r in traits.values(): rows += list(r)
        for k, v in rows: p.add(k, v)
        keys = tuple(k for k in traits if k != 'disp') + (('fdisp', 'mdisp') if 'disp' in traits else ())
        return self.indfactory('sp%d' % self.n, list(pheno) if 'quant' in traits else [], p, trait_keys=keys, seed=seed, metaibm_nemo_params=mnp)

##****************************** trait parameter rows (every species of a process uses the same loci layout) ******************************
def quant_rows(case='continuous', mut=0.0, ve=0.0, extra=()):
    if case == 'three_axes':
        return [('quanti_traits', 3), ('quanti_loci', 40), ('quanti_loci_per_trait', '{{14, 13, 13}}'), ('quanti_allele_model', 'continuous'), ('quanti_init_model', 2),
                ('quanti_init_value', '{{0.3, -0.2, 0.1}}'), ('quanti_init_variance', '{{0.002, 0.002, 0.002}}'), ('quanti_mutation_rate', mut),
                ('quanti_mutation_variance', 0.01), ('quanti_environmental_variance', '{{%r, %r, %r}}' % (ve, ve, ve))] + list(extra)
    rows = [('quanti_traits', 2), ('quanti_loci', 40), ('quanti_mutation_rate', mut), ('quanti_environmental_variance', '{{%r, %r}}' % (ve, ve))]
    cont = [('quanti_allele_model', 'continuous'), ('quanti_init_model', 2), ('quanti_init_value', '{{0.3, -0.2}}'), ('quanti_init_variance', '{{0.002, 0.002}}'), ('quanti_mutation_variance', 0.01)]
    if case == 'continuous':     rows += cont
    elif case == 'diallelic':    rows += [('quanti_allele_model', 'diallelic'), ('quanti_allele_value', 0.01)]
    elif case == 'asymmetrical': rows += [('quanti_allele_model', 'diallelic_asymmetrical'), ('quanti_allele_value', '{{0.01}}')]   # matrix form: {0, a}
    elif case == 'asymmetrical_scalar': rows += [('quanti_allele_model', 'diallelic_asymmetrical'), ('quanti_allele_value', 0.01)]  # scalar form (Nemo stores +a / -a)
    elif case == 'pleiotropy':   rows += cont + [('quanti_pleiotropy', 'full')]
    elif case == 'per_trait':    rows += cont + [('quanti_loci_per_trait', '{{30, 10}}')]
    elif case == 'mapped':       rows += cont + [('quanti_recombination_rate', 0.05)]
    else: raise KeyError(case)
    return rows + list(extra)

def ntrl_rows(n_all=4, mut=0.0): return [('ntrl_loci', 12), ('ntrl_all', n_all), ('ntrl_mutation_rate', mut), ('ntrl_mutation_model', 1)]
def delet_rows(mut=0.05, model=1, scale=None, init=0.3, s=0.05, h=0.3):
    rows = [('delet_loci', 30), ('delet_mutation_rate', mut), ('delet_effects_mean', s), ('delet_dominance_mean', h), ('delet_fitness_model', model), ('delet_init_freq', init)]
    return rows + ([('delet_fitness_scaling_factor', scale)] if scale is not None else [])
def table_str(t): return '{' + ''.join('{%s}' % ','.join(repr(x) for x in row) for row in t) + '}'
def dmi_rows(mut=0.1, table=None): return [('dmi_loci', 8), ('dmi_mutation_rate', mut), ('dmi_genot_table', table_str(table or DMI16))]
def disp_rows(mut=0.0, mean=0.1, fem=0.37, mal=0.81): return [('disp_mutation_rate', mut), ('disp_mutation_mean', mean), ('disp_init_rate_fem', fem), ('disp_init_rate_mal', mal)]
def wolb_rows(rate=1.0): return [('wolbachia_transmission_rate', rate)]
def all_traits(mut=1.0, ve=0.0):
    ''' the 7-trait species used by several groups; mut scales every mutation rate (0 = none) '''
    return {'quant': quant_rows('continuous', mut=0.01 * mut, ve=ve), 'ntrl': ntrl_rows(4, 0.05 * mut), 'delet': delet_rows(0.2 * mut), 'dmi': dmi_rows(0.2 * mut),
            'disp': disp_rows(0.2 * mut, 0.3), 'wolb': wolb_rows(1.0)}

##****************************** oracles ******************************
def n_loci(ind, key):
    ''' locus count straight from nemo_params (independent of the reader's getTraitPrototype().get_locus_number()) '''
    return int(ind.indfactory.nemo_params_dict[key + '_loci'])

def raw_table(np, ind, key, method, trait_dim=0):
    ''' (n, 2) table by direct TTrait calls; quant positions from Nemo's own locus table; wolb -> bool '''
    t = ind.get_trait(key)
    if key == 'wolb': return bool(t.get_allele(0, 0))
    if key == 'quant': pos = [t.get_locus_seq_pos(l, trait_dim) for l in range(t.get_num_loci(trait_dim))]
    elif key in ('fdisp', 'mdisp'): pos = [0]
    else: pos = range(n_loci(ind, key))
    f = getattr(t, method)
    return np.array([[f(q, 0), f(q, 1)] for q in pos], dtype=float if method == 'get_allele_value' else int).reshape(-1, 2)

def stored(ind, key, d=0):
    ''' value at the TTrait::getValue() address (sandbox oracle, read-only); None when Nemo stores nothing (ntrl) '''
    v = ind.get_trait(key).getValue()
    if v is None: return None
    if key == 'wolb': return bool(ctypes.c_bool.from_address(int(v)).value)
    return ctypes.c_double.from_address(int(v) + 8 * d).value

def delet_formula(st, s, h, model, scale=1.0):
    ''' W from the (n, 2) state table: model 1 (1-s)^hmz (1-hs)^htz * scale, model 2 1 - s(hmz + h htz) * scale; RETURN (W, hmz, htz) '''
    hmz = int(((st[:, 0] == 1) & (st[:, 1] == 1)).sum()); htz = int((st[:, 0] != st[:, 1]).sum())
    if model == 1: return math.pow(1.0 - s, hmz) * math.pow(1.0 - h * s, htz) * scale, hmz, htz
    return 1.0 - s * (hmz + h * htz) * scale, hmz, htz

def dmi_table16(table):
    rows = [list(r) for r in table]
    if len(rows[0]) == 9: rows = [[r[i] for i in DMI9_INDEX] for r in rows]
    return rows

def dmi_formula(st, table16):
    ''' W = prod over locus pairs (1 + T[pair][4 cA + cB]), c = 2 strand0 + strand1; rows recycled over pairs; RETURN (W, codes) '''
    W, codes = 1.0, []
    for pair in range(len(st) // 2):
        cA = 2 * int(st[2 * pair, 0]) + int(st[2 * pair, 1]); cB = 2 * int(st[2 * pair + 1, 0]) + int(st[2 * pair + 1, 1])
        W *= 1 + table16[pair % len(table16)][4 * cA + cB]; codes.append((pair, 4 * cA + cB))
    return W, codes

def raises(f, *a, **kw):
    ''' RETURN: 'ExcType: message' or None '''
    try: f(*a, **kw); return None
    except Exception as e: return '%s: %s' % (type(e).__name__, e)

def founders(sp, n=20):
    return [sp.make_founder('female' if j % 2 else 'male') for j in range(n)]

def generations(sp, pop, gens, size=40):
    ''' random mating from pop (mother = female); RETURN: every individual of every generation, plus the last generation '''
    every = list(pop)
    for _ in range(gens):
        F = [y for y in pop if y.gender == 'female']; M = [y for y in pop if y.gender == 'male']
        pop = [sp.make_offspring(random.choice(F), random.choice(M), 'female' if j % 2 else 'male') for j in range(size)]
        every += pop
    return every, pop

def orientation(np, kids, table):
    ''' per locus: column 1 from the mother and column 0 from the father? RETURN (mismatch rate as read, mismatch rate if the columns were swapped) '''
    ok = sw = n = 0
    for k, mo, fa in kids:
        K, Mo, Fa = table(k), table(mo), table(fa)
        for l in range(len(K)):
            n += 1; ok += not (K[l, 1] in Mo[l] and K[l, 0] in Fa[l]); sw += not (K[l, 0] in Mo[l] and K[l, 1] in Fa[l])
    return ok / n, sw / n

def dump_species(c):
    ''' the 7-trait species with every mutation on; RETURN: all individuals (founders + 3 generations + 100 kids) '''
    random.seed(11); c.np.random.seed(11)
    sp = c.species(all_traits(1.0, ve=0.5), seed=21)
    pop = founders(sp, 20)
    for j, y in enumerate(pop):
        if y.gender == 'female' and j % 4 == 1: y.get_trait('wolb').set_allele_value(0, 0, 1.0)
    every, last = generations(sp, pop, 3)
    F = [y for y in last if y.gender == 'female']; M = [y for y in last if y.gender == 'male']
    every += [sp.make_offspring(random.choice(F), random.choice(M), 'female' if j % 2 else 'male') for j in range(100)]
    return sp, every

##****************************** groups ******************************
def g_routing(c):
    np = c.np
    sp = c.species(all_traits(0.5), seed=3)
    pop = founders(sp, 10)
    pop[1].get_trait('wolb').set_allele_value(0, 0, 1.0)
    inds = pop + [sp.make_offspring(pop[1], pop[0], 'female') for _ in range(10)]
    y = inds[0]
    # get_trait_alleles: every valid pair
    for d in range(2):
        L = y.get_trait('quant').get_num_loci(d)
        v = [z.get_trait_alleles('quant', 'get_allele_value', trait_dim=d) for z in inds]; s = [z.get_trait_alleles('quant', 'get_allele', trait_dim=d) for z in inds]
        c.check('routing/alleles/quant/get_allele_value/axis%d' % d, all(isinstance(a, np.ndarray) and a.dtype == np.float64 and a.shape == (L, 2) for a in v), 'np.ndarray float64 shape (%d, 2)' % L)
        c.check('routing/alleles/quant/get_allele/axis%d' % d, all(isinstance(a, np.ndarray) and a.dtype.kind == 'i' and a.shape == (L, 2) for a in s), 'np.ndarray int shape (%d, 2)' % L)
    for k in ('ntrl', 'delet', 'dmi'):
        a = [z.get_trait_alleles(k, 'get_allele') for z in inds]
        c.check('routing/alleles/%s/get_allele' % k, all(isinstance(x, np.ndarray) and x.dtype.kind == 'i' and x.shape == (n_loci(y, k), 2) for x in a), 'np.ndarray int shape (%d, 2)' % n_loci(y, k))
    for k in ('fdisp', 'mdisp'):
        a = [z.get_trait_alleles(k, 'get_allele_value') for z in inds]
        c.check('routing/alleles/%s/get_allele_value' % k, all(isinstance(x, np.ndarray) and x.dtype == np.float64 and x.shape == (1, 2) for x in a), 'np.ndarray float64 shape (1, 2)')
    w = [z.get_trait_alleles('wolb', 'get_allele') for z in inds]
    c.check('routing/alleles/wolb/get_allele', all(type(x) is bool for x in w) and set(w) == {True, False}, 'Python bool, both values seen')
    # get_trait_value: every trait
    c.check('routing/value/quant', all(type(z.get_trait_value('quant', trait_dim=d)) is float for z in inds for d in range(2)), 'Python float per axis')
    for k in DBL: c.check('routing/value/%s' % k, all(type(z.get_trait_value(k)) is float for z in inds), 'Python float')
    c.check('routing/value/wolb', all(type(z.get_trait_value('wolb')) is bool for z in inds) and {z.get_trait_value('wolb') for z in inds} == {True, False}, 'Python bool, both values seen')
    c.check('routing/value/ntrl', all(z.get_trait_value('ntrl') is None for z in inds), 'None (native getValue() returns NULL)')
    # invalid pairs, unknown keys / methods
    invalid = [(k, m) for k in KEYS for m in ('get_allele', 'get_allele_value') if (k, m) not in VALID_PAIRS]
    errs = {p: raises(y.get_trait_alleles, *p) for p in invalid + [('ntral', 'get_allele'), ('quant', 'get_alleles'), ('quant', None), (None, 'get_allele'), ('Quant', 'get_allele_value')]}
    c.check('routing/alleles/invalid-pairs', len(invalid) == 6 and all(e and e.startswith('ValueError') and "'get_allele_value': quant/fdisp/mdisp; 'get_allele': quant/ntrl/delet/dmi/wolb" in e for e in errs.values()),
            '%d invalid pairs + 5 unknown key/method -> ValueError naming the valid pairs; e.g. %s' % (len(invalid), errs[('wolb', 'get_allele_value')]))
    errs = {k: raises(y.get_trait_value, k) for k in ('ntral', 'Quant', None, '', 'disp')}
    c.check('routing/value/unknown-key', all(e and e.startswith('ValueError') and "trait_key = 'quant'/'ntrl'/'delet'/'dmi'/'fdisp'/'mdisp'/'wolb'" in e for e in errs.values()),
            "5 unknown keys (incl. 'disp', '') -> ValueError listing the 7 keys; e.g. %s" % errs['disp'])
    # traits absent from a species (recorded, not judged: traits_index raises KeyError)
    q_only = c.species({'quant': quant_rows('continuous')}, seed=4); z = q_only.make_founder('female')
    c.check('routing/absent-trait/quant-only-species', True, 'species without delet: get_trait_value(delet) -> %s; get_trait_alleles(delet, get_allele) -> %s'
            % (raises(z.get_trait_value, 'delet'), raises(z.get_trait_alleles, 'delet', 'get_allele')), note=True)
    nq = c.species({'ntrl': ntrl_rows(4), 'delet': delet_rows()}, pheno=(), seed=5); z2 = nq.make_founder('male')
    c.check('routing/no-quant-species/readers-work', type(z2.get_trait_value('delet')) is float and z2.get_trait_value('ntrl') is None and z2.get_trait_alleles('delet', 'get_allele').shape == (30, 2),
            'species with traits_num 0 (ntrl + delet): both readers work on its traits')
    c.check('routing/no-quant-species/quant', True, "get_trait_value('quant') -> %s; get_trait_alleles('quant', 'get_allele_value') -> %s"
            % (raises(z2.get_trait_value, 'quant'), raises(z2.get_trait_alleles, 'quant', 'get_allele_value')), note=True)
    # recycled individual
    dead = inds[-1]; sp.recycle(dead)
    e = [raises(dead.get_trait_value, k, **({'trait_dim': 0} if k == 'quant' else {})) for k in KEYS] + [raises(dead.get_trait_alleles, k, m) for k, m in VALID_PAIRS]
    c.check('routing/recycled', all(x and x.startswith('RuntimeError') for x in e), '7 get_trait_value + 8 get_trait_alleles calls on a recycled individual -> %s' % e[0])

def g_kwargs(c):
    np = c.np
    sp = c.species(all_traits(0.5), seed=6)
    y = sp.make_founder('female')
    bad, calls = [], []
    def expect(tag, f, *a, **kw):
        e = raises(f, *a, **kw); calls.append(tag)
        if not (e and e.startswith('TypeError')): bad.append('%s -> %s' % (tag, e))
        return e
    # get_trait_value
    ex = expect('value(quant, trait_dims=1)', y.get_trait_value, 'quant', trait_dims=1)
    expect('value(quant, trait=1)', y.get_trait_value, 'quant', trait=1)
    expect('value(quant, 1)', y.get_trait_value, 'quant', 1)
    expect('value(quant, trait_dim=0, x=1)', y.get_trait_value, 'quant', trait_dim=0, x=1)
    for k in KEYS[1:]:
        expect('value(%s, trait_dim=1)' % k, y.get_trait_value, k, trait_dim=1)
        expect('value(%s, x=1)' % k, y.get_trait_value, k, x=1)
    # get_trait_alleles
    for k, m in VALID_PAIRS:
        if k == 'quant':
            expect('alleles(quant, %s, trait_dims=1)' % m, y.get_trait_alleles, k, m, trait_dims=1)
            expect('alleles(quant, %s, 1)' % m, y.get_trait_alleles, k, m, 1)
            expect('alleles(quant, %s, trait_dim=0, x=1)' % m, y.get_trait_alleles, k, m, trait_dim=0, x=1)
        else:
            expect('alleles(%s, %s, trait_dim=1)' % (k, m), y.get_trait_alleles, k, m, trait_dim=1)
            expect('alleles(%s, %s, x=1)' % (k, m), y.get_trait_alleles, k, m, x=1)
    c.check('kwargs/stray-or-misspelled -> TypeError', not bad, '%d calls (every branch of both readers: stray trait_dim, misspelled trait_dims, trait=, positional trait_dim, extra x=); e.g. %s'
            % (len(calls), ex) if not bad else '; '.join(bad[:4]))
    c.check('kwargs/invalid-pair-before-kwargs', (raises(y.get_trait_alleles, 'ntrl', 'get_allele_value', trait_dim=1) or '').startswith('ValueError'),
            'an invalid pair with kwargs -> ValueError (routing first)')
    c.check('kwargs/default-axis-0', y.get_trait_value('quant') == y.get_trait_value('quant', trait_dim=0) and
            all(np.array_equal(y.get_trait_alleles('quant', m), y.get_trait_alleles('quant', m, trait_dim=0)) for m in ('get_allele_value', 'get_allele')),
            'no trait_dim == trait_dim=0 in both readers')
    c.check('kwargs/axis-1-differs', y.get_trait_value('quant', trait_dim=1) != y.get_trait_value('quant', trait_dim=0) and
            not np.array_equal(y.get_trait_alleles('quant', 'get_allele_value', trait_dim=1), y.get_trait_alleles('quant', 'get_allele_value', trait_dim=0)),
            'trait_dim=1 reads the other axis')
    kw = {'trait_dim': 1}; y.get_trait_value('quant', **kw); y.get_trait_alleles('quant', 'get_allele_value', **kw)
    c.check('kwargs/caller-dict-unchanged', kw == {'trait_dim': 1}, "kwargs.pop('trait_dim') works on the function's own dict")

def g_quant_layout(c, case):
    np = c.np
    random.seed(2); np.random.seed(2)
    mnp = {'quanti_init_freq': 0.5} if case in ('diallelic', 'asymmetrical') else None
    dims = 3 if case == 'three_axes' else 2
    pheno = tuple('t%d' % i for i in range(dims))
    sp = c.species({'quant': quant_rows(case)}, pheno=pheno, seed=7, mnp=mnp)
    pop = founders(sp, 20)
    kids = [(sp.make_offspring(pop[1], pop[0], 'female' if j % 2 else 'male'), pop[1], pop[0]) for j in range(60)]
    inds = pop + [k for k, _, _ in kids]
    t = pop[0].get_trait('quant')
    loci = [t.get_num_loci(d) for d in range(dims)]
    raw_ok = all(np.array_equal(z.get_trait_alleles('quant', m, trait_dim=d), raw_table(np, z, 'quant', m, d)) for z in inds for d in range(dims) for m in ('get_allele_value', 'get_allele'))
    c.check('quant_layout/%s/alleles==raw' % case, raw_ok, '%d individuals x %d axes x 2 methods; loci per axis %s' % (len(inds), dims, loci))
    shape_ok = all(z.get_trait_alleles('quant', 'get_allele_value', trait_dim=d).shape == (loci[d], 2) for z in inds for d in range(dims))
    want_total = [40] * dims if case == 'pleiotropy' else ([30, 10] if case == 'per_trait' else ([14, 13, 13] if case == 'three_axes' else [20, 20]))
    c.check('quant_layout/%s/shape' % case, shape_ok and loci == want_total, 'shape (loci_d, 2), loci per axis %s (expected %s)' % (loci, want_total))
    val_ok = all(z.get_trait_value('quant', trait_dim=d) == z.get_indi_phenotype_ls()[d] == z.phenotype_set[pheno[d]] == stored(z, 'quant', d) for z in inds for d in range(dims))
    c.check('quant_layout/%s/value==legacy==stored' % case, val_ok, 'get_trait_value(quant, trait_dim=d) == get_indi_phenotype_ls()[d] == phenotype_set[t_d] == stored _phenotypes[d]')
    add_ok = all(abs(z.get_trait_alleles('quant', 'get_allele_value', trait_dim=d).sum() - z.get_trait('quant').get_genotypic_value(d)) < 1e-12 and
                 abs(z.get_trait('quant').get_genotypic_value(d) - z.get_trait_value('quant', trait_dim=d)) < 1e-12 for z in inds for d in range(dims))
    c.check('quant_layout/%s/sum==G==P' % case, add_ok, 'additive model, Ve 0: sum of the allele table == G == P on every axis')
    if case in ('diallelic', 'asymmetrical'):
        a = 0.01
        f = (lambda st: np.where(st == 0, a, -a)) if case == 'diallelic' else (lambda st: np.where(st == 0, 0.0, a))
        sv_ok = all(np.allclose(z.get_trait_alleles('quant', 'get_allele_value', trait_dim=d), f(z.get_trait_alleles('quant', 'get_allele', trait_dim=d)), atol=1e-15) for z in inds for d in range(dims))
        both = set(np.concatenate([z.get_trait_alleles('quant', 'get_allele', trait_dim=d).ravel() for z in pop for d in range(dims)]).tolist())
        c.check('quant_layout/%s/value-by-state' % case, sv_ok and both == {0, 1}, 'states {0, 1} both present; value = %s by state' % ('+a / -a' if case == 'diallelic' else '0 / a (quanti_allele_value {{0.01}})'))
        if case == 'asymmetrical':   # the scalar form of quanti_allele_value (same layout, so the same process is fine)
            sc = c.species({'quant': quant_rows('asymmetrical_scalar')}, seed=8, mnp=mnp).make_founder('female')
            pairs = sorted({(int(a), float(b)) for d in range(dims) for a, b in zip(sc.get_trait_alleles('quant', 'get_allele', trait_dim=d).ravel(), sc.get_trait_alleles('quant', 'get_allele_value', trait_dim=d).ravel())})
            c.check('quant_layout/asymmetrical/scalar-allele-value', True, 'diallelic_asymmetrical with scalar quanti_allele_value 0.01: Nemo stores state -> value %s (+a / -a, not 0 / a); the readers return exactly that (sum == G: %s)'
                    % (pairs, all(abs(sc.get_trait_alleles('quant', 'get_allele_value', trait_dim=d).sum() - sc.get_trait('quant').get_genotypic_value(d)) < 1e-12 for d in range(dims))), note=True)
        tab = lambda z: np.concatenate([z.get_trait_alleles('quant', 'get_allele', trait_dim=d) for d in range(dims)])
    else:
        c.check('quant_layout/%s/states-all-0' % case, all(not z.get_trait_alleles('quant', 'get_allele', trait_dim=d).any() for z in pop for d in range(dims)),
                'continuous allele model: get_allele reads 0 everywhere (documented in the helper)')
        tab = lambda z: np.concatenate([z.get_trait_alleles('quant', 'get_allele_value', trait_dim=d) for d in range(dims)])
    a_, b_ = orientation(np, kids, tab)
    c.check('quant_layout/%s/strands' % case, a_ == 0 and b_ > 0, 'mutation 0: column 1 maternal, column 0 paternal; mismatch as read %.3f, if swapped %.3f' % (a_, b_))
    if case in ('continuous', 'mapped'):    # recombination along the reader's locus order: maternal-origin switches between adjacent loci of axis 0
        sw = n = 0
        for k, mo, fa in kids:
            K, M = k.get_trait_alleles('quant', 'get_allele_value', trait_dim=0), mo.get_trait_alleles('quant', 'get_allele_value', trait_dim=0)
            orig = [0 if K[l, 1] == M[l, 0] else 1 for l in range(len(K))]
            sw += sum(orig[l] != orig[l + 1] for l in range(len(orig) - 1)); n += len(orig) - 1
        r = sw / n
        ok = (0.4 < r < 0.6) if case == 'continuous' else (r < 0.12)
        c.check('quant_layout/%s/recombination' % case, ok, 'adjacent-locus switch rate of the maternal strand %.3f (%s)' % (r, 'free recombination, expected ~0.5' if case == 'continuous' else 'quanti_recombination_rate 0.05, expected ~0.05'))

def g_quant_genetics(c):
    np = c.np
    random.seed(3); np.random.seed(3)
    k = 0.3
    dom = c.species({'quant': quant_rows('continuous', extra=[('quanti_dominance_effects', k)])}, seed=8)
    pop = founders(dom, 20); inds = pop + [dom.make_offspring(pop[1], pop[0], 'female') for _ in range(40)]
    def g_dom(z, d):
        a = z.get_trait_alleles('quant', 'get_allele_value', trait_dim=d)
        return float(sum(a0 + a1 + k * abs(a1 - a0) for a0, a1 in a))
    ok = all(abs(z.get_trait_value('quant', trait_dim=d) - z.get_trait('quant').get_genotypic_value(d)) < 1e-12 and abs(z.get_trait_value('quant', trait_dim=d) - g_dom(z, d)) < 1e-12 for z in inds for d in range(2))
    gap = min(abs(z.get_trait_value('quant', trait_dim=d) - z.get_trait_alleles('quant', 'get_allele_value', trait_dim=d).sum()) for z in inds for d in range(2))
    c.check('quant_genetics/dominance', ok and gap > 1e-6, 'quanti_dominance_effects %g, Ve 0: value == G == sum(a0 + a1 + k|a1 - a0|) on every axis; min |value - additive sum| = %.2e > 0' % (k, gap))
    ve = 0.5
    env = c.species({'quant': quant_rows('continuous', ve=ve)}, seed=9)
    big = founders(env, 4000)
    r = np.array([[z.get_trait_value('quant', trait_dim=d) - z.get_trait('quant').get_genotypic_value(d) for d in range(2)] for z in big])
    m, v = r.mean(axis=0), r.var(axis=0, ddof=1)
    se_m, se_v = math.sqrt(ve / len(big)), ve * math.sqrt(2.0 / (len(big) - 1))
    ok = all(abs(m) < 4 * se_m) and all(abs(v - ve) < 4 * se_v) and all(z.get_trait_value('quant', trait_dim=d) == stored(z, 'quant', d) for z in big[:200] for d in range(2))
    c.check('quant_genetics/Ve', ok, '4000 founders, Ve %g: value == stored P; P - G mean %s (4 se = %.3f), variance %s (4 se = %.3f)' % (ve, m.round(4).tolist(), 4 * se_m, v.round(4).tolist(), 4 * se_v))

def g_ntrl(c):
    np = c.np
    random.seed(4); np.random.seed(4)
    for n_all in (4, 2):                        # byte (ntrl_all > 2) and bitstring (ntrl_all 2) implementations
        sp = c.species({'quant': quant_rows('continuous'), 'ntrl': ntrl_rows(n_all)}, seed=10 + n_all)
        pop = founders(sp, 20)
        kids = [(sp.make_offspring(pop[1], pop[0], 'male'), pop[1], pop[0]) for _ in range(50)]
        inds = pop + [k for k, _, _ in kids]
        T = lambda z: z.get_trait_alleles('ntrl', 'get_allele')
        seen = set(np.concatenate([T(z) for z in pop]).ravel().tolist())
        c.check('ntrl/all%d/alleles==raw' % n_all, all(np.array_equal(T(z), raw_table(np, z, 'ntrl', 'get_allele')) and T(z).shape == (12, 2) for z in inds) and seen == set(range(n_all)),
                'int (12, 2) == raw for %d individuals; states seen %s' % (len(inds), sorted(seen)))
        a_, b_ = orientation(np, kids, T)
        c.check('ntrl/all%d/strands' % n_all, a_ == 0 and b_ > 0, 'mutation 0: mismatch as read %.3f, if swapped %.3f' % (a_, b_))
        c.check('ntrl/all%d/value-None' % n_all, all(z.get_trait_value('ntrl') is None and stored(z, 'ntrl') is None for z in inds), "get_trait_value('ntrl') is None (getValue() NULL)")

def g_delet(c):
    np = c.np
    random.seed(5); np.random.seed(5)
    s, h = 0.05, 0.3
    for model in (1, 2):
        for scale in (None, 0.9):
            sp = c.species({'quant': quant_rows('continuous'), 'delet': delet_rows(0.05, model=model, scale=scale, init=0.3, s=s, h=h)}, seed=20 + model)
            every, _ = generations(sp, founders(sp, 20), 3)
            res = [(z.get_trait_value('delet'), delet_formula(z.get_trait_alleles('delet', 'get_allele'), s, h, model, 1.0 if scale is None else scale), stored(z, 'delet')) for z in every]
            ok = all(abs(v - f[0]) < 1e-12 and v == st for v, f, st in res)
            hmz = [f[1] for _, f, _ in res]; htz = [f[2] for _, f, _ in res]
            c.check('delet/model%d/scale%s/value==formula==stored' % (model, scale), ok and len(set(v for v, _, _ in res)) > 10,
                    '%d individuals (founders delet_init_freq 0.3 + 3 generations, mutation 0.05): max |value - formula| %.1e; hmz %d..%d, htz %d..%d; W %.4f..%.4f'
                    % (len(res), max(abs(v - f[0]) for v, f, _ in res), min(hmz), max(hmz), min(htz), max(htz), min(v for v, _, _ in res), max(v for v, _, _ in res)))
    sp = c.species({'quant': quant_rows('continuous'), 'delet': delet_rows(0.2, init=0.3)}, seed=25)
    pop = founders(sp, 20)
    T = lambda z: z.get_trait_alleles('delet', 'get_allele')
    bad = bad_sw = 0
    for j in range(10):
        mo, fa = pop[2 * j + 1], pop[2 * j]; Mo, Fa = T(mo), T(fa)
        for _ in range(50):
            K = T(sp.make_offspring(mo, fa, 'female'))
            for l in range(30):
                bad += (K[l, 1] == 0 and 0 not in Mo[l]) + (K[l, 0] == 0 and 0 not in Fa[l])
                bad_sw += (K[l, 0] == 0 and 0 not in Mo[l]) + (K[l, 1] == 0 and 0 not in Fa[l])
    c.check('delet/strands', bad == 0 and bad_sw > 0, 'no back mutation (a 0 is only inherited): violations as read %d, if swapped %d (500 kids x 30 loci)' % (bad, bad_sw))
    c.check('delet/alleles==raw', all(np.array_equal(T(z), raw_table(np, z, 'delet', 'get_allele')) for z in pop), 'int (30, 2) == raw')

def g_dmi(c):
    np = c.np
    random.seed(6); np.random.seed(6)
    for tag, table in (('table16x3rows', DMI16), ('table9x1row', DMI9)):
        sp = c.species({'quant': quant_rows('continuous'), 'dmi': dmi_rows(0.2, table)}, seed=30)
        every, _ = generations(sp, founders(sp, 20), 6)
        t16 = dmi_table16(table)
        res = [(z.get_trait_value('dmi'),) + dmi_formula(z.get_trait_alleles('dmi', 'get_allele'), t16) + (stored(z, 'dmi'),) for z in every]
        ok = all(abs(v - f) < 1e-12 and v == st for v, f, _, st in res)
        codes = sorted({cd for _, _, cs, _ in res for cd in cs})
        cover = {p: sorted({cd for q, cd in codes if q == p}) for p in range(4)}
        c.check('dmi/%s/value==formula==stored' % tag, ok, '%d individuals (6 generations, mutation 0.2, 4 locus pairs): max |value - formula| %.1e; W %.4f..%.4f; genotype codes seen per pair %s'
                % (len(res), max(abs(v - f) for v, f, _, _ in res), min(r[0] for r in res), max(r[0] for r in res), {p: len(v) for p, v in cover.items()}))
        c.data('dmi/%s/codes' % tag, cover)
    sp = c.species({'quant': quant_rows('continuous'), 'dmi': dmi_rows(0.1)}, seed=31)
    every, last = generations(sp, founders(sp, 40), 6)
    T = lambda z: z.get_trait_alleles('dmi', 'get_allele')
    hit = hit_sw = n = 0
    for mo in [z for z in every if z.gender == 'female']:
        for fa in [z for z in every if z.gender == 'male'][:5]:
            Mo, Fa = T(mo), T(fa)
            loci = [l for l in range(8) if Mo[l, 0] == Mo[l, 1] and Fa[l, 0] == Fa[l, 1] and Mo[l, 0] != Fa[l, 0]]
            if not loci: continue
            K = T(sp.make_offspring(mo, fa, 'female'))
            for l in loci:
                n += 1; hit += (K[l, 1] == Mo[l, 0] and K[l, 0] == Fa[l, 0]); hit_sw += (K[l, 0] == Mo[l, 0] and K[l, 1] == Fa[l, 0])
            if n > 2000: break
        if n > 2000: break
    c.check('dmi/strands', n > 200 and hit / n > 0.75 and hit_sw / n < 0.1, 'opposite-homozygote parents (%d loci, mutation 0.1): column 1 maternal and column 0 paternal %.3f, swapped %.3f' % (n, hit / n, hit_sw / n))

def g_disp(c):
    np = c.np
    random.seed(7); np.random.seed(7)
    sp = c.species({'quant': quant_rows('continuous'), 'disp': disp_rows(0.0, 0.1, 0.37, 0.81)}, seed=40)
    pop = founders(sp, 6)
    c.check('disp/founder-init', all(z.get_trait_value('fdisp') == 0.37 and z.get_trait_value('mdisp') == 0.81 and np.array_equal(z.get_trait_alleles('fdisp', 'get_allele_value'), [[0.37, 0.37]]) for z in pop),
            'disp_init_rate_fem 0.37 / _mal 0.81: fdisp value 0.37 and alleles (0.37, 0.37), mdisp 0.81, both genders')
    mo, fa = pop[1], pop[0]
    t = mo.get_trait('fdisp'); before = mo.get_trait_value('fdisp')
    t.set_allele_value(0, 0, 0.11); t.set_allele_value(0, 1, 0.12)
    stale = mo.get_trait_value('fdisp'); t.set_value(); fresh = mo.get_trait_value('fdisp')
    c.check('disp/set_value-writes-get_value-reads', stale == before and abs(fresh - 0.115) < 1e-15 and np.array_equal(mo.get_trait_alleles('fdisp', 'get_allele_value'), [[0.11, 0.12]]),
            'after set_allele_value the alleles reader sees (0.11, 0.12) at once, get_trait_value stays %.4f until set_value(), then %.4f' % (stale, fresh))
    for z, vals in ((mo, {'mdisp': (0.21, 0.22)}), (fa, {'fdisp': (0.77, 0.78), 'mdisp': (0.67, 0.68)})):
        for k, (a0, a1) in vals.items():
            tt = z.get_trait(k); tt.set_allele_value(0, 0, a0); tt.set_allele_value(0, 1, a1); tt.set_value()
    kids = [(sp.make_offspring(mo, fa, 'female' if j % 2 else 'male'), mo, fa) for j in range(60)]
    for k in ('fdisp', 'mdisp'):
        T = lambda z, k=k: z.get_trait_alleles(k, 'get_allele_value')
        a_, b_ = orientation(np, kids, T)          # orientation() asks 'column 1 maternal'; disp is reversed, so expect the swapped reading to be exact
        c.check('disp/%s/strands' % k, b_ == 0 and a_ > 0, 'mutation 0: column 0 maternal, column 1 paternal (reversed vs the other traits); mismatch if column 1 were maternal %.3f, as documented %.3f' % (a_, b_))
    sp2 = c.species({'quant': quant_rows('continuous'), 'disp': disp_rows(0.5, 0.3, 0.5, 0.5)}, seed=41)
    every, _ = generations(sp2, founders(sp2, 20), 4)
    for k in ('fdisp', 'mdisp'):
        vals = [z.get_trait_value(k) for z in every]
        ok = all(abs(v - z.get_trait_alleles(k, 'get_allele_value').mean()) < 1e-15 and v == stored(z, k) for v, z in zip(vals, every))
        al = np.concatenate([z.get_trait_alleles(k, 'get_allele_value').ravel() for z in every])
        c.check('disp/%s/mutation' % k, ok and al.min() >= 0.0 and al.max() <= 1.0 and len(set(vals)) > 20,
                '%d individuals, mutation 0.5 mean 0.3: value == mean of the alleles == stored; alleles in [%.3f, %.3f] (fix 640b9da)' % (len(every), al.min(), al.max()))

def g_wolb(c):
    np = c.np
    random.seed(8); np.random.seed(8)
    sp = c.species({'quant': quant_rows('continuous'), 'wolb': wolb_rows(1.0)}, seed=50)
    pop = founders(sp, 20)
    for j, z in enumerate(pop):
        if z.gender == 'female' and j % 4 == 1: z.get_trait('wolb').set_allele_value(0, 0, 1.0)
        if z.gender == 'male' and j % 4 == 0: z.get_trait('wolb').set_allele_value(0, 0, 1.0)
    F = [z for z in pop if z.gender == 'female']; M = [z for z in pop if z.gender == 'male']
    trios = []
    for _ in range(300):
        mo, fa = random.choice(F), random.choice(M); trios.append((sp.make_offspring(mo, fa, 'female' if random.random() < 0.5 else 'male'), mo, fa))
    inds = pop + [k for k, _, _ in trios]
    c.check('wolb/value==alleles==stored', all(type(z.get_trait_value('wolb')) is bool and z.get_trait_value('wolb') == z.get_trait_alleles('wolb', 'get_allele') == stored(z, 'wolb') for z in inds),
            '%d individuals, %d infected' % (len(inds), sum(z.get_trait_value('wolb') for z in inds)))
    mat = all(k.get_trait_value('wolb') == mo.get_trait_value('wolb') for k, mo, fa in trios)
    pat = [(k, mo, fa) for k, mo, fa in trios if fa.get_trait_value('wolb') and not mo.get_trait_value('wolb')]
    c.check('wolb/maternal-rate1', mat and pat and not any(k.get_trait_value('wolb') for k, _, _ in pat),
            'transmission 1: kid infected iff its mother is (300 kids); %d kids of infected father x uninfected mother: none infected' % len(pat))
    cl = [(sp.make_clone(mo), mo) for mo in F for _ in range(3)]
    c.check('wolb/clone', all(x.get_trait_value('wolb') == mo.get_trait_value('wolb') for x, mo in cl), 'make_clone keeps the mother state (%d clones)' % len(cl))
    z = [y for y in pop if not y.get_trait_value('wolb')][0]
    z.get_trait('wolb').set_allele_value(0, 0, 1.0); on = z.get_trait_value('wolb'); z.get_trait('wolb').set_allele_value(0, 0, 0.0); off = z.get_trait_value('wolb')
    c.check('wolb/set_allele_value', on is True and off is False, 'set_allele_value(0, 0, 1) -> True, (0, 0, 0) -> False (no set_value() needed: the flag is the value)')
    sp2 = c.species({'quant': quant_rows('continuous'), 'wolb': wolb_rows(0.5)}, seed=51)
    mo, fa = sp2.make_founder('female'), sp2.make_founder('male')
    mo.get_trait('wolb').set_allele_value(0, 0, 1.0)
    kids = [sp2.make_offspring(mo, fa, 'female') for _ in range(4000)]
    fr = sum(k.get_trait_value('wolb') for k in kids) / len(kids); se = math.sqrt(0.25 / len(kids))
    c.check('wolb/rate0.5', abs(fr - 0.5) < 4 * se, 'transmission 0.5: infected fraction among 4000 kids of an infected mother %.4f (4 se = %.4f)' % (fr, 4 * se))

def g_stored(c):
    np = c.np
    sp, every = dump_species(c)
    bad = []
    for z in every:
        for d in range(2):
            if z.get_trait_value('quant', trait_dim=d) != stored(z, 'quant', d): bad.append('quant%d' % d)
        for k in DBL + ('wolb',):
            if z.get_trait_value(k) != stored(z, k): bad.append(k)
        if z.get_trait_value('ntrl') is not None: bad.append('ntrl')
        for k, m in VALID_PAIRS:
            if k == 'quant':
                if not all(np.array_equal(z.get_trait_alleles(k, m, trait_dim=d), raw_table(np, z, k, m, d)) for d in range(2)): bad.append('alleles quant ' + m)
            elif k == 'wolb':
                if z.get_trait_alleles(k, m) != raw_table(np, z, k, m): bad.append('alleles wolb')
            elif not np.array_equal(z.get_trait_alleles(k, m), raw_table(np, z, k, m)): bad.append('alleles ' + k)
    c.check('stored/value==stored', not [b for b in bad if not b.startswith('alleles')], '%d individuals (7 traits, every mutation on, Ve 0.5): get_trait_value == stored for quant x2, delet, dmi, fdisp, mdisp, wolb; ntrl None' % len(every))
    c.check('stored/alleles==raw', not [b for b in bad if b.startswith('alleles')], '%d individuals x 8 valid pairs (quant x2 axes)' % len(every))
    c.check('stored/infected', 0 < sum(z.get_trait_value('wolb') for z in every) < len(every), '%d of %d infected' % (sum(z.get_trait_value('wolb') for z in every), len(every)))

def snapshot(np, z):
    ''' every value and allele table of one individual, JSON-friendly '''
    out = {'P': [z.get_trait_value('quant', trait_dim=d) for d in range(z.traits_num)]}
    for k in DBL + ('wolb', 'ntrl'): out[k] = z.get_trait_value(k)
    for k, m in VALID_PAIRS:
        if k == 'quant':
            for d in range(z.traits_num): out['%s/%s/%d' % (k, m, d)] = z.get_trait_alleles(k, m, trait_dim=d).tolist()
        else:
            a = z.get_trait_alleles(k, m); out['%s/%s' % (k, m)] = a if isinstance(a, bool) else a.tolist()
    return out

def g_lifecycle(c):
    np = c.np
    random.seed(9); np.random.seed(9)
    A = c.species(all_traits(1.0), seed=60)
    B = c.species(all_traits(0.0), seed=61)
    everyA, popA = generations(A, founders(A, 20), 2)
    popB = founders(B, 20)
    for z in popA[:10] + popB[:10]:
        if z.gender == 'female': z.get_trait('wolb').set_allele_value(0, 0, 1.0)
    cps = [z.make_copy() for z in popA]
    c.check('lifecycle/make_copy', all(snapshot(np, z) == snapshot(np, x) for z, x in zip(popA, cps)), '%d copies: every value and allele table equal (delet / dmi recomputed by operator=, disp copied)' % len(cps))
    for x in cps: x.recycle()
    clB = [(B.make_clone(z), z) for z in popB]
    c.check('lifecycle/make_clone-mutation0', all(snapshot(np, x) == snapshot(np, z) for x, z in clB), '%d clones without mutation == their mother' % len(clB))
    clA = [(A.make_clone(z), z) for z in popA]
    fresh = [(x, x.make_copy()) for x, _ in clA]
    c.check('lifecycle/make_clone-mutation', all(snapshot(np, x) == snapshot(np, y) for x, y in fresh),
            '%d clones with mutation: stored values current (== a fresh make_copy, which re-runs set_value); delet W differs from the mother in %d'
            % (len(clA), sum(x.get_trait_value('delet') != z.get_trait_value('delet') for x, z in clA)))
    for _, y in fresh: y.recycle()
    # recycle -> reuse from the pool
    old = popA[:20]
    old_addr = {int(z.nemo_indi_obj.this): z.get_trait_value('quant') for z in old}
    A.purge_recycling_pool()                    # the pool is FIFO: empty it first, so the next 20 births reuse exactly these 20 objects
    for z in old: z.recycle()
    rest = popA[20:]
    F = [z for z in rest if z.gender == 'female']; M = [z for z in rest if z.gender == 'male']
    born = [A.make_offspring(random.choice(F), random.choice(M), 'female') for _ in range(20)]
    reused = [z for z in born if int(z.nemo_indi_obj.this) in old_addr]
    fresh_ok = all(z.get_trait_value(k) == stored(z, k) for z in born for k in DBL + ('wolb',)) and all(z.get_trait_value('quant') == z.get_indi_phenotype_ls()[0] for z in born) and \
               all(abs(z.get_trait_value('fdisp') - z.get_trait_alleles('fdisp', 'get_allele_value').mean()) < 1e-15 for z in born) and \
               all(np.array_equal(z.get_trait_alleles(k, 'get_allele'), raw_table(np, z, k, 'get_allele')) for z in born for k in ('ntrl', 'delet', 'dmi'))
    c.check('lifecycle/recycled-memory-reused', len(reused) == 20 and fresh_ok and all(z.get_trait_value('quant') != old_addr[int(z.nemo_indi_obj.this)] for z in reused),
            '%d of 20 newborns built on recycled Nemo objects: values == stored / legacy / disp mean, alleles == raw, none equals the previous occupant' % len(reused))
    e = [raises(z.get_trait_value, 'delet') for z in old] + [raises(z.get_trait_alleles, 'delet', 'get_allele') for z in old]
    c.check('lifecycle/old-wrappers', all(x and x.startswith('RuntimeError') for x in e), '20 recycled wrappers: both readers -> RuntimeError (their Nemo objects now hold newborns)')
    alive = [z for z in everyA + born + [x for x, _ in clA] + popB + [x for x, _ in clB] if z.is_alive]
    for z in alive: z.indfactory.recycle(z)
    c.check('lifecycle/minted==recycled', all(f.minted_num == f.recycled_num for f in (A, B)), 'after recycling every survivor: A %s, B %s (minted, recycled)' % ((A.minted_num, A.recycled_num), (B.minted_num, B.recycled_num)))

def g_out_of_range(c, case):
    ''' NOTE only: an out-of-range trait_dim (no guard by design; Nemo's asserts are compiled out) -- one probe per process '''
    sp = c.species({'quant': quant_rows('continuous')}, seed=70)
    z = sp.make_founder('female')
    call = {'value_dim2': lambda: z.get_trait_value('quant', trait_dim=2), 'value_dim-1': lambda: z.get_trait_value('quant', trait_dim=-1),
            'alleles_values_dim2': lambda: z.get_trait_alleles('quant', 'get_allele_value', trait_dim=2), 'alleles_states_dim2': lambda: z.get_trait_alleles('quant', 'get_allele', trait_dim=2),
            'alleles_dim-1': lambda: z.get_trait_alleles('quant', 'get_allele_value', trait_dim=-1)}[case]
    print('@@ ' + json.dumps({'group': c.group, 'id': 'out_of_range/%s/started' % case, 'status': 'NOTE', 'detail': 'started'}), flush=True)
    np = c.np
    try:
        r = call()
        if hasattr(r, 'shape'):
            t = z.get_trait('quant'); m = 'get_allele' if 'states' in case else 'get_allele_value'
            same = [d for d in range(2) if np.array_equal(r, z.get_trait_alleles('quant', m, trait_dim=d))]
            pos = [t.get_locus_seq_pos(l, 2) for l in range(t.get_num_loci(2))]
            desc = 'returned array shape %s (get_num_loci(2) = %d, positions %s..%s); equal to axis %s; first rows %s' % (r.shape, t.get_num_loci(2), min(pos), max(pos), same or 'none', r[:2].tolist())
        else:
            desc = 'returned %r (axis 0 = %r, axis 1 = %r)' % (r, z.get_trait_value('quant', trait_dim=0), z.get_trait_value('quant', trait_dim=1))
    except Exception as e: desc = 'raised %s: %s' % (type(e).__name__, e)
    c.check('out_of_range/%s' % case, True, desc, note=True)

def g_dump(c):
    ''' allele tables + legacy readers of the 7-trait species: works on HEAD individual.py and on pynemo without get_value; RETURN via @@DATA '''
    np = c.np
    sp, every = dump_species(c)
    rows = []
    for z in every:
        r = {'id': z.individual_id, 'gender': z.gender, 'P': [repr(x) for x in z.get_indi_phenotype_ls()], 'pset': [repr(v) for v in z.phenotype_set.values()]}
        for k, m in VALID_PAIRS:
            if k == 'quant':
                for d in range(2): r['%s/%s/%d' % (k, m, d)] = z.get_trait_alleles(k, m, trait_dim=d).tolist()
            else:
                a = z.get_trait_alleles(k, m); r['%s/%s' % (k, m)] = a if isinstance(a, bool) else a.tolist()
        rows.append(r)
    c.data('dump/md5', hashlib.md5(json.dumps(rows).encode()).hexdigest()); c.data('dump/n', len(rows))
    c.data('dump/has_get_trait_value', hasattr(every[0], 'get_trait_value'))
    c.data('dump/value_delet', raises(every[0].get_trait_value, 'delet') if hasattr(every[0], 'get_trait_value') else 'no get_trait_value')
    c.check('dump/ran', True, '%d individuals' % len(rows), note=True)

def g_perf(c):
    np = c.np
    sp = c.species(all_traits(0.5), seed=80)
    z = sp.make_founder('female')
    out = []
    for label, f in [('value quant', lambda: z.get_trait_value('quant', trait_dim=1)), ('value delet', lambda: z.get_trait_value('delet')), ('value wolb', lambda: z.get_trait_value('wolb')),
                     ('value ntrl', lambda: z.get_trait_value('ntrl')), ('legacy get_indi_phenotype_ls', lambda: z.get_indi_phenotype_ls()),
                     ('alleles quant values (20 loci)', lambda: z.get_trait_alleles('quant', 'get_allele_value', trait_dim=1)), ('alleles delet (30 loci)', lambda: z.get_trait_alleles('delet', 'get_allele')),
                     ('alleles fdisp', lambda: z.get_trait_alleles('fdisp', 'get_allele_value')), ('alleles wolb', lambda: z.get_trait_alleles('wolb', 'get_allele'))]:
        n = 20000 if label.startswith(('value', 'legacy', 'alleles wolb', 'alleles fdisp')) else 2000
        t0 = time.perf_counter()
        for _ in range(n): f()
        out.append('%s %.2f us' % (label, (time.perf_counter() - t0) / n * 1e6))
    c.check('perf/per-call', True, '; '.join(out), note=True)

##****************************** parent side ******************************
def run_child(group, root=PROJECT_ROOT, py=sys.executable, timeout=600):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    t0 = time.time()
    try:
        r = subprocess.run([str(py), '-I', '-B', str(Path(__file__).resolve()), '--child', group, str(root)], capture_output=True, text=True, env=env, cwd=tempfile.gettempdir(), timeout=timeout)
        raw, code = r.stdout + r.stderr, r.returncode
    except subprocess.TimeoutExpired as e:
        raw, code = ((e.stdout or b'').decode() if isinstance(e.stdout, bytes) else (e.stdout or '')), 'timeout %ds' % timeout
    checks = [json.loads(l[3:]) for l in raw.splitlines() if l.startswith('@@ ')]
    data = {d['key']: d['value'] for d in (json.loads(l[7:]) for l in raw.splitlines() if l.startswith('@@DATA '))}
    return checks, data, raw, code, time.time() - t0

def head_copy():
    ''' metaibm/ + extension/ of the working tree with HEAD's individual.py, in a temp dir '''
    d = Path(tempfile.mkdtemp(prefix='ntr_head_'))
    for pkg in ('metaibm', 'extension'): shutil.copytree(PROJECT_ROOT / pkg, d / pkg, ignore=shutil.ignore_patterns('__pycache__'))
    head = subprocess.run(['git', 'show', 'HEAD:metaibm/individual.py'], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True).stdout
    (d / 'metaibm' / 'individual.py').write_text(head)
    return d

def md5(p): return hashlib.md5(Path(p).read_bytes()).hexdigest()

def main(groups):
    import pynemo
    so = [f for f in os.listdir(os.path.dirname(pynemo.__file__)) if f.startswith('_nemo')][0]
    print('python %s | pynemo %s (%s md5 %s) | individual.py md5 %s' % (sys.executable, os.path.dirname(pynemo.__file__), so, md5(Path(pynemo.__file__).parent / so), md5(PROJECT_ROOT / 'metaibm' / 'individual.py')))
    results, t_all = [], time.time()
    def add(group, checks, code, raw):
        results.extend(checks)
        if code != 0: results.append({'group': group, 'id': '%s/child-exit' % group, 'status': 'FAIL', 'detail': 'exit %s: %s' % (code, raw[-600:].replace('\n', ' '))})
    for g in groups:
        if g == 'out_of_range':
            for case in OOR_CASES:
                ch, _, raw, code, dt = run_child('out_of_range:' + case, timeout=60)
                ch = [x for x in ch if not x['id'].endswith('/started')]
                if code != 0:        # a crash or a hang is the recorded behaviour here, not a test failure
                    ch.append({'group': 'out_of_range', 'id': 'out_of_range/%s' % case, 'status': 'NOTE', 'detail': 'child ended with %s (%s)' % (code, 'SIGSEGV' if code in (-11, 139) else raw.strip().splitlines()[-1][:200] if raw.strip() else 'no output')})
                results.extend(ch)
                print('  %-24s %5.1fs  %s' % ('out_of_range:' + case, dt, ' '.join(x['status'][0] for x in ch)), flush=True)
            continue
        if g == 'regression_head':
            hd = head_copy()
            try:
                ch, dh, rh, kh, t1 = run_child('dump', root=hd); cw, dw, rw, kw, t2 = run_child('dump')
                add(g, [], kh, rh); add(g, [], kw, rw)
                ok = dh.get('dump/md5') and dh.get('dump/md5') == dw.get('dump/md5')
                results.append({'group': g, 'id': 'regression_head/alleles+legacy', 'status': 'PASS' if ok else 'FAIL',
                                'detail': 'HEAD individual.py (has get_trait_value: %s) vs working tree (%s), same seeds, %s individuals: md5 %s / %s'
                                % (dh.get('dump/has_get_trait_value'), dw.get('dump/has_get_trait_value'), dw.get('dump/n'), str(dh.get('dump/md5'))[:12], str(dw.get('dump/md5'))[:12])})
            finally: shutil.rmtree(hd, ignore_errors=True)
            print('  %-24s %5.1fs  %s' % (g, t1 + t2, results[-1]['status'][0]), flush=True)
            continue
        if g == 'regression_env':
            if not PY_ALT.exists():
                results.append({'group': g, 'id': 'regression_env/skipped', 'status': 'NOTE', 'detail': '%s not found' % PY_ALT}); continue
            c5, d5, r5, k5, t1 = run_child('dump', py=PY_ALT); c6, d6, r6, k6, t2 = run_child('dump')
            add(g, [], k5, r5); add(g, [], k6, r6)
            ok = d5.get('dump/md5') and d5.get('dump/md5') == d6.get('dump/md5')
            results.append({'group': g, 'id': 'regression_env/alleles+legacy test5==test6', 'status': 'PASS' if ok else 'FAIL',
                            'detail': 'pynemo of metaibm-nemo-test5 vs this env, same seeds, %s individuals: md5 %s / %s' % (d6.get('dump/n'), str(d5.get('dump/md5'))[:12], str(d6.get('dump/md5'))[:12])})
            results.append({'group': g, 'id': 'regression_env/dependency', 'status': 'NOTE',
                            'detail': "get_trait_value('delet') on test5: %s; on this env: %s" % (d5.get('dump/value_delet') or 'no error', d6.get('dump/value_delet') or 'no error')})
            print('  %-24s %5.1fs  %s' % (g, t1 + t2, results[-2]['status'][0]), flush=True)
            continue
        ch, _, raw, code, dt = run_child(g)
        add(g, ch, code, raw)
        print('  %-24s %5.1fs  %s' % (g, dt, ' '.join(x['status'][0] for x in ch)), flush=True)
    n = {s: sum(1 for x in results if x['status'] == s) for s in ('PASS', 'FAIL', 'NOTE')}
    print('\n| status | id | detail |\n|---|---|---|')
    for x in results: print('| %s | `%s` | %s |' % (x['status'], x['id'], x['detail'].replace('|', '\\|')))
    print('\nPASS %d  FAIL %d  NOTE %d  (%.0fs)' % (n['PASS'], n['FAIL'], n['NOTE'], time.time() - t_all))
    return 1 if n['FAIL'] else 0

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--child':
        group = sys.argv[2]
        ctx = Ctx(group, Path(sys.argv[3]))
        try:
            name, _, arg = group.partition(':')
            f = globals()['g_' + name]
            f(ctx, arg) if arg else f(ctx)
        finally: ctx.close()
    else:
        sys.exit(main([g for g in sys.argv[1:]] or list(GROUPS)))
