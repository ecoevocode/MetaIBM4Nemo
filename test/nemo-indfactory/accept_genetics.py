# genotype_table / phenotype / inheritance read back correctly through individual.py
from common import nemo_params, indfactory

sp = indfactory('sp_gen', ['pheno1', 'pheno2'], nemo_params('sp_gen', 0.3, mutation_rate=0.0, env_variance=0.0), seed=7)
print(sp)

# 1. phenotype == sum of both allele values over the loci of that axis (additive, no mutation, Ve = 0)
founders = [sp.make_founder('female' if i % 2 == 0 else 'male') for i in range(10)]
max_err = 0.0
for indi in founders:
    for axis in range(sp.traits_num):
        table = indi.genotype_table(axis)
        assert len(table) == 20, len(table)
        max_err = max(max_err, abs(indi.get_indi_phenotype_ls()[axis] - sum(a + b for a, b in table)))
print('max |phenotype - sum(alleles)| = %.3e' % max_err)
assert max_err < 1e-9

# 2. phenotype_set mirrors get_indi_phenotype_ls and is read-only
ps = founders[0].phenotype_set
assert list(ps.values()) == founders[0].get_indi_phenotype_ls() and list(ps) == ['pheno1', 'pheno2']
try: ps['pheno1'] = 1.0; raise AssertionError('phenotype_set writable')
except TypeError: pass

# 3. inheritance: every offspring locus carries one allele from each parent (mutation rate 0)
mother, father = founders[0], founders[1]
votes = {0: 0, 1: 0}
for _ in range(30):
    child = sp.make_offspring(mother, father, 'female')
    for axis in range(sp.traits_num):
        for (c0, c1), (m0, m1), (f0, f1) in zip(child.genotype_table(axis), mother.genotype_table(axis), father.genotype_table(axis)):
            if c1 in (m0, m1) and c0 in (f0, f1): votes[1] += 1
            elif c0 in (m0, m1) and c1 in (f0, f1): votes[0] += 1
            else: raise AssertionError('allele not inherited from parents')
    child.recycle()
print('loci with maternal allele on strand 1 / strand 0:', votes[1], '/', votes[0])
print('ACCEPT GENETICS PASSED')
