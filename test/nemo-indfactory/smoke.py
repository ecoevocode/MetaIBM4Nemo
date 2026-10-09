import hashlib
from common import nemo_params, indfactory

# 1. boot two species in one process, same seed
sp1 = indfactory('sp1', ['pheno1', 'pheno2'], nemo_params('sp1', 0.2), seed=42)
sp2 = indfactory('sp2', ['pheno1', 'pheno2'], nemo_params('sp2', 0.6), seed=42)
print(sp1); print(sp2); print(repr(sp1))

# 2. founders: phenotypes centred on each species' quanti_init_value, parent IDs are sentinels
founders = {sp: [sp.make_founder('female' if i % 2 == 0 else 'male') for i in range(20)] for sp in (sp1, sp2)}
for sp, target in ((sp1, 0.2), (sp2, 0.6)):
    mean = sum(indi.get_indi_phenotype_ls()[0] for indi in founders[sp]) / 20
    print('%s founder mean phenotype[0] = %.4f (target %.1f)' % (sp.species_id, mean, target))
    assert abs(mean - target) < 0.1, (sp.species_id, mean)
    assert all(indi.mother_id >= 2**63 for indi in founders[sp])
print(founders[sp1][0])

# 3. offspring: pedigree copied from the parents
mother, father = founders[sp1][0], founders[sp1][1]
child = sp1.make_offspring(mother, father, 'female', home_habitat=3)
assert (child.mother_id, child.father_id) == (mother.individual_id, father.individual_id)
assert child.nemo_indi_obj.getHome() == 3
print(child)

# 4. IDs unique across both species (the ID counter survived both boots)
alive = founders[sp1] + founders[sp2] + [child]
ids = [indi.individual_id for indi in alive]
assert len(set(ids)) == len(ids), 'ID collision'
print('ids: min %d, max %d, n %d' % (min(ids), max(ids), len(ids)))

# 5. guards: wrong species, recycled individual, double recycle
try: sp1.make_offspring(founders[sp2][0], father, 'male'); raise AssertionError('cross-species parent accepted')
except RuntimeError as e: print('ok, refused:', e)
child.recycle()
print(child); print(repr(child))
for action in (lambda: child.gender, lambda: child.recycle(), lambda: sp1.make_offspring(child, father, 'male')):
    try: action(); raise AssertionError('recycled individual still usable')
    except RuntimeError: pass

# 6. recycling-pool reuse + counters
reborn = sp1.make_offspring(mother, father, 'male')
assert reborn.individual_id not in ids, 'reused memory got an old ID'
assert (sp1.minted_num, sp1.recycled_num) == (22, 1), (sp1.minted_num, sp1.recycled_num)
sp1.purge_recycling_pool()

digest = hashlib.sha1(repr([indi.get_indi_phenotype_ls() for indi in founders[sp1] + founders[sp2] + [reborn]]).encode()).hexdigest()
print('phenotype digest:', digest)
print('SMOKE TEST PASSED')
