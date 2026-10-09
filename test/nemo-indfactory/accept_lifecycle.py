# many generations of birth + recycle -- IDs, counters, recycling-pool reuse, speed
import random, time
from common import nemo_params, indfactory

random.seed(1)
sp = indfactory('sp_life', ['pheno1', 'pheno2'], nemo_params('sp_life', 0.2), seed=11)
N, GENS = 200, 30
pop = [sp.make_founder('female' if i % 2 == 0 else 'male') for i in range(N)]
seen_ids = {indi.individual_id for indi in pop}
recycled_ptrs, reused, newborn_total = set(), 0, 0
t0 = time.time()
for gen in range(1, GENS + 1):
    females = [indi for indi in pop if indi.gender == 'female']
    males = [indi for indi in pop if indi.gender == 'male']
    offspring = [sp.make_offspring(random.choice(females), random.choice(males), random.choice(('female', 'male'))) for _ in range(N)]
    for child in offspring:
        reused += int(child.nemo_indi_obj.this) in recycled_ptrs
        newborn_total += 1
        assert child.individual_id not in seen_ids, 'ID reused: %d' % child.individual_id
        seen_ids.add(child.individual_id)
    for parent in pop:
        recycled_ptrs.add(int(parent.nemo_indi_obj.this))
        parent.recycle()
    pop = offspring
    assert sp.minted_num - sp.recycled_num == N, (gen, sp.minted_num, sp.recycled_num)
elapsed = time.time() - t0
print('%d gens x %d births: %.2f s (%.1f us per birth incl. recycle)' % (GENS, N, elapsed, 1e6 * elapsed / (GENS * N)))
print('minted %d, recycled %d' % (sp.minted_num, sp.recycled_num))
print('newborns placed in recycled memory: %d / %d' % (reused, newborn_total))
assert reused >= newborn_total - N, 'recycling pool not reused'
sp.purge_recycling_pool()
print('ACCEPT LIFECYCLE PASSED')
