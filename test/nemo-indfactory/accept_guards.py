# input guards, and what a FAILED boot does to the global ID counter
from common import nemo_params, indfactory

sp1 = indfactory('sp1', ['pheno1', 'pheno2'], nemo_params('sp1', 0.2), seed=3)
before = [sp1.make_founder('female') for _ in range(5)]
print('ids before failed boots:', [i.individual_id for i in before], ' counter:', indfactory.get_current_individual_id())

def expect(exc_type, fn, label):
    try: fn(); print('!! %s: no error' % label)
    except exc_type as e: print('ok %-36s %s: %s' % (label, exc_type.__name__, e))

expect(ValueError, lambda: sp1.make_founder('f'), 'bad gender')
expect(ValueError, lambda: indfactory('bad_a', ['pheno1'], nemo_params('bad_a', 0.2), seed=3), 'pheno names != quanti_traits')
print('counter after a failure BEFORE _init:', indfactory.get_current_individual_id())
expect(ValueError, lambda: indfactory('bad_c', ['pheno1', 'pheno2'], nemo_params('bad_c', 0.2), trait_keys=('quant', 'ntrl'), seed=3), 'trait_keys asks for unbuilt ntrl')
print('counter after a failure AFTER _init:', indfactory.get_current_individual_id())
for big in (70000, -1):
    try: sp1.make_founder('male', home_habitat=big); print('!! home_habitat=%d accepted' % big)
    except Exception as e: print('ok home_habitat=%d -> %s' % (big, type(e).__name__))

after = [sp1.make_founder('male') for _ in range(5)]
print('ids after failed boots:', [i.individual_id for i in after])
print('ID collision with earlier individuals:', sorted({i.individual_id for i in before} & {i.individual_id for i in after}) or 'none')
