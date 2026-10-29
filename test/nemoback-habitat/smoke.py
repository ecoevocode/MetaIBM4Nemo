# Smoke (< 30 s): imports + one call of every NemoBackend habitat function, invariants after each. run_all.py stops if this fails.
import sys, time, traceback
from common import seed_all, make_species, make_hab, check_invariants, recycle_everything, save

T0 = time.time()
STEPS, CURRENT = [], ['import']


def main():
    seed_all(0)
    sp = make_species('sp1', 0.2, seed=1)
    hab = make_hab(L=6, K=10, global_index=3)

    def step(name, fn):
        CURRENT[0] = name
        r = fn()
        STEPS.append({'step': name, 'return': repr(r), 'invariants': check_invariants(hab, sp, name)})
        print('[smoke] %-58s -> %r' % (name, r))
        return r

    step('hab_initialize', lambda: hab.hab_initialize({'sp1': sp}, [[0.2, 0.2]]))
    step('hab_dead_selection', lambda: hab.hab_dead_selection(0.1, 0.5, 'niche_gaussian'))
    step('hab_sex_reproduce_mutate_into_offspring_pool', lambda: hab.hab_sex_reproduce_mutate_into_offspring_pool(1.0))
    assert hab.offspring_pool[0].nemo_indi_obj.getHome() == 3
    step('hab_local_germinate_from_offspring_and_immigrant_pool', lambda: hab.hab_local_germinate_from_offspring_and_immigrant_pool())
    step('hab_clear_up_offspring_pool + immigrant_pool', lambda: hab.hab_clear_up_offspring_pool() + hab.hab_clear_up_immigrant_pool())
    hab.hab_dead_selection(0.1, 0.5, 'niche_gaussian')
    step('reproduce', lambda: hab.hab_sex_reproduce_mutate_into_offspring_pool(1.0))
    step('hab_local_germinate_from_offspring_immigrant_and_dormancy_pool', lambda: hab.hab_local_germinate_from_offspring_immigrant_and_dormancy_pool())
    step('hab_dormancy_process_from_offspring_pool_to_dormancy_pool', hab.hab_dormancy_process_from_offspring_pool_to_dormancy_pool)
    step('reproduce', lambda: hab.hab_sex_reproduce_mutate_into_offspring_pool(1.0))
    step('hab_dormancy_process_from_offspring_pool_and_immigrant_pool', hab.hab_dormancy_process_from_offspring_pool_and_immigrant_pool)
    minted = sp.minted_num
    step('hab_sex_reproduce_calculation_into_offspring_marker_pool', lambda: hab.hab_sex_reproduce_calculation_into_offspring_marker_pool('patch1', 1.0))
    assert sp.minted_num == minted
    hab.hab_dead_selection(0.1, 0.5, 'niche_gaussian')
    step('hab_sexual_reprodece_germinate', lambda: hab.hab_sexual_reprodece_germinate(1.0))
    dorm = list(hab.dormancy_pool)
    step('habitat_disturbance_process', hab.habitat_disturbance_process)
    assert hab.indi_num == 0 and len(hab.dormancy_pool) == len(dorm) and all(a is b for a, b in zip(hab.dormancy_pool, dorm))
    CURRENT[0] = 'germinate(is_remove=False)'
    try: hab.hab_local_germinate_from_offspring_and_immigrant_pool(is_remove=False); raise AssertionError('no NotImplementedError')
    except NotImplementedError: pass
    step('teardown (recycle_everything)', lambda: recycle_everything(hab))
    assert sp.minted_num == sp.recycled_num
    sp.purge_recycling_pool()


if __name__ == '__main__':
    out = {'steps': STEPS}
    try:
        main(); out['ok'] = True
    except BaseException:
        out.update(ok=False, failed_step=CURRENT[0], traceback=traceback.format_exc())
        print('SMOKE FAILED at %s\n%s' % (CURRENT[0], out['traceback']))
    out['seconds'] = round(time.time() - T0, 2)
    save('smoke', out)
    if not out['ok']: sys.exit(1)
    print('SMOKE PASSED (%.2f s)' % out['seconds'])
