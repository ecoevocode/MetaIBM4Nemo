# Smoke (< 30 s): imports + one call of every NemoBackend patch function on real Nemo individuals, invariants after each. run_all.py stops if this fails.
import sys, time, traceback, warnings
import numpy as np
from common import make_patch, habs, check_invariants, recycle_everything, nemo_species, seed_all, save

warnings.filterwarnings('ignore', category=RuntimeWarning); warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
T0 = time.time()
STEPS, CURRENT = [], ['import']


def main():
    seed_all(0)
    F = {'sp1': nemo_species('sp1', 0.2, seed=1), 'sp2': nemo_species('sp2', 0.5, seed=2)}
    fs = list(F.values())
    np.random.seed(0)
    p = make_patch(3, L=5, K=10, gidx0=3, means=[(0.2, 0.2), (0.5, 0.5), (0.2, 0.2)])

    def step(name, fn):
        CURRENT[0] = name
        r = fn()
        STEPS.append({'step': name, 'return': repr(r), 'invariants': check_invariants(p, fs, name)})
        print('[smoke] %-70s -> %r' % (name, r))
        return r

    step('patch_initialize', lambda: p.patch_initialize(F, [[0.2, 0.2], [0.5, 0.5]]))
    assert all(x.nemo_indi_obj.getHome() == h.global_index for h in habs(p) for row in h.set['microsite_individuals'] for x in row)
    step('patch_dead_selection', lambda: p.patch_dead_selection(0.1, 0.5, 'niche_gaussian'))
    step('patch_sex_reproduce_mutate_into_offspring_pool', lambda: p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None))
    assert all(k.nemo_indi_obj.getHome() == h.global_index for h in habs(p) for k in h.offspring_pool)
    step('patch_dispersal_within_from_offspring_pool_to_immigrant_pool', lambda: p.patch_dispersal_within_from_offspring_pool_to_immigrant_pool(0.5))
    step('patch_local_germinate_from_offspring_and_immigrant_pool', p.patch_local_germinate_from_offspring_and_immigrant_pool)
    step('patch_dormancy_process_from_offspring_pool_and_immigrant_pool', p.patch_dormancy_process_from_offspring_pool_and_immigrant_pool)
    step('patch_dead_selection', lambda: p.patch_dead_selection(0.3, 0.5, 'niche_gaussian'))
    step('patch_sex_reproduce_mutate_into_offspring_pool', lambda: p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None))
    step('patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool', lambda: p.patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool(0.5))
    step('patch_local_germinate_from_offspring_immigrant_and_dormancy_pool', p.patch_local_germinate_from_offspring_immigrant_and_dormancy_pool)
    step('patch_clear_up_offspring_pool', p.patch_clear_up_offspring_pool)
    step('patch_clear_up_immigrant_pool', p.patch_clear_up_immigrant_pool)
    step('patch_dead_selection', lambda: p.patch_dead_selection(0.3, 0.5, 'niche_gaussian'))
    step('patch_sex_reproduce_mutate_into_offspring_pool', lambda: p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None))
    step('patch_dipersal_within_from_offspring_and_dormancy_pool', lambda: p.patch_dipersal_within_from_offspring_and_dormancy_pool(0.5))
    step('patch_local_germinate_from_offspring_and_dormancy_pool', p.patch_local_germinate_from_offspring_and_dormancy_pool)
    step('patch_dormancy_process_from_offspring_pool_to_dormancy_pool', p.patch_dormancy_process_from_offspring_pool_to_dormancy_pool)
    minted = sum(f.minted_num for f in fs)
    step('patch_sex_reproduce_calculation_into_offspring_marker_pool', lambda: p.patch_sex_reproduce_calculation_into_offspring_marker_pool(1.0))
    step('patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool', lambda: p.patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool(0.5, True))
    assert sum(f.minted_num for f in fs) == minted
    step('patch_clear_up_offspring_marker_pool', p.patch_clear_up_offspring_marker_pool)
    step('patch_clear_up_immigrant_marker_pool', p.patch_clear_up_immigrant_marker_pool)
    step('patch_dead_selection', lambda: p.patch_dead_selection(0.3, 0.5, 'niche_gaussian'))
    step('patch_sexual_birth_germinate', lambda: p.patch_sexual_birth_germinate(1.0, None, None))
    dorm = [list(h.dormancy_pool) for h in habs(p)]
    step('patch_disturbance_process', p.patch_disturbance_process)
    assert all(h.indi_num == 0 and len(h.dormancy_pool) == len(d) and all(a is b for a, b in zip(h.dormancy_pool, d)) for h, d in zip(habs(p), dorm))
    step('patch_initialize (re-initialize, dormancy kept)', lambda: p.patch_initialize(F, [[0.2, 0.2], [0.5, 0.5]]))
    step('patch_clear_up_dormancy_pool', p.patch_clear_up_dormancy_pool)
    step('patch_clear_up_microsite_individuals', p.patch_clear_up_microsite_individuals)
    CURRENT[0] = 'dispersal(is_remove=False)'
    for fn in (p.patch_dispersal_within_from_offspring_pool_to_immigrant_pool, p.patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool, p.patch_dipersal_within_from_offspring_and_dormancy_pool):
        try: fn(0.5, False); raise AssertionError('no NotImplementedError')
        except NotImplementedError: pass
    step('teardown (recycle_everything)', lambda: recycle_everything(p))
    assert all(f.minted_num == f.recycled_num for f in fs)
    for f in fs: f.purge_recycling_pool()


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
