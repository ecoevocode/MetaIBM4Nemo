''' test-non-quanti-neutral-selection -- a species WITHOUT any Nemo trait lives through complete time steps of one habitat under neutral selection.

    Species: one indfactory, trait_keys=(), pheno_name_ls=() -- no genetic rows in nemo_params: no quanti, no phenotype, no other trait.
    Founders: filled by hand (make_founder + add_individual, as fill_founders() in examples/quanti-examples/example1), so no species_2_phenotype_ls is needed.
    One time step: hab_dead_selection(method='neutral_uniform') -> reproduce into offspring_pool -> local germination -> recycle the ungerminated offspring.
    Checked for reproduce_mode = 'asexual', 'sexual' ('mix' is not tested):
      (1) the founders fill the habitat; bookkeeping is consistent; no individual carries a trait
      (2) every time step: bookkeeping is consistent, the habitat is full again, every resident is alive and of this species;
          asexual: every offspring is a clone (father_id 0), sexual: none is
      (3) deaths over all time steps match base_dead_rate (binomial, 4 sd)
      (4) niche_gaussian on this species raises KeyError('quant') before it kills or recycles anyone
      (5) teardown: minted == recycled
    Run: python test-non-quanti-neutral-selection.py   (conda env metaibm-nemo-test6); exit code 1 if any check fails.
'''
import random
import sys
import tempfile
from itertools import product
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]      # test/test-nemoback-habitat/non-quanti-selection -> the project root
sys.path.insert(0, str(PROJECT_ROOT))
from metaibm.habitat import habitat                     # noqa: E402
from metaibm.indfactory import indfactory               # noqa: E402

from pynemo import _nemo                                # noqa: E402
from pynemo.params import SimParameters                 # noqa: E402

SEED = 1
LENGTH, WIDTH = 20, 20
TIME_STEPS = 10
BASE_DEAD_RATE, FITNESS_WID = 0.1, 0.5
ASEXUAL_BIRTH_RATE, SEXUAL_BIRTH_RATE = 0.5, 1
NEMO_CACHE = tempfile.mkdtemp(prefix='nemo_cache_')     # Nemo's boot logs stay out of the repo

FAILED = []


def check(ok, label):
    print('%s  %s' % ('PASS' if ok else 'FAIL', label), flush=True)
    if not ok: FAILED.append(label)


def traitless_species(species_id):
    ''' one species without any Nemo trait: only the rows Nemo needs to boot; RETURN: indfactory '''
    params = SimParameters()
    rows = [('replicates', 1), ('generations', 1), ('patch_number', 1), ('patch_capacity', 1), ('aging', 1),
            ('run_mode', 'silent'), ('root_dir', NEMO_CACHE), ('logfile', NEMO_CACHE + '/nemo.log'), ('filename', species_id)]
    for key, value in rows: params.add(key, value)
    return indfactory(species_id, (), params, trait_keys=(), seed=SEED)


def make_habitat():
    ''' one 20 x 20 habitat with one environment type: neutral selection never reads it '''
    return habitat(hab_name='h1', hab_index=0, hab_location=(0, 0), num_env_types=1, env_types_name=('env1',), mean_env_ls=[0.5], var_env_ls=[0.0],
                   length=LENGTH, width=WIDTH, global_index=0)


def fill_founders(h_object, factory, reproduce_mode):
    ''' every microsite gets a founder (as examples/quanti-examples/example1 fill_founders): asexual -> all female, else a random gender '''
    for row, col in product(range(h_object.length), range(h_object.width)):
        gender = 'female' if reproduce_mode == 'asexual' else np.random.choice(('male', 'female'))
        h_object.add_individual(factory.make_founder(gender, home_habitat=h_object.global_index), row, col)


def residents(h_object):
    ''' RETURN: {(row, col): individual} of the occupied microsites '''
    return {(r, c): ind for r, row in enumerate(h_object.set['microsite_individuals']) for c, ind in enumerate(row) if ind is not None}


def bookkeeping_errors(h_object):
    ''' grid, indi_num, occupied / empty lists and species_category must describe the same individuals; RETURN: list of str (empty = consistent) '''
    grid = residents(h_object)
    errors = []
    if h_object.indi_num != len(grid): errors.append('indi_num %d != %d on the grid' % (h_object.indi_num, len(grid)))
    if set(h_object.occupied_site_pos_ls) != set(grid) or len(h_object.occupied_site_pos_ls) != len(grid): errors.append('occupied_site_pos_ls != grid')
    if set(h_object.empty_site_pos_ls) != set(product(range(h_object.length), range(h_object.width))) - set(grid): errors.append('empty_site_pos_ls != free microsites')
    if len(h_object.empty_site_pos_ls) + len(h_object.occupied_site_pos_ls) != h_object.size: errors.append('empty + occupied != size')
    category = {(sp_id, gender, pos) for sp_id, by_gender in h_object.species_category.items() for gender, pos_ls in by_gender.items() for pos in pos_ls}
    if category != {(ind.species_id, ind.gender, pos) for pos, ind in grid.items()}: errors.append('species_category != grid (species, gender, position)')
    return errors


def reproduce(h_object, reproduce_mode):
    ''' RETURN: offspring born into offspring_pool '''
    if reproduce_mode == 'asexual': return h_object.hab_asex_reproduce_mutate_into_offspring_pool(ASEXUAL_BIRTH_RATE)
    return h_object.hab_sex_reproduce_mutate_into_offspring_pool(SEXUAL_BIRTH_RATE)


def run_mode(reproduce_mode):
    print('\n== reproduce_mode = %s' % reproduce_mode, flush=True)
    random.seed(SEED); np.random.seed(SEED)
    species_id = 'sp1_%s' % reproduce_mode
    factory = traitless_species(species_id)
    h_object = make_habitat()

    # (1) founders
    fill_founders(h_object, factory, reproduce_mode)
    founders = residents(h_object)
    check(len(founders) == h_object.size and h_object.indi_num == h_object.size and not h_object.empty_site_pos_ls, '(1) the founders fill all %d microsites' % h_object.size)
    check(not bookkeeping_errors(h_object), '(1) bookkeeping consistent after filling %s' % bookkeeping_errors(h_object))
    check(factory.traits_index == {} and all(ind.traits_num == 0 for ind in founders.values()), '(1) no trait: traits_index = %r, traits_num = 0' % (factory.traits_index,))
    if reproduce_mode == 'asexual': check(all(ind.gender == 'female' for ind in founders.values()), '(1) asexual: every founder is female')

    # (2) complete time steps
    exposed = died = 0
    for t in range(TIME_STEPS):
        before = h_object.indi_num
        dead = h_object.hab_dead_selection(BASE_DEAD_RATE, FITNESS_WID, 'neutral_uniform')
        after_selection = h_object.indi_num
        born = reproduce(h_object, reproduce_mode)
        clones = sum(ind.father_id == 0 for ind in h_object.offspring_pool)
        germinated = h_object.hab_local_germinate_from_offspring_and_immigrant_pool(is_remove=True)
        ungerminated = h_object.hab_clear_up_offspring_pool()
        exposed += before; died += dead

        step_errors = bookkeeping_errors(h_object)
        if after_selection != before - dead: step_errors.append('indi_num after selection %d != %d - %d' % (after_selection, before, dead))
        if clones != (born if reproduce_mode == 'asexual' else 0): step_errors.append('%d clones among %d %s offspring' % (clones, born, reproduce_mode))
        if germinated != min(dead, born): step_errors.append('germinated %d != min(dead %d, born %d)' % (germinated, dead, born))
        if ungerminated != born - germinated: step_errors.append('ungerminated %d != born %d - germinated %d' % (ungerminated, born, germinated))
        if h_object.offspring_pool: step_errors.append('offspring_pool not empty after clear-up')
        if born >= dead and h_object.indi_num != h_object.size: step_errors.append('habitat not full again: %d / %d' % (h_object.indi_num, h_object.size))
        for ind in residents(h_object).values():
            if not ind.is_alive: step_errors.append('a resident was recycled'); break
            if ind.species_id != species_id: step_errors.append('a resident of another species'); break
            if reproduce_mode == 'asexual' and ind.gender != 'female': step_errors.append('a male in an asexual run'); break
        check(not step_errors, '(2) t=%d  dead=%3d born=%3d (clones %3d) germinated=%3d ungerminated=%3d  n=%d %s' % (t, dead, born, clones, germinated, ungerminated, h_object.indi_num, step_errors))

    # (3) neutral deaths: binomial(exposed, base_dead_rate)
    expected, sd = exposed * BASE_DEAD_RATE, np.sqrt(exposed * BASE_DEAD_RATE * (1 - BASE_DEAD_RATE))
    check(abs(died - expected) <= 4 * sd, '(3) deaths %d over %d exposures, expected %.0f +- %.0f (4 sd)' % (died, exposed, expected, 4 * sd))

    # (4) niche_gaussian needs a phenotype: KeyError('quant') before anyone dies
    n_before, recycled_before = h_object.indi_num, factory.recycled_num
    try:
        h_object.hab_dead_selection(BASE_DEAD_RATE, FITNESS_WID, 'niche_gaussian')
        error = None
    except KeyError as e:
        error = e
    check(error is not None and error.args == ('quant',), "(4) niche_gaussian -> KeyError('quant'), got %r" % (error,))
    check(h_object.indi_num == n_before and factory.recycled_num == recycled_before and not bookkeeping_errors(h_object), '(4) ... and nobody died or was recycled')

    # (5) teardown
    h_object.hab_clear_up_microsite_individuals()
    h_object.hab_clear_up_offspring_pool(); h_object.hab_clear_up_immigrant_pool(); h_object.hab_clear_up_dormancy_pool()
    check(factory.minted_num == factory.recycled_num, '(5) teardown: minted %d == recycled %d' % (factory.minted_num, factory.recycled_num))
    factory.purge_recycling_pool()


if __name__ == '__main__':
    _nemo.cvar.SILENT_RUN = True
    for mode in ('asexual', 'sexual'):
        run_mode(mode)
    print('\n%s: %d check(s) failed' % ('FAILED' if FAILED else 'ALL PASSED', len(FAILED)) + (''.join('\n  - ' + f for f in FAILED)), flush=True)
    sys.exit(1 if FAILED else 0)
