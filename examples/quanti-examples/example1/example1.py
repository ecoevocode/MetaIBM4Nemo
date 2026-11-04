''' example1 -- four species, four patches, MetaIBM's life cycle on Nemo's genetics (metaibm-nemo).
    Founders are built from Nemo's genetic parameters alone -- no allele is edited by hand; each species lives alone in its own patch.
    Each species = one spN.json (its values, then a description of the species); quanti_parameters.md explains every row.
    One species = one metaibm.indfactory (= one Nemo boot sim + IndFactory): Nemo mints, breeds, clones and recycles the individuals in C++,
    MetaIBM keeps the landscape, the selection and the bookkeeping.
    Needs pynemo with Individual.make_copy, IndFactory.makeClonedOffsprg and TTrait.set_allele_bit (e.g. the conda env metaibm-nemo-test3).
'''
import json
import os
import random
from itertools import product

import numpy as np

import bootstrap_metaibm                    # puts the project root on sys.path: the metaibm + extension packages; pynemo comes from the conda env
import metaibm
from metaibm.indfactory import indfactory

from pynemo import _nemo
from pynemo.params import SimParameters

HERE = os.path.dirname(os.path.abspath(__file__))
NEMO_CACHE, METAIBM_OUT = 'nemo_cache', 'metaibm_output'    # relative to HERE: main() chdirs there first

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
_UNSET = '_UNSET'   # placeholder in sp*.json: the row stays in the file, the key never reaches Nemo.
                    # Not null: Nemo would read a bare key as the BOOL flag `k 1`.

def load_species(species_id):
    ''' one species' parameters from <species_id>.json (its description is not read; quanti_parameters.md explains every row)
        RETURN: nemo_params (SimParameters: non_genetic + quanti rows), metaibm_nemo_params (dict) -- rows set to _UNSET are left out '''
    with open(os.path.join(HERE, species_id + '.json')) as f:
        cfg = json.load(f)

    nemo_params = SimParameters()
    for k, v in list(cfg['non_genetic'].items()) + list(cfg['quanti'].items()):
        if v != _UNSET:
            nemo_params.add(k, v)

    metaibm_nemo_params = {k: v for k, v in cfg['metaibm_nemo_params'].items() if v != _UNSET}
    return nemo_params, metaibm_nemo_params

''' MetaIBM Landscape - one metacommunity: four patches, one 30 x 30 habitat each '''
def build_metacommunity():
    ''' RETURN: an empty metacommunity -- environments drawn, no individuals yet '''
    meta_obj = metaibm.metacommunity(metacommunity_name='metacommunity')

    patch_obj_1 = metaibm.patch(patch_name='patch1', patch_index=0, location=(0, 0))
    patch_obj_2 = metaibm.patch(patch_name='patch2', patch_index=1, location=(0, 1))
    patch_obj_3 = metaibm.patch(patch_name='patch3', patch_index=2, location=(1, 0))
    patch_obj_4 = metaibm.patch(patch_name='patch4', patch_index=3, location=(1, 1))

    # global_index: unique across the whole landscape, it becomes the Nemo homepatch of every individual born in that habitat (0-65535)
    patch_obj_1.add_habitat(hab_name='h1', hab_index=0, hab_location=(0, 0), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.3, 0.3], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=0)
    patch_obj_2.add_habitat(hab_name='h2', hab_index=0, hab_location=(0, 1), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.5, 0.5], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=1)
    patch_obj_3.add_habitat(hab_name='h3', hab_index=0, hab_location=(1, 0), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.7, 0.7], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=2)
    patch_obj_4.add_habitat(hab_name='h4', hab_index=0, hab_location=(1, 1), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.9, 0.9], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=3)

    meta_obj.add_patch(patch_name='patch1', patch_object=patch_obj_1)
    meta_obj.add_patch(patch_name='patch2', patch_object=patch_obj_2)
    meta_obj.add_patch(patch_name='patch3', patch_object=patch_obj_3)
    meta_obj.add_patch(patch_name='patch4', patch_object=patch_obj_4)
    return meta_obj

def fill_founders(meta_obj, patch_2_factory, reproduce_mode):
    ''' one species per patch: patch_2_factory = {patch_id: the indfactory that mints that patch's founders};
        each founder comes out of make_founder() as its species' genetics make it (its spN.json, + quanti_init_freq for sp1) -- nothing is set by hand '''
    for patch_id, patch_obj in meta_obj.set.items():
        for h_object in patch_obj.set.values():
            for row, col in product(range(h_object.length), range(h_object.width)):
                gender = 'female' if reproduce_mode == 'asexual' else np.random.choice(('male', 'female'))
                founder = patch_2_factory[patch_id].make_founder(gender, home_habitat=h_object.global_index)
                h_object.add_individual(founder, row, col)

def save_snapshot(meta_obj, columns, time_step, sp_file, pheno_files):
    ''' one landscape snapshot into each file: species ids, then one file per phenotype '''
    mode = 'w' if time_step == 0 else 'a'    # 'w' writes the header, 'a' appends a gzip member
    index = ['time_step%d' % time_step]
    meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_sp_id_values(), file_name=sp_file, index=index, columns=columns, mode=mode)
    for trait_name, pheno_file in pheno_files.items():
        meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_phenotype_values(trait_name=trait_name), file_name=pheno_file, index=index, columns=columns, mode=mode)

''' MetaIBM life cycle: selection -> reproduction -> germination (every step recycles its dead and its ungerminated offspring into Nemo) '''
def meta_reproduce(meta_obj, reproduce_mode, asexual_birth_rate, sexual_birth_rate):
    ''' whichever mode the run is configured for, into the offspring_pool; RETURN: offspring born '''
    if   reproduce_mode == 'asexual': meta_obj.meta_asex_reproduce_mutate_into_offspring_pool(asexual_birth_rate)                    # indfactory.make_clone()
    elif reproduce_mode == 'sexual':  meta_obj.meta_sex_reproduce_mutate_into_offspring_pool(sexual_birth_rate)                      # indfactory.make_offspring()
    elif reproduce_mode == 'mix':     meta_obj.meta_mix_reproduce_mutate_into_offspring_pool(asexual_birth_rate, sexual_birth_rate)  # both, split by fitness at selection
    else: raise ValueError("reproduce_mode must be 'asexual', 'sexual' or 'mix', got %r" % (reproduce_mode,))
    return meta_obj.meta_offspring_pool_individual_num()

def recycle_everything(meta_obj, factories):
    ''' teardown: recycle the last generation, then free each species' recycling pool '''
    meta_obj.meta_clear_up_microsite_individuals()
    meta_obj.meta_clear_up_offspring_and_immigrant_pool()
    meta_obj.meta_clear_up_dormancy_pool()
    for factory in factories:
        assert factory.minted_num == factory.recycled_num, '%s: minted %d != recycled %d' % (factory.species_id, factory.minted_num, factory.recycled_num)   # every individual went back to Nemo
        factory.purge_recycling_pool()
    print('done: minted == recycled for %s' % ', '.join('%s (%d)' % (f.species_id, f.minted_num) for f in factories), flush=True)

def main(reproduce_mode='sexual', all_time_step=2000, seed=None):
    ''' one run: boot the four species, fill the landscape with founders, run the life cycle all_time_step times, recycle everything
        INPUT:  reproduce_mode = 'asexual' (indfactory.make_clone) | 'sexual' (indfactory.make_offspring) | 'mix' (both); all_time_step = int;
                seed = int (seeds Python's random, numpy and every indfactory: the same seed reruns the same script identically) | None (a time-based seed)
        RETURN: the output folder -- the run lives on disk, not in memory '''
    os.chdir(HERE)
    os.makedirs(NEMO_CACHE, exist_ok=True)      # printLogHeader() opens `logfile` raw, no mkdir of its own
    os.makedirs(METAIBM_OUT, exist_ok=True)

    ''' MetaIBM parameters (reproduce_mode and all_time_step are main()'s arguments) '''
    base_dead_rate, fitness_wid = 0.1, 0.5
    asexual_birth_rate, sexual_birth_rate = 0.5, 1          # <= 1: parents are sampled without replacement
    record_every, print_every = 10, 100                     # a landscape snapshot every record_every steps, a progress line every print_every steps

    ''' MetaIBM output stream - one row per recorded step, one column per microsite '''
    sp_file = os.path.join(METAIBM_OUT, 'sp_dis_over_time-%s.csv.gz'%(reproduce_mode))
    pheno_files = {'phenotype-1': os.path.join(METAIBM_OUT, 'pheno1_dis_over_time-%s.csv.gz'%(reproduce_mode)),
                   'phenotype-2': os.path.join(METAIBM_OUT, 'pheno2_dis_over_time-%s.csv.gz'%(reproduce_mode))}

    ''' Nemo parameters '''
    # 2 traits per species; founder phenotypes (0.1, 0.1) / (0.3, 0.3) / (0.5, 0.5) / (0.7, 0.7) built from each species' genetics alone (sp1.json - sp4.json):
    # sp1 diallelic + quanti_init_freq, sp2 diallelic_asymmetrical + dominance, sp3 continuous + epistasis, sp4 continuous_in_place + full pleiotropy
    params1, metaibm_params1 = load_species('sp1')
    params2, metaibm_params2 = load_species('sp2')
    params3, metaibm_params3 = load_species('sp3')
    params4, metaibm_params4 = load_species('sp4')

    ''' Indfactory backended by Nemo: one species = one boot sim = one trait prototype = one Nemo IndFactory '''
    _nemo.cvar.SILENT_RUN = True                # the first boot prints its seed line before it reads run_mode 'silent'; process-global, stays on
    random.seed(seed); np.random.seed(seed)     # MetaIBM draws from both (habitat.py uses Python's random); None: fresh seeds
    
    sp1_factory = indfactory('sp1', pheno_name_ls=('phenotype-1', 'phenotype-2'), nemo_params=params1, trait_keys=('quant',), metaibm_nemo_params=metaibm_params1, seed=seed)
    sp2_factory = indfactory('sp2', pheno_name_ls=('phenotype-1', 'phenotype-2'), nemo_params=params2, trait_keys=('quant',), metaibm_nemo_params=metaibm_params2, seed=seed)
    sp3_factory = indfactory('sp3', pheno_name_ls=('phenotype-1', 'phenotype-2'), nemo_params=params3, trait_keys=('quant',), metaibm_nemo_params=metaibm_params3, seed=seed)
    sp4_factory = indfactory('sp4', pheno_name_ls=('phenotype-1', 'phenotype-2'), nemo_params=params4, trait_keys=('quant',), metaibm_nemo_params=metaibm_params4, seed=seed)
    factories = (sp1_factory, sp2_factory, sp3_factory, sp4_factory)

    meta_obj = build_metacommunity()
    fill_founders(meta_obj, {'patch1': sp1_factory, 'patch2': sp2_factory, 'patch3': sp3_factory, 'patch4': sp4_factory}, reproduce_mode)
    columns = list(meta_obj.columns_patch_habitat_microsites_id())

    for time_step in range(all_time_step):
        alive = meta_obj.get_meta_individual_num()
        meta_obj.meta_dead_selection(base_dead_rate, fitness_wid, method='niche_gaussian')         # the dead are recycled inside
        dead = alive - meta_obj.get_meta_individual_num()
        born = meta_reproduce(meta_obj, reproduce_mode, asexual_birth_rate, sexual_birth_rate)
        meta_obj.meta_local_germinate_from_offspring_and_immigrant_pool(is_remove=True)
        ungerminated = meta_obj.meta_clear_up_offspring_and_immigrant_pool()                      # recycles the offspring left without a microsite

        if time_step % record_every == 0:
            save_snapshot(meta_obj, columns, time_step, sp_file, pheno_files)

        if time_step % print_every == 0:
            print('t=%4d  n=%5d  dead=%4d born=%4d ungerminated=%4d'% (time_step, meta_obj.get_meta_individual_num(), dead, born, ungerminated), flush=True)

    recycle_everything(meta_obj, factories)
    return os.path.join(HERE, METAIBM_OUT)


if __name__ == '__main__':
    main(reproduce_mode='asexual', all_time_step=1000, seed=1)
    main(reproduce_mode='sexual', all_time_step=1000, seed=1)
