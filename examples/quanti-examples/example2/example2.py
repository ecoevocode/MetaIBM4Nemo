''' example2 -- mainland-island: example1's four-patch landscape, now empty at the start and colonized from a mainland that holds one species (sp1).
    The species uses Nemo's 'continuous' allele model with a purely polygenic architecture (20 het loci +-0.0125).
    One species = one metaibm.indfactory (= one Nemo boot sim + IndFactory); individuals are minted, bred, cloned and recycled in C++ by Nemo,
    while MetaIBM keeps the landscape, the selection and the bookkeeping. Needs pynemo with Individual.make_copy and IndFactory.makeClonedOffsprg
    (e.g. the conda env metaibm-nemo-test3).
'''
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
_UNSET = object()   # placeholder: the parameter keeps its line in the table but never reaches Nemo.
                    # Not None -- SimParameters.add(k, None) writes a bare key (pynemo/params.py:42)
                    # and the parser rewrites that to `k 1` (src/paramsparser.cc:1470-1473), which is
                    # how a BOOL flag is spelled. Nemo's "not set" is the key being absent.

def nemo_params(sp_log_name, p1=0.5, p2=0.5):
    ''' one species' Nemo config: 
        (1) non-genetic, only there to build the boot sim, and 
        (2) the genetics, which is what MetaIBM reads back off the prototype. 
    '''
    params = SimParameters()

    ''' (1) non-genetic: builds the boot sim, which nothing here steps -- measured bit-identical either way '''
    non_genetic = [
        ('replicates', 1),                                           # OFF REQ FIX  boot sim only; without it the consistency check fails, no record is built, and no prototype comes out; e.g. 1
        ('generations', 1),                                          # OFF REQ FIX  boot sim only; same failure mode as replicates; e.g. 1
        ('patch_number', 1),                                         # OFF REQ FIX  boot sim only; without it 'param "patch_number" is missing!'; e.g. 1
        ('patch_capacity', 1),                                       # OFF REQ FIX  boot sim only; without it 'population parameters are not properly set!' + TMatrix overflow; e.g. 1
        ('aging', 1),                                                # OFF REQ FIX  boot sim only; To avoid NEMO's warning >>>> WARNING::Metapop::setPopulation: generation 0 is empty!!
        ('run_mode', 'silent'),                                      # OFF OPT FIX  file placement + console only; 'silent' turns SILENT_RUN on (no 'setting random seed' / 'outputs:' lines; warnings muted too, ***ERROR*** still shows), the .log files are still written. Measured vs 'overwrite' with seed 42: founders identical, logs differ only in start time and this line; e.g. 'overwrite' to see Nemo's warnings again
        ('root_dir', NEMO_CACHE),                                    # OFF OPT FIX  file placement only, but the directory IS created at boot; e.g. 'nemo_cache'
        ('logfile', NEMO_CACHE + '/nemo.log'),                       # OFF OPT FIX  file placement only; No effect in metaibm-nemo; e.g. 'nemo_cache/nemo.log'
        ('filename', sp_log_name),                                   # OFF REQ VAR! names this species' .log -- the only one that is per species; sharing it does not error, the species just overwrite each other's log; e.g. 'sp1'
    ]
    for k, v in non_genetic:
        params.add(k, v)

    ''' (2) the genetics: trait_keys=('quant',) '''
    quanti = [
        # -- architecture ------------------------------------------------------------------------ #
        ('quanti_traits', 2),                                        # ON  REQ FIX* quantitative traits per individual.
        ('quanti_loci', 40),                                         # ON  REQ FIX  loci summed over all traits
        ('quanti_loci_per_trait', '{{20, 20}}'),                     # ON  OPT VAR  loci of each trait, one row of quanti_traits columns
        # -- allele model: continuous (a mutation ADDS to the allele value) ---------------------- #
        ('quanti_allele_model', 'continuous'),                       # ON  REQ VAR  Nemo's continuous allele model
        ('quanti_allele_value', _UNSET),                             # ON  OPT VAR  di-allelic models only, and mandatory there: a scalar a gives the alleles +a / -a, or a locus-specific matrix; e.g. 0.01
        ('quanti_diallele_datatype', _UNSET),                        # ON  OPT VAR  di-allelic only: how Nemo stores the sequence, 'bit' (default) | 'byte'; same individuals either way, except next to quanti_dominance_effects, where 'bit' gives a different phenotype -- dominance on -> 'byte'; e.g. 'byte'
        # -- founder initialisation -------------------------------------------------------------- #
        ('quanti_init_model', 2),                                    # ON  OPT VAR  2 Gaussian jitter per allele
        ('quanti_init_value', '{{%g, %g}}' % (p1, p2)),              # ON  OPT VAR  genotypic-value centre per trait
        ('quanti_init_variance', '{{0.0014, 0.0014}}'),              # ON  OPT VAR  initial genotypic VARIANCE per trait
        ('quanti_init_freq', _UNSET),                                # OFF - but MetaIBM rebuilt it for diallelic(-asymmetric) models without pleiotropy in indfactory.py: a parameter of Nemo's LCE quanti_init, and MetaIBM4Nemo never reach Nemo's LCE. 
        ('quanti_init_trait_values', _UNSET),                        # OFF - the other half of LCE quanti_init: per-patch trait centres, never read in metaibm-nemo; the per-species centre is quanti_init_value, which IS read at founding
        # -- genotype to phenotype --------------------------------------------------------------- #
        ('quanti_environmental_variance', '{{0.000625, 0.000625}}'), # ON  OPT VAR  Ve, non-heritable part: genotypic value + N(0, sqrt(Ve))
        ('quanti_heritability', _UNSET),                             # OFF - a target heritability; only Nemo's LCEs (quanti_modifier, breed_quanti) use it, to reset Ve, and no LCE runs in metaibm-nemo -- set Ve itself with quanti_environmental_variance
        ('quanti_heritability_setTime', _UNSET),                     # OFF - how often (in generations) LCE breed_quanti re-applies quanti_heritability; never read in metaibm-nemo
        ('quanti_heritability_isBroad', _UNSET),                     # OFF - switches the heritability target from narrow (Va) to broad sense (Vg); read only by the quanti stat handler, which never runs in metaibm-nemo
        ('quanti_modifier_at_generation', _UNSET),                   # OFF - LCE quanti_modifier: the generation(s), one row, at which it resets Ve to reach quanti_heritability; no LCE runs in metaibm-nemo
        ('quanti_modifier_at_age', _UNSET),                          # OFF - LCE quanti_modifier: the age class Va is measured on, 'adults' | 'offspring' | 'all'
        ('quanti_modifier_is_permanent', _UNSET),                    # OFF - LCE quanti_modifier: turns quanti_modifier_at_generation into a period (1 = every generation)
        # -- mutation ---------------------------------------------------------------------------- #
        ('quanti_mutation_rate', 0.0001),                            # ON  REQ VAR  per ALLELE COPY per generation
        ('quanti_mutation_variance', 0.0025),                        # ON  OPT VAR  VARIANCE of one mutation's effect (continuous adds N(0, sqrt(v)))
        # -- pleiotropy -------------------------------------------------------------------------- #
        ('quanti_pleiotropy', _UNSET),                               # ON  OPT VAR  unset/'no' (default) one locus one trait
        ('quanti_mutation_correlation', _UNSET),                     # ON  OPT VAR  needs pleiotropy
        ('quanti_mutation_matrix', _UNSET),                          # ON  OPT VAR  variance-COVARIANCE matrix of mutation effects
        ('quanti_pleio_matrix', _UNSET),                             # ON  OPT VAR  loci-by-traits 0/1 connectivity
        # -- dominance --------------------------------------------------------------------------- #
        ('quanti_dominance_effects', _UNSET),                        # ON  OPT VAR  k: locus value a1 + a2 + k |a2 - a1| (k 0 additive, 1 the bigger allele dominant, > 1 overdominance), scalar or matrix (rows = traits, cols = loci, recycled); di-allelic models: also write quanti_diallele_datatype 'byte'; e.g. 0.3
        ('quanti_dominance_mean', _UNSET),                           # OFF - dead in Nemo: declared, read nowhere; draw random k with a config macro instead, quanti_dominance_effects {{rnorm(n, mean, sd)}}
        ('quanti_dominance_sd', _UNSET),                             # OFF - dead in Nemo: the other half of quanti_dominance_mean
        # -- epistasis --------------------------------------------------------------------------- #
        ('quanti_epistatic_effects', _UNSET),                        # ON  OPT VAR  one coefficient per locus pair: G += e_ij g_i g_j; this parameter switches epistasis on. Without a network it is safe only when one trait owns every locus (quanti_traits 1 or full pleiotropy); e.g. {{0.5, -0.5, 0.2}}
        ('quanti_epistatic_network', _UNSET),                        # ON  OPT VAR  the interacting loci, one network per row: M loci -> M(M-1)/2 coefficients, in the order (1,2)(1,3)..(2,3)..; ids are 1-based and counted within ONE trait, the same network acts on every trait; e.g. {{1, 2, 3}}
        # -- genetic map: process-global, keyed by trait type ------------------------------------ #
        ('quanti_recombination_rate', 0.5),                          # ON  OPT FIX* free recombination (= 0.5)
        ('quanti_genetic_map_resolution', _UNSET),                   # ON  OPT FIX  map grain in cM (default 1), read only once a map is active; process-global: Nemo keeps the smallest over all species and rescales the whole map, so keep it identical; e.g. 0.1
        ('quanti_chromosome_num_locus', _UNSET),                     # ON  OPT FIX  loci per chromosome, one row; cannot switch the map on alone -- read next to a matrix quanti_recombination_rate or quanti_random_genetic_map; e.g. {{10, 10, 10, 10}}
        ('quanti_genetic_map', _UNSET),                              # ON  OPT FIX  locus positions in cM, one row per chromosome, totalling quanti_loci; switches the map on, not together with quanti_recombination_rate; every species recombines on the FIRST registrant's map; e.g. 40 loci on 2 chromosomes: '{' + '{0, 5, ... , 95}' * 2 + '}'
        ('quanti_random_genetic_map', _UNSET),                       # ON  OPT FIX  chromosome lengths in cM, loci scattered uniformly; switches the map on, shared like quanti_genetic_map (the first registrant's map wins); e.g. {{100, 100}}
        # -- output: Nemo's FileServices, forbiddened by MetaIBM ---------------------------------- #
        ('quanti_output', _UNSET),                                   # FORBIDDEN BY METAIBM: ON  - Nemo's trait file: it only ever sees the boot sim's empty Metapop, never MetaIBM's individuals, and a bad option or a missing quanti_logtime kills the boot. Options: genotypes | snp_genotypes | snp_id | plink
        ('quanti_logtime', _UNSET),                                  # FORBIDDEN BY METAIBM: OFF - how often quanti_output writes (generations); read only next to it, where it is mandatory
        ('quanti_dir', _UNSET),                                      # FORBIDDEN BY METAIBM: OFF - directory of the quanti output files, relative to root_dir
        ('quanti_freq_output', _UNSET),                              # FORBIDDEN BY METAIBM: ON  - per-locus allele frequencies, di-allelic only (inert on continuous models); same empty Metapop, and fatal without quanti_freq_logtime
        ('quanti_freq_logtime', _UNSET),                             # FORBIDDEN BY METAIBM: OFF - how often quanti_freq_output writes; read only next to it, where it is mandatory
        ('quanti_ohta_output', _UNSET),                              # FORBIDDEN BY METAIBM: ON  - Ohta's LD statistics per locus pair; fatal on anything but a di-allelic, single-trait, multi-patch setup
        ('quanti_ohta_logtime', _UNSET),                             # FORBIDDEN BY METAIBM: OFF - how often quanti_ohta_output writes; read only next to it, where it is mandatory
    ]
    for k, v in quanti:
        if v is not _UNSET:
            params.add(k, v)

    return params

''' MetaIBM Landscape - the islands: example1's metacommunity, four patches, one 30 x 30 habitat each; no founders, the species arrives from the mainland '''
def build_metacommunity():
    ''' RETURN: an empty metacommunity -- environments drawn, no individuals yet '''
    meta_obj = metaibm.metacommunity(metacommunity_name='metacommunity')

    patch_obj_1 = metaibm.patch(patch_name='patch0', patch_index=0, location=(0, 0))
    patch_obj_2 = metaibm.patch(patch_name='patch1', patch_index=1, location=(0, 1))
    patch_obj_3 = metaibm.patch(patch_name='patch2', patch_index=2, location=(1, 0))
    patch_obj_4 = metaibm.patch(patch_name='patch3', patch_index=3, location=(1, 1))

    # global_index: becomes the Nemo homepatch of every individual born in that habitat (0-65535); one per island (0-3), never 65535, which marks the mainland
    # hab_index: the index WITHIN the patch (0 .. hab_num-1), so 0 for every one-habitat island -- patch.migrants_to_patch_into_habs_immigrant_pool indexes its per-habitat split with it
    patch_obj_1.add_habitat(hab_name='h0', hab_index=0, hab_location=(0, 0), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.1, 0.1], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=0)
    patch_obj_2.add_habitat(hab_name='h1', hab_index=0, hab_location=(0, 1), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.3, 0.3], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=1)
    patch_obj_3.add_habitat(hab_name='h2', hab_index=0, hab_location=(1, 0), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.7, 0.7], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=2)
    patch_obj_4.add_habitat(hab_name='h3', hab_index=0, hab_location=(1, 1), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.9, 0.9], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=3)

    meta_obj.add_patch(patch_name='patch0', patch_object=patch_obj_1)
    meta_obj.add_patch(patch_name='patch1', patch_object=patch_obj_2)
    meta_obj.add_patch(patch_name='patch2', patch_object=patch_obj_3)
    meta_obj.add_patch(patch_name='patch3', patch_object=patch_obj_4)
    return meta_obj

''' MetaIBM Landscape - the mainland: the source pool, one patch with four 30 x 30 habitats of sp1 (3600 individuals), each habitat at the species' optimum (= quanti_init_value) '''
def build_mainland():
    ''' RETURN: an empty mainland -- environments drawn, no individuals yet '''
    mainland_obj = metaibm.metacommunity(metacommunity_name='mainland')

    patch_obj = metaibm.patch(patch_name='mainland', patch_index=0, location=(0, 0))

    # global_index 65535 (the top of Nemo's homepatch range) for all four habitats: homepatch 65535 just says 'born on the mainland'; a colonizer is moved, not re-born, so it keeps it
    patch_obj.add_habitat(hab_name='h0', hab_index=0, hab_location=(0, 0), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.5, 0.5], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=65535)
    patch_obj.add_habitat(hab_name='h1', hab_index=1, hab_location=(0, 1), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.5, 0.5], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=65534)
    patch_obj.add_habitat(hab_name='h2', hab_index=2, hab_location=(1, 0), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.5, 0.5], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=65533)
    patch_obj.add_habitat(hab_name='h3', hab_index=3, hab_location=(1, 1), num_env_types=2, env_types_name=('env1', 'env2'), mean_env_ls=[0.5, 0.5], var_env_ls=[0.025**2, 0.025**2], length=30, width=30, global_index=65532)

    mainland_obj.add_patch(patch_name='mainland', patch_object=patch_obj)
    return mainland_obj

def configure_heterozygous_founder(indi_obj):
    ''' cryptic founder: every locus heterozygous c + a / c - a (c = 0.0125, a = 0.02, sign drawn per locus), plus a small jitter per allele (Vg(F0) = 0.00125). '''
    quant_trait = indi_obj.nemo_indi_obj.getTrait(indi_obj.indfactory.traits_index['quant'])
    for trait_dim in range(2):
        pos_ls = [quant_trait.get_locus_seq_pos(loc, trait_dim) for loc in range(20)]
        signs = np.random.choice((-1.0, 1.0), size=20)
        eta = np.random.normal(0.0, np.sqrt(0.00125 / 40.0), size=(20, 2))
        for loc, pos in enumerate(pos_ls):
            quant_trait.set_allele_value(pos, 0, float(0.0125 + signs[loc] * 0.02 + eta[loc, 0]))
            quant_trait.set_allele_value(pos, 1, float(0.0125 - signs[loc] * 0.02 + eta[loc, 1]))
    quant_trait.set_value()

def fill_founders(mainland_obj, factory, reproduce_mode):
    ''' every mainland microsite gets a founder minted by factory (the one species' indfactory) '''
    for patch_id, patch_obj in mainland_obj.set.items():
        for h_id, h_object in patch_obj.set.items():
            for row, col in product(range(h_object.length), range(h_object.width)):
                gender = 'female' if reproduce_mode == 'asexual' else np.random.choice(('male', 'female'))
                founder = factory.make_founder(gender, home_habitat=h_object.global_index)
                configure_heterozygous_founder(founder)
                h_object.add_individual(founder, row, col)

def save_snapshot(meta_obj, columns, time_step, sp_file, pheno_files):
    ''' one landscape snapshot into each file: species ids, then one file per phenotype '''
    mode = 'w' if time_step == 0 else 'a'    # 'w' writes the header, 'a' appends a gzip member
    index = ['time_step%d' % time_step]
    meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_sp_id_values(), file_name=sp_file, index=index, columns=columns, mode=mode)
    for trait_name, pheno_file in pheno_files.items():
        meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_phenotype_values(trait_name=trait_name), file_name=pheno_file, index=index, columns=columns, mode=mode)

''' MetaIBM life cycle on the mainland: selection -> birth straight into its empty microsites (no pools, nothing left to recycle) '''
def mainland_birth(mainland_obj, reproduce_mode, asexual_birth_rate, sexual_birth_rate):
    ''' whichever mode the run is configured for; RETURN: offspring born '''
    alive = mainland_obj.get_meta_individual_num()
    if   reproduce_mode == 'asexual': mainland_obj.meta_mainland_asexual_birth_mutate_germinate(asexual_birth_rate)                    # indfactory.make_clone()
    elif reproduce_mode == 'sexual':  mainland_obj.meta_mainland_sexual_birth_mutate_germinate(sexual_birth_rate)                     # indfactory.make_offspring()
    elif reproduce_mode == 'mix':     mainland_obj.meta_mainland_mixed_birth_mutate_germinate(asexual_birth_rate, sexual_birth_rate) # both, split by fitness at selection
    else: raise ValueError("reproduce_mode must be 'asexual', 'sexual' or 'mix', got %r" % (reproduce_mode,))
    return mainland_obj.get_meta_individual_num() - alive

''' MetaIBM life cycle on the islands: selection -> reproduction -> dispersal among islands -> colonization -> germination (every step recycles its dead and its ungerminated offspring into Nemo) '''
def meta_reproduce(meta_obj, reproduce_mode, asexual_birth_rate, sexual_birth_rate):
    ''' whichever mode the run is configured for, into the offspring_pool; RETURN: offspring born '''
    if   reproduce_mode == 'asexual': meta_obj.meta_asex_reproduce_mutate_into_offspring_pool(asexual_birth_rate)                    # indfactory.make_clone()
    elif reproduce_mode == 'sexual':  meta_obj.meta_sex_reproduce_mutate_into_offspring_pool(sexual_birth_rate)                     # indfactory.make_offspring()
    elif reproduce_mode == 'mix':     meta_obj.meta_mix_reproduce_mutate_into_offspring_pool(asexual_birth_rate, sexual_birth_rate) # both, split by fitness at selection
    else: raise ValueError("reproduce_mode must be 'asexual', 'sexual' or 'mix', got %r" % (reproduce_mode,))
    return meta_obj.meta_offspring_pool_individual_num()

def disperse_among_islands(meta_obj, disp_among_rate):
    ''' a share disp_among_rate of every island's offspring_pool leaves for the other islands' immigrant_pools (uniform kernel: each other island equally likely),
        where it competes with the local offspring at germination; is_remove=True: an emigrant leaves its source pool (moved, not copied). RETURN: emigrants '''
    meta_obj.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(disp_among_rate, method='uniform', is_remove=True)   # 'aomng' is metaibm's own spelling
    return meta_obj.meta_immigrant_pool_individual_num()     # the immigrant pools are empty before this step: cleared at the end of every time step

def colonize(meta_obj, mainland_obj, reproduce_mode, propagules_rain_num):
    ''' propagules from random occupied mainland microsites land on empty island microsites -- straight in, ahead of the local offspring's germination;
        is_remove=True: a colonizer leaves the mainland (moved, not copied), whose next birth refills its site. RETURN: colonizers '''
    alive = meta_obj.get_meta_individual_num()
    if reproduce_mode == 'asexual': meta_obj.meta_colonize_from_propagules_rains(mainland_obj, propagules_rain_num, is_remove=True)               # one at a time: a lone clone founds a population
    else:                           meta_obj.pairwise_sexual_colonization_from_prpagules_rains(mainland_obj, propagules_rain_num, is_remove=True) # (female, male) of one species into one habitat: a lone sexual colonizer could never breed, and under 'mix' only individuals with survival_rate >= 0.897 clone
    return meta_obj.get_meta_individual_num() - alive

def recycle_everything(meta_objs, factories):
    ''' teardown: recycle the last generation of the islands and the mainland, then free each species' recycling pool '''
    for meta_obj in meta_objs:
        meta_obj.meta_clear_up_microsite_individuals()
        meta_obj.meta_clear_up_offspring_and_immigrant_pool()
        meta_obj.meta_clear_up_dormancy_pool()
    for factory in factories:
        assert factory.minted_num == factory.recycled_num, '%s: minted %d != recycled %d' % (factory.species_id, factory.minted_num, factory.recycled_num)   # every individual went back to Nemo
        factory.purge_recycling_pool()
    print('done: minted == recycled for %s' % ', '.join('%s (%d)' % (f.species_id, f.minted_num) for f in factories), flush=True)

def main(reproduce_mode='sexual', all_time_step=2000, seed=None):
    ''' one run: boot the species, fill the mainland with founders, run the life cycle all_time_step times, recycle everything
        INPUT:  reproduce_mode = 'asexual' (indfactory.make_clone) | 'sexual' (indfactory.make_offspring) | 'mix' (both); all_time_step = int;
                seed = int (seeds Python's random, numpy and every indfactory: the same seed reruns the same script identically) | None (a time-based seed)
        RETURN: the output folder -- the run lives on disk, not in memory '''
    os.chdir(HERE)
    os.makedirs(NEMO_CACHE, exist_ok=True)      # printLogHeader() opens `logfile` raw, no mkdir of its own
    os.makedirs(METAIBM_OUT, exist_ok=True)

    ''' MetaIBM parameters (reproduce_mode and all_time_step are main()'s arguments) '''
    base_dead_rate, fitness_wid = 0.1, 0.5
    asexual_birth_rate, sexual_birth_rate = 0.5, 1          # <= 1: parents are sampled without replacement
    disp_among_rate = 0.1                                   # share of every island's offspring that disperses to the other islands
    record_every, print_every = 10, 100                   # a landscape snapshot every record_every steps, a progress line every print_every steps

    ''' MetaIBM output stream - the islands only: one row per recorded step, one column per microsite '''
    sp_file = os.path.join(METAIBM_OUT, 'sp_dis_over_time-%s.csv.gz'%(reproduce_mode))
    pheno_files = {'phenotype-1': os.path.join(METAIBM_OUT, 'pheno1_dis_over_time-%s.csv.gz'%(reproduce_mode)),
                   'phenotype-2': os.path.join(METAIBM_OUT, 'pheno2_dis_over_time-%s.csv.gz'%(reproduce_mode))}

    ''' Nemo parameters '''
    # 2 traits, 20 diploid loci (40 allele copies) per trait; mainland founder phenotypes = (0.5, 0.5), canalized cryptic heterozygous founders (configure_heterozygous_founder):
    # sp1: ('continuous', 20 polygenic het loci +-0.0125): E[G] = 0.50, Vg(F0) = 0.00125 -> Vg_seg(F1) ~ 0.00688
    params1 = nemo_params('sp1')

    ''' Indfactory backended by Nemo: one species = one boot sim = one trait prototype = one Nemo IndFactory '''
    _nemo.cvar.SILENT_RUN = True                # the first boot prints its seed line before it reads run_mode 'silent'; process-global, stays on
    random.seed(seed); np.random.seed(seed)     # MetaIBM draws from both (habitat.py uses Python's random); None: fresh seeds

    sp1_factory = indfactory('sp1', pheno_name_ls=('phenotype-1', 'phenotype-2'), nemo_params=params1, trait_keys=('quant',), seed=seed)
    factories = (sp1_factory,)

    mainland_obj = build_mainland()
    fill_founders(mainland_obj, sp1_factory, reproduce_mode)
    meta_obj = build_metacommunity()
    columns = list(meta_obj.columns_patch_habitat_microsites_id())
    propagules_rain_num = 10

    for time_step in range(all_time_step):
        mainland_obj.meta_dead_selection(base_dead_rate, fitness_wid, method='niche_gaussian')    # the dead are recycled inside
        mainland_born = mainland_birth(mainland_obj, reproduce_mode, asexual_birth_rate, sexual_birth_rate)

        alive = meta_obj.get_meta_individual_num()
        meta_obj.meta_dead_selection(base_dead_rate, fitness_wid, method='niche_gaussian')        # the dead are recycled inside
        dead = alive - meta_obj.get_meta_individual_num()
        born = meta_reproduce(meta_obj, reproduce_mode, asexual_birth_rate, sexual_birth_rate)
        dispersed = disperse_among_islands(meta_obj, disp_among_rate)
        colonized = colonize(meta_obj, mainland_obj, reproduce_mode, propagules_rain_num)
        meta_obj.meta_local_germinate_from_offspring_and_immigrant_pool(is_remove=True)
        ungerminated = meta_obj.meta_clear_up_offspring_and_immigrant_pool()                      # recycles the offspring left without a microsite

        if time_step % record_every == 0:
            save_snapshot(meta_obj, columns, time_step, sp_file, pheno_files)

        if time_step % print_every == 0:
            print('t=%4d  n=%5d  dead=%4d born=%4d dispersed=%3d colonized=%2d ungerminated=%4d  mainland n=%5d born=%4d'% (time_step, meta_obj.get_meta_individual_num(), dead, born, dispersed, colonized, ungerminated, mainland_obj.get_meta_individual_num(), mainland_born), flush=True)

    recycle_everything((meta_obj, mainland_obj), factories)
    return os.path.join(HERE, METAIBM_OUT)


if __name__ == '__main__':
    main(reproduce_mode='asexual', all_time_step=1000, seed=1)
    main(reproduce_mode='sexual', all_time_step=1000, seed=1)
