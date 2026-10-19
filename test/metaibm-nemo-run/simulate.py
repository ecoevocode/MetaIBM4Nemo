#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Full sexual simulations of the experiments/*.py interfaces on the NemoBackend core (metacommunity -> patch -> habitat -> indfactory -> individual).
Nothing in metaibm/, extension/ or experiments/ is edited: each experiment's landscape generator and time loop is copied here with only these interface changes:
  * reproduce_mode='sexual' runs the migrated sexual functions (meta_mainland_sexual_birth_mutate_germinate, meta_sex_reproduce_calculation_with_offspring_marker_pool;
    asex / mix are not migrated); meta_initialize(indfactories, species_2_phenotype_ls, reproduce_mode) with one Nemo indfactory per species
    (traits_num quanti traits, geno_len loci PER TRAIT: Nemo quanti_loci = geno_len x traits_num)
  * every habitat gets a distinct global_index (Nemo homepatch); --steps (default 2000), CSV rows every --save-every steps, invariant checks every --check-every steps
    (one owner, alive, minted - recycled == held, microsite bookkeeping; via test/nemoback-metacommunity/common.py), teardown at the end (minted == recycled)
  * the global-habitat-network dispersal is not tested for now (grfe uses the experiment's commented patch-level alternative)
Sub-commands (outputs in runs/<name>/ unless --out):
    python simulate.py model   [opts]                         experiments/model.py: 1 trait, 2 mainlands, 100 islands x 1 habitat 32x32, environment +0.1 / 100 steps
    python simulate.py sloss   [--patch-num 16] [--hetero] [--schedule marker|pool --dormancy-total N] [opts]     experiments/model-sloss.py
    python simulate.py grfe    [--patch-num 16] [opts]        experiments/model-sloss-GRFE.py (landscape from the experiments' CSV files, read only)
    python simulate.py disturbance [--steps 300 --dist 0.01]  sloss-p16 with a high disturbance rate; every disturbance checked (recycling, untouched patches, invariants)
    python simulate.py all     [--steps 2000 --save-every 10 --max-procs 6 --only NAME ...]   every run of JOBS in parallel (at most --max-procs at once), then RUN_REPORT.md
    python simulate.py report                                  RUN_REPORT.md from runs/*/progress.json (runs may still be going)
"""
import os, sys, ast, json, time, random, argparse, tempfile, hashlib, platform, subprocess, datetime, importlib.util
from collections import defaultdict

sys.dont_write_bytecode = True
os.environ.setdefault('MPLBACKEND', 'Agg')
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
EXPERIMENTS = os.path.join(REPO, 'experiments')
RUNS_DIR = os.path.join(HERE, 'runs')
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
sys.path.insert(0, REPO)

_spec = importlib.util.spec_from_file_location('mc_common', os.path.join(REPO, 'test', 'nemoback-metacommunity', 'common.py'))
mc_common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mc_common)                      # imports metaibm (+ the patch / habitat suites' helpers): owner / leak invariants, teardown

import numpy as np                                       # noqa: E402
import pandas as pd                                      # noqa: E402
import matplotlib.pyplot as plt                          # noqa: E402
from metaibm.indfactory import indfactory                # noqa: E402
from metaibm.patch import patch                          # noqa: E402
from metaibm.metacommunity import metacommunity          # noqa: E402
from pynemo.params import SimParameters                  # noqa: E402



# ====================================================================== run harness: Nemo species, global_index, monitor
NEMO_CACHE = tempfile.mkdtemp(prefix='nemo_run_cache_')


# ------------------------------------------------------------------ species = Nemo indfactories
def nemo_params(name, init_ls, loci=20, mutation_rate=0.0001, mutation_variance=0.01, init_variance=0.002025, env_variance=0.025 ** 2):
    ''' one species' boot parameters: len(init_ls) quanti traits (= traits_num), `loci` loci PER TRAIT (= geno_len), founders start at init_ls;
        environmental variance = pheno_var ** 2 (v3 drew phenotype = mean(genotype) + N(0, pheno_var)).
        Nemo's quanti_loci is the TOTAL over all traits (it warns and rounds when it is not a multiple of quanti_traits), so it is loci * len(init_ls). '''
    fmt = lambda v: '{{%s}}' % ', '.join('%g' % x for x in v)
    rows = [('replicates', 1), ('generations', 1), ('patch_number', 1), ('patch_capacity', 1), ('aging', 1),
            ('run_mode', 'overwrite'), ('root_dir', NEMO_CACHE), ('logfile', NEMO_CACHE + '/nemo.log'), ('filename', name),
            ('quanti_traits', len(init_ls)), ('quanti_loci', loci * len(init_ls)), ('quanti_allele_model', 'continuous'),
            ('quanti_init_model', 2), ('quanti_init_value', fmt(init_ls)), ('quanti_init_variance', fmt([init_variance] * len(init_ls))),
            ('quanti_mutation_rate', mutation_rate), ('quanti_mutation_variance', mutation_variance),
            ('quanti_environmental_variance', fmt([env_variance] * len(init_ls)))]
    p = SimParameters()
    for k, v in rows: p.add(k, v)
    return p


def make_indfactories(species_2_phenotype_ls, pheno_names_ls, seed, loci=20, mutation_rate=0.0001, pheno_var=0.025):
    ''' species 'sp<k+1>' <-> species_2_phenotype_ls[k] (the convention of habitat.hab_initialize) '''
    return {'sp%d' % (k + 1): indfactory('sp%d' % (k + 1), list(pheno_names_ls), nemo_params('sp%d' % (k + 1), ph, loci=loci, mutation_rate=mutation_rate, env_variance=pheno_var ** 2), seed=seed + k)
            for k, ph in enumerate(species_2_phenotype_ls)}


class GlobalIndex():
    ''' hands out distinct habitat.global_index values (Nemo homepatch, 0..65535) across mainland(s) + metacommunity '''
    def __init__(self, start=0): self.next = start
    def __call__(self):
        self.next += 1
        assert self.next - 1 <= 65535
        return self.next - 1


# ------------------------------------------------------------------ run monitor
class Monitor():
    ''' per-phase timing, periodic owner / leak invariants (mc_common.check_world), progress.json, teardown balance '''
    def __init__(self, goal_path, name, metas, factories, all_time_step, check_every=100):
        self.goal_path, self.name, self.metas, self.F, self.T = goal_path, name, metas, factories, all_time_step
        self.fs = list(factories.values())
        self.check_every = check_every
        self.phase = defaultdict(float); self.calls = defaultdict(int)
        self.t0 = time.time(); self.rows = []; self.checks = 0
        self.progress_file = os.path.join(goal_path, 'progress.json')

    def __call__(self, name, fn, *a, **kw):
        t = time.perf_counter()
        r = fn(*a, **kw)
        self.phase[name] += time.perf_counter() - t; self.calls[name] += 1
        return r

    def rss_mb(self):
        return mc_common.rss_mb()

    def step_done(self, t, force=False):
        if not (force or t == 0 or (t + 1) % self.check_every == 0): return
        tc = time.perf_counter()
        mc_common.check_world(self.metas, self.fs, '%s t=%d' % (self.name, t)); self.checks += 1
        self.phase['(invariant checks)'] += time.perf_counter() - tc
        row = {'t': t, 'elapsed_s': round(time.time() - self.t0, 1), 'rss_mb': round(self.rss_mb(), 1),
               'alive': {m.metacommunity_name: m.get_meta_individual_num() for m in self.metas},
               'minted': sum(f.minted_num for f in self.fs), 'recycled': sum(f.recycled_num for f in self.fs),
               'species_alive_metacommunity': self.species_counts(self.metas[-1])}
        self.rows.append(row)
        self.write(status='running')
        print('[%s] t=%d/%d %.0fs rss=%.0fMB alive=%s' % (self.name, t + 1, self.T, row['elapsed_s'], row['rss_mb'], row['alive']), flush=True)

    @staticmethod
    def species_counts(m):
        c = defaultdict(int)
        for p in m.set.values():
            for h in p.set.values():
                for s, d in h.species_category.items(): c[s] += sum(len(v) for v in d.values())
        return dict(sorted(c.items(), key=lambda kv: int(kv[0][2:])))

    def write(self, status, extra=None):
        d = {'name': self.name, 'status': status, 'all_time_step': self.T, 'rows': self.rows, 'phase_seconds': {k: round(v, 1) for k, v in sorted(self.phase.items(), key=lambda kv: -kv[1])},
             'calls': dict(self.calls), 'invariant_checks': self.checks}
        if extra: d.update(extra)
        tmp = self.progress_file + '.tmp'
        with open(tmp, 'w') as fh: json.dump(d, fh, indent=1, default=str)
        os.replace(tmp, self.progress_file)

    def finish(self, extra=None):
        t = time.perf_counter()
        n = mc_common.teardown(self.metas, self.fs)
        bal = {f.species_id: [f.minted_num, f.recycled_num] for f in self.fs}
        for f in self.fs: f.purge_recycling_pool()
        e = dict(extra or {}); e.update(teardown_recycled=n, teardown_balance=bal, teardown_ok=all(a == b for a, b in bal.values()), total_s=round(time.time() - self.t0, 1))
        self.phase['(teardown)'] += time.perf_counter() - t
        self.write(status='done', extra=e)
        return e


def write_logger(log_info, is_logging=False, logger_file=None):
    if is_logging == True: print(log_info, file=logger_file)
    else: print(log_info)


# ====================================================================== experiments/model-sloss.py (sloss)
def new_round(_float, _len=0):
    if isinstance(_float, float):
        if str(_float)[::-1].find('.') <= _len: return (_float)
        if str(_float)[-1] == '5': return (round(float(str(_float)[:-1] + '6'), _len))
        else: return (round(_float, _len))
    else: return (round(_float, _len))


def sloss_generating_empty_metacommunity(meta_name, patch_num, patch_location_ls, hab_num, hab_length, hab_width, dormancy_pool_max_size,
                                   x_axis_environment_values_ls, y_axis_environment_values_ls, environment_types_num, environment_types_name,
                                   environment_variation_ls, is_heterogeneity=False, gidx=None):
    ''' as experiments/model-sloss.py, plus global_index=gidx() for every habitat '''
    log_info = 'generating empty metacommunity ... \n'
    meta_object = metacommunity(metacommunity_name=meta_name)
    patch_num_x_axis, patch_num_y_axis = int(np.sqrt(patch_num)), int(np.sqrt(patch_num))
    hab_num_x_axis, hab_num_y_axis = int(np.sqrt(hab_num)), int(np.sqrt(hab_num))
    for i in range(0, patch_num):
        patch_name = 'patch%d' % (i + 1)
        location = patch_location_ls[i]
        p = patch(patch_name, i, location)
        patch_x_loc, patch_y_loc = location[0], location[1]
        for j in range(hab_num):
            habitat_name = 'h%s' % str(j)
            hab_x_loc, hab_y_loc = patch_x_loc * hab_num_x_axis + j // hab_num_y_axis, patch_y_loc * hab_num_y_axis + j % hab_num_y_axis
            hab_location = (hab_x_loc, hab_y_loc)
            if is_heterogeneity == False:
                x_axis_environment_means_value = x_axis_environment_values_ls[int(hab_x_loc * 4 / (patch_num_x_axis * hab_num_x_axis)) % 4]
                y_axis_environment_means_value = y_axis_environment_values_ls[int(hab_y_loc * 4 / (patch_num_y_axis * hab_num_y_axis)) % 4]
            elif is_heterogeneity == True:
                x_axis_min, x_axis_max = x_axis_environment_values_ls[0], x_axis_environment_values_ls[1]
                y_axis_min, y_axis_max = y_axis_environment_values_ls[0], y_axis_environment_values_ls[1]
                x_step = (x_axis_max - x_axis_min) / (patch_num_x_axis * hab_num_x_axis)
                y_step = (y_axis_max - y_axis_min) / (patch_num_y_axis * hab_num_y_axis)
                x_axis_environment_means_value = x_axis_min + 1 / 2 * x_step + x_step * hab_x_loc
                y_axis_environment_means_value = y_axis_min + 1 / 2 * y_step + y_step * hab_y_loc
            p.add_habitat(hab_name=habitat_name, hab_index=j, hab_location=hab_location, num_env_types=environment_types_num, env_types_name=environment_types_name,
                          mean_env_ls=[x_axis_environment_means_value, y_axis_environment_means_value], var_env_ls=environment_variation_ls, length=hab_length, width=hab_width,
                          dormancy_pool_max_size=dormancy_pool_max_size, global_index=gidx())
            log_info += '%s: %s, %s, %s, %s: x_axis_environment_means_value=%s, y_axis_environment_means_value=%s\n' % (meta_object.metacommunity_name, patch_name, str(location), habitat_name, str(hab_location), str(x_axis_environment_means_value), str(y_axis_environment_means_value))
        meta_object.add_patch(patch_name=patch_name, patch_object=p)
    return meta_object, log_info


def sloss_generate_empty_mainland(meta_name, patch_num, patch_location_ls, hab_num, hab_length, hab_width, dormancy_pool_max_size,
                            micro_environment_values_ls, macro_environment_values_ls, environment_types_num, environment_types_name, environment_variation_ls, gidx=None):
    ''' as experiments/model-sloss.py, plus global_index=gidx() for every habitat '''
    log_info = 'generating empty mainland ... \n'
    meta_object = metacommunity(metacommunity_name=meta_name)
    for i in range(0, patch_num):
        patch_name = 'patch%d' % (i)
        location = patch_location_ls[i]
        patch_x_loc, patch_y_loc = location[0], location[1]
        p = patch(patch_name, i, location)
        hab_num_length, hab_num_width = int(np.sqrt(hab_num)), int(np.sqrt(hab_num))
        for j in range(hab_num):
            hab_location = (patch_x_loc * hab_num_length + j // hab_num_width, patch_y_loc * hab_num_width + j % hab_num_width)
            micro_environment_mean_value = micro_environment_values_ls[j // hab_num_width]
            macro_environment_mean_value = macro_environment_values_ls[j % hab_num_width]
            p.add_habitat(hab_name='h%s' % str(j), hab_index=j, hab_location=hab_location, num_env_types=environment_types_num, env_types_name=environment_types_name,
                          mean_env_ls=[micro_environment_mean_value, macro_environment_mean_value], var_env_ls=environment_variation_ls, length=hab_length, width=hab_width,
                          dormancy_pool_max_size=dormancy_pool_max_size, global_index=gidx())
            log_info += '%s: %s, %s, h%s, %s: micro_environment_mean_value=%s, macro_environment_mean_value=%s \n' % (meta_object.metacommunity_name, patch_name, str(location), j, str(hab_location), str(micro_environment_mean_value), str(macro_environment_mean_value))
        meta_object.add_patch(patch_name=patch_name, patch_object=p)
    return meta_object, log_info


def sloss_main(rep, patch_num, is_heterogeneity, reproduce_mode, total_disp_among_rate, disp_within_rate, patch_dist_rate, goal_path=None, all_time_step=2000, among='patch',
         seed=0, save_every=1, check_every=100, plots=True, name='sloss', schedule='marker', dormancy_total=0):
    if goal_path == None: goal_path = os.getcwd()
    assert reproduce_mode == 'sexual', 'only the sexual branch is migrated'
    np.random.seed(seed); __import__('random').seed(seed)
    all_time_start = time.time()

    meta_length, meta_width = 320, 320
    patch_num_x_axis, patch_num_y_axis = int(np.sqrt(patch_num)), int(np.sqrt(patch_num))
    patch_location_ls = [(i, j) for i in range(patch_num_x_axis) for j in range(patch_num_y_axis)]
    hab_num_in_patch = 16 if patch_num == 1 else 4
    hab_num_x_axis, hab_num_y_axis = int(np.sqrt(hab_num_in_patch)), int(np.sqrt(hab_num_in_patch))
    hab_length, hab_width = int(meta_length / patch_num_y_axis / hab_num_y_axis), int(meta_width / patch_num_x_axis / hab_num_x_axis)

    environment_types_num = 2
    environment_types_name = ('x_axis_environment', 'y_axis_environment')
    environment_variation_ls = [0.025, 0.025]
    total_dormancy_pool_max_size = dormancy_total            # experiment: 0 (only the pool schedule uses the dormancy pool)
    if is_heterogeneity == False: x_axis_environment_values_ls, y_axis_environment_values_ls = [0.2, 0.4, 0.6, 0.8], [0.2, 0.4, 0.6, 0.8]
    elif is_heterogeneity == True: x_axis_environment_values_ls, y_axis_environment_values_ls = [0.1, 0.9], [0.1, 0.9]

    base_dead_rate = 0.1
    fitness_wid = 0.5
    sexual_birth_rate = 1
    mutation_rate = 0.0001

    colonize_rate = 0.001
    propagules_rain_num = 6400 * colonize_rate

    traits_num = 2
    pheno_names_ls = ('x_axis_phenotype', 'y_axis_phenotype')
    pheno_var_ls = (0.025, 0.025)
    geno_len_ls = (20, 20)
    species_2_phenotype_ls = [[i / 10, j / 10] for i in range(2, 10, 2) for j in range(2, 10, 2)]

    is_logging = True
    logger_file = open(goal_path + '/' + 'rep=%d-logger.log' % (rep), 'a')

    write_in_logger_info = ''
    gidx = GlobalIndex(0)
    mainland, log_info = sloss_generate_empty_mainland(meta_name='mainland', patch_num=1, patch_location_ls=[(0, 0)], hab_num=16, hab_length=20, hab_width=20, dormancy_pool_max_size=0,
                                                 micro_environment_values_ls=[0.2, 0.4, 0.6, 0.8], macro_environment_values_ls=[0.2, 0.4, 0.6, 0.8],
                                                 environment_types_num=environment_types_num, environment_types_name=environment_types_name, environment_variation_ls=environment_variation_ls, gidx=gidx)
    write_in_logger_info += log_info
    gidx.next = 100
    meta_obj, log_info = sloss_generating_empty_metacommunity('metacommunity', patch_num, patch_location_ls, hab_num_in_patch, hab_length, hab_width, int(new_round(total_dormancy_pool_max_size / (patch_num * hab_num_in_patch))),
                                                        x_axis_environment_values_ls, y_axis_environment_values_ls, environment_types_num, environment_types_name, environment_variation_ls, is_heterogeneity, gidx=gidx)
    write_in_logger_info += log_info

    indfactories = make_indfactories(species_2_phenotype_ls, pheno_names_ls, seed=1000 * seed + 1, loci=geno_len_ls[0], mutation_rate=mutation_rate, pheno_var=pheno_var_ls[0])
    write_in_logger_info += mainland.meta_initialize(indfactories, species_2_phenotype_ls, reproduce_mode)
    write_logger(write_in_logger_info, is_logging, logger_file)
    mon = Monitor(goal_path, name, [mainland, meta_obj], indfactories, all_time_step, check_every)

    columns_patch_id, columns_habitat_id, columns_mocrosite_id = meta_obj.columns_patch_habitat_microsites_id()
    columns = [columns_patch_id, columns_habitat_id, columns_mocrosite_id]
    files = {k: goal_path + '/' + 'rep=%d-%s_all_time.csv.gz' % (rep, k) for k in ('meta_species_distribution', 'meta_x_axis_phenotype', 'meta_y_axis_phenotype')}
    meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_optimum_sp_id_val(base_dead_rate, fitness_wid, species_2_phenotype_ls), file_name=files['meta_species_distribution'], index=['optimun_sp_id_values'], columns=columns, mode='w')
    meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsite_environment_values(environment_name='x_axis_environment'), file_name=files['meta_x_axis_phenotype'], index=['x_axis_environment_values'], columns=columns, mode='w')
    meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsite_environment_values(environment_name='y_axis_environment'), file_name=files['meta_y_axis_phenotype'], index=['y_axis_environment_values'], columns=columns, mode='w')

    for time_step in range(all_time_step):
        write_in_logger_info = 'time_step=%d \n' % time_step
        ''' dead selection process in mainland and metacommunity '''
        write_in_logger_info += mon('mainland.meta_dead_selection', mainland.meta_dead_selection, base_dead_rate, fitness_wid, method='niche_gaussian')
        write_in_logger_info += mon('mainland.meta_mainland_sexual_birth_mutate_germinate', mainland.meta_mainland_sexual_birth_mutate_germinate, sexual_birth_rate, mutation_rate, pheno_var_ls)
        write_in_logger_info += mon('meta_dead_selection', meta_obj.meta_dead_selection, base_dead_rate, fitness_wid, method='niche_gaussian')
        if schedule == 'pool':                                   # the experiment's commented individual-pool lines, sexual versions
            ''' reproduction process '''
            write_in_logger_info += mon('meta_sex_reproduce_mutate_into_offspring_pool', meta_obj.meta_sex_reproduce_mutate_into_offspring_pool, sexual_birth_rate, mutation_rate, pheno_var_ls)
            ''' dispersal processes '''
            write_in_logger_info += mon('pairwise_sexual_colonization_from_prpagules_rains', meta_obj.pairwise_sexual_colonization_from_prpagules_rains, mainland, propagules_rain_num)
            write_in_logger_info += mon('meta_dispersal_within_patch_from_offspring_to_immigrant_pool', meta_obj.meta_dispersal_within_patch_from_offspring_to_immigrant_pool, disp_within_rate)
            write_in_logger_info += mon('dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool', meta_obj.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool, total_disp_among_rate, method='exponential', rho=2)
            ''' germination processes '''
            write_in_logger_info += mon('meta_local_germinate_from_offspring_immigrant_and_dormancy_pool', meta_obj.meta_local_germinate_from_offspring_immigrant_and_dormancy_pool)
            ''' dormancy process '''
            write_in_logger_info += mon('meta_dormancy_process_from_offspring_pool_and_immigrant_pool', meta_obj.meta_dormancy_process_from_offspring_pool_and_immigrant_pool)
            ''' disturbance process '''
            write_in_logger_info += mon('meta_disturbance_process_in_patches', meta_obj.meta_disturbance_process_in_patches, patch_dist_rate)
            ''' eliminating offspring and immigrant pool '''
            mon('meta_clear_up_offspring_and_immigrant_pool', meta_obj.meta_clear_up_offspring_and_immigrant_pool)
            mon('(logging)', write_logger, write_in_logger_info, is_logging, logger_file)
            if (time_step + 1) % save_every == 0 or time_step == all_time_step - 1:
                mon('(data: sp_id csv.gz)', lambda: meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_sp_id_values(), file_name=files['meta_species_distribution'], index=['time_step%d' % time_step], columns=columns, mode='a'))
                mon('(data: phenotype csv.gz)', lambda: meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_phenotype_values(trait_name='x_axis_phenotype'), file_name=files['meta_x_axis_phenotype'], index=['time_step%d' % time_step], columns=columns, mode='a'))
                mon('(data: phenotype csv.gz)', lambda: meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_phenotype_values(trait_name='y_axis_phenotype'), file_name=files['meta_y_axis_phenotype'], index=['time_step%d' % time_step], columns=columns, mode='a'))
            mon.step_done(time_step)
            continue
        ''' reproduction process '''
        write_in_logger_info += mon('meta_sex_reproduce_calculation_with_offspring_marker_pool', meta_obj.meta_sex_reproduce_calculation_with_offspring_marker_pool, sexual_birth_rate)
        ''' dispersal processes '''
        write_in_logger_info += mon('pairwise_sexual_colonization_from_prpagules_rains', meta_obj.pairwise_sexual_colonization_from_prpagules_rains, mainland, propagules_rain_num)
        write_in_logger_info += mon('meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool', meta_obj.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool, disp_within_rate)
        if among == 'patch':
            write_in_logger_info += mon('dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool', meta_obj.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool, total_disp_among_rate, method='exponential', rho=2)
        else:
            write_in_logger_info += mon('dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool', meta_obj.dispersal_among_patches_in_global_habitat_network_from_offspring_marker_pool_to_immigrant_marker_pool, total_disp_among_rate, method='exponential', rho=2)
        ''' germination processes '''
        write_in_logger_info += mon('meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool', meta_obj.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool, mutation_rate, pheno_var_ls)
        ''' disturbance process '''
        write_in_logger_info += mon('meta_disturbance_process_in_patches', meta_obj.meta_disturbance_process_in_patches, patch_dist_rate)
        ''' eliminating offspring and immigrant (marker) pool '''
        mon('meta_clear_up_offspring_marker_and_immigrant_marker_pool', meta_obj.meta_clear_up_offspring_marker_and_immigrant_marker_pool)

        mon('(logging)', write_logger, write_in_logger_info, is_logging, logger_file)
        if (time_step + 1) % save_every == 0 or time_step == all_time_step - 1:
            mon('(data: sp_id csv.gz)', lambda: meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_sp_id_values(), file_name=files['meta_species_distribution'], index=['time_step%d' % time_step], columns=columns, mode='a'))
            mon('(data: phenotype csv.gz)', lambda: meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_phenotype_values(trait_name='x_axis_phenotype'), file_name=files['meta_x_axis_phenotype'], index=['time_step%d' % time_step], columns=columns, mode='a'))
            mon('(data: phenotype csv.gz)', lambda: meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_phenotype_values(trait_name='y_axis_phenotype'), file_name=files['meta_y_axis_phenotype'], index=['time_step%d' % time_step], columns=columns, mode='a'))
        mon.step_done(time_step)

    if plots:
        kw = dict(sub_row=patch_num_y_axis, sub_col=patch_num_x_axis, hab_num_x_axis_in_patch=hab_num_x_axis, hab_num_y_axis_in_patch=hab_num_y_axis, hab_y_len=hab_length, hab_x_len=hab_width)
        mon('(plots)', meta_obj.meta_show_species_distribution, vmin=1, vmax=16, cmap=plt.get_cmap('tab20'), file_name=goal_path + '/' + 'rep=%d-time_step=%d-metacommunity_sp_dis.jpg' % (rep, time_step), **kw)
        mon('(plots)', meta_obj.meta_show_species_phenotype_distribution, trait_name='x_axis_phenotype', cmap=plt.get_cmap('Blues'), file_name=goal_path + '/' + 'rep=%d-time_step=%d-metacommunity_sp_x_axis_phenotype_dis.jpg' % (rep, time_step), **kw)
        mon('(plots)', meta_obj.meta_show_species_phenotype_distribution, trait_name='y_axis_phenotype', cmap=plt.get_cmap('Greens'), file_name=goal_path + '/' + 'rep=%d-time_step=%d-metacommunity_sp_y_axis_phenotype_dis.jpg' % (rep, time_step), **kw)
        mon('(plots)', meta_obj.meta_show_two_environment_distribution, environment1_name='x_axis_environment', environment2_name='y_axis_environment', mask_loc1='upper', mask_loc2='lower',
            cmap1=plt.get_cmap('Blues'), cmap2=plt.get_cmap('Greens'), file_name=goal_path + '/' + 'rep=%d-metacommunity_environment.jpg' % (rep), **kw)
    write_logger('总模拟运行时间：%.8s s \n' % (time.time() - all_time_start), is_logging, logger_file)
    logger_file.close()
    mon.step_done(all_time_step - 1, force=True)
    return mon.finish({'files_mb': {os.path.basename(f): round(os.path.getsize(f) / 2 ** 20, 1) for f in files.values()}})


# ====================================================================== experiments/model.py (model)
def model_generating_empty_metacommunity(meta_name, patch_num, patch_location_ls, hab_num, hab_length, hab_width, dormancy_pool_max_size,
                                   environment_mean_value, environment_types_num, environment_types_name, environment_variation_ls, gidx=None):
    ''' as experiments/model.py, plus global_index=gidx() for every habitat '''
    log_info = 'generating empty metacommunity ... \n'
    meta_object = metacommunity(metacommunity_name=meta_name)
    hab_num_x_axis, hab_num_y_axis = int(np.sqrt(hab_num)), int(np.sqrt(hab_num))
    for i in range(patch_num):
        patch_name = 'patch%d' % (i + 1)
        location = patch_location_ls[i]
        p = patch(patch_name, i, location)
        patch_x_loc, patch_y_loc = location[0], location[1]
        for j in range(hab_num):
            habitat_name = 'h%s' % str(j + 1)
            hab_location = (patch_x_loc * hab_num_x_axis + j // hab_num_y_axis, patch_y_loc * hab_num_y_axis + j % hab_num_y_axis)
            p.add_habitat(hab_name=habitat_name, hab_index=j, hab_location=hab_location, num_env_types=environment_types_num, env_types_name=environment_types_name,
                          mean_env_ls=[environment_mean_value], var_env_ls=environment_variation_ls, length=hab_length, width=hab_width, dormancy_pool_max_size=dormancy_pool_max_size,
                          global_index=gidx())
            log_info += '%s: %s, %s, %s, %s: environment_mean_value=%s\n' % (meta_object.metacommunity_name, patch_name, str(location), habitat_name, str(hab_location), str(environment_mean_value))
        meta_object.add_patch(patch_name=patch_name, patch_object=p)
    return meta_object, log_info


def model_main(reproduce_mode, patch_dist_rate, mutation_rate, environment_mean_value, delta_mean_ls, rep=0, goal_path=None, all_time_step=2000, seed=0, save_every=1, check_every=100, plots=True, name='model'):
    assert reproduce_mode == 'sexual', 'only the sexual branch is migrated'
    np.random.seed(seed); __import__('random').seed(seed)
    all_time_start = time.time()
    time_step = -1

    patch_num = 100
    patch_num_x_axis, patch_num_y_axis = int(np.sqrt(patch_num)), int(np.sqrt(patch_num))
    patch_location_ls = [(i, j) for i in range(patch_num_x_axis) for j in range(patch_num_y_axis)]
    meta_length, meta_width = 320, 320
    hab_num_in_patch = 1
    hab_num_x_axis, hab_num_y_axis = int(np.sqrt(hab_num_in_patch)), int(np.sqrt(hab_num_in_patch))
    hab_length, hab_width = int(meta_length / patch_num_y_axis / hab_num_y_axis), int(meta_width / patch_num_x_axis / hab_num_x_axis)

    environment_types_num = 1
    environment_types_name = ['environment']
    environment_variation_ls = [0.025]
    start_change_time, end_change_time = time_step + 100, all_time_step - 100
    change_time_step = 100
    env_name_ls = ['environment']
    delta_var_ls = [0]

    species_2_phenotype_ls = [[0.0], [1.0]]
    pheno_names_ls = ['phenotype']
    pheno_var_ls = [0.025]
    geno_len_ls = [20]

    dispersal_amomg_rate = 0.0001
    propagules_rain_num = 10
    base_dead_rate = 0.1
    fitness_wid = 0.5
    sexual_birth_rate = 1

    if goal_path == None: goal_path = os.getcwd()
    is_logging = True
    logger_file = open(goal_path + '/' + 'rep=%d-logger.log' % (rep), 'a')

    write_in_logger_info = ''
    gidx = GlobalIndex(100)
    meta_object, log_info = model_generating_empty_metacommunity('metacommunity', patch_num, patch_location_ls, hab_num_in_patch, hab_length, hab_width, 0,
                                                           environment_mean_value, environment_types_num, environment_types_name, environment_variation_ls, gidx=gidx)
    write_in_logger_info += log_info
    mainland1, log_info = model_generating_empty_metacommunity('mainland1', 1, [(0, 0)], 1, 100, 100, 0, 0.0, 1, ['environment'], [0.025], gidx=GlobalIndex(0))
    write_in_logger_info += log_info
    mainland2, log_info = model_generating_empty_metacommunity('mainland2', 1, [(0, 0)], 1, 100, 100, 0, 1.0, 1, ['environment'], [0.025], gidx=GlobalIndex(1))
    write_in_logger_info += log_info

    indfactories = make_indfactories(species_2_phenotype_ls, pheno_names_ls, seed=1000 * seed + 1, loci=geno_len_ls[0], mutation_rate=mutation_rate, pheno_var=pheno_var_ls[0])
    write_in_logger_info += mainland1.meta_initialize(indfactories, species_2_phenotype_ls, reproduce_mode)
    write_in_logger_info += mainland2.meta_initialize(indfactories, species_2_phenotype_ls, reproduce_mode)
    write_logger(write_in_logger_info, is_logging, logger_file)
    mon = Monitor(goal_path, name, [mainland1, mainland2, meta_object], indfactories, all_time_step, check_every)

    columns_patch_id, columns_habitat_id, columns_mocrosite_id = meta_object.columns_patch_habitat_microsites_id()
    columns = [columns_patch_id, columns_habitat_id, columns_mocrosite_id]
    f_sp, f_ph = goal_path + '/' + 'rep=%d-meta_species_distribution_all_time.csv.gz' % (rep), goal_path + '/' + 'rep=%d-meta_phenotype_all_time.csv.gz' % (rep)
    meta_object.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_object.get_meta_microsites_optimum_sp_id_val(base_dead_rate, fitness_wid, species_2_phenotype_ls), file_name=f_sp, index=['optimun_sp_id_values'], columns=columns, mode='w')
    meta_object.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_object.get_meta_microsite_environment_values(environment_name='environment'), file_name=f_ph, index=['environment'], columns=columns, mode='w')

    env_changes = 0
    for time_step in range(all_time_step):
        write_in_logger_info = 'time_step=%d \n' % time_step
        write_in_logger_info += mon('mainland1.meta_dead_selection', mainland1.meta_dead_selection, base_dead_rate, fitness_wid, method='niche_gaussian')
        write_in_logger_info += mon('mainland1.meta_mainland_sexual_birth_mutate_germinate', mainland1.meta_mainland_sexual_birth_mutate_germinate, sexual_birth_rate, mutation_rate, pheno_var_ls)
        write_in_logger_info += mon('mainland2.meta_dead_selection', mainland2.meta_dead_selection, base_dead_rate, fitness_wid, method='niche_gaussian')
        write_in_logger_info += mon('mainland2.meta_mainland_sexual_birth_mutate_germinate', mainland2.meta_mainland_sexual_birth_mutate_germinate, sexual_birth_rate, mutation_rate, pheno_var_ls)
        write_in_logger_info += mon('meta_dead_selection', meta_object.meta_dead_selection, base_dead_rate, fitness_wid, method='niche_gaussian')
        write_in_logger_info += mon('meta_sex_reproduce_calculation_with_offspring_marker_pool', meta_object.meta_sex_reproduce_calculation_with_offspring_marker_pool, sexual_birth_rate)
        if (time_step + 1) % change_time_step == 0 and start_change_time <= time_step <= end_change_time:
            write_in_logger_info += mon('meta_offset_environmental_values', meta_object.meta_offset_environmental_values, env_name_ls, delta_mean_ls, delta_var_ls); env_changes += 1
        write_in_logger_info += mon('meta_colonize_from_propagules_rains', meta_object.meta_colonize_from_propagules_rains, mainland1, propagules_rain_num)
        write_in_logger_info += mon('meta_colonize_from_propagules_rains', meta_object.meta_colonize_from_propagules_rains, mainland2, propagules_rain_num)
        write_in_logger_info += mon('dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool', meta_object.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool, dispersal_amomg_rate)
        write_in_logger_info += mon('meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool', meta_object.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool, mutation_rate, pheno_var_ls)
        write_in_logger_info += mon('meta_disturbance_process_in_patches', meta_object.meta_disturbance_process_in_patches, patch_dist_rate)
        mon('meta_clear_up_offspring_marker_and_immigrant_marker_pool', meta_object.meta_clear_up_offspring_marker_and_immigrant_marker_pool)
        mon('(logging)', write_logger, write_in_logger_info, is_logging, logger_file)
        if (time_step + 1) % save_every == 0 or time_step == all_time_step - 1:
            mon('(data: sp_id csv.gz)', lambda: meta_object.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_object.get_meta_microsites_individuals_sp_id_values(), file_name=f_sp, index=['time_step%d' % time_step], columns=columns, mode='a'))
            mon('(data: phenotype csv.gz)', lambda: meta_object.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_object.get_meta_microsites_individuals_phenotype_values(trait_name='phenotype'), file_name=f_ph, index=['time_step%d' % time_step], columns=columns, mode='a'))
        mon.step_done(time_step)

    write_logger('all the simulation time：%.8s s \n' % (time.time() - all_time_start), is_logging, logger_file)
    logger_file.close()
    if plots:
        mon('(plots)', meta_object.meta_show_species_distribution, sub_row=10, sub_col=10, hab_num_x_axis_in_patch=1, hab_num_y_axis_in_patch=1, hab_y_len=hab_width, hab_x_len=hab_length, vmin=1, vmax=2, cmap=plt.get_cmap('Paired'), file_name=goal_path + '/' + 'rep=%d-time_step=%d-metacommunity_sp_dis.jpg' % (rep, time_step))
        mon('(plots)', meta_object.meta_show_species_phenotype_distribution, trait_name='phenotype', sub_row=10, sub_col=10, hab_num_x_axis_in_patch=1, hab_num_y_axis_in_patch=1, hab_y_len=hab_width, hab_x_len=hab_length, cmap=plt.get_cmap('Paired'), file_name=goal_path + '/' + 'rep=%d-time_step=%d-metacommunity_sp_x_axis_phenotype_dis.jpg' % (rep, time_step))
    mon.step_done(all_time_step - 1, force=True)
    return mon.finish({'env_changes': env_changes, 'final_env_mean': meta_object.set['patch1'].set['h1'].mean_env_ls, 'files_mb': {os.path.basename(f): round(os.path.getsize(f) / 2 ** 20, 1) for f in (f_sp, f_ph)}})


# ====================================================================== experiments/model-sloss-GRFE.py (grfe)
def _load_patch_habitat_layout(layout_csv_path, patch_num):
    df = pd.read_csv(layout_csv_path)
    return df[df['patch_num'] == patch_num].copy()


def _read_env_multiindex_csv(env_csv):
    df = pd.read_csv(env_csv, header=[0, 1])
    df.columns = pd.MultiIndex.from_tuples(df.columns)
    return df


def _extract_environment_grid(env_df, patch_num, is_same_heterogeneity):
    key = (f'patch_num={patch_num}', str(bool(is_same_heterogeneity)))
    sub = env_df.loc[:, [('meta', 'x_index'), ('meta', 'y_index'), key]].copy()
    sub.columns = ['x_index', 'y_index', 'value']
    sub['x_index'], sub['y_index'] = sub['x_index'].astype(int), sub['y_index'].astype(int)
    return sub


def read_landscape_configuration(layout_csv, env1_csv, env2_csv, patch_num, is_same_heterogeneity):
    df_layout = _load_patch_habitat_layout(layout_csv_path=layout_csv, patch_num=patch_num)
    env1 = _extract_environment_grid(_read_env_multiindex_csv(env1_csv), patch_num, is_same_heterogeneity)
    env2 = _extract_environment_grid(_read_env_multiindex_csv(env2_csv), patch_num, is_same_heterogeneity)
    env = env1.merge(env2, on=['x_index', 'y_index'], how='inner')
    env.columns = ['x_index', 'y_index', 'env_axis1', 'env_axis2']
    return df_layout.merge(env, left_on=['habitat_x_location', 'habitat_y_location'], right_on=['x_index', 'y_index'], how='left')


def grfe_generate_empty_metacommunity(meta_name, layout_csv, env1_csv, env2_csv, patch_num, is_same_heterogeneity, num_env_types, env_types_name, var_env_ls, hab_length, hab_width, dormancy_pool_max_size, gidx=None):
    ''' as experiments/model-sloss-GRFE.py, plus global_index=gidx() for every habitat '''
    log_info = 'generating empty metacommunity ... \n'
    layout = read_landscape_configuration(layout_csv, env1_csv, env2_csv, patch_num, is_same_heterogeneity)
    meta_obj = metacommunity(metacommunity_name=meta_name)
    reg = layout[['patch_id', 'patch_index', 'patch_location_x', 'patch_location_y']].drop_duplicates().sort_values(['patch_index', 'patch_id'])
    for _, pr in reg.iterrows():
        patch_name = pr['patch_id']
        patch_location = (float(pr['patch_location_x']), float(pr['patch_location_y']))
        patch_object = patch(patch_name, int(pr['patch_index']), patch_location)
        for _, hr in layout[layout['patch_id'] == patch_name].copy().sort_values(['habitat_index', 'habitat_id']).iterrows():
            mean_env_ls = [float(hr['env_axis1']), float(hr['env_axis2'])]
            patch_object.add_habitat(hab_name=str(hr['habitat_id']), hab_index=int(hr['habitat_index']), hab_location=(int(hr['habitat_x_location']), int(hr['habitat_y_location'])),
                                     num_env_types=num_env_types, env_types_name=env_types_name, mean_env_ls=mean_env_ls, var_env_ls=var_env_ls,
                                     length=hab_length, width=hab_width, dormancy_pool_max_size=dormancy_pool_max_size, global_index=gidx())
            log_info += '%s: %s, %s, %s: env_axis1=%s, env_axis2=%s\n' % (meta_obj.metacommunity_name, patch_name, str(patch_location), hr['habitat_id'], mean_env_ls[0], mean_env_ls[1])
        meta_obj.add_patch(patch_name=patch_name, patch_object=patch_object)
    return meta_obj, log_info


def grfe_main(rep, patch_num, is_same_heterogeneity, reproduce_mode, total_disp_among_rate, disp_within_rate, patch_dist_rate, goal_path=None, all_time_step=2000, seed=0, save_every=1, check_every=100, plots=True, name='grfe'):
    assert reproduce_mode == 'sexual', 'only the sexual branch is migrated'
    if goal_path == None: goal_path = os.getcwd()
    np.random.seed(seed); __import__('random').seed(seed)
    all_time_start = time.time()
    patch_num_x_axis, patch_num_y_axis = int(np.sqrt(patch_num)), int(np.sqrt(patch_num))
    hab_num_in_patch = int(1024 / patch_num)
    hab_num_x_axis, hab_num_y_axis = int(np.sqrt(hab_num_in_patch)), int(np.sqrt(hab_num_in_patch))
    hab_length, hab_width = 10, 10
    environment_types_num = 2
    environment_types_name = ('env_axis1', 'env_axis2')
    environment_variation_ls = [0.025, 0.025]
    base_dead_rate, fitness_wid, sexual_birth_rate, mutation_rate = 0.1, 0.5, 1, 0.0001
    propagules_rain_num = 6400 * 0.001
    pheno_names_ls = ('phenotype_axis1', 'phenotype_axis2')
    pheno_var_ls, geno_len_ls = (0.025, 0.025), (20, 20)
    species_2_phenotype_ls = [[i / 10, j / 10] for i in range(2, 10, 2) for j in range(2, 10, 2)]

    is_logging = True
    logger_file = open(goal_path + '/' + 'rep=%d-logger.log' % (rep), 'a')
    write_in_logger_info = ''
    gidx = GlobalIndex(0)
    mainland, log_info = sloss_generate_empty_mainland(meta_name='mainland', patch_num=1, patch_location_ls=[(0, 0)], hab_num=16, hab_length=20, hab_width=20, dormancy_pool_max_size=0,
                                                           micro_environment_values_ls=[0.2, 0.4, 0.6, 0.8], macro_environment_values_ls=[0.2, 0.4, 0.6, 0.8],
                                                           environment_types_num=environment_types_num, environment_types_name=environment_types_name, environment_variation_ls=environment_variation_ls, gidx=gidx)
    write_in_logger_info += log_info
    gidx.next = 100
    t_build = time.time()
    meta_obj, log_info = grfe_generate_empty_metacommunity('metacommunity', os.path.join(EXPERIMENTS, 'patch_habitat_layouts.csv'), os.path.join(EXPERIMENTS, '32x32_habitats_env1.csv'),
                                                      os.path.join(EXPERIMENTS, '32x32_habitats_env2.csv'), patch_num, is_same_heterogeneity, environment_types_num, environment_types_name,
                                                      environment_variation_ls, hab_length, hab_width, 0, gidx=gidx)
    build_s = time.time() - t_build
    write_in_logger_info += log_info
    indfactories = make_indfactories(species_2_phenotype_ls, pheno_names_ls, seed=1000 * seed + 1, loci=geno_len_ls[0], mutation_rate=mutation_rate, pheno_var=pheno_var_ls[0])
    write_in_logger_info += mainland.meta_initialize(indfactories, species_2_phenotype_ls, reproduce_mode)
    write_logger(write_in_logger_info, is_logging, logger_file)
    mon = Monitor(goal_path, name, [mainland, meta_obj], indfactories, all_time_step, check_every)
    mon.phase['(build metacommunity incl. global habitat distance matrix)'] = build_s

    columns = list(meta_obj.columns_patch_habitat_microsites_id())
    files = {k: goal_path + '/' + 'rep=%d-%s_all_time.csv.gz' % (rep, k) for k in ('meta_species_distribution', 'meta_x_axis_phenotype', 'meta_y_axis_phenotype')}
    meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_optimum_sp_id_val(base_dead_rate, fitness_wid, species_2_phenotype_ls), file_name=files['meta_species_distribution'], index=['optimun_sp_id_values'], columns=columns, mode='w')
    meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsite_environment_values(environment_name='env_axis1'), file_name=files['meta_x_axis_phenotype'], index=['x_axis_environment_values'], columns=columns, mode='w')
    meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsite_environment_values(environment_name='env_axis2'), file_name=files['meta_y_axis_phenotype'], index=['y_axis_environment_values'], columns=columns, mode='w')

    for time_step in range(all_time_step):
        write_in_logger_info = 'time_step=%d \n' % time_step
        write_in_logger_info += mon('mainland.meta_dead_selection', mainland.meta_dead_selection, base_dead_rate, fitness_wid, method='niche_gaussian')
        write_in_logger_info += mon('mainland.meta_mainland_sexual_birth_mutate_germinate', mainland.meta_mainland_sexual_birth_mutate_germinate, sexual_birth_rate, mutation_rate, pheno_var_ls)
        write_in_logger_info += mon('meta_dead_selection', meta_obj.meta_dead_selection, base_dead_rate, fitness_wid, method='niche_gaussian')
        write_in_logger_info += mon('meta_sex_reproduce_calculation_with_offspring_marker_pool', meta_obj.meta_sex_reproduce_calculation_with_offspring_marker_pool, sexual_birth_rate)
        write_in_logger_info += mon('pairwise_sexual_colonization_from_prpagules_rains', meta_obj.pairwise_sexual_colonization_from_prpagules_rains, mainland, propagules_rain_num)
        write_in_logger_info += mon('meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool', meta_obj.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool, disp_within_rate)
        write_in_logger_info += mon('dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool', meta_obj.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool, total_disp_among_rate, method='exponential', rho=0.2)   # the experiment's commented alternative: global-habitat-network not tested for now
        write_in_logger_info += mon('meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool', meta_obj.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool, mutation_rate, pheno_var_ls)
        write_in_logger_info += mon('meta_disturbance_process_in_patches', meta_obj.meta_disturbance_process_in_patches, patch_dist_rate)
        mon('meta_clear_up_offspring_marker_and_immigrant_marker_pool', meta_obj.meta_clear_up_offspring_marker_and_immigrant_marker_pool)
        mon('(logging)', write_logger, write_in_logger_info, is_logging, logger_file)
        if (time_step + 1) % save_every == 0 or time_step == all_time_step - 1:
            mon('(data: sp_id csv.gz)', lambda: meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_sp_id_values(), file_name=files['meta_species_distribution'], index=['time_step%d' % time_step], columns=columns, mode='a'))
            mon('(data: phenotype csv.gz)', lambda: meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_phenotype_values(trait_name='phenotype_axis1'), file_name=files['meta_x_axis_phenotype'], index=['time_step%d' % time_step], columns=columns, mode='a'))
            mon('(data: phenotype csv.gz)', lambda: meta_obj.meta_distribution_data_all_time_to_csv_gz(dis_data_all_time=meta_obj.get_meta_microsites_individuals_phenotype_values(trait_name='phenotype_axis2'), file_name=files['meta_y_axis_phenotype'], index=['time_step%d' % time_step], columns=columns, mode='a'))
        mon.step_done(time_step)

    if plots:
        kw = dict(sub_row=patch_num_y_axis, sub_col=patch_num_x_axis, hab_num_x_axis_in_patch=hab_num_x_axis, hab_num_y_axis_in_patch=hab_num_y_axis, hab_y_len=hab_length, hab_x_len=hab_width)
        mon('(plots)', meta_obj.meta_show_species_distribution, vmin=1, vmax=16, cmap=plt.get_cmap('tab20'), file_name=goal_path + '/' + 'rep=%d-time_step=%d-metacommunity_sp_dis.jpg' % (rep, time_step), **kw)
        mon('(plots)', meta_obj.meta_show_species_phenotype_distribution, trait_name='phenotype_axis1', cmap=plt.get_cmap('turbo'), file_name=goal_path + '/' + 'rep=%d-time_step=%d-metacommunity_sp_x_axis_phenotype_dis.jpg' % (rep, time_step), **kw)
        mon('(plots)', meta_obj.meta_show_species_phenotype_distribution, trait_name='phenotype_axis2', cmap=plt.get_cmap('turbo'), file_name=goal_path + '/' + 'rep=%d-time_step=%d-metacommunity_sp_y_axis_phenotype_dis.jpg' % (rep, time_step), **kw)
        mon('(plots)', meta_obj.meta_show_environment_distribution, environment_name='env_axis1', mask_loc=None, cmap=plt.get_cmap('turbo'), file_name=goal_path + '/' + 'rep=%d-metacommunity_axis1_environment.jpg' % (rep), **kw)
        mon('(plots)', meta_obj.meta_show_environment_distribution, environment_name='env_axis2', mask_loc=None, cmap=plt.get_cmap('turbo'), file_name=goal_path + '/' + 'rep=%d-metacommunity_axis2_environment.jpg' % (rep), **kw)
    write_logger('总模拟运行时间：%.8s s \n' % (time.time() - all_time_start), is_logging, logger_file)
    logger_file.close()
    mon.step_done(all_time_step - 1, force=True)
    return mon.finish({'files_mb': {os.path.basename(f): round(os.path.getsize(f) / 2 ** 20, 1) for f in files.values()}})


# ====================================================================== high-disturbance short run (disturbance)
def disturbance_snapshot(p):
    return {'sites': {(h.name, i, j): id(h.set['microsite_individuals'][i][j]) for h in p.set.values() for (i, j) in h.occupied_site_pos_ls},
            'pools': {(h.name, q): [id(x) for x in getattr(h, q)] for h in p.set.values() for q in mc_common.POOLS},
            'markers': {(h.name, q): list(getattr(h, q)) for h in p.set.values() for q in mc_common.MARKER_POOLS}}


def disturbance_main(out, steps, dist, seed):
    np.random.seed(seed); random.seed(seed)
    gidx = GlobalIndex(0)
    env = dict(environment_types_num=2, environment_types_name=('x_axis_environment', 'y_axis_environment'), environment_variation_ls=[0.025, 0.025])
    mainland, _ = sloss_generate_empty_mainland('mainland', 1, [(0, 0)], 16, 20, 20, 0, [0.2, 0.4, 0.6, 0.8], [0.2, 0.4, 0.6, 0.8], gidx=gidx, **env)
    gidx.next = 100
    meta, _ = sloss_generating_empty_metacommunity('metacommunity', 16, [(i, j) for i in range(4) for j in range(4)], 4, 40, 40, 0, [0.2, 0.4, 0.6, 0.8], [0.2, 0.4, 0.6, 0.8],
                                                       env['environment_types_num'], env['environment_types_name'], env['environment_variation_ls'], False, gidx=gidx)
    species_2_phenotype_ls = [[i / 10, j / 10] for i in range(2, 10, 2) for j in range(2, 10, 2)]
    F = make_indfactories(species_2_phenotype_ls, ('x_axis_phenotype', 'y_axis_phenotype'), seed=1000 * seed + 1)
    fs = list(F.values())
    mainland.meta_initialize(F, species_2_phenotype_ls, 'sexual')
    world = [mainland, meta]
    events, recol, t0 = [], {}, time.time()
    for t in range(steps):
        mainland.meta_dead_selection(0.1, 0.5, method='niche_gaussian')
        mainland.meta_mainland_sexual_birth_mutate_germinate(1, 0.0001, (0.025, 0.025))
        meta.meta_dead_selection(0.1, 0.5, method='niche_gaussian')
        meta.meta_sex_reproduce_calculation_with_offspring_marker_pool(1)
        meta.pairwise_sexual_colonization_from_prpagules_rains(mainland, 6.4)
        meta.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(0.1)
        meta.dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool(0.1, method='exponential', rho=2)
        meta.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(0.0001, (0.025, 0.025))
        # ---- instrumented disturbance
        before = {pn: disturbance_snapshot(p) for pn, p in meta.set.items()}
        objs = {pn: [x for h in p.set.values() for x in mc_common.micro(h) + h.offspring_pool + h.immigrant_pool] for pn, p in meta.set.items()}
        dorm = {pn: {h.name: list(h.dormancy_pool) for h in p.set.values()} for pn, p in meta.set.items()}
        m0, r0 = {f.species_id: f.minted_num for f in fs}, {f.species_id: f.recycled_num for f in fs}
        log = meta.meta_disturbance_process_in_patches(dist)
        hit = ast.literal_eval(log[log.index('occurred in ') + len('occurred in '):].strip())
        mc_common.check_world(world, fs, 'after disturbance t=%d' % t)
        exp_rec = {f.species_id: 0 for f in fs}
        for pn in hit:
            for x in objs[pn]: exp_rec[x.species_id] += 1
        assert {f.species_id: f.recycled_num - r0[f.species_id] for f in fs} == exp_rec, 't=%d recycled per factory != individuals of the disturbed patches' % t
        assert all(f.minted_num == m0[f.species_id] for f in fs), 't=%d disturbance minted'
        for pn, p in meta.set.items():
            if pn in hit:
                assert all(not x.is_alive and x.nemo_indi_obj is None for x in objs[pn]), 't=%d %s: a disturbed individual is still alive' % (t, pn)
                assert all(h.indi_num == 0 and not h.offspring_pool and not h.immigrant_pool and not h.offspring_marker_pool and not h.immigrant_marker_pool for h in p.set.values()), pn
                assert all([id(x) for x in h.dormancy_pool] == [id(x) for x in dorm[pn][h.name]] for h in p.set.values()), 'dormancy pool touched'
            else:
                assert disturbance_snapshot(p) == before[pn], 't=%d %s: undisturbed patch changed' % (t, pn)
        emptied = sorted(pn for pn, p in meta.set.items() if before[pn]['sites'] and not any(h.indi_num for h in p.set.values()))
        assert emptied == sorted(pn for pn in hit if before[pn]['sites']), (t, emptied, hit)
        for pn in hit:
            events.append({'t': t, 'patch': pn, 'individuals_recycled': len(objs[pn]), 'markers_cleared': sum(len(v) for v in before[pn]['markers'].values())})
            recol[(t, pn)] = []
        meta.meta_clear_up_offspring_marker_and_immigrant_marker_pool()
        for (te, pn), ls in recol.items():
            if t - te in (1, 5, 10, 25, 50, 100): ls.append((t - te, sum(h.indi_num for h in meta.set[pn].set.values())))
        if (t + 1) % 50 == 0:
            print('[disturbance] t=%d/%d %.0fs events=%d islands=%d rss=%.0fMB' % (t + 1, steps, time.time() - t0, len(events), meta.get_meta_individual_num(), mc_common.rss_mb()), flush=True)
    n = mc_common.teardown(world, fs)
    res = {'steps': steps, 'patch_dist_rate': dist, 'expected_events': 16 * steps * dist, 'events': len(events), 'individuals_recycled_by_disturbance': sum(e['individuals_recycled'] for e in events),
           'markers_cleared': sum(e['markers_cleared'] for e in events), 'per_event': events, 'recolonization (steps after, individuals)': {'%s@%d' % (pn, te): ls for (te, pn), ls in recol.items()},
           'teardown_recycled': n, 'teardown_balance': {f.species_id: [f.minted_num, f.recycled_num] for f in fs}, 'seconds': round(time.time() - t0, 1)}
    json.dump(res, open(os.path.join(out, 'disturbance.json'), 'w'), indent=1)
    print('DONE events=%d (expected %.0f), recycled by disturbance %d, markers cleared %d, teardown minted == recycled: %s' % (
        len(events), res['expected_events'], res['individuals_recycled_by_disturbance'], res['markers_cleared'], all(a == b for a, b in res['teardown_balance'].values())))


# ====================================================================== all: parallel launcher + RUN_REPORT.md
def git(*a): return subprocess.run(['git'] + list(a), capture_output=True, text=True, cwd=REPO).stdout.rstrip()
def md5(rel): return hashlib.md5(open(os.path.join(REPO, rel), 'rb').read()).hexdigest()
def esc(s): return str(s if s is not None else '').replace('|', '/').replace('\n', ' ')
def table(head, rows): return ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)] + ['| ' + ' | '.join(esc(c) for c in r) + ' |' for r in rows]

# name, sub-command + args, interface (the runs of `simulate.py all`); 'disturbance' takes its own --steps / --dist
JOBS = [
    ('model', ['model'], 'experiments/model.py: 1 trait, 2 mainlands (100x100), 100 islands x 1 habitat 32x32, two colonizations / step, environment +0.1 every 100 steps'),
    ('sloss-p16', ['sloss', '--patch-num', '16'], 'experiments/model-sloss.py: 16 species, mainland 16 x 20x20, 16 patches x 4 habitats 40x40 (discrete environment), pairwise colonization'),
    ('sloss-p64-hetero', ['sloss', '--patch-num', '64', '--hetero'], 'experiments/model-sloss.py: 64 patches x 4 habitats 20x20, continuous environment (is_heterogeneity=True)'),
    ('sloss-p1', ['sloss', '--patch-num', '1'], 'experiments/model-sloss.py: 1 patch x 16 habitats 80x80 (among-patch dispersal takes the patch_num < 2 branch; slow: marker birth re-pairs the 6400 residents per birth)'),
    ('sloss-pool-p16', ['sloss', '--patch-num', '16', '--schedule', 'pool', '--dormancy-total', '10240'], 'experiments/model-sloss.py with its commented individual-pool lines (sexual): offspring pool -> pairwise colonization -> within / among dispersal to immigrant pools -> germinate O+I+D -> dormancy (10240 total, 160 / habitat) -> disturbance -> clear-ups'),
    ('grfe-p16', ['grfe', '--patch-num', '16'], 'experiments/model-sloss-GRFE.py: landscape from the CSV layouts, 16 patches x 64 habitats 10x10, among-patch marker dispersal (rho=0.2)'),
    ('disturbance-p16', ['disturbance'], 'sloss-p16 landscape, patch_dist_rate 0.01, 300 steps, every disturbance checked')]
RUNS = [(n, a, i) for n, a, i in JOBS if a[0] != 'disturbance']
FILES = ['metaibm/metacommunity.py', 'metaibm/patch.py', 'metaibm/habitat.py', 'metaibm/indfactory.py', 'metaibm/individual.py', 'extension/global_habitat_network.py']


def launch(names, steps, save_every, check_every, max_procs):
    """ start every selected job (at most max_procs at once), wait, write runs/returncodes.json """
    os.makedirs(RUNS_DIR, exist_ok=True)
    todo, procs, rc = [j for j in JOBS if not names or j[0] in names], {}, {}
    while todo or len(rc) < len(procs):
        while todo and sum(1 for n in procs if n not in rc) < max_procs:
            name, args, _ = todo.pop(0)
            out = os.path.join(RUNS_DIR, name); os.makedirs(out, exist_ok=True)
            for f in os.listdir(out): os.remove(os.path.join(out, f))
            extra = ['--out', out] + (['--steps', str(steps), '--save-every', str(save_every), '--check-every', str(check_every), '--name', name] if args[0] != 'disturbance' else [])
            log = open(os.path.join(out, 'stdout.log'), 'w')
            procs[name] = (subprocess.Popen([sys.executable, os.path.abspath(__file__)] + args + extra, stdout=log, stderr=subprocess.STDOUT, cwd=HERE, env=ENV), log, time.time())
            print('started %-18s pid %d' % (name, procs[name][0].pid), flush=True)
        for n, (p, log, t0) in procs.items():
            if n not in rc and p.poll() is not None:
                rc[n] = (p.returncode, round(time.time() - t0, 1)); log.close()
                print('finished %-18s rc=%d after %.0f s' % (n, p.returncode, rc[n][1]), flush=True)
        time.sleep(5)
    json.dump(rc, open(os.path.join(RUNS_DIR, 'returncodes.json'), 'w'), indent=1)
    return rc

def report():
    rc = json.load(open(os.path.join(RUNS_DIR, 'returncodes.json'))) if os.path.exists(os.path.join(RUNS_DIR, 'returncodes.json')) else {}
    L = ['# MetaIBM-Nemo full simulations (sexual, experiments/*.py interfaces)', '', 'Generated %s by `test/metaibm-nemo-run/simulate.py report`.' % datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), '',
         '- python `%s`, numpy `%s`, `%s`; HEAD `%s` (`%s`), branch `%s`' % (sys.version.split()[0], __import__('numpy').__version__, platform.platform(), git('rev-parse', '--short', 'HEAD'), git('log', '-1', '--format=%s'), git('rev-parse', '--abbrev-ref', 'HEAD')), '']
    L += table(('file', 'md5', 'state'), [('`%s`' % f, '`%s`' % md5(f), 'modified (working tree)' if git('status', '--short', '--', f) else 'committed') for f in FILES])
    L += ['', 'Nothing in metaibm/, extension/ or experiments/ was edited: simulate.py copies each experiment\'s landscape generator and time loop, with only these interface changes: '
          'sexual functions instead of the (not migrated) mix / asex ones, `meta_initialize(indfactories, species_2_phenotype_ls, reproduce_mode)` with one Nemo indfactory per species, '
          'a distinct `global_index` per habitat, CSV rows every `--save-every` steps (the experiments write every step: ~1 GB per run at 102400 sites), invariant checks every 100 steps, teardown at the end. '
          'The global-habitat-network dispersal is not tested for now (GRFE uses the experiment\'s commented patch-level alternative).', '']
    rows, details = [], []
    for name, extra, iface in RUNS:
        pf = os.path.join(RUNS_DIR, name, 'progress.json')
        if not os.path.exists(pf): rows.append((name, 'simulate.py ' + ' '.join(extra), 'not started', '', '', '', '', '', '')); continue
        d = json.load(open(pf)); r = d['rows'][-1] if d['rows'] else {}
        st = d['status'] if name not in rc else ('done' if rc[name][0] == 0 and d['status'] == 'done' else 'stopped by user (SIGTERM, rc=-15)' if rc[name][0] == -15 else 'FAILED rc=%s' % rc[name][0])
        if d['status'].startswith('stopped'): st = 'stopped by user (SIGTERM, rc=-15)'
        steps_done = r.get('t', -1) + 1
        sim_rows = d['rows'][:-1] if d['status'] == 'done' and len(d['rows']) > 1 else d['rows']      # the last row is taken after the final plots
        rss = [x['rss_mb'] for x in sim_rows]
        sps = [x for x in d['rows'] if x['t'] >= 99]
        sec_step = (sps[-1]['elapsed_s'] - sps[0]['elapsed_s']) / max(1, sps[-1]['t'] - sps[0]['t']) if len(sps) > 1 else None
        rows.append((name, 'simulate.py ' + ' '.join(extra), st, '%d / %d' % (steps_done, d['all_time_step']), '%.2f' % sec_step if sec_step else '', '%s' % r.get('alive'),
                     '%d' % len([v for v in r.get('species_alive_metacommunity', {}).values() if v]), '%s -> %s (max %s)%s' % (rss[0], rss[-1], max(rss), '; after plots %s' % d['rows'][-1]['rss_mb'] if sim_rows is not d['rows'] else '') if rss else '',
                     ('OK %s' % (sum(a for a, b in d['teardown_balance'].values()))) if d.get('teardown_ok') else ('FAIL' if 'teardown_ok' in d else '-')))
        tail = open(os.path.join(RUNS_DIR, name, 'stdout.log')).read().splitlines()[-25:] if os.path.exists(os.path.join(RUNS_DIR, name, 'stdout.log')) else []
        details += ['', '### %s' % name, '', '- interface: %s' % iface, '- invariant checks (one owner, alive, minted - recycled == held, microsite bookkeeping): %d, all passed%s' % (d.get('invariant_checks', 0), '' if 'FAIL' not in st else ' up to the failure'),
                    '- files: %s; total %s s' % (d.get('files_mb'), d.get('total_s')), '- phase seconds: %s' % ', '.join('%s %s' % kv for kv in list(d['phase_seconds'].items())[:10])]
        if d.get('env_changes') is not None: details.append('- environment changes: %s, final island mean %s' % (d['env_changes'], d.get('final_env_mean')))
        details += ['', '| t | elapsed s | RSS MB | alive | minted | recycled |', '|---|---|---|---|---|---|'] + ['| %d | %s | %s | %s | %d | %d |' % (x['t'] + 1, x['elapsed_s'], x['rss_mb'], esc(x['alive']), x['minted'], x['recycled']) for x in d['rows'][::2] + d['rows'][-1:]]
        if 'FAIL' in st: details += ['', '```'] + tail + ['```']
    L += ['## Runs', ''] + table(('run', 'command', 'status', 'steps', 's / step (after t=100)', 'alive at last check', 'species on islands', 'RSS MB first -> last simulation check (max)', 'teardown minted == recycled (individuals)'), rows)
    dj = os.path.join(RUNS_DIR, 'disturbance-p16', 'disturbance.json')
    if os.path.exists(dj):
        z = json.load(open(dj)); rc_ = z['recolonization (steps after, individuals)']
        mean_at = lambda lag: round(sum(x[1] for ls in rc_.values() for x in ls if x[0] == lag) / max(1, sum(1 for ls in rc_.values() for x in ls if x[0] == lag)))
        L += ['', '## High-disturbance short run (`simulate.py disturbance`: sloss-p16 landscape, patch_dist_rate %s, %d steps)' % (z['patch_dist_rate'], z['steps']), '',
              '- %d disturbances (expected %.0f), every one checked: the disturbed patches\' individuals (microsites + offspring / immigrant pools) all recycled with the exact per-factory count, their marker pools emptied, dormancy kept, every other patch untouched, full invariants after the call' % (z['events'], z['expected_events']),
              '- recycled by disturbance %d individuals (mean %.0f per event), %d markers cleared; teardown minted == recycled: %s; %s s' % (z['individuals_recycled_by_disturbance'], z['individuals_recycled_by_disturbance'] / max(1, z['events']), z['markers_cleared'], all(a == b for a, b in z['teardown_balance'].values()), z['seconds']),
              '- recolonization of a disturbed patch (6400 sites), mean individuals after 1 / 5 / 10 / 25 / 50 / 100 steps: %s' % ' / '.join(str(mean_at(l)) for l in (1, 5, 10, 25, 50, 100))]
    L += ['', '## Details'] + details
    L += ['', '## Repository state', '', '```', git('status', '--short'), '```']
    open(os.path.join(HERE, 'RUN_REPORT.md'), 'w').write('\n'.join(L) + '\n')
    print('report:', os.path.join(HERE, 'RUN_REPORT.md'))


# ====================================================================== command line
def _sim_opts(ap, name):
    ap.add_argument('--out', default=None); ap.add_argument('--steps', type=int, default=2000); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--save-every', type=int, default=1); ap.add_argument('--check-every', type=int, default=100); ap.add_argument('--no-plots', action='store_true'); ap.add_argument('--name', default=name)
    ap.add_argument('--dist', type=float, default=0.00001)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    m = sub.add_parser('model'); _sim_opts(m, 'model'); m.add_argument('--mutation', type=float, default=0.0001); m.add_argument('--env', type=float, default=0.0); m.add_argument('--delta', type=float, default=0.1)
    s = sub.add_parser('sloss'); _sim_opts(s, 'sloss'); s.add_argument('--patch-num', type=int, default=16); s.add_argument('--hetero', action='store_true'); s.add_argument('--among', default='patch')
    s.add_argument('--disp-among', type=float, default=0.1); s.add_argument('--disp-within', type=float, default=0.1)
    s.add_argument('--schedule', choices=('marker', 'pool'), default='marker'); s.add_argument('--dormancy-total', type=int, default=0)
    g = sub.add_parser('grfe'); _sim_opts(g, 'grfe'); g.add_argument('--patch-num', type=int, default=16); g.add_argument('--not-same-heterogeneity', action='store_true')
    g.add_argument('--disp-among', type=float, default=0.01); g.add_argument('--disp-within', type=float, default=0.1)
    d = sub.add_parser('disturbance'); d.add_argument('--out', default=None); d.add_argument('--steps', type=int, default=300); d.add_argument('--dist', type=float, default=0.01); d.add_argument('--seed', type=int, default=0)
    a_ = sub.add_parser('all'); a_.add_argument('--steps', type=int, default=2000); a_.add_argument('--save-every', type=int, default=10); a_.add_argument('--check-every', type=int, default=100)
    a_.add_argument('--max-procs', type=int, default=6); a_.add_argument('--only', nargs='*', default=[])
    sub.add_parser('report')
    a = ap.parse_args()
    if a.cmd in ('model', 'sloss', 'grfe', 'disturbance'):
        a.out = a.out or os.path.join(RUNS_DIR, getattr(a, 'name', 'disturbance-p16')); os.makedirs(a.out, exist_ok=True)
    if a.cmd == 'model':
        r = model_main('sexual', a.dist, a.mutation, a.env, [a.delta], rep=0, goal_path=a.out, all_time_step=a.steps, seed=a.seed, save_every=a.save_every, check_every=a.check_every, plots=not a.no_plots, name=a.name)
    elif a.cmd == 'sloss':
        r = sloss_main(rep=0, patch_num=a.patch_num, is_heterogeneity=a.hetero, reproduce_mode='sexual', total_disp_among_rate=a.disp_among, disp_within_rate=a.disp_within, patch_dist_rate=a.dist,
                       goal_path=a.out, all_time_step=a.steps, among=a.among, seed=a.seed, save_every=a.save_every, check_every=a.check_every, plots=not a.no_plots, name=a.name,
                       schedule=a.schedule, dormancy_total=a.dormancy_total)
    elif a.cmd == 'grfe':
        r = grfe_main(0, a.patch_num, not a.not_same_heterogeneity, 'sexual', a.disp_among, a.disp_within, a.dist, goal_path=a.out, all_time_step=a.steps, seed=a.seed,
                      save_every=a.save_every, check_every=a.check_every, plots=not a.no_plots, name=a.name)
    elif a.cmd == 'disturbance':
        r = disturbance_main(a.out, a.steps, a.dist, a.seed)
    elif a.cmd == 'all':
        t = time.time(); rc = launch(a.only, a.steps, a.save_every, a.check_every, a.max_procs)
        print('all finished in %.0f s: %s' % (time.time() - t, rc)); report(); return
    else:
        report(); return
    print('DONE', r)


if __name__ == '__main__':
    main()
