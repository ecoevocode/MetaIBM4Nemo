# -*- coding: utf-8 -*-
"""
Split from metacommunity_IBM.py.
Compatibility-preserving class module.
"""

import numpy as np
import random
import re

from .habitat import habitat


class patch():
    def __init__(self, patch_name, patch_index, location):
        self.name = patch_name
        self.index = patch_index
        self.set: dict[str, habitat] = {}            # {hab_name: habitat}
        self.hab_num = 0
        self.hab_id_ls = []
        self.location = location
        
    def get_data(self):
        output = {}
        for key, value in self.set.items():
            output[key]=value.set
        return output

    def add_habitat(self, hab_name, hab_index, hab_location, num_env_types, env_types_name, mean_env_ls, var_env_ls, length, width, dormancy_pool_max_size=0, global_index=0):
        ''' create a habitat and add it to the patch
        INPUT: (as habitat.__init__); global_index = int 0-65535 -- unique habitat index across the metacommunity, used as Nemo's homepatch (default 0) '''
        h_object = habitat(hab_name, hab_index, hab_location, num_env_types, env_types_name, mean_env_ls, var_env_ls, length, width, dormancy_pool_max_size, global_index=global_index)
        self.set[hab_name] = h_object
        self.hab_id_ls.append(hab_name)
        self.hab_num += 1
        
    def patch_initialize(self, indfactories, species_2_phenotype_ls, reproduce_mode='sexual'):
        ''' NemoBackend: fill every microsite of every habitat in the patch with founders of the species matched to each habitat's mean environment (if occupied, residents recycled first, pools untouched)
        INPUT: indfactories = {species_id: indfactory}; species_2_phenotype_ls = species' founder phenotypes, species 'sp<k+1>' <-> index k; reproduce_mode = 'sexual' (default) | 'asexual' (all founders female) '''
        for h_id, h_object in self.set.items():
            h_object.hab_initialize(indfactories, species_2_phenotype_ls, reproduce_mode)   # NemoBackend Methods: indfactory.make_founder() <--> IndFactory::makeNewIndividual() + Individual::create_first_gen()

    def patch_dead_selection(self, base_dead_rate, fitness_wid, method):
        ''' NemoBackend: dead selection in every habitat of the patch -- each individual survives with survival_rate(phenotype, micro-environment); the dead are recycled
        INPUT: base_dead_rate (d); fitness_wid (w); method = 'niche_gaussian' | 'neutral_uniform'
        RETURN: number dead '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_dead_selection(base_dead_rate, fitness_wid, method)     # NemoBackend Methods: recycling the dead individuals
        return counter

    def get_patch_individual_num(self):
        num = 0
        for key, value in self.set.items():
            num += value.indi_num
        return num

    def get_patch_empty_sites_ls(self):
        ''' return patch_empty_pos_ls as [(h_id, len_id, wid_id)] '''
        patch_empty_pos_ls = []
        for h_id, h_object in self.set.items():
            empty_site_pos_ls = h_object.empty_site_pos_ls
            for site_pos in empty_site_pos_ls:
                site_pos = (h_id, ) + site_pos
                patch_empty_pos_ls.append(site_pos)
        return patch_empty_pos_ls

    def patch_empty_sites_num(self):
        ''' return the number of empty microsite in the patches '''
        return len(self.get_patch_empty_sites_ls())

    def get_patch_pairwise_empty_sites_ls(self):
        ''' return patch_empty_pos_ls as [((h_id, len_id, wid_id), (h_id, len_id, wid_id))...] '''
        patch_pairwise_empty_pos_ls = []
        for h_id, h_object in self.set.items():
            pairwise_empty_sites_pos_ls = h_object.get_hab_pairwise_empty_site_pos_ls()
            for (empty_site_1_pos, empty_site_2_pos) in pairwise_empty_sites_pos_ls:
                empty_site_1_pos = (h_id, ) + empty_site_1_pos
                empty_site_2_pos = (h_id, ) + empty_site_2_pos
                patch_pairwise_empty_pos_ls.append((empty_site_1_pos, empty_site_2_pos))
        return patch_pairwise_empty_pos_ls
    
    def get_patch_pairwise_occupied_sites_ls(self):
        ''' return patch_empty_pos_ls as [((h_id, len_id, wid_id), (h_id, len_id, wid_id))...] '''
        patch_pairwise_occupied_pos_ls = []
        for h_id, h_object in self.set.items():
            pairwise_occupied_sites_pos_ls = h_object.get_hab_pairwise_occupied_site_pos_ls()
            for (occupied_site_1_pos, occupied_site_2_pos) in pairwise_occupied_sites_pos_ls:
                occupied_site_1_pos = (h_id, ) + occupied_site_1_pos
                occupied_site_2_pos = (h_id, ) + occupied_site_2_pos
                patch_pairwise_occupied_pos_ls.append((occupied_site_1_pos, occupied_site_2_pos))
        return patch_pairwise_occupied_pos_ls
    
    def get_patch_occupied_sites_ls(self):
        ''' return patch_occupied_pos_ls as [(h_id, len_id, wid_id)] '''
        patch_occupied_pos_ls = []
        for h_id, h_object in self.set.items():
            occupied_site_pos_ls = h_object.occupied_site_pos_ls
            for site_pos in occupied_site_pos_ls:
                site_pos = (h_id, ) + site_pos
                patch_occupied_pos_ls.append(site_pos)
        return patch_occupied_pos_ls
    
    def get_patch_offspring_marker_pool(self):
        ''' get the combination of all the offspring marker pool in all habitat in the patch 
        return as a list of offspring marker. '''
        patch_offsprings_marker_pool = []
        for h_id, h_object in self.set.items():
            patch_offsprings_marker_pool += h_object.offspring_marker_pool
        return patch_offsprings_marker_pool
    
    def get_patch_offspring_pool(self):
        ''' get the combination of all the offspring pool in all habitat in the patch 
        return as a list of offspring individual object. '''
        patch_offsprings_pool = []
        for h_id, h_object in self.set.items():
            patch_offsprings_pool += h_object.offspring_pool
        return patch_offsprings_pool
    
    def get_patch_dormancy_pool(self):
        ''' get the combination of all the dormancy pool in all habitat in the patch 
        return as a list of dormancy individual object. '''
        patch_dormancy_pool = []
        for h_id, h_object in self.set.items():
            patch_dormancy_pool += h_object.dormancy_pool
        return patch_dormancy_pool
    
    def get_patch_offspring_and_dormancy_pool(self):
        ''' get the combination of all the offspring_pool_and_dormancy_pool in all habitat in the patch 
        return as a list of offspring_and_dormancy individual object. '''
        patch_offspring_and_dormancy_pool = []
        for h_id, h_object in self.set.items():
            patch_offspring_and_dormancy_pool += h_object.offspring_pool
            patch_offspring_and_dormancy_pool += h_object.dormancy_pool
        return patch_offspring_and_dormancy_pool

    def sample_offspring_without_replacement(self, num, pool_name=('offspring_pool', )):
        ''' Sample `num` offspring (individual_objs or offspring_markers) across all habitats of this patch 
        and REMOVE them from their source pool, then return a Python-list of sampled offspring (objs or markers).
        pool_name could be ('offspring_pool', 'dormancy_pool') or ('offspring_marker_pool') '''
        idx_ls = [(h_object, name, ObjOrMarker) for h_object in self.set.values() for name in pool_name for ObjOrMarker in getattr(h_object, name)]
        picked_ls = random.sample(idx_ls, min(num, len(idx_ls)))
        for h_object, name, ObjOrMarker in picked_ls: 
            getattr(h_object, name).remove(ObjOrMarker)
        return [ObjOrMarker for _, _, ObjOrMarker in picked_ls]
    
    def get_patch_microsites_individals_sp_id_values(self):
        ''' get species_id, phenotypes distribution in the patch as values set '''
        values_set = np.array([], dtype=int)
        for h_id, h_object in self.set.items():
            hab_len = h_object.length
            hab_wid = h_object.width
            for row in range(hab_len):
                for col in range(hab_wid):
                    individual_object = h_object.set['microsite_individuals'][row][col]
                    if individual_object ==None:
                        values_set = np.append(values_set, np.nan)
                    else:
                        species_id = individual_object.species_id
                        species_id_value = int(re.findall(r"\d+",species_id)[0])
                        values_set = np.append(values_set, species_id_value)
        values_set = values_set.reshape(self.hab_num, h_object.size)
        return values_set
    
    def get_patch_microsites_individals_phenotype_values(self, trait_name):
        ''' get species_id, phenotypes distribution in the patch as values set '''
        values_set = []
        for h_id, h_object in self.set.items():
            hab_len = h_object.length
            hab_wid = h_object.width
            for row in range(hab_len):
                for col in range(hab_wid):
                    individual_object = h_object.set['microsite_individuals'][row][col]
                    if individual_object ==None:
                        values_set.append(np.nan)
                    else:
                        phenotype = individual_object.phenotype_set[trait_name]
                        values_set.append(phenotype)
        values_set = np.array(values_set).reshape(self.hab_num, h_object.size)
        return values_set
    
    def get_patch_microsites_environment_values(self, environment_name):
        ''' get microsite environment values distribution in the patch as values set '''
        values_array = np.array([])
        for h_id, h_object in self.set.items():
            hab_environment_values_array = h_object.set[environment_name] 
            values_array = np.append(values_array, hab_environment_values_array) # dimension of the return of np.append() is always in dim=1
        values_array = values_array.reshape(self.hab_num, h_object.size)
        return values_array

# ********************* changing environmental values (offset) ********************************
    def patch_offset_environmental_values(self, env_name_ls, delta_mean_ls, delta_var_ls=None):
        ''' multi mean (var) values of env_name_ls in habitats + or - delta_mean_ls (or delta_val_ls) across the patches '''
        if delta_var_ls == None: delta_var_ls = [0 for i in range(len(env_name_ls))]
        for h_id, h_object in self.set.items():
            h_object.hab_offset_environment_values(env_name_ls, delta_mean_ls, delta_var_ls)
        return 0
  
##********************* for patch in mainland only and do not contains dormant bank *******************************####
    def patch_asexual_birth_germinate(self, asexual_birth_rate):
        ''' NemoBackend: asexual birth directly into empty microsites in every habitat of the patch, without competition between local offspring and immigrants
        INPUT: asexual_birth_rate
        RETURN: number germinated '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_asexual_reprodece_germinate(asexual_birth_rate)
        return counter
    
    def patch_sexual_birth_germinate(self, sexual_birth_rate):
        ''' NemoBackend: sexual birth directly into empty microsites in every habitat of the patch, without competition between local offspring and immigrants
        INPUT: sexual_birth_rate
        RETURN: number germinated '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_sexual_reprodece_germinate(sexual_birth_rate)
        return counter
    
    def patch_mixed_birth_germinate(self, asexual_birth_rate, sexual_birth_rate):
        ''' NemoBackend: mixed birth (asexual + sexual) directly into empty microsites in every habitat of the patch
        INPUT: asexual_birth_rate, sexual_birth_rate
        RETURN: number germinated '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_mixed_reproduce_germinate(asexual_birth_rate, sexual_birth_rate)
        return counter
    
#** calculating the offspring num and generating offspring markers but do not reproduce actually until local germination process afterward after dispersal process.*
    def patch_asex_reproduce_calculation_into_offspring_marker_pool(self, asexual_birth_rate):
        ''' used only when dormancy process do not occur '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_asex_reproduce_calculation_into_offspring_marker_pool(self.name, asexual_birth_rate)
        return counter
    
    def patch_sex_reproduce_calculation_into_offspring_marker_pool(self, sexual_birth_rate):
        ''' used only when dormancy process do not occur '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_sex_reproduce_calculation_into_offspring_marker_pool(self.name, sexual_birth_rate)
        return counter
    
    def patch_mix_reproduce_calculation_into_offspring_marker_pool(self, asexual_birth_rate, sexual_birth_rate):
        ''' used only when dormancy process do not occur '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_mix_reproduce_calculation_into_offspring_marker_pool(self.name, asexual_birth_rate, sexual_birth_rate)
        return counter
        
#********************************************************************************************************************#
#************************ all hab in patch reproduce into offspring pool ********************************************#
    def patch_asex_reproduce_mutate_into_offspring_pool(self, asexual_birth_rate):
        ''' NemoBackend: asexual birth into offspring_pool in every habitat of the patch '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_asex_reproduce_mutate_into_offspring_pool(asexual_birth_rate)
        return counter
    
    def patch_sex_reproduce_mutate_into_offspring_pool(self, sexual_birth_rate):
        ''' NemoBackend: sexual birth into offspring_pool in every habitat of the patch '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_sex_reproduce_mutate_into_offspring_pool(sexual_birth_rate)
        return counter
    
    def patch_mix_reproduce_mutate_into_offspring_pool(self, asexual_birth_rate, sexual_birth_rate):
        ''' NemoBackend: mixed birth (asexual + sexual) into offspring_pool in every habitat of the patch '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_mix_reproduce_mutate_into_offspring_pool(asexual_birth_rate, sexual_birth_rate)
        return counter
#****************************************** dispersal among patch into habs offspring pool *******************************************************
    def split(self, ls, n):
        ''' split a ls into n pieces '''
        k, m = divmod(len(ls), n)
        divided_ls = list(ls[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))
        random.shuffle(divided_ls)
        return divided_ls

    def migrants_to_patch_into_habs_immigrant_pool(self, migrants_indi_obj_to_patch_ls):
        divide_into_habs_indi_obj_ls = self.split(migrants_indi_obj_to_patch_ls, self.hab_num)
        for h_id, h_object in self.set.items():
            h_object.immigrant_pool += divide_into_habs_indi_obj_ls[h_object.index]
        return 0
    
    def migrants_marker_to_patch_into_habs_immigrant_marker_pool(self, migrants_to_j_offspring_marker_ls):
        divide_into_habs_offs_marker_ls = self.split(migrants_to_j_offspring_marker_ls, self.hab_num)
        for h_id, h_object in self.set.items():
            h_object.immigrant_marker_pool += divide_into_habs_offs_marker_ls[h_object.index]
        return 0
#****************************************** dispersal within patch ******************************************************
    def patch_matrix_around(self, matrix):   
        ''' Stochastic Rounding: fractional parts as changes (per-row: systematic sampling & sum kept); NaN -> 0 '''
        a = np.where(np.isnan(matrix), 0.0, np.asarray(matrix, dtype=float))
        base = np.floor(a); frac = a - base
        c  = np.round(np.cumsum(frac, axis=1), 9)
        c0 = np.hstack([np.zeros((a.shape[0], 1)), c[:, :-1]])
        u  = 1.0 - np.random.uniform(0, 1, (a.shape[0], 1))
        matrix[:, :] = base + np.floor(c - u + 1) - np.floor(c0 - u + 1)
        return matrix

    def get_patch_dispersal_within_rate_matrix(self, disp_within_rate):
        dispersal_within_rate_matrix = np.matrix(np.ones((self.hab_num, self.hab_num)))*disp_within_rate/(self.hab_num-1)
        dispersal_within_rate_matrix[np.diag_indices_from(dispersal_within_rate_matrix)] = 1-disp_within_rate  # elements in diagonal_line = 1-m
        return dispersal_within_rate_matrix

    def get_patch_habs_offspring_num_matrix(self):
        habs_offspring_num_matrix = np.matrix(np.zeros((self.hab_num, self.hab_num)))
        for h_id, h_object in self.set.items():
            offspring_num = len(h_object.offspring_pool)
            habs_offspring_num_matrix[h_object.index, h_object.index] = offspring_num
        return habs_offspring_num_matrix
    
    def get_patch_habs_offspring_marker_num_matrix(self):
        habs_offspring_marker_num_matrix = np.matrix(np.zeros((self.hab_num, self.hab_num)))
        for h_id, h_object in self.set.items():
            offspring_marker_num = len(h_object.offspring_marker_pool)
            habs_offspring_marker_num_matrix[h_object.index, h_object.index] = offspring_marker_num
        return habs_offspring_marker_num_matrix
    
    def get_patch_habs_dormancy_num_matrix(self):
        habs_dormancy_num_matrix = np.matrix(np.zeros((self.hab_num, self.hab_num)))
        for h_id, h_object in self.set.items():
            dormancy_num = len(h_object.dormancy_pool)
            habs_dormancy_num_matrix[h_object.index, h_object.index] = dormancy_num
        return habs_dormancy_num_matrix

    def get_patch_habs_empty_sites_num_matrix(self):
        habs_empty_sites_num_matrix = np.matrix(np.zeros((self.hab_num, self.hab_num)))
        for h_id, h_object in self.set.items():
            empty_site_num = len(h_object.empty_site_pos_ls)
            habs_empty_sites_num_matrix[h_object.index, h_object.index] = empty_site_num
        return habs_empty_sites_num_matrix
    
    def get_habs_emigrants_matrix(self, disp_within_rate):
        habs_offspring_num_matrix = self.get_patch_habs_offspring_num_matrix()
        habs_dormancy_num_matrix = self.get_patch_habs_dormancy_num_matrix()
        dispersal_within_rate_matrix = self.get_patch_dispersal_within_rate_matrix(disp_within_rate)
        return (habs_offspring_num_matrix + habs_dormancy_num_matrix) * dispersal_within_rate_matrix
    
    def get_habs_immigrants_matrix(self, disp_within_rate):
        habs_emigrants_matrix = self.get_habs_emigrants_matrix(disp_within_rate)
        habs_empty_sites_num_matrix = self.get_patch_habs_empty_sites_num_matrix()
        return (habs_emigrants_matrix/habs_emigrants_matrix.sum(axis=0))*habs_empty_sites_num_matrix
    
    def get_dispersal_within_num_matrix(self, disp_within_rate):
        habs_emigrants_matrix = self.get_habs_emigrants_matrix(disp_within_rate)
        habs_immigrants_matrix = self.get_habs_immigrants_matrix(disp_within_rate)
        return self.patch_matrix_around(np.minimum(habs_emigrants_matrix, habs_immigrants_matrix))

    def patch_dipersal_within_from_offspring_and_dormancy_pool(self, disp_within_rate, is_remove=True):
        ''' NemoBackend: dispersal within the patch -- individuals from habitat i's offspring + dormancy pools germinate directly into empty microsites of habitat j (i != j); migrants left without an empty site are recycled
        INPUT: disp_within_rate (m_within); is_remove = True (default): sampling without replacement, migrants leave their source pool / False: sampling with replacement across destinations, the original stays in its source pool and a copy disperses (individual.make_copy())
        RETURN: number germinated (after dispersal) '''
        disp_within_num_matrix = self.get_dispersal_within_num_matrix(disp_within_rate)
        counter = 0
        for j in range(self.hab_num):
            h_j_id = self.hab_id_ls[j]
            h_j_object = self.set[h_j_id]
            h_j_empty_site_ls = h_j_object.empty_site_pos_ls
            migrants_to_j_indi_object_ls = []
            
            for i in range(self.hab_num):
                h_i_id = self.hab_id_ls[i]
                h_i_object = self.set[h_i_id]
                if i==j:
                    continue
                else:
                    migrants_i_j_num = int(disp_within_num_matrix[i, j])
                    if is_remove: migrants_to_j_indi_object_ls += h_i_object.sample_offspring_without_replacement(migrants_i_j_num, pool_name=('offspring_pool', 'dormancy_pool'))
                    else: migrants_to_j_indi_object_ls += [indi_object.make_copy() for indi_object in random.sample(h_i_object.offspring_pool + h_i_object.dormancy_pool, migrants_i_j_num)]   # NemoBackend Method: individual.make_copy(), the originals stay in their source pools
                    
            random.shuffle(h_j_empty_site_ls)
            random.shuffle(migrants_to_j_indi_object_ls)
            for (len_id, wid_id), migrants_object in list(zip(h_j_empty_site_ls, migrants_to_j_indi_object_ls)):
                self.set[h_j_id].add_individual(indi_object = migrants_object, len_id=len_id, wid_id=wid_id)
                counter += 1
                migrants_to_j_indi_object_ls.remove(migrants_object)

            for migrants_object in migrants_to_j_indi_object_ls:                      # no empty site left: already out of its source pool (is_remove=False: a copy, the original stays in its source pool), dies
                migrants_object.recycle()                                             # NemoBackend Methods: recycling individuals
        return counter
    
    def patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool(self, disp_within_rate, is_remove=True):
        offspring_marker_dispersal_matrix =  self.patch_matrix_around(self.get_patch_habs_offspring_marker_num_matrix()*self.get_patch_dispersal_within_rate_matrix(disp_within_rate))
        counter = 0
        for j in range(self.hab_num):
            h_j_id = self.hab_id_ls[j]
            h_j_object = self.set[h_j_id]
            migrants_to_j_offspring_marker_ls = []
            
            for i in range(self.hab_num):
                h_i_id = self.hab_id_ls[i]
                h_i_object = self.set[h_i_id]
                if i==j:
                    continue
                else:
                    migrants_i_j_num = int(offspring_marker_dispersal_matrix[i, j])
                    if is_remove: migrants_to_j_offspring_marker_ls += h_i_object.sample_offspring_without_replacement(migrants_i_j_num, pool_name=('offspring_marker_pool',))
                    else: migrants_to_j_offspring_marker_ls += random.sample(h_i_object.offspring_marker_pool, migrants_i_j_num)
                    
            random.shuffle(migrants_to_j_offspring_marker_ls)
            h_j_object.immigrant_marker_pool += migrants_to_j_offspring_marker_ls
            counter += len(migrants_to_j_offspring_marker_ls)
        return counter
    
    def patch_dispersal_within_from_offspring_pool_to_immigrant_pool(self, disp_within_rate, is_remove=True):
        ''' NemoBackend: dispersal within the patch: offspring from habitat i's offspring_pool move into habitat j's immigrant_pool (i != j)
        INPUT: disp_within_rate (m_within); is_remove = True (default): sampling without replacement, migrants leave their source pool / False: sampling with replacement across destinations, the original stays in its source pool and a copy disperses (individual.make_copy())
        RETURN: number of migrants '''
        offspring_dispersal_matrix = self.patch_matrix_around(self.get_patch_habs_offspring_num_matrix()*self.get_patch_dispersal_within_rate_matrix(disp_within_rate))
        counter = 0
        for j in range(self.hab_num):
            h_j_id = self.hab_id_ls[j]
            h_j_object = self.set[h_j_id]
            migrants_to_j_indi_object_ls = []
            
            for i in range(self.hab_num):
                h_i_id = self.hab_id_ls[i]
                h_i_object = self.set[h_i_id]
                if i==j: 
                    continue
                else:
                    migrants_i_j_num = int(offspring_dispersal_matrix[i, j])
                    if is_remove: migrants_to_j_indi_object_ls += h_i_object.sample_offspring_without_replacement(migrants_i_j_num)
                    else: migrants_to_j_indi_object_ls += [indi_object.make_copy() for indi_object in random.sample(h_i_object.offspring_pool, migrants_i_j_num)]   # NemoBackend Method: individual.make_copy(), the originals stay in their source pool
                    
            random.shuffle(migrants_to_j_indi_object_ls)
            h_j_object.immigrant_pool += migrants_to_j_indi_object_ls
            counter += len(migrants_to_j_indi_object_ls)
        return counter
    
    def patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool(self, disp_within_rate, is_remove=True):
        ''' NemoBackend: dispersal within the patch (m_within): individuals from habitat i's offspring + dormancy pools move into habitat j's immigrant_pool (i != j)
        INPUT: disp_within_rate (m_within); is_remove = True (default): sampling without replacement, migrants leave their source pool / False: sampling with replacement across destinations, the original stays in its source pool and a copy disperses (individual.make_copy())
        RETURN: number of migrants '''
        offspring_and_dormancy_dispersal_matrix = self.patch_matrix_around(self.get_habs_emigrants_matrix(disp_within_rate))
        counter = 0
        for j in range(self.hab_num):
            h_j_id = self.hab_id_ls[j]
            h_j_object = self.set[h_j_id]
            migrants_to_j_indi_object_ls = []
            
            for i in range(self.hab_num):
                h_i_id = self.hab_id_ls[i]
                h_i_object = self.set[h_i_id]
                if i==j:
                    continue
                else:
                    migrants_i_j_num = int(offspring_and_dormancy_dispersal_matrix[i, j])
                    if is_remove: migrants_to_j_indi_object_ls += h_i_object.sample_offspring_without_replacement(migrants_i_j_num, pool_name=('offspring_pool', 'dormancy_pool'))
                    else: migrants_to_j_indi_object_ls += [indi_object.make_copy() for indi_object in random.sample(h_i_object.offspring_pool + h_i_object.dormancy_pool, migrants_i_j_num)]   # NemoBackend Method: individual.make_copy(), the originals stay in their source pools
                    
            random.shuffle(migrants_to_j_indi_object_ls)
            h_j_object.immigrant_pool += migrants_to_j_indi_object_ls
            counter += len(migrants_to_j_indi_object_ls)
        return counter  
        
#******************** local germination in all habitat in the patch *****************************#
    def patch_local_germinate_from_offspring_and_dormancy_pool(self, is_remove=True):
        ''' NemoBackend: germinate individuals from offspring_pool + dormancy_pool into empty microsites in every habitat of the patch
        INPUT: is_remove = True (default): germinated individuals leave their pool / False: a copy germinates, the original stays in its pool
        RETURN: number germinated '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_local_germinate_from_offspring_and_dormancy_pool(is_remove)
        return counter
    
    def patch_local_germinate_from_offspring_and_immigrant_pool(self, is_remove=True):
        ''' NemoBackend: same as patch_local_germinate_from_offspring_and_dormancy_pool(), from offspring_pool + immigrant_pool '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_local_germinate_from_offspring_and_immigrant_pool(is_remove)
        return counter
            
    def patch_local_germinate_from_offspring_immigrant_and_dormancy_pool(self, is_remove=True):
        ''' NemoBackend: same as patch_local_germinate_from_offspring_and_dormancy_pool(), from offspring_pool + immigrant_pool + dormancy_pool '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_local_germinate_from_offspring_immigrant_and_dormancy_pool(is_remove)
        return counter

#********************************** dormancy process in the patch **********************************
    def patch_dormancy_process_from_offspring_pool_to_dormancy_pool(self):
        ''' NemoBackend: dormancy process in every habitat of the patch -- offspring enter the dormancy pool (capacity per habitat); immigrants and everyone else dropped from the pools are recycled
        RETURN: (survival_num, eliminate_num, new_dormancy_num, dormancy_num) summed over the habitats '''
        survival_counter = 0
        eliminate_counter = 0
        new_dormancy_counter = 0
        all_dormancy_num = 0
        for h_id, h_object in self.set.items():
            survival_num, eliminate_num, new_dormancy_num, dormancy_num = h_object.hab_dormancy_process_from_offspring_pool_to_dormancy_pool() # NemoBackend Methods: recycling the dropped
            survival_counter += survival_num
            eliminate_counter += eliminate_num
            new_dormancy_counter += new_dormancy_num
            all_dormancy_num += dormancy_num
        return survival_counter, eliminate_counter, new_dormancy_counter, all_dormancy_num
    
    def patch_dormancy_process_from_offspring_pool_and_immigrant_pool(self):
        ''' NemoBackend: dormancy process in every habitat of the patch -- offspring and immigrants enter the dormancy pool (capacity per habitat; priority offspring > immigrants > old dormants); everyone dropped is recycled
        RETURN: (survival_num, eliminate_num, new_dormancy_num, dormancy_num) summed over the habitats '''
        survival_counter = 0
        eliminate_counter = 0
        new_dormancy_counter = 0
        all_dormancy_num = 0
        for h_id, h_object in self.set.items():
            survival_num, eliminate_num, new_dormancy_num, dormancy_num = h_object.hab_dormancy_process_from_offspring_pool_and_immigrant_pool() # NemoBackend Methods: recycling the dropped
            survival_counter += survival_num
            eliminate_counter += eliminate_num
            new_dormancy_counter += new_dormancy_num
            all_dormancy_num += dormancy_num
        return survival_counter, eliminate_counter, new_dormancy_counter, all_dormancy_num    
    
# ******************** clearing up microsites, offspring pool and immigrant pool and dormancy pool ***************************************************#
    def patch_clear_up_microsite_individuals(self):
        ''' NemoBackend: recycle every individual in the microsites of every habitat of the patch and reset the microsite bookkeeping (pools untouched)
        RETURN: number recycled '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_clear_up_microsite_individuals()   # NemoBackend Methods: recycling individuals
        return counter

    def patch_clear_up_offspring_pool(self):
        ''' NemoBackend: recycle every ungerminated offspring in every habitat of the patch, then empty offspring_pool
        RETURN: number recycled '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_clear_up_offspring_pool()          # NemoBackend Methods: recycling individuals
        return counter

    def patch_clear_up_immigrant_pool(self):
        ''' NemoBackend: recycle every ungerminated immigrant in every habitat of the patch, then empty immigrant_pool
        RETURN: number recycled '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_clear_up_immigrant_pool()          # NemoBackend Methods: recycling individuals
        return counter

    def patch_clear_up_dormancy_pool(self):
        ''' NemoBackend: recycle every dormant individual in every habitat of the patch, then empty dormancy_pool -- teardown only
        (disturbance keeps the dormancy pool; never call it inside the dormancy processes, where _recycle_dropped already recycles)
        RETURN: number recycled '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_clear_up_dormancy_pool()           # NemoBackend Methods: recycling individuals
        return counter

    def patch_clear_up_offspring_marker_pool(self):
        ''' empty offspring_marker_pool in every habitat of the patch (markers are the tags of the to-be-born: nothing to recycle)
        RETURN: number of markers cleared '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_clear_up_offspring_marker_pool()
        return counter

    def patch_clear_up_immigrant_marker_pool(self):
        ''' empty immigrant_marker_pool in every habitat of the patch (markers are the tags of the to-be-born: nothing to recycle)
        RETURN: number of markers cleared '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.hab_clear_up_immigrant_marker_pool()
        return counter
    
    def patch_disturbance_process(self):
        ''' NemoBackend: disturbance in every habitat of the patch -- individuals in the microsites and in the offspring / immigrant pools are recycled (dormancy_pool survives), marker pools are emptied
        RETURN: number of individuals recycled (both microsites and pool) + number of markers cleared '''
        counter = 0
        for h_id, h_object in self.set.items():
            counter += h_object.habitat_disturbance_process()          # NemoBackend Methods: recycling individuals
        return counter
            
#********************************** data generating for saving ******************************************
    def get_patch_microsites_optimum_sp_id_value_array(self, d, w, species_2_phenotype_ls):
        ''''''
        values_array = np.array([], dtype=int)
        for h_id, h_object in self.set.items():
            hab_len = h_object.length
            hab_wid = h_object.width
            mean_env_tuple = h_object.mean_env_ls
            
            hab_opt_survival_rate = 0
            hab_opt_sp_id_val = np.nan
            for phenotype_ls in species_2_phenotype_ls:
                survival_rate = h_object.survival_rate(d=d, phenotype_ls=phenotype_ls, env_val_ls=mean_env_tuple, w=w)
                if survival_rate > hab_opt_survival_rate:
                    hab_opt_survival_rate = survival_rate
                    hab_opt_sp_id_val = species_2_phenotype_ls.index(phenotype_ls)+1

            hab_sp_id_val_array = np.ones(hab_len*hab_wid, dtype=int)*hab_opt_sp_id_val
            values_array = np.append(values_array, hab_sp_id_val_array)
        
        values_array = values_array.reshape(self.hab_num, h_object.size)
        return values_array
