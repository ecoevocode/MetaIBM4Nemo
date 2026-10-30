# -*- coding: utf-8 -*-
"""
Split from metacommunity_IBM.py.
Compatibility-preserving class module.
"""

import numpy as np
import random
import math
import copy

class habitat():
    def __init__(self, hab_name, hab_index, hab_location, num_env_types, env_types_name, mean_env_ls, var_env_ls, length, width, dormancy_pool_max_size=0, global_index=0):
        '''
        int num_env_types is the number of environment types in the habitat.
        env_types_name is the list of names of env_types.
        list mean_env_ls is the tuple of mean environment values; the len(mean_env_ls)=num_env_types.
        list var_env_ls is the list of variation of enviroment distribution in the habitat.
        int length is the length of the habitat.
        int width is the width of the habitat.
        int size is the the number of microsites within a habitat.
        '''
        self.name = hab_name
        self.index = hab_index            # unique index within patches
        self.global_index = global_index  # unique index across global metacommunity, used as 'homepatch' of Nemo's individual
        self.location = hab_location 
        self.num_env_types = num_env_types
        self.env_types_name = env_types_name
        self.mean_env_ls = list(mean_env_ls)   # own copy: hab_reset_environment_values() writes it in place
        self.var_env_ls = list(var_env_ls)     # own copy: callers often pass one shared list to every habitat
        self.length = length
        self.width = width
        self.size = length*width
        self.set = {}                     # self.data_set={} # to be improved
        self.indi_num = 0
        
        self.offspring_pool = []
        self.immigrant_pool = []
        self.dormancy_pool = []
        
        self.offspring_marker_pool = []   # used only without dormancy process, marker = (patch_id, h_id, reproduction_mode)
        self.immigrant_marker_pool = []   # used only without dormancy process, marker = (patch_id, h_id, reproduction_mode)
        
        self.species_category = {}
        self.occupied_site_pos_ls = []
        
        self.empty_site_pos_ls = [(i, j) for i in range(length) for j in range(width)]
        self.dormancy_pool_max_size = dormancy_pool_max_size
        
        self.reproduction_mode_threhold = 0.897
        self.asexual_parent_pos_ls = []                           # If an individual can fit its environment condition well, it goes through asexual reproduction. only for mix
        self.species_category_for_sexual_parents_pos = {}         # If an individual can not fit its environmet condition, it goes through sexual reproduction. only for mix

        for index in range(0, len(mean_env_ls)):
            mean_e_index = self.mean_env_ls[index]
            var_e_index = self.var_env_ls[index]
            name_e_index = self.env_types_name[index]
            microsite_e_values = np.random.normal(loc=0, scale=np.sqrt(var_e_index), size=(self.length, self.width)) + mean_e_index
            self.set[name_e_index] = microsite_e_values

        microsite_individuals = self._empty_microsite_grid()
        self.set['microsite_individuals'] = microsite_individuals
        
    def __str__(self):
        return str(self.set)
    
#************************** chaning environmental values (reset & offset) ******************************************* 
    def hab_reset_environment_values(self, env_name_ls, env_mean_val_ls, env_var_val_ls):
        ''' env_name_ls: a list of env name; env_mean_val_ls: a list of new mean values of the changing env; env_var_val_ls: a list of new var values of the changing env '''
        
        for name_i, mean_i, var_i in zip(env_name_ls, env_mean_val_ls, env_var_val_ls):
            env_index = self.env_types_name.index(name_i) # name to index
            self.mean_env_ls[env_index] = mean_i          # setting new mean value
            self.var_env_ls[env_index] = var_i            # setting new var value
            
            new_microsite_e_values = np.random.normal(loc=0, scale=np.sqrt(var_i), size=(self.length, self.width)) + mean_i # set the new env landscape of the habitat
            self.set[name_i] = new_microsite_e_values # update all the microsites env val in the habitat
            
        return 0
           
    def hab_offset_environment_values(self, env_name_ls, delta_mean_ls, delta_var_ls=None):
        ''' env_name_ls: a list of the name of changing environment;
            delta_mean_ls: a list of changing (mean) values of the changing environment, can be both + or - 
            delta_var_ls: a list of changing (var) values of the changing environment, can be both + or - '''
        
        if delta_var_ls == None: delta_var_ls = [0 for i in range(len(env_name_ls))]
        
        new_env_mean_ls = [v + d for v, d in zip(self.mean_env_ls, delta_mean_ls)]
        new_env_var_ls = [v + d for v, d in zip(self.var_env_ls, delta_var_ls)]
        
        self.hab_reset_environment_values(env_name_ls=env_name_ls, env_mean_val_ls=new_env_mean_ls, env_var_val_ls=new_env_var_ls)
        
        return 0
    
#****************************************************************************************
    def _empty_microsite_grid(self):        
        """Return an empty microsite grid with shape [length][width]."""        
        return [[None for _ in range(self.width)] for _ in range(self.length)]

    def add_individual(self, indi_object, len_id, wid_id):
        if self.set['microsite_individuals'][len_id][wid_id] != None:
            print('the microsite in the habitat is occupied.')
        else:
            self.set['microsite_individuals'][len_id][wid_id] = indi_object
            self.empty_site_pos_ls.remove((len_id, wid_id))
            self.occupied_site_pos_ls.append((len_id, wid_id))
            self.indi_num +=1

            if indi_object.species_id in self.species_category.keys():
                if indi_object.gender in self.species_category[indi_object.species_id].keys():
                    self.species_category[indi_object.species_id][indi_object.gender].append((len_id, wid_id))
                else:
                    self.species_category[indi_object.species_id][indi_object.gender] = [(len_id, wid_id)]             
            else:
                self.species_category[indi_object.species_id] = {indi_object.gender:[(len_id, wid_id)]}
                
    def del_individual(self, len_id, wid_id):
        ''' the individual leaves the microsite, not necessarily dead (no recycle()) '''
        if self.set['microsite_individuals'][len_id][wid_id] == None:
            print('the microsite in the habitat is empty.')
        else:
            indi_object = self.set['microsite_individuals'][len_id][wid_id]
            self.set['microsite_individuals'][len_id][wid_id] = None
            self.empty_site_pos_ls.append((len_id, wid_id))
            self.occupied_site_pos_ls.remove((len_id, wid_id))
            self.indi_num -=1 
            self.species_category[indi_object.species_id][indi_object.gender].remove((len_id, wid_id))    
        
    def get_hab_pairwise_empty_site_pos_ls(self):
        ''' return as [((len_id, wid_id), (len_id, wid_id)) ...]'''
        hab_pairwise_empty_sites_pos_ls = []
        if len(self.empty_site_pos_ls) < 2:
            return hab_pairwise_empty_sites_pos_ls
        else:
            empty_sites_pos_ls = copy.deepcopy(self.empty_site_pos_ls) 
            random.shuffle(empty_sites_pos_ls)
            for i in range(0, len(empty_sites_pos_ls)-1, 2):
                empty_site_1_pos = empty_sites_pos_ls[i]
                empty_site_2_pos = empty_sites_pos_ls[i+1]
                
                hab_pairwise_empty_sites_pos_ls.append((empty_site_1_pos, empty_site_2_pos))
            return hab_pairwise_empty_sites_pos_ls    
    
    def get_hab_pairwise_occupied_site_pos_ls(self):
        ''' return as [((len_id, wid_id), (len_id, wid_id)) ...]'''
        hab_pairwise_occupied_sites_pos_ls = []
        if self.indi_num < 2:
            return hab_pairwise_occupied_sites_pos_ls
        else:
            species_category = copy.deepcopy(self.species_category)
            for sp_id, sp_id_val in species_category.items():
                try:
                    sp_id_female_ls = sp_id_val['female']
                    random.shuffle(sp_id_female_ls)
                except:
                    continue
                try:
                    sp_id_male_ls = sp_id_val['male']
                    random.shuffle(sp_id_male_ls)
                except:
                    continue
                try:
                    sp_id_pairwise_occupied_pos_ls = list(zip(sp_id_female_ls, sp_id_male_ls))
                    hab_pairwise_occupied_sites_pos_ls += sp_id_pairwise_occupied_pos_ls
                except:
                    continue
            return hab_pairwise_occupied_sites_pos_ls

    def sample_offspring_without_replacement(self, num, pool_name=('offspring_pool', )):
        ''' Sample `num` offspring (individual_objs or offspring_markers) of this habitat
        and REMOVE them from their source pool, then return a Python-list of sampled offspring (objs or markers).
        pool_name could be ('offspring_pool', 'dormancy_pool') or ('offspring_marker_pool') '''

        idx_ls = [(name, ObjOrMarker) for name in pool_name for ObjOrMarker in getattr(self, name)]
        num = min(num, len(idx_ls))
        picked_ls = random.sample(idx_ls, num)
        for name, ObjOrMarker in picked_ls: 
            getattr(self, name).remove(ObjOrMarker)
        return [ObjOrMarker for _, ObjOrMarker in picked_ls]

    def hab_initialize(self, indfactories, species_2_phenotype_ls, reproduce_mode='sexual'):
        ''' NemoBackend: reset the microsites (residents recycled, pools untouched), then fill every microsite with a founder of the species matched to this habitat's mean environment
        INPUT: indfactories {species_id: indfactory}, species_2_phenotype_ls, reproduce_mode 'sexual' | 'asexual'
        RETURN: 0 '''
        species_id = 'sp%d' % (species_2_phenotype_ls.index(self.mean_env_ls) + 1)
        indfactory = indfactories[species_id]
        self.hab_clear_up_microsite_individuals()    # NemoBackend Methods: recycle residents, empty all microsites (re-initializing mid-simulation is rare, but safe)
        for row in range(self.length):
            for col in range(self.width):
                if reproduce_mode in ('sexual', 'mix'): gender = random.sample(('male', 'female'), 1)[0]
                elif reproduce_mode == 'asexual': gender = 'female'
                self.add_individual(indfactory.make_founder(gender, self.global_index), row, col)
        return 0 
    
    def get_microsite_env_val_ls(self, len_id, wid_id):
        ''' return a list of environment value of all the environment type in the order of env_types_name '''
        env_val_ls = []
        for env_name in self.env_types_name:
            env_val = self.set[env_name][len_id][wid_id]
            env_val_ls.append(env_val)
        return env_val_ls
    
    def survival_rate(self, d, phenotype_ls, env_val_ls, w=0.5, method='niche_gaussian'):
        #d is the baseline death rate responding to the disturbance strength.
        #phenotype_ls is a list of phenotype of each trait.
        #env_val_ls is a list of environment value responding to the environment type.
        #w is the width of the fitness function.
        
        
        if method == 'niche_gaussian':
            survival_rate = (1-d)
            for index in range(len(phenotype_ls)):
                ei = phenotype_ls[index]               # individual phenotype of a trait 
                em = env_val_ls[index]                 # microsite environment value of a environment type
                survival_rate *= math.exp((-1)*math.pow(((ei-em)/w),2))
                
        elif method == 'neutral_uniform':
            survival_rate = 1-d

        return survival_rate

    def hab_dead_selection(self, base_dead_rate, fitness_wid, method):
        ''' NemoBackend '''
        self.asexual_parent_pos_ls = []                           # If an individual can fit its environment condition well, it goes through asexual reproduction.
        self.species_category_for_sexual_parents_pos = {}         # If an individual can not fit its environmet condition, it goes through sexual reproduction.
        counter = 0
        for row in range(self.length):
            for col in range(self.width):
                env_val_ls = self.get_microsite_env_val_ls(row, col)
                
                if self.set['microsite_individuals'][row][col] != None:
                    individual_object = self.set['microsite_individuals'][row][col]
                    phenotype_ls = individual_object.get_indi_phenotype_ls()
                    survival_rate = self.survival_rate(d=base_dead_rate, phenotype_ls=phenotype_ls, env_val_ls=env_val_ls, w=fitness_wid, method=method)
                    
                    if survival_rate < np.random.uniform(0,1,1)[0]:
                        self.del_individual(len_id=row, wid_id=col)           # bookkeeping first (del_individual reads gender)
                        individual_object.recycle()                           # set individual.is_alive = False; move the nemo_individual into Nemo's recycling-pool
                        counter += 1
                    else:
                        if survival_rate >= self.reproduction_mode_threhold: 
                            self.asexual_parent_pos_ls.append((row, col))     # the individual fits its local environment
                        else:
                            if individual_object.species_id in self.species_category_for_sexual_parents_pos.keys():
                                if individual_object.gender in self.species_category_for_sexual_parents_pos[individual_object.species_id].keys():
                                    self.species_category_for_sexual_parents_pos[individual_object.species_id][individual_object.gender].append((row, col))
                                else:
                                    self.species_category_for_sexual_parents_pos[individual_object.species_id][individual_object.gender] = [(row, col)]
                            else:
                                self.species_category_for_sexual_parents_pos[individual_object.species_id] = {individual_object.gender:[(row, col)]}         
                else:
                    continue
        return counter
#*************** for habitat for mainland only and do not contains dormant bank ******************************************
    def hab_asex_reproduce_mutate_with_num(self, num):
        ''' NemoBackend: asexual birth of num offspring; INPUT: num = number of offspring; RETURN: list of offspring indi_obj '''
        list_of_offspring = []
        mother_pos_ls = random.sample(self.occupied_site_pos_ls, num)
        for row, col in mother_pos_ls:
            mother_indi_obj = self.set['microsite_individuals'][row][col]
            list_of_offspring.append(mother_indi_obj.indfactory.make_clone(mother_indi_obj, home_habitat=self.global_index)) # Nemobackended Method: indfactory::makeClonedOffsprg(), home_patch = home_habitat
        return list_of_offspring
    
    def hab_sexual_pairwise_parents_ls(self):
        pair_parents_ls = []
        for sp_id, sp_id_val in self.species_category.items():
            try:
                sp_id_female_ls = sp_id_val['female']
            except:
                continue
            try:
                sp_id_male_ls = sp_id_val['male']
            except:
                continue
            
            random.shuffle(sp_id_female_ls) #list of individuals location in habitat, i.e., (len_id, wid_id)
            random.shuffle(sp_id_male_ls)   #random sample of pairwise parents in sexual reproduction
            
            pair_parents_ls += list(zip(sp_id_female_ls, sp_id_male_ls))
        return pair_parents_ls
    
    def hab_sexual_pairwise_parents_num(self):
        return len(self.hab_sexual_pairwise_parents_ls())  
    
    def hab_sex_reproduce_mutate_with_num(self, num):
        ''' NemoBackend: sexual birth of num offspring from random parent pairs; inheritance, recombination and mutation by Nemo's makeOffsprg() (mutation and Ve set by the species' Nemo parameters)
        INPUT: num = number of offspring
        RETURN: list of offspring indi_obj '''
        list_of_offspring = []
        pairwise_parents_pos_ls = random.sample(self.hab_sexual_pairwise_parents_ls(), num)     # TODO(offspring_num): now num <= number of parent pairs (sampled without replacement: at most one offspring per pair); lift this so sexual_birth_rate is a true birth rate (offspring per pair, may be > 1)
        for female_pos, male_pos in pairwise_parents_pos_ls:
            female_indi_obj = self.set['microsite_individuals'][female_pos[0]][female_pos[1]]
            male_indi_obj = self.set['microsite_individuals'][male_pos[0]][male_pos[1]]
            gender = random.sample(('male', 'female'), 1)[0]
            list_of_offspring.append(female_indi_obj.indfactory.make_offspring(female_indi_obj, male_indi_obj, gender, home_habitat=self.global_index))
        return list_of_offspring
    
    def hab_mixed_sexual_pairwise_parents_ls(self):
        pair_parents_ls = []
        for sp_id, sp_id_val in self.species_category_for_sexual_parents_pos.items():
            try:
                sp_id_female_ls = sp_id_val['female']
            except:
                continue
            try:
                sp_id_male_ls = sp_id_val['male']
            except:
                continue
            
            random.shuffle(sp_id_female_ls) #list of individuals location in habitat, i.e., (len_id, wid_id)
            random.shuffle(sp_id_male_ls)   #random sample of pairwise parents in sexual reproduction
            
            pair_parents_ls += list(zip(sp_id_female_ls, sp_id_male_ls))
        return pair_parents_ls
    
    def hab_mixed_sexual_pairwse_parents_num(self):
        return len(self.hab_mixed_sexual_pairwise_parents_ls())
    
    def hab_mixed_asexual_parent_num(self):
        return len(self.asexual_parent_pos_ls)    
    
    def hab_mix_asex_reproduce_mutate_with_num(self, num):
        ''' NemoBackend: mixed mode, asexual part; INPUT: num = number of asexual offspring; RETURN: list of asexual offspring indi_obj '''
        list_of_offspring = []
        mother_pos_ls = random.sample(self.asexual_parent_pos_ls, num)
        for row, col in mother_pos_ls:
            mother_indi_obj = self.set['microsite_individuals'][row][col]
            list_of_offspring.append(mother_indi_obj.indfactory.make_clone(mother_indi_obj, home_habitat=self.global_index))
        return list_of_offspring
    
    def hab_mix_sex_reproduce_mutate_with_num(self, num):
        ''' NemoBackend: mixed mode, sexual part; INPUT: num = number of sexual offspring; RETURN: list of sexual offspring indi_obj '''
        list_of_offspring = []
        pairwise_parents_pos_ls = random.sample(self.hab_mixed_sexual_pairwise_parents_ls(), num)
        for female_pos, male_pos in pairwise_parents_pos_ls:
            female_indi_obj = self.set['microsite_individuals'][female_pos[0]][female_pos[1]]
            male_indi_obj   = self.set['microsite_individuals'][male_pos[0]][male_pos[1]]
            gender = random.sample(('male', 'female'), 1)[0]
            list_of_offspring.append(female_indi_obj.indfactory.make_offspring(female_indi_obj, male_indi_obj, gender, home_habitat=self.global_index))
        return list_of_offspring
##******************************* birth and local germination *************************************#####    
    def hab_asexual_reprodece_germinate(self, asexual_birth_rate):
        ''' NemoBackend: asexual birth directly into empty microsites (without competition between local offspring and immigrants)
        INPUT: asexual_birth_rate
        RETURN: number germinated (num <= empty sites) '''
        counter = 0
        num = min(len(self.empty_site_pos_ls), int(self.indi_num * asexual_birth_rate))                          # offspring_num <= empty sites
        hab_offsprings_for_germinate = self.hab_asex_reproduce_mutate_with_num(num)                              # NemoBackend Methods

        random.shuffle(self.empty_site_pos_ls)
        random.shuffle(hab_offsprings_for_germinate)

        for pos, indi_object in list(zip(self.empty_site_pos_ls, hab_offsprings_for_germinate)):
            len_id = pos[0]
            wid_id = pos[1]
            self.add_individual(indi_object, len_id, wid_id)
            counter += 1
        return counter
        
    def hab_sexual_reprodece_germinate(self, sexual_birth_rate):
        ''' NemoBackend: sexual birth directly into empty microsites (without competition between local offspring and immigrants)
        INPUT: sexual_birth_rate
        RETURN: number germinated (num <= empty sites) '''
        counter = 0
        num = min(len(self.empty_site_pos_ls), int(self.hab_sexual_pairwise_parents_num() * sexual_birth_rate)) # offspring_num <= empty sites
        hab_offsprings_for_germinate = self.hab_sex_reproduce_mutate_with_num(num)                              # NemoBackended methods
    
        random.shuffle(self.empty_site_pos_ls)
        random.shuffle(hab_offsprings_for_germinate)
        
        for pos, indi_object in list(zip(self.empty_site_pos_ls, hab_offsprings_for_germinate)):
            len_id = pos[0]
            wid_id = pos[1]
            self.add_individual(indi_object, len_id, wid_id)
            counter += 1
        return counter
    
    def hab_mixed_reproduce_germinate(self, asexual_birth_rate, sexual_birth_rate):
        ''' NemoBackend: mixed birth directly into empty microsites (asexual + sexual);
        INPUT: asexual_birth_rate, sexual_birth_rate
        RETURN: number germinated '''
        counter = 0
        empty_sites_pos_ls = self.empty_site_pos_ls
        empty_sites_num = len(empty_sites_pos_ls)
        
        asex_offs_expectation_num = int(np.around(self.hab_mixed_asexual_parent_num() * asexual_birth_rate))
        sex_offs_expectation_num = int(np.around(self.hab_mixed_sexual_pairwse_parents_num() * sexual_birth_rate))
        
        total = asex_offs_expectation_num + sex_offs_expectation_num
        asex_num = min(asex_offs_expectation_num, int(np.around(empty_sites_num * asex_offs_expectation_num / total))) if total > 0 else 0
        sex_num  = min(sex_offs_expectation_num, empty_sites_num - asex_num)
        
        hab_offsprings_for_germinate = self.hab_mix_asex_reproduce_mutate_with_num(asex_num) + self.hab_mix_sex_reproduce_mutate_with_num(sex_num) # NemoBackend Methods

        random.shuffle(self.empty_site_pos_ls)
        random.shuffle(hab_offsprings_for_germinate)

        for pos, indi_object in list(zip(empty_sites_pos_ls, hab_offsprings_for_germinate)):
            len_id = pos[0]
            wid_id = pos[1]
            self.add_individual(indi_object, len_id, wid_id)
            counter += 1
        return counter

#*****************  reproduction into offspring_marker_pool process *********************************
    def hab_asex_reproduce_calculation_into_offspring_marker_pool(self, patch_name, asexual_birth_rate):
        ''' 
        Calculating the offspring num and generating offspring markers 
        but do not reproduce actually until local germination process afterward after dispersal process 
        patch_id: the patch of the habitat belong with, (i.e., patch_object.patch_id)
        The function will be called by the the patch_object of the habitat belong with
        '''
        self.offspring_marker_pool = []
        offspring_num = self.indi_num * asexual_birth_rate
        offspring_num_int = int(offspring_num)                    #整数部分
        offspring_num_dem = offspring_num - offspring_num_int     #小数部分
        
        offspring_marker_ls = [(patch_name, self.name, 'asexual') for _ in range(offspring_num_int)]
        if offspring_num_dem > np.random.uniform(0,1,1)[0]:
            one_offspring_marker_ls = [(patch_name, self.name, 'asexual')]
        else:
            one_offspring_marker_ls = []
        self.offspring_marker_pool = offspring_marker_ls + one_offspring_marker_ls
        return len(self.offspring_marker_pool)

    def hab_sex_reproduce_calculation_into_offspring_marker_pool(self, patch_name, sexual_birth_rate):
        ''' 
        Calculating the offspring num and generating offspring markers 
        but do not reproduce actually until local germination process afterward after dispersal process 
        patch_id: the patch of the habitat belong with, (i.e., patch_object.patch_id)
        The function will be called by the the patch_object of the habitat belong with
        '''
        self.offspring_marker_pool = []
        offspring_num = self.hab_sexual_pairwise_parents_num() * sexual_birth_rate
        offspring_num_int = int(offspring_num)                    #整数部分
        offspring_num_dem = offspring_num - offspring_num_int     #小数部分
        
        offspring_marker_ls = [(patch_name, self.name, 'sexual') for _ in range(offspring_num_int)]
        if offspring_num_dem > np.random.uniform(0,1,1)[0]:
            one_offspring_marker_ls = [(patch_name, self.name, 'sexual')]
        else:
            one_offspring_marker_ls = []
        self.offspring_marker_pool = offspring_marker_ls + one_offspring_marker_ls
        return len(self.offspring_marker_pool)
    
    def hab_mix_reproduce_calculation_into_offspring_marker_pool(self, patch_name, asexual_birth_rate, sexual_birth_rate):
        ''' 
        Calculating the offspring num and generating offspring markers 
        but do not reproduce actually until local germination process afterward after dispersal process 
        patch_id: the patch of the habitat belong with, (i.e., patch_object.patch_id)
        The function will be called by the the patch_object of the habitat belong with
        '''
        self.offspring_marker_pool = []
        asex_offspring_num = len(self.asexual_parent_pos_ls) * asexual_birth_rate
        asex_offspring_num_int = int(asex_offspring_num)
        asex_offspring_num_dem = asex_offspring_num - asex_offspring_num_int
        
        asex_offspring_marker_ls = [(patch_name, self.name, 'mix_asexual') for _ in range(asex_offspring_num_int)]
        if asex_offspring_num_dem > np.random.uniform(0,1,1)[0]:
            asex_one_offspring_marker_ls = [(patch_name, self.name, 'mix_asexual')]
        else:
            asex_one_offspring_marker_ls = []
            
        sex_offspring_num = self.hab_mixed_sexual_pairwse_parents_num() * sexual_birth_rate
        sex_offspring_num_int = int(sex_offspring_num)
        sex_offspring_num_dem = sex_offspring_num - sex_offspring_num_int
        
        sex_offspring_marker_ls = [(patch_name, self.name, 'mix_sexual') for _ in range(sex_offspring_num_int)]
        if sex_offspring_num_dem > np.random.uniform(0,1,1)[0]:
            sex_one_offspring_marker_ls = [(patch_name, self.name, 'mix_sexual')]
        else:
            sex_one_offspring_marker_ls = []
        self.offspring_marker_pool = asex_offspring_marker_ls + asex_one_offspring_marker_ls + sex_offspring_marker_ls + sex_one_offspring_marker_ls
            
        return len(self.offspring_marker_pool)
        
#*************************** reprodction into offsprings pool processes *************************************
    def hab_asex_reproduce_mutate_into_offspring_pool(self, asexual_birth_rate):
        ''' NemoBackend: asexual birth into offspring_pool (expected offspring = indi_num x asexual_birth_rate; the fractional part is born by chance)
        INPUT: asexual_birth_rate
        RETURN: len(offspring_pool) '''
        self.hab_clear_up_offspring_pool()                                                                 # NemoBackend Method: recycle leftover offspring instead of dropping them
        offspring_num = self.indi_num * asexual_birth_rate
        num = int(offspring_num + np.random.uniform(0,1,1)[0])                                             # stochastic rounding: the fractional part is born by chance
        self.offspring_pool = self.hab_asex_reproduce_mutate_with_num(num)                                 # NemoBackend Method: return a list of offspring individual object
        return len(self.offspring_pool)
    
    def hab_sex_reproduce_mutate_into_offspring_pool(self, sexual_birth_rate):
        ''' NemoBackend: sexual birth into offspring_pool (expected offspring = parent_pairs x sexual_birth_rate; the fractional part is born by chance)
        INPUT: sexual_birth_rate
        RETURN: len(offspring_pool) '''
        self.hab_clear_up_offspring_pool()                          # NemoBackend Method: recycle leftover offspring instead of dropping them
        offspring_num = self.hab_sexual_pairwise_parents_num() * sexual_birth_rate
        num = int(offspring_num + np.random.uniform(0,1,1)[0])      # stochastic rounding: the fractional part is born by chance
        self.offspring_pool = self.hab_sex_reproduce_mutate_with_num(num) # NemoBackend Method: return a list of offspring individual object
        return len(self.offspring_pool)

    def hab_mix_reproduce_mutate_into_offspring_pool(self, asexual_birth_rate, sexual_birth_rate):
        ''' NemoBackend: mixed birth into offspring_pool (asexual: asexual parents x asexual_birth_rate; sexual: parent pairs x sexual_birth_rate; each fractional part is born by chance)
        INPUT: asexual_birth_rate, sexual_birth_rate
        RETURN: len(offspring_pool) '''
        self.hab_clear_up_offspring_pool()                                              # NemoBackend Method: recycle leftover offspring instead of dropping them
        asex_offspring_num = self.hab_mixed_asexual_parent_num() * asexual_birth_rate
        asex_num = int(asex_offspring_num + np.random.uniform(0,1,1)[0])                # stochastic rounding: the fractional part is born by chance
        sex_offspring_num = self.hab_mixed_sexual_pairwse_parents_num() * sexual_birth_rate
        sex_num = int(sex_offspring_num + np.random.uniform(0,1,1)[0])                  # stochastic rounding: the fractional part is born by chance
        self.offspring_pool = self.hab_mix_asex_reproduce_mutate_with_num(asex_num) + self.hab_mix_sex_reproduce_mutate_with_num(sex_num) # NemoBackend Methods
        return len(self.offspring_pool)

#*********************************** local germination processes ************************************************
    def _germinate_from_pools(self, pool_name_ls, is_remove=True):
        ''' NemoBackend: germinate individuals randomly chosen from the given pools into empty microsites
        INPUT: pool_name_ls = pool names in the original order, e.g. ('offspring_pool', 'immigrant_pool', 'dormancy_pool'); is_remove = True (default): germinated individuals leave their pool / False: a copy germinates, the original stays in its pool
        RETURN: number germinated '''
        hab_empty_pos_ls = self.empty_site_pos_ls
        hab_pools_with_labels = [(indi_object, pool_name) for pool_name in pool_name_ls for indi_object in getattr(self, pool_name)]

        random.shuffle(hab_empty_pos_ls)
        random.shuffle(hab_pools_with_labels)

        counter = 0
        for (row_id, col_id), (indi_object, pool_label) in list(zip(hab_empty_pos_ls, hab_pools_with_labels)):
            if is_remove:
                self.add_individual(indi_object=indi_object, len_id=row_id, wid_id=col_id)
                getattr(self, pool_label).remove(indi_object)
            else:
                self.add_individual(indi_object=indi_object.make_copy(), len_id=row_id, wid_id=col_id)    # NemoBackend Method: individual.make_copy() <--> Individual::operator=
            counter += 1
        return counter

    def hab_local_germinate_from_offspring_and_dormancy_pool(self, is_remove=True):
        ''' NemoBackend: germinate individuals from offspring_pool + dormancy_pool into empty microsites
        INPUT: is_remove = True (default): germinated individuals leave their pool / False: a copy germinates, the original stays in its pool
        RETURN: number germinated '''
        return self._germinate_from_pools(('offspring_pool', 'dormancy_pool'), is_remove)

    def hab_local_germinate_from_offspring_and_immigrant_pool(self, is_remove=True):
        ''' NemoBackend: same as hab_local_germinate_from_offspring_and_dormancy_pool(), from offspring_pool + immigrant_pool '''
        return self._germinate_from_pools(('offspring_pool', 'immigrant_pool'), is_remove)

    def hab_local_germinate_from_offspring_immigrant_and_dormancy_pool(self, is_remove=True):
        ''' NemoBackend: same as hab_local_germinate_from_offspring_and_dormancy_pool(), from offspring_pool + immigrant_pool + dormancy_pool '''
        return self._germinate_from_pools(('offspring_pool', 'immigrant_pool', 'dormancy_pool'), is_remove)
    
#************************** dormancy process in the habitat *********************************************************
    def _recycle_dropped(self, before_ls, after_ls):
        ''' NemoBackend: recycle every individual of before_ls that is not in after_ls (dropped by MetaIBM's pool bookkeeping)
        INPUT: before_ls = individuals in offspring / immigrant / dormancy pools before the bookkeeping; after_ls = individuals in those pools after it
        RETURN: number recycled '''
        after_set = set(after_ls)                                                       # identity-based lookup (individual has no __eq__ / __hash__)
        dropped_ls = [indi_obj for indi_obj in before_ls if indi_obj not in after_set]  # keep before_ls order; a duplicate shows up as a double recycle
        for indi_obj in dropped_ls:
            indi_obj.recycle()                                                          # NemoBackend Methods
        return len(dropped_ls)

    def hab_dormancy_process_from_offspring_pool_to_dormancy_pool(self):
        ''' NemoBackend: offspring enter the dormancy pool (capacity dormancy_pool_max_size); immigrants and everyone else dropped from the pools are recycled
        RETURN: (survival_in_dormancy_pool_num, eliminate_from_dormancy_pool_num, offspring_num, len(dormancy_pool)) '''
        before_ls = self.offspring_pool + self.immigrant_pool + self.dormancy_pool

        offspring_num = len(self.offspring_pool)
        dormancy_num = len(self.dormancy_pool)
        
        if self.dormancy_pool_max_size == 0:
            eliminate_from_dormancy_pool_num = dormancy_num  # 如果中途改变 dormancy_pool_max_size 的值
            survival_in_dormancy_pool_num = 0
            offspring_num = 0
            self.offspring_pool = []
            self.immigrant_pool = []
            self.dormancy_pool = []
        
        elif offspring_num > self.dormancy_pool_max_size:
            eliminate_from_dormancy_pool_num = dormancy_num   # when full: dormancy_num = self.dormancy_pool_max_size
            survival_in_dormancy_pool_num = 0
            offspring_num = self.dormancy_pool_max_size
            self.dormancy_pool = random.sample(self.offspring_pool, self.dormancy_pool_max_size)
            self.offspring_pool = []
            self.immigrant_pool = []
        
        elif offspring_num + dormancy_num <= self.dormancy_pool_max_size:
            eliminate_from_dormancy_pool_num = 0
            survival_in_dormancy_pool_num = dormancy_num
            self.dormancy_pool += self.offspring_pool
            self.offspring_pool = []
            self.immigrant_pool = []
        
        else:
            eliminate_from_dormancy_pool_num = offspring_num + dormancy_num - self.dormancy_pool_max_size
            survival_in_dormancy_pool_num = self.dormancy_pool_max_size - offspring_num
            survival_dormancy_pool = random.sample(self.dormancy_pool, survival_in_dormancy_pool_num)
            self.dormancy_pool = survival_dormancy_pool + self.offspring_pool
            self.offspring_pool = []
            self.immigrant_pool = []

        after_ls = self.offspring_pool + self.immigrant_pool + self.dormancy_pool
        self._recycle_dropped(before_ls, after_ls)              # NemoBackend Method
        return survival_in_dormancy_pool_num, eliminate_from_dormancy_pool_num, offspring_num, len(self.dormancy_pool)
    
    def hab_dormancy_process_from_offspring_pool_and_immigrant_pool(self):
        ''' NemoBackend: offspring and immigrants enter the dormancy pool (capacity dormancy_pool_max_size; priority offspring > immigrants > old dormants); everyone dropped is recycled
        RETURN: (survival_in_dormancy_pool_num, eliminate_from_dormancy_pool_num, new_dormancy_num, len(dormancy_pool)) '''
        before_ls = self.offspring_pool + self.immigrant_pool + self.dormancy_pool

        offspring_num = len(self.offspring_pool)
        immigrant_num = len(self.immigrant_pool)
        dormancy_num = len(self.dormancy_pool)
    
        if self.dormancy_pool_max_size < offspring_num:
            eliminate_from_dormancy_pool_num = dormancy_num     # when full: dormancy_num = self.dormancy_pool_max_size
            survival_in_dormancy_pool_num = 0
            new_dormancy_num = self.dormancy_pool_max_size
            self.dormancy_pool = random.sample(self.offspring_pool, self.dormancy_pool_max_size)
            self.offspring_pool = []
            self.immigrant_pool = []
            
        elif offspring_num <= self.dormancy_pool_max_size < (offspring_num + immigrant_num):
            eliminate_from_dormancy_pool_num = dormancy_num     # when full: dormancy_num = self.dormancy_pool_max_size
            survival_in_dormancy_pool_num = 0
            new_dormancy_num = self.dormancy_pool_max_size
            self.dormancy_pool = self.offspring_pool + random.sample(self.immigrant_pool, self.dormancy_pool_max_size-offspring_num)
            self.offspring_pool = []
            self.immigrant_pool = []
        
        elif (offspring_num + immigrant_num) <= self.dormancy_pool_max_size < (offspring_num + immigrant_num + dormancy_num):
            eliminate_from_dormancy_pool_num = offspring_num + immigrant_num + dormancy_num - self.dormancy_pool_max_size   # when full: dormancy_num = self.dormancy_pool_max_size
            survival_in_dormancy_pool_num = dormancy_num - eliminate_from_dormancy_pool_num                                 # when full: dormancy_num = self.dormancy_pool_max_size
            new_dormancy_num = offspring_num + immigrant_num
            self.dormancy_pool = self.offspring_pool + self.immigrant_pool + random.sample(self.dormancy_pool, survival_in_dormancy_pool_num)
            self.offspring_pool = []
            self.immigrant_pool = []
            
        else:
            eliminate_from_dormancy_pool_num = 0
            survival_in_dormancy_pool_num = dormancy_num
            new_dormancy_num = offspring_num + immigrant_num
            self.dormancy_pool += self.immigrant_pool + self.offspring_pool
            self.offspring_pool = []
            self.immigrant_pool = []

        after_ls = self.offspring_pool + self.immigrant_pool + self.dormancy_pool
        self._recycle_dropped(before_ls, after_ls)                   # NemoBackend Method
        return survival_in_dormancy_pool_num, eliminate_from_dormancy_pool_num, new_dormancy_num, len(self.dormancy_pool)
    
# ******************** clearing up offspring_and_immigrant_pool when dormancy pool do not run *******************************************************************#
    def hab_clear_up_microsite_individuals(self):
        ''' NemoBackend: recycle every individual in the microsites and reset the microsite bookkeeping (pools untouched)
        RETURN: number recycled '''
        dead_ls = [indi_obj for row in self.set['microsite_individuals'] for indi_obj in row if indi_obj is not None]

        self.set['microsite_individuals'] = self._empty_microsite_grid()
        self.empty_site_pos_ls = [(i, j) for i in range(self.length) for j in range(self.width)]
        self.occupied_site_pos_ls = []
        self.indi_num = 0
        self.species_category = {}
        self.asexual_parent_pos_ls = []
        self.species_category_for_sexual_parents_pos = {}

        for indi_obj in dead_ls:
            indi_obj.recycle()                                                          # NemoBackend Methods
        return len(dead_ls)

    def hab_clear_up_offspring_pool(self):
        ''' NemoBackend: recycle every ungerminated offspring, then empty offspring_pool
        RETURN: number recycled '''
        ungerminated_ls = self.offspring_pool
        self.offspring_pool = []
        for indi_obj in ungerminated_ls:
            indi_obj.recycle()
        return len(ungerminated_ls)

    def hab_clear_up_immigrant_pool(self):
        ''' NemoBackend: recycle every ungerminated immigrant, then empty immigrant_pool
        RETURN: number recycled '''
        ungerminated_ls = self.immigrant_pool
        self.immigrant_pool = []
        for indi_obj in ungerminated_ls:
            indi_obj.recycle()
        return len(ungerminated_ls)

    def hab_clear_up_dormancy_pool(self):
        ''' NemoBackend: recycle every dormant individual, then empty dormancy_pool -- teardown only, e.g. releasing a finished run's Nemo memory before the next run in the same process
        (disturbance keeps the dormancy pool; never call it inside the dormancy processes, where _recycle_dropped already recycles)
        RETURN: number recycled '''
        dormant_ls = self.dormancy_pool
        self.dormancy_pool = []
        for indi_obj in dormant_ls:
            indi_obj.recycle()
        return len(dormant_ls)

    def hab_clear_up_offspring_marker_pool(self):
        ''' empty offspring_marker_pool (markers are tuples: nothing to recycle)
        RETURN: number of markers cleared '''
        num = len(self.offspring_marker_pool)
        self.offspring_marker_pool = []
        return num

    def hab_clear_up_immigrant_marker_pool(self):
        ''' empty immigrant_marker_pool (markers are tuples: nothing to recycle)
        RETURN: number of markers cleared '''
        num = len(self.immigrant_marker_pool)
        self.immigrant_marker_pool = []
        return num

    def habitat_disturbance_process(self):
        ''' NemoBackend: disturbance kills every individual in the microsites and in the offspring / immigrant pools (dormancy_pool survives); marker pools are emptied
        RETURN: number of individuals recycled + number of markers cleared '''
        return (self.hab_clear_up_microsite_individuals() + self.hab_clear_up_offspring_pool() + self.hab_clear_up_immigrant_pool() + self.hab_clear_up_offspring_marker_pool() + self.hab_clear_up_immigrant_marker_pool())