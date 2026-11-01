# -*- coding: utf-8 -*-
''' indfactory -- one species = one boot sim = one trait prototype = one Nemo IndFactory. The only module of MetaIBM-Nemo that imports pynemo.
'''

import numpy as np
from pynemo import nemo
from pynemo.sim import NemoSim

from .individual import individual

class indfactory():
    """ One species = one boot sim = one trait prototype = one Nemo IndFactory """

    ''' Class Attributes '''
    _nemo_id_probe = nemo.Individual()         # access to getcurrentID()/setCurrentID() on Nemo's static Individual::currentID; created once at import (takes one ID)

    def __init__(self, species_id, pheno_name_ls, nemo_params, trait_keys=('quant',), seed=None, metaibm_nemo_params=None):

        ''' MetaIBM Attributes '''
        self.species_id = species_id
        self.pheno_name_ls = pheno_name_ls      # quanti-trait axis name
        self.traits_num = len(pheno_name_ls)    # quanti-trait axis num

        ''' Nemo Attributes '''
        self.trait_keys = trait_keys
        self.nemo_params = nemo_params          # pynemo.params.SimParameters: record only, do not edit after boot
        self.nemo_params_dict = dict(nemo_params._params)    # CACHE of nemo_params: e.g., {'quanti_traits': '2', 'quanti_allele_model': 'diallelic'}
        
        self.traits_index = {}                  # CACHE of IndFactory::getTraitIndex(), avoiding frequent SWIG call; e.g., {'ntrl': 0, 'quant': 1}, sorted by type string
        self._nemo_simulator = None             # pynemo.sim.NemoSim (Python wrapper): boot fixture, never stepped
        self._nemo_indfactory = None            # pynemo.nemo.IndFactory (C++ objects)

        ''' MetaIBM-Nemo Attributes '''
        self.metaibm_nemo_params = dict(metaibm_nemo_params or {})   # dict: Nemo parameters that MetaIBM re-implements itself (Nemo never sees them)
        self._quanti_seq_pos = None                                  # CACHE: list of int, sequence position of quanti locus k (quanti_init_freq order, trait by trait)
        
        ''' Individual Counters (for leak checks) '''
        self.minted_num = 0                     # MetaIBM-Nemo individuals created so far
        self.recycled_num = 0                   # MetaIBM-Nemo individuals recycled so far

        self._boot(seed)

    ##***************************** display ***********************************
    def __str__(self):
        metaibm_str = 'species_id=%s, pheno_names=%s, minted=%d, recycled=%d' % (str(self.species_id), str(tuple(self.pheno_name_ls)), self.minted_num, self.recycled_num)
        nemo_str = 'traits_index=%s, seed=%d' % (str(self.traits_index), self._nemo_simulator._seed)
        return 'indfactory\n' + '  metaibm-info : %s\n' % metaibm_str + '  nemo-info    : %s' % nemo_str

    def __repr__(self):
        return '<indfactory species_id=%r traits=%r>' % (self.species_id, tuple(sorted(self.traits_index)))
    
    ##********** Nemo's static Individual::currentID Interface *************
    @classmethod
    def get_current_individual_id(cls):
        ''' read Nemo's global individual ID counter -- process-wide, shared by all species
            RETURN: int -- the ID the next Nemo individual will get '''
        return cls._nemo_id_probe.getcurrentID()

    @classmethod
    def set_current_individual_id(cls, value):
        ''' set Nemo's global individual ID counter -- makePrototype() zeroes it, so _boot() restores it
            INPUT:  value = int -- the ID the next Nemo individual will get '''
        cls._nemo_id_probe.setCurrentID(value)

    ##***************************** boot ***********************************
    def _check_quanti_traits_num(self):
        ''' (quant) traits_num must equal the number of quanti trait axes (0 without 'quant'), otherwise get_phenotype(i) reads past the last axis '''
        quanti_traits = int(self.nemo_params_dict.get('quanti_traits', 0)) if 'quant' in self.trait_keys else 0
        if self.traits_num != quanti_traits or ('quant' in self.trait_keys and quanti_traits < 1): raise ValueError('species %r: trait_keys=%r, %d pheno names, %d quanti trait axes' % (self.species_id, self.trait_keys, self.traits_num, quanti_traits))

    def _check_trait_prototypes_built(self):
        ''' trait_keys must equal the traits Nemo built from nemo_params -- Nemo silently skips a trait whose mandatory parameters are missing '''
        built_traits = {trait for trait in ('delet', 'dmi', 'fdisp', 'mdisp', 'ntrl', 'quant', 'wolb') if self._nemo_simulator._pop.getTraitPrototype(trait) is not None}
        if built_traits != set(self.trait_keys): raise ValueError('species %r: trait_keys=%r, but nemo_params built %r' % (self.species_id, self.trait_keys, sorted(built_traits)))

    def _boot(self, seed):
        ''' boot this species' sim and build its IndFactory from the prototypes in trait_keys; save the static Individual::currentID first, restore it last -- makePrototype() zeroes it, and colliding IDs reuse crossovers
            INPUT:  seed = int or None -- None: random_seed in nemo_params, else a time-based seed '''
        self._check_quanti_traits_num()
        saved_id = self.get_current_individual_id()
        try:
            self._nemo_simulator = NemoSim(self.nemo_params, seed=seed)   # params -> .ini -> parsed records; resolve the seed; empty Metapop + SimRunner -- no prototypes, RAND not seeded yet
            self._nemo_simulator._init()                                  # seed RAND -> build trait prototypes (+ LCEs) -> write log file, stays at gen 0; NEVER step_through(): _teardown() -> clearPrototype()

            self._check_trait_prototypes_built()
            proto_traits = {trait_key: self._nemo_simulator._pop.getTraitPrototype(trait_key) for trait_key in self.trait_keys}     # IndFactory::_protoTraits

            self._nemo_indfactory = nemo.IndFactory()
            self._nemo_indfactory.makePrototype(proto_traits)             # zeroes the static Individual::currentID; builds the proto-individual + _TraitsIndex (sorted by type)
            self.traits_index = {trait_key: self._nemo_indfactory.getTraitIndex(trait_key) for trait_key in self.trait_keys}
        finally:
            self.set_current_individual_id(saved_id)                      # restore the static Individual::currentID back to the saved_id (even if error was raise above)

    ##**************************** birth ************************************
    @staticmethod
    def metaibm2nemo_sex(gender):
        ''' convert MetaIBM's gender into Nemo's sex_t
            INPUT:  gender = 'female' or 'male'
            RETURN: nemo.FEM or nemo.MAL '''
        if gender not in ('female', 'male'): raise ValueError("gender must be 'female' or 'male', got %r" % (gender,))
        return nemo.FEM if gender == 'female' else nemo.MAL       # FEM=1, MAL=0 (int)

    def _wrap_nemo_individual(self, nemo_indi_obj):
        ''' wrap a newborn Nemo individual into a MetaIBM individual of this species, and count it in minted_num
            INPUT:  nemo_indi_obj = pynemo.nemo.Individual -- fresh from this IndFactory
            RETURN: individual '''
        self.minted_num += 1      # counter of MetaIBM-Nemo individuals created so far
        return individual(species_id=self.species_id, pheno_names_ls=self.pheno_name_ls, nemo_indi_obj=nemo_indi_obj, indfactory=self)

    def _require_own(self, indi_obj):
        ''' check the individual is alive and minted by this indfactory -- Nemo checks neither (parents from another IndFactory may carry another trait layout)
            INPUT:  indi_obj = individual '''
        indi_obj._require_alive()
        if indi_obj.indfactory is not self: raise RuntimeError('individual of species %r does not belong to the indfactory of species %r' % (indi_obj.species_id, self.species_id))

    def _use_quanti_init_freq(self):
        ''' whether make_founder() redraws the quanti alleles from metaibm_nemo_params['quanti_init_freq']; RETURN: bool '''
        if self.metaibm_nemo_params.get('quanti_init_freq') is None: return False
        if self.nemo_params_dict.get('quanti_allele_model') in ('continuous', 'continuous_in_place', 'continuous_HC') or self.nemo_params_dict.get('quanti_pleiotropy') in ('full', 'f', 'var', 'variable', 'v'):
            raise ValueError('species %r: quanti_init_freq requires a di-allelic model without pleiotropy' % (self.species_id,))
        return self.nemo_params_dict.get('quanti_allele_model') in ('diallelic', 'diallelic_HC', 'diallelic_asymmetrical') 

    def _set_quanti_init_alleles(self, nemo_indi_obj):
        ''' redraw every quanti allele of a fresh founder on its own: allele 0 with P = metaibm_nemo_params['quanti_init_freq'][k], else allele 1 (TTQuanti::set_allele_bit); then recompute the phenotype
            INPUT:  nemo_indi_obj = pynemo.nemo.Individual -- fresh from create_first_gen(), not wrapped yet '''
        nemo_trait = nemo_indi_obj.getTrait(self.traits_index['quant'])
        if self._quanti_seq_pos is None:
            self._quanti_seq_pos = [nemo_trait.get_locus_seq_pos(locus, trait_dim) for trait_dim in range(self.traits_num) for locus in range(nemo_trait.get_num_loci(trait_dim))]
        is_allele_1 = np.random.random((2, len(self._quanti_seq_pos))) >= np.asarray(self.metaibm_nemo_params['quanti_init_freq'], dtype=float)      # (strand, locus) bool
        for strand in (0, 1):
            for k, pos in enumerate(self._quanti_seq_pos):
                nemo_trait.set_allele_bit(pos, strand, bool(is_allele_1[strand, k]))
        nemo_trait.set_value()                                                     # phenotype = new genotype + a new Ve draw

    def make_founder(self, gender, home_habitat=0):
        ''' one founder with random genotypes drawn from this species' nemo_params (e.g. quanti_init_value)
            INPUT:  gender = 'female' or 'male'; home_habitat = int 0-65535 -- habitat.global_index of its birthplace (as Nemo's homepatch)
            RETURN: individual -- parent IDs are sentinels counting down from ULONG_MAX '''
        nemo_indi_obj = self._nemo_indfactory.makeNewIndividual(None, None, self.metaibm2nemo_sex(gender), home_habitat)   # no parents
        nemo_indi_obj.create_first_gen()
        if self._use_quanti_init_freq(): self._set_quanti_init_alleles(nemo_indi_obj)
        return self._wrap_nemo_individual(nemo_indi_obj)

    def make_clone(self, mother_indi_obj, home_habitat=0):
        ''' one asexual offspring: a copy of the mother that only mutates (no inheritance, no recombination) -- Nemo's makeClonedOffsprg (= LCE_Breed_base::breed_cloning)
            INPUT:  mother_indi_obj = individual of this species, either gender (the clone keeps it); home_habitat = int 0-65535 -- habitat.global_index where the clone is born (as Nemo's homepatch)
            RETURN: individual -- new ID, age 0, mother_id = the mother's ID, father_id = 0, pedigree class 4; alleles = the mother's + Nemo mutations; phenotype recomputed by Nemo '''
        self._require_own(mother_indi_obj)
        return self._wrap_nemo_individual(self._nemo_indfactory.makeClonedOffsprg(mother_indi_obj.nemo_indi_obj, home_habitat))

    def make_offspring(self, mother_indi_obj, father_indi_obj, gender, home_habitat=0):
        ''' one sexual offspring: inheritance, recombination and mutation all done in C++ by Nemo's makeOffsprg() -- replaces MetaIBM's deepcopy + recombine + mutate
            INPUT:  mother_indi_obj, father_indi_obj = individual of this species; gender = 'female' or 'male' -- the offspring's; home_habitat = int 0-65535 -- habitat.global_index where the offspring is born (as Nemo's homepatch)
            RETURN: individual '''
        self._require_own(mother_indi_obj)
        self._require_own(father_indi_obj)
        return self._wrap_nemo_individual(self._nemo_indfactory.makeOffsprg(mother_indi_obj.nemo_indi_obj, father_indi_obj.nemo_indi_obj, self.metaibm2nemo_sex(gender), home_habitat))

    ##***************************** life recycle ***********************************
    def recycle(self, indi_obj):
        ''' return an individual's nemo_indi_obj to Nemo's recycling-pool and mark it dead -- call it after MetaIBM's bookkeeping (del_individual reads gender)
            INPUT:  indi_obj = individual of this species '''
        self._require_own(indi_obj)
        nemo_indi_obj = indi_obj.nemo_indi_obj
        indi_obj.is_alive = False
        indi_obj.nemo_indi_obj = None
        self._nemo_indfactory.recycle(nemo_indi_obj)
        self.recycled_num += 1

    def purge_recycling_pool(self):
        ''' free the C++ memory of every individual in the recycling-pool '''
        self._nemo_indfactory.purgeRecyclingPOOL()
    

    
