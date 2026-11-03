# -*- coding: utf-8 -*-
''' individual -- a MetaIBM-Nemo individual whose genetics live in a Nemo `Individual`. 

ATTENTION:
1. Do NOT give this class an __eq__ / __hash__ -- keep identity comparison (the core relies on list.remove and `in`).
2. Never keep nemo_indi_obj or a trait from get_trait() outside this class: after recycle() Nemo reuses that memory for a newborn.
3. MetaIBM bookkeeping first, recycle() last: del_individual() reads gender, which is unreadable after recycle().
'''
import numpy as np
from types import MappingProxyType

class individual():
    def __init__(self, species_id, pheno_names_ls, nemo_indi_obj, indfactory):

        ''' MetaIBM Attributes '''
        self.species_id = species_id
        self.pheno_names_ls = pheno_names_ls
        self.traits_num = len(pheno_names_ls)
        self.age = 0

        ''' Nemo Attributes '''
        self.nemo_indi_obj = nemo_indi_obj         # pynemo.nemo.Individual (C++ objects)
        self.indfactory = indfactory               # the indfactory that manages it
        self.is_alive = True

    def _require_alive(self):
        ''' check whether the nemo_indi_obj is already in the Nemo's recycling-pool '''
        if not self.is_alive: raise RuntimeError('individual (species_id=%r) was already recycled: it can no longer be read, bred from or recycled again' % (self.species_id,))

    ##******************* display ***********************
    def __str__(self):
        metaibm_str = 'species_id=%s, age=%d, pheno_names=%s' % (str(self.species_id), self.age, str(tuple(self.pheno_names_ls)))
        if not self.is_alive: return 'individual (recycled)\n' + '  metaibm-info : %s\n' % metaibm_str + '  nemo-info    : no data after recycle()'
        pheno_str = ', '.join('%.4f' % value for value in self.get_indi_phenotype_ls()) if self.traits_num else 'none'   # no 'quant' -> no phenotype to read
        pedigree_str = 'founder' if self.mother_id >= 2**63 else 'clone of %d' % self.mother_id if self.father_id == 0 else 'mother %d, father %d' % (self.mother_id, self.father_id)
        return ('individual\n' + '  metaibm-info : %s\n' % metaibm_str + '  nemo-info    : id=%d, gender=%s, phenotype=(%s), %s' % (self.individual_id, self.gender, pheno_str, pedigree_str))

    def __repr__(self):
        if not self.is_alive: return '<individual species_id=%r recycled>' % (self.species_id,)
        return '<individual species_id=%r id=%d>' % (self.species_id, self.individual_id)

    ##***************** Nemo‘s interface ******************
    @property
    def individual_id(self):
        ''' RETURN: Nemo's ID of this individual object: assigned at birth, unique across one simulation '''
        self._require_alive()
        return self.nemo_indi_obj.getID()

    @property
    def mother_id(self):
        ''' RETURN: an int copied at birth, valid even after the mother is recycled (never use getMother()); founders: sentinels near ULONG_MAX; clones: the parent '''
        self._require_alive()
        return self.nemo_indi_obj.getMotherID()

    @property
    def father_id(self):
        ''' RETURN: an int copied at birth, valid even after the father is recycled (never use getFather()); founders: sentinels near ULONG_MAX; clones: 0 '''
        self._require_alive()
        return self.nemo_indi_obj.getFatherID()

    def get_trait(self, trait_key='quant'):
        ''' get one Nemo trait of this individual; slots are sorted by type string, so the slot is looked up via the indfactory
            INPUT:  trait_key = 'quant', 'ntrl', 'delet', 'dmi', 'wolb', 'fdisp' or 'mdisp' -- Nemo trait type
            RETURN: pynemo.nemo.TTrait -- a live view, not a copy; never keep it after recycle() '''
        self._require_alive()
        return self.nemo_indi_obj.getTrait(self.indfactory.traits_index[trait_key])

    # - - - - - - - - trait values readers - - - - - - - - #
    def get_trait_value(self, trait_key, **kwargs):
        ''' RETURN: value of one trait as stored by Nemo, in its own type; kwargs: trait_dim=i for 'quant' '''
        if trait_key == 'quant':                              return self.get_trait('quant').get_phenotype(kwargs.pop('trait_dim', 0), **kwargs)   # float: P of axis trait_dim
        elif trait_key in ('delet', 'dmi', 'fdisp', 'mdisp'): return self.get_trait(trait_key).get_value(**kwargs)        # float
        elif trait_key == 'wolb':                             return self.get_trait(trait_key).get_wolb_value(**kwargs)   # bool
        elif trait_key == 'ntrl':                             return self.get_trait(trait_key).getValue(**kwargs)         # None: Nemo stores no value
        else: raise ValueError("get_trait_value(%r) unsupported; trait_key = 'quant'/'ntrl'/'delet'/'dmi'/'fdisp'/'mdisp'/'wolb'; quant takes trait_dim=i" % trait_key)

    # - - - - - - - - allele readers - - - - - - - - #
    def get_trait_alleles(self, trait_key, method, **kwargs):
        ''' RETURN: alleles of one trait via TTrait::get_allele() or get_allele_value(); kwargs: trait_dim=i for 'quant' '''
        if trait_key == 'quant' and method == 'get_allele_value':              return self._get_quanti_allele_values_array(**kwargs)            # np.ndarray float,     shape=(n_loci, 2)
        elif trait_key == 'quant' and method == 'get_allele':                  return self._get_quanti_allele_states_array(**kwargs)            # np.ndarray int (0/1), shape=(n_loci, 2)
        elif trait_key in ('ntrl', 'delet', 'dmi') and method == 'get_allele': return self._get_allele_states_array(trait_key, **kwargs)        # np.ndarray int,       shape=(n_loci, 2)
        elif trait_key in ('fdisp', 'mdisp') and method == 'get_allele_value': return self._get_disp_allele_values_array(trait_key, **kwargs)   # np.ndarray float,     shape=(1, 2)
        elif trait_key == 'wolb' and method == 'get_allele':                   return self._get_wolb_infection_state(**kwargs)                  # bool value,           shape=()
        else: raise ValueError("get_trait_alleles(%r, %r) unsupported; 'get_allele_value': quant/fdisp/mdisp; 'get_allele': quant/ntrl/delet/dmi/wolb; quant takes trait_dim=i" % (trait_key, method))

    def _get_quanti_allele_values_array(self, trait_dim=0):
        ''' allele values at each locus of quanti trait axis trait_dim (0 ... traits_num-1)
            RETURN: np.ndarray of float, shape (n_loci, 2) -- row = locus, column = strand (0 paternal, 1 maternal) '''
        trait = self.get_trait('quant')
        positions = [trait.get_locus_seq_pos(locus, trait_dim) for locus in range(trait.get_num_loci(trait_dim))]   # loci of this axis -> sequence positions
        return np.array([(trait.get_allele_value(pos, 0), trait.get_allele_value(pos, 1)) for pos in positions], dtype=float).reshape(-1, 2)

    def _get_quanti_allele_states_array(self, trait_dim=0):
        ''' allelic states (0/1) at each locus of quanti trait axis trait_dim; di-allelic only (continuous: all 0)
            RETURN: np.ndarray of int, shape (n_loci, 2) -- row = locus, column = strand (0 paternal, 1 maternal) '''
        trait = self.get_trait('quant')
        positions = [trait.get_locus_seq_pos(locus, trait_dim) for locus in range(trait.get_num_loci(trait_dim))]   # loci of this axis -> sequence positions
        return np.array([(trait.get_allele(pos, 0), trait.get_allele(pos, 1)) for pos in positions], dtype=int).reshape(-1, 2)

    def _get_allele_states_array(self, trait_key):
        ''' allelic states (0-based) at each locus of trait_key = 'ntrl', 'delet' or 'dmi' (not 'quant')
            RETURN: np.ndarray of int, shape (n_loci, 2) -- row = locus, column = strand (0 paternal, 1 maternal) '''
        trait = self.get_trait(trait_key)
        num_loci = self.indfactory._nemo_indfactory.getTraitPrototype(trait_key).get_locus_number()                 # one sequence: locus == position
        return np.array([(trait.get_allele(locus, 0), trait.get_allele(locus, 1)) for locus in range(num_loci)], dtype=int).reshape(-1, 2)

    def _get_disp_allele_values_array(self, trait_key='fdisp'):
        ''' allele values (dispersal rates, 0 ... 1) at the single locus of trait_key = 'fdisp' or 'mdisp' (every individual carries both)
            RETURN: np.ndarray of float, shape (1, 2) -- column = strand (0 maternal, 1 paternal: reversed vs the other traits) '''
        trait = self.get_trait(trait_key)
        return np.array([(trait.get_allele_value(0, 0), trait.get_allele_value(0, 1))], dtype=float)

    def _get_wolb_infection_state(self):
        ''' RETURN: bool -- Wolbachia infection (no loci), inherited from the mother; set by Nemo's inoculation LCE, which MetaIBM does not run '''
        return bool(self.get_trait('wolb').get_allele(0, 0))

    ##*************** MetaIBM legacy interface (rebuilt by reading Nemo) ***************
    @property
    def gender(self):
        ''' RETURN: 'female' or 'male' -- read from Nemo, fixed at birth '''
        self._require_alive()
        return 'female' if self.nemo_indi_obj.isFemale() else 'male'

    @property
    def phenotype_set(self):
        ''' ONLY 'quant': read live from Nemo; the i-th quanti trait axis is named pheno_names_ls[i]
            RETURN: {pheno_name: phenotype} -- writing to it raises TypeError '''
        trait = self.get_trait('quant')
        return MappingProxyType({pheno_name: trait.get_phenotype(i) for i, pheno_name in enumerate(self.pheno_names_ls)})

    def get_indi_phenotype_ls(self):
        ''' ONLY 'quant': phenotypes in trait order -- what selection reads
            RETURN: list of float, length traits_num '''
        trait = self.get_trait('quant')
        return [trait.get_phenotype(i) for i in range(self.traits_num)]

    ##***************** life cycle events ******************
    def make_copy(self):
        ''' RETURN: an exact copy of this individual with a new ID -- calling Nemo's Individual::operator= (pynemo Individual.make_copy) '''
        self._require_alive()
        nemo_copy_indi_obj = self.indfactory._nemo_indfactory.getNewIndividual().make_copy(self.nemo_indi_obj)  # Individual::operator=
        copy_indi_obj = self.indfactory._wrap_nemo_individual(nemo_copy_indi_obj)                               # minted_num += 1
        copy_indi_obj.age = self.age
        return copy_indi_obj

    def recycle(self):
        ''' return the nemo_indi_obj to the recycling-pool of its indfactory '''
        self._require_alive()
        self.indfactory.recycle(self)
