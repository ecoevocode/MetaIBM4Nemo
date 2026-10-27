# -*- coding: utf-8 -*-
''' individual -- a MetaIBM-Nemo individual whose genetics live in a Nemo `Individual`. 

ATTENTION:
1. Do NOT give this class an __eq__ / __hash__ -- keep identity comparison (the core relies on list.remove and `in`).
2. Never keep nemo_indi_obj or a trait from get_trait() outside this class: after recycle() Nemo reuses that memory for a newborn.
3. MetaIBM bookkeeping first, recycle() last: del_individual() reads gender, which is unreadable after recycle().
'''

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

    ##***************** Nemo-only interface ******************
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

    def genotype_table(self, trait_dim=0):
        ''' QUANTI ONLY: allele values on both chromosome copies at each locus of one quanti trait axis
            INPUT:  trait_dim -- index of the trait axis, 0 ... traits_num-1
            RETURN: list of (strand0, strand1) -- one tuple (of allele values) per locus '''
        trait = self.get_trait('quant')
        positions = [trait.get_locus_seq_pos(locus, trait_dim) for locus in range(trait.get_num_loci(trait_dim))] # get loci_num of this axis, then their positions in the allele sequence
        return [(trait.get_allele_value(pos, 0), trait.get_allele_value(pos, 1)) for pos in positions]

    def allele_table(self, trait_key, num_loci):
        ''' allele identities on both chromosome copies at each locus of given trait_key (for quanti use genotype_table)
            INPUT:  trait_key = 'ntrl', 'delet' or 'dmi'; num_loci = int -- must equal ntrl_loci, delet_loci or dmi_loci
            RETURN: list of (strand0, strand1) -- one tuple (of allele identities) per locus '''
        trait = self.get_trait(trait_key)
        return [(trait.get_allele(locus, 0), trait.get_allele(locus, 1)) for locus in range(num_loci)]

    ##*************** Nemo-MetaIBM interface ***************
    @property
    def gender(self):
        ''' RETURN: 'female' or 'male' -- read from Nemo, fixed at birth '''
        self._require_alive()
        return 'female' if self.nemo_indi_obj.isFemale() else 'male'

    @property
    def phenotype_set(self):
        ''' read-only, read live from Nemo; the i-th quanti trait axis is named pheno_names_ls[i]
            RETURN: {pheno_name: phenotype} -- writing to it raises TypeError '''
        trait = self.get_trait('quant')
        return MappingProxyType({pheno_name: trait.get_phenotype(i) for i, pheno_name in enumerate(self.pheno_names_ls)})

    def get_indi_phenotype_ls(self):
        ''' phenotypes in trait order -- what selection reads
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
