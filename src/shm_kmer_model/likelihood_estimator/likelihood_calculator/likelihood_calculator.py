#!/usr/bin/env python2

from __future__ import print_function

import numpy as np
import pandas as pd

import mismatch_finder.mismatch_finder as mismatch_finder
import shm_kmer_model.shm_kmer_model as shm_kmer_model

class LikehoodCalculator(object):
    def __init__(self, model, kmer_len=5, model_mode=shm_kmer_model.ModelMode.Both,
            mismatch_finder=mismatch_finder.NoKNeighboursMismatchFinder(5)):
        self.model = model
        self.kmer_len = kmer_len
        self.half_kmer_len = self.kmer_len // 2
        self.mismatch_finder = mismatch_finder
        self.model_mode = model_mode

    def __calculate_subst_likelihood(self, source_fivemer, mutated_nucl):
        return np.log(self.model.expect_subst_prob(source_fivemer, mutated_nucl))

    def __calculate_mutation_likelihood(self, source_fivemer):
        return np.log(self.model.expect_mut_prob(source_fivemer))

    def calculate_likelihood(self, source, destination):
        assert len(source) == len(destination)
        mismatch_positions = self.mismatch_finder.find_mismatch_positions(source, destination)
        likelihood = 0
        for mismatch_position in mismatch_positions:
            source_fivemer = source[(mismatch_position - self.half_kmer_len) : \
                                    (mismatch_position + self.half_kmer_len + 1)]
            if self.model_mode != shm_kmer_model.ModelMode.Mutation:
                likelihood += self.__calculate_subst_likelihood(source_fivemer,
                              destination[mismatch_position])
            if self.model_mode != shm_kmer_model.ModelMode.Substitution:
                likelihood += self.__calculate_mutation_likelihood(source_fivemer)
        return likelihood
