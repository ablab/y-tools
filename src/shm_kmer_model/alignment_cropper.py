#!/usr/bin/env python2

from __future__ import print_function
from abc import ABCMeta, abstractmethod


class AbstractAlignmentCropper:
    __metaclass__ = ABCMeta

    @abstractmethod
    def crop(self, alignment):
        pass


class UptoLastReliableKMerAlignmentCropper(AbstractAlignmentCropper):
    def crop(self, alignment):
        """ Crop the tail """
        start_last_reliable_k_mer = alignment.read.__len__() - self.k_mer_len
        while alignment.read.seq[start_last_reliable_k_mer:
                                 start_last_reliable_k_mer + self.k_mer_len] != \
                alignment.germline_seq.seq[start_last_reliable_k_mer:
                                           start_last_reliable_k_mer + self.k_mer_len] and \
                start_last_reliable_k_mer >= 0:
            start_last_reliable_k_mer -= 1

        alignment.read.seq = alignment.read.seq[:start_last_reliable_k_mer + self.k_mer_len]
        alignment.germline_seq.seq = alignment.germline_seq.seq[:start_last_reliable_k_mer +
                                                                self.k_mer_len]

        """ Crop the head """
        start_first_reliable_k_mer = 0
        while alignment.read.seq[start_first_reliable_k_mer:
                                 start_first_reliable_k_mer + self.k_mer_len] != \
                alignment.germline_seq.seq[start_first_reliable_k_mer:
                                           start_first_reliable_k_mer + self.k_mer_len] and \
                start_first_reliable_k_mer + self.k_mer_len <= len(alignment.read.seq):
            start_first_reliable_k_mer += 1

        alignment.read.seq = alignment.read.seq[start_first_reliable_k_mer:]
        alignment.germline_seq.seq = alignment.germline_seq.seq[start_first_reliable_k_mer:]

        return alignment

    def __init__(self, config):
        self.k_mer_len = config.k_mer_len
