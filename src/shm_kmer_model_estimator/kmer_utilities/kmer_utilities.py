import itertools


def nucl_bases():
    return ['A', 'C', 'G', 'T']


def kmer_names(kmer_len=5):
    """ Simple routine to get all kmer_names over DNA alphabet. """
    bases = nucl_bases()
    return [''.join(p) for p in itertools.product(bases, repeat=kmer_len)]


def central_nucl_indexes(kmer_len=5):
    """ Simple routine to get indexes in the alphabet
    of a central nucl in each kmer. """
    bases = nucl_bases()
    n_nucl = len(bases)
    half_kmer_len = kmer_len // 2
    indexes = [(x % n_nucl**(half_kmer_len + 1)) //
               n_nucl**half_kmer_len
               for x in xrange(len(nucl_bases())**kmer_len)]
    return indexes


def nonmutated_index(kmer, kmer_len=5):
    """ Simple routine to get index in the alphabet
    of a central nucl in a kmer. """
    nonmutated_nucl = kmer[kmer_len // 2]
    return nucl_bases().index(nonmutated_nucl)
