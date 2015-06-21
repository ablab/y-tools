def lfactorial(x):
    """
    Compute log-factorial using log-gamma function
    """
    from math import lgamma

    return lgamma(x + 1)


def anonnamedtuple(field_names, *args, **kwargs):
    """
    Make anonymous named tuple (like make_tuple in C++11)
    """
    from collections import namedtuple

    return namedtuple("_", field_names)(*args, **kwargs)


def groupby_dict(iterable, key=None, lazy=False, key_as_index=False,
                 raw_dict=False):
    """
    Analogue of itertools.groupby, but using dict instead of sorted list
    """
    from collections import defaultdict

    if key is None:
        key = lambda x: x

    d = defaultdict(list)
    if not key_as_index:
        for e in iterable:
            d[key(e)].append(e)
    else:
        for k, e in zip(key, iterable):
            d[k].append(e)

    if raw_dict:
        return d
    else:
        return d.iteritems() if lazy else list(d.iteritems())


def count_multiplicity(reads):
    gd = groupby_dict(reads, reads, key_as_index=True)
    reads = []
    multiplicity = []
    for read, l in gd:
        reads.append(read)
        multiplicity.append(len(l))

    return reads, multiplicity


def hamming_dist(s1, s2):
    """
    Ordinary Hamming distance between two strings
    """
    s1, s2 = s1.strip(), s2.strip()
    assert(len(s1) == len(s2))
    return sum([c1 != c2 for c1, c2 in zip(s1, s2)])

def levenshtein_dist(s1, s2):
    import Levenshtein
    s1, s2 = s1.strip(), s2.strip()
    return Levenshtein.distance(s1, s2)


def hamming_matrix(reads):
    import numpy as np

    N = len(reads)
    m = np.zeros((N, N), dtype=int)

    for i in range(N):
        for j in range(i):
            dist = hamming_dist(reads[i], reads[j])
            m[i, j] = m[j, i] = dist

    return m


def levenshtein_matrix(reads):
    import numpy as np

    N = len(reads)
    m = np.zeros((N, N), dtype=int)

    for i in range(N):
        for j in range(i):
            dist = levenshtein_dist(reads[i], reads[j])
            m[i, j] = m[j, i] = dist

    return m


def hamming_graph_knuth(reads, tau=1, **kwargs):
    """
    Construct hamming(tau) graph using Knuth's algorithm
    """
    import igraph as ig
    from collections import defaultdict
    from Levenshtein import hamming

    l = len(reads[0])

    for read in reads:
        assert(len(read) == l)

    piece_len = int(l / (tau + 1))
    piece_len_last = l - piece_len * tau

    piece_lens = [piece_len] * tau + [piece_len_last]

    g = ig.Graph(len(reads))
    g.vs["reads"] = reads

    for attr_name, attr_data in kwargs.iteritems():
        g.vs[attr_name] = attr_data

    for j in range(tau + 1):
        sets = defaultdict(list)
        for i in range(len(reads)):
            substr = reads[i][j*piece_len:j*piece_len + piece_lens[j]]
            sets[substr].append(i)

        for v_list in sets.itervalues():
            for i1 in range(len(v_list)):
                for i2 in range(i1 + 1, len(v_list)):
                    v1, v2 = v_list[i1], v_list[i2]
                    read1, read2 = reads[v1], reads[v2]
                    d = hamming(read1, read2)
                    if d <= tau and g.get_eid(v1, v2, error=False) == -1:
                        g.add_edge(v1, v2, weight=d)

    return g


def hamming_graph_naive(reads, tau=1, **kwargs):
    """
    Construct hamming(tau) graph using naive O(N**2 d) algorithm
    """
    import igraph as ig
    import numpy as np
    from Levenshtein import hamming

    l = len(reads[0])

    for read in reads:
        assert(len(read) == l)

    N = len(reads)
    m = np.zeros((N, N), dtype=int)

    for i in range(N):
        for j in range(i):
            dist = hamming(reads[i], reads[j])
            m[i, j] = m[j, i] = dist if dist <= tau else 0

    # Be careful! Zero elements are not interpreted as zero-length edges
    g = ig.Graph.Weighted_Adjacency(m.tolist(),
                                    mode="UNDIRECTED",
                                    attr="weight",
                                    loops=False)

    g.vs["read"] = reads

    for attr_name, attr_data in kwargs.iteritems():
        g.vs[attr_name] = attr_data

    return g


def hamming_graph(reads, tau=1, multiplicity=None, use_naive=False, **kwargs):
    """
    Hamming graph(tau) construction wrapper
    """
    l = len(reads[0])

    for read in reads:
        assert(len(read) == l)

    if multiplicity is None:
        multiplicity = [1] * len(reads)

    if use_naive or tau + 1 > l:
        g = hamming_graph_naive(reads, tau=tau,
                                multiplicity=multiplicity, **kwargs)
    else:
        g = hamming_graph_knuth(reads, tau=tau,
                                multiplicity=multiplicity, **kwargs)

    return g


def levenshtein_graph(reads, tau=1, multiplicity=None):
    """
    Construct Levenshtein(tau) graph using naive O(N**2 d) algorithm
    """
    import igraph as ig
    import numpy as np

    N = len(reads)
    m = np.zeros((N, N), dtype=int)

    for i in range(N):
        for j in range(i):
            dist = levenshtein_dist(reads[i], reads[j])
            m[i, j] = m[j, i] = dist if dist <= tau else 0

    # Be careful! Zero elements are not interpreted as zero-length edges
    g = ig.Graph.Weighted_Adjacency(m.tolist(),
                                    mode="UNDIRECTED",
                                    attr="weight",
                                    loops=False)

    g.vs["read"] = reads

    if multiplicity is None:
        multiplicity = [1] * N
    g.vs["multiplicity"] = multiplicity

    return g


def findDuplicates(l):
    """
    Find and count duplicates in list using dict
    """
    from collections import defaultdict, namedtuple

    element_count = defaultdict(int)
    for e in l:
        element_count[e] += 1

    uniques, counts = [list(_) for _ in zip(*element_count.iteritems())]

    return namedtuple("elements_counts",
                      ["elements", "counts"])(uniques, counts)


def count_triplet(s, t, align=0):
    return sum([s[i:(i+3)] == t for i in range(align, len(s) - 3 + 1, 3)])


def stop_codons_count(s, align=0):
    """
    Stop codons (from wiki):
        - TAG ("amber")
        - TAA ("ochre")
        - TGA ("opal" or "umber")
    """
    cs = [count_triplet(s, t, align=align) for t in ["TAG", "TAA", "TGA"]]
    return sum(cs)


def optimal_align_using_stop_codons(reads):
    import numpy as np
    return np.argmin([sum([stop_codons_count(read, align=align)
                           for read in reads])
                      for align in [0, 1, 2]])


def stop_codon_indices(s):
    res = []
    for i in range(len(s)-3+1):
        if s[i:i+3] in ["TAG", "TAA", "TGA"]:
            res.append(i)
    return res

def stop_codon_indices_all(reads):
    res = []
    for read in reads:
        res.extend(stop_codon_indices(read))

    return res
