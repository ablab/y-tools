#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys

from itertools import izip

from ig_report_supernodes import smart_open, parse_size

if __name__ == "__main__":
    parser = ArgumentParser(description="Construct repertoire from TrieCompressor output")
    parser.add_argument("--input-reads", "-i",
                        type=str,
                        help="file with aligned RepSeq sequences")
    parser.add_argument("--input-map", "-m",
                        type=str,
                        help="file with map reported by TrieCompressor")
    parser.add_argument("--input-compressed", "-c",
                        type=str,
                        help="file with compressed reads reported by TrieCompressor")
    parser.add_argument("--output-repertoire", "-r",
                        type=str,
                        help="output file with repertoire sequences")
    parser.add_argument("--output-rcm", "-R",
                        type=str,
                        help="output file with repertoire RCM")

    args = parser.parse_args()

    print "Construct repertoire from TrieCompressor output..."
    print "Command line: %s" % " ".join(sys.argv)

    # Fix ids
    with smart_open(args.input_compressed) as fin, smart_open(args.output_repertoire, "w") as fout:
        for i, record in enumerate(SeqIO.parse(fin, "fasta")):
            id = record.description
            size = parse_size(id)
            record.id = record.description = "cluster___%d___size___%d" % (i, size)
            SeqIO.write(record, fout, "fasta")

    with smart_open(args.input_reads) as fin_reads, smart_open(args.input_map) as fin_map, smart_open(args.output_rcm, "w") as fout_rcm:
        for record, cluster in izip(SeqIO.parse(fin_reads, "fasta"), fin_map):
            id = record.description
            cluster = cluster.strip()
            fout_rcm.write("%s\t%s\n" % (id, cluster))

    print "Done"
