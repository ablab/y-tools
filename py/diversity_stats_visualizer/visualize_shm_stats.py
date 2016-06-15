import os
import sys
import operator
import warnings

import matplotlib as mplt
mplt.use('Agg')

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class AlignedRead:
    def _parse_line(self, line):
        splits = line.split()
        if len(splits) != 6:
            print "Header line " + line + " is not valid"
            sys.exit(1)
        self.read_name = splits[0][len("Read_name:"):]
        self.read_len = int(splits[1][len("Read_length:"):])
        self.gene_name = splits[2][len("Gene_name:"):]
        self.gene_len = int(splits[3][len("Gene_length:"):])
        self.segment = splits[4][len("Segment:"):]
        self.chain_type = splits[5][len("Chain_type:"):]

    def __init__(self, line):
        self._parse_line(line)

    def __hash__(self):
        return hash(self.read_name) * hash(self.gene_name)

    def __eq__(self, other):
        return self.read_name == other.read_name and self.gene_name == other.gene_name

    def is_variable(self):
        return self.segment == 'V'

class SHM:
    def __init__(self, line):
        splits = line.split()
        if len(splits) != 9:
            print "Invalid SHM line " + line
        self.type = splits[0]
        self.read_pos = int(splits[1])
        self.gene_pos = int(splits[2])
        self.read_nucl = splits[3]
        self.gene_nucl = splits[4]
        self.read_aa = splits[5]
        self.gene_aa = splits[6]
        self.synonymous = False
        if int(splits[7]) == 1:
            self.synonymous = True
        self.to_stop_codon = False
        if int(splits[8]) == 1:
            self.to_stop_codon = True

    def is_deletion(self):
        return self.type == 'D'

    def is_insertion(self):
        return self.type == 'I'

class SHMs:
    def __init__(self, df_fname):
        self.shm_dict = dict()
        fhandler = open(df_fname, "r")
        lines = fhandler.readlines()
        current_read = ""
        for i in range(1, len(lines)):
            l = lines[i].strip()
            if l[:len("Read_name:")] == "Read_name:":
                current_read = AlignedRead(l)
                self.shm_dict[current_read] = list()
            else:
                self.shm_dict[current_read].append(SHM(l))
        print str(len(self.shm_dict)) + " records were extracted from " + df_fname

    def __len__(self):
        return self.shm_dict

    def __getitem__(self, item):
        return self.shm_dict[item]

    def __iter__(self):
        for it in self.shm_dict:
            yield it

def output_shms_pos(all_shms_pos, colors):
    pos = []
    labels = []
    cols = []
    for isotype in all_shms_pos:
        if len(all_shms_pos[isotype]) > 0:
            pos.append(all_shms_pos[isotype])
            labels.append(str(isotype))
            cols.append(colors[isotype])
    plt.hist(pos, bins= 30, color = cols, alpha = .75, label = labels)
    plt.legend(loc = 'upper center', ncol = len(pos), fontsize = 16)
    plt.xlabel("Relative position of SHM in V gene segment", fontsize = 16)
    plt.ylabel("# SHMs", fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)

def output_num_shms(num_all_shms, colors):
    nums = []
    cols = []
    for isotype in num_all_shms:
        if len(num_all_shms[isotype]) > 0:
            sns.distplot(num_all_shms[isotype], hist = False, label = str(isotype), color = colors[str(isotype)])
            #nums.append(num_shms[isotype])
            #labels.append(str(isotype))
            #cols.append(colors[isotype])
    #plt.hist(nums, bins= 25, color = cols, alpha = .5, label = labels)
    #plt.legend(loc = 'upper center', ncol = len(nums))
    plt.xlabel("#SHM in V gene segment", fontsize = 16)
    plt.ylabel("# sequences", fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlim(0, 150)
    plt.legend(loc = 'upper center', ncol = len(num_all_shms), fontsize = 16)

def visualize_v_mutations_stats(shms_df, output_fname):
    all_shms_pos = {'IGH': [], 'IGK': [], 'IGL': []}
    num_all_shms = {'IGH': [], 'IGK': [], 'IGL': []}
    colors = {'IGH': 'b', 'IGK': 'g', 'IGL': 'r'}
    for it in shms_df:
        if not it.is_variable():
            continue
        read_shms = shms_df[it]
        for shm in read_shms:
            all_shms_pos[it.chain_type].append(float(shm.read_pos) / float(it.gene_len))
        num_all_shms[it.chain_type].append(len(shms_df[it]))
    plt.figure(1, figsize=(9, 6))
    plt.subplot(211)
    output_shms_pos(all_shms_pos, colors)
    plt.subplot(212)
    output_num_shms(num_all_shms, colors)
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()
    print "Distribution of SHMs in V was written to " + output_fname

def get_aa_list():
    return ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def aa_is_valid(aa):
    return aa != '*' and aa != '-' and aa != 'X'

def visualize_aa_substitution_matrix(shms_df, output_fname):
    dict_aa = dict()
    for it in shms_df:
        read_shms = shms_df[it]
        prev_pos = -1
        for shm in read_shms:
            if prev_pos / 3 != shm.read_pos / 3:
                if aa_is_valid(shm.gene_aa) and aa_is_valid(shm.read_aa):
                    aa_pair = shm.gene_aa + shm.read_aa
                    if not aa_pair in dict_aa:
                        dict_aa[aa_pair] = 0
                    dict_aa[aa_pair] += 1
            prev_pos = shm.read_pos
    aa_list = get_aa_list()
    aa_freq = []
    for i in range(0, len(aa_list)):
        aa_freq.append([0] * len(aa_list))
    for aa_pair in dict_aa:
        aa_freq[aa_list.index(aa_pair[1])][aa_list.index(aa_pair[0])] = dict_aa[aa_pair]
    fig, ax = plt.subplots()
    sns.heatmap(aa_freq, cmap = plt.cm.jet, xticklabels = aa_list, yticklabels = aa_list, square = True, ax = ax)
    ax.tick_params(labelsize = 14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12, rotation='horizontal')
    plt.xlabel("From", fontsize = 14)
    plt.ylabel("To", fontsize = 14, rotation='horizontal')
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()
    print "Amino acid substitution heatmap was written to " + output_fname

def visualize_special_shm_positions(shm_df, output_fname):
    synonymous_pos = []
    stop_codon_pos = []
    deletion_pos = []
    insertion_pos = []
    for it in shm_df:
        read_shms = shm_df[it]
        for shm in read_shms:
            relative_pos = float(shm.read_pos) / float(it.read_len)
            if shm.synonymous:
                synonymous_pos.append(relative_pos)
            elif shm.to_stop_codon:
                stop_codon_pos.append(relative_pos)
            elif shm.is_deletion():
                deletion_pos.append(relative_pos)
            elif shm.is_insertion():
                insertion_pos.append(relative_pos)
    pos = []
    labels = []
    colors = []
    plt.figure(figsize=(12, 9))
    if len(synonymous_pos) > 100:
        pos.append(synonymous_pos)
        labels.append('Synonymous')
        colors.append('r')
        #sns.distplot(synonymous_pos, hist = False, label = "Synonymous SHMs", color = 'r')
    if len(stop_codon_pos) > 100:
        pos.append(stop_codon_pos)
        labels.append('Stop codon')
        colors.append('g')
        #sns.distplot(stop_codon_pos, hist = False, label = "Stop codon SHMs", color = 'g')
    if len(deletion_pos) > 100:
        pos.append(deletion_pos)
        labels.append('Deletions')
        colors.append('b')
        #sns.distplot(deletion_pos, hist = False, label = "Deletion SHMs", color = 'b')
    if len(insertion_pos) > 100:
        pos.append(insertion_pos)
        labels.append('Insertions')
        colors.append('orange')
        #sns.distplot(insertion_pos, hist = False, label = "Insertion SHMs", color = 'orange')
    plt.hist(pos, color = colors, label= labels, bins = 30)
    plt.xlim(0, 1)
    plt.legend(loc = 'upper center', ncol = len(pos), fontsize = 12, bbox_to_anchor=(0.5, -0.07))
    plt.xlabel("Relative position on read", fontsize = 14)
    plt.ylabel("# SHMs", fontsize = 14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()
    print "Distribution of special SHMs was written to " + output_fname

def main(shm_df_fname, output_dir):
    print "== Output SHMs statistics"
    shm_df = SHMs(shm_df_fname)
    visualize_v_mutations_stats(shm_df, os.path.join(output_dir, "v_mutations_distribution.pdf"))
    visualize_aa_substitution_matrix(shm_df, os.path.join(output_dir, "aa_substitutions.pdf"))
    visualize_special_shm_positions(shm_df, os.path.join(output_dir, "special_shms_positions.pdf"))

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Invalid input parameters"
        print "python visualize_shm_stats.py shm_df.txt output_dir"
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])