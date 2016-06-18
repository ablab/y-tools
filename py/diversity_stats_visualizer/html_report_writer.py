import os
import sys
import warnings

import pandas as pd

import visualize_shm_stats
import visualize_cdr_stats


class HTMLReportWriter:
    def _WriteHeader(self, level, text, align):
        self.fhandler.write("<h" + str(level) + " align = " + align + ">" + text + "</h" + str(level) + ">\n")

    def _GetAlignment(self, align_to_center = True):
        if not align_to_center:
            return "left"
        return "center"

    def _WriteTableCell(self, elem):
        self.fhandler.write("<td>" + str(elem) + "</td>")

    def _WriteTableRow(self, row):
        self.fhandler.write("<tr>\n")
        for r in row:
            self._WriteTableCell(r)
        self.fhandler.write("</tr>\n")

    def __init__(self, output_fname):
        self.fhandler = open(output_fname, "w")

    def WriteH1(self, text, align_to_center = True):
        self._WriteHeader(1, text, self._GetAlignment(align_to_center))

    def WriteH2(self, text, align_to_center = True):
        self._WriteHeader(2, text, self._GetAlignment(align_to_center))

    def WriteH3(self, text, align_to_center = True):
        self._WriteHeader(3, text, self._GetAlignment(align_to_center))

    # width in percent
    def WriteImage(self, path_to_image, width = 60):
        self.fhandler.write("<p align = center>\n")
        self.fhandler.write("<image src = " + path_to_image + " width = " + str(width) + "%></image>\n")
        self.fhandler.write("</p>\n")

    def WriteTable(self, col_names, row_names, values, width = 60):
        if len(values) == 0:
            return
        if len(row_names) != len(values):
            print "# rows in table and # row names are not consistent"
            sys.exit(1)
        if len(col_names) != len(values[0]):
            print "# columns in table and # column names are not consistent"
            sys.exit(1)
        self.fhandler.write("<table width = " + str(width) + "% align = center>\n")
        header_row = [""]
        header_row.extend(col_names)
        self._WriteTableRow(header_row)
        for i in range(0, len(row_names)):
            cur_row = [row_names[i]]
            cur_row.extend(values[i])
            self._WriteTableRow(cur_row)
        self.fhandler.write("</table>\n")

    def WriteEmptyLine(self):
        self.fhandler.write("<br>\n")

    def WriteHorizontalLine(self, width = 100):
        self.fhandler.write("<hr width=" + str(width) + "%>\n")

    def CloseFile(self):
        self.fhandler.close()

#######################################################################################################################
def ComputeLocusCharacteristics(vj_df, locus):
    stats = []
    locus_df = vj_df.loc[vj_df['Chain_type'] == locus]
    stats.append(len(locus_df["Read_name"]))
    stats.append(len([s for s in locus_df['Productive'] if s == 1]))
    stats.append(len([s for s in locus_df['Has_stop_codon'] if s == 1]))
    stats.append(len([s for s in locus_df['In-frame'] if s == 0]))
    return stats

def ComputeGeneralCharacteristicsTable(vj_df):
    col_names = ["IGH", "IGK", "IGL"]
    row_names = ["# sequences", "# productive sequences", "# sequences with stop codon", "# out-of-frame sequences"]
    table = []
    for r in row_names:
        table.append([0] * len(col_names))
    for i in range(0, len(col_names)):
        locus_stats = ComputeLocusCharacteristics(vj_df, col_names[i])
        for j in range(0, len(locus_stats)):
            table[j][i] = locus_stats[j]
    return table, row_names, col_names

def WriteGeneralCharacteristics(html_writer, vj_df, images_dict):
    html_writer.WriteH1("General characteristics")
    table, row_names, col_names = ComputeGeneralCharacteristicsTable(vj_df)
    html_writer.WriteTable(col_names, row_names, table)
    if "vj_heatmap" in images_dict:
        html_writer.WriteH2("Heatmap of VJ hit concentrations")
        html_writer.WriteImage(images_dict["vj_heatmap"], width=70)

#######################################################################################################################
def ComputeGeneralSHMCharacteristics(shm_df):
    col_names = ['V segment', 'J segment']
    row_names = ['# SHMs', '# substitutions', '# deletions', '# insertions', '# synonymous SHMs', '# stop codon SHMs']
    row_index_dict = dict()
    for i in range(0, len(row_names)):
        row_index_dict[row_names[i]] = i
    table = []
    for r in row_names:
        table.append([0] * len(col_names))
    for read in shm_df:
        shms = shm_df[read]
        index_j = 0
        if not read.is_variable():
            index_j = 1
        for shm in shms:
            table[row_index_dict["# SHMs"]][index_j] += 1
            if shm.is_substitution():
                table[row_index_dict["# substitutions"]][index_j] += 1
            if shm.is_deletion():
                table[row_index_dict["# deletions"]][index_j] += 1
            if shm.is_insertion():
                table[row_index_dict["# insertions"]][index_j] += 1
            if shm.synonymous:
                table[row_index_dict["# synonymous SHMs"]][index_j] += 1
            if shm.to_stop_codon:
                table[row_index_dict["# stop codon SHMs"]][index_j] += 1
    return row_names, col_names, table

def WriteSHMCharacteristics(html_writer, shm_df, images_dict):
    html_writer.WriteH1("SHM characteristics")
    row_names, col_names, table = ComputeGeneralSHMCharacteristics(shm_df)
    html_writer.WriteTable(col_names, row_names, table)
    if "aa_substitutions" in images_dict:
        html_writer.WriteH2("Heatmap of amino acid substitutions:")
        html_writer.WriteImage(images_dict["aa_substitutions"], width = 60)
    if "general_v_mutations" in images_dict:
        html_writer.WriteH2("Distribution of SHM positions in V gene segment (above) and "
                            "distribution of SHM number in V gene segment per sequence (below):")
        html_writer.WriteImage(images_dict["general_v_mutations"], width = 80)
    if "synonymous_shms" in images_dict:
        html_writer.WriteH2("Distribution of synonymous SHM positions in V gene segment:")
        html_writer.WriteImage(images_dict["synonymous_shms"], width = 60)
    if "special_shms" in images_dict:
        html_writer.WriteH2("Distribution of special SHM (stop codon, deletions, insertions) positions in V gene segment:")
        html_writer.WriteImage(images_dict["special_shms"], width = 70)

#######################################################################################################################
def ComputeLocusCDRCharacteristics(vj_df, locus):
    locus_df = vj_df.loc[vj_df['Chain_type'] == locus]
    nums = list()
    cdrs = ['CDR1', 'CDR2', 'CDR3']
    for cdr in cdrs:
        region_seq = list(locus_df[cdr + "_nucls"])
        seq_dict = dict()
        for seq in region_seq:
            if seq not in seq_dict:
                seq_dict[seq] = 0
            seq_dict[seq] += 1
        max_abun = 0
        for seq in seq_dict:
            max_abun = max(max_abun, seq_dict[seq])
        nums.extend([len(seq_dict), max_abun])
    return nums

def ComputeGeneralCDRCharacteristicsTable(vj_df):
    col_names = ["IGH", "IGK", "IGL"]
    row_names = ["# CDR1s", "max CDR1 abundance",
                 "# CDR2s", "max CDR2 abundance",
                 "# CDR3s", "max CDR3 abundance"]
    table = []
    for r in row_names:
        table.append([0] * len(col_names))
    for i in range(0, len(col_names)):
        nums = ComputeLocusCDRCharacteristics(vj_df, col_names[i])
        for j in range(0, len(nums)):
            table[j][i] = nums[j]
    return table, row_names, col_names

def WriteGeneralCDRCharacteristics(html_writer, vj_df):
    table, row_names, col_names = ComputeGeneralCDRCharacteristicsTable(vj_df)
    html_writer.WriteTable(col_names, row_names, table)

def GetLengthAbundanceLargestGroup(locus_df, locus, cdr):
    region_seq = list(locus_df[cdr + "_nucls"])
    max_group = visualize_cdr_stats.get_region_largest_group(region_seq)
    if len(max_group) == 0:
        return 0, 0
    return len(max_group[0]), round(float(len(max_group)) / float(len(region_seq)) * 100, 2)

def WriteCDRPlots(html_writer, vj_df, images_dict):
    loci = ['IGH', 'IGK', 'IGL']
    cdrs = ['CDR1', 'CDR2', 'CDR3']
    for l in loci:
        locus_df = vj_df.loc[vj_df['Chain_type'] == l]
        for cdr in cdrs:
            if l + "_" + cdr + "_length" in images_dict:
                html_writer.WriteH2(l + " " + cdr + " length distribution:")
                html_writer.WriteImage(images_dict[l + "_" + cdr + "_length"], 65)
            else:
                continue
            max_len, max_abun = GetLengthAbundanceLargestGroup(locus_df, l, cdr)
            if l + "_" + cdr + "_nucls" in images_dict:
                html_writer.WriteH2("Most " + l + " " + cdr + "s (" + str(max_abun) + "%) have length " + str(max_len) + " nt")
                html_writer.WriteH3("Distribution of nucleotide abundance per position for " + l + " " + cdr +
                                    "s of length " + str(max_len) + " nt:")
                html_writer.WriteImage(images_dict[l + "_" + cdr + "_nucls"], 70)
                if l + "_" + cdr + "_aa" in images_dict:
                    html_writer.WriteH3(l + " " + cdr + "s of length " + str(max_len) + " nt are in-frame. "
                                                                                   "Plot of the most abundant amino "
                                                                                   "acids per position:")
                    html_writer.WriteImage(images_dict[l + "_" + cdr + "_aa"], 70)
        html_writer.WriteHorizontalLine(70)

def WriteCDRCharacteristics(html_writer, vj_df, images_dict):
    html_writer.WriteH1("CDR characteristics")
    WriteGeneralCDRCharacteristics(html_writer, vj_df)
    WriteCDRPlots(html_writer, vj_df, images_dict)

#######################################################################################################################

def create_html(output_html_fname, vj_df, shm_df, images_dict):
    print "Annotation report will be written to " + output_html_fname
    print "Printing general characteristics of the repertoire"
    html_writer = HTMLReportWriter(output_html_fname)
    WriteGeneralCharacteristics(html_writer, vj_df, images_dict)
    html_writer.WriteHorizontalLine()
    print "Printing SHM characteristics"
    WriteSHMCharacteristics(html_writer, shm_df, images_dict)
    html_writer.WriteHorizontalLine()
    print "Printing CDR characteristics"
    WriteCDRCharacteristics(html_writer, vj_df, images_dict)
    html_writer.CloseFile()
    print "Annotation report was written to " + output_html_fname

#######################################################################################################################
def add_cdr_plots(image_dict, plots_dir):
    inner_plots_dir = "plots/cdr_plots"
    if not os.path.exists(os.path.join(plots_dir, "cdr_plots")):
        return
    cdr_dict = os.path.join(plots_dir, "cdr_plots")
    loci = ['IGH', 'IGK', 'IGL']
    cdrs = ['CDR1', 'CDR2', 'CDR3']
    for l in loci:
        for cdr in cdrs:
            file_dict = {l + "_" + cdr + "_length" : l + "_" + cdr + "_length.png",
                         l + "_" + cdr + "_nucls" : l + "_" + cdr + "_nucls.png",
                         l + "_" + cdr + "_aa" : l + "_" + cdr + "_aa.png"}
            for f in file_dict:
                if os.path.exists(os.path.join(cdr_dict, file_dict[f])):
                    image_dict[f] = os.path.join(inner_plots_dir, file_dict[f])

def add_shm_plots(image_dict, plots_dir):
    inner_plots_dir = "plots"
    if os.path.exists(os.path.join(plots_dir, "aa_substitutions.png")):
        image_dict["aa_substitutions"] = os.path.join(inner_plots_dir, "aa_substitutions.png")
    if os.path.exists(os.path.join(plots_dir, "v_mutations_distribution.png")):
        image_dict["general_v_mutations"] = os.path.join(inner_plots_dir, "v_mutations_distribution.png")
    if os.path.exists(os.path.join(plots_dir, "synonymous_shms_positions.png")):
        image_dict["synonymous_shms"] = os.path.join(inner_plots_dir, "synonymous_shms_positions.png")
    if os.path.exists(os.path.join(plots_dir, "special_shms_positions.png")):
        image_dict["special_shms"] = os.path.join(inner_plots_dir, "special_shms_positions.png")


def create_image_dict(plots_dir):
    image_dict = dict()
    if os.path.exists(os.path.join(plots_dir, "vj_heatmap.png")):
        image_dict["vj_heatmap"] = os.path.join("plots", "vj_heatmap.png")
    add_cdr_plots(image_dict, plots_dir)
    add_shm_plots(image_dict, plots_dir)
    return image_dict

#######################################################################################################################

def main(argv):
    warnings.filterwarnings('ignore')
    if len(argv) != 5:
        print "Invalid input"
        print "python html_report_writer.py cdr_details.txt shm_details.txt plots_dir output.html"
        sys.exit(1)
    vj_df_fname = argv[1]
    if not os.path.exists(vj_df_fname):
        print "CDR details file " + vj_df_fname + " was not found"
        sys.exit(1)
    vj_df = pd.read_table(vj_df_fname, delim_whitespace = True)
    shm_df_fname = argv[2]
    if not os.path.exists(shm_df_fname):
        print "SHM details file " + shm_df_fname + " was not found"
        sys.exit(1)
    shm_df = visualize_shm_stats.SHMs(shm_df_fname)
    plots_dir = argv[3]
    if not os.path.exists(plots_dir):
        print "Plot directory " + plots_dir + " does not exist"
        sys.exit(1)
    image_dict = create_image_dict(plots_dir)
    output_html_fname = argv[4]
    create_html(output_html_fname, vj_df, shm_df, image_dict)

if __name__ == "__main__":
    main(sys.argv)

