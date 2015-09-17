#!/usr/bin/env python

import os
import shutil
import sys
import logging

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
spades_src = os.path.join(home_directory, "src/python_pipeline/")

import dense_subgraph_finder

sys.path.append(spades_src)
import support

#######################################################################################
#           Error messages
#######################################################################################

def ErrorMessagePrepareCfg(log):
    log.info("Probably you forgot to prepare IgRepertoireConstructor. Please follow instructions:")
    log.info("  (1) remove build/ directory")
    log.info("  (2) type command \'./prepare cfg\' to check all dependencies")
    log.info("  (3) type \'make\' to compile IgRepertoireConstrictor")
    log.info("  (4) rerun IgRepertoireConstructor")

def SupportInfo(log):
    log.info("\nIn case you have troubles running IgRepertoireConstructor, "
             "you can write to igtools_support@googlegroups.com.")
    log.info("Please provide us with ig_repertoire_constructor.log file from the output directory.")

#######################################################################################
#           Binary routines
#######################################################################################
class PhaseNames:
    def __init__(self):
        self.__vj_alignment = 'vj_alignment'
        self.__trie_compressor = 'trie_compressor'
        self.__graph_construction = 'graph_constructor'
        self.__dsf = 'dsf'
        self.__consensus_constructor = 'consensus_constructor'
        self.__remove_low_abundance_reads = 'remove_low_abundance_reads'
        self.__phase_order = [self.__vj_alignment,
                              self.__trie_compressor,
                              self.__graph_construction,
                              self.__dsf,
                              self.__consensus_constructor,
                              self.__remove_low_abundance_reads]
        self.__long_names = {'vj_alignment' : 'VJ Alignment',
                             'trie_compressor' : 'Trie Compressor',
                             'graph_constructor' : 'Graph Constructor',
                             'dsf' : 'Dense Subgraph Finder',
                             'consensus_constructor' : 'Consensus Constructor',
                             'remove_low_abundance_reads' : 'Low Abundant Clusters Remover'}

    def __iter__(self):
        for sname in self.__phase_order:
            yield sname

    def __len__(self):
        return len(self.__phase_order)

    def GetPhaseNameBy(self, index):
        return self.__phase_order[index]

    def GetPhaseIndex(self, phase_name):
        for i in range(len(self)):
            if self.GetPhaseNameBy(i) == phase_name:
                return i
        return -1

    def PhaseIsVJAlignment(self, phase_name):
        return phase_name == self.__vj_alignment

    def GetVJAlignmentLongName(self):
        return self.__long_names[self.__vj_alignment]

    def PhaseIsTrieCompressor(self, phase_name):
        return phase_name == self.__trie_compressor

    def GetTrieCompressorLongName(self):
        return self.__long_names[self.__trie_compressor]

    def PhaseIsGraphConstructor(self, phase_name):
        return phase_name == self.__graph_construction

    def GetGraphConstructionLongName(self):
        return self.__long_names[self.__graph_construction]

    def PhaseIsDSF(self, phase_name):
        return phase_name == self.__dsf

    def GetDSFLongName(self):
        return self.__long_names[self.__dsf]

    def PhaseIsConsensusConstructor(self, phase_name):
        return phase_name == self.__consensus_constructor

    def GetConsensusConstructorLongName(self):
        return self.__long_names[self.__consensus_constructor]

    def PhaseIsRemoveLowAbundanceReads(self, phase_name):
        return phase_name == self.__remove_low_abundance_reads

    def GetRemoveLowAbundanceReadsName(self):
        return self.__long_names[self.__remove_low_abundance_reads]

###########
class IgRepConConfig:
    def __initBinaryPaths(self):
        self.path_to_vj_aligner = os.path.join(home_directory, 'build/release/bin/ig_kplus_vj_finder')
        self.run_vj_aligner = os.path.join(home_directory, 'build/release/bin/./ig_kplus_vj_finder')
        self.path_to_trie_compressor = os.path.join(home_directory, 'build/release/bin/ig_trie_compressor')
        self.run_trie_compressor = os.path.join(home_directory, 'build/release/bin/./ig_trie_compressor')
        self.path_to_graph_constructor = os.path.join(home_directory, 'build/release/bin/ig_swgraph_construct')
        self.run_graph_constructor = os.path.join(home_directory, 'build/release/bin/./ig_swgraph_construct')
        self.path_to_consensus_constructor = os.path.join(home_directory, 'build/release/bin/ig_consensus_finder')
        self.run_consensus_constructor = os.path.join(home_directory, 'build/release/bin/./ig_consensus_finder')
        self.run_rcm_recoverer = os.path.join(home_directory, 'src/ig_quast_tool/rcm_recoverer.py')
        self.run_remove_low_abundance_reads = os.path.join(home_directory, 'src/ig_quast_tool/ig_remove_low_abundance_reads.py')
        self.run_report_supernodes = os.path.join(home_directory, 'src/ig_quast_tool/ig_report_supernodes.py')
        self.path_to_dsf = os.path.join(home_directory, 'build/release/bin/dense_sgraph_finder')
        self.path_to_germline = os.path.join(home_directory, "build/release/bin/germline")

    def __init__(self):
        self.__initBinaryPaths()
        self.cluster_size_limit = 5

    def CheckBinaries(self, log):
        phase_names = PhaseNames()
        if not os.path.exists(self.path_to_vj_aligner):
            log.info("ERROR: Binary file of " + phase_names.GetVJAlignmentLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_trie_compressor):
            log.info("ERROR: Binary file of " + phase_names.GetTrieCompressorLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_graph_constructor):
            log.info("ERROR: Binary file of " + phase_names.GetGraphConstructionLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_dsf):
            log.info("ERROR: Binary file of " + phase_names.GetDSFLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_consensus_constructor) or not os.path.exists(self.run_rcm_recoverer):
            log.info("ERROR: Binary file of " + phase_names.GetConsensusConstructorLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.run_remove_low_abundance_reads):
            log.info("ERROR: Binary file of " + phase_names.GetRemoveLowAbundanceReadsName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)

class IgRepConIO:
    def __initVJFinderOutput(self, output_dir):
        self.vj_finder_output = os.path.join(output_dir, "vj_finder")
        self.cropped_reads = os.path.join(self.vj_finder_output, "cropped.fa")
        self.bad_reads = os.path.join(self.vj_finder_output, "bad.fa")
        self.vj_alignment_info = os.path.join(self.vj_finder_output, "add_info.csv")

    def __initCompressorOutput(self, output_dir):
        self.compressed_reads = os.path.join(output_dir, "compressed.fa")
        self.map_file = os.path.join(output_dir, "map.txt")
        self.supernodes_file = os.path.join(output_dir, "supernodes.fa")

    def __initDSFOutput(self, output_dir):
        self.dsf_output = os.path.join(output_dir, "dense_sgraph_finder")
        self.dense_sgraph_decomposition = os.path.join(self.dsf_output, 'dense_subgraphs.txt')

    def __initFinalOutput(self, output_dir):
        self.final_clusters_fa = os.path.join(output_dir, 'final_repertoire.fa')
        self.final_rcm = os.path.join(output_dir, 'final_repertoire.rcm')

    def __init__(self, output_dir, log):
        self.__log = log
        self.__initVJFinderOutput(output_dir)
        self.__initCompressorOutput(output_dir)
        self.sw_graph = os.path.join(output_dir, "sw.graph")
        self.__initDSFOutput(output_dir)
        self.__initFinalOutput(output_dir)
        self.final_stripped_clusters_fa = os.path.join(output_dir, 'final_repertoire_stripped.fa')

    def CheckCroppedReadsExistance(self):
        if not os.path.exists(self.cropped_reads):
            self.__log.info("ERROR: File containing cropped reads was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckBadReadsExistance(self):
        if not os.path.exists(self.bad_reads):
            self.__log.info("ERROR: File containing contaminated reads (not Ig-Seq) was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckVJAlignmentInfoExistance(self):
        if not os.path.exists(self.vj_alignment_info):
            self.__log.info("ERROR: File containing VJ alignment info was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckCompressedReadsExistance(self):
        if not os.path.exists(self.compressed_reads):
            self.__log.info("ERROR: File containing compressed reads was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckCroppedCompressedMapExistance(self):
        if not os.path.exists(self.map_file):
            self.__log.info("ERROR: File containing map from cropped reads to compressed reads was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckSupernodesExistance(self):
        if not os.path.exists(self.supernodes_file):
            self.__log.info("ERROR: File containing supernodes was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckSWGraphExistance(self):
        if not os.path.exists(self.sw_graph):
            self.__log.info("ERROR: File containing Smith-Waterman graph was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckDenseSubgraphDecompositionExistance(self):
        if not os.path.exists(self.dense_sgraph_decomposition):
            self.__log("ERROR: File containing dense subgraph decomposition was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckFinalClustersExistance(self):
        if not os.path.exists(self.final_clusters_fa):
            self.__log("ERROR: File containing clusters of final repertoire was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckFinalRCMExistance(self):
        if not os.path.exists(self.final_rcm):
            self.__log("ERROR: File containing RCM of final repertoire was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckFinalStrippedClustersExistance(self):
        if not os.path.exists(self.final_stripped_clusters_fa):
            self.__log("ERROR: File containing stripped clusters of final repertoire was not found")
            SupportInfo(self.__log)
            sys.exit(1)

#######################################################################################
#           Phases
#######################################################################################
class Phase:
    def __init__(self, long_name, log):
        self._long_name = long_name
        self._log = log

    def PrintStartMessage(self):
        self._log.info("==== " + self._long_name + " starts\n")

    def Run(self):
        print "This method should be overloaded"

    def PrintFinishMessage(self):
        self._log.info("\n==== " + self._long_name + " finished")

###########
class VJAlignmentPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetVJAlignmentLongName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        if not os.path.exists(self.__params.reads):
            self._log.info("ERROR: Input reads " + self.__params.reads + " were not found")
            SupportInfo(self._log)
            sys.exit(1)

    def __CheckOutputExistance(self):
        self.__params.io.CheckCroppedReadsExistance()
        self.__params.io.CheckBadReadsExistance()
        self.__params.io.CheckVJAlignmentInfoExistance()

    def Run(self):
        self.__CheckInputExistance()
        self.__params.vj_finder_output = os.path.join(self.__params.output, "vj_finder")
        command_line = IgRepConConfig().run_vj_aligner + " -i " + self.__params.reads + \
                       " -o " + self.__params.io.vj_finder_output + \
                       " --db-directory " + IgRepConConfig().path_to_germline + \
                       " -t " + str(self.__params.num_threads) + \
                       " -C " + self.__params.chain + \
                       " --organism " + self.__params.organism
        if self.__params.no_pseudogenes:
            command_line += " --no-pseudogenes"
        support.sys_call(command_line, self._log)

    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files: ")
        self._log.info("  * Cropped reads were written to " + self.__params.io.cropped_reads)
        self._log.info("  * Contaminated (not Ig-Seq) reads were written to " + self.__params.io.bad_reads)
        self._log.info("  * VJ alignment output was written to " + self.__params.io.vj_alignment_info)

###########
class TrieCompressionPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetTrieCompressorLongName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        self.__params.io.CheckCroppedReadsExistance()

    def __CheckOutputExistance(self):
        self.__params.io.CheckCompressedReadsExistance()
        self.__params.io.CheckCroppedCompressedMapExistance()
        self.__params.io.CheckSupernodesExistance()

    def Run(self):
        self.__CheckInputExistance()
        command_line = IgRepConConfig().run_trie_compressor + " -i " + self.__params.io.cropped_reads + \
                       " -o " + self.__params.io.compressed_reads + " -m " + self.__params.io.map_file
        support.sys_call(command_line, self._log)
        command_line = "%s %s %s --limit=%d" % (IgRepConConfig().run_report_supernodes,
                                                self.__params.io.compressed_reads,
                                                self.__params.io.supernodes_file,
                                                5)
        support.sys_call(command_line, self._log)

    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files:")
        self._log.info("  * Compressed reads were written to " + self.__params.io.compressed_reads)
        self._log.info("  * Supernode sequences were written to " + self.__params.io.supernodes_file)

###########
class GraphConstructionPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetGraphConstructionLongName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        self.__params.io.CheckCompressedReadsExistance()

    def __CheckOutputExistance(self):
        self.__params.io.CheckSWGraphExistance()

    def Run(self):
        self.__CheckInputExistance()
        command_line = IgRepConConfig().run_graph_constructor + " -i " + self.__params.io.compressed_reads + \
                       " -o " + self.__params.io.sw_graph + " -t " + str(self.__params.num_threads) + \
                       " --tau=" + str(self.__params.max_mismatches)
        support.sys_call(command_line, self._log)

    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files:")
        self._log.info("  * File containing Smith-Waterman graph was written to " + self.__params.io.sw_graph)

###########
class DSFPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetDSFLongName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        self.__params.io.CheckSWGraphExistance()

    def __GetDSFParams(self):
        dsf_params = ['-g', self.__params.io.sw_graph, '-o', self.__params.io.dsf_output, '-t',
                      str(self.__params.num_threads)]
        if self.__params.create_trivial_decomposition:
            dsf_params.append('--create-triv-dec')
        if self.__params.save_aux_files:
            dsf_params.append('--save-aux-files')
        return dsf_params

    def __CheckOutputExistance(self):
        self.__params.io.CheckDenseSubgraphDecompositionExistance()

    def Run(self):
        self.__CheckInputExistance()
        dense_subgraph_finder.main(self.__GetDSFParams(), self.__params.log_filename)

    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files:")
        self._log.info("  * File containing dense subgraph decomposition was written to " +
                       self.__params.io.dense_sgraph_decomposition)

###########
class ConsensusConstructionPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetConsensusConstructorLongName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        self.__params.io.CheckCompressedReadsExistance()
        self.__params.io.CheckDenseSubgraphDecompositionExistance()
        self.__params.io.CheckCroppedReadsExistance()
        self.__params.io.CheckCroppedCompressedMapExistance()

    def __CheckOutputExistance(self):
        self.__params.io.CheckFinalClustersExistance()
        self.__params.io.CheckFinalRCMExistance()

    def Run(self):
        self.__CheckInputExistance()
        command_line = "%s -i %s -c %s -q %s -o %s" % (IgRepConConfig().run_rcm_recoverer,
                                                       self.__params.io.cropped_reads,
                                                       self.__params.io.map_file,
                                                       self.__params.io.dense_sgraph_decomposition,
                                                       self.__params.io.final_rcm)
        support.sys_call(command_line, self._log)
        command_line = IgRepConConfig().run_consensus_constructor + \
                       " -i " + self.__params.io.cropped_reads + \
                       " -R " + self.__params.io.final_rcm + \
                       " -o " + self.__params.io.final_clusters_fa + \
                       " -H " + " -t " + str(self.__params.num_threads)
        support.sys_call(command_line, self._log)


    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files:")
        self._log.info("  * File containing clusters of final repertoire was written to " +
                       self.__params.io.final_clusters_fa)
        self._log.info("  * File containing read-cluster map of final repertoire was written to " +
                       self.__params.io.final_rcm)

class RemoveLowAbundanceReadsPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetRemoveLowAbundanceReadsName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        self.__params.io.CheckFinalClustersExistance()

    def __CheckOutputExistance(self):
        self.__params.io.CheckFinalStrippedClustersExistance()

    def Run(self):
        self.__CheckInputExistance()
        command_line = "%s %s %s --limit=%d" % (IgRepConConfig().run_remove_low_abundance_reads,
                                                self.__params.io.final_clusters_fa,
                                                self.__params.io.final_stripped_clusters_fa,
                                                IgRepConConfig().cluster_size_limit)
        support.sys_call(command_line, self._log)


    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files:")
        self._log.info("  * File containing stripped clusters of final repertoire was written to " +
                       self.__params.io.final_stripped_clusters_fa)

###########
class PhaseFactory:
    def __init__(self, phase_names, params, log):
        self.__phase_names = phase_names
        self.__entry_point = params.entry_point
        self.__params = params
        self.__log = log

    def __CreatePhaseByName(self, phase_name):
        if self.__phase_names.PhaseIsVJAlignment(phase_name):
            return VJAlignmentPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsTrieCompressor(phase_name):
            return TrieCompressionPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsGraphConstructor(phase_name):
            return GraphConstructionPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsDSF(phase_name):
            return DSFPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsConsensusConstructor(phase_name):
            return ConsensusConstructionPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsRemoveLowAbundanceReads(phase_name):
            return RemoveLowAbundanceReadsPhase(self.__params, self.__log)

    def CreatePhases(self):
        phase_list = list()
        first_phase_index = self.__phase_names.GetPhaseIndex(self.__entry_point)
        if first_phase_index == -1:
            self.__log.info("Incorrect name of entry-point")
            sys.exit(1)
        for i in range(first_phase_index, len(self.__phase_names )):
            phase_list.append(self.__CreatePhaseByName(self.__phase_names.GetPhaseNameBy(i)))
        return phase_list

############
class PhaseManager:
    def __init__(self, phase_factory, params, log):
        self.__params = params
        self.__log = log
        self.__phase_factory = phase_factory
        self.__phases = self.__phase_factory.CreatePhases()

    def __RunSinglePhase(self, phase_index):
            self.__phases[phase_index].PrintStartMessage()
            self.__phases[phase_index].Run()
            self.__phases[phase_index].PrintOutputFiles()
            self.__phases[phase_index].PrintFinishMessage()

    def __PrintPhaseDelimeter(self):
        self.__log.info("\n============================================\n")

    def Run(self):
        self.__RunSinglePhase(0)
        for i in range(1, len(self.__phases) - 1):
            self.__PrintPhaseDelimeter()
            self.__RunSinglePhase(i)
        if len(self.__phases) != 1:
            self.__PrintPhaseDelimeter()
            self.__RunSinglePhase(len(self.__phases) - 1)

#######################################################################################
#           IO routines
#######################################################################################

def CreateLogger():
    log = logging.getLogger('ig_repertoire_constructor')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log

def ParseCommandLineParams():
    from src.python_add.argparse_ext import ArgumentHiddenParser
    parser = ArgumentHiddenParser(description="IgRepertoireConstructor: an algorithm for construction of "
                                              "antibody repertoire from immunosequencing data",
                                  epilog="""
    In case you have troubles running IgRepertoireConstructor, you can write to igtools_support@googlegroups.com.
    Please provide us with ig_repertoire_constructor.log file from the output directory.
                                  """,
                                  add_help=False)

    req_args = parser.add_argument_group("Input")
    input_args = req_args.add_mutually_exclusive_group(required=True)
    input_args.add_argument("-s", "--reads",
                            type=str,
                            default="", # FIXME This is only for ace's version of python. Locally it works great w/o it
                            help="Immunosequencing reads in FASTQ format")
    input_args.add_argument("--test",
                            action="store_const",
                            const="test_dataset/merged_reads.fastq",
                            dest="reads",
                            help="Running of test dataset")

    out_args = parser.add_argument_group("Output")
    out_args.add_argument("-o", "--output",
                          type=str,
                          default="igrepcon_test",
                          help="Output directory [default: \"%(default)s\"]")

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-t", "--threads",
                               type=int,
                               default=16,
                               dest="num_threads",
                               help="Thread number [default: %(default)d]")
    optional_args.add_argument("--tau",
                               type=int,
                               default=4,
                               dest="max_mismatches",
                               help="Maximum allowed mismatches between identical error-prone reads "
                                    "[default: %(default)d]")

    optional_args.add_argument("-h", "--help",
                               action="help",
                               help="Show this help message and exit")

    vj_align_args = parser.add_argument_group("VJ alignment arguments")
    vj_align_args.add_argument("--no-pseudogenes",
                               action="store_const",
                               const=True,
                               dest="no_pseudogenes",
                               help="Do not use pseudogenes along with normal gene segments for VJ alignment [default: False]")
    vj_align_args.add_argument("--organism",
                               type=str,
                               default="human",
                               dest="organism",
                               help="Organism (human and mouse only are supported for this moment) [default: %(default)s]")
    vj_align_args.add_argument("-C", "--chain",
                               type=str,
                               dest="chain",
                               default="all",
                               help="Ig chain type: all (for both heavy and light chains)/ heavy / light (for both kappa and lambda chains) "
                                    "/ lambda / kappa [default: %(default)s]")

    dev_args = parser.add_argument_group("_Developer arguments")
    dev_args.add_argument("-f", "--min-fillin",
                          type=float,
                          default=0.6,
                          help="_Minimum edge fill-in of dense subgraphs [default: %(default)2.1f]")
    dev_args.add_argument('--entry-point',
                          type=str,
                          default=PhaseNames().GetPhaseNameBy(0),
                          help="_Continue from the given stage [default: %(default)s]")
    dev_args.add_argument("--create-triv-dec",
                          action="store_const",
                          const=True,
                          dest="create_trivial_decomposition",
                          help='_Creating decomposition according to connected components [default: False]')
    dev_args.add_argument("--save-aux-files",
                          action="store_const",
                          const=True,
                          dest="save_aux_files",
                          help="_Saving auxiliary files: subgraphs in GRAPH format and their decompositions "
                                    "[default: False]")

    ods_args = dev_args.add_mutually_exclusive_group(required=False)
    ods_args.add_argument("--help-hidden", "-H",
                          action="help_hidden",
                          help="_Show hidden help")
    parser.set_defaults(config_dir="configs",
                        config_file="config.info")
    params = parser.parse_args()
    return parser, params

def CheckParamsCorrectness(parser, params, log):
    if not "output" in params or params.output == "":
        log.info("ERROR: Output directory (-o) was not specified\n")
        parser.print_help()
        sys.exit(-1)
    if not "reads" in params or params.reads == "":
        log.info("ERROR: Reads (-s) were not specified\n")
        parser.print_help()
        sys.exit(-1)
    if not os.path.exists(params.reads):
        log.info("ERROR: File with reads " + params.reads + " were not found\n")
        parser.print_help()
        sys.exit(-1)
    if not os.path.isabs(params.reads):
        params.reads = os.path.abspath(params.reads)

def PrepareOutputDir(params):
    if params.entry_point == "vj_alignment" and os.path.exists(params.output):
        shutil.rmtree(params.output)
    if not os.path.isdir(params.output):
        os.makedirs(params.output)

def PrintParams(params, log):
    log.info("IgRepertoireConstructor parameters:")
    log.info("  Input reads:\t\t\t" + params.reads)
    log.info("  Output directory:\t\t" + params.output)
    log.info("  Number of threads:\t\t" + str(params.num_threads))
    log.info("  Maximal number of mismatches:\t" + str(params.max_mismatches))
    log.info("  Entry point:\t\t\t" + params.entry_point)

def CreateFileLogger(params, log):
    params.log_filename = os.path.join(params.output, "ig_repertoire_constructor.log")
    if os.path.exists(params.log_filename):
        os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
    log.addHandler(log_handler)
    log.info("Log will be written to " + params.log_filename + "\n")

def PrintCommandLine(log):
    command_line = "Command line: " + " ".join(sys.argv)
    log.info("\n" + command_line + "\n")

#######################################################################################
#           Main
#######################################################################################
def main():
    binary_config = IgRepConConfig()
    log = CreateLogger()
    binary_config.CheckBinaries(log)
    parser, params = ParseCommandLineParams()
    CheckParamsCorrectness(parser, params, log)
    PrepareOutputDir(params)
    CreateFileLogger(params, log)
    PrintParams(params, log)
    PrintCommandLine(log)
    params.io = IgRepConIO(params.output, log)

    try:
        ig_phase_factory = PhaseFactory(PhaseNames(), params, log)
        ig_repertoire_constructor = PhaseManager(ig_phase_factory, params, log)
        ig_repertoire_constructor.Run()
        log.info("\nThank you for using IgRepertoireConstructor!")
    except (KeyboardInterrupt):
        log.info("\nIgRepertoireConstructor was interrupted!")
    except Exception:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught.")
            SupportInfo(log)
    except BaseException:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught.")
            SupportInfo(log)

    log.info("Log was written to " + params.log_filename)

if __name__ == '__main__':
    main()
