#!/usr/bin/env python2

from simulate import *
from aimquast_impl import get_clusters_sizes
import os.path
from joblib import Parallel, delayed

path_to_aimquast = igrec_dir + "/aimquast.py"


def run_and_quast_all(input_reads,
                      ideal_repertoire_fa,
                      ideal_repertoire_rcm,
                      out_dir,
                      threads=4,
                      rerun_mixcr=True,
                      do_not_run=False):
    import os.path
    import shutil

    mkdir_p(out_dir)

    class IgReCRun:

        def __init__(self,
                     name,
                     tau=4,
                     fillin=0.6,
                     min_sread_size=5,
                     max_votes=1000500,
                     trivial=False,
                     loci="all",
                     additional_args=" "):
            self.name = name
            self.tau = tau
            self.fillin = fillin
            self.trivial = trivial
            self.loci = loci
            self.min_sread_size = min_sread_size
            self.additional_args = additional_args
            self.max_votes = max_votes

        def run(self):
            additional_args = " --no-alignment --max-votes=%d" % self.max_votes
            if self.trivial:
                additional_args += " --create-triv-dec "
            additional_args += " " + self.additional_args
            run_igrec(input_reads,
                      tau=self.tau,
                      additional_args=additional_args,
                      min_fillin=self.fillin,
                      threads=threads,
                      loci=self.loci,
                      min_sread_size=self.min_sread_size,
                      output_dir=out_dir + "/" + self.name + "/")

    igrec_runs = []
    # igrec_runs.append(IgReCRun("igrec_trivial", trivial=True))
    igrec_runs.append(IgReCRun("igrec"))
    # igrec_runs.append(IgReCRun("igrec_nomsns", min_sread_size=1005000))
    # igrec_runs.append(IgReCRun("igrec_msns2", min_sread_size=2))
    # igrec_runs.append(IgReCRun("igrec_msns3", min_sread_size=3))
    igrec_runs.append(IgReCRun("igrec_vote", max_votes=1))
    # igrec_runs.append(IgReCRun("igrec_vote2", max_votes=2))
    # igrec_runs.append(IgReCRun("igrec_trivial_tau3", tau=3, trivial=True))
    # igrec_runs.append(IgReCRun("igrec_trivial_tau2", tau=2, trivial=True))
    # igrec_runs.append(IgReCRun("igrec_trivial_tau1", tau=1, trivial=True))

    # igrec_runs.append(IgReCRun("igrec", additional_args="--debug"))
    # igrec_runs.append(IgReCRun("igrec_split", additional_args="--no-equal-compression --debug"))
    # igrec_runs.append(IgReCRun("igrec_tau3", tau=3))
    # igrec_runs.append(IgReCRun("igrec_split_tau3", tau=3, additional_args=" --no-equal-compression --debug"))
    # igrec_runs.append(IgReCRun("igrec_tau2", tau=2))
    # igrec_runs.append(IgReCRun("igrec_tau1", tau=1))

    # igrec_runs.append(IgReCRun("igrec_tau3_msns2", tau=3, min_sread_size=2))

    # igrec_runs.append(IgReCRun("igrec_tau3_vote", tau=3, max_votes=1))
    # igrec_runs.append(IgReCRun("igrec_tau3_vote2", tau=3, max_votes=2))
    # igrec_runs.append(IgReCRun("igrec_tau2_msns2", tau=2, min_sread_size=2))
    # igrec_runs.append(IgReCRun("igrec_tau1_msns2", tau=1, min_sread_size=2))

    # # igrec_runs.append(IgReCRun("igrec_msns3", min_sread_size=3))
    # # igrec_runs.append(IgReCRun("igrec_msns4", min_sread_size=4))
    # # igrec_runs.append(IgReCRun("igrec_f03", fillin=0.3))
    # # igrec_runs.append(IgReCRun("igrec_f075", fillin=0.75))
    # # igrec_runs.append(IgReCRun("igrec_f09", fillin=0.9))
    # # igrec_runs.append(IgReCRun("igrec_f095", fillin=0.95))
    # igrec_runs.append(IgReCRun("igrec_tau3", tau=3))
    # igrec_runs.append(IgReCRun("igrec_tau3_msns3", tau=3, min_sread_size=3))
    # igrec_runs.append(IgReCRun("igrec_tau3_f03", tau=3, fillin=0.3))
    # # igrec_runs.append(IgReCRun("igrec_tau3_f075", tau=3, fillin=0.75))
    # # igrec_runs.append(IgReCRun("igrec_tau3_f09", tau=3, fillin=0.9))
    # igrec_runs.append(IgReCRun("igrec_tau3_f095", tau=3, fillin=0.95))
    # igrec_runs.append(IgReCRun("igrec_tau2_msns2", tau=2, min_sread_size=2))
    # # igrec_runs.append(IgReCRun("igrec_tau2", tau=2))
    # # igrec_runs.append(IgReCRun("igrec_tau2_f03", tau=2, fillin=0.3))
    # # igrec_runs.append(IgReCRun("igrec_tau2_f075", tau=2, fillin=0.75))
    # # igrec_runs.append(IgReCRun("igrec_tau2_f09", tau=2, fillin=0.9))
    # igrec_runs.append(IgReCRun("igrec_tau2_f095", tau=2, fillin=0.95))
    # igrec_runs.append(IgReCRun("igrec_tau2_f099", tau=2, fillin=0.99))
    # # igrec_runs.append(IgReCRun("igrec_tau1", tau=1))
    # # igrec_runs.append(IgReCRun("igrec_tau1_f03", tau=1, fillin=0.3))
    # # igrec_runs.append(IgReCRun("igrec_tau1_f075", tau=1, fillin=0.75))
    # # igrec_runs.append(IgReCRun("igrec_tau1_f09", tau=1, fillin=0.9))
    # # igrec_runs.append(IgReCRun("igrec_tau1_f095", tau=1, fillin=0.95))

    if not do_not_run:
        for run in igrec_runs:
            run.run()

        if rerun_mixcr or not os.path.isfile(out_dir + "/mixcr/final_repertoire.fa"):
            run_mixcr2(input_reads, threads=threads, output_dir=out_dir + "/mixcr/", loci="all")

        run_presto(input_reads, output_dir=out_dir + "/presto/")

        mkdir_p(out_dir + "/supernode")
        shutil.copy(out_dir + "/" + igrec_runs[0].name + "/supernode_repertoire.fa",
                    out_dir + "/supernode/final_repertoire.fa")
        shutil.copy(out_dir + "/" + igrec_runs[0].name + "/supernode_repertoire.rcm",
                    out_dir + "/supernode/final_repertoire.rcm")

    kinds = [run.name for run in igrec_runs] + ["supernode", "mixcr"]

    for kind in kinds:
        args = {"ideal_repertoire_fa": ideal_repertoire_fa,
                "ideal_repertoire_rcm": ideal_repertoire_rcm,
                "input_reads": input_reads,
                "out_dir": out_dir,
                "kind": kind}
        cmd = path_to_aimquast + " -s %(input_reads)s -r %(ideal_repertoire_fa)s -R %(ideal_repertoire_rcm)s -c %(out_dir)s/%(kind)s/final_repertoire.fa -o %(out_dir)s/%(kind)s/aimquast --no-reference-free -F png,pdf --rcm-based --reference-free" % args

        rcm = "%(out_dir)s/%(kind)s/final_repertoire.rcm" % args
        if os.path.isfile(rcm):
            cmd += " -C %s" % rcm

        os.system(cmd)


if __name__ == "__main__":
    mkdir_p("various_error_rate")
    mkdir_p("var_err_rate_real")

    ig_simulator_output_dir = "/tmp/ig_simulator"
    datasets = [ig_simulator_output_dir + "/final_repertoire.fasta",
                igrec_dir + "/var_err_rate_real/error_free_reads.fa.gz"]
    output_dirs = [igrec_dir + "/various_error_rate", igrec_dir + "/var_err_rate_real"]

    if True:
        run_ig_simulator(ig_simulator_output_dir,
                         chain="HC", num_bases=100, num_mutated=1000, reprtoire_size=5000)

        try:
            convert_abvitro_to_repertoire("/Jake/data/input/ImmunoSeq/AbVitro/flu_time_course/FV/assembled_umis/21_assemble_combined.fastq",
                                          igrec_dir + "/var_err_rate_real/flu_repertoire.fa.gz")
            multiplex_repertoire(igrec_dir + "/var_err_rate_real/flu_repertoire.fa.gz",
                                 igrec_dir + "/var_err_rate_real/error_free_reads.fa.gz")
        except BaseException as ex:
            print ex
            print "Cannot multiplex reperoire, file not found"

        for dataset, output_dir in zip(datasets, output_dirs):
            if not os.path.isfile(dataset):
                continue
            simulate_data_wo_errors(dataset,
                                    output_dir + "/data")

    # lambdas = [0, 0.0625, 0.125, 0.25, 0.375, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4]
    lambdas = [0, 0.0625, 0.125, 0.25, 0.375, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
    for output_dir in output_dirs:
        if not os.path.isfile(output_dir + "/data/error_free_reads.fa.gz"):
            continue
        min_error_interval = [0, 1] if "real" not in output_dir else [0]
        for min_error in min_error_interval:
            def JOB(error_rate):
                out_dir = output_dir + "/errate_%0.4f" % error_rate if not min_error else output_dir + "/errate_%0.4f_woans" % error_rate
                simulate_data_from_dir(output_dir + "/data",
                                       out_dir,
                                       error_rate=error_rate,
                                       seed=0,
                                       min_error=min_error,
                                       erroneous_site_len=300)

                sizes = get_clusters_sizes(output_dir + "/data/ideal_final_repertoire.fa.gz")
                print "Reference consists of %d clusters" % len(sizes)
                print "Reference consists of %d large (>=5) clusters" % len([size for size in sizes if size >= 5])

                run_and_quast_all(out_dir + "/merged_reads.fa.gz",
                                  output_dir + "/data/ideal_final_repertoire.fa.gz",
                                  output_dir + "/data/ideal_final_repertoire.rcm", out_dir,
                                  rerun_mixcr=True)

            import multiprocessing
            n_jobs = 1 if multiprocessing.cpu_count() <= 16 else 4
            Parallel(n_jobs=n_jobs)(delayed(JOB)(error_rate) for error_rate in lambdas)
