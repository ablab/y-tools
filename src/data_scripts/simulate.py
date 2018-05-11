#!/usr/bin/env python2
import shutil

import sys
import os
import subprocess


ROOT_PATH = os.path.abspath(os.path.join(__file__, os.path.pardir, os.path.pardir, os.path.pardir))
CONFIG_PATHS = os.path.abspath(os.path.join(__file__, os.path.pardir, "simulation_configs"))
BIN_PATH = os.path.join(ROOT_PATH, "build/release/bin")
IG_SIMULATOR_BIN = os.path.join(BIN_PATH, "ig_simulator")
IG_SIMULATOR_PATH = os.path.join(ROOT_PATH, "ig_simulator.py")
PCR_SIMULATOR_PATH = os.path.join(BIN_PATH, "simulate_barcoded")

METAROOTS_BY_STRATEGY = {
    "naive": "20000",
    "plasma": "100",
}


def simulate_all(output_dir):
    simulate(output_dir, "plasma")
    simulate(output_dir, "naive")


def simulate(output_dir, strategy_name):
    strategy_output_dir = os.path.join(output_dir, strategy_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    simulation_output_dir = os.path.join(strategy_output_dir, "simulation")
    subprocess.check_call([IG_SIMULATOR_PATH,
                           "--config", os.path.join(CONFIG_PATHS, strategy_name + "_config.info"),
                           "--output", simulation_output_dir,
                           "-n", METAROOTS_BY_STRATEGY[strategy_name]])
    simulated_reads = os.path.join(simulation_output_dir, "filtered_pool.fasta")

    simulate_with_pcr(simulated_reads, strategy_name, strategy_output_dir)

    simulate_no_pcr(simulated_reads, simulation_output_dir, strategy_name, strategy_output_dir)


def simulate_no_pcr(simulated_reads, simulation_output_dir, strategy_name, strategy_output_dir):
    # all simulations without PCR should have all reads repeated several times (5 at least)
    # replacing with zero rate PCR now
    pcr_output_dir = os.path.join(strategy_output_dir, "pcr0")
    subprocess.call([PCR_SIMULATOR_PATH,
                     "--input-file", simulated_reads,
                     "--output-dir", pcr_output_dir,
                     "--pcr-error1", "0.0",
                     "--pcr-error2", "0.0", ])
    no_pcr_dir = os.path.join(strategy_output_dir, strategy_name + "_no_pcr")
    prepare_test_input(
        target_dir = no_pcr_dir,
        reads_path = os.path.join(pcr_output_dir, "amplified.fasta"),
        repertoire_path = simulated_reads,
        rcm_file = os.path.join(pcr_output_dir, "amplified_to_orig.rcm")
    )

    # id_rcm_file = os.path.join(simulation_output_dir, "identity.rcm")
    # generate_identity_rcm(id_rcm_file, simulated_reads)
    # no_pcr_dir = os.path.join(strategy_output_dir, strategy_name + "_no_pcr")
    # prepare_test_input(
    #     no_pcr_dir,
    #     simulated_reads,
    #     simulated_reads,
    #     id_rcm_file
    # )


def simulate_with_pcr(simulated_reads, strategy_name, strategy_output_dir):
    pcr_output_dir = os.path.join(strategy_output_dir, "pcr")
    subprocess.call([PCR_SIMULATOR_PATH,
                     "--input-file", simulated_reads,
                     "--output-dir", pcr_output_dir,
                     "--pcr-error1", "0.003",
                     "--pcr-error2", "0.003", ])
    pcr_dir = os.path.join(strategy_output_dir, strategy_name + "_with_pcr")
    prepare_test_input(
        target_dir = pcr_dir,
        reads_path = os.path.join(pcr_output_dir, "amplified.fasta"),
        repertoire_path = simulated_reads,
        rcm_file = os.path.join(pcr_output_dir, "amplified_to_orig.rcm")
    )


def generate_identity_rcm(output_path, input_path):
    with open(output_path, "w") as out_file, open(input_path, "r") as in_file:
        current = 0
        for line in in_file:
            if not line.strip().startswith(">"):
                continue
            out_file.write(line.strip()[1:] + "\t" + str(current) + "\n")
            current += 1


def generate_repertoire(repertoire_path, rcm_file, compressed):
    from Bio import SeqIO
    from collections import defaultdict

    with open(repertoire_path) as repertoire_file:
        records = list(SeqIO.parse(repertoire_file, "fasta"))

    cluster_size = defaultdict(int)
    for line in open(rcm_file):
        columns = line.strip().split("\t")
        if len(columns) != 2:
            continue
        _, cluster_idx = columns
        cluster_size[int(cluster_idx)] += 1

    for idx in range(len(records)):
        records[idx].id = "cluster___" + str(idx) + "___size___" + str(cluster_size[idx])
        records[idx].description = ""


    tmp_path = compressed + ".tmp"
    with open(tmp_path, "w") as tmp_file:
        for record in records:
            SeqIO.write(record, tmp_file, "fasta")

    compress(tmp_path, compressed)
    os.remove(tmp_path)


def prepare_test_input(target_dir, reads_path, repertoire_path, rcm_file):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    compress(reads_path, os.path.join(target_dir, "input_reads.fa.gz"))
    # compress(repertoire_file, os.path.join(target_dir, "repertoire.fa.gz"))
    generate_repertoire(repertoire_path, rcm_file, os.path.join(target_dir, "repertoire.fa.gz"))
    shutil.copy2(rcm_file, os.path.join(target_dir, "repertoire.rcm"))


def compress(original, compressed):
    import gzip
    with open(original, "rb") as infile, gzip.open(compressed, "wb") as outfile:
        shutil.copyfileobj(infile, outfile)


if __name__ == '__main__':
    output_dir = sys.argv[1]
    simulate_all(output_dir)