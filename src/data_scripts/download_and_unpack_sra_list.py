#!/usr/bin/env python2

import subprocess
import sys
import os

CSV_DELIMITER = ","
DATASET_ID_COLUMN_NAME = "Run"
PATH_TO_ASPERA_SCRIPT = "/home/sbankevich/scripts/aspera.sh"


def download_and_unpack_sra_list(csv_file, destination):
    if not os.path.exists(destination):
        os.makedirs(destination)
    dataset_ids = get_dataset_ids(csv_file)
    for dataset_id in dataset_ids:
        print "Processing dataset", dataset_id
        sys.stdout.flush()
        dataset_destination = os.path.join(destination, dataset_id)
        print "Downloading to", dataset_destination
        sys.stdout.flush()
        subprocess.call([PATH_TO_ASPERA_SCRIPT, dataset_id, os.path.join(destination, dataset_id)])
        dataset_sra_file = os.path.join(destination, dataset_id + ".sra")
        print "Running fastq-dump on", dataset_sra_file
        sys.stdout.flush()
        subprocess.check_call(["fastq-dump", "--split-files", dataset_sra_file, "--gzip"], cwd = destination)
        fastq_dump_out = os.path.join(destination, dataset_id + "_1.fastq.gz")
        check_file_existence(fastq_dump_out)
        fastq_file = os.path.join(destination, dataset_id + ".fastq.gz")
        os.rename(fastq_dump_out, fastq_file)
        os.remove(dataset_sra_file)
        print "Result written to", fastq_dump_out
        sys.stdout.flush()


def get_dataset_ids(file_path):
    csv_file = open(file_path)
    dataset_id_column = get_dataset_id_column(csv_file.readline().strip())
    result = []
    for line in csv_file.readlines():
        dataset_id = line.strip().split(CSV_DELIMITER)[dataset_id_column]
        result.append(dataset_id)
    return result


def get_dataset_id_column(line):
    return line.split(CSV_DELIMITER).index(DATASET_ID_COLUMN_NAME)


def check_file_existence(filename, message="", log=None, dipspades=False):
    filename = os.path.abspath(filename)
    if not os.path.isfile(filename):
        print "file not found: %s (%s)" % (filename, message)
        sys.exit(1)
    return filename


if __name__ == '__main__':
    download_and_unpack_sra_list(sys.argv[1], sys.argv[2])
