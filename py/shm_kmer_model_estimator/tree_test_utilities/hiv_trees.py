#!/usr/bin/env python2

from __future__ import print_function

import os
import re

import numpy as np

from special_utils.cd_manager import cd
from special_utils.largest_files_in_dir import n_largest_files
from shm_kmer_model import cab_shm_model, yale_model
import shm_kmer_model_estimator.standard_model_estimations as \
    standard_model_estimations
from chains import chains

from disk_memoize.disk_memoize import memoize_to_disk 


def get_hiv_trees_path(chain_type=chains.Chains.IGH, n_largest_trees=5,
        data_path_prefix='/Sid/abzikadze/antevolo_paper_results/antevolo/HIV/antevolo_backup'):
    chain_type = chain_type.name
    return get_hiv_trees_path__(chain_type, n_largest_trees, data_path_prefix)


@memoize_to_disk
def get_hiv_trees_path__(chain_type, n_largest_trees, data_path_prefix):
    paths = []
    chain_path = os.path.join(data_path_prefix, chain_type)
    hiv_ind_names = os.listdir(chain_path)
    print(hiv_ind_names)
    for ind_name in hiv_ind_names:
        clonal_trees = os.path.join(chain_path, ind_name, 'clonal_trees')
        paths += n_largest_files(clonal_trees, n_largest_trees)
    return paths
