#!/usr/bin/env python3
# coding: utf-8

import RNA
import numpy as np

from sklearn import preprocessing
min_max_scaler = preprocessing.MinMaxScaler()

import subprocess
import matplotlib.pyplot as plt
from PIL import Image
from IPython.display import SVG, display
from collections import Counter
from collections import defaultdict

import difflib
import sys
import os
import random
import string
import time

# from helper import print_moves

import pandas as pd

sys.path.append('../')
from pretty_print_path import print_moves
import findpath_librna
import findpath


# find a sequence + structure pair as example (good sw-multiplier scaling)

sys.path.append('./indirect')
import helper
target_count = 1000
x = 300
min_bp_dist = 14

def detect_local_minimum(fc, structure):
    # perform gradient walk from sample to determine direct local minimum
    pt = RNA.IntVector(RNA.ptable(structure))
    fc.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
    return RNA.db_from_ptable(list(pt))

# check confirmations for merging
s_list = []
while (len(s_list) < target_count):

    # RNAsubopt wrapper
    sequence, s1, s2 = helper.generate_structures(x)

    # overwrite s1 s2 to local minimum
    fc = RNA.fold_compound(sequence)
    s1 = detect_local_minimum(fc, s1)
    s2 = detect_local_minimum(fc, s2)
    bp_dist = RNA.bp_distance(s1, s2)

    if bp_dist < min_bp_dist:
        continue

    # ensure that we only process paths which are not trivial to process - if we get the best result at a 
    # low search width, its probably not a good candidate for training

    sws = [0.5, 1, 2]
    results = []
    for sw in sws:
        fp = findpath.findpath_single(sequence, s1, s2, search_width_multiplier=sw, mp=True)
        result = fp.get_en()/100.0
        results.append(result)
        
    # if results[0] != results[1] and results[1] != results[2]:
    if results[0] != results[2]:
        
        # print(f"sequence = '{sequence}'")
        # print(f"s1       = '{s1}'")
        # print(f"s2       = '{s2}'")
        # print (results)
        
        s_list.append((sequence, s1, s2))

        if len(s_list)%100 == 0:
            print (len(s_list), end=' ')

            # save seqs every 10 iterations.
            df = pd.DataFrame(s_list, columns=['sequence', 's1', 's2']).set_index('sequence')
            filename = f'./dataset_{x}_large.csv'
            df.to_csv(filename)

        

    # break