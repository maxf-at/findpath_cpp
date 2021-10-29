#!/usr/bin/env python3
# coding: utf-8

# generic Python libraries
import numpy as np
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
from PIL import Image
from IPython.display import SVG, display
from collections import Counter
from collections import defaultdict
from dataclasses import dataclass
import itertools
import difflib
import sys
import os
import random
import string
import time
import logging
import re

# coloredlogs
import coloredlogs

import RNA
from helper import p_to_s, print_moves
from findpath_i_findmoves import find_moves

from findpath_i_sorted import find_path as find_path_sorted
from findpath_i_unsorted import find_path as find_path_unsorted

if __name__ == '__main__':


    search_width = 50
    Verbose = True
    # Debug = True
    Debug = False

    sws = [2]


    o_filename = r"indirect_local_min_100_multiple_sections_min10.csv"


    filename = r"./results/" + o_filename

    df = pd.read_csv(filename)
    elements = len(df.index)


    sequences = []
    s1s = []
    s2s = []
    indices = []
    seq_lengths = []

    all_indirect_moves = []
    
    all_en1 = []
    all_en2 = []

    en1_runtimes = []
    en2_runtimes = []


    for index, row in df.iterrows():

        # if index != 1:
        #     continue

        sequence = row.sequence
        s1 = row.s1
        s2 = row.s2  


        # indirect_moves = row.indirect_moves.split(",")
        # indirect_moves = [int(re.sub("[^0-9]", '', i)) for i in indirect_moves]

        indirect_moves = re.findall(r'\d+', row.indirect_moves)
        indirect_moves = [int(i) for i in indirect_moves]

        indirect_moves = [(indirect_moves[i*2], indirect_moves[i*2+1]) for i in range(int(len(indirect_moves)/2))  ]
        
        if len(indirect_moves) == 0:
            continue


        print (indirect_moves)


        indices.append(index)
        sequences.append(sequence)
        seq_lengths.append(len(sequence))
        s1s.append(s1)
        s2s.append(s2)



        percent_complete = 100-(elements-index)/elements*100
        bar_length = 20
        filled_length = int(percent_complete/100*bar_length)
        rest = bar_length - filled_length
        bar = "█" * filled_length + '_' * rest

        print(
            f'Computing |{bar}| {percent_complete:.1f}% complete {index} ')

        sequence = row.sequence
        s1 = row.s1
        s2 = row.s2  

        print(f"sequence = \"{sequence}\"")
        print(f"s1 = \"{s1}\"")
        print(f"s2 = \"{s2}\"")

        # print (indirect_moves)

        fc = RNA.fold_compound(sequence)
        s1_eval = round(fc.eval_structure(s1), 2)



        start_findpath = time.time()
        # different search width multipliers
        for swm in sws:
            search_width = RNA.bp_distance(s1, s2) * swm
            paths = find_path_unsorted(sequence, s1, s2, add_moves=indirect_moves,
                            search_width=search_width, Debug=Debug, Verbose=Verbose)            
            saddle_en = paths[0][0]
            all_en1.append(saddle_en)
        en1_runtimes.append(round(time.time()-start_findpath, 4))

        
        start_findpath = time.time()
        # different search width multipliers
        for swm in sws:
            search_width = RNA.bp_distance(s1, s2) * swm
            paths = find_path_sorted(sequence, s1, s2, add_moves=indirect_moves,
                            search_width=search_width, Debug=Debug, Verbose=Verbose)            
            saddle_en = paths[0][0]
            all_en2.append(saddle_en)
        en2_runtimes.append(round(time.time()-start_findpath, 4))




    data = [indices, sequences, s1s, s2s, seq_lengths, all_indirect_moves, all_en1, all_en2, en1_runtimes, en2_runtimes]

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = ["i", "sequence", "s1", "s2", "seq_length", "indirect_moves", "en1", "en2", "en1_runtimes", "en2_runtimes"]

    prefix = "indirect_results_"

    savefile = r"./results/" + prefix + o_filename
    df.to_csv(savefile)
    print(df)