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

# coloredlogs
import coloredlogs

import RNA
from helper import p_to_s, print_moves
from findpath_i_sorted import find_path
from findpath_i_findmoves import find_moves


if __name__ == '__main__':

     
    search_width = 50
    Verbose = True
    # Debug = True
    Debug = False


    o_filename = r"local_min_100_multiple_sections_min10.csv"
    # o_filename = r"local_min_200_multiple_sections_min10.csv"
    # o_filename = r"local_min_300_multiple_sections_min10.csv"
    # o_filename = r"local_min_400_multiple_sections_min10.csv"
    # o_filename = r"local_min_500_multiple_sections_min10.csv"
    # o_filename = r"local_min_600_multiple_sections_min10.csv"
    # o_filename = r"local_min_800_multiple_sections_min10.csv"
    # o_filename = r"local_min_1000_multiple_sections_min10.csv"

    filename = r"./sample_seqs/" + o_filename

    df = pd.read_csv(filename)
    elements = len(df.index)


    sequences = []
    s1s = []
    s2s = []
    indices = []
    seq_lengths = []

    all_indirect_moves = []
    all_en = []

    for index, row in df.iterrows():

        if index > 90:
            break

        # if index != 1:
        #     continue
    

        sequence = row.sequence
        s1 = row.s1
        s2 = row.s2  

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



        print(f"sequence = \"{sequence}\"")
        print(f"s1 = \"{s1}\"")
        print(f"s2 = \"{s2}\"")

        fc = RNA.fold_compound(sequence)
        s1_eval = round(fc.eval_structure(s1), 2)

        search_width = RNA.bp_distance(s1, s2) * 2

        # break

        candidates = find_moves(sequence, s1, s2, search_width=500, Verbose=Verbose)
        s, count, indirect_moves, search_width, moves = candidates[0]
        
        all_indirect_moves.append(list(indirect_moves).copy())
        all_en.append(s)

        print (indirect_moves)
        print_moves(sequence, s1, s2, moves, Verbose=Verbose)

        # break


    df = pd.DataFrame([indices, sequences, s1s, s2s, seq_lengths, all_indirect_moves, all_en]).transpose()
    df.columns = ["i", "sequence", "s1", "s2", "seq_length", "indirect_moves", "saddle_en"]

    prefix = "indirect_"
    savefile = r"./results/" + prefix + o_filename
    df.to_csv(savefile)

    
    print(df)