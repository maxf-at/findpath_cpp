#!/usr/bin/env python3
# coding: utf-8

# eval.c benchmark: line 356: PUBLIC int vrna_eval_move_pt(vrna_fold_compound_t
# eval.c vrna_eval_loop_pt_v line 270

# vrna_eval_loop_pt = vrna_eval_loop_pt_v



import RNA
import numpy as np

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


import coloredlogs

import pandas as pd


filename = r"../sample_seqs/eval_move_dataset.csv"


df = pd.read_csv(filename)
elements = len(df)


time_elapsed = 0

# look at 2

for index, row in df.iterrows():


    # if index>=100:
    #     break
    
    

    sequence = row.sequence
    s1 = row.s1
    i = row.i_move
    j = row.j_move

    fc = RNA.fold_compound(sequence)
    pt = RNA.ptable(s1)

    fc.eval_structure_pt(pt)

    for _ in range(1):
        start = time.time()
        en = fc.eval_move_pt(pt, i, j)
        end = time.time()-start
        time_elapsed += end


    percent_complete = 100-(elements-index-1)/elements*100
    bar_length = 20
    filled_length = int(percent_complete/100*bar_length)
    rest = bar_length - filled_length
    bar = "█" * filled_length + '_' * rest

    print(
        f'\rComputing |{bar}| {percent_complete:.1f}% complete', end="\r")


print()
print("runtime:", time_elapsed)


# init dataframe, then save data as csv
# data = [runtimes_merge_outer, runtimes_merge_inner, runtimes_merge_join, runtimes_regular, results_merge, results_cpp]  
  
# df = pd.DataFrame(data)

# # transpose
# df = df.transpose() 

# df.columns = ["runtimes_merge_outer", "runtimes_merge_inner", "runtimes_merge_join", "runtimes_regular", "results_merge", "results_regular"]

# print (df)
# filename = f'../documents/benchmarks/local_min_300_2_sections.csv'
# filename = f'../documents/benchmarks/local_min_400_2_sections.csv'
# filename = f'../documents/benchmarks/local_min_500_2_sections.csv'
# df.to_csv(filename)



# print ("cpp runtime ext", runtimes_sep_ext) # post-processing networkx

# 127.82060360908508 <=

# new opt a/b
# better 8
# worse 4
# new opt a/b
# better 5
# worse 9


# return b->opt - a->opt; a/b
# better 5 / 7
# worse 8 / 8


#return b->curr_en - a->curr_en;
# better 7
# worse 9