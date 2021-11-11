#!/usr/bin/env python3
# coding: utf-8


import RNA

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
from IPython.display import SVG, display

from collections import Counter
from collections import defaultdict
import subprocess
import difflib
import sys
import os
import random
import string
import time

sys.path.append('../')
from pretty_print_path import print_moves
import findpath_librna
import findpath


# import feature_generation
from features import ij_distance, new_move_dist, plt_moves, config_distance, balance_in_all_things, return_shift_moves
from process_features import fp_call, find_moves, process


sequence = 'AGGUGUUAUAAAAGUGCAGAGCCAAGUCACCGCCUUGCCGUCGCUUCCGCUAUGCCUCAAAAGCGGAAUAUCGAGCGAGCGAUUGUGAUUCGAAUGGUCG'
s1       = '.(((................))).((((((....(((((((((.(((((((.(......).)))))))...))).)).))))..))))))..........'
s2       = '...(((.........)))..((((((((((.(((((((((....(((((((.(......).)))))))...)).))))).).).))))))....))))..'
search_width_multiplier = 2

def adjust_pt(pt, i, j):
    pt_adj = pt.copy()
    if i<0:
        pt_adj[-i] = 0
        pt_adj[-j] = 0
    else:
        pt_adj[i] = j
        pt_adj[j] = i
    return pt_adj

def structure_evaluation(sequence, pt1, pt2):
    avail_moves = []
    fp_results = []
    # pt1 = list(RNA.ptable(s))
    # pt2 = list(RNA.ptable(s2))
    for pos, (i,j) in enumerate(find_moves(pt1, pt2)):    
        next_en = fc.eval_move_pt(pt1, i, j)
        # mark where we found our move
        avail_moves.append((i, j, next_en))

    avail_moves.sort(key=lambda x: x[2])

    if len(avail_moves) == 0:
        return 0.0

    for pos, (i,j, en) in enumerate(avail_moves):  
        pt = adjust_pt(pt1, i, j)
        snew = RNA.db_from_ptable(pt)
        # if i > 0:
        #     snew = s[:i-1] + "(" + s[i:j-1] + ")" + s[j:]
        # if i < 0:
        #     snew = s[:-i-1] + "." + s[-i:-j-1] + "." + s[-j:]

        ptnew = list(RNA.ptable(snew))
        result_new, path = fp_call(sequence, snew, s2, search_width_multiplier)
        fp_results.append(result_new)

    best_result = np.argmin(fp_results)/len(fp_results)
    return best_result



fc = RNA.fold_compound(sequence)
pt1 = list(RNA.ptable(s1))
pt2 = list(RNA.ptable(s2))
structure_evaluation(sequence, pt1, pt2)


sE, path = fp_call(sequence, s1, s2, search_width_multiplier)

ptlast = list(RNA.ptable(s1))
for i, j, en in path:

    if i==0: 
        continue

    pt = adjust_pt(ptlast, i, j)
    # print (i, RNA.db_from_ptable(pt))

    eval = structure_evaluation(sequence, pt, pt2)

    s = RNA.db_from_ptable(pt)
    print (i, eval)
    

    ptlast = pt



