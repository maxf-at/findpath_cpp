#!/usr/bin/env python3
# coding: utf-8

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

import pathfinder_python
import pathfinder
import pathfinder_cpp
import helper
import merge_composition
import merge_recursive
import merge_dp

import pandas as pd

# generate X random sequences ~ x nt for merge join

# target_count = 1000
# min_bp_dist = 60
# max_bp_dist = 120
# x = 200
# target_count = 500
# min_bp_dist = 40
# max_bp_dist = 80
# x = 150
# target_count = 400
# min_bp_dist = 10
# max_bp_dist = 30
# x = 60

# indirect generation
# target_count = 100
# min_bp_dist = 10
# max_bp_dist = 20
# x = 70
# target_count = 100
# min_bp_dist = 20
# max_bp_dist = 30
# x = 100
# target_count = 100
# min_bp_dist = 40
# max_bp_dist = 60
# x = 150


target_count = 1000
min_bp_dist = 10
max_bp_dist = float("inf")
x = 200

Verbose = False

s_list = []
results = []


def find_moves(s_ptable, t_ptable):
    """
    generator function, yields possible structures 1 move away
    from the original structure by finding fitting i and j with
    RNA pair and loop tables
    s_ptable: current ptable
    t_ptable: s2 end ptable
    """

    ls = RNA.loopidx_from_ptable(s_ptable)

    for i in range(len(s_ptable)):
        if i == 0:
            continue

        if s_ptable[i] == 0 and t_ptable[i] > i:
            j = t_ptable[i]
            # found j has to be empty and currently on the same loop as i
            if s_ptable[j] == 0 and ls[i]==ls[j]:
                return i, j
        # test for bp removal: i has to be paired with a different j in s2
        j = s_ptable[i]
        # dont remove things which are present in s2
        if s_ptable[i] > i and s_ptable[i] != s_ptable[j] and\
                s_ptable[i] != t_ptable[i] and s_ptable[j] != t_ptable[j]:
            return -i, -j


# check confirmations for merging

def detect_local_minimum(fc, structure):
    # perform gradient walk from sample to determine direct local minimum
    pt = RNA.IntVector(RNA.ptable(structure))
    fc.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
    return RNA.db_from_ptable(list(pt))



useless_sections = 0

while (len(results) < target_count):
    
    # RNAsubopt wrapper
    sequence, s1, s2 = helper.generate_structures(x)

    # overwrite s1 s2 to local minimum
    fc = RNA.fold_compound(sequence)
    s1 = detect_local_minimum(fc, s1)
    s2 = detect_local_minimum(fc, s2)

    bp_dist = RNA.bp_distance(s1, s2)
    if bp_dist < min_bp_dist:
        continue
    if bp_dist > max_bp_dist:
        continue
    l = merge_composition.merge_check(sequence, s1, s2, Debug=False)

    # ignore the case where sections provide no benefit
    if len(l) == 1 or l==[0,x]:
        useless_sections += 1
        # continue


    pt1 = RNA.ptable(s1)
    pt2 = RNA.ptable(s2)
    
    i_move, j_move = find_moves(pt1, pt2)

    # print (sequence)
    # print (s1)
    # print (s2)
    # print (l, len(s_list), len(l), bp_dist)
    # print (i_move, j_move)

    # break
    results.append((sequence, s1, i_move, j_move))


df = pd.DataFrame(results, columns=['sequence','s1', 'i_move', 'j_move']).set_index('sequence') 
filename = f'../sample_seqs/eval_move_dataset.csv'
df.to_csv(filename)

print (df)
