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

from sklearn.cluster import DBSCAN

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

def next_vicinity(lastmove, avail_moves):

    """
    percentage of moves which can continue in the vicinity of the last move
    """
    lasti, lastj = lastmove
    moves = len(avail_moves)
    cntr = 0

    for a in avail_moves:
        i=a[0]
        j=a[1]      
        if lasti + 1 == i or lasti - 1 == i or\
           lastj + 1 == j or lastj - 1 == j:
           cntr += 1
    return cntr

def last_vs_best(lastmove, bestmove, lensequence):
    lasti = abs(lastmove[0])
    lastj = abs(lastmove[1])
    i = abs(bestmove[0])
    j = abs(bestmove[1])
    dist = abs(lasti - i) + abs(lastj - j)
    return dist/(lensequence*2) #, dist, lensequence


def cluster_moves(inputmoves, absolute_values=False):
    """
    cluster moves with DBscan. needs to include the last move as well...
    """

    if absolute_values:
        m = [(abs(i[0]), abs(i[1])) for i in inputmoves]
        moves = np.array(m)
    else:
        m = [(i[0], i[1]) for i in inputmoves]
        moves = np.array(m)

    # print ("clustering", moves)
    clustering = DBSCAN(eps=2, min_samples=1).fit(moves)
    
    # generate clustering dict (move-tuple:label, ...)
    return {i[0]:i[1] for i in zip(m, clustering.labels_)}, clustering.labels_


def structure_evaluation(fc, pt1, pt2, path, lastmove):
    avail_moves = []
    fp_results = []
    en_contrib = []

    available_add = set()
    available_delete = set()
    lensequence = pt1[0]

    # pt1 = list(RNA.ptable(s))
    # pt2 = list(RNA.ptable(s2))
    for pos, (i,j) in enumerate(find_moves(pt1, pt2)):    
        next_en = fc.eval_move_pt(pt1, i, j)
        # mark where we found our move    
        # map energies somehow to [0,1] such that the network will understand.  
        en = np.interp(next_en/100, [-10,10], [0,1]) 
        en_contrib.append(en)
        avail_moves.append((i, j, next_en))

        if i>0:
            available_add.add((i,j))
        elif i<0:
            available_delete.add((i,j))

    if len(avail_moves)==0:
        add_delete = 0
    else:
        add_delete = len(available_delete)/len(avail_moves)

    # print ('best', en_contrib[0], 'avg', en_mean)

    avail_moves.sort(key=lambda x: x[2]) # best move at [0]


    vic = next_vicinity(lastmove, avail_moves)
    vic_best = next_vicinity(lastmove, avail_moves[0:1])

    avail_moves = [(i[0], i[1]) for i in avail_moves]
    avail_moves_abs = [(abs(i[0]), abs(i[1])) for i in avail_moves] # this is sorted as well, best move at [0]

    # absolute_moves+=[lastmove]
    # print ("avail+last", avail_moves)

    # last_index = -1
    # best_index = absolute_moves.index((abs(absolute_moves[0][0]), abs(absolute_moves[0][1])))
    # print ("best", best_index)

    p = [(0, 0)] + [(i[0], i[1]) for i in path]
    clustered_moves, _ = cluster_moves(p, absolute_values=False)
    clustered_moves_abs, _ = cluster_moves(p, absolute_values=True)

    bestmove = avail_moves[0]
    bestmove_abs = avail_moves_abs[0]

    lastmove_abs = abs(lastmove[0]), abs(lastmove[1])

    # last move vicinity is available, but not among the best energy options:

    unique_found = set()
    foundpos = -1
    for pos, m in enumerate(avail_moves):
        l = clustered_moves[m]
        if l == clustered_moves[lastmove]:
            foundpos = len(unique_found) # cant add pos here

        if l not in unique_found:
            unique_found.add(l)
    if foundpos == -1:
        foundpos = len(unique_found)
    foundpos /= len(unique_found)

    unique_found_abs = set()
    foundpos_abs = -1
    for pos, m in enumerate(avail_moves_abs):
        l = clustered_moves_abs[m]
        if l == clustered_moves_abs[lastmove_abs]:
            foundpos_abs = len(unique_found_abs) # cant add pos here

        if l not in unique_found_abs:
            unique_found_abs.add(l)
    if foundpos_abs == -1:
        foundpos_abs = len(unique_found_abs)
    foundpos_abs /= len(unique_found_abs)

    # print (foundpos, last_is_best, last_is_best_abs)
    # print ("bestmove", bestmove, clustered_moves[bestmove], bestmove_abs, clustered_moves_abs[bestmove_abs])
    # print ("lastmove", lastmove, clustered_moves[lastmove], lastmove_abs, clustered_moves_abs[lastmove_abs])

    # distlast = last_vs_best(clustered_moves)

    unique_labels = 1 / len(unique_found)
    # unique_moves = 1/unique_labels

    distlast = last_vs_best(lastmove, bestmove, lensequence)
    unique_moves = 0

    # print ("unique", unique_labels, unique_moves)

    count_available = len(avail_moves)/len(path)


    for pos, (i,j) in enumerate(avail_moves):  
        pt = adjust_pt(pt1, i, j)
        snew = RNA.db_from_ptable(pt)
        # if i > 0:
        #     snew = s[:i-1] + "(" + s[i:j-1] + ")" + s[j:]
        # if i < 0:
        #     snew = s[:-i-1] + "." + s[-i:-j-1] + "." + s[-j:]

        ptnew = list(RNA.ptable(snew))
        result_new, path = fp_call(sequence, snew, s2, search_width_multiplier)
        fp_results.append(result_new)

    if len(avail_moves) == 0:
        label = 0
        en_mean = 0
        en_std = 0
        best_en = 0

    else:
        label = np.argmin(fp_results)/len(fp_results)

        en_mean = np.mean(en_contrib)
        en_std = np.std(en_contrib)
        en_contrib.sort(key=lambda x: x)
        en_diff = en_contrib[0] / en_mean # best en, relative

        best_en = en_contrib[0]

   
    return label, en_mean, en_std, best_en, foundpos, foundpos_abs, distlast, unique_labels



filename_samples = f'./dataset_100_large.csv'
samples_df = pd.read_csv(filename_samples)
sequence = ''
s1 = ''
s2 = ''
search_width_multiplier = 6

df = pd.DataFrame()
all_results = []



for index, row in samples_df.iterrows():
    # if index != 90:
    #     continue

    sequence = row.sequence
    s1 = row.s1
    s2 = row.s2


    bp_dist = RNA.bp_distance(s1, s2)

    # print (index)
    # print (sequence)
    # print (s1)
    # print (s2)
    # print (bp_dist)

    # random_move = random.choice(range(bp_dist))

    fc = RNA.fold_compound(sequence)
    pt1 = list(RNA.ptable(s1))
    pt2 = list(RNA.ptable(s2))
    # structure_evaluation(sequence, pt1, pt2)

    # sE, path = fp_call(sequence, s1, s2, search_width_multiplier)
    fp = findpath.findpath_single(sequence, s1, s2, search_width_multiplier=search_width_multiplier, mp=False)
    sE = fp.get_en()/100.0
    path = fp.get_path()

    # print (path)

    p = [(0, 0)] + [(i[0], i[1]) for i in path]
    clustered_moves, _ = cluster_moves(p, absolute_values=False)
    clustered_moves_abs, _ = cluster_moves(p, absolute_values=True)

    # print (clustered_moves)

    random_moves = random.choices([i[0] for i in path], k=3)
    # print (random_moves)
    # print (len(path))

    ptlast = list(RNA.ptable(s1))
    pt = list(RNA.ptable(s1))

    lastmove = 0, 0

    all_evals = 0
    results = []

    # find a random structure in the path
    # ignore the last moves, we're not interested in meaningless decisions
    for i, j, en in path:
        if i==0: 
            continue
        # if i not in random_moves:
        #     continue

        # print (i, RNA.db_from_ptable(pt))

        if i in random_moves:
            label, en_mean, en_std, best_en, foundpos, foundpos_abs, distlast, unique_labels = structure_evaluation(fc, pt, pt2, path, lastmove)
            

            # pred label with old data
            # data = {'en_mean': [en_mean],
            # 'en_std': en_std,
            # 'best_en': best_en,
            # # 'vic': vic,
            # # 'vic_best': vic_best,
            # # 'unique_moves': unique_moves,
            # # 'distlast': distlast
            # }
            # test_features = pd.DataFrame.from_dict(data)  
            # y = reloaded_model.predict(test_features)[0][0]
            # pred_label = round(y,4)
            
            
            # results.append((sequence, RNA.db_from_ptable(pt), s2, i, pred_label, label, en_mean, en_std, best_en, foundpos, foundpos_abs, distlast, unique_labels))
            results.append((sequence, RNA.db_from_ptable(pt), s2, i, j, label, en_mean, en_std, best_en, foundpos, foundpos_abs, distlast, unique_labels))

            all_evals += label       
            # print (i, eval)  

        pt = adjust_pt(ptlast, i, j)

        lastmove = i, j
        ptlast = pt

    if all_evals != 0:
        all_results += results


    # break

    if index%10 == 0:        
        print (index, end=" ")
        # print (all_results)
        df = pd.DataFrame(all_results, columns=['sequence', 's', 's2', 'i', 'j', 'target', 'en_mean', 'en_std', 'best_en', 'foundpos', 'foundpos_abs', 'distlast', 'unique_labels'])# .set_index('sequence')
        filename = f'./move_diversity4.csv'
        df.to_csv(filename)

    if index == 800:
        break

print (df)


