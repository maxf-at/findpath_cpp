#!/usr/bin/env python3
# coding: utf-8


import time
import os
import sys
import subprocess
import concurrent.futures
import pandas as pd

import RNA

sys.path.append('../')
import findpath_librna
import findpath




# o_filename = r"local_min_100_multiple_sections_min10.csv"
# o_filename = r"local_min_200_multiple_sections_min10.csv"
# o_filename = r"local_min_300_multiple_sections_min10.csv"
# o_filename = r"local_min_400_multiple_sections_min10.csv"
# o_filename = r"local_min_500_multiple_sections_min10.csv"
# o_filename = r"local_min_600_multiple_sections_min10.csv"
# o_filename = r"local_min_800_multiple_sections_min10.csv"
# o_filename = r"local_min_1000_multiple_sections_min10.csv"

# barriers direct
# o_filename = "5000x_100-20.csv"
# o_filename = "1000x_100-20.csv"
# o_filename = "1000x_150-25.csv"
# o_filename = "1000x_200-30.csv"

# w/o barriers
# o_filename = "local_min_100x1000_multiple_sections_min10.csv"
# o_filename = "local_min_200x1000_multiple_sections_min10.csv"
o_filename = "local_min_300x1000_multiple_sections_min10.csv"


filename = r"./sample_seqs/" + o_filename
df = pd.read_csv(filename)
elements = len(df.index)
print("processing:", filename)







# new api
def launch_new_fp(sequence, s1, s2, swm=None, mp=True):
    start_findpath = time.time()
    result = findpath.init_single_findpath(sequence, s1, s2, swm, mp)
    end_findpath = round(time.time()-start_findpath, 4)
    result = round(result/100.0,2)
    return end_findpath, result
def launch_new_merge_fp(sequence, s1, s2, swm=None, mp=True):
    start_findpath = time.time()
    fp = findpath.findpath_class(sequence, mp)
    fp.init(s1, s2, swm)
    result = fp.get_en()
    end_findpath = round(time.time()-start_findpath, 4)
    result = round(result/100.0,2)
    return end_findpath, result

def launch_new_merge_ext_fp(sequence, s1, s2, swm=None, mp=True):
    start_findpath = time.time()
    result = findpath.init_merge_ext_findpath(sequence, s1, s2, swm, mp)
    end_findpath = round(time.time()-start_findpath, 4)
    result = round(result/100.0,2)
    return end_findpath, result

def launch_new_merge_mfe_fp(sequence, s1, s2, swm=None, mp=True):
    start_findpath = time.time()
    result = findpath.init_mfe_findpath(sequence, s1, s2, swm, mp)
    end_findpath = round(time.time()-start_findpath, 4)
    result = round(result/100.0,2)
    return end_findpath, result

sequences = []
s1s = []
s2s = []
search_width_multipliers = []
indices = []
seq_lengths = []
new_fp_runtimes = []
new_fp_results = []
py_runtimes = []
py_results = []
new_merge_runtimes = []
new_merge_results = []
new_merge_ext_runtimes = []
new_merge_ext_results = []
new_merge_mfe_runtimes = []
new_merge_mfe_results = []
bp_dist = []
all_barriers_max_en = []
all_moves = []
all_indirect_moves = []
old_merge_results = []


for index, row in df.iterrows():


    percent_complete = 100-(elements-index)/elements*100
    bar_length = 20
    filled_length = int(percent_complete/100*bar_length)
    rest = bar_length - filled_length
    bar = "???" * filled_length + '_' * rest


    sequence = row.sequence
    s1 = row.s1
    s2 = row.s2

    # w/barriers
    # i = row.i
    # barriers_max_en = row.barriers_max_en
    # moves = row.moves
    # indirect_moves = row.indirect_moves

    # w/o barriers
    i = index
    barriers_max_en = 0
    moves = []
    indirect_moves = []

    print(
        f'\rComputing |{bar}| {percent_complete:.1f}% complete {i} ', end="\r")

    fc = RNA.fold_compound(sequence)
    s1_eval = round(fc.eval_structure(s1), 2)

    # search width scaling
    # sws = [1,1.5,2,3,4,6,8,10,12,14,16,18,20,30,40]
    # sws = [2]
    # sws = [1]

#     sws = [1,2,3,4,6,8,10,15,20]
    # sws = [1,2,3,4,6,8,10,15,20,30,40,60]
    # sws = [0.5, 0.75, 1,2,4,8,16,32] # up to 200
    # sws = [0.5, 0.75, 1,2,4,8] # up to 200
    sws = [0.5, 1,2,4,7,10] # up to 200
    # sws = [1,2,3,4,6]


    
    if index > 100: 
        break


    for search_width_multiplier in sws:

        # print(index, "current sw:", search_width_multiplier)

        indices.append(i)
        sequences.append(sequence)
        seq_lengths.append(len(sequence))
        s1s.append(s1)
        s2s.append(s2)
        search_width_multipliers.append(search_width_multiplier)
        bp_dist.append(RNA.bp_distance(s1, s2))
        search_width = int(RNA.bp_distance(s1, s2)*search_width_multiplier)

        # if search_width == 0:
            # print ("0!")
        # print (search_width_multiplier, search_width, RNA.bp_distance(s1, s2))
        # print ()


        all_barriers_max_en.append(barriers_max_en)
        all_moves.append(moves)
        all_indirect_moves.append(indirect_moves)
        
        # print ("new fp")
        time_fp, result = launch_new_fp(
            sequence, s2, s1, swm=search_width_multiplier, mp=True)
        new_fp_runtimes.append(time_fp)
        new_fp_results.append(round(result,2))

        # interior loop merge
        # print ("int loop merge")
        time_fp, result = launch_new_merge_fp(
            sequence, s1, s2, swm=search_width_multiplier, mp=True)
        new_merge_runtimes.append(time_fp)
        new_merge_results.append(round(result,2))

        # time_fp, result = launch_new_merge_ext_fp(
        #     sequence, s1, s2, swm=search_width_multiplier, mp=True)
        # new_merge_ext_runtimes.append(time_fp)
        # new_merge_ext_results.append(round(result-s1_eval,2))

        # mfe
        # print ("mfe")
        time_fp, result = launch_new_merge_mfe_fp(
            sequence, s1, s2, swm=search_width_multiplier, mp=True)
        new_merge_mfe_runtimes.append(time_fp)
        new_merge_mfe_results.append(round(result,2))

        # regular findpath
        # print ("regular findpath")
        start_findpath3 = time.time()
        result = findpath_librna.pathfinder(
            sequence, s1, s2, search_width=search_width, verbose=False)
        time_fp = round(time.time()-start_findpath3, 4)
        py_runtimes.append(time_fp)
        py_results.append(round(result,2))

   
        old_result = 0

        old_merge_results.append(old_result)



data = [indices, sequences, s1s, s2s, seq_lengths, search_width_multipliers, bp_dist, new_fp_runtimes, new_fp_results,
        py_runtimes, py_results, new_merge_runtimes, new_merge_results, new_merge_ext_runtimes,
        new_merge_ext_results, new_merge_mfe_runtimes, new_merge_mfe_results,
        all_barriers_max_en, all_moves, all_indirect_moves, old_merge_results]

df_save = pd.DataFrame(data)
df_save = df_save.transpose()
df_save.columns = ["i", "sequence", "s1", "s2", "seq_length", "search_width_multiplier", "bp_dist", "new_fp_runtimes", "new_fp_results",
        "py_runtimes", "py_results", "new_merge_runtimes", "new_merge_results", "new_merge_ext_runtimes",
        "new_merge_ext_results", "new_merge_mfe_runtimes", "new_merge_mfe_results",
        "barriers_max_en", "moves", "indirect_moves", "old_merge_results"]

# prefix = "sw_scaling_fasthash_"
# prefix = "sw_scaling_"
prefix = "sw_scaling_final_"


savefile = r"./results/" + prefix + o_filename
df_save.to_csv(savefile)
print(df_save)
