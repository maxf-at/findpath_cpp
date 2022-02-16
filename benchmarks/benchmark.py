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


# sys.path.append('../pathfinder/src')
# import merge_recursive
# import merge_composition
# import pathfinder


# o_filename = r"local_min_100_multiple_sections_min10.csv"
# o_filename = r"local_min_200_multiple_sections_min10.csv"
# o_filename = r"local_min_300_multiple_sections_min10.csv"
# o_filename = r"local_min_400_multiple_sections_min10.csv"
# o_filename = r"local_min_500_multiple_sections_min10.csv"
# o_filename = r"local_min_600_multiple_sections_min10.csv"
# o_filename = r"local_min_800_multiple_sections_min10.csv"
# o_filename = r"local_min_1000_multiple_sections_min10.csv"


# o_filename = "local_min_100x1000_multiple_sections_min10.csv"
# o_filename = "local_min_200x1000_multiple_sections_min10.csv"
o_filename = "local_min_300x1000_multiple_sections_min10.csv"

filename = r"./sample_seqs/" + o_filename

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


for index, row in df.iterrows():

    # if index!=5:
    #     continue
    # if index>2:
    #     break
    # if index != 16:
    #     continue
    #

    percent_complete = 100-(elements-index)/elements*100
    bar_length = 20
    filled_length = int(percent_complete/100*bar_length)
    rest = bar_length - filled_length
    bar = "█" * filled_length + '_' * rest

    print(
        f'\rComputing |{bar}| {percent_complete:.1f}% complete {index} ', end="\r")

    sequence = row.sequence
    s1 = row.s1
    s2 = row.s2

    # print(f"sequence = \"{sequence}\"")
    # print(f"s1 = \"{s1}\"")
    # print(f"s2 = \"{s2}\"")

    fc = RNA.fold_compound(sequence)
    s1_eval = round(fc.eval_structure(s1), 2)

    # search width scaling
    # sws = [1,1.5,2,3,4,6,8,10] #,12,14,16,18,20,30,40]
    sws = [2]
    # sws = [1]

    for search_width_multiplier in sws:

        # print(index, "current sw:", search_width_multiplier)

        indices.append(index)
        sequences.append(sequence)
        seq_lengths.append(len(sequence))
        s1s.append(s1)
        s2s.append(s2)
        search_width_multipliers.append(search_width_multiplier)
        bp_dist.append(RNA.bp_distance(s1, s2))
        search_width = int(RNA.bp_distance(s1, s2)*search_width_multiplier)

        time_fp, result = launch_new_fp(
            sequence, s2, s1, swm=search_width_multiplier, mp=True)
        new_fp_runtimes.append(time_fp)
        new_fp_results.append(round(result-s1_eval,2))

        time_fp, result = launch_new_merge_fp(
            sequence, s1, s2, swm=search_width_multiplier, mp=True)
        new_merge_runtimes.append(time_fp)
        new_merge_results.append(round(result-s1_eval,2))

        # time_fp, result = launch_new_merge_ext_fp(
        #     sequence, s1, s2, swm=search_width_multiplier, mp=True)
        # new_merge_ext_runtimes.append(time_fp)
        # new_merge_ext_results.append(round(result-s1_eval,2))
        new_merge_ext_runtimes.append(0)
        new_merge_ext_results.append(0)


        time_fp, result = launch_new_merge_mfe_fp(
            sequence, s1, s2, swm=search_width_multiplier, mp=True)
        new_merge_mfe_runtimes.append(time_fp)
        new_merge_mfe_results.append(round(result-s1_eval,2))

        start_findpath3 = time.time()
        result = findpath_librna.pathfinder(
            sequence, s1, s2, search_width=search_width, verbose=False)
        time_fp = round(time.time()-start_findpath3, 4)
        py_runtimes.append(time_fp)
        py_results.append(round(result-s1_eval,2))


# new_fp_mp_runtimes = []
# new_fp_mp_results = []
# new_fp_runtimes = []
# new_fp_results = []
# cpp_orig_runtimes = []
# cpp_orig_results = []
# py_runtimes = []
# py_results = []

data = [indices, sequences, s1s, s2s, seq_lengths, search_width_multipliers, bp_dist, new_fp_runtimes, new_fp_results,
        py_runtimes, py_results, new_merge_runtimes, new_merge_results, new_merge_ext_runtimes,
        new_merge_ext_results, new_merge_mfe_runtimes, new_merge_mfe_results]

df = pd.DataFrame(data)
df = df.transpose()
df.columns = ["i", "sequence", "s1", "s2", "seq_length", "search_width_multiplier", "bp_dist", "new_fp_runtimes", "new_fp_results",
        "py_runtimes", "py_results", "new_merge_runtimes", "new_merge_results", "new_merge_ext_runtimes",
        "new_merge_ext_results", "new_merge_mfe_runtimes", "new_merge_mfe_results"]

prefix = "final2_"

savefile = r"./results/" + prefix + o_filename
df.to_csv(savefile)
print(df)

# print ("new fp")
# print (new_fp_results)
# print ("merge")
# print (new_merge_results)