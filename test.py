#!/usr/bin/env python3
# coding: utf-8

import time
import subprocess
import RNA
import pandas as pd

import findpath
import findpath_librna





# o_filename = r"local_min_100_multiple_sections_min10.csv"
# o_filename = r"local_min_200_multiple_sections_min10.csv"
o_filename = r"local_min_300_multiple_sections_min10.csv"
# o_filename = r"local_min_400_multiple_sections_min10.csv"
# o_filename = r"local_min_500_multiple_sections_min10.csv"
# o_filename = r"local_min_600_multiple_sections_min10.csv"
# o_filename = r"local_min_800_multiple_sections_min10.csv"

filename = r"./benchmarks/sample_seqs/" + o_filename

df = pd.read_csv(filename)
elements = len(df.index)
print("processing:", filename)


search_width_multiplier = 2
reference_fp_results = [3.1, 6.6, 9.6, 2.5, 15.3, 13.4, 15.3, 1.8, 23.5, 6.4, 15.7, 1.8, 14.0, 20.7, 17.3, 14.2, 14.0, 11.7, 9.3, 3.5, 1.7, 11.9, 6.5, 15.66, 3.1, 10.3, 7.2, 24.0, 11.3, 11.3, 6.5, 8.3, 14.7, 3.8, 11.6, 8.6, 2.6, 7.7, 30.7, 7.7, 10.3, 6.3, 16.3, 5.7, 5.6, 12.0, 10.6, 14.1, 10.4, 0.6, 8.9, 9.7, 12.4, 15.1, 8.0, 8.1, 9.2, 2.0, 6.7, 5.5, 18.7, 13.3, 5.4, 8.1, 9.7, 26.2, 11.2, 7.1, 12.5, 19.5, 16.6, 4.8, 12.1, 10.98, 16.8, 10.0, 4.0, 13.2, 24.7, 24.8, 18.1, 26.81, 1.0, 7.4, 6.2, 17.4, 1.5, 4.3, 9.9, 5.9, 9.93, 23.2, 13.2, 25.5, 13.73, 6.9, 16.8, 22.0, 10.3, 15.9]
reference_merge_results = [2.3, 6.6, 9.6, 2.5, 14.7, 13.4, 15.3, 1.8, 23.5, 6.4, 15.7, 1.8, 14.0, 20.7, 17.3, 14.2, 12.8, 10.5, 9.3, 3.5, 1.7, 11.9, 6.5, 15.66, 3.1, 8.8, 7.2, 24.0, 11.3, 10.6, 6.3, 8.3, 14.7, 3.8, 11.6, 8.6, 2.6, 7.0, 30.7, 7.6, 10.1, 6.3, 15.6, 5.7, 5.6, 12.0, 9.9, 14.1, 9.5, 0.6, 8.9, 9.7, 12.1, 15.1, 8.0, 8.0, 9.2, 2.0, 6.7, 5.5, 18.7, 13.3, 5.4, 8.1, 9.7, 26.2, 10.6, 7.1, 9.8, 17.9, 16.0, 4.8, 12.1, 10.98, 16.8, 10.0, 4.0, 13.2, 24.7, 24.8, 18.3, 26.81, 1.0, 6.6, 6.2, 17.4, 1.5, 4.3, 9.9, 5.9, 9.2, 23.2, 13.2, 25.5, 13.73, 6.9, 15.0, 22.0, 10.3, 15.9]
reference_merge_mfe_results = []

#sw 1
search_width_multiplier = 1
if o_filename == "local_min_200_multiple_sections_min10.csv":
    reference_fp_results = [3.1, 6.6, 9.6, 2.5, 15.3, 13.4, 16.7, 1.8, 24.0, 6.4, 17.2, 1.8, 14.0, 21.0, 18.7, 14.3, 16.1, 10.5, 9.4, 3.5, 1.7, 11.9, 6.5, 15.66, 3.1, 11.5, 7.2, 26.6, 11.3, 11.3, 6.5, 8.3, 14.7, 3.8, 12.9, 9.3, 2.6, 7.7, 31.9, 7.8, 10.3, 6.3, 16.3, 5.7, 5.6, 12.0, 10.8, 15.1, 10.4, 0.6, 8.9, 9.7, 12.4, 15.7, 8.0, 8.6, 9.2, 2.0, 6.7, 5.5, 17.8, 13.7, 5.4, 8.1, 9.7, 26.2, 11.6, 7.1, 14.0, 19.5, 16.6, 4.8, 12.2, 10.98, 16.8, 10.5, 4.0, 13.2, 24.7, 24.8, 18.4, 30.71, 1.0, 7.9, 6.2, 19.5, 1.5, 4.3, 9.9, 5.9, 9.93, 24.4, 13.2, 25.5, 13.73, 6.9, 16.8, 23.2, 10.3, 15.9]
    reference_merge_results = [2.3, 6.6, 9.6, 2.5, 14.7, 13.4, 16.7, 1.8, 24.0, 6.4, 17.2, 1.8, 14.0, 21.0, 18.7, 14.2, 16.7, 10.5, 9.3, 3.5, 1.7, 11.9, 6.5, 15.66, 3.1, 8.8, 7.2, 26.6, 11.3, 10.6, 6.3, 8.3, 14.7, 3.8, 12.9, 8.6, 2.6, 7.0, 31.9, 7.6, 10.1, 6.3, 15.6, 5.7, 5.6, 12.0, 9.9, 15.1, 10.4, 0.6, 8.9, 9.7, 12.1, 15.7, 8.0, 8.0, 9.2, 2.0, 6.7, 5.5, 17.8, 13.9, 5.4, 8.1, 9.7, 26.2, 10.9, 7.1, 11.6, 17.9, 16.0, 4.8, 12.2, 10.98, 16.8, 10.5, 4.0, 13.2, 24.7, 24.8, 18.5, 30.71, 1.0, 6.6, 6.2, 19.5, 1.5, 4.3, 9.9, 5.9, 9.2, 24.4, 13.2, 25.5, 13.9, 6.9, 16.0, 23.2, 10.3, 15.9]
if o_filename == "local_min_300_multiple_sections_min10.csv":
    reference_fp_results = [7.8, 12.6, 15.7, 7.7, 25.1, 5.3, 3.0, 9.1, 9.6, 3.6, 3.4, 13.23, 27.8, 8.6, 12.8, 11.6, 14.3, 24.0, 13.1, 7.3, 8.5, 11.83, 7.6, 32.4, 17.0, 11.4, 12.3, 11.5, 1.4, 10.9, 19.4, 17.7, 16.1, 11.8, 11.3, 17.0, 8.2, 22.4, 17.6, 16.9, 23.2, 9.8, 29.0, 1.2, 10.8, 8.0, 6.2, 16.5, 18.7, 1.1, 2.3, 15.14, 6.6, 29.5, 6.5, 6.2, 5.2, 25.6, 8.6, 6.3, 15.0, 7.4, 4.3, 9.6, 24.8, 17.4, 42.6, 7.8, 8.5, 3.1, 25.63, 13.3, 29.4, 23.4, 10.1, 8.4, 15.8, 16.6, 3.1, 9.5, 4.6, 4.6, 17.8, 26.7, 6.6, 22.3, 32.1, 4.8, 21.1, 22.3, 1.5, 24.1, 17.9, 19.5, 3.8, 17.7, 3.7, 14.3, 26.5, 14.6]
    reference_merge_results = [7.8, 12.6, 14.8, 7.7, 23.8, 4.7, 3.0, 9.1, 9.6, 3.6, 3.7, 13.23, 27.8, 8.5, 9.6, 10.5, 10.9, 23.4, 13.1, 6.7, 8.1, 11.33, 7.1, 32.3, 13.0, 6.4, 12.2, 11.4, 1.4, 9.6, 19.4, 17.4, 14.8, 11.7, 9.1, 17.0, 8.2, 22.4, 16.1, 16.0, 17.4, 6.7, 29.0, 1.2, 10.8, 5.8, 6.2, 16.5, 15.4, 1.1, 2.3, 14.1, 5.5, 28.4, 6.5, 6.2, 4.2, 24.3, 8.6, 6.3, 15.0, 7.4, 4.3, 9.6, 25.96, 17.4, 42.6, 6.3, 8.3, 3.1, 22.73, 13.3, 31.9, 25.1, 10.1, 8.4, 15.8, 16.3, 3.1, 5.8, 4.6, 4.6, 14.9, 26.0, 6.6, 18.4, 34.9, 4.8, 21.1, 22.3, 1.5, 22.8, 17.9, 19.5, 3.8, 17.7, 3.7, 13.2, 26.5, 14.7]



reference_merge_mfe_results = []


start = time.time()


def launch_new_fp(sequence, s1, s2, swm=None, mp=True):
    start_findpath = time.time()
    result = findpath.init_single_findpath(sequence, s1, s2, swm, mp)
    end_findpath = round(time.time()-start_findpath, 4)
    result = round(result/100.0,2)
    return end_findpath, result
def launch_new_merge_fp(sequence, s1, s2, swm=None, mp=True):
    start_findpath = time.time()
    fp = findpath.findpath_class(sequence, mp)
    result = fp.init(s1, s2, swm)
    # result = findpath.init_merge_findpath(sequence, s1, s2, swm, mp)
    end_findpath = round(time.time()-start_findpath, 4)
    result = round(result/100.0,2)

    path = fp.return_path()
    print (path) 
    
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

all_offset_values = 0

for index, row in df.iterrows():

    # if index!=5:
    #     continue

    sequence = row.sequence
    s1 = row.s1
    s2 = row.s2
    fc = RNA.fold_compound(sequence)
    s1_eval = round(fc.eval_structure(s1), 2)
    
    print(index, "current sw:", search_width_multiplier)

    time_fp, result = launch_new_fp(
        sequence, s2, s1, swm=search_width_multiplier, mp=True)

    result = round(result-s1_eval,2)
    ref_result = reference_fp_results[index]
    offset = result - ref_result
    print (index, "reference result:", ref_result, "now:", result, "offset:", offset)
    all_offset_values += offset
    # if offset != 0.0:
        # print(f"sequence = \"{sequence}\"")
        # print(f"s1 = \"{s1}\"")
        # print(f"s2 = \"{s2}\"")
    #     break

    time_fp, result = launch_new_merge_fp(
        sequence, s2, s1, swm=search_width_multiplier, mp=True)

    result = round(result-s1_eval,2)
    ref_result = reference_merge_results[index]
    offset = result - ref_result
    print (index, "reference result:", ref_result, "now:", result, "offset:", offset)
    # if offset != 0.0:
        # print(f"sequence = \"{sequence}\"")
        # print(f"s1 = \"{s1}\"")
        # print(f"s2 = \"{s2}\"")
    #     break
    all_offset_values += offset

end = round(time.time()-start, 4)

print ("all deviations:", all_offset_values, "time:", end)







