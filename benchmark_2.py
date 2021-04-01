#!/usr/bin/env python3
# coding: utf-8

import pathfinder
import time
import os
import sys
import subprocess
import concurrent.futures
import pandas as pd

import RNA



sys.path.append('../pathfinder/src')

import merge_recursive
import merge_composition

# os.chdir("../../pathfinder_cpp")


# o_filename = r"local_min_100_multiple_sections_min10.csv"
o_filename = r"local_min_200_multiple_sections_min10.csv"
# o_filename = r"local_min_300_multiple_sections_min10.csv"
# o_filename = r"local_min_400_multiple_sections_min10.csv"
# o_filename = r"local_min_500_multiple_sections_min10.csv"

filename = r"./sample_seqs/" + o_filename

single_1_runtimes = []
single_1_results = []

single_2_runtimes = []
single_2_results = []

mp_runtimes = []
mp_results = []

new_merge_runtimes = []
new_merge_results = []

old_merge_runtimes = []
old_merge_results = []

cpp_orig_runtimes = []
cpp_orig_results = []

py_runtimes = []
py_results = []

bp_dist = []


df = pd.read_csv(filename)
elements = len(df.index)


print ("processing:", filename)

def launch_fp(findpath_cmd, sequence, s1, s2, search_width=None, search_width_multiplier=None, rprint=False):
    start_findpath = time.time()

    if search_width:
        cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width} '
    if search_width_multiplier:
        cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width_multiplier} '

    # cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width} '

    
    # ignore stderr
    DEVNULL = open(os.devnull, 'wb')
    # result = subprocess.check_output(cmd, shell=True, encoding="utf8", stderr=DEVNULL)    
    result = subprocess.check_output(cmd, shell=True, encoding="utf8")

    time_fp = round(time.time()-start_findpath,4)


    if (rprint):
        print (result)



    result = result.splitlines()[-1]
    # print (result)
    if result == "Graph object being created":
        result = 0
    else:
        result = float(result)
    
    return time_fp, result



# def launch_fp(sequence, s1, s2, findpath_cmd, search_width=None, search_width_multiplier=None):
#     start_findpath = time.time()
#     if search_width:
#         cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width} '
#     if search_width_multiplier:
#         cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width_multiplier} '

#     result1 = subprocess.check_output(cmd, shell=True, encoding="utf8")
#     end_findpath = round(time.time()-start_findpath,4)
#     return end_findpath, result1


def launch_fp_mp(findpath_cmd, sequence, s1, s2, search_width, rprint=False):
    start_findpath = time.time()
    
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        result = executor.map(launch_fp, [findpath_cmd, findpath_cmd], [sequence, sequence], [s1, s2], [s2, s1], [search_width, search_width])


    result = list(result)   
    
    result = min(list(result)[0][1],list(result)[1][1])
    
    time_fp = round(time.time()-start_findpath,4)
    
    return time_fp, result


for index, row in df.iterrows():

    # if index>=1:
    #     break

    #300.23
    if index != 23:
        continue    

    # if index != 17:
    #     continue    
    #    
    percent_complete = 100-(elements-index)/elements*100       
    bar_length = 20
    filled_length = int(percent_complete/100*bar_length)
    rest = bar_length - filled_length
    bar = "█" * filled_length + '_' * rest

    print(f'\rComputing |{bar}| {percent_complete:.1f}% complete {index} ', end = "\r")
    
    sequence = row.sequence
    s1 = row.s1
    s2 = row.s2

    print (f"sequence = \"{sequence}\"")
    print (f"s1 = \"{s1}\"")
    print (f"s2 = \"{s2}\"")
    
    bp_dist.append(RNA.bp_distance(s1, s2))

    search_width_multiplier = 20

    search_width = RNA.bp_distance(s1, s2)*search_width_multiplier
    # search_width = RNA.bp_distance(s1, s2)*2

    
    time_fp, result = launch_fp("./fp_single_test", sequence, s1, s2, search_width=search_width)
#     print ("time elapsed single s1s2 thread:", time_fp, result)
    single_1_runtimes.append(time_fp)
    single_1_results.append(result)    
    
    time_fp, result = launch_fp("./fp_single_test", sequence, s2, s1, search_width=search_width)
#     print ("time elapsed single s2s1 thread:", time_fp, result)
    single_2_runtimes.append(time_fp)
    single_2_results.append(result)    

    time_fp, result = launch_fp("./fp_multi_test", sequence, s1, s2, search_width=search_width)
#     print ("time elapsed inbuilt mp:", time_fp, result) 
    mp_runtimes.append(time_fp)
    mp_results.append(result)   

    time_fp, result = launch_fp("./main", sequence, s1, s2, search_width_multiplier=search_width_multiplier)
#     print ("time elapsed inbuilt mp:", time_fp, result) 
    new_merge_runtimes.append(time_fp)
    new_merge_results.append(result)   


    sections = merge_composition.merge_check(sequence, s1, s2, Debug=False)
    result = merge_recursive.recursive_merge(sequence, s1, s2, sections, search_width_multiplier=search_width_multiplier,
                        Verbose=False, Debug=False, plot_graph=False, new=True,
                        print_sections=True, max_rec_depth=99)
    old_merge_runtimes.append(result.runtime_merge)
    old_merge_results.append(result.max_en)   


    time_fp, result = launch_fp("./findpath_orig", sequence, s1, s2, search_width=search_width)
#     print ("time elapsed orig c++:", time_fp, result)
    cpp_orig_runtimes.append(time_fp)
    cpp_orig_results.append(result)   
    
    start_findpath3 = time.time()
    result = pathfinder.pathfinder(sequence, s1, s2, search_width=search_width, verbose=False).max_en
    time_fp = round(time.time()-start_findpath3,4)    
    py_runtimes.append(time_fp)
    py_results.append(result)   
    
  
data = [single_1_runtimes, single_1_results, single_2_runtimes, single_2_results, mp_runtimes, mp_results, cpp_orig_runtimes, cpp_orig_results, new_merge_runtimes, new_merge_results, old_merge_runtimes, old_merge_results, py_runtimes, py_results, bp_dist]  
  
df = pd.DataFrame(data)
df = df.transpose() 
df.columns = ["single_1_runtimes", "single_1_results", "single_2_runtimes", "single_2_results", "mp_runtimes", "mp_results", "cpp_orig_runtimes", "cpp_orig_results", "new_merge_runtimes", "new_merge_results", "old_merge_runtimes", "old_merge_results", "py_runtimes", "py_results", "bp_dist"] 

savefile = r"./results/" + "3" + o_filename
df.to_csv(savefile)
print (df)