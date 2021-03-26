#!/usr/bin/env python3
# coding: utf-8

import pathfinder
import time
import os
import subprocess
import concurrent.futures
import pandas as pd

import RNA

o_filename = r"local_min_300_multiple_sections_min10.csv"
# o_filename = r"local_min_400_multiple_sections_min10.csv"
# o_filename = r"local_min_500_multiple_sections_min10.csv"

filename = r"./sample_seqs/" + o_filename

single_1_runtimes = []
single_1_results = []

single_2_runtimes = []
single_2_results = []

mp_runtimes = []
mp_results = []

cpp_orig_runtimes = []
cpp_orig_results = []

py_runtimes = []
py_results = []


df = pd.read_csv(filename)
elements = len(df.index)


print ("processing:", filename)

def launch_fp(findpath_cmd, sequence, s1, s2, search_width):
    start_findpath = time.time()
    cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width} '
    # cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width} -v'
    DEVNULL = open(os.devnull, 'wb')
    result = subprocess.check_output(cmd, shell=True, encoding="utf8", stderr=DEVNULL)    
           
    time_fp = round(time.time()-start_findpath,4)
    
    result = float(result)
    
    return time_fp, result


def launch_fp_mp(findpath_cmd, sequence, s1, s2, search_width):
    start_findpath = time.time()
    
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        result = executor.map(launch_fp, [findpath_cmd, findpath_cmd], [sequence, sequence], [s1, s2], [s2, s1], [search_width, search_width])
    # result = executor.map(pathfinder_cpp.find_path, [sequence, sequence, sequence, sequence], [s1, s2, s1, s2], [s2, s1, s2, s1], [search_width, search_width,search_width, search_width])
    # result = executor.map(pathfinder_cpp.find_path, [sequence, sequence, sequence, sequence, sequence, sequence], [s1, s2, s1, s2, s1, s2], [s2, s1, s2, s1, s2, s1], [search_width, search_width, search_width, search_width, search_width, search_width])
    result = list(result)   
    
    result = min(list(result)[0][1],list(result)[1][1])
    
    time_fp = round(time.time()-start_findpath,4)
    
    return time_fp, result


for index, row in df.iterrows():

    # if index>=2:
    #     break
    # if index != 7:
    #     continue       
    percent_complete = 100-(elements-index)/elements*100       
    bar_length = 20
    filled_length = int(percent_complete/100*bar_length)
    rest = bar_length - filled_length
    bar = "█" * filled_length + '_' * rest

    print(f'\rComputing |{bar}| {percent_complete:.1f}% complete {index} ', end = "\r")
    
    sequence = row.sequence
    s1 = row.s1
    s2 = row.s2
    
    search_width = RNA.bp_distance(s1, s2)*2

    
    time_fp, result = launch_fp("./fp_single_test", sequence, s1, s2, search_width)
#     print ("time elapsed single s1s2 thread:", time_fp, result)
    single_1_runtimes.append(time_fp)
    single_1_results.append(result)    
    
    time_fp, result = launch_fp("./fp_single_test", sequence, s2, s1, search_width)
#     print ("time elapsed single s2s1 thread:", time_fp, result)
    single_2_runtimes.append(time_fp)
    single_2_results.append(result)    

    time_fp, result = launch_fp("./fp_multi_test", sequence, s1, s2, search_width)
#     print ("time elapsed inbuilt mp:", time_fp, result) 
    mp_runtimes.append(time_fp)
    mp_results.append(result)   

    time_fp, result = launch_fp("./findpath_orig", sequence, s1, s2, search_width)
#     print ("time elapsed orig c++:", time_fp, result)
    cpp_orig_runtimes.append(time_fp)
    cpp_orig_results.append(result)   
    
    start_findpath3 = time.time()
    result = pathfinder.pathfinder(sequence, s1, s2, search_width=search_width, verbose=False).max_en
    time_fp = round(time.time()-start_findpath3,4)    
    py_runtimes.append(time_fp)
    py_results.append(result)   
    
  
data = [single_1_runtimes, single_1_results, single_2_runtimes, single_2_results, mp_runtimes, mp_results, cpp_orig_runtimes, cpp_orig_results, py_runtimes, py_results]  
  
df = pd.DataFrame(data)
df = df.transpose() 
df.columns = ["single_1_runtimes", "single_1_results", "single_2_runtimes", "single_2_results", "mp_runtimes", "mp_results", "cpp_orig_runtimes", "cpp_orig_results", "py_runtimes", "py_results"] 

savefile = r"./results/" + "_" + o_filename
df.to_csv(savefile)
print (df)