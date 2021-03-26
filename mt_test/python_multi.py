#!/usr/bin/env python3
# coding: utf-8

import time
import concurrent.futures
import subprocess

# launching ./single concurrently

start_findpath = time.time()

def single(dummy=""):
    cmd = "./single"
    result = subprocess.check_output(cmd, shell=True, encoding="utf8")    
    print (result)


# single()
time_1thread = round(time.time()-start_findpath,4)

start_findpath = time.time()

with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
    result = executor.map(single, ["./main", "./main"])


result = list(result)
time_2threads = round(time.time()-start_findpath,4)


# print ("1 thread: ", time_1thread)
print ("2 threads:", time_2threads)