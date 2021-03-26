#!/usr/bin/env python3
# coding: utf-8

import time
import subprocess
import numpy as np
import concurrent.futures


import RNA

import pathfinder




# sequence = "AUGAACAGAUGGUACCUCGCGGCGGGGCCUGACCCAUCUACAGUUUGUGCAG"
# s1 =       "((((((((((((..(((((...)))))......))))))...))))))...."
# s2 =       ".((..((((..(((((....))..(((.....)))...)))..))))..))."

# sequence = "UGGGUCCACGCAAGCCGGCACGCGAAAGGGCGGCUAUGUACAGCUCUUAAAACCACCAGAGGUUUAGUGAUCACUGAGGCUUGUUCGCAAAUCACUGCAAUUAGAUAUGACUCACGAUAUGGGGCACGGUGCAUACAUACAAGAAUUACUUGGCACCCCCGCGGAAAGUAGUCUCCUUGUACGCAUAUCGCUACAAAAUCGGGGAACUUCCAUACAUCCUGUUAUGUGAAAUGCGACGUAGAAGCCACACACUCGUCAUUGACAGAUCUACGAAUUAUGGGUCUACUAAGAAAUUUCACAUUCCGUAGAUUCUUGUACAAGGAUAGUGACUCCGCAAUUAGGACUGAACAAUAACUUGGGCGUAAAGUAAGCGCAGGACUAGAGGAAGGCUGGUGCGCGA"
# s1 =       ".(.(((((.(..(((((.(.(......)).))))).(((.(((.((((((((((......))))..(((((((((((.((.((.........))..))..))))...))).))))......(((...(((((.......((((.....))))))))))))(((((...((.(.(((((((((((((.((((((((..((..(((....)))....))..))).).)))).))))..((((((..........(((((...))).))))))))......((((((((...(((........))).))))))))..))))))))).).))..)))))..))))))))).)))....).))))).)........(((.(.(((((.......))))).))))."
# s2 =       "((.(((((.(..((((.((.(.......))))))).(((.(((.((((((((((......))))).((((((((((..((................))...)))...))).))))......(((...(((((....................))))))))(((((...(..(.(((((((((((((.((((.(((......(((....)))........)))...)))).))))((((.((..........)))))).......((((((((....((((((((....)))....)).)))..))))))))...))))))))).)...).)))))...)))))))).)))....).))))).))......((((...(((((.......))))))))).."

# this one is better with findpath_cpp
# sequence = "GCCGCUAGAUCGAGAUGUUUAGGAUAUGCAAACAGGUCUGGGGGGUCAUGUCUUUAAAGGUCAUCAUGAGGCUCCGCCAAGACUAGAGUCCGGCCGAGCCUACCACGGCUCAUCAUGGCGCACCCUCUUAUUUUGUCAACUUGCUCUUAUCCGAACGACAACUAUCACUUGGUUACUCUGCAAUGAUUACGUACCUAGUUGGGGACGACAAUUGCUCAAUCGUUCGCGCCGUACGUGUGGCGAAGGCCUUAUGAUGAUCAAAGAGACUGGCCAGCGCCCGGAAGUACGCCGAACGGAGCG"
# s1 =       ".((.(((((((....(((((.(......)))))))))))))))(((((.((((((...(((((((((((((((.(((((........((.((((.(((((......)))))((...((((((...........))......))))...)).(((((((.(((((.....)))))..((......))...(((.(((....))).)))............))))))).)))).))...)))))..)))))))))))))))..)))))))))))..((((((.............))).)))"
# s2 =       ".((.(((((((....(((.........)))....)))))))))(((((.((((((...(((((((((((((((.((((((((..((((....((((((((......))))).....))).....))))..))).((((..((((......((((..((......))..)))).......))))))))(((((((...(.((.((((((...........)))))).))).))))))))))))..)))))))))))))))..)))))))))))..((((((.............))).)))"

# 300 nt, 190 bp_dist, 2 sections
sequence = "CAAUUGCGCCGCAUAAAACCAGAGAUUCUACCCUCAAUCGGUUCUUAAGACGUACUGCGCGUUUCACCAGACCACAAUGCAGGGCGGCACCGUUAGGCAACACAACGAGACUACUCAUGCACAUAAGGAAGGUUAUCGCCAUAGACAUGGCGCGGCAGCGCAGAAUGUUUAAAUCUAAAUCUGGUAUGGGAGGCGUGCCCGUUGGUAUGAAGAAAUUUGCUGGGAGAAAAAGUCUAAGGCCUUGAAUCCGGCGGGUCUUAAUACUUACCUACAAAAUCAUCAGGCUGUACUUCCUGUAUC"
# s1 =       "......(((((.......((..(((((((.(((.(((....(((((....(((((((((((((((.(((.((((...(((..((((....))))..))).......(((....)))...........((......))....(((((((.((((....))))...)))))))..........)))).))))))))))).....))))))))))))..)))..))))))....))))..)).........)))))(((..........)))...........((((........))))...."
# s2 =       "....((((((...........(((........)))....(((..(((((((.....(((((((((.(((.((((........(((((.(((...........((..(((....))).))...........)))..))))).(((((((.((((....))))...)))))))..........)))).))))))))))))(((((((..(.(((......(((.(((.......)))..)))))).)..)))))))))))))).....))).............))).)))..........."


# 600 nt, 2 sections
# sequence = "UUGUGUGAUGUGCAUCGGUAGCUCGGUAAAUACGUUCGCCGGGCUUGAAGCUGUAUCCGUUGGAAUCCAAGUUGGGAGUUAUAACGUAAAACGGCAUUCCCGAUAGUCUCUCAGGGGACGAUCGAUGUUCAGGCGGCUGCGAUUUUGAACGCCGUGAAGCGUGAUUAAACCGGUCGACUUUUGCUCCUUGAAGUUUUCAUCUAUCACCUGUACUUGCGGCCACCACCUCAUUUGCCCGUUGGCCGACUCAUCGCGUGGUAGGGGUUGUCGACUGUGGGAGUUUUCAGAUUACGUCCUAACCAACCGUAAAAGAACCUCCGGUAAGAGUAAAAGGAGGGUACAGGGUACAGUAAAGACGAUAGGACAGCCGUGUGUUAGGGCCGGAAGCCGCGGGCAAAUUGACGGGACAGCACCUAGACAUACUCAAAUUAGGCCUAGCCACUAUUUGGCCGCACGCUUGUAUCAGAUUGGUCAAAGUACACUCCAUUUCUAAAGAUUUACCCCUCUAGCCCCUCCAGUGAGUGGCGCUGUCCGUUCAUAAAACUCCGUACACUCAUAUAGUCGAUGGCGAACGCGAUGUGGUGUAACCGAUUUACGUUA"
# s1 =       ".....(((((((.((((((((((((((..........)))))))).((((((((...(((((...(((......)))....)))))....))))).))).((((.(((((....))))).))))....(((((((....((....))..)))).)))((.((((.((.....(..((((((........))))))..)...)))))))).......((((((.(.............).))))))(((((((((((.......(((((((((((((.((((.......(((((...........))))).....(((((.............)))))((((.((....(((......((..(((((((.((((((((((((.(....((((((.....))).)))..).)).)))).)))))).(((((..(((...)))...)))))(((((.((((..........(((............)))...(((.((........)).))).......)))).))))))))))))..)).....))))))))))))).)))))))))))))..)))))))).)))...)))))).)))))))"
# s2 =       ".....(((((((.((((((((((((((..........))))))))....((((((((((((((..(((......)))......((((...(((((.....((((.(((((....))))).))))..(((((((((....))..))))))))))))...))))........((((((((....(((......))).....((((((((.........((((((.(.............).)))))).........))))))))....))))))))(((((((.........)).))))).)))))..........(((((.............))))).....)))))))))..........(((((((.(((((((..((((((((((.((((.(((((.....((......(((((...........))))).....))...))))).))))..))))........))))))..)))))))...........(((........))).(((.(((....))).))))))))))...............(((((.(((...((((.((.....)))))))))))))))))))).)))))))"

# random 600 nt example
sequence = "AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAAUUUUAUUAAAAGUCCAAGUUGGACUGACAAAACGCGUGCGGUGUCCUAGGGAUUGGUGGCAUAACCAGCGGUUUAAAAGCUGUGUAUAUCCGCAGCAAAUCACCGGAAAGCGGCGUUAUUAGCACCACAAAUUGAUGGUUGGUACGAGUACAAUUGCGCCGCAUAAAACCAGAGAUUCUACCCUCAAUCGGUUCUUAAGACGUACUGCGCGUUUCACCAGACCACAAUGCAGGGCGGCACCGUUAGGCAACACAACGAGACUACUCAUGCACAUAAGGAAGGUUAUCGCCAUAGACAUGGCGCGGCAGCGCAGAAUGUUUAAAUCUAAAUCUGGUAUGGGAGGCGUGCCCGUUGGUAUGAAGAAAUUUGCUGGGAGAAAAAGUCUAAGGCCUUGAAUCCGGCGGGUCUUAAUACUUACCUACAAAAUCAUCAGGCUGUACUUCCUGUAUC"
s1 =       "........(((((((((((((((((..((((.((.((((((((((...((.((((((((((....(((.((((((.......)))))).....)))....))))..(((((...(((..(((((((....(((......((((.................((((((...))))))........((((((((((((.(((((((((((...........((((......((((((......)))))).....))))....(((((((((((.(.((((((......))..)))).).)))).....))))))).........(((........))).))))))))))).)))..)))))))))....)))).......))))))).)))..))).))))).))))..)).))..))).))).)))).))..))))..((((((....)))))).....)))).)))))............))))))))((((((((.(((.(.((((.........((((..(((((.....((.((((((((((....))..)))))))))).))).))..)))))))).).))))))..)))))....."
s2 =       ".............((((((((((((.(((((.((.((((((((((...((.((((((((((...(((..((((((.......))))))....))).....))))..(((((...(((..(((((((....(((......((((.................((((((...))))))........((((((((((((.(((((((((((((.....))..((((......((((((......)))))).....))))....((((((((((((..((((........)))))))))...........))))))).........(((........))).))))))))))).)))..)))))))))....)))).......))))))).)))..))).))))).))))..)).))..))).))).)))).))..)))...((((((....))))))))...)))).))))))))............(((((((((((((.(((.((.(((((.((((..(((((((((.......(((...)))......))))))))))))).))))).))..............)))))..)))))))))))"

#100 nt with 2 inner sections
# sequence = "CGCAUCUCUUUAGGGUAUGAAAUGUUAUAUGCUACGGGAACAAUGCCGACCUUCGGAGACCUAAGGAAUACGUCUUUCGAGCGGAAGGAUUCCUCGUUCA"
# s1 =       ".(((((((..(((.((((((....)))))).))).))))....)))(((...)))..(((...((((....(((((((.....))))))))))).))).."
# s2 =       "....((.((.(((.((((((....)))))).))).))))......((((...)))).(((...(((((....((((((.....))))))))))).))).."



search_width = RNA.bp_distance(s1, s2)*2

# search_width = RNA.bp_distance(s1, s2)*4

# search_width = 4

# findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_f_minmax_debug"
# findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_orig"
# findpath_cmd = "/mnt/c/Users/Max/Google\ Drive/arbeit/cppfinder/cppfinder/findpath_orig_th"


# 60.0 works with findpath_f_minmax but not with minmax2...

# findpath_cmd = "/scratch/maxf/cppfinder/findpath_f"
# findpath_cmd = "/scratch/maxf/cppfinder/findpath_f_debug"
# findpath_cmd = "/scratch/maxf/cppfinder/findpath_f_all"

def launch_fp(sequence, s1, s2, search_width, findpath_cmd):
    start_findpath = time.time()
    cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width} '
    # cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width} -v'
    result1 = subprocess.check_output(cmd, shell=True, encoding="utf8")
    end_findpath = round(time.time()-start_findpath,4)
    return end_findpath, result1




start_c = time.time()
# with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
#     result = executor.map(launch_main, [sequence, sequence], [s1, s2], [s2, s1], [search_width, search_width])
# result = list(result)
result = launch_fp(sequence, s1, s2, search_width, "./main")
end_c = round(time.time()-start_c,4)

# print ("new:", result)

result = result[1].split("\n")
for row in result:
    print (row)


# result = result[0]
# end_findpath, result1 = launch_main(sequence, s1, s2, search_width)




# findpath_cmd = "./findpath_orig"
# start_findpath2 = time.time()
# cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width} '
# # cmd = f'printf "{sequence}\n{s1}\n{s2}" | {findpath_cmd} -m {search_width} -v'
# result2 = subprocess.check_output(cmd, shell=True, encoding="utf8")
# end_findpath2 = round(time.time()-start_findpath2,4)



start_findpath3 = time.time()
orig = pathfinder.pathfinder(sequence, s1, s2, search_width=search_width, verbose=False).max_en
# print ("orig pathfinder:", orig, round (time.time()-start,4), "s")
end_findpath3 = round(time.time()-start_findpath3,4)


# print (result)
print ("time elapsed th:", end_c)
# print ("time elapsed: no_th", end_findpath,  result1)
# print ("findpath_orig elapsed:", end_findpath2, result2)
print ("python bindings:", end_findpath3, orig)