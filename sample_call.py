#!/usr/bin/env python3
# coding: utf-8

import time
import RNA

import findpath
import findpath_librna

from pretty_print_path import print_moves

# random 600 nt example
sequence = "AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAAUUUUAUUAAAAGUCCAAGUUGGACUGACAAAACGCGUGCGGUGUCCUAGGGAUUGGUGGCAUAACCAGCGGUUUAAAAGCUGUGUAUAUCCGCAGCAAAUCACCGGAAAGCGGCGUUAUUAGCACCACAAAUUGAUGGUUGGUACGAGUACAAUUGCGCCGCAUAAAACCAGAGAUUCUACCCUCAAUCGGUUCUUAAGACGUACUGCGCGUUUCACCAGACCACAAUGCAGGGCGGCACCGUUAGGCAACACAACGAGACUACUCAUGCACAUAAGGAAGGUUAUCGCCAUAGACAUGGCGCGGCAGCGCAGAAUGUUUAAAUCUAAAUCUGGUAUGGGAGGCGUGCCCGUUGGUAUGAAGAAAUUUGCUGGGAGAAAAAGUCUAAGGCCUUGAAUCCGGCGGGUCUUAAUACUUACCUACAAAAUCAUCAGGCUGUACUUCCUGUAUC"
s1 =       "........(((((((((((((((((..((((.((.((((((((((...((.((((((((((....(((.((((((.......)))))).....)))....))))..(((((...(((..(((((((....(((......((((.................((((((...))))))........((((((((((((.(((((((((((...........((((......((((((......)))))).....))))....(((((((((((.(.((((((......))..)))).).)))).....))))))).........(((........))).))))))))))).)))..)))))))))....)))).......))))))).)))..))).))))).))))..)).))..))).))).)))).))..))))..((((((....)))))).....)))).)))))............))))))))((((((((.(((.(.((((.........((((..(((((.....((.((((((((((....))..)))))))))).))).))..)))))))).).))))))..)))))....."
s2 =       ".............((((((((((((.(((((.((.((((((((((...((.((((((((((...(((..((((((.......))))))....))).....))))..(((((...(((..(((((((....(((......((((.................((((((...))))))........((((((((((((.(((((((((((((.....))..((((......((((((......)))))).....))))....((((((((((((..((((........)))))))))...........))))))).........(((........))).))))))))))).)))..)))))))))....)))).......))))))).)))..))).))))).))))..)).))..))).))).)))).))..)))...((((((....))))))))...)))).))))))))............(((((((((((((.(((.((.(((((.((((..(((((((((.......(((...)))......))))))))))))).))))).))..............)))))..)))))))))))"

# 100 nt with 2 inner sections
# sequence = "CGCAUCUCUUUAGGGUAUGAAAUGUUAUAUGCUACGGGAACAAUGCCGACCUUCGGAGACCUAAGGAAUACGUCUUUCGAGCGGAAGGAUUCCUCGUUCA"
# s1 = ".(((((((..(((.((((((....)))))).))).))))....)))(((...)))..(((...((((....(((((((.....))))))))))).))).."
# s2 = "....((.((.(((.((((((....)))))).))).))))......((((...)))).(((...(((((....((((((.....))))))))))).))).."

# sequence = "UCUACUAUUCCGGCUUGACAUAAAUAUCGAGUGCUCGACCGCUAUUAUGGUACUUUCCAGCGUUUUGAUUGGUGGAUAAUAUCCCCCAAAAACGCGAGUC"
# s1 = "............(((((..........)))))((((..((........)).........((((((...((((.((((...)))).))))))))))))))."
# s2 = "..((((((..(((.(((((((.........))).))))))).....)))))).....(..(((((...((((.((((...)))).))))))))).)...."

# sequence = "UGAAGACCCAUUGAGUAAAA"
# s1       = "(((((((.....)))).)))"
# s2       = "(((((.........)).)))"
# s2       = "(((((((....)).)).)))"

# settings
search_width_multiplier = 4
mp = False
# mp = True
Verbose = True

# for debug
search_width = RNA.bp_distance(s1, s2)*search_width_multiplier
print("search width:", search_width, "bp_dist:",
      RNA.bp_distance(s1, s2), "seq_length:", len(sequence))
fc = RNA.fold_compound(sequence)
s1_eval = fc.eval_structure(s1)
bp_dist = RNA.bp_distance(s1, s2)



# single section findpath
print("~~~~~~~~~~~")
print(f'single section')
start_findpath = time.time()

# call with search width multiplier:
fp = findpath.findpath_single(sequence, s1, s2, search_width_multiplier=search_width_multiplier, mp=mp)

# call with fixed search width integer:
# fp = findpath.findpath_single(sequence, s1, s2, search_width=50, mp=mp)

# call with added model detail dictionary:
# fp = findpath.findpath_single(sequence, s1, s2, search_width_multiplier=search_width_multiplier, 
#     mp=mp, model_details={"temperature": -14.5,})

runtime = time.time()-start_findpath
result = fp.get_en()/100.0
path = fp.get_path()
print_moves(sequence, s1, s2, path, convert_to_float=True)
print(f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')


# ---------------

# model dict schema
# "noLP":        int,
# "logML":       int,
# "temperature": float,
# "dangles":     int,
# "special_hp":  int,
# "noGU":        int,
# "noGUclosure": int,

# example: set model details
model_details = {
    "noLP": 1,
    "logML": 0,
    "temperature": 37.0,
    "dangles": 2,
    "special_hp": 1,
    "noGU": 0,
    "noGUclosure": 0,
}

# set only T
model_details = {
    "temperature": 37.0,
}

# # first merge findpath call
# print("~~~~~~~~~~~")
# print(f'first merge findpath')
# start_findpath = time.time()
# fp = findpath.findpath_class(sequence, mp=mp, model_details=model_details)
# fp.init(s1, s2, search_width_multiplier)

# runtime = time.time()-start_findpath
# result = fp.get_en()/100.0
# path = fp.get_path()
# print_moves(sequence, s1, s2, path, convert_to_float=True)
# # print(fp.get_sections())
# print(f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')



# # second merge findpath call
# print("~~~~~~~~~~~")
# print(f'second merge findpath')
# start_findpath = time.time()
# # fp = findpath.findpath_class(sequence, mp=mp, model_details=model_details)
# fp.init(s1, s2, search_width_multiplier)

# runtime = time.time()-start_findpath
# result = fp.get_en()/100.0
# path = fp.get_path()
# print_moves(sequence, s1, s2, path, convert_to_float=True)
# # print(fp.get_sections())
# print(f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')
