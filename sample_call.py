#!/usr/bin/env python3
# coding: utf-8

import time
import RNA

import findpath
import findpath_librna

from pretty_print_path import print_moves

# random 600 nt example
# sequence = "AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAAUUUUAUUAAAAGUCCAAGUUGGACUGACAAAACGCGUGCGGUGUCCUAGGGAUUGGUGGCAUAACCAGCGGUUUAAAAGCUGUGUAUAUCCGCAGCAAAUCACCGGAAAGCGGCGUUAUUAGCACCACAAAUUGAUGGUUGGUACGAGUACAAUUGCGCCGCAUAAAACCAGAGAUUCUACCCUCAAUCGGUUCUUAAGACGUACUGCGCGUUUCACCAGACCACAAUGCAGGGCGGCACCGUUAGGCAACACAACGAGACUACUCAUGCACAUAAGGAAGGUUAUCGCCAUAGACAUGGCGCGGCAGCGCAGAAUGUUUAAAUCUAAAUCUGGUAUGGGAGGCGUGCCCGUUGGUAUGAAGAAAUUUGCUGGGAGAAAAAGUCUAAGGCCUUGAAUCCGGCGGGUCUUAAUACUUACCUACAAAAUCAUCAGGCUGUACUUCCUGUAUC"
# s1 =       "........(((((((((((((((((..((((.((.((((((((((...((.((((((((((....(((.((((((.......)))))).....)))....))))..(((((...(((..(((((((....(((......((((.................((((((...))))))........((((((((((((.(((((((((((...........((((......((((((......)))))).....))))....(((((((((((.(.((((((......))..)))).).)))).....))))))).........(((........))).))))))))))).)))..)))))))))....)))).......))))))).)))..))).))))).))))..)).))..))).))).)))).))..))))..((((((....)))))).....)))).)))))............))))))))((((((((.(((.(.((((.........((((..(((((.....((.((((((((((....))..)))))))))).))).))..)))))))).).))))))..)))))....."
# s2 =       ".............((((((((((((.(((((.((.((((((((((...((.((((((((((...(((..((((((.......))))))....))).....))))..(((((...(((..(((((((....(((......((((.................((((((...))))))........((((((((((((.(((((((((((((.....))..((((......((((((......)))))).....))))....((((((((((((..((((........)))))))))...........))))))).........(((........))).))))))))))).)))..)))))))))....)))).......))))))).)))..))).))))).))))..)).))..))).))).)))).))..)))...((((((....))))))))...)))).))))))))............(((((((((((((.(((.((.(((((.((((..(((((((((.......(((...)))......))))))))))))).))))).))..............)))))..)))))))))))"

#100 nt with 2 inner sections
sequence = "CGCAUCUCUUUAGGGUAUGAAAUGUUAUAUGCUACGGGAACAAUGCCGACCUUCGGAGACCUAAGGAAUACGUCUUUCGAGCGGAAGGAUUCCUCGUUCA"
s1 =       ".(((((((..(((.((((((....)))))).))).))))....)))(((...)))..(((...((((....(((((((.....))))))))))).))).."
s2 =       "....((.((.(((.((((((....)))))).))).))))......((((...)))).(((...(((((....((((((.....))))))))))).))).."


search_width_multiplier = 4

# mp = False
mp = True
Verbose = True


search_width = RNA.bp_distance(s1, s2)*search_width_multiplier
print ("search width:", search_width, "bp_dist:", RNA.bp_distance(s1, s2), "seq_length:", len(sequence))
fc = RNA.fold_compound(sequence)
s1_eval = fc.eval_structure(s1)
bp_dist = RNA.bp_distance(s1, s2)


# single section findpath call
start_findpath = time.time()
result = findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp)
result = result/100.0
runtime = time.time()-start_findpath
print ("~~~~~~~~~~~")
print (f'single section')
print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')


# merge findpath call
start_findpath = time.time()
fp = findpath.findpath_class(sequence, mp)
result = fp.init(s1, s2, search_width_multiplier)
result = result/100.0
runtime = time.time()-start_findpath
print ("~~~~~~~~~~~")
print (f'merge findpath')
print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')

# return paths from merge call
path = fp.return_path()
print (path) 
print_moves(sequence, s1, s2, path, convert_to_float=True)
fp.return_sections()
