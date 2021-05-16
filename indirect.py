#!/usr/bin/env python3
# coding: utf-8

import time
import RNA

import findpath
import findpath_librna

sequence = "GGGCGCGGUUCGCCCUCCGCUAAAUGCGGAAGAUAAAUUGUGUCU"
s1 = "(((((.....)))))(((((.....)))))(((((.....)))))"
s2 = "((((((((((.....(((((.....))))).....))))))))))"
add_moves =[6, 12, 7, 11, 5, 13]

# sequence = 'AAACGUAAAAGUGAGUACCAUCCGAUAGUCUACGUGGUGCUUUAAUCAAACACCUAACGGAUCAUCGUACGCCGUCAGCAUAAAAUGAUCACGGUUUCAGCUUUGGACAGGGCGUUCCAC'
# s1       = '............(((((((((............)))))))))......(((.(((..((((((((..((.((.....)).))..)))))).))((.(((....))))).))).)))....'
# s2       = '...((((...(((((((((((..((...))...))))))))).......))...)).)).........((((((((..((.....(((........)))....)))))..))))).....'
# add_moves = [58, 67, 59, 66, 57, 68]

# sequence = 'UAAACUCCUGGUGACAUAGAUAUAUUGGAUUGCAACUACUCGUCGGUCCGGUUGGCGUUCAACCCGCGAUAAGUUGGCGCUUGUCCCACUCCGUAAACCUGUGUCUCUCAGGCGGUUACC'
# s1       = '(((.(((((((.(((((((...(((.(((.((..((....((((((((((((((.....))).))).)))....)))))...))..)).))))))...)))))))..))))).)))))..'
# s2       = '..(((((((((.(((((((...(((.(((((((........).))))))(((((.....)))))...(((((((....))))))).......)))...))))))).).)))).))))...'
# add_moves = [-25, -93, -24, -94, -23, -95, 26, 93]

# sequence = 'AGCUUAAUGAUCCACGGGUUCGCUGCACUGCAGCUUUUCACGGACUCUGUGCCGUUUGGAGUGACCAGGUCUCUCCCGCUUCGCGCGAGGUGGGUCCACGCAAGCCGGCACGCGAAAGGG'
# s1       = '...............(((((((((((...))))).......))))))(((((((((((..((((((...((((...(((...))).))))..))).))).))))).))))))........'
# s2       = '........((.((...)).))(((((...)))))(((((........(((((((((((..(((....((.(((.(((((.....))).)).)))))))).)))).))))))).)))))..'
# add_moves = [22, 41, 8, 97, 12, 119, 69, 96, -63, -97, 13, 118]
# add_moves = [68, 97, 12, 119, -61, -99, 69, 96, -63, -97, -62, -98, 13, 118]
# # add_moves = [68, 97, 69, 96, -63, -97, -62, -98]

# sequence = 'AGCGGAGGGGUCUAGUGGUCACAUCCUACGAGUUCGGGGCACCGAACAAAUCGAAGUCUACCGAGAUAGGGCGUUCGUUAGGCGGACGCUAGCCCAUCAUUUCUCAUGGUUCAUUGAUAC'
# s1       = '....((.(((.((((((..(....(((((((((((((....)))))).............((.......))...)))).))).)..))))))))).))...((....))...........'
# s2       = '...((.(((((((..((((.((.((......((((((....)))))).....)).))..))))))))..(((((((((...)))))))))..))).))...((..(((...))).))...'
# add_moves = [20, 85, 18, 87, 105, 116, 61, 70, 106, 115, 104, 117, 19, 86, 62, 69]
# add_moves = [20, 85, 105, 116, 61, 70, 106, 115, 104, 117, 19, 86, 28, 79, 62, 69]



# print (sequence)
# print (s1)
# print (s2)
print ("~~~~~~~~~~~")

search_width_multiplier = 4
# search_width_multiplier = 100

search_width = RNA.bp_distance(s1, s2)*search_width_multiplier
print ("search width:", search_width, "bp_dist:", RNA.bp_distance(s1, s2), "seq_length:", len(sequence))


fc = RNA.fold_compound(sequence)
s1_eval = fc.eval_structure(s1)
bp_dist = RNA.bp_distance(s1, s2)


mp = False
mp = True
Verbose = True

# warmup
# findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp)


# direct
start_findpath = time.time()
result = findpath.init_single_findpath_i(sequence, s1, s2, search_width_multiplier, [], mp, Verbose=Verbose)
result = result/100.0
runtime = time.time()-start_findpath
print ("~~~~~~~~~~~")
print (f'single section')
print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')

# indirect
start_findpath = time.time()
result = findpath.init_single_findpath_i(sequence, s1, s2, search_width_multiplier, add_moves, mp, Verbose=Verbose)
result = result/100.0
runtime = time.time()-start_findpath
print ("~~~~~~~~~~~")
print (f'single section')
print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')
