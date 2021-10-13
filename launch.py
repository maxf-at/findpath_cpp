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

#100 nt with 2 inner sections
# sequence = "CGCAUCUCUUUAGGGUAUGAAAUGUUAUAUGCUACGGGAACAAUGCCGACCUUCGGAGACCUAAGGAAUACGUCUUUCGAGCGGAAGGAUUCCUCGUUCA"
# s1 =       ".(((((((..(((.((((((....)))))).))).))))....)))(((...)))..(((...((((....(((((((.....))))))))))).))).."
# s2 =       "....((.((.(((.((((((....)))))).))).))))......((((...)))).(((...(((((....((((((.....))))))))))).))).."

# sequence = "GUUGGGGUAGGGGCCCGACUAAUAGCAGUUCUUUGACAAUUCUUCUGCGUUAUUCAUUUUGAUAACAAUUAUUAUAUAAGUGCUGGAAAGCCAACCAUAUCGGCCUUAAUCCCCAACAGAAAACCUUCACGAGUGGGGUCGGUUUCGCAAUUCCUUGCGUGAGAGCCGAGGUUGGUUACAUGAGAAUAACCACGCCGAAC"
# s1 = "(((((((..((((((.((.((((.((((................)))))))).))((((..((((......))))..))))..(((....)))........))))))...))))))).......((.((((...((((..............)))).))))))...((..((.(((((........)))))))..))..."
# s2 = "((((((((((.((((.........((((................))))((((((......)))))).............(((.(((....)))..)))...)))))))..)))))))......(((........))).((((((((((((....)))))...)))))))(((((((((........)))))).)))...."

# sequence = "GGUAUGCUGCCUAUAAAGUAAGUGCACUGAUCUACAAGCAUAAGCGGGACGCUUCGUGAGUACAUUCAGA"
# s1 = "....((((((.......)))...)))((((..(((...((((((((...))))).))).)))...))))."
# s2 = ".(((((((........)))..)))).((((..(((..((..(((((...))))).))..)))...))))."

# this example requires suboptimal paths for merging
# S: -2.50, barrier: 1.90, runtime: 0.0023 s
sequence = "UGAAGACCCAUUGAGUAACGACACCGCACGGCGCAUGGCGUCAGAGUAGCACUGCCUCGU"
s1 =       "....(((((((.(.(...((....)).....).))))).)))((.((......))))..."
s2 =       "....(((((((...((..((........))..)))))).))).(((..((...))))).."

# no sections 300
# sequence = "UCACGACACCCCUCAACUAUAACAACGGUCCGUACAUACUAGCCCUGCAAUGGAACGGGCAGGGCCAGCACAGGUGGGGCGCCCGCUUGGGGGAUCAAAUGUGUGAUCCAGAUACUUUAGACGCGUGCAGAACUUUUUAGAUCGAUCAGUGGGAACAGGCAUUGAUUAUGAAAUCAAUUAGGGGGUUUAGGACCGCACCACAAACUGCGGGAGGGCACGCUUUGGUUCCUGUGUUACGCUAAUCCUCUAGCCACGGAGGGCUUCUUCGUACAAUGAUUGGGUUACCAGGGUUCCAGUGUG"
# s1 = ".........(((((.(((........)))............(((((((..........))))))).....((((((((...))))))))..((((((......))))))......(((((((.(.(........(((((....((((((((........)))))))).))))).......).).))))))).(((((........))))))))))((((((..((.(((((.((.((.((((((........((((((.....)))))).....)))))))).)))))))..))))))))"
# s2 = ".(((...((((........((((..(((((.((((...((.(((((((..........))))))).))..((((((((...))))))))((((((.....(((((((.(((..((..((((.((((((.(((((((((((((.((((((((........)))))))).....)))...))))))))))....(((((........)))))....))))))))))))..))).)))))))..)))))).((((......)))).....))))...))))).))))...))))....))).."

# 3 ext loops
# sequence = "AACCCGCGAUAAGUUGGCGCUUGUCCCACUCCGUAAACCUGUGUCUCUCAGGCGGUUACCCGAUAGAAGGCAGUAGGAUGUAUCACCCCCCACCGACUCCACUAUACAACGAAACCGGCCCUUUGUGAAUAAACGCUCUCAGUUAGAUUGAGCGGCGCAACAGAACCAGAUCGAUUCCCGGGGCAGUAGUGUCCUACUGCCCAACCCGGUAGAGGGAUCGAGUGGUCUAGGCGAAUGAUACGAGAGGCUCACUGGGGACGAUGGGAGUGAUCUACUCGAUGUUGCGUCCGAUCACAUCCA"
# s1 =       "....(((.(.....).))).(((((((............((.(((((((.(((((((.....((((..((.(((.((.((..........)))).))))).)))).......))))).)).............(((((...........))))).(((...(((.(((..(((((((((((((((((((....))))))))....))).....)))))))).))))))..)))........))))))).))...)))))))...(((((...)))))(((((...........))))).."
# s2 =       "....(((.(.....).))).(((((((............((.(((((((.(((((((...((((((..((.(((.((.((..........)))).))))).))))....)).)))).))).((((....))))(((((...........))))).(((...((.((((..(((((((((((((((((((....)))))))))....)).....)))))))).))))))..)))........))))))).))...)))))))((((.((((((....((......))...)))))).))))"

sequence = "UGAAGACCCAUUGAGUAAAA"
s1       = "(((((((.....)))).)))"
s2       = "(((((.........)).)))"

# sequence = "AAGAAGACCUCAAUCGAAUCACGGGCAAGUCCGACGAGGAACGCCUAGGCGAGGUGAUCGGCCCGAUCUUAAUGUAGGAUCCCCGGAGUCGCAUGACGACAGCUUAAUGUUCGUCCAGGGGGCAUACCCUUGGUGACUGUAAGCCGUGCCUGGUCCUUUUCUCGAAUGAGUCCACAGAUUAGCAAAUUUAAAAAGUGCGG"
# s1 = ".....((.((((.((((....(((......)))..(((((....(((((((.(((.(((((((((((((((...)))))))...)).)))).)).))....(((((..((.((.(((((((.....)))))))))))..)))))).))))))))))))...)))).)))))).........(((...........))).."
# s2 = ".....((.((((.((((....(((......)))..((((...(((.(((((.(((...(((...(((((((...))))))).))).((((((..(((((((......))).))))((((((.....))))))))))))....))).))))))))))))...)))).)))))).........(((..((....)).))).."

# sequence = "UAAAAUGAUCACGGUUUCAGCUUUGGACAGGGCGUUCCACUAAACUCCUGGUGACAUAGAUAUAUUGGAUUGCAACUACUCGUCGGUCCGGUUGGCGUUCAACCCGCGAUAAGUUGGCGCUUGUCCCACUCCGUAAACCUGUGUCUCUCAGGCGGUUACCCGAUAGAAGGCAGUAGGAUGUAUCACCCCCCACCGACUCC"
# s1 = ".....(((........))).....(((..(((.(((......))))))(((((.....((((((((..(((((..((..(((((((...(((((.....)))))...(((((((....))))))).......((.((((...((((...))))))))))))))).)))))))))..))))))))......)))))..)))"
# s2 = ".....(((........)))...((((...(((.((((((((.(((((((((.(((((((.......(((((((........).))))))(((((.....)))))...(((((((....))))))).............))))))).).)))).))))..((.......)).))).)))......)))))...))))...."


# 4
search_width_multiplier = 4
# mp = False
mp = True
Verbose = True


search_width = RNA.bp_distance(s1, s2)*search_width_multiplier
print ("search width:", search_width, "bp_dist:", RNA.bp_distance(s1, s2), "seq_length:", len(sequence))
fc = RNA.fold_compound(sequence)
s1_eval = fc.eval_structure(s1)
bp_dist = RNA.bp_distance(s1, s2)

# print ("bp dist:", bp_dist, len(sequence))


# warmup
# findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp)

# en_limit = -15191

start_findpath = time.time()
result = findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp)
result = result/100.0
runtime = time.time()-start_findpath
print ("~~~~~~~~~~~")
print (f'single section')
print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')

start_findpath = time.time()
fp = findpath.findpath_single(sequence, s1, s2, search_width_multiplier=search_width_multiplier, mp=mp)
# fp = findpath.findpath_single(sequence, s1, s2, search_width=50, mp=mp)
runtime = time.time()-start_findpath
result = fp.get_en()/100.0
# path = fp.get_path()
# print_moves(sequence, s1, s2, path, convert_to_float=True)
print(f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')


# start_findpath = time.time()
# fp = findpath.findpath_class(sequence, mp)
# fp.init(s1, s2, search_width_multiplier)
# result = fp.get_en()/100.0
# runtime = time.time()-start_findpath
# print ("~~~~~~~~~~~")
# print (f'merge findpath')
# print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')


# start_findpath = time.time()
# fp.init(s1, s2, search_width_multiplier)
# result = fp.get_en()/100.0
# runtime = time.time()-start_findpath
# print ("~~~~~~~~~~~")
# print (f'cached merge findpath call')
# print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')


# path = fp.return_path()
# print (path) 
# print_moves(sequence, s1, s2, path, convert_to_float=True)
# fp.return_sections()

# start_findpath = time.time()
# result = findpath.init_mfe_findpath(sequence, s1, s2, search_width_multiplier, True)
# result = result/100.0
# runtime = time.time()-start_findpath
# print ("~~~~~~~~~~~")
# print (f'merge findpath (+MFE)')
# print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')

# start_findpath = time.time()
# result = findpath.init_merge_ext_findpath(sequence, s1, s2, search_width_multiplier, True)
# result = result/100.0
# runtime = time.time()-start_findpath
# print ("~~~~~~~~~~~")
# print (f'merge findpath (+ext loops)')
# print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')

# start_findpath = time.time()
# result = findpath_librna.pathfinder(sequence, s1, s2, search_width=search_width_multiplier*bp_dist)
# runtime = time.time()-start_findpath
# print ("~~~~~~~~~~~")
# print (f'librna findpath (orig SWIG bindings)')
# print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')

# start_findpath = time.time()
# result = findpath.init_vrna_findpath(sequence, s1, s2, search_width_multiplier, True)
# result = result/100.0
# runtime = time.time()-start_findpath
# print ("~~~~~~~~~~~")
# print (f'librna findpath (PYBIND11)')
# print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')



# caching
# runtimes = []
# fp = findpath.findpath_class(sequence, mp)
# for i in range(10):
#     start_findpath = time.time()    
#     result = fp.init(s1, s2, search_width_multiplier)
#     print (result)
#     runtimes.append(time.time()-start_findpath)
# print (f'min runtime: {min(runtimes):2.4f} s / total runtime: {sum(runtimes):2.4f}')

# find the average runtime
# runtimes = []
# for i in range(10):
#     start_findpath = time.time()    
#     result = findpath.init_merge_findpath(sequence, s1, s2, search_width_multiplier, True)
#     runtimes.append(time.time()-start_findpath)
# print (f'min runtime: {min(runtimes):2.4f} s / total runtime: {sum(runtimes):2.4f}')
# runtimes = []
# for i in range(10):
#     start_findpath = time.time()
#     findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp)
#     runtimes.append(time.time()-start_findpath)
# print (f'min runtime: {min(runtimes):2.4f} s / total runtime: {sum(runtimes):2.4f}')
