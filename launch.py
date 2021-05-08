#!/usr/bin/env python3
# coding: utf-8

import time
import subprocess
import RNA

import findpath
import findpath_librna

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
# sequence = "UGAAGACCCAUUGAGUAACGACACCGCACGGCGCAUGGCGUCAGAGUAGCACUGCCUCGU"
# s1 =       "....(((((((.(.(...((....)).....).))))).)))((.((......))))..."
# s2 =       "....(((((((...((..((........))..)))))).))).(((..((...))))).."


# no sections 300
# sequence = "UCACGACACCCCUCAACUAUAACAACGGUCCGUACAUACUAGCCCUGCAAUGGAACGGGCAGGGCCAGCACAGGUGGGGCGCCCGCUUGGGGGAUCAAAUGUGUGAUCCAGAUACUUUAGACGCGUGCAGAACUUUUUAGAUCGAUCAGUGGGAACAGGCAUUGAUUAUGAAAUCAAUUAGGGGGUUUAGGACCGCACCACAAACUGCGGGAGGGCACGCUUUGGUUCCUGUGUUACGCUAAUCCUCUAGCCACGGAGGGCUUCUUCGUACAAUGAUUGGGUUACCAGGGUUCCAGUGUG"
# s1 = ".........(((((.(((........)))............(((((((..........))))))).....((((((((...))))))))..((((((......))))))......(((((((.(.(........(((((....((((((((........)))))))).))))).......).).))))))).(((((........))))))))))((((((..((.(((((.((.((.((((((........((((((.....)))))).....)))))))).)))))))..))))))))"
# s2 = ".(((...((((........((((..(((((.((((...((.(((((((..........))))))).))..((((((((...))))))))((((((.....(((((((.(((..((..((((.((((((.(((((((((((((.((((((((........)))))))).....)))...))))))))))....(((((........)))))....))))))))))))..))).)))))))..)))))).((((......)))).....))))...))))).))))...))))....))).."

# 3 ext loops
# sequence = "AACCCGCGAUAAGUUGGCGCUUGUCCCACUCCGUAAACCUGUGUCUCUCAGGCGGUUACCCGAUAGAAGGCAGUAGGAUGUAUCACCCCCCACCGACUCCACUAUACAACGAAACCGGCCCUUUGUGAAUAAACGCUCUCAGUUAGAUUGAGCGGCGCAACAGAACCAGAUCGAUUCCCGGGGCAGUAGUGUCCUACUGCCCAACCCGGUAGAGGGAUCGAGUGGUCUAGGCGAAUGAUACGAGAGGCUCACUGGGGACGAUGGGAGUGAUCUACUCGAUGUUGCGUCCGAUCACAUCCA"
# s1 =       "....(((.(.....).))).(((((((............((.(((((((.(((((((.....((((..((.(((.((.((..........)))).))))).)))).......))))).)).............(((((...........))))).(((...(((.(((..(((((((((((((((((((....))))))))....))).....)))))))).))))))..)))........))))))).))...)))))))...(((((...)))))(((((...........))))).."
# s2 =       "....(((.(.....).))).(((((((............((.(((((((.(((((((...((((((..((.(((.((.((..........)))).))))).))))....)).)))).))).((((....))))(((((...........))))).(((...((.((((..(((((((((((((((((((....)))))))))....)).....)))))))).))))))..)))........))))))).))...)))))))((((.((((((....((......))...)))))).))))"

# sequence = "CUAUCACCUGUACUUGCGGCCACCACCUCAUUUGCCCGUUGGCCGACUCAUCGCGUGGUAGGGGUUGUCGACUGUGGGAGUUUUCAGAUUACGUCCUAACCAACCGUAAAAGAACCUCCGGUAAGAGUAAAAGGAGGGUACAGGGUACAGUAAAGACGAUAGGACAGCCGUGUGUUAGGGCCGGAAGCCGCGGGCAAAUU"
# s1 = ".......((((....))))..........(((((((((.((((((.(((..((..((((((((((.(((((((.....))))....))).)).)))).))))..))........(((((.............))))).((((.(((...((....)).........))).))))...))).))...)))))))))))))."
# s2 = ".......((((....))))..........(((((((((((((((.......((((((((..(..((((((((((((((.((.........)).)))......(((.(.......(((((.............)))))....).)))))))).....))))))..).))))))))....)))))......))))))))))."

# sequence = "CACAUGGGAAGCAAUGUAGCUCUCUGAUUCAGCUCUGCUACCAAAUAAUUGGCCAGGUCCAGCGCGAUCUCUCGGUUGGGCAACUAUGUCAGUGGAGGGU"
# s1 = "...((.((((((......))).))).))...........(((.....((((((...(((((((.(((....))))))))))......))))))....)))"
# s2 = ".....(((.(((......)))))).......(((((.((((........((((..(((...((.(((((....))))).)).)))..)))))))))))))"


# sequence = "CACAACGAGACUACUCAUGCACAUAAGGAAGGUUAUCGCCAUAGACAUGGCGCGGCAGCGCAGAAUGUUUAAAUCUAAAUCUGGUAUGGGAGGCGUGCCCGUUGGUAUGAAGAAAUUUGCUGGGAGAAAAAGUCUAAGGCCUUGAAUCCGGCGGGUCUUAAUACUUACCUACAAAAUCAUCAGGCUGUACUUCCUGUAUC"
# s1 = "......(((....))).........((((((......(((.(((((((.((((....))))...)))))))...........((((((.....)))))).((.(((((.((((..(((((((((.....(((.(....).)))...))))))))))))).))))).))..............)))....))))))....."
# s2 = "......(((....)))((((.....((((.(((.((.(((.(((((((...(((....)))...)))))))...........((((((.....)))))).((.(((((.((((..(((((((((.....(((((...)).)))...))))))))))))).))))).))..............))).)))))))))))))."

# sequence = "UGAAGACCCAUUGAGUAAAA"
# s1       = "(((((((.....)))).)))"
# s2       = "(((((.........)).)))"

# search width m. 1 reference librna: -31.4 / single: -33.10 / merge fp: -32.50
# sequence = "UAAGGAAACUGAUGAGGGCAAAGUCUCUUCAUUGGCGCAAAACGGGGUAGAUGUCGGCCUGCAGGUGAUGAUGUAUCCGUCCGGCGAUCAACGCCUAAUUCGCAUCUCUUUAGGGUAUGAAAUGUUAUAUGCUACGGGAACAAUGCCGACCUUCGGAGACCUAAGGAAUACGUCUUUCGAGCGGAAGGAUUCCUCGUUCA"
# s1 = ".........(((((((((....((((((...(((((((((...((.((.((((((((.(((((........))))...).)))))).)).)).))...)).)).((((..(((.((((((....)))))).))).)))).....))))).....))))))...........(((((((.....)))))))))))))).))"
# s2 = "..(((((.((((((((((.......))))))))))........((.((.(((((((((.((((........))))...).))))).))).)).))......((((.(((((((.((((((....)))))).))).))))...))))(..(((((((((((...........))).)))))).))..)..)))))......"


# print (sequence)
# print (s1)
# print (s2)
print ("~~~~~~~~~~~")

search_width_multiplier = 4
# search_width_multiplier = 10

search_width = RNA.bp_distance(s1, s2)*search_width_multiplier
print ("search width:", search_width, "bp_dist:", RNA.bp_distance(s1, s2), "seq_length:", len(sequence))



fc = RNA.fold_compound(sequence)
s1_eval = fc.eval_structure(s1)


bp_dist = RNA.bp_distance(s1, s2)

# print ("bp dist:", bp_dist, len(sequence))

# mp = False
mp = True

findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp)

# en_limit = -15191

start_findpath = time.time()
# result = findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp, en_limit)
result = findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp)
result = result/100.0
runtime = time.time()-start_findpath
print ("~~~~~~~~~~~")
print (f'single section')
print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')


# start_findpath = time.time()
# result = findpath.init_multi_findpath(sequence, s1, [s2], search_width_multiplier, mp)
# result = result/100.0
# runtime = time.time()-start_findpath
# print ("~~~~~~~~~~~")
# print (f'multi section')
# print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')




start_findpath = time.time()
fp = findpath.findpath_class(sequence, mp)
result = fp.init(s1, s2, search_width_multiplier)
result = result/100.0
runtime = time.time()-start_findpath
print ("~~~~~~~~~~~")
print (f'merge findpath')
print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')


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



# runtimes = []
# for i in range(10):
#     start_findpath = time.time()
#     findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp)
#     runtimes.append(time.time()-start_findpath)
# print (f'min runtime: {min(runtimes):2.4f} s / total runtime: {sum(runtimes):2.4f}')



# runtimes = []
# fp = findpath.findpath_class(sequence, mp)

# for i in range(10):
#     start_findpath = time.time()
    
#     result = fp.init(s1, s2, search_width_multiplier)
#     print (result)

#     runtimes.append(time.time()-start_findpath)
# print (f'min runtime: {min(runtimes):2.4f} s / total runtime: {sum(runtimes):2.4f}')


# runtimes = []
# for i in range(10):
#     start_findpath = time.time()    
#     result = findpath.init_merge_findpath(sequence, s1, s2, search_width_multiplier, True)
#     runtimes.append(time.time()-start_findpath)
# print (f'min runtime: {min(runtimes):2.4f} s / total runtime: {sum(runtimes):2.4f}')