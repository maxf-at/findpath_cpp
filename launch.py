#!/usr/bin/env python3
# coding: utf-8

import time
import subprocess
import RNA

import findpath
import findpath_librna

# sequence = "AUGAACAGAUGGUACCUCGCGGCGGGGCCUGACCCAUCUACAGUUUGUGCAG"
# s1 =       "((((((((((((..(((((...)))))......))))))...))))))...."
# s2 =       ".((..((((..(((((....))..(((.....)))...)))..))))..))."

# sequence = "AACGCCGCGCCGGCCAGGCUGCCAGUUGUCUUAUUGAGUUGGCCAUCUACGUAAUACUUUUGUAGCCGACGACUAUAACG"
# s1 =       "(((...(((((.....)).)))..)))...((((...((((((....((((.((....))))))))))))....)))).."
# s2 =       "...((.(((((.....))).))..)).(((.......((((((....((((.........)))))))))))))......."


# sequence = "UGGGUCCACGCAAGCCGGCACGCGAAAGGGCGGCUAUGUACAGCUCUUAAAACCACCAGAGGUUUAGUGAUCACUGAGGCUUGUUCGCAAAUCACUGCAAUUAGAUAUGACUCACGAUAUGGGGCACGGUGCAUACAUACAAGAAUUACUUGGCACCCCCGCGGAAAGUAGUCUCCUUGUACGCAUAUCGCUACAAAAUCGGGGAACUUCCAUACAUCCUGUUAUGUGAAAUGCGACGUAGAAGCCACACACUCGUCAUUGACAGAUCUACGAAUUAUGGGUCUACUAAGAAAUUUCACAUUCCGUAGAUUCUUGUACAAGGAUAGUGACUCCGCAAUUAGGACUGAACAAUAACUUGGGCGUAAAGUAAGCGCAGGACUAGAGGAAGGCUGGUGCGCGA"
# s1 =       ".(.(((((.(..(((((.(.(......)).))))).(((.(((.((((((((((......))))..(((((((((((.((.((.........))..))..))))...))).))))......(((...(((((.......((((.....))))))))))))(((((...((.(.(((((((((((((.((((((((..((..(((....)))....))..))).).)))).))))..((((((..........(((((...))).))))))))......((((((((...(((........))).))))))))..))))))))).).))..)))))..))))))))).)))....).))))).)........(((.(.(((((.......))))).))))."
# s2 =       "((.(((((.(..((((.((.(.......))))))).(((.(((.((((((((((......))))).((((((((((..((................))...)))...))).))))......(((...(((((....................))))))))(((((...(..(.(((((((((((((.((((.(((......(((....)))........)))...)))).))))((((.((..........)))))).......((((((((....((((((((....)))....)).)))..))))))))...))))))))).)...).)))))...)))))))).)))....).))))).))......((((...(((((.......))))))))).."

# this one is better with findpath_cpp
# sequence = "GCCGCUAGAUCGAGAUGUUUAGGAUAUGCAAACAGGUCUGGGGGGUCAUGUCUUUAAAGGUCAUCAUGAGGCUCCGCCAAGACUAGAGUCCGGCCGAGCCUACCACGGCUCAUCAUGGCGCACCCUCUUAUUUUGUCAACUUGCUCUUAUCCGAACGACAACUAUCACUUGGUUACUCUGCAAUGAUUACGUACCUAGUUGGGGACGACAAUUGCUCAAUCGUUCGCGCCGUACGUGUGGCGAAGGCCUUAUGAUGAUCAAAGAGACUGGCCAGCGCCCGGAAGUACGCCGAACGGAGCG"
# s1 =       ".((.(((((((....(((((.(......)))))))))))))))(((((.((((((...(((((((((((((((.(((((........((.((((.(((((......)))))((...((((((...........))......))))...)).(((((((.(((((.....)))))..((......))...(((.(((....))).)))............))))))).)))).))...)))))..)))))))))))))))..)))))))))))..((((((.............))).)))"
# s2 =       ".((.(((((((....(((.........)))....)))))))))(((((.((((((...(((((((((((((((.((((((((..((((....((((((((......))))).....))).....))))..))).((((..((((......((((..((......))..)))).......))))))))(((((((...(.((.((((((...........)))))).))).))))))))))))..)))))))))))))))..)))))))))))..((((((.............))).)))"

# 300 nt, 190 bp_dist, 2 sections
# good example sw_multiplier scaling... 
# sequence = "CAAUUGCGCCGCAUAAAACCAGAGAUUCUACCCUCAAUCGGUUCUUAAGACGUACUGCGCGUUUCACCAGACCACAAUGCAGGGCGGCACCGUUAGGCAACACAACGAGACUACUCAUGCACAUAAGGAAGGUUAUCGCCAUAGACAUGGCGCGGCAGCGCAGAAUGUUUAAAUCUAAAUCUGGUAUGGGAGGCGUGCCCGUUGGUAUGAAGAAAUUUGCUGGGAGAAAAAGUCUAAGGCCUUGAAUCCGGCGGGUCUUAAUACUUACCUACAAAAUCAUCAGGCUGUACUUCCUGUAUC"
# s1 =       "......(((((.......((..(((((((.(((.(((....(((((....(((((((((((((((.(((.((((...(((..((((....))))..))).......(((....)))...........((......))....(((((((.((((....))))...)))))))..........)))).))))))))))).....))))))))))))..)))..))))))....))))..)).........)))))(((..........)))...........((((........))))...."
# s2 =       "....((((((...........(((........)))....(((..(((((((.....(((((((((.(((.((((........(((((.(((...........((..(((....))).))...........)))..))))).(((((((.((((....))))...)))))))..........)))).))))))))))))(((((((..(.(((......(((.(((.......)))..)))))).)..)))))))))))))).....))).............))).)))..........."


# 600 nt, 2 sections
# sequence = "UUGUGUGAUGUGCAUCGGUAGCUCGGUAAAUACGUUCGCCGGGCUUGAAGCUGUAUCCGUUGGAAUCCAAGUUGGGAGUUAUAACGUAAAACGGCAUUCCCGAUAGUCUCUCAGGGGACGAUCGAUGUUCAGGCGGCUGCGAUUUUGAACGCCGUGAAGCGUGAUUAAACCGGUCGACUUUUGCUCCUUGAAGUUUUCAUCUAUCACCUGUACUUGCGGCCACCACCUCAUUUGCCCGUUGGCCGACUCAUCGCGUGGUAGGGGUUGUCGACUGUGGGAGUUUUCAGAUUACGUCCUAACCAACCGUAAAAGAACCUCCGGUAAGAGUAAAAGGAGGGUACAGGGUACAGUAAAGACGAUAGGACAGCCGUGUGUUAGGGCCGGAAGCCGCGGGCAAAUUGACGGGACAGCACCUAGACAUACUCAAAUUAGGCCUAGCCACUAUUUGGCCGCACGCUUGUAUCAGAUUGGUCAAAGUACACUCCAUUUCUAAAGAUUUACCCCUCUAGCCCCUCCAGUGAGUGGCGCUGUCCGUUCAUAAAACUCCGUACACUCAUAUAGUCGAUGGCGAACGCGAUGUGGUGUAACCGAUUUACGUUA"
# s1 =       ".....(((((((.((((((((((((((..........)))))))).((((((((...(((((...(((......)))....)))))....))))).))).((((.(((((....))))).))))....(((((((....((....))..)))).)))((.((((.((.....(..((((((........))))))..)...)))))))).......((((((.(.............).))))))(((((((((((.......(((((((((((((.((((.......(((((...........))))).....(((((.............)))))((((.((....(((......((..(((((((.((((((((((((.(....((((((.....))).)))..).)).)))).)))))).(((((..(((...)))...)))))(((((.((((..........(((............)))...(((.((........)).))).......)))).))))))))))))..)).....))))))))))))).)))))))))))))..)))))))).)))...)))))).)))))))"
# s2 =       ".....(((((((.((((((((((((((..........))))))))....((((((((((((((..(((......)))......((((...(((((.....((((.(((((....))))).))))..(((((((((....))..))))))))))))...))))........((((((((....(((......))).....((((((((.........((((((.(.............).)))))).........))))))))....))))))))(((((((.........)).))))).)))))..........(((((.............))))).....)))))))))..........(((((((.(((((((..((((((((((.((((.(((((.....((......(((((...........))))).....))...))))).))))..))))........))))))..)))))))...........(((........))).(((.(((....))).))))))))))...............(((((.(((...((((.((.....)))))))))))))))))))).)))))))"

# random 600 nt example
# sequence = "AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAAUUUUAUUAAAAGUCCAAGUUGGACUGACAAAACGCGUGCGGUGUCCUAGGGAUUGGUGGCAUAACCAGCGGUUUAAAAGCUGUGUAUAUCCGCAGCAAAUCACCGGAAAGCGGCGUUAUUAGCACCACAAAUUGAUGGUUGGUACGAGUACAAUUGCGCCGCAUAAAACCAGAGAUUCUACCCUCAAUCGGUUCUUAAGACGUACUGCGCGUUUCACCAGACCACAAUGCAGGGCGGCACCGUUAGGCAACACAACGAGACUACUCAUGCACAUAAGGAAGGUUAUCGCCAUAGACAUGGCGCGGCAGCGCAGAAUGUUUAAAUCUAAAUCUGGUAUGGGAGGCGUGCCCGUUGGUAUGAAGAAAUUUGCUGGGAGAAAAAGUCUAAGGCCUUGAAUCCGGCGGGUCUUAAUACUUACCUACAAAAUCAUCAGGCUGUACUUCCUGUAUC"
# s1 =       "........(((((((((((((((((..((((.((.((((((((((...((.((((((((((....(((.((((((.......)))))).....)))....))))..(((((...(((..(((((((....(((......((((.................((((((...))))))........((((((((((((.(((((((((((...........((((......((((((......)))))).....))))....(((((((((((.(.((((((......))..)))).).)))).....))))))).........(((........))).))))))))))).)))..)))))))))....)))).......))))))).)))..))).))))).))))..)).))..))).))).)))).))..))))..((((((....)))))).....)))).)))))............))))))))((((((((.(((.(.((((.........((((..(((((.....((.((((((((((....))..)))))))))).))).))..)))))))).).))))))..)))))....."
# s2 =       ".............((((((((((((.(((((.((.((((((((((...((.((((((((((...(((..((((((.......))))))....))).....))))..(((((...(((..(((((((....(((......((((.................((((((...))))))........((((((((((((.(((((((((((((.....))..((((......((((((......)))))).....))))....((((((((((((..((((........)))))))))...........))))))).........(((........))).))))))))))).)))..)))))))))....)))).......))))))).)))..))).))))).))))..)).))..))).))).)))).))..)))...((((((....))))))))...)))).))))))))............(((((((((((((.(((.((.(((((.((((..(((((((((.......(((...)))......))))))))))))).))))).))..............)))))..)))))))))))"

#100 nt with 2 inner sections
# sequence = "CGCAUCUCUUUAGGGUAUGAAAUGUUAUAUGCUACGGGAACAAUGCCGACCUUCGGAGACCUAAGGAAUACGUCUUUCGAGCGGAAGGAUUCCUCGUUCA"
# s1 =       ".(((((((..(((.((((((....)))))).))).))))....)))(((...)))..(((...((((....(((((((.....))))))))))).))).."
# s2 =       "....((.((.(((.((((((....)))))).))).))))......((((...)))).(((...(((((....((((((.....))))))))))).))).."

# 300 nt fails?
# sequence = "GUUGGGGUAGGGGCCCGACUAAUAGCAGUUCUUUGACAAUUCUUCUGCGUUAUUCAUUUUGAUAACAAUUAUUAUAUAAGUGCUGGAAAGCCAACCAUAUCGGCCUUAAUCCCCAACAGAAAACCUUCACGAGUGGGGUCGGUUUCGCAAUUCCUUGCGUGAGAGCCGAGGUUGGUUACAUGAGAAUAACCACGCCGAAC"
# s1 = "(((((((..((((((.((.((((.((((................)))))))).))((((..((((......))))..))))..(((....)))........))))))...))))))).......((.((((...((((..............)))).))))))...((..((.(((((........)))))))..))..."
# s2 = "((((((((((.((((.........((((................))))((((((......)))))).............(((.(((....)))..)))...)))))))..)))))))......(((........))).((((((((((((....)))))...)))))))(((((((((........)))))).)))...."




# dist0 , fix!
# sequence = "GGUAUGCUGCCUAUAAAGUAAGUGCACUGAUCUACAAGCAUAAGCGGGACGCUUCGUGAGUACAUUCAGA"
# s1 = "....((((((.......)))...)))((((..(((...((((((((...))))).))).)))...))))."
# s2 = ".(((((((........)))..)))).((((..(((..((..(((((...))))).))..)))...))))."

sequence = "UGAAGACCCAUUGAGUAACGACACCGCACGGCGCAUGGCGUCAGAGUAGCACUGCCUCGU"
s1 =       "....(((((((.(.(...((....)).....).))))).)))((.((......))))..."
s2 =       "....(((((((...((..((........))..)))))).))).(((..((...))))).."
s2 =       "....(((((((...((..((........))..)))))).)))((.((......))))..."

# no sections 300
# sequence = "UCACGACACCCCUCAACUAUAACAACGGUCCGUACAUACUAGCCCUGCAAUGGAACGGGCAGGGCCAGCACAGGUGGGGCGCCCGCUUGGGGGAUCAAAUGUGUGAUCCAGAUACUUUAGACGCGUGCAGAACUUUUUAGAUCGAUCAGUGGGAACAGGCAUUGAUUAUGAAAUCAAUUAGGGGGUUUAGGACCGCACCACAAACUGCGGGAGGGCACGCUUUGGUUCCUGUGUUACGCUAAUCCUCUAGCCACGGAGGGCUUCUUCGUACAAUGAUUGGGUUACCAGGGUUCCAGUGUG"
# s1 = ".........(((((.(((........)))............(((((((..........))))))).....((((((((...))))))))..((((((......))))))......(((((((.(.(........(((((....((((((((........)))))))).))))).......).).))))))).(((((........))))))))))((((((..((.(((((.((.((.((((((........((((((.....)))))).....)))))))).)))))))..))))))))"
# s2 = ".(((...((((........((((..(((((.((((...((.(((((((..........))))))).))..((((((((...))))))))((((((.....(((((((.(((..((..((((.((((((.(((((((((((((.((((((((........)))))))).....)))...))))))))))....(((((........)))))....))))))))))))..))).)))))))..)))))).((((......)))).....))))...))))).))))...))))....))).."

# 3 ext loops
# sequence = "AACCCGCGAUAAGUUGGCGCUUGUCCCACUCCGUAAACCUGUGUCUCUCAGGCGGUUACCCGAUAGAAGGCAGUAGGAUGUAUCACCCCCCACCGACUCCACUAUACAACGAAACCGGCCCUUUGUGAAUAAACGCUCUCAGUUAGAUUGAGCGGCGCAACAGAACCAGAUCGAUUCCCGGGGCAGUAGUGUCCUACUGCCCAACCCGGUAGAGGGAUCGAGUGGUCUAGGCGAAUGAUACGAGAGGCUCACUGGGGACGAUGGGAGUGAUCUACUCGAUGUUGCGUCCGAUCACAUCCA"
# s1 =       "....(((.(.....).))).(((((((............((.(((((((.(((((((.....((((..((.(((.((.((..........)))).))))).)))).......))))).)).............(((((...........))))).(((...(((.(((..(((((((((((((((((((....))))))))....))).....)))))))).))))))..)))........))))))).))...)))))))...(((((...)))))(((((...........))))).."
# s2 =       "....(((.(.....).))).(((((((............((.(((((((.(((((((...((((((..((.(((.((.((..........)))).))))).))))....)).)))).))).((((....))))(((((...........))))).(((...((.((((..(((((((((((((((((((....)))))))))....)).....)))))))).))))))..)))........))))))).))...)))))))((((.((((((....((......))...)))))).))))"

# sequence = "AUCAUAAUUACUAUAACGCUCUGCCGGGGGGAUGGCCGUCUUAACCACUGGCGAUAAACCCCCCAUCAAUUUGCUAAGGCUAGGGUUUUUUGUGUAAGCCGUCAUGAAGUAAAACCGUUUAACCCGGGUUUGGAGGUAUCCGGUGCGGAGCGACCCGAAAUCACUGAGGGUCCAGACGAUCGUGCUCACUGACUUACUUC"
# s1 = ".................((((((((((((((((.((((..........)))).))...)))))).......((((..(((...(((((.......))))))))....))))(((((..........)))))(((....))).).)))))))(((((...........)))))(((..((......)).)))........."
# s2 = ".................(((((((.((((((((.((((..........)))).))...))))))........((....))..((...(((((((........))))))).....)).......((((((........)))))).)))))))(((((...........))))).........((........))......."
#
# sequence = "CGCAUCUCUUUAGGGUAUGAAAUGUUAUAUGCUACGGGAACAAUGCCGACCUUCGGAGACCUAAGGAAUACGUCUUUCGAGCGGAAGGAUUCCUCGUUCA"
# s1 = ".((((.(((((((.((((((....)))))).))).))))...)))).....(((((((((.((.....)).))))).)))).((((...))))......."
# s2 = ".((((.(((((((.((((((....)))))).))).))))...))))((...(((((((((...........))).))))))))(((.((....)).)))."

# sequence = "UAAGGUCCCACACCCCCACGGGAUUUAUGCCUCAUGCGAGACUCGAUUAUCUGACGCUUAGGAGCCUAAGUGGCGGAGUCUCGCGCUGAACCGAGCACAG"
# s1 = "..(((((((..........))))))).........(((((((.......((((.(((((((....))))))).)))))))))))(((......)))...."
# s2 = "..((((......(((....)))......))))...(((((((.......((((.(((((((....))))))).)))))))))))(((......)))...."

# sequence = "CUAUCACCUGUACUUGCGGCCACCACCUCAUUUGCCCGUUGGCCGACUCAUCGCGUGGUAGGGGUUGUCGACUGUGGGAGUUUUCAGAUUACGUCCUAACCAACCGUAAAAGAACCUCCGGUAAGAGUAAAAGGAGGGUACAGGGUACAGUAAAGACGAUAGGACAGCCGUGUGUUAGGGCCGGAAGCCGCGGGCAAAUU"
# s1 = ".......((((....))))..........(((((((((.((((((.(((..((..((((((((((.(((((((.....))))....))).)).)))).))))..))........(((((.............))))).((((.(((...((....)).........))).))))...))).))...)))))))))))))."
# s2 = ".......((((....))))..........(((((((((((((((.......((((((((..(..((((((((((((((.((.........)).)))......(((.(.......(((((.............)))))....).)))))))).....))))))..).))))))))....)))))......))))))))))."

# s1 = ".......((((....))))..........(((((((((((((((.......(((.((((........(((((((.((..((.........))..)).......((.........(((((.............)))))......))..)))).....))).......)))).)))....)))))......))))))))))."

sequence = "CACAUGGGAAGCAAUGUAGCUCUCUGAUUCAGCUCUGCUACCAAAUAAUUGGCCAGGUCCAGCGCGAUCUCUCGGUUGGGCAACUAUGUCAGUGGAGGGU"
s1 = "...((.((((((......))).))).))...........(((.....((((((...(((((((.(((....))))))))))......))))))....)))"
s2 = ".....(((.(((......)))))).......(((((.((((........((((..(((...((.(((((....))))).)).)))..)))))))))))))"


search_width_multiplier = 0.000001
search_width = RNA.bp_distance(s1, s2)*search_width_multiplier

print ("sw:", search_width, "bp_dist:", RNA.bp_distance(s1, s2))

cmd = ""
if search_width:
    cmd = f'./main "{sequence}" "{s1}" "{s2}"'
if search_width_multiplier:
    cmd = f'./main "{sequence}" "{s1}" "{s2}"'

fc = RNA.fold_compound(sequence)
s1_eval = fc.eval_structure(s1)


bp_dist = RNA.bp_distance(s1, s2)

# print ("bp dist:", bp_dist, len(sequence))



start_findpath = time.time()
result = findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, True)
# result = findpath.init_merge_findpath(sequence, s1, s2, search_width_multiplier, True)
result = result/100.0
runtime = time.time()-start_findpath
print ("~~~~~~~~~~~")
print ("findpath single", runtime)
print (result, "barrier:", result-s1_eval)


# start_findpath = time.time()
# result = findpath.init_mfe_findpath(sequence, s1, s2, search_width_multiplier, True)
# result = result/100.0
# runtime = time.time()-start_findpath
# print ("~~~~~~~~~~~")
# print ("mfe merge", runtime)
# print (result, "barrier:", result-s1_eval)



# start_findpath = time.time()
# result = findpath_librna.pathfinder(sequence, s1, s2, search_width=search_width_multiplier*bp_dist)
# runtime = time.time()-start_findpath
# print ("~~~~~~~~~~~")
# print ("librna SWIG bindings:", runtime)
# print (result, "barrier:", result-s1_eval)




# result = findpath.init_merge_ext_findpath(sequence, s1, s2, search_width_multiplier, True)
# end_findpath = round(time.time()-start_findpath,4)
# result = result/100.0
# print (result, "barrier:", result-s1_eval)