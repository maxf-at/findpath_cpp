#!/usr/bin/env python3
# coding: utf-8

import time
import subprocess
import RNA

import findpath
import findpath_librna


sequence = "CGGGAGCGAGUAUGGACCGGGUUGAAUUGGUACCUCGGCCGCCCUCUGGGGGUAGCCGACGGGCGGCUUCAGCCGGGCCC"
s1 = ".............((.((.(((((((.(.((.(((((((.(((((...))))).))))).)))).).))))))).)).))"
s2 = "((....)).....((.((.(((((((...........(((((((((.(.(.....)))).)))))))))))))).))))."

intermediates = [

'((....)).....((.((.(((((((...........((((((((...))))).(((....))))))))))))).)))).',
'((....)).....((.((.(((((((...........((((((((...))))...((....))))))))))))).)))).',
'((....)).....((.((.(((((((...........(((((((....)).)...........))))))))))).)))).',
'((....)).....((.((.(((((((...........(((((((.....)))...((....))))))))))))).)))).',
'((....)).....((.((.(((((((...........((((((.......))...((....))))))))))))).)))).',
'((....)).....((.((.(((((((...........(((.((((...))))..(((....))))))))))))).)))).',
'((....)).....((.((.(((((((...........((.(((((...))))).(((....))).))))))))).)))).',
'((....)).....((.((.(((((((...........(.((((((...))))).(((....))).))))))))).)))).',

# '................((.(((((((..............((((((...((....))...))).)))))))))).))...',
# '................((.(((((((..............((((.(((..((...))..)))).)))))))))).))...',
# '................((.(((((((..............((((.(((.((....))..)))).)))))))))).))...',
# '................((.(((((((..............((((((.(..((...))..)))).)))))))))).))...',
# '................((.(((((((..............((((((.(.((....))..)))).)))))))))).))...',
# '................((.(((((((..............(((((((..((....)).).))).)))))))))).))...',
# '................((.(((((((..............(((((((.(.(....)).).))).)))))))))).))...',
# '................((.(((((((...........(((.((....)).)))..............))))))).))...',
# '................((.(((((((...........(((.(((...))))))..............))))))).))...',
# '................((.(((((((......((.((...(((((....)))..))...))))....))))))).))...',
# '................((.(((((((......((.((...((.(((...)))..))...))))....))))))).))...',
# '................((.(((((((............................((((.....))))))))))).))...',
# '................((.(((((((......((((((.......))))))................))))))).))...',
# '................((.(((((((.((............(((....)))............))..))))))).))...',
# '((....)).....((.((.(((((((...((.(.(((((.(((((...))))).))))).).))...))))))).)))).',
# '((....)).....((.((.(((((((...(..(.(((((.(((((...))))).))))).)..)...))))))).)))).',
# '((....)).....((.((.(((((((..(.(.(((((((.(((((...))))).))))..))).).)))))))).)))).',
# '((....)).....((.((.((((((((((...(.(((((.(((((...))))).))))).)..))).))))))).)))).',
# '((....)).....((.((.(((((((..(.(.(((.(((.(((((...))))).)))...))).).)))))))).)))).',
# '((....)).....((.((.(((((((..(.(.((.((((.(((((...))))).))...)))).).)))))))).)))).',
# '((....)).....((.((.(((((((..............(((((...))))).((((.....))))))))))).)))).',
# '((....)).....((.((.(((((((...........(.((((((...))))).(((....))).))))))))).)))).',
# '((....)).....((.((.(((((((...........((.(((((...))))).(((....))).))))))))).)))).',
# '((....)).....((.((.(((((((...........((((((((...))))).(((....))))))))))))).)))).',
# '((....)).....((.((.(((((((...........(((.((((...))))..(((....))))))))))))).)))).',
# '((....)).....((.((.(((((((...........((((((((...))))...((....))))))))))))).)))).',
# '((....)).....((.((.(((((((...........(((((((.....)))...((....))))))))))))).)))).',
# '((....)).....((.((.(((((((...........((((((.......))...((....))))))))))))).)))).',
# '((....)).....((.((.(((((((...........(((((((....)).)...........))))))))))).)))).',
# '((....)).....((.((.(((((((...........((((.(((....((....))...)))))))))))))).)))).',
# '((....)).....((.((.(((((((...........(((((.((....((....))...)))))))))))))).)))).',
# '((....)).....((.((.(((((((...........((((((......((....))....))))))))))))).)))).',
# '((....)).....((.((.(((((((...........(((((((.....((....))...)))))))))))))).)))).',
# '((....)).....((.((.(((((((...........(((((((.(...((....)))..)))))))))))))).)))).',
# '((....)).....((.((.(((((((...........(((((((((...((....)))).)))))))))))))).)))).',
# '((....)).....((.((.(((((((...........(((((((((.((......)))).)))))))))))))).)))).',
]


# sequence = "CUGGGUCGUGGUGGCCUUUUAGAUACGAUUCACGAACGUAGCACGUUUCGGUCUCCGGAGACGCAAUGAUCUCGAGGGUA"
# s1 = ".(((((((((.(((....)))..))))))))).((.(((.((..(((((((...))))).)))).))).))........."
# s2 = ".(((((((((.............)))))))))(((((((...))))).)).((((.(((..........))).))))..."
# intermediates = [
# '.((((((((((...))........))))))))(((((((...)))).)))..(((..((((........))))..)))..', 
# '.((((((((((...))........))))))))(((((((...))).))))..(((..((((........))))..)))..', 
# '.((((((((((...))........)))))))).((((((...))))))....(((..((((........))))..)))..', 
# '.(((((((((.(((....)))..)))))))))(((((((...))))).))..(((..((((........))))..)))..', 
# '.(((((((((.(((....)))..)))))))))(((((((...)))).)))..(((..((((........))))..)))..', 
# '.(((((((((.(((....)))..)))))))))(((((((...))).))))..(((..((((........))))..)))..', 
# '.(((((((((.(((....)))..))))))))).((((((...))))))....(((..((((........))))..)))..', 
# '.(((((((((.....((...)).)))))))))(((((((...)))).)))..(((..((((........))))..)))..', 
# '.(((((((((.............)))))))))(((((((...))))).))..(((..((((........))))..)))..', 
# '.(((((((((.............)))))))))(((((((...)))).)))..(((..((((........))))..)))..', 
# '.(((((((((.............)))))))))(((((((...))).))))..(((..((((........))))..)))..', 
# '.(((((((((.............))))))))).((((((...))))))....(((..((((........))))..)))..', 
# '.((((((((.((..........))))))))))(((((((...))))).))..(((..((((........))))..)))..', 
# '.((((((((.((..........))))))))))(((((((...)))).)))..(((..((((........))))..)))..', 
# '.((((((((.((..........))))))))))(((((((...))).))))..(((..((((........))))..)))..', 
# '.((((((((.((..........)))))))))).((((((...))))))....(((..((((........))))..)))..', 
# ]

sequence = 'AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAA'
s1       = '.............((((((..((...(..(((.....)))..).))..))))))(((((((....(((.((((((.......)))))).....)))....)))))))(((......((.((.........)).))......)))......'
s2       = '......(((.((((((...(((....))).(((((((((......))))))...)))((...(((....((((((.......))))))...))).)))))))).)))(((......((.((.........)).))......)))......'
intermediates = ['......(((.((((((.(((((....)))((((((((((......))))))...))))....(((....((((((.......))))))...))))).)))))).)))(((......((.((.........)).))......)))......', '......(((.((((((...(((....)))((((((((((......))))))...))))....(((....((((((.......))))))...)))...)))))).)))(((......((.((.........)).))......)))......', '......(((.((((((.............((((((((((......))))))...))))....(((....((((((.......))))))...)))...)))))).)))(((......((.((.........)).))......)))......', '..........((((((...(((....)))((((((((((......))))))...))))....(((....((((((.......))))))...)))...))))))....(((......((.((.........)).))......)))......']


sequence = 'GGAACAAGCAGCUGCAACAAGUUUGAACGGUUUCUCCUGCAGAAUGAGUGGUACCCUAACAUGCUAGGCUGCACUGGAAGCAAUGUCCCAUUGACCACUCCGCUCUCCCCUCUGGGAGUCGUGGGCCUGACCCCGAACUCACAUGAUCUA'
s1       = '.......((....))....((((((...((((...(((((.....(((((((..((((......)))).(((.......)))...........)))))))....(((((....)))))..)))))...)))).))))))...........'
s2       = '..............((...((((((....((((((..(((((..((.(((.((...)).)))..))..)))))..))))))...(((.......((((...((((((......)))))).))))....)))..))))))....)).....'
intermediates = ['.......((....))..((((((((...((((...(((((.((..(((((((.......(((((......))).))....(((((...))))))))))))....(((((....))))))))))))...)))).))))))....)).....', '.......((....))....((((((...((((...(((((.((..(((((((...(((...(((......))).)))...(((((...))))))))))))....(((((....))))))))))))...)))).))))))...........', '.......((....))....((((((...((((...(((((.((..(((((((...((..(((((......))).))..))(((((...))))))))))))....(((((....))))))))))))...)))).))))))...........', '.......((....))....((((((...((((...(((((.....(((((((.........(((......))).......(((((...))))))))))))....(((((....)))))..)))))...)))).))))))...........', '.................((((((((...((((...(((((.((..(((((((...((..(((((......))).))..))(((((...))))))))))))....(((((....))))))))))))...)))).))))))....)).....']

sequence = 'AUCAUUGCCUAAGGGACAUUUCCACGUGUAAGAAACCCGGAAACAUCGUUUGGCUUAGAAUGAACAGAUAGAUUGCCUACAACAAAAUACUGAUGACGGGGGAGCCGUGACAUGGCGGCCACCGCUAGUCGGCCUACGCUUUAAAACAUA'
s1       = '.((((.(((.(((.((..(((((..((.......))..)))))..)).))))))..........(((.((..((..........)).)))))))))((..((.((((.....((((((...))))))..)))))).))............'
s2       = '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).(((.....(((.......)))....)))..((((.((..(((((......)))))..)).)..)))(((....)))..........'
intermediates = ['.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................', '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).........(((....))).......(((((..(((((..(((((...)))))..)).)))...)))))..................', '.((((((((.(((.((..(((((..(.((.....))).)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................', '.((((((((.(((.((..(((((...............)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................', '.((((((((.(((.((..(((((...............)))))..)).))))))....)))))..........................(((((..(((((..(((((...)))))..)).)))...)))))..................', '......(((.(((.((..(((((..((.......))..)))))..)).))))))..............(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................']


sequence = 'AAAACGGCUAGACGUUUAUGACUGGCGAUAAUUAUUUGCAUGAUGGAUCCAAUCCAAUGCAGUCGAGUGAAGUACCGUCGUAGUAGCUAUAACCACCUUAUUCUUGCCCCCUGCAAGCGCUAGAUACUGAUCCAUGCUAAUAAAGUAGUC'
s1       = '.....(((((.....(((((((.((.....(((..((((((..(((((...)))))))))))...)))......))))))))))))))...((.((.((((((((((.....))))).((..(((....)))...)).))))).)).)).'
s2       = '.(((((......)))))....((((((........((((((..(((((...)))))))))))....(((.(((..(......)..))).....)))......(((((.....))))))))))).....(((...((((.....)))))))'
intermediates = ['.(((((......)))))((..((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).)))))))).........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).))))))..(((.......(((.....)))))).', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).))))))....((...))((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).)))))).....(((...((((.....)))))))', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).))))))...........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..(((((.((....((.((((...)))).))...)))))))((((.....)))).))))))...........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))..((((.(((....((.((((...)))).))...)))))))((((.....)))).))))))...........((((.....))))...', '.(((((......)))))....((((((........((((((..(((((...)))))))))))...(((((.((...((.((((...)))).)).)).)))))(((((.....)))))))))))...........((((.....))))...', '.(((((......))))).....(((((........((((((..(((((...)))))))))))..((((((.((...((.((((...)))).)).)).))))))((((.....)))).)))))(((....)))..((((.....))))...']

sequence = 'AUCAUUGCCUAAGGGACAUUUCCACGUGUAAGAAACCCGGAAACAUCGUUUGGCUUAGAAUGAACAGAUAGAUUGCCUACAACAAAAUACUGAUGACGGGGGAGCCGUGACAUGGCGGCCACCGCUAGUCGGCCUACGCUUUAAAACAUA'
s1 = '.((((.(((.(((.((..(((((..((.......))..)))))..)).))))))..........(((.((..((..........)).)))))))))((..((.((((.....((((((...))))))..)))))).))............'
s2 = '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).(((.....(((.......)))....)))..((((.((..(((((......)))))..)).)..)))(((....)))..........'
intermediates = [
    '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    # '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).........(((....))).......(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '.((((((((.(((.((..(((((..(.((.....))).)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '.((((((((.(((.((..(((((...............)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '.((((((((.(((.((..(((((...............)))))..)).))))))....)))))..........................(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '......(((.(((.((..(((((..((.......))..)))))..)).))))))..............(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '.((((((((.(((.((..(((((..(.((.....))).)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '.((((((((.(((.((..(((((...............)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '.((((((((.(((.((..(((((...............)))))..)).))))))....)))))..........................(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '......(((.(((.((..(((((..((.......))..)))))..)).))))))..............(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '.((((((((.(((.((..(((((..((.......))..)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '.((((((((.(((.((..(((((..(.((.....))).)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '.((((((((.(((.((..(((((...............)))))..)).))))))....))))).....(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '.((((((((.(((.((..(((((...............)))))..)).))))))....)))))..........................(((((..(((((..(((((...)))))..)).)))...)))))..................',
    '......(((.(((.((..(((((..((.......))..)))))..)).))))))..............(((.....)))..........(((((..(((((..(((((...)))))..)).)))...)))))..................',

]


# sequence = "UGAAGACCCAUUGAGUAAAA"
# s1       = "(((((((.....)))).)))"
# s2       = "(((((.........)).)))"
# intermediates = [
#     "(((((.........)).)))",
#     "((((((.......))).)))",
#     "((((((((...))))).)))"]


# print (sequence)
# print (s1)
# print (s2)
print ("~~~~~~~~~~~")




target_s = s1
# target_s = s2

# mp = False
mp = True
search_width_multiplier = 2
target_en = findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp)
print (f'S target: {target_en:2.2f}')


# target_en = 999
search_width_multiplier = 2

start_findpath = time.time()
result = findpath.init_multi_findpath(sequence, target_s, intermediates, search_width_multiplier*2, target_en, mp)
result = result/100.0
runtime = time.time()-start_findpath
print ("~~~~~~~~~~~")
print (f'multi section')
print (f'S: {result:2.2f}, runtime: {runtime:2.4f} s')

search_width_multiplier = 2

print ("comparison individual findpath calls")

fp = findpath.findpath_class(sequence, mp)
    

start_findpath = time.time()
for dest in intermediates:
    result = fp.init(target_s, dest, search_width_multiplier)
    # print (result)
    # result = findpath.init_single_findpath(sequence, target_s, dest, search_width_multiplier, mp, target_en)
    # result = findpath.init_single_findpath(sequence, target_s, dest, search_width_multiplier, mp)
    print (f'{dest}, S: {result:2.2f}, bp_dist: {RNA.bp_distance(target_s, dest)}, sw: {search_width_multiplier*RNA.bp_distance(target_s, dest)}')
runtime = time.time()-start_findpath
print (f'runtime: {runtime:2.4f} s')

# findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp)


# start_findpath = time.time()
# result = findpath.init_single_findpath(sequence, s1, s2, search_width_multiplier, mp)
# result = result/100.0
# runtime = time.time()-start_findpath
# print ("~~~~~~~~~~~")
# print (f'single section')
# print (f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}, runtime: {runtime:2.4f} s')

