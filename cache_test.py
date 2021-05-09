#!/usr/bin/env python3
# coding: utf-8

import time
import subprocess
import RNA

import findpath
import findpath_librna

search_width_multiplier = 4
mp = True
start_findpath = time.time()

sequence = 'UCCGACAUUAAGACACACCAGGGUCUCGUAUCCCUAGGGUAAGGUACGCGCGGACCGGCCAAUCGGGUAUUGCUGCAAACUAUGGCAAUAGUGCAUAGGUUCAGACGAAGUACGGGUGGAUAUUUGUAGCCAGUAUGCUGGGUCUCCGGG'

fp = findpath.findpath_class(sequence, mp)

s1       = '((((..........((.((((((........)))).))))..........))))((((...(((.(((((.((((((((((((.(((....)))))))(((((..((.....))..))))).))))))))....))))).)))..)))).'
s2       = '((((..........((.((((((........)))).))))..........))))((((....((((((((((((((((((((((((....)).)))))(((((..((.....))..))))).)))))))).))))).))))....)))).'
result = fp.init(s1, s2, search_width_multiplier)
print (result)

s1       = '((((..........((.((((((........)))).))))..........))))((((...(((.(((((.((((((((((((.(((....)))))))(((((..((.....))..))))).))))))))....))))).)))..)))).'
s2       = '((((....((....((.((((((........)))).))))....))....))))((((....((((((((((((((((((((((((....)).)))))(((((..((.....))..))))).))))))).)))))).))))....)))).'
result = fp.init(s1, s2, search_width_multiplier)
print (result)


print ("~~~~~~~~~~~")

