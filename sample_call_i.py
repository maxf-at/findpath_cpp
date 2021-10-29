#!/usr/bin/env python3
# coding: utf-8

import time
import RNA

import findpath
import findpath_librna

from pretty_print_path import print_moves

sequence = "GGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAAUUUUAUUAAAAGUCCAAGUUGGACUGACAAAACGCGUGCGGUGUCCUAGG"
s1       = "...(((.((((..(((((.(((...................)))................((((((...))))))........))))).)).).).)))."
s2       = "((...(((((...(((((.(((...................)))................((((((...))))))........))))).)))))))...."

add_moves = [8, 92, 9, 91]


# Address boundary error
sequence = "CUCACGAUAUGGGGCACGGUGCAUACAUACAAGAAUUACUUGGCACCCCCGCGGAAAGUAGUCUCCUUGUACGCAUAUCGCUACAAAAUCGGGGAACUUC"
s1 = "(((.((((.((.(((...((((.((((....(((.((((((.((......)).))..)))))))...)))).))))...))).))..)))))))......"
s2 = "...((....((.((...(((((.......((((.....))))))))).)).))....))((((((((((((.((.....)))))))....)))).))).."

add_moves = [53, 66, 54, 65, 10, 85, 19, 78, 14, 81, 20, 77, 15, 80, 55, 64, 18, 79]
add_moves = [53, 66, 54, 65, 10, 85, 19, 78, 14, 81, 20, 77, 15, 80, 55, 64]


sequence = "GGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAAUUUUAUUAAAAGUCCAAGUUGGACUGACAAAACGCGUGCGGUGUCCUAGG"
s1 = "...(((.((((..(((((.(((...................)))................((((((...))))))........))))).)).).).)))."
s2 = "((...(((((...(((((.(((...................)))................((((((...))))))........))))).)))))))...."

sequence = "GGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAAUUUUAUUAAAAGUCCAAGUUGGACUGACAAAACGCGUGCGGUGUCCUAGG"
s1 = "...(((.((((..(((((.(((...................)))................((((((...))))))........))))).)).).).)))."
s2 = "((...(((((...(((((.(((...................)))................((((((...))))))........))))).)))))))...."

add_moves = [6, 97, 8, 94, 5, 98, 4, 99, 12, 90]


# settings
search_width_multiplier = 2
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




# call with search width multiplier:
fp = findpath.findpath_single_i(sequence, s1, s2, add_moves=add_moves, search_width_multiplier=search_width_multiplier, mp=mp)


result = fp.get_en()/100.0
path = fp.get_path()

print (path)

print_moves(sequence, s1, s2, path, convert_to_float=True)
print(f'S: {result:2.2f}, barrier: {result-s1_eval:2.2f}')

