#!/usr/bin/env python3
# coding: utf-8


import RNA
import time
from dataclasses import dataclass

import matplotlib.pyplot as plt



def pathfinder(sequence, s1, s2, verbose=False, output=False, search_width = None, section = None):

    # prune off sections which will currently not be regarded / later merged
    if section:
        if len(section)==2:
            start = section[0]
            end = section[1]
            s1 = '.'*start + s1[start:end] + '.'*(len(s1)-end)
            s2 = '.'*start + s2[start:end] + '.'*(len(s1)-end)
        if len(section)==4:
            # print ('section 4', section)
            start = section[0]
            mid1 = section[1]
            mid2 = section[2]
            end = section[3]
            s1 = '.'*start + s1[start:mid1] + '.'*(mid2-mid1) + s1[mid2:end] + '.'*(len(s1)-end)
            s2 = '.'*start + s2[start:mid1] + '.'*(mid2-mid1) + s2[mid2:end] + '.'*(len(s1)-end)
            # print (s1)
            # print (s2)

    start = time.time()
    # Set model details (by Stefan)
    md = RNA.md()
    md.temperature = 37.0
    # md.temperature = 22.0
    # How to treat "dangling end" energies for bases adjacent to helices in free ends and multi-loops.
    # md.dangles = 0
    md.dangles = 2
    md.logML = 0
    md.special_hp = True
    md.noGU = False
    md.noGUclosure = False
    fc = RNA.fold_compound(sequence, md)
    # fc = RNA.fold_compound(sequence)

    bpd = RNA.bp_distance(s1, s2)

    if verbose:
        print(f"Base-pair distance: {bpd}")


    max_energy = None  # Specify upper bound for barrier energy (kcal/mol).
    
    # search_width = None  # Adjust upper bound for findpath search.

    if search_width is None:
        search_width = 2 * bpd

    paths = fc.path_findpath(s1, s2, width=search_width)

    sE = float('-inf')     # saddle energy
    s_pos = 0
    en_min = float('inf')  # minimum, assuming we get below E1

    for i, step in enumerate(paths):
        current_structure = step.s
        if step.en > sE:
            sE = step.en
            s_pos = i

    return sE


if __name__ == '__main__':
    # 150 nt seq with 2 frozen parts, 1 non-frozen mid section
    sequence = "CGUAAAGUAAGCGCAGGACUAGAGGAAGGCUGGUGCGCGAUAGGCCCCCGAUCAUUGAUCACGUCAAGGACACUCUUGUUGAUGAAGCUAAGGGGCGCCGCCUCCCAUCGCUUUGAAUGACGGCGCAAUGAGGCCCGGAUAUAACUCGGG"
    s1 = "......((...(((.(.(((((.......))))).))))....)).(((((((...))))....(((((....)))))..((.........(((((((((.(...((......))...).)))))).......))).........)))))"
    s2 = "..(.......((((...(((((.......)))))))))(((.((...)).))).(((((...))))).)...((((((((.....))).)))))((((((.(...((......))...).)))))).......(((((.......)))))"

    path = pathfinder(sequence, s1, s2, verbose=True)
    print (path)
    # plt = plot_paths(path)
    # plt.show()

    # print(path)

    # sE, s_pos, b1, b2, e1, e2, moves_str, moves_en = pathfinder(sequence, s1, s2)
    # print(moves_str)
