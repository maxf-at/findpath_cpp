#!/usr/bin/env python3
# coding: utf-8


import RNA
import time
from dataclasses import dataclass

import matplotlib.pyplot as plt

def str_compare(s1, s2):
    for i, (c1, c2) in enumerate(zip(s1, s2)):
        if c1 != c2:
            if c1 == ".":
                yield i
            else:
                yield -i

def pathfinder(sequence, s1, s2, verbose=False, output=False, search_width = None, return_paths=False, return_s=False):

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

    moves_pos = []

    current_structure = ""

    for i, step in enumerate(paths):    

        if return_paths:
            moves = [i for i in str_compare('.'+current_structure, '.'+step.s)]
            if moves==[]:
                if return_s:
                    moves_pos.append((0,0,step.en, step.s))
                else:
                    moves_pos.append((0,0,step.en))
            else:
                if return_s:
                    moves_pos.append((moves[0],moves[1],step.en, step.s))
                else:
                    moves_pos.append((moves[0],moves[1],step.en))

        current_structure = step.s

        if step.en > sE:
            sE = step.en
            s_pos = i

    if return_paths:
        return sE, moves_pos

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
