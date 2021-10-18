#!/usr/bin/env python3
# coding: utf-8

# generic Python libraries
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from PIL import Image
from IPython.display import SVG, display
from collections import Counter
from collections import defaultdict
from dataclasses import dataclass
import itertools
import difflib
import sys
import os
import random
import string
import time
import logging

# coloredlogs
import coloredlogs

import RNA
from helper import p_to_s, print_moves
from findpath_i_sorted import find_path

sys.path.append('../../')
import findpath_librna
import findpath


processed_stacks = set()

def find_stack(sequence, s1, s2, path, moves):
    """
    find compensation stacks to lower the energy of e.g. saddle points
    """

    # print ("launch find_stack")

    fc = RNA.fold_compound(sequence)

    for max_pos in range(len(path)):

        # print ("start find stack", max_pos, len(path), moves)

        if max_pos <= 1:
            A = path[0]
            B = path[2]
            current_ptable = path[1]
        elif max_pos >= len(path)-1:
            A = path[-3]
            B = path[-1]
            current_ptable = path[-2]
        else:
            A = path[max_pos-1]
            B = path[max_pos+1]
            current_ptable = path[max_pos]

        # workaround...
        A[0] = len(A)-1
        B[0] = len(B)-1
        current_ptable[0] = len(A)-1

        A_str = p_to_s(A)
        B_str = p_to_s(B)
        C_str = p_to_s(current_ptable)

        # quick workaround - ignore previously computed stacks
        if (A_str, B_str, C_str) in processed_stacks:
            continue
        processed_stacks.add((A_str, B_str, C_str))

        # energy of current ptable
        base_energy = fc.eval_structure_pt(current_ptable)/100


        # compatible = A.copy()
        # offset_loop_table = A.copy()

        loop_A = RNA.loopidx_from_ptable(A)
        loop_B = RNA.loopidx_from_ptable(B)
        current_loop = RNA.loopidx_from_ptable(current_ptable)


        allowed = ["AU", "UA", "GC", "CG", "GU", "UG"]

        # generate indirect candidates for the current structure
        candidates = []
        for i in range(1, len(A)):
            # if not compatible[i]: continue
            # offset = offset_loop_table[i]

            if current_ptable[i] != 0:
                continue

            for j in range(i, len(A)):
                # requirements: minimum length between i and j
                if i+4 > j:
                    continue
                if current_ptable[j] != 0:
                    continue

                # if not compatible[j]: continue
                # if offset_loop_table[j]-offset != 0: continue
                # if current_loop[i] != current_loop[j]: continue

                # incompatible loops (with respect to bordering structures)
                # if loop_A[i] != loop_A[j]:
                #     continue
                # if loop_B[i] != loop_B[j]:
                #     continue

                # i and j have to be on the same loop
                if current_loop[i] != current_loop[j]:
                    continue

                # only consider allowed base pairs
                seq_i = sequence[i-1]
                seq_j = sequence[j-1]
                if seq_i+seq_j not in allowed:
                    continue

                # dont consider candidates which are in the primary move set
                if (i, j) in moves:
                    continue
                if (-i, -j) in moves:
                    continue

                candidates.append((i, j))

        def moves_to_en(moves, Verbose=False):
            # Verbose = True
            new_ptable = current_ptable.copy()
            for i, j in moves:
                new_ptable[i] = j
                new_ptable[j] = i
            if Verbose:
                print(p_to_s(new_ptable), end=" ")
                print(fc.eval_structure_pt(new_ptable)/100)
            return fc.eval_structure_pt(new_ptable)/100

        # if max_pos == 13:
        #     print (max_pos, p_to_s(current_ptable))
        #     for candidate in candidates:
        #         moves_to_en([candidate], Verbose=True)

        # candidates_with_energies = [[moves_to_en([(x[0], x[1])]), [(x[0], x[1])]] for x in candidates]
        # candidates_with_energies.sort(key=lambda x: x[0])

        # input is a list of tuples which can potentially added [(1, 63), (1, 65), (1, 67), (1, 69), (7, 11), (7, 12)]
        # lets find out which tuple combinations have potential to lower the energy
        combinations = [
            [moves_to_en([(x[0], x[1])]), [(x[0], x[1])]] for x in candidates]
        all_candidates = []

        # combinations.sort(key=lambda x: x[0])
        # combinations = combinations[:50]

        # try: rnasubopt for candidate generation
        # cmd = f'printf "{sequence}\n{C_str}\n" | RNAsubopt -C −−enforceConstraint -e 2'
        # result = subprocess.check_output(cmd, shell=True, encoding="utf8")
        # subopt_list = result.split("\n")
        # subopt_structure = subopt_list[1].split()[0]
        # subopt_ptable = RNA.ptable_from_string(subopt_structure)
        # candidates = []
        # for i in range(1, len(current_ptable)):
        #     if current_ptable[i] == subopt_ptable[i]: continue
        #     if subopt_ptable[i] == 0: continue
        #     j = subopt_ptable[i]
        #     if i>j: continue # closing bracket

        #     if current_loop[i] != current_loop[j]: continue
        #     if (i,j) in moves: continue
        #     if (-i,-j) in moves: continue
        #     candidates.append((i, j))

        # combinations = [[moves_to_en([(x[0],x[1])]),[(x[0],x[1])]] for x in candidates]
        # all_candidates = []

        iteration = 0
        while iteration < 10:
            iteration += 1
            # print (test)
            all_candidates += combinations
            next_combinations = []
            for current_energy, candidate_list in combinations:
                for next_candidate in candidates:
                    next_candidate_list = candidate_list.copy()

                    if next_candidate in candidate_list:
                        continue
                    start = candidate_list[-1][0]
                    end = candidate_list[-1][1]
                    # next candidate has to fit into the last one
                    if next_candidate[0] <= start or next_candidate[1] >= end:
                        continue
                    next_candidate_list += [next_candidate]
                    next_combinations.append(
                        [moves_to_en(next_candidate_list), next_candidate_list])

            if next_combinations == []:
                break
            next_combinations.sort()
            combinations = next_combinations[0:20]

        all_candidates.sort()

        # print ("bfs done")

        # for en, candidate in all_candidates[:4]:
        # print (max_pos, candidate, en, base_energy, candidate)
        # if en<base_energy:
        #     print (max_pos, candidate, en, base_energy, candidate)

        # if max_pos == 13:
        #     # print (candidates)
        #     for en, candidate in all_candidates[:10]:
        #         moves_to_en(candidate, Verbose=True)
        #         yield set(candidate)
        #         # print (max_pos, candidate, en, base_energy, candidate)

        if all_candidates != [] and all_candidates[0][0] < base_energy:
            # print (p_to_s(current_ptable))
            # print ("c:", all_candidates[1][-1], moves_to_en(all_candidates[1][-1], Verbose=True))
            # print ()
            # i = all_candidates[0][1][0][0]
            # j = all_candidates[0][1][0][1]
            # print ("candidate:", all_candidates[0][0]-base_energy, all_candidates[0][1], (i,j) in moves, (-i,-j) in moves)
            yield set(all_candidates[0][-1])

        # if all_candidates != [] and len(all_candidates)>1 and all_candidates[1][0]<base_energy:
            # # print ("c:", all_candidates[1][-1], moves_to_en(all_candidates[1][-1], Verbose=True))
            # yield set(all_candidates[1][-1])

    return


def find_moves(sequence, s1, s2, search_width=0, sw=0, add_moves=[], Verbose=True):

    # take direct paths and feed it into the indirct move finder

    if search_width==0:
        search_width = RNA.bp_distance(s1, s2) * sw
    if sw==0:
        sw = search_width / RNA.bp_distance(s1, s2)



    # direct_paths = find_path(sequence, s1, s2, add_moves=add_moves,
    #                   search_width=search_width, Debug=False, Verbose=Verbose)
    # max_en = direct_paths[0][0]
    # max_en_d = max_en


    add_moves2 = []
    for (i,j) in add_moves:
        add_moves2.append(i)
        add_moves2.append(j)

    # max_en_i = findpath.init_single_findpath_i(sequence, s1, s2, sw, add_moves2, True, Verbose=True)  

    fp = findpath.findpath_single_i(sequence, s1, s2, search_width_multiplier=sw, mp=True)

    max_en = max_en_d = fp.get_en()/100.0
    direct_paths = [fp.get_path()]


    print (fp.get_path())

    print ("max en d:", max_en)

    indirect_candidates = set()

    for moves in direct_paths:               

        # if s < max_en:
        #     continue

        # print_moves(sequence, s1, s2, moves, Verbose=Verbose)
        # max_en = s        
                
        current_ptable = list(RNA.ptable(s1))
        ptables = [current_ptable.copy()]
        for pos, (i, j, en) in enumerate(moves):
            if (i,j) == (0,0):
                continue
            if i < 0:
                current_ptable[-i] = 0
                current_ptable[-j] = 0
            else:
                current_ptable[i] = j
                current_ptable[j] = i
            ptables.append(current_ptable.copy())
        

        for indirect_moves in find_stack(sequence, s1, s2, ptables, moves):
            
            if add_moves != [] or add_moves != set():
                # if already previous indirect moves are present, add them to
                # the newly found moves
                indirect_moves = indirect_moves | set(add_moves)

            indirect_candidates.add(tuple(indirect_moves))

            # indirect_candidates |= indirect_moves
        
    # fp indirect with our best candidates 
       
    print (len(indirect_candidates))

    # iterate over all indirect paths, return the one with the lowest energy which also
    # has the lower amout of indirect moves

    candidates = []

    for indirect_moves in indirect_candidates:

        print (indirect_moves)

        # indirect_paths = find_path(sequence, s1, s2, add_moves=indirect_moves,
        #                 search_width=search_width, Debug=False, Verbose=Verbose)
     
        add_moves2 = []
        for (i,j) in indirect_moves:
            add_moves2.append(i)
            add_moves2.append(j)


        if len(add_moves2) > 16:
            continue

        
        print ("try", add_moves2)
        try:
            fp = findpath.findpath_single_i(sequence, s1, s2, add_moves=add_moves2, search_width_multiplier=sw, mp=True)
            s = fp.get_en()/100.0
            indirect_paths = [fp.get_path()]
        except:
            s = max_en_d
            indirect_paths = []

        print ("try2")
        # print (s, add_moves2)

        for moves in indirect_paths:            
            
            # check if our indirect path has actually a lower saddle energy
            if s >= max_en_d:
                indirect_moves = ()

            moves = [(i[0], i[1]) for i in moves]

            count = len(moves)
            candidates.append([s, count, indirect_moves, search_width, moves])
            
            # print (s)
            # print_moves(sequence, s1, s2, moves, Verbose=Verbose)

    candidates.sort(key=lambda x: (x[0], x[1])) # sort by max_en, count

    return candidates





if __name__ == '__main__':

    sequence   = "CCAGCGUAUUAGUUAUGGCCUGGAGGUAGAAGCGUUAGAGCAAUACUUCUACAGAGACCACGUGAGGUAG"
    s1         = "((((..((.......))..))))..(((((((.............))))))).................."
    s2         = "((((((((.....))).)).)))..(((((((.((.......)).)))))))....(((......))).."

    # Tabu paper example
    sequence = "CGCGACGGCUACGCGACGGCAAUGCCGUUGCGAAGCCGUCGCGAUC"
    s1 = "(((((((((..............))))))))).............."
    s2 = "...........(((((((((..............)))))))))..."

    # rna2dfold example
    # sequence = "GGGCGCGGUUCGCCCUCCGCUAAAUGCGGAAGAUAAAUUGUGUCU"
    # s1 = "(((((.....)))))(((((.....)))))(((((.....)))))"
    # s2 = "((((((((((.....(((((.....))))).....))))))))))"

    # dsrA
    sequence = "ACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUAAGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUU"
    s1 = "..(((((((.....)))))))...(((((((((((((.......))))))).)))))).((((((((((.....))))))))))."
    s2 = "..(((((((.....)))))))................(...(((((....)))))...)((((((((((.....))))))))))."

    sequence = "GGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAAUUUUAUUAAAAGUCCAAGUUGGACUGACAAAACGCGUGCGGUGUCCUAGG"
    s1 = "...(((.((((..(((((.(((...................)))................((((((...))))))........))))).)).).).)))."
    s2 = "((...(((((...(((((.(((...................)))................((((((...))))))........))))).)))))))...."


    search_width = 500
    Verbose = True
    # Debug = True
    Debug = False

    # add_moves = [(9, 15), (36, 42)]



    candidates = find_moves(sequence, s1, s2, sw=2, Verbose=True)
    s, count, indirect_moves, search_width, moves = candidates[0]
    print (indirect_moves)
    print_moves(sequence, s1, s2, moves, Verbose=Verbose)


    processed_stacks = set()


    # second iteration

    # candidates = find_moves(sequence, s1, s2, search_width=500, Verbose=True, add_moves=[(9, 15)])
    # s, count, indirect_moves, search_width, moves = candidates[0]
    # print (indirect_moves)
    # print_moves(sequence, s1, s2, moves, Verbose=Verbose)









    # paths = find_path(sequence, s1, s2, add_moves=add_moves,
    #                   search_width=search_width, Debug=Debug, Verbose=Verbose)
    # for s, sw, mode, moves in paths:
    #     print_moves(sequence, s1, s2, moves, Verbose=Verbose)
    #     break
    # print ([(i[0], i[1]) for i in paths])



    # print (se)
    # print('orig findpath:')

    # pathfinder_result = pathfinder.pathfinder(
    #     sequence, s1, s2, search_width=search_width, verbose=Verbose)

    # print(pathfinder_result.max_en)
