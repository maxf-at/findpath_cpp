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

# custom
import RNA
from helper import p_to_s, print_moves



def print_d(*args):
    """
    debug print for coloredlogs library
    """
    msg = ""
    for m in args:
        msg += str(m)
        msg += " "
    logging.debug(msg)


# https://www.reddit.com/r/Python/comments/27crqg/making_defaultdict_create_defaults_that_are_a/
class key_dependent_dict(defaultdict):
    def __init__(self, f_of_x):
        super().__init__(None)  # base class doesn't get a factory
        self.f_of_x = f_of_x  # save f(x)

    def __missing__(self, key):  # called when a default needed
        ret = self.f_of_x(key)  # calculate default value
        self[key] = ret  # and install it in the dict
        return ret


@dataclass
class Intermediate:
    p_table:      list
    saddle_e:     int
    current_e:    int
    mode:         bool  # 1 = fwd, 0 = bwd
    moves:        list
    energies:     list
    opt:          float
    add_moves:    list
    # distance:     int


class Fp_class():
    def __init__(self, sequence, s1, s2, add_moves=[]):

        self.fc = RNA.fold_compound(sequence)
        self.sequence = sequence
        self.s1 = s1
        self.s2 = s2

        # energy caching
        self.evals = key_dependent_dict(
            lambda x: round(self.fc.eval_structure(x), 2))
        self.p_tables = key_dependent_dict(
            lambda x: RNA.ptable(x))
        self.l_tables = lambda x: RNA.loopidx_from_ptable(x)
        self.bp_dist = lambda x, y: RNA.bp_distance(x, y)

        # additional class variables for indirect paths
        self.processed_stacks = set()
        self.moves_add = add_moves




    def find_moves(self, s_ptable, t_ptable):
        """
        generator function, yields possible structures 1 move away
        from the original structure by finding fitting i and j with
        RNA pair and loop tables
        s_ptable: current ptable
        t_ptable: s2 end ptable
        """
        # loop table
        ls = self.l_tables(s_ptable)

        # extra indirect moves: i and j have to be unpaired, and on the same loop
        for i, j in self.moves_add:
            if i > 0 and s_ptable[i] == 0 and s_ptable[j] == 0 and ls[i] == ls[j]:
                yield i, j

        for i in range(len(s_ptable)):
            if i == 0:
                continue

            if s_ptable[i] == 0 and t_ptable[i] > i:
                j = t_ptable[i]
                # found j has to be empty and currently on the same loop as i
                if s_ptable[j] == 0 and ls[i] == ls[j]:
                    yield i, j
            # test for bp removal: i has to be paired with a different j in s2
            j = s_ptable[i]
            # dont remove things which are present in s2
            if s_ptable[i] > i and s_ptable[i] != s_ptable[j] and\
                    s_ptable[i] != t_ptable[i] and s_ptable[j] != t_ptable[j]:
                yield -i, -j

    def find_path_once(self, s1, s2, max_energy, width, mode=True, sort_min=False, Debug=False, Verbose=False):
        """
        main findpath algorithm (bounded BFS)
        """

        if Debug:
            coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s %(message)s'
            coloredlogs.install(level='DEBUG')

        # caching
        evals = self.evals
        p_tables = self.p_tables
        l_tables = self.l_tables
        bp_dist = self.bp_dist
        # s1 = fp_class.s1
        # s2 = fp_class.s2

        e1 = self.evals[s1]
        e2 = self.evals[s2]

        runtimes = 0
        current_bp = 0
        current_bp_end = bp_dist(s1, s2)      # current_bp_end = 4
        # paths = [(list(p_tables[s1]), float("-inf"), evals[s1], [])]

        s1_en = round(evals[s1], 2)
        end_p_table = self.p_tables[s2]
        init_intermediate = Intermediate(p_table=list(p_tables[s1]), mode=mode, saddle_e=float(
            "-inf"), current_e=evals[s1], moves=[(0, 0, s1_en)], energies=[0], opt=0, add_moves=[])
        # initial path start with 1 intermediate
        init_path = [init_intermediate]
        # paths start with 1 initial path
        paths = [init_path]

        # dont stop at current_bp_end, consider potential indirect moves
        while (current_bp != current_bp_end+2*len(self.moves_add)):

            # collect all new paths here (next iteration)
            collect_paths = []

            for current_path in paths:

                collect_paths_per_move = []
                add_moves = 0
                delete_moves = 0

                current_p_table = current_path[-1].p_table
                current_e = current_path[-1].current_e
                current_s = current_path[-1].saddle_e
                current_string = p_to_s(current_p_table)
                current_moves = current_path[-1].moves
                current_energies = current_path[-1].energies

                # "try_moves"
                for i, j in self.find_moves(current_p_table, end_p_table):

                    if i<0:
                        delete_moves += 1
                    else:
                        add_moves += 1

                    current_add_moves = current_path[-1].add_moves.copy()

                    if (i, j) in current_add_moves:
                        continue  # this optional move is already in the path

                    if (i, j) in self.moves_add:
                        current_add_moves.append((i, j))

                    # next energy calculations
                    next_e = self.fc.eval_move(
                        current_string, i, j) + current_e
                    next_e = round(next_e, 2)

                    next_p_table = current_p_table.copy()
                    if i < 0:
                        next_p_table[-i] = 0
                        next_p_table[-j] = 0
                    else:
                        next_p_table[i] = j
                        next_p_table[j] = i

                    # next saddle energy
                    next_s = round(max(current_s, next_e), 2)

                    # if this move is valid... append to list
                    if next_s <= max_energy:

                        next_moves = current_moves.copy()
                        next_moves.append((i, j, next_e))

                        next_energies = current_energies.copy()
                        if next_e < s1_en:
                            next_energies.append(round(next_e-s1_en, 2))
                        else:
                            next_energies.append(0)

                        # unused?
                        en_moves = [x[2] for x in next_moves]

                        new_intermediate = Intermediate(p_table=next_p_table, mode=mode, saddle_e=next_s, current_e=next_e,
                                                        moves=next_moves, energies=next_energies, opt=[], add_moves=current_add_moves)

                        new_path = current_path.copy() + [new_intermediate]
                        collect_paths_per_move.append(new_path)
                        # collect_paths.append(new_path)



                collect_paths_per_move.sort(key=lambda x: x[-1].saddle_e)

                # if current_bp == 3:
                #     collect_paths_per_move = collect_paths_per_move[0:2]
                # if current_bp >=4 and current_bp <=6:
                #     collect_paths_per_move = collect_paths_per_move[0:2]
                # if current_bp >=12:
                #     collect_paths_per_move = collect_paths_per_move[0:3]

                # print (current_bp, add_moves, delete_moves)

                # if add_moves>delete_moves:
                #     collect_paths_per_move = collect_paths_per_move[0:3]

                # collect_paths_per_move = collect_paths_per_move[0:5]

                collect_paths += collect_paths_per_move



            # first sorting step
            collect_paths.sort(key=lambda x: (x[-1].p_table, x[-1].saddle_e))

            last_ptable = []
            # last_ptable = collect_paths[-1][0].p_table
            # last_ptable = collect_paths[0][-1].p_table
            print_d("sort done", last_ptable, init_intermediate.p_table)

            # remove duplicates ptables
            if current_bp+1 != current_bp_end:
                for i in range(len(collect_paths)):
                    c_ptable = collect_paths[i][-1].p_table
                    if c_ptable == last_ptable:
                        # set saddle energy high
                        collect_paths[i][-1].saddle_e = 999
                    else:
                        last_ptable = c_ptable
                i = 0
                while True:
                    if i > len(collect_paths)-1:
                        break
                    if collect_paths[i][-1].saddle_e == 999:
                        collect_paths.pop(i)
                        continue
                    i += 1

            # second sorting step
            collect_paths.sort(key=lambda x: (x[-1].saddle_e, x[-1].current_e))

            print("iteration", current_bp, len(collect_paths))

            # discard paths
            # if current_bp == 6:                 
            #     width = 1

            info = [0.125, 0.1875, 0.0, 0.058823529411764705, 0.0, 0.0, 0.0, 0.07692307692307693, 0.08333333333333333, 0.16666666666666666, 0.0, 0.16666666666666666, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            
            # if info[current_bp] == 0:
            #     width = 1


            collect_paths = collect_paths[:width]

            # return valid paths if we're past the minimum bp_dist
            if current_bp >= current_bp_end-1:# or current_bp==30:
                for i in range(len(collect_paths)):
                    if collect_paths[i][-1].p_table == list(end_p_table):
                        # print("y")
                        yield collect_paths[i][-1]
                    # yield collect_paths[i][-1]


            # next iteration
            paths = collect_paths
            current_bp += 1

        # return remaining paths
        if paths:
            for path in paths:
                if path:
                    yield path[-1]


def find_path(sequence, s1, s2, add_moves=[], results=1, indirect_iterations=2, search_width=1000, Debug=False, Verbose=False):
    """
    indirect findpath, main function

    settings:
    indirect_iterations, default value 2 (means 1 direct pass, 1 indirect pass)

    """

    fp_class = Fp_class(sequence, s1, s2, add_moves)
    max_energy = float("inf")
   

    accelerator = 1.0;  # unused
    last_iteration = search_width # search width for final pass
    iterations = [search_width]

    # iterations = [4, 20, 100, 500]
    # the first iterations should be somewhere between 2 and 16
    while (last_iteration > 16):
        next_iteration = int(last_iteration / 5.0 * accelerator)        
        iterations.insert(0, next_iteration)
        last_iteration = next_iteration

    saddle_en1 = max_energy
    saddle_en2 = max_energy

    paths = []

    iterations = [search_width]

    for search_width in iterations:

        # fwd path
        for path in fp_class.find_path_once(
                s1, s2, max_energy, search_width, mode=True, sort_min=False, Debug=Debug, Verbose=Verbose):

                if path.saddle_e < saddle_en1:
                    saddle_en1 = path.saddle_e
                current_moves = []
                for pos, (i, j, e) in enumerate(path.moves):
                    current_moves.append((i, j))                
                paths.append([path.saddle_e, search_width, True, current_moves])

        # bwd path
        # for path in fp_class.find_path_once(
        #         s2, s1, max_energy, search_width, mode=False, sort_min=False, Debug=Debug, Verbose=Verbose):
        #         if path.saddle_e < saddle_en1:
        #             saddle_en1 = path.saddle_e
        #         current_moves = []
        #         # flip backwards path
        #         for pos, (i, j, e) in enumerate(path.moves):
        #             if i == 0: 
        #                 continue
        #             current_moves.insert(0, (-i, -j))                  
        #         current_moves.insert(0, (0, 0))  
        #         paths.append([path.saddle_e, search_width, False, current_moves])

        if max_energy > saddle_en1:
            max_energy = saddle_en1
        if max_energy > saddle_en2:
            max_energy = saddle_en2               

    # best path to 0
    paths.sort(key = lambda x: x[0])

    return paths


if __name__ == '__main__':

    # various random sequences

    # # 60.1 inner
    # sequence = "GCCAACAAACCGUGAUGGGCUAUGUUUGAUUCAUCCUAUUUAUGUUUUUCGAGAUGCGCG"
    # s1       = "...........(((..(((..((.....))...)))....)))................."
    # s2       = "...........((((((((..(((.......)))))))).)))................."

    sequence   = "CCAGCGUAUUAGUUAUGGCCUGGAGGUAGAAGCGUUAGAGCAAUACUUCUACAGAGACCACGUGAGGUAG"
    s1         = "((((..((.......))..))))..(((((((.............))))))).................."
    s2         = "((((((((.....))).)).)))..(((((((.((.......)).)))))))....(((......))).."

    sequence = 'AGGUGUUAUAAAAGUGCAGAGCCAAGUCACCGCCUUGCCGUCGCUUCCGCUAUGCCUCAAAAGCGGAAUAUCGAGCGAGCGAUUGUGAUUCGAAUGGUCG'
    s1       = '.(((................))).((((((....(((((((((.(((((((.(......).)))))))...))).)).))))..))))))..........'
    s2       = '...(((.........)))..((((((((((.(((((((((....(((((((.(......).)))))))...)).))))).).).))))))....))))..'


    section = ()
    search_width = 60
    Verbose = True
    # Debug = True
    Debug = False


    add_moves = []
    indirect_iterations = 1
    paths = find_path(sequence, s1, s2, indirect_iterations=indirect_iterations, add_moves=add_moves,
                      search_width=search_width, Debug=Debug, Verbose=Verbose)

    for s, sw, mode, moves in paths:
        print_moves(sequence, s1, s2, moves, Verbose=Verbose)
        break

    print ([(i[0], i[1]) for i in paths])
