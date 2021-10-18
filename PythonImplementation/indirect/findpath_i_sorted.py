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
import pathfinder


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
    distance:     int


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
        current_bp_end = bp_dist(s1, s2)     
        
        s1_en = round(evals[s1], 2)
        end_p_table = self.p_tables[s2]

        init_intermediate = Intermediate(p_table=list(p_tables[s1]), mode=mode, saddle_e=s1_en, current_e=evals[s1], moves=[(0, 0, s1_en)], energies=[0], opt=0, add_moves=[], distance=current_bp_end)
        
        # initial path start with 1 intermediate
        init_path = [init_intermediate]
        # paths start with 1 initial path
        paths = [init_path]


        # Ivos variation: define the largest distance, 
        # next iteration will only ever iterate over the longest distances

        longest_distance = current_bp_end


        # dont stop at current_bp_end, consider potential indirect moves
        while (len(paths)!=0):

            # copy all paths from the queue into 'paths' which have the longest distance


            # iterate over all queue items for longest distance first

            longest_distance = 0

            for p in paths:
                if p[-1].distance > longest_distance:
                    longest_distance = p[-1].distance

            # all elements in paths, which dont have the longest distance, get onto the queue.
            # from the queue, all items with the longest distance get to paths


            # print ('longest d', longest_distance, 'paths:', len(paths), 'queue:', len(queue),            
            #  'avg', np.mean(test), np.max(test), np.min(test))
         


            next_longest_distance = 0
            
            # collect all new paths here (next iteration)
            collect_paths = []
            ignored = []

            for current_path in paths:

                # print ('collect', current_path[-1].distance, longest_distance)
                if current_path[-1].distance != longest_distance:
                    ignored.append(current_path.copy())
                    continue

                current_distance = current_path[-1].distance
                # if current_distance < longest_distance:
                #     print ('ignore', current_distance)
                # else:
                #     print ("accept", current_distance)

                current_p_table = current_path[-1].p_table
               
                current_s = current_path[-1].saddle_e
                current_e = current_path[-1].current_e
                current_string = p_to_s(current_p_table)
                current_moves = current_path[-1].moves
                current_energies = current_path[-1].energies

                # "try_moves"
                for i, j in self.find_moves(current_p_table, end_p_table):

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

                        distance_s2 = RNA.bp_distance(RNA.db_from_ptable(next_p_table), s2)

                        if distance_s2 > next_longest_distance:
                            next_longest_distance = distance_s2

                        # print (distance_s2, end=" ")

                        new_intermediate = Intermediate(p_table=next_p_table, mode=mode, saddle_e=next_s, current_e=next_e,
                                                        moves=next_moves, energies=next_energies, opt=[],
                                                        add_moves=current_add_moves, distance=distance_s2)

                        new_path = current_path.copy() + [new_intermediate]
                        collect_paths.append(new_path)

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

            # discard paths
            collect_paths = collect_paths[:width]

            # return valid paths if we're past the minimum bp_dist
            if current_bp >= current_bp_end-1:
                for i in range(len(collect_paths)):
                    if collect_paths[i][-1].p_table == list(end_p_table):
                        yield collect_paths[i][-1]

            # next iteration

            # for the next iteration, copy only those intermediates into paths, which have 
            # the longest distance. the rest gets appended to a queue.

            paths = []
            for path in collect_paths:
                paths.append(path)

            # add postponed moves back into the paths queue
            for path in ignored:
                # print("add back", path[-1].distance)
                paths.append(path)
                               

            # paths = collect_paths
            current_bp += 1
            longest_distance = next_longest_distance
            # end while loop, get to the next distance class

        # return remaining paths
        if paths:
            for path in paths:
                if path:
                    yield path[-1]


def find_path(sequence, s1, s2, add_moves=[], search_width=1000, Debug=False, Verbose=False):
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
        for path in fp_class.find_path_once(
                s2, s1, max_energy, search_width, mode=False, sort_min=False, Debug=Debug, Verbose=Verbose):

                if path.saddle_e < saddle_en1:
                    saddle_en1 = path.saddle_e
                current_moves = []

                # flip backwards path
                for pos, (i, j, e) in enumerate(path.moves):
                    if i == 0: 
                        continue
                    current_moves.insert(0, (-i, -j))                  
                current_moves.insert(0, (0, 0))  
                paths.append([path.saddle_e, search_width, False, current_moves])

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

    # sequence = "CACGCUGGGGAAUUCCGACGAUAAAUAACCCUUCGCGCCCAUAGGAACGUACACAAAUCUGAAGCUUACG"
    # s1       = "..(((.(((...................)))...)))((....))..((((...............))))"
    # s2       = "..(((.((((..................))))..)))................(...............)"

    # sequence = "GCCGCCUGAGCCAUUUACAAACUAUCGGAUACAGUGUGUCUUAGAUGUUGUGGCAACAUUAGAUAAGGUA"
    # s1       = "(((..(((((((((..(((..(((..((((((...))))))))).))).)))))....))))....)))."
    # s2       = "...((((..(((((..(((..(((...(((.......))).))).))).)))))...........))))."

    # sequence   = "UAGGGGUGCACUAAAGCUGGUAUCCCCUAUGAGUGGAUAAAUGAUACAGGUCACCCUACGACAUAUACGC"
    # s1         = "(((((.((.(((........((((((......).))))).........))))))))))............"
    # s2         = "..((((((..((........(((((.........)))))........))..)))))).((.......))."

    sequence   = "CCAGCGUAUUAGUUAUGGCCUGGAGGUAGAAGCGUUAGAGCAAUACUUCUACAGAGACCACGUGAGGUAG"
    s1         = "((((..((.......))..))))..(((((((.............))))))).................."
    s2         = "((((((((.....))).)).)))..(((((((.((.......)).)))))))....(((......))).."

    # sequence   = "UCUUGAACCCAUGGCGUUUUCAACCGACAUCCUGCUCCCGCAAUCACCUAGGUUAAGGGUUCUUCAAGGA"
    # s1         = "((((((((((.((((.........(.......(((....)))........))))).))))...))))))."
    # s2         = "....((((((.(((((((...)))........(((....))).........)))).))))))........"

    # sequence = "UCACUGAGGCUUGUUCGCAAAUCACUGCAAUUAGAUAUGACUCACGAUAUGGGGCACGGUGCAUACAUAC"
    # s1       = ".(((((.(.((((((((....((...............))....))).))))).).)))))........."
    # s2       = ".(((((..((((..(((..(.(((.((........))))).)..)))....)))).)))))........."

    # ~~~

    # Tabu paper example
    # sequence = "CGCGACGGCUACGCGACGGCAAUGCCGUUGCGAAGCCGUCGCGAUC"
    # s1 = "(((((((((..............))))))))).............."
    # s2 = "...........(((((((((..............)))))))))..."

    # rna2dfold example
    # sequence = "GGGCGCGGUUCGCCCUCCGCUAAAUGCGGAAGAUAAAUUGUGUCU"
    # s1 = "(((((.....)))))(((((.....)))))(((((.....)))))"
    # s2 = "((((((((((.....(((((.....))))).....))))))))))"

    # EApath supplement

    # rb2
    # sequence = "GGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCGUAAACUCUUGAUUAUGAAGUCUGUCGCUUUAUCCGAAAUUUUAUAAAGAGAAGACUCAUGAAU"
    # s1       = "............((((((.........))))))........((((((.......)))))).((((((.((((.((.((((((..........)))))).)))))))))))).."
    # s2       = "(((((((((...((((((.........))))))........((((((.......))))))..))))))))).........................................."

    # rb3
    # sequence = "CUCUUAUCAAGAGAGGUGGAGGGACUGGCCCGAUGAAACCCGGCAACCAGCCUUAGGGCAUGGUGCCAAUUCCUGCAGCGGUUUCGCGUUGAAAGAUGAGAGAUUCUUGUAGUCUCUUCUUUUAGCGAAGGGACUUUUUUU"
    # s1       = "((((((((....(.(((...(((.....)))......))))(((..((((((....))).))).)))........(((((......)))))...)))))))).........((((((((.......))))))))......."
    # s2       = "............(.(((...(((.....)))......))))(((..((((((....))).))).)))........(((((......)))))(((((.((((...........)))).)))))..................."

    # rb5 -26.40 kcal/mol vs -23.60 kcal/mol direct
    # sequence = "GGGAAUAUAAUAGGAACACUCAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGACUAUGGGUGAGCAAUGGAACCGCACGUGUACGGUUUUUUGUGAUAUCAGCAUUGCUUGCUCUUUAUUUGAGCGGGCAAUGCUUUUUUUAUUCUCAUAACGGAGGUAGACAGGAUGGAUCCACUGA"
    # s1       = "................((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))........(((((........)))))............((((((((((((((.......))))))))))))))..........................................."
    # s2       = ".....................(((...(.(((((.......))))).)........((((((.......))))))..)))((((((((((((.(((((........)))))..............))))))))))))................................................................"

    # dsrA
    # sequence = "ACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUAAGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUU"
    # s1 = "..(((((((.....)))))))...(((((((((((((.......))))))).)))))).((((((((((.....))))))))))."
    # s2 = "..(((((((.....)))))))................(...(((((....)))))...)((((((((((.....))))))))))."

    # HIV leader
    # sequence = "GGUCUCUCUGGUUAGACCAGAUUUGAGCCUGGGAGCUCUCUGGCUAACUGGGAACCCACUGCUUAAGCCUCAAUAAAGCUUGCCUUGAGUGCUUCAAGUAGUGUGUGCCCGUCUGUUGUGUGACUCUGGUAACUAGAGAUCCCUCAGACCCUUUUAGUCAGUGUGGAAAAUCUCUAGCAGUGGCGCCCGAACAGGGACUUGAAAGCGAAAGGGAAACCAGAGGAGCUCUCUCGACGCAGGACUCGGCUUGCUGAAGCGCGCACGGCAAGAGGCGAGGGGA"
    # s1       = "((..((.(..(((((.(((((...((((......))))))))))))))..)))..)).((.(((..(((((........(((((.((.((((((((.((((.....(((.((((((((.((..(((((((.....(((....)))...(((((((.((((.((((((.....)))).)).)))).(((.....)))..........)))))))..)))))))..)).....))))...))))..))))))))))))))).)).))))))))))))).))."
    # s2       = "((..((.(..(((((.(((((...((((......))))))))))))))..)))..)).((.(((..(((((........(((((.((.(((((((.(((((.....(.(((...(((.((((.(((((((..((.(((....))))).(((((((.((((.((((((.....)))).)).))))((((.....))).)........)))))))..)))))))(((....)))..)))).))).)))))))))))))))).)).))))))))))))).))."

    # barriers
    # ((((((((((...(.(((((.....))))).)...)))))))))) (-17.70) L0002
    # ((((((((((.....(((((.....))))).....)))))))))) (-17.70) I
    # (((((((((......(((((.....)))))......))))))))) (-16.60) I
    # ((((((((.......(((((.....))))).......)))))))) (-15.00) I
    # (((((((........(((((.....)))))........))))))) (-14.40) I
    # (((((((...)....(((((.....))))).........)))))) ( -8.90) I
    # (((((((...))...(((((.....)))))..........))))) ( -7.00) S
    # (((((((...)))..(((((.....)))))...........)))) (-10.40) L0020
    # .((((((...)))..(((((.....)))))...........))). ( -9.80) I
    # ((.((((...)))..(((((.....)))))...........))). ( -6.20) S
    # ((..(((...)))..(((((.....)))))............)). ( -6.30) L0120
    # (...(((...)))..(((((.....))))).............). ( -6.00) S
    # (...(((...))).)(((((.....)))))............... ( -7.70) I
    # ((..(((...)))))(((((.....)))))............... (-10.50) L0019
    # (((..((...)))))(((((.....)))))............... (-10.40) S
    # ((((..(...)))))(((((.....)))))............... (-11.40) I
    # (((((.....)))))(((((.....)))))............... (-16.80) L0003
    # (((((.....)))))(((((.....))))).(...........). (-12.70) S
    # (((((.....)))))(((((.....))))).((.........)). (-13.90) I
    # (((((.....)))))(((((.....))))).(((.......))). (-15.90) I
    # (((((.....)))))(((((.....))))).((((.....)))). (-17.20) I
    # (((((.....)))))(((((.....)))))(((((.....))))) (-18.10) L0001

    # s1, s2 = s2, s1

    section = ()
    search_width = 50
    Verbose = True
    # Debug = True
    Debug = False

    # add_moves = [(39, 45), (38, 46)]
    # indirect_iterations = 1

    # tabu paper example
    # add_moves = [(1, 43)]
    # add_moves = [(1, 43), (2, 42)]
    # add_moves = [(1, 43), (2, 42), (3, 41)]
    # add_moves = [(1, 43), (2, 42), (3, 41), (4, 40)]
    # indirect_iterations = 1

    add_moves = [(9, 15), (36, 42)]

    # add_moves = []
    # indirect_iterations = 2
    indirect_iterations = 1
    paths = find_path(sequence, s1, s2, indirect_iterations=indirect_iterations, add_moves=add_moves,
                      search_width=search_width, Debug=Debug, Verbose=Verbose)


    for s, sw, mode, moves in paths:
        print_moves(sequence, s1, s2, moves, Verbose=Verbose)
        break

    print ([(i[0], i[1]) for i in paths])



    # print (se)
    print('orig findpath:')

    pathfinder_result = pathfinder.pathfinder(
        sequence, s1, s2, search_width=search_width, verbose=Verbose)

    print(pathfinder_result.max_en)
    # print (pathfinder.pathfinder(sequence, s2, s1, section=section, search_width=search_width, verbose=Verbose).sE)


def detect_local_minimum(fc, structure):
    # perform gradient walk from sample to determine direct local minimum
    pt = RNA.IntVector(RNA.ptable(structure))
    fc.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
    return RNA.db_from_ptable(list(pt))


def get_neighbors(fc, db=None, pt=None):
    """
    """
    if pt is None:
        pt = RNA.ptable(db)
    else:
        assert db == None
    nbrs = []
    for move in fc.neighbors(pt):
        npt = list(pt)
        if move.is_removal():
            npt[-move.pos_3] = 0
            npt[-move.pos_5] = 0
            ndb = RNA.db_from_ptable(npt)
        elif move.is_insertion():
            npt[move.pos_3] = move.pos_5
            npt[move.pos_5] = move.pos_3
        else:
            rlog.warning(f"Are you using shift moves?")
            rlog.warning(f" shift = {move.is_shift()}")
            rlog.warning(f" pos3 = {move.pos_3}")
            rlog.warning(f" pos5 = {move.pos_5}")
            raise NotImplementedError('Are you using shift moves?')
        dG = fc.eval_move_pt(pt, move.pos_5, move.pos_3)
        if db:
            nbrs.append([RNA.db_from_ptable(npt), dG])
        else:
            nbrs.append([npt, dG])
        return nbrs
