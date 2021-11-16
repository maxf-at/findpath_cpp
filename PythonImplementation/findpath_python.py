#!/usr/bin/env python3
# coding: utf-8

# generic Python libraries
import numpy as np
import pandas as pd
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

from sklearn.cluster import DBSCAN
import tensorflow as tf
from tensorflow import feature_column
from tensorflow.keras import layers

reloaded_model = tf.keras.models.load_model('dnn_model')

# import feature_generation
from features import ij_distance, new_move_dist, plt_moves, config_distance, balance_in_all_things, return_shift_moves
from process_features import fp_call, find_moves, process



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








def next_vicinity(lastmove, avail_moves):

    """
    percentage of moves which can continue in the vicinity of the last move
    """
    lasti = lastmove[0]
    lastj = lastmove[1]
    moves = len(avail_moves)
    cntr = 0

    for a in avail_moves:
        i=a[0]
        j=a[1]      
        if lasti + 1 == i or lasti - 1 == i or\
           lastj + 1 == j or lastj - 1 == j:
           cntr += 1
    return cntr


def last_vs_best(lastmove, avail_moves, lensequence):
    lasti = abs(lastmove[0])
    lastj = abs(lastmove[1])
    i = abs(avail_moves[0][0])
    j = abs(avail_moves[0][1])

    dist = abs(lasti - i) + abs(lastj - j)
    return dist/(lensequence*2) #, dist, lensequence


def cluster_moves(moves):
    moves = np.array([(abs(i[0]), abs(i[1])) for i in moves])
    clustering = DBSCAN(eps=2, min_samples=1).fit(moves)
    unique_labels = len(set(clustering.labels_))
    # normalize to 1
    return 1/unique_labels


# def structure_evaluation(sequence, s, s2):
def structure_evaluation(fc, pt1, pt2, lastmove):
    avail_moves = []
    fp_results = []
    en_contrib = []

    available_add = set()
    available_delete = set()
    lensequence = pt1[0]

    # pt1 = list(RNA.ptable(s))
    # pt2 = list(RNA.ptable(s2))
    for pos, (i,j) in enumerate(find_moves(pt1, pt2)):    
        next_en = fc.eval_move_pt(pt1, i, j)
        # mark where we found our move    
        # map energies somehow to [0,1] such that the network will understand.  
        en = np.interp(next_en/100, [-10,10], [0,1]) 
        en_contrib.append(en)
        avail_moves.append((i, j, next_en))

        if i>0:
            available_add.add((i,j))
        elif i<0:
            available_delete.add((i,j))

    if len(avail_moves)==0:
        add_delete = 0
    else:
        add_delete = len(available_delete)/len(avail_moves)

    # print ('best', en_contrib[0], 'avg', en_mean)

    avail_moves.sort(key=lambda x: x[2])


    vic = next_vicinity(lastmove, avail_moves)
    vic_best = next_vicinity(lastmove, avail_moves[0:1])

    distlast = last_vs_best(lastmove, avail_moves, lensequence)

    unique_moves = cluster_moves(avail_moves)

    # print ("vicinity", vic, lastmove, "->", avail_moves)

    # print ('avail add', available_add)
    # print ('avail del', available_delete)

    if len(avail_moves) == 0:
        en_mean = 0
        en_std = 0
        best_en = 0


    else:
        en_mean = np.mean(en_contrib)
        en_std = np.std(en_contrib)
        en_contrib.sort(key=lambda x: x)
        en_diff = en_contrib[0] / en_mean # best en, relative

        best_en = en_contrib[0]


    # compare with prediction
    # data = {'en_mean': [en_mean],
    #         'en_std': en_std,
    #         'best_en': best_en,
    #         'vic': vic,
    #         'vic_best': vic_best,
    #         'unique_moves': unique_moves,
    #         'distlast': distlast}
    # test_features = pd.DataFrame.from_dict(data)  

    return (en_mean, en_std, best_en, vic, vic_best, unique_moves, distlast)

    # y = reloaded_model.predict(test_features)[0][0]
    y = 1
    
    distlast = round(y,4)

    
    return distlast











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

    def find_path_once(self, s1, s2, max_energy, width, mode=True, sort_min=False,\
        Debug=False, Verbose=False, limit=False):
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

        eval_counter = 0

        # dont stop at current_bp_end, consider potential indirect moves
        while (current_bp != current_bp_end+2*len(self.moves_add)):

            # collect all new paths here (next iteration)
            collect_paths = []
            collect_paths2 = []
            collect_seval = []

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

                avail_moves = []

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

                    avail_moves.append((i,j, self.fc.eval_move(current_string, i, j)))
                    eval_counter += 1

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

                avail_moves.sort(key=lambda x: x[2])

                lastmove = collect_paths_per_move[0][-1].moves[-2]
                selected = collect_paths_per_move[0][-1].moves[-1]

                seval = structure_evaluation(self.fc, current_p_table, end_p_table, lastmove)
                collect_seval.append(seval)

                # if seval < 0.05:
                #     collect_paths_per_move = collect_paths_per_move[0:1]

                collect_paths2.append(collect_paths_per_move)

            collect_seval = pd.DataFrame(collect_seval, columns=['en_mean', 'en_std',\
                 'best_en', 'vic', 'vic_best', 'unique_moves', 'distlast'])
            
            if limit:
                y = reloaded_model.predict(collect_seval).T[0]
            
            # print (collect_seval)
            # print (y)

            if (current_bp == 19 or current_bp == 20) and limit:
                print_d(y)

            for i in range(len(collect_seval)):
                
                
                # if limit and y[i] < 0.05:
                # if limit and y[i] < 0.05 and (current_bp < 19 or current_bp > 20):
                #     collect_paths2[i] = collect_paths2[i][0:1]
                
                # if limit and y[i] < 0.05 and y[i] > 0.00040:
                #     collect_paths2[i] = collect_paths2[i][0:1]
                
                collect_paths += collect_paths2[i]


            # first sorting step
            collect_paths.sort(key=lambda x: (x[-1].p_table, x[-1].saddle_e))


            last_ptable = []
            # print_d("sort done", last_ptable, init_intermediate.p_table)

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

            print_d("iteration", current_bp, len(collect_paths), "evals:", eval_counter)

            # discard paths
            # if current_bp == 6:                 
            #     width = 1

            


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


def find_path(sequence, s1, s2, add_moves=[], results=1, indirect_iterations=2,\
    search_width=1000, Debug=False, Verbose=False, limit=False):
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
                s1, s2, max_energy, search_width, mode=True, sort_min=False,\
                Debug=Debug, Verbose=Verbose, limit=limit):

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
    search_width = 30

    sequence = "CCCGUUUCAGCGUCUGGGUCGUAAUCGUCUUGUCUCUGGGCGAAUCGAUAACUGAUCCGUCACUACAGCCGAAGCCAAAGCCAAUGACGGUCCCACCUAU"
    s1 = "...(((((.((((.((((((((.((((..(((((....)))))..)))).)).))))))..))....)).)))))....(((......)))........."
    s2 = "..(((....)))..(((((.(..((((((......(((((((.((((.....)))).))))....))).....((....))....))))))..))))))."
    search_width = 50

    sequence = 'AGGUGUUAUAAAAGUGCAGAGCCAAGUCACCGCCUUGCCGUCGCUUCCGCUAUGCCUCAAAAGCGGAAUAUCGAGCGAGCGAUUGUGAUUCGAAUGGUCG'
    s1       = '.(((................))).((((((....(((((((((.(((((((.(......).)))))))...))).)).))))..))))))..........'
    s2       = '...(((.........)))..((((((((((.(((((((((....(((((((.(......).)))))))...)).))))).).).))))))....))))..'
    # search_width = 30

    sequence = "GACCGACUUACUACCCGCAGACGCAAAGUAAUCAAGACACCGUAUACCAAUGCCGUCCUAUAUUGGGAAGCCAGUUAAAGCGUACUAUUUCGGACCAAGA"
    # s1       = "..((((..........((....))..............((.((((....)))).)).....(((((....)))))..............))))......."
    s1       = "..((((.(((((....((....))..))))).......((.((((....)))).))..((.(((((....)))))))............))))......."
    s2       = "..((((..((.(((..((.((((((..(((..(........)..)))...)).))))....(((((....)))))....))))).))..))))......."


    section = ()
    search_width = 100
    Verbose = True
    Debug = True
    # Debug = False

    limit=False
    limit=0.003

    # paths = find_path(sequence, s1, s2, indirect_iterations=1, add_moves=[],
    #                   search_width=search_width, Debug=Debug, Verbose=Verbose, limit=limit)
    # for s, sw, mode, moves in paths:
    #     print_moves(sequence, s1, s2, moves, Verbose=Verbose)
    #     break

    # print ([(i[0], i[1]) for i in paths])

    filename_samples = f'./dataset_100_large.csv'
    samples_df = pd.read_csv(filename_samples)
    for index, row in samples_df.iterrows():
        # if index != 809:
        #     continue
        if index != 13:
            continue
        # if index == 20:
        #     break
        sequence = row.sequence
        s1 = row.s1
        s2 = row.s2
        en1, en2 = 0, 0

        output_path = False
        output_path = True

        # limit=False
        # paths = find_path(sequence, s1, s2, indirect_iterations=1, add_moves=[],
        #                 search_width=search_width, Debug=Debug, Verbose=Verbose, limit=limit)
        # for s, sw, mode, moves in paths:
        #     en1 = print_moves(sequence, s1, s2, moves, Verbose=output_path)
        #     break
        en1 = 0

        limit=0.05
        paths = find_path(sequence, s1, s2, indirect_iterations=1, add_moves=[],
                        search_width=search_width, Debug=Debug, Verbose=Verbose, limit=limit)
        for s, sw, mode, moves in paths:
            en2 = print_moves(sequence, s1, s2, moves, Verbose=output_path)
            break

        print(f"sequence = \"{sequence}\"")
        print(f"s1 = \"{s1}\"")
        print(f"s2 = \"{s2}\"")
        print (RNA.bp_distance(s1, s2))

        print ('regular', en1, 'with limit', en2)
