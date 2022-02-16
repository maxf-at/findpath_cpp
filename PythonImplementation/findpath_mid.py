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

sys.path.append('../')
# from pretty_print_path import print_moves
import findpath_librna
import findpath
import findpath_helix
# import helper

from sklearn.cluster import DBSCAN
from scipy.spatial.distance import cdist
# import tensorflow as tf
# from tensorflow import feature_column
# from tensorflow.keras import layers

# reloaded_model = tf.keras.models.load_model('dnn_model4')

# # import feature_generation
# from features import ij_distance, new_move_dist, plt_moves, config_distance, balance_in_all_things, return_shift_moves
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
    cluster_counter:    list
    lastaddcluster: int
    lastaddleft: int
    lastdelcluster: int
    lastdelleft: int

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

    # return (en_mean, en_std, best_en, vic, vic_best, unique_moves, distlast)
    return (en_mean, en_std, best_en)


    # y = reloaded_model.predict(test_features)[0][0]
    y = 1
    
    distlast = round(y,4)

    
    return distlast



def cluster_moves(inputmoves, absolute_values=False):
    """
    cluster moves with DBscan. needs to include the last move as well...
    """

    if absolute_values:
        m = [(abs(i[0]), abs(i[1])) for i in inputmoves]
        moves = np.array(m)
    else:
        m = [(i[0], i[1]) for i in inputmoves]
        moves = np.array(m)
    # print ("clustering", moves)
    clustering = DBSCAN(eps=2, min_samples=1).fit(moves)    
    # generate clustering dict (move-tuple:label, ...)
    return {i[0]:i[1] for i in zip(m, clustering.labels_)}, clustering.labels_

from process_features import fp_call, find_moves, process

def adjust_pt(pt, i, j):
    pt_adj = pt.copy()
    if i<0:
        pt_adj[-i] = 0
        pt_adj[-j] = 0
    else:
        pt_adj[i] = j
        pt_adj[j] = i
    return pt_adj


def cluster_eval(i,j, clustered_moves, avail_clusters): # +global vars...
    

    current_label = clustered_moves[(i,j)]
    members = avail_clusters[current_label]

    clusterpos = members.index((i,j))
    adj_clusterpos = members.index((i,j))
    if len(members)-1 != 0:
        adj_clusterpos /= (len(members)-1)
    adj_clusterpos = abs(abs(adj_clusterpos-0.5)-0.5)*2

    if i<0:
        # del move. helix removal from the side
        # print ("adj", adj_clusterpos)
        return adj_clusterpos
    
    return 0







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
        Debug=False, Verbose=False, limit=False, hopp=False):
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
        midpoint = int(current_bp_end/2)


        # clustering
        # if limit:
        fp = findpath.findpath_single(sequence, s1, s2, search_width_multiplier=1, mp=True)
        p = fp.get_path()      
        # p = [(0, 0)] + [(i[0], i[1]) for i in p]
        p = [(i[0], i[1]) for i in p]
        p.sort(key=lambda x: (x[0], x[1]))
        clustered_moves, clabels = cluster_moves(p, absolute_values=False)


        init_clusters = defaultdict(list)
        add_clusters = defaultdict(list)
        del_clusters = dict()

        for i,j in p:
            if i==0: continue
            current_label = clustered_moves[(i,j)]
            if i<0:
                init_clusters[current_label].append((i,j))
            else:
                add_clusters[current_label].append((i,j))

        for a,b in init_clusters.items():
            del_clusters[a] = len(b)

            # print_d (i, j, current_label)

        # print (clustered_moves)
        
        # init_clusters = [0] * (clabels[-1]+1)
        # print (init_clusters)
        # print (add_clusters)

        lastaddcluster = -1
        lastaddleft = -1
        lastdelcluster = -1
        lastdelleft = -1

        s1_en = round(evals[s1], 2)
        end_p_table = self.p_tables[s2]
        init_intermediate = Intermediate(p_table=list(p_tables[s1]), mode=mode, saddle_e=float(
            "-inf"), current_e=evals[s1], moves=[(0, 0, s1_en)], energies=[0], opt=0, cluster_counter=init_clusters,
                                                        lastaddcluster=lastaddcluster, lastaddleft=lastaddleft, 
                                                        lastdelcluster=lastdelcluster, lastdelleft=lastdelleft                                                        
                                                        )
        # initial path start with 1 intermediate
        init_path = [init_intermediate]
        # paths start with 1 initial path
        paths = [init_path]

        eval_counter = 0
        sorting_couter = 0





        # dont stop at current_bp_end, consider potential indirect moves
        while (current_bp != current_bp_end+2*len(self.moves_add)):

            # collect all new paths here (next iteration)
            collect_paths = []
            collect_paths2 = []
            collect_seval = []
            collect_moves = []

            for current_path in paths:

                collect_paths_per_move = []
                delete_moves = 0

                current_p_table = current_path[-1].p_table
                current_e = current_path[-1].current_e
                current_s = current_path[-1].saddle_e
                current_string = p_to_s(current_p_table)
                current_moves = current_path[-1].moves
                current_energies = current_path[-1].energies

                current_clusters = current_path[-1].cluster_counter

                avail_moves = []


                # simple helix model: best move per cluster (greedy step IN a cluster)
                best_per_cluster = defaultdict(lambda: np.inf)
                best_per_cluster_id = dict()
                for i, j in self.find_moves(current_p_table, end_p_table):
                    # if i<0:
                    #     continue

                    next_cluster = clustered_moves[(i,j)]
                    e = self.fc.eval_move(current_string, i, j)
                    if best_per_cluster[next_cluster]>e:
                        best_per_cluster[next_cluster] = e
                        best_per_cluster_id[next_cluster] = i


                # cluster diversity

                # compare init clusters with current c.


                # diversity = dict()

                # print (current_clusters)


                # print (best_per_cluster_id)

                # "try_moves"
                for i, j in self.find_moves(current_p_table, end_p_table):

                    if i<0:
                        delete_moves += 1


                    # next energy calculations
                    next_e = self.fc.eval_move(
                        current_string, i, j) + current_e
                    next_e = round(next_e, 2)

                    next_cluster = clustered_moves[(i,j)]
                    next_clusters = current_clusters.copy()
                    

                    # print (i,j, next_cluster)
                    next_clusters[next_cluster] = current_clusters[next_cluster].copy() #deepcopy
                    next_clusters[next_cluster].sort(key=lambda x: x[0])




                    diversity = np.zeros(shape=(len(add_clusters)+len(del_clusters)+1))
                    for a,b in next_clusters.items():
                        # if a not in del_clusters:
                        #     continue
                        diversity[a] = len(b)#/del_clusters[a]
                        # if diversity[a] != 1:
                            # diversity[a] = 0
                        if a in del_clusters:
                            diversity[a]/=del_clusters[a]
                        else:
                            diversity[a]/=len(add_clusters[a])

                    # print (current_bp, diversity, len(add_clusters)+len(del_clusters))


                    # continue with current cluster if possible

                    # if hopp and i<0 and lastdelcluster!=-1 and lastdelleft!=0:
                    #     if lastdelcluster != next_cluster: 
                    #         continue

                        # print (current_bp, (i,j), "last", lastaddcluster, "current", next_cluster)

                    # print (next_clusters, a, b)

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
                                                        moves=next_moves, energies=next_energies, opt=diversity, cluster_counter=next_clusters,
                                                        lastaddcluster=False, lastaddleft=False, 
                                                        lastdelcluster=False, lastdelleft=False                                                        
                                                        )

                        new_path = current_path.copy() + [new_intermediate]


                        collect_paths.append(new_path)



           



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

            # print_d("iteration", current_bp, len(collect_paths), "evals:", eval_counter)

            # discard paths
            # if current_bp == 6:                 
            #     width = 1

            


            # if info[current_bp] == 0:
            #     width = 1

            rest = collect_paths[width:]
            collect_paths = collect_paths[:width]

            if limit:
                if current_bp==midpoint:
                    # print ("midpoint", current_bp)

                    for i in range(len(collect_paths)):
                        c_ptable = collect_paths[i][-1].p_table
                        cse = collect_paths[i][-1].saddle_e
                        cs = RNA.db_from_ptable(c_ptable)

                        cmoves = collect_paths[i][-1].moves

                        entest, pathtest, _ = find_path(sequence, s1, cs, indirect_iterations=1, add_moves=[],
                                        search_width=width, Debug=Debug, Verbose=Verbose, limit=False, hopp=hopp)

                        cmoves2 = pathtest[0][3]


                        se1 = print_moves(sequence, s1, cs, cmoves, Verbose=False)                        
                        se2 = print_moves(sequence, s1, cs, cmoves2, Verbose=False)  

                        if se2<se1:
                            # print (i, cse, cs, "se1se2", se1, se2)
                            # print (cmoves, se1)
                            # print (cmoves2, se2)
                            # print_moves(sequence, s1, cs, cmoves, Verbose=True)    
                            # print_moves(sequence, s1, cs, cmoves2, Verbose=True)  

                            # replace path
                            collect_paths[i][-1].saddle_e = se2
                            collect_paths[i][-1].moves = pathtest[0][3]



            # return valid paths if we're past the minimum bp_dist
            if current_bp >= current_bp_end-1:# or current_bp==30:
                for i in range(len(collect_paths)):
                    if collect_paths[i][-1].p_table == list(end_p_table):
                        # print("y")
                        yield collect_paths[i][-1], eval_counter


                    # yield collect_paths[i][-1]


            # next iteration
            paths = collect_paths
            current_bp += 1

        # print ('sorting', limit, sorting_couter)
        # print ('eval counter', eval_counter, limit, sorting_couter, "maxen:", max_energy)
        # return remaining paths
        if paths:
            for path in paths:
                if path:
                    yield path[-1], eval_counter


def find_path(sequence, s1, s2, add_moves=[], results=1, indirect_iterations=2,\
    search_width=1000, Debug=False, Verbose=False, limit=False, simple=False, hopp=False):
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

    if simple: iterations = [search_width]


    evf = 0
    evb = 0

    for search_width in iterations:

        # fwd path
        for path, eval_counter in fp_class.find_path_once(
                s1, s2, max_energy, search_width, mode=True, sort_min=False,\
                Debug=Debug, Verbose=Verbose, limit=limit, hopp=hopp):

                if path.saddle_e < saddle_en1:
                    saddle_en1 = path.saddle_e
                current_moves = []
                for pos, (i, j, e) in enumerate(path.moves):
                    current_moves.append((i, j, e))                
                paths.append([path.saddle_e, search_width, True, current_moves])

                if eval_counter>evf:
                    evf = eval_counter
                # print (path.saddle_e, max_energy)

        # bwd path
        if not simple:
            for path, eval_counter in fp_class.find_path_once(
                    s2, s1, max_energy, search_width, mode=False, sort_min=False,\
                    Debug=Debug, Verbose=Verbose, limit=limit, hopp=hopp):
                    if path.saddle_e < saddle_en1:
                        saddle_en1 = path.saddle_e
                    current_moves = []
                    # flip backwards path
                    for pos, (i, j, e) in enumerate(path.moves):
                        if i == 0: 
                            continue

                        etest = path.moves[pos-1][2]
                        current_moves.insert(0, (-i, -j, etest))                  
                    current_moves.insert(0, (0, 0, e))  
                    paths.append([path.saddle_e, search_width, False, current_moves])

                    if eval_counter>evb:
                        evb = eval_counter

        if max_energy > saddle_en1:
            max_energy = saddle_en1
        if max_energy > saddle_en2:
            max_energy = saddle_en2               
        # print (path.saddle_e, max_energy)
    # best path to 0
    paths.sort(key = lambda x: x[0])

    return max_energy, paths, evf+evb
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

    # sequence = "CCCGUUUCAGCGUCUGGGUCGUAAUCGUCUUGUCUCUGGGCGAAUCGAUAACUGAUCCGUCACUACAGCCGAAGCCAAAGCCAAUGACGGUCCCACCUAU"
    # s1 = "...(((((.((((.((((((((.((((..(((((....)))))..)))).)).))))))..))....)).)))))....(((......)))........."
    # s2 = "..(((....)))..(((((.(..((((((......(((((((.((((.....)))).))))....))).....((....))....))))))..))))))."
    # search_width = 50

    # sequence = 'AGGUGUUAUAAAAGUGCAGAGCCAAGUCACCGCCUUGCCGUCGCUUCCGCUAUGCCUCAAAAGCGGAAUAUCGAGCGAGCGAUUGUGAUUCGAAUGGUCG'
    # s1       = '.(((................))).((((((....(((((((((.(((((((.(......).)))))))...))).)).))))..))))))..........'
    # s2       = '...(((.........)))..((((((((((.(((((((((....(((((((.(......).)))))))...)).))))).).).))))))....))))..'
    # # search_width = 30

    # sequence = "GACCGACUUACUACCCGCAGACGCAAAGUAAUCAAGACACCGUAUACCAAUGCCGUCCUAUAUUGGGAAGCCAGUUAAAGCGUACUAUUUCGGACCAAGA"
    # # s1       = "..((((..........((....))..............((.((((....)))).)).....(((((....)))))..............))))......."
    # s1       = "..((((.(((((....((....))..))))).......((.((((....)))).))..((.(((((....)))))))............))))......."
    # s2       = "..((((..((.(((..((.((((((..(((..(........)..)))...)).))))....(((((....)))))....))))).))..))))......."

    # sequence = "UUCUAGGAACCACCCAACGUCCUAAUAGCUCGAUGGAUGAUCAGGAUGCAUAUCACGCGACGUCCCGAACCUAUACAGACCGGGCAGAACACAGCCACUU"
    # s1 = ".................((((...((((.(((..(((((.((.((((....))).)..))))))))))..))))...)).))(((........)))...."
    # s2 = "...(((((.(........))))))...(((((.((..((.((.((((((.........).))))).))...))..))...)))))..............."
    # s2 = "...(((((.(........))))))..........................................................(((........)))...."


    # sequence = "AUCUGACCUCAGCGGGAGGAGUUUACUGGACUAUAAUAUUGCCCUCCAUGCAACCCAAAAUCUCCAAUGAUUCCAACUAGCGCGAUCAUCACCAAGGAUC"
    # s1 = "...(((.(.(.((((..((((((...((((........((((.......)))).........))))..))))))..)).))).).)))............"
    # s2 = "...........((((((((.......((.....)).......))))).)))........((((...((((((.(.......).))))))......))))."

    # sequence = "GAAAUGGUCCUACUCCGUAGCUUAAUGAUCCGAAGUCAGUAUCUGAUACACUCUUCAGGCCGGUAGCCCACAUCGCGCAGCCUUGUCGAUGAUCCACGCU"
    # s1 = "...((((.......)))).(((...((..((((((...((((...))))...)))).)).))..)))...(((((.(((....))))))))........."
    # s2 = "....(((..(((((........((.((((.....)))).))(((((........)))))..))))).)))(((((.(((....))))))))........."

    # sequence = "CAAUUGCGCCGCAUAAAACCAGAGAUUCUACCCUCAAUCGGUUCUUAAGACGUACUGCGCGUUUCACCAGACCACAAUGCAGGGCGGCACCGUUAGGCAA"
    # s1 =       ".......(((((.........((((....(((.......))))))).....................................))))).((....))..."
    # s2 =       ".......(((((.........((((....(((.......))))))).((((((.....))))))...................))))).((....))..."

    search_width = 20

    # sws = 20
    sws = 1
    # sws = False

    Verbose = True

    limit=False
    limit=0.003


    Debug = False
    # Debug = 1
    # Debug = 1
    # Debug = 8
    # Debug = -1

    # paths = find_path(sequence, s1, s2, indirect_iterations=1, add_moves=[],
    #                   search_width=search_width, Debug=Debug, Verbose=Verbose, limit=limit)
    # for s, sw, mode, moves in paths:
    #     print_moves(sequence, s1, s2, moves, Verbose=Verbose)
    #     break

    # print ([(i[0], i[1]) for i in paths])

    better = []
    worse = []
    collect = []

    filename_samples = 'dataset_100_large.csv'
    # filename_samples = 'dataset_200_large.csv'
    # filename_samples = 'dataset_300_large.csv'

    # better [16, 25, 31, 57, 64, 66, 76, 78, 81, 91, 99, 128, 138, 146, 152, 165, 168, 171, 173, 175, 197, 203, 230, 231, 232, 237, 239, 242, 256, 270, 272, 283, 292, 303, 311, 319, 321, 331, 347, 364, 375, 377, 382, 396, 401, 409, 420, 423, 424, 436, 440, 451, 453, 463, 468, 469, 470, 475, 479, 481, 489, 490, 500, 518, 527, 532, 535, 545, 562, 583, 589, 595, 613, 615, 625, 631, 640, 646, 651, 675, 680, 682, 721, 725, 731, 739, 755, 761, 763, 773, 793, 799, 801, 803, 806, 814, 817, 818, 821, 832, 844, 850, 855, 858, 864, 870, 873, 882, 895, 912, 913, 925, 931, 947, 948, 950, 966, 974, 978, 981, 992]
    # worse [80, 297, 367, 616, 637, 693, 743, 777, 815, 904, 915]
    # ! 693, 815, 915

    # better [2, 5, 9, 23, 25, 29, 35, 36, 41, 45, 51, 52, 55, 56, 57, 60, 61, 62, 64, 67, 71, 76, 77, 79, 81, 83, 88, 90, 91, 95, 96]
    # worse [0, 3, 4, 7, 8, 11, 12, 14, 15, 17, 19, 22, 24, 26, 30, 31, 34, 38, 40, 42, 44, 47, 48, 54, 58, 65, 66, 68, 69, 70, 72, 78, 82, 84, 86, 87, 92, 94, 97, 99]

    samples_df = pd.read_csv(filename_samples)
    for index, row in samples_df.iterrows():
        

        if index == 1000:
            break

        if index < 750:
            continue

        # if index != 417:
        #     continue
        # if bp_dist > 20:
        #     continue

        # if index == 20:
        #     break


        if not Debug:
            sequence = row.sequence
            s1 = row.s1
            s2 = row.s2
        if Debug:
            if Debug != -1:
                sequence = row.sequence
                s1 = row.s1
                s2 = row.s2
                if index!=Debug: continue

            if Debug == -1:
                if index==1:
                    break


        # sequence = "CAAUUGCGCCGCAUAAAACCAGAGAUUCUACCCUCAAUCGGUUCUUAAGACGUACUGCGCGUUUCACCAGACCACAAUGCAGGGCGGCACCGUUAGGCAA"
        # s1 =       ".......(((((.........((((....(((.......))))))).....................................))))).((....))..."
        # s2 =       ".......(((((.........((((....(((.......))))))).((((((.....))))))...................))))).((....))..."

        # sequence = "AACGGGGGCUUCAACUCGCUCAGAAUCAGCGGUAUAGAUAUCCGGGUAGCGGCUUAAAGCAGCACUUUACCAUCGAGGGGGCAAGGAACACUAGCCGACU"
        # s1       = "..((((...((..((.((((.......))))))..))...))))(((((.((((......))).).)))))....((..(((.((.....)).)))..))"
        # s1       = "..((((....((...(((((.......)))))....))..))))(((((.(.(........)..).)))))....((..(((.((.....)).)))..))"
        # s2       = "..((((....((...(((((.......)))))....))..))))((((.(.((.....)).).....))))....((..(((.((.....)).)))..))"

        # sequence = "UUCUAGGAACCACCCAACGUCCUAAUAGCUCGAUGGAUGAUCAGGAUGCAUAUCACGCGACGUCCCGAACCUAUACAGACCGGGCAGAACACAGCCACUU"
        # s1 = "..................................................................................(((........)))...."
        # s2 = ".................................((.....((.((((((.........).))))).)).......)).....(((........)))...."


        bp_dist = RNA.bp_distance(s1, s2)
        if sws:
            search_width = int(bp_dist * sws)
        
        print (index, bp_dist, Debug)

        en1, en2 = 0, 0

        output_path = False
        # output_path = True

        limit=False        
        hopp=False

        start = time.time()
        en1, paths1, c1 = find_path(sequence, s1, s2, indirect_iterations=1, add_moves=[],
                        search_width=search_width, Debug=Debug, Verbose=Verbose, limit=limit, hopp=hopp)
        t1 = time.time()-start

        if Debug:
            for s, sw, mode, moves in paths1:
                print_moves(sequence, s1, s2, moves, Verbose=output_path)
                mp = print_moves(sequence, s1, s2, moves, midpoint=True)
                print ("mid", mp)
                break
        # en1 = 0

        limit=0.05
        hopp=False

        start = time.time()
        en2, paths2, c2 = find_path(sequence, s1, s2, indirect_iterations=1, add_moves=[],
                        search_width=search_width, Debug=Debug, Verbose=Verbose, limit=limit, hopp=hopp)
        t2 = time.time()-start

        # print ("en2", en2, paths2)

        if Debug:
            for s, sw, mode, moves in paths2:

                print ("here", moves)

                print_moves(sequence, s1, s2, moves, Verbose=True)
                break



        # print (moves)
        if en2 < en1:
            better.append(index)
        if en1 < en2:
            worse.append(index)
        collect.append((en1, c1, t1, en2, c2, t2))

        print(f"sequence = \"{sequence}\"")
        print(f"s1 = \"{s1}\"")
        print(f"s2 = \"{s2}\"")
        print (bp_dist)

        # print ('regular', en1, 'with limit', en2, "limit2", en3)
        print ('regular', en1, 'with limit', en2)
    
    print ("better", better)
    print ("worse", worse)

    df = pd.DataFrame(collect)
    df.to_csv(f"./results/{filename_samples}_midfp4")
    print (df)