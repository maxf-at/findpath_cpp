#!/usr/bin/env python3
# coding: utf-8


import RNA
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# from PIL import Image
# from IPython.display import SVG, display
# import seaborn as sns

from collections import Counter
from collections import defaultdict
import subprocess
import difflib
import sys
import os
import random
import string
import time

# import tensorflow as tf
# from tensorflow import feature_column
# from tensorflow.keras import layers
from sklearn.model_selection import train_test_split
from sklearn.cluster import DBSCAN
from sklearn import preprocessing
min_max_scaler = preprocessing.MinMaxScaler()

sys.path.append('../')
from pretty_print_path import print_moves
import findpath_librna
import findpath
import helper


# dataset_200_large.csv: 368 447 [504, 623, 806] high

sw = 40
# seqid = 806
seqid = 0


# seqid = 1920 # 100 set

startid = 0

# filename_samples = f'./dataset_100_large.csv'
filename_samples = f'./dataset_200_large.csv'

# debug = False
debug = True

counter = []

def printd(*args):
    print (args)

samples_df = pd.read_csv(filename_samples)
for index, row in samples_df.iterrows():
    if index < startid: 
        continue
    if debug and index != seqid:
        continue
    sequence = row.sequence
    s1 = row.s1
    s2 = row.s2
    # s1 = row.s2
    # s2 = row.s1

    # s1 = ".((((((((.((.(((((.......((((......))))((((((((.......))))))))........))))))).)))))).))............."
    # s2 = "......((((((((((.......................((((((((.......))))))))...............))))).)))))............"

    # 504 problem
    # sequence = "CUGAGCAGGUCAACACAGGCUUCAGCAGGGUGGCAUCAAUGAUGGAGACUACAAAAUGUCCUCUAGGGAGAGGAAAUAACUAGUAGGAAUUGUACAGGACCGAUGUUGUCACUCCAAGGUCCGGCAAUCACGAAUCGAUAAUUGGUCUCGAAGGUCAAUUUUUUCUUGUUGCUCAUUCACUCGCUCUGAAUUACUAGGAC"
    # s1 = "((((((..((....))..))).)))((((((((((((...))))....((((......((((((....))))))........))))..(((((...(((((................))))).)))))..............((((((...)))))).......................))))))))............"
    # s2 = "((((((..((....))..))).)))((((((((((((...))))....((((......((((((....))))))........))))..(((((...(((((................))))).)))))...........................((((........)))).........))))))))............"

    # 806 problem
    # sequence = "GAGAGUCUGGAAUACGCCUCACAGUUCAGUGAGCUGUAAACUUUAUGCGGAUUUCUUUAAGAGUCACAUUGUUGAGCUCCAUACUUCGAUGCGUAGAACGCGAGGAUGGGGGGCGCUAUUCGGCCUCGAGGUAUCUGUCGGCGUCUGUAGCAGCGCACCUAGGGAGUAAAGAACGUAACCCGACGGUGGCUUAAAUACGC"
    # s1 = "(((.(((..(((((((((((..((((((((((..((...(((((...............)))))))..))))))))))((((.((((...((((...)))))))))))))))))).)))))))))))(((.((.((((((((((((.((.(..(........)..))).)).)))....))))))))).)))........"
    # s2 = "(((.(((..(((((((((((...................(((((...............)))))..............((((.((((...((((...)))))))))))))))))).)))))))))))(((.((.((((((((((((...((..(........)..))..)).)))....))))))))).)))........"
    # s2 = "(((.(((..(((((((((((..............((...(((((...............)))))))............((((.((((...((((...)))))))))))))))))).)))))))))))(((.((.((((((((((((.((.(..(........)..))).)).)))....))))))))).)))........"



    # print(f"sequence = \"{sequence}\"")
    # print(f"s1 = \"{s1}\"")
    # print(f"s2 = \"{s2}\"")
    # print (RNA.bp_distance(s1, s2))


    # fp = findpath.findpath_single(sequence, s1, s2, search_width_multiplier=sw, mp=True)
    # result = fp.get_en()/100.0
    # path = fp.get_path()
    

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


    def cluster_eval(i,j): # +global vars...
        

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



    # sequence = row.sequence
    # s1 = row.s1
    # s2 = row.s2

    bp_dist = RNA.bp_distance(s1, s2)
    fc = RNA.fold_compound(sequence)
    pt1 = list(RNA.ptable(s1))
    pt2 = list(RNA.ptable(s2))

    shifts = pt1.copy()
    ptmid = pt1.copy()
    # which move / clusters are shift moves?
    for it, (a,b) in enumerate(zip(pt1, pt2)):
        if it==0: continue
        if a!=b:
            ptmid[it] = 0

        if (a==0 and b!=0) or (a!=0 and b==0) or (a==0 and b==0):
            shifts[it] = 0
        elif (a < it and b > it) or (a > it and b < it):
            shifts[it] = 1
        else:
            shifts[it] = 1

    smid = RNA.db_from_ptable(ptmid)



    fp = findpath.findpath_single(sequence, s1, s2, search_width_multiplier=sw, mp=True)
    sE = fp.get_en()/100.0
    path = fp.get_path()


    if debug:
        printd(f"sequence = \"{sequence}\"")
        printd(f"s1 = \"{s1}\"")
        printd(f"s2 = \"{s2}\"")
        printd (RNA.bp_distance(s1, s2))
        print_moves(sequence, s1, s2, path, convert_to_float=True)
    

    p = [(0, 0)] + [(i[0], i[1]) for i in path]

    p.sort(key=lambda x: (x[0], x[1]))


    clustered_moves, clabels = cluster_moves(p, absolute_values=False)
    clustered_moves_abs, clabels_abs = cluster_moves(p, absolute_values=True)

    clusters = clabels[-1]+1
    used_clusters = [0]*(clusters)

    if debug:
        printd (index)
    # print (path)
        print (clustered_moves, len(clustered_moves))
        print (clabels, len(clabels), clusters, used_clusters)

    #shift clusters
    shift_clusters = set()
    for (a,b),c in clustered_moves.items():        
        if shifts[abs(a)]==1 or shifts[abs(b)]==1:
            shift_clusters.add(c)

    print ("shifts", shifts)
    print (shift_clusters)

    lastcluster = -1
    current_cluster = []

    starter_cluster = defaultdict(list)

    for (a,b),c in clustered_moves.items(): 
        
        if c != lastcluster and current_cluster != []:
            print (c, a,b, current_cluster)
            
            # del -> add
            if current_cluster[0][0]<0:

                ptend = ptmid.copy()
                for i,j in current_cluster:
                    ptend[-i] = -j
                    ptend[-j] = -i
                send = RNA.db_from_ptable(ptend)

                fp = findpath.findpath_single(sequence, smid, send, search_width_multiplier=sw, mp=True)
                tpath = fp.get_path()
                e = fp.get_en()
                

                pt = ptmid.copy()
                tsE = fc.eval_structure_pt(pt)
                # compare with greedy path
                while current_cluster != []:
                    
                    en_min = np.inf
                    en_pos = -1
                    for pos, (i,j) in enumerate(current_cluster):
                        next_en = fc.eval_move_pt(pt, -i, -j)
                        en = np.interp(next_en/100, [-10,10], [0,1]) 
                        if en<en_min:
                            en_min = en
                            en_pos = pos, -i, -j
                    

                    pos, i, j = en_pos
                    pt = adjust_pt(pt, i, j)
                    en = fc.eval_structure_pt(pt)
                    # print (i,j,en, RNA.db_from_ptable(pt))
                    if en>tsE:
                        tsE = en
                    current_cluster.pop(pos)

                if tsE!=e:
                    print_moves(sequence, smid, send, tpath, convert_to_float=True)
                    print ("greedy:", tsE, e)
                    quit()

                # print ("greedy:", tsE, e, c)
                starter_cluster[c-1].append(tpath[-1]) # last iteration c was correcet
                if debug:
                    print_moves(sequence, smid, send, tpath, convert_to_float=True)

            
            current_cluster = []

        lastcluster = c
        current_cluster.append((a,b))

    # print (starter_cluster)


    """
    # cluster findpaths test
    lastcluster = -1
    temppath = []
    for (a,b),c in clustered_moves.items():

        if a>0: continue
        
        if c!=lastcluster:            

            if temppath != []:                  
                
                print (RNA.db_from_ptable(pt1))
                print (RNA.db_from_ptable(pt))
                print (temppath)

                fp = findpath.findpath_single(sequence, s1, RNA.db_from_ptable(pt), search_width_multiplier=sw, mp=True)
                #  = fp.get_en()/100.0
                tpath = fp.get_path()
                print_moves(sequence, s1, RNA.db_from_ptable(pt), tpath, convert_to_float=True)

            pt = pt1.copy()
            lastcluster = c
            temppath = []

        temppath.append((a,b))
        print (a,b, c, temppath)
        pt = adjust_pt(pt, a,b)
    # cluster findpaths test end
    """

    for e, (i,j, cen) in enumerate(path):
        # move on along the predefined path

        if e+1==len(path):
            break

        if e!=0:
            # print ("adj1", RNA.db_from_ptable(ptlast), lastmove[0], lastmove[1])
            pt = adjust_pt(ptlast, i, j)
            # print ("adj2", RNA.db_from_ptable(pt))
        else:
            lastmove = 0,0
            pt = pt1.copy()

        nexti = path[e+1][0]
        nextj = path[e+1][1]
        # print ("iteration", e, i, j, nextmove, used_clusters)
        # print ("s  :", RNA.db_from_ptable(pt), lastmove[0], lastmove[1])

        lastmove = i, j
        ptlast = pt.copy()

        avail_moves = []
        fp_results = []
        en_contrib = []
        avail_clusters = defaultdict(list)
        avail_dict = dict()
        available_add = set()
        available_delete = set()
        lensequence = pt1[0]

        for pos, (i,j) in enumerate(find_moves(pt, pt2)):    
            next_en = fc.eval_move_pt(pt, i, j)
            # mark where we found our move    
            # map energies somehow to [0,1] such that the network will understand.  
            en = np.interp(next_en/100, [-10,10], [0,1]) 
            en_contrib.append(en)
            avail_moves.append((i, j, next_en))
            avail_dict[(i,j)] = next_en

            if i>0:
                available_add.add((i,j))
            elif i<0:
                available_delete.add((i,j))

            current_label = clustered_moves[(i,j)]
            # current_label = clabels[pos]
            avail_clusters[current_label].append((i,j))

            # print ("iteration!", e, pos, i, j, current_label)

        if len(avail_moves)==0:
            add_delete = 0
        else:
            add_delete = len(available_delete)/len(avail_moves)

        # print ('best', en_contrib[0], 'avg', en_mean)

        avail_moves.sort(key=lambda x: x[2]) # best move at [0]

        # print (sE)
        # print (clustered_moves)
        # print (clabels)
        # print ("avail", avail_dict)
        count_available = len(avail_moves)/len(path)
        # print (e, len(avail_dict), avail_clusters)

        ce = cluster_eval(nexti, nextj)

        # if nexti<0:
        #     print ((nexti, nextj), ce)

        # check for alternatives if our choice is not obvious
        if ce !=0:

            current_label = clustered_moves[(nexti,nextj)]
            members = avail_clusters[current_label]
            # for (i, j, nexten) in avail_moves:

            ptnew = adjust_pt(pt, nexti, nextj)
            snew = RNA.db_from_ptable(ptnew)
            res0, _ = fp_call(sequence, snew, s2, sw)

            i,j = members[0]
            ptnew = adjust_pt(pt, i, j)
            snew = RNA.db_from_ptable(ptnew)
            res1, _ = fp_call(sequence, snew, s2, sw)

            i,j = members[-1]
            ptnew = adjust_pt(pt, i, j)
            snew = RNA.db_from_ptable(ptnew)
            res2, _ = fp_call(sequence, snew, s2, sw)

            if res0<min(res1,res2) and min(res1,res2)>sE:
                printd ("err", res0, res1, res2, members)
                
                for (i,j) in members:
                    ptnew = adjust_pt(pt, i, j)
                    snew = RNA.db_from_ptable(ptnew)
                    result_new, _ = fp_call(sequence, snew, s2, sw)
                    # print ((i,j), result_new, avail_dict[(i,j)])
                    fp_results.append(result_new)
                
                counter.append(index)
                # quit()

            else:
                # pass
                printd ("helix removal ok", res0, res1, res2, "/ alt min:", min(res1,res2), "/ sE:", sE)

        current_label = clustered_moves[(nexti,nextj)]
        members = avail_clusters[current_label]

        # helix remove

        # helix extension idea: greedy, adjacent move. 
        if nexti < 0:

            en_contrib = []
            en_min = np.inf
            en_pos = -1
            for i,j in members:
                next_en = fc.eval_move_pt(pt, i, j)
                en = np.interp(next_en/100, [-10,10], [0,1]) 
                if en<en_min:
                    en_min = en
                    en_pos = i,j
                en_contrib.append(en)
                
            if used_clusters[current_label]==0:
                print ("new helix removal", (nexti,nextj), "vs", en_pos, en_pos==(nexti,nextj))
            else:
                print ("extend helix removal", (nexti,nextj), "vs", en_pos, en_pos==(nexti,nextj))
            print ("shift", current_label, current_label in shifts)

            if en_pos!=(nexti,nextj):

                ptnew = adjust_pt(pt, nexti, nextj)
                snew = RNA.db_from_ptable(ptnew)
                res0, _ = fp_call(sequence, snew, s2, sw)

                i,j = en_pos
                ptnew = adjust_pt(pt, i, j)
                snew = RNA.db_from_ptable(ptnew)
                res1, test0 = fp_call(sequence, snew, s2, sw)

                if used_clusters[current_label]==0:
                    i = -starter_cluster[current_label][0][0]
                    j = -starter_cluster[current_label][0][1]
                    ptnew = adjust_pt(pt, i, j)
                    snew = RNA.db_from_ptable(ptnew)
                    res2, test0 = fp_call(sequence, snew, s2, sw)
                else:
                    res2 = np.inf

                # print ("helix removal", (nexti,nextj), "vs", en_pos, "/",  en_pos==(nexti,nextj))


                # if res0<res1 and res1>sE:
                if res1>sE and res2>sE and used_clusters[current_label]!=0:
                    print  (index, "helix ext error", en_pos, res1, (i,j), res2, "/", (nexti, nextj), res0, sE)
                    print_moves(sequence, snew, s2, test0, convert_to_float=True)
                    quit()
                else:
                    print ("ok", en_pos, res1, (i,j), res2, res0, sE)


        # helix add
        # if nexti > 0 and used_clusters[current_label]==0:

        #     en_contrib = []
        #     en_min = np.inf
        #     en_pos = -1
        #     for i,j in members:
        #         next_en = fc.eval_move_pt(pt, i, j)
        #         en = np.interp(next_en/100, [-10,10], [0,1]) 
        #         if en<en_min:
        #             en_min = en
        #             en_pos = i,j
        #         en_contrib.append(en)
                
        #     print ("add new helix", (nexti,nextj), "vs", en_pos, en_pos==(nexti,nextj))

        #     # 1920 error -5.8 -4.8 -5.8

        #     if en_pos!=(nexti,nextj):

        #         ptnew = adjust_pt(pt, nexti, nextj)
        #         snew = RNA.db_from_ptable(ptnew)
        #         res0, _ = fp_call(sequence, snew, s2, sw)

        #         i,j = en_pos
        #         ptnew = adjust_pt(pt, i, j)
        #         snew = RNA.db_from_ptable(ptnew)
        #         res1, test0 = fp_call(sequence, snew, s2, sw)

        #         # if res0<res1 and res1>sE:
        #         if res1>sE:
        #             print  (index, "error", res0, res1, sE, (nexti, nextj), en_pos)
        #             print_moves(sequence, snew, s2, test0, convert_to_float=True)
        #             quit()
        #         else:
        #             print ("ok", res0, res1, sE)


        elif nexti > 0:

            en_contrib = []
            en_min = np.inf
            en_pos = -1
            for i,j in members:
                next_en = fc.eval_move_pt(pt, i, j)
                en = np.interp(next_en/100, [-10,10], [0,1]) 
                if en<en_min:
                    en_min = en
                    en_pos = i,j
                en_contrib.append(en)
                
            if used_clusters[current_label]==0:
                print ("add new helix", (nexti,nextj), "vs", en_pos, en_pos==(nexti,nextj))
            else:
                print ("helix extension", (nexti,nextj), "vs", en_pos, en_pos==(nexti,nextj))

            if en_pos!=(nexti,nextj):

                ptnew = adjust_pt(pt, nexti, nextj)
                snew = RNA.db_from_ptable(ptnew)
                res0, _ = fp_call(sequence, snew, s2, sw)

                i,j = en_pos
                ptnew = adjust_pt(pt, i, j)
                snew = RNA.db_from_ptable(ptnew)
                res1, test0 = fp_call(sequence, snew, s2, sw)

                # if res0<res1 and res1>sE:
                if res1>sE:
                    print  (index, "error", res0, res1, sE, (nexti, nextj), en_pos)
                    print_moves(sequence, snew, s2, test0, convert_to_float=True)
                    quit()
                else:
                    print ("ok", res0, res1, sE)


        used_clusters[current_label] += 1
        print (clustered_moves, nexti, nextj)
        break


            # members = [members[0]] + [members[1]] # first and last


    # break
print (counter)