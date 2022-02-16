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

sw = 10
# seqid = 806

seqid = 46
# seqid = 171
# seqid = 2820
# seqid = 2187
# seqid = 4
# seqid = 8192

# seqid = 1188
# seqid = 3770
# seqid = 1733 #<--!
seqid = 14
# seqid = 9050

# 1188 ADD extend error -6.2 -6.2 -6.3 (62, 88) (62, 88) (65, 84)
# 1476 ADD extend error -9.0 -9.0 -9.5 (37, 45) (37, 45) (22, 85)
# 1733 ADD extend error -10.6 -10.6 -11.0 (57, 71) (57, 71) (48, 87)
# 2579 ADD extend error -13.1 -13.1 -14.7 (44, 64) (44, 64) (38, 73)
# 3191 ADD extend error -10.7 -10.7 -11.2 (4, 62) (4, 62) (66, 95)
# 3191 ADD extend error -10.8 -10.8 -11.2 (66, 95) (66, 95) (6, 60)
# 3736 ADD extend error -9.0 -9.0 -10.5 (18, 34) (18, 34) (13, 39)
# 3770 ADD extend error -5.1 -5.1 -5.2 (81, 87) (81, 87) (22, 61)
# 4122 ADD extend error -9.2 -9.2 -9.4 (7, 82) (7, 82) (21, 44)
# 4561 ADD extend error -6.5 -6.5 -6.6 (14, 80) (14, 80) (6, 89)
# 6466 ADD extend error -27.0 -27.0 -27.3 (12, 19) (13, 18) (40, 66)
# 7197 ADD extend error -9.6 -9.6 -9.9 (14, 66) (14, 66) (36, 46)
# 7493 ADD extend error -14.2 -14.2 -15.3 (30, 73) (30, 73) (32, 71)
# 7579 ADD extend error -16.2 -16.2 -16.8 (30, 95) (29, 96) (82, 91)
# 8532 ADD extend error -14.0 -14.0 -14.3 (13, 35) (13, 35) (18, 29)
# 9050 ADD extend error -3.7 -3.7 -4.3 (76, 91) (76, 91) (70, 97)

# seqid = 1920 # 100 set

startid = 0


filename_samples = f'./dataset_100_large.csv'
# filename_samples = f'./dataset_200_large.csv'

debug = False
debug = True

# counter = []

counter_add_begin = []
counter_add_extension = []
counter_add_adj = []

counter_add_hopp = []
counter_del_hopp = []

counter_del_begin = []
counter_del_begin_pred = []
counter_del_extension = []
counter_del_adj = []

def printd(*args):
    if debug: print (args)

samples_df = pd.read_csv(filename_samples)


for index, row in samples_df.iterrows():
    if index < startid: 
        continue
    if debug and index != seqid:
        continue
    
    if index == 20:
        break

    
    percent_complete = 100-(len(samples_df)-index)/len(samples_df)*100
    bar_length = 20
    filled_length = int(percent_complete/100*bar_length)
    rest = bar_length - filled_length
    bar = "█" * filled_length + '_' * rest
    print(
        f'\rComputing |{bar}| {percent_complete:.1f}% complete {index} ', end="\r")
    

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

    # 171
    # sequence = "GAAAUGAUACCUCAUUGCAGCUCCAGCGCGCUGAGACACUGUGAGGAUGGGGGGCACAAUUGCUCCUUGUCGAGAGGUCGCAUGGUAGGUCACUUGACUA"
    # s2 = "..........((((((.((((........)))).)).....))))(((((((((((....)))))))))))....((((........))))........."
    # s1 = "..........((((((.((((........)))).)).....))))(((((((((((....)))))))))))................((((....))))."

    #1920#
    # sequence = "UAGGUGCUUGUGACAGGUGGGCACGUUUCGCACAUAUGUUCUGGCGGCUA"
    # s1       = "(((.((.......................(((....))).....)).)))"
    # s2       = "(((.(((((((.....)))))).......(((....)))......).)))"

    # non-adj 6407
    # sequence = "UUCUCCCUCUUGUUCCUUGAACAGAUCUCAUGGCGGGGUGGCUAGCAGCCUAACUAUCUCACGUGGUCCAGGUCUGCCGUCCCAACUGGACAUGCCGGCG"
    # s1 = "..((((((((.((((...)))))))......)).)))..(((.((..((((.(((((.....)))))..)))))))))......................"
    # s2 = ".......(((.((((...))))))).........((((((((.((..((((.(((((.....)))))..))))))))))))))................."

    # non-adj example add
    # sequence = "GCAUUGAUUACGCAACUUCAGAGUGAGUCCCAGGAACGUCCCUUGCG"
    # s1       = "((.((((..........)))).))........(((...)))......"
    # s2       = "...((((..........)))).(((((.....(((...))))))))."

    # sequence = 'UGUGUGUUGGGAAACUUGGAAGCUAUUACAGUAUGGUCUA'
    # s1       = '.....(((....)))......(((((......)))))...'
    # s2       = '...((........)).((((..((((......))))))))'

    # sequence = 'ACUAAUGGUUUUUCCUUGCGUAAUUCCUGAGGGAGAUUGUCUGACCCGGCAUAGGCGACG'
    # s1       = '.........(((((((((.(.....).)))))))))(((((((........)))))))..'
    # s2       = '.........(((((((((.(.....).)))))))))..(((.(.((.......)))))).'

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


        
    
    # sorted moves (for clustering)
    p = [(0, 0)] + [(i[0], i[1]) for i in path]
    if debug:
        printd(f"sequence = \"{sequence}\"")
        printd(f"s1 = \"{s1}\"")
        printd(f"s2 = \"{s2}\"")
        printd (RNA.bp_distance(s1, s2))
        print_moves(sequence, s1, s2, path, convert_to_float=True)
        print (p)


    p.sort(key=lambda x: (x[0], x[1]))
    clustered_moves, clabels = cluster_moves(p, absolute_values=False)
    clustered_moves_abs, clabels_abs = cluster_moves(p, absolute_values=True)

    # m = [(i[0], i[1]) for i in inputmoves]
    clustered_moves_i = {i[0][0]:i[1] for i in zip(p, clabels)}
    clustered_moves_j = {i[0][1]:i[1] for i in zip(p, clabels)}

    clusters = clabels[-1]+1
    used_clusters = [0]*(clusters)
    last_cluster_move = defaultdict(list)


    # if debug:
    #     printd (index)
    #     print (clustered_moves, len(clustered_moves))
    #     print (clabels, len(clabels), clusters, used_clusters)





    lastcluster = -1
    current_cluster = []
    starter_cluster = defaultdict(list)

    cluster_list = defaultdict(list)
    interacting_clusters = defaultdict(set)

    # interacting clusters
    for (a,b),c in clustered_moves.items(): 
        cluster_list[c].append((a,b))
        # print ("current", c, (a,b))
        if -a in clustered_moves_i:
            # print ("interact", clustered_moves_i[-a])
            interacting_clusters[c].add(clustered_moves_i[-a])
        if -b in clustered_moves_i:
            # print ("interact", clustered_moves_i[-b])
            interacting_clusters[c].add(clustered_moves_i[-b])
        if -a in clustered_moves_j:
            # print ("interact", clustered_moves_j[-a])
            interacting_clusters[c].add(clustered_moves_j[-a])
        if -b in clustered_moves_j:
            # print ("interact", clustered_moves_j[-b])
            interacting_clusters[c].add(clustered_moves_j[-b])

    # print (interacting_clusters)
    # print (cluster_list)

    # create an initial path for every helix and check where it starts / ends. can we predict a helix greedily?

    for (a,b),c in clustered_moves.items(): 
        
        #  new cluster
        if c != lastcluster and current_cluster != []:
            # print ("current cluster", c, current_cluster)
            # print ()

            # for this method, only look at delete clusters. 
            if current_cluster[0][0]<0:





                ptend = ptmid.copy()
                for i,j in current_cluster:
                    ptend[-i] = -j
                    ptend[-j] = -i

                # interacting clusters?
                ptmid2 = ptmid.copy()
                # print ("interact", interacting_clusters[c-1])
                for related_cluster in interacting_clusters[c-1]:
                    for k, l in cluster_list[related_cluster]:
                        ptmid2[k] = l
                        ptmid2[l] = k
                        # print (k, l)
                        


                send = RNA.db_from_ptable(ptend)
                smid = RNA.db_from_ptable(ptmid2)

                fp = findpath.findpath_single(sequence, smid, send, search_width_multiplier=sw, mp=True)
                tpath = fp.get_path()
                e = fp.get_en()
                

                pt = ptmid.copy()
                tsE = fc.eval_structure_pt(pt)

                # compare with greedy path
                greedy_path = []
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
                    
                    greedy_path.append((i,j))
                    pt = adjust_pt(pt, i, j)


                    en = fc.eval_structure_pt(pt)
                    # print (i,j,en, RNA.db_from_ptable(pt))
                    if en>tsE:
                        tsE = en
                    current_cluster.pop(pos)

                # if tsE!=e:
                #     print_moves(sequence, smid, send, tpath, convert_to_float=True)
                #     print ("greedy:", tsE, e)
                #     quit()

                # print ("greedy:", tsE, e, c)
                # starter_cluster[c-1].append(tpath[-1]) # last iteration c was correcet


                starter_cluster[c-1].append(greedy_path[-1])

                # if debug:
                #     print ("greedy path", greedy_path)
                #     print_moves(sequence, smid, send, tpath, convert_to_float=True)

            
            current_cluster = []

        lastcluster = c
        current_cluster.append((a,b))

    # quit()

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

    lastaddcluster = -1
    lastdelcluster = -1

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

        current_label = clustered_moves[(nexti,nextj)]
        members = avail_clusters[current_label]
        # print ("iteration", e, nexti, nextj, "label", current_label, members)

        # this checks if the optimal result can be found deleting helices with only removing the end of a helix stack
        if ce !=0 and used_clusters[current_label]!=0:

            # for (i, j, nexten) in avail_moves:

            ptnew = adjust_pt(pt, nexti, nextj)
            snew = RNA.db_from_ptable(ptnew)
            res0, _ = fp_call(sequence, snew, s2, sw)

            i,j = members[0]
            ptnew = adjust_pt(pt, i, j)
            snew = RNA.db_from_ptable(ptnew)
            res1, test1 = fp_call(sequence, snew, s2, sw)

            i,j = members[-1]
            ptnew = adjust_pt(pt, i, j)
            snew = RNA.db_from_ptable(ptnew)
            res2, test2 = fp_call(sequence, snew, s2, sw)

            if res0<min(res1,res2) and min(res1,res2)>sE:
                printd ("adj. err removal, valid:",  members[0],  members[-1], "chosen", (nexti,nextj))
                printd ("adj. err removal", i, j, res0, res1, res2, members)
                if debug: print_moves(sequence, snew, s2, test1, convert_to_float=True)
                
                counter_del_adj.append(index)
                # quit()

            else:
                # pass
                printd ("helix removal ok", res0, res1, res2, "/ alt min:", min(res1,res2), "/ sE:", sE)

        # current_label = clustered_moves[(nexti,nextj)]
        # members = avail_clusters[current_label]

        # helix remove

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
                printd ("removal: new helix", (nexti,nextj), "vs", en_pos, en_pos==(nexti,nextj))
            else:
                printd ("removal: extend helix", (nexti,nextj), "vs", en_pos, en_pos==(nexti,nextj))
            # print ("shift", current_label, current_label in shifts)

            # the greedy choice differs from the actual choice. the initial step is the most critical one. 
            if en_pos!=(nexti,nextj):

                i,j = en_pos
                ptnew = adjust_pt(pt, i, j)
                snew = RNA.db_from_ptable(ptnew)
                res1, test1 = fp_call(sequence, snew, s2, sw)
                # print_moves(sequence, snew, s2, test1, convert_to_float=True)

                if used_clusters[current_label]==0:
                    i = -starter_cluster[current_label][0][0]
                    j = -starter_cluster[current_label][0][1]
                    ptnew = adjust_pt(pt, i, j)
                    snew = RNA.db_from_ptable(ptnew)
                    res2, test0 = fp_call(sequence, snew, s2, sw)
                else:
                    res2 = np.inf

                # res2 = np.inf

                # print ("helix removal", (nexti,nextj), "vs", en_pos, "/",  en_pos==(nexti,nextj))




                if res1>sE and used_clusters[current_label]!=0:
                    # print  (index, "helix ext error", en_pos, res1, (i,j), res2, "/", (nexti, nextj), res0, sE)
                    # if debug: print_moves(sequence, snew, s2, test0, convert_to_float=True)                    
                    counter_del_extension.append(index)

                # take the greedy and predicted starting positions (2 options as starting pos)
                if res1>sE and res2>sE and used_clusters[current_label]==0:
                    # if debug:
                        # print ("new pred helix error", (nexti, nextj), "greedy:", en_pos, "pred", starter_cluster[current_label][0])
                    counter_del_begin_pred.append(index)
                
                # what happens if we only start with the greedy option (not very good...)
                if res1>sE and used_clusters[current_label]==0:
                    # if debug: 
                        # print ("new helix error, greedy:", res1, en_pos, res2, starter_cluster[current_label][0], "sE:", sE)
                        # print_moves(sequence, snew, s2, test1, convert_to_float=True)     
                    counter_del_begin.append(index)
                else:
                    printd ("ok", en_pos, res1, (i,j), res2, sE)
                

            # if lastdelcluster > -1:
            #     print ("DEL", lastdelcluster, lastdelcluster != current_label, len(avail_clusters[lastdelcluster]))
            #     print()



            if lastdelcluster > -1 and lastdelcluster != current_label and len(avail_clusters[lastdelcluster])!=0:
                if avail_clusters[current_label][0][0] < 0 and avail_clusters[lastdelcluster][0][0] < 0:
                    
                    # theres still moves left in the last cluster - cant we just complete it?
                    
                    printd("DEL different")
            
                    lmembers = avail_clusters[lastdelcluster]
                    en_min = np.inf
                    lpos = -1
                    for i,j in lmembers:
                        next_en = fc.eval_move_pt(pt, i, j)
                        en = np.interp(next_en/100, [-10,10], [0,1]) 
                        if en<en_min:
                            en_min = en
                            lpos = i,j

                    i,j = lpos
                    ptnew = adjust_pt(pt, i, j)
                    snew = RNA.db_from_ptable(ptnew)
                    res1, test0 = fp_call(sequence, snew, s2, sw)

                    if res1>sE:
                        print  (index, "DEL extend error", res1, res1, sE, (nexti, nextj), en_pos, lpos)
                        counter_del_hopp.append(index)
                        # quit()
                    else:
                        printd  ("DEL extend OK", lpos, "vs", en_pos)
                # else:
                #     printd  ("ADD extend OK", lpos, "vs", en_pos)

            # as long as a cluster is not used up...
            # if used_clusters[current_label]==0: 
            # if len(avail_clusters[current_label])>0:
            lastdelcluster = current_label



        # helix addition: check if greedy = best result?

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
                printd ("add new helix", (nexti,nextj), "vs", en_pos, en_pos==(nexti,nextj))
            else:
                printd ("helix extension", (nexti,nextj), "vs", en_pos, en_pos==(nexti,nextj))

                # 1 below min, 1 above move
                min_ext = min(last_cluster_move[current_label])-1
                max_ext = max(last_cluster_move[current_label])+1


                # find all adjacent moves
                adjacent = set()
                for s in last_cluster_move[current_label]:
                    adjacent.add(abs(s)+1)
                    adjacent.add(abs(s)-1)

                # is the greedy move adjacent to existing helix moves? (only valid for helix extension)
                
                nextig = abs(en_pos[0])
                if abs(nextig) not in adjacent:

                    # last move not adjecent - try if the greedy move leads to the best result
                    # see if greedy move would have been fine
                    # we migh have started the helix in a weird starting position, check if one adjacent move leads us to the goal
                    
                    printd ("fail", nextig, last_cluster_move[current_label], "exp:", min_ext, max_ext, en_pos, (nexti,nextj))
                    
                    fail = 0
                    for i,j in members:
                        if abs(i) in adjacent:
                            # print ("try", i)

                            ptnew = adjust_pt(pt, i, j)
                            snew = RNA.db_from_ptable(ptnew)
                            res0, test0 = fp_call(sequence, snew, s2, sw)
                            if debug:
                                print ("fail", res0, sE)
                                print_moves(sequence, snew, s2, test0, convert_to_float=True)
                            if res0<=sE: # one adjacent move at least is ok
                                fail = 1
                    
                    if fail == 0:
                        counter_add_adj.append(index)

            # check if we stay inside a add cluster, can we finish properly?

            if lastaddcluster and lastaddcluster != current_label and len(avail_clusters[lastaddcluster])!=0:
                if avail_clusters[current_label][0][0] > 0 and avail_clusters[lastaddcluster][0][0] > 0:
                    
                    # theres still moves left in the last cluster - cant we just complete it?
                    
                    printd("ADD different")
            
                    lmembers = avail_clusters[lastaddcluster]
                    en_min = np.inf
                    lpos = -1
                    for i,j in lmembers:
                        next_en = fc.eval_move_pt(pt, i, j)
                        en = np.interp(next_en/100, [-10,10], [0,1]) 
                        if en<en_min:
                            en_min = en
                            lpos = i,j

                    i,j = lpos
                    ptnew = adjust_pt(pt, i, j)
                    snew = RNA.db_from_ptable(ptnew)
                    res1, test0 = fp_call(sequence, snew, s2, sw)

                    if res1>sE:
                        print  (index, "ADD extend error", res1, res1, sE, (nexti, nextj), en_pos, lpos)
                        counter_add_hopp.append(index)
                        # quit()
                    else:
                        printd  ("ADD extend OK", lpos, "vs", en_pos)
                # else:
                #     printd  ("ADD extend OK", lpos, "vs", en_pos)

            # as long as a cluster is not used up...
            # if used_clusters[current_label]==0: 
            # if len(avail_clusters[current_label])>0:
            lastaddcluster = current_label



            # if the next i is not in the first energy position of the cluster (non-greedy choice was computed by findpath)
            if en_pos!=(nexti,nextj):

                i,j = en_pos
                ptnew = adjust_pt(pt, i, j)
                snew = RNA.db_from_ptable(ptnew)
                res1, test0 = fp_call(sequence, snew, s2, sw)

                # if res0<res1 and res1>sE:
                if res1>sE:
                    printd  (index, "error", res1, res1, sE, (nexti, nextj), en_pos)
                    # print_moves(sequence, snew, s2, test0, convert_to_float=True)
                    if used_clusters[current_label]==0:
                        printd ("error: add new helix", (nexti,nextj), "vs", en_pos, en_pos==(nexti,nextj))
                        if debug: print_moves(sequence, snew, s2, test0, convert_to_float=True)
                        counter_add_begin.append(index)
                    else:
                        counter_add_extension.append(index)
                
                else: # Findpath had a non-greedy choice, but the greedy choice gets us also to the optimal result
                    printd ("ok", res1, res1, sE)




        used_clusters[current_label] += 1
        last_cluster_move[current_label].append(abs(nexti))




            # members = [members[0]] + [members[1]] # first and last


print ("error counter")

counter_add_begin = list(set(counter_add_begin))
counter_add_extension = list(set(counter_add_extension))
counter_add_adj = list(set(counter_add_adj))

counter_del_begin = list(set(counter_del_begin))
counter_del_extension = list(set(counter_del_extension))
counter_del_adj = list(set(counter_del_adj))
counter_add_hopp = list(set(counter_add_hopp))

print("counter_add_begin", len(counter_add_begin), counter_add_begin)
print("counter_add_extension", len(counter_add_extension), counter_add_extension)
print("counter_add_adj", len(counter_add_adj), counter_add_adj)

print("counter_add_hopp", len(counter_add_hopp), counter_add_hopp)
print("counter_del_hopp", len(counter_del_hopp), counter_del_hopp)

print("counter_del_begin", len(counter_del_begin), counter_del_begin)
print("counter_del_begin_pred", len(counter_del_begin_pred), counter_del_begin_pred)
print("counter_del_extension", len(counter_del_extension), counter_del_extension)
print("counter_del_adj", len(counter_del_adj), counter_del_adj)