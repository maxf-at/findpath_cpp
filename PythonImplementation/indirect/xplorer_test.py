#!/usr/bin/env python3
# coding: utf-8

from operator import __truediv__, truediv
import time
import subprocess
import sys

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json

from dataclasses import dataclass

from sklearn.metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS, pairwise_distances
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA


import hdbscan

import findpath
import pathfinder
import pathfinder_i_2
import pathfinder_i_cache


import RNA
# import RNAxplorer


# coloredlogs
import coloredlogs
import logging

def print_d(*args):
    """
    debug print for coloredlogs library
    """
    msg = ""
    for m in args:
        msg += str(m)
        msg += " "
    logging.info(msg)


# @dataclass(frozen=True) # autohash
@dataclass
class structure_node:
    __slots__ = ('__dict__', 'structure', 'bp_dist_s1', 'bp_dist_s2', 'en')
    structure: str
    bp_dist_s1: int
    bp_dist_s2: int
    en: int

    max_en_s1: int = 0
    max_en_s2: int = 0

    def __str__(self):
        return f'{self.structure} -> {self.en:2.2f}'

    def __eq__(self, other):
        "Checks if equal."
        return self.structure == other.structure

    def __hash__(self):
        """Make the instances hashable."""
        return hash(self.structure)




def callRNAxplorer(seq, ref_struct1, ref_struct2, s1, s2, n=500):
    """
    Call RNAxplorer using ref_struct1, ref_struct2 as
    structures that are made equally probable within
    the ensemble. Tell RNAxplorer to generate n samples

    returns a Matrix2D object which contains all unique samples and the inputdata of the xplorer call.
    """
    RNA.cvar.uniq_ML = 1
    RNA.cvar.betaScale = 1

    # md = RNA.md("global")
    md = RNA.md()
    md.circ = 0
    md.uniq_ML = 1
    md.compute_bpp = 0
    md.betaScale = 1


    vc = RNA.fold_compound(seq, md, RNA.OPTION_MFE | RNA.OPTION_PF)
    extended_options = ""

    # samples = list (k,l,energy,[strucutres])
    tmpSamples = RNAxplorer.estimate_landscape(
        vc, ref_struct1, ref_struct2, n, extended_options)

    return_structures = set()

    # convert tmpSamples to a list with one structure per k,l-cell (min. energy for example).
    xplorerSamples = []
    for s in tmpSamples:
        uniqueStructures = list(set(s[3]))
        minEnergy = float("inf")
        minStructure = ""
        for structure in uniqueStructures:

            energy = int(RNA.energy_of_struct(seq, structure)*100)
            # return_structures.add((structure,energy))

            return_structures.add(structure_node(structure, RNA.bp_distance(structure, s1), RNA.bp_distance(structure, s2), energy))

            if energy < minEnergy:
                minStructure = structure
                minEnergy = energy
        xplorerSamples.append((s[0], s[1], minEnergy, minStructure))

    return list(return_structures)



def detect_local_minimum(fc, structure):
    # perform gradient walk from sample to determine direct local minimum
    pt = RNA.IntVector(RNA.ptable(structure))
    fc.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
    return RNA.db_from_ptable(list(pt))





def detour(pt, pt1, pt2):
    # fix: should also include missing bps 

    detour_i = set()
    moves_pt1 = set()
    moves_pt2 = set()

    for i in range(1,pt[0]):
        if pt[i]>0 and i<pt[i]:
            if pt[i]!=pt1[i] and pt[i]!=pt2[i]:
                detour_i.add((i, pt[i]))

        if (i < pt1[i] and pt1[i]==pt2[i] and pt[i]==0):
            # print ("delete", (-i, -pt1[i]))
            detour_i.add((-i, -pt1[i]))

        # distance to pt1
        if pt[i]!=pt1[i]:
            if i<pt1[i]:
                moves_pt1.add((-i, -pt1[i]))
            if i<pt[i]:
                moves_pt1.add((i, pt[i]))
        # distance to pt2
        if pt[i]!=pt2[i]:
            if i<pt2[i]:
                moves_pt2.add((-i, -pt1[i]))
            if i<pt2[i]:
                moves_pt2.add((i, pt[i]))

    return detour_i, moves_pt1, moves_pt2


def direct_indirect_filter(pt, pt1, pt2):

    direct_pt = pt.copy()

    for i in range(1,pt[0]):
        if pt[i]>0 and i<pt[i]:
            if pt[i]!=pt1[i] and pt[i]!=pt2[i]:
                direct_pt[i] = 0
                direct_pt[pt[i]] = 0


    return direct_pt



def stabilize(s, detour_i, fc):

    s_eval = fc.eval_structure(s)

    pt = list(RNA.ptable_from_string(s))

    for i in detour_i:
        j = pt[i]
        pt[i] = 0
        pt[j] = 0
    
    s_new = RNA.db_from_ptable(pt)
    s_eval_new = fc.eval_structure(s_new)

    # print (s, s_eval)
    # print (s_new, s_eval_new)

    return s_eval_new-s_eval



def common_structure(pt1, pt2):
    bp_common = []
    pt_common = pt1.copy()
    for i in range(1,pt1[0]):
        if (pt1[i] == pt2[i] and pt1[i] != 0):
            pt_common[i] = pt1[i]
        else:
            pt_common[i] = 0



def remove_incompatible_bps(pt, l, i, j):
    i = abs(i)
    j = abs(j)
    # print ("adj1", RNA.db_from_ptable(pt))

    adjusted_bps = 0
    if (l[i] != l[j] or pt[i] != 0 or pt[j] != 0):
        
        if pt[i] != 0:
            b = pt[i]
            pt[b] = 0
            pt[i] = 0
            adjusted_bps += 1

        if pt[j] != 0:
            b = pt[j]
            pt[b] = 0
            pt[j] = 0
            adjusted_bps += 1
        
        # print ("adj2", RNA.db_from_ptable(pt), i, j)

        # update loop table
        l = RNA.loopidx_from_ptable(pt)
        
        # if we're still incompatible, move in
        if (l[i] != l[j]):      
            if (i+2 == j-1):
                adjusted_bps += remove_incompatible_bps(pt, l, i+1, j)
            else:
                adjusted_bps += remove_incompatible_bps(pt, l, i+1, j-1)
                
    return adjusted_bps


def cluster_intermediates(sequence, s1, s2, intermediates, Debug=False):


    if Debug:
        coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s \033[96m %(message)s \033[0m'
        # coloredlogs.install(level='VERBOSE')
        coloredlogs.install(level='INFO')


    fc = RNA.fold_compound(sequence)
    s1_eval = fc.eval_structure(s1)
    s2_eval = fc.eval_structure(s2)




    intermediates2 = []
    for structure in intermediates:

        # print_d(structure)
        # structure = detect_local_minimum(fc, structure)
        # print_d(structure)
        # print_d("~~~")

        energy = int(RNA.energy_of_struct(sequence, structure)*100)
        intermediates2.append(structure_node(structure, RNA.bp_distance(structure, s1), RNA.bp_distance(structure, s2), energy))
    
    intermediates2.sort(key=lambda x: x.en)
    intermediates = [i.structure for i in intermediates2]

    # limit intermediates to 1/2 or less
    # intermediates = intermediates[0:int(len(intermediates)/5)]

    # for s in intermediates2[0:30]:
    #     print_d (s)




    search_width_multiplier = 2
    # search_width_multiplier = 10


    # print (f"sequence = '{sequence}'")
    # print (f"s1       = '{s1}'")
    # print (f"s2       = '{s2}'")

    start_time = time.time()
    fp_call = findpath.findpath_class(sequence, True)
    max_en = fp_call.init(s1, s2, search_width_multiplier)
    print (f'direct path: {time.time()-start_time:2.4f} s')
    max_en_direct = max_en

    start_time = time.time()
    # cmd = f'RNAxplorer -M RSH -n 200 --sequence {sequence}'
    # result = subprocess.check_output(cmd, shell=True, encoding="utf8")
    xplorer_time = time.time()-start_time
    print (f'RNAxplorer: {xplorer_time:2.4f} s')

    fp_time = 0

    # compute the direct path (compute distance metrics)    
    bp_dist_total = RNA.bp_distance(s1, s2)
    Verbose = False
    direct_path = pathfinder.pathfinder(sequence, s1, s2, verbose=Verbose, search_width=bp_dist_total*search_width_multiplier)


    # print (direct_path.paths[0])

    direct_structures = set([s1])
    direct_structures_en = []
    s = s1

    for i, j, en in direct_path.paths[0]:
        if i > 0:
            s = s[:i-1] + "(" + s[i:j-1] + ")" + s[j:]
        if i < 0:
            s = s[:-i-1] + "." + s[-i:-j-1] + "." + s[-j:]
        direct_structures.add(s)
        direct_structures_en.append((s,en))
    
    direct_structures_en.sort(key=lambda x:-x[1])

    worst_s =  [x[0] for x in direct_structures_en[0:5]]

    # print (direct_structures_en)
    # print (worst_s)





    # return


    max_en_float = round(max_en/100.0,2)
    # print(max_en_float)

    return_max_en = max_en_float


    # print("max_en:", max_en, max_en_float)

    intermediates = list(set(intermediates))
    intermediates.sort()

    pt1 = list(RNA.ptable_from_string(s1))
    pt2 = list(RNA.ptable_from_string(s2))

    direct_moves = set()
    indirect_moves = set()

    for i in range(1,pt1[0]):
        if pt1[i]!=pt2[i]:
            if i<pt1[i]:
                direct_moves.add((-i, -pt1[i]))
            if i<pt2[i]:
                direct_moves.add((i, pt2[i]))



    
    candidates1 = []
    candidates2 = []

    best_en_2 = 9999

    all_moves = set()


    # every intermediate gets its move set, store in dict
    for i, intermediate in enumerate(intermediates):

        pt = RNA.ptable_from_string(intermediate)
        # print (intermediate, pt)
        det, moves_pt1, moves_pt2 = detour(pt, pt1, pt2)
        
        dist1 = RNA.bp_distance(intermediate, s1)
        dist2 = RNA.bp_distance(intermediate, s2)

        # print (intermediate[0:20], i, "detour:", det, dist1, dist2)
        
        # add current moves to all_moves
        all_moves |= moves_pt1
        all_moves |= moves_pt2     
        indirect_moves |= det

        i_list = []

        max_bp_dist1 = float("-inf")
        max_bp_dist2 = float("-inf")
        max_member_en = float("inf")
        max_s = ""


        intermediate = intermediates[i]
        dist1 = RNA.bp_distance(intermediate, s1)
        dist2 = RNA.bp_distance(intermediate, s2)

        lowest_bp = 9999
        # lowest bp dist to direct path intermediates
        for s_i in direct_structures:
            a = RNA.bp_distance(s_i, intermediate)
            if a < lowest_bp:
                lowest_bp = a


        pt = list(RNA.ptable_from_string(intermediate))
        det, m1, m2 = detour(pt, pt1, pt2)

        en_improvement = []
        en_improvement_inc = []
        add_moves = det.copy()

        # this compares the current intermediate structure to highest
        # energy structures, and verifies if potential moves can "stabilize" a saddle point
        
        for s in worst_s:                        
            temp_pt = list(RNA.ptable_from_string(s))
            
            # print ("worst s", s, add_moves)


            l = RNA.loopidx_from_ptable(temp_pt)
            adjusted_bps = 0
            en_base = fc.eval_structure_pt(temp_pt)
            for i,j in add_moves:
                adjusted_bps += remove_incompatible_bps(temp_pt, l, i, j)
                
            for i,j in add_moves:
                if i>0:
                    temp_pt[i] = j
                    temp_pt[j] = i
                else:
                    temp_pt[-i] = 0
                    temp_pt[-j] = 0

            s = RNA.db_from_ptable(temp_pt)

            en_adj = int(fc.eval_structure(s)*100)
            # en_adj = fc.eval_structure_pt(temp_pt)
            
            en_improvement.append(en_adj-en_base)
            en_improvement_inc.append(adjusted_bps)

        en_improvement = np.mean(en_improvement)
        en_improvement_inc = np.mean(en_improvement_inc)

        add_moves = []
        for (i,j) in det:
            add_moves.append(i)
            add_moves.append(j)

        # detour distance (distance compared to direct path intermediates)
        total_detour = len(det)
        if total_detour < 1: continue # we're on the direct path then

        direct_pt = direct_indirect_filter(pt, pt1, pt2)
        current_en = fc.eval_structure_pt(pt)
        en_direct_pt = fc.eval_structure_pt(direct_pt)
        en_offset = current_en - en_direct_pt

        if dist1 > max_bp_dist1:
            max_bp_dist1 = dist1
        if dist2 > max_bp_dist2:
            max_bp_dist2 = dist2

        if current_en < max_member_en:
            max_member_en = current_en
            max_s = intermediate

        # continue

        # debug - all with 2
        search_width_multiplier = 2
        max_en_s1 = findpath.init_single_findpath(sequence, s1, intermediate, search_width_multiplier, True)
        max_en_s2 = findpath.init_single_findpath(sequence, intermediate, s2, search_width_multiplier, True)
        combined_en = max(max_en_s1, max_en_s2)
        if combined_en < best_en_2:
            best_en_2 = combined_en


        search_width_multiplier = 0.2
        start_time = time.time()
        max_en_s1 = findpath.init_single_findpath(sequence, s1, intermediate, search_width_multiplier, True)
        max_en_s2 = findpath.init_single_findpath(sequence, intermediate, s2, search_width_multiplier, True)
        fp_time += (time.time()-start_time)
        
        
        combined_en = max(max_en_s1, max_en_s2)

        candidates1.append((combined_en, intermediate, pt, add_moves))

        max_en_i = 0
        candidates2.append((max_en_i, intermediate, pt, add_moves))


        

        

        if combined_en < return_max_en:
            return_max_en = combined_en

        



        i_list.append(intermediate)

        total_detour = (dist1 + dist2 - bp_dist_total) / 2

        # print_d (det)
        # print_d (intermediate[0:100], i, "en_i", en_improvement, "/", en_improvement_inc, "detour:", total_detour, "bp dist to d-path:", lowest_bp, "2:", combined_en, "i:", max_en_i, "stabilize en:", en_offset)
        # print_d (intermediate[0:100], i, "en_i", en_improvement, "/", en_improvement_inc, "best:", best_improvement, "stabilize en:", en_offset, "detour:", total_detour, "bp dist to d-path:", lowest_bp, "max_en", max_en, "current_en", current_en, "result:", combined_en)

        # print_d (intermediate[0:110], i, "detour:", total_detour, dist1+dist2, combined_en, "(", max_en_s1, max_en_s2, ")", "en:", current_en, en_offset)
        # print_d (intermediate[0:130], i, "detour:", total_detour, "bp dist to d-path:", lowest_bp, "2:", combined_en, "i:", max_en_i, "stabilize en:", en_offset)
        # print_d (intermediate[0:100], i, "detour:", dist1, dist2, dist1+dist2, combined_en, "(", max_en_s1, max_en_s2, ")", "en:", current_en)

        # print_d (intermediate, "bp_d", total_detour, "max_en", max_en, "current_en", current_en, "result:", combined_en)






    # return

    search_width_multiplier = 100

    sws = [0.5, 1, 1.5, 2]
    remaining_candidates = [50,25,10,1]

    sws = [1, 2]
    remaining_candidates = [10,1]

    candidates1.sort(key=lambda x: x[0])
    candidates2.sort(key=lambda x: x[0])

    # print (candidates1[0])

    Verbose = True
    # Verbose = False

    max_en_1 = 999

    for sw, limit in zip(sws, remaining_candidates):

        candidates = candidates1.copy()
        candidates1 = []

        # print_d("sw:", sw)

        for i, (en, intermediate, pt, add_moves) in enumerate(candidates[0:limit]):        
            
            start_time = time.time()
            max_en_s1 = findpath.init_single_findpath(sequence, s1, intermediate, sw, True)
            max_en_s2 = findpath.init_single_findpath(sequence, intermediate, s2, sw, True)
            fp_time += (time.time()-start_time)

            Verbose = False
            
            # break
            combined_en = max(max_en_s1, max_en_s2)        
            # print_d(i, en, "->", combined_en, intermediate)
            
            candidates1.append((combined_en, intermediate, pt, add_moves))
            if combined_en < max_en_1:
                max_en_1 = combined_en

        candidates1.sort(key=lambda x: x[0])

    print_d ("gradient:", max_en_1, "absolute best:", best_en_2, "direct", max_en_direct, fp_time)

    max_en_2 = 999


    return return_max_en, max_en_1, max_en_2







def init(sequence, s1, s2, search_width_multiplier, Debug=False):


    if Debug:
        coloredlogs.DEFAULT_LOG_FORMAT = '%(levelname)s \033[96m %(message)s \033[0m'
        # coloredlogs.install(level='VERBOSE')
        coloredlogs.install(level='INFO')


    fc = RNA.fold_compound(sequence)
    s1_eval = fc.eval_structure(s1)
    s2_eval = fc.eval_structure(s2)


    max_en = findpath.init_single_findpath(
        sequence, s1, s2, search_width_multiplier, True)
    max_en_float = round(max_en/100.0,2)


    bp_dist_total = RNA.bp_distance(s1, s2)

    # intermediates = callRNAxplorer(sequence, s1, s2, s1, s2, 4)
    # intermediates += callRNAxplorer(sequence, s1, s2, s1, s2, 2)
    # intermediates += callRNAxplorer(sequence, s1, s2, s1, s2, 2)

    cmd = f'RNAxplorer -M RSH -n 200 --sequence {sequence}'
    result = subprocess.check_output(cmd, shell=True, encoding="utf8")
    result = result.split("\n")
    intermediates = []
    for line in result:
        line = line.split()
        if len(line) != 4: continue
        structure = line[1]
        energy = int(RNA.energy_of_struct(sequence, structure)*100)
        intermediates.append(structure_node(structure, RNA.bp_distance(structure, s1), RNA.bp_distance(structure, s2), energy))

    # print_d ("RNAxplorer done")
    # intermediates.sort(key=lambda x: x.en)

    # intermediates = intermediates[0:int(len(intermediates)/3)]

    # for s in intermediates[0:10]:
    #     print_d (s)


    intermediate_structures = [i.structure for i in intermediates]

    Debug = True
    # Debug = False


    new_en = cluster_intermediates(sequence, s1, s2, intermediate_structures, Debug=Debug)


    print_d ("best result direct:", max_en_float, "indirect:", new_en)

    return max_en_float, new_en, intermediate_structures








if __name__ == "__main__":

    # sequence = 'UAAAAUGAUCACGGUUUCAGCUUUGGACAGGGCGUUCCACUAAACUCCUGGUGACAUAGAUAUAUUGGAUUGCAACUACUCGUCGGUCCGGUUGGCGUUC'
    # s1       = '..........(((...((((((..((((..((((...(((((......)))))...........(((.....))).....)))).)))))))))))))..'
    # s2       = '....(((.(((((((....))).....(((((.(((......)))))))))))))))......(((((((((..((.....)))))))))))........'

    # simple 60 nt ex / hard...
    sequence = "GAAAGUUUUCAAGAAAGUAAUAUAGUGUCCCAGUCAGGGGGUAGUGGCGCGUUCUGGACGCAUAGUAUUC"
    s1 = "(((....))).........(((((((((.((((...(.(.((....)).).).)))))))).)).))).."
    s2 = ".................((.(((.((((.((((..(.(..((....)).).).))))))))))).))..."

    # 80 nt
    # sequence = "GGAAGCCGGCGAGGCAGUACCAUUAUAUAGUUUGUCUUCCAAGAAUGGGUACGACCGCGGGACCGUUCGGUUAUCGUCUG"
    # s1 = "((.((((((((.((..((((((((....((.....))......))))).)))..)).)).......)))))).))....."
    # s2 = "((.((((((.(((((((.((.........)))))))))))..((((((...((....))...)))))))))).))....."

    # 80 nt, S: -16.60 simple detour, -17.70 better one
    # sequence = "CGGGAGCGAGUAUGGACCGGGUUGAAUUGGUACCUCGGCCGCCCUCUGGGGGUAGCCGACGGGCGGCUUCAGCCGGGCCC"
    # s1 = ".............((.((.(((((((.(.((.(((((((.(((((...))))).))))).)))).).))))))).)).))"
    # s2 = "((....)).....((.((.(((((((...........(((((((((.(.(.....)))).)))))))))))))).))))."

    # 80 nt, -7.20 direct, -8.30 simple detour, -8.40
    # sequence = "CUGGGUCGUGGUGGCCUUUUAGAUACGAUUCACGAACGUAGCACGUUUCGGUCUCCGGAGACGCAAUGAUCUCGAGGGUA"
    # s1 = ".(((((((((.(((....)))..))))))))).((.(((.((..(((((((...))))).)))).))).))........."
    # s2 = ".(((((((((.............)))))))))(((((((...))))).)).((((.(((..........))).))))..."
    
    # 150 nt, -20.10 kcal/mol direct, indirect -24.50
    # sequence = 'AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAA'
    # s1       = '.............((((((..((...(..(((.....)))..).))..))))))(((((((....(((.((((((.......)))))).....)))....)))))))(((......((.((.........)).))......)))......'
    # s2       = '......(((.((((((...(((....))).(((((((((......))))))...)))((...(((....((((((.......))))))...))).)))))))).)))(((......((.((.........)).))......)))......'

    # sequence = 'AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAA'
    # s1       = '.....(((((.....))))).......(((....(((((......))))))))(((((..(((....((((((((.......)))...............(((...(((.......)))))).........)))))....))).))))).'
    # s2       = '.......((..(((((.....((((.............(((...((.((((...(((((((........((((((.......))))))............)))))))....)))).)).)))........)).))..))).))..))...'

    # init(sequence, s1, s2, True)


    Debug = True

    # Debug = False

    # filename = f'./sample_seqs/indirect_new_150.csv'
    filename = f'./sample_seqs/indirect_new_120.csv'

    df = pd.read_csv(filename)

    start = time.time()

    for index, row in df.iterrows():
        # break
        # 1: 9.3
        # 2: -24.7
        # 18: -23.1

    # 5 check out

        # if index!=1:
        #     continue
        if index>30:
            break


        sequence = row.sequence
        s1 = row.s1
        s2 = row.s2
        e1 = int(row.e1*100)
        e2 = int(row.e2*100)

        intermediate_structures = row.intermediate_structures.replace('\'', '"')
        intermediate_structures = json.loads(intermediate_structures)

        if Debug:
            print (f"sequence = '{sequence}'")
            print (f"s1       = '{s1}'")
            print (f"s2       = '{s2}'")

        # init(sequence, s1, s2, True)
        new_en_all, en1, en2 = cluster_intermediates(sequence, s1, s2, intermediate_structures, Debug=Debug)

        # en1 = en1 - e2
        # en2 = en2 - e2

        print ("i:", index, "/", e1, e2, "/", new_en_all, en1, en2)

    # max_en_s1 = findpath.init_single_findpath(sequence, s1, s2, 2, True)
    # print ("direct:", max_en_s1)

    print ("time elapsed:", round(time.time()-start,4))


