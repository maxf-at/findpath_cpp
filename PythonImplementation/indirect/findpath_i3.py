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

    def find_stack(self, path, moves):
        """
        find compensation stacks to lower the energy of e.g. saddle points
        """

        print ("launch find_stack")

        sequence = self.sequence

        # if max_pos == 0: return set()
        # if max_pos == len(path): return set()

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
            if (A_str, B_str, C_str) in self.processed_stacks:
                continue
            self.processed_stacks.add((A_str, B_str, C_str))

            # energy of current ptable
            base_energy = self.fc.eval_structure_pt(current_ptable)/100

            # compatible = A.copy()
            # offset_loop_table = A.copy()

            loop_A = self.l_tables(A)
            loop_B = self.l_tables(B)
            current_loop = self.l_tables(current_ptable)

            # for i in range(1, len(A)):
            # if loop_A[i] == loop_B[i] and A[i]==0 and B[i]==0:
            # offset_loop_table[i] = loop_A[i] - loop_B[i]
            # if loop_A[i] == loop_B[i] and current_ptable[i] == 0:
            #     compatible[i] = True
            # else:
            #     compatible[i] = False

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
                    print(self.fc.eval_structure_pt(new_ptable)/100)
                return self.fc.eval_structure_pt(new_ptable)/100

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
            # cmd = f'printf "{sequence}\n{C_str}\n" | RNAsubopt -C ??????enforceConstraint -e 2'
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
                # self.moves_add.remove((i,j))
                yield i, j
            # if i<0 and ls[i]==ls[j]:
                # print ("try", i, j)
                # self.moves_add.remove((i,j))
                # yield i,j

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
            "-inf"), current_e=evals[s1], moves=[(0, 0, s1_en)], energies=[0], opt=0, add_moves=[], distance=current_bp_end)
        # initial path start with 1 intermediate
        init_path = [init_intermediate]
        # paths start with 1 initial path
        paths = [init_path]


        # Ivos variation: define the largest distance, 
        # next iteration will only ever iterate over the longest distances

        longest_distance = current_bp_end

        queue = []

        # dont stop at current_bp_end, consider potential indirect moves
        while (len(paths)!=0):

            # copy all paths from the queue into 'paths' which have the longest distance


            # iterate over all queue items for longest distance first

            longest_distance = 0
            for p in queue:
                if p[-1].distance > longest_distance:
                    longest_distance = p[-1].distance
            for p in paths:
                if p[-1].distance > longest_distance:
                    longest_distance = p[-1].distance

            print ('longest d', longest_distance, 'paths:', len(paths), 'queue:', len(queue))

            # all elements in paths, which dont have the longest distance, get onto the queue.
            # from the queue, all items with the longest distance get to paths

            if longest_distance==14:
                print('before')
                for p in paths:
                    print(p[-1].distance)

            for p in queue:
                if p[-1].distance == longest_distance:
                    # print ('add back')
                    paths.append(p)
                    queue.remove(p)

            for p in paths:
                if p[-1].distance != longest_distance:
                    # print ('add back')
                    queue.append(p)
                    paths.remove(p)
                    if longest_distance==14:
                        print('remove', len(paths))
                else:
                    if longest_distance==14:
                        print('stay', p[-1].distance, len(paths))
                    
            
            if longest_distance==14:
                for i in range(len(paths)):

                    print (i, len(paths), paths[i][-1].distance)

            new_paths = [i for i in paths if i[-1].distance == longest_distance]
            



            if longest_distance==14:
                print('after')
                for p in new_paths:
                    print(p[-1].distance)

            test = []
            for p in paths:
                test.append(p[-1].distance)

                if longest_distance==14:
                    print(p[-1].distance)

            print ('longest d', longest_distance, 'paths:', len(paths), 'queue:', len(queue),
             'avg', np.mean(test), np.max(test), np.min(test))
         



            next_longest_distance = 0
            
            # collect all new paths here (next iteration)
            collect_paths = []

            for current_path in paths:

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
                               


            # paths = collect_paths
            current_bp += 1

            longest_distance = next_longest_distance
            # end while loop, get to the next distance class


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

    current_width = search_width  # for the time being, only 1 pass

    indirect_moves = [add_moves]
    best_max_en = float("inf")
    best_max_pos = 0
    best_path = False
    best_indirect_moves = []
    best_ptables = []
    best_indirect_move_count = 0

    for iteration in range(indirect_iterations):

        # first pass is direct, then iterate over different indirect moves

        next_indirect_moves = []

        for indirect_move in indirect_moves:

            if indirect_move == None:
                indirect_move = []  # first iteration only

            indirect_move = set(best_indirect_moves) | set(indirect_move)
            if Verbose: print("Iteration", iteration,
                  "launching findpath with addtional moves:", indirect_move)

            # add a set of optional indirect moves to the direct move set, e.g. add {(6, 12), (5, 13)}
            fp_class.moves_add = indirect_move
            # main findpath_once function call with supplied settings
            path_generator = fp_class.find_path_once(
                s1, s2, max_energy, current_width, mode=True, sort_min=False, Debug=Debug, Verbose=Verbose)

            # iterate over all currently found paths. somepaths might not be optimal
            for path in path_generator:
                # path is an intermediate struct, containing saddle_e, moves, etc.

                if path.saddle_e < max_energy:
                    max_energy = path.saddle_e

                e_0 = fp_class.evals[s1]
                e_1 = fp_class.evals[s2]
                current_ptable = list(fp_class.p_tables[s1])

                ptables = []  # list of ptables according to moves, should not be needed later

                # workaround convert moves to list of ptables, del this later
                current_used_indirect_moves = []  # which indirect moves are actually used
                current_moves = []
                for pos, (i, j, e) in enumerate(path.moves):
                    current_moves.append((i, j))
                    if (i, j) in indirect_move:
                        current_used_indirect_moves.append((i, j))
                    if i < 0:
                        current_ptable[-i] = 0
                        current_ptable[-j] = 0
                    else:
                        current_ptable[i] = j
                        current_ptable[j] = i
                    ptables.append(current_ptable.copy())

                # save best path for print output / next iteration
                if max_energy < best_max_en or \
                   (max_energy == best_max_en and len(current_used_indirect_moves) < best_indirect_move_count):
                   # or save best path with same energy, but fewer indirect moves: take this one instead

                    best_indirect_move_count = len(current_used_indirect_moves)
                    best_max_en = max_energy
                    best_path = path.moves.copy()
                    best_ptables = ptables.copy()
                    best_indirect_moves = current_used_indirect_moves.copy()
                    barrier = max_energy-e_0
                    if Verbose: print(
                        f"New best result: {max_energy:6.2f} kcal/mol | B: {barrier:6.2f} kcal/mol | E[start]:{e_0:6.2f} E[end]:{e_1:6.2f} | additional moves: {current_used_indirect_moves}")

                if iteration+1 != indirect_iterations:  # dont find new indirect moves in the last iteration

                    for current_indirect_moves in fp_class.find_stack(ptables, current_moves):
                        # print ("find move", current_indirect_moves, next_indirect_moves)
                        if current_indirect_moves not in next_indirect_moves:
                            next_indirect_moves.append(current_indirect_moves)

        # print path during last iteration
        if iteration+1 == indirect_iterations and Verbose:
            print_moves(sequence, s1, s2, best_path)
        # prepare for next indirect iteration
        else:
            # print ("found indirect moves:")
            # print (next_indirect_moves)
            indirect_moves = next_indirect_moves

        # best_path = [(0, 0), (-31, -45), (-35, -41), (-34, -42), (-33, -43),
        #  (-32, -44), (-5, -11), (-4, -12), (6, 12), (7, 11),
        #  (-3, -13), (5, 13), (-2, -14), (-1, -15), (4, 42), (3, 43), (2, 44), (1, 45),
        #  (-5, -13),(5, 41), (-6, -12), (6, 40), (-7, -11), (7, 39),
        #  (8, 38), (9, 37), (10, 36)]


    # post processing - check every indirect move - 
    # confirm if they're actually useful (linear time)

    best_path_indirect_moves = []
    # check which indirect moves were actually used
    filtered = [(i[0], i[1]) for i in best_path]
    for i, j in filtered:
        if (-i, -j) in filtered and (-i, -j) not in best_path_indirect_moves and (i, j) != (0, 0):
            best_path_indirect_moves.append((i, j))    

    indirect_moves = best_path_indirect_moves.copy()
    for i,j in best_path_indirect_moves:
        current_test = []
        for m,n,e in best_path:
            # if (-i, -j) in filtered: continue
            if (m, n) in filtered and (-m, -n) in filtered: continue
            current_test.append((m,n))
        test_en = print_moves(sequence, s1, s2, current_test, Verbose=False)
        if test_en == best_max_en:            
            filtered.remove((i,j))
            filtered.remove((-i,-j))
            indirect_moves.remove((i,j))


    # code here to not just return nothing
    return best_max_en, filtered, indirect_moves


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
    search_width = 30
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

    add_moves = {(9, 15)}

    # add_moves = []
    # indirect_iterations = 2
    indirect_iterations = 1
    paths = find_path(sequence, s1, s2, indirect_iterations=indirect_iterations, add_moves=add_moves,
                      search_width=search_width, Debug=Debug, Verbose=Verbose)


    


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
