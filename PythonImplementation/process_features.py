#!/usr/bin/env python3
# coding: utf-8

# process features

import sys
import time

import numpy as np
import matplotlib.pyplot as plt
import RNA

from sklearn import preprocessing
min_max_scaler = preprocessing.MinMaxScaler()

sys.path.append('../')
from pretty_print_path import print_moves
import findpath_librna
import findpath


from features import ij_distance, new_move_dist, plt_moves, config_distance, balance_in_all_things, return_shift_moves

def find_moves(s_ptable, t_ptable):
    """
    generator function, yields possible structures 1 move away
    from the original structure by finding fitting i and j with
    RNA pair and loop tables
    s_ptable: current ptable
    t_ptable: s2 end ptable
    """
    # loop table
    ls = RNA.loopidx_from_ptable(s_ptable)

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


def fp_call(sequence, s1, s2, search_width_multiplier = 20):    
    fp = findpath.findpath_single(sequence, s1, s2, search_width_multiplier=search_width_multiplier, mp=True)
    result = fp.get_en()/100.0
    path = fp.get_path()
    # return result, path
    return result, path


def process(sequence, s1, s2, search_width_multiplier = 20):

    s = s1
    lasts = s
    lasti = None
    lastj = None
    pt2 = list(RNA.ptable(s2))
    fc = RNA.fold_compound(sequence)

    result, path = fp_call(sequence, s1, s2, search_width_multiplier)
    all_moves = [(abs(i[0]), abs(i[1])) for i in path]

    prev_move_dist = np.zeros((20))   

    # print ("best result", result, "swm:", search_width_multiplier)
    # print_moves(sequence, s1, s2, path)

    results = []
    vec_results = []

    for e, (a,b, en) in enumerate(path):
        if (a,b) == (0,0):
            continue  

        

        # check where we can go, compare with our best move. 
        pt = list(RNA.ptable(s))

        # check available moves, save them, sort them    
        avail_moves = []
        ij_moves = []
        found_pos = None

        for pos, (i,j) in enumerate(find_moves(pt, pt2)):    
            next_en = fc.eval_move_pt(pt, i, j)
            # mark where we found our move
            found = (i,j) == (a,b)
            avail_moves.append((i, j, next_en, found))
            ij_moves.append((abs(i),abs(j)))

        # sort moves independent of delete insert moves

        ij_moves.sort(key=lambda x: (x[0], x[1]))


        i_shift_moves, j_shift_moves = return_shift_moves(all_moves)
        # print (i_shift_moves, ij_moves)

        avail_moves.sort(key=lambda x: x[2])
        found_list = [x[3] for x in avail_moves]
        en_list = np.array([[x[2] for x in avail_moves]])
        en_list_scaled = min_max_scaler.fit_transform(en_list.T).T[0]

        

        # find where our move is after sorting
        found_pos = found_list.index(True)
        rel_pos = found_pos * 1.0 / len(found_list)

        # print (e, a,b, 'found at pos:', found_pos, 'of', len(avail_moves), ':',  1-rel_pos)
       
        # for every move we take we have to run a new findpath, see if this move will yield the ideal result
        
        t1 = 0

        for pos, (i,j, en, found) in enumerate(avail_moves):
            
            t0 = time.time()
            
            if i > 0:
                snew = s[:i-1] + "(" + s[i:j-1] + ")" + s[j:]
            if i < 0:
                snew = s[:-i-1] + "." + s[-i:-j-1] + "." + s[-j:]
            ptnew = list(RNA.ptable(snew))

            result_new, path = fp_call(sequence, snew, s2, search_width_multiplier)

            if result_new <= result:
                pos_result = 1
            else:
                pos_result = 0

            if found: found = "<-- taken"
            else: found = ""

            this_move = (abs(i), abs(j))
            last_move = (lasti, lastj)

            t1 += time.time()- t0

            if lasti:
                # print (this_move, last_move, ij_moves)
                ijd, thisclose, lastclose = ij_distance(last_move, this_move, ij_moves)
            else:
                ijd, thisclose, lastclose = 0, 0, 0

            cd = config_distance(pt, this_move)
            
            if abs(i) in i_shift_moves:
                i_shift = 1
            else:
                i_shift = 0
            if abs(j) in j_shift_moves:
                j_shift = 1
            else:
                j_shift = 0

            if i<0:
                insert_or_delete = 0
            else:
                insert_or_delete = 1

            balance = balance_in_all_things(s1, s2, snew) / len(all_moves)

            # move_dist = new_move_dist(abs(i), abs(j), prev_move_dist)  


            # ignore results if we dont have past results
            # ignore if we're not making meaningful decisions
            if lasti and len(avail_moves) > 3:
                results.append((snew, i, j, pos_result, ijd, thisclose, lastclose, cd, en_list_scaled[pos],\
                i_shift, j_shift, insert_or_delete, balance))                
                vec_results.append(prev_move_dist.copy())
                
                # print (f'{snew}, {i}, {j}, {result_new}/{pos_result}: {ijd:2.2f},\
                #          {thisclose:2.2f}, {lastclose:2.2f}, {cd}, {en}, {en_list_scaled[pos]:2.2f} {found}')
                # print (move_dist)


        # if e==63:
        #     print (avail_moves)

        # update s for the next iteration
        


        lasts = s

        lasti = abs(a)
        lastj = abs(b)

        prev_move_dist = new_move_dist(lasti, lastj, prev_move_dist)

        if a > 0:
            s = s[:a-1] + "(" + s[a:b-1] + ")" + s[b:]
        if a < 0:
            s = s[:-a-1] + "." + s[-a:-b-1] + "." + s[-b:]
        
        # print (t1, lasti, lastj, prev_move_dist)

        # if e>1:
        #     break
    
    return results, vec_results

        # if not en:
        #     en = round(fc.eval_structure(s), 2)


        # break