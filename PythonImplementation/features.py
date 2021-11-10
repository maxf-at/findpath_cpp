#!/usr/bin/env python3
# coding: utf-8

# contains all sorts of feature engineering functions

from collections import Counter
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import RNA

def ij_distance(last_move, this_move, ij_moves):
    # how far is the last move away from the current move.
    # it is likely, that the next move is close to the last one
    # there are better distance metrices probably...

    # ij move list is supposed to be sorted, find indices
    
    ijmoves = ij_moves + [last_move]
    ijmoves.sort(key=lambda x: (x[0], x[1]))

    pos_old = ijmoves.index(last_move)
    pos_new = ijmoves.index(this_move)

    distance = abs(pos_old-pos_new)/len(ijmoves)
        
    # moves left in vicinity out of total moves
    thisi, thisj = this_move
    lasti, lastj = last_move
    thisclose = 0
    lastclose = 0

    for i, j in ij_moves:
        if (abs(i-thisi) < 5) and (abs(j-thisj) < 5):
            thisclose += 1
        if (abs(i-lasti) < 5) and (abs(j-lastj) < 5):
            lastclose += 1
    
    thisclose /= len(ij_moves)
    lastclose /= len(ij_moves)

    # print ("thisclose", thisclose, lastclose, distance)
    return distance, thisclose, lastclose

    print (distance)


# distribution of moves

def pdf(x, target):
    mean = target
    std = np.std(x)/8.0
    y_out = 1/(std * np.sqrt(2 * np.pi)) * np.exp( - (x - mean)**2 / (2 * std**2))
    return y_out

def calc_move_dist(i, j):
    # overlaps 2 normal distributions for our 2 move positions
    x = np.arange(-50, 50, 5)
    y = pdf(x, i-50)*4
    y += pdf(x, j-50)*4
    return y

def new_move_dist(i, j, prev):
    if i==0:
        return prev
    if np.mean(prev) != 0:
        prev  = (prev  - prev.min()) / (prev.max() - prev.min()) * 0.25
    new = calc_move_dist(i, j)
    return new + prev

def plt_moves(a):
    # Plotting the bell-shaped curve
    y = a
    x = np.arange(1,100,5)
    plt.style.use('seaborn')
    plt.figure(figsize = (6, 6))
    plt.plot(x, y, color = 'black',
            linestyle = 'dashed')    
    plt.scatter( x, y, marker = 'o', s = 25, color = 'red')
    plt.show()

# still needed?

def config_distance(pt, move):
    """
    are we extending / removing the outside/inside layer of a loop or adding something in the middle?
    """
    i = move[0]
    j = move[1]
    points = 0

    # if we're extending from outside to inside, the position i+1 and j-1 should be ideally unpaired
    # inside to outside: i-1 and j+1 should be ideally unpaired

    if i>0:
        # print ("add") 
        # outside/inside paired?        
        if j+1 < pt[0] and i-1 > 0: # outside - boundary check
            if pt[i-1] == j+1:
                points += 1
        if pt[i+1] == j-1:
            points += 1
    if i<0:
        # print ("del")
        i, j = -i, -j
        # outside/inside paired?
        if j+1 < pt[0] and i-1 > 0: # outside - boundary check
            if pt[i-1] == j+1:
                points += 1
        if pt[i+1] == j-1:
            points += 1
        elif pt[i+1] == 0 and pt[j-1] == 0:
            pass

    if points == 2:
        points = 0
    return points
    print (pt[i], pt[j])


def balance_in_all_things(s1, s2, s):
    """
    lets compare the current structure with s1 and s2. It can be assumed that a refolding
    path will not first delete all pairs and then add new ones, the process goes hand in hand:
    there's always a balance between paired and unpaired base pairs.
    """

    paired_s1 = len([len(i) for i in s1 if i=='(']) * 2
    paired_s2 = len([len(i) for i in s2 if i=='(']) * 2

    paired_s = len([len(i) for i in s if i=='(']) * 2

    s1_dist = abs(paired_s1 - paired_s)
    s2_dist = abs(paired_s2 - paired_s)

    return min(s1_dist, s2_dist)


def return_shift_moves(ij_moves):

    # ij_moves = [(abs(i[0]), abs(i[1])) for i in path]
    # print (ij_moves)
    i_dict = defaultdict(list)
    j_dict = defaultdict(list)
    i_list = [i for i in i_dict.keys()]

    # def shift_move():
    for i, j in ij_moves:
        i_dict[i].append(j)
        j_dict[j].append(i)

    i_shift = [i for i in list(i_dict) if len(i_dict[i]) > 1]
    j_shift = [i for i in list(j_dict) if len(j_dict[i]) > 1]    
    return i_shift, j_shift

