#!/usr/bin/env python3
# coding: utf-8

import RNA

def print_moves(sequence, s1, s2, moves, move_color='\033[93m', Verbose = True, exclude=None, include=None, convert_to_float=False):

    """
    print a folding path with colour coding
    
    moves have to contain i,j or optionally i,j,en
    e.g. [(0, 0, -15.4), (-2, -24, -13.2), (-3, -23, -11.6)]

    without verbosity, this just returns max_en
    """

    # print (moves)

    # from stackexchange...
    class c:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        CYAN = '\033[96m'
        GREEN = '\033[92m'
        YELLOW = '\033[93m'
        RED = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'

    if Verbose: print(sequence)
    s = s1

    fc = RNA.fold_compound(sequence)
    e1 = en = round(fc.eval_structure(s), 2)
    max_en = float("-inf")

    output_rows = []
    moves_i_j = [(x[0],x[1]) for x in moves]

    # preprocessing - generate strings & energies if required
    for move in moves:
        en = False
        if len(move) == 2:
            i, j = move
        if len(move) == 3:
            i, j, en = move
        if i > 0:
            s = s[:i-1] + "(" + s[i:j-1] + ")" + s[j:]
        if i < 0:
            s = s[:-i-1] + "." + s[-i:-j-1] + "." + s[-j:]
        if not en:
            en = round(fc.eval_structure(s), 2)
        
        if convert_to_float:
            en = en / 100.0
        
        e2 = en
        if en > max_en:
            max_en = en

        output_rows.append((s, i, j, en))

    for s, i, j, en in output_rows:

        # print initial row with move (0,0)
        if i == 0:
            info = f'{move_color}[{i:4}, {j:4} ]{c.ENDC} {en:6.2f}'
            if Verbose: print(f"{s} {info}")
            continue

        pos_i = abs(i)-1
        pos_j = abs(j)-1

        # if a move has an identical inverse copy (1,2) <-> (-1,-2)
        # it is automatically an indirect move - these are colored in red 
        if (-i, -j) in moves_i_j: # indirect move
            colored_s = s[0:pos_i] + c.RED + c.BOLD + s[pos_i] + c.ENDC +\
                s[pos_i+1:pos_j] + c.RED + c.BOLD + \
                s[pos_j] + c.ENDC + s[pos_j+1:]
            info = f'{c.RED}[{i:4}, {j:4} ]{c.ENDC}'
        else:  # direct move
            colored_s = s[0:pos_i] + move_color + c.BOLD + s[pos_i] + c.ENDC +\
                s[pos_i+1:pos_j] + move_color + c.BOLD + \
                s[pos_j] + c.ENDC + s[pos_j+1:]
            info = f'{move_color}[{i:4}, {j:4} ]{c.ENDC}'

        # mark pos x
        # x = 10
        # colored_s = colored_s[0:x] + c.CYAN + colored_s[x:x+1] + c.ENDC + colored_s[x+1:]

        if en == max_en:
            info += f' {c.RED}{c.BOLD}{en:6.2f}{c.ENDC}'
        else:
            info += f' {en:6.2f}'

        if Verbose:
            if include != None:
                if abs(i) in include:
                    print(f"{info}")
            elif exclude != None:
                if abs(i) not in exclude:
                    print(f"{info}")
            else:
                print(f"{colored_s} {info}")

    barrier = max_en - e1
    if Verbose: print(
        f"S: {max_en:6.2f} kcal/mol | B: {barrier:6.2f} kcal/mol | E[start]:{e1:6.2f} E[end]:{e2:6.2f}")
    
    return max_en