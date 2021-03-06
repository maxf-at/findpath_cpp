/* gcc -fopenmp -g3 -DTEST_FINDPATH findpath.c -o FINDpath -lRNA -lm -I../ -L./

        g++ findpath.cpp -o3 -std=c++17 -o findpath -lRNA -lm -DTEST_FINDPATH -DEBUG_FINDPATH

 */

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <algorithm>

#include <future>
#include <thread>

// #include <execution>

#include "./icecream.hpp"
using icecream::ic;

using namespace std;

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "findpath_B.hpp"

// extern "C" {
// #include <limits.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>

// #include "ViennaRNA/cofold.h"
// #include "ViennaRNA/datastructures/basic.h"
// #include "ViennaRNA/fold.h"
// #include "ViennaRNA/fold_vars.h"
// #include "ViennaRNA/landscape/findpath.h"
// #include "ViennaRNA/model.h"
// #include "ViennaRNA/params/basic.h"
// #include "ViennaRNA/utils/basic.h"
// #include "ViennaRNA/utils/strings.h"
// #include "ViennaRNA/utils/structures.h"

// }

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#ifdef _OPENMP
#include <omp.h>
#endif

#endif

#define LOOP_EN

#define PATH_DIRECT_FINDPATH 1U

/**
 *  @brief
 */
typedef struct move {
    int i; /* i,j>0 insert; i,j<0 delete */
    int j;
    int when; /* 0 if still available, else resulting distance from start */
    int E;
} move_t;

/**
 *  @brief
 */
typedef struct intermediate {
    short*  pt;      /**<  @brief  pair table */
    int     Sen;     /**<  @brief  saddle energy so far */
    int     curr_en; /**<  @brief  current energy */
    move_t* moves;   /**<  @brief  remaining moves to target */
} intermediate_t;

struct vrna_path_options_s {
    unsigned int type;
    unsigned int method;

    int width;
};

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

// PRIVATE int BP_dist;
// PRIVATE move_t* path = NULL;
// PRIVATE int     path_fwd; /* 1: s1->s2, else s2 -> s1 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PRIVATE vrna_fold_compound_t* backward_compat_compound = NULL;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
 */
// #pragma omp threadprivate(BP_dist, path, path_fwd, backward_compat_compound)
#pragma omp threadprivate(backward_compat_compound)

#endif

#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE move_t* copy_moves(move_t* mvs, int bp_dist);

PRIVATE int compare_ptable(const void* A, const void* B);

PRIVATE int compare_energy(const void* A, const void* B);

PRIVATE int compare_moves_when(const void* A, const void* B);

PRIVATE void free_intermediate(intermediate_t* i);

#ifdef TEST_FINDPATH

/* TEST_FINDPATH, COFOLD */
PRIVATE void usage(void);

#endif

auto find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2, int current_search_width,
                    int maxE, bool direction, int search_width) -> std::vector<sorted_path>;

PRIVATE int try_moves(vrna_fold_compound_t* vc, intermediate_t c, int maxE, intermediate_t* next,
                      int dist, int bp_dist);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void vrna_path_free(vrna_path_t* path)
{
    vrna_path_t* tmp = path;

    if (tmp) {
        if (tmp->type == VRNA_PATH_TYPE_DOT_BRACKET) {
            while (tmp->s) {
                free(tmp->s);
                tmp++;
            }
        } else if (tmp->type == VRNA_PATH_TYPE_MOVES) {
            while ((tmp->move.pos_5)) {
                vrna_move_list_free(tmp->move.next);
                tmp++;
            }
        }

        free(path);
    }
}

PUBLIC int vrna_path_findpath_saddle(vrna_fold_compound_t* vc, const char* s1, const char* s2,
                                     int width)
{
    std::vector<sorted_path> obj;

    short* pt1 = vrna_ptable(s1);
    short* pt2 = vrna_ptable(s2);

    obj = custom_findpath_method(vc, pt1, pt2, width, INT_MAX - 1, false);

    return obj[0].max_en;
}

// auto custom_findpath_method(vrna_fold_compound_t* vc, const char* s1, const char* s2, int width,
//                             int maxE) -> std::vector<sorted_path>

// auto custom_findpath_method(vrna_fold_compound_t* vc, short* pt1, short* pt2, int width, int
// maxE)
//     -> std::vector<sorted_path>

// // ub = upper bound?
// {
//     int current_search_width;
//     // short * ptr, *pt1, *pt2;
//     short* ptr;

//     move_t* bestpath = NULL;
//     int     dir;

//     std::vector<sorted_path> all_paths;

//     path_fwd = dir = 0;

//     // new private bool
//     bool direction = false;  // 1: s1->s2, else s2 -> s1

//     // pt1 = vrna_ptable(s1);
//     // pt2 = vrna_ptable(s2);

//     current_search_width = 16;
//     // current_search_width = 1;
//     // do while: at least executed once.
//     do {
//         int saddleE;
//         path_fwd  = !path_fwd;
//         direction = !direction;

//         if (current_search_width > width) current_search_width = width;

//         if (path) free(path);

//         // the compiler should take care of copy elision here...
//         std::vector<sorted_path> current_paths;
//         current_paths = find_path_once(vc, pt1, pt2, current_search_width, maxE, direction,
//         width);

//         // did we actually receive a path / paths?
//         if (current_paths.size() > 0) {
//             saddleE = current_paths[0].max_en;

//             // for (auto &move:current_paths[0].moves) {
//             //     IC(move.i, move.j);

//             // }
//             // IC(direction, current_search_width, current_paths.size(), saddleE);

//         } else {
//             saddleE = maxE;
//         }

//         if (saddleE < maxE) {
//             maxE = saddleE;
//             if (bestpath) free(bestpath);

//             bestpath = path;
//             path     = NULL;
//             dir      = path_fwd;
//         } else {
//             free(path);
//             path = NULL;
//         }

//         ptr = pt1;
//         pt1 = pt2;
//         pt2 = ptr;
//         current_search_width *= 2;

//         // concatenate vectors
//         std::move(current_paths.begin(), current_paths.end(), std::back_inserter(all_paths));

//     } while (current_search_width < 2 * width);

//     IC(all_paths.size(), current_search_width, all_paths[0].max_en);

//     // current workaround...
//     all_paths[0].max_en = maxE;

//     /* (re)set some globals */
//     path     = bestpath;
//     path_fwd = dir;

//     free(pt1);
//     free(pt2);

//     return all_paths;
//     // return maxE;
// }

auto custom_findpath_method(vrna_fold_compound_t* vc, short* pt1, short* pt2, int width, int maxE,
                            bool init_direction) -> std::vector<sorted_path>

// ub = upper bound?
{
    // short * ptr, *pt1, *pt2;
    short* ptr;

    move_t* bestpath = NULL;
    int     dir;

    std::vector<sorted_path> all_paths;

    bool path_fwd = dir = 0;

    // new private bool

    bool direction = false;  // true: s1 -> s2, false: s2 -> s1

    if (not init_direction) {
        direction = true;
        ptr       = pt1;
        pt1       = pt2;
        pt2       = ptr;
    }

    // pt1 = vrna_ptable(s1);
    // pt2 = vrna_ptable(s2);

    int              last_iteration = width;
    std::vector<int> iterations{width};

    float accelerator = 1.0;

    while (last_iteration > 16) {
        const int next_iteration = int(last_iteration / 5.0 * accelerator);
        iterations.push_back(next_iteration);
        last_iteration = next_iteration;

        accelerator *= accelerator;

        if (accelerator > 2) accelerator = 2;
    }

    IC(iterations);

    int current_search_width;
    // current_search_width = 1;
    // do while: at least executed once.

    for (vector<int>::reverse_iterator it = iterations.rbegin(); it != iterations.rend(); ++it) {
        current_search_width = *it;
        // IC(current_search_width);

        int saddleE;
        path_fwd  = !path_fwd;
        direction = !direction;

        // overshoot
        // if (current_search_width > width) current_search_width = width;

        // if (path) free(path);

        // the compiler should take care of copy elision here...

        std::vector<sorted_path> current_paths;
        current_paths = find_path_once(vc, pt1, pt2, current_search_width, maxE, direction, width);

        // std::future<std::vector<sorted_path>> ret1 =
        //     std::async(std::launch::async, &find_path_once, vc, pt1, pt2, current_search_width,
        //                maxE, direction, width);
        // std::future<std::vector<sorted_path>> ret2 =
        //     std::async(std::launch::async, &find_path_once, vc, pt1, pt2, current_search_width,
        //                maxE, not direction, width);

        // std::vector<sorted_path> current_paths  = ret1.get();

        // std::vector<sorted_path> current_paths2 = ret2.get();

        // did we actually receive a path / paths?

        // #ifdef IC_DEBUG
        // IC("r", direction, current_search_width, current_paths.size(), saddleE);
        // #endif

        if (current_paths.size() > 0) {
            saddleE = current_paths[0].max_en;

            // for (auto &move:current_paths[0].moves) {
            //     IC(move.i, move.j);

            // }
            // IC(direction, current_search_width, current_paths.size(), saddleE);

        } else {
            saddleE = maxE;
        }

        if (saddleE < maxE) {
            maxE = saddleE;
            // if (bestpath) free(bestpath);

            // bestpath = path;
            // path     = NULL;
            // dir      = path_fwd;
            // } else {
            //     free(path);
            //     path = NULL;
        }

        ptr = pt1;
        pt1 = pt2;
        pt2 = ptr;
        // current_search_width *= 4;

        // concatenate vectors
        std::move(current_paths.begin(), current_paths.end(), std::back_inserter(all_paths));
    }

    // while (current_search_width < 4 * width);

    // #ifdef IC_DEBUG
    IC(all_paths.size(), current_search_width, all_paths[all_paths.size() - 1].max_en);
    // #endif

    // current workaround...
    // all_paths[0].max_en = maxE;

    /* (re)set some globals */
    // path     = bestpath;
    // path_fwd = dir;

    free(pt1);
    free(pt2);

    return all_paths;
    // return maxE;
}

PUBLIC int vrna_path_findpath_saddle_ub(vrna_fold_compound_t* vc, const char* s1, const char* s2,
                                        int width, int maxE)
{  // dummy
    return maxE;
}

PUBLIC vrna_path_t* vrna_path_findpath(vrna_fold_compound_t* fc, const char* s1, const char* s2,
                                       int width)
{
    struct vrna_path_options_s* opt;
    vrna_path_t*                path;

    opt  = vrna_path_options_findpath(width, VRNA_PATH_TYPE_DOT_BRACKET);
    path = vrna_path_direct_ub(fc, s1, s2, INT_MAX - 1, opt);

    free(opt);

    return path;
}

PUBLIC vrna_path_t* vrna_path_findpath_ub(vrna_fold_compound_t* fc, const char* s1, const char* s2,
                                          int width, int maxE)
{
    struct vrna_path_options_s* opt;
    vrna_path_t*                path;

    opt  = vrna_path_options_findpath(width, VRNA_PATH_TYPE_DOT_BRACKET);
    path = vrna_path_direct_ub(fc, s1, s2, maxE, opt);

    free(opt);

    return path;
}

PUBLIC struct vrna_path_options_s* vrna_path_options_findpath(int width, unsigned int type)
{
    struct vrna_path_options_s* options =
        (struct vrna_path_options_s*)vrna_alloc(sizeof(struct vrna_path_options_s));

    options->type   = type;
    options->method = PATH_DIRECT_FINDPATH;
    options->width  = width;

    return options;
}

PUBLIC void vrna_path_options_free(struct vrna_path_options_s* options) { free(options); }

PRIVATE vrna_path_t* findpath_method(vrna_fold_compound_t* fc, const char* s1, const char* s2,
                                     int width, int maxE, unsigned int return_type)
{
    // this does not work any longer... delete
    move_t* path     = NULL;
    int     BP_dist  = 0;
    int     path_fwd = 0;

    int          E, d;
    float        last_E;
    vrna_path_t* route = NULL;

    E = vrna_path_findpath_saddle_ub(fc, s1, s2, width, maxE);

    /* did we find a better path than one with saddle maxE? */
    if (E < maxE) {
        route = (vrna_path_t*)vrna_alloc((BP_dist + 2) * sizeof(vrna_path_t));

        qsort(path, BP_dist, sizeof(move_t), compare_moves_when);

        switch (return_type) {
            case VRNA_PATH_TYPE_MOVES:

                // #ifdef IC_DEBUG
                IC("return type:", return_type);
                // #endif

                if (path_fwd) {
                    last_E = vrna_eval_structure(fc, s1);
                    for (d = 0; d < BP_dist; d++) {
                        route[d].type = return_type;
                        route[d].move = vrna_move_init(path[d].i, path[d].j);
                        route[d].en   = (path[d].E / 100.0) - last_E;
                        last_E        = path[d].E / 100.0;
                    }

                    route[BP_dist].type = return_type;
                    route[BP_dist].move = vrna_move_init(0, 0);
                } else {
                    last_E = vrna_eval_structure(fc, s2);
                    for (d = 0; d < BP_dist; d++) {
                        route[BP_dist - d - 2].type = return_type;
                        route[BP_dist - d - 2].move = vrna_move_init(path[d].i, path[d].j);
                        route[BP_dist - d - 2].en   = last_E - (path[d].E / 100.0);
                        last_E                      = path[d].E / 100;
                    }

                    route[BP_dist].type = return_type;
                    route[BP_dist].move = vrna_move_init(0, 0);
                }

                break;

            case VRNA_PATH_TYPE_DOT_BRACKET:
                /* fall through */

            default:
                if (path_fwd) {
                    /* memorize start of path */
                    route[0].type = return_type;
                    route[0].s    = strdup(s1);
                    route[0].en   = vrna_eval_structure(fc, s1);

                    for (d = 0; d < BP_dist; d++) {
                        int i, j;
                        route[d + 1].type = return_type;
                        route[d + 1].s    = strdup(route[d].s);
                        i                 = path[d].i;
                        j                 = path[d].j;
                        if (i < 0) {
                            /* delete */
                            route[d + 1].s[(-i) - 1] = route[d + 1].s[(-j) - 1] = '.';
                        } else {
                            route[d + 1].s[i - 1] = '(';
                            route[d + 1].s[j - 1] = ')';
                        }

                        route[d + 1].en = path[d].E / 100.0;
                    }
                } else {
                    /* memorize start of path */
                    route[0].type     = return_type;
                    route[BP_dist].s  = strdup(s2);
                    route[BP_dist].en = vrna_eval_structure(fc, s2);

                    for (d = 0; d < BP_dist; d++) {
                        int i, j;
                        route[BP_dist - d - 1].type = return_type;
                        route[BP_dist - d - 1].s    = strdup(route[BP_dist - d].s);
                        i                           = path[d].i;
                        j                           = path[d].j;
                        if (i < 0) {
                            /* delete */
                            route[BP_dist - d - 1].s[(-i) - 1] =
                                route[BP_dist - d - 1].s[(-j) - 1] = '.';
                        } else {
                            route[BP_dist - d - 1].s[i - 1] = '(';
                            route[BP_dist - d - 1].s[j - 1] = ')';
                        }

                        route[BP_dist - d - 1].en = path[d].E / 100.0;
                    }
                }

                break;
        }

#if _DEBUG_FINDPATH_
        if (return_type & VRNA_PATH_TYPE_DOT_BRACKET) {
            fprintf(stderr, "\n%s\n%s\n%s\n\n", seq, s1, s2);
            for (d = 0; d <= BP_dist; d++) fprintf(stderr, "%s %6.2f\n", route[d].s, route[d].en);
            fprintf(stderr, "%d\n", *num_entry);
        }

#endif
    }

    free(path);
    path = NULL;

    return route;
}

PUBLIC vrna_path_t* vrna_path_direct(vrna_fold_compound_t* fc, const char* s1, const char* s2,
                                     struct vrna_path_options_s* options)
{
    return vrna_path_direct_ub(fc, s1, s2, INT_MAX - 1, options);
}

PUBLIC vrna_path_t* vrna_path_direct_ub(vrna_fold_compound_t* fc, const char* s1, const char* s2,
                                        int maxE, struct vrna_path_options_s* options)
{
    int                         E, d;
    struct vrna_path_options_s* o;
    vrna_path_t*                route = NULL;

    /* we default to findpath method */
    o = options ? options : vrna_path_options_findpath(10, VRNA_PATH_TYPE_DOT_BRACKET);

    switch (o->method) {
        case PATH_DIRECT_FINDPATH:
        /* fall through */
        default: route = findpath_method(fc, s1, s2, o->width, maxE, o->type); break;
    }

    if (!options) vrna_path_options_free(o);

    return route;
}

#ifdef TEST_FINDPATH

PUBLIC void print_path(const char* seq, const char* struc)
{
    int   d;
    char* s;

    s = strdup(struc);
    if (cut_point == -1) {
        printf("%s\n%s\n", seq, s);
    }
    /* printf("%s\n%s %6.2f\n", seq, s, vrna_eval_structure_simple(seq,s)); */
    else {
        char *pstruct, *pseq;
        pstruct = vrna_cut_point_insert(s, cut_point);
        pseq    = vrna_cut_point_insert(seq, cut_point);
        printf("%s\n%s\n", pseq, pstruct);
        /* printf("%s\n%s %6.2f\n", pseq, pstruct, vrna_eval_structure_simple(seq,s)); */
        free(pstruct);
        free(pseq);
    }

    qsort(path, BP_dist, sizeof(move_t), compare_moves_when);
    for (d = 0; d < BP_dist; d++) {
        int i, j;
        i = path[d].i;
        j = path[d].j;
        if (i < 0) {
            /* delete */
            s[(-i) - 1] = s[(-j) - 1] = '.';
        } else {
            s[i - 1] = '(';
            s[j - 1] = ')';
        }

        /* printf("%s %6.2f - %6.2f\n", s, vrna_eval_structure_simple(seq,s), path[d].E/100.0); */
    }
    free(s);
}

int main(int argc, char* argv[])
{
    int          buffer = 81;
    char         line[1000];
    char *       seq, *s1, *s2;
    int          E, maxkeep = 1000;
    int          verbose = 0, i;
    vrna_path_t *route, *r;

    for (i = 1; i < argc; i++) {
        switch (argv[i][1]) {
            case 'm':
                if (strcmp(argv[i], "-m") == 0) sscanf(argv[++i], "%d", &maxkeep);

                break;
            case 'v': verbose = !strcmp(argv[i], "-v"); break;
            case 'd':
                if (strcmp(argv[i], "-d") == 0) sscanf(argv[++i], "%d", &dangles);

                break;
            default: usage();
        }
    }

    cut_point = -1;

    // line      = vrna_read_line(stdin);

    fgets(line, 1000, stdin);
    strtok(line, "\n");
    seq = vrna_cut_point_remove(line, &cut_point);
    // free(line);
    // line  = vrna_read_line(stdin);
    fgets(line, 1000, stdin);
    strtok(line, "\n");
    s1 = vrna_cut_point_remove(line, &cut_point);
    // free(line);
    // line  = vrna_read_line(stdin);
    fgets(line, 1000, stdin);
    strtok(line, "\n");
    s2 = vrna_cut_point_remove(line, &cut_point);
    // free(line);

    // int                   maxE;
    // char*                 sequence;
    // vrna_fold_compound_t* vc;
    // vrna_md_t             md, *md_p;
    // vc = NULL;
    // set_model_details(&md);

    // if (backward_compat_compound) {
    //     if (!strcmp(seq, backward_compat_compound->sequence)) {
    //         /* check if sequence is the same as before */
    //         md.window_size = backward_compat_compound->length;
    //         md.max_bp_span = backward_compat_compound->length;
    //         md_p           = &(backward_compat_compound->params->model_details);
    //         if (!memcmp(&md, md_p, sizeof(vrna_md_t))) /* check if model_details are the same as
    //         before */
    //             vc = backward_compat_compound;         /* re-use previous vrna_fold_compound_t */
    //     }
    // }

    // if (!vc) {
    //     vrna_fold_compound_free(backward_compat_compound);

    //     sequence = vrna_cut_point_insert(seq, cut_point);

    //     backward_compat_compound = vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);

    //     free(sequence);
    // }
    // maxE = vrna_path_findpath_saddle(vc, s1, s2, width);

    vrna_fold_compound_t* fc = nullptr;
    vrna_md_t             md;
    set_model_details(&md);
    // md.window_size = backward_compat_compound->length;
    // md.max_bp_span = backward_compat_compound->length;
    fc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);

    E = vrna_path_findpath_saddle(fc, s1, s2, maxkeep);

    // E = find_saddle(seq, s1, s2, maxkeep);

    if (verbose) {
        if (path_fwd)
            print_path(seq, s1);
        else
            print_path(seq, s2);

        free(path);
        path  = NULL;
        route = get_path(seq, s1, s2, maxkeep);
        for (r = route; r->s; r++) {
            if (cut_point == -1) {
                printf("%s %6.2f\n", r->s, r->en);
                /* printf("%s %6.2f - %6.2f\n", r->s, vrna_eval_structure_simple(seq,r->s), r->en);
                 */
            } else {
                char* pstruct;
                pstruct = vrna_cut_point_insert(r->s, cut_point);
                printf("%s %6.2f\n", pstruct, r->en);
                /* printf("%s %6.2f - %6.2f\n", pstruct, vrna_eval_structure_simple(seq,r->s),
                 * r->en); */
                free(pstruct);
            }

            free(r->s);
        }
        free(route);
    }
    printf("%6.2f\n", E / 100.);
    free(seq);
    free(s1);
    free(s2);
    return EXIT_SUCCESS;
}

static void usage(void) { vrna_message_error("usage: findpath.c  [-m depth] [-d[0|1|2]] [-v]"); }

#endif

/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE int try_moves(vrna_fold_compound_t* vc, intermediate_t c, int maxE, intermediate_t* next,
                      int dist, int bp_dist)
{
    int *   loopidx, len, num_next = 0, en, oldE;
    move_t* mv;
    short*  pt;

    // len     = c.pt[0];
    // loopidx = vrna_loopidx_from_ptable(c.pt);
    // oldE    = c.Sen;
    // for (mv = c.moves; mv->i != 0; mv++) {
    //     int i, j;
    //     if (mv->when > 0) continue;

    //     i  = mv->i;
    //     j  = mv->j;
    //     pt = (short*)vrna_alloc(sizeof(short) * (len + 1));
    //     memcpy(pt, c.pt, (len + 1) * sizeof(short));
    //     if (j < 0) {
    //         /*it's a delete move */
    //         pt[-i] = 0;
    //         pt[-j] = 0;
    //     } else {
    //         /* insert move */
    //         if ((loopidx[i] == loopidx[j]) && /* i and j belong to same loop */
    //             (pt[i] == 0) && (pt[j] == 0)  /* ... and are unpaired */
    //         ) {
    //             pt[i] = j;
    //             pt[j] = i;
    //         } else {
    //             free(pt);
    //             continue; /* llegal move, try next; */
    //         }
    //     }

    // this does not work without LOOP_EN

    len     = c.pt[0];
    loopidx = vrna_loopidx_from_ptable(c.pt);
    oldE    = c.Sen;
    for (mv = c.moves; mv->i != 0; mv++) {
        int i, j;

        auto source_1 = 0;
        auto source_2 = 0;
        auto dest_1   = 0;
        auto dest_2   = 0;

        if (mv->when > 0) continue;

        i = mv->i;
        j = mv->j;
        // pt = (short*)vrna_alloc(sizeof(short) * (len + 1));
        // memcpy(pt, c.pt, (len + 1) * sizeof(short));
        if (j < 0) {
            /*it's a delete move */
            source_1 = -i;
            dest_1   = 0;
            source_2 = -j;
            dest_2   = 0;

        } else {
            /* insert move */
            if ((loopidx[i] == loopidx[j]) &&    /* i and j belong to same loop */
                (c.pt[i] == 0) && (c.pt[j] == 0) /* ... and are unpaired */
            ) {
                source_1 = i;
                dest_1   = j;
                source_2 = j;
                dest_2   = i;

            } else {
                // free(pt);
                continue; /* llegal move, try next; */
            }
        }

        // pt = (short*)vrna_alloc(sizeof(short) * (len + 1));
        // memcpy(pt, c.pt, (len + 1) * sizeof(short));

#ifdef LOOP_EN
        en = c.curr_en + vrna_eval_move_pt(vc, c.pt, i, j);
#else
        en = vrna_eval_structure_pt(vc, pt);
#endif
        // this used to be en < maxE
        if (en <= maxE) {
            pt = (short*)vrna_alloc(sizeof(short) * (len + 1));
            memcpy(pt, c.pt, (len + 1) * sizeof(short));

            pt[source_1] = dest_1;
            pt[source_2] = dest_2;

            // if (en < maxE) {
            next[num_next].Sen     = (en > oldE) ? en : oldE;
            next[num_next].curr_en = en;
            next[num_next].pt      = pt;
            mv->when               = dist;
            mv->E                  = en;
            next[num_next++].moves = copy_moves(c.moves, bp_dist);
            mv->when               = 0;
        }
        // else {
        //     free(pt);
        // }
    }
    free(loopidx);
    return num_next;
}

// removed global variales
// removed global direction, path

auto find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2, int current_search_width,
                    int maxE, bool direction, int final_search_width) -> std::vector<sorted_path>
{
    move_t*         mlist;
    int             i, len, d, dist = 0, result;
    short*          pt;
    intermediate_t *current, *next;  // array of current, array of next

    len = (int)pt1[0];
    pt  = vrna_ptable_copy(pt1);

    mlist = (move_t*)vrna_alloc(sizeof(move_t) * len); /* bp_dist < n */

    // generate all possible moves
    // this weirdly also calculates bp_dist...

    for (i = 1; i <= len; i++) {
        if (pt[i] != pt2[i]) {
            if (i < pt[i]) {
                /* need to delete this pair */
                mlist[dist].i      = -i;
                mlist[dist].j      = -pt[i];
                mlist[dist++].when = 0;
            }

            if (i < pt2[i]) {
                /* need to insert this pair */
                mlist[dist].i      = i;
                mlist[dist].j      = pt2[i];
                mlist[dist++].when = 0;
            }
        }
    }

    int bp_dist = dist;  // this is now local...

    // #ifdef IC_DEBUG
    IC(direction, dist, current_search_width, maxE, std::this_thread::get_id());
    // #endif

    current = (intermediate_t*)vrna_alloc(sizeof(intermediate_t) * (current_search_width + 1));
    current[0].pt  = pt;
    current[0].Sen = current[0].curr_en = vrna_eval_structure_pt(vc, pt);
    current[0].moves                    = mlist;
    next = (intermediate_t*)vrna_alloc(sizeof(intermediate_t) * (dist * current_search_width + 1));

    std::vector<sorted_path> all_paths;

    for (d = 1; d <= dist; d++) {
        /* go through the distance classes */
        int             c, u, num_next = 0;
        intermediate_t* cc;

        for (c = 0; current[c].pt != NULL; c++)
            // this fills up teh next array of intermediates
            num_next += try_moves(vc, current[c], maxE, next + num_next, d, bp_dist);
        if (num_next == 0) {
            // case where we don't find any moves -> abort
            for (cc = current; cc->pt != NULL; cc++) free_intermediate(cc);
            current[0].Sen = INT_MAX;
            break;
        }

        /* remove duplicates via sort|uniq
         * if this becomes a bottleneck we can use a hash instead */
        qsort(next, num_next, sizeof(intermediate_t), compare_ptable);

        // this shrinks the next array which makes the following sorting step faster
        for (u = 0, c = 1; c < num_next; c++) {
            if (memcmp(next[u].pt, next[c].pt, sizeof(short) * len) != 0) {
                next[++u] = next[c];

            } else {
                // new part - dont just delete valid duplicates
                // we dont get here during the greedy first iteration
                // only consider last fwd & bwd pass

                if (d >= dist && next[c].Sen <= maxE &&
                    current_search_width >= final_search_width) {
                    // if (next[c].Sen <= maxE && current_search_width>=search_width){
                    move_t* temp_moves = copy_moves(next[c].moves, bp_dist);
                    qsort(temp_moves, bp_dist, sizeof(move_t), compare_moves_when);

                    sorted_path current_path(dist);
                    // sorted_path_2 tests(dist);

                    // qsort(path, bp_dist, sizeof(move_t), compare_moves_when);
                    for (d = 0; d < bp_dist; d++) {
                        int i, j;
                        i = temp_moves[d].i;
                        j = temp_moves[d].j;

                        if (direction) {
                            // fwd path
                            current_path.moves[d] = {i, j, 0};
                        } else {
                            // bwd path, fill vector back to front
                            current_path.moves[dist - d - 1] = {-i, -j, 0};
                        }
                    }
                    current_path.max_en = next[c].Sen;
                    all_paths.push_back(current_path);
                }

                free_intermediate(next + c);
            }
        }

        // IC(num_next);

        num_next = u + 1;
        qsort(next, num_next, sizeof(intermediate_t), compare_energy);

        // next is now again reduced to size current, means replace current with next.
        /* free the old stuff */
        for (cc = current; cc->pt != NULL; cc++) free_intermediate(cc);
        for (u = 0; u < current_search_width && u < num_next; u++) current[u] = next[u];
        for (; u < num_next; u++) free_intermediate(next + u);
        num_next = 0;
    }
    free(next);

    // generate temporary path object...
    move_t* path = nullptr;

    path   = current[0].moves;
    result = current[0].Sen;
    free(current[0].pt);
    free(current);

    if (path) {
        // IC("path", direction);

        // std::vector<sorted_move> current_path;

        // preallocate bp_dist elements in vector
        sorted_path current_path(dist);
        // sorted_path_2 tests(dist);

        qsort(path, bp_dist, sizeof(move_t), compare_moves_when);
        for (d = 0; d < bp_dist; d++) {
            int i, j;
            i = path[d].i;
            j = path[d].j;

            if (direction) {
                // fwd path
                current_path.moves[d] = {i, j, 0};
            } else {
                // bwd path, fill vector back to front
                current_path.moves[dist - d - 1] = {-i, -j, 0};
            }
        }

        current_path.max_en = result;
        all_paths.push_back(current_path);

        // IC(all_paths[0].moves.size(), all_paths[0].moves[0].i);
        free(path);
    }

    // IC(all_paths.size(), result, current_search_width);
    return all_paths;
    // return result;
}

PRIVATE void free_intermediate(intermediate_t* i)
{
    free(i->pt);
    free(i->moves);
    i->pt    = NULL;
    i->moves = NULL;
    i->Sen   = INT_MAX;
}

PRIVATE int compare_ptable(const void* A, const void* B)
{
    intermediate_t *a, *b;
    int             c;

    a = (intermediate_t*)A;
    b = (intermediate_t*)B;

    c = memcmp(a->pt, b->pt, a->pt[0] * sizeof(short));
    if (c != 0) return c;

    if ((a->Sen - b->Sen) != 0) return a->Sen - b->Sen;

    return a->curr_en - b->curr_en;
}

PRIVATE int compare_energy(const void* A, const void* B)
{
    intermediate_t *a, *b;

    a = (intermediate_t*)A;
    b = (intermediate_t*)B;

    if ((a->Sen - b->Sen) != 0) return a->Sen - b->Sen;

    return a->curr_en - b->curr_en;
}

PRIVATE int compare_moves_when(const void* A, const void* B)
{
    move_t *a, *b;

    a = (move_t*)A;
    b = (move_t*)B;

    return a->when - b->when;
}

PRIVATE move_t* copy_moves(move_t* mvs, int bp_dist)
{
    move_t* new_2;

    new_2 = (move_t*)vrna_alloc(sizeof(move_t) * (bp_dist + 1));
    memcpy(new_2, mvs, sizeof(move_t) * (bp_dist + 1));
    return new_2;
}

/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC void free_path(vrna_path_t* path) { vrna_path_free(path); }

PUBLIC int find_saddle(const char* seq, const char* s1, const char* s2, int width)
{
    int                   maxE;
    char*                 sequence;
    vrna_fold_compound_t* vc;
    vrna_md_t             md, *md_p;

    vc = NULL;
    set_model_details(&md);

    if (backward_compat_compound) {
        if (!strcmp(seq, backward_compat_compound->sequence)) {
            /* check if sequence is the same as before */
            md.window_size = backward_compat_compound->length;
            md.max_bp_span = backward_compat_compound->length;
            md_p           = &(backward_compat_compound->params->model_details);
            if (!memcmp(&md, md_p,
                        sizeof(vrna_md_t)))    /* check if model_details are the same as before */
                vc = backward_compat_compound; /* re-use previous vrna_fold_compound_t */
        }
    }

    if (!vc) {
        vrna_fold_compound_free(backward_compat_compound);

        sequence = vrna_cut_point_insert(seq, cut_point);

        backward_compat_compound = vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);

        free(sequence);
    }

    maxE = vrna_path_findpath_saddle(vc, s1, s2, width);

    return maxE;
}

PUBLIC vrna_path_t* get_path(const char* seq, const char* s1, const char* s2, int maxkeep)
{
    vrna_path_t*          route    = NULL;
    char*                 sequence = NULL;
    vrna_fold_compound_t* vc       = NULL;
    vrna_md_t             md, *md_p;

    set_model_details(&md);

    if (backward_compat_compound) {
        if (!strcmp(seq, backward_compat_compound->sequence)) {
            /* check if sequence is the same as before */
            md.window_size = backward_compat_compound->length;
            md.max_bp_span = backward_compat_compound->length;
            md_p           = &(backward_compat_compound->params->model_details);
            if (!memcmp(&md, md_p,
                        sizeof(vrna_md_t)))    /* check if model_details are the same as before */
                vc = backward_compat_compound; /* re-use previous vrna_fold_compound_t */
        }
    }

    if (!vc) {
        vrna_fold_compound_free(backward_compat_compound);

        sequence = vrna_cut_point_insert(seq, cut_point);

        backward_compat_compound = vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);

        free(sequence);
    }

    route = vrna_path_findpath(vc, s1, s2, maxkeep);

    return route;
}

#endif
