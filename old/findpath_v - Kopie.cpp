/* gcc -fopenmp -g3 -DTEST_FINDPATH findpath.c -o FINDpath -lRNA -lm -I../ -L./
 */

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <iomanip>
#include <cstring>

#include <algorithm>

// #include <execution>


#include "./icecream.hpp"
using icecream::ic;




using namespace std;

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

extern "C" {
	#include <limits.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>

	#include "ViennaRNA/cofold.h"
	#include "ViennaRNA/datastructures/basic.h"
	#include "ViennaRNA/fold.h"
	#include "ViennaRNA/fold_vars.h"
	#include "ViennaRNA/landscape/findpath.h"
	#include "ViennaRNA/model.h"
	#include "ViennaRNA/params/basic.h"
	#include "ViennaRNA/utils/basic.h"
	#include "ViennaRNA/utils/strings.h"
	#include "ViennaRNA/utils/structures.h"

}


#define def auto

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
typedef struct move
{
	int i; /* i,j>0 insert; i,j<0 delete */
	int j;
	int when; /* 0 if still available, else resulting distance from start */
	int E;
} move_t;

/**
 *  @brief
 */
typedef struct intermediate
{
	short *pt;	   /**<  @brief  pair table */
	int Sen;	   /**<  @brief  saddle energy so far */
	int curr_en;   /**<  @brief  current energy */
	move_t *moves; /**<  @brief  remaining moves to target */
} intermediate_t;

struct vrna_path_options_s
{
	unsigned int type;
	unsigned int method;

	int width;
};

double w_time;

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
PRIVATE int BP_dist;
PRIVATE move_t *path = NULL;
PRIVATE int path_fwd; /* 1: s1->s2, else s2 -> s1 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PRIVATE vrna_fold_compound_t *backward_compat_compound = NULL;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as
 * threadprivate
 */
#pragma omp threadprivate(BP_dist, path, path_fwd, backward_compat_compound)

#endif

#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE move_t *copy_moves(move_t *mvs);

PRIVATE int compare_ptable(const void *A, const void *B);

PRIVATE int compare_energy(const void *A, const void *B);

PRIVATE int compare_moves_when(const void *A, const void *B);

PRIVATE void free_intermediate(intermediate_t *i);

#ifdef TEST_FINDPATH

/* TEST_FINDPATH, COFOLD */
PRIVATE void usage(void);

#endif

PRIVATE int find_path_once(vrna_fold_compound_t *vc, short *pt1, short *pt2,
						   int maxl, int maxE);

PRIVATE int try_moves(vrna_fold_compound_t *vc, intermediate_t c, int maxE,
					  std::vector<intermediate_t> &next_2, int outer_num_next, int dist);

PRIVATE int try_moves(vrna_fold_compound_t *vc, intermediate_t c, int maxE,
					  std::array<intermediate_t,50000> &next_2, int outer_num_next, int dist);

PRIVATE int try_moves(vrna_fold_compound_t *vc, intermediate_t c, int maxE,
					  intermediate_t *next, int dist);





	// c = memcmp(a->pt, b->pt, a->pt[0] * sizeof(short));
	// fprintf(stderr, "compare c=%d, a= %d %d, b = %d %d", c, a->Sen, b->Sen, a->curr_en, b->curr_en);
	// if (c != 0){
	// 	fprintf(stderr, "return %d\n", c);
	// 	return c;
	// }
	// if ((a->Sen - b->Sen) != 0){
	// 	fprintf(stderr, "return %d\n", a->Sen - b->Sen);
	// 	return a->Sen - b->Sen;		
	// 	}

	// fprintf(stderr, "return %d\n", a->curr_en - b->curr_en);
	// return a->curr_en - b->curr_en;


auto compare_ptables_n(const intermediate_t &a, const intermediate_t &b) -> bool
{
	// Compares the first num bytes of the block of memory pointed
	// returning zero if they all match or a value different from zero representing which is greater if they do not.
	auto c = memcmp(a.pt, b.pt, a.pt[0] * sizeof(short));
	// fprintf(stderr, "compare c=%d, a= %d %d, b = %d %d", c, a.Sen, b.Sen, a.curr_en, b.curr_en);
	if (c != 0){
		// fprintf(stderr, "return false %d\n", c);
		if (c>0)
			return false;
		else
		{
			return true;
		}		
	}
	auto s_diff = a.Sen - b.Sen;
	if (s_diff != 0){
		// fprintf(stderr, "return sd %d\n", a.Sen - b.Sen);
		if (s_diff > 0)
			return false;
		else
			return true;		
		}
	// fprintf(stderr, "return ed %d\n", a.curr_en - b.curr_en);
	auto e_diff = a.curr_en - b.curr_en;
	if (e_diff > 0)
			return true;
		else
			return false;	
}

auto compare_energy_n(const intermediate_t &a, const intermediate_t &b) -> bool
{
	// compare saddle energy or current energy
	auto s_diff = a.Sen - b.Sen;
	if (s_diff != 0){
		// fprintf(stderr, "return sd %d\n", a.Sen - b.Sen);
		if (s_diff > 0)
			return false;
		else
			return true;		
		}
	// fprintf(stderr, "return ed %d\n", a.curr_en - b.curr_en);
	auto e_diff = a.curr_en - b.curr_en;
	if (e_diff > 0)
			return true;
		else
			return false;	
}

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void vrna_path_free(vrna_path_t *path)
{
	vrna_path_t *tmp = path;

	if (tmp)
	{
		if (tmp->type == VRNA_PATH_TYPE_DOT_BRACKET)
		{
			while (tmp->s)
			{
				free(tmp->s);
				tmp++;
			}
		}
		else if (tmp->type == VRNA_PATH_TYPE_MOVES)
		{
			while ((tmp->move.pos_5))
			{
				vrna_move_list_free(tmp->move.next);
				tmp++;
			}
		}

		free(path);
	}
}

PUBLIC int vrna_path_findpath_saddle(vrna_fold_compound_t *vc, const char *s1,
									 const char *s2, int width)
{
	return vrna_path_findpath_saddle_ub(vc, s1, s2, width, INT_MAX - 1);
}

PUBLIC int vrna_path_findpath_saddle_ub(vrna_fold_compound_t *vc,
										const char *s1, const char *s2,
										int width, int maxE)
{
	int maxl;
	short *ptr, *pt1, *pt2;
	move_t *bestpath = NULL;
	int dir;

	path_fwd = dir = 0;
	pt1 = vrna_ptable(s1);
	pt2 = vrna_ptable(s2);

	maxl = 1;
	do
	{
		int saddleE;
		// switch forward / backward
		path_fwd = !path_fwd;
		if (maxl > width)
			maxl = width;

		if (path)
			free(path);

		saddleE = find_path_once(vc, pt1, pt2, maxl, maxE);

		

		if (saddleE < maxE)
		{
			maxE = saddleE;
			if (bestpath)
				free(bestpath);

			bestpath = path;
			path = NULL;
			dir = path_fwd;
		}
		else
		{
			free(path);
			path = NULL;
		}

		ptr = pt1;
		pt1 = pt2;
		pt2 = ptr;
		maxl *= 2;
	} while (maxl < 2 * width);

	/* (re)set some globals */
	path = bestpath;
	path_fwd = dir;

	free(pt1);
	free(pt2);

	return maxE;
}

PUBLIC vrna_path_t *vrna_path_findpath(vrna_fold_compound_t *fc, const char *s1,
									   const char *s2, int width)
{
	struct vrna_path_options_s *opt;
	vrna_path_t *path;

	opt = vrna_path_options_findpath(width, VRNA_PATH_TYPE_DOT_BRACKET);
	path = vrna_path_direct_ub(fc, s1, s2, INT_MAX - 1, opt);

	free(opt);

	return path;
}

PUBLIC vrna_path_t *vrna_path_findpath_ub(vrna_fold_compound_t *fc,
										  const char *s1, const char *s2,
										  int width, int maxE)
{
	struct vrna_path_options_s *opt;
	vrna_path_t *path;

	opt = vrna_path_options_findpath(width, VRNA_PATH_TYPE_DOT_BRACKET);
	path = vrna_path_direct_ub(fc, s1, s2, maxE, opt);

	free(opt);

	return path;
}

PUBLIC struct vrna_path_options_s *
vrna_path_options_findpath(int width, unsigned int type)
{
	struct vrna_path_options_s *options =
		(struct vrna_path_options_s *)vrna_alloc(
			sizeof(struct vrna_path_options_s));

	options->type = type;
	options->method = PATH_DIRECT_FINDPATH;
	options->width = width;

	return options;
}

PUBLIC void vrna_path_options_free(struct vrna_path_options_s *options)
{
	free(options);
}

PRIVATE vrna_path_t *findpath_method(vrna_fold_compound_t *fc, const char *s1,
									 const char *s2, int width, int maxE,
									 unsigned int return_type)
{
	// returns a route (vrna_path_t)
	// reconstructs i,j pairs to generate a string path
	
	int E, d;
	float last_E;
	vrna_path_t *route = nullptr;

	E = vrna_path_findpath_saddle_ub(fc, s1, s2, width, maxE);

	/* did we find a better path than one with saddle maxE? */
	if (E < maxE)
	{
		route = (vrna_path_t *)vrna_alloc((BP_dist + 2) * sizeof(vrna_path_t));

		qsort(path, BP_dist, sizeof(move_t), compare_moves_when);

		switch (return_type)
		{
		case VRNA_PATH_TYPE_MOVES: // unused?
			if (path_fwd)
			{
				last_E = vrna_eval_structure(fc, s1);
				for (d = 0; d < BP_dist; d++)
				{
					route[d].type = return_type;
					route[d].move = vrna_move_init(path[d].i, path[d].j);
					route[d].en = (path[d].E / 100.0) - last_E;
					last_E = path[d].E / 100.0;
				}

				route[BP_dist].type = return_type;
				route[BP_dist].move = vrna_move_init(0, 0);
			}
			else
			{
				last_E = vrna_eval_structure(fc, s2);
				for (d = 0; d < BP_dist; d++)
				{
					route[BP_dist - d - 2].type = return_type;
					route[BP_dist - d - 2].move = vrna_move_init(path[d].i, path[d].j);
					route[BP_dist - d - 2].en = last_E - (path[d].E / 100.0);
					last_E = path[d].E / 100;
				}

				route[BP_dist].type = return_type;
				route[BP_dist].move = vrna_move_init(0, 0);
			}

			break;

		case VRNA_PATH_TYPE_DOT_BRACKET:
			/* fall through */

		default:
			if (path_fwd)
			{
				/* memorize start of path */
				route[0].type = return_type;
				route[0].s = strdup(s1);
				route[0].en = vrna_eval_structure(fc, s1);

				for (d = 0; d < BP_dist; d++)
				{
					int i, j;
					route[d + 1].type = return_type;
					route[d + 1].s = strdup(route[d].s);
					i = path[d].i;
					j = path[d].j;
					if (i < 0)
					{
						/* delete */
						route[d + 1].s[(-i) - 1] = route[d + 1].s[(-j) - 1] = '.';
					}
					else
					{
						route[d + 1].s[i - 1] = '(';
						route[d + 1].s[j - 1] = ')';
					}

					route[d + 1].en = path[d].E / 100.0;
				}
			}
			else
			{
				/* memorize start of path */
				route[0].type = return_type;
				route[BP_dist].s = strdup(s2);
				route[BP_dist].en = vrna_eval_structure(fc, s2);

				for (d = 0; d < BP_dist; d++)
				{
					int i, j;
					route[BP_dist - d - 1].type = return_type;
					route[BP_dist - d - 1].s = strdup(route[BP_dist - d].s);
					i = path[d].i;
					j = path[d].j;
					if (i < 0)
					{
						/* delete */
						route[BP_dist - d - 1].s[(-i) - 1] =
							route[BP_dist - d - 1].s[(-j) - 1] = '.';
					}
					else
					{
						route[BP_dist - d - 1].s[i - 1] = '(';
						route[BP_dist - d - 1].s[j - 1] = ')';
					}
					route[BP_dist - d - 1].en = path[d].E / 100.0;
				}
			}

			break;
		}

		// does the same as main verbose 
		// fprintf(stderr, "\n%s\n%s\n\n", s1, s2);
		// for (d = 0; d <= BP_dist; d++)
		// 	fprintf(stderr, "%s %6.2f\n", route[d].s, route[d].en);
		// // fprintf(stderr, "%d\n", *num_entry);

#if _DEBUG_FINDPATH_
		if (return_type & VRNA_PATH_TYPE_DOT_BRACKET)
		{
			fprintf(stderr, "\n%s\n%s\n%s\n\n", seq, s1, s2);
			for (d = 0; d <= BP_dist; d++)
				fprintf(stderr, "%s %6.2f\n", route[d].s, route[d].en);
			fprintf(stderr, "%d\n", *num_entry);
		}

#endif
	}

	free(path);
	path = NULL;

	return route;
}

PUBLIC vrna_path_t *vrna_path_direct(vrna_fold_compound_t *fc, const char *s1,
									 const char *s2,
									 struct vrna_path_options_s *options)
{
	return vrna_path_direct_ub(fc, s1, s2, INT_MAX - 1, options);
}


PUBLIC vrna_path_t *vrna_path_direct_ub(vrna_fold_compound_t *fc,
										const char *s1, const char *s2,
										int maxE,
										struct vrna_path_options_s *options)
{
	int E, d;
	struct vrna_path_options_s *o;
	vrna_path_t *route = NULL;

	/* we default to findpath method */
	o = options ? options
				: vrna_path_options_findpath(10, VRNA_PATH_TYPE_DOT_BRACKET);

	switch (o->method)
	{
	case PATH_DIRECT_FINDPATH:
	/* fall through */
	default:
		route = findpath_method(fc, s1, s2, o->width, maxE, o->type);
		break;
	}

	if (!options)
		vrna_path_options_free(o);

	return route;
}

// #ifdef TEST_FINDPATH

PUBLIC void print_path(const char *seq, const char *struc)
{
	int d;
	char *s;

	s = strdup(struc);
	if (cut_point == -1)
	{
		printf("%s\n%s\n", seq, s);
	}
	/* printf("%s\n%s %6.2f\n", seq, s, vrna_eval_structure_simple(seq,s)); */
	else
	{
		char *pstruct, *pseq;
		pstruct = vrna_cut_point_insert(s, cut_point);
		pseq = vrna_cut_point_insert(seq, cut_point);
		printf("%s\n%s\n", pseq, pstruct);
		/* printf("%s\n%s %6.2f\n", pseq, pstruct,
	 * vrna_eval_structure_simple(seq,s)); */
		free(pstruct);
		free(pseq);
	}

	qsort(path, BP_dist, sizeof(move_t), compare_moves_when);
	for (d = 0; d < BP_dist; d++)
	{
		int i, j;
		i = path[d].i;
		j = path[d].j;
		if (i < 0)
		{
			/* delete */
			s[(-i) - 1] = s[(-j) - 1] = '.';
		}
		else
		{
			s[i - 1] = '(';
			s[j - 1] = ')';
		}

		/* printf("%s %6.2f - %6.2f\n", s, vrna_eval_structure_simple(seq,s),
	 * path[d].E/100.0); */
	}
	free(s);
}


int main(int argc, char *argv[])
{
	// vrna_message_error("usage: findpath.c  [-m depth] [-d[0|1|2]] [-v]");

	printf("main function\n");

	ic.line_wrap_width(990);
    // ic.prefix("ic: ");
	ic.prefix("\033[1;36mDEBUG: \033[0m");	
    ic.show_c_string(false);
	ic.disable();

	// debug

	// char *line,
	char *seq, *s1, *s2;
	float E;
	int maxkeep = 1000;
	maxkeep = 500; // 30 is enough

	int verbose = 0;
	int i = 0;

	vrna_path_t *route, *r;

	for (i = 1; i < argc; i++)
	{
		switch (argv[i][1])
		{
		case 'm':
			if (strcmp(argv[i], "-m") == 0)
				printf("m\n");
			sscanf(argv[++i], "%d", &maxkeep);

			break;
		case 'v':
			printf("v\n");
			verbose = !strcmp(argv[i], "-v");
			break;
		case 'd':
			printf("d\n");
			if (strcmp(argv[i], "-d") == 0)
				sscanf(argv[++i], "%d", &dangles);

			break;
		default:
			verbose = 1;
			// usage();
		}
	}

	cut_point = -1;
	// cut_point = 10;

	// # sequence = "CCUCCCAUCGCUUUGAAUGACGGCGCAAUGAGGCCCGGAUAUAACUCGGGGGAACAAGCAGCUGCAACAAGUUUGAACGGUUUCUCCUGCAGAAUGAGUGGUACCCUAACAUGCUAGGCUGCACUGGAAGCAAUGUCCCAUUGACCACUC";
    // # s1 =       "((((....((((.........))))....))))(((.((......)).)))..........(((((........(((....)))...)))))...(((((((......(((((((((......))...))).)))).......)))))))";
    // # s2 =       "((((.....(((.........))).....))))...........((((((((((....(((((......)))).)......)))))))).))...(((((((.........(.((((......)))).).(((((...))))))))))))";


	
	printf("seq\n");
	// line      = vrna_read_line(stdin);
	// line = "AAAA";
	// string l = "AUCAAGCUUUUGUAGCUAAACCUACGCGGGUUUGGCGUAGGGGGU";
	// string l = "CUUCCUCUCUGUCUCUCCUUGNNNCCAUGGCUAUAAGNNNCUCGAGGGUACCAGCAU";
	// string l = "AACGCCGCGCCGGCCAGGCUGCCAGUUGUCUUAUUGAGUUGGCCAUCUACGUAAUACUUUUGUAGCCGACGACUAUAACG";
	string l = "AUGAACUUUCGUCCGAUCCCUGGAAAGGAGAUUGCUAGGGCCUAGAUCCACACAGACUAGACUGGCCGUUGGAUGAUAAAGGUUGCUGAUCUCAGAGUAAGCGUAACACCCGUAUACGACUAGUUCAAUGUUUUAGCUAAAGGAUCAGCCGGAUGCCAGUCCAUAGCUCAUUGGUUUCGGUGUGUGUGUCUUAGGAUUCU";
	// string l = "UGGUCGGUCCUUAACCUAGACCUAGUACUUUUGACGUAGCCUGCUUUUAAACUUACAACUCAUGUUCAAC";



	printf("seq\n");
	char* line = const_cast<char*>(l.c_str());
	seq = vrna_cut_point_remove(line, &cut_point);
	// free(line);

	printf("s1\n");
	// line  = vrna_read_line(stdin);
	// line = "(())";
	// l = ".....((((((...(((((((((....)))))))))...))))))";
	// l = ".........((((((((...(...)..........((...)).))))).....))).";
	// l = "(((...(..............)..)))...((((...((((((....((((.((....))))))))))))....))))..";
	l = "((((....)))).....((((((.((.....)).))))))(((((((.((((((.....((((((.((((((.((.((..(.(((((........))))).).)).)).)))...))).))))))(((((....(((((..((((..((.....))..)))).)))))))))).......)))))).)))).))).....";
	// l = "(((((((........)).)))).).....(..(((((........................)))))..).";
	line = const_cast<char*>(l.c_str());
	s1 = vrna_cut_point_remove(line, &cut_point);
	// free(line);

	printf("s2\n");
	// line  = vrna_read_line(stdin);
	// line = "()()";
	// l = "....((((.....))))...(((((((.......)))))))....";
	// l = "...((((..(......((.((...).).)).....)(...)..))))..........";
	// l = "(((...(..............)..)))..........((((((....((((.........))))))))))..........";
	l = ".....((((.((..(((((((....))).)))))).))))((.(((..((((((.(((((((((((...(((((..........((((((((.((...(((((((((...((....))....)))..))))))...))...)))))))).........)))))..)))....))))).))).)))))))))..)).....";
	// l = ".((((((.......))..))))..(((..((((....((....))..))))..))).((....)).....";
	line = const_cast<char*>(l.c_str());
	s2 = vrna_cut_point_remove(line, &cut_point);
	// free(line);

   



	// E = find_saddle(seq, s1, s2, maxkeep);
	// printf("saddle_energy = %6.2f\n", E / 100.);

	if (verbose)
	{
		if (path_fwd)
			print_path(seq, s1);
		else
			print_path(seq, s2);

		free(path);
		path = NULL;
		// main function call
		// ? vrna_path_findpath(vc, s1, s2, maxkeep);
		// ? ->

		E = INT_MIN;

		route = get_path(seq, s1, s2, maxkeep);
		for (r = route; r->s; r++)
		{
			if (cut_point == -1)
			{
				if (r->en > E)
					E = r->en;

				printf("%s %6.2f\n", r->s, r->en);
				// printf("%s %6.2f - %6.2f\n", r->s,
				// vrna_eval_structure_simple(seq,r->s), r->en);
			}
			else
			{
				char *pstruct;
				pstruct = vrna_cut_point_insert(r->s, cut_point);
				printf("%s %6.2f\n", pstruct, r->en);
				// printf("%s %6.2f - %6.2f\n", pstruct,
				// vrna_eval_structure_simple(seq,r->s), r->en);
				free(pstruct);
			}

			free(r->s);
		}

		printf("saddle_energy = %6.2f\n", E);

		free(route);
	}

	free(seq);
	free(s1);
	free(s2);
	return EXIT_SUCCESS;
}

static void usage(void)
{
	vrna_message_error("usage: findpath.c  [-m depth] [-d[0|1|2]] [-v]");
}

// #endif

/*
 #################################
 # STATIC helper functions below #
 #################################
*/

// PRIVATE int try_moves(vrna_fold_compound_t *vc, intermediate_t c, int maxE,
// 					  std::vector<intermediate_t> &next_2, int outer_num_next, int dist)

// PRIVATE int try_moves(vrna_fold_compound_t *vc, intermediate_t c, int maxE,
// 					  std::array<intermediate_t,50000> &next_2, int outer_num_next, int dist)
PRIVATE int try_moves(vrna_fold_compound_t *vc, intermediate_t c, int maxE,
					  intermediate_t *next,  int dist)
{

	// fprintf(stderr, "start try_moves, maxE=%d, dist=%d, c.saddle=%d, %d, c.en=%d \n", maxE, dist, c.Sen, *c.pt, c.curr_en);

	int *loopidx, len, num_next = 0, en, oldE;
	move_t *mv;
	short *pt;

	len = c.pt[0];
	loopidx = vrna_loopidx_from_ptable(c.pt);
	oldE = c.Sen;
	// create new p table for moves
	for (mv = c.moves; mv->i != 0; mv++)
	{
		int i, j;
		if (mv->when > 0)
			continue;

		i = mv->i;
		j = mv->j;
		pt = (short *)vrna_alloc(sizeof(short) * (len + 1));
		memcpy(pt, c.pt, (len + 1) * sizeof(short));
		if (j < 0)
		{
			/*it's a delete move */
			pt[-i] = 0;
			pt[-j] = 0;
		}
		else
		{
			/* insert move */
			if ((loopidx[i] == loopidx[j]) && /* i and j belong to same loop */
				(pt[i] == 0) && (pt[j] == 0)  /* ... and are unpaired */
			)
			{
				pt[i] = j;
				pt[j] = i;
			}
			else
			{
				free(pt);
				continue; /* illegal move, try next; */
			}
		}

#ifdef LOOP_EN
		en = c.curr_en + vrna_eval_move_pt(vc, c.pt, i, j);
#else
		en = vrna_eval_structure_pt(vc, pt);
#endif
		if (en < maxE)
		{
			next[num_next].Sen = (en > oldE) ? en : oldE;
			next[num_next].curr_en = en;
			next[num_next].pt = pt;

			// next_2[num_next+outer_num_next].Sen = (en > oldE) ? en : oldE;
			// next_2[num_next+outer_num_next].curr_en = en;
			// next_2[num_next+outer_num_next].pt = pt;

			mv->when = dist;
			mv->E = en;
			next[num_next].moves = copy_moves(c.moves);
			// next_2[num_next+outer_num_next].moves = copy_moves(c.moves);
			mv->when = 0;

			num_next++;
		}
		else
		{
			free(pt);
		}
	}
	free(loopidx);
	return num_next;
}

PRIVATE int find_path_once(vrna_fold_compound_t *vc, short *pt1, short *pt2,
						   int maxl, int maxE)
{

	// calls try_moves - iterate over BP_dist below

	IC("fp once");
	IC(maxE, maxl, *pt1, *pt2);


	move_t *mlist;
	int i, len, d, dist = 0, result;
	short *pt;
	// intermediate_t *current, *next;
	intermediate_t *current;

	len = (int)pt1[0];
	pt = vrna_ptable_copy(pt1);

	mlist = (move_t *)vrna_alloc(sizeof(move_t) * len); /* bp_dist < n */

	// find fitting i,j, save into mlist / distance class

	for (i = 1; i <= len; i++)
	{
		if (pt[i] != pt2[i])
		{
			if (i < pt[i])
			{
				/* need to delete this pair */
				mlist[dist].i = -i;
				mlist[dist].j = -pt[i];
				mlist[dist++].when = 0;
			}

			if (i < pt2[i])
			{
				/* need to insert this pair */
				mlist[dist].i = i;
				mlist[dist].j = pt2[i];
				mlist[dist++].when = 0;
			}
		}
	}

	BP_dist = dist;
	current = (intermediate_t *)vrna_alloc(sizeof(intermediate_t) * (maxl + 1));
	current[0].pt = pt;
	current[0].Sen = current[0].curr_en = vrna_eval_structure_pt(vc, pt);
	current[0].moves = mlist;

	// initialize array of paths
	// next =
	// 	(intermediate_t *)vrna_alloc(sizeof(intermediate_t) * (dist * maxl + 1));

	// std::vector<intermediate_t> next_2;
	// next_2.resize(dist * maxl + 1);

	std::array<intermediate_t, 50000> next_2;

	// intermediate_t *next_2;
	// next_2 =
	// 	(intermediate_t *)vrna_alloc(sizeof(intermediate_t) * (dist * maxl + 1));

	// fprintf(stderr, "call try_moves for %d %d \n", dist, current->Sen);

	// iterate over BP_dist (always e.g. 26 times)
	for (d = 1; d <= dist; d++)
	{
		/* go through the distance classes */
		int c, u, num_next = 0;
		intermediate_t *cc;

		// fprintf(stderr, "bp_dist = %d / go through distance classes to %d \n", d, dist);

		for (c = 0; current[c].pt != NULL; c++){
		// for (c = 0; c >= 10; c++)


			IC(d, c, num_next, "bp_dist");
			// debug

			std::string temp_str = "(((...(((((.....)).)))..)))...((((...((((((....((((.((....))))))))))))....))))..";
			for (int m = 1; m<80; m++)
				{
				int p1 = current[c].pt[m];
				if (p1 == 0)
				{
					temp_str[m] = '.';			
					}
				else if (p1 > m)
					{
					temp_str[m] = '(';
					}
				else 
					{
					temp_str[m] = ')';
					}					
				}

			// float test_E = vrna_eval_structure_pt(vc, current[c].pt) / 100.0;
			// IC(temp_str, test_E, c);
			// fprintf(stderr, "before %d  \n", num_next);
			// num_next += try_moves(vc, current[c], maxE, next_2, num_next, d);
			num_next += try_moves(vc, current[c], maxE, next_2.begin()+num_next, d);

			// fprintf(stderr, "received %d  \n", num_next);

			// num_next += try_moves(vc, current[c], 999 , next + num_next, d);
		}
		if (num_next == 0)
		{
			IC("break", num_next);
			for (cc = current; cc->pt != NULL; cc++)
				free_intermediate(cc);
			current[0].Sen = INT_MAX;
			break;
		}


		


		IC ("received candidates:", num_next);

		/* remove duplicates via sort|uniq
	 	 * if this becomes a bottleneck we can use a hash instead */
		
		// quicksort compare_ptable
		// std::sort();


		// typedef struct intermediate
		// {
		// 	short *pt;	   /**<  @brief  pair table */
		// 	int Sen;	   /**<  @brief  saddle energy so far */
		// 	int curr_en;   /**<  @brief  current energy */
		// 	move_t *moves; /**<  @brief  remaining moves to target */
		// } intermediate_t;		
		// qsort params:
		// a base pointer, which points to the start of a contiguous block of data,
		// the number of elements in the block,
		// the size of one data member, and
		// a function that compares two data values.
		// std::sort(next, next+num_next, compare_ptable);
		
		// next+num_next;
		
		auto t1 = std::chrono::high_resolution_clock::now();
		
		// std::sort(std::execution::par_unseq, next, next+num_next, compare_ptables);

		// std::sort(next, next+num_next, compare_ptables_n);

		std::sort(next_2.begin(), next_2.begin()+num_next, compare_ptables_n);

		// qsort(next, num_next, sizeof(intermediate_t), compare_ptable);



		// fprintf(stderr, "2received candidates: num_next = %d  \n", num_next);

		// int a;
		// std::string line;

		//for (u = 0, c = 1; c < num_next; c++)
		for (u = 0, c = 0; c < num_next; c++)
		{		
			// char temp_str[] = "(((...(((((.....)).)))..)))...((((...((((((....((((.((....))))))))))))....))))..";
			std::string temp_str = "(((...(((((.....)).)))..)))...((((...((((((....((((.((....))))))))))))....))))..";
			for (int m = 1; m<80; m++)
				{
				int p1 = next_2[u].pt[m];
				if (p1 == 0)
				{
					temp_str[m] = '.';			
					}
				else if (p1 > m)
					{
					temp_str[m] = '(';
					}
				else 
					{
					temp_str[m] = ')';
					}					
				}
			// float test_E = vrna_eval_structure_pt(vc, next[u].pt) / 100.0;
		
			if (c == 0) {
				// no need to duplicate check the first entry
				continue;			
			}
		
			// mem compare u vs c
			if (memcmp(next_2[u].pt, next_2[c].pt, sizeof(short) * len) != 0)
			{
				// IC(temp_str, test_E, "A");
				++u;
				// next[u] = next[c];
				next_2[u] = next_2[c];
			}
			else{
				// remove duplicate
				// IC(temp_str, test_E, "D");
				
				// free_intermediate(next + c);

				// free(next_2[c].pt); `???
				// free(next_2[c].moves);
				next_2[c].pt = nullptr;
				next_2[c].moves = nullptr;
				next_2[c].Sen = INT_MAX;

			}
				
		}
		num_next = u + 1;


		
		// quicksort compare_energy

		// std::sort(next, next+num_next, compare_energy_n);
		
		intermediate_t* addr = next_2.begin();

		std::sort(addr, addr+num_next, compare_energy_n);

		// qsort(next, num_next, sizeof(intermediate_t), compare_energy);
		
		// if (d==50)
		// {
		// 	printf("%d\n", next[0].Sen);
		// 	printf("%d\n", next_2[0].Sen);

		// 	std::string line;
		// 	std::getline(std::cin, line);
		// }
		


		// new
		for (u = 0, c = 0; c < num_next; c++)
		{
			std::string temp_str = "(((...(((((.....)).)))..)))...((((...((((((....((((.((....))))))))))))....))))..";
			for (int m = 1; m<80; m++)
				{
				int p1 = next_2[c].pt[m];
				if (p1 == 0)
				{
					temp_str[m] = '.';			
					}
				else if (p1 > m)
					{
					temp_str[m] = '(';
					}
				else 
					{
					temp_str[m] = ')';
					}					
				}
			
			// float test_E = vrna_eval_structure_pt(vc, next[c].pt) / 100.0;
			// float test_E = vrna_eval_structure(vc, temp_str);
			// IC(temp_str, test_E, "S");
		}

		

		/* free the old stuff */
		for (cc = current; cc->pt != NULL; cc++)
			free_intermediate(cc);
		for (u = 0; u < maxl && u < num_next; u++){
			
			IC(u, "current=next");
			// fprintf(stderr, "current %d = next \n", u);	
			current[u] = next_2[u];
			}
		for (; u < num_next; u++)
		{
			next_2[u].pt = nullptr;
			next_2[u].moves = nullptr;
			next_2[u].Sen = INT_MAX;
			// free_intermediate(next + u);
		}
			




		num_next = 0;

		auto t2 = std::chrono::high_resolution_clock::now();
		auto c_time = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
		w_time += c_time;


	}

	// free(next);


	// start
	ic.enable();
	IC(w_time/1000000);
	ic.disable();

	path = current[0].moves;
	result = current[0].Sen;
	free(current[0].pt);
	free(current);
	return result;
}

PRIVATE void free_intermediate(intermediate_t *i)
{
	free(i->pt);
	free(i->moves);
	i->pt = NULL;
	i->moves = NULL;
	i->Sen = INT_MAX;
}

PRIVATE int compare_ptable(const void *A, const void *B)
{
	intermediate_t *a, *b;
	int c;
	a = (intermediate_t *)A;
	b = (intermediate_t *)B;
	c = memcmp(a->pt, b->pt, a->pt[0] * sizeof(short));
	// fprintf(stderr, "compare c=%d, a= %d %d, b = %d %d", c, a->Sen, b->Sen, a->curr_en, b->curr_en);
	if (c != 0){
		// fprintf(stderr, "return %d\n", c);
		return c;
	}
	if ((a->Sen - b->Sen) != 0){
		// fprintf(stderr, "return %d\n", a->Sen - b->Sen);
		return a->Sen - b->Sen;		
		}
	// fprintf(stderr, "return %d\n", a->curr_en - b->curr_en);
	return a->curr_en - b->curr_en;
}



PRIVATE int compare_energy(const void *A, const void *B)
{
	intermediate_t *a, *b;

	a = (intermediate_t *)A;
	b = (intermediate_t *)B;

	// compare saddle energy or current energy

	if ((a->Sen - b->Sen) != 0)
		return a->Sen - b->Sen;

	return a->curr_en - b->curr_en;
}

PRIVATE int compare_moves_when(const void *A, const void *B)
{
	move_t *a, *b;

	a = (move_t *)A;
	b = (move_t *)B;

	return a->when - b->when;
}

PRIVATE move_t *copy_moves(move_t *mvs)
{
	move_t *new2;

	new2 = (move_t *)vrna_alloc(sizeof(move_t) * (BP_dist + 1));

	memcpy(new2, mvs, sizeof(move_t) * (BP_dist + 1));
	
	return new2;
}

/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC void free_path(vrna_path_t *path)
{
	vrna_path_free(path);
}

PUBLIC int find_saddle(const char *seq, const char *s1, const char *s2,
					   int width)
{
	int maxE;
	char *sequence;
	vrna_fold_compound_t *vc;
	vrna_md_t md, *md_p;

	vc = NULL;
	set_model_details(&md);

	if (backward_compat_compound)
	{
		if (!strcmp(seq, backward_compat_compound->sequence))
		{
			/* check if sequence is the same as before */
			md.window_size = backward_compat_compound->length;
			md.max_bp_span = backward_compat_compound->length;
			md_p = &(backward_compat_compound->params->model_details);
			if (!memcmp(&md, md_p, sizeof(vrna_md_t))) /* check if model_details are
													the same as before */
				vc =
					backward_compat_compound; /* re-use previous vrna_fold_compound_t */
		}
	}

	if (!vc)
	{
		vrna_fold_compound_free(backward_compat_compound);

		sequence = vrna_cut_point_insert(seq, cut_point);

		backward_compat_compound = vc =
			vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);

		free(sequence);
	}

	maxE = vrna_path_findpath_saddle(vc, s1, s2, width);

	return maxE;
}

PUBLIC vrna_path_t *get_path(const char *seq, const char *s1, const char *s2,
							 int maxkeep)
{
	vrna_path_t *route = NULL;
	char *sequence = NULL;
	vrna_fold_compound_t *vc = NULL;
	vrna_md_t md, *md_p;

	set_model_details(&md);

	if (backward_compat_compound)
	{
		if (!strcmp(seq, backward_compat_compound->sequence))
		{
			/* check if sequence is the same as before */
			md.window_size = backward_compat_compound->length;
			md.max_bp_span = backward_compat_compound->length;
			md_p = &(backward_compat_compound->params->model_details);
			if (!memcmp(&md, md_p, sizeof(vrna_md_t))) /* check if model_details are
													the same as before */
				vc =
					backward_compat_compound; /* re-use previous vrna_fold_compound_t */
		}
	}

	if (!vc)
	{
		vrna_fold_compound_free(backward_compat_compound);

		sequence = vrna_cut_point_insert(seq, cut_point);

		backward_compat_compound = vc =
			vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);

		free(sequence);
	}

	route = vrna_path_findpath(vc, s1, s2, maxkeep);

	return route;
}

#endif
