

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

// extern "C" {

// }

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

// #include <ViennaRNA/fold_compound.h>
// #include <ViennaRNA/landscape/paths.h>

}



struct sorted_move {
    int i; /* i,j>0 insert; i,j<0 delete */
    int j;
    // int when; /* 0 if still available, else resulting distance from start */
    int E;

    bool operator==(const sorted_move &other) const
    { return (i == other.i
                && j == other.j);
    }

};

struct sorted_path {
    std::vector<sorted_move> moves;
    int32_t                  max_en;

    // struct initializer, preallocate vector
    sorted_path(int init_size) {
        moves.resize(init_size);
        }
};

struct merge_path {
    std::vector<sorted_move> moves;
    int32_t                  current_s;
    int32_t                  current_en;
    short*                  current_ptable;

    int current_G1_node;
    int current_G2_node;

    // struct initializer, preallocate vector
    // sorted_path(int i2) {
    //     moves.resize(i2);
    //     }
};



using sorted_paths = std::vector<sorted_path>;

// class sorted_path
// {
//     // int size;
//     public:
//     //    sorted_path_2(double r) { radius = r; }
//        std::vector<sorted_move> moves;
//        int32_t                  max_en;
//        sorted_path(int init_size) {
//            moves.resize(init_size);
//            }

       
       
// };

// auto custom_findpath_method(vrna_fold_compound_t* fc, const char* s1, const char* s2, int width, int maxE)
//     -> std::vector<sorted_path>;

// function declarations...

// auto custom_findpath_method(vrna_fold_compound_t* fc, const char* s1, const char* s2, int width, int maxE) -> std::vector<sorted_path>;
// auto custom_findpath_method(vrna_fold_compound_t* fc, short* pt1, short* pt2, int width, int maxE) -> std::vector<sorted_path>;
auto custom_findpath_method(vrna_fold_compound_t* fc, short* pt1, short* pt2, int width, int maxE, bool init_direction) -> std::vector<sorted_path>;

