

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




/**
 *  @brief
 */
struct sorted_move {
    int i; /* i,j>0 insert; i,j<0 delete */
    int j;
    // int when; /* 0 if still available, else resulting distance from start */
    int E;
};

struct sorted_path {
    std::vector<sorted_move> path;
    int32_t                  max_en;
};

// auto custom_findpath_method(vrna_fold_compound_t* fc, const char* s1, const char* s2, int width, int maxE)
//     -> std::vector<sorted_path>;


std::vector<sorted_path> custom_findpath_method(vrna_fold_compound_t* fc, const char* s1, const char* s2, int width, int maxE);

