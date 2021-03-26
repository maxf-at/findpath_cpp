

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <algorithm>

// extern "C" {
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/landscape/paths.h>
// }

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

auto custom_findpath_method(vrna_fold_compound_t* fc, const char* s1, const char* s2, int width, int maxE)
    -> std::vector<sorted_path>;



