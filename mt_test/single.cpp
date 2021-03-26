

// g++ single.cpp -o3 -std=c++20 -o single -lRNA -lm -I./ -I/usr/include/ViennaRNA -DEBUG_FINDPATH
// -fno-lto -pthread  -fopenmp

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <unordered_map>

// replace this with std::format (C++23?)
// #define FMT_HEADER_ONLY
// #include "fmt/color.h"
// #include "fmt/format.h"

#include <algorithm>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include <chrono>
#include <future>
#include <thread>

// #include <execution>

extern "C" {
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cofold.h"
#include "datastructures/basic.h"
#include "fold.h"
#include "fold_vars.h"
#include "landscape/findpath.h"
#include "model.h"
#include "params/basic.h"
#include "utils/basic.h"
#include "utils/strings.h"
#include "utils/structures.h"
}

auto test_func() -> int
{
    // # recursive 90 nt example
    const char* seq = const_cast<char*>(
        "CUCCGUUCGGCACAGUGGGAUUCAGACUUCCUGCCGCUGGGAGAAACGCGCGGUUCGGGGUGUAAUCAUGGUUUAAGCCUCCGAAGCCC"
        "A");
    const char* s1 = const_cast<char*>(
        "((((...((((...((((((........)))))).))))))))........((((((((((..(((....)))...))).))))).))."
        ".");
    const char* s2 = const_cast<char*>(
        "((((...(((((.....(((........))))))))...))))........(((((((((..(.(........).)..))))).))))."
        ".");

    // ((((...(((((.....(((........))))))))...))))........(((((((((..(.(........).)..))))).)))).

    // -1 -42

    std::cout << "entering test_func, thread id: " << std::this_thread::get_id() << "\n";

    vrna_fold_compound_t* fc = nullptr;
    vrna_md_t             md1;
    set_model_details(&md1);

    fc          = vrna_fold_compound(seq, &md1, VRNA_OPTION_EVAL_ONLY);
    short* pt_1 = vrna_ptable(s1);
    short* pt_2 = vrna_ptable(s2);

    int iterations = 1000000;
    int en         = 0;

    int i = -1;
    int j = -41;

    for (int i = 0; i < iterations; i++) {
        en += vrna_eval_structure_pt(fc, pt_1);
    }

    return en;
}

auto main(int argc, char** argv) -> int
{
    std::cout << "entering main function\n";

    auto start = std::chrono::high_resolution_clock::now();

    int en0 = test_func();

    auto                          finish  = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "single thread time elapsed:" << elapsed.count() << "\n";

}