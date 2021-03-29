
// g++ main.cpp -o3 -std=c++20 -o main findpath_B.cpp -lRNA -lm -I./ -I./ViennaRNA/ -DEBUG_FINDPATH
// -fno-lto

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <unordered_map>

// replace this with std::format (C++23?)
#define FMT_HEADER_ONLY
#include "fmt/color.h"
#include "fmt/format.h"

#include <algorithm>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include <chrono>

#include <future>
#include <thread>

// #include <execution>

#include "./icecream.hpp"
using icecream::ic;
using namespace std;

#include "findpath_B.hpp"

#include "s_graph.hpp"

#include "findpath_header.hpp"

auto main(int argc, char** argv) -> int
{
    // cout << "Hello World! ";
    // cout << "I'm a C++ program";

    ic.line_wrap_width(990);
    ic.prefix("\033[1;36mDEBUG: \033[0m");
    ic.show_c_string(false);

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
                // default: usage();
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

    vrna_fold_compound_t* fc1 = nullptr;
    vrna_md_t             md1;
    // set model params
    set_model_details(&md1);

    vrna_fold_compound_t* fc2 = nullptr;
    vrna_md_t             md2;
    // set model params
    set_model_details(&md2);

    // fc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);

    // result = custom_findpath_method(fc, s1, s2, maxkeep, INT_MAX - 1);
    // cout << result[0].max_en;

    // init graph

    // seq = const_cast<char*>(
    //     "AGUUCAGCUACUUACACUGCGGCCUACGGGGGAAUCACUGACCGCGUUAUACCCUAUAUUCGUAUAUGCACGUACAAUACCAACAGUAGG"
    //     "AAUACUUUCUGGGGAGUAGAGUAGAAAUCUUGGCCCCCGACUUGACUGUUAAUCCUCUAC");
    // s1 = const_cast<char*>(
    //     ".((...(((...........)))..))(((((..(((......(((.(((((.........))))))))...(((..((((..(((.((("
    //     "....))).)))..).)))..))).......))))))))(((......)))..........");
    // s2 = const_cast<char*>(
    //     "...........(((....(((((...((((((..(((......((.((((.(((((.((((.((..((...(((...)))...)).)).)"
    //     ")))......))))).)))).))........)))))))))....).)))))))........");

    // // ex2
    // seq = const_cast<char*>(
    //     "AGGGGAAAGUUAUCGAACAGACCUUGUUUGGCUUAUGGCACUUGGGUAGAUUUUCUCUACCGUUGUCAUGAAAACGCACCGCAAGAGAAG"
    //     "UCACGAACCUGCGCCAGCGCACUACAGCUUGACUUUGGGAACUGAGUGGCGUAGUGACGG");
    // s1 = const_cast<char*>(
    //     "...((...(((.(((((((.....))))))).((((((((....((((((.....))))))..)))))))).)))...)).........("
    //     "((((.....(((((((((........)))..((((.(....).)))))))))))))))..");
    // s2 = const_cast<char*>(
    //     "...((...((..(((.......(((((...((((((((((..((.(((((.....))))))).))))))))....))...)))))....."
    //     "...)))....)).))..((((((((.(((..((((........))))))))))))).)).");

    // // # recursive 90 nt example
    // seq = const_cast<char*>(
    //     "CUCCGUUCGGCACAGUGGGAUUCAGACUUCCUGCCGCUGGGAGAAACGCGCGGUUCGGGGUGUAAUCAUGGUUUAAGCCUCCGAAGCCC"
    //     "A");
    // s1 = const_cast<char*>(
    //     "((((...((((...((((((........)))))).))))))))........((((((((((..(((....)))...))).))))).))."
    //     ".");
    // s2 = const_cast<char*>(
    //     "((((...(((((.....(((........))))))))...))))........(((((((((..(.(........).)..))))).))))."
    //     ".");

    // fc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);
    auto start = std::chrono::high_resolution_clock::now();

    // findpath fp(seq, s1, s2);

    auto                          finish   = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = finish - start;

    // fmt::print("init: {}\n{}\n{}\n", fp.seq, fp.s1, fp.s2);

    // maxkeep     = 128;

    fc1         = vrna_fold_compound(seq, &md1, VRNA_OPTION_EVAL_ONLY);
    short* pt_1 = vrna_ptable(s1);
    short* pt_2 = vrna_ptable(s2);

    fc2           = vrna_fold_compound(seq, &md2, VRNA_OPTION_EVAL_ONLY);
    short* pt_1_2 = vrna_ptable(s1);
    short* pt_2_2 = vrna_ptable(s2);

    start = std::chrono::high_resolution_clock::now();

    // sorted_paths result;
    std::vector<sorted_path> result = custom_findpath_method(fc1, pt_1, pt_2, maxkeep, INT_MAX - 1, false);



    // std::future<std::vector<sorted_path>> ret1 = std::async(
    //     std::launch::async, custom_findpath_method, fc1, pt_1, pt_2, maxkeep, INT_MAX - 1, false);
    // std::future<std::vector<sorted_path>> ret2 = std::async(
    //     std::launch::async, custom_findpath_method, fc2, pt_1_2, pt_2_2, maxkeep, INT_MAX - 1, true);
    // std::vector<sorted_path> result = ret1.get();
    // std::vector<sorted_path> result2 = ret2.get();

    // // combine vectors
    // std::move(result2.begin(), result2.end(), std::back_inserter(result));



    // best path to [0] (lowest max_en)
    sort(result.begin(), result.end(), 
        [](const auto & a, const auto & b) -> bool
    { 
        return a.max_en < b.max_en; 
    });

    // for (const auto res : result){
    //     IC(res.max_en);
    //     print_moves(res, fc1, s1, false);
    // }


    // int TIMES = 2;
    // std::future<sorted_paths> res1;
    // for (int i = 0; i < TIMES; ++i) { 
    //     res1 = std::async(std::launch::async, &custom_findpath_method, fc2, pt_1_2, pt_2_2, maxkeep, INT_MAX - 1);        
    //     fc1 = fc2;
    //     pt_1 = pt_1_2;
    //     pt_2 = pt_2_2;
    //      }

    // std::thread t1(custom_findpath_method, fc1, pt_1, pt_2, maxkeep, INT_MAX - 1);
    // std::thread t2(custom_findpath_method, fc2, pt_2_2, pt_1_2, maxkeep, INT_MAX - 1);
    // t1.join();
    // t2.join();



    // result;

    // std::vector<sorted_path> current_paths2 = ret2.get();

    finish                                 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = finish - start;

    // auto number_of_results = result.size();

    // print_moves(result[number_of_results - 1], fc1, s1, false);
    print_moves(result[0], fc1, s1, false);

    // cout << result[0].max_en;

    std::chrono::duration<double> elapsed3 = finish - start;

    // fmt::print("done 1: {}\n", elapsed1.count());
    // fmt::print("done 2: {}\n", elapsed2.count());

    // fmt::print("done 3: {}\n", elapsed3.count());

    // find_moves
}