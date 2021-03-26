

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
        // std::this_thread::sleep_for(std::chrono::nanoseconds(100000));
    }

    return en;
}

auto system_call(){
    system("./single");
    return 1;

}

auto primes() {

    int i, num = 1, primes = 0;

    while (num <= 100000) { 
        i = 2; 
        while (i <= num) { 
            if(num % i == 0)
                break;
            i++; 
        }
        if (i == num)
            primes++;
        num++;
    }

    return 1;

}

auto main(int argc, char** argv) -> int
{
    std::cout << "entering main function\n";

    auto start = std::chrono::high_resolution_clock::now();

    // int en0 = test_func();
    // int en0 = system_call();
    // int en0 = primes();


    auto                          finish  = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    // std::cout << "1 thread: " << elapsed.count() << "\n";


    // with 2 threads

    start = std::chrono::high_resolution_clock::now();

    std::future<int> ret1 = std::async(std::launch::async, test_func);
    std::future<int> ret2 = std::async(std::launch::async, test_func);
    // std::future<int> ret1 = std::async(std::launch::async, system_call);
    // std::future<int> ret2 = std::async(std::launch::async, system_call);
    // std::future<int> ret1 = std::async(std::launch::async, primes);
    // std::future<int> ret2 = std::async(std::launch::async, primes);

    int en1 = ret1.get();
    int en2 = ret2.get();

    finish  = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << "2 threads: " << elapsed.count() << "\n";


}