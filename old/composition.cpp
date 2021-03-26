#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <algorithm>
#include <any>
#include <array>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <vector>

// #include <execution>

#include "./icecream.hpp"

using icecream::ic;
// using namespace std;

// #ifdef HAVE_CONFIG_H
// #include "config.h"
// #endif

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

auto p_to_s(short* p_table) -> std::string
{
    // """ptable --> string"""
    std::string s = "";
    for (int i = 1; i < p_table[0]; i++) {
        if (p_table[i] == 0)
            s += '.';
        else if (p_table[i] > i)
            s += '(';
        else
            s += ')';
    }
    return s;
}

auto sections_to_string(std::vector<std::any>& results) -> std::string
{
    std::string s = "[";
    int         i = 0;
    // this is ... questionable
    for (const auto& r : results) {
        if (r.type() == typeid(int)) s += std::to_string(std::any_cast<int>(r));

        if (r.type() == typeid(std::vector<std::any>)) {
            std::vector<std::any> k = std::any_cast<std::vector<std::any>>(r);
            s += sections_to_string(k);
        }
        if (i + 1 != results.size()) s += ", ";
        i++;
    }
    s += ']';
    return s;
}

typedef struct section {
    int start;
    int end;
    // std::vector<section>;
    std::vector<section> nested_sections{};

    std::ostream& operator<<(std::ostream& out) { return out << '[' << start << " " << end << ']'; }

} sections;

// print funcion for sections
std::ostream& operator<<(std::ostream& os, const section& s)
{
    os << '[' << s.start << ", ";
    for (const auto& r : s.nested_sections) {
        os << r; // recursive call
        os << ", ";
    }
    os << s.end << ']';
    return os;
}

auto main(int argc, char** argv) -> int
{
    printf("main function\n");

    ic.line_wrap_width(990);
    ic.prefix("\033[1;36mDEBUG: \033[0m");
    ic.show_c_string(false);
    // ic.disable();

    std::string input_sequence = "AACGGGUGGGUACUCCCUGGUAAAGCCCGAGUCGAGACAUUGUCAUAUGUAUGAGAUUCCUUUGUUGUUGGUCGGCUGGG";
    std::string input_s1       = "..((((((((....))).......)))))((((((.....(.((((....)))).)...((........))))))))...";
    std::string input_s2       = "...((((...((((....))))..)))).(((((((((....((..........)).......))).....))))))...";

    char* sequence = const_cast<char*>(input_sequence.c_str());
    char* s1       = const_cast<char*>(input_s1.c_str());
    char* s2       = const_cast<char*>(input_s2.c_str());

    short *ptr, *p_table1, *p_table2;
    // move_t *bestpath = NULL;
    // int dir;
    // path_fwd = dir = 0;
    vrna_md_t             md;
    vrna_fold_compound_t* fc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
    p_table1                 = vrna_ptable(s1);
    p_table2                 = vrna_ptable(s2);

    // for (int i = 0; i < input_sequence.length(); i++) { std::cout << pt1[i] << " "; }

    int    len            = p_table1[0];
    short* inner_p_table1 = (short*)vrna_alloc(sizeof(short) * (len + 1));
    short* inner_p_table2 = (short*)vrna_alloc(sizeof(short) * (len + 1));
    short* outer_p_table1 = (short*)vrna_alloc(sizeof(short) * (len + 1));
    short* outer_p_table2 = (short*)vrna_alloc(sizeof(short) * (len + 1));

    char *inner_s1, *inner_s2, *outer_s1, *outer_s2;

    int min_pos = 0;
    int max_pos = p_table1[0];

    std::vector<std::tuple<float, int, int>> candidates;
    // std::vector<std::any>                    results{0};
    sections results;
    results.start = 1;

    for (int i = min_pos + 2; i < max_pos - 1; i++) {
        // iterate over both pair tables, keep 1 nt spacer between start and end
        // compared to strings, everything here has a +1 offset
        int p1 = p_table1[i];
        int p2 = p_table2[i];

        // IC("s", i, p1, p2);
        if (p1 == 0 or p2 == 0) continue;   // some of these continue statements may be redundant
        if (p1 != p2 or i >= p1) continue;  // opening brackets only

        int last_i = i - 1;
        int j      = p1;
        int next_j = p1 + 1;

        // if (p_table1[j + 1] != p_table2[j + 1]) continue;
        if (p_table1[next_j] != p_table2[next_j]) continue;

        // check which compatible sections have the highest potential for recursion

        if (p_table1[last_i] == 0 and p_table2[last_i] != 0) continue;
        if (p_table1[next_j] == 0 and p_table2[next_j] != 0) continue;
        // IC(last_i, i, j, next_j);

        // std::cout << p1 << "/" << p2 << " ";

        // last_i has to be linked to last_j

        // IC(p_table1[p_table1[last_i]], p_table1[next_j + 1]);
        // last_i has to be linked with last_j: )( instead of ()
        // this guarantees that i & j are correctly linked
        if (p_table1[p_table1[last_i]] != p_table1[next_j]) continue;

        // std::cout << "s1: " << s1[last_i - 1] << s1[i - 1] << s1[j - 1] << s1[next_j - 1] << "\n";
        // std::cout << "s2: " << s2[last_i - 1] << s2[i - 1] << s2[j - 1] << s2[next_j - 1] << "\n";

        // memcpy(inner_p_table1, p_table1, (len + 1) * sizeof(short));
        inner_p_table1 = vrna_ptable_copy(p_table1);
        inner_p_table2 = vrna_ptable_copy(p_table2);
        outer_p_table1 = vrna_ptable_copy(p_table1);
        outer_p_table2 = vrna_ptable_copy(p_table2);

        for (int k = 1; k < len; k++) {
            if (k < i or k > j + 1) {
                inner_p_table1[k] = 0;
                inner_p_table2[k] = 0;
            } else {
                outer_p_table1[k] = 0;
                outer_p_table2[k] = 0;
            }
        }

        inner_s1 = vrna_db_from_ptable(inner_p_table1);
        inner_s2 = vrna_db_from_ptable(inner_p_table2);
        outer_s1 = vrna_db_from_ptable(outer_p_table1);
        outer_s2 = vrna_db_from_ptable(outer_p_table2);

        int inner_bp_dist = vrna_bp_distance(inner_s1, inner_s2);
        int outer_bp_dist = vrna_bp_distance(outer_s1, outer_s2);
        if (std::min(inner_bp_dist, outer_bp_dist) < 1) continue;

        const float inner_size = j - i + 1;
        const float outer_size = max_pos - min_pos - inner_size;

        const float optimize = std::abs(0.6 - (inner_size / outer_size));

        IC(last_i, i, j, next_j, inner_size, outer_size, optimize);
        std::cout << inner_s1 << "\n";
        std::cout << inner_s2 << "\n";
        std::cout << outer_s1 << "\n";
        std::cout << outer_s2 << "\n";

        const auto tuple = std::make_tuple(optimize, i, j);
        candidates.push_back(tuple);

        // bp_dist_inner = RNA.bp_distance(inner_s1, inner_s2)
        // bp_dist_outer = RNA.bp_distance(outer_s1, outer_s2)

        // # out of bounds
        // if next_j>len(s1) or last_i < 0:
        //     continue
    }

    // sort candidatess by first column
    std::sort(candidates.begin(), candidates.end(),
              [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });

    std::vector<int> in_use(len, 0);
    ;
    // int i,j;
    // float opt;
    for (const auto& [opt, i, j] : candidates) {
        bool ignore = false;
        // generate in use vector - 1 for stretches where recursion might be viable [00011110001111]
        for (int k = i; k < j; k++) {
            if (in_use[k] != 0) {
                ignore = true;
                break;
            }
        }
        if (ignore) continue;
        for (int k = i; k < j; k++) { in_use[k] = 1; }

        // recursion...

        sections current_result;
        current_result.start = i;
        current_result.end   = j;
        results.nested_sections.push_back(current_result);

        // std::vector<std::any> current_result;
        // current_result.push_back(i);
        // current_result.push_back(j);
        // results.push_back(current_result);

        IC(i, j, opt);
    }

    // results.push_back(len);  // final element is lenght
    results.end = len;

    // for (const auto &r : results)
    //     {
    //         if (r.type() == typeid(int))
    //             std::cout << std::any_cast<int>(r) << "\n";

    //         if (r.type() == typeid(std::vector<std::any>))
    //             std::cout << "recursion" << "\n";
    //     }
    //
    std::cout << results << "\n";

    // for i, (p1, p2) in enumerate(zip(ptables_s1[1:], ptables_s2[1:])):
    //     # if i==0: continue
    //     if i < min_pos or i > max_pos: continue

    //     if p1==p2 and p1 > i:
    //         curr_lvl += 1
    //     elif p1==p2 and p1 < i:
    //         curr_lvl -= 1

    // IC(pt1);

    // const char* s3  = "..............(((((((((....))))))))).........";
    // short*      pt3 = vrna_ptable(s3);
}