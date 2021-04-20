
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <future>
#include <thread>

#define FMT_HEADER_ONLY
#include <fmt/color.h>
#include <fmt/format.h>

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

#include "ViennaRNA/mfe.h"

#include "ViennaRNA/model.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/structures.h"
}

// various data structures
struct sorted_move;
struct sorted_path;
struct merge_path;
struct int_loops;

// exterior & interior loop related functions
auto bp_distance(short* pt_a, short* pt_b, int min_pos, int max_pos) -> std::tuple<int, int>;
auto bp_distance(short* pt_a, short* pt_b) -> int;

// ----

struct sorted_move {
    int  i;
    int  j;
    int  E;
    bool operator==(const sorted_move& other) const { return (i == other.i && j == other.j); }
};

struct sorted_path {
    std::vector<sorted_move> moves;
    int32_t                  max_en;
    // struct initializer, preallocate vector
    sorted_path(int init_size) { moves.resize(init_size); }
};

struct merge_path {
    std::vector<sorted_move> moves;
    int32_t                  current_s;
    int32_t                  current_en;
    short*                   current_ptable;

    int current_G1_bp_dist;
    int current_G2_bp_dist;

    int current_G1_node;
    int current_G2_node;

    // struct initializer, preallocate vector
    // sorted_path(int i2) {
    //     moves.resize(i2);
    //     }

    bool operator==(const merge_path& current) const
    {
        // this vs. current
        if (memcmp(current.current_ptable, this->current_ptable,
                   sizeof(short) * (this->current_ptable[0] + 1))) {
            return false;
        }
        return true;
    }
};

// interior loop data structure
struct int_loops {
    int                    start;
    int                    end;
    int                    bp_dist;
    std::vector<int_loops> nested_sections{}; // recursion

    // constructors
    int_loops(int a, int b)
    {
        start = a;
        end   = b;
    }
    int_loops(int a, int b, int c)
    {
        start   = a;
        end     = b;
        bp_dist = c;
    }
};

// print funcion for interior loop sections
auto operator<<(std::ostream& os, const int_loops& s) -> std::ostream&
{
    os << '[' << s.start << ", ";
    for (const auto& r : s.nested_sections) {
        os << r;  // recursive ostream call
        os << ", ";
    }
    os << s.end << ']';
    return os;
}

// exterior loops are just multiple interior loops
using ext_loops = std::vector<int_loops>;

// print funcion for exterior loops
auto operator<<(std::ostream& os, const ext_loops& s) -> std::ostream&
{
    int i = 0;
    for (const auto& int_loop : s) {
        i++;
        os << int_loop;
        if (i != s.size()) { os << ", "; }
    }
    os << "\n";
    return os;
}

auto find_interior_loops(short* pt_1, short* pt_2, int min_pos, int max_pos)
    -> std::vector<int_loops>
{
    // todo: this should be done much simpler...

    std::vector<std::tuple<int, int, int, int, float>> candidates;

    // go from min to max pos, find all nested sections and respective basepair distances
    for (int i = min_pos + 1; i < max_pos - 1; i++) {
        const int j1 = pt_1[i];
        const int j2 = pt_2[i];

        if (i < min_pos or i + 1 > max_pos) continue;
        if (j1 == 0 and j2 == 0) continue;
        if (j1 != j2 or i > j1) continue;

        const int j = j1;

        // inner / outer basepair distance (between ptable1 and 2 with i and j as limits)
        auto [inner_bp_dist, outer_bp_dist] = bp_distance(pt_1, pt_2, i, j);
        // outer_bp_dist = bp_distance(outer_p_table1, outer_p_table2);

        if (std::min(inner_bp_dist, outer_bp_dist) < 1) continue;

        const float inner_size = j - i;
        const float outer_size = max_pos - min_pos - inner_size;
        const float optimize   = std::abs(0.7 - (inner_size / outer_size));

        // IC(i, j, inner_bp_dist, outer_bp_dist, optimize);
        candidates.emplace_back(i, j, inner_bp_dist, outer_bp_dist, optimize);
    }

    // minimize the outer basepair distance - this maximizes recursive potential
    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
        if (std::get<3>(a) == std::get<3>(b)) return std::get<4>(a) < std::get<4>(b);
        return std::get<3>(a) < std::get<3>(b);
    });

    std::vector<int_loops> nested_sections{};
    int                    bp_dist_left = 99999;

    // for (const auto& candidate : candidates) {
    for (const auto& [i, j, inner_bp_dist, outer_bp_dist, optimize] : candidates) {
        // IC(nested_sections, i, j, bp_dist_left);

        if (bp_dist_left - inner_bp_dist < 2) continue;  // nothing left

        // section i to j can't be part of the already existing section
        // dont add something twice
        bool ignore = false;
        for (const auto& section : nested_sections) {
            if ((i > section.start and j < section.end) or (i < section.start and j > section.end))
                ignore = true;
        }
        if (ignore) continue;

        if (bp_dist_left == 99999)
            bp_dist_left = outer_bp_dist;
        else
            bp_dist_left -= inner_bp_dist;

        int_loops current_section(i, j, inner_bp_dist);

        // recursive call
        current_section.nested_sections = find_interior_loops(pt_1, pt_2, i + 1, j - 1);
        nested_sections.push_back(current_section);

        // IC(i, j);
    }

    return nested_sections;
}

auto find_exterior_loops(short* pt_1, short* pt_2) -> std::vector<int_loops>
{
    ext_loops ext_loops{};
    int       i = 1;

    // generate loop table like indices, but with redundancy
    int current_l1 = 0;
    int current_l2 = 0;

    // max loop level
    int max_l1 = 0;
    int max_l2 = 0;

    // const int len = pt_1[0];
    // for (int i=1, i<=len+1, i++){

    for (int j = 1; j <= pt_1[0]; j++) {
        // fmt::print ("pt: {} {} \n", j, pt_1[j]);

        if (pt_1[j] > j) {  // opening bracket
            current_l1 += 1;
            max_l1 = current_l1;
        } else if (pt_1[j] > 0 and pt_1[j] < j) {
            current_l1 -= 1;
        }
        if (pt_2[j] > j) {  // opening bracket
            current_l2 += 1;
            max_l2 = current_l2;
        } else if (pt_2[j] > 0 and pt_2[j] < j) {
            current_l2 -= 1;
        }

        if ((max_l1 != 0 and max_l2 != 0) and (current_l1 < 1 and current_l2 < 1)) {
            
            // append and start a new ext loop
            
            max_l1 = 0;
            max_l2 = 0;
            const auto [outer_bp_dist, _] = bp_distance(pt_1, pt_2, i, j);

            // recursively find all interior loops
            int_loops current_section(i, j, outer_bp_dist);
            if (outer_bp_dist > 0) {
                current_section.nested_sections = find_interior_loops(pt_1, pt_2, i + 1, j - 1);
            }
            ext_loops.emplace_back(current_section);

            // std::cout << "size" << i << " " << j << " " << outer_bp_dist << "\n";

            i = j + 1;
        }
    }

    // std::cout << "size" << ext_loops.size() << "\n";
    // if (ext_loops.size() > 0) { std::cout << "size" << ext_loops << "\n"; }

    return ext_loops;
}

auto bp_distance(short* pt_a, short* pt_b, int min_pos, int max_pos) -> std::tuple<int, int>
{
    /*
    basepair distance between 2 pairing tables, see vrna_bp_distance(), with boundaries
     */
    int inner_bp_dist = 0;
    int outer_bp_dist = 0;
    int length        = pt_a[0];
    for (int i = 1; i <= length; i++)
        if (pt_a[i] != pt_b[i]) {
            if (pt_a[i] > i) {
                if (i > min_pos and i < max_pos)
                    inner_bp_dist++;
                else
                    outer_bp_dist++;
            }
            // both can happen at once...
            if (pt_b[i] > i) {
                if (i > min_pos and i < max_pos)
                    inner_bp_dist++;
                else
                    outer_bp_dist++;
            }
        }
    return {inner_bp_dist, outer_bp_dist};
}

auto bp_distance(short* pt_a, short* pt_b) -> int
{
    /*
    basepair distance between 2 pairing tables, see vrna_bp_distance()
     */
    int bp_dist = 0;
    int length  = pt_a[0];
    for (int i = 1; i <= length; i++)
        if (pt_a[i] != pt_b[i]) {
            if (pt_a[i] > i) { bp_dist++; }
            // both can happen at once...
            if (pt_b[i] > i) { bp_dist++; }
        }
    return bp_dist;
}

auto merge_pairing_table(short* pt_1, short* pt_2) -> short*
{
    // if pt_2 has base pairs which are not present in pt_1, add them
    for (int i = 1; i <= pt_1[0]; i++) {
        if (pt_1[i] == 0 and pt_2 != 0) { pt_1[i] = pt_2[i]; }
    }
    return pt_1;
}