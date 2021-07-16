
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <deque>
#include <set>

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
auto bp_distance(std::vector<short> pt_a, std::vector<short> pt_b, int min_pos, int max_pos) -> std::tuple<int, int>;
auto bp_distance(std::vector<short> pt_a, std::vector<short> pt_b) -> int;

// ----

struct sorted_move {
    int  i;
    int  j;
    int  E;
    bool operator==(const sorted_move& other) const { return (i == other.i && j == other.j); }
};

struct sorted_path {
    std::vector<sorted_move> moves;
    int                      max_en;
    std::string              destination;

    // struct initializer, preallocate vector
    sorted_path(int init_size) { moves.resize(init_size); }
    sorted_path(int init_size, int max_en) : max_en{max_en} { moves.resize(init_size); }
    sorted_path(int init_size, int max_en, std::string destination)
            : max_en{max_en}, destination{destination}
    {
        moves.resize(init_size);
    }
};

struct merge_path {
    std::vector<sorted_move> moves;
    sorted_move*             move_ptr;

    int    current_s;
    int    current_en;
    short* current_ptable;

    size_t s_hash;

    int current_G1_bp_dist;
    int current_G2_bp_dist;

    int current_G1_node;
    int current_G2_node;

    int last_index;
    // int last_G1_node;
    // int last_G2_node;
    int i_move;
    int j_move;

    // struct initializer, preallocate vector
    // sorted_path(int i2) {
    //     moves.resize(i2);
    //     }

    bool operator==(const merge_path& current) const
    {
        // this vs. current
        // if (memcmp(current.current_ptable, this->current_ptable,
        //            sizeof(short) * (this->current_ptable[0] + 1))) {
        //     return false;
        // }
        // return true;

        return current.s_hash == this->s_hash;
    }
};

// interior loop data structure
struct int_loops {
    int                    start;
    int                    end;
    int                    bp_dist;
    std::vector<int_loops> nested_sections{};  // recursion

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

auto find_interior_loops(std::vector<short> pt_1, std::vector<short> pt_2, int min_pos, int max_pos)
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

auto find_exterior_loops(std::vector<short> pt_1, std::vector<short> pt_2) -> std::vector<int_loops>
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

            max_l1                        = 0;
            max_l2                        = 0;
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

auto bp_distance(std::vector<short> pt_a, std::vector<short> pt_b, int min_pos, int max_pos) -> std::tuple<int, int>
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

auto bp_distance(std::vector<short> pt_a, std::vector<short> pt_b) -> int
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





// legacy

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




// auto merge_pairing_table(std::vector<short> pt_1, std::vector<short> pt_2) -> std::vector<short>
// {
//     // if pt_2 has base pairs which are not present in pt_1, add them
//     for (int i = 1; i <= pt_1[0]; i++) {
//         if (pt_1[i] == 0 and pt_2 != 0) { pt_1[i] = pt_2[i]; }
//     }
//     return pt_1;
// }


auto merge_pairing_table(short* pt_1, short* pt_2) -> short*
{
    // if pt_2 has base pairs which are not present in pt_1, add them
    for (int i = 1; i <= pt_1[0]; i++) {
        if (pt_1[i] == 0 and pt_2 != 0) { pt_1[i] = pt_2[i]; }
    }
    return pt_1;
}

auto print_moves(const auto& path, vrna_fold_compound_t* fc, const short* pt1,
                 bool show_path = true)
{
    // this print path function is unused (?)
    auto const& moves = path.moves;

    short* pt = vrna_ptable_copy(pt1);

    char* s1 = vrna_db_from_ptable(pt);

    // short* pt;
    // pt = vrna_ptable(s1);

    float en     = vrna_eval_structure_pt(fc, pt) / 100.0;
    float max_en = float(-INT_MAX);

    if (show_path) fmt::print("{} {:7.2f} ({:4}/{:4})\n", s1, en, 0, 0);

    for (auto const& move : moves) {
        std::string insert1, insert2;

        if (move.i == 0) { continue; }  // unused indirect move

        if (move.j < 0) {
            /*it's a delete move */
            pt[-move.i] = 0;
            pt[-move.j] = 0;
            insert1 = insert2 = fmt::format(fmt::emphasis::bold | fg(fmt::color::red), ".");

        } else {
            pt[move.i] = move.j;
            pt[move.j] = move.i;
            insert1 = insert2 = "#";
        }

        const char* s = vrna_db_from_ptable(pt);
        en            = vrna_eval_structure_pt(fc, pt) / 100.0;
        if (en > max_en) { max_en = en; }

        if (show_path) fmt::print("{} {:7.2f} ({:4}/{:4})\n", s, en, move.i, move.j);
    }

    fmt::print("S: {:6.2f} kcal/mol\n", max_en);

    free(pt);
    free(s1);
}

// internal function structure_utils.c, unchanged
int extract_pairs(short* pt, const char* structure, const char* pair)
{
    const char*  ptr;
    char         open, close;
    short*       stack;
    unsigned int i, j, n;
    int          hx;

    n     = (unsigned int)pt[0];
    stack = (short*)vrna_alloc(sizeof(short) * (n + 1));

    open  = pair[0];
    close = pair[1];

    for (hx = 0, i = 1, ptr = structure; (i <= n) && (*ptr != '\0'); ptr++, i++) {
        if (*ptr == open) {
            stack[hx++] = i;
        } else if (*ptr == close) {
            j = stack[--hx];

            if (hx < 0) {
                vrna_message_warning(
                    "%s\nunbalanced brackets '%2s' found while extracting base pairs", structure,
                    pair);
                free(stack);
                return 0;
            }

            pt[i] = j;
            pt[j] = i;
        }
    }

    free(stack);

    if (hx != 0) {
        vrna_message_warning("%s\nunbalanced brackets '%2s' found while extracting base pairs",
                             structure, pair);
        return 0;
    }

    return 1; /* success */
}

// pairing table as std::vector
auto ptable_from_string(std::string s) -> std::vector<short>
{
    char         pairs[3];
    unsigned int i, n;
    unsigned int options = VRNA_BRACKETS_RND;

    n = s.length();
    std::vector<short> pt(n+2);

    if (n > SHRT_MAX) {
        vrna_message_warning(
            "vrna_ptable_from_string: "
            "Structure too long to be converted to pair table (n=%d, max=%d)",
            n, SHRT_MAX);
        return pt;
    }

    // pt    = (short*)vrna_alloc(sizeof(short) * (n + 2));
    pt[0] = (short)n;

    if ((options & VRNA_BRACKETS_RND) && (!extract_pairs(pt.data(), s.c_str(), "()"))) {
        return pt;
    }

    if ((options & VRNA_BRACKETS_ANG) && (!extract_pairs(pt.data(), s.c_str(), "<>"))) {
        return pt;
    }

    if ((options & VRNA_BRACKETS_CLY) && (!extract_pairs(pt.data(), s.c_str(), "{}"))) {
        return pt;
    }

    if ((options & VRNA_BRACKETS_SQR) && (!extract_pairs(pt.data(), s.c_str(), "[]"))) {
        return pt;
    }

    if (options & VRNA_BRACKETS_ALPHA) {
        for (i = 65; i < 91; i++) {
            pairs[0] = (char)i;
            pairs[1] = (char)(i + 32);
            pairs[2] = '\0';
            if (!extract_pairs(pt.data(), s.c_str(), pairs)) {
                return pt;
            }
        }
    }

    return pt;
}

// Robert Jenkins' 32 bit integer hash function
uint32_t int_hash(uint32_t a)
{
    a = (a + 0x7ed55d16) + (a << 12);
    a = (a ^ 0xc761c23c) ^ (a >> 19);
    a = (a + 0x165667b1) + (a << 5);
    a = (a + 0xd3a2646c) ^ (a << 9);
    a = (a + 0xfd7046c5) + (a << 3);
    a = (a ^ 0xb55a4f09) ^ (a >> 16);
    return a;
}

/*  Written in 2015 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>

/* This is a fixed-increment version of Java 8's SplittableRandom generator
   See http://dx.doi.org/10.1145/2714064.2660195 and
   http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

   It is a very fast generator passing BigCrush, and it can be useful if
   for some reason you absolutely want 64 bits of state. */

// static uint64_t x; /* The state can be seeded with any value. */

uint64_t int_hash_64(uint64_t x)
{
    uint64_t z = (x += 0x9e3779b97f4a7c15);
    z          = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z          = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

#define mix(h)                        \
    ({                                \
        (h) ^= (h) >> 23;             \
        (h) *= 0x2127599bf4325c37ULL; \
        (h) ^= (h) >> 47;             \
    })

uint64_t fasthash64(const void* buf, size_t len, uint64_t seed)
{
    const uint64_t       m   = 0x880355f21e6d1965ULL;
    const uint64_t*      pos = (const uint64_t*)buf;
    const uint64_t*      end = pos + (len / 8);
    const unsigned char* pos2;
    uint64_t             h = seed ^ (len * m);
    uint64_t             v;

    while (pos != end) {
        v = *pos++;
        h ^= mix(v);
        h *= m;
    }

    pos2 = (const unsigned char*)pos;
    v    = 0;

    switch (len & 7) {
        case 7: v ^= (uint64_t)pos2[6] << 48;
        case 6: v ^= (uint64_t)pos2[5] << 40;
        case 5: v ^= (uint64_t)pos2[4] << 32;
        case 4: v ^= (uint64_t)pos2[3] << 24;
        case 3: v ^= (uint64_t)pos2[2] << 16;
        case 2: v ^= (uint64_t)pos2[1] << 8;
        case 1:
            v ^= (uint64_t)pos2[0];
            h ^= mix(v);
            h *= m;
    }

    return mix(h);
}

size_t hash_c_string(const char* p, size_t s)
{
    size_t       result = 0;
    const size_t prime  = 31;
    for (size_t i = 0; i < s; ++i) { result = p[i] + (result * prime); }
    return result;
}

size_t hash_pt(const short* p, short s)
{
    size_t       result = 0;
    const size_t prime  = 31;
    for (short i = 0; i < s; ++i) { result = p[i] + (result * prime); }
    return result;
}

#include <string_view>

static size_t hash_cstr(const char* s, int len)
{
    // C++ 17 magic
    // https://stackoverflow.com/questions/34597260/stdhash-value-on-char-value-and-not-on-memory-address

    // fmt::print ("hash: {} / {} {} {} \n", s, s[0], sizeof(char), std::strlen(s));
    // return std::hash<std::string_view>()(std::string_view(s, s[0]*sizeof(char)));
    // return std::hash<std::string_view>()(std::string_view(s, std::strlen(s))); // does not work
    return std::hash<std::string_view>()(std::string_view(s, len));
}