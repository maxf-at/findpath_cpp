
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

// #include <execution>

#include "./icecream.hpp"
using icecream::ic;
using namespace std;

#include "findpath_B.hpp"

#include "s_graph.hpp"

#include "findpath_header.hpp"

int compare_ptable(const void* A, const void* B)
{
    merge_path *a, *b;
    int         c;

    a = (merge_path*)A;
    b = (merge_path*)B;

    c = memcmp(a->current_ptable, b->current_ptable, a->current_ptable[0] * sizeof(short));
    if (c != 0) return c;

    if ((a->current_s - b->current_s) != 0) return a->current_s - b->current_s;

    return a->current_en - b->current_en;
}

PRIVATE int compare_energy(const void* A, const void* B)
{
    merge_path *a, *b;

    a = (merge_path*)A;
    b = (merge_path*)B;

    if ((a->current_s - b->current_s) != 0) return a->current_s - b->current_s;

    return a->current_en - b->current_en;
}

// auto print_moves(sorted_path& path, vrna_fold_compound_t* fc, const char* s1)
auto print_moves(const auto& path, vrna_fold_compound_t* fc, const char* s1)
{
    auto const& moves = path.moves;

    // std::format("{} {}!", "Hello", "world", "something"); // OK, produces "Hello world!"

    short* pt;
    pt = vrna_ptable(s1);

    std::cout << fmt::format("Hello!\n");

    float en     = vrna_eval_structure_pt(fc, pt) / 100.0;
    float max_en = -INT_MAX;

    fmt::print("{} {:7.2f} ({:4}/{:4})\n", s1, en, 0, 0);

    for (auto const& move : moves) {
        std::string insert1, insert2;

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

        // std::string str2 = str.substr (3,5);

        const char* s = vrna_db_from_ptable(pt);
        en            = vrna_eval_structure_pt(fc, pt) / 100.0;
        if (en > max_en) { max_en = en; }

        fmt::print("{} {:7.2f} ({:4}/{:4})\n", s, en, move.i, move.j);

        // printf("%.*s", std::abs(move.i), s);
        // fmt::print("{}", insert1);
        // // printf("%.*s", std::abs(move.i-move.j), s + std::abs(move.i-move.j));
        // printf("%.*s", std::abs(move.i-move.j-3), s + std::abs(move.i) + 1);
        // fmt::print("{}", insert2);
        // // printf("%.*s", std::abs(move.i-move.j-3), s + std::abs(move.j) + 1);
        // printf("%s", s + strlen(s) - std::abs(move.i));
    }

    fmt::print("S: {:6.2f} kcal/mol\n", max_en);
}

auto available_edges(const auto& G1_node, const auto& G2_node, int G1_G2, auto& next_paths,
                     const auto& current_path, vrna_fold_compound_t* fc, int max_en)
{
    for (const auto& current_edge : G1_node.out_edges) {
        const auto& i = current_edge.i;
        const auto& j = current_edge.j;
        auto        current_en =
            current_path.current_en + vrna_eval_move_pt(fc, current_path.current_ptable, i, j);

        if (current_en <= max_en) {
            // generate new pairing tables
            // auto current_ptable = vrna_ptable_copy(pt_1);
            merge_path new_path;
            if (G1_G2 == 1) {
                new_path.current_G1_node = current_edge.destination;
                new_path.current_G2_node = current_path.current_G2_node;
            } else {  // reversed
                new_path.current_G1_node = current_path.current_G1_node;
                new_path.current_G2_node = current_edge.destination;
            }

            new_path.current_en     = current_en;
            new_path.current_s      = std::max(current_en, current_path.current_s);
            new_path.current_ptable = vrna_ptable_copy(current_path.current_ptable);
            new_path.moves          = current_path.moves;
            new_path.moves.push_back({i, j, current_en});

            if (j < 0) {  // delete a basepair
                new_path.current_ptable[-i] = 0;
                new_path.current_ptable[-j] = 0;
            } else {  // add a basepair
                new_path.current_ptable[i] = j;
                new_path.current_ptable[j] = i;
            }

            next_paths.push_back(new_path);
            // auto en = vrna_eval_structure_pt(fc, new_path.current_ptable);
            // std::string s = vrna_db_from_ptable(new_path.current_ptable);
            // fmt::print("{} d={} G1 {} {} {} {}\n", s, d, i, j, en, current_en);
        }
    }
}

auto merge_once(auto G1, auto G2, short* pt_1, int s1_en, short* pt_2, int s2_en, bool direction,
                int total_bp_dist, int max_en, int merge_search_width, vrna_fold_compound_t* fc)
{
    // this is the start of findpath_once, assuming we have maxE etc...

    int current_i_node = 0;
    int current_j_node = 0;

    if (!direction) {
        current_i_node = G1.bp_dist;
        current_i_node = G2.bp_dist;
    }

    // init single empty path
    std::vector<merge_path> all_paths;

    merge_path current_path;
    current_path.current_G1_node = 0;
    current_path.current_G2_node = 0;
    current_path.current_en      = s1_en;
    current_path.current_s       = s1_en;
    current_path.current_ptable  = pt_1;

    all_paths.push_back(current_path);

    const auto& G1_node = G1.node_list[current_path.current_G1_node];
    const auto& G2_node = G2.node_list[current_path.current_G2_node];

    if (direction) { const auto& edges = G1_node.out_edges; }

    // fmt::print("moves {} {}\n", G1_node.out_edges[0].i, G1_node.out_edges[0].j);
    // fmt::print("moves0 {} {}\n", G1.node_list[0].out_edges[0].i, G1.node_list[0].out_edges[0].j);
    // fmt::print("moves1 {} {}\n", G1.node_list[1].out_edges[0].i, G1.node_list[1].out_edges[0].j);
    // fmt::print("moves2 {} {}\n", G1.node_list[2].out_edges[0].i, G1.node_list[2].out_edges[0].j);
    // fmt::print("moves3 {} {}\n", G1.node_list[3].out_edges[0].i, G1.node_list[3].out_edges[0].j);

    // G1 moves
    fmt::print("start\n");

    // total_bp_dist = 4;
    for (int d = 1; d <= total_bp_dist; d++) {
        std::vector<merge_path> next_paths;

        for (const auto current_path : all_paths) {
            const auto& G1_node = G1.node_list[current_path.current_G1_node];
            const auto& G2_node = G2.node_list[current_path.current_G2_node];
            // fill up next_paths - either move 1 step on G1 or G2
            available_edges(G1_node, G2_node, 1, next_paths, current_path, fc, max_en);
            available_edges(G2_node, G1_node, 2, next_paths, current_path, fc, max_en);
        }
        // fmt::print("size before sort: {} \n", next_paths.size());

        if (next_paths.size() == 0) {
            std::vector<merge_path> empty_paths;
            return empty_paths;
        }

        std::qsort(&next_paths[0], next_paths.size(), sizeof(merge_path), compare_ptable);

        // for (const auto& next_path : next_paths) {
        //     for (const auto& move : next_path.moves) { fmt::print("({} / {}) ", move.i, move.j);
        //     } fmt::print(" d={}, sorted: {} {}\n", d, next_path.current_en, next_path.current_s);
        // }

        // compare current with next pairing table...

        // optimize findpath: this sorting / try moves without ptable creation
        // this is probably better than in findpath? does not guarantee dropping...

        int len = next_paths[0].current_ptable[0];

        if (d != total_bp_dist) {
            for (int u = 0, c = 1; c < next_paths.size(); c++) {
                if (memcmp(next_paths[u].current_ptable, next_paths[c].current_ptable,
                           sizeof(short) * len) != 0) {
                    u++;
                    // next_paths[u] = next_paths[c];
                } else {
                    // next_paths[c].current_s = INT_MAX;
                    // next_paths[c].current_en = INT_MAX;
                    next_paths[c].current_s  = 999;
                    next_paths[c].current_en = 999;
                    // fmt::print("drop {}, {} \n", u, c);

                    // free_intermediate(next + c);
                }
            }
        }

        std::qsort(&next_paths[0], next_paths.size(), sizeof(merge_path), compare_energy);

        if (next_paths.size() > merge_search_width) { next_paths.resize(merge_search_width); }

        // last iteration
        // if (d == total_bp_dist) {
        //     int i = 0;
        //     for (const auto& next_path : next_paths) {
        //         ++i;
        //         if (i > 10) break;
        //         for (const auto& move : next_path.moves) {
        //             fmt::print("({} / {}) ", move.i, move.j);
        //         }
        //         fmt::print(" d={}, sorted: {} {}\n", d, next_path.current_en,
        //         next_path.current_s);
        //     }
        // }

        all_paths = next_paths;
    }

    return all_paths;
}

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

    vrna_fold_compound_t* fc = nullptr;
    vrna_md_t             md;
    // set model params
    set_model_details(&md);

    // fc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);
    sorted_paths result;

    // result = custom_findpath_method(fc, s1, s2, maxkeep, INT_MAX - 1);
    // cout << result[0].max_en;

    // init graph

    seq = const_cast<char*>(
        "AGUUCAGCUACUUACACUGCGGCCUACGGGGGAAUCACUGACCGCGUUAUACCCUAUAUUCGUAUAUGCACGUACAAUACCAACAGUAGG"
        "AAUACUUUCUGGGGAGUAGAGUAGAAAUCUUGGCCCCCGACUUGACUGUUAAUCCUCUAC");
    s1 = const_cast<char*>(
        ".((...(((...........)))..))(((((..(((......(((.(((((.........))))))))...(((..((((..(((.((("
        "....))).)))..).)))..))).......))))))))(((......)))..........");
    s2 = const_cast<char*>(
        "...........(((....(((((...((((((..(((......((.((((.(((((.((((.((..((...(((...)))...)).)).)"
        ")))......))))).)))).))........)))))))))....).)))))))........");

    // ex2
    seq = const_cast<char*>(
        "AGGGGAAAGUUAUCGAACAGACCUUGUUUGGCUUAUGGCACUUGGGUAGAUUUUCUCUACCGUUGUCAUGAAAACGCACCGCAAGAGAAG"
        "UCACGAACCUGCGCCAGCGCACUACAGCUUGACUUUGGGAACUGAGUGGCGUAGUGACGG");
    s1 = const_cast<char*>(
        "...((...(((.(((((((.....))))))).((((((((....((((((.....))))))..)))))))).)))...)).........("
        "((((.....(((((((((........)))..((((.(....).)))))))))))))))..");
    s2 = const_cast<char*>(
        "...((...((..(((.......(((((...((((((((((..((.(((((.....))))))).))))))))....))...)))))....."
        "...)))....)).))..((((((((.(((..((((........))))))))))))).)).");

    // # recursive 90 nt example
    seq = const_cast<char*>(
        "CUCCGUUCGGCACAGUGGGAUUCAGACUUCCUGCCGCUGGGAGAAACGCGCGGUUCGGGGUGUAAUCAUGGUUUAAGCCUCCGAAGCCC"
        "A");
    s1 = const_cast<char*>(
        "((((...((((...((((((........)))))).))))))))........((((((((((..(((....)))...))).))))).))."
        ".");
    s2 = const_cast<char*>(
        "((((...(((((.....(((........))))))))...))))........(((((((((..(.(........).)..))))).))))."
        ".");

    // fc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);
    findpath fp(seq, s1, s2);

    // fmt::print("init: {}\n{}\n{}\n", fp.seq, fp.s1, fp.s2);

    maxkeep = 128;

    // seq = "GCCGUGUGCUCUGGCUCUAGACUCUAACACAGGGCGGUUGUCCCUAUUACUGGUACUCGU";
    // s1  = "(((.((((....(........).....)))).))).........................";
    // s2  = "(((.((((..(((....))).......)))).))).........................";
    seq = const_cast<char*>(
        "AGUUCAGCUACUUACACUGCGGCCUACGGGGGAAUCACUGACCGCGUUAUACCCUAUAUUCGUAUAUGCACGUACAAUACCAACAGUAGG"
        "AAUACUUUCUGGGGAGUAGAGUAGAAAUCUUGGCCCCCGACUUGACUGUUAAUCCUCUAC");
    s1 = const_cast<char*>(
        ".((...(((...........)))..))(((((..((....................................................."
        "................................)))))))(((......)))..........");
    s2 = const_cast<char*>(
        "...........(((....(((((...((((((..((....................................................."
        "................................))))))))....).)))))))........");
    fc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);

    auto start = std::chrono::high_resolution_clock::now();
    result     = custom_findpath_method(fc, s1, s2, maxkeep, INT_MAX - 1);
    cout << result[0].max_en;

    // print_moves(result[0], fc, s1);
    // print_moves(result[1], fc, s1);
    // print_moves(result[2], fc, s1);
    // print_moves(result[3], fc, s1);

    s_graph G1;
    G1.fc      = fc;
    G1.pt_1    = vrna_ptable(s1);
    G1.pt_2    = vrna_ptable(s2);
    G1.bp_dist = vrna_bp_distance(s1, s2);
    G1.add_paths(result);

    auto                          finish   = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = finish - start;

    // G1.info();

    // s1     = ".......................................((.((.(....))).))....";
    // s2     = "...................................((.((.((........))))))...";
    s1 = const_cast<char*>(
        "....................................(......(((.(((((.........))))))))...(((..((((..(((.(("
        "(....))).)))..).)))..))).......).............................");
    s2 = const_cast<char*>(
        "....................................(......((.((((.(((((.((((.((..((...(((...)))...)).))."
        "))))......))))).)))).))........).............................");

    start = std::chrono::high_resolution_clock::now();
    // fc     = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);
    result = custom_findpath_method(fc, s1, s2, maxkeep, INT_MAX - 1);
    cout << result[0].max_en;

    // print_moves(result[0], fc, s1);
    // print_moves(result[1], fc, s1);
    // print_moves(result[2], fc, s1);
    // print_moves(result[3], fc, s1);

    s_graph G2;
    G2.fc      = fc;
    G2.pt_1    = vrna_ptable(s1);
    G2.pt_2    = vrna_ptable(s2);
    G2.bp_dist = vrna_bp_distance(s1, s2);
    G2.add_paths(result);

    finish                                 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = finish - start;
    // G2.info();

    // print_moves(result[0], fc, s1);
    // print_moves(result[1], fc, s1);
    // print_moves(result[3], fc, s1);

    start = std::chrono::high_resolution_clock::now();

    // step 1 merge pairing tables (maybe don't do this here...)

    short* merged_pt_1 = vrna_ptable_copy(G1.pt_1);
    short* merged_pt_2 = vrna_ptable_copy(G1.pt_2);

    for (int i = 1; i < merged_pt_1[0]; i++) {
        // fmt::print("i {}\n", i);
        if ((merged_pt_1[i] == 0 and G2.pt_1 != 0) or (merged_pt_1[i] != 0 and G2.pt_1 == 0)) {
            merged_pt_1[i] = G2.pt_1[i];
        }
        if ((merged_pt_2[i] == 0 and G2.pt_2 != 0) or (merged_pt_2[i] != 0 and G2.pt_2 == 0)) {
            merged_pt_2[i] = G2.pt_2[i];
        }
    }

    int s1_en = vrna_eval_structure_pt(fc, merged_pt_1);
    int s2_en = vrna_eval_structure_pt(fc, merged_pt_2);

    fmt::print("{} {} {}\n", vrna_db_from_ptable(merged_pt_1), G1.bp_dist, s1_en);
    fmt::print("{} {} {}\n", vrna_db_from_ptable(merged_pt_2), G2.bp_dist, s2_en);

    int total_bp_dist = G1.bp_dist + G2.bp_dist;

    int max_en = INT_MAX;
    // int  max_en             = -1390;
    int  merge_search_width = 100;
    bool direction          = true;

    auto all_paths = merge_once(G1, G2, merged_pt_1, s1_en, merged_pt_2, s2_en, direction,
                                total_bp_dist, max_en, merge_search_width, fc);

    // for (const auto& move : all_paths[0].moves) { fmt::print("({} / {}) ", move.i, move.j); }
    // fmt::print("sorted: {} {}\n", all_paths[0].current_en, all_paths[0].current_s);
    // s1 =
    // ".((...(((...........)))..))(((((..(((......(((.(((((.........))))))))...(((..((((..(((.(("
    //      "(....))).)))..).)))..))).......))))))))(((......)))..........";
    // print_moves(all_paths[0], fc, s1);

    // for (const auto &current_edge : G2_node.out_edges) {
    //     const auto &i = current_edge.i;
    //     const auto &j = current_edge.j;
    //     auto en = current_en_1 + vrna_eval_move_pt(fc, pt_1, i, j);
    //     fmt::print("{} {} {}\n", i,j, en);
    // }

    finish                                 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed3 = finish - start;

    fmt::print("done 1: {}\n", elapsed1.count());
    fmt::print("done 2: {}\n", elapsed2.count());
    fmt::print("done 3: {}\n", elapsed3.count());

    // find_moves
}