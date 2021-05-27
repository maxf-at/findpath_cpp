// this includes the basic elementary findpath algorithm

// #i

#include <pybind11/pybind11.h>
namespace py = pybind11;

// #include "./quadsort/quadsort.hpp"
#include "./sort/timsort.hpp"

// basic findpath structs
struct move_t {
    // short i; /* i,j>0 insert; i,j<0 delete */
    // short j;
    short when; /* 0 if still available, else resulting distance from start */
    // int E;
};
struct move_ij {
    short i; /* i,j>0 insert; i,j<0 delete */
    short j;
};

struct intermediate_t {
    char*  s;
    short* moves;
    short  length;
    int    saddle_en; /**<  @brief  saddle energy so far */
    int    curr_en;   /**<  @brief  current energy */

    int last_id;
    int move_id;  // which move will take next turn

    int move_delete;
    int move_i;
    int move_j;

    size_t s_hash;
};

class single_findpath
{
   private:
    static auto try_moves(vrna_fold_compound_t* vc, intermediate_t c, int maxE,
                          intermediate_t* next, int dist, int bp_dist, short* temp_pt, short* stack)
        -> int;
    static auto find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2,
                               const int current_search_width, int maxE, bool direction,
                               int final_search_width) -> std::vector<sorted_path>;
    static auto findpath_method(vrna_fold_compound_t* vc, short* pt1, short* pt2, int width,
                                int maxE, bool init_direction) -> std::vector<sorted_path>;
    static auto free_intermediate(intermediate_t* i) -> void;
    static auto compare_ptable(const void* A, const void* B) -> int;
    static auto compare_energy(const void* A, const void* B) -> int;

   public:
    auto init(vrna_fold_compound_t* fc, short* pt1, short* pt2, int final_search_width,
              bool mp = true, int en_limit = INT_MAX - 1) -> std::vector<sorted_path>;
    single_findpath(){};
};

inline auto single_findpath::init(vrna_fold_compound_t* fc, short* pt1, short* pt2,
                                  int final_search_width, bool mp, int en_limit)
    -> std::vector<sorted_path>
{
    // fmt::print("{}\n", vrna_db_from_ptable(pt1));
    // fmt::print("{}\n", vrna_db_from_ptable(pt2));

    // mp = false;
    // char* temp_s = vrna_db_from_ptable(pt1);
    // short* temp_pt = vrna_ptable_copy(pt1);
    // encode_pt(temp_pt);

    auto result = findpath_method(fc, vrna_ptable_copy(pt1), vrna_ptable_copy(pt2),
                                  final_search_width, en_limit, mp);

    // best path to [0] (lowest max_en)

    std::sort(result.begin(), result.end(),
              [](const auto& a, const auto& b) -> bool { return a.max_en < b.max_en; });

    // fmt::print("hello fmt\n");
    // fmt::print ("s1 = '{}'\n", vrna_db_from_ptable(pt1));
    // fmt::print ("s2 = '{}'\n", vrna_db_from_ptable(pt2));
    // for (const auto& path : result) {
    //     for (const auto m : path.moves) { fmt::print("({}/{}) ", m.i, m.j); }
    //         fmt::print("| max_en: {}\n", path.max_en);

    // }

    return result;
};

inline auto single_findpath::findpath_method(vrna_fold_compound_t* fc, short* pt1, short* pt2,
                                             int final_search_width, int max_en, bool mp)
    -> std::vector<sorted_path>

{
    short*                   temp_pt;
    std::vector<sorted_path> all_paths;

    // new private bool
    bool direction = true;  // true: s1 -> s2, false: s2 -> s1

    short* pt1_bwd = vrna_ptable_copy(pt1);
    short* pt2_bwd = vrna_ptable_copy(pt2);

    int              last_iteration = final_search_width;
    std::vector<int> iterations{final_search_width};

    // if we have a given energy limit, only take the final search width pass
    // if (max_en != INT_MAX - 1) {
    // construct the list of iterations, e.g. [400, 80, 16]
    float accelerator = 1.0;  // this might not be useful...
    while (last_iteration > 16) {
        const int next_iteration = int(last_iteration / 5.0 * accelerator);
        iterations.push_back(next_iteration);
        last_iteration = next_iteration;
        // accelerator *= accelerator;
        // if (accelerator > 2) accelerator = 2;
    }
    // }

    // iterate back to front [16, 80, 400]
    for (std::vector<int>::reverse_iterator it = iterations.rbegin(); it != iterations.rend();
         ++it) {
        int                      current_search_width = *it;
        int                      saddle_en_1          = max_en;
        int                      saddle_en_2          = max_en;
        std::vector<sorted_path> fwd_paths;
        std::vector<sorted_path> bwd_paths;

        if (mp) {
            std::future<std::vector<sorted_path>> ret1 =
                std::async(std::launch::async, &find_path_once, fc, pt1, pt2, current_search_width,
                           max_en, direction, final_search_width);
            std::future<std::vector<sorted_path>> ret2 =
                std::async(std::launch::async, &find_path_once, fc, pt2_bwd, pt1_bwd,
                           current_search_width, max_en, not direction, final_search_width);
            fwd_paths = ret1.get();
            bwd_paths = ret2.get();
        } else {
            fwd_paths = find_path_once(fc, pt1, pt2, current_search_width, max_en, direction,
                                       final_search_width);
            bwd_paths = find_path_once(fc, pt2_bwd, pt1_bwd, current_search_width, max_en,
                                       not direction, final_search_width);
        }

        if (fwd_paths.size() > 0) {
            // for (const auto m : fwd_paths[0].moves) { fmt::print("{} {} / ", m.i, m.j); }
            std::move(fwd_paths.begin(), fwd_paths.end(),
                      std::back_inserter(all_paths));  // concatenate paths
            saddle_en_1 = fwd_paths[0].max_en;
        }

        if (bwd_paths.size() > 0) {
            // std::cout << "\n";
            // for (const auto m : bwd_paths[0].moves) { fmt::print("{} {} / ", m.i, m.j); }
            std::move(bwd_paths.begin(), bwd_paths.end(),
                      std::back_inserter(all_paths));  // concatenate paths
            saddle_en_2 = bwd_paths[0].max_en;
        }

        // std::cout << "\n" << current_search_width << " / " << saddle_en_1 << " / " <<
        // saddle_en_2
        // << " / "
        //           << max_en << "\n";

        // set max_en to minimum of both passes
        if (saddle_en_1 < max_en) { max_en = saddle_en_1; }
        if (saddle_en_2 < max_en) { max_en = saddle_en_2; }
    }

    free(pt1);
    free(pt2);
    free(pt1_bwd);
    free(pt2_bwd);

    return all_paths;
}

// private functions below

auto ptable_from_string(short* pt, short* loop, const char* s, short* stack)
{
    /*
    this function combines: (utils/structure_utils.c)
        vrna_loopidx_from_ptable(const short *pt)
        vrna_ptable_from_string(structure, VRNA_BRACKETS_***);

    the stack will be overwritten with every call, it needs to be allocated beforehand.

    input:  structure as string (char* s)
    output: pairing table: short* pt
            loop table:    short* loop
    */

    const char*  ptr;
    unsigned int i, j, n;
    int          hx = 0;
    int          l  = 0;
    int          nl = 0;

    n = (unsigned int)pt[0];

    // reset everything to unpaired
    memset(pt + 1, 0, sizeof(short) * pt[0]);
    memset(loop + 1, 0, sizeof(short) * pt[0]);
    // whatever currently in the stack is does not matter...

    const char open  = '(';
    const char close = ')';

    // fmt::print("before:{}\n", vrna_db_from_ptable(pt));
    for (hx = 0, i = 1, ptr = s; (i <= n) && (*ptr != '\0'); ptr++, i++) {
        if (*ptr == open) {
            nl++;
            l           = nl;
            stack[hx++] = i;
        }

        loop[i] = l;

        if (*ptr == close) {
            j = stack[--hx];

            if (hx > 0) {
                l = loop[stack[hx - 1]]; /* index of enclosing loop   */
            } else {
                l = 0; /* external loop has index 0 */
            }
            // no time for exception handling...
            // if (hx < 0) {
            //     vrna_message_warning(
            //         "%s\nunbalanced brackets found while extracting base pairs", s);
            // }
            pt[i] = j;
            pt[j] = i;
        }
    }
    pt[0]   = n;
    loop[0] = nl;
}

inline auto single_findpath::try_moves(vrna_fold_compound_t* vc, intermediate_t c, int maxE,
                                       intermediate_t* next, int dist, int bp_dist, short* temp_pt,
                                       short* stack) -> int
{
    return 1;
}

// size_t hash_c_string(const char* p)
// {
//     size_t       result = 0;
//     const size_t prime  = 31;
//     for (size_t i = 0; i < 599; ++i) { result = p[i] + (result * prime); }

//     // much slower
//     // for(; *p!= '\0'; ++p) {
//     //     result = *p + (result * prime);
//     // }

//     return result;
// }

uint32_t hash_string(const char* s, size_t size)
{
    uint32_t hash = 0;

    for (size_t i = 0; i < size; ++i) {
        hash += s[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }

    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);

    return hash;
}

// removed global direction, path

inline auto single_findpath::find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2,
                                            const int current_search_width, int maxE,
                                            bool direction, int final_search_width)
    -> std::vector<sorted_path>
{
    move_t* mlist;

    int i, d;

    int bp_dist = 0;
    int len     = static_cast<int>(pt1[0]);

    short* temp_pt = vrna_ptable_copy(pt1);
    // temp_pt[0] = pt1[0];

    short* stack     = vrna_ptable_copy(pt1);
    short* temp_loop = vrna_ptable_copy(pt1);

    // mlist = (move_t*)vrna_alloc(sizeof(move_t) * len); /* bp_dist < n */

    // fmt::print("init: maxE: {} sw: {} dir: {}\n", maxE, current_search_width, direction);

    // generate all possible moves
    // this also calculates bp_dist...

    // generate move list. uses push backs, we don't know bp_dist (vector size) yet
    // std::vector<move_t> move_list{};
    // std::vector<short> move_list{};

    std::vector<move_ij> move_list{};

    std::vector<size_t> h_list{};

    // move_list.push_back({0,0,0});

    for (i = 1; i <= len; i++) {
        if (pt1[i] != pt2[i]) {
            if (i < pt1[i]) {
                /* need to delete this pair */
                // move_list.push_back(0);
                move_list.push_back({static_cast<short>(-i), static_cast<short>(-pt1[i])});

                size_t move_hash;
                move_hash = int_hash_64(-i);
                move_hash *= (-pt1[i]);
                move_hash += (-pt1[i]) * 13;
                h_list.push_back(move_hash);

                bp_dist++;
            }

            if (i < pt2[i]) {
                /* need to insert this pair */
                // move_list.push_back(0);
                move_list.push_back({static_cast<short>(i), static_cast<short>(pt2[i])});

                size_t move_hash;
                move_hash = int_hash_64(i);
                move_hash *= (pt2[i]);
                move_hash += (pt2[i]) * 13;
                h_list.push_back(move_hash);

                bp_dist++;
            }
        }
    }
    // move_list.push_back(0);
    move_list.push_back({0, 0});

    // memory pool for structures and moves
    // we're using a double buffer such that current results don't get immediately overwritten in
    // the next iteration (the pointers from current vector need to be valid for the next vector)

    std::vector<std::vector<char>> structure_storage(2 * current_search_width + 4);
    for (auto& col : structure_storage) {
        col.reserve(len + 1);
    }  // column space allocation so that we can memcpy into it

    std::vector<std::vector<short>> move_storage(2 * current_search_width + 4);
    for (auto& col : move_storage) { col.reserve(bp_dist + 1); }

    std::vector<intermediate_t> current_vec(current_search_width + 1);
    intermediate_t*             current = current_vec.data();

    std::vector<intermediate_t> next_vec(bp_dist * current_search_width + 1);
    intermediate_t*             next = next_vec.data();

    // Initialization data point

    char* init_s = vrna_db_from_ptable(pt1);

    current[0].last_id = 0;
    current[0].s_hash  = 0;

    current[0].s         = init_s;
    current[0].length    = pt1[0];
    current[0].saddle_en = current[0].curr_en = vrna_eval_structure_pt(vc, pt1);

    // every intermediate saves the move index, which references i and j
    std::vector<short> empty_move_index(bp_dist + 1, 0);
    current[0].moves = empty_move_index.data();

    int c = 0;

    std::vector<sorted_path> all_paths;

    for (d = 1; d <= bp_dist; d++) {
        /* go through the distance classes */
        int             u, num_next = 0;
        intermediate_t* cc;

        // memory pool double buffer offset
        int pool_offset = 0;
        if (d % 2 == 0) {
            // pool_offset = bp_dist * current_search_width + 2;
            pool_offset = current_search_width + 2;
        }

        // fmt::print("d: {} \n", d);

        // find valid moves (extracted from original try moves function)
        for (c = 0; current[c].s != nullptr; c++) {
            int en;
            // move_t* mv;

            // update temporary pairing table & loop table to match the current structure string
            ptable_from_string(temp_pt, temp_loop, current[c].s, stack);
            // int* loopidx = vrna_loopidx_from_ptable(temp_pt);

            // fmt::print("d: {}, c: {}, s: {}\n", d, c, current[c].s);

            for (int a = 0; a < move_list.size(); a++) {
                // fmt::print ("d: {}, a: {}, num_next: {}, c: {}\n", d, a, num_next, c);

                const int move_index = current[c].moves[a];
                if (move_index > 0) {
                    continue;
                    // the move at index a was already taken
                }

                const int i = move_list[a].i;
                const int j = move_list[a].j;

                int source_1;
                int source_2;
                int dest_1;
                int dest_2;

                if (j < 0) {  // bp deletion
                    source_1 = -i;
                    dest_1   = 0;
                    source_2 = -j;
                    dest_2   = 0;

                } else {                                         // bp insertion
                    if ((temp_loop[i] == temp_loop[j]) and       // i and j belong to same loop
                        (temp_pt[i] == 0) and (temp_pt[j] == 0)  // and are unpaired
                    ) {
                        source_1 = i;
                        dest_1   = j;
                        source_2 = j;
                        dest_2   = i;

                    } else {
                        // free(pt);
                        continue; /* llegal move, try next; */
                    }
                }

                // this does not work with LOOP_EN
                en = current[c].curr_en + vrna_eval_move_pt(vc, temp_pt, i, j);
                // en = vrna_eval_structure_pt(vc, pt);

                // this used to be en < maxE
                if (en <= maxE) {
                    // take the current string, depending on if it is a delete or add move,
                    // adjust the string accordingly, calculate hash value, and return it as it was
                    // before

                    if (dest_1 == 0) {  // del
                        current[c].s[source_1 - 1] = '.';
                        current[c].s[source_2 - 1] = '.';
                    } else {  // add
                        current[c].s[source_1 - 1] = '(';
                        current[c].s[source_2 - 1] = ')';
                    }

                    // precalculated hash values
                    next[num_next].s_hash = current[c].s_hash + h_list[a];

                    // with hash calculation for every string
                    // next[num_next].s_hash      = hash_cstr(current[c].s, len-1);

                    next[num_next].move_delete = dest_1;
                    next[num_next].move_i      = source_1;
                    next[num_next].move_j      = source_2;

                    if (dest_1 == 0) {  // undo del
                        current[c].s[source_1 - 1] = '(';
                        current[c].s[source_2 - 1] = ')';
                    } else {  // undo add
                        current[c].s[source_1 - 1] = '.';
                        current[c].s[source_2 - 1] = '.';
                    }

                    next[num_next].length = current[c].length;
                    next[num_next].saddle_en =
                        (en > current[c].saddle_en) ? en : current[c].saddle_en;
                    next[num_next].curr_en = en;

                    // set last_id (from which node we came from)
                    next[num_next].last_id = c;
                    // set move_id to current bp dist (further down after sorting)
                    next[num_next].move_id = a;

                    // fmt::print ("d: {}, c: {}, a: {}\n", d, c, a);
                    num_next++;
                }
            }
        }

        if (num_next == 0) {
            // case where we don't find any moves -> abort
            for (cc = current; cc->s != nullptr; cc++) free_intermediate(cc);
            current[0].saddle_en = INT_MAX;
            break;
        }

        /* remove duplicates via sort|uniq
         * if this becomes a bottleneck we can use a hash instead */
        std::qsort(next, num_next, sizeof(intermediate_t), compare_ptable);

        // std::stable_sort(next, next + num_next, [](const auto& a, const auto& b) -> bool {
        //     if (a.s_hash < b.s_hash) { return true; }
        //     if (a.s_hash > b.s_hash) { return false; }
        //     if (a.saddle_en < b.saddle_en) { return true; }
        //     if (a.saddle_en > b.saddle_en) { return false; }
        //     if (a.curr_en < b.curr_en) { return true; }
        //     return false;
        // });

        // std::stable_sort(next, next + num_next, [](const auto& a, const auto& b) -> bool {
        //     if (a.curr_en < b.curr_en) { return true; }
        //     if (a.curr_en > b.curr_en) { return false; }
        //     if (a.s_hash < b.s_hash) { return true; }
        //     if (a.s_hash > b.s_hash) { return false; }
        //     if (a.saddle_en < b.saddle_en) { return true; }
        //     return false;
        // });

        bool flag = true;

        // dont delete duplicates at the end
        if (d == bp_dist and current_search_width >= final_search_width) {
            flag = false;
            for (u = 0, c = 1; c < num_next; c++) { next[++u] = next[c]; }
            num_next = u + 1;
        }

        // this deletes duplicates and shrinks the next array which makes the following sorting step
        // faster
        if (d <= bp_dist and flag) {
            for (u = 0, c = 1; c < num_next; c++) {
                // fmt::print ("d: {}, c: {}, hash: {}, en1: {}, en2: {}\n", d, c, next[c].s_hash,
                // next[c].saddle_en, next[c].curr_en);

                if (next[u].s_hash != next[c].s_hash) {
                    next[++u] = next[c];

                } else {
                    // fmt::print("F\n");
                    free_intermediate(next + c);
                }
            }
            num_next = u + 1;
        }

        // py::print("C:", u, "/", c, sizeof(intermediate_t));
        // fmt::print ("C: {} / {}\n", u, c);

        // if (d <= bp_dist and flag) {
        //     auto last_hash  = next[1].s_hash;
        //     int  best_en    = next[1].saddle_en;
        //     int  best_index = 1;

        //     for (u = 0, c = 1; c < num_next; c++) {
        //         if (next[c].s_hash != last_hash) {
        //             // we now now the best structure
        //             next[u] = next[best_index];
        //             u++;

        //             // reset for next structure
        //             last_hash  = next[c].s_hash;
        //             best_en    = next[c].saddle_en;
        //             best_index = c;

        //         } else {
        //             // we have 2+ consecutive structures with potentially different energies
        //             if (next[c].saddle_en < best_en) {
        //                 best_en    = next[c].saddle_en;
        //                 best_index = c;
        //             }
        //         }
        //     }
        // }

        num_next = u + 1;
        // std::qsort(next, num_next, sizeof(intermediate_t), compare_energy);

        // std::stable_sort(next, next + num_next, [](const auto& a, const auto& b) -> bool {
        //     if (a.saddle_en < b.saddle_en) { return true; }
        //     if (a.saddle_en > b.saddle_en) { return false; }
        //     if (a.curr_en < b.curr_en) { return true; }
        //     return false;
        // });

        std::nth_element(next, next + current_search_width + 1, next + num_next,
                         [](const auto& a, const auto& b) -> bool {
                             // return a.saddle_en < b.saddle_en;
                             if (a.saddle_en < b.saddle_en) { return true; }
                             if (a.saddle_en > b.saddle_en) { return false; }
                             if (a.curr_en < b.curr_en) { return true; }
                             return false;
                         });
                         
        auto max_sort = std::min(current_search_width + 1, num_next);
        std::stable_sort(next, next + max_sort, [](const auto& a, const auto& b) -> bool {
            if (a.saddle_en < b.saddle_en) { return true; }
            if (a.saddle_en > b.saddle_en) { return false; }
            if (a.curr_en < b.curr_en) { return true; }
            return false;
        });

        // auto max_sort = std::min(current_search_width+1, num_next);
        // std::partial_sort(next, next+max_sort, next + num_next, [](const auto& a, const auto& b)
        // -> bool {
        //     // return a.saddle_en < b.saddle_en;
        //     if (a.saddle_en < b.saddle_en) { return true; }
        //     if (a.saddle_en > b.saddle_en) { return false; }
        //     if (a.curr_en < b.curr_en) { return true; }
        //     return false;
        // });

        // update moves & structures for next iteration
        for (u = 0; u < current_search_width and u < num_next; u++) {
            // fmt::print ("d: {}/{}, u: {}, hash: {}, en1: {}, en2: {}\n", d, bp_dist, u,
            // next[u].s_hash, next[u].saddle_en, next[u].curr_en);

            auto last_id = next[u].last_id;
            auto move_id = next[u].move_id;

            // fmt::print("u: {}, last_id: {}, move_id: {} \n", u, last_id, move_id);

            // fmt::print("u: {}, last_id: {}, move_id: {} / {} {}\n", u, last_id, move_id,
            //            current[last_id].moves == nullptr, current[last_id].saddle_en);

            // whatever currently in next[u].moves is needs to be discarded -
            // allocate memory, overwrite with moves from last node, then update

            next[u].moves =
                move_storage[u + pool_offset].data();  // set move pointer (valid for 2 iterations)
            memcpy(next[u].moves, current[last_id].moves, sizeof(short) * (bp_dist + 1));
            next[u].moves[move_id] = d;  // update current move id to current basepair distance

            // update structure

            auto dest_1   = next[u].move_delete;
            auto source_1 = next[u].move_i;
            auto source_2 = next[u].move_j;

            memcpy(structure_storage[u + pool_offset].data(), current[last_id].s,
                   (current[last_id].length) * sizeof(char));
            next[u].s = structure_storage[u + pool_offset].data();  // set structure pointer

            if (dest_1 == 0) {
                next[u].s[source_1 - 1] = '.';
                next[u].s[source_2 - 1] = '.';
            } else {
                next[u].s[source_1 - 1] = '(';
                next[u].s[source_2 - 1] = ')';
            }
            // fmt::print("done2\n");
            // fmt::print("u: {} {} {}\n", next[u].moves[0], next[u].moves[1], next[u].moves[2]);
        }

        /* free the old stuff */
        for (cc = current; cc->s != nullptr; cc++) { free_intermediate(cc); }

        // adjust current for the next iteration
        for (u = 0; u < current_search_width and u < num_next; u++) {
            current[u] = next[u];

            // test
        }

        for (; u < num_next; u++) { free_intermediate(next + u); }
        num_next = 0;
    }

    // dummy path
    if (bp_dist == 0) { c = 1; }

    for (int index = 0; index < c; index++) {
        if (not current[index].moves) { continue; }

        // fmt::print("E2 {} {} {}\n", index, c, current[index].moves[0]);

        // allocate bp_dist elements in current_path vector and set max_en
        auto& current_path = all_paths.emplace_back(bp_dist, current[index].saddle_en);

        // std::qsort(current[0].moves, bp_dist, sizeof(move_t), compare_moves_when);
        for (d = 0; d < bp_dist; d++) {
            int current_dist = current[index].moves[d] - 1;

            int i, j;
            i = move_list[d].i;
            j = move_list[d].j;

            if (direction) {
                // fwd path
                current_path.moves[current_dist] = {i, j, 0};
            } else {
                // bwd path, fill vector back to front
                current_path.moves[bp_dist - current_dist - 1] = {-i, -j, 0};
            }
        }
    }

    free(init_s);
    free(temp_pt);
    free(temp_loop);
    free(stack);

    return all_paths;
}

inline auto single_findpath::free_intermediate(intermediate_t* i) -> void
{
    i->s         = nullptr;
    i->moves     = nullptr;
    i->saddle_en = INT_MAX;
}

inline auto single_findpath::compare_ptable(const void* A, const void* B) -> int
{
    const intermediate_t* a = (const intermediate_t*)A;
    const intermediate_t* b = (const intermediate_t*)B;

    /*
    Identical structures are grouped by identical energy & structure hash.
    We can sort first by s_hash or curr_en: since comparing 32 bit integers
    is cheaper than 64 bit ints, curr_en is compared first.

    (ascending or descending ordering is irrelevant here)
    */

    if (a->curr_en > b->curr_en) { return -1; }
    if (a->curr_en < b->curr_en) { return 1; }
    if (a->s_hash > b->s_hash) { return -1; }
    if (a->s_hash < b->s_hash) { return 1; }

    // if structures are identical sort by saddle energy: lowest energy to [0]
    if (a->saddle_en < b->saddle_en) { return -1; }
    if (a->saddle_en > b->saddle_en) { return 1; }
    return 0;
}

// inline auto single_findpath::compare_ptable(const void* A, const void* B) -> int
// {
//     intermediate_t *a, *b;
//     int             c;

//     a = (intermediate_t*)A;
//     b = (intermediate_t*)B;

//     // if both structures are not identical, sort them according to the hash value
//     // if (a->s_hash != b->s_hash) {
//         if (a->s_hash > b->s_hash) {
//             return 1;
//         } else {
//             return -1;
//         }

//     // }

//     // return 1;
//     // if ((a->saddle_en - b->saddle_en) != 0) return a->saddle_en - b->saddle_en;
//     // return a->curr_en - b->curr_en;

// }

inline auto single_findpath::compare_energy(const void* A, const void* B) -> int
{
    const intermediate_t* a = (const intermediate_t*)A;
    const intermediate_t* b = (const intermediate_t*)B;

    if ((a->saddle_en - b->saddle_en) != 0) return a->saddle_en - b->saddle_en;

    return a->curr_en - b->curr_en;
}
