// this includes the basic elementary findpath algorithm

// // basic findpath structs
// struct move_ij {
//     short i; /* i,j>0 insert; i,j<0 delete */
//     short j;
//     short when; /* 0 if still available, else resulting distance from start */
//     // int E;
// };

// struct intermediate_multi {
//     // short*  pt; /**<  @brief  pair table */
//     char*   s;
//     short   length;
//     int     saddle_en; /**<  @brief  saddle energy so far */
//     int     curr_en;   /**<  @brief  current energy */
//     move_ij* moves;     /**<  @brief  remaining moves to target */
// };

// idea: https://stackoverflow.com/questions/35826416/why-no-default-hash-for-c-pod-structs

// // basic findpath structs
// struct move_s {
//     short i; /* i,j>0 insert; i,j<0 delete */
//     short j;
//     short when; /* 0 if still available, else resulting distance from start */
//                 // int E;

//     bool const operator==(const move_s& o) const { return i == o.i && j == o.j; }
//     bool const operator<(const move_s& o) const { return i < o.i || (i == o.i && j < o.j); }
// };
// // The specialized hash function for `unordered_map` keys
// struct hash_fn {
//     std::size_t operator()(const move_s& node) const
//     {
//         std::size_t h1 = std::hash<short>()(node.i);
//         std::size_t h2 = std::hash<short>()(node.j);

//         return h1 ^ h2;
//     }
// };

struct intermediate_multi {
    // short*  pt; /**<  @brief  pair table */
    char*  s;
    short* moves;
    short* available_targets;
    short  length;
    int    last_move_index;
    int    moves_available;

    int saddle_en; /**<  @brief  saddle energy so far */
    int curr_en;   /**<  @brief  current energy */
    // move_t* moves;     /**<  @brief  remaining moves to target */
};

class multi_findpath
{
   private:
    // vrna_fold_compound_t* fc;
    // short*                pt1;
    // short*                pt2;
    // int                   final_search_width;

    // static auto compare_ptable(const void* A, const void* B) -> int;
    // static auto compare_energy(const void* A, const void* B) -> int;
    // static auto compare_moves_when(const void* A, const void* B) -> int;
    // static auto copy_moves(move_ij* mvs, int bp_dist) -> move_ij*;

    std::vector<std::vector<move_ij>> move_lists{};
    int                               len_move_lists;

    std::vector<move_ij> combined_moves{};
    std::vector<int>     move_list_index{};

    std::vector<int> s_code{};
    std::set<int>    s_code_set;

    int bp_dist_min = 999999;
    int bp_dist_max = 0;

    std::vector<std::vector<bool>> allowed_moves{};
    std::vector<std::vector<bool>> move_matrix{};

    auto find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2,
                        const int current_search_width, int maxE, bool direction,
                        int final_search_width) -> std::vector<sorted_path>;

    static auto free_intermediate(intermediate_multi* i) -> void;
    static auto compare_ptable(const void* A, const void* B) -> int;
    static auto compare_energy(const void* A, const void* B) -> int;
    // static auto compare_moves_when(const void* A, const void* B) -> int;
    // static auto copy_moves(move_t* mvs, int bp_dist) -> move_t*;

   public:
    // auto init(vrna_fold_compound_t* fc, short* pt1, short* pt2, int final_search_width,
    //           bool mp = true) -> std::vector<sorted_path>;

    multi_findpath(vrna_fold_compound_t* fc, std::string s, std::vector<std::string> destinations,
                   float search_width_multiplier, int en_limit);
};

// struct move_ij {
//     int x, y;
// };

multi_findpath::multi_findpath(vrna_fold_compound_t* fc, std::string s,
                               std::vector<std::string> destinations, float search_width_multiplier,
                               int en_limit)
{
    auto start = std::chrono::system_clock::now();

    // fmt::print("reference: {}\n", s);

    int    len = s.size();
    short* pt  = vrna_ptable(s.c_str());

    auto equal = [](const move_ij& c1, const move_ij& c2) {
        if (c1.i == c2.i and c1.j == c2.j) { return true; }
        return false;
    };
    auto hash = [](const move_ij& n) {
        size_t res = 17;
        res        = res * 31 + std::hash<int>()(n.i);
        res        = res * 31 + std::hash<int>()(n.j);
        return res;
    };

    // std::vector<int>                  move_list_index(destinations.size(), 0);

    // std::unordered_map<move_ij, int, decltype(hash), decltype(equal)> combined_moves;
    int len_moves = 0;
    move_list_index.resize(destinations.size());

    s_code.resize(destinations.size());

    for (auto& current_s : destinations) {
        static int d       = 0;
        int        bp_dist = 0;
        // fmt::print("s: {}\n", current_s);
        auto&  current_move_list = move_lists.emplace_back();
        short* current_pt        = vrna_ptable(current_s.c_str());

        for (int i = 1; i <= len; i++) {
            if (pt[i] != current_pt[i]) {
                if (i < pt[i]) {
                    /* need to delete this pair */
                    current_move_list.push_back(
                        {static_cast<short>(-i), static_cast<short>(-pt[i])});
                    combined_moves.push_back({static_cast<short>(-i), static_cast<short>(-pt[i])});
                    s_code[d] -= i * 7;
                    s_code[d] -= pt[i] * 3;
                    bp_dist++;
                }

                if (i < current_pt[i]) {
                    /* need to insert this pair */
                    current_move_list.push_back(
                        {static_cast<short>(i), static_cast<short>(current_pt[i])});
                    combined_moves.push_back(
                        {static_cast<short>(i), static_cast<short>(current_pt[i])});
                    s_code[d] += i * 7;
                    s_code[d] += current_pt[i] * 3;
                    bp_dist++;
                }
            }
        }
        // sort current path
        std::sort(current_move_list.begin(), current_move_list.end(),
                  [](const auto& a, const auto& b) -> bool {
                      if (a.i != b.i) { return a.i < b.i; }
                      return a.j < b.j;
                  });

        if (current_move_list.size() < bp_dist_min) { bp_dist_min = current_move_list.size(); }

        // fmt::print("max: {} {} {}\n",bp_dist_max, current_move_list.size(),
        // current_move_list.size() > bp_dist_max);
        if (current_move_list.size() > bp_dist_max) { bp_dist_max = current_move_list.size(); }

        s_code_set.insert(s_code[d]);
        // fmt::print("d : {} / code: {} / ", d, s_code[d]);
        // for (auto const& m : current_move_list) { fmt::print("({}/{}) ", m.i, m.j); }
        // fmt::print("\n");
        free(current_pt);
        d++;
    }

    std::sort(combined_moves.begin(), combined_moves.end(),
              [](const auto& a, const auto& b) -> bool {
                  if (a.i != b.i) { return a.i < b.i; }
                  return a.j < b.j;
              });

    // sort unique
    combined_moves.erase(std::unique(combined_moves.begin(), combined_moves.end(), equal),
                         combined_moves.end());

    // print all available moves
    // for (auto const& m : combined_moves) { fmt::print("({}/{}) ", m.i, m.j); }
    // fmt::print("\n");

    // build a 2d array with (size of combined moves)² to keep track which moves allowed where

    // std::vector<std::vector<bool>> move_matrix(move_lists.size());

    len_move_lists = move_lists.size();

    move_matrix.resize(move_lists.size());
    for (auto& row : move_matrix) { row.resize(combined_moves.size()); }

    // std::vector<std::vector<bool>> allowed_moves(combined_moves.size());

    allowed_moves.resize(combined_moves.size());
    for (auto& row : allowed_moves) { row.resize(combined_moves.size()); }

    // check which moves every destination requires
    // fill up a m x n matrix first
    for (int i = 0; i < combined_moves.size(); i++) {
        const auto& outer_move = combined_moves[i];

        // fmt::print("current move ({}/{}): \n", outer_move.i, outer_move.j);
        for (int j = 0; j < move_lists.size(); j++) {
            int&  index      = move_list_index[j];
            auto& inner_move = move_lists[j][index];

            // fmt::print("{}: ({}/{}) ({}) ", j, inner_move.i, inner_move.j, index);
            if (equal(inner_move, outer_move)) {
                move_matrix[j][i] = true;
                index++;  // move_list_index[j]++
                // inner_move = move_lists[j][index];
                // fmt::print("available\n");
                // fmt::print("{}: ({}/{}) ({}) \n", j, inner_move.i, inner_move.j, index);
            }
            // else: skip move, [j][i] remains false
        }
    }

    // print m x n matrix
    // fmt::print("A\n");
    // for (auto& row : move_matrix) {
    //     for (auto m : row) {
    //         if (m)
    //             fmt::print("1, ");
    //         else
    //             fmt::print("0, ");
    //     }
    //     fmt::print("\n");
    // }
    // fmt::print("\n");

    // multiply A^T with A (boolean matrix multiplication)
    int n = combined_moves.size();
    // https://stackoverflow.com/questions/29588384/boolean-matrix-multiplication-algorithm
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            bool value = false;
            for (int m = 0; m < move_lists.size(); m++) {
                // value |= move_matrix[i][m] && move_matrix[j][m];
                value |= move_matrix[m][i] && move_matrix[m][j];

                if (value) break;  // early out
            }
            allowed_moves[i][j] = value;
        }
    }

    // print m x m matrix
    // fmt::print("A^TA\n");
    // for (auto& row : allowed_moves) {
    //     for (auto m : row) {
    //         if (m)
    //             fmt::print("1, ");
    //         else
    //             fmt::print("0, ");
    //     }
    //     fmt::print("\n");
    // }
    // fmt::print("\n");

    // allowed & disallowed vector indices (int last move; check last move, set when to -1)

    short* pt1 = vrna_ptable(s.c_str());
    short* pt2 = vrna_ptable(destinations[1].c_str());

    // int  current_search_width = 1500;
    // int  final_search_width   = 1500;
    // int  current_search_width = 150;
    // int  final_search_width   = 150;
    // int  current_search_width = 40;
    // int  final_search_width   = 40;
    int current_search_width = static_cast<int>(bp_dist_max * 1.0 * search_width_multiplier);
    int final_search_width   = static_cast<int>(bp_dist_max * 1.0 * search_width_multiplier);

    // int maxE = -1260;
    // int  maxE                 = -1600;

    bool direction = true;

    fmt::print("bp dist min: {}\n", bp_dist_min);
    fmt::print("bp dist max: {}\n", bp_dist_max);
    fmt::print("search width: {}\n", current_search_width);

    auto result =
        find_path_once(fc, pt1, pt2, current_search_width, en_limit, direction, final_search_width);

    std::sort(result.begin(), result.end(),
              [](const auto& a, const auto& b) -> bool { return a.max_en < b.max_en; });

    // int bp_dist = 2;
    // s_graph G_inner{fc, pt1, pt2, bp_dist, result};
    // G_inner.display_path();

    fmt::print("~~~~~~~~~~\n");

    for (const auto& path : result) { fmt::print("{} -> {}\n", path.destination, path.max_en); }

    auto end     = std::chrono::system_clock::now();
    auto elapsed = end - start;
    fmt::print("runtime int: {}\n", elapsed.count() / 1000000000.0);
    return;
}

inline auto multi_findpath::find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2,
                                           const int current_search_width, int maxE, bool direction,
                                           int final_search_width) -> std::vector<sorted_path>
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

    // move_list.push_back({0,0,0});

    // for (i = 1; i <= len; i++) {
    //     if (pt1[i] != pt2[i]) {
    //         if (i < pt1[i]) {
    //             /* need to delete this pair */
    //             // move_list.push_back(0);
    //             move_list.push_back({static_cast<short>(-i), static_cast<short>(-pt1[i])});
    //             bp_dist++;
    //         }

    //         if (i < pt2[i]) {
    //             /* need to insert this pair */
    //             // move_list.push_back(0);
    //             move_list.push_back({static_cast<short>(i), static_cast<short>(pt2[i])});
    //             bp_dist++;
    //         }
    //     }
    // }
    // move_list.push_back({0, 0});

    move_list = combined_moves;
    // move_list.push_back({0, 0});

    bp_dist = combined_moves.size();  // individual paths might have a lower bp_dist

    // short* stack = new short[200];
    // fmt::print("input: {}\n", vrna_db_from_ptable(pt2));
    // fmt::print("input: {}\n", vrna_db_from_ptable(pt1));
    // fmt::print("input: {}\n", vrna_db_from_ptable(temp_pt));
    // fmt::print("\n");

    // memory pool for structures and moves
    // we're using a double buffer such that current results don't get immediately overwritten in
    // the next iteration (the pointers from current vector need to be valid for the next vector)

    std::vector<std::vector<char>> s_pool(2 * bp_dist * current_search_width + 1);
    for (auto& col : s_pool) {
        col.reserve(len);  // len+1
    }                      // column space allocation so that we can memcpy into it

    std::vector<std::vector<short>> move_pool(2 * bp_dist * current_search_width + 1);
    for (auto& col : move_pool) { col.reserve(bp_dist + 1); }

    std::vector<std::vector<short>> available_targets(2 * bp_dist * current_search_width + 1);
    for (auto& col : available_targets) { col.reserve(len_move_lists); }

    std::vector<intermediate_multi> current_vec((current_search_width + 1) * len_move_lists);
    intermediate_multi*             current = current_vec.data();

    std::vector<intermediate_multi> next_vec(bp_dist * current_search_width + 1);
    intermediate_multi*             next = next_vec.data();

    // Initialization data point

    char* init_s = vrna_db_from_ptable(pt1);

    current[0].s         = init_s;
    current[0].length    = pt1[0];
    current[0].saddle_en = current[0].curr_en = vrna_eval_structure_pt(vc, pt1);

    // every intermediate saves the move index, which references i and j
    std::vector<short> empty_move_index(bp_dist + 1, 0);
    current[0].moves = empty_move_index.data();

    // all moves available at start
    std::vector<short> full_available_targets(len_move_lists, 1);
    current[0].available_targets = full_available_targets.data();

    current[0].last_move_index = -1;
    current[0].moves_available = combined_moves.size();

    // ptable_from_string(temp_pt, temp_loop, current[0].s, stack);

    // next =
    //     (intermediate_multi*)vrna_alloc(sizeof(intermediate_multi) * (bp_dist *
    //     current_search_width + 1));

    std::vector<sorted_path> all_paths;

    for (d = 1; d <= bp_dist; d++) {
        /* go through the distance classes */
        int                 c, u, num_next = 0;
        intermediate_multi* cc;

        // memory pool double buffer offset
        int pool_offset = 0;
        if (d % 2 == 0) { pool_offset = bp_dist * current_search_width + 2; }

        // fmt::print("d: {} \n", d);

        // find valid moves (extracted from original try moves function)
        for (c = 0; current[c].s != nullptr; c++) {
            int en;

            // update temporary pairing table & loop table to match the current structure string
            ptable_from_string(temp_pt, temp_loop, current[c].s, stack);

            int all_move_indices = 0;
            // extra pass to filter out all move branches which are invalid
            for (int a = 0; a < move_list.size(); a++) {
                const int move_index = current[c].moves[a];
                if (move_index != 0) {
                    if (move_index > 0) {
                        // if we make an error, we're just adding an additional sequence to the
                        // result vector
                        all_move_indices += move_list[a].i * 7;
                        all_move_indices += move_list[a].j * 3;
                    }
                    continue;
                    // the move at index a was already taken
                }
                const int last_move_index = current[c].last_move_index;

                if (last_move_index >= 0 and not allowed_moves[last_move_index][a]) {
                    const int allowed_move = allowed_moves[last_move_index][a];
                    // fmt::print("not allowed: {} / last_move: {} / c: {}\n", a, last_move_index,
                    // c);
                    current[c].moves[a] = -1;
                    current[c].moves_available--;
                }

                // if (current[c].moves[a] > 0) {
                // all_move_indices += current[c].moves[a];  // del
                // }
            }

            // if (c >= 400) {
            //     fmt::print("c: {} / 1:{} / 2:{} / 3:{} \n", c, current[c].available_targets[0],
            //                current[c].available_targets[1], current[c].available_targets[2]);
            // }

            // if (d >= bp_dist_min) {
            //     fmt::print("d: {} / {} / max_en: {}\n", d, current[c].s, current[c].saddle_en);
            //     // test
            // }

            // fmt::print("d: {}, c:{}, moves av: {} / code: {}:{}\n", d, c,
            //            current[c].moves_available, all_move_indices,
            //            s_code_set.contains(all_move_indices));

            // if (current[c].moves_available == 0) {
            if (current[c].moves_available == 0 or (d >= bp_dist_min and all_move_indices != 0 and
                                                    s_code_set.contains(all_move_indices))) {
                fmt::print("{} return: d: {}, c:{}, moves av: {} / code: {}\n", current[c].s, d, c,
                           current[c].moves_available, all_move_indices);

                std::string destination(current[c].s);

                auto& current_path =
                    all_paths.emplace_back(bp_dist, current[c].saddle_en, destination);

                // std::qsort(current[0].moves, bp_dist, sizeof(move_t), compare_moves_when);

                for (int current_d = 0; current_d < bp_dist; current_d++) {
                    int i, j;

                    int current_dist = current[c].moves[current_d] - 1;

                    i = move_list[current_d].i;
                    j = move_list[current_d].j;

                    if (current_dist < 0) { continue; }

                    // fmt::print("d: {}, cd: {}, i:{}, j:{}\n", current_d, current_dist, i, j);

                    if (direction) {
                        // fwd path
                        current_path.moves[current_dist] = {i, j, 0};
                    } else {
                        // bwd path, fill vector back to front
                        current_path.moves[bp_dist - current_dist - 1] = {-i, -j, 0};
                    }
                }

                continue;
            }

            for (int a = 0; a < move_list.size(); a++) {
                const int move_index = current[c].moves[a];
                if (move_index != 0) {
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
                    // fmt::print(
                    //     "d: {} / move_index: {} / i: {} / j: {} / en: {} / av. moves: {} / acc: "
                    //     "{}\n",
                    //     d, a, i, j, en, current[c].moves_available, all_move_indices);

                    // reset everything to unpaired
                    // memset(s_pool[num_next+pool_offset].data(), 0, sizeof(char) * pt[0]);

                    memcpy(s_pool[num_next + pool_offset].data(), current[c].s,
                           (current[c].length) * sizeof(char));
                    // memcpy(s_pool[num_next + pool_offset], current[c].s,
                    //        (current[c].length) * sizeof(char));

                    if (dest_1 == 0) {
                        s_pool[num_next + pool_offset][source_1 - 1] = '.';
                        s_pool[num_next + pool_offset][source_2 - 1] = '.';
                    } else {
                        s_pool[num_next + pool_offset][source_1 - 1] = '(';
                        s_pool[num_next + pool_offset][source_2 - 1] = ')';
                    }

                    next[num_next].s =
                        s_pool[num_next + pool_offset].data();  // set structure pointer
                    // next[num_next].s    = s_pool[num_next + pool_offset];

                    next[num_next].length = current[c].length;

                    next[num_next].saddle_en =
                        (en > current[c].saddle_en) ? en : current[c].saddle_en;
                    next[num_next].curr_en = en;

                    next[num_next].moves_available = current[c].moves_available - 1;
                    next[num_next].last_move_index = a;
                    // next[num_next].pt      = pt;
                    // mv->when = d;
                    // mv->E                  = en;

                    // next[num_next].moves = (move_t*)vrna_alloc(sizeof(move_t) * (bp_dist +
                    // 1)); memcpy(next[num_next].moves, current[c].moves, sizeof(move_t) *
                    // (bp_dist + 1));

                    memcpy(move_pool[num_next + pool_offset].data(), current[c].moves,
                           sizeof(short) * (bp_dist + 1));
                    next[num_next].moves =
                        move_pool[num_next + pool_offset].data();  // set move pointer

                    move_pool[num_next + pool_offset][a] = d;

                    memcpy(available_targets[num_next + pool_offset].data(),
                           current[c].available_targets, sizeof(short) * (len_move_lists));
                    next[num_next].available_targets =
                        available_targets[num_next + pool_offset].data();  // set move pointer

                    // next[num_next].moves = copy_moves(current[c].moves, bp_dist);
                    // mv->when = 0;

                    num_next++;
                }
            }

            // free(loopidx);
            // return num_next;
        }

        // fmt::print("B\n");

        if (num_next == 0) {
            // case where we don't find any moves -> abort
            // for (cc = current; cc->s != nullptr; cc++) free_intermediate(cc);
            // current[0].saddle_en = INT_MAX;
            break;
        }

        /* remove duplicates via sort|uniq
         * if this becomes a bottleneck we can use a hash instead */
        std::qsort(next, num_next, sizeof(intermediate_multi), compare_ptable);

        // this shrinks the next array which makes the following sorting step faster
        for (u = 0, c = 1; c < num_next; c++) {
            // if (memcmp(next[u].pt, next[c].pt, sizeof(short) * len) != 0) {
            if (memcmp(next[u].s, next[c].s, sizeof(char) * len) != 0) {
                next[++u] = next[c];

            } else {
                // new part - save multiple result paths

                // if (d >= bp_dist && next[c].saddle_en <= maxE &&
                //     current_search_width >= final_search_width) {
                if (next[c].moves_available == 0 and next[c].saddle_en <= maxE &&
                    current_search_width >= final_search_width) {
                    // fmt::print("return\n");

                    std::string destination = next[c].s;

                    // allocate bp_dist elements in current_path vector and set max_en
                    auto& current_path =
                        all_paths.emplace_back(bp_dist, next[c].saddle_en, destination);
                    // sorted_path current_path(bp_dist, next[c].saddle_en);

                    for (int current_d = 0; current_d < bp_dist; current_d++) {
                        int i, j;

                        int current_dist = next[c].moves[current_d] - 1;

                        i = move_list[current_d].i;
                        j = move_list[current_d].j;

                        if (current_dist < 0) { continue; }

                        // fmt::print("d: {}, cd: {}, i:{}, j:{}\n", current_d, current_dist, i, j);

                        if (direction) {
                            // fwd path
                            current_path.moves[current_dist] = {i, j, 0};
                        } else {
                            // bwd path, fill vector back to front
                            current_path.moves[bp_dist - current_dist - 1] = {-i, -j, 0};
                        }
                    }

                    // current_path.max_en = next[c].saddle_en;
                    // all_paths.push_back(std::move(current_path));
                    // free(temp_moves);
                }

                // fmt::print("F\n");

                free_intermediate(next + c);
                // fmt::print("G\n");
            }
        }

        // IC(num_next);

        num_next = u + 1;
        std::qsort(next, num_next, sizeof(intermediate_multi), compare_energy);

        // next is now again reduced to size current, means replace current with next.
        /* free the old stuff */
        for (cc = current; cc->s != nullptr; cc++) { free_intermediate(cc); }

        // set next to current for next iteration

        std::vector<short> available_target_counter(len_move_lists, 0);

        for (u = 0; u < current_search_width && u < num_next; u++) {
            current[u] = next[u];

            // move invalidation

            // auto last_move = current[u].last_move_index;
            // if (last_move < 0) { continue; }

            // for (int b = 0; b < len_move_lists; b++) {
            //     // if a target is already marked unavailable, we don't need to check it again
            //     if (current[u].available_targets[b] == 0) { continue; }
            //     // fmt::print("last move {} target {}\n", last_move, b);
            //     // fmt::print("size 1 {} size2 {}\n", move_matrix[0].size(),move_matrix.size());

            //     auto invalidate = move_matrix[b][last_move];

            //     if (invalidate == 0) {
            //         current[u].available_targets[b] = 0;
            //         // fmt::print("last move {} invalidates target {}\n", last_move, b);
            //         // fmt::print("{} {} {} {}\n", current[u].available_targets[0],
            //         //    current[u].available_targets[1], current[u].available_targets[2],
            //         //    current[u].available_targets[3]);
            //     } else {
            //         available_target_counter[b]++;
            //     }
            // }
        }

        // find extra moves deep in next...
        // int offset = 0;
        // for (int b = 0; b < len_move_lists; b++) {
        //     fmt::print("d: {} target {}: {} counter out of {}\n", d, b, available_target_counter[b],
        //                num_next);

        //     if (available_target_counter[b] < 180) {
        //         auto start = u;

        //         for (; start < num_next; start++) {
        //             auto last_move = next[start].last_move_index;

        //             // auto invalidate = move_matrix[b][last_move];

        //             if (next[start].available_targets[b] == 1 and move_matrix[b][last_move] == 1) {
        //                 // fmt::print("d: {} target {}: available {}/{}\n", d, b, u, start);
        //                 available_target_counter[b]++;

        //                 // todo: check if this invalidates other things...
        //                 // set start to u;
        //                 current[u] = next[start];
        //                 next[start] = next[num_next-1]; // random data entry
        //                 u++;
        //             }

        //             if (available_target_counter[b] > 180) { break; }
        //         }
        //     }
        // }

        // fmt::print("C\n");
        for (; u < num_next; u++) { free_intermediate(next + u); }
        // fmt::print("D\n");
        num_next = 0;
    }

    // free(next);

    if (current[0].moves) {
        // allocate bp_dist elements in current_path vector and set max_en

        std::string destination = current[0].s;

        auto& current_path = all_paths.emplace_back(bp_dist, current[0].saddle_en, destination);

        // std::qsort(current[0].moves, bp_dist, sizeof(move_t), compare_moves_when);

        for (int current_d = 0; current_d < bp_dist; current_d++) {
            int i, j;

            int current_dist = current[0].moves[current_d] - 1;

            i = move_list[current_d].i;
            j = move_list[current_d].j;

            if (current_dist < 0) { continue; }

            // fmt::print("d: {}, cd: {}, i:{}, j:{}\n", current_d, current_dist, i, j);

            if (direction) {
                // fwd path
                current_path.moves[current_dist] = {i, j, 0};
            } else {
                // bwd path, fill vector back to front
                current_path.moves[bp_dist - current_dist - 1] = {-i, -j, 0};
            }
        }

        // for (d = 0; d < bp_dist; d++) {
        //     int current_dist = current[0].moves[d] - 1;

        //     int i, j;
        //     i = move_list[d].i;
        //     j = move_list[d].j;

        //     // i = 0;
        //     // j = 0;

        //     // fmt::print("d={}, i={}, j={}, current_dist={}\n", d, i, j, current_dist);

        //     if (direction) {
        //         // fwd path
        //         current_path.moves[current_dist] = {i, j, 0};
        //         // current_path.moves[d] = {i, j, 0};
        //     } else {
        //         // bwd path, fill vector back to front
        //         current_path.moves[bp_dist - current_dist - 1] = {-i, -j, 0};
        //         // current_path.moves[bp_dist - d - 1] = {-i, -j, 0};
        //     }
        // }

        // current_path.max_en = current[0].saddle_en;
        // all_paths.push_back(std::move(current_path));
    }

    free(init_s);
    free(temp_pt);
    free(temp_loop);
    free(stack);

    // fmt::print("F\n");

    return all_paths;
}

inline auto multi_findpath::free_intermediate(intermediate_multi* i) -> void
{
    i->s         = nullptr;
    i->moves     = nullptr;
    i->saddle_en = INT_MAX;
}

inline auto multi_findpath::compare_ptable(const void* A, const void* B) -> int
{
    intermediate_multi *a, *b;
    int                 c;

    a = (intermediate_multi*)A;
    b = (intermediate_multi*)B;

    c = memcmp(a->s, b->s, a->length * sizeof(char));
    // c = strcmp (a->s, b->s);
    // c = memcmp(a->pt, b->pt, a->pt[0] * sizeof(short));

    if (c != 0) {
        // int d = memcmp(a->s, b->s, 299 * sizeof(char));

        // if (d != c) { fmt::print("error {} {}\n", c, d); }

        // I have no idea why a minus needs to be there...🤦‍♂️
        return -c;
    }
    // same structures, c==0
    // if (memcmp(a->s, b->s, 299 * sizeof(char)) != 0) { fmt::print("error\n"); }

    if ((a->saddle_en - b->saddle_en) != 0) return a->saddle_en - b->saddle_en;

    return a->curr_en - b->curr_en;
}

inline auto multi_findpath::compare_energy(const void* A, const void* B) -> int
{
    intermediate_multi *a, *b;

    a = (intermediate_multi*)A;
    b = (intermediate_multi*)B;

    if ((a->saddle_en - b->saddle_en) != 0) return a->saddle_en - b->saddle_en;

    return a->curr_en - b->curr_en;
}