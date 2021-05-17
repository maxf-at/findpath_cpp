// this includes the basic elementary findpath algorithm


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
    // short*  pt; /**<  @brief  pair table */
    char*  s;
    size_t s_hash;
    
    short* moves;
    short  length;
    int    saddle_en; /**<  @brief  saddle energy so far */
    int    curr_en;   /**<  @brief  current energy */




    int last_id;
    int move_id;  // which move will take next turn

    int move_delete;
    int move_i;
    int move_j;


};

class single_findpath
{
   private:
    // vrna_fold_compound_t* fc;
    // short*                pt1;
    // short*                pt2;
    // int                   final_search_width;

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
    static auto compare_moves_when(const void* A, const void* B) -> int;
    static auto copy_moves(move_t* mvs, int bp_dist) -> move_t*;

   public:
    auto init(vrna_fold_compound_t* fc, short* pt1, short* pt2, int final_search_width,
              bool mp = true, int en_limit = INT_MAX - 1) -> std::vector<sorted_path>;
    void test();
    single_findpath(){};
};

inline void single_findpath::test() {}

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

    // move_list.push_back({0,0,0});

    for (i = 1; i <= len; i++) {
        if (pt1[i] != pt2[i]) {
            if (i < pt1[i]) {
                /* need to delete this pair */
                // move_list.push_back(0);
                move_list.push_back({static_cast<short>(-i), static_cast<short>(-pt1[i])});
                bp_dist++;
            }

            if (i < pt2[i]) {
                /* need to insert this pair */
                // move_list.push_back(0);
                move_list.push_back({static_cast<short>(i), static_cast<short>(pt2[i])});
                bp_dist++;
            }
        }
    }
    // move_list.push_back(0);
    move_list.push_back({0, 0});

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
        col.reserve(len + 1);
    }  // column space allocation so that we can memcpy into it

    std::vector<std::vector<short>> move_pool(2 * bp_dist * current_search_width + 1);
    for (auto& col : move_pool) { col.reserve(bp_dist + 1); }

    std::vector<intermediate_t> current_vec(current_search_width + 1);
    intermediate_t*             current = current_vec.data();

    std::vector<intermediate_t> next_vec(bp_dist * current_search_width + 1);
    intermediate_t*             next = next_vec.data();

    // Initialization data point

    char* init_s = vrna_db_from_ptable(pt1);

    current[0].s         = init_s;
    current[0].length    = pt1[0];
    current[0].saddle_en = current[0].curr_en = vrna_eval_structure_pt(vc, pt1);

    current[0].s_hash = 42; // initialize first structure hash with a totally random integer

    // every intermediate saves the move index, which references i and j
    std::vector<short> empty_move_index(bp_dist + 1, 0);
    current[0].moves = empty_move_index.data();

    // ptable_from_string(temp_pt, temp_loop, current[0].s, stack);

    // next =
    //     (intermediate_t*)vrna_alloc(sizeof(intermediate_t) * (bp_dist * current_search_width +
    //     1));

    std::vector<sorted_path> all_paths;

    int c=0;

    for (d = 1; d <= bp_dist; d++) {

        std::unordered_map<size_t, int> testmap;

        /* go through the distance classes */
        int             u, num_next = 0;
        intermediate_t* cc;

        // memory pool double buffer offset
        int pool_offset = 0;
        if (d % 2 == 0) { pool_offset = bp_dist * current_search_width + 2; }

        // fmt::print("{} \n", d);

        // find valid moves (extracted from original try moves function)
        for (c = 0; current[c].s != nullptr; c++) {
            int en;
            // move_t* mv;

            // update temporary pairing table & loop table to match the current structure string
            ptable_from_string(temp_pt, temp_loop, current[c].s, stack);
            // int* loopidx = vrna_loopidx_from_ptable(temp_pt);

            for (int a = 0; a < move_list.size(); a++) {
                // for (mv = current[c].moves; mv->i != 0; mv++) {

                const int move_index = current[c].moves[a];
                if (move_index > 0) {
                    continue;
                    // the move at index a was already taken
                }

                const int i = move_list[a].i;
                const int j = move_list[a].j;

                // fmt::print("{} {} {} {}\n", i,j,when,d);

                // mv = &current[c].moves[a];

                // int i, j;

                int source_1;
                int source_2;
                int dest_1;
                int dest_2;

                // if (mv->when > 0) continue;

                // i = mv->i;
                // j = mv->j;

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
                    // pt = (short*)vrna_alloc(sizeof(short) * (len + 1));
                    // memcpy(pt, temp_pt, (len + 1) * sizeof(short));

                    // pt[source_1] = dest_1;
                    // pt[source_2] = dest_2;

                    // next[num_next].s = (char*)vrna_alloc(sizeof(char) * (current[c].length));
                    // memcpy(next[num_next].s, current[c].s, (current[c].length) *
                    // sizeof(char));

                    // if (dest_1 == 0) {
                    //     next[num_next].s[source_1 - 1] = '.';
                    //     next[num_next].s[source_2 - 1] = '.';
                    // } else {
                    //     next[num_next].s[source_1 - 1] = '(';
                    //     next[num_next].s[source_2 - 1] = ')';
                    // }

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
                    
                    
                    size_t move_hash;
                    move_hash = int_hash(source_1);
                    move_hash *= source_2;
                    move_hash += source_2*13;

                    const size_t hash = current[c].s_hash + move_hash;

                    next[num_next].s_hash = hash;

                    // next[num_next].moves = (move_t*)vrna_alloc(sizeof(move_t) * (bp_dist +
                    // 1)); memcpy(next[num_next].moves, current[c].moves, sizeof(move_t) *
                    // (bp_dist + 1));

                    memcpy(move_pool[num_next + pool_offset].data(), current[c].moves,
                           sizeof(move_t) * (bp_dist + 1));
                    next[num_next].moves =
                        move_pool[num_next + pool_offset].data();  // set move pointer

                    move_pool[num_next + pool_offset][a] = d;

                    // next[num_next].moves = copy_moves(current[c].moves, bp_dist);
                    // mv->when = 0;

                    num_next++;
                }
            }

            // free(loopidx);
            // return num_next;
        }

        if (num_next == 0) {
            // case where we don't find any moves -> abort
            // fmt::print ("abort\n");
            for (cc = current; cc->s != nullptr; cc++) free_intermediate(cc);
            current[0].saddle_en = INT_MAX;
            break;
        }

        // fmt::print("B\n");

        // remove duplicates via sort|uniq
        std::qsort(next, num_next, sizeof(intermediate_t), compare_ptable);

        // this shrinks the next array which makes the following sorting step faster
        
        // if (d<=bp_dist) 



        bool flag = true;

        if (d == bp_dist && current_search_width >= final_search_width) {

            flag = false;
            // u = num_next-1;
            c = num_next;
            // num_next = u + 1;

            // std::qsort(next, num_next, sizeof(intermediate_t), compare_energy);

            /* free the old stuff */
            // for (cc = current; cc->s != nullptr; cc++) { free_intermediate(cc); }
            // for (u = 0; u < current_search_width && u < num_next; u++) { current[u] = next[u]; }


            // break;


        }

        // this shrinks the next array which makes the following sorting step faster
        if (d <= bp_dist and flag) {

            for (u = 0, c = 1; c < num_next; c++) {
                // if (memcmp(next[u].pt, next[c].pt, sizeof(short) * len) != 0) {
                if (memcmp(next[u].s, next[c].s, sizeof(char) * len) != 0) {
                    next[++u] = next[c];

                } else {

                    // fmt::print("F\n");
                    free_intermediate(next + c);
                }
            }
            num_next = u + 1;

        }

        std::qsort(next, num_next, sizeof(intermediate_t), compare_energy);

        // next is now again reduced to size current, means replace current with next.
        /* free the old stuff */
        for (cc = current; cc->s != nullptr; cc++) { free_intermediate(cc); }

        // fmt::print("D\n");

        for (u = 0; u < current_search_width && u < num_next; u++) { current[u] = next[u]; }
        for (; u < num_next; u++) { free_intermediate(next + u); }

    }

    // else 
    if (c==0){
        c = 1;
    }
    // c = 1;
        
    for (int index=0; index<c; index++){


        
        if (not current[index].moves) {
           continue; 
           
        }

        // fmt::print("E2 {} {} {}\n", index, c, current[0].moves[0]);

        // allocate bp_dist elements in current_path vector and set max_en
        auto& current_path = all_paths.emplace_back(bp_dist, current[index].saddle_en);

        // std::qsort(current[0].moves, bp_dist, sizeof(move_t), compare_moves_when);
        for (d = 0; d < bp_dist; d++) {
            int current_dist = current[index].moves[d] - 1;

            int i, j;
            i = move_list[d].i;
            j = move_list[d].j;

            // i = 0;
            // j = 0;

            // fmt::print("d={}, i={}, j={}, current_dist={}\n", d, i, j, current_dist);

            if (direction) {
                // fwd path
                current_path.moves[current_dist] = {i, j, 0};
                // current_path.moves[d] = {i, j, 0};
            } else {
                // bwd path, fill vector back to front
                current_path.moves[bp_dist - current_dist - 1] = {-i, -j, 0};
                // current_path.moves[bp_dist - d - 1] = {-i, -j, 0};
            }
        }

    }

    // fmt::print("E3, size:{}\n", all_paths.size());
    // int in = 0;
    // for (const auto & current_path : all_paths){
    //     fmt::print("{} bp_dist: {} / max_en: {}\n", in, bp_dist, current[index].saddle_en );
    //     for (const auto & move : current_path.moves){
    //         fmt::print("({}/{}), ", move.i, move.j);
    //     }
    //     fmt::print("\n");
    //     in++;
    // }





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
    intermediate_t *a, *b;
    int             c;

    a = (intermediate_t*)A;
    b = (intermediate_t*)B;


    // if both structures are not identical, sort them according to the hash value
    if (a->s_hash != b->s_hash) {       
        if (a->s_hash > b->s_hash) {
            return 1;
        }
        else{
            return -1;
        }
    }

    // c = memcmp(a->s, b->s, a->length * sizeof(char));
    // if (c != 0) {
    //     return -c;
    // }

    // same structures, c==0
    if (memcmp(a->s, b->s,  a->length * sizeof(char)) != 0) { fmt::print("error\n"); }

    if ((a->saddle_en - b->saddle_en) != 0) return a->saddle_en - b->saddle_en;

    return a->curr_en - b->curr_en;
}

inline auto single_findpath::compare_energy(const void* A, const void* B) -> int
{
    intermediate_t *a, *b;

    a = (intermediate_t*)A;
    b = (intermediate_t*)B;

    if ((a->saddle_en - b->saddle_en) != 0) return a->saddle_en - b->saddle_en;

    return a->curr_en - b->curr_en;
}

inline auto single_findpath::compare_moves_when(const void* A, const void* B) -> int
{
    move_t *a, *b;

    a = (move_t*)A;
    b = (move_t*)B;

    return a->when - b->when;
}

inline auto single_findpath::copy_moves(move_t* mvs, int bp_dist) -> move_t*
{
    move_t* new_2;

    new_2 = (move_t*)vrna_alloc(sizeof(move_t) * (bp_dist + 1));
    memcpy(new_2, mvs, sizeof(move_t) * (bp_dist + 1));
    return new_2;
}