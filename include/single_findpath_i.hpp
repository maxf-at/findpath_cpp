// this includes the basic elementary findpath algorithm



class single_findpath_i
{
   private:


    // std::vector<int> add_moves;

    static auto find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2,
                               const int current_search_width, int maxE, bool direction,
                               int final_search_width, std::vector<int> add_moves) -> std::vector<sorted_path>;
    static auto findpath_method(vrna_fold_compound_t* vc, short* pt1, short* pt2, int width,
                                int maxE, bool mp, std::vector<int> add_moves) -> std::vector<sorted_path>;
    static auto free_intermediate(intermediate_t* i) -> void;
    static auto compare_ptable(const void* A, const void* B) -> int;
    static auto compare_energy(const void* A, const void* B) -> int;
    static auto compare_moves_when(const void* A, const void* B) -> int;
    static auto ptable_from_string(short* pt, short* loop, const char* s, short* stack) -> void;

    vrna_fold_compound_t* fc;
    std::string sequence; 
    std::string s1; 
    std::string s2;
    std::vector<sorted_path> result; // results are stored here

   public:

    // bool Verbose = false;

    auto init(vrna_fold_compound_t* fc, short* pt1, short* pt2, std::vector<int> add_moves = {}, int final_search_width = 2,
              bool mp = true, int en_limit = INT_MAX - 1, bool Verbose = false) -> std::vector<sorted_path>;

    single_findpath_i(){};

    // result getters
    auto get_en() -> int {
        if (result.size() > 0) {
            return result[0].max_en;
        } else {
            return INT_MAX - 1;
        }
    }
    auto get_path() -> std::vector<std::tuple<int,int, int>> {
        
        // get list of moves, e.g.
        // [(0, 0, -18), (-33, -99, -17), (-34, -98, -14), ... ]

        std::vector<std::tuple<int,int, int>> path;
        auto const& moves = result[0].moves;
        short* pt = vrna_ptable(s1.c_str());
        float en     = vrna_eval_structure_pt(fc, pt);
        
        path.push_back({0, 0, en});

        for (auto const& move : moves) {
            if (move.i == 0)
                break; // any trailing (0,0) indicate we're done (unused indirect moves)

            if (move.j < 0) {
                pt[-move.i] = 0;
                pt[-move.j] = 0;
            } else {
                pt[move.i] = move.j;
                pt[move.j] = move.i;
            }
            en            = vrna_eval_structure_pt(fc, pt);
            path.push_back({move.i, move.j, en});
        }

        free(pt);
        return path;
    }

    // alternative constructor for Python export with combined init
    single_findpath_i(std::string sequence, std::string s1, std::string s2, std::vector<int> add_moves = {}, int search_width = 0, float search_width_multiplier=0,
                    bool mp = true, int en_limit = INT_MAX - 1, const py::dict& model_details = {}) : sequence{sequence}, s1{s1}, s2{s2} {
        


        vrna_md_t md;

        vrna_md_set_default(&md); // copy global settings
        // set_model_details(&md);

        if (model_details.contains("noLP")){ 
            // Only consider canonical structures, i.e. no 'lonely' base pairs. 
            md.noLP = model_details["noLP"].cast<int>();
        }
        if (model_details.contains("logML")){ 
            // Use logarithmic scaling for multiloops. 
            md.logML = model_details["logML"].cast<int>();
        }
        if (model_details.contains("temperature")){
            // The temperature used to scale the thermodynamic parameters. 
            // py::print("set T to", model_details["temperature"]);
            md.temperature = model_details["temperature"].cast<double>();
        }
        if (model_details.contains("dangles")){
            // Specifies the dangle model used in any energy evaluation (0,1,2 or 3) 
            md.dangles = model_details["dangles"].cast<int>();
        }
        if (model_details.contains("special_hp")){
            // Include special hairpin contributions for tri, tetra and hexaloops. 
            md.special_hp = model_details["special_hp"].cast<int>();
        }
        if (model_details.contains("noGU")){
            // Do not allow GU pairs. 
            md.noGU = model_details["noGU"].cast<int>();
        }
        if (model_details.contains("noGUclosure")){
            // Do not allow loops to be closed by GU pair. 
            md.noGUclosure = model_details["noGUclosure"].cast<int>();
        }





        fc = vrna_fold_compound(sequence.c_str(), &md, VRNA_OPTION_EVAL_ONLY);

        short* pt1 = vrna_ptable(s1.c_str());
        short* pt2 = vrna_ptable(s2.c_str());

        int bp_dist            = vrna_bp_distance(s1.c_str(), s2.c_str());

        // if we have a valid search width multiplier or nothing selected, use it:
        if (search_width_multiplier>0 or search_width==0){

            if (search_width_multiplier==0) {
                search_width_multiplier=2;
            }

            search_width = bp_dist * search_width_multiplier;
            if (search_width < 2) { search_width = 2; }
        }
        
        en_limit = INT_MAX - 1; // workaround?

        result = init(fc, pt1, pt2, add_moves, search_width, mp, en_limit, false);

        // fmt::print("init: {} {}\n", search_width, result.size());

        // free(pt1);
        // free(pt2);

    };


};


inline auto single_findpath_i::init(vrna_fold_compound_t* fc, short* pt1, short* pt2, std::vector<int> add_moves,
                                  int final_search_width, bool mp, int en_limit, bool Verbose)
    -> std::vector<sorted_path>
{
    // fmt::print("{}\n", vrna_db_from_ptable(pt1));
    // fmt::print("{}\n", vrna_db_from_ptable(pt2));

    // mp = false;
    // char* temp_s = vrna_db_from_ptable(pt1);
    // short* temp_pt = vrna_ptable_copy(pt1);
    // encode_pt(temp_pt);

    auto result = findpath_method(fc, vrna_ptable_copy(pt1), vrna_ptable_copy(pt2),
                                  final_search_width, en_limit, mp, add_moves);

    // best path to [0] (lowest max_en)

    std::sort(result.begin(), result.end(),
              [](const auto& a, const auto& b) -> bool { return a.max_en < b.max_en; });

    // fmt::print("hello fmt\n");
    // fmt::print ("s1 = '{}'\n", vrna_db_from_ptable(pt1));
    // fmt::print ("s2 = '{}'\n", vrna_db_from_ptable(pt2));

    if (Verbose) {
        for (const auto& path : result) {
            print_moves(path, fc, pt1);
            break;
            // for (const auto m : path.moves) { fmt::print("({}/{}) ", m.i, m.j); }
            //     fmt::print("| max_en: {}\n", path.max_en);
        }
    }

    return result;
};

inline auto single_findpath_i::findpath_method(vrna_fold_compound_t* fc, short* pt1, short* pt2,
                                             int final_search_width, int max_en, bool mp, std::vector<int> add_moves)
    -> std::vector<sorted_path>

{
    short*                   temp_pt;
    std::vector<sorted_path> all_paths;

    // mp                 = false;
    // final_search_width = 100;

    // new private bool
    bool direction = true;  // true: s1 -> s2, false: s2 -> s1

    short* pt1_bwd = vrna_ptable_copy(pt1);
    short* pt2_bwd = vrna_ptable_copy(pt2);

    int              last_iteration = final_search_width;
    std::vector<int> iterations{final_search_width};

    // if we have a given energy limit, only take the final search width pass

    // construct the list of iterations, e.g. [400, 80, 16]
    float accelerator = 1.0;  // this might not be useful...
    while (last_iteration > 16) {
        const int next_iteration = int(last_iteration / 5.0 * accelerator);
        iterations.push_back(next_iteration);
        last_iteration = next_iteration;
        // accelerator *= accelerator;
        // if (accelerator > 2) accelerator = 2;
    }

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
                           max_en, direction, final_search_width, add_moves);
            std::future<std::vector<sorted_path>> ret2 =
                std::async(std::launch::async, &find_path_once, fc, pt2_bwd, pt1_bwd,
                           current_search_width, max_en, not direction, final_search_width, add_moves);
            fwd_paths = ret1.get();
            bwd_paths = ret2.get();
        } else {
            fwd_paths = find_path_once(fc, pt1, pt2, current_search_width, max_en, direction,
                                       final_search_width, add_moves);
            bwd_paths = find_path_once(fc, pt2_bwd, pt1_bwd, current_search_width, max_en,
                                       not direction, final_search_width, add_moves);
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

inline auto single_findpath_i::ptable_from_string(short* pt, short* loop, const char* s, short* stack) -> void
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





// removed global direction, path

inline auto single_findpath_i::find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2,
                                            const int current_search_width, int maxE,
                                            bool direction, int final_search_width, std::vector<int> add_moves)
    -> std::vector<sorted_path>
{


    int i, d;

    int bp_dist = 0;
    int len     = static_cast<int>(pt1[0]);

    short* temp_pt = vrna_ptable_copy(pt1);
    // temp_pt[0] = pt1[0];

    short* stack     = vrna_ptable_copy(pt1);
    short* temp_loop = vrna_ptable_copy(pt1);



    // fmt::print("init: maxE: {} sw: {} dir: {}\n", maxE, current_search_width, direction);

    // generate all possible moves
    // this also calculates bp_dist...

    // generate move list. uses push backs, we don't know bp_dist (vector size) yet

    // std::vector<short> move_list{};

    std::vector<move_ij> move_list{};

    std::vector<size_t> h_list{};

    // move_list.push_back({0,0,0});

    size_t end_hash = 0;

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
                end_hash += move_hash;

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
                end_hash += move_hash;

                bp_dist++;
            }
        }
    }

    for (int i=0; i<add_moves.size(); i+=2) {
        move_list.push_back({static_cast<short>(add_moves[i]), static_cast<short>(add_moves[i+1])});
        move_list.push_back({static_cast<short>(-add_moves[i]), static_cast<short>(-add_moves[i+1])}); // undo extra move
        // fmt::print ("add {} {} \n", add_moves[i], add_moves[i+1]);
    }


    move_list.push_back({0, 0});

    int bp_dist_min = bp_dist;
    bp_dist += add_moves.size();

    // short* stack = new short[200];
    // fmt::print("input: {}\n", vrna_db_from_ptable(pt2));
    // fmt::print("input: {}\n", vrna_db_from_ptable(pt1));
    // fmt::print("input: {}\n", vrna_db_from_ptable(temp_pt));
    // fmt::print("bp_dist {}\n", bp_dist);

    // memory pool for structures and moves
    // we're using a double buffer such that current results don't get immediately overwritten in
    // the next iteration (the pointers from current vector need to be valid for the next vector)

    // std::vector<std::vector<char>> s_pool(2 * bp_dist * current_search_width + 1);
    // for (auto& col : s_pool) {
    //     col.reserve(len + 1);
    // }  // column space allocation so that we can memcpy into it

    std::vector<std::vector<char>> structure_storage(2 * current_search_width + 4);
    for (auto& col : structure_storage) {
        col.reserve(len + 1);
    }  // column space allocation so that we can memcpy into it

    // std::vector<std::vector<short>> move_pool(2 * bp_dist * current_search_width + 1);
    // for (auto& col : move_pool) { col.reserve(bp_dist + 1); }

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

    // ptable_from_string(temp_pt, temp_loop, current[0].s, stack);

    // next =
    //     (intermediate_t*)vrna_alloc(sizeof(intermediate_t) * (bp_dist * current_search_width +
    //     1));

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
 

            // update temporary pairing table & loop table to match the current structure string
            ptable_from_string(temp_pt, temp_loop, current[c].s, stack);
            // int* loopidx = vrna_loopidx_from_ptable(temp_pt);

            // fmt::print("d: {}, c: {}, s: {} / {}\n", d, c, current[c].s, current[c].s_hash);

            for (int a = 0; a < move_list.size(); a++) {
                // for (mv = current[c].moves; mv->i != 0; mv++) {

                // fmt::print ("d: {}, a: {}, num_next: {}, c: {}\n", d, a, num_next, c);

                const int move_index = current[c].moves[a];
                if (move_index > 0) {
                    continue;
                    // the move at index a was already taken
                }

                const int i = move_list[a].i;
                const int j = move_list[a].j;

                // mv = &current[c].moves[a];

                // int i, j;

                int source_1;
                int source_2;
                int dest_1;
                int dest_2;

                // if (mv->when > 0) continue;

                // i = mv->i;
                // j = mv->j;

                if (j < 0 and temp_pt[-i] == -j) {  // bp deletion
                    source_1 = -i;
                    dest_1   = 0;
                    source_2 = -j;
                    dest_2   = 0;

                } else {                                              // bp insertion
                    if (j > 0 and (temp_loop[i] == temp_loop[j]) and  // i and j belong to same loop
                        (temp_pt[i] == 0) and (temp_pt[j] == 0)       // and are unpaired
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

                // fmt::print("d: {} / eval {}/{} (move:{})\n", d, i,j,a);

                // this does not work with LOOP_EN
                en = current[c].curr_en + vrna_eval_move_pt(vc, temp_pt, i, j);
                // en = vrna_eval_structure_pt(vc, pt);

                // this used to be en < maxE
                if (en <= maxE) {
                    // pt = (short*)vrna_alloc(sizeof(short) * (len + 1));
                    // memcpy(pt, temp_pt, (len + 1) * sizeof(short));

                    // pt[source_1] = dest_1;
                    // pt[source_2] = dest_2;

                    // memcpy(s_pool[num_next + pool_offset].data(), current[c].s,
                    //        (current[c].length) * sizeof(char));
                    // if (dest_1 == 0) {
                    //     s_pool[num_next + pool_offset][source_1 - 1] = '.';
                    //     s_pool[num_next + pool_offset][source_2 - 1] = '.';
                    // } else {
                    //     s_pool[num_next + pool_offset][source_1 - 1] = '(';
                    //     s_pool[num_next + pool_offset][source_2 - 1] = ')';
                    // }
                    // next[num_next].s =
                    //     s_pool[num_next + pool_offset].data();  // set structure pointer

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

                    // size_t move_hash;
                    // move_hash = 17*std::abs(source_1);
                    // // move_hash += (move_hash << 10);
                    // // move_hash ^= (move_hash >> 6);
                    // move_hash += 31*std::abs(source_2);
                    // // move_hash += (move_hash << 8);
                    // // move_hash ^= (move_hash >> 4);

                    // next[num_next].s_hash = current[c].s_hash + move_hash;

                    // size_t move_hash;
                    // move_hash = int_hash(i);
                    // move_hash *= j;
                    // move_hash += j * 13;
                    // next[num_next].s_hash = current[c].s_hash + move_hash;

                    // size_t move_hash;
                    // move_hash = int_hash_64(i);
                    // move_hash *= j;
                    // move_hash += j * 13;
                    // next[num_next].s_hash = current[c].s_hash + move_hash;

                    next[num_next].s_hash = current[c].s_hash + h_list[a];

                    // next[num_next].s_hash      = hash_cstr(current[c].s, len-1);
                    // next[num_next].s_hash = fasthash64(current[c].s, sizeof(char)*(len-1), 42);

                    // fmt::print ("d: {}, c: {}, a: {}, s:{}\n", d, c, a, current[c].s);

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

                    // next[num_next].s    = s_pool[num_next + pool_offset];

                    next[num_next].length = current[c].length;

                    next[num_next].saddle_en =
                        (en > current[c].saddle_en) ? en : current[c].saddle_en;
                    next[num_next].curr_en = en;
                    // next[num_next].pt      = pt;
                    // mv->when = d;
                    // mv->E                  = en;

                    // memcpy(move_pool[num_next + pool_offset].data(), current[c].moves,
                    //        sizeof(short) * (bp_dist + 1));
                    // next[num_next].moves =
                    //     move_pool[num_next + pool_offset].data();  // set move pointer

                    // move_pool[num_next + pool_offset][a] = d;

                    next[num_next].last_id = c;  // set last_id (from which node we came from)
                    next[num_next].move_id =
                        a;  // set move_id to current bp dist (further down after sorting)

                    // next[num_next].moves = copy_moves(current[c].moves, bp_dist);
                    // mv->when = 0;

                    num_next++;
                }
            }

            // free(loopidx);
            // return num_next;
        }
        // fmt::print("A1 {}\n", num_next);

        if (num_next == 0) {
            // case where we don't find any moves -> abort
            for (cc = current; cc->s != nullptr; cc++) free_intermediate(cc);
            current[0].saddle_en = INT_MAX;
            break;
        }

        // fmt::print("A2\n");

        /* remove duplicates via sort|uniq
         * if this becomes a bottleneck we can use a hash instead */
        std::qsort(next, num_next, sizeof(intermediate_t), compare_ptable);

        // std::sort(next, next + num_next, [](const auto& a, const auto& b) -> bool {
        //     //  return a.max_en < b.max_en;
        //     if (a.s_hash != b.s_hash) {
        //         return a.s_hash > b.s_hash;
        //     }
        //     if (a.saddle_en != b.saddle_en) return a.saddle_en < b.saddle_en;

        //     return a.curr_en < b.curr_en;
        // });

        bool flag = true;

        if (d == bp_dist && current_search_width >= final_search_width) {
            flag = false;

            for (u = 0, c = 1; c < num_next; c++) { next[++u] = next[c]; }

            num_next = u + 1;
        }

        // this shrinks the next array which makes the following sorting step faster
        if (d <= bp_dist and flag) {
            for (u = 0, c = 1; c < num_next; c++) {
                if (next[u].s_hash != next[c].s_hash) {
                    next[++u] = next[c];

                } else {
                    // fmt::print("F\n");
                    free_intermediate(next + c);
                }
            }
            num_next = u + 1;
        }

        num_next = u + 1;
        std::qsort(next, num_next, sizeof(intermediate_t), compare_energy);

        // next is now again reduced to size current, means replace current with next.

        // fmt::print("d: {}/{} update structures / n: {}\n", d, bp_dist, num_next);

        // if (d == bp_dist) {
        //     for (u = 0; u < num_next; u++) {
        //         auto last_id = next[u].last_id;
        //         auto move_id = next[u].move_id;

        //         fmt::print("C {} {} {} {}\n", u, last_id, current[last_id].moves == nullptr,
        //                    current[last_id].saddle_en);
        //     }
        // }

        // update moves & structures for next iteration
        for (u = 0; u < current_search_width and u < num_next; u++) {
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

            // fmt::print("A {} {}\n", last_id, current[last_id].s);

            memcpy(structure_storage[u + pool_offset].data(), current[last_id].s,
                   (current[last_id].length) * sizeof(char));

            // fmt::print("done\n");

            next[u].s = structure_storage[u + pool_offset].data();  // set structure pointer

            if (dest_1 == 0) {
                next[u].s[source_1 - 1] = '.';
                next[u].s[source_2 - 1] = '.';
            } else {
                next[u].s[source_1 - 1] = '(';
                next[u].s[source_2 - 1] = ')';
            }

            if (d >= bp_dist_min and next[u].s_hash == end_hash) {
                // fmt::print("found extra, d = {}\n", d);
                if (not next[u].moves) { continue; }
                // allocate d (current dist) elements in current_path vector and set max_en


                auto& current_path = all_paths.emplace_back(bp_dist, next[u].saddle_en);
                for (int dist = 0; dist < bp_dist; dist++) {
                    int current_dist = next[u].moves[dist] - 1;
                    if (current_dist == -1) {
                        // unused indirect move
                        continue;
                    }
                    int i, j;
                    i = move_list[dist].i;
                    j = move_list[dist].j;
                    // fmt::print("dist={}, i={}, j={}, current_dist={} / ({})\n", dist, i, j,
                    //            current_dist, d);

                    if (direction) {
                        // fwd path
                        current_path.moves[current_dist] = {i, j, 0};
                    } else {
                        // bwd path, fill vector back to front
                        current_path.moves[bp_dist - current_dist - 1] = {-i, -j, 0};
                    }
                }





            }

            // fmt::print("u: {} {} {}\n", next[u].moves[0], next[u].moves[1],
            // next[u].moves[2]);
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

    // fmt::print("E1 {}\n", c);

    for (int index = 0; index < c; index++) {
        if (not current[index].moves) { continue; }

        // fmt::print("E2 {} {} {}\n", index, c, current[index].moves[0]);

        // allocate bp_dist elements in current_path vector and set max_en

        auto& current_path = all_paths.emplace_back(bp_dist, current[index].saddle_en);
        for (d = 0; d < bp_dist; d++) {
            int current_dist = current[index].moves[d] - 1;
            int i, j;
            i = move_list[d].i;
            j = move_list[d].j;
            // fmt::print("d={}, i={}, j={}, current_dist={}\n", d, i, j, current_dist);

            if (direction) {
                // fwd path
                current_path.moves[current_dist] = {i, j, 0};
            } else {
                // bwd path, fill vector back to front
                current_path.moves[bp_dist - current_dist - 1] = {-i, -j, 0};
            }
        }

    }

    // fmt::print("End, return: {} paths / end_hash: {}\n", all_paths.size(), end_hash);

    free(init_s);
    free(temp_pt);
    free(temp_loop);
    free(stack);

    return all_paths;
}

inline auto single_findpath_i::free_intermediate(intermediate_t* i) -> void
{
    i->s         = nullptr;
    i->moves     = nullptr;
    i->saddle_en = INT_MAX;
}

inline auto single_findpath_i::compare_ptable(const void* A, const void* B) -> int
{
    intermediate_t *a, *b;
    int             c;

    a = (intermediate_t*)A;
    b = (intermediate_t*)B;

    // if both structures are not identical, sort them according to the hash value
    if (a->s_hash != b->s_hash) {
        if (a->s_hash > b->s_hash) {
            return 1;
        } else {
            return -1;
        }

        return a->s_hash < b->s_hash;

        return a->s_hash - b->s_hash;
    }
    // same structures, c==0
    // if (memcmp(a->s, b->s, 299 * sizeof(char)) != 0) { fmt::print("error\n"); }

    if ((a->saddle_en - b->saddle_en) != 0) return a->saddle_en - b->saddle_en;

    return a->curr_en - b->curr_en;
}

// inline auto single_findpath::compare_ptable(const void* A, const void* B) -> int
// {
//     intermediate_t *a, *b;
//     int             c;

//     a = (intermediate_t*)A;
//     b = (intermediate_t*)B;

//     c = memcmp(a->s, b->s, a->length * sizeof(char));
//     // c = strcmp (a->s, b->s);
//     // c = memcmp(a->pt, b->pt, a->pt[0] * sizeof(short));

//     if (c != 0) {
//         // int d = memcmp(a->s, b->s, 299 * sizeof(char));

//         // if (d != c) { fmt::print("error {} {}\n", c, d); }

//         // I have no idea why a minus needs to be there...🤦‍♂️
//         return -c;
//     }
//     // same structures, c==0
//     // if (memcmp(a->s, b->s, 299 * sizeof(char)) != 0) { fmt::print("error\n"); }

//     if ((a->saddle_en - b->saddle_en) != 0) return a->saddle_en - b->saddle_en;

//     return a->curr_en - b->curr_en;
// }

inline auto single_findpath_i::compare_energy(const void* A, const void* B) -> int
{
    intermediate_t *a, *b;

    a = (intermediate_t*)A;
    b = (intermediate_t*)B;

    if ((a->saddle_en - b->saddle_en) != 0) return a->saddle_en - b->saddle_en;

    return a->curr_en - b->curr_en;
}
