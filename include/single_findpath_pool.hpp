// this includes the basic elementary findpath algorithm

// #include <common.hpp>
#include <MemoryPool/robin_hood.h>

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

// move_t* copy_moves(move_t* mvs, int bp_dist);

// int compare_ptable(const void* A, const void* B);

// int compare_energy(const void* A, const void* B);

// int compare_moves_when(const void* A, const void* B);

// void free_intermediate(intermediate_t* i);

// auto find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2, int current_search_width,
//                     int maxE, bool direction, int search_width) -> std::vector<sorted_path>;

// int try_moves(vrna_fold_compound_t* fc, intermediate_t c, int maxE, intermediate_t* next, int
// dist,
//               int bp_dist);

// basic findpath structs
struct move_t {
    short i; /* i,j>0 insert; i,j<0 delete */
    short j;
    short when; /* 0 if still available, else resulting distance from start */
    // int E;
};

struct intermediate_t {
    short*  pt; /**<  @brief  pair table */
    int     pt_link;
    int     Sen;     /**<  @brief  saddle energy so far */
    int     curr_en; /**<  @brief  current energy */
    move_t* moves;   /**<  @brief  remaining moves to target */
};

class single_findpath
{
   private:
    // vrna_fold_compound_t* fc;
    // short*                pt1;
    // short*                pt2;
    // int                   final_search_width;

    static auto try_moves(vrna_fold_compound_t* vc, intermediate_t c, int maxE,
                          intermediate_t* next, int dist, int bp_dist, std::vector<std::vector<short>> &ptable_pool) -> int;
    static auto find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2,
                               int current_search_width, int maxE, bool direction,
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
              bool mp = true) -> std::vector<sorted_path>;
    void test();
    single_findpath(){};
};

inline void single_findpath::test() {}

inline auto single_findpath::init(vrna_fold_compound_t* fc, short* pt1, short* pt2,
                                  int final_search_width, bool mp) -> std::vector<sorted_path>
{
    // fmt::print("{}\n", vrna_db_from_ptable(pt1));
    // fmt::print("{}\n", vrna_db_from_ptable(pt2));
    fmt::print("hello fmt\n");
    mp = false;

    auto result = findpath_method(fc, vrna_ptable_copy(pt1), vrna_ptable_copy(pt2),
                                  final_search_width, INT_MAX - 1, mp);

    // best path to [0] (lowest max_en)
    std::sort(result.begin(), result.end(),
              [](const auto& a, const auto& b) -> bool { return a.max_en < b.max_en; });

    // for (const auto m : result[0].moves) { fmt::print("{} {}\n", m.i, m.j); }
    // fmt::print("{}\n", result[0].max_en);

    // fmt::print("test\n");
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

    float accelerator = 1.0;  // this might not be useful...

    // construct the list of iterations, e.g. [400, 80, 16]
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
                           max_en, direction, final_search_width);
            std::future<std::vector<sorted_path>> ret2 =
                std::async(std::launch::async, &find_path_once, fc, pt2_bwd, pt1_bwd,
                           current_search_width, max_en, not direction, final_search_width);
            fwd_paths = ret1.get();
            bwd_paths = ret2.get();
        } else {
            fwd_paths = find_path_once(fc, pt1, pt2, current_search_width, max_en, direction,
                                       final_search_width);
            // bwd_paths = find_path_once(fc, pt2_bwd, pt1_bwd, current_search_width, max_en,
            //                            not direction, final_search_width);
        }

        if (fwd_paths.size() > 0) {
            // for (const auto m : fwd_paths[0].moves) { fmt::print("{} {} / ", m.i, m.j); }
            std::move(fwd_paths.begin(), fwd_paths.end(),
                      std::back_inserter(all_paths));  // concatenate paths
            saddle_en_1 = fwd_paths[0].max_en;
        }

        // if (bwd_paths.size() > 0) {
        //     // std::cout << "\n";
        //     // for (const auto m : bwd_paths[0].moves) { fmt::print("{} {} / ", m.i, m.j); }
        //     std::move(bwd_paths.begin(), bwd_paths.end(),
        //               std::back_inserter(all_paths));  // concatenate paths
        //     saddle_en_2 = bwd_paths[0].max_en;
        // }

        // std::cout << "\n" << current_search_width << " / " << saddle_en_1 << " / " << saddle_en_2
        // << " / "
        //           << max_en << "\n";

        // set max_en to minimum of both passes
        if (saddle_en_1 < max_en) { max_en = saddle_en_1; }
        // if (saddle_en_2 < max_en) { max_en = saddle_en_2; }
    }

    free(pt1);
    free(pt2);
    free(pt1_bwd);
    free(pt2_bwd);

    return all_paths;
}

// private functions below
robin_hood::unordered_map<unsigned int, std::string> my_unordered_map; 

inline auto single_findpath::try_moves(vrna_fold_compound_t* vc, intermediate_t c, int maxE,
                                       intermediate_t* next, int dist, int bp_dist, std::vector<std::vector<short>> &ptable_pool) -> int
{
    int *   loopidx, len, num_next = 0, en, oldE;
    move_t* mv;
    short*  pt;

    // std::unordered_map<std::string, unsigned int> my_unordered_map;
    

    // unordered_map<const char*, const char*> umap;
    // umap["banana"] = "yellow";

    // this does not work with LOOP_EN
    len     = c.pt[0];
    loopidx = vrna_loopidx_from_ptable(c.pt);
    oldE    = c.Sen;
    for (mv = c.moves; mv->i != 0; mv++) {
        int i, j;

        int source_1 = 0;
        int source_2 = 0;
        int dest_1   = 0;
        int dest_2   = 0;

        if (mv->when > 0) continue;

        i = mv->i;
        j = mv->j;
        // pt = (short*)vrna_alloc(sizeof(short) * (len + 1));
        // memcpy(pt, c.pt, (len + 1) * sizeof(short));
        if (j < 0) {
            /*it's a delete move */
            source_1 = -i;
            dest_1   = 0;
            source_2 = -j;
            dest_2   = 0;

        } else {
            /* insert move */
            if ((loopidx[i] == loopidx[j]) &&    /* i and j belong to same loop */
                (c.pt[i] == 0) && (c.pt[j] == 0) /* ... and are unpaired */
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

        en = c.curr_en + vrna_eval_move_pt(vc, c.pt, i, j);
        // en = vrna_eval_structure_pt(vc, pt);

        // this used to be en < maxE
        if (en <= maxE) {
            pt = (short*)vrna_alloc(sizeof(short) * (len + 1));
            memcpy(pt, c.pt, (len + 1) * sizeof(short));

            


            pt[source_1] = dest_1;
            pt[source_2] = dest_2;

            my_unordered_map[num_next] = "test";
            // memcpy(ptable_pool[num_next].data(), pt, (len + 1) * sizeof(short));


            // if (en < maxE) {
            next[num_next].Sen     = (en > oldE) ? en : oldE;
            next[num_next].curr_en = en;
            next[num_next].pt      = pt;
            mv->when               = dist;
            // mv->E                  = en;
            next[num_next++].moves = copy_moves(c.moves, bp_dist);
            mv->when               = 0;
        }
    }
    free(loopidx);
    return num_next;
}

// removed global direction, path

inline auto single_findpath::find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2,
                                            int current_search_width, int maxE, bool direction,
                                            int final_search_width) -> std::vector<sorted_path>
{
    move_t*         mlist;
    int             i, len, d, dist = 0, result;
    short*          pt;
    intermediate_t *current, *next;  // array of current, array of next

    len = (int)pt1[0];
    pt  = vrna_ptable_copy(pt1);

    mlist = (move_t*)vrna_alloc(sizeof(move_t) * len); /* bp_dist < n */

    // generate all possible moves
    // this also calculates bp_dist...

    for (i = 1; i <= len; i++) {
        if (pt[i] != pt2[i]) {
            if (i < pt[i]) {
                /* need to delete this pair */
                mlist[dist].i      = -i;
                mlist[dist].j      = -pt[i];
                mlist[dist++].when = 0;
            }

            if (i < pt2[i]) {
                /* need to insert this pair */
                mlist[dist].i      = i;
                mlist[dist].j      = pt2[i];
                mlist[dist++].when = 0;
            }
        }
    }

    int bp_dist = dist;  // this is now local...

    // #ifdef IC_DEBUG
    // IC(direction, dist, current_search_width, maxE, std::this_thread::get_id());
    // #endif

    // std::vector<short> vpt(pt[0] + 1);

    std::vector<std::vector<short>> ptable_pool(current_search_width+bp_dist+1, std::vector<short>(len + 1));

    // vector<vector<int>> table(3, vector<int>(4));

    // ptable_pool.resize(10000);

    // ptable_pool[0] = {};
    memcpy(ptable_pool[0].data(), pt, (len + 1) * sizeof(short));

    // memcpy(ptable_pool[10000].data(), pt, (len + 1) * sizeof(short));



    // ptable_pool[0] = vpt;
    // std::cout << ptable_pool.size() << "\n";
    // ptable_pool[1] = vpt;
    // std::cout << ptable_pool.size() << "\n";

    // for (const short &a : ptable_pool[0]){
    //     std::cout << a << " ";
    // }

    // std::cout << std::endl;

    current = (intermediate_t*)vrna_alloc(sizeof(intermediate_t) * (current_search_width + 1));
    current[0].pt  = pt;
    current[0].Sen = current[0].curr_en = vrna_eval_structure_pt(vc, pt);
    current[0].moves                    = mlist;
    next = (intermediate_t*)vrna_alloc(sizeof(intermediate_t) * (dist * current_search_width + 1));

    std::vector<sorted_path> all_paths;

    for (d = 1; d <= dist; d++) {
        /* go through the distance classes */
        int             c, u, num_next = 0;
        intermediate_t* cc;

        for (c = 0; current[c].pt != NULL; c++)
            // this fills up teh next array of intermediates
            num_next += try_moves(vc, current[c], maxE, next + num_next, d, bp_dist, ptable_pool);
        if (num_next == 0) {
            // case where we don't find any moves -> abort
            for (cc = current; cc->pt != NULL; cc++) free_intermediate(cc);
            current[0].Sen = INT_MAX;
            break;
        }

        /* remove duplicates via sort|uniq
         * if this becomes a bottleneck we can use a hash instead */
        std::qsort(next, num_next, sizeof(intermediate_t), compare_ptable);

        // this shrinks the next array which makes the following sorting step faster
        for (u = 0, c = 1; c < num_next; c++) {
            if (memcmp(next[u].pt, next[c].pt, sizeof(short) * len) != 0) {
                next[++u] = next[c];

            } else {
                // new part - save multiple result paths
                // we dont get here during first passes...
                // only consider last fwd & bwd pass

                if (d >= dist && next[c].Sen <= maxE &&
                    current_search_width >= final_search_width) {
                    // if (next[c].Sen <= maxE && current_search_width>=search_width){
                    move_t* temp_moves = copy_moves(next[c].moves, bp_dist);
                    std::qsort(temp_moves, bp_dist, sizeof(move_t), compare_moves_when);

                    // preallocate bp_dist elements in vector
                    sorted_path current_path(dist);

                    for (d = 0; d < bp_dist; d++) {
                        int i, j;
                        i = temp_moves[d].i;
                        j = temp_moves[d].j;

                        if (direction) {
                            // fwd path
                            current_path.moves[d] = {i, j, 0};
                        } else {
                            // bwd path, fill vector back to front
                            current_path.moves[dist - d - 1] = {-i, -j, 0};
                        }
                    }
                    current_path.max_en = next[c].Sen;
                    all_paths.push_back(current_path);
                }

                free_intermediate(next + c);
            }
        }

        // IC(num_next);

        num_next = u + 1;
        std::qsort(next, num_next, sizeof(intermediate_t), compare_energy);

        // next is now again reduced to size current, means replace current with next.
        /* free the old stuff */
        for (cc = current; cc->pt != NULL; cc++) free_intermediate(cc);
        for (u = 0; u < current_search_width && u < num_next; u++) current[u] = next[u];
        for (; u < num_next; u++) free_intermediate(next + u);
        num_next = 0;
    }
    free(next);

    // generate temporary path object...
    move_t* path = nullptr;

    path   = current[0].moves;
    result = current[0].Sen;
    free(current[0].pt);
    free(current);

    if (path) {
        // IC("path", direction);

        // preallocate bp_dist elements in vector
        sorted_path current_path(dist);

        std::qsort(path, bp_dist, sizeof(move_t), compare_moves_when);
        for (d = 0; d < bp_dist; d++) {
            int i, j;
            i = path[d].i;
            j = path[d].j;

            if (direction) {
                // fwd path
                current_path.moves[d] = {i, j, 0};
            } else {
                // bwd path, fill vector back to front
                current_path.moves[dist - d - 1] = {-i, -j, 0};
            }
        }

        current_path.max_en = result;
        all_paths.push_back(current_path);

        // IC(all_paths[0].moves.size(), all_paths[0].moves[0].i);
        free(path);
    }

    // IC(all_paths.size(), result, current_search_width);
    return all_paths;
}

inline auto single_findpath::free_intermediate(intermediate_t* i) -> void
{
    free(i->pt);
    free(i->moves);
    i->pt    = NULL;
    i->moves = NULL;
    i->Sen   = INT_MAX;
}

inline auto single_findpath::compare_ptable(const void* A, const void* B) -> int
{
    intermediate_t *a, *b;
    int             c;

    a = (intermediate_t*)A;
    b = (intermediate_t*)B;

    c = memcmp(a->pt, b->pt, a->pt[0] * sizeof(short));
    if (c != 0) return c;

    if ((a->Sen - b->Sen) != 0) return a->Sen - b->Sen;

    return a->curr_en - b->curr_en;
}

inline auto single_findpath::compare_energy(const void* A, const void* B) -> int
{
    intermediate_t *a, *b;

    a = (intermediate_t*)A;
    b = (intermediate_t*)B;

    if ((a->Sen - b->Sen) != 0) return a->Sen - b->Sen;

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