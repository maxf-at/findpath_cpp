// this includes the basic elementary findpath algorithm

// #include <common.hpp>

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

// basic findpath
struct move_t {
    int i; /* i,j>0 insert; i,j<0 delete */
    int j;
    int when; /* 0 if still available, else resulting distance from start */
    int E;
};

struct intermediate_t {
    short*  pt;      /**<  @brief  pair table */
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

    auto try_moves(vrna_fold_compound_t* vc, intermediate_t c, int maxE, intermediate_t* next,
                   int dist, int bp_dist) -> int;
    auto find_path_once(vrna_fold_compound_t* vc, short* pt1, short* pt2, int current_search_width,
                        int maxE, bool direction, int final_search_width)
        -> std::vector<sorted_path>;
    auto findpath_method(vrna_fold_compound_t* vc, short* pt1, short* pt2, int width, int maxE,
                         bool init_direction) -> std::vector<sorted_path>;
    auto free_intermediate(intermediate_t* i) -> void;
    static auto compare_ptable(const void* A, const void* B) -> int;
    static auto compare_energy(const void* A, const void* B) -> int;
    static auto compare_moves_when(const void* A, const void* B) -> int;
    auto        copy_moves(move_t* mvs, int bp_dist) -> move_t*;

   public:
    auto init(vrna_fold_compound_t* fc, short* pt1, short* pt2, int final_search_width)
        -> std::vector<sorted_path>;
    void test();
    single_findpath(){};
};

inline void single_findpath::test() {}

inline auto single_findpath::init(vrna_fold_compound_t* fc, short* pt1, short* pt2,
                                  int final_search_width) -> std::vector<sorted_path>
{
    // fmt::print("{}\n", vrna_db_from_ptable(pt1));
    // fmt::print("{}\n", vrna_db_from_ptable(pt2));
    // fmt::print("hello fmt\n");

    auto result = findpath_method(fc, vrna_ptable_copy(pt1), vrna_ptable_copy(pt2),
                                  final_search_width, INT_MAX - 1, false);

    // auto result2 = findpath_method(fc, vrna_ptable_copy(pt1), vrna_ptable_copy(pt2),
    //                                final_search_width, INT_MAX - 1, true);
    // std::move(result2.begin(), result2.end(), std::back_inserter(result));

    // best path to [0] (lowest max_en)
    sort(result.begin(), result.end(),
         [](const auto& a, const auto& b) -> bool { return a.max_en < b.max_en; });


    // for (const auto m : result[0].moves) { fmt::print("{} {}\n", m.i, m.j); }
    fmt::print("{}\n", result[0].max_en);

    // fmt::print("test\n");
    return result;
};

inline auto single_findpath::findpath_method(vrna_fold_compound_t* vc, short* pt1, short* pt2,
                                             int width, int maxE, bool init_direction)
    -> std::vector<sorted_path>

// ub = upper bound?
{
    // short * ptr, *pt1, *pt2;
    short* ptr;

    move_t* bestpath = NULL;
    int     dir;

    std::vector<sorted_path> all_paths;

    bool path_fwd = dir = 0;

    // new private bool

    bool direction = false;  // true: s1 -> s2, false: s2 -> s1

    if (not init_direction) {
        direction = true;
        ptr       = pt1;
        pt1       = pt2;
        pt2       = ptr;
    }

    // pt1 = vrna_ptable(s1);
    // pt2 = vrna_ptable(s2);

    int              last_iteration = width;
    std::vector<int> iterations{width};

    float accelerator = 1.0;

    while (last_iteration > 16) {
        const int next_iteration = int(last_iteration / 5.0 * accelerator);
        iterations.push_back(next_iteration);
        last_iteration = next_iteration;

        accelerator *= accelerator;

        if (accelerator > 2) accelerator = 2;
    }

    int current_search_width;
    // current_search_width = 1;
    // do while: at least executed once.

    for (std::vector<int>::reverse_iterator it = iterations.rbegin(); it != iterations.rend();
         ++it) {
        current_search_width = *it;
        // IC(current_search_width);

        int saddleE;
        path_fwd  = !path_fwd;
        direction = !direction;

        // overshoot
        // if (current_search_width > width) current_search_width = width;

        // if (path) free(path);

        // the compiler should take care of copy elision here...

        std::vector<sorted_path> current_paths;
        current_paths = find_path_once(vc, pt1, pt2, current_search_width, maxE, direction, width);

        // std::future<std::vector<sorted_path>> ret1 =
        //     std::async(std::launch::async, &find_path_once, vc, pt1, pt2, current_search_width,
        //                maxE, direction, width);
        // std::future<std::vector<sorted_path>> ret2 =
        //     std::async(std::launch::async, &find_path_once, vc, pt1, pt2, current_search_width,
        //                maxE, not direction, width);

        // std::vector<sorted_path> current_paths  = ret1.get();

        // std::vector<sorted_path> current_paths2 = ret2.get();

        // did we actually receive a path / paths?

        // #ifdef IC_DEBUG
        // IC("r", direction, current_search_width, current_paths.size(), saddleE);
        // #endif

        if (current_paths.size() > 0) {
            saddleE = current_paths[0].max_en;

        } else {
            saddleE = maxE;
        }

        if (saddleE < maxE) {
            maxE = saddleE;
            // if (bestpath) free(bestpath);
        }

        ptr = pt1;
        pt1 = pt2;
        pt2 = ptr;
        // current_search_width *= 4;

        // concatenate vectors
        std::move(current_paths.begin(), current_paths.end(), std::back_inserter(all_paths));
    }

    // while (current_search_width < 4 * width);

    // #ifdef IC_DEBUG
    // IC(all_paths.size(), current_search_width, all_paths[all_paths.size() - 1].max_en);
    // #endif

    // current workaround...
    // all_paths[0].max_en = maxE;

    /* (re)set some globals */
    // path     = bestpath;
    // path_fwd = dir;

    free(pt1);
    free(pt2);

    return all_paths;
}

// private functions below

inline auto single_findpath::try_moves(vrna_fold_compound_t* vc, intermediate_t c, int maxE,
                                       intermediate_t* next, int dist, int bp_dist) -> int
{
    int *   loopidx, len, num_next = 0, en, oldE;
    move_t* mv;
    short*  pt;

    // this does not work with LOOP_EN
    len     = c.pt[0];
    loopidx = vrna_loopidx_from_ptable(c.pt);
    oldE    = c.Sen;
    for (mv = c.moves; mv->i != 0; mv++) {
        int i, j;

        auto source_1 = 0;
        auto source_2 = 0;
        auto dest_1   = 0;
        auto dest_2   = 0;

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

            // if (en < maxE) {
            next[num_next].Sen     = (en > oldE) ? en : oldE;
            next[num_next].curr_en = en;
            next[num_next].pt      = pt;
            mv->when               = dist;
            mv->E                  = en;
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
            num_next += try_moves(vc, current[c], maxE, next + num_next, d, bp_dist);
        if (num_next == 0) {
            // case where we don't find any moves -> abort
            for (cc = current; cc->pt != NULL; cc++) free_intermediate(cc);
            current[0].Sen = INT_MAX;
            break;
        }

        /* remove duplicates via sort|uniq
         * if this becomes a bottleneck we can use a hash instead */
        qsort(next, num_next, sizeof(intermediate_t), compare_ptable);

        // this shrinks the next array which makes the following sorting step faster
        for (u = 0, c = 1; c < num_next; c++) {
            if (memcmp(next[u].pt, next[c].pt, sizeof(short) * len) != 0) {
                next[++u] = next[c];

            } else {
                // new part - dont just delete valid duplicates
                // we dont get here during the greedy first iteration
                // only consider last fwd & bwd pass

                if (d >= dist && next[c].Sen <= maxE &&
                    current_search_width >= final_search_width) {
                    // if (next[c].Sen <= maxE && current_search_width>=search_width){
                    move_t* temp_moves = copy_moves(next[c].moves, bp_dist);
                    qsort(temp_moves, bp_dist, sizeof(move_t), compare_moves_when);

                    sorted_path current_path(dist);
                    // sorted_path_2 tests(dist);

                    // qsort(path, bp_dist, sizeof(move_t), compare_moves_when);
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
        qsort(next, num_next, sizeof(intermediate_t), compare_energy);

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

        // std::vector<sorted_move> current_path;

        // preallocate bp_dist elements in vector
        sorted_path current_path(dist);
        // sorted_path_2 tests(dist);

        qsort(path, bp_dist, sizeof(move_t), compare_moves_when);
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
    // return result;
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