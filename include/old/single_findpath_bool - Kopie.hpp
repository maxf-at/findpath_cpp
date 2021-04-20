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

#include <FastMemcpy/FastMemcpy_Avx.h>
// #include <FastMemcpy/FastMemcpy.h>

// #include <assert.h>
// #include <immintrin.h>
// #include <cstdint>
// /* ... */
// void fastMemcpy(void *pvDest, void *pvSrc, size_t nBytes) {
//   assert(nBytes % 32 == 0);
//   assert((intptr_t(pvDest) & 31) == 0);
//   assert((intptr_t(pvSrc) & 31) == 0);
//   const __m256i *pSrc = reinterpret_cast<const __m256i*>(pvSrc);
//   __m256i *pDest = reinterpret_cast<__m256i*>(pvDest);
//   int64_t nVects = nBytes / sizeof(*pSrc);
//   for (; nVects > 0; nVects--, pSrc++, pDest++) {
//     const __m256i loaded = _mm256_stream_load_si256(pSrc);
//     _mm256_stream_si256(pDest, loaded);
//   }
//   _mm_sfence();
// }

double t_a      = 0;
double t_b      = 0;
double t_c      = 0;
double t_d      = 0;
double t_e      = 0;
double memory_a = 0;
double memory_b = 0;

// basic findpath structs
// struct move_t {
//     int i; /* i,j>0 insert; i,j<0 delete */
//     int j;
//     int when; /* 0 if still available, else resulting distance from start */
//     int E;
// };

// struct move_t {
//     short i; /* i,j>0 insert; i,j<0 delete */
//     short j;
//     short when; /* 0 if still available, else resulting distance from start */
//     short E;
// };

struct move_t {
    int16_t i; /* i,j>0 insert; i,j<0 delete */
    int16_t j;
    int16_t when; /* 0 if still available, else resulting distance from start */
    int16_t E;
};

// struct move_t {
//     int i; /* i,j>0 insert; i,j<0 delete */
//     int j;
//     int when; /* 0 if still available, else resulting distance from start */
//     int E;
// };

// struct intermediate_t {
//     short*  pt;      /**<  @brief  pair table */
//     int     Sen;     /**<  @brief  saddle energy so far */
//     int     curr_en; /**<  @brief  current energy */
//     move_t* moves;   /**<  @brief  remaining moves to target */
// };

struct intermediate_t {
    // std::vector<std::vector<bool>>* pt_bool;
    // int pt_id;

    // std::vector<bool> m_allFalse;
    short   pt_a[201];
    short*  pt;      /**<  @brief  pair table */
    int     Sen;     /**<  @brief  saddle energy so far */
    int     curr_en; /**<  @brief  current energy */
    move_t* moves;   /**<  @brief  remaining moves to target */
};

struct intermediate_new {
    // std::vector<std::vector<bool>>* pt_bool;
    // int pt_id;

    // std::vector<bool> m_allFalse;
    short   pt_a[201];
    // short*  pt;      /**<  @brief  pair table */
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
                          intermediate_t* next, int dist, int bp_dist) -> int;
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
    // fmt::print("hello fmt\n");

    auto result = findpath_method(fc, vrna_ptable_copy(pt1), vrna_ptable_copy(pt2),
                                  final_search_width, INT_MAX - 1, mp);

    // fmt::print("runtimes: A: {}\n", t_a / 1000);
    // fmt::print("runtimes: B: {}\n", t_b / 1000);
    // fmt::print("runtimes: C: {}\n", t_c / 1000);
    // fmt::print("runtimes: D: {}\n", t_d / 1000);
    // fmt::print("runtimes: E: {}\n", t_e / 1000);

    // fmt::print("{}:[{}, {}, {}, {}, {}, {}, {}],\n", final_search_width, t_a, t_b, t_c, t_d, t_e,
    // memory_a, sizeof(move_t));

    fmt::print("{}:[{}, {}, {}],\n", final_search_width, t_a / 1000000,
               (t_b + t_c + t_d + t_e) / 1000000, memory_a / 1000000, memory_b / 1000000);

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

    float accelerator = 1.0;

    // construct the list of iterations, e.g. [400, 80, 16]
    while (last_iteration > 16) {
        const int next_iteration = int(last_iteration / 5.0 * accelerator);
        iterations.push_back(next_iteration);
        last_iteration = next_iteration;

        accelerator *= accelerator;

        if (accelerator > 2) accelerator = 2;
    }

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

        // std::cout << "\n" << current_search_width << " / " << saddle_en_1 << " / " << saddle_en_2
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

// inline auto single_findpath::findpath_method(vrna_fold_compound_t* fc, short* pt1, short* pt2,
//                                              int final_search_width, int max_en, bool mp)
//     -> std::vector<sorted_path>

// {

//     short* temp_pt;

//     std::vector<sorted_path> all_paths;

//     // new private bool
//     bool direction = true;  // true: s1 -> s2, false: s2 -> s1

//     // if (not init_direction) {
//     //     direction = true;
//     //     temp_pt   = pt1;
//     //     pt1       = pt2;
//     //     pt2       = temp_pt;
//     // }

//     // pt1 = vrna_ptable(s1);
//     // pt2 = vrna_ptable(s2);

//     int              last_iteration = final_search_width;
//     std::vector<int> iterations{final_search_width};

//     float accelerator = 1.0;

//     while (last_iteration > 16) {
//         const int next_iteration = int(last_iteration / 5.0 * accelerator);
//         iterations.push_back(next_iteration);
//         last_iteration = next_iteration;

//         accelerator *= accelerator;

//         if (accelerator > 2) accelerator = 2;
//     }

//     int current_search_width;

//         // std::future<std::vector<sorted_path>> ret1 =
//         //     std::async(std::launch::async, &find_path_once, vc, pt1, pt2,
//         current_search_width,
//         //                maxE, direction, width);
//         // std::future<std::vector<sorted_path>> ret2 =
//         //     std::async(std::launch::async, &find_path_once, vc, pt1, pt2,
//         current_search_width,
//         //                maxE, not direction, width);

//         // std::vector<sorted_path> current_paths  = ret1.get();
//         // std::vector<sorted_path> current_paths2 = ret2.get();

//     for (std::vector<int>::reverse_iterator it = iterations.rbegin(); it != iterations.rend();
//          ++it) {

//         current_search_width = *it;
//         int saddle_en_1;
//         int saddle_en_2;
//         std::vector<sorted_path> current_paths;

//         current_paths =
//             find_path_once(fc, pt1, pt2, current_search_width, max_en, direction,
//             final_search_width);
//         if (current_paths.size() > 0) {
//             std::move(current_paths.begin(), current_paths.end(),
//                       std::back_inserter(all_paths));  // concatenate paths
//             saddle_en_1 = current_paths[0].max_en;
//         } else {
//             saddle_en_1 = max_en;
//         }

//         direction = not direction;
//         temp_pt   = pt1;
//         pt1       = pt2;
//         pt2       = temp_pt;

//         current_paths =
//             find_path_once(fc, pt1, pt2, current_search_width, max_en, not direction,
//             final_search_width);
//         if (current_paths.size() > 0) {
//             std::move(current_paths.begin(), current_paths.end(),
//                       std::back_inserter(all_paths));  // concatenate paths
//             saddle_en_2 = current_paths[0].max_en;
//         } else {
//             saddle_en_2 = max_en;
//         }

//         direction = not direction;
//         temp_pt   = pt1;
//         pt1       = pt2;
//         pt2       = temp_pt;

//         // std::cout << current_search_width << " / " << saddle_en_1 << " / " << saddle_en_2 << "
//         /

//         // set max_en to minimum of both passes
//         if (saddle_en_1 < max_en) { max_en = saddle_en_1; }
//         if (saddle_en_2 < max_en) { max_en = saddle_en_2; }
//     }

//     free(pt1);
//     free(pt2);

//     return all_paths;
// }

// inline auto single_findpath::findpath_method(vrna_fold_compound_t* vc, short* pt1, short* pt2,
//                                              int width, int maxE, bool init_direction)
//     -> std::vector<sorted_path>

// // ub = upper bound?
// {
//     // short * ptr, *pt1, *pt2;
//     short* ptr;

//     move_t* bestpath = NULL;
//     int     dir;

//     std::vector<sorted_path> all_paths;

//     bool path_fwd = dir = 0;

//     // new private bool

//     bool direction = false;  // true: s1 -> s2, false: s2 -> s1

//     if (not init_direction) {
//         direction = true;
//         ptr       = pt1;
//         pt1       = pt2;
//         pt2       = ptr;
//     }

//     // pt1 = vrna_ptable(s1);
//     // pt2 = vrna_ptable(s2);

//     int              last_iteration = width;
//     std::vector<int> iterations{width};

//     float accelerator = 1.0;

//     while (last_iteration > 16) {
//         const int next_iteration = int(last_iteration / 5.0 * accelerator);
//         iterations.push_back(next_iteration);
//         last_iteration = next_iteration;

//         accelerator *= accelerator;

//         if (accelerator > 2) accelerator = 2;
//     }

//     int current_search_width;
//     // current_search_width = 1;
//     // do while: at least executed once.

//     for (std::vector<int>::reverse_iterator it = iterations.rbegin(); it != iterations.rend();
//          ++it) {
//         current_search_width = *it;
//         // IC(current_search_width);

//         int saddleE;
//         path_fwd  = !path_fwd;
//         direction = !direction;

//         // overshoot
//         // if (current_search_width > width) current_search_width = width;

//         // if (path) free(path);

//         // the compiler should take care of copy elision here...

//         std::vector<sorted_path> current_paths;
//         current_paths = find_path_once(vc, pt1, pt2, current_search_width, maxE, direction,
//         width);

//         // std::future<std::vector<sorted_path>> ret1 =
//         //     std::async(std::launch::async, &find_path_once, vc, pt1, pt2,
//         current_search_width,
//         //                maxE, direction, width);
//         // std::future<std::vector<sorted_path>> ret2 =
//         //     std::async(std::launch::async, &find_path_once, vc, pt1, pt2,
//         current_search_width,
//         //                maxE, not direction, width);

//         // std::vector<sorted_path> current_paths  = ret1.get();

//         // std::vector<sorted_path> current_paths2 = ret2.get();

//         // did we actually receive a path / paths?

//         // #ifdef IC_DEBUG
//         // IC("r", direction, current_search_width, current_paths.size(), saddleE);
//         // #endif

//         if (current_paths.size() > 0) {
//             saddleE = current_paths[0].max_en;

//         } else {
//             saddleE = maxE;
//         }

//         if (saddleE < maxE) {
//             maxE = saddleE;
//             // if (bestpath) free(bestpath);
//         }

//         ptr = pt1;
//         pt1 = pt2;
//         pt2 = ptr;
//         // current_search_width *= 4;

//         // concatenate vectors
//         std::move(current_paths.begin(), current_paths.end(), std::back_inserter(all_paths));
//     }

//     // while (current_search_width < 4 * width);

//     // #ifdef IC_DEBUG
//     // IC(all_paths.size(), current_search_width, all_paths[all_paths.size() - 1].max_en);
//     // #endif

//     // current workaround...
//     // all_paths[0].max_en = maxE;

//     /* (re)set some globals */
//     // path     = bestpath;
//     // path_fwd = dir;

//     free(pt1);
//     free(pt2);

//     return all_paths;
// }

// private functions below

auto ptable_to_bool(short* pt, std::vector<bool>& pt_bool)
{
    // pt_bool.assign(200, false);
    pt_bool.resize(200);

    for (int i = 1; i <= pt[0]; i++) {
        // const int j = i * 2;
        if (pt[i] > i) {
            pt_bool[i * 2]     = true;
            pt_bool[i * 2 + 1] = true;
        } else if (pt[i] < i and pt[i] != 0) {
            pt_bool[i * 2]     = true;
            pt_bool[i * 2 + 1] = false;
        }

        // fmt::print("{}>{}/{} ", i, pt[i], pt_bool[i]);
    }
}

auto print_ptable_bool(std::vector<bool>& pt_bool) -> void
{
    std::string s = "";
    // the first 2 bits are unused
    for (int a = 2; a <= pt_bool.size(); a += 2) {
        if (pt_bool[a]) {
            if (pt_bool[a + 1]) {
                s += '(';
            } else {
                s += ')';
            }
        } else {
            s += '.';
        }
    }
    std::cout << s << std::endl;
    // return s;
}

auto bool_to_ptable(std::vector<bool>& pt_bool)
{
    std::string s = "";
    // the first 2 bits are unused
    for (int a = 2; a <= pt_bool.size(); a += 2) {
        if (pt_bool[a]) {
            if (pt_bool[a + 1]) {
                s += '(';
            } else {
                s += ')';
            }
        } else {
            s += '.';
        }
    }
    // std::cout << s.size() << std::endl;
    return vrna_ptable_from_string(s.c_str(), 100);

    // from extract_pairs function

    // const char*  ptr;
    // char         open, close;
    // short*       stack;
    // unsigned int i, j, n;
    // int          hx;

    // n     = (unsigned int)pt[0];
    // stack = (short*)vrna_alloc(sizeof(short) * (n + 1));

    // open  = pair[0];
    // close = pair[1];

    // for (hx = 0, i = 1, ptr = structure; (i <= n) && (*ptr != '\0'); ptr++, i++) {
    //     if (*ptr == open) {
    //         stack[hx++] = i;
    //     } else if (*ptr == close) {
    //         j = stack[--hx];

    //         if (hx < 0) {
    //             vrna_message_warning(
    //                 "%s\nunbalanced brackets '%2s' found while extracting base pairs", structure,
    //                 pair);
    //             free(stack);
    //             return 0;
    //         }

    //         pt[i] = j;
    //         pt[j] = i;
    //     }
    // }

    // free(stack);

    // if (hx != 0) {
    //     vrna_message_warning("%s\nunbalanced brackets '%2s' found while extracting base pairs",
    //                          structure, pair);
    //     return 0;
    // }

    // return 1; /* success */
}

// global
// std::vector<std::vector<bool>> bool_vector(200000);

inline auto single_findpath::try_moves(vrna_fold_compound_t* vc, intermediate_t c, int maxE,
                                       intermediate_t* next, int dist, int bp_dist) -> int
{
    int *   loopidx, len, num_next = 0, en, oldE;
    move_t* mv;
    short*  pt;

    // this does not work with LOOP_EN
    len = c.pt_a[0];

    loopidx = vrna_loopidx_from_ptable(&c.pt_a[0]);

    // short* temp = (short*)vrna_alloc(sizeof(short) * (len + 1));
    // temp        = bool_to_ptable(c.pt_bool);
    // free(temp);

    oldE = c.Sen;
    for (mv = c.moves; mv->i != 0; mv++) {
        int i, j;

        short source_1 = 0;
        short source_2 = 0;
        short dest_1   = 0;
        short dest_2   = 0;

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
            if ((loopidx[i] == loopidx[j]) &&        /* i and j belong to same loop */
                (c.pt_a[i] == 0) && (c.pt_a[j] == 0) /* ... and are unpaired */
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

        en = c.curr_en + vrna_eval_move_pt(vc, &c.pt_a[0], i, j);

        // en = vrna_eval_structure_pt(vc, pt);

        // this used to be en < maxE
        if (en <= maxE) {
            auto start = std::chrono::system_clock::now();

            // pt = (short*)vrna_alloc(sizeof(short) * (len + 1));
            // memcpy_fast(pt, c.pt, (len + 1) * sizeof(short));  // !!!

            // new part
            memcpy_fast(next[num_next].pt_a, &c.pt_a[0], (len + 1) * sizeof(short));
            next[num_next].pt_a[source_1] = dest_1;
            next[num_next].pt_a[source_2] = dest_2;

            auto end     = std::chrono::system_clock::now();
            auto elapsed = end - start;
            memory_a += elapsed.count();

            // pt[source_1] = dest_1;
            // pt[source_2] = dest_2;

            // if (en < maxE) {
            next[num_next].Sen     = (en > oldE) ? en : oldE;
            next[num_next].curr_en = en;
            next[num_next].pt      = &next[num_next].pt_a[0];
            mv->when               = dist;
            mv->E                  = en;

            // std::copy(c.pt_bool.begin(), c.pt_bool.end(),
            //           std::back_inserter(bool_vector[num_next]));
            // if (dest_1 == 0) {
            //     bool_vector[num_next][source_1] = false;
            //     bool_vector[num_next][source_2] = false;
            //     // next[num_next].pt_bool[-source_1] = false;
            //     // next[num_next].pt_bool[-source_2] = false;
            // } else {
            //     bool_vector[num_next][source_1] = true;
            //     bool_vector[num_next][source_2] = true;
            // }

            // next[num_next].pt_bool.assign(1, true);
            // ptable_to_bool(pt, next[num_next].pt_bool);

            start = std::chrono::system_clock::now();

            next[num_next++].moves = copy_moves(c.moves, bp_dist);  // !!!

            end     = std::chrono::system_clock::now();
            elapsed = end - start;
            memory_b += elapsed.count();

            mv->when = 0;
        }
    }
    free(loopidx);
    return num_next;
}

// removed global direction, path

auto compare_bool_ptables(const void* A, const void* B) -> int
{
    std::vector<bool>*a, *b;
    int               c;

    a = (std::vector<bool>*)A;
    b = (std::vector<bool>*)B;

    // std::vector<bool> a1 = *(std::vector<bool>*)A;
    // std::vector<bool> b1 = *(std::vector<bool>*)B;

    // std::cout << "compare " << *a[0] << " / " << *b[0] << "\n";

    // return a1 == b1;

    c = memcmp(a, b, 50);  // &a[0] *
    if (c != 0) {
        std::cout << "return" << c << "\n";
        return c;
    }

    std::cout << "return 5 \n";

    return 0;
    // if ((a->Sen - b->Sen) != 0) return a->Sen - b->Sen;

    // return a->curr_en - b->curr_en;
}

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
    current[0].pt = pt;

    // current[0].pt_a  = vrna_ptable_copy(pt1);

    memcpy_fast(current[0].pt_a, pt, sizeof(short) * 200);  // new, old, size

    // std::cout << "c/c " << sizeof(current[0]) << "/" << vrna_eval_structure_pt(vc, pt) << " "
    //           << vrna_eval_structure_pt(vc, &current[0].pt_a[0]) << "\n";

    current[0].Sen = current[0].curr_en = vrna_eval_structure_pt(vc, pt);
    current[0].moves                    = mlist;
    next = (intermediate_t*)vrna_alloc(sizeof(intermediate_t) * (dist * current_search_width + 1));

    std::vector<std::vector<bool>> bool_vector(dist * current_search_width + 1);
    // bool_vector.resize(dist * current_search_width + 1);

    for (int a = 0; a < dist * current_search_width + 1; a++) {
        bool_vector[a].resize(200);
        bool_vector[a][0] = true;
    }

    // current[0].pt_bool = &bool_vector;

    // for (int a=0; a<dist * current_search_width + 1; a++){
    //     next[a].pt_bool.resize(200);
    // }

    // ptable_to_bool(pt, &current[0].pt_bool);
    // print_ptable_bool(current[0].pt_bool);

    std::vector<sorted_path> all_paths;

    int current_elements = 1;

    for (d = 1; d <= dist; d++) {
        /* go through the distance classes */
        int             c, u, num_next = 0;
        intermediate_t* cc;

        auto start = std::chrono::system_clock::now();

        // short* temp_pt = vrna_ptable_copy(pt);

        // for (c = 0; current[c].pt != NULL; c++) {
        for (c = 0; c < current_elements; c++) {
            // this fills up teh next array of intermediates

            // current[c].pt_bool.assign(200, false);

            // ptable_to_bool(temp_pt, current[c].pt_bool);

            num_next += try_moves(vc, current[c], maxE, next + num_next, d, bp_dist);
        }

        auto end     = std::chrono::system_clock::now();
        auto elapsed = end - start;
        t_a += elapsed.count();
        start = std::chrono::system_clock::now();

        if (num_next == 0) {
            // case where we don't find any moves -> abort

            // todo...

            // for (cc = current; cc->pt != NULL; cc++) {
            //     free_intermediate(cc);
            //      }

            current[0].Sen = INT_MAX;
            break;
        }

        /* remove duplicates via sort|uniq
         * if this becomes a bottleneck we can use a hash instead */
        std::qsort(next, num_next, sizeof(intermediate_t), compare_ptable);

        // std::qsort(&bool_vector[0], num_next, sizeof(bool_vector), compare_bool_ptables);

        // c = memcmp(a->pt, b->pt, a->pt[0] * sizeof(short));
        // if (c != 0) return c;

        // if ((a->Sen - b->Sen) != 0) return a->Sen - b->Sen;

        // return a->curr_en - b->curr_en;

        // std::sort(bool_vector.begin(), bool_vector.begin()+(num_next),
        //         [](const auto& a, const auto& b) -> bool {
        //                 bool c;
        //                 c = memcmp(&a, &b, a[0] * 20);
        //                 if (c != 0) return false;
        //                 return true;
        //             // return a==b;
        //              });

        // sort for unique pairing tables

        // compare_bool_ptables);

        end     = std::chrono::system_clock::now();
        elapsed = end - start;
        t_b += elapsed.count();
        start = std::chrono::system_clock::now();

        // this shrinks the next array which makes the following sorting step faster

        int redundant_entries = 0;

        for (u = 0, c = 1; c < num_next; c++) {
            // if memcmp is 0, both are identical
            // if (memcmp(next[u].pt, next[c].pt, sizeof(short) * len) != 0) {
            if (memcmp(&next[u].pt_a[0], &next[c].pt_a[0], sizeof(short) * len) != 0) {
                // both consecutive elements are different
                u++;
                next[u] = next[c];

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

                // redundant_entries++;

                free_intermediate(next + c);
            }
        }

        // std::cout << "red: " << redundant_entries << " " << c << "\n";

        end     = std::chrono::system_clock::now();
        elapsed = end - start;
        t_c += elapsed.count();
        start = std::chrono::system_clock::now();

        // IC(num_next);

        num_next = u + 1;
        std::qsort(next, num_next, sizeof(intermediate_t), compare_energy);

        end     = std::chrono::system_clock::now();
        elapsed = end - start;
        t_d += elapsed.count();
        start = std::chrono::system_clock::now();

        // next is now again reduced to size current, means replace current with next.
        /* free the old stuff */
        // for (cc = current; cc->pt != NULL; cc++) {
        
        for (int d=0; d < current_elements; d++) {
            // empty the current vector
            free_intermediate(current+d);
        }

        // std::cout << "cc c" << cc << "/" << current_elements << "\n";

        // copy next vector to current one
        current_elements = 0;
        for (u = 0; u < current_search_width && u < num_next; u++) {
            current[u] = next[u];
            current_elements++;
        }
        for (; u < num_next; u++) free_intermediate(next + u);
        num_next = 0;

        end     = std::chrono::system_clock::now();
        elapsed = end - start;
        t_e += elapsed.count();
    }
    free(next);

    // generate temporary path object...
    move_t* path = nullptr;

    path   = current[0].moves;
    result = current[0].Sen;

    // free(current[0].pt);
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
    // free(i->pt);
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

    // c = memcmp(a->pt, b->pt, a->pt[0] * sizeof(short));
    c = memcmp(&a->pt_a[0], &b->pt_a[0], 200 * sizeof(short));

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
    // memcpy(new_2, mvs, sizeof(move_t) * (bp_dist + 1));
    memcpy_fast(new_2, mvs, sizeof(move_t) * (bp_dist + 1));
    return new_2;
}