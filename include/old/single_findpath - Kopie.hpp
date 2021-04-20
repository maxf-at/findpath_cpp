// this includes the basic elementary findpath algorithm

// #include <common.hpp>

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */



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

move_t* copy_moves(move_t* mvs, int bp_dist);
int compare_ptable(const void* A, const void* B);
int compare_energy(const void* A, const void* B);
int compare_moves_when(const void* A, const void* B);
void free_intermediate(intermediate_t* i);

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

    // auto free_intermediate(intermediate_t* i) -> void;
    // static auto compare_ptable(const void* A, const void* B) -> int;
    // static auto compare_energy(const void* A, const void* B) -> int;
    // static auto compare_moves_when(const void* A, const void* B) -> int;
    // auto        copy_moves(move_t* mvs, int bp_dist) -> move_t*;

   public:
    auto init(vrna_fold_compound_t* fc, short* pt1, short* pt2, int final_search_width)
        -> std::vector<sorted_path>;
    void test();
    single_findpath(){};
};

