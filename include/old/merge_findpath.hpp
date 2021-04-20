

class findpath
{
   private:
    int  test = 0;
    void init_pt();
    // std::vector<int_loops> find_interior_loop(int start, int end);
    s_graph process_ext_loops(ext_loops ext_loops, short* pt_1, short* pt_2);
    // void      process_ext_loops(ext_loops ext_loops, short* pt_1, short* pt_2);
    s_graph   process_int_loops(int_loops current_sections, short* pt_1, short* pt_2);
    int_loops all_sections{0, 0};

   public:
    float                 search_width_multiplier;
    char *                seq, *s1, *s2;
    short *               pt_1, *pt_2;
    vrna_fold_compound_t* fc = nullptr;
    // std::tuple<int, int>  bp_distance(short* pt_a, short* pt_b, int min_pos = 0,
    //                                   int max_pos = INT_MAX / 2);

    auto init_int() -> s_graph;
    auto init_ext() -> s_graph;

    // Constructors with fold compound or sequence
    findpath(vrna_fold_compound_t* init_fc, char* a_s1, char* a_s2, float sw);
    findpath(char* seq, char* s1, char* s2, float sw);
    findpath(vrna_fold_compound_t* a_fc, short* a_pt1, short* a_pt2, float sw);
    // s_graph(int size); //Constructor of the class
};

// Constructors
findpath::findpath(vrna_fold_compound_t* a_fc, char* a_s1, char* a_s2, float sw)
{

    this->fc = a_fc;
    // fc = a_fc;
    s1                      = a_s1;
    s2                      = a_s2;
    search_width_multiplier = sw;

    // IC(1, search_width_multiplier);

    init_pt();
}

findpath::findpath(vrna_fold_compound_t* a_fc, short* a_pt1, short* a_pt2, float sw)
{

    this->fc = a_fc;
    // fc = a_fc;
    pt_1                    = a_pt1;
    pt_2                    = a_pt2;
    search_width_multiplier = sw;

    s1 = vrna_db_from_ptable(pt_1);
    s2 = vrna_db_from_ptable(pt_2);
}

findpath::findpath(char* a_seq, char* a_s1, char* a_s2, float sw)
{
    seq                     = a_seq;
    s1                      = a_s1;
    s2                      = a_s2;
    search_width_multiplier = sw;

    // IC(1, search_width_multiplier);

    // set model params
    vrna_md_t md;
    set_model_details(&md);
    fc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);

    init_pt();
}

void findpath::init_pt()
{
    pt_1 = vrna_ptable(s1);
    pt_2 = vrna_ptable(s2);
}

auto findpath::init_int() -> s_graph
{

    all_sections.start = 1;
    all_sections.end   = pt_1[0];

    all_sections.bp_dist         = vrna_bp_distance(s1, s2);
    all_sections.nested_sections = find_interior_loops(
        pt_1, pt_2, 1, pt_1[0]);  // recursively add all nested sections (if available)

    // std::cout << "Sections:" << all_sections << "\n";
    auto G = process_int_loops(all_sections, pt_1, pt_2);
    
    return G;
}


auto findpath::init_ext() -> s_graph
{

    auto ext_loops = find_exterior_loops(pt_1, pt_2);
    auto G = process_ext_loops(ext_loops, pt_1, pt_2);

    // G.display_path(true);
    return G;
}




// quicksort compare functions (Ivos code)

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
auto print_moves(const auto& path, vrna_fold_compound_t* fc, const char* s1, bool show_path = true)
{
    auto const& moves = path.moves;

    short* pt;
    pt = vrna_ptable(s1);

    // std::cout << fmt::format("Hello!\n");

    float en     = vrna_eval_structure_pt(fc, pt) / 100.0;
    float max_en = float(-INT_MAX);

    if (show_path) fmt::print("{} {:7.2f} ({:4}/{:4})\n", s1, en, 0, 0);

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

        if (show_path) fmt::print("{} {:7.2f} ({:4}/{:4})\n", s, en, move.i, move.j);

        // printf("%.*s", std::abs(move.i), s);
        // fmt::print("{}", insert1);
        // // printf("%.*s", std::abs(move.i-move.j), s + std::abs(move.i-move.j));
        // printf("%.*s", std::abs(move.i-move.j-3), s + std::abs(move.i) + 1);
        // fmt::print("{}", insert2);
        // // printf("%.*s", std::abs(move.i-move.j-3), s + std::abs(move.j) + 1);
        // printf("%s", s + strlen(s) - std::abs(move.i));
    }

    // fmt::print("S: {:6.2f} kcal/mol\n", max_en);
    fmt::print("{:6.2f}\n", max_en);
}

auto available_edges(const auto& input_node, int G1_G2, auto& next_paths, const auto& current_path,
                     vrna_fold_compound_t* fc, int max_en, bool direction)
{
    // we use ingoing or outgoing edges for graph traversal, depending on direction

    if (direction) {
        for (const auto& current_edge : input_node.out_edges) {
            const auto& i = current_edge.i;
            const auto& j = current_edge.j;
            auto        current_en =
                current_path.current_en + vrna_eval_move_pt(fc, current_path.current_ptable, i, j);

            if (current_en <= max_en) {
                merge_path new_path;
                // are we currently processing the G1 or G2 node?
                if (G1_G2 == 1) {  // one move on G1
                    new_path.current_G1_node    = current_edge.destination;
                    new_path.current_G2_node    = current_path.current_G2_node;
                    new_path.current_G1_bp_dist = current_path.current_G1_bp_dist + 1;
                    new_path.current_G2_bp_dist = current_path.current_G2_bp_dist;

                } else {  // one move in graph G2
                    new_path.current_G1_node    = current_path.current_G1_node;
                    new_path.current_G2_node    = current_edge.destination;
                    new_path.current_G1_bp_dist = current_path.current_G1_bp_dist;
                    new_path.current_G2_bp_dist = current_path.current_G2_bp_dist + 1;
                }

                new_path.current_en = current_en;
                new_path.current_s  = std::max(current_en, current_path.current_s);
                // generate new pairing tables
                new_path.current_ptable = vrna_ptable_copy(current_path.current_ptable);
                new_path.moves          = current_path.moves;
                new_path.moves[new_path.current_G1_bp_dist + new_path.current_G2_bp_dist - 1] = {
                    i, j, current_en};
                // new_path.moves.push_back({i, j, current_en});

                if (j < 0) {  // delete a basepair
                    new_path.current_ptable[-i] = 0;
                    new_path.current_ptable[-j] = 0;
                } else {  // add a basepair
                    new_path.current_ptable[i] = j;
                    new_path.current_ptable[j] = i;
                }

                next_paths.push_back(new_path);
                // auto en = vrna_eval_structure_pt(fc, new_path.current_ptable);
                std::string s = vrna_db_from_ptable(new_path.current_ptable);
                // fmt::print("{} FW {} {} {}\n", s, i, j, current_en);
            }
        }
    } else {
        // traverse graph backwards - this has a few changes compared to the fwd version, lots of
        // duplication here...
        for (const auto& current_edge : input_node.in_edges) {
            const auto& i = current_edge.i;
            const auto& j = current_edge.j;

            // IC("r", i, j);

            auto current_en =
                current_path.current_en + vrna_eval_move_pt(fc, current_path.current_ptable, i, j);

            if (current_en <= max_en) {
                merge_path new_path;
                // are we currently processing the G1 or G2 node?
                if (G1_G2 == 1) {  // one move on G1
                    new_path.current_G1_node    = current_edge.destination;
                    new_path.current_G2_node    = current_path.current_G2_node;
                    new_path.current_G1_bp_dist = current_path.current_G1_bp_dist - 1;
                    new_path.current_G2_bp_dist = current_path.current_G2_bp_dist;

                } else {  // one move in graph G2
                    new_path.current_G1_node    = current_path.current_G1_node;
                    new_path.current_G2_node    = current_edge.destination;
                    new_path.current_G1_bp_dist = current_path.current_G1_bp_dist;
                    new_path.current_G2_bp_dist = current_path.current_G2_bp_dist - 1;
                }

                new_path.current_en = current_en;
                new_path.current_s  = std::max(current_en, current_path.current_s);
                // generate new pairing tables
                new_path.current_ptable = vrna_ptable_copy(current_path.current_ptable);
                new_path.moves          = current_path.moves;
                new_path.moves[new_path.current_G1_bp_dist + new_path.current_G2_bp_dist] = {
                    -i, -j, current_en};
                // new_path.moves.push_back({i, j, current_en});

                if (j < 0) {  // delete a basepair
                    new_path.current_ptable[-i] = 0;
                    new_path.current_ptable[-j] = 0;
                } else {  // add a basepair
                    new_path.current_ptable[i] = j;
                    new_path.current_ptable[j] = i;
                }

                next_paths.push_back(new_path);
                // auto en = vrna_eval_structure_pt(fc, new_path.current_ptable);
                std::string s = vrna_db_from_ptable(new_path.current_ptable);
                // if (G1_G2 == 1)
                // fmt::print("{} BW G1 {}->{} / {} {} {}\n", s, current_path.current_G1_node,
                // current_edge.destination, i, j, current_en); else fmt::print("{} BW G2 {}->{} /
                // {} {} {}\n", s, current_path.current_G2_node, current_edge.destination, i, j,
                // current_en);
            }
        }
        // IC("B");
    }
}

auto merge_once(auto G1, auto G2, short* pt_1, int s1_en, short* pt_2, int s2_en, bool direction,
                int total_bp_dist, int max_en, int merge_search_width, vrna_fold_compound_t* fc)
{
    // this is the start of findpath_once, assuming we have maxE etc...

    int current_i_node = 0;
    int current_j_node = 0;

    // direction = not direction;

    std::vector<merge_path> all_paths;

    // init the first path element, depending on direction
    merge_path current_path;
    current_path.current_G1_node    = 0;
    current_path.current_G2_node    = 0;
    current_path.current_G1_bp_dist = 0;
    current_path.current_G2_bp_dist = 0;
    current_path.current_ptable     = pt_1;
    current_path.current_en         = s1_en;
    current_path.current_s          = s1_en;

    if (not direction) {
        current_path.current_G1_node    = G1.bp_dist;
        current_path.current_G2_node    = G2.bp_dist;
        current_path.current_G1_bp_dist = G1.bp_dist;
        current_path.current_G2_bp_dist = G2.bp_dist;
        current_path.current_ptable     = pt_2;
        current_path.current_en         = s2_en;
        current_path.current_s          = s2_en;
    }

    current_path.moves.resize(G1.bp_dist + G2.bp_dist);  // fill move vector with zero tuples

    all_paths.push_back(current_path);

    const auto& G1_node = G1.node_list[current_path.current_G1_node];
    const auto& G2_node = G2.node_list[current_path.current_G2_node];

    // if (direction) { const auto& edges = G1_node.out_edges; }

    // G1 moves
    fmt::print("start merging, direction: {}, max_en: {} / distances: {} / {} \n", direction,
               max_en, G1.bp_dist, G2.bp_dist);
    // fmt::print("{} \n", vrna_db_from_ptable(all_paths[0].current_ptable));

    for (int d = 1; d <= total_bp_dist; d++) {
        // we will transfer all current paths (in all_paths) into next_paths for the next iteration
        std::vector<merge_path> next_paths;

        // reference?

        for (const auto current_path : all_paths) {
            // IC(current_path.current_G1_node, current_path.current_G2_node);
            // fmt::print("cp, direction: {}, {} \n", direction, d);
            // fmt::print("{} \n", vrna_db_from_ptable(current_path.current_ptable));

            // for (const auto& m : current_path.moves) { IC(m.i, m.j, m.E); }

            const auto& G1_node = G1.node_list[current_path.current_G1_node];
            const auto& G2_node = G2.node_list[current_path.current_G2_node];
            // fill up next_paths - either move 1 step on G1 or G2
            available_edges(G1_node, 1, next_paths, current_path, fc, max_en, direction);
            available_edges(G2_node, 2, next_paths, current_path, fc, max_en, direction);
        }
        // fmt::print("size before sort: {} \n", next_paths.size());

        // IC(d, next_paths.size(), 2);

        if (next_paths.size() == 0) {
            std::vector<merge_path> empty_paths{};

            merge_path element;
            element.current_s = max_en;

            empty_paths.push_back(element);

            return empty_paths;
        }

        // sort for unique pairing tables
        std::qsort(&next_paths[0], next_paths.size(), sizeof(merge_path), compare_ptable);

        // for (const auto& next_path : next_paths) {
        //     for (const auto& move : next_path.moves) { fmt::print("({} / {}) ", move.i, move.j);
        //     } fmt::print(" d={}, sorted: {} {}\n", d, next_path.current_en, next_path.current_s);
        // }

        // compare current with next pairing table...

        int len = next_paths[0].current_ptable[0];

        // remove duplicates (the struct has an equal operator, referencing ptables)
        // the best path is always at the first position, the rest gets pruned off
        next_paths.erase(std::unique(next_paths.begin(), next_paths.end()), next_paths.end());

        // if (d != total_bp_dist) {
        //     for (int u = 0, c = 1; c < next_paths.size(); c++) {
        //         if (memcmp(next_paths[u].current_ptable, next_paths[c].current_ptable,
        //                    sizeof(short) * len) != 0) {
        //             u++;
        //             // next_paths[u] = next_paths[c];
        //         } else {
        //             // next_paths[c].current_s = INT_MAX;
        //             // next_paths[c].current_en = INT_MAX;
        //             next_paths[c].current_s  = 999;
        //             next_paths[c].current_en = 999;
        //             // fmt::print("drop {}, {} \n", u, c);

        //             // free_intermediate(next + c);
        //         }
        //     }
        // }

        std::qsort(&next_paths[0], next_paths.size(), sizeof(merge_path), compare_energy);

        // IC(d, next_paths.size(), 3);

        // not needed?
        // if (next_paths.size() == 0) {
        //     std::vector<merge_path> empty_paths;
        //     return empty_paths;
        // }

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
    // IC("se");
    return all_paths;
}

auto merge_method(auto& G1, auto& G2, short* pt_1, int s1_en, short* pt_2, int s2_en,
                  int total_bp_dist, int max_en, int final_merge_search_width,
                  vrna_fold_compound_t* fc)
{
    // IC("method");

    // we might merge constant sections, ignore them
    if (G1.bp_dist == 0 ){
        // fmt::print("ignore merge \n");
        G2.pt_1 = pt_1;
        G2.pt_2 = pt_2;        
        return G2;    
    }
    if (G2.bp_dist == 0 ){
        // fmt::print("ignore merge \n");
        G1.pt_1 = pt_1;
        G1.pt_2 = pt_2;        
        return G1;    
    }

// or G2.bp_dist == 0


    if (final_merge_search_width < 17) { final_merge_search_width = 17; }
    int  current_merge_search_width = 16;
    bool direction                  = true;

    int              last_iteration = final_merge_search_width;
    std::vector<int> iterations{final_merge_search_width};

    while (last_iteration > 16) {
        const int next_iteration = int(last_iteration / 2.0);
        iterations.push_back(next_iteration);
        last_iteration = next_iteration;
    }

    // IC(iterations);
    // iterations = {1};

    std::vector<merge_path> all_paths{};

    // for (const auto& current_merge_search_width : iterations | std::views::reverse) {
    for (std::vector<int>::reverse_iterator it = iterations.rbegin(); it != iterations.rend();
         ++it) {
        current_merge_search_width = *it;

        direction = true;
        std::vector<merge_path> fwd_paths =
            merge_once(G1, G2, pt_1, s1_en, pt_2, s2_en, direction, total_bp_dist, max_en,
                       current_merge_search_width, fc);

        direction = false;
        std::vector<merge_path> bwd_paths =
            merge_once(G1, G2, pt_1, s1_en, pt_2, s2_en, direction, total_bp_dist, max_en,
                       current_merge_search_width, fc);

        // IC(current_merge_search_width, fwd_paths[0].current_s, bwd_paths[0].current_s);

        // current_merge_search_width = 100;
        max_en = std::max(fwd_paths[0].current_s, bwd_paths[0].current_s);

        // concatenate result vectors
        std::move(fwd_paths.begin(), fwd_paths.end(), std::back_inserter(all_paths));
        std::move(bwd_paths.begin(), bwd_paths.end(), std::back_inserter(all_paths));
    }

    // move best path (lowest saddle energy) to [0]
    std::sort(all_paths.begin(), all_paths.end(),
              [](const auto& a, const auto& b) -> bool { return a.current_s < b.current_s; });

    s_graph G_merged{fc, pt_1, pt_2, total_bp_dist, all_paths};
    G_merged.max_en = all_paths[0].current_s;

    return G_merged;
    // return all_paths;
}

s_graph findpath::process_ext_loops(ext_loops ext_loops, short* pt_1, short* pt_2)
// void findpath::process_ext_loops(ext_loops ext_loops, short* pt_1, short* pt_2)
{
    // return;

    bool Verbose = false;
    // bool Verbose                       = true;
    int merge_search_width_multiplier = 1;


    bool    first_run = true;
    s_graph last_G;
    short*  last_pt_1;
    short*  last_pt_2;
    int     last_bp_dist = -1;

    // s_graph temp;

    for (int i = 0; i < ext_loops.size(); i++) {
        // if (i==1){
        //     continue;
        // }

        const auto& current_ext_loop = ext_loops[i];
        // const auto& next_ext_loop = ext_loops[i+1];

        // std::cout << current_ext_loop << " / "
        //           << "\n";
        // fmt::print("process ext loops {}\n", i);

        short* current_pt_1 = vrna_ptable_copy(pt_1);
        short* current_pt_2 = vrna_ptable_copy(pt_2);

        // adjust pairing table
        for (int i = 1; i <= pt_1[0]; i++) {
            if (i < current_ext_loop.start or i > current_ext_loop.end) {
                current_pt_1[i] = 0;
                current_pt_2[i] = 0;
            }
        }
        int current_bp_dist = bp_distance(current_pt_1, current_pt_2);

        // if (last_bp_dist == 0 or current_bp_dist == 0) {
        //     fmt::print("ignore: {} / {}\n", last_bp_dist, current_bp_dist);
        // }

        // ext loops without further interior loop recursion
        // const int search_width = current_bp_dist * search_width_multiplier;
        // single_findpath test;
        // auto            result = test.init(fc, current_pt_1, current_pt_2, search_width, true);
        // s_graph current_G {fc, current_pt_1, current_pt_2, current_bp_dist, result};
        // current_G.max_en = result[0].max_en;

        s_graph current_G = process_int_loops(current_ext_loop, current_pt_1, current_pt_2);
        current_G.bp_dist = current_bp_dist;

        // fmt::print("current bp dist: {}\n", current_bp_dist);

        if (first_run) {
            // temp = current_G;

            last_G       = current_G;
            last_pt_1    = current_pt_1;
            last_pt_2    = current_pt_2;
            first_run    = false;
            last_bp_dist = current_bp_dist;
            continue;
        }

        // last_G    = current_G;
        // last_pt_1 = current_pt_1;
        // last_pt_2 = current_pt_2;
        // first_run = false;
        // continue;

        int merge_bp_dist      = last_G.bp_dist + current_G.bp_dist;
        int max_en             = INT_MAX;
        int merge_search_width = merge_bp_dist * merge_search_width_multiplier;

        // merge pairing tables (outer and current inner)
        short* merged_pt_1 = vrna_ptable_copy(current_pt_1);
        short* merged_pt_2 = vrna_ptable_copy(current_pt_2);

        for (int i = 1; i <= merged_pt_1[0]; i++) {
            // fmt::print("i {}\n", i);
            if ((merged_pt_1[i] == 0 and last_pt_1 != 0) or
                (merged_pt_1[i] != 0 and last_pt_1 == 0)) {
                merged_pt_1[i] = last_pt_1[i];
            }
            if ((merged_pt_2[i] == 0 and last_pt_2 != 0) or
                (merged_pt_2[i] != 0 and last_pt_2 == 0)) {
                merged_pt_2[i] = last_pt_2[i];
            }
        }

        // current_G.display_path();
        // current_G.info();
        // last_G.display_path();
        // last_G.info();

        int s1_en = vrna_eval_structure_pt(fc, merged_pt_1);
        int s2_en = vrna_eval_structure_pt(fc, merged_pt_2);

        auto start = std::chrono::high_resolution_clock::now();

        auto G_merged = merge_method(current_G, last_G, merged_pt_1, s1_en, merged_pt_2, s2_en,
                                      merge_bp_dist, max_en, merge_search_width, fc);

        auto                          finish  = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        // IC(elapsed.count());

        // fmt::print("merge time: {} \n", elapsed.count());

        // outer_bp_dist += G_inner.bp_dist;
        last_pt_1    = merged_pt_1;
        last_pt_2    = merged_pt_2;
        last_bp_dist = current_bp_dist;

        // s_graph G_merged{fc, merged_pt_1, merged_pt_2, merge_bp_dist, all_paths};
        // G_merged.max_en = all_paths[0].current_s;
        last_G          = G_merged;

        // fmt::print("display merged path:\n");
        // last_G.display_path();
        // fmt::print("merged: {} {}\n", all_paths.size(), all_paths[0].current_s);
        // G_outer.info();
    }

    return last_G;
}

s_graph findpath::process_int_loops(int_loops current_sections, short* pt_1, short* pt_2)
{
    bool Verbose = false;
    // bool Verbose = true;

    // int  search_width_multiplier = 1;
    // float search_width_multiplier       = 2;
    int merge_search_width_multiplier = 1;

    if (current_sections.nested_sections.size() == 0) {
        // top of the recursion tree - findpath for a single section with start and end

        const int search_width = current_sections.bp_dist * search_width_multiplier;
        // const int search_width = 200;

        single_findpath test;
        auto            result = test.init(fc, pt_1, pt_2, search_width, true);

        s_graph G_inner{fc, pt_1, pt_2, current_sections.bp_dist, result};
        G_inner.max_en = result[0].max_en;

        // s_graph G_inner;
        // G_inner.fc      = fc;
        // G_inner.pt_1    = pt_1;
        // G_inner.pt_2    = pt_2;
        // G_inner.bp_dist = current_sections.bp_dist;
        // G_inner.add_paths(result);

        // G_inner.info();

        if (Verbose) print_moves(result[0], fc, vrna_db_from_ptable(pt_1), true);

        return G_inner;
    }

    // else: generate the outer pairing table and substract all inner sections, where we will
    // recurse into
    // IC("s", current_sections);
    // IC(current_sections.start, current_sections.end, current_sections.bp_dist);

    int outer_bp_dist = current_sections.bp_dist;

    short* outer_pt_1 = vrna_ptable_copy(pt_1);
    short* outer_pt_2 = vrna_ptable_copy(pt_2);

    for (int i = 1; i <= pt_1[0]; i++) {
        if (i < current_sections.start or i > current_sections.end) {
            outer_pt_1[i] = 0;
            outer_pt_2[i] = 0;
        }
    }

    // iterate over nested inner sections to adjust outer section
    for (const auto& nested_section : current_sections.nested_sections) {
        outer_bp_dist -= nested_section.bp_dist;

        // IC("r", nested_section.start, nested_section.end);

        for (int i = 1; i <= pt_1[0]; i++) {
            if (i > nested_section.start and i < nested_section.end) {
                outer_pt_1[i] = 0;
                outer_pt_2[i] = 0;
            }
        }
    }

    const int search_width = outer_bp_dist * search_width_multiplier;
    // const int search_width = 200;

    single_findpath test;
    auto            result = test.init(fc, outer_pt_1, outer_pt_2, search_width, true);

    s_graph G_outer{fc, outer_pt_1, outer_pt_2, outer_bp_dist, result};

    // G_outer.info();

    if (Verbose) print_moves(result[0], fc, vrna_db_from_ptable(outer_pt_1), true);

    // fmt::print("outer section 1:\n{}\n", vrna_db_from_ptable(outer_pt_1));
    // fmt::print("outer section 2:\n{}\n", vrna_db_from_ptable(outer_pt_2));

    // iterate over nested inner sections to adjust outer section
    for (const auto& nested_section : current_sections.nested_sections) {
        short* inner_pt_1 = vrna_ptable_copy(pt_1);
        short* inner_pt_2 = vrna_ptable_copy(pt_2);
        for (int i = 1; i <= pt_1[0]; i++) {
            if (i < nested_section.start or i > nested_section.end) {
                inner_pt_1[i] = 0;
                inner_pt_2[i] = 0;
            }
        }
        // fmt::print("inner section 1:\n{}\n", vrna_db_from_ptable(inner_pt_1));
        // fmt::print("inner section 2:\n{}\n", vrna_db_from_ptable(inner_pt_2));

        // recursive call to build up the current inner graph
        s_graph G_inner = process_int_loops(nested_section, inner_pt_1, inner_pt_2);

        // merge pairing tables (outer and current inner)
        short* merged_pt_1 = vrna_ptable_copy(outer_pt_1);
        short* merged_pt_2 = vrna_ptable_copy(outer_pt_2);

        for (int i = 1; i < merged_pt_1[0]; i++) {
            // fmt::print("i {}\n", i);
            if ((merged_pt_1[i] == 0 and G_inner.pt_1 != 0) or
                (merged_pt_1[i] != 0 and G_inner.pt_1 == 0)) {
                merged_pt_1[i] = G_inner.pt_1[i];
            }
            if ((merged_pt_2[i] == 0 and G_inner.pt_2 != 0) or
                (merged_pt_2[i] != 0 and G_inner.pt_2 == 0)) {
                merged_pt_2[i] = G_inner.pt_2[i];
            }
        }

        int s1_en = vrna_eval_structure_pt(fc, merged_pt_1);
        int s2_en = vrna_eval_structure_pt(fc, merged_pt_2);

        // fmt::print("starting merge:\n{} {} {}\n", vrna_db_from_ptable(merged_pt_1),
        // G_outer.bp_dist,
        //    s1_en);
        // fmt::print("{} {} {}\n", vrna_db_from_ptable(merged_pt_2), G_inner.bp_dist, s2_en);

        int total_bp_dist = G_inner.bp_dist + G_outer.bp_dist;

        // start here the merge supervisor gradually increasing search width

        int max_en = INT_MAX;
        // int  max_en             = -1390;

        int merge_search_width = total_bp_dist * merge_search_width_multiplier;

        auto start = std::chrono::high_resolution_clock::now();

        auto G_merged = merge_method(G_outer, G_inner, merged_pt_1, s1_en, merged_pt_2, s2_en,
                                      total_bp_dist, max_en, merge_search_width, fc);

        auto                          finish  = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        // IC(elapsed.count());

        fmt::print("merged: {} \n", G_merged.max_en);
        // IC("merged");

        outer_bp_dist += G_inner.bp_dist;
        outer_pt_1 = merged_pt_1;
        outer_pt_2 = merged_pt_2;

        // s_graph G_new{fc, merged_pt_1, merged_pt_2, outer_bp_dist, all_paths};
        // G_new.max_en = all_paths[0].current_s;
        G_outer      = G_merged;

        // IC("merged2");

        // if (Verbose) print_moves(all_paths[0], fc, vrna_db_from_ptable(merged_pt_1), true);

        // print_moves(all_paths[0], fc, vrna_db_from_ptable(merged_pt_1), false);

        // break;
    }

    // IC(current_sections.nested_sections[0].start, current_sections.nested_sections[0].end,
    //    current_sections.nested_sections[0].bp_dist);
    // IC(outer_bp_dist);
    if (Verbose) { G_outer.display_path(); }

    return G_outer;
}
