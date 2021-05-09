

class findpath
{
   private:
    int  test = 0;
    void init_pt();

    // fp cache

    auto single_findpath_cache(vrna_fold_compound_t* fc, short* pt1, short* pt2, int search_width,
                               bool mp);

    auto      process_ext_loops(ext_loops ext_loops, short* pt_1, short* pt_2) -> s_graph;
    auto      process_int_loops(int_loops current_sections, short* pt_1, short* pt_2) -> s_graph;
    int_loops all_sections{0, 0};

    std::map<size_t, s_graph> G_cache;
    bool                                cache = true;

   public:
    float search_width_multiplier;
    bool  mp;

    // todo: delete these...
    // char *                seq, *s1, *s2;
    // short *               pt_1, *pt_2;
    vrna_fold_compound_t* fc;

    auto init(std::string s1, std::string s2, float sw = 2) -> s_graph;
    auto init_python(std::string s1, std::string s2, float sw = 2) -> int;

    auto init_ext(std::string s1, std::string s2, float sw = 2) -> s_graph;

    // various constructors (with fold compound or sequence)

    findpath(std::string sequence, bool mp = true);

    // findpath(vrna_fold_compound_t* init_fc, char* a_s1, char* a_s2, float sw);
    // findpath(char* seq, char* s1, char* s2, float sw);
    // findpath(vrna_fold_compound_t* a_fc, short* a_pt1, short* a_pt2, float sw);
};

findpath::findpath(std::string sequence, bool mp) : mp{mp}
{
    vrna_md_t md;
    set_model_details(&md);
    fc = vrna_fold_compound(sequence.c_str(), &md, VRNA_OPTION_EVAL_ONLY);

    // s1          = vrna_db_from_ptable(pt_1);
    // s2          = vrna_db_from_ptable(pt_2);
}

auto findpath::init_python(std::string s1, std::string s2, float sw) -> int
{
    // s1 = vrna_db_from_ptable(pt_1);
    // s2 = vrna_db_from_ptable(pt_2);

    search_width_multiplier = sw;
    short* pt1              = vrna_ptable(s1.c_str());
    short* pt2              = vrna_ptable(s2.c_str());

    // initialise section structure
    all_sections.start = 1;
    all_sections.end   = pt1[0];

    all_sections.bp_dist         = vrna_bp_distance(s1.c_str(), s2.c_str());
    all_sections.nested_sections = find_interior_loops(
        pt1, pt2, 1, pt1[0]);  // recursively add all nested sections (if available)

    std::cout << "Sections:" << all_sections << "\n";
    // fmt::print ("s1: {}\n", vrna_db_from_ptable(pt1));
    // fmt::print ("s2: {}\n", vrna_db_from_ptable(pt2));

    // recursively process nested sections
    auto G = process_int_loops(all_sections, pt1, pt2);
    return G.max_en;
}

auto findpath::init(std::string s1, std::string s2, float sw) -> s_graph
{
    // s1 = vrna_db_from_ptable(pt_1);
    // s2 = vrna_db_from_ptable(pt_2);

    search_width_multiplier = sw;
    short* pt1              = vrna_ptable(s1.c_str());
    short* pt2              = vrna_ptable(s2.c_str());

    // initialise section structure
    all_sections.start = 1;
    all_sections.end   = pt1[0];

    all_sections.bp_dist         = vrna_bp_distance(s1.c_str(), s2.c_str());
    all_sections.nested_sections = find_interior_loops(
        pt1, pt2, 1, pt1[0]);  // recursively add all nested sections (if available)

    // std::cout << "Sections:" << all_sections << "\n";
    // fmt::print ("s1: {}\n", vrna_db_from_ptable(pt1));
    // fmt::print ("s2: {}\n", vrna_db_from_ptable(pt2));

    // recursively process nested sections
    auto G = process_int_loops(all_sections, pt1, pt2);
    return G;
}

auto findpath::init_ext(std::string s1, std::string s2, float sw) -> s_graph
{
    search_width_multiplier = sw;
    short* pt1              = vrna_ptable(s1.c_str());
    short* pt2              = vrna_ptable(s2.c_str());

    // variant with exterior loops
    auto ext_loops = find_exterior_loops(pt1, pt2);
    auto G         = process_ext_loops(ext_loops, pt1, pt2);

    // G.display_path(true);
    return G;
}

// quicksort compare functions (mostly origianl findpath code)

auto compare_ptable(const void* A, const void* B) -> int
{
    merge_path *a, *b;
    int         c;

    a = (merge_path*)A;
    b = (merge_path*)B;

    // c = memcmp(a->current_ptable, b->current_ptable, a->current_ptable[0] * sizeof(short));
    // if (c != 0) return c;

    // if both structures are not identical, sort them according to the hash value
    if (a->s_hash != b->s_hash) {
        if (a->s_hash > b->s_hash) {
            return 1;
        } else {
            return -1;
        }
    }

    if ((a->current_s - b->current_s) != 0) return a->current_s - b->current_s;

    return a->current_en - b->current_en;
}

auto compare_energy(const void* A, const void* B) -> int
{
    merge_path *a, *b;

    a = (merge_path*)A;
    b = (merge_path*)B;

    if ((a->current_s - b->current_s) != 0) return a->current_s - b->current_s;

    return a->current_en - b->current_en;
}

auto print_moves(const auto& path, vrna_fold_compound_t* fc, const char* s1, bool show_path = true)
{
    // this print path function is unused (?)
    auto const& moves = path.moves;

    short* pt;
    pt = vrna_ptable(s1);

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

        const char* s = vrna_db_from_ptable(pt);
        en            = vrna_eval_structure_pt(fc, pt) / 100.0;
        if (en > max_en) { max_en = en; }

        if (show_path) fmt::print("{} {:7.2f} ({:4}/{:4})\n", s, en, move.i, move.j);
    }

    // fmt::print("S: {:6.2f} kcal/mol\n", max_en);
    fmt::print("{:6.2f}\n", max_en);
}

auto available_edges(const auto& input_node, int G1_G2, auto& next_paths, const auto& current_path,
                     int path_index, vrna_fold_compound_t* fc, int max_en, bool direction,
                     int pool_offset, int total_bp_dist,
                     std::vector<std::vector<sorted_move>>& move_storage,
                     std::vector<std::vector<short>>&       ptable_storage, const s_graph &G1, const s_graph &G2) -> void
{
    // this is the try_moves equivalent (regular findpath)
    // given:  an input node (either G1 or G2) + direction
    // output: possible moves for the next move
    //         everything gets appended to next_paths

    // we use ingoing or outgoing edges for graph traversal, depending on direction
    if (direction) {
        for (const auto& current_edge : input_node.out_edges) {
            int i = current_edge.i;
            int j = current_edge.j;

            // auto        current_en =
            //     current_path.current_en + vrna_eval_move_pt(fc, current_path.current_ptable, i,
            //     j);
            auto current_en = current_path.current_en + current_edge.en;

            if (current_en <= max_en) {
                auto& new_path = next_paths.emplace_back();

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

                const int n = next_paths.size();
                // new_path.move_ptr = move_storage[n + pool_offset].data();
                // memcpy(new_path.move_ptr, current_path.move_ptr,
                //        sizeof(sorted_move) * (total_bp_dist));
                // new_path.move_ptr[new_path.current_G1_bp_dist + new_path.current_G2_bp_dist - 1]
                // = {
                //     i, j, current_en};

                // generate new pairing tables
                // new_path.current_ptable = vrna_ptable_copy(current_path.current_ptable);

                // const int len           = current_path.current_ptable[0];
                // new_path.current_ptable = ptable_storage[n + pool_offset].data();
                // memcpy(new_path.current_ptable, current_path.current_ptable,
                //        sizeof(short) * (len + 1));

                // if (j < 0) {  // delete a basepair
                //     new_path.current_ptable[-i] = 0;
                //     new_path.current_ptable[-j] = 0;
                // } else {  // add a basepair
                //     new_path.current_ptable[i] = j;
                //     new_path.current_ptable[j] = i;
                // }

                new_path.s_hash = G1.node_list[new_path.current_G1_node].s_hash +
                    G2.node_list[new_path.current_G2_node].s_hash;

                new_path.i_move     = i;
                new_path.j_move     = j;
                new_path.last_index = path_index;
                // new_path.last_bp = new_path.current_G1_bp_dist + new_path.current_G2_bp_dist - 1;
                // fmt::print("{} FW {} {} {}\n", s, i, j, current_en);
            }
        }
    } else {
        // traverse graph backwards - this has a few changes compared to the fwd version, lots of
        // duplication here...
        for (const auto& current_edge : input_node.in_edges) {
            int i = current_edge.i;
            int j = current_edge.j;

            // auto current_en =
            // current_path.current_en + vrna_eval_move_pt(fc, current_path.current_ptable, i, j);
            auto current_en = current_path.current_en - current_edge.en;

            if (current_en <= max_en) {
                auto& new_path = next_paths.emplace_back();

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

                int n = next_paths.size();
                // new_path.move_ptr = move_storage[n + pool_offset].data();
                // memcpy(new_path.move_ptr, current_path.move_ptr,
                //        sizeof(sorted_move) * (total_bp_dist));
                // new_path.move_ptr[new_path.current_G1_bp_dist + new_path.current_G2_bp_dist] = {
                //     -i, -j, current_en};
                // fmt::print("n: {}, offset: {}, cptr: {}\n", n, pool_offset, current_path.move_ptr
                // == nullptr);

                // generate new pairing tables
                // new_path.current_ptable = vrna_ptable_copy(current_path.current_ptable);
                // const int len           = current_path.current_ptable[0];
                // new_path.current_ptable = ptable_storage[n + pool_offset].data();
                // memcpy(new_path.current_ptable, current_path.current_ptable,
                //        sizeof(short) * (len + 1));

                // if (j < 0) {  // delete a basepair
                //     new_path.current_ptable[-i] = 0;
                //     new_path.current_ptable[-j] = 0;
                // } else {  // add a basepair
                //     new_path.current_ptable[i] = j;
                //     new_path.current_ptable[j] = i;
                // }

                new_path.s_hash = G1.node_list[new_path.current_G1_node].s_hash +
                    G2.node_list[new_path.current_G2_node].s_hash;

                new_path.i_move     = -i;
                new_path.j_move     = -j;
                new_path.last_index = path_index;
            }
        }
    }
}

auto merge_once(const s_graph G1, const s_graph G2, short* pt, int en_start, bool direction,
                int total_bp_dist, int max_en, int merge_search_width, vrna_fold_compound_t* fc)
    -> std::vector<merge_path>
{
    // merge 2 given graphs with a given search width and a given direction

    // fmt::print("paths: {} {} / bp dist: {}\n", G1.paths_count, G2.paths_count, total_bp_dist);

    // std::vector<std::vector<sorted_move>> move_storage(
    //     2 * (G1.paths_count + G2.paths_count) * merge_search_width + 2);
    // for (auto& col : move_storage) { col.reserve(total_bp_dist); }
    // std::vector<std::vector<short>> ptable_storage(
    //     2 * (G1.paths_count + G2.paths_count) * merge_search_width + 2);
    // for (auto& col : ptable_storage) { col.reserve(pt[0] + 1); }

    std::vector<std::vector<sorted_move>> move_storage(
        2 * merge_search_width + 2);
    for (auto& col : move_storage) { col.reserve(total_bp_dist); }
    std::vector<std::vector<short>> ptable_storage(
        2 * merge_search_width + 2);
    for (auto& col : ptable_storage) { col.reserve(pt[0] + 1); }


    // std::vector<std::vector<sorted_move>> move_storage;

    std::vector<sorted_move> init_move(G1.bp_dist + G2.bp_dist);

    std::vector<merge_path> all_paths;

    auto& current_path = all_paths.emplace_back();

    // init first path element, depending on direction
    if (direction) {
        current_path.current_G1_node    = 0;
        current_path.current_G2_node    = 0;
        current_path.current_G1_bp_dist = 0;
        current_path.current_G2_bp_dist = 0;
        current_path.current_ptable     = pt;
        current_path.current_en         = en_start;
        current_path.current_s          = en_start;
    } else {
        current_path.current_G1_node    = G1.bp_dist;
        current_path.current_G2_node    = G2.bp_dist;
        current_path.current_G1_bp_dist = G1.bp_dist;
        current_path.current_G2_bp_dist = G2.bp_dist;
        current_path.current_ptable     = pt;
        current_path.current_en         = en_start;
        current_path.current_s          = en_start;
    }

    current_path.move_ptr = init_move.data();

    const auto& G1_node = G1.node_list[current_path.current_G1_node];
    const auto& G2_node = G2.node_list[current_path.current_G2_node];

    // until we find s2
    for (int d = 1; d <= total_bp_dist; d++) {
        // iterate over all current paths, the next iteration will be stored in next_paths
        std::vector<merge_path> next_paths;

        // memory pool double buffer offset
        int pool_offset = 0;
        if (d % 2 == 0) {
            // pool_offset = (G1.paths_count + G2.paths_count) * merge_search_width + 2;
            pool_offset = merge_search_width + 2;
        }

        // for (const auto &current_path : all_paths) {
        for (int i = 0; i < all_paths.size(); i++) {
            const auto& current_path = all_paths[i];
            // fmt::print("cp, direction: {}, {} \n", direction, d);
            // fmt::print("{} \n", vrna_db_from_ptable(current_path.current_ptable));

            const auto& G1_node = G1.node_list[current_path.current_G1_node];
            const auto& G2_node = G2.node_list[current_path.current_G2_node];

            // fill up next_paths - either move 1 step on G1 or G2
            available_edges(G1_node, 1, next_paths, current_path, i, fc, max_en, direction,
                            pool_offset, total_bp_dist, move_storage, ptable_storage, G1, G2);
            available_edges(G2_node, 2, next_paths, current_path, i, fc, max_en, direction,
                            pool_offset, total_bp_dist, move_storage, ptable_storage, G1, G2);
        }

        // fmt::print("size before sort: {} \n", next_paths.size());

        if (next_paths.size() == 0) {
            // nothing found below max_en

            // merge_path element;
            // element.current_s = max_en;
            // all_paths = {element};
            // break;

            std::vector<merge_path> empty_paths{};
            merge_path              element;
            element.current_s = max_en;
            empty_paths.push_back(element);
            return empty_paths;
        }

        // sort for unique pairing tables
        std::qsort(next_paths.data(), next_paths.size(), sizeof(merge_path), compare_ptable);

        // remove duplicates (the struct has an equal operator, referencing ptables)
        // the best path is always at the first position, the rest gets pruned off
        next_paths.erase(std::unique(next_paths.begin(), next_paths.end()), next_paths.end());

        // sort best energies to [0]
        std::qsort(next_paths.data(), next_paths.size(), sizeof(merge_path), compare_energy);

        // resize vector if needed (this drops all many entries), we're done
        if (next_paths.size() > merge_search_width) { next_paths.resize(merge_search_width); }

        // update moves & pairing tables for the next iteration
        for (int i = 0; i < next_paths.size(); i++) {
            auto&     new_path   = next_paths[i];
            const int last_index = new_path.last_index;
            const int i_move     = new_path.i_move;
            const int j_move     = new_path.j_move;

            const auto& current_path = all_paths[last_index];

            // fmt::print("d: {} n: {}, offset: {}, cptr: {} / i {} j {} / last n: {}\n", d, i,
            //            pool_offset, current_path.move_ptr == nullptr, i_move, j_move, last_index);

            // update moves
            new_path.move_ptr = move_storage[i + pool_offset].data();
            memcpy(new_path.move_ptr, current_path.move_ptr, sizeof(sorted_move) * (total_bp_dist));

            if (direction) {
                new_path.move_ptr[new_path.current_G1_bp_dist + new_path.current_G2_bp_dist - 1] = {
                    i_move, j_move, new_path.current_en};
            } else {
                new_path.move_ptr[new_path.current_G1_bp_dist + new_path.current_G2_bp_dist] = {
                    i_move, j_move, new_path.current_en};
            }

            // fmt::print("d: {} n: {}, offset: {}, cptr: {} / i {} j {} / last n: {}\n", d, i,
            //            pool_offset, current_path.current_ptable == nullptr, i_move, j_move, last_index);

            // update pairing tables
            const int len           = current_path.current_ptable[0];
            new_path.current_ptable = ptable_storage[i + pool_offset].data();
            memcpy(new_path.current_ptable, current_path.current_ptable,
                    sizeof(short) * (len + 1));

            if (j_move < 0) {  // delete a basepair
                new_path.current_ptable[-i_move] = 0;
                new_path.current_ptable[-j_move] = 0;
            } else {  // add a basepair
                new_path.current_ptable[i_move] = j_move;
                new_path.current_ptable[j_move] = i_move;
            }

        }

        all_paths = next_paths;
    }

    // copy moves from the temporary move storage pool into the vector
    // first, init a new vector and allocate space
    all_paths[0].moves = {};
    all_paths[0].moves.resize(G1.bp_dist + G2.bp_dist);

    if (all_paths[0].move_ptr != nullptr) {
        sorted_move* dest = all_paths[0].moves.data();
        sorted_move* src  = all_paths[0].move_ptr;
        memcpy(dest, src, sizeof(sorted_move) * (total_bp_dist));
    }

    all_paths[0].current_ptable = nullptr;

    return all_paths;
}

auto merge_method(auto& G1, auto& G2, short* pt_1, int s1_en, short* pt_2, int s2_en,
                  int total_bp_dist, int max_en, int final_merge_search_width,
                  vrna_fold_compound_t* fc, bool merge_constant_sections)
{
    // merge 2 graphs, increase search widths, consider both directions, etc.

    // workaround for ext_loops version
    if (merge_constant_sections) {
        // we might merge constant sections, ignore them
        if (G1.bp_dist == 0) {
            // fmt::print("ignore merge \n");
            G2.pt_1 = pt_1;
            G2.pt_2 = pt_2;
            return G2;
        }
        if (G2.bp_dist == 0) {
            // fmt::print("ignore merge \n");
            G1.pt_1 = pt_1;
            G1.pt_2 = pt_2;
            return G1;
        }
    }

    // fix this...
    total_bp_dist = G1.bp_dist + G2.bp_dist;

    // generate search widths vector [400, 84, 16]

    if (final_merge_search_width < 17) { final_merge_search_width = 17; }
    int  current_merge_search_width = 16;
    bool direction                  = true;

    int             last_iteration = final_merge_search_width;
    std::deque<int> iterations{final_merge_search_width};

    // todo: optimize this
    while (last_iteration > 16) {
        const int next_iteration = int(last_iteration / 2.0);
        iterations.push_front(next_iteration);
        last_iteration = next_iteration;
    }

    // all paths will be stored here...
    std::vector<merge_path> all_paths{};

    // iterate search width vector from back to front
    for (const int current_merge_search_width : iterations) {
        direction                         = true;
        std::vector<merge_path> fwd_paths = merge_once(
            G1, G2, pt_1, s1_en, direction, total_bp_dist, max_en, current_merge_search_width, fc);
        direction                         = false;
        std::vector<merge_path> bwd_paths = merge_once(
            G1, G2, pt_2, s2_en, direction, total_bp_dist, max_en, current_merge_search_width, fc);

        // std::future<std::vector<merge_path>> ret1 =
        //     std::async(std::launch::async, &merge_once, G1, G2, pt_1, s1_en, true, total_bp_dist,
        //                max_en, current_merge_search_width, fc);
        // auto                                 G1_copy = G1;
        // auto                                 G2_copy = G2;
        // std::future<std::vector<merge_path>> ret2 =
        //     std::async(std::launch::async, &merge_once, G1_copy, G2_copy, pt_2, s2_en, false,
        //                total_bp_dist, max_en, current_merge_search_width, fc);
        // std::vector<merge_path> fwd_paths = ret1.get();
        // std::vector<merge_path> bwd_paths = ret2.get();

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

    // postprocessing: transform paths into graph without redundancy
    max_en = all_paths[0].current_s;

    s_graph G_merged{fc, pt_1, pt_2, total_bp_dist, all_paths, max_en};

    return G_merged;
}

auto findpath::process_ext_loops(ext_loops ext_loops, short* pt_1, short* pt_2) -> s_graph
{
    // this needs some work...

    bool Verbose = false;
    // bool Verbose                       = true;
    int merge_search_width_multiplier = 1;

    bool    first_run = true;
    s_graph last_G;
    short*  last_pt_1;
    short*  last_pt_2;
    int     last_bp_dist = -1;

    // iterate over available ext loops

    for (int i = 0; i < ext_loops.size(); i++) {
        const auto& current_ext_loop = ext_loops[i];

        std::cout << "current ext loop:" << current_ext_loop << " / "
                  << "\n";
        fmt::print("process ext loops {}\n", i);

        short* current_pt_1 = vrna_ptable_copy(pt_1);
        short* current_pt_2 = vrna_ptable_copy(pt_2);

        // fmt::print("s1 {}\n", vrna_db_from_ptable(pt_1));
        // fmt::print("s2 {}\n", vrna_db_from_ptable(pt_2));

        // adjust pairing table
        for (int i = 1; i <= pt_1[0]; i++) {
            if (i < current_ext_loop.start or i > current_ext_loop.end) {
                current_pt_1[i] = 0;
                current_pt_2[i] = 0;
            }
        }
        int current_bp_dist = bp_distance(current_pt_1, current_pt_2);
        // current_G.bp_dist = current_bp_dist;

        // fmt::print("s1 {}\n", vrna_db_from_ptable(current_pt_1));
        // fmt::print("s2 {}\n", vrna_db_from_ptable(current_pt_2));

        // ext loops without further interior loop recursion
        // const int search_width = current_bp_dist * search_width_multiplier;
        // single_findpath test;
        // auto            result = test.init(fc, current_pt_1, current_pt_2, search_width, true);
        // s_graph current_G {fc, current_pt_1, current_pt_2, current_bp_dist, result};
        // current_G.max_en = result[0].max_en;

        s_graph current_G = process_int_loops(current_ext_loop, current_pt_1, current_pt_2);
        current_G.bp_dist = current_bp_dist;

        // fmt::print("ext: current bp dist: {}\n", current_bp_dist);

        if (first_run) {
            last_G       = current_G;
            last_pt_1    = current_pt_1;
            last_pt_2    = current_pt_2;
            first_run    = false;
            last_bp_dist = current_bp_dist;
            continue;
        }

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
                                     merge_bp_dist, max_en, merge_search_width, fc, false);

        auto                          finish  = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        // IC(elapsed.count());

        // fmt::print("s1 {}\n", vrna_db_from_ptable(merged_pt_1));
        // fmt::print("s2 {}\n", vrna_db_from_ptable(merged_pt_2));
        // fmt::print("s1 {}\n", vrna_db_from_ptable(G_merged.pt_1));
        // fmt::print("s2 {}\n", vrna_db_from_ptable(G_merged.pt_2));
        // fmt::print("merge time: {} \n", elapsed.count());

        last_pt_1    = merged_pt_1;
        last_pt_2    = merged_pt_2;
        last_bp_dist = current_bp_dist;
        last_G       = G_merged;
        // fmt::print("display merged path:\n");
        // last_G.display_path();
        // fmt::print("merged: {} {}\n", all_paths.size(), all_paths[0].current_s);
        // G_outer.info();
    }

    return last_G;
}

auto findpath::process_int_loops(int_loops current_sections, short* pt_1, short* pt_2) -> s_graph
{
    // recursive processing of sections - findpath and merge calls below

    // bool Verbose = false;
    bool Verbose = true;

    fmt::print ("start\n");
    int bp_dist = 0;

    size_t moves_hash = 0;

    std::vector<move_ij> move_list{};

    for (int i = 1; i <= pt_1[0]; i++) {
        if (pt_1[i] != pt_2[i]) {
            if (i < pt_1[i]) {  // need to delete this pair
                move_list.push_back({static_cast<short>(-i), static_cast<short>(-pt_1[i])});

                size_t current_hash = int_hash(-i);
                current_hash *= (-pt_1[i]);
                current_hash += (-pt_1[i]) * 13;
                moves_hash += current_hash;
                bp_dist++;
            }
            if (i < pt_2[i]) {  // need to insert this pair
                move_list.push_back({static_cast<short>(i), static_cast<short>(pt_2[i])});

                size_t current_hash = int_hash(i);
                current_hash *= (pt_2[i]);
                current_hash += (pt_2[i]) * 13;
                moves_hash += current_hash;
                bp_dist++;
            }
        }
    }
    move_list.push_back({0, 0});

    if (bp_dist != 0 and G_cache.contains(moves_hash)) {
        fmt::print("G avail: {} / {}\n", moves_hash, bp_dist);
        
        s_graph temp = G_cache[moves_hash];
        // temp.info();
        // temp.display_path();

        fmt::print("return\n");
        // s_graph temp;

        temp.pt_1 = pt_1;
        temp.pt_2 = pt_2;
        fmt::print("return2\n");

        return temp;
        // G_cache[moves_hash].info();
        // return G_cache[moves_hash];
    }

    // continue with recursive processing of sections...



    // int  search_width_multiplier = 1;
    // float search_width_multiplier       = 2;
    int merge_search_width_multiplier = 1;

    if (current_sections.nested_sections.size() == 0) {
        // we're at the top of the recursion tree - findpath for a single section with start and end

        fmt::print ("top\n");

        const int search_width = current_sections.bp_dist * search_width_multiplier;

        single_findpath fp_call;
        auto            result = fp_call.init(fc, pt_1, pt_2, search_width, true);

        // postprocess paths into graph
        s_graph G_inner{fc, pt_1, pt_2, current_sections.bp_dist, result};
        G_inner.max_en = result[0].max_en;

        // G_inner.info();
        if (Verbose) print_moves(result[0], fc, vrna_db_from_ptable(pt_1), true);

        G_inner.info();
        if (cache) { G_cache[moves_hash] = G_inner; }
        return G_inner;
    }

    fmt::print ("B\n");

    // else: generate the outer pairing table and substract all inner sections, where we will
    // recurse into

    int outer_bp_dist = current_sections.bp_dist;

    // generate pairing table for outer sections  from X to Y: [X, [39, 64], [124, 134], Y]
    short* outer_pt_1 = vrna_ptable_copy(pt_1);
    short* outer_pt_2 = vrna_ptable_copy(pt_2);
    for (int i = 1; i <= pt_1[0]; i++) {
        if (i < current_sections.start or i > current_sections.end) {
            outer_pt_1[i] = 0;
            outer_pt_2[i] = 0;
        }
    }

    // iterate over nested inner sections to remove everything for the outer pairing
    // table (e.g. remove [39, 64] and [124, 134])
    for (const auto& nested_section : current_sections.nested_sections) {
        outer_bp_dist -= nested_section.bp_dist;

        for (int i = 1; i <= pt_1[0]; i++) {
            if (i > nested_section.start and i < nested_section.end) {
                outer_pt_1[i] = 0;
                outer_pt_2[i] = 0;
            }
        }
    }

    // fmt::print ("s1: {}\n", vrna_db_from_ptable(outer_pt_1));
    // fmt::print ("s2: {}\n", vrna_db_from_ptable(outer_pt_2));

    // findpath call for outer section
    const int search_width = outer_bp_dist * search_width_multiplier;

    // intermission:

    size_t moves_hash_outer = 0;
    for (int i = 1; i <= outer_pt_1[0]; i++) {
        if (outer_pt_1[i] != pt_2[i]) {
            if (i < outer_pt_1[i]) {  // need to delete this pair
                size_t current_hash = int_hash(-i);
                current_hash *= (-outer_pt_1[i]);
                current_hash += (-outer_pt_1[i]) * 13;
                moves_hash_outer += current_hash;
            }
            if (i < outer_pt_2[i]) {  // need to insert this pair
                size_t current_hash = int_hash(i);
                current_hash *= (outer_pt_2[i]);
                current_hash += (outer_pt_2[i]) * 13;
                moves_hash_outer += current_hash;
            }
        }
    }

    s_graph G_outer;

    if (cache and G_cache.contains(moves_hash_outer)) {
        fmt::print ("G avail o.: {} \n", moves_hash_outer);
        G_outer = G_cache[moves_hash_outer];
    } else {
        single_findpath fp_call;
        auto            result = fp_call.init(fc, outer_pt_1, outer_pt_2, search_width, true);
        // postprocessing: outer paths into graph
        G_outer = s_graph {fc, outer_pt_1, outer_pt_2, outer_bp_dist, result};

        if (Verbose) print_moves(result[0], fc, vrna_db_from_ptable(outer_pt_1), true);
    }
    // G_outer.info();

    // iterate over nested inner sections and merge them one by one into the outer section
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

        int max_en             = INT_MAX;
        int merge_search_width = total_bp_dist * merge_search_width_multiplier;

        auto start = std::chrono::high_resolution_clock::now();

        // G_inner.bp_dist = b;
        int c           = bp_distance(outer_pt_1, outer_pt_2);
        G_outer.bp_dist = c;

        // this is the only merge call
        auto G_merged = merge_method(G_outer, G_inner, merged_pt_1, s1_en, merged_pt_2, s2_en,
                                     total_bp_dist, max_en, merge_search_width, fc, true);

        auto                          finish  = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;

        // fmt::print("merged: {} \n", G_merged.max_en);

        // adjustments to merge the next inner section into the outer section + current inner
        // section
        outer_bp_dist += G_inner.bp_dist;
        outer_pt_1 = merged_pt_1;
        outer_pt_2 = merged_pt_2;
        G_outer    = G_merged;

        // if (Verbose) print_moves(all_paths[0], fc, vrna_db_from_ptable(merged_pt_1), true);
        // print_moves(all_paths[0], fc, vrna_db_from_ptable(merged_pt_1), false);
    }

    if (Verbose) { G_outer.display_path(); }

    if (cache) { G_cache[moves_hash] = G_outer; }
    return G_outer;
}
