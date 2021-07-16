

class findpath
{
   private:
    int  test = 0;
    void init_pt();

    // fp cache

    auto single_findpath_cache(vrna_fold_compound_t* fc, short* pt1, short* pt2, int search_width,
                               bool mp);

    auto      process_ext_loops(ext_loops ext_loops, short* pt_1, short* pt_2) -> s_graph;
    auto      process_int_loops(int_loops current_sections, std::vector<short> pt1,
                                std::vector<short> pt2) -> s_graph;
    int_loops all_sections{0, 0};

    s_graph                             G;  // result graph
    std::unordered_map<size_t, s_graph> G_cache;
    bool                                cache = true;

   public:
    float search_width_multiplier;
    bool  mp;

    // todo: delete these...
    // char *                seq, *s1, *s2;
    // short *               pt_1, *pt_2;
    vrna_fold_compound_t* fc;

    auto init(std::string s1, std::string s2, float sw = 2) -> s_graph;
    auto init_python(std::string s1, std::string s2, float sw = 2) -> void;

    // result getters
    auto get_path() -> std::vector<std::tuple<int, int, int>>;
    auto get_sections() -> void;
    auto get_en() -> int; 

    auto init_ext(std::string s1, std::string s2, float sw = 2) -> s_graph;

    // various constructors (with fold compound or sequence)

    findpath(std::string sequence, bool mp = true, const py::dict &model_details = {});

    ~findpath() { vrna_fold_compound_free(fc); }
};

findpath::findpath(std::string sequence, bool mp, const py::dict &model_details) : mp{mp}
{
    vrna_md_t md;
    
    vrna_md_set_default(&md); // copy global settings
    // set_model_details(&md);
    

    // key and values of the dict are pybind11 objects, with included cast function to convert them as needed.

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

    // s1          = vrna_db_from_ptable(pt_1);
    // s2          = vrna_db_from_ptable(pt_2);
}

auto findpath::init_python(std::string s1, std::string s2, float sw) -> void
{
    // class member function, exported to Python

    search_width_multiplier = sw;
    std::vector<short> pt1  = ptable_from_string(s1);
    std::vector<short> pt2  = ptable_from_string(s2);

    // initialise section structure
    all_sections.start = 1;
    all_sections.end   = pt1[0];

    all_sections.bp_dist = vrna_bp_distance(s1.c_str(), s2.c_str());
    // recursively add all nested sections (if available)
    all_sections.nested_sections = find_interior_loops(pt1, pt2, 1, pt1[0]);

    // short* pt_1 = pt1.data();
    // short* pt_2 = pt2.data();

    // recursively process nested sections
    G = process_int_loops(all_sections, pt1, pt2);

    // return G.max_en;
}

auto findpath::get_en() -> int {
    return G.max_en;
}

auto findpath::get_path() -> std::vector<std::tuple<int, int, int>>
{
    std::vector<std::tuple<int, int, int>> path = G.return_path();

    return path;
    // class member function, exported to Python
}

auto findpath::get_sections() -> void { std::cout << all_sections << "\n"; }

auto findpath::init(std::string s1, std::string s2, float sw) -> s_graph
{
    // class member function, exported to Python

    search_width_multiplier = sw;
    std::vector<short> pt1  = ptable_from_string(s1);
    std::vector<short> pt2  = ptable_from_string(s2);

    // initialise section structure
    all_sections.start = 1;
    all_sections.end   = pt1[0];

    all_sections.bp_dist = vrna_bp_distance(s1.c_str(), s2.c_str());
    // recursively add all nested sections (if available)
    all_sections.nested_sections = find_interior_loops(pt1, pt2, 1, pt1[0]);

    // short* pt_1 = pt1.data();
    // short* pt_2 = pt2.data();

    // recursively process nested sections
    G = process_int_loops(all_sections, pt1, pt2);

    return G;
}

auto findpath::init_ext(std::string s1, std::string s2, float sw) -> s_graph
{
    search_width_multiplier = sw;
    std::vector<short> pt1  = ptable_from_string(s1);
    std::vector<short> pt2  = ptable_from_string(s2);
    // variant with exterior loops
    auto   ext_loops = find_exterior_loops(pt1, pt2);
    short* pt_1      = pt1.data();
    short* pt_2      = pt2.data();
    auto   G         = process_ext_loops(ext_loops, pt_1, pt_2);
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
                     std::vector<std::vector<short>>& ptable_storage, const s_graph& G1,
                     const s_graph& G2) -> void
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

    std::vector<std::vector<sorted_move>> move_storage(2 * merge_search_width + 2);
    for (auto& col : move_storage) { col.reserve(total_bp_dist); }
    std::vector<std::vector<short>> ptable_storage(2 * merge_search_width + 2);
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

        // fmt::print("size after sort: {} \n", next_paths.size());

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
            //    pool_offset, current_path.move_ptr == nullptr, i_move, j_move, last_index);

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
            //            pool_offset, current_path.current_ptable == nullptr, i_move, j_move,
            //            last_index);

            // update pairing tables
            const int len           = current_path.current_ptable[0];
            new_path.current_ptable = ptable_storage[i + pool_offset].data();
            memcpy(new_path.current_ptable, current_path.current_ptable, sizeof(short) * (len + 1));

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

auto merge_method(auto& G1, auto& G2, std::vector<short> pt1, int s1_en, std::vector<short> pt2, int s2_en,
                  int total_bp_dist, int max_en, int final_merge_search_width,
                  vrna_fold_compound_t* fc, bool merge_constant_sections)
{
    // merge 2 graphs, increase search widths, consider both directions, etc.

    // workaround for ext_loops version
    if (merge_constant_sections) {
        // we might merge constant sections, ignore them
        if (G1.bp_dist == 0) {
            // fmt::print("ignore merge \n");
            G2.pt1 = pt1;
            G2.pt2 = pt2;
            return G2;
        }
        if (G2.bp_dist == 0) {
            // fmt::print("ignore merge \n");
            G1.pt1 = pt1;
            G1.pt2 = pt2;
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
        // direction                         = true;
        // std::vector<merge_path> fwd_paths = merge_once(
        //     G1, G2, pt_1, s1_en, direction, total_bp_dist, max_en, current_merge_search_width,
        //     fc);
        // direction                         = false;
        // std::vector<merge_path> bwd_paths = merge_once(
        //     G1, G2, pt_2, s2_en, direction, total_bp_dist, max_en, current_merge_search_width,
        //     fc);

        std::future<std::vector<merge_path>> ret1 =
            std::async(std::launch::async, &merge_once, G1, G2, pt1.data(), s1_en, true, total_bp_dist,
                       max_en, current_merge_search_width, fc);
        auto                                 G1_copy = G1;
        auto                                 G2_copy = G2;
        std::future<std::vector<merge_path>> ret2 =
            std::async(std::launch::async, &merge_once, G1_copy, G2_copy, pt2.data(), s2_en, false,
                       total_bp_dist, max_en, current_merge_search_width, fc);
        std::vector<merge_path> fwd_paths = ret1.get();
        std::vector<merge_path> bwd_paths = ret2.get();

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

    // std::vector<short> pt1(total_bp_dist+2);
    // std::vector<short> pt2(total_bp_dist+2);

    // pt1


    s_graph G_merged{fc, pt1, pt2, total_bp_dist, all_paths, max_en};

    return G_merged;
}

auto findpath::process_ext_loops(ext_loops ext_loops, short* pt_1, short* pt_2) -> s_graph
{
    // this needs some work...

    s_graph test;
    return test;

    // bool Verbose = false;
    // // bool Verbose                       = true;
    // int merge_search_width_multiplier = 1;

    // bool    first_run = true;
    // s_graph last_G;
    // short*  last_pt_1;
    // short*  last_pt_2;
    // int     last_bp_dist = -1;

    // // iterate over available ext loops

    // for (int i = 0; i < ext_loops.size(); i++) {
    //     const auto& current_ext_loop = ext_loops[i];

    //     // std::cout << "current ext loop:" << current_ext_loop << " / "
    //     //           << "\n";
    //     // fmt::print("process ext loops {}\n", i);

    //     short* current_pt_1 = vrna_ptable_copy(pt_1);
    //     short* current_pt_2 = vrna_ptable_copy(pt_2);

    //     // fmt::print("s1 {}\n", vrna_db_from_ptable(pt_1));
    //     // fmt::print("s2 {}\n", vrna_db_from_ptable(pt_2));

    //     // adjust pairing table
    //     for (int i = 1; i <= pt_1[0]; i++) {
    //         if (i < current_ext_loop.start or i > current_ext_loop.end) {
    //             current_pt_1[i] = 0;
    //             current_pt_2[i] = 0;
    //         }
    //     }
    //     int current_bp_dist = bp_distance(current_pt_1, current_pt_2);
    //     // current_G.bp_dist = current_bp_dist;

    //     // fmt::print("s1 {}\n", vrna_db_from_ptable(current_pt_1));
    //     // fmt::print("s2 {}\n", vrna_db_from_ptable(current_pt_2));

    //     // ext loops without further interior loop recursion
    //     // const int search_width = current_bp_dist * search_width_multiplier;
    //     // single_findpath test;
    //     // auto            result = test.init(fc, current_pt_1, current_pt_2, search_width,
    //     true);
    //     // s_graph current_G {fc, current_pt_1, current_pt_2, current_bp_dist, result};
    //     // current_G.max_en = result[0].max_en;

    //     s_graph current_G = process_int_loops(current_ext_loop, current_pt_1, current_pt_2);
    //     current_G.bp_dist = current_bp_dist;

    //     // fmt::print("ext: current bp dist: {}\n", current_bp_dist);

    //     if (first_run) {
    //         last_G       = current_G;
    //         last_pt_1    = current_pt_1;
    //         last_pt_2    = current_pt_2;
    //         first_run    = false;
    //         last_bp_dist = current_bp_dist;
    //         continue;
    //     }

    //     int merge_bp_dist      = last_G.bp_dist + current_G.bp_dist;
    //     int max_en             = INT_MAX;
    //     int merge_search_width = merge_bp_dist * merge_search_width_multiplier;

    //     // merge pairing tables (outer and current inner)
    //     short* merged_pt_1 = vrna_ptable_copy(current_pt_1);
    //     short* merged_pt_2 = vrna_ptable_copy(current_pt_2);

    //     for (int i = 1; i <= merged_pt_1[0]; i++) {
    //         // fmt::print("i {}\n", i);
    //         if ((merged_pt_1[i] == 0 and last_pt_1 != 0) or
    //             (merged_pt_1[i] != 0 and last_pt_1 == 0)) {
    //             merged_pt_1[i] = last_pt_1[i];
    //         }
    //         if ((merged_pt_2[i] == 0 and last_pt_2 != 0) or
    //             (merged_pt_2[i] != 0 and last_pt_2 == 0)) {
    //             merged_pt_2[i] = last_pt_2[i];
    //         }
    //     }

    //     // current_G.display_path();
    //     // current_G.info();
    //     // last_G.display_path();
    //     // last_G.info();

    //     int s1_en = vrna_eval_structure_pt(fc, merged_pt_1);
    //     int s2_en = vrna_eval_structure_pt(fc, merged_pt_2);

    //     auto start = std::chrono::high_resolution_clock::now();

    //     auto G_merged = merge_method(current_G, last_G, merged_pt_1, s1_en, merged_pt_2, s2_en,
    //                                  merge_bp_dist, max_en, merge_search_width, fc, false);

    //     auto                          finish  = std::chrono::high_resolution_clock::now();
    //     std::chrono::duration<double> elapsed = finish - start;
    //     // IC(elapsed.count());

    //     // fmt::print("s1 {}\n", vrna_db_from_ptable(merged_pt_1));
    //     // fmt::print("s2 {}\n", vrna_db_from_ptable(merged_pt_2));
    //     // fmt::print("s1 {}\n", vrna_db_from_ptable(G_merged.pt_1));
    //     // fmt::print("s2 {}\n", vrna_db_from_ptable(G_merged.pt_2));
    //     // fmt::print("merge time: {} \n", elapsed.count());

    //     last_pt_1    = merged_pt_1;
    //     last_pt_2    = merged_pt_2;
    //     last_bp_dist = current_bp_dist;
    //     last_G       = G_merged;
    //     // fmt::print("display merged path:\n");
    //     // last_G.display_path();
    //     // fmt::print("merged: {} {}\n", all_paths.size(), all_paths[0].current_s);
    //     // G_outer.info();
    // }

    // return last_G;
}

auto findpath::process_int_loops(int_loops current_sections, std::vector<short> pt1,
                                 std::vector<short> pt2) -> s_graph
{
    // recursive processing of sections - findpath and merge calls below

    bool Verbose = false;
    // bool Verbose = true;

    if (Verbose) fmt::print("start\n");
    int bp_dist = 0;

    size_t moves_hash = 0;

    std::vector<move_ij> move_list{};

    for (int i = 1; i <= pt1[0]; i++) {
        if (pt1[i] != pt2[i]) {
            if (i < pt1[i]) {  // need to delete this pair
                move_list.push_back({static_cast<short>(-i), static_cast<short>(-pt1[i])});

                size_t current_hash = int_hash(-i);
                current_hash *= (-pt1[i]);
                current_hash += (-pt1[i]) * 13;
                moves_hash += current_hash;
                bp_dist++;
            }
            if (i < pt2[i]) {  // need to insert this pair
                move_list.push_back({static_cast<short>(i), static_cast<short>(pt2[i])});

                size_t current_hash = int_hash(i);
                current_hash *= (pt2[i]);
                current_hash += (pt2[i]) * 13;
                moves_hash += current_hash;
                bp_dist++;
            }
        }
    }
    move_list.push_back({0, 0});

    if (bp_dist != 0 and G_cache.contains(moves_hash)) {
        // fmt::print("G avail: {} / {}\n", moves_hash, bp_dist);

        s_graph temp = G_cache[moves_hash];
        // temp.info();
        // temp.display_path();

        // fmt::print("return\n");
        // s_graph temp;

        temp.pt1 = pt1;
        temp.pt2 = pt2;
        // fmt::print("return2\n");

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

        // fmt::print ("top\n");

        int search_width = current_sections.bp_dist * search_width_multiplier;
        if (search_width == 0) { search_width = 1; }
        single_findpath fp_call;
        auto            result = fp_call.init(fc, pt1.data(), pt2.data(), search_width, true);

        // postprocess paths into graph
        // s_graph G_inner{fc, pt_1, pt_2, current_sections.bp_dist, result};
        s_graph G_inner{fc, pt1, pt2, current_sections.bp_dist, result, result[0].max_en};

        // G_inner.max_en = result[0].max_en;

        // G_inner.info();
        if (Verbose) print_moves(result[0], fc, pt1.data(), true);

        // G_inner.info();
        if (cache) { G_cache[moves_hash] = G_inner; }

        if (Verbose) fmt::print("return G_inner\n");

        return G_inner;
    }

    // fmt::print ("B\n");

    // else: generate the outer pairing table and substract all inner sections, where we will
    // recurse into

    int outer_bp_dist = current_sections.bp_dist;

    // generate pairing table for outer sections  from X to Y: [X, [39, 64], [124, 134], Y]
    // short* outer_pt_1 = vrna_ptable_copy(pt_1);
    // short* outer_pt_2 = vrna_ptable_copy(pt_2);

    // this is a memcpy = deepcopy
    std::vector<short> outer_pt1 = pt1;
    std::vector<short> outer_pt2 = pt2;

    for (int i = 1; i <= pt1[0]; i++) {
        if (i < current_sections.start or i > current_sections.end) {
            outer_pt1[i] = 0;
            outer_pt2[i] = 0;
        }
    }

    // iterate over nested inner sections to remove everything for the outer pairing
    // table (e.g. remove [39, 64] and [124, 134])
    for (const auto& nested_section : current_sections.nested_sections) {
        outer_bp_dist -= nested_section.bp_dist;

        for (int i = 1; i <= pt1[0]; i++) {
            if (i > nested_section.start and i < nested_section.end) {
                outer_pt1[i] = 0;
                outer_pt2[i] = 0;
            }
        }
    }

    if (Verbose) fmt::print("s1: {}\n", vrna_db_from_ptable(outer_pt1.data()));
    if (Verbose) fmt::print("s2: {}\n", vrna_db_from_ptable(outer_pt2.data()));

    // findpath call for outer section
    int search_width = outer_bp_dist * search_width_multiplier;

    // intermission:

    size_t moves_hash_outer = 0;
    for (int i = 1; i <= outer_pt1[0]; i++) {
        if (outer_pt1[i] != pt2[i]) {
            if (i < outer_pt1[i]) {  // need to delete this pair
                size_t current_hash = int_hash(-i);
                current_hash *= (-outer_pt1[i]);
                current_hash += (-outer_pt1[i]) * 13;
                moves_hash_outer += current_hash;
            }
            if (i < outer_pt2[i]) {  // need to insert this pair
                size_t current_hash = int_hash(i);
                current_hash *= (outer_pt2[i]);
                current_hash += (outer_pt2[i]) * 13;
                moves_hash_outer += current_hash;
            }
        }
    }

    s_graph G_outer;

    if (cache and G_cache.contains(moves_hash_outer)) {
        // fmt::print ("G avail o.: {} \n", moves_hash_outer);
        G_outer = G_cache[moves_hash_outer];
    } else {
        if (search_width == 0) { search_width = 1; }
        single_findpath fp_call;
        auto result = fp_call.init(fc, outer_pt1.data(), outer_pt2.data(), search_width, true);
        // postprocessing: outer paths into graph
        G_outer = s_graph{fc, outer_pt1, outer_pt2, outer_bp_dist, result, result[0].max_en};

        if (Verbose) print_moves(result[0], fc, outer_pt1.data(), true);
    }
    // G_outer.info();

    // iterate over nested inner sections and merge them one by one into the outer section
    for (const auto& nested_section : current_sections.nested_sections) {
        // short* inner_pt_1 = vrna_ptable_copy(pt_1);
        // short* inner_pt_2 = vrna_ptable_copy(pt_2);
        std::vector<short> inner_pt1 = pt1;
        std::vector<short> inner_pt2 = pt2;

        for (int i = 1; i <= pt1[0]; i++) {
            if (i < nested_section.start or i > nested_section.end) {
                inner_pt1[i] = 0;
                inner_pt2[i] = 0;
            }
        }
        if (Verbose) fmt::print("inner section 1:\n{}\n", vrna_db_from_ptable(inner_pt1.data()));
        if (Verbose) fmt::print("inner section 2:\n{}\n", vrna_db_from_ptable(inner_pt2.data()));

        // recursive call to build up the current inner graph
        s_graph G_inner = process_int_loops(nested_section, inner_pt1, inner_pt2);

        // merge pairing tables (outer and current inner)
        // short* merged_pt_1 = vrna_ptable_copy(outer_pt_1);
        // short* merged_pt_2 = vrna_ptable_copy(outer_pt_2);
        std::vector<short> merged_pt1 = outer_pt1;
        std::vector<short> merged_pt2 = outer_pt2;

        // for (int i = 1; i < merged_pt1[0]; i++) {
        //     // fmt::print("i {}\n", i);
        //     if ((merged_pt1[i] == 0 and G_inner.pt_1 != 0) or
        //         (merged_pt1[i] != 0 and G_inner.pt_1 == 0)) {
        //         merged_pt1[i] = G_inner.pt_1[i];
        //     }
        //     if ((merged_pt2[i] == 0 and G_inner.pt_2 != 0) or
        //         (merged_pt2[i] != 0 and G_inner.pt_2 == 0)) {
        //         merged_pt2[i] = G_inner.pt_2[i];
        //     }
        // }

        // if ptB has base pairs which are not present in ptA, add them
        for (int i = 1; i < merged_pt1[0]; i++) {
            // fmt::print("i {}\n", i);
            if (merged_pt1[i] == 0 and inner_pt1[i] != 0) { merged_pt1[i] = inner_pt1[i]; }
            if (merged_pt2[i] == 0 and inner_pt2[i] != 0) { merged_pt2[i] = inner_pt2[i]; }
        }

        // fmt::print("inner1 {}\n", vrna_db_from_ptable(inner_pt1.data()));
        // fmt::print("outer1 {}\n", vrna_db_from_ptable(outer_pt1.data()));
        // fmt::print("merge1 {}\n", vrna_db_from_ptable(merged_pt1.data()));
        // fmt::print("inner2 {}\n", vrna_db_from_ptable(inner_pt2.data()));
        // fmt::print("outer2 {}\n", vrna_db_from_ptable(outer_pt2.data()));
        // fmt::print("merge2 {}\n", vrna_db_from_ptable(merged_pt2.data()));


        // for (int i = 1; i < merged_pt1[0]; i++) { fmt::print("{} ", pt2test[i]); }
        // fmt::print("\n");

        int s1_en = vrna_eval_structure_pt(fc, merged_pt1.data());
        int s2_en = vrna_eval_structure_pt(fc, merged_pt2.data());

        // fmt::print("to2p\n");

        // fmt::print("starting merge:\n{} {} {}\n", vrna_db_from_ptable(merged_pt1.data()),
        //            G_outer.bp_dist, s1_en);
        // fmt::print("{} {} {}\n", vrna_db_from_ptable(merged_pt2.data()), G_inner.bp_dist, s2_en);

        int total_bp_dist = G_inner.bp_dist + G_outer.bp_dist;

        int max_en             = INT_MAX;
        int merge_search_width = total_bp_dist * merge_search_width_multiplier;

        auto start = std::chrono::high_resolution_clock::now();

        // G_inner.bp_dist = b;
        int c           = bp_distance(outer_pt1, outer_pt2);
        G_outer.bp_dist = c;

        // this is the only merge call
        auto G_merged = merge_method(G_outer, G_inner, merged_pt1, s1_en, merged_pt2,
                                     s2_en, total_bp_dist, max_en, merge_search_width, fc, true);

        auto                          finish  = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;

        if (Verbose) fmt::print("merged: {} \n", G_merged.max_en);

        // free(inner_pt_1);
        // free(inner_pt_2);
        // free(outer_pt_1);
        // free(outer_pt_2);

        // adjustments to merge the next inner section into the outer section +
        // current inner

        // section
        outer_bp_dist += G_inner.bp_dist;
        outer_pt1 = merged_pt1;
        outer_pt2 = merged_pt2;
        G_outer   = G_merged;

        // if (Verbose) print_moves(all_paths[0], fc, vrna_db_from_ptable(merged_pt_1), true);
        // print_moves(all_paths[0], fc, vrna_db_from_ptable(merged_pt_1), false);
    }

    if (Verbose) { G_outer.display_path(); }

    if (cache) { G_cache[moves_hash] = G_outer; }
    return G_outer;
}
