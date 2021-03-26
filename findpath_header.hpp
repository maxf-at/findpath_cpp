

struct section {
    int                  start;
    int                  end;
    int                  bp_dist;
    std::vector<section> nested_sections{};

    // constructor
    section(int a, int b)
    {
        start = a;
        end   = b;
    }
    section(int a, int b, int c)
    {
        start   = a;
        end     = b;
        bp_dist = c;
    }

    std::ostream& operator<<(std::ostream& out) { return out << '[' << start << " " << end << ']'; }
};

// print funcion for section
std::ostream& operator<<(std::ostream& os, const section& s)
{
    os << '[' << s.start << ", ";
    for (const auto& r : s.nested_sections) {
        os << r;  // recursive ostream call
        os << ", ";
    }
    os << s.end << ']';
    return os;
}

class findpath
{
   private:
    int                  test = 0;
    void                 init_pt();
    void                 init_sections();
    std::vector<section> find_interior_loop(int start, int end);
    s_graph              findpath_sections(section current_sections, short* pt_1, short* pt_2);
    section              all_sections{0, 0};

   public:
    char *                seq, *s1, *s2;
    short *               pt_1, *pt_2;
    vrna_fold_compound_t* fc = nullptr;
    std::tuple<int, int>  bp_distance(short* pt_a, short* pt_b, int min_pos = 0,
                                      int max_pos = INT_MAX / 2);

    // Constructors with fold compound or sequence
    findpath(vrna_fold_compound_t* init_fc, char* a_s1, char* a_s2);
    findpath(char* seq, char* s1, char* s2);
    // s_graph(int size); //Constructor of the class
};

// Constructors
findpath::findpath(vrna_fold_compound_t* a_fc, char* a_s1, char* a_s2)
{
    std::cout << "Findpath object being created\n";
    this->fc = a_fc;
    // fc = a_fc;
    s1 = a_s1;
    s2 = a_s2;

    init_pt();
    init_sections();

    // graphSize = size;
    // cout << "Graph object created"<<endl;
}

findpath::findpath(char* a_seq, char* a_s1, char* a_s2)
{
    seq = a_seq;
    s1  = a_s1;
    s2  = a_s2;

    // set model params
    vrna_md_t md;
    set_model_details(&md);
    fc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);

    init_pt();
    init_sections();

    // std::cout << "Findpath seq object done\n";
}

std::tuple<int, int> findpath::bp_distance(short* pt_a, short* pt_b, int min_pos, int max_pos)
{
    /*
    basepair distance between 2 pairing tables (see vrna_bp_distance() )
     */

    // int dist = 0;

    int inner_bp_dist = 0;
    int outer_bp_dist = 0;

    int length = pt_a[0];

    for (int i = 1; i <= length; i++)
        if (pt_a[i] != pt_b[i]) {
            if (pt_a[i] > i) {
                if (i > min_pos and i < max_pos)
                    inner_bp_dist++;
                else
                    outer_bp_dist++;
            }
            // both can happen at once...
            if (pt_b[i] > i) {
                if (i > min_pos and i < max_pos)
                    inner_bp_dist++;
                else
                    outer_bp_dist++;
            }
        }

    return {inner_bp_dist, outer_bp_dist};
}

void findpath::init_pt()
{
    pt_1 = vrna_ptable(s1);
    pt_2 = vrna_ptable(s2);
}

void findpath::init_sections()
{
    std::cout << "Findpath seq object done\n";

    // find potential sections to split up the findpath process
    // section all_sections(1, pt_1[0]);  // initial large section

    all_sections.start = 1;
    all_sections.end   = pt_1[0];

    all_sections.bp_dist = vrna_bp_distance(s1, s2);
    all_sections.nested_sections =
        find_interior_loop(1, pt_1[0]);  // recursively add all nested sections (if available)

    std::cout << "Sections:" << all_sections << "\n";

    findpath_sections(all_sections, pt_1, pt_2);
}

std::vector<section> findpath::find_interior_loop(int min_pos, int max_pos)
{
    // IC(min_pos, max_pos);

    std::vector<std::tuple<int, int, int, int, float>> candidates;

    for (int i = min_pos + 1; i < max_pos - 1; i++) {
        const int j1 = pt_1[i];
        const int j2 = pt_2[i];

        if (i < min_pos or i + 1 > max_pos) continue;
        if (j1 == 0 and j2 == 0) continue;
        if (j1 != j2 or i > j1) continue;

        const int j = j1;

        // short* inner_p_table1 = vrna_ptable_copy(pt_1);
        // short* inner_p_table2 = vrna_ptable_copy(pt_2);
        // short* outer_p_table1 = vrna_ptable_copy(pt_1);
        // short* outer_p_table2 = vrna_ptable_copy(pt_2);
        // for (int k = 1; k < pt_1[0]; k++) {
        //     if (k < i or k > j + 1) {
        //         inner_p_table1[k] = 0;
        //         inner_p_table2[k] = 0;
        //     } else {
        //         outer_p_table1[k] = 0;
        //         outer_p_table2[k] = 0;
        //     }
        // }
        // char* inner_s1 = vrna_db_from_ptable(inner_p_table1);
        // char* inner_s2 = vrna_db_from_ptable(inner_p_table2);
        // char* outer_s1 = vrna_db_from_ptable(outer_p_table1);
        // char* outer_s2 = vrna_db_from_ptable(outer_p_table2);
        // int inner_bp_dist = vrna_bp_distance(inner_s1, inner_s2);
        // int outer_bp_dist = vrna_bp_distance(outer_s1, outer_s2);
        // IC(i, j, inner_bp_dist, outer_bp_dist);

        // inner / outer basepair distance (between ptable1 and 2 with i and j as limits)
        auto [inner_bp_dist, outer_bp_dist] = bp_distance(pt_1, pt_2, i, j);
        // outer_bp_dist = bp_distance(outer_p_table1, outer_p_table2);

        if (std::min(inner_bp_dist, outer_bp_dist) < 1) continue;

        const float inner_size = j - i;
        const float outer_size = max_pos - min_pos - inner_size;
        const float optimize   = std::abs(0.7 - (inner_size / outer_size));

        // IC(i, j, inner_bp_dist, outer_bp_dist, optimize);

        candidates.emplace_back(i, j, inner_bp_dist, outer_bp_dist, optimize);

        // break;
    }

    // minimize the outer basepair distance - this maximizes recursive potential
    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
        if (get<3>(a) == get<3>(b)) return get<4>(a) < get<4>(b);
        return get<3>(a) < get<3>(b);
    });

    std::vector<section> nested_sections{};
    int                  bp_dist_left = 9999;

    // for (const auto& candidate : candidates) {
    for (const auto& [i, j, inner_bp_dist, outer_bp_dist, optimize] : candidates) {
        // IC(nested_sections, i, j, bp_dist_left);

        if (bp_dist_left - inner_bp_dist < 2) continue;  // nothing left

        // section i to j can't be part of the already existing section
        // dont add something twice
        bool ignore = false;
        for (const auto& section : nested_sections) {
            if ((i > section.start and j < section.end) or (i < section.start and j > section.end))
                ignore = true;
        }
        if (ignore) continue;

        if (bp_dist_left == 9999)
            bp_dist_left = outer_bp_dist;
        else
            bp_dist_left -= inner_bp_dist;

        section current_section(i, j, inner_bp_dist);

        // recursive call
        current_section.nested_sections = {};
        // current_section.nested_sections = find_interior_loop(i + 1, j - 1);

        nested_sections.push_back(current_section);

        // IC(i, j);
    }

    return nested_sections;
}

// merge routine

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

    // std::format("{} {}!", "Hello", "world", "something"); // OK, produces "Hello world!"

    short* pt;
    pt = vrna_ptable(s1);

    // std::cout << fmt::format("Hello!\n");

    float en     = vrna_eval_structure_pt(fc, pt) / 100.0;
    float max_en = -INT_MAX;

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

auto available_edges(const auto& G1_node, const auto& G2_node, int G1_G2, auto& next_paths,
                     const auto& current_path, vrna_fold_compound_t* fc, int max_en)
{
    for (const auto& current_edge : G1_node.out_edges) {
        const auto& i = current_edge.i;
        const auto& j = current_edge.j;
        auto        current_en =
            current_path.current_en + vrna_eval_move_pt(fc, current_path.current_ptable, i, j);

        if (current_en <= max_en) {
            // generate new pairing tables
            // auto current_ptable = vrna_ptable_copy(pt_1);
            merge_path new_path;
            if (G1_G2 == 1) {
                new_path.current_G1_node = current_edge.destination;
                new_path.current_G2_node = current_path.current_G2_node;
            } else {  // reversed
                new_path.current_G1_node = current_path.current_G1_node;
                new_path.current_G2_node = current_edge.destination;
            }

            new_path.current_en     = current_en;
            new_path.current_s      = std::max(current_en, current_path.current_s);
            new_path.current_ptable = vrna_ptable_copy(current_path.current_ptable);
            new_path.moves          = current_path.moves;
            new_path.moves.push_back({i, j, current_en});

            if (j < 0) {  // delete a basepair
                new_path.current_ptable[-i] = 0;
                new_path.current_ptable[-j] = 0;
            } else {  // add a basepair
                new_path.current_ptable[i] = j;
                new_path.current_ptable[j] = i;
            }

            next_paths.push_back(new_path);
            // auto en = vrna_eval_structure_pt(fc, new_path.current_ptable);
            // std::string s = vrna_db_from_ptable(new_path.current_ptable);
            // fmt::print("{} d={} G1 {} {} {} {}\n", s, d, i, j, en, current_en);
        }
    }
}

auto merge_once(auto G1, auto G2, short* pt_1, int s1_en, short* pt_2, int s2_en, bool direction,
                int total_bp_dist, int max_en, int merge_search_width, vrna_fold_compound_t* fc)
{
    // this is the start of findpath_once, assuming we have maxE etc...

    int current_i_node = 0;
    int current_j_node = 0;

    if (!direction) {
        current_i_node = G1.bp_dist;
        current_i_node = G2.bp_dist;
    }

    // init single empty path
    std::vector<merge_path> all_paths;

    merge_path current_path;
    current_path.current_G1_node = 0;
    current_path.current_G2_node = 0;
    current_path.current_en      = s1_en;
    current_path.current_s       = s1_en;
    current_path.current_ptable  = pt_1;

    all_paths.push_back(current_path);

    const auto& G1_node = G1.node_list[current_path.current_G1_node];
    const auto& G2_node = G2.node_list[current_path.current_G2_node];

    if (direction) { const auto& edges = G1_node.out_edges; }

    // fmt::print("moves {} {}\n", G1_node.out_edges[0].i, G1_node.out_edges[0].j);
    // fmt::print("moves0 {} {}\n", G1.node_list[0].out_edges[0].i, G1.node_list[0].out_edges[0].j);
    // fmt::print("moves1 {} {}\n", G1.node_list[1].out_edges[0].i, G1.node_list[1].out_edges[0].j);
    // fmt::print("moves2 {} {}\n", G1.node_list[2].out_edges[0].i, G1.node_list[2].out_edges[0].j);
    // fmt::print("moves3 {} {}\n", G1.node_list[3].out_edges[0].i, G1.node_list[3].out_edges[0].j);

    // G1 moves
    fmt::print("start\n");

    // total_bp_dist = 4;
    for (int d = 1; d <= total_bp_dist; d++) {
        std::vector<merge_path> next_paths;

        for (const auto current_path : all_paths) {
            const auto& G1_node = G1.node_list[current_path.current_G1_node];
            const auto& G2_node = G2.node_list[current_path.current_G2_node];
            // fill up next_paths - either move 1 step on G1 or G2
            available_edges(G1_node, G2_node, 1, next_paths, current_path, fc, max_en);
            available_edges(G2_node, G1_node, 2, next_paths, current_path, fc, max_en);
        }
        // fmt::print("size before sort: {} \n", next_paths.size());

        if (next_paths.size() == 0) {
            std::vector<merge_path> empty_paths;
            return empty_paths;
        }

        std::qsort(&next_paths[0], next_paths.size(), sizeof(merge_path), compare_ptable);

        // for (const auto& next_path : next_paths) {
        //     for (const auto& move : next_path.moves) { fmt::print("({} / {}) ", move.i, move.j);
        //     } fmt::print(" d={}, sorted: {} {}\n", d, next_path.current_en, next_path.current_s);
        // }

        // compare current with next pairing table...

        // optimize findpath: this sorting / try moves without ptable creation
        // this is probably better than in findpath? does not guarantee dropping...

        int len = next_paths[0].current_ptable[0];

        if (d != total_bp_dist) {
            for (int u = 0, c = 1; c < next_paths.size(); c++) {
                if (memcmp(next_paths[u].current_ptable, next_paths[c].current_ptable,
                           sizeof(short) * len) != 0) {
                    u++;
                    // next_paths[u] = next_paths[c];
                } else {
                    // next_paths[c].current_s = INT_MAX;
                    // next_paths[c].current_en = INT_MAX;
                    next_paths[c].current_s  = 999;
                    next_paths[c].current_en = 999;
                    // fmt::print("drop {}, {} \n", u, c);

                    // free_intermediate(next + c);
                }
            }
        }

        std::qsort(&next_paths[0], next_paths.size(), sizeof(merge_path), compare_energy);

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

    return all_paths;
}

s_graph findpath::findpath_sections(section current_sections, short* pt_1, short* pt_2)
{

    bool Verbose = false;
    int search_width_multiplier = 20;
    // int search_width_multiplier = 2;

    if (current_sections.nested_sections.size() == 0) {
        // top of the recursion tree - findpath for a single section with start and end

        const int search_width = current_sections.bp_dist*search_width_multiplier;
        // const int search_width = 200;

        // auto result = custom_findpath_method(fc, outer_pt_1, outer_pt_2, search_width, INT_MAX -
        // 1);
        auto result = custom_findpath_method(fc, vrna_ptable_copy(pt_1), vrna_ptable_copy(pt_2),
                                             search_width, INT_MAX - 1, false);

        auto result2 = custom_findpath_method(fc, vrna_ptable_copy(pt_1), vrna_ptable_copy(pt_2),
                                              search_width, INT_MAX - 1, true);
        std::move(result2.begin(), result2.end(), std::back_inserter(result));
        // best path to [0] (lowest max_en)
        sort(result.begin(), result.end(),
             [](const auto& a, const auto& b) -> bool { return a.max_en < b.max_en; });

        s_graph G_inner;
        G_inner.fc      = fc;
        G_inner.pt_1    = pt_1;
        G_inner.pt_2    = pt_2;
        G_inner.bp_dist = current_sections.bp_dist;
        G_inner.add_paths(result);

        // G_inner.info();
        if (Verbose) print_moves(result[0], fc, vrna_db_from_ptable(pt_1), true);

        return G_inner;
    }

    // else: generate the outer pairing table and substract all inner sections, where we will
    // recurse into

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

    // auto result = custom_findpath_method(fc, outer_pt_1, outer_pt_2, search_width, INT_MAX - 1);
    auto result =
        custom_findpath_method(fc, vrna_ptable_copy(outer_pt_1), vrna_ptable_copy(outer_pt_2),
                               search_width, INT_MAX - 1, false);
    auto result2 =
        custom_findpath_method(fc, vrna_ptable_copy(outer_pt_1), vrna_ptable_copy(outer_pt_2),
                               search_width, INT_MAX - 1, true);
    std::move(result2.begin(), result2.end(), std::back_inserter(result));
    // best path to [0] (lowest max_en)
    sort(result.begin(), result.end(),
         [](const auto& a, const auto& b) -> bool { return a.max_en < b.max_en; });

    s_graph G_outer;
    G_outer.fc      = fc;
    G_outer.pt_1    = outer_pt_1;
    G_outer.pt_2    = outer_pt_2;
    G_outer.bp_dist = outer_bp_dist;
    G_outer.add_paths(result);

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
        s_graph G_inner = findpath_sections(nested_section, inner_pt_1, inner_pt_2);

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

        // fmt::print("starting merge:\n{} {} {}\n", vrna_db_from_ptable(merged_pt_1), G_outer.bp_dist,
                //    s1_en);
        // fmt::print("{} {} {}\n", vrna_db_from_ptable(merged_pt_2), G_inner.bp_dist, s2_en);

        int total_bp_dist = G_inner.bp_dist + G_outer.bp_dist;

        // start here the merge supervisor gradually increasing search width

        int max_en = INT_MAX;
        // int  max_en             = -1390;
        int  merge_search_width = 100;
        bool direction          = true;

        auto all_paths = merge_once(G_outer, G_inner, merged_pt_1, s1_en, merged_pt_2, s2_en,
                                    direction, total_bp_dist, max_en, merge_search_width, fc);

        fmt::print("merged: {}\n", all_paths.size());

        outer_bp_dist += G_inner.bp_dist;
        outer_pt_1 = merged_pt_1;
        outer_pt_2 = merged_pt_2;

        s_graph G_new;
        G_outer         = G_new;
        G_outer.fc      = fc;
        G_outer.pt_1    = merged_pt_1;
        G_outer.pt_2    = merged_pt_2;
        G_outer.bp_dist = outer_bp_dist;
        G_outer.add_paths(all_paths);

        // for (const auto p : all_paths) {
        //     print_moves(p, fc, s1, false);
        // }

        if (Verbose) print_moves(all_paths[0], fc, vrna_db_from_ptable(merged_pt_1), true);

        print_moves(all_paths[0], fc, vrna_db_from_ptable(merged_pt_1), false);

        // break;
    }

    

    // IC(current_sections.nested_sections[0].start, current_sections.nested_sections[0].end,
    //    current_sections.nested_sections[0].bp_dist);
    // IC(outer_bp_dist);

    return G_outer;
}

