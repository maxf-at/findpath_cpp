

struct s_edge {
    int i;
    int j;
    int destination;  // references id from s_node
    int en;
};

struct s_node {
    // this is the node object for graph traversal
    int                 id;
    size_t s_hash;
    

    std::vector<s_edge> out_edges = {};
    std::vector<s_edge> in_edges  = {};
    s_node(int id) : id{id} {};  // construct with id (node_count + 1)

    // ~s_node() {
    //      }
};

class s_graph
{
   private:
    void add_edge(int node1, int node2);
    void add_paths(const auto& paths);  // sorted_paths or merged_paths template function
    void add_node(short* p_table, int i, int j, int en);

    // this dictionary is used to populate the node_list vector
    // during the graph initialization, then it is no longer used.
    std::unordered_map<std::string, int> node_map;

    int  node_count      = 0;
    int  last_node       = 0;
    bool last_existing   = false;
    int  current_bp_dist = 0;

   public:
    // node list for graph traversal (only this is used for the merge algorithm)
    std::vector<s_node> node_list = {};

    vrna_fold_compound_t* fc;
    short*                pt_1;
    short*                pt_2;
    int                   bp_dist;
    int                   max_en = -INT_MAX / 2;  // referencing the best path

    int paths_count = 0; // how many individual paths were added to this graph

    void reset();
    void info();
    void display_path(bool result_only = false);

    // s_graph(); //Constructor of the class
    s_graph(vrna_fold_compound_t* fc, short* pt_1, short* pt_2, int bp_dist, const auto& paths);
    s_graph(vrna_fold_compound_t* fc, short* pt_1, short* pt_2, int bp_dist, const auto& paths, int max_en);
    s_graph();  // dummy init for empty graph


    ~s_graph() {
        // free (pt_1);
        // free (pt_2);
    }

};

// Constructor
s_graph::s_graph(vrna_fold_compound_t* fc, short* pt_1, short* pt_2, int bp_dist, const auto& paths)
        : fc{fc}, pt_1{pt_1}, pt_2{pt_2}, bp_dist{bp_dist}
{
    // std::cout << "Graph object being created\n";
    add_paths(paths);
}

// Constructor
s_graph::s_graph(vrna_fold_compound_t* fc, short* pt_1, short* pt_2, int bp_dist, const auto& paths, int max_en)
        : fc{fc}, pt_1{pt_1}, pt_2{pt_2}, bp_dist{bp_dist}, max_en{max_en}
{
    // std::cout << "Graph object being created\n";
    add_paths(paths);
}


s_graph::s_graph()
{
    // dummy init for empty graph
}

void s_graph::reset()
{
    // unused
    fc      = nullptr;
    pt_1    = nullptr;
    pt_2    = nullptr;
    bp_dist = 0;
    node_list.clear();
    node_map.clear();
}

void s_graph::add_node(short* p_table, int i, int j, int en)
{
    // generating a unique hashable string, consisting of the structure
    // and the current basepair distance concatenated
    std::string structure = vrna_db_from_ptable(p_table);
    structure += std::to_string(current_bp_dist);

    

    int  current_node;
    bool current_existing;

    // check if we can reuse an existing node in the graph for the current structure
    if (node_map.contains(structure)) {
        current_existing = true;
        current_node     = node_map[structure];
        // std::cout << "avail. node" << structure << " : " << current_node << endl;

    } else {
        // std::cout << "adding node" << structure << " : " << node_count << endl;

        // add a new node for the current structure
        current_existing    = false;
        current_node        = node_count;
        node_map[structure] = node_count;

        // this creates a new node object (s_node struct) initialized with id=node_count
        node_list.push_back(node_count);
        node_count++;
    }

    node_list[current_node].s_hash = std::hash<std::string_view>()(std::string_view(structure));


    current_bp_dist++;

    // if there's nothing to do...
    if (last_node == current_node or current_node == 0 or (current_existing and last_existing))
    // std::cout << "n edge " << last_node << " : " << current_node << endl;
    {
        last_node     = current_node;  // for case current_node == 0
        last_existing = current_existing;
        return;
    }
   

    // std::cout << "edge " << last_node << " : " << current_node << endl;

    // connect our new node with the rest of the graph
    node_list[last_node].out_edges.push_back({i, j, current_node, en});   // fwd edges
    node_list[current_node].in_edges.push_back({-i, -j, last_node, en});  // bwd edges

    last_node     = current_node;
    last_existing = current_existing;
}

void s_graph::add_paths(const auto& paths)
{
    // add multiple paths to the graph

    for (auto const& path : paths) {
        current_bp_dist            = 0;
        paths_count++;

        auto const& moves          = path.moves;
        short*      current_ptable = vrna_ptable_copy(pt_1);
        
        // fmt::print("S: {} / ({} {})\n", vrna_db_from_ptable(current_ptable), 0,0);

        add_node(current_ptable, 0, 0, 0);  // init graph at S1

        for (auto const& move : moves) {

            int en = vrna_eval_move_pt(fc, current_ptable, move.i, move.j);
            
            if (move.j < 0) {  // delete a basepair
                current_ptable[-move.i] = 0;
                current_ptable[-move.j] = 0;
            } else {  // add a basepair
                current_ptable[move.i] = move.j;
                current_ptable[move.j] = move.i;
            }

            // fmt::print("S: {} / ({} {})\n", vrna_db_from_ptable(current_ptable), move.i, move.j);

            add_node(current_ptable, move.i, move.j, en);
        }

        // path_number++;
    }
}

void s_graph::info()
{
    // for debug:
    // display all nodes / edges present in the graph
    for (auto const& node : node_list) {
        for (auto const& edge : node.out_edges) {
            std::cout << "node: " << node.id;  // << std::endl;
            std::cout << " ( " << edge.i << "  / " << edge.j << " )  dest: " << edge.destination
                      << std::endl;
        }
        if (node.out_edges.size() == 0)
            std::cout << "node: " << node.id << " final node\n";  // << std::endl;
    }

    auto current_ptable = vrna_ptable_copy(pt_1);  // short*

    std::string s  = vrna_db_from_ptable(current_ptable);
    float       en = vrna_eval_structure_pt(fc, current_ptable) / 100.0;
}

void s_graph::display_path(bool result_only)
{
    // prints the best path
    auto current_ptable = vrna_ptable_copy(pt_1);

    std::string s  = vrna_db_from_ptable(current_ptable);
    float       en = vrna_eval_structure_pt(fc, current_ptable) / 100.0;

    if (not result_only) { fmt::print("{} {:7.2f} ({:4}/{:4})\n", s, en, 0, 0); }

    // graph traversal
    int    iter         = 0;
    s_node current_node = node_list[iter];

    float max_en = -9999;

    while (current_node.out_edges.size() != 0) {
        // the first edge always leads to the best result since we added them sorted to the graph
        const auto& m = current_node.out_edges[0];

        if (m.j < 0) {
            current_ptable[-m.i] = 0;
            current_ptable[-m.j] = 0;
        } else {
            current_ptable[m.i] = m.j;
            current_ptable[m.j] = m.i;
        }
        s  = vrna_db_from_ptable(current_ptable);
        en = vrna_eval_structure_pt(fc, current_ptable) / 100.0;
        if (en > max_en) { max_en = en; }

        if (not result_only) { fmt::print("{} {:7.2f} ({:4}/{:4})\n", s, en, m.i, m.j); }

        // fmt::print("ij: {} / {} \n", m.i, m.j);
        current_node = node_list[m.destination];
    }
    // fmt::print("S:{:7.2f} best path in G\n", max_en);
    fmt::print("{:7.2f}\n", max_en);

    free(current_ptable);
}

// unused below

// graph definitions

// struct s_move {
//     int i; /* i,j>0 insert; i,j<0 delete */
//     int j;
// };

// specialized hash function for unordered_map keys
// struct hash_fn {
//     std::size_t operator()(const sorted_move& m) const
//     {
//         std::size_t h1 = std::hash<int>()(m.i);
//         std::size_t h2 = std::hash<int>()(m.j);

//         return h1 ^ h2;
//     }
// };

// https://stackoverflow.com/questions/42701688/using-an-unordered-map-with-arrays-as-keys
// struct array_hash {
//     std::size_t operator()(const short* a) const {
//         std::size_t h = 0;

//         // for (auto e : a) {
//         // }

//         auto e = a[0];
//         h ^= std::hash<int>{}(e)  + 0x9e3779b9 + (h << 6) + (h >> 2);

//         IC("h", a[0], a[1], a[2], h);

//         return h;
//     }
// };

// auto array_hasher(const short* p_table) -> int
// {
//     // std::size_t h = 0;
//     int iterations = p_table[0];
//     int h          = 0;

//     for (int i = 1; i < iterations + 1; i++) {
//         h ^= std::hash<int>{}(p_table[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
//     }
//     // IC("h", p_table[0], sizeof(p_table), p_table[1], p_table[2], h);
//     return h;
// }

// auto hash1  = [](const sorted_move& m) { return m.i + 10 * m.j; };
// auto equal1 = [](const sorted_move& m1, const sorted_move& m2) {
//     return m1.i == m2.i && m1.j == m2.j;
// };
// using move_set = std::unordered_set<sorted_move, decltype(hash1), decltype(equal1)>;
// //
// https://stackoverflow.com/questions/50888127/how-can-i-use-an-unordered-set-with-a-custom-struct
