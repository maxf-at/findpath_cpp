

// graph definitions

// struct s_move {
//     int i; /* i,j>0 insert; i,j<0 delete */
//     int j;
// };

// specialized hash function for unordered_map keys
struct hash_fn {
    std::size_t operator()(const sorted_move& m) const
    {
        std::size_t h1 = std::hash<int>()(m.i);
        std::size_t h2 = std::hash<int>()(m.j);

        return h1 ^ h2;
    }
};

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

auto array_hasher(const short* p_table) -> int
{
    // std::size_t h = 0;
    int iterations = p_table[0];
    int h          = 0;

    for (int i = 1; i < iterations + 1; i++) {
        h ^= std::hash<int>{}(p_table[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    // IC("h", p_table[0], sizeof(p_table), p_table[1], p_table[2], h);
    return h;
}

auto hash1  = [](const sorted_move& m) { return m.i + 10 * m.j; };
auto equal1 = [](const sorted_move& m1, const sorted_move& m2) {
    return m1.i == m2.i && m1.j == m2.j;
};
using move_set = std::unordered_set<sorted_move, decltype(hash1), decltype(equal1)>;
// https://stackoverflow.com/questions/50888127/how-can-i-use-an-unordered-set-with-a-custom-struct

struct s_node;

struct s_edge {
    // Data data;
    // sorted_move move;
    int i;
    int j;
    int destination;  // references s_node::id
};

struct s_node {
    // this is the node object for graph traversal
    int                 id;
    std::vector<s_edge> out_edges = {};
    std::vector<s_edge> in_edges  = {};
    s_node(int id_s) { id = id_s; };  // construct with id

    // ~s_node() { 
    //     IC("de");
    //     out_edges.clear();
    //     in_edges.clear();        
    //      }
};

class s_graph
{
   private:
    // std::vector<int> edge1Array;
    // std::vector<int> edge2Array;
    void add_edge(int node1, int node2);
    // void deleteEdge(int node1, int node2);

    // this dictionary is used to populate the node_list
    // during the graph initialization, then it is no longer used.
    std::unordered_map<std::string, int> node_map;

    int  node_count    = 0;
    int  last_node     = 0;
    bool last_existing = false;

   public:
    void reset();
    void add_paths(const auto& paths); // sorted_paths or merged_paths template function
    void add_node(short* p_table, int i, int j);
    void info();

    std::vector<s_node> node_list = {};

    int current_bp_dist = 0;

    vrna_fold_compound_t* fc;
    short*                pt_1;
    short*                pt_2;
    int                   bp_dist;

    // void addUEdge(int node1, int node2);
    // int connectedComponents(void);
    // int checkNode(int node);
    // int checkEdge(int node1 , int node2);
    // void deleteNode(int node);
    // void deleteUEdge(int node1, int node2);
    // void printNodeList(void);
    // s_graph(); //Constructor of the class
    s_graph();
    // s_graph(int size); //Constructor of the class
};

// Constructor
s_graph::s_graph()
{
    std::cout << "Graph object being created\n";
    // graphSize = size;
    // cout << "Graph object created"<<endl;
}

void s_graph::reset()
{
    fc      = nullptr;
    pt_1    = nullptr;
    pt_2    = nullptr;
    bp_dist = 0;

    node_list.clear();

    node_map.clear();

    // node_count    = 0;
    // last_node     = 0;
    // last_existing = false;
}

void s_graph::add_node(short* p_table, int i, int j)
{
    // generating a unique hashable string, consisting of the structure
    // and the current basepair distance concatenated
    std::string structure = vrna_db_from_ptable(p_table);
    structure += std::to_string(current_bp_dist);

    int  current_node;
    bool current_existing;

    // bool std::unordered_map::contains( const Key& key ) const t;
    if (node_map.contains(structure)) {
        current_existing = true;
        current_node     = node_map[structure];
        // std::cout << "avail. node" << structure << " : " << current_node << endl;
        current_bp_dist++;
    } else {
        // std::cout << "adding node" << structure << " : " << node_count << endl;
        current_existing = false;
        current_node     = node_count;

        node_map[structure] = node_count;

        // this does not add an integer to the node list,
        // it creates a new node object initialized with id=node_count
        node_list.push_back(node_count);

        node_count++;
        current_bp_dist++;
    }

    if (last_node == current_node or current_node == 0 or (current_existing and last_existing))
    // std::cout << "n edge " << last_node << " : " << current_node << endl;
    {
        last_node     = current_node;  // for case current_node == 0
        last_existing = current_existing;
        return;
    }

    // if (current_existing and last_existing) {
    //     std::cout << "existing link" << endl;
    // }

    // std::cout << "edge " << last_node << " : " << current_node << endl;

    node_list[last_node].out_edges.push_back({i, j, current_node});   // fwd edges
    node_list[current_node].in_edges.push_back({-i, -j, last_node});  // bwd edges

    last_node     = current_node;
    last_existing = current_existing;
    // node_list.push_back(node);
}

void s_graph::add_paths(const auto& paths)
{
    for (auto const& path : paths) {
        current_bp_dist            = 0;
        auto const& moves          = path.moves;
        short*      current_ptable = vrna_ptable_copy(pt_1);

        add_node(current_ptable, 0, 0);  // init

        for (auto const& move : moves) {
            if (move.j < 0) {  // delete a basepair
                current_ptable[-move.i] = 0;
                current_ptable[-move.j] = 0;
            } else {  // add a basepair
                current_ptable[move.i] = move.j;
                current_ptable[move.j] = move.i;
            }
            add_node(current_ptable, move.i, move.j);
        }
    }
}

void s_graph::info()
{
    // nodes debug
    for (auto const& node : node_list) {
        for (auto const& edge : node.out_edges) {
            std::cout << "node: " << node.id;  // << endl;
            std::cout << " ( " << edge.i << "  / " << edge.j << " )  dest: " << edge.destination <<
            endl;
        }
        if (node.out_edges.size()==0)
            std::cout << "node: " << node.id << " final node\n" ;  // << endl;
    }

    auto current_ptable = vrna_ptable_copy(pt_1);  // short*

    std::string s  = vrna_db_from_ptable(current_ptable);
    float       en = vrna_eval_structure_pt(fc, current_ptable) / 100.0;

    // fmt::print("{} {:7.2f} ({:4}/{:4})\n", s, en, 0, 0);

    // // random graph traversal
    // int    iter         = 0;
    // s_node current_node = node_list[iter];  // dont copy this

    // while (current_node.out_edges.size() != 0) {
    //     // if (current_node.out_edges.size() == 0) break;

    //     // fmt::print("nodes: {} \n", current_node.out_edges.size());

    //     int next_node = current_node.out_edges.size() - 1;

    //     const auto& m = current_node.out_edges[next_node];

    //     if (m.j < 0) {
    //         /*it's a delete move */
    //         current_ptable[-m.i] = 0;
    //         current_ptable[-m.j] = 0;
    //     } else {
    //         current_ptable[m.i] = m.j;
    //         current_ptable[m.j] = m.i;
    //     }
    //     s  = vrna_db_from_ptable(current_ptable);
    //     en = vrna_eval_structure_pt(fc, current_ptable) / 100.0;
    //     fmt::print("{} {:7.2f} ({:4}/{:4})\n", s, en, m.i, m.j);

    //     // fmt::print("ij: {} / {} \n", m.i, m.j);
    //     current_node = node_list[m.destination];
    // }



}

