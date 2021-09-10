#ifndef TOPO_GRAPH_H
#define TOPO_GRAPH_H

#include <assert.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <vector>
#include <set>
#include <map>
#include <tuple>

#include <utility>      // std::pair

#include <algorithm>
#include <iterator>

#include "topo_topo.h"
#include "topo_extinfo.h"

// ---------------------------------------------------------
// Graph based objects.
// ---------------------------------------------------------
// Adapted by keeping a similar interface to the boost graph library.
// A set of graph related objects. 
// These are to be used in searching certain relations such as fundamental 
// cycles given a graph.
// It is not the most secure (vertices and edge descriptions can be changed 
// easily) or efficient code. Std map and set were used for  
// vertex adjacency lists and edge lists.

class vd_graph;
class cell_graph;
class path_adder;
class cell_path;

typedef int vd_Vertex;

// A cell path is a nonintersecting cell sequence connecting two 3cells, not
// already connected by a 2cell. The first element denotes the 3cells, and the
// second element contains the list of paths connecting them. The second 
// element of the path list contains indices of the nonintersecting paths.
typedef std::pair<int,int > cell_cp;

void dummy_graph_stop();

class vd_Edge_desc {
  private:
    int U;
    int V;

  public:
    bool operator<(const vd_Edge_desc& e2)const;
    bool operator==(const vd_Edge_desc& e2)const;
    int source () const;
    int target () const;

    // Default constructor:
    vd_Edge_desc();
    vd_Edge_desc(int U_in, int V_in);

    // Copy constructor:
    vd_Edge_desc(const vd_Edge_desc &obj);
    vd_Edge_desc& operator=( const vd_Edge_desc& obj );
    friend std::ostream& operator<<(std::ostream& s, const vd_Edge_desc& e);
};

class vd_Edge {
  private:
  public:
    vd_Edge_desc first;
    bool second;
    // Default constructor:
    vd_Edge() : first(-1,-1), second(false) {}

    vd_Edge(int U, int V, bool state = false) : first(U, V), second(state) {}

    // Copy constructor:
    vd_Edge(const std::pair<vd_Edge_desc, bool>& obj) : 
                                first(obj.first.source(), obj.first.target()),
                                  second(obj.second) {}
    vd_Edge(const vd_Edge &obj) : 
                                first(obj.first.source(), obj.first.target()),
                                  second(obj.second) {}
    vd_Edge& operator=( const vd_Edge& obj ) {
      first = obj.first;
      second = obj.second;
      return *this;
    }
    vd_Edge& operator=( const std::pair<vd_Edge_desc, bool>& obj ) {
      first = obj.first;
      second = obj.second;
      return *this;
    }
    bool operator==(const vd_Edge& e2)const {
      return first == e2.first;
    }
    ~vd_Edge() {}
};

// Check if two edges have at least a joint vertex.
bool chk_edge_joint(vd_Edge &e1, vd_Edge &e2);

// Assuming the edges are joint (share a vertex), join the edges. Joining an
// edge with itself creates an edge starting and ending at the target vertex of
// the edge.
vd_Edge join_edges(vd_Edge &e1, vd_Edge &e2);

bool chk_discon_circ(std::vector<vd_Edge> temp);

void sort_edge (std::vector<vd_Edge>* edge_list, 
                                          int left = -1, int right = -1);
int partition_edge_list(std::vector<vd_Edge>* edge_list, int left, int right);

class s_graph {
  private:
  public:
    std::vector<int > first;
    std::vector<int > second;

    // Default constructor:
    s_graph() : first(0), second(0) {}

    // Copy constructor:
    s_graph(const s_graph &obj) : first(obj.first),
                                  second(obj.second) {}
    s_graph& operator=( const s_graph& obj ) {
      first = obj.first;
      second = obj.second;
      return *this;
    }
    ~s_graph() {}
};

class circuit {
  private:
  public:
    // first.first: Local 3cell ids
    // first.second: Local 3cell ids
    // second.first: Local 3cell and 2cell ids for the first disjoint graph
    // second.first: Local 3cell and 2cell ids for the second disjoint graph
    s_graph first;
    s_graph second;
    // Default constructor:
    circuit() : first(), second() {}

    // Copy constructor:
    circuit(const circuit &obj) : first(obj.first),
                                  second(obj.second) {}
    circuit& operator=( const circuit& obj ) {
      first = obj.first;
      second = obj.second;
      return *this;
    }
    ~circuit() {}
};

class crc {
  private:
  public:
    std::vector<vd_Edge> first;
    s_graph second;
    // Default constructor:
    crc() : first(), second() {}

    // Copy constructor:
    crc(const crc &obj) : first(obj.first),
                                  second(obj.second) {}
    crc& operator=( const crc& obj ) {
      first = obj.first;
      second = obj.second;
      return *this;
    }
    ~crc() {}
};

// Used in 2cell insertion. Using the nonintersecting list, create a list of 
// ngon_path to store all possible combinations non-intersecting paths.
class ngon_path {
  private:
    cell_cp first;
    std::vector< std::vector< int > > second;
  public:
    // Default constructor:
    ngon_path() : first(std::make_pair(-1,-1)), 
                  second(0, std::vector< int >(0)) {}

    // Copy constructor:
    ngon_path(const ngon_path &obj) : first(obj.first) {
      for(int i = 0; i < second.size(); i++) {
        second.at(i).clear();
      }
      second.resize(obj.second.size());
      for(int i = 0; i < second.size(); i++) {
        second.at(i) = obj.second.at(i);
      }
    }
    ngon_path& operator=( const ngon_path& obj ) {
      first = obj.first;
      for(int i = 0; i < second.size(); i++) {
        second.at(i).clear();
      }
      second.resize(obj.second.size());
      for(int i = 0; i < second.size(); i++) {
        second.at(i) = obj.second.at(i);
      }
      return *this;
    }

    ~ngon_path() {
      for(int i = 0; i < second.size(); i++) {
        second.at(i).clear();
      }
      second.clear();
    }

};

// Used in 2cell insertion. Seperate the path 2cell and 3cell lists and express 
// them in gmi notation. The first list of ngon_gmi stores 2cells and 3cells of 
// each path seperately. Second list stores the ngons. 
// ngon_gmi.cp: cell_cp
// ngon_gmi.cells[i].first: ith path 3cell list
// ngon_gmi.cells[i].second: ith path 2cell list
// ngon_gmi.ngons[j]: jth ngon path list
class ngon_gmi{
  public:
    cell_cp cp;
    // The 2- and 3cells of each path.
    std::vector< std::pair< std::vector< int >, std::vector< int > > > cells;
    // The ngon list, for each ngon possible.
    std::vector< std::vector< int > > ngons;

    // Flag to signify a boundary 0cell.
    bool bd;
    // Holds the index of the current ngon.
    int curr;

    // Default constructor
    ngon_gmi();
    // Copy constructor:
    ngon_gmi(const ngon_gmi &obj);
    ngon_gmi& operator=( const ngon_gmi& obj );

    void clear();

    // Shift cells by 1 or -1, depending on gmi or index notation.
    void shift_cells(int shift);

    // Destructor
    ~ngon_gmi();

    // Print the contents of the ngon_gmi
    void print();

};


class cell_path {
  public:
    std::pair<int,int> first;
    // First component of s_graph is the path nodes, second is the non 
    // intersecting paths. 
    std::vector<s_graph> second;

    ~cell_path();
    cell_path();

    // Copy constructor:
    cell_path(const cell_path &obj);
    cell_path& operator=( const cell_path& obj );
    friend std::ostream& operator<<(std::ostream& s, const cell_path& e);

};


class vd_graph {
  private:
    // reserve_sz is the total expected size of vertices.
    // For faster filling of containers.
    int reserve_sz;

    bool directed;
    std::set<int> v_id;
    //std::map<int, int> index;
    // Adjacency lists.
    // In general, this should allow for repeated elements. 
    std::vector<std::set<int> > adj_vec;
    std::map<int, int > adj_map;
    std::multiset<vd_Edge_desc > edge_list;

    // No copy constructor
    vd_graph(const vd_graph &obj);

  public:
    // Returns the iterators to the vertex range of the graph vertices.
    std::pair<std::set<int>::iterator, 
              std::set<int>::iterator > vertices(); 

    // Returns the iterators to the vertex range of the vertices adjacent to V.
    std::pair<std::set<int>::iterator, 
              std::set<int>::iterator > adjacent_vertices(int V); 

    // Returns the iterators to the edge range.
    std::pair<std::multiset<vd_Edge_desc >::iterator, 
              std::multiset<vd_Edge_desc >::iterator > edge_it(); 

    // Source and target of given edge.
    int source(vd_Edge* e) const;
    int target(vd_Edge* e) const;

    // Add-remove, check vertex:
    void add_vertex(int id);
    bool check_vertex(int id);
    void remove_vertex(int id);

    // Add-remove edge.
    bool remove_edge(int v1, int v2);
    bool remove_edge_dir(int v1, int v2);

    vd_Edge add_edge(int v1, int v2);
    vd_Edge add_edge_dir(int v1, int v2);

    bool check_edge(int v1, int v2);
    // Get edge. Non existing edges will have boolean part false.
    vd_Edge get_edge(int v1, int v2);
    vd_Edge get_edge_dir(int v1, int v2);

    // Assuming given node has exactly two connections, given a node that is 
    // connected return the other connection and assert crash otherwise.
    int get_node_other(int n2, int n_end);

    // Get the number of vertices, edges, graph type (undirected, directed).
    int num_vertices();
    int num_edges();
    bool dir_type();

    // Change the graph type. Resets all connections.
    void set_dir(bool dir);

    // Default constructor.
    vd_graph();
    // Copy constructor.
    vd_graph& operator=(  vd_graph& obj );
    //vd_graph(const vd_graph &obj);
    // Clear the content of the graph.
    void reserve_v(int sz);
    void clear();

    // Test function to print the current graph.
    void print_adj();
    void print_graph(const char* dotname = "output/graph.dot", 
                                          std::vector<vd_Edge>* hl = NULL);
    void print_graph_p(std::map<vd_Edge_desc, bool>& hl_map, std::map<int, std::string>& color_map, std::map<vd_Edge_desc, int>& path_map, const char* dotname = "output/graph.dot");
    void print_graph(std::vector<std::string>* lbls, 
                          const char* dotname);

};

// ---------------------------------------------------------
// Paton elementary circuit finder and Gibbs algorithm.
// ---------------------------------------------------------

// Sorting and comparison functions.
bool edge_comp (const vd_Edge& e1, const vd_Edge& e2);
// Assume sorted. If all edges are the same, they will be sorted the same.
// The comparison should hold.
bool set_edge_comp (std::vector<vd_Edge> e1, std::vector<vd_Edge> e2);

// Paton, K.
// Communications of the ACM. Vol 12/9,1969 September
// Given a set of vertices V and set of edges E on a graph,
// To find a set of fundamental cycles in the graph construct a tree graph in the
// following fashion:
// Let T be the set of vertices already on the tree, initially empty. 
// Let X be set of vertices not on the tree.
// 0. T = 0, X = V
// 1. Take a vertex v from X, put it on the tree and in T, s.t. T = [v], X = V
// 2. Z = (T intersect X).  
// 3. Consider each edge z-w in E for each z in Z.
//    - If w in T, find the unique path between z-w and join the edge to form a 
//      fundamental cycle.
//    - If w not in T, add z-w to the tree, and w to T. Remove z-w from E.
// 4. When all edges of z on E are exhausted, remove z from X. 
class PatonFinder {

  private:
    // Pointer to the original graph. During the course of the algorithm g_init
    // should be valid.
    vd_graph* g_init;
    // Copy of *g_init.
    vd_graph g;
    // Paton algorithm related objects.
    // The graph storing the tree in Paton algorithm.
    vd_graph tree;

    // Vertex containers for tree search:
    std::set<int> T;
    std::set<int> V;
    std::set<int>::iterator T_it;

    // The vertex list of elementary circuits. The ordering represents the 
    // circuit ordering.
    std::vector< std::vector<int> > paths;

    // Path, expressed as cell edge connection :
    // The edge vector of elementary circuits. There is an associated sorting. 
    std::vector< std::vector<vd_Edge> > paths_e;

    // Path expressed by its cells:
    std::vector< std::vector<int> > paths_eset;

    // All unique edge circuits obtained by Gibbs algorithm.
    std::vector< std::vector<vd_Edge> > S;
    // The associated subgraphs.

    // A container for the current circuit disjoint graphs that are to be grouped  
    // into two disjoint graphs using the comp_conn function.
    std::vector< std::vector <int > > disj_frag;
    // A container for circuit edge lists paired with disjoint graphs.
    std::vector< std::pair< std::vector<vd_Edge>, s_graph > > circ_tup;

    // Keeps the indices of Vertices. TODO I'm not sure how necessary.
    std::map<int, int> index;

    // Iterators.
    std::pair<std::set<int>::const_iterator,
                std::set<int>::const_iterator> vp;

    std::set<int>::const_iterator ai;
    std::set<int>::const_iterator ai_end;

    std::set<int>::const_iterator ai_int;
    std::set<int>::const_iterator ai_end_int;

    // c_graph is referred to when collecting disjoint graphs into two groups. 
    // If it is NULL, the two disjoint graphs are not populated.
    cell_graph* c_graph;

    // Copy constructor and assignment private.
    PatonFinder& operator=(  PatonFinder& obj );
    PatonFinder(const PatonFinder& obj );
  public:

    // Sorting and comparison functions.

    //bool edge_comp (const vd_Edge& e1, const vd_Edge& e2);

    // Assume sorted. If all edges are the same, they will be sorted the same.
    // The comparison should hold.
    //bool set_edge_comp (std::vector<vd_Edge> e1, std::vector<vd_Edge> e2);

    //void sort_edge (std::vector<vd_Edge>* edge_list);
    int partition_edge_set(std::vector<std::vector<vd_Edge> >* e_set, 
                                                      int left, int right);

    void sort_edge_sets (std::vector<std::vector<vd_Edge> >* e_set, 
                                                  int left = -1, int right = -1);
    //void sort_edge_sets (std::vector<std::vector<vd_Edge> >* e_set);

    void isect_edge (std::vector<vd_Edge>* e1, std::vector<vd_Edge>* e2, 
                     std::vector<vd_Edge>* out);

    void diff_edge (std::vector<vd_Edge>* e1, std::vector<vd_Edge>* e2, 
                     std::vector<vd_Edge>* out);

    void merge_edge (std::vector<vd_Edge>* e1, std::vector<vd_Edge>* e2, 
                     std::vector<vd_Edge>* out);

    void merge_edge_set (std::vector<std::vector<vd_Edge> >* e1, 
                         std::vector<std::vector<vd_Edge> >* e2, 
                     std::vector<std::vector<vd_Edge> >* out);

    // Default constructor.
    PatonFinder();
    PatonFinder(vd_graph* g_in);
    void reload_graph(vd_graph* g_in);

    // PATON RELATED:
    // Elementary circuit detection, recursive one does not work as the 
    // iterator cannot be used in multiple recursions.
    bool find_nonrecursive(int v1, int v2);

    // Using the Paton algorithm, find all elementary circuits.
    void find_paton();

    // GIBBS RELATED:

    // Given two vector list of cells and a function that tests connectivity of 
    // cells in these lists, return true if they have at least a connection.
    // Mainly used in collecting the disjoint graphs obtained by removing the 
    // circuit into at most two groups.
    bool chk_conn_list(std::vector<int>* c1_list,
      std::vector<int>* c2_list);

    // Create a subgraph for each path.
    void create_sub(cell_graph* const c_graph);

    // Using the Gibbs algorithm, find all unique circuits.
    void find_gibbs();

    // Convert the elementary circuit lists into edge circuits lists.
    // These are used in finding circuits by Gibbs algorithm.
    void vert_2_edge();

    // GET THE CIRCUITS:
    // Given a graph, first find the elementary circuits by Paton algorithm.
    // Using the Gibbs algorithm, find all circuits as edge lists. 
    // Convert the edge lists into vertex lists.
    void find_circuits();

    // Return a pointer to the circuits.
    std::vector< std::pair< std::vector<vd_Edge>, s_graph > >* get_circuits();

    // Clear the contents of the containers.
    void clear();

    // Default destructor.
    ~PatonFinder();
};

// ---------------------------------------------------------
// Functions for doing graph based search on a cellbase object.
// ---------------------------------------------------------
class cell_graph {
  private:
    struct cell_base* cb;

    // ID of the central 0cell.
    int cell_ctr;
    // 2-3 adjacency graph.
    vd_graph g;
    // 1-2 adjacency graph.
    vd_graph g_12;

    std::map<int, int> cell_type;
    std::vector<int> cell_lbls;
    //std::map <int, std::string > cell_lbls;
    PatonFinder PF;
    //std::vector< std::pair< std::vector<vd_Edge>, s_graph > > circ_tup;
    std::vector< crc > circ_tup;

    // The possible paths and non-intersecting lists.
    std::vector< cell_path > path;
    // Current path index, used to limit path calculations to a single 3cell 
    // couple.
    int path_act;

    // Used in DFS algorithm. By exhaustively going over the graph, find all
    // eligible paths between two 3cells. 
    // Path_curr stores the current trace of cells. 
    // visited keeps track of the current visited cells. 
    // cell1_flag denotes that the current path is around a 1cell and 
    // cell_1act is the 1cell index. 
    // It is possible that a series of 2cells might be joined around more than
    // one 1cell. 
    int cell_start;
    int cell_end;

    std::vector< int > path_curr;
    int path_curr_id;

    int c2_curr;
    int c3_curr;

    std::map< int, bool > visited;
    bool cell1_flag; // 0: 2nd 2cell, check if joint. If not remain 0. Any 
                     // following match is problematic. 
                     // 1: All following 2cell must join at the same 1cell.
    std::set<int> cell1_act;
    std::set<int> c1_int;

    // ---------------------------------------------------------
    // Variables for detecting insertion at corner 0cells.
    // ---------------------------------------------------------
    // Number of exterior 3cells.
    int ext_sz;

    // Flag to store whether exterior is included or not.
    bool calc_corner;
    // Flag to store whether the 0cell is exterior or not.
    bool cell_ext;
    bool cell_ext_corner;
/*
    // 2c exterior group membership and for the active 3c couple, path group.
    // Corner 1cell flag
    std::map< int, bool > ext_1c;
    // Corner 1cell 1st exterior 2cell neigbor
    std::map< int, int > ext_1c1;
    // Corner 1cell 2nd exterior 2cell neigbor
    std::map< int, int > ext_1c2;
    // Ext 2cell flag
    std::map< int, bool > ext_2c;
    // Ext 2cell part tag
    std::map< int, int > ext_2c_part;
    // Path exterior flag
    std::map< int, int > ext_path;

    // Clear the variables.
    void clear_ext_map();
*/
    // ---------------------------------------------------------
    // Call the routines of the PatonFinder.
    // ---------------------------------------------------------
    // Treat the ciruits and paths for graphs containing multiple exterior 
    // cells.
    void treat_circuits();
    void treat_paths_act();

    void treat_graphs();

    // Check circuit list for elements with multiple circuits. There shouldn't
    // be any.
    void chk_discon_circs();
    // Find the circuits within the loaded graph. Get the disjoint graphs 
    // separated by the circuits.
    void find_circuits();
    // Set the non-intersecting lists of all paths.
    void find_ext_nonintsct();
    void find_nonintsct();
    // Find the disconnected 3cells and the accompanying paths within the 
    // loaded graph.
    bool chk_1cell();
    // The elements of the input list of 1cells that are adjacent to the 0cell
    // are kept. The others are removed.
    // bool chk_1cell_adj(struct ent_conn* e_in);

    // Check if the current path has reached the end 3cell. If so, add it to 
    // the path list.
    bool try_path();

    // Appends the current 2- and 3-cells to the current path.
    void path_append();
    // Remove the last 2- and 3-cells from the current path.
    void path_remove();

    void trace_path();
    void find_paths();


    // Copy constructor and assignment private.
    cell_graph(const cell_graph &obj);
    cell_graph& operator=( const cell_graph& obj );

  public:

    struct ent_conn cells3;
    struct ent_conn cells2;
    struct ent_conn cells1;

    int get_0c();
    bool get_ext();

    void find_slice(std::vector< 
            std::pair<std::vector<int >, std::vector<int > > >* path_cells, 
            std::vector< std::pair< std::pair<int,int>, 
                         std::vector<std::vector<int > > > >* slice_cells);

    // ---------------------------------------------------------
    // Check connections of different topologies.
    // ---------------------------------------------------------
    // Convert the id of a 23 adjacency graph 2cell to a 12 adjacency graph.
    int conv_c2_12(int c2);
    // Given two 2cells in 2/3-cell adjacency graph, check if they share a  
    // common 1cell adjacency. Used within disjoint cell grouping inside 
    // PatonFinder.
    bool comp_conn_12(int c21, int c22);

    // ---------------------------------------------------------
    // Load the graphs
    // ---------------------------------------------------------
    // By using a given cell and adjacency list.
    void load_cells(int tag_0cell, std::vector<std::vector<int > >* cells, 
                       std::vector<std::vector<std::vector<int > > >* adj);
    // Get the adjacency maps in the neighborhood of 0cell, in topo indices.
    void get_23adj(int tag_0cell);

    int node_nbr(int n, int i, int j);
    void load_tiles(int n, int m);

    // Communicate with the graph externally. Only for test purposes.
    void remove_edge(int v1, int v2);

    // Print the topology into a dot file or to the terminal.
    void print_adj();
    void print_graph(const char* dotname = "output/graph.dot");

    void print_graph_p(const char* dotname = "output/graph.dot", 
     std::vector<std::pair<std::vector<int >, std::vector<int > > >* path_cells = NULL); 

    void print_graph_lbl(const char* dotname = "output/graph.dot");

    // ---------------------------------------------------------
    // Retrieve information
    // ---------------------------------------------------------
    // Get the node cell type. Each node(on 23 or 12 adjacency map) has a cell
    // type. Return the cell type of the node.
    int g23type(int node_id);
    int g12type(int node_id);
    int g23id(int node_id);
    int g12id(int node_id);

    int gn23(int c_dim, int c_id);
    int gn12(int c_dim, int c_id);
    // ---------------------------------------------------------
    // Return the pointer to the circuits and the disjoint graphs.
    // ---------------------------------------------------------

    // Return a pointer to the circuits.
    //std::vector< std::pair< std::vector<vd_Edge>, s_graph > > * get_circuits();
    std::vector< crc > * get_circuits();
    // Return a pointer to the paths.
    std::vector< cell_path >* get_paths();
    void get_path_ngon_gmi(int path_id, ngon_gmi* ng);
    void get_path_ngon(int path_id, ngon_gmi* ng);

    // Print the circuits and the associated subgraphs.
    void print_circ();

    // ---------------------------------------------------------
    // Circuit classification and reordering.
    // ---------------------------------------------------------

    // After obtaining the circuits using the PatonFinder, sort the lists
    // and group the circuits 

    bool cont_3cell_and (const s_graph* circ_graph);
    bool cont_3cell_or (const s_graph* circ_graph);
    bool cont_23cell (const s_graph* circ_graph);

    bool cont_3cell (const std::vector<int>* disj_graph);

    // Check both disjoint graphs, if both of them are non-empty, return true.
    bool chk_nondegen (const s_graph* circ_graph);
    bool chk_nondegen (const std::vector<int>* disj_graph);

    bool circ_tup_comp (const crc& circ_1, const crc& circ_2);

    // Compare based on disjoint graphs. If they are the same, return true.
    // Used in treating corner 0cell circuits, where same circuit may be
    // repeated.
    bool circ_tup_comp_disj (const crc& circ_1, const crc& circ_2);

    //void quicksort_circ_tup (std::vector< crc >* circ_tup_list);
    void quicksort_circ_tup (int left, int right);
    int partition (int left, int right);
    // Circuit classification based on the cell dimensions on the disjoint 
    // graphs generated by removing the circuit from the cell complex.
    // Store the number of non-degenerate cases in an ordered circ_tup list.
    // 0: Both disjoint graphs contain at least a 3cell.
    // 1: Both disjoint graphs contain at least a 2cell or 3cell.
    // Rest of the circuits: One disjoint graph contains no cell.
    int circ_class[2];
    void get_circ_class(int cell0, std::vector<int>* pos);
    int get_circ_type(int c_in);

    // Return the number of n-cells around the 0cell.
    int get_3cell_sz();
    int get_2cell_sz();
    int get_1cell_sz();

    // ---------------------------------------------------------
    // TODO may be replaced if a better volume preservation scheme is used.
    // ---------------------------------------------------------
    // For the moment, corner 1cells are kept for preserving volume. Set the 
    // flag for trying 2cell insertions between exterior and 3cells connected
    // to the exterior.
    void set_calc_corner(bool fix);

    // ---------------------------------------------------------
    // Reset the cell_graph except for the c_base connection.
    // ---------------------------------------------------------
    void clear();
    void load_cb(struct cell_base* c_base, bool corner_flag = false);

    // ---------------------------------------------------------
    // Constructor.
    // ---------------------------------------------------------
    cell_graph();
    cell_graph(struct cell_base* c_base, bool corner_flag = false);
    //cell_graph();
    ~cell_graph();
};


// ---------------------------------------------------------
// Circuit and path information stored at each 0cell. Cell insertion check.
// ---------------------------------------------------------
class cell_ins_chk {
  private:

    struct cell_base* cb;
    cell_graph cg;

    std::vector<int > circ_class;

    // [0_cell, circuit] first(circuit 3/2cells) second(disjoint graphs) 
    std::vector< std::vector< circuit > > circ_tup;
    // We could also store non intersecting path lists for each path.
    std::vector< std::vector< cell_path > > path;

    // Check for the loading of the cell_base.
    bool load_flag;

    // Check for calculating insertions at corner 0cells differently.
    bool calc_corner;

    // Copy constructor and assignment private.
    cell_ins_chk(const cell_ins_chk &obj);
    cell_ins_chk& operator=( const cell_ins_chk& obj );

  public:
    // TODO a method to call the pointer to these cell vectors.
    std::vector<ent_conn > cells3;
    std::vector<ent_conn > cells2;
    std::vector<ent_conn > cells1;
    // ---------------------------------------------------------
    // Constructor.
    // ---------------------------------------------------------
    cell_ins_chk();
    cell_ins_chk(struct cell_base* c_base, bool corner_flag);

    void update();

    void load_cb(struct cell_base* c_base, bool corner_flag = false);

    // ---------------------------------------------------------
    // Update path and circuit information.
    // ---------------------------------------------------------
    // Given a 0cell, derive the path and circuit information from the cell 
    // graph.
    // ---------------------------------------------------------
    void update_path(int cell0);
    void update_path_gmi(int cell0);

    void clear_path(int cell0);

    // ---------------------------------------------------------
    // Information retrival.
    // ---------------------------------------------------------
    // Return the pointer to the cellbase structure.
    cell_base* get_cbase();

    // Returns the circuit type, 0 for both containing at least a 3cell, 
    // 1 for at least 2cells, 
    // other cases for degenerate.
    int get_circ_type(int cell0, int circ);
    int get_circ_type_gmi(int cell0, int circ);
    // Return number of 0cell circuits.
    int get_circ_sz();
    int get_circ_sz(int cell0);
    int get_circ_sz_gmi(int cell0);

    circuit* get_circ(int cell0, int circ);
    circuit* get_circ_gmi(int cell0, int circ);

    // Given the circuit, return the actual cell indices in the cell_base 
    // structure. 
    void get_circ_topo(ent_conn* c2, ent_conn* c2_circ, ent_conn* c3_circ, 
                      circuit* circ_in,
                      int g_id, int cell_id);

    // Use gmi indices. circ_in is in local graph indices.
    void get_circ_topo_gmi(ent_conn* c2, ent_conn* c2_circ, ent_conn* c3_circ, 
                      circuit* circ_in,
                      int g_id, int cell_id);

    void get_circ_topo_dis(std::vector<ent_conn>* cs, 
                               ent_conn* c2_circ, ent_conn* c3_circ, 
                      circuit* circ_in,
                      int g_id, int cell_id);

    // Use gmi indices. circ_in is in local graph indices.
    void get_circ_topo_dis_gmi(std::vector<ent_conn>* cs,
                               ent_conn* c2_circ, ent_conn* c3_circ, 
                      circuit* circ_in,
                      int g_id, int cell_id);


    // Given the circuit, return the actual 2cell indices of a disjoint graph in
    // the cell_base structure. 
    void get_circ_2c_disj(ent_conn* c2, 
                      circuit* circ_in,
                      int g_id, int cell_id);

    // Using the gmi indices.
    void get_circ_2c_disj_gmi(ent_conn* c2,
                      circuit* circ_in,
                      int g_id, int cell_id);

    // Given the 0cell, return the number of disconnected 3cell couples. 
    int get_path_sz(int cell0);

    // Given the 0cell, return the all 3cells member of a disconnected 3cell 
    // couple. Mainly used in mesh preconditioning step of cell insertion.
    std::vector<int> get_path_3cells(int cell0);

    // Given the 0cell and disconnected 3cell couple id, return a pointer to 
    // the cell_path. 
    cell_path* get_path(int cell0, int path_id);

    // Return the 3cells to be connected, in gmi indices. 
    std::pair<int, int> get_path_3c_gmi(int cell0, int path_id);
    void get_path_ngon_gmi(int cell0, int path_id, ngon_gmi* ng);

    // Given the 0cell(gmi), return all 3cells(cbase) bounded by the 0cell.
    std::vector<int> get_3cells(int cell0);

    // Given the 0cell(gmi), return all 3cells(gmi) bounded by the 0cell.
    std::vector<int> get_3cells_gmi(int cell0);

    // Find the slices given a 0-cell and path.
    void find_slice(int c0, std::vector< 
            std::pair<std::vector<int >, std::vector<int > > >* path_cells, 
            std::vector< std::pair< std::pair<int,int>, 
                         std::vector<std::vector<int > > > >* slice_cells);

    // Find the slices given a 0-cell and path.
    void find_slice(path_adder* pa);
    // Find the slices given a 0-cell and circuit.
    void find_circuit(path_adder* pa); 

    // Replace the indices from an input vector using the graph indices.
    void repl_index(int dim, int c0, std::vector<int>* c_in);

    // For the isotropic case, return the 0cells around which new cells can be 
    // inserted. Quick way of correcting the starting microstructure.
    //std::vector<int> ret_ins();
    std::vector<int> ret_ins_gmi();

    // For isotropic case, check if cell is insertible based on Euler 
    // characteristic.
    bool insertible(int tag_0cell, bool isotropic);

    // ---------------------------------------------------------
    // TODO may be replaced if a better volume preservation scheme is used.
    // ---------------------------------------------------------
    // For the moment, corner 1cells are kept for preserving volume. Set the 
    // flag for trying 2cell insertions between exterior and 3cells connected
    // to the exterior.
    void set_calc_corner(bool fix);

/*
    // Given the 0cell and the disconnected 3cell couple index, return the 
    // pointer to the path.
    void trace_non_int(int cell0, int path_id, std::vector<int>* non_int);

    // Given the 0cell and the disconnected 3cell couple index, return the 
    // pointer to the path.
    void get_path(int cell0, int path_id, std::vector<ngon_path>* ngons);
    void get_path_gmi(int cell0, int path_id, std::vector<ngon_path>* ngons);

    // Given the path, return the actual cell indices in the cell_base 
    // structure. 
    void get_path_topo(std::vector<ent_conn>* c2_path, 
                  std::vector<ent_conn>* c3_path, cell_cp c3_dj, 
                  ngon_path* path_in, int cell_id);

    // Use gmi indices. circ_in is in local graph indices.
    void get_path_topo_gmi(std::vector<ent_conn>* c2_path, 
                  std::vector<ent_conn>* c3_path, cell_cp c3_dj, 
                  ngon_path* path_in, int cell_id);
*/

    // ---------------------------------------------------------
    // Output related.
    // ---------------------------------------------------------
    // Print the circuits and paths of the 0cells.
    void print_circ();
    void print_path();

    void print_circ(int cell);
    void print_path(int i);

    // Print all the 0cell topologies into dot files.
    void print_graph();

    // ---------------------------------------------------------
    // Destructor.
    // ---------------------------------------------------------
    void clear();
    ~cell_ins_chk();
};

// ---------------------------------------------------------
// Ngon tracer. Used in extracting ngons from paths of a given 0cell.
// ---------------------------------------------------------

class ngon_tracer {
  private:
    ngon_gmi* ng;
    cell_ins_chk* c_ins;
    cell_path* c_path;

    // Copy constructor and assignment private.
    ngon_tracer(const ngon_tracer &obj);
    ngon_tracer& operator=( const ngon_tracer& obj );

  public:
    ngon_tracer ();
    ngon_tracer (cell_ins_chk* c_ins_in, ngon_gmi* ng_in);

    void get_path_ngon_gmi(int cell0, int path_id);
    void trace_path(std::vector<int> in_curr, std::vector<int> ngon_curr);

};

// Cell_graph version of ngon_tracer
class ngon_tracer_g {
  private:
    ngon_gmi* ng;
    cell_graph* cg;
    cell_path* c_path;

    // Copy constructor and assignment private.
    ngon_tracer_g(const ngon_tracer_g &obj);
    ngon_tracer_g& operator=( const ngon_tracer_g& obj );

  public:
    ngon_tracer_g (cell_graph* cg_in, ngon_gmi* ng_in);

    void get_path_ngon(int path_id);
    void get_path_ngon_gmi(int path_id);
    void trace_path(std::vector<int> in_curr, std::vector<int> ngon_curr);

};

// This object is fed into the cell_graph to extract the cells 
class path_adder {
  private:
    int c0;
    int p_id;
    int ng_id;

    // Copy constructor and assignment private.
    path_adder(const path_adder &obj);
    path_adder& operator=( const path_adder& obj );

  public:
    cell_adder* ca;
    cell_base* cb;

    path_adder(struct cell_base* c_base, cell_adder* c_add, int cell0);
    void reload(struct cell_base* c_base, cell_adder* c_add, int cell0); 
    void set_path(int path_id, int ngon_id); 
    void set_circ(int circ_id); 

    int get_c0();
    int get_pid();
    int get_ng();

    void clear();
};

#endif
