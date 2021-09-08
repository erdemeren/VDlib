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

#include <algorithm>
#include <iterator>

#include <boost/config.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/subgraph.hpp>

//#include <boost/graph/johnson_all_pairs_shortest.hpp>

#include "topo_topo.h"

// ---------------------------------------------------------
// Commonly used graph types.
// ---------------------------------------------------------

typedef boost::adjacency_list<boost::setS, boost::vecS,
    boost::undirectedS> graph;

// Somehow subgraph requires these additional property options.
typedef boost::subgraph < boost::adjacency_list<boost::vecS, boost::vecS,
    boost::undirectedS, boost::property<boost::vertex_color_t, int>, boost::property<boost::edge_index_t, int> > > s_graph;

typedef boost::graph_traits<graph>::vertex_iterator v_it;
typedef boost::property_map<graph, boost::vertex_index_t>::type IndexMap;

typedef boost::graph_traits<graph>::vertex_descriptor Vertex;
typedef std::pair<graph::edge_descriptor, bool> Edge;

// ---------------------------------------------------------
// Paton elementary circuit finder and Gibbs algorithm.
// ---------------------------------------------------------

class PatonFinder {

  private:
    // Pointer to the original graph. During the course of the algorithm g_init
    // should be valid.
    graph* g_init;
    // Copy of *g_init.
    graph g;
    // Paton algorithm related objects.
    // The graph storing the tree in Paton algorithm.
    graph tree;
    std::set<Vertex> T;
    std::set<Vertex> V;
    std::set<Vertex>::iterator T_it;

    // The vertex list of elementary circuits. The ordering represents the 
    // circuit ordering.
    std::vector< std::vector<Vertex> > paths;
    // The edge vector of elementary circuits. There is an associated sorting. 
    std::vector< std::vector<Edge> > paths_e;
    std::vector< std::vector<Vertex> > paths_eset;

    // All unique edge circuits obtained by Gibbs algorithm.
    std::vector< std::vector<Edge> > S;
    // The associated subgraphs.
    std::vector<std::vector<Vertex > > s_graphs;

    // Keeps the indices of Vertices. TODO I'm not sure how necessary.
    IndexMap index;

    // Iterators.
    std::pair<v_it, v_it> vp;

    boost::graph_traits<graph>::adjacency_iterator ai;
    boost::graph_traits<graph>::adjacency_iterator ai_end;

    boost::graph_traits<graph>::adjacency_iterator ai_int;
    boost::graph_traits<graph>::adjacency_iterator ai_end_int;

  public:

    // Sorting and comparison functions.

    bool edge_comp (Edge e1, Edge e2);

    // Assume sorted. If all edges are the same, they will be sorted the same.
    // The comparison should hold.
    bool set_edge_comp (std::vector<Edge> e1, std::vector<Edge> e2);

    void sort_edge (std::vector<Edge>* edge_list);

    void sort_edge_sets (std::vector<std::vector<Edge> >* e_set);

    void isect_edge (std::vector<Edge>* e1, std::vector<Edge>* e2, 
                     std::vector<Edge>* out);

    void diff_edge (std::vector<Edge>* e1, std::vector<Edge>* e2, 
                     std::vector<Edge>* out);

    void merge_edge (std::vector<Edge>* e1, std::vector<Edge>* e2, 
                     std::vector<Edge>* out);

    void merge_edge_set (std::vector<std::vector<Edge> >* e1, 
                         std::vector<std::vector<Edge> >* e2, 
                     std::vector<std::vector<Edge> >* out);

    // Default constructor.
    PatonFinder(graph* g_in);

    // Elementary circuit detection, recursive one does not work as the 
    // iterator cannot be used in multiple recursions.
    bool find_nonrecursive(Vertex v1, Vertex v2);

    // Create a subgraph for each path.
    void create_sub();

    // Convert the elementary circuit lists into edge circuits lists.
    // These are used in finding circuits by Gibbs algorithm.
    void vert_2_edge();

    // Convert the edge circuits lists found by Gibbs algorithm into vertex
    // lists.
    void edge_2_vert();

    // Using the Gibbs algorithm, find all unique circuits.
    void find_gibbs();

    // Using the Paton algorithm, find all elementary circuits.
    void find_paton();

    // Given a graph, first find the elementary circuits by Paton algorithm.
    // Using the Gibbs algorithm, find all circuits as edge lists. 
    // Convert the edge lists into vertex lists.
    void find_circuits();

    // Return a pointer to the circuits.
    const std::vector< std::vector<Edge> >* get_circuits();
    // Return a pointer to the associated subgraphs.
    const std::vector< std::vector<Vertex> >* get_disjoint();


    // Default destructor.
    ~PatonFinder() {};
};

// ---------------------------------------------------------
// Functions for doing graph based search on a cellbase object.
// ---------------------------------------------------------
struct cell_graph {
  private:
    struct cell_base* cb;
    struct ent_conn cells3;
    struct ent_conn cells2;
    graph g;
    std::vector<int> cell_type;
    std::vector<std::string> cell_lbls;
    //std::map <int, std::string > cell_lbls;
    PatonFinder* PF;
    const std::vector< std::vector<Edge> >* S;
    const std::vector< std::vector<Vertex> >* s_graphs;

    // Default constructor not allowed.
    cell_graph();

    // ---------------------------------------------------------
    // Call the routines of the PatonFinder.
    // ---------------------------------------------------------

    // Find the circuits within the loaded graph. Get the disjoint graphs 
    // separated by the circuits.
    void find_circuits();

  public:

    // ---------------------------------------------------------
    // Load the graph different topologies.
    // ---------------------------------------------------------
    void get_23adj(int tag_0cell);
    void print_graph();
    // ---------------------------------------------------------
    // Constructor.
    // ---------------------------------------------------------
    cell_graph(struct cell_base* c_base);
    ~cell_graph();
};

#endif
