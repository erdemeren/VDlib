#include <assert.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <functional>

#include <vector>
#include <set>
#include <map>
#include <tuple>

#include <utility>      // std::pair

#include <algorithm>
#include <iterator>

#include <cstring>

#include "topo_topo.h"
#include "topo_graph.h"

// Directed graph: Reversed edges are not equivalent.
// Indirected graph: When adding an edge, smallest index is always added 
// first, so comparing cross vertices (e.g. U1 vs V2) is not necessary. 

bool vd_Edge_desc::operator<(const vd_Edge_desc& e2)const {
  bool u1u2 = (this->U == e2.U);
  bool v1v2 = (this->V == e2.V);
  //std::cout << this->U << " - "<< this->V << " vs ";
  //std::cout << e2.U << " - "<< e2.V << std::endl;

  if (u1u2) {
    if (v1v2) {
      //std::cout << "Same edge"<< std::endl;
      return false;
    }
    else {
      return (this->V < e2.V);
    }
  }
  else {
    return (this->U < e2.U);
  }
}

bool vd_Edge_desc::operator==(const vd_Edge_desc& e2)const {

  bool u1u2 = (this->U == e2.U);
  bool v1v2 = (this->V == e2.V);
  return (u1u2 and v1v2);
}


int vd_Edge_desc::source () const {
  return U;
}

int vd_Edge_desc::target() const {
  return V;
}

// Default constructor:
vd_Edge_desc::vd_Edge_desc() {
  U = 0;
  V = 0;
}

vd_Edge_desc::vd_Edge_desc(int U_in, int V_in)
   :  U(U_in), V(V_in) {
}

// Copy constructor:
vd_Edge_desc::vd_Edge_desc(const vd_Edge_desc &obj) 
  : U(obj.U), V(obj.V) {
}

vd_Edge_desc& vd_Edge_desc::operator=( const vd_Edge_desc& obj ) {

  U = obj.U;
  V = obj.V;
  return *this;
}

bool chk_edge_joint(vd_Edge &e1, vd_Edge &e2) {
  return (e1.first.source() == e2.first.source() or
          e1.first.source() == e2.first.target() or
          e1.first.target() == e2.first.source() or
          e1.first.target() == e2.first.target());
}

// Assuming the edges are joint (share a vertex), join the edges. Joining an
// edge with itself creates an edge starting and ending at the target vertex of
// the edge.
vd_Edge join_edges(vd_Edge &e1, vd_Edge &e2) {
  if(e1.first.source() == e2.first.source()) {
    if(e1.first.target() < e2.first.target())
      return vd_Edge(e1.first.target(), e2.first.target(), e1.second);
    else
      return vd_Edge(e2.first.target(), e1.first.target(), e1.second);
  }
  else if(e1.first.source() == e2.first.target()) {
    if(e1.first.target() < e2.first.source())
      return vd_Edge(e1.first.target(), e2.first.source(), e1.second);
    else
      return vd_Edge(e2.first.source(), e1.first.target(), e1.second);
  }
  else if(e1.first.target() == e2.first.source()) {
    if(e1.first.source() < e2.first.target())
      return vd_Edge(e1.first.source(), e2.first.target(), e1.second);
    else
      return vd_Edge(e2.first.target(), e1.first.source(), e1.second);
  }
  else if(e1.first.target() == e2.first.target()) {
    if(e1.first.source() < e2.first.source())
      return vd_Edge(e1.first.source(), e2.first.source(), e1.second);
    else
      return vd_Edge(e2.first.source(), e1.first.source(), e1.second);
  }
  return vd_Edge(-1,-1,false);
}

// Quicksort implementation for an edge list.
void sort_edge(std::vector<vd_Edge>* edge_list, int left, int right) {
  if(right == -1 and left == -1) {
    left = 0;
    right = edge_list->size()-1;
  }

  // left > right
  if (left < right) {
    int part = partition_edge_list(edge_list, left, right);
    sort_edge(edge_list, left, part - 1);
    sort_edge(edge_list, part + 1, right);
  }
}

int partition_edge_list(std::vector<vd_Edge>* edge_list, int left, int right) {
  vd_Edge pivot = edge_list->at(right);
  // move the mid point value to the front.
  int i = left-1;
  int j = left;
  for (; j < right; j++) {
    if(edge_comp(edge_list->at(j),pivot)) {
      i++;
      std::swap(edge_list->at(i), edge_list->at(j));
    }
  }
  std::swap(edge_list->at(i+1),edge_list->at(right));
  return i + 1;
}

// Check if a given list of edges is not composed of a single circuit.
bool chk_discon_circ(std::vector<vd_Edge> temp) {

  bool fail = false;
  while(!fail and !temp.empty()) {
    if(temp.size() > 1) {
      if(!chk_edge_joint(temp.at(0), temp.at(1)))
        fail = true;
      else {
        temp.at(0) = join_edges(temp.at(0), temp.at(1));
        temp.at(1) = temp.back();
        temp.pop_back();
        sort_edge(&temp);
      }
    }
    else {
      if(temp.at(0).first.source() != temp.at(0).first.target() )
        fail = true;
      temp.clear();
    }
  }
  return fail;
}

std::ostream& operator<<(std::ostream& s, 
                     const vd_Edge_desc& e) {
  s << "[" << e.source() << ", " << e.target() << "]";
  return s;
}

// ---------------------------------------------------------
// ngon_gmi.
// ---------------------------------------------------------

ngon_gmi::ngon_gmi() : cp(std::make_pair(-1,-1)), 
  cells(0, std::make_pair(std::vector< int > (0), std::vector< int > (0))),
  ngons(0, std::vector< int > (0)) {
}

// Copy constructor:
ngon_gmi::ngon_gmi(const ngon_gmi& obj) : cp(std::make_pair(-1,-1)), 
  cells(0, std::make_pair(std::vector< int > (0), std::vector< int > (0))),
  ngons(0, std::vector< int > (0)) {
  clear();

  ngons.resize(obj.ngons.size());
  for(int i = 0; i < ngons.size(); i++) {
    ngons.at(i) = obj.ngons.at(i);
  }

  cells.resize(obj.cells.size());
  for(int i = 0; i < cells.size(); i++) {
    cells.at(i).first = obj.cells.at(i).first;
    cells.at(i).second = obj.cells.at(i).second;
  }
}

ngon_gmi& ngon_gmi::operator=( const ngon_gmi& obj ) {
  clear();

  cp = obj.cp;

  ngons.resize(obj.ngons.size());
  for(int i = 0; i < ngons.size(); i++) {
    ngons.at(i) = obj.ngons.at(i);
  }

  cells.resize(obj.cells.size());
  for(int i = 0; i < cells.size(); i++) {
    cells.at(i).first = obj.cells.at(i).first;
    cells.at(i).second = obj.cells.at(i).second;
  }
  return *this;
}

void ngon_gmi::clear() {
  for(int i = 0; i < ngons.size(); i++) {
    ngons.at(i).clear();
  }
  ngons.clear();

  for(int i = 0; i < cells.size(); i++) {
    cells.at(i).first.clear();
    cells.at(i).second.clear();
  }
  cells.clear();

  bd = false;
  curr = 0;
}

void ngon_gmi::shift_cells(int shift) {
  for(int i = 0; i < cells.size(); i++) {
    for(int j = 0; j < cells.at(i).first.size(); j++) {
      assert(cells.at(i).first.at(j) != 0);
      if(cells.at(i).first.at(j) > -1)
        cells.at(i).first.at(j) = cells.at(i).first.at(j) + shift;
    }
    for(int j = 0; j < cells.at(i).second.size(); j++) {
      cells.at(i).second.at(j) = cells.at(i).second.at(j) + shift;
    }
  }
}

ngon_gmi::~ngon_gmi() {
  clear();
}

void ngon_gmi::print() {
  if(bd)
    std::cout << "Boundary 0cell!" << std::endl;

  std::cout << "3cell couple " << std::endl;
  std::cout << cp.first << " " << cp.second << std::endl;

  std::cout << "Paths " << std::endl;
  for(int i = 0; i < cells.size(); i++) {
    std::cout << "\t " << i << ": 3cells" << std::endl;
    for(int j = 0; j < cells[i].first.size(); j++) {
      std::cout << "\t\t" << cells[i].first[j] << " " << std::endl;
    }

    std::cout << "\t " << i << ": 2cells" << std::endl;
    for(int j = 0; j < cells[i].second.size(); j++) {
      std::cout << "\t\t" << cells[i].second[j] << " " << std::endl;
    }
  }

  std::cout << "Ngons " << std::endl;
  for(int i = 0; i < ngons.size(); i++) {
    std::cout << "\t " << i << ": " << std::endl;
    for(int j = 0; j < ngons[i].size(); j++) {
      std::cout << "\t\t" << ngons[i][j] << " " << std::endl;
    }

  }
}


cell_path::cell_path() : first(std::make_pair(-1,-1)), 
  second(0) {
}

cell_path::~cell_path() {
  for(int i = 0; i < second.size(); i++) {
    second.at(i).first.clear();
    second.at(i).second.clear();
  }
  second.clear();
}

// Copy constructor:
cell_path::cell_path(const cell_path &obj) : first(std::make_pair(-1,-1)), 
  second(0) {
  second.resize(obj.second.size());
  for(int i = 0; i < second.size(); i++) {
    second.at(i).first = obj.second.at(i).first;
    second.at(i).second = obj.second.at(i).second;
  }
  first = obj.first;
}

cell_path& cell_path::operator=( const cell_path& obj ) {
  second.resize(obj.second.size());
  for(int i = 0; i < second.size(); i++) {
    second.at(i).first = obj.second.at(i).first;
    second.at(i).second = obj.second.at(i).second;
  }
  first = obj.first;
  return *this;
}


std::ostream& operator<<(std::ostream& s, 
                     const cell_path& cp) {
  s << "Couple " << cp.first.first << " - " << cp.first.second << ":\n";
  s << "Paths:\n";
  for(int i = 0; i < cp.second.size(); i++) {
    s << i << ": ";
    for(int j = 0; j < cp.second.at(i).first.size(); j++) {
      s << cp.second.at(i).first.at(j) << " ";
    }
    s << "\n\tNon-intersecting:";
    for(int j = 0; j < cp.second.at(i).second.size(); j++) {
      s << cp.second.at(i).second.at(j) << " ";
    }
    s << "\n";
  }
  return s;
}


// ---------------------------------------------------------
// Graph based objects.
// ---------------------------------------------------------


std::pair<std::set<int>::const_iterator, std::set<int>::const_iterator> 
            vd_graph::vertices() {
  std::set<int>::const_iterator v_beg = v_id.begin();
  std::set<int>::const_iterator v_end = v_id.end();
  return std::make_pair(v_beg, v_end);
}

std::pair<std::set<int>::const_iterator, std::set<int>::const_iterator> 
            vd_graph::adjacent_vertices(int V) {

  std::set<int>::const_iterator v_beg = adj_vec.at(adj_map[V]).begin();
  std::set<int>::const_iterator v_end = adj_vec.at(adj_map[V]).end();
  return std::make_pair(v_beg, v_end);
}

std::pair<std::multiset<vd_Edge_desc >::iterator, 
          std::multiset<vd_Edge_desc >::iterator > vd_graph::edge_it() {
  std::multiset<vd_Edge_desc>::const_iterator e_beg = edge_list.begin();
  std::multiset<vd_Edge_desc>::const_iterator e_end = edge_list.end();
  return std::make_pair(e_beg, e_end);
}

int vd_graph::source(vd_Edge* e) const {
  //if (edge_list.find(e->first) != edge_list.end())
    return e->first.source();
}

int vd_graph::target(vd_Edge* e) const {
  //if (edge_list.find(e->first) != edge_list.end())
    return e->first.target();
}

// Add-remove, check vertex:
void vd_graph::add_vertex(int id) {
  if (v_id.find(id) == v_id.end()) {
    v_id.insert(id);
    std::set<int > v_set{};
    adj_vec.push_back(v_set);

    adj_map[id] = adj_vec.size()-1;
  }
}

bool vd_graph::check_vertex(int id) {
  if (v_id.find(id) != v_id.end())
    return true;
  return false;
}

// Add-remove, check vertex:
void vd_graph::remove_vertex(int id) {
  if (check_vertex(id)) {

    std::map<int, int>::iterator it;

    int id_temp = -1;
    for(it = adj_map.begin(); it != adj_map.end(); it++) {
      if(it->second == adj_vec.size()-1) {
        id_temp = it->first;
        break;
      }
    }
    adj_map[id_temp] = adj_map[id];
    adj_vec.at(adj_map[id]) = adj_vec.back();
    adj_vec.pop_back();

    v_id.erase(id);
  }
}


// Add-remove edge.
bool vd_graph::remove_edge(int v1, int v2) {
  if (directed) 
    return remove_edge_dir(v1,v2);
  else {
    vd_Edge_desc e(v1,v2);
    if (v1 < v2)
      e = vd_Edge_desc(v1, v2);
    else
      e = vd_Edge_desc(v2, v1);

    std::multiset<vd_Edge_desc>::iterator e_it = edge_list.find(e);
    if (e_it != edge_list.end()) {
      edge_list.erase(e_it);
      // If no other edge between vertices exists, remove their adjacency.
      e_it = edge_list.find(e);
      assert(e_it == edge_list.end());
      adj_vec.at(adj_map[v1]).erase(v2);
      adj_vec.at(adj_map[v2]).erase(v1);
      return true;
    }

    return false;
  }
}

bool vd_graph::remove_edge_dir(int v1, int v2) {
  vd_Edge_desc e(v1,v2);

  std::multiset<vd_Edge_desc>::iterator e_it = edge_list.find(e);
  if (e_it != edge_list.end()) {
    edge_list.erase(e_it);
    // If no other edge between vertices exists, remove their adjacency.
    e_it = edge_list.find(e);
    if (e_it == edge_list.end()) {
      adj_vec.at(adj_map[v1]).erase(v2);
    }
    return true;
  }
  return false;
}


vd_Edge vd_graph::add_edge(int v1, int v2) {

  if (!check_vertex(v1))
    add_vertex(v1);
  if (!check_vertex(v2))
    add_vertex(v2);

  if (directed)
    return add_edge_dir(v1,v2);
  else {
    vd_Edge_desc e(v1,v2);
    if (v1 < v2)
      e = vd_Edge_desc(v1, v2);
    else
      e = vd_Edge_desc(v2, v1);

    // Make sure that the vertices are adjacent correctly or become adjacent.
    if (adj_vec.at(adj_map[v1]).find(v2) != adj_vec.at(adj_map[v1]).end())
      assert(adj_vec.at(adj_map[v2]).find(v1) != 
                                      adj_vec.at(adj_map[v2]).end());
    else {
      assert(adj_vec.at(adj_map[v2]).find(v1) == 
                                      adj_vec.at(adj_map[v2]).end());
      adj_vec.at(adj_map[v1]).insert(v2);
      adj_vec.at(adj_map[v2]).insert(v1);
    }
    edge_list.insert(e);
    return std::make_pair(e, true);
  }
}

vd_Edge vd_graph::add_edge_dir(int v1, int v2) {
  vd_Edge_desc e(v1,v2);

  // Make sure that the vertices are adjacent correctly or become adjacent.
  if (adj_vec.at(adj_map[v1]).find(v2) == adj_vec.at(adj_map[v1]).end())
    adj_vec.at(adj_map[v1]).insert(v2);

  edge_list.insert(e);
  return std::make_pair(e, true);
}

// Undirected graph check edge.
bool vd_graph::check_edge(int v1, int v2) {
  if (v2 < v1) {
    int temp = v1; 
    v1 = v2;
    v2 = temp;
  }
  if (check_vertex(v1)) {
    if (check_vertex(v2)) {
      std::set<int>::iterator it;
      it = adj_vec.at(adj_map[v1]).find(v2);
      if(it != adj_vec.at(adj_map[v1]).end()) {
        it = adj_vec.at(adj_map[v2]).find(v1);
        assert(it != adj_vec.at(adj_map[v2]).end());
        vd_Edge_desc e(v1, v2);
        std::multiset<vd_Edge_desc>::iterator e_it = edge_list.find(e);
        assert(e_it != edge_list.end());
        return true;
      }
    }
  }
  return false;
}

// Get edge. Non existing edges will have boolean part false.
vd_Edge vd_graph::get_edge(int v1, int v2) {
  if(directed)
    return get_edge_dir(v1,v2);
  else {
    vd_Edge_desc e(v1,v2);
    if (v1 < v2)
      e = vd_Edge_desc(v1, v2);
    else
      e = vd_Edge_desc(v2, v1);

    if (edge_list.find(e) != edge_list.end())
      return std::make_pair(e, true);

    return std::make_pair(e, false);
  }
}

vd_Edge vd_graph::get_edge_dir(int v1, int v2) {
  vd_Edge_desc e(v1,v2);

  if (edge_list.find(e) != edge_list.end())
    return std::make_pair(e, true);

  return std::make_pair(e, false);
}

// Assuming given node has exactly two connections, given a node that is 
// connected return the other connection and assert crash otherwise.
int vd_graph::get_node_other(int n2, int n_end) {
  std::set<int>::const_iterator ai_int;
  std::set<int>::const_iterator ai_end_int;
  std::tie(ai_int, ai_end_int) = adjacent_vertices(n2);
  assert(std::distance(ai_int, ai_end_int) == 2);

  // Check the first element. If the first matches, return the second.
  if(*ai_int == n_end) {
    ++ai_int;
    return *ai_int;
  }
  // Assert the second mathes the first. Return the first.
  else {
    int temp = *ai_int;
    ++ai_int;
    assert(*ai_int == n_end);
    return temp;
  }
}


// Get the number of vertices, edges, graph type (undirected, directed).
int vd_graph::num_vertices() {
  return v_id.size();
}

int vd_graph::num_edges() {
  return edge_list.size();
}

bool vd_graph::dir_type() {
  return directed;
}

void vd_graph::set_dir(bool dir) {
  clear();
  directed = dir;
}


// Default constructor.
vd_graph::vd_graph() : directed(false), v_id{}, reserve_sz(0),
                      adj_vec(0, std::set<int>{}), adj_map{}, edge_list{} {
}

vd_graph& vd_graph::operator=( vd_graph& obj ) {
  clear();

  directed = obj.directed;
  v_id = obj.v_id;

  reserve_sz = obj.reserve_sz;
  std::set<int>::const_iterator id_beg = obj.v_id.begin();
  adj_map = obj.adj_map;

  adj_vec.resize(obj.adj_vec.size());
  for(int i = 0; i < adj_vec.size(); i++) {
    adj_vec.at(i) = obj.adj_vec.at(i);
  }

  std::multiset<vd_Edge_desc>::const_iterator e_beg = obj.edge_list.begin();

  // Iterate over the adjacency map, copy each vertex adjacency.
  while (e_beg != obj.edge_list.end()) {
    edge_list.insert(*e_beg);
    e_beg++;
  }
  return *this;
}

// Clear the content of the graph.
void vd_graph::reserve_v(int sz) {
  reserve_sz = sz;
  adj_vec.reserve(reserve_sz);
}

// Clear the content of the graph.
void vd_graph::clear() {
  directed = false;
  v_id.clear();
  //index.clear();
  adj_map.clear();
  for(int i = 0; i < adj_vec.size(); i++)
    adj_vec.at(i).clear();
  adj_vec.clear();

  edge_list.clear();
}

// Test function to print the current graph.
void vd_graph::print_adj() {
  std::set<int>::const_iterator ai_int;
  std::set<int>::const_iterator ai_end_int;

  for (std::set<int>::iterator v_it = v_id.begin(); 
                               v_it != v_id.end(); v_it++) {
    std::cout << *v_it << ": ";
    for (std::tie(ai_int, ai_end_int) = adjacent_vertices(*v_it); 
                              ai_int != ai_end_int; ++ai_int) {
      std::cout << *ai_int << ", ";
    }
      std::cout << std::endl;
  }

  std::cout << "Edges: " << std::endl;
  for (  std::multiset<vd_Edge_desc >::iterator e_it = edge_list.begin();
                                           e_it != edge_list.end();
                                           e_it++) {
    std::cout << e_it->source() << " - " << e_it->target() << std::endl;
  }
}

// Test function to print the current graph. hl is used to highlight edges.
void vd_graph::print_graph(const char* dotname, std::vector<vd_Edge>* hl) {
  std::ofstream fout(dotname);
  std::map<vd_Edge_desc, bool> hl_map{};
  if(hl != NULL) {
    for (int i = 0; i < hl->size(); i++) {
      hl_map[hl->at(i).first] = true;
    }
  }

  if (directed)
    fout << "digraph A {\n";
  else
    fout << "graph A {\n";

  fout << "  rankdir=LR\n"
    << "ratio=\"fill\"\n"
    << "edge[style=\"bold\"]\n" << "node[shape=\"circle\"]\n";

  std::multiset<vd_Edge_desc >::iterator it = edge_list.begin();
  std::multiset<vd_Edge_desc >::iterator end = edge_list.end();

  if (directed) {
    for (; it != end; ++it) {
      fout << it->source() << " -> " << it->target() 
           << (hl_map[*it] ? " [color=\"red\"]" : "")
           << "\n";
    }
  }
  else {
    for (; it != end; ++it) {
      fout << it->source() << " -- " << it->target()
           << (hl_map[*it] ? " [color=\"red\"]" : "")
           << "\n";
    }
  }
  fout << "}\n";
  fout.close();
  //std::ofstream fout2("output/graph2.dot");
  //write_graphviz(fout2, g);
}

// Test function to print the current graph. hl is used to highlight edges.
void vd_graph::print_graph_p(std::map<vd_Edge_desc, bool>& hl_map, std::map<int, std::string>& color_map, std::map<vd_Edge_desc, int>& path_map, const char* dotname) {
  std::ofstream fout(dotname);

  if (directed)
    fout << "digraph A {\n";
  else
    fout << "graph A {\n";

  fout << "  rankdir=LR\n"
    << "ratio=\"fill\"\n"
    << "edge[style=\"bold\"]\n" << "node[shape=\"circle\"]\n";

  std::multiset<vd_Edge_desc >::iterator it = edge_list.begin();
  std::multiset<vd_Edge_desc >::iterator end = edge_list.end();

  if (directed) {
    for (; it != end; ++it) {
      fout << it->source() << " -> " << it->target() 
           << (hl_map[*it] ? color_map[path_map[*it]] : "")
           << "\n";
    }
  }
  else {
    for (; it != end; ++it) {
      fout << it->source() << " -- " << it->target()
           << (hl_map[*it] ? color_map[path_map[*it]] : "")
           << "\n";
    }
  }
  fout << "}\n";
  fout.close();
  //std::ofstream fout2("output/graph2.dot");
  //write_graphviz(fout2, g);
}
// Print the current graph using labels. 
void vd_graph::print_graph(std::vector<std::string>* lbls, 
                          const char* dotname = "output/graph.dot" ) {
  assert(lbls->size() == v_id.size());
  std::ofstream fout(dotname);

  if (directed)
    fout << "digraph A {\n";
  else
    fout << "graph A {\n";

  fout << "  rankdir=LR\n"
    << "ratio=\"fill\"\n"
    << "edge[style=\"bold\"]\n" << "node[shape=\"circle\"]\n";

  std::multiset<vd_Edge_desc >::iterator it = edge_list.begin();
  std::multiset<vd_Edge_desc >::iterator end = edge_list.end();


  if (directed) {
    for (; it != end; ++it)
      fout << lbls->at(it->source()) << " -> " 
           << lbls->at(it->target()) << "\n";
  }
  else {
    for (; it != end; ++it)
      fout << lbls->at(it->source()) << " -- " 
           << lbls->at(it->target()) << "\n";
  }

  fout << "}\n";
  fout.close();

  //std::ofstream fout2("output/graph2.dot");
  //write_graphviz(fout2, g);
}


// ---------------------------------------------------------
// Paton elementary circuit finder and Gibbs algorithm.
// ---------------------------------------------------------

bool edge_comp (const vd_Edge& e1, const vd_Edge& e2) { 
  return e1.first < e2.first;
}

//bool PatonFinder::edge_comp (const vd_Edge& e1, const vd_Edge& e2) { 
//  return e1.first < e2.first;
//}

// Assume sorted. If all edges are the same, they will be sorted the same.
// The comparison should hold.
bool set_edge_comp (std::vector<vd_Edge> e1, std::vector<vd_Edge> e2) { 

  // This is used to sort edge sets such that the shorter sets will be 
  // considered "less".
  if (e1.size() < e2.size()) {
    return true;
  }
  else if (e1.size() > e2.size()) {
    return false;
  }

  std::vector<vd_Edge>::iterator it_e1;
  std::vector<vd_Edge>::iterator it_e2;

  it_e1 = e1.begin(); 
  it_e2 = e2.begin(); 
  while (it_e1 != e1.end() and it_e2 != e2.end()) {
    int u1 = it_e1->first.source();
    int u2 = it_e2->first.source();
    int v1 = it_e1->first.target();
    int v2 = it_e2->first.target();

    bool u1u2 = (u1 == u2);
    bool u1v2 = (u1 == v2);
    bool v1u2 = (v1 == u2);
    bool v1v2 = (v1 == v2);

    if ((u1u2 and v1v2) or (u1v2 and v1u2)) {
      it_e1++;
      it_e2++;
    }
    else {
      return edge_comp (*it_e1, *it_e2);
    }
  }
  if (it_e1 != e1.end()) {
    return false;
  }
  else if (it_e2 != e2.end()) {
    return true;
  }

  return false;
}

void PatonFinder::sort_edge_sets(std::vector<std::vector<vd_Edge> >* e_set, int left, int right) {
  if(right == -1 and left == -1) {
    left = 0;
    right = e_set->size()-1;
  }
  //std::cout << "Sorting edge set " << e_set << std::endl;
  //std::cout << left << " " << right << std::endl;

  if (left < right) {
    int part = partition_edge_set(e_set, left, right);
    sort_edge_sets(e_set, left, part - 1);
    sort_edge_sets(e_set, part + 1, right);
  }
}
//Function to determine the partitions
// partitions the array and returns the middle subscript
int PatonFinder::partition_edge_set(std::vector<std::vector<vd_Edge> >* e_set, int left, int right) {
  std::vector<vd_Edge> pivot(0);
  pivot = e_set->at(right);
  // move the mid point value to the front.
  int i = left-1;
  int j = left;
  for (; j < right; j++) {
    if(set_edge_comp(e_set->at(j), pivot)) {
      i++;
      std::swap(e_set->at(i), e_set->at(j));
    }
  }
  std::swap(e_set->at(i+1),e_set->at(right));
  return i + 1;
}

void PatonFinder::isect_edge (std::vector<vd_Edge>* e1, 
                              std::vector<vd_Edge>* e2, 
                 std::vector<vd_Edge>* out) {
  set_intersection(e1->begin(),e1->end(), e2->begin(),e2->end(), 
      std::inserter(*out,out->begin()), 
                edge_comp);
            //std::bind(&PatonFinder::edge_comp, this, 
            //std::placeholders::_1, std::placeholders::_2));
}

void PatonFinder::diff_edge (std::vector<vd_Edge>* e1, 
                             std::vector<vd_Edge>* e2, 
                 std::vector<vd_Edge>* out) {
  set_difference(e1->begin(),e1->end(),
                    e2->begin(),e2->end(), 
                    std::inserter(*out,out->begin()), 
                    edge_comp);
            //std::bind(&PatonFinder::edge_comp, this, 
            //std::placeholders::_1, std::placeholders::_2));
}
void PatonFinder::merge_edge (std::vector<vd_Edge>* e1, 
                              std::vector<vd_Edge>* e2, 
                 std::vector<vd_Edge>* out) {
  std::merge(e1->begin(),e1->end(),
                    e2->begin(),e2->end(), 
                std::inserter(*out,out->end()), 
                edge_comp);

            //std::bind(&PatonFinder::edge_comp, this, 
            //std::placeholders::_1, std::placeholders::_2));
}

void PatonFinder::merge_edge_set (std::vector<std::vector<vd_Edge> >* e1, 
                 std::vector<std::vector<vd_Edge> >* e2, 
                 std::vector<std::vector<vd_Edge> >* out) {
  //std::cout << "Merging e sets" << e1 << " " << e2 << std::endl;
  std::merge(e1->begin(),e1->end(),
                e2->begin(),e2->end(), 
                    std::inserter(*out,out->end()), 
                    set_edge_comp);
            //std::bind(&PatonFinder::set_edge_comp, this, 
            //std::placeholders::_1, std::placeholders::_2));
}

PatonFinder::PatonFinder() :
    g(), tree(),
    T{}, V{}, 
    paths(0, std::vector<int>(0) ), paths_e(0, std::vector<vd_Edge>(0)),
    paths_eset(0, std::vector<int>(0)), S(0, std::vector<vd_Edge>(0)),
    disj_frag(0, std::vector<int>(0)),
    circ_tup(0, std::make_pair(std::vector<vd_Edge>(0), s_graph()) ),
    index{} {
}

PatonFinder::PatonFinder(vd_graph* g_in) :
    g(), tree(),
    T(), V(), 
    paths(0, std::vector<int>(0) ), paths_e(0, std::vector<vd_Edge>(0)),
    paths_eset(0, std::vector<int>(0)), S(0, std::vector<vd_Edge>(0)),
    disj_frag(0, std::vector<int>(0)),
    circ_tup(0, std::make_pair(std::vector<vd_Edge>(0), s_graph()) ),
    index{} {

  g_init = g_in;
  g = *g_in;
  //index = get(boost::vertex_index, g);
  //std::cout << g.num_vertices() << std::endl;

  if (g.num_vertices() > 0) {
    for (vp = g.vertices(); vp.first != vp.second; ++vp.first) {
      //std::cout << *vp.first <<  " ";
      V.insert(*vp.first);
    }
    --vp.first;
    T.insert(*vp.first);
    tree.add_vertex(*vp.first);
    //std::cout << std::endl;
  }
  else {
    //std::cout << "The graph is empty" << std::endl;
  }
}

void PatonFinder::reload_graph(vd_graph* g_in) {
  clear();

  g_init = g_in;
  g = *g_init;
  //index = get(boost::vertex_index, g);
  //std::cout << g.num_vertices() << std::endl;

  if (g.num_vertices() > 0) {
    for (vp = g.vertices(); vp.first != vp.second; ++vp.first) {
      //std::cout << *vp.first <<  " ";
      V.insert(*vp.first);
    }
    // Insert the last vertex as the root of the tree.
    --vp.first;
    T.insert(*vp.first);
    //std::cout << std::endl;
  }
  else {
    //std::cout << "The graph is empty" << std::endl;
  }
}

// Elementary circuit detection.
bool PatonFinder::find_nonrecursive(int v1, int v2) {

  std::vector<int> path(0);
  std::map<int, int> predecessor{};
  std::map<int, bool> visited{};

  std::vector<int> cand(0);
  predecessor.clear();
  cand.clear();
  cand.push_back(v1);

  visited[v1] = true;

  // Visit every accessible vertex, tag them as visited and record their
  // predecessor. At any given level, each vertex will be connected to 
  // its closest neighbor.
  while (!cand.empty()) {
    for (std::tie(ai_int, ai_end_int) = tree.adjacent_vertices(cand.front()); ai_int != ai_end_int; ++ai_int) {
      if (!visited[*ai_int]) {
        cand.push_back(*ai_int);
        visited[*ai_int] = true;
        predecessor[*ai_int] = cand.front();
      }
    }
    cand.erase(cand.begin());
  }

  //if(visited[v2]){
    //std::cout << "Visited " << v2 << std::endl;
  //}

  // Collect the vertices starting from v2.
  path.push_back(v2);
  while (path.back() != v1) {
    path.push_back(predecessor[path.back()]);
    //std::cout << " " << path.back() << std::endl;
  }

  // Vector to set conversion. Could be faster when combining different 
  // circuits.
  std::set<int> path_set(std::make_move_iterator(path.begin()),
                std::make_move_iterator(path.end()));

  std::set<int> intersect{};
  set_intersection(V.begin(),V.end(),path_set.begin(),path_set.end(),
              std::inserter(intersect,intersect.begin()));

  //std::cout << std::endl;
  //while(!intersect.empty()) {
    //std::cout << " " << *intersect.begin() << std::endl;
  //  intersect.erase(*intersect.begin());
  //}
  //std::cout << std::endl;
  paths.push_back(path);

  return true;
}


void PatonFinder::find_paton() {
  std::set<int> X(std::make_move_iterator(V.begin()),
                std::make_move_iterator(V.end()));

  std::set<int> intersect{};
  set_intersection(X.begin(),X.end(),T.begin(),T.end(),
              std::inserter(intersect,intersect.begin()));

  //std::cout << "Size: " << intersect.size();
  //std::cout << "max: " << intersect.max_size();
  //std::cout << std::endl;
  //int map_nbr = 0;

  std::map<int, int> ord_map{};
  std::map<int, bool> used_map{};
  int last_id = 1;
  assert(T.size() == 1);
  ord_map[last_id] = *T.begin();

  while(!intersect.empty()) {
    int id_c = ord_map[last_id];
    if(!used_map[id_c]) {
      std::tie(ai, ai_end) = g.adjacent_vertices(id_c);
      while (ai != ai_end) {
      //for (std::tie(ai, ai_end) = g.adjacent_vertices(id_c); ai != ai_end; ++ai) {
        if (T.find(*ai) == T.end()) {
          T.insert(*ai);
          intersect.insert(*ai);

          tree.add_vertex(*ai);
          tree.add_edge(id_c, *ai);
          last_id = last_id + 1;
          ord_map[last_id] = *ai;
        }
        else {
          find_nonrecursive(id_c, *ai);
        }
        std::set<int>::const_iterator ai_curr = ai;
        ++ai;
        g.remove_edge(id_c, *ai_curr);
        //std::cout << "After: g: " << std::endl;
        //g.print_adj();
        //std::cout << "After: tree: " << std::endl;
        //tree.print_adj();

      }
      used_map[id_c] = true;
      //std::cout << std::endl;
      X.erase(X.find(id_c));
      intersect.erase(intersect.find(id_c));
    }
    else {
      last_id = last_id - 1;
    }

  }

  for (int i = 0; i < paths.size(); i++) {
    std::cout << "Cycle " << i << ": ";
    for (int j = 0; j < paths.at(i).size(); ++j) {
      std::cout << paths[i][j] << " ";
    }
    std::cout << std::endl;
  }

}


// Convert the elementary circuit lists into edge circuits lists.
// These are used in finding circuits by gibbs algorithm.
void PatonFinder::vert_2_edge() {
  g = *g_init;

  for (int i = 0; i < paths_e.size(); i++)
    paths_e.at(i).clear();
  paths_e.clear();

  paths_e.resize(paths.size());
  for (int i = 0; i < paths_e.size(); i++) {
    paths_e.at(i).resize(paths.at(i).size());
  }

  for (int i = 0; i < paths_e.size(); i++) {
    for (int j = 0; j < paths.at(i).size()-1; j++) {
      paths_e.at(i).at(j) = g.get_edge(paths.at(i).at(j), paths.at(i).at(j+1));
    }
    paths_e.at(i).back() = g.get_edge(paths.at(i).at(0), paths.at(i).back());
  }

  for (int i = 0; i < paths_e.size(); i++) {
    sort_edge (&paths_e.at(i));
  }

}

// Given two vector list of cells and a function that tests connectivity of 
// cells in these lists, return true if they have at least a connection.
// Mainly used in collecting the disjoint graphs obtained by removing the  
// circuit into at most two groups.
// The 2cells have the graph ids and ids lower than c3_nbr are discarded.
bool PatonFinder::chk_conn_list(std::vector<int>* c1_list,
      std::vector<int>* c2_list) {

  for(int i = 0; i < c1_list->size(); i++) {
    for(int j = 0; j < c2_list->size(); j++) {
      if(c_graph->comp_conn_12(c1_list->at(i), c2_list->at(j)) )
        return true;
    }
  }
  return false;
}

// Create a subgraph for each path. 
// For each path: 
// For each edge on the path remove the edges connecting the source and target
// vertices. 
// Using the burn algorithm, find disjoint subgraphs. 
// Reload the graph.
// Print out the disjoint subgraphs and paths.
void PatonFinder::create_sub(cell_graph* const c_graph_in) {
  g = *g_init;

  c_graph = c_graph_in;
  circ_tup.resize(S.size());

  for (int i = 0; i < S.size(); i++) {
    //std::cout << "Circuit " << i <<": ";
    circ_tup.at(i).first.resize(S.at(i).size());

    for (int j = 0; j < S.at(i).size(); j++) {
      int u=g.source(&S.at(i).at(j));
      int v=g.target(&S.at(i).at(j));
      //std::cout << u << " - " << v << ", " ;
      circ_tup.at(i).first.at(j) = S.at(i).at(j);
    }
    //std::cout << std::endl;
  }

  // Stores the disjoints graphs as couples.

  for (int i = 0; i < circ_tup.size(); i++) {
    circ_tup.at(i).second.first.reserve(g.num_vertices());
    circ_tup.at(i).second.second.reserve(g.num_vertices());
  }

  // Remove all edges connecting the vertices on the path and the rest of the
  // original graph. Burn the subgraphs.
  for (int i = 0; i < S.size(); i++) {

    // To check if the vertex is burnable. No vertex on path or already burned
    // is burnable.
    std::map<int, bool> burned{};
    // After the vertices of the circuit are removed, the remaining vertices of 
    // the graph are grouped in a number of disjoint graphs. 
    // Burn the vertices of the same trial with the same burn id.
    std::map<int, int> burn_id{};

    for (int j = 0; j < S.at(i).size(); j++) {
      int u=g.source(&S.at(i).at(j));
      int v=g.target(&S.at(i).at(j));
      if (!burned[u]) {
        //std::cout << "Removing " << u << " on circuit." << std::endl;
        //for (std::tie(ai, ai_end) = g.adjacent_vertices(u); ai != ai_end; ++ai) {
        std::tie(ai, ai_end) = g.adjacent_vertices(u);
         while (ai != ai_end) {
          std::set<int>::const_iterator ai_curr = ai;
          ++ai;
          g.remove_edge(u, *ai_curr);
        }
        burned[u] = true;
        burn_id[u] = -1;
      }
      if (!burned[v]) {
        //std::cout << "Removing " << v << " on circuit." << std::endl;
        //for (std::tie(ai, ai_end) = g.adjacent_vertices(v); ai != ai_end; ++ai) {
        //  g.remove_edge(v, *ai);
        std::tie(ai, ai_end) = g.adjacent_vertices(v);
         while (ai != ai_end) {
          std::set<int>::const_iterator ai_curr = ai;
          ++ai;
          g.remove_edge(v, *ai_curr);
        }
        burned[v] = true;
        burn_id[v] = -1;
      }
    }

    // Flip true if all passed vertices are burned.
    bool passed = false;
    // id of the current burning disjoint graph.
    int id = 0;

    while(!passed) {
      passed = true;
      // Burn algorithm.
      // Find an unburned vertex on the graph. 
      std::vector<int> burning(0);
      burning.clear();
      burning.reserve(g.num_vertices());

      //std::cout << "Burning the first subgraph..." << std::endl;
      // Find a vertex on the subgraph.
      for (vp = g.vertices(); vp.first != vp.second; ++vp.first) {
        if (!burned[*vp.first]) {
          burning.push_back(*vp.first);
          burned[*vp.first] = true;
          burn_id[*vp.first] = id;
          //std::cout << "Burning " << *vp.first << " in " << id << std::endl;
          vp.first = vp.second;
          --vp.first;
          passed = false;
        }
      }

      // Burn the first subgraph.
      while(!burning.empty()) {
        int v = *burning.begin();
        //std::cout << "Neighbors of " << v << " in " << id << std::endl;
        for (std::tie(ai, ai_end) = g.adjacent_vertices(v); 
                                                      ai != ai_end; ++ai) {
          if (!burned[*ai]) {
            //std::cout << "Burning " << *ai << std::endl;
            burning.push_back(*ai);
            burned[*ai] = true;
            burn_id[*ai] = id;
          }
        }
        //std::cout << "Removing " << *burning.begin() << std::endl;
        //circ_tup.at(i).second.first.push_back(burning[0]);
        burning.erase(burning.begin());
      }

      // If an unburned patch is found, increase the patch id for the next run.
      if(!passed)
        id++;
    }

    //std::cout << id << " patches. " << std::endl;

    disj_frag.clear();
    disj_frag.resize(id);
    for (int j = 0; j < id; j++) {
      disj_frag.at(j).reserve(g.num_vertices());
    }

    for (vp = g.vertices(); vp.first != vp.second; ++vp.first) {
      //std::cout << *vp.first << " in " << burn_id[*vp.first] << std::endl;
      if(burn_id[*vp.first] > -1)
        disj_frag[burn_id[*vp.first]].push_back(*vp.first);
    }

    // If there are more than one disj_frag,
    while (disj_frag.size() > 2) {
      bool found = false;
      bool found2 = false;
      for (int j = 0; j < disj_frag.size() - 1; j++) {
        if(chk_conn_list(&disj_frag.back(), &disj_frag.at(j))) {
          std::copy(disj_frag.back().begin(), disj_frag.back().end(),
           std::inserter(disj_frag.at(j), disj_frag.at(j).end()));
          disj_frag.back().clear();
          disj_frag.pop_back();
          j = disj_frag.size();
          found = true;
        }
      }
      // Switch the first one with the last one.
      if(!found) {
        assert(!found2);
        found2 = true;
        std::vector<int> temp(0);
        temp = disj_frag.back();
        disj_frag.back() = disj_frag.at(0);
        disj_frag.at(0) = temp;
      }
    }
    // Check if the last two disj_frag are connected, as well:
    if(disj_frag.size() == 2) {
      bool found = false;
      bool found2 = false;
      if(chk_conn_list(&disj_frag.back(), &disj_frag.at(0))) {
        std::copy(disj_frag.back().begin(), disj_frag.back().end(),
         std::inserter(disj_frag.at(0), disj_frag.at(0).end()));
        disj_frag.back().clear();
        disj_frag.pop_back();
        found = true;
      }
    }
    // TODO is it ever NULL?
    if (c_graph != NULL and id > 0) {
      circ_tup.at(i).second.first.reserve(g.num_vertices());
      circ_tup.at(i).second.second.reserve(g.num_vertices());

      std::copy(disj_frag.at(0).begin(), disj_frag.at(0).end(),
                 std::inserter(circ_tup.at(i).second.first,
                                circ_tup.at(i).second.first.begin()));

      if (disj_frag.size() == 2) {
        std::copy(disj_frag.at(1).begin(), disj_frag.at(1).end(),
                 std::inserter(circ_tup.at(i).second.second,
                                  circ_tup.at(i).second.second.begin()));
      }
    }
    // Reload the graph for the next path.
    g = *g_init;
  }

  // Print the paths and associated subgraphs.
  /*
  for (int i = 0; i < S.size(); i++) {
    std::cout << i << "th circuit: " << std::endl;
    for (int j = 0; j < S[i].size(); j++) {
      std::cout << S[i][j].first << " ["
                << S[i][j].second << "] ";
    }
    std::cout << std::endl;

    std::cout << "Subgraph 1_tup: " << std::endl;
    for (int j = 0; j < circ_tup[i].second.first.size(); j++) {
      std::cout << circ_tup[i].second.first[j] << " ";
    }
    std::cout << std::endl;

    std::cout << "Subgraph 2_tup: " << std::endl;
    for (int j = 0; j < circ_tup[i].second.second.size(); j++) {
      std::cout << circ_tup[i].second.second[j] << " ";
    }
    std::cout << std::endl;
  }
  */
}

// Using the gibbs algorithm, find all unique circuits.
// TODO FIX the hackish solution for removal of non-circuit edge lists from the
// candidate list...
void PatonFinder::find_gibbs() {
  std::vector< std::vector<vd_Edge> > S_temp
    (0, std::vector<vd_Edge>(0));
  std::vector< std::vector<vd_Edge> > Q
    (0, std::vector<vd_Edge>(0));
  std::vector< std::vector<vd_Edge> > Q_temp
    (0, std::vector<vd_Edge>(0));
  std::vector< std::vector<vd_Edge> > Q_temp2
    (0, std::vector<vd_Edge>(0));
  std::vector< std::vector<vd_Edge> > R
    (0, std::vector<vd_Edge>(0));
  std::vector< std::vector<vd_Edge> > R_star
    (0, std::vector<vd_Edge>(0));
 
  std::vector<vd_Edge> diff(0);
  std::vector<vd_Edge> diff_join(0);
  std::vector<vd_Edge> intersect(0);

  std::map<int,int> self_int{};

  if (!paths_e.empty()) {
    S.push_back(paths_e.at(0));
    Q.push_back(paths_e.at(0));
  }

  //std::cout << "S-Q sizes" << S.size() << Q.size() << std::endl;

  std::vector<vd_Edge>::iterator it_e;

  for (int i = 1; i < paths_e.size(); i++) {

    for (int j = 0; j < S_temp.size(); j++) {
      S_temp.at(j).clear();
    }
    S_temp.clear();
    for (int j = 0; j < Q_temp.size(); j++) {
      Q_temp.at(j).clear();
    }
    Q_temp.clear();
    for (int j = 0; j < Q_temp2.size(); j++) {
      Q_temp2.at(j).clear();
    }
    Q_temp2.clear();
    for (int j = 0; j < R.size(); j++) {
      R.at(j).clear();
    }
    R.clear();
    for (int j = 0; j < R_star.size(); j++) {
      R_star.at(j).clear();
    }
    R_star.clear();

/*
    std::cout << "paths_e[" << i << "] = ";
    for (it_e = paths_e[i].begin(); 
         it_e != paths_e[i].end(); ++it_e) {
      std::cout << it_e->first << " ";
    }
    std::cout << std::endl;
*/
    // Compare with T in Q.
    for (int j = 0; j < Q.size(); j++) {
      diff.clear();
      diff_join.clear();
      intersect.clear();

      // path-T
      diff_edge(&paths_e.at(i),&Q.at(j),&diff);

      sort_edge (&diff);

      // T-path
      diff_edge(&Q.at(j),&paths_e.at(i),&intersect);

      sort_edge (&intersect);

      // (T U path) - (path n T)
      merge_edge(&diff,&intersect,&diff_join);

      sort_edge (&diff_join);

      //std::cout << std::endl;
      intersect.clear();
      // T n path
      isect_edge(&paths_e.at(i), &Q.at(j), &intersect);

      if(intersect.empty()) {
        R_star.push_back(diff_join);
      }
      else {
        //std::cout << "Diff_join:" << diff_join.size() << std::endl;
        R.push_back(diff_join);
      }
    } // for_2 ends
/*
    std::cout << "R-size " << R.size() << std::endl;
    for (int iii = 0; iii < R.size(); iii++) {
      std::cout << "R-Cycle " << iii << ", " << "size: " << R[iii].size() 
                << std::endl;
      for (it_e = R[iii].begin(); 
           it_e != R[iii].end(); ++it_e) {
        std::cout << it_e->first << " ";
      }
      std::cout << std::endl;
    }
*/
    std::vector< std::vector<vd_Edge> >::iterator it_R;
    it_R = std::unique (R.begin(), R.end()); 
    R.resize(std::distance(R.begin(),it_R));

    it_R = std::unique (R_star.begin(), R_star.end()); 
    R_star.resize(std::distance(R_star.begin(),it_R));
/*
    std::cout << "R-size " << R.size() << std::endl;
    for (int iii = 0; iii < R.size(); iii++) {
      std::cout << "R-Cycle " << iii << ", " << "size: " << R[iii].size() 
                << std::endl;
      for (it_e = R[iii].begin(); 
           it_e != R[iii].end(); ++it_e) {
        std::cout << it_e->first << " ";
      }
      std::cout << std::endl;
    }
*/
    // Track if v is moved to R_star.
    std::vector<bool> in_star(0);

    in_star.resize(R.size());
    std::fill(in_star.begin(), in_star.end(), false); 

    for (int u = 0; u < R.size(); u++) {
      self_int.clear();
      for (int s_it = 0; s_it < R.at(u).size(); s_it++) {
        self_int[R.at(u).at(s_it).first.source()] = self_int[R.at(u).at(s_it).first.source()] + 1;
        self_int[R.at(u).at(s_it).first.target()] = self_int[R.at(u).at(s_it).first.target()] + 1;
        if(self_int[R.at(u).at(s_it).first.source()] > 2 or self_int[R.at(u).at(s_it).first.target()] > 2)
          in_star[u] = true;
      }

      for (int v = u+1; v < R.size(); v++) {
        if (in_star[u]) {
          v = R.size();
        }
        else if (in_star[v]) {
          v++;
        }
        else {
          intersect.clear();
          isect_edge(&R[u], &R[v], &intersect);
/*
          std::cout << "R[" << u << "]:" << std::endl;
          for(int i_u = 0; i_u < R[u].size(); i_u++)
            std::cout << R[u].at(i_u).first << " ";
          std::cout << std::endl;

          std::cout << "R[" << v << "]:" << std::endl;
          for(int i_v = 0; i_v < R[v].size(); i_v++)
            std::cout << R[v].at(i_v).first << " ";
          std::cout << std::endl;

          std::cout << "R[" << u << "], R[" << v << "], int " 
                    << R[u].size() << " " << R[v].size() << " " << intersect.size() 
                    << std::endl;
*/
          if (intersect.size() == R.at(u).size() ) {
            //std::cout << "v " << v << std::endl;
            in_star[v] = true;
          }
          else if (intersect.size() == R.at(v).size() ) {
            //std::cout << "u " << u << std::endl;
            in_star[u] = true;
            v = R.size();
          }
        }
      }// for_3 ends
      // TODO this is a hackish attempt at removing edge lists not composed of 
      // a single circuit. Rather than revising this section, I'm suppressing the
      // symptom at the moment, as this only creates a problem of extra edgelists
      // that are not circuits and they are resolved by this step.
      if(!in_star[u]) {
        if(chk_discon_circ(R.at(u)))
        in_star[u] = true;
      }

    }// for_2 ends

    for (int u = R.size()-1; u > -1; u--) {
      if (in_star[u]) {
        //std::cout << "Popping back.." << u << std::endl;
        R_star.push_back(R.at(u));
        std::swap(R.at(u), R.back());
        R.pop_back();
      }
    }// for_2 ends

    //std::cout << "S-size " << S.size() << ", R-size " << R.size() 
    //<< std::endl;

    // Sort sets
    for (int iii = 0; iii < S.size(); iii++)
      sort_edge(&S.at(iii));
    for (int iii = 0; iii < R.size(); iii++)
      sort_edge(&R.at(iii));
    sort_edge_sets (&S);
    sort_edge_sets (&R);

    S_temp.reserve(S.size()+R.size());
    merge_edge_set(&S,&R,&S_temp);

    for (int iii = 0; iii < S_temp.size(); iii++)
      sort_edge(&S_temp.at(iii));
    sort_edge_sets (&S_temp);

    for (int iii = 0; iii < Q.size(); iii++)
      sort_edge(&Q.at(iii));
    for (int iii = 0; iii < R_star.size(); iii++)
      sort_edge(&R_star.at(iii));
    sort_edge_sets (&Q);
    sort_edge_sets (&R_star);

    Q_temp.reserve(Q.size()+R_star.size());

    merge_edge_set(&Q,&R_star,&Q_temp);

    Q.clear();
    S.clear();

    std::vector<std::vector<vd_Edge> > set_temp
            (0, std::vector<vd_Edge>(0) );
    set_temp.clear();
    set_temp.push_back(paths_e.at(i));

    //Q.reserve(Q_temp.size()+R.size());
    //merge_edge_set(&Q_temp,&R,&Q);

    for (int iii = 0; iii < Q_temp.size(); iii++)
      sort_edge(&Q_temp.at(iii));
    sort_edge_sets (&Q_temp);

    Q_temp2.reserve(Q_temp.size()+R.size());

    merge_edge_set(&Q_temp,&R,&Q_temp2);

    for (int iii = 0; iii < Q_temp2.size(); iii++)
      sort_edge(&Q_temp2.at(iii));
    sort_edge_sets (&Q_temp2);

    Q.reserve(Q_temp2.size()+set_temp.size());

    merge_edge_set(&Q_temp2,&set_temp,&Q);

    for (int iii = 0; iii < Q.size(); iii++)
      sort_edge(&Q.at(iii));
    sort_edge_sets (&Q);

    S.reserve(S_temp.size()+set_temp.size());

    merge_edge_set(&S_temp,&set_temp,&S);

    for (int iii = 0; iii < S.size(); iii++)
      sort_edge(&S.at(iii));
    sort_edge_sets (&S);

    //std::cout << "S-size " << S.size() << ", Q-size " << Q.size() 
    //<<std::endl;

    it_R = std::unique (S.begin(), S.end()); 
    S.resize(std::distance(S.begin(),it_R));
    it_R = std::unique (Q.begin(), Q.end()); 
    Q.resize(std::distance(Q.begin(),it_R));

    //std::cout << "S-size " << S.size() << ", Q-size " << Q.size() 
    //<<std::endl;
  } // for_1 ends
/*
  for (int u = 0; u < S.size(); u++) {
    for (int v = u+1; v < S.size(); v++) {
      intersect.clear();
      isect_edge(&S[u], &S[v], &intersect);
      std::cout << "S[" << u << "]:" << std::endl;
      for(int i_u = 0; i_u < S.at(u).size(); i_u++)
        std::cout << S.at(u).at(i_u).first << " ";
      std::cout << std::endl;

      std::cout << "R[" << v << "]:" << std::endl;
      for(int i_v = 0; i_v < S.at(v).size(); i_v++)
        std::cout << S.at(v).at(i_v).first << " ";
      std::cout << std::endl;

      std::cout << "S[" << u << "], S[" << v << "], int " 
                << S.at(u).size() << " " << S.at(v).size() << " " << intersect.size() 
                << std::endl;

      if (intersect.size() == S.at(u).size() or 
          intersect.size() == S.at(v).size() ) {
        std::cout << "Intersect same size!!!"  << std::endl;
      }
    }// for_3 ends
    assert(!chk_discon_circ(S.at(u)));
  }  
*/

} // find_gibbs ends


// Given a graph, first find the elementary circuits by Paton algorithm.
// Using the gibbs algorithm, find all circuits as edge lists. 
// Convert the edge lists into vertex lists.
void PatonFinder::find_circuits() {
  find_paton();
  vert_2_edge();

/*
  std::cout << "Paton found " << paths_e.size() << std::endl;

  for (int i = 0; i < paths_e.size(); i++) {
    for (int j = 0; j < paths_e[i].size(); j++) {
      std::cout << paths_e[i][j].first << " ["
                << paths_e[i][j].second << "] ";
    }
    std::cout << std::endl;
  }
*/

  //std::set<vd_Edge >::iterator it_e;
  find_gibbs();
}

// Return a pointer to the list with the circuit and subgraph tuples.
//std::vector< crc > *
std::vector< std::pair< std::vector<vd_Edge>, s_graph > > *
          PatonFinder::get_circuits() {
  return &circ_tup;
}

void PatonFinder::clear() {
    g.clear();
    tree.clear();

    T.clear();
    V.clear();

    for(int i = 0; i < paths.size(); i++) {
      paths.at(i).clear();
    }
    paths.clear();

    for(int i = 0; i < paths_e.size(); i++) {
      paths_e.at(i).clear();
    }
    paths_e.clear();

    for(int i = 0; i < paths_eset.size(); i++) {
      paths_eset.at(i).clear();
    }
    paths_eset.clear();

    for(int i = 0; i < S.size(); i++) {
      S.at(i).clear();
    }
    S.clear();

    for(int i = 0; i < disj_frag.size(); i++) {
      disj_frag.at(i).clear();
    }
    disj_frag.clear();

    for(int i = 0; i < circ_tup.size(); i++) {
      circ_tup.at(i).first.clear();
      circ_tup.at(i).second.first.clear();
      circ_tup.at(i).second.second.clear();
    }
    circ_tup.clear();

    index.clear();
}

PatonFinder::~PatonFinder() {
  clear();

}

// ---------------------------------------------------------
// Functions for doing graph based search on a cellbase object.
// ---------------------------------------------------------

// ---------------------------------------------------------
// ---------------------------------------------------------
// cell_graph object:
// ---------------------------------------------------------
// Default constructor private.
// ---------------------------------------------------------
cell_graph::cell_graph() :
    g(), g_12(),
    PF(),
    cell_ctr(0), cell_type{}, cell_lbls(0),
    circ_tup(0),
    path(0),

    path_act(0), cell_start(0), cell_end(0),
    path_curr(0), path_curr_id(0),
    c2_curr(0), c3_curr(0),

    visited{}, cell1_flag(false), cell1_act{}, c1_int{},
    ext_sz(0), calc_corner(false),
    cell_ext(false), cell_ext_corner(false),
    cells3(), cells2(), cells1() {
  circ_class[0] = 0;
  circ_class[1] = 0;
}


// ---------------------------------------------------------
// Call the routines of the PatonFinder.
// ---------------------------------------------------------
// Convert the id of a 23 adjacency graph 2cell to a 12 adjacency graph.
int cell_graph::conv_c2_12(int c2) {
  if(c2 > -1)
    return c2 - cells3.conn.size();
  else
    return c2 + ext_sz;
}

// Given two 2cells in 2/3-adj graph, check if they share a common 1cell adjacency. 
// Used within disjoint cell grouping inside PatonFinder.
bool cell_graph::comp_conn_12(int c21, int c22) {
  // Both cells should be 2cells in the 2/3-adj graph.
  if ((c21 > -1 and c21 < (int)(cells3.conn.size())) or (c22 > -1 and c22 < (int)(cells3.conn.size())))
    return false;
  else if ((c21 < 0 and c21 > - ext_sz - 1) or (c22 < 0 and c22 > - ext_sz - 1))
    return false;

  // The indices on the 1/2-adj graph.
  int c21_12 = conv_c2_12(c21);
  int c22_12 = conv_c2_12(c22);

  std::set<int> c12_int{};

  std::set<int>::const_iterator c21_st;
  std::set<int>::const_iterator c21_end;

  std::set<int>::const_iterator c22_st;
  std::set<int>::const_iterator c22_end;

  std::tie(c21_st, c21_end) = g_12.adjacent_vertices(
                                       c21_12);
  std::tie(c22_st, c22_end) = g_12.adjacent_vertices(
                                       c22_12);

  c12_int.clear();
  std::set_intersection(c21_st,c21_end, c22_st, c22_end,
                          std::inserter(c12_int,c12_int.begin()));
  if (c12_int.size() > 0) {
    //std::cout << "2c" << cells2.conn.at(c21_12)+1 
    //  << " and 2c" << cells2.conn.at(c22_12)+1 << "connected over ";
    //for(c21_st = c12_int.begin(); c21_st != c12_int.end(); c21_st++) 
    //  std::cout << " 1c" << cells1.conn[*c21_st-cells2.conn.size()]+1;
    //std::cout << std::endl;
    return true;
  }
  return false;
}

// Treat the ciruits and paths for graphs containing multiple exterior cells.
void cell_graph::treat_circuits() {
  //vd_Edge_desc null_edge(-1,-1);

  //reduce the edges

  for(int i = 0; i < circ_tup.size(); i++) {
    int skip = 0;
    for(int j = 0; j < circ_tup.at(i).first.size(); j++) {
      int n1 = circ_tup.at(i).first.at(j).first.source();
      int n2 = circ_tup.at(i).first.at(j).first.target();

      if(n1 < 0 and n2 < 0) {
        skip = skip + 1;
      }
      else {
        if(n1 < 0)
          n1 = -1;
        if(n2 < 0)
          n2 = -1;
        circ_tup.at(i).first.at(j).first = vd_Edge_desc(n1, n2);
        circ_tup.at(i).first.at(j-skip) = circ_tup.at(i).first.at(j);
      }
    }
    circ_tup.at(i).first.resize(circ_tup.at(i).first.size() - skip);
  }

  // Check the circuit for multiply repeated exterior. If so, empty the disj.
  // graphs. The circuit is going to be classified degenerate.
  for(int i = 0; i < circ_tup.size(); i++) {
    bool flip1 = false;
    bool flip2 = false;
    bool degen = false;
    bool all_virt = true;
    bool ext_circ = false;

    int ext_count = 0;
    //std::cout << "circ " << i << std::endl;
    for(int j = 0; j < circ_tup.at(i).first.size(); j++) {
      //std::cout << "edge " << j << std::endl;
      int n1 = circ_tup.at(i).first.at(j).first.source();
      int n2 = circ_tup.at(i).first.at(j).first.target();
      if(n1 < 0 or n2 < 0) {
        ext_count = ext_count + 1;
        //std::cout << "ext_count " << ext_count << std::endl;
        if(flip2 or ext_count > 2) {
          degen = true;
          j = circ_tup.at(i).first.size();
        }

        if(flip1)
          flip2 = true;
        else
          flip1 = true;
        if(flip1 or flip2)
          ext_circ = true;
      }
      else
        all_virt = false;
    }
    if(all_virt)
      degen = true;
    if(degen) {
      circ_tup.at(i).first.clear();
      circ_tup.at(i).second.first.clear();
      circ_tup.at(i).second.second.clear();
    }
    else {
      // Only one exterior 3c is allowed.
      bool ext3c = false;
      // If the exterior is on the circuit, remove from disjoint
      // graphs.
      if(ext_circ)
        ext3c = true;
      int skip = 0;
      for(int j = 0; j < circ_tup.at(i).second.first.size(); j++) {
        int n1 = circ_tup.at(i).second.first.at(j);
        // Circuits that detach 2cells are degenerate. 
        if((cell_type[n1] == -2) or (ext3c and cell_type[n1] == -3)) {
          skip = skip + 1;
        }
        else {
          if(cell_type[n1] == -3) {
            circ_tup.at(i).second.first.at(j) = -1;
            ext3c = true;
          }
          circ_tup.at(i).second.first.at(j-skip) =
                                 circ_tup.at(i).second.first.at(j);
        }
      }
      circ_tup.at(i).second.first.resize(
                        circ_tup.at(i).second.first.size()-skip);

      skip = 0;
      for(int j = 0; j < circ_tup.at(i).second.second.size(); j++) {
        int n1 = circ_tup.at(i).second.second.at(j);
        // Circuits that detach only virtual 2cells are degenerate. 
        if((cell_type[n1] == -2) or (ext3c and cell_type[n1] == -3)) {
          skip = skip + 1;
        }
        else {
          if(cell_type[n1] == -3) {
            circ_tup.at(i).second.second.at(j) = -1;
            ext3c = true;
          }
          circ_tup.at(i).second.second.at(j-skip) =
                                 circ_tup.at(i).second.second.at(j);
        }
      }
      circ_tup.at(i).second.second.resize(
                        circ_tup.at(i).second.second.size()-skip);

    }
    std::sort(circ_tup.at(i).second.first.begin(), 
            circ_tup.at(i).second.first.end());
    std::sort(circ_tup.at(i).second.second.begin(), 
            circ_tup.at(i).second.second.end());
  }

  for(int i = 0; i < circ_tup.size(); i++) {
    for(int j = i+1; j < circ_tup.size(); j++) {
      if(circ_tup_comp_disj(circ_tup.at(i), circ_tup.at(j))) {

        circ_tup.at(j).first.clear();
        circ_tup.at(j).second.first.clear();
        circ_tup.at(j).second.second.clear();
      }
    }
  }

}

void cell_graph::treat_paths_act() {

  int skip_out = 0;
  for(int i = 0; i < path.at(path_act).second.size(); i++) {
    bool flip1 = false;
    bool flip2 = false;

    bool degen = false;

    int skip_in = 0;
    if(path.at(path_act).second.at(i).first.at(0) < 0) {
      assert(path.at(path_act).second.at(i).first.at(0) == -1);
      flip1 = true;
    }
    for(int j = 1; j < path.at(path_act).second.at(i).first.size(); j++) {
      //std::cout << "f1 " << flip1 << " f2 " << flip2 << std::endl;
      //std::cout << "s_in " << skip_in << " s_out " << skip_out << std::endl;
      //std::cout << path.at(path_act).second.at(i).first.at(j) << std::endl;
      if(flip2) {
        if(path.at(path_act).second.at(i).first.at(j) < 0) {
          j = path.at(path_act).second.at(i).first.size();
          skip_out = skip_out + 1;
          degen = true;
        }
        else
          path.at(path_act).second.at(i).first.at(j - skip_in) = 
          path.at(path_act).second.at(i).first.at(j);
      }
      else if(flip1) {
        if(path.at(path_act).second.at(i).first.at(j) < 0) {
          skip_in = skip_in + 1;
        }
        else {
          flip2 = true;
          path.at(path_act).second.at(i).first.at(j - skip_in) = 
          path.at(path_act).second.at(i).first.at(j);
        }
      }
      else {
        if(path.at(path_act).second.at(i).first.at(j) < 0) {
          path.at(path_act).second.at(i).first.at(j) = -1;
          flip1 = true;
        }
        else {
          path.at(path_act).second.at(i).first.at(j - skip_in) = 
          path.at(path_act).second.at(i).first.at(j);
        }
      }
    }
    path.at(path_act).second.at(i).first.resize(
                path.at(path_act).second.at(i).first.size() - skip_in);

    //std::cout << "After " << path.at(path_act).second.at(i).first.size() 
    //          << " skip " << skip_in
    //          << std::endl;
    //for(int j = 0; j < path.at(path_act).second.at(i).first.size(); j++) {
    //  std::cout << path.at(path_act).second.at(i).first.at(j) << " ";
    //}
    //std::cout << std::endl;

    if(!degen and path.at(path_act).second.at(i).first.size() == 1 and
      path.at(path_act).second.at(i).first.at(0) == -1) {
      skip_out = skip_out + 1;
      degen = true;
    }
    //std::cout << "degen " << degen << std::endl;
    if(!degen) {
      //std::cout << "Replacing " << i-skip_out << " " << i << std::endl;
      path.at(path_act).second.at(i-skip_out).first =
                                  path.at(path_act).second.at(i).first;
      path.at(path_act).second.at(i-skip_out).second = 
                                  path.at(path_act).second.at(i).second;
    }
  }
  if(skip_out > 0) {
    //std::cout << "Replacing after " 
    //          << path.at(path_act).second.size() - skip_out 
    //          << " New size " << path.at(path_act).second.size()-skip_out 
    //          << std::endl;
    for(int i = path.at(path_act).second.size() - skip_out; 
            i < path.at(path_act).second.size(); i++) {
      path.at(path_act).second.at(i).first.clear();
      path.at(path_act).second.at(i).second.clear();
    }
    path.at(path_act).second.resize(path.at(path_act).second.size()-skip_out);

  }

  // Also, remove duplicates.
  for(int i = 0; i < path.at(path_act).second.size(); i++) {
    skip_out = 0;
    for(int j = i + 1; j < path.at(path_act).second.size(); j++) {
      if(path.at(path_act).second.at(i).first ==
         path.at(path_act).second.at(j).first) {
        skip_out = skip_out + 1;
      }
      else {
        path.at(path_act).second.at(j-skip_out).first =
                                    path.at(path_act).second.at(j).first;
        path.at(path_act).second.at(j-skip_out).second = 
                                    path.at(path_act).second.at(j).second;

      }
    }
    for(int i = path.at(path_act).second.size() - skip_out; 
            i < path.at(path_act).second.size(); i++) {
      path.at(path_act).second.at(i).first.clear();
      path.at(path_act).second.at(i).second.clear();
    }
    path.at(path_act).second.resize(path.at(path_act).second.size()-skip_out);
  }

}

// Treat the graphs containing multiple exterior cells after detection of 
// circuits and paths.
void cell_graph::treat_graphs() {
  std::set<int>::const_iterator ai_st; 
  std::set<int>::const_iterator ai_end;

  std::set<int>::const_iterator af_st; 
  std::set<int>::const_iterator af_end;
  // g23.
  for (std::tie(ai_st, ai_end) = g.vertices(); 
                    ai_st != ai_end; ++ai_st) {
    int ci_id = *ai_st;

    bool all_clear = false;
    while(!all_clear) {
      all_clear = true;
      for (std::tie(af_st, af_end) = g.adjacent_vertices(ci_id); 
                        af_st != af_end; ++af_st) {
        int cf_id = *af_st;
        if(ci_id < 0 and cf_id < 0) {
          if(g.check_edge(ci_id, cf_id)) {
            g.remove_edge(ci_id, cf_id);
            all_clear = false;
            //af_st = af_end;
            //--af_st;
            break;
          }
        }
        else {
          if(ci_id < -1) {
            if(!g.check_edge(-1, cf_id))
              g.add_edge(-1, cf_id);
            if(g.check_edge(ci_id, cf_id)) {
              g.remove_edge(ci_id, cf_id);
              all_clear = false;
              //af_st = af_end;
              //--af_st;
              break;
            }
          }
          if(cf_id < -1) {
            if(!g.check_edge(-1, ci_id))
              g.add_edge(-1, ci_id);
            if(g.check_edge(ci_id, cf_id)) {
              g.remove_edge(ci_id, cf_id);
              all_clear = false;
              //af_st = af_end;
              //--af_st;
              break;
            }
          }
        }
      }
    }
  }

  bool all_clear = false;
  while(!all_clear) {
    all_clear = true;
    for (std::tie(ai_st, ai_end) = g.vertices(); ai_st != ai_end; ++ai_st) {
      int ci_id = *ai_st;
      if(ci_id < -1) {
        g.remove_vertex(ci_id);
        std::tie(ai_st, ai_end) = g.vertices();
        all_clear = false;
      }
    }
  }

  // g12.
  for (std::tie(ai_st, ai_end) = g_12.vertices(); 
                    ai_st != ai_end; ++ai_st) {
    int ci_id = *ai_st;
    for (std::tie(af_st, af_end) = g_12.adjacent_vertices(ci_id); 
                      af_st != af_end; ++af_st) {
      int cf_id = *af_st;
      if(ci_id < 0 or cf_id < 0) {        
        g_12.remove_edge(ci_id, cf_id);
      }
    }
  }
  all_clear = false;
  while(!all_clear) {
    all_clear = true;
    for (std::tie(ai_st, ai_end) = g_12.vertices(); ai_st != ai_end; ++ai_st) {
      int ci_id = *ai_st;
      if(ci_id < 0) {
        g_12.remove_vertex(ci_id);
        std::tie(ai_st, ai_end) = g_12.vertices();
        all_clear = false;
      }
    }
  }

}

void cell_graph::chk_discon_circs() {
  std::vector<vd_Edge> temp(0);
  for(int i = 0; i < circ_tup.size(); i++) {
    temp = circ_tup.at(i).first;
    assert(!chk_discon_circ(temp));
  }
}

// Find the circuits within the loaded graph. Get the disjoint graphs 
// separated by the circuits.
void cell_graph::find_circuits() {
  PF.reload_graph(&g);
  PF.find_circuits();
  PF.create_sub(this);

  std::vector< std::pair< std::vector<vd_Edge>, s_graph > >* pf_circ;
  pf_circ = PF.get_circuits();

  circ_tup.resize(pf_circ->size());
  for(int i = 0; i < pf_circ->size(); i++) {
    circ_tup.at(i).first = pf_circ->at(i).first;
    circ_tup.at(i).second.first = pf_circ->at(i).second.first;
    circ_tup.at(i).second.second = pf_circ->at(i).second.second;
  }

  if(calc_corner and cell_ext_corner)
    treat_circuits();
  // Every circuit on the circuit list should have a single circuit on it.
  chk_discon_circs();

  //print_circ();
  //std::cout << "Sorting the circuits" << std::endl;
  quicksort_circ_tup (0, circ_tup.size()-1);

  // Get the number of degenerate and non-degenerate circuits. Assume sorted.
  int i = 0;
  while(i < circ_tup.size() and cont_23cell (&circ_tup.at(i).second)) {
    i++;
  }

  circ_class[0] = i-1;

  while(i < circ_tup.size() and chk_nondegen (&circ_tup.at(i).second)) {
    i++;
  }
  circ_class[1] = i-1;

  //std::cout << "Number of non-degenerate cases: 3cell: " << circ_class[0]
  //<< " 2cell: " << circ_class[1] << std::endl;

  //print_circ();

  //std::ofstream fout2("output/23adj.dot");
  //write_graphviz(fout2, g, boost::make_label_writer(&cell_lbls[0]));

}

// Set the non-intersecting lists of all paths.
void cell_graph::find_nonintsct() {
  // Reserve the non-intersecting lists.
  for (int i = 0; i < path.at(path_act).second.size(); i++) {
    path.at(path_act).second.at(i).second.clear();
    path.at(path_act).second.at(i).second.reserve(path.at(path_act).second.size());
  }

  // For couples of paths, compare the visited cells. If there is a match, 
  // skip.
  // If no match were found between the paths, add to the list for both paths.
  // The non-intersecting path lists are sorted.
  for (int i = 0; i < path.at(path_act).second.size(); i++) {
    for (int j = i+1; j < path.at(path_act).second.size(); j++) {
      bool match = false;
      //std::cout << "p " << i << " " << j << std::endl;
      for (int ik = 1; ik < path.at(path_act).second.at(i).first.size()-1; ik++) {
        for (int jk = 1; jk < path.at(path_act).second.at(j).first.size()-1; jk++) {
          //std::cout << "\t"
          //          << path.at(path_act).second.at(i).first.at(ik) << " "
          //          << path.at(path_act).second.at(j).first.at(jk)
          //          << std::endl;
          if (path.at(path_act).second.at(i).first.at(ik) 
              == path.at(path_act).second.at(j).first.at(jk)) {
            match = true;
            ik = path.at(path_act).second.at(i).first.size(); 
            jk = path.at(path_act).second.at(j).first.size(); 
          }
        }
      }
      if (!match) {
        //std::cout << "no match, adding " << std::endl;
        path.at(path_act).second.at(i).second.push_back(j);
        path.at(path_act).second.at(j).second.push_back(i);
      }
    }
  }
/*
  for (int i = 0; i < path[path_act].second.size(); i++) {
    std::cout << "Cell couple " << i << std::endl;
    for (int j = 0; j < path[path_act].second[i].second.size(); j++) {
      std::cout << path[path_act].second[i].second[j] << ", ";
    }
    std::cout << std::endl;
  }
*/
}

// Check the path length and 1cell_flag. If the path length is longer than 3 
// (more than one 2cell) 1cell_flag will remain the same. 
// If 1cell_flag is true, than do a check over the 2cells, to find out 
// candidate 1cells.
bool cell_graph::chk_1cell() {
  cell1_act.clear();
  c1_int.clear();

  cell1_flag = false;

  // For path size 1, there are no 2cells on path, yet. So, clearing the 1cell 
  // list is enough.

  // If there is at least a 2cells on path, The 1cells touching the 0cell are 
  // to be added to the cell1_act list.
  // For the first two 2cell, check if they share 1cells. Set the 1cell_flag 
  // accordingly.
  if (path_curr.size() >= 1) {

    if (path_curr.size() >= 3) {

      std::set<int>::const_iterator c1_st;
      std::set<int>::const_iterator c1_end;
      std::tie(c1_st, c1_end) = g_12.adjacent_vertices(
                                              path_curr.at(1)-cells3.conn.size());

      for (; c1_st != c1_end; c1_st++) {
        cell1_act.insert(*c1_st);
      }

      // There are at least two 2cells, which should be checked for joint 1cell
      // adjacencies. If the path had been updated correctly, they should either
      // all joint at a number of 1cells or not join for subsequent 2cell 
      // couples.
      if (path_curr.size() >= 5) {
        std::tie(c1_st, c1_end) = g_12.adjacent_vertices(
                                              path_curr.at(3)-cells3.conn.size());

        c1_int.clear();
        std::set_intersection(cell1_act.begin(),cell1_act.end(),
                           c1_st, c1_end, 
                          std::inserter(c1_int,c1_int.begin()));

        // If the first two 2cells intersect, all other should intersect as well.
        if (c1_int.size() > 0) {
          cell1_act.clear();

          std::copy(c1_int.begin(), c1_int.end(),
            std::inserter(cell1_act,cell1_act.begin()));

          // Get the intersection of current 1cells and the 1cells of the next
          // 2cell on the path. 
          for(int i = 5; i < path_curr.size(); i=i+2) {
            std::tie(c1_st, c1_end) = g_12.adjacent_vertices(
                                              path_curr.at(i)-cells3.conn.size());

            c1_int.clear();
            std::set_intersection(cell1_act.begin(),cell1_act.end(),
                           c1_st, c1_end, 
                          std::inserter(c1_int,c1_int.begin()));

            assert(c1_int.size() > 0);

            cell1_act.clear();

            std::copy(c1_int.begin(), c1_int.end(),
              std::inserter(cell1_act,cell1_act.begin()));
          }

          cell1_flag = true;
        }

        // If the they don't, no 2cell may intersect. Store the 1cell adj of the
        // last 2cell.
        else {
          std::tie(c1_st, c1_end) = g_12.adjacent_vertices(
                           path_curr[path_curr.size()-2]-cells3.conn.size());

          cell1_act.clear();
          for (; c1_st != c1_end; c1_st++) {
            cell1_act.insert(*c1_st);
          }
          cell1_flag = false;
        }
      }

    }
  }

  return cell1_flag;
}

int cell_graph::get_circ_type(int c_in) {
  if(c_in > circ_class[1])
    return 2;
  else if(c_in > circ_class[0])
    return 1;

  return 0;
}

// Check if the current path has reached the end 3cell. If so, add it to the 
// path list.
bool cell_graph::try_path() {
  int c3_next_id = -1;
  if(c3_curr < 0)
    c3_next_id = -1;
  else
    c3_next_id = cells3.conn.at(c3_curr);
  //std::cout << std::endl;
  //std::cout << "End vertex: " << c3_next_id 
  //          << ", current vertex: "<< cell_end << std::endl;

  if (c3_next_id == cell_end) {
    s_graph s_curr;
    s_curr.first.reserve(path_curr.size()+2);
    for(int i = 0; i < path_curr.size(); i++) {
      s_curr.first.push_back(path_curr.at(i));
    }
/*
    std::cout << "Path found: type " << cell1_flag << std::endl;
    for(int i = 0; i < s_curr.first.size(); i++) {
      if (s_curr.first[i] < cells3.conn.size())
        std::cout << "3c" << cells3.conn.at(s_curr.first[i]) + 1 << ", ";
      else
        std::cout << "2c" 
              << cells2.conn.at(s_curr.first[i]-cells3.conn.size()) + 1 << ", ";
    }
    std::cout << std::endl;
*/
    path.at(path_act).second.push_back(s_curr);

    return true;
  } // If path end check.
  return false;
}

// Appends the current 2- and 3-cells to the current path.
void cell_graph::path_append() {
  path_curr.push_back(c2_curr);
  path_curr.push_back(c3_curr);
  visited[c2_curr] = true;
  visited[c3_curr] = true;
}

// Remove the last 2- and 3-cells from the current path.
void cell_graph::path_remove() {
  assert(path_curr.size() >= 3);
  path_curr.pop_back();
  path_curr.pop_back();
  visited[c2_curr] = false;
  visited[c3_curr] = false;
}

// Trace all possible paths by recursively visiting neighboring unvisited 
// vertices of the last vertex on path.
// Extend path, recursively call trace_path. When you finish a line, exit one 
// recursion. Remove the entended part. Do a 1cell check.
// After each iteration, check if the 2cell and the following 3cell are 
// appropriate(both unvisited, 3cell exists).
// If the current length is 1 simply add the 2cell and the following 3cell.
// If the current length is 3, in addition, do a 1cell check.
// For the following extensions, check 1cell_flag.
// If true, check if intersecting 1cells exist. If not, skip.
// Else, check if intersecting 1cells exist. If so, skip.
void cell_graph::trace_path() {
  // Get the adjacencies of the last 3cell on the path. Initially, this is one
  // the disconnected 3cells.
  int c3_last = path_curr.back();
  std::set<int>::const_iterator ai_st; 
  std::set<int>::const_iterator ai_end;
  std::tie(ai_st, ai_end) = g.adjacent_vertices(c3_last);

/*
  std::cout << c3_last << " 3c_curr";

  if(c3_last == -1)
    std::cout << "Ext" <<  std::endl;
  else
    std::cout << cells3.conn.at(c3_last)+1 <<  std::endl;
*/
  // Iterate over the neighboring 2cells.
  //std::cout << "Going over the neighbors " << std::endl;
  while (ai_st != ai_end) {
    c2_curr = *ai_st;
    //int c2_id = cells2.conn.at(c2_curr-cells3.conn.size());

    //std::cout << c2_curr << " 2c_curr"
    //        << c2_id+1 << "[" << visited[*ai_st] <<  "]" <<  std::endl;

    std::set<int>::const_iterator c2_st; 
    std::set<int>::const_iterator c2_end;
    std::tie(c2_st, c2_end) = g.adjacent_vertices(c2_curr);

    // Check if the current 2cell is used.
    if (!visited[c2_curr]) {
      int c2_neigh = std::distance(c2_st,c2_end);
      // Get the following 3cell by getting the 3cell neighbors of the 2cell.

      //std::cout << "Neighbors: " << c2_neigh << std::endl;

      //assert((c2_neigh == 1) or (c2_neigh == 2));
      // With the exterior added as a 3cell, all 2cells should have exactly  
      // two adjacencies. 
      assert((c2_neigh == 2));

      int c3_temp;
      if (*c2_st == c3_last) {
        c2_st++;
        assert(c3_last != *c2_st);
        c3_curr = *c2_st;
      }
      else {
        c3_curr = *c2_st;
      }

      c3_temp = c3_curr;
/*
        if(c3_temp == -1)
          std::cout << "Ext" <<  std::endl;
        else
          std::cout << c3_curr << " 3c_curr"
              << cells3.conn.at(c3_curr)+1 << "[" 
              << visited[*ai_st] <<  "]" <<  std::endl;
*/

      if (!visited[c3_curr]) {
        // 1cell adjacencies of the current 2cell:
        std::tie(c2_st, c2_end) = g.adjacent_vertices(c2_curr);

        path_append();
        // The length of path_curr is 3, so this configures the cell1_act
        // for a single 2cell.
        //chk_1cell();

        // Add to the path.
        if(try_path()) {
        }
        else {
          trace_path();
          c2_curr = *ai_st;
          c3_curr = c3_temp;
        }
        path_remove();
        //std::cout << "\\2c" << c2_id+1 << "\\3c" 
        //          << c3_id+1 << ", ";
      } // If 3cell not visited check. 
    } // if visited ends 
    ai_st++;
  } // while ends 

} // trace_paths ends

// TODO if the 0cell is on simulating cell boundary, check 3cells not connected
// to the exterior. These can be connected to the exterior by a 2cell insertion.
void cell_graph::find_paths() {
  path.clear();
  int cell3_nbr = cells3.conn.size();
  std::map<vd_Edge_desc, bool > c3_map{};
  std::cout << "Loading paths" << std::endl;
  if ( cell3_nbr > 0) {
    // TODO not a hard limit, and finding a hard limit shouldn't be possible.
    //if(calc_corner and cell_ext_corner)
    //  path.reserve((cell3_nbr+ext_sz) * (cell3_nbr+ext_sz-1) );
    if(cell_ext)
      path.reserve((cell3_nbr+1) * cell3_nbr );
    else
      path.reserve(cell3_nbr * (cell3_nbr-1) );

    for (int i = 0; i < cell3_nbr; i++) {
      std::set<int>::const_iterator ai_st; 
      std::set<int>::const_iterator ai_end;
      // Going over the 2cell adjacencies of the 3cell, get the other 3cells.
      for (std::tie(ai_st, ai_end) = g.adjacent_vertices(i); 
                        ai_st != ai_end; ++ai_st) {
        int n3 = g.get_node_other(*ai_st, i);
        if(n3 < 0)
          c3_map[vd_Edge_desc(n3, i)] = true;
        else
          c3_map[vd_Edge_desc(i, n3)] = true;
      }
      for (int j = i+1; j < cell3_nbr; j++) {
        //std::cout << std::endl;
        // If the 3cells are not connected by a 2cell, or if the 2cell is not
        // bounded by the current 0cell, these 3cells can be connected. 
        if (!c3_map[vd_Edge_desc(i, j)]) {
          cell_path cp_1;
          cell_cp coup(cells3.conn.at(i), cells3.conn.at(j));
          std::cout << "Couple " << cells3.conn.at(i) 
                    << " " <<cells3.conn.at(j) 
                    << " not connected." << std::endl;
          cp_1.first = coup;
          path.push_back(cp_1);
        }
      }
      // Also, check the exterior. If there is a bounding 2cell bounding 
      // only this 3cell, it is a candidate.
      // If the 0cell is connected to a side 1cell, it is still a candidate,
      // even if the exterior is connected.
      if(calc_corner and cell_ext_corner) {
        cell_path cp_1;
        cell_cp coup(-1, cells3.conn.at(i));
        std::cout << "Corner 0cell: Couple " << cells3.conn.at(i) + 1
                  << " " << -1
                  << " is considered not connected." << std::endl;
        cp_1.first = coup;
        path.push_back(cp_1);
      }
      else if(cell_ext) {
        if (!c3_map[vd_Edge_desc(-1, i)]) {
          cell_path cp_1;
          cell_cp coup(-1, cells3.conn.at(i));
          std::cout << "Couple " << cells3.conn.at(i) + 1
                    << " " << -1
                    << " not connected." << std::endl;
          cp_1.first = coup;
          path.push_back(cp_1);
        }
      }
    }
  }

  // We have the possible 3cell couples to connect by a 2cell generation.
  // Now, a DFS algorithm to find the eligible paths between the couples.

  // Longest possible path would at most include all cells once.
  path_curr.clear();
  if(cell_ext)
    path_curr.reserve(g.num_vertices());
  else
    path_curr.reserve(cells3.conn.size()+cells2.conn.size());

  if (path.size() > 0) {
    for (path_act = 0; path_act < path.size(); path_act++) {
      // Actual ids of the cells.
      cell_start = path.at(path_act).first.first;
      cell_end = path.at(path_act).first.second;
      // Local id of the current cell.
      int cell_curr = cells3.find_ent(cell_start);
      //std::cout << cell_start + 1 << " [" << cell_curr << "] " 
      //          << cell_end + 1 << std::endl; 
      path_curr.clear();
      path_curr.push_back(cell_curr);
      visited.clear();
      visited[cell_curr] = true;
      trace_path();
      visited[cell_curr] = false;
      if(calc_corner and cell_ext_corner)
        treat_paths_act();
//      if(calc_corner and cb->get_0c_corner_gmi(cell_ctr))
//        find_ext_nonintsct();
//      else
      find_nonintsct();
      //std::cout << path_act << "Path size: " 
      //          << path.at(path_act).second.size() << std::endl;

    }
  }

}

int cell_graph::get_0c() {
  return cell_ctr;
}

bool cell_graph::get_ext() {
  return cell_ext;
}

// Given path 2-cell ids, return slice 2-cells and their paths.
void cell_graph::find_slice(std::vector<
            std::pair<std::vector<int >,std::vector<int > > >* path_cells, 
    std::vector< std::pair< std::pair<int,int>, 
                 std::vector<std::vector<int > > > >* slice_cells) {
  // Burn paths. Go over 2-cells until finding one unburned. 
  // Add the first one to the on_fire list.
  // Put the 1-cell and 3-cell adjacencies of the 2-cell to the next on fire 
  // list.
  // Burn and clear the previous list.
  // Alternate until cell on fire is non-empty.
  std::vector<int> fire_2cell(0);
  std::vector<int> fire_3cell(0);
  std::vector<int> fire_1cell(0);

  std::map<int, int> path_id{};
  std::map<int, int> slice_id_23{};
  std::map<int, int> slice_id_12{};
  std::map<int, bool> burned_23{};
  std::map<int, bool> burned_12{};

  std::set<int>::const_iterator ai_st; 
  std::set<int>::const_iterator ai_end;

  //for(int i = 0; i < cells2.conn.size(); i++) {
  //  std::cout << " 2c " << cells2.conn.at(i);
  //}
  //for(int i = 0; i < cells3.conn.size(); i++) {
  //  std::cout << " 3c " << cells3.conn.at(i);
  //}

  //std::cout << "Path ";
  for(int i = 0; i < path_cells->size(); i++) {
    std::cout << i;
    for(int j = 0; j < path_cells->at(i).second.size(); j++) {
      //c_base id
      int c_id = path_cells->at(i).second.at(j);
      //std::cout << " 2c " << c_id + 1;
      //local id
      c_id = cells2.find_ent(c_id);
      assert(c_id > -1);
      burned_12[c_id] = true;
      burned_23[c_id+cells3.conn.size()] = true;
      path_id[c_id+cells3.conn.size()] = i+1;
      //std::cout << " [2c " << cell_lbls[c_id+cells3.conn.size()] << "] ";
    }
    std::cout << std::endl;
    for(int j = 0; j < path_cells->at(i).first.size(); j++) {
      //c_base id
      int c_id = path_cells->at(i).first.at(j);
      std::cout << " 3c " << c_id + 1;
      if(c_id == -1) {
        //std::cout << " Ext ";
      }
      else {
        //local id
        c_id = cells3.find_ent(c_id);
        //std::cout << " [3c " << cell_lbls[c_id] << "] ";
      }

      path_id[c_id] = i+1;
      burned_23[c_id] = true;
    }
    std::cout << std::endl;
  }

  for(int i = 0; i < slice_cells->size(); i++) {
    for(int j = 0; j < slice_cells->at(i).second.size(); j++) {
      slice_cells->at(i).second.at(j).clear();
      slice_cells->at(i).second.at(j).clear();
      slice_cells->at(i).second.at(j).clear();
    }
    slice_cells->at(i).second.clear();
  }
  slice_cells->clear();

  int slice = 1;
  bool all_burned = false;
  while(all_burned == false) {
    //std::cout << "All burned " << all_burned << std::endl;

    fire_2cell.clear();
    fire_3cell.clear();
    fire_1cell.clear();

    fire_2cell.reserve(cells2.conn.size());
    fire_3cell.reserve(cells3.conn.size());
    fire_1cell.reserve(cells1.conn.size());
    // Starting from a 2-cell, the final cell to burn would be a 1-cell.
    // First do 2-cells.
    all_burned = true;
    for(int i = 0; i < cells2.conn.size(); i++) {
      int c_id = i+cells3.conn.size();

      if(!burned_23[c_id]) {
        burned_23[c_id] = true;
        slice_id_23[c_id] = slice;

        burned_12[c_id-cells3.conn.size()] = true;
        slice_id_12[c_id-cells3.conn.size()] = slice;

        fire_2cell.push_back(c_id);
        all_burned = false;
        i = cells2.conn.size();
      }
    }
    if(all_burned == false) {
      bool all_ablaze = false;

      while((!all_ablaze)) {
/*
        std::cout << "2c_sz " << fire_2cell.size()
                << " empty " << fire_2cell.empty()
                << " 3c_sz " << fire_3cell.size()
                << " empty " << fire_3cell.empty()
                << " 1c_sz " << fire_1cell.size()
                << " empty " << fire_1cell.empty()
                << " !all_ablaze " << !all_ablaze
                << std::endl;
*/
        // First loop over 2cell adjacencies.
        all_ablaze = true;
        while(!fire_2cell.empty()) {
          int cell2_curr = fire_2cell.back();
          fire_2cell.pop_back();
          //std::cout << "popping 2c" 
          //          << cells2.conn.at(cell2_curr-cells3.conn.size()) + 1 << " ";
          for (std::tie(ai_st, ai_end) = g.adjacent_vertices(cell2_curr); 
                        ai_st != ai_end; ++ai_st) {
            if (!burned_23[*ai_st]) {
              fire_3cell.push_back(*ai_st);
              burned_23[*ai_st] = true;
              slice_id_23[*ai_st] = slice;

              all_ablaze = false;
            }
          }
          for (std::tie(ai_st, ai_end) = 
                     g_12.adjacent_vertices(cell2_curr-cells3.conn.size()); 
                        ai_st != ai_end; ++ai_st) {
            if (!burned_12[*ai_st]) {
              fire_1cell.push_back((*ai_st));
              burned_12[*ai_st] = true;
              slice_id_12[*ai_st] = slice;

              all_ablaze = false;
            }
          }
        }

        // Next loop over 1- and 3-cell adjacencies.
        while(!fire_3cell.empty()) {
          int cell3_curr = fire_3cell.back();
          fire_3cell.pop_back();

          //std::cout << "popping 3c" << cell3_curr << " ";
          for (std::tie(ai_st, ai_end) = g.adjacent_vertices(cell3_curr); 
                        ai_st != ai_end; ++ai_st) {
            if (!burned_23[*ai_st]) {
              fire_2cell.push_back(*ai_st);

              burned_23[*ai_st] = true;
              slice_id_23[*ai_st] = slice;

              burned_12[(*ai_st) - cells3.conn.size()] = true;
              slice_id_12[(*ai_st) - cells3.conn.size()] = slice;

              all_ablaze = false;
            }
          }
        }

        while(!fire_1cell.empty() ) {
          int cell1_curr = fire_1cell.back();
          fire_1cell.pop_back();
          //std::cout << "popping 1c" 
          //          << cells1.conn.at(cell1_curr-cells2.conn.size()) + 1 << " ";
          for (std::tie(ai_st, ai_end) = 
                     g_12.adjacent_vertices(cell1_curr); 
                        ai_st != ai_end; ++ai_st) {
            if (!burned_12[*ai_st]) {
              fire_2cell.push_back((*ai_st) + cells3.conn.size());

              burned_23[(*ai_st) + cells3.conn.size()] = true;
              slice_id_23[(*ai_st) + cells3.conn.size()] = slice;

              burned_12[*ai_st] = true;
              slice_id_12[*ai_st] = slice;

              all_ablaze = false;
            }
          }
        }

      } // While fire_empty ends
      slice++;
      //std::cout << "Slice next " << slice << std::endl;
    } // If ends
    //std::cout << "All burned " << all_burned << std::endl;
  }
  for(std::tie(ai_st, ai_end) = g.vertices(); ai_st != ai_end; ++ai_st) {
    if(*ai_st == -1)
      std::cout << "Ext ";
    else if(*ai_st < (int)(cells3.conn.size()))
      std::cout << "3c" << cell_lbls[*ai_st];
    else
      std::cout << "2c" << cell_lbls[*ai_st];

    std::cout << " burn: " << burned_23[*ai_st]
              << " slice: " << slice_id_23[*ai_st]
              << " path: " << path_id[*ai_st]

              << std::endl;
  }
  for(std::tie(ai_st, ai_end) = g_12.vertices(); ai_st != ai_end; ++ai_st) {
    if(*ai_st < (int)(cells2.conn.size()))
      std::cout << "2c" << cell_lbls[(*ai_st)+cells3.conn.size()];
    else
      std::cout << "1c" << cells1.conn.at((*ai_st)-cells2.conn.size());

    std::cout << " burn: " << burned_12[*ai_st]
              << " slice: " << slice_id_12[*ai_st]
              << std::endl;
  }

  slice_cells->resize(slice-1);

  for(int i = 0; i < slice_cells->size(); i++) {
    slice_cells->at(i).second.resize(3);
    slice_cells->at(i).second.at(0).reserve(cells1.conn.size());
    slice_cells->at(i).second.at(1).reserve(cells2.conn.size());
    slice_cells->at(i).second.at(2).reserve(cells3.conn.size());
  }

  for(std::tie(ai_st, ai_end) = g_12.vertices(); ai_st != ai_end; ++ai_st) {
    if(*ai_st >= (int)(cells2.conn.size())) {
      if(slice_id_12[*ai_st] == 0) {
        std::cout << "cell(1) " << *ai_st << " " << (*ai_st)-cells2.conn.size() 
                  << " is the only cell on its slice. Adding new slice for it."
                  << std::endl;

        slice_cells->push_back(std::make_pair(std::make_pair(0,0), 
                 std::vector<std::vector<int > > (3, std::vector<int >(0)) ));
        //slice_cells->back().second.at(0).push_back((*ai_st)-cells2.conn.size());
        slice_id_12[*ai_st] = slice_cells->size();
      }
    }
  }

  // Obtain the path adjacencies of the slices by looking at the slice 1cell
  // adjacencies.
  std::cout << "Obtain the path adjacencies of slices." << std::endl;

  for(std::tie(ai_st, ai_end) = g_12.vertices(); ai_st != ai_end; ++ai_st) {

    if(*ai_st >= (int)(cells2.conn.size())) {
      assert(slice_id_12[*ai_st] > 0);
      std::cout << "cell(1) " << *ai_st << " " << (*ai_st)-cells2.conn.size() 
                << " slice " << slice_id_12[*ai_st] << std::endl;

      int slice_curr = slice_id_12[*ai_st] - 1;
      slice_cells->at(slice_curr).second.at(0).
                        push_back((*ai_st)-cells2.conn.size());

      std::set<int>::const_iterator a1_st; 
      std::set<int>::const_iterator a1_end;

      for(std::tie(a1_st, a1_end) = g_12.adjacent_vertices(*ai_st); 
                  a1_st != a1_end; ++a1_st) {

        // The 2cell adjacent to the slice 1cell is on a path.
        if(*a1_st < (int)(cells2.conn.size()) and 
            slice_id_23[*a1_st + cells3.conn.size()] == 0) {
          //std::cout << *a1_st << " 2c" << cells2.conn[*a1_st] << " ";
          int path_curr = path_id[(*a1_st) + cells3.conn.size()];
          //std::cout << " path" << path_curr << std::endl;
          if(path_curr != 0) {
            if(slice_cells->at(slice_curr).first.first == 0) {
              slice_cells->at(slice_curr).first.first = path_curr;
            }
            else if (slice_cells->at(slice_curr).first.first != path_curr and
                      slice_cells->at(slice_curr).first.second == 0) {
              slice_cells->at(slice_curr).first.second = path_curr;
            }
          }
        }
      }
    }
  }

  for(std::tie(ai_st, ai_end) = g.vertices(); ai_st != ai_end; ++ai_st) {
    if(slice_id_23[*ai_st] > 0) {
      if(*ai_st >= (int)(cells3.conn.size())) {
        std::cout << "cell(2) " << *ai_st << " " << (*ai_st)-cells3.conn.size() 
                  << " slice " << slice_id_23[*ai_st] << std::endl;

        int slice_curr = slice_id_23[*ai_st] - 1;
        slice_cells->at(slice_curr).second.at(1).
                          push_back((*ai_st)-cells3.conn.size());
      }
      else {
        std::cout << "cell(3) " << *ai_st << " " << (*ai_st) 
                  << " slice " << slice_id_23[*ai_st] << std::endl;
        int slice_curr = slice_id_23[*ai_st] - 1;
        slice_cells->at(slice_curr).second.at(2).
                          push_back((*ai_st));
      }
    }
  }

  for(int i = 0; i < slice_cells->size(); i++) {
    std::cout << "Slice " << i
              << " path1: " << slice_cells->at(i).first.first
              << " path2: " << slice_cells->at(i).first.second
              << std::endl;
  }

}

// Reorder the circuits and subgraphs, based on the content of the subgraphs.
// If both of the subgraphs contain no cell, it is degenerate.
// If one of the subgraphs contains no cell, it is semi-degenerate.
// If the disjoint graphs both contain at least a 2cell, 

// ---------------------------------------------------------
// Custom label functions to differentiate cells of different dimensions.
// ---------------------------------------------------------


// ---------------------------------------------------------
// Load the graph different topologies.
// ---------------------------------------------------------

// Load the 2- and 3-cells from the cellbase. Boost graph handles both kind 
// as vertices. The vertices are ordered such that the 3cells come first, with 
// the same order in the cell-complex and these are followed by the 2cells. 
// The distinction between the 2- and 3-cells is stored in a vector array, 
// cell_type. 
// If exterior exists, it will be the last entry in cells->at(2) list.
void cell_graph::load_cells(int tag_0cell, 
                       std::vector<std::vector<int > >* cells, 
                       std::vector<std::vector<std::vector<int > > >* adj) {
  clear();

  cell_ctr = tag_0cell - 1;

  int cell3_nbr = 0;

  ext_sz = 0;
  if (cb->get_cell_ext(0, cell_ctr)) {
    cell_ext = true;
    if (calc_corner and cb->get_0c_corner(cell_ctr)) {
      cell_ext_corner = true;
    }
  }

  if (calc_corner and cb->get_0c_corner(cell_ctr)) {
    for(int i = 0; i < cells->at(2).size(); i++) {
      if(cells->at(2).at(i) > -1)
        cell3_nbr = cell3_nbr + 1;
      else{
        i = cells->at(2).size();
      }
    }
    ext_sz = cells->at(2).size() - cell3_nbr;
    cells3.conn.reserve(cell3_nbr);
    cells2.conn.reserve(cells->at(1).size() - ext_sz);
  }
  else if(cell_ext) {
    ext_sz = 1;
    cell3_nbr = cells->at(2).size() - ext_sz;
  }
  else {
    cell3_nbr = cells->at(2).size();
  }
  std::copy(cells->at(2).begin(), cells->at(2).end() - ext_sz,
                   std::inserter(cells3.conn, cells3.conn.begin()));
  if (calc_corner and cb->get_0c_corner(cell_ctr)) {
    std::copy(cells->at(1).begin(), cells->at(1).end() - ext_sz,
                   std::inserter(cells2.conn, cells2.conn.begin()));
    std::copy(cells->at(0).begin(), cells->at(0).end() - 1,
                   std::inserter(cells1.conn, cells1.conn.begin()));
  }
  else {
    std::copy(cells->at(1).begin(), cells->at(1).end(),
                   std::inserter(cells2.conn, cells2.conn.begin()));
    std::copy(cells->at(0).begin(), cells->at(0).end(),
                   std::inserter(cells1.conn, cells1.conn.begin()));
  }

  int cell2_nbr = cells2.conn.size();
  int cell1_nbr = cells1.conn.size();
  if (calc_corner and cb->get_0c_corner(cell_ctr)) {
    cell1_nbr = cell1_nbr - 1;
  }

  //cell1_act.reserve(cell1_nbr);
  //cell_lbls.reserve(cell3_nbr+cell2_nbr);
  cell_lbls.resize(cell3_nbr+cell2_nbr);

  for (int i = 0; i < cells->at(2).size(); i++) {
    if(cells->at(2).at(i) < 0)
      cell_type[cells->at(2).at(i)] = -3;
    else
      cell_type[cells->at(2).at(i)] = 3;
    if(cells->at(2).at(i) > 0)
      cell_lbls.at(i) = cells->at(2).at(i);
  }

  for (int i = 0; i < cells->at(1).size(); i++) {
    if(cells->at(1).at(i) < 0)
      cell_type[cells->at(1).at(i)-ext_sz] = -2;
    else
      cell_type[cells->at(1).at(i)+cell3_nbr] = 2;
    if(cells->at(1).at(i) > 0)
      cell_lbls.at(i+cell3_nbr) = cells->at(1).at(i);
  }
  g.reserve_v(cell3_nbr+cell2_nbr+ext_sz*2);

  for (int i = 0; i < cells->at(2).size(); i++) {
    if(cells->at(2).at(i) < 0)
      g.add_vertex(cells->at(2).at(i));
    else
      g.add_vertex(i);
  }

  for (int i = 0; i < cells->at(1).size(); i++) {
    if(cells->at(1).at(i) < 0)
      g.add_vertex(-ext_sz+cells->at(1).at(i));
    else
      g.add_vertex(cell3_nbr+i);
  }

  for (int i = 0; i < cells->at(2).size(); i++) {
    int c3_curr = cells->at(2).at(i);
    for (int j = 0; j < adj->at(1).at(i).size(); j++) {
      int c2_curr = adj->at(1).at(i).at(j);

      if(c3_curr < 0) {
        if(c2_curr < 0) {
          g.add_edge(c3_curr, -ext_sz+c2_curr);
        }
        else {
          g.add_edge(c3_curr, cell3_nbr+c2_curr);
        }
      }
      else {
        if(adj->at(1).at(i).at(j) < 0) {
          g.add_edge(i, -ext_sz+c2_curr);
        }
        else {
          g.add_edge(i, cell3_nbr+c2_curr);
        }
      }
    }
  }

  // g_12 contains the 1cell-2cell adjacency map. First cell2_nbr indices 
  // correspond to the 2cells, and the subsequent ones correspond to 1cells.
  if (calc_corner and cb->get_0c_corner(cell_ctr)) {
    // Add the exterior 1c
    g_12.add_vertex(-ext_sz);
  }
  for (int i = 0; i < cells1.conn.size(); i++) {
    g_12.add_vertex(cell2_nbr + i);
  }

  for (int i = 0; i < cells->at(1).size(); i++) {
    g_12.add_vertex(cells->at(1).at(i));
    for (int j = 0; j < adj->at(0).at(i).size(); j++) {
      int c2_curr = cells->at(1).at(i);
      int c1_curr = adj->at(0).at(i).at(j);
      if(c1_curr < 0) {
        g_12.add_edge(c2_curr, -ext_sz+c1_curr);
      }
      else {
        g_12.add_edge(c2_curr, cell2_nbr+c1_curr);
      }
    }
  }

  print_graph("output/cells.dot");

  find_circuits();
  //print_circ();
  find_paths();
  treat_graphs();
  //std::cout << "Number of cell couples, 2-3adj: " << path.size() << std::endl;
}

void cell_graph::get_23adj(int tag_0cell) {

  clear();

  cell_ctr = tag_0cell;

  cb->get_conn_dim(3, 0, tag_0cell, &cells3);
  cb->get_conn_dim(2, 0, tag_0cell, &cells2);
  cb->get_conn_dim(1, 0, tag_0cell, &cells1);

  int cell3_nbr = cells3.conn.size();
  int cell2_nbr = cells2.conn.size();
  int cell1_nbr = cells1.conn.size();

  //cell1_act.reserve(cell1_nbr);
  //cell_lbls.reserve(cell3_nbr+cell2_nbr);
  cell_lbls.resize(cell3_nbr+cell2_nbr);
  //cell_type.reserve(cell3_nbr+cell2_nbr);
  //cell_type.resize(cell3_nbr+cell2_nbr);

  for (int i = 0; i < cells3.conn.size(); i++) {
    cell_type[i] = 3;
    //cell_lbls[i] = "3c";
    //std::string str3 = std::to_string(cells3.conn.at(i)+1);
    //cell_lbls[i].append(str3);
    cell_lbls[i] = cells3.conn.at(i)+1;
  }

  g.reserve_v(cell3_nbr+cell2_nbr+1);

  for (int i = 0; i < cells2.conn.size(); i++) {
    cell_type[cell3_nbr+i] = 2;
    //cell_lbls[cell3_nbr+i] = "2c";
    //std::string str2 = std::to_string(cells2.conn.at(i)+1);
    //cell_lbls[cell3_nbr+i].append(str2);
    cell_lbls[cell3_nbr+i] = cells2.conn.at(i)+1;
  }

  for (int i = 0; i < cell3_nbr+cell2_nbr; i++)
    g.add_vertex(i);

  for (int i = 0; i < cells3.conn.size(); i++) {
    struct ent_conn e_con;
    cb->get_conn(3, cells3.conn.at(i), &e_con);
    for (int j = 0; j < cells2.conn.size(); j++) {
      if (e_con.chk_ent(cells2.conn.at(j)))
        g.add_edge(i, cell3_nbr+j);
    }
  }

  // g_12 contains the 1cell-2cell adjacency map. First cell2_nbr indices 
  // correspond to the 2cells, and the subsequent ones correspond to 1cells.
  for (int i = 0; i < cell2_nbr+cell1_nbr; i++)
    g_12.add_vertex(i);

  for (int i = 0; i < cells2.conn.size(); i++) {
    struct ent_conn e_con;
    cb->get_conn(2, cells2.conn.at(i), &e_con);
    for (int j = 0; j < cells1.conn.size(); j++) {
      if (e_con.chk_ent(cells1.conn.at(j)))
        g_12.add_edge(i, cell2_nbr+j);
    }
  }

  if (cb->get_cell_ext(0, cell_ctr)) {
    cell_ext = true;
    ext_sz = 0;
    if (calc_corner and cb->get_0c_corner(cell_ctr)) {
      cell_ext_corner = true;

      struct ent_conn e1;
      std::map<int, int> c3_mark{};
      std::map<int, bool> c2_flag{};
      std::map<int, int> c2_mark{};

      int c2_v_sz = 0;
      // For each exterior c2, add a virtual external c3
      for (int i = 0; i < cells2.conn.size(); i++) {
        if (cb->get_cell_ext(2, cells2.conn.at(i))) {
          //std::cout << "One bounded 3cell, Ext must exist" << std::endl; 
          ext_sz = ext_sz + 1;
          g.add_vertex(-ext_sz);
          cell_type[-ext_sz] = -3;
          g.add_edge(-ext_sz, cell3_nbr+i);
          c3_mark[cells2.conn.at(i)] = -ext_sz;
        }
      }

      // Going over exterior c2, for each exterior c1 bounding the c2, add a 
      // virtual external c2. Connect it to the exterior c3.
      for (int i = 0; i < cells2.conn.size(); i++) {
        if (cb->get_cell_ext(2, cells2.conn.at(i))) {
          cb->get_conn(2, cells2.conn.at(i), &e1);
          cb->rem_conn(0, cell_ctr, 1, &e1);
          int ext3_curr = c3_mark[cells2.conn.at(i)];
          for (int j = 0; j < e1.conn.size(); j++) {
            if(!c2_flag[e1.conn.at(j)]) {
              c2_v_sz = c2_v_sz + 1;
              c2_mark[e1.conn.at(j)] = -c2_v_sz;
              c2_flag[e1.conn.at(j)] = true;
              g.add_vertex(-ext_sz + c2_mark[e1.conn.at(j)]);
              cell_type[-ext_sz + c2_mark[e1.conn.at(j)]] = -2;
              g_12.add_vertex(-c2_v_sz);

              std::vector<int>::iterator e1_it = std::find(e1.conn.begin(),
                                                          e1.conn.end(), 
                                                          e1.conn.at(j));
              assert(e1_it != e1.conn.end());
              g_12.add_edge(cell2_nbr + std::distance(e1_it, e1.conn.end()), 
                                                                -c2_v_sz);
            }
            g.add_edge(-ext_sz + c2_mark[e1.conn.at(j)], ext3_curr);
          }
        }
      }

      g_12.add_vertex(-c2_v_sz - 1);
      for (int i = 0; i < c2_v_sz; i++) {
        g_12.add_edge(-i-1, -c2_v_sz - 1);
      }

    }
    else {
      ext_sz = 1;
      g.add_vertex(-ext_sz);
      cell_type[-ext_sz] = -3;
      for (int i = 0; i < cells2.conn.size(); i++) {
        if (cb->get_cell_ext(2, cells2.conn.at(i))) {
          //std::cout << "One bounded 3cell, Ext must exist" << std::endl; 
          g.add_edge(-ext_sz, cell3_nbr+i);
        }
      }
    }

    //std::cout << "Boundary 0cell, external grain added" << std::endl;
  }

  print_graph("output/graph.dot");

  find_circuits();
  //print_circ();
  find_paths();
  treat_graphs();
  //std::cout << "Number of cell couples, 2-3adj: " << path.size() << std::endl;

}

// Trial functions:

int cell_graph::node_nbr(int n, int i, int j) {
  return j*n+i;
}

void cell_graph::load_tiles(int n, int m) {
  clear();

  std::vector<std::string> lbls(0, std::string(""));

  lbls.reserve(n*m);
  lbls.resize(n*m);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      int index = node_nbr(n,i,j);
      g.add_vertex(index);
      lbls[index] = std::to_string(i);
      std::string str2 = std::to_string(j);
      lbls[index].append("-");
      lbls[index].append(str2);
    }
  }

  for (int i = 0; i < n-1; i++) {
    for (int j = 0; j < m; j++) {
      g.add_edge(node_nbr(n, i, j), node_nbr(n, i+1, j));
    }
  }

  // Collect vertical edges
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m-1; j++) {
      g.add_edge(node_nbr(n, i, j), node_nbr(n, i, j+1));
    }
  }

  find_circuits(); 
}

void cell_graph::remove_edge(int v1, int v2) {
  g.remove_edge(v1, v2);
}

// Test function to print the current graph.
void cell_graph::print_adj() {
  g.print_adj();

}

void cell_graph::print_graph(const char* dotname) {
  g.print_graph(dotname);
}

// Test function to print the current graph. pathcells is a list of 2-cells andd 
// 3-cells associated with individual paths. It used to highlight edges.
void cell_graph::print_graph_p(const char* dotname, 
     std::vector<std::pair<std::vector<int >, std::vector<int > > >* path_cells) {
  std::ofstream fout(dotname);
  std::map<vd_Edge_desc, bool> hl_map{};
  std::map<vd_Edge_desc, int> path_map{};
  std::map<int, std::string> color_map{};
  // Simple way of a colormap. It is not a good colormap.
  // Count the number of groups of 3 paths.
  // For groups of 3 paths, loop over the three and assign the current cardinal
  // color the value 1-step*loop, where loop is the id of the current 3 paths.
  // Assign the next cardinal color the value step*loop.
  // step = 1/loop_max + 1;
  assert(path_cells != NULL);
  int cycles = path_cells->size()/3 + 1;
  double step = 1./cycles;
  int c3_nbr = cells3.conn.size();

  for (int i = 0; i < path_cells->size(); i++) {
    int loop_curr = i / 3;
    int cycle_curr = (i % 3);
    double shift = loop_curr * step;
    double val1 = 1 - shift;
    int i2 = (i+1)%3;

    std::stringstream ss;
    ss << " [color=\"";
    for (int j = 0; j < 3; j++) {
      if(j == 0)
        ss << "green";
      else if(j == 1)
        ss << "red";
      else
        ss << "blue";

      if(j == cycle_curr)
        ss << ";" << val1;
      else if((j + 1)%3 == cycle_curr)
        ss << ";" << shift;
      else
        ss << ";" << 0;
      if(j < 2)
        ss << ":";
    }
    ss << "\"]";
    std::string tmp = ss.str();
    color_map[i] = tmp;
  }

  if(path_cells->size() == 1) {
    for (int k = 0; k < path_cells->at(0).first.size() - 1; k++) {
      int c3_curr = path_cells->at(0).first.at(k);
      int c2_curr = path_cells->at(0).second.at(k) + c3_nbr;
      vd_Edge_desc edge_curr(c3_curr, c2_curr);
      path_map[edge_curr] = 0;
      hl_map[edge_curr] = true;

      c2_curr = path_cells->at(0).second.at(k+1) + c3_nbr;
      vd_Edge_desc edge_curr2(c3_curr, c2_curr);
      path_map[edge_curr2] = 0;
      hl_map[edge_curr2] = true;
    }
    int c3_curr = path_cells->at(0).first.back();
    int c2_curr = path_cells->at(0).second.back() + c3_nbr;
    vd_Edge_desc edge_curr(c3_curr, c2_curr);
    path_map[edge_curr] = 0;
    hl_map[edge_curr] = true;

    c2_curr = path_cells->at(0).second.at(0) + c3_nbr;
    vd_Edge_desc edge_curr2(c3_curr, c2_curr);
    path_map[edge_curr2] = 0;
    hl_map[edge_curr2] = true;
  }
  else {
    for (int i = 0; i < path_cells->size(); i++) {
      for (int k = 0; k < path_cells->at(i).first.size() - 1; k++) {
        int c3_curr = path_cells->at(i).first.at(k);
        int c2_curr = path_cells->at(i).second.at(k) + c3_nbr;
        vd_Edge_desc edge_curr(c3_curr, c2_curr);
        path_map[edge_curr] = i;
        hl_map[edge_curr] = true;

        c3_curr = path_cells->at(i).first.at(k+1);
        vd_Edge_desc edge_curr2(c3_curr, c2_curr);
        path_map[edge_curr2] = i;
        hl_map[edge_curr2] = true;
      }
    }
  }
  g.print_graph_p(hl_map, color_map, path_map, dotname);
}


//void cell_graph::print_graph_lbl(const char* dotname = "output/graph.dot") {
//  g.print_graph(&cell_lbls, dotname);
//}

// Print the current graph using labels. 
void cell_graph::print_graph_lbl(const char* dotname) {
  //assert(cell_lbls.size() == g.num_vertices() or
  //       cell_lbls.size()+1 == g.num_vertices());
  //assert(cell_type.size() == g.num_vertices() or
  //       cell_type.size()+1 == g.num_vertices());

  std::multiset<vd_Edge_desc >::iterator it;
  std::multiset<vd_Edge_desc >::iterator end;

  //std::tie(it, end) = g.edge_it();

  //std::cout << "lbl: " << cell_lbls.size() 
  //          << " type: " << cell_type.size() 
  //          << " num_vert: " << g.num_vertices()
  //          << std::endl;
  //for (; it != end; ++it) {
  //    std::cout << it->source() << "--" << it->target() << std::endl;
  //}
  std::tie(it, end) = g.edge_it();

  std::ofstream fout(dotname);

  fout << "graph A {\n";
  fout << "  rankdir=LR\n"
    << "ratio=\"fill\"\n"
    << "edge[style=\"bold\"]\n" << "node[shape=\"circle\"]\n";

  std::set<int>::const_iterator ai_st; 
  std::set<int>::const_iterator ai_end;

  for(std::tie(ai_st, ai_end) = g.vertices(); ai_st != ai_end; ++ai_st) {
    if (cell_type[*ai_st] == -3) {
      fout << "\""  << "3Ext" << -*ai_st  
           << "\"" << "[shape=\"circle\"]\n";
    }
    else if (cell_type[*ai_st] == -2) {
      fout << "\""  << "2Ext" << -*ai_st << "\"" 
           << "[style=filled fillcolor=blue width=0.3 shape=\"square\"]\n";
    }
    else if (cell_type[*ai_st] == 3) {
      fout << "\""  << cell_type[*ai_st] << "c" << cell_lbls[*ai_st]  
           << "\"" << "[shape=\"circle\"]\n";
    }
    else {
      assert(cell_type[*ai_st] == 2);
      fout << "\""  << cell_type[*ai_st] << "c" << cell_lbls[*ai_st]  << "\"" 
           << "[style=filled fillcolor=blue width=0.3 shape=\"square\"]\n";
    }
  }

  if (cell_ext) {

    for (; it != end; ++it) {
      //std::cout << it->source() << "--" << it->target() << std::endl;

      fout << "\"";

      if (cell_type[it->source()] == -3) {
        fout << "3Ext" << std::abs(it->source()) << "\"";
      }
      else if (cell_type[it->source()] == -2) {
        fout << "2Ext" << std::abs(it->source()) << "\"";
      }
      else {
        assert(cell_type[it->source()] == 2 or cell_type[it->source()] == 3);
        fout << cell_type[it->source()] << "c"
             << cell_lbls.at(it->source()) << "\"";
      }

      fout << " -- ";

      if (cell_type[it->target()] == -3) {
        fout << "\"" << "3Ext" << std::abs(it->target()) << "\"";
      }
      else if (cell_type[it->target()] == -2) {
        fout << "\"" << "2Ext" << std::abs(it->target()) << "\"";
      }
      else {
        assert(cell_type[it->target()] == 2 or cell_type[it->target()] == 3);
        fout << "\"" << cell_type[it->target()] << "c"
             << cell_lbls.at(it->target()) << "\"";
      }

      fout << "\n";
    }

  }

  else {
    for (; it != end; ++it) {
      fout << "\"" << cell_type[it->source()] << "c" 
             << cell_lbls.at(it->source()) << "\"" << " -- " 
           << "\"" << cell_type[it->target()] << "c" 
             << cell_lbls.at(it->target()) << "\"" << "\n";
    }
  }

  fout << "}\n";
  fout.close();

  //std::ofstream fout2("output/graph2.dot");
  //write_graphviz(fout2, g);
}

int cell_graph::g23type(int node_id) {
  return cell_type[node_id];
}

int cell_graph::g12type(int node_id) {
  return cell_type[node_id];
}

int cell_graph::g23id(int node_id) {
  if(cell_type[node_id])
    return cell_type[node_id];
  else
    return cell_type[node_id];
}

int cell_graph::g12id(int node_id) {
  return cell_type[node_id];
}

// ---------------------------------------------------------
// Return the pointer to the circuits and the disjoint graphs.
// ---------------------------------------------------------

// Return a pointer to the circuits.
std::vector< crc >* cell_graph::get_circuits() {
  return &circ_tup;
}
// Return a pointer to the paths.
std::vector< cell_path >* cell_graph::get_paths() {
  //std::cout << "Path size " << path.size() << std::endl;
  return &path;
}

void cell_graph::get_path_ngon(int path_id, ngon_gmi* ng) {
  ngon_tracer_g n_tr(this, ng);
  n_tr.get_path_ngon(path_id);

}

void cell_graph::get_path_ngon_gmi(int path_id, ngon_gmi* ng) {
  ngon_tracer_g n_tr(this, ng);
  n_tr.get_path_ngon_gmi(path_id);
}

// Print the circuits and the associated subgraphs.
// Edges are actually sorted such that U is always a 3cell and V is always a 
// 2cell. TODO
void cell_graph::print_circ() {

  std::cout << "0cell" << cell_ctr + 1 << std::endl;
  std::cout << " " << circ_tup.size() << std::endl;
  for (int i = 0; i < circ_tup.size(); i++) {
    //std::cout << i << "th path: " << std::endl;
    std::cout << i << "th path: " << circ_tup.at(i).first.size() << std::endl;
    for (int j = 0; j < circ_tup.at(i).first.size(); j++) {
      int U = circ_tup.at(i).first.at(j).first.source();
      int V = circ_tup.at(i).first.at(j).first.target();

      if (U < (int)(cells3.conn.size())) {
        if (U == -1)
          std::cout << "[Ext" << ", ";
        else
          std::cout << "[3c" << cells3.conn.at(U) + 1 << ", ";
      }
      else
        std::cout << "[2c" 
              << cells2.conn.at(U-cells3.conn.size()) + 1 << ", ";

      if (V < (int)(cells3.conn.size())) {
        if (V == -1)
          std::cout << "[Ext" << ", ";
        else
          std::cout << "3c" << cells3.conn.at(V) + 1 << "], ";
      }
      else
        std::cout << "2c" 
              << cells2.conn.at(V-cells3.conn.size()) + 1 << "], ";

      //std::cout << circ_tup.at(i).first[j].second << "] ";
    }
    std::cout << std::endl;

    std::cout << "Subgraph 1_tup: " << std::endl;
    std::cout << "sz: " << circ_tup.at(i).second.first.size() << std::endl;
    for (int j = 0; j < circ_tup.at(i).second.first.size(); j++) {
      int cell_id = circ_tup.at(i).second.first[j];
      if (cell_id < (int)(cells3.conn.size())) {
        if (cell_id == -1)
          std::cout << "Ext" << ", ";
        else
          std::cout << "3c" << cells3.conn.at(cell_id) + 1 << ", ";
      }
      else
        std::cout << "2c" 
              << cells2.conn.at(cell_id-cells3.conn.size()) + 1 << ", ";
    }
    std::cout << std::endl;

    std::cout << "Subgraph 2_tup: " << std::endl;
    std::cout << "sz: " << circ_tup.at(i).second.second.size() << std::endl;
    for (int j = 0; j < circ_tup.at(i).second.second.size(); j++) {
      int cell_id = circ_tup.at(i).second.second[j];
      if (cell_id < (int)(cells3.conn.size())) {
        if (cell_id == -1)
          std::cout << "Ext" << ", ";
        else
          std::cout << "3c" << cells3.conn.at(cell_id) + 1 << ", ";
      }
      else
        std::cout << "2c" 
              << cells2.conn.at(cell_id-cells3.conn.size()) + 1 << ", ";
    }
    std::cout << std::endl;
  }

}


// ---------------------------------------------------------
// Circuit classification and reordering.
// ---------------------------------------------------------

// Check both disjoint graphs, if both of them contain a 3cell, return true.
bool cell_graph::cont_3cell_and (const s_graph* circ_graph) { 
  return (cont_3cell(&circ_graph->first) and cont_3cell(&circ_graph->second));
}

bool cell_graph::cont_3cell_or (const s_graph* circ_graph) { 
  return (cont_3cell(&circ_graph->first) or cont_3cell(&circ_graph->second));
}

// Check if both disjoint graphs contain 2cells and at least one contains a 
// 3cell.
bool cell_graph::cont_23cell (const s_graph* circ_graph) { 
  if (chk_nondegen(circ_graph)) {
    bool b1 = cont_3cell(&circ_graph->first);
    bool b2 = cont_3cell(&circ_graph->second);
    return (b1 or b2);
  }
  return false;
}
//bool cell_graph::cont_23cell (const s_graph* circ_graph) { 
//  if (chk_nondegen(circ_graph))
//    return (cont_3cell(&circ_graph->first) or cont_3cell(&circ_graph->second));
//  return false;
//}

// Check the disjoint graph, if it contains at least a 3cell, return true.
bool cell_graph::cont_3cell (const std::vector<int >* disj_graph) { 
  //std::cout << "Sz " << disj_graph->size() << " " << cells3.conn.size() << std::endl;
  for (int i = 0; i < disj_graph->size(); i++) {
    //std::cout << "\t" << disj_graph->at(i) << std::endl;
    //if (disj_graph->at(i) < (int)(cells3.conn.size())) {
    if (std::abs(cell_type[disj_graph->at(i)]) == 3) {
      return true;
    }
  }
  return false;
}

// Check both disjoint graphs, if one of them is non-empty, return true.
bool cell_graph::chk_nondegen (const s_graph* circ_graph) { 
  return ((!circ_graph->first.empty() and !circ_graph->second.empty()));
}

// Check both disjoint graphs, if one of them is non-empty, return true.
bool cell_graph::chk_nondegen (const std::vector<int >* disj_graph) { 
  return (!disj_graph->empty() );
}

// Used in weak ordering. Returns false, if the first circuit is degenerate, 
// i.e. contains at least an empty disjoint graph, the second circuit contains a
// 3cell or the first doesn't and if both circuits are equivalent.
bool cell_graph::circ_tup_comp (const crc& circ_1, const crc& circ_2) { 

  // Tuples with disjoint graphs at least one containing 3cells should come 
  // first.
  // If one of the disjoint graphs is empty, it should come last.
  //return e1.first < e2.first;
  bool c3c_1 = cont_23cell (&circ_1.second );
  bool c3c_2 = cont_23cell (&circ_2.second );
  bool c2c_1 = chk_nondegen (&circ_1.second );
  bool c2c_2 = chk_nondegen (&circ_2.second );

  // ? - 2+
  if (c2c_2) {

    // 2+ - 2+
    if (c2c_1) {

      // 2+ - 3
      if (c3c_2) {
        // 3 - 3
        if (c3c_1) {
          return false;
        }
        // 2 - 3
        else {
          return false;
        }
      }
      // 3 - 2
      else if (c3c_1) {
        return true;
      }
      // 2 - 2
      else {
        return false;
      }
    }

    // 0 - 2
    else {
      return false;
    }
  }

  // ? - 0
  else {
    return true;
  }
}

//
bool cell_graph::circ_tup_comp_disj (const crc& circ_1, const crc& circ_2) { 

  // Skip degenerate circuits.
  if(circ_1.second.first.size() == 0 or circ_1.second.second.size() == 0 or
     circ_2.second.first.size() == 0 or circ_2.second.second.size() == 0)
    return false;
  // Compare the disjoint graph: 
  // 1st disjoint of the 1st circuit with the 1st disjoint of the 2nd circuit
  // 1st disjoint of the 1st circuit with the 2nd disjoint of the 2nd circuit
  bool try11 = true;
  bool try12 = true;
  if(circ_1.second.first.size() == circ_2.second.first.size()) {
    for(int i = 0; i < circ_1.second.first.size(); i++) {
      if(circ_1.second.first.at(i) != circ_2.second.first.at(i)) {
        try11 = false;
        i = circ_1.second.first.size();
      }
    }
    try12 = false;
  }
  else {
    try11 = false;
  }
  if(!try11 and circ_1.second.first.size() == circ_2.second.second.size()) {
    for(int i = 0; i < circ_1.second.first.size(); i++) {
      if(circ_1.second.first.at(i) != circ_2.second.second.at(i)) {
        try12 = false;
        i = circ_1.second.first.size();
      }
    }
  }
  else {
    try12 = false;
  }

  if(!try11 and !try12)
    return false;

  if(try11) {
    if(circ_1.second.second.size() != circ_2.second.second.size())
      return false;
    for(int i = 0; i < circ_1.second.second.size(); i++) {
      if(circ_1.second.second.at(i) != circ_2.second.second.at(i)) {
        return false;
      }
    }
  }
  else {
    if(circ_1.second.second.size() != circ_2.second.first.size());
      return false;
    for(int i = 0; i < circ_1.second.second.size(); i++) {
      if(circ_1.second.second.at(i) != circ_2.second.first.at(i)) {
        return false;
      }
    }
  }
  return true;

}

// Quick Sort Functions for Ascending Order 
// (2 Functions) 
void cell_graph::quicksort_circ_tup(int left, int right) {
  // left > right
  if (left < right) {
    int part = partition(left, right);
    quicksort_circ_tup(left, part - 1);
    quicksort_circ_tup(part + 1, right);
  }
}
//Function to determine the partitions
// partitions the array and returns the middle subscript
int cell_graph::partition(int left, int right) {
  crc pivot = circ_tup[right];
  // move the mid point value to the front.
  int i = left-1;
  int j = left;
  for (; j < right; j++) {
    if(!circ_tup_comp(pivot, circ_tup.at(j))) {
      i++;
      std::swap(circ_tup.at(i), circ_tup.at(j));
    }
  }
  std::swap(circ_tup.at(i+1),circ_tup.at(right));
  return i + 1;
}

// Circuit classification based on the cell dimensions on the disjoint graphs
// generated by removing the circuit from the cell complex.
// crc_both_3: Both disjoint graphs contain at least a 3cell.
// crc_2: Both disjoint graphs contain at least a 2cell or 3cell.
// crc_degen: One disjoint graph contains no cell.
enum crc_type {crc_both_3, crc_2, crc_degen, 
               END};

// Stores the circuit type to a vector that stores the information. 
// TODO Not a good application.
void cell_graph::get_circ_class(int cell0, std::vector<int>* pos) {
  pos->at(2*cell0) = circ_class[0];
  pos->at(2*cell0+1) = circ_class[1];
}

int cell_graph::get_3cell_sz() {
  return cells3.conn.size();
}

int cell_graph::get_2cell_sz() {
  return cells2.conn.size();
}

int cell_graph::get_1cell_sz() {
  return cells1.conn.size();
}

void cell_graph::set_calc_corner(bool fix) {
  calc_corner = fix;
}

// ---------------------------------------------------------
// Reset the cell_graph except for the c_base connection.
// ---------------------------------------------------------
void cell_graph::clear() {

  cells3.conn.clear();
  cells2.conn.clear();
  cells1.conn.clear();

  cell_type.clear();
  cell_lbls.clear();
  
  g.clear();
  g_12.clear();

  path.clear();
  path_curr.clear();

  cell1_act.clear();
  c1_int.clear();

  visited.clear();

  for(int i = 0; i < circ_tup.size(); i++) {
    circ_tup.at(i).first.clear();
    circ_tup.at(i).second.first.clear();
    circ_tup.at(i).second.second.clear();
  }
  circ_tup.clear();

  cell_ext = false;
  cell_ext_corner = false;
  //clear_ext_map();
}

void cell_graph::load_cb(struct cell_base* c_base, bool corner_flag) {
  clear();
  cb = c_base;

  calc_corner = corner_flag;
}

// ---------------------------------------------------------
// Constructor.
// ---------------------------------------------------------
//cell_graph::cell_graph() {
cell_graph::cell_graph(struct cell_base* c_base, bool corner_flag) :
    g(), g_12(),
    PF(),
    cell_ctr(0), cell_type{}, cell_lbls(0),
    circ_tup(0), path(0),

    path_act(0), cell_start(0), cell_end(0),
    path_curr(0), path_curr_id(0),
    c2_curr(0), c3_curr(0),

    visited{}, cell1_flag(false), cell1_act(), c1_int(),
    ext_sz(0), calc_corner(false),
    cell_ext(false), cell_ext_corner(false),
    cells3(), cells2(), cells1() {

  circ_class[0] = 0;
  circ_class[1] = 0;

  //clear();
  cb = c_base;

  calc_corner = corner_flag;

  PF.reload_graph(&g);
}

cell_graph::~cell_graph() {
  clear();
  //delete PF;
}

void cell_ins_chk::update_path(int cell0) {
  std::cout << "Updating 0cell " << cell0+1 << std::endl; 
  //std::cout << "0cell " << cell0+1 << " nbr of 0cells " 
  //          << cb->get_sz(0) << " circ_tup.size " << circ_tup.size()
  //          << std::endl;
  assert(cb->get_sz(0) > cell0);
  assert(cb->get_sz(0) == circ_tup.size());

  clear_path(cell0);

  if (!cb->is_free(0, cell0)) {
    std::cout << "Loading the 0cell..." << cell0+1 << std::endl;

    //std::cout << "Number of couples before: " << cg->get_paths()->size() 
    //    << std::endl;

    cg.get_23adj(cell0);

    //std::cout << "Cell graph obtained..." << std::endl;
    cb->get_conn_dim(3, 0, cell0, &cells3.at(cell0));
    cb->get_conn_dim(2, 0, cell0, &cells2.at(cell0));
    cb->get_conn_dim(1, 0, cell0, &cells1.at(cell0));

    std::vector< crc > * circ = cg.get_circuits();

    // TODO copy the source and target entities. Than find the unique ones.
    //std::cout << "Copying the circuit cells..." << std::endl;
    circ_tup.at(cell0).clear();
    circ_tup.at(cell0).resize(circ->size());
    for (int i = 0; i < circ->size(); i++) {
      //std::cout << "Circuit " << i << std::endl;
      // The circuit entities, to be taken from vd_Edge objects.
      std::map<int, int> map_1{};
      std::map<int, int> map_2{};

      circ_tup.at(cell0).at(i).first.first.reserve(circ->at(i).first.size());
      circ_tup.at(cell0).at(i).first.second.reserve(circ->at(i).first.size());

      for (int j = 0; j < circ->at(i).first.size(); j++) {
        int U = circ->at(i).first.at(j).first.source();
        int V = circ->at(i).first.at(j).first.target();
        map_1[U] = 0;
        map_1[V] = 0;
        map_2[U] = 0;
        map_2[V] = 0;
      }

      for (int j = 0; j < circ->at(i).first.size(); j++) {
        int U = circ->at(i).first.at(j).first.source();
        int V = circ->at(i).first.at(j).first.target();

        //std::cout << "U " << U << " V " << V << std::endl;
        //if (U == -1)
        //  std::cout << "3c" << "Ext ";
        //else
        //  std::cout << "3c" << cells3[cell0].conn.at(U) + 1;

        //std::cout << " 2c" << cells2[cell0].conn.at(V-cells3[cell0].conn.size()) + 1 
        //          << std::endl;

        //std::cout << "map_1[U] " << map_1[U] << " map_1[V] " << map_1[V] 
        //          << " map_2[U] " << map_2[U] << " map_2[V] " << map_2[V] 
        //          << std::endl;

        if(map_1[U] == 0) {
          map_1[U] = V+1;
        }
        // There can be at most two connections. Both connections are filled 
        // after visiting both neighboring cells.
        else {
          assert(map_2[U] == 0);
          map_2[U] = V+1;
        }
        if(map_1[V] == 0) {
          if(U == -1)
            map_1[V] = -1;
          else
            map_1[V] = U+1;
        }
        else {
          assert(map_2[V] == 0);
          if(U == -1)
            map_2[V] = -1;
          else
            map_2[V] = U+1;
        }
      }
      if(circ->at(i).first.size() > 0) {
        int last = circ->at(i).first.at(0).first.target();
        int next;
        if(map_1[last] == -1)
          next = -1;
        else
          next = map_1[last]-1;

        circ_tup.at(cell0).at(i).first.second.push_back(last);
        //std::cout << "First " << last << std::endl;
        // Assume circuit is closed.
        while(next != circ->at(i).first.at(0).first.target()) {
          //std::cout << "Last " << last << "Next " << next << std::endl;
          if(next < (int)(cells3.at(cell0).conn.size())) {
            circ_tup.at(cell0).at(i).first.first.push_back(next);
          }
          else {
            circ_tup.at(cell0).at(i).first.second.push_back(next);
          }
          //std::cout << "map_1[" << next << "] " << map_1[next] 
          //          << "map_2[" << next << "] " << map_2[next] 
          //          << std::endl;

          int next_curr;
          if(map_1[next] == -1)
            next_curr = -1;
          else
            next_curr = map_1[next]-1;

          if(last == next_curr) {
            last = next;
            if(map_2[next] == -1)
              next = -1;
            else
              next = map_2[next] - 1;
          }
          else {
            last = next;
            if(map_1[next] == -1)
              next = -1;
            else
              next = map_1[next] - 1;
          }

        }

        //std::cout << "The circuit cells: " << std::endl;

        //for (int j = 0; j < circ_tup[cell0][i].first.first.size(); j++) {
        //  std::cout << "id: " << circ_tup[cell0][i].first.first[j] << "[";
        //  if(circ_tup[cell0][i].first.first[j] == -1)
        //    std::cout << "Ext] ";
        //  else {
        //    int c3_curr = cells3[cell0].conn.at(
        //                      circ_tup[cell0][i].first.first[j]);
        //    std::cout << "3c" << c3_curr+1 << "] ";
        //  }
        //}
        //for (int j = 0; j < circ_tup[cell0][i].first.second.size(); j++) {
        //  std::cout << "id: " << circ_tup[cell0][i].first.second[j] << "[";
        //  int c2_curr = cells2[cell0].conn.at(
        //     circ_tup[cell0][i].first.second[j)-cells3[cell0].conn.size()];
        //  std::cout << "2c" << c2_curr+1 << "] ";
        //}
        //std::cout << std::endl;
      }
      // The disjoint graph entities, to be copied.
      circ_tup.at(cell0).at(i).second = circ->at(i).second;
    }

    std::vector< cell_path >* paths = cg.get_paths();

    //path[cell0];
    path.at(cell0).clear();
    path.at(cell0).resize(paths->size());

    //std::cout << "Number of couples after: " << paths->size() << std::endl;
    for (int i = 0; i < paths->size(); i++) {
      path.at(cell0).at(i).first = paths->at(i).first;
      path.at(cell0).at(i).second.resize(paths->at(i).second.size());
      for (int j = 0; j < paths->at(i).second.size(); j++) {
        path.at(cell0).at(i).second.at(j) = paths->at(i).second.at(j);
      }
    }
    //print_path(cell0);

    cg.get_circ_class(cell0, &circ_class);
    
  }
  else {
    std::cout << cell0+1 << " is freed already. " << std::endl;
  }

}

// Shift gmi indices to the graph indices.
void cell_ins_chk::update_path_gmi(int cell0) {
  update_path(cell0-1);
}

void cell_ins_chk::clear_path(int cell0) {

  cells3.at(cell0).conn.clear();
  cells2.at(cell0).conn.clear();
  cells1.at(cell0).conn.clear();
  circ_class.at(2*cell0) = 0;
  circ_class.at(2*cell0+1) = 0;
  circ_tup.at(cell0).clear();
  path.at(cell0).clear();

}

// Return the pointer to the cellbase structure.
cell_base* cell_ins_chk::get_cbase() {
  return cb;
}


// Returns the circuit type, 0 for both containing at least a 3cell, 
// 1 for at least 2cells, 
// other cases for degenerate.
// Circuits are sorted such that the 3-cell containing nondegenerate circuits are
// followed by nondegenerate circuits, which are followed by the degenerate ones.
// circ_class contains the end position of the first and second kind.
int cell_ins_chk::get_circ_type(int cell0, int circ) {
  // Cell0 index is not available.
  if (cells3.size() < cell0 or cell0 < 0)
    return 2;

  // Circuit index is not available or degenerate:
  if (circ > circ_class.at(2*cell0+1))
    return 2;

  // 2cell type:
  else if (circ > circ_class.at(2*cell0))
    return 1;

  // 3cell type:
  else
    return 0;
}

int cell_ins_chk::get_circ_type_gmi(int cell0, int circ) {
  // Cell0 index is not available.
  return get_circ_type(cell0-1,circ);
}

int cell_ins_chk::get_circ_sz() {
  return circ_tup.size();
}

int cell_ins_chk::get_circ_sz(int cell0) {
  return circ_tup.at(cell0).size();
}

int cell_ins_chk::get_circ_sz_gmi(int cell0) {
  return get_circ_sz(cell0-1);
}

// Return the circuit circ in the neighborhood of cell0. The circuit and 
// disjoint graph indices are cell_base indices. cell0 has a cell_base index.
circuit* cell_ins_chk::get_circ(int cell0, int circ) {
  return &circ_tup.at(cell0).at(circ);
}

// Return the circuit circ in the neighborhood of cell0. The circuit and 
// disjoint graph indices are cell_base indices. cell0 has a gmi index.
circuit * cell_ins_chk::get_circ_gmi(int cell0, int circ) {
  return &circ_tup.at(cell0-1).at(circ);
}

// Given the circuit, return the actual cell indices in the cell_base 
// structure. 
void cell_ins_chk::get_circ_topo(ent_conn* c2, 
                                 ent_conn* c2_circ, ent_conn* c3_circ, 
                                 circuit* circ_in, int g_id, int cell_id) {

  c2->conn.clear();
  c2_circ->conn.clear();
  c3_circ->conn.clear();

  assert(cell_id < cells3.size());
  //std::cout << "Size of circ_in: "<< circ_in->first.size() 
  //    << ", size of cell3: " << cells3[cell_id].conn.size() << std::endl;
  // For the circuit, simply convert the local indices to the global indices.
  for(int i=0; i < circ_in->first.first.size(); i++) {

    if (circ_in->first.first.at(i) != -1) {
      //std::cout << circ_in->first[i] << " ";
      int c3_cur = cells3.at(cell_id).conn.at(circ_in->first.first.at(i) );
      //std::cout << c3_cur << std::endl;
      c3_circ->add_ent(c3_cur);
    }
    else
      c3_circ->add_ent(-2);
  }

  for(int i=0; i < circ_in->first.second.size(); i++) {
    int c2_cur = cells2.at(cell_id).conn.at(circ_in->first.second.at(i) - 
                                              cells3.at(cell_id).conn.size());
    c2_circ->add_ent(c2_cur);
  }

  // For the disjoint graphs, use the selected graph given as input.

  std::vector<int>* g_dis;
  if (g_id == 0) {
    g_dis = &circ_in->second.first;
  }
  else 
    g_dis = &circ_in->second.second;

  struct ent_conn e_neigh;

  // Going over the disjoint graph entities, get the 2cells, which are adjacent 
  // to the path 3cells.
  for(int i=0; i < g_dis->size(); i++) {

    if (g_dis->at(i) >= (int)(cells3.at(cell_id).conn.size())) {
      int c2_cur = cells2.at(cell_id)
                      .conn.at(g_dis->at(i) - cells3.at(cell_id).conn.size());
      cb->get_conn_dim(3, 2, c2_cur, &e_neigh);

      for(int j = 0; j < e_neigh.conn.size(); j++) {
        if (c3_circ->chk_ent(e_neigh.conn.at(j))) {
          c2->add_ent(c2_cur);
          j = e_neigh.conn.size();
        }
      }
    }

  }

}

void cell_ins_chk::get_circ_topo_gmi(ent_conn* c2, 
                                 ent_conn* c2_circ, ent_conn* c3_circ, 
                                 circuit* circ_in, int g_id, int cell_id) {

  get_circ_topo(c2, c2_circ, c3_circ, circ_in, g_id, cell_id-1);
  for (int i = 0; i < c2->conn.size(); i++)
    c2->conn.at(i)++;

  for (int i = 0; i < c3_circ->conn.size(); i++)
    c3_circ->conn.at(i)++;

  for (int i = 0; i < c2_circ->conn.size(); i++)
    c2_circ->conn.at(i)++;
}

// Given the circuit, return the actual cell indices in the cell_base 
// structure of the 3cells and the 2cells, which are adjacent to either only 
// 3cells on the path, or to only 3cells in the disjoint graph. 
// The edges on these cells are to be moved to the new vertex by the disc 
// expansion.
// So, these edges can be used to find the disc.
void cell_ins_chk::get_circ_topo_dis(std::vector<ent_conn>* cs, 
                                 ent_conn* c2_circ, ent_conn* c3_circ, 
                                 circuit* circ_in, int g_id, int cell_id) {

  for(int i = 0; i < cs->size(); i++)
    cs->at(i).conn.clear();
  cs->clear();

  c2_circ->conn.clear();
  c3_circ->conn.clear();

  cs->resize(3);

  assert(cell_id < cells3.size());

  //std::cout << "Size of circ_in: "<< circ_in->first.size() 
  //    << ", size of cell3: " << cells3[cell_id].conn.size() << std::endl;
  // For the circuit, simply convert the local indices to the global indices.
  //std::cout << "Circuit: " << std::endl;
  for(int i=0; i < circ_in->first.first.size(); i++) {

    if (circ_in->first.first.at(i) != -1) {
      //std::cout << circ_in->first.first[i] << " ";
      int c3_cur = cells3.at(cell_id).conn.at(circ_in->first.first.at(i) );
      //std::cout << "3c " << c3_cur + 1 << std::endl;
      c3_circ->add_ent(c3_cur);
    }
    else {
      c3_circ->add_ent(-2);
    }
  }

  for(int i=0; i < circ_in->first.second.size(); i++) {
    int c2_cur = cells2.at(cell_id).conn.at(circ_in->first.second.at(i) - 
                                            cells3.at(cell_id).conn.size());
    //std::cout << "2c " << c2_cur + 1 << std::endl;
    c2_circ->add_ent(c2_cur);
  }

  // For the disjoint graphs, use the selected graph given as input.

  std::vector<int>* g_dis;
  if (g_id == 0) {
    g_dis = &circ_in->second.first;
  }
  else 
    g_dis = &circ_in->second.second;

  struct ent_conn e_neigh;

  for(int i=0; i < g_dis->size(); i++) {
    // 2cell
    if (g_dis->at(i) >= (int)(cells3.at(cell_id).conn.size())) {
      int c2_cur = cells2.at(cell_id).conn.at(g_dis->at(i) - cells3.at(cell_id).conn.size());
      //std::cout << "2c" << c2_cur + 1 << " ";
      cs->at(1).add_ent(c2_cur);
    }
    // 3cell
    else {
      if(g_dis->at(i) == -1)
        cs->at(2).add_ent(-2);
      else
        cs->at(2).add_ent(cells3.at(cell_id).conn.at(g_dis->at(i)));
      //std::cout << "3c" << cells3[cell_id].conn.at(g_dis->at(i))  + 1 << " ";
    }
  }

  // Going over the 2cells, check for every bounding 1cell, adjacent to the 
  // 0 cell, add it to cs->at(0). 
  std::map <int, bool> c1_trial{};

  for(int i=0; i < cs->at(1).conn.size(); i++) {
    cb->get_conn_12(cell_id, cs->at(1).conn.at(i), &e_neigh);
    for(int j=0; j < e_neigh.conn.size(); j++) {

      if(!c1_trial[e_neigh.conn.at(j)]) {
        cs->at(0).add_ent(e_neigh.conn.at(j));
        //std::cout << "1c" << e_neigh.conn.at(j) + 1 << " ";
        c1_trial[e_neigh.conn.at(j)] = true;
      }
    }
  }

}

void cell_ins_chk::get_circ_topo_dis_gmi(std::vector<ent_conn>* cs,
                                 ent_conn* c2_circ, ent_conn* c3_circ, 
                                 circuit* circ_in, int g_id, int cell_id) {

  get_circ_topo_dis(cs, c2_circ, c3_circ, circ_in, g_id, cell_id-1);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < cs->at(i).conn.size(); j++)
      cs->at(i).conn.at(j)++;
  }

  for (int i = 0; i < c3_circ->conn.size(); i++)
    c3_circ->conn.at(i)++;

  for (int i = 0; i < c2_circ->conn.size(); i++)
    c2_circ->conn.at(i)++;
}

// Given the circuit, return the actual 2cell indices of a disjoint graph in 
// the cell_base structure. 
void cell_ins_chk::get_circ_2c_disj(ent_conn* c2, circuit* circ_in, 
                                    int g_id, int cell_id) {

  std::vector<int>* g_dis;
  if (g_id == 0) {
    g_dis = &circ_in->second.first;
  }
  else 
    g_dis = &circ_in->second.second;

  struct ent_conn e_neigh;

  // Going over the disjoint graph entities, get the 2cells, which are adjacent 
  // to the path 3cells.
  for(int i=0;i<g_dis->size(); i++) {

    if (g_dis->at(i) >= (int)(cells3.at(cell_id).conn.size())) {
      int c2_cur = cells2.at(cell_id).conn.at(g_dis->at(i) - cells3.at(cell_id).conn.size());
      c2->add_ent(c2_cur);
    }

  }

}

// Using the gmi indices.
void cell_ins_chk::get_circ_2c_disj_gmi(ent_conn* c2, circuit* circ_in, 
                int g_id, int cell_id) {

  get_circ_2c_disj(c2, circ_in, g_id, cell_id-1);
  for (int i = 0; i < c2->conn.size(); i++)
    c2->conn.at(i)++;
}


// Given the 0cell, return the number of disconnected 3cell couples. 
int cell_ins_chk::get_path_sz(int cell0) {
  assert (path.size() > cell0);
  return path.at(cell0).size();
}

// Given the 0cell(gmi), return the all 3cells member of a disconnected 3cell 
// couple. Mainly used in mesh preconditioning step of cell insertion.
std::vector<int> cell_ins_chk::get_path_3cells(int cell0) {
  std::vector<int> e_con(0);
  e_con.reserve(cells3.at(cell0-1).conn.size());

  std::cout << "Getting path 3cells of 0cell" << cell0 << std::endl;
  for (int i = 0; i < get_path_sz(cell0-1); i++) {
    std::pair<int, int> i_gmi = path.at(cell0-1).at(i).first;
    i_gmi.first = i_gmi.first + 1;
    i_gmi.second = i_gmi.second + 1;
    std::cout << "\t3c" << i_gmi.first << "/3c" << i_gmi.second << std::endl;
    if(std::find(e_con.begin(), e_con.end(), i_gmi.first) == e_con.end()) {
      e_con.push_back(i_gmi.first);
      std::cout << "\t3c" << i_gmi.first << std::endl;
    }
    if(std::find(e_con.begin(), e_con.end(), i_gmi.second) == e_con.end()) {
      e_con.push_back(i_gmi.second);
      std::cout << "\t3c" << i_gmi.second << std::endl;
    }
  }
  return e_con;
}

// Given the 0cell and disconnected 3cell couple id, return a pointer to 
// the cell_path. 
cell_path* cell_ins_chk::get_path(int cell0, int path_id) {
  return &path.at(cell0).at(path_id);
}

// Return the 3cells to be connected, in gmi indices. 
// TODO ugly std::pair copying.
std::pair<int, int> cell_ins_chk::get_path_3c_gmi(int cell0, int path_id) {
  std::pair<int, int> i_gmi = path.at(cell0-1).at(path_id).first;
  i_gmi.first = i_gmi.first + 1;
  i_gmi.second = i_gmi.second + 1;
  return i_gmi;
}

void cell_ins_chk::get_path_ngon_gmi(int cell0, int path_id, ngon_gmi* ng) {
  ngon_tracer n_tr(this, ng);
  n_tr.get_path_ngon_gmi(cell0, path_id);

}

// Given the 0cell(gmi), return all 3cells(cbase) bounded by the 0cell.
std::vector<int> cell_ins_chk::get_3cells(int cell0) {
  return cells3.at(cell0-1).conn;
}

// Given the 0cell(gmi), return all 3cells(gmi) bounded by the 0cell.
std::vector<int> cell_ins_chk::get_3cells_gmi(int cell0) {
  std::vector<int> ent(0);
  ent = cells3.at(cell0-1).conn;
  for(int i = 0; i < ent.size(); i++)
    ent.at(i)++;
  return ent;
}

void cell_ins_chk::find_slice(int c0, std::vector< 
          std::pair<std::vector<int >, std::vector<int > > >* path_cells, 
          std::vector< std::pair< std::pair<int,int>, 
                       std::vector<std::vector<int > > > >* slice_cells) {

  cg.get_23adj(c0);
  cg.find_slice(path_cells, slice_cells);
}

// TODO this only inserts a 2cell.
void cell_ins_chk::find_slice(path_adder* pa) {

  int c0 = pa->get_c0();
  int p_id = pa->get_pid();
  int ng_id = pa->get_ng();

  //std::cout << "p_id " << p_id << " ng_id " << ng_id << std::endl;

  ngon_gmi ng;
  get_path_ngon_gmi(c0+1, p_id, &ng);
  //ng.print();

  std::vector<std::pair<std::vector<int >, std::vector<int > > > path_cells
     (0, std::make_pair(std::vector<int >(0), std::vector<int >(0) ) );

  std::vector< std::pair< std::pair<int,int>, 
                          std::vector<std::vector<int > > > > slice_cells
     (0, std::make_pair(
            std::make_pair(0,0), 
            std::vector<std::vector<int > > (0, 
                        std::vector<int >(0) )) );

  slice_cells.resize(0);
  path_cells.resize(0);

  path_cells.resize(ng.ngons.at(ng_id).size());
  for(int i = 0; i < ng.ngons.at(ng_id).size(); i++ ) {
    int path_id = ng.ngons.at(ng_id).at(i);

    //c_base id
    // 3c
    path_cells.at(i).first.reserve(ng.cells.at(path_id).first.size());
    for(int j = 0; j < ng.cells.at(path_id).first.size(); j++) {
      path_cells.at(i).first.push_back(ng.cells.at(path_id).first.at(j) - 1);
    }

    // 2c
    path_cells.at(i).second.reserve(ng.cells.at(path_id).second.size());
    for(int j = 0; j < ng.cells.at(path_id).second.size(); j++) {
      path_cells.at(i).second.push_back(ng.cells.at(path_id).second.at(j) - 1);
    }
  }

  cg.get_23adj(c0);
  cg.find_slice(&path_cells, &slice_cells);

  pa->ca->dim = 2;
  pa->ca->tag = cb->use_free(2);
  pa->ca->tag_0 = c0;

  std::vector<int> new_1cell(0);
  new_1cell.reserve(ng.ngons.at(ng_id).size());
  for(int i = 0; i < ng.ngons.at(ng_id).size(); i++)
    new_1cell.push_back(cb->use_free(1));

  //c_base id
  std::vector<int> new_0cell(0);
  new_0cell.reserve(ng.ngons.at(ng_id).size());
  new_0cell.push_back(c0);

  for(int i = 0; i < ng.ngons.at(ng_id).size()-1; i++)
    new_1cell.push_back(cb->use_free(0));

  std::vector<std::vector<int> > c0s(0, std::vector<int>(0));
  c0s.resize(3);
  for(int i = 0; i < c0s.size(); i++) {
    c0s.at(i).resize(2);
    c0s.at(i).at(0) = -1;
    c0s.at(i).at(1) = -1;
  }

  std::vector<int> e_cell(0);

  for(int i = 0; i < new_0cell.size(); i++) {
    // local id
    std::vector<int> c1_adj(0);
    c1_adj = slice_cells.at(i).second.at(0);

    std::cout << "upper 1c(local) of " << new_0cell.at(i) + 1
              << std::endl;

    for(int j = 0; j < c1_adj.size(); j++) {
      std::cout << "1c" << c1_adj.at(j)+1 << " ";
    }
    std::cout << std::endl;

    //c_base id
    repl_index(1, c0, &c1_adj);

    std::cout << "upper 1c(local) of " << new_0cell.at(i) + 1
              << std::endl;
    for(int j = 0; j < c1_adj.size(); j++) {
      std::cout << "1c" << c1_adj.at(j)+1 << " ";
    }
    std::cout << std::endl;

    pa->ca->add_cell(0, new_0cell.at(i), &e_cell, &c1_adj);
    // TODO Adds the same 0cell twice.
    std::cout << "path_id 1 " << slice_cells.at(i).first.first - 1
              << " path_id 2 " << slice_cells.at(i).first.second - 1
              << std::endl;

    // path_id
    int cc = slice_cells.at(i).first.first - 1;

    if(c0s.at(cc).at(0) == -1)
      c0s.at(cc).at(0) = new_0cell.at(i);
    else
      c0s.at(cc).at(1) = new_0cell.at(i);

    std::cout << "path_id 1 " << cc << " 0c(1) " << c0s.at(cc).at(0) + 1
              << " 0c(2) " << c0s.at(cc).at(1) + 1
              << std::endl;

    cc = slice_cells.at(i).first.second - 1;
    if(c0s.at(cc).at(0) == -1)
      c0s.at(cc).at(0) = new_0cell.at(i);
    else
      c0s.at(cc).at(1) = new_0cell.at(i);

    std::cout << "path_id 2 " << cc << " 0c(1) " << c0s.at(cc).at(0) + 1
              << " 0c(2) " << c0s.at(cc).at(1) + 1
              << std::endl;
  }

  // Path id of the 1cells.
  for(int i = 0; i < new_1cell.size(); i++) {
    std::cout << "1c" << new_1cell.at(i) + 1 << " is being added to ";
    for(int j = 0; j < path_cells.at(i).second.size(); j++) {
      std::cout << "2c" << path_cells.at(i).second.at(j) + 1 << " ";
    }
    std::cout << std::endl;

    // Slice id of the 0cells.
    //c_base id
    pa->ca->add_cell(1, new_1cell.at(i), &c0s.at(i), &path_cells.at(i).second);
  }

  e_cell.clear();
  e_cell.reserve(2);
  if(ng.cp.first > -1)
    e_cell.push_back(ng.cp.first-1);
  if(ng.cp.second > -1)
    e_cell.push_back(ng.cp.second-1);

  pa->ca->add_cell(2, pa->ca->tag, &new_1cell, &e_cell);
}

void cell_ins_chk::find_circuit(path_adder* pa) {

  int c0 = pa->get_c0();
  int circ_id = pa->get_pid();

  std::vector<ent_conn> cs(0, ent_conn());
  cs.resize(3);
  circuit* circ_tup;

  circ_tup = get_circ_gmi(c0+1, circ_id);
  ent_conn c2_circ, c3_circ;

  get_circ_topo_dis_gmi(&cs, &c2_circ, &c3_circ, circ_tup, 0, c0+1); 

  std::vector<std::pair<std::vector<int >, std::vector<int > > > path_cells
      (0, std::make_pair(std::vector<int >(0), std::vector<int >(0)));

  std::vector< std::pair< std::pair<int,int>, 
               std::vector<std::vector<int > > > > slice_cells
     (0, std::make_pair(
            std::make_pair(0,0), 
            std::vector<std::vector<int > > (0, 
                        std::vector<int >(0) )) );

  slice_cells.resize(0);
  path_cells.resize(1);

  path_cells.at(0).first.reserve(c3_circ.conn.size());
  for(int j = 0; j < c3_circ.conn.size(); j++) {
    path_cells.at(0).first.push_back(c3_circ.conn.at(j) - 1);
  }

  path_cells.at(0).second.reserve(c2_circ.conn.size());
  for(int j = 0; j < c2_circ.conn.size(); j++) {
    path_cells.at(0).second.push_back(c2_circ.conn.at(j) - 1);
  }

  cg.get_23adj(c0);
  cg.find_slice(&path_cells, &slice_cells);

  // Get the new cells.
  pa->ca->dim = 1;
  pa->ca->tag = cb->use_free(1);
  pa->ca->tag_0 = c0;

  //c_base id
  std::vector<int> new_0cell(0);
  new_0cell.reserve(2);
  new_0cell.push_back(c0);
  new_0cell.push_back(cb->use_free(0));

  std::vector<int> e_cell(0);

  for(int i = 0; i < new_0cell.size(); i++) {
    std::vector<int> c1_adj(0);
    c1_adj = slice_cells.at(i).second.at(0);

    std::cout << "upper 1c(local) of " << new_0cell.at(i) + 1
              << std::endl;

    for(int j = 0; j < c1_adj.size(); j++) {
      std::cout << "1c" << c1_adj.at(j)+1 << " ";
    }
    std::cout << std::endl;

    //c_base id
    repl_index(1, c0, &c1_adj);

    std::cout << "upper 1c(local) of " << new_0cell.at(i) + 1
              << std::endl;
    for(int j = 0; j < c1_adj.size(); j++) {
      std::cout << "1c" << c1_adj.at(j)+1 << " ";
    }
    std::cout << std::endl;

    pa->ca->add_cell(0, new_0cell.at(i), &e_cell, &c1_adj);
    // TODO Adds the same 0cell twice.
    std::cout << "path_id 1 " << slice_cells.at(i).first.first - 1
              << " path_id 2 " << slice_cells.at(i).first.second - 1
              << std::endl;
  }

  //c_base id
  std::cout << "1c" << pa->ca->tag+1 << " is being added to ";
  for(int j = 0; j < c2_circ.conn.size(); j++) {
    c2_circ.conn.at(j) = c2_circ.conn.at(j) - 1;

    std::cout << "2c" << c2_circ.conn.at(j) << " ";
  }
  std::cout << std::endl;

  pa->ca->add_cell(1, pa->ca->tag, &new_0cell, &c2_circ.conn);
}

void cell_ins_chk::repl_index(int dim, int c0, std::vector<int>* c_in) {
  if(dim == 1) {
    for(int i = 0; i < c_in->size(); i++) {
      c_in->at(i) = cells1.at(c0).conn.at(c_in->at(i));
    }
  }
  else if(dim == 2) {
    for(int i = 0; i < c_in->size(); i++) {
      c_in->at(i) = cells2.at(c0).conn.at(c_in->at(i));
    }
  }
  else if(dim == 3) {
    for(int i = 0; i < c_in->size(); i++) {
      c_in->at(i) = cells3.at(c0).conn.at(c_in->at(i));
    }
  }
}

std::vector<int> cell_ins_chk::ret_ins_gmi() {
  std::vector<int> cell_ins(0);
  cell_ins.reserve(cells3.size());

  std::cout << "Finding insertible 0cells" << std::endl;
  for(int i = 0; i < cells3.size(); i++) {
    std::cout << "0c" << i+1 << std::endl;
    int c3_nbr = 4;
    if(cb->get_0c_corner_gmi(i+1)) {
      int ext_2 = 0;
      for(int j = 0; j < cells2.at(i).conn.size(); j++)
        ext_2 = ext_2 + 1;
      c3_nbr = 4 - ext_2;
    }
    else if(cb->chk_0cell_ext_gmi(i+1))
    //if(cb->chk_0cell_ext_gmi(i+1))
      c3_nbr = 3;
      //c3_nbr = 5;
    if((int)(cells3.at(i).conn.size()) > c3_nbr) {
      std::cout << "More than " << c3_nbr << " 3-cell" << std::endl;
      for(int j = 0; j < cells3.at(i).conn.size(); j ++)
        std::cout << "3c" << cells3.at(i).conn.at(j) + 1 << " ";
      std::cout << std::endl;
      cell_ins.push_back(i+1);
    }
    else {
      bool c1_pass = true;
      for(int j = 0; j < cells1.at(i).conn.size(); j++) {
        if(cb->get_conn_sz(1, cells1.at(i).conn.at(j)) == 1) {
          std::cout << "1c" << cells1.at(i).conn.at(j) + 1 << "bounded only by"
                    << " 0c" << i+1
                    << std::endl;
          c1_pass = false;
          cell_ins.push_back(i+1);
          j = cells1.at(i).conn.size();
        }
      }
      if(c1_pass) {
        std::vector<ent_conn> e_con(0);
        e_con.resize(cells2.at(i).conn.size());
        for(int j = 0; j < cells2.at(i).conn.size(); j++) {
          cb->get_conn_dim(3, 2, cells2.at(i).conn.at(j), &e_con.at(j));
          std::sort(e_con.at(j).conn.begin(), e_con.at(j).conn.end());
          std::cout << "2c" << cells2.at(i).conn.at(j) + 1 << ": ";
          for(int k = 0; k < e_con.at(j).conn.size(); k++) {
            std::cout << "3c" << e_con.at(j).conn.at(k) + 1 << " ";
          }
          std::cout << std::endl;
        }
        for(int j = 0; j < cells2.at(i).conn.size(); j++) {
          for(int k = j+1; k < cells2.at(i).conn.size(); k++) {
            if(e_con.at(j).conn.size() == 1 or e_con.at(k).conn.size() == 1) {
/*
              // Boundary 2cells have the same adjacencies.
              if(e_con.at(j).conn.size() == 1 and 
                 e_con.at(k).conn.size() == 1) {
                if(e_con.at(j).conn.at(0) == e_con.at(k).conn.at(0)) {
                  std::cout << "Upward 3-cell ext" << std::endl;
                  cell_ins.push_back(i+1);
                  j = cells2.at(i).conn.size();
                  k = j;
                }
              }
*/
            }
            // Spurious 1-cell
            else if(e_con.at(j).conn == e_con.at(k).conn) {
              std::cout << "Upward 3-cell the same" << std::endl;
              cell_ins.push_back(i+1);
              j = cells2.at(i).conn.size();
              k = j;
            }
          }
        }
      }
    }
  }

  return cell_ins;
}

bool cell_ins_chk::insertible(int tag_0cell, bool isotropic) {
  if (isotropic) {
    struct ent_conn c_1;
    cb->get_conn_lower(1, 2, &cells2.at(tag_0cell-1), &c_1);
    cb->rem_conn(0, tag_0cell-1, 1, &c_1);
    int n_c3 = cells3.at(tag_0cell-1).conn.size();
    int n_c2 = cells2.at(tag_0cell-1).conn.size();
    int n_c1 = c_1.conn.size();
    int g_1 = n_c1;
    // 4f_2 − 2f_1 = 4f_3 + g_1 − 4,
    // 4f_2 − 2f_1 - 4f_3 + 4 - g_1 = 0,
    int euler_char = 4*n_c2 - 2*n_c1 - 4*n_c3 + 4 - g_1;
    std::cout << "c3 " << n_c3 << " c2 " << n_c2 
              << " c1 " << n_c1
              << std::endl;
    std::cout << "Euler characteristic " << euler_char << std::endl;
    if ( euler_char == 0) 
      return false;
    else
      return true;
  }
  return true;
}

void cell_ins_chk::set_calc_corner(bool fix) {
  calc_corner = fix;
  cg.set_calc_corner(fix);
}

cell_ins_chk::cell_ins_chk() :
    cg(), circ_class(0),
    circ_tup(0, std::vector< circuit >(0)),
    path(0, std::vector< cell_path >(0)),
    load_flag(false), calc_corner(false),
    cells3(0), cells2(0), cells1(0) {
}

cell_ins_chk::cell_ins_chk(struct cell_base* c_base, bool corner_flag) :
    cg(), circ_class(0),
    circ_tup(0, std::vector< circuit >(0)),
    path(0, std::vector< cell_path >(0)),
    load_flag(false), calc_corner(false),
    cells3(0), cells2(0), cells1(0) {
  //std::cout << "Loading the cell insert checker..." << std::endl;

  cb = c_base;

  calc_corner = corner_flag;
  //cg = new cell_graph(cb);
  cg.load_cb(cb, corner_flag);
  update();
}

void cell_ins_chk::update() {

  std::cout << "Updating the cell graph..." << std::endl;

  clear();

  int n_0 = cb->get_sz(0);
  circ_tup.resize(n_0);
  path.resize(n_0);

  cells3.resize(n_0);
  cells2.resize(n_0);
  cells1.resize(n_0);

  circ_class.resize(n_0*2);

  std::cout << "Containers resized..." << std::endl;

  for (int i = 0; i< n_0; i++) {
    update_path(i);
  }
  std::cout << "Paths loaded into the cell insert structure..." << std::endl;
  //print_circ();
  //print_path();
}

void cell_ins_chk::load_cb(struct cell_base* c_base, bool corner_flag) {
  clear();
  cb = c_base;

  cg.load_cb(cb, corner_flag);
  update();
  load_flag = true;
}

// ---------------------------------------------------------
// Output related.
// ---------------------------------------------------------
// Print the circuits and paths of the 0cells.

void cell_ins_chk::print_circ() {

  for (int i = 0; i < circ_tup.size(); i++) {
    std::cout << "\n" << "Printing circuits of 0cell " << i << std::endl;

    for (int ii = 0; ii < circ_tup[i].size(); ii++) {
      std::cout << "\n" << "Circuit " << ii << std::endl;
      for (int j = 0; j < circ_tup[i][ii].first.first.size(); j++) {
        int U = circ_tup.at(i)[ii].first.first.at(j);
        if (U == -1)
          std::cout << "Ext" << ", ";
        else
          std::cout << "3c" << cells3.at(i).conn.at(U) + 1 << ", ";

      }
      for (int j = 0; j < circ_tup.at(i)[ii].first.second.size(); j++) {
        int U = circ_tup.at(i)[ii].first.second.at(j);
        std::cout << "2c" 
              << cells2.at(i).conn.at(U-cells3.at(i).conn.size()) + 1 << ", ";
      }

      std::cout << std::endl;

      std::cout << "Subgraph 1_tup: " << std::endl;
      for (int j = 0; j < circ_tup.at(i)[ii].second.first.size(); j++) {
        int cell_nbr = circ_tup.at(i)[ii].second.first.at(j);
        if (cell_nbr == -1)
          std::cout << "Ext" << ", ";
        else if (cell_nbr < cells3.at(i).conn.size())
          std::cout << "3c" << cells3.at(i).conn.at(cell_nbr) + 1 << ", ";
        else
          std::cout << "2c" 
                << cells2.at(i).conn.at(cell_nbr-cells3.at(i).conn.size()) + 1 << ", ";
      }
      std::cout << std::endl;

      std::cout << "Subgraph 2_tup: " << std::endl;
      for (int j = 0; j < circ_tup[i][ii].second.second.size(); j++) {
        int cell_nbr = circ_tup[i][ii].second.second[j];
        if (cell_nbr == -1)
          std::cout << "Ext" << ", ";
        else if (cell_nbr < cells3.at(i).conn.size())
          std::cout << "3c" << cells3.at(i).conn.at(cell_nbr) + 1 << ", ";
        else
          std::cout << "2c" 
                << cells2.at(i).conn.at(cell_nbr-cells3.at(i).conn.size()) + 1 << ", ";
      }
      std::cout << std::endl;
    }

  }
}


void cell_ins_chk::print_path() {
  std::cout << "Printing paths: " << std::endl;
  for (int i = 0; i < path.size(); i++) {
    print_path(i);
  }
}


void cell_ins_chk::print_circ(int cell) {
  assert(cell < circ_tup.size());

  for (int ii = 0; ii < circ_tup[cell].size(); ii++) {
    std::cout << "\n" << "Circuit " << ii << std::endl;
    for (int j = 0; j < circ_tup[cell][ii].first.first.size(); j++) {
      int U = circ_tup[cell][ii].first.first.at(j);
      if (U == -1)
        std::cout << "Ext" << ", ";
      else 
        std::cout << "3c" << cells3.at(cell).conn.at(U) + 1 << ", ";
    }
    for (int j = 0; j < circ_tup[cell][ii].first.second.size(); j++) {
      int U = circ_tup[cell][ii].first.second.at(j);
      std::cout << "2c" 
            << cells2.at(cell).conn.at(U-cells3.at(cell).conn.size()) + 1 << ", ";
    }
    std::cout << std::endl;

    std::cout << "Subgraph 1_tup: " << std::endl;
    for (int j = 0; j < circ_tup[cell][ii].second.first.size(); j++) {
      int cell_nbr = circ_tup[cell][ii].second.first[j];
      if (cell_nbr == -1)
        std::cout << "Ext" << ", ";
      else if (cell_nbr < cells3.at(cell).conn.size())
        std::cout << "3c" << cells3.at(cell).conn.at(cell_nbr) + 1 << ", ";
      else
        std::cout << "2c" 
          << cells2.at(cell).conn.at(cell_nbr-cells3.at(cell).conn.size()) + 1 << ", ";
    }
    std::cout << std::endl;

    std::cout << "Subgraph 2_tup: " << std::endl;
    for (int j = 0; j < circ_tup[cell][ii].second.second.size(); j++) {
      int cell_nbr = circ_tup[cell][ii].second.second[j];
      if (cell_nbr == -1)
        std::cout << "Ext" << ", ";
      else if (cell_nbr < cells3.at(cell).conn.size())
        std::cout << "3c" << cells3.at(cell).conn.at(cell_nbr) + 1 << ", ";
      else
        std::cout << "2c" 
         << cells2.at(cell).conn.at(cell_nbr-cells3.at(cell).conn.size()) + 1 << ", ";
    }
    std::cout << std::endl;
  }
}

void cell_ins_chk::print_path(int i) {
  std::cout << "Printing paths of 0Cell " << i+1 << " sz " << path.at(i).size() 
            << std::endl;
  for (int ii = 0; ii < path.at(i).size(); ii++) {
    std::cout << path[i][ii].first.first + 1 << "-" 
              << path[i][ii].first.second + 1 << std::endl;
    std::cout << "Paths" << std::endl;
    for (int j = 0; j < path[i][ii].second.size(); j++) {
      std::cout << j << ": " << std::endl;
      for (int k = 0; k < path[i][ii].second[j].first.size(); k++) {
        int cell_nbr = path[i][ii].second[j].first[k];

        if (cell_nbr == -1)
          std::cout << "Ext" << ", ";
        else if (cell_nbr < cells3.at(i).conn.size())
          std::cout << "3c" << cells3.at(i).conn.at(cell_nbr) + 1 << ", ";
        else
          std::cout << "2c" 
                << cells2.at(i).conn.at(cell_nbr-cells3.at(i).conn.size()) + 1 << ", ";
      }
      std::cout << "Non-intersecting" << std::endl;
      for (int k = 0; k < path[i][ii].second[j].second.size(); k++) {
        std::cout << path[i][ii].second[j].second[k] << ", ";
      }
      std::cout << std::endl;
    }
  }
}

// Print all the 0cell topologies into dot files.
void cell_ins_chk::print_graph() {

  std::cout << "Printing graphs into dot files..." << std::endl;

  for (int i = 0; i < circ_tup.size(); i++) {

    if(!cb->is_free(0, i)) {
      std::cout << i << "th 0cell out of "<< circ_tup.size() << std::endl;
      cg.get_23adj(i);

      std::stringstream ss;
      ss << "./output/23adj" << "_" << i+1 << ".dot";
      std::string tmp = ss.str();
      const char* cstr = tmp.c_str();

      cg.print_graph_lbl(cstr);
    }
  }
}

void cell_ins_chk::clear() {

  cg.clear();
  circ_class.clear();

  for(int i = 0; i < circ_tup.size(); i++) {
    for(int j = 0; j < circ_tup.at(i).size(); j++) {
      circ_tup.at(i).at(j).first.first.clear();
      circ_tup.at(i).at(j).first.second.clear();
      circ_tup.at(i).at(j).second.first.clear();
      circ_tup.at(i).at(j).second.second.clear();
    }
    circ_tup.at(i).clear();
  }
  circ_tup.clear();

  for(int i = 0; i < path.size(); i++) {
    path.at(i).clear();
  }
  path.clear();

  cells3.clear();
  cells2.clear();
  cells1.clear();
}

cell_ins_chk::~cell_ins_chk() {
  std::cout << "Destroying cell insert checker...";
  std::cout << "cell_ins_chk: " << this;
  std::cout << " Path: " << &path << " Cb: " << cb << " cg: " << &cg << 
              " circ tup: " << &circ_tup << " circ_clas: " << &circ_class <<
              " cells2: " << &cells2 << " cells3: " << &cells3 <<
              std::endl;
  //delete cg;
  clear();
  std::cout << "Destroyed" << std::endl;

}

// ---------------------------------------------------------
// Ngon tracer. Used in extracting ngons from paths of a given 0cell.
// ---------------------------------------------------------

ngon_tracer::ngon_tracer(cell_ins_chk* c_ins_in, ngon_gmi* ng_in) :
    ng(NULL), c_ins(NULL), c_path(NULL) {
  c_ins = c_ins_in;
  ng = ng_in;
}

void ngon_tracer::get_path_ngon_gmi(int cell0, int path_id) {

  std::cout << "Obtaining the ngons for 0cell" << cell0 << ", " << "couple "
            << path_id << std::endl;

  ng->clear();

  c_path = c_ins->get_path(cell0-1, path_id);
  ng->cells.resize(c_path->second.size());

  if(c_ins->get_cbase()->chk_0cell_ext_gmi(cell0)) {
    ng->bd = true;
  }
  else {
    ng->bd = false;
  }

  // There are at most 2^(n-1) ngons, where n is the number of paths. With 
  // intersecting paths, this number is less. Usually, the number of paths is
  // low.

  int g = 1 << (c_path->second.size() - 1);
  std::cout << "Maximum number of ngon combinations: " << g << std::endl;
  ng->ngons.reserve(g);

  for(int i = 0; i < ng->cells.size(); i++) {
    // 3cells
    ng->cells.at(i).first.reserve(c_path->second.at(i).first.size());
    // 2cells
    ng->cells.at(i).second.reserve(c_path->second.at(i).first.size());

    std::cout << "Checking cell_nbrs: " << std::endl;
    for(int j = 0; j < c_path->second.at(i).first.size(); j++) {
      int cell_nbr = c_path->second.at(i).first.at(j);

      std::cout << cell_nbr << " ";
      if (cell_nbr == -1) {
        std::cout << "(Ext)" << " ";
        ng->cells.at(i).first.push_back(-1);
      }
      else if (cell_nbr < (int)(c_ins->cells3.at(cell0-1).conn.size()))
        ng->cells.at(i).first.push_back(c_ins->cells3.at(cell0-1).conn.at(cell_nbr) + 1);
      else
        ng->cells.at(i).second.push_back(c_ins->cells2.at(cell0-1).conn.at(cell_nbr-c_ins->cells3.at(cell0-1).conn.size()) + 1);

    }

    std::vector<int> ngon_curr(0);
    ngon_curr.push_back(i);
    std::vector<int> in_curr(0);
    in_curr = c_path->second.at(i).second;

    for(int j = 0; j < i; j++) {
      std::vector<int>::iterator newEnd = std::remove(in_curr.begin(),
                                                      in_curr.end(), j);
      if(newEnd != in_curr.end())
        in_curr.erase(newEnd, in_curr.end());
    }
    trace_path(in_curr, ngon_curr);
  }

  ng->cp = c_ins->get_path_3c_gmi(cell0, path_id);

}


void ngon_tracer::trace_path(std::vector<int> in_curr, 
                              std::vector<int> ngon_curr) {
  if (!in_curr.empty()) {
    for(int i = 0; i < in_curr.size(); i++) {
      // Non-intersecting list of the next path.
      std::vector<int> in_path(0);
      in_path = c_path->second[in_curr[i]].second;

      for(int j = 0; j < in_curr.at(i); j++) {
        std::vector<int>::iterator newEnd = std::remove(in_path.begin(),
                                                        in_path.end(), j);
        if(newEnd != in_path.end())
          in_path.erase(newEnd, in_path.end());
      }

      std::vector<int> in_next(0);
      std::set_intersection( in_path.begin(),in_path.end(),
                             in_curr.begin(),in_curr.end(),
                          std::inserter(in_next,in_next.begin()));

      std::vector<int> ngon_next(0);
      ngon_next = ngon_curr;
      ngon_next.push_back(in_curr.at(i));
      ng->ngons.push_back(ngon_next);

      trace_path(in_next, ngon_next);
    }
  }
}

// ngon_tracer for cell_graph

ngon_tracer_g::ngon_tracer_g (cell_graph* cg_in, ngon_gmi* ng_in) : 
  ng(NULL), cg(NULL), c_path(NULL) {
  cg = cg_in;
  ng = ng_in;
}

void ngon_tracer_g::get_path_ngon(int path_id) {

  std::cout << "Obtaining the ngons for 0cell" << cg->get_0c() + 1
            << ", " << "couple "
            << path_id << std::endl;

  ng->clear();

  c_path = &cg->get_paths()->at(path_id);
  ng->cells.resize(c_path->second.size());

  if(cg->get_ext()) {
    ng->bd = true;
  }
  else {
    ng->bd = false;
  }

  // There are at most 2^(n-1) ngons, where n is the number of paths. With 
  // intersecting paths, this number is less. Usually, the number of paths is
  // low.

  int g = 1 << (c_path->second.size() - 1);
  std::cout << "Maximum number of ngon combinations: " << g << std::endl;
  ng->ngons.reserve(g);

  for(int i = 0; i < ng->cells.size(); i++) {
    // 3cells
    ng->cells.at(i).first.reserve(c_path->second.at(i).first.size());
    // 2cells
    ng->cells.at(i).second.reserve(c_path->second.at(i).first.size());

    std::cout << "Checking cell_nbrs: " << std::endl;
    for(int j = 0; j < c_path->second.at(i).first.size(); j++) {
      int cell_nbr = c_path->second.at(i).first.at(j);

      std::cout << cell_nbr << " ";
      if (cell_nbr == -1) {
        std::cout << "Ext" << ", ";
        ng->cells.at(i).first.push_back(-1);
      }
      else if (cell_nbr < (int)(cg->cells3.conn.size()))
        ng->cells.at(i).first.push_back(cg->cells3.conn.at(cell_nbr));
      else
        ng->cells.at(i).second.push_back(cg->cells2.conn.at(cell_nbr-
                                              cg->cells3.conn.size()));
    }

    std::vector<int> ngon_curr(0);
    ngon_curr.push_back(i);
    std::vector<int> in_curr(0);
    in_curr = c_path->second.at(i).second;

    for(int j = 0; j < i; j++) {
      std::vector<int>::iterator newEnd = std::remove(in_curr.begin(),
                                                      in_curr.end(), j);
      if(newEnd != in_curr.end())
        in_curr.erase(newEnd, in_curr.end());
    }
    trace_path(in_curr, ngon_curr);
  }

  std::pair<int, int> i_gmi = c_path->first;
  i_gmi.first = i_gmi.first;
  i_gmi.second = i_gmi.second;

  ng->cp = i_gmi;

}

void ngon_tracer_g::get_path_ngon_gmi(int path_id) {
  get_path_ngon(path_id);
  // There are at most 2^(n-1) ngons, where n is the number of paths. With 
  // intersecting paths, this number is less. Usually, the number of paths is
  // low.

  for(int i = 0; i < ng->cells.size(); i++) {
    for(int j = 0; j < ng->cells.at(i).first.size(); j++) {
      ng->cells.at(i).first.at(j) = ng->cells.at(i).first.at(j) + 1;
    }
    for(int j = 0; j < ng->cells.at(i).second.size(); j++) {
      ng->cells.at(i).second.at(j) = ng->cells.at(i).second.at(j) + 1;
    }
  }
}

void ngon_tracer_g::trace_path(std::vector<int> in_curr, 
                              std::vector<int> ngon_curr) {
  if (!in_curr.empty()) {
    for(int i = 0; i < in_curr.size(); i++) {
      // Non-intersecting list of the next path.
      std::vector<int> in_path(0);
      in_path = c_path->second[in_curr[i]].second;

      for(int j = 0; j < in_curr.at(i); j++) {
        std::vector<int>::iterator newEnd = std::remove(in_path.begin(),
                                                        in_path.end(), j);
        if(newEnd != in_path.end())
          in_path.erase(newEnd, in_path.end());
      }

      std::vector<int> in_next(0);
      std::set_intersection( in_path.begin(),in_path.end(),
                             in_curr.begin(),in_curr.end(),
                          std::inserter(in_next,in_next.begin()));

      std::vector<int> ngon_next(0);
      ngon_next = ngon_curr;
      ngon_next.push_back(in_curr.at(i));
      ng->ngons.push_back(ngon_next);

      trace_path(in_next, ngon_next);
    }
  }
}


path_adder::path_adder(struct cell_base* c_base, 
                                cell_adder* c_add, int cell0) :
    c0(0), p_id(0), ng_id(0), ca(NULL), cb(NULL) {
  ca = c_add;
  cb = c_base;
  c0 = cell0;

  ca->clear();

  p_id = -1;
  ng_id = -1;
}
void path_adder::reload(struct cell_base* c_base, 
                                cell_adder* c_add, int cell0) {
  ca = c_add;
  cb = c_base;
  c0 = cell0;

  ca->clear();

  p_id = -1;
  ng_id = -1;
}

void path_adder::set_path(int path_id, int ngon_id) {
  ca->dim = 2;
  p_id = path_id;
  ng_id = ngon_id;
}

void path_adder::set_circ(int circ_id) {
  ca->dim = 1;
  p_id = circ_id;
}

int path_adder::get_c0() {
  return c0;
}

int path_adder::get_pid() {
  return p_id;
}

int path_adder::get_ng() {
  return ng_id;
}

void path_adder::clear() {

  ca->clear();
}



