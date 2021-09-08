
#include <assert.h>

#include <fstream>
#include <iostream>
#include <iomanip>

#include <vector>

#include <boost/config.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
//#include <boost/graph/johnson_all_pairs_shortest.hpp>

#include "topo_topo.h"
#include "topo_graph.h"

// ---------------------------------------------------------
// Paton elementary circuit finder and Gibbs algorithm.
// ---------------------------------------------------------

bool PatonFinder::edge_comp (Edge e1, Edge e2) { 
  Vertex u1=source(e1.first,g);
  Vertex u2=source(e2.first,g);
  Vertex v1=target(e1.first,g);
  Vertex v2=target(e2.first,g);

  bool u1u2 = (index[u1] == index[u2]);
  bool u1v2 = (index[u1] == index[v2]);
  bool v1u2 = (index[v1] == index[u2]);
  bool v1v2 = (index[v1] == index[v2]);

  if (u1u2) {
    if (v1v2) {
      return false;
    }
    else {
      return (index[v1]<index[v2]);
    }
  }
  else if (u1v2) {
    if (v1u2) {
      return false;
    }
    else {
      return (index[u1]<index[u2]);
    }
  }
  else {
    return (index[u1]<index[u2]);
  }
}

// Assume sorted. If all edges are the same, they will be sorted the same.
// The comparison should hold.
bool PatonFinder::set_edge_comp (std::vector<Edge> e1, std::vector<Edge> e2) { 

  // This is used to sort edge sets such that the shorter sets will be 
  // considered "less".
  if (e1.size() < e2.size()) {
    return true;
  }
  else if (e1.size() > e2.size()) {
    return false;
  }

  std::vector<Edge>::iterator it_e1;
  std::vector<Edge>::iterator it_e2;

  it_e1 = e1.begin(); 
  it_e2 = e2.begin(); 
  while (it_e1 != e1.end() and it_e2 != e2.end()) {
    Vertex u1=source(it_e1->first,g);
    Vertex u2=source(it_e2->first,g);
    Vertex v1=target(it_e1->first,g);
    Vertex v2=target(it_e2->first,g);

    bool u1u2 = (index[u1] == index[u2]);
    bool u1v2 = (index[u1] == index[v2]);
    bool v1u2 = (index[v1] == index[u2]);
    bool v1v2 = (index[v1] == index[v2]);

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

void PatonFinder::sort_edge (std::vector<Edge>* edge_list) { 
  std::sort(edge_list->begin(), edge_list->end(), 
      std::bind(&PatonFinder::edge_comp, this, 
      std::placeholders::_1, std::placeholders::_2) );
}

void PatonFinder::sort_edge_sets (std::vector<std::vector<Edge> >* e_set) { 
  std::sort(e_set->begin(), e_set->end(), 
      std::bind(&PatonFinder::set_edge_comp, this, 
      std::placeholders::_1, std::placeholders::_2) );
}

void PatonFinder::isect_edge (std::vector<Edge>* e1, std::vector<Edge>* e2, 
                 std::vector<Edge>* out) {
  set_intersection(e1->begin(),e1->end(), e2->begin(),e2->end(), 
      std::inserter(*out,out->begin()), 
            std::bind(&PatonFinder::edge_comp, this, 
            std::placeholders::_1, std::placeholders::_2));
}

void PatonFinder::diff_edge (std::vector<Edge>* e1, std::vector<Edge>* e2, 
                 std::vector<Edge>* out) {
  set_difference(e1->begin(),e1->end(),
                    e2->begin(),e2->end(), 
                    std::inserter(*out,out->begin()), 
            std::bind(&PatonFinder::edge_comp, this, 
            std::placeholders::_1, std::placeholders::_2));
}
void PatonFinder::merge_edge (std::vector<Edge>* e1, std::vector<Edge>* e2, 
                 std::vector<Edge>* out) {
  std::merge(e1->begin(),e1->end(),
                    e2->begin(),e2->end(), 
                std::inserter(*out,out->end()), 
            std::bind(&PatonFinder::edge_comp, this, 
            std::placeholders::_1, std::placeholders::_2));
}

void PatonFinder::merge_edge_set (std::vector<std::vector<Edge> >* e1, 
                 std::vector<std::vector<Edge> >* e2, 
                 std::vector<std::vector<Edge> >* out) {
  std::merge(e1->begin(),e1->end(),
                e2->begin(),e2->end(), 
                    std::inserter(*out,out->end()), 
            std::bind(&PatonFinder::set_edge_comp, this, 
            std::placeholders::_1, std::placeholders::_2));
}

PatonFinder::PatonFinder(graph* g_in) {
  g_init = g_in;
  g = *g_in;
  index = get(boost::vertex_index, g);

  for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
    std::cout << *vp.first <<  " ";
    V.insert(*vp.first);
  }
  --vp.first;
  T.insert(*vp.first);
  std::cout << std::endl;
}

    // Elementary circuit detection, recursive one does not work as the 
    // iterator cannot be used in multiple recursions.
bool PatonFinder::find_nonrecursive(Vertex v1, Vertex v2) {

  std::vector<Vertex> path;
  std::map<Vertex, Vertex> predecessor;
  std::map<Vertex, bool> visited;
  //path.push_back(v1);

  std::vector<Vertex> cand;
  predecessor.clear();
  cand.clear();
  cand.push_back(v1);

  visited[v1] = true;

  // Visit every accessible vertex, tag them as visited and record their
  // predecessor. At any given level, each vertex will be connected to 
  // its closest neighbor.
  while (!cand.empty()) {
    for (std::tie(ai_int, ai_end_int) = boost::adjacent_vertices(cand.front(), tree); ai_int != ai_end_int; ++ai_int) {
      if (!visited[*ai_int]) {
        cand.push_back(*ai_int);
        visited[*ai_int] = true;
        predecessor[*ai_int] = cand.front();
      }
    }
    cand.erase(cand.begin());
  }

  if(visited[v2]){
    //std::cout << "Visited " << v2 << std::endl;
  }

  // Collect the vertices starting from v2.
  path.push_back(v2);
  while (path.back() != v1) {
    path.push_back(predecessor[path.back()]);
    //std::cout << " " << path.back() << std::endl;
  }

  // Vector to set conversion. Could be faster when combining different 
  // circuits.
  std::set<Vertex> path_set(std::make_move_iterator(path.begin()),
                std::make_move_iterator(path.end()));

  std::set<Vertex> intersect;
  set_intersection(V.begin(),V.end(),path_set.begin(),path_set.end(),
              std::inserter(intersect,intersect.begin()));

  //std::cout << std::endl;
  while(!intersect.empty()) {
    //std::cout << " " << *intersect.begin() << std::endl;
    intersect.erase(*intersect.begin());
  }
  //std::cout << std::endl;
  paths.push_back(path);

  return true;
}


// Create a subgraph for each path. 
// For each path: 
// For each edge on the path remove the edges connecting the source and target
// vertices. 
// Using the burn algorithm, find disjoint subgraphs. 
// Reload the graph.
// Print out the disjoint subgraphs and paths.
void PatonFinder::create_sub() {
  g = *g_init;

  // Stores the disjoints graphs as couples.
  s_graphs.reserve(2*S.size());
  s_graphs.resize(2*S.size());
  for (int i = 0; i < s_graphs.size(); i++) {
    s_graphs[i].reserve(num_vertices(g));
  }

  // Remove all edges connecting the vertices on the path and the rest of the
  // original graph. Burn the subgraphs.
  for (int i = 0; i < S.size(); i++) {

    // To check if the vertex is burnable. No vertex on path or already burned
    // is burnable.
    std::map<Vertex, bool> burned;

    for (int j = 0; j < S[i].size(); j++) {
      Vertex u=source(S[i][j].first,g);
      Vertex v=target(S[i][j].first,g);
      if (!burned[u]) {
        for (std::tie(ai, ai_end) = boost::adjacent_vertices(u, g); ai != ai_end; ++ai) {
          remove_edge(u, *ai, g);
        }
        burned[u] = true;
      }
      if (!burned[v]) {
        for (std::tie(ai, ai_end) = boost::adjacent_vertices(v, g); ai != ai_end; ++ai) {
          remove_edge(v, *ai, g);
        }
        burned[v] = true;
      }
    }

    // Burn algorithm.
    // Find an unburned vertex on the graph. 
    std::vector<Vertex> burning;
    burning.clear();
    burning.reserve(num_vertices(g));

    //std::cout << "Burning the first subgraph..." << std::endl;
    // Find a vertex on the subgraph.
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
      if (!burned[*vp.first]) {
        //s_graphs[2*i].push_back(*vp.first);
        burning.push_back(*vp.first);
        burned[*vp.first] = true;
        vp.first = vp.second-1;
      }
    }

    // Burn the first subgraph.
    while(!burning.empty()) {
      Vertex v = *burning.begin();
      for (std::tie(ai, ai_end) = boost::adjacent_vertices(v, g); ai != ai_end; ++ai) {
        if (!burned[*ai]) {
          //s_graphs[2*i].push_back(*ai);
          burning.push_back(*ai);
          burned[*ai] = true;
        }
      }
      //std::cout << "Removing " << *burning.begin() << std::endl;
      s_graphs[2*i].push_back(burning[0]);
      burning.erase(burning.begin());
    }

    burning.clear();
    // Find a vertex on the subgraph.
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
      if (!burned[*vp.first]) {
        burning.push_back(*vp.first);
        burned[*vp.first] = true;
        vp.first = vp.second-1;
      }
    }

    //std::cout << "Burning the second subgraph..." << std::endl;
    // Burn the second subgraph.
    while(!burning.empty()) {
      Vertex v = *burning.begin();
      for (std::tie(ai, ai_end) = boost::adjacent_vertices(v, g); ai != ai_end; ++ai) {
        if (!burned[*ai]) {
          burning.push_back(*ai);
          burned[*ai] = true;
        }
      }
      //std::cout << "Removing " << *burning.begin() << std::endl;
      s_graphs[2*i+1].push_back(burning[0]);
      burning.erase(burning.begin());
    }

    // Reload the graph for the next path.
    g = *g_init;
  }

  // Print the paths and associated subgraphs.
  for (int i = 0; i < S.size(); i++) {
    std::cout << i << "th path: " << std::endl;
    for (int j = 0; j < S[i].size(); j++) {
      std::cout << S[i][j].first << " ["
                << S[i][j].second << "] ";
    }
    std::cout << std::endl;

    std::cout << "Subgraph 1: " << std::endl;
    for (int j = 0; j < s_graphs[2*i].size(); j++) {
      std::cout << s_graphs[2*i][j] << " ";
    }
    std::cout << std::endl;

    std::cout << "Subgraph 2: " << std::endl;
    for (int j = 0; j < s_graphs[2*i+1].size(); j++) {
      std::cout << s_graphs[2*i+1][j] << " ";
    }
    std::cout << std::endl;
  }

}

// Convert the elementary circuit lists into edge circuits lists.
// These are used in finding circuits by gibbs algorithm.
void PatonFinder::vert_2_edge() {
  g = *g_init;
  paths_e.clear();
  paths_e.reserve(paths.size());
  paths_e.resize(paths.size());
  for (int i = 0; i < paths_e.size(); i++) {
    paths_e[i].reserve(paths[i].size());
    paths_e[i].resize(paths[i].size());
  }

  for (int i = 0; i < paths_e.size(); i++) {
    for (int j = 0; j < paths[i].size()-1; j++) {
      paths_e[i][j] = boost::edge(paths[i][j], paths[i][j+1], g);
    }
    paths_e[i].back() = boost::edge(paths[i][0], paths[i].back(), g);
  }

  for (int i = 0; i < paths_e.size(); i++) {
    sort_edge (&paths_e[i]);
  }

/*
  paths_eset.reserve(paths.size());
  paths_eset.resize(paths.size());

  for (int i = 0; i < paths_e.size(); i++) {
    std::copy(paths_e[i].begin(), paths_e[i].end(), 
             std::inserter(paths_eset[i], paths_eset[i].end() ));

  }
*/
}

// These are used in finding circuits by gibbs algorithm.
void PatonFinder::edge_2_vert() {
}

    // Using the gibbs algorithm, find all unique circuits.
void PatonFinder::find_gibbs() {
  std::vector< std::vector<Edge> > S_temp;
  std::vector< std::vector<Edge> > Q;
  std::vector< std::vector<Edge> > Q_temp;
  std::vector< std::vector<Edge> > Q_temp2;
  std::vector< std::vector<Edge> > R;
  std::vector< std::vector<Edge> > R_star;
 
  std::vector<Edge> diff;
  std::vector<Edge> diff_join;
  std::vector<Edge> intersect;

  S.push_back(paths_e[0]);
  Q.push_back(paths_e[0]);

  std::vector<Edge>::iterator it_e;

  for (int i = 1; i < paths_e.size(); i++) {
    S_temp.clear();
    Q_temp.clear();
    Q_temp2.clear();
    R.clear();
    R_star.clear();

    // Compare with T in Q.
    for (int j = 0; j < Q.size(); j++) {
      diff.clear();
      diff_join.clear();
      intersect.clear();

      // path-T
      diff_edge(&paths_e[i],&Q[j],&diff);
      sort_edge (&diff);

      // T-path
      diff_edge(&Q[j],&paths_e[i],&intersect);
      sort_edge (&intersect);

      // (T U path) - (path n T)
      merge_edge(&diff,&intersect,&diff_join);
      sort_edge (&diff_join);

      //std::cout << std::endl;
      intersect.clear();
      // T n path
      isect_edge(&paths_e[i], &Q[j], &intersect);

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
      std::cout << "R-Cycle " << iii << ", " << "size: " << R[iii].size() << std::endl;
      for (it_e = R[iii].begin(); 
           it_e != R[iii].end(); ++it_e) {
        std::cout << it_e->first << " ";
      }
      std::cout << std::endl;
    }
*/
    // Track if v is moved to R_star.
    std::vector<bool> in_star;


    in_star.resize(R.size());
    std::fill(in_star.begin(), in_star.end(), false); 

    for (int u = 0; u < R.size(); u++) {
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

          if (intersect.size() == R[u].size() ) {
            //std::cout << "v " << v << std::endl;
            in_star[v] = true;
          }
          else if (intersect.size() == R[v].size() ) {
            //std::cout << "u " << u << std::endl;
            in_star[u] = true;
            v = R.size();
          }
        }
      }// for_3 ends
    }// for_2 ends

    for (int u = R.size()-1; u > -1; u--) {
      if (in_star[u]) {
        //std::cout << "Popping back.." << u << std::endl;
        R_star.push_back(R[u]);
        std::swap(R[u], R.back());
        R.pop_back();
      }
    }// for_2 ends

    //std::cout << "S-size " << S.size() << ", R-size " << R.size() 
    //<<std::endl;

    // Sort sets
    for (int iii = 0; iii < S.size(); iii++)
      sort_edge(&S.at(iii));
    for (int iii = 0; iii < R.size(); iii++)
      sort_edge(&R.at(iii));
    sort_edge_sets (&S);
    sort_edge_sets (&R);

    //std::cout << "S-size " << S.size() << ", R-size " << R.size() 
    //<<std::endl;

/*
    for (int iii = 0; iii < S.size(); iii++) {
      //std::cout << "S-Cycle " << iii << ", " << "size: " 
      //<< S[iii].size() << std::endl;
      for (it_e = S[iii].begin(); 
           it_e != S[iii].end(); ++it_e) {
        std::cout << it_e->first << " ";
      }
      std::cout << std::endl;
    }

    for (int iii = 0; iii < R.size(); iii++) {
      std::cout << "R-Cycle " << iii << ", " << "size: " << R[iii].size() << std::endl;
      for (it_e = R[iii].begin(); it_e != R[iii].end(); ++it_e) {
        std::cout << it_e->first << " ";
      }
      std::cout << std::endl;
    }
*/

    merge_edge_set(&S,&R,&S_temp);

/*
    for (int iii = 0; iii < S_temp.size(); iii++) {
      std::cout << "Stemp-Cycle " << iii << ", " << "size: " 
      << S_temp[iii].size() << std::endl;
      for (it_e = S_temp[iii].begin(); it_e != S_temp[iii].end(); ++it_e) {
        std::cout << it_e->first << " ";
      }
      std::cout << std::endl;
    }
*/

    for (int iii = 0; iii < S_temp.size(); iii++)
      sort_edge(&S_temp.at(iii));
    sort_edge_sets (&S_temp);

    for (int iii = 0; iii < Q.size(); iii++)
      sort_edge(&Q.at(iii));
    for (int iii = 0; iii < R_star.size(); iii++)
      sort_edge(&R_star.at(iii));
    sort_edge_sets (&Q);
    sort_edge_sets (&R_star);
    merge_edge_set(&Q,&R_star,&Q_temp);

    Q.clear();
    S.clear();

    std::vector<std::vector<Edge> > set_temp;
    set_temp.clear();
    set_temp.push_back(paths_e[i]);

    merge_edge_set(&Q_temp,&R,&Q);

    for (int iii = 0; iii < Q_temp.size(); iii++)
      sort_edge(&Q_temp.at(iii));
    sort_edge_sets (&Q_temp);

    merge_edge_set(&Q_temp,&R,&Q_temp2);

    for (int iii = 0; iii < Q_temp2.size(); iii++)
      sort_edge(&Q_temp2.at(iii));
    sort_edge_sets (&Q_temp2);

    merge_edge_set(&Q_temp2,&set_temp,&Q);

    for (int iii = 0; iii < Q.size(); iii++)
      sort_edge(&Q.at(iii));
    sort_edge_sets (&Q);

    merge_edge_set(&S_temp,&set_temp,&S);

    for (int iii = 0; iii < S.size(); iii++)
      sort_edge(&S.at(iii));
    sort_edge_sets (&S);

  } // for_1 ends

/*
  for (int iii = 0; iii < S.size(); iii++) {
    std::cout << "S-Cycle " << iii << ", " << "size: " << 
    S[iii].size() << std::endl;
    for (it_e = S[iii].begin(); it_e != S[iii].end(); ++it_e) {
      std::cout << it_e->first << " ";
    }
    std::cout << std::endl;
  }
*/
} // find_gibbs ends

void PatonFinder::find_paton() {
  std::set<Vertex> X = V;
  tree.clear();

  std::set<Vertex> intersect;
  set_intersection(X.begin(),X.end(),T.begin(),T.end(),
              std::inserter(intersect,intersect.begin()));

  //std::cout << "Size: " << intersect.size();
  //std::cout << "max: " << intersect.max_size();
  //std::cout << std::endl;

  while(!intersect.empty()) {
    //std::cout << "Checking " << *intersect.begin() << std::endl;
    for (std::tie(ai, ai_end) = boost::adjacent_vertices(*intersect.begin(), g); ai != ai_end; ++ai) {
      //std::cout << *ai <<  " ";
      if (T.find(*ai) == T.end()) {
        T.insert(*ai);
        boost::add_edge(*intersect.begin(), *ai, tree);
      }
      else {
        //std::cout << "Cycle Found\n" << std::endl;
        find_nonrecursive(*intersect.begin(), *ai);
        //std::cout << " " << *intersect.begin() << std::endl;
      }
      remove_edge(*intersect.begin(), *ai, g);
    }
    //std::cout << std::endl;
    X.erase(X.find(*intersect.begin()));
    intersect.clear();
    set_intersection(X.begin(),X.end(),T.begin(),T.end(),
              std::inserter(intersect,intersect.begin()));
  }

  for (int i = 0; i < paths.size(); i++) {
    std::cout << "Cycle " << i << ": ";
    for (int j = 0; j < paths.at(i).size(); ++j) {
      std::cout << paths[i][j] << " ";
    }
    std::cout << std::endl;
  }

  // Detect the cycles:
  // Compare the elements of cycles. If there is a match, the rest of the
  // matches should be continuous in both cycles. 
  // Record the indices of the first and the last entities. In case there 
  // is a break in matches in a cycle and the first entity was matched, 
  // check the end entities in the cycle.
  // Remove the matching parts, join the rest. The cycles should be in 
  // the same sense.
  // m3 00 01 m1 m2    00 01 m1 m2 m3
  // 02 03 m3 m2 m1 => 02 03 m3 m2 m1 => 00 01 02 03
  for (int i = 0; i < paths.size()-1; i++) {
    for (int j = i+1; j < paths.size(); j++) {
      intersect.clear();
      set_intersection(paths[i].begin(),paths[i].end(),
                       paths[j].begin(),paths[j].end(), 
                       std::inserter(intersect,intersect.begin()));
      if(!intersect.empty()) {
      }
    }
  }

  std::ofstream fout("figs/tree.dot");
  write_graphviz(fout, tree);
}

    // Given a graph, first find the elementary circuits by Paton algorithm.
    // Using the gibbs algorithm, find all circuits as edge lists. 
    // Convert the edge lists into vertex lists.
void PatonFinder::find_circuits() {
  find_paton();
  vert_2_edge();

  for (int i = 0; i < paths_e.size(); i++) {
    for (int j = 0; j < paths_e[i].size(); j++) {
      std::cout << paths_e[i][j].first << " ["
                << paths_e[i][j].second << "] ";
    }
    std::cout << std::endl;
  }

  std::set<Edge >::iterator it_e;
/*
  for (int i = 0; i < paths_eset.size(); i++) {
    std::cout << "Cycle " << i << ", " << "size: " << 
    paths_eset[i].size() << std::endl;
    for (it_e = paths_eset[i].begin(); it_e != paths_eset[i].end(); ++it_e) {
      std::cout << it_e->first << " ";
    }
    std::cout << std::endl;
  }

  Edge e_cur = *paths_eset[0].begin();
  Vertex u=source(e_cur.first,g);
  Vertex v=target(e_cur.first,g);
  std::cout << "Internal: " << e_cur.first << " " << u << ", " << v;
  std::cout << std::endl;
*/

  find_gibbs();
}


// Return a pointer to the circuits.
const std::vector< std::vector<Edge> >* PatonFinder::get_circuits() {
  return &S;
}
// Return a pointer to the associated subgraphs.
const std::vector< std::vector<Vertex> >* PatonFinder::get_disjoint() {
  return &s_graphs;
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
cell_graph::cell_graph(){
}


// ---------------------------------------------------------
// Call the routines of the PatonFinder.
// ---------------------------------------------------------

// Find the circuits within the loaded graph. Get the disjoint graphs 
// separated by the circuits.

void cell_graph::find_circuits() {
  PF = new PatonFinder (&g);
  PF->find_circuits();
  PF->create_sub();

  S = PF->get_circuits();
  s_graphs = PF->get_disjoint();

  std::ofstream fout2("output/23adj.dot");
  write_graphviz(fout2, g, boost::make_label_writer(&cell_lbls[0]));

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
void cell_graph::get_23adj(int tag_0cell) {
  cb->get_conn_dim(3, 0, tag_0cell-1, &cells3);
  cb->get_conn_dim(2, 0, tag_0cell-1, &cells2);

  int cell3_nbr = cells3.n_fill;
  int cell2_nbr = cells2.n_fill;

  g.clear();
  cell_lbls.clear();
  cell_type.clear();
  cell_lbls.reserve(cell3_nbr+cell2_nbr);
  cell_lbls.resize(cell3_nbr+cell2_nbr);
  cell_type.reserve(cell3_nbr+cell2_nbr);
  cell_type.resize(cell3_nbr+cell2_nbr);

  for (int i = 0; i < cells3.n_fill; i++) {
    cell_type[i] = 3;
    cell_lbls[i] = "3c";
    std::string str2 = std::to_string(cells3.conn[i]+1);
    cell_lbls[i].append(str2);
  }

  for (int i = 0; i < cells2.n_fill; i++) {
    cell_type[i] = 2;
    cell_lbls[cell3_nbr+i] = "2c";
    std::string str2 = std::to_string(cells2.conn[i]+1);
    cell_lbls[cell3_nbr+i].append(str2);
  }

  for (int i = 0; i < cells3.n_fill; i++) {
    struct ent_conn e_con;
    cb->get_conn(3, cells3.conn[i], &e_con);
    for (int j = 0; j < cells2.n_fill; j++) {
      if (e_con.chk_ent(cells2.conn[j]))
        boost::add_edge(i, cell3_nbr+j, g);
    }
  }
  find_circuits();
}

// Test function to print the current graph.
void cell_graph::print_graph() {
  std::ofstream fout("output/graph.dot");
  fout << "graph A {\n"
    << "  rankdir=LR\n"
    << "ratio=\"fill\"\n"
    << "edge[style=\"bold\"]\n" << "node[shape=\"circle\"]\n";

  graph::edge_iterator it, end;
  std::tie(it, end) = boost::edges(g);

  for (std::tie(it, end) = edges(g); it != end; ++it)
    fout << source(*it, g) << " -- " << target(*it, g) << "\n";

  fout << "}\n";
  //std::ofstream fout2("output/graph2.dot");
  //write_graphviz(fout2, g);
}

// ---------------------------------------------------------
// Return the pointer to the circuits and the disjoint graphs.
// ---------------------------------------------------------

// Return a pointer to the circuits.
const std::vector< std::vector<Edge> >* cell_graph::get_circuits() {
  return &S;
}
// Return a pointer to the associated subgraphs.
const std::vector< std::vector<Vertex> >* cell_graph::get_disjoint() {
  return &s_graphs;
}

// ---------------------------------------------------------
// Constructor.
// ---------------------------------------------------------
cell_graph::cell_graph(struct cell_base* c_base) {
  cb = c_base;
  g.clear();
}

cell_graph::~cell_graph() {
  //delete PF;
}

