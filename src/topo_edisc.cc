#include <vector>
#include <map>
#include <deque>
#include <algorithm>    
#include <functional>
#include <cstring>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>
#include <gmi_null.h>

#include <ma.h>

#include "topo_extinfo.h"
#include "topo_geom.h" // Also contains pi.
#include "topo_field.h"

#include "topo_manip.h"
#include "topo_disp.h"

#include "topo_topo.h"
#include "topo_graph.h"

#include "topo_lens.h"
#include "topo_glens.h"

#include "topo_edisc.h"

#include "topo_entlist.h"

#include "topo_tmesh.h"

#include "topo_energy.h"
#include "topo_ma.h"
#include "topo_write.h"

// Given the vertex, return the adjacent triangles.
int lookup_tet_surf [4][3] = {{0,1,3},{0,2,1},{0,3,2},{1,2,3}};
// Given the vertex, return the adjacent edges.
int lookup_tet_edge [4][3] = {{0,2,3},{0,1,4},{1,2,5},{3,4,5}};
// Given the triangle, return the adjacent vertices.
int lookup_tet_surf_v [4][3] = {{0,2,1},{0,1,3},{1,2,3},{0,3,2}};
// Given the vertex, return the triangle across.
int lookup_tet_surf_x [4] = {2, 3, 1, 0};
// Given the triangle, return the vertex across.
int lookup_tet_x_surf [4] = {3, 2, 0, 1};
// Given the vertex, return the edges across within the tetrahedron.
int lookup_v_tet_e_x [4][3] = {{1,4,5},{2,5,3},{0,4,3},{0,1,2}};
// Given the vertex, and an edge bounded by the vertex, return the triangle 
// across within the tetrahedron.
int lookup_v_tet_e_t_x [4][6] = {{3,-1,1,0,-1,-1},{2,1,-1,-1,0,-1}, 
                                 {-1,3,2,-1,-1,0},{-1,-1,-1,2,3,1}};
// Given the vertex, and a tri bounded by the vertex, return the two tris 
// across the tri within the tetrahedron.
int lookup_v_tet_t_t_x [4][4][2] = {{{1,3},{0,3},{-1,-1},{0,1}},
                                    {{1,2},{0,2},{0,1},{-1,-1}},
                                    {{2,3},{-1,-1},{0,3},{0,2}},
                                    {{-1,-1},{2,3},{1,3},{1,2}}};
// Given the vertex, return the adjacent edges within the triangle.
int lookup_tri_ed [3][2] = {{2,0},{0,1},{1,2}};
// Given the vertex, return the edges across within the triangle.
int lookup_v_tri_e_x [3] = {1,2,0};
// Given the edge, return the vertex across within the triangle.
int lookup_v_tri_x_e [3] = {2,0,1};
// Given a vertex on a tri, return the indices of the other two vertices.
int lookup_v_tri_x_x [3][2] = {{1, 2}, {0, 2}, {0, 1}};
// Given the vertex, and an edge return the other edge adjacent to the vertex 
// within the triangle.
int lookup_v_tri_n_e [3][3] = {{2,-1,0},{1,0,-1},{-1,2,1}};

// Given two vertex indices, return the index of the other one. Not necessary 
// but cleaner.
int lookup_triverts2 [3][3] = {{-1,2,1},{2,-1,0},{1,0,-1}};

// Given two triangle indices, find the edge not contained within one of the
// triangles.
int lookup_tri_edges [4][4] = {{-1,5,3,4},{5,-1,2,1},{3,2,-1,0},{4,1,0,-1}};


// ----------------------------------------------------------
// Operations used in finding triangles of discs.
// ----------------------------------------------------------

// Given a mesh, a center vertex and two edges that are bounded by this vertex,
// find the shortest path of triangles passing through a single 3cell. Assumes
// the edges to be connected over entities belonging to that 3cell.

vd_disc& vd_disc::operator=( const vd_disc& obj ) {

  tri = obj.tri;
  edge = obj.edge;
  v_pos = obj.v_pos; 
  //elem_top = obj.elem_top;
  //elem_bottom = obj.elem_bottom;
  broken = obj.broken;
  t_em = obj.t_em;
  t2_em = obj.t2_em;
  e_em = obj.e_em;
  e_t_map1 = obj.e_t_map1;
  e_t_map2 = obj.e_t_map2;
  t_e_map1 = obj.t_e_map1;
  t_e_map2 = obj.t_e_map2;
  return *this;
}

vd_disc::vd_disc () : m(NULL), v_ctr(NULL), 
    tri(0), edge(0), v_pos(0, apf::Vector3(0,0,0)), 
    //elem_top(0), elem_bottom(0),
    broken(false), t_em(0), e_em(0), t2_em(0), 
    e_end1(NULL), e_end2(NULL),
    e_t_map1{}, e_t_map2{}, t_e_map1{}, t_e_map2{} {

}

// Print the content of the disc.
void vd_disc::print_cont() {
  std::cout << "Vd_disc: Mesh: " << m << std::endl;
  std::cout << "Disc is " << broken << std::endl;

  std::cout << "Tri: ";
  for (int i = 0; i < tri.size(); i++) {
    std::cout << tri.at(i) << " ";
  }
  std::cout << std::endl;

  std::cout << "Edge: ";
  for (int i = 0; i < edge.size(); i++) {
    std::cout << edge.at(i) << " ";
  }
  std::cout << std::endl;

  //std::cout << "Elem_top: ";
  //for (int i = 0; i < elem_top.size(); i++) {
  //  std::cout << elem_top.at(i) << " ";
  //}
  //std::cout << std::endl;

  //std::cout << "Elem_bottom: ";
  //for (int i = 0; i < elem_bottom.size(); i++) {
  //  std::cout << elem_bottom.at(i) << " ";
  //}
  //std::cout << std::endl;

  std::cout << "Cell_memb of tri: ";
  for (int i = 0; i < t_em.size(); i++) {
    std::cout << t_em.at(i) << " ";
  }
  std::cout << std::endl;

  std::cout << "Cell_memb of tri-2: ";
  for (int i = 0; i < t2_em.size(); i++) {
    std::cout << t2_em.at(i) << " ";
  }
  std::cout << std::endl;

  std::cout << "Cell_memb of edge: ";
  for (int i = 0; i < e_em.size(); i++) {
    std::cout << e_em.at(i) << " ";
  }
  std::cout << std::endl;
}

// Reset the disc.
void vd_disc::clear() {

  tri.clear();
  edge.clear();
  v_pos.clear();

  broken = false;

  //elem_top.clear();
  //elem_bottom.clear();
  t_em.clear();
  e_em.clear();

  t2_em.clear();

  e_t_map1.clear();
  e_t_map2.clear();
  t_e_map1.clear();
  t_e_map2.clear();

}

vd_disc::~vd_disc() {
  clear();
}

// Copy the triangles and edges without sorting. They only need to be unique.
void cp_tri(const std::vector<apf::MeshEntity*>* es_tri, 
                                    vd_disc* disc) {

  apf::Up up;

  disc->broken = false;

  disc->tri.clear();
  disc->edge.clear();

  if(es_tri->size() > 0) {
    int tri_sz = es_tri->size();
    disc->tri.resize(tri_sz);
    disc->edge.reserve(2*tri_sz);

    // If the number of unique edges is equal to that of the triangles, the disc
    // is not broken. If they are +1, the disc is broken. If +2, the disc is two
    // fins. 

    for(int i = 0; i < es_tri->size(); i++) {
      disc->tri.at(i) = es_tri->at(i);
      apf::Downward down;
      apf::Downward v;
      disc->m->getDownward(es_tri->at(i), 1, down);
      disc->m->getDownward(es_tri->at(i), 0, v);
      int v0 = findIn(v, 3, disc->v_ctr);

      disc->edge.push_back(down[lookup_tri_ed[v0][0]]);
      disc->edge.push_back(down[lookup_tri_ed[v0][1]]);
    }

    std::sort(disc->edge.begin(), disc->edge.end());
    std::vector<apf::MeshEntity*>::iterator it;
    it = std::unique (disc->edge.begin(), disc->edge.end());
    disc->edge.resize(std::distance(disc->edge.begin(),it));

    if(disc->edge.size() > disc->tri.size() )
      disc->broken = true;

    disc->t_em.resize(disc->tri.size());
    disc->e_em.resize(disc->edge.size());
    disc->t2_em.resize(disc->edge.size());

    //disc->elem_top.resize(disc->tri.size());
    //disc->elem_bottom.resize(disc->tri.size());

    disc->v_pos.resize(disc->edge.size());

    for(int i = 0; i < disc->v_pos.size(); i++) {
      apf::MeshEntity* v_curr = getEdgeVertOppositeVert(disc->m, 
                                        disc->edge.at(i), disc->v_ctr);
      disc->m->getPoint(v_curr, 0, disc->v_pos.at(i));
      //disc->v_pos.at(i) = disc->v_pos.at(i) - v_pos;
    }
    for (int j = 0; j < disc->edge.size(); j++) {
      apf::MeshEntity* e_curr = disc->edge.at(j);
      disc->m->getUp(e_curr, up);
      bool found_1 = false;
      bool found_2 = false;
      for (int k = 0; k < up.n; k++) {
        std::vector<apf::MeshEntity*>::iterator it_e;
        it_e = std::find(disc->tri.begin(), 
                         disc->tri.end(), up.e[k]);
        // If on the disc, replace it with the new edge for the top disc.
        if(it_e != disc->tri.end()) {
          if(!found_1) {
            found_1 = true;
            disc->e_t_map1[e_curr] = std::distance(disc->tri.begin(), it_e)+1;
            if(disc->t_e_map1[*it_e] == 0) {
              disc->t_e_map1[*it_e] = j+1;
            }
            else {
              assert(disc->t_e_map2[*it_e] == 0);
              disc->t_e_map2[*it_e] = j+1;
            }
          }
          else {
            assert(!found_2);
            found_2 = true;
            disc->e_t_map2[e_curr] = std::distance(disc->tri.begin(), it_e)+1;
            if(disc->t_e_map1[*it_e] == 0) {
              disc->t_e_map1[*it_e] = j+1;
            }
            else {
              assert(disc->t_e_map2[*it_e] == 0);
              disc->t_e_map2[*it_e] = j+1;
            }
            k = up.n;
          }
        }
      }
    }
  }
}

/*
vd_fin::vd_fin (apf::Mesh2* m_in, std::vector<int>* path_in, 
            apf::MeshEntity* tri_in, apf::MeshEntity* edge_in) {
  m = m_in;
  path = path_in;
  // By definition, it will be broken.
  broken = true;
  tri_first = tri_in;
  edge_first = edge_in;

  // TODO Not sure if it is necessary, performance related. Usage of set? 
  tri.push_back(tri_in);
  // The disc edges.
  edge.push_back(edge_in);

  // This section crawls to find the disc.

  apf::ModelEntity* mdl = m->toModel(tri_in);

  int type_ctr = m->getModelType(mdl);
  int tag_ctr = m->getModelTag(mdl);

}
*/
// Vd elens is a container for the discs used in lens expansion.
vd_elens::vd_elens () :
  m(NULL), discs(0), init_flag(false), null_flag(false), 
  es_vert(0), es_edge(0),  es_surf(0),  es_elem(0), vert_ctr(NULL),
  pos(0,0,0), vert_ctr_new(0), edge_ctr_new_t(0),  edge_ctr_new_b(0),  
  vert_sp_new(0), c1(0, std::vector<int >(0)),
  c1_mult(0),
  edge_sp_new(0, std::vector<apf::MeshEntity*>(0)),
  tri_disc_new(0, std::vector<apf::MeshEntity*>(0)), c3_f(-1), c3_e(-1),

  e_3c_t_flag(false), e_3c_b_flag(false), e_3c_t(NULL), e_3c_b(NULL), 
  v_2c(NULL),
  new_cell0_id(0), new_cell1_id(0), new_cell2_id(-1), c_base_curr(NULL),

  vt(0, apf::Vector3(0,0,0)), split_vt(0, apf::Vector3(0,0,0)),
  slices(0), slice_map{}, 
  create_map(0, std::map<apf::MeshEntity*, apf::MeshEntity* >{}),
  tri_tet_map(0, std::map<apf::MeshEntity*, apf::MeshEntity*> {}),
  tri_map(0, std::map<apf::MeshEntity*, apf::MeshEntity* >{}),
  split_map(0, std::map<apf::MeshEntity*, apf::MeshEntity* >{}),
  create_map_rest{} {
}

vd_elens::vd_elens (apf::Mesh2* m_in) :
  m(NULL), discs(0), init_flag(false), null_flag(false), 
  es_vert(0), es_edge(0),  es_surf(0),  es_elem(0), vert_ctr(NULL),
  pos(0,0,0), vert_ctr_new(0), edge_ctr_new_t(0),  edge_ctr_new_b(0),  
  vert_sp_new(0), c1(0, std::vector<int >(0)),
  edge_sp_new(0, std::vector<apf::MeshEntity*>(0)),
  tri_disc_new(0, std::vector<apf::MeshEntity*>(0)), c3_f(-1), c3_e(-1),

  e_3c_t_flag(false), e_3c_b_flag(false), e_3c_t(NULL), e_3c_b(NULL), 
  v_2c(NULL),
  new_cell0_id(0), new_cell1_id(0), new_cell2_id(-1), c_base_curr(NULL),

  vt(0, apf::Vector3(0,0,0)), split_vt(0, apf::Vector3(0,0,0)),
  slices(0), slice_map{}, 
  create_map(0, std::map<apf::MeshEntity*, apf::MeshEntity* >{}),
  tri_tet_map(0, std::map<apf::MeshEntity*, apf::MeshEntity*> {}),
  tri_map(0, std::map<apf::MeshEntity*, apf::MeshEntity* >{}),
  split_map(0, std::map<apf::MeshEntity*, apf::MeshEntity* >{}),
  create_map_rest{} {

  m = m_in;
}

void vd_elens::set_m(apf::Mesh2* m_in) {
  clear();
  m = m_in;
}

// ----------------------------------------------------------
// Operations used in finding triangles of discs.
// ----------------------------------------------------------

// Given a mesh, a center vertex and two edges that are bounded by this vertex,
// find the shortest path of triangles passing through a single 3cell. Assumes
// the edges to be connected over entities belonging to that 3cell.
bool vd_elens::vd_mesh_find_short(apf::MeshEntity* e_1, apf::MeshEntity* e_2,
                           std::vector<apf::MeshEntity*>* es_tri, int c3_id) {

  apf::Up up;
  apf::Downward v;
  apf::Downward e;

  es_tri->reserve(es_tri->size()+es_surf.size());

  // The map of previous entities.
  std::map<apf::MeshEntity*,apf::MeshEntity* > e_prev{};
  std::map<apf::MeshEntity*,bool> e_burn{};

  // The sum of the internal angles of the central vertex corner of the passed 
  // triangles. At each edge, for each adjacent edge not visited, add the new 
  // edge and add set the sum of the internal angle so far. For those visited 
  // before, compare the current to the the value before. If lower, replace the
  // previous edge. 
  // While the end edge has not been found or there are branches with sum_ang 
  // less than the end edge, try branches. If the minimum of the current 
  // branches is more than the end edge, the shortest path has been found.
  std::map<apf::MeshEntity*,double> int_ang{};
  std::map<apf::MeshEntity*,double> sum_ang{};
  double ang_min = 6.28;
  double ang_end = 6.28;

  for(int i = 0; i < es_surf.size(); i++) {
    int_ang[es_surf.at(i)] = vd_int_angle(m, es_surf.at(i), vert_ctr);
  }

  std::vector<std::vector <apf::MeshEntity* > > e_fire(2, 
                                          std::vector <apf::MeshEntity* >(0));

  // Reserve space for all edges.
  m->getUp(vert_ctr, up);
  e_fire[0].reserve(up.n);
  e_fire[1].reserve(up.n);
  e_fire[0].push_back(e_1);

  sum_ang[e_1] = 0.;

  bool found = false;
  bool faulty = false;
  int i_curr = 0;
  int i_next = 1;

  apf::MeshEntity* e_curr;
  apf::MeshEntity* e_next;

  while((!found or ang_min < ang_end) and e_fire[i_curr].size() > 0) {

    // A triangle corner angle is always less than pi.
    double ang_min_temp = ang_min + 3.15;
    while(e_fire[i_curr].size() > 0) {

      e_curr = e_fire[i_curr].back();

      // Get the triangles adjacent to the current edge.
      m->getUp(e_curr, up);
      for(int i = 0; i < up.n; i++) {
        if (!e_burn[up.e[i]]) {
          m->getDownward(up.e[i], 0, v);
          m->getDownward(up.e[i], 1, e);

          int v1 = findIn(v, 3, vert_ctr);
          int e1 = findIn(e, 3, e_curr);

          // Get the other edge adjacent to the center vertex.
          if(lookup_tri_ed[v1][0] == e1) {
            e_next = e[lookup_tri_ed[v1][1]];
          }
          else {
            e_next = e[lookup_tri_ed[v1][0]];
          }
          //std::cout << "Tri " << up.e[i]
          //         << std::endl;
          //for (int k = 0; k < 3; k++) {
          //  m->getDownward(e[k], 0, v);
          //  std::cout << "Edge " << e[k] << " v1 " << v[0] << " v2 " << v[1]
          //            << std::endl;
          //  e1 = findIn(v, 3, vert_ctr);
          //  assert(e1 > -1);
          //}

          // Check if the edge is burned.

          // The edge can't belong to another disc.
          if (!(slice_map[e_next] < 0)) {
            // Visiting the edge for the first time.
            if (!e_burn[e_next]) {
              // Destination reached.
              sum_ang[e_next] = sum_ang[e_curr]+int_ang[up.e[i]];
              if (e_next == e_2) {
                e_burn[e_next] = true;
                e_prev[e_next] = up.e[i];
                e_prev[up.e[i]] = e_curr;

                found = true;
                i = up.n;
              }
              // Check if the edge is the end edge or member of a 3cell.
              else {
                int c_dim = m->getModelType(m->toModel(e_next));
                if(c_dim == 3) {
                  if(c3_id == 0) {
                    if(sum_ang[e_next] < ang_end) {
                      e_fire[i_next].push_back(e_next);

                      e_burn[e_next] = true;
                      e_prev[e_next] = up.e[i];
                      e_prev[up.e[i]] = e_curr;
                    }
                  }
                  else if(c3_id == m->getModelTag(m->toModel(e_next)) ) {
                    if(sum_ang[e_next] < ang_end) {
                      e_fire[i_next].push_back(e_next);

                      e_burn[e_next] = true;
                      e_prev[e_next] = up.e[i];
                      e_prev[up.e[i]] = e_curr;
                    }
                  }
                }
              }

              if (e_next == e_2) {
                ang_end = sum_ang[e_next];
              }
              if(sum_ang[e_next] < ang_min_temp) {
                ang_min_temp = sum_ang[e_next];
              }
            }
            else {
              // If this branch is geometrically shorter, replace the previous
              // entities.
              if(sum_ang[e_curr]+int_ang[up.e[i]] < sum_ang[e_next]) {

                e_prev[e_next] = up.e[i];
                e_prev[up.e[i]] = e_curr;

                if (e_next == e_2) {
                  ang_end = sum_ang[e_next];
                }
                else {
                  // In order to update the nodes on the path, add the node back
                  // to the list.
                  e_fire[i_next].push_back(e_next);
                }

                if(sum_ang[e_next] < ang_min_temp) {
                  ang_min_temp = sum_ang[e_next];
                }
              }
            }

          }
          else {
          }
          e_burn[up.e[i]] = true;

        }

      }

      //if(found) {
      //  e_fire[i_curr].clear();
      //  e_fire[i_next].clear();
      //}
      //else {
      e_fire[i_curr].pop_back();
      //}
    }
    ang_min = ang_min_temp;

    i_next = i_curr;
    i_curr = (i_curr + 1) % 2;

    if(!found and e_fire[i_curr].empty() ) {
      faulty = true;
      found = true;
    }
  }

  if(faulty)
    found = false;

  if(found) {
    e_next = e_2;
    // Follow the line of predecessing entities back to the starting edge.
    while(e_next != e_1) {
      if(m->getType(e_next) == 2) {
        es_tri->push_back(e_next);
      }
      assert(e_prev[e_prev[e_next]] != e_next);
      e_next = e_prev[e_next];
    }
  }

  return found;
}

// Check the slice tetrahedro for negative volume entities after the insertion.
long int vd_elens::vd_chk_neg_slice() {

  apf::MeshIterator* it = m->begin(3);

  long int n = 0;
  double meas = 0;
  apf::MeshEntity* elem;
  while(elem = m->iterate(it)) {


  //cell_id = m->getModelTag(vert_ctr_em);

    apf::MeshElement* ee = createMeshElement(m, elem);
    meas = measure(ee);
    destroyMeshElement(ee);

    if (meas < 0 or std::fabs(meas) < std::numeric_limits<double>::min()) {
      apf::ModelEntity* mdl = m->toModel(elem);
      std::cout << "Elem " << elem << " " << m->getModelType(mdl) << "c"
            << m->getModelTag(mdl) << " is inverted with vol "
            << meas << std::endl;
      n++;
    }
  }
  m->end(it);

  return n;
}


// Print the content of the lens.
void vd_elens::print_cont() {

  for (int i = 0; i < discs.size(); i++) {
    std::cout << "Disc " << i << ":" << std::endl;
    discs.at(i).print_cont();
  }
}

// Reset the contents.
void vd_elens::clear_old_ent() {

  es_vert.clear();
  es_edge.clear();
  es_surf.clear();
  es_elem.clear();
}

void vd_elens::clear_discs() {
  // Clear the old entity discs
  for (int i = 0; i < discs.size(); i++)
    discs.at(i).clear();
  discs.clear();

  // Clear slices associated with the discs. Two per disc.
  slices.clear();

}

// Clear new vertices
void vd_elens::clear_new_vert() {

  vert_ctr_new.clear();
  vert_sp_new.clear();
  vt.clear();
  split_vt.clear();
}

// Clear create maps
void vd_elens::clear_create_maps() {

  // Clear the maps from old to new entities
  for (int i = 0; i < create_map.size(); i++)
    create_map.at(i).clear();
  create_map.clear();

  // Clear the maps from old to new entities
  for (int i = 0; i < tri_tet_map.size(); i++)
    tri_tet_map.at(i).clear();
  tri_tet_map.clear();

  for (int i = 0; i < tri_map.size(); i++)
    tri_map.at(i).clear();
  for (int i = 0; i < split_map.size(); i++)
    split_map.at(i).clear();

  split_map.clear();
  tri_map.clear();
}

// Reset the lens.
void vd_elens::clear() {
  // Clear old entities
  clear_old_ent();
  clear_discs();

  // Clear topology
  for (int i = 0; i < c1.size(); i++)
    c1.at(i).clear();
  c1.clear();
  c1_mult.clear();

  // Clear old entity slice map.
  slice_map.clear();

  clear_new_vert();

  clear_create_maps();

  // Clear new top and bottom entities
  edge_ctr_new_t.clear();
  edge_ctr_new_b.clear();


  for (int i = 0; i < edge_sp_new.size(); i++)
    edge_sp_new.at(i).clear();
  for (int i = 0; i < tri_disc_new.size(); i++)
    tri_disc_new.at(i).clear();

  new_cell0_id.clear();
  new_cell1_id.clear();

  create_map_rest.clear();

  init_flag = false;
}

void vd_elens::null_slice_map() {
  for(int i = 0; i < es_elem.size(); i++) {
    slice_map[es_elem.at(i)] = 0;
  }
  for(int i = 0; i < es_surf.size(); i++) {
    slice_map[es_surf.at(i)] = 0;
  }
  for(int i = 0; i < es_edge.size(); i++) {
    slice_map[es_edge.at(i)] = 0;
  }
  for(int i = 0; i < es_vert.size(); i++) {
    slice_map[es_vert.at(i)] = 0;
  }
  null_flag = true;

  //for(int i = 0; i < slices.size(); i++) {
  //  slices.at(i) = 0;
  //}

  slices.clear();

}

vd_elens::~vd_elens () {
  clear();
}

void vd_disc_cut::clear() {

  c2_edge_new.clear();
  c2_edge_flag.clear();

  c2_cuts.clear();
  c2_edges.clear();
  c2_cut_pos.clear();
  ent_burn.clear();
  ents.clear();

}
vd_disc_cut::vd_disc_cut() : 
  circ_tup(NULL), ng_pt(NULL), mdl_curr(NULL),
  c2_edge_new(0), c3_mdl(0), c2_edge_flag(0), int_tol(0.1),
  c2_cuts(0), c2_edges(0), c2_cut_pos(0), ents(0), ent_burn{} {
}
vd_disc_cut::~vd_disc_cut() {
  clear();
}

void vd_edisc::clear_map() {
  preid_map.clear();
  actid_map.clear();
  main2pre_map.clear();
  trial2pre_map.clear();

  pre2main_map.clear();
  pre2trial_map.clear();

  wing_map.clear();
}

void vd_edisc::clear_ext() {

  ext_new.clear();
  ext_slice.clear();
  ext_pc.clear();
  ext_cor.clear();
  ext_0cell = false;

  for(int i = 0; i < ext_proj_dir.size(); i++) {
    ext_proj_dir.at(i).clear();
  }
  ext_proj_dir.clear();
}

void vd_edisc::clear() {
  clear_ext();

  c2_edge.clear();

  for(int i = 0; i < elem_slice.size(); i++) {
    elem_slice.at(i).clear();
  }
  elem_slice.clear();

  elem_burn.clear();
  tri_burn.clear();

  elem_ofire.clear();

  es_vert.clear();
  es_edge.clear();
  es_surf.clear();
  es_elem.clear();

  es_vert_act.clear();
  es_edge_act.clear();
  es_surf_act.clear();
  es_elem_act.clear();

  vert.clear();
  vert_trial.clear();

  shell2_map.clear();
  shell2_list.clear();

  w1.clear();
  w2.clear();
  e1.clear();
  e2.clear();
  w1_exp.clear();
  w2_exp.clear();

  clear_vel_map();
  clear_map();
}


// Change the active mesh in a centralized manner.
void vd_edisc::act_precond() {
  assert(precond_ex);
  m_act = m_precond;
  vert_ctr_act = vert_ctr_precond;
  init_elens();

  //c_base_init->copy_local(*c_base, m_main->getModelTag(vert_ctr_em) - 1);
  //c_base_init = c_base;
  c_base_act = c_base;

  e_list_act->change_mesh(m_act, c_base_act);

  f_calc.reload_cb(c_base_act);
  f_calc.refresh(m_act, c_base_act, e_list_act);

  e_lens.c_base_curr = c_base_act;
}

void vd_edisc::act_trial(){
  assert(trial_ex);
  m_act = m_trial;
  vert_ctr_act = vert_ctr_cp;
  init_elens();

  //c_base_init->copy_local(*c_base, m_main->getModelTag(vert_ctr_em) - 1);
  c_base_act = c_base;

  e_list_act->change_mesh(m_act, c_base_act);

  f_calc.reload_cb(c_base_act);
  f_calc.refresh(m_act, c_base_act, e_list_act);

  e_lens.c_base_curr = c_base_act;
}

void vd_edisc::act_main(){
  assert(m_ex);
  m_act = m_main;
  vert_ctr_act = vert_ctr;
  init_elens();

  c_base_act = c_base;

  e_list_act->change_mesh(m_act, c_base_act);
  f_calc.reload_cb(c_base_act);
  f_calc.refresh(m_act, c_base_act, e_list_act);

  e_lens.c_base_curr = c_base_act;
}


// ---------------------------------------------------
// 1cell insertion with the current circuit.

// Constructor.
vd_edisc::vd_edisc(apf::Mesh2* m_in, cell_base* c, vd_entlist* el_in) 
  :   
    f_calc(), vd_cd(NULL), vd_3c(NULL), e_lens(), ng(), 
//e_sh_save(), 
//    f_calc(), vd_cd(), e_sh_save(), e_lens(), ng(), 
    spur(false), skip_en(false), calc_ext(false), calc_corner(false), ext_0cell(false),
    cell_flag(false),
    mov_flag(false), sub_vtk(true), verb_flag(true), drag_flag(false),
    shell2_map{}, wing_map{}, ext_new{}, ext_slice{}, ext_pc{}, ext_cor{},
    ext_proj_type{},
    wg_tag(WG_TYPE::EVOLVE),
    elem_burn{}, tri_burn{}, 
    slice_tris(0, std::vector<apf::MeshEntity* >(0)),
    slice_edges(0, std::vector<apf::MeshEntity* >(0)),
    preid_map{}, actid_map{}, main2pre_map{}, trial2pre_map{},
    pre2main_map{}, pre2trial_map{}, 
    e_list(NULL),
    shell2_list(0),
    sh_old(-1,-1), sh_ext(false),
    ext_proj_dir(0, std::vector<apf::Vector3> (0, apf::Vector3(0, 0, 0))),
    isotropic(true),
    m_ex(false), precond_ex(false), trial_ex(false), m_load(false),
    precond_load(false), trial_load(false),
    trial_type(0), trial_curr(0), ngon_curr(0),
    disc_load(false), pre_flag(false), ins_flag(false),

    circ_min(0), path_min(0), ngon_min(0), circ_sz(0), pt_sz(0),

    pos_old(0, 0, 0), v_avg(0, 0, 0),

    exp_dir_cos_th(0.001),
    energy_th(0), sh_used(false),
    vel_max(1), en_min_1c(0), en_min_2c(0), dt_curr(1), dt_max(1), t_total(0),

    cells3(), cells2(),

    cell_id(0), len_sh(-1), len_sp(-1), l_min(-1), len_edge(1),
    rho_rat(16), fudge_factor(1), rat_init(10),

    c3_edge(0), c2_edge(0),
    elem_slice(0, std::vector < apf::MeshEntity* > (0)),
    elem_ofire(0),

    es_vert(0), es_edge(0), es_surf(0), es_elem(0),
    es_vert_act(0), es_edge_act(0), es_surf_act(0), es_elem_act(0),
    vert(0), vert_trial(0),
    w1(0), w2(0, std::vector<double>(0)),
    w1_exp(0, std::vector<double>(0)), w2_exp(0, std::vector<std::vector<double> >(0, std::vector<double>(0))),
    e1(0), e2(0, std::vector<double>(0)) {

  c_base = c;

  //c_base_init = new cell_base(c_base->get_sz(0), 
  //                            c_base->get_sz(1),
  //                            c_base->get_sz(2),
  //                            c_base->get_sz(3));

  //c_base_init->print_ent();

  c_base_act = c_base;

  //c_ins = cins;
  m_main = m_in;
  e_list = el_in;

  // The empty model used in generating the trial meshes.
  //gmi_model* empty_mdl = gmi_load(".null");

  gmi_register_null();

  m_precond = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);
  m_trial = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);

  f_calc.reload_cb(c_base_act);

  isotropic = false;

  m_ex = true;
  precond_ex = true;
  trial_ex = true;

  m_load = true;
  precond_load = false;
  trial_load = false;

  //energy_func_tri = energy_func_tri_trial;
  //energy_func_tet = energy_func_tet_trial;

  cell_id = 1;
  trial_type = 0;
  trial_curr = 0;

  dt_max = 0;
  vd_cd = new vd_cell_det();
  vd_3c = new vd_3c_det();

  e_list_act = new vd_entlist();
}

// Destructor.
vd_edisc::~vd_edisc() {
  //clear_ext();
  clear();
  //delete c_base;
  if (precond_load) {
    //apf::destroyMesh(m_precond);
  }
  if (trial_load)
    apf::destroyMesh(m_trial);

  //delete c_base_init;
  delete vd_cd;
  delete vd_3c;

  delete e_list_act;
}

// Set the vertex belonging to the selected 0cell. Normally, this is to be 
// called internally. 
// This also creates the vertices in the trial mesh.
void vd_edisc::set_vertex(apf::MeshEntity* vertex) {
  vert_ctr = vertex;

  vert_ctr_em = m_main->toModel(vert_ctr);
  std::cout << "Vertex " << m_main->getModelType(vert_ctr_em)
            << "cell"<< m_main->getModelTag(vert_ctr_em) << std::endl;

  m_main->getPoint(vert_ctr, 0, pos_old);

  c_base->get_conn_dim_gmi(3, 0, cell_id, &cells3);
  c_base->get_conn_dim_gmi(2, 0, cell_id, &cells2);

  pre_flag = true;
  set_trial();
  if(!pre_flag) {
    std::cout << "Precondition not successful." << std::endl;
    return;
  }
  vtk_precond();

  //std::cout << "Verifying m_precond" << std::endl;
  //m_act->verify();
  //std::cout << "Initial:" << m_main << ", Precond:" << m_precond << std::endl;

  dt_max = calc_dt();
  f_calc.vdparam.adj_dt(dt_max);

  energy = 0;
  energy = calc_energy();
  //std::cout << "E per characteristic time length " << f_calc.vdparam.dt 
  std::cout << "E is " << energy << "." << std::endl;
}

// TODO revert these to function pointers. This was done to get rid of a 
// symptom where creating a new vd_edisc was causing stack overflows.
// The trial energy functions.
double vd_edisc::energy_func_tri(apf::Mesh2* m, apf::MeshEntity* tri) {

  int ent_type = m->getType(tri);
  int d = m->typeDimension[ent_type];
  int m_type = m->getModelType(m->toModel(tri));
  //std::cout << "Entity being measured: " << d << ", " << ent_type << ", "
  //          << m->getModelType(m->toModel(tri)) << ", "
  //          << m->getModelTag(m->toModel(tri)) << std::endl;

  //std::cout << "Area is " << meas << std::endl;

  if (d == 2 and m_type == 2) {
    //double meas = std::fabs(vd_area_out(e_lens.m, tri).getLength());
    double meas = std::fabs(vd_area_out(e_lens.m, tri).getLength())
                                            *f_calc.gam2(e_lens.m, tri);
    return meas;
  }
  //std::cout << "Returning 0" << std::endl;

  return 0;
}

double vd_edisc::energy_func_tet(apf::Mesh2* m, apf::MeshEntity* tet) {
  //int ent_type = m->getType(tet);
  //int d = m->typeDimension[ent_type];

  //if (d == 3 and m->getModelType(m->toModel(tet)) == 3) {
  //  apf::MeshElement* ee = createMeshElement(m, tet);
  //  return measure(ee);
  //}
  return 0;
}

// Copy entities from the original mesh.
// Copy entities from the original mesh to the preconditioned mesh.
// This one uses buildelement, apf implementation with actual mds support.
apf::MeshEntity* vd_edisc::copyElement_pre(
    apf::MeshEntity* original,
    apf::BuildCallback* cb)
{
  // Get the type in the original mesh.
  int type = m_main->getType(original);
  if (type==0)  {
    int e1 = findIn(&es_vert, es_vert.size(), original);
    assert(e1 > -1);
    // Return the corresponding vertex in trial mesh:
    return vert.at(e1);
  }

  int d = apf::Mesh::typeDimension[type];
  apf::Downward down;
  int nd = m_main->getDownward(original,d-1,down);

  //apf::Downward vert_down; // Downward vertices of the current entity.
  //m->getDownward(original,0,vert_down);

  for (int i=0; i < nd; ++i)
    down[i] = copyElement_pre(down[i],cb);

  nd = m_main->getDownward(original,0,down);
  for (int i=0; i < nd; ++i) {
    int e1 = findIn(&es_vert, es_vert.size(), down[i]);
    assert(e1 > -1);
    // Return the corresponding vertex in trial mesh:
    down[i] = vert.at(e1);
  }

  // The topology object in the actual mesh:
  apf::ModelEntity* mdl = m_main->toModel(original);
  int type_ori = m_main->getModelType(mdl);
  int tag_ori = m_main->getModelTag(mdl);

  //std::cout << type << "-Ent " << original 
  //          << " " << type_ori << "c" << tag_ori
  //          << std::endl;
  // The topology object in the trial mesh:
  mdl = m_precond->findModelEntity(type_ori,tag_ori);

  apf::MeshEntity* ent;
  ent = buildElement(m_precond, mdl, type, down);

  //vd_print_vert(m_precond, ent);

  return ent;
}

// Copy entities from the preconditioned mesh to the trial mesh. 
apf::MeshEntity* vd_edisc::copyElement(
    apf::MeshEntity* original,
    apf::BuildCallback* cb)
{
  // Get the type in the original mesh.
  int type = m_precond->getType(original);
  if (type==0)
  {
    std::vector<apf::MeshEntity*>::iterator it;
    it = std::find (vert.begin(), vert.end(), original);
    assert(it != vert.end());
    // Return the corresponding vertex in trial mesh:
    return vert_trial[std::distance(vert.begin(), it)];
  }
  int d = apf::Mesh::typeDimension[type];
  apf::Downward down;
  int nd = m_precond->getDownward(original,d-1,down);

  //apf::Downward vert_down; // Downward vertices of the current entity.
  //m->getDownward(original,0,vert_down);

  for (int i=0; i < nd; ++i)
    down[i] = copyElement(down[i],cb);

  nd = m_precond->getDownward(original,0,down);

  //apf::Downward vert_down; // Downward vertices of the current entity.
  //m->getDownward(original,0,vert_down);

  for (int i=0; i < nd; ++i) {
    int e1 = findIn(&vert, vert.size(), down[i]);
    assert(e1 != -1);
    down[i] = vert_trial.at(e1);
  }

  // The topology object in the actual mesh:
  apf::ModelEntity* mdl = m_precond->toModel(original);
  int type_ori = m_precond->getModelType(mdl);
  int tag_ori = m_precond->getModelTag(mdl);


  //std::cout << type << "-Ent " << original << " " 
  //          << type_ori << "c" << tag_ori  << std::endl;
  // The topology object in the trial mesh:
  mdl = m_trial->findModelEntity(type_ori,tag_ori);

  apf::MeshEntity* ent;
  ent = buildElement(m_trial, mdl, type, down, cb);
  //vd_print_vert(m_trial, ent);

  return ent;
}

// Copy element from preconditioned mesh to the trial mesh by replacing the given
// vertex with the replacement vertex. p2t is the map of vertices from 
// preconditioned to the trial mesh.
apf::MeshEntity* vd_edisc::copyElement_rpl(
    apf::MeshEntity* original,
    apf::MeshEntity* v_orig,
    apf::MeshEntity* v_repl,
    std::map<apf::MeshEntity*,apf::MeshEntity*> &p2t,
    std::map<apf::MeshEntity*,apf::MeshEntity*> &t2p,
    apf::BuildCallback* cb) {
  // Get the type in the original mesh.
  int type = m_precond->getType(original);
  if (type==0) {
    return NULL;
  }
  int d = apf::Mesh::typeDimension[type];
  apf::Downward down;
  int nd = m_precond->getDownward(original,d-1,down);

  //apf::Downward vert_down; // Downward vertices of the current entity.
  //m->getDownward(original,0,vert_down);

  for (int i=0; i < nd; ++i)
    down[i] = copyElement_rpl(down[i], v_orig, v_repl, p2t, t2p, cb);

  nd = m_precond->getDownward(original, 0, down);

  //apf::Downward vert_down; // Downward vertices of the current entity.
  //m->getDownward(original,0,vert_down);

  for (int i=0; i < nd; ++i) {
    if(down[i] == v_orig)
      down[i] = v_repl;
    else
      down[i] = p2t[down[i]];
    assert(down[i] != NULL);
  }

  // The topology object in the actual mesh:
  apf::ModelEntity* mdl = m_precond->toModel(original);
  int type_ori = m_precond->getModelType(mdl);
  int tag_ori = m_precond->getModelTag(mdl);

  mdl = m_trial->findModelEntity(type_ori,tag_ori);

  apf::MeshEntity* ent;
  ent = buildElement(m_trial, mdl, type, down, cb);
  p2t[original] = ent;
  t2p[ent] = original;
  return ent;
}

// Get the entity list around the new vertices.
void vd_edisc::crt_new_vert_list() {

  e_lens.vert_ctr_new.clear();
  e_lens.vert_sp_new.clear();
  e_lens.vt.clear();
  e_lens.split_vt.clear();
  double z[3] = {0,0,0};

  if(trial_type == 1) {
    e_lens.vert_ctr_new.resize(2);
    e_lens.edge_ctr_new_t.resize(1);
    e_lens.edge_ctr_new_b.resize(1);
    e_lens.vert_sp_new.resize(1);

    e_lens.vt.resize(2);
    e_lens.split_vt.resize(1);

    e_lens.vt.at(0).fromArray(z);
    e_lens.vt.at(1).fromArray(z);
    e_lens.split_vt.at(0).fromArray(z);

    e_lens.new_cell0_id.clear();
    e_lens.new_cell0_id.resize(2);
    e_lens.new_cell1_id.clear();
    e_lens.new_cell1_id.resize(1);

  }
  else {
    assert(trial_type == 2);

    e_lens.vert_ctr_new.resize(ng.ngons.at(ng.curr).size());
    e_lens.edge_ctr_new_t.resize(ng.ngons.at(ng.curr).size());
    e_lens.edge_ctr_new_b.resize(ng.ngons.at(ng.curr).size());
    e_lens.vert_sp_new.resize(ng.ngons.at(ng.curr).size());

    e_lens.vt.resize(ng.ngons.at(ng.curr).size());
    e_lens.split_vt.resize(ng.ngons.at(ng.curr).size());

    for(int i = 1; i < e_lens.vert_ctr_new.size(); i++) {
      e_lens.vt.at(i).fromArray(z);
    }
    for(int i = 1; i < e_lens.vert_sp_new.size(); i++) {
      e_lens.split_vt.at(i).fromArray(z);
    }

    e_lens.new_cell0_id.clear();
    e_lens.new_cell0_id.resize(ng.ngons.at(ng.curr).size());
    e_lens.new_cell1_id.clear();
    e_lens.new_cell1_id.resize(ng.ngons.at(ng.curr).size());
  }
}

void vd_edisc::crt_new_vert() {
  if(trial_type == 1) {

    std::cout << "Old vert " << e_lens.vert_ctr << std::endl;

    apf::ModelEntity* vert_ctr_new_em = e_lens.m->toModel(e_lens.vert_ctr);
    e_lens.vert_ctr_new.at(0) = e_lens.m->createVert(vert_ctr_new_em);
    std::cout << "New vert " << e_lens.vert_ctr_new.at(0)
              << " " << e_lens.m->getModelType(vert_ctr_new_em)
              << "c" << e_lens.m->getModelTag(vert_ctr_new_em)
              << std::endl;

    vert_ctr_new_em = e_lens.m->findModelEntity(0, e_lens.new_cell0_id.at(1));
    e_lens.vert_ctr_new.at(1) = e_lens.m->createVert(vert_ctr_new_em);

    std::cout << "New vert " << e_lens.vert_ctr_new.at(1)
              << " " << e_lens.m->getModelType(vert_ctr_new_em)
              << "c" << e_lens.m->getModelTag(vert_ctr_new_em)
              << std::endl;

    e_lens.create_map[0][vert_ctr_act] = e_lens.vert_ctr_new.at(0);
    e_lens.create_map[1][vert_ctr_act] = e_lens.vert_ctr_new.at(1);

    assert(e_lens.new_cell1_id.at(0) > 0);
    apf::ModelEntity* edge_ctr_new_em = e_lens.m->findModelEntity(1, 
                                                  e_lens.new_cell1_id.at(0));

    e_lens.vert_sp_new.at(0) = e_lens.m->createVert(edge_ctr_new_em);

    std::cout << "New vert_sp " << e_lens.vert_sp_new.at(0)
              << " " << e_lens.m->getModelType(edge_ctr_new_em)
              << "c" << e_lens.m->getModelTag(edge_ctr_new_em)
              << std::endl;
  }
  else {
    assert(trial_type == 2);

    std::cout << "Old vert " << e_lens.vert_ctr << std::endl;

    apf::ModelEntity* vert_ctr_new_em = e_lens.m->toModel(e_lens.vert_ctr);
    e_lens.vert_ctr_new.at(0) = e_lens.m->createVert(vert_ctr_new_em);

    std::cout << "New vert " << e_lens.vert_ctr_new.at(0)
              << " " << e_lens.m->getModelType(vert_ctr_new_em)
              << "c" << e_lens.m->getModelTag(vert_ctr_new_em)
              << std::endl;

    for(int i = 1; i < e_lens.vert_ctr_new.size(); i++) {
      vert_ctr_new_em = e_lens.m->findModelEntity(0, 
                                        e_lens.new_cell0_id.at(i));
      e_lens.vert_ctr_new.at(i) = e_lens.m->createVert(vert_ctr_new_em);

      std::cout << "New vert " << e_lens.vert_ctr_new.at(i)
                << " " << e_lens.m->getModelType(vert_ctr_new_em)
                << "c" << e_lens.m->getModelTag(vert_ctr_new_em)
                << std::endl;
    }

    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      assert(e_lens.new_cell1_id.at(i) > 0);
      vert_ctr_new_em = e_lens.m->findModelEntity(1, 
                                        e_lens.new_cell1_id.at(i));
      e_lens.vert_sp_new.at(i) = e_lens.m->createVert(vert_ctr_new_em);

      std::cout << "New vert_sp " << e_lens.vert_sp_new.at(i)
                << " " << e_lens.m->getModelType(vert_ctr_new_em)
                << "c" << e_lens.m->getModelTag(vert_ctr_new_em)
                << std::endl;

    }

    for(int i = 0; i < e_lens.discs.size(); i++) {
      int v1 = e_lens.slices.at(2*i)-1;
      int v2 = e_lens.slices.at(2*i+1)-1;
      e_lens.create_map[2*i][vert_ctr_act] = e_lens.vert_ctr_new.at(v1);
      e_lens.create_map[2*i+1][vert_ctr_act] = e_lens.vert_ctr_new.at(v2);
    }
  }
}

void vd_edisc::crt_new_cell() {
  if(trial_type == 1) {
    e_lens.new_cell1_id.at(0) = e_lens.c_base_curr->use_free_gmi(1);
    assert(e_lens.new_cell1_id.at(0) > 0);

    e_lens.new_cell0_id.at(0) = e_lens.m->getModelTag(
                                       e_lens.m->toModel(e_lens.vert_ctr));
    e_lens.new_cell0_id.at(1) = e_lens.c_base_curr->use_free_gmi(0);

    assert(e_lens.new_cell0_id.at(0) > 0);
    assert(e_lens.new_cell0_id.at(1) > 0);

    //apf::ModelEntity* edge_ctr_new_em = e_lens.m->findModelEntity(1, 
    //                                              e_lens.new_cell1_id.at(0));

    //e_lens.vert_sp_new.at(0) = e_lens.m->createVert(edge_ctr_new_em);

  }
  else {
    assert(trial_type == 2);

    e_lens.new_cell0_id.at(0) = m_act->getModelTag(
                            e_lens.m->toModel(e_lens.vert_ctr));

    for(int i = 1; i < e_lens.vert_ctr_new.size(); i++) {
      e_lens.new_cell0_id.at(i) = e_lens.c_base_curr->use_free_gmi(0);
      assert(e_lens.new_cell0_id.at(i) > 0);
    }

    assert(e_lens.new_cell0_id.at(0) > 0);

    for(int i = 0; i < e_lens.discs.size(); i++) {
      e_lens.new_cell1_id.at(i) = e_lens.c_base_curr->use_free_gmi(1);
      assert(e_lens.new_cell1_id.at(i) > 0);
    }
    e_lens.new_cell2_id = e_lens.c_base_curr->use_free_gmi(2);
    assert(e_lens.new_cell2_id > 0);
  }
  for(int i = 0; i < e_lens.new_cell0_id.size(); i++) {
    std::cout << e_lens.new_cell0_id.at(i) << " ";
  }
  std::cout << std::endl;
}

void vd_edisc::asgn_new_vert_ext() {
  if(trial_type == 1) {
    ext_new[e_lens.vert_ctr_new.at(0)] = (ext_slice[0] or ext_pc[0]);
    ext_new[e_lens.vert_ctr_new.at(1)] = (ext_slice[1] or ext_pc[0]);
    ext_new[e_lens.vert_sp_new.at(0)] = ext_pc[0];
  }
  else {
    assert(trial_type == 2);

    std::vector<int > paths(0);
    paths.resize(e_lens.slices.size());
    for(int i = 0; i < e_lens.slices.size()/2; i++) {
      int path1_id = e_lens.slices.at(2*i)-1;
      int path2_id = e_lens.slices.at(2*i+1)-1;

      if(paths.at(2*path1_id) == 0)
        paths.at(2*path1_id) = i+1;
      else
        paths.at(2*path1_id+1) = i+1;

      if(paths.at(2*path2_id) == 0)
        paths.at(2*path2_id) = i+1;
      else
        paths.at(2*path2_id+1) = i+1;

      //paths.at(e_lens.slices.at(2*i)-1) = i+1;
      //paths.at(e_lens.slices.at(2*i+1)-1) = i+1;
    }
/*
    int p1 = paths.at(0)-1;
    int p2 = paths.at(1)-1;

    ext_new[e_lens.vert_ctr_new.at(0)] = (ext_pc[p1] or ext_pc[p2] 
                                                          or ext_slice[0]);
*/
    for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      int p1 = paths.at(2*i)-1;
      int p2 = paths.at(2*i+1)-1;

      ext_new[e_lens.vert_ctr_new.at(i)] = (ext_pc[p1] or ext_pc[p2] 
                                                          or ext_slice[i]);
    }

    for(int i = 0; i < e_lens.discs.size(); i++) {
      int v1 = e_lens.slices.at(2*i)-1;
      int v2 = e_lens.slices.at(2*i+1)-1;
      ext_new[e_lens.vert_sp_new.at(i)] = ext_pc[i];
    }
  }
}

// Get the entity list around the new vertices.
void vd_edisc::get_ent_set_new() {
  e_lens.es_vert = e_lens.vert_ctr_new;
  vd_set_up(e_lens.m, &e_lens.es_vert, &e_lens.es_edge);
  vd_set_up(e_lens.m, &e_lens.es_edge, &e_lens.es_surf);
  vd_set_up(e_lens.m, &e_lens.es_surf, &e_lens.es_elem);
  vd_set_down(e_lens.m, &e_lens.es_elem, &e_lens.es_surf);
  vd_set_down(e_lens.m, &e_lens.es_surf, &e_lens.es_edge);
  vd_set_down(e_lens.m, &e_lens.es_edge, &e_lens.es_vert);
}

bool vd_edisc::burn_slice(int slice) {
  // The containers for downward lists used in checking for center vertex 
  // adjacency and obtaining the triangle adjacencies of the active element.
  apf::Downward down_b;
  apf::Downward v_b;

  apf::Up up_b;

  // Check if this slice contains the exterior. If it does, burn the unburned
  // tets adjacent to the current disc in this color.
  bool hit_ext = false;
  //std::cout << "Burning slice: " << slice+1 << std::endl;

  while(elem_ofire.size() > 0) {
    //std::cout << "Elem " << elem_ofire.front() << std::endl;
    e_lens.m->getDownward(elem_ofire.front(), 2, down_b);

    //std::cout << "\t Tri: " << std::endl;

    bool all_burned = true;
    for(int i = 0; i < 4; i++) {
      //std::cout << "\t\t" << down_b[i] << " ";

      if(!tri_burn[down_b[i]]) {
        all_burned = false;

        e_lens.m->getDownward(down_b[i], 0, v_b);
        int e1 = findIn(v_b, 3, vert_ctr_act);

        //std::cout << "v_ctr: " << e1 << " not burned, ";

        e_lens.m->getUp(down_b[i], up_b);

        // If this is a boundary triangle, skip.
        if (e1 < -1) {
        //  std::cout << " v_ctr missing." << std::endl;
        }

        else if (up_b.n == 1) {
          hit_ext = true;
        //  std::cout << "boundary." << std::endl;
        }
        // Else, try the adjacent tet.
        else {
          assert(up_b.n == 2);
          bool skip = false;

          apf::MeshEntity* elem_next;
          // The first is burned.
          if(elem_burn[up_b.e[0]]) {
            // The second is burned.
            if(elem_burn[up_b.e[1]]) {
              // skip
              skip = true;
              //std::cout << "Both elem burned." << std::endl;
            }
            // The second is not burned.
            else {
              assert(up_b.e[1] != elem_ofire.front());
              elem_next = up_b.e[1];
              //std::cout << "Elem next: " << elem_next;
            }
          }
          // The first is not burned.
          else {
            assert(up_b.e[0] != elem_ofire.front());
            elem_next = up_b.e[0];
            //std::cout << "Elem next: " << elem_next;
          }

          if (!skip) {
            e_lens.m->getDownward(elem_next, 0, v_b);

            int e1 = findIn(v_b, 4, vert_ctr_act);
            // Next element belongs to the slice. Set it to burning.
            if (e1 > -1) {
              //std::cout << " burning... " << std::endl;
              elem_slice.at(slice).push_back(elem_next);
              elem_ofire.push_back(elem_next);
              elem_burn[elem_next] = true;
            }
            else {
              //std::cout << " no v_ctr... " << std::endl;
              elem_burn[elem_next] = true;
            }
          }
          // Triangle burned.
        }
        tri_burn[down_b[i]] = true;
      }
      // All neighbours of the element is tried, so remove it.
      //else {
      //  std::cout << " burned " << std::endl;
      //}
    }
    elem_ofire.pop_front();
  }

  //std::cout << "Slice " << slice+1 << " " << "nbr of ent:" 
  //          << elem_slice.at(slice).size() << std::endl;
  //for(int i = 0; i < elem_slice.at(slice).size(); i++) {
  //  std::cout << elem_slice.at(slice][i) << " ["
  //            << findIn(e_lens.es_elem.e, e_lens.es_elem.n, 
  //                                    elem_slice.at(slice][i)) << "] ";
  //}
  //std::cout << std::endl;
  return hit_ext;
}

// The locations of the discs and fins are arbitrary the resulting slices
// can contain arbitrarily shaped triangles belonging to path and circuit 
// 2cells. For consistency, do not consider these triangles in velocity 
// calculation. 
void vd_edisc::burn_slice_wings(int disc_curr) {
  if(trial_type == 3)
    return;
  for (int i = 0; i < e_lens.discs.at(disc_curr).edge.size(); i++) {
    apf::MeshEntity* e_curr = e_lens.discs.at(disc_curr).edge.at(i);
    apf::ModelEntity* mdl = e_lens.m->toModel(e_curr);

    if(e_lens.m->getModelType(mdl) == 2) {
      std::deque<apf::MeshEntity*> edge_fire;
      std::deque<apf::MeshEntity*> surf_fire;

      edge_fire.push_back(e_curr);
      wing_map[e_curr] = true;

      while(edge_fire.size() > 0) {
        apf::Up up;
        apf::Downward down;
        apf::Downward d_vert;

        while(edge_fire.size() > 0) {
          e_lens.m->getUp(edge_fire.front(), up);

          for(int j = 0; j < up.n; j++) {
            if(e_lens.m->toModel(up.e[j]) == mdl and !wing_map[up.e[j]]) {
              surf_fire.push_back(up.e[j]);
              wing_map[up.e[j]] = true;
            }
          }
          edge_fire.pop_front();
        }

        while(surf_fire.size() > 0) {
          e_lens.m->getDownward(surf_fire.front(), 1, down);
          e_lens.m->getDownward(surf_fire.front(), 0, d_vert);

          int v1 = findIn(d_vert, 3, vert_ctr_act);
          assert(v1 > -1);
          int e1 = lookup_tri_ed[v1][0];
          int e2 = lookup_tri_ed[v1][1];

          if(wing_map[ down[e1] ]) {
            assert(!wing_map[ down[e2] ]);
            if(e_lens.m->toModel(down[e2]) == mdl) {
              edge_fire.push_back(down[e2]);
              wing_map[down[e2]] = true;
            }
          }
          else {
            if(e_lens.m->toModel(down[e1]) == mdl) {
              edge_fire.push_back(down[e1]);
              wing_map[down[e1]] = true;
            }
          }
          surf_fire.pop_front();
        }
      }
    }
  }
}

// Burn trial;
void vd_edisc::burn_trial() {

  wing_map.clear();

  for(int i = 0; i < elem_slice.size(); i++)
    elem_slice.at(i).clear();
  elem_slice.clear();

  elem_burn.clear();

  elem_ofire.clear();
  //elem_ofire.reserve(es_elem.n);

  assert(trial_type >= 1 or trial_type <= 3);

  if(trial_type == 1) {
    //std::cout << "1cell insertion, burning..." << std::endl;
    elem_slice.resize(2);
    elem_slice[0].reserve(e_lens.es_elem.size());
    elem_slice[1].reserve(e_lens.es_elem.size());
  }

  else if(trial_type == 2) {
    //std::cout << "2cell insertion, burning..." << std::endl;
    elem_slice.resize(ng.ngons.at(ng.curr).size());
    for(int i = 0; i < elem_slice.size(); i++)
      elem_slice.at(i).reserve(e_lens.es_elem.size());
  }
  else {
    //std::cout << "3cell preconditioning, burning..." << std::endl;
    elem_slice.resize(2);
    elem_slice[0].reserve(e_lens.es_elem.size());
    elem_slice[1].reserve(e_lens.es_elem.size());
  }

  //std::cout << "nbr of elem: " << e_lens.es_elem.n << std::endl;
  //for (int i = 0; i < e_lens.es_elem.n; i++) {
  //  std::cout << e_lens.es_elem.e[i] << " ";
  //}
  //std::cout << std::endl;

  int slice_curr = 0;

  // Assume that exterior happens once. So one of the slices is complete.
  // The other will have at least one unburned triangle with a single upper 
  // adjacency.
  for (int i = 0; i < e_lens.discs.size(); i++) {
    // Find the first triangle bounded by two tetrahedra. Try to burn the slices
    // using the two tetrahedra. Associate the first unassociated disc slice 
    // with the current slice. If one of them is burned already, associate
    // the disc with that slice.
    // Looking at the discs, get the disc slice colors. The top color is the 
    // path color for each disc. The bottom color may be obtained during 
    // burning, by comparing disc colors (slices) to current slice color. If 
    // the top slice color is different color the bottom one. Assert that it  
    // is not colored beforehand.
    bool found = false;
    apf::Up up;
    for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {
      e_lens.discs.at(i).m->getUp(e_lens.discs.at(i).tri.at(j), up);

      if(up.n == 2) {
        j = e_lens.discs.at(i).tri.size();
        found = true;
      }
    }

    //std::cout << "Disc " << i+1 << std::endl; 
    //std::cout << "up[0] " << up.e[0] << " " << elem_burn[up.e[0]] 
    //          << " " << e_lens.slice_map[up.e[0]] 
    //          << std::endl;
    //std::cout << "up[1] " << up.e[1] << " " << elem_burn[up.e[1]]
    //          << " " << e_lens.slice_map[up.e[1]] 
    //          << std::endl;

    if(trial_type != 3)
      assert(found);

    if (trial_type == 3) {
      if(!found) {
        e_lens.slices.at(2*i) = 2;

        elem_slice.at(1).push_back(up.e[0]);
        elem_ofire.push_back(up.e[0]);
        elem_burn[up.e[0]] = true;

        burn_slice(1);
        color_slice(1);

        assert(elem_ofire.empty());

        e_lens.slices.at(2*i+1) = 1;
        elem_slice.at(0).resize(0);
        //elem_slice.at(slice_curr).push_back(up.e[0]);
        //elem_ofire.push_back(up.e[0]);
        //elem_burn[up.e[0]] = true;

        burn_slice(0);
        color_slice(0);

        assert(elem_ofire.empty());
      }
/*
shattered slices should be collectible, there can be more than one piece of 
the slice. Also, check for 1 cell and 2cell insertion. Revise the part with
finding, as well
*/
      else {
        e_lens.slices.at(2*i+1) = 2;

        if(e_lens.m->toModel(up.e[0]) == vd_3c->get_mdl(3, trial_curr)) {
        //if(e_lens.m->toModel(up.e[0]) == 
        //      e_lens.m->findModelEntity(3, trial_curr)) {

          elem_slice.at(1).push_back(up.e[0]);
          elem_ofire.push_back(up.e[0]);
          elem_burn[up.e[0]] = true;
        }
        else {
          elem_slice.at(1).push_back(up.e[1]);
          elem_ofire.push_back(up.e[1]);
          elem_burn[up.e[1]] = true;
        }
        burn_slice(1);
        color_slice(1);
        assert(elem_ofire.empty());

        slice_curr = 0;
        e_lens.slices.at(2*i) = 1;
        for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {
          e_lens.discs.at(i).m->getUp(e_lens.discs.at(i).tri.at(j), up);

          for (int k = 0; k < up.n; k++) {
            if(e_lens.m->toModel(up.e[k]) == vd_3c->get_mdl(3, trial_curr)) {
            //if(e_lens.m->toModel(up.e[k]) == 
            //      e_lens.m->findModelEntity(3, trial_curr)) {
              assert(elem_burn[up.e[k]]);
            }
            else {
              if(!elem_burn[up.e[k]] ) {
                elem_slice.at(slice_curr).push_back(up.e[k]);
                elem_ofire.push_back(up.e[k]);
                elem_burn[up.e[k]] = true;
                burn_slice(slice_curr);
                color_slice(slice_curr);
                assert(elem_ofire.empty());
              }
            }
          }
        }
      }

    }

    else {
      // Assume exterior is not multiply connected. If exterior is encountered,
      // burn all remaining tets using the exterior slice.
      bool hit_ext1 = false;
      bool hit_ext2 = false;
      int slice_ext = -1;
      if(!elem_burn[up.e[0]]) {
        assert(slice_curr < elem_slice.size());
        //std::cout << "Slice " << slice_curr+1 << std::endl;

        if(e_lens.slices.at(2*i) == 0) {
          e_lens.slices.at(2*i) = slice_curr+1;
        }
        else {
          assert(e_lens.slices.at(2*i+1) == 0);
          e_lens.slices.at(2*i+1) = slice_curr+1;
        }

        elem_slice.at(slice_curr).push_back(up.e[0]);

        elem_ofire.push_back(up.e[0]);
        elem_burn[up.e[0]] = true;

        hit_ext1 = burn_slice(slice_curr);
        if(hit_ext1) {
          slice_ext = slice_curr;
          ext_slice[slice_curr] = true;
        }
        color_slice(slice_curr);

        assert(elem_ofire.empty());
        slice_curr++;
      }
      else {
        assert(e_lens.slice_map[up.e[0]] > 0);
        e_lens.slices.at(2*i) = e_lens.slice_map[up.e[0]];
      }

      if (!elem_burn[up.e[1]]) {
        assert(slice_curr < elem_slice.size());
        //std::cout << "Slice " << slice_curr+1 << std::endl;

        if(e_lens.slices.at(2*i) == 0) {
          e_lens.slices.at(2*i) = slice_curr+1;
        }
        else {
          assert(e_lens.slices.at(2*i+1) == 0);
          e_lens.slices.at(2*i+1) = slice_curr+1;
        }

        elem_slice.at(slice_curr).push_back(up.e[1]);

        elem_ofire.push_back(up.e[1]);
        elem_burn[up.e[1]] = true;

        hit_ext2 = burn_slice(slice_curr);
        if(hit_ext2) {
          //assert(!hit_ext1);
          slice_ext = slice_curr;
          ext_slice[slice_curr] = true;
        }
        color_slice(slice_curr);

        elem_ofire.clear();
        slice_curr++;
      }

      else {
        assert(e_lens.slice_map[up.e[1]] > 0);
        e_lens.slices.at(2*i+1) = e_lens.slice_map[up.e[1]];
      }
      if(hit_ext1 or hit_ext2) {
        assert(slice_ext > -1);
        for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {
          e_lens.discs.at(i).m->getUp(e_lens.discs.at(i).tri.at(j), up);

          assert(up.n > 0);
          if(!elem_burn[up.e[0]]) {
            if(up.n == 2)
              assert(elem_burn[up.e[1]]);
            elem_slice.at(slice_ext).push_back(up.e[0]);

            elem_ofire.push_back(up.e[0]);
            elem_burn[up.e[0]] = true;

            burn_slice(slice_ext);
            color_slice(slice_ext);
          }
          else if(up.n == 2 and !elem_burn[up.e[1]]) {
            assert(elem_burn[up.e[0]]);
            elem_slice.at(slice_ext).push_back(up.e[1]);

            elem_ofire.push_back(up.e[1]);
            elem_burn[up.e[1]] = true;

            burn_slice(slice_ext);
            color_slice(slice_ext);
          }
        }

      }
    }
    burn_slice_wings(i);
  }


  //for (int i = 0; i < e_lens.discs.size(); i++) {
  //  std::cout << "disc " << i 
  //            << " slices[0] " << e_lens.slices[2*i] 
  //            << " slices[1] " << e_lens.slices[2*i+1] << std::endl;
  //}

  // Provide VTK output for the slices.
  //for (int i = 0; i < elem_slice.size(); i++) {
    // Write a vtk for the burning entities after burn_slice.
  //  std::vector<std::vector<apf::MeshEntity*> > es_burn 
  //                                (4, std::vector<apf::MeshEntity*> (0) );

  //  es_burn[3].resize(elem_slice.at(i).size());
  //  for (int j = 0; j < elem_slice.at(i).size(); j++) {
  //    es_burn[3].at(j) = elem_slice.at(i).at(j);
  //  }
  //  vd_set_down(e_lens.m, &es_burn.at(3), &es_burn.at(2));
  //  vd_set_down(e_lens.m, &es_burn.at(2), &es_burn.at(1));
  //  vd_set_down(e_lens.m, &es_burn.at(1), &es_burn.at(0));

  //  char s[50];
  //  sprintf(s, "./output/after_burn%d", i+1);
  //  vd_save_vtk_set(e_lens.m, &es_burn, s);
  //}

}

// Collect triangles of the slice, belonging to a 2cell. Also, if the bounding
// disc of the slice is not separated during the lens expansion, they also 
// belong to the slice. 
// Skip the triangles belonging to cells on paths or circuits.
void vd_edisc::collect_slice_tris(
                    std::vector<std::vector<apf::MeshEntity*> > * es_tri) {

  for (int i = 0; i < es_tri->size(); i++)
    es_tri->at(i).clear();
  es_tri->resize(e_lens.vert_ctr_new.size());
  for (int i = 0; i < es_tri->size(); i++)
    es_tri->at(i).reserve(e_lens.es_surf.size());

  for (int i = 0; i < e_lens.es_surf.size(); i++) {
    int slice_curr = e_lens.slice_map[e_lens.es_surf.at(i)];
    if(slice_curr > 0) {
      apf::ModelEntity* mdl = e_lens.m->toModel(e_lens.es_surf.at(i));
      if(e_lens.m->getModelType(mdl) == 2 
          and !wing_map[e_lens.es_surf.at(i)]) {
        es_tri->at(slice_curr - 1).push_back(e_lens.es_surf.at(i));
        //std::cout << e_lens.es_surf.at(i) << " ";
      }
    }
  }
  //std::cout << std::endl;

  for (int i = 0; i < es_tri->size(); i++) {
    //assert(es_tri->at(i).size() > 0);
    assert(e_lens.c1.at(i).size() > 0);
    if(e_lens.c1.at(i).size() == 1) {
      //assert(es_tri->at(i).size() == 0);

      std::cout << "Single 1cell edge on slice." << std::endl;
      bool found = false;
      for (int j = 0; j < e_lens.es_edge.size(); j++) {
        int slice_curr = e_lens.slice_map[e_lens.es_edge.at(j)];
        if(slice_curr - 1 == i) {
          apf::ModelEntity* mdl = e_lens.m->toModel(e_lens.es_edge.at(j));
          if(e_lens.m->getModelType(mdl) == 1) {
            assert(e_lens.c1.at(i).at(0) == e_lens.m->getModelTag(mdl));

            //Removed assert(es_tri->at(i).size() == 0); above, reset the 
            // slice triangles instead...
            es_tri->at(slice_curr - 1).clear();
            apf::Up up;
            e_lens.m->getUp(e_lens.es_edge.at(j), up);
            for (int k = 0; k < up.n; k++) {
              int slice_tri = e_lens.slice_map[up.e[k]];
              assert(slice_tri == slice_curr);
              mdl = e_lens.m->toModel(up.e[k]);
              if(e_lens.m->getModelType(mdl) == 2)
                es_tri->at(slice_curr - 1).push_back(up.e[k]);
            }
            j = e_lens.es_edge.size();
            found = true;
          }
        }
      }
      assert(found);
    }
    std::vector<std::vector<apf::MeshEntity*> > e_set
                   (4, std::vector<apf::MeshEntity*> (0) );

    e_set.at(2) = es_tri->at(i);
    vd_set_up(e_lens.m, &e_set[2], &e_set[3]);
    vd_set_down(e_lens.m, &e_set[3], &e_set[2]);
    vd_set_down(e_lens.m, &e_set[2], &e_set[1]);
    vd_set_down(e_lens.m, &e_set[1], &e_set[0]);

    char s[50];
    sprintf(s, "./output/tri_slice%d", i);

    if(sub_vtk)
      vd_save_vtk_set(e_lens.m, &e_set, s);
  }
  //std::cout << "\tSlice tri " << es_tri->size() << std::endl;

}

// Collect triangles of the slice, belonging to a 2cell. Also, if the bounding
// disc of the slice is not separated during the lens expansion, they also 
// belong to the slice. 
// Skip the triangles belonging to cells on paths or circuits.
void vd_edisc::collect_slice_tri(std::vector<apf::MeshEntity*>* es_tri, 
                                  int slice) {

  es_tri->clear();
  es_tri->reserve(e_lens.es_surf.size());
  for (int i = 0; i < e_lens.es_surf.size(); i++) {
    if(e_lens.slice_map[e_lens.es_surf.at(i)] == slice) {
      apf::ModelEntity* mdl = e_lens.m->toModel(e_lens.es_surf.at(i));
      if(e_lens.m->getModelType(mdl) == 2) {
        es_tri->push_back(e_lens.es_surf.at(i));
        //std::cout << e_lens.es_surf.at(i) << " ";
      }
    }
  }
  //char s[50];
  //sprintf(s, "./output/tri_slice%d", slice);

  //vd_save_vtk_ent(e_lens.m, es_tri, s);
  //std::cout << "\tSlice tri " << es_tri->size() << std::endl;

}

// Collect triangles of the path, belonging to a 2cell. 
void vd_edisc::collect_path_tri(std::vector<apf::MeshEntity*>* es_tri, 
                                  int path_id) {
  assert(trial_type == 2);
  std::vector<apf::MeshEntity*> e_2c(0);
  std::vector<apf::MeshEntity*> s_2c(0);

  std::pair<std::vector<int>, std::vector<int> >* ng_pt =
                              &ng.cells.at(ng.ngons.at(ng.curr).at(path_id));
  e_2c.resize(ng_pt->second.size());

  for(int k = 0; k < ng_pt->second.size(); k++) {
    e_2c.at(k) = c2_edge.at(ng_pt->second.at(k));
  }

  vd_set_up(m_act, &e_2c, &s_2c);

  es_tri->reserve(s_2c.size());
  for (int i = 0; i < s_2c.size(); i++) {
    apf::ModelEntity* mdl = e_lens.m->toModel(s_2c.at(i));
    if(e_lens.m->getModelType(mdl) == 2 ) {
      //std::cout << s_2c.at(i) << " ";
      es_tri->push_back(s_2c.at(i));
    }
  }
  //std::cout << std::endl;

  //std::vector<std::vector<apf::MeshEntity*> > e_set
  //               (4, std::vector<apf::MeshEntity*> (0) );


  //e_set.at(2) = *es_tri;
  //vd_set_up(e_lens.m, &e_set[2], &e_set[3]);
  //vd_set_down(e_lens.m, &e_set[3], &e_set[2]);
  //vd_set_down(e_lens.m, &e_set[2], &e_set[1]);
  //vd_set_down(e_lens.m, &e_set[1], &e_set[0]);

  //char s[50];
  //sprintf(s, "./output/tri_path%d", path_id);
  //vd_save_vtk_set(e_lens.m, &e_set, s);

  //for(int k = 0; k < 4; k++) {
  //  e_set.at(k).clear();
  //}
  //std::cout << "\tSlice tri " << es_tri->size() << std::endl;

}

// Collect the discs of the slices. 
// For an interior 0cell, there should number of ngon sides slides. 
// FIX, each 2cell trial (3cell couples) can be tried for different ngons.
// So at each collect_disc_path, different ngons should be tried.
// 

// TODO Actually, as the mesh is preconditioned, the disc triangles can be 
// selected on the path. Would that be helpful or useful? How to move the 
// vertices?
void vd_edisc::collect_disc_path(int ngon_id) {

  trial_type = 2;
  ng.curr = ngon_id;

  init_elens();
  e_lens.c_base_curr = c_base_act;

  //assert(e_lens.init_flag);
  //assert(e_lens.null_flag);

  collect_disc();
  color_discs();

  //for(int i = 0; i < ng.ngons.at(ng.curr).size(); i++) {
  //  std::cout << "ext_path[" << i << "] " << ext_pc[i] << " ";
  //}
  //std::cout << std::endl;
  //for(int i = 0; i < ng.ngons.at(ng.curr).size(); i++) {
  //  std::cout << "ext_slice[" << i << "] " << ext_slice[i] << " ";
  //}
  //std::cout << std::endl;

  //for (int i = 0; i < e_lens.discs.size(); i++) {
  //  std::cout << "Disc triangle: " << e_lens.discs.at(i).tri.size() << std::endl;
  //}

}

void vd_edisc::precond_map() {
  vd_cd->set_cb_flag(true);
  vd_cd->reload(m_act, c_base_act, vert_ctr_act);

  preid_map = vd_cd->get_vert_c1_map();
}

void vd_edisc::trans_map_pre2act() {
  apf::MeshIterator* it = m_act->begin(1);

  apf::MeshEntity* edge;
  apf::MeshEntity* v_curr;
  apf::ModelEntity* mdl;
  apf::Downward down;

  actid_map.clear();
  std::cout << "Transferring disjoint cell ids " << std::endl;
  while(edge = m_act->iterate(it)) {

    mdl = m_act->toModel(edge);
    int cell_type = m_act->getModelType(mdl);

    if (cell_type == 1) {
      m_act->getDownward(edge, 0, down);
      int v1 = findIn(down, 2, vert_ctr_act);
      if(v1 > -1) {
        v1 = (v1+1) % 2;
        v_curr = down[v1];
        if(m_act == m_trial)
          actid_map[v_curr] = preid_map[trial2pre_map[v_curr]];
        else if(m_act == m_main)
          actid_map[v_curr] = preid_map[main2pre_map[v_curr]];
        //std::cout << v_curr << " " << actid_map[v_curr] << std::endl;
      }
    }
  }
  m_act->end(it);
}

void vd_edisc::expand_disc_path() {

  // TODO simplify. Two checks: Is the current 0c exterior? Is the ext_shell
  // algorithm being used. If both, shell needs to be updated before 
  // find_ext_dir upd_cell changes the topology information of the 0cell.
  bool shell_ext = false;
  bool cell_ext = false;
  if(calc_ext and e_lens.c_base_curr->get_cell_ext_gmi(0, cell_id)) {
    cell_ext = true;
  }

  if(cell_ext) {
    if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
      shell_ext = true;
    }
  }

  double t1;
  // TIME
  double dt_collect = 0;
  double dt_color_disc = 0;
  double dt_reload_inv = 0;

  double dt_crt_ent = 0;

  double dt_calc_vel_init = 0;
  double dt_upd_pos_init = 0;

  double dt_relax = 0;
  double dt_vtk_rest = 0;
  double dt_tot = 0;

  double t0 = PCU_Time();
  double t_init = t0;

  ins_flag = true;
  // Create the new vertices:

  std::cout << "m_act " << m_act << std::endl;

  crt_new_vert_list();

  // Collect the affected cells.
  collect_cell();
  coll_proj_map();

  // TIME
  t1 = PCU_Time();
  dt_collect = t1 - t0;

  calc_ctr_vel(false);
  //if(calc_ext and e_lens.c_base_curr->get_cell_ext_gmi(0, cell_id))
    fix_vel_new();

  calc_sp_vel(false);

  bool skip_ins = false;
/*
  if(detect_inv()) {
    std::cout << "disc inversion, trying to fix (2cell)" << std::endl;
    vd_disc_cut vd_cut;
    if(!do_the_cuts(&vd_cut)) {
      //skip_ins = true;
      //ins_flag = false;
      std::cout << "ins_flag: disc inverting 2cell, try insertion anyways" << std::endl;
    }
    else
      std::cout << "fixed" << std::endl;

    init_elens();
    e_lens.c_base_curr = c_base_act;

    reload_edges();
    overwrite_pc(&vd_cut);

    collect_disc_path(ng.curr);

    crt_new_vert_list();

    collect_cell();
    coll_proj_map();
  }
*/
  // TIME
  t0 = PCU_Time();
  dt_reload_inv = t0 - t1;

  crt_new_cell();
  crt_new_vert();
  asgn_new_vert_ext();


  apf::MeshEntity* v_edge[3];

  std::cout << "Creating the new path 1cell edges" << std::endl;
  for(int i = 0; i < e_lens.discs.size(); i++) {
    int v1 = e_lens.slices.at(2*i)-1;
    int v2 = e_lens.slices.at(2*i+1)-1;

    v_edge[0] = e_lens.vert_ctr_new.at(v1);
    v_edge[1] = e_lens.vert_sp_new.at(i);
    v_edge[2] = e_lens.vert_ctr_new.at(v2);

    apf::ModelEntity* edge_ctr_new_em = m_act->findModelEntity(1,
                                                 e_lens.new_cell1_id.at(i));
    e_lens.edge_ctr_new_t.at(i) = buildElement(m_act, edge_ctr_new_em, 
                                             apf::Mesh::EDGE, v_edge, 0);
    e_lens.edge_ctr_new_b.at(i) = buildElement(m_act, edge_ctr_new_em, 
                                            apf::Mesh::EDGE, &v_edge[1], 0);

    int mdl_type = m_act->getModelType(edge_ctr_new_em);
    int mdl_tag = m_act->getModelTag(edge_ctr_new_em);

    std::cout << "Edge " << e_lens.edge_ctr_new_t.at(i) 
              << "Edge " << e_lens.edge_ctr_new_b.at(i) 
              << mdl_type << "c" << mdl_tag
              << std::endl; 
  }

  // TIME
  t1 = PCU_Time();
  dt_crt_ent = t1 - t0;

  // Expand the disc:
  expand_lens();

  e_lens.m->acceptChanges();

  e_lens.m->getPoint(vert_ctr_act, 0, e_lens.pos);

  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_old);
  }
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_old);
  }

  // Check the motion of the new 0cell vertices.
  std::cout << "\t" << "dt: " << f_calc.vdparam.dt << std::endl;
  std::cout << "\t" << "pos_ctr = " << e_lens.pos  << ";" << std::endl;

  std::cout << "get_ext() " << f_calc.get_ext() 
            << " ext_0cell " << ext_0cell
            << std::endl;

  std::cout << "pos_ctr_new = " << e_lens.pos << std::endl;

  fill_lens();
  fill_lens_elem();

  e_lens.m->acceptChanges();
  recreate_slices();

  e_lens.m->acceptChanges();

  // Recreate the slice entities:

  fill_void();

  e_lens.m->setPoint(e_lens.v_2c, 0, pos_old);

  e_lens.m->acceptChanges();
  //vd_chk_neg_sgn(e_lens.m);

  // Destroy the old entities:
  destroy_ent();

  f_calc.refresh_mesh();

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    f_calc.vd_att_fields(m_act, e_lens.vert_ctr_new.at(i));
  }
  for(int i = 0; i < e_lens.discs.size(); i++) {
    f_calc.vd_att_fields(m_act, e_lens.vert_sp_new.at(i));
  }

  //e_lens.pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
  //              vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;

  //e_lens.m->setPoint(e_lens.v_2c, 0, e_lens.pos);

  // TIME
  t0 = PCU_Time();
  dt_crt_ent = t0 - t1;

  int trial_temp = trial_type;
  trial_type = 2223;
  /*
  if (m_act == m_trial)
    vtk_trial();
  else if (m_act == m_precond)
    vtk_precond();
  */
  vtk_mesh();
  trial_type = trial_temp;

  // Get the entity list.
  get_ent_set_new();

  std::cout << "get_ext() " << f_calc.get_ext() 
            << " ext_0cell " << ext_0cell
            << std::endl;

  std::cout << "pos_ctr_new = " << e_lens.pos << std::endl;

  e_lens.m->acceptChanges();

  upd_cell();

  if(skip_ins) {
    ins_flag = false;
    std::cout << "ins_flag: skip_ins" << std::endl;
    for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_old);
    }
    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_old);
    }
    e_lens.m->setPoint(e_lens.v_2c, 0, pos_old);
  }

  //if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL)
  //  update_shell_gmi();

  //if(shell_ext) {
  //  if(ins_flag)
  //    update_shell_pos();
  //}

  find_ctr_vel();
  fix_vel_new();

  find_pos_non_inv_init();
  // TIME

  t1 = PCU_Time();
  dt_calc_vel_init = t1 - t0;

  trial_temp = trial_type;
  trial_type = 223;
  /*
  if (m_act == m_trial)
    vtk_trial();
  else if (m_act == m_precond)
    vtk_precond();
  */
  vtk_mesh();
  trial_type = trial_temp;

  if(vd_chk_neg_sgn(e_lens.m) > 0) {
    ins_flag = false;
    skip_ins = true;
    std::cout << "ins_flag: element inversion 2cell" << std::endl;
  }

  // TIME
  t0 = PCU_Time();
  dt_upd_pos_init = t0 - t1;

  std::cout << "pow_dis(" << w1.size() + ng.curr + 1 
            << ") = " << calc_energy_diss_rate_sing()
            << ";" << std::endl;
  if(ins_flag) {
    find_ext_dir();
  }

  // TIME
  t1 = PCU_Time();
  dt_relax = t1 - t0;

  //std::cout << m_act << " " << m_main << " " << m_precond << " "<< m_trial
  //          << std::endl;
  trial_temp = trial_type;
  trial_type = 22;

  vtk_mesh();

  trial_type = trial_temp;

  std::cout << "Check for negative entities. " << std::endl;
  // Check for negative elements.
  if(!ins_flag or vd_chk_neg_sgn(e_lens.m) > 0) {
    for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_old);
    }
    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_old);
    }
    ins_flag = false;

    if(e_lens.m == m_trial and !skip_en)
      w2.at(trial_curr).at(ng.curr) = 0.1;
    std::cout << "ins_flag: negative element 2cell" << std::endl;
  }
  else
    ins_flag = true;

  e_lens.m->acceptChanges();

  //assert(vd_chk_neg_sgn(e_lens.m) < 1);
  if(ins_flag) {

    f_calc.vd_upd_vel_field(e_lens.m, e_lens.v_2c);
    for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr_new.at(i), &f_calc);
    }

    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_sp_new.at(i), &f_calc);
    }
  }
  // The motion of vertices

  //vtk_trial();
  vtk_mesh();

  // TIME
  double t_final = PCU_Time();
  dt_vtk_rest = t_final - t1;

  dt_tot = t_final - t_init;
  std::cout << " dt_collect " << dt_collect
            << " dt_color_disc " << dt_color_disc
            << " dt_reload_inv " << dt_reload_inv
            << " dt_crt_ent " << dt_crt_ent
            << " dt_calc_vel_init " << dt_calc_vel_init
            << " dt_upd_pos_init " << dt_upd_pos_init
            << " dt_relax " << dt_relax
            << " dt_vtk_rest " << dt_vtk_rest
            << " dt_tot " << dt_tot
            << std::endl;

}


// Color the slice map. Disc entities are colored -i-1, where i is the index of 
// the disc. The slice entities are colored i+1, where i is the slice number.
// TODO this should change to incorporate 2cell insertion, as well.
void vd_edisc::color_disc(int i) {
  //std::cout << "Coloring the disc " << i << std::endl;

  e_lens.null_flag = false;
  //std::cout << "Coloring the the disc entities." << std::endl;
  //std::cout << "Tris " << std::endl;
  for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {
    e_lens.slice_map[e_lens.discs.at(i).tri.at(j)] = -i-1;
    //std::cout << e_lens.discs.at(i).tri.at(j) << " " << -i-1  << std::endl;
  }
  //std::cout << "Edges " << std::endl;
  for (int j = 0; j < e_lens.discs.at(i).edge.size(); j++) {
    e_lens.slice_map[e_lens.discs.at(i).edge.at(j)] = -i-1;
    //std::cout << e_lens.discs.at(i).edge.at(j) << " " << -i-1  << std::endl;
  }

  for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {
    tri_burn[e_lens.discs.at(i).tri.at(j)] = true;
  }

  //std::cout << "Coloring the slices..." << std::endl;

}

// Color the slice map. Disc entities are colored -i-1, where i is the index of 
// the disc. The slice entities are colored i+1, where i is the slice number.
// TODO this should change to incorporate 2cell insertion, as well.
void vd_edisc::color_slice(int i) {
  
  //std::cout << "vert_ctr_act " << vert_ctr_act 
  //          << " e_lens.vert_ctr " << e_lens.vert_ctr << std::endl; 
  e_lens.null_flag = false;

  for (int j = 0; j < elem_slice.at(i).size(); j++) {

    e_lens.slice_map[elem_slice.at(i).at(j)] = i+1;
    //std::cout << "Elem " << elem_slice.at(i).at(j) << ", slice" << i+1 
    //          << std::endl;

    apf::Downward down;
    apf::Downward v;
    // Triangles
    e_lens.m->getDownward(elem_slice.at(i).at(j), 2, down);
    for (int k = 0; k < 4; k++) {
      //std::cout << "Tri " << down[k] << ", slice ";
      if(e_lens.slice_map[down[k]] == 0) {
        e_lens.m->getDownward(down[k], 0, v);
        if(findIn(v,3,e_lens.vert_ctr) > -1) {
          e_lens.slice_map[down[k]] = i+1;
        }
        //std::cout << "F " << e_lens.slice_map[down[k]] << ", ";
      }
      else {
        //std::cout << "NF " << e_lens.slice_map[down[k]] << ", ";
      }
    }
    //std::cout << std::endl;

    // Edges
    e_lens.m->getDownward(elem_slice.at(i).at(j), 1, down);
    for (int k = 0; k < 6; k++) {
      //std::cout << "Edge " << down[k] << std::endl;
      if(e_lens.slice_map[down[k]] == 0) {
        e_lens.m->getDownward(down[k], 0, v);

        //vd_print_vert(e_lens.m, down[k]);
        if(findIn(v,2,vert_ctr_act) > -1) {
          e_lens.slice_map[down[k]] = i+1;
        }
        //std::cout << "\t slice ";
        //std::cout << "F " << e_lens.slice_map[down[k]] << ", ";
      }
      else {
        //std::cout << "NF " << e_lens.slice_map[down[k]] << ", ";
      }
    }
    //std::cout << std::endl;

  }

  //std::cout << "Elems: " << std::endl;
  //for(int i = 0; i < e_lens.es_elem.size(); i++) {
    //std::cout << e_lens.es_elem.at(i) << ", " 
    //              << e_lens.slice_map[e_lens.es_elem.at(i)] << std::endl;
  //}

  //std::cout << "Tris: " << std::endl;
  //for(int i = 0; i < e_lens.es_surf.size(); i++) {
  //  std::cout << e_lens.es_surf.at(i) << ", " 
  //                << e_lens.slice_map[e_lens.es_surf.at(i)] << std::endl;
  //}

  //std::cout << "Edges: " << std::endl;
  //for(int i = 0; i < e_lens.es_edge.size(); i++) {
  //  std::cout << e_lens.es_edge.at(i) << ", " 
  //                << e_lens.slice_map[e_lens.es_edge.at(i)] << std::endl;
  //}

}

// ---------------------------------------------------
// Try 1cell and 2cell insertions. If one minimum energy configuration is found
// apply to the actual mesh.
// <cell_dim, cell_id>, cell_dim = -1 if not successful 
std::pair<int, int> vd_edisc::try_insert() {
  if(!pre_flag) {
    return std::make_pair(-1,-1);
  }

  en_min_1c = 0.1;
  en_min_2c = 0.1;

  // TIME
  double t0 = PCU_Time();
  double t1;
  double t_1cell;
  double t_2cell;
  if(wg_tag == WG_TYPE::TRI or wg_tag == WG_TYPE::EDGE)
    calc_max_diss_wg();
  if(wg_tag == WG_TYPE::SPHERE)
    calc_max_diss_trial_wg();
  else {
    try_1cell();
    t1 = PCU_Time();
    t_1cell = t1 - t0;
    try_2cell();
    t0 = PCU_Time();
    t_2cell = t0 - t1;
  }
  
  std::cout << "Cell_id " << cell_id 
//            << "Trials took: 1c: " << t_1cell
//            << " 2c: " << t_2cell
            << " energy_1c: " << en_min_1c
            << " energy_2c: " << en_min_2c
            << " energy_th: " << energy_th
            << std::endl;

  // en_min should be <= 0.
  if(pre_flag and en_min_1c < -std::numeric_limits<double>::epsilon()
   and en_min_1c < energy_th - std::numeric_limits<double>::epsilon() and
      en_min_1c - en_min_2c < -std::numeric_limits<double>::epsilon()) {
    std::cout << "Minimum energy 1cell insertion " << circ_min 
              << " en_min " << en_min_1c << std::endl;

    trial_type = 1;

    act_main();
    precond_mesh();

    collect_pc();

    // This updates the cell complex, as well.
    insert_1cell(circ_min);
// TODO delegate this process to vd_sim. Instead, return the ids of the 0cells
// to be modified. If the cell number changes, the cell_base structure and the
// cell_ins_chk should be reinitialized.
// Keep an update list for 1,2,3-cells. Make the list accessible to the vd_sim.
// After each mesh is modified, collect the lists. Update the topology at once.
    if(ins_flag) {
      std::cout << "Successful insertion " << "1c" 
                << e_lens.new_cell1_id.at(0) << std::endl;

      evolve_ins();

      for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        vd_cell_val val_vert(e_lens.m, e_lens.c_base_curr, 
                                  e_lens.vert_ctr_new.at(i));
        assert(val_vert.vd_vert_valid());
      }

      for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
        vd_cell_val val_vert(e_lens.m, e_lens.c_base_curr, 
                                  e_lens.vert_sp_new.at(i));
        assert(val_vert.vd_vert_valid());
      }

      //chk_ma_new();
      return std::make_pair(1, e_lens.new_cell1_id.at(0));
    }
    else {

      std::cout << "Unsuccessful insertion " << "1c" 
                << e_lens.new_cell1_id.at(0) << std::endl;
      find_ctr_vel();
      //if(calc_ext and e_lens.c_base_curr->get_cell_ext_gmi(0, cell_id))
        fix_vel_new();

      find_max_disp();
      move_ctr();
      e_lens.pos = vd_get_center(e_lens.m, &e_lens.vert_ctr_new);
      e_lens.m->setPoint(e_lens.vert_sp_new.at(0), 0, e_lens.pos);

      e_list->refresh();
      vd_glens g_lens2(m_main, e_lens.c_base_curr, e_list);
      g_lens2.set_field_calc(&f_calc, false);

      g_lens2.set_verify(false);
      g_lens2.set_inv(false);
      g_lens2.set_precond(false);

      std::pair<int, int> new_0c;

      new_0c = g_lens2.col_cell(1, e_lens.new_cell1_id.at(0), cell_id);
      assert(new_0c.first == 0 and new_0c.second == cell_id);

      if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
        ext_shell* e_sh = f_calc.get_e_sh();

        //if(ext_0cell) {
        e_sh->set_shell(0, new_0c.second - 1, sh_old, sh_ext);
      }

      //c_ins->update_path_gmi(new_0c.second);
      //if(new_0c.second != e_lens.new_cell0_id.at(0))
      //  c_ins->update_path_gmi(e_lens.new_cell0_id.at(0));

      assert(vd_chk_neg_sgn(e_lens.m) == 0);

      //chk_ma_new();
      return std::make_pair(-1,-1);
    }
  }
  else if (pre_flag and en_min_2c < -std::numeric_limits<double>::epsilon()
   and en_min_2c < energy_th - std::numeric_limits<double>::epsilon()) {
    std::cout << "Minimum energy 2cell insertion " << path_min << " "
              << ngon_min << " en_min " << en_min_2c << std::endl;
    trial_type = 2;

    // This updates the cell complex, as well.
    act_main();

    //e_list.refresh();
    precond_mesh();

    //e_list.refresh();
    collect_pc();

    insert_2cell(path_min, ngon_min);
    //e_list.refresh();

    if(ins_flag) {
      std::cout << "Successful insertion " << "2c" 
                << e_lens.new_cell2_id << std::endl;

      evolve_ins();

      for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        vd_cell_val val_vert(e_lens.m, e_lens.c_base_curr, 
                                  e_lens.vert_ctr_new.at(i));
        assert(val_vert.vd_vert_valid());
      }

      for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
        vd_cell_val val_vert(e_lens.m, e_lens.c_base_curr, 
                                  e_lens.vert_sp_new.at(i));
        assert(val_vert.vd_vert_valid());
      }

      //chk_ma_new();
      return std::make_pair(2, e_lens.new_cell2_id);
    }
    else {
      std::cout << "Unsuccessful insertion " << "2c" 
                << e_lens.new_cell2_id << std::endl;
      find_ctr_vel();
      //if(calc_ext and e_lens.c_base_curr->get_cell_ext_gmi(0, cell_id))
        fix_vel_new();

      find_max_disp();
      move_ctr();

      find_sp_vel();
      //if(calc_ext and e_lens.c_base_curr->get_cell_ext_gmi(0, cell_id))
        fix_vel_new();

      find_max_disp();
      move_sp();
      e_lens.pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
                    vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;
      //}
      e_lens.m->setPoint(e_lens.v_2c, 0, e_lens.pos);

      e_list->refresh();
      vd_glens g_lens2(m_main, e_lens.c_base_curr, e_list);
      g_lens2.set_field_calc(&f_calc, false);

      g_lens2.set_verify(false);
      g_lens2.set_inv(false);
      g_lens2.set_precond(false);

      std::pair<int, int> new_0c;

      new_0c = g_lens2.col_cell(2, e_lens.new_cell2_id, cell_id);
      assert(new_0c.first == 0 and new_0c.second == cell_id);

      if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
        ext_shell* e_sh = f_calc.get_e_sh();

        //if(ext_0cell) {
        e_sh->set_shell(0, new_0c.second - 1, sh_old, sh_ext);
      }

      //c_ins->update_path_gmi(new_0c.second);
      //if(new_0c.second != e_lens.new_cell0_id.at(0))
      //  c_ins->update_path_gmi(e_lens.new_cell0_id.at(0));

      //chk_ma_new();
      return std::make_pair(-1,-1);
    }
  }
  return std::make_pair(-1,-1);
}

// ---------------------------------------------------
// 1cell insertion
void vd_edisc::try_1cell() {

  trial_type = 1;
  reload_trial();

  //std::cout << "Energy before " << calc_energy() << std::endl;

  std::cout << "Circ_sz: " << circ_sz << std::endl;

  //c_ins->print_circ(cell_id-1);
  // Store the energy change associated with the insertion.
  w1.resize(circ_sz);
  e1.resize(circ_sz);
  w1_exp.resize(circ_sz);
 
  circ_min = 0;
  en_min_1c = 0.1;
  for (int circ = 0; circ < circ_sz; circ++) {
    w1.at(circ) = 0.1;
  }

  for (int circ = 0; circ < circ_sz; circ++) {
    int circ_ty = get_circ_type(circ);
    std::cout << "Circuit " << circ << " type: " << circ_ty << std::endl;
    if (circ_ty < 2) {

      insert_1cell(circ);

      e_lens.c_base_curr->coll_cell_gmi(1, e_lens.new_cell1_id.at(0), cell_id);
      if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
        //f_calc.set_e_sh(&e_sh_save);
        for (int i = 0; i < e_lens.new_cell0_id.size(); i++)
          f_calc.get_e_sh()->set_shell(0, e_lens.new_cell0_id.at(i) - 1, 
                                         sh_old, false);
        f_calc.get_e_sh()->set_shell(1, e_lens.new_cell1_id.at(0) - 1, 
                                         sh_old, false);
        f_calc.get_e_sh()->set_shell(0, cell_id - 1, sh_old, sh_ext);
      }
      reload_trial();
    }
    else {
      w1.at(circ) = 0.1;
    }
/*
    // Is it possible that a collapsed configuration is locally minimal in a 
    // isotropic case? TODO In an isotropic case, the split should continue 
    // until euler number is satisfied.
    if(splt_all) {
      if(circ == 0) {
        en_min = w.at(circ);
        circ_min = circ;
      }
      else if(w.at(circ) < en_min) {
        en_min = w.at(circ);
        circ_min = circ;
      }
    }
    else if(w.at(circ) < en_min) {
*/
    if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
      if(sh_used) {
        energy_th = w1.at(circ);
        en_min_1c = w1.at(circ);
        circ_min = circ;
        std::cout << " sh_used." << std::endl;
      }
      else if(w1.at(circ) < en_min_1c) {
        en_min_1c = w1.at(circ);
        circ_min = circ;
      }
    }
    else if(w1.at(circ) < en_min_1c) {
      en_min_1c = w1.at(circ);
      circ_min = circ;
    }
    e_lens.clear();

  }

  for (int i = 0; i < w1.size(); i++) {
    //std::cout << "Energy of circ " << i << ": " << w1.at(i) << std::endl;
    std::cout << "E(i," << i+1 << ")= " << w1.at(i) << ";" << std::endl;
  }


/*
*/

}

// ---------------------------------------------------
// 2cell insertion
// 
void vd_edisc::try_2cell() {

  // Store the energy change associated with the 2-cell insertion connecting 
  // the 3-cell couple. The 3-cell couple index and the ngon combination are 
  // kept.

  trial_type = 2;
  std::cout << "Trying 2cell insertions " << std::endl;
  std::cout << "Cell_id " << cell_id 
            << " paths: " << pt_sz << std::endl;
  refresh_e();

  for(trial_curr = 0; trial_curr < pt_sz; trial_curr++) {

    std::cout << trial_curr << std::endl;
    // c_base->print_ent();

    // Collect the relevant cells for the current 3cell couple to be joined.
    ng.clear();
    //c_ins->print_path(cell_id-1);
    get_path_ngon(trial_curr, &ng);

    ng.print();

    // Store the energy change associated with the 2-cell insertion connecting 
    // the 3-cell couple. The 3-cell couple index and the ngon combination are 
    // kept.
    w2.at(trial_curr).resize(ng.ngons.size());
    for(int i = 0; i < ng.ngons.size(); i++) {
      w2.at(trial_curr).at(i) = 0.1;
    }
    e2.at(trial_curr).resize(ng.ngons.size());
    for(int i = 0; i < ng.ngons.size(); i++) {
      e2.at(trial_curr).at(i) = 0.;
    }
    int step = 10;
    w2_exp.at(trial_curr).resize(ng.ngons.size());
    for(int i = 0; i < ng.ngons.size(); i++) {
      w2_exp.at(trial_curr).at(i).resize(step);
    }

    // For each ngon to combine the 3-cell couple, try the insertion. Calculate
    // the energy and compare to the least energy insertion.
    for(int i = 0; i < ng.ngons.size(); i++) {
//      if(ng.ngons.at(i).size() != 2) {
      reload_trial();

      ng.clear();
      get_path_ngon(trial_curr, &ng);

      collect_disc_path(i);
      expand_disc_path();
      m_act->acceptChanges();
      //chk_ma_new();

      if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
        if(sh_used) {
          energy_th = w2.at(trial_curr).at(i);
          ngon_min = i;
          path_min = trial_curr;
          std::cout << "path_min " << path_min << " ngon_min " 
                    << ngon_min << " sh_used." << std::endl;
        }
        else if(w2.at(trial_curr).at(i) < en_min_2c) {
          en_min_2c = w2.at(trial_curr).at(i);
          ngon_min = i;
          path_min = trial_curr;
          std::cout << "path_min " << path_min << " ngon_min " 
                    << ngon_min << std::endl;
        }
      }
      else if(w2.at(trial_curr).at(i) < en_min_2c) {
        en_min_2c = w2.at(trial_curr).at(i);
        ngon_min = i;
        path_min = trial_curr;
        std::cout << "path_min " << path_min << " ngon_min " 
                  << ngon_min << std::endl;
      }
      e_lens.clear();

      e_lens.c_base_curr->coll_cell_gmi(2, e_lens.new_cell2_id, cell_id);
      if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
        //f_calc.set_e_sh(&e_sh_save);
        for (int i = 0; i < e_lens.new_cell0_id.size(); i++)
          f_calc.get_e_sh()->set_shell(0, e_lens.new_cell0_id.at(i) - 1, 
                                         sh_old, false);
        for (int i = 0; i < e_lens.new_cell1_id.size(); i++)
          f_calc.get_e_sh()->set_shell(1, e_lens.new_cell1_id.at(i) - 1, 
                                         sh_old, false);
        f_calc.get_e_sh()->set_shell(2, e_lens.new_cell2_id - 1, sh_old, false);

        f_calc.get_e_sh()->set_shell(0, cell_id - 1, sh_old, sh_ext);
      }
//      }
//      else {
//        std::cout << "Skipping digon insertion. " << std::endl;
//      }
    }

    for (int i = 0; i < w2.at(trial_curr).size(); i++) {
      //std::cout << "Energy of ngon " << i << ": " << w2.at(i) << std::endl;
      std::cout << "E(i," << i+1+w1.size() << ")= " << w2.at(trial_curr).at(i) 
                << ";" << std::endl;
    }

  }
  // TODO should be after each 2cell ngon trial. 
}

// 2cell insertion interface.
void vd_edisc::insert_2cell(int c3_cp_id, int ngon_id) {

  trial_curr = c3_cp_id;
  get_path_ngon(c3_cp_id, &ng);
  ng.curr = ngon_id;
  //e_lens.null_slice_map();
  collect_disc_path(ngon_id);

  expand_disc_path();
  e_lens.m->acceptChanges();
  //e_lens.m->verify();

}

// Check both disjoint graphs, if both of them contain a 3cell, return 0, if  
// one of them contains a 3cell, return the index+1, else return -1.
int vd_edisc::cont_3cell_min (s_graph* circ_graph) { 
  bool first = false;
  for (int i = 0; i < circ_graph->first.size(); i++) {
    if (circ_graph->first.at(i) < cells3.conn.size()) {
      first = true;
      i = circ_graph->first.size();
    }
  }

  bool second = false;
  for (int i = 0; i < circ_graph->second.size(); i++) {
    if (circ_graph->second.at(i) < cells3.conn.size()) {
      second = true;
      i = circ_graph->second.size();
    }
  }

  if (first) {
    if (second)
      return 0;
    else
      return 1;
  }
  else if (second)
    return 2;
  else
    return -1;
}

// Insert an edge inside the 3cell. Used in preconditioning step to treat
// disconnected 3cell couples.
void vd_edisc::init_elens() {
  //std::cout << "Trial type " << trial_type << std::endl;

  e_lens.clear();

  e_lens.set_m(m_act);
  e_lens.vert_ctr = vert_ctr_act;

  init_ent_set(&e_lens.es_vert, vert_ctr_act);
  vd_set_up(e_lens.m, &e_lens.es_vert, &e_lens.es_edge);
  vd_set_up(e_lens.m, &e_lens.es_edge, &e_lens.es_surf);
  vd_set_up(e_lens.m, &e_lens.es_surf, &e_lens.es_elem);

  e_lens.null_slice_map();

  if(trial_type == 0) {
    e_lens.init_flag = false;
  }

  else if(trial_type == 1) {
    e_lens.discs.resize(1);
    e_lens.discs.at(0).clear();
    e_lens.discs.at(0).m = m_act;
    e_lens.discs.at(0).v_ctr = vert_ctr_act;

    e_lens.create_map.resize(2);
    e_lens.tri_tet_map.resize(2);
    e_lens.slices.resize(2);
    e_lens.slices.at(0) = 0;
    e_lens.slices.at(1) = 0;

    e_lens.split_map.resize(1);
    e_lens.tri_map.resize(2);

    e_lens.slice_map.clear();
    e_lens.init_flag = true;
  }

  else if(trial_type == 2) {
    //ng.curr creates problems

    //std::cout << "ng.curr " << ng.curr << std::endl;
    e_lens.discs.resize(ng.ngons.at(ng.curr).size());
    e_lens.create_map.resize(2*ng.ngons.at(ng.curr).size());
    e_lens.tri_tet_map.resize(2*ng.ngons.at(ng.curr).size());
    e_lens.slices.resize(2*ng.ngons.at(ng.curr).size());

    e_lens.split_map.resize(ng.ngons.at(ng.curr).size());
    e_lens.tri_map.resize(2*ng.ngons.at(ng.curr).size());

    // Initialize the slices as -1.
    for(int j = 0; j < e_lens.slices.size(); j++) {
      e_lens.slices.at(j) = 0;
    }

    for(int j = 0; j < e_lens.discs.size(); j++) {
      e_lens.discs.at(j).clear();
      e_lens.discs.at(j).m = m_act;
      e_lens.discs.at(j).v_ctr = vert_ctr_act;
    }

    tri_burn.clear();
    e_lens.slice_map.clear();
    e_lens.init_flag = true;

  }
  else if (trial_type == 3) {
    e_lens.discs.resize(1);
    e_lens.discs.at(0).clear();
    e_lens.discs.at(0).m = m_act;
    e_lens.discs.at(0).v_ctr = vert_ctr_act;

    e_lens.create_map.resize(2);
    e_lens.tri_tet_map.resize(2);
    e_lens.slices.resize(2);
    e_lens.slices.at(0) = 0;
    e_lens.slices.at(1) = 0;

    e_lens.split_map.resize(1);
    e_lens.tri_map.resize(2);
    e_lens.init_flag = true;
  }
  sh_used = false;

  //vd_print_ent(m_act);
  cell_flag = false;
  clear_vel_map();
}


// Insert an edge inside the 3cell. Used in preconditioning step to treat
// disconnected 3cell couples.
void vd_edisc::reset_elens_ent() {

  e_lens.clear_old_ent();

  for (int i = 0; i < e_lens.discs.size(); i++)
    e_lens.discs.at(i).clear();
  e_lens.slice_map.clear();

  e_lens.clear_create_maps();

  assert(trial_type == 1 or trial_type == 2);

  init_ent_set(&e_lens.es_vert, vert_ctr_act);
  vd_set_up(e_lens.m, &e_lens.es_vert, &e_lens.es_edge);
  vd_set_up(e_lens.m, &e_lens.es_edge, &e_lens.es_surf);
  vd_set_up(e_lens.m, &e_lens.es_surf, &e_lens.es_elem);

  e_lens.null_slice_map();

  if(trial_type == 0) {
    e_lens.init_flag = false;
  }

  else if(trial_type == 1) {
    e_lens.discs.resize(1);
    e_lens.discs.at(0).clear();
    e_lens.discs.at(0).m = m_act;
    e_lens.discs.at(0).v_ctr = vert_ctr_act;

    e_lens.create_map.resize(2);
    e_lens.tri_tet_map.resize(2);
    e_lens.slices.resize(2);
    e_lens.slices.at(0) = 0;
    e_lens.slices.at(1) = 0;

    e_lens.split_map.resize(1);
    e_lens.tri_map.resize(2);

    e_lens.slice_map.clear();
    e_lens.init_flag = true;
  }

  else if(trial_type == 2) {
    //ng.curr creates problems

    //std::cout << "ng.curr " << ng.curr << std::endl;
    e_lens.discs.resize(ng.ngons.at(ng.curr).size());
    e_lens.create_map.resize(2*ng.ngons.at(ng.curr).size());
    e_lens.tri_tet_map.resize(2*ng.ngons.at(ng.curr).size());
    e_lens.slices.resize(2*ng.ngons.at(ng.curr).size());

    e_lens.split_map.resize(ng.ngons.at(ng.curr).size());
    e_lens.tri_map.resize(2*ng.ngons.at(ng.curr).size());

    // Initialize the slices as -1.
    for(int j = 0; j < e_lens.slices.size(); j++) {
      e_lens.slices.at(j) = 0;
    }

    for(int j = 0; j < e_lens.discs.size(); j++) {
      e_lens.discs.at(j).clear();
      e_lens.discs.at(j).m = m_act;
      e_lens.discs.at(j).v_ctr = vert_ctr_act;
    }

    tri_burn.clear();
    e_lens.slice_map.clear();
    e_lens.init_flag = true;

  }
}

void vd_edisc::collect_disc() {

  assert(e_lens.init_flag);
  assert(e_lens.null_flag);

  if(trial_type == 1) {
    std::cout << "0cell" << cell_id << ", circ:" << trial_curr << std::endl;
    //c_ins->print_circ(cell_id-1);

    circuit* circ_tup;
    circ_tup = get_circ(trial_curr);

    get_disc_edges(circ_tup, &e_lens.discs.at(0));

  }

  else if(trial_type == 2) {
    std::cout << "0cell" << cell_id << ", ngon:" << ng.curr << std::endl;

    std::vector <std::vector<apf::MeshEntity*> > es_tri
                 (0, std::vector<apf::MeshEntity*> (0) );

    es_tri.resize(ng.ngons.at(ng.curr).size());

    vd_entlist_v ent_list(e_lens.m, e_lens.vert_ctr, e_lens.c_base_curr);
    //ent_list.print();

    apf::MeshEntity* e_3;

    int id_curr = ng.cells.at(ng.ngons.at(ng.curr).at(0)).first.at(0);
    //e_lens.c3_f = ng.cells.at(ng.ngons.at(ng.curr).at(0)).first.at(0);
    if (id_curr == -1) {
    //if (e_lens.c3_f == -1) {
      e_lens.c3_f = -1;
      e_lens.e_3c_t_flag = false;

      for(int j = 0; j < ng.ngons.at(ng.curr).size(); j++)
        ext_pc[j] = true;
    }
    else {
      apf::ModelEntity* mdl_curr = vd_cd->get_mdl(3, id_curr);
      e_lens.c3_f = e_lens.m->getModelTag(mdl_curr);
      // Find the innermost edge of the 3cell.
      //e_lens.e_3c_t = find_edge_clst_int(
      //                &ent_list.e.at(3).at(e_lens.c3_f-1).at(1), e_lens.c3_f);

      std::vector<apf::MeshEntity*> tet_3c(0);
      std::vector<apf::MeshEntity*> es_surf(0);
      std::vector<apf::MeshEntity*> es_edge(0);
      tet_3c = vd_cd->get_dim_ent(3, id_curr);

      vd_set_down(e_lens.m, &tet_3c, &es_surf);
      vd_set_down(e_lens.m, &es_surf, &es_edge);

      int edge_nbr = 0;
      for(int j = 0; j < es_edge.size(); j++) {
        apf::ModelEntity* mdl = e_lens.m->toModel(es_edge.at(j));

        if(mdl_curr == mdl) {
          apf::Downward d_v;
          e_lens.m->getDownward(es_edge.at(j), 0, d_v);

          int j_v = findIn(d_v, 2, e_lens.vert_ctr);
          if(j_v > -1) {            
            e_lens.e_3c_t = es_edge.at(j);
            edge_nbr = edge_nbr + 1;
          }
        }
        if(edge_nbr > 1)
          j = es_edge.size();
      }
      assert(edge_nbr == 1);
      //e_lens.e_3c_t = ent_list.e.at(3).at(e_lens.c3_f-1).at(1).at(0);
      e_lens.e_3c_t_flag = true;
    }
    std::cout << "3cf " << e_lens.c3_f;

    id_curr = ng.cells.at(ng.ngons.at(ng.curr).at(0)).first.back();
    if (id_curr == -1) {
    //if (e_lens.c3_e == -1) {
      e_lens.c3_e = -1;
      e_lens.e_3c_b_flag = false;

      for(int j = 0; j < ng.ngons.at(ng.curr).size(); j++)
        ext_pc[j] = true;
    }
    else {
      apf::ModelEntity* mdl_curr = vd_cd->get_mdl(3, id_curr);
      //e_lens.c3_e = ng.cells.at(ng.ngons.at(ng.curr).at(0)).first.back();
      e_lens.c3_e = e_lens.m->getModelTag(mdl_curr);

      std::vector<apf::MeshEntity*> tet_3c(0);
      std::vector<apf::MeshEntity*> es_surf(0);
      std::vector<apf::MeshEntity*> es_edge(0);
      tet_3c = vd_cd->get_dim_ent(3, id_curr);

      vd_set_down(e_lens.m, &tet_3c, &es_surf);
      vd_set_down(e_lens.m, &es_surf, &es_edge);

      int edge_nbr = 0;
      for(int j = 0; j < es_edge.size(); j++) {
        apf::ModelEntity* mdl = e_lens.m->toModel(es_edge.at(j));

        if(mdl_curr == mdl) {
          apf::Downward d_v;
          e_lens.m->getDownward(es_edge.at(j), 0, d_v);

          int j_v = findIn(d_v, 2, e_lens.vert_ctr);
          if(j_v > -1) {            
            e_lens.e_3c_b = es_edge.at(j);
            edge_nbr = edge_nbr + 1;
          }
        }
        if(edge_nbr > 1)
          j = es_edge.size();
      }
      assert(edge_nbr == 1);

      //e_lens.e_3c_b = ent_list.e.at(3).at(e_lens.c3_e-1).at(1).at(0);
      e_lens.e_3c_b_flag = true;
    }

    std::cout << " 3ce " << e_lens.c3_e 
              << std::endl;

  //TODO find the edges across the vertex at the end of the edge e_3,  
  // belonging to the triangles bounding the tets surrounding the e_3.

    for(int j = 0; j < ng.ngons.at(ng.curr).size(); j++) {

      es_tri.at(j).clear();

      // If the initial 3cell is not the exterior, the first edge inside the
      // starting 3cell. This is only necessary once, but this is easier to 
      // follow.
      // int c3_f = ng.cells.at(ng.ngons.at(ng.curr).at(j)).first.at(0);

      if (e_lens.c3_f == -1) {
      }
      else {
        int c2_f = ng.cells.at(ng.ngons.at(ng.curr).at(j)).second.at(0);
        apf::MeshEntity* e_temp;

        e_temp = c2_edge.at(c2_f);

        //std::cout << "Flag: " << e_lens.e_3c_t_flag << " ";
        //std::cout << "First edge, 3c" << e_lens.c3_f << " , 2c" << c2_f 
        //          << std::endl;
        //std::cout << "First edge, " << e_lens.e_3c_t << " , 2c edge " << e_temp 
        //          << std::endl;


        //std::cout << "Edge3 " << e_3 << " , Edge2 " << e_temp << std::endl;

        assert(e_lens.vd_mesh_find_short(e_lens.e_3c_t, e_temp, 
                                               &es_tri.at(j), e_lens.c3_f));
      }

      // Iterate over the cells of the current path. 
      std::pair<std::vector<int>, std::vector<int> >* ng_pt =
                                     &ng.cells.at(ng.ngons.at(ng.curr).at(j));
      for(int k = 1; k < ng_pt->first.size()-1; k++) {
        int c3_curr = ng_pt->first.at(k);

        int c2_fr = ng_pt->second.at(k-1);
        int c2_sc = ng_pt->second.at(k);

        apf::MeshEntity* e_temp;
        apf::MeshEntity* e_temp2;
        //std::cout << "2c" << c2_fr << " 3c" << c3_curr 
        //          << " , 2c" << c2_sc << std::endl;
        //e_temp = ent_list.e.at(2).at(c2_fr-1).at(1).at(0);
        //e_temp2 = ent_list.e.at(2).at(c2_sc-1).at(1).at(0);
        e_temp = c2_edge.at(c2_fr);
        e_temp2 = c2_edge.at(c2_sc);
        //std::cout << "Edge2 " << e_temp << " , Edge2 " << e_temp2 
        //          << std::endl;

        if(c3_curr == -1) {
          ext_pc[j] = true;
          //k = ng_pt->first.size();
          // Rest of the triangles are to be attached to the final disc.
        }
        else {
          apf::ModelEntity* mdl_curr = vd_cd->get_mdl(3, c3_curr);
          assert(e_lens.vd_mesh_find_short(e_temp, e_temp2, 
                         &es_tri.at(j), e_lens.m->getModelTag(mdl_curr)));
        }
      }

      if (e_lens.c3_e == -1) {
      }
      else {
        int c2_e = ng_pt->second.back();
        apf::MeshEntity* e_temp;
        //e_temp = ent_list.e.at(2).at(c2_e-1).at(1).at(0);
        e_temp = c2_edge.at(c2_e);
        //std::cout << "Flag: " << e_lens.e_3c_b_flag << " ";
        //std::cout << "Last edge, 3c" << e_lens.c3_e 
        //          << " , 2c" << c2_e << std::endl;
        //std::cout << "Last edge, " << e_lens.e_3c_b << " , 2c edge " << e_temp 
        //          << std::endl;

        //std::cout << "Edge3 " << e_3 << " , Edge2 " << e_temp << std::endl;

        assert(e_lens.vd_mesh_find_short(e_lens.e_3c_b, e_temp, 
                                                &es_tri.at(j), e_lens.c3_e));
      }
      cp_tri(&es_tri.at(j), &e_lens.discs.at(j));
    }
  }
  //TODO check if the 3c has disjoint components. In that case, insert an edge
  // for each component.
  else if (trial_type == 3) {
/*
    std::vector<apf::MeshEntity*> es_1(0);
    std::vector<apf::MeshEntity*> es_2(0);
    std::vector<apf::MeshEntity*> es_3(0);

    vd_set_up(e_lens.m, e_lens.vert_ctr, &es_1);
    vd_set_up(e_lens.m, &es_1, &es_2);
    vd_set_up(e_lens.m, &es_2, &es_3);

    std::vector<apf::MeshEntity*> es_tet(0);
    std::vector<apf::MeshEntity*> es_surf(0);

    //es_tet = vd_3c->get_tet_id(trial_curr);

    es_tet.reserve(es_3.size());
    for(int i = 0; i < es_3.size(); i++) {
      apf::ModelEntity* mdl = e_lens.m->toModel(es_3.at(i));
      int m_type = e_lens.m->getModelType(mdl);
      int m_tag = e_lens.m->getModelTag(mdl);
      if(m_type == 3 and m_tag == trial_curr) {
        es_tet.push_back(es_3.at(i));
        //std::cout << es_3.at(i) 
        //          << m_type << "c" << m_tag << std::endl;
      }
    }
*/
    std::vector<apf::MeshEntity*> es_tet(0);
    std::vector<apf::MeshEntity*> es_tri(0);
    es_tet = vd_3c->get_tet_id(trial_curr);
    es_tri = coll_bound_tri(es_tet);

    cp_tri(&es_tri, &e_lens.discs.at(0));
  }
}

// Check the given set of triangles for members that 
// belong to a 3cell, 
// bound the center vertex and
// have edges adjacent to central vertex with lower dimensional cell membership. 
bool vd_edisc::chk_span_surf(std::vector<apf::MeshEntity*> &es_surf) {
  apf::Downward down;
  for(int i = 0; i < es_surf.size(); i++) {
    apf::MeshEntity* t_curr = es_surf.at(i);
    apf::ModelEntity* mdl = m_act->toModel(t_curr);
    int c2_type = m_act->getModelType(mdl);

    if(c2_type == 3) {
      m_act->getDownward(t_curr, 0, down);
      int v1 = findIn(down, 3, vert_ctr_act);

      if (v1 > -1) {
        m_act->getDownward(t_curr, 1, down);

        mdl = m_act->toModel(down[lookup_tri_ed[v1][0]]);
        int c1_type = m_act->getModelType(mdl);

        if(c1_type < c2_type) {
          mdl = m_act->toModel(down[lookup_tri_ed[v1][1]]);
          c1_type = m_act->getModelType(mdl);
          if(c1_type < c2_type) {
            return true;
          }
        }

        apf::MeshEntity* edge_cross = down[lookup_v_tri_e_x[v1]];
        mdl = m_act->toModel(edge_cross);
        c1_type = m_act->getModelType(mdl);
        m_act->getDownward(edge_cross, 0, down);

        if(c1_type > m_act->getModelType(m_act->toModel(down[0])) ) {
          if(c1_type > m_act->getModelType(m_act->toModel(down[1])) ) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

// Given a set of tetrahedra belonging to a disjoint 3stratum, collect the 
// bounding 2stratum triangles touching the central vertex. Calculate the 
// anglecosine weighted sum of area normals. Also, return false if the anglecosine
// between the direction and one of plane normals is lower than a threshold.
bool vd_edisc::coll_bound_tri_dir(std::vector<apf::MeshEntity*> &es_tet,
                                  apf::Vector3 &dir, double ang_th) {
  if(es_tet.size() == 0)
    return false;
  apf::ModelEntity* mdl_3c = m_act->toModel(es_tet.at(0));
  std::map<apf::MeshEntity*, apf::MeshEntity*> tri_2_tet_map{};
  std::map<apf::MeshEntity*, apf::Vector3> tet_dir_map{};
  apf::Up up;

  apf::Vector3 n(0,0,0);
  apf::Vector3 m(0,0,0);

  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_tri(0);
  apf::Downward d_v;

  vd_set_down(e_lens.m, &es_tet, &es_surf);
  es_tri.reserve(es_surf.size());

  for(int i = 0; i < es_tet.size(); i++) {
    tet_dir_map[es_tet.at(i)] = vd_get_pos(m_act, es_tet.at(i));
  }

  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* mdl = e_lens.m->toModel(es_surf.at(i));
    int m_type = e_lens.m->getModelType(mdl);
    int m_tag = e_lens.m->getModelTag(mdl);

    e_lens.m->getDownward(es_surf.at(i), 0, d_v);
    int j_v = findIn(d_v, 3, vert_ctr_act);
    if(m_type == 2 and j_v > -1) {
      es_tri.push_back(es_surf.at(i));
      //std::cout << e_lens.es_surf.at(i) 
      //          << m_type << "c" << m_tag << std::endl;
      m_act->getUp(es_surf.at(i), up);
      if(m_act->toModel(up.e[0]) == mdl_3c) {
        assert(findIn(&es_tet, es_tet.size(), up.e[0]) > -1);
        tri_2_tet_map[es_surf.at(i)] = up.e[0];
        tet_dir_map[es_surf.at(i)] = tet_dir_map[up.e[0]] 
                                    - vd_get_pos(m_act, es_surf.at(i));
      }
      else {
        assert(findIn(&es_tet, es_tet.size(), up.e[1]) > -1);
        assert(up.n == 2 and m_act->toModel(up.e[1]) == mdl_3c);
        tri_2_tet_map[es_surf.at(i)] = up.e[1];
        tet_dir_map[es_surf.at(i)] = tet_dir_map[up.e[1]] 
                                    - vd_get_pos(m_act, es_surf.at(i));
      }
    }
  }
  vd_disc disc_bound;
  disc_bound.m = m_act;
  disc_bound.v_ctr = vert_ctr_act;
  cp_tri(&es_tri, &disc_bound);

  apf::Vector3 norm_p(0,0,0);
  apf::Vector3 rel_dir0(0,0,0);
  apf::Vector3 rel_dir1(0,0,0);
  apf::Vector3 rel_dir2(0,0,0);

  for(int i = 0; i < disc_bound.tri.size(); i++) {
    apf::MeshEntity* tri = disc_bound.tri.at(i);
    int e_1;
    int e_2;

    e_1 = disc_bound.t_e_map1[tri] - 1;
    e_2 = disc_bound.t_e_map2[tri] - 1;

    // Relative directions of the edges preceeding and succeeding the current 
    // edge on the disc.
    rel_dir1 = disc_bound.v_pos.at(e_1) - pos_old;
    rel_dir2 = disc_bound.v_pos.at(e_2) - pos_old;
    rel_dir1 = norm_0(rel_dir1);
    rel_dir2 = norm_0(rel_dir2);
    m = cross(rel_dir1, rel_dir2);
    if(m * tet_dir_map[tri] < - std::numeric_limits<double>::min())
      m = m * (-1);

    double acos = std::acos(rel_dir1*rel_dir2);
    n = n + m*acos*0.5;
  }

  // Calc m:
  m = apf::Vector3(0,0,0);
  for (int i = 0; i < disc_bound.edge.size(); i++) {
    apf::MeshEntity* e_curr = disc_bound.edge.at(i);
    // Edge is not exterior:
    if(disc_bound.e_t_map1[e_curr] != 0 and 
       disc_bound.e_t_map2[e_curr] != 0) {
      int t_1 = disc_bound.e_t_map1[e_curr];
      int t_2 = disc_bound.e_t_map2[e_curr];
      apf::MeshEntity* tri_1 = disc_bound.tri.at(t_1 - 1);
      apf::MeshEntity* tri_2 = disc_bound.tri.at(t_2 - 1);

      int e_1;
      int e_2;

      if(disc_bound.t_e_map1[tri_1] != i+1) {
        e_1 = disc_bound.t_e_map1[tri_1] - 1;
      }
      else {
        assert(disc_bound.t_e_map1[tri_1] == i+1);
        e_1 = disc_bound.t_e_map2[tri_1] - 1;
      }
      if(disc_bound.t_e_map1[tri_2] != i+1) {
        e_2 = disc_bound.t_e_map1[tri_2] - 1;
      }
      else {
        assert(disc_bound.t_e_map1[tri_2] == i+1);
        e_2 = disc_bound.t_e_map2[tri_2] - 1;
      }
      // Relative directions of the edges preceeding and succeeding the current 
      // edge on the disc.
      rel_dir0 = disc_bound.v_pos.at(i) - pos_old;
      rel_dir1 = disc_bound.v_pos.at(e_1) - pos_old;
      rel_dir2 = disc_bound.v_pos.at(e_2) - pos_old;

      rel_dir0 = norm_0(rel_dir0);
      rel_dir1 = norm_0(rel_dir1);
      rel_dir2 = norm_0(rel_dir2);

      double acos = std::acos(rel_dir1*rel_dir0) + std::acos(rel_dir2*rel_dir0);
      m = m + rel_dir0*acos*0.5;
    }
  }
  std::cout << "n = [" << n[0] << ", " << n[1] << ", " << n[2] << "]; "
            << " m = [" << m[0] << ", " << m[1] << ", " << m[2] << "]; "
            << std::endl;
  if(n.getLength() > m.getLength())
    dir = norm_0(n);
  else
    dir = norm_0(m);
  assert(!std::isnan(dir.getLength()));

  if(es_tet.size() == 1) {
    assert((vd_get_pos(e_lens.m, es_tet.at(0)) - pos_old)*dir > 
                            std::numeric_limits<double>::min());

  }
  for(int i = 0; i < disc_bound.tri.size(); i++) {
    apf::MeshEntity* tri = disc_bound.tri.at(i);
    int e_1;
    int e_2;

    e_1 = disc_bound.t_e_map1[tri] - 1;
    e_2 = disc_bound.t_e_map2[tri] - 1;

    // Relative directions of the edges preceeding and succeeding the current 
    // edge on the disc.
    rel_dir1 = disc_bound.v_pos.at(e_1) - pos_old;
    rel_dir2 = disc_bound.v_pos.at(e_2) - pos_old;
    rel_dir1 = norm_0(rel_dir1);
    rel_dir2 = norm_0(rel_dir2);
    m = cross(rel_dir1, rel_dir2);
    if(m * tet_dir_map[tri] <- std::numeric_limits<double>::min())
      m = m * (-1);
    if(m*dir > ang_th)
      return false;
  }
  return true;
}

// Given a set of tetrahedra belonging to a disjoint 3stratum, collect the 
// bounding 2stratum triangles.
std::vector<apf::MeshEntity*> vd_edisc::coll_bound_tri(
                          std::vector<apf::MeshEntity*> &es_tet) {
  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_tri(0);
  apf::Downward d_v;

  vd_set_down(e_lens.m, &es_tet, &es_surf);
  es_tri.reserve(es_surf.size());

  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* mdl = e_lens.m->toModel(es_surf.at(i));
    int m_type = e_lens.m->getModelType(mdl);
    int m_tag = e_lens.m->getModelTag(mdl);

    e_lens.m->getDownward(es_surf.at(i), 0, d_v);
    int j_v = findIn(d_v, 3, vert_ctr_act);
    if(m_type == 2 and j_v > -1) {
      es_tri.push_back(es_surf.at(i));
      //std::cout << e_lens.es_surf.at(i) 
      //          << m_type << "c" << m_tag << std::endl;
    }
  }
  return es_tri;
}

void vd_edisc::color_discs() {

  if(trial_type == 1) {

    // The disc(s) divide the entities into segments.
    // Given the discs, for every entity adjacent to the central vertex, 
    // find the segment it belongs to. For the disc entities, they will be 
    // associated with an old and a new modelentity and two vertices. 

    //e_lens.null_slice_map();

    tri_burn.clear();
    color_disc(0);

    burn_trial();

  }
  else if(trial_type == 2) {

    // On the disc:
    // For each triangle, look at the neighboring tetrahedra. The tetrahedra  
    // belonging to the cone is the top tetrahedra, into which the lens will 
    // expand to is stored.
    // The other tet is stored in the bottom list.
    std::cout << "Discs size " << e_lens.discs.size() << std::endl;
    for(int j = 0; j < e_lens.discs.size(); j++) {

      for (int i = 0; i < e_lens.discs.at(j).tri.size(); i++) {
        apf::Up up;
        e_lens.discs.at(j).m->getUp(e_lens.discs.at(j).tri.at(i), up);
        assert(up.n == 2);

        // 3cell model.
        apf::ModelEntity* mdl = e_lens.discs.at(j).m->toModel(up.e[0]);
        e_lens.discs.at(j).t_em.at(i) = mdl;
      }

      // For each edge, look at the cell membership. Assuming precondtioning,  
      // it must belong to a 2cell or 3cell, and the associated lens triangle  
      // belongs to the cell same cell. 
      for (int i = 0; i < e_lens.discs.at(j).edge.size(); i++) {
        apf::ModelEntity* e_m = 
                e_lens.discs.at(j).m->toModel(e_lens.discs.at(j).edge.at(i));
        int type_e = e_lens.m->getModelType(e_m);
        int tag_e = e_lens.m->getModelTag(e_m);

        //std::cout << disc_in->edge.at(i) << " " << type_e << "c" << tag_e 
        //          << std::endl;

        // assert(0 < type_e and type_e < 3);

        e_lens.discs.at(j).e_em.at(i) = e_m;
        assert(type_e > 1);
        e_lens.discs.at(j).t2_em.at(i) = e_lens.discs.at(j).e_em.at(i);
      }

      color_disc(j);
    }

    if (e_lens.c3_f != -1) {
      e_lens.slice_map[e_lens.e_3c_t] = -e_lens.discs.size()-1;
      //std::cout << "Top " << e_lens.e_3c_t << " " 
      //          << e_lens.slice_map[e_lens.e_3c_t] << std::endl;
    }

    if (e_lens.c3_e != -1) {
      e_lens.slice_map[e_lens.e_3c_b] = -e_lens.discs.size()-1;
      //std::cout << "Bot " << e_lens.e_3c_b << " " 
      //          << e_lens.slice_map[e_lens.e_3c_b] << std::endl;

    }

    //for(int j = 0; j < ng.ngons.at(ng.curr).size(); j++) {
    //  std::vector<std::vector<apf::MeshEntity*> > e_set
    //             (4, std::vector<apf::MeshEntity*> (0) );
      
    //  e_set.at(2) = e_lens.discs.at(j).tri;
    //  vd_set_up(e_lens.m, &e_set.at(2), &e_set.at(3));
    //  vd_set_down(e_lens.m, &e_set.at(3), &e_set.at(2));
    //  vd_set_down(e_lens.m, &e_set.at(2), &e_set.at(1));
    //  vd_set_down(e_lens.m, &e_set.at(1), &e_set.at(0));

    //  char s[50];
    //  sprintf(s, "./output/tri_disc%d", j);
    //  vd_save_vtk_set(e_lens.m, &e_set, s);

    //  for(int k = 0; k < 4; k++) {
    //    e_set.at(k).clear();
    //  }
    //}

    burn_trial();
  }

  else if(trial_type == 3) {
    e_lens.slice_map.clear();

    tri_burn.clear();
    color_disc(0);

    burn_trial();

    //std::cout << "Slices: "
    //          << e_lens.slices.at(0) << " "
    //          << e_lens.slices.at(1) << std::endl;

    apf::ModelEntity* mdl = vd_3c->get_mdl(3, trial_curr);
    //apf::ModelEntity* mdl = e_lens.m->findModelEntity(3, trial_curr);
    for (int i = 0; i < e_lens.discs.at(0).tri.size(); i++) {
      e_lens.discs.at(0).t_em.at(i) = mdl;
    }
    for (int i = 0; i < e_lens.discs.at(0).edge.size(); i++) {
      e_lens.discs.at(0).e_em.at(i) = mdl;
      e_lens.discs.at(0).t2_em.at(i) = mdl;
    }
  }
  e_lens.null_flag = false;
}

// Calculate the velocities of the new 0cell vertices using the Mason algorithm.
void vd_edisc::calc_ctr_vel(bool fix) {

  apf::Field* vel_field = e_lens.m->findField("velocity_field");
  std::vector<std::vector<apf::MeshEntity*> > tri_slice
            (e_lens.vert_ctr_new.size(), std::vector<apf::MeshEntity*> (0) );

  e_lens.vt.clear();
  e_lens.vt.resize(e_lens.vert_ctr_new.size());

  collect_slice_tris(&tri_slice);
  //for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
  //  collect_slice_tri(&tri_slice.at(i), i+1);
  //}

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    f_calc.vd_upd_vel_field(e_lens.m, vert_ctr_act, &tri_slice.at(i), fix);
    //vd_upd_vel_field_mason(e_lens.m, vert_ctr_act, &tri_slice.at(i), &f_calc);
    apf::getVector(vel_field, vert_ctr_act, 0, e_lens.vt.at(i));
    //std::cout << "calc_ctr[" << i << "] " << e_lens.vt.at(i) << std::endl;
    //apf::setVector(vel_field, e_lens.vert_ctr_new.at(i), 0, e_lens.vt.at(i));
  }
  //apf::Vector3 v_avg = calc_v_avg();
  //for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
  //  apf::getVector(vel_field, vert_ctr_act, 0, e_lens.vt.at(i));
  //  e_lens.vt.at(i) = e_lens.vt.at(i) - v_avg;
  //  std::cout << "calc_ctr_rel[" << i << "] " << e_lens.vt.at(i) 
  //            << std::endl;
  //}

  for (int i = 0; i < tri_slice.size(); i++) {
    tri_slice.at(i).clear();
  }
  tri_slice.clear();
}

// Calculate the velocities of the new 1cell vertices using the Mason algorithm.
void vd_edisc::calc_sp_vel(bool fix) {

  apf::Field* vel_field = e_lens.m->findField("velocity_field");
  std::vector<std::vector<apf::MeshEntity*> > tri_path
                 (0, std::vector<apf::MeshEntity*> (0) );
  e_lens.split_vt.clear();
  e_lens.split_vt.resize(ng.ngons.at(ng.curr).size());
  tri_path.resize(ng.ngons.at(ng.curr).size());

  for (int i = 0; i < ng.ngons.at(ng.curr).size(); i++) {
    collect_path_tri(&tri_path.at(i), i);
  }

  for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    f_calc.vd_upd_vel_field(e_lens.m, vert_ctr_act, &tri_path.at(i), fix);
    //vd_upd_vel_field_mason(e_lens.m, vert_ctr_act, &tri_path.at(i), &f_calc);
    apf::getVector(vel_field, vert_ctr_act, 0, e_lens.split_vt.at(i));
    //std::cout << "calc_sp[" << i << "] " << e_lens.split_vt.at(i) 
    //          << std::endl;
    //apf::setVector(vel_field, e_lens.vert_sp_new.at(i), 0, 
    //                                        e_lens.split_vt.at(i));
  }

  // Remove the velocity component in the direction of relative motion of the 
  // center vertices to expand laterally, and the average component of split 
  // vertices to expand outwardly. 
  //apf::Vector3 v_sp_avg(0,0,0);
  apf::Vector3 v_ctr_dir(0,0,0);

  for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    int v1 = e_lens.slices.at(2*i)-1;
    int v2 = e_lens.slices.at(2*i+1)-1;
    v_ctr_dir = norm_0(e_lens.vt.at(v1) - e_lens.vt.at(v2) );
    e_lens.split_vt.at(i) = e_lens.split_vt.at(i) 
                          - v_ctr_dir*(e_lens.split_vt.at(i)*v_ctr_dir);
    //v_sp_avg = v_sp_avg + e_lens.split_vt.at(i);
  }
  //v_sp_avg = v_sp_avg/e_lens.vert_sp_new.size();
  //for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
  //  e_lens.split_vt.at(i) = e_lens.split_vt.at(i) - v_sp_avg;
  //}

  //apf::Vector3 v_avg = calc_v_avg();
  //for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
  //  apf::getVector(vel_field, vert_ctr_act, 0, e_lens.split_vt.at(i));
  //  e_lens.split_vt.at(i) = e_lens.split_vt.at(i) - v_avg;
  //  std::cout << "calc_sp_rel[" << i << "] " << e_lens.split_vt.at(i) 
  //            << std::endl;
  //}

  for (int i = 0; i < tri_path.size(); i++) {
    tri_path.at(i).clear();
  }
  tri_path.clear();
}

void vd_edisc::clear_vel_map() {

  ext_proj_type.clear();
  for(int i = 0; i < ext_proj_dir.size(); i++) {
    ext_proj_dir.at(i).clear();
  }
  ext_proj_dir.clear();
}

void vd_edisc::coll_proj_map() {
/*
use slices
each slice will have either none, one or two exterior cells on its path, 
each slice will have a number of side edges. collect the directions of the side 
edges if there are any
if a new 1cell has two slices that are corner, use their projection dirs
if one has a lower number of projection dirs, use that one
*/
  clear_vel_map();
  ext_proj_dir.resize(e_lens.vert_ctr_new.size());

  for(int i = 0; i < e_lens.es_edge.size(); i++) {
    apf::MeshEntity* ent_c = e_lens.es_edge.at(i);
    //std::cout << "Check slice id of " << ent_c << std::endl;
    int slice_id = e_lens.slice_map[ent_c];
    if(slice_id > 0) {
      apf::ModelEntity* m_edge = e_lens.m->toModel(ent_c);
      int m_dim = e_lens.m->getModelType(m_edge);
      int m_id = e_lens.m->getModelTag(m_edge);
      if(m_dim == 1) {
        if(e_lens.c_base_curr->get_1c_corner_gmi(m_id)) {
          ext_proj_type[slice_id] = 1;
          ext_proj_dir.at(slice_id-1).push_back(get_edge_dir(e_lens.m, ent_c));
        }
      }
    }
  }

  for(int i = 0; i < e_lens.es_surf.size(); i++) {
    apf::MeshEntity* ent_c = e_lens.es_surf.at(i);
    //std::cout << "Check slice id of " << ent_c << std::endl;
    int slice_id = e_lens.slice_map[ent_c];
    //std::cout << slice_id << std::endl;

    if(slice_id > 0) {
      if(ext_proj_type[slice_id] == 0) {
        apf::ModelEntity* m_surf = e_lens.m->toModel(ent_c);
        int m_dim = e_lens.m->getModelType(m_surf);
        int m_id = e_lens.m->getModelTag(m_surf);
        if(m_dim == 2 and e_lens.c_base_curr->get_cell_ext_gmi(2, m_id)) {
          ext_proj_type[slice_id] = 2;
          ext_proj_dir.at(slice_id-1).push_back(
                                  norm_0(vd_area_out(e_lens.m, ent_c)));
        }
      }
    }
  }
}

apf::Vector3 vd_edisc::calc_v_avg() {
  apf::Vector3 v_avg(0,0,0);

  if(!ext_0cell or 
     (trial_type == 2 and (e_lens.c3_f == -1 or e_lens.c3_e == -1))) {
    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      //vd_upd_vel_field_mason(e_lens.m, e_lens.vert_ctr_new.at(i), &f_calc);
      v_avg = v_avg + e_lens.vt.at(i);
    }
    v_avg = v_avg/e_lens.vert_ctr_new.size();
  }
  else {
    bool sp_flag = false;
    bool ctr_flag = false;
    int v_id = -1;
    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      if(ext_pc[i]) {
        sp_flag = true;
        v_id = i;
        i = e_lens.vert_sp_new.size();
      }
    }
    if(!sp_flag) {
      for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        if(ext_slice[i]) {
          ctr_flag = true;
          v_id = i;
          i = e_lens.vert_ctr_new.size();
        }
      }
    }
    assert(sp_flag or ctr_flag);
    if(ctr_flag) {
      v_avg = e_lens.vt.at(v_id);
    }
    else {
      int v1 = e_lens.slices.at(2*v_id)-1;
      int v2 = e_lens.slices.at(2*v_id+1)-1;
      v_avg = e_lens.vt.at(v1)/2 + e_lens.vt.at(v2)/2;
    }
  }
  return v_avg;
}

apf::Vector3 vd_edisc::calc_v_avg(std::map<apf::MeshEntity*, bool> sk_map) {
  apf::Vector3 v_avg(0,0,0);

  if(!ext_0cell or 
     (trial_type == 2 and (e_lens.c3_f == -1 or e_lens.c3_e == -1))) {
    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      //vd_upd_vel_field_mason(e_lens.m, e_lens.vert_ctr_new.at(i), &f_calc);
      v_avg = v_avg + e_lens.vt.at(i);
    }
    v_avg = v_avg/e_lens.vert_ctr_new.size();
  }
  else {
    bool sp_flag = false;
    bool ctr_flag = false;
    int v_id = -1;
    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      if(sk_map[e_lens.vert_sp_new.at(i)]) {
        sp_flag = true;
        v_id = i;
        i = e_lens.vert_sp_new.size();
      }
    }
    if(!sp_flag) {
      for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        if(sk_map[e_lens.vert_ctr_new.at(i)]) {
          ctr_flag = true;
          v_id = i;
          i = e_lens.vert_ctr_new.size();
        }
      }
    }
    assert(sp_flag or ctr_flag);
    if(ctr_flag) {
      v_avg = e_lens.vt.at(v_id);
    }
    else {
      int v1 = e_lens.slices.at(2*v_id)-1;
      int v2 = e_lens.slices.at(2*v_id+1)-1;
      v_avg = e_lens.vt.at(v1)/2 + e_lens.vt.at(v2)/2;
    }
  }
  return v_avg;
}

void vd_edisc::fix_vel_new() {
  apf::Vector3 v_avg(0,0,0);
  apf::Vector3 v_res(0,0,0);
  apf::Vector3 v_miss(0,0,0);
  apf::Vector3 v_temp(0,0,0);

//first fix the exterior velocities, than calc average and add the missing parts
//v_avg should be the velocity of the exterior vertex, belonging to the dominant shell

  int ext_count = 0;

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    if(ext_proj_type[i+1] == 1) {
      v_res = e_lens.vt.at(i);
      for(int j = 0; j < ext_proj_dir.at(i).size(); j++) {
        assert(ext_proj_dir.at(i).at(j).getLength() > 
                              std::numeric_limits<double>::epsilon());
        e_lens.vt.at(i) = ext_proj_dir.at(i).at(j)*
                                  (e_lens.vt.at(i)*ext_proj_dir.at(i).at(j));
      }
      v_miss = v_miss - e_lens.vt.at(i) + v_res;
      ext_count = ext_count + 1;
    }
    else if(ext_proj_type[i+1] == 2) {
      v_res = e_lens.vt.at(i);
      e_lens.vt.at(i) = e_lens.vt.at(i) - ext_proj_dir.at(i).at(0)*
                                  (e_lens.vt.at(i)*ext_proj_dir.at(i).at(0));
      v_miss = v_miss - e_lens.vt.at(i) + v_res;
      ext_count = ext_count + 1;
    }
  }

  ext_count = e_lens.vert_ctr_new.size() - ext_count;
  if(ext_count > 0) {
    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      if(ext_proj_type[i+1] == 0) {
        e_lens.vt.at(i) = e_lens.vt.at(i) + v_miss/ext_count;
      }
    }
  }

  if(trial_type == 2) {
    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      int v1 = e_lens.slices.at(2*i)-1;
      int v2 = e_lens.slices.at(2*i+1)-1;
      if(ext_pc[i]) {
        if(ext_proj_type[v1+1] == 1 and ext_proj_type[v2+1] == 1) {
          ext_cor[i] = true;
          assert(ext_proj_dir.at(v1).size() == 1 or 
                 ext_proj_dir.at(v2).size() == 1);
          if(ext_proj_dir.at(v1).size() == 1) {
            e_lens.split_vt.at(i) = ext_proj_dir.at(v1).at(0)*
                             (e_lens.split_vt.at(i)*ext_proj_dir.at(v1).at(0));

            v_temp =  ext_proj_dir.at(v1).at(0)*
                             (v_avg*ext_proj_dir.at(v1).at(0));
            e_lens.split_vt.at(i) = e_lens.split_vt.at(i) - v_temp;
          }
          if(ext_proj_dir.at(v2).size() == 1) {
            e_lens.split_vt.at(i) = ext_proj_dir.at(v2).at(0)*
                             (e_lens.split_vt.at(i)*ext_proj_dir.at(v2).at(0));
            v_temp =  ext_proj_dir.at(v2).at(0)*
                             (v_avg*ext_proj_dir.at(v2).at(0));
            e_lens.split_vt.at(i) = e_lens.split_vt.at(i) - v_temp;
          }
        }
        else if (ext_proj_type[v1+1] == 2 and ext_proj_type[v2+1] == 2) {
          e_lens.split_vt.at(i) = e_lens.split_vt.at(i) 
                            - ext_proj_dir.at(v1).at(0)*
                             (e_lens.split_vt.at(i)*ext_proj_dir.at(v1).at(0));

          v_temp = v_temp - ext_proj_dir.at(v1).at(0)*
                             (v_avg*ext_proj_dir.at(v1).at(0));
          e_lens.split_vt.at(i) = e_lens.split_vt.at(i) - v_temp;
        }
      }
      else if(e_lens.c3_f == -1 or e_lens.c3_e == -1) {
        assert(ext_proj_type[v1+1] == 2 and ext_proj_type[v2+1] == 2);

        e_lens.split_vt.at(i) = e_lens.split_vt.at(i) 
                                - ext_proj_dir.at(v1).at(0)*
                           (e_lens.split_vt.at(i)*ext_proj_dir.at(v1).at(0));
        v_temp = v_avg - ext_proj_dir.at(v1).at(0)*
                           (v_avg*ext_proj_dir.at(v1).at(0));
        e_lens.split_vt.at(i) = e_lens.split_vt.at(i) - v_temp;
      }
      else
        e_lens.split_vt.at(i) = e_lens.split_vt.at(i) - v_avg;
    }
  }
  v_avg = calc_v_avg();

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.vt.at(i) = e_lens.vt.at(i) - v_avg;
  }

}

void vd_edisc::find_ctr_vel() {
  apf::Field* vel_field = e_lens.m->findField("velocity_field");
  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    f_calc.vd_upd_vel_field(e_lens.m,  e_lens.vert_ctr_new.at(i), drag_flag);
    apf::setVector(vel_field, e_lens.vert_ctr_new.at(i), 0, e_lens.vt.at(i));
  }
}

void vd_edisc::find_sp_vel() {
  apf::Field* vel_field = e_lens.m->findField("velocity_field");
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    f_calc.vd_upd_vel_field(e_lens.m,  e_lens.vert_sp_new.at(i), drag_flag);
    apf::setVector(vel_field, e_lens.vert_sp_new.at(i), 0, 
                                              e_lens.split_vt.at(i));
  }
}

void vd_edisc::find_max_disp() {
  vel_max = 0;
  double v_curr;

  apf::Field* vel_field = e_lens.m->findField("velocity_field");
  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    v_curr = e_lens.vt.at(i).getLength();
    if(v_curr > vel_max)
      vel_max = v_curr;
  }
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    v_curr = e_lens.split_vt.at(i).getLength();
    if(v_curr > vel_max)
      vel_max = v_curr;
  }
  assert(!std::isnan(vel_max));
  if (vel_max < std::numeric_limits<double>::epsilon() ) {
    vel_max = 1;
  }
  std::cout << "vel_max " << vel_max 
            << " len_sp " << len_sp << std::endl;

  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    std::cout << i << " " << e_lens.vt.at(i) << std::endl;
    e_lens.vt.at(i) = e_lens.vt.at(i)/vel_max*len_sp/rat_init;
  }
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.split_vt.at(i) = e_lens.split_vt.at(i)/vel_max*len_sp/rat_init;
  }
}

// Evolve the inserted shape until all vertices are expanding. This is achieved
// by sequentially moving the vertices, calculating the center of the 
// expanding sphere projecting onto the expanding sphere, until the condition 
// is reached. 
void vd_edisc::evolve_ins() {

  assert(len_sp > std::numeric_limits<double>::epsilon());
  // TODO Dirty way of overcoming the jumping over the minimum energy position
  // due to zero mean velocity expansion. A better approach should be by 
  // setting the objective as finding the minimum energy positions for the 
  // given insertion and calculating a mean radius of expansion to use in 
  // calculating the energy rate. This would be necessary for anisotropic 
  // properties.
  double rho = len_sp/rho_rat/fudge_factor;
  double disp_th = rho/rho_rat/1.2/fudge_factor;

  double rad_tol = disp_th/6/fudge_factor;

  double r_th_min = rho/200/fudge_factor;
  double r_th_max = 2*rho/fudge_factor;

  // Displacement threshold. If the maximum motion is less than this amount
  // stop iteration.
  double d_th = disp_th/2400/fudge_factor;

  double r_min = rho;
  double r_max = 0;

  double ang_th;
  if(trial_type == 1)
    ang_th = 2*pi/18;
  else
    ang_th = 2*pi/e_lens.vert_ctr_new.size()/6;

  //double dir_tol = rho*ang_th;

  double max_ang = 6.14/f_calc.vdparam.dt;

  double v_r_curr = 0;
  double v_t_curr = 0;
  double r_curr = 0;

  double v_r_max = 0;

  int counter_sz = 50;
  double rat_t = 3;

  std::cout << "rho " << rho << " len_sp " << len_sp << std::endl;
  std::cout << "ang_th " << ang_th << " max_ang " << max_ang << std::endl;

  int iter_sz = 4;

  double dt_glob = f_calc.vdparam.dt;

  apf::Field* vel_field = e_lens.m->findField("velocity_field");

  apf::Vector3 pos_curr(0,0,0);
  apf::Vector3 r_hat(0,0,0);
  apf::Vector3 t_hat(0,0,0);
  apf::Vector3 v_temp(0,0,0);

  std::vector<apf::Vector3> v_r_ctr(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_r_sp(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> v_r_ctr_old(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_r_sp_old(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> v_t_ctr(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_t_sp(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> v_t_ctr_old(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_t_sp_old(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> v_ctr_old(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_sp_old(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> v_ctr(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_sp(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> r_ctr(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> r_sp(0, apf::Vector3(0,0,0));

  int v_sz = e_lens.vert_ctr_new.size()+e_lens.vert_sp_new.size();
  int v_ctr_sz = e_lens.vert_ctr_new.size();
  int v_sp_sz = e_lens.vert_sp_new.size();

  r_ctr.resize(e_lens.vert_ctr_new.size());
  r_sp.resize(e_lens.vert_sp_new.size());

  v_r_ctr.resize(e_lens.vert_ctr_new.size());
  v_t_ctr.resize(e_lens.vert_ctr_new.size());

  v_r_sp.resize(e_lens.vert_sp_new.size());
  v_t_sp.resize(e_lens.vert_sp_new.size());

  v_r_ctr_old.resize(e_lens.vert_ctr_new.size());
  v_r_sp_old.resize(e_lens.vert_sp_new.size());

  v_t_ctr_old.resize(e_lens.vert_ctr_new.size());
  v_t_sp_old.resize(e_lens.vert_sp_new.size());

  v_ctr.resize(e_lens.vert_ctr_new.size());
  v_sp.resize(e_lens.vert_sp_new.size());

  v_ctr_old.resize(e_lens.vert_ctr_new.size());
  v_sp_old.resize(e_lens.vert_sp_new.size());

  //The exterior vertex will be positioned at the old location, so the closing
  // check will be positive. Skip that vertex.
  //std::map<apf::MeshEntity*, bool> skip_dist{};
  //skip_dist = calc_skip();

  double z[3] = {0,0,0};

  double dt_th;
  //double dt_inv;
  double dt_t = -1;
  double dt_r = -1;

/*
  std::vector<apf::MeshEntity*> ent_col(0);
  if(trial_type == 1)
    vd_find_ent_topo(e_lens.m, e_lens.vert_sp_new.at(0), &ent_col, 
                                               e_lens.new_cell1_id.at(0), 1);
  else
    vd_find_ent_topo(e_lens.m, e_lens.v_2c, &ent_col, 
                                               e_lens.new_cell2_id, 2);

  fix_skip(&ent_col, skip_dist);
*/

  proj_around_ctr(rho);

  int mov_tag = 0;
  bool closing = false;
  bool all_pos = false;

  int iter = 0;

  std::vector<std::vector<apf::MeshEntity*> > ent(4, 
                                          std::vector<apf::MeshEntity*>(0));
  if(trial_type == 1) {
    ent.at(0).resize(3);
    ent.at(0).at(0) = e_lens.vert_ctr_new.at(0);
    ent.at(0).at(1) = e_lens.vert_ctr_new.at(1);
    ent.at(0).at(2) = e_lens.vert_sp_new.at(0);
  }
  else {
    assert(trial_type == 2);
    ent.at(0).resize(e_lens.vert_ctr_new.size()*2+1);
    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++)
      ent.at(0).at(i) = e_lens.vert_ctr_new.at(i);
    for (int i = 0; i < e_lens.vert_sp_new.size(); i++)
      ent.at(0).at(e_lens.vert_ctr_new.size()+i) = e_lens.vert_sp_new.at(i);
    ent.at(0).at(e_lens.vert_ctr_new.size()*2) = e_lens.v_2c;
  }
  for(int dim = 0; dim < 3; dim++) {
    vd_set_up(e_lens.m, &ent.at(dim), &ent.at(dim+1));
  }

  std::vector<apf::MeshEntity*> tri(0);
  tri.reserve(ent.at(3).size());
  apf::Downward d_t;
  apf::Downward d_v;
  int lookup_ts [4] = {2, 3, 1, 0};

  for(int i = 0; i < ent.at(3).size(); i++) {
    e_lens.m->getDownward(ent.at(3).at(i), 0, d_v);
    bool found1 = false;
    bool found2 = false;

    apf::MeshEntity* tri_temp;
    for(int j = 0; j < 4; j++) {
      int i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[j]);
      if(i1 > -1) {
        if(!found1) {
          e_lens.m->getDownward(ent.at(3).at(i), 2, d_t);
          tri_temp = d_t[lookup_ts[j]];
          found1 = true;
        }
        else {
          found2 = true;
          j = 4;
        }
      }
    }
    assert(found1);
    if(!found2) {
      tri.push_back(tri_temp);
    }
  }

  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    //vd_upd_vel_field_mason(e_lens.m, e_lens.vert_ctr_new.at(i), &f_calc);
    f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr_new.at(i), drag_flag);
    apf::getVector(vel_field, e_lens.vert_ctr_new.at(i), 0, e_lens.vt.at(i));
  }

  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    //vd_upd_vel_field_mason(e_lens.m, e_lens.vert_ctr_new.at(i), &f_calc);
    f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_sp_new.at(i), drag_flag);
    apf::getVector(vel_field, e_lens.vert_sp_new.at(i), 0, 
                                                    e_lens.split_vt.at(i));
  }

  double rho_st = rho;
  while((!all_pos and iter < iter_sz) and rho > rho_st/16) {

    t_total = 0;

    double dt_t_curr = -1;
    double dt_r_curr = -1;

    apf::Vector3 v_avg(0,0,0);
    // The components of the average velocities not removed from the vertices.
    int ext_count = 0;
    apf::Vector3 v_miss(0,0,0);
    apf::Vector3 v_res(0,0,0);

    double dt_inv_ctr = f_calc.find_min_t2(e_lens.m, &e_lens.vert_ctr_new);
    double dt_inv_sp = f_calc.find_min_t2(e_lens.m, &e_lens.vert_sp_new);

    double dt = std::min(dt_inv_ctr, dt_inv_sp);

    if(dt < 0)
      dt = f_calc.vdparam.dt;
    else
      dt = dt/9;

    apf::Vector3 pos_curr(0,0,0);

    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, pos_curr);
      pos_curr = pos_curr + e_lens.vt.at(i)*dt;
      e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_curr);
    }
    if(trial_type == 1) {
      e_lens.pos = vd_get_center(e_lens.m, &e_lens.vert_ctr_new);
      e_lens.m->setPoint(e_lens.vert_sp_new.at(0), 0, e_lens.pos);
    }
    else {
      assert(trial_type == 2);
      if(e_lens.vert_ctr_new.size() == 2) {
        for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
          e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, pos_curr);
          pos_curr = pos_curr + e_lens.split_vt.at(i)*dt;
          e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_curr);
        }
        e_lens.pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
                vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;
        e_lens.m->setPoint(e_lens.v_2c, 0, e_lens.pos);
      }
      else {
        for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
          int v1 = e_lens.slices.at(2*i)-1;
          int v2 = e_lens.slices.at(2*i+1)-1;

          apf::Vector3 temp1;
          apf::Vector3 temp2;

          e_lens.m->getPoint(e_lens.vert_ctr_new.at(v1), 0, temp1);
          e_lens.m->getPoint(e_lens.vert_ctr_new.at(v2), 0, temp2);
          e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, temp1/2+temp2/2);
        }
        e_lens.pos = vd_get_center(e_lens.m, &e_lens.vert_ctr_new);
        e_lens.m->setPoint(e_lens.v_2c, 0, e_lens.pos);

      }
    }
    t_total = t_total + dt;

    double r_max = -1;
    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      if(r_max < 0 or r_max < r_curr)
        r_max = r_curr;
    }
    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      if(r_max < 0 or r_max < r_curr)
        r_max = r_curr;

    }
    double len_max;

    std::tie(rho, len_max) = tri_dist_sphere(e_lens.m, e_lens.pos, &tri);

    proj_around_ctr(rho/4);

    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      //vd_upd_vel_field_mason(e_lens.m, e_lens.vert_ctr_new.at(i), &f_calc);
      f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr_new.at(i), drag_flag);
      apf::getVector(vel_field, e_lens.vert_ctr_new.at(i), 0, e_lens.vt.at(i));
    }

    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      //vd_upd_vel_field_mason(e_lens.m, e_lens.vert_ctr_new.at(i), &f_calc);
      f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_sp_new.at(i), drag_flag);
      apf::getVector(vel_field, e_lens.vert_sp_new.at(i), 
                                                   0, e_lens.split_vt.at(i));
    }
    iter = iter + 1;
  }
  if(sub_vtk)
    vd_save_vtk_vert(e_lens.m, &e_lens.vert_sp_new, "output/inserted");
}

void vd_edisc::move_ctr() {

  apf::Field* vel_field = e_lens.m->findField("velocity_field");
  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    //e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, e_lens.pos+e_lens.vt.at(i)*f_calc.vdparam.dt);
    //e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, e_lens.pos + 
    //                                   norm_0(e_lens.vt.at(i))*len_sp/rat_init);
    e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, e_lens.pos+
                                                        e_lens.vt.at(i));
    std::cout << "Vert_ctr " << e_lens.vert_ctr_new.at(i)  << ";" 
              << std::endl;
    std::cout << "\t" << "disp_ctr(" << i + 1 << ",:) = " 
              << e_lens.vt.at(i)  
              << ";" << std::endl;
    //std::cout << "\t" << "new_pos_ctr(" << i + 1 << ",:) = "
              //<< e_lens.pos+e_lens.vt.at(i)*f_calc.vdparam.dt  << ";" << std::endl;
    //          << e_lens.pos+ e_lens.vt.at(i) << ";" << std::endl;

    assert (sqrt(e_lens.vt.at(i)*e_lens.vt.at(i)) < len_sp);
  }
}

// Calculate the velocities of the new 1cell vertices using the Mason algorithm.
void vd_edisc::move_sp() {

  apf::Field* vel_field = e_lens.m->findField("velocity_field");
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    //e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, e_lens.pos
    //                                      +e_lens.split_vt.at(i)*f_calc.vdparam.dt);
    //std::cout << "Vert_sp " << e_lens.vert_sp_new.at(i) << std::endl;
    if(!ext_new[e_lens.vert_sp_new.at(i)]) {
      //std::cout << "Not ext" << std::endl;
      //e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, e_lens.pos
      //                          +norm_0(e_lens.split_vt.at(i))*len_sp/rat_init);
      e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, e_lens.pos
                                + e_lens.split_vt.at(i));
      //std::cout << "\t" << "disp_sp(" << i + 1 << ",:) = " 
      //          << e_lens.split_vt.at(i) << ";"
      //          << std::endl;
      //std::cout << "\t" << "new_pos_sp(" << i + 1 << ",:) = " 
      //          << e_lens.pos + e_lens.split_vt.at(i) 
      //          << std::endl << ";";
      assert (sqrt(e_lens.split_vt.at(i)*e_lens.split_vt.at(i)) < len_sp);
    }

    else {
      //std::cout << "Ext" << std::endl;

      int v1 = e_lens.slices.at(2*i)-1;
      int v2 = e_lens.slices.at(2*i+1)-1;

      apf::Vector3 temp1;
      apf::Vector3 temp2;

      e_lens.m->getPoint(e_lens.vert_ctr_new.at(v1), 0, temp1);
      e_lens.m->getPoint(e_lens.vert_ctr_new.at(v2), 0, temp2);

      //std::cout << "Vert_sp " << e_lens.vert_sp_new.at(i) << std::endl;
      e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, (temp1+temp2)/2);
    }
  }
  if(e_lens.vert_ctr_new.size() > 2) {
    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      int v1 = e_lens.slices.at(2*i)-1;
      int v2 = e_lens.slices.at(2*i+1)-1;

      apf::Vector3 temp1;
      apf::Vector3 temp2;

      e_lens.m->getPoint(e_lens.vert_ctr_new.at(v1), 0, temp1);
      e_lens.m->getPoint(e_lens.vert_ctr_new.at(v2), 0, temp2);

      //std::cout << "Vert_sp " << e_lens.vert_sp_new.at(i) << std::endl;
      e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, (temp1+temp2)/2);
    }
  }
}

void vd_edisc::fix_skip(std::vector<apf::MeshEntity*>* ent_col, 
                              std::map<apf::MeshEntity*, bool> skip_map) {
  if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL and ext_0cell) {
    if(trial_type == 1) {
      if(skip_map[e_lens.vert_ctr_new.at(0)] or 
         skip_map[e_lens.vert_ctr_new.at(1)]) {
        apf::Downward down;
        for(int i = ent_col->size()-1; i > -1; i--) {
          int d_sz = e_lens.m->getDownward(ent_col->at(i), 0, down);
          for(int j = 0; j < d_sz; j++) {
            if(skip_map[down[j]]) {
              j = d_sz;
              ent_col->erase(ent_col->begin()+i);
            }
          }
        }
      }
    }
    else {
      apf::Downward down;
      for(int i = ent_col->size()-1; i > -1; i--) {
        int d_sz = e_lens.m->getDownward(ent_col->at(i), 0, down);
        for(int j = 0; j < d_sz; j++) {
          if(skip_map[down[j]]) {
            j = d_sz;
            ent_col->erase(ent_col->begin()+i);
          }
        }
      }
    }
  }
}

std::map<apf::MeshEntity*, bool> vd_edisc::calc_skip() {
  std::map<apf::MeshEntity*, bool> skip_dist;
  if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL and 
                                                  ext_0cell) {
    ext_shell* e_sh = f_calc.get_e_sh();
    if(trial_type == 1) {
      if(e_sh->chk_shell(trial_type, e_lens.new_cell1_id.at(0)-1) ) {
        shell sh_1c = e_sh->get_shell(trial_type, e_lens.new_cell1_id.at(0)-1);

        assert(e_sh->chk_shell(0, e_lens.new_cell0_id.at(0)-1) and 
               e_sh->chk_shell(0, e_lens.new_cell0_id.at(1)-1) );
        shell sh_0c1 = e_sh->get_shell(0, e_lens.new_cell0_id.at(0)-1);
        shell sh_0c2 = e_sh->get_shell(0, e_lens.new_cell0_id.at(1)-1);

        if(sh_0c1.dim < sh_1c.dim) {
          assert(sh_0c1.dim < sh_0c2.dim);
          skip_dist[e_lens.vert_ctr_new.at(0)] = true;
        }
        else if(sh_0c2.dim < sh_1c.dim) {
          skip_dist[e_lens.vert_ctr_new.at(1)] = true;
        }
        else
          skip_dist[e_lens.vert_sp_new.at(0)] = true;
      }
      else {
        if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(0)-1)) {
          if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(1)-1)) {
            shell sh_0c1 = e_sh->get_shell(0, e_lens.new_cell0_id.at(0)-1);
            shell sh_0c2 = e_sh->get_shell(0, e_lens.new_cell0_id.at(1)-1);
            if(sh_0c1.dim < sh_0c2.dim)
              skip_dist[e_lens.vert_ctr_new.at(0)] = true;
            else
              skip_dist[e_lens.vert_ctr_new.at(1)] = true;
          }
          else {
            skip_dist[e_lens.vert_ctr_new.at(0)] = true;
          }
        }
        else {
          assert(e_sh->chk_shell(0, e_lens.new_cell0_id.at(1)-1));
          skip_dist[e_lens.vert_ctr_new.at(1)] = true;
        }
      }

    }
    else {
      assert(trial_type == 2);
      bool found = false;
      int id = -1;
      int type = -1;
      shell temp;
      shell sh_2c;

      if(e_sh->chk_shell(trial_type, e_lens.new_cell2_id-1)) {
        sh_2c = e_sh->get_shell(trial_type, e_lens.new_cell2_id-1);

        for(int i = 0; i < e_lens.new_cell1_id.size(); i++) {
          if(e_sh->chk_shell(1, e_lens.new_cell1_id.at(i)-1 ) ) {
            shell sh_1c = e_sh->get_shell(1, e_lens.new_cell1_id.at(i)-1 );
            if(found) {
              assert(!(sh_1c.dim < temp.dim));
            }
            else {
              if(sh_1c.dim < sh_2c.dim) {
                found = true;
                temp = sh_1c;
                id = i;
                type = 1;
              }
            }
          }
        }

        for(int i = 0; i < e_lens.new_cell0_id.size(); i++) {
          if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(i)-1 ) ) {
            shell sh_0c = e_sh->get_shell(0, e_lens.new_cell0_id.at(i)-1 );
            if(found) {
              if(sh_0c.dim < temp.dim) {
                temp = sh_0c;
                id = i;
                type = 0;
              }
            }
            else {
              if(sh_0c.dim < sh_2c.dim) {
                found = true;
                temp = sh_0c;
                id = i;
                type = 0;
              }
            }
          }
        }
      }
      else {
        for(int i = 0; i < e_lens.new_cell1_id.size(); i++) {
          if(e_sh->chk_shell(1, e_lens.new_cell1_id.at(i)-1 ) ) {
            shell sh_1c = e_sh->get_shell(1, e_lens.new_cell1_id.at(i)-1 );
            if(found) {
              assert(!(sh_1c.dim < temp.dim));
            }
            else {
              found = true;
              temp = sh_1c;
              id = i;
              type = 1;
            }
          }
        }

        for(int i = 0; i < e_lens.new_cell0_id.size(); i++) {
          if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(i)-1 ) ) {
            shell sh_0c = e_sh->get_shell(0, e_lens.new_cell0_id.at(i)-1 );
            if(found) {
              if(sh_0c.dim < temp.dim) {
                temp = sh_0c;
                id = i;
                type = 0;
              }
            }
            else {
              found = true;
              temp = sh_0c;
              id = i;
              type = 0;
            }
          }
        }
      }

      if(found) {
        if(type == 0) {
          skip_dist[e_lens.vert_ctr_new.at(id)] = true;
        }
        else
          skip_dist[e_lens.vert_sp_new.at(id)] = true;
      }
    }
  }
  return skip_dist;
}

// Calculate the highest multiplier for velocity field that would invert any 
// tetrahedra. Position the new vertices such that nothing inverts. Used before
// the lenses are fully expanded.
void vd_edisc::find_pos_non_inv_init() {
  apf::Field* vel_field = e_lens.m->findField("velocity_field");
/*
  double t_min = -1;
  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr_new.at(i), drag_flag);
    apf::getVector(vel_field, e_lens.vert_ctr_new.at(i), 0, e_lens.vt.at(i));
  }
  for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.split_vt.at(0) = apf::Vector3(0,0,0);
  }
  t_min = f_calc.find_min_t2(e_lens.m, &e_lens.vert_ctr_new);

  if(t_min < -std::numeric_limits<double>::min())
    t_min = 0;
  t_min = t_min/rat_init;

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, e_lens.vt.at(i)*t_min 
                                                                    + pos_old);
  }
*/
  std::vector<apf::MeshEntity*> e_e(0);
  std::vector<apf::MeshEntity*> e_s(0);
  std::vector<apf::MeshEntity*> e_t(0);

  apf::Vector3 t_pos(0,0,0);
  apf::Vector3 v_pos(0,0,0);
  apf::Vector3 sn(0,0,0);

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr_new.at(i), drag_flag);
    apf::getVector(vel_field, e_lens.vert_ctr_new.at(i), 0, e_lens.vt.at(i));
    e_lens.vt.at(i) = norm_0(e_lens.vt.at(i));
  }
  for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.split_vt.at(0) = apf::Vector3(0,0,0);
  }

  vd_set_up(e_lens.m, &e_lens.vert_sp_new, &e_e);
  vd_set_up(e_lens.m, &e_e, &e_s);

  std::map<apf::MeshEntity*, bool> tri_lens{};
  for(int i = 0; i < e_s.size(); i++) {
    tri_lens[e_s.at(i)] = true;
  }

  double t_min = -1;
  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    //init_ent_set(&e_v, vert->at(j));
    vd_set_up(e_lens.m, e_lens.vert_ctr_new.at(i), &e_e);
    vd_set_up(e_lens.m, &e_e, &e_s);
    vd_set_up(e_lens.m, &e_s, &e_t);

    //std::cout << "v_pos " << v_pos << " dn " << dn << std::endl;
    for(int j = 0; j < e_t.size(); j++) {
      apf::Downward d_t;
      apf::Downward d_v;
      e_lens.m->getDownward(e_t.at(j), 2, d_t);
      e_lens.m->getDownward(e_t.at(j), 0, d_v);
      int v_id = findIn(d_v, 4, e_lens.vert_ctr_new.at(i));
      assert(v_id > -1);

      apf::MeshEntity* t_curr = d_t[lookup_tet_surf_x[v_id]];
      if(!tri_lens[t_curr]) {
        sn = vd_area_out(e_lens.m, t_curr, 0);

        t_pos = getLinearCentroid(e_lens.m, t_curr);
        t_pos = t_pos - pos_old;

        assert(t_pos.getLength() > std::numeric_limits<double>::min());
        sn = norm_0(sn);

        if(sn*t_pos < - std::numeric_limits<double>::min()
            and std::fabs(sn*t_pos) > std::numeric_limits<double>::min())
          sn = sn*(-1);

        double dist = t_pos*sn;
        double v = e_lens.vt.at(i)*sn;
        if(v > std::numeric_limits<double>::min() 
            and std::fabs(v) > std::numeric_limits<double>::min()) {
          double t = dist/v;
          t = fabs(t);
          if (t_min < 0. or t < t_min) {
            t_min = t;
          }
        }
      }
    }
  }
  if(t_min < -std::numeric_limits<double>::min())
    t_min = 0;
  t_min = t_min/rat_init;

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, e_lens.vt.at(i)*t_min 
                                                                    + pos_old);
  }
  t_min = t_min*rat_init;
  if(trial_type == 1 or (trial_type == 2 and e_lens.vert_sp_new.size() > 2)) {
    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      int v1 = e_lens.slices.at(2*i)-1;
      int v2 = e_lens.slices.at(2*i+1)-1;
      e_lens.m->getPoint(e_lens.vert_ctr_new.at(v1), 0, t_pos);
      e_lens.m->getPoint(e_lens.vert_ctr_new.at(v2), 0, v_pos);
      e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, t_pos/2
                                                          +v_pos/2);
    }
  }
  else {
    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_sp_new.at(i), drag_flag);
      apf::getVector(vel_field, e_lens.vert_sp_new.at(i), 0, 
                                                e_lens.split_vt.at(i));
    }

    vd_set_up(e_lens.m, e_lens.v_2c, &e_e);
    vd_set_up(e_lens.m, &e_e, &e_s);

    std::map<apf::MeshEntity*, bool> tri_void{};
    for(int i = 0; i < e_s.size(); i++) {
      tri_void[e_s.at(i)] = true;
    }

    double t_min = -1;
    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      //init_ent_set(&e_v, vert->at(j));
      vd_set_up(e_lens.m, e_lens.vert_sp_new.at(i), &e_e);
      vd_set_up(e_lens.m, &e_e, &e_s);
      vd_set_up(e_lens.m, &e_s, &e_t);

      //std::cout << "v_pos " << v_pos << " dn " << dn << std::endl;
      for(int j = 0; j < e_t.size(); j++) {
        apf::Downward d_t;
        apf::Downward d_v;
        e_lens.m->getDownward(e_t.at(j), 2, d_t);
        e_lens.m->getDownward(e_t.at(j), 0, d_v);
        int v_id = findIn(d_v, 4, e_lens.vert_sp_new.at(i));
        assert(v_id > -1);

        apf::MeshEntity* t_curr = d_t[lookup_tet_surf_x[v_id]];
        if(!tri_void[t_curr]) {
          sn = vd_area_out(e_lens.m, t_curr, 0);

          t_pos = getLinearCentroid(e_lens.m, t_curr);
          t_pos = t_pos - pos_old;

          assert(t_pos.getLength() > std::numeric_limits<double>::min());
          sn = norm_0(sn);

          if(sn*t_pos < - std::numeric_limits<double>::min()
              and std::fabs(sn*t_pos) > std::numeric_limits<double>::min())
            sn = sn*(-1);

          double dist = t_pos*sn;
          double v = e_lens.split_vt.at(i)*sn;
          if(v > std::numeric_limits<double>::min() 
              and std::fabs(v) > std::numeric_limits<double>::min()) {
            double t = dist/v;
            t = fabs(t);
            if (t_min < 0. or t < t_min) {
              t_min = t;
            }
          }
        }
      }
    }
    if(t_min < -std::numeric_limits<double>::min())
      t_min = 0;
    t_min = t_min/rat_init;

    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, e_lens.split_vt.at(i)*t_min 
                                                                      + pos_old);
    }

  }
  if(trial_type == 2) {
    apf::Vector3 v_2c_pos(0,0,0);
    v_2c_pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
              vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;

    e_lens.m->setPoint(e_lens.v_2c, 0, v_2c_pos);
  }
}


// Calculate the highest multiplier for vector field that would invert any 
// tetrahedra. Position the new vertices such that nothing inverts.
void vd_edisc::find_pos_non_inv() {
  std::vector<apf::MeshEntity*> e_e(0);
  std::vector<apf::MeshEntity*> e_s(0);
  std::vector<apf::MeshEntity*> e_t(0);

  apf::Vector3 t_pos(0,0,0);
  apf::Vector3 v_pos(0,0,0);
  apf::Vector3 sn(0,0,0);

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, v_pos);
    e_lens.vt.at(i) = v_pos - pos_old;
    e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_old);
  }
  for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, v_pos);
    e_lens.split_vt.at(i) = v_pos - pos_old;
    e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_old);
  }

  vd_set_up(e_lens.m, &e_lens.vert_sp_new, &e_e);
  vd_set_up(e_lens.m, &e_e, &e_s);

  std::map<apf::MeshEntity*, bool> tri_lens{};
  for(int i = 0; i < e_s.size(); i++) {
    tri_lens[e_s.at(i)] = true;
  }

  double t_min = -1;
  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    //init_ent_set(&e_v, vert->at(j));
    vd_set_up(e_lens.m, e_lens.vert_ctr_new.at(i), &e_e);
    vd_set_up(e_lens.m, &e_e, &e_s);
    vd_set_up(e_lens.m, &e_s, &e_t);

    //std::cout << "v_pos " << v_pos << " dn " << dn << std::endl;
    for(int j = 0; j < e_t.size(); j++) {
      apf::Downward d_t;
      apf::Downward d_v;
      e_lens.m->getDownward(e_t.at(j), 2, d_t);
      e_lens.m->getDownward(e_t.at(j), 0, d_v);
      int v_id = findIn(d_v, 4, e_lens.vert_ctr_new.at(i));
      assert(v_id > -1);

      apf::MeshEntity* t_curr = d_t[lookup_tet_surf_x[v_id]];
      if(!tri_lens[t_curr]) {
        sn = vd_area_out(e_lens.m, t_curr, 0);

        t_pos = getLinearCentroid(e_lens.m, t_curr);
        t_pos = t_pos - pos_old;

        assert(t_pos.getLength() > std::numeric_limits<double>::min());
        sn = norm_0(sn);

        if(sn*t_pos < - std::numeric_limits<double>::min()
            and std::fabs(sn*t_pos) > std::numeric_limits<double>::min())
          sn = sn*(-1);

        double dist = t_pos*sn;
        double v = e_lens.vt.at(i)*sn;
        if(v > std::numeric_limits<double>::min() 
            and std::fabs(v) > std::numeric_limits<double>::min()) {
          double t = dist/v;
          t = fabs(t);
          if (t_min < 0 or t < t_min) {
            t_min = t;
          }
        }
      }
    }
  }
  if(t_min < -std::numeric_limits<double>::min())
    t_min = 0;
  t_min = t_min/rat_init;

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, e_lens.vt.at(i)*t_min 
                                                                    + pos_old);
  }
  if(trial_type == 2 and e_lens.vert_sp_new.size() > 2) {
    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      int v1 = e_lens.slices.at(2*i)-1;
      int v2 = e_lens.slices.at(2*i+1)-1;
      e_lens.m->getPoint(e_lens.vert_ctr_new.at(v1), 0, t_pos);
      e_lens.m->getPoint(e_lens.vert_ctr_new.at(v2), 0, v_pos);
      e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, t_pos/2
                                                          +v_pos/2);
    }
  }
  else {
    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, e_lens.split_vt.at(i)*t_min 
                                                                      + pos_old);
    }
  }
  if(trial_type == 2) {
    apf::Vector3 v_2c_pos(0,0,0);
    v_2c_pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
              vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;

    e_lens.m->setPoint(e_lens.v_2c, 0, v_2c_pos);
  }
}

double vd_edisc::calc_dt_ext_dir(std::vector<apf::Vector3> &v_ctr, 
                               std::vector<apf::Vector3> &v_sp, 
                               std::vector<apf::Vector3> &r_ctr, 
                               std::vector<apf::Vector3> &r_sp,
                               std::vector<apf::Vector3> &v_r_ctr, 
                               std::vector<apf::Vector3> &v_r_ctr_old, 
                               std::vector<apf::Vector3> &v_t_ctr, 
                               std::vector<apf::Vector3> &v_t_ctr_old, 
                               std::vector<apf::Vector3> &v_r_sp, 
                               std::vector<apf::Vector3> &v_r_sp_old, 
                               std::vector<apf::Vector3> &v_t_sp, 
                               std::vector<apf::Vector3> &v_t_sp_old,
                               double &r_min, double &r_max,
                               double &r_th_min, double &r_th_max, 
                               double &ang_th, double rat_t) {
  apf::Vector3 pos_curr(0,0,0);
  apf::Vector3 r_hat(0,0,0);
  apf::Vector3 t_hat(0,0,0);
  apf::Vector3 v_temp(0,0,0);

  apf::Field* vel_field = e_lens.m->findField("velocity_field");

  double v_r_curr = 0;
  double v_t_curr = 0;
  double r_curr = 0;
  double v_r_max = 0;

  double dt_th;
  //double dt_inv;
  double dt_t = -1;
  double dt_r = -1;
  double dt_t_curr = -1;
  double dt_r_curr = -1;

  e_lens.pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
                vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;

  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    apf::getVector(vel_field, e_lens.vert_ctr_new.at(i), 0, v_temp);
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, pos_curr);
    v_ctr.at(i) = v_temp;
    r_ctr.at(i) = pos_curr - e_lens.pos;

    r_hat = norm_0(r_ctr.at(i));
    r_curr = r_ctr.at(i).getLength();

    v_r_ctr.at(i) = r_hat*(v_temp*r_hat);
    v_t_ctr.at(i) = v_temp - v_r_ctr.at(i);
    t_hat = norm_0(v_t_ctr.at(i));

    v_r_curr = (v_r_ctr.at(i) + v_r_ctr_old.at(i))*r_hat/2;
    v_t_curr = (v_t_ctr.at(i) + v_t_ctr_old.at(i))*t_hat/2;

    if(r_curr < r_min)
      r_min = r_curr;
    if(r_curr > r_max)
      r_max = r_curr;

    if(std::fabs(v_r_curr) > std::numeric_limits<double>::epsilon()) {
      if(v_r_curr > 0)
        dt_r_curr = std::fabs((r_th_max*1.2 - r_curr)/v_r_curr/rat_t);
      else
        dt_r_curr = std::fabs((r_curr - r_th_min/2)/v_r_curr/rat_t);
    }
    else
      dt_r_curr = -1;

    if(std::fabs(v_t_curr) > std::numeric_limits<double>::epsilon())
      dt_t_curr = std::fabs(ang_th*r_curr/v_t_curr/rat_t);
    else
      dt_t_curr = -1;

    if(dt_r_curr > 0 and (dt_r_curr < dt_r or dt_r < 0))
      dt_r = dt_r_curr;
    if(dt_r_curr > 0 and (dt_t_curr < dt_t or dt_t < 0))
      dt_t = dt_t_curr;
  }

  if(trial_type == 2) {
    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      apf::getVector(vel_field, e_lens.vert_sp_new.at(i), 0, v_temp);

      if(ext_new[e_lens.vert_sp_new.at(i)]) {

        int v1 = e_lens.slices.at(2*i)-1;
        int v2 = e_lens.slices.at(2*i+1)-1;

        apf::Vector3 vel_1(0,0,0);
        apf::getVector(vel_field, e_lens.vert_ctr_new.at(v1), 0, vel_1);
        apf::getVector(vel_field, e_lens.vert_ctr_new.at(v2), 0, v_temp);
        v_temp = (v_temp + vel_1)/2;
      }
      //else
      //  v_temp = v_temp - v_avg;
      //v_temp = v_temp - v_avg;
      apf::setVector(vel_field, e_lens.vert_sp_new.at(i), 0, v_temp);
      v_sp.at(i) = v_temp;

      e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, pos_curr);
      r_sp.at(i) = pos_curr - e_lens.pos;

      r_hat = norm_0(r_sp.at(i));
      r_curr = r_sp.at(i).getLength();

      v_r_sp.at(i) = r_hat*(v_temp*r_hat);
      v_t_sp.at(i) = v_temp - v_r_sp.at(i);
      t_hat = norm_0(v_t_sp.at(i));

      v_r_curr = (v_r_sp.at(i) + v_r_sp_old.at(i))*r_hat/2;
      v_t_curr = (v_t_sp.at(i) + v_t_sp_old.at(i))*t_hat/2;

      if(r_curr < r_min)
        r_min = r_curr;
      if(r_curr > r_max)
        r_max = r_curr;

      if(std::fabs(v_r_curr) > std::numeric_limits<double>::epsilon()) {
        if(v_r_curr > 0)
          dt_r_curr = std::fabs((r_th_max*1.2 - r_curr)/v_r_curr/rat_t);
        else
          dt_r_curr = std::fabs((r_curr - r_th_min/2)/v_r_curr/rat_t);
      }
      else
        dt_r_curr = -1;

      if(std::fabs(v_t_curr) > std::numeric_limits<double>::epsilon())
        dt_t_curr = std::fabs(ang_th*r_curr/v_t_curr/rat_t);
      else
        dt_t_curr = -1;

      if(dt_r_curr > 0 and (dt_r_curr < dt_r or dt_r < 0))
        dt_r = dt_r_curr;
      if(dt_r_curr > 0 and (dt_t_curr < dt_t or dt_t < 0))
        dt_t = dt_t_curr;
    }
  }

  if(dt_r < 0) {
    if(dt_t < 0)
      dt_th = -1;
    else
      dt_th = dt_t;
  }
  else {
    if(dt_t < 0)
      dt_th = dt_r;
    else
      dt_th = std::min(dt_t, dt_r);
  }
  return dt_th;
}


void vd_edisc::move_vert_ext(std::vector<apf::Vector3> &v_ctr, 
                               std::vector<apf::Vector3> &v_ctr_old, 
                               std::vector<apf::Vector3> &v_sp,
                               std::vector<apf::Vector3> &v_sp_old, 
                               double &d_max) {
  apf::Vector3 pos_curr(0,0,0);
  apf::Vector3 v_temp(0,0,0);
  double r_curr = 0;
  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, pos_curr);

    v_temp = (v_ctr.at(i)+v_ctr_old.at(i))/2;
    //std::cout << "x_ctr("<< i << ")= " << pos_curr << std::endl;
    //std::cout << "v_ctr("<< i << ")= " << v_temp << std::endl;
    r_curr = v_temp.getLength()*f_calc.vdparam.dt;
    //std::cout << "r("<< i << ")= " << r_curr << std::endl;
    if(d_max < r_curr)
      d_max = r_curr;
    pos_curr = pos_curr + v_temp*f_calc.vdparam.dt;

    e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_curr);
  }

  if(trial_type == 2) {
    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, pos_curr);

      //pos_curr = pos_curr + (v_sp.at(i)+v_sp_old.at(i))/2*f_calc.vdparam.dt;
      v_temp = (v_sp.at(i)+v_sp_old.at(i))/2;
      //std::cout << "v_sp("<< i << ")= " << v_temp << std::endl;

      r_curr = v_temp.getLength()*f_calc.vdparam.dt;
      //std::cout << "r("<< i << ")= " << r_curr << std::endl;
      if(d_max < r_curr)
        d_max = r_curr;

      if(e_lens.vert_sp_new.size() > 2) {
        int v1 = e_lens.slices.at(2*i)-1;
        int v2 = e_lens.slices.at(2*i+1)-1;
        e_lens.m->getPoint(e_lens.vert_ctr_new.at(v1), 0, pos_curr);
        e_lens.m->getPoint(e_lens.vert_ctr_new.at(v2), 0, v_temp);
        e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_curr/2
                                                            +v_temp/2);
      }
      else {
        pos_curr = pos_curr + v_temp*f_calc.vdparam.dt;
        e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_curr);
      }
      //r_hat = norm_0(r_sp.at(i));
      //r_curr = r_sp.at(i).getLength();
      //std::cout << "r_curr " << r_curr << std::endl;
    }

    apf::Vector3 v_2c_pos(0,0,0);
    v_2c_pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
              vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;

    e_lens.m->setPoint(e_lens.v_2c, 0, v_2c_pos);

  }
  else {
    assert(trial_type == 1);
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(0), 0, v_temp);
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(1), 0, pos_curr);
    e_lens.m->setPoint(e_lens.vert_sp_new.at(0), 0, pos_curr/2
                                                        +v_temp/2);       
  }
}

// Calculate the positions of the new 0cell vertices using the Mason algorithm,
// by iterating until the direction of motions of vertices converge.
void vd_edisc::find_ext_dir() {
  assert(len_sp > std::numeric_limits<double>::epsilon());
  // TODO Dirty way of overcoming the jumping over the minimum energy position
  // due to zero mean velocity expansion. A better approach should be by 
  // setting the objective as finding the minimum energy positions for the 
  // given insertion and calculating a mean radius of expansion to use in 
  // calculating the energy rate. This would be necessary for anisotropic 
  // properties.
  double rho = len_sp/rho_rat/fudge_factor;
  double disp_th = rho/rho_rat/1.2/fudge_factor;

  double rad_tol = disp_th/6/fudge_factor;

  double r_th_min = rho/100/fudge_factor;
  double r_th_max = 2*rho/fudge_factor;

  // Displacement threshold. If the maximum motion is less than this amount
  // stop iteration.
  double d_th = disp_th/2400/fudge_factor;

  double r_min = rho;
  double r_max = 0;

  double ang_th;
  if(trial_type == 1)
    ang_th = 2*pi/18;
  else
    ang_th = 2*pi/e_lens.vert_ctr_new.size()/6;

  //double dir_tol = rho*ang_th;

  double max_ang = 6.14/f_calc.vdparam.dt;

  double v_r_curr = 0;
  double v_t_curr = 0;
  double r_curr = 0;

  double v_r_max = 0;

  int counter_sz = 50;
  double rat_t = 3;

  std::cout << "rho " << rho << " len_sp " << len_sp << std::endl;
  std::cout << "ang_th " << ang_th << " max_ang " << max_ang << std::endl;

  int iter_sz = 4;

  double dt_glob = f_calc.vdparam.dt;

  apf::Field* vel_field = e_lens.m->findField("velocity_field");

  apf::Vector3 pos_curr(0,0,0);
  apf::Vector3 r_hat(0,0,0);
  apf::Vector3 t_hat(0,0,0);
  apf::Vector3 v_temp(0,0,0);

  std::vector<apf::Vector3> v_r_ctr(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_r_sp(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> v_r_ctr_old(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_r_sp_old(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> v_t_ctr(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_t_sp(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> v_t_ctr_old(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_t_sp_old(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> v_ctr_old(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_sp_old(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> v_ctr(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> v_sp(0, apf::Vector3(0,0,0));

  std::vector<apf::Vector3> r_ctr(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> r_sp(0, apf::Vector3(0,0,0));

  int v_sz = e_lens.vert_ctr_new.size()+e_lens.vert_sp_new.size();
  int v_ctr_sz = e_lens.vert_ctr_new.size();
  int v_sp_sz = e_lens.vert_sp_new.size();

  r_ctr.resize(e_lens.vert_ctr_new.size());
  r_sp.resize(e_lens.vert_sp_new.size());

  v_r_ctr.resize(e_lens.vert_ctr_new.size());
  v_t_ctr.resize(e_lens.vert_ctr_new.size());

  v_r_sp.resize(e_lens.vert_sp_new.size());
  v_t_sp.resize(e_lens.vert_sp_new.size());

  v_r_ctr_old.resize(e_lens.vert_ctr_new.size());
  v_r_sp_old.resize(e_lens.vert_sp_new.size());

  v_t_ctr_old.resize(e_lens.vert_ctr_new.size());
  v_t_sp_old.resize(e_lens.vert_sp_new.size());

  v_ctr.resize(e_lens.vert_ctr_new.size());
  v_sp.resize(e_lens.vert_sp_new.size());

  v_ctr_old.resize(e_lens.vert_ctr_new.size());
  v_sp_old.resize(e_lens.vert_sp_new.size());

  //The exterior vertex will be positioned at the old location, so the closing
  // check will be positive. Skip that vertex.
  //std::map<apf::MeshEntity*, bool> skip_dist{};
  //skip_dist = calc_skip();

  double z[3] = {0,0,0};
  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    v_r_ctr.at(i).fromArray(z);
    v_t_ctr.at(i).fromArray(z);
    v_r_ctr_old.at(i).fromArray(z);
    v_t_ctr_old.at(i).fromArray(z);
    v_ctr_old.at(i).fromArray(z);
    v_ctr.at(i).fromArray(z);
  }

  for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    v_t_sp.at(i).fromArray(z);
    v_r_sp.at(i).fromArray(z);
    v_t_sp_old.at(i).fromArray(z);
    v_r_sp_old.at(i).fromArray(z);
    v_sp_old.at(i).fromArray(z);
    v_sp.at(i).fromArray(z);
  }

  double dt_th;
  //double dt_inv;
  double dt_t = -1;
  double dt_r = -1;

  proj_around_ctr(rho);

  std::tie(r_min, r_max) = vd_find_sph_bound(); 

  if(r_min < r_th_min) {
    if (r_min < std::numeric_limits<double>::epsilon() ) {
      proj_around_ctr(rho);
    }
    else // Project the minimum radius vertex to twice the min radius threshold.
      proj_around_ctr(rho*2*r_th_min/r_min);
  }


  double e_i_old = calc_energy_lens();
  double e_f_old = e_i_old;

  double e_i_new = 2*e_i_old;
  double e_f_new = 2*e_f_old;

  // The iteration will end if the percentile change in initial and final 
  // energies are both below. 
  // If there is some singularity, energy will be zero
  // TODO with non-positive energy terms, the energy can actually become positive.
  double e_th;
  if(std::fabs(e_i_old) > std::numeric_limits<double>::epsilon())
    e_th = 1/2.1;
//    e_th = std::fabs((e_i_new-e_i_old)/e_i_old)/2.1;
  else
    e_th = 1;

  int mov_tag = 0;
  bool closing = false;

  int iter = 0;

  while(std::fabs(e_i_old) > std::numeric_limits<double>::epsilon() and
        std::fabs(e_f_old) > std::numeric_limits<double>::epsilon() and
        std::fabs((e_i_new-e_i_old)/e_i_old) > e_th and 
        std::fabs((e_f_new-e_f_old)/e_f_old) > e_th and
        iter < iter_sz) {

    t_total = 0;

    double dt_t_curr = -1;
    double dt_r_curr = -1;

    //apf::Vector3 v_avg(0,0,0);
    // The components of the average velocities not removed from the vertices.
    int ext_count = 0;

    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      //vd_upd_vel_field_mason(e_lens.m, e_lens.vert_ctr_new.at(i), &f_calc);
      f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr_new.at(i), drag_flag);
      apf::getVector(vel_field, e_lens.vert_ctr_new.at(i), 0, e_lens.vt.at(i));
    }
    if(trial_type == 2) {
      for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
        //vd_upd_vel_field_mason(e_lens.m, e_lens.vert_sp_new.at(i), &f_calc);
        f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_sp_new.at(i), drag_flag);
      }
    }

    dt_th = calc_dt_ext_dir(v_ctr, v_sp, r_ctr, r_sp, v_r_ctr, v_r_ctr_old, 
            v_t_ctr, v_t_ctr_old, v_r_sp, v_r_sp_old, v_t_sp, v_t_sp_old,
                    r_min, r_max, r_th_min, r_th_max, ang_th, rat_t);
/*
    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      apf::getVector(vel_field, e_lens.vert_ctr_new.at(i), 0, v_temp);
      //std::cout << "v_ctr("<< i << ")= " << v_temp << std::endl;
      e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, pos_curr);
      r_ctr.at(i) = pos_curr-e_lens.pos;
      r_hat = norm_0(r_ctr.at(i));

      r_curr = r_ctr.at(i).getLength();
      //if(!skip_dist[e_lens.vert_ctr_new.at(i)] and r_curr < r_th_min) {
        //std::cout << "r_curr < r_th_min " << r_curr 
        //          << " " << r_th_min << std::endl;
      //  closing = true;
      //}
      std::tie(r_min, r_max) = vd_find_sph_bound(); 
      if(r_min < r_th_min) {
        //std::cout << "r_curr < r_th_min " << r_curr 
        //          << " " << r_th_min << std::endl;
        closing = true;
      }

      v_r_ctr.at(i) = r_hat*(v_temp*r_hat);
      v_t_ctr.at(i) = v_temp - v_r_ctr.at(i);
      t_hat = norm_0(v_t_ctr.at(i));

      v_r_curr = (v_r_ctr.at(i) + v_r_ctr_old.at(i))*r_hat/2;
      v_t_curr = (v_t_ctr.at(i) + v_t_ctr_old.at(i))*t_hat/2;

      if(std::fabs(v_r_curr) > std::numeric_limits<double>::epsilon()) {
        if(v_r_curr > 0)
          dt_r_curr = std::fabs((r_th_max*1.2 - r_curr)/v_r_curr/rat_t);
        else
          dt_r_curr = std::fabs((r_curr - r_th_min/2)/v_r_curr/rat_t);
      }
      else
        dt_r_curr = -1;

      if(std::fabs(v_t_curr) > std::numeric_limits<double>::epsilon())
        dt_t_curr = std::fabs(ang_th*r_curr/v_t_curr/rat_t);
      else
        dt_t_curr = -1;

      if(dt_r_curr > 0 and (dt_r_curr < dt_r or dt_r < 0))
        dt_r = dt_r_curr;
      if(dt_r_curr > 0 and (dt_t_curr < dt_t or dt_t < 0))
        dt_t = dt_t_curr;
    }

    //and e_lens.vert_sp_new.size() == 2
    if(trial_type == 2) {
      for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
        //vd_upd_vel_field_mason(e_lens.m, e_lens.vert_sp_new.at(i), &f_calc);
        f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_sp_new.at(i), drag_flag);
        apf::getVector(vel_field, e_lens.vert_sp_new.at(i), 0, v_temp);

        if(e_lens.vert_sp_new.size() > 2 and
           ext_new[e_lens.vert_sp_new.at(i)]) {
          int v1 = e_lens.slices.at(2*i)-1;
          int v2 = e_lens.slices.at(2*i+1)-1;

          apf::Vector3 vel_1;
          apf::getVector(vel_field, e_lens.vert_ctr_new.at(v1), 0, vel_1);
          apf::getVector(vel_field, e_lens.vert_ctr_new.at(v2), 0, v_temp);
          v_temp = (v_temp + vel_1)/2;
        }
        //else
        //  v_temp = v_temp - v_avg;

        //v_temp = v_temp - v_avg;
        apf::setVector(vel_field, e_lens.vert_sp_new.at(i), 0, v_temp);
        //std::cout << "v_sp("<< i << ")= " << v_temp << std::endl;
        e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, pos_curr);
        r_sp.at(i) = pos_curr-e_lens.pos;
        r_hat = norm_0(r_sp.at(i));

        r_curr = r_sp.at(i).getLength();
        //if(!skip_dist[e_lens.vert_sp_new.at(i)] and r_curr < r_th_min) {
          //std::cout << "r_curr < r_th_min " << r_curr 
          //          << " " << r_th_min << std::endl;
        //  closing = true;
        //}
        std::tie(r_min, r_max) = vd_find_sph_bound(); 
        if(r_min < r_th_min) {
          //std::cout << "r_curr < r_th_min " << r_curr 
          //          << " " << r_th_min << std::endl;
          closing = true;
        }

        v_r_sp.at(i) = r_hat*(v_temp*r_hat);
        v_t_sp.at(i) = v_temp - v_r_sp.at(i);
        t_hat = norm_0(v_t_sp.at(i));

        v_r_curr = (v_r_sp.at(i) + v_r_sp_old.at(i))*r_hat/2;
        v_t_curr = (v_t_sp.at(i) + v_t_sp_old.at(i))*t_hat/2;

        if(std::fabs(v_r_curr) > std::numeric_limits<double>::epsilon()) {
          if(v_r_curr > 0)
            dt_r_curr = std::fabs((r_th_max*1.2 - r_curr)/v_r_curr/rat_t);
          else
            dt_r_curr = std::fabs((r_curr - r_th_min/2)/v_r_curr/rat_t);
        }
        else
          dt_r_curr = -1;

        if(std::fabs(v_t_curr) > std::numeric_limits<double>::epsilon())
          dt_t_curr = std::fabs(ang_th*r_curr/v_t_curr/rat_t);
        else
          dt_t_curr = -1;

        if(dt_r_curr > 0 and (dt_r_curr < dt_r or dt_r < 0))
          dt_r = dt_r_curr;
        if(dt_r_curr > 0 and (dt_t_curr < dt_t or dt_t < 0))
          dt_t = dt_t_curr;
      }

    }

    r_min = rho;
    r_max = 0;

    if(dt_r < 0) {
      if(dt_t < 0)
        dt_th = -1;
      else
        dt_th = dt_t;
    }
    else {
      if(dt_t < 0)
        dt_th = dt_r;
      else
        dt_th = std::min(dt_t, dt_r);
    }
*/
    if(dt_th > std::numeric_limits<double>::epsilon())
      f_calc.vdparam.adj_dt(dt_th);
    else {
      f_calc.vdparam.adj_dt(dt_glob);
      std::cout << "Using dt_glob" << std::endl;
    }

    std::cout << "r_th_min " << r_th_min << " r_th_max " << r_th_max 
              << std::endl;
    std::cout << "r_min " << r_min << std::endl;

    double d_max = 2*d_th;
    std::cout << "d_max " << d_max << " d_th " << d_th 
              << std::endl;
    int counter = 0;
    // The zero mean velocity expansion may cause the energy to increase.
    // Expand as long as the energy is decreasing.
    double del_e = -1;
    while(r_min > r_th_min and r_max < r_th_max and d_max > d_th 
                           and counter < counter_sz
                           and del_e < 0) {
      double en_curr = calc_energy_lens();
      //while(r_min > r_th_min and r_max < r_th_max) {
      //std::cout << "rho " << rho << std::endl;
      // cancel the average velocity of the central vertices, so that biased 
      // initial conditions doesn't lead to migration away from the old vertex 
      // position.
      dt_r = -1;
      dt_t = -1;

      for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, pos_curr);

        //vd_upd_vel_field_mason(e_lens.m, e_lens.vert_ctr_new.at(i), &f_calc);
        f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr_new.at(i), drag_flag);
        apf::getVector(vel_field, e_lens.vert_ctr_new.at(i), 0, 
                                                            e_lens.vt.at(i));
      }
      if(trial_type == 2) {
        for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
          //vd_upd_vel_field_mason(e_lens.m, e_lens.vert_sp_new.at(i), &f_calc);

          f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_sp_new.at(i), drag_flag);
        }
      }
      dt_th = calc_dt_ext_dir(v_ctr, v_sp, r_ctr, r_sp, v_r_ctr, v_r_ctr_old, 
            v_t_ctr, v_t_ctr_old, v_r_sp, v_r_sp_old, v_t_sp, v_t_sp_old,
                    r_min, r_max, r_th_min, r_th_max, ang_th, rat_t);
/*
      for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        apf::getVector(vel_field, e_lens.vert_ctr_new.at(i), 0, v_temp);
        r_ctr.at(i) = pos_curr-e_lens.pos;
        v_ctr.at(i) = v_temp;

        r_hat = norm_0(r_ctr.at(i));
        r_curr = r_ctr.at(i).getLength();

        v_r_ctr.at(i) = r_hat*(v_temp*r_hat);
        v_t_ctr.at(i) = v_temp - v_r_ctr.at(i);
        t_hat = norm_0(v_t_ctr.at(i));

        v_r_curr = (v_r_ctr.at(i) + v_r_ctr_old.at(i))*r_hat/2;
        v_t_curr = (v_t_ctr.at(i) + v_t_ctr_old.at(i))*t_hat/2;

        if(r_curr < r_min)
          r_min = r_curr;
        if(r_curr > r_max)
          r_max = r_curr;

        if(std::fabs(v_r_curr) > std::numeric_limits<double>::epsilon()) {
          if(v_r_curr > 0)
            dt_r_curr = std::fabs((r_th_max*1.2 - r_curr)/v_r_curr/rat_t);
          else
            dt_r_curr = std::fabs((r_curr - r_th_min/2)/v_r_curr/rat_t);
        }
        else
          dt_r_curr = -1;

        if(std::fabs(v_t_curr) > std::numeric_limits<double>::epsilon())
          dt_t_curr = std::fabs(ang_th*r_curr/v_t_curr/rat_t);
        else
          dt_t_curr = -1;

        if(dt_r_curr > 0 and (dt_r_curr < dt_r or dt_r < 0))
          dt_r = dt_r_curr;
        if(dt_r_curr > 0 and (dt_t_curr < dt_t or dt_t < 0))
          dt_t = dt_t_curr;
      }

      if(trial_type == 2) {
        for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
          apf::getVector(vel_field, e_lens.vert_sp_new.at(i), 0, v_temp);

          if(ext_new[e_lens.vert_sp_new.at(i)]) {

            int v1 = e_lens.slices.at(2*i)-1;
            int v2 = e_lens.slices.at(2*i+1)-1;

            apf::Vector3 vel_1(0,0,0);
            apf::getVector(vel_field, e_lens.vert_ctr_new.at(v1), 0, vel_1);
            apf::getVector(vel_field, e_lens.vert_ctr_new.at(v2), 0, v_temp);
            v_temp = (v_temp + vel_1)/2;
          }
          //else
          //  v_temp = v_temp - v_avg;
          //v_temp = v_temp - v_avg;
          apf::setVector(vel_field, e_lens.vert_sp_new.at(i), 0, v_temp);
          v_sp.at(i) = v_temp;

          e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, pos_curr);
          r_sp.at(i) = pos_curr-e_lens.pos;

          r_hat = norm_0(r_sp.at(i));
          r_curr = r_sp.at(i).getLength();

          v_r_sp.at(i) = r_hat*(v_temp*r_hat);
          v_t_sp.at(i) = v_temp - v_r_sp.at(i);
          t_hat = norm_0(v_t_sp.at(i));

          v_r_curr = (v_r_sp.at(i) + v_r_sp_old.at(i))*r_hat/2;
          v_t_curr = (v_t_sp.at(i) + v_t_sp_old.at(i))*t_hat/2;

          if(r_curr < r_min)
            r_min = r_curr;
          if(r_curr > r_max)
            r_max = r_curr;

          if(std::fabs(v_r_curr) > std::numeric_limits<double>::epsilon()) {
            if(v_r_curr > 0)
              dt_r_curr = std::fabs((r_th_max*1.2 - r_curr)/v_r_curr/rat_t);
            else
              dt_r_curr = std::fabs((r_curr - r_th_min/2)/v_r_curr/rat_t);
          }
          else
            dt_r_curr = -1;

          if(std::fabs(v_t_curr) > std::numeric_limits<double>::epsilon())
            dt_t_curr = std::fabs(ang_th*r_curr/v_t_curr/rat_t);
          else
            dt_t_curr = -1;

          if(dt_r_curr > 0 and (dt_r_curr < dt_r or dt_r < 0))
            dt_r = dt_r_curr;
          if(dt_r_curr > 0 and (dt_t_curr < dt_t or dt_t < 0))
            dt_t = dt_t_curr;
        }
      }

      if(dt_r < 0) {
        if(dt_t < 0)
          dt_th = -1;
        else
          dt_th = dt_t;
      }
      else {
        if(dt_t < 0)
          dt_th = dt_r;
        else
          dt_th = std::min(dt_t, dt_r);
      }
*/
      if(dt_th > std::numeric_limits<double>::epsilon())
        f_calc.vdparam.adj_dt(dt_th);
      else {
        f_calc.vdparam.adj_dt(dt_glob);
        std::cout << "Using dt_glob" << std::endl;
      }
      t_total = t_total + f_calc.vdparam.dt;
      move_vert_ext(v_ctr, v_ctr_old, v_sp, v_sp_old, d_max);

      del_e = calc_energy_lens() - en_curr;
      // If energy is larger, revert changes and exit.
      if(del_e > 0) {
        for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
          e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, pos_curr);

          v_temp = (v_ctr.at(i)+v_ctr_old.at(i))/2;
          pos_curr = pos_curr - v_temp*f_calc.vdparam.dt;

          e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_curr);
        }

        if(trial_type == 2) {
          for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
            e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, pos_curr);

            v_temp = (v_sp.at(i)+v_sp_old.at(i))/2;

            if(e_lens.vert_sp_new.size() > 2) {
              int v1 = e_lens.slices.at(2*i)-1;
              int v2 = e_lens.slices.at(2*i+1)-1;
              e_lens.m->getPoint(e_lens.vert_ctr_new.at(v1), 0, pos_curr);
              e_lens.m->getPoint(e_lens.vert_ctr_new.at(v2), 0, v_temp);
              e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_curr/2
                                                                  +v_temp/2);
            }
            else {
              pos_curr = pos_curr - v_temp*f_calc.vdparam.dt;
              e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_curr);
            }
          }

          apf::Vector3 v_2c_pos(0,0,0);
          v_2c_pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
                    vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;

          e_lens.m->setPoint(e_lens.v_2c, 0, v_2c_pos);

        }
        else {
          assert(trial_type == 1);
          e_lens.m->getPoint(e_lens.vert_ctr_new.at(0), 0, v_temp);
          e_lens.m->getPoint(e_lens.vert_ctr_new.at(1), 0, pos_curr);
          e_lens.m->setPoint(e_lens.vert_sp_new.at(0), 0, pos_curr/2
                                                            +v_temp/2);
        }

      }

      for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        v_r_ctr_old.at(i) = v_r_ctr.at(i);
        v_t_ctr_old.at(i) = v_t_ctr.at(i);
        v_ctr_old.at(i) = v_ctr.at(i);
      }

      for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
        v_r_sp_old.at(i) = v_r_sp.at(i);
        v_t_sp_old.at(i) = v_t_sp.at(i);
        v_sp_old.at(i) = v_sp.at(i);
      }
      //apf::Vector3 ctr_pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) 
      //              + vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;
      //std::cout << "trial_type " << trial_type << std::endl;
      std::tie(r_min, r_max) = vd_find_sph_bound(); 

      //std::cout << "r_min " << r_min << " r_max " << r_max << std::endl;
      //if(dt_inv < dt_th/10) {
      //  f_calc.vdparam.adj_dt(dt_inv);
      //}
      counter++;
      if(mov_flag)
        save_vtk_mov(e_lens.m, mov_tag, "./outputtemp/series");
      mov_tag++;
    }
    std::cout << "r_min " << r_min << " r_max " << r_max << std::endl;
    if(counter == counter_sz) {
      std::cout << "Ext_loop: Takes too long to converge." << std::endl;
      e_i_old = e_i_new;
      e_f_old = e_f_new;
      e_f_new = calc_energy_lens();
      if(e_f_new > energy) {
        ins_flag = false;
        std::cout << "energy_flag: energy is increasing" << std::endl;
      }
    }
    else if(del_e >= 0) {
      std::cout << "Ext_loop: energy reached a minimum." << std::endl;
      e_i_old = e_i_new;
      e_f_old = e_f_new;
      e_f_new = calc_energy_lens();
      if(e_f_new > energy) {
        ins_flag = false;
        std::cout << "energy_flag: energy is increasing" << std::endl;
      }
    }

    else if(d_max < d_th) {
      std::cout << "Ext_loop: d_max " << d_max << " < " << d_th << " d_th " 
                << std::endl;
      iter = iter_sz;
      e_i_old = e_i_new;
      e_f_old = e_f_new;

      e_f_new = calc_energy_lens();
      if(e_f_new > energy) {
        ins_flag = false;
        std::cout << "energy_flag: energy is increasing" << std::endl;
      }

      //e_f_new = calc_energy_lens() + energy;
      e_f_new = calc_energy_lens();
      //r_th_min = 2*r_max;
      proj_around_ctr(rho);

      //apf::Vector3 ctr_pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) 
      //              + vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;
      std::tie(r_min, r_max) = vd_find_sph_bound(); 

      if(r_min < r_th_min) {
        if (r_min < std::numeric_limits<double>::epsilon() ) {
          proj_around_ctr(rho);
        }
        else
          proj_around_ctr(rho*2*r_th_min/r_min);
      }
      //e_i_new = calc_energy_lens() + energy;
      e_i_new = calc_energy_lens();
      iter++;
    }
    else if(r_min < r_th_min) {
      std::cout << "Ext_loop: closing" << std::endl;
      closing = true;
      iter = iter_sz;
      e_i_old = e_i_new;
      e_f_old = e_f_new;
    }
    else if (r_max > r_th_max) {
      std::cout << "Ext_loop: hit outer sphere" << std::endl;
      e_i_old = e_i_new;
      e_f_old = e_f_new;

      //e_f_new = calc_energy_lens() + energy;
      e_f_new = calc_energy_lens();

      proj_around_ctr(rho);

      apf::Vector3 ctr_pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
                    vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;

      std::tie(r_min, r_max) = vd_find_sph_bound(); 

      if(r_min < r_th_min*2/3) {
        if (r_min < std::numeric_limits<double>::epsilon() ) {
          proj_around_ctr(rho);
        }
        else
          proj_around_ctr(rho*2*r_th_min/r_min);
      }
      //e_i_new = calc_energy_lens() + energy;
      e_i_new = calc_energy_lens();
      iter++;
    }

    std::cout << "e_i_old " << e_i_old
              << " e_i_new " << e_i_new
              << " e_f_old " << e_f_old
              << " e_f_new " << e_f_new
              << std::endl;

    if(trial_type == 1) {
      e1.at(trial_curr) = e_i_new;
    }
    else {
      assert(trial_type == 2);
      e2.at(trial_curr).at(ng.curr) = e_i_new;
    }
  }
/*
  if(r_max < rho) {
    std::cout << "r_max " << r_max << " rho " << rho << std::endl;
    ins_flag = false;
    std::cout << "ins_flag: distance less than rho" << std::endl;
  }
*/
  if(closing) {
    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_old);
    }
    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_old);
    }
    if(trial_type == 2)
      e_lens.m->setPoint(e_lens.v_2c, 0, pos_old);
    ins_flag = false;
    std::cout << "ins_flag: closing" << std::endl;
  }
  if (e_lens.m == m_trial) {
    if (ins_flag) {
      proj_around_ctr(rho);

      for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr_new.at(i), drag_flag);
      }

      for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
        f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_sp_new.at(i), drag_flag);
      }
  /*
      for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        v_temp = (v_ctr.at(i)+v_ctr_old.at(i))/2;
        apf::setVector(vel_field, e_lens.vert_ctr_new.at(i), 0, v_temp);
      }
      for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
        v_temp = (v_sp.at(i)+v_sp_old.at(i))/2;
        apf::setVector(vel_field, e_lens.vert_sp_new.at(i), 0, v_temp);
      }
  */
      calc_energy_diss_rate();

      vd_entlist_v e_list_v(e_lens.m, &e_lens.vert_sp_new, e_lens.c_base_curr);

      int new_cell_id;
      if(trial_type == 1) {
        new_cell_id = e_lens.new_cell1_id.at(0);
      }
      else {
        assert(trial_type == 2);
        new_cell_id = e_lens.new_cell2_id;
      }

      find_pos_non_inv();
      double roc = f_calc.calc_roc(e_lens.m, &e_list_v.e.at(trial_type)
                                 .at(new_cell_id-1).at(trial_type));

      std::cout << "Rate of change of insertion " << roc << std::endl;

      if(!skip_en) {
        if(roc > std::numeric_limits<double>::epsilon()) {
          if(trial_type == 1) {
            w1.at(trial_curr) = calc_energy_diss_rate_sing();
          }
          else {
            assert(trial_type == 2);
            w2.at(trial_curr).at(ng.curr) = calc_energy_diss_rate_sing();
          }
          // Due to rotations or migration, some surrounding tets might invert
          // despite not closing. If that is the case, position the new vertices
          // such that no tet inverts.
        }
        else {
          if(trial_type == 1) {
            w1.at(trial_curr) = 0.1;
          }
          else {
            assert(trial_type == 2);
            w2.at(trial_curr).at(ng.curr) = 0.1;
          }
        }
      }
    }
    else {
      if(!skip_en) {
        if(trial_type == 1) {
          w1.at(trial_curr) = 0.1;
        }
        else {
          assert(trial_type == 2);
          w2.at(trial_curr).at(ng.curr) = 0.1;
        }
      }
    }
  }
  else {
    if (ins_flag) {
      //proj_around_ctr(rho);
      // Due to rotations or migration, some surrounding tets might invert
      // despite not closing. If that is the case, position the new vertices
      // such that no tet inverts.
      //if(vd_chk_neg_sgn(e_lens.m)) {
        find_pos_non_inv();
      //  assert(vd_chk_neg_sgn(e_lens.m));
      //}

    }
  }
  //else {
  //  rho = len_sp/4;
  //  proj_around_ctr(rho);
  //}
  f_calc.vdparam.adj_dt(dt_glob);
}

// Find the minimum and maximum distances from the 1stratum centers for vertices
// bounding 1strata and 2stratum center for vertices bounding the 2stratum.
std::pair<double, double> vd_edisc::vd_find_sph_bound() {
  apf::Vector3 p_ctr(0, 0, 0);
  double r_min = -1;
  double r_max = -1;
  double r_curr;
  std::vector<apf::Vector3> temp_ctr(e_lens.vert_ctr_new.size(), 
                                      apf::Vector3(0,0,0));

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, temp_ctr.at(i));
  }
  if(trial_type == 1) {
    r_curr = (temp_ctr.at(0) - temp_ctr.at(1)).getLength()/2;
    if(r_min < 0 or r_min > r_curr)
      r_min = r_curr;
    if(r_max < 0 or r_max < r_curr)
      r_max = r_curr;
  }
  else {
    assert(trial_type == 2);
    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      int v1 = e_lens.slices.at(2*i)-1;
      int v2 = e_lens.slices.at(2*i+1)-1;

      r_curr = (temp_ctr.at(v1) - temp_ctr.at(v2)).getLength()/2;
      if(r_min < 0 or r_min > r_curr)
        r_min = r_curr;
    }

    p_ctr = vd_get_center(e_lens.m, &e_lens.vert_ctr_new)/2 + 
            vd_get_center(e_lens.m, &e_lens.vert_sp_new)/2;

    for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      r_curr = (temp_ctr.at(i) - p_ctr).getLength();
      if(r_min < 0 or r_min > r_curr)
        r_min = r_curr;
      if(r_max < 0 or r_max < r_curr)
        r_max = r_curr;
    }

    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, temp_ctr.at(i));
      r_curr = (temp_ctr.at(i) - p_ctr).getLength();
      if(r_min < 0 or r_min > r_curr)
        r_min = r_curr;
      if(r_max < 0 or r_max < r_curr)
        r_max = r_curr;
    }
  }
  return std::make_pair(r_min, r_max);
}

double vd_edisc::get_r_max(std::vector<double> &r_ctr, std::vector<double> &r_sp) {
  double r_max = 0;
  apf::Vector3 temp(0,0,0);
  apf::Vector3 temp_r(0,0,0);

  assert(!isvecnan(pos_old));

  //std::cout << "pos_old " << pos_old << std::endl; 

  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, temp);

    assert(!isvecnan(temp));

    temp = temp - pos_old;
    r_ctr.at(i) = temp.getLength();
    if(r_ctr.at(i) > r_max)
      r_max = r_ctr.at(i);

    //std::cout << "v_ctr " << temp << std::endl; 
    //std::cout << "r_max " << r_max << " r_ctr " << r_ctr.at(i) << std::endl; 
  }
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, temp);

    assert(!isvecnan(temp));

    temp = temp - pos_old;
    r_sp.at(i) = temp.getLength();
    if(r_sp.at(i) > r_max)
      r_max = r_sp.at(i);
    //std::cout << "v_sp " << temp << std::endl; 
    //std::cout << "r_max " << r_max << " r_sp " << r_sp.at(i) << std::endl; 
  }

  if (r_max < std::numeric_limits<double>::epsilon() ) {
    r_max = 1;
    //std::cout << "All new vertices at pos_old, r_max = " << r_max 
    //          << std::endl; 
    itisnan();
  }
  return r_max;
}

void vd_edisc::proj_around_ctr(double rho) {
  std::vector<double> r_ctr(e_lens.vert_ctr_new.size());
  std::vector<double> r_sp(e_lens.vert_sp_new.size());
  double r_max = get_r_max(r_ctr, r_sp);

  apf::Vector3 temp(0,0,0);
  apf::Vector3 temp_r(0,0,0);
  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    if(calc_ext) {
      if(ext_new[e_lens.vert_ctr_new.at(i)]) {
        //std::cout << "Projecting exterior 0cell vertex " 
        //          << e_lens.vert_ctr_new.at(i)
        //          << std::endl; 
        e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, temp);

        //std::cout << "temp " << temp << std::endl; 
        //std::cout << "pos_old " << pos_old << std::endl; 
        //std::cout << "(- ).norm() " << norm_0(temp - pos_old) << std::endl; 
        //std::cout << "rho " << rho << " r_ctr " << r_ctr.at(i) 
        //          << " r_max " << r_max << std::endl; 
        temp_r = norm_0(temp - pos_old)*rho*r_ctr.at(i)/r_max;

        //std::cout << "Projected motion " << temp_r << std::endl;

        apf::Vector3 volm_pres(0,0,0);

        volm_pres = f_calc.get_vec_special(e_lens.m, 
                              e_lens.vert_ctr_new.at(i), temp_r);

        assert(!isvecnan(pos_old+volm_pres));
        e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_old + volm_pres);
        //std::cout << "Volm preserving " << volm_pres << std::endl;
        //std::cout << "Position " << pos_old + volm_pres << std::endl;
      }
      else {
        //std::cout << "Projecting interior 0cell vertex " 
        //          << e_lens.vert_ctr_new.at(i)
        //          << std::endl; 
          vd_proj_v_sphere(e_lens.m, e_lens.vert_ctr_new.at(i), pos_old,
                                                         rho*r_ctr.at(i)/r_max);
          e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, temp);
          assert(!isvecnan(temp));
        //std::cout << "Position " << temp << std::endl;
      }
    }
    else {
      if(!ext_new[e_lens.vert_ctr_new.at(i)]) {
        //std::cout << "Projecting interior 0cell vertex " 
        //          << e_lens.vert_ctr_new.at(i)
        //          << std::endl; 
          vd_proj_v_sphere(e_lens.m, e_lens.vert_ctr_new.at(i), pos_old,
                                                         rho*r_ctr.at(i)/r_max);
          e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, temp);
          assert(!isvecnan(temp));
        //std::cout << "Position " << temp << std::endl;
      }
    }
  }
  if(trial_type == 2) {
    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      if(calc_ext) {
        if(ext_new[e_lens.vert_sp_new.at(i)]) {
          //std::cout << "Projecting exterior 1cell vertex " 
          //          << e_lens.vert_sp_new.at(i)
          //          << std::endl; 
          e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, temp);
          temp_r = norm_0(temp - pos_old)*rho*r_sp.at(i)/r_max;

          //std::cout << "Projected motion " << temp_r << std::endl; 
          apf::Vector3 volm_pres(0,0,0);

          volm_pres = f_calc.get_vec_special(e_lens.m, 
                              e_lens.vert_sp_new.at(i), temp_r);

          assert(!isvecnan(pos_old+volm_pres));
          e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_old + volm_pres);
          //std::cout << "Volm preserving " << volm_pres << std::endl; 
          //std::cout << "Position " << pos_old + volm_pres << std::endl;
        }
        else {
          //std::cout << "Projecting interior 1cell vertex " 
          //          << e_lens.vert_sp_new.at(i)
          //          << std::endl; 
          vd_proj_v_sphere(e_lens.m, e_lens.vert_sp_new.at(i), pos_old,
                                                         rho*r_sp.at(i)/r_max);
          e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, temp);
          assert(!isvecnan(temp));
          //std::cout << "Position " << temp << std::endl;
        }
      }
      else if (!ext_new[e_lens.vert_sp_new.at(i)]) {
        //std::cout << "Projecting interior 1cell vertex " 
        //          << e_lens.vert_sp_new.at(i)
        //          << std::endl; 
        vd_proj_v_sphere(e_lens.m, e_lens.vert_sp_new.at(i), pos_old,
                                                       rho*r_sp.at(i)/r_max);
        e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, temp);
          assert(!isvecnan(temp));
        //std::cout << "Position " << temp << std::endl;
      }
    }
    if(f_calc.get_ext() and ext_0cell) {
      e_lens.pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
                    vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;

      assert(!isvecnan(e_lens.pos));
      if(f_calc.get_e_sh()->chk_shell(2, e_lens.new_cell2_id - 1))
        e_lens.pos = f_calc.get_e_sh()->find_int_pos(2, e_lens.new_cell2_id - 1,
                                                               e_lens.pos);
      e_lens.m->setPoint(e_lens.v_2c, 0, e_lens.pos);
      //std::cout << "Position " << e_lens.pos << std::endl;
    }
    else {
      e_lens.pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
                    vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;
      assert(!isvecnan(e_lens.pos));
      e_lens.m->setPoint(e_lens.v_2c, 0, e_lens.pos);
    }
  }
  else {
    assert(trial_type == 1);
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(0), 0, temp);
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(1), 0, temp_r);

    e_lens.pos = vd_get_center(e_lens.m, &e_lens.vert_ctr_new);
    e_lens.m->setPoint(e_lens.vert_sp_new.at(0), 0, temp/2+temp_r/2);  
  }
}

// Insert an edge inside the 3cell. Used in preconditioning step to treat
// disconnected 3cell couples.
// TODO this should be per disjoint 3cell. Also, for topologies such as a bead
// touching a surface, this will detach the bead! A better way would be 
// detecting the sets of 3cell tets that will be used to insert an edge.
// Also, the assumption that there is a single 3cell edge needs to be relaxed.
void vd_edisc::insert_3celledge(int cell3_curr) {
  if(cell3_curr == -1)
    return;

  trial_curr = cell3_curr;

  init_elens();
  e_lens.c_base_curr = c_base_act;

  collect_disc();
  color_discs();

  apf::Vector3 vt(0,0,0);
  apf::Vector3 vb(0,0,0);
  apf::Vector3 v_tet(0,0,0);

  double V_total = 0;
  double dist = -1;


  //vd_3c_det v3_cd(e_lens.m, e_lens.c_base_curr, e_lens.vert_ctr);
/*
Consider the convex hull composed of the center vertex, the end vertices of the edges on the disc and a reflection of the center vertex through the plane of the end vertices. Calculate the convex hull. If all vertices are on the hull, pick vd_find_dist_w(e_lens.m, &tet_3c, vert_ctr_act) as the direction.
If there are some vertices inside the hull, 
if 1, pick a point interior in the direction of n1/2+n2/2, n1, n2 normals of adjacent tris
if 2, pick the center of these two and check if interior ((adjacent tris side check. If not, preconditioning not possible)
If more than 3:
  Pick two, construct a line between them
   Pick a third, using the adjacent tris, check if a part of the line remains between the planes. If not, not feasible. 
  If possible, mark these intersections.
  If more vertices, pick a new one. Check the planes of the adjacent tris for intersections with the lines between the last vertex and the intersections from the previous line.
Recursively check for all vertices.
If a feasible intersection is available for the last vertex, pick the center of the vertex and the center of the last feasible line.

Have a tolerance when calculating intersections, so that the center is not too close to the triangle planes...

Start:
collect center vertex, edge end vertices, reflection of the center vertex
find the convex hull
assert the top and bottom vertices are in convex hull
find the interior vertices, if any
  If any, find the feasible interior position, if possible
*/
/*
  for(int i = 0; i < e_lens.es_elem.size(); i++) {
    apf::ModelEntity* mdl = e_lens.m->toModel(e_lens.es_elem.at(i));
    int m_type = e_lens.m->getModelType(mdl);
    int m_tag = e_lens.m->getModelTag(mdl);
    if(m_type == 3 and m_tag == trial_curr) {
      tet_3c.push_back(e_lens.es_elem.at(i));
    }
  }
*/
  std::vector<apf::MeshEntity*> tet_3c(0);
  tet_3c = vd_3c->get_tet_id(cell3_curr);

  // Better alternative to vd_find_dist_w, consider the intersections to the 
  // triangles across.
  std::vector<apf::MeshEntity*> tri(0);
  tri.reserve(tet_3c.size());
  apf::Downward d_t;
  apf::Downward d_v;
  int lookup_ts [4] = {2, 3, 1, 0};

  for(int i = 0; i < tet_3c.size(); i++) {
    m_act->getDownward(tet_3c.at(i), 0, d_v);
    m_act->getDownward(tet_3c.at(i), 2, d_t);

    int i1 = findIn(d_v, 4, vert_ctr_act);
    assert(i1 > -1);

    tri.push_back(d_t[lookup_ts[i1]]);
  }
  double len_main = 0;
  //double len_max = 0;
  //std::tie(len_main, len_max) = tri_dist_sphere(m_act, pos_old, &tri);

  if(!coll_bound_tri_dir(tet_3c, vb)){
    std::cout << "The direction " << vb 
              << " almost in plane with one of the bounding triangles." 
              << std::endl;
    pre_flag = false;
    return;
  }
  len_main = vd_dist_v_x_pl_dir(m_act, vert_ctr_act, tri, vb);

  //vb = vd_find_dist_w(e_lens.m, &tet_3c, vert_ctr_act); 
  std::cout << "vb before " << vb << std::endl;
  vb = norm_0(vb)*len_main;
  //assert(vb.getLength() < len_sh);

  // TODO actually, find the directions the vertices should move. Scale the 
  // motions using the largest motion. 
  std::cout << "vt of top slice " << vt << std::endl;
  std::cout << "vb of bottom slice " << vb << std::endl;
  double norm = sqrt((vt-vb)*(vt-vb));

  e_lens.vt.clear();
  e_lens.vt.resize(2);

  std::cout << "norm: " << norm 
            << ", len_sh: " << len_sh << std::endl;

  double z[3] = {0,0,0};
  e_lens.vt.at(0).fromArray(z);

  //e_lens.vt.at(1) = vb*len_sh*5/6/norm;
  //e_lens.vt.at(1) = vb/norm*dist_min*5/6;
  e_lens.vt.at(1) = vb;

  std::cout << "vt: " << e_lens.vt.at(0) 
            << ", vb: " << e_lens.vt.at(1) << std::endl;

  e_lens.m->getPoint(vert_ctr_act, 0, e_lens.pos);

  if(detect_inv()) {
    pre_flag = false;
    return;
  }

  e_lens.vert_ctr_new.resize(2);
  e_lens.edge_ctr_new_t.resize(1);
  e_lens.edge_ctr_new_b.resize(1);
  e_lens.vert_sp_new.resize(1);

  e_lens.vert_ctr_new.at(0) = vert_ctr_act;
  std::cout << "Old vert " << vert_ctr_act
            << " " << e_lens.m->getModelType(e_lens.m->toModel(vert_ctr_act))
            << "c" << e_lens.m->getModelTag(e_lens.m->toModel(vert_ctr_act))
            << std::endl;

  //apf::ModelEntity* vert_ctr_new_em = e_lens.m->findModelEntity(3, trial_curr);
  apf::ModelEntity* vert_ctr_new_em = vd_3c->get_mdl(3, trial_curr);

  // The first one is the actual vertex, second one is the new vertex, 
  // the internal entities after the expansion will be associated with the 
  // new vertex.
  // 3cell vertex vert_ctr_new.at(1).
  e_lens.vert_ctr_new.at(1) = e_lens.m->createVert(vert_ctr_new_em);
  std::cout << "New vert " << e_lens.vert_ctr_new.at(1)
            << " " << e_lens.m->getModelType(vert_ctr_new_em)
            << "c" << e_lens.m->getModelTag(vert_ctr_new_em)
            << std::endl;

  e_lens.vert_sp_new.at(0) = e_lens.m->createVert(vert_ctr_new_em);
  std::cout << "New vert_sp " << e_lens.vert_sp_new.at(0)
            << " " << e_lens.m->getModelType(vert_ctr_new_em)
            << "c" << e_lens.m->getModelTag(vert_ctr_new_em)
            << std::endl;


  apf::MeshEntity* v_edge[3];
  v_edge[0] = e_lens.vert_ctr_new.at(0);
  v_edge[1] = e_lens.vert_sp_new.at(0);
  v_edge[2] = e_lens.vert_ctr_new.at(1);

  //std::cout << "Old center vertex " << e_lens.vert_ctr << std::endl;
  //std::cout << "Edge vertices " << e_lens.vert_ctr_new.at(0) << " "
  //                              << e_lens.vert_sp_new.at(0) << " " 
  //                              << e_lens.vert_ctr_new.at(1) << " "
  //          << std::endl;


  e_lens.edge_ctr_new_t.at(0) = buildElement(e_lens.m, vert_ctr_new_em,
                                               apf::Mesh::EDGE, v_edge, 0);
  e_lens.edge_ctr_new_b.at(0) = buildElement(e_lens.m, vert_ctr_new_em, 
                                               apf::Mesh::EDGE, &v_edge[1], 0);

  //vd_print_vert(e_lens.m, e_lens.edge_ctr_new_t.at(0));
  //vd_print_vert(e_lens.m, e_lens.edge_ctr_new_t.at(0));

  //std::cout << "New edges " << e_lens.edge_ctr_new_t.at(0) << " "
  //                              << e_lens.edge_ctr_new_b.at(0) << " "
  //          << std::endl;

  int mdl_type = m_act->getModelType(vert_ctr_new_em);
  int mdl_tag = m_act->getModelTag(vert_ctr_new_em);

  std::cout << "Edge_t " << e_lens.edge_ctr_new_t.at(0) << " "
            << mdl_type << "c" << mdl_tag
            << " edge_b " << e_lens.edge_ctr_new_b.at(0) << " "
            << mdl_type << "c" << mdl_tag
            << std::endl; 

  e_lens.create_map[0][vert_ctr_act] = e_lens.vert_ctr_new.at(0);
  e_lens.create_map[1][vert_ctr_act] = e_lens.vert_ctr_new.at(1);

  e_lens.m->setPoint(e_lens.vert_ctr_new.at(0), 0, e_lens.pos);
  e_lens.m->setPoint(e_lens.vert_ctr_new.at(1), 0, e_lens.pos+e_lens.vt.at(1)*5/6);
  e_lens.m->setPoint(e_lens.vert_sp_new.at(0), 0, 
                                e_lens.pos+e_lens.vt.at(1)*4/6);

  std::cout << "p+vt: " << e_lens.pos+e_lens.vt.at(0)
            << ", p+vb: " << e_lens.pos+e_lens.vt.at(1)*4/6 << std::endl;

  // Expand the disc:
  expand_lens();
  fill_lens();
  fill_3c_lens_elem();

  // Recreate the slice entities:
  recreate_slices();

  e_lens.m->getPoint(e_lens.vert_ctr, 0, e_lens.pos);

  // Destroy the old entities:
  destroy_ent();
  e_lens.m->acceptChanges();

  f_calc.vd_att_fields(e_lens.m, e_lens.vert_sp_new.at(0));
  f_calc.vd_att_fields(e_lens.m, e_lens.vert_ctr_new.at(0));
  f_calc.vd_att_fields(e_lens.m, e_lens.vert_ctr_new.at(1));

  vd_chk_neg_sgn(e_lens.m);

  //int trial_temp = trial_type;
  //trial_type = 3;
  /*
  if (e_lens.m == m_trial)
    vtk_trial();
  else if (e_lens.m == m_precond) {
    vtk_precond();
  }
  */
  vtk_mesh();

  //trial_type = trial_temp;

}

void vd_edisc::collect_tris_wg(
                std::vector<std::pair<std::vector<int >, std::vector<int > > >*
                                                                  path_cells,
                std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells) {

  for(int i = 0; i < slice_tris.size(); i++) {
    slice_tris.at(i).clear();
  }
  slice_tris.resize(slice_cells->size());

  for(int i = 0; i < slice_tris.size(); i++) {
    slice_tris.at(i).reserve(e_lens.es_surf.size());
  }
  std::vector<apf::MeshEntity*> temp(0);

  for(int slice = 0; slice < slice_cells->size(); slice++) {
    for(int i = 0; i < slice_cells->at(slice).second.at(1).size(); i++) {
      int c2_local = slice_cells->at(slice).second.at(1).at(i);
      temp = vd_cd->get_dim_ent(2, c2_local);
      for(int j = 0; j < temp.size(); j++) {
        slice_tris.at(slice).push_back(temp.at(j));
      }
    }
    // Add the triangles associated with the path 2cells to all slices.
    int path_curr = slice_cells->at(slice).first.first;
    if(path_curr == 0) {
      assert(slice_cells->at(slice).first.second == 0);
      for(int i = 0; i < path_cells->at(path_curr).second.size(); i++) {
        int c2_local = path_cells->at(path_curr).second.at(i);
        temp = vd_cd->get_dim_ent(2, c2_local);
        for(int j = 0; j < temp.size(); j++) {
          slice_tris.at(slice).push_back(temp.at(j));
        }
      }
    }
    else {
      path_curr = path_curr - 1;
      for(int i = 0; i < path_cells->at(path_curr).second.size(); i++) {
        int c2_local = path_cells->at(path_curr).second.at(i);
        temp = vd_cd->get_dim_ent(2, c2_local);
        for(int j = 0; j < temp.size(); j++) {
          slice_tris.at(slice).push_back(temp.at(j));
        }
      }
      path_curr = slice_cells->at(slice).first.second - 1;
      for(int i = 0; i < path_cells->at(path_curr).second.size(); i++) {
        int c2_local = path_cells->at(path_curr).second.at(i);
        temp = vd_cd->get_dim_ent(2, c2_local);
        for(int j = 0; j < temp.size(); j++) {
          slice_tris.at(slice).push_back(temp.at(j));
        }
      }
    }
  }
}

void vd_edisc::collect_edges_wg(
                std::vector<std::pair<std::vector<int >, std::vector<int > > >*
                                                                  path_cells,
                std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells) {
  apf::Downward d_v;
  apf::Downward d_e;

  for(int i = 0; i < slice_edges.size(); i++) {
    slice_edges.at(i).clear();
  }
  slice_edges.resize(slice_cells->size());

  for(int i = 0; i < slice_edges.size(); i++) {
    slice_edges.at(i).reserve(e_lens.es_edge.size());
  }
  std::vector<apf::MeshEntity*> temp(0);

  for(int slice = 0; slice < slice_cells->size(); slice++) {
    for(int dim = 0; dim < 3; dim++) {
      for(int i = 0; i < slice_cells->at(slice).second.at(dim).size(); i++) {
        int c_local = slice_cells->at(slice).second.at(dim).at(i);
        apf::ModelEntity* mdl = vd_cd->get_mdl(dim + 1, c_local);

        temp = vd_cd->get_dim_ent(dim+1, c_local);
        if(dim == 0) { //1stratum
          for(int j = 0; j < temp.size(); j++) {
            slice_edges.at(slice).push_back(temp.at(j));
          }
        }
        else if (dim == 1) {
          for(int j = 0; j < temp.size(); j++) {
            e_lens.m->getDownward(temp.at(j), 0, d_v);
            e_lens.m->getDownward(temp.at(j), 1, d_e);

            int i_v = findIn(d_v, 3, e_lens.vert_ctr);
            for(int k = 0; k < 2; k++) {
              int k_curr = lookup_tri_ed[i_v][k];
              if(e_lens.m->toModel(d_e[k_curr]) == mdl)
                slice_edges.at(slice).push_back(d_e[k_curr]);
            }
          }
        }
        else if (dim == 2){
          for(int j = 0; j < temp.size(); j++) {
            e_lens.m->getDownward(temp.at(j), 0, d_v);
            e_lens.m->getDownward(temp.at(j), 1, d_e);

            int i_v = findIn(d_v, 4, e_lens.vert_ctr);
            for(int k = 0; k < 3; k++) {
              int k_curr = lookup_tet_edge[i_v][k];
              if(e_lens.m->toModel(d_e[k_curr]) == mdl)
                slice_edges.at(slice).push_back(d_e[k_curr]);
            }
          }
        }
      }
    }
    std::sort(slice_edges.at(slice).begin(), slice_edges.at(slice).end());
    std::vector<apf::MeshEntity*>::iterator it;
    it = std::unique (slice_edges.at(slice).begin(), slice_edges.at(slice).end());
    slice_edges.at(slice).resize(std::distance(slice_edges.at(slice).begin(),
                                                                           it));
/*
    // Actually this is not necessary as internal edges bounding coplanar 2-
    // stratum triangles have no effect for constant energy case.
    // Add the edges associated with the path cells to all slices.
    int path_curr = slice_cells->at(slice).first.first;
    // Circuit
    if(path_curr == 0) {
      assert(slice_cells->at(slice).first.second == 0);
      for(int i = 0; i < path_cells->at(path_curr).first.size(); i++) {
        int c3_local = path_cells->at(path_curr).first.at(i);
        apf::ModelEntity* mdl = vd_cd->get_mdl(3, c3_local);
        temp = vd_cd->get_dim_ent(3, c3_local);
        for(int j = 0; j < temp.size(); j++) {
          e_lens.m->getDownward(temp.at(j), 0, d_v);
          e_lens.m->getDownward(temp.at(j), 1, d_e);

          int i_v = findIn(d_v, 4, e_lens.vert_ctr);
          for(int k = 0; k < 3; k++) {
            int k_curr = lookup_tet_edge[i_v][k];
            if(e_lens.m->toModel(d_e[k_curr]) == mdl)
              slice_edges.at(slice).push_back(d_e[k_curr]);
          }
        }
      }
      for(int i = 0; i < path_cells->at(path_curr).second.size(); i++) {
        int c2_local = path_cells->at(path_curr).second.at(i);
        apf::ModelEntity* mdl = vd_cd->get_mdl(2, c2_local);
        temp = vd_cd->get_dim_ent(2, c2_local);
        for(int j = 0; j < temp.size(); j++) {
          e_lens.m->getDownward(temp.at(j), 0, d_v);
          e_lens.m->getDownward(temp.at(j), 1, d_e);

          int i_v = findIn(d_v, 3, e_lens.vert_ctr);
          for(int k = 0; k < 2; k++) {
            int k_curr = lookup_tri_ed[i_v][k];
            if(e_lens.m->toModel(d_e[k_curr]) == mdl)
              slice_edges.at(slice).push_back(d_e[k_curr]);
          }
        }
      }
    }

    // Set of paths
    else {
      path_curr = path_curr - 1;
      for(int i = 0; i < path_cells->at(path_curr).second.size(); i++) {
        int c2_local = path_cells->at(path_curr).second.at(i);
        apf::ModelEntity* mdl = vd_cd->get_mdl(2, c2_local);
        temp = vd_cd->get_dim_ent(2, c2_local);
        for(int j = 0; j < temp.size(); j++) {
          e_lens.m->getDownward(temp.at(j), 0, d_v);
          e_lens.m->getDownward(temp.at(j), 1, d_e);

          int i_v = findIn(d_v, 3, e_lens.vert_ctr);
          for(int k = 0; k < 2; k++) {
            int k_curr = lookup_tri_ed[i_v][k];
            if(e_lens.m->toModel(d_e[k_curr]) == mdl)
              slice_edges.at(slice).push_back(d_e[k_curr]);
          }
        }
      }
      path_curr = slice_cells->at(slice).first.second - 1;
      for(int i = 0; i < path_cells->at(path_curr).first.size(); i++) {
        int c3_local = path_cells->at(path_curr).first.at(i);
        apf::ModelEntity* mdl = vd_cd->get_mdl(3, c3_local);
        temp = vd_cd->get_dim_ent(3, c3_local);
        for(int j = 0; j < temp.size(); j++) {
          e_lens.m->getDownward(temp.at(j), 0, d_v);
          e_lens.m->getDownward(temp.at(j), 1, d_e);

          int i_v = findIn(d_v, 3, e_lens.vert_ctr);
          for(int k = 0; k < 3; k++) {
            int k_curr = lookup_tet_edge[i_v][k];
            if(e_lens.m->toModel(d_e[k_curr]) == mdl)
              slice_edges.at(slice).push_back(d_e[k_curr]);
          }
        }
      }
    }
    std::sort(slice_edges.at(slice).begin(), slice_edges.at(slice).end());
    it = std::unique (slice_edges.at(slice).begin(), slice_edges.at(slice).end());
    slice_edges.at(slice).resize(std::distance(slice_edges.at(slice).begin(),
                                                                           it));
*/
  }
}

void vd_edisc::update_circpath_energies(
       std::vector<std::pair<std::vector<int >, std::vector<int > > >*
                                                              path_cells) {
  apf::Field* gam2_field = e_lens.m->findField("gam2");
  apf::Field* d2_field = e_lens.m->findField("d2");

  for(int i = 0; i < e_lens.es_surf.size(); i++) {
    apf::setScalar(gam2_field, e_lens.es_surf.at(i), 0, f_calc.vdparam.vd_gam);
    apf::setScalar(d2_field, e_lens.es_surf.at(i), 0, 1/f_calc.vdparam.mob);
  }

  if(wg_tag == WG_TYPE::TRI) {
    std::vector<apf::MeshEntity*> temp(0);
    for(int path_curr = 0; path_curr < path_cells->size(); path_curr++) {
      for(int i = 0; i < path_cells->at(path_curr).second.size(); i++) {
        int c2_local = path_cells->at(path_curr).second.at(i);
        temp = vd_cd->get_dim_ent(2, c2_local);
        for(int j = 0; j < temp.size(); j++) {
          apf::setScalar(gam2_field, temp.at(j), 0, f_calc.vdparam.vd_gam/2);
          apf::setScalar(d2_field, e_lens.es_surf.at(i), 0, 1/f_calc.vdparam.mob/2);
        }
      }
    }
  }
}

bool vd_edisc::chk_inv_wg() {
  collect_disc();
  color_discs();

  std::vector<apf::Vector3> pos(4, apf::Vector3(0,0,0));
  apf::Downward down;
  apf::Downward d_tri;

  for (int i = 0; i < e_lens.discs.size(); i++) {
    int v1 = e_lens.slices.at(2*i)-1;
    int v2 = e_lens.slices.at(2*i+1)-1;
    apf::MeshEntity* vert_1 = e_lens.vert_ctr_new.at(v1);
    apf::MeshEntity* vert_2 = e_lens.vert_ctr_new.at(v2);

    // The inside tetrahedra. The vertex ordering is taken from upward adjacency
    // tets. Replace the old vert_ctr with vert_sp_new and the slice vertex with
    // the corresponding vert_ctr_new.
    apf::Up up;
    for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {
      apf::MeshEntity* tri_curr = e_lens.discs.at(i).tri.at(j);
      e_lens.m->getUp(tri_curr, up);
      assert(up.n == 2);
      assert(e_lens.slice_map[up.e[0]] != e_lens.slice_map[up.e[1]]);

      // Top tet:
      e_lens.m->getDownward(up.e[0], 0, down);
      e_lens.m->getDownward(up.e[0], 2, d_tri);
      int j_tri = findIn(d_tri, 4, tri_curr);

      // Replace the vertices:
      int v_slice_id = lookup_tet_x_surf[j_tri];
      int j_v = findIn(down, 4, vert_ctr_act);
      assert(j_v > -1);

      int slice_curr = e_lens.slice_map[up.e[0]];
      assert(slice_curr > 0);

      for(int k = 0; k < 4; k++) {
        e_lens.m->getPoint(down[k], 0, pos.at(k));
      }

      // TODO For some reason this is not the correct ordering?
      int i_shift = 0;
      if(slice_curr == e_lens.slices.at(2*i)) {
        pos.at(v_slice_id) = pos_old + norm_0(e_lens.vt.at(slice_curr - 1))*len_sh;
      }
      else {
        assert(slice_curr == e_lens.slices.at(2*i+1));
        pos.at(v_slice_id) = pos_old + norm_0(e_lens.vt.at(slice_curr - 1))*len_sh;
        i_shift = 1;
      }
      if(vd_volume_tet(&pos) < std::numeric_limits<double>::min())
        return false;
    }
  }
  return true;
}
// Calculate the maximum dissipation rate insertion without modifying the mesh. 
void vd_edisc::calc_max_diss_trial_wg() {
  //init_elens();
  apf::Field* vel_field = e_lens.m->findField("velocity_field");

  trial_type = 1;
  std::vector<apf::Vector3> force_ctr(2, apf::Vector3(0,0,0));

  std::vector<std::pair<std::vector<int >, std::vector<int > > > path_cells
      (0, std::make_pair(std::vector<int>(0), std::vector<int>(0)));
  std::vector< std::pair< std::pair<int,int>, 
                         std::vector<std::vector<int > > > > slice_cells
      (0, std::make_pair(std::make_pair(0,0), std::vector<std::vector<int > > (0, 
                                                      std::vector<int > (0))));

  w1.resize(circ_sz);
  e1.resize(circ_sz);


  int step = 10;
  w1_exp.resize(circ_sz);
  for(int i = 0; i < circ_sz; i++)
    w1_exp.at(i).resize(step);

  apf::Vector3 vel_ctr(0,0,0);
  apf::Vector3 f_ctr(0,0,0);
  f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr, false);
  apf::getVector(vel_field, e_lens.vert_ctr, 0, vel_ctr);
  f_ctr = f_calc.vd_calc_force(e_lens.m, e_lens.vert_ctr, false);
  if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
    ext_shell* e_sh = f_calc.get_e_sh();
    if(e_sh->chk_shell(0, cell_id - 1) ) {
      vel_ctr = e_sh->find_para_dir(0, cell_id - 1, vel_ctr);
      f_ctr = e_sh->find_para_dir(0, cell_id - 1, f_ctr);
    }
  }
  double diss_base = (vel_ctr*f_ctr)*(-1);
  std::cout << "p_diss_ctr = " << diss_base << std::endl;

  for(trial_curr = 0; trial_curr < vd_cd->get_circ_sz(); trial_curr++) {
    if(vd_cd->get_circ_type(trial_curr) < 2) {
      ins_flag = true;

      double p_diss = 0;

      init_elens();
      vd_cd->find_slice_circ(trial_curr, &path_cells, &slice_cells);
      collect_cell_wg(&slice_cells);

      crt_new_vert_list();
      crt_new_cell();

      e_lens.slices.at(0) = 1;
      e_lens.slices.at(1) = 2;

      e_lens.vt.resize(2);
      e_lens.split_vt.resize(1);

      // Update the energies for circuits and path triangles to be half.
      update_circpath_energies(&path_cells);

      reload_trial_wg(&path_cells, &slice_cells);

      upd_cell_wg();

      e_lens.c_base_curr->coll_cell_gmi(1, e_lens.new_cell1_id.at(0), cell_id);
      if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
        //f_calc.set_e_sh(&e_sh_save);
        for (int i = 0; i < e_lens.new_cell0_id.size(); i++)
          f_calc.get_e_sh()->set_shell(0, e_lens.new_cell0_id.at(i) - 1, 
                                         sh_old, false);
        f_calc.get_e_sh()->set_shell(1, e_lens.new_cell1_id.at(0) - 1, 
                                         sh_old, false);
        f_calc.get_e_sh()->set_shell(0, cell_id - 1, sh_old, sh_ext);
      }
    }
  }

  trial_type = 2;
  // TODO fix the circular calls required to restore local ids of the ng object
  // after replacing them with gmi ids. There should be a separation between the
  // two,
  // ng with local ids is used in extracting the entities.
  // ng with gmi ids is used to update the stratification.
  refresh_e();

  for(trial_curr = 0; trial_curr < pt_sz; trial_curr++) {
    // Collect the relevant cells for the current 3cell couple to be joined.
    ng.clear();
    get_path_ngon(trial_curr, &ng);
    ng.print();
    int ngon_sz = ng.ngons.size();

    w2.at(trial_curr).resize(ngon_sz);
    e2.at(trial_curr).resize(ngon_sz);
    w2_exp.at(trial_curr).resize(ngon_sz);

    for(int ngons_curr = 0; ngons_curr < ngon_sz; ngons_curr++) {
      ng.clear();
      get_path_ngon(trial_curr, &ng);
      ng.curr = ngons_curr;

      vd_cd->set_ng_couple(trial_curr);

      ins_flag = true;
      double p_diss = 0;
      init_elens();

      if(ng.cp.first != -1)
        e_lens.e_3c_t_flag = true;
      else
        e_lens.e_3c_t_flag = false;

      if(ng.cp.second != -1)
        e_lens.e_3c_b_flag = true;
      else
        e_lens.e_3c_b_flag = false;

      vd_cd->find_slice_ng(ng.curr, &path_cells, &slice_cells);
      collect_cell_wg(&slice_cells);

      force_ctr.resize(slice_cells.size());

      crt_new_vert_list();
      crt_new_cell();

      for(int i = 0; i < slice_cells.size(); i++) {
        int path_curr = slice_cells.at(i).first.first - 1;
        if(e_lens.slices.at(2*path_curr) == 0)
          e_lens.slices.at(2*path_curr) = i + 1;

        else {
          assert(e_lens.slices.at(2*path_curr + 1) == 0);
          e_lens.slices.at(2*path_curr + 1) = i + 1;
        }
        path_curr = slice_cells.at(i).first.second - 1;
        if(e_lens.slices.at(2*path_curr) == 0)
          e_lens.slices.at(2*path_curr) = i + 1;

        else {
          assert(e_lens.slices.at(2*path_curr + 1) == 0);
          e_lens.slices.at(2*path_curr + 1) = i + 1;
        }
      }

      e_lens.vt.resize(ng.ngons.at(ng.curr).size());
      e_lens.split_vt.resize(ng.ngons.at(ng.curr).size());

      // Update the energies for circuits and path triangles to be half.
      update_circpath_energies(&path_cells);
      reload_trial_wg(&path_cells, &slice_cells);

      upd_cell_wg();

      e_lens.c_base_curr->coll_cell_gmi(2, e_lens.new_cell2_id, cell_id);
      if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
        //f_calc.set_e_sh(&e_sh_save);
        for (int i = 0; i < e_lens.new_cell0_id.size(); i++)
          f_calc.get_e_sh()->set_shell(0, e_lens.new_cell0_id.at(i) - 1, 
                                         sh_old, false);
        for (int i = 0; i < e_lens.new_cell1_id.size(); i++)
          f_calc.get_e_sh()->set_shell(1, e_lens.new_cell1_id.at(i) - 1, 
                                         sh_old, false);
        f_calc.get_e_sh()->set_shell(2, e_lens.new_cell2_id - 1, sh_old, false);

        f_calc.get_e_sh()->set_shell(0, cell_id - 1, sh_old, sh_ext);
      }
    }
  }

  int shift = 0;
  for (int i = 0; i < w1.size(); i++) {
    //std::cout << "Energy of ngon " << i << ": " << w2.at(i) << std::endl;
    std::cout << "E(i," << i+1+shift << ")= " << w1.at(i) 
              << ";" << std::endl;
  }
  shift = w1.size();
  for (trial_curr = 0; trial_curr < pt_sz; trial_curr++) {
    for (int i = 0; i < w2.at(trial_curr).size(); i++) {
      //std::cout << "Energy of ngon " << i << ": " << w2.at(i) << std::endl;
      std::cout << "E(i," << i+1+shift << ")= " << w2.at(trial_curr).at(i) 
                << ";" << std::endl;
    }
    shift = shift + w2.at(trial_curr).size();
  }

}


// Calculate the maximum dissipation rate insertion without modifying the mesh. 
void vd_edisc::calc_max_diss_wg() {
  //init_elens();
  apf::Field* vel_field = e_lens.m->findField("velocity_field");

  trial_type = 1;
  std::vector<apf::Vector3> force_ctr(2, apf::Vector3(0,0,0));

  std::vector<std::pair<std::vector<int >, std::vector<int > > > path_cells
      (0, std::make_pair(std::vector<int>(0), std::vector<int>(0)));
  std::vector< std::pair< std::pair<int,int>, 
                         std::vector<std::vector<int > > > > slice_cells
      (0, std::make_pair(std::make_pair(0,0), std::vector<std::vector<int > > (0, 
                                                      std::vector<int > (0))));

  w1.resize(circ_sz);
  e1.resize(circ_sz);
  w1_exp.resize(circ_sz);
  int step = 10;
  for(int i = 0; i < circ_sz; i++)
    w1_exp.at(i).resize(step);

  apf::Vector3 vel_ctr(0,0,0);
  apf::Vector3 f_ctr(0,0,0);
  f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr, false);
  apf::getVector(vel_field, e_lens.vert_ctr, 0, vel_ctr);
  f_ctr = f_calc.vd_calc_force(e_lens.m, e_lens.vert_ctr, false);
  if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
    ext_shell* e_sh = f_calc.get_e_sh();
    if(e_sh->chk_shell(0, cell_id - 1) ) {
      vel_ctr = e_sh->find_para_dir(0, cell_id - 1, vel_ctr);
      f_ctr = e_sh->find_para_dir(0, cell_id - 1, f_ctr);
    }
  }
  double diss_base = (vel_ctr*f_ctr)*(-1);
  std::cout << "p_diss_ctr = " << diss_base << std::endl;

  for(trial_curr = 0; trial_curr < vd_cd->get_circ_sz(); trial_curr++) {
    if(vd_cd->get_circ_type(trial_curr) < 2) {
      ins_flag = true;

      double p_diss = 0;

      init_elens();
      vd_cd->find_slice_circ(trial_curr, &path_cells, &slice_cells);
      collect_cell_wg(&slice_cells);

      if(wg_tag == WG_TYPE::TRI) {
        collect_tris_wg(&path_cells, &slice_cells);
      }
      if(wg_tag == WG_TYPE::EDGE) {
        collect_edges_wg(&path_cells, &slice_cells);
      }
      else {
        collect_tris_wg(&path_cells, &slice_cells);
      }

      crt_new_vert_list();
      crt_new_cell();

      e_lens.slices.at(0) = 1;
      e_lens.slices.at(1) = 2;

      e_lens.vt.resize(2);
      e_lens.split_vt.resize(1);

      // Update the energies for circuits and path triangles to be half.
      update_circpath_energies(&path_cells);

      // Don't fix the exterior velocities, yet.
      //calc_ctr_vel_wg();
      if(wg_tag == WG_TYPE::TRI) {
        for(int i = 0; i < slice_tris.size(); i++) {
          f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr, &slice_tris.at(i), 
                                                                        false);
          apf::getVector(vel_field, e_lens.vert_ctr, 0, e_lens.vt.at(i));
        }
        for(int i = 0; i < slice_tris.size(); i++) {
          force_ctr.at(i) = f_calc.vd_calc_force_tri(e_lens.m, e_lens.vert_ctr, 
                                                      &slice_tris.at(i), false);
        }
      }
      else if(wg_tag == WG_TYPE::EDGE) {
        for(int i = 0; i < slice_edges.size(); i++) {
          f_calc.vd_upd_vel_field_edge(e_lens.m, e_lens.vert_ctr, 
                                                             &slice_edges.at(i), 
                                                                        false);
          apf::getVector(vel_field, e_lens.vert_ctr, 0, e_lens.vt.at(i));
        }
        for(int i = 0; i < slice_edges.size(); i++) {
          force_ctr.at(i) = f_calc.vd_calc_force_edge(e_lens.m, e_lens.vert_ctr, 
                                                      &slice_edges.at(i), false);
        }


      }
      else {
        for(int i = 0; i < slice_tris.size(); i++) {
          f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr, &slice_tris.at(i), 
                                                                        false);
          apf::getVector(vel_field, e_lens.vert_ctr, 0, e_lens.vt.at(i));
        }
        for(int i = 0; i < slice_tris.size(); i++) {
          force_ctr.at(i) = f_calc.vd_calc_force_tri(e_lens.m, e_lens.vert_ctr, 
                                                      &slice_tris.at(i), false);
        }
      }

      upd_cell_wg();

      if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
        ext_shell* e_sh = f_calc.get_e_sh();
        for(int i = 0; i < slice_tris.size(); i++) {
          if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(i) - 1) ) {
            e_lens.vt.at(i) = 
              e_sh->find_para_dir(0, e_lens.new_cell0_id.at(i) - 1, 
                                                        e_lens.vt.at(i));
            force_ctr.at(i) = 
              e_sh->find_para_dir(0, e_lens.new_cell0_id.at(i) - 1, 
                                                        force_ctr.at(i));
          }
        }
      }
      // Also, check for inversions.
      //ins_flag = (ins_flag and chk_inv_wg());

      //calc_ctr_force_wg(force_ctr);
      if(ins_flag) {
        for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
          p_diss = p_diss - force_ctr.at(i)*e_lens.vt.at(i);
        }
        w1.at(trial_curr) = p_diss;
        std::cout << "p_diss(" << trial_curr << ") = " 
                  << p_diss << ";" << std::endl;
        if(w1.at(trial_curr) < en_min_1c) {
          en_min_1c = w1.at(trial_curr);
          circ_min = trial_curr;
          std::cout << "circ_min " << circ_min << std::endl;
        }
      }
      else
        w1.at(trial_curr) = 0.1;

      e_lens.c_base_curr->coll_cell_gmi(1, e_lens.new_cell1_id.at(0), cell_id);
      if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
        //f_calc.set_e_sh(&e_sh_save);
        for (int i = 0; i < e_lens.new_cell0_id.size(); i++)
          f_calc.get_e_sh()->set_shell(0, e_lens.new_cell0_id.at(i) - 1, 
                                         sh_old, false);
        f_calc.get_e_sh()->set_shell(1, e_lens.new_cell1_id.at(0) - 1, 
                                         sh_old, false);
        f_calc.get_e_sh()->set_shell(0, cell_id - 1, sh_old, sh_ext);
      }
    }
  }

  trial_type = 2;
  // TODO fix the circular calls required to restore local ids of the ng object
  // after replacing them with gmi ids. There should be a separation between the
  // two,
  // ng with local ids is used in extracting the entities.
  // ng with gmi ids is used to update the stratification.
  refresh_e();

  for(trial_curr = 0; trial_curr < pt_sz; trial_curr++) {
    // Collect the relevant cells for the current 3cell couple to be joined.
    ng.clear();
    get_path_ngon(trial_curr, &ng);
    ng.print();
    int ngon_sz = ng.ngons.size();

    w2.at(trial_curr).resize(ngon_sz);
    e2.at(trial_curr).resize(ngon_sz);
    w2_exp.at(trial_curr).resize(ngon_sz);

    for(int ngons_curr = 0; ngons_curr < ngon_sz; ngons_curr++) {
      ng.clear();
      get_path_ngon(trial_curr, &ng);
      ng.curr = ngons_curr;

      vd_cd->set_ng_couple(trial_curr);

      ins_flag = true;
      double p_diss = 0;
      init_elens();

      if(ng.cp.first != -1)
        e_lens.e_3c_t_flag = true;
      else
        e_lens.e_3c_t_flag = false;

      if(ng.cp.second != -1)
        e_lens.e_3c_b_flag = true;
      else
        e_lens.e_3c_b_flag = false;

      vd_cd->find_slice_ng(ng.curr, &path_cells, &slice_cells);
      collect_cell_wg(&slice_cells);
      if(wg_tag == WG_TYPE::TRI) {
        collect_tris_wg(&path_cells, &slice_cells);
      }
      if(wg_tag == WG_TYPE::EDGE) {
        collect_edges_wg(&path_cells, &slice_cells);
      }
      else {
        collect_tris_wg(&path_cells, &slice_cells);
      }

      force_ctr.resize(slice_cells.size());

      crt_new_vert_list();
      crt_new_cell();

      for(int i = 0; i < slice_cells.size(); i++) {
        int path_curr = slice_cells.at(i).first.first - 1;
        if(e_lens.slices.at(2*path_curr) == 0)
          e_lens.slices.at(2*path_curr) = i + 1;

        else {
          assert(e_lens.slices.at(2*path_curr + 1) == 0);
          e_lens.slices.at(2*path_curr + 1) = i + 1;
        }
        path_curr = slice_cells.at(i).first.second - 1;
        if(e_lens.slices.at(2*path_curr) == 0)
          e_lens.slices.at(2*path_curr) = i + 1;

        else {
          assert(e_lens.slices.at(2*path_curr + 1) == 0);
          e_lens.slices.at(2*path_curr + 1) = i + 1;
        }
      }

      e_lens.vt.resize(ng.ngons.at(ng.curr).size());
      e_lens.split_vt.resize(ng.ngons.at(ng.curr).size());

      // Update the energies for circuits and path triangles to be half.
      update_circpath_energies(&path_cells);
      // Don't fix the exterior velocities, yet.
      //calc_ctr_vel_wg();

      if(wg_tag == WG_TYPE::TRI) {
        for(int i = 0; i < slice_tris.size(); i++) {
          f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr, &slice_tris.at(i), 
                                                                        false);
          apf::getVector(vel_field, e_lens.vert_ctr, 0, e_lens.vt.at(i));
        }
        for(int i = 0; i < slice_tris.size(); i++) {
          force_ctr.at(i) = f_calc.vd_calc_force_tri(e_lens.m, e_lens.vert_ctr, 
                                                      &slice_tris.at(i), false);
        }
      }
      else if(wg_tag == WG_TYPE::EDGE) {
        for(int i = 0; i < slice_edges.size(); i++) {
          f_calc.vd_upd_vel_field_edge(e_lens.m, e_lens.vert_ctr, 
                                                             &slice_edges.at(i), 
                                                                        false);
          apf::getVector(vel_field, e_lens.vert_ctr, 0, e_lens.vt.at(i));
        }
        for(int i = 0; i < slice_edges.size(); i++) {
          force_ctr.at(i) = f_calc.vd_calc_force_edge(e_lens.m, e_lens.vert_ctr, 
                                                      &slice_edges.at(i), false);
        }
      }
      else {
        for(int i = 0; i < slice_tris.size(); i++) {
          f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr, &slice_tris.at(i), 
                                                                        false);
          apf::getVector(vel_field, e_lens.vert_ctr, 0, e_lens.vt.at(i));
        }
        for(int i = 0; i < slice_tris.size(); i++) {
          force_ctr.at(i) = f_calc.vd_calc_force_tri(e_lens.m, e_lens.vert_ctr, 
                                                      &slice_tris.at(i), false);
        }
      }

      upd_cell_wg();

      if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
        ext_shell* e_sh = f_calc.get_e_sh();
        for(int i = 0; i < slice_tris.size(); i++) {
          if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(i) - 1) ) {
            e_lens.vt.at(i) = 
              e_sh->find_para_dir(0, e_lens.new_cell0_id.at(i) - 1, 
                                                        e_lens.vt.at(i));
            force_ctr.at(i) = 
              e_sh->find_para_dir(0, e_lens.new_cell0_id.at(i) - 1, 
                                                        force_ctr.at(i));
          }
        }
      }
      // Also, check for inversions.
      //ins_flag = (ins_flag and chk_inv_wg());

      if(ins_flag) {
        //calc_ctr_force_wg(force_ctr);
        for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
          p_diss = p_diss - force_ctr.at(i)*e_lens.vt.at(i);
        }
        w2.at(trial_curr).at(ng.curr) = p_diss;
        std::cout << "p_diss(" << circ_sz + ng.curr << ") = " 
                  << p_diss << ";" << std::endl;
        if(w2.at(trial_curr).at(ng.curr) < en_min_2c) {
          en_min_2c = w2.at(trial_curr).at(ng.curr);
          ngon_min = ng.curr;
          path_min = trial_curr;
          std::cout << "path_min " << path_min << " ngon_min " 
                    << ngon_min << std::endl;
        }
      }
      else
        w2.at(trial_curr).at(ng.curr) = 0.1;

      e_lens.c_base_curr->coll_cell_gmi(2, e_lens.new_cell2_id, cell_id);
      if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
        //f_calc.set_e_sh(&e_sh_save);
        for (int i = 0; i < e_lens.new_cell0_id.size(); i++)
          f_calc.get_e_sh()->set_shell(0, e_lens.new_cell0_id.at(i) - 1, 
                                         sh_old, false);
        for (int i = 0; i < e_lens.new_cell1_id.size(); i++)
          f_calc.get_e_sh()->set_shell(1, e_lens.new_cell1_id.at(i) - 1, 
                                         sh_old, false);
        f_calc.get_e_sh()->set_shell(2, e_lens.new_cell2_id - 1, sh_old, false);

        f_calc.get_e_sh()->set_shell(0, cell_id - 1, sh_old, sh_ext);
      }
    }
  }

  int shift = 0;
  for (int i = 0; i < w1.size(); i++) {
    //std::cout << "Energy of ngon " << i << ": " << w2.at(i) << std::endl;
    std::cout << "E(i," << i+1+shift << ")= " << w1.at(i) 
              << ";" << std::endl;
  }
  shift = w1.size();
  for (trial_curr = 0; trial_curr < pt_sz; trial_curr++) {
    for (int i = 0; i < w2.at(trial_curr).size(); i++) {
      //std::cout << "Energy of ngon " << i << ": " << w2.at(i) << std::endl;
      std::cout << "E(i," << i+1+shift << ")= " << w2.at(trial_curr).at(i) 
                << ";" << std::endl;
    }
    shift = shift + w2.at(trial_curr).size();
  }
}

// Calculate the velocities without generating new vertices.
void vd_edisc::calc_ctr_vel_wg() {

  apf::Field* vel_field = e_lens.m->findField("velocity_field");
  std::vector<std::vector<apf::MeshEntity*> > tri_slice
                 (e_lens.slices.size(), std::vector<apf::MeshEntity*> (0) );

  e_lens.vt.clear();
  e_lens.vt.resize(e_lens.slices.size());

  for (int i = 0; i < e_lens.slices.size(); i++) {
    collect_slice_tri(&tri_slice.at(i), i+1);
  }
  e_lens.m->getPoint(vert_ctr_act, 0, e_lens.pos);

  for(int i = 0; i < e_lens.slices.size(); i++) {
    f_calc.vd_upd_vel_field(e_lens.m, vert_ctr_act, &tri_slice.at(i), false);
    //vd_upd_vel_field_mason(e_lens.m, vert_ctr_act, 
    //                                    &tri_slice.at(i), &f_calc);
    apf::getVector(vel_field, vert_ctr_act, 0, e_lens.vt.at(i));
    //std::cout << "calc_ctr[" << i << "] " << e_lens.vt.at(i) << std::endl;
  }

  for (int i = 0; i < tri_slice.size(); i++) {
    tri_slice.at(i).clear();
  }
  tri_slice.clear();

}

// Calculate the velocities without generating new vertices.
void vd_edisc::calc_sp_vel_wg() {

  apf::Field* vel_field = e_lens.m->findField("velocity_field");
  std::vector<std::vector<apf::MeshEntity*> > tri_path
                 (0, std::vector<apf::MeshEntity*> (0) );

  e_lens.split_vt.clear();
  e_lens.split_vt.resize(ng.ngons.at(ng.curr).size());
  tri_path.resize(ng.ngons.at(ng.curr).size());

  for (int i = 0; i < ng.ngons.at(ng.curr).size(); i++) {
    collect_path_tri(&tri_path.at(i), i);
  }

  for(int i = 0; i < tri_path.size(); i++) {
    f_calc.vd_upd_vel_field(e_lens.m, vert_ctr_act, &tri_path.at(i), false);
    //vd_upd_vel_field_mason(e_lens.m, vert_ctr_act, &tri_path.at(i), &f_calc);
    apf::getVector(vel_field, vert_ctr_act, 0, e_lens.split_vt.at(i));
    //std::cout << "calc_sp[" << i << "] " << e_lens.split_vt.at(i) 
    //          << std::endl;
  }

  for (int i = 0; i < tri_path.size(); i++) {
    tri_path.at(i).clear();
  }
  tri_path.clear();
}

// Going over all slices of all possible insertions, find the dt such that
// the maximum of vertex motion is less than the half of the edge lengths.
double vd_edisc::calc_dt() {

  //double max_vel = 0;
  //double vel_curr;
  double min_t = -1;

  trial_type = 1;

  std::vector<apf::Vector3> a_pos(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> area(0, apf::Vector3(0,0,0));

  vd_get_apos(e_lens.m, e_lens.vert_ctr, &e_lens.es_elem, &area, &a_pos);

  for(trial_curr = 0; trial_curr < vd_cd->get_circ_sz(); trial_curr++) {
    if(vd_cd->get_circ_type(trial_curr) < 2) {
      init_elens();
      e_lens.c_base_curr = c_base_act;

      collect_disc();
      color_discs();
  /*
      calc_ctr_vel_wg();

      for(int i = 0; i < e_lens.slices.size(); i++) {
        vel_curr = e_lens.vt.at(i)*e_lens.vt.at(i);
        std::cout << "\tvel_curr^2: " << vel_curr << std::endl; 
        if(vel_curr > max_vel)
          max_vel = vel_curr;
      }
  */
      calc_ctr_vel_wg();

      double t_ctr = vd_find_min_t_v(e_lens.m, e_lens.vert_ctr, 
                                          &area, &a_pos, &e_lens.vt);

      std::cout << "t_ctr " << t_ctr
                << " t_min " << min_t << std::endl;

      if(min_t < 0)
        min_t = t_ctr;
      else {
        min_t = std::min(t_ctr, min_t);
      }
    }
    if(std::isnan(min_t))
      min_t = f_calc.vdparam.dt;

  }

  trial_type = 2;

  for(trial_curr = 0; trial_curr < pt_sz; trial_curr++) {

    // Collect the relevant cells for the current 3cell couple to be joined.
    get_path_ngon(trial_curr, &ng);

    for(int i = 0; i < ng.ngons.size(); i++) {
      init_elens();
      // init_elens called internally.
      collect_disc_path(i);
/*
      calc_ctr_vel_wg();

      for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        vel_curr = e_lens.vt.at(i)*e_lens.vt.at(i);
        std::cout << "\tvel_curr^2: " << vel_curr << std::endl; 
        if(vel_curr > max_vel)
          max_vel = vel_curr;
      }
*/
      calc_ctr_vel_wg();
      calc_sp_vel_wg();

      double t_ctr = vd_find_min_t_v(e_lens.m, e_lens.vert_ctr, 
                                          &area, &a_pos, &e_lens.vt);
      if(min_t < 0) {
        min_t = t_ctr;
      }

      double t_sp = vd_find_min_t_v(e_lens.m, e_lens.vert_ctr, 
                                          &area, &a_pos, &e_lens.split_vt);

      std::cout << "t_ctr " << t_ctr
                << " t_sp " << t_sp
                << " t_min " << min_t << std::endl;
      if(t_ctr < t_sp)
        min_t = std::min(t_ctr, min_t);
      else
        min_t = std::min(t_sp, min_t);
    }
  }
  std::cout << "dt_min: " << min_t << " " << f_calc.vdparam.dt << std::endl;
  if (min_t < 0 or std::isnan(min_t)) 
    return f_calc.vdparam.dt;
  else
    return std::min(min_t/2, f_calc.vdparam.dt);

}


// 1cell insertion with the current circuit.
void vd_edisc::insert_1cell(int circ) {
  //assert(trial_ex and trial_load);

  // TODO simplify. Two checks: Is the current 0c exterior? Is the ext_shell
  // algorithm being used. If both, shell needs to be updated before 
  // find_ext_dir. upd_cell changes the topology information of the 0cell.
  bool shell_ext = false;
  bool cell_ext = false;
  if(calc_ext and e_lens.c_base_curr->get_cell_ext_gmi(0, cell_id)) {
    cell_ext = true;
  }

  if(cell_ext) {
    if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
      shell_ext = true;
    }
  }

  trial_type = 1;
  trial_curr = circ;

  
  // TODO Assert statement for circuit number.
  assert(trial_curr < circ_sz);

  // TIME
  double t1;
  double dt_init = 0;
  double dt_collect = 0;
  double dt_color_disc = 0;
  double dt_reload_inv = 0;

  double dt_crt_ent = 0;

  double dt_calc_vel_init = 0;
  double dt_upd_pos_init = 0;

  double dt_relax = 0;
  double dt_vtk_rest = 0;
  double dt_tot = 0;

  double t0 = PCU_Time();
  double t_init = t0;

  init_elens();

  // TIME
  t1 = PCU_Time();

  e_lens.c_base_curr = c_base_act;

  collect_disc();

  // TIME
  t0 = PCU_Time();
  dt_collect = t0 - t1;

  //print_ent_set(m_act, &e_lens.es_surf);

  // Disc slices, i denotes slice id, -1 denotes no slice adjacency(boundary
  // disc)
  // When regenerating the slice entities, look into the entity slice_map.
  // If it is negative, it is a disc entity. Find the associated slices value.
  // If it is 2*i, it belongs to the old disc, so use create_map[2*i][entity]
  // to find the new corresponding entity. 
  // If it is 2*i+1, it belongs to the new disc, so use
  //  create_map[2*i+1][entity] to find the new entity. 

  //for (int i = 0; i < e_lens.discs.size(); i++) {
  //  std::cout << "Disc " << i << std::endl; 
  //  e_lens.discs.at(i).print_cont();
  //}

  ins_flag = true;
  if(e_lens.discs.at(0).tri.size() > 0) {

    color_discs();
    // TIME
    t1 = PCU_Time();
    dt_color_disc = t1 - t0;

    std::cout << "ext_slice[0] " << ext_slice[0]
              << " ext_slice[1] " << ext_slice[1]
              << " ext_pc[0] " << ext_pc[0]
              << std::endl;

    crt_new_vert_list();

    // Collect the affected cells.
    collect_cell();
    coll_proj_map();

    calc_ctr_vel(false);
    //if(calc_ext and e_lens.c_base_curr->get_cell_ext_gmi(0, cell_id))
      fix_vel_new();

    bool skip_ins = false;
/*
    if(detect_inv()) {

      std::cout << "disc inversion, trying to fix" << std::endl;
      vd_disc_cut vd_cut;
      if(!do_the_cuts(&vd_cut)) {
      std::cout << "ins_flag: disc inverting 1cell, try insertion anyways" << std::endl;
        //skip_ins = true;
        //ins_flag = false;
      }
      else
        std::cout << "fixed" << std::endl;

      init_elens();
      e_lens.c_base_curr = c_base_act;

      //reset_elens_ent();
      reload_edges();
      overwrite_pc(&vd_cut);
      collect_disc();
      color_discs();

      crt_new_vert_list();

      collect_cell();
      coll_proj_map();
    }
*/
    // TIME
    t0 = PCU_Time();
    dt_reload_inv = t0 - t1;

    crt_new_cell();
    crt_new_vert();

    asgn_new_vert_ext();

    apf::MeshEntity* v_edge[3];
    v_edge[0] = e_lens.vert_ctr_new.at(0);
    v_edge[1] = e_lens.vert_sp_new.at(0);
    v_edge[2] = e_lens.vert_ctr_new.at(1);

    //std::cout << "Old center vertex " << e_lens.vert_ctr << std::endl;
    //std::cout << "Edge vertices " << e_lens.vert_ctr_new[0] << " "
    //                              << e_lens.vert_sp_new[0] << " " 
    //                              << e_lens.vert_ctr_new[1] << " "
    //          << std::endl;
    apf::ModelEntity* edge_ctr_new_em = m_act->findModelEntity(1,
                                                   e_lens.new_cell1_id.at(0));
    e_lens.edge_ctr_new_t.at(0) = buildElement(m_act, edge_ctr_new_em, 
                                              apf::Mesh::EDGE, v_edge, 0);
    e_lens.edge_ctr_new_b.at(0) = buildElement(m_act, edge_ctr_new_em, 
                                              apf::Mesh::EDGE, &v_edge[1], 0);

    int mdl_type = m_act->getModelType(edge_ctr_new_em);
    int mdl_tag = m_act->getModelTag(edge_ctr_new_em);

    std::cout << "Edge " << e_lens.edge_ctr_new_t.at(0) << " "
              << mdl_type << "c" << mdl_tag << "\n"
              << "Edge " << e_lens.edge_ctr_new_b.at(0) << " "
              << mdl_type << "c" << mdl_tag
              << std::endl; 

    //vd_print_vert(m_act, e_lens.edge_ctr_new_t.at(0));
    //vd_print_vert(m_act, e_lens.edge_ctr_new_b.at(0));

    // Check the motion of the new vertices.
    apf::Field* vel_field = e_lens.m->findField("velocity_field");

    apf::Vector3 point;
    apf::Vector3 dn;

    e_lens.m->getPoint(vert_ctr_act, 0, e_lens.pos);
    
    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0,  pos_old);
    }
    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_old);
    }
    e_lens.vt.at(0) = apf::Vector3(0,0,0);
    e_lens.vt.at(1) = apf::Vector3(0,0,0);
    e_lens.split_vt.at(0) = apf::Vector3(0,0,0);

    // Expand the disc:
    expand_lens();
    e_lens.m->acceptChanges();

    fill_lens();
    fill_lens_elem();

    // TIME
    t1 = PCU_Time();
    dt_crt_ent = t1 - t0;

    e_lens.m->acceptChanges();
    // Recreate the slice entities:
    recreate_slices();

    e_lens.m->acceptChanges();
    // Destroy the old entities:
    destroy_ent();

    f_calc.refresh_mesh();

    f_calc.vd_att_fields(e_lens.m, e_lens.vert_ctr_new.at(0));
    f_calc.vd_att_fields(e_lens.m, e_lens.vert_ctr_new.at(1));
    f_calc.vd_att_fields(e_lens.m, e_lens.vert_sp_new.at(0));

    e_lens.m->acceptChanges();

    get_ent_set_new();

    // TIME
    t0 = PCU_Time();
    dt_calc_vel_init = t0 - t1;

    upd_cell();

    //if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL)
    //  update_shell_gmi();
    //if(shell_ext) {
    //  update_shell_pos();
    //}

    find_ctr_vel();
    fix_vel_new();

    find_pos_non_inv_init();
    //find_pos_non_inv();

    if(m_act == m_trial) {
      int trial_temp = trial_type;
      trial_type = 113;
      //vtk_trial();
      vtk_mesh();
      trial_type = trial_temp;
    }

    if(vd_chk_neg_sgn(e_lens.m) > 0) {
      skip_ins = true;
      std::cout << "ins_flag: neg tet, skip_ins" << std::endl;
    }

    if(skip_ins) {
      ins_flag = false;
      std::cout << "ins_flag: skip_ins" << std::endl;
      for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_old);
      }
      for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
        e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_old);
      }
    }

    // TIME
    t1 = PCU_Time();
    dt_upd_pos_init = t1 - t0;

    std::cout << "pow_dis(" << trial_curr + 1 << ") = " 
              << calc_energy_diss_rate_sing() << ";" << std::endl;
    if(ins_flag) {
      find_ext_dir();
    }

    // TIME
    t0 = PCU_Time();
    dt_relax = t0 - t1;

    if(m_act == m_trial) {
      int trial_temp = trial_type;
      trial_type = 11;
      //vtk_trial();
      vtk_mesh();
      trial_type = trial_temp;
    }
    else {
      int trial_temp = trial_type;
      trial_type = 111;
      //vtk_trial();
      vtk_mesh();
      trial_type = trial_temp;
    }

    // Check for negative elements.
    if(!ins_flag or vd_chk_neg_sgn(e_lens.m) > 0) {
      for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_old);
      }
      for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
         e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_old);
      }
      ins_flag = false;
      if(e_lens.m == m_trial and !skip_en)
        w1.at(trial_curr) = 0.1;
      std::cout << "ins_flag: neg element" << std::endl;
    }
    else {
      ins_flag = true;
    }
    //vtk_trial();
    vtk_mesh();

    e_lens.m->acceptChanges();
    //e_lens.m->verify();
  }
  else {
    std::cout << "Disc is empty!" << std::endl;
    disc_load = false;
    ins_flag = false;
    if(e_lens.m == m_trial and !skip_en)
      w1.at(trial_curr) = 0.1;
    std::cout << "ins_flag: empty disc" << std::endl;
  }
  // TIME
  double t_final = PCU_Time();
  dt_tot = t_final - t_init;
  std::cout << " dt_collect " << dt_collect
            << " dt_color_disc " << dt_color_disc
            << " dt_reload_inv " << dt_reload_inv
            << " dt_crt_ent " << dt_crt_ent
            << " dt_calc_vel_init " << dt_calc_vel_init
            << " dt_upd_pos_init " << dt_upd_pos_init
            << " dt_relax " << dt_relax
            << " dt_vtk_rest " << dt_vtk_rest
            << " dt_tot " << dt_tot
            << std::endl;
}

// Not very useful, as the inversion could happen after generation of new 
// entities and fixing the velocities.
void vd_edisc::assert_inv_pre() {
  apf::Up up;
  apf::Downward d_v;
  apf::Downward d_t;

  apf::Vector3 a_curr;
  apf::Vector3 v_pos;
  apf::Vector3 t_pos;

  std::cout << "Detect_inv." << std::endl;

  for (int i = 0; i < e_lens.discs.size(); i++) {
    for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {
      apf::MeshEntity* e_curr = e_lens.discs.at(i).tri.at(j);
      e_lens.m->getUp(e_curr, up);

      a_curr = vd_area_out(e_lens.m, e_curr);
      assert(a_curr.getLength() > std::numeric_limits<double>::epsilon());
      t_pos = getLinearCentroid(e_lens.m, e_curr);

      std::cout << "t_pos " << t_pos << " a_curr " << a_curr << std::endl;

      for (int k = 0; k < up.n; k++) {
        e_lens.m->getDownward(up.e[k], 0, d_v);
        e_lens.m->getDownward(up.e[k], 2, d_t);

        int j_t = findIn(d_t, 4, e_curr);
        int j_v = lookup_tet_x_surf[j_t];

        e_lens.m->getPoint(d_v[j_v], 0, v_pos);

        std::cout << "v_pos " << v_pos << std::endl;

        v_pos = v_pos - t_pos;
        assert(v_pos.getLength() > std::numeric_limits<double>::epsilon());

        double dir = v_pos*a_curr;
        if(dir < 0)
          a_curr = a_curr*(-1);

        int slice_id = e_lens.slice_map[up.e[k]] - 1;

        dir = e_lens.vt.at(slice_id)*a_curr;
        std::cout << "v_pos " << v_pos << " a_curr " << a_curr 
                  << " vt " << e_lens.vt.at(slice_id) << " dir " << dir
                  << std::endl;

        assert(!(dir < std::numeric_limits<double>::epsilon()) );
      }
    }
  }
}


// This also uses the old entities...
bool vd_edisc::detect_inv() {
  std::cout << "Checking for inversion " << std::endl;
  for (int i = 0; i < e_lens.discs.size(); i++) {
    int v1 = e_lens.slices.at(2*i)-1;
    int v2 = e_lens.slices.at(2*i+1)-1;
    apf::Vector3 temp1(0,0,0);
    apf::Vector3 temp2(0,0,0);

    //e_lens.m->getPoint(e_lens.vert_ctr_new.at(v1), 0, temp1);
    //e_lens.m->getPoint(e_lens.vert_ctr_new.at(v2), 0, temp2);
    //std::cout << i << " " 
    //          << e_lens.vt.at(v1) << " "  << e_lens.vt.at(v2) << std::endl;

    apf::Vector3 mot_rel(0,0,0);
    mot_rel = e_lens.vt.at(v1) - e_lens.vt.at(v2);

    //if(detect_inv(i, temp1, e_lens.vt.at(v1))) {
    //if(detect_inv(i, temp1, mot_rel)) {
    //  //std::cout << "Inverting " << e_lens.vt.at(v1) << std::endl;
    //  std::cout << "Inverting " << mot_rel << std::endl;
    //  return true;
    //}

    //if(detect_inv(i, temp2, e_lens.vt.at(v2))) {
    if(detect_inv(i, pos_old, mot_rel)) {
      //std::cout << "Inverting " << e_lens.vt.at(v2) << std::endl;
      std::cout << "Inverting " << mot_rel << std::endl;
      return true;
    }
  }
  return false;
}

bool vd_edisc::detect_inv(int disc_id, apf::Vector3 v_pos, 
                                            apf::Vector3 mot_dir) {

  apf::Up up;

  apf::Vector3 norm_p(0,0,0);
  apf::Vector3 rel_dir1(0,0,0);
  apf::Vector3 rel_dir2(0,0,0);

  // Detect flip for the disc, top and bottom.
  for (int j = 0; j < e_lens.discs.at(disc_id).edge.size(); j++) {
    apf::MeshEntity* e_curr = e_lens.discs.at(disc_id).edge.at(j);
    // Edge is not exterior:
    if(e_lens.discs.at(disc_id).e_t_map1[e_curr] != 0 and 
       e_lens.discs.at(disc_id).e_t_map2[e_curr] != 0) {
      int t_1 = e_lens.discs.at(disc_id).e_t_map1[e_curr];
      int t_2 = e_lens.discs.at(disc_id).e_t_map2[e_curr];
      apf::MeshEntity* tri_1 = e_lens.discs.at(disc_id).tri.at(t_1 - 1);
      apf::MeshEntity* tri_2 = e_lens.discs.at(disc_id).tri.at(t_2 - 1);

      int e_1;
      int e_2;

      if(e_lens.discs.at(disc_id).t_e_map1[tri_1] != j+1) {
        e_1 = e_lens.discs.at(disc_id).t_e_map1[tri_1] - 1;
      }
      else {
        assert(e_lens.discs.at(disc_id).t_e_map1[tri_1] == j+1);
        e_1 = e_lens.discs.at(disc_id).t_e_map2[tri_1] - 1;
      }
      if(e_lens.discs.at(disc_id).t_e_map1[tri_2] != j+1) {
        e_2 = e_lens.discs.at(disc_id).t_e_map1[tri_2] - 1;
      }
      else {
        assert(e_lens.discs.at(disc_id).t_e_map1[tri_2] == j+1);
        e_2 = e_lens.discs.at(disc_id).t_e_map2[tri_2] - 1;
      }
      // Relative directions of the edges preceeding and succeeding the current 
      // edge on the disc.
      rel_dir1 = e_lens.discs.at(disc_id).v_pos.at(e_1) - v_pos;
      rel_dir2 = e_lens.discs.at(disc_id).v_pos.at(e_2) - v_pos;
      rel_dir1 = norm_0(rel_dir1);
      rel_dir2 = norm_0(rel_dir2);
      apf::Vector3 j_dir = e_lens.discs.at(disc_id).v_pos.at(j) - v_pos;
      j_dir = norm_0(j_dir);
      mot_dir = norm_0(mot_dir);
      //if(mot_dir.getLength() < std::numeric_limits<double>::epsilon())
      norm_p = cross(j_dir, mot_dir);
      norm_p = norm_0(norm_p);
      double dir;
      double len = mot_dir.getLength();
      if(!(len < std::numeric_limits<double>::epsilon()) )
        dir = (norm_p*rel_dir1)*(norm_p*rel_dir2);
      else
        dir = -1;
      if(len > std::numeric_limits<double>::epsilon() and 
          (std::fabs(mot_dir*rel_dir1) < exp_dir_cos_th or
         std::fabs(mot_dir*rel_dir2) < exp_dir_cos_th )) {
        std::cout << "Motion inplane with one of the triangles. "
                  << "Edge " << j << " norm_p " << norm_p
                  << " rel_dir1 " << rel_dir1
                  << " rel_dir2 " << rel_dir2
                  << " dir " << dir
                  << std::endl;
        return true;
      }

      if(dir > std::numeric_limits<double>::min()) {
        std::cout << "Edge " << j << " norm_p " << norm_p
                  << " rel_dir1 " << rel_dir1
                  << " rel_dir2 " << rel_dir2
                  << " dir " << dir
                  << std::endl;
        return true;
      }
    }
  }

  return false;
}

void vd_edisc::expand_lens() {
  for (int i = 0; i < e_lens.discs.size(); i++) {
    // Center vertex of the old disc.
    apf::MeshEntity* v_d_cp = e_lens.create_map[2*i][vert_ctr_act];
    // Create the edges, old disc
    apf::Downward down;

    //std::cout << "Top disc, or the disc with the old 0cell, vertex " 
    //          << v_d_cp << std::endl;

    for (int j = 0; j < e_lens.discs.at(i).edge.size(); j++) {
      e_lens.m->getDownward(e_lens.discs.at(i).edge.at(j), 0, down);
      //if (down[0] == vert_ctr_act)
      //  down[0] = v_d_cp;
      //else
      //  down[1] = v_d_cp;

      int j_v = findIn(down, 2, vert_ctr_act);
      assert(j_v > -1);
      down[j_v] = v_d_cp;

      apf::ModelEntity* mdl;
      if(trial_type == 3)
        mdl = e_lens.m->toModel(e_lens.discs.at(i).edge.at(j));
      else
        mdl = e_lens.discs.at(i).e_em.at(j);

      apf::MeshEntity* edge_new = buildElement(e_lens.m, mdl, 
                                                  apf::Mesh::EDGE, down, 0);
      e_lens.create_map[2*i][e_lens.discs.at(i).edge.at(j)] = edge_new;
      //std::cout << "old edge " << e_lens.discs.at(i).edge.at(j) << " " 
      //          << e_lens.m->getModelType(mdl) << "c"
      //          << e_lens.m->getModelTag(mdl)
      //          << " new edge " << edge_new << std::endl;
      //vd_print_vert(e_lens.m, e_lens.discs.at(i).edge.at(j));
      //vd_print_vert(e_lens.m, edge_new);
    }

    // Create the triangles, old disc
    for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {
      e_lens.m->getDownward(e_lens.discs.at(i).tri.at(j), 0, down);
      //for (int j_v = 0; j_v < 3; j_v++) {
      //  if (down[j_v] == vert_ctr_act)
      //    down[j_v] = v_d_cp;
      //}
      int j_v = findIn(down, 3, vert_ctr_act);
      assert(j_v > -1);
      down[j_v] = v_d_cp;

/*
      e_lens.m->getDownward(e_lens.discs.at(i).tri.at(j), 1, down);
      for (int k = 0; k < 3; k++) {
        std::vector<apf::MeshEntity*>::iterator it_e;
        it_e = std::find(e_lens.discs.at(i).edge.begin(), 
                         e_lens.discs.at(i).edge.end(), 
                        down[k]);
        // If on the disc, replace it with the new edge for the top disc.
        if(it_e != e_lens.discs.at(i).edge.end())
          down[k] = e_lens.create_map[2*i][down[k]];
      }
*/
      apf::ModelEntity* mdl;
      if(trial_type == 3)
        mdl = e_lens.m->toModel(e_lens.discs.at(i).tri.at(j));
      else
        mdl = e_lens.discs.at(i).t_em.at(j);

      apf::MeshEntity* tri_new = buildElement(e_lens.m, mdl, 
                                              apf::Mesh::TRIANGLE, down, 0);
      e_lens.create_map[2*i][e_lens.discs.at(i).tri.at(j)] = tri_new;

      //std::cout << "old tri " << e_lens.discs.at(i).tri.at(j) << " " 
      //          << e_lens.m->getModelType(mdl) << "c"
      //          << e_lens.m->getModelTag(mdl)
      //          << " new tri " << tri_new << std::endl;

      //vd_print_vert(e_lens.m, e_lens.discs.at(i).tri.at(j));
      //vd_print_vert(e_lens.m, tri_new);
    }

    v_d_cp = e_lens.create_map[2*i+1][vert_ctr_act];

    //std::cout << "Bottom disc, or the disc with the new 0cell, vertex " 
    //          << v_d_cp << std::endl;

    // Create the edges, new disc
    for (int j = 0; j < e_lens.discs.at(i).edge.size(); j++) {
      e_lens.m->getDownward(e_lens.discs.at(i).edge.at(j), 0, down);
      //if (down[0] == vert_ctr_act)
      //  down[0] = v_d_cp;
      //else
      //  down[1] = v_d_cp;

      int j_v = findIn(down, 2, vert_ctr_act);
      assert(j_v > -1);
      down[j_v] = v_d_cp;

      apf::ModelEntity* mdl = e_lens.discs.at(i).t2_em.at(j);
      apf::MeshEntity* edge_new = buildElement(e_lens.m, mdl, 
                                                  apf::Mesh::EDGE, down, 0);
      e_lens.create_map[2*i+1][e_lens.discs.at(i).edge.at(j)] = edge_new;
      //std::cout << "old edge " << e_lens.discs.at(i).edge.at(j) << " " 
      //          << e_lens.m->getModelType(mdl) << "c"
      //          << e_lens.m->getModelTag(mdl)
      //          << " new edge " << edge_new << std::endl;

      //vd_print_vert(e_lens.m, e_lens.discs.at(i).edge.at(j));
      //vd_print_vert(e_lens.m, edge_new);
    }

    // Create the triangles, new disc
    for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {

      e_lens.m->getDownward(e_lens.discs.at(i).tri.at(j), 0, down);
      //for (int j_v = 0; j_v < 3; j_v++) {
      //  if (down[j_v] == vert_ctr_act)
      //    down[j_v] = v_d_cp;
      //}
      int j_v = findIn(down, 3, vert_ctr_act);
      assert(j_v > -1);
      down[j_v] = v_d_cp;
/*
      e_lens.m->getDownward(e_lens.discs.at(i).tri.at(j), 1, down);
      for (int k = 0; k < 3; k++) {
        std::vector<apf::MeshEntity*>::iterator it_e;
        it_e = std::find(e_lens.discs.at(i).edge.begin(),
                         e_lens.discs.at(i).edge.end(),
                        down[k]);
        if(it_e != e_lens.discs.at(i).edge.end())
          down[k] = e_lens.create_map[2*i+1][down[k]];
      }
*/
      apf::ModelEntity* mdl = e_lens.discs.at(i).t_em.at(j);
      apf::MeshEntity* tri_new = buildElement(e_lens.m, mdl, 
                                              apf::Mesh::TRIANGLE, down, 0);
      e_lens.create_map[2*i+1][e_lens.discs.at(i).tri.at(j)] = tri_new;

      //std::cout << "old tri " << e_lens.discs.at(i).tri.at(j) << " " 
      //          << e_lens.m->getModelType(mdl) << "c"
      //          << e_lens.m->getModelTag(mdl)
      //          << " new tri " << tri_new << std::endl;
      //vd_print_vert(e_lens.m, e_lens.discs.at(i).tri.at(j));
      //vd_print_vert(e_lens.m, tri_new);
    }

  }

}

// Fill the insides of the expanded lenses.
void vd_edisc::fill_lens() {
  apf::Downward down;
  apf::Downward down_e;
  apf::MeshEntity* e_temp;

  std::vector<apf::MeshEntity* > tri_new(0);
  for (int i = 0; i < e_lens.discs.size(); i++) {

    int v1 = e_lens.slices.at(2*i)-1;
    int v2 = e_lens.slices.at(2*i+1)-1;
    apf::MeshEntity* vert_1 = e_lens.vert_ctr_new.at(v1);
    apf::MeshEntity* vert_2 = e_lens.vert_ctr_new.at(v2);

    // Recreate the edges using the vertex joining the two lenses.
    // The central disc edge and triangles.
    std::vector<apf::MeshEntity* > edge_d_new(e_lens.discs.at(i).edge.size());
    std::vector<apf::MeshEntity* > tri_d_new(e_lens.discs.at(i).tri.size());

    // The central disc.
    for (int j = 0; j < e_lens.discs.at(i).edge.size(); j++) {
      e_lens.m->getDownward(e_lens.discs.at(i).edge.at(j), 0, down);
      if (down[0] == vert_ctr_act)
        down[0] = e_lens.vert_sp_new.at(i);
      else
        down[1] = e_lens.vert_sp_new.at(i);

      apf::ModelEntity* mdl = e_lens.discs.at(i).e_em.at(j);
      edge_d_new.at(j) = buildElement(e_lens.m, mdl, 
                                            apf::Mesh::EDGE, down, 0);

      //std::cout << "old edge " << e_lens.discs.at(i).edge.at(j) << " " 
      //          << e_lens.m->getModelType(mdl) << "c"
      //          << e_lens.m->getModelTag(mdl)
      //          << " new edge " << edge_d_new.at(j) << std::endl;

      //vd_print_vert(e_lens.m, e_lens.discs.at(i).edge.at(j));
      //vd_print_vert(e_lens.m, edge_d_new.at(j));
      e_lens.split_map[i][e_lens.discs.at(i).edge.at(j)] = edge_d_new.at(j);
    }

    // Create the triangles, central disc
    for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {
/*
      e_lens.m->getDownward(e_lens.discs.at(i).tri.at(j), 1, down);
      for (int k = 0; k < 3; k++) {
        std::vector<apf::MeshEntity*>::iterator it_e;
        it_e = std::find(e_lens.discs.at(i).edge.begin(), 
                         e_lens.discs.at(i).edge.end(), 
                        down[k]);
        // If on the disc, replace it with the new edge for the central disc.
        if(it_e != e_lens.discs.at(i).edge.end()) {
          int e1 = std::distance(e_lens.discs.at(i).edge.begin(), it_e);
          down[k] = edge_d_new.at(e1);
        }
      }
*/
      e_lens.m->getDownward(e_lens.discs.at(i).tri.at(j), 0, down);
      for (int k = 0; k < 3; k++) {
        if (down[k] == vert_ctr_act)
          down[k] = e_lens.vert_sp_new.at(i);
      }

      apf::ModelEntity* mdl = e_lens.discs.at(i).t_em.at(j);
      tri_d_new.at(j) = buildElement(e_lens.m, mdl, 
                                          apf::Mesh::TRIANGLE, down, 0);

      e_lens.split_map.at(i)[e_lens.discs.at(i).tri.at(j)] = tri_d_new.at(j);
      //std::cout << "old tri " << e_lens.discs.at(i).tri.at(j) << " " 
      //          << e_lens.m->getModelType(mdl) << "c"
      //          << e_lens.m->getModelTag(mdl)
      //          << " disc tri " << tri_d_new.at(j) << std::endl;

      //std::cout << "old tri " << e_lens.discs.at(i).tri.at(j) 
      //          << " new tri " <<  << std::endl;

      //vd_print_vert(e_lens.m, e_lens.discs.at(i).tri.at(j));
      //vd_print_vert(e_lens.m, tri_d_new.at(j));
    }

    //std::cout << "Disc " << i << std::endl;
    tri_new.resize(e_lens.discs.at(i).edge.size());

    // The top.
    // The inside triangles. TODO problem here, edge model is not compliant
    // with 3cell edge insertion.
    //std::cout << "Filling Inside " << std::endl;
    for (int j = 0; j < e_lens.discs.at(i).edge.size(); j++) {

      e_lens.m->getDownward(e_lens.discs.at(i).edge.at(j), 0, down);
      if (down[0] == vert_ctr_act)
        down[0] = e_lens.vert_sp_new.at(i);
      else
        down[1] = e_lens.vert_sp_new.at(i);

      down[2] = vert_1;

      apf::ModelEntity* mdl;
      if(trial_type == 3) {
        mdl = vd_3c->get_mdl(3, trial_curr);
        //mdl = e_lens.m->findModelEntity(trial_type, trial_curr);
      }
      else {
        mdl = e_lens.m->toModel(e_lens.discs.at(i).edge.at(j));
      }
      tri_new.at(j) = buildElement(e_lens.m, mdl, 
                                             apf::Mesh::TRIANGLE, down, 0);

      e_lens.tri_map[2*i][e_lens.discs.at(i).edge.at(j)] = tri_new.at(j);
      //std::cout << "Tri " << tri_new.at(j) 
      //          << e_lens.m->getModelType(mdl) << "c"
      //          << e_lens.m->getModelTag(mdl)
      //          << std::endl;

      //vd_print_vert(e_lens.m, tri_new.at(j));

      //apf::Downward vert; 
      //e_lens.m->getDownward(tri_new.at(j), 0, vert);
      //for(int k = 0; k < 3; k++) 
      //  std::cout << vert[k] << ", ";

      //std::cout << std::endl;

    }

    // The bottom.
    // The inside triangles.
    //std::cout << "Filling Inside " << std::endl;
    for (int j = 0; j < e_lens.discs.at(i).edge.size(); j++) {

      e_lens.m->getDownward(e_lens.discs.at(i).edge.at(j), 0, down);
      if (down[0] == vert_ctr_act)
        down[0] = e_lens.vert_sp_new.at(i);
      else
        down[1] = e_lens.vert_sp_new.at(i);

      down[2] = vert_2;

      apf::ModelEntity* mdl;
      if(trial_type == 3) {
        mdl = vd_3c->get_mdl(3, trial_curr);
        //mdl = e_lens.m->findModelEntity(trial_type, trial_curr);
      }
      else {
        mdl = e_lens.m->toModel(e_lens.discs.at(i).edge.at(j));
      }
      tri_new.at(j) = buildElement(e_lens.m, mdl, 
                                              apf::Mesh::TRIANGLE, down, 0);

      e_lens.tri_map[2*i+1][e_lens.discs.at(i).edge.at(j)] = tri_new.at(j);
      //std::cout << "Tri " << tri_new.at(j) 
      //          << e_lens.m->getModelType(mdl) << "c"
      //          << e_lens.m->getModelTag(mdl)
      //          << std::endl;

      //vd_print_vert(e_lens.m, tri_new.at(j));

      //apf::Downward vert; 
      //e_lens.m->getDownward(tri_new.at(j), 0, vert);
      //for(int k = 0; k < 3; k++) 
      //  std::cout << vert[k] << ", ";

      //std::cout << std::endl;

    }
  }
}

apf::MeshEntity* vd_edisc::recreate_inv_elem(apf::MeshEntity* tet, 
                                              apf::ModelEntity* mdl) {
  apf::Downward down;
  double meas;
  apf::MeshEntity* e_temp;
  e_lens.m->getDownward(tet, 0, down);

  apf::MeshElement* ee = createMeshElement(e_lens.m, tet);
  meas = measure(ee);
  destroyMeshElement(ee);

  if(meas < 0. 
     and std::fabs(meas) > std::numeric_limits<double>::min()) {
  //if(vd_volume_tet(e_lens.m, down) < 0) {
    e_lens.m->destroy(tet);
    e_temp = down[1];
    down[1] = down[2];
    down[2] = e_temp;
    tet = buildElement(e_lens.m, mdl, apf::Mesh::TET, down, 0);
    e_lens.m->getDownward(tet, 0, down);

    ee = createMeshElement(e_lens.m, tet);
    meas = measure(ee);
    destroyMeshElement(ee);
    if(std::fabs(meas) > std::numeric_limits<double>::min())
      assert(meas > 0);

    //assert(vd_volume_tet(e_lens.m, down) > 0);
    //vd_print_vert(e_lens.m, e_void);
  }
  return tet;
}

//The adjacencies of some of the newly generated tets inside the lens seem to be false or incomplete.. Some of these tetrahedra are created from triangles adjacent to the top and bottom 3c edges and one of the void elements touching the same edge

// Fill the insides of the 3cell lenses. Used in insert_3celledge.
// TODO the first part seems superfluous: Just use the first upper adjacency for
// the ordering, shift the ordering and replace the vertex for the second one.
// Extends to disc triangles with one upper adjacency.
void vd_edisc::fill_3c_lens_elem() {
  apf::Downward down;
  apf::Downward d_tri;

  for (int i = 0; i < e_lens.discs.size(); i++) {

    int v1 = e_lens.slices.at(2*i)-1;
    int v2 = e_lens.slices.at(2*i+1)-1;
    apf::MeshEntity* vert_1 = e_lens.vert_ctr_new.at(v1);
    apf::MeshEntity* vert_2 = e_lens.vert_ctr_new.at(v2);

    // The inside tetrahedra. The vertex ordering is taken from upward adjacency
    // tets. Replace the old vert_ctr with vert_sp_new and the slice vertex with
    // the corresponding vert_ctr_new.
    apf::Up up;
    for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {
      apf::MeshEntity* tri_curr = e_lens.discs.at(i).tri.at(j);
      e_lens.m->getUp(tri_curr, up);
      if(up.n == 2) {
        assert(e_lens.slice_map[up.e[0]] != e_lens.slice_map[up.e[1]]);
        // Top tet:
        e_lens.m->getDownward(up.e[0], 0, down);
        e_lens.m->getDownward(up.e[0], 2, d_tri);
        int j_tri = findIn(d_tri, 4, tri_curr);

        // Replace the vertices:
        int v_slice_id = lookup_tet_x_surf[j_tri];
        int j_v = findIn(down, 4, vert_ctr_act);
        assert(j_v > -1);

        int slice_curr = e_lens.slice_map[up.e[0]];
        assert(slice_curr > 0);

        // TODO For some reason this is not the correct ordering?
        down[j_v] = e_lens.vert_sp_new.at(i);
        int i_shift = 0;
        if(slice_curr == e_lens.slices.at(2*i)) {
          down[v_slice_id] = vert_1;
        }
        else {
          assert(slice_curr == e_lens.slices.at(2*i+1));
          down[v_slice_id] = vert_2;
          i_shift = 1;
        }

        apf::ModelEntity* mdl;
        mdl = vd_3c->get_mdl(3, trial_curr);

        apf::MeshEntity* elem_new;
        elem_new = buildElement(e_lens.m, mdl, apf::Mesh::TET, down, 0);

        e_lens.tri_tet_map.at(2*i+i_shift)[e_lens.discs.at(i).tri.at(j)] = elem_new;

        // Bottom tet:
        e_lens.m->getDownward(up.e[1], 0, down);
        e_lens.m->getDownward(up.e[1], 2, d_tri);
        j_tri = findIn(d_tri, 4, tri_curr);

        // Replace the vertices:
        v_slice_id = lookup_tet_x_surf[j_tri];
        j_v = findIn(down, 4, vert_ctr_act);
        assert(j_v > -1);

        slice_curr = e_lens.slice_map[up.e[1]];
        assert(slice_curr > 0);

        down[j_v] = e_lens.vert_sp_new.at(i);

        i_shift = 0;
        if(slice_curr == e_lens.slices.at(2*i)) {
          down[v_slice_id] = vert_1;
        }
        else {
          assert(slice_curr == e_lens.slices.at(2*i+1));
          down[v_slice_id] = vert_2;
          i_shift = 1;
        }

        mdl = vd_3c->get_mdl(3, trial_curr);

        elem_new = buildElement(e_lens.m, mdl, apf::Mesh::TET, down, 0);
        e_lens.tri_tet_map.at(2*i+i_shift)[e_lens.discs.at(i).tri.at(j)] = elem_new;
      }
      else {
        assert(up.n == 1);
        // Top tet:
        e_lens.m->getDownward(up.e[0], 0, down);
        e_lens.m->getDownward(up.e[0], 2, d_tri);
        int j_tri = findIn(d_tri, 4, tri_curr);

        // Replace the vertices:
        int v_slice_id = lookup_tet_x_surf[j_tri];
        int j_v = findIn(down, 4, vert_ctr_act);
        assert(j_v > -1);

        int slice_curr = e_lens.slice_map[up.e[0]];
        assert(slice_curr > 0);

        // TODO For some reason this is not the correct ordering?
        down[j_v] = e_lens.vert_sp_new.at(i);
        int i_shift = 0;
        if(slice_curr == e_lens.slices.at(2*i)) {
          down[v_slice_id] = vert_1;
        }
        else {
          assert(slice_curr == e_lens.slices.at(2*i+1));
          down[v_slice_id] = vert_2;
          i_shift = 1;
        }

        apf::ModelEntity* mdl;
        mdl = vd_3c->get_mdl(3, trial_curr);

        apf::MeshEntity* elem_new;
        elem_new = buildElement(e_lens.m, mdl, apf::Mesh::TET, down, 0);

        e_lens.tri_tet_map.at(2*i+i_shift)[e_lens.discs.at(i).tri.at(j)] = elem_new;
        // Bottom tet, reverse the vertex order and replace the new center vertex:
        down[v_slice_id] = e_lens.vert_sp_new.at(i);
        if(i_shift == 0) {
          down[j_v] = vert_2;
        }
        else {
          down[j_v] = vert_1;
        }

        mdl = vd_3c->get_mdl(3, trial_curr);
        elem_new = buildElement(e_lens.m, mdl, apf::Mesh::TET, down, 0);

        e_lens.tri_tet_map.at(2*i+i_shift)[e_lens.discs.at(i).tri.at(j)] = elem_new;


      }

    }
  }
}

// Fill the insides of the expanded lenses.
void vd_edisc::fill_lens_elem() {
  apf::Downward down;
  apf::Downward d_tri;

  for (int i = 0; i < e_lens.discs.size(); i++) {

    int v1 = e_lens.slices.at(2*i)-1;
    int v2 = e_lens.slices.at(2*i+1)-1;
    apf::MeshEntity* vert_1 = e_lens.vert_ctr_new.at(v1);
    apf::MeshEntity* vert_2 = e_lens.vert_ctr_new.at(v2);

    // The inside tetrahedra. The vertex ordering is taken from upward adjacency
    // tets. Replace the old vert_ctr with vert_sp_new and the slice vertex with
    // the corresponding vert_ctr_new.
    apf::Up up;
    for (int j = 0; j < e_lens.discs.at(i).tri.size(); j++) {
      apf::MeshEntity* tri_curr = e_lens.discs.at(i).tri.at(j);
      e_lens.m->getUp(tri_curr, up);
      assert(up.n == 2);
      assert(e_lens.slice_map[up.e[0]] != e_lens.slice_map[up.e[1]]);

      // Top tet:
      e_lens.m->getDownward(up.e[0], 0, down);
      e_lens.m->getDownward(up.e[0], 2, d_tri);
      int j_tri = findIn(d_tri, 4, tri_curr);

      // Replace the vertices:
      int v_slice_id = lookup_tet_x_surf[j_tri];
      int j_v = findIn(down, 4, vert_ctr_act);
      assert(j_v > -1);

      int slice_curr = e_lens.slice_map[up.e[0]];
      assert(slice_curr > 0);

      // TODO For some reason this is not the correct ordering?
      down[j_v] = e_lens.vert_sp_new.at(i);
      int i_shift = 0;
      if(slice_curr == e_lens.slices.at(2*i)) {
        down[v_slice_id] = vert_1;
      }
      else {
        assert(slice_curr == e_lens.slices.at(2*i+1));
        down[v_slice_id] = vert_2;
        i_shift = 1;
      }

      apf::ModelEntity* mdl;
      mdl = e_lens.m->toModel(e_lens.discs.at(i).tri.at(j));

      apf::MeshEntity* elem_new;
      elem_new = buildElement(e_lens.m, mdl, apf::Mesh::TET, down, 0);

      e_lens.tri_tet_map.at(2*i+i_shift)[e_lens.discs.at(i).tri.at(j)] = elem_new;

      // Sanity check for the downward triangle adjacencies:
/*
      e_lens.m->getDownward(elem_new, 2, d_tri);
      int i_edge = e_lens.discs.at(i).t_e_map1[e_lens.discs.at(i).tri.at(j)] - 1;
      apf::MeshEntity* e_adj = e_lens.discs.at(i).edge.at(i_edge);
      e_adj = e_lens.tri_map[2*i+i_shift][e_adj];
      int i_san = findIn(d_tri, 4, e_adj);
      assert(i_san > -1);

      i_edge = e_lens.discs.at(i).t_e_map2[e_lens.discs.at(i).tri.at(j)] - 1;
      e_adj = e_lens.discs.at(i).edge.at(i_edge);
      e_adj = e_lens.tri_map[2*i+i_shift][e_adj];
      i_san = findIn(d_tri, 4, e_adj);
      assert(i_san > -1);

      e_adj = e_lens.split_map.at(i)[e_lens.discs.at(i).tri.at(j)];
      i_san = findIn(d_tri, 4, e_adj);
      assert(i_san > -1);
      
      e_adj = e_lens.create_map.at(2*i+i_shift)[e_lens.discs.at(i).tri.at(j)];
      i_san = findIn(d_tri, 4, e_adj);
      assert(i_san > -1);
*/
      // Bottom tet:
      e_lens.m->getDownward(up.e[1], 0, down);
      e_lens.m->getDownward(up.e[1], 2, d_tri);
      j_tri = findIn(d_tri, 4, tri_curr);

      // Replace the vertices:
      v_slice_id = lookup_tet_x_surf[j_tri];
      j_v = findIn(down, 4, vert_ctr_act);
      assert(j_v > -1);

      slice_curr = e_lens.slice_map[up.e[1]];
      assert(slice_curr > 0);

      down[j_v] = e_lens.vert_sp_new.at(i);

      i_shift = 0;
      if(slice_curr == e_lens.slices.at(2*i)) {
        down[v_slice_id] = vert_1;
      }
      else {
        assert(slice_curr == e_lens.slices.at(2*i+1));
        down[v_slice_id] = vert_2;
        i_shift = 1;
      }

      mdl = e_lens.m->toModel(e_lens.discs.at(i).tri.at(j));

      elem_new = buildElement(e_lens.m, mdl, apf::Mesh::TET, down, 0);
      e_lens.tri_tet_map.at(2*i+i_shift)[e_lens.discs.at(i).tri.at(j)] = elem_new;

      // Sanity check for the downward triangle adjacencies:
/*
      e_lens.m->getDownward(elem_new, 2, d_tri);
      i_edge = e_lens.discs.at(i).t_e_map1[e_lens.discs.at(i).tri.at(j)] - 1;
      e_adj = e_lens.discs.at(i).edge.at(i_edge);
      e_adj = e_lens.tri_map[2*i+i_shift][e_adj];
      i_san = findIn(d_tri, 4, e_adj);
      assert(i_san > -1);

      i_edge = e_lens.discs.at(i).t_e_map2[e_lens.discs.at(i).tri.at(j)] - 1;
      e_adj = e_lens.discs.at(i).edge.at(i_edge);
      e_adj = e_lens.tri_map[2*i+i_shift][e_adj];
      i_san = findIn(d_tri, 4, e_adj);
      assert(i_san > -1);

      e_adj = e_lens.split_map.at(i)[e_lens.discs.at(i).tri.at(j)];
      i_san = findIn(d_tri, 4, e_adj);
      assert(i_san > -1);
      
      e_adj = e_lens.create_map.at(2*i+i_shift)[e_lens.discs.at(i).tri.at(j)];
      i_san = findIn(d_tri, 4, e_adj);
      assert(i_san > -1);
*/
    }
  }
}

bool vd_edisc::chk_ma_swap_lens() {
  MaSwap3Dcheck swap_chk(e_lens.m);

  if(trial_type == 1 or trial_type == 3) {
    for(int i = 0; i < e_lens.discs.at(0).edge.size(); i++) {
      apf::MeshEntity* e = e_lens.discs.at(0).edge.at(i);

      apf::MeshEntity* e_new = e_lens.create_map.at(0)[e];
      apf::ModelEntity* mdl = e_lens.m->toModel(e_new);
      if(e_lens.m->getModelType(mdl) == 3) {
        if(!swap_chk.run(e_new)) {
          return false;
        }
      }
      e_new = e_lens.create_map.at(1)[e];
      mdl = e_lens.m->toModel(e_new);
      if(e_lens.m->getModelType(mdl) == 3) {
        if(!swap_chk.run(e_new)) {
          return false;
        }
      }
      e_new = e_lens.split_map.at(0)[e];
      mdl = e_lens.m->toModel(e_new);
      if(e_lens.m->getModelType(mdl) == 3) {
        if(!swap_chk.run(e_new)) {
          return false;
        }
      }
    }
  }
  else {
    assert(trial_type == 2);
    for(int i = 0; i < e_lens.discs.size(); i++) {
      for(int j = 0; j < e_lens.discs.at(i).edge.size(); j++) {
        apf::MeshEntity* e = e_lens.discs.at(0).edge.at(j);
        if(e != e_lens.e_3c_t or e != e_lens.e_3c_t) {
          apf::MeshEntity* e_new = e_lens.create_map.at(2*i)[e];
          apf::ModelEntity* mdl = e_lens.m->toModel(e_new);
          if(e_lens.m->getModelType(mdl) == 3) {
            if(!swap_chk.run(e_new)) {
              return false;
            }
          }
          e_new = e_lens.create_map.at(2*i + 1)[e];
          mdl = e_lens.m->toModel(e_new);
          if(e_lens.m->getModelType(mdl) == 3) {
            if(!swap_chk.run(e_new)) {
              return false;
            }
          }
          e_new = e_lens.split_map.at(i)[e];
          mdl = e_lens.m->toModel(e_new);
          if(e_lens.m->getModelType(mdl) == 3) {
            if(!swap_chk.run(e_new)) {
              return false;
            }
          }
        }
      }
    }
  }
  return true;
}

// Create the 2cell triangles and fill the openings inside the 3cells getting 
// connected.
void vd_edisc::fill_void() {
  // Create a new edge for both connecting 3cells.
  apf::MeshEntity* e_3c_t_new;
  apf::MeshEntity* e_3c_b_new;

  // For each new center vertex, create an edge.
  std::vector<apf::MeshEntity*> e_ctr(0);
  // For each split vertex, create an edge.
  std::vector<apf::MeshEntity*> e_sp(0);

  // The new triangles of the connecting 3-cells.
  std::vector<apf::MeshEntity*> tri_f_ctr(0);
  std::vector<apf::MeshEntity*> tri_e_ctr(0);
  std::vector<apf::MeshEntity*> tri_f_sp(0);
  std::vector<apf::MeshEntity*> tri_e_sp(0);
  // The triangles of the new 2cell.
  std::vector<apf::MeshEntity*> tri_2(0);

  apf::Downward down;
  apf::Up up;

  apf::MeshElement* ee;
  double meas = 0.0;

  apf::ModelEntity* c2c = e_lens.m->findModelEntity(2, e_lens.new_cell2_id);
  std::cout << "New 2c " << e_lens.new_cell2_id << " " 
            << e_lens.m->getModelType(c2c) << "c"
            << e_lens.m->getModelTag(c2c) << std::endl;

  // Create a new vertex on the 2cell center.
  e_lens.v_2c = e_lens.m->createVert(c2c);

  //v_pos = v_pos/(e_lens.vert_ctr_new.size()+e_lens.vert_sp_new.size());
  //std::cout << "Created the 2cell vertex " << v_2c << " at " << v_pos 
  //          << std::endl;
  //e_lens.m->setPoint(e_lens.v_2c, 0, v_pos);
  e_lens.m->setPoint(e_lens.v_2c, 0, e_lens.pos);

  // Create the new 3cell edges if the 3cells are interior.
  apf::MeshEntity* v3_t;
  apf::MeshEntity* v3_b;
  if (e_lens.e_3c_t_flag) {
    apf::ModelEntity* mdl = e_lens.m->toModel(e_lens.e_3c_t);
    e_lens.m->getDownward(e_lens.e_3c_t, 0, down);
    if (down[0] == e_lens.vert_ctr) {
      down[0] = e_lens.v_2c;
      v3_t = down[1];
    }
    else {
      assert(down[1] == e_lens.vert_ctr);
      down[1] = e_lens.v_2c;
      v3_t = down[0];
    }

    e_3c_t_new = buildElement(e_lens.m, mdl, apf::Mesh::EDGE, down, 0);
    //std::cout << "e_3c_t_new " << vd_meas_ent(e_lens.m, e_3c_t_new) 
    //          << std::endl;

    std::cout << "e_3c_t_old " << e_lens.e_3c_t << " " 
                << e_lens.m->getModelType(mdl) << "c"
                << e_lens.m->getModelTag(mdl)
                << " e_3c_t_new " << e_3c_t_new << std::endl;

    //vd_print_vert(e_lens.m, e_3c_t_new);

  }
  if (e_lens.e_3c_b_flag) {
    apf::ModelEntity* mdl = e_lens.m->toModel(e_lens.e_3c_b);
    e_lens.m->getDownward(e_lens.e_3c_b, 0, down);
    if (down[0] == e_lens.vert_ctr) {
      down[0] = e_lens.v_2c;
      v3_b = down[1];
    }
    else {
      assert(down[1] == e_lens.vert_ctr);
      down[1] = e_lens.v_2c;
      v3_b = down[0];
    }
    e_3c_b_new = buildElement(e_lens.m, mdl, apf::Mesh::EDGE, down, 0);

    std::cout << "e_3c_b_old " << e_lens.e_3c_b << " " 
                << e_lens.m->getModelType(mdl) << "c"
                << e_lens.m->getModelTag(mdl)
                << " e_3c_b_new " << e_3c_b_new << std::endl;

    //vd_print_vert(e_lens.m, e_3c_b_new);
  }

  // Create the new edges.
  e_ctr.resize(e_lens.vert_ctr_new.size());
  e_sp.resize(e_lens.vert_sp_new.size());
  tri_f_ctr.resize(e_lens.vert_ctr_new.size());
  tri_e_ctr.resize(e_lens.vert_ctr_new.size());
  tri_f_sp.resize(e_lens.vert_sp_new.size());
  tri_e_sp.resize(e_lens.vert_sp_new.size());

  tri_2.resize(e_lens.discs.size()*2);

  down[0] = e_lens.v_2c;
  //std::cout << "Creating the 2cell edges " << std::endl;

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    down[1] = e_lens.vert_ctr_new.at(i);
    e_ctr.at(i) = buildElement(e_lens.m, c2c, apf::Mesh::EDGE, down, 0);
    std::cout << "e_ctr " << e_ctr[i] << std::endl;
    vd_print_vert(e_lens.m, e_ctr[i]);
    std::cout << "e_ctr(" << i << ") " << vd_meas_ent(e_lens.m, e_ctr.at(i)) 
              << std::endl;
    std::cout << "e_ctr(" << i << ") " << e_ctr.at(i) << " " 
                << e_lens.m->getModelType(c2c) << "c"
                << e_lens.m->getModelTag(c2c) << std::endl;

  }
  for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    down[1] = e_lens.vert_sp_new.at(i);
    e_sp.at(i) = buildElement(e_lens.m, c2c, apf::Mesh::EDGE, down, 0);
    std::cout << "e_sp " << e_sp[i] << std::endl;
    vd_print_vert(e_lens.m, e_sp[i]);
    std::cout << "e_sp(" << i << ") " << vd_meas_ent(e_lens.m, e_sp.at(i)) 
              << std::endl;

    std::cout << "e_sp(" << i << ") " << e_sp.at(i) << " " 
                << e_lens.m->getModelType(c2c) << "c"
                << e_lens.m->getModelTag(c2c) << std::endl;
  }

  // Create the new 3-cell triangles.
  if (e_lens.e_3c_t_flag) {
    //int c3_f = ng.cells.at(ng.ngons.at(ng.curr).at(0)).first.at(0);
    //apf::ModelEntity* mdl = e_lens.m->findModelEntity(3, c3_f);
    apf::ModelEntity* mdl = e_lens.m->findModelEntity(3, e_lens.c3_f);

    for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      down[0] = v3_t;
      down[1] = e_lens.v_2c;
      down[2] = e_lens.vert_ctr_new.at(i);
      tri_f_ctr.at(i) = buildElement(e_lens.m, mdl, 
                                      apf::Mesh::TRIANGLE, down, 0);

      std::cout << "tri_f_ctr(" << i << ") " 
                << vd_meas_ent(e_lens.m, tri_f_ctr.at(i)) 
                << std::endl;
      std::cout << "tri_f_ctr(" << i << ") " << tri_f_ctr.at(i) << " " 
                << e_lens.m->getModelType(mdl) << "c"
                << e_lens.m->getModelTag(mdl) << std::endl;

    }

    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      down[0] = v3_t;
      down[1] = e_lens.v_2c;
      down[2] = e_lens.vert_sp_new.at(i);
      tri_f_sp.at(i) = buildElement(e_lens.m, mdl, 
                                        apf::Mesh::TRIANGLE, down, 0);

      std::cout << "tri_f_sp(" << i << ") " 
                << vd_meas_ent(e_lens.m, tri_f_sp.at(i)) 
                << std::endl;
      std::cout << "tri_f_sp(" << i << ") " << tri_f_sp.at(i) << " " 
                << e_lens.m->getModelType(mdl) << "c"
                << e_lens.m->getModelTag(mdl) << std::endl;

    }
  }

  if (e_lens.e_3c_b_flag) {
    apf::ModelEntity* mdl = e_lens.m->findModelEntity(3, e_lens.c3_e);

    for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      down[0] = v3_b;
      down[1] = e_lens.v_2c;
      down[2] = e_lens.vert_ctr_new.at(i);
      tri_e_ctr.at(i) = buildElement(e_lens.m, mdl, 
                                      apf::Mesh::TRIANGLE, down, 0);

      std::cout << "tri_e_ctr(" << i << ") " 
                << vd_meas_ent(e_lens.m, tri_e_ctr.at(i)) 
                << std::endl;

      std::cout << "tri_e_ctr(" << i << ") " << tri_e_ctr.at(i) << " " 
                << e_lens.m->getModelType(mdl) << "c"
                << e_lens.m->getModelTag(mdl) << std::endl;
    }

    for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      down[0] = v3_b;
      down[1] = e_lens.v_2c;
      down[2] = e_lens.vert_sp_new.at(i);
      tri_e_sp.at(i) = buildElement(e_lens.m, mdl, 
                                      apf::Mesh::TRIANGLE, down, 0);

      std::cout << "tri_e_sp(" << i << ") " 
                << vd_meas_ent(e_lens.m, tri_e_sp.at(i)) 
                << std::endl;

      std::cout << "tri_e_sp(" << i << ") " << tri_e_sp.at(i) << " " 
                << e_lens.m->getModelType(mdl) << "c"
                << e_lens.m->getModelTag(mdl) << std::endl;
    }
  }

  // Create the new triangles.
  for(int i = 0; i < e_lens.discs.size(); i++) {
    int v1 = e_lens.slices.at(2*i)-1;
    int v2 = e_lens.slices.at(2*i+1)-1;

    down[0] = e_lens.vert_ctr_new.at(v1);
    down[1] = e_lens.vert_sp_new.at(i);
    down[2] = e_lens.v_2c;

    tri_2.at(2*i) = buildElement(e_lens.m, c2c, 
                                  apf::Mesh::TRIANGLE, down, 0);

    std::cout << "tri_2 [" << 2*i << "] " << tri_2.at(2*i) << std::endl;

    std::cout   << vd_meas_ent(e_lens.m, tri_2.at(2*i)) 
                << std::endl;

    std::cout << "tri_2 [" << 2*i << "] " << tri_2.at(2*i) << " "
              << e_lens.m->getModelType(c2c) << "c"
              << e_lens.m->getModelTag(c2c) << std::endl;

    down[0] = e_lens.vert_ctr_new.at(v2);
    down[1] = e_lens.vert_sp_new.at(i);
    down[2] = e_lens.v_2c;

    tri_2.at(2*i+1) = buildElement(e_lens.m, c2c, 
                                   apf::Mesh::TRIANGLE, down, 0);

    std::cout << "tri_2 [" << 2*i+1 << "] " << tri_2.at(2*i+1) << std::endl;

    std::cout   << vd_meas_ent(e_lens.m, tri_2.at(2*i+1)) 
                << std::endl;

    std::cout << "tri_2 [" << 2*i+1 << "] " << tri_2.at(2*i+1)
              << e_lens.m->getModelType(c2c) << "c"
              << e_lens.m->getModelTag(c2c) << std::endl;

  }

  vd_disc disc_new;
  disc_new.m = e_lens.m;
  disc_new.v_ctr = e_lens.v_2c;
  cp_tri(&tri_2, &disc_new);

  // The disc cannot be broken.
  assert(!disc_new.broken);

  // Create the new elements.
  // Top
  // Create the new 3-cell triangles.
  std::vector<apf::MeshEntity*> tet_void(0);

  // Find a disc triangle such that the vertices at the ends of bounding edges
  // on the discs and the top vertex form a triangle that has a tet adjacent to
  // it (it will be a lens tet). Get the vertex ordering of that tet and switch
  // the remaining vertex with the e_lens.v_2c. Going along the disc, follow
  // the same ordering to generate the tets.
  bool flip = false;

  apf::MeshEntity* e_curr = NULL;
  apf::MeshEntity* e_next = NULL;
  apf::MeshEntity* tri_curr = NULL;
  apf::MeshEntity* e_void = NULL;
  int i_st = -1;

  int i_3c = 0;
  int i1 = 1;
  int i2 = 2;
  int i_2c = 3;

  apf::ModelEntity* mdl_t = NULL;
  apf::ModelEntity* mdl_b = NULL;
  apf::MeshEntity* v_1 = NULL;
  apf::MeshEntity* v_2 = NULL;

  if (e_lens.e_3c_t_flag)
    mdl_t = e_lens.m->findModelEntity(3, e_lens.c3_f);
  if (e_lens.e_3c_b_flag)
    mdl_b = e_lens.m->findModelEntity(3, e_lens.c3_e);

  if (e_lens.e_3c_t_flag) {
    for(int i = 0; i < disc_new.tri.size(); i++) {
      tri_curr = disc_new.tri.at(i);

      e_curr = disc_new.edge.at(disc_new.t_e_map1[tri_curr] - 1);
      e_next = disc_new.edge.at(disc_new.t_e_map2[tri_curr] - 1);

      v_1 = getEdgeVertOppositeVert(e_lens.m, e_curr, e_lens.v_2c);

      v_2 = getEdgeVertOppositeVert(e_lens.m, e_next, e_lens.v_2c);

      down[0] = v3_t;
      down[i1] = v_1;
      down[i2] = v_2;

      apf::MeshEntity* tri_leanto = findElement(e_lens.m, apf::Mesh::TRIANGLE,
                                                                     down);
      if(tri_leanto) {

        e_lens.m->getUp(tri_leanto, up);
        assert(up.n < 2);
        if(up.n == 1) {
          // Flip the indices of 3c and 2c vertices:
          e_lens.m->getDownward(up.e[0], 0, down);
          i_2c = findIn(down, 4, v3_t);
          i1 = findIn(down, 4, v_1);
          i2 = findIn(down, 4, v_2);
          i_3c = 6 - i_2c - i1 - i2;

          assert(i1 > -1 and i2 > -1 and i_3c > -1 and i_2c > -1);

          i_st = i;
          i = disc_new.tri.size();
        }
      }
    }
  }
  else {
    assert(e_lens.e_3c_b_flag);
    for(int i = 0; i < disc_new.tri.size(); i++) {
      tri_curr = disc_new.tri.at(i);

      e_curr = disc_new.edge.at(disc_new.t_e_map1[tri_curr] - 1);
      e_next = disc_new.edge.at(disc_new.t_e_map2[tri_curr] - 1);

      v_1 = getEdgeVertOppositeVert(e_lens.m, e_curr, e_lens.v_2c);
      v_2 = getEdgeVertOppositeVert(e_lens.m, e_next, e_lens.v_2c);

      down[0] = v3_b;
      down[i1] = v_1;
      down[i2] = v_2;

      apf::MeshEntity* tri_leanto = findElement(e_lens.m, apf::Mesh::TRIANGLE,
                                                                     down);
      if(tri_leanto) {

        e_lens.m->getUp(tri_leanto, up);
        assert(up.n < 2);
        if(up.n == 1) {
          e_lens.m->getDownward(up.e[0], 0, down);
          // Flip the indices of 3c and 2c vertices:
          i_2c = findIn(down, 4, v3_b);
          i1 = findIn(down, 4, v_1);
          i2 = findIn(down, 4, v_2);
          i_3c = 6 - i_2c - i1 - i2;

          assert(i1 > -1 and i2 > -1 and i_3c > -1 and i_2c > -1);

          i_st = i;
          i = disc_new.tri.size();
        }
      }
    }
  }
  assert(i_st > -1);

  apf::Downward d_tri;

  tet_void.reserve(4*e_lens.discs.size());
  apf::MeshEntity* e_init = e_curr;

  if (e_lens.e_3c_t_flag) {
    for(int i = 0; i < disc_new.tri.size(); i++) {
      down[i_3c] = v3_t;
      down[i1] = v_1;
      down[i2] = v_2;
      down[i_2c] = e_lens.v_2c;

      e_void = buildElement(e_lens.m, mdl_t, apf::Mesh::TET, down, 0);
      tet_void.push_back(e_void);

      // Sanity check for downward triangles:
      // e_lens.m->getDownward(e_void, 2, d_tri);

      if (e_lens.e_3c_b_flag) {
        down[i_2c] = v3_b;
        down[i1] = v_1;
        down[i2] = v_2;
        down[i_3c] = e_lens.v_2c;
        e_void = buildElement(e_lens.m, mdl_b, apf::Mesh::TET, down, 0);
        tet_void.push_back(e_void);
      }

      // Collect the next edges and vertices:
      e_curr = e_next;

      if(tri_curr == disc_new.tri.at(disc_new.e_t_map1[e_curr] - 1)) {
        tri_curr = disc_new.tri.at(disc_new.e_t_map2[e_curr] - 1);
      }
      else {
        assert(tri_curr == disc_new.tri.at(disc_new.e_t_map2[e_curr] - 1));
        tri_curr = disc_new.tri.at(disc_new.e_t_map1[e_curr] - 1);
      }
      if(disc_new.edge.at(disc_new.t_e_map1[tri_curr] - 1) == e_curr) {
        e_next = disc_new.edge.at(disc_new.t_e_map2[tri_curr] - 1);
      }
      else {
        assert(disc_new.edge.at(disc_new.t_e_map2[tri_curr] - 1) == e_curr);
        e_next = disc_new.edge.at(disc_new.t_e_map1[tri_curr] - 1);
      }
      v_1 = getEdgeVertOppositeVert(e_lens.m, e_curr, e_lens.v_2c);
      v_2 = getEdgeVertOppositeVert(e_lens.m, e_next, e_lens.v_2c);

    }
  }
  else {
    for(int i = 0; i < disc_new.tri.size(); i++) {
      down[i_3c] = v3_b;
      down[i1] = v_1;
      down[i2] = v_2;
      down[i_2c] = e_lens.v_2c;

      e_void = buildElement(e_lens.m, mdl_b, apf::Mesh::TET, down, 0);
      tet_void.push_back(e_void);
      // Collect the next edges and vertices:
      e_curr = e_next;

      if(tri_curr == disc_new.tri.at(disc_new.e_t_map1[e_curr] - 1)) {
        tri_curr = disc_new.tri.at(disc_new.e_t_map2[e_curr] - 1);
      }
      else {
        assert(tri_curr == disc_new.tri.at(disc_new.e_t_map2[e_curr] - 1));
        tri_curr = disc_new.tri.at(disc_new.e_t_map1[e_curr] - 1);
      }
      if(disc_new.edge.at(disc_new.t_e_map1[tri_curr] - 1) == e_curr) {
        e_next = disc_new.edge.at(disc_new.t_e_map2[tri_curr] - 1);
      }
      else {
        assert(disc_new.edge.at(disc_new.t_e_map2[tri_curr] - 1) == e_curr);
        e_next = disc_new.edge.at(disc_new.t_e_map1[tri_curr] - 1);
      }
      v_1 = getEdgeVertOppositeVert(e_lens.m, e_curr, e_lens.v_2c);
      v_2 = getEdgeVertOppositeVert(e_lens.m, e_next, e_lens.v_2c);
    }
  }
  assert(e_curr == e_init);
  f_calc.vd_att_fields(e_lens.m, e_lens.v_2c);
}

int vd_edisc::get_circ_type(int circ_in) {
  return vd_cd->get_circ_type(circ_in);
/*
  if(spur)
    vd_cd->get_circ_type(circ_in);
  else
    c_ins->get_circ_type_gmi(cell_id, circ_in);
*/
}

circuit* vd_edisc::get_circ(int circ_in) {
  return vd_cd->get_circ(circ_in);
/*
  if(spur)
    return vd_cd->get_circ(circ_in);
  else
    return c_ins->get_circ_gmi(cell_id, circ_in);
*/
}

void vd_edisc::get_circ_topo_gmi(std::vector<int>* c2_list, std::vector<int>* c3_list, circuit* circ) {
  vd_cd->get_circ_topo_gmi(c2_list, c3_list, circ);

}

void vd_edisc::get_path_ngon(int c3_cp, ngon_gmi* ng) {
  vd_cd->get_path_ngon(c3_cp, ng);
/*
  if(spur)
    vd_cd->get_path_ngon_gmi(c3_cp, ng);
  else
    c_ins->get_path_ngon_gmi(cell_id, c3_cp, ng);
*/
}

void vd_edisc::conv_path_ngon_gmi(ngon_gmi* ng) {
  vd_cd->conv_path_ngon_gmi(ng);
}

// Collect the affected 1cells on slices.
void vd_edisc::collect_cell() {
  std::cout << "Collecting cells 1cells" << std::endl;
  if(trial_type == 1)
    e_lens.c1.resize(2);
  if(trial_type == 2)
    e_lens.c1.resize(e_lens.discs.size());

  ent_conn* c1_adj = new ent_conn();
  e_lens.c_base_curr->get_conn_dim_gmi(1, 0, cell_id, c1_adj);
  std::map<int, bool> c1_map{};
  std::map<int, int> c1_map_id{};

  for(int i = 0; i < e_lens.c1.size(); i++) {
    e_lens.c1.at(i).reserve(c1_adj->conn.size());
  }
  e_lens.c1_mult.reserve(c1_adj->conn.size());

  // Check the slice of each edge, add to the slice 1cell list if not in the 
  // list.
  for(int i = 0; i < e_lens.es_edge.size(); i++) {
    apf::MeshEntity* ent_c = e_lens.es_edge.at(i);
    //std::cout << "Check slice id of " << ent_c << std::endl;
    int slice_id = e_lens.slice_map[ent_c];
    //std::cout << slice_id << std::endl;

    if(slice_id > 0) {
      apf::ModelEntity* mdl = e_lens.m->toModel(ent_c);
      int dim = e_lens.m->getModelType(mdl);
      int tag = e_lens.m->getModelTag(mdl);
      if(dim == 1) {
        if(!c1_map[tag]) {
          //std::cout << slice_id << ": " << tag << " " << std::endl;
          e_lens.c1.at(slice_id-1).push_back(tag);
          c1_map[tag] = true;
          c1_map_id[tag] = slice_id;
        }
        else {
          if(c1_map_id[tag] != slice_id) {
            e_lens.c1.at(slice_id-1).push_back(tag);
            e_lens.c1_mult.push_back(tag);
          }
        }
      }
    }
  }
  for(int i = 0; i < e_lens.c1.size(); i++) {
    std::sort(e_lens.c1.at(i).begin(), e_lens.c1.at(i).end());
  }
  delete c1_adj;
  //std::cout << std::endl;
}

// Collect the affected 1cells on slices, without relying on mesh information.
void vd_edisc::collect_cell_wg(std::vector< std::pair< std::pair<int,int>, 
                         std::vector<std::vector<int > > > >* slice_cells) {
  std::cout << "Collecting cells 1cells" << std::endl;
  if(trial_type == 1)
    e_lens.c1.resize(2);
  if(trial_type == 2)
    e_lens.c1.resize(e_lens.discs.size());

  std::map<int, bool> c1_map{};
  for(int slice = 0; slice < slice_cells->size(); slice++) {
    e_lens.c1.at(slice).resize(slice_cells->at(slice).second.at(0).size());
    for(int i = 0; i < slice_cells->at(slice).second.at(0).size(); i++) {
      int c1_disj_curr = slice_cells->at(slice).second.at(0).at(i);
      int c1_curr = e_lens.m->getModelTag(vd_cd->get_mdl(1, c1_disj_curr));
      e_lens.c1.at(slice).at(i) = c1_curr;

      if(c1_map[c1_curr])
        e_lens.c1_mult.push_back(c1_curr);
      c1_map[c1_curr] = true;

    }
  }
}

// Update the cell_base structure.
void vd_edisc::upd_cell() {
  std::cout << "Updating the cell base." << std::endl;
  ent_conn* e_con = new ent_conn();
  ent_conn* e_con2 = new ent_conn();
  e_con2->resize(e_lens.discs.size());

  e_con->clear();
  // Empty the upward adjacency lists:
  if(trial_type == 1) {
    for(int i = 0; i < e_lens.new_cell0_id.size(); i++)
      e_lens.c_base_curr->set_conn_up_gmi(0, e_lens.new_cell0_id.at(i), e_con);
    e_lens.c_base_curr->set_conn_up_gmi(1, e_lens.new_cell1_id.at(0), e_con);
    e_lens.c_base_curr->set_conn_gmi(1, e_lens.new_cell1_id.at(0), e_con);
  }
  else {
    assert(trial_type == 2);
    for(int i = 0; i < e_lens.new_cell0_id.size(); i++)
      e_lens.c_base_curr->set_conn_up_gmi(0, e_lens.new_cell0_id.at(i), e_con);
    for(int i = 0; i < e_lens.new_cell1_id.size(); i++) {
      e_lens.c_base_curr->set_conn_up_gmi(1, e_lens.new_cell1_id.at(i), e_con);
      e_lens.c_base_curr->set_conn_gmi(1, e_lens.new_cell1_id.at(i), e_con);
    }
    e_lens.c_base_curr->set_conn_up_gmi(2, e_lens.new_cell2_id, e_con);
    e_lens.c_base_curr->set_conn_gmi(2, e_lens.new_cell2_id, e_con);
  }

  // These are the c1 bounded by a single c0 and that belong to multiple slices.
  // Repeat the c0 so that the replacement can be assigned twice.
  for(int i = 0; i < e_lens.c1_mult.size(); i++) {
    e_con->conn.clear();
    e_con->conn.resize(2);
    e_con->conn.at(0) = cell_id;
    e_con->conn.at(1) = cell_id;
    e_lens.c_base_curr->set_conn_gmi(1, e_lens.c1_mult.at(i), e_con);
  }

  for(int i = 0; i < e_lens.c1.size(); i++) {
    for(int j = 0; j < e_lens.c1.at(i).size(); j++) {
      std::cout << "c[" << i << "][" << j << "]=" 
                << e_lens.c1.at(i).at(j) << std::endl;
      // Replace the trial 0cell id with the 0cell corresponding to the
      // slice. A 1cell connected to the 0cell is bounded by either one 0cell 
      // or two 0cells, unless it is repeated in multiple slices.
      e_lens.c_base_curr->get_conn_gmi(1, e_lens.c1.at(i).at(j), e_con);
      if(e_con->conn.size() == 1) {
        assert(e_con->conn.at(0) == cell_id);
        e_con->conn.at(0) = e_lens.new_cell0_id.at(i);
      }
      else {
        if(e_con->conn.at(0) == cell_id)
          e_con->conn.at(0) = e_lens.new_cell0_id.at(i);
        else {
          assert(e_con->conn.at(1) == cell_id);
          e_con->conn.at(1) = e_lens.new_cell0_id.at(i);
        }
      }
      e_lens.c_base_curr->set_conn_gmi(1, e_lens.c1.at(i).at(j), e_con);
    }
  }

  for(int i = 0; i < e_lens.discs.size(); i++) {
    e_con->clear();
    e_con->resize(2);
    int v1 = e_lens.slices.at(2*i)-1;
    int v2 = e_lens.slices.at(2*i+1)-1;
    e_con->conn.at(0) = e_lens.new_cell0_id.at(v1);
    e_con->conn.at(1) = e_lens.new_cell0_id.at(v2);

    e_lens.c_base_curr->set_conn_gmi(1, e_lens.new_cell1_id.at(i), e_con);
    std::cout << "1c" << e_lens.new_cell1_id.at(i) 
              << " 0c" << e_con->conn.at(0)
              << " 0c" << e_con->conn.at(1)
              << std::endl;

    e_con2->conn.at(i) = e_lens.new_cell1_id.at(i);
  }

  // TODO Convert these to use the cell_adder object.
  // TODO These should be registered to the e_lens.
  // Given the slices, collect the 1cells to associate with the slice.
  // This is less complex in the case of the 2cell insertion, where the slice
  // 1cells cannot be simply obtained from the paths.
  // If a 1cell is inserted, update the 2cells of the circuit. Insert the 
  // 1cell.
  // If a 2cell is inserted, update the 2cells of the paths, update the 3cell
  // couple. Insert the 2cell.

  if(trial_type == 1) {
    // Update the 0/1-cell and the 1/2-cell adjacencies.
    // Get the 1cells in both subgraphs adjacent to the old 0cell. Reattach 
    // the groups of 1cells to the corresponding vert_ctr_new.
    // Get the 2cells on the circuit. Add the new 1cell to their adjacency 
    // lists.
    circuit* circ_tup;

    circ_tup = get_circ(trial_curr);
    ent_conn* c2_circ = new ent_conn();
    ent_conn* c3_circ = new ent_conn();

    get_circ_topo_gmi(&c2_circ->conn, &c3_circ->conn, circ_tup); 

    // In case a 2stratum is has two disjoint components around the 0stratum 
    // vertex, it can be repeated twice on the circuit.
    std::map<int, bool> c2_map{};
    // Flag the new 1cell and 0cells as exterior if one circuit 2cell is on the
    // exterior.
    bool exterior = false;

    for(int i = 0; i < c2_circ->conn.size(); i++) {
      e_lens.c_base_curr->get_conn_gmi(2, c2_circ->conn.at(i), e_con);

      if(!c2_map[c2_circ->conn.at(i)]) {
        std::cout << "2c" << c2_circ->conn.at(i)
                << " added 1c" << e_lens.new_cell1_id.at(0)
                << std::endl;

        e_con->add_ent(e_lens.new_cell1_id.at(0));
        e_lens.c_base_curr->set_conn_gmi(2, c2_circ->conn.at(i), e_con);
        c2_map[c2_circ->conn.at(i)] = true;
        //if(e_lens.c_base_curr->chk_cell_ext_gmi(2, c2_circ->conn.at(i)))
        if(e_lens.c_base_curr->get_cell_ext_gmi(2, c2_circ->conn.at(i)))
          exterior = true;
      }
    }
    if(exterior) {
      e_lens.c_base_curr->set_ext_gmi(1, e_lens.new_cell1_id.at(0), true);
      for(int i = 0; i < e_lens.new_cell0_id.size(); i++)
        e_lens.c_base_curr->set_ext_gmi(0, e_lens.new_cell0_id.at(i), true);
    }
    // Check if any of the 1strata shares the same 2strata as the new one.
    std::sort(c3_circ->conn.begin(), c3_circ->conn.end());
    std::sort(c2_circ->conn.begin(), c2_circ->conn.end());
    bool spur_1c = false;
    for(int i = 0; i < e_lens.c1.size(); i++) {
      for(int j = 0; j < e_lens.c1.at(i).size(); j++) {
        e_lens.c_base_curr->get_conn_dim_gmi(2, 1, e_lens.c1.at(i).at(j), 
                                                                    e_con);
        std::sort(e_con->conn.begin(), e_con->conn.end());
        if(c2_circ->conn == e_con->conn) {
          j = e_lens.c1.at(i).size();
          spur_1c = true;
        }
      }
      if(spur_1c)
        i = e_lens.c1.size();
    }
    if(spur_1c) {
      std::cout << "ins_flag: Spurious insertion at the exterior." << std::endl;
      ins_flag = false;
    }
    if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL)
      determine_shell_gmi();
    delete c2_circ;
    delete c3_circ;
  }

  if(trial_type == 2) {
    // Update the 0/1-cell, 1/2-cell and 2/3-cell adjacencies.
    // Get the 1cells in both slices adjacent to the old 0cell. Reattach 
    // the groups of 1cells to the corresponding vert_ctr_new.
    // Get the 2cells on the paths. Add the corresponding new 1cell to their 
    // adjacency lists. Add the new 2cell to the adjacency list of both of  
    // the 3cell couples.
    conv_path_ngon_gmi(&ng);
    // Obtain the list of available shells. These will be used to determine the
    // shell of the new 2cell by comparing dissipation rates.
    if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL and ext_0cell) {
      ext_shell* e_sh = f_calc.get_e_sh();
      bool ext_3c = (ng.cp.first == -1 or ng.cp.second == -1);
      if(ext_3c) {
        std::vector<std::vector<int> > sh2_cell_list(0, std::vector<int> (0));
        std::map<int, int> sh2id_map{};
        ent_conn c2_local;
        ent_conn sh2_local;
        e_lens.c_base_curr->get_conn_dim_gmi(2, 0, cell_id, &c2_local);
        if(sh_old.dim == 2) {
          shell2_list.resize(1);
          shell2_list.at(0) = sh_old.id;
          // Remove from the list of 2cells for each shell the 2cells bounding  
          // the current 3cell.
          if(ng.cp.first == -1)
            e_lens.c_base_curr->get_conn_gmi(3, ng.cp.second, &c2_local);
          else
            e_lens.c_base_curr->get_conn_gmi(3, ng.cp.first, &c2_local);
          sh2_local.clear();
          sh2_local.reserve(c2_local.conn.size());
          for(int j = 0; j < c2_local.conn.size(); j++) {
            int ng_id = c2_local.conn.at(j);
            if(e_lens.c_base_curr->get_cell_ext_gmi(2, ng_id) ) {
              assert(e_sh->chk_shell(2, ng_id - 1) );
              sh2_local.conn.push_back(e_sh->get_shell(2, ng_id - 1).id);
            }
          }
          for(int i = 0; i < shell2_list.size(); i++) {
            rem_ents(&shell2_list, &sh2_local.conn);
          }
        }
        else {
          shell2_list.clear();
          e_sh->sh_base.get_conn_dim_gmi(2, sh_old.dim, sh_old.id, &sh2_local);

          sh2_cell_list.resize(sh2_local.conn.size());
          // For each shell, collect 2cells attached to it and the 0cell.
          for(int i = 0; i < sh2_local.conn.size(); i++) {
            sh2id_map[sh2_local.conn.at(i)] = i;
            sh2_cell_list.at(i).reserve(c2_local.conn.size());
          }
          // Obtain the list of 2cells for each shell.
          for(int i = 0; i < ng.ngons.at(ng.curr).size(); i++) {
            bool found = false;
            for(int j = 0; j < ng.cells.at(ng.ngons.at(ng.curr).at(i))
                                                          .second.size(); j++) {
              int ng_id = ng.cells.at(ng.ngons.at(ng.curr).at(i)).second.at(j);
              if(e_lens.c_base_curr->get_cell_ext_gmi(2, ng_id) ) {
                assert(!found);
                assert(e_sh->chk_shell(2, ng_id - 1) );
                shell sh_curr = e_sh->get_shell(2, ng_id - 1);
                sh2_cell_list.at(sh2id_map[sh_curr.id]).push_back(ng_id);
                found = true;
              }
            }
          }
          // Remove from the list of 2cells for each shell the 2cells bounding  
          // the current 3cell.
          if(ng.cp.first == -1)
            e_lens.c_base_curr->get_conn_gmi(3, ng.cp.second, &c2_local);
          else
            e_lens.c_base_curr->get_conn_gmi(3, ng.cp.first, &c2_local);
          for(int i = 0; i < sh2_cell_list.size(); i++) {
            rem_ents(&sh2_cell_list.at(i), &c2_local.conn);
          }
          for(int i = 0; i < sh2_cell_list.size(); i++) {
            if(sh2_cell_list.at(i).size() > 0)
              shell2_list.push_back(sh2_local.conn.at(i));
          }
        }
        std::sort(shell2_list.begin(), shell2_list.end());
        std::vector<int>::iterator it;
        it = std::unique (shell2_list.begin(), shell2_list.end());
        shell2_list.resize(std::distance(shell2_list.begin(),it));

        std::cout << "Available 2shells are ";
        for(int i = 0; i < shell2_list.size(); i++) {
          std::cout << shell2_list.at(i) << " ";
        }
        std::cout << std::endl;
      }
    }

    for(int i = 0; i < ng.ngons.at(ng.curr).size(); i++) {
      // In case there are disjoint sets of entities belonging to the same 
      // 2stratum, add the new 1stratum once.
      std::map<int, bool> c2_map {};
      int p_id = ng.ngons.at(ng.curr).at(i);
      for(int j = 0; j < ng.cells.at(p_id).second.size(); j++) {
        int ng_id = ng.cells.at(p_id).second.at(j);
        e_lens.c_base_curr->get_conn_gmi(2, ng_id, e_con);

        if(!c2_map[ng_id]) {

          e_con->add_ent(e_lens.new_cell1_id.at(i));

          std::cout << "2c" << ng_id
                << " added 1c" << e_lens.new_cell1_id.at(i)
                << std::endl;

          e_lens.c_base_curr->set_conn_gmi(2, ng_id, e_con);
          c2_map[ng_id] = true;
        }
      }
      // Check if the exterior is on the path.
      std::vector<int>* c3_vec = &ng.cells.at(p_id).first;
      std::vector<int>::iterator c3_it;
      c3_it = std::find(c3_vec->begin(), c3_vec->end(), -1);

      if(c3_it != c3_vec->end()) {
        e_lens.c_base_curr->set_ext_gmi(1, e_lens.new_cell1_id.at(i), true);
      }
    }

    e_lens.c_base_curr->set_conn_gmi(2, e_lens.new_cell2_id, e_con2);

    if(e_lens.e_3c_t_flag) {
      e_lens.c_base_curr->get_conn_gmi(3, ng.cp.first, e_con);
      e_con->add_ent(e_lens.new_cell2_id);
      e_lens.c_base_curr->set_conn_gmi(3, ng.cp.first, e_con);

      std::cout << "3c" << ng.cp.first
              << " added 2c" << e_lens.new_cell2_id
              << std::endl;
    }

    if(e_lens.e_3c_b_flag) {
      e_lens.c_base_curr->get_conn_gmi(3, ng.cp.second, e_con);
      e_con->add_ent(e_lens.new_cell2_id);
      e_lens.c_base_curr->set_conn_gmi(3, ng.cp.second, e_con);

      std::cout << "3c" << ng.cp.second
              << " added 2c" << e_lens.new_cell2_id
              << std::endl;
    }

    if(!e_lens.e_3c_b_flag or !e_lens.e_3c_t_flag) {
      for(int i = 0; i < e_lens.new_cell1_id.size(); i++)
        e_lens.c_base_curr->set_ext_gmi(1, e_lens.new_cell1_id.at(i), true);
      for(int i = 0; i < e_lens.new_cell0_id.size(); i++)
        e_lens.c_base_curr->set_ext_gmi(0, e_lens.new_cell0_id.at(i), true);
      e_lens.c_base_curr->set_ext_gmi(2, e_lens.new_cell2_id, true);
    }

    if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL)
      determine_shell_gmi();
  }
  cell_flag = true;

  delete e_con;
  delete e_con2;
}

// Update the cell_base structure.
void vd_edisc::upd_cell_wg() {
  std::cout << "Updating the cell base." << std::endl;
  ent_conn* e_con = new ent_conn();
  ent_conn* e_con2 = new ent_conn();
  e_con2->resize(e_lens.discs.size());

  e_con->clear();
  // Empty the upward adjacency lists:
  if(trial_type == 1) {
    for(int i = 0; i < e_lens.new_cell0_id.size(); i++)
      e_lens.c_base_curr->set_conn_up_gmi(0, e_lens.new_cell0_id.at(i), e_con);
    e_lens.c_base_curr->set_conn_up_gmi(1, e_lens.new_cell1_id.at(0), e_con);
    e_lens.c_base_curr->set_conn_gmi(1, e_lens.new_cell1_id.at(0), e_con);
  }
  else {
    assert(trial_type == 2);
    for(int i = 0; i < e_lens.new_cell0_id.size(); i++)
      e_lens.c_base_curr->set_conn_up_gmi(0, e_lens.new_cell0_id.at(i), e_con);
    for(int i = 0; i < e_lens.new_cell1_id.size(); i++) {
      e_lens.c_base_curr->set_conn_up_gmi(1, e_lens.new_cell1_id.at(i), e_con);
      e_lens.c_base_curr->set_conn_gmi(1, e_lens.new_cell1_id.at(i), e_con);
    }
    e_lens.c_base_curr->set_conn_up_gmi(2, e_lens.new_cell2_id, e_con);
    e_lens.c_base_curr->set_conn_gmi(2, e_lens.new_cell2_id, e_con);
  }

  // These are the c1 bounded by a single c0 and that belong to multiple slices.
  // Repeat the c0 so that the replacement can be assigned.
  for(int i = 0; i < e_lens.c1_mult.size(); i++) {
    e_lens.c_base_curr->get_conn_gmi(1, e_lens.c1_mult.at(i), e_con);
    assert(e_con->conn.size() == 1 and e_con->conn.at(0) == cell_id);
    e_con->conn.resize(2);
    e_con->conn.at(1) = cell_id;
    e_lens.c_base_curr->set_conn_gmi(1, e_lens.c1_mult.at(i), e_con);
  }

  for(int i = 0; i < e_lens.c1.size(); i++) {
    for(int j = 0; j < e_lens.c1.at(i).size(); j++) {
      std::cout << "c[" << i << "][" << j << "]=" 
                << e_lens.c1.at(i).at(j) << std::endl;
      // Replace the trial 0cell id with the 0cell corresponding to the
      // slice. A 1cell connected to the 0cell is bounded by either one 0cell 
      // or two 0cells, unless it is repeated in multiple slices.
      e_lens.c_base_curr->get_conn_gmi(1, e_lens.c1.at(i).at(j), e_con);
      if(e_con->conn.size() == 1) {
        assert(e_con->conn.at(0) == cell_id);
        e_con->conn.at(0) = e_lens.new_cell0_id.at(i);
      }
      else {
        if(e_con->conn.at(0) == cell_id)
          e_con->conn.at(0) = e_lens.new_cell0_id.at(i);
        else {
          assert(e_con->conn.at(1) == cell_id);
          e_con->conn.at(1) = e_lens.new_cell0_id.at(i);
        }
      }
      e_lens.c_base_curr->set_conn_gmi(1, e_lens.c1.at(i).at(j), e_con);
    }
  }

  for(int i = 0; i < e_lens.discs.size(); i++) {
    e_con->clear();
    e_con->resize(2);
    int v1 = e_lens.slices.at(2*i)-1;
    int v2 = e_lens.slices.at(2*i+1)-1;
    e_con->conn.at(0) = e_lens.new_cell0_id.at(v1);
    e_con->conn.at(1) = e_lens.new_cell0_id.at(v2);

    e_lens.c_base_curr->set_conn_gmi(1, e_lens.new_cell1_id.at(i), e_con);
    std::cout << "1c" << e_lens.new_cell1_id.at(i) 
              << " 0c" << e_con->conn.at(0)
              << " 0c" << e_con->conn.at(1)
              << std::endl;

    e_con2->conn.at(i) = e_lens.new_cell1_id.at(i);
  }

  // TODO Convert these to use the cell_adder object.
  // TODO These should be registered to the e_lens.
  // Given the slices, collect the 1cells to associate with the slice.
  // This is less complex in the case of the 2cell insertion, where the slice
  // 1cells cannot be simply obtained from the paths.
  // If a 1cell is inserted, update the 2cells of the circuit. Insert the 
  // 1cell.
  // If a 2cell is inserted, update the 2cells of the paths, update the 3cell
  // couple. Insert the 2cell.

  if(trial_type == 1) {
    // Update the 0/1-cell and the 1/2-cell adjacencies.
    // Get the 1cells in both subgraphs adjacent to the old 0cell. Reattach 
    // the groups of 1cells to the corresponding vert_ctr_new.
    // Get the 2cells on the circuit. Add the new 1cell to their adjacency 
    // lists.
    circuit* circ_tup;

    circ_tup = get_circ(trial_curr);
    ent_conn* c2_circ = new ent_conn();
    ent_conn* c3_circ = new ent_conn();

    get_circ_topo_gmi(&c2_circ->conn, &c3_circ->conn, circ_tup); 

    // In case a 2stratum is has two disjoint components around the 0stratum 
    // vertex, it can be repeated twice on the circuit.
    std::map<int, bool> c2_map{};
    // Flag the new 1cell and 0cells as exterior if one circuit 2cell is on the
    // exterior.
    bool exterior = false;

    for(int i = 0; i < c2_circ->conn.size(); i++) {
      e_lens.c_base_curr->get_conn_gmi(2, c2_circ->conn.at(i), e_con);

      if(!c2_map[c2_circ->conn.at(i)]) {
        std::cout << "2c" << c2_circ->conn.at(i)
                << " added 1c" << e_lens.new_cell1_id.at(0)
                << std::endl;

        e_con->add_ent(e_lens.new_cell1_id.at(0));
        e_lens.c_base_curr->set_conn_gmi(2, c2_circ->conn.at(i), e_con);
        c2_map[c2_circ->conn.at(i)] = true;
        //if(e_lens.c_base_curr->chk_cell_ext_gmi(2, c2_circ->conn.at(i)))
        if(e_lens.c_base_curr->get_cell_ext_gmi(2, c2_circ->conn.at(i)))
          exterior = true;
      }
    }
    if(exterior) {
      e_lens.c_base_curr->set_ext_gmi(1, e_lens.new_cell1_id.at(0), true);
      for(int i = 0; i < e_lens.new_cell0_id.size(); i++)
        e_lens.c_base_curr->set_ext_gmi(0, e_lens.new_cell0_id.at(i), true);
    }
    // Check if any of the 1strata shares the same 2strata as the new one.
    std::sort(c3_circ->conn.begin(), c3_circ->conn.end());
    std::sort(c2_circ->conn.begin(), c2_circ->conn.end());
    bool spur_1c = false;
    for(int i = 0; i < e_lens.c1.size(); i++) {
      for(int j = 0; j < e_lens.c1.at(i).size(); j++) {
        e_lens.c_base_curr->get_conn_dim_gmi(2, 1, e_lens.c1.at(i).at(j), 
                                                                    e_con);
        std::sort(e_con->conn.begin(), e_con->conn.end());
        if(c2_circ->conn == e_con->conn) {
          j = e_lens.c1.at(i).size();
          spur_1c = true;
        }
      }
      if(spur_1c)
        i = e_lens.c1.size();
    }
    if(spur_1c) {
      std::cout << "ins_flag: Spurious insertion at the exterior." << std::endl;
      ins_flag = false;
    }
    if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL)
      determine_shell_gmi_wg();
    delete c2_circ;
    delete c3_circ;
  }

  if(trial_type == 2) {
    // Update the 0/1-cell, 1/2-cell and 2/3-cell adjacencies.
    // Get the 1cells in both slices adjacent to the old 0cell. Reattach 
    // the groups of 1cells to the corresponding vert_ctr_new.
    // Get the 2cells on the paths. Add the corresponding new 1cell to their 
    // adjacency lists. Add the new 2cell to the adjacency list of both of  
    // the 3cell couples.
    conv_path_ngon_gmi(&ng);

    // Obtain the list of available shells. These will be used to determine the
    // shell of the new 2cell by comparing dissipation rates.
    if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL and ext_0cell) {
      ext_shell* e_sh = f_calc.get_e_sh();
      bool ext_3c = (ng.cp.first == -1 or ng.cp.second == -1);
      if(ext_3c) {
        std::vector<std::vector<int> > sh2_cell_list(0, std::vector<int> (0));
        std::map<int, int> sh2id_map{};
        ent_conn c2_local;
        ent_conn sh2_local;
        e_lens.c_base_curr->get_conn_dim_gmi(2, 0, cell_id, &c2_local);
        if(sh_old.dim == 2) {
          shell2_list.resize(1);
          shell2_list.at(0) = sh_old.id;
          // Remove from the list of 2cells for each shell the 2cells bounding  
          // the current 3cell.
          if(ng.cp.first == -1)
            e_lens.c_base_curr->get_conn_gmi(3, ng.cp.second, &c2_local);
          else
            e_lens.c_base_curr->get_conn_gmi(3, ng.cp.first, &c2_local);
          sh2_local.clear();
          sh2_local.reserve(c2_local.conn.size());
          for(int j = 0; j < c2_local.conn.size(); j++) {
            int ng_id = c2_local.conn.at(j);
            if(e_lens.c_base_curr->get_cell_ext_gmi(2, ng_id) ) {
              assert(e_sh->chk_shell(2, ng_id - 1) );
              sh2_local.conn.push_back(e_sh->get_shell(2, ng_id - 1).id);
            }
          }
          for(int i = 0; i < shell2_list.size(); i++) {
            rem_ents(&shell2_list, &sh2_local.conn);
          }
        }
        else {
          e_sh->sh_base.get_conn_dim_gmi(2, sh_old.dim, sh_old.id, &sh2_local);

          sh2_cell_list.resize(sh2_local.conn.size());
          // For each shell, collect 2cells attached to it and the 0cell.
          for(int i = 0; i < sh2_local.conn.size(); i++) {
            sh2id_map[sh2_local.conn.at(i)] = i;
            sh2_cell_list.at(i).reserve(c2_local.conn.size());
          }
          // Obtain the list of 2cells for each shell.
          for(int i = 0; i < ng.ngons.at(ng.curr).size(); i++) {
            bool found = false;
            for(int j = 0; j < ng.cells.at(ng.ngons.at(ng.curr).at(i))
                                                          .second.size(); j++) {
              int ng_id = ng.cells.at(ng.ngons.at(ng.curr).at(i)).second.at(j);
              if(e_lens.c_base_curr->get_cell_ext_gmi(2, ng_id) ) {
                assert(!found);
                assert(e_sh->chk_shell(2, ng_id - 1) );
                shell sh_curr = e_sh->get_shell(2, ng_id - 1);
                sh2_cell_list.at(sh2id_map[sh_curr.id]).push_back(ng_id);
                found = true;
              }
            }
          }
          // Remove from the list of 2cells for each shell the 2cells bounding  
          // the current 3cell.
          if(ng.cp.first == -1)
            e_lens.c_base_curr->get_conn_gmi(3, ng.cp.second, &c2_local);
          else
            e_lens.c_base_curr->get_conn_gmi(3, ng.cp.first, &c2_local);
          for(int i = 0; i < sh2_cell_list.size(); i++) {
            rem_ents(&sh2_cell_list.at(i), &c2_local.conn);
          }
          for(int i = 0; i < sh2_cell_list.size(); i++) {
            if(sh2_cell_list.at(i).size() > 0)
              shell2_list.push_back(sh2_local.conn.at(i));
          }
        }
        std::sort(shell2_list.begin(), shell2_list.end());
        std::vector<int>::iterator it;
        it = std::unique (shell2_list.begin(), shell2_list.end());
        shell2_list.resize(std::distance(shell2_list.begin(),it));

        std::cout << "Available 2shells are ";
        for(int i = 0; i < shell2_list.size(); i++) {
          std::cout << shell2_list.at(i) << " ";
        }
        std::cout << std::endl;
      }
    }

    for(int i = 0; i < ng.ngons.at(ng.curr).size(); i++) {
      // In case there are disjoint sets of entities belonging to the same 
      // 2stratum, add the new 1stratum once.
      std::map<int, bool> c2_map {};
      int p_id = ng.ngons.at(ng.curr).at(i);
      for(int j = 0; j < ng.cells.at(p_id).second.size(); j++) {
        int ng_id = ng.cells.at(p_id).second.at(j);
        e_lens.c_base_curr->get_conn_gmi(2, ng_id, e_con);

        if(!c2_map[ng_id]) {

          e_con->add_ent(e_lens.new_cell1_id.at(i));

          std::cout << "2c" << ng_id
                << " added 1c" << e_lens.new_cell1_id.at(i)
                << std::endl;

          e_lens.c_base_curr->set_conn_gmi(2, ng_id, e_con);
          c2_map[ng_id] = true;
        }
      }
      // Check if the exterior is on the path.
      std::vector<int>* c3_vec = &ng.cells.at(p_id).first;
      std::vector<int>::iterator c3_it;
      c3_it = std::find(c3_vec->begin(), c3_vec->end(), -1);

      if(c3_it != c3_vec->end()) {
        e_lens.c_base_curr->set_ext_gmi(1, e_lens.new_cell1_id.at(i), true);
      }
    }

    e_lens.c_base_curr->set_conn_gmi(2, e_lens.new_cell2_id, e_con2);

    if(e_lens.e_3c_t_flag) {
      e_lens.c_base_curr->get_conn_gmi(3, ng.cp.first, e_con);
      e_con->add_ent(e_lens.new_cell2_id);
      e_lens.c_base_curr->set_conn_gmi(3, ng.cp.first, e_con);

      std::cout << "3c" << ng.cp.first
              << " added 2c" << e_lens.new_cell2_id
              << std::endl;
    }

    if(e_lens.e_3c_b_flag) {
      e_lens.c_base_curr->get_conn_gmi(3, ng.cp.second, e_con);
      e_con->add_ent(e_lens.new_cell2_id);
      e_lens.c_base_curr->set_conn_gmi(3, ng.cp.second, e_con);

      std::cout << "3c" << ng.cp.second
              << " added 2c" << e_lens.new_cell2_id
              << std::endl;
    }

    if(!e_lens.e_3c_b_flag or !e_lens.e_3c_t_flag) {
      for(int i = 0; i < e_lens.new_cell1_id.size(); i++)
        e_lens.c_base_curr->set_ext_gmi(1, e_lens.new_cell1_id.at(i), true);
      for(int i = 0; i < e_lens.new_cell0_id.size(); i++)
        e_lens.c_base_curr->set_ext_gmi(0, e_lens.new_cell0_id.at(i), true);
      e_lens.c_base_curr->set_ext_gmi(2, e_lens.new_cell2_id, true);
    }

    if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL)
      determine_shell_gmi_wg();
  }
  cell_flag = true;

  delete e_con;
  delete e_con2;
}


void vd_edisc::recreate_slices() {
  // Iterate over the edges, which do not belong to a disc.
  // Recreate them and store them in the e_lens.create_map_rest.
  //std::cout << "Recreating the slice entities," 
  //          << " vert_ctr is " << vert_ctr_act
  //          << std::endl;
  apf::Downward down;
  apf::MeshEntity* e_temp;

  apf::MeshElement* ee;
  //double meas = 0.0;

  for (int i = 0; i < e_lens.es_edge.size(); i++) {
    int slice_id = e_lens.slice_map[e_lens.es_edge.at(i)];
    //std::cout << "Edge: " << e_lens.es_edge.e[i] <<
    //             ", slice " << slice_id << std::endl;
    // If the entity is on a slice, not on a disc, recreate it:
    if( slice_id > 0) {
      // Get the downward adjacencies, check if they are on a disc.
      e_lens.m->getDownward(e_lens.es_edge.at(i), 0, down);
      // Find the central vertex. Use slice map, as it will be extended for the
      // triangle and tet rebuilding.
      // The slice id of the first vertex.
      int vert_id = 0;
      int on_slice = e_lens.slice_map[down[vert_id]];
      int disc_id;
      // If it is on a disc, 
      //std::cout << "Vert " << down[0] << " id " << 0 
      //          << ", on_slice " << e_lens.slice_map[down[0]]
      //          << "Vert " << down[1] << " id " << 1 
      //          << ", on_slice " << e_lens.slice_map[down[1]] 
      //          << ",to be replaced by " 
      //          << e_lens.vert_ctr_new[slice_id-1]<< std::endl;

      if(down[vert_id] == vert_ctr_act) {
      }
      // If it is on a slice, try the other vertex. It should be the central 
      // one. 
      else {
        vert_id = 1;
      }

      down[vert_id] = e_lens.vert_ctr_new.at(slice_id-1);

      apf::ModelEntity* mdl = e_lens.m->toModel(e_lens.es_edge.at(i));
      apf::MeshEntity* edge_new = buildElement(e_lens.m, mdl, 
                                            apf::Mesh::EDGE, down, 0);
      //std::cout << "old edge " << e_lens.es_edge.at(i)
      //          << " new edge, " << edge_new << " "
      //        << e_lens.m->getModelType(mdl) << "c"
      //        << e_lens.m->getModelTag(mdl) << std::endl;

      //vd_print_vert(e_lens.m, e_lens.es_edge.e[i]);
      //vd_print_vert(e_lens.m, edge_new);
      e_lens.create_map_rest[e_lens.es_edge.at(i)] = edge_new;
    }
    //else {
      //std::cout << "skipped" << std::endl;
    //}

  }

  // Recreate the triangles.
  //std::cout << "Rebuilding the tris" << std::endl;
  for (int i = 0; i < e_lens.es_surf.size(); i++) {
    int slice_id = e_lens.slice_map[e_lens.es_surf.at(i)];
    //std::cout << "Tri " << e_lens.es_surf.at(i) 
    //          << " Slice id " << slice_id << std::endl;
    // If the entity is on a slice, not on a disc, recreate it:
    if( slice_id > 0) {
      // Get the downward adjacencies, check if they are on a disc.
      apf::Downward down;
      e_lens.m->getDownward(e_lens.es_surf.at(i), 0, down);
      int v_id = findIn(down, 3, vert_ctr_act);

      assert(v_id > -1);

      down[v_id] = e_lens.vert_ctr_new.at(slice_id-1);
      apf::ModelEntity* mdl = e_lens.m->toModel(e_lens.es_surf.at(i));
      apf::MeshEntity* tri_new = buildElement(e_lens.m, mdl, 
                                      apf::Mesh::TRIANGLE, down, 0);
      e_lens.create_map_rest[e_lens.es_surf.at(i)] = tri_new;
      //std::cout << "old tri " << e_lens.es_surf.at(i)
      //          << " new tri, " << tri_new
      //        << e_lens.m->getModelType(mdl) << "c"
      //        << e_lens.m->getModelTag(mdl) << std::endl;

      //vd_print_vert(e_lens.m, e_lens.es_surf.e[i]);
      //vd_print_vert(e_lens.m, tri_new);

    }

  }

  // Recreate the tets.
  //std::cout << "Recreating slice tets" << std::endl;

  for (int i = 0; i < e_lens.es_elem.size(); i++) {
    int slice_id = e_lens.slice_map[e_lens.es_elem.at(i)];
    //std::cout << "Slice id " << slice_id << std::endl;
    // If the entity is on a slice, not on a disc, recreate it:
    assert( slice_id > 0);
    // Create the new tet.
    apf::ModelEntity* mdl = e_lens.m->toModel(e_lens.es_elem.at(i));

    apf::Downward down;
    e_lens.m->getDownward(e_lens.es_elem.at(i), 0, down);
    int v_id = findIn(down, 4, vert_ctr_act);
    down[v_id] = e_lens.vert_ctr_new.at(slice_id-1);

    if(std::fabs(vd_volume_tet(e_lens.m, down)/6) 
                                < std::numeric_limits<double>::min()) {
      std::cout << "ins_flag: Tet volume is too small, inaccurate calc. "
                << std::endl;
      //ins_flag = false;
    }

    apf::MeshEntity* tet_new;
    tet_new = buildElement(e_lens.m, mdl, apf::Mesh::TET, down, 0);

    //tet_new = recreate_inv_elem(tet_new, mdl);


    //vd_print_vert(e_lens.m, e_lens.es_elem.e[i]);
    //vd_print_vert(e_lens.m, tet_new);


    //e_lens.m->getDownward(tet_new, 0, down);
    //std::cout << " Tet old " << e_lens.es_elem.at(i)
    //          << " tet new " << tet_new
    //          << " meas " << meas
    //          << " Volume is " << vd_volume_tet(e_lens.m, down) << std::endl;
    //if(vd_volume_tet(e_lens.m, down) < 0) {
    //if(meas < 0.) {
    //  std::cout << "Inverted ! ";
    //}

    //std::cout << "new surr tet " << tet_new << std::endl;
    e_lens.create_map_rest[e_lens.es_elem.at(i)] = tet_new;
    e_lens.slice_map[tet_new] = slice_id;
  }
  //std::cout << "Recreation successful." << std::endl;
}


// Destroy the old entities. Basically, most of the old entities have to be 
// destroyed. If one of the new vertices is the same as the old one, some of 
// the edges and triangles might be the same. 
// Starting from highest dimension, apply the following going over all 
// dimensions.
// For the entities around the 0cell vertex, check if they lie on a disc or 
// slice. If on a disc, check if any of the new disc entities are the same as 
// the old one. If not, destroy the entity. 
// If on a slice, check if the new entity is the same as the old one. If not,
// destroy.
void vd_edisc::destroy_ent() {
  std::cout << "Destroying old entities." << std::endl;
  std::cout << "n_elem: " << e_lens.es_elem.size()
            << " n_surf: " << e_lens.es_surf.size()
            << " n_edge: " << e_lens.es_edge.size()
            << std::endl;

  for (int i = 0; i < e_lens.es_elem.size(); i++) {
    apf::MeshEntity* ent_c = e_lens.es_elem.at(i);
    int slice_id = e_lens.slice_map[ent_c];

    //std::cout << ent_c << " " << e_lens.create_map_rest[ent_c] << std::endl; 
    if(slice_id > 0) {
      if(ent_c == e_lens.create_map_rest[ent_c]) {
      }
      else {
        //std::cout << "Destroyed" << std::endl; 
        e_lens.m->destroy(ent_c);
      }
    }

    else {
      assert(slice_id != 0);
    }

  }


  std::cout << "Elem destroyed." << std::endl;

  for (int i = 0; i < e_lens.es_surf.size(); i++) {
    apf::MeshEntity* ent_c = e_lens.es_surf.at(i);
    int slice_id = e_lens.slice_map[ent_c];
    if( slice_id < 0) {
      int disc_id = -slice_id-1;
      //std::cout << ent_c << " " 
      //          << e_lens.create_map[2*disc_id+1][ent_c] << " "
      //          << e_lens.create_map[2*disc_id][ent_c] << std::endl; 
      if(ent_c == e_lens.create_map[2*disc_id][ent_c]
         or ent_c == e_lens.create_map[2*disc_id+1][ent_c]) {
      }
      else {
        //std::cout << "Destroyed" << std::endl; 
        e_lens.m->destroy(ent_c);
      }

    }

    else if(slice_id > 0) {
      //std::cout << ent_c << " " << e_lens.create_map_rest[ent_c] 
      //          << std::endl; 
      if(ent_c == e_lens.create_map_rest[ent_c]) {
      }
      else {
        //std::cout << "Destroyed" << std::endl; 
        e_lens.m->destroy(ent_c);
      }
    }

    else {
      assert(slice_id != 0);
    }
  }

  std::cout << "Surf destroyed." << std::endl;

  for (int i = 0; i < e_lens.es_edge.size(); i++) {
    apf::MeshEntity* ent_c = e_lens.es_edge.at(i);
    int slice_id = e_lens.slice_map[ent_c];
    //std::cout << "Edge " << ent_c << " slice " << slice_id << std::endl;

    if( slice_id == -e_lens.discs.size()-1) {
      bool found = false;
      for(int j = 0; j < e_lens.discs.size()*2; j++) {
        if(ent_c == e_lens.create_map[j][ent_c]) {
          found = true;
          j = e_lens.discs.size();
        }
      }
      if(found) {
      }
      else {
        //std::cout << "Destroyed" << std::endl; 
        e_lens.m->destroy(ent_c);
      }
    }
    else if( slice_id < 0) {
      int disc_id = -slice_id-1;
      //std::cout << ent_c << " " 
      //          << e_lens.create_map[2*disc_id+1][ent_c] << " "
      //          << e_lens.create_map[2*disc_id][ent_c] << std::endl; 
      if(ent_c == e_lens.create_map[2*disc_id][ent_c]
         or ent_c == e_lens.create_map[2*disc_id+1][ent_c]) {
      }
      else {
        e_lens.m->destroy(ent_c);
      }
    }
    else if(slice_id > 0) {
      //std::cout << ent_c << " " << e_lens.create_map_rest[ent_c] << std::endl;
      if(ent_c == e_lens.create_map_rest[ent_c]) {
      }
      else {
        //std::cout << "Destroyed" << std::endl; 
        e_lens.m->destroy(ent_c);
      }
    }

    else {
      assert(slice_id != 0);
    }
  }

  std::cout << "Edge destroyed." << std::endl;

  //std::cout << "Vert old " << vert_ctr_act << std::endl;

  bool found = false;
  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    //std::cout << "Vert new " << e_lens.vert_ctr_new.at(i) << std::endl;
    if(vert_ctr_act == e_lens.vert_ctr_new.at(i))
      found = true;
  }
  if(!found){
    //std::cout << "Destroyed" << std::endl; 
    e_lens.m->destroy(vert_ctr_act);
  }
}

// Using the edge and 3cell model list, connect the edges of the disc or fin. 
void vd_edisc::get_disc_edges(circuit* circ_in, vd_disc* disc_in) {

  // TODO Get the actual 3cell indices. It is still confusing because the 2-
  // cell indices used belong to the graph...
  std::vector<int> c2_list(0);
  std::vector<int> c3_list(0);
  get_circ_topo_gmi(&c2_list, &c3_list, circ_in);

  ext_pc[0] = (std::find(c3_list.begin(), c3_list.end(), -1) != c3_list.end());
  //std::cout << "Getting the circuit disc" << std::endl;
  //for(int i = 0; i < c3_circ->conn.size(); i++) {
  //  std::cout << "3c" << c3_circ->conn.at(i) << " ";
  //}
  //for(int i = 0; i < c2_circ->conn.size(); i++) {
  //  std::cout << "2c" << c2_circ->conn.at(i) << " ";
  //}

  //vd_entlist_v ent_list(m_act, vert_ctr_act, c_base);
  //ent_list.print();

  disc_in->clear();
  disc_in->m = m_act;

  disc_in->v_ctr = vert_ctr_act;

  std::vector<apf::MeshEntity*> es_tri(0);

  apf::MeshEntity* e_2f;
  apf::MeshEntity* e_2e;

  tri_burn.clear();
  // Going over the path(circuit), connect two edges on subsequent 2cells. 
  // The mesh is preconditioned, so the edges will have a shortest path that 
  // goes through the 3cell connecting the 2cells.
  for(int i = 0; i < circ_in->first.first.size()-1; i++) {
    std::cout << "graph ids " << circ_in->first.second.at(i) << " " 
              << circ_in->first.second.at(i+1) << std::endl;

    e_2f = c2_edge.at(circ_in->first.second.at(i));
    e_2e = c2_edge.at(circ_in->first.second.at(i+1));
    int c2_f = m_act->getModelTag(m_act->toModel(e_2f));
    //int c3_f = circ_in->first.first.at(i);
    int c3_f = c3_list.at(i);
    int c2_e = m_act->getModelTag(m_act->toModel(e_2e));

    std::cout << "2c" << c2_f << " 3c" << c3_f << " 2c" << c2_e << std::endl;
    if (c3_f != -1) {
      std::cout << "e_2f " << e_2f << " e_2e " << e_2e << std::endl;
      assert(e_lens.vd_mesh_find_short(e_2f, e_2e, &es_tri, c3_f));
    }
  }

  e_2f = c2_edge.at(circ_in->first.second.back());
  e_2e = c2_edge.at(circ_in->first.second.at(0));
  int c2_f = m_act->getModelTag(m_act->toModel(e_2f));
  //int c3_f = circ_in->first.first.back();
  int c3_f = c3_list.back();
  int c2_e = m_act->getModelTag(m_act->toModel(e_2e));

  std::cout << "2c" << c2_f << " 3c" << c3_f << " 2c" << c2_e << std::endl;
  if (c3_f != -1) {
    std::cout << "e_2f " << e_2f << " e_2e " << e_2e << std::endl;
    assert(e_lens.vd_mesh_find_short(e_2f, e_2e, &es_tri, c3_f));
  }

  std::vector<std::vector<apf::MeshEntity*> > e_set
                 (4, std::vector<apf::MeshEntity*> (0) );

  e_set.at(2) = es_tri;
  vd_set_up(disc_in->m, &e_set.at(2), &e_set.at(3));
  vd_set_down(disc_in->m, &e_set.at(3), &e_set.at(2));
  vd_set_down(disc_in->m, &e_set.at(2), &e_set.at(1));
  vd_set_down(disc_in->m, &e_set.at(1), &e_set.at(0));

  if(sub_vtk)
    vd_save_vtk_set(disc_in->m, &e_set, "output/tri_disc");

  cp_tri(&es_tri, disc_in);

  //TODO This is obsolete, as the path or circuit goes through the 3cell, such
  // that both upper adjacencies have the same topology. Still, it works.
  // On the disc:
  // For each triangle, look at the neighboring tetrahedra. The tetrahedra  
  // belonging to the cone is the top tetrahedra, into which the lens will 
  // expand to is stored.
  // The other tet is stored in the bottom list.

  for (int i = 0; i < disc_in->tri.size(); i++) {
    apf::Up up;
    disc_in->m->getUp(disc_in->tri.at(i), up);
    assert(up.n == 2);

    // 3cell model.
    apf::ModelEntity* mdl = disc_in->m->toModel(up.e[0]);
    disc_in->t_em.at(i) = mdl;
  }

  // For each edge, look at the cell membership. Assuming preconditioning, it 
  // must belong to a 2cell or 3cell, and the associated lens triangle belongs 
  // to the cell same cell. 
  for (int i = 0; i < disc_in->edge.size(); i++) {
    apf::ModelEntity* e_m = disc_in->m->toModel(disc_in->edge.at(i));
    int type_e = m_act->getModelType(e_m);
    int tag_e = m_act->getModelTag(e_m);

    //std::cout << disc_in->edge.at(i) << " " << type_e << "c" << tag_e 
    //          << std::endl;

    // assert(0 < type_e and type_e < 3);

    disc_in->e_em.at(i) = e_m;
    assert(type_e > 1);
    disc_in->t2_em.at(i) = disc_in->e_em.at(i);
  }
}

// This always operates on the actual mesh m_main.
bool vd_edisc::set_0cell(int cell, bool spur_in) {
  assert(m_ex and m_load);
  spur = spur_in;

  cell_id = cell;

  std::vector<apf::MeshEntity*> es_vert_ctr(0);
  //vd_print_ent(m);
  vd_find_ent_geom(m_main, &es_vert_ctr, cell_id, 0);

  // TODO
  // As a check, assert that the cell actually exists and a mesh entity exists.
  std::cout << "Setting 0cell" << cell_id << ", v_sz:" 
            << es_vert_ctr.size() << std::endl;

  if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
    ext_shell* e_sh = f_calc.get_e_sh();

    if(c_base->get_cell_ext_gmi(0, cell)) {
      assert(e_sh->chk_shell(0, cell - 1));
      sh_ext = true;
      sh_old = e_sh->get_shell(0, cell - 1);
    }
    else {
      sh_ext = false;
      sh_old = shell(-1,-1);
    }
  }
  else {
    sh_ext = false;
    sh_old = shell(-1,-1);
  }
  std::cout << "Shell[" << sh_ext << "] " 
            << sh_old.dim << "sh" << sh_old.id << std::endl;

  if (es_vert_ctr.size() == 1) {
    set_vertex(es_vert_ctr.at(0));

    //std::vector<apf::MeshEntity*> ev(0);
    //std::vector<apf::MeshEntity*> ee(0);
    //std::vector<apf::MeshEntity*> es(0);
    //std::vector<apf::MeshEntity*> et(0);

    //init_ent_set(&ev, es_vert_ctr.at(0));

    // Get the tetrahedra:
    //vd_set_up(m_main, &ev, &ee);
    //vd_set_up(m_main, &ee, &es);
    //vd_set_up(m_main, &es, &et);
    //vd_set_down(m_main, &et, &es);
    //vd_set_down(m_main, &es, &ee);
    //vd_set_down(m_main, &ee, &ev);
/*
    for(int i = 0; i < ev.size(); i++) {
      int type = m_main->getModelType(m_main->toModel(ev.at(i) ));
      int tag = m_main->getModelTag(m_main->toModel(ev.at(i) ));
      apf::Vector3 pos;
      m_main->getPoint(ev.at(i), 0, pos);
      std::cout << "\t" << "Vert " 
              << ev.at(i) << ", " << type << "c" << tag << " ";
      std::cout << pos << std::endl;
    }
*/
    return true;
  }
  return false;

  // Copy the cell_base and reserve the actual cell_base.
}

// Set preconditioned mesh, initial configuration.
void vd_edisc::set_trial() {

  assert(m_ex and m_load);
  assert(precond_ex);

  if(precond_load) {
    apf::destroyMesh(m_precond);

    //m_precond = apf::makeEmptyMdsMesh(m->getModel(), 3, false);
    m_precond = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);
    precond_load = false;
  }

  int type_ctr = m_main->getModelType(vert_ctr_em);
  int tag_ctr = m_main->getModelTag(vert_ctr_em);

  trial_type = 0;

  // Get all surrounding entities:

  es_vert.clear();
  es_edge.clear();
  es_surf.clear();
  es_elem.clear();

  apf::Up up;
  apf::Downward v;

  init_ent_set(&es_vert, vert_ctr);

  // Get the tetrahedra:
  vd_set_up(m_main, &es_vert, &es_edge);
  vd_set_up(m_main, &es_edge, &es_surf);
  vd_set_up(m_main, &es_surf, &es_elem);

  vd_set_down(m_main, &es_elem, &es_surf);
  vd_set_down(m_main, &es_surf, &es_edge);

  if(sub_vtk) {
    std::vector<apf::MeshEntity*> es_vert_vtk(0);
    std::vector<apf::MeshEntity*> es_edge_vtk(0);
    std::vector<apf::MeshEntity*> es_surf_vtk(0);
    std::vector<apf::MeshEntity*> es_elem_vtk(0);

    init_ent_set(&es_vert_vtk, vert_ctr);

    // Get the tetrahedra:
    vd_set_up(m_main, &es_vert_vtk, &es_edge_vtk);
    vd_set_up(m_main, &es_edge_vtk, &es_surf_vtk);
    vd_set_up(m_main, &es_surf_vtk, &es_elem_vtk);

    // Save the entities before the transfer.
    vd_rem_tag(m_main);
    vd_tag_mesh(m_main);
    vd_tag_set(m_main, &es_elem_vtk, "Vertex_ent");

    apf::writeVtkFiles("./output/before_main", m_main);

    // Save the entities before the transfer.
    vd_save_vtk_vert(m_main, vert_ctr, "./output/before_main");
  }

  // Get the lower dimensional entities:
  vd_set_down(m_main, &es_elem, &es_surf);
  vd_set_down(m_main, &es_surf, &es_edge);
  vd_set_down(m_main, &es_edge, &es_vert);

  // These are to be copied to the new trial mesh:

  vert.resize(es_vert.size());

  std::cout << "Copied the surrounding entities, creating the vertices"
            << std::endl;

  main2pre_map.clear();
  pre2main_map.clear();
  // Copy the vertices to the new trial mesh.
  for (int i = 0; i < es_vert.size(); i++) {
    //std::cout << "Vert " << es_vert.at(i) << "["
    //          << m_main->toModel(es_vert.at(i)) << "] ," ;
    int type = m_main->getModelType(m_main->toModel(es_vert.at(i)));
    int tag = m_main->getModelTag(m_main->toModel(es_vert.at(i)));
    apf::ModelEntity* vert_em = m_precond->findModelEntity(type,tag);
    vert.at(i) = m_precond->createVert(vert_em);

    apf::MeshEntity* edge;
    if(es_vert.at(i) != vert_ctr and 
       vd_find_edge(m_main, es_vert.at(i), vert_ctr, &edge) and 
       m_main->getModelType(m_main->toModel(edge)) == 1) {
      std::cout << m_main->getModelType(m_main->toModel(edge)) << "c"
                << m_main->getModelTag(m_main->toModel(edge))
                << std::endl;
      main2pre_map[es_vert.at(i)] = vert.at(i);
      pre2main_map[vert.at(i)] = es_vert.at(i);
      std::cout << "Vert main " << es_vert.at(i) << " "
                << " precond " << vert.at(i) << std::endl;

    }

    apf::Vector3 vert_pos(0,0,0);
    m_main->getPoint(es_vert.at(i), 0, vert_pos);
    //std::cout<< vert_pos << std::endl;
    m_precond->setPoint(vert.at(i), 0, vert_pos);
  }

  // Copy the elements recursively:
  for (int i = 0; i < es_elem.size(); i++) {
    copyElement_pre(es_elem.at(i),0);
  }

  m_precond->acceptChanges();

  //for (int i = 0; i < 4; i++) {
    //apf::MeshIterator* it_e = m_precond->begin(i);
    //std::cout << i << " dim Iterator." << std::endl;
    //apf::MeshEntity* e;

    //std::cout << "Created ents: ";
    //while ((e = m_precond->iterate(it_e))) {
      //std::cout << e << ", ";
      //vd_ent_deg_adj(m_precond, e, 1);
    //}
    //std::cout << std::endl;

    //m_precond->end(it_e);
  //}

  std::vector<apf::MeshEntity*> es_vert_cp(0);
  // Get the central vertex in the (not-yet) preconditioned mesh.
  vd_find_ent_geom(m_precond, &es_vert_cp, cell_id, 0);
  assert(es_vert_cp.size() == 1);
  vert_ctr_precond = es_vert_cp.at(0);

  precond_load = true;

  act_precond();

  //vd_entlist_v ent_list(m_act, vert_ctr_act, c_base);
  //ent_list.print();

  if(sub_vtk)
    vd_save_vtk_vert(m_act, vert_ctr_act, "./output/before_precond");

  f_calc.vd_att_fields(m_precond);

  std::vector<apf::MeshEntity*> es_v_ctr(0);
  std::vector<apf::MeshEntity*> es_e_ctr(0);
  std::vector<apf::MeshEntity*> es_s_ctr(0);
  std::vector<apf::MeshEntity*> es_ee_ctr(0);

  init_ent_set(&es_v_ctr, vert_ctr_precond);

  // Get the tetrahedra:
  vd_set_up(m_precond, &es_v_ctr, &es_e_ctr);
  vd_set_up(m_precond, &es_e_ctr, &es_s_ctr);
  vd_set_up(m_precond, &es_s_ctr, &es_ee_ctr);

  if(len_sh < 0) {
    std::cout << "len_sh not set" << std::endl; 
    len_sh = vd_minmax_meas(m_precond, &es_ee_ctr, vert_ctr_precond, true)*3/4;
  }

  std::cout << "after initializing m_precond, len_sh = " << len_sh 
            << std::endl;

  vd_set_down(m_precond, &es_ee_ctr, &es_s_ctr);
  vd_set_down(m_precond, &es_s_ctr, &es_e_ctr);

  precond_mesh();
  m_precond->acceptChanges();
  if(pre_flag) {

    f_calc.vd_att_fields(m_precond);

    /////////////////////////////////////////////////////
    // Reload the vert list used in trial mesh generation.

    if(sub_vtk)
      vd_save_vtk_vert(m_act, vert_ctr_act, "./output/after_precond");

    std::vector<std::vector<apf::MeshEntity*> > es_pre
                 (4, std::vector<apf::MeshEntity*> (0) );

    init_ent_set(&es_pre[0], vert_ctr_act);

    // Get the tetrahedra:
    for(int i = 0; i < 3; i++)
      vd_set_up(m_act, &es_pre[i], &es_pre[i+1]);
    for(int i = 3; i > 0; i--)
      vd_set_down(m_act, &es_pre[i], &es_pre[i-1]);

    if(sub_vtk)
      vd_save_vtk_set(m_act, &es_pre, "output/pre_vert");

    vert.clear();
    vert.reserve(es_pre[0].size());

    apf::MeshIterator* it_e = m_act->begin(0);
    apf::MeshEntity* e;
    while ((e = m_act->iterate(it_e))) {
      vert.push_back(e);
      //std::cout << e << ", ";
      //vd_ent_deg_adj(m_precond, e, 1);
    }
    //std::cout << std::endl;

    m_act->end(it_e);

    init_ent_set(&es_v_ctr, vert_ctr_precond);

    // Get the tetrahedra:
    vd_set_up(m_precond, &es_v_ctr, &es_e_ctr);
    vd_set_up(m_precond, &es_e_ctr, &es_s_ctr);
    vd_set_up(m_precond, &es_s_ctr, &es_ee_ctr);

    std::cout << "after initializing m_precond, len_sp = "
              << len_sp << std::endl;

    vd_set_down(m_precond, &es_ee_ctr, &es_s_ctr);
    vd_set_down(m_precond, &es_s_ctr, &es_e_ctr);


  /*
    // Update the list of vertices to be copied to the new trial mesh.
    for (int i = 0; i < es_pre[0].n; i++) {
      vert[i] = es_pre[0].e[i];

      std::cout << "Vert " << es_pre[0].e[i] << "["
                << m_act->toModel(es_pre[0].e[i]) << "] ," ;

      apf::Vector3 vert_pos;
      m_act->getPoint(es_pre[0].e[i], 0, vert_pos);
      std::cout<< vert_pos << std::endl;
    }
  */
    collect_pc();
  }
}

bool vd_edisc::edge_comp(apf::MeshEntity* e1, apf::MeshEntity* e2, apf::Mesh2* m) {
  apf::ModelEntity* mdl = m->toModel(e1);
  int c1_type = m_act->getModelType(mdl);

  mdl = m->toModel(e2);

  // Highest cell membership first.
  if(c1_type != m_act->getModelType(mdl))
    return c1_type > m_act->getModelType(mdl);

  apf::MeshElement* ee = createMeshElement(m, e1);
  double len1 = measure(ee);
  destroyMeshElement(ee);

  ee = createMeshElement(m, e2);
  double len2 = measure(ee);
  destroyMeshElement(ee);

  // Widest element first.
  return len1 > len2;
}

void vd_edisc::sort_edge(std::vector<apf::MeshEntity*>* edge_list) { 
  std::sort(edge_list->begin(), edge_list->end(), 
      std::bind(&vd_edisc::edge_comp, this, 
      std::placeholders::_1, std::placeholders::_2, m_act) );
}
/*

void vd_edisc::sort_edge(std::vector<apf::MeshEntity*>* edge_list, int left, int right) {
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

int vd_edisc::partition_edge_list(std::vector<apf::MeshEntity*>* edge_list, int left, int right) {
  apf::MeshEntity* pivot = edge_list->at(right);
  // move the mid point value to the front.
  int i = left-1;
  int j = left;
  for (; j < right; j++) {
    if(edge_comp(edge_list->at(j),pivot, m_act)) {
      i++;
      std::swap(edge_list->at(i), edge_list->at(j));
    }
  }
  std::swap(edge_list->at(i+1),edge_list->at(right));
  return i + 1;
}
*/

void vd_edisc::chk_ma_new() {
  std::vector<apf::MeshEntity*> verts(0);
  if(trial_type == 2)
    verts.reserve(e_lens.vert_ctr_new.size() + e_lens.vert_sp_new.size()+1);
  else
    verts.reserve(e_lens.vert_ctr_new.size() + e_lens.vert_sp_new.size());

  std::copy(verts.begin(), verts.end(), std::back_inserter(e_lens.vert_ctr_new));
  verts.insert(verts.end(), e_lens.vert_sp_new.begin(), 
                              e_lens.vert_sp_new.end() );
  if(trial_type == 2)
    verts.push_back(e_lens.v_2c);

  if(e_lens.m == m_main)
    chk_ma_swap(e_lens.m);
  else
    chk_ma_swap(e_lens.m, verts);
}


// Precondition the mesh for 1cell and 2cell insertions.
// Precondition the input mesh for 2cell insertion. Preconditioning  
// involves making sure that there is an internal edge for each 2cell and  
// 3cell. This is to be called once on the trial mesh to prepare for the 
// trials and once for the actual mesh for splitting.

void vd_edisc::precond_mesh() {
  // Obtain the discs that seperate the 3cell couple tetrahedra from the rest 
  // of the entities. 
  // Expand the disc by an internal edge insertion for each 3cell. This will 
  // ensure that there will be at least a tetrahedron for each path. 

  // In order to make sure that there are at least two entities per 2cell and
  // a tet in a 3cell for each 2cell adjacency of that 3cell, count the number 
  // of entities in each 2- and 3-cell. Divide the entities until there are 
  // at least enough number of entities in each 2- and 3-cell.

  //std::vector<int> c3_list = c_ins->get_path_3cells(cell_id);
  ent_conn* e_adj = new ent_conn();
  c_base_act->get_conn_0d_gmi(cell_id, 3, e_adj);
  std::vector<int> c3_list(0);
  c3_list = e_adj->conn;
  delete e_adj;

  //std::vector<int> c3_list = c_ins->get_3cells_gmi(cell_id);
  std::cout << "3cells " << std::endl;
  trial_type = 3;

  init_elens();
  e_lens.c_base_curr = c_base_act;

  apf::Downward dw;
  std::cout << "Volumes before " << std::endl;
  for (int i = 0; i < e_lens.es_elem.size(); i++) {
    m_act->getDownward(e_lens.es_elem.at(i), 0, dw);
    //std::cout << " Tet old " << e_lens.es_elem.e[i] 
    //          << vd_volume_tet(m_act, dw) << std::endl;
    if(vd_volume_tet(m_act, dw) < 0) {
      std::cout << "Inverted ! ";
    }
  }
  double l_3cell;

  apf::ModelEntity* c3_g;
  bool precond_3c = true;
  // Reload the disjoint sets of tetrahedra until all of them have an edge 
  // inserted.
  while(precond_3c) {
    precond_3c = false;
    vd_3c->reload(m_act, c_base_act, vert_ctr_act);
    for (int i = 0; i < vd_3c->get_3c_nbr(); i++) {
      c3_g = vd_3c->get_mdl(3, i);

      std::vector<apf::MeshEntity*> tet_3c(0);
      tet_3c = vd_3c->get_tet_id(i);
      std::vector<apf::MeshEntity*> es_surf(0);
      std::vector<apf::MeshEntity*> es_edge(0);

      vd_set_down(e_lens.m, &tet_3c, &es_surf);
      vd_set_down(e_lens.m, &es_surf, &es_edge);

      int edge_nbr = 0;
      for(int j = 0; j < es_edge.size(); j++) {
        apf::ModelEntity* mdl = e_lens.m->toModel(es_edge.at(j));

        if(c3_g == mdl) {
          apf::Downward d_v;
          e_lens.m->getDownward(es_edge.at(j), 0, d_v);

          int j_v = findIn(d_v, 2, vert_ctr_act);
          if(j_v > -1) {
            edge_nbr = edge_nbr + 1;
          }
        }
        if(edge_nbr > 1)
          j = es_edge.size();
      }
      if(edge_nbr != 1) {
        std::cout << "Inserting edge inside 3c" << m_act->getModelTag(c3_g)
                  << std::endl;
        insert_3celledge(i);
        m_act->acceptChanges();

        //chk_ma_new();

        if(pre_flag) {
          //init_elens();
          //e_lens.c_base_curr = c_base_act;
          //ent_list.change_mesh(m_act, vert_ctr_act, c_base);
          //ent_list.print();
          init_elens();
          e_lens.c_base_curr = c_base_act;
          l_3cell = vd_minmax_meas(m_act, &e_lens.es_elem, 
                                                    vert_ctr_act, true)*5/6;
          std::cout << "after 3cell" << m_act->getModelTag(c3_g)
                    <<" edge insertions, l_3cell = " 
                    << l_3cell << std::endl;
          i = vd_3c->get_3c_nbr();
          precond_3c = true;
        }
        else {
          i = vd_3c->get_3c_nbr();
          precond_3c = false;
        }
      }
      // Check for triangles spanning through the disjoint 3cell, bounded by
      // lower dim edges adjacent to the center vertex or vertices at the end of
      // these edges. 
      else if(chk_span_surf(es_surf)) {
        std::cout << "Inserting edge inside 3c" << m_act->getModelTag(c3_g)
                  << std::endl;
        insert_3celledge(i);
        m_act->acceptChanges();

        //chk_ma_new();

        if(pre_flag) {
          //init_elens();
          //e_lens.c_base_curr = c_base_act;
          //ent_list.change_mesh(m_act, vert_ctr_act, c_base);
          //ent_list.print();
          init_elens();
          e_lens.c_base_curr = c_base_act;
          l_3cell = vd_minmax_meas(m_act, &e_lens.es_elem, 
                                                    vert_ctr_act, true)*5/6;
          std::cout << "after 3cell" << m_act->getModelTag(c3_g)
                    <<" edge insertions, l_3cell = " 
                    << l_3cell << std::endl;
          i = vd_3c->get_3c_nbr();
          precond_3c = true;
        }
        else {
          i = vd_3c->get_3c_nbr();
          precond_3c = false;
        }
      }
    }
  }

/*
  for (int i = 0; i < c3_list.size(); i++) {
    std::cout << "Inserting edge inside 3c" << c3_list.at(i) << std::endl;
    insert_3celledge(c3_list.at(i));
    m_act->acceptChanges();
    if(pre_flag) {
      //init_elens();
      //e_lens.c_base_curr = c_base_act;
      //ent_list.change_mesh(m_act, vert_ctr_act, c_base);
      //ent_list.print();
      init_elens();
      e_lens.c_base_curr = c_base_act;
      l_3cell = vd_minmax_meas(m_act, &e_lens.es_elem, 
                                                vert_ctr_act, true)*5/6;
      std::cout << "after 3cell" << c3_list.at(i) 
                <<" edge insertions, l_3cell = " 
                << l_3cell << std::endl;
    }
    else
      i = c3_list.size();
  }
*/
  if(pre_flag) {
    std::cout << "Volumes after 3cell " << std::endl;
    assert(vd_chk_neg_sgn(m_act) == 0);

    if(sub_vtk)
      vd_save_vtk_vert(m_act, vert_ctr_act, "./output/after_3cell");
    init_elens();
    e_lens.c_base_curr = c_base_act;

    vd_entlist_v ent_list(m_act, vert_ctr_act, c_base);

    std::vector<apf::MeshEntity* > tri_2_sp(0);
    collect_split_tri(&tri_2_sp, &ent_list);

    vd_bipy* tri_split = new vd_bipy(m_act, c_base, &f_calc);

    apf::Vector3 temp(0,0,0);
    for(int i = 0; i < tri_2_sp.size(); i++) {
      std::cout << "tri_2_sp " << tri_2_sp.at(i) << std::endl; 
      temp = vd_get_pos(m_act, tri_2_sp.at(i));

      tri_split->load_tri(tri_2_sp.at(i));
      tri_split->split_bipy();
      m_act->acceptChanges();
      vd_chk_neg_sgn(m_act);
      apf::MeshEntity* ev_sp = tri_split->get_vert_ctr();
      //f_calc.vd_att_fields(m_act, ev_sp);

      m_act->setPoint(ev_sp, 0, pos_old + (temp - pos_old)*1.2);
    }
    delete tri_split;
    m_act->acceptChanges();

    ent_list.change_mesh(m_act, vert_ctr_act, c_base);
    // init_elens
    init_elens();
    e_lens.c_base_curr = c_base_act;
    if(sub_vtk)
      vd_save_vtk_vert(m_act, vert_ctr_act, "./output/after_cross");

    assert(e_lens.init_flag);

    vd_rem_tag(m_act);

    // Split the edges bounded by the central vertex and project the new vertices
    // on a sphere.
    std::pair<double, double> l_sph = vert_dist_sphere(m_act, vert_ctr_act);
    if(len_sh > std::numeric_limits<double>::min())
      l_min = std::min(l_sph.first*5/6, len_sh*3/4*5/6);

    std::vector<apf::MeshEntity* > ev(0);
    std::vector<apf::MeshEntity* > ee(0);

    vd_set_up(m_act, vert_ctr_act, &ee);
    ev.resize(ee.size());

    //double l_min = vd_minmax_meas(m_act, &ee)*5/6;
    //double l_min = vd_minmax_meas(m_act, &e_lens.es_elem, 
    //                                          vert_ctr_act, true)*5/6;

    //if(l_min < std::numeric_limits<double>::min())

    vd_lens* lens_sph = new vd_lens(m_act, c_base, &f_calc);

    for (int i = 0; i < ee.size(); i++) {
      //std::cout << "e2_sp " << ee.at(i) << std::endl; 

      apf::MeshEntity* v_curr;
      apf::MeshEntity* v_temp;
      int c_type = m_act->getModelType(m_act->toModel(ee.at(i)));
      apf::Downward down;
      m_act->getDownward(ee.at(i), 0, down);
      int v1 = findIn(down, 2, vert_ctr_act);
      assert(v1 > -1);
      v1 = (v1+1) % 2;

      if(c_type == 1) {
        v_curr = down[v1];
        //std::cout << c_type << "c"
        //          << m_act->getModelTag(m_act->toModel(ee.at(i)))
        //          << " " << v_curr << std::endl;
      }

      lens_sph->load_edge(ee.at(i));
      lens_sph->split_lens(
                  getEdgeVertOppositeVert(m_act, ee.at(i), vert_ctr_act));
      m_act->acceptChanges();

      if(c_type == 1) {
        if(m_act == m_precond) {
          //std::cout << "preid_map[old] " << preid_map[v_curr]
          //          << " new_vert " << lens_sph->get_vert_ctr();

          //preid_map[lens_sph->get_vert_ctr()] = preid_map[v_curr];
          v_temp = pre2main_map[v_curr];
          pre2main_map[lens_sph->get_vert_ctr()] = v_temp;
          main2pre_map[v_temp] = lens_sph->get_vert_ctr();

          //std::cout << " pre2main_map[old] " << pre2main_map[v_curr]
          //          << " main2pre_map[v_temp] " << main2pre_map[v_temp]
          //          << std::endl;

        }
        else {
          assert(m_act != m_trial);
          assert(m_act == m_main);
          v_temp = main2pre_map[v_curr];
          pre2main_map[v_temp] = lens_sph->get_vert_ctr();
          main2pre_map[lens_sph->get_vert_ctr()] = v_temp;

          //std::cout << " main2pre_map[v_curr] " << main2pre_map[v_curr]
          //          << " pre2main_map[v_temp] " << pre2main_map[v_temp]
          //          << std::endl;
        }
      }

      vd_chk_neg_sgn(m_act);
      //f_calc.vd_att_fields(m_act, lens_sph->get_vert_ctr());

      ev.at(i) = lens_sph->get_vert_ctr();
    }

    delete lens_sph;

    std::cout << "Volumes after split " << std::endl;
    vd_chk_neg_sgn(m_act);
    if(sub_vtk)
      vd_save_vtk_vert(m_act, vert_ctr_act, "./output/after_split");

    apf::Vector3 ctr_pos(0,0,0);
    m_act->getPoint(vert_ctr_act, 0, ctr_pos);

    for (int i = 0; i < ev.size(); i++) {
      vd_proj_v_sphere(m_act, ev.at(i), ctr_pos, l_min);
      vd_chk_neg_sgn(m_act);
    }

    if(sub_vtk)
      vd_save_vtk_vert(m_act, vert_ctr_act, "./output/after_sphere");

    //vd_entlist_v ent_list(m_act, vert_ctr_act, c_base);
    //ent_list.print();

    // For each 3cell, get the 2cell adjacencies. For each triangle of each 
    // 2cell adjacency, there should be a tet belonging to the 3cell.

    // Instead, check every tet triangles. If a tet has two triangle faces 
    // belonging to 2cells in the cells2 list, mark the edge not contained by
    // any of the two triangles. 
    // If a tet belongs to a 3cell with a single tet inside, that tet should be
    // split into 3 tets. 
    // Split the marked edges.

    m_act->acceptChanges();

    init_elens();
    e_lens.c_base_curr = c_base_act;

    std::cout << "Volumes after sphere " << std::endl;
    vd_chk_neg_sgn(m_act);

    //ent_list.change_mesh(m_act, vert_ctr_act, c_base);
    //ent_list.change_mesh(m_act, vert_ctr_act, c_base);
    //ent_list.print();
    //len_sp = vd_minmax_meas(m_act, &e_lens.es_elem, vert_ctr_act, true)*3/4;
    len_sp = len_sh*3/4;
    len_edge = vd_minmax_meas(m_act, &e_lens.es_edge, true)/3*4;

    std::cout << "after preconditioning m_precond, len_sp = " 
              << len_sp << ", len_edge = " << len_edge << std::endl;
  }
}

void vd_edisc::collect_split_tri(std::vector<apf::MeshEntity* >* tri_2_sp, 
                                               vd_entlist_v* ent_list) {

  apf::Downward down;

  tri_2_sp->clear();
  tri_2_sp->reserve(ent_list->es.at(2).size());

  std::cout << "Preconditioning the mesh " << std::endl; 

  // Check the edges of each triangle around the center vertex, which are also 
  // adjacent to the central vertex. If both of them belong to a lower 
  // dimensional cell, the edge across should be split.
  for(int i = 0; i < ent_list->es.at(2).size(); i++) {
    apf::ModelEntity* mdl = m_act->toModel(ent_list->es.at(2).at(i));
    int c2_type = m_act->getModelType(mdl);

    m_act->getDownward(ent_list->es.at(2).at(i), 0, down);
    int v1 = findIn(down, 3, vert_ctr_act);

    if (v1 > -1) {
      m_act->getDownward(ent_list->es.at(2).at(i), 1, down);

      mdl = m_act->toModel(down[lookup_tri_ed[v1][0]]);
      int c1_type = m_act->getModelType(mdl);

      if(c1_type < c2_type) {
        mdl = m_act->toModel(down[lookup_tri_ed[v1][1]]);
        int c1_type = m_act->getModelType(mdl);
        if(c1_type < c2_type) {
          tri_2_sp->push_back(ent_list->es.at(2).at(i));
          //std::cout << "splitting " << down[lookup_v_tri_e_x[v1]] 
          //          << std::endl; 
        }
      }

      // Also, over all triangles neighboring the center vertex, check the 
      // counter edge. If both of the end vertices belong to cells different
      // (lower dim) than the edge cell, split the triangle, as well.
      apf::MeshEntity* edge_cross = down[lookup_v_tri_e_x[v1]];
      mdl = m_act->toModel(edge_cross);
      c1_type = m_act->getModelType(mdl);
      m_act->getDownward(edge_cross, 0, down);
/*
      std::cout 
        << "\te_cross " << edge_cross << " "
        << m_act->getModelType(m_act->toModel(edge_cross))
        << "c" << m_act->getModelTag(m_act->toModel(edge_cross))
        << " v1 " << down[0] << " "
        << m_act->getModelType(m_act->toModel(down[0]))
           << "c" << m_act->getModelTag(m_act->toModel(down[0]))
        << " v2 "  << down[1] << " "
        << m_act->getModelType(m_act->toModel(down[1]))
            << "c" << m_act->getModelTag(m_act->toModel(down[1]))
        << std::endl;
*/
      if(c1_type > m_act->getModelType(m_act->toModel(down[0])) ) {
        if(c1_type > m_act->getModelType(m_act->toModel(down[1])) ) {
          tri_2_sp->push_back(ent_list->es.at(2).at(i));
        }
      }
    }
  }

  //std::cout << "Tris to split before unique " << std::endl; 
  for(int i = 0; i < tri_2_sp->size(); i++) {
    apf::ModelEntity* mdl = m_act->toModel(tri_2_sp->at(i));
    int c1_type = m_act->getModelType(mdl);
    int c1_tag = m_act->getModelTag(mdl);
    std::cout << tri_2_sp->at(i) << "("<< c1_type << "c"<< c1_tag << ") "; 
  }

  std::sort(tri_2_sp->begin(), tri_2_sp->end());
  std::vector<apf::MeshEntity*>::iterator it;
  it = std::unique (tri_2_sp->begin(), tri_2_sp->end());
  tri_2_sp->resize(std::distance(tri_2_sp->begin(),it));

  //std::cout << "Tris to split for preconditioning " << std::endl; 
  for(int i = 0; i < tri_2_sp->size(); i++) {
    apf::ModelEntity* mdl = m_act->toModel(tri_2_sp->at(i));
    int c1_type = m_act->getModelType(mdl);
    int c1_tag = m_act->getModelTag(mdl);
    std::cout << tri_2_sp->at(i) << "("<< c1_type << "c"<< c1_tag << ") "; 
  }
  std::cout << std::endl;
}

void vd_edisc::copy_trial_wg(apf::MeshEntity* v_new, 
                             std::vector<apf::MeshEntity*> &e_set, 
                             std::map<apf::MeshEntity*,apf::MeshEntity* > &p2t, 
                             std::map<apf::MeshEntity*,apf::MeshEntity* > &t2p) {

  for(int i = 0; i < e_set.size(); i++) {
    copyElement_rpl(e_set.at(i), e_lens.vert_ctr, v_new, p2t, t2p, 0);
  }

}

void vd_edisc::create_trial_vert_wg(std::vector<apf::MeshEntity*> &e_set, 
                             std::map<apf::MeshEntity*,apf::MeshEntity* > &p2t, 
                             std::map<apf::MeshEntity*,apf::MeshEntity* > &t2p) {

  apf::Downward down;
  for(int i = 0; i < e_set.size(); i++) {
    apf::MeshEntity* v_old;
    e_lens.m->getDownward(e_set.at(i), 0, down);
    if(down[0] == e_lens.vert_ctr) {
      v_old = down[1];
    }
    else {
      assert(down[1] == e_lens.vert_ctr);
      v_old = down[0];
    }
    apf::ModelEntity* vert_em = e_lens.m->toModel(v_old);
    int type = e_lens.m->getModelType(vert_em);
    int tag = e_lens.m->getModelTag(vert_em);
    vert_em = m_trial->findModelEntity(type, tag);
    apf::Vector3 pos(0,0,0);

    e_lens.m->getPoint(v_old, 0, pos);
    apf::MeshEntity* v_new = m_trial->createVert(vert_em);
    m_trial->setPoint(v_new, 0, pos);
    p2t[v_old] = v_new;
    t2p[v_new] = v_old;
  }

}

void vd_edisc::copy_slice_ents_wg(int slice, 
            std::vector<apf::MeshEntity*> & v_ctr_new,
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_ctr,
                std::map<apf::MeshEntity*, int> &s_map,
                std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells) {
  std::vector<apf::MeshEntity*> temp(0);
  std::vector<std::vector<apf::MeshEntity*> > e_slice(3, 
                                      std::vector<apf::MeshEntity*> (0));

  e_slice.at(0).reserve(e_lens.es_edge.size());
  e_slice.at(1).reserve(e_lens.es_surf.size());
  e_slice.at(2).reserve(e_lens.es_elem.size());

  apf::Downward down;
  for(int dim = 0; dim < 3; dim++) {
    for(int i = 0; i < slice_cells->at(slice).second.at(dim).size(); i++) {
      int c_local = slice_cells->at(slice).second.at(dim).at(i);
      temp = vd_cd->get_dim_ent(dim+1, c_local);
      std::copy(temp.begin(), temp.end(), std::back_inserter(e_slice.at(dim)));
    }
  }
  // Collect all the edges on the slice bounded by the center vertex.
  std::vector<apf::MeshEntity*>::iterator it;
  for(int dim = 2; dim > 0; dim--) {
    vd_set_down(e_lens.m, &e_slice.at(dim), &temp);
    std::copy(temp.begin(), temp.end(), std::back_inserter(e_slice.at(dim-1)));
    std::sort(e_slice.at(dim-1).begin(), e_slice.at(dim-1).end());
    it = std::unique (e_slice.at(dim-1).begin(), e_slice.at(dim-1).end());
    e_slice.at(dim-1).resize(std::distance(e_slice.at(dim-1).begin(),it));
  }
  temp.clear();
  temp.reserve(e_slice.at(0).size());
  for(int i = 0; i < e_slice.at(0).size(); i++) {
    e_lens.m->getDownward(e_slice.at(0).at(i), 0, down);
    if(down[0] == e_lens.vert_ctr or down[1] == e_lens.vert_ctr)
      temp.push_back(e_slice.at(0).at(i));
  }
  create_trial_vert_wg(temp, p2t_ctr.at(slice), t2p_ctr.at(slice));

  vd_set_up(e_lens.m, &temp, &e_slice.at(1));
  vd_set_up(e_lens.m, &e_slice.at(1), &e_slice.at(2));

  for(int i = 0; i < e_slice.at(2).size(); i++) {
    s_map[e_slice.at(2).at(i)] = slice + 1;
  }
  // Copy the disc end vertices, edges and tris:
  if(slice_cells->at(slice).first.first == 0) {
    create_trial_vert_wg(e_lens.discs.at(0).edge, 
                                  p2t_ctr.at(slice), t2p_ctr.at(slice));

    copy_trial_wg(v_ctr_new.at(slice), e_lens.discs.at(0).edge, 
                                    p2t_ctr.at(slice),t2p_ctr.at(slice));
    copy_trial_wg(v_ctr_new.at(slice), e_lens.discs.at(0).tri, 
                                    p2t_ctr.at(slice),t2p_ctr.at(slice));
  }
  else {
    int d1 = slice_cells->at(slice).first.first - 1;
    int d2 = slice_cells->at(slice).first.second - 1;
    create_trial_vert_wg(e_lens.discs.at(d1).edge, 
                                  p2t_ctr.at(slice), t2p_ctr.at(slice));
    create_trial_vert_wg(e_lens.discs.at(d2).edge, 
                                  p2t_ctr.at(slice), t2p_ctr.at(slice));

    copy_trial_wg(v_ctr_new.at(slice), e_lens.discs.at(d1).edge, 
                                    p2t_ctr.at(slice),t2p_ctr.at(slice));
    copy_trial_wg(v_ctr_new.at(slice), e_lens.discs.at(d1).tri, 
                                    p2t_ctr.at(slice),t2p_ctr.at(slice));

    copy_trial_wg(v_ctr_new.at(slice), e_lens.discs.at(d2).edge, 
                                    p2t_ctr.at(slice),t2p_ctr.at(slice));
    copy_trial_wg(v_ctr_new.at(slice), e_lens.discs.at(d2).tri, 
                                    p2t_ctr.at(slice),t2p_ctr.at(slice));
  }

  // Copy the slice edges, tris and tets:
  for(int dim = 0; dim < 3; dim++) {
    copy_trial_wg(v_ctr_new.at(slice), e_slice.at(dim), 
                                    p2t_ctr.at(slice), t2p_ctr.at(slice));
  }
}

// Recreate the path edges, tris for a disc on 1cell vertices.
void vd_edisc::create_path_int_p_wg(int path, int slice,
                    std::vector<apf::MeshEntity*> & v_sp_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_sp,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_sp,
                    std::vector<std::map<int, apf::MeshEntity*> > &p2sv_map) {
  apf::Downward down;

  // 1-cell sphere, the vertex corresponding to the 0-cell vertex is going to be
  // 1cell as the sphere represents the infinitesimal neigborhood:
  down[0] = v_sp_new.at(path);
  int tag = e_lens.new_cell1_id.at(path);
  apf::ModelEntity* edge_em = m_trial->findModelEntity(1,tag);
  down[1] = m_trial->createVert(edge_em);
  m_trial->setPoint(down[1], 0, pos_old);

  p2sv_map.at(path)[slice] = down[1];
  buildElement(m_trial, edge_em, apf::Mesh::EDGE, down, 0);
  
  // Split disc edges and internal triangles:
  for(int i = 0; i < e_lens.discs.at(path).edge.size(); i++) {
    apf::MeshEntity* e_curr = e_lens.discs.at(path).edge.at(i);
    edge_em = e_lens.m->toModel(e_curr);
    int type = e_lens.m->getModelType(edge_em);
    tag = e_lens.m->getModelTag(edge_em);
    edge_em = m_trial->findModelEntity(type,tag);

    e_lens.m->getDownward(e_curr, 0, down);
    int v1 = findIn(down, 2, e_lens.vert_ctr);
    assert(v1 > -1);

    down[v1] = p2sv_map.at(path)[slice];
    int j_curr = (v1+1)%2;
    down[j_curr] = p2t_sp.at(path)[down[j_curr]];
    assert(down[j_curr]);
    // Map the copy. Mainly used to obtain the copy of disconnected 3-cell edges
    // on the 1-cell spheres.
    p2t_sp.at(path)[e_curr] = buildElement(m_trial, edge_em, apf::Mesh::EDGE, 
                                                                      down, 0);
    t2p_sp.at(path)[p2t_sp.at(path)[e_curr]] = e_curr;
    down[2] = v_sp_new.at(path);
    buildElement(m_trial, edge_em, apf::Mesh::TRIANGLE, down, 0);
  }
}

// Recreate the path edges, tris for a disc on slices connected to the paths.
// The disc edges and the corresponding internal triangles are generated.
// The other edges bounding the internal triangles are already generated when 
// copying slices.
void vd_edisc::create_path_int_s_wg(int path, int slice,
                    std::vector<apf::MeshEntity*> & v_ctr_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_ctr,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map) {
  apf::Downward down;
  // 0-cell sphere:
  down[0] = v_ctr_new.at(slice);
  int tag = e_lens.new_cell1_id.at(path);
  apf::ModelEntity* edge_em = m_trial->findModelEntity(1,tag);
  down[1] = m_trial->createVert(edge_em);
  m_trial->setPoint(down[1], 0, pos_old);

  s2pv_map.at(slice)[path] = down[1];
  buildElement(m_trial, edge_em, apf::Mesh::EDGE, down, 0);
  // Split disc edges and internal triangles:
  for(int i = 0; i < e_lens.discs.at(path).edge.size(); i++) {
    apf::MeshEntity* e_curr = e_lens.discs.at(path).edge.at(i);
    edge_em = e_lens.m->toModel(e_curr);
    int type = e_lens.m->getModelType(edge_em);
    tag = e_lens.m->getModelTag(edge_em);
    edge_em = m_trial->findModelEntity(type,tag);

    e_lens.m->getDownward(e_curr, 0, down);
    int v1 = findIn(down, 2, e_lens.vert_ctr);
    assert(v1 > -1);

    down[v1] = s2pv_map.at(slice)[path];
    int j_curr = (v1+1)%2;
    down[j_curr] = p2t_ctr.at(slice)[down[j_curr]];
    assert(down[j_curr]);

    // These edges already have a copy on the slice. 
    buildElement(m_trial, edge_em, apf::Mesh::EDGE, down, 0);
    down[2] = v_ctr_new.at(slice);
    buildElement(m_trial, edge_em, apf::Mesh::TRIANGLE, down, 0);
  }

}

void vd_edisc::create_path_int_s_tet_wg(int path, apf::MeshEntity* e_curr,
                    apf::MeshEntity* tet, 
                    std::vector<apf::MeshEntity*> & v_ctr_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_ctr,
            std::map<apf::MeshEntity*, int> &s_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map,
                    std::vector<std::vector<apf::MeshEntity*> > &s_tets) {
  apf::Downward down;
  apf::Downward d_t;

  // Given two indices, return the other two:
  int ind4others[4][4][2] = {{{-1, -1}, {2, 3}, {1, 3}, {1, 2}},
                         {{2, 3}, {-1, -1}, {0, 3}, {0, 2}}, 
                         {{1, 3}, {0, 3}, {-1, -1}, {0, 1}},
                         {{1 ,2}, {0, 2}, {0, 1}, {-1, -1}}};

  int s_id = s_map[tet] - 1;

  apf::ModelEntity* tri_em = e_lens.m->toModel(e_curr);
  int type = e_lens.m->getModelType(tri_em);
  int tag = e_lens.m->getModelTag(tri_em);
  tri_em = m_trial->findModelEntity(type,tag);

  e_lens.m->getDownward(e_curr, 0, down);
  int v1 = findIn(down, 3, e_lens.vert_ctr);
  assert(v1 > -1);

  down[v1] = s2pv_map.at(s_id)[path];
  for(int j = 1; j < 3; j++) {
    int j_curr = (v1+j)%3;
    down[j_curr] = p2t_ctr.at(s_id)[down[j_curr]];
    assert(down[j_curr]);
  }

  buildElement(m_trial, tri_em, apf::Mesh::TRIANGLE, down, 0);

  e_lens.m->getDownward(tet, 0, down);
  e_lens.m->getDownward(tet, 2, d_t);
  v1 = findIn(down, 4, e_lens.vert_ctr);
  int t1 = findIn(d_t, 4, e_curr);
  assert(v1 > -1 and t1 > -1);
  int v2 = lookup_tet_x_surf[t1];
  // Invert the vertex ordering: v1 to be replaced with the new 1-cell vertex
  // and v2 to be replaced with new 0-cell vertex.
  down[v1] = s2pv_map.at(s_id)[path];
  down[v2] = v_ctr_new.at(s_id);
  for(int j = 0; j < 2; j++) {
    int j_curr = ind4others[v1][v2][j];
    down[j_curr] = p2t_ctr.at(s_id)[down[j_curr]];
    assert(down[j_curr]);
  }

  apf::MeshEntity* tet_new = buildElement(m_trial, tri_em, apf::Mesh::TET, down, 0);
  s_tets.at(s_id).push_back(tet_new);
}

apf::MeshEntity* vd_edisc::create_path_int_p_tet_wg(int path, 
                                        apf::MeshEntity* e_curr,
                    apf::MeshEntity* tet, 
                    std::vector<apf::MeshEntity*> & v_sp_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_sp,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_sp,
            std::map<apf::MeshEntity*, int> &s_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &p2sv_map) {
  apf::Downward down;
  apf::Downward d_t;

  // Given two indices, return the other two:
  int ind4others[4][4][2] = {{{-1, -1}, {2, 3}, {1, 3}, {1, 2}},
                         {{2, 3}, {-1, -1}, {0, 3}, {0, 2}}, 
                         {{1, 3}, {0, 3}, {-1, -1}, {0, 1}},
                         {{1 ,2}, {0, 2}, {0, 1}, {-1, -1}}};

  int s_id = s_map[tet] - 1;

  apf::ModelEntity* tri_em = e_lens.m->toModel(e_curr);
  int type = e_lens.m->getModelType(tri_em);
  int tag = e_lens.m->getModelTag(tri_em);
  tri_em = m_trial->findModelEntity(type,tag);

  e_lens.m->getDownward(e_curr, 0, down);
  int v1 = findIn(down, 3, e_lens.vert_ctr);
  assert(v1 > -1);

  down[v1] = p2sv_map.at(path)[s_id];
  for(int j = 1; j < 3; j++) {
    int j_curr = (v1+j)%3;
    down[j_curr] = p2t_sp.at(path)[down[j_curr]];
    assert(down[j_curr]);
  }

  buildElement(m_trial, tri_em, apf::Mesh::TRIANGLE, down, 0);

  e_lens.m->getDownward(tet, 0, down);
  e_lens.m->getDownward(tet, 2, d_t);
  v1 = findIn(down, 4, e_lens.vert_ctr);
  int t1 = findIn(d_t, 4, e_curr);
  assert(v1 > -1 and t1 > -1);
  int v2 = lookup_tet_x_surf[t1];
  // v1 to be replaced with the new 0-cell vertex
  // and v2 to be replaced with new 1-cell vertex.
  down[v1] = v_sp_new.at(path);
  down[v2] = p2sv_map.at(path)[s_id];
  for(int j = 0; j < 2; j++) {
    int j_curr = ind4others[v1][v2][j];
    down[j_curr] = p2t_sp.at(path)[down[j_curr]];
    assert(down[j_curr]);
  }

  return buildElement(m_trial, tri_em, apf::Mesh::TET, down, 0);
}


// Create the 3-cell tetrahedra for all spheres.
// Using the vertex ordering of the tet on slices, replace the following:
// the center vertex by the 2-cell vertex or it's copy
// the vertex across triangle by the 0-cell vertex or it's copy
// the 3-cell vertex by the copy of the 3-cell vertex
// the last vertex by the 1-cell vertex 
void vd_edisc::create_3c_tet_wg(int path, apf::MeshEntity* e_curr,
                    apf::MeshEntity* tet, 
                    std::vector<apf::MeshEntity*> & v_ctr_new,
                    std::vector<apf::MeshEntity*> & v_sp_new, 
                    apf::MeshEntity* v2c_new, apf::MeshEntity* v3c,
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_sp,
            std::map<apf::MeshEntity*,apf::MeshEntity* > &p2c_2c,
                    std::vector<std::map<int, apf::MeshEntity*> > &p2sv_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map,
                    std::vector<apf::MeshEntity*> &s2cv_map,
                    std::map<int, apf::MeshEntity*> &c2sv_map,
                    std::vector<apf::MeshEntity*> &p2cv_map,
                    std::map<int, apf::MeshEntity*> &c2pv_map,
                    std::map<apf::MeshEntity*, int> &s_map,
                    std::vector<std::vector<apf::MeshEntity*> > &s_tets) {
  apf::Downward down;
  apf::Downward d_t;

  int s_id = s_map[tet] - 1;

  apf::ModelEntity* tet_em = e_lens.m->toModel(tet);
  int tag = e_lens.m->getModelTag(tet_em);
  tet_em = m_trial->findModelEntity(3, tag);

  e_lens.m->getDownward(tet, 0, down);
  e_lens.m->getDownward(tet, 2, d_t);
  int v1 = findIn(down, 4, e_lens.vert_ctr);
  int t1 = findIn(d_t, 4, e_curr);
  assert(v1 > -1 and t1 > -1);
  int v2 = lookup_tet_x_surf[t1];
  int v3 = findIn(down, 4, v3c);

  int v4 = 6 - v1 - v2 - v3;

  // 0-cell sphere
  down[v1] = s2cv_map.at(s_id);
  down[v2] = v_ctr_new.at(s_id);
  down[v3] = p2t_ctr.at(s_id)[v3c];
  down[v4] = s2pv_map.at(s_id)[path];
  apf::MeshEntity* temp = buildElement(m_trial, tet_em, apf::Mesh::TET, down, 0);
  s_tets.at(s_id).push_back(temp);
  // 1-cell sphere
  down[v1] = p2cv_map.at(path);
  down[v2] = p2sv_map.at(path)[s_id];
  down[v3] = p2t_sp.at(path)[v3c];
  down[v4] = v_sp_new.at(path);
  buildElement(m_trial, tet_em, apf::Mesh::TET, down, 0);
  // 2-cell sphere
  down[v1] = v2c_new;
  down[v2] = c2sv_map.at(s_id);
  down[v3] = p2c_2c[v3c];
  down[v4] = c2pv_map.at(path);
  buildElement(m_trial, tet_em, apf::Mesh::TET, down, 0);
}



void vd_edisc::create_path_s_wg(int path, apf::MeshEntity* e_curr,
                    std::vector<apf::MeshEntity*> & v_ctr_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_ctr,
            std::map<apf::MeshEntity*, int> &s_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map,
                    std::vector<std::vector<apf::MeshEntity*> > &s_tets) {

  apf::Up up;
  e_lens.m->getUp(e_curr, up);
  create_path_int_s_tet_wg(path, e_curr, up.e[0], v_ctr_new, p2t_ctr, t2p_ctr,
                          s_map, s2pv_map, s_tets);
  if(up.n == 2)
    create_path_int_s_tet_wg(path, e_curr, up.e[1], v_ctr_new, p2t_ctr, t2p_ctr,
                          s_map, s2pv_map, s_tets);
}

void vd_edisc::create_path_p_wg(int path, int tri_id, apf::MeshEntity* e_curr,
                    std::vector<apf::MeshEntity*> & v_sp_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_sp,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_sp,
            std::map<apf::MeshEntity*, int> &s_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &p2sv_map,
                    std::vector<std::vector<apf::MeshEntity*> > &p_tets_s1,
                    std::vector<std::vector<apf::MeshEntity*> > &p_tets_s2) {

  apf::Up up;
  e_lens.m->getUp(e_curr, up);
  p_tets_s1.at(path).at(tri_id) = create_path_int_p_tet_wg(path, e_curr, up.e[0], 
                                              v_sp_new, p2t_sp, t2p_sp,
                          s_map, p2sv_map);
  if(up.n == 2)
    p_tets_s2.at(path).at(tri_id) = create_path_int_p_tet_wg(path, e_curr, 
                      up.e[1], v_sp_new, p2t_sp, t2p_sp, s_map, p2sv_map);
}

// Recreate the path entities for all discs on all slices connected to the paths.
void vd_edisc::copy_path_disc_wg(int path, 
                    std::vector<apf::MeshEntity*> & v_ctr_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_ctr,
                    std::vector<apf::MeshEntity*> & v_sp_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_sp,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_sp,
            std::map<apf::MeshEntity*, int> &s_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &p2sv_map,
                    std::vector<std::vector<apf::MeshEntity*> > &p_tets_s1,
                    std::vector<std::vector<apf::MeshEntity*> > &p_tets_s2,
                    std::vector<std::vector<apf::MeshEntity*> > &s_tets ) {


  //1cell
  create_trial_vert_wg(e_lens.discs.at(path).edge, p2t_sp.at(path), t2p_sp.at(path));
  copy_trial_wg(v_sp_new.at(path), e_lens.discs.at(path).edge, 
                                    p2t_sp.at(path),t2p_sp.at(path));
  copy_trial_wg(v_sp_new.at(path), e_lens.discs.at(path).tri, 
                                    p2t_sp.at(path),t2p_sp.at(path));

  int s1 = e_lens.slices.at(2*path) - 1;
  int s2 = e_lens.slices.at(2*path+1) - 1;

  // The 1cell vertex and edge and split edges and triangles for the first slice:
  create_path_int_s_wg(path, s1, v_ctr_new, p2t_ctr, t2p_ctr, s2pv_map);
  // The 1cell vertex and edge and split edges and triangles for the second slice:
  create_path_int_s_wg(path, s2, v_ctr_new, p2t_ctr, t2p_ctr, s2pv_map);

  // The 1cell vertex and edge and split edges and triangles for the first slice:
  create_path_int_p_wg(path, s1, v_sp_new, p2t_sp, t2p_sp, p2sv_map);
  // The 1cell vertex and edge and split edges and triangles for the second slice:
  create_path_int_p_wg(path, s2, v_sp_new, p2t_sp, t2p_sp, p2sv_map);

  // Split disc triangles and internal tets, from triangles only for upper 
  // adjacency tets belonging to the slice:
  for(int i = 0; i < e_lens.discs.at(path).tri.size(); i++) {
    apf::MeshEntity* e_curr = e_lens.discs.at(path).tri.at(i);
    create_path_s_wg(path, e_curr, v_ctr_new, p2t_ctr, t2p_ctr, s_map, s2pv_map, s_tets);
    create_path_p_wg(path, i, e_curr, v_sp_new, p2t_sp, t2p_sp, s_map, p2sv_map,
                                                  p_tets_s1, p_tets_s2);
  }
}

// Update the positions of the vertex copies on the spheres until non-radial  
// forces vanish.
void vd_edisc::evolve_equi_wg(
                    std::vector<apf::MeshEntity*> & v_ctr_new, 
                    std::vector<apf::MeshEntity*> & v_sp_new, 
                    apf::MeshEntity* v2c_new,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &p2sv_map,
                    std::vector<apf::MeshEntity*> &s2cv_map,
                    std::map<int, apf::MeshEntity*> &c2sv_map,
                    std::vector<apf::MeshEntity*> &p2cv_map,
                    std::map<int, apf::MeshEntity*> &c2pv_map,
                    std::vector<apf::Vector3> &vel_ctr,
                    std::vector<apf::Vector3> &force_ctr,
                    std::vector<apf::Vector3> &vel_sp,
                    std::vector<apf::Vector3> &force_sp,
                    apf::Vector3 &vel_2c,
                    apf::Vector3 &force_2c) {
  apf::Field* vel_field = m_trial->findField("velocity_field");

  // Ratio of the maximum tangential motion length to the sphere radius. 
  double rad_rat = 10.;
  int iter_sz = 50; 
  double rat_th = 10e-3;
  apf::Vector3 v_r(0,0,0);
  apf::Vector3 v_t(0,0,0);
  apf::Vector3 pos_cp(0,0,0);

  // The maximum tangential magnitude of the velocity at the start.
  double mag_t_st = 0;
  // The maximum tangential magnitude of the velocity. It will correspond to 
  // l_min/rad_rat*|v_r|/|v_t| of displacement.
  double mag_t_max = 0;
  double r2t_rat = 1;

  std::vector<std::map<int, apf::Vector3> > s2pv_v_r(s2pv_map.size(), std::map<int, apf::Vector3>{});
  std::vector<std::map<int, apf::Vector3> > p2sv_v_r(p2sv_map.size(), std::map<int, apf::Vector3>{});
  std::map<int, apf::Vector3> s2cv_v_r{};
  std::map<int, apf::Vector3> c2sv_v_r{};
  std::map<int, apf::Vector3> p2cv_v_r{};
  std::map<int, apf::Vector3> c2pv_v_r{};

  std::vector<std::map<int, apf::Vector3> > s2pv_v_t(s2pv_map.size(), std::map<int, apf::Vector3>{});
  std::vector<std::map<int, apf::Vector3> > p2sv_v_t(p2sv_map.size(), std::map<int, apf::Vector3>{});
  std::map<int, apf::Vector3> s2cv_v_t{};
  std::map<int, apf::Vector3> c2sv_v_t{};
  std::map<int, apf::Vector3> p2cv_v_t{};
  std::map<int, apf::Vector3> c2pv_v_t{};

  for(int i = 0; i < v_ctr_new.size(); i++) {
    f_calc.vd_upd_vel_field(m_trial, v_ctr_new.at(i), false);
    apf::getVector(vel_field, v_ctr_new.at(i), 0, vel_ctr.at(i));
    force_ctr.at(i) = f_calc.vd_calc_force(m_trial, v_ctr_new.at(i), false);
  }

  for(int i = 0; i < v_sp_new.size(); i++) {
    f_calc.vd_upd_vel_field(m_trial, v_sp_new.at(i), false);
    apf::getVector(vel_field, v_sp_new.at(i), 0, vel_sp.at(i));

    int s1 = e_lens.slices.at(2*i) - 1;
    int s2 = e_lens.slices.at(2*i+1) - 1;
    vel_sp.at(i) = vel_sp.at(i) + vel_ctr.at(s1)/2 + vel_ctr.at(s2)/2;
    force_sp.at(i) = f_calc.vd_calc_force(m_trial, v_sp_new.at(i), false);
  }

//v2c positioning can be done based on relative velocities of 1-cell vertices. This way the 2-cell will open up in the case of digon insertions, as well.
  // Update positions on the 1cell vertices, check for inversion
  if(trial_type == 2) {
    apf::Vector3 v_mean(0,0,0);
    for(int i = 0; i < vel_ctr.size(); i++)
      v_mean = v_mean + vel_ctr.at(i);
    v_mean = v_mean/vel_ctr.size();

    f_calc.vd_upd_vel_field(m_trial, v2c_new, false);
    apf::getVector(vel_field, v2c_new, 0, vel_2c);
    vel_2c = vel_2c + v_mean;
    force_2c = f_calc.vd_calc_force(m_trial, v2c_new, false);
  }

  for(int path = 0; path < v_sp_new.size(); path++) {
    int s1 = e_lens.slices.at(2*path) - 1;
    int s2 = e_lens.slices.at(2*path+1) - 1;

    // The velocity of the 0cell vertex copy on the 1cell vertex sphere.
    v_t = vel_ctr.at(s1) - vel_sp.at(path);

    apf::MeshEntity* v_cp = p2sv_map.at(path)[s1];
    m_trial->getPoint(v_sp_new.at(path), 0, pos_old);
    m_trial->getPoint(v_cp, 0, pos_cp);

    // 0cell vertex copies:
    pos_cp = pos_cp - pos_old;
    pos_cp = norm_0(pos_cp);
    v_r = pos_cp*(v_t*pos_cp);
    v_t = v_t - v_r;
    p2sv_v_r.at(path)[s1] = v_r;
    p2sv_v_t.at(path)[s1] = v_t;
    s2pv_v_r.at(s1)[path] = v_r*(-1);
    s2pv_v_t.at(s1)[path] = v_t*(-1);
    double mag_curr = v_t.getLength();
    if(mag_curr > mag_t_max) {
      mag_t_max = mag_curr;
      if(v_r.getLength() > v_t.getLength()) {
        r2t_rat = v_t.getLength()/v_r.getLength();
      }
      else
        r2t_rat = 1;
    }

    v_t = vel_ctr.at(s2) - vel_sp.at(path);
    v_cp = p2sv_map.at(path)[s2];
    m_trial->getPoint(v_cp, 0, pos_cp);

    pos_cp = pos_cp - pos_old;
    pos_cp = norm_0(pos_cp);
    v_r = pos_cp*(v_t*pos_cp);
    v_t = v_t - v_r;
    p2sv_v_r.at(path)[s2] = v_r;
    p2sv_v_t.at(path)[s2] = v_t;
    s2pv_v_r.at(s2)[path] = v_r*(-1);
    s2pv_v_t.at(s2)[path] = v_t*(-1);

    mag_curr = v_t.getLength();
    if(mag_curr > mag_t_max) {
      mag_t_max = mag_curr;
      if(v_r.getLength() > v_t.getLength()) {
        r2t_rat = v_t.getLength()/v_r.getLength();
      }
      else
        r2t_rat = 1;
    }

    // 2cell vertex copies:
    v_t = vel_2c - vel_sp.at(path);
    v_cp = p2cv_map[path];
    m_trial->getPoint(v_cp, 0, pos_cp);
    pos_cp = pos_cp - pos_old;
    pos_cp = norm_0(pos_cp);
    v_r = pos_cp*(v_t*pos_cp);
    v_t = v_t - v_r;
    p2cv_v_r[path] = v_r;
    p2cv_v_t[path] = v_t;
    s2cv_v_r[path] = v_r*(-1);
    s2cv_v_t[path] = v_t*(-1);

    mag_curr = v_t.getLength();
    if(mag_curr > mag_t_max) {
      mag_t_max = mag_curr;
      if(v_r.getLength() > v_t.getLength()) {
        r2t_rat = v_t.getLength()/v_r.getLength();
      }
      else
        r2t_rat = 1;
    }

  }

  if(trial_type == 2) {
    for(int i = 0; i < v_ctr_new.size(); i++) {
      apf::MeshEntity* v_cp = c2sv_map[i];
      m_trial->getPoint(v_cp, 0, pos_cp);

      // The velocity of the 0cell vertex copy on the 2cell vertex sphere.
      v_t = vel_ctr.at(i) - vel_2c;
      pos_cp = pos_cp - pos_old;
      pos_cp = norm_0(pos_cp);
      v_r = pos_cp*(v_t*pos_cp);
      v_t = v_t - v_r;

      c2sv_v_r[i] = v_r;
      c2sv_v_t[i] = v_t;
      s2cv_v_r[i] = v_r*(-1);
      s2cv_v_t[i] = v_t*(-1);
      double mag_curr = v_t.getLength();
      if(mag_curr > mag_t_max) {
        mag_t_max = mag_curr;
        if(v_r.getLength() > v_t.getLength()) {
          r2t_rat = v_t.getLength()/v_r.getLength();
        }
        else
          r2t_rat = 1;
      }
    }
  }

  // Repositioning:
  // 0- and 1-cell vertices:
  double scale = 1/mag_t_max/rad_rat*r2t_rat;
  for(int path = 0; path < v_sp_new.size(); path++) {
    int s1 = e_lens.slices.at(2*path) - 1;
    int s2 = e_lens.slices.at(2*path+1) - 1;

    // s1
    apf::MeshEntity* v_cp = p2sv_map.at(path)[s1];
    m_trial->getPoint(v_cp, 0, pos_cp);

    double mag_curr = p2sv_v_t.at(path)[s1].getLength();
    apf::Vector3 temp(0,0,0);
    temp = norm_0(p2sv_v_t.at(path)[s1])*mag_curr*scale;
    m_trial->setPoint(v_cp, 0, pos_cp + temp);
    vd_proj_v_sphere(m_trial, v_cp, pos_old, l_min);

    v_cp = s2pv_map.at(s1)[path];
    m_trial->getPoint(v_cp, 0, pos_cp);
    m_trial->setPoint(v_cp, 0, pos_cp - temp);
    vd_proj_v_sphere(m_trial, v_cp, pos_old, l_min);

    // s2
    v_cp = p2sv_map.at(path)[s2];
    m_trial->getPoint(v_cp, 0, pos_cp);

    mag_curr = p2sv_v_t.at(path)[s2].getLength();
    temp = norm_0(p2sv_v_t.at(path)[s1])*mag_curr*scale;
    m_trial->setPoint(v_cp, 0, pos_cp + temp);
    vd_proj_v_sphere(m_trial, v_cp, pos_old, l_min);

    v_cp = s2pv_map.at(s2)[path];
    m_trial->getPoint(v_cp, 0, pos_cp);
    m_trial->setPoint(v_cp, 0, pos_cp - temp);
    vd_proj_v_sphere(m_trial, v_cp, pos_old, l_min);
  }
  // 2-cell vertices:

  if(trial_type == 2) {
    for(int i = 0; i < v_ctr_new.size(); i++) {
      apf::MeshEntity* v_cp = c2sv_map[i];
      m_trial->getPoint(v_cp, 0, pos_cp);

      double mag_curr = c2sv_v_r[i].getLength();
      apf::Vector3 temp(0,0,0);
      temp = norm_0(c2sv_v_r[i])*mag_curr*scale;
      m_trial->setPoint(v_cp, 0, pos_cp + temp);
      vd_proj_v_sphere(m_trial, v_cp, pos_old, l_min);

      v_cp = s2cv_map[i];
      m_trial->getPoint(v_cp, 0, pos_cp);
      m_trial->setPoint(v_cp, 0, pos_cp - temp);
      vd_proj_v_sphere(m_trial, v_cp, pos_old, l_min);
    }
    for(int i = 0; i < v_sp_new.size(); i++) {
      apf::MeshEntity* v_cp = c2pv_map[i];
      m_trial->getPoint(v_cp, 0, pos_cp);

      double mag_curr = c2pv_v_r[i].getLength();
      apf::Vector3 temp(0,0,0);
      temp = norm_0(c2pv_v_r[i])*mag_curr*scale;
      m_trial->setPoint(v_cp, 0, pos_cp + temp);
      vd_proj_v_sphere(m_trial, v_cp, pos_old, l_min);

      v_cp = p2cv_map[i];
      m_trial->getPoint(v_cp, 0, pos_cp);
      m_trial->setPoint(v_cp, 0, pos_cp - temp);
      vd_proj_v_sphere(m_trial, v_cp, pos_old, l_min);
    }
  }
/*
At the start make sure all tets(new tets around 1cell vertices?) are positive volume
Find the relative velocities. Find the tangential component. Find the maximum of the tangential velocity magnitudes. Scale the velocities such that the maximum tangential velocity vertex will move by the (minimum of 1 or |vel_t|/|vel_r|for the maximum tangential velocity vertex)*radius/rad_rat
Iterate while |vel_t|/|vel_r| for maximum vel_t vertex is above th or iter < iter_sz

If the tets invert at any point, discard insertion
*/
}
/*

// Update the positions of the vertices on the sphere until non-radial forces 
// vanish.
void vd_edisc::evolve_equi_wg(apf::MeshEntity* v_ctr) {
  double mult_rat = 5.;
  double end_rat = 10e-4;
  std::map<apf::MeshEntity*, bool> tet_skip{};

  std::vector<apf::MeshEntity*> edges(0);
  std::vector<apf::MeshEntity*> v_other(0);
  std::vector<apf::Vector3> f_other(0);
  vd_set_up(m_trial, v_ctr, &edges);
  v_other.resize(edges.size());
  f_other.resize(edges.size());
  for(int i = 0; i < edges.size(); i++) {
    v_other.at(i) = getEdgeVertOppositeVert(m_trial, edges.at(i), v_ctr);
  }
  double mag_st = 0;
  for(int i = 0; i < v_other.size(); i++) {
    f_other.at(i) = f_calc.vd_calc_force(m_trial, v_other.at(i), false);
    double mag_curr = f_other.at(i).getLength();
    if(mag_curr > mag_st)
      mag_st = mag_curr;
  }
  double mult_st = vd_find_min_mult(m_trial, &v_other, tet_skip, f_other)
                                                              /mult_rat;
  if(mag_st < std::numeric_limits<double>::min())
    mag_st = 1;
  if(mult_st < std::numeric_limits<double>::min())
    mult_st = 1;

  apf::Vector3 pos_rel(0,0,0);
  double rat_max = 1;
  double mult = mult_st;
  while(rat_max > end_rat and mult > mult_st*end_rat) {
    double mag_max = 0;
    for(int i = 0; i < v_other.size(); i++) {
      f_other.at(i) = f_calc.vd_calc_force(m_trial, v_other.at(i), false);
      m_trial->getPoint(v_other.at(i), 0, pos_rel);
      pos_rel = norm_0(pos_rel - pos_old);
      f_other.at(i) = f_other.at(i) - pos_rel*(f_other.at(i)*pos_rel);
      double mag_curr = f_other.at(i).getLength();
      if(mag_curr > mag_max)
        mag_max = mag_curr;
    }
    rat_max = mag_max/mag_st;
    if(rat_max > end_rat) {
      mult = vd_find_min_mult(m_trial, &v_other, tet_skip, f_other)
                                                                  /mult_rat;
      for(int i = 0; i < v_other.size(); i++) {
        m_trial->getPoint(v_other.at(i), 0, pos_rel);
        pos_rel = pos_rel + f_other.at(i)*mult;
        m_trial->setPoint(v_other.at(i), 0, pos_rel);
        vd_proj_v_sphere(m_trial, v_other.at(i), pos_old, l_min);
      }
    }
  }
}
*/
// Update the positions on the spheres by updating the forces and velocities.
void vd_edisc::upd_pos_sph_wg(std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells,
                    std::vector<apf::MeshEntity*> & v_ctr_new, 
                    std::vector<apf::MeshEntity*> & v_sp_new, 
                    apf::MeshEntity* v2c_new,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &p2sv_map,
                    std::vector<apf::MeshEntity*> &s2cv_map,
                    std::map<int, apf::MeshEntity*> &c2sv_map,
                    std::vector<apf::MeshEntity*> &p2cv_map,
                    std::map<int, apf::MeshEntity*> &c2pv_map,
                    std::vector<apf::Vector3> &vel_ctr,
                    std::vector<apf::Vector3> &force_ctr,
                    std::vector<apf::Vector3> &vel_sp,
                    std::vector<apf::Vector3> &force_sp,
                    apf::Vector3 &vel_2c,
                    apf::Vector3 &force_2c,
                    std::vector<std::vector<apf::MeshEntity*> >& p_tets_s1,
                    std::vector<std::vector<apf::MeshEntity*> >& p_tets_s2,
                    std::vector<std::vector<apf::MeshEntity*> > &s_tets) {

  apf::Field* vel_field = m_trial->findField("velocity_field");
  assert(l_min > std::numeric_limits<double>::min());
  
  for(int i = 0; i < v_ctr_new.size(); i++) {
    f_calc.vd_upd_vel_field(m_trial, v_ctr_new.at(i), false);
    apf::getVector(vel_field, v_ctr_new.at(i), 0, vel_ctr.at(i));
    force_ctr.at(i) = f_calc.vd_calc_force(m_trial, v_ctr_new.at(i), false);
  }
  if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
    ext_shell* e_sh = f_calc.get_e_sh();
    for(int i = 0; i < v_ctr_new.size(); i++) {
      if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(i) - 1) ) {
        vel_ctr.at(i) = 
          e_sh->find_para_dir(0, e_lens.new_cell0_id.at(i) - 1, vel_ctr.at(i));
        force_ctr.at(i) = 
          e_sh->find_para_dir(0, e_lens.new_cell0_id.at(i) - 1, force_ctr.at(i));
      }
    }
  }
  for(int path = 0; path < v_sp_new.size(); path++) {
    int s1 = e_lens.slices.at(2*path) - 1;
    int s2 = e_lens.slices.at(2*path+1) - 1;

    apf::MeshEntity* v_ctr_cp = p2sv_map.at(path)[s1];
    apf::Vector3 pos_new(0,0,0);
    pos_new = norm_0(vel_ctr.at(s1))*(l_min) + pos_old;
    m_trial->setPoint(v_ctr_cp, 0, pos_new);

    v_ctr_cp = p2sv_map.at(path)[s2];
    pos_new = norm_0(vel_ctr.at(s2))*(l_min) + pos_old;
    m_trial->setPoint(v_ctr_cp, 0, pos_new);
  }

  for(int i = 0; i < v_sp_new.size(); i++) {
    f_calc.vd_upd_vel_field(m_trial, v_sp_new.at(i), false);
    apf::getVector(vel_field, v_sp_new.at(i), 0, vel_sp.at(i));
    force_sp.at(i) = f_calc.vd_calc_force(m_trial, v_sp_new.at(i), false);
  }
  if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
    ext_shell* e_sh = f_calc.get_e_sh();
    for(int i = 0; i < v_sp_new.size(); i++) {
      if(e_sh->chk_shell(1, e_lens.new_cell1_id.at(i) - 1) ) {
        vel_sp.at(i) = 
          e_sh->find_para_dir(1, e_lens.new_cell1_id.at(i) - 1, vel_sp.at(i));
        force_sp.at(i) = 
          e_sh->find_para_dir(1, e_lens.new_cell1_id.at(i) - 1, force_sp.at(i));
      }
    }
  }
  apf::Vector3 dir_rel(0,0,0);
  if(slice_cells->at(0).first.first == 0) {
    apf::MeshEntity* v_sp_cp = s2pv_map.at(0)[0];
    dir_rel = norm_0(vel_ctr.at(1) - vel_ctr.at(0))*(l_min);
    m_trial->setPoint(v_sp_cp, 0, pos_old + dir_rel);
    v_sp_cp = s2pv_map.at(1)[0];
    m_trial->setPoint(v_sp_cp, 0, pos_old - dir_rel);
  }
  else {
    for(int slice = 0; slice < v_ctr_new.size(); slice++) {
      int d1 = slice_cells->at(slice).first.first - 1;
      int d2 = slice_cells->at(slice).first.second - 1;

      int s1 = e_lens.slices.at(2*d1) - 1;
      int s2 = e_lens.slices.at(2*d1+1) - 1;
      if(s1 != slice) {
        s2 = s1;
        s1 = slice;
      }
      apf::MeshEntity* v_sp_cp = s2pv_map.at(slice)[d1];
      dir_rel = norm_0(vel_ctr.at(s2) - vel_ctr.at(s1))*(l_min);
      m_trial->setPoint(v_sp_cp, 0, pos_old + dir_rel);

      s1 = e_lens.slices.at(2*d2) - 1;
      s2 = e_lens.slices.at(2*d2+1) - 1;
      if(s1 != slice) {
        s2 = s1;
        s1 = slice;
      }
      v_sp_cp = s2pv_map.at(slice)[d2];
      dir_rel = norm_0(vel_ctr.at(s2) - vel_ctr.at(s1))*(l_min);
      m_trial->setPoint(v_sp_cp, 0, pos_old + dir_rel);
    }
  }
  for(int i = 0; i < vel_sp.size(); i++) {
    for(int j = 0; j < p_tets_s1.at(i).size(); j++) {
      if(vd_volume_tet(m_trial, p_tets_s1.at(i).at(j)) < 
                            std::numeric_limits<double>::min()) {
        ins_flag = false;
        return;
      }
    }
    for(int j = 0; j < p_tets_s2.size(); j++) {
      if(vd_volume_tet(m_trial, p_tets_s2.at(i).at(j)) < 
                            std::numeric_limits<double>::min()) {
        ins_flag = false;
        return;
      }
    }
  }

//v2c positioning can be done based on relative velocities of 1-cell vertices. This way the 2-cell will open up in the case of digon insertions, as well.
  // Update positions on the 1cell vertices, check for inversion
  if(trial_type == 2) {
    apf::Vector3 v_mean(0,0,0);
    for(int i = 0; i < vel_ctr.size(); i++)
      v_mean = v_mean + vel_ctr.at(i);
    v_mean = v_mean/vel_ctr.size();

    for(int i = 0; i < vel_ctr.size(); i++) {
      apf::MeshEntity* v_ctr_cp = c2sv_map[i];
      dir_rel = norm_0(vel_ctr.at(i) - v_mean)*(l_min);
      m_trial->setPoint(v_ctr_cp, 0, pos_old + dir_rel);
      apf::MeshEntity* v_2c_cp = s2cv_map.at(i);
      m_trial->setPoint(v_2c_cp, 0, pos_old - dir_rel);
    }
    for(int i = 0; i < vel_sp.size(); i++) {
      apf::MeshEntity* v_sp_cp = c2pv_map[i];
      dir_rel = norm_0(vel_sp.at(i) - v_mean)*(l_min);
      m_trial->setPoint(v_sp_cp, 0, pos_old + dir_rel);
      apf::MeshEntity* v_2c_cp = p2cv_map.at(i);
      m_trial->setPoint(v_2c_cp, 0, pos_old - dir_rel);
    }

    f_calc.vd_upd_vel_field(m_trial, v2c_new, false);
    apf::getVector(vel_field, v2c_new, 0, vel_2c);
    force_2c = f_calc.vd_calc_force(m_trial, v2c_new, false);

    if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
      ext_shell* e_sh = f_calc.get_e_sh();
      if(e_sh->chk_shell(2, e_lens.new_cell2_id - 1) ) {
        vel_2c = e_sh->find_para_dir(2, e_lens.new_cell2_id - 1, vel_2c);
        force_2c = e_sh->find_para_dir(2, e_lens.new_cell2_id - 1, force_2c);
      }
    }
  }

  if(sub_vtk) {

    for(int slice = 0; slice < v_ctr_new.size(); slice++) {
      std::string vtk_name = "./output/sub/v_trial_exp_" + std::to_string(trial_type) + "_" + std::to_string(trial_curr); 
      if(trial_type == 2)
        vtk_name = vtk_name + "_" + std::to_string(ng.curr);
      vtk_name = vtk_name + "_s" + std::to_string(slice);
      vd_save_vtk_vert(m_trial, v_ctr_new.at(slice), vtk_name.c_str());
      vtk_name = "./output/v_trial_exp"+ std::to_string(slice);
      vd_save_vtk_vert(m_trial, v_ctr_new.at(slice), vtk_name.c_str());
    }

    for (int path = 0; path < v_sp_new.size(); path++) {
      std::string vtk_name = "./output/sub/v_trial_exp_" + std::to_string(trial_type) + "_" + std::to_string(trial_curr); 
      if(trial_type == 2)
        vtk_name = vtk_name + "_" + std::to_string(ng.curr);
      vtk_name = vtk_name + "_p" + std::to_string(path);
      vd_save_vtk_vert(m_trial, v_sp_new.at(path), vtk_name.c_str());
    }
  }
/*
  // Evolve the surrounding vertices.
  for(int i = 0; i < v_ctr_new.size(); i++) {
    evolve_equi_wg(v_ctr_new.at(i));
  }
  for(int i = 0; i < v_sp_new.size(); i++) {
    evolve_equi_wg(v_sp_new.at(i));
  }

  for(int slice = 0; slice < v_ctr_new.size(); slice++) {
    std::string vtk_name = "./output/sub/v_trial_evol_" + std::to_string(trial_type) + "_" + std::to_string(trial_curr); 
    if(trial_type == 2)
      vtk_name = vtk_name + "_" + std::to_string(ng.curr);
    vtk_name = vtk_name + "_s" + std::to_string(slice);
    vd_save_vtk_vert(m_trial, v_ctr_new.at(slice), vtk_name.c_str());
  }

  for (int path = 0; path < v_sp_new.size(); path++) {
    std::string vtk_name = "./output/sub/v_trial_evol_" + std::to_string(trial_type) + "_" + std::to_string(trial_curr); 
    if(trial_type == 2)
      vtk_name = vtk_name + "_" + std::to_string(ng.curr);
    vtk_name = vtk_name + "_p" + std::to_string(path);
    vd_save_vtk_vert(m_trial, v_sp_new.at(path), vtk_name.c_str());
  }
*/
  for(int i = 0; i < v_ctr_new.size(); i++) {
    f_calc.vd_upd_vel_field(m_trial, v_ctr_new.at(i), false);
    apf::getVector(vel_field, v_ctr_new.at(i), 0, vel_ctr.at(i));
    force_ctr.at(i) = f_calc.vd_calc_force(m_trial, v_ctr_new.at(i), false);
  }
  for(int i = 0; i < v_sp_new.size(); i++) {
    f_calc.vd_upd_vel_field(m_trial, v_sp_new.at(i), false);
    apf::getVector(vel_field, v_sp_new.at(i), 0, vel_sp.at(i));
    force_sp.at(i) = f_calc.vd_calc_force(m_trial, v_sp_new.at(i), false);
  }
  if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
    ext_shell* e_sh = f_calc.get_e_sh();
    for(int i = 0; i < v_ctr_new.size(); i++) {
      if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(i) - 1) ) {
        vel_ctr.at(i) = 
          e_sh->find_para_dir(0, e_lens.new_cell0_id.at(i) - 1, vel_ctr.at(i));
        force_ctr.at(i) = 
          e_sh->find_para_dir(0, e_lens.new_cell0_id.at(i) - 1, force_ctr.at(i));
      }
    }
    for(int i = 0; i < v_sp_new.size(); i++) {
      if(e_sh->chk_shell(1, e_lens.new_cell1_id.at(i) - 1) ) {
        vel_sp.at(i) = 
          e_sh->find_para_dir(1, e_lens.new_cell1_id.at(i) - 1, vel_sp.at(i));
        force_sp.at(i) = 
          e_sh->find_para_dir(1, e_lens.new_cell1_id.at(i) - 1, force_sp.at(i));
      }
    }
  }
/*
Actually consider the local environments as additional degrees of freedom.
After initial expansion allow the vertices on the sphere to relax.They will reach equilibrium configuration (flat boundries for 2cell vertices, dihedral angles for 1cell vertices. Calculate rate of change of tets afterwards.
But an important question... is the equilibrium also reached for the vertex center of the sphere. The centers of the spheres will act like suddenly unloaded springs. If the resultant motion is closing the newly generated tets, than it is closing.
The tets to calculate the rates of change are not the 1-cell sphere ones, but the copies of the tets on the 0-cells for 1-stratum insertion. For 2-stratum insertion, the copies of the void ones on 1-cell spheres
*/
  for(int i = 0; i < vel_ctr.size(); i++) {
    for(int j = 0; j < s_tets.at(i).size(); j++) {
      //if(vd_volume_tet(m_trial, p_tets_s1.at(i).at(j)) < 
      if(vd_volume_tet(m_trial, s_tets.at(i).at(j)) < 
                            std::numeric_limits<double>::min()
          or
          vd_calc_roc(m_trial, s_tets.at(i).at(j), NULL) < 
                            std::numeric_limits<double>::min()) {
        ins_flag = false;
        return;
      }
    }
  }
/*
  for(int i = 0; i < vel_sp.size(); i++) {
    for(int j = 0; j < p_tets_s1.at(i).size(); j++) {
      //if(vd_volume_tet(m_trial, p_tets_s1.at(i).at(j)) < 
      if(f_calc.calc_roc(m_trial, p_tets_s1.at(i).at(j)) < 
                            std::numeric_limits<double>::min()) {
        ins_flag = false;
        return;
      }
    }
    for(int j = 0; j < p_tets_s2.size(); j++) {
      //if(vd_volume_tet(m_trial, p_tets_s2.at(i).at(j)) < 
      if(f_calc.calc_roc(m_trial, p_tets_s2.at(i).at(j)) < 
                            std::numeric_limits<double>::min()) {
        ins_flag = false;
        return;
      }
    }
  }
*/
}

  // Expand the 2-cell void:
void vd_edisc::expand_2c_void_wg(
      std::vector<std::pair<std::vector<int >, std::vector<int > > >* path_cells,
                    std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells,
                    std::vector<apf::MeshEntity*> & v_ctr_new,
                    std::vector<apf::MeshEntity*> & v_sp_new, 
                    apf::MeshEntity* v2c_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_sp,
            std::map<apf::MeshEntity*,apf::MeshEntity* > &p2c_2c,
                    std::vector<std::map<int, apf::MeshEntity*> > &p2sv_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map,
                    std::vector<apf::MeshEntity*> &s2cv_map,
                    std::map<int, apf::MeshEntity*> &c2sv_map,
                    std::vector<apf::MeshEntity*> &p2cv_map,
                    std::map<int, apf::MeshEntity*> &c2pv_map,
                    std::map<apf::MeshEntity*, int> &s_map,
                    std::vector<std::vector<apf::MeshEntity*> > &s_tets) {
  if(trial_type == 2) {
    apf::Downward down;
    int tag = e_lens.new_cell2_id;
    apf::ModelEntity* em_2c = m_trial->findModelEntity(2,tag);
    // 2-cell edges between the 2-cell vertex and 0- and 1-cell vertices on
    // all spheres:
    for(int slice = 0; slice < slice_cells->size(); slice++) {
      // 0-cell copy:
      down[0] = v_ctr_new.at(slice);
      apf::MeshEntity* vc2_cp = m_trial->createVert(em_2c);
      down[1] = vc2_cp;
      s2cv_map.at(slice) = down[1];
      m_trial->setPoint(vc2_cp, 0, pos_old);
      buildElement(m_trial, em_2c, apf::Mesh::EDGE, down, 0);

      // 2-cell copy:
      down[0] = v2c_new;
      apf::MeshEntity* v_ctr_cp = m_trial->createVert(em_2c);
      down[1] = v_ctr_cp;
      c2sv_map[slice] = down[1];
      m_trial->setPoint(v_ctr_cp, 0, pos_old);
      buildElement(m_trial, em_2c, apf::Mesh::EDGE, down, 0);
    }
    for(int path = 0; path < e_lens.discs.size(); path++) {
      // 1-cell copy:
      down[0] = v_sp_new.at(path);
      apf::MeshEntity* vc2_cp = m_trial->createVert(em_2c);
      down[1] = vc2_cp;
      p2cv_map.at(path) = down[1];
      m_trial->setPoint(vc2_cp, 0, pos_old);
      buildElement(m_trial, em_2c, apf::Mesh::EDGE, down, 0);

      // 2-cell copy:
      down[0] = v2c_new;
      apf::MeshEntity* v_sp_cp = m_trial->createVert(em_2c);
      down[1] = v_sp_cp;
      c2pv_map[path] = down[1];
      m_trial->setPoint(v_sp_cp, 0, pos_old);
      buildElement(m_trial, em_2c, apf::Mesh::EDGE, down, 0);
    }
    // 2-cell edges from the 2-cell vertex copy to the 1-cell copies on the 0-
    // cell spheres:
    for(int slice = 0; slice < slice_cells->size(); slice++) {
      int d1 = slice_cells->at(slice).first.first - 1;
      int d2 = slice_cells->at(slice).first.second - 1;

      down[0] = s2cv_map.at(slice);
      down[1] = s2pv_map.at(slice)[d1];

      buildElement(m_trial, em_2c, apf::Mesh::EDGE, down, 0);
      down[1] = s2pv_map.at(slice)[d2];
      buildElement(m_trial, em_2c, apf::Mesh::EDGE, down, 0);
    }

    // 2-cell edges from the 2-cell vertex copy to the 0-cell copies on the 1-
    // cell spheres:
    for(int path = 0; path < e_lens.discs.size(); path++) {
      int s1 = e_lens.slices.at(2*path) - 1;
      int s2 = e_lens.slices.at(2*path+1) - 1;

      down[0] = p2cv_map.at(path);
      down[1] = p2sv_map.at(path)[s1];

      buildElement(m_trial, em_2c, apf::Mesh::EDGE, down, 0);
      down[1] = p2sv_map.at(path)[s2];
      buildElement(m_trial, em_2c, apf::Mesh::EDGE, down, 0);
    }

    // 2-cell edges from the 0-cell copies to the 1-cell copies on the 2-
    // cell spheres:
    for(int path = 0; path < e_lens.discs.size(); path++) {
      int s1 = e_lens.slices.at(2*path) - 1;
      int s2 = e_lens.slices.at(2*path+1) - 1;

      down[0] = c2pv_map[path];
      down[1] = c2sv_map[s1];
      buildElement(m_trial, em_2c, apf::Mesh::EDGE, down, 0);
      down[1] = c2sv_map[s2];
      buildElement(m_trial, em_2c, apf::Mesh::EDGE, down, 0);
    }
    /////
    int c3_top = path_cells->at(0).first.at(0);
    int c3_bot = path_cells->at(0).first.back();
    apf::MeshEntity* c3_t_e = NULL;
    apf::MeshEntity* c3_b_e = NULL;
    apf::MeshEntity* c3_t_v = NULL;
    apf::MeshEntity* c3_b_v = NULL;
    apf::ModelEntity* mdl_c3_t = NULL;
    apf::ModelEntity* mdl_c3_b = NULL;
    apf::Vector3 pos(0,0,0);
    // Find the 3-cell edges on paths, their copies exist on 0-cell and 1-cell 
    // spheres. Assert existence
    // Copy both 3-cell edge to the 2-cell sphere.
    if(c3_top != -1) {
      mdl_c3_t = vd_cd->get_mdl(3, c3_top);
      int tag = e_lens.m->getModelTag(mdl_c3_t);
      mdl_c3_t = m_trial->findModelEntity(3, tag);

      c3_t_e = c3_edge.at(c3_top);
      c3_t_v = getEdgeVertOppositeVert(e_lens.m, c3_t_e, e_lens.vert_ctr);
      for(int path = 0; path < e_lens.discs.size(); path++) {
        assert(p2t_sp.at(path)[c3_t_e]);
      }
      for(int slice = 0; slice < slice_cells->size(); slice++) {
        assert(p2t_ctr.at(slice)[c3_t_e]);
      }
      e_lens.m->getPoint(c3_t_v, 0, pos);
      apf::MeshEntity* vc3_cp = m_trial->createVert(mdl_c3_t);
      m_trial->setPoint(vc3_cp, 0, pos);
      p2c_2c[c3_t_v] = vc3_cp;

      down[0] = v2c_new;
      down[1] = vc3_cp;
      buildElement(m_trial, mdl_c3_t, apf::Mesh::EDGE, down, 0);
    }
    if(c3_bot != -1) {
      mdl_c3_b = vd_cd->get_mdl(3, c3_bot);
      int tag = e_lens.m->getModelTag(mdl_c3_b);
      mdl_c3_b = m_trial->findModelEntity(3, tag);

      c3_b_e = c3_edge.at(c3_bot);
      c3_b_v = getEdgeVertOppositeVert(e_lens.m, c3_b_e, e_lens.vert_ctr);
      for(int path = 0; path < e_lens.discs.size(); path++) {
        assert(p2t_sp.at(path)[c3_b_e]);
      }
      for(int slice = 0; slice < slice_cells->size(); slice++) {
        assert(p2t_ctr.at(slice)[c3_b_e]);
      }
      e_lens.m->getPoint(c3_b_v, 0, pos);
      apf::MeshEntity* vc3_cp = m_trial->createVert(mdl_c3_b);
      m_trial->setPoint(vc3_cp, 0, pos);
      p2c_2c[c3_b_v] = vc3_cp;

      down[0] = v2c_new;
      down[1] = vc3_cp;
      buildElement(m_trial, mdl_c3_b, apf::Mesh::EDGE, down, 0);
    }

    // Form the 2-cell triangles for 0-, 1-cell and 2-cell spheres
    for(int path = 0; path < e_lens.discs.size(); path++) {
      int s1 = e_lens.slices.at(2*path) - 1;
      int s2 = e_lens.slices.at(2*path+1) - 1;
      // s1
      // 0-cell copies:
      down[0] = v_ctr_new.at(s1);
      down[1] = s2pv_map.at(s1)[path];
      down[2] = s2cv_map.at(s1);
      buildElement(m_trial, em_2c, apf::Mesh::TRIANGLE, down, 0);

      // 1-cell copies:
      down[0] = p2sv_map.at(path)[s1];
      down[1] = v_sp_new.at(path);
      down[2] = p2cv_map.at(s1);
      buildElement(m_trial, em_2c, apf::Mesh::TRIANGLE, down, 0);

      // 2-cell copies:
      down[0] = c2sv_map[s1];
      down[1] = c2pv_map[path];
      down[2] = v2c_new;
      buildElement(m_trial, em_2c, apf::Mesh::TRIANGLE, down, 0);
      // s2
      // 0-cell copies:
      down[0] = v_ctr_new.at(s2);
      down[1] = s2pv_map.at(s2)[path];
      down[2] = s2cv_map.at(s2);
      buildElement(m_trial, em_2c, apf::Mesh::TRIANGLE, down, 0);

      // 1-cell copies:
      down[0] = p2sv_map.at(path)[s2];
      down[1] = v_sp_new.at(path);
      down[2] = p2cv_map.at(s2);
      buildElement(m_trial, em_2c, apf::Mesh::TRIANGLE, down, 0);

      // 2-cell copies:
      down[0] = c2sv_map[s2];
      down[1] = c2pv_map[path];
      down[2] = v2c_new;
      buildElement(m_trial, em_2c, apf::Mesh::TRIANGLE, down, 0);
    }

    apf::Up up;
    if(c3_top != -1) {
      e_lens.m->getUp(c3_t_e, up);
      std::vector<apf::MeshEntity*> tri_3c(e_lens.discs.size(), NULL);
      for(int path = 0; path < e_lens.discs.size(); path++) {
        bool found = false;
        for(int i = 0; i < up.n; i++) {
          for(int j = 0; j < e_lens.discs.at(path).tri.size(); j++) {
            if(up.e[i] == e_lens.discs.at(path).tri.at(j)) {
              assert(!found);
              assert(e_lens.m->toModel(up.e[i]) == vd_cd->get_mdl(3, c3_top));
              found = true;
              tri_3c.at(path) = up.e[i];
            }
          }
        }
      }
      // Loop over paths, create 3c tets for all spheres.
      for(int path = 0; path < e_lens.discs.size(); path++) {
        e_lens.m->getUp(tri_3c.at(path), up);

        create_3c_tet_wg(path, tri_3c.at(path), up.e[0], v_ctr_new, v_sp_new, 
                    v2c_new, c3_t_v, p2t_ctr, p2t_sp,
                    p2c_2c, p2sv_map, s2pv_map, s2cv_map, c2sv_map,
                    p2cv_map, c2pv_map, s_map, s_tets);
        assert(up.n == 2);
        create_3c_tet_wg(path, tri_3c.at(path), up.e[1], v_ctr_new, v_sp_new, 
                    v2c_new, c3_t_v, p2t_ctr, p2t_sp,
                    p2c_2c, p2sv_map, s2pv_map, s2cv_map, c2sv_map,
                    p2cv_map, c2pv_map, s_map, s_tets);
      }
    }
    if(c3_bot != -1) {
      e_lens.m->getUp(c3_b_e, up);
      std::vector<apf::MeshEntity*> tri_3c(e_lens.discs.size(), NULL);
      for(int path = 0; path < e_lens.discs.size(); path++) {
        bool found = false;
        for(int i = 0; i < up.n; i++) {
          for(int j = 0; j < e_lens.discs.at(path).tri.size(); j++) {
            if(up.e[i] == e_lens.discs.at(path).tri.at(j)) {
              assert(!found);
              assert(e_lens.m->toModel(up.e[i]) == vd_cd->get_mdl(3, c3_bot));
              found = true;
              tri_3c.at(path) = up.e[i];
            }
          }
        }
      }
      // Loop over paths, create 3c tets for all spheres.
      for(int path = 0; path < e_lens.discs.size(); path++) {
        e_lens.m->getUp(tri_3c.at(path), up);

        create_3c_tet_wg(path, tri_3c.at(path), up.e[0], v_ctr_new, v_sp_new, 
                    v2c_new, c3_b_v, p2t_ctr, p2t_sp,
                    p2c_2c, p2sv_map, s2pv_map, s2cv_map, c2sv_map,
                    p2cv_map, c2pv_map, s_map, s_tets);
        assert(up.n == 2);
        create_3c_tet_wg(path, tri_3c.at(path), up.e[1], v_ctr_new, v_sp_new, 
                    v2c_new, c3_b_v, p2t_ctr, p2t_sp,
                    p2c_2c, p2sv_map, s2pv_map, s2cv_map, c2sv_map,
                    p2cv_map, c2pv_map, s_map, s_tets);
      }
    }
  }
}

// Create a local neighborhood for each vertex to calculate dissipation rate for
// and infinitesimal insertion. 
// TODO Update the local fields for copied entities and the assign correct fields
// for newly generated ones. Doesn't create any issues for constant boundary 
// energy, but would create issues for things like stored deformation energy.
void vd_edisc::reload_trial_wg(
                std::vector<std::pair<std::vector<int >, std::vector<int > > >*
                                                                  path_cells,
                std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells) {
  assert(m_ex and m_load);
  assert(precond_ex and precond_load);

  //std::cout << "Mesh cond: "<< trial_ex << " " << trial_load << std::endl;
  //vd_print_ent(m_precond);

  if(!trial_ex) {
    assert(!trial_load);

    m_trial = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);
    trial_ex = true;
    std::cout << "Created empty mesh, " << trial_load << std::endl;
  }

  if(trial_load) {
    apf::destroyMesh(m_trial);
    m_trial = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);

    trial_load = false;
    std::cout << "Recreate empty mesh, " << trial_load << std::endl;
  }

  std::cout << "Initial:" << m_main << ", Precond:" << m_precond
            << ", Trial:" << m_trial << std::endl;

  // Entity maps from m_precond to m_trial for each new 0- and 1-cell 
  // and vica versa
  int v_nbr_ctr = e_lens.new_cell0_id.size();
  int v_nbr_sp = e_lens.new_cell1_id.size();

  std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > p2t_ctr(v_nbr_ctr, 
                          std::map<apf::MeshEntity*,apf::MeshEntity* > {});
  std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > t2p_ctr(v_nbr_ctr, 
                          std::map<apf::MeshEntity*,apf::MeshEntity* > {});
  std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > p2t_sp(v_nbr_sp, 
                          std::map<apf::MeshEntity*,apf::MeshEntity* > {});
  std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > t2p_sp(v_nbr_sp, 
                          std::map<apf::MeshEntity*,apf::MeshEntity* > {});

  // m_precond to m_trial for 2-cell vertex and vice versa.
  std::map<apf::MeshEntity*,apf::MeshEntity* > p2c_2c {};
  std::map<apf::MeshEntity*,apf::MeshEntity* > c2p_2c {};

  // The slice map of the preconditioned mesh tetrahedra. Used in determining
  // the vertex ordering for the tets surrounding new 1-cell edges.
  std::map<apf::MeshEntity*, int> s_map{};
  // The map from path id to vertex for each slice. Each slice will be at most 
  // associated with two paths. Used to quickly update the positions of the
  // associated copies of each new 0-cell and 1-cell in other local neighborhoods.
  std::vector<std::map<int, apf::MeshEntity*> > s2pv_map(
            e_lens.new_cell0_id.size(), std::map<int, apf::MeshEntity*> {});
  // The map from slice id to vertex for each path. Each path will be associated 
  // with two slices.
  std::vector<std::map<int, apf::MeshEntity*> > p2sv_map(
            e_lens.new_cell1_id.size(), std::map<int, apf::MeshEntity*> {});

  // Maps to 2cell vertex and vice versa
  std::vector<apf::MeshEntity*> s2cv_map(e_lens.new_cell0_id.size());
  std::map<int, apf::MeshEntity*> c2sv_map{};
  std::vector<apf::MeshEntity*> p2cv_map(e_lens.new_cell1_id.size());
  std::map<int, apf::MeshEntity*> c2pv_map{};

  std::vector<apf::MeshEntity*> v_ctr_new(e_lens.new_cell0_id.size());
  std::vector<apf::MeshEntity*> v_sp_new(e_lens.new_cell1_id.size());
  apf::MeshEntity* v2c_new;

  std::vector<std::vector<apf::MeshEntity*> > p_tets_s1(e_lens.new_cell1_id.size(), std::vector<apf::MeshEntity*>(0));
  std::vector<std::vector<apf::MeshEntity*> > p_tets_s2(e_lens.new_cell1_id.size(), std::vector<apf::MeshEntity*>(0));

  std::vector<std::vector<apf::MeshEntity*> > s_tets(e_lens.new_cell0_id.size(), std::vector<apf::MeshEntity*>(0));

  // For each slice,
  // Create a new 0cell vertex, positioned at old position
  for(int i = 0; i < e_lens.new_cell0_id.size(); i++) {
    int tag = e_lens.new_cell0_id.at(i);
    apf::ModelEntity* vert_em = m_trial->findModelEntity(0,tag);
    v_ctr_new.at(i) = m_trial->createVert(vert_em);
    m_trial->setPoint(v_ctr_new.at(i), 0, pos_old);
  }
  for(int i = 0; i < e_lens.new_cell1_id.size(); i++) {
    int tag = e_lens.new_cell1_id.at(i);
    apf::ModelEntity* vert_em = m_trial->findModelEntity(1,tag);
    v_sp_new.at(i) = m_trial->createVert(vert_em);
    m_trial->setPoint(v_sp_new.at(i), 0, pos_old);
  }
  if(trial_type == 2) {
    int tag = e_lens.new_cell2_id;
    apf::ModelEntity* vert_em = m_trial->findModelEntity(2,tag);
    v2c_new = m_trial->createVert(vert_em);
    m_trial->setPoint(v2c_new, 0, pos_old);
  }

  collect_disc();
  for(int path = 0; path < e_lens.discs.size(); path++) {
    p_tets_s1.at(path).resize(e_lens.discs.at(path).tri.size());
    p_tets_s2.at(path).resize(e_lens.discs.at(path).tri.size());
  }

  for(int slice = 0; slice < slice_cells->size(); slice++) {
    if(slice_cells->at(0).first.first == 0) {
      s_tets.at(slice).reserve(e_lens.discs.at(0).tri.size());
    }
    else {
      int d1 = slice_cells->at(slice).first.first - 1;
      int d2 = slice_cells->at(slice).first.second - 1;

      s_tets.at(slice).reserve(e_lens.discs.at(d1).tri.size() + e_lens.discs.at(d2).tri.size() + 4);
    }
  }

  // Copy the slice entities on the m_precond to the m_trial by replacing the 
  // center vertex with the new 0cell vertex
  for(int slice = 0; slice < slice_cells->size(); slice++) {
    copy_slice_ents_wg(slice, v_ctr_new, p2t_ctr, t2p_ctr, s_map, slice_cells);
    m_trial->acceptChanges();
    std::string vtk_name = "./output/v_trial"+ std::to_string(slice);
    if(sub_vtk)
      vd_save_vtk_vert(m_trial, v_ctr_new.at(slice), vtk_name.c_str());
  }


  // Expand the paths:
  for(int path = 0; path < e_lens.discs.size(); path++) {
    copy_path_disc_wg(path, v_ctr_new, p2t_ctr, t2p_ctr, 
                            v_sp_new, p2t_sp, t2p_sp,
                     s_map, s2pv_map, p2sv_map, p_tets_s1, p_tets_s2, s_tets);
    m_trial->acceptChanges();

    int s1 = e_lens.slices.at(2*path) - 1;
    int s2 = e_lens.slices.at(2*path+1) - 1;
    if(sub_vtk) {
      std::string vtk_name = "./output/v_trial"+ std::to_string(s1);
      vd_save_vtk_vert(m_trial, v_ctr_new.at(s1), vtk_name.c_str());
      vtk_name = "./output/v_trial"+ std::to_string(s2);
      vd_save_vtk_vert(m_trial, v_ctr_new.at(s2), vtk_name.c_str());
    }
  }
  expand_2c_void_wg(path_cells, slice_cells, v_ctr_new, v_sp_new, 
                    v2c_new, p2t_ctr, p2t_sp,
                    p2c_2c, p2sv_map, s2pv_map, s2cv_map, c2sv_map,
                    p2cv_map, c2pv_map, s_map, s_tets);

  f_calc.vd_att_fields(m_trial);

  std::vector<apf::Vector3> vel_ctr(v_ctr_new.size(), 
                                          apf::Vector3(0,0,0));
  std::vector<apf::Vector3> force_ctr(v_ctr_new.size(), 
                                          apf::Vector3(0,0,0));
  std::vector<apf::Vector3> vel_sp(v_sp_new.size(), 
                                          apf::Vector3(0,0,0));
  std::vector<apf::Vector3> force_sp(v_sp_new.size(), 
                                          apf::Vector3(0,0,0));
  apf::Vector3 vel_2c(0,0,0);
  apf::Vector3 force_2c(0,0,0);

  upd_pos_sph_wg(slice_cells, v_ctr_new, v_sp_new, v2c_new, s2pv_map, p2sv_map,
                 s2cv_map, c2sv_map, p2cv_map, c2pv_map,
                    vel_ctr, force_ctr, vel_sp, force_sp,
                    vel_2c, force_2c, p_tets_s1, p_tets_s2, s_tets);

  double p_diss = 0;
  if(ins_flag) {
    if(trial_type == 1) {
      for (int i = 0; i < v_ctr_new.size(); i++) {
        p_diss = p_diss - force_ctr.at(i)*vel_ctr.at(i);
      }
      for (int i = 0; i < v_sp_new.size(); i++) {
        p_diss = p_diss - force_sp.at(i)*vel_sp.at(i);
      }
      w1.at(trial_curr) = p_diss;
    }
    else {
      assert(trial_type == 2);
      double p_diss = 0;
      for (int i = 0; i < v_ctr_new.size(); i++) {
        p_diss = p_diss - force_ctr.at(i)*vel_ctr.at(i);
      }
      for (int i = 0; i < v_sp_new.size(); i++) {
        p_diss = p_diss - force_sp.at(i)*vel_sp.at(i);
      }
      //p_diss = p_diss - force_2c*vel_2c;
      w2.at(trial_curr).at(ng.curr) = p_diss;
    }
  }
  else {
    if(trial_type == 1) {
      w1.at(trial_curr) = 0.1;
    }
    else {
      assert(trial_type == 2);
      w2.at(trial_curr).at(ng.curr) = 0.1;
    }
  }


  trial_load = true;
}


// Destroy the trial mesh and reload.
void vd_edisc::reload_trial() {
  assert(m_ex and m_load);
  assert(precond_ex and precond_load);

  //std::cout << "Mesh cond: "<< trial_ex << " " << trial_load << std::endl;
  //vd_print_ent(m_precond);

  if(!trial_ex) {
    assert(!trial_load);

    m_trial = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);
    trial_ex = true;
    std::cout << "Created empty mesh, " << trial_load << std::endl;
  }

  if(trial_load) {
    apf::destroyMesh(m_trial);
    m_trial = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);

    trial_load = false;
    std::cout << "Recreate empty mesh, " << trial_load << std::endl;
  }

  std::cout << "Initial:" << m_main << ", Precond:" << m_precond
            << ", Trial:" << m_trial << std::endl;
  // Copy the vertices to the new trial mesh.

  vert_trial.clear();
  vert_trial.reserve(vert.size());

  std::cout << "Copying the vertices from the preconditioned mesh." 
            << std::endl;

  apf::MeshIterator* it_e = m_precond->begin(0);
  //std::cout << "Iterator." << std::endl;
  apf::MeshEntity* e;

  trial2pre_map.clear();
  pre2trial_map.clear();
  while (e = m_precond->iterate(it_e)) {
    //std::cout << "Vert " << e << "["<< m_precond->toModel(e) << "] ," ;

    int type = m_precond->getModelType(m_precond->toModel(e));
    int tag = m_precond->getModelTag(m_precond->toModel(e));
    apf::ModelEntity* vert_em = m_trial->findModelEntity(type,tag);

    vert_trial.push_back(m_trial->createVert(vert_em));

    trial2pre_map[vert_trial.back()] = e;
    pre2trial_map[e] = vert_trial.back();
    //std::cout << "Vert precond " << e << " "
    //          << " trial " << vert_trial.back() << std::endl;

    apf::Vector3 vert_pos(0,0,0);
    m_precond->getPoint(e, 0, vert_pos);
    //std::cout<< vert_trial.back() << ", " << vert_pos << std::endl;
    m_trial->setPoint(vert_trial.back(), 0, vert_pos);

  }
  m_precond->end(it_e);

  //std::cout << vert_trial.size() << " " << vert.size() << std::endl;
  assert(vert_trial.size() == vert.size());
  //std::cout << "Copied." << std::endl;

  // Going over the tets of the preconditioned mesh, copy them recursively to 
  // the trial mesh.
  it_e = m_precond->begin(3);

  while (e = m_precond->iterate(it_e)) {
    copyElement(e,0);
  }
  m_precond->end(it_e);

  m_trial->acceptChanges();
  f_calc.vd_att_fields(m_trial);
  // m_trial->verify();

  // Make sure that the cell complex has enough free cells.
  // set_free_cells();

  std::vector<apf::MeshEntity*> es_vert_cp(0);
  // Get the central vertex in the trial mesh.
  vd_find_ent_geom(m_trial, &es_vert_cp, cell_id, 0);
  assert(es_vert_cp.size() == 1);
  vert_ctr_cp = es_vert_cp.at(0);

  f_calc.vdparam.adj_dt(dt_max);

  act_trial();
  trial_load = true;

  //collect_pc();
  reload_edges();
}

// Collect paths and circuits.
void vd_edisc::collect_pc() {

  //clear_path();
  clear_edges();

  if(m_act == m_precond)
    precond_map();
  else {
    trans_map_pre2act();
    vd_cd->set_vert_c1_map(actid_map);
    vd_cd->set_cb_flag(false);
    vd_cd->reload(m_act, c_base_act, vert_ctr_act);
  }
  c2_edge = vd_cd->get_c2_edges();
  c3_edge = vd_cd->get_c3_edges();

  //std::cout << "2cell edges are ";
  //for(int i = 0; i < c2_edge.size(); i++)
  //  std::cout << c2_edge.at(i) << " ";
  //std::cout << std::endl;

  circ_sz = vd_cd->get_circ_sz();
  pt_sz = vd_cd->get_path_sz();

  clear_ext();
  ext_0cell = vd_cd->get_ext();
}

// Collect paths and circuits.
void vd_edisc::overwrite_pc(vd_disc_cut* vd_cut) {

  //clear_path();
  clear_edges();

  assert(m_act != m_precond);
  trans_map_pre2act();
  vd_cd->set_vert_c1_map(actid_map);
  vd_cd->set_cb_flag(false);

  vd_cd->reload(m_act, c_base_act, vert_ctr_act);
  set_the_edges(vd_cut);

  c2_edge = vd_cd->get_c2_edges();
  c3_edge = vd_cd->get_c3_edges();

  //std::cout << "2cell edges are ";
  //for(int i = 0; i < c2_edge.size(); i++)
  //  std::cout << c2_edge.at(i) << " ";
  //std::cout << std::endl;

  circ_sz = vd_cd->get_circ_sz();
  pt_sz = vd_cd->get_path_sz();

  ext_new.clear();
  ext_slice.clear();
  ext_pc.clear();
  ext_cor.clear();
  ext_0cell = false;

  ext_0cell = vd_cd->get_ext();
}

void vd_edisc::refresh_e() {

  w1.resize(circ_sz);
  e1.resize(circ_sz);

  int step = 10;
  w1_exp.resize(circ_sz);
  for(int i = 0; i < circ_sz; i++)
    w1_exp.at(i).resize(step);

  w2.resize(pt_sz);
  e2.resize(pt_sz);
  w2_exp.resize(pt_sz);

  for(trial_curr = 0; trial_curr < pt_sz; trial_curr++) {
    // Collect the relevant cells for the current 3cell couple to be joined.
    ng.clear();
    get_path_ngon(trial_curr, &ng);
    ng.print();
    int ngon_sz = ng.ngons.size();

    w2.at(trial_curr).resize(ngon_sz);
    e2.at(trial_curr).resize(ngon_sz);
    w2_exp.at(trial_curr).resize(ngon_sz);
  }
}

bool vd_edisc::burn_the_cut_c2(int k, vd_disc_cut* vd_cut, vd_plane* pl) {

  int c2_curr;

  if(trial_type == 1) {
    c2_curr = vd_cut->circ_tup->first.second.at(k);
  }
  else {
    c2_curr = vd_cut->ng_pt->second.at(k);
  }

  vd_inter v_int;

  apf::Downward d_e;
  apf::Downward d_v;

  vd_cut->ents = vd_cd->get_dim_ent(2, c2_curr);

  apf::ModelEntity* mdl_curr = e_lens.m->toModel(vd_cut->ents.at(0));
  std::cout << "Trying to cut 2c" << e_lens.m->getModelTag(mdl_curr)
            << std::endl;

  bool cut_once = false;
  bool cut_twice = false;
  int cut_curr = -1;
  bool fail = false;
  for(int e = 0; e < vd_cut->ents.size(); e++) {
    e_lens.m->getDownward(vd_cut->ents.at(e), 0, d_v);
    e_lens.m->getDownward(vd_cut->ents.at(e), 1, d_e);
    int v1 = findIn(d_v, 3, e_lens.vert_ctr);
    pl_int_edge(e_lens.m, d_e[lookup_v_tri_e_x[v1]], pl, &v_int);

    if(v_int.cut == 1) {
      if(cut_once) {
        cut_twice = true;
        fail = true;
      }
      else {
        cut_once = true;
        vd_cut->c2_edges.at(k) = d_e[lookup_v_tri_e_x[v1]];
        vd_cut->c2_cuts.at(k) = v_int.cut;
        vd_cut->c2_cut_pos.at(k) = v_int.pos;
        cut_curr = e;
      }
    }
    else if(v_int.cut < 0) {
      if(cut_twice) {
        fail = true;
      }
      else if(cut_once) {
        if(vd_cut->c2_cuts.at(k) < 0) {

          if(vd_cut->c2_edges.at(k) != d_e[lookup_v_tri_e_x[v1]])
            fail = true;
          cut_twice = true;
        }
        else
          fail = true;
      }
      else {
        cut_once = true;

        vd_cut->c2_cuts.at(k) = v_int.cut;
        vd_cut->c2_cut_pos.at(k) = v_int.pos;
        cut_curr = e;

        e_lens.m->getDownward(d_e[lookup_v_tri_e_x[v1]], 0, d_v);
        apf::MeshEntity* e_bound;
        int v_id = (vd_cut->c2_cuts.at(k)+1)*(-1);
        assert(vd_find_edge(e_lens.m, d_v[v_id], e_lens.vert_ctr, &e_bound));
        if(e_lens.m->toModel(e_bound) != mdl_curr)
          fail = true;
        else
          vd_cut->c2_edge_new.at(c2_curr) = e_bound;

        vd_cut->c2_edges.at(k) = e_bound;
      }
    }
    else {
      if(!cut_once) {
        vd_cut->c2_cuts.at(k) = v_int.cut;
      }
    }
  }
  if(!cut_once)
    fail = true;

  if(fail) {
    std::cout << "Cutting not successful." << std::endl;
    return false;
  }
  else {
    assert(cut_curr > -1);
    std::cout << "Cutting type" << vd_cut->c2_cuts.at(k) 
              << std::endl;
    if(vd_cut->c2_cuts.at(k) == 1) {
      vd_lens* lens_sph = new vd_lens(e_lens.m, e_lens.c_base_curr);
      lens_sph->load_edge(vd_cut->c2_edges.at(k));
      lens_sph->split_lens();
      e_lens.m->acceptChanges();
      vd_chk_neg_sgn(e_lens.m);

      f_calc.vd_att_fields(e_lens.m, lens_sph->get_vert_ctr());

      apf::MeshEntity* v_new = lens_sph->get_vert_ctr();
      e_lens.m->setPoint(v_new, 0, vd_cut->c2_cut_pos.at(k));
      assert(vd_find_edge(e_lens.m, v_new, e_lens.vert_ctr, 
           &vd_cut->c2_edge_new.at(c2_curr) ));

      vd_cut->c2_edge_flag.at(c2_curr) = true;
      delete lens_sph;
    }
    else {
      assert(vd_cut->c2_cuts.at(k) < 0);
    }

    return true;
  }
}

// Either return NULL, an edge or tri. 
apf::MeshEntity* vd_edisc::burn_c3_edge(apf::MeshEntity* e_curr, 
                                     vd_disc_cut* vd_cut, vd_plane* pl) {
  int e_type = e_lens.m->getType(e_curr);
  assert(e_type == apf::Mesh::EDGE);

  apf::MeshElement* ee = createMeshElement(e_lens.m, e_curr);
  double meas = measure(ee);
  destroyMeshElement(ee);
  assert(meas > std::numeric_limits<double>::min());

  std::vector<apf::MeshEntity*> tri_set(0);
  std::vector<apf::MeshEntity*> tet_set(0);

  vd_set_up(e_lens.m, e_curr, &tri_set);
  vd_set_up(e_lens.m, &tri_set, &tet_set);

  bool cut_once = false;
  bool cut_twice = false;
  bool fail = false;

  apf::MeshEntity* e_next = NULL;
  int cut_type = 0;

  apf::Downward d_v;
  apf::Downward d_e;
  apf::Downward d_t;

  vd_inter v_int;

  for(int e = 0; e < tet_set.size(); e++) {
    std::cout << tet_set.at(e) << " " << vd_cut->ent_burn[tet_set.at(e)]
             << std::endl;
    if(!vd_cut->ent_burn[tet_set.at(e)]) {
      if(vd_cut->mdl_curr == e_lens.m->toModel(tet_set.at(e))) {
        assert(!(cut_type == 1 and cut_twice));
        e_lens.m->getDownward(tet_set.at(e), 0, d_v);
        e_lens.m->getDownward(tet_set.at(e), 1, d_e);
        e_lens.m->getDownward(tet_set.at(e), 2, d_t);
        int v1_tet = findIn(d_v, 4, e_lens.vert_ctr);
        int e1_tet = findIn(d_e, 6, e_curr);
        assert(v1_tet > -1 and e1_tet > -1);
        int i1_tet = lookup_v_tet_e_t_x[v1_tet][e1_tet];
        assert(i1_tet > -1);

        std::cout << "v1 " << v1_tet << " e1 " << e1_tet
                 << " t1 " << i1_tet
                 << std::endl;

        e_lens.m->getDownward(d_t[i1_tet], 0, d_v);
        e_lens.m->getDownward(d_t[i1_tet], 1, d_e);
        int v1_tri = findIn(d_v, 3, e_lens.vert_ctr);

        pl_int_edge(e_lens.m, d_e[lookup_v_tri_e_x[v1_tri]], pl, &v_int);

        if(v_int.cut == 1) {
          if(cut_once) {
            cut_twice = true;
            fail = true;
          }
          else {
            cut_once = true;
            cut_type = 1;
            e_next = d_t[i1_tet];
          }
        }
        else if(v_int.cut < 0) {
          if(cut_twice) {
            fail = true;
          }
          if(cut_once) {
            if(cut_type < 0)
              cut_twice = true;
            else
              fail = true;
          }
          else {
            cut_once = true;
            cut_type = v_int.cut;

            e_lens.m->getDownward(d_e[lookup_v_tri_e_x[v1_tri]], 0, d_v);
            int v_id = (v_int.cut+1)*(-1);
            vd_find_edge(e_lens.m, d_v[v_id], e_lens.vert_ctr, &e_next);
          }
        }
      }
      vd_cut->ent_burn[tet_set.at(e)] = true;
    }
  }
  if(!cut_once)
    fail = true;
  if(fail)
    return NULL;
  else
    return e_next;
}

apf::MeshEntity* vd_edisc::burn_c3_tri(apf::MeshEntity* e_curr, 
                                     vd_disc_cut* vd_cut, vd_plane* pl) {
  int e_type = e_lens.m->getType(e_curr);
  assert(e_type == apf::Mesh::TRIANGLE);

  apf::MeshElement* ee = createMeshElement(e_lens.m, e_curr);
  double meas = measure(ee);
  destroyMeshElement(ee);
  assert(meas > std::numeric_limits<double>::min());

  std::vector<apf::MeshEntity*> tet_set(0);

  vd_set_up(e_lens.m, e_curr, &tet_set);

  bool cut_once = false;
  bool cut_twice = false;
  bool fail = false;

  int cut_type = 0;

  apf::Downward d_v;
  apf::Downward d_e;
  apf::Downward d_t;

  vd_inter v_int;

  apf::MeshEntity* e_next = NULL;

  for(int e = 0; e < tet_set.size(); e++) {
    if(!vd_cut->ent_burn[tet_set.at(e)]) {

      if(vd_cut->mdl_curr == e_lens.m->toModel(tet_set.at(e)) ) {
        assert(!cut_once);
        e_lens.m->getDownward(tet_set.at(e), 0, d_v);
        e_lens.m->getDownward(tet_set.at(e), 2, d_t);
        int v1_tet = findIn(d_v, 4, e_lens.vert_ctr);
        int t1_tet = findIn(d_t, 4, e_curr);
        assert(v1_tet > -1 and t1_tet > -1);
        int i1_tet = lookup_v_tet_t_t_x[v1_tet][t1_tet][0];
        int i2_tet = lookup_v_tet_t_t_x[v1_tet][t1_tet][1];
        assert(i1_tet > -1 and i2_tet > -1);

        e_lens.m->getDownward(d_t[i1_tet], 0, d_v);
        e_lens.m->getDownward(d_t[i1_tet], 1, d_e);
        int v1_tri = findIn(d_v, 3, e_lens.vert_ctr);

        pl_int_edge(e_lens.m, d_e[lookup_v_tri_e_x[v1_tri]], pl, &v_int);

        if(v_int.cut == 1) {
          cut_once = true;
          cut_type = 1;
          e_next = d_t[i1_tet];
        }
        else if(v_int.cut < 0) {
          cut_once = true;
          cut_type = v_int.cut;

          e_lens.m->getDownward(d_e[lookup_v_tri_e_x[v1_tri]], 0, d_v);
          int v_id = (v_int.cut+1)*(-1);
          vd_find_edge(e_lens.m, d_v[v_id], e_lens.vert_ctr, &e_next);
        }
        else {
          e_lens.m->getDownward(d_t[i2_tet], 0, d_v);
          e_lens.m->getDownward(d_t[i2_tet], 1, d_e);
          int v1_tri = findIn(d_v, 3, e_lens.vert_ctr);

          pl_int_edge(e_lens.m, d_e[lookup_v_tri_e_x[v1_tri]], pl, &v_int);
          if(v_int.cut != 1)
            return NULL;

          cut_once = true;
          cut_type = 1;
          e_next = d_t[i1_tet];
        }
      }
      vd_cut->ent_burn[tet_set.at(e)] = true;
    }
  }

  if(fail)
    return NULL;
  else
    return e_next;
}

// Around the central vertex, given a 2cell edge, going over the 2cell edges 
// and tris, find the bounding 1cell edges. Assert found. 
std::pair<apf::MeshEntity*, apf::MeshEntity*> vd_edisc::get_bound
                                          (apf::MeshEntity* edge_c2) {

  apf::Up up;
  apf::Downward v;
  apf::Downward e;

  std::pair<apf::MeshEntity*, apf::MeshEntity*> bounds;
  bool found1 = false;
  bool found2 = false;

  std::vector<std::vector <apf::MeshEntity* > > e_fire(2, 
                                          std::vector <apf::MeshEntity* >(0));
  // Reserve space for all edges.
  e_lens.m->getUp(e_lens.vert_ctr, up);
  e_fire.at(0).reserve(up.n);
  e_fire.at(1).reserve(up.n);

  e_fire.at(0).push_back(edge_c2);

  apf::MeshEntity* e_curr;
  apf::MeshEntity* e_next;

  std::map<apf::MeshEntity*,bool> e_burn{};
  e_burn[edge_c2] = true;

  apf::ModelEntity* mdl_c2 = e_lens.m->toModel(edge_c2);
  while(!(found1 and found2) ) {
    while(e_fire.at(0).size() > 0) {

      e_curr = e_fire.at(0).back();

      e_lens.m->getUp(e_curr, up);
      for(int i = 0; i < up.n; i++) {
        apf::ModelEntity* mdl_curr = e_lens.m->toModel(up.e[i]);
        if (mdl_curr == mdl_c2 and !e_burn[up.e[i]]) {
          e_fire.at(1).push_back(up.e[i]);
        }
        e_burn[up.e[i]] = true;
      }
      e_fire.at(0).pop_back();
    }

    while(e_fire.at(1).size() > 0) {

      e_curr = e_fire.at(1).back();
      e_lens.m->getDownward(e_curr, 0, v);
      e_lens.m->getDownward(e_curr, 1, e);

      int v1 = findIn(v, 3, e_lens.vert_ctr);
      assert(v1 > -1);

      if(!e_burn[e[lookup_tri_ed[v1][0]]] ) {
        assert(e_burn[e[lookup_tri_ed[v1][1]]]);
        e_next = e[lookup_tri_ed[v1][0]];
      }
      else {
        e_next = e[lookup_tri_ed[v1][1]];
      }

      apf::ModelEntity* mdl_curr = e_lens.m->toModel(e_next);
      if (mdl_curr == mdl_c2) {
        e_fire.at(0).push_back(e_next);
      }
      else {
        if(!found1) {
          bounds.first = e_next;
          found1 = true;
        }
        else {
          assert(!found2);
          bounds.second = e_next;
          found2 = true;
        }
      }
      e_burn[e_next] = true;

      e_fire.at(1).pop_back();
    }
  }

  assert(found1 and found2);
  return bounds;
}

bool vd_edisc::burn_the_cut_c3(int k, vd_disc_cut* vd_cut) {

  vd_inter v_int;
  int c3_curr;
  int c2_f;
  int c2_e;

  if(trial_type == 1) {
    assert(k < vd_cut->circ_tup->first.first.size());
    if(k < vd_cut->circ_tup->first.first.size()-1) {
      c3_curr = vd_cut->circ_tup->first.first.at(k);
      c2_f = vd_cut->circ_tup->first.second.at(k);
      c2_e = vd_cut->circ_tup->first.second.at(k+1);
      if(c3_curr == -1)
        return true;
    }
    else {
      c3_curr = vd_cut->circ_tup->first.first.back();
      c2_f = vd_cut->circ_tup->first.second.back();
      c2_e = vd_cut->circ_tup->first.second.at(0);
      if(c3_curr == -1)
        return true;
    }
  }
  else {
    assert(k < vd_cut->ng_pt->first.size()-1 and k > 0);
    c3_curr = vd_cut->ng_pt->first.at(k);
    c2_f = vd_cut->ng_pt->second.at(k-1);
    c2_e = vd_cut->ng_pt->second.at(k);
    if(c3_curr == -1)
      return true;
  }

  std::cout << "c2_f " << c2_f
            << " c3 " << c3_curr 
            << " c2_e " << c2_e 
            << std::endl;
  vd_plane pl;
  pl.pos = pos_old;
  pl.norm = vd_get_e_plane(e_lens.m, e_lens.vert_ctr, c2_edge.at(c2_f), 
                                                    c2_edge.at(c2_e));

  //Check if the plane passes between the bounding edges of c2_f and c2_e.
  std::pair<apf::MeshEntity*, apf::MeshEntity*> bounds;
  bounds = get_bound(c2_edge.at(c2_f));

  apf::Vector3 p1(0,0,0);
  apf::Vector3 p2(0,0,0);
  e_lens.m->getPoint(getEdgeVertOppositeVert(e_lens.m, bounds.first, 
                      e_lens.vert_ctr), 0, p1);
  e_lens.m->getPoint(getEdgeVertOppositeVert(e_lens.m, bounds.second, 
                      e_lens.vert_ctr), 0, p2);

  pl_int_2pts(e_lens.m, p1, p2, &pl, &v_int);
  bool cut1 = (v_int.cut == 1);

  bounds = get_bound(c2_edge.at(c2_e));
  e_lens.m->getPoint(getEdgeVertOppositeVert(e_lens.m, bounds.first, 
                      e_lens.vert_ctr), 0, p1);
  e_lens.m->getPoint(getEdgeVertOppositeVert(e_lens.m, bounds.second, 
                      e_lens.vert_ctr), 0, p2);

  pl_int_2pts(e_lens.m, p1, p2, &pl, &v_int);
  bool cut2 = (v_int.cut == 1);

  std::cout << "normal plane between edges " << pl.norm << std::endl;
  if(!cut1 or !cut2) {
    std::cout << "Cutting not successful." << std::endl;
    return false;
  }

  vd_cut->mdl_curr = vd_cut->c3_mdl.at(c3_curr);
  vd_cut->ents = vd_cd->get_dim_ent(3, c3_curr);

  std::vector<apf::MeshEntity*> edge_cut(0);

  edge_cut.reserve(vd_cut->ents.size()*2);

  assert(vd_cut->ents.size() > 0);
  apf::ModelEntity* mdl_end = e_lens.m->toModel(c2_edge.at(c2_e));
  // Start cutting from the first c2 edge.
  // Then either cut through tris, or if cutting through an edge, cut from an
  // edge, until the end c2 is reached.

  apf::MeshEntity* e_curr = c2_edge.at(c2_f);

  apf::Downward d_v;
  apf::Downward d_e;

  bool end = false;
  bool fail = false;
  while(!end) {
    int e_type = e_lens.m->getType(e_curr);
    std::cout << e_type << "-dim " << e_curr << std::endl;
    if(e_type == apf::Mesh::EDGE) {
      e_curr = burn_c3_edge(e_curr, vd_cut, &pl);
    }
    else {
      assert(e_type == apf::Mesh::TRIANGLE);
      e_lens.m->getDownward(e_curr, 0, d_v);
      e_lens.m->getDownward(e_curr, 1, d_e);
      int v1 = findIn(d_v, 3, e_lens.vert_ctr);

      edge_cut.push_back(d_e[lookup_v_tri_e_x[v1]]);

      e_curr = burn_c3_tri(e_curr, vd_cut, &pl);
    }

    if(!e_curr) {
      std::cout << "Cutting through 3c" 
                << e_lens.m->getModelTag(vd_cut->mdl_curr)
                << " not successful" << std::endl;
      end = true;
      fail = true;
    }
    else {
      if(e_lens.m->toModel(e_curr) == mdl_end)
        end = true;
    }
  }

  if(fail) {
    std::cout << "Cutting not successful." << std::endl;
    return false;
  }
  else {
    std::cout << "Cutting type 3c"
              << std::endl;
    for(int i = 0; i < edge_cut.size(); i++) {
      if(pl_int_edge(e_lens.m, edge_cut.at(i), pl.pos, pl.norm, 0.05) == 1) {
        pl_int_edge(e_lens.m, edge_cut.at(i), &pl, &v_int);
        vd_lens* lens_sph = new vd_lens(e_lens.m, e_lens.c_base_curr);
        lens_sph->load_edge(edge_cut.at(i));
        lens_sph->split_lens();
        e_lens.m->acceptChanges();
        vd_chk_neg_sgn(e_lens.m);
        f_calc.vd_att_fields(e_lens.m, lens_sph->get_vert_ctr());

        apf::MeshEntity* v_new = lens_sph->get_vert_ctr();
        e_lens.m->setPoint(v_new, 0, v_int.pos);
        delete lens_sph;
      }
    }
    return true;
  }
}

bool vd_edisc::burn_the_cut_c3_old(int k, vd_disc_cut* vd_cut, vd_plane* pl) {

  vd_inter v_int;
  int c3_curr;
  int c2_f;
  int c2_e;

  if(trial_type == 1) {
    assert(k < vd_cut->circ_tup->first.first.size());
    if(k < vd_cut->circ_tup->first.first.size()-1) {
      c3_curr = vd_cut->circ_tup->first.first.at(k);
      c2_f = vd_cut->circ_tup->first.second.at(k);
      c2_e = vd_cut->circ_tup->first.second.at(k+1);
      if(c3_curr == -1)
        return true;
    }
    else {
      c3_curr = vd_cut->circ_tup->first.first.back();
      c2_f = vd_cut->circ_tup->first.second.back();
      c2_e = vd_cut->circ_tup->first.second.at(0);
      if(c3_curr == -1)
        return true;
    }
  }
  else {
    assert(k < vd_cut->ng_pt->first.size()-1 and k > 0);
    c3_curr = vd_cut->ng_pt->first.at(k);
    c2_f = vd_cut->ng_pt->second.at(k-1);
    c2_e = vd_cut->ng_pt->second.at(k);
    if(c3_curr == -1)
      return true;
  }

  vd_cut->mdl_curr = vd_cut->c3_mdl.at(c3_curr);
  vd_cut->ents = vd_cd->get_dim_ent(3, c3_curr);

  std::vector<apf::MeshEntity*> edge_cut(0);

  edge_cut.reserve(vd_cut->ents.size()*2);

  assert(vd_cut->ents.size() > 0);
  apf::ModelEntity* mdl_end = e_lens.m->toModel(vd_cut->c2_edge_new.at(c2_e));
  // Start cutting from the first c2 edge.
  // Then either cut through tris, or if cutting through an edge, cut from an
  // edge, until the end c2 is reached.

  apf::MeshEntity* e_curr = vd_cut->c2_edge_new.at(c2_f);

  apf::Downward d_v;
  apf::Downward d_e;

  bool end = false;
  bool fail = false;
  while(!end) {
    int e_type = e_lens.m->getType(e_curr);
    std::cout << e_type << "-dim " << e_curr << std::endl;
    if(e_type == apf::Mesh::EDGE) {
      e_curr = burn_c3_edge(e_curr, vd_cut, pl);
    }
    else {
      assert(e_type == apf::Mesh::TRIANGLE);
      e_lens.m->getDownward(e_curr, 0, d_v);
      e_lens.m->getDownward(e_curr, 1, d_e);
      int v1 = findIn(d_v, 3, e_lens.vert_ctr);

      // Cut the edge across the triangle intersecting the plane.
      edge_cut.push_back(d_e[lookup_v_tri_e_x[v1]]);

      e_curr = burn_c3_tri(e_curr, vd_cut, pl);
    }

    if(!e_curr) {
      end = true;
      fail = true;
    }
    else {
      if(e_lens.m->toModel(e_curr) == mdl_end)
        end = true;
    }
  }

  if(fail) {
    std::cout << "Cutting not successful." << std::endl;
    return false;
  }
  else {
    std::cout << "Cutting type 3c"
              << std::endl;
    for(int i = 0; i < edge_cut.size(); i++) {
      pl_int_edge(e_lens.m, edge_cut.at(i), pl, &v_int);

      vd_lens* lens_sph = new vd_lens(e_lens.m, e_lens.c_base_curr);
      lens_sph->load_edge(edge_cut.at(i));
      lens_sph->split_lens();
      e_lens.m->acceptChanges();
      vd_chk_neg_sgn(e_lens.m);
      f_calc.vd_att_fields(m_act, lens_sph->get_vert_ctr());

      apf::MeshEntity* v_new = lens_sph->get_vert_ctr();
      e_lens.m->setPoint(v_new, 0, v_int.pos);
      delete lens_sph;
    }
    return true;
  }
}

// TODO Placeholder for another possible disc fixing function.
// Try to cut singular creases if they exist and cutting fixes inversion.
// Currently, if a single plane cannot cut through all 3cells on the disc in  
// a continuous fashion, the insertion fails. If this creates problems with
// physical insertions, it should be mitigated. The rationale behind skipping
// it for now is the evolution tends to simply the geometry for the physical
// insertion and the curvature would allow for more volume in the 3cells on the
// circuit of the physical insertion.
bool vd_edisc::cut_sing_crease(vd_disc_cut* vd_cut) {
  return false;
}

bool vd_edisc::do_the_cuts(vd_disc_cut* vd_cut) {

  assert(trial_type == 1 or trial_type == 2);

  vd_cut->c2_edge_new.clear();
  vd_cut->c2_edge_new.resize(vd_cd->get_2c_nbr());

  vd_cut->c2_edge_flag.clear();
  vd_cut->c2_edge_flag.resize(vd_cd->get_2c_nbr());

  vd_cut->c3_mdl.resize(vd_cd->get_3c_nbr());
  for(int i = 0; i < vd_cut->c3_mdl.size(); i++) {
    vd_cut->c3_mdl.at(i) = vd_cd->get_mdl(3, i);
  }

  for(int i = 0; i < vd_cut->c2_edge_flag.size(); i++)
    vd_cut->c2_edge_flag.at(i) = false;

  if(trial_type == 1) {
    vd_cut->c2_cuts.clear();
    vd_cut->c2_edges.clear();
    vd_cut->c2_cut_pos.clear();

    vd_cut->circ_tup = get_circ(trial_curr);

    vd_cut->c2_cuts.resize(vd_cut->circ_tup->first.second.size());
    vd_cut->c2_edges.resize(vd_cut->circ_tup->first.second.size());
    vd_cut->c2_cut_pos.resize(vd_cut->circ_tup->first.second.size());

    //for(int k = 0; k < vd_cut->circ_tup->first.second.size(); k++) {
    //  if(!burn_the_cut_c2(k, vd_cut, &pl))
    //    return false;
    //}

    if(!cut_sing_crease(vd_cut)) {
      // Otherwise try to cut through a single plane for all 3cells on the 
      // circuit.
      for(int k = 0; k < vd_cut->circ_tup->first.first.size(); k++) {
        if(!burn_the_cut_c3(k, vd_cut))
          return false;
      }
    }
  }
  else {
    for (int i = 0; i < e_lens.discs.size(); i++) {
      int v1 = e_lens.slices.at(2*i)-1;
      int v2 = e_lens.slices.at(2*i+1)-1;

      std::cout << i << " " 
                << e_lens.vt.at(v1) << " "  << e_lens.vt.at(v2) << std::endl;

      apf::Vector3 mot_rel = e_lens.vt.at(v1) - e_lens.vt.at(v2);
      mot_rel = norm_0(mot_rel);

      if(detect_inv(i, pos_old, mot_rel)) {
        // Iterate over the cells of the current path. 
        // The directions of the top and bottom 3cell edges define an expansion
        // plane.
        // If the projection of the disjoint 2cell entities intersects 
        // itself, the insertion is skipped. It requires spherical insertion. 
        // This is TODO
        // Check if the current disjoint 2cell entities are cut by the plane
        // normal to mot_rel. 
        // If the plane normal to mot_rel does not cut the entities find the 
        // bounding planes defined by the bounding edges and the expansion 
        // plane normal.
        // Going over all disjoint cells on the path, collect if the plane cuts
        // and the bounds. If all cut, collect the edges to be sliced.
        // If there are more than one non-intersecting disjoint cell, check
        // if the bounds overlap. If they do, check if there is a plane that
        // doesn't invert entities w.r.t. mot_rel.
        // This is not guaranteed to work. 
        // Even if all checks are passed it is possible that something inverts 
        // during find_ext_dir. The inversion check in the end should revert
        // things if that happens.
        vd_cut->ng_pt = &ng.cells.at(ng.ngons.at(ng.curr).at(i));

        vd_cut->c2_cuts.clear();
        vd_cut->c2_edges.clear();
        vd_cut->c2_cut_pos.clear();
        vd_cut->c2_cuts.resize(vd_cut->ng_pt->second.size());
        vd_cut->c2_edges.resize(vd_cut->ng_pt->second.size());
        vd_cut->c2_cut_pos.resize(vd_cut->ng_pt->second.size());

        //vd_plane pl;
        //pl.norm = mot_rel;
        //pl.pos = pos_old;

        //for(int k = 0; k < vd_cut->ng_pt->second.size(); k++) {
        //  if(!burn_the_cut_c2(k, vd_cut, &pl))
        //    return false;
        //}

        for(int k = 1; k < vd_cut->ng_pt->first.size()-1; k++) {
          if(!burn_the_cut_c3(k, vd_cut))
            return false;
        }
      }
    }

  }
  return true;
}

// Collect paths and circuits.
void vd_edisc::set_the_edges(vd_disc_cut* vd_cut) {
  for(int i = 0; i < vd_cd->get_2c_nbr(); i++) {
    if(vd_cut->c2_edge_flag.at(i)) {
      vd_cd->set_c2_edge(i, vd_cut->c2_edge_new.at(i));
    }
  }
}

void vd_edisc::reload_edges() {
  clear_edges();

  trans_map_pre2act();
  vd_cd->set_cb_flag(false);
  vd_cd->set_vert_c1_map(actid_map);
  vd_cd->reload(m_act, c_base_act, vert_ctr_act);
  c2_edge = vd_cd->get_c2_edges();
  c3_edge = vd_cd->get_c3_edges();

  clear_ext();
  ext_0cell = vd_cd->get_ext();
}

void vd_edisc::clear_edges() {
  c2_edge.clear();
}

/*
void vd_edisc::clear_path() {

  std::vector<std::pair<std::vector<int>, std::vector<int> > > circ_id;
  std::vector<std::pair<std::vector<int>, std::vector<int> > > path_id;

  std::vector<std::vector<int> > ngons;

}
*/

// Lens and disc ordering. These operations are used to find and prepare the 
// disc to extend the lens.
// For a petal(2cell) removal, find the triangles belonging to the 2cell.
// After finding the triangles, order the triangles and extract the edges 
// bounding the triangles. 
// For each edge, find the surrounding tet lens. The lens is ordered with the 
// vertex on one side and all lenses are ordered in the same sense. 
// Find the intersections of the sequential tet lenses. Join the non-
// intersecting parts. 
// Get the disc (the side triangles of the joint lens containing the vertex).

// For a 3cell containing disjoint graph removal, find the 2cells in the 
// disjoint graph connecting to the circuit 3cells. Find the triangles on 
// these 2cells. Order the triangles. Check if two sequential triangles 
// share a tet. If so, split the tet. These triangles will form the disc.

// For each disc triangle, find the adjacent tetrahedra. Store the 3cell 
// membership of the tet belonging to the circuit 3cell.
// Extend the disc.
// The newly created disc and tets will belong to the stored 3cells.
/*
// Copy the triangles and edges without sorting. They only need to be unique.
void vd_edisc::cp_tri(const std::vector<apf::MeshEntity*>* es_tri, 
                                    vd_disc* disc) {

  apf::Up up;

  disc->broken = false;

  disc->tri.clear();
  disc->edge.clear();

  if(es_tri->size() > 0) {
    int tri_sz = es_tri->size();
    disc->tri.resize(tri_sz);
    disc->edge.reserve(2*tri_sz);

    // If the number of unique edges is equal to that of the triangles, the disc
    // is not broken. If they are +1, the disc is broken. If +2, the disc is two
    // fins. 

    for(int i = 0; i < es_tri->size(); i++) {
      disc->tri.at(i) = es_tri->at(i);
      apf::Downward down;
      apf::Downward v;
      disc->m->getDownward(es_tri->at(i), 1, down);
      disc->m->getDownward(es_tri->at(i), 0, v);
      int v0 = findIn(v, 3, disc->v_ctr);

      disc->edge.push_back(down[lookup_tri_ed[v0][0]]);
      disc->edge.push_back(down[lookup_tri_ed[v0][1]]);
    }

    std::sort(disc->edge.begin(), disc->edge.end());
    std::vector<apf::MeshEntity*>::iterator it;
    it = std::unique (disc->edge.begin(), disc->edge.end());
    disc->edge.resize(std::distance(disc->edge.begin(),it));

    if(disc->edge.size() > disc->tri.size() )
      disc->broken = true;

    disc->t_em.resize(disc->tri.size());
    disc->e_em.resize(disc->edge.size());
    disc->t2_em.resize(disc->edge.size());

    //disc->elem_top.resize(disc->tri.size());
    //disc->elem_bottom.resize(disc->tri.size());

    disc->v_pos.resize(disc->edge.size());

    for(int i = 0; i < disc->v_pos.size(); i++) {
      apf::MeshEntity* v_curr = getEdgeVertOppositeVert(disc->m, 
                                        disc->edge.at(i), disc->v_ctr);
      disc->m->getPoint(v_curr, 0, disc->v_pos.at(i));
      //disc->v_pos.at(i) = disc->v_pos.at(i) - v_pos;
    }
    for (int j = 0; j < disc->edge.size(); j++) {
      apf::MeshEntity* e_curr = disc->edge.at(j);
      e_lens.m->getUp(e_curr, up);
      bool found_1 = false;
      bool found_2 = false;
      for (int k = 0; k < up.n; k++) {
        std::vector<apf::MeshEntity*>::iterator it_e;
        it_e = std::find(disc->tri.begin(), 
                         disc->tri.end(), up.e[k]);
        // If on the disc, replace it with the new edge for the top disc.
        if(it_e != disc->tri.end()) {
          if(!found_1) {
            found_1 = true;
            disc->e_t_map1[e_curr] = std::distance(disc->tri.begin(), it_e)+1;
            if(disc->t_e_map1[*it_e] == 0) {
              disc->t_e_map1[*it_e] = j+1;
            }
            else {
              assert(disc->t_e_map2[*it_e] == 0);
              disc->t_e_map2[*it_e] = j+1;
            }
          }
          else {
            assert(!found_2);
            found_2 = true;
            disc->e_t_map2[e_curr] = std::distance(disc->tri.begin(), it_e)+1;
            if(disc->t_e_map1[*it_e] == 0) {
              disc->t_e_map1[*it_e] = j+1;
            }
            else {
              assert(disc->t_e_map2[*it_e] == 0);
              disc->t_e_map2[*it_e] = j+1;
            }
            k = up.n;
          }
        }
      }
    }
  }
}
*/
// TODO correct rebuilding for 2cell insertion requires:
// - we have a number of paths

// - designate the new 0cells, associated with the paths,
// - designate the new 1cells, associated with the paths,
// - designate the new 2cell

// - Add the triangles belonging to the 2cell, and adding the tetrahedra 
//   inside the 3cells for each triangle.

// - sequential disc expansion: 
//   - select the first path, and find the associated disc,
//   - expand the disc using the designated 0cell and 1cell associated with the
//     path. The center edges of these discs are the ones already generated with
//     the 2cell.
//   - select the next path, by finding one with at least one of the vertices of 
//     its 0cells formed (if we go sequentially, only the last one will have 
//     both formed).
//   - Keep track of the newly added 3cell edges and 1cell edges. ??


// ------------------------------------------------

// Calculate the energy change. 
// This calculates rate of energy change of the current mesh.
// Trial mesh may require addition of fields for energy related 
// calculations. This should be governed by the vd_sim object.

double vd_edisc::calc_energy() {

  //std::cout << " Preconditioned mesh E calculation" << std::endl;
  vd_entlist_v ent_list(m_act, vert_ctr_act, c_base);
  //ent_list.print();

  double e_tot = 0;

  std::vector<std::vector<apf::MeshEntity*> > es_ent
      (4, std::vector<apf::MeshEntity*>(0));

  es_ent.at(0).clear();

  init_ent_set(&es_ent.at(0), vert_ctr_act);
  
  // Collect the upper adjacencies of these vertices.
  for (int i = 0; i < 3; i++) {
    vd_set_up(m_act, &es_ent[i], &es_ent[i+1]);    
  }

  //std::cout << "Volumetric entities: " << es_ent[3].n << std::endl;
  // Calculate the volumetric and surface energies.
  for (int i = 0; i < es_ent[3].size(); i++) {
    //e_tot = e_tot + (*energy_func_tet)(m_act, es_ent.at(3).at(i));
    e_tot = e_tot + energy_func_tet(m_act, es_ent.at(3).at(i));
  }

  //std::cout << "Surface entities: " << es_ent.at(2).size() << std::endl;
  for (int i = 0; i < es_ent[2].size(); i++) {
    //std::cout << es_ent.at(2).at(i) << "[" << 
    //m_act->getModelType(m_act->toModel(es_ent.at(2).at(i))) << "c" << 
    //m_act->getModelTag(m_act->toModel(es_ent.at(2).at(i))) << "]" 
    //  << std::endl;

    e_tot = e_tot + energy_func_tri(m_act, es_ent.at(2).at(i));
    //e_tot = e_tot + (*energy_func_tri)(m_act, es_ent.at(2).at(i));
    //std::cout << e_tot << " Surf area, "
    //   << vd_area_out(m_act, es_ent.at(2).at(i), 0) << std::endl;

  }
  //return e_tot / ins_length;
  //return e_tot / f_calc.vdparam.dt;
  return e_tot;
}

// Calculate the energy after the lens expansion.
double vd_edisc::calc_energy_lens() {

  //vd_entlist_v ent_list(e_lens.m, &e_lens.vert_ctr_new, e_lens.c_base_curr);
  //ent_list.print();

  double e_tot = 0;

  std::vector<std::vector<apf::MeshEntity*> > es_ent
      (4, std::vector<apf::MeshEntity*>(0));


  es_ent.reserve(e_lens.vert_ctr_new.size()+e_lens.vert_sp_new.size());
  // Collect the new center vertices.
  //std::cout << "Lens calculation" << std::endl;
  //std::cout << "Vertices:" << std::endl;
  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    //e_lens.m->getModelTag(e_lens.m->toModel(es_ent[2].e[i])) << std::endl;
    es_ent.at(0).push_back(e_lens.vert_ctr_new.at(i));
    //std::cout << e_lens.vert_ctr_new[i] << " ["
    //  << e_lens.m->getModelType(e_lens.m->toModel(e_lens.vert_ctr_new[i])) 
    //<<"c"<< e_lens.m->getModelTag(e_lens.m->toModel(e_lens.vert_ctr_new[i])) 
    //<< "], ";
  }
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    //e_lens.m->getModelTag(e_lens.m->toModel(es_ent[2].e[i])) << std::endl;
    es_ent.at(0).push_back(e_lens.vert_sp_new.at(i));
    //std::cout << e_lens.vert_sp_new[i] << " ["
    //<< e_lens.m->getModelType(e_lens.m->toModel(e_lens.vert_sp_new[i])) 
    //<<"c"<< e_lens.m->getModelTag(e_lens.m->toModel(e_lens.vert_sp_new[i])) 
    //<< "], ";
  }
  //std::cout << std::endl;
  
  // Collect the upper adjacencies of these vertices.
  for (int i = 0; i < 3; i++) {
    vd_set_up(e_lens.m, &es_ent[i], &es_ent[i+1]);    
  }

  //std::cout << "Volumetric entities: " << es_ent[3].size() << std::endl;
  // Calculate the volumetric and surface energies.
  for (int i = 0; i < es_ent[3].size(); i++) {
    e_tot = e_tot + energy_func_tet(e_lens.m, es_ent.at(3).at(i));
    //e_tot = e_tot + (*energy_func_tet)(e_lens.m, es_ent.at(3).at(i));
  }

  //std::cout << "Surface entities: " << es_ent[2].size() << std::endl;
  for (int i = 0; i < es_ent[2].size(); i++) {
    //std::cout << es_ent.at(2).at(i) << "[" << 
    //e_lens.m->getModelType(e_lens.m->toModel(es_ent.at(2).at(i))) << "c" << 
    //e_lens.m->getModelTag(e_lens.m->toModel(es_ent.at(2).at(i))) << "]" 
    //  << std::endl;

    e_tot = e_tot + energy_func_tri(e_lens.m, es_ent.at(2).at(i));
    //e_tot = e_tot + (*energy_func_tri)(e_lens.m, es_ent.at(2).at(i));
    //std::cout << e_tot << " Surf area, "
    //   << vd_area_out(e_lens.m, es_ent.at(2).at(i), 0) << std::endl;

  }
  //return e_tot / ins_length  - energy;
  std::cout << "Energy calculation, "
            << " vd_par.dt " << f_calc.vdparam.dt
            << " e_tot " << e_tot
            << " energy " << energy
            << std::endl;

  return e_tot  - energy;
}

// Using the force on each vertex and their velocities, calculate the energy
// dissipation rate at singular configuration.
double vd_edisc::calc_energy_diss_rate_sing() {

  //apf::Vector3 v_avg(0,0,0);
  apf::Vector3 v_res(0,0,0);
  apf::Vector3 v_miss(0,0,0);
  apf::Vector3 v_temp(0,0,0);

  apf::Field* vel_field = e_lens.m->findField("velocity_field");

  std::vector<apf::Vector3> force_ctr(e_lens.vert_ctr_new.size(), 
                                                  apf::Vector3(0,0,0));
  std::vector<apf::Vector3> force_sp(e_lens.vert_sp_new.size(), 
                                                  apf::Vector3(0,0,0));

  std::map<apf::MeshEntity*, apf::Vector3> pos_ori{};
  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, v_temp);
    pos_ori[e_lens.vert_ctr_new.at(i)] = v_temp;
  }
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, v_temp);
    pos_ori[e_lens.vert_sp_new.at(i)] = v_temp;
  }
  if(trial_type == 2) {
    e_lens.m->getPoint(e_lens.v_2c, 0, v_temp);
    pos_ori[e_lens.v_2c] = v_temp;
  }

  double p_diss = 0;
  v_temp = apf::Vector3(0,0,0);
  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_old);
  }
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_old);
  }

  if(trial_type == 2) {
    e_lens.m->setPoint(e_lens.v_2c, 0, pos_old);
  }

  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr_new.at(i), drag_flag);
  }

  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_sp_new.at(i), drag_flag);
  }

  if(trial_type == 2) {
    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      apf::getVector(vel_field, e_lens.vert_sp_new.at(i), 0, v_temp);

      if(e_lens.vert_sp_new.size() > 2 and
         ext_new[e_lens.vert_sp_new.at(i)]) {
        int v1 = e_lens.slices.at(2*i)-1;
        int v2 = e_lens.slices.at(2*i+1)-1;

        apf::Vector3 vel_1;
        apf::getVector(vel_field, e_lens.vert_ctr_new.at(v1), 0, vel_1);
        apf::getVector(vel_field, e_lens.vert_ctr_new.at(v2), 0, v_temp);
        v_temp = (v_temp + vel_1)/2;
      }
      //else
      //  v_temp = v_temp - v_avg;

      //v_temp = v_temp - v_avg;
      apf::setVector(vel_field, e_lens.vert_sp_new.at(i), 0, v_temp);
    }
  }


  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    apf::getVector(vel_field, e_lens.vert_ctr_new.at(i), 0, 
                                                  e_lens.vt.at(i));
  }

  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    apf::getVector(vel_field, e_lens.vert_sp_new.at(i), 0, 
                                                  e_lens.split_vt.at(i));
  }

  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    force_ctr.at(i) = f_calc.vd_calc_force(e_lens.m, e_lens.vert_ctr_new.at(i));
  }

  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    force_sp.at(i) = f_calc.vd_calc_force(e_lens.m, e_lens.vert_sp_new.at(i));
  }

  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    p_diss = p_diss - force_ctr.at(i)*e_lens.vt.at(i);
  }

  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    p_diss = p_diss - force_sp.at(i)*e_lens.split_vt.at(i);
  }
  std::cout << "rad(1,:) = " << 0 
            << "; p_diss(1,:) = " << p_diss << ";" << std::endl;

  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, 
                              pos_ori[e_lens.vert_ctr_new.at(i)]);
  }
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, 
                              pos_ori[e_lens.vert_sp_new.at(i)]);
  }
  if(trial_type == 2) {
    e_lens.m->setPoint(e_lens.v_2c, 0, pos_ori[e_lens.v_2c]);
  }

  return p_diss;
}

// Using the force on each vertex and their velocities, calculate the energy
// dissipation rate at singular configuration.
void vd_edisc::calc_energy_diss_rate() {

  //apf::Vector3 v_avg(0,0,0);
  apf::Vector3 v_res(0,0,0);
  apf::Vector3 v_miss(0,0,0);
  apf::Vector3 v_temp(0,0,0);

  apf::Field* vel_field = e_lens.m->findField("velocity_field");

  std::vector<apf::Vector3> force_ctr(e_lens.vert_ctr_new.size(), 
                                                  apf::Vector3(0,0,0));
  std::vector<apf::Vector3> force_sp(e_lens.vert_sp_new.size(), 
                                                  apf::Vector3(0,0,0));

  std::map<apf::MeshEntity*, apf::Vector3> pos_ori{};
  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, v_temp);
    pos_ori[e_lens.vert_ctr_new.at(i)] = v_temp;
  }
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, v_temp);
    pos_ori[e_lens.vert_sp_new.at(i)] = v_temp;
  }
  if(trial_type == 2) {
    e_lens.m->getPoint(e_lens.v_2c, 0, v_temp);
    pos_ori[e_lens.v_2c] = v_temp;
  }

  int step = 10;
  if(trial_type == 1) {
    w1_exp.at(trial_curr).resize(step);
  }
  else {
    w2_exp.at(trial_curr).at(ng.curr).resize(step);
  }

  std::vector<double> r_ctr(e_lens.vert_ctr_new.size());
  std::vector<double> r_sp(e_lens.vert_sp_new.size());
  double r_max = get_r_max(r_ctr, r_sp);

  for (int st = step - 1; st > -1; st--) {
    double p_diss = 0;
    v_temp = apf::Vector3(0,0,0);
    if(st == 0) {
      for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
        e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos_old);
      }
      for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
        e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos_old);
      }
      if(trial_type == 2) {
        e_lens.m->setPoint(e_lens.v_2c, 0, pos_old);
      }
    }
    else {
      proj_around_ctr(r_max/step*st);
    }

    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_ctr_new.at(i), drag_flag);
    }

    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      f_calc.vd_upd_vel_field(e_lens.m, e_lens.vert_sp_new.at(i), drag_flag);
    }

    if(trial_type == 2) {
      for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
        apf::getVector(vel_field, e_lens.vert_sp_new.at(i), 0, v_temp);

        if(e_lens.vert_sp_new.size() > 2 and
           ext_new[e_lens.vert_sp_new.at(i)]) {
          int v1 = e_lens.slices.at(2*i)-1;
          int v2 = e_lens.slices.at(2*i+1)-1;

          apf::Vector3 vel_1;
          apf::getVector(vel_field, e_lens.vert_ctr_new.at(v1), 0, vel_1);
          apf::getVector(vel_field, e_lens.vert_ctr_new.at(v2), 0, v_temp);
          v_temp = (v_temp + vel_1)/2;
        }
        //else
        //  v_temp = v_temp - v_avg;

        //v_temp = v_temp - v_avg;
        apf::setVector(vel_field, e_lens.vert_sp_new.at(i), 0, v_temp);
      }
    }


    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      apf::getVector(vel_field, e_lens.vert_ctr_new.at(i), 0, 
                                                    e_lens.vt.at(i));
    }

    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      apf::getVector(vel_field, e_lens.vert_sp_new.at(i), 0, 
                                                    e_lens.split_vt.at(i));
    }

    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      force_ctr.at(i) = f_calc.vd_calc_force(e_lens.m, e_lens.vert_ctr_new.at(i));
    }

    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      force_sp.at(i) = f_calc.vd_calc_force(e_lens.m, e_lens.vert_sp_new.at(i));
    }

    for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      p_diss = p_diss - force_ctr.at(i)*e_lens.vt.at(i);
    }

    for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
      p_diss = p_diss - force_sp.at(i)*e_lens.split_vt.at(i);
    }
    if(trial_type == 1) {
      w1_exp.at(trial_curr).at(st) = p_diss;
    }
    else {
      w2_exp.at(trial_curr).at(ng.curr).at(st) = p_diss;
    }

  }

  for (int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, 
                              pos_ori[e_lens.vert_ctr_new.at(i)]);
  }
  for (int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, 
                              pos_ori[e_lens.vert_sp_new.at(i)]);
  }
  if(trial_type == 2) {
    e_lens.m->setPoint(e_lens.v_2c, 0, pos_ori[e_lens.v_2c]);
  }

}

// Calculate the energy around a given set of vertices.
double vd_edisc::calc_energy_vertices(apf::Mesh2* m_in,
                                     std::vector<apf::MeshEntity*>* v_in) {

  //vd_entlist_v ent_list(e_lens.m, &e_lens.vert_ctr_new, e_lens.c_base_curr);
  //ent_list.print();

  double e_tot = 0;

  std::vector<std::vector<apf::MeshEntity*> > es_ent
      (4, std::vector<apf::MeshEntity*>(0));

  es_ent.at(0).reserve(v_in->size());

  // Collect the new center vertices.
  //std::cout << "Lens calculation" << std::endl;
  //std::cout << "Vertices:" << std::endl;
  // TODO remnant of Entity_set simply replace this with std copy 
  for (int i = 0; i < v_in->size(); i++) {
    //e_lens.m->getModelTag(e_lens.m->toModel(es_ent[2].e[i])) << std::endl;
    es_ent.at(0).push_back(v_in->at(i));
    //std::cout << e_lens.vert_ctr_new[i] << " ["
    //  << e_lens.m->getModelType(e_lens.m->toModel(e_lens.vert_ctr_new[i])) 
    //<<"c"<< e_lens.m->getModelTag(e_lens.m->toModel(e_lens.vert_ctr_new[i])) 
    //<< "], ";
  }
  //std::cout << std::endl;
  
  // Collect the upper adjacencies of these vertices.
  for (int i = 0; i < 3; i++) {
    vd_set_up(m_in, &es_ent[i], &es_ent[i+1]);    
  }

  //std::cout << "Volumetric entities: " << es_ent[3].n << std::endl;
  // Calculate the volumetric and surface energies.
  for (int i = 0; i < es_ent[3].size(); i++) {
    e_tot = e_tot + energy_func_tet(m_in, es_ent.at(3).at(i));
    //e_tot = e_tot + (*energy_func_tet)(m_in, es_ent.at(3).at(i));
  }

  //std::cout << "Surface entities: " << es_ent[2].n << std::endl;
  for (int i = 0; i < es_ent[2].size(); i++) {
    //std::cout << es_ent[2].e[i] << "[" << 
    //e_lens.m->getModelType(e_lens.m->toModel(es_ent[2].e[i])) << "c" << 
    //e_lens.m->getModelTag(e_lens.m->toModel(es_ent[2].e[i])) << "]" 
    //  << std::endl;

    //e_tot = e_tot + (*energy_func_tri)(m_in, es_ent.at(2).at(i));
    e_tot = e_tot + energy_func_tri(m_in, es_ent.at(2).at(i));
    //std::cout << e_tot << " Surf area, "
    //   << vd_area_out(e_lens.m, es_ent[2].e[i], 0) << std::endl;

  }
  //return e_tot / ins_length  - energy;
  return e_tot / f_calc.vdparam.dt;

  //return e_tot  - energy;
}

void vd_edisc::set_len_sh(double len_in) {
  len_sh = len_in;
}

void vd_edisc::set_rho_rat(double rat_in, double fudge_in) {
  rho_rat = rat_in;
  fudge_factor = fudge_in;
}

double vd_edisc::get_len_sp() {
  return len_sp;
}

void vd_edisc::set_field_calc(const field_calc& FC) {
  f_calc = FC;
  calc_ext = f_calc.get_ext();
  calc_corner = f_calc.get_corner();

  //if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
  //  e_sh_save = *f_calc.get_e_sh();
  //}
  f_calc.set_integ_type(INTEG_TYPE::EULER);
}

// Set the flag for skipping energy dissipation rate calculation. Used in
// applying a predetermined topological change.
void vd_edisc::set_en_skip(bool skip_in) {
  skip_en = skip_in;
}

// Set the simulation parameters.
// TODO f_calc is already passed. Copy from vdparam of f_calc.
void vd_edisc::set_vdpar(vd_param par_in) {
  f_calc.vdparam = par_in;
}

void vd_edisc::set_isotropic(bool iso_in) {
  isotropic = iso_in;
}

void vd_edisc::set_calc_ext(bool fix) {
  calc_ext = fix;
  vd_cd->set_calc_ext(fix);
}

void vd_edisc::set_calc_corner(bool fix) {
  assert(calc_ext);
  calc_corner = fix;
  vd_cd->set_calc_corner(fix);
}

void vd_edisc::set_proj(PROJ_TYPE PROJ) {
  int proj_flag = (int)PROJ;
  assert(PROJ < PROJ_TYPE::END);

  if(PROJ > PROJ_TYPE::FIXED) {
    calc_ext = true;
    vd_cd->set_calc_ext(calc_ext);
  }

  if(PROJ > PROJ_TYPE::FIXED and 
     PROJ < PROJ_TYPE::EXT_SHELL) {
    calc_corner = true;
    vd_cd->set_calc_corner(calc_corner);
  }

  if(PROJ == PROJ_TYPE::EXT_SHELL) {
    calc_corner = true;
    vd_cd->set_calc_corner(calc_corner);
  }

}

// Visual output of the trial mesh. 
void vd_edisc::vtk_mesh() {
  if(!sub_vtk)
    return;
  safe_mkdir("./output");

  std::stringstream ss;
  ss << "./output/trial";

  vd_rem_tag(m_act);
  vd_tag_mesh(m_act);

  if(m_act == m_precond) {
    ss << "0cell" << cell_id << "_precond";
    std::string tmp = ss.str();   
    const char* cstr = tmp.c_str();

    safe_mkdir(cstr);
    if (trial_type == 3)
      ss << "/ins"<< trial_type << "cell_" << trial_curr;
    else if (trial_type == 1 or trial_type == 11 or trial_type == 113)
      ss << "/ins"<< trial_type << "cell_" << trial_curr;
    else 
      ss << "/ins"<< trial_type
         << "cell_" << trial_curr << "_" << ng.curr;
  }
/*
  else if(m_act == m_trial) {
    if (trial_type == 1 or trial_type == 11 or trial_type == 113)
      ss << "0cell" << cell_id << "/ins"<< trial_type << "cell_" << trial_curr;
    else 
      ss << "0cell" << cell_id << "/ins"<< trial_type
         << "cell_" << trial_curr << "_" << ng.curr;
  }
*/
  else {
    //assert(m_act == m_main);
    ss << "0cell" << cell_id;

    std::string tmp = ss.str();   
    const char* cstr = tmp.c_str();

    safe_mkdir(cstr);
    if (trial_type == 3)
      ss << "/ins"<< trial_type << "cell_" << trial_curr;
    else if (trial_type == 1 or trial_type == 11 or trial_type == 113)
      ss << "/ins"<< trial_type << "cell_" << trial_curr;
    else 
      ss << "/ins"<< trial_type
         << "cell_" << trial_curr << "_" << ng.curr;
  }
  std::string tmp = ss.str();   
  const char* cstr = tmp.c_str();
  apf::writeVtkFiles(cstr, m_act);
}

void vd_edisc::vtk_trial() {
  if(!sub_vtk)
    return;

  std::stringstream ss;
  ss << "./output/trial";

  vd_rem_tag(m_trial);
  vd_tag_mesh(m_trial);

  if (trial_type == 1 or trial_type == 11 or trial_type == 113)
    ss << "0cell" << cell_id << "/ins"<< trial_type << "cell_" << trial_curr;
  else if (trial_type == 2 or trial_type == 22 or trial_type == 223)
    ss << "0cell" << cell_id << "/ins"<< trial_type
       << "cell_" << trial_curr << "_" << ng.curr;
  else
    ss << "0cell" << cell_id << "/ins"<< trial_type
       << "cell_" << trial_curr;

  std::string tmp = ss.str();   
  const char* cstr = tmp.c_str();
  apf::writeVtkFiles(cstr, m_trial);
}

void vd_edisc::vtk_precond() {
  if(!sub_vtk)
    return;

  std::stringstream ss;
  ss << "./output/trial";

  vd_rem_tag(m_precond);
  vd_tag_mesh(m_precond);

  ss << "0cell" << cell_id << "_precond";

  std::string tmp = ss.str();   
  const char* cstr = tmp.c_str();
  apf::writeVtkFiles(cstr, m_precond);
}

void vd_edisc::set_sub_vtk_flag(bool vtk_flag) {
  sub_vtk = vtk_flag;
}

void vd_edisc::set_verbose_flag(bool flag_in) {
  verb_flag = flag_in;
}

// Given a list of 2shells, find the joint shell.
shell vd_edisc::get_2shell_joint(std::vector<int>* shell_2_ids) {

  ext_shell* e_sh = f_calc.get_e_sh();

  std::sort(shell_2_ids->begin(), shell_2_ids->end());
  std::vector<int>::iterator it;
  it = std::unique (shell_2_ids->begin(), shell_2_ids->end());
  shell_2_ids->resize(std::distance(shell_2_ids->begin(),it));
  // Belongs to 0shell.
  if(shell_2_ids->size() == 3) {
    ent_conn* e0 = new ent_conn();
    ent_conn* e1 = new ent_conn();
    ent_conn* e2 = new ent_conn();
    e_sh->sh_base.get_conn_dim_gmi(0, 2, shell_2_ids->at(0), e0);
    e_sh->sh_base.get_conn_dim_gmi(0, 2, shell_2_ids->at(1), e1);
    e0->chk_intsct(e1, e2);
    e_sh->sh_base.get_conn_dim_gmi(0, 2, shell_2_ids->at(2), e1);
    e1->chk_intsct(e2, e0);
    assert(e0->conn.size() == 1);
    int sh_id = e0->conn.at(0);
    delete e0;
    delete e1;
    delete e2;

    return shell(0, sh_id);
  }
  // Belongs to 1shell.
  else if(shell_2_ids->size() == 2) {
    ent_conn* e0 = new ent_conn();
    ent_conn* e1 = new ent_conn();
    ent_conn* e2 = new ent_conn();
    e_sh->sh_base.get_conn_dim_gmi(1, 2, shell_2_ids->at(0), e0);
    e_sh->sh_base.get_conn_dim_gmi(1, 2, shell_2_ids->at(1), e1);
    e0->chk_intsct(e1, e2);
    assert(e2->conn.size() == 1);
    int sh_id = e2->conn.at(0);
    delete e0;
    delete e1;
    delete e2;

    return shell(1, sh_id);
  }
  // Belongs to 2shell.
  else {
    assert(shell_2_ids->size() == 1);
    return shell(2, shell_2_ids->at(0));
  }
}

std::vector<int> vd_edisc::get_shell_id_list(int dim, std::vector<int>* tag) {

  ext_shell* e_sh = f_calc.get_e_sh();
  std::vector<int> sh_list(0);
  shell sh_temp;
  sh_list.reserve(tag->size());

  std::map<int, bool> sh_chk{};

  for(int i = 0; i < tag->size(); i++) {
    if(e_lens.c_base_curr->get_cell_ext_gmi(dim, tag->at(i)) ) {
      assert(e_sh->chk_shell(dim, tag->at(i) - 1));
      sh_temp = e_sh->get_shell(dim, tag->at(i) - 1);
      if(!sh_chk[sh_temp.id]) {
        sh_chk[sh_temp.id] = true;
        sh_list.push_back(sh_temp.id);
      }
    }
  }
  std::sort(sh_list.begin(), sh_list.end()); 
  return sh_list;
}

void vd_edisc::assign_01sh() {

  ext_shell* e_sh = f_calc.get_e_sh();
  ent_conn* e_up = new ent_conn();
  std::vector<int> sh2_adj(0);
  // Determine and assign the shells:
  // 1 strata
  for(int j = 0; j < e_lens.new_cell1_id.size(); j++) {
    sh2_adj.clear();
    e_lens.c_base_curr->get_conn_dim_gmi(2, 1, 
                              e_lens.new_cell1_id.at(j), e_up);
    sh2_adj.reserve(e_up->conn.size());
    std::cout << "Available shells:";
    for(int k = 0; k < e_up->conn.size(); k++) {
      //if(e_lens.c_base_curr->get_cell_ext_gmi(2, e_up->conn.at(k)) ) {
        //assert(e_sh->chk_shell(2, e_up->conn.at(k) - 1) );
      if(e_sh->chk_shell(2, e_up->conn.at(k) - 1) ) {
        int id_curr = e_sh->get_shell(2, e_up->conn.at(k) - 1).id;
        sh2_adj.push_back(id_curr);
        std::cout << " 2sh" << id_curr;
      }
    }
    std::cout << std::endl;

    // External, either 1shell or 2shell.
    if(sh2_adj.size() > 0) {
      shell sh_curr = get_2shell_joint(&sh2_adj);
      e_sh->set_shell(1, e_lens.new_cell1_id.at(j) - 1, sh_curr);
      std::cout << "\t1c" << e_lens.new_cell1_id.at(j) << " "
                << sh_curr.dim << "sh" << sh_curr.id << " "
                << " wins!"
                << std::endl;
    }
    // Internal, no shell.
    else {
      e_sh->set_shell(1, e_lens.new_cell1_id.at(j) - 1, sh_old, false);
      std::cout << "\t1c" << e_lens.new_cell1_id.at(j) << " "
                << " Interior!"
                << std::endl;
    }
  }

  for(int j = 0; j < e_lens.new_cell0_id.size(); j++) {
    e_lens.c_base_curr->get_conn_dim_gmi(2, 0, 
                              e_lens.new_cell0_id.at(j), e_up);
    sh2_adj.clear();
    sh2_adj.reserve(e_up->conn.size());
    for(int k = 0; k < e_up->conn.size(); k++) {
      //if(e_lens.c_base_curr->get_cell_ext_gmi(2, e_up->conn.at(k)) ) {
        //assert(e_sh->chk_shell(2, e_up->conn.at(k) - 1) );
      if(e_sh->chk_shell(2, e_up->conn.at(k) - 1) ) {
        sh2_adj.push_back(e_sh->get_shell(2, e_up->conn.at(k) - 1).id);
      }
    }
    // External, either 1shell or 2shell.
    if(sh2_adj.size() > 0) {
      shell sh_curr = get_2shell_joint(&sh2_adj);
      e_sh->set_shell(0, e_lens.new_cell0_id.at(j) - 1, sh_curr);
      std::cout << "\t0c" << e_lens.new_cell0_id.at(j) << " "
                << sh_curr.dim << "sh" << sh_curr.id
                << " wins!"
                << std::endl;
    }
    // Internal, no shell.
    else {
      e_sh->set_shell(0, e_lens.new_cell0_id.at(j) - 1, sh_old, false);
      std::cout << "\t0c" << e_lens.new_cell0_id.at(j) << " "
                << " Interior!"
                << std::endl;
    }
  }
  delete e_up;
}

void vd_edisc::assign_2sh(int sh_2id) {
  ext_shell* e_sh = f_calc.get_e_sh();

  shell sh_curr = shell(2, sh_2id);
  e_sh->set_shell(2, e_lens.new_cell2_id - 1, sh_curr);
  std::cout << "2c" << e_lens.new_cell2_id << " "
          << sh_curr.dim << "sh" << sh_curr.id
          << " trial..."
          << std::endl;
  assign_01sh();
}

void vd_edisc::reset_ext() {
  if(ext_0cell) {
    e_lens.c_base_curr->set_ext_gmi(0, cell_id, true);
    if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
      if(sh_old.dim == 0 or sh_old.dim == 1) {
        e_lens.c_base_curr->set_cor_ext_0c_gmi(cell_id, true);
      }
    }
  }
}

// Based on the shell memberships, update the exterior labels of cells.
void vd_edisc::update_ext() {
  if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {

    ext_shell* e_sh = f_calc.get_e_sh();
    // Cell exterior and corner status update.
    for(int j = 0; j < e_lens.new_cell0_id.size(); j++) {
      if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(j) - 1)) {
        shell sh_curr = e_sh->get_shell(0, e_lens.new_cell0_id.at(j) - 1);
        if(sh_curr.dim == 0)
          e_lens.c_base_curr->set_cor_ext_0c_gmi(e_lens.new_cell0_id.at(j), 
                                                                    true);
        e_lens.c_base_curr->set_ext_gmi(0, e_lens.new_cell0_id.at(j), true);
      }
      else {
        e_lens.c_base_curr->set_cor_ext_0c_gmi(e_lens.new_cell0_id.at(j), 
                                                                    false);
        e_lens.c_base_curr->set_ext_gmi(0, e_lens.new_cell0_id.at(j), false);

      }
    }
    for(int j = 0; j < e_lens.new_cell1_id.size(); j++) {
      if(e_sh->chk_shell(1, e_lens.new_cell1_id.at(j) - 1)) {
        shell sh_curr = e_sh->get_shell(1, e_lens.new_cell1_id.at(j) - 1);
        if(sh_curr.dim == 1)
          e_lens.c_base_curr->set_cor_ext_gmi(e_lens.new_cell1_id.at(j), 
                                                                    true);
        e_lens.c_base_curr->set_ext_gmi(1, e_lens.new_cell1_id.at(j), true);
      }
      else {
        e_lens.c_base_curr->set_cor_ext_gmi(e_lens.new_cell1_id.at(j), 
                                                                    false);
        e_lens.c_base_curr->set_ext_gmi(1, e_lens.new_cell1_id.at(j), false);
      }
    }

    if(trial_type == 2 and e_sh->chk_shell(2, e_lens.new_cell2_id - 1))
      e_lens.c_base_curr->set_ext_gmi(2, e_lens.new_cell2_id, true);
    else
      e_lens.c_base_curr->set_ext_gmi(2, e_lens.new_cell2_id, false);
  }
}

bool vd_edisc::chk_bound_spur() {
  assert(trial_type == 2);
  int spur_nbr = 0;
  for(int i = 0; i < e_lens.new_cell1_id.size(); i++) {
    if(e_lens.c_base_curr->chk_spur(1, e_lens.new_cell1_id.at(i) - 1))
      spur_nbr = spur_nbr + 1;
  }
  ext_shell* e_sh = f_calc.get_e_sh();

  bool multiple = false;
  std::map<int, bool> sh1_map{};
  for(int i = 0; i < e_lens.new_cell1_id.size(); i++) {
    if(e_sh->chk_shell(1, e_lens.new_cell1_id.at(i) - 1)) {
      shell sh_curr = e_sh->get_shell(1, e_lens.new_cell1_id.at(i) - 1);
      if(sh_curr.dim == 1) {
        if(sh1_map[sh_curr.id])
          multiple = true;
        else
          sh1_map[sh_curr.id] = true;
      }
    }
  }

  return (spur_nbr == 1 or multiple);
}


void vd_edisc::write_diss_csv(const char* outfile) {
  std::ofstream ofs;
  ofs.open(outfile, std::ofstream::out | std::ofstream::trunc);

  csvfile csvdiss(outfile);
  for(int i = 0; i < w1.size(); i++) {
    csvdiss << w1.at(i);
    csvdiss << endrow;
  }
  for(int i = 0; i < w2.size(); i++) {
    for(int j = 0; j < w2.at(i).size(); j++) {
      csvdiss << w2.at(i).at(j);
      csvdiss << endrow;
    }
  }
}

void vd_edisc::write_diss_exp_csv(const char* outfile) {
  std::ofstream ofs;
  ofs.open(outfile, std::ofstream::out | std::ofstream::trunc);

  csvfile csvdiss(outfile);
  for(int i = 0; i < w1_exp.size(); i++) {
    for(int j = 0; j < w1_exp.at(i).size(); j++)
      csvdiss << w1_exp.at(i).at(j);
    csvdiss << endrow;
  }
  for(int i = 0; i < w2_exp.size(); i++) {
    for(int j = 0; j < w2_exp.at(i).size(); j++) {
      for(int k = 0; k < w2_exp.at(i).at(j).size(); k++)
        csvdiss << w2_exp.at(i).at(j).at(k);
      csvdiss << endrow;
    }
  }
}

void vd_edisc::write_ei_csv(const char* outfile) {
  std::ofstream ofs;
  ofs.open(outfile, std::ofstream::out | std::ofstream::trunc);

  csvfile csvdiss(outfile);
  for(int i = 0; i < e1.size(); i++) {
    csvdiss << e1.at(i);
    csvdiss << endrow;
  }
  for(int i = 0; i < e2.size(); i++) {
    for(int j = 0; j < e2.at(i).size(); j++) {
      csvdiss << e2.at(i).at(j);
      csvdiss << endrow;
    }
  }
}

// The updated shell assignment. Considers unused shells and dissipation rates
// using different shells when determining shell assignment.
void vd_edisc::determine_shell_gmi() {
  ext_shell* e_sh = f_calc.get_e_sh();
  if(trial_type == 1) {
    // Exterior shell determination of the new strata.
    circuit* circ_tup;

    circ_tup = get_circ(trial_curr);

    if(ext_0cell) {
      assign_01sh();

    }
    else {
      e_sh->set_shell(1, e_lens.new_cell1_id.at(0) - 1, sh_old, false);
      e_sh->set_shell(0, e_lens.new_cell0_id.at(0) - 1, sh_old, false);
      e_sh->set_shell(0, e_lens.new_cell0_id.at(1) - 1, sh_old, false);
      std::cout << " All interior!" << std::endl;
    }
    update_ext();
  }
  else {
    assert(trial_type == 2);
    // Do all shell assignments and compare energies.
    if(ext_0cell) {
      bool ext_3c = (ng.cp.first == -1 or ng.cp.second == -1);
      ent_conn* e_up = new ent_conn();
      std::vector<int> sh2_adj(0);
      if(ext_3c) {
        if(shell2_list.size() == 0) {
          std::cout << "ins_flag: No 2shell available." 
                    << std::endl;
          ins_flag = false;
        }
        else if(!chk_bound_spur() and shell2_list.size() == 1) {
          // Assign the min shell.
          assign_2sh(shell2_list.at(0));
        }
        else {
          // The id of the max magnitude dissipation rate shell assignment.
          int sh2_id = -1;
          double diss_max = 0;

          for(int i = 0; i < shell2_list.size(); i++) {
            assign_2sh(shell2_list.at(i));
            update_ext();
            // Check if only a single bounding 1-cell is spurious. If so,
            // the insertion is spurious and should be skipped.
            if(chk_bound_spur()) {
              std::cout << "2sh" << shell2_list.at(i) << " is migration only."
                        << std::endl;
            }
            else {
              double diss_curr = calc_energy_diss_rate_sing();
              std::cout << "2sh" << shell2_list.at(i) << " p_diss = " 
                        << diss_curr << std::endl;
              if(diss_curr < diss_max) {
                diss_max = diss_curr;
                sh2_id = i;
              }
            }
            // Reset the shells:
            for(int j = 0; j < e_lens.new_cell0_id.size(); j++)
              e_sh->set_shell(0, e_lens.new_cell0_id.at(j) - 1, sh_old, false);
            e_sh->set_shell(0, cell_id - 1, sh_old);
            for(int j = 0; j < e_lens.new_cell1_id.size(); j++)
              e_sh->set_shell(1, e_lens.new_cell1_id.at(j) - 1, sh_old, false);
            e_sh->set_shell(2, e_lens.new_cell2_id - 1, sh_old, false);
          }
          // Assign the min shell.
          if(sh2_id == -1) {
            std::cout << "ins_flag: No 2shell selected." 
                      << std::endl;
            ins_flag = false;
          }
          else
            assign_2sh(shell2_list.at(sh2_id));
        }
      }
      // The 2cell is not exterior but some 0-1 cells are. There can be only one
      // assignment.
      else {
        e_sh->set_shell(2, e_lens.new_cell2_id - 1, sh_old, false);
        std::cout << "\t2c" << e_lens.new_cell2_id << " "
                  << " Interior!"
                  << std::endl;

        assign_01sh();
      }

      delete e_up;
    }
    // No exterior, everything interior.
    else {
      for(int i = 0; i < e_lens.new_cell1_id.size(); i++) {
        e_sh->set_shell(1, e_lens.new_cell1_id.at(i) - 1, sh_old, false);
      }
      for(int i = 0; i < e_lens.new_cell0_id.size(); i++) {
        e_sh->set_shell(0, e_lens.new_cell0_id.at(i) - 1, sh_old, false);
      }
      e_sh->set_shell(2, e_lens.new_cell2_id - 1, sh_old, false);
      std::cout << " All interior!" << std::endl;
    }
    update_ext();
  }
}


// The updated shell assignment. Considers unused shells and dissipation rates
// using different shells when determining shell assignment.
void vd_edisc::determine_shell_gmi_wg() {
  reset_ext();
  ext_shell* e_sh = f_calc.get_e_sh();
  if(trial_type == 1) {
    // Exterior shell determination of the new strata.
    circuit* circ_tup;

    circ_tup = get_circ(trial_curr);

    if(ext_0cell) {
      assign_01sh();

    }
    else {
      e_sh->set_shell(1, e_lens.new_cell1_id.at(0) - 1, sh_old, false);
      e_sh->set_shell(0, e_lens.new_cell0_id.at(0) - 1, sh_old, false);
      e_sh->set_shell(0, e_lens.new_cell0_id.at(1) - 1, sh_old, false);
      std::cout << " All interior!" << std::endl;
    }
    update_ext();
  }
  else {
    assert(trial_type == 2);
    // Do all shell assignments and compare energies.
    if(ext_0cell) {
      bool ext_3c = (ng.cp.first == -1 or ng.cp.second == -1);
      ent_conn* e_up = new ent_conn();
      std::vector<int> sh2_adj(0);
      if(ext_3c) {
        if(shell2_list.size() == 0) {
          std::cout << "ins_flag: No 2shell available." 
                    << std::endl;
          ins_flag = false;
        }
        else if(!chk_bound_spur() and shell2_list.size() == 1) {
          // Assign the min shell.
          assign_2sh(shell2_list.at(0));
        }
        else {
          // The id of the max magnitude dissipation rate shell assignment.
          int sh2_id = -1;
          double diss_max = 0;

          std::vector<apf::Vector3> force_ctr(e_lens.vert_ctr_new.size(), 
                                                  apf::Vector3(0,0,0));
          std::vector<apf::Vector3> force_ctr_temp(e_lens.vert_ctr_new.size(), 
                                                  apf::Vector3(0,0,0));
          std::vector<apf::Vector3> vel_ctr_temp(e_lens.vert_ctr_new.size(), 
                                                  apf::Vector3(0,0,0));
          for(int i = 0; i < slice_tris.size(); i++) {
            force_ctr.at(i) = f_calc.vd_calc_force_tri(e_lens.m, e_lens.vert_ctr, 
                                                        &slice_tris.at(i), false);
          }

          for(int i = 0; i < shell2_list.size(); i++) {
            update_ext();
            assign_2sh(shell2_list.at(i));

            if(f_calc.get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
              ext_shell* e_sh = f_calc.get_e_sh();
              for(int i = 0; i < slice_tris.size(); i++) {
                if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(i) - 1) ) {
                  vel_ctr_temp.at(i) = 
                    e_sh->find_para_dir(0, e_lens.new_cell0_id.at(i) - 1, 
                                                              e_lens.vt.at(i));
                  force_ctr_temp.at(i) = 
                    e_sh->find_para_dir(0, e_lens.new_cell0_id.at(i) - 1, 
                                                              force_ctr.at(i));
                }
              }
            }
            update_ext();
            if(chk_bound_spur()) {
              std::cout << "2sh" << shell2_list.at(i) << " is migration only."
                        << std::endl;
            }
            else {
              double diss_curr = 0;
              for(int i = 0; i < e_lens.new_cell0_id.size(); i++) {
                diss_curr = diss_curr - vel_ctr_temp.at(i)*force_ctr_temp.at(i);
              }

              std::cout << "2sh" << shell2_list.at(i) << " p_diss = " 
                        << diss_curr << std::endl;
              if(diss_curr < diss_max) {
                diss_max = diss_curr;
                sh2_id = i;
              }
            }

            // Reset the shells:
            for(int j = 0; j < e_lens.new_cell0_id.size(); j++)
              e_sh->set_shell(0, e_lens.new_cell0_id.at(j) - 1, sh_old, false);
            e_sh->set_shell(0, cell_id - 1, sh_old);
            for(int j = 0; j < e_lens.new_cell1_id.size(); j++)
              e_sh->set_shell(1, e_lens.new_cell1_id.at(j) - 1, sh_old, false);
            e_sh->set_shell(2, e_lens.new_cell2_id - 1, sh_old, false);
            reset_ext();
          }
          // Assign the min shell.

          if(sh2_id == -1) {
            std::cout << "ins_flag: No 2shell selected." 
                      << std::endl;
            ins_flag = false;
          }
          else
            assign_2sh(shell2_list.at(sh2_id));
        }
      }
      // The 2cell is not exterior but some 0-1 cells are. There can be only one
      // assignment.
      else {
        e_sh->set_shell(2, e_lens.new_cell2_id - 1, sh_old, false);
        std::cout << "\t2c" << e_lens.new_cell2_id << " "
                  << " Interior!"
                  << std::endl;

        assign_01sh();
      }

      delete e_up;
    }
    // No exterior, everything interior.
    else {
      for(int i = 0; i < e_lens.new_cell1_id.size(); i++) {
        e_sh->set_shell(1, e_lens.new_cell1_id.at(i) - 1, sh_old, false);
      }
      for(int i = 0; i < e_lens.new_cell0_id.size(); i++) {
        e_sh->set_shell(0, e_lens.new_cell0_id.at(i) - 1, sh_old, false);
      }
      e_sh->set_shell(2, e_lens.new_cell2_id - 1, sh_old, false);
      std::cout << " All interior!" << std::endl;
    }
    update_ext();
  }
}
/*
void vd_edisc::update_shell_pos() {

  ext_shell* e_sh = f_calc.get_e_sh();
  shell sh_temp;
  shell sh_temp2;

  std::map<apf::MeshEntity*, bool> sk_map{};
  sk_map = calc_skip();
  bool sp_flag = false;
  bool ctr_flag = false;
  int v_id = -1;
  for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    if(sk_map[e_lens.vert_sp_new.at(i)]) {
      sp_flag = true;
      v_id = i;
      i = e_lens.vert_sp_new.size();
    }
  }
  if(!sp_flag) {
    for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      if(sk_map[e_lens.vert_ctr_new.at(i)]) {
        ctr_flag = true;
        v_id = i;
        i = e_lens.vert_ctr_new.size();
      }
    }
  }

  // Before projecting, shift the vertices such that the skipped vertex is on 
  // pos_old.
  apf::Vector3 v_shift;
  apf::Vector3 pos1;
  apf::Vector3 pos2;
  //assert(sp_flag or ctr_flag);
  if(ctr_flag) {
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(v_id), 0, v_shift);
  }
  else if(sp_flag) {
    e_lens.m->getPoint(e_lens.vert_sp_new.at(v_id), 0, v_shift);
  }
  else {
    v_shift = vd_get_center(e_lens.m, &e_lens.vert_ctr_new);
  }
  v_shift = pos_old - v_shift;

  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, pos1);
    e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, pos1+v_shift);
  }

  for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {
    e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, pos1);
    e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, pos1+v_shift);
  }

  if(trial_type == 2) {
    e_lens.m->getPoint(e_lens.v_2c, 0, pos1);
    e_lens.m->setPoint(e_lens.v_2c, 0, pos1+v_shift);
  }

  std::vector<apf::Vector3> shifts_ctr(0, apf::Vector3(0,0,0));
  shifts_ctr.resize(e_lens.vert_ctr_new.size());
  double z[3] = {0,0,0};
  for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
    shifts_ctr.at(i).fromArray(z);
  }

  if(trial_type == 2) {
    for(int i = 0; i < e_lens.vert_ctr_new.size(); i++) {
      if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(i) - 1)) {
        sh_temp = e_sh->get_shell(0, e_lens.new_cell0_id.at(i) - 1);
        if(sh_temp.dim == 1) {
          apf::Up up;
          e_lens.m->getUp(e_lens.vert_ctr_new.at(i), up);

          bool found = false;
          for(int j = 0; j < up.n; j++) {
            apf::ModelEntity* mdl;
            mdl = e_lens.m->toModel(up.e[j]);
            int c_dim = e_lens.m->getModelType(mdl);
            int c_id = e_lens.m->getModelTag(mdl);
            if(e_sh->chk_shell(c_dim, c_id - 1)) {
              sh_temp2 = e_sh->get_shell(c_dim, c_id - 1);
              if(sh_temp.dim == sh_temp2.dim and 
                 sh_temp.id == sh_temp2.id) {

                apf::MeshEntity* v_other = getEdgeVertOppositeVert(e_lens.m, 
                                  up.e[j], e_lens.vert_ctr_new.at(i));
                e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, pos1);
                e_lens.m->getPoint(v_other, 0, v_shift);
                shifts_ctr.at(i) = pos1*(-0.01)+v_shift*0.01;
                e_lens.m->setPoint(e_lens.vert_ctr_new.at(i), 0, 
                                                    pos1*0.99+v_shift*0.01);

                found = true;
                j = up.n;
              }
            }
          }
        }
      }
    }
  }

  for(int dim = 0; dim < 3; dim++) {
    for(int i = 0; i < e_lens.new_cell0_id.size(); i++) {
      if(e_sh->chk_shell(0, e_lens.new_cell0_id.at(i) - 1)) {
        sh_temp = e_sh->get_shell(0, e_lens.new_cell0_id.at(i) - 1);
        if(sh_temp.dim == dim) {
          apf::Vector3 disp_curr;
          apf::Vector3 v_pos_old;
          e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, v_pos_old);
          f_calc.corr_pos(e_lens.m, e_lens.vert_ctr_new.at(i));
          e_lens.m->getPoint(e_lens.vert_ctr_new.at(i), 0, disp_curr);
          shifts_ctr.at(i) = shifts_ctr.at(i) + disp_curr - v_pos_old;
        }
      }
    }
  }

  for(int i = 0; i < e_lens.vert_sp_new.size(); i++) {

    int v1 = e_lens.slices.at(2*i)-1;
    int v2 = e_lens.slices.at(2*i+1)-1;
    apf::Vector3 pos_sp;
    e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, pos_sp);
    e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, 
                    pos_sp + shifts_ctr.at(v1)/2 + shifts_ctr.at(v2)/2);

    if(e_sh->chk_shell(1, e_lens.new_cell1_id.at(i) - 1)) {
      sh_temp = e_sh->get_shell(1, e_lens.new_cell1_id.at(i) - 1);

      //apf::Vector3 pos1;
      //apf::Vector3 pos2;
      //apf::Vector3 pos_sp;
      //apf::Vector3 pos_avg;
      //e_lens.m->getPoint(e_lens.vert_ctr_new.at(v1), 0, pos1);
      //e_lens.m->getPoint(e_lens.vert_ctr_new.at(v2), 0, pos2);
      //e_lens.m->getPoint(e_lens.vert_sp_new.at(i), 0, pos_sp);
      //pos_avg = (pos1 + pos2)/2;
      //pos2 = norm_0(pos2 - pos1);
      //pos_avg = pos_avg - pos_sp;
      //e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, 
      //                pos_sp + shifts_ctr.at(v1)/2 + shifts_ctr.at(v2)/2);

      //e_lens.m->setPoint(e_lens.vert_sp_new.at(i), 0, 
      //                                  pos_sp + pos2*(pos_avg*pos2));
      f_calc.corr_pos(e_lens.m, e_lens.vert_sp_new.at(i));
    }
  }

  if(trial_type == 2) {
    e_lens.pos = (vd_get_center(e_lens.m, &e_lens.vert_ctr_new) + 
                  vd_get_center(e_lens.m, &e_lens.vert_sp_new))/2;
    //}
    e_lens.m->setPoint(e_lens.v_2c, 0, e_lens.pos);
  }
  // TODO Think of a fix in the case of a ngon insertion with 4 or more
  // 1cells, where the shell projection can cause triangles to collapse
  // if prior to the projection some of the vertices lie outside the 
  // simulating cell.
}
*/
