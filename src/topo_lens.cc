#include <vector>
#include <algorithm>    
#include <cstring>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_extinfo.h"
#include "topo_geom.h" // Also contains pi.

#include "topo_topo.h"
#include "topo_geom.h"

#include "topo_ma.h"

#include "topo_lens.h"

/* Topology manipulator:  */
#include "topo_manip.h"

// Lookups for entity ordering, given the order of the collapsing edge.
// The elements are considered in a ringlike fashion, so only a couple 
// of entities of each kind are processed.
// Depending on the vertex to be removed, the order of the indices might
// change. Assuming the first second one is to be removed:
int lookup_vert_col [6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
// Only check the ccw face, with the kept vertex on the top. First 2 ccw, last 
// 2 cw:
//int lookup_edge_mer [6][4] = {{2,1,3,4},{0,2,4,5},{3,5,0,1},{0,4,2,5},{1,5,0,3},{2,3,1,4}};
int lookup_edge_mer [6][4] = {{2,1,4,3},{0,2,5,4},{3,5,1,0},{0,4,5,2},{1,5,3,0},{2,3,4,1}};

int lookup_edge_remain [6] = {5,3,4,1,2,0};
int lookup_surf_mer [6][2] = {{3,2},{1,3},{1,2},{0,2},{0,3},{0,1}};
int lookup_surf_col [6][2] = {{0,1},{0,2},{3,0},{1,3},{2,1},{3,2}};
// int lookup_vert_out [6][2] = {{2,3},{0,3},{3,1},{1,2},{2,0},{0,1}};

// enum for the flag keeping.
//enum {split_6, split_5, split_4, split_3, split_2, split_1, split_0,
//keep_0, keep_1, keep_2, keep_3, keep_4, keep_5, keep_6} merg_flags;
//int merg_index = [6,5,4,3,2,1,0,0,1,2,3,4,5,6];
enum merg_flags {split_1, split_0, keep_0, keep_1, 
                 END};
const char* merg_flag_str[] = {"sp_1", "sp_0", "kp_0", "kp_1"};
int merg_index[] = {1,0,0,1};

// Rebuild the merging entities for the central disc. 
apf::MeshEntity* vd_lens::build_repl(apf::MeshEntity* ent_in, apf::ModelEntity* mdl_in) {

  int i1;
  apf::Downward vert;

  int ent_type = m->getType(ent_in);
  int d = m->typeDimension[ent_type];

  m->getDownward(ent_in, 0, vert);
/*
  std::cout << d << "-ent " << ent_in;
  for(int i = 0; i < d+1; i++) {
    std::cout << " v" << i << " " << vert[i] << " ";
  }
  std::cout << std::endl;
*/
  i1 = findIn(vert, d+1, vert_mer.at(0));
  if(i1 > -1) {
    vert[i1] = vert_ctr;
    //std::cout << vert[i1] << " replaced by " << vert_ctr << std::endl;
    i1 = findIn(vert, d+1, vert_mer.at(1));
    assert(i1 == -1);
  }
  else {
    i1 = findIn(vert, d+1, vert_mer.at(1));
    //std::cout << vert[i1] << " replaced by " << vert_ctr << std::endl;
    assert(i1 > -1);
    vert[i1] = vert_ctr;
  }

  //std::cout << ent_in << " " << m->getModelType(m->toModel(ent_in)) << "c" 
  //          << m->getModelTag(m->toModel(ent_in));

  apf::MeshEntity* ent_n = buildElement(m, mdl_in, ent_type, vert);
  if(d == 3) {
    assert(vd_volume_tet_sign(m, vert));
    //std::cout << "New tet " << ent_n << " " 
    //        << m->getModelType(m->toModel(ent_n))
    //        << "c" << m->getModelTag(m->toModel(ent_n))
    //        << std::endl;
  }
  //std::cout << " replaced by " << ent_n << " "
  //          << m->getModelType(mdl_in) << "c" << m->getModelTag(mdl_in)
  //          << std::endl;

  assert(ent_n != ent_in);
  return ent_n;
}

// Rebuild the splitting entities for the central disc. 
apf::MeshEntity* vd_lens::build_repl(apf::MeshEntity* ent_in, apf::MeshEntity* v_repl) {

  int i1;
  apf::Downward vert;

  int ent_type = m->getType(ent_in);
  int d = m->typeDimension[ent_type];

  apf::ModelEntity* mdl_in = m->toModel(ent_in);

  //std::cout << d << "-ent " << ent_in;

  m->getDownward(ent_in, 0, vert);

  //for(int i = 0; i < d+1; i++) {
  //  std::cout << " v" << i << " " << vert[i] << " ";
  //}
  //std::cout << std::endl;

  i1 = findIn(vert, d+1, v_repl);

  //std::cout << " vert " << vert[i1] << " replaced by " 
  //          << vert_ctr << std::endl;

  assert(i1 > -1);
  vert[i1] = vert_ctr;

  //std::cout << m->getModelType(mdl_in) << "c" << m->getModelTag(mdl_in)
  //          << " replaced by ";

  apf::MeshEntity* ent_n = buildElement(m, mdl_in, ent_type, vert);
  if(d == 3) {
    assert(vd_volume_tet_sign(m, vert));
    //std::cout << "New tet " << ent_n << " " 
    //        << m->getModelType(m->toModel(ent_n))
    //        << "c" << m->getModelTag(m->toModel(ent_n))
    //        << std::endl;
  }
  //std::cout << ent_n << std::endl;

  assert(ent_n != ent_in);
  return ent_n;
}

void vd_lens::l_clear() {
  dummy_clear_stop();

  surf_col.clear();
  elem_col.clear();
  vert_mer.clear();
  type_vertex = -1;
  tag_vertex = -1;

}

// Go over the collapsing edges. Get the entities around them. For merging 
// entities if one is collapsing, the other should as well.
// If they are not collapsing, assign them merging and find the lowest cell dim
// cell one belongs to. If they are same dim, the cell should be the same.
void vd_lens::split_ent() {

  apf::Downward d_edge;
  apf::Downward d_surf;
  apf::Downward vert;

  apf::Downward dv;
  apf::MeshEntity* v_temp;
  apf::MeshEntity* ent_n;

  apf::MeshElement* ee;

  int ent_type;
  int i1;

  double meas; 

  edge_new.first = build_repl(edge_col, vert_mer.at(0));
  edge_new.second = build_repl(edge_col, vert_mer.at(1));

  for(int j = 0; j < elem_col.size(); j++) {
    m->getDownward(elem_col.at(j), 1, d_edge);
    m->getDownward(elem_col.at(j), 2, d_surf);

    int i1 = findIn(d_edge, 6, edge_col);
    int e1 = lookup_edge_mer[i1][0];
    int e2 = lookup_edge_mer[i1][1];
    int e3 = lookup_edge_mer[i1][2];
    int e4 = lookup_edge_mer[i1][3];

    int s1 = lookup_surf_mer[i1][0];
    int s2 = lookup_surf_mer[i1][1];
    int s3 = lookup_surf_col[i1][0];
    int s4 = lookup_surf_col[i1][1];

    apf::ModelEntity* mdl_col1 = m->toModel(d_surf[s3]);
    apf::ModelEntity* mdl_col2 = m->toModel(d_surf[s4]);

    apf::ModelEntity* mdl_col3 = m->toModel(elem_col.at(j));

    build_repl(d_edge[e1], mdl_col1);
    build_repl(d_edge[e3], mdl_col2);
    build_repl(d_surf[s1], mdl_col3);

    build_repl(d_surf[s3], vert_mer.at(0));
    build_repl(d_surf[s3], vert_mer.at(1));

    build_repl(d_surf[s4], vert_mer.at(0));
    build_repl(d_surf[s4], vert_mer.at(1));

    ent_n = build_repl(elem_col.at(j), vert_mer.at(0));
/*
    ee = createMeshElement(m, ent_n);
    meas = measure(ee);
    destroyMeshElement(ee);

    if(meas < 0.) {
    //if(vd_volume_tet(e_lens.m, down) < 0) {
      m->getDownward(ent_n, 0, dv);
      m->destroy(ent_n);
      v_temp = dv[1];
      dv[1] = dv[2];
      dv[2] = v_temp;

      ent_n = buildElement(m, m->toModel(elem_col.at(j)), apf::Mesh::TET, dv);

      ee = createMeshElement(m, ent_n);
      meas = measure(ee);
      destroyMeshElement(ee);
      assert(meas > 0);

      //assert(vd_volume_tet(e_lens.m, down) > 0);
      //vd_print_vert(e_lens.m, e_void);
    }
*/
    ent_n = build_repl(elem_col.at(j), vert_mer.at(1));
/*
    ee = createMeshElement(m, ent_n);
    meas = measure(ee);
    destroyMeshElement(ee);

    if(meas < 0.) {
    //if(vd_volume_tet(e_lens.m, down) < 0) {
      m->getDownward(ent_n, 0, dv);
      m->destroy(ent_n);
      v_temp = dv[1];
      dv[1] = dv[2];
      dv[2] = v_temp;

      ent_n = buildElement(m, m->toModel(elem_col.at(j)), apf::Mesh::TET, dv);

      ee = createMeshElement(m, ent_n);
      meas = measure(ee);
      destroyMeshElement(ee);
      assert(meas > 0);

      //assert(vd_volume_tet(e_lens.m, down) > 0);
      //vd_print_vert(e_lens.m, e_void);
    }

*/
  }
  for(int i = 0; i < elem_col.size(); i++) {
    //printf("Element split %d: %p, ", i, elem_col.at(i));
    m->destroy(elem_col.at(i));
  }
  //printf("Split elements collapsed.\n");

  for(int i = 0; i < surf_col.size(); i++) {
    //printf("Surface split %d: %p, ", i, surf_col.at(i));
    m->destroy(surf_col.at(i));
  }
  //std::cout << std::endl;
  //printf("Split surfaces collapsed.\n");

  m->destroy(edge_col);
  //printf("Edge destroyed.\n");
  m->acceptChanges();
}

void vd_lens::split_lens(apf::MeshEntity* v_cls) {

  if (save_vtk) {
    vd_rem_tag(m);
    vd_tag_mesh(m);

    std::stringstream ss;
    ss << vtk_name << "Edge-" << type_vertex << "c" << tag_vertex << "before";
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    //apf::writeVtkFiles(cstr, m);
    vd_save_vtk_vert(m, &vert_mer, cstr);
  }

  // Create a center vertex. Recreate the entities with the top and bottom vertices replaced by
  // the new vertex. This will double the number of elements. Destroy the collapsing entities.
  //printf("Check vertices: %p, %p against %p.\n", vert_mer[0], vert_mer[1], vert_ctr);

  // If a specific end vertex is given, move the new vertex close to that vertex.
  if(v_cls == NULL)
    m->setPoint(vert_ctr, 0, midpoint);
  else {
    apf::Downward v_d;
    m->getDownward(edge_col, 0, v_d);
    assert(findIn(v_d, 2, v_cls) > -1);

    apf::Vector3 v_cls_pos;
    m->getPoint(v_cls, 0, v_cls_pos);

    v_cls_pos = (v_cls_pos*9 + midpoint)/10;
    m->setPoint(vert_ctr, 0, v_cls_pos);
  }

  split_ent();
  //m->verify();

  vd_rem_tag(m);
  vd_tag_mesh(m);

  if(f_calc != NULL)
    f_calc->vd_att_fields(m, vert_ctr);

  vd_mesh_bad_adj(m);

  if (save_vtk) {
    vd_rem_tag(m);
    vd_tag_mesh(m);

    std::stringstream ss;
    ss << vtk_name << "Edge-" << type_vertex << "c" << tag_vertex;
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    //apf::writeVtkFiles(cstr, m);
    vd_save_vtk_vert(m, &vert_mer, cstr);
  }

  //vd_print_ent(m);

}

void vd_lens::load_edge(apf::MeshEntity* edge_in) {

  assert(m->getType(edge_in) == apf::Mesh::EDGE);

  edge_col = edge_in;

  vd_set_up(m, edge_col, &surf_col);
  vd_set_up(m, &surf_col, &elem_col);

  vd_set_down(m, edge_col, &vert_mer);

  vert_ctr_em = m->toModel(edge_col);
  //if(f_calc != NULL)
  //  vert_ctr = f_calc->create_v(m, vert_ctr_em);
  //else
    vert_ctr = m->createVert(vert_ctr_em);

  //std::cout << "New vert " << vert_ctr
  //          << " " << m->getModelType(vert_ctr_em)
  //          << "c" << m->getModelTag(vert_ctr_em)
  //          << std::endl;

  type_vertex = m->getModelType(vert_ctr_em);
  tag_vertex = m->getModelTag(vert_ctr_em);

  apf::Vector3 v1;
  apf::Vector3 v2;

  m->getPoint(vert_mer.at(0), 0, v1);
  m->getPoint(vert_mer.at(1), 0, v2);

  midpoint = (v1+v2)/2;
}


vd_lens::vd_lens(apf::Mesh2* msh, struct cell_base* c, field_calc* f_in) :
    m(NULL), mdl(NULL), c_base(NULL), f_calc(NULL),
    save_vtk(false), vtk_name("./output/"), smb_name("./tempmesh/temp"),
    edge_new(NULL, NULL),
    type_vertex(0), tag_vertex(0), 
    vert_mer(0), surf_col(0), elem_col(0), 
    midpoint(0, 0, 0) {

  m = msh;
  mdl = m->getModel();
  c_base = c;

  f_calc = f_in;
}

std::pair<int, int> vd_lens::get_cell() {
  return std::make_pair(type_vertex, tag_vertex);
}

void vd_lens::set_field_calc(field_calc* f_in) {
  f_calc = f_in;
}

// Set the output file prefixes.
void vd_lens::set_files(const char* vtkFile, const char* meshFile) {
  vtk_name = vtkFile;
  smb_name = meshFile;
}

std::pair<apf::MeshEntity*, apf::MeshEntity*> vd_lens::get_edge_new() {
  return edge_new;
}

void vd_lens::turn_save_vtk(bool onoff) {
  save_vtk = onoff;
}

// Used after collapsing internal structure within the lens to get the new 
// vertex.
apf::MeshEntity* vd_lens::get_vert_ctr() {
  return vert_ctr;
}

vd_lens::~vd_lens() {
  l_clear();
}

void vd_split_edge(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* edge,  field_calc* f_calc) {
  int ent_type = m->getType(edge);
  int d = m->typeDimension[ent_type];
  assert(d == 1);

  vd_lens* lens_split = new vd_lens(m, c_base, f_calc);

  lens_split->load_edge(edge);
  lens_split->split_lens();
  m->acceptChanges();

  if(f_calc != NULL)
    f_calc->vd_att_fields(m, lens_split->get_vert_ctr());

  delete lens_split;
}

void split_edge_target_len(apf::Mesh2* m, cell_base* cb, field_calc* f_calc, apf::MeshEntity* edge, double len, int max_split) {
  std::vector<apf::MeshEntity*> vert_mer(2);
  std::vector<apf::MeshEntity*> surf_split(0);
  std::vector<apf::MeshEntity*> edge_split(0);
  std::pair<apf::MeshEntity*, apf::MeshEntity*> edge_new{};

  vd_lens lens_split(m, cb, f_calc);

  vd_set_down(m, edge, &vert_mer);
  vd_set_up(m, edge, &surf_split);
  edge_split.resize(surf_split.size());

  // Collect the edges bounding the triangles and joining the splitting edge on the remaining length.
  for(int i = 0; i < surf_split.size(); i++)
    edge_split.at(i) = vd_find_esjv(m, edge, surf_split.at(i), vert_mer.at(1));

  apf::Vector3 p0 = vd_get_pos(m, vert_mer.at(0));
  apf::Vector3 p1 = vd_get_pos(m, vert_mer.at(1));
  p1 = (p1 - p0);
  double len_curr = p1.getLength();
  int iter = 0;
  while(len_curr > len) {
    if(iter == max_split) {
      for(int i = 0; i < surf_split.size(); i++) {
        lens_split.load_edge(edge_split.at(i));
        lens_split.split_lens();
      }
      vd_set_up(m, edge, &surf_split);
      edge_split.resize(surf_split.size());

      // Collect the edges bounding the triangles and joining the splitting edge on the remaining length.
      for(int i = 0; i < surf_split.size(); i++)
        edge_split.at(i) = vd_find_esjv(m, edge, surf_split.at(i), vert_mer.at(1));
      iter = 0;
    }
 
    lens_split.load_edge(edge);
    lens_split.split_lens();
    apf::MeshEntity* v_ctr = lens_split.get_vert_ctr();

    p1 = p0 + p1*(len/len_curr);
    m->setPoint(v_ctr, 0, p1);
    // edge_new.first/second are associated with vert_mer.at(1/0)
    // (replace vert_mer.at(0) to generate edge_new.first, so edge_new.first is bounded by v_ctr_new and vert_mer.at(1))
    // position the new vertex close to vert_mer.at(0) and split edge_new.first
    edge_new = lens_split.get_edge_new();
    edge = edge_new.first;

    v_ctr = vert_mer.at(0);
    vd_set_down(m, edge, &vert_mer);
    assert(vert_mer.at(0) != v_ctr and vert_mer.at(1) != v_ctr);
    p0 = vd_get_pos(m, vert_mer.at(0));
    p1 = vd_get_pos(m, vert_mer.at(1));
    p1 = (p1 - p0);
    len_curr = p1.getLength();
    iter++;
  }
  return;
}

void split_edge_target_nbr(apf::Mesh2* m, cell_base* cb, field_calc* f_calc, apf::MeshEntity* edge, int nbr, int max_split) {
  double len = vd_meas_ent(m, edge)/nbr;
  split_edge_target_len(m, cb, f_calc, edge, len, max_split);
}

void vd_bipy::l_clear() {

  vert_surr.clear();
  vert_surr.resize(3);
  tets.clear();
  tets.reserve(2);

  edges.clear();
  edges.reserve(2);

  type_vertex = -1;
  tag_vertex = -1;

  edge_tri.clear();
  edge_tri.resize(3);
  tri_tri.clear();
  tri_tri.resize(3);
}

// 
void vd_bipy::split_ent() {

  // Given the vertex, return the triangle across.
  int lookup_tet_surf_x [4] = {2, 3, 1, 0};
  // Given the triangle, return the vertex across.
  int lookup_tet_x_surf [4] = {3, 2, 0, 1};

  apf::ModelEntity* m_tri;
  apf::ModelEntity* m_tet;

  apf::Downward d_v;
  apf::Downward d_temp;
  apf::Downward d_t;
  apf::MeshEntity* temp;

  d_v[0] = vert_ctr;

  for(int i = 0; i < vert_surr.size(); i++) {
    d_v[1] = vert_surr.at(i);
    apf::MeshEntity* ent_n = buildElement(m, vert_ctr_em, apf::Mesh::EDGE, d_v);
    std::cout << "New edge " << ent_n << " " 
              << m->getModelType(m->toModel(ent_n))
              << "c" << m->getModelTag(m->toModel(ent_n))
              << std::endl;

  }

  m->getDownward(tri_col, 0, d_v);
  for(int i = 0; i < 3; i++) {
    temp = d_v[i];
    d_v[i] = vert_ctr;
    apf::MeshEntity* ent_n = buildElement(m, vert_ctr_em, apf::Mesh::TRIANGLE,
                                                                           d_v);
    std::cout << "New triangle " << ent_n << " " 
              << m->getModelType(m->toModel(ent_n))
              << "c" << m->getModelTag(m->toModel(ent_n))
              << std::endl;
    d_v[i] = temp;
  }

  for(int i = 0; i < tets.size(); i++) {
    m_tet = m->toModel(tets.at(i));
    assert(m->getModelType(m_tet) == 3);

    m->getDownward(tets.at(i), 0, d_v);
    m->getDownward(tets.at(i), 2, d_t);
    int t1 = findIn(d_t, 4, tri_col);
    assert(t1 > -1);
    d_temp[0] = d_v[lookup_tet_x_surf[t1]];
    d_temp[1] = vert_ctr;
    edges.at(i) = buildElement(m, m_tet, apf::Mesh::EDGE, d_temp);

    std::cout << "New edge " << edges.at(i) << " " 
              << m->getModelType(m->toModel(edges.at(i)))
              << "c" << m->getModelTag(m->toModel(edges.at(i)))
              << std::endl;

    for(int j = 0; j < vert_surr.size(); j++) {
      d_temp[2] = vert_surr.at(j);
      apf::MeshEntity* ent_n = buildElement(m, m_tet, apf::Mesh::TRIANGLE, 
                                                                    d_temp);

      std::cout << "New triangle " << ent_n << " " 
              << m->getModelType(m->toModel(ent_n))
              << "c" << m->getModelTag(m->toModel(ent_n))
              << std::endl;
    }

    for(int j = 0; j < vert_surr.size(); j++) {
      int i1 = findIn(d_v, 4, vert_surr.at(j));
      assert(i1 > -1);

      temp = d_v[i1];

      d_v[i1] = vert_ctr;
      assert(vd_volume_tet_sign(m, d_v));
      apf::MeshEntity* ent_n = buildElement(m, m_tet, apf::Mesh::TET, d_v);
      //std::cout << "New tet " << ent_n << " " 
      //        << m->getModelType(m->toModel(ent_n))
      //        << "c" << m->getModelTag(m->toModel(ent_n))
      //        << std::endl;
      d_v[i1] = temp;
    }
  }

  for(int i = 0; i < tets.size(); i++) {
    //printf("Element split %d: %p, ", i, elem_col.at(i));
    m->destroy(tets.at(i));
  }
  //printf("Split elements collapsed.\n");

  m->destroy(tri_col);

  m->acceptChanges();
}

void vd_bipy::split_bipy(apf::MeshEntity* v_cls) {
  if (save_vtk) {
    vd_rem_tag(m);
    vd_tag_mesh(m);

    std::stringstream ss;
    ss << vtk_name << "Edge-" << type_vertex << "c" << tag_vertex << "before";
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    //apf::writeVtkFiles(cstr, m);
    vd_save_vtk_vert(m, &vert_surr, cstr);
  }

  // Create a center vertex. Recreate the entities with the top and bottom vertices replaced by
  // the new vertex. This will double the number of elements. Destroy the collapsing entities.
  //printf("Check vertices: %p, %p against %p.\n", vert_mer[0], vert_mer[1], vert_ctr);

  // If a specific end vertex is given, move the new vertex close to that vertex.
  if(v_cls == NULL)
    m->setPoint(vert_ctr, 0, midpoint);
  else {
    apf::Downward v_d;
    m->getDownward(tri_col, 0, v_d);
    assert(findIn(v_d, 3, v_cls) > -1);

    apf::Vector3 v_cls_pos(0,0,0);
    m->getPoint(v_cls, 0, v_cls_pos);

    v_cls_pos = (v_cls_pos*9 + midpoint)/10;
    m->setPoint(vert_ctr, 0, v_cls_pos);
  }

  split_ent();

  //m->verify();

  vd_rem_tag(m);
  vd_tag_mesh(m);

  if(f_calc != NULL)
    f_calc->vd_att_fields(m, vert_ctr);

  vd_mesh_bad_adj(m);

  if (save_vtk) {
    vd_rem_tag(m);
    vd_tag_mesh(m);

    std::stringstream ss;
    ss << vtk_name << "Edge-" << type_vertex << "c" << tag_vertex;
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    //apf::writeVtkFiles(cstr, m);
    vd_save_vtk_vert(m, &vert_surr, cstr);
  }

  //vd_print_ent(m);

}

void vd_bipy::load_tri(apf::MeshEntity* tri_in) {

  assert(m->getType(tri_in) == apf::Mesh::TRIANGLE);

  tri_col = tri_in;

  vd_set_up(m, tri_col, &tets);

  vd_set_down(m, tri_col, &edges);
  vd_set_down(m, &edges, &vert_surr);

  edges.resize(tets.size());

  vert_ctr_em = m->toModel(tri_col);
  //if(f_calc != NULL)
  //  vert_ctr = f_calc->create_v(m, vert_ctr_em);
  //else
    vert_ctr = m->createVert(vert_ctr_em);

  //std::cout << "New vert " << vert_ctr
  //          << " " << m->getModelType(vert_ctr_em)
  //          << "c" << m->getModelTag(vert_ctr_em)
  //          << std::endl;

  type_vertex = m->getModelType(vert_ctr_em);
  tag_vertex = m->getModelTag(vert_ctr_em);

  midpoint = vd_get_center(m, &vert_surr);
}


vd_bipy::vd_bipy(apf::Mesh2* msh, struct cell_base* c, field_calc* f_in) :
    m(NULL), mdl(NULL), c_base(NULL), f_calc(NULL),
    save_vtk(false), vtk_name("./output/"), smb_name("./tempmesh/temp"),

    type_vertex(0), tag_vertex(0), 
    vert_surr(3), edge_tri(3), tri_tri(3), tets(1), edges(1),
    midpoint(0, 0, 0) {

  m = msh;
  mdl = m->getModel();
  c_base = c;
  f_calc = f_in;
}

std::pair<int, int> vd_bipy::get_cell() {
  return std::make_pair(type_vertex, tag_vertex);
}

void vd_bipy::set_field_calc(field_calc* f_in) {
  f_calc = f_in;
}

// Set the output file prefixes.
void vd_bipy::set_files(const char* vtkFile, const char* meshFile) {
  vtk_name = vtkFile;
  smb_name = meshFile;
}

void vd_bipy::turn_save_vtk(bool onoff) {
  save_vtk = onoff;
}

// Used after collapsing internal structure within the bipyramid to get the new 
// vertex.
apf::MeshEntity* vd_bipy::get_vert_ctr() {
  return vert_ctr;
}

vd_bipy::~vd_bipy() {
  l_clear();
}

void vd_sp_tet::l_clear() {

  vert_surr.clear();
  vert_surr.resize(4);

  type_vertex = -1;
  tag_vertex = -1;
}

// 
void vd_sp_tet::split_ent() {

  apf::Downward d_v;
  apf::Downward d_e;
  apf::Downward d_temp;
  apf::MeshEntity* temp;

  d_v[0] = vert_ctr;

  for(int i = 0; i < vert_surr.size(); i++) {
    d_v[1] = vert_surr.at(i);
    apf::MeshEntity* ent_n = buildElement(m, vert_ctr_em, apf::Mesh::EDGE, d_v);
    std::cout << "New edge " << ent_n << " " 
              << m->getModelType(m->toModel(ent_n))
              << "c" << m->getModelTag(m->toModel(ent_n))
              << std::endl;
  }

  m->getDownward(tet_col, 1, d_e);
  for(int i = 0; i < 6; i++) {
    m->getDownward(d_e[i], 0, d_v);
    d_v[2] = vert_ctr;
    apf::MeshEntity* ent_n = buildElement(m, vert_ctr_em, apf::Mesh::TRIANGLE,
                                                                           d_v);
    std::cout << "New triangle " << ent_n << " " 
              << m->getModelType(m->toModel(ent_n))
              << "c" << m->getModelTag(m->toModel(ent_n))
              << std::endl;
  }

  m->getDownward(tet_col, 0, d_v);
  for(int i = 0; i < 4; i++) {
    temp = d_v[i];
    d_v[i] = vert_ctr;
    assert(vd_volume_tet_sign(m, d_v));
    apf::MeshEntity* ent_n = buildElement(m, vert_ctr_em, apf::Mesh::TET, d_v);
    //std::cout << "New tet " << ent_n << " " 
    //          << m->getModelType(m->toModel(ent_n))
    //          << "c" << m->getModelTag(m->toModel(ent_n))
    //          << std::endl;

    d_v[i] = temp;
  }

  m->destroy(tet_col);

  m->acceptChanges();
}

void vd_sp_tet::split_tet(apf::MeshEntity* v_cls) {
  if (save_vtk) {
    vd_rem_tag(m);
    vd_tag_mesh(m);

    std::stringstream ss;
    ss << vtk_name << "Edge-" << type_vertex << "c" << tag_vertex << "before";
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    //apf::writeVtkFiles(cstr, m);
    vd_save_vtk_vert(m, &vert_surr, cstr);
  }

  // Create a center vertex. Recreate the entities with the top and bottom vertices replaced by
  // the new vertex. This will double the number of elements. Destroy the collapsing entities.
  //printf("Check vertices: %p, %p against %p.\n", vert_mer[0], vert_mer[1], vert_ctr);

  // If a specific end vertex is given, move the new vertex close to that vertex.
  if(v_cls == NULL)
    m->setPoint(vert_ctr, 0, midpoint);
  else {
    apf::Downward v_d;
    m->getDownward(tet_col, 0, v_d);
    assert(findIn(v_d, 3, v_cls) > -1);

    apf::Vector3 v_cls_pos;
    m->getPoint(v_cls, 0, v_cls_pos);

    v_cls_pos = (v_cls_pos*9 + midpoint)/10;
    m->setPoint(vert_ctr, 0, v_cls_pos);
  }

  split_ent();

  //m->verify();

  vd_rem_tag(m);
  vd_tag_mesh(m);

  if(f_calc != NULL)
    f_calc->vd_att_fields(m, vert_ctr);

  vd_mesh_bad_adj(m);

  if (save_vtk) {
    vd_rem_tag(m);
    vd_tag_mesh(m);

    std::stringstream ss;
    ss << vtk_name << "Edge-" << type_vertex << "c" << tag_vertex;
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    //apf::writeVtkFiles(cstr, m);
    vd_save_vtk_vert(m, &vert_surr, cstr);
  }

  //vd_print_ent(m);

}

void vd_sp_tet::load_tet(apf::MeshEntity* tet_in) {

  assert(m->getType(tet_in) == apf::Mesh::TET);

  tet_col = tet_in;

  std::vector<apf::MeshEntity*> tris(0);
  std::vector<apf::MeshEntity*> edges(0);

  vd_set_down(m, tet_col, &tris);
  vd_set_down(m, &tris, &edges);
  vd_set_down(m, &edges, &vert_surr);

  vert_ctr_em = m->toModel(tet_col);
  //if(f_calc != NULL)
  //  vert_ctr = f_calc->create_v(m, vert_ctr_em);
  //else
    vert_ctr = m->createVert(vert_ctr_em);

  //std::cout << "New vert " << vert_ctr
  //          << " " << m->getModelType(vert_ctr_em)
  //          << "c" << m->getModelTag(vert_ctr_em)
  //          << std::endl;

  type_vertex = m->getModelType(vert_ctr_em);
  tag_vertex = m->getModelTag(vert_ctr_em);

  midpoint = vd_get_center(m, &vert_surr);
}


vd_sp_tet::vd_sp_tet(apf::Mesh2* msh, struct cell_base* c, field_calc* f_in) :
    m(NULL), mdl(NULL), c_base(NULL), f_calc(NULL),
    save_vtk(false), vtk_name("./output/"), smb_name("./tempmesh/temp"),

    type_vertex(0), tag_vertex(0), 
    vert_surr(3),
    midpoint(0, 0, 0) {

  m = msh;
  mdl = m->getModel();
  c_base = c;
  f_calc = f_in;
}

std::pair<int, int> vd_sp_tet::get_cell() {
  return std::make_pair(type_vertex, tag_vertex);
}

void vd_sp_tet::set_field_calc(field_calc* f_in) {
  f_calc = f_in;
}

// Set the output file prefixes.
void vd_sp_tet::set_files(const char* vtkFile, const char* meshFile) {
  vtk_name = vtkFile;
  smb_name = meshFile;
}

void vd_sp_tet::turn_save_vtk(bool onoff) {
  save_vtk = onoff;
}

// Used after collapsing internal structure within the bipyramid to get the new 
// vertex.
apf::MeshEntity* vd_sp_tet::get_vert_ctr() {
  return vert_ctr;
}

vd_sp_tet::~vd_sp_tet() {
  l_clear();
}
