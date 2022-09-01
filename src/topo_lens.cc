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

void vd_split_edge_2c(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* tri, apf::MeshEntity* vert, field_calc* f_calc) {
  apf::ModelEntity* mdl = m->toModel(tri);
  assert(mdl == m->toModel(vert));
  int c_type = m->getModelType(mdl);
  int c_tag = m->getModelTag(mdl);
  assert(c_type == 2);

  std::vector<apf::MeshEntity*> edge_tri(0);
  std::vector<apf::MeshEntity*> edge_vert(0);
  vd_set_down(m, tri, &edge_tri);
  vd_set_up(m, vert, &edge_vert);
  vd_keep_cell(m, &edge_tri, c_tag, c_type);
  vd_keep_cell(m, &edge_vert, c_tag, c_type);

  vd_remove_set(&edge_vert, &edge_tri);

  if(edge_vert.size() == 0) {
    std::cout << "No edges fit the criterion." << std::endl;
    return;
  }

  apf::MeshEntity* edge = NULL;
  double len_max = -1.;
  for(int i = 0; i < edge_vert.size(); i++) {
    double len_temp = vd_meas_ent(m, edge_vert.at(i));
    if(edge == NULL or len_temp > len_max)
      edge = edge_vert.at(i);
  }

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

void vd_split_tet(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* tet,  field_calc* f_calc) {
  int ent_type = m->getType(tet);
  int d = m->typeDimension[ent_type];
  assert(d == 3);

  vd_sp_tet* tet_split = new vd_sp_tet(m, c_base, f_calc);

  tet_split->load_tet(tet);
  tet_split->split_tet();
  m->acceptChanges();

  if(f_calc != NULL)
    f_calc->vd_att_fields(m, tet_split->get_vert_ctr());

  delete tet_split;
}

void vd_split_tet_vert(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* vert,  field_calc* f_calc) {
  int ent_type = m->getType(vert);
  int d = m->typeDimension[ent_type];
  assert(d == 0);

  std::vector<apf::MeshEntity* > edge(0);
  std::vector<apf::MeshEntity* > tri(0);
  std::vector<apf::MeshEntity* > tet(0);
  vd_set_up(m, vert, &edge);
  vd_set_up(m, &edge, &tri);
  vd_set_up(m, &tri, &tet);

  vd_sp_tet* tet_split = new vd_sp_tet(m, c_base, f_calc);
  for(int i = 0; i < tet.size(); i++) {
    tet_split->load_tet(tet.at(i));
    tet_split->split_tet();
    m->acceptChanges();
    if(f_calc != NULL)
      f_calc->vd_att_fields(m, tet_split->get_vert_ctr());
  }

  delete tet_split;
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
    //assert(vd_volume_tet_sign(m, d_v));
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

///////////////////////////////////////
// vd_lens_col
///////////////////////////////////////

// Rebuild the splitting entities for the central disc. 
apf::MeshEntity* vd_lens_col::build_repl(std::vector<apf::MeshEntity*> &vert_in, apf::ModelEntity* mdl_in) {

  apf::Downward vert;

  for(int i = 0; i < vert_in.size(); i++)
    vert[i] = vert_in.at(i);
  apf::MeshEntity* ent_n(NULL);
  int d = vert_in.size() - 1;
  if(d == 1)
    ent_n = buildElement(m, mdl_in, apf::Mesh::EDGE, vert);
  else if(d == 2)
    ent_n = buildElement(m, mdl_in, apf::Mesh::TRIANGLE, vert);
  else {
    assert(d == 3);
    ent_n = buildElement(m, mdl_in, apf::Mesh::TET, vert);
    assert(vd_volume_tet_sign(m, vert));
  }
  return ent_n;
}

// Rebuild the splitting entities for the central disc. 
apf::MeshEntity* vd_lens_col::build_repl(apf::MeshEntity* ent_in, apf::ModelEntity* mdl_in, apf::MeshEntity* v_repl) {

  int i1;
  apf::Downward vert;

  int ent_type = m->getType(ent_in);
  int d = m->typeDimension[ent_type];

  //apf::ModelEntity* mdl_in = m->toModel(ent_in);

  m->getDownward(ent_in, 0, vert);
  i1 = findIn(vert, d+1, v_repl);

  if(i1 == -1)
    return NULL;
  vert[i1] = vert_ctr;

  apf::MeshEntity* ent_n = buildElement(m, mdl_in, ent_type, vert);
  if(d == 3) {
    assert(vd_volume_tet_sign(m, vert));
  }

  assert(ent_n != ent_in);
  return ent_n;
}

void vd_lens_col::l_clear() {

  clear();
  type_vertex = -1;
  tag_vertex = -1;

}

bool vd_lens_col::chk_ents(apf::MeshEntity* ec, apf::MeshEntity* e1, apf::MeshEntity* e2) {
  apf::ModelEntity* mdl1 = m->toModel(e1);
  apf::ModelEntity* mdl2 = m->toModel(e2);
  apf::ModelEntity* mdlc = m->toModel(ec);

  int d_1 = m->getModelType(mdl1);
  int id_1 = m->getModelTag(mdl1);
  int d_2 = m->getModelType(mdl2);
  int id_2 = m->getModelTag(mdl2);
  int d_c = m->getModelType(mdlc);
  int id_c = m->getModelTag(mdlc);

  bool same_1c = (mdl1 == mdlc);
  bool same_2c = (mdl2 == mdlc);
  bool bound_1c = c_base->chk_conn_d_gmi(d_1, id_1, d_c, id_c);
  bool bound_2c = c_base->chk_conn_d_gmi(d_2, id_2, d_c, id_c);

  bool valid_curr = (same_1c and (same_2c or bound_2c)) or 
                    (bound_1c and same_2c);

  std::cout << "e1: " << e1 << " " << d_1 << "c" << id_1
            << " e2: " << e2 << " " << d_2 << "c" << id_2
            << " ec: " << ec << " " << d_c << "c" << id_c
            << " valid_mdl: " << valid_curr << std::endl;

  if(f_calc and f_calc->get_proj() == (int) PROJ_TYPE::EXT_SHELL) {
    ext_shell* e_sh = f_calc->get_e_sh();
    cell_base* s_base = &e_sh->sh_base;
    shell sh_1(3, 1);
    shell sh_2(3, 1);

    if(e_sh->chk_shell(d_1, id_1 - 1))
      sh_1 = e_sh->get_shell(d_1, id_1 - 1);
    if(e_sh->chk_shell(d_2, id_2 - 1))
      sh_2 = e_sh->get_shell(d_2, id_2 - 1);
    std::cout << "e1: " << e1 << " " << sh_1.dim << "sh" << sh_1.id
              << " e2: " << e2 << " " << sh_2.dim << "sh" << sh_2.id
              << " valid_mdl: " << valid_curr << std::endl;

    bool sh_valid = ((sh_1.dim == 3 or sh_2.dim == 3) or sh_1 == sh_2 or s_base->chk_conn_d_gmi(sh_1.dim, sh_1.id, sh_2.dim, sh_2.id));
    valid_curr = valid_curr and sh_valid;
  }
  std::cout << "valid_shell: " << valid_curr << std::endl;

  return valid_curr;
}

bool vd_lens_col::chk_cell() {

  // The right handed pair of edges given axis in the v_i -> v_j direction
  // tet_2v_2e[i][j][0] and tet_2v_2e[i][j][0]
  int tet_2v_2e[4][4][2] = {{{-1, -1}, {3, 4}, {0, 1}, {2, 5}}, 
                            {{1, 2}, {-1, -1}, {4, 5}, {0, 3}}, 
                            {{3, 5}, {0, 2}, {-1, -1}, {1, 4}}, 
                            {{0, 4}, {1, 5}, {2, 3}, {-1, -1}}};
  // The right handed triangle given axis in the v_i -> v_j direction
  // tet_2v_t[i][j]
  int tet_2v_t[4][4] = {{-1, 1, 0, 3}, 
                        {0, -1, 2, 1}, 
                        {3, 0, -1, 2}, 
                        {1, 2, 3, -1}};

  // In a tet, given two vertices, ids of the two tris not adjacent to both 
  // vertices.
  int lookup_surf_n_gg [4][4][2] = {{{-1,-1},{2,3},{1,2},{0,2}},
                                    {{2,3},{-1,-1},{1,3},{0,3}},
                                    {{1,2},{1,3},{-1,-1},{0,1}},
                                    {{0,2},{0,3},{0,1},{-1,-1}}};

  apf::Downward dv;
  apf::Downward de;
  apf::Downward ds;

  std::map<apf::MeshEntity*, bool> e_chk{};

  for(int i = 0; i < elem_col.size(); i ++) {
    m->getDownward(elem_col.at(i), 0, dv);
    m->getDownward(elem_col.at(i), 1, de);
    m->getDownward(elem_col.at(i), 2, ds);
    int v1 = findIn(dv, 4, vert_mer.at(0));
    int v2 = findIn(dv, 4, vert_mer.at(1));

    int e1 = tet_2v_2e[v1][v2][0];
    int e2 = tet_2v_2e[v1][v2][1];
    int e3 = tet_2v_2e[v2][v1][0];
    int e4 = tet_2v_2e[v2][v1][1];

    int s1 = lookup_surf_n_gg[v1][v2][0];
    int s2 = lookup_surf_n_gg[v1][v2][1];

    int sc1 = tet_2v_t[v1][v2];
    int sc2 = tet_2v_t[v2][v1];

    if(!e_chk[de[e1]]) {
      bool valid_curr = chk_ents(ds[sc1], de[e1], de[e2]);
      std::cout << "e1: " << de[e1]
                << " e2: " << de[e2]
                << std::endl;
      valid.first = valid.first and valid_curr;
      e_chk[de[e1]] = true;
      e_chk[de[e2]] = true;
    }
    if(!e_chk[de[e3]]) {
      bool valid_curr = chk_ents(ds[sc2], de[e3], de[e4]);
      std::cout << "e3: " << de[e3]
                << " e4: " << de[e4]
                << std::endl;
      valid.first = valid.first and valid_curr;
      e_chk[de[e3]] = true;
      e_chk[de[e4]] = true;
    }
    if(!e_chk[ds[s1]]) {
      bool valid_curr = chk_ents(elem_col.at(i), ds[s1], ds[s2]);
      std::cout << "s1: " << ds[s1]
                << " s2: " << ds[s2]
                << std::endl;
      valid.first = valid.first and valid_curr;
      e_chk[ds[s1]] = true;
      e_chk[ds[s2]] = true;
    }
  }
  return valid.first;
}

bool vd_lens_col::chk_inv() {
  //std::vector<apf::MeshEntity*> edge_temp(0);
  //std::vector<apf::MeshEntity*> surf_temp(0);
  //std::vector<apf::MeshEntity*> elem_temp(0);

  //vd_set_up(m, &vert_mer, &edge_temp);
  //vd_set_up(m, &edge_temp, &surf_temp);
  //vd_set_up(m, &surf_temp, &elem_temp);

  //vd_remove_set(&elem_temp, &elem_col);

  apf::Vector3 v1(0,0,0);
  apf::Vector3 v2(0,0,0);

  m->getPoint(vert_mer.at(0), 0, v1);
  m->getPoint(vert_mer.at(1), 0, v2);

  m->setPoint(vert_mer.at(0), 0, midpoint);
  m->setPoint(vert_mer.at(1), 0, midpoint);

  vd_save_vtk_ent(m, edge_col, "./output/edge_col_inv");

  bool inv = false;
  //for(int i = 0; i < elem_temp.size(); i++) {
  //  if(!vd_volume_tet_sign(m, elem_temp.at(i) ) ) {
  for(int i = 0; i < elem_sur.size(); i++) {
    if(!vd_volume_tet_sign(m, elem_sur.at(i) ) ) {
      inv = true;
      i = elem_sur.size();
    }
  }
  m->setPoint(vert_mer.at(0), 0, v1);
  m->setPoint(vert_mer.at(1), 0, v2);

  return inv;
}

void vd_lens_col::chk_adj() {
  std::vector<apf::MeshEntity*> edges_1(0);
  std::vector<apf::MeshEntity*> edges_2(0);

  std::vector<apf::MeshEntity*> temp_s(0);
  std::vector<apf::MeshEntity*> temp_e(0);
  apf::Downward dv;

  vd_set_down(m, &elem_col, &temp_s);
  vd_set_down(m, &temp_s, &temp_e);

  vd_set_up(m, vert_mer.at(0), &edges_1);
  vd_set_up(m, vert_mer.at(1), &edges_2);

  vd_remove_set(&edges_1, &temp_e);
  vd_remove_set(&edges_2, &temp_e);

  std::map<apf::MeshEntity*, int> v_count{};

  for(int i = 0; i < edges_1.size(); i++) {
    m->getDownward(edges_1.at(i), 0, dv);
    if(dv[0] == vert_mer.at(0) ) {
      v_count[dv[1]] = v_count[dv[1]] + 1;
    }
    else if(dv[0] == vert_mer.at(1) ) {
      v_count[dv[1]] = v_count[dv[1]] + 1;
    }
    else if(dv[1] == vert_mer.at(0) ) {
      v_count[dv[0]] = v_count[dv[0]] + 1;
    }
    else {
      assert(dv[1] == vert_mer.at(1) );
      v_count[dv[0]] = v_count[dv[0]] + 1;
    }
  }

  for(int i = 0; i < edges_2.size(); i++) {
    m->getDownward(edges_2.at(i), 0, dv);
    if( (dv[0] == vert_mer.at(0) and v_count[dv[1]] == 1) or
        (dv[0] == vert_mer.at(1) and v_count[dv[1]] == 1) or
        (dv[1] == vert_mer.at(0) and v_count[dv[0]] == 1) or
        (dv[1] == vert_mer.at(1) and v_count[dv[0]] == 1) ) {
      std::cout << "invalid adj, e2: " << edges_2.at(i)
                << " dv[0]: " << dv[0] << " " << v_count[dv[0]]
                << " dv[1]: " << dv[1] << " " << v_count[dv[1]]
                << std::endl;
      valid.first = false;
      return;
    }
/*
    if(dv[0] == vert_mer.at(0) ) {
      if(v_count[dv[1]] == 1)
        return false;
    }
    else if(dv[0] == vert_mer.at(1) ) {
      if(v_count[dv[1]] == 1)
        return false;
    }
    else if(dv[1] == vert_mer.at(0) ) {
      if(v_count[dv[0]] == 1)
        return false;
    }
    else {
      assert(dv[1] == vert_mer.at(1) );
      if(v_count[dv[0]] == 1)
        return false;
    }
*/
  }
  valid.first = true;
  return;
}

void vd_lens_col::fix_adj() {
  std::vector<apf::MeshEntity*> edges_1(0);
  std::vector<apf::MeshEntity*> edges_2(0);

  std::vector<apf::MeshEntity*> temp_s(0);
  std::vector<apf::MeshEntity*> temp_e(0);

  vd_set_down(m, &elem_col, &temp_s);
  vd_set_down(m, &temp_s, &temp_e);

  vd_set_up(m, vert_mer.at(0), &edges_1);
  vd_set_up(m, vert_mer.at(1), &edges_2);

  vd_remove_set(&edges_1, &temp_e);
  vd_remove_set(&edges_2, &temp_e);

  std::map<apf::MeshEntity*, int> v_count{};
  apf::Downward dv;
  for(int i = 0; i < edges_1.size(); i++) {
    m->getDownward(edges_1.at(i), 0, dv);
    if(dv[0] == vert_mer.at(0) ) {
      v_count[dv[1]] = v_count[dv[1]] + 1;
    }
    else if(dv[0] == vert_mer.at(1) ) {
      v_count[dv[1]] = v_count[dv[1]] + 1;
    }
    else if(dv[1] == vert_mer.at(0) ) {
      v_count[dv[0]] = v_count[dv[0]] + 1;
    }
    else {
      assert(dv[1] == vert_mer.at(1) );
      v_count[dv[0]] = v_count[dv[0]] + 1;
    }
  }

  bool fixed = false;
  field_calc* f_calc_temp = NULL;
  if(m->findField("velocity_field"))
    f_calc_temp = f_calc;
  for(int i = 0; i < edges_2.size(); i++) {
    m->getDownward(edges_2.at(i), 0, dv);
    if( (dv[0] == vert_mer.at(0) and v_count[dv[1]] == 1) or
        (dv[0] == vert_mer.at(1) and v_count[dv[1]] == 1) or
        (dv[1] == vert_mer.at(0) and v_count[dv[0]] == 1) or
        (dv[1] == vert_mer.at(1) and v_count[dv[0]] == 1) ) {
      std::cout << "invalid adj, e2: " << edges_2.at(i)
                << " dv[0]: " << dv[0] << " " << v_count[dv[1]]
                << " dv[1]: " << dv[1] << " " << v_count[dv[0]]
                << std::endl;
      vd_split_edge(m, c_base, edges_2.at(i), f_calc_temp);
      fixed = true;
    }
/*
    if(dv[0] == vert_mer.at(0) ) {
      if(v_count[dv[1]] == 1) {
        vd_split_edge(m, c_base, edges_2.at(i), f_calc_temp);
        fixed = true;
      }
    }
    else if(dv[0] == vert_mer.at(1) ) {
      if(v_count[dv[1]] == 1) {
        vd_split_edge(m, c_base, edges_2.at(i), f_calc_temp);
        fixed = true;
      }
    }
    else if(dv[1] == vert_mer.at(0) ) {
      if(v_count[dv[0]] == 1) {
        vd_split_edge(m, c_base, edges_2.at(i), f_calc_temp);
        fixed = true;
      }
    }
    else {
      assert(dv[1] == vert_mer.at(1) );
      if(v_count[dv[0]] == 1) {
        vd_split_edge(m, c_base, edges_2.at(i), f_calc_temp);
        fixed = true;
      }
    }
*/
  }
  valid.second = fixed;
}

// Check the stratum and shell memberships of the merging entities.
// Get the stratum membership and position of the final vertex.
// Check and fix unintended merging entities by splitting one of the 
// associated edges.
// Remove the collapsing entities, recreate the surrounding entities 
// adjacent to the vertex to be deleted by replacing the vertex with the 
// remaining vertex. 
void vd_lens_col::col_ent() {
  apf::MeshEntity* ent_n(NULL);

  m->setPoint(vert_ctr, 0, midpoint);

  for(int i = 0; i < elem_col.size(); i++) {
    m->destroy(elem_col.at(i));
  }
  for(int i = 0; i < elem_sur.size(); i++) {
    m->destroy(elem_sur.at(i));
  }

  for(int i = 0; i < surf_col.size(); i++) {
    m->destroy(surf_col.at(i));
  }
  for(int i = 0; i < surf_sur.size(); i++) {
    m->destroy(surf_sur.at(i));
  }

  for(int i = 0; i < edge_sur.size(); i++) {
    m->destroy(edge_sur.at(i));
  }

  for(int i = 0; i < surf_mer.size(); i++) {
    m->destroy(surf_mer.at(i));
  }

  for(int i = 0; i < edge_mer.size(); i++) {
    m->destroy(edge_mer.at(i));
  }

  m->destroy(edge_col);
  m->destroy(vert_other);

  m->acceptChanges();

  for(int i = 0; i < edge_mer_vert.size(); i++) {
    ent_n = build_repl(edge_mer_vert.at(i), edge_mer_mdl.at(i));
    std::cout << "edge_mer " << ent_n << std::endl;
  }
  for(int i = 0; i < surf_mer_vert.size(); i++) {
    ent_n = build_repl(surf_mer_vert.at(i), surf_mer_mdl.at(i));
    std::cout << "surf_mer " << ent_n << std::endl;
  }

  for(int i = 0; i < edge_sur_vert.size(); i++) {
    ent_n = build_repl(edge_sur_vert.at(i), edge_sur_mdl.at(i));
    std::cout << "edge_sur " << ent_n << std::endl;
  }
  for(int i = 0; i < surf_sur_vert.size(); i++) {
    ent_n = build_repl(surf_sur_vert.at(i), surf_sur_mdl.at(i));
    std::cout << "surf_sur " << ent_n << std::endl;
  }
  elem_sur_new.resize(elem_sur_vert.size());
  for(int i = 0; i < elem_sur_vert.size(); i++) {
    elem_sur_new.at(i) = build_repl(elem_sur_vert.at(i), elem_sur_mdl.at(i));
    std::cout << "elem_sur " << elem_sur_new.at(i) << std::endl;
  }

  m->acceptChanges();
  m->verify();

}

std::pair<bool, bool> vd_lens_col::col_lens(apf::MeshEntity* v_cls) {

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

  if(valid.first) {
    chk_adj();
    if(!valid.first) {
      std::vector<apf::MeshEntity*> edge_col_set(1, edge_col);
      // Try to fix merging ents outside lens, otherwise not valid.
      fix_adj();
      if(valid.second)
        valid.first = true;
      vd_set_up(m, &vert_mer, &edge_sur);
      vd_set_up(m, &edge_sur, &surf_sur);
      vd_set_up(m, &surf_sur, &elem_sur);

      std::vector<apf::MeshEntity*> surf_temp(0);
      std::vector<apf::MeshEntity*> edge_temp(0);
      vd_set_down(m, &elem_col, &surf_temp);
      vd_set_down(m, &surf_temp, &edge_temp);

      vd_remove_set(&elem_sur, &elem_col);
      vd_remove_set(&surf_sur, &surf_temp);
      vd_remove_set(&edge_sur, &edge_temp);

      //vd_remove_set(&elem_sur, &elem_col);
      //vd_remove_set(&surf_sur, &surf_col);
      //vd_remove_set(&edge_sur, &edge_col_set);
    }
  }

  if(!valid.first)
    return valid;

  apf::ModelEntity* mdl1 = m->toModel(vert_mer.at(0));
  apf::ModelEntity* mdl2 = m->toModel(vert_mer.at(1));

  int d_1 = m->getModelType(mdl1);
  int id_1 = m->getModelTag(mdl1);
  int d_2 = m->getModelType(mdl2);
  int id_2 = m->getModelTag(mdl2);

  if(v_cls != NULL and !(prs_bound and (d_1 < d_2 or d_2 < d_1)) ) {
    apf::Downward v_d;
    m->getDownward(edge_col, 0, v_d);
    assert(findIn(v_d, 2, v_cls) > -1);

    apf::Vector3 v_cls_pos;
    m->getPoint(v_cls, 0, v_cls_pos);

    v_cls_pos = (v_cls_pos*9 + midpoint)/10;

    if(f_calc and f_calc->get_proj() == (int) PROJ_TYPE::EXT_SHELL) {
      ext_shell* e_sh = f_calc->get_e_sh();
      cell_base* s_base = &e_sh->sh_base;

      shell sh_1;

      if(e_sh->chk_shell(type_vertex, tag_vertex - 1)) {
        sh_1 = e_sh->get_shell(type_vertex, tag_vertex - 1);
        midpoint = e_sh->find_int_pos(sh_1, v_cls_pos);
      }
    }
  }

  if(chk_inv()) {
    valid.first = false;
    return valid;
  }
  collect_merge();
  col_ent();
  //m->verify();

  vd_rem_tag(m);
  vd_tag_mesh(m);

  if(f_calc != NULL and m->findField("velocity_field"))
    f_calc->vd_att_fields(m, vert_ctr);

  if (save_vtk) {
    vd_rem_tag(m);
    vd_tag_mesh(m);

    std::stringstream ss;
    ss << vtk_name << "Edge-" << type_vertex << "c" << tag_vertex;
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    vd_save_vtk_vert(m, &vert_mer, cstr);
  }
  return valid;
}

std::pair<bool, bool> vd_lens_col::load_edge(apf::MeshEntity* edge_in) {

  assert(m->getType(edge_in) == apf::Mesh::EDGE);

  clear();
  edge_col = edge_in;
  vd_save_vtk_ent(m, edge_col, "./output/edge_col");

  vd_set_up(m, edge_col, &surf_col);
  vd_set_up(m, &surf_col, &elem_col);

  vd_set_down(m, edge_col, &vert_mer);

  apf::ModelEntity* mdl1 = m->toModel(vert_mer.at(0));
  apf::ModelEntity* mdl2 = m->toModel(vert_mer.at(1));
  apf::ModelEntity* mdlc = m->toModel(edge_col);

  int d_1 = m->getModelType(mdl1);
  int id_1 = m->getModelTag(mdl1);
  int d_2 = m->getModelType(mdl2);
  int id_2 = m->getModelTag(mdl2);

  int d_c = m->getModelType(mdlc);
  int id_c = m->getModelTag(mdlc);

  bool same_1c = (mdl1 == mdlc);
  bool same_2c = (mdl2 == mdlc);
  bool bound_1c = c_base->chk_conn_d_gmi(d_1, id_1, d_c, id_c);
  bool bound_2c = c_base->chk_conn_d_gmi(d_2, id_2, d_c, id_c);

  valid.first = (same_1c and (same_2c or bound_2c)) or 
                    (bound_1c and same_2c);

  //valid = c_base->chk_conn_d_gmi(d_1, id_1, d_2, id_2);

  if(mdl1 == mdl2 or d_1 < d_2) {
    vert_ctr_em = mdl1;
    vert_ctr = vert_mer.at(0);
    vert_other = vert_mer.at(1);
  }
  else if (d_2 < d_1) {
    vert_ctr_em = mdl2;
    vert_ctr = vert_mer.at(1);
    vert_other = vert_mer.at(0);
  }
  else {
    assert(d_1 == d_2);
    vert_ctr_em = mdl1;
    vert_ctr = vert_mer.at(0);
    vert_other = vert_mer.at(1);
  }
  std::vector<apf::MeshEntity*> edge_col_set(1, edge_col);

  std::vector<apf::MeshEntity*> edge_temp(0);
  std::vector<apf::MeshEntity*> surf_temp(0);
  std::vector<apf::MeshEntity*> elem_temp(0);

  vd_set_up(m, &vert_mer, &edge_sur);
  vd_set_up(m, &edge_sur, &surf_sur);
  vd_set_up(m, &surf_sur, &elem_sur);

  vd_set_down(m, &elem_col, &surf_temp);
  vd_set_down(m, &surf_temp, &edge_temp);

  //vd_remove_set(&surf_sur, &surf_col);
  //vd_remove_set(&edge_sur, &edge_col_set);
  vd_remove_set(&elem_sur, &elem_col);
  vd_remove_set(&surf_sur, &surf_temp);
  vd_remove_set(&edge_sur, &edge_temp);

  type_vertex = m->getModelType(vert_ctr_em);
  tag_vertex = m->getModelTag(vert_ctr_em);

  apf::Vector3 v1(0,0,0);
  apf::Vector3 v2(0,0,0);

  m->getPoint(vert_mer.at(0), 0, v1);
  m->getPoint(vert_mer.at(1), 0, v2);

  midpoint = (v1+v2)/2;
  if(f_calc and f_calc->get_proj() == (int) PROJ_TYPE::EXT_SHELL) {

    ext_shell* e_sh = f_calc->get_e_sh();
    cell_base* s_base = &e_sh->sh_base;

    shell sh_1;

    if(e_sh->chk_shell(type_vertex, tag_vertex - 1)) {
      sh_1 = e_sh->get_shell(type_vertex, tag_vertex - 1);
      midpoint = e_sh->find_int_pos(sh_1, midpoint);
    }
  }
  if(prs_bound and (d_1 < d_2 or d_2 < d_1))
    m->getPoint(vert_ctr, 0, midpoint);

  chk_cell();
  return valid;
}

void vd_lens_col::add_merge(apf::MeshEntity* e1, apf::MeshEntity* e2) {
  apf::Downward dv;

  apf::ModelEntity* mdl1 = m->toModel(e1);
  apf::ModelEntity* mdl2 = m->toModel(e2);
  apf::ModelEntity* mdl_min = NULL;

  int ent_type = m->getType(e1);
  int d = m->typeDimension[ent_type];
  std::vector<apf::MeshEntity*> vert(d + 1);

  int d_1 = m->getModelType(mdl1);
  int id_1 = m->getModelTag(mdl1);
  int d_2 = m->getModelType(mdl2);
  int id_2 = m->getModelTag(mdl2);

  if(mdl1 == mdl2 or d_1 < d_2) {
    mdl_min = mdl1;
  }
  else if (d_2 < d_1) {
    mdl_min = mdl2;
  }
  m->getDownward(e1, 0, dv);
  for(int i = 0; i < d + 1; i++)
    vert.at(i) = dv[i];

  int v1 = findIn(dv, d + 1, vert_other);
  if(v1 != -1) {
    vert.at(v1) = vert_ctr;
  }
  if(d == 1) {
    edge_mer_vert.push_back(vert);
    edge_mer_mdl.push_back(mdl_min);
  }
  else {
    assert(d == 2);
    surf_mer_vert.push_back(vert);
    surf_mer_mdl.push_back(mdl_min);
  }
}

void vd_lens_col::add_surr(apf::MeshEntity* e) {
  apf::Downward dv;

  apf::ModelEntity* mdl = m->toModel(e);

  int ent_type = m->getType(e);
  int d = m->typeDimension[ent_type];
  std::vector<apf::MeshEntity*> vert(d + 1);

  m->getDownward(e, 0, dv);
  for(int i = 0; i < d + 1; i++)
    vert.at(i) = dv[i];

  int v1 = findIn(dv, d + 1, vert_other);
  if(v1 != -1) {
    vert.at(v1) = vert_ctr;
    assert(findIn(dv, d + 1, vert_ctr) == -1);
  }
  if(d == 1) {
    edge_sur_vert.push_back(vert);
    edge_sur_mdl.push_back(mdl);
  }
  else if (d == 2) {
    surf_sur_vert.push_back(vert);
    surf_sur_mdl.push_back(mdl);
  }
  else {
    assert(d == 3);
    elem_sur_vert.push_back(vert);
    elem_sur_mdl.push_back(mdl);
  }
}

void vd_lens_col::clear() {
  for(int i = 0; i < edge_mer_vert.size(); i ++) {
    edge_mer_vert.at(i).clear();
  }
  edge_mer_vert.clear();
  for(int i = 0; i < surf_mer_vert.size(); i ++) {
    surf_mer_vert.at(i).clear();
  }
  surf_mer_vert.clear();

  for(int i = 0; i < edge_sur_vert.size(); i ++) {
    edge_sur_vert.at(i).clear();
  }
  edge_sur_vert.clear();
  for(int i = 0; i < surf_sur_vert.size(); i ++) {
    surf_sur_vert.at(i).clear();
  }
  surf_sur_vert.clear();
  for(int i = 0; i < elem_sur_vert.size(); i ++) {
    elem_sur_vert.at(i).clear();
  }
  elem_sur_vert.clear();
  elem_sur_new.clear();

  elem_sur_mdl.clear();
  surf_sur_mdl.clear();
  edge_sur_mdl.clear();
  edge_mer_mdl.clear();
  surf_mer_mdl.clear();

  surf_sur.clear();
  edge_sur.clear();

  vert_mer.clear();
  elem_col.clear();
  surf_col.clear();
  surf_mer.clear();
  edge_mer.clear();

}

void vd_lens_col::collect_merge() {
  // The right handed pair of edges given axis in the v_i -> v_j direction
  // tet_2v_2e[i][j][0] and tet_2v_2e[i][j][0]
  int tet_2v_2e[4][4][2] = {{{-1, -1}, {3, 4}, {0, 1}, {2, 5}}, 
                            {{1, 2}, {-1, -1}, {4, 5}, {0, 3}}, 
                            {{3, 5}, {0, 2}, {-1, -1}, {1, 4}}, 
                            {{0, 4}, {1, 5}, {2, 3}, {-1, -1}}};

  // In a tet, given two vertices, ids of the two tris not adjacent to both 
  // vertices.
  int lookup_surf_n_gg [4][4][2] = {{{-1,-1},{2,3},{1,2},{0,2}},
                                    {{2,3},{-1,-1},{1,3},{0,3}},
                                    {{1,2},{1,3},{-1,-1},{0,1}},
                                    {{0,2},{0,3},{0,1},{-1,-1}}};
  apf::Downward dv;
  apf::Downward de;
  apf::Downward ds;

  std::map<apf::MeshEntity*, bool> e_chk{};

  edge_mer_vert.reserve(4*elem_col.size());
  surf_mer_vert.reserve(2*elem_col.size());
  edge_mer_mdl.reserve(4*elem_col.size());
  surf_mer_mdl.reserve(2*elem_col.size());

  elem_sur_vert.reserve(elem_sur.size());
  elem_sur_mdl.reserve(elem_sur.size());
  surf_sur_vert.reserve(surf_sur.size());
  surf_sur_mdl.reserve(surf_sur.size());
  edge_sur_vert.reserve(edge_sur.size());
  edge_sur_mdl.reserve(edge_sur.size());

  vd_set_up(m, edge_col, &surf_col);
  vd_set_up(m, &surf_col, &elem_col);

  std::vector<apf::MeshEntity*> edge_col_temp(1, edge_col);
  vd_set_down(m, &surf_col, &edge_mer);
  vd_remove_set(&edge_mer, &edge_col_temp);
  vd_set_down(m, &elem_col, &surf_mer);
  vd_remove_set(&surf_mer, &surf_col);

  vd_set_down(m, edge_col, &vert_mer);

  std::vector<apf::MeshEntity*> edge_temp(0);
  std::vector<apf::MeshEntity*> surf_temp(0);

  vd_set_up(m, &vert_mer, &edge_sur);
  vd_set_up(m, &edge_sur, &surf_sur);
  vd_set_up(m, &surf_sur, &elem_sur);
  vd_set_down(m, &elem_col, &surf_temp);
  vd_set_down(m, &surf_temp, &edge_temp);

  vd_remove_set(&elem_sur, &elem_col);
  vd_remove_set(&surf_sur, &surf_temp);
  vd_remove_set(&edge_sur, &edge_temp);

  for(int i = 0; i < elem_col.size(); i ++) {
    m->getDownward(elem_col.at(i), 0, dv);
    m->getDownward(elem_col.at(i), 1, de);
    m->getDownward(elem_col.at(i), 2, ds);
    int v1 = findIn(dv, 4, vert_mer.at(0));
    int v2 = findIn(dv, 4, vert_mer.at(1));

    int e1 = tet_2v_2e[v1][v2][0];
    int e2 = tet_2v_2e[v1][v2][1];
    int e3 = tet_2v_2e[v2][v1][0];
    int e4 = tet_2v_2e[v2][v1][1];

    int s1 = lookup_surf_n_gg[v1][v2][0];
    int s2 = lookup_surf_n_gg[v1][v2][1];

    if(!e_chk[de[e1]]) {
      add_merge(de[e1], de[e2]);
      e_chk[de[e1]] = true;
      e_chk[de[e2]] = true;
    }
    if(!e_chk[de[e3]]) {
      add_merge(de[e3], de[e4]);
      e_chk[de[e3]] = true;
      e_chk[de[e4]] = true;
    }

    if(!e_chk[ds[s1]]) {
      add_merge(ds[s1], ds[s2]);
      e_chk[ds[s1]] = true;
      e_chk[ds[s2]] = true;
    }
  }
  for(int i = 0; i < elem_sur.size(); i ++) {
    add_surr(elem_sur.at(i));
  }
  for(int i = 0; i < surf_sur.size(); i ++) {
    add_surr(surf_sur.at(i));
  }
  for(int i = 0; i < edge_sur.size(); i ++) {
    add_surr(edge_sur.at(i));
  }

}

void vd_lens_col::assign_surr_map(std::map<apf::MeshEntity*, apf::MeshEntity*>& map_in) {
  for(int i = 0; i < elem_sur.size(); i++) {
    map_in[elem_sur.at(i)] = elem_sur_new.at(i);    
  }
}

void vd_lens_col::assign_col_map(std::map<apf::MeshEntity*, bool>& map_in) {
  for(int i = 0; i < elem_col.size(); i++) {
    map_in[elem_col.at(i)] = true;
  }
}

vd_lens_col::vd_lens_col(apf::Mesh2* msh, struct cell_base* c, field_calc* f_in) :
    m(NULL), mdl(NULL), c_base(NULL), f_calc(NULL),
    save_vtk(false), vtk_name("./output/"), smb_name("./tempmesh/temp"),
    edge_col(NULL), vert_ctr(NULL), vert_ctr_em(NULL), vert_other(NULL),
    type_vertex(0), tag_vertex(0), 
    edge_sur(0), surf_sur(0), elem_sur(0), vert_mer(0), 
    surf_col(0), elem_col(0), elem_sur_new(0),
    surf_mer(0), edge_mer(0),
    edge_sur_vert(0, std::vector<apf::MeshEntity*> (2)),
    surf_sur_vert(0, std::vector<apf::MeshEntity*> (3)),
    elem_sur_vert(0, std::vector<apf::MeshEntity*> (4)),
    elem_sur_mdl(0), surf_sur_mdl(0), edge_sur_mdl(0),

    edge_mer_vert(0, std::vector<apf::MeshEntity*> (2)),
    surf_mer_vert(0, std::vector<apf::MeshEntity*> (3)),
    edge_mer_mdl(0), surf_mer_mdl(0),
    midpoint(0, 0, 0), valid(true, false), prs_bound(false) {

  m = msh;
  mdl = m->getModel();
  c_base = c;

  f_calc = f_in;
}

std::pair<int, int> vd_lens_col::get_cell() {
  return std::make_pair(type_vertex, tag_vertex);
}

void vd_lens_col::set_field_calc(field_calc* f_in) {
  f_calc = f_in;
}

// Set the output file prefixes.
void vd_lens_col::set_files(const char* vtkFile, const char* meshFile) {
  vtk_name = vtkFile;
  smb_name = meshFile;
}

void vd_lens_col::turn_save_vtk(bool onoff) {
  save_vtk = onoff;
}

// Used after collapsing internal structure within the lens to get the new 
// vertex.
apf::MeshEntity* vd_lens_col::get_vert_ctr() {
  return vert_ctr;
}

void vd_lens_col::preserve_bound(bool flag_in) {
  prs_bound = flag_in;
}

vd_lens_col::~vd_lens_col() {
  l_clear();
}

std::pair<bool, bool> vd_col_edge(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* edge, std::map<apf::MeshEntity*, apf::MeshEntity*>& tet_surr_map, std::map<apf::MeshEntity*, bool>& tet_col_map, field_calc* f_calc, bool preserv_bound = false) {
  int ent_type = m->getType(edge);
  int d = m->typeDimension[ent_type];
  assert(d == 1);

  vd_lens_col* lens_col = new vd_lens_col(m, c_base, f_calc);

  std::pair<bool, bool> valid(false, false);
  if(lens_col->load_edge(edge).first) {
    lens_col->preserve_bound(preserv_bound);
    valid = lens_col->col_lens();
    m->acceptChanges();
    if(f_calc != NULL and m->findField("velocity_field"))
      f_calc->vd_att_fields(m, lens_col->get_vert_ctr());

    lens_col->assign_surr_map(tet_surr_map);
    lens_col->assign_col_map(tet_col_map);
  }

  delete lens_col;
  return valid;
}

std::pair<bool, bool> vd_col_edge(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* edge, field_calc* f_calc, bool preserv_bound) {
  int ent_type = m->getType(edge);
  int d = m->typeDimension[ent_type];
  assert(d == 1);

  vd_lens_col* lens_col = new vd_lens_col(m, c_base, f_calc);

  std::pair<bool, bool> valid(false, false);
  if(lens_col->load_edge(edge).first) {
    lens_col->preserve_bound(preserv_bound);
    valid = lens_col->col_lens();
    m->acceptChanges();
    if(f_calc != NULL and m->findField("velocity_field"))
      f_calc->vd_att_fields(m, lens_col->get_vert_ctr());
  }

  delete lens_col;
  return valid;
}

typedef std::pair<double, int> e_len_type;
bool e_len_comp ( const e_len_type& l, const e_len_type& r) { return l.first < r.first; }


// Try to collapse the tet.
// Return true if possible.
// If l_min > sqrt(A_min*2)/10, where l_min is the shortest edge length and A 
// is the smallest triangle area, split the longest edge of the triangle and 
// try to collapse the edge bisecting the triangle. 
// Otherwise, find the shortest collapsable edge of a tet and collapse if 
// possible.
std::pair<bool, bool> vd_col_edge_tet(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* tet, std::map<apf::MeshEntity*, apf::MeshEntity*>& tet_surr_map, std::map<apf::MeshEntity*, bool>& tet_col_map, field_calc* f_calc, bool preserv_bound) {
  int ent_type = m->getType(tet);
  int d = m->typeDimension[ent_type];
  assert(d == 3);

  double vol_tet = vd_meas_ent(m, tet);

  std::vector<e_len_type> t_area(4);
  std::vector<e_len_type> e_len(6);
  apf::Downward dv;
  apf::Downward de;
  apf::Downward ds;
  m->getDownward(tet, 1, de);
  m->getDownward(tet, 2, ds);
  for(int i = 0; i < 6; i++) {
    e_len.at(i).first = vd_meas_ent(m, de[i]);
    e_len.at(i).second = i;
  }
  std::sort(e_len.begin(), e_len.end(), e_len_comp);

  for(int i = 0; i < 4; i++) {
    t_area.at(i).first = vd_meas_ent(m, ds[i]);
    t_area.at(i).second = i;
  }
  std::sort(t_area.begin(), t_area.end(), e_len_comp);

  std::pair<bool, bool> valid(false, false);
  if(e_len.at(0).first > std::sqrt(t_area.at(0).first*2)/10) {
    m->getDownward(ds[t_area.at(0).second], 1, de);
    m->getDownward(ds[t_area.at(0).second], 0, dv);
    e_len.resize(3);
    for(int i = 0; i < 3; i++) {
      e_len.at(i).first = vd_meas_ent(m, de[i]);
      e_len.at(i).second = i;
    }
    std::sort(e_len.begin(), e_len.end(), e_len_comp);

    vd_lens* lens_split = new vd_lens(m, c_base, NULL);

    int i_max = e_len.at(2).second;
    apf::MeshEntity* vert1 = dv[(i_max + 1) % 3];

    lens_split->load_edge(de[i_max]);
    lens_split->split_lens();
    m->acceptChanges();
    apf::MeshEntity* vert2 = lens_split->get_vert_ctr();
    if(f_calc != NULL and m->findField("velocity_field"))
      f_calc->vd_att_fields(m, vert2);

    delete lens_split;

    apf::MeshEntity* edge_new(NULL);
    assert(vd_find_edge(m, vert1, vert2, &edge_new));    

    vd_lens_col* lens_col = new vd_lens_col(m, c_base, f_calc);
    valid = lens_col->load_edge(edge_new);
    lens_col->preserve_bound(preserv_bound);

    if(valid.first) {
      valid = lens_col->col_lens();
      m->acceptChanges();
      if(f_calc != NULL and m->findField("velocity_field"))
        f_calc->vd_att_fields(m, lens_col->get_vert_ctr());
    }
    valid.second = true;
    delete lens_col;
  }
  else if(e_len.at(0).first > std::cbrt(vol_tet*6)/10) {

    vd_bipy* lens_split = new vd_bipy(m, c_base, NULL);

    int i_max = t_area.at(3).second;
    int lookup_tet_x_surf [4] = {3, 2, 0, 1};
    apf::MeshEntity* vert1 = dv[lookup_tet_x_surf [i_max]];

    lens_split->load_tri(ds[i_max]);
    lens_split->split_bipy();

    m->acceptChanges();
    apf::MeshEntity* vert2 = lens_split->get_vert_ctr();
    if(f_calc != NULL and m->findField("velocity_field"))
      f_calc->vd_att_fields(m, vert2);

    delete lens_split;

    apf::MeshEntity* edge_new(NULL);
    assert(vd_find_edge(m, vert1, vert2, &edge_new));    

    vd_lens_col* lens_col = new vd_lens_col(m, c_base, f_calc);
    valid = lens_col->load_edge(edge_new);
    lens_col->preserve_bound(preserv_bound);
    if(valid.first) {
      valid = lens_col->col_lens();
      m->acceptChanges();
      if(f_calc != NULL and m->findField("velocity_field"))
        f_calc->vd_att_fields(m, lens_col->get_vert_ctr());
    }
    valid.second = true;
    delete lens_col;

  }
  else {
    for(int i = 0; i < 6; i++) {
      vd_lens_col* lens_col = new vd_lens_col(m, c_base, f_calc);

      valid = std::pair<bool, bool> (false, false);
      if(lens_col->load_edge(de[e_len.at(i).second]).first) {
        lens_col->preserve_bound(preserv_bound);
        valid = lens_col->col_lens();
        if(valid.first) {
          m->acceptChanges();
          if(f_calc != NULL and m->findField("velocity_field"))
            f_calc->vd_att_fields(m, lens_col->get_vert_ctr());

          lens_col->assign_surr_map(tet_surr_map);
          lens_col->assign_col_map(tet_col_map);
          i = 6;
        }
      }
      delete lens_col;
    }
  }

  return valid;
}
