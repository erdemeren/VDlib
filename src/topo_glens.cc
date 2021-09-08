#include <vector>
#include <algorithm>    
#include <functional>
#include <cstring>

#include <math.h>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>
#include <gmi_null.h>

#include "topo_extinfo.h"
#include "topo_geom.h" // Also contains pi.

#include "topo_topo.h"

#include "topo_ma.h"

/* Topology manipulator:  */
#include "topo_manip.h"
#include "topo_energy.h"

#include "topo_tmesh.h"

#include "topo_field.h"

#include "topo_glens.h"
#include "topo_chull.h"


// Lookups for entity ordering, given the order of the collapsing edge.
// The elements are considered in a ringlike fashion, so only a couple 
// of entities of each kind are processed.
// Depending on the vertex to be removed, the order of the indices might
// change. Assuming the first second one is to be removed:
int lookup_vert_col_g [6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
// Only check the ccw face, with the kept vertex on the top. First 2 ccw, last 
// 2 cw:
//int lookup_edge_mer_g [6][4] = {{2,1,3,4},{0,2,4,5},{3,5,0,1},{0,4,2,5},{1,5,0,3},{2,3,1,4}};
int lookup_edge_mer_g [6][4] = {{2,1,4,3},{0,2,5,4},{3,5,1,0},{0,4,5,2},{1,5,3,0},{2,3,4,1}};

int lookup_edge_remain_g [6] = {5,3,4,1,2,0};
int lookup_surf_mer_g [6][2] = {{3,2},{1,3},{1,2},{0,2},{0,3},{0,1}};
int lookup_surf_col_g [6][2] = {{0,1},{0,2},{3,0},{1,3},{2,1},{3,2}};
int lookup_vert_out_g [6][2] = {{2,3},{0,3},{3,1},{1,2},{2,0},{0,1}};
// The edge id bounded by vertices i and j in a tet.
int lookup_vert_edge_g [4][4] = {{-1,0,2,3},{0,-1,1,4},{2,1,-1,5},{3,4,5,-1}};
// Used in finding the shared triangles of two surrounding tets. Also pick any 
// edge connecting the triangles.
// Triangle opposite vertex on tet.
int lookup_surf_opp_g [4] = {2,3,1,0};

// In a tet, id of tris next to a vertex.
int lookup_surf_n_g [4][3] = {{0,1,3}, {0,1,2}, {0,2,3}, {1,2,3}};
// In a tet, given two vertices, ids of the two tris not adjacent to both 
// vertices.
int lookup_surf_n_gg [4][4][2] = {{{-1,-1},{2,3},{1,2},{0,2}},
                                  {{2,3},{-1,-1},{1,3},{0,3}},
                                  {{1,2},{1,3},{-1,-1},{0,1}},
                                  {{0,2},{0,3},{0,1},{-1,-1}}};
// In a tet, id of vertices adjacent a tri.
int lookup_v_n_s [4][3] = {{0,1,2}, {0,1,3}, {1,2,3}, {0,2,3}};
// In a tet, id of edge across edge.
int lookup_edge_x_e [6] = {5,3,4,1,2,0};

// In a tet, ids of the vertices other than given vertex.
int lookup_v_other [4][3] = {{1,2,3}, {0,2,3}, {0,1,3}, {0,1,2}};

// The edges next to each vertex on the tet.
int lookup_edge_next_g [4][3] = {{0,2,3},{0,1,4},{1,2,5},{3,4,5}};
// The edges next to vertices bounding edge i, on the tet.
int lookup_edge_next_e [6][4] = {{1,2,3,4},{0,2,4,5},{0,1,3,5},{0,2,4,5},
                                 {0,1,3,5},{1,2,3,4}};

// The edges across the vertex on each triangle on the tet.
int lookup_edge_x_g [4][3] = {{1,4,5},{2,3,5},{0,3,4},{0,1,2}};
// The ids of edges in lookup_edge_next_g list. For each edge across the vertex
// on each triangle on the tet. Used to get the index to be used in vector 
// lists.
int lookup_edge_x_bound_id [4][3][2] = {{{0,1},{0,2},{1,2}},
                                        {{0,1},{0,2},{1,2}},
                                        {{0,1},{1,2},{0,2}},
                                        {{0,1},{1,2},{0,2}}};

// On a triangle, given two vertices, gives the id of the two edges not joining
// the vertices.
int lookup_tri_edge_n_g [3][3][2] = {{{-1,-1}, {1,2}, {0,1}},
                                 {{1,2}, {-1,-1}, {0,2}},
                                 {{0,1}, {0,2}, {-1,-1}}};
// On a triangle, given a vertex, gives the id of the two edges next to the  
// vertex.
int lookup_tri_edge_v_g [3][2] = {{0, 2}, {0,1}, {1,2}};
// On a triangle, given a vertex, gives the id of the edge across the vertex.
int lookup_tri_edge_x_g [3] = {1, 2, 0};
// enum for the flag keeping.
//enum {split_6, split_5, split_4, split_3, split_2, split_1, split_0,
//keep_0, keep_1, keep_2, keep_3, keep_4, keep_5, keep_6} merg_flags;
//int merg_index = [6,5,4,3,2,1,0,0,1,2,3,4,5,6];
enum merg_flags_g {split_1, split_0, keep_0, keep_1, keep_2, keep_3, col_fl,
                   END};
const char* merg_flag_str_g[] = {"sp_1", "sp_0", "kp_0", "kp_1", "kp_2", 
                               "kp_3", "col"};
int merg_index_g[] = {1,0,0,1,2,3,0};


vd_relax_cont::vd_relax_cont(apf::Mesh2* m_in, cell_base* cb_in, 
                             field_calc* f_calc_in,
              int dim_in, int id_in,
              std::vector<apf::MeshEntity*> &v_sp_set_in,
              std::vector<std::vector<apf::MeshEntity*> > &ent_in,
              std::vector<apf::MeshEntity*> &et_in, 
              std::map<apf::MeshEntity*, int> &e_map_in,
              bool ext_cell_in, 
              bool upd_shell_in,
              apf::Vector3 midpoint_in, double r_in, double v_th_tol_in) : 
                  // Initialization of mesh related information 
                  mesh(m_in), c_base(cb_in),
                  f_calc(f_calc_in), ext_cell(ext_cell_in),

                  dim_col(dim_in), id_col(id_in),
                  upd_shell(upd_shell_in),
                  midpoint(0,0,0),
                  r_edge_min(r_in), 
                  ent (4, std::vector<apf::MeshEntity*>(0)),
                  // Initialization of entity vectors and maps 
                  ee(0), es(0), et(0),
                  e_map{},
                  tet_surr{}, tet_skip{}, tet_inv{}, tet_i1{}, 
                  v_skip{}, t_nbr{},
                  vol_t(et_in.size()),
                  tet_v (et_in.size(), std::vector<apf::MeshEntity*>(4)),
                  tet_v_id (et_in.size(), std::vector<int>(4)),

                  v_tet (v_sp_set_in.size(), std::vector<int>(0)),
                  v_vert (v_sp_set_in.size(), std::vector<int>(0)),

                  pos_ori{}, v_pos{}, v_disp{}, 

                  A(et_in.size(), 0),
                  B(v_sp_set_in.size(), 
                     std::vector<apf::Vector3> (0, apf::Vector3(0,0,0)) ),
                  C(et_in.size(), 0),
                  exp_pw_tot(v_sp_set_in.size(), 0),

                  tagnumbering(NULL),
                  tetset(NULL),
                  f_field(NULL),
                  // Conjugate gradient initialization
                  x0(v_sp_set_in.size(), apf::Vector3(0,0,0)),
                  x1(v_sp_set_in.size(), apf::Vector3(0,0,0)),
                  g0(v_sp_set_in.size(), apf::Vector3(0,0,0)),
                  g1(v_sp_set_in.size(), apf::Vector3(0,0,0)),
                  dir(v_sp_set_in.size(), apf::Vector3(0,0,0)),
                  ndir(v_sp_set_in.size(), apf::Vector3(0,0,0)),
                  vol_scale(0), vol_max(0), vol_min(-1), 
                  vol_neg(0),
                  epsilon(std::numeric_limits<double>::min()*10e10),
                  v_th_tol(v_th_tol_in), v_th(0.01), w(100),

                  inv_flag(false),
                  iter_limit(400),
                  iter(0), g_norm(0), 
                  xa(0), xb(0), xc(0), xd(0), xe(0), xu(0),
                  fa(0), fb(0), fc(0), fd(0), fe(0), fu(0),
                  p(0), q(0), r(0), s(0), t(0), m(0), tol(0), tol2(0),
                  inv_quad_step(false),
                  Phi0(0), Phi1(0)
{
  midpoint = midpoint_in;
  v_sp_set = v_sp_set_in;
  et = et_in;
  e_map = e_map_in;
  assert(ent_in.size() == 4);
  for(int i = 0; i < 4; i++)
    ent.at(i) = ent_in.at(i);

  f_field = mesh->findField("f_restore");
  /////////////////////////////////////////////
  // Collect adjacency related information:  //
  /////////////////////////////////////////////

  apf::Downward d_v;
  apf::Downward d_e;
  apf::Vector3 temp(0,0,0);
  std::vector<apf::Vector3> pts (4, apf::Vector3(0,0,0));
  std::vector<apf::MeshEntity*>::iterator it;

  for(int i = 0; i < v_sp_set.size(); i++) {
    t_nbr[v_sp_set.at(i)] = 0;
  }

  for(int i = 0; i < et.size(); i++) {
    mesh->getDownward(et.at(i), 0, d_v);
    for(int j = 0; j < 4; j++) {
      tet_v.at(i).at(j) = d_v[j];
      mesh->getPoint(d_v[j], 0, temp);
      v_pos[d_v[j]] = temp;

      it = std::find (v_sp_set.begin(), v_sp_set.end(), tet_v.at(i).at(j));
      if(it != v_sp_set.end()) {
        int i1 = std::distance(v_sp_set.begin(), it);
        tet_v_id.at(i).at(j) = i1;
        t_nbr[v_sp_set.at(i1)] = t_nbr[v_sp_set.at(i1)] + 1;
        v_disp[v_sp_set.at(i1)] = apf::Vector3(0,0,0);
      }
      else
        tet_v_id.at(i).at(j) = -1;
    }
  }

  // Reserve the vertex-tetrahedra lists:
  for(int i = 0; i < v_sp_set.size(); i++) {
    v_tet.at(i).reserve(t_nbr[v_sp_set.at(i)]);
    // In case exterior, there could be one more vertex:
    v_vert.at(i).reserve(t_nbr[v_sp_set.at(i)] + 1);

    mesh->getPoint(v_sp_set.at(i), 0, temp);
    pos_ori[v_sp_set.at(i)] = temp;
    x0.at(i) = temp;
  }

  shell sh_curr;
  v_skip.clear();

  if(upd_shell and f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
    if(ext_cell) {
      for(int i = 0; i < v_sp_set.size(); i++) {
        apf::ModelEntity* mdl = mesh->toModel(v_sp_set.at(i));
        int c_dim = mesh->getModelType(mdl);
        int c_tag = mesh->getModelTag(mdl);
        if(c_dim < 3 and c_base->get_cell_ext_gmi(c_dim, c_tag)) {
          v_skip[v_sp_set.at(i)] = true;
        }
      }
    }
  }

  // Fill the vertex-tetrahedra and vertex-vertex lists for split vertices:
  for(int i = 0; i < et.size(); i++) {
    for(int j = 0; j < 4; j++) {
      if(tet_v_id.at(i).at(j) > -1)
        v_tet.at(tet_v_id.at(i).at(j)).push_back(i);
    }
    for(int j = 0; j < 4; j++) {
      int i1 = tet_v_id.at(i).at(j);
      if(i1 > -1) {
        for(int k = j+1; k < 4; k++) {
          int i2 = tet_v_id.at(i).at(k);
          if(i2 > -1) {
            v_vert.at(i1).push_back(i2);
            v_vert.at(i2).push_back(i1);
          }
        }
      }
    }
  }

  // Remove the repeated entities in the vertex-vertex lists for split vertices:
  for(int i = 0; i < v_sp_set.size(); i++) {
    std::sort(v_vert.at(i).begin(), v_vert.at(i).end());
    std::vector<int>::iterator it;
    it = std::unique (v_vert.at(i).begin(), v_vert.at(i).end());
    v_vert.at(i).resize(std::distance(v_vert.at(i).begin(),it));
  }

  for (int i = 0; i < v_sp_set.size(); i++) {
    B.at(i).resize(v_tet.at(i).size());
  }

  // Tag inverting
  //apf::Numbering* a_f;
  if (mesh->findNumbering("inv_tag")) {
    //printf("grain_memb tag already exists.\n");
    tagnumbering = mesh->findNumbering("inv_tag");
    apf::destroyNumbering(tagnumbering);
    tagnumbering = apf::createNumbering(mesh, "inv_tag", apf::getConstant(3),1);
  }
  else {
    tagnumbering = apf::createNumbering(mesh, "inv_tag", apf::getConstant(3),1);
  }
  //apf::Numbering* a_f;
  if (mesh->findNumbering("tet_tag")) {
    //printf("grain_memb tag already exists.\n");
    tetset = mesh->findNumbering("tet_tag");
    apf::destroyNumbering(tetset);
    tetset = apf::createNumbering(mesh, "tet_tag", apf::getConstant(3),1);
  }
  else {
    tetset = apf::createNumbering(mesh, "tet_tag", apf::getConstant(3),1);
  }
  //if (mesh->findNumbering("a_f")) {
    //printf("grain_memb tag already exists.\n");
  //  a_f = mesh->findNumbering("a_f");
  //  apf::destroyNumbering(a_f);
  //  a_f = apf::createNumbering(mesh, "a_f", apf::getConstant(3),1);
  //}
  //else {
  //  a_f = apf::createNumbering(mesh, "a_f", apf::getConstant(3),1);
  //}

  apf::MeshEntity* e;

  apf::MeshIterator* e_it = mesh->begin(3);
  while ((e = mesh->iterate(e_it))) 
  {
    apf::number(tagnumbering, e, 0, 0, 2);
    apf::number(tetset, e, 0, 0, 0);
    //apf::number(a_f, e, 0, 0, 0);
  }
  mesh->end(e_it);
  
  /////////////////////////////////////////////
  // Collect volume related information:     //
  /////////////////////////////////////////////
  for (int i = 0; i < et.size(); i++) {
    if(e_map[et.at(i)] != 1) {
      mesh->getDownward(et.at(i), 0, d_v);
      for(int j = 0; j < 4; j++) {
        mesh->getPoint(d_v[j], 0, pts.at(j));
      }
      vol_t.at(i) = vd_volume_tet(&pts);

      if(std::fabs(vol_t.at(i)) < vol_min or
         vol_min < - std::numeric_limits<double>::min())
        vol_min = std::fabs(vol_t.at(i));
    }
  }
  // Try total of volumes:
  vol_scale = 0;

  if(dim_col == 3) {
    apf::ModelEntity* mdl = mesh->findModelEntity(dim_col, id_col);
    for (int i = 0; i < et.size(); i++) {
      apf::ModelEntity* mdl_curr = mesh->toModel(et.at(i));
      if(mdl_curr == mdl) {
        assert(e_map[et.at(i)] == 1);
      }
    }
    for (int i = 0; i < ent.at(3).size(); i++) {
      apf::ModelEntity* mdl_curr = mesh->toModel(ent.at(2).at(i));
      if(mdl_curr == mdl) {
        assert(e_map[ent.at(2).at(i)] == 1);
        apf::number(tetset, et.at(i), 0, 0, 1);
      }
    }
  }

  for (int i = 0; i < et.size(); i++) {
    apf::number(tetset, et.at(i), 0, 0, 4);

    if(e_map[et.at(i)] != 1) {
      int ent_type = mesh->getType(et.at(i));
      int d = mesh->typeDimension[ent_type];
      assert(d == 3); 

      mesh->getDownward(et.at(i), 0, d_v);
      mesh->getDownward(et.at(i), 1, d_e);

      std::vector<apf::MeshEntity*>::iterator it;
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;

      int i1 = -1;
      int i2 = -1;
      int i4 = -1;
      //int i3 = -1;
      for (int j = 0; j < 4; j++) {
        it = std::find (ent.at(0).begin(), ent.at(0).end(), d_v[j]);
        //int i_curr = std::distance(ent.at(0).begin(), it);
        if(it != ent.at(0).end()) {
          //std::cout << "d_v[" << j << "] = " << d_v[j] << std::endl;

          if(!found2) {
            if(!found1) {
              found1 = true;
              i1 = j;
              //std::cout << "\tfirst" << std::endl;
            }
            else {
              found2 = true;
              i2 = j;
              //std::cout << "\tsecond" << std::endl;
            }
          }
          else {
            found3 = true;
            //std::cout << "\tthird" << std::endl;
          }
        }
        else {
          i4 = j;
        }
      }

      if(found2) {
        apf::number(tagnumbering, et.at(i), 0, 0, 1);
        apf::number(tetset, et.at(i), 0, 0, 1);
        tet_skip[et.at(i)] = true;
      }
      else if(found1) {
        tet_i1[et.at(i)] = i1;
        // Store the vertices, such that 0th vertex is the merging vertex and
        // the others are on the exterior surface.
        for(int j = 0; j < 4; j++) {
          tet_v.at(i).at(j) = d_v[j];
          mesh->getPoint(d_v[j], 0, pts.at(j));
        }
        //pts.at(0) = midpoint;
        pts.at(i1) = midpoint;
        vol_t.at(i) = vd_volume_tet(&pts);
        if(vol_t.at(i) < epsilon) {
          inv_flag = true;
          tet_inv[et.at(i)] = true;
          //std::cout << et.at(i) << " inverting."
          //          << " vol_t: " << vol_t.at(i)
          //          << std::endl;
          apf::number(tagnumbering, et.at(i), 0, 0, 3);
          if(vol_t.at(i) < vol_neg)
            vol_neg = vol_t.at(i);
        }
        else
          apf::number(tagnumbering, et.at(i), 0, 0, 2);

        if(std::fabs(vol_t.at(i)) > vol_max)
          vol_max = vol_t.at(i);

        apf::number(tetset, et.at(i), 0, 0, 2);
      }
      else {
        tet_surr[et.at(i)] = true;
        // TODO this is currently not used. It should possibly be aiming at 
        // collecting the id of the non-split vertex of the surrounding tet, but
        // is not doing that at the moment. tet_i1 is not being used for 
        // surrounding tets, either...
        tet_i1[et.at(i)] = i4;
        for(int j = 0; j < 4; j++) {
          tet_v.at(i).at(j) = d_v[j];
          mesh->getPoint(d_v[j], 0, pts.at(j));
        }
        vol_t.at(i) = vd_volume_tet(&pts);
        if(std::fabs(vol_t.at(i)) > vol_max)
          vol_max = vol_t.at(i);
        if(vol_t.at(i) < std::numeric_limits<double>::min()) {
          std::cout << "tet_surr " << i << " vol " << vol_t.at(i) << std::endl;
          apf::number(tagnumbering, et.at(i), 0, 0, 3);
        }
        else
          apf::number(tagnumbering, et.at(i), 0, 0, 2);

        apf::number(tetset, et.at(i), 0, 0, 3);
      }
    }
    else {
      vol_t.at(i) = vd_volume_tet(mesh, et.at(i));
      apf::number(tagnumbering, et.at(i), 0, 0, 1);
      apf::number(tetset, et.at(i), 0, 0, 1);
      tet_skip[et.at(i)] = true;
    }

    vol_scale = vol_scale + vol_t.at(i);
  }
  //vol_scale = getMedianEntSize_hist(vol_t);
  //if(vol_scale < std::numeric_limits<double>::min())
  //  vol_scale = vol_max;
  if(vol_scale < std::numeric_limits<double>::min()) {
    std::cout << "Total volume is close to 0." << std::endl;
    vol_scale = 1;
  }

  v_th = std::fabs(v_th_tol*(vol_min+epsilon)/vol_scale);

  for (int i = 0; i < et.size(); i++) {
    if(e_map[et.at(i)] != 1 and !tet_skip[et.at(i)] and !tet_inv[et.at(i)]) {
      if(vol_t.at(i) < v_th*vol_scale) {
        //std::cout << et.at(i) << " the volume " << vol_t.at(i) 
        //          << " is smaller than" << v_th*vol_scale
        //          << std::endl;
        apf::number(tagnumbering, et.at(i), 0, 0, 3);
      }
    }
  }
  for (int i = 0; i < ent.at(3).size(); i++) {
    apf::number(tagnumbering, ent.at(3).at(i), 0, 0, 1);
    apf::number(tetset, ent.at(3).at(i), 0, 0, 1);
  }

  apf::writeVtkFiles("./output/Glens_inv", mesh);

  for (int i = 0; i < ent.at(0).size(); i++) {
    v_pos[ent.at(0).at(i)] = midpoint;
  }

  w = 80*vol_scale/(std::fabs(vol_neg) + epsilon);
  double v_th_i = v_th;
}

// Calculate the terms in the softmax function. B is the 1/6 of the sum of 
// the inward pointing normals of the triangles adjacent to vertex v.
apf::Vector3 vd_relax_cont::calc_B(std::vector<apf::MeshEntity*> &v_set, 
        apf::MeshEntity* v, std::map<apf::MeshEntity*, apf::Vector3> &pos) {

  int v_sense[4] = {1, -1, 1, -1};
  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find (v_set.begin(), v_set.end(), v);
  int i1 = std::distance(v_set.begin(), it);

  apf::Vector3 temp(0,0,0);

  //for(int i = 0; i < 4; i++) {
  //  std::cout << "pos(" << i+1 << ",:) = ["
  //            << pos[v_set.at(i)][0] << ", " 
  //            << pos[v_set.at(i)][1] << ", " 
  //            << pos[v_set.at(i)][2] << "]; " 
  //            << std::endl;
  //}

  int ids[3] = {1, 2, 3};
  apf::Vector3 B_curr(0,0,0);
  for(int i = 0; i < 3; i++) {
    int i2 = (i1+ids[i]) %4;
    int i3 = (i1+ids[(i+1)%3]) %4;
    //temp = vd_cross(pos[v_set.at(i2)] - pos[v_set.at(i1)], 
    //                 pos[v_set.at(i3)] - pos[v_set.at(i1)]);
    //std::cout << "i2 = " << i2+1 << ";";
    //std::cout << "i3 = " << i2+3 << ";";
    //std::cout << "t_n(" << i+1 << ",:) = ["
    //          << temp[0] << ", " 
    //          << temp[1] << ", " 
    //          << temp[2] << "]; " 
    //          << std::endl;
    //B_curr = B_curr + temp;

    //temp = (pos[v_set.at(i2)] + pos[v_set.at(i3)] + pos[v_set.at(i1)])/3;
    //std::cout << "t_pos(" << i+1 << ",:) = ["
    //          << temp[0] << ", " 
    //          << temp[1] << ", " 
    //          << temp[2] << "]; " 
    //          << std::endl;

     B_curr = B_curr + vd_cross(pos[v_set.at(i2)] - pos[v_set.at(i1)], 
                      pos[v_set.at(i3)] - pos[v_set.at(i1)]);

  }
  //std::cout << "B_curr = ["
  //          << B_curr[0] << ", " 
  //          << B_curr[1] << ", " 
  //          << B_curr[2] << "]; " 
  //          << std::endl;
  return B_curr/6*v_sense[i1];
}

// The rate of change of volume wrt. alfa. If v is specified, only consider
// the residual at that vertex. iith vertex is fixed.
double vd_relax_cont::calc_C(std::vector<apf::MeshEntity*> &v, 
                      std::map<apf::MeshEntity*, apf::Vector3> &pos, int ii, 
                                            apf::MeshEntity* v_r) {

  int v_sense[4] = {1, -1, 1, -1};

  std::vector<apf::Vector3> r(4, apf::Vector3(0,0,0));
  double C_curr = 0;

  assert(f_field);
  for(int i = 0; i < 4; i++) {
    if(i != ii) {
      apf::getVector(f_field, v.at(i), 0, r.at(i));
    }
  }

  int ids[3] = {1, 2, 3};
  if(v_r != NULL) {
    std::vector<apf::MeshEntity*>::iterator it;
    it = std::find (v.begin(), v.end(), v_r);
    assert(it != v.end());
    int i1 = std::distance(v.begin(), it);
    // Shift in ids array:
    int id_shift = (i1 + 3 - ii) %4;
    int i2 = (ii+ids[(id_shift+1)%3])%4;
    int i3 = (ii+ids[(id_shift+2)%3])%4;
    C_curr = C_curr + vd_trip_prod(pos[v.at(i3)] - pos[v.at(ii)], 
                     r.at(i1), pos[v.at(i2)] - pos[v.at(ii)]);
  }
  else {
    for(int i = 0; i < 3; i++) {
      int i1 = (ii+ids[(i)%3]) %4;
      int i2 = (ii+ids[(i+1)%3]) %4;
      int i3 = (ii+ids[(i+2)%3]) %4;
      //C_curr = C_curr + vd_trip_prod(pos[v.at(i3)] - pos[v.at(i0)], 
      //                 pos[v.at(i1)] - pos[v.at(i0)], 
      //                 pos[v.at(i2)] - pos[v.at(i0)]);
      C_curr = C_curr + vd_trip_prod(pos[v.at(i3)] - pos[v.at(ii)], 
                       r.at(i1), pos[v.at(i2)] - pos[v.at(ii)]);

    }
  }
  return C_curr/6*v_sense[ii];
}

bool vd_relax_cont::update_vol_inv_grad() {
  double min_vol = v_th*vol_scale;
  std::vector<apf::Vector3> pts(4, apf::Vector3(0,0,0));

  inv_flag = false;
  vol_neg = vol_min;

  for (int i = 0; i < et.size(); i++) {
    if(tet_surr[et.at(i)]) {
      tet_inv[et.at(i)] = false;
      for(int j = 0; j < 4; j++) {
        mesh->getPoint(tet_v.at(i).at(j), 0, pts.at(j));
      }
      vol_t.at(i) = vd_volume_tet(&pts);
      if(vol_t.at(i) < std::numeric_limits<double>::min() and
         vol_t.at(i) < vol_neg)
        vol_neg = vol_t.at(i);

      if(vol_t.at(i) < min_vol) {
        inv_flag = true;
        tet_inv[et.at(i)] = true;
        //std::cout << "tet_surr " << i << " inverting."
        //          << " vol_t: " << vol_t.at(i)
        //          << std::endl;
        apf::number(tagnumbering, et.at(i), 0, 0, 3);
      }
      else
        apf::number(tagnumbering, et.at(i), 0, 0, 2);
      if(std::fabs(vol_t.at(i)) < vol_min or
         vol_min < - std::numeric_limits<double>::min())
        vol_min = std::fabs(vol_t.at(i));
      //if(std::fabs(vol_t.at(i)) > vol_scale)
      //  vol_scale = vol_t.at(i);
    }

    else if(!tet_skip[et.at(i)] and e_map[et.at(i)] != 1) {
      tet_inv[et.at(i)] = false;
      for(int j = 0; j < 4; j++) {
        mesh->getPoint(tet_v.at(i).at(j), 0, pts.at(j));
      }
      pts.at(tet_i1[et.at(i)]) = midpoint;
      vol_t.at(i) = vd_volume_tet(&pts);
      //if(vol_t.at(i) < std::numeric_limits<double>::min()) {
      if(vol_t.at(i) < min_vol) {
        if(vol_t.at(i) < std::numeric_limits<double>::min() and
            vol_t.at(i) < vol_neg)
          vol_neg = vol_t.at(i);

        inv_flag = true;
        tet_inv[et.at(i)] = true;
        //std::cout << "tet " << i << " inverting."
        //          << " vol_t: " << vol_t.at(i)
        //          << std::endl;
        apf::number(tagnumbering, et.at(i), 0, 0, 3);
      }
      else {
        apf::number(tagnumbering, et.at(i), 0, 0, 2);
      }
      //if(std::fabs(vol_t.at(i)) > vol_scale)
      //  vol_scale = vol_t.at(i);
      if(std::fabs(vol_t.at(i)) < vol_min or
         vol_min < - std::numeric_limits<double>::min())
        vol_min = std::fabs(vol_t.at(i));

    }
    else {
      apf::number(tagnumbering, et.at(i), 0, 0, 1);
    }
    assert(!std::isnan(vol_t.at(i)));
  }

  return inv_flag;
}

// Calculate the potential.
double vd_relax_cont::calc_energy() {
  for (int i = 0; i < et.size(); i++) {
    if(e_map[et.at(i)] != 1) {
      A.at(i) = vol_t.at(i)/vol_scale;
      //if(A.at(i) > 10)
      //  A.at(i) = 10; 
      A.at(i) = std::exp(-w*A.at(i));
      //apf::number(a_f, et.at(i), 0, 0, A.at(i));
    }
  }

  for (int i = 0; i < v_sp_set.size(); i++) {
    exp_pw_tot.at(i) = 0;

    bool skipped = true;
    //std::cout << v_sp_set.at(i) << std::endl;
    for (int j = 0; j < v_tet.at(i).size(); j++) {
      int t1 = v_tet.at(i).at(j);
      apf::MeshEntity* tet_curr = et.at(t1);
      if(e_map[tet_curr] != 1 and !tet_skip[tet_curr]) {
        exp_pw_tot.at(i) = exp_pw_tot.at(i) + A.at(t1);
      }
    }
    exp_pw_tot.at(i) = exp_pw_tot.at(i) + 1;
  }

  double Phi = 0;
  for (int i = 0; i < v_sp_set.size(); i++) {
    double Phi_local = std::log(exp_pw_tot.at(i))/w;
    Phi = Phi + Phi_local;
  }
  assert(!std::isnan(Phi));
  return Phi;
}

void vd_relax_cont::calc_Cs() {
  for (int i = 0; i < et.size(); i++) {
    if(e_map[et.at(i)] != 1 and !tet_skip[et.at(i)]) {
      C.at(i) = calc_C(tet_v.at(i), v_pos, tet_i1[et.at(i)])/vol_scale;
    }
    else
      C.at(i) = 0;
  }
}

// TODO this works well enough, but it doesn't consider the cross terms
// grad_i phi_l and only considers grad_i phi_i. See SM of the method paper.
void vd_relax_cont::upd_grad(std::vector<apf::Vector3> &g_curr) {

  for (int i = 0; i < et.size(); i++) {
    if(e_map[et.at(i)] != 1) {
      A.at(i) = vol_t.at(i)/vol_scale;
      //if(A.at(i) > 10)
      //  A.at(i) = 10; 
      A.at(i) = std::exp(-w*A.at(i));

      assert(!std::isinf(A.at(i)));
      assert(!std::isnan(A.at(i)));
      //apf::number(a_f, et.at(i), 0, 0, A.at(i));
    }
  }

  for (int i = 0; i < v_sp_set.size(); i++) {
    g_curr.at(i) = apf::Vector3(0,0,0);
    exp_pw_tot.at(i) = 0;

    //bool skipped = true;
    //std::cout << v_sp_set.at(i) << std::endl;
    for (int j = 0; j < v_tet.at(i).size(); j++) {
      int t1 = v_tet.at(i).at(j);
      apf::MeshEntity* tet_curr = et.at(v_tet.at(i).at(j));
      if(e_map[tet_curr] != 1 and !tet_skip[tet_curr]) {
        //skipped = false;
        B.at(i).at(j) = calc_B(tet_v.at(t1), v_sp_set.at(i), v_pos)/vol_scale;
        //B.at(i).at(j) = norm_0(B.at(i).at(j));

        exp_pw_tot.at(i) = exp_pw_tot.at(i) + A.at(t1);
        //vd_print_vert(mesh, et.at(t1));
      }
      else {
        B.at(i).at(j) = apf::Vector3(0,0,0);
      }
    }
    exp_pw_tot.at(i) = exp_pw_tot.at(i) + 1;

    assert(!std::isinf(exp_pw_tot.at(i)));

    for (int j = 0; j < v_tet.at(i).size(); j++) {
      int t1 = v_tet.at(i).at(j);
      apf::MeshEntity* tet_curr = et.at(t1);
      if(e_map[tet_curr] != 1 and !tet_skip[tet_curr]) {
        g_curr.at(i) = g_curr.at(i) + B.at(i).at(j)*A.at(t1);
      }
    }
    //if(!skipped)
      g_curr.at(i) = g_curr.at(i)/exp_pw_tot.at(i);

    if(v_skip[v_sp_set.at(i)] and f_calc->chk_vert_special(mesh, v_sp_set.at(i))) {
      // TODO Don't use the get_vec_special as it nullifies the motion of shell 2stratum vertices. If possible, find a better solution using get_vec_special.
      apf::ModelEntity* mdl = mesh->toModel(v_sp_set.at(i));
      int em_type = mesh->getModelType(mdl);
      int em_tag = mesh->getModelTag(mdl);

      ext_shell* e_sh;
      e_sh = f_calc->get_e_sh();
      shell sh = e_sh->get_shell(em_type, em_tag-1);
      g_curr.at(i) = e_sh->find_para_dir(em_type, em_tag-1, g_curr.at(i));

      //g_curr.at(i) = f_calc->get_vec_special(mesh, v_sp_set.at(i), g_curr.at(i));
    }
    apf::setVector(f_field, v_sp_set.at(i), 0, g_curr.at(i));
    //std::cout << " f_i: " << g_curr << std::endl;
  }
}

void vd_relax_cont::shift_v_sp_pos(std::vector<apf::Vector3> &x_curr,
                    double x_in) {
  apf::Vector3 pos(0,0,0);
  //apf::Vector3 temp(0,0,0);
  for(int i = 0; i < v_sp_set.size(); i++) {
    //2
    //pos = ndir.at(i)*x_in;
    //assert(!std::isnan(pos.getLength()));
    //pos = x_curr.at(i) + pos;
    //assert(!std::isnan(pos.getLength()));
    //mesh->setPoint(v_sp_set.at(i), 0, pos);

    //3
    if (f_calc->chk_vert_special(mesh, v_sp_set.at(i)))
      pos = f_calc->get_vec_special(mesh, v_sp_set.at(i), ndir.at(i)*x_in);
    else
      pos = ndir.at(i)*x_in;

    assert(!std::isnan(pos.getLength()));
    pos = x_curr.at(i) + pos;
    assert(!std::isnan(pos.getLength()));
    mesh->setPoint(v_sp_set.at(i), 0, pos);

    //1
    //mesh->getPoint(v_sp_set.at(i), 0, temp);
    //temp = f_calc->get_vec_special(mesh, v_sp_set.at(i), x_curr.at(i) + 
    //                                          ndir.at(i)*x_in - temp) + temp;
    //4
    //temp = f_calc->fix_disp_bound(mesh, v_sp_set.at(i), pos, temp);
    //assert(!std::isnan(temp.getLength()));
    //pos = pos + temp;
    //pos = f_calc->fix_vert_bound(mesh, v_sp_set.at(i), pos + temp);
    //assert(!std::isnan(pos.getLength()));
    //mesh->setPoint(v_sp_set.at(i), 0, pos);
    //mesh->setPoint(v_sp_set.at(i), 0, temp);
  }
  update_vol_inv_grad();
}

// Revert the split vertices to their original positions.
void vd_relax_cont::revert_ori() {
  for (int i = 0; i < v_sp_set.size(); i++) {
    mesh->setPoint(v_sp_set.at(i), 0, pos_ori[v_sp_set.at(i)]);
  }
}

bool vd_relax_cont::relax() {
  // Scale it such that most negative volume tetrahedron has A about e^(2)
  // w * v_th = -2*v_th_tol*vol_min/vol_neg
  // w *- vol_t.at(i)/vol_scale = 10/vol_neg * vol_t.at(i)

  double STEP_FIRST = r_edge_min/100;
  //inv_flag = false;

  //inv_flag = update_vol_inv_grad();
  if(!inv_flag) {
    std::cout<< "No inversions remain!" << std::endl;
    return inv_flag;
  }
  //w = 80*vol_scale/vol_neg;
  //if(w*v_th_i > 1)
  //  v_th = std::fabs(v_th_tol*vol_min/vol_scale);
  //else
  //  v_th = v_th_i;

//review the equations of energy and gradient. try manually setting the xb and others to see what makes you jump out of the local minimum, for that's what seems to be happening
  Phi0 = calc_energy();
  upd_grad(g0);

  for (int i = 0; i < v_sp_set.size(); i++) {
    dir.at(i) = g0.at(i)*(-1);
  }
  g_norm = 0;
  for (int i = 0; i < v_sp_set.size(); i++) {
    g_norm = g_norm + dir.at(i) * dir.at(i);
  }
  g_norm = std::sqrt(g_norm);

  for (int i = 0; i < v_sp_set.size(); i++) {
    ndir.at(i) = dir.at(i) / g_norm;
  }

  for (int b = 0; b < MAX_CG_ITER; ++b) {
    ////////////////////////////////////////////////////////////
    // bounds the line search, with ax < bx < cx and fa > fb < fc
    ////////////////////////////////////////////////////////////
    //w = 80*vol_scale/vol_neg;
    //v_th = std::fabs(v_th_tol*vol_min/vol_scale);
    //double v_th_i = v_th;

    // The gradient becomes get reduced so much that the convergence is not 
    // achieved. The weight is rescaled and gradient and energy are recalculated.
    // Below statement is equivalent to abs(vol_neg) < abs(vol_neg_0) and ... but 
    // we are not keeping track of vol_neg.
    // w = 8*vol_scale/(std::fabs(vol_neg) + epsilon);
    // w = 80*vol_scale/(std::fabs(vol_neg) + epsilon);
    if(vol_neg < 2*vol_neg_0 or 
     (vol_neg < -std::numeric_limits<double>::min() and vol_neg > vol_neg_0/10)) {
      w = 80*vol_scale/(std::fabs(vol_neg) + epsilon);
      vol_neg_0 = vol_neg;

      Phi0 = calc_energy();
      upd_grad(g0);

      for (int i = 0; i < v_sp_set.size(); i++) {
        dir.at(i) = g0.at(i)*(-1);
      }
      g_norm = 0;
      for (int i = 0; i < v_sp_set.size(); i++) {
        g_norm = g_norm + dir.at(i) * dir.at(i);
      }
      g_norm = std::sqrt(g_norm);

      if(std::fabs(g_norm) < std::numeric_limits<double>::min())
        g_norm = 1;

      for (int i = 0; i < v_sp_set.size(); i++) {
        ndir.at(i) = dir.at(i) / g_norm;
      }
    }

    xa = 0.;
    fa = Phi0;

    xb = STEP_FIRST;
    shift_v_sp_pos(x0, xb);
    fb = calc_energy();

    while (fb > fa && xb > EPS) {
      // decrease step until energy is decreasing along ndir
      xb /= 10.;
      shift_v_sp_pos(x0, xb);
      fb = calc_energy();
    }

    // try inverse quadratic interpolation
    p = 0;
    for (int i = 0; i < v_sp_set.size(); i++) {
      p = p - g0.at(i)*ndir.at(i);
    }
    p = p * xb;
    q = (fb - fa) + p;
    inv_quad_step = false;
    if (q > EPS) {
      // parabola is concave up, find the minimum
      xc = (p * xb) / (2. * q);
      if (xc > (MAX_BND_STEP + 1.) * xb) {
        // maximum step length
        inv_quad_step = true;
        // MAX_BND_STEP could be too large when xb is not reduced.
        //if (xb/STEP_FIRST < 0.01)
        //  xc = (100 + 1.) * xb;
        //else
        //  xc = (MAX_BND_STEP + 1.) * xb;
        xc = (MAX_BND_STEP + 1.) * xb;

        shift_v_sp_pos(x0, xc);
        fc = calc_energy();
      } else if (xc > (BND_MAG + 1.) * xb) {
        // normal step
        inv_quad_step = true;
        shift_v_sp_pos(x0, xc);
        fc = calc_energy();
        if (fc < fb) {
          // try to step past minimum
          SHIFT2(xa, xb, xc);
          xc = xb + BND_MAG * (xb - xa);
          SHIFT2(fa, fb, fc);
          shift_v_sp_pos(x0, xc);
          fc = calc_energy();
        }
      } else if (xc > xa + SQRT_EPS && xc < xb - SQRT_EPS) {
        // minimum falls in (ax, bx)
        shift_v_sp_pos(x0, xc);
        fc = calc_energy();
        if (fc < fb) {
          // found bracket, all done
          inv_quad_step = true; 
          std::swap(xb, xc);
          std::swap(fb, fc);
        }
      }
    }
    if (!inv_quad_step) {
      // quadratic interpolation failed, conservative step
      xc = (BND_MAG + 1.) * xb;
      shift_v_sp_pos(x0, xc);
      fc = calc_energy();
    }
    while (fc < fb) {
      // try inverse quadratic interpolation
      p = xc - xb;
      q = xa - xb;
      r = (fa - fb) * p;
      s = (fc - fb) * q;
      t = r - s;
      inv_quad_step = false;
      if (t > EPS) {
        // parabola is concave up, find the minimum
        xd = xb + (r * p - s * q) / (2. * t);
        if (xd > xc + MAX_BND_STEP * p) {
          // maximum step length
          inv_quad_step = true;
          xd = xc + MAX_BND_STEP * p;

          shift_v_sp_pos(x0, xd);
          fd = calc_energy();
        } else if (xd > xc + BND_MAG * p) {
          // normal step
          inv_quad_step = true;
          shift_v_sp_pos(x0, xd);
          fd = calc_energy();
          if (fd < fc) {
            // try to step past minimum
            SHIFT3(xa, xb, xc, xd);
            xd = xc + BND_MAG * (xc - xb);
            SHIFT3(fa, fb, fc, fd);
            shift_v_sp_pos(x0, xd);
            fd = calc_energy();
          }
        } else if (xd > xb + SQRT_EPS && xd < xc - SQRT_EPS) {
          // minimum falls in (ax, bx)
          shift_v_sp_pos(x0, xd);
          fd = calc_energy();
          if (fd < fc) {
            // found bracket, all done
            inv_quad_step = true;
            SHIFT2(xa, xb, xd);
            SHIFT2(fa, fb, fd);
            break;
          }
        }
      }
      if (!inv_quad_step) {
        // quadratic interpolation failed, conservative step
        xd = xc + BND_MAG * p;
        shift_v_sp_pos(x0, xd);
        fd = calc_energy();
      }
      // bookkeeping for next iteration
      SHIFT3(xa, xb, xc, xd);
      SHIFT3(fa, fb, fc, fd);
    }

    ////////////////////////////////////////////////////////////
    // Brent's method to find minimum along search direction, a translation
    // of the ALGOL 60 algorithm on page 79 of R. P. Brent, Algorithms for
    // Minimization Without Derivatives, 1973 with minor modifications. The
    // author gave permission to use this algorithm in private communication.
    ////////////////////////////////////////////////////////////
    // use values from bounding the line search
    if (fc < fa) {
      xd = xc;
      xe = xa;
      fd = fc;
      fe = fa;
    } else {
      xd = xa;
      xe = xc;
      fd = fa;
      fe = fc;
    }
    t = (xb < 0.5 * (xa + xc) ? xc : xa) - xb;
    s = PHI_SQ_INV * t;
    int c = 0;
    for (c = 0; c < MAX_LS_ITER; ++c) {
      m = 0.5 * (xa + xc);
      tol  = SQRT_EPS * (fabs(xb) + 1.);
      tol2 = 2. * tol;
      // check stopping criterion
      if (fabs(xb - m) > tol2 - 0.5 * (xc - xa)) {
        inv_quad_step = false;
        if (fabs(t) > tol) {
          // inverse quadratic interpolation
          p = (xb - xd) * (fb - fe);
          q = (xb - xe) * (fb - fd);
          r = (xb - xe) * q - (xb - xd) * p;
          q = 2. * (q - p);
          if (q > 0.) { r = -r; }
          q = fabs(q);
          SHIFT2(p, t, s);
          // mistake in ALGOL 60 routine, second condition is inverted
          if (fabs(r) < fabs(0.5 * q * p) && r > q * (xa - xb) && r < q * (xc - xb)) {
            // take inverse quadratic interpolation step
            inv_quad_step = true;
            s = r / q;
            xu = xb + s;
            // f should not be evaluated too close to xa or xc
            if (xu - xa < tol2 || xc - xu < tol2) {
              s = (xb < m ? tol : -tol);
            }
          }
        }
        if (!inv_quad_step) {
          // interpolation failed, take golden section step
          t = (xb < m ? xc : xa) - xb;
          s = PHI_SQ_INV * t;
        }

        // f should not be evaluated too close to xb
        xu = xb + (fabs(s) >= tol ? s : (s > 0. ? tol : -tol));

        shift_v_sp_pos(x0, xu);
        fu = calc_energy();
        // bookkeeping for next iteration
        if (fu <= fb) {
          if (xu < xb) { xc = xb; } else { xa = xb; }
          SHIFT3(xe, xd, xb, xu);
          SHIFT3(fe, fd, fb, fu);
        } else {
          if (xu < xb) { xa = xu; } else { xc = xu; }
          if (fu <= fd || xd == xb) {
            SHIFT2(xe, xd, xu);
            SHIFT2(fe, fd, fu);
          } else if (fu <= fe || xe == xb || xe == xd) {
            xe = xu;
            fe = fu;
          }
        }
      }
      else {
        for (int i = 0; i < v_sp_set.size(); i++) {
          x1.at(i) = x0.at(i) + ndir.at(i) * xb;
        }

        // found minimum, apply change and update energy
        Phi1 = fb;
        break;
      }
    }
    if (c == MAX_LS_ITER)
      std::cerr << "RunSimulation: max LS iteration reached" << std::endl;

    ////////////////////////////////////////////////////////////
    // Conjugate gradient
    ////////////////////////////////////////////////////////////
    // check energy convergence
    //if (fabs(Phi1 - Phi0) < E_TOL) {
    //  x0 = x1;
    //}

    // check gradient convergence
    shift_v_sp_pos(x1, 0);
    if(!(w < -8*vol_scale/vol_neg or w > -160*vol_scale/vol_neg)) {
      upd_grad(g1);

      // Direction update given by Y.H Dai, C.X. Kou, SIAM Journal of
      // Optimization, v 23, p 296-320, 2013
      for (int i = 0; i < v_sp_set.size(); i++) {
        g0.at(i) = g0.at(i) - g1.at(i);
      }
      for (int i = 0; i < v_sp_set.size(); i++) {
        x0.at(i) = x0.at(i) - x1.at(i);
      }

      double B0 = 0;
      double B1 = 0;
      double B0_a = 0;
      double B0_b = 0;
      double B0_c = 0;

      double B1_a = 0;
      double B1_b = 0;
      for (int i = 0; i < v_sp_set.size(); i++) {
        B0_a = B0_a + g0.at(i)*g0.at(i);
        B0_b = B0_a + g1.at(i)*x0.at(i);
        B0_c = B0_a + g0.at(i)*x0.at(i);
      }
      if(std::fabs(B0_c) < std::numeric_limits<double>::min())
        B0_c = 1;
      B0 = B0_a*B0_b/B0_c;

      for (int i = 0; i < v_sp_set.size(); i++) {
        B0_a = B0_a + g0.at(i)*g1.at(i);
        B0_b = B0_a + dir.at(i)*g0.at(i);

        B1_a = B0_a + g1.at(i)*dir.at(i);
        B1_b = B0_a + dir.at(i)*dir.at(i);
      }
      if(std::fabs(B0_b) < std::numeric_limits<double>::min())
        B0_b = 1;

      B0 = (B0_a - B0)/B0_b;

      if(std::fabs(B1_b) < std::numeric_limits<double>::min())
        B1_b = 1;

      B1 = B1_a / (2. * B1_b);

      double mult = fmax(B0, B1);
      for (int i = 0; i < v_sp_set.size(); i++) {
        dir.at(i) = dir.at(i) * mult - g1.at(i);
      }

      g_norm = 0;
      for (int i = 0; i < v_sp_set.size(); i++) {
        g_norm = g_norm + dir.at(i) * dir.at(i);
      }
      g_norm = std::sqrt(g_norm);

      if(std::fabs(g_norm) < std::numeric_limits<double>::min())
        g_norm = 1;

      for (int i = 0; i < v_sp_set.size(); i++) {
        ndir.at(i) = dir.at(i) / g_norm;
      }

      // prepare for next iteration
      for (int i = 0; i < v_sp_set.size(); i++) {
        x0.at(i) = x1.at(i);
      }
      shift_v_sp_pos(x0, 0);

      Phi0 = Phi1;
      for (int i = 0; i < v_sp_set.size(); i++) {
        g0.at(i) = g1.at(i);
      }

    }
    inv_flag = false;
    //v_th = v_th_tol*vol_min/vol_scale;
    //vol_scale = 0;
    inv_flag = update_vol_inv_grad();

    apf::writeVtkFiles("./output/Glens_step", mesh);

    if(!inv_flag) {
      std::cout<< "No inversions remain!" << std::endl;
      b = MAX_CG_ITER;
    }
    else
      std::cout<< "Threshold volume:" << v_th*vol_scale << std::endl;

    if(std::fabs(vol_neg) > 4*vol_scale) {
      std::cout<< "Negative volume elements are larger than the starting volume, relaxation using the v_th is likely unfeasible. Stopping relaxation. "
               << " v_th: " << v_th
               << " vol_neg: " << vol_neg
               << " vol_scale: " << vol_scale
               << std::endl;
      break;
    }

  }

  if(inv_flag) {
    std::cout << dim_col << "c" << id_col <<  ": Relaxation was not successful!" << std::endl;
    revert_ori();
  }
  else
    std::cout << dim_col << "c" << id_col <<  ": Relaxation was successful!" << std::endl;

  // End conjugate gradient method
  return inv_flag;
}

vd_relax_cont::~vd_relax_cont() {
}
////////////////////////////////////////
// vd_glens
////////////////////////////////////////
/*
void vd_glens::get_v_sp_pos(std::map<MeshEntity*, apf::MeshEntity*> &e2e_map_in, 
                  std::map<MeshEntity*, apf::Vector3> &e_sp_pos_trial_in) {
  e2e_map = e2e_map_in;
  e_sp_pos_trial = e_sp_pos_trial_in;
  calc_sp_pos(false);
}

bool vd_glens::set_v_sp_pos(std::vector<apf::MeshEntity*> &ee) {
  for(int i = 0; i < ee.size; i++) {
    m->setPoint(sp_e2v_map[ee.at(i)], 0, e_sp_pos[ee.at(i)]);
  }
}
*/
// Update the merging entities.
void vd_glens::set_merg_map(apf::MeshEntity* e1, apf::MeshEntity* e2) {

  apf::ModelEntity* mdl_1 = m->toModel(e1);
  apf::ModelEntity* mdl_2 = m->toModel(e2);
/*
  std::cout << m->typeDimension[m->getType(e1)] << "-ents "
            << e1 << " emap[" << e_map[e1] << "] ["
            << m->getModelType(mdl_1) << "c" << m->getModelTag(mdl_1) << "] "
            << e2 << " emap[" << e_map[e2] << "] ["
            << m->getModelType(mdl_2) << "c" << m->getModelTag(mdl_2)
            << "] are merging." << std::endl;
*/
  if(e_map[e1] == 2) {
    mdl_1 = m_map[e1];
  }
  else
    e_map[e1] = 2;
  if(e_map[e2] == 2) {
    mdl_2 = m_map[e2];
  }
  else
    e_map[e2] = 2;
/*
  std::cout << "Updated "
            << m->getModelType(mdl_1) << "c" << m->getModelTag(mdl_1) << " "
            << m->getModelType(mdl_2) << "c" << m->getModelTag(mdl_2)
            << " are merging." << std::endl;
*/
  apf::ModelEntity* mdl_curr;
  if(m->getModelType(mdl_1) < m->getModelType(mdl_2))
    mdl_curr = mdl_1;
  else if(m->getModelType(mdl_2) < m->getModelType(mdl_1))
    mdl_curr = mdl_2;
  else {
    assert(mdl_1 == mdl_2);
    mdl_curr = mdl_1;
  }
/*
  std::cout << "New model is "
            << m->getModelType(mdl_curr) << "c" << m->getModelTag(mdl_curr)
            << std::endl;
*/
  m_map[e1] = mdl_curr;
  m_map[e2] = mdl_curr;
}

// A top down approach to collecting merging entities
// Collect the vertices across the merging edges. Bundle the merging edges by
// joint vertices across. Find the model for the created entities by finding
// the lowest dimensional cell in a bundle. Assert there is a single one.
void vd_glens::set_merg() {
  std::vector<apf::MeshEntity*> elem(0);
  std::vector<apf::MeshEntity*> surf(0);
  apf::Downward d_vert;
  apf::Downward d_edge;
  // The entity across the merging vertices in a given entity 
  // (edge to vertex (tri to edge) across)
  std::map<apf::MeshEntity*, apf::MeshEntity*> e_x_map{};
  // The count of d+1 dimensional merging entities bounded by the d 
  // dimensional entity across the merging vertices.
  std::map<apf::MeshEntity*, int> e_merg_ct{};
  // The lowest dimensional model
  std::map<apf::MeshEntity*, apf::ModelEntity*> m_x_map{};
  apf::ModelEntity* mdl;
  int c_type = -1;
  int c_tag = -1;

  // Index of the vertex across. Map this edge to the vertex. All edges
  // mapped to the same vertex will be merging to a single edge.
  for(int i = 0; i < ent.at(1).size(); i++) {
    m->getDownward(ent.at(1).at(i), 0, d_vert);

    int i1 = findIn(&ent.at(0), ent.at(0).size(), d_vert[0]);
    if(i1 > -1) {
      i1 = findIn(&ent.at(0), ent.at(0).size(), d_vert[1]);
      if(i1 > -1) {
        e_map[ent.at(1).at(i)] = 1;

        vd_set_up(m, ent.at(1).at(i), &surf);
        vd_set_up(m, &surf, &elem);

        for(int j = 0; j < elem.size(); j++) {
          //std::cout << "Tet " << elem.at(j) << std::endl;
          e_map[elem.at(j)] = 1;
        }
        for(int j = 0; j < surf.size(); j++) {
          e_map[surf.at(j)] = 1;
        }
      }
      else {
        e_x_map[ent.at(1).at(i)] = d_vert[1];
        e_merg_ct[d_vert[1]] = e_merg_ct[d_vert[1]] + 1;

        if(m_x_map[d_vert[1]] == NULL)
          m_x_map[d_vert[1]] = m->toModel(ent.at(1).at(i));
        else {
          c_type = m->getModelType(m_x_map[d_vert[1]]);
          mdl = m->toModel(ent.at(1).at(i));
          if(c_type > m->getModelType(mdl))
            m_x_map[d_vert[1]] = mdl;
        }
      }
    }
    else {
      e_x_map[ent.at(1).at(i)] = d_vert[0];
      e_merg_ct[d_vert[0]] = e_merg_ct[d_vert[0]] + 1;

      i1 = findIn(&ent.at(0), ent.at(0).size(), d_vert[1]);
      assert(i1 > -1);

      if(m_x_map[d_vert[0]] == NULL)
        m_x_map[d_vert[0]] = m->toModel(ent.at(1).at(i));
      else {
        c_type = m->getModelType(m_x_map[d_vert[0]]);
        mdl = m->toModel(ent.at(1).at(i));
        if(c_type > m->getModelType(mdl))
          m_x_map[d_vert[0]] = mdl;
      }
    }
  }

  for(int i = 0; i < ent.at(1).size(); i++) {
    m->getDownward(ent.at(1).at(i), 0, d_vert);

    //std::cout << ent.at(1).at(i) << " "
    //          << vd_print_pos(m, ent.at(1).at(i)) 
    //          << " " << d_vert[0] << " " << d_vert[1]
    //          << std::endl;

    if(e_map[ent.at(1).at(i)] != 1) {
      mdl = m->toModel(ent.at(1).at(i));
      apf::MeshEntity* v_x = e_x_map[ent.at(1).at(i)];

      c_type = m->getModelType(m_x_map[v_x]);
      int type_curr = m->getModelType(mdl);
      assert(type_curr >= c_type);
      if(c_type == type_curr)
        assert(mdl == m_x_map[v_x]);
      m_map[ent.at(1).at(i)] = m_x_map[v_x];
      //std::cout << "m_old " << m->getModelType(mdl) << "c" 
      //                      << m->getModelTag(mdl)
      //          << " m_new " << m->getModelType(m_x_map[v_x]) << "c" 
      //                      << m->getModelTag(m_x_map[v_x])
      //          << std::endl;
      if(e_merg_ct[v_x] > 1) {
        e_map[ent.at(1).at(i)] = 2;
        //std::cout << "Vert " << v_x << std::endl;
      }
    }
  }

  for(int i = 0; i < ent.at(2).size(); i++) {
    if(e_map[ent.at(2).at(i)] != 1) {
      m->getDownward(ent.at(2).at(i), 0, d_vert);
      bool found1 = false;
      int i0 = -1;
      for(int j = 0; j < 3; j++) {
        int i1 = findIn(&ent.at(0), ent.at(0).size(), d_vert[j]);
        if(i1 > -1) {
          if(found1) {
            e_map[ent.at(2).at(i)] = 1;
            j = 3;
          }
          else {
            i0 = j;
            found1 = true;
          }
        }
      }
      if(e_map[ent.at(2).at(i)] != 1 and found1) {
        m->getDownward(ent.at(2).at(i), 1, d_edge);
        // Index of the edge across. Map this triangle to the edge. All triangles
        // mapped to the same edge will be merging to a single triangle.
        int e1 = (i0+1) % 3;

        e_x_map[ent.at(2).at(i)] = d_edge[e1];
        e_merg_ct[d_edge[e1]] = e_merg_ct[d_edge[e1]] + 1;

        if(m_x_map[d_edge[e1]] == NULL)
          m_x_map[d_edge[e1]] = m->toModel(ent.at(2).at(i));
        else {
          c_type = m->getModelType(m_x_map[d_edge[e1]]);
          mdl = m->toModel(ent.at(2).at(i));
          if(c_type > m->getModelType(mdl))
            m_x_map[d_edge[e1]] = mdl;
        }
      }
    }
  }
  for(int i = 0; i < ent.at(2).size(); i++) {
    m->getDownward(ent.at(2).at(i), 1, d_edge);

    //std::cout << ent.at(2).at(i) << " "
    //          << vd_print_pos(m, ent.at(2).at(i)) 
    //          << " " << d_edge[0] 
    //          << " " << d_edge[1] 
    //          << " " << d_edge[2] 
    //          << std::endl;
    if(e_map[ent.at(2).at(i)] != 1) {
      mdl = m->toModel(ent.at(2).at(i));

      apf::MeshEntity* e_x = e_x_map[ent.at(2).at(i)];

      c_type = m->getModelType(m_x_map[e_x]);
      int type_curr = m->getModelType(mdl);
      assert(type_curr >= c_type);
      if(c_type == type_curr)
        assert(mdl == m_x_map[e_x]);
      m_map[ent.at(2).at(i)] = m_x_map[e_x];
      //std::cout << "m_old " << m->getModelType(mdl) << "c" 
      //                      << m->getModelTag(mdl)
      //          << " m_new " << m->getModelType(m_x_map[e_x]) << "c" 
      //                      << m->getModelTag(m_x_map[e_x])
      //          << std::endl;

      if(e_merg_ct[e_x] > 1) {
        //std::cout << "Edge " << e_x << " "
        //          << ent.at(2).at(i) << " "
        //          << vd_print_pos(m, ent.at(2).at(i)) << std::endl;
        e_map[ent.at(2).at(i)] = 2;
      }
    }
  }

}

// Go over the collapsing edges. Get the entities around them. For merging 
// entities if one is collapsing, the other should as well.
// If they are not collapsing, assign them merging and find the lowest cell dim
// cell one belongs to. If they are same dim, the cell should be the same.
void vd_glens::set_merg2() {
  for(int i = 0; i < ent.at(1).size(); i++) {
    apf::Downward d_vert;
    m->getDownward(ent.at(1).at(i), 0, d_vert);
    int i1 = findIn(&ent.at(0), ent.at(0).size(), d_vert[0]);
    if(i1 > -1) {
      i1 = findIn(&ent.at(0), ent.at(0).size(), d_vert[1]);
      if(i1 > -1)
        e_map[ent.at(1).at(i)] = 1;
    }
  }
  for(int i = 0; i < ent.at(2).size(); i++) {
    apf::Downward d_vert;
    m->getDownward(ent.at(2).at(i), 0, d_vert);
    bool found1 = false;
    for(int j = 0; j < 3; j++) {
      int i1 = findIn(&ent.at(0), ent.at(0).size(), d_vert[j]);
      if(i1 > -1) {
        if(found1) {
          e_map[ent.at(2).at(i)] = 1;
          j = 3;
        }
        else {
          found1 = true;
        }
      }
    }
  }

  for(int i = 0; i < ent.at(1).size(); i++) {
    if(e_map[ent.at(1).at(i)] == 1) {
      std::vector<apf::MeshEntity* > surf(0);
      std::vector<apf::MeshEntity* > elem(0);

      vd_set_up(m, ent.at(1).at(i), &surf);
      vd_set_up(m, &surf, &elem);

      apf::Downward d_edge;
      apf::Downward d_surf;

      for(int j = 0; j < elem.size(); j++) {
        //std::cout << "Tet " << elem.at(j) << std::endl;
        e_map[elem.at(j)] = 1;

        m->getDownward(elem.at(j), 1, d_edge);
        m->getDownward(elem.at(j), 2, d_surf);

        int i1 = findIn(d_edge, 6, ent.at(1).at(i));
        // e1-e2 and e3-e4 merging edge couples:
        int e1 = lookup_edge_mer_g[i1][0];
        int e2 = lookup_edge_mer_g[i1][1];
        int e3 = lookup_edge_mer_g[i1][2];
        int e4 = lookup_edge_mer_g[i1][3];

        // s1-s2 merging triangle couple and s3-s4 are collapsing:
        int s1 = lookup_surf_mer_g[i1][0];
        int s2 = lookup_surf_mer_g[i1][1];
        int s3 = lookup_surf_col_g[i1][0];
        int s4 = lookup_surf_col_g[i1][1];

        if(e_map[d_edge[e1]] == 1 or e_map[d_edge[e2]] == 1) {
          //assert(e_map[d_edge[e2]] == 1);
          e_map[d_edge[e1]] = 1;
          e_map[d_edge[e2]] = 1;
        }
        else {
          set_merg_map(d_edge[e1], d_edge[e2]);
        }

        if(e_map[d_edge[e3]] == 1 or e_map[d_edge[e4]] == 1) {
          e_map[d_edge[e3]] = 1;
          e_map[d_edge[e4]] = 1;
          //assert(e_map[d_edge[e4]] == 1);
        }
        else {
          set_merg_map(d_edge[e3], d_edge[e4]);
        }

        if(e_map[d_surf[s1]] == 1 or e_map[d_surf[s2]] == 1) {
          e_map[d_surf[s1]] = 1;
          e_map[d_surf[s2]] = 1;
          //assert(e_map[d_surf[s2]] == 1);
        }
        else {
          set_merg_map(d_surf[s1], d_surf[s2]);
        }

        e_map[d_surf[s3]] = 1;
        e_map[d_surf[s4]] = 1;
      }
    }
  }
}

// Load the entities surrounding the collapsing edge.
bool vd_glens::load_cell() {
  if (not c_flag) {
    return false;
  }

  ent.resize(4);
  v_list.resize(3);
  m_list.resize(3);

  //vd_find_ent_geom(m, &ent.at(cell_col_dim), cell_col_id, cell_col_dim);
  ent.at(cell_col_dim) = e_list->e.at(cell_col_dim).at(cell_col_id-1).at(cell_col_dim);

  for(int i = 0; i < ent.at(cell_col_dim).size(); i++) {
    e_map[ent.at(cell_col_dim).at(i)] = 1;
  }
  for(int dim = cell_col_dim; dim > 1; dim--) {
    vd_set_down(m, &ent.at(dim), &ent.at(dim-1));
    for(int i = 0; i < ent.at(dim-1).size(); i++) {
      e_map[ent.at(dim-1).at(i)] = 1;
    }
  }

  vd_set_down(m, &ent.at(1), &ent.at(0));
  ext_c0_vec.reserve(ent.at(0).size());

  apf::Vector3 zero(0,0,0);

  if(tag_0c == -1) {
    for(int i = 0; i < ent.at(0).size(); i++) {
      e_map[ent.at(0).at(i)] = 2;
      v_map[ent.at(0).at(i)] = true;
      apf::ModelEntity* mdl_curr = m->toModel(ent.at(0).at(i));
      if(m->getModelType(mdl_curr) == 0) {
        c0_flag = true;
        vert_ctr_em = mdl_curr;
        vert_ctr = ent.at(0).at(i);
        int c_id = m->getModelTag(mdl_curr);
        if(ext_cell and c_base->get_cell_ext_gmi(0, c_id)) {
          ext_c0_vec.push_back(zero);
          m->getPoint(ent.at(0).at(i), 0, ext_c0_vec.back());
        }
      }
    }
    if(c0_flag) {
      v_map[vert_ctr] = false;
      k_map[vert_ctr] = true;
    }
  }
  else {
    bool found = false;
    for(int i = 0; i < ent.at(0).size(); i++) {
      e_map[ent.at(0).at(i)] = 2;
      v_map[ent.at(0).at(i)] = true;
      apf::ModelEntity* mdl_curr = m->toModel(ent.at(0).at(i));
      if(m->getModelType(mdl_curr) == 0) { 
        if(m->getModelTag(mdl_curr) == tag_0c) {
          c0_flag = true;
          assert(!found);
          found = true;
          vert_ctr_em = mdl_curr;
          vert_ctr = ent.at(0).at(i);
        }
        int c_id = m->getModelTag(mdl_curr);
        if(ext_cell and c_base->get_cell_ext_gmi(0, c_id)) {
          ext_c0_vec.push_back(zero);
          m->getPoint(ent.at(0).at(i), 0, ext_c0_vec.back());
        }
      }
    }
    assert(c0_flag);
    v_map[vert_ctr] = false;
    k_map[vert_ctr] = true;
  }

  if (save_vtk) {
    vd_rem_tag(m);
    vd_tag_mesh(m);

    std::stringstream ss;
    ss << vtk_name << "Glens_"<< cell_col_dim << "c" << cell_col_id;
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    //apf::writeVtkFiles(cstr, m);
    vd_save_vtk_vert(m, &ent.at(0), cstr);
  }

  //TODO split here. Go up once, if it's not collapsing, it is to be split.
  // Going up in the following step will update the entities. 
  createfields();
  if(!corr_lens_new()) {
    destroyfields();
    return false;
  }
  destroyfields();

  e_map.clear();
  m_map.clear();

  for(int i = 0; i < ent.at(0).size(); i++) {
    e_map[ent.at(0).at(i)] = 2;
  }

  for(int dim = 0; dim < 3; dim++) {
    vd_set_up(m, &ent.at(dim), &ent.at(dim+1));
  }
  set_merg();

  // Label the surrounding entities.
  for(int dim = 1; dim < 4; dim++) {
    for(int i = 0; i < ent.at(dim).size(); i++) {
      if(e_map[ent.at(dim).at(i)] == 0)
        e_map[ent.at(dim).at(i)] = 3;
    }
  }

  if(inv_flag)
    assert(!chk_inversion());

  apf::ModelEntity* mdl_vert;
  if(!c0_flag) {
    int c_nbr;
    std::map<apf::ModelEntity*, bool> mdl_map{};
    for(int dim = 1; dim < 4; dim++) {
      c_nbr = 0;
      for(int i = 0; i < ent.at(dim).size(); i++) {
        apf::ModelEntity* mdl_curr = m->toModel(ent.at(dim).at(i));
        std::cout << dim << "-ent " << ent.at(dim).at(i) << " "
                  << m->getModelType(mdl_curr) << "c"
                  << m->getModelTag(mdl_curr)
                  << std::endl;
        if(e_map[ent.at(dim).at(i)] > 1) {
          if(m->getModelType(mdl_curr) == dim and !mdl_map[mdl_curr]) {
            c_nbr++;
            mdl_vert = mdl_curr;
            mdl_map[mdl_curr] = true;
          }
        }
      }
      assert(c_nbr < 2);
      if(c_nbr == 1) {
        dim = 4;
      }
    }
    assert(c_nbr > 0);
  }

  if(!c0_flag) {
    vert_ctr_em = mdl_vert;
    vert_ctr = m->createVert(vert_ctr_em);
    //vert_ctr = f_calc->create_v(m, vert_ctr_em);

    std::cout << "New vert " << vert_ctr
              << " " << m->getModelType(vert_ctr_em)
              << "c" << m->getModelTag(vert_ctr_em)
              << std::endl;

  }

  //m->setPoint(vert_ctr, 0, midpoint);

  type_vertex = m->getModelType(vert_ctr_em);
  tag_vertex = m->getModelTag(vert_ctr_em);

  if(!c0_flag)
    std::cout << "New center vertex " << vert_ctr;
  else
    std::cout << "Center vertex " << vert_ctr << " emap[" << e_map[vert_ctr]
              << "] " << std::endl;
//There is no such vertex!

  std::cout << type_vertex << "c" << tag_vertex
            << std::endl;

  return true;
}
/*
apf::MeshEntity* vd_glens::vd_rebuildElement(apf::MeshEntity* original, 
                                                apf::BuildCallback* cb) {
  int ent_type = m->getType(original);
  int d = m->typeDimension[ent_type];

  apf::Downward vert;
  m->getDownward(original, 0, vert);

  bool repl = false;
  bool has_0c = false;
  for(int i = 0; i < d+1; i++) {
    if(v_map[vert[i]]) {
      assert(!repl);
      vert[i] = vert_ctr;
      repl = true;
    }
    else if(vert[i] == vert_ctr) {
      has_0c = true;
    }
  }

  apf::MeshEntity* ent_n;

  //std::cout << d << "-ent " << original << " ";
  //std::cout << m->getModelType(m->toModel(original)) 
  //            << "c" << m->getModelTag(m->toModel(original));
  // If merging entity
  if(e_map[original] == 2) {
    if(has_0c)
      assert(m_map[original] == m->toModel(original));
    ent_n = buildElement(m, m_map[original], ent_type, vert);
    if(ent_n == original)
      k_map[ent_n] = true;
    //std::cout << " replaced by " << m->getModelType(m_map[original]) 
    //          << "c" << m->getModelTag(m_map[original]);
  }

  // If surrounding entity
  else {
    assert(e_map[original] == 3);

    ent_n = buildElement(m, m->toModel(original), ent_type, vert);
    if(ent_n == original)
      k_map[ent_n] = true;
  }
  //std::cout << " " << d << "-ent " << ent_n << std::endl;

  return ent_n;
}
*/
void vd_glens::get_vert() {
  apf::Downward down;

  for(int dim = 1; dim < 4; dim++) {
    int e_nbr = 0;
    for(int i = 0; i < ent.at(dim).size(); i++) {
      if(e_map[ent.at(dim).at(i)] > 1) {
        e_nbr++;
      }
    }
    m_list.at(dim-1).reserve(e_nbr);
    v_list.at(dim-1).reserve(e_nbr);

    std::vector<apf::MeshEntity*> emp_list(0);
    for(int i = 0; i < ent.at(dim).size(); i++) {
      if(e_map[ent.at(dim).at(i)] > 1) {
        int e_sz = m->getDownward(ent.at(dim).at(i), 0, down);

        bool repl = false;
        bool has_0c = false;
        for(int j = 0; j < e_sz; j++) {
          // These entities should be either merging or surrounding.
          if(v_map[down[j]]) {
            assert(!repl);
            down[j] = vert_ctr;
            repl = true;
          }
          else if(down[j] == vert_ctr) {
            has_0c = true;
          }
        }
        //Recreating merging entities shouldn't cause any problem, but check.
        //Something fishy is going on, as the newly generated surrounding entities are not using vert_ctr!
        if(e_map[ent.at(dim).at(i)] == 2) {
          m_list.at(dim-1).push_back(m_map[ent.at(dim).at(i)]);
        }
        else
          m_list.at(dim-1).push_back(m->toModel(ent.at(dim).at(i)));

        v_list.at(dim-1).push_back(emp_list);
        v_list.at(dim-1).back().resize(e_sz);
        for(int j = 0; j < dim+1; j++) {
          v_list.at(dim-1).back().at(j) = down[j];
        }
      }
    }
  }
}

void vd_glens::recreate_ent() {
  get_vert();
  destroy_ent();

  m->setPoint(vert_ctr, 0, midpoint);

  apf::Downward down;
  for(int i = 0; i < m_list.at(0).size(); i++) {
    //std::cout << m->getModelType(m_list.at(0).at(i)) <<"c" 
    //          << m->getModelTag(m_list.at(0).at(i)) << std::endl << "\t";
    for(int j = 0; j < v_list.at(0).at(i).size(); j++) {
      down[j] = v_list.at(0).at(i).at(j);
      //std::cout << down[j] << " ";
    }
    //std::cout << std::endl; 
    apf::MeshEntity* ent_n = buildElement(m, m_list.at(0).at(i),
                        apf::Mesh::EDGE, down);
    //vd_print_down(m, ent_n);
  }
  for(int i = 0; i < m_list.at(1).size(); i++) {
    //std::cout << dim << "-ent " << ent.at(dim).at(i);
    //std::cout << m->getModelType(m_list.at(1).at(i)) <<"c" 
    //          << m->getModelTag(m_list.at(1).at(i)) << std::endl << "\t";
    for(int j = 0; j < v_list.at(1).at(i).size(); j++) {
      down[j] = v_list.at(1).at(i).at(j);
      //std::cout << down[j] << " ";
    }
    //std::cout << std::endl; 
    apf::MeshEntity* ent_n = buildElement(m, m_list.at(1).at(i),
                        apf::Mesh::TRIANGLE, down);
    //vd_print_down(m, ent_n);
  }
  for(int i = 0; i < m_list.at(2).size(); i++) {
    //std::cout << dim << "-ent " << ent.at(dim).at(i);
    //std::cout << m->getModelType(m_list.at(2).at(i)) <<"c" 
    //          << m->getModelTag(m_list.at(2).at(i)) << std::endl << "\t";
    for(int j = 0; j < v_list.at(2).at(i).size(); j++) {
      down[j] = v_list.at(2).at(i).at(j);
      //std::cout << down[j] << " ";
    }
    //std::cout << std::endl; 
    apf::Downward dv;
    apf::MeshEntity* v_temp;
    apf::MeshEntity* ent_n = buildElement(m, m_list.at(2).at(i),
                        apf::Mesh::TET, down);
    //vd_print_down(m, ent_n);

    // Assuming precondition is not used only in the case of an unsuccessful 
    // insertion. A small neighborhood may lead to very small volume elements.
    if(precond_flag) {
      double meas = vd_volume_tet(m, down);
      //apf::MeshElement* ee = createMeshElement(m, ent_n);
      //double meas = measure(ee);
      //destroyMeshElement(ee);

      assert(meas > std::numeric_limits<double>::min());
    }
/*
    if(meas < std::numeric_limits<double>::min()) {
    //if(vd_volume_tet(m, down) < 0) {
      m->getDownward(ent_n, 0, dv);
      m->destroy(ent_n);
      v_temp = dv[1];
      dv[1] = dv[2];
      dv[2] = v_temp;

      ent_n = buildElement(m, m_list.at(2).at(i), apf::Mesh::TET, dv);

      ee = createMeshElement(m, ent_n);
      double meas = measure(ee);
      destroyMeshElement(ee);
      assert(meas > 0);

      //assert(vd_volume_tet(m, down) > 0);
      //vd_print_vert(m, e_void);
    }
*/
  }
}

void vd_glens::destroy_ent() {
  for(int dim = 3; dim > 0; dim--) {
    for(int i = 0; i < ent.at(dim).size(); i++) {
      //std::cout << dim << "-ent " << ent.at(dim).at(i) << std::endl;
      m->destroy(ent.at(dim).at(i));
    }
  }
  for(int i = 0; i < ent.at(0).size(); i++) {
    //std::cout << 0 << "-ent " << ent.at(0).at(i) << std::endl;
    if(v_map[ent.at(0).at(i)])
      m->destroy(ent.at(0).at(i));
    else {
      //std::cout << "\tKept" << std::endl;
    }
  }
  m->acceptChanges();
}

// TODO Move the merging vertices to the center. Mainly used in debugging.
void vd_glens::set_merg_ctr() {
  for (int i = 0; i < ent.at(0).size(); i++) {
    m->setPoint(ent.at(0).at(i), 0, midpoint);
  }
}

// TODO Move the merging vertices to the center. Mainly used in debugging.
void vd_glens::get_merg_ori() {
  apf::Vector3 temp(0,0,0);
  for (int i = 0; i < ent.at(0).size(); i++) {
    m->getPoint(ent.at(0).at(i), 0, temp);
    pos_merg[ent.at(0).at(i)] = temp;
  }
}

// TODO Move the merging vertices to the center. Mainly used in debugging.
void vd_glens::set_merg_ori() {
  apf::Vector3 temp(0,0,0);
  for (int i = 0; i < ent.at(0).size(); i++) {
    temp = pos_merg[ent.at(0).at(i)];
    m->setPoint(ent.at(0).at(i), 0, temp);
  }
}

// Check if the surrounding elements become inverted after the collapse.
// Print the inverting elements. Return true if inverted.
bool vd_glens::chk_inversion() {

  bool inv = false;
/*
  std::cout << "Inversion before: " << std::endl;
  for (int i = 0; i < ent.at(3).size(); i++) {
    if(e_map[ent.at(3).at(i)] == 3) {
      double meas_ee = vd_volume_tet(m, ent.at(3).at(i));
      if(meas_ee < std::numeric_limits<double>::min()) {

        int type1 = m->getModelType(m->toModel(ent.at(3).at(i)));
        int tag1 = m->getModelTag(m->toModel(ent.at(3).at(i)));
        printf("Element %p %dcell%d inverts, vol: %e.\n", 
                                                (void*)ent.at(3).at(i), 
                                                type1,tag1,meas_ee);
        std::cout << getLinearCentroid(m, ent.at(3).at(i)) << std::endl;
      }
        inv = true;

    }
  }
*/
  std::vector<apf::Vector3> pos(0, apf::Vector3(0,0,0));
  pos.resize(ent.at(0).size());

  for (int i = 0; i < ent.at(0).size(); i++) {
    m->getPoint(ent.at(0).at(i), 0, pos.at(i));
    m->setPoint(ent.at(0).at(i), 0, midpoint);
  }

  std::cout << "Inversion after: " << std::endl;

  inv = false;
  for (int i = 0; i < ent.at(3).size(); i++) {
    if(e_map[ent.at(3).at(i)] == 3) {
      //apf::MeshElement* ee = createMeshElement(m, ent.at(3).at(i));
      //double meas_ee = measure(ee);
      double meas_ee = vd_volume_tet(m, ent.at(3).at(i));
      if(meas_ee < std::numeric_limits<double>::min()) {

        int type1 = m->getModelType(m->toModel(ent.at(3).at(i)));
        int tag1 = m->getModelTag(m->toModel(ent.at(3).at(i)));
        printf("Element %p %dcell%d inverts, vol: %e.\n", 
                                                (void*)ent.at(3).at(i), 
                                                type1,tag1,meas_ee);
        std::cout << getLinearCentroid(m, ent.at(3).at(i)) << std::endl;
        inv = true;
      }
    }
  }

  // Revert the positions.
  for (int i = 0; i < ent.at(0).size(); i++) {
    m->setPoint(ent.at(0).at(i), 0, pos.at(i));
  }
  return inv;
}

bool vd_glens::edge_comp(apf::MeshEntity* e1, apf::MeshEntity* e2, apf::Mesh2* m) {
  apf::ModelEntity* mdl = m->toModel(e1);
  int c1_type = m->getModelType(mdl);

  mdl = m->toModel(e2);

  // Highest cell membership first.
  if(c1_type != m->getModelType(mdl))
    return c1_type > m->getModelType(mdl);

  apf::MeshElement* ee = createMeshElement(m, e1);
  double len1 = measure(ee);
  destroyMeshElement(ee);

  ee = createMeshElement(m, e2);
  double len2 = measure(ee);
  destroyMeshElement(ee);

  // Widest element first.
  return len1 > len2;
}

void vd_glens::sort_edge(std::vector<apf::MeshEntity*>* edge_list) { 
  std::sort(edge_list->begin(), edge_list->end(), 
      std::bind(&vd_glens::edge_comp, this, 
      std::placeholders::_1, std::placeholders::_2, m) );
}

void vd_glens::extract_conn_map(cell_base* conn_map, 
                         std::vector<apf::MeshEntity*> ee,
                         std::vector<apf::MeshEntity*> es_in,
                         std::vector<apf::MeshEntity*> et) { 
  conn_map->clear();
  // Cell ids of all entities. 
  apf::Downward d_v;
  apf::Downward d_e;
  apf::MeshEntity* e_curr;

  std::map<apf::MeshEntity*, int> ent_id_map{};
  std::vector<std::map<int, apf::MeshEntity*> > id_ent_map(3, 
                                    std::map<int, apf::MeshEntity*> {});
  // Going over the edges, go over the adjacent non-collapsing triangles, 
  // and check if the exterior angle is smaller than pi. 
  // If so, split the edge across. Continue until no exterior angle smaller 
  // than pi exists for such triangles.
  // Going over the collapsing tets, check pairs of edges across 
  // (2 pairs intersecting the sphere)
  // The first one to be split among these should be one of the closest pair.
  // If this still has an exterior angle is smaller than pi, split the truss
  // afterwards and project the vertex onto the sphere
  // (the edge upwards of the split vertices of the two edges).
  conn_map->add_free(0, ee.size());
  conn_map->add_free(1, es_in.size());
  conn_map->add_free(2, et.size());

  // Loop over non-collapsing edges
  // ent_id_map, id_ent_map, int_pt, vert_pt
  for (int i = 0; i < ee.size(); i++) {
    int id = conn_map->use_free_gmi(0);
    assert(id > 0);
    ent_id_map[ee.at(i)] = id;
    id_ent_map.at(0)[id] = ee.at(i);
  }
  // Loop over triangles not belonging to the collapsing cell.
  for (int i = 0; i < es_in.size(); i++) {
      //std::cout << "e2_sp " << ee.at(i) << std::endl; 
    if(e_map[es_in.at(i)] != 1) {
      int id = conn_map->use_free_gmi(1);
      assert(id > 0);
      ent_id_map[es_in.at(i)] = id;
      id_ent_map.at(1)[id] = es_in.at(i);

      m->getDownward(es_in.at(i), 0, d_v);
      m->getDownward(es_in.at(i), 1, d_e);

      bool found1 = false;
      bool found2 = false;

      int i1 = -1;
      int i2 = -1;
      for(int j = 0; j < 3; j++) {

        int i_curr = findIn(&ent.at(0), ent.at(0).size(), d_v[j]);
        if(i_curr > -1) {
          assert(!found2);
          if(found1) {
            found2 = true;
            i2 = j;
            assert(vd_find_edge(m, d_v[i1], d_v[i2], &e_curr));
          }
          else {
            found1 = true;
            i1 = j;
          }
        }
      }

      assert(found1);

      int e1 = -1;
      int e2 = -1;
      int e3 = -1;

      if(found2) {
        e1 = findIn(d_e, 3, e_curr);
        e2 = (e1 + 1) % 3;
        e3 = (e2 + 1) % 3;

        assert((lookup_tri_edge_n_g[i1][i2][0] == e2 and
                lookup_tri_edge_n_g[i1][i2][1] == e3) or
               (lookup_tri_edge_n_g[i1][i2][0] == e3 and
                lookup_tri_edge_n_g[i1][i2][1] == e2));
      }
      else {
        bool found_e1 = false;
        bool found_e2 = false;
        for(int j = 0; j < 3; j++) {

          int e_curr = findIn(&ee, ee.size(), d_e[j]);
          if(e_curr > -1) {
            if(found_e1) {
              e3 = j;
              found_e2 = true;
              j = 3;
            }
            else {
              e2 = j;
              found_e1 = true;
            }
          }
        }
        assert((lookup_tri_edge_v_g[i1][0] == e2 and
                lookup_tri_edge_v_g[i1][1] == e3) or
               (lookup_tri_edge_v_g[i1][0] == e3 and
                lookup_tri_edge_v_g[i1][1] == e2));
      }

      int id_1 = ent_id_map[d_e[e2]];
      int id_2 = ent_id_map[d_e[e3]];
      assert(id_1 > 0 and id_2 > 0);

      ent_conn e_con(2);
      e_con.add_ent(id_1);
      e_con.add_ent(id_2);
      conn_map->set_conn_gmi(1, id, &e_con);
      std::cout << id_1 << " " << id_2 << "; ";
    }
  }
  std::cout << "];" << std::endl;
  // Loop over tets not belonging to the collapsing cell.
  apf::Vector3 pos1(0,0,0);
  apf::Vector3 pos2(0,0,0);
  apf::Vector3 l_dir(0,0,0);
  apf::Vector3 l_pos(0,0,0);

  for (int i = 0; i < et.size(); i++) {

    if(e_map[et.at(i)] != 1) {
      int id = conn_map->use_free_gmi(2);
      assert(id > 0);
      ent_id_map[et.at(i)] = id;
      id_ent_map.at(2)[id] = et.at(i);

      apf::Downward d_v;
      apf::Downward d_t;

      m->getDownward(et.at(i), 0, d_v);
      m->getDownward(et.at(i), 2, d_t);

      bool found1 = false;
      bool found2 = false;
      bool found3 = false;

      int i1 = -1;
      int i2 = -1;
      int i3 = -1;
      int i_not = -1;
      for(int j = 0; j < 4; j++) {

        int i_curr = findIn(&ent.at(0), ent.at(0).size(), d_v[j]);
        if(i_curr > -1) {
          assert(!found3);
          if(found2) {
            found3 = true;
            i3 = j;
          }
          else {
            if(found1) {
              found2 = true;
              i2 = j;
            }
            else {
              found1 = true;
              i1 = j;
            }
          }
        }
        else {
          i_not = j;
        }
      }
      ent_conn e_con(4);
      // Three vertices found: Add all tris other than the one bounded by
      // the vertices.
      if(found3) {
        for(int j = 0; j < 3; j++) {
          int id_curr = ent_id_map[d_t[lookup_surf_n_g[i_not][j]]];
          e_con.add_ent(id_curr);
        }
      }
      // Two vertices found: Add all tris.
      else if(found2) {
        for(int j = 0; j < 4; j++) {
          int id_curr = ent_id_map[d_t[j]];
          e_con.add_ent(id_curr);
        }
      }
      // One vertex found: Add all tris other than the one across the vertex.
      else {
        for(int j = 0; j < 3; j++) {
          int id_curr = ent_id_map[d_t[lookup_surf_n_g[i1][j]]];
          e_con.add_ent(id_curr);
        }
      }
      //vd_print_vert(m, et.at(i));
      conn_map->set_conn_gmi(2, id, &e_con);
    }
  }
}

// Given a triangle and a list of intersection points lying on the triangle,
// refine the triangle and the subsequent triangles recursively.
void vd_glens::refine_tri(apf::MeshEntity* tri, std::vector<apf::Vector3> pts) {
  if(pts.size() == 0)
    return;
  std::vector<std::vector<apf::Vector3> > pts_div(3, 
                          std::vector<apf::Vector3> (0, apf::Vector3(0,0,0)));
  apf::Downward d_v;
  apf::Downward d_v2;
  m->getDownward(tri, 0, d_v);

  // Take the points closest to the incenter as the intersection point.
  // Divide the others into the resulting triangles.
  apf::Vector3 i_ctr(0,0,0);

  for(int i = 0; i < pts.size(); i++) {
    i_ctr = i_ctr + pts.at(i);
  }
  i_ctr = i_ctr/pts.size();
  //i_ctr = inctr_ent(m, tri);

  if(pts.size() > 1) {
    pts_div.at(0).reserve(pts.size()-1);
    pts_div.at(1).reserve(pts.size()-1);
    pts_div.at(2).reserve(pts.size()-1);

    double dist_min = -1;
    double dist_curr = -1;
    int i_div = -1;
    for(int i = 0; i < pts.size(); i++) {
      dist_curr = (pts.at(i) - i_ctr).getLength();
      if(dist_min < 0 or dist_curr < dist_min) {
        i_div = i;
      }
    }
    i_ctr = pts.at(i_div);
    pts.erase(pts.begin()+i_div);
  }
  else
    i_ctr = pts.at(0);

  bipy_split->load_tri(tri);
  bipy_split->split_bipy();
  apf::MeshEntity* v_sp = bipy_split->get_vert_ctr();

  //m->setPoint(v_sp, 0, proj_pt/2+l/2);
  m->setPoint(v_sp, 0, i_ctr);
  std::cout << "v_new(1,:) = ["
            << i_ctr[0] << "," 
            << i_ctr[1] << "," 
            << i_ctr[2] << "];" << std::endl;
  assert(vd_chk_neg_vert(m, v_sp) == 0);

  f_calc->vd_att_fields(m, v_sp);

  if(pts.size() > 1) {
    std::vector<apf::MeshEntity*> ee(0);
    std::vector<apf::MeshEntity*> es(0);

    vd_set_up(m, v_sp, &ee);
    vd_set_up(m, &ee, &es);

    for(int i = es.size()-1; i > -1; i--) {
      bool found1 = false;
      bool found2 = false;
      m->getDownward(es.at(i), 0, d_v2);
      for(int j = 0; j < 3; j++) {
        if(findIn(d_v, 3, d_v2[j]) > -1) {
          if(found1) {
            found2 = true;
            j = 3;
          }
          else {
            found1 = true;
          }
        }
      }
      if(!found2) {
        es.erase(es.begin()+i);
      }
    }

    assert(es.size() == 3);
    for(int i = 0; i < pts.size(); i++) {
      bool found = false;
      for(int j = 0; j < es.size(); j++) {
        m->getDownward(es.at(j), 0, d_v2);
        if(pt_int_tri(m, es.at(j), pts.at(i), 0.05)) {
          pts_div.at(j).push_back(pts.at(j));
          found = true;
          j = es.size();
        }
      }
      if(!found) {
        std::cout << "Pts too close, skipping." << std::endl;
      }
    }
    
    for(int j = 0; j < es.size(); j++) {
      if(pts_div.at(j).size() > 0)
        refine_tri(es.at(j), pts_div.at(j));
    }
  }
  // Assign points
}

void vd_glens::get_int_pos2(std::vector<apf::MeshEntity*> ee, 
                           std::map<apf::MeshEntity*, apf::Vector3> &mer_pt,
                           std::map<apf::MeshEntity*, apf::Vector3> &int_pt,
                           std::map<apf::MeshEntity*, apf::Vector3> &vert_pt,
                           double r_cvx) {
  // Loop over non-collapsing edges
  // ent_id_map, id_ent_map, int_pt, vert_pt
  for (int i = 0; i < ee.size(); i++) {
    //std::cout << "e2_sp " << ee.at(i) << std::endl; 
    apf::Downward down;

    apf::MeshEntity* v_in;
    apf::MeshEntity* v_out;

    apf::Vector3 l(0,0,0);
    apf::Vector3 o(0,0,0);
    apf::Vector3 oc(0,0,0);

    m->getDownward(ee.at(i), 0, down);

    std::vector<apf::MeshEntity*>::iterator it;
    it = std::find (ent.at(0).begin(), ent.at(0).end(), down[0]);
    int i1 = std::distance(ent.at(0).begin(), it);

    if(it != ent.at(0).end()) {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);

      //assert(it == ent.at(0).end());
      if(it != ent.at(0).end()) {
        std::cout << "Edge shared by two merging vertices" << std::endl;
      }
      v_in = ent.at(0).at(i1);
      v_out = down[1];
    }
    else {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);
      i1 = std::distance(ent.at(0).begin(), it);

      assert(it != ent.at(0).end());
      v_in = ent.at(0).at(i1);
      v_out = down[0];
    }
    m->getPoint(v_in, 0, o);
    //o = midpoint;
    m->getPoint(v_out, 0, l);

    mer_pt[ee.at(i)] = o;
    vert_pt[ee.at(i)] = l;

    //std::cout << "c_pos " << c_pos << " "
    //          << " v_in " << o << " "
    //          << " v_out " << l
    //          << std::endl;
    double split_len = (l - o).getLength()/2;
    if(r_cvx < - std::numeric_limits<double>::min() or split_len < r_cvx) {
      int_pt[ee.at(i)] = o/2+l/2;
    }
    else {
      l = l - o;
      l = norm_0(l);
      int_pt[ee.at(i)] = o+l*r_cvx;
    }

    std::cout << "vert(" << i+1 << ",:) = [" << vert_pt[ee.at(i)][0] << ", "
              << vert_pt[ee.at(i)][1] << ", " << vert_pt[ee.at(i)][2] << "];"
              << std::endl;
    std::cout << "int(" << i+1 << ",:) = [" << int_pt[ee.at(i)][0] << ", "
              << int_pt[ee.at(i)][1] << ", " << int_pt[ee.at(i)][2] << "];"
              << std::endl;
  }
}

void vd_glens::get_int_pos(std::vector<apf::MeshEntity*> ee, 
                           std::map<apf::MeshEntity*, apf::Vector3> &mer_pt,
                           std::map<apf::MeshEntity*, apf::Vector3> &int_pt,
                           std::map<apf::MeshEntity*, apf::Vector3> &vert_pt,
                           double r_cvx) {
  // Loop over non-collapsing edges
  // ent_id_map, id_ent_map, int_pt, vert_pt
  for (int i = 0; i < ee.size(); i++) {
    //std::cout << "e2_sp " << ee.at(i) << std::endl; 
    apf::Downward down;

    apf::MeshEntity* v_in;
    apf::MeshEntity* v_out;

    apf::Vector3 l(0,0,0);
    apf::Vector3 o(0,0,0);
    apf::Vector3 oc(0,0,0);

    m->getDownward(ee.at(i), 0, down);

    std::vector<apf::MeshEntity*>::iterator it;
    it = std::find (ent.at(0).begin(), ent.at(0).end(), down[0]);
    int i1 = std::distance(ent.at(0).begin(), it);

    if(it != ent.at(0).end()) {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);

      //assert(it == ent.at(0).end());
      if(it != ent.at(0).end()) {
        std::cout << "Edge shared by two merging vertices" << std::endl;
      }
      v_in = ent.at(0).at(i1);
      v_out = down[1];
    }
    else {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);
      i1 = std::distance(ent.at(0).begin(), it);

      assert(it != ent.at(0).end());
      v_in = ent.at(0).at(i1);
      v_out = down[0];
    }
    m->getPoint(v_in, 0, o);
    //o = midpoint;
    m->getPoint(v_out, 0, l);

    mer_pt[ee.at(i)] = o;
    vert_pt[ee.at(i)] = l;

    //std::cout << "c_pos " << c_pos << " "
    //          << " v_in " << o << " "
    //          << " v_out " << l
    //          << std::endl;

    if(r_cvx < - std::numeric_limits<double>::min()) {
      int_pt[ee.at(i)] = o/2+l/2;
    }
    else {
      l = l - o;
      l = norm_0(l);
      int_pt[ee.at(i)] = o+l*r_cvx;
    }

    std::cout << "vert(" << i+1 << ",:) = [" << vert_pt[ee.at(i)][0] << ", "
              << vert_pt[ee.at(i)][1] << ", " << vert_pt[ee.at(i)][2] << "];"
              << std::endl;
    std::cout << "int(" << i+1 << ",:) = [" << int_pt[ee.at(i)][0] << ", "
              << int_pt[ee.at(i)][1] << ", " << int_pt[ee.at(i)][2] << "];"
              << std::endl;
  }
}

// Split the entitiesby positioning the new vertex at the center.
void vd_glens::split_ctr(std::vector<apf::MeshEntity*> &ee,
                         std::vector<apf::MeshEntity*> &es_in,
                         std::vector<apf::MeshEntity*> &et,
              std::vector<apf::MeshEntity*> &v_sp_set, 
              std::map<apf::MeshEntity*, apf::Vector3> &int_pt, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_sp_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_mer_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_ori_dir, 
              double r_edge_min) {

  apf::Downward d_v;
  apf::Downward d_e;
  apf::Downward d_t;

  std::vector<apf::Vector3> verts(ee.size(), apf::Vector3(0,0,0));

  for (int i = 0; i < ee.size(); i++) {
    verts.at(i) = int_pt[ee.at(i)];
  }

  std::cout << "Splitting tets: " << std::endl;

  v_sp_set.reserve(et.size() + es_in.size() + ee.size());

  for (int i = 0; i < et.size(); i++) {
    if(e_map[et.at(i)] != 1) {
      int ent_type = m->getType(et.at(i));
      int d = m->typeDimension[ent_type];
      assert(d == 3); 
      apf::Vector3 o(0,0,0);
      apf::Vector3 l(0,0,0);

      m->getDownward(et.at(i), 0, d_v);
      m->getDownward(et.at(i), 1, d_e);
      m->getDownward(et.at(i), 2, d_t);
      std::cout << et.at(i)
                << std::endl;

      std::vector<apf::MeshEntity*>::iterator it;
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;

      apf::MeshEntity* v_other = NULL;

      int i1 = -1;
      int i2 = -1;
      //int i3 = -1;
      for (int j = 0; j < 4; j++) {
        it = std::find (ent.at(0).begin(), ent.at(0).end(), d_v[j]);
        //int i_curr = std::distance(ent.at(0).begin(), it);
        if(it != ent.at(0).end()) {
          o = o + vd_get_pos(m, d_v[j]);
          std::cout << "d_v[" << j << "] = " << d_v[j] << std::endl;

          if(!found2) {
            if(!found1) {
              found1 = true;
              i1 = j;
              std::cout << "\tfirst" << std::endl;
            }
            else {
              found2 = true;
              i2 = j;
              std::cout << "\tsecond" << std::endl;
            }
          }
          else {
            found3 = true;
            //i3 = j;
            std::cout << "\tthird" << std::endl;
          }
        }
        else
          v_other = d_v[j];
      }
      assert(found1);
      assert(v_other != NULL);

      //l = inctr_ent(m, et.at(i));
      if(found3) {
        //v.resize(3);
        //m->getPoint(d_v[i1], 0, v.at(0));
        //m->getPoint(d_v[i2], 0, v.at(1));
        //m->getPoint(d_v[i3], 0, v.at(2));
        //o = inctr_ent(&v);
        o = o/3;
        l = vd_get_pos(m, v_other);
      }
      else if(found2) {
        //v.resize(2);
        //m->getPoint(d_v[i1], 0, v.at(0));
        //m->getPoint(d_v[i2], 0, v.at(1));
        //o = inctr_ent(&v);
        o = o/2;
        l = vd_get_pos(m, d_e[lookup_edge_x_e[lookup_vert_edge_g[i1][i2]]]);
      }
      else {
        l = vd_get_pos(m, d_t[lookup_surf_opp_g[i1]]);
      }

      if(found2) {
        tet_split->load_tet(et.at(i));
        tet_split->split_tet();
        apf::MeshEntity* v_sp = tet_split->get_vert_ctr();

        v_sp_set.push_back(v_sp);
        v_ori_dir[v_sp] = norm_0(l - o);
        if(r_edge_min > std::numeric_limits<double>::min())
          v_sp_pos[v_sp] = o + v_ori_dir[v_sp]*r_edge_min;
        else
          v_sp_pos[v_sp] = l/2 + o/2;
        v_mer_pos[v_sp] = o;

        m->setPoint(v_sp, 0, v_sp_pos[v_sp]);
        assert(vd_chk_neg_vert(m, v_sp) == 0);

        if (f_calc != NULL) {
          f_calc->vd_att_fields(m, v_sp);
        }

      }
    }
  }
  vd_save_vtk_vert(m, &ent.at(0), "./output/Glens_tets");

  std::cout << "Splitting tris: " << std::endl;
  for (int i = 0; i < es_in.size(); i++) {
    int ent_type = m->getType(es_in.at(i));
    int d = m->typeDimension[ent_type];
    assert(d == 2); 
    if(e_map[es_in.at(i)] != 1) {
      apf::Vector3 o(0,0,0);
      apf::Vector3 l(0,0,0);

      m->getDownward(es_in.at(i), 0, d_v);
      std::cout << es_in.at(i)
                << std::endl;

      std::vector<apf::MeshEntity*>::iterator it;
      bool found1 = false;
      bool found2 = false;

      apf::MeshEntity* v_other = NULL;

      //l = inctr_ent(m, es_in.at(i));
      int i1 = -1;
      int i2 = -1;
      for (int j = 0; j < 3; j++) {
        it = std::find (ent.at(0).begin(), ent.at(0).end(), d_v[j]);
        //int i_curr = std::distance(ent.at(0).begin(), it);
        if(it != ent.at(0).end()) {
          o = o + vd_get_pos(m, d_v[j]);
          std::cout << "d_v[" << j << "] = " << d_v[j] << std::endl;

          if(!found1) {
            found1 = true;
            i1 = j;
            std::cout << "\tfirst" << std::endl;
          }
          else {
            found2 = true;
            i2 = j;
            std::cout << "\tsecond" << std::endl;
          }
        }
        else
          v_other = d_v[j];
      }
      assert(found1);
      assert(v_other != NULL);

      if(found2) {
        o = o/2;
        l = vd_get_pos(m, v_other);
      }
      else {
        l = vd_get_pos(m, d_v[(i1 + 1) % 3]);
        l = l/2 + vd_get_pos(m, d_v[(i1 + 2) % 3])/2;
      }
      if(found2) {

        bipy_split->load_tri(es_in.at(i));
        bipy_split->split_bipy();
        apf::MeshEntity* v_sp = bipy_split->get_vert_ctr();

        v_sp_set.push_back(v_sp);
        v_ori_dir[v_sp] = norm_0(l - o);
        if(r_edge_min > std::numeric_limits<double>::min())
          v_sp_pos[v_sp] = o + v_ori_dir[v_sp]*r_edge_min;
        else
          v_sp_pos[v_sp] = l/2 + o/2;
        v_mer_pos[v_sp] = o;

        m->setPoint(v_sp, 0, v_sp_pos[v_sp]);

        assert(vd_chk_neg_vert(m, v_sp) == 0);

        if (f_calc != NULL) {
          f_calc->vd_att_fields(m, v_sp);
        }
      }
    }
  }
  for (int i = 0; i < ee.size(); i++) {
    //std::cout << "e2_sp " << ee.at(i) << std::endl; 

    apf::Downward down;

    apf::MeshEntity* v_in;
    apf::MeshEntity* v_out;

    apf::Vector3 l(0,0,0);
    apf::Vector3 o(0,0,0);
    apf::Vector3 oc(0,0,0);

    m->getDownward(ee.at(i), 0, down);

    std::vector<apf::MeshEntity*>::iterator it;
    it = std::find (ent.at(0).begin(), ent.at(0).end(), down[0]);
    int i1 = std::distance(ent.at(0).begin(), it);

    if(it != ent.at(0).end()) {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);

      assert(it == ent.at(0).end());
      v_in = ent.at(0).at(i1);
      v_out = down[1];
    }
    else {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);
      i1 = std::distance(ent.at(0).begin(), it);

      assert(it != ent.at(0).end());
      v_in = ent.at(0).at(i1);
      v_out = down[0];
    }

    m->getPoint(v_in, 0, o);
    //o = midpoint;
    m->getPoint(v_out, 0, l);

    lens_split->load_edge(ee.at(i));
    lens_split->split_lens();
    v_in = lens_split->get_vert_ctr();

    v_sp_set.push_back(v_in);
    v_ori_dir[v_in] = norm_0(l - o);
    if(r_edge_min > std::numeric_limits<double>::min())
      v_sp_pos[v_in] = o + v_ori_dir[v_in]*r_edge_min;
    else
      v_sp_pos[v_in] = l/2 + o/2;

    v_mer_pos[v_in] = o;


    std::cout << "Inversion before repositioning " << std::endl;
    assert(vd_chk_neg_vert(m, v_in) == 0);

    if (f_calc != NULL) {
      f_calc->vd_att_fields(m, v_in);
    }

    //m->setPoint(v_in, 0, o+l*d);
    m->setPoint(v_in, 0, v_sp_pos[v_in]);
    //split_verts.at(i) = v_in;

    std::cout << "Inversion after repositioning " << std::endl;
    assert(vd_chk_neg_vert(m, v_in) == 0);

    if(vrfy_msh)
      m->verify();
  }

  vd_save_vtk_vert(m, &ent.at(0), "./output/Glens_edges");
}

// Test whether all edges are intersected by the convex hull generated.
bool vd_glens::test_cvx(int dim, int cell_id) {

  cell_col_dim = dim; 
  cell_col_id = cell_id;
  l_clear();

  ent.resize(4);
  v_list.resize(3);
  m_list.resize(3);

  //vd_find_ent_geom(m, &ent.at(cell_col_dim), cell_col_id, cell_col_dim);
  ent.at(cell_col_dim) = e_list->e.at(cell_col_dim).at(cell_col_id-1).at(cell_col_dim);

  for(int i = 0; i < ent.at(cell_col_dim).size(); i++) {
    e_map[ent.at(cell_col_dim).at(i)] = 1;
  }
  for(int dim = cell_col_dim; dim > 1; dim--) {
    vd_set_down(m, &ent.at(dim), &ent.at(dim-1));
    for(int i = 0; i < ent.at(dim-1).size(); i++) {
      e_map[ent.at(dim-1).at(i)] = 1;
    }
  }

  vd_set_down(m, &ent.at(1), &ent.at(0));
  ext_c0_vec.reserve(ent.at(0).size());

  apf::Vector3 zero(0,0,0);

  struct ent_conn e0;
  c_base->get_conn_dim_gmi(0, dim, cell_id, &e0);
  ext_c0.clear();
  for(int i = 0; i < e0.conn.size(); i++) {
    if(c_base->get_cell_ext_gmi(0, e0.conn.at(i)) ) {
      ext_cell = true;
      ext_c0.push_back(e0.conn.at(i));
    }
  }

  for(int i = 0; i < ent.at(0).size(); i++) {
    e_map[ent.at(0).at(i)] = 2;
    v_map[ent.at(0).at(i)] = true;
    apf::ModelEntity* mdl_curr = m->toModel(ent.at(0).at(i));
    if(m->getModelType(mdl_curr) == 0) {
      c0_flag = true;
      vert_ctr_em = mdl_curr;
      vert_ctr = ent.at(0).at(i);
      int c_id = m->getModelTag(mdl_curr);
      if(ext_cell and c_base->get_cell_ext_gmi(0, c_id)) {
        ext_c0_vec.push_back(zero);
        m->getPoint(ent.at(0).at(i), 0, ext_c0_vec.back());
      }
    }
  }
  if(c0_flag) {
    v_map[vert_ctr] = false;
    k_map[vert_ctr] = true;
  }

  if(ext_cell) {
    assert(ext_c0_vec.size() > 0);
    for(int i = 0; i < ext_c0_vec.size(); i++) {
      midpoint = midpoint + ext_c0_vec.at(i);
    }
    midpoint = midpoint/ext_c0_vec.size();

    ext_shell* e_sh;
    sh_min.dim = -1;
    if(upd_shell and f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
      e_sh = f_calc->get_e_sh();

      ent_conn e_conn;
      c_base->get_conn_dim_gmi(0, cell_col_dim, cell_col_id, &e_conn);

      for(int i = 0; i < e_conn.conn.size(); i++) {
        if(c_base->get_cell_ext_gmi(0, e_conn.conn.at(i))) {
          shell sh_0curr = e_sh->get_shell(0, e_conn.conn.at(i)-1);
          if(sh_min.dim == -1 or sh_0curr.dim < sh_min.dim)
            sh_min = sh_0curr;
        }
      }
      if(tag_0c != -1)
        assert(sh_min == e_sh->get_shell(0, tag_0c - 1));
    }

    if(upd_shell and f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
      assert(sh_min.dim > -1);
      midpoint = f_calc->get_e_sh()->find_int_pos(sh_min, midpoint);
    }
  }
  else
    midpoint = vd_get_center(m, &ent.at(0));
    //midpoint = vd_get_center_e(m, &ent.at(cell_col_dim));

  std::vector<apf::MeshEntity*> ee(0);
  std::vector<apf::MeshEntity*> he(0);
  std::vector<apf::MeshEntity*> es(0);
  std::vector<apf::MeshEntity*> es_in(0);
  std::vector<apf::MeshEntity*> et(0);

  vd_set_up(m, &ent.at(0), &ee);
  vd_set_up(m, &ee, &es_in);
  vd_set_up(m, &es_in, &et);
  vd_set_down(m, &et, &es);

  vd_remove_set(&es, &es_in);

  for (int i = 0; i < es.size(); i++) {
    apf::Downward down;
    m->getDownward(es.at(i), 0, down);
    for (int j = 0; j < 3; j++) {
      assert(findIn(&ent.at(0), ent.at(0).size(), down[j]) == -1);
    }
  }

  // Remove the collapsing edges.
  for (int i = 0; i < ent.at(1).size(); i++) {
    assert(e_map[ent.at(1).at(i)] == 1);

    std::vector<apf::MeshEntity*>::iterator it;
    it = std::find (ee.begin(), ee.end(), ent.at(1).at(i));
    assert(it != ee.end());
    ee.erase(it);
  }

  apf::Downward d_v;
  apf::Vector3 pos_in(0,0,0);
  apf::Vector3 pos_out(0,0,0);
  double r_edge_min = -1;

  // Consider external edges that are collapsing:
  std::map<apf::MeshEntity*, bool> e_col_map{};

  for (int i = 0; i < ee.size(); i++) {
    m->getDownward(ee.at(i), 0, d_v);
    int i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[0]);
    //bool merg = false;
    if(i1 == -1) {
      i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[1]);
      assert(i1 > -1);

      m->getPoint(d_v[0], 0, pos_out);
      m->getPoint(d_v[1], 0, pos_in);
    }
    else {
      i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[1]);
      if(i1 != -1) {
        //std::cout << "Cannot collapse, entity outside cell bounded by merging vertices." << std::endl;
        e_col_map[ee.at(i)] = true;
        std::cout << "Collapsing edge outside collapsing cell. Check for problems after preconditioning." << std::endl;
        //return false;
        //std::cout << "Entity outside cell bounded by merging vertices." << std::endl;
        //merg = true;
      }
      m->getPoint(d_v[1], 0, pos_out);
      m->getPoint(d_v[0], 0, pos_in);
    }
    if(!e_col_map[ee.at(i)]) {
      double meas = (pos_out-pos_in).getLength();
      //if(merg)
      //  meas = meas/2;
      if (r_edge_min < 0 or r_edge_min > meas) {
        r_edge_min = meas;
      }
    }
  }

  //for (int i = 0; i < ent.at(1).size(); i++) {
  //  double meas = vd_meas_ent(m, ent.at(1).at(i));
  //  if (r_edge_min < 0 or r_edge_min > meas) {
  //    r_edge_min = meas;
  //  }
  //}
  double r_cvx = r_edge_min/cvx_rat;

  int count = 0;
  for (int i = 0; i < ee.size()-count; i++) {
    if(e_col_map[ee.at(i)]) {
      apf::MeshEntity* e_temp = ee.at(i);
      ee.at(i) = ee.at(ee.size()-count-1);
      ee.at(ee.size()-count-1) = e_temp;
      count = count + 1;
      i = i - 1;
    }
  }
  ee.resize(ee.size()-count);

  std::map<apf::MeshEntity*, apf::Vector3> int_pt{};
  std::map<apf::MeshEntity*, apf::Vector3> vert_pt{};
  std::map<apf::MeshEntity*, apf::Vector3> mer_pt{};
  get_int_pos(ee, mer_pt, int_pt, vert_pt, r_cvx);

  std::vector<apf::Vector3> verts(ee.size(), apf::Vector3(0,0,0));
  // Include the exterior vertices to the convex hull, as well.
  if(ext_cell) {
    verts.resize(ee.size() + ext_c0_vec.size());
  }
  else {
    verts.resize(ee.size());
  }

  for (int i = 0; i < ee.size(); i++) {
    if(!e_col_map[ee.at(i)])
      verts.at(i) = int_pt[ee.at(i)];
  }

  if(ext_cell) {
    for (int i = 0; i < ext_c0_vec.size(); i++) {
      verts.at(ee.size() + i) = ext_c0_vec.at(i);
    }
  }

  std::vector<int> hull_id(0);
  std::vector<std::vector<int> > f_id(0, std::vector<int>(0));
  std::vector<apf::Vector3 > f_norm(0, apf::Vector3(0,0,0));
  std::vector<std::vector<int> > e_id(0, std::vector<int>(0));
  cvx_hull chull;
  hull_id = chull.find_chull_vertices(&verts);
  std::sort(hull_id.begin(), hull_id.end());
  chull.get_face_id(&f_id);
  chull.get_edge_id(&e_id);

  f_norm.resize(f_id.size());
  for(int i = 0; i < f_id.size(); i++) {
    f_norm.at(i) = vd_area_out_n(verts.at(f_id.at(i).at(0)), 
                                 verts.at(f_id.at(i).at(1)),
                                 verts.at(f_id.at(i).at(2)));
  }

  vd_plane pl;
  vd_inter v_int;

  for(int i = 0; i < ee.size(); i++) {
    bool found = false;
    if(!e_col_map[ee.at(i)]) {

      for(int j = 0; j < f_id.size(); j++) {
        pl.pos = verts.at(f_id.at(j).at(0));
        pl.norm = f_norm.at(j);
        pl_int_2pts(m, mer_pt[ee.at(i)], vert_pt[ee.at(i)], &pl, &v_int);

        if(v_int.cut == 1) {
          found = true;
          j = f_id.size();
        }
      }
      if(!found) {
        std::cout << dim << "c" << cell_id 
                  << ": No face of the convex hull intersects the edge. "
                  << std::endl;
        return false;
      }
    }
  }
  return true;
}

void vd_glens::split_cvx(std::vector<apf::MeshEntity*> &ee,
                         std::vector<apf::MeshEntity*> &es_in,
                         std::vector<apf::MeshEntity*> &et,
              std::vector<apf::MeshEntity*> &v_sp_set, 
              std::map<apf::MeshEntity*, apf::Vector3> &int_pt, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_sp_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_mer_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_ori_dir) {

  apf::Downward d_v;
  apf::Downward d_e;
  apf::Downward d_t;

  std::vector<apf::Vector3> verts(ee.size(), apf::Vector3(0,0,0));
  // Include the exterior vertices to the convex hull, as well.
  if(ext_cell) {
    verts.resize(ee.size() + ext_c0_vec.size());
  }
  else {
    verts.resize(ee.size());
  }

  for (int i = 0; i < ee.size(); i++) {
    verts.at(i) = int_pt[ee.at(i)];
  }

  if(ext_cell) {
    for (int i = 0; i < ext_c0_vec.size(); i++) {
      verts.at(ee.size() + i) = ext_c0_vec.at(i);
      std::cout << "int(" << ee.size()+i << ",:) = ["
                << ext_c0_vec.at(i)[0] << ","
                << ext_c0_vec.at(i)[1] << ","
                << ext_c0_vec.at(i)[2] << "];"
                << std::endl;
    }
  }

  std::vector<int> hull_id(0);
  std::vector<std::vector<int> > f_id(0, std::vector<int>(0));
  std::vector<std::vector<int> > e_id(0, std::vector<int>(0));
  cvx_hull chull;
  hull_id = chull.find_chull_vertices(&verts);
  std::sort(hull_id.begin(), hull_id.end());
  chull.get_face_id(&f_id);
  chull.get_edge_id(&e_id);

  //std::cout << "Splitting tets: " << std::endl;

  v_sp_set.reserve(et.size() + es_in.size() + ee.size());

  for (int i = 0; i < ee.size(); i++) {
    //std::cout << "e2_sp " << ee.at(i) << std::endl; 

    apf::Downward down;

    apf::MeshEntity* v_in;
    apf::MeshEntity* v_out;

    apf::Vector3 l(0,0,0);
    apf::Vector3 o(0,0,0);
    apf::Vector3 oc(0,0,0);

    m->getDownward(ee.at(i), 0, down);

    std::vector<apf::MeshEntity*>::iterator it;
    it = std::find (ent.at(0).begin(), ent.at(0).end(), down[0]);
    int i1 = std::distance(ent.at(0).begin(), it);

    if(it != ent.at(0).end()) {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);

      assert(it == ent.at(0).end());
      v_in = ent.at(0).at(i1);
      v_out = down[1];
    }
    else {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);
      i1 = std::distance(ent.at(0).begin(), it);

      assert(it != ent.at(0).end());
      v_in = ent.at(0).at(i1);
      v_out = down[0];
    }

    m->getPoint(v_in, 0, o);
    //o = midpoint;
    m->getPoint(v_out, 0, l);

    apf::Vector3 proj_pt(0, 0, 0);
    int f_j = -1;
    //if((it = std::find (he.begin(), he.end(), ee.at(i))) != he.end()) {
    //  int i_curr = std::distance(he.begin(), it);
    //  proj_pt = verts.at(i_curr);
    //}
    // Either find the intersection with the hull, or the end vertex of the
    // edge is already in the hull.
    bool found_cvx = false;
    bool skip_cvx = false;

    if(std::find (hull_id.begin(), hull_id.end(), i) != hull_id.end()) {
      skip_cvx = true;
      proj_pt = verts.at(i);
    }
    else {
      for (int j = 0; j < f_id.size(); j++) {
        int ii0 = f_id.at(j).at(0);
        int ii1 = f_id.at(j).at(1);
        int ii2 = f_id.at(j).at(2);

        // End of the edge is on cvx hull.
        if(pt_int_tri(verts.at(ii0), verts.at(ii1), verts.at(ii2), 
                      verts.at(i), 0.01)) {
          found_cvx = true;
          skip_cvx = true;
          f_j = j+1;
          j = f_id.size();
          proj_pt = verts.at(i);
        }
        else if(tri_int_edge(verts.at(ii0), verts.at(ii1), 
                        verts.at(ii2), o, l, -10e-4, &proj_pt)) {
          found_cvx = true;
          f_j = j+1;
          j = f_id.size();
        }
      }
      //assert(found_cvx);
      if(!found_cvx) {
        proj_pt = verts.at(i);
      }
    }
    double t0 = PCU_Time();

    lens_split->load_edge(ee.at(i));
    lens_split->split_lens();
    v_in = lens_split->get_vert_ctr();

    if (f_calc != NULL) {
      f_calc->vd_att_fields(m, v_in);
    }

    double t1 = PCU_Time();
    std::cout << "Split took " << t1 - t0 << "s." << std::endl;

    apf::setScalar(ch_f, v_in, 0, f_j);
    apf::setScalar(sp_f, v_in, 0, 1);

    v_sp_set.push_back(v_in);
    v_ori_dir[v_in] = norm_0(l - o);
    // If the proj_pt is closer to the vertex outside by 0.1 of the distance
    // between the vertices, position the split vertex at 0.9*l + 0.1*o.
    double l1 = (l - proj_pt).getLength();
    double l2 = (l - o).getLength()*0.1;
    if(l1 < l2)
      proj_pt = l*0.9 + o*0.1;
    v_sp_pos[v_in] = proj_pt;
    v_mer_pos[v_in] = o;

    std::cout << "on the hull: " << skip_cvx << std::endl;
    std::cout << "Inversion before repositioning " << std::endl;
    assert(vd_chk_neg_vert(m, v_in) == 0);

    f_calc->vd_att_fields(m, v_in);

    //m->setPoint(v_in, 0, o+l*d);
    m->setPoint(v_in, 0, proj_pt);
    //split_verts.at(i) = v_in;

    std::cout << "Inversion after repositioning " << std::endl;
    assert(vd_chk_neg_vert(m, v_in) == 0);

    if(vrfy_msh)
      m->verify();
  }

  vd_save_vtk_vert(m, &ent.at(0), "./output/Glens_edges");

}

// Imposing r_cvx as the smallest edge length everywhere is too constrictive.
// Instead, take r_cvx as a fraction of longest edge length. If an edge is shorter
// than this, use the half of the original edge length.
void vd_glens::split_noncvx2(std::vector<apf::MeshEntity*> &ee,
                         std::vector<apf::MeshEntity*> &es_in,
                         std::vector<apf::MeshEntity*> &et,
              std::vector<apf::MeshEntity*> &v_sp_set, 
              std::map<apf::MeshEntity*, apf::Vector3> &int_pt, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_sp_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_mer_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_ori_dir, double r_cvx) {

  apf::Downward d_v;
  apf::Downward d_e;
  apf::Downward d_t;

  std::vector<apf::Vector3> verts(ee.size(), apf::Vector3(0,0,0));
  // Include the exterior vertices to the convex hull, as well.
  if(ext_cell) {
    verts.resize(ee.size() + ext_c0_vec.size());
  }
  else {
    verts.resize(ee.size());
  }

  for (int i = 0; i < ee.size(); i++) {
    verts.at(i) = int_pt[ee.at(i)];
  }

  if(ext_cell) {
    for (int i = 0; i < ext_c0_vec.size(); i++) {
      verts.at(ee.size() + i) = ext_c0_vec.at(i);
      std::cout << "int(" << ee.size()+i << ",:) = ["
                << ext_c0_vec.at(i)[0] << ","
                << ext_c0_vec.at(i)[1] << ","
                << ext_c0_vec.at(i)[2] << "];"
                << std::endl;
    }
  }

  //std::cout << "Splitting tets: " << std::endl;

  v_sp_set.reserve(et.size() + es_in.size() + ee.size());

  for (int i = 0; i < ee.size(); i++) {
    //std::cout << "e2_sp " << ee.at(i) << std::endl; 

    apf::Downward down;

    apf::MeshEntity* v_in;
    apf::MeshEntity* v_out;

    apf::Vector3 l(0,0,0);
    apf::Vector3 o(0,0,0);
    apf::Vector3 oc(0,0,0);

    m->getDownward(ee.at(i), 0, down);

    std::vector<apf::MeshEntity*>::iterator it;
    it = std::find (ent.at(0).begin(), ent.at(0).end(), down[0]);
    int i1 = std::distance(ent.at(0).begin(), it);

    if(it != ent.at(0).end()) {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);

      assert(it == ent.at(0).end());
      v_in = ent.at(0).at(i1);
      v_out = down[1];
    }
    else {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);
      i1 = std::distance(ent.at(0).begin(), it);

      assert(it != ent.at(0).end());
      v_in = ent.at(0).at(i1);
      v_out = down[0];
    }

    m->getPoint(v_in, 0, o);
    //o = midpoint;
    m->getPoint(v_out, 0, l);

    double t0 = PCU_Time();
    lens_split->load_edge(ee.at(i));
    lens_split->split_lens();
    v_in = lens_split->get_vert_ctr();

    if (f_calc != NULL) {
      f_calc->vd_att_fields(m, v_in);
    }

    double t1 = PCU_Time();
    std::cout << "Split took " << t1 - t0 << "s." << std::endl;

    v_sp_set.push_back(v_in);
    v_ori_dir[v_in] = l - o;
    double split_len = v_ori_dir[v_in].getLength()/2;
    v_ori_dir[v_in] = norm_0(v_ori_dir[v_in]);
    apf::Vector3 proj_pt(0, 0, 0);
    if(split_len < r_cvx)
      proj_pt = o + v_ori_dir[v_in]*split_len;
    else {
      proj_pt = o + v_ori_dir[v_in]*r_cvx;
      // If the proj_pt is closer to the vertex outside by 0.1 of the distance
      // between the vertices, position the split vertex at 0.9*l + 0.1*o.
      double l1 = (l - proj_pt).getLength();
      double l2 = (l - o).getLength()*0.1;
      if(l1 < l2)
        proj_pt = l*0.9 + o*0.1;
    }
    v_sp_pos[v_in] = proj_pt;
    v_mer_pos[v_in] = o;

    std::cout << "Inversion before repositioning " << std::endl;
    assert(vd_chk_neg_vert(m, v_in) == 0);

    f_calc->vd_att_fields(m, v_in);

    //m->setPoint(v_in, 0, o+l*d);
    m->setPoint(v_in, 0, proj_pt);
    //split_verts.at(i) = v_in;

    std::cout << "Inversion after repositioning " << std::endl;
    assert(vd_chk_neg_vert(m, v_in) == 0);

    if(vrfy_msh)
      m->verify();
  }

  vd_save_vtk_vert(m, &ent.at(0), "./output/Glens_edges");

}

void vd_glens::split_noncvx(std::vector<apf::MeshEntity*> &ee,
                         std::vector<apf::MeshEntity*> &es_in,
                         std::vector<apf::MeshEntity*> &et,
              std::vector<apf::MeshEntity*> &v_sp_set, 
              std::map<apf::MeshEntity*, apf::Vector3> &int_pt, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_sp_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_mer_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_ori_dir, double r_cvx) {

  apf::Downward d_v;
  apf::Downward d_e;
  apf::Downward d_t;

  std::vector<apf::Vector3> verts(ee.size(), apf::Vector3(0,0,0));
  // Include the exterior vertices to the convex hull, as well.
  if(ext_cell) {
    verts.resize(ee.size() + ext_c0_vec.size());
  }
  else {
    verts.resize(ee.size());
  }

  for (int i = 0; i < ee.size(); i++) {
    verts.at(i) = int_pt[ee.at(i)];
  }

  if(ext_cell) {
    for (int i = 0; i < ext_c0_vec.size(); i++) {
      verts.at(ee.size() + i) = ext_c0_vec.at(i);
      std::cout << "int(" << ee.size()+i << ",:) = ["
                << ext_c0_vec.at(i)[0] << ","
                << ext_c0_vec.at(i)[1] << ","
                << ext_c0_vec.at(i)[2] << "];"
                << std::endl;
    }
  }

  //std::cout << "Splitting tets: " << std::endl;

  v_sp_set.reserve(et.size() + es_in.size() + ee.size());

  for (int i = 0; i < ee.size(); i++) {
    //std::cout << "e2_sp " << ee.at(i) << std::endl; 

    apf::Downward down;

    apf::MeshEntity* v_in;
    apf::MeshEntity* v_out;

    apf::Vector3 l(0,0,0);
    apf::Vector3 o(0,0,0);
    apf::Vector3 oc(0,0,0);

    m->getDownward(ee.at(i), 0, down);

    std::vector<apf::MeshEntity*>::iterator it;
    it = std::find (ent.at(0).begin(), ent.at(0).end(), down[0]);
    int i1 = std::distance(ent.at(0).begin(), it);

    if(it != ent.at(0).end()) {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);

      assert(it == ent.at(0).end());
      v_in = ent.at(0).at(i1);
      v_out = down[1];
    }
    else {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), down[1]);
      i1 = std::distance(ent.at(0).begin(), it);

      assert(it != ent.at(0).end());
      v_in = ent.at(0).at(i1);
      v_out = down[0];
    }

    m->getPoint(v_in, 0, o);
    //o = midpoint;
    m->getPoint(v_out, 0, l);

    double t0 = PCU_Time();
    lens_split->load_edge(ee.at(i));
    lens_split->split_lens();
    v_in = lens_split->get_vert_ctr();

    if (f_calc != NULL) {
      f_calc->vd_att_fields(m, v_in);
    }

    double t1 = PCU_Time();
    std::cout << "Split took " << t1 - t0 << "s." << std::endl;

    v_sp_set.push_back(v_in);
    v_ori_dir[v_in] = norm_0(l - o);
    apf::Vector3 proj_pt(0, 0, 0);
    proj_pt = o + v_ori_dir[v_in]*r_cvx;

    // If the proj_pt is closer to the vertex outside by 0.1 of the distance
    // between the vertices, position the split vertex at 0.9*l + 0.1*o.
    double l1 = (l - proj_pt).getLength();
    double l2 = (l - o).getLength()*0.1;
    if(l1 < l2)
      proj_pt = l*0.9 + o*0.1;
    v_sp_pos[v_in] = proj_pt;
    v_mer_pos[v_in] = o;

    std::cout << "Inversion before repositioning " << std::endl;
    assert(vd_chk_neg_vert(m, v_in) == 0);

    f_calc->vd_att_fields(m, v_in);

    //m->setPoint(v_in, 0, o+l*d);
    m->setPoint(v_in, 0, proj_pt);
    //split_verts.at(i) = v_in;

    std::cout << "Inversion after repositioning " << std::endl;
    assert(vd_chk_neg_vert(m, v_in) == 0);

    if(vrfy_msh)
      m->verify();
  }

  vd_save_vtk_vert(m, &ent.at(0), "./output/Glens_edges");

}

// Check if all surrounding tris are non-merging during collapse. If not, 
// preconditioning is necessary.
// Mark each surrounding tri by the edge across. For merging tris sharing the same
// edge, in the tri list sharing the edge, check if there is at least one tri
// that they share a collapsing tet with. 
// Also, merging entities should cause no additional topological transitions.
// For merging edges, find the lowest stratum membership. Check that it is either
// adjacent or same. 
bool vd_glens::chk_merg() {

  std::vector<std::vector<apf::MeshEntity*> > ee(4, 
                                    std::vector<apf::MeshEntity*>(0));
  ee.at(cell_col_dim) = e_list->e.at(cell_col_dim).at(cell_col_id-1).at(cell_col_dim);
  // Map of surrounding triangles to the edge across the merging vertex.
  std::map<apf::MeshEntity*, bool> t_surr{};
  std::map<apf::MeshEntity*, bool> t_merg{};
  std::map<apf::MeshEntity*, bool> t_col{};
  std::map<apf::MeshEntity*, apf::MeshEntity*> t_ex{};
  // Map of triangle pairs connected over collapsing tets.
  //std::map<std::pair<apf::MeshEntity*, apf::MeshEntity*>, bool> t_pair{};

  for(int dim = cell_col_dim; dim > 0; dim--)
    vd_set_down(m, &ee.at(dim), &ee.at(dim-1));
  vd_set_up(m, &ee.at(0), &ee.at(1));
  std::vector<apf::MeshEntity*> es(0);
  std::vector<apf::MeshEntity*> et(0);
  vd_set_up(m, &ee.at(1), &es);
  vd_set_up(m, &es, &et);

  apf::Downward d_v;
  apf::Downward d_e;
  for(int i = 0; i < es.size(); i++) {
    m->getDownward(es.at(i), 0, d_v);
    int j1 = -1;
    bool found1 = false;
    bool found2 = false;
    for(int j = 0; j < 3; j++) {
      if(e_map[d_v[j]] == 2) {
        if(!found1) {
          found1 = true;
          j1 = j;
        }
        else {
          if(!found2) {
            found2 = true;
            j = 3;
            t_col[es.at(i)] = true;
          }
        }
      }
    }
    if(!found2) {
      assert(found1);
      m->getDownward(es.at(i), 1, d_e);
      t_ex[es.at(i)] = d_e[lookup_tri_edge_x_g[j1]];
      t_surr[es.at(i)] = true;
    }
  }
  // For collapsing tets with two merging vertices mark the merging triangle 
  // couples. For surrounding triangles, count the number of times an edge across
  // a vertex is encountered. If it is more than once, the triangle should be a
  // merging triangle. Otherwise the collapse is not valid without 
  // preconditioning.
  for(int i = 0; i < et.size(); i++) {
    m->getDownward(et.at(i), 0, d_v);
    int j1 = -1;
    int j2 = -1;
    bool found1 = false;
    bool found2 = false;
    bool found3 = false;
    for(int j = 0; j < 4; j++) {
      if(e_map[d_v[j]] == 2) {
        if(!found1) {
          found1 = true;
          j1 = j;
        }
        else {
          if(!found2) {
            found2 = true;
            j2 = j;
          }
          else {
            found3 = true;
            j = 4;
          }
        }
      }
    }
    if(!found3 and found2) {
      m->getDownward(et.at(i), 2, d_e);
      int t1 = lookup_surf_n_gg[j1][j2][0];
      int t2 = lookup_surf_n_gg[j1][j2][0];
      //t_pair[std::make_pair<d_e[t1], d_e[t2]>] = true;
      //t_pair[std::make_pair<d_e[t2], d_e[t1]>] = true;
      t_merg[d_e[t1]] = true;
      t_merg[d_e[t2]] = true;
    }
  }

  // collect the counter edge lists of merging tris.
  // for a counter edge, the tris should have at least one merging tri in the 
  // list. 

  int e_count = 0;
  std::map<apf::MeshEntity*, int> e_t_count{};
  std::map<apf::MeshEntity*, bool> e_added{};
  std::map<apf::MeshEntity*, int> e_id{};
  std::vector<std::vector<apf::MeshEntity*> > e_t_list(0, 
                                      std::vector<apf::MeshEntity*>(0));
  for(int i = 0; i < es.size(); i++) {
    if(!t_col[es.at(i)] and t_surr[es.at(i)]) {
      apf::MeshEntity* e_x = t_ex[es.at(i)];

      if(e_t_count[e_x] == 0)
        e_count = e_count + 1;
      e_t_count[e_x] = e_t_count[e_x] + 1;
    }
  }
  for(int i = 0; i < es.size(); i++) {
    if(!t_col[es.at(i)] and t_surr[es.at(i)]) {
      apf::MeshEntity* e_x = t_ex[es.at(i)];
      if(e_t_count[e_x] > 1) {
        if(!t_merg[es.at(i)]) {
          return true;
        }
      }
    }
  }

  // Also, for non-collapsing edges check for disallowed merges.
  e_t_list.reserve(e_count);
  for(int i = 0; i < es.size(); i++) {
    if(!t_col[es.at(i)] and t_surr[es.at(i)]) {
      apf::MeshEntity* e_x = t_ex[es.at(i)];
      if(!e_added[e_x]) {
        e_id[e_x] = e_t_list.size();
        e_t_list.push_back(std::vector<apf::MeshEntity*>(0));
        e_t_list.back().reserve(e_t_count[e_x]);
      }
      e_added[e_x] = true;
    }
  }
  for(int i = 0; i < es.size(); i++) {
    if(!t_col[es.at(i)] and t_surr[es.at(i)]) {
      apf::MeshEntity* e_x = t_ex[es.at(i)];
      int e1 = e_id[e_x];
      e_t_list.at(e1).push_back(es.at(i));
    }
  }
  for(int i = 0; i < e_t_list.size(); i++) {
    int d_min = 4;
    int id_min = 0;
    for(int j = 0; j < e_t_list.at(i).size(); j++) {
      apf::ModelEntity* mdl1 = m->toModel(e_t_list.at(i).at(j));
      int d1 = m->getModelType(mdl1);
      int id1 = m->getModelTag(mdl1);
      if(d1 < d_min) {
        d_min = d1;
        id_min = id1;
      }
    }
    for(int j = 0; j < e_t_list.at(i).size(); j++) {
      apf::ModelEntity* mdl1 = m->toModel(e_t_list.at(i).at(j));
      int d1 = m->getModelType(mdl1);
      int id1 = m->getModelTag(mdl1);
      if(d1 == d_min) {
        if(id1 != id_min)
          return true;
      }
      else {
        if(!c_base->chk_conn_d_gmi(d_min, id_min, d1, id1))
          return true;
      }
    }
  }
  return false;
}

// Check if all surrounding tets are non-inverting during collapse. If so, no 
// preconditioning necessary.
bool vd_glens::chk_precond(int dim, int cell_id) {

  precond_flag = false;

  cell_col_dim = dim; 
  cell_col_id = cell_id;
  l_clear();

  ent.resize(4);
  v_list.resize(3);
  m_list.resize(3);

  //vd_find_ent_geom(m, &ent.at(cell_col_dim), cell_col_id, cell_col_dim);
  ent.at(cell_col_dim) = e_list->e.at(cell_col_dim).at(cell_col_id-1).at(cell_col_dim);

  for(int i = 0; i < ent.at(cell_col_dim).size(); i++) {
    e_map[ent.at(cell_col_dim).at(i)] = 1;
  }
  for(int dim = cell_col_dim; dim > 1; dim--) {
    vd_set_down(m, &ent.at(dim), &ent.at(dim-1));
    for(int i = 0; i < ent.at(dim-1).size(); i++) {
      e_map[ent.at(dim-1).at(i)] = 1;
    }
  }

  struct ent_conn e0;
  c_base->get_conn_dim_gmi(0, dim, cell_id, &e0);
  ext_c0.clear();
  for(int i = 0; i < e0.conn.size(); i++) {
    if(c_base->get_cell_ext_gmi(0, e0.conn.at(i)) ) {
      ext_cell = true;
      ext_c0.push_back(e0.conn.at(i));
    }
  }

  vd_set_down(m, &ent.at(1), &ent.at(0));
  ext_c0_vec.reserve(ent.at(0).size());

  apf::Vector3 zero(0,0,0);

  for(int i = 0; i < ent.at(0).size(); i++) {
    e_map[ent.at(0).at(i)] = 2;
    v_map[ent.at(0).at(i)] = true;
    apf::ModelEntity* mdl_curr = m->toModel(ent.at(0).at(i));
    if(m->getModelType(mdl_curr) == 0) {
      c0_flag = true;
      vert_ctr_em = mdl_curr;
      vert_ctr = ent.at(0).at(i);
      int c_id = m->getModelTag(mdl_curr);
      if(ext_cell and c_base->get_cell_ext_gmi(0, c_id)) {
        ext_c0_vec.push_back(zero);
        m->getPoint(ent.at(0).at(i), 0, ext_c0_vec.back());
      }
    }
  }
  if(c0_flag) {
    v_map[vert_ctr] = false;
    k_map[vert_ctr] = true;
  }

  if(ext_cell) {
    assert(ext_c0_vec.size() > 0);
    for(int i = 0; i < ext_c0_vec.size(); i++) {
      midpoint = midpoint + ext_c0_vec.at(i);
    }
    midpoint = midpoint/ext_c0_vec.size();

    ext_shell* e_sh;
    sh_min.dim = -1;
    if(upd_shell and f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
      e_sh = f_calc->get_e_sh();
      if(ext_cell) {
        ent_conn e_conn;
        c_base->get_conn_dim_gmi(0, cell_col_dim, cell_col_id, &e_conn);

        for(int i = 0; i < e_conn.conn.size(); i++) {
          if(c_base->get_cell_ext_gmi(0, e_conn.conn.at(i))) {
            shell sh_0curr = e_sh->get_shell(0, e_conn.conn.at(i)-1);
            if(sh_min.dim == -1 or sh_0curr.dim < sh_min.dim)
              sh_min = sh_0curr;
          }
        }
        if(tag_0c != -1)
          assert(sh_min == e_sh->get_shell(0, tag_0c - 1));
      }
    }

    if(upd_shell and f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
      assert(sh_min.dim > -1);
      midpoint = f_calc->get_e_sh()->find_int_pos(sh_min, midpoint);
    }
  }
  else
    midpoint = vd_get_center(m, &ent.at(0));
    //midpoint = vd_get_center_e(m, &ent.at(cell_col_dim));

  // Collect the non-collapsing tetrahedra and check for inversion.

  apf::Downward d_v;
  apf::Downward d_s;
  std::vector<std::vector<apf::MeshEntity*> > ee(4, 
                                    std::vector<apf::MeshEntity*>(0));
  ee.at(cell_col_dim) = e_list->e.at(cell_col_dim).at(cell_col_id-1).at(cell_col_dim);
  for(int dim = cell_col_dim; dim > 0; dim--)
    vd_set_down(m, &ee.at(dim), &ee.at(dim-1));
  for(int dim = 0; dim < 3; dim++)
    vd_set_up(m, &ee.at(dim), &ee.at(dim+1));

  std::map<apf::MeshEntity*, bool> tet_skip{};
  for (int i = 0; i < ee.at(3).size(); i++) {
    m->getDownward(ee.at(3).at(i), 0, d_v);

    std::vector<apf::MeshEntity*>::iterator it;
    bool found1 = false;
    bool found2 = false;
    bool found3 = false;

    int i1 = -1;
    for (int j = 0; j < 4; j++) {
      it = std::find (ent.at(0).begin(), ent.at(0).end(), d_v[j]);
      //int i_curr = std::distance(ent.at(0).begin(), it);
      if(it != ent.at(0).end()) {
        //std::cout << "d_v[" << j << "] = " << d_v[j] << std::endl;

        if(!found2) {
          if(!found1) {
            found1 = true;
            i1 = j;
          }
          else {
            found2 = true;
          }
        }
        else {
          found3 = true;
        }
      }
    }
    if(found2)
      tet_skip[ee.at(3).at(i)] = true;
  }

  std::map<apf::MeshEntity*, apf::Vector3> v_pos{};
  apf::Vector3 pos(0,0,0);
  apf::Vector3 n_area(0,0,0);
  for (int i = 0; i < ee.at(0).size(); i++) {
    m->getPoint(ee.at(0).at(i), 0, pos);
    v_pos[ee.at(0).at(i)] = pos;
  }

  for (int i = 0; i < ee.at(3).size(); i++) {
    if(!tet_skip[ee.at(3).at(i)]) {

      m->getDownward(ee.at(3).at(i), 0, d_v);
      m->getDownward(ee.at(3).at(i), 2, d_s);

      std::vector<apf::MeshEntity*>::iterator it;
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;

      int i1 = -1;
      for (int j = 0; j < 4; j++) {
        it = std::find (ent.at(0).begin(), ent.at(0).end(), d_v[j]);
        //int i_curr = std::distance(ent.at(0).begin(), it);
        if(it != ent.at(0).end()) {
          //std::cout << "d_v[" << j << "] = " << d_v[j] << std::endl;

          if(!found2) {
            if(!found1) {
              found1 = true;
              i1 = j;
            }
            else {
              found2 = true;
            }
          }
          else {
            found3 = true;
          }
        }
      }
      assert(!found2 and found1);
      int t1 = lookup_surf_opp_g[i1];
      pos = vd_get_pos(m, d_s[t1]);
      n_area = norm_0(vd_area_out_n(m, d_s[t1]));
      double side1 = n_area*(v_pos[d_v[i1]] - pos);
      double side2 = n_area*(midpoint - pos);
      if(!(side1*side2 > std::numeric_limits<double>::min())) {
        precond_flag = true;
        i = ee.at(3).size();
      }
    }
  }
  if(!precond_flag) {
    precond_flag = precond_flag or chk_merg();
  }
  return precond_flag;
}

// Check and correct the lens. First fixes the merging lens entities, so that
// no disallowed topological changes are taking place. Next, for each edge
// radiating from the merging vertices, estimate the position of the split 
// vertices. Relax the positions of these vertices, such that no element is 
// inverting.

// Two rectifiers, one for volume one for edge length. The minimum distance
// from the actual vertex positions can at most be 1/10th of the original edge 
// length. The minimum volume should be 1/1000th of the original volume. 
bool vd_glens::corr_lens_new() {

  vd_save_vtk_vert(m, &ent.at(0), "./output/Glens_before");

  midpoint = apf::Vector3(0,0,0);

  apf::Downward d_v;
  apf::Downward d_e;
  apf::Downward d_t;
  // Fix the edges bounded by the merging vertices.
  // Collect the edges.
  std::vector<apf::MeshEntity*> ee(0);
  std::vector<apf::MeshEntity*> he(0);
  std::vector<apf::MeshEntity*> es(0);
  std::vector<apf::MeshEntity*> es_in(0);
  std::vector<apf::MeshEntity*> et(0);

  vd_set_up(m, &ent.at(0), &ee);
  vd_set_up(m, &ee, &es_in);
  vd_set_up(m, &es_in, &et);
  vd_set_down(m, &et, &es);

  vd_remove_set(&es, &es_in);

  for (int i = 0; i < es.size(); i++) {
    apf::Downward down;
    m->getDownward(es.at(i), 0, down);
    for (int j = 0; j < 3; j++) {
      assert(findIn(&ent.at(0), ent.at(0).size(), down[j]) == -1);             
    }
  }

  // Remove the collapsing edges.
  for (int i = 0; i < ent.at(1).size(); i++) {
    assert(e_map[ent.at(1).at(i)] == 1);

    std::vector<apf::MeshEntity*>::iterator it;
    it = std::find (ee.begin(), ee.end(), ent.at(1).at(i));
    assert(it != ee.end());
    ee.erase(it);
  }
  std::map<apf::MeshEntity*, bool> e_col_map{};

  for (int i = 0; i < ee.size(); i++) {
    m->getDownward(ee.at(i), 0, d_v);
    int i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[0]);
    //bool merg = false;
    if(i1 == -1) {
      i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[1]);
      assert(i1 > -1);
    }
    else {
      i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[1]);
      if(i1 != -1) {
        //std::cout << "Cannot collapse, entity outside cell bounded by merging vertices." << std::endl;
        e_col_map[ee.at(i)] = true;
        std::cout << "Collapsing edge outside collapsing cell. Check for problems after preconditioning." << std::endl;
        //return false;
        //std::cout << "Entity outside cell bounded by merging vertices." << std::endl;
        //merg = true;
      }
    }
  }
  // Consider exterior edges that are collapsing:
  ent.at(1).reserve(ent.at(1).size() + e_col_ext.size());
  for(int i = 0; i < ee.size(); i++) {
    if(e_map[ee.at(i)] == 1)
      ent.at(1).push_back(ee.at(i));
  }

  int count = 0;
  for (int i = 0; i < ee.size()-count; i++) {
    if(e_col_map[ee.at(i)]) {
      apf::MeshEntity* e_temp = ee.at(i);
      ee.at(i) = ee.at(ee.size()-count-1);
      ee.at(ee.size()-count-1) = e_temp;
      count = count + 1;
      i = i - 1;
    }
  }
  ee.resize(ee.size()-count);

  if(ext_cell) {
    assert(ext_c0_vec.size() > 0);
    for(int i = 0; i < ext_c0_vec.size(); i++) {
      midpoint = midpoint + ext_c0_vec.at(i);
    }
    midpoint = midpoint/ext_c0_vec.size();
    if(upd_shell and f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
      assert(sh_min.dim > -1);
      midpoint = f_calc->get_e_sh()->find_int_pos(sh_min, midpoint);
    }
  }
  else
    midpoint = vd_get_center_e(m, &ent.at(cell_col_dim));
    //midpoint = vd_get_center_e(m, &ent.at(cell_col_dim));

  std::cout << "c_pos " << midpoint << std::endl;

  // Find the maximum of the distance to the merging vertices.
  double r_edge_min = -1;
  double r_edge_max = 0;

  apf::Vector3 pos_in(0,0,0);
  apf::Vector3 pos_out(0,0,0);

  for (int i = 0; i < ee.size(); i++) {
    m->getDownward(ee.at(i), 0, d_v);
    int i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[0]);
    if(i1 == -1) {
      i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[1]);
      assert(i1 > -1);

      m->getPoint(d_v[0], 0, pos_out);
      m->getPoint(d_v[1], 0, pos_in);
    }
    else {
      i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[1]);
      //assert(i1 == -1);
      assert(i1 == -1);

      m->getPoint(d_v[1], 0, pos_out);
      m->getPoint(d_v[0], 0, pos_in);
    }
    double meas = (pos_out-pos_in).getLength();
    if (r_edge_min < 0 or r_edge_min > meas) {
      r_edge_min = meas;
    }
    if (r_edge_max < meas) {
      r_edge_max = meas;
    }
  }
  //for (int i = 0; i < ent.at(1).size(); i++) {
  //  double meas = vd_meas_ent(m, ent.at(1).at(i));
  //  if (r_edge_min < 0 or r_edge_min > meas) {
  //    r_edge_min = meas;
  //  }
  //}

  double r_cvx = r_edge_min/cvx_rat;
  //double r_cvx = r_edge_max/cvx_rat2;

  apf::Vector3 pos_curr(0, 0, 0);

  // Check if preconditioning is needed:
  //chk_precond();
  if(precond_flag) {

    // New approach:
    // Check and fix convexity around intersection points. 
    // Edge end vertex position, and intersection point with the sphere.
    std::map<apf::MeshEntity*, apf::Vector3> int_pt{};
    std::map<apf::MeshEntity*, apf::Vector3> vert_pt{};

    // The list of new vertices and the direction that they can relax in.
    std::vector<apf::MeshEntity*> v_sp_set(0);
    std::map<apf::MeshEntity*, apf::Vector3> v_sp_pos{};
    std::map<apf::MeshEntity*, apf::Vector3> v_ori_dir{};
    std::map<apf::MeshEntity*, apf::Vector3> v_mer_pos{};

    //get_int_pos(ee, int_pt, vert_pt);
    //split_ctr(ee, es_in, et, v_sp_set, int_pt, v_sp_pos, v_mer_pos, v_ori_dir, r_cvx);
    get_int_pos(ee, v_mer_pos, int_pt, vert_pt, r_cvx);
    // TODO Convex hull sometimes fails with edges having undefined faces. The 
    // debugging is intricate so I'm skipping convex hull to see how relaxation 
    // alone fixes things.
    //split_cvx(ee, es_in, et, v_sp_set, int_pt, v_sp_pos, v_mer_pos, v_ori_dir);
    split_noncvx(ee, es_in, et, v_sp_set, int_pt, v_sp_pos, v_mer_pos, v_ori_dir, r_cvx);

    vd_set_up(m, &ent.at(0), &ee);
    vd_set_up(m, &ee, &es_in);
    vd_set_up(m, &es_in, &et);
    vd_set_down(m, &et, &es);

    vd_remove_set(&es, &es_in);

    // Remove the collapsing edges.
    for (int i = 0; i < ent.at(1).size(); i++) {
      assert(e_map[ent.at(1).at(i)] == 1);

      std::vector<apf::MeshEntity*>::iterator it;
      it = std::find (ee.begin(), ee.end(), ent.at(1).at(i));
      assert(it != ee.end());
      ee.erase(it);
    }
    double t0 = PCU_Time();
    bool relaxed = false;
    //if(sp_trial_load) {
      //bool relaxed = relax_inv_grad(v_sp_set, r_edge_min);
      relaxed = relax_inv_conjgrad(v_sp_set, r_edge_min);
    //}
    //else {
    //  relaxed = set_v_sp_pos(ee);
    //}
    double t1 = PCU_Time();
    std::cout << "Relax took " << t1 - t0 << "s." << std::endl;
    //relax_inv_dir(et, es_in, ee, v_sp_set, v_sp_pos, v_ori_dir, v_mer_pos,
    //                                                      r_edge_min);
    return relaxed;
  }
  return true;
/*
calc the volumes of the non collapsing tets
  tag inverting
calc m
while any inverting tet
  calc virtual force
    find max virtual force
  find the motion along force to fix the minimum inversion
    for inverting volumes, del v_i / del m_j = v_curr + 0.1*v_start
*/
}

// Calculate the terms in the softmax function. B is the 1/6 of the sum of the 
// inward pointing normals of the triangles adjacent to vertex v.
apf::Vector3 vd_glens::calc_B(std::vector<apf::MeshEntity*> &v_set, 
        apf::MeshEntity* v, std::map<apf::MeshEntity*, apf::Vector3> &pos) {

  int v_sense[4] = {1, -1, 1, -1};
  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find (v_set.begin(), v_set.end(), v);
  int i1 = std::distance(v_set.begin(), it);

  int ids[3] = {1, 2, 3};
  apf::Vector3 B(0,0,0);
  for(int i = 0; i < 3; i++) {
    int i2 = (i1+ids[i]) %4;
    int i3 = (i1+ids[(i+1)%3]) %4;
    B = B + vd_cross(pos[v_set.at(i2)] - pos[v_set.at(i1)], 
                     pos[v_set.at(i3)] - pos[v_set.at(i1)]);
  }
  return B/6*v_sense[i1];
}

// The rate of change of volume wrt. alfa. If v is specified, only consider
// the residual at that vertex. iith vertex is fixed.
double vd_glens::calc_C(std::vector<apf::MeshEntity*> &v, 
                      std::map<apf::MeshEntity*, apf::Vector3> &pos, int ii, 
                                            apf::MeshEntity* v_r) {

  int v_sense[4] = {1, -1, 1, -1};

  std::vector<apf::Vector3> r(4, apf::Vector3(0,0,0));
  double C = 0;

  assert(f_field);
  for(int i = 0; i < 4; i++) {
    if(i != ii) {
      apf::getVector(f_field, v.at(i), 0, r.at(i));
    }
  }

  int ids[3] = {1, 2, 3};
  if(v_r != NULL) {
    std::vector<apf::MeshEntity*>::iterator it;
    it = std::find (v.begin(), v.end(), v_r);
    assert(it != v.end());
    int i1 = std::distance(v.begin(), it);
    // Shift in ids array:
    int id_shift = (i1 + 3 - ii) %4;
    int i2 = (ii+ids[(id_shift+1)%3])%4;
    int i3 = (ii+ids[(id_shift+2)%3])%4;
    C = C + vd_trip_prod(pos[v.at(i3)] - pos[v.at(ii)], 
                     r.at(i1), pos[v.at(i2)] - pos[v.at(ii)]);
  }
  else {
    for(int i = 0; i < 3; i++) {
      int i1 = (ii+ids[(i)%3]) %4;
      int i2 = (ii+ids[(i+1)%3]) %4;
      int i3 = (ii+ids[(i+2)%3]) %4;
      //C = C + vd_trip_prod(pos[v.at(i3)] - pos[v.at(i0)], 
      //                 pos[v.at(i1)] - pos[v.at(i0)], 
      //                 pos[v.at(i2)] - pos[v.at(i0)]);
      C = C + vd_trip_prod(pos[v.at(i3)] - pos[v.at(ii)], 
                       r.at(i1), pos[v.at(i2)] - pos[v.at(ii)]);

    }
  }
  return C/6*v_sense[ii];
}

bool vd_glens::relax_inv_grad(std::vector<apf::MeshEntity*> &v_sp_set,
                              double r_edge_min) {
  int iter_limit = 400;
  apf::Downward d_v;
  apf::Downward d_e;
  std::vector<apf::MeshEntity*>::iterator it;

  std::vector<apf::MeshEntity*> ee(0);
  std::vector<apf::MeshEntity*> es(0);
  std::vector<apf::MeshEntity*> et(0);

  vd_set_up(m, &v_sp_set, &ee);
  vd_set_up(m, &ee, &es);
  vd_set_up(m, &es, &et);

  std::map<apf::MeshEntity*, int> t_nbr{};
  for(int i = 0; i < v_sp_set.size(); i++) {
    t_nbr[v_sp_set.at(i)] = 0;
  }

  std::map<apf::MeshEntity*, apf::Vector3> pos_ori{};
  std::map<apf::MeshEntity*, apf::Vector3> v_pos{};
  std::map<apf::MeshEntity*, apf::Vector3> v_disp{};
  apf::Vector3 temp(0,0,0);

  std::vector<std::vector<apf::MeshEntity*> > tet_v (et.size(),
                                       std::vector<apf::MeshEntity*>(4));
  std::vector<std::vector<int> > tet_v_id (et.size(), std::vector<int>(4));

  for(int i = 0; i < et.size(); i++) {
    m->getDownward(et.at(i), 0, d_v);
    for(int j = 0; j < 4; j++) {
      tet_v.at(i).at(j) = d_v[j];
      m->getPoint(d_v[j], 0, temp);
      v_pos[d_v[j]] = temp;

      it = std::find (v_sp_set.begin(), v_sp_set.end(), tet_v.at(i).at(j));
      if(it != v_sp_set.end()) {
        int i1 = std::distance(v_sp_set.begin(), it);
        tet_v_id.at(i).at(j) = i1;
        t_nbr[v_sp_set.at(i1)] = t_nbr[v_sp_set.at(i1)] + 1;
        v_disp[v_sp_set.at(i1)] = apf::Vector3(0,0,0);
      }
      else
        tet_v_id.at(i).at(j) = -1;
    }
  }

  // Adjacent tets of split vertices on the split surface.
  std::vector<std::vector<int> > v_tet (v_sp_set.size(), std::vector<int>(0));
  // Adjacent split vertices of split vertices on the split surface.
  std::vector<std::vector<int> > v_vert (v_sp_set.size(), std::vector<int>(0));
  for(int i = 0; i < v_sp_set.size(); i++) {
    v_tet.at(i).reserve(t_nbr[v_sp_set.at(i)]);
    // In case exterior, there could be one more vertex:
    v_vert.at(i).reserve(t_nbr[v_sp_set.at(i)] + 1);

    m->getPoint(v_sp_set.at(i), 0, temp);
    pos_ori[v_sp_set.at(i)] = temp;
  }

  shell sh_curr;
  std::map<apf::MeshEntity*, bool> v_skip{};

  if(upd_shell and f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
    if(ext_cell) {
      for(int i = 0; i < v_sp_set.size(); i++) {
        apf::ModelEntity* mdl = m->toModel(v_sp_set.at(i));
        int c_dim = m->getModelType(mdl);
        int c_tag = m->getModelTag(mdl);
        if(c_dim < 3 and c_base->get_cell_ext_gmi(c_dim, c_tag)) {
          v_skip[v_sp_set.at(i)] = true;
        }
      }
    }
  }


  for(int i = 0; i < et.size(); i++) {
    for(int j = 0; j < 4; j++) {
      if(tet_v_id.at(i).at(j) > -1)
        v_tet.at(tet_v_id.at(i).at(j)).push_back(i);
    }
    for(int j = 0; j < 4; j++) {
      int i1 = tet_v_id.at(i).at(j);
      if(i1 > -1) {
        for(int k = j+1; k < 4; k++) {
          int i2 = tet_v_id.at(i).at(k);
          if(i2 > -1) {
            v_vert.at(i1).push_back(i2);
            v_vert.at(i2).push_back(i1);
          }
        }
      }
    }
  }

  for(int i = 0; i < v_sp_set.size(); i++) {
    std::sort(v_vert.at(i).begin(), v_vert.at(i).end());
    std::vector<int>::iterator it;
    it = std::unique (v_vert.at(i).begin(), v_vert.at(i).end());
    v_vert.at(i).resize(std::distance(v_vert.at(i).begin(),it));
  }

  std::map<apf::MeshEntity*, bool> tet_surr{};
  std::map<apf::MeshEntity*, bool> tet_skip{};
  std::map<apf::MeshEntity*, bool> tet_inv{};
  std::map<apf::MeshEntity*, int> tet_i1{};

  std::vector<apf::Vector3> pts (4, apf::Vector3(0,0,0));
  // Starting volume
  std::vector<double> vol_t (et.size());

  // Tag inverting
  apf::Numbering* tagnumbering;
  apf::Numbering* tetset;
  //apf::Numbering* a_f;
  if (m->findNumbering("inv_tag")) {
    //printf("grain_memb tag already exists.\n");
    tagnumbering = m->findNumbering("inv_tag");
    apf::destroyNumbering(tagnumbering);
    tagnumbering = apf::createNumbering(m, "inv_tag", apf::getConstant(3),1);
  }
  else {
    tagnumbering = apf::createNumbering(m, "inv_tag", apf::getConstant(3),1);
  }
  //apf::Numbering* a_f;
  if (m->findNumbering("tet_tag")) {
    //printf("grain_memb tag already exists.\n");
    tetset = m->findNumbering("tet_tag");
    apf::destroyNumbering(tetset);
    tetset = apf::createNumbering(m, "tet_tag", apf::getConstant(3),1);
  }
  else {
    tetset = apf::createNumbering(m, "tet_tag", apf::getConstant(3),1);
  }
  //if (m->findNumbering("a_f")) {
    //printf("grain_memb tag already exists.\n");
  //  a_f = m->findNumbering("a_f");
  //  apf::destroyNumbering(a_f);
  //  a_f = apf::createNumbering(m, "a_f", apf::getConstant(3),1);
  //}
  //else {
  //  a_f = apf::createNumbering(m, "a_f", apf::getConstant(3),1);
  //}

  apf::MeshEntity* e;

  apf::MeshIterator* e_it = m->begin(3);
  while ((e = m->iterate(e_it))) 
  {
    apf::number(tagnumbering, e, 0, 0, 2);
    apf::number(tetset, e, 0, 0, 0);
    //apf::number(a_f, e, 0, 0, 0);
  }
  m->end(e_it);

  double vol_scale = 0;
  double vol_max = 0;
  double vol_min = -1;
  // Most negative volume, used to scale w, st. w*vol_t/vol_scale = 100;
  double vol_neg = 0;
  // The volume threshold aimed as percentage of the original tet volume.
  double v_th_tol = 0.005;
  double v_th = 0.01;
  double w = 100;
  double w_m = 20;

  // The relative amount of change in the parameter m allowed in one step, in 
  // relation to the current one.
  double tol_step = 0.1;
  double tol_inv = 10;

  // The tolerances for convergence of the step multiplier.
  double tol_alfa = 0.01;

  for (int i = 0; i < et.size(); i++) {
    if(e_map[et.at(i)] != 1) {
      // Store the vertices, such that 0th vertex is the merging vertex and
      // the others are on the exterior surface.
      m->getDownward(et.at(i), 0, d_v);
      for(int j = 0; j < 4; j++) {
        m->getPoint(d_v[j], 0, pts.at(j));
      }
      vol_t.at(i) = vd_volume_tet(&pts);

      if(std::fabs(vol_t.at(i)) < vol_min or
         vol_min < - std::numeric_limits<double>::min())
        vol_min = std::fabs(vol_t.at(i));
    }
  }
  // Try total of volumes:
  vol_scale = 0;

  for (int i = 0; i < et.size(); i++) {
    apf::number(tetset, et.at(i), 0, 0, 1);

    if(e_map[et.at(i)] != 1) {
      int ent_type = m->getType(et.at(i));
      int d = m->typeDimension[ent_type];
      assert(d == 3); 

      m->getDownward(et.at(i), 0, d_v);
      m->getDownward(et.at(i), 1, d_e);

      std::vector<apf::MeshEntity*>::iterator it;
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;

      int i1 = -1;
      int i2 = -1;
      int i4 = -1;
      //int i3 = -1;
      for (int j = 0; j < 4; j++) {
        it = std::find (ent.at(0).begin(), ent.at(0).end(), d_v[j]);
        //int i_curr = std::distance(ent.at(0).begin(), it);
        if(it != ent.at(0).end()) {
          //std::cout << "d_v[" << j << "] = " << d_v[j] << std::endl;

          if(!found2) {
            if(!found1) {
              found1 = true;
              i1 = j;
              //std::cout << "\tfirst" << std::endl;
            }
            else {
              found2 = true;
              i2 = j;
              //std::cout << "\tsecond" << std::endl;
            }
          }
          else {
            found3 = true;
            //std::cout << "\tthird" << std::endl;
          }
        }
        else {
          i4 = j;
        }
      }
      if(found2) {
        apf::number(tagnumbering, et.at(i), 0, 0, 1);

        tet_skip[et.at(i)] = true;
      }
      else if(found1) {
        tet_i1[et.at(i)] = i1;
        // Store the vertices, such that 0th vertex is the merging vertex and
        // the others are on the exterior surface.
        for(int j = 0; j < 4; j++) {
          tet_v.at(i).at(j) = d_v[j];
          m->getPoint(d_v[j], 0, pts.at(j));
        }
        //pts.at(0) = midpoint;
        pts.at(i1) = midpoint;
        vol_t.at(i) = vd_volume_tet(&pts);
        if(vol_t.at(i) < std::numeric_limits<double>::min()) {
          tet_inv[et.at(i)] = true;
          //std::cout << et.at(i) << " inverting."
          //          << " vol_t: " << vol_t.at(i)
          //          << std::endl;
          apf::number(tagnumbering, et.at(i), 0, 0, 3);
          if(vol_t.at(i) < vol_neg)
            vol_neg = vol_t.at(i);
        }

        if(std::fabs(vol_t.at(i)) > vol_max)
          vol_max = vol_t.at(i);
      }
      else {
        tet_surr[et.at(i)] = true;
        tet_i1[et.at(i)] = i4;
        for(int j = 0; j < 4; j++) {
          tet_v.at(i).at(j) = d_v[j];
          m->getPoint(d_v[j], 0, pts.at(j));
        }
        vol_t.at(i) = vd_volume_tet(&pts);
        v_pos[d_v[i4]] = pts.at(i4);
        if(std::fabs(vol_t.at(i)) > vol_max)
          vol_max = vol_t.at(i);
      }
    }
    else {
      vol_t.at(i) = vd_volume_tet(m, et.at(i));
      apf::number(tagnumbering, et.at(i), 0, 0, 1);
    }

    vol_scale = vol_scale + vol_t.at(i);
  }
  //vol_scale = getMedianEntSize_hist(vol_t);
  //if(vol_scale < std::numeric_limits<double>::min())
  //  vol_scale = vol_max;


  for (int i = 0; i < et.size(); i++) {
    if(e_map[et.at(i)] != 1 and !tet_skip[et.at(i)] and !tet_inv[et.at(i)]) {
      if(vol_t.at(i) < v_th*vol_scale) {
        std::cout << et.at(i) << " the volume " << vol_t.at(i) 
                  << " is smaller than" << v_th*vol_scale
                  << std::endl;
        apf::number(tagnumbering, et.at(i), 0, 0, 3);
      }
    }
  }
  for (int i = 0; i < ent.at(3).size(); i++) {
    apf::number(tagnumbering, ent.at(3).at(i), 0, 0, 1);
    apf::number(tetset, ent.at(3).at(i), 0, 0, 2);
  }

  apf::writeVtkFiles("./output/Glens_inv", m);

  for (int i = 0; i < ent.at(0).size(); i++) {
    v_pos[ent.at(0).at(i)] = midpoint;
  }
  if(vol_neg > - std::numeric_limits<double>::min())
    vol_neg = -vol_scale;

  std::vector<double> A(et.size(), 0);
  std::vector<double> C(et.size(), 0);
  std::vector<std::vector<apf::Vector3> > B(v_sp_set.size(), 
                   std::vector<apf::Vector3> (0, apf::Vector3(0,0,0)) );
  std::vector<double> exp_pw_tot(v_sp_set.size(), 0);

  for (int i = 0; i < v_sp_set.size(); i++) {
    B.at(i).resize(v_tet.at(i).size());
  }

  // Scale it such that most negative volume tetrahedron has A about e^(2)
  // w * v_th = -2*v_th_tol*vol_min/vol_neg
  // w *- vol_t.at(i)/vol_scale = 10/vol_neg * vol_t.at(i)
  w = -80*vol_scale/vol_neg;
  v_th = std::fabs(v_th_tol*vol_min/vol_scale);
  double v_th_i = v_th;

  bool inv_flag = true;
  int iter = 0;
  while(iter < iter_limit and inv_flag) {
    iter = iter + 1;

    w = -80*vol_scale/vol_neg;
    if(w*v_th_i > 1)
      v_th = std::fabs(v_th_tol*vol_min/vol_scale);
    else
      v_th = v_th_i;

    inv_flag = false;
    //v_th = v_th_tol*vol_min/vol_scale;

    apf::Vector3 force(0,0,0);
    for (int i = 0; i < et.size(); i++) {
      if(e_map[et.at(i)] != 1) {
        A.at(i) = v_th - vol_t.at(i)/vol_scale;
        //if(A.at(i) > 10)
        //  A.at(i) = 10; 
        A.at(i) = std::exp(w*A.at(i));
        //apf::number(a_f, et.at(i), 0, 0, A.at(i));
      }
    }

    for (int i = 0; i < v_sp_set.size(); i++) {
      force = apf::Vector3(0,0,0);
      exp_pw_tot.at(i) = 0;

      bool skipped = true;
      //std::cout << v_sp_set.at(i) << std::endl;
      for (int j = 0; j < v_tet.at(i).size(); j++) {
        apf::MeshEntity* tet_curr = et.at(v_tet.at(i).at(j));
        if(e_map[tet_curr] != 1 and !tet_skip[tet_curr]) {
          skipped = false;
          int t1 = v_tet.at(i).at(j);
          B.at(i).at(j) = calc_B(tet_v.at(t1), v_sp_set.at(i), v_pos)/vol_scale;
          B.at(i).at(j) = norm_0(B.at(i).at(j));
          //exp_pw.at(j) = w*(v_th - B.at(j)*v_disp[v_sp_set.at(i)]/vol_scale);
          //exp_pw.at(j) = std::exp(w*A.at(v_tet.at(i).at(j)));
          exp_pw_tot.at(i) = exp_pw_tot.at(i) + A.at(t1);
          //vd_print_vert(m, et.at(t1));
          //std::cout << "\tA(" << j+1 << ") =" << A.at(t1) << ";" 
          //          << std::endl;
          //std::cout << "\tB(" << j+1 << ") = [" 
          //          << B.at(i).at(j)[0] << "," << B.at(i).at(j)[1] << "," 
          //          << B.at(i).at(j)[2] << ";"
          //          << std::endl;
        }
        else {
          B.at(i).at(j) = apf::Vector3(0,0,0);
        }
      }
      exp_pw_tot.at(i) = exp_pw_tot.at(i) + 1;
      for (int j = 0; j < v_tet.at(i).size(); j++) {
        int t1 = v_tet.at(i).at(j);
        apf::MeshEntity* tet_curr = et.at(t1);
        if(e_map[tet_curr] != 1 and !tet_skip[tet_curr]) {
          force = force - B.at(i).at(j)*A.at(t1);
        }
      }
      if(!skipped)
        force = force/exp_pw_tot.at(i);
      apf::setVector(f_field, v_sp_set.at(i), 0, force);

      //std::cout << " f_i: " << force << std::endl;
    }
 
    // Transfer the force due to potentials on neighboring vertices.
/*
    for (int i = 0; i < v_sp_set.size(); i++) {
      for (int j = 0; j < v_vert.at(i).size(); j++) {
        int i2 = v_vert.at(i).at(j);
        for(int k = 0; k < v_tet.at(i).size(); k++) {
          int t1 = v_tet.at(i).at(k);
          apf::MeshEntity* tet_curr = et.at(t1);
          if(e_map[tet_curr] != 1 and !tet_skip[tet_curr]) {
            for(int l = 0; l < v_tet.at(i2).size(); l++) {
              int t2 = v_tet.at(i2).at(l);
              if(t1 == t2) {
                apf::getVector(f_field, v_sp_set.at(i), 0, force);
                //std::cout << v_sp_set.at(i) << "-" << v_sp_set.at(i2)
                //          << "(" << tet_curr << ")" << "\n"
                //          << " f_i: " << force << std::endl;
                force = force - B.at(i).at(k)*A.at(t1)/exp_pw_tot.at(i2);
                apf::setVector(f_field, v_sp_set.at(i), 0, force);
                //std::cout << " f_f: " << force
                //          << std::endl;
              }
            }
          }
        }
      }
    }
*/
    // Calculate alfa
    for (int i = 0; i < et.size(); i++) {
      if(e_map[et.at(i)] != 1 and !tet_skip[et.at(i)]) {
        C.at(i) = calc_C(tet_v.at(i), v_pos, tet_i1[et.at(i)])/vol_scale;
      }
      else
        C.at(i) = 0;
    }

    double alfa = -1;
    for (int i = 0; i < v_sp_set.size(); i++) {

      double dPhida = 0;
      double Phi = 0;
      for (int j = 0; j < v_tet.at(i).size(); j++) {
        int i_tet = v_tet.at(i).at(j);
        apf::MeshEntity* tet_curr = et.at(i_tet);
        if(e_map[tet_curr] != 1 and !tet_skip[tet_curr]) {
          double C_curr = calc_C(tet_v.at(i_tet), v_pos, tet_i1[et.at(i_tet)], 
                                              v_sp_set.at(i))/vol_scale;
          //dPhida = dPhida + A.at(i_tet)*C.at(i_tet);
          dPhida = dPhida - A.at(i_tet)*C_curr;
        }
      }
      dPhida = dPhida/exp_pw_tot.at(i);
      Phi = std::log(exp_pw_tot.at(i))/w;
/*
      for (int j = 0; j < v_vert.at(i).size(); j++) {
        int i2 = v_vert.at(i).at(j);
        for(int k = 0; k < v_tet.at(i).size(); k++) {
          int t1 = v_tet.at(i).at(k);
          apf::MeshEntity* tet_curr = et.at(t1);
          if(e_map[tet_curr] != 1 and !tet_skip[tet_curr]) {
            for(int l = 0; l < v_tet.at(i2).size(); l++) {
              int t2 = v_tet.at(i2).at(l);
              if(t1 == t2) {
                std::cout << v_sp_set.at(i) << "-" << v_sp_set.at(i2)
                          << "(" << tet_curr << ")" << "\n"
                          << " dPhida_i: " << dPhida << std::endl;
                dPhida = dPhida - A.at(t1)*C.at(t1)/exp_pw_tot.at(i2);
                std::cout << " dPhida_f: " << dPhida << std::endl;
              }
            }
          }
        }
      }
*/

      //std::cout << v_sp_set.at(i) 
      //          << " Phi: " << Phi
      //          << " dPhida: " << dPhida
      //          << std::endl;
    }

    
    for (int i = 0; i < et.size(); i++) {
      // The multiplier to reach half of current volume: 
      // For small tets with faster moving vertices, this will be smaller.
      // For larger tets this will be larger. The smallest is picked.

      if(e_map[et.at(i)] != 1 and !tet_skip[et.at(i)]) {
        double alfa_curr = 0;
        if(vol_t.at(i)/vol_scale > v_th) {
          for(int j = 0; j < 4; j++) {
            if(j != tet_i1[et.at(i)]) {
              double C_curr = calc_C(tet_v.at(i), v_pos, tet_i1[et.at(i)], 
                                              tet_v.at(i).at(j));
              //                                tet_v.at(i).at(j))/vol_scale;
              if(alfa_curr < std::numeric_limits<double>::min()
                  and C_curr < -std::numeric_limits<double>::min()
                  and std::fabs(vol_t.at(i) > vol_min*1.1)) {
                  //alfa_curr = -vol_t.at(i)/C_curr/vol_scale;
                  alfa_curr = -vol_t.at(i)/C_curr;

                if(std::fabs(alfa_curr) > std::numeric_limits<double>::min()) {
                  if(alfa_curr < alfa or 
                     alfa < - std::numeric_limits<double>::min()) {
                    alfa = alfa_curr;
                    //std::cout << et.at(i) 
                    //          << " V: " << vol_t.at(i)/vol_scale
                    //          << " C: " << C_curr
                    //          << std::endl;
                  }
                }
              }
            }
          }
        }
        //if(vol_t.at(i)/vol_scale > v_th) {
        //  if(C.at(i) < -std::numeric_limits<double>::min())
        //    alfa_curr = -vol_t.at(i)/C.at(i)/vol_scale;
         //else if (C.at(i) > std::numeric_limits<double>::min())
         //   alfa_curr = vol_t.at(i)/C.at(i)/vol_scale;
        //}
        //else if(vol_t.at(i) < - std::numeric_limits<double>::min()) {
        //  if(C.at(i) < -std::numeric_limits<double>::min())
        //    alfa_curr = vol_t.at(i)/C.at(i)/vol_scale/10;
        //}
        //if(std::fabs(alfa_curr) > std::numeric_limits<double>::min()) {
        //  if(alfa_curr < alfa or alfa < - std::numeric_limits<double>::min()) {
        //    alfa = alfa_curr;
        //    std::cout << et.at(i) 
        //              << " V: " << vol_t.at(i)/vol_scale
        //              << " C: " << C.at(i)
        //              << std::endl;
        //  }
        //}
      }
    }

    for (int i = 0; i < v_sp_set.size(); i++) {
      force = apf::Vector3(0,0,0);
      apf::getVector(f_field, v_sp_set.at(i), 0, force);
      if(!v_skip[v_sp_set.at(i)]) {
        //v_disp[v_sp_set.at(i)] = force*alfa/2 + v_disp[v_sp_set.at(i)]/2;
        v_disp[v_sp_set.at(i)] = force*alfa/4;
        //std::cout << v_sp_set.at(i)
        //          << " pos: " << v_pos[v_sp_set.at(i)]
        //          << " disp: " << v_disp[v_sp_set.at(i)]
        //          << std::endl;
      }
      else {
        if (f_calc->chk_vert_special(m, v_sp_set.at(i)))
          force = f_calc->get_vec_special(m, v_sp_set.at(i), force);
        v_disp[v_sp_set.at(i)] = force*alfa/4;
      }
      v_pos[v_sp_set.at(i)] = v_pos[v_sp_set.at(i)] + v_disp[v_sp_set.at(i)];
      m->setPoint(v_sp_set.at(i), 0, v_pos[v_sp_set.at(i)]);

    }

/*

    // Calculate alfa
    for (int i = 0; i < v_sp_set.size(); i++) {
      //Calculate alfa_curr
      for (int j = 0; j < v_tet.at(i).size(); j++) {
        apf::MeshEntity* tet_curr = et.at(v_tet.at(i).at(j));
        if(e_map[tet_curr] != 1 and !tet_skip[tet_curr]) {
          a.at(j) = B.at(j)*force;
          b.at(j) = std::exp(w*A.at(j));
        }
      }
      // Start with alfa = 0, find the zero by Newton-Rhapson of f:
      // f = sum(a*b**exp(alfa)^a)
      // alfa = ln(gam)
      double alfa_curr = 0;
      double gam_curr = 1;
      double gam_next = 2;
      if(!skipped) {
        while(std::fabs(std::log1p(gam_next) - std::log1p(gam_curr)/std::log1p(gam_curr)) > tol_alfa) {
          double f = 0;
          double f_prime = 0;
          for (int j = 0; j < v_tet.at(i).size(); j++) {
            apf::MeshEntity* tet_curr = et.at(v_tet.at(i).at(j));
            if(e_map[tet_curr] != 1 and !tet_skip[tet_curr]) {
              f = f + a.at(j)*b.at(j)*std::pow(gam_curr, a.at(j));
              f_prime = f_prime + a.at(j)*a.at(j)*b.at(j)
                              *std::pow(gam_curr, a.at(j)-1);
            }
          }
          assert(f_prime > std::numeric_limits<double>::min());
          gam_next = gam_curr - f/f_prime;
        }
        assert(gam_next > 1 - std::numeric_limits<double>::epsilon());
        alfa_curr = std::log1p(gam_next);
        assert(alfa_curr > - std::numeric_limits<double>::min());
      }
      if(alfa_curr < alfa or alfa < - std::numeric_limits<double>::min()) {
        alfa = alfa_curr;
      }
    }
    std::cout << "alfa: " << alfa << std::endl;
 
//gradient descent
    double force_max = 0;
    double force_min = -1;
    for (int i = 0; i < v_sp_set.size(); i++) {
      apf::getVector(f_field, v_sp_set.at(i), 0, force);
      double force_mag = force.getLength();
      if(force_mag > force_max)
        force_max = force_mag;
      if(force_min < -std::numeric_limits<double>::min() or 
         force_mag < force_min)
        force_min = force_mag;
    }
    for (int i = 0; i < v_sp_set.size(); i++) {
      double mot_mult_curr = -1;
      apf::getVector(f_field, v_sp_set.at(i), 0, force);
      force = force/force_max;
      apf::setVector(f_field, v_sp_set.at(i), 0, force);
      // A motion proportional to force, in the direction of the displacement
      // inverts a tetrahedron by a motion of
      // tol_inv*mot_mult_curr*force*v_disp
      std::cout << "mult for " << v_sp_set.at(i) << std::endl;
      mot_mult_curr = vd_dist_v_x(m, v_sp_set.at(i), force, tet_inv);
      if(mot_mult_curr < std::numeric_limits<double>::min()) {
        mot_mult_curr = mot_mult;
      }
      else if(mot_mult_curr < m_th) {
        force = apf::Vector3(0,0,0);
        apf::setVector(f_field, v_sp_set.at(i), 0, force);
      }
      else {
        if(force.getLength() < std::numeric_limits<double>::min())
          mot_mult_curr = 1;
        else
          mot_mult_curr = mot_mult_curr/force.getLength()/tol_inv;
      }
      if(mot_mult_curr > std::numeric_limits<double>::min() and
         (mot_mult < - std::numeric_limits<double>::min() 
                      or mot_mult_curr < mot_mult)) {
        mot_mult = mot_mult_curr;
      }
      std::cout << v_sp_set.at(i) 
                << " pos " << v_pos[v_sp_set.at(i)]
                << " mot_mult " << mot_mult_curr
                << " f_rel: " << force 
                << std::endl;
    }

    for (int i = 0; i < v_sp_set.size(); i++) {
      apf::getVector(f_field, v_sp_set.at(i), 0, force);
      v_disp[v_sp_set.at(i)] = force*mot_mult/2 + v_disp[v_sp_set.at(i)]/2;
      std::cout << v_sp_set.at(i)
                << " pos: " << v_pos[v_sp_set.at(i)]
                << " disp: " << v_disp[v_sp_set.at(i)]
                << std::endl;
      v_pos[v_sp_set.at(i)] = v_pos[v_sp_set.at(i)] + v_disp[v_sp_set.at(i)];
      m->setPoint(v_sp_set.at(i), 0, v_pos[v_sp_set.at(i)]);
    }
*/

    //vol_scale = 0;
    inv_flag = update_vol_inv_grad(et, vol_t, tet_inv, tet_surr, tet_skip, 
                                              tet_i1, tet_v,
                                                  vol_neg, vol_min,
                                                  v_th*vol_scale);

    apf::writeVtkFiles("./output/Glens_step", m);

    if(!inv_flag) {
      std::cout<< "No inversions remain!" << std::endl;
    }
    else
      std::cout<< "Threshold volume:" << v_th*vol_scale << std::endl;

    if(vol_neg > - std::numeric_limits<double>::min())
      vol_neg = -vol_scale;
  }
  if(inv_flag) {
    std::cout << "Relaxation was not successful!" << std::endl;
    for (int i = 0; i < v_sp_set.size(); i++) {
      m->setPoint(v_sp_set.at(i), 0, pos_ori[v_sp_set.at(i)]);
    }
  }

  return (!inv_flag);
}

// Calculate the potential.
double vd_glens::calc_energy(std::vector<apf::MeshEntity*> &v_sp_set,
                        std::vector<apf::MeshEntity*> &et,
                        std::vector<double> &C,
                  std::vector<std::vector<apf::MeshEntity*> > &tet_v,
                  std::vector<double> &exp_pw_tot,
                  std::map<apf::MeshEntity*, apf::Vector3> &v_pos,
                  std::map<apf::MeshEntity*, bool> &tet_skip,
                  std::map<apf::MeshEntity*, bool> &v_skip,
                  std::map<apf::MeshEntity*, int> &tet_i1,
                  double vol_scale, double w) {

  for (int i = 0; i < et.size(); i++) {
    if(e_map[et.at(i)] != 1 and !tet_skip[et.at(i)]) {
      C.at(i) = calc_C(tet_v.at(i), v_pos, tet_i1[et.at(i)])/vol_scale;
    }
    else
      C.at(i) = 0;
  }

  double Phi = 0;
  for (int i = 0; i < v_sp_set.size(); i++) {
    //double dPhida = 0;
    //for (int j = 0; j < v_tet.at(i).size(); j++) {
    //  int i_tet = v_tet.at(i).at(j);
    //  apf::MeshEntity* tet_curr = et.at(i_tet);
    //  if(e_map[tet_curr] != 1 and !tet_skip[tet_curr]) {
    //    double C_curr = calc_C(tet_v.at(i_tet), v_pos, tet_i1[et.at(i_tet)], 
    //                                        v_sp_set.at(i))/vol_scale;
    //    //dPhida = dPhida + A.at(i_tet)*C.at(i_tet);
    //    dPhida = dPhida - A.at(i_tet)*C_curr;
    //  }
    //}
    //dPhida = dPhida/exp_pw_tot.at(i);
    double Phi_local = std::log(exp_pw_tot.at(i))/w;
    Phi = Phi + Phi_local;
  }
  return Phi;
}

////////////////////////////////////////////////
bool vd_glens::chk_aperture(std::vector<apf::MeshEntity*> &v_sp_set) {

  std::vector<apf::MeshEntity*> ee(0);
  std::vector<apf::MeshEntity*> es(0);

  vd_set_up(m, &v_sp_set, &ee);
  vd_set_up(m, &ee, &es);

  std::vector<apf::MeshEntity*> t_3c(0);
  t_3c.reserve(es.size());

  apf::Downward d_v;
  // Check that all bounding vertices are split vertices.
  for(int i = 0; i < es.size(); i++) {
    m->getDownward(es.at(i), 0, d_v);
    bool found = true;
    for(int j = 0; j < 3; j++) {
      int i1 = findIn(&v_sp_set, v_sp_set.size(), d_v[j]);
      if(i1 == -1)
        found = false;
    }
    if(found) {
      assert(m->getModelType(m->toModel(es.at(i))) == 3);
      t_3c.push_back(es.at(i));
    }
  }

  apf::Downward d_e;
  std::map<apf::MeshEntity*, int> e_t_nbr{};
  std::map<apf::MeshEntity*, apf::MeshEntity*> e_t_map1{};
  std::map<apf::MeshEntity*, apf::MeshEntity*> e_t_map2{};

  bool disc_like = true;
  for(int i = 0; i < t_3c.size(); i++) {
    m->getDownward(t_3c.at(i), 1, d_e);
    for(int j = 0; j < 3; j++) {
      e_t_nbr[d_e[j]] = e_t_nbr[d_e[j]] + 1;
      if(!e_t_map1[d_e[j]])
        e_t_map1[d_e[j]] = t_3c.at(i);
      else if(!e_t_map2[d_e[j]])
        e_t_map2[d_e[j]] = t_3c.at(i);
      else {
        disc_like = false;
        i = t_3c.size();
      }
      if(e_t_nbr[d_e[j]] > 2) {
        disc_like = false;
        i = t_3c.size();
      }
    }
  }

  assert(disc_like);

  // Store the current ordering of vertices. 
  // If an orientation of any vertex pair is used, the ordering of the vertices
  // should be flipped. For a surface composed of edge sharing triangles, 
  // this should yield a consistent arrangement. 
  std::map<std::pair<int,int>, bool> e_used{};
  std::vector<int> v_id(3);
  std::map<apf::MeshEntity*, bool> t_burned{};
  std::map<apf::MeshEntity*, int> t_3c_map{};
  std::map<int, int> c3_count{};
  std::vector<std::vector<apf::MeshEntity*> > burning(2, 
                                      std::vector<apf::MeshEntity*> (0));
  burning.at(0).reserve(t_3c.size());
  burning.at(1).reserve(t_3c.size());
  int id_curr = 0;
  int id_next = 1;
  // The number of disjoint sets.
  int set_nbr = 0;
  for(int i = 0; i < t_3c.size(); i++) {
    if(!t_burned[t_3c.at(i)]) {
      set_nbr = set_nbr + 1;
      burning.at(id_curr).clear();
      burning.at(id_curr).push_back(t_3c.at(i));
      t_burned[burning.at(id_curr).back()] = true;
      t_3c_map[burning.at(id_curr).back()] = set_nbr;
      c3_count[set_nbr] = 1;
      apf::ModelEntity* mdl_curr = m->toModel(t_3c.at(i));
      while(!burning.at(id_curr).empty()) {
        apf::MeshEntity* t_curr = burning.at(id_curr).back();
        m->getDownward(t_curr, 1, d_e);

        // Collect the next triangles:
        for(int j = 0; j < 3; j++) {
          if(m->toModel(d_e[j]) == mdl_curr) {
            apf::MeshEntity* t_next = NULL;
            if(e_t_map1[d_e[j]] != t_curr) {
              assert(e_t_map2[d_e[j]] == t_curr);
              t_next = e_t_map1[d_e[j]];
            }
            else {
              assert(e_t_map2[d_e[j]] != NULL);
              t_next = e_t_map2[d_e[j]];
            }
            if(t_next != NULL and !t_burned[t_next]) {
              burning.at(id_next).push_back(t_next);
              t_burned[t_next] = true;
              t_3c_map[t_next] = set_nbr;
              c3_count[set_nbr] = c3_count[set_nbr] + 1;
            }
          }
        }
        burning.at(id_curr).pop_back();
        if(burning.at(id_curr).empty()) {
          id_next = id_curr;
          id_curr = (id_next + 1) % 2;
        }
      }
    }
  }

  std::vector<std::vector<apf::MeshEntity*> > c3_t(set_nbr, 
                                      std::vector<apf::MeshEntity*> (0));
  for(int i = 0; i < c3_t.size(); i++) {
    c3_t.at(i).reserve(c3_count[i+1]);
  }
  for(int i = 0; i < t_3c.size(); i++) {
    int set_curr = t_3c_map[t_3c.at(i)];
    c3_t.at(set_curr - 1).push_back(t_3c.at(i));
  }

  // Estimate the f-number (N) of the aperture area:
  // N = f/r, where r = sqrt(A/pi) and f is the average distance of the hull 
  // triangles.
  // Set a maximum N difference allowed to prevent stretched grains.

  double N_th = 2;
  double N_max = 0.;
  int i_max = -1;
  // Instead of 
  for(int i = 0; i < c3_t.size(); i++) {
    // Find the total area vector. Find the mean position of the associated 
    // merging vertices across. Find the change in aperture by their motion to 
    // the new midpoint.
    apf::Vector3 A_tot(0,0,0);
    apf::Vector3 v_pos_mean(0,0,0);
    apf::Vector3 t_pos_mean(0,0,0);
    for(int j = 0; j < c3_t.at(i).size(); j++) {
      apf::MeshEntity* tri_curr = c3_t.at(i).at(j);
      apf::Vector3 v_pos = get_pos_v_x_t_in_list(m, tri_curr, ent.at(0));
      apf::Vector3 t_pos = vd_get_pos(m, tri_curr);
      apf::Vector3 A_out = vd_area_out_n(m, tri_curr);
      if(A_out*(v_pos - t_pos) < -std::numeric_limits<double>::min())
        A_out = A_out*(-1);
      A_tot = A_tot + A_out;
      v_pos_mean = v_pos_mean + v_pos;
      t_pos_mean = t_pos_mean + t_pos;
    }
    v_pos_mean = v_pos_mean/c3_t.at(i).size();
    t_pos_mean = t_pos_mean/c3_t.at(i).size();

    double r_equi = std::sqrt(A_tot.getLength()/PI_L);
    A_tot = norm_0(A_tot);
    double f_old = A_tot*(v_pos_mean - t_pos_mean);
    double f_new = A_tot*(midpoint - t_pos_mean);

    double N_diff = std::fabs(f_old/r_equi - f_new/r_equi);
    if(std::isnan(N_diff))
      return false;
    if(N_diff > N_max) {
      N_max = N_diff;
      i_max = i;
    }
  }

/*
// TODO
//Is f-number a better measure or weighted aperture area (inverse of f-number)
//Minimum bound 4pi for the total of area.
  double N_th = 2;
  double N_max = 0.;
  int i_max = -1;
  for(int i = 0; i < c3_t.size(); i++) {
    double f_avg = 0;
    double A_tot = 0;


    for(int j = 0; j < c3_t.at(i).size(); j++) {
      double A_curr = vd_meas_ent(m, c3_t.at(i).at(j));
      double dist_curr = (vd_get_pos(m, c3_t.at(i).at(j)) - midpoint).getLength();
      A_tot = A_tot + A_curr;
      f_avg = f_avg + dist_curr*A_curr;
    }
    f_avg = f_avg/A_tot;
    double r_equi = std::sqrt(A_tot/PI_L);
    double N_curr = f_avg/r_equi;
    if(std::isnan(N_curr))
      return false;
    if(N_curr > N_max) {
      N_max = N_curr;
      i_max = i;
    }
  }
*/
  if(N_max > N_th)
    return false;
  else
    return true;
}


bool vd_glens::relax_inv_conjgrad(std::vector<apf::MeshEntity*> &v_sp_set,
                              double r_edge_min) {
  apf::Downward d_v;
  apf::Downward d_e;
  std::vector<apf::MeshEntity*>::iterator it;

  std::vector<apf::MeshEntity*> ee(0);
  std::vector<apf::MeshEntity*> es(0);
  std::vector<apf::MeshEntity*> et(0);

  vd_set_up(m, &v_sp_set, &ee);
  vd_set_up(m, &ee, &es);
  vd_set_up(m, &es, &et);

  vd_relax_cont rc(m, c_base, f_calc, cell_col_dim, cell_col_id, 
                             v_sp_set, ent, et, e_map,
                             ext_cell, upd_shell, midpoint, r_edge_min);

  bool inv_flag = rc.relax();
  bool aperture_OK = chk_aperture(v_sp_set);
  std::cout << "Aperture: " << aperture_OK << std::endl;
  return (!inv_flag and aperture_OK);
}


bool vd_glens::update_vol_inv_grad(std::vector<apf::MeshEntity*> &et, 
                                   std::vector<double> &vol_t,
                                   std::map<apf::MeshEntity*, bool> &tet_inv,
                                   std::map<apf::MeshEntity*, bool> &tet_surr,
                                   std::map<apf::MeshEntity*, bool> &tet_skip,
                                   std::map<apf::MeshEntity*, int> &tet_i1,
                           std::vector<std::vector<apf::MeshEntity*> > &tet_v,
                   double &vol_neg, double &vol_min, double min_vol) {
  std::vector<apf::Vector3> pts(4, apf::Vector3(0,0,0));

  apf::Numbering* tagnumbering = m->findNumbering("inv_tag");
  assert(tagnumbering);

  bool inv_flag = false;
  vol_neg = vol_min;

  for (int i = 0; i < et.size(); i++) {
    if(tet_surr[et.at(i)]) {
      for(int j = 0; j < 4; j++) {
        m->getPoint(tet_v.at(i).at(j), 0, pts.at(j));
      }
      vol_t.at(i) = vd_volume_tet(&pts);
      if(vol_t.at(i) < std::numeric_limits<double>::min() and
         vol_t.at(i) < vol_neg)
        vol_neg = vol_t.at(i);

      if(vol_t.at(i) < min_vol) {
        inv_flag = true;
        tet_inv[et.at(i)] = true;
        std::cout << "tet " << i << " inverting."
                  << " vol_t: " << vol_t.at(i)
                  << std::endl;
        apf::number(tagnumbering, et.at(i), 0, 0, 3);
      }
      else
        apf::number(tagnumbering, et.at(i), 0, 0, 2);
      if(std::fabs(vol_t.at(i)) < vol_min or
         vol_min < - std::numeric_limits<double>::min())
        vol_min = std::fabs(vol_t.at(i));
      //if(std::fabs(vol_t.at(i)) > vol_scale)
      //  vol_scale = vol_t.at(i);
    }

    else if(!tet_skip[et.at(i)] and e_map[et.at(i)] != 1) {
      tet_inv[et.at(i)] = false;
      for(int j = 0; j < 4; j++) {
        m->getPoint(tet_v.at(i).at(j), 0, pts.at(j));
      }
      pts.at(tet_i1[et.at(i)]) = midpoint;
      vol_t.at(i) = vd_volume_tet(&pts);
      //if(vol_t.at(i) < std::numeric_limits<double>::min()) {
      if(vol_t.at(i) < min_vol) {
        if(vol_t.at(i) < std::numeric_limits<double>::min() and
            vol_t.at(i) < vol_neg)
          vol_neg = vol_t.at(i);

        inv_flag = true;
        tet_inv[et.at(i)] = true;
        std::cout << "tet " << i << " inverting."
                  << " vol_t: " << vol_t.at(i)
                  << std::endl;
        apf::number(tagnumbering, et.at(i), 0, 0, 3);
      }
      else {
        apf::number(tagnumbering, et.at(i), 0, 0, 2);
      }
      //if(std::fabs(vol_t.at(i)) > vol_scale)
      //  vol_scale = vol_t.at(i);
      if(std::fabs(vol_t.at(i)) < vol_min or
         vol_min < - std::numeric_limits<double>::min())
        vol_min = std::fabs(vol_t.at(i));

    }
    else {
      apf::number(tagnumbering, et.at(i), 0, 0, 1);
    }
  }
  return inv_flag;
}

vd_glens::vd_glens(apf::Mesh2* msh, struct cell_base* c, vd_entlist* el_in) :
    m(NULL), mdl(NULL), c_base(NULL), 
    bipy_split(NULL), tet_split(NULL), lens_split(NULL), f_calc(NULL),
    vd_par(),
    sh_min(-1,-1), 
    vrfy_msh(true), inv_flag(true), precond_flag(true), shell_flag(true), 
    save_vtk(true), save_vtk_sub(false), vtk_name("./output/"), smb_name("./tempmesh/temp"), 
    mod_c_base(true),
    e_list(NULL),
    cvx_rat(2),
    cvx_rat2(10),
    cell_col_dim(0), cell_col_id(0), c_flag(true), tag_0c(-1),
    ext_cell(false), 
    ext_c0(0), 
    ext_c0_vec(0, apf::Vector3(0,0,0)),
    type_vertex(0), tag_vertex(0), 
    ch_f(NULL), sp_f(NULL), f_field(NULL),
    c0_flag(false), upd_shell(false),

    ent(0, std::vector<apf::MeshEntity*>(0)),
    m_list(0, std::vector<apf::ModelEntity*>(0)),
    v_list(0, std::vector<std::vector<apf::MeshEntity*> > (0, 
           std::vector<apf::MeshEntity*>(0) ) ),
    e_col_ext(0),
    pos_merg{}, 
    e_map(), v_map(), k_map(), m_map(), 
    midpoint(0, 0, 0) {

  m = msh;
  mdl = m->getModel();
  c_base = c;
  e_list = el_in;

  lens_split = new vd_lens(m, c, f_calc);
  bipy_split = new vd_bipy(m, c, f_calc);
  tet_split = new vd_sp_tet(m, c, f_calc);

  lens_split->turn_save_vtk(save_vtk_sub);
  lens_split->set_files("./output/sub/collapse", "./tempmesh/sub/temp");

  bipy_split->turn_save_vtk(save_vtk_sub);
  bipy_split->set_files("./output/sub/splitpy", "./tempmesh/sub/temp");
  tet_split->turn_save_vtk(save_vtk_sub);
  tet_split->set_files("./output/sub/splittet", "./tempmesh/sub/temp");

}

void vd_glens::l_clear() {

  for(int i = 0; i < ent.size(); i++)
    ent.at(i).clear();
  ent.clear();

  for(int i = 0; i < v_list.size(); i++) {
    for(int j = 0; j < v_list.at(i).size(); j++) {
      v_list.at(i).at(j).clear();
    }
    v_list.at(i).clear();
  }
  v_list.clear();

  for(int i = 0; i < m_list.size(); i++)
    m_list.at(i).clear();
  m_list.clear();

  e_map.clear();
  v_map.clear();
  m_map.clear();
  k_map.clear();

  ext_c0.clear();
  ext_c0_vec.clear();
  ext_cell = false;

  e_col_ext.clear();
}

// Set the output file prefixes.
void vd_glens::set_files(const char* vtkFile, const char* meshFile) {
  vtk_name = vtkFile;
  smb_name = meshFile;
}

// Set 
std::pair<int, int> vd_glens::set_cell(int dim, int cell_id) {

  cell_col_dim = dim;
  cell_col_id = cell_id;
  mdl_col = m->findModelEntity(cell_col_dim, cell_col_id);

  if(c_base->is_free(cell_col_dim, cell_col_id-1)) {
    c_flag = false;
  }
  else
    c_flag = true;

  c0_flag = false;

  l_clear();

  struct ent_conn e0;
  c_base->get_conn_dim_gmi(0, dim, cell_id, &e0);


  ext_c0.reserve(e0.conn.size());

  for(int i = 0; i < e0.conn.size(); i++) {
    if(c_base->get_cell_ext_gmi(0, e0.conn.at(i)) ) {
      ext_cell = true;
      ext_c0.push_back(e0.conn.at(i));
    }
  }

  ext_shell* e_sh;
  sh_min.dim = -1;
  if(upd_shell and f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
    e_sh = f_calc->get_e_sh();
    if(ext_cell) {
      ent_conn e_conn;
      c_base->get_conn_dim_gmi(0, dim, cell_id, &e_conn);

      for(int i = 0; i < e_conn.conn.size(); i++) {
        if(c_base->get_cell_ext_gmi(0, e_conn.conn.at(i))) {
          shell sh_0curr = e_sh->get_shell(0, e_conn.conn.at(i)-1);
          if(sh_min.dim == -1 or sh_0curr.dim < sh_min.dim)
            sh_min = sh_0curr;
        }
      }
    }
  }

  if(load_cell()) {
    return std::make_pair(dim, cell_id);
  }
  else {
    return std::make_pair(-2, -1);
  }
}


// After the collapse, smooth the boundaries by moving the vertices at the ends
// of the edges connected to the central vertex.
// TODO this step can use conjugate gradient approach with energy and energy 
// dissipation rate calculations.
void vd_glens::evolve_boundaries() {
  std::vector<apf::MeshEntity*> ee(0);
  std::vector<apf::MeshEntity*> ev(0);
  apf::Vector3 temp(0,0,0);
  apf::Field* vel_field = m->findField("velocity_field");

  vd_save_vtk_vert(m, vert_ctr, "output/collapsed_before");

  if (f_calc != NULL) {
    f_calc->vd_att_fields(m, vert_ctr);
  }

  vd_set_up(m, vert_ctr, &ee);
  ev.resize(ee.size());

  for(int i = 0; i < ee.size(); i++) {
    ev.at(i) = apf::getEdgeVertOppositeVert(m, ee.at(i), vert_ctr);
  }

  if (f_calc != NULL) {
    for(int i = 0; i < ev.size(); i++) {
      f_calc->vd_att_fields(m, ev.at(i));
    }
  }

  double mag_st = -1;
  double mag_temp = -1;
  for (int i = 0; i < ev.size(); i++) {
    f_calc->vd_upd_vel_field(m, ev.at(i));
    apf::getVector(vel_field, ev.at(i), 0, temp);
      mag_temp = temp.getLength();
    if(mag_st < -std::numeric_limits<double>::min() or mag_st < mag_temp)
      mag_st = mag_temp;
  }
  mag_temp = mag_st;
  double dt_st = f_calc->find_min_t2(m, &ev);
  double dt_inv = dt_st/9;

  //double dt_tot = 0;

  // Relaxation is terminated either if the maximum magnitude of the velocities is
  // below mag_st/rat_term, time step that can cause an inversion is below 
  // dt_st/rat_dt.
  double rat_term = 10;
  double rat_dt = 10;
  int iter_sz = 100;
  int iter = 0;
  while(iter < iter_sz and mag_temp > mag_st/rat_term and dt_inv > dt_st/rat_dt) {

    apf::Vector3 pos_curr(0,0,0);

    //dt_tot = dt_tot + dt_inv;
    for (int i = 0; i < ev.size(); i++) {
      m->getPoint(ev.at(i), 0, pos_curr);
      apf::getVector(vel_field, ev.at(i), 0, temp);
      pos_curr = pos_curr + temp*dt_inv;
      m->setPoint(ev.at(i), 0, pos_curr);
    }

    mag_temp = -1;
    for (int i = 0; i < ev.size(); i++) {
      f_calc->vd_upd_vel_field(m, ev.at(i));
      apf::getVector(vel_field, ev.at(i), 0, temp);
      double mag_curr = temp.getLength();
      if(mag_temp < -std::numeric_limits<double>::min() or mag_temp < mag_curr)
        mag_temp = mag_curr;
    }
    dt_inv = f_calc->find_min_t2(m, &ev)/9;
    if(dt_inv < -std::numeric_limits<double>::min())
      dt_inv = dt_st/9;
    iter = iter + 1;
  }
  vd_save_vtk_vert(m, vert_ctr, "output/collapsed_after");
}

// Return -1 if no cell is collapsed.
std::pair<int, int> vd_glens::col_cell(int dim, int cell_id, int tag_0c_in) {

  //vd_pert_cell_equi(m, dim, cell_id);
  //vd_shr_cell(m, dim, cell_id, float sp_size);
  //vd_adapt(m);
  if(vrfy_msh)
    m->verify();

  if(tag_0c_in != -1)
    tag_0c = tag_0c_in;
  else
    tag_0c = -1;


  std::pair<int, int> cell_coll = set_cell(dim, cell_id);
  if(cell_coll.first == -2)
    return cell_coll;
  //if(inv_flag)
  //  assert(!chk_inversion());

  // Recreates the merged and surrounding entities and destroys old entities.
  //get_vert();
  //destroy_ent();
  recreate_ent();
  m->acceptChanges();
  f_calc->refresh_mesh();

  vd_mesh_bad_adj(m);
  assert(vd_chk_neg(m) == 0);

  if(mod_c_base) {
    if(c0_flag) {
      ext_shell* e_sh;
      if(upd_shell and f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL)
        e_sh = f_calc->get_e_sh();

      if(ext_cell) {
        if(upd_shell and f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
          shell sh_false;
          for(int i = 1; i < cell_col_dim; i++) {
            sh_false.dim = i;
            sh_false.dim = 1;

            ent_conn e_conn;
            c_base->get_conn_dim_gmi(i, cell_col_dim, cell_col_id, &e_conn);
            for(int j = 0; j < e_conn.conn.size(); j++) {
              std::cout << "Removed shell of " << i << "c" << e_conn.conn.at(j)
                        << std::endl;
              e_sh->set_shell(i, e_conn.conn.at(j)-1, sh_false, false);
            }
          }
        }
      }

      c_base->coll_cell_gmi(cell_col_dim, cell_col_id, tag_vertex);

      if(ext_cell) {
        c_base->set_ext_gmi(0, tag_vertex, true);
        if(upd_shell and f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL) {
          std::cout << "0c" << tag_vertex
                    << " " << sh_min.dim << "sh" << sh_min.id
                    << std::endl;
          e_sh->set_shell(0, tag_vertex-1, sh_min);
          f_calc->corr_pos(m, vert_ctr);
        }
      }
    }
    else
      c_base->coll_cell_gmi(cell_col_dim, cell_col_id, -1);
  }

  m->acceptChanges();
  //m->verify();

  if (save_vtk) {
    vd_rem_tag(m);
    vd_tag_mesh(m);

    std::stringstream ss;
    ss << vtk_name << dim << "cell" << cell_id;
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    apf::writeVtkFiles(cstr, m);
  }

  //m->acceptChanges();
  //m->verify();
  if(mod_c_base) {
    vd_cell_val val_vert(m, c_base, vert_ctr);
    assert(val_vert.vd_vert_valid());
  }

  evolve_boundaries();

  if(c0_flag)
    return std::make_pair(0, tag_vertex);
  else
    return std::make_pair(-1, tag_vertex);
}

std::pair<int, int> vd_glens::get_cell() {
  return std::make_pair(type_vertex, tag_vertex);
}

void vd_glens::save_glens_vtk(const char* vtk_filename) {

  vd_rem_tag(m);
  vd_tag_mesh(m);

  //vd_tag_glens();
  apf::writeVtkFiles(vtk_filename, m);
}

void vd_glens::turn_save_vtk(bool onoff) {
  save_vtk = onoff;
}

void vd_glens::turn_save_vtk_sub(bool onoff) {
  save_vtk_sub = onoff;
}

void vd_glens::set_verify(bool onoff) {
  vrfy_msh = onoff;
}

void vd_glens::set_precond(bool onoff) {
  precond_flag = onoff;
}

// Used after collapsing internal structure within the lens to get the new 
// vertex.
apf::MeshEntity* vd_glens::get_vert_ctr() {
  return vert_ctr;
}

void vd_glens::createfields() {
  ch_f = vd_att_vs_field(m, "cvx");
  sp_f = vd_att_vs_field(m, "v_sp");
  f_field = vd_att_vv_field(m, "f_restore");
}

void vd_glens::destroyfields() {
  assert(ch_f != NULL);
  assert(sp_f != NULL);
  assert(f_field != NULL);

  apf::destroyField(ch_f);
  apf::destroyField(sp_f);
  apf::destroyField(f_field);

  ch_f = NULL;
  sp_f = NULL;
  f_field = NULL;
}

void vd_glens::set_field_calc(field_calc* f_calc_in, bool on_off) {
  upd_shell = on_off;
  f_calc = f_calc_in;
  lens_split->set_field_calc(f_calc);
  bipy_split->set_field_calc(f_calc);
  tet_split->set_field_calc(f_calc);
}

void vd_glens::set_inv(bool on_off) {
  inv_flag = on_off;
}

//void vd_glens::calc_sp_pos(bool on_off) {
//  sp_trial_load = on_off;
//}

// Used in attaching arbitrary tags.
/*
// Given a mesh object, create element and element volume membership tags.
// findNumbering is problematic, so destroy and rewrite.
void vd_glens::vd_tag_glens() {
  int meshDimension = m->getDimension();

  apf::Numbering* tagnumbering;
  if (m->findNumbering("lens_memb")) {
    printf("lens_memb tag already exists.\n");
    tagnumbering = m->findNumbering("lens_memb");
    apf::destroyNumbering(tagnumbering);
    tagnumbering = apf::createNumbering(m, "lens_memb", apf::getConstant(meshDimension),1);
  }
  else {
    tagnumbering = apf::createNumbering(m, "lens_memb", apf::getConstant(meshDimension),1);
  }

  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(meshDimension);

  while ((e = m->iterate(it)))
  {
    apf::number(tagnumbering, e, 0, 0, 0);
  }
  m->end(it);

  for (int i = 0; i < edge_nbr; i++) {
    for (int j = 0; j < elem_col.at(i).size(); j++)
      apf::number(tagnumbering, elem_col.at(i).at(j), 0, 0, 1);
  }


  if (m->findNumbering("surr_memb")) {
    printf("surr_memb tag already exists.\n");
    tagnumbering = m->findNumbering("surr_memb");
    apf::destroyNumbering(tagnumbering);
    tagnumbering = apf::createNumbering(m, "surr_memb", apf::getConstant(meshDimension),1);
  }
  else {
    tagnumbering = apf::createNumbering(m, "surr_memb", apf::getConstant(meshDimension),1);
  }

  it = m->begin(meshDimension);

  while ((e = m->iterate(it)))
  {
    apf::number(tagnumbering, e, 0, 0, 9);
  }
  m->end(it);

  for (int i = 0; i < surr_elem.size(); i++) {
    apf::number(tagnumbering, surr_elem.at(i), 0, 0, 7);
  }

  for (int i = 0; i < edge_nbr; i++) {
    for (int j = 0; j < elem_col.at(i).size(); j++)
      apf::number(tagnumbering, elem_col.at(i).at(j), 0, 0, i+1);
  }

	if (ent_d == 3)
    apf::number(tagnumbering, ent_col, 0, 0, 0);

}
void vd_glens::vd_start_tag() {
  apf::Numbering* tagnumbering;
  if (m->findNumbering("tag_ent")) {
    printf("tag_ent tag already exists.\n");
    tagnumbering = m->findNumbering("tag_ent");
    apf::destroyNumbering(tagnumbering);
    tagnumbering = apf::createNumbering(m, "tag_ent", apf::getConstant(3),1);
  }
  else {
    tagnumbering = apf::createNumbering(m, "tag_ent", apf::getConstant(3),1);
  }

  tag_flag = true;

}

bool vd_glens::vd_tag_ent_glens(apf::MeshEntity* ent, int tag) {
  if (tag_flag) {
    apf::Numbering* tagnumbering = m->findNumbering("tag_ent");
    apf::number(tagnumbering, ent, 0, 0, tag);
    return true;
  }

  return false;
}
*/
// Set the simulation parameters.
void vd_glens::set_vdpar(vd_param par_in) {
  vd_par = par_in;
}

void vd_glens::set_mod_c_base(bool on_off) {
  mod_c_base = on_off;
}

vd_glens::~vd_glens() {
  l_clear();
  delete lens_split;
  delete bipy_split;
  delete tet_split;
}


// vd_glens_trial
vd_glens_trial::vd_glens_trial(apf::Mesh2* msh, struct cell_base* c, vd_entlist* el_in, field_calc* f_calc_in) :  m(NULL), m_trial(NULL), c_base(NULL),
                      e_list(NULL), f_calc(NULL), trial_load(false),
                      precond_flag(true), vrfy_msh(true),
                      c_id(-1), c_dim(-1), es(0), 
                      es_ent(0, std::vector<apf::MeshEntity*> (0)), vert_t(0)
  {

  m = msh;
  e_list = el_in;
  trial_load = false;
  f_calc = f_calc_in;

  es_ent.resize(4);

  clear();

  c_base = c;
}

vd_glens_trial::~vd_glens_trial() {
  clear();
  clear_ent();
  es_ent.clear();
}

void vd_glens_trial::set_cell(int cell_dim, int cell_id) {

  c_dim = cell_dim;
  c_id = cell_id;

  assert(c_dim > 0 and c_dim < 4);

  es = e_list->e.at(cell_dim).at(cell_id-1).at(cell_dim);
  //vd_find_ent_geom(m, &es, cell_id, cell_dim, cell_dim);

  assert(es.size() > 0);
  
  clear();
  //m_precond = apf::makeEmptyMdsMesh(m->getModel(), 3, false);
  gmi_register_null();
  m_trial = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);

  copy_ent();

  trial_load = true;
}

std::pair<int, int> vd_glens_trial::coll_cell(int cell_dim, int cell_id, int tag_0c_in) {
  set_cell(cell_dim, cell_id);
  vd_entlist e_list_trial(m_trial, c_base);

  vd_glens g_lens(m_trial, c_base, &e_list_trial);
  g_lens.set_field_calc(f_calc);
  g_lens.set_verify(false);
  g_lens.set_precond(precond_flag);
  f_calc->vd_att_fields(m_trial);

  f_calc->refresh(m_trial, c_base, &e_list_trial);

  std::pair<int, int> cell_new = std::make_pair(-2, 0);

  std::cout << "Trial collapse" << cell_new.first 
                                << "c" << cell_new.second << std::endl;
  g_lens.set_mod_c_base(false);
  g_lens.set_verify(vrfy_msh);
  cell_new = g_lens.col_cell(cell_dim, cell_id, tag_0c_in);

  f_calc->refresh(m, c_base, e_list);
  if(cell_new.first != -2) {
    std::cout << "Cell collapse" << cell_new.first 
                                  << "c" << cell_new.second << std::endl;

    vd_glens g_lens2(m, c_base, e_list);
    g_lens2.set_field_calc(f_calc);
    g_lens2.set_precond(precond_flag);
    cell_new = g_lens2.col_cell(cell_dim, cell_id, tag_0c_in);
    //g_lens2.get_v_sp_pos(std::map<MeshEntity*, apf::MeshEntity*> &e2e_map_in, 
    //                  std::map<MeshEntity*, apf::Vector3> &e_sp_pos_trial_in);
  }
  return cell_new;
}
void vd_glens_trial::set_precond(bool on_off) {
  precond_flag = on_off;
}

void vd_glens_trial::set_verify(bool on_off) {
  vrfy_msh = on_off;
}

// Copy entities from the original mesh.
// Copy entities from the original mesh to the preconditioned mesh.
// This one uses buildelement, apf implementation with actual mds support.
apf::MeshEntity* vd_glens_trial::copyElement_pre(
    apf::MeshEntity* original,
    apf::BuildCallback* cb) {
  // Get the type in the original mesh.
  int type = m->getType(original);
  if (type==0)  {
    int e1 = findIn(&es_ent.at(0), es_ent.at(0).size(), original);
    assert(e1 > -1);
    // Return the corresponding vertex in trial mesh:
    return vert_t.at(e1);
  }

  int d = apf::Mesh::typeDimension[type];
  apf::Downward down;
  int nd = m->getDownward(original,d-1,down);

  //apf::Downward vert_down; // Downward vertices of the current entity.
  //m->getDownward(original,0,vert_down);

  for (int i=0; i < nd; ++i)
    down[i] = copyElement_pre(down[i],cb);

  nd = m->getDownward(original,0,down);
  for (int i=0; i < nd; ++i) {
    int e1 = findIn(&es_ent.at(0), es_ent.at(0).size(), down[i]);
    assert(e1 > -1);
    // Return the corresponding vertex in trial mesh:
    down[i] = vert_t.at(e1);
  }

  // The topology object in the actual mesh:
  apf::ModelEntity* mdl = m->toModel(original);
  int type_ori = m->getModelType(mdl);
  int tag_ori = m->getModelTag(mdl);

  //std::cout << type << "-Ent " << original 
  //          << " " << type_ori << "c" << tag_ori
  //          << std::endl;
  // The topology object in the trial mesh:
  mdl = m_trial->findModelEntity(type_ori,tag_ori);

  apf::MeshEntity* ent;
  ent = buildElement(m_trial, mdl, type, down);

  //vd_print_vert(m_precond, ent);

  return ent;
}

void vd_glens_trial::copy_ent() {

  apf::Downward v;

  es_ent.at(c_dim) = es;

  for(int i = c_dim; i > 0; i--) {
    vd_set_down(m, &es_ent.at(i), &es_ent.at(i-1));
  }

  for(int i = 0; i < 3; i++) {
    vd_set_up(m, &es_ent.at(i), &es_ent.at(i+1));
  }

  for(int i = 3; i > 0; i--) {
    vd_set_down(m, &es_ent.at(i), &es_ent.at(i-1));
  }

  // Save the entities before the transfer.
  vd_rem_tag(m);
  vd_tag_mesh(m);
  vd_tag_set(m, &es_ent.at(3), "Vertex_ent");

  apf::writeVtkFiles("./output/before_col_trial", m);

  // These are to be copied to the new trial mesh:
  vert_t.resize(es_ent.at(0).size());

  // Copy the vertices to the new trial mesh.
  for (int i = 0; i < es_ent.at(0).size(); i++) {
    //std::cout << "Vert " << es_vert.at(i) << "["
    //          << m_main->toModel(es_vert.at(i)) << "] ," ;
    int type = m->getModelType(m->toModel(es_ent.at(0).at(i)));
    int tag = m->getModelTag(m->toModel(es_ent.at(0).at(i)));
    apf::ModelEntity* vert_em = m_trial->findModelEntity(type,tag);
    vert_t.at(i) = m_trial->createVert(vert_em);

    apf::Vector3 vert_pos(0,0,0);
    m->getPoint(es_ent.at(0).at(i), 0, vert_pos);
    //std::cout<< vert_pos << std::endl;
    m_trial->setPoint(vert_t.at(i), 0, vert_pos);
  }

  // Copy the elements recursively:
  for (int i = 0; i < es_ent.at(3).size(); i++) {
    copyElement_pre(es_ent.at(3).at(i),0);
  }

  m_trial->acceptChanges();
}

void vd_glens_trial::clear_ent() {
  for(int i = 0; i < es_ent.size(); i++) {
    es_ent.at(i).clear();
  }
}

void vd_glens_trial::clear() {
  if(trial_load)
    apf::destroyMesh(m_trial);
  trial_load = false;
  clear_ent();
}


/*



double vd_glens_trial::get_energy() {

  double e_tot = 0;

  apf::MeshEntity* ent;
  apf::MeshIterator* it = m_trial->begin(3);
  while(ent = m_trial->iterate(it)) {
    e_tot = e_tot + vd_e.en_tet(m_trial, ent);
  }
  m_trial->end(it);

  it = m_trial->begin(2);
  while(ent = m_trial->iterate(it)) {
    e_tot = e_tot + vd_e.en_tri(m_trial, ent);
  }
  m_trial->end(it);

  return e_tot;
}

void vd_glens_trial::clear_ent() {
  for(int i = 0; i < es_ent.size(); i++) {
    es_ent.at(i).clear();
  }
}

void vd_glens_trial::clear() {
  if(trial_load)
    apf::destroyMesh(m_trial);
  clear_ent();
}


void vd_glens_trial::reload(apf::Mesh2* msh, struct cell_base* c, 
                                                    vd_entlist* el_in) {

  m = msh;

  trial_load = false;

  c_base_in = c;
  *c_base = *c;
  e_list = el_in;
}


*/

// vd_glens_trial
vd_glens_trial_ins::vd_glens_trial_ins(apf::Mesh2* msh, struct cell_base* c, vd_entlist* el_in, field_calc* f_calc_in) :  m(NULL), m_trial(NULL), c_base(NULL),
                      e_list(NULL), f_calc(NULL), trial_load(false),
                      precond_flag(true), 
                      c_id(-1), c_dim(-1), es(0), 
                      es_ent(0, std::vector<apf::MeshEntity*> (0)), vert_t(0)
  {

  m = msh;
  e_list = el_in;
  trial_load = false;
  f_calc = f_calc_in;

  es_ent.resize(4);

  clear();

  c_base = c;
}

vd_glens_trial_ins::~vd_glens_trial_ins() {
  clear();
  clear_ent();
  es_ent.clear();
}

void vd_glens_trial_ins::set_cell(int cell_dim, int cell_id) {

  c_dim = cell_dim;
  c_id = cell_id;

  assert(c_dim > 0 and c_dim < 4);

  es = e_list->e.at(cell_dim).at(cell_id-1).at(cell_dim);
  //vd_find_ent_geom(m, &es, cell_id, cell_dim, cell_dim);

  assert(es.size() > 0);
  
  clear();
  //m_precond = apf::makeEmptyMdsMesh(m->getModel(), 3, false);
  gmi_register_null();
  m_trial = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);

  copy_ent();

  trial_load = true;
}

std::pair<int, int> vd_glens_trial_ins::coll_cell(int cell_dim, int cell_id, int tag_0c_in) {
  set_cell(cell_dim, cell_id);
  vd_entlist e_list_trial(m_trial, c_base);

  vd_glens g_lens(m_trial, c_base, &e_list_trial);
  g_lens.set_field_calc(f_calc);
  g_lens.set_verify(false);
  g_lens.set_precond(precond_flag);
  f_calc->vd_att_fields(m_trial);

  f_calc->refresh(m_trial, c_base, &e_list_trial);

  std::pair<int, int> cell_new = std::make_pair(-2, 0);

  std::cout << "Trial collapse" << cell_new.first 
                                << "c" << cell_new.second << std::endl;
  g_lens.set_mod_c_base(false);
  cell_new = g_lens.col_cell(cell_dim, cell_id, tag_0c_in);

  f_calc->refresh(m, c_base, e_list);
  if(cell_new.first != -2) {
    std::cout << "Cell collapse" << cell_new.first 
                                  << "c" << cell_new.second << std::endl;

    vd_glens g_lens2(m, c_base, e_list);
    g_lens2.set_field_calc(f_calc);
    g_lens2.set_precond(precond_flag);
    cell_new = g_lens2.col_cell(cell_dim, cell_id, tag_0c_in);
    //g_lens2.get_v_sp_pos(std::map<MeshEntity*, apf::MeshEntity*> &e2e_map_in, 
    //                  std::map<MeshEntity*, apf::Vector3> &e_sp_pos_trial_in);
  }
  return cell_new;
}
void vd_glens_trial_ins::set_precond(bool on_off) {
  precond_flag = on_off;
}

// Copy entities from the original mesh.
// Copy entities from the original mesh to the preconditioned mesh.
// This one uses buildelement, apf implementation with actual mds support.
apf::MeshEntity* vd_glens_trial_ins::copyElement_pre(
    apf::MeshEntity* original,
    apf::BuildCallback* cb) {
  // Get the type in the original mesh.
  int type = m->getType(original);
  if (type==0)  {
    int e1 = findIn(&es_ent.at(0), es_ent.at(0).size(), original);
    assert(e1 > -1);
    // Return the corresponding vertex in trial mesh:
    return vert_t.at(e1);
  }

  int d = apf::Mesh::typeDimension[type];
  apf::Downward down;
  int nd = m->getDownward(original,d-1,down);

  //apf::Downward vert_down; // Downward vertices of the current entity.
  //m->getDownward(original,0,vert_down);

  for (int i=0; i < nd; ++i)
    down[i] = copyElement_pre(down[i],cb);

  nd = m->getDownward(original,0,down);
  for (int i=0; i < nd; ++i) {
    int e1 = findIn(&es_ent.at(0), es_ent.at(0).size(), down[i]);
    assert(e1 > -1);
    // Return the corresponding vertex in trial mesh:
    down[i] = vert_t.at(e1);
  }

  // The topology object in the actual mesh:
  apf::ModelEntity* mdl = m->toModel(original);
  int type_ori = m->getModelType(mdl);
  int tag_ori = m->getModelTag(mdl);

  //std::cout << type << "-Ent " << original 
  //          << " " << type_ori << "c" << tag_ori
  //          << std::endl;
  // The topology object in the trial mesh:
  mdl = m_trial->findModelEntity(type_ori,tag_ori);

  apf::MeshEntity* ent;
  ent = buildElement(m_trial, mdl, type, down);

  //vd_print_vert(m_precond, ent);

  return ent;
}

void vd_glens_trial_ins::copy_ent() {

  apf::Downward v;

  es_ent.at(c_dim) = es;

  for(int i = c_dim; i > 0; i--) {
    vd_set_down(m, &es_ent.at(i), &es_ent.at(i-1));
  }

  for(int i = 0; i < 3; i++) {
    vd_set_up(m, &es_ent.at(i), &es_ent.at(i+1));
  }

  for(int i = 3; i > 0; i--) {
    vd_set_down(m, &es_ent.at(i), &es_ent.at(i-1));
  }

  // Save the entities before the transfer.
  vd_rem_tag(m);
  vd_tag_mesh(m);
  vd_tag_set(m, &es_ent.at(3), "Vertex_ent");

  apf::writeVtkFiles("./output/before_col_trial", m);

  // These are to be copied to the new trial mesh:
  vert_t.resize(es_ent.at(0).size());

  // Copy the vertices to the new trial mesh.
  for (int i = 0; i < es_ent.at(0).size(); i++) {
    //std::cout << "Vert " << es_vert.at(i) << "["
    //          << m_main->toModel(es_vert.at(i)) << "] ," ;
    int type = m->getModelType(m->toModel(es_ent.at(0).at(i)));
    int tag = m->getModelTag(m->toModel(es_ent.at(0).at(i)));
    apf::ModelEntity* vert_em = m_trial->findModelEntity(type,tag);
    vert_t.at(i) = m_trial->createVert(vert_em);

    apf::Vector3 vert_pos(0,0,0);
    m->getPoint(es_ent.at(0).at(i), 0, vert_pos);
    //std::cout<< vert_pos << std::endl;
    m_trial->setPoint(vert_t.at(i), 0, vert_pos);
  }

  // Copy the elements recursively:
  for (int i = 0; i < es_ent.at(3).size(); i++) {
    copyElement_pre(es_ent.at(3).at(i),0);
  }

  m_trial->acceptChanges();
}

void vd_glens_trial_ins::clear_ent() {
  for(int i = 0; i < es_ent.size(); i++) {
    es_ent.at(i).clear();
  }
}

void vd_glens_trial_ins::clear() {
  if(trial_load)
    apf::destroyMesh(m_trial);
  trial_load = false;
  clear_ent();
}

