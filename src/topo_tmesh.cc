#include <vector>
#include <map>
#include <deque>
#include <algorithm>    
#include <functional>

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

#include "topo_topo.h"
#include "topo_graph.h"

#include "topo_entlist.h"

#include "topo_tmesh.h"

// The cell detector object that works with expansions around 0cells with
// disjoint entity sets belonging to the same cell. This extracts the disjoint
// cells and tags them uniquely and extracts a 2/3-cell adjacency graph. 
vd_cell_det::vd_cell_det() : 
    c0_tag(0), c_dim(1), d_ent(1), id(0),
    calc_ext(false), calc_corner(false), ext_0cell(false), ext_corner(false),
    ext_2cell{}, ext_sz(0), c3_nbr(0), cb_write_flag(true),
    cg(), ent_list(),
    ng_local(),
    vert_c1_map{},
    c2_c3real1(0), c2_c3real2(0),
    v_0cell(NULL), mdl_curr(NULL),
    cells(3, std::vector<int > (0) ),
    adj(2, std::vector<std::vector<int > > (0, std::vector<int> (0) ) ),
    e_2c(0), c3_edge{}, c2_dir{}, c3_dir{}, circ_tup(0), path(0),
    ent_sz(4), burn_mdl(0, std::map<int, apf::ModelEntity*>{} ), burn_ent{},
    ent_cell{}, ent_ofire(2, std::vector<apf::MeshEntity* >(0)) {

  cb_dj = new cell_base(0, 0, 0, 0);
  //cb_dj->print_ent();
}

vd_cell_det::vd_cell_det(apf::Mesh2* m_in, struct cell_base* cb_in, apf::MeshEntity* v_0cell_in) : 
    c0_tag(0), c_dim(1), d_ent(1), id(0),
    calc_ext(false), calc_corner(false), ext_0cell(false), ext_corner(false),
    ext_2cell{}, ext_sz(0), c3_nbr(0), cb_write_flag(true),
    cg(), ent_list(),
    vert_c1_map{},
    c2_c3real1(0), c2_c3real2(0),
    v_0cell(NULL), mdl_curr(NULL),
    cells(3, std::vector<int > (0) ),
    adj(2, std::vector<std::vector<int > > (0, std::vector<int> (0) ) ),
    e_2c(0), c3_edge{}, c2_dir{}, c3_dir{}, circ_tup(0), path(0),
    ent_sz(4), burn_mdl(0, std::map<int, apf::ModelEntity*>{} ), burn_ent{},
    ent_cell{}, ent_ofire(2, std::vector<apf::MeshEntity* >(0)) {

  cb = cb_in;
  m = m_in;
  v_0cell = v_0cell_in;
  assert(m->getModelType(m->toModel(v_0cell)) == 0);
  c0_tag = m->getModelTag(m->toModel(v_0cell));

  calc_ext = false;
  calc_corner = false;
  cb_write_flag = true;

  ent_sz.resize(4);

  cb_dj = new cell_base(0, 0, 0, 0);

  ent_list.change_mesh(m, v_0cell, cb);
  //ent_list.print();
  load_graph();
}

// TODO CAN two disjoint sets of the same 2stratum can cause problems...
// vert_c1_map associates the new vertices on the trial mesh with the ones already
// associated on the preconditioned mesh... So they should be matching
// Using the ids of the bounding disjoint cells, return the id of dim-cell
int vd_cell_det::find_cell_dj(int dim, apf::MeshEntity* ent, std::vector<int> bound) {
  ent_conn e_con;

  for (int c_id = 0; c_id < cb_dj->get_sz(dim); c_id++) {
    bool found = true;
    cb_dj->get_conn(dim, c_id, &e_con);

    if(e_con.conn.size() == bound.size()) {
      std::sort(e_con.conn.begin(), e_con.conn.end());
      for (int i = 0; i < e_con.conn.size(); i++) {
        if(e_con.conn.at(i) != bound.at(i)) {
          found = false;
          i = e_con.conn.size();
        }
      }
      if(found) {

        if(dim == 2) {
          apf::Up up;
          m->getUp(ent, up);
          assert(up.n > 0);
          assert(c2_c3real1.at(c_id) > 0);
          int c3_1 = c2_c3real1.at(c_id);
          int c3_2 = c2_c3real2.at(c_id);
          //std::cout << "c3_1 " << c3_1 << " " << " c3_2 " << c3_2 
          //          << std::endl;
          //for (int i = 0; i < up.n; i++) {
          //  std::cout << m->getModelTag(m->toModel(up.e[i])) << " ";
          //}
          //std::cout << std::endl;
          if(c3_2) {
            if(up.n == 2) {
              int c3_r1 = m->getModelTag(m->toModel(up.e[0]));
              int c3_r2 = m->getModelTag(m->toModel(up.e[1]));
              //std::cout << "c3_r1 " << c3_r1 << " "
              //          << "c3_r2 " << c3_r2 << " "
              //          << std::endl;
              if((c3_1 == c3_r1 and c3_2 == c3_r2) or
                 (c3_1 == c3_r2 and c3_2 == c3_r1) ) {
                return c_id;
              }
            }
          }
          else {
            if(up.n == 1) {
              int c3_r = m->getModelTag(m->toModel(up.e[0]));
              if(c3_1 == c3_r)
                return c_id;
            }
          }
        }
        else
          return c_id;
      }
    }
  }
  return -1;
}

// Going over disjoint sets of dim-dimensional entities, assign their lower 
// dimensional adjacencies with the same stratum the same id. 
void vd_cell_det::burn_lower() {
  ent_conn e_con;
  apf::Downward d_e;
  apf::Downward d_v;
  std::vector<apf::MeshEntity*> ents(0);
  for (int dim = 2; dim < 4; dim++) {
    for (int c_id = 0; c_id < cb_dj->get_sz(dim); c_id++) {
      apf::ModelEntity* mdl = get_mdl(dim, c_id);
      ents = get_dim_ent(dim, c_id);
      for (int i = 0; i < ents.size(); i++) {
        for (int dim_d = 1; dim_d < dim; dim_d++) {
          m->getDownward(ents.at(i), dim_d, d_e);
          for (int j = 0; j < dim+1; j++) {
            apf::ModelEntity* mdl_curr = m->toModel(d_e[j]);
            if(mdl == mdl_curr) {
              m->getDownward(ents.at(i), 0, d_v);
              int i1 = findIn(d_v, dim+1, v_0cell);
              if(i1 > -1) {
                ent_cell[ents.at(i)] = c_id;
              }
            }
          }
        }
      }
    }
  }
}

// TODO cg seems to assume the cells list to be ordered. This is a bug, but sort
// it out for now by reordering the cells before feeding into the cell_graph 
// object.
void vd_cell_det::reorder_cells() {
  std::vector<int> adj_temp(0);
  for (int dim = 2; dim < 4; dim++) {
    for (int i = 0; i < cells.at(dim-1).size(); i++) {
      if(cells.at(dim-1).at(i) != i) {
        int i_act = cells.at(dim-1).at(i);
        adj_temp = adj.at(dim-2).at(i_act);
        adj.at(dim-2).at(i_act) = adj.at(dim-2).at(i);
        adj.at(dim-2).at(i) = adj_temp;
        cells.at(dim-1).at(i_act) = cells.at(dim-1).at(i);
        cells.at(dim-1).at(i) = i;
        i = i - 1;
      }
    }
  }
}


// Return the disjoint set id of the entity.
int vd_cell_det::get_id(apf::MeshEntity* ent) {
  return ent_cell[ent];
}

void vd_cell_det::burn_new() {

  apf::MeshEntity* ent;
  cb_dj->clear();
  struct ent_conn e_con;

  cb_dj->add_free(0, 1);
  int c0_id = cb_dj->use_free(0);
  assert(c0_id == 0);

  c_dim = 1;
  d_ent = 1;

  id = 0;

  vert_c1_map.clear();

  apf::Downward down;
  apf::MeshEntity* v_curr;

  for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
    for (int i = 0; i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
      ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

      if(!burn_ent[ent]) {
        burn_cell(c_dim, id, ent);
        cb_dj->add_free(1,1);
        int c1_id = cb_dj->use_free(1);
        assert(id == c1_id);

        //std::cout << ent << " " 
        //        << vd_get_pos(m, ent) - vd_get_pos(m, v_0cell)
        //        << std::endl;
        //std::cout << "1c" << c_id+1 << std::endl;
        //std::cout << "1c_dj" << id << std::endl;

        e_con.clear();
        e_con.conn.push_back(c0_id);
        cb_dj->set_conn(1, c1_id, &e_con);

        m->getDownward(ent, 0, down);
        if(down[0] == v_0cell)
          v_curr = down[1];
        else
          v_curr = down[0];
        //std::cout << m->getModelType(m->toModel(ent)) << "c"
        //          << m->getModelTag(m->toModel(ent))
        //          << " " << v_curr << " " << id << std::endl;
        vert_c1_map[v_curr] = id + 1;

        id++;
      }
    }
  }

  c_dim = 2;
  d_ent = 2;

  id = 0;

  for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
    for (int i = 0; i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
      ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

      if(!burn_ent[ent]) {
        burn_cell(c_dim, id, ent);
        cb_dj->add_free(2,1);
        int c2_id = cb_dj->use_free(2);
        assert(id == c2_id);

        e_con.clear();
        e_con.conn = adj.at(c_dim-2).back();
        cb_dj->set_conn(2, c2_id, &e_con);

        //std::cout << "2c" << c_id + 1 << std::endl;
        //std::cout << "bound ";
        for (int j = 0; j < e_con.conn.size(); j++) {
          std::cout << e_con.conn.at(j) << " ";
        }
        std::cout << std::endl;

        //std::cout << "2c_dj" << id << std::endl;
        id++;
      }
    }
  }
  //coll_c2_dir();

  c_dim = 3;
  d_ent = 3;


  c2_c3real1.clear();
  c2_c3real2.clear();
  c2_c3real1.resize(cb_dj->get_sz(2));
  c2_c3real2.resize(cb_dj->get_sz(2));

  id = 0;

  //std::map<int, int> c2_c3map1;
  //std::map<int, int> c2_c3map2;
  // Loop over c3 tets. Burn disjoint groups of tets. Assuming there is a 
  // single c3 edge per disjoint group of entities, find the edge and 
  // associate it with the c3.
  for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
    int id_start = id;

    for (int i = 0; i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
      ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

      if(!burn_ent[ent]) {
        burn_cell(c_dim, id, ent);
        cb_dj->add_free(3,1);
        int c3_id = cb_dj->use_free(3);
        assert(id == c3_id);

        e_con.clear();
        e_con.conn = adj.at(c_dim-2).back();
        cb_dj->set_conn(3, c3_id, &e_con);

        //std::cout << "3c" << c_id + 1 << std::endl;
        //std::cout << "bound ";
        //for (int j = 0; j < e_con.conn.size(); j++) {
        //  std::cout << e_con.conn.at(j) << " ";
        //}
        //std::cout << std::endl;

        //std::cout << "3c_dj" << id << std::endl;

        int c3_real = m->getModelTag(m->toModel(ent));
        for(int j = 0; j < e_con.conn.size(); j++) {
/*
          if(!c2_c3map1[e_con.conn.at(j)]) {
            std::cout << "2c" << e_con.conn.at(j) << " 3c" << c3_real
                      << std::endl;
            c2_c3map1[e_con.conn.at(j)] = c3_real;
          }
          else {
            if(!c2_c3map2[e_con.conn.at(j)]) {
              std::cout << "2c" << e_con.conn.at(j) << " 3c" << c3_real
                        << std::endl;
              c2_c3map2[e_con.conn.at(j)] = c3_real;
            }
          }
*/
          if(!c2_c3real1.at(e_con.conn.at(j))) {
            //std::cout << "2c" << e_con.conn.at(j) << " 3c" << c3_real
            //          << std::endl;
            c2_c3real1.at(e_con.conn.at(j)) = c3_real;
          }
          else {
            assert(!c2_c3real2.at(e_con.conn.at(j)));
            //std::cout << "2c" << e_con.conn.at(j) << " 3c" << c3_real
            //          << std::endl;
            c2_c3real2.at(e_con.conn.at(j)) = c3_real;
          }
        }
        std::vector<apf::MeshEntity*>* tri_list;
        tri_list = &ent_list.e.at(c_dim).at(c_id).at(2);
        for(int j = 0; j < tri_list->size(); j++) {
          if(burn_ent[tri_list->at(j)] and ent_cell[tri_list->at(j)] == id ) {
            apf::Downward d_edge;
            apf::Downward d_vert;
            m->getDownward(tri_list->at(j), 1, d_edge);
            m->getDownward(tri_list->at(j), 0, d_vert);
            int v1 = findIn(d_vert, 3, v_0cell);
            assert(v1 > -1);
            int lu_tri_ed [3][2] = {{2,0},{0,1},{1,2}};
            if(m->toModel(d_edge[lu_tri_ed[v1][0]]) == m->toModel(ent)) {
              c3_edge[id] = d_edge[lu_tri_ed[v1][0]];
            }
            else {
              assert(m->toModel(d_edge[lu_tri_ed[v1][1]]) == m->toModel(ent));
              c3_edge[id] = d_edge[lu_tri_ed[v1][1]];
            }
            j = tri_list->size();
          }
        }

        id++;
      }
    }
    // If the c3 tets are joint, assert the number of c3 edges are less than
    // 2.
    if(id_start == id - 1)
      assert(ent_list.e.at(c_dim).at(c_id).at(1).size() < 2);
  }
}
//TODO Somehow the disjoint 3stratum and the 2stratum connections are switched when two disjoint 3strata of the same 3stratum are connected to the same 3stratum by two disjoint 2strata belonging to the same 2stratum...
void vd_cell_det::burn_from_cb() {

  apf::Downward down;
  apf::MeshEntity* v_curr;
  apf::MeshEntity* ent;

  std::vector<int> bound(0);

  c_dim = 1;
  d_ent = 1;

  for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
    for (int i = 0; i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
      ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

      if(!burn_ent[ent]) {
        m->getDownward(ent, 0, down);
        if(down[0] == v_0cell)
          v_curr = down[1];
        else
          v_curr = down[0];

        //std::cout << ent << " " 
        //        << vd_get_pos(m, ent) - vd_get_pos(m, v_0cell)
        //        << std::endl;

        id = vert_c1_map[v_curr] - 1;
        assert(id > -1);

        burn_cell(c_dim, id, ent);
        //std::cout << "1c" << c_id + 1 << std::endl;
        //std::cout << "1c_dj" << id << std::endl;

      }
    }
  }

  c_dim = 2;
  d_ent = 2;
/*
  for (int c_id = 0; c_id < cb_dj->get_sz(c_dim); c_id++) {
    std::cout << "c2_c3_real1 ";
    std::cout << c2_c3real1.at(c_id) << " ";
    std::cout << c2_c3real2.at(c_id) << " ";
    std::cout << std::endl;
    ent_conn e_con;
    cb_dj->get_conn(c_dim, c_id, &e_con);
    std::cout << "bound ";
    for (int j = 0; j < e_con.conn.size(); j++) {
      std::cout << e_con.conn.at(j) << " ";
    }
    std::cout << std::endl;
  }
*/
  for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
    for (int i = 0; i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
      ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

      if(!burn_ent[ent]) {
        std::vector<int> bound = burn_bound(c_dim, ent);

        apf::Up up;
        m->getUp(ent, up);

        int id = find_cell_dj(c_dim, ent, bound);
        assert(id > -1);
        burn_cell(c_dim, id, ent);

        //std::cout << "2c" << c_id + 1 << std::endl;
        //if(up.n > 0)
        //  std::cout << m->getModelTag(m->toModel(up.e[0])) << " ";
        //if(up.n == 2)
        //  std::cout << m->getModelTag(m->toModel(up.e[1])) << std::endl;

        //std::cout << "bound ";
        //for (int j = 0; j < bound.size(); j++) {
        //  std::cout << bound.at(j) << " ";
        //}
        //std::cout << std::endl;

        //std::cout << "2c_dj" << id << std::endl;

      }
    }
  }

  c_dim = 3;
  d_ent = 3;

  // Loop over c3 tets. Burn disjoint groups of tets. Assuming there is a 
  // single c3 edge per disjoint group of entities, find the edge and 
  // associate it with the c3.
  for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
    for (int i = 0; i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
      ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

      if(!burn_ent[ent]) {
        std::vector<int> bound = burn_bound(c_dim, ent);

        int id = find_cell_dj(c_dim, ent, bound);
        assert(id > -1);
        burn_cell(c_dim, id, ent);

        //std::cout << "3c" << c_id + 1 << std::endl;
        //std::cout << "bound ";
        //for (int j = 0; j < bound.size(); j++) {
        //  std::cout << bound.at(j) << " ";
        //}
        //std::cout << std::endl;
        //std::cout << "3c_dj" << id << std::endl;

        std::vector<apf::MeshEntity*>* tri_list;
        tri_list = &ent_list.e.at(c_dim).at(c_id).at(2);
        // Find the disjoint c3 edge by finding a tri on the disjoint c3, and
        // finding the edge adjacent to it with the same c3 membership.
        for(int j = 0; j < tri_list->size(); j++) {
          if(burn_ent[tri_list->at(j)] and ent_cell[tri_list->at(j)] == id ) {
            apf::Downward d_edge;
            apf::Downward d_vert;
            m->getDownward(tri_list->at(j), 1, d_edge);
            m->getDownward(tri_list->at(j), 0, d_vert);
            int v1 = findIn(d_vert, 3, v_0cell);
            assert(v1 > -1);
            int lu_tri_ed [3][2] = {{2,0},{0,1},{1,2}};
            if(m->toModel(d_edge[lu_tri_ed[v1][0]]) == m->toModel(ent)) {
              c3_edge[id] = d_edge[lu_tri_ed[v1][0]];
            }
            else {
              assert(m->toModel(d_edge[lu_tri_ed[v1][1]]) == m->toModel(ent));
              c3_edge[id] = d_edge[lu_tri_ed[v1][1]];
            }
            j = tri_list->size();
          }
        }

      }
    }
  }

}

void vd_cell_det::load_graph() {
  for(int i = 0; i < 4; i++) {
    ent_sz.at(i) = ent_list.es.at(i).size();
  }

  clear_mesh();
  if(cb_write_flag) {
    cg.load_cb(cb, calc_corner);
    if(cb->get_cell_ext_gmi(0, c0_tag))
      ext_0cell = true;

    if(cb->get_0c_corner_gmi(c0_tag))
      ext_corner = true;
    cells.at(0).reserve(cb->get_sz(1));
    cells.at(1).reserve(cb->get_sz(2));
    cells.at(2).reserve(cb->get_sz(3));

    adj.at(0).reserve(cb->get_sz(2));
    adj.at(1).reserve(cb->get_sz(3));

    burn_new();
    reorder_cells();
  }
  else
    burn_from_cb();

  burn_lower();

  if(cb_write_flag) {
    c3_nbr = cells.at(2).size();
    if(ext_0cell){
      if(calc_corner and ext_corner) {

        std::map<int, int> c3_mark;
        std::map<int, bool> c2_flag;
        std::map<int, int> c2_mark;

        ext_sz = 0;
        int c2_v_sz = 0;
        // For each disjoint set of triangles on the boundary, insert a new ext  
        // 3c. Add the adjancency list for the exterior.
        adj.at(1).reserve(adj.at(1).size()+cells.at(1).size());
        cells.at(2).reserve(cells.at(2).size()+cells.at(1).size());

        for(int i = 0; i < cells.at(1).size(); i++) {
          if(ext_2cell[cells.at(1).at(i)]) {

            ext_sz = ext_sz + 1;
            c3_mark[cells.at(1).at(i)] = -ext_sz;

            cells.at(2).push_back(-ext_sz);
            adj.at(1).resize(adj.at(1).size() + 1);
            adj.at(1).back().reserve(3);
            adj.at(1).back().push_back(cells.at(1).at(i));
          }
        }
        // Add the exterior 1cell. 
        cells.at(0).push_back(-1);
        // For each disj. set of tri. on the bound. insert a new ext 2c for the 
        // bounding 1c. Add the bounding 1c to the adj. list of the new ext 2c.
        adj.at(0).reserve(adj.at(0).size()+ext_sz);
        cells.at(1).reserve(cells.at(1).size()+ext_sz);

        for(int i = 0; i < cells.at(1).size(); i++) {
          if(ext_2cell[cells.at(1).at(i)]) {
            int ext3_curr = c3_mark[cells.at(1).at(i)];
            for (int j = 0; j < adj.at(0).at(i).size(); j++) {
              if(!c2_flag[adj.at(0).at(i).at(j)]) {
                c2_v_sz = c2_v_sz + 1;
                c2_mark[adj.at(0).at(i).at(j)] = -c2_v_sz;
                c2_flag[adj.at(0).at(i).at(j)] = true;
                cells.at(1).push_back(-c2_v_sz);

                adj.at(0).resize(adj.at(0).size() + 1);
                adj.at(0).back().reserve(2);
                adj.at(0).back().push_back(adj.at(0).at(i).at(j));
                adj.at(0).back().push_back(-1);

              }
              adj.at(1).at(c3_nbr - ext3_curr - 1).push_back(
                                                c2_mark[adj.at(0).at(i).at(j)]);
            }
          }
        }
      }

      else {
        ext_sz = 1;
        cells.at(2).push_back(-1);
        adj.at(1).resize(adj.at(1).size() + 1);
        for(int i = 0; i < cells.at(1).size(); i++) {
          if(ext_2cell[cells.at(1).at(i)]) {
            //std::cout << cells.at(1).at(i) << " " << std::endl;
            adj.at(1).back().push_back(cells.at(1).at(i));
          }
        }
        //std::cout << std::endl;
      }
    }
  }

  std::cout << "Cells " << std::endl;
  for(int i = 0; i < cells.size(); i++) {
    for(int j = 0; j < cells.at(i).size(); j++) {
      if(cells.at(i).at(j) < 0)
        std::cout << i+1 << "Ext " << cells.at(i).at(j) << " ";
      else
        std::cout << i+1 << "c" << cells.at(i).at(j) << " ";
    }
    std::cout << std::endl;
  }

  if(cb_write_flag) {
    cg.load_cells(c0_tag, &cells, &adj);
    update_path();
  }
  coll_2c_edges();
  compare_cells();
}

vd_cell_det::~vd_cell_det() {
  clear();

  ent_list.clear();
  vert_c1_map.clear();
  c2_c3real1.clear();
  c2_c3real2.clear();

  delete cb_dj;
}

void vd_cell_det::compare_cells() {
  ent_conn e1;
  ent_conn e2;
  ent_conn e_temp;
  cb->get_conn_dim_gmi(1, 0, c0_tag, &e1);
  cb->get_conn_dim_gmi(2, 0, c0_tag, &e2);

  ent_conn e1_dj;
  ent_conn e2_dj;
  cb_dj->get_conn_dim(1, 0, 0, &e1_dj);
  cb_dj->get_conn_dim(2, 0, 0, &e2_dj);

  std::vector<int>::iterator it_1c;

  std::map<int, int> comp_nbr1{};
  std::map<int, int> comp_nbr2{};

  // The count the number of disjoint sets belonging to the same 1cell.
  // Assert that the disjoint set also belongs to a cell bounded by the 0cell.
  for(int i = 0; i < e1_dj.conn.size(); i++) {
    int c1_curr = m->getModelTag(get_mdl(1, e1_dj.conn.at(i)));
    it_1c = std::find(e1.conn.begin(), e1.conn.end(), c1_curr);
    assert(it_1c != e1.conn.end());
    comp_nbr1[c1_curr] = comp_nbr1[c1_curr] + 1;

    // The cell of the disjoint set of entities is bounded by the 0cell.
    cb->get_conn_gmi(1, c1_curr, &e_temp);
    it_1c = std::find(e_temp.conn.begin(), e_temp.conn.end(), c0_tag);
    assert(it_1c != e_temp.conn.end());
  }

  // The count the number of disjoint sets belonging to the same 2cell.
  // Assert that the disjoint set also belongs to a cell bounded by the 0cell.
  for(int i = 0; i < e2_dj.conn.size(); i++) {
    int c2_curr = m->getModelTag(get_mdl(2, e2_dj.conn.at(i)));
    it_1c = std::find(e2.conn.begin(), e2.conn.end(), c2_curr);
    assert(it_1c != e2.conn.end());
    comp_nbr2[c2_curr] = comp_nbr2[c2_curr] + 1;
  }

  // For 2cells belonging to disjoint entities, get the bounding disjoint 1cells.
  // If the actual 1cell has two disjoint parts, assert either the 2cell has 
  // multiple disjoint components or is only bounded by the 1cell.
  for(int i = 0; i < e2_dj.conn.size(); i++) {
    cb_dj->get_conn(2, e2_dj.conn.at(i), &e_temp);
    int c2_curr = m->getModelTag(get_mdl(2, e2_dj.conn.at(i)));
    for(int j = 0; j < e_temp.conn.size(); j++) {
      int c1_curr = m->getModelTag(get_mdl(1, e_temp.conn.at(j)));
      //if(comp_nbr1[c1_curr] > 1) {
      // TODO the following is wrong. A disjoint 2cell is either is in a single 
      // piece and connected to two copies of the same 1cell or is in n 
      // pieces and connected to n-1 two copies of the same 1cell.
        //assert(comp_nbr2[c2_curr] > 1 or e_temp.conn.size() == 1);
      //}
    }
  }
}

int vd_cell_det::get_3c_nbr() {
  return c3_nbr;
}

int vd_cell_det::get_3c_nbr_ext() {
  if(calc_corner and ext_corner) {
    return c3_nbr - 1;
  }
  else
    return c3_nbr;
}

int vd_cell_det::get_2c_nbr() {
  return cells.at(1).size();
}

// Get the list of same dim entities belonging to the given disjoint cell.
std::vector<apf::MeshEntity*> vd_cell_det::get_dim_ent(int dim, int id) {
  assert(dim < 4 and dim > 0);
  std::vector<apf::MeshEntity*> ents(0);
  if(ent_list.es.size() > 0) {
    ents.reserve(ent_list.es.at(dim).size());
    for(int i = 0; i < ent_list.es.at(dim).size(); i++) {
      apf::MeshEntity* ent = ent_list.es.at(dim).at(i);
      apf::ModelEntity* mdl = m->toModel(ent);

      if(m->getModelType(mdl) == dim and ent_cell[ent] == id)
        ents.push_back(ent);
    }
  }
  return ents;
}

void vd_cell_det::coll_2c_edges() {

  apf::Vector3 temp;
  apf::Vector3 v_pos;

  m->getPoint(v_0cell, 0, v_pos);

  e_2c.clear();
  e_2c.resize(cells.at(1).size());

  std::vector<int> c3_adj_1(0);
  std::vector<int> c3_adj_2(0);
  std::vector<double> c3_ang(0);

  c3_adj_1.resize(cells.at(1).size());
  c3_adj_2.resize(cells.at(1).size());
  c3_ang.resize(cells.at(1).size());

  for(int i = 0; i < cells.at(1).size(); i++) {
    c3_adj_1.at(i) = -1;
    c3_adj_2.at(i) = -1;
    c3_ang.at(i) = 0;
  }

  // TODO pick the edge closest to both bounded 3cell edges.
  for(int i = 0; i < cells.at(2).size(); i++) {
    if(cells.at(2).at(i) > -1) {
      for(int j = 0; j < adj.at(1).at(i).size(); j++) {
        if(c3_adj_1.at(adj.at(1).at(i).at(j)) == -1)
          c3_adj_1.at(adj.at(1).at(i).at(j)) = i;
        else {
          assert(c3_adj_2.at(adj.at(1).at(i).at(j)) == -1);
          c3_adj_2.at(adj.at(1).at(i).at(j)) = i;
        }
      }
    }
  }

  // TODO pick the edge closest to both bounded 3cell edges.
  for(int i = 0; i < ent_list.es.at(1).size(); i++) {
    apf::MeshEntity* ent = ent_list.es.at(1).at(i);
    mdl_curr = m->toModel(ent);
    int c2_id = ent_cell[ent];
    //std::cout << ent << " " << m->getModelType(mdl_curr) << "c"
    //          << m->getModelTag(mdl_curr) << " " << c2_id
    //          << std::endl;
    if (m->getModelType(mdl_curr) == 2) {
      assert(c3_adj_1.at(c2_id) != -1);
      double ang_curr = std::fabs(vd_int_angle(m, c3_edge[c3_adj_1.at(c2_id)],
                                                    ent, v_0cell));

      if(c3_adj_2.at(ent_cell[ent]) != -1) {
        ang_curr = ang_curr + 
                    std::fabs(vd_int_angle(m, c3_edge[c3_adj_2.at(c2_id)],
                                                    ent, v_0cell));
      }
      if(std::fabs(c3_ang.at(c2_id)) < std::numeric_limits<double>::epsilon() 
         or ang_curr < c3_ang.at(c2_id) ) {
        c3_ang.at(c2_id) = ang_curr;
        e_2c.at(ent_cell[ent]) = ent;
      }
    }
  }

  /*
        apf::MeshEntity* v_curr = getEdgeVertOppositeVert(m, ent, v_0cell);
        m->getPoint(v_curr, 0, temp);
        temp = norm_0(temp - v_pos);
        double ang_curr = temp*c2_dir[ent_cell[ent]];
        if(!c2_flag[ent_cell[ent]]) {
          e_2c.at(ent_cell[ent]) = ent;
          c2_flag[ent_cell[ent]] = true;
          c2_ang[ent_cell[ent]] = ang_curr;
        }
        else {
          if(ang_curr < c2_ang[ent_cell[ent]]) {
            e_2c.at(ent_cell[ent]) = ent;
            c2_ang[ent_cell[ent]] = ang_curr;
          }
        }
  */

//check the 3cell connections, have the 3cell interior directions
//get the smallest internal angle edges
}

std::vector<apf::MeshEntity*> vd_cell_det::get_c2_edges() {
  return e_2c;
}

std::vector<apf::MeshEntity*> vd_cell_det::get_c3_edges() {
  std::vector<apf::MeshEntity*> temp(c3_nbr);
  for(int i = 0; i < c3_nbr; i++)
    temp.at(i) = c3_edge[i];
  return temp;
}

void vd_cell_det::set_c2_edge(int ent_cell_curr, apf::MeshEntity* e_curr) {
  e_2c.at(ent_cell_curr) = e_curr;
}

void vd_cell_det::get_path_ngon(int c3_cp, ngon_gmi* ng) {
  cg.get_path_ngon(c3_cp, ng);
  //shift_ngon_gmi(ng);
}

void vd_cell_det::conv_path_ngon_gmi(ngon_gmi* ng) {
  std::cout << "Converting ngon paths" << std::endl;
  int cells3_nbr;
  if(ext_0cell)
    cells3_nbr = cells.at(2).size()-1;
  else
    cells3_nbr = cells.at(2).size();

  //std::cout << "cells3_nbr " << cells3_nbr << std::endl;
  for(int i = 0; i < ng->cells.size(); i++) {
    for(int j = 0; j < ng->cells.at(i).second.size(); j++) {
      std::cout << "Local 2c" << ng->cells.at(i).second.at(j) << " ";
      //apf::ModelEntity* mdl = burn_mdl.at(1)[ng->cells.at(i).second.at(j)-cells3_nbr];
      apf::ModelEntity* mdl = burn_mdl.at(1)[ng->cells.at(i).second.at(j)];
      ng->cells.at(i).second.at(j) = m->getModelTag(mdl);
      //ng->cells.at(i).second.at(j) = ng->cells.at(i).second.at(j)-cells3_nbr;
      std::cout << "2c" << ng->cells.at(i).second.at(j) << std::endl;
      //ng->cells.at(i).second.at(j) = ng->cells.at(i).second.at(j);
      //std::cout << " new 2c" << ng->cells.at(i).second.at(j) << std::endl;
    }
  }
  if(ng->cp.first != -1) {
    apf::ModelEntity* mdl = burn_mdl.at(2)[ng->cp.first];
    ng->cp.first = m->getModelTag(mdl);
  }
  if(ng->cp.second != -1) {
    apf::ModelEntity* mdl = burn_mdl.at(2)[ng->cp.second];
    ng->cp.second = m->getModelTag(mdl);
  }
}
int vd_cell_det::conv_path_2c_gmi(int n_id) {
  return m->getModelTag(burn_mdl.at(1)[n_id]);
}
int vd_cell_det::conv_path_3c_gmi(int n_id) {
  if(n_id != -1) {
    return m->getModelTag(burn_mdl.at(2)[n_id]);
  }
  return -1;
}

void vd_cell_det::shift_ngon_gmi(ngon_gmi* ng) {

  std::cout << "Shifting ngon paths" << std::endl;
  int cells3_nbr;
  if(ext_0cell)
    cells3_nbr = cells.at(2).size()-1;
  else
    cells3_nbr = cells.at(2).size();

  std::cout << "cells3_nbr: " << cells3_nbr << std::endl;

  for(int i = 0; i < ng->cells.size(); i++) {
    for(int j = 0; j < ng->cells.at(i).first.size(); j++) {
      std::cout << "3c" << ng->cells.at(i).first.at(j) << std::endl;
      apf::ModelEntity* mdl = burn_mdl.at(2)[ng->cells.at(i).first.at(j)];
      ng->cells.at(i).first.at(j) = m->getModelTag(mdl);
      std::cout << " new 3c" << ng->cells.at(i).first.at(j) << std::endl;
    }
    for(int j = 0; j < ng->cells.at(i).second.size(); j++) {
      //apf::ModelEntity* mdl = burn_mdl.at(1)[ng->cells.at(i).second.at(j)-cells3_nbr];
      //ng->cells.at(i).second.at(j) = m->getModelTag(mdl);
      //ng->cells.at(i).second.at(j) = ng->cells.at(i).second.at(j)-cells3_nbr;
      std::cout << "2c" << ng->cells.at(i).second.at(j) << std::endl;
      //ng->cells.at(i).second.at(j) = ng->cells.at(i).second.at(j);
      //std::cout << " new 2c" << ng->cells.at(i).second.at(j) << std::endl;
    }
  }
}

int vd_cell_det::get_circ_type(int c_in) {
  return cg.get_circ_type(c_in);
}

bool vd_cell_det::get_ext() {
  return ext_0cell;
}

void vd_cell_det::set_calc_ext(bool fix) {
  calc_ext = fix;
}

void vd_cell_det::set_calc_corner(bool fix) {
  calc_corner = fix;
  cg.set_calc_corner(fix);
}

void vd_cell_det::set_proj(PROJ_TYPE PROJ) {
  int proj_flag = (int)PROJ;
  assert(proj_flag < (int)PROJ_TYPE::END);

  if(proj_flag > (int)PROJ_TYPE::FIXED)
    calc_ext = true;
  else
    calc_ext = false;

  if(proj_flag > (int)PROJ_TYPE::FIXED and 
     proj_flag < (int)PROJ_TYPE::END) {
    calc_corner = true;
    cg.set_calc_corner(calc_corner);
  }
  else {
    calc_corner = false;
    cg.set_calc_corner(false);
  }
}

int vd_cell_det::get_circ_sz() {
  return cg.get_circuits()->size();
}
int vd_cell_det::get_path_sz() {
  return cg.get_paths()->size();
}

circuit * vd_cell_det::get_circ(int t_curr) {
  return &circ_tup.at(t_curr);
}

// For corner vertices, remove the virtual exterior cells before collecting mesh
// entities on the paths and circuits.
void vd_cell_det::treat_cells() {
  for(int i = 1; i < cells.size(); i++) {
    for(int j = 0; j < ext_sz; j++) {
      cells.at(i).pop_back();

      adj.at(i-1).back().clear();
      adj.at(i-1).pop_back();
    }
  }
  cells.at(0).pop_back();

  // Add a single exterior.
  cells.at(2).push_back(-1);
  adj.at(1).resize(adj.at(1).size() + 1);
  adj.at(1).back().reserve(cells.at(1).size());
  for(int i = 0; i < cells.at(1).size(); i++) {
    if(ext_2cell[cells.at(1).at(i)]) {
      //std::cout << cells.at(1).at(i) << " " << std::endl;
      adj.at(1).back().push_back(cells.at(1).at(i));
    }
  }

}

void vd_cell_det::update_path() {
  std::cout << "Updating 0cell " 
            << c0_tag << std::endl; 

  clear_path();

  // Graph id.
  std::vector< crc > * circ = cg.get_circuits();
  circ_tup.resize(circ->size());

  if(calc_corner and ext_corner)
    treat_cells();

  int cells3_nbr;
  if(ext_0cell)
    cells3_nbr = cells.at(2).size()-1;
  else
    cells3_nbr = cells.at(2).size();

  for (int i = 0; i < circ->size(); i++) {
    if(circ->at(i).first.size() > 0) {
  //std::cout << "Circuit " << i << std::endl;
      // The circuit entities, to be taken from vd_Edge objects.
      std::map<int, int> map_1;
      std::map<int, int> map_2;

      circ_tup.at(i).first.first.reserve(circ->at(i).first.size());
      circ_tup.at(i).first.second.reserve(circ->at(i).first.size());

      for (int j = 0; j < circ->at(i).first.size(); j++) {
        int U = circ->at(i).first.at(j).first.source();
        int V = circ->at(i).first.at(j).first.target();
        map_1[U] = 0;
        map_1[V] = 0;
        map_2[U] = 0;
        map_2[V] = 0;
      }

      // Iterate over the edges.
      for (int j = 0; j < circ->at(i).first.size(); j++) {
        int U = circ->at(i).first.at(j).first.source();
        int V = circ->at(i).first.at(j).first.target();
        //std::cout << "U " << U << " V " << V << std::endl;
        //if (U == -1)
        //  std::cout << "3c" << "Ext ";
        //else
        //  std::cout << "3c" <<  cells.at(2).at(U);

        //std::cout << " 2c" <<  cells.at(1).at(V-cells3_nbr) 
        //          << std::endl;

        //std::cout << "map_1[U] " << map_1[U] << " map_1[V] " << map_1[V] 
        //          << " map_2[U] " << map_2[U] << " map_2[V] " << map_2[V] 
        //          << std::endl;

        if(map_1[U] == 0) {
          map_1[U] = V;
        }
        // There can be at most two connections. Both connections are filled 
        // after visiting both neighboring cells.
        else {
          assert(map_2[U] == 0);
          map_2[U] = V;
        }
        if(map_1[V] == 0) {
          if(U == -1)
            map_1[V] = -1;
          else
            map_1[V] = U;
        }
        else {
          assert(map_2[V] == 0);
          if(U == -1)
            map_2[V] = -1;
          else
            map_2[V] = U;
        }
      }

      int last = circ->at(i).first.at(0).first.target();
      int next;
      if(map_1[last] == -1)
        next = -1;
      else
        next = map_1[last];

      assert(last >= cells3_nbr);
      circ_tup.at(i).first.second.push_back(last-cells3_nbr);
      //std::cout << "First " << last << std::endl;
      // Assume circuit is closed.
      while(next != circ->at(i).first.at(0).first.target()) {
        //std::cout << "Last " << last << " Next " << next << std::endl;
        if(next < cells3_nbr) {
          circ_tup.at(i).first.first.push_back(next);
        }
        else {
          circ_tup.at(i).first.second.push_back(next-cells3_nbr);
        }
        //std::cout << "map_1[" << next << "] " << map_1[next] 
        //          << "map_2[" << next << "] " << map_2[next] 
        //          << std::endl;

        int next_curr;
        if(map_1[next] == -1)
          next_curr = -1;
        else
          next_curr = map_1[next];

        if(last == next_curr) {
          last = next;
          if(map_2[next] == -1)
            next = -1;
          else
            next = map_2[next];
        }
        else {
          last = next;
          if(map_1[next] == -1)
            next = -1;
          else
            next = map_1[next];
        }
      }

      // The disjoint graph entities, to be copied.
      circ_tup.at(i).second = circ->at(i).second;
    }
  }

  std::vector< cell_path >* paths = cg.get_paths();
  //path[cell0];

  path.clear();
  path.resize(paths->size());

  //std::cout << "Number of couples after: " << paths->size() << std::endl;
  for (int i = 0; i < paths->size(); i++) {
    path.at(i).first = paths->at(i).first;
    path.at(i).second.resize(paths->at(i).second.size());
    for (int j = 0; j < paths->at(i).second.size(); j++) {
      path.at(i).second.at(j) = paths->at(i).second.at(j);
    }
  }

/*
  for(int dim = 3; dim > 0; dim--) {
    for(int i = 0; i < cells.at(dim-1).size(); i++) {
      int id_curr = cells.at(dim-1).at(i);
      if(id_curr == -1) {
        std::cout << dim << "c" << id_curr << " Ext";
      }
      else {
        std::cout << dim << "c" << id_curr
                  << " " << m->getModelTag(burn_mdl.at(dim-1)[id_curr]);
      }
      std::cout << std::endl;
    }
  }
*/
}

// Return bounding disjoint strata of given stratum.
std::vector<int> vd_cell_det::get_bound_ids(int dim, int id) {
  return adj.at(dim-2).at(id);
}

void vd_cell_det::set_ng_couple(int p_id) {
  cg.get_path_ngon(p_id, &ng_local);
}

void vd_cell_det::find_slice_ng(int ng_id, 
            std::vector<std::pair<std::vector<int >, std::vector<int > > >*
                                                              path_cells,
            std::vector< std::pair< std::pair<int,int>, 
                         std::vector<std::vector<int > > > >* slice_cells) {

  for(int i = 0; i < slice_cells->size(); i++) {
    for(int j = 0; j < slice_cells->at(i).second.size(); j++) {
      slice_cells->at(i).second.at(j).clear();
      slice_cells->at(i).second.at(j).clear();
      slice_cells->at(i).second.at(j).clear();
    }
    slice_cells->at(i).second.clear();
  }
  slice_cells->clear();

  //ng.print();
  for(int i = 0; i < path_cells->size(); i++) {
    path_cells->at(i).first.clear();
    path_cells->at(i).second.clear();
  }

  path_cells->resize(ng_local.ngons.at(ng_id).size());
  for(int i = 0; i < ng_local.ngons.at(ng_id).size(); i++ ) {
    int path_id = ng_local.ngons.at(ng_id).at(i);

    // 3c
    path_cells->at(i).first.reserve(ng_local.cells.at(path_id).first.size());
    for(int j = 0; j < ng_local.cells.at(path_id).first.size(); j++) {
      path_cells->at(i).first.push_back(ng_local.cells.at(path_id).first.at(j));
    }

    // 2c
    path_cells->at(i).second.reserve(ng_local.cells.at(path_id).second.size());
    for(int j = 0; j < ng_local.cells.at(path_id).second.size(); j++) {
      path_cells->at(i).second.push_back(ng_local.cells.at(path_id).second.at(j));
    }
  }
  cg.find_slice(path_cells, slice_cells);
}

// Return the indices of real local ids 1-, 2- and, 3-cells.
void vd_cell_det::find_slice_circ(int circ_id, 
            std::vector<std::pair<std::vector<int >, std::vector<int > > >*
                                                              path_cells,
            std::vector< std::pair< std::pair<int,int>, 
                         std::vector<std::vector<int > > > >* slice_cells) {

  for(int i = 0; i < slice_cells->size(); i++) {
    for(int j = 0; j < slice_cells->at(i).second.size(); j++) {
      slice_cells->at(i).second.at(j).clear();
    }
    slice_cells->at(i).second.clear();
  }
  slice_cells->clear();
  // path_cells
  for(int i = 0; i < path_cells->size(); i++) {
    path_cells->at(i).first.clear();
    path_cells->at(i).second.clear();
  }
  path_cells->resize(1);

  // 3c
  path_cells->at(0).first.reserve(circ_tup.at(circ_id).first.first.size());
  for(int j = 0; j < circ_tup.at(circ_id).first.first.size(); j++) {
    path_cells->at(0).first.push_back(circ_tup.at(circ_id).first.first.at(j));
  }

  // 2c
  path_cells->at(0).second.reserve(circ_tup.at(circ_id).first.second.size());
  for(int j = 0; j < circ_tup.at(circ_id).first.second.size(); j++) {
    path_cells->at(0).second.push_back(circ_tup.at(circ_id).first.second.at(j));
  }

  // 2 slices
  slice_cells->resize(2);
  //1-, 2-, 3-cells
  slice_cells->at(0).second.resize(3);
  slice_cells->at(1).second.resize(3);

  // 3-cells
  slice_cells->at(0).second.at(2).reserve(circ_tup.at(circ_id).first.first.size());
  // 2-cells
  slice_cells->at(0).second.at(1).reserve(circ_tup.at(circ_id).first.second.size());
  // 1-cells. Add all bounding 1-cells of the 2-cells on the slice, sort and find 
  // unique ones.
  slice_cells->at(0).second.at(0).reserve(circ_tup.at(circ_id).first.second.size() * 2);

  // First disjoint graph:
  std::vector<int>::iterator it_2c;
  for(int i = 0; i < circ_tup.at(circ_id).second.first.size(); i++) {
    int c_id = circ_tup.at(circ_id).second.first.at(i);
    if(c_id < c3_nbr and c_id > - ext_sz - 1) {
      slice_cells->at(0).second.at(2).push_back(c_id);
    }
    else {
      c_id = c_id - c3_nbr;
      slice_cells->at(0).second.at(1).push_back(c_id);
      it_2c = std::find(cells.at(1).begin(), cells.at(1).end(), c_id);
      assert(it_2c != cells.at(1).end());
      c_id = std::distance(cells.at(1).begin(), it_2c);
      for(int j = 0; j < adj.at(0).at(c_id).size(); j++) {
        slice_cells->at(0).second.at(0).push_back(adj.at(0).at(c_id).at(j));
      }
    }
  }
  // Get the unique disjoint 1cells:
  std::sort(slice_cells->at(0).second.at(0).begin(), 
            slice_cells->at(0).second.at(0).end());

  it_2c = std::unique(slice_cells->at(0).second.at(0).begin(), 
            slice_cells->at(0).second.at(0).end());
  int un_ln = std::distance(slice_cells->at(0).second.at(0).begin(), it_2c);
  slice_cells->at(0).second.at(0).resize(un_ln);

  // Second disjoint graph
  for(int i = 0; i < circ_tup.at(circ_id).second.second.size(); i++) {
    int c_id = circ_tup.at(circ_id).second.second.at(i);
    if(c_id < c3_nbr and c_id > -ext_sz - 1) {
      slice_cells->at(1).second.at(2).push_back(c_id);
    }
    else {
      c_id = c_id - c3_nbr;
      slice_cells->at(1).second.at(1).push_back(c_id);
      it_2c = std::find(cells.at(1).begin(), cells.at(1).end(), c_id);
      assert(it_2c != cells.at(1).end());
      c_id = std::distance(cells.at(1).begin(), it_2c);
      for(int j = 0; j < adj.at(0).at(c_id).size(); j++) {
        slice_cells->at(1).second.at(0).push_back(adj.at(0).at(c_id).at(j));
      }
    }
  }
  // Get the unique disjoint 1cells:
  std::sort(slice_cells->at(1).second.at(0).begin(), 
            slice_cells->at(1).second.at(0).end());

  it_2c = std::unique(slice_cells->at(1).second.at(0).begin(), 
            slice_cells->at(1).second.at(0).end());
  un_ln = std::distance(slice_cells->at(1).second.at(0).begin(), it_2c);
  slice_cells->at(1).second.at(0).resize(un_ln);

  for(int i = 0; i < slice_cells->at(0).second.at(0).size(); i++) {
    for(int j = 0; j < slice_cells->at(1).second.at(0).size(); j++) {
      assert(slice_cells->at(0).second.at(0).at(i) !=
              slice_cells->at(1).second.at(0).at(j));
    }
  }

}


// From local id, return the model object.
apf::ModelEntity* vd_cell_det::get_mdl(int dim_in, int id_in) {
  assert(id_in > -1 and id_in < cells.at(dim_in-1).size());
  return burn_mdl.at(dim_in-1)[id_in];
}

void vd_cell_det::coll_c2_dir() {
  apf::Vector3 v_pos;
  m->getPoint(v_0cell, 0, v_pos);

  std::vector<std::vector<apf::MeshEntity*> > tri
      (0, std::vector<apf::MeshEntity*>(0) );
  tri.resize(cells.at(1).size());
  for(int i = 0; i < tri.size(); i++) {
    tri.at(i).reserve(ent_list.es.at(2).size());
  }
  for(int i = 0; i < ent_list.es.at(2).size(); i++) {
    apf::MeshEntity* e_curr = ent_list.es.at(2).at(i);
    apf::ModelEntity* mdl_next = m->toModel(e_curr);
    if(m->getModelType(mdl_next) == 2) {
      tri.at(ent_cell[e_curr]).push_back(e_curr);
    }
  }

  for(int i = 0; i < tri.size(); i++) {
    c2_dir[i] = vd_find_dist_w_tri(m, &tri.at(i), v_0cell);
    c2_dir[i] = norm_0(c2_dir[i] - v_pos);
  }

  for(int i = 0; i < tri.size(); i++) {
    tri.at(i).clear();
  }
  tri.clear();
}

void vd_cell_det::coll_c3_dir() {
  apf::Vector3 v_pos;
  m->getPoint(v_0cell, 0, v_pos);

  std::vector<std::vector<apf::MeshEntity*> > tet
      (0, std::vector<apf::MeshEntity*>(0) );

  tet.resize(cells.at(2).size());
  for(int i = 0; i < tet.size(); i++) {
    tet.at(i).reserve(ent_list.es.at(3).size());
  }
  for(int i = 0; i < ent_list.es.at(3).size(); i++) {
    apf::MeshEntity* e_curr = ent_list.es.at(3).at(i);
    apf::ModelEntity* mdl_next = m->toModel(e_curr);
    if(m->getModelType(mdl_next) == 3) {
      tet.at(ent_cell[e_curr]).push_back(e_curr);
    }
  }

  for(int i = 0; i < tet.size(); i++) {
    c3_dir[i] = vd_find_dist_w(m, &tet.at(i), v_0cell);
    c3_dir[i] = norm_0(c3_dir[i] - v_pos);
  }

  if(ext_0cell) {
    int tet_sz = 0;
    for(int i = 0; i < tet.size(); i++)
      tet_sz = tet_sz + tet.at(i).size();
    tet.at(0).reserve(tet_sz);
    for(int i = 1; i < tet.size(); i++)
      tet.at(0).insert(tet.at(0).end(), tet.at(i).begin(), tet.at(i).end());
    c3_dir[-1] = vd_find_dist_w(m, &tet.at(0), v_0cell);
    c3_dir[-1] = norm_0(c3_dir[-1] - v_pos)*(-1);
  }

  for(int i = 0; i < tet.size(); i++) {
    tet.at(i).clear();
  }
  tet.clear();
}

apf::Vector3 vd_cell_det::get_c2_dir(int c_id) {
  return c2_dir[c_id];
}

apf::Vector3 vd_cell_det::get_c3_dir(int c_id) {
  return c3_dir[c_id];
}

void vd_cell_det::set_vert_c1_map(std::map<apf::MeshEntity*, int> map_in) {
  vert_c1_map = map_in;
}

std::map<apf::MeshEntity*, int> vd_cell_det::get_vert_c1_map() {
  return vert_c1_map;
}

void vd_cell_det::set_cb_flag(bool on_off) {
  cb_write_flag = on_off;
}

void vd_cell_det::burn_cell(int c_dim, int id, apf::MeshEntity* ent) {

  int type = m->getType(ent);
  int d = apf::Mesh::typeDimension[type];
  assert(d == c_dim);

  // Add the current disjoint cell segment to the list.

  if(cb_write_flag) {
    cells.at(c_dim-1).push_back(id);
  }

  // Associate the current disjoint cell segment with the current model.
  mdl_curr = m->toModel(ent);
  burn_mdl.at(c_dim-1)[id] = mdl_curr;

  ent_cell[ent] = id;
  burn_ent[ent] = true;

  if(c_dim == 1) {
    ent_ofire.at(0).clear();
    ent_ofire.at(1).clear();
    return;
  }

  if(cb_write_flag) {
    if(c_dim > 1) {
      adj.at(c_dim-2).resize(adj.at(c_dim-2).size() + 1);
      adj.at(c_dim-2).back().reserve(cb->get_sz(c_dim-1));
    }

    if(c_dim == 2 and cb->get_cell_ext_gmi(2, m->getModelTag(mdl_curr)))
      ext_2cell[id] = true;
  }

  apf::ModelEntity* mdl_next;
  apf::MeshEntity* ent_curr;
  apf::MeshEntity* ent_next;

  apf::Downward down;
  apf::Downward vert;
  apf::Up up;

  std::vector<int>::iterator it_l;


  // Burn the first entity of the cell dimension.
  ent_ofire.at(0).clear();
  ent_ofire.at(1).clear();

  assert(ent_sz.at(c_dim) >= 0);
  assert(ent_sz.at(c_dim-1) >= 0);
  ent_ofire.at(0).reserve(ent_sz.at(c_dim));
  ent_ofire.at(1).reserve(ent_sz.at(c_dim-1));

  ent_ofire.at(0).push_back(ent);

  bool all_burned = false;
  while(!all_burned) {
    all_burned = true;
    while(ent_ofire.at(0).size() > 0) {
      ent_curr = ent_ofire.at(0).back();
      //std::cout << ent_curr << " " 
      //          << vd_get_pos(m, ent_curr) - vd_get_pos(m, v_0cell)
      //          << std::endl;

      int e_sz = m->getDownward(ent_curr, c_dim-1, down);
      for(int i = 0; i < e_sz; i++) {
        ent_next = down[i];

        int e_sz2 = m->getDownward(ent_next, 0, vert);
        if(findIn(vert, e_sz2, v_0cell) > -1) {
          mdl_next = m->toModel(ent_next);
          int c_next = m->getModelType(mdl_next);

          if(burn_ent[ent_next]) {
            if(cb_write_flag) {
              int i_next = ent_cell[ent_next];
              if(c_next == c_dim - 1) {
                it_l = std::find(adj.at(c_dim-2).back().begin(), 
                            adj.at(c_dim-2).back().end(), i_next);
                if(it_l == adj.at(c_dim-2).back().end())
                  adj.at(c_dim-2).back().push_back(i_next);
              }
            }
          }
          else {
            if(c_next == c_dim) {
              assert(mdl_next == mdl_curr);

              ent_ofire.at(1).push_back(ent_next);
              burn_ent[ent_next] = true;
              ent_cell[ent_next] = id;

            }
            all_burned = false;
          }
        }
      }

      ent_ofire.at(0).pop_back();
    }

    while(ent_ofire.at(1).size() > 0) {
      ent_curr = ent_ofire.at(1).back();

      m->getUp(ent_curr, up);
      for(int i = 0; i < up.n; i++) {
        ent_next = up.e[i];

        int e_sz2 = m->getDownward(ent_next, 0, vert);
        // TODO all should be finding
        if(findIn(vert, e_sz2, v_0cell) > -1) {
          mdl_next = m->toModel(ent_next);
          int c_next = m->getModelType(mdl_next);

          if(burn_ent[ent_next]) {
            int i_next = ent_cell[ent_next];
            assert(c_next != c_dim - 1);
/*
            if(c_next == c_dim - 1) {
              it_l = std::find(adj.at(c_dim-2).back().begin(), 
                          adj.at(c_dim-2).back().end(), i_next);
              if(it_l == adj.at(c_dim-2).back().end())
                adj.at(c_dim-2).back().push_back(i_next);
            }
*/
          }
          else {
            // TODO is it not possible to have entities belonging to other 
            // cells?
            if(c_next == c_dim) {
              assert(mdl_next == mdl_curr);

              ent_ofire.at(0).push_back(ent_next);
              burn_ent[ent_next] = true;
              ent_cell[ent_next] = id;

              all_burned = false;
            }
          }
        }
      }

      ent_ofire.at(1).pop_back();
    }
  }
  if(cb_write_flag)
    std::sort(adj.at(c_dim-2).back().begin(), adj.at(c_dim-2).back().end());
}

// Burn the disjoint cell to find the bounding disjoint cells.
std::vector<int> vd_cell_det::burn_bound(int c_dim, apf::MeshEntity* ent) {

  std::vector<int> bound(0);

  int type = m->getType(ent);
  int d = apf::Mesh::typeDimension[type];
  assert(d == c_dim);

  // Associate the current disjoint cell segment with the current model.
  mdl_curr = m->toModel(ent);

  std::map<apf::MeshEntity*, bool> burn_ent_bound;

  burn_ent_bound[ent] = true;

  if(c_dim == 1) {
    ent_ofire.at(0).clear();
    ent_ofire.at(1).clear();
    bound.push_back(0);
    return bound;
  }

  apf::ModelEntity* mdl_next;
  apf::MeshEntity* ent_curr;
  apf::MeshEntity* ent_next;

  apf::Downward down;
  apf::Downward vert;
  apf::Up up;

  std::vector<int>::iterator it_l;

  // Burn the first entity of the cell dimension.
  ent_ofire.at(0).clear();
  ent_ofire.at(1).clear();

  assert(ent_sz.at(c_dim) >= 0);
  assert(ent_sz.at(c_dim-1) >= 0);
  ent_ofire.at(0).reserve(ent_sz.at(c_dim));
  ent_ofire.at(1).reserve(ent_sz.at(c_dim-1));
  bound.reserve(ent_sz.at(c_dim-1));

  ent_ofire.at(0).push_back(ent);

  bool all_burned = false;
  while(!all_burned) {
    //get_3c_nbr();
    all_burned = true;
    while(ent_ofire.at(0).size() > 0) {
      ent_curr = ent_ofire.at(0).back();

      int e_sz = m->getDownward(ent_curr, c_dim-1, down);
      for(int i = 0; i < e_sz; i++) {
        ent_next = down[i];

        int e_sz2 = m->getDownward(ent_next, 0, vert);
        if(findIn(vert, e_sz2, v_0cell) > -1) {
          mdl_next = m->toModel(ent_next);
          int c_next = m->getModelType(mdl_next);

          if(burn_ent[ent_next]) {
            int i_next = ent_cell[ent_next];
            if(c_next == c_dim - 1) {
              it_l = std::find(bound.begin(), bound.end(), i_next);
              if(it_l == bound.end())
                bound.push_back(i_next);
            }
          }
          else if(burn_ent_bound[ent_next]) {
          }
          else {
            if(c_next == c_dim) {
              assert(mdl_next == mdl_curr);

              ent_ofire.at(1).push_back(ent_next);
              burn_ent_bound[ent_next] = true;
            }
            all_burned = false;
          }
        }
      }

      //get_3c_nbr();
      ent_ofire.at(0).pop_back();
    }

    while(ent_ofire.at(1).size() > 0) {
      ent_curr = ent_ofire.at(1).back();

      m->getUp(ent_curr, up);
      for(int i = 0; i < up.n; i++) {
        ent_next = up.e[i];

        int e_sz2 = m->getDownward(ent_next, 0, vert);
        if(findIn(vert, e_sz2, v_0cell) > -1) {
          mdl_next = m->toModel(ent_next);
          int c_next = m->getModelType(mdl_next);

          if(burn_ent_bound[ent_next]) {
            assert(c_next != c_dim - 1);
          }
          else {
            if(c_next == c_dim) {
              assert(mdl_next == mdl_curr);

              ent_ofire.at(0).push_back(ent_next);
              burn_ent_bound[ent_next] = true;

              all_burned = false;
            }
          }
        }
      }
      //get_3c_nbr();
      ent_ofire.at(1).pop_back();
    }
  }
  //get_3c_nbr();

  std::sort(bound.begin(), bound.end());
  return bound;
}

void vd_cell_det::get_circ_topo_gmi(std::vector<int>* c2_list, std::vector<int>* c3_list, circuit* circ) {

  c3_list->resize(circ->first.first.size());
  c2_list->resize(circ->first.second.size());

  //std::cout << "Size of circ_in: "<< circ_in->first.size() 
  //    << ", size of cell3: " << cells3[cell_id].conn.size() << std::endl;
  // For the circuit, simply convert the local indices to the global indices.
  //std::cout << "Circuit: " << std::endl;
  for(int i = 0; i < circ->first.first.size(); i++) {

    if (circ->first.first.at(i) != -1) {
      apf::ModelEntity* mdl = burn_mdl.at(2)[circ->first.first.at(i)];
      //std::cout << circ_in->first.first[i] << " ";
      int c3_cur = m->getModelTag(mdl);
      //std::cout << "3c " << c3_cur + 1 << std::endl;
      c3_list->at(i) = c3_cur;
    }
    else {
      c3_list->at(i) = -1;
    }
  }

  for(int i = 0; i < circ->first.second.size(); i++) {
    apf::ModelEntity* mdl = burn_mdl.at(1)[circ->first.second.at(i)];
    int c2_cur = m->getModelTag(mdl);
    //std::cout << "3c " << c3_cur + 1 << std::endl;
    c2_list->at(i) = c2_cur;
    //std::cout << circ_in->first.first[i] << " ";
    //c2_list->at(i) = circ_in->at(circ).first.second.at(i);
  }

}

void vd_cell_det::get_circ_topo_dis_gmi(std::vector<ent_conn>* cs, 
                                                circuit* circ, int g_id) {

  for(int i = 0; i < cs->size(); i++)
    cs->at(i).conn.clear();
  cs->clear();

  cs->resize(3);
  cs->at(0).conn.reserve(cells.at(0).size());
  cs->at(1).conn.reserve(cells.at(1).size());
  cs->at(2).conn.reserve(cells.at(2).size());

  // For the disjoint graphs, use the selected graph given as input.
  std::vector<int>* g_dis;
  if (g_id == 0) {
    g_dis = &circ->second.first;
  }
  else
    g_dis = &circ->second.second;

  struct ent_conn e_neigh;

  for(int i = 0; i < g_dis->size(); i++) {
    // 2cell
    if (g_dis->at(i) >= (int)(cells.at(2).size())) {
      apf::ModelEntity* mdl = burn_mdl.at(1)[g_dis->at(i) - cells.at(2).size()];
      int c2_cur = m->getModelTag(mdl);
      //std::cout << "3c " << c3_cur + 1 << std::endl;
      cs->at(1).conn.push_back(c2_cur);
    }
    else if (g_dis->at(i) != -1) {
      apf::ModelEntity* mdl = burn_mdl.at(2)[g_dis->at(i)];
      //std::cout << circ->first.first[i] << " ";
      int c3_cur = m->getModelTag(mdl);
      //std::cout << "3c " << c3_cur + 1 << std::endl;
      cs->at(2).conn.push_back(c3_cur);
    }
    else {
      cs->at(2).conn.push_back(-1);
    }
  }

  // Going over the 2cells, check for every bounding 1cell, adjacent to the 
  // 0 cell, add it to cs->at(0). 
  std::map <int, bool> c1_trial;

  for(int i=0; i < cs->at(1).conn.size(); i++) {
    cb->get_conn_12(c0_tag, cs->at(1).conn.at(i), &e_neigh);
    for(int j=0; j < e_neigh.conn.size(); j++) {
      if(!c1_trial[e_neigh.conn.at(j)]) {
        cs->at(0).add_ent(e_neigh.conn.at(j));
        //std::cout << "1c" << e_neigh.conn.at(j) + 1 << " ";
        c1_trial[e_neigh.conn.at(j)] = true;
      }
    }
  }
}

/*
// Convert the cell graph information of circuits and paths into 
// gmi indices.
void vd_cell_det::update_path() {
  std::vector< cell_path >* paths = cg.get_paths();

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

}
*/
void vd_cell_det::clear_path() {

  for(int j = 0; j < circ_tup.size(); j++) {
    circ_tup.at(j).first.first.clear();
    circ_tup.at(j).first.second.clear();
    circ_tup.at(j).second.first.clear();
    circ_tup.at(j).second.second.clear();
  }
  circ_tup.clear();

  path.clear();
}

// Clear the topology and graph related containers and objects.
// Called when reloading vd_cell_det for a new 0stratum vertex.
void vd_cell_det::clear_topo() {
  cb_write_flag = true;
  cg.clear();
  c2_c3real1.clear();
  c2_c3real2.clear();
  for(int i = 0; i < cells.size(); i++) {
    cells.at(i).clear();
  }
  cells.clear();
  for(int i = 0; i < adj.size(); i++) {
    for(int j = 0; j < adj.at(i).size(); j++) {
      adj.at(i).at(j).clear();
    }
    adj.at(i).clear();
  }
  adj.clear();

  vert_c1_map.clear();
  c2_dir.clear();
  c3_dir.clear();
  clear_path();
}

// Clear the mesh related containers.
// Called when reloading vd_cell_det using the same topology but copied mesh.
void vd_cell_det::clear_mesh() {
  for(int i = 0; i < ent_sz.size(); i++) {
    ent_sz.at(i) = 0;
  }

  e_2c.clear();
  c3_edge.clear();
  for(int i = 0; i < burn_mdl.size(); i++) {
    burn_mdl.at(i).clear();
  }
  burn_mdl.clear();

  burn_mdl.resize(3);
  burn_ent.clear();
  ent_cell.clear();
  for(int i = 0; i < ent_ofire.size(); i++) {
    ent_ofire.at(i).clear();
  }
  ent_ofire.clear();

  ent_ofire.resize(2);
}

void vd_cell_det::clear() {
  clear_mesh();
  clear_topo();
}

void vd_cell_det::reload(apf::Mesh2* m_in, struct cell_base* cb_in, apf::MeshEntity* v_0cell_in) {

  cb = cb_in;
  m = m_in;

  v_0cell = v_0cell_in;
  assert(m->getModelType(m->toModel(v_0cell)) == 0);
  c0_tag = m->getModelTag(m->toModel(v_0cell));

  ent_list.change_mesh(m, v_0cell, cb);
  //ent_list.print();
  load_graph();
}

void vd_spur_mesh::get_vert() {
  int ent_type = m->getType(ent_c);
  d_ent = m->typeDimension[ent_type];

  mdl = m->toModel(ent_c);

  std::vector<apf::MeshEntity*> v_in(0);
  if(d_ent == 0) {
    v_in.push_back(ent_c);
  }
  else {
    assert(d_ent == 1);
    apf::Downward down;
    m->getDownward(ent_c, 0, down);
    v_in.reserve(2);
    v_in.push_back(down[0]);
    v_in.push_back(down[1]);
  }

  ent_list.change_mesh(m, &v_in, cb);
}

vd_spur_mesh::vd_spur_mesh(apf::Mesh2* m_in, cell_base* cb_in, apf::MeshEntity* ent_c_in) :
    ent_c(NULL), d_ent(0), mdl(NULL), cell_repl(0), c_nbr(0), m(NULL),
    cb(NULL), ent_list() {
  ent_c = ent_c_in;
  cb = cb_in;
  m = m_in;

  get_vert();
  //ent_list.print();
}
void vd_spur_mesh::reload(apf::Mesh2* m_in, cell_base* cb_in,
                                              apf::MeshEntity* ent_c_in) {
  ent_c = ent_c_in;
  cb = cb_in;
  m = m_in;

  get_vert();
  //ent_list.print();
}

std::pair<int,int> vd_spur_mesh::get_cell_repl() {

  c_nbr = 0;
  int c_dim = m->getModelType(mdl) + 1;

  for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
    if(ent_list.e.at(c_dim).at(c_id).at(c_dim).size() > 0) {
      return std::make_pair(-1,-1);
    }
  }

  c_dim++;

  for(; c_dim < 4; c_dim++) {

    cell_repl = 0;

    for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
      if(ent_list.e.at(c_dim).at(c_id).at(c_dim).size() > 0) {
        c_nbr++;
        cell_repl = c_id + 1;
      }
    }
    assert(c_nbr < 2);
    if(c_nbr == 1)
      return std::make_pair(c_dim, cell_repl);
  }

  // A entity not connected to any other higher dimensional entity is called.
  assert(c_nbr > 0);
}

// vd_3c_det

vd_3c_det::vd_3c_det() : 
    c0_tag(0), c_dim(1), d_ent(1), id(0),
    c3_nbr(0),
    m(NULL), cb(NULL), cb_dj(NULL),
    ent_list(),
    v_0cell(NULL), mdl_curr(NULL),
    cells(0, std::vector<int > (0) ),
    adj(0, std::vector<std::vector<int > > (0, std::vector<int> (0) ) ),
    ent_sz(4), burn_mdl(0, std::map<int, apf::ModelEntity*>{} ), burn_ent{},
    c3_ents(0, std::vector<apf::MeshEntity*>(0)),
    ent_cell{}, ent_ofire(2, std::vector<apf::MeshEntity* >(0)) {

  cb_dj = new cell_base(0, 0, 0, 0);
  //cb_dj->print_ent();
}

vd_3c_det::vd_3c_det(apf::Mesh2* m_in, struct cell_base* cb_in, apf::MeshEntity* v_0cell_in) : 
    c0_tag(0), c_dim(1), d_ent(1), id(0),
    c3_nbr(0),
    m(NULL), cb(NULL), cb_dj(NULL),
    ent_list(),
    v_0cell(NULL), mdl_curr(NULL),
    cells(0, std::vector<int > (0) ),
    adj(0, std::vector<std::vector<int > > (0, std::vector<int> (0) ) ),
    ent_sz(4), burn_mdl(0, std::map<int, apf::ModelEntity*>{} ), burn_ent{},
    ent_cell{}, ent_ofire(2, std::vector<apf::MeshEntity* >(0)) {

  cb = cb_in;
  m = m_in;
  v_0cell = v_0cell_in;
  int center_dim = m->getModelType(m->toModel(v_0cell));
  assert(center_dim == 0 or center_dim == 1);
  c0_tag = m->getModelTag(m->toModel(v_0cell));

  ent_sz.resize(4);

  cb_dj = new cell_base(0, 0, 0, 0);

  ent_list.change_mesh(m, v_0cell, cb);
  //ent_list.print();
  load_graph();
}

void vd_3c_det::reload(apf::Mesh2* m_in, struct cell_base* cb_in, apf::MeshEntity* v_0cell_in) {

  cb = cb_in;
  m = m_in;

  v_0cell = v_0cell_in;
  int center_dim = m->getModelType(m->toModel(v_0cell));
  assert(center_dim == 0 or center_dim == 1);
  c0_tag = m->getModelTag(m->toModel(v_0cell));


  ent_list.change_mesh(m, v_0cell, cb);
  //ent_list.print();
  load_graph();
}

void vd_3c_det::load_graph() {
  for(int i = 0; i < 4; i++) {
    ent_sz.at(i) = ent_list.es.at(i).size();
  }

  clear();

  cells.at(0).reserve(cb->get_sz(1));
  cells.at(1).reserve(cb->get_sz(2));
  cells.at(2).reserve(cb->get_sz(3));

  adj.at(0).reserve(cb->get_sz(2));
  adj.at(1).reserve(cb->get_sz(3));

  burn_new();
  c3_nbr = cells.at(2).size();

  for(int i = 0; i < c3_ents.size(); i++) {
    c3_ents.at(i).clear();
  }

  c3_ents.resize(get_3c_nbr());
  for(int i = 0; i < get_3c_nbr(); i++) {
    c3_ents.at(i).reserve(ent_list.es.at(3).size());
  }

  for(int i = 0; i < ent_list.es.at(3).size(); i++) {
    apf::MeshEntity* ent = ent_list.es.at(3).at(i);
    c3_ents.at(ent_cell[ent]).push_back(ent);
  }
}

void vd_3c_det::burn_new() {

  apf::MeshEntity* ent;
  cb_dj->clear();
  struct ent_conn e_con;

  cb_dj->add_free(0, 1);
  int c0_id = cb_dj->use_free(0);
  assert(c0_id == 0);

  c_dim = 1;
  d_ent = 1;

  id = 0;

  apf::Downward down;
  apf::MeshEntity* v_curr;

  for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
    for (int i = 0; i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
      ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

      if(!burn_ent[ent]) {
        burn_cell(c_dim, id, ent);
        cb_dj->add_free(1,1);
        int c1_id = cb_dj->use_free(1);
        assert(id == c1_id);

        e_con.clear();
        e_con.conn.push_back(c0_id);
        cb_dj->set_conn(1, c1_id, &e_con);

        m->getDownward(ent, 0, down);
        if(down[0] == v_0cell)
          v_curr = down[1];
        else
          v_curr = down[0];
        id++;
      }
    }
  }

  c_dim = 2;
  d_ent = 2;

  id = 0;

  for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
    for (int i = 0; i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
      ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

      if(!burn_ent[ent]) {
        burn_cell(c_dim, id, ent);
        cb_dj->add_free(2,1);
        int c2_id = cb_dj->use_free(2);
        assert(id == c2_id);

        e_con.clear();
        e_con.conn = adj.at(c_dim-2).back();
        cb_dj->set_conn(2, c2_id, &e_con);

        id++;
      }
    }
  }
  //coll_c2_dir();

  c_dim = 3;
  d_ent = 3;

  id = 0;

  //std::map<int, int> c2_c3map1;
  //std::map<int, int> c2_c3map2;
  // Loop over c3 tets. Burn disjoint groups of tets. Assuming there is a 
  // single c3 edge per disjoint group of entities, find the edge and 
  // associate it with the c3.
  for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
    int id_start = id;

    for (int i = 0; i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
      ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

      if(!burn_ent[ent]) {
        burn_cell(c_dim, id, ent);
        cb_dj->add_free(3,1);
        int c3_id = cb_dj->use_free(3);
        assert(id == c3_id);

        e_con.clear();
        e_con.conn = adj.at(c_dim-2).back();
        cb_dj->set_conn(3, c3_id, &e_con);
        id++;
      }
    }
  }


}

void vd_3c_det::burn_cell(int c_dim, int id, apf::MeshEntity* ent) {

  int type = m->getType(ent);
  int d = apf::Mesh::typeDimension[type];
  assert(d == c_dim);

  // Add the current disjoint cell segment to the list.
  cells.at(c_dim-1).push_back(id);

  // Associate the current disjoint cell segment with the current model.
  mdl_curr = m->toModel(ent);
  burn_mdl.at(c_dim-1)[id] = mdl_curr;

  ent_cell[ent] = id;
  burn_ent[ent] = true;

  if(c_dim == 1) {
    ent_ofire.at(0).clear();
    ent_ofire.at(1).clear();
    return;
  }

  if(c_dim > 1) {
    adj.at(c_dim-2).resize(adj.at(c_dim-2).size() + 1);
    adj.at(c_dim-2).back().reserve(cb->get_sz(c_dim-1));
  }

  apf::ModelEntity* mdl_next;
  apf::MeshEntity* ent_curr;
  apf::MeshEntity* ent_next;

  apf::Downward down;
  apf::Downward vert;
  apf::Up up;

  std::vector<int>::iterator it_l;

  // Burn the first entity of the cell dimension.
  ent_ofire.at(0).clear();
  ent_ofire.at(1).clear();

  assert(ent_sz.at(c_dim) >= 0);
  assert(ent_sz.at(c_dim-1) >= 0);
  ent_ofire.at(0).reserve(ent_sz.at(c_dim));
  ent_ofire.at(1).reserve(ent_sz.at(c_dim-1));

  ent_ofire.at(0).push_back(ent);

  bool all_burned = false;
  while(!all_burned) {
    all_burned = true;
    while(ent_ofire.at(0).size() > 0) {
      ent_curr = ent_ofire.at(0).back();

      int e_sz = m->getDownward(ent_curr, c_dim-1, down);
      for(int i = 0; i < e_sz; i++) {
        ent_next = down[i];

        int e_sz2 = m->getDownward(ent_next, 0, vert);
        if(findIn(vert, e_sz2, v_0cell) > -1) {
          mdl_next = m->toModel(ent_next);
          int c_next = m->getModelType(mdl_next);

          if(burn_ent[ent_next]) {
            int i_next = ent_cell[ent_next];
            if(c_next == c_dim - 1) {
              it_l = std::find(adj.at(c_dim-2).back().begin(), 
                          adj.at(c_dim-2).back().end(), i_next);
              if(it_l == adj.at(c_dim-2).back().end())
                adj.at(c_dim-2).back().push_back(i_next);
            }
          }
          else {
            if(c_next == c_dim) {
              assert(mdl_next == mdl_curr);

              ent_ofire.at(1).push_back(ent_next);
              burn_ent[ent_next] = true;
              ent_cell[ent_next] = id;

            }
            all_burned = false;
          }
        }
      }

      ent_ofire.at(0).pop_back();
    }

    while(ent_ofire.at(1).size() > 0) {
      ent_curr = ent_ofire.at(1).back();

      m->getUp(ent_curr, up);
      for(int i = 0; i < up.n; i++) {
        ent_next = up.e[i];

        int e_sz2 = m->getDownward(ent_next, 0, vert);
        // TODO all should be finding
        if(findIn(vert, e_sz2, v_0cell) > -1) {
          mdl_next = m->toModel(ent_next);
          int c_next = m->getModelType(mdl_next);

          if(burn_ent[ent_next]) {
            int i_next = ent_cell[ent_next];
            assert(c_next != c_dim - 1);
/*
            if(c_next == c_dim - 1) {
              it_l = std::find(adj.at(c_dim-2).back().begin(), 
                          adj.at(c_dim-2).back().end(), i_next);
              if(it_l == adj.at(c_dim-2).back().end())
                adj.at(c_dim-2).back().push_back(i_next);
            }
*/
          }
          else {
            // TODO is it not possible to have entities belonging to other 
            // cells?
            if(c_next == c_dim) {
              assert(mdl_next == mdl_curr);

              ent_ofire.at(0).push_back(ent_next);
              burn_ent[ent_next] = true;
              ent_cell[ent_next] = id;

              all_burned = false;
            }
          }
        }
      }

      ent_ofire.at(1).pop_back();
    }
  }
  std::sort(adj.at(c_dim-2).back().begin(), adj.at(c_dim-2).back().end());
}

// Get the list of same dim entities belonging to the given disjoint cell.
std::vector<apf::MeshEntity*> vd_3c_det::get_dim_ent(int dim, int id) {
  assert(dim < 4 and dim > 0);
  std::vector<apf::MeshEntity*> ents(0);
  if(ent_list.es.size() > 0) {
    ents.reserve(ent_list.es.at(dim).size());
    for(int i = 0; i < ent_list.es.at(dim).size(); i++) {
      apf::MeshEntity* ent = ent_list.es.at(dim).at(i);
      apf::ModelEntity* mdl = m->toModel(ent);

      if(m->getModelType(mdl) == dim and ent_cell[ent] == id)
        ents.push_back(ent);
    }
  }
  return ents;
}

// Return bounding disjoint strata of given stratum.
std::vector<int> vd_3c_det::get_bound_ids(int dim, int id) {
  return adj.at(dim-2).at(id);
}

// Get the list of same dim entities belonging to the given disjoint cell.
std::vector<apf::MeshEntity*> vd_3c_det::get_tet_id(int id) {
  assert(id < c3_nbr);
  return c3_ents.at(id);
}

int vd_3c_det::get_3c_nbr() {
  return c3_nbr;
}

// From local id, return the model object.
apf::ModelEntity* vd_3c_det::get_mdl(int dim_in, int id_in) {
  assert(id_in > -1 and id_in < cells.at(dim_in-1).size());
  return burn_mdl.at(dim_in-1)[id_in];
}

// From ent, return the id.
int vd_3c_det::get_id(apf::MeshEntity* ent) {

  return ent_cell[ent];
}

void vd_3c_det::clear() {
  for(int i = 0; i < cells.size(); i++) {
    cells.at(i).clear();
  }
  cells.clear();

  for(int i = 0; i < adj.size(); i++) {
    for(int j = 0; j < adj.at(i).size(); j++) {
      adj.at(i).at(j).clear();
    }
    adj.at(i).clear();
  }
  adj.clear();

  adj.resize(2);
  cells.resize(3);

  for(int i = 0; i < ent_ofire.size(); i++) {
    ent_ofire.at(i).clear();
  }
  ent_ofire.clear();

  ent_ofire.resize(2);

  for(int i = 0; i < burn_mdl.size(); i++) {
    burn_mdl.at(i).clear();
  }
  burn_mdl.clear();

  for(int i = 0; i < burn_mdl.size(); i++) {
    burn_mdl.at(i).clear();
  }
  burn_mdl.clear();

  burn_mdl.resize(3);

  for(int i = 0; i < ent_sz.size(); i++) {
    ent_sz.at(i) = 0;
  }

  burn_ent.clear();
  ent_cell.clear();
}

vd_3c_det::~vd_3c_det() {
  clear();

  ent_list.clear();

  delete cb_dj;
}


// The cell replacer replaces the cell membership of the entities with the
// new cell. 
vd_cell_repl::vd_cell_repl(apf::Mesh2* m_in) :
  m(NULL), c_map(), e_map(), repl(), 
  dest_ent(0, std::vector<apf::MeshEntity*> (0)), 
  down_list(0, std::vector<std::vector<apf::MeshEntity*> > 
           (0, std::vector<apf::MeshEntity*> (0) )),
  mdl_list(0, std::vector<apf::ModelEntity*> (0)) {

  m = m_in;
  clear();
}

// Update the mesh in the object.
void vd_cell_repl::update_m(apf::Mesh2* m_in) {
  m = m_in;
  clear();
}

// Sets the cells.
void vd_cell_repl::set_cell(std::vector<c_elem>* in_cells) {
  apf::ModelEntity* mdl_i;
  apf::ModelEntity* mdl_f;

  clear();

  for(int i = 0; i < in_cells->size(); i++) {
    mdl_i = m->findModelEntity(in_cells->at(i).first.first,
                             in_cells->at(i).first.second);
    mdl_f = m->findModelEntity(in_cells->at(i).second.first,
                             in_cells->at(i).second.second);
    c_map[mdl_i] = mdl_f;
    repl[mdl_i] = true;

    int dim_e = m->getModelType(mdl_i); 
    int tag_e = m->getModelTag(mdl_i);
    std::cout << dim_e << "c" << tag_e << " "
                << " repl " << repl[mdl_i] << " ";

    dim_e = m->getModelType(c_map[mdl_i]); 
    tag_e = m->getModelTag(c_map[mdl_i]);
    std::cout << " new_repl "
              << dim_e << "c" << tag_e << " "
              << std::endl;

  }
  for(int i = 0; i < in_cells->size(); i++) {
    mdl_i = m->findModelEntity(in_cells->at(i).first.first,
                             in_cells->at(i).first.second);
    // The assigned cell might be spurious as well. So assign the ultimate 
    // replacing cell.
    while(repl[c_map[mdl_i]] == true) {
      c_map[mdl_i] = c_map[c_map[mdl_i]];
    }
  }

}

// Replaces the cells.
void vd_cell_repl::repl_cell() {

  int const types[4] =
    {apf::Mesh::VERTEX
    ,apf::Mesh::EDGE
    ,apf::Mesh::TRIANGLE
    ,apf::Mesh::TET};

  apf::Downward down;
  apf::Downward down_l;

  apf::ModelEntity* mdl_i;
  apf::ModelEntity* mdl_f;
  apf::ModelEntity* mdl_v;

  apf::MeshEntity* ent;
  apf::MeshEntity* ent_n;

  apf::Vector3 pos;

  int e_sz;

  int count[4];
  for(int i = 0; i < 4; i++)
    count[i] = 0;

  // Count the vertices to be destroyed.
  apf::MeshIterator* it = m->begin(0);

  while(ent = m->iterate(it)) {
    mdl_i = m->toModel(ent);
    if(repl[mdl_i]) {
      count[0]++;
    }
  }
  m->end(it);

  dest_ent.at(0).reserve(count[0]);

  // Create new vertices, add the old ones to be destroyed list.
  it = m->begin(0);

  while(ent = m->iterate(it)) {
    mdl_i = m->toModel(ent);
    if(repl[mdl_i]) {
      ent_n = m->createVert(c_map[mdl_i]);

      if (m->findField("velocity_field")) {
        vd_add_0_vel_field(m, ent_n);
      }

      assert(ent_n != ent);
      m->getPoint(ent, 0, pos);
      m->setPoint(ent_n, 0, pos);
      e_map[ent] = ent_n;
      dest_ent.at(0).push_back(ent);
    }
  }
  m->end(it);

  // Count entities to be replaced.
  for(int dim = 1; dim < 4; dim++) {
    it = m->begin(dim);
    
    while(ent = m->iterate(it)) {
      bool e_repl = false;

      mdl_i = m->toModel(ent);

      e_sz = m->getDownward(ent, 0, down);
      for(int j = 0; j < e_sz; j++) {
        mdl_v = m->toModel(down[j]);

        if(repl[mdl_v]) {
          e_repl = true;
        }
      }

      for(int dim_l = 1; dim_l < dim; dim_l++) {
        e_sz = m->getDownward(ent, dim_l, down_l);
        for(int j = 0; j < e_sz; j++) {
          mdl_v = m->toModel(down_l[j]);
          if(repl[mdl_v]) {
            e_repl = true;
          }
        }
      }


      if(repl[mdl_i] or e_repl) {
        e_repl = true;
        count[dim]++;
      }
/*
      int dim_e = m->getModelType(m->toModel(ent)); 
      int tag_e = m->getModelTag(m->toModel(ent));
      std::cout << "ent_old " << ent << " " << dim_e << "c" << tag_e << " " 
                << std::endl;
      m->getDownward(ent, dim-1, down);
      for(int j = 0; j < dim+1; j++) {
        dim_e = m->getModelType(m->toModel(down[j])); 
        tag_e = m->getModelTag(m->toModel(down[j]));
        std::cout << dim-1 << "-ent " << down[j] << " " 
                  << dim_e << "c" << tag_e << " "; 
      }
      std::cout << std::endl; 
*/
    }
    //std::cout << "count[" << dim << "] " << count[dim];
    m->end(it);

    dest_ent.at(dim).resize(count[dim]);
    down_list.at(dim-1).resize(count[dim]);
    mdl_list.at(dim-1).resize(count[dim]);
    count[dim] = 0;
  }

  // Populate the model and vertex information of the new entities to be 
  // created. Add the old ones to the to be destroyed list.
  for(int dim = 3; dim > 0; dim--) {
    it = m->begin(dim);
    
    while(ent = m->iterate(it)) {
      bool e_repl = false;
      int type = m->getType(ent);

      mdl_i = m->toModel(ent);
/*
      int dim_e = m->getModelType(mdl_i); 
      int tag_e = m->getModelTag(mdl_i);
      std::cout << dim << "-ent " << ent << " "
                << dim_e << "c" << tag_e << " "
                << " repl " << repl[mdl_i] << " "
                << std::endl << "\t";
*/
      e_sz = m->getDownward(ent, 0, down);
      for(int j = 0; j < e_sz; j++) {
        mdl_v = m->toModel(down[j]);
/*
        dim_e = m->getModelType(mdl_v); 
        tag_e = m->getModelTag(mdl_v);
        std::cout << "vert " << down[j] << " "
                  << dim_e << "c" << tag_e << " "
                  << " repl " << repl[mdl_v] << " ";
*/
        if(repl[mdl_v]) {
          e_repl = true;
          down[j] = e_map[down[j]];
        }
      }
      for(int dim_l = 1; dim_l < dim; dim_l++) {
        e_sz = m->getDownward(ent, dim_l, down_l);
        for(int j = 0; j < e_sz; j++) {
          mdl_v = m->toModel(down_l[j]);
/*
          dim_e = m->getModelType(mdl_v); 
          tag_e = m->getModelTag(mdl_v);
          std::cout << "ent " << down_l[j] << " "
                    << dim_e << "c" << tag_e << " "
                    << " repl " << repl[mdl_v] << " ";
*/
          if(repl[mdl_v]) {
            e_repl = true;
          }
        }
      }

      if(repl[mdl_i] or e_repl) {
        e_repl = true;
        down_list.at(dim-1).at(count[dim]).resize(dim+1);
        for(int j = 0; j < e_sz; j++) {
          down_list.at(dim-1).at(count[dim]).at(j) = down[j];
/*
          dim_e = m->getModelType(m->toModel(down[j])); 
          tag_e = m->getModelTag(m->toModel(down[j]));
          std::cout << "vert_new " << down[j] << " "
                    << dim_e << "c" << tag_e << " ";
*/
        }
        //std::cout << std::endl;

        if(repl[mdl_i]) {
          mdl_list.at(dim-1).at(count[dim]) = c_map[mdl_i];
        }
        else {
          mdl_list.at(dim-1).at(count[dim]) = mdl_i;
        }
/*
        dim_e = m->getModelType(mdl_list.at(dim-1).at(count[dim])); 
        tag_e = m->getModelTag(mdl_list.at(dim-1).at(count[dim]));
        std::cout << "ent_new " << dim_e << "c" << tag_e << " " << std::endl;
*/
        dest_ent.at(dim).at(count[dim]) = ent;
        count[dim]++;
      }
      //std::cout << std::endl;

    }
    //std::cout << "count[" << dim << "] " << count[dim];
    m->end(it);
  }

  // Destroy the old entities and vertices.
  for(int dim = 3; dim > -1; dim--) {
    for(int i = 0; i < dest_ent.at(dim).size(); i++) {
/*
      int dim_e = m->getModelType(m->toModel(dest_ent.at(dim).at(i))); 
      int tag_e = m->getModelTag(m->toModel(dest_ent.at(dim).at(i)));

      std::cout << "Destroying " << dim << "-ent " 
                << dest_ent.at(dim).at(i) << " " 
                << dim_e << "c" << tag_e << " " 
                << std::endl;
*/
      m->destroy(dest_ent.at(dim).at(i));
    }
  }

  // Create the new entities.
  for(int dim = 1; dim < 4; dim++) {

    for(int i = 0; i < mdl_list.at(dim-1).size(); i++) {
      for(int j = 0; j < dim+1; j++) {
        down[j] = down_list.at(dim-1).at(i).at(j);
      }
      ent_n = buildElement(m, mdl_list.at(dim-1).at(i), types[dim], down);
/*
      int dim_e = m->getModelType(m->toModel(ent_n)); 
      int tag_e = m->getModelTag(m->toModel(ent_n));
      std::cout << "ent_new " << ent_n << " " << dim_e << "c" << tag_e << " " 
                << std::endl;
      m->getDownward(ent_n, dim-1, down);
      for(int j = 0; j < dim+1; j++) {
        dim_e = m->getModelType(m->toModel(down[j])); 
        tag_e = m->getModelTag(m->toModel(down[j]));
        std::cout << dim-1 << "-ent " << down[j] << " " 
                  << dim_e << "c" << tag_e << " "; 
      }
      std::cout << std::endl; 
*/

    }
  }
  m->acceptChanges();
  std::cout << "Cell entities replaced." << std::endl;
}

// Clear the map.
void vd_cell_repl::clear() {

  c_map.clear();
  repl.clear();
  e_map.clear();

  for(int i = 0; i < dest_ent.size(); i++) {
    dest_ent.at(i).clear();
  }
  dest_ent.clear();
  dest_ent.resize(4);

  for(int i = 0; i < down_list.size(); i++) {
    for(int j = 0; j < down_list.at(i).size(); j++) {
      down_list.at(i).at(j).clear();
    }
    down_list.at(i).clear();
  }
  down_list.clear();
  down_list.resize(3);

  for(int i = 0; i < mdl_list.size(); i++) {
    mdl_list.at(i).clear();
  }
  mdl_list.clear();
  mdl_list.resize(3);

}

vd_cell_val::vd_cell_val() : 
    cell_dim(0), cell_tag(0), c_dim(1), d_ent(1), id(0), ent_list(),
    v_ctr(NULL), mdl_curr(NULL), valid(true),
    cells(0, std::vector<int > (0) ),
    adj(0, std::vector<std::vector<int > > (0, std::vector<int> (0) ) ),
    ent_sz(4), burn_mdl(0, std::map<int, apf::ModelEntity*>{} ), burn_ent{},
    ent_cell{}, ent_ofire(2, std::vector<apf::MeshEntity* >(0)) {

  cb_dj = new cell_base(0, 0, 0, 0);
  //cb_dj->print_ent();
}

void vd_cell_val::burn_cell(int c_dim, int id, apf::MeshEntity* ent) {

  int type = m->getType(ent);
  int d = apf::Mesh::typeDimension[type];
  assert(d == c_dim);

  // Add the current disjoint cell segment to the list.
  cells.at(c_dim-1).push_back(id);

  // Associate the current disjoint cell segment with the current model.
  mdl_curr = m->toModel(ent);
  burn_mdl.at(c_dim-1)[id] = mdl_curr;

  ent_cell[ent] = id;
  burn_ent[ent] = true;

  if(c_dim == 1) {
    ent_ofire.at(0).clear();
    ent_ofire.at(1).clear();
    return;
  }

  // TODO pushback rather than resize? 
  if(c_dim > 1) {
    adj.at(c_dim-2).resize(adj.at(c_dim-2).size() + 1);
    adj.at(c_dim-2).back().reserve(cb->get_sz(c_dim-1));
  }

  apf::ModelEntity* mdl_next;
  apf::MeshEntity* ent_curr;
  apf::MeshEntity* ent_next;

  apf::Downward down;
  apf::Downward vert;
  apf::Up up;

  std::vector<int>::iterator it_l;


  // Burn the first entity of the cell dimension.
  ent_ofire.at(0).clear();
  ent_ofire.at(1).clear();

  assert(ent_sz.at(c_dim) >= 0);
  assert(ent_sz.at(c_dim-1) >= 0);
  ent_ofire.at(0).reserve(ent_sz.at(c_dim));
  ent_ofire.at(1).reserve(ent_sz.at(c_dim-1));

  ent_ofire.at(0).push_back(ent);

  bool all_burned = false;
  while(!all_burned) {
    all_burned = true;
    while(ent_ofire.at(0).size() > 0) {
      ent_curr = ent_ofire.at(0).back();

      int e_sz = m->getDownward(ent_curr, c_dim-1, down);
      for(int i = 0; i < e_sz; i++) {
        ent_next = down[i];

        int e_sz2 = m->getDownward(ent_next, 0, vert);
        if(findIn(vert, e_sz2, v_ctr) > -1) {
          mdl_next = m->toModel(ent_next);
          int c_next = m->getModelType(mdl_next);

          if(burn_ent[ent_next]) {
            int i_next = ent_cell[ent_next];
            if(c_next == c_dim - 1) {
              it_l = std::find(adj.at(c_dim-2).back().begin(), 
                          adj.at(c_dim-2).back().end(), i_next);
              if(it_l == adj.at(c_dim-2).back().end())
                adj.at(c_dim-2).back().push_back(i_next);
            }
          }
          else {
            if(c_next == c_dim) {
              assert(mdl_next == mdl_curr);

              ent_ofire.at(1).push_back(ent_next);
              burn_ent[ent_next] = true;
              ent_cell[ent_next] = id;

            }
            all_burned = false;
          }
        }
      }

      ent_ofire.at(0).pop_back();
    }

    while(ent_ofire.at(1).size() > 0) {
      ent_curr = ent_ofire.at(1).back();

      m->getUp(ent_curr, up);
      for(int i = 0; i < up.n; i++) {
        ent_next = up.e[i];

        int e_sz2 = m->getDownward(ent_next, 0, vert);
        int e_i = findIn(vert, e_sz2, v_ctr);
        // TODO all should be finding
        if(e_i > -1) {
          mdl_next = m->toModel(ent_next);
          int c_next = m->getModelType(mdl_next);

          if(burn_ent[ent_next]) {
            int i_next = ent_cell[ent_next];
            assert(c_next != c_dim - 1);
/*
            if(c_next == c_dim - 1) {
              it_l = std::find(adj.at(c_dim-2).back().begin(), 
                          adj.at(c_dim-2).back().end(), i_next);
              if(it_l == adj.at(c_dim-2).back().end())
                adj.at(c_dim-2).back().push_back(i_next);
            }
*/
          }
          else {
            // TODO is it not possible to have entities belonging to other 
            // cells?
            if(c_next == c_dim) {
              assert(mdl_next == mdl_curr);

              ent_ofire.at(0).push_back(ent_next);
              burn_ent[ent_next] = true;
              ent_cell[ent_next] = id;

              all_burned = false;
            }
          }
        }
      }

      ent_ofire.at(1).pop_back();
    }
  }
  std::sort(adj.at(c_dim-2).back().begin(), adj.at(c_dim-2).back().end());
}


vd_cell_val::vd_cell_val(apf::Mesh2* m_in, struct cell_base* cb_in, apf::MeshEntity* v_ctr_in) : 
    cell_dim(0), cell_tag(0), c_dim(1), d_ent(1), id(0), ent_list(),
    v_ctr(NULL), mdl_curr(NULL), valid(true),
    cells(0, std::vector<int > (0) ),
    adj(0, std::vector<std::vector<int > > (0, std::vector<int> (0) ) ),
    ent_sz(4), burn_mdl(0, std::map<int, apf::ModelEntity*>{} ), burn_ent{},
    ent_cell{}, ent_ofire(2, std::vector<apf::MeshEntity* >(0)) {

  cb = cb_in;
  m = m_in;
  v_ctr = v_ctr_in;
  cell_dim = m->getModelType(m->toModel(v_ctr));
  cell_tag = m->getModelTag(m->toModel(v_ctr));

  ent_sz.resize(4);

  cb_dj = new cell_base(0, 0, 0, 0);

  ent_list.change_mesh(m, v_ctr, cb);
  //ent_list.print();
  load_graph();
}

void vd_cell_val::load_graph() {

  for(int i = 0; i < 4; i++) {
    ent_sz.at(i) = ent_list.es.at(i).size();
  }

  clear();

  cells.at(0).reserve(cb->get_sz(1));
  cells.at(1).reserve(cb->get_sz(2));
  cells.at(2).reserve(cb->get_sz(3));

  adj.at(0).reserve(cb->get_sz(2));
  adj.at(1).reserve(cb->get_sz(3));

  valid = true;

  burn_new();

  std::cout << "Cells " << std::endl;
  for(int i = 0; i < cells.size(); i++) {
    for(int j = 0; j < cells.at(i).size(); j++) {
      if(cells.at(i).at(j) < 0)
        std::cout << i+1 << "Ext " << cells.at(i).at(j) << " ";
      else
        std::cout << i+1 << "c" << cells.at(i).at(j) << " ";
    }
    std::cout << std::endl;
  }

}

int vd_cell_val::find_cell_dj(int dim, apf::MeshEntity* ent, std::vector<int> bound) {
}

void vd_cell_val::burn_new() {
  apf::MeshEntity* ent;
  cb_dj->clear();
  struct ent_conn e_con;

  if(cell_dim == 0) {
    cb_dj->add_free(0, 1);
    int c0_id = cb_dj->use_free(0);
    assert(c0_id == 0);
  }

  c_dim = 1;
  d_ent = 1;

  id = 0;

  apf::Downward down;
  apf::MeshEntity* v_curr;

  if(cell_dim == 0) {
    for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
      if(ent_list.e.at(c_dim).at(c_id).at(d_ent).size() > 1) {
        ent_conn e_con;
        cb->get_conn(1, c_id, &e_con);
        if(e_con.conn.size() != 1) {
          valid = false;
          std::cout << "Invalid" << cell_dim << "c" << cell_tag 
                    << " is connected to "
                    << "1c" << c_id + 1 << "multiple times" 
                    << std::endl;
        }
      }
      for (int i = 0; 
              i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
        ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

        if(!burn_ent[ent]) {
          burn_cell(c_dim, id, ent);
          cb_dj->add_free(1,1);
          int c1_id = cb_dj->use_free(1);
          assert(id == c1_id);

          e_con.clear();
          e_con.conn.push_back(0);
          cb_dj->set_conn(1, c1_id, &e_con);

          m->getDownward(ent, 0, down);
          if(down[0] == v_ctr)
            v_curr = down[1];
          else
            v_curr = down[0];
          //std::cout << m->getModelType(m->toModel(ent)) << "c"
          //          << m->getModelTag(m->toModel(ent))
          //          << " " << v_curr << " " << id << std::endl;

          id++;
        }
      }
    }
  }
  else {
    for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
      if(ent_list.e.at(c_dim).at(c_id).at(d_ent).size() > 0) {
        if((cell_dim != c_dim) and cell_tag != c_id+1) {
          valid = false;
          std::cout << "Invalid" << cell_dim << "c" << cell_tag 
                    << " is connected to "
                    << "1c" << c_id + 1 << std::endl;
        }
      }
    }
  }

  c_dim = 2;
  d_ent = 2;

  id = 0;

  if(cell_dim < 2) {
    for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
      for (int i = 0; 
              i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
        ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

        if(!burn_ent[ent]) {
          burn_cell(c_dim, id, ent);
          cb_dj->add_free(2,1);
          int c2_id = cb_dj->use_free(2);
          assert(id == c2_id);


          e_con.clear();
          e_con.conn = adj.at(c_dim-2).back();
          cb_dj->set_conn(2, c2_id, &e_con);

          id++;
        }
      }
    }
  }
  else {
    for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
      if(ent_list.e.at(c_dim).at(c_id).at(d_ent).size() > 0) {
        if((cell_dim != c_dim) and cell_tag != c_id+1) {
          valid = false;
          std::cout << "Invalid" << cell_dim << "c" << cell_tag 
                    << " is connected to "
                    << "2c" << c_id + 1 << std::endl;
        }
      }
    }
  }

  c_dim = 3;
  d_ent = 3;

  id = 0;

  //std::map<int, int> c2_c3map1;
  //std::map<int, int> c2_c3map2;

  if(cell_dim < 3) {
    for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
      int id_start = id;

      for (int i = 0; 
              i < ent_list.e.at(c_dim).at(c_id).at(d_ent).size(); i++) {
        ent = ent_list.e.at(c_dim).at(c_id).at(d_ent).at(i);

        if(!burn_ent[ent]) {
          burn_cell(c_dim, id, ent);
          cb_dj->add_free(3,1);
          int c3_id = cb_dj->use_free(3);
          assert(id == c3_id);


          e_con.clear();
          e_con.conn = adj.at(c_dim-2).back();
          cb_dj->set_conn(3, c3_id, &e_con);
          id++;
        }
      }
    }
  }
  else {
    for (int c_id = 0; c_id < cb->get_sz(c_dim); c_id++) {
      if(ent_list.e.at(c_dim).at(c_id).at(d_ent).size() > 0) {
        if((cell_dim != c_dim) and cell_tag != c_id+1) {
          valid = false;
          std::cout << "Invalid" << cell_dim << "c" << cell_tag 
                    << " is connected to "
                    << "3c" << c_id + 1 << std::endl;
        }
      }
    }
  }

}

std::vector<int> vd_cell_val::burn_bound(int c_dim, apf::MeshEntity* ent) {
}

vd_cell_val::~vd_cell_val() {
  clear();
}

void vd_cell_val::clear() {

  for(int i = 0; i < cells.size(); i++) {
    cells.at(i).clear();
  }
  cells.clear();

  for(int i = 0; i < adj.size(); i++) {
    for(int j = 0; j < adj.at(i).size(); j++) {
      adj.at(i).at(j).clear();
    }
    adj.at(i).clear();
  }
  adj.clear();

  adj.resize(2);
  cells.resize(3);

  for(int i = 0; i < ent_ofire.size(); i++) {
    ent_ofire.at(i).clear();
  }
  ent_ofire.clear();

  ent_ofire.resize(2);

  for(int i = 0; i < burn_mdl.size(); i++) {
    burn_mdl.at(i).clear();
  }
  burn_mdl.clear();

  for(int i = 0; i < burn_mdl.size(); i++) {
    burn_mdl.at(i).clear();
  }
  burn_mdl.clear();

  burn_mdl.resize(3);

  for(int i = 0; i < ent_sz.size(); i++) {
    ent_sz.at(i) = 0;
  }

  burn_ent.clear();
  ent_cell.clear();
}

void vd_cell_val::reload(apf::Mesh2* m_in, struct cell_base* cb_in,
                                              apf::MeshEntity* v_ctr_in) {
  cb = cb_in;
  m = m_in;
  v_ctr = v_ctr_in;
  cell_dim = m->getModelTag(m->toModel(v_ctr));
  cell_tag = m->getModelTag(m->toModel(v_ctr));

  ent_sz.resize(4);

  cb_dj = new cell_base(0, 0, 0, 0);

  ent_list.change_mesh(m, v_ctr, cb);
  //ent_list.print();
  load_graph();
}

std::vector<apf::MeshEntity*> vd_cell_val::get_dim_ent(int dim, int id) {
  assert(dim < 4 and dim > 0);
  std::vector<apf::MeshEntity*> ents(0);
  ents.reserve(ent_list.es.at(dim).size());
  for(int i = 0; i < ent_list.es.at(dim).size(); i++) {
    apf::MeshEntity* ent = ent_list.es.at(dim).at(i);
    apf::ModelEntity* mdl = m->toModel(ent);

    if(m->getModelType(mdl) == dim and ent_cell[ent] == id)
      ents.push_back(ent);
  }
  return ents;
}

apf::ModelEntity* vd_cell_val::get_mdl(int dim_in, int id_in) {
  assert(id_in > -1 and id_in < cells.at(dim_in-1).size());
  return burn_mdl.at(dim_in-1)[id_in];
}

bool vd_cell_val::vd_vert_valid() {
  return valid;
}

