#include <limits>       /* numeric_limits */
#include <functional>
#include <cstring>

#include <vector>
#include <map>
#include <deque>
#include <algorithm>    

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include "ma.h"
#include "maShape.h"
#include "maSize.h"

#include <PCU.h>

#include <gmi.h>
#include <gmi_mesh.h>

/* Topology manipulator:  */
#include "topo_manip.h"
#include "topo_disp.h"
/* Topology write:  */
#include "topo_write.h"

#include "topo_topo.h"

#include "topo_graph.h"

#include "topo_lens.h"

#include "topo_vd.h"

#include "topo_geom.h"

#include "topo_edisc.h"

#include "topo_ma.h"

#include "topo_tmesh.h"

#include "topo_tess.h"

#include "topo_feat.h"

/*
// Here is a template for meshadapt:
class Linear : public ma::IsotropicFunction
{
  public:
    Linear(ma::Mesh* m, double sz)
    {
      sizing = sz;

      mesh = m;
      average = ma::getAverageEdgeLength(m);
    }
    virtual double getValue(ma::Entity* v)
    {
      return average*sizing;
    }
  private:
    ma::Mesh* mesh;
    double average;
    double sizing;
};
*/
// Edge id across an edge within a tet.
int lookup_t_e_x_e [6] = {5, 3, 4, 1, 2, 0};

// Given an edge within a tet return adjacent triangles.
int lookup_t_e_n_t [6][2] = {{0,1},{0,2},{3,0},{1,3},{2,1},{3,2}};
// Given an edge within a tet return vertices across the edge and on the 
// triangles in lookup_t_e_n_t.
int lookup_t_e_n_g [6][2] = {{2,3},{0,3},{3,1},{1,2},{2,0},{0,1}};

// Given two triangle indices, find the edge not contained within one of the
// triangles.
int lookup_tt_e [4][4] = {{-1,5,3,4},{5,-1,2,1},{3,2,-1,0},{4,1,0,-1}};

// Given the tri, return the edges across within the tetrahedron.
int lookup_t_tet_e_x [4][3] = {{3,4,5},{1,2,5},{0,2,3},{0,1,4}};

// Used in adapt_prob_edge_sp. For a low quality tet, this gives which entity
// should be used to increase quality by the split and later by a quality 
// improvement operation.
// Given the edge e_i, find the triangles that are joined by the edge. Find the 
// e1 across e_i. For the triangles, find the other tets t2 and t3. Find edges 
// e2 and e3 across e_i on t2 and t3. Split e1, e2, and e3.
// 0 Split nothing. Internal tet, should be fixed by MA.
// 1 0th tri: Longest edge
// 2 1st tri: Longest edge
// 3 Edge across the the tris 0 and 1: edge 5
// 4 2nd tri: Longest edge
// 5 Edge across the the tris 0 and 2: edge 3
// 6 Edge across the the tris 1 and 2: edge 2
// 7 Split 3rd tri
// 8 3rd tri: Longest edge
// 9 Edge across the the tris 0 and 3: edge 4
// 10 Edge across the the tris 1 and 3: edge 1
// 11 Split 2nd tri
// 12 Edge across the the tris 2 and 3: edge 0
// 13 Split 1st tri
// 14 Split 0th tri
// 15 Split nothing. the last tet in a volume, strangely very low quality
// 6 + n, n = [0, 3]. Split tri_{n}
// n, n < 6. Split edge_{n}
// -2, split nothing
// -1, split longest edge
//int tri_split_map [16] = {-2, -1, -1, 5, 
//                          -1, 3, 2, 9, 
//                          -1, 4, 1, 8,
//                           0, 7, 6, -2};
// A different scheme. When single boundary triangle, split the boundary triangle.
int tri_split_map [16] = {-2, 6, 7, 5, 
                           8, 3, 2, 9, 
                           9, 4, 1, 8,
                           0, 7, 6, -2};
// -1: Nothing - longest edge
// 0: Split edge
// 1: Split tri
// 1: Split tri if largest area
//int tri_split_type [16] = {-1, -1, -1, 0, 
//                          -1, 0, 0, 1, 
//                          -1, 0, 0, 1,
//                           0, 1, 1, -1};
int tri_split_type [16] = {-1, 2, 2, 0, 
                            2, 0, 0, 1, 
                            2, 0, 0, 1,
                            0, 1, 1, -1};
int tri_single_2c [16] =  {0, 1, 1, 0, 
                           1, 0, 0, 1, 
                           0, 0, 0, 0,
                           0, 0, 0, 0};
int tri_split_shift [4] = {1, 2, 4, 8};

// Similar to tri split type, used to identify which edge to use to determine 
// the edges to be split for the three tetrahedra.
int edge_split_type [16] = {-1, 2, 2, 0, 
                            2, 0, 0, 1, 
                            2, 0, 0, 1,
                            0, 1, 1, -1};


class Ref_Bound : public ma::IsotropicFunction
{
  public:
    Ref_Bound(ma::Mesh* m, double sz)
    {
      sizing = sz;

      mesh = m;
      average = ma::getAverageEdgeLength(m);
    }
    virtual double getValue(ma::Entity* v)
    {
      int type = mesh->getModelType(mesh->toModel(v));
      if(type < 3)
        return average*sizing/3;
      else
        return average*sizing;
    }
  private:
    ma::Mesh* mesh;
    double average;
    double sizing;
};

bool vd_sim::check_ma() {

  apf::MeshEntity* elem;
  apf::MeshIterator* it = m->begin(3);
  while(elem = m->iterate(it)) {
    apf::MeshElement* ee = createMeshElement(m, elem);
    double meas = measure(ee);
    destroyMeshElement(ee);
    if (meas < std::numeric_limits<double>::min() ) {
      printf("Volume %p is %f. Reverting changes and reducing time step.\n",
                                                     (void*)elem, meas);
      std::cout << getLinearCentroid(m, elem) << std::endl;
      return true;
    }
  }
  m->end(it);
  return false;
}

bool vd_sim::check_ma(std::vector<apf::MeshEntity*>* verts) {

  std::vector<apf::MeshEntity*> edges(0);
  std::vector<apf::MeshEntity*> tris(0);
  std::vector<apf::MeshEntity*> tets(0);

  vd_set_up(m, verts, &edges);
  vd_set_up(m, &edges, &tris);
  vd_set_up(m, &tris, &tets);

  apf::MeshEntity* elem;
  for(int i = 0; i < tets.size(); i++) {
    apf::MeshElement* ee = createMeshElement(m, tets.at(i));
    double meas = measure(ee);
    destroyMeshElement(ee);
    if (meas < std::numeric_limits<double>::min() ) {
      printf("Volume %p is %f. Reverting changes and reducing time step.\n",
                                                     (void*)tets.at(i), meas);
      std::cout << getLinearCentroid(m, tets.at(i)) << std::endl;
      return true;
    }
  }
  return false;
}

bool vd_sim::check_ma_sgn() {

  apf::Downward down;

  apf::MeshEntity* elem;
  apf::MeshIterator* it = m->begin(3);
  while(elem = m->iterate(it)) {
    m->getDownward(elem, 0, down);
    bool sgn = vd_volume_tet_sign(m, down);
    if (!sgn) {
      double meas = vd_volume_tet(m, down);
      apf::ModelEntity* mdl = m->toModel(elem);
      std::cout << "Elem " << elem << " " << m->getModelType(mdl) << "c"
            << m->getModelTag(mdl) << " is inverted with vol "
            << meas << std::endl;
      std::cout << getLinearCentroid(m, elem) << std::endl;
      return true;
    }
  }
  m->end(it);
  return false;
}

bool vd_sim::check_ma_sgn(std::vector<apf::MeshEntity*>* verts) {

  apf::Downward down;

  std::vector<apf::MeshEntity*> edges(0);
  std::vector<apf::MeshEntity*> tris(0);
  std::vector<apf::MeshEntity*> tets(0);

  vd_set_up(m, verts, &edges);
  vd_set_up(m, &edges, &tris);
  vd_set_up(m, &tris, &tets);

  apf::MeshEntity* elem;
  for(int i = 0; i < tets.size(); i++) {
    m->getDownward(tets.at(i), 0, down);
    bool sgn = vd_volume_tet_sign(m, down);
    if (!sgn) {
      double meas = vd_volume_tet(m, down);
      apf::ModelEntity* mdl = m->toModel(tets.at(i));
      std::cout << "Elem " << tets.at(i) << " " << m->getModelType(mdl) << "c"
            << m->getModelTag(mdl) << " is inverted with vol "
            << meas << std::endl;
      std::cout << getLinearCentroid(m, tets.at(i)) << std::endl;
      return true;
    }
  }
  return false;
}

void vd_sim::iter_step() {

  double iter_time = time_curr;

  //vd_att_disp_field(m);
  f_calc.vd_att_fields(m);

  if(ext_gam_flag == EXT_GAM_TYPE::ZERO) {
    int sz = c_base->get_sz(2);
    for(int i = 0; i < sz; i++) {
      if(!c_base->is_free(2,i) and c_base->get_cell_ext(2, i)) {
        f_calc.upd_gam2(m, e_list.e.at(2).at(i).at(2), 0);
      }
    }
  }

  f_calc.vd_calc_vel(m);

  //double vel_max = vd_find_max_vel(m);
  //f_calc.vdparam.adj_dt(std::min(dt_set, ma::getAverageEdgeLength(m)/vel_max));
  //f_calc.vdparam.adj_dt(std::min(f_calc.find_min_t2(m)/9, dt_set));
  f_calc.vdparam.adj_dt(std::min(f_calc.find_min_t2(m)/t_rat, dt_set));

  if(f_calc.vdparam.dt < std::numeric_limits<double>::min())
    f_calc.vdparam.adj_dt(dt_set);

  printf("Time: %.2f, End: %.2f:\n", time_curr, time_end);
  bool sub_suc = false;
  // TODO apf field has an operation to apply displacement field and revert it.
  double sub_param = 1;
  f_calc.set_drag_glob(getAverageEntSize(m, 2));

  for (int i = 0; i < iter_sub; i++) {
    printf("Trial %d.\n", i);

    f_calc.vd_calc_vel(m);

    std::stringstream ss;
    ss <<"output/before_trial" << i;

    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    save_vtk_name(cstr);

    //vd_tet_vel(m);
    //vd_find_max_vel(m);
    t_sub = f_calc.find_min_t2(m);
    // double t_t = f_calc.find_min_t_tet(m);
    //f_calc.vdparam.adj_dt(std::min(vd_find_min_t(m), f_calc.vdparam.dt));
    f_calc.vdparam.adj_dt(std::min(t_sub/t_rat/sub_param, dt_set));
    if(f_calc.vdparam.dt < std::numeric_limits<double>::min())
      f_calc.vdparam.adj_dt(dt_set);
    printf("dt: %.8f, dt_set: %.8f, t_sub: %.8f.\n", f_calc.vdparam.dt, 
                                                    dt_set, t_sub);
    if(!fixed_time and t_sub > 3*dt_set)
      dt_set = t_sub;

    vd_apply_vel_field(m, f_calc.vdparam.dt);
    //displaceMesh(Mesh2* m, Field* d, double factor)

    m->verify();
    ss.str("");

    ss.clear();
    ss <<"output/after_trial" << i;
    tmp = ss.str();
    cstr = tmp.c_str();

    save_vtk_name(cstr);

    //PCU_Thrd_Barrier();
    bool chk_ma = check_ma_sgn();
    //PCU_Thrd_Barrier();

    if (PCU_Or(chk_ma)) {
      vd_apply_vel_field(m, -f_calc.vdparam.dt);
      sub_param = sub_param*2;
      //PCU_Thrd_Barrier();
    }
    else {
      time_curr = time_curr + f_calc.vdparam.dt;
      sub_suc = true;
    }
    vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_ITER);
    vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
    if(extract_flag) {
      extract_data();
      //extract_MS();
    }
    vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
    vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_ITER);
  }
  assert(sub_suc);
}

// Move only a select set of boundaries. Mainly used to circumnavigate the  
// gsl memory issue for large meshes. 
// TODO Better incorporate into the overall code structure. Ideally, have an 
// input interface for every option including the way motion is handled.
void vd_sim::iter_cells() {

  double iter_time = time_curr;

  //vd_att_disp_field(m);
  f_calc.vd_att_fields(m);

  if(ext_gam_flag == EXT_GAM_TYPE::ZERO) {
    int sz = c_base->get_sz(2);
    for(int i = 0; i < sz; i++) {
      if(!c_base->is_free(2,i) and c_base->get_cell_ext(2, i)) {
        f_calc.upd_gam2(m, e_list.e.at(2).at(i).at(2), 0);
      }
    }
  }
  std::vector<apf::MeshEntity*> v_2move(0);
  get_v2move(v_2move);

  f_calc.vd_att_fields(m);
  f_calc.vd_calc_vel(m, &v_2move);

  //double vel_max = vd_find_max_vel(m);
  //f_calc.vdparam.adj_dt(std::min(dt_set, ma::getAverageEdgeLength(m)/vel_max));
  //f_calc.vdparam.adj_dt(std::min(f_calc.find_min_t2(m)/9, dt_set));
  f_calc.vdparam.adj_dt(std::min(f_calc.find_min_t2(m, &v_2move)/t_rat, dt_set));

  if(f_calc.vdparam.dt < std::numeric_limits<double>::min())
    f_calc.vdparam.adj_dt(dt_set);

  printf("Time: %.2f, End: %.2f:\n", time_curr, time_end);
  bool sub_suc = false;
  // TODO apf field has an operation to apply displacement field and revert it.
  double sub_param = 1;
  f_calc.set_drag_glob(getAverageEntSize(m, 2));

  for (int i = 0; i < iter_sub; i++) {
    printf("Trial %d.\n", i);

    f_calc.vd_att_fields(m);
    f_calc.vd_calc_vel(m, &v_2move);

    std::stringstream ss;
    ss <<"output/before_trial" << i;

    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    save_vtk_name(cstr);

    //vd_tet_vel(m);
    //vd_find_max_vel(m);
    t_sub = f_calc.find_min_t2(m, &v_2move);
    // double t_t = f_calc.find_min_t_tet(m);
    //f_calc.vdparam.adj_dt(std::min(vd_find_min_t(m), f_calc.vdparam.dt));
    f_calc.vdparam.adj_dt(std::min(t_sub/t_rat/sub_param, dt_set));
    if(f_calc.vdparam.dt < std::numeric_limits<double>::min())
      f_calc.vdparam.adj_dt(dt_set);
    printf("dt: %.8f, dt_set: %.8f, t_sub: %.8f.\n", f_calc.vdparam.dt, 
                                                    dt_set, t_sub);
    if(!fixed_time and t_sub > 3*dt_set)
      dt_set = t_sub;

    vd_apply_vel_field(m, f_calc.vdparam.dt);
    //displaceMesh(Mesh2* m, Field* d, double factor)

    m->verify();
    ss.str("");

    ss.clear();
    ss <<"output/after_trial" << i;
    tmp = ss.str();
    cstr = tmp.c_str();

    save_vtk_name(cstr);

    //PCU_Thrd_Barrier();
    bool chk_ma = check_ma_sgn();
    //PCU_Thrd_Barrier();

    if (PCU_Or(chk_ma)) {
      vd_apply_vel_field(m, -f_calc.vdparam.dt);
      sub_param = sub_param*2;
      //PCU_Thrd_Barrier();
    }
    else {
      time_curr = time_curr + f_calc.vdparam.dt;
      sub_suc = true;
    }

    vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_ITER);

    double elapsed_temp = vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_AUX);
    vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
    if(extract_flag) {
      extract_data();
      //extract_MS();
    }
    double elapsed = vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
    std::cout << "extract took " << elapsed-elapsed_temp << std::endl;
    vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_ITER);

  }
  assert(sub_suc);
}

bool vd_sim::chk_vert_val(apf::MeshEntity* vert) {
  vd_cell_val val_vert(m, c_base, vert);
  return val_vert.vd_vert_valid();
}

void vd_sim::shr_axis(apf::Vector3 ax, double shr) {

  vd_shr_axis(m, ax, shr);
  if(proj_flag == (int) PROJ_TYPE::EXT_SHELL) {

    std::map<int, bool> shell_map{};
    for(int i = 0; i < c_base->get_sz(0); i++) {
      if(c_base->get_cell_ext(0, i) and e_sh.chk_shell(0, i)) {
        shell sh = e_sh.get_shell(0, i);
        if(sh.dim == 0 and !shell_map[sh.id]) {
          shell_map[sh.id] = true;

          e_sh.set_pos(sh, vd_shr_vec(e_sh.get_shell_pos(sh), ax, shr));
          //e_sh.set_dir(sh, vd_shr_vec(get_shell_dir(sh), ax, shr));
        }
      }
    }
    shell_map.clear();

    for(int i = 0; i < c_base->get_sz(1); i++) {
      if(c_base->get_cell_ext(1, i) and e_sh.chk_shell(1, i)) {
        shell sh = e_sh.get_shell(1, i);
        if(sh.dim == 1 and !shell_map[sh.id]) {
          shell_map[sh.id] = true;

          e_sh.set_pos(sh, vd_shr_vec(e_sh.get_shell_pos(sh), ax, shr));
          e_sh.set_dir(sh, norm_0(vd_shr_vec(e_sh.get_shell_dir(sh), ax, shr)));
        }
      }
    }

    shell_map.clear();

    for(int i = 0; i < c_base->get_sz(2); i++) {
      if(c_base->get_cell_ext(2, i) and e_sh.chk_shell(2, i)) {
        shell sh = e_sh.get_shell(2, i);
        if(sh.dim == 2 and !shell_map[sh.id]) {
          shell_map[sh.id] = true;

          e_sh.set_pos(sh, vd_shr_vec(e_sh.get_shell_pos(sh), ax, shr));
          e_sh.set_dir(sh, norm_0(
                              vd_shr_vec(e_sh.get_shell_dir(sh), ax, 1/shr)));

        }
      }
    }
  }
}

// Given a same dim list of 0-,1-,2-cells and an axis, take the weighted 
// average of the velocity. Evolve the vertices of the boundary cells until
// the weighted average of the relative velocity is below a threshold of the  
// weighted average.
void vd_sim::evolve_bound_ax_conv(std::vector<std::pair<int, int> > * cells, 
                          apf::Vector3 ax, double th) {
  f_calc.vd_att_fields(m);

  if(th > 1 or th < - std::numeric_limits<double>::min())
    th = 1;

  if(cells->size() == 0)
    return;
  int dim = cells->at(0).first;
  int v_nbr = dim + 1;
  for(int i = 1; i < cells->size(); i++) {
    if(dim != cells->at(i).first)
      return;
  }
  e_list.refresh();

  apf::Field* vel_field = m->findField("velocity_field");

  std::vector<apf::MeshEntity*> verts(0);
  std::vector<apf::MeshEntity*> ents(0);
  std::vector<std::vector<apf::MeshEntity*> > ent_verts(0, 
                              std::vector<apf::MeshEntity*>(0));
  std::vector<apf::MeshEntity*> temp_vert(0);
  // Collect the entities and their downward vertices:
  int e_nbr = 0;
  for(int i = 0; i < cells->size(); i++) {
    int id = cells->at(i).second - 1;
    std::vector<apf::MeshEntity*> * temp = &e_list.e.at(dim).at(id).at(dim);
    e_nbr = e_nbr + temp->size();
  }

  ents.reserve(e_nbr);
  ent_verts.reserve(e_nbr);
  for(int i = 0; i < cells->size(); i++) {
    int id = cells->at(i).second - 1;
    std::vector<apf::MeshEntity*> * temp = &e_list.e.at(dim).at(id).at(dim);
    for(int j = 0; j < temp->size(); j++) {
      vd_set_down(m, temp->at(j), &temp_vert, dim);
      ents.push_back(temp->at(j));
      ent_verts.push_back(temp_vert);
    }
  }
  vd_set_down(m, &ents, &verts, dim);

  // Update the velocities:
  f_calc.vd_upd_vel_field(m, &verts);

  apf::Vector3 v_temp(0,0,0);
  // Calculate the weighted average velocity component along the axis:
  double w_tot = 0;
  double v_avg = 0;
  for(int i = 0; i < ents.size(); i++) {
    double w_curr = vd_meas_ent(m, ents.at(i));
    w_tot = w_tot + w_curr;

    double comp = 0;
    for(int j = 0; j < v_nbr; j++) {
      apf::getVector(vel_field, ent_verts.at(i).at(j), 0, v_temp);
      comp = comp + v_temp*ax;
    }
    v_avg = comp*w_curr/v_nbr;
  }
  v_avg = v_avg/w_tot;

  // Calculate the relative velocities:
  for(int i = 0; i < verts.size(); i++) {
    apf::getVector(vel_field, verts.at(i), 0, v_temp);
    v_temp = v_temp - ax*v_avg;
    if(proj_flag == (int) PROJ_TYPE::EXT_SHELL and 
                              f_calc.chk_vert_special(m, verts.at(i)) ) {
      v_temp = f_calc.get_vec_special(m, verts.at(i), v_temp);
    }
    apf::setVector(vel_field, verts.at(i), 0, v_temp);
  }

  double v_r_avg = 0;
  for(int i = 0; i < ents.size(); i++) {
    double w_curr = vd_meas_ent(m, ents.at(i));

    double comp = 0;
    for(int j = 0; j < v_nbr; j++) {
      apf::getVector(vel_field, ent_verts.at(i).at(j), 0, v_temp);
      comp = comp + v_temp*ax;
    }
    v_r_avg = comp*w_curr/v_nbr;
  }
  v_r_avg = v_r_avg/w_tot;
  // Do at least one cycle.
  double err = v_avg/v_avg;

  double mult_start = f_calc.find_min_t2(m, &verts)/20;
  double mult = mult_start;
  // Update the relative velocities and evolve:
  while(err > th and mult > mult_start/1000) {
    // Move:
    double mult = f_calc.find_min_t2(m, &verts)/20;
    vd_apply_vel_field(m, &verts, mult);

    // Update the velocity fields:
    f_calc.vd_upd_vel_field(m, &verts);

    // Calculate the weighted average velocity component along the axis:
    w_tot = 0;
    v_avg = 0;
    for(int i = 0; i < ents.size(); i++) {
      double w_curr = vd_meas_ent(m, ents.at(i));
      w_tot = w_tot + w_curr;

      double comp = 0;
      for(int j = 0; j < v_nbr; j++) {
        apf::getVector(vel_field, ent_verts.at(i).at(j), 0, v_temp);
        comp = comp + v_temp*ax;
      }
      v_avg = comp*w_curr/v_nbr;
    }
    v_avg = v_avg/w_tot;

    // Calculate the relative velocities:
    for(int i = 0; i < verts.size(); i++) {
      apf::getVector(vel_field, verts.at(i), 0, v_temp);
      v_temp = v_temp - ax*v_avg;
      if(proj_flag == (int) PROJ_TYPE::EXT_SHELL and 
                                f_calc.chk_vert_special(m, verts.at(i)) ) {
        v_temp = f_calc.get_vec_special(m, verts.at(i), v_temp);
      }
      apf::setVector(vel_field, verts.at(i), 0, v_temp);
    }

    v_r_avg = 0;
    for(int i = 0; i < ents.size(); i++) {
      double w_curr = vd_meas_ent(m, ents.at(i));

      double comp = 0;
      for(int j = 0; j < v_nbr; j++) {
        apf::getVector(vel_field, ent_verts.at(i).at(j), 0, v_temp);
        comp = comp + v_temp*ax;
      }
      v_r_avg = comp*w_curr/v_nbr;
    }
    v_r_avg = v_r_avg/w_tot;
    // Do at least one cycle.
    double err = std::fabs((v_avg - v_r_avg)/v_avg);
    std::cout << "v_avg: " << v_avg
              << " v_r_avg: " << v_r_avg
              << " err: " << err
              << " t: " << mult
              << std::endl;

    save_vtk_name("./output/evolve");

    if(mult < mult_start/100) {
      adapt();
      //e_list.refresh();
      // Collect the entities and their downward vertices:
      e_nbr = 0;
      for(int i = 0; i < cells->size(); i++) {
        int id = cells->at(i).second - 1;
        std::vector<apf::MeshEntity*> * temp = &e_list.e.at(dim).at(id).at(dim);
        e_nbr = e_nbr + temp->size();
      }
      ents.clear();
      for(int i = 0; i < ent_verts.size(); i++)
        ent_verts.at(i).clear();
      ent_verts.clear();

      ents.reserve(e_nbr);
      ent_verts.reserve(e_nbr);
      for(int i = 0; i < cells->size(); i++) {
        int id = cells->at(i).second - 1;
        std::vector<apf::MeshEntity*> * temp = &e_list.e.at(dim).at(id).at(dim);
        for(int j = 0; j < temp->size(); j++) {
          vd_set_down(m, temp->at(j), &temp_vert, dim);
          ents.push_back(temp->at(j));
          ent_verts.push_back(temp_vert);
        }
      }
      vd_set_down(m, &ents, &verts, dim);

    }
  }
}

// Given a 1- or 2-cell, and a position, calculate the weighted average 
// radial velocity. Evolve the vertices of the boundary cells until
// the weighted average of the relative velocity is below a threshold of the  
// weighted average.
void vd_sim::evolve_bound_rad(std::pair<int, int> * cell, 
                          apf::Vector3 pos, double th) {
}

// Given a same dim list of 0-,1-,2-cells and an axis, evolve the boundary 
// vertices until the component of the velocity in the axis direction for all
// vertices for subsequent steps falls below the threshold.
void vd_sim::evolve_bound_ax(std::vector<std::pair<int, int> > * cells, 
                          apf::Vector3 ax, double th) {
  f_calc.vd_att_fields(m);

  if(th > 1 or th < - std::numeric_limits<double>::min())
    th = 1;

  if(cells->size() == 0)
    return;
  int dim = cells->at(0).first;
  int v_nbr = dim + 1;
  for(int i = 1; i < cells->size(); i++) {
    if(dim != cells->at(i).first)
      return;
  }
  e_list.refresh();

  apf::Field* vel_field = m->findField("velocity_field");

  apf::Vector3 v_temp(0,0,0);
  std::vector<apf::MeshEntity*> verts(0);
  std::vector<double> vert_vel(0);
  std::vector<apf::MeshEntity*> ents(0);
  std::vector<std::vector<apf::MeshEntity*> > ent_verts(0, 
                              std::vector<apf::MeshEntity*>(0));
  std::vector<apf::MeshEntity*> temp_vert(0);
  // Collect the entities and their downward vertices:
  int e_nbr = 0;
  for(int i = 0; i < cells->size(); i++) {
    int id = cells->at(i).second - 1;
    std::vector<apf::MeshEntity*> * temp = &e_list.e.at(dim).at(id).at(dim);
    e_nbr = e_nbr + temp->size();
  }

  ents.reserve(e_nbr);
  ent_verts.reserve(e_nbr);
  for(int i = 0; i < cells->size(); i++) {
    int id = cells->at(i).second - 1;
    std::vector<apf::MeshEntity*> * temp = &e_list.e.at(dim).at(id).at(dim);
    for(int j = 0; j < temp->size(); j++) {
      vd_set_down(m, temp->at(j), &temp_vert, dim);
      ents.push_back(temp->at(j));
      ent_verts.push_back(temp_vert);
    }
  }
  vd_set_down(m, &ents, &verts, dim);
  vert_vel.resize(verts.size());

  // Update the velocities:
  f_calc.vd_upd_vel_field(m, &verts);

  // Calculate the relative velocities:
  for(int i = 0; i < verts.size(); i++) {
    apf::getVector(vel_field, verts.at(i), 0, v_temp);
    vert_vel.at(i) = ax*v_temp;
  }

  double mult_start = f_calc.find_min_t2(m, &verts)/20;
  double mult = mult_start;
  // Update the relative velocities and evolve:
  double rat_max = 1;

  while(rat_max > th and mult > mult_start/1000) {
    // Move:
    double mult = f_calc.find_min_t2(m, &verts)/20;
    vd_apply_vel_field(m, &verts, mult);

    // Update the velocity fields:
    f_calc.vd_upd_vel_field(m, &verts);
    // Calculate the relative velocities:
    rat_max = 0;
    for(int i = 0; i < verts.size(); i++) {
      apf::getVector(vel_field, verts.at(i), 0, v_temp);
      double proj = ax*v_temp;
      double rat_curr;
      if(std::fabs(proj) > std::numeric_limits<double>::min())
        rat_curr = std::fabs((proj - vert_vel.at(i))/vert_vel.at(i));
      else
        rat_curr = 0;
      vert_vel.at(i) = proj;
      if(rat_curr > rat_max)
        rat_max = rat_curr;
    }

    // Do at least one cycle.
    std::cout << "rat_max: " << rat_max
              << " t: " << mult
              << std::endl;

    save_vtk_name("./output/evolve");

    if(mult < mult_start/100) {
      adapt();
      //e_list.refresh();
      // Collect the entities and their downward vertices:
      e_nbr = 0;
      for(int i = 0; i < cells->size(); i++) {
        int id = cells->at(i).second - 1;
        std::vector<apf::MeshEntity*> * temp = &e_list.e.at(dim).at(id).at(dim);
        e_nbr = e_nbr + temp->size();
      }
      ents.clear();
      for(int i = 0; i < ent_verts.size(); i++)
        ent_verts.at(i).clear();
      ent_verts.clear();

      ents.reserve(e_nbr);
      ent_verts.reserve(e_nbr);
      for(int i = 0; i < cells->size(); i++) {
        int id = cells->at(i).second - 1;
        std::vector<apf::MeshEntity*> * temp = &e_list.e.at(dim).at(id).at(dim);
        for(int j = 0; j < temp->size(); j++) {
          vd_set_down(m, temp->at(j), &temp_vert, dim);
          ents.push_back(temp->at(j));
          ent_verts.push_back(temp_vert);
        }
      }
      vd_set_down(m, &ents, &verts, dim);
      for(int i = 0; i < verts.size(); i++) {
        apf::getVector(vel_field, verts.at(i), 0, v_temp);
        double proj = ax*v_temp;
        double rat_curr;
        if(std::fabs(proj) > std::numeric_limits<double>::min())
          rat_curr = std::fabs((proj - vert_vel.at(i))/vert_vel.at(i));
        else
          rat_curr = 0;
        vert_vel.at(i) = proj;
      }
      f_calc.vd_upd_vel_field(m, &verts);
    }
  }
}

// Given a 1- or 2-cell, and a position, calculate the weighted average 
// radial velocity. Evolve the vertices of the boundary cells until
// the weighted average of the relative velocity is below a threshold of the  
// weighted average.
void vd_sim::evolve_bound_rad_conv(std::pair<int, int> * cell, 
                          apf::Vector3 pos, double th) {
}


/*
double vd_sim::shr_cell(int cell_dim, int cell_id, double sp_size) {

  assert(!check_ma_sgn());
  assert(sp_size < 1);

  std::vector<apf::MeshEntity*> vert_b(0);
  std::vector<apf::MeshEntity*> vert(0);
  std::vector<apf::MeshEntity*>* es_cell = 
                        &e_list.e.at(cell_dim).at(cell_id-1).at(cell_dim);

  vd_set_down(m, es_cell, &vert_b, cell_dim);
  vert.reserve(vert.size());
  if(c_base->get_cell_bound_ext_gmi(cell_dim, cell_id)) {
    for (int i = 0; i < vert_b.size(); i ++) {
      apf::ModelEntity* mdl = m->toModel(vert_b.at(i));
      int c_dim = m->getModelType(mdl);
      int c_id = m->getModelTag(mdl);
      if(c_dim < 3 and c_base->get_cell_ext_gmi(c_dim, c_id))
        vert.push_back(vert_b.at(i));
    }
  }
  else
    vert = vert_b;

  if(proj_flag == (int)PROJ_TYPE::EXT_SHELL) {
    for (int i = 0; i < vert_b.size(); i ++) {
      f_calc.corr_pos(m, vert_b.at(i));
    }
  }

  apf::Vector3 cell_center = vd_get_center(m, &vert);

  std::vector<apf::Vector3> point(vert_b.size(), apf::Vector3(0,0,0));
  std::vector<apf::Vector3> disp(vert_b.size(), apf::Vector3(0,0,0));

  apf::Field* vec_field = vd_att_vv_field(m, "vec_field");

  // Displace the domain vertices:
  for (int i = 0; i < vert_b.size(); i ++) {
    m->getPoint(vert_b.at(i), 0, point.at(i));    
    disp.at(i) = (cell_center-point.at(i));
    if (f_calc.chk_vert_special(m, vert_b.at(i))) {
      disp.at(i) = f_calc.get_vec_special(m, vert_b.at(i), disp.at(i));
    }
    apf::setVector(vec_field, vert_b.at(i), 0, disp.at(i));
  }

  double mult = vd_find_min_mult(m, &vert_b, "vec_field");

  mult = std::min(mult/5, sp_size);

  // Displace the domain vertices:
  for (int i = 0; i < vert_b.size(); i ++) {
    m->setPoint(vert_b.at(i), 0, point.at(i) + disp.at(i)*mult);
  }
  m->removeField(vec_field);

  assert(!check_ma_sgn(&vert_b));

  e_list.refresh();

  return sp_size;

}
*/

// Check if the cell is still shrinking after shrinking slightly. 
// TODO, a more stable time iteration scheme could be used instead. 
bool vd_sim::chk_cell_shr(int cell_dim, int cell_id) {
  // TODO this check could be more efficient if it considers only vertices 
  // associated with the stratum in question.
  assert(!check_ma_sgn());
  double sp_size = 0.1;
  //assert(sp_size < 1);

  std::vector<apf::MeshEntity*> vert_b(0);
  std::vector<apf::MeshEntity*> vert(0);
  std::vector<apf::MeshEntity*>* es_cell = 
                        &e_list.e.at(cell_dim).at(cell_id-1).at(cell_dim);

  // The collapsing triangles and tets.
  std::vector<apf::MeshEntity*> tri_col(0);
  std::vector<apf::MeshEntity*> tet_col(0);

  std::map<apf::MeshEntity*, bool> tet_skip {};
  vd_set_up(m, &e_list.e.at(cell_dim).at(cell_id-1).at(1), &tri_col);
  vd_set_up(m, &tri_col, &tet_col);

  for(int i = 0; i < tet_col.size(); i++) {
    tet_skip[tet_col.at(i)] = true;
  }

  vd_set_down(m, es_cell, &vert_b, cell_dim);
  vert.reserve(vert.size());
  if(c_base->get_cell_bound_ext_gmi(cell_dim, cell_id)) {
    for (int i = 0; i < vert_b.size(); i ++) {
      apf::ModelEntity* mdl = m->toModel(vert_b.at(i));
      int c_dim = m->getModelType(mdl);
      int c_id = m->getModelTag(mdl);
      if(c_dim < 3 and c_base->get_cell_ext_gmi(c_dim, c_id))
        vert.push_back(vert_b.at(i));
    }
  }
  else
    vert = vert_b;

  if(proj_flag == (int)PROJ_TYPE::EXT_SHELL) {
    for (int i = 0; i < vert_b.size(); i ++) {
      f_calc.corr_pos(m, vert_b.at(i));
    }
  }

  apf::Vector3 cell_center = vd_get_center(m, &vert);

  std::vector<apf::Vector3> point(vert_b.size(), apf::Vector3(0,0,0));
  std::vector<apf::Vector3> disp(vert_b.size(), apf::Vector3(0,0,0));

  apf::Field* vec_field = vd_att_vv_field(m, "vec_field");

  // Displace the domain vertices:
  for (int i = 0; i < vert_b.size(); i ++) {
    m->getPoint(vert_b.at(i), 0, point.at(i));    
    disp.at(i) = (cell_center-point.at(i));
    if (f_calc.chk_vert_special(m, vert_b.at(i))) {
      disp.at(i) = f_calc.get_vec_special(m, vert_b.at(i), disp.at(i));
    }
    apf::setVector(vec_field, vert_b.at(i), 0, disp.at(i));
  }

  double mult = vd_find_min_mult(m, &vert_b, tet_skip, "vec_field")/2;
  mult = std::min(mult, 1. - 0.5*sp_size);
  std::cout << "Shrinking by " << mult << std::endl;
  //mult = std::min(mult/2, sp_size);
  // Force the multiplier to be less than sp_size. If not, continue to evolve.
  // If the cell is shrinking, it will become smaller and smaller until it will
  // become more easy for it to be shrunk and possible for it to collapse.
  if((1-mult) * sm_map.at(cell_dim-1)[cell_id] < len_trans/ratio_col_sh) {
  //if((1-mult) < sp_size) {
    // Displace the domain vertices:
    for (int i = 0; i < vert_b.size(); i ++) {
      m->setPoint(vert_b.at(i), 0, point.at(i) + disp.at(i)*mult);
    }
    m->removeField(vec_field);
    for (int i = 0; i < vert_b.size(); i ++) {
      f_calc.vd_upd_vel_field(m, vert_b.at(i));
    }

    double roc = f_calc.calc_roc(m, 
                 &e_list.e.at(cell_dim).at(cell_id-1).at(cell_dim));
    std::cout << "After shrinking: Rate of change " << roc << std::endl;

    // Displace the domain vertices:
    for (int i = 0; i < vert_b.size(); i ++) {
      m->setPoint(vert_b.at(i), 0, point.at(i));
    }
    for (int i = 0; i < vert_b.size(); i ++) {
      f_calc.vd_upd_vel_field(m, vert_b.at(i));
    }

    std::cout << "Shr mult " << mult << " roc after shr " << roc << std::endl;
    return (roc < -std::numeric_limits<double>::min());
  }
  else {
    std::cout << "Above collapse threshold." << std::endl;
    m->removeField(vec_field);
    return false;
  }
}

void vd_sim::set_glens_param(vd_glens& g_lens) {
  g_lens.set_field_calc(&f_calc);
  g_lens.turn_save_vtk(sub_vtk);
}

void vd_sim::set_glens_trial_param(vd_glens_trial& g_trial, bool precond_flag) {
  g_trial.set_precond(precond_flag);
  g_trial.set_verify(false);
}

// Collapse a cell given in gmi notation.
std::pair<int, int> vd_sim::col_cell(int cell_dim, int cell_id) {
  // First, get the list of the 0cells that require their path and circuit lists
  // updated, after the collapse.
  std::pair<int, int> tag_cell = std::make_pair(-1, -1);

  if(!(c_base->is_free(cell_dim, cell_id-1)) ) {
    if(load_flag) {
      struct ent_conn e_0cell;
      c_base->get_conn_dim_gmi(0, cell_dim, cell_id, &e_0cell);
      std::cout << "Merging 0cells, ";
      for (int i = 0; i < e_0cell.conn.size(); i++) {
        std::cout << e_0cell.conn.at(i) << ", ";
      }
      std::cout << std::endl;

/*      
      apf::Field* vel_field = m->findField("velocity_field");
      apf::destroyField(vel_field);

      if(ad_type == ADAPT_TYPE::ADAPT_STEP_1CELL) {
        adapt_col(cell_dim, cell_id);
      }
*/
      // Try to coarsen the neighborhood around the collapsing cell and try
      // collapsing again. If not skip.
      e_list.refresh();
      f_calc.vd_del_fields(m);
      adapt_coarsen_cell(cell_dim, cell_id);
      //e_list.refresh();
      f_calc.vd_att_fields(m);

      //If the chk_cell_thresh is satisfied, no need to shrink...
      //double shr_sz = 0.2;
      //for(int i = 0; i < 2; i++) {
      //  shr_sz = shr_cell(cell_dim, cell_id, shr_sz);
      //  if(shr_sz < 0.05)
      //    i = 2;
      //}
      std::vector<apf::MeshEntity*> tets(0);

      vd_get_tet(m, e_list.e.at(cell_dim).at(cell_id-1).at(cell_dim), tets);

      std::vector<std::vector<apf::MeshEntity*> > ent_col(3-cell_dim + 1, 
                                          std::vector<apf::MeshEntity*>(0));
      ent_col.at(0) = e_list.e.at(cell_dim).at(cell_id-1).at(cell_dim);
      for(int dim = 0; dim < 3 - cell_dim; dim++) {
        vd_set_up(m, &ent_col.at(dim), &ent_col.at(dim+1));
      }
      vd_remove_set(&tets, &ent_col.back());

      apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
      adapt_mark_0c_min();
      Step_ns sf(m, get_adapt_ln());
      ma::Input* in = ma::configure(m, &sf);

      bool high_q = true;
      for(int i = 0; i < tets.size(); i++) {
        double q_temp = measureTetQuality(m, in->sizeField, tets.at(i));
        if(q_temp < q_th*q_th) {
          high_q = false;
          i = tets.size();
        }
      }
      if(!high_q) {
        rad_map_cc.at(cell_dim-1)[cell_id] = true;
        std::cout << "Low quality tets... " << std::endl;
      }
      delete in;
      apf::destroyField(field_step);

      vd_glens g_lens2(m, c_base, &e_list);
      set_glens_param(g_lens2);
      bool precond_flag = g_lens2.chk_precond(cell_dim, cell_id);
      if(high_q or !precond_flag) {

        //if(!g_lens2.test_cvx(cell_dim, cell_id)) {
        //  rad_map_cc.at(cell_dim-1)[cell_id] = true;
        //  std::cout << "Cannot collapse, no cvx is found. " << std::endl;
        //}
        //else {
          std::cout << "Trying to collapse " << cell_dim << "c" << cell_id << std::endl;

          std::vector<ent_conn> e_low(cell_dim - 1, ent_conn());

          for(int dim_low = 1; dim_low < cell_dim; dim_low++) {
            c_base->get_conn_dim_gmi(dim_low, cell_dim, cell_id, 
                                                      &e_low.at(dim_low-1));
          }

          vd_glens_trial g_trial(m, c_base, &e_list, &f_calc);
          set_glens_trial_param(g_trial, precond_flag);
          tag_cell = g_trial.coll_cell(cell_dim, cell_id);

          //tag_cell = g_lens2.col_cell(cell_dim, cell_id);
          // Currently preconditioning could fail if relaxation takes too long.
          // In that case, the collapse is skipped for this iteration.
          bool skipped = false;
          if(tag_cell.first == -2) {
            rad_map_cc.at(cell_dim - 1)[cell_id] = true;

            ent_conn* e_l = new ent_conn();
            for(int dim_low = cell_dim - 1; dim_low > 0; dim_low--) {
              c_base->get_conn_dim_gmi(dim_low, cell_dim, cell_id, e_l);
              for(int i = 0; i < e_l->conn.size(); i++) {
                rad_map_cc.at(dim_low - 1)[e_l->conn.at(i)] = true;
              }
            }
            delete e_l;
            skipped = true;
          }
          if(!skipped) {
            //c_base->print_ent();

            // TODO when doing things in parallel, the correct thing is to 
            // apply the topological modifications at this level.
            f_calc.vd_att_fields(m);

            for (int i = 0; i < e_0cell.conn.size(); i++) {
              std::cout << "Updating paths 0cell" << e_0cell.conn.at(i) << std::endl;
              //c_ins.update_path_gmi(e_0cell.conn.at(i));
            }
            for (int i = 0; i < e_0cell.conn.size(); i++) {
              ins_rad_map[e_0cell.conn.at(i)] = -1;
            }
            for (int i = 0; i < e_0cell.conn.size(); i++)
              tag_0cell_stable[e_0cell.conn.at(i)] = false;

            m->acceptChanges();
            reload_mesh();
            m->verify();

            e_list.refresh();
            if(tag_cell.first == 0) {

              ent_conn* e_l = new ent_conn();
              for(int dim_high = 1; dim_high < 4; dim_high++) {
                c_base->get_conn_dim_gmi(dim_high, 0, tag_cell.second, e_l);
                for(int i = 0; i < e_l->conn.size(); i++) {
                  tag_shell_chk.at(dim_high-1)[e_l->conn.at(i)] = false;
                  tag_shell_map.at(dim_high-1)[e_l->conn.at(i)] = false;
                }
              }
              delete e_l;

              assert(e_list.e.at(0).at(tag_cell.second - 1).at(0).size() == 1);
              apf::MeshEntity* v_ctr = e_list.e.at(0).at(tag_cell.second - 1).at(0).at(0);

              std::vector<apf::MeshEntity*> ee(0);
              std::vector<apf::MeshEntity*> ev(0);

              vd_set_up(m, v_ctr, &ee);
              ev.resize(ee.size());

              for(int i = 0; i < ee.size(); i++) {
                ev.at(i) = apf::getEdgeVertOppositeVert(m, ee.at(i), v_ctr);
              }

              apf::Vector3 temp(0,0,0);
              for(int step = 0; step < 10; step++) {
                // Displace the domain vertices:
                for (int j = 0; j < ev.size(); j++) {
                  f_calc.vd_upd_vel_field(m, ev.at(j));
                }
                double mult = f_calc.find_min_t2(m, &ev)/20;
                vd_apply_vel_field(m, &ev, mult);
              }

              // Replacement of spurious cells after collapse:
              if(c_base->fix_spur(tag_cell.first, tag_cell.second-1, true)) {
                c_base->process_fix_list();

                sort_celem(&c_base->fix_list);

                std::cout << "To be replaced " << std::endl;
                for(int i = 0; i < c_base->fix_list.size(); i++) {
                  std::cout << c_base->fix_list.at(i).first.first << "c"
                            << c_base->fix_list.at(i).first.second << " "
                            << c_base->fix_list.at(i).second.first << "c"
                            << c_base->fix_list.at(i).second.second << std::endl;
                }
                vd_cell_repl v_repl(m);
                //v_repl.set_cell(&fix_new);
                v_repl.set_cell(&c_base->fix_list);
                v_repl.repl_cell();
                std::cout << "Replaced " << std::endl;
                //c_base->print_ent();
              }

              save_vtk_name("output/after_replace");

            }
            t_below_map.at(cell_dim - 1)[cell_id] = -1;
            for(int dim_low = 1; dim_low < cell_dim; dim_low++) {
              for(int i = 0; i < e_low.at(dim_low - 1).conn.size(); i++)
                t_below_map.at(dim_low - 1)[e_low.at(dim_low - 1).conn.at(i)] = -1;
            }

          }
        //}
      }
    }
  }

  return tag_cell;
}


// Given a set of tetrahedra belonging to a disjoint 3stratum, collect the 
// bounding 2stratum triangles touching the central vertex. Calculate the 
// anglecosine weighted sum of area normals. Also, return false if the anglecosine
// between the direction and one of plane normals is lower than a threshold.
bool vd_sim::coll_bound_tri_dir(apf::MeshEntity* vert, 
                                  std::vector<apf::MeshEntity*> &es_tet,
                                  apf::Vector3 &dir, apf::Vector3 pos, 
                                  double ang_th) {

  apf::ModelEntity* mdl_3c = m->toModel(es_tet.at(0));
  std::map<apf::MeshEntity*, apf::MeshEntity*> tri_2_tet_map{};
  std::map<apf::MeshEntity*, apf::Vector3> tet_dir_map{};
  apf::Up up;

  apf::Vector3 nn(0,0,0);
  apf::Vector3 mm(0,0,0);

  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_tri(0);
  apf::Downward d_v;

  vd_set_down(m, &es_tet, &es_surf);
  es_tri.reserve(es_surf.size());

  for(int i = 0; i < es_tet.size(); i++) {
    tet_dir_map[es_tet.at(i)] = vd_get_pos(m, es_tet.at(i));
  }

  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* mdl = m->toModel(es_surf.at(i));
    int m_type = m->getModelType(mdl);
    int m_tag = m->getModelTag(mdl);

    m->getDownward(es_surf.at(i), 0, d_v);
    int j_v = findIn(d_v, 3, vert);
    if(m_type == 2 and j_v > -1) {
      es_tri.push_back(es_surf.at(i));
      //std::cout << e_lens.es_surf.at(i) 
      //          << m_type << "c" << m_tag << std::endl;
      m->getUp(es_surf.at(i), up);
      if(m->toModel(up.e[0]) == mdl_3c) {
        assert(findIn(&es_tet, es_tet.size(), up.e[0]) > -1);
        tri_2_tet_map[es_surf.at(i)] = up.e[0];
        tet_dir_map[es_surf.at(i)] = tet_dir_map[up.e[0]] 
                                    - vd_get_pos(m, es_surf.at(i));
      }
      else {
        assert(findIn(&es_tet, es_tet.size(), up.e[1]) > -1);
        assert(up.n == 2 and m->toModel(up.e[1]) == mdl_3c);
        tri_2_tet_map[es_surf.at(i)] = up.e[1];
        tet_dir_map[es_surf.at(i)] = tet_dir_map[up.e[1]] 
                                    - vd_get_pos(m, es_surf.at(i));
      }
    }
  }
  vd_disc disc_bound;
  disc_bound.m = m;
  disc_bound.v_ctr = vert;
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
    rel_dir1 = disc_bound.v_pos.at(e_1) - pos;
    rel_dir2 = disc_bound.v_pos.at(e_2) - pos;
    rel_dir1 = norm_0(rel_dir1);
    rel_dir2 = norm_0(rel_dir2);
    mm = cross(rel_dir1, rel_dir2);
    if(mm * tet_dir_map[tri] < - std::numeric_limits<double>::min())
      mm = mm * (-1);

    double acos = std::acos(rel_dir1*rel_dir2);
    nn = nn + mm*acos*0.5;
  }

  // Calc m:
  mm = apf::Vector3(0,0,0);
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
      rel_dir0 = disc_bound.v_pos.at(i) - pos;
      rel_dir1 = disc_bound.v_pos.at(e_1) - pos;
      rel_dir2 = disc_bound.v_pos.at(e_2) - pos;

      rel_dir0 = norm_0(rel_dir0);
      rel_dir1 = norm_0(rel_dir1);
      rel_dir2 = norm_0(rel_dir2);

      double acos = std::acos(rel_dir1*rel_dir0) + std::acos(rel_dir2*rel_dir0);
      mm = mm + rel_dir0*acos*0.5;
    }
  }
  std::cout << "n = [" << nn[0] << ", " << nn[1] << ", " << nn[2] << "]; "
            << " m = [" << mm[0] << ", " << mm[1] << ", " << mm[2] << "]; "
            << std::endl;
  if(nn.getLength() > mm.getLength())
    dir = norm_0(nn);
  else
    dir = norm_0(mm);
  assert(!std::isnan(dir.getLength()));

  if(es_tet.size() == 1) {
    if((vd_get_pos(m, es_tet.at(0)) - pos)*dir < 
                            std::numeric_limits<double>::min())
      return false;

  }
  for(int i = 0; i < disc_bound.tri.size(); i++) {
    apf::MeshEntity* tri = disc_bound.tri.at(i);
    int e_1;
    int e_2;

    e_1 = disc_bound.t_e_map1[tri] - 1;
    e_2 = disc_bound.t_e_map2[tri] - 1;

    // Relative directions of the edges preceeding and succeeding the current 
    // edge on the disc.
    rel_dir1 = disc_bound.v_pos.at(e_1) - pos;
    rel_dir2 = disc_bound.v_pos.at(e_2) - pos;
    rel_dir1 = norm_0(rel_dir1);
    rel_dir2 = norm_0(rel_dir2);
    mm = cross(rel_dir1, rel_dir2);
    if(mm * tet_dir_map[tri] <- std::numeric_limits<double>::min())
      mm = mm * (-1);
    if(mm*dir > ang_th)
      return false;
  }
  return true;
}

// Given a 0cell vertex, find the minimum of the distances towards the disjoint
// set of entities belonging to a 3cell around that vertex. 
double vd_sim::dist_3c_min(apf::MeshEntity* vert) {
  int lookup_ts [4] = {2, 3, 1, 0};
  double len_min = sim_len;

  vd_3c_det* vd_3c = new vd_3c_det();
  vd_3c->reload(m, c_base, vert);
  int c3_nbr = vd_3c->get_3c_nbr();

  apf::Vector3 vb(0,0,0);
  apf::Vector3 pos(0,0,0);
  m->getPoint(vert, 0, pos);

  for (int i = 0; i < c3_nbr; i++) {

    std::vector<apf::MeshEntity*> tet_3c(0);
    tet_3c = vd_3c->get_tet_id(i);
    assert(tet_3c.size() > 0);
    // Better alternative to vd_find_dist_w, consider the intersections to the 
    // triangles across.
    std::vector<apf::MeshEntity*> tri(0);
    tri.reserve(tet_3c.size());
    apf::Downward d_t;
    apf::Downward d_v;
    int lookup_ts [4] = {2, 3, 1, 0};

    for(int i = 0; i < tet_3c.size(); i++) {
      m->getDownward(tet_3c.at(i), 0, d_v);
      m->getDownward(tet_3c.at(i), 2, d_t);

      int i1 = findIn(d_v, 4, vert);
      assert(i1 > -1);

      tri.push_back(d_t[lookup_ts[i1]]);
    }

    if(coll_bound_tri_dir(vert, tet_3c, vb, pos)) {
      double len_curr = vd_dist_v_x_pl_dir(m, vert, tri, vb);
      if(len_curr < len_min)
        len_min = len_curr;
    }
    else {
      std::cout << "The disc of the tetrahedra are ill-posed." << std::endl;
      len_min = -1;
    }
  }

  delete vd_3c;
  return len_min;
}


// Try inserting cells around the given 0cell. For the isentropic case, it also
// checks for the Euler characteristic.
bool vd_sim::ins_cell(int tag_0cell) {
  // TODO this check could be more efficient if it only considers the vertex.
  assert(!check_ma_sgn());

  if(!ins_flag)
    return false;
  // TODO update when including exterior motion.
  if(!calc_ext and c_base->chk_0cell_ext_gmi(tag_0cell))
    return false;
  set_free_cells();

  assert(load_flag);

  bool ins_res = false;
  vd_edisc* e_d = new vd_edisc(get_mesh(), get_c_base(), &e_list);
  set_edisc_param(e_d);

  //get_c_ins()->print_graph();

  int c0_sz = c_base->get_sz(0);
  int c1_sz = c_base->get_sz(1);
  int c2_sz = c_base->get_sz(2);
  int c3_sz = c_base->get_sz(3);

  bool prob_2c = c_base->chk_2c_gmi(tag_0cell);

  std::cout << "0c" << tag_0cell <<  "Problematic " << prob_2c << std::endl;

  e_list.refresh();
  std::vector<apf::MeshEntity*>* es_vert_ctr = 
                              &e_list.e.at(0).at(tag_0cell-1).at(0);
  //vd_print_ent(m);
  if (es_vert_ctr->size() == 1) {
    std::vector<apf::MeshEntity*> es_e_ctr(0);
    std::vector<apf::MeshEntity*> es_s_ctr(0);
    std::vector<apf::MeshEntity*> es_ee_ctr(0);
    apf::Vector3 midpoint(0, 0, 0);
    m->getPoint(es_vert_ctr->at(0), 0, midpoint);

    vd_set_up(m, es_vert_ctr, &es_e_ctr);
    vd_set_up(m, &es_e_ctr, &es_s_ctr);
    vd_set_up(m, &es_s_ctr, &es_ee_ctr);

    std::vector<apf::MeshEntity*> tri(0);
    tri.reserve(es_ee_ctr.size());
    apf::Downward d_t;
    apf::Downward d_v;
    int lookup_ts [4] = {2, 3, 1, 0};

    for(int i = 0; i < es_ee_ctr.size(); i++) {
      m->getDownward(es_ee_ctr.at(i), 0, d_v);
      m->getDownward(es_ee_ctr.at(i), 2, d_t);

      int i1 = findIn(d_v, 4, es_vert_ctr->at(0));
      assert(i1 > -1);

      tri.push_back(d_t[lookup_ts[i1]]);
    }

    double len_main = 0;
    double len_max = 0;
    std::tie(len_main, len_max) = tri_dist_sphere(m, midpoint, &tri);
    // Find the minimum distances to the triangle planes across for each disjoint
    // set of 3cell tets along the direction obtained using the 2cell triangles
    // bounding the disjoint set of tets.
    double len_min = dist_3c_min(es_vert_ctr->at(0));

    double ins_len = len_trans/ratio_ins;
    std::cout << "Precond sphere " << len_main 
              << " l_3c_min " << len_min
              << " ins_rad " << ins_rad_map[tag_0cell]
              << " topo len ins " << ins_len << std::endl;
    if(len_main < min_ins_rad)
      min_ins_rad = len_main;
    // TODO The insertion used to require a neighborhood around the 0c vertex
    // with a preconditioning sphere with a radius larger than the topological
    // length scale, or a ratio of it. This is not always achievable and causes
    // the simulation to get stuck. To improve, removed the restriction.
    // This should be revisited in the future, as it is closely related to the
    // mesh quality and time step length and may contribute to ping pong like
    // behavior. 

    //double ins_len = len_main/1.5;
    e_d->set_len_sh(ins_len);
    if (len_main < ins_rad_map[tag_0cell] and 
        len_main < ins_len/2)
      std::cout << "Precond sphere is reducing and small!" << std::endl;

    if(ins_rad_map[tag_0cell] < -std::numeric_limits<double>::min())
      ins_rad_map[tag_0cell] = len_min;
    else
      ins_rad_map[tag_0cell] = (len_min + ins_rad_map[tag_0cell])/2;
    //ins_rad_map[tag_0cell] = len_main;

    // Smallest radius of the spheres that can be fit inside the triangular hull
    // len_main
    // Is smaller than a multiple of the longest edge insertible inside disjoint 
    // sets of 3cell tetrahedra
    // len_min.
    // TODO which doesn't make very much sense? Instead try len_min > ins_len*2
    // which requires the 3-cell edges inserted to be longer than the topological 
    // length scale.
    // And len_main larger than topological length scale or 
    // Precond sphere is reducing and small

    // if (len_main < len_min*10 and
    //if (len_min > ins_len*1.1 and
    if (
        (len_min > ins_len*1.1 or (len_min < ins_rad_map[tag_0cell] and 
        len_min < ins_len/2)) and
        //(len_main > ins_len*1.1 or (len_main < ins_rad_map[tag_0cell] and 
        //len_main < ins_len/2)) and
        //len_main < ins_len/2 and ins_len/4 < len_main)) and
    //if (
        insertible(tag_0cell) and e_d->set_0cell(tag_0cell, prob_2c) ) {
      //vd_print_ent(m);

      std::pair<int, int> new_cell = e_d->try_insert();
      if (new_cell.first > -1) {
        if(proj_flag == (int)PROJ_TYPE::EXT_SHELL) {
          f_calc.corr_pos(m);
          ent_conn e_down;
          for(int dim_l = new_cell.first-1; dim_l > 0; dim_l--) {
            c_base->get_conn_dim_gmi(dim_l, new_cell.first, new_cell.second, 
                                                                        &e_down);
            for(int i = 0; i < e_down.conn.size(); i++) {
              tag_shell_chk.at(dim_l-1)[e_down.conn.at(i)] = false;
              tag_shell_map.at(dim_l-1)[e_down.conn.at(i)] = false;
            }

          }
          for(int dim_h = new_cell.first + 1; dim_h < 4; dim_h++) {
            c_base->get_conn_dim_gmi(dim_h, new_cell.first, new_cell.second, 
                                                                        &e_down);
            for(int i = 0; i < e_down.conn.size(); i++) {
              tag_shell_chk.at(dim_h-1)[e_down.conn.at(i)] = false;
              tag_shell_map.at(dim_h-1)[e_down.conn.at(i)] = false;
            }
          }
          tag_shell_chk.at(new_cell.first-1)[new_cell.second] = false;
          tag_shell_map.at(new_cell.first-1)[new_cell.second] = false;
        }

        //std::cout << "Check c_ins circ sz" << c_ins.get_circ_sz() << std::endl;

/*
        std::cout << c0_sz << " " << c_base->get_sz(0) << " "
                  << c1_sz << " " << c_base->get_sz(1) << " "
                  << c2_sz << " " << c_base->get_sz(2) << " "
                  << c3_sz << " " << c_base->get_sz(3) << " "
                  << std::endl;
        if(c0_sz != c_base->get_sz(0) or c1_sz != c_base->get_sz(1)
           or c2_sz != c_base->get_sz(2) or c3_sz != c_base->get_sz(3)) {
          std::cout << "Updating c_ins." << std::endl;
          c_ins.update();
        }
*/
        struct ent_conn e0;
        c_base->get_conn_dim_gmi(0, new_cell.first, new_cell.second, &e0);
        //if(sync_cins()) {
        //}
        //else {
        //  for(int i = 0; i < e0.conn.size(); i++)
        //    c_ins.update_path_gmi(e0.conn.at(i));
        //}

        m->acceptChanges();
        reload_mesh();
        m->verify();

        chk_ma_swap(m);

        //vd_print_ent(m);
        if(ad_flag)
          adapt();
        save_vtk_name("output/after_adapt");

        if(c_base->fix_spur(new_cell.first, new_cell.second-1, true)) {
          c_base->process_fix_list();

          sort_celem(&c_base->fix_list);

          std::cout << "To be replaced " << std::endl;
          for(int i = 0; i < c_base->fix_list.size(); i++) {
            std::cout << c_base->fix_list.at(i).first.first << "c"
                      << c_base->fix_list.at(i).first.second << " "
                      << c_base->fix_list.at(i).second.first << "c"
                      << c_base->fix_list.at(i).second.second << std::endl;
            if(c_base->fix_list.at(i).first.first == 0)
              e0.rem_ent(c_base->fix_list.at(i).first.second);
          }
          vd_cell_repl v_repl(m);
          //v_repl.set_cell(&fix_new);
          v_repl.set_cell(&c_base->fix_list);
          v_repl.repl_cell();
          std::cout << "Replaced " << std::endl;
          //c_base->print_ent();
        }

        for(int i = 0; i < e0.conn.size(); i++) {
          tag_0cell_ins.push_back(e0.conn.at(i));
          tag_0cell_stable[e0.conn.at(i)] = false;
        }

        for (int i = 0; i < e0.conn.size(); i++) {
          ins_rad_map[e0.conn.at(i)] = -1;
        }

        save_vtk_name("output/after_replace");

        reload_mesh();

        double r_equi = 0;
        for(int dim_low = 1; dim_low < new_cell.first; dim_low++) {
          c_base->get_conn_dim_gmi(dim_low, new_cell.first, new_cell.second, &e0);
          for(int i = 0; i < e0.conn.size(); i++) {
            chk_cell_th(dim_low, e0.conn.at(i), r_equi);
            if(r_equi < min_cell.at(dim_low - 1))
              min_cell.at(dim_low - 1) = r_equi;
          }
        }
        m->verify();

        //assert(!c_base->chk_spur(new_cell.first, new_cell.second-1));

        save_vtk_name("output/after_ins");

        //get_c_ins()->update();
        //time_curr = time_curr + 0.0001*save_count;
        f_calc.vd_att_fields(m);
        ins_res = true;
      }
      f_calc.vd_att_fields(m);
      ins_res = false;
    }
    //ins_rad_map[tag_0cell] = -1;
  }

  delete e_d;
  return ins_res;
}

// Transfer the simulation parameters to the vd_edisc object.
void vd_sim::set_edisc_param(vd_edisc* e_d) {
  e_d->set_proj((PROJ_TYPE)proj_flag);

  e_d->set_vdpar(f_calc.vdparam);
  e_d->set_field_calc(f_calc);
  e_d->set_sub_vtk_flag(sub_vtk);
  e_d->set_verbose_flag(ed_verb_flag);
}

bool vd_sim::insertible(int tag_0cell) {
/*
  if (isotropic) {
    struct ent_conn cells1;
    c_base->get_conn_lower(1, 2, &c_ins.cells2.at(tag_0cell-1), &cells1);
    c_base->rem_conn(0, tag_0cell-1, 1, &cells1);
    int n_c3 = c_ins.cells3.at(tag_0cell-1).conn.size();
    int n_c2 = c_ins.cells2.at(tag_0cell-1).conn.size();
    int n_c1 = cells1.conn.size();
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
*/
  return true;
}

// If c_ins and c_base are not synchronized, reload c_ins. Assuming each 
// change in the cell complex is through vd_sim, and updates both c_ins
// and c_base, only compare the total number of 0cells.
/*
bool vd_sim::sync_cins() {
  int c0_sz = c_base->get_sz(0);
  if(c0_sz != c_ins.get_circ_sz()) {
    std::cout << "c0_sz: " << c0_sz << ", c_ins c0_sz: "
              << c_ins.get_circ_sz() << std::endl;
    std::cout << "Updating c_ins." << std::endl;
    c_ins.load_cb(c_base, calc_corner);
    return true;
  }
  else
    assert(c_ins.get_cbase() == c_base);

  return false;
}
*/
// Default constructor:
vd_sim::vd_sim() : 
  f_calc(VEL_TYPE::MASON),  
  e_list(), e_sh(),
  //c_ins(), 
  m(NULL), c_base(NULL), 
  mesh_flag(false), load_flag(false), first_time(true),
  save_mov(false), sub_vtk(true), ins_flag(true), ed_verb_flag(false), 
  col_flag(true), extract_flag(false), 
  //outputMS("./MS.csv"),
  //outputROC("./ROC.csv"),
  //outputVOL("./VOL.csv"),
  ext_opts(""),
  isotropic(true),
  calc_ext(false), calc_corner(false),
  proj_flag((int)PROJ_TYPE::FIXED), 
  integ_flag(INTEG_TYPE::EULER),
  ext_gam_flag(EXT_GAM_TYPE::CONSTANT),
  sim_sz(1), sim_len(-1), lay_rat(4),
  t_rat(20), roc_tol(-1./15),
// cell_meas_max(1), cell_meas_min(1),
  avg_cell(3), min_cell(3), min_ins_rad(-1),
               rad_map(3, std::map<int, bool>{}), 
               rad_map_th(3, std::map<int, bool>{}), 
               rad_map_cc(3, std::map<int, bool>{}), 
               sm_map(3, std::map<int, double>{}), 
               meas_map(3, std::map<int, double>{}), 
               meas_below_map(3, std::map<int, double>{}), 
               t_below_map(3, std::map<int, double>{}), 
               ins_rad_map{},
  smb_name("./tempmesh/temp"), vtk_name("./outputiter/after"),

  adapt_len(2), ad_param(1), 
  coarse_th(0.5), split_th(1.5),

  len_f(true), len_col_f(false),
  ad_flag(true), ad_type(ADAPT_TYPE::ADAPT_1CELL),
  ad_fix_shape(true),
  cells_ref(0, std::make_pair(-1,-1)), bound_len(1), 
  correct(false), edge_min(0.03), surf_min(0.003), vol_min(0.003),
  len_topo(-1), 
  len_trans(-1), tag_0cell_stable{},
  tag_shell_chk(3, std::map<int, bool>{}), 
  tag_shell_map(3, std::map<int, bool>{}), 
  tag_0cell_ins(0),

  dt_adapt(0),
  fixed_time(false), time_curr(0), time_end(0), dt_set(0.001), t_sub(0.0),
  t_sub_limit(10e-8), limit_flag(true), 
  flag_g_nbr(false), cond_g_nbr(0),
  flag_s_nbr(false), cond_s_nbr(0),
  flag_l_nbr(false), cond_l_nbr(0),

  dt_vtk_period(-1), t_save_last(0.),
  iter_sub(10),
  ratio_ent(2.5), ratio_ins(10), ratio_col_sh(16), ratio_col(32),
  q_th(10e-05), ln_min(-1), vol_th(10e-15),
  vd_tim((int)T_COST_TYPE::END), time_cost_str("./output/t_cost"),
  cost_t_accum(true),
  f_v2move(false), c_2move(0, std::make_pair(0,0)),
  bound_manual_count(0), bound_manual_count_limit(5)
// Strict prevention to stop ping-pong-like insertions 
// r_th_min > len_trans/ratio_col
// len_trans/ratio_ins*3/rho_rat/200/4 > len_trans/ratio_col
// => ratio_col/4266 > ratio_ins
{
  //set_euler_char(true);
}

void vd_sim::set_mesh(const char* modelFile) {

  // Load the mesh. set cellbase
  if (load_flag) {
    printf("VDSIM: Mesh already loaded, destroying first.\n");
    clean_mesh();
  }
  m = apf::makeEmptyMdsMesh(gmi_load(modelFile), 3, false);
  vd_create_tess(m, modelFile);
  mesh_flag = true;

  vd_rem_tag(m);
  vd_tag_mesh(m);

  apf::writeVtkFiles("./output/meshset", m);
  m->writeNative("./tempmesh/temp_load.smb");

  if(first_time) {
    c_base = new cell_base(m->getModel());
    reload_len();
  }
  e_list.change_mesh(m, c_base);

  if(first_time and proj_flag == (int)PROJ_TYPE::EXT_SHELL) {
    collect_ext();
    c_base->set_ext_spur(true);
  }

  save_vtk_name("./output/mesh_created");

  f_calc.reload_cb(c_base);
  //sync_cins();
  //g_lens = new vd_glens(m, c_base);

  //vd_att_disp_field(m);
  //vd_att_vel_field(m);
  f_calc.vd_att_fields(m);
  f_calc.refresh(m, c_base, &e_list);

  load_flag = true;
  printf("VDSIM: Mesh loaded.\n");
  set_free_cells();

  sim_sz = std::cbrt(vd_tot_volm(m));
  first_time = false;

  get_length_scale();
}

void vd_sim::set_mesh(const char* modelFile, const char* meshFile) {

  // Load the mesh. set cellbase
  if (load_flag) {
    printf("VDSIM: Mesh already loaded, destroying first.\n");
    clean_mesh();
  }
  m = apf::loadMdsFromGmsh(gmi_load(modelFile), meshFile);
  mesh_flag = true;

  if(first_time) {
    c_base = new cell_base(m->getModel());
    reload_len();
  }
  e_list.change_mesh(m, c_base);
  if(first_time and proj_flag == (int)PROJ_TYPE::EXT_SHELL) {
    collect_ext();
    c_base->set_ext_spur(true);
  }

  f_calc.reload_cb(c_base);

  //sync_cins();
  //c_ins.load_cb(c_base);
  //g_lens = new vd_glens(m, c_base);

  //vd_att_disp_field(m);
  //vd_att_vel_field(m);
  f_calc.vd_att_fields(m);
  f_calc.refresh(m, c_base, &e_list);

  load_flag = true;

  printf("VDSIM: Mesh loaded.\n");
  set_free_cells();

  sim_sz = std::cbrt(vd_tot_volm(m));
  first_time = false;

  get_length_scale();
}

// TODO not compliant with the above function.
void vd_sim::set_mesh_smb(const char* modelFile, const char* meshFile) {

  if (load_flag) {
    printf("VDSIM: Mesh already loaded, destroying first.\n");
    clean_mesh();
  }

  m = apf::loadMdsMesh(gmi_load(modelFile), meshFile);
  mesh_flag = true;

  if(first_time) {
    c_base = new cell_base(m->getModel());
    reload_len();
  }
  e_list.change_mesh(m, c_base);

  if(first_time and proj_flag == (int)PROJ_TYPE::EXT_SHELL) {
    collect_ext();
    c_base->set_ext_spur(true);
  }

  save_vtk_name("./output/mesh_loaded");
  f_calc.reload_cb(c_base);

  //sync_cins();
  //c_ins.load_cb(c_base);
  //g_lens = new vd_glens(m, c_base);

  //vd_att_disp_field(m);
  //vd_att_vel_field(m);
  f_calc.vd_att_fields(m);
  f_calc.refresh(m, c_base, &e_list);

  load_flag = true;
  printf("VDSIM: Mesh loaded.\n");
  set_free_cells();

  sim_sz = std::cbrt(vd_tot_volm(m));
  first_time = false;

  get_length_scale();

}


cell_base* vd_sim::get_c_base() {

  if (load_flag)
    return c_base;

  return NULL;
}
/*
cell_ins_chk* vd_sim::get_c_ins() {

  //std::cout << "Load flag " << load_flag << std::endl;
  if (load_flag)
    return &c_ins;

  return NULL;
}
*/

// Check if the connectivity of the 0cell allows insertion (assume constant GBP).
bool vd_sim::ret_ins_gmi(int c0_id) {
  int c0_id_cb = c0_id - 1;
  bool insert = false;
  ent_conn* cells1 = new ent_conn();
  ent_conn* cells2 = new ent_conn();
  ent_conn* cells3 = new ent_conn();
  ent_conn* e_2sh = new ent_conn();
  c_base->get_conn_dim(3, 0, c0_id_cb, cells3);
  c_base->get_conn_dim(2, 0, c0_id_cb, cells2);
  vd_3c_det vd_3c;
  assert(e_list.e.at(0).at(c0_id_cb).at(0).size() == 1);
  vd_3c.reload(m, c_base, e_list.e.at(0).at(c0_id_cb).at(0).at(0));
  int c3_nbr = vd_3c.get_3c_nbr();
  if(c_base->chk_0cell_ext_gmi(c0_id)) {
    // If the number of 2shells and 3strata is larger than 4, candidate: 
    if(proj_flag == (int)PROJ_TYPE::EXT_SHELL) {
      assert(e_sh.chk_shell(0, c0_id_cb));
      shell sh_0c = e_sh.get_shell(0, c0_id_cb);
      if(sh_0c.dim == 2)
        c3_nbr = 1 + cells3->conn.size();
      else {
        e_sh.sh_base.get_conn_dim_gmi(2, sh_0c.dim, sh_0c.id, e_2sh);
        c3_nbr = e_2sh->conn.size() + c3_nbr;
        std::vector<int> sh2_bound(0);
        std::map<int, bool> sh2_map{};
        for(int c3_curr = 0; c3_curr < vd_3c.get_3c_nbr(); c3_curr++) {
          std::vector<int> c2_bound(0);
          c2_bound = vd_3c.get_bound_ids(3, c3_curr);
          sh2_bound.reserve(c2_bound.size());
          for(int j = 0; j < c2_bound.size(); j++) {
            int c2_curr = m->getModelTag(vd_3c.get_mdl(2, c2_bound.at(j)));
            if(c_base->get_cell_ext_gmi(2, c2_curr) ) {
              assert(e_sh.chk_shell(2, c2_curr - 1));
              shell sh_2c = e_sh.get_shell(2, c2_curr - 1);
              if(!sh2_map[sh_2c.id]) {
                sh2_bound.push_back(sh_2c.id);
                sh2_map[sh_2c.id] = true;
              }
            }
          }

        }
        if(sh2_bound.size() < e_2sh->conn.size())
          c3_nbr = c3_nbr + 1;
      }
    }
  }
  if(4 < c3_nbr) {
  //if((int)(cells3->conn.size()) > c3_nbr) {
    std::cout << "More than " << c3_nbr << " 3-cell" << std::endl;
    for(int j = 0; j < cells3->conn.size(); j ++)
      std::cout << "3c" << cells3->conn.at(j) + 1 << " ";
    std::cout << std::endl;
    insert = true;
  }
  else {
    bool c1_pass = true;
    for(int j = 0; j < cells1->conn.size(); j++) {
      if(c_base->get_conn_sz(1, cells1->conn.at(j)) == 1) {
        std::cout << "1c" << cells1->conn.at(j) + 1 
                  << "bounded only by"
                  << " 0c" << c0_id
                  << std::endl;
        c1_pass = false;
        //tag_0cell_ins.push_back(c0_id);
        j = cells1->conn.size();
        insert = true;
      }
    }
    if(c1_pass) {
      std::vector<ent_conn> e_con(0);
      e_con.resize(cells2->conn.size());
      for(int j = 0; j < cells2->conn.size(); j++) {
        c_base->get_conn_dim(3, 2, cells2->conn.at(j), &e_con.at(j));
        std::sort(e_con.at(j).conn.begin(), e_con.at(j).conn.end());
        std::cout << "2c" << cells2->conn.at(j) + 1 << ": ";
        for(int k = 0; k < e_con.at(j).conn.size(); k++) {
          std::cout << "3c" << e_con.at(j).conn.at(k) + 1 << " ";
        }
        std::cout << std::endl;
      }
      for(int j = 0; j < cells2->conn.size(); j++) {
        for(int k = j+1; k < cells2->conn.size(); k++) {
          if(e_con.at(j).conn.size() == 1 or e_con.at(k).conn.size() == 1) {
/*
            // Boundary 2cells have the same adjacencies.
            if(e_con.at(j).conn.size() == 1 and 
               e_con.at(k).conn.size() == 1) {
              if(e_con.at(j).conn.at(0) == e_con.at(k).conn.at(0)) {
                std::cout << "Upward 3-cell ext" << std::endl;
                tag_0cell_ins.push_back(c0_id);
                j = cells2->conn.size();
                k = j;
              }
            }
*/
          }
          // Spurious 1-cell
          else if(e_con.at(j).conn == e_con.at(k).conn) {
            std::cout << "Upward 3-cell the same" << std::endl;
            //tag_0cell_ins.push_back(c0_id);
            j = cells2->conn.size();
            k = j;
            insert = true;
          }
        }
      }
    }
  }
  if(!insert)
    tag_0cell_stable[c0_id] = true;

  delete cells1;
  delete cells2;
  delete cells3;
  delete e_2sh;
  return insert;
}

// This should be based on disjoint sets of entities. Count each disjoint 3-
// strata. For each disjoint 3strata check the following.
void vd_sim::ret_ins_gmi() {
  e_list.refresh();

  tag_0cell_ins.clear();
  tag_0cell_ins.reserve(c_base->get_sz(0));

  std::cout << "Finding insertible 0cells" << std::endl;
  for(int i = 0; i < c_base->get_sz(0); i++) {
    if(!c_base->is_free(0, i)) {

      if(!tag_0cell_stable[i+1] and ret_ins_gmi(i + 1))
        tag_0cell_ins.push_back(i+1);
    }
  }

}

void vd_sim::set_end(float time_e) {
  time_end = time_e;
}

void vd_sim::corr_0cell() {

  for (int i = 0; i < c_base->get_sz(0); i++) {
    ins_rad_map[i+1] = -1;
  }
  //tag_0cell_ins = c_ins.ret_ins_gmi();
  ret_ins_gmi();

  if(tag_0cell_ins.size() > 0)
    correct = true;

  if(correct) {
    for (int i = 0; i < tag_0cell_ins.size(); i++) {
      //time_curr = time_end + 0.001*i;
      //time_curr = time_curr + 0.0001*i;
      std::cout << "Trying inserting around 0cell" << tag_0cell_ins.at(i)
                << std::endl;
      //c_ins.print_path(tag_0cell_ins.at(i)-1);
      // If insertion is successful, check the same 0cell for more insertions.
      // Append the new 0cells to the tag_0cell_ins list.
      while(ins_cell(tag_0cell_ins.at(i))) {
        m->acceptChanges();
        m->verify();
      }
      m->acceptChanges();
      m->verify();
      if(save_mov)
        save_vtk_mov();
    }

    tag_0cell_ins.clear();
  }

}

void vd_sim::corr_1cell() {

  for(int i = 0; i < c_base->get_sz(1); i++) {
    if(!c_base->chk_1cell_ext_c(i) and 
       (!c_base->is_free(1, i) and c_base->fix_spur(1, i, true)) ) {
    //if((!c_base->is_free(1, i) and c_base->fix_spur(1, i, true)) ) {
      c_base->process_fix_list();
      sort_celem(&c_base->fix_list);

      std::cout << "To be replaced " << std::endl;
      for(int j = 0; j < c_base->fix_list.size(); j++) {
        std::cout << c_base->fix_list.at(j).first.first << "c"
                  << c_base->fix_list.at(j).first.second << " "
                  << c_base->fix_list.at(j).second.first << "c"
                  << c_base->fix_list.at(j).second.second << std::endl;
      }

      vd_cell_repl v_repl(m);
      //v_repl.set_cell(&fix_new);
      v_repl.set_cell(&c_base->fix_list);
      v_repl.repl_cell();
      std::cout << "Replaced " << std::endl;
      //c_base->print_ent();

      reload_mesh();
    }
  }

}

// TODO this doesn't need to be called at each time iteration by storing 0-cells
// that can be used for insertions and only update those affected by collapses
// and insertions.
// Although in the general case, this should be more carefully handled.
void vd_sim::find_0cell() {

  //tag_0cell_ins = c_ins.ret_ins_gmi();
  ret_ins_gmi();
  if(tag_0cell_ins.size() > 0)
    correct = true;

}

void vd_sim::start_sim() {
  chk_ma_swap(m);
  f_calc.vd_att_fields(m);

  if(ad_flag) {
    adapt();
  }

  tag_0cell_ins.reserve(c_base->get_sz(0));

  //corr_1cell();
  get_length_scale();
  vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_INS);
  //double t0 = PCU_Time();
  if(ins_flag)
    corr_0cell();
  double elapsed = vd_tim.get_time((int)T_COST_TYPE::T_COST_FLAG_INS);
  double elapsed_temp;
  std::cout << "Insertions took " << elapsed << std::endl;

  f_calc.vd_att_fields(m);

  if(f_v2move) {
    std::vector<apf::MeshEntity*> v_2move(0);
    get_v2move(v_2move);
    f_calc.vd_calc_vel(m, &v_2move);
  }
  else
    f_calc.vd_calc_vel(m);

  if(extract_flag) {
    extract_data();
    //extract_MS();
  }

  while(time_end > time_curr) {
    if(!cost_t_accum)
      vd_tim.reset_elapsed();

    // 0cells to check for cell insertion.
    tag_0cell_ins.reserve(c_base->get_sz(0));
    dt_adapt = 0;
    //if (not correct) {
      assert(!check_ma_sgn());

      if(ad_flag) {
        adapt();
        // TODO This was a sanity check to identify that it was indeed 
        // meshAdapt that was introducing bad elements after starting with 
        // mesh generated by SCOREC using a tess file. It is potentially 
        // spurious now.
        assert(!check_ma_sgn());
      }

      vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_ITER);
      if(f_v2move)
        iter_cells();
      else
        iter_step();
      elapsed = vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_ITER);
      std::cout << "Iter_step took " << elapsed << std::endl;

      len_col_f = false;

      //vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
      //e_list.refresh();
      //elapsed_temp = vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_AUX);
      //elapsed = vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
      //std::cout << "e_list.refresh took " << elapsed-elapsed_temp << std::endl;

      vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
      if(len_f and !len_col_f) {
        clear_rad_map();
        avg_cell = upd_cell_rad_sim();
      }

      elapsed_temp = vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_AUX);
      elapsed = vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
      std::cout << "Update rad took " << elapsed-elapsed_temp << std::endl;

      vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_ITER);
      if(proj_flag == (int)PROJ_TYPE::EXT_SHELL) {
        f_calc.corr_pos(m);
      }
      elapsed_temp = vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_ITER);
      elapsed = vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_ITER);
      std::cout << "corr_pos took " << elapsed-elapsed_temp << std::endl;

      if(ad_flag)
        adapt();
      vd_rem_tag(m);
      vd_tag_mesh(m);
      if(!col_flag and limit_flag and t_sub < t_sub_limit) {
        std::cout << "breaking: t_sub is " << t_sub << std::endl;
        break;
      }
      if(flag_g_nbr) {
        int g_nbr = c_base->get_sz(3) - c_base->get_free_sz(3);
        if(g_nbr <= cond_g_nbr) {
          std::cout << "Remaining grains " << g_nbr << std::endl;
          break;
        }
      }
      if(flag_s_nbr) {
        int s_nbr = c_base->get_sz(3) - c_base->get_free_sz(3);
        if(s_nbr <= cond_s_nbr) {
          std::cout << "Remaining surfaces " << s_nbr << std::endl;
          break;
        }
      }
      if(flag_l_nbr) {
        int l_nbr = c_base->get_sz(3) - c_base->get_free_sz(3);
        if(l_nbr <= cond_l_nbr) {
          std::cout << "Remaining lines " << l_nbr << std::endl;
          break;
        }
      }

    //}
    //iter_step();

    cell_2_col c2c(-1,-1);
    if(col_flag)
      c2c.first = 1;
    else
      c2c.first = -1;

    vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
    get_length_scale();

    correct = false;
    if(ins_flag)
      find_0cell();

    elapsed_temp = vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_AUX);
    elapsed = vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
    std::cout << "len_topo and find_0cell took " << elapsed-elapsed_temp << std::endl;

    vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
    int save_count = 1;
    if(!ad_flag) {
      clear_rad_map();
      avg_cell = upd_cell_rad_sim();
    }
    elapsed_temp = vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_AUX);
    elapsed = vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_AUX);
    std::cout << "cell rad took " << elapsed-elapsed_temp << std::endl;

    vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_COL);
    while(c2c.first > 0) {
      e_list.refresh();
      c2c.first = 0;
      // To reduce the computations, map checked cells.
      c2c = chk_small_cell_rad();
      if (c2c.first > 0) {
        //e_list.refresh();
        //vd_glens g_lens2(m, c_base, &e_list);
        //g_lens2.set_field_calc(&f_calc);
        //reload_mesh();
        //if(g_lens2.test_cvx(c2c.first, c2c.second)) {
          std::pair<int, int> tag_cell = col_cell(c2c.first, c2c.second);
          f_calc.vd_att_fields(m);

          if(tag_cell.first != -2) {
            if(ad_flag) {
              adapt();
            }
            e_list.refresh();

            if(sub_vtk)
              save_vtk_name("output/after_adapt_col");
            chk_ma_swap(m);
          }
          // Assume collapse fixes spurious entities already.

          if(tag_cell.first == 0)
            tag_0cell_ins.push_back(tag_cell.second);

          //e_list.refresh();

          correct = true;
          //time_curr = time_end + 0.001*save_count;
          //time_curr = time_curr + 0.0001*save_count;
          save_count++;
          if(save_mov)
            save_vtk_mov();
        //}
        //else {
        //  std::cout << "A convex hull cannot be created." << std::endl;
        //  rad_map_cc.at(c2c.first - 1)[c2c.second] = true;
        //}
      }
    }
    elapsed_temp = vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_COL);
    elapsed = vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_COL);
    std::cout << "collapses took " << elapsed-elapsed_temp << std::endl;

    vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_INS);
    if(correct) {
      if(ad_flag)
        adapt();
      for (int i = 0; i < tag_0cell_ins.size(); i++) {
        //time_curr = time_end + 0.001*save_count;
        //time_curr = time_curr + 0.0001*save_count;
        save_count++;
        std::cout << "Trying inserting around 0cell" << tag_0cell_ins.at(i)
                  << std::endl;
        //c_ins.print_path(tag_0cell_ins.at(i)-1);
        // If insertion is successful, check the same 0cell for more 
        // insertions. Append the new 0cells to the tag_0cell_ins list.
        //adapt();
        while(ins_cell(tag_0cell_ins.at(i))) {
          m->acceptChanges();
          m->verify();
          if(ad_flag)
            adapt();
        }
        m->acceptChanges();
        m->verify();
        if(save_mov)
          save_vtk_mov();
      }

      tag_0cell_ins.clear();
    }
    elapsed_temp = vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_INS);
    elapsed = vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_INS);
    std::cout << "insertions took " << elapsed-elapsed_temp << std::endl;

    if(chk_save_vtk_interval()) {
      if(save_mov)
        save_vtk_mov();
      else
        save_vtk();
    }
    std::cout << "Time taken: \ndt_it: " << vd_tim.get_elapsed_sum()
              << " dt_iter: " << vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_ITER)
              << " dt_aux: " << vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_AUX)
              << " dt_col: " << vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_COL)
              << " dt_ins: " << vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_INS)
              << " dt_adapt(intrinsic to the other times): " << vd_tim.get_elapsed((int)T_COST_TYPE::T_COST_FLAG_ADAPT)
              << std::endl;
    if(!cost_t_accum) {
      write_time_cost();
    }
  }
  if(save_mov)
    save_vtk_mov();
  else
    save_vtk();
  if(cost_t_accum) {
    write_time_cost();
  }
}

// Mesh save and load functions:
void vd_sim::save_mesh(char const* FileName) {

  std::stringstream ss;
  ss << FileName << ".dmg";
  std::string tmp = ss.str();
  const char* cstr = tmp.c_str();

  c_base->vd_write_dmg(cstr);

  ss.str("");

  ss.clear();

  ss << FileName << ".smb";
  tmp = ss.str();
  cstr = tmp.c_str();

  m->writeNative(cstr);
  e_list.refresh();
  //m = apf::loadMdsFromGmsh(gmi_load(modelFiletemp), meshFile);

}

void vd_sim::reload_len() {
  if(sim_len < std::numeric_limits<double>::min()) {
    apf::Vector3 min(0,0,0);
    apf::Vector3 max(0,0,0);
    apf::Vector3 p(0,0,0);
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* v;
    assert(v = m->iterate(it));
    m->getPoint(v, 0, p);
    min = p;
    max = p;

    while (v = m->iterate(it)) {
      m->getPoint(v, 0, p);
      if(min[0] > p[0])
        min[0] = p[0];
      if(min[1] > p[1])
        min[1] = p[1];
      if(min[2] > p[2])
        min[2] = p[2];

      if(max[0] < p[0])
        max[0] = p[0];
      if(max[1] < p[1])
        max[1] = p[1];
      if(max[2] < p[2])
        max[2] = p[2];
    }
    m->end(it);
    sim_len = (max-min).getLength();
  }
}

void vd_sim::check_model() {
  for(int dim = 0; dim < 4; dim++) {
    for(int i = 0; i < c_base->get_sz(dim); i++) {
      assert(m->findModelEntity(dim, i + 1));
    }
  }

}

void vd_sim::reload_mesh() {

  if (load_flag) {

    c_base->vd_write_dmg("./temp_save.dmg");
    m->writeNative("./temp_save.smb");

    //clean_mesh();
    std::cout << "Mesh reloaded." << std::endl;
    set_mesh_smb("./temp_save.dmg", "./temp_save.smb");

    f_calc.vd_att_fields(m);

    check_model();
  }

}

// Print the model object. 
void vd_sim::print_mesh_topo() {

  struct gmi_model* mdl = m->getModel();

  std::cout << "gmi entities:" << std::endl;
  for(int i = 0; i < 4; i++) {
    struct gmi_iter* g_it = gmi_begin(mdl, i);
    struct gmi_ent* g_ent;

    while ((g_ent = gmi_next(mdl, g_it))) {
      std::cout << "dim:" << i << " " 
                << gmi_dim(mdl,g_ent) << "c" 
                << gmi_tag(mdl,g_ent) << std::endl;

      if(i < 3) {
        struct gmi_set* s = gmi_adjacent(mdl, g_ent, i+1);
        std::cout << "Up ";
        for (int i = 0 ; i < s->n; i++) {
          std::cout << gmi_tag(mdl, s->e[i]) << " ";
        }
        gmi_free_set(s);
        std::cout << std::endl;
      }
      if(i > 0) {
        struct gmi_set* s = gmi_adjacent(mdl, g_ent, i-1);
        std::cout << "Down ";
        for (int i = 0 ; i < s->n; i++) {
          std::cout << gmi_tag(mdl, s->e[i]) << " ";
        }
        gmi_free_set(s);
        std::cout << std::endl;
      }
    }
    gmi_end(mdl, g_it);
  }

  c_base->print_ent();
}

void vd_sim::set_save_vtk_interval(double dt_period_in) {
  dt_vtk_period = dt_period_in;
  t_save_last = time_curr;
}

bool vd_sim::chk_save_vtk_interval() {
  if(dt_vtk_period < std::numeric_limits<double>::min())
    return true;

  double temp = time_curr - t_save_last;
  if(temp > dt_vtk_period) {
    temp = std::fmod(temp, dt_vtk_period);
    t_save_last = time_curr - temp;
    return true;
  }
  return false;
}

void vd_sim::save_vtk() {

  vd_rem_tag(m);
  vd_tag_mesh(m);

  std::stringstream ss;
  ss.precision(2);
  ss << vtk_name << std::fixed << time_curr*10000;
  std::string tmp = ss.str();
  const char* cstr = tmp.c_str();

  apf::writeVtkFiles(cstr, m);

  ss.str("");
  ss.clear();
  ss << smb_name << std::fixed << time_curr*10000;

  std::string::size_type dot_pos = 1;
  tmp = ss.str();
  dot_pos = tmp.find ('.', dot_pos);
  assert(dot_pos != std::string::npos);

  tmp = tmp.substr(0, dot_pos) + "_" + 
                            tmp.substr(dot_pos+1, tmp.length()) + ".smb";
  
  cstr = tmp.c_str();

  //m->writeNative(smb_name.c_str());
  m->writeNative(cstr);

  e_list.refresh();
  //ss.str("");
  //ss.clear();
  //ss << "./tempmesh/temp.dmg";
  //tmp = ss.str();
  tmp = tmp.substr(0, tmp.length()-3) + "dmg";
  cstr = tmp.c_str();

  c_base->vd_write_dmg(cstr);
  vd_rem_tag(m);

}

// Generates a vtk per grain. Used in generating movies in a hackish way.
// TODO this can be safely discarded now that Paraview python interface allows
// selection of individual grains.
void vd_sim::save_vtk_mov() {

  vd_rem_tag(m);
  vd_tag_mesh(m);

  // Adding this part to generate a video or gif of the grain growth 
  // simulation in a hackish way.
  std::vector<std::vector<apf::MeshEntity*> > es_set_save(4, 
                                  std::vector<apf::MeshEntity*>(0));
  for (int i = 0; i < c_base->get_sz(3); i++) {
    es_set_save.at(3) = e_list.e.at(3).at(i).at(3);
    for (int j = 3; j > 0; j--) {
      vd_set_down(m, &es_set_save.at(j), &es_set_save.at(j-1));
    }

    std::stringstream ss;
    ss.precision(2);
    ss << vtk_name << "cell" << i+1 << "t"
       << std::fixed << time_curr*10000;
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    vd_save_vtk_set(m, &es_set_save, cstr);
  }
  vd_rem_tag(m);

}

void vd_sim::save_vtk_name(const char* vtk_out, bool smb_out) {

  vd_rem_tag(m);
  vd_tag_mesh(m);
  apf::writeVtkFiles(vtk_out, m);

  if(smb_out) {
    std::stringstream ss;

    ss << smb_name << ".smb";
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    m->writeNative(cstr);

    c_base->vd_write_dmg("./tempmesh/temp_iter.dmg");
    m->writeNative("./tempmesh/temp_iter.smb");
    e_list.refresh();
  }
}

apf::Mesh2* vd_sim::get_mesh() {
  return m;
}

void vd_sim::verify() {
  m->verify();
}

std::pair<double,double> vd_sim::get_length_scale() {
  //return ma::getAverageEdgeLength(m)/ratio_ent;
  // This number should be related to refinement (and either directly or 
  // indirectly through the refinement, to the number of cells.)

  //if(len_f and !len_col_f) {
  //  return avg_cell.at(0)/adapt_len/ad_param/ratio_ent;
  //}
  //else
  //  return getMedianEntSize_hist(m, 1, 100)/ratio_ent;
  std::vector<double> c1_len(0);
  ln_min = sm_map.at(0)[1];
  c1_len.reserve(c_base->get_sz(1));
  for(int i = 0; i < c_base->get_sz(1); i++) {
    if(!c_base->is_free(1, i)) {
      c1_len.push_back(sm_map.at(0)[i+1]);
      if(ln_min > sm_map.at(0)[i + 1])
        ln_min = sm_map.at(0)[i + 1];
    }
  }

  if(len_f and !len_col_f) {
    avg_cell = upd_cell_rad_sim();
  }
  if(ad_type == ADAPT_TYPE::ADAPT_1CELL) {
    //return adapt_len*avg_cell.at(0)/ad_param/ma::getAverageEdgeLength(m);
    //return adapt_len*ad_param*avg_cell.at(0);
    //return avg_cell.at(0)/adapt_len/ad_param;
    if(c1_len.size() > 10)
      len_topo = getMedianEntSize_hist(c1_len, 100);
    else {
      len_topo = 0;
      for(int i = 0; i < c1_len.size(); i++)
        len_topo = len_topo + c1_len.at(i);
      len_topo = len_topo/c1_len.size();
    }
  }
  else {
    assert(ad_type == ADAPT_TYPE::ADAPT_3CELL or
           ad_type == ADAPT_TYPE::ADAPT_STEP or
           ad_type == ADAPT_TYPE::ADAPT_STEP_1CELL or
           ad_type == ADAPT_TYPE::ADAPT_BOUND or
           ad_type == ADAPT_TYPE::ADAPT_CURVE);

    //return adapt_len*std::cbrt(avg_cell.at(2))
    //                  /ad_param/ma::getAverageEdgeLength(m);
    //return adapt_len*ad_param*std::cbrt(avg_cell.at(2));
    //return avg_cell.at(2)/adapt_len/ad_param;
    std::vector<double> c3_vol(0);
    c3_vol.reserve(c_base->get_sz(3));
    for(int i = 0; i < c_base->get_sz(3); i++) {
      if(!c_base->is_free(3, i)) {
        c3_vol.push_back(sm_map.at(2)[i+1]);
      }
    }
    if(c3_vol.size() > 4)
      len_topo = getMedianEntSize_hist(c3_vol, 100);
    else {
      len_topo = 0;
      for(int i = 0; i < c3_vol.size(); i++)
        len_topo = len_topo + c3_vol.at(i);
      len_topo = len_topo/c3_vol.size();
    }
  }

  len_trans = len_topo/ratio_ent;

  q_th = calc_good_vol(m, ln_min, len_trans);

  std::cout << "len_topo is " << len_topo << " len_trans is " << len_trans
            << " q_th is " << q_th
            << std::endl;
  return std::make_pair(len_topo, len_trans);
}

// Update the normal distance field for the vertices of the given cell, from
// the one lower dimensional entities on the bounding strata.
void vd_sim::upd_dist_field(int c_dim, int c_id, 
                            std::map<apf::MeshEntity*, double> &v_dist) {
  if(c_dim == 0)
    return;

  int v_sz = c_dim + 1;

}

// Assuming an equiaxed ball of appropriate dimension of radius len_topo*2, 
// check whether the length/area/volume associated with the cell is smaller
// than that associated with the ball.
// TODO change the name of chk_cell_thresh or this one, as it is confusing.
// Also, chk_cell_thresh seems outdated, so remove if not to be recycled.
bool vd_sim::chk_cell_th(int d, int cell_id, double &r_equi) {
  double w_th = 0;
  if(d == 1)
    w_th = 2*len_trans;
  else if(d == 2)
    w_th = len_trans*len_trans*PI_L;
  else {
    assert(d == 3);
    w_th = std::pow(len_trans, 3.)*PI_L*4./3.;
  }
  double w_tot = vd_meas_set(m, &e_list.e.at(d).at(cell_id-1).at(d));
  //double w_tot = 0;
  //for(int i = 0; i < e_list.e.at(d).at(cell_id-1).at(d).size(); i++) {
  //  w_tot = w_tot + vd_meas_ent(m, e_list.e.at(d).at(cell_id-1).at(d).at(i));
  //}
  if(d == 1)
    r_equi = w_tot/2;
  else if(d == 2) {
    //r_equi = std::pow(w_tot/PI_L, 0.5);
    r_equi = std::pow(w_tot/8., 0.5);
  }
  else {
    assert(d == 3);
    //r_equi = std::pow(w_tot/PI_L/4.*3., 1.0/3);
    r_equi = std::pow(w_tot/16., 1.0/3);
  }
  meas_map.at(d-1)[cell_id] = w_tot;
  sm_map.at(d-1)[cell_id] = r_equi;

  if(w_tot > len_trans/ratio_col_sh)
    t_below_map.at(d-1)[cell_id] = -1;

  if(sm_map.at(d-1)[cell_id] > len_trans/ratio_col_sh) {
    rad_map.at(d-1)[cell_id] = true;
  }
  else
    rad_map.at(d-1)[cell_id] = false;


  if(w_tot > w_th)
    return true;
  else
    return false;
}

std::vector<double> vd_sim::upd_cell_rad_sim() {
  //double average_ent = ma::getAverageEdgeLength(m);
  double average_ent = getAverageEntSize(m, 1);
  f_calc.set_drag_glob(getAverageEntSize(m,2));

  std::cout << "Average, edge: " << average_ent << std::endl;
  double cell_meas_max;
  double cell_meas_min;

  //std::vector<double> avg_cell(0);
  //avg_cell.resize(3);

  std::vector<apf::MeshEntity*> ev(0);
  std::vector<apf::MeshEntity*> ev_e(0);
  std::vector<apf::MeshEntity*> ev_d(0);

  gmi_model* mdl = m->getModel();

  struct gmi_iter* it;
  struct gmi_ent* e;

  for (int d = 3; d > 0; d--) {

    min_cell.at(d-1) = -1;

    it = gmi_begin(mdl, d);
    int counter = 0;
    while ((e = gmi_next(mdl, it))) {
      int cell_id = gmi_tag(mdl, e);

      if(e_list.e.at(d).at(cell_id-1).at(d).size() > 0) {
        //std::cout << d << "c" << cell_id << std::endl;
        double cell_meas_min;
        chk_cell_th(d, cell_id, cell_meas_min);
        if(cell_meas_min < min_cell.at(d-1) or min_cell.at(d-1) < 0)
          min_cell.at(d-1) = cell_meas_min;

        counter++;
        //std::cout << cell_meas_min << " counter " << counter << std::endl;
        avg_cell.at(d-1) = avg_cell.at(d-1) + cell_meas_min;
      }
    }
    gmi_end(mdl, it);

    avg_cell.at(d-1) = avg_cell.at(d-1)/counter;
    std::cout << "Average " << d << "c " << avg_cell.at(d-1) << std::endl;
  }
  len_col_f = true;
  if(c_base->get_sz(1) - c_base->get_free_sz(1) == 0)
    min_ins_rad = get_adapt_ln();
  else if(min_ins_rad < std::numeric_limits<double>::min())
    min_ins_rad = min_cell.at(0);
  return avg_cell;

}


// For stepwise graded adaptation, mark the vertices.
void vd_sim::adapt_mark_step() {

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");

  apf::MeshEntity* e;
  apf::ModelEntity* mdl;

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it))) {
    mdl = m->toModel(e);
    if(m->getModelType(mdl) < 3)
      apf::setScalar(field_step, e, 0, 1);
    else
      apf::setScalar(field_step, e, 0, 3);
  }
  m->end(it);

  apf::Up up;
  it = m->begin(0);
  while ((e = m->iterate(it))) {
    mdl = m->toModel(e);
    if(m->getModelType(mdl) < 3) {
      m->getUp(e, up);
      for(int i = 0; i < up.n; i++) {
        apf::MeshEntity* v2 = apf::getEdgeVertOppositeVert(m, up.e[i], e);
        if(m->getModelType(mdl) == 3) {
          apf::setScalar(field_step, v2, 0, 3);
        }
      }
    }
  }
  m->end(it);

}



// In addition to stepwise graded adaptation, mark the vertices connected to 1c
// vertices such that the minimum 1c has two edges and the surrounding edges 
// are of similar length.
// The ratio of the minimum 1c length to the reference length should by larger
// than 3. 
void vd_sim::adapt_mark_1c_min() {

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");

  apf::MeshEntity* e;
  apf::ModelEntity* mdl;

  double ref_len = get_adapt_ln();
  std::cout << "ref_len " << ref_len << std::endl;

  double e_weight = 1;

  apf::Up up;
  apf::Downward d_v;

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it))) {
    apf::setScalar(field_step, e, 0, 1);
  }
  m->end(it);

  it = m->begin(1);
  while ((e = m->iterate(it))) {
    m->getDownward(e, 0, d_v);
    mdl = m->toModel(e);
    int c_type = m->getModelType(mdl);
    int c_tag = m->getModelTag(mdl);
    e_weight = 1/sm_map.at(c_type - 1)[c_tag];
    if(c_type == 1) {
      //if(e_weight < 1)
      //  e_weight = e_weight*e_weight;
      e_weight = e_weight/16;
    }
    else if(c_type == 3) {
      //if(e_weight < 1)
      //  e_weight = e_weight*e_weight;
      e_weight = e_weight*3;
    }
    if(c_type < 3 and e_weight < 1) {

      double w1 = apf::getScalar(field_step, d_v[0], 0);
      if(w1 > e_weight) {
        apf::setScalar(field_step, d_v[0], 0, e_weight);
        //std::cout << d_v[0] << " "
        //          << m->getModelType(m->toModel(d_v[0])) << "c"
        //          << m->getModelTag(m->toModel(d_v[0]))
        //          << " " << e_weight << std::endl;
      }
      w1 = apf::getScalar(field_step, d_v[1], 0);
      if(w1 > e_weight) {
        apf::setScalar(field_step, d_v[1], 0, e_weight);
        //std::cout << d_v[1] << " "
        //          << m->getModelType(m->toModel(d_v[1])) << "c"
        //          << m->getModelTag(m->toModel(d_v[1]))
        //          << " " << e_weight << std::endl;
      }
    }
    else {
      double w1 = apf::getScalar(field_step, d_v[0], 0);
      if(w1 > 1) {
        apf::setScalar(field_step, d_v[0], 0, 1);
        //std::cout << d_v[0] << " "
        //          << m->getModelType(m->toModel(d_v[0])) << "c"
        //          << m->getModelTag(m->toModel(d_v[0]))
        //          << " " << 1 << std::endl;
      }
      //else {
        //std::cout << d_v[0] << " "
        //          << m->getModelType(m->toModel(d_v[0])) << "c"
        //          << m->getModelTag(m->toModel(d_v[0]))
        //          << " " << w1 << std::endl;
      //}
      w1 = apf::getScalar(field_step, d_v[1], 0);
      if(w1 > 1) {
        apf::setScalar(field_step, d_v[1], 0, 1);
        //std::cout << d_v[1] << " "
        //          << m->getModelType(m->toModel(d_v[1])) << "c"
        //          << m->getModelTag(m->toModel(d_v[1]))
        //          << " " << 1 << std::endl;
      }
      //else {
      //  std::cout << d_v[1] << " "
      //            << m->getModelType(m->toModel(d_v[1])) << "c"
      //            << m->getModelTag(m->toModel(d_v[1]))
      //            << " " << w1 << std::endl;
      //}
    }
  }
  m->end(it);
/*
  it = m->begin(1);
  while ((e = m->iterate(it))) {
    m->getDownward(e, 0, d_v);

    double w1 = apf::getScalar(field_step, d_v[0], 0);
    double w2 = apf::getScalar(field_step, d_v[1], 0);

    // e_weight = w1/2+w2/2;
    e_weight = vd_meas_ent(m, e)/ref_len;
    if(e_weight < 1) {
      e_weight = 0.5;
    }
    else {
      e_weight = exp(-e_weight);
    }

    if(w1 < w2) {
      apf::setScalar(field_step, d_v[1], 0, (1-e_weight)*w2 + e_weight*w1);
    }
    else if (w2 < w1) {
      apf::setScalar(field_step, d_v[0], 0, (1-e_weight)*w1 + e_weight*w2);
    }
    //std::cout << d_v[0] << " "
    //          << m->getModelType(m->toModel(d_v[0])) << "c"
    //          << m->getModelTag(m->toModel(d_v[0]))
    //          << " " << w1 << " " 
    //          << d_v[1] << " "
    //          << m->getModelType(m->toModel(d_v[1])) << "c"
    //          << m->getModelTag(m->toModel(d_v[1]))
    //          << " " << w2
    //          << std::endl;

  }
  m->end(it);

  std::cout << "Rescaling weights based on cell membership" << std::endl;
  it = m->begin(0);
  while ((e = m->iterate(it))) {
    mdl = m->toModel(e);
    int c_type = m->getModelType(mdl);
    int c_tag = m->getModelTag(mdl);
    double w1 = apf::getScalar(field_step, e, 0);
    if(c_type > 0)
      apf::setScalar(field_step, e, 0, w1*c_type);
    //std::cout << e << " "
    //          << c_type << "c"
    //          << c_tag
    //          << " " << w1 << " " << w1*c_type
    //          << std::endl;

  }
  m->end(it);
*/
}

// Mark the vertices of the collapsing cell with the weight 1.
void vd_sim::adapt_mark_0c_col(int d_col, int t_col) {
  apf::Field* field_step = m->findField("adapt_step");
  assert(field_step);
  std::vector<std::vector<apf::MeshEntity*> > ent_d(d_col, 
                                      std::vector<apf::MeshEntity*>(0));
  vd_set_down(m, &e_list.e.at(d_col).at(t_col-1).at(d_col), &ent_d.at(d_col-1));

  for(int d = d_col-1; d > 0; d--) {
    vd_set_down(m, &ent_d.at(d), &ent_d.at(d-1));
  }
  double temp = 0;
  std::cout << "Setting adapt step. " << std::endl;
  for(int i = 0; i < ent_d.at(0).size(); i++) {
    apf::setScalar(field_step, ent_d.at(0).at(i), 0, 10);
    //temp = apf::getScalar(field_step, ent_d.at(0).at(i), 0);
    //std::cout << ent_d.at(0).at(i) << " " << temp << std::endl;
  }

}

// Mark the vertices based on the distances to the 0cell vertices, bounding the
// cell they belong to.
void vd_sim::adapt_mark_0c_min() {

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");

  apf::ModelEntity* mdl;

  double ref_len = get_adapt_ln();
  std::cout << "ref_len " << ref_len << std::endl;

  double e_weight = 1/ref_len;
  apf::MeshEntity* e;

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it))) {
    apf::setScalar(field_step, e, 0, e_weight);
  }
  m->end(it);
/*
  std::map<int, apf::Vector3> pos_0c{};
  apf::Vector3 temp(0,0,0);
  for(int i = 0; i < e_list.e.at(0).size(); i++) {
    if(e_list.e.at(0).at(i).at(0).size() == 1) {
      m->getPoint(e_list.e.at(0).at(i).at(0).at(0), 0, temp);
      pos_0c[i] = temp;
    }
  }
  // Going over each cell of d > 0, check each vertex for distances to bounding
  // 0cells.

  ent_conn* e_0c = new ent_conn();
  for(int d = 1; d < e_list.e.size(); d++) {
    for(int i = 0; i < e_list.e.at(d).size(); i++) {
      c_base->get_conn_dim(0, d, i, e_0c);
      for(int j = 0; j < e_list.e.at(d).at(i).at(0).size(); j++) {
        double dist = -1;
        double dist_curr = -1;
        m->getPoint(e_list.e.at(d).at(i).at(0).at(j), 0, temp);
        for(int k = 0; k < e_0c->conn.size(); k++) {
          dist_curr = (temp - pos_0c[e_0c->conn.at(k)]).getLength();
          if(dist < 0 or dist_curr < dist)
            dist = dist_curr;
        }
        if(dist < 0) {
          apf::setScalar(field_step, e_list.e.at(d).at(i).at(0).at(j), 0, 
                                                   sm_map.at(d)[i+1]/ref_len);
          //apf::setScalar(field_step, e_list.e.at(d).at(i).at(0).at(j), 0, 
          //                                         sm_map.at(d)[i+1]);
        }
        else {
          apf::setScalar(field_step, e_list.e.at(d).at(i).at(0).at(j), 0, 
                                                                dist/ref_len);
          //apf::setScalar(field_step, e_list.e.at(d).at(i).at(0).at(j), 0, 
          //                                                      dist);
        }
      }
    }
  }

  apf::Up up;
  apf::MeshEntity* v_curr;
  apf::MeshEntity* v_other;
  for(int i = 0; i < e_list.e.at(0).size(); i++) {
    // 0c being used
    if(e_list.e.at(0).at(i).at(0).size() == 1) {
      v_curr = e_list.e.at(0).at(i).at(0).at(0);
      c_base->get_conn_dim(0, 0, i, e_0c);
      double dist = -1;
      double dist_curr = -1;

      // 0c is connected to other 0cells over 1cells
      if(e_0c->conn.size() > 0) {
        for(int k = 0; k < e_0c->conn.size(); k++) {
          dist_curr = (pos_0c[i] - pos_0c[e_0c->conn.at(k)]).getLength();
          //std::cout << "0c" << e_0c->conn.at(k)+1 << " "
          //          << dist_curr << std::endl;

          if(dist < 0 or dist_curr < dist)
            dist = dist_curr;
        }
      }
      else {
        dist = 1;
      }
      assert(dist > 0);
      //std::cout << "0c" << i+1 << " " << v_curr << " " << dist << std::endl;
      apf::setScalar(field_step, v_curr, 0, dist/ref_len);
      //apf::setScalar(field_step, v_curr, 0, dist);

    }
  }
  delete e_0c;
*/
}

// Curvature based mesh refinement, mark the vertices.
// Top-down calculation, should be efficient.
// H_i = 1/2*sum_{j in N_i} cot(alpha_ij) + cot(beta_ij)(len e_j)
// where Ni is set of edges adjacent to vertex i, alpha_ij and beta_ij are the 
// angles of the surface triangles adjcent to edge j, across the edge j.
// Loop over surface triangles, calculate the interior angles, add the cot(angle) 
// to the field_cot on the edge. 
// Loop over the edges, 
// if edge belongs to 2-stratum: add 1/2 of field_cot*len_e to the vertex
// if edge belongs to 1-stratum: add 1/2n of field_cot*len_e to the vertex, where
// n is the number of adjacent 3-strata of the vertex.
void vd_sim::adapt_mark_curv() {
  int lookup_tri_ed [3][2] = {{2,0},{0,1},{1,2}};

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  apf::Field* field_curv = vd_att_vs_field(m, "curv_mean");
  apf::Field* field_cot = vd_att_es_field(m, "cot_field");

  apf::ModelEntity* mdl;
  apf::Downward d_v;
  apf::Downward d_e;

  std::vector<double> ang(3, 0.);
  //std::vector<apf::vector3> (3, apf::Vector3(0,0,0));
  apf::Vector3 temp_pos(0,0,0);

  double ref_len = get_adapt_ln();
  std::cout << "ref_len " << ref_len << " bound_len " << bound_len << std::endl;

  double e_weight = 1./ref_len;
  double b_weight = 1./bound_len;
  apf::MeshEntity* e;

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it))) {
    apf::setScalar(field_step, e, 0, e_weight);
    apf::setScalar(field_curv, e, 0, 0.);
  }
  m->end(it);

  std::map<apf::MeshEntity*, apf::Vector3> ori_map{};
  std::map<apf::MeshEntity*, apf::Vector3> pos_map{};

  for(int dim = 1; dim < 3; dim++) {
    int sz = c_base->get_sz(dim);
    for(int c_id = 0; c_id < sz; c_id++) {
      if(!c_base->is_free(dim, c_id)) {
        std::vector<apf::MeshEntity*>* ents = &e_list.e.at(dim).at(c_id).at(1);
        for(int i = 0; i < ents->size(); i++) {
          ori_map[ents->at(i)] = get_edge_dir(m, ents->at(i));
          pos_map[ents->at(i)] = vd_get_pos(m, ents->at(i));
        }
      }
    }
  }

  it = m->begin(2);
  while ((e = m->iterate(it))) {
    if(m->getModelType(m->toModel(e)) == 2) {
      m->getDownward(e, 0, d_v);
      m->getDownward(e, 1, d_e);

      for(int v1 = 0; v1 < 3; v1++) {
        int e1 = lookup_tri_ed[v1][0];
        int e2 = lookup_tri_ed[v1][1];

        m->getPoint(d_v[v1], 0, temp_pos);

        double coef = 1.;
        temp_pos = norm_0(pos_map[d_e[e1]] - temp_pos);
        if(temp_pos*ori_map[d_e[e1]] < -std::numeric_limits<double>::min())
          coef = -coef;
        temp_pos = norm_0(pos_map[d_e[e2]] - temp_pos);
        if(temp_pos*ori_map[d_e[e2]] < -std::numeric_limits<double>::min())
          coef = -coef;
        double ang_cos = ori_map[d_e[e2]]*ori_map[d_e[e1]]*coef;

        double ang_curr = std::acos(std::min(std::max(ang_cos,-1.0),1.0));
        double ang_sin = std::sin(ang_curr);

        if(std::fabs(ang_sin) > std::numeric_limits<double>::min()) { 
          e1 = 3 - e1 + e2;
          double cot_curr = apf::getScalar(field_cot, d_e[e1], 0);
          apf::setScalar(field_cot, d_e[e1], 0, cot_curr + ang_cos/ang_sin);
        }
      }
    }
  }
  m->end(it);

  for(int dim = 1; dim < 3; dim++) {
    int sz = c_base->get_sz(dim);
    for(int c_id = 0; c_id < sz; c_id++) {
      if(!c_base->is_free(1, c_id)) {
        std::vector<apf::MeshEntity*>* ents = &e_list.e.at(dim).at(c_id).at(1);

        for(int i = 0; i < ents->size(); i++) {
          apf::MeshEntity* e_curr = ents->at(i);
          m->getDownward(e_curr, 0, d_v);
          double len = vd_meas_ent(m, e_curr);
          double cot_curr = apf::getScalar(field_cot, e_curr, 0);
          double curv_mean = apf::getScalar(field_curv, d_v[0], 0);

          apf::setScalar(field_curv, d_v[0], 0, curv_mean + cot_curr*len);
        }
      }
    }
  }

  for(int dim = 0; dim < 3; dim++) {
    int sz = c_base->get_sz(dim);
    for(int c_id = 0; c_id < sz; c_id++) {
      if(!c_base->is_free(dim, c_id)) {

        ent_conn* eup = new ent_conn();
        c_base->get_conn_dim(3, 1, c_id, eup);
        int n = eup->conn.size();
        delete eup;

        std::vector<apf::MeshEntity*>* ents = &e_list.e.at(dim).at(c_id).at(0);
        for(int i = 0; i < ents->size(); i++) {
          double curv_mean = apf::getScalar(field_curv, d_v[0], 0);
          if(n/curv_mean > ref_len)
            apf::setScalar(field_step, d_v[0], 0, n/curv_mean);
        }
      }
    }
  }
  apf::destroyField(field_cot);

  if(sub_vtk)
    save_vtk_name("./output/after_curv");
}


// Mark the vertices based on cell membership.
void vd_sim::adapt_mark_bound_min() {

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");

  apf::ModelEntity* mdl;

  double ref_len = get_adapt_ln();
  std::cout << "ref_len " << ref_len << " bound_len " << bound_len << std::endl;

  double e_weight = 1./ref_len;
  double b_weight = 1./bound_len;
  apf::MeshEntity* e;

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it))) {
    apf::setScalar(field_step, e, 0, e_weight);
  }
  m->end(it);
  for(int i = 0; i < cells_ref.size(); i++) {
    int dim_curr = cells_ref.at(i).first;
    int id_curr = cells_ref.at(i).second - 1;
    std::vector<apf::MeshEntity*>* v_set = &e_list.e.at(dim_curr).
                                                        at(id_curr).at(0);
    for(int j = 0; j < v_set->size(); j++) {
      apf::setScalar(field_step, v_set->at(j), 0, b_weight);
    }
  }
}

// Mark the vertices based on the normal distances to the 2strata.
void vd_sim::adapt_mark_2c_dist(std::map<apf::MeshEntity*, double> &v_dist,
std::vector<double> &dist_c3) {
/*
  double ref_len = get_adapt_ln();

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");

  // The tetrahedra burning front:
  std::vector<std::vector<apf::MeshEntity*> > front(c_base->get_sz(3), 
                                  std::vector<apf::MeshEntity*>(0));
  // The minimum distance from any bounding 2strata, for each grain.
  dist_c3.clear();
  dist_c3.resize(c_base->get_sz(3));

  // Candidate tets, the adjacencies of the current front. They are added if
  // unburnt and not already in the front.
  std::vector<apf::MeshEntity*> tet_c(0);
  std::vector<apf::MeshEntity*> tris(0);
  // Once distances over all possible tets are exhausted, burn the vertex.
  std::map<apf::MeshEntity*, bool> v_burn{};
  // The stratum id of the minimum distance boundary. -1 for 1 and 0 strata,
  // id for 2 strata.
  std::map<apf::MeshEntity*, bool> v_2c{};

  v_dist.clear();
  // The approximate normal distance for a given vertex from the closest 2
  // stratum.
  std::map<apf::MeshEntity*, double> t_dist{};
  // If all vertex adjacencies of a tet are burnt, burn the tet.
  std::map<apf::MeshEntity*, bool> t_burn{};
  // The unburnt vertex associated with the tet.
  std::map<apf::MeshEntity*, apf::MeshEntity*> tv_map{};
  std::map<apf::MeshEntity*, int> tv_id{};

  // The percent difference of the maximum distance being burned to the  
  // minimum distance being burned in the current loop.  
  double dist_tol = 0.1;

  bool v_per_3c = false;
  while(!v_per_3c) {
    v_per_3c = true;
    for(int c3 = 0; c3 < e_list.e.at(3).size(); c3++) {
      if(e_list.e.at(3).at(c3).at(3).size() > 0 and
         e_list.e.at(3).at(c3).at(0).size() == 0) {
        v_per_3c = false;
        c3 = e_list.e.at(3).size();
      }
    }
    if(!v_per_3c) {
      adapt_split();
      save_vtk_name("./output/after_split");
    }
    else {
      std::cout << "All c3 have a vertex." << std::endl;
    }
  }

  // Maximum size of entities is the number of tets per grain.
  for(int c3 = 0; c3 < e_list.e.at(3).size(); c3++) {
    front.at(c3).reserve(e_list.e.at(3).at(c3).at(3).size());
  }

  int nbr = 0;
  for(int c2 = 0; c2 < e_list.e.at(2).size(); c2++) {
    if(e_list.e.at(2).at(c2).at(2).size() > 0) {
      nbr = nbr + e_list.e.at(2).at(c2).at(2).size();
    }
  }
  tris.reserve(nbr);
  for(int c2 = 0; c2 < e_list.e.at(2).size(); c2++) {
    if(e_list.e.at(2).at(c2).at(2).size() > 0) {
      for(int i = 0; i < e_list.e.at(2).at(c2).at(2).size(); i++) {
        apf::MeshEntity* t_curr = e_list.e.at(2).at(c2).at(2).at(i);
        tris.push_back(t_curr);
      }
    }
  }

  vd_set_up(m, &tris, &tet_c);
  for(int i = 0; i < tet_c.size(); i++) {
    apf::ModelEntity* mdl = m->toModel(tet_c.at(i));
    int c3 = m->getModelTag(mdl)-1;
    front.at(c3).push_back(tet_c.at(i));
  }

  // Burn the boundary vertices
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;
  while(vert = m->iterate(it)) {
    v_dist[vert] = sim_len;
  }
  m->end(it);

  for(int d = 0; d < 3; d++) {
    for(int c = 0; c < e_list.e.at(d).size(); c++) {
      if(e_list.e.at(d).at(c).at(0).size() > 0) {
        for(int i = 0; i < e_list.e.at(d).at(c).at(0).size(); i++) {
          apf::MeshEntity* v_curr = e_list.e.at(d).at(c).at(0).at(i);
          v_dist[v_curr] = 0;
          v_burn[v_curr] = true;
          if(d == 2)
            v_2c[v_curr] = c;
          else
            v_2c[v_curr] = -1;
        }
      }
    }
  }

  apf::Downward d_v;
  apf::Downward d_e;
  apf::Downward d_t;
  apf::Downward d_t2;
  apf::Up up;

  // Vertex to adjacent triangles:
  int lookup_tet_surf [4][3] = {{0,1,3},{0,2,1},{0,3,2},{1,2,3}};
  // Vertex to surface across
  int lookup_tet_surf_x [4] = {2, 3, 1, 0};

  for(int c3 = 0; c3 < front.size(); c3++) {
    int burnt = 0;
    for(int i = 0; i < front.at(c3).size() - burnt; i++) {
      apf::MeshEntity* t_curr = front.at(c3).at(i);
      m->getDownward(t_curr, 0, d_v);
      bool found = false;
      for(int j = 0; j < 4; j++) {
        if(!v_burn[d_v[j]]) {
          assert(!found);
          found = true;
          tv_map[t_curr] = d_v[j];
          tv_id[t_curr] = j;
        }
      }
      if(!found) {
        // Last unburnt tet:
        int i_ub = front.at(c3).size() - burnt - 1;
        assert(!t_burn[front.at(c3).at(i_ub)]);
        if(burnt > 0)
          assert(t_burn[front.at(c3).at(i_ub+1)]);

        front.at(c3).at(i) = front.at(c3).at(i_ub);
        front.at(c3).at(i_ub) = t_curr;
        t_burn[t_curr] = true;

        i = i-1;
        burnt = burnt + 1;
      }
    }
    front.at(c3).resize(front.at(c3).size() - burnt);
  }

  for(int c3 = 0; c3 < front.size(); c3++) {
    std::map<apf::MeshEntity*, apf::Vector3> t_pos{};
    std::map<apf::MeshEntity*, apf::Vector3> t_dir{};
    std::map<apf::MeshEntity*, apf::Vector3> e_pos{};
    std::map<apf::MeshEntity*, apf::Vector3> e_dir{};
    std::map<apf::MeshEntity*, apf::Vector3> v_pos{};

    std::map<apf::MeshEntity*, apf::MeshEntity*> e_t_map1{};
    std::map<apf::MeshEntity*, apf::MeshEntity*> e_t_map2{};

    ent_conn* edown = new ent_conn();
    c_base->get_conn(3, c3, edown);

    // Collect the vertex and edge adjacencies of tris bounding the 3strata,
    // and the triangle adjacencies of the edges bounding the 3strata.
    std::map<std::vector<apf::MeshEntity*> > t_v_adj{};
    std::vector<std::vector<apf::MeshEntity*> > t_e_adj{};
    for(int i = 0; i < edown->conn.size(); i++) {
      int c2 = edown->conn.at(i);
      int t_sz = e_list.e.at(2).at(c2).at(2).size();
      for(int j = 0; j < t_sz; j++) {
        apf::MeshEntity* t_curr = e_list.e.at(2).at(c2).at(2).at(j);
        m->getDownward(t_curr, 0, d_v);
        std::vector<apf::MeshEntity*> es(0);
        copy_ent_set(&es, d_v, 3);
        t_v_adj[t_curr] = es;

        m->getDownward(t_curr, 1, d_e);
        copy_ent_set(&es, d_e, 3);
        t_e_adj[t_curr] = es;
        for(int k = 0; k < 3; k++) {
          if(e_t_map1[es.at(k)] == NULL)
            e_t_map1[es.at(k)] = t_curr;
          else {
            assert(e_t_map2[es.at(k)] == NULL);
            e_t_map2[es.at(k)] = t_curr;
          }
        }

Collect the positions of the vertices, edges, tris
Collect the directions
For tris, direction that points into the 3strata
      }
    }

    delete edown;

    dist_c3.at(c3) = sim_len;
    while(front.at(c3).size() > 0) {
      // The end burnt number of tets are burnt and not to be considered.
      int burnt = 0;
      double dist_min = sim_len;
      for(int i = 0; i < front.at(c3).size() - burnt; i++) {
        apf::MeshEntity* t_curr = front.at(c3).at(i);
        assert(!t_burn[t_curr]);
        m->getDownward(t_curr, 0, d_v);
        bool found = false;
        for(int j = 0; j < 4; j++) {
          if(!v_burn[d_v[j]]) {
            assert(!found);
            found = true;
            tv_map[t_curr] = d_v[j];
            tv_id[t_curr] = j;
          }
        }
        assert(found);
        t_dist[t_curr] = approx_dist(m, v_dist, d_v, tv_map[t_curr]);
        if(t_dist[t_curr] < v_dist[tv_map[t_curr]]) {
          v_dist[tv_map[t_curr]] = t_dist[t_curr];
          if(v_dist[tv_map[t_curr]] < dist_min)
            dist_min = v_dist[tv_map[t_curr]];
        }
      }
      int added = 0;
      double dist_th = dist_min*(1+dist_tol);
      for(int i = 0; i < front.at(c3).size() - burnt - added; i++) {
        apf::MeshEntity* t_curr = front.at(c3).at(i);
        if(v_dist[tv_map[t_curr]] < dist_th) {
          v_burn[tv_map[t_curr]] = true;
        }
      }

      for(int i = 0; i < front.at(c3).size() - burnt - added; i++) {
        apf::MeshEntity* t_curr = front.at(c3).at(i);
        assert(tv_id[t_curr] > -1);
        assert(tv_map[t_curr] != NULL);

        if(v_dist[tv_map[t_curr]] < dist_th) {
          m->getDownward(t_curr, 2, d_t);
          int v_id = tv_id[t_curr];
          for(int j = 0; j < 3; j++) {
            int t_id = lookup_tet_surf[v_id][j];
            apf::MeshEntity* tri_curr = d_t[t_id];

            m->getUp(tri_curr, up);
            if(up.n == 2) {
              apf::MeshEntity* t_next;
              if(up.e[0] == t_curr)
                t_next = up.e[1];
              else
                t_next = up.e[0];
              if(!t_burn[t_next]) {
                int t1 = findIn(&front.at(c3), front.at(c3).size(), t_next);
                if(t1 == -1) {
                  m->getDownward(t_next, 0, d_v);
                  bool found1 = false;
                  bool found2 = false;
                  for(int j = 0; j < 4; j++) {
                    if(!v_burn[d_v[j]]) {
                      if(found1)
                        found2 = true;
                      else
                        found1 = true;
                    }
                  }
                  assert(!found2);
                  if(found1) {
                    if(burnt > 0) {
                      int i_b = front.at(c3).size() - burnt;
                      assert(!t_burn[front.at(c3).at(i_b-1)]);
                      assert(t_burn[front.at(c3).at(i_b)]);
                      front.at(c3).push_back(front.at(c3).at(i_b));
                      front.at(c3).at(i_b) = t_next;
                      //m->getDownward(t_next, 2, d_t2);
                      //int t2 = findIn(d_t2, 4, tri_curr);
                      //assert(t2 > -1);
                      //tv_map[t_next] = d_v[lookup_tet_surf_x[t2]];
                      //t_dist[t_next] = approx_dist(m, v_dist, d_v, 
                      //                              tv_map[t_curr]);
                    }
                    else
                      front.at(c3).push_back(t_next);

                    tv_id[t_next] = -1;
                    tv_map[t_next] = NULL;

                    added = added + 1;
                  }
                  else {
                    t_burn[t_next] = true;
                  }
                }
              }
            }
          }
          t_burn[t_curr] = true;
          int i_ub = front.at(c3).size() - burnt - added - 1;
          int i_ad = front.at(c3).size() - burnt - 1;
          if(i_ub != i)
            assert(!t_burn[front.at(c3).at(i_ub)]);
          if(added > 0)
            assert(!t_burn[front.at(c3).at(i_ub+1)]);
          if(burnt > 0)
            assert(t_burn[front.at(c3).at(i_ub+added+1)]);

          front.at(c3).at(i) = front.at(c3).at(i_ub);
          front.at(c3).at(i_ub) = front.at(c3).at(i_ad);
          front.at(c3).at(i_ad) = t_curr;
          burnt = burnt + 1;
          i = i-1;
        }
      }
      front.at(c3).resize(front.at(c3).size() - burnt);
    }
    // Actually, going over the c3 edges, check for vertices belonging to c3
    // Mark vertex having less distance as false.
    // Next, going over the vertices, 

    std::map<apf::MeshEntity*, bool> local_max{};
    for(int i = 0; i < e_list.e.at(3).at(c3).at(0).size(); i++) {
      local_max[e_list.e.at(3).at(c3).at(0).at(i)] = true;
    }

    for(int i = 0; i < e_list.e.at(3).at(c3).at(1).size(); i++) {
      apf::MeshEntity* e_curr = e_list.e.at(3).at(c3).at(1).at(i);
      m->getDownward(e_curr, 0, d_v);
      apf::ModelEntity* mdl = m->toModel(d_v[0]);
      int d_0 = m->getModelType(mdl);
      mdl = m->toModel(d_v[1]);
      int d_1 = m->getModelType(mdl);
      if(d_0 == 3 and d_1 == 3) {
        if(v_dist[d_v[0]] < v_dist[d_v[1]]) {
          local_max[d_v[0]] = false;
        }
        else if(v_dist[d_v[1]] < v_dist[d_v[0]]) {
          local_max[d_v[1]] = false;
        }
      }
    }
    for(int i = 0; i < e_list.e.at(3).at(c3).at(0).size(); i++) {
      apf::MeshEntity* v_curr = e_list.e.at(3).at(c3).at(0).at(i);
      if(local_max[v_curr]) {
        if(dist_c3.at(c3) > v_dist[v_curr])
          dist_c3.at(c3) = v_dist[v_curr];
      }
    }
  }
  double EdgeMaxBndry = get_lay_thk();
  double min_c3 = sim_len;
  for(int i = 0; i < dist_c3.size(); i++) {
    if(dist_c3.at(i) < min_c3)
      min_c3 = dist_c3.at(i);
  }
  if(min_c3 > EdgeMaxBndry)
    EdgeMaxBndry = min_c3;

  it = m->begin(0);
  while(vert = m->iterate(it)) {
    if(m->getModelType(m->toModel(vert)) < 3)
      apf::setScalar(field_step, vert, 0, 1/EdgeMaxBndry);
    else
      apf::setScalar(field_step, vert, 0, 1/v_dist[vert]);
  }
  m->end(it);
*/
/*
actually, burn by distance, or larger tets will cover more volume and by current algorithm lock distances of certain vertices
the front should be tets
for each grain keep a list of front tets
burn by layers of distance. find the minimum distance in the front tets
pick say 110% threshold of that distance. burn the tets below the threshold
add the unburned tets adjacent to the burning tets into the front
discard the burned tets by a merge sort fashion (move them to the end of the list, resize the list)
The layer thickness should be half of the minimum, roughly the length of 2 edges, for each individual grain

So, rather than calculating the small cell size, just run this burning algorithm and calculating the layer widths for each grain, grade them individually

I think this will also detect local constrictions. If a tet doesn't have any burnable adjacent tets, 
*/
}

void vd_sim::adapt_iso() {
  f_calc.del_tri_fields(m);

  assert(ad_type < ADAPT_TYPE::END);
  double adapt_ln = get_adapt_ln();

  double m_len =  min_cell.at(0)/2;

  if(ad_type < ADAPT_TYPE::ADAPT_STEP) {
    //Ref_Bound sf(m, adapt_ln);
    Linear sf(m, adapt_ln);
    //Linear_0cell sf(m, c_base, adapt_ln);
    ma::Input* in = ma::configure(m, &sf);

    std::vector<apf::MeshEntity*> tets(0);
    double valid = get_low_q(m, in->sizeField, tets, q_th);
    double good = get_good_q();

    good = std::max(valid*1.1, good);
    //double good = 0.1;

    in->shouldRunPreZoltan = true;
    in->shouldRunMidParma = true;
    in->shouldRunPostParma = true;
    in->shouldRefineLayer = true;

    in->shouldFixShape = ad_fix_shape;
    in->validQuality = valid;
    in->goodQuality = good;
    //in->goodQuality = 0.2;
    double t0 = PCU_Time();
    ma::adapt(in);
    double t1 = PCU_Time();
    dt_adapt = dt_adapt + (t1 - t0);
  }
  else {
    adapt_mark_step();
    Step sf(m, adapt_ln);
    //Linear_0cell sf(m, c_base, adapt_ln);
    ma::Input* in = ma::configure(m, &sf);

    std::vector<apf::MeshEntity*> tets(0);
    double valid = get_low_q(m, in->sizeField, tets, q_th);
    double good = get_good_q();

    good = std::max(valid*1.1, good);
    //double good = 0.1;

    in->shouldRunPreZoltan = true;
    in->shouldRunMidParma = true;
    in->shouldRunPostParma = true;
    in->shouldRefineLayer = true;

    in->shouldFixShape = ad_fix_shape;
    in->validQuality = valid;
    in->goodQuality = good;
    //in->goodQuality = 0.2;
    double t0 = PCU_Time();
    ma::adapt(in);
    double t1 = PCU_Time();
    dt_adapt = dt_adapt + (t1 - t0);
    apf::Field* a_f = m->findField("adapt_step");
    if(a_f)
      apf::destroyField(a_f);

  }

  m->acceptChanges();
  m->verify();
  f_calc.vd_att_fields(m);

  //save_vtk();

}

double vd_sim::get_lay_thk() {
  return get_adapt_ln();
}

double vd_sim::get_adapt_ln() {
/*
  if(ad_type == ADAPT_TYPE::ADAPT_1CELL) {
    std::vector<double> c1_len(0);
    c1_len.reserve(c_base->get_sz(1));
    for(int i = 0; i < c_base->get_sz(1); i++) {
      if(!c_base->is_free(1, i)) {
        c1_len.push_back(sm_map.at(0)[i+1]);
      }
    }
    return getMedianEntSize_hist(c1_len, 100)/adapt_len/ad_param*16;
    //return adapt_len*avg_cell.at(0)/ad_param/ma::getAverageEdgeLength(m);
    //return adapt_len*ad_param*avg_cell.at(0);
    //return avg_cell.at(0)/adapt_len/ad_param;
  }
  else {
    assert(ad_type == ADAPT_TYPE::ADAPT_3CELL or
           ad_type == ADAPT_TYPE::ADAPT_STEP or
           ad_type == ADAPT_TYPE::ADAPT_STEP_1CELL);

    //return adapt_len*std::cbrt(avg_cell.at(2))
    //                  /ad_param/ma::getAverageEdgeLength(m);
    //return adapt_len*ad_param*std::cbrt(avg_cell.at(2));
    //return avg_cell.at(2)/adapt_len/ad_param;
    std::vector<double> c3_vol(0);
    c3_vol.reserve(c_base->get_sz(3));
    for(int i = 0; i < c_base->get_sz(3); i++) {
      if(!c_base->is_free(3, i)) {
        c3_vol.push_back(sm_map.at(2)[i+1]);
      }
    }
    return getMedianEntSize_hist(c3_vol, 100)/adapt_len/ad_param*16;
  }
*/
  return len_topo/adapt_len/ad_param;
}

// Given a cell, coarsen the edges emanating from bounding vertices within a
// sphere with a radius twice the equivalent radius of the given cell. 
// Given a cell, imagine a sphere with 2* the equivalent radius of the cell 
// and centered at the midpoint.
// Coarsen all edges with both their end points within the sphere as a 
// preconditioning before collapse. 
// This could affect the evolution of adjacent cells if the coarsening changes
// the mesh considerably. Still, this should be uncommon as small cells tend to
// be flat and have flat bounding strata.
void vd_sim::adapt_coarsen_cell(int cell_dim, int cell_id) {
  e_list.refresh();

  int rad_curr = sm_map.at(cell_dim-1)[cell_id];
  apf::Vector3 midpoint(0,0,0);
  apf::Vector3 disp(0,0,0);
  midpoint = get_cell_midpoint(cell_dim, cell_id);

  // Collect all bounding 0cells. For these 0cells, collect all the bounded 3cells
  // Collect the bounding cells of these 3cells. 
  // The vertices outside the sphere or belonging to the non-collapsing bounding 
  // cells are the end vertices of the coarsening cavity.
  std::vector<std::map<int, bool> > end_c(3, std::map<int,bool>{});
  std::vector<std::map<int, bool> > adj_c(3, std::map<int,bool>{});

  std::map<apf::MeshEntity*, bool> burned{};

  c_base->get_end_gmi(cell_dim, cell_id, end_c, adj_c);

  apf::Downward d_v;
  std::vector<apf::MeshEntity*> verts(0);
  std::vector<apf::MeshEntity*> edges(0);
  std::vector<apf::MeshEntity*> ents(0);

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  Step_ns sf(m, get_adapt_ln());
  ma::Input* in = ma::configure(m, &sf);

  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  repl_sz_field(in, m, MA_SIZE_TYPE::EDGE_COL);
  ModelEdgeCollapse* ref_0c = (ModelEdgeCollapse*) in->sizeField;

  ents = e_list.e.at(cell_dim).at(cell_id - 1).at(cell_dim);
  vd_set_down(m, &ents, &edges, cell_dim - 1);
  vd_set_down(m, &edges, &verts);

  for(int i = 0; i < edges.size(); i++) {
    ref_0c->coarse_map[edges.at(i)] = true;
    burned[edges.at(i)] = true;
  }
  for(int i = 0; i < verts.size(); i++) {
    burned[verts.at(i)] = true;
  }

  bool all_burned = false;
  while(verts.size() > 0) {
    all_burned = true;
    vd_set_up(m, &verts, &edges);
    verts.clear();
    verts.reserve(edges.size());

    for(int i = 0; i < edges.size(); i++) {
      if(!burned[edges.at(i)]) {
        m->getDownward(edges.at(i), 0, d_v);
        if(burned[d_v[0]] and !burned[d_v[1]]) {
          apf::MeshEntity* v_next = d_v[1];
          apf::ModelEntity* mdl = m->toModel(v_next);
          int d_next = m->getModelType(mdl);
          int id_next = m->getModelTag(mdl);
          if(d_next == 3 or !end_c.at(d_next)[id_next]) {
            m->getPoint(v_next, 0, disp);
            disp = disp - midpoint;
            if(disp.getLength() < rad_curr) {
              verts.push_back(v_next);
            }
          }
          burned[v_next] = true;
        }
        else if(burned[d_v[1]] and !burned[d_v[0]]) {
          apf::MeshEntity* v_next = d_v[0];
          apf::ModelEntity* mdl = m->toModel(v_next);
          int d_next = m->getModelType(mdl);
          int id_next = m->getModelTag(mdl);
          if(d_next == 3 or !end_c.at(d_next)[id_next]) {
            m->getPoint(v_next, 0, disp);
            disp = disp - midpoint;
            if(disp.getLength() < rad_curr) {
              verts.push_back(v_next);
            }
          }
          burned[v_next] = true;
        }
      }
      ref_0c->coarse_map[edges.at(i)] = true;
    }
  }

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = false;

  in->shouldCoarsen = true;
  in->shouldCoarsenLayer = true;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);


  if(sub_vtk)
    save_vtk_name("./output/adapt_coarsen");

  apf::destroyField(field_step);
  e_list.refresh();
}

// Given the low quality tets, set the edge belonging to the lowest dimensional 
// cell to be coarsened. Only coarsen.
void vd_sim::adapt_prob_edge(std::vector<apf::MeshEntity*> &tet) {

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  Step_ns sf(m, get_adapt_ln());
  ma::Input* in = ma::configure(m, &sf);

  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  repl_sz_field(in, m, MA_SIZE_TYPE::EDGE_COL);
  ModelEdgeCollapse* ref_0c = (ModelEdgeCollapse*) in->sizeField;

  for(int i = 0; i < tet.size(); i++) {
    int e1 = 0;
    int d = 3;
    apf::Downward d_e;
    m->getDownward(tet.at(i), 1, d_e);
    for(int j = 0; j < 6; j++) {
      int type = m->getModelType(m->toModel(d_e[j]));
      if(type < d) {
        e1 = j;
        d = type;
      }
    }
    ref_0c->coarse_map[d_e[e1]] = true;
  }

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = false;

  in->shouldCoarsen = true;
  in->shouldCoarsenLayer = true;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  if(sub_vtk)
    save_vtk_name("./output/adapt_prob");

}

// If adapt_prob_edge doesn't reduce the number of low quality tets, try marking
// all edges adjacent to the vertices to be coarsened.
void vd_sim::adapt_prob_edge_all(std::vector<apf::MeshEntity*> &tet) {

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  Step_ns sf(m, get_adapt_ln());
  ma::Input* in = ma::configure(m, &sf);

  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  delete in->sizeField;
  ModelEdgeCollapse* ref_0c = new ModelEdgeCollapse(m);
  in->sizeField = ref_0c;

  std::vector<apf::MeshEntity*> verts(0);
  std::vector<apf::MeshEntity*> edges(0);

  vd_set_down(m, &tet, &verts, 3);
  vd_set_up(m, &verts, &edges);
  for(int i = 0; i < edges.size(); i++) {
    ref_0c->coarse_map[edges.at(i)] = true;
  }

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = false;

  in->shouldCoarsen = true;
  in->shouldCoarsenLayer = true;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  if(sub_vtk)
    save_vtk_name("./output/adapt_prob_all");

}

// If adapt_prob_edge doesn't reduce the number of low quality tets, try marking
// all edges adjacent to the vertices belonging to the same cell as the edge 
// belonging to the tet with the lowest dim cell membership.
void vd_sim::adapt_prob_edge_low(std::vector<apf::MeshEntity*> &tet) {

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  Step_ns sf(m, get_adapt_ln());
  ma::Input* in = ma::configure(m, &sf);

  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  delete in->sizeField;
  ModelEdgeCollapse* ref_0c = new ModelEdgeCollapse(m);
  in->sizeField = ref_0c;

  std::vector<apf::MeshEntity*> verts(0);
  std::vector<apf::MeshEntity*> edges(0);

  for(int i = 0; i < tet.size(); i++) {
    int e1 = 0;
    int d = 3;
    apf::Downward d_e;
    m->getDownward(tet.at(i), 1, d_e);
    for(int j = 0; j < 6; j++) {
      int type = m->getModelType(m->toModel(d_e[j]));
      if(type < d) {
        e1 = j;
        d = type;
      }
    }
    //ref_0c->coarse_map[d_e[e1]] = true;
    apf::ModelEntity* mdl_curr = m->toModel(d_e[e1]);
    vd_set_down(m, d_e[e1], &verts);
    vd_set_up(m, &verts, &edges);
    for(int j = 0; j < edges.size(); j++) {
      if(m->toModel(edges.at(j)) == mdl_curr)
        ref_0c->coarse_map[edges.at(j)] = true;
    }
  }

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = false;

  in->shouldCoarsen = true;
  in->shouldCoarsenLayer = true;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  if(sub_vtk)
    save_vtk_name("./output/adapt_prob_low");

}

void vd_sim::adapt_prob_edge_sp_col(std::vector<apf::MeshEntity*> &tet) {
  apf::Up up;
  apf::Downward d_s;
  apf::Downward d_e;
  apf::Downward d_e3;
  apf::Downward d_e4;
  apf::Downward d_v;

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  double adapt_ln = get_adapt_ln();
  double m_len =  min_cell.at(0)*2;
  double q_frac = 1.1;

  Step_ns sf(m, adapt_ln);
  ma::Input* in = ma::configure(m, &sf);

//old solution and mix of this one
  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  repl_sz_field(in, m, MA_SIZE_TYPE::EDGE_SPLIT);

  std::vector<apf::MeshEntity*> tets_v(0);
  std::vector<apf::MeshEntity*> verts(0);
  std::vector<apf::MeshEntity*> edges(0);
  vd_set_down(m, &tet, &verts, 3);

  ModelEdgeSplit* ref_0c = (ModelEdgeSplit*) in->sizeField;

  //for(int i = 0; i < edges.size(); i++)
  //  ref_0c->split_map[edges.at(i)] = true;
  if(sub_vtk)
    vd_save_vtk_ent(m, &tet, "./output/adapt_prob_sp_ents");

  std::map<apf::MeshEntity*, apf::MeshEntity*> tri_tet_map{};
  std::map<apf::MeshEntity*, int> tet_neigh_nbr{};
  for(int i = 0; i < tet.size(); i++) {
    m->getDownward(tet.at(i), 2, d_s);
    for(int j = 0; j < 4; j++) {
      if(!tri_tet_map[d_s[j]]) {
        tri_tet_map[d_s[j]] = tet.at(i);
      }
      else {
        tet_neigh_nbr[tri_tet_map[d_s[j]]] = 
                                 tet_neigh_nbr[tri_tet_map[d_s[j]]] + 1;
        tet_neigh_nbr[tet.at(i)] = tet_neigh_nbr[tet.at(i)] + 1;
      }
    }
  }
  int zero_tet = 0;
  int sing_tet = 0;
  for(int i = 0; i < tet.size() - zero_tet; i++) {
    if(tet_neigh_nbr[tet.at(i)] == 0) {
      apf::MeshEntity* temp = tet.at(i);
      int subs_id = tet.size() - zero_tet - 1;
      tet.at(i) = tet.at(subs_id);
      tet.at(subs_id) = temp;
      zero_tet = zero_tet + 1;
      i = i - 1;
    }
  }

  for(int i = 0; i < tet.size() - zero_tet - sing_tet; i++) {
    if(tet_neigh_nbr[tet.at(i)] == 1) {
      apf::MeshEntity* temp = tet.at(i);
      int subs_id = tet.size() - zero_tet - sing_tet - 1;
      tet.at(i) = tet.at(subs_id);
      tet.at(subs_id) = temp;
      sing_tet = sing_tet + 1;
      i = i - 1;
    }
  }

  // Sort edges by length. Starting from the shortest edge, try until one the 
  // first edge can be collapsed.
  std::map<int, bool> map_col_id{};
  std::map<apf::MeshEntity*, bool> map_col{};
  std::map<apf::MeshEntity*, apf::MeshEntity*> map_sur{};

  std::pair<bool, bool> valid(false, false);
  for(int i = tet.size() - zero_tet; i < tet.size(); i++) {
    if(map_col_id[i] ) {
    }
    else {
      valid = vd_col_edge_tet(m, c_base, tet.at(i), map_sur, map_col, &f_calc, true);
      if(valid.second)
        break;
      if(valid.first) {
        map_col_id[i] = true;
        for(int j = 0; j < i; j++) {
          if(map_sur[tet.at(j)]) {
            tet.at(j) = map_sur[tet.at(j)];
          }
          if(map_col[tet.at(j)]) {
            map_col_id[j] = true;
          }
        }
        map_sur.clear();
        map_col.clear();
      }
    }
  }
  if(valid.second)
    tet.clear();
  else {
    tet.resize(tet.size() - zero_tet);
    valid = std::pair<bool, bool>(false, false);

    for(int i = 0; i < tet.size(); i++) {
      if(map_col_id[i] ) {

      }
      else {
        valid = vd_col_edge_tet(m, c_base, tet.at(i), map_sur, map_col, &f_calc, true);
        if(valid.second)
          break;
        if(valid.first) {
          map_col_id[i] = true;
          for(int j = 0; j < i; j++) {
            if(map_sur[tet.at(j)]) {
              tet.at(j) = map_sur[tet.at(j)];
            }
            if(map_col[tet.at(j)]) {
              map_col_id[j] = true;
            }
          }
          map_sur.clear();
          map_col.clear();
          break;
        }
      }
    }
  }

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = true;

  in->shouldCoarsen = false;
  in->shouldCoarsenLayer = false;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  if(sub_vtk)
    save_vtk_name("./output/adapt_prob_sp");
}
// This preconditions problematic low quality tetrahedra that remain un-fixed by
// the meshadapt.
// Splitting the edge e1 across the boundary tris doesn't seem to fix things
// as the tet adjacencies of the triangles bounded by the edge can't be merged
// with the interior tets obtained after the split. Instead, find the edge
// joining the boundary tris. Also for the tet adjacencies of the other two
// tris. On these tets, find the edges across the edge e1. Splitting these
// three edges should yield two couples of tets, one of each pair low quality
// both can be improved by a edge flip operation by meshadapt.
// This is applied first before resorting to other split operations.
void vd_sim::adapt_prob_edge_sp_bound(std::vector<apf::MeshEntity*> &tet) {
  apf::Up up;
  apf::Downward d_s;
  apf::Downward d_e;
  apf::Downward d_e3;
  apf::Downward d_e4;
  apf::Downward d_v;

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  double adapt_ln = get_adapt_ln();
  double m_len =  min_cell.at(0)*2;
  double q_frac = 1.1;

  Step_ns sf(m, adapt_ln);
  ma::Input* in = ma::configure(m, &sf);

//old solution and mix of this one
  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  repl_sz_field(in, m, MA_SIZE_TYPE::EDGE_SPLIT);

  std::vector<apf::MeshEntity*> tets_v(0);
  std::vector<apf::MeshEntity*> verts(0);
  std::vector<apf::MeshEntity*> edges(0);
  vd_set_down(m, &tet, &verts, 3);

  ModelEdgeSplit* ref_0c = (ModelEdgeSplit*) in->sizeField;

  //for(int i = 0; i < edges.size(); i++)
  //  ref_0c->split_map[edges.at(i)] = true;
  if(sub_vtk)
    vd_save_vtk_ent(m, &tet, "./output/adapt_prob_sp_ents");

/*
Check for tets with a single low quality tet neighbor.

*/

/*
check triangles first
If no triangles, check edges
The central idea is to find two(three) interior triangles such that the line passing through middle point of the edge(triangle) across and the vertex on the other tets for all triangles intersect the triangle. 

This might be too strong a condition... It might be mesh adapt moves the vertices such that even without satisfying this condition, low q tet is resolved.
So first try finding two interior triangles.

Depending on the number of boundary entity found:
4 tris -> nothing, collapsing
3 tris -> split the tri across
2 tris -> find the edge joining them. Find the other triangle tets. For 3 tets, split edge across the joining edge.
1 tri -> try to find two triangles, s.t. the line between the middle point of the edge not bounding them on the tet and the vertex across the triangle on the tets sharing the triangles both intersect their respective triangles. If found, find the edge joining them. Find the other triangle tets. For 3 tets, split edge across the joining edge.
          If not found, check the lines between the middle point of the bounding 
          triangle and the vertices on the tets sharing their respective tris.
          If they all satisfy, split the triangle.

If the above is not satisfied:
Some of the edges cannot be flipped.
6 edges -> nothing, likely will disappear.
3 edges -> 3 configurations. Same tri, same vertex, trail.

If the above is not satisfied, skip.
*/
  //m->getDownward(tet.at(0), 1, d_e);
  //vd_col_edge(m, c_base, d_e[0], &f_calc);
  //tet.clear();

  for(int i = 0; i < tet.size(); i++) {
    int d = 3;
    m->getDownward(tet.at(i), 2, d_s);
    m->getDownward(tet.at(i), 1, d_e);
    int tri_state = 0;
/*
    // Check if the given tet has an exterior 2-stratum vertex bounding an 
    // exterior 2-stratum triangle. 
    if(apf::MeshEntity* vert_b = get_tet_bound_vert(m, tet.at(i))) {
    // Check if the tetrahedron lies between two boundary triangles, one belonging
    // to the exterior. If so, split the boundary vertex bounding the exterior
    // triangle. TODO this is very prone to crashes if multiple tetrahedra 
    //if(apf::MeshEntity* tri_2c = get_tet_bound_tri(m, tet.at(i))) {
    //  apf::MeshEntity* vert_b = get_tri_vert_bound(m, tri_2c);
      std::vector<apf::MeshEntity*> edge_b(0);
      std::vector<apf::MeshEntity*> tri_b(0);
      std::vector<apf::MeshEntity*> tet_b(0);

      vd_set_up(m, vert_b, &edge_b);
      vd_set_up(m, &edge_b, &tri_b);
      vd_set_up(m, &tri_b, &tet_b);

      for(int j = 0; j < tet_b.size(); j++) {
        if(tet_b.at(j) != tet.at(i)) {
          int j1 = findIn(&tet, tet.size(), tet_b.at(j));
          if(j1 > i) {
            tet.at(j1) = tet.back();
            tet.resize(tet.size() - 1);
          }
        }
      }
      vd_split_tet_vert(m, c_base, vert_b, NULL);
    }
    else {
*/
      for(int j = 0; j < 4; j++) {
        apf::ModelEntity* mdl = m->toModel(d_s[j]);
        int type_e = m->getModelType(mdl);
        if(type_e == 2)
          tri_state = tri_state + tri_split_shift[j];
      }
      int e1 = -1;
      if(tri_split_type[tri_state] == 0) {
        // Also check for maximum edge length.
        double l_max = 0;
        int j_max = -1;
        for(int j = 0; j < 6; j++) {

          double l_curr = vd_meas_ent(m, d_e[j]);
          if(l_curr > l_max) {
            l_max = l_curr;
            j_max = j;
          }
        }
        if(l_max > m_len and
           j_max != -1) {
          e1 = tri_split_map[tri_state];
        }
      }
  /*
      // Single boundary triangle. From the remaining edges, determine the longest 
      // internal edge. 
      else if(tri_split_type[tri_state] == 2) {
        int t1 = tri_split_map[tri_state] - 6;
        //lookup_tt_e
        double l_max = 0;
        int j_max = -1;
        for(int j = 0; j < 3; j++) {
          e1 = lookup_t_tet_e_x[t1][j];
          apf::ModelEntity* mdl = m->toModel(d_e[e1]);
          int type_e = m->getModelType(mdl);
          if(type_e == 3) {
            double l_curr = vd_meas_ent(m, d_e[e1]);
            if(l_curr > l_max) {
              l_max = l_curr;
              j_max = e1;
            }
          }
        }
        if(l_max > min_ins_rad*2 and
           j_max != -1) {
          e1 = j_max;
        }
      }
      // Single internal triangle. Don't modify.
      else if(tri_split_type[tri_state] == 1) {

      }
      else {
        // Find the longest edge e2. Get the edge across(e1). 
        // Get the adjacent tris.
        // Find the edges e3 and e4 across e_i on the tets adjacent to the tris.
        // Split e2,e3,e4.
        double l_max = 0;
        int j_max = -1;
        for(int j = 0; j < 6; j++) {
          apf::ModelEntity* mdl = m->toModel(d_e[j]);
          int type_e = m->getModelType(mdl);
          if(type_e == 3) {
            double l_curr = vd_meas_ent(m, d_e[j]);
            if(l_curr > l_max) {
              l_max = l_curr;
              j_max = j;
            }
          }
        }
        if(l_max > min_ins_rad*2 and
           j_max != -1) {
          e1 = lookup_t_e_x_e[j_max];
          //e1 = j_max;
        }
      }
  */
      if(e1 != -1) {
        m->getDownward(tet.at(i), 0, d_v);
        int e2 = lookup_t_e_x_e[e1];
        int t1 = lookup_t_e_n_t[e1][0];
        int t2 = lookup_t_e_n_t[e1][1];

        // Vertices on the main tet that bound the edge to be split.
        int v1 = lookup_t_e_n_g[e1][0];
        int v2 = lookup_t_e_n_g[e1][1];

        apf::MeshEntity* vert_1 = d_v[v1];
        apf::MeshEntity* vert_2 = d_v[v2];

        // First tri:
        m->getUp(d_s[t1], up);
        assert(up.n == 2);
        apf::MeshEntity* tet_next = up.e[0];
        if(up.e[0] == tet.at(i))
          tet_next = up.e[1];
        m->getDownward(tet_next, 1, d_e3);
        int e1_next = findIn(d_e3, 6, d_e[e1]);
        assert(e1_next > -1);
        int e3 = lookup_t_e_x_e[e1_next];

        // Second tri:
        m->getUp(d_s[t2], up);
        assert(up.n == 2);
        tet_next = up.e[0];
        if(up.e[0] == tet.at(i))
          tet_next = up.e[1];
        m->getDownward(tet_next, 1, d_e4);
        e1_next = findIn(d_e4, 6, d_e[e1]);
        assert(e1_next > -1);
        int e4 = lookup_t_e_x_e[e1_next];
  /*
        // Do not split tets if vertices bounding e2 and e3 belong to strata not 
        // bounding each other.
        // Do not split tets where the improved tet quality is not larger by a  
        // multiple. This is to prevent explosive refinement of tets between two 
        // close boundaries... 

        apf::MeshEntity* v_oth3 = apf::getEdgeVertOppositeVert(m, d_e3[e3],
                                                                       vert_1);
        apf::MeshEntity* v_oth4 = apf::getEdgeVertOppositeVert(m, d_e4[e4], 
                                                                       vert_2);
        apf::ModelEntity* mdl_v = m->toModel(vert_1);
        int d_sp1 = m->getModelType(mdl_v);
        int c_sp1 = m->getModelTag(mdl_v);
        mdl_v = m->toModel(vert_2);
        int d_sp2 = m->getModelType(mdl_v);
        int c_sp2 = m->getModelTag(mdl_v);

        mdl_v = m->toModel(v_oth3);
        int d_sp3 = m->getModelType(mdl_v);
        int c_sp3 = m->getModelTag(mdl_v);
        mdl_v = m->toModel(v_oth4);
        int d_sp4 = m->getModelType(mdl_v);
        int c_sp4 = m->getModelTag(mdl_v);

        bool cell_ok = c_base->chk_conn_d_gmi(d_sp1, c_sp1, d_sp3, c_sp3) and 
                       c_base->chk_conn_d_gmi(d_sp2, c_sp2, d_sp4, c_sp4);

        bool cell_tri_ok = true;
        if(tri_split_type[tri_state] == 2) {
          mdl_v = m->toModel(d_e[e2]);
          int d_e_sp2 = m->getModelType(mdl_v);
          int c_e_sp2 = m->getModelTag(mdl_v);
          mdl_v = m->toModel(d_e3[e3]);
          int d_e_sp3 = m->getModelType(mdl_v);
          int c_e_sp3 = m->getModelTag(mdl_v);
          mdl_v = m->toModel(d_e4[e4]);
          int d_e_sp4 = m->getModelType(mdl_v);
          int c_e_sp4 = m->getModelTag(mdl_v);
          cell_tri_ok = !(d_e_sp2 == d_e_sp3 or
                          d_e_sp2 == d_e_sp4);
        }

        if(cell_tri_ok or cell_ok) {
  */
  /*
          double q_curr = measureTetQuality(m, ref_0c, tet.at(i));
          std::vector<apf::Vector3> v(4, apf::Vector3(0,0,0));
          std::vector<apf::Vector3> v_new(4, apf::Vector3(0,0,0));
          for(int j = 0; j < 4; j++) {
            m->getPoint(d_v[j], 0, v.at(j));
            v_new.at(j) = v.at(j);
          }
          double q_init = calc_q(m, v);

          apf::Vector3 pos_e2(0,0,0);
          apf::Vector3 pos_e3(0,0,0);
          apf::Vector3 pos_e4(0,0,0);

          pos_e2 = vd_get_pos(m, d_e[e2]);
          pos_e3 = vd_get_pos(m, d_e3[e3]);
          pos_e4 = vd_get_pos(m, d_e4[e4]);
          // Check the quality of the new tets generated by triangle t1.
          // Shift the position corresponding to v1 or v2 to the center of e1.
          // Looping over the other v(not v1/v2), shift the position to the center 
          // of e3. Check quality.
          bool q_ok = true;
          v_new.at(v2) = pos_e2;
          int j = 0;
          while(j < 3 and q_ok) {
            int v_next = (v2 + 1 + j) % 4;
            if(v_next != v1) {
              v_new.at(v_next) = pos_e3;
              double q_curr = calc_q(m, v_new);
              q_ok = (q_ok and q_curr > q_frac*q_init);
              v_new.at(v_next) = v.at(v_next);
            }
            j = j + 1;
          }
          v_new.at(v2) = v.at(v2);
          v_new.at(v1) = pos_e2;
          j = 0;
          while(j < 3 and q_ok) {
            int v_next = (v1 + 1 + j) % 4;
            v_new.at(v_next) = pos_e4;
            if(v_next != v2) {
              double q_curr = calc_q(m, v_new);
              q_ok = (q_ok and q_curr > q_frac*q_init);
              v_new.at(v_next) = v.at(v_next);
            }
            j = j + 1;
          }
          if(q_ok) {
  */
            ref_0c->split_map[d_e[e2]] = true;
            ref_0c->split_map[d_e3[e3]] = true;
            ref_0c->split_map[d_e4[e4]] = true;
  //        }
  //      }
    
    }
  }

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = true;

  in->shouldCoarsen = false;
  in->shouldCoarsenLayer = false;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  if(sub_vtk)
    save_vtk_name("./output/adapt_prob_sp");
}

// Tets with a single boundary triangle cannot be fixed by just splitting. This 
// fixes them by first splitting the triangle and for the resulting tets, marking 
// the interior edge of the triangular bipyramid bounding the three tets generated
// by splitting the problematic tet to be collapsed.
void vd_sim::adapt_prob_edge_col_bound_manual(std::vector<apf::MeshEntity*> &tet) {
  apf::Up up;
  apf::Downward d_s;
  apf::Downward d_e;
  apf::Downward d_e3;
  apf::Downward d_e4;
  apf::Downward d_v;

  int lookup_tet_x_s [4] = {3, 2, 0, 1};

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  double adapt_ln = get_adapt_ln();
  double m_len =  min_cell.at(0)/2;
  double q_frac = 1.1;

  Step_ns sf(m, adapt_ln);
  ma::Input* in = ma::configure(m, &sf);

  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  repl_sz_field(in, m, MA_SIZE_TYPE::EDGE_COL);
  ModelEdgeCollapse* ref_0c = (ModelEdgeCollapse*) in->sizeField;


  // The single boundary triangle tet lists, to be handled separately:
  std::vector<apf::MeshEntity*> split_tri(0);
  std::vector<apf::MeshEntity*> split_tet(0);
  std::vector<apf::MeshEntity*> split_vert(0);
  split_tet.reserve(tet.size());
  split_tri.reserve(tet.size());
  split_vert.reserve(tet.size());

  if(sub_vtk)
    vd_save_vtk_ent(m, &tet, "./output/adapt_prob_sp_ents");

  split_vert.reserve(tet.size());
  for(int i = 0; i < tet.size(); i++) {
    // Check if the given tet has an exterior 2-stratum vertex bounding an 
    // exterior 2-stratum triangle. 
    if(apf::MeshEntity* vert_b = get_tet_bound_vert(m, tet.at(i))) {
      split_vert.push_back(vert_b);
    }
  }
  {
    std::vector<apf::MeshEntity*> edges_v(0);
    std::vector<apf::MeshEntity*> tris_v(0);

    std::vector<apf::MeshEntity*>::iterator it;
    it = std::unique(split_vert.begin(), split_vert.end());
    split_vert.resize(std::distance(split_vert.begin(),it));
    vd_set_up(m, &split_vert, &edges_v);
    vd_set_up(m, &edges_v, &tris_v);
    vd_set_up(m, &tris_v, &split_tet);
  }
  vd_remove_set(&tet, &split_tet);
  for(int i = 0; i < split_tet.size(); i++) {
    vd_split_tet(m, c_base, split_tet.at(i), NULL);
  }
  split_tet.clear();

  for(int i = 0; i < tet.size(); i++) {
    int d = 3;
    m->getDownward(tet.at(i), 2, d_s);
    m->getDownward(tet.at(i), 1, d_e);
    int tri_state = 0;

    for(int j = 0; j < 4; j++) {
      apf::ModelEntity* mdl = m->toModel(d_s[j]);
      int type_e = m->getModelType(mdl);
      if(type_e == 2)
        tri_state = tri_state + tri_split_shift[j];
    }
    // Single boundary triangle. From the remaining edges, determine the longest 
    // internal edge. 
    if(tri_split_type[tri_state] == 2) {
      int t1 = tri_split_map[tri_state] - 6;
      //lookup_tt_e
      double l_max = 0;
      int j_max = -1;
      for(int j = 0; j < 3; j++) {
        int e1 = lookup_t_tet_e_x[t1][j];
        apf::ModelEntity* mdl = m->toModel(d_e[e1]);
        int type_e = m->getModelType(mdl);
        if(type_e == 3) {
          double l_curr = vd_meas_ent(m, d_e[e1]);
          if(l_curr > l_max) {
            l_max = l_curr;
            j_max = e1;
          }
        }
      }
      if(l_max > min_ins_rad*2 and
         j_max != -1) {
        m->getDownward(tet.at(i), 0, d_v);
        split_tet.push_back(tet.at(i));
        split_tri.push_back(d_s[t1]);
        int v1 = lookup_tet_x_s[t1];
        split_vert.push_back(d_v[v1]);
      }
    }
  }

  vd_bipy* bipy_split = new vd_bipy(m, c_base);
  std::map<apf::MeshEntity*, bool> split_tri_map{};
  for(int i = 0; i < split_tet.size(); i++) {
    if(!split_tri_map[split_tri.at(i)]) {
      bipy_split->load_tri(split_tri.at(i));
      bipy_split->split_bipy();

      apf::MeshEntity* e_col;
      apf::MeshEntity* v_sp = bipy_split->get_vert_ctr();
      vd_find_edge(m, split_vert.at(i), v_sp, &e_col);
      ref_0c->coarse_map[e_col] = true;
      split_tri_map[split_tri.at(i)] = true;
    }
  }

  delete bipy_split;

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = false;

  in->shouldCoarsen = true;
  in->shouldCoarsenLayer = true;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  if(sub_vtk)
    save_vtk_name("./output/adapt_prob");
}

// Splitting the edge e1 across the boundary tris doesn't seem to fix things
// as the tet adjacencies of the triangles bounded by the edge can't be merged
// with the interior tets obtained after the split. Instead, find the edge
// joining the boundary tris. Also for the tet adjacencies of the other two
// tris. On these tets, find the edges across the edge e1. Splitting these
// three edges should yield two couples of tets, one of each pair low quality
// both can be improved by a edge flip operation by meshadapt.
// This is applied first before resorting to other split operations.
void vd_sim::adapt_prob_edge_sp_bound_manual(std::vector<apf::MeshEntity*> &tet) {
  apf::Up up;
  apf::Downward d_s;
  apf::Downward d_e;
  apf::Downward d_e3;
  apf::Downward d_e4;
  apf::Downward d_v;

  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  double adapt_ln = get_adapt_ln();
  double m_len =  min_cell.at(0)*2;
  double q_frac = 1.1;

  Step_ns sf(m, adapt_ln);
  ma::Input* in = ma::configure(m, &sf);

//old solution and mix of this one
  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  repl_sz_field(in, m, MA_SIZE_TYPE::EDGE_SPLIT);
  ModelEdgeSplit* ref_0c = (ModelEdgeSplit*) in->sizeField;

  if(sub_vtk)
    vd_save_vtk_ent(m, &tet, "./output/adapt_prob_sp_ents");

  for(int i = 0; i < tet.size(); i++) {
    int d = 3;
    m->getDownward(tet.at(i), 2, d_s);
    m->getDownward(tet.at(i), 1, d_e);
    int tri_state = 0;


    if(apf::MeshEntity* vert_b = get_tet_bound_vert(m, tet.at(i))) {
    // Check if the tetrahedron lies between two boundary triangles, one belonging
    // to the exterior. If so, split the boundary vertex bounding the exterior
    // triangle. TODO this is very prone to crashes if multiple tetrahedra 
    //if(apf::MeshEntity* tri_2c = get_tet_bound_tri(m, tet.at(i))) {
    //  apf::MeshEntity* vert_b = get_tri_vert_bound(m, tri_2c);
      std::vector<apf::MeshEntity*> edge_b(0);
      std::vector<apf::MeshEntity*> tri_b(0);
      std::vector<apf::MeshEntity*> tet_b(0);

      vd_set_up(m, vert_b, &edge_b);
      vd_set_up(m, &edge_b, &tri_b);
      vd_set_up(m, &tri_b, &tet_b);

      for(int j = 0; j < tet_b.size(); j++) {
        if(tet_b.at(j) != tet.at(i)) {
          int j1 = findIn(&tet, tet.size(), tet_b.at(j));
          if(j1 > i) {
            tet.at(j1) = tet.back();
            tet.resize(tet.size() - 1);
          }
        }
      }
      vd_split_tet_vert(m, c_base, vert_b, NULL);
    }
    else {

      for(int j = 0; j < 4; j++) {
        apf::ModelEntity* mdl = m->toModel(d_s[j]);
        int type_e = m->getModelType(mdl);
        if(type_e == 2)
          tri_state = tri_state + tri_split_shift[j];
      }
      int e1 = -1;
      if(tri_split_type[tri_state] == 0) {
        e1 = tri_split_map[tri_state];
      }
  /*
      // Single boundary triangle. From the remaining edges, determine the longest 
      // internal edge.  
      // Currently handled by adapt_prob_col_bound_manual
      else if(tri_split_type[tri_state] == 2) {
        int t1 = tri_split_map[tri_state] - 6;
        //lookup_tt_e
        double l_max = 0;
        int j_max = -1;
        for(int j = 0; j < 3; j++) {
          e1 = lookup_t_tet_e_x[t1][j];
          apf::ModelEntity* mdl = m->toModel(d_e[e1]);
          int type_e = m->getModelType(mdl);
          if(type_e == 3) {
            double l_curr = vd_meas_ent(m, d_e[e1]);
            if(l_curr > l_max) {
              l_max = l_curr;
              j_max = e1;
            }
          }
        }
        if(l_max > m_len and
           j_max != -1) {
          e1 = j_max;
        }
      }
  */
      // Single internal triangle. Don't modify.
      else if(tri_split_type[tri_state] == 1) {

      }
      else {
        // Find the longest edge e2. Get the edge across(e1). 
        // Get the adjacent tris.
        // Find the edges e3 and e4 across e_i on the tets adjacent to the tris.
        // Split e2,e3,e4.
        double l_max = 0;
        int j_max = -1;
        for(int j = 0; j < 6; j++) {
          apf::ModelEntity* mdl = m->toModel(d_e[j]);
          int type_e = m->getModelType(mdl);
          if(type_e == 3) {
            double l_curr = vd_meas_ent(m, d_e[j]);
            if(l_curr > l_max) {
              l_max = l_curr;
              j_max = j;
            }
          }
        }
        if(l_max > m_len and
           j_max != -1) {
          e1 = lookup_t_e_x_e[j_max];
          //e1 = j_max;
        }
      }

      if(e1 != -1) {
        m->getDownward(tet.at(i), 0, d_v);
        int e2 = lookup_t_e_x_e[e1];
        int t1 = lookup_t_e_n_t[e1][0];
        int t2 = lookup_t_e_n_t[e1][1];

        // Vertices on the main tet that bound the edge to be split.
        int v1 = lookup_t_e_n_g[e1][0];
        int v2 = lookup_t_e_n_g[e1][1];

        apf::MeshEntity* vert_1 = d_v[v1];
        apf::MeshEntity* vert_2 = d_v[v2];

        // First tri:
        m->getUp(d_s[t1], up);
        assert(up.n == 2);
        apf::MeshEntity* tet_next = up.e[0];
        if(up.e[0] == tet.at(i))
          tet_next = up.e[1];
        m->getDownward(tet_next, 1, d_e3);
        int e1_next = findIn(d_e3, 6, d_e[e1]);
        assert(e1_next > -1);
        int e3 = lookup_t_e_x_e[e1_next];

        // Second tri:
        m->getUp(d_s[t2], up);
        assert(up.n == 2);
        tet_next = up.e[0];
        if(up.e[0] == tet.at(i))
          tet_next = up.e[1];
        m->getDownward(tet_next, 1, d_e4);
        e1_next = findIn(d_e4, 6, d_e[e1]);
        assert(e1_next > -1);
        int e4 = lookup_t_e_x_e[e1_next];

        // Do not split tets if vertices bounding e2 and e3 belong to strata not 
        // bounding each other.
        // Do not split tets where the improved tet quality is not larger by a  
        // multiple. This is to prevent explosive refinement of tets between two 
        // close boundaries... 

        apf::MeshEntity* v_oth3 = apf::getEdgeVertOppositeVert(m, d_e3[e3],
                                                                       vert_1);
        apf::MeshEntity* v_oth4 = apf::getEdgeVertOppositeVert(m, d_e4[e4], 
                                                                       vert_2);
        apf::ModelEntity* mdl_v = m->toModel(vert_1);
        int d_sp1 = m->getModelType(mdl_v);
        int c_sp1 = m->getModelTag(mdl_v);
        mdl_v = m->toModel(vert_2);
        int d_sp2 = m->getModelType(mdl_v);
        int c_sp2 = m->getModelTag(mdl_v);

        mdl_v = m->toModel(v_oth3);
        int d_sp3 = m->getModelType(mdl_v);
        int c_sp3 = m->getModelTag(mdl_v);
        mdl_v = m->toModel(v_oth4);
        int d_sp4 = m->getModelType(mdl_v);
        int c_sp4 = m->getModelTag(mdl_v);

        bool cell_ok = c_base->chk_conn_d_gmi(d_sp1, c_sp1, d_sp3, c_sp3) and 
                       c_base->chk_conn_d_gmi(d_sp2, c_sp2, d_sp4, c_sp4);

        bool cell_tri_ok = true;
        if(tri_split_type[tri_state] == 2) {
          mdl_v = m->toModel(d_e[e2]);
          int d_e_sp2 = m->getModelType(mdl_v);
          int c_e_sp2 = m->getModelTag(mdl_v);
          mdl_v = m->toModel(d_e3[e3]);
          int d_e_sp3 = m->getModelType(mdl_v);
          int c_e_sp3 = m->getModelTag(mdl_v);
          mdl_v = m->toModel(d_e4[e4]);
          int d_e_sp4 = m->getModelType(mdl_v);
          int c_e_sp4 = m->getModelTag(mdl_v);
          cell_tri_ok = !(d_e_sp2 == d_e_sp3 or
                          d_e_sp2 == d_e_sp4);
        }

        if(cell_tri_ok or cell_ok) {
          double q_curr = measureTetQuality(m, ref_0c, tet.at(i));
          std::vector<apf::Vector3> v(4, apf::Vector3(0,0,0));
          std::vector<apf::Vector3> v_new(4, apf::Vector3(0,0,0));
          for(int j = 0; j < 4; j++) {
            m->getPoint(d_v[j], 0, v.at(j));
            v_new.at(j) = v.at(j);
          }
          double q_init = calc_q(m, v);

          apf::Vector3 pos_e2(0,0,0);
          apf::Vector3 pos_e3(0,0,0);
          apf::Vector3 pos_e4(0,0,0);

          pos_e2 = vd_get_pos(m, d_e[e2]);
          pos_e3 = vd_get_pos(m, d_e3[e3]);
          pos_e4 = vd_get_pos(m, d_e4[e4]);
          // Check the quality of the new tets generated by triangle t1.
          // Shift the position corresponding to v1 or v2 to the center of e1.
          // Looping over the other v(not v1/v2), shift the position to the center 
          // of e3. Check quality.
          bool q_ok = true;
          v_new.at(v2) = pos_e2;
          int j = 0;
          while(j < 3 and q_ok) {
            int v_next = (v2 + 1 + j) % 4;
            if(v_next != v1) {
              v_new.at(v_next) = pos_e3;
              double q_curr = calc_q(m, v_new);
              q_ok = (q_ok and q_curr > q_frac*q_init);
              v_new.at(v_next) = v.at(v_next);
            }
            j = j + 1;
          }
          v_new.at(v2) = v.at(v2);
          v_new.at(v1) = pos_e2;
          j = 0;
          while(j < 3 and q_ok) {
            int v_next = (v1 + 1 + j) % 4;
            v_new.at(v_next) = pos_e4;
            if(v_next != v2) {
              double q_curr = calc_q(m, v_new);
              q_ok = (q_ok and q_curr > q_frac*q_init);
              v_new.at(v_next) = v.at(v_next);
            }
            j = j + 1;
          }
          if(q_ok) {

            ref_0c->split_map[d_e[e2]] = true;
            ref_0c->split_map[d_e3[e3]] = true;
            ref_0c->split_map[d_e4[e4]] = true;
          }
        }
      }
    }
  }

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = true;

  in->shouldCoarsen = false;
  in->shouldCoarsenLayer = false;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  if(sub_vtk)
    save_vtk_name("./output/adapt_prob_sp");
}

// Given the low quality tets, set the edge across the shortest edge 
// and the edges bounded by vertices belonging to lower dimensional strata than 
// the edge of the tets to be split.
// Only split.
void vd_sim::adapt_prob_edge_sp(std::vector<apf::MeshEntity*> &tet) {
  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  double adapt_ln = get_adapt_ln();
  double m_len =  min_cell.at(0)/2;

  Step_ns sf(m, adapt_ln);
  ma::Input* in = ma::configure(m, &sf);

//old solution and mix of this one
  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  repl_sz_field(in, m, MA_SIZE_TYPE::EDGE_SPLIT);

  std::vector<apf::MeshEntity*> verts(0);
  std::vector<apf::MeshEntity*> edges(0);
  vd_set_down(m, &tet, &verts, 3);
  vd_set_up(m, &verts, &edges);

  ModelEdgeSplit* ref_0c = (ModelEdgeSplit*) in->sizeField;

  //for(int i = 0; i < edges.size(); i++)
  //  ref_0c->split_map[edges.at(i)] = true;
  if(sub_vtk)
    vd_save_vtk_ent(m, &tet, "./output/adapt_prob_sp_ents");

  //vd_bipy* bipy_split = new vd_bipy(m, c_base, &f_calc);
  vd_bipy* bipy_split = new vd_bipy(m, c_base);
  std::vector<apf::MeshEntity*> tris(0);
  tris.reserve(tet.size());


  apf::Up up;
  apf::Downward d_s;
  apf::Downward d_e;
  apf::Downward d_v;
  apf::ModelEntity* mdl;
  int type_e;
  int type_v0;
  int type_v1;

  // Loop over tets marking triangles not belonging to a boundary as seen. 
  // If there is a single 2cell triangle, mark the tet. 
  // If a tet marked as such shares a single 3cell triangle with any other low 
  // quality tet, mark the edge across the edge connecting the shared triangle and
  // the 2cell triangle to be refined.
  // TODO it may be possible to have a lens-like set of pancake elements exists 
  // and if for some reason mesh-adapt doesn't fix those, also consider elements
  // with two neighbours.
  // unfurl the low quality tets in successive iterations. 
  // Given a tri on a tet, return the indices of the edges not bounding the tri.
  // Note:
  // The old algorithm where the longest edge not belonging to the 2cell triangle 
  // was split was causing excessive refinement when a cluster of low quality 
  // tetrahedra are tried to be fixed.

  int lookup_tri_e_other [4][3] = {{3,4,5}, {1,2,4}, {0,2,3}, {0,1,4}};
  std::map<apf::MeshEntity*, int> tri_id {};
  std::map<apf::MeshEntity*, int> tri_count {};
  std::map<apf::MeshEntity*, bool> tet_mark {};
  // If the minimum edge length is less than a fraction of the topological 
  // edgelength, don't refine.
  int sk_sz = 0;
  int tet_sz = tet.size();

  for(int i = 0; i < tet_sz - sk_sz; i++) {
    apf::MeshEntity* tet_curr = tet.at(i);
    m->getDownward(tet_curr, 1, d_e);
    for(int j = 0; j < 6; j++) {
      if(vd_meas_ent(m, d_e[j]) < len_topo/4) {
        tet.at(i) = tet.at(tet_sz - sk_sz - 1);
        tet.at(tet_sz - sk_sz - 1) = tet_curr;
        sk_sz = sk_sz + 1;
        i = i - 1;
        j = 6;
      }
    }
  }
  tet.resize(tet_sz - sk_sz);

  for(int i = 0; i < tet.size(); i++) {
    int d = 3;
    int t_id = -1;
    m->getDownward(tet.at(i), 2, d_s);
    int tri_state = 0;
    for(int j = 0; j < 4; j++) {
      mdl = m->toModel(d_s[j]);
      type_e = m->getModelType(mdl);
      if(type_e == 2) {
        tri_state = tri_state + tri_split_shift[j];
        tri_id[tet.at(i)] = j;
      }
      else
        tri_count[d_s[j]] = tri_count[d_s[j]] + 1;
    }
    if(tri_single_2c[tri_state] == 1)
      tet_mark[tet.at(i)] = true;
  }

  for(int i = 0; i < tet.size(); i++) {
    if(tet_mark[tet.at(i)]) {
      m->getDownward(tet.at(i), 2, d_s);
      m->getDownward(tet.at(i), 1, d_e);
      int count = 0;
      int t_id = -1;
      for(int j = 0; j < 4; j++) {
        if(tri_count[d_s[j]] == 2) {
          count = count + 1;
          t_id = j;
        }
      }
      if(count == 1) {
        int tri_state = tri_split_shift[tri_id[tet.at(i)]] + tri_split_shift[t_id];
        apf::MeshEntity* e1 = d_e[tri_split_map[tri_state]];
        ref_0c->split_map[e1] = true;
      }
      // TODO unless a lens-like set of pancake elements exists the above should
      // unfurl the low quality tets in successive iterations. 
      //else if(count == 2) {
      //}
    }
  }
/*
  for(int i = 0; i < tet.size(); i++) {
    int d = 3;
    m->getDownward(tet.at(i), 2, d_s);
    for(int j = 0; j < 4; j++) {
      mdl = m->toModel(d_s[j]);
      type_e = m->getModelType(mdl);
      if(type_e == 2)
        tri_state = tri_state + tri_split_shift[j];
    }

    int e1 = 0;
    if(tri_split_type[tri_state] == -1) {

    for(int j = 0; j < 4; j++) {
      if(tri_mark[d_s[j]]) {
        mdl = m->toModel(d_s[j]);
        type_e = m->getModelType(mdl);
        // If boundary tri, refining in this manner won't help. 
        if(type_e != 2) {
          m->getUp(d_s[j], up);
          assert(up.n == 2);

          m->getDownward(up.e[0], 1, d_e);
          m->getDownward(up.e[0], 2, d_v);
          int s1 = findIn(d_v, 4, d_s[j]);
          bool found1 = false;
          apf::MeshEntity* e1;
          double len_max = 0;
          int id_max = -1;
          for(int k = 0; k < 3; k++) {
            int id_curr = lookup_tri_e_other[s1][k];
            mdl = m->toModel(d_e[id_curr]);
            type_e = m->getModelType(mdl);
            if(type_e == 3) {
              found1 = true;
              double len_curr = vd_meas_ent(m, d_e[id_curr]);
              if(len_curr > len_max) {
                id_max = id_curr;
                len_max = len_curr;
              }
            }
          }
          if(found1) {
            e1 = d_e[id_max];
          }
          m->getDownward(up.e[1], 1, d_e);
          m->getDownward(up.e[1], 2, d_v);
          s1 = findIn(d_v, 4, d_s[j]);
          bool found2 = false;
          len_max = 0;
          id_max = -1;
          for(int k = 0; k < 3; k++) {
            int id_curr = lookup_tri_e_other[s1][k];
            mdl = m->toModel(d_e[id_curr]);
            type_e = m->getModelType(mdl);
            if(type_e == 3) {
              found2 = true;
              double len_curr = vd_meas_ent(m, d_e[id_curr]);
              if(len_curr > len_max) {
                id_max = id_curr;
                len_max = len_curr;
              }
            }
          }
          if(found1 and found2) {
vertex sharing also
            ref_0c->split_map[d_e[id_max]] = true;
            ref_0c->split_map[e1] = true;
          }
        }
      }
      else
        tri_mark[d_s[j]] = true;
    }
  }
*/

  for(int i = 0; i < tet.size(); i++) {
    int d = 3;
    m->getDownward(tet.at(i), 2, d_s);
    m->getDownward(tet.at(i), 1, d_e);
    int tri_state = 0;

    double area_max = 0;
    int j_max = -1;

    for(int j = 0; j < 4; j++) {
      mdl = m->toModel(d_s[j]);
      type_e = m->getModelType(mdl);
      if(type_e == 2)
        tri_state = tri_state + tri_split_shift[j];

      double area_curr = vd_area_out_n(m, d_s[j]).getLength();
      if(area_curr > area_max) {
        j_max = j;
        area_max = area_curr;
      }
    }

    int e1 = 0;
    if(tri_split_type[tri_state] == 0) {
      e1 = tri_split_map[tri_state];
      if(vd_chk_vol_valid(m, tet.at(i), d_e[e1]))
        ref_0c->split_map[d_e[e1]] = true;
    }
    else if(tri_split_type[tri_state] == 1) {
      e1 = tri_split_map[tri_state] - 6;
      if(vd_chk_vol_valid(m, tet.at(i), d_s[e1]))
        tris.push_back(d_s[e1]);
    }
    else if(tri_split_type[tri_state] == 2) {
      e1 = tri_split_map[tri_state] - 6;
      if(e1 == j_max and vd_chk_vol_valid(m, tet.at(i), d_s[e1]))
        tris.push_back(d_s[e1]);
    }

    double temp = vd_meas_ent(m, d_e[0]);

    for(int j = 1; j < 6; j++) {
      mdl = m->toModel(d_e[j]);
      type_e = m->getModelType(mdl);
      m->getDownward(d_e[j], 0, d_v);
      mdl = m->toModel(d_v[0]);
      type_v0 = m->getModelType(mdl);
      mdl = m->toModel(d_v[1]);
      type_v1 = m->getModelType(mdl);
      if(type_v0 < type_e and type_v1 < type_e and
         vd_chk_vol_valid(m, tet.at(i), d_e[j]) )
        ref_0c->split_map[d_e[j]] = true;

      double curr = vd_meas_ent(m, d_e[j]);
      if(curr > temp) {
        temp = curr;
        e1 = j;
      }
    }
    //ref_0c->split_map[d_e[lookup_t_e_x_e[e1]]] = true;
    //ref_0c->split_map[d_e[e1]] = false;
    //if(vd_chk_vol_valid(m, tet.at(i), d_e[e1]))
    //  ref_0c->split_map[d_e[e1]] = true;
  }
  std::vector<apf::MeshEntity*>::iterator it;
  it = std::unique(tris.begin(), tris.end());
  tris.resize(std::distance(tris.begin(),it));

  std::vector<apf::MeshEntity*> up_tet(0);
  for(int i = 0; i < tris.size(); i++) {
    bipy_split->load_tri(tris.at(i));
    bipy_split->split_bipy();
  }

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = true;

  in->shouldCoarsen = false;
  in->shouldCoarsenLayer = false;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);


  delete bipy_split;

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  if(sub_vtk)
    save_vtk_name("./output/adapt_prob_sp");
}


// Given the low quality tets, set the edge across the shortest edge 
// and the edges bounded by vertices belonging to lower dimensional strata than 
// the edge of the tets to be split.
// Only split.
void vd_sim::adapt_prob_edge_sp2(std::vector<apf::MeshEntity*> &tet) {
  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  Step_ns sf(m, get_adapt_ln());
  ma::Input* in = ma::configure(m, &sf);

//old solution and mix of this one
  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  repl_sz_field(in, m, MA_SIZE_TYPE::EDGE_SPLIT);

  std::vector<apf::MeshEntity*> verts(0);
  std::vector<apf::MeshEntity*> edges(0);
  vd_set_down(m, &tet, &verts, 3);
  vd_set_up(m, &verts, &edges);

  ModelEdgeSplit* ref_0c = (ModelEdgeSplit*) in->sizeField;

  //for(int i = 0; i < edges.size(); i++)
  //  ref_0c->split_map[edges.at(i)] = true;
  if(sub_vtk)
    vd_save_vtk_ent(m, &tet, "./output/adapt_prob_sp_ents");

  apf::Downward d_e;
  apf::Downward d_v;
  apf::ModelEntity* mdl;
  int type_e;
  int type_v0;
  int type_v1;
  for(int i = 0; i < tet.size(); i++) {
    int d = 3;
    m->getDownward(tet.at(i), 1, d_e);
    double temp = vd_meas_ent(m, d_e[0]);
    int e1 = 0;

    for(int j = 1; j < 6; j++) {
      mdl = m->toModel(d_e[j]);
      type_e = m->getModelType(mdl);
      m->getDownward(d_e[j], 0, d_v);
      mdl = m->toModel(d_v[0]);
      type_v0 = m->getModelType(mdl);
      mdl = m->toModel(d_v[1]);
      type_v1 = m->getModelType(mdl);
      if(type_v0 < type_e and type_v1 < type_e)
        ref_0c->split_map[d_e[j]] = true;

      double curr = vd_meas_ent(m, d_e[j]);
      if(curr > temp) {
        temp = curr;
        e1 = j;
      }
    }
    //ref_0c->split_map[d_e[lookup_t_e_x_e[e1]]] = true;
    //ref_0c->split_map[d_e[e1]] = false;
    ref_0c->split_map[d_e[e1]] = true;
  }

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = true;

  in->shouldCoarsen = false;
  in->shouldCoarsenLayer = false;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  if(sub_vtk)
    save_vtk_name("./output/adapt_prob_sp");
}

void vd_sim::adapt_prob_edge_sp3(std::vector<apf::MeshEntity*> &tet) {
  apf::Field* field_step = vd_att_vs_field(m, "adapt_step");
  Step_ns sf(m, get_adapt_ln());
  ma::Input* in = ma::configure(m, &sf);

  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  delete in->sizeField;
  ModelEdgeSplit* ref_0c = new ModelEdgeSplit(m);
  in->sizeField = ref_0c;

  for(int i = 0; i < tet.size(); i++) {
    int e1 = 0;
    int d = 3;
    apf::Downward d_e;
    m->getDownward(tet.at(i), 1, d_e);
    double temp = 0;
    for(int j = 0; j < 6; j++) {
      double curr = vd_meas_ent(m, d_e[j]);
      if(curr > temp) {
        temp = curr;
        e1 = j;
      }
    }
    ref_0c->split_map[d_e[e1]] = true;
  }

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = true;

  in->shouldCoarsen = false;
  in->shouldCoarsenLayer = false;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  if(sub_vtk)
    save_vtk_name("./output/adapt_prob_sp");
}
/*
void vd_sim::fix_low_q() {
  Step_ns sf(m, adapt_ln);
  ma::Input* in = ma::configure(m, &sf);

  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  delete in->sizeField;
  ModelEdgeRefinerDist* ref_0c = new ModelEdgeRefinerDist(m);
  ref_0c->set_target_th(coarse_th, split_th);
  in->sizeField = ref_0c;

  double ref_len = get_adapt_ln();
  double m_len =  min_cell.at(0)/2;

  //double good = (calc_good_q(m, adapt_ln, m_len)+0.2)/2;
  double good = (calc_good_q(m, std::min(adapt_ln, min_ins_rad), m_len)+0.2)/2;
  //std::pair<apf::MeshEntity*, double> valid = calc_valid_q(m, ref_0c);
  std::vector<apf::MeshEntity*> tets(0);
  double valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  //double good = 0.2;
  if(valid < q_th) {
    adapt_prob_edge(tets);
  }
  //int sz_temp = tets.size();
  //valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  //if(tets.size() == sz_temp) {
  //  adapt_prob_edge_low(tets);
  //  valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  //  if(valid < q_th) {
  //    adapt_prob_edge_sp(tets);
  //  }
  //}
  else if(valid < q_th) {
    adapt_prob_edge_sp(tets);
  }

}
*/
// For each tet, check if bounded by an exterior triangle, if not check if 
// bounded by a 2-stratum triangle. If found, return the vertex 
void vd_sim::get_low_q_verts(std::vector<apf::MeshEntity*>& tets,
                             std::vector<apf::MeshEntity*>& verts) {
}

void vd_sim::fix_low_q_split_tet() {

  double t0 = PCU_Time();
  //f_calc.del_tri_fields(m);
  if(ad_type == ADAPT_TYPE::ADAPT_BOUND)
    adapt_mark_bound_min();
  else if(ad_type == ADAPT_TYPE::ADAPT_CURVE)
    adapt_mark_curv();
  else
    adapt_mark_0c_min();

  double adapt_ln = get_adapt_ln();
  Step_ns sf(m, adapt_ln);
  ma::Input* in = ma::configure(m, &sf);

  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  delete in->sizeField;
  ModelEdgeRefinerDist* ref_0c = new ModelEdgeRefinerDist(m);
  ref_0c->set_target_th(coarse_th, split_th);
  in->sizeField = ref_0c;

  double m_len =  min_cell.at(0)/2;

  //std::pair<apf::MeshEntity*, double> valid = calc_valid_q(m, ref_0c);       
  std::vector<apf::MeshEntity*> tets(0);
  std::vector<apf::MeshEntity*> verts(0);
  double valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  double good = get_good_q();
  good = std::max(valid*1.1, good);

  std::cout << "Minimum quality is " << valid
            << ", good quality is " << good 
            << std::endl;
  if(sub_vtk)
    vd_save_vtk_ent(m, &tets, "./output/adapt_prob_sp_ents");

  get_low_q_verts(tets, verts);
  for(int i = 0; i < verts.size(); i++) {
    vd_split_tet_vert(m, c_base, verts.at(i), NULL);
  }

  if(ad_type == ADAPT_TYPE::ADAPT_BOUND)
    adapt_mark_bound_min();
  else if(ad_type == ADAPT_TYPE::ADAPT_CURVE)
    adapt_mark_curv();
  else
    adapt_mark_0c_min();

  in = ma::configure(m, &sf);
  delete in->sizeField;
  ref_0c = new ModelEdgeRefinerDist(m);
  ref_0c->set_target_th(coarse_th, split_th);
  in->sizeField = ref_0c;

  //in->shouldRunPreZoltan = true;
  //in->shouldRunMidParma = true;
  //in->shouldRunPostParma = true;
  //in->shouldRefineLayer = true;

  //in->shouldCoarsenLayer = true;
  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = false;
  in->shouldCoarsenLayer = false;
  in->shouldCoarsen = false;

  in->shouldFixShape = ad_fix_shape;
  in->goodQuality = good;
  in->validQuality = valid;
  //in->validQuality = 0;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  ma::adapt(in);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

}

// TODO adaptation function should be unified across all adaptation related 
// functions in vd_sim to prevent emergence of different metric criteria.
// In addition, the structure of the code should reflect that for easier 
// maintenance.
void vd_sim::fix_low_q() {

  double t0 = PCU_Time();
  //f_calc.del_tri_fields(m);
  if(ad_type == ADAPT_TYPE::ADAPT_BOUND)
    adapt_mark_bound_min();
  else if(ad_type == ADAPT_TYPE::ADAPT_CURVE)
    adapt_mark_curv();
  else
    adapt_mark_0c_min();

  double adapt_ln = get_adapt_ln();
  Step_ns sf(m, adapt_ln);
  ma::Input* in = ma::configure(m, &sf);

  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  delete in->sizeField;
  ModelEdgeRefinerDist* ref_0c = new ModelEdgeRefinerDist(m);
  ref_0c->set_target_th(coarse_th, split_th);
  in->sizeField = ref_0c;

  double m_len =  min_cell.at(0)/2;

  //std::pair<apf::MeshEntity*, double> valid = calc_valid_q(m, ref_0c);
  std::vector<apf::MeshEntity*> tets(0);
  double valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  double good = get_good_q();
  good = std::max(valid*1.1, good);
  delete in;
  if(valid < q_th) {
    adapt_prob_edge(tets);
  }
  std::cout << "Minimum quality is " << valid
            << ", good quality is " << good 
            << std::endl;

  if(ad_type == ADAPT_TYPE::ADAPT_BOUND)
    adapt_mark_bound_min();
  else if(ad_type == ADAPT_TYPE::ADAPT_CURVE)
    adapt_mark_curv();
  else
    adapt_mark_0c_min();

  in = ma::configure(m, &sf);
  delete in->sizeField;
  ref_0c = new ModelEdgeRefinerDist(m);
  ref_0c->set_target_th(coarse_th, split_th);
  in->sizeField = ref_0c;

  //in->shouldRunPreZoltan = true;
  //in->shouldRunMidParma = true;
  //in->shouldRunPostParma = true;
  //in->shouldRefineLayer = true;

  //in->shouldCoarsenLayer = true;
  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = false;
  in->shouldCoarsenLayer = false;
  in->shouldCoarsen = false;

  in->shouldFixShape = ad_fix_shape;
  in->goodQuality = good;
  in->validQuality = valid;
  //in->validQuality = 0;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  ma::adapt(in);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  //int sz_temp = tets.size();
  //valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  //if(tets.size() == sz_temp) {
  //  adapt_prob_edge_low(tets);
  //  valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  //  if(valid < q_th) {
  //    adapt_prob_edge_sp(tets);
  //  }
  //}

  // adapt_sp_bound
  if(ad_type == ADAPT_TYPE::ADAPT_BOUND)
    adapt_mark_bound_min();
  else if(ad_type == ADAPT_TYPE::ADAPT_CURVE)
    adapt_mark_curv();
  else
    adapt_mark_0c_min();

  in = ma::configure(m, &sf);
  delete in->sizeField;
  ref_0c = new ModelEdgeRefinerDist(m);
  ref_0c->set_target_th(coarse_th, split_th);
  in->sizeField = ref_0c;

  valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  //good = std::max(valid*1.1, good);
  delete in;
  if(valid < q_th) {
    // Refine all edges upward to the bounding vertices of the problematic tets
    // Run the general adaptation.
    //adapt_prob_edge_sp(tets);
    if(tets.size() < 10)
      adapt_prob_edge_sp_col(tets);
    else
      adapt_prob_edge_sp_bound(tets);
  }
  std::cout << "Minimum quality is " << valid
            << ", good quality is " << good 
            << std::endl;

  if(ad_type == ADAPT_TYPE::ADAPT_BOUND)
    adapt_mark_bound_min();
  else if(ad_type == ADAPT_TYPE::ADAPT_CURVE)
    adapt_mark_curv();
  else
    adapt_mark_0c_min();

  in = ma::configure(m, &sf);
  delete in->sizeField;
  ref_0c = new ModelEdgeRefinerDist(m);
  ref_0c->set_target_th(coarse_th, split_th);
  in->sizeField = ref_0c;

  //in->shouldRunPreZoltan = true;
  //in->shouldRunMidParma = true;
  //in->shouldRunPostParma = true;
  //in->shouldRefineLayer = true;

  //in->shouldCoarsenLayer = true;
  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = false;
  in->shouldCoarsenLayer = false;
  in->shouldCoarsen = false;

  in->shouldFixShape = ad_fix_shape;
  in->goodQuality = good;
  in->validQuality = valid;
  //in->validQuality = 0;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 1;

  ma::adapt(in);

  a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  // TODO manual operation:
  bool bound_manual = false;
  if(bound_manual ) {
     //or bound_manual_count > bound_manual_count_limit) {

    if(ad_type == ADAPT_TYPE::ADAPT_BOUND)
      adapt_mark_bound_min();
    else if(ad_type == ADAPT_TYPE::ADAPT_CURVE)
      adapt_mark_curv();
    else
      adapt_mark_0c_min();

    in = ma::configure(m, &sf);
    delete in->sizeField;
    ref_0c = new ModelEdgeRefinerDist(m);
    ref_0c->set_target_th(coarse_th, split_th);
    in->sizeField = ref_0c;

    valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
    //good = std::max(valid*1.1, good);
    delete in;
    if(valid < q_th) {
      // Refine all edges upward to the bounding vertices of the problematic tets
      // Run the general adaptation.
      //adapt_prob_edge_sp(tets);
      adapt_prob_edge_col_bound_manual(tets);

      std::cout << "Minimum quality is " << valid
                << ", good quality is " << good 
                << std::endl;

      if(ad_type == ADAPT_TYPE::ADAPT_BOUND)
        adapt_mark_bound_min();
      else if(ad_type == ADAPT_TYPE::ADAPT_CURVE)
        adapt_mark_curv();
      else
        adapt_mark_0c_min();

      in = ma::configure(m, &sf);
      delete in->sizeField;
      ref_0c = new ModelEdgeRefinerDist(m);
      ref_0c->set_target_th(coarse_th, split_th);
      in->sizeField = ref_0c;

      //in->shouldRunPreZoltan = true;
      //in->shouldRunMidParma = true;
      //in->shouldRunPostParma = true;
      //in->shouldRefineLayer = true;

      //in->shouldCoarsenLayer = true;
      in->shouldRunPreZoltan = false;
      in->shouldRunMidParma = false;
      in->shouldRunPostParma = false;
      in->shouldRefineLayer = false;
      in->shouldCoarsenLayer = false;
      in->shouldCoarsen = false;

      in->shouldFixShape = ad_fix_shape;
      in->goodQuality = good;
      in->validQuality = valid;
      //in->validQuality = 0;
      in->maximumEdgeRatio = 12;
      in->maximumIterations = 1;

      ma::adapt(in);

      a_f = m->findField("adapt_step");
      if(a_f)
        apf::destroyField(a_f);
    }
    if(ad_type == ADAPT_TYPE::ADAPT_BOUND)
      adapt_mark_bound_min();
    else if(ad_type == ADAPT_TYPE::ADAPT_CURVE)
      adapt_mark_curv();
    else
      adapt_mark_0c_min();

    in = ma::configure(m, &sf);
    delete in->sizeField;
    ref_0c = new ModelEdgeRefinerDist(m);
    ref_0c->set_target_th(coarse_th, split_th);
    in->sizeField = ref_0c;

    valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
    //good = std::max(valid*1.1, good);
    delete in;
    if(valid < q_th) {
      // Refine all edges upward to the bounding vertices of the problematic tets
      // Run the general adaptation.
      //adapt_prob_edge_sp(tets);
      adapt_prob_edge_sp_bound_manual(tets);

      std::cout << "Minimum quality is " << valid
                << ", good quality is " << good 
                << std::endl;

      if(ad_type == ADAPT_TYPE::ADAPT_BOUND)
        adapt_mark_bound_min();
      else if(ad_type == ADAPT_TYPE::ADAPT_CURVE)
        adapt_mark_curv();
      else
        adapt_mark_0c_min();

      in = ma::configure(m, &sf);
      delete in->sizeField;
      ref_0c = new ModelEdgeRefinerDist(m);
      ref_0c->set_target_th(coarse_th, split_th);
      in->sizeField = ref_0c;

      //in->shouldRunPreZoltan = true;
      //in->shouldRunMidParma = true;
      //in->shouldRunPostParma = true;
      //in->shouldRefineLayer = true;

      //in->shouldCoarsenLayer = true;
      in->shouldRunPreZoltan = false;
      in->shouldRunMidParma = false;
      in->shouldRunPostParma = false;
      in->shouldRefineLayer = false;
      in->shouldCoarsenLayer = false;
      in->shouldCoarsen = false;

      in->shouldFixShape = ad_fix_shape;
      in->goodQuality = good;
      in->validQuality = valid;
      //in->validQuality = 0;
      in->maximumEdgeRatio = 12;
      in->maximumIterations = 1;

      ma::adapt(in);

      a_f = m->findField("adapt_step");
      if(a_f)
        apf::destroyField(a_f);
    }
    bound_manual_count = 0;
  }
  else
    bound_manual_count = bound_manual_count + 1;
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);
}
/*
// TODO adaptation function should be unified across all adaptation related 
// functions in vd_sim to prevent emergence of different metric criteria.
// In addition, the structure of the code should reflect that for easier 
// maintenance.
void vd_sim::fix_low_q() {
  //f_calc.del_tri_fields(m);
  
  adapt_mark_0c_min();
  double adapt_ln = get_adapt_ln();
  Step_ns sf(m, adapt_ln);
  ma::Input* in = ma::configure(m, &sf);

  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  delete in->sizeField;
  ModelEdgeRefinerDist* ref_0c = new ModelEdgeRefinerDist(m);
  ref_0c->set_target_th(coarse_th, split_th);
  in->sizeField = ref_0c;

  double m_len =  min_cell.at(0)/2;

  //double good = (calc_good_q(m, adapt_ln, m_len)+0.2)/2;
  //std::pair<apf::MeshEntity*, double> valid = calc_valid_q(m, ref_0c);
  std::vector<apf::MeshEntity*> tets(0);
  double valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  double good = std::min(0.1, calc_good_q(m, adapt_ln, m_len)*50);
  //double good = std::min(0.1, calc_good_q(m, std::min(adapt_ln, min_ins_rad), m_len)*50);
  good = std::max(valid*1.1, good);
  //double good = 0.1;

  //double good = 0.2;
  delete in;
  if(valid < q_th) {
    adapt_prob_edge(tets);
  }
  //int sz_temp = tets.size();
  //valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  //if(tets.size() == sz_temp) {
  //  adapt_prob_edge_low(tets);
  //  valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  //  if(valid < q_th) {
  //    adapt_prob_edge_sp(tets);
  //  }
  //}

  adapt_mark_0c_min();

  in = ma::configure(m, &sf);
  delete in->sizeField;
  ref_0c = new ModelEdgeRefinerDist(m);
  ref_0c->set_target_th(coarse_th, split_th);
  in->sizeField = ref_0c;

  valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  //good = std::max(valid*1.1, good);
  delete in;
  if(valid < q_th) {
    // Refine all edges upward to the bounding vertices of the problematic tets
    // Run the general adaptation.
    adapt_prob_edge_sp(tets);
  }
  std::cout << "Minimum quality is " << valid
            << ", good quality is " << good 
            << std::endl;

  adapt_mark_0c_min();
  in = ma::configure(m, &sf);
  delete in->sizeField;
  ref_0c = new ModelEdgeRefinerDist(m);
  ref_0c->set_target_th(coarse_th, split_th);
  in->sizeField = ref_0c;

  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;

  in->shouldCoarsenLayer = true;

  in->shouldFixShape = ad_fix_shape;
  in->goodQuality = good;
  in->validQuality = valid;
  //in->validQuality = 0;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 3;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

}
*/
/*
// TODO adaptation function should be unified across all adaptation related 
// functions in vd_sim to prevent emergence of different metric criteria.
// In addition, the structure of the code should reflect that for easier 
// maintenance.
void vd_sim::fix_low_q_bad() {
  //f_calc.del_tri_fields(m);
  ma::Input* in;
  
  double adapt_ln = get_adapt_ln();
  if(ad_type == ADAPT_TYPE::ADAPT_BOUND) {
    adapt_ln = std::min(adapt_ln, bound_len*4);
    Linear sf(m, adapt_ln);
    in = ma::configure(m, &sf);
  }
  else {
    adapt_mark_0c_min();
    Step_ns sf(m, adapt_ln);
    in = ma::configure(m, &sf);
  }

  // Use a different sizefield to split each 1stratum into two edges 
  // at least.
  // Doing only this causes time step issues. The quality of elements 
  // at the boundaries should be high. 
  delete in->sizeField;

  double m_len =  min_cell.at(0)/2;

  //double good = (calc_good_q(m, adapt_ln, m_len)+0.2)/2;
  //std::pair<apf::MeshEntity*, double> valid = calc_valid_q(m, ref_0c);
  std::vector<apf::MeshEntity*> tets(0);
  double valid = 1;
  double good = 1;
  
  if(ad_type == ADAPT_TYPE::ADAPT_BOUND) {
    ModelEdgeRefinerVarying* ref_0c_vary = new ModelEdgeRefinerVarying(m);
    ref_0c_vary->set_target_th(coarse_th, split_th);
    ref_0c_vary->set_all_cell(adapt_ln);
    if(cells_ref.size() != 0) {
      ref_0c_vary->set_bound_cell(bound_len, &cells_ref);
    }
    else
      ref_0c_vary->set_bound_cell(bound_len, NULL);
    in->sizeField = ref_0c_vary;
    valid = get_low_q(m, ref_0c_vary, tets, q_th)*0.99;
  }
  else {
    adapt_mark_0c_min();
    ModelEdgeRefinerDist* ref_0c = new ModelEdgeRefinerDist(m);
    ref_0c->set_target_th(coarse_th, split_th);
    in->sizeField = ref_0c;
    valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  }
  good = std::min(0.1, calc_good_q(m, adapt_ln, m_len)*50);
  //double good = std::min(0.1, calc_good_q(m, std::min(adapt_ln, min_ins_rad), m_len)*50);
  good = std::max(valid*1.1, good);
  //double good = 0.1;

  //double good = 0.2;
  delete in;
  if(valid < q_th) {
    adapt_prob_edge(tets);
  }
  //int sz_temp = tets.size();
  //valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  //if(tets.size() == sz_temp) {
  //  adapt_prob_edge_low(tets);
  //  valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  //  if(valid < q_th) {
  //    adapt_prob_edge_sp(tets);
  //  }
  //}

  if(ad_type == ADAPT_TYPE::ADAPT_BOUND) {
    Linear sf(m, adapt_ln);
    in = ma::configure(m, &sf);
  }
  else {
    Step_ns sf(m, adapt_ln);
    in = ma::configure(m, &sf);
  }

  delete in->sizeField;

  if(ad_type == ADAPT_TYPE::ADAPT_BOUND) {
    ModelEdgeRefinerVarying* ref_0c_vary = new ModelEdgeRefinerVarying(m);
    ref_0c_vary->set_target_th(coarse_th, split_th);
    ref_0c_vary->set_all_cell(adapt_ln);
    if(cells_ref.size() != 0) {
      ref_0c_vary->set_bound_cell(bound_len, &cells_ref);
    }
    else
      ref_0c_vary->set_bound_cell(bound_len, NULL);
    in->sizeField = ref_0c_vary;
    valid = get_low_q(m, ref_0c_vary, tets, q_th)*0.99;
  }
  else {
    adapt_mark_0c_min();
    ModelEdgeRefinerDist* ref_0c = new ModelEdgeRefinerDist(m);
    ref_0c->set_target_th(coarse_th, split_th);
    in->sizeField = ref_0c;
    valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
  }

  //good = std::max(valid*1.1, good);
  delete in;
  if(valid < q_th) {
    // Refine all edges upward to the bounding vertices of the problematic tets
    // Run the general adaptation.
    adapt_prob_edge_sp(tets);
  }
  std::cout << "Minimum quality is " << valid
            << ", good quality is " << good 
            << std::endl;

  if(ad_type == ADAPT_TYPE::ADAPT_BOUND) {
    Linear sf(m, adapt_ln);
    in = ma::configure(m, &sf);
  }
  else {
    Step_ns sf(m, adapt_ln);
    in = ma::configure(m, &sf);
  }

  delete in->sizeField;
  
  if(ad_type == ADAPT_TYPE::ADAPT_BOUND) {
    ModelEdgeRefinerVarying* ref_0c_vary = new ModelEdgeRefinerVarying(m);
    ref_0c_vary->set_target_th(coarse_th, split_th);
    ref_0c_vary->set_all_cell(adapt_ln);
    if(cells_ref.size() != 0) {
      ref_0c_vary->set_bound_cell(bound_len, &cells_ref);
    }
    else
      ref_0c_vary->set_bound_cell(bound_len, NULL);
    in->sizeField = ref_0c_vary;
  }
  else {
    adapt_mark_0c_min();
    ModelEdgeRefinerDist* ref_0c = new ModelEdgeRefinerDist(m);
    ref_0c->set_target_th(coarse_th, split_th);
    in->sizeField = ref_0c;
  }

  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;

  in->shouldCoarsenLayer = true;

  in->shouldFixShape = ad_fix_shape;
  in->goodQuality = good;
  in->validQuality = valid;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 3;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  Linear sf(m, adapt_ln);
  in = ma::configure(m, &sf);

  valid = get_low_q(m, in->sizeField, tets, q_th)*0.99;
  vd_save_vtk_ent(m, &tet, "./output/adapt_prob_sp_ents_after");
  delete in;
}
*/

double vd_sim::get_good_q() {
  //double adapt_ln = get_adapt_ln();
  //double m_len =  min_cell.at(0)/2;
  //return std::min(0.1, calc_good_q(m, adapt_ln, m_len)*5e4);
  //return (calc_good_q(m, adapt_ln, m_len)+0.2)/2;
  return 0.1;
}

void vd_sim::adapt(bool avg_flag) {
  vd_tim.set_time((int)T_COST_TYPE::T_COST_FLAG_ADAPT);

  //f_calc.del_tri_fields(m);
  f_calc.vd_del_fields(m);

  //Ref_Bound sf(m, adapt_len);
  e_list.refresh();
  if(avg_flag and len_f and !len_col_f) {
    clear_rad_map();
    avg_cell = upd_cell_rad_sim();
    get_length_scale();
  }
  if(ad_type == ADAPT_TYPE::ADAPT_1CELL or 
     ad_type == ADAPT_TYPE::ADAPT_3CELL) {

    double adapt_ln = get_adapt_ln();
    Linear sf(m, adapt_ln);

    ma::Input* in = ma::configure(m, &sf);
    in->shouldRunPreZoltan = true;
    in->shouldRunMidParma = true;
    in->shouldRunPostParma = true;
    in->shouldRefineLayer = true;

    in->shouldFixShape = ad_fix_shape;

    std::vector<apf::MeshEntity*> tets(0);
    double m_len =  min_cell.at(0)/2;
    double valid = get_low_q(m, in->sizeField, tets, q_th);
    double good = get_good_q();
    good = std::max(valid*1.1, good);
    //double good = 0.1;

    in->goodQuality = good;
    //in->goodQuality = 0.10;
    //in->validQuality = 0.2;

    double t0 = PCU_Time();
    ma::adapt(in);
    double t1 = PCU_Time();
    dt_adapt = dt_adapt + (t1 - t0);
  }
  else if(ad_type == ADAPT_TYPE::ADAPT_BOUND) {
    double adapt_ln = get_adapt_ln();
    double m_len =  min_cell.at(0)/2;
    //double adapt_ln_mod = std::min(adapt_ln, bound_len*4);
    ModelEdgeRefinerVarying* ref_0c = new ModelEdgeRefinerVarying(m);
    ref_0c->set_target_th(coarse_th, split_th);

    std::vector<apf::MeshEntity*> tets(0);
    m_len =  min_cell.at(0)/2;
    double valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
    delete ref_0c;
    if(valid < q_th) {
      for(int i = 0; i < 1; i++) {
        fix_low_q();
        //fix_low_q_split_tet();
      }
    }

    adapt_mark_bound_min();
    Step_ns sf(m, bound_len);
    //Linear sf(m, adapt_ln);

    ma::Input* in = ma::configure(m, &sf);

    ModelEdgeRefinerDist* ref_dist = new ModelEdgeRefinerDist(m);
    ref_dist->set_target_th(coarse_th, split_th);

    delete in->sizeField;
    in->sizeField = ref_dist;
/*
    ref_0c = (ModelEdgeRefinerVarying*) in->sizeField;
    delete in->sizeField;
    ref_0c = new ModelEdgeRefinerVarying(m);
    ref_0c->set_all_cell(adapt_ln_mod);
    if(cells_ref.size() != 0) {
      ref_0c->set_bound_cell(bound_len, &cells_ref);
    }
    else
      ref_0c->set_bound_cell(bound_len, NULL);

    in->sizeField = ref_0c;
*/

    tets.clear();
    valid = get_low_q(m, in->sizeField, tets, q_th);
    double good = get_good_q();
    good = std::max(valid*1.1, good);

    std::cout << "Minimum quality is " << valid
              << ", good quality is " << good 
              << std::endl;

    in->shouldRunPreZoltan = true;
    in->shouldRunMidParma = true;
    in->shouldRunPostParma = true;
    in->shouldRefineLayer = true;

    in->shouldCoarsenLayer = true;

    in->shouldFixShape = ad_fix_shape;
    in->goodQuality = good;
    in->validQuality = valid;
    in->maximumEdgeRatio = 12;
    in->maximumIterations = 1;

    double t0 = PCU_Time();
    ma::adapt(in);
    double t1 = PCU_Time();
    dt_adapt = dt_adapt + (t1 - t0);

    apf::Field* a_f = m->findField("adapt_step");
    if(a_f)
      apf::destroyField(a_f);

    if(sub_vtk)
      save_vtk_name("./output/adapt_iter");
  }

  else if(ad_type == ADAPT_TYPE::ADAPT_CURVE) {
    double adapt_ln = get_adapt_ln();
    double m_len =  min_cell.at(0)/2;
    //double adapt_ln_mod = std::min(adapt_ln, bound_len*4);
    ModelEdgeRefinerVarying* ref_0c = new ModelEdgeRefinerVarying(m);
    ref_0c->set_target_th(coarse_th, split_th);

    std::vector<apf::MeshEntity*> tets(0);
    m_len =  min_cell.at(0)/2;
    double valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
    delete ref_0c;
    if(valid < q_th) {
      for(int i = 0; i < 1; i++) {
        fix_low_q();
        //fix_low_q_split_tet();
      }
    }

    adapt_mark_curv();
    Step_ns sf(m, bound_len);

    ma::Input* in = ma::configure(m, &sf);

    ModelEdgeRefinerDist* ref_dist = new ModelEdgeRefinerDist(m);
    ref_dist->set_target_th(coarse_th, split_th);

    delete in->sizeField;
    in->sizeField = ref_dist;

    tets.clear();
    valid = get_low_q(m, in->sizeField, tets, q_th);
    double good = get_good_q();
    good = std::max(valid*1.1, good);

    std::cout << "Minimum quality is " << valid
              << ", good quality is " << good 
              << std::endl;

    in->shouldRunPreZoltan = true;
    in->shouldRunMidParma = true;
    in->shouldRunPostParma = true;
    in->shouldRefineLayer = true;

    in->shouldCoarsenLayer = true;

    in->shouldFixShape = ad_fix_shape;
    in->goodQuality = good;
    in->validQuality = valid;
    in->maximumEdgeRatio = 12;
    in->maximumIterations = 1;

    double t0 = PCU_Time();
    ma::adapt(in);
    double t1 = PCU_Time();
    dt_adapt = dt_adapt + (t1 - t0);

    apf::Field* a_f = m->findField("adapt_step");
    if(a_f)
      apf::destroyField(a_f);

    if(sub_vtk)
      save_vtk_name("./output/adapt_iter");
  }
  else {
    double adapt_ln = get_adapt_ln();
    double m_len =  min_cell.at(0)/2;
    //double good = (calc_good_q(m, adapt_ln, m_len)+0.2)/2;

    if(ad_type == ADAPT_TYPE::ADAPT_STEP) {
      adapt_mark_step();
      Step_ns sf(m, adapt_ln);
      ma::Input* in = ma::configure(m, &sf);

      std::vector<apf::MeshEntity*> tets(0);
      double valid = get_low_q(m, in->sizeField, tets, q_th);
      double good = get_good_q();

      good = std::max(valid*1.1, good);
      //double good = 0.1;

      in->shouldRunPreZoltan = true;
      in->shouldRunMidParma = true;
      in->shouldRunPostParma = true;
      in->shouldRefineLayer = true;

      in->shouldFixShape = ad_fix_shape;
      in->goodQuality = good;
      //in->goodQuality = 0.10;

      double t0 = PCU_Time();
      ma::adapt(in);
      double t1 = PCU_Time();
      dt_adapt = dt_adapt + (t1 - t0);
      apf::Field* a_f = m->findField("adapt_step");
      if(a_f)
        apf::destroyField(a_f);
    }

    else {
      assert(ad_type == ADAPT_TYPE::ADAPT_STEP_1CELL);
      //adapt_mark_1c_min();
      std::map<apf::MeshEntity*, double> v_dist{};

      adapt_mark_0c_min();
      ModelEdgeRefinerDist* ref_0c = new ModelEdgeRefinerDist(m);
      ref_0c->set_target_th(coarse_th, split_th);

      std::vector<apf::MeshEntity*> tets(0);
      m_len =  min_cell.at(0)/2;
      double valid = get_low_q(m, ref_0c, tets, q_th)*0.99;
      delete ref_0c;
      if(valid < q_th) {
        for(int i = 0; i < 2; i++) {
          fix_low_q();
	        //fix_low_q_split_tet();
        }
      }

      adapt_mark_0c_min();
      Step_ns sf(m, adapt_ln);
      ma::Input* in = ma::configure(m, &sf);
      delete in->sizeField;
      ref_0c = new ModelEdgeRefinerDist(m);
      ref_0c->set_target_th(coarse_th, split_th);

      in->sizeField = ref_0c;

      m_len =  min_cell.at(0)/2;
      valid = get_low_q(m, in->sizeField, tets, q_th);
      //double good = (calc_good_q(m, adapt_ln, m_len)+0.2)/2;
      double good = get_good_q();

      good = std::max(valid*1.1, good);
      //double good = 0.1;

      in->shouldRunPreZoltan = true;
      in->shouldRunMidParma = true;
      in->shouldRunPostParma = true;
      in->shouldRefineLayer = true;

      in->shouldCoarsenLayer = true;

      in->shouldFixShape = ad_fix_shape;
      in->goodQuality = good;
      in->validQuality = valid;
      //in->validQuality = 0;
      in->maximumEdgeRatio = 12;
      in->maximumIterations = 3;

      double t0 = PCU_Time();
      ma::adapt(in);
      double t1 = PCU_Time();
      dt_adapt = dt_adapt + (t1 - t0);
//TODO fixElementShapes actually removes split 1c edges, to improve quality
// if in->shouldFixShape == true
// ma::adapt is self contained. Unless a new adaptation function is written
// by making the ma headers available in the library in->maximumEdgeRatio and
// in->validQuality and in->goodQuality should be adjusted. But this will
// affect other elements, as well. Maximum edge ratio could be based on the 
// median edge length and the shortest 1 stratum length.
      apf::Field* a_f = m->findField("adapt_step");
      if(a_f)
        apf::destroyField(a_f);

      m->acceptChanges();
      m->verify();
      e_list.refresh();

      if(sub_vtk)
        save_vtk_name("./output/adapt_iter");

    }
  }
  m->acceptChanges();
  m->verify();

  e_list.refresh();
  f_calc.vd_att_fields(m);
  //save_vtk();
  vd_tim.add_time((int)T_COST_TYPE::T_COST_FLAG_ADAPT);
}

// TODO This should focus only near the collapsing stratum. Right now it adapts
// the whole mesh, which might be inefficient.
void vd_sim::adapt_col(int dim_col, int tag_col) {
  f_calc.del_tri_fields(m);

  assert(ad_type == ADAPT_TYPE::ADAPT_STEP_1CELL);
  double adapt_ln = get_adapt_ln();
  double m_len =  min_cell.at(0)/2;

  //adapt_mark_1c_min();
  adapt_mark_0c_min();
  adapt_mark_0c_col(dim_col, tag_col);
  if(sub_vtk)
    save_vtk_name("./output/before_adapt");

  Step_ns sf(m, adapt_ln);
  ma::Input* in = ma::configure(m, &sf);

  delete in->sizeField;
  ModelEdgeRefinerDist* ref_0c = new ModelEdgeRefinerDist(m);
  ref_0c->set_target_th(coarse_th, split_th);

  ref_0c->coarse_map.at(dim_col-1)[tag_col] = true;
  if(dim_col == 3) {
    ent_conn* edown = new ent_conn();
    c_base->get_conn_gmi(dim_col, tag_col, edown);
    for(int i = 0; i < edown->conn.size(); i++) {
      ref_0c->coarse_map.at(dim_col-2)[edown->conn.at(i)] = true;
    }
    c_base->get_conn_dim_gmi(1, dim_col, tag_col, edown);
    for(int i = 0; i < edown->conn.size(); i++) {
      ref_0c->coarse_map.at(dim_col-3)[edown->conn.at(i)] = true;
    }
    delete edown;
  }
  else if(dim_col == 2) {
    ent_conn* edown = new ent_conn();
    c_base->get_conn_gmi(dim_col, tag_col, edown);
    for(int i = 0; i < edown->conn.size(); i++) {
      ref_0c->coarse_map.at(dim_col-2)[edown->conn.at(i)] = true;
    }
    delete edown;
  }

  in->sizeField = ref_0c;

  double valid = 1;
  apf::MeshEntity* elem;
  apf::MeshIterator* it = m->begin(3);
  while(elem = m->iterate(it)) {
    double q_temp = measureTetQuality(m, ref_0c, elem);
    if(q_temp < valid)
      valid = q_temp;
  }
  m->end(it);
  std::cout << "Minimum quality is " << valid << std::endl;

  double good = get_good_q();

  good = std::max(valid*1.1, good);
  //double good = 0.1;

  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;

  in->shouldCoarsenLayer = true;
  in->shouldCoarsen = true;

  in->shouldFixShape = ad_fix_shape;
  in->goodQuality = good;
  //in->goodQuality = 0.03;
  in->validQuality = valid;
  //in->validQuality = 0;
  in->maximumEdgeRatio = 12;
  in->maximumIterations = 3;

  double t0 = PCU_Time();
  ma::adapt(in);
  double t1 = PCU_Time();
  dt_adapt = dt_adapt + (t1 - t0);

  f_calc.vd_att_fields(m);

  apf::Field* a_f = m->findField("adapt_step");
  if(a_f)
    apf::destroyField(a_f);

  m->acceptChanges();
  m->verify();

  e_list.refresh();

  //save_vtk();

}

void vd_sim::clean_mesh() {

  if (load_flag) {

    std::cout << "VD_SIM: Cleaning up... ";
    std::cout << "Deleting c_base... ";
    //std::cout << "Deleting g_lens... ";
    //delete g_lens;

    //std::cout << "Clearing c_ins... ";
    //c_ins.clear();

    //c_base = NULL;
    f_calc.off_cb();

    //g_lens = NULL;
    std::cout << "VD_SIM: Cleaned up." << std::endl;

    load_flag = false;
    mesh_flag = false;

    vd_rem_tag(m);
    m->destroyNative();
    apf::destroyMesh(m);
  }

}

void vd_sim::clean_up() {

  if (load_flag) {

    std::cout << "VD_SIM: Cleaning up... ";
    std::cout << "Deleting c_base... ";
    if(!first_time)
      delete c_base;
    //c_ins.clear();

    c_base = NULL;
    f_calc.off_cb();

    //g_lens = NULL;
    std::cout << "VD_SIM: Cleaned up." << std::endl;

    load_flag = false;
    mesh_flag = false;
    first_time = true;

    vd_rem_tag(m);
    m->destroyNative();
    apf::destroyMesh(m);
  }

}

// By considering the possible cell split results, make sure there are enough 
// free cells in the cell_base object.
void vd_sim::set_free_cells() {
  int free_count = 10;
  bool c_flag = false;
  if (c_base->get_free_sz(0) < free_count) {
    c_base->add_free(0, 3*free_count-c_base->get_free_sz(0));
    c_flag = true;
  }
  if (c_base->get_free_sz(1) < free_count) {
    c_base->add_free(1, 3*free_count-c_base->get_free_sz(1));
    c_flag = true;
  }
  if (c_base->get_free_sz(2) < free_count) {
    c_base->add_free(2, 3*free_count-c_base->get_free_sz(2));
    c_flag = true;
  }
/*
  if (c_base->get_free_sz(3) < 3) {
    c_base->add_free(3, 6-c_base->get_free_sz(3));
    c_flag = true;
  }
*/
  if (c_flag)
    reload_mesh();

}

void vd_sim::get_v2move(std::vector<apf::MeshEntity*>& v_2move) {
  int v_count = 0;
  for(int i = 0; i < c_2move.size(); i++) {
    int dim = c_2move.at(i).first;
    int c_id = c_2move.at(i).second - 1;
    v_count = v_count + e_list.e.at(dim).at(c_id).at(0).size();
  }
  v_2move.clear();
  v_2move.reserve(v_count);

  for(int i = 0; i < c_2move.size(); i++) {
    int dim = c_2move.at(i).first;
    int c_id = c_2move.at(i).second - 1;
    v_2move.insert(v_2move.end(), e_list.e.at(dim).at(c_id).at(0).begin(), 
                                  e_list.e.at(dim).at(c_id).at(0).end());
  }
}

void vd_sim::set_c2move(std::vector<std::pair<int, int> >& c_2move_in, 
                                                      bool f_v2move_in) {
  c_2move = c_2move_in;
  f_v2move = f_v2move_in;
}

// opts is a vector of strings in the form ...;d id;...
// where d and id are the dimension and id of the strata to be moved.
void opts2vecintpair(std::vector<std::string>& opts, 
                                      std::vector<std::pair<int,int>>& vipair) {

  vipair.clear();
  vipair.reserve(opts.size());

  for(int i = 0; i < opts.size(); i++) {
    std::string temp("");
    std::string::size_type sz;
    temp = opts.at(i);

    int c_dim = std::stoi(temp, &sz);
    temp = temp.substr(sz);
    int c_id = std::stoi(temp, &sz);
    vipair.push_back(std::make_pair(c_dim, c_id));
  }
}


void vd_sim::set_adapt_file(std::string ad_opts) {
  std::vector<std::vector<std::string> > opts(0, std::vector<std::string > (0, std::string("")));
  ReadNames(ad_opts.c_str(), ";", opts);
  set_adapt_opts(opts.at(0));
}

// Used with adaptation scheme ADAPT_TYPE::ADAPT_BOUND. These cells are refined 
// ADAPT_TYPE::X; ad_flag; ad_len; ad_param; additional_opt1; additional_opt2;
// Examples:
// ADAPT_STEP_1CELL; 1; 1.; 1.; 0; 0; 0;
// ADAPT_STEP_1CELL; 0; 1.; 1.; 0; 1.; 0;
// ADAPT_BOUND; 0; 1.; 1.; ./mshfiles/ref_list.txt; ./mshfiles/move_list.txt;
void vd_sim::set_adapt_opts(std::vector<std::string>& opts) {
  assert(opts.size() == 7);

  // ad_type:
  if(opts.at(0) == "ADAPT_1CELL") {
    set_adapt_type(ADAPT_TYPE::ADAPT_1CELL);
  }
  else if(opts.at(0) == "ADAPT_3CELL") {
    set_adapt_type(ADAPT_TYPE::ADAPT_3CELL);
  }
  else if(opts.at(0) == "ADAPT_BOUND") {
    set_adapt_type(ADAPT_TYPE::ADAPT_BOUND);
  }
  else if(opts.at(0) == "ADAPT_CURVE") {
    set_adapt_type(ADAPT_TYPE::ADAPT_CURVE);
  }
  else if(opts.at(0) == "ADAPT_STEP") {
    set_adapt_type(ADAPT_TYPE::ADAPT_STEP);
  }
  else if(opts.at(0) == "ADAPT_STEP_1CELL") {
    set_adapt_type(ADAPT_TYPE::ADAPT_STEP_1CELL);
  }
  else {
    set_adapt_type(ADAPT_TYPE::ADAPT_STEP_1CELL);
  }
  // ad_flag:
  set_ad_flag(!(opts.at(1).compare(std::string("0")) == 0));
  // ad_ratio:
  if(std::strcmp(opts.at(2).c_str(), "0") != 0) {
    set_adapt(atof(opts.at(2).c_str()));
  }
  else
    set_adapt(1.);
  // ad_param:
  if(opts.at(3).compare(std::string("0")) != 0) {
    set_adapt_param(atof(opts.at(3).c_str()));
  }
  else
    set_adapt_param(1.);
  // cells to refine:
  if(opts.at(4).compare(std::string("0")) != 0) {
    assert(ad_type == ADAPT_TYPE::ADAPT_BOUND);
    std::vector<std::vector<std::string> > temp_str_vec(0,
                                   std::vector<std::string> (0,""));
    std::vector<std::pair<int,int> > cells_ref(0, std::make_pair(0,0));

    ReadNames(opts.at(4).c_str(), ";", temp_str_vec);
    opts2vecintpair(temp_str_vec.at(0), cells_ref);
    set_adapt_bound_cells(std::stod(opts.at(5)), cells_ref);
  }
  // cells to move:
  if(opts.at(6).compare(std::string("0")) != 0) {
    assert(ad_type == ADAPT_TYPE::ADAPT_BOUND);
    std::vector<std::vector<std::string> > temp_str_vec(0,
                                   std::vector<std::string> (0,""));
    std::vector<std::pair<int,int> > cells_move (0, std::make_pair(0,0));

    ReadNames(opts.at(6).c_str(), ";", temp_str_vec);
    opts2vecintpair(temp_str_vec.at(0), cells_move);

    set_c2move(cells_move, true);
  }
}

// Used with adaptation scheme ADAPT_TYPE::ADAPT_BOUND. These cells are refined 
void vd_sim::set_adapt_bound_cells(double len, std::vector<std::pair<int, int> > &cells_ref_in) {
  bound_len = len;
  cells_ref = cells_ref_in;
}

void vd_sim::ref_mesh_iso(double len) {

  double temp = adapt_len;
  adapt_len = len;
  adapt_iso();
  adapt_len = temp;

}

void vd_sim::set_adapt_param(double param_in) {
  adapt_len = param_in;
  std::cout << "adapt_len multiplier is set to " << adapt_len <<"."<< std::endl;
}

void vd_sim::set_adapt_type(ADAPT_TYPE ad_in) {
  ad_type = ad_in;
  std::cout << "adapt_type is " << (int)ad_in << "." << std::endl;
}

void vd_sim::set_field_calc(const field_calc& FC) {

  f_calc = FC;

  if(load_flag)
    f_calc.reload_cb(c_base);

  set_vec_sp_calc((PROJ_TYPE)f_calc.get_proj());
  f_calc.refresh(m, c_base, &e_list);
}

void vd_sim::set_field_calc(VEL_TYPE vel_in) {

  f_calc = field_calc(vel_in);

  if(load_flag)
    f_calc.reload_cb(c_base);

  f_calc.set_vec_sp_calc((PROJ_TYPE)proj_flag);
  f_calc.refresh(m, c_base, &e_list);
}

void vd_sim::set_vec_sp_calc(PROJ_TYPE PROJ) {

  proj_flag = (int)PROJ;
  assert(proj_flag < (int)PROJ_TYPE::END);

  if(proj_flag > (int)PROJ_TYPE::FIXED) {
    calc_ext = true;
  }
  else
    calc_ext = false;

  if(proj_flag > (int)PROJ_TYPE::FIXED and 
     proj_flag < (int)PROJ_TYPE::END)
    calc_corner = true;

  //c_ins.set_calc_corner(calc_corner);

  if(proj_flag == (int)PROJ_TYPE::EXT_SHELL) {
    f_calc.set_vec_sp_calc(PROJ, &e_sh);
    collect_ext();
    c_base->set_ext_spur(true);
  }
  else
    f_calc.set_vec_sp_calc(PROJ);

}

// Set whether the exterior vertices are moving or not.
void vd_sim::set_calc_ext(bool calc_ext_in) {

  calc_ext = calc_ext_in;
  f_calc.set_calc_ext(calc_ext);

}

void vd_sim::set_integ_type(INTEG_TYPE IT) {
  integ_flag = IT;
  f_calc.set_integ_type(IT);
}

// Geometrically detect the side 1cells. Used to collect topology information 
// to preserve exterior geometry.
void vd_sim::collect_ext() {
  assert(mesh_flag);
  shell_burner sb(&e_sh, c_base, m, &e_list);
  sb.collect_ext();
}

void vd_sim::set_time(double time_in) {
  time_curr = time_in;
}

double vd_sim::get_time() {
  return time_curr;
}

void vd_sim::set_dt(double dt_in, bool fixed_in) {
  dt_set = dt_in;
  fixed_time = fixed_in;
}

void vd_sim::set_iter_sub(int sub_nbr) {
  iter_sub = sub_nbr;
}

void vd_sim::set_time_cost_str(const std::string str_in, bool accumulate) {
  time_cost_str = str_in;
  cost_t_accum = accumulate;
}

void vd_sim::set_time(std::vector<std::string>& opts) {
  std::vector<std::string> out(0, std::string(""));
  std::string delim("::");
  for(int i = 0; i < opts.size(); i++) {
    split_str_delim(opts.at(i), delim, out);
    assert(out.size() == 2);
    if(out.at(0) == "TIME_START") {
      time_curr = std::stod(out.at(1));
    }
    else if(out.at(0) == "TIME_STEP") {
      dt_set = std::stod(out.at(1));
    }
    else if(out.at(0) == "TIME_END") {
      time_end = std::stod(out.at(1));
    }
    else if(out.at(0) == "TIME_FIXED") {
      fixed_time = std::stoi(out.at(1));
    }
    else if(out.at(0) == "TIME_SUB_LIMIT") {
      limit_flag = true;
      t_sub_limit = std::stod(out.at(1));
    }
    else if(out.at(0) == "CONDITION_G_NBR") {
      flag_g_nbr = true;
      cond_g_nbr = std::stoi(out.at(1));
    }
    else if(out.at(0) == "CONDITION_S_NBR") {
      flag_s_nbr = true;
      cond_s_nbr = std::stoi(out.at(1));
    }
    else if(out.at(0) == "CONDITION_L_NBR") {
      flag_l_nbr = true;
      cond_l_nbr = std::stoi(out.at(1));
    }
    else if(out.at(0) == "VTK_INTERVAL") {
      set_save_vtk_interval(std::stod(out.at(1)));
    }
  }
}

void vd_sim::write_time_cost() {
  csvfile csvtime(time_cost_str);
  csvtime << time_curr;
  for(int i = 0; i < vd_tim.elapsed.size(); i++)
      csvtime << vd_tim.elapsed.at(i);
  csvtime << endrow;
}

void vd_sim::set_euler_char(bool on_off) {
  isotropic = on_off;
}

void vd_sim::set_ins_flag(bool on_off) {
  ins_flag = on_off;
}

void vd_sim::set_col_flag(bool on_off) {
  col_flag = on_off;
}

void vd_sim::set_col_rat(double ratio_col_sh_in, double ratio_col_in) {
  ratio_col_sh = ratio_col_sh_in;
  ratio_col = ratio_col_in;
}

void vd_sim::set_ins_rat(double ratio_ins_in) {
  ratio_ins = ratio_ins_in;
}

void vd_sim::set_adapt(double param_in) {
  ad_param = param_in;
}

void vd_sim::set_ad_flag(bool on_off) {
  ad_flag = on_off;
}

// Controls adaptations functions that fix shape to actually fix shape or not. By default this should be true.
void vd_sim::set_ad_fix_flag(bool on_off) {
  ad_fix_shape = on_off;
}

void vd_sim::set_ad_th(double coarse_in, double split_in) {
  coarse_th = coarse_in;
  split_th = split_in;
}

/*
void vd_sim::set_extract(bool on_off, std::vector<vd_ext_opts>* opts) {
  extract_flag = on_off;
  //m_opts = *opts;
}

void vd_sim::set_extract_MS(std::string MS_i, std::string ROC_i, 
                            std::string VOL_i, bool on_off) {
  outputMS = MS_i;
  outputROC = ROC_i;
  outputVOL = VOL_i;
  extract_flag = on_off;
}

void vd_sim::set_extract_cells(int i, const std::vector<c_p>* c_top,
                              const std::vector<std::vector<c_p> >* c_bound) {
  //m_opts.at(i).set_extract_cells(c_top, c_bound);
}
*/

// Set equations and energy functions. Multiple entries for same parameter 
// override old ones.
// Default (VEL_TYPE::MASON PROJ_TYPE::FIXED EXT_GAM_TYPE::ZERO)

void vd_sim::set_equations(std::vector<std::string>& opts) {
  bool flag_v = false;
  bool flag_p = false;
  bool flag_eg = false;
  bool flag_int = false;
  bool flag_drag = false;
  double drag_rat = 1000.;

  std::vector<std::string> out(0, std::string(""));
  std::string delim("::");
  for(int i = 0; i < opts.size(); i++) {
    split_str_delim(opts.at(i), delim, out);
    assert(out.size() == 2);
    if(out.at(0) == "VEL_TYPE") {
      set_field_calc(conv_str2vel(out.at(1)));
      flag_v = true;
    }
    else if(out.at(0) == "PROJ_TYPE") {
      set_vec_sp_calc(conv_str2proj(out.at(1)));
      flag_p = true;
    }
    else if(out.at(0) == "EXT_GAM_TYPE") {
      set_ext_gam_flag(conv_str2extgam(out.at(1)));
      flag_eg = true;
    }
    else if(out.at(0) == "INTEG_TYPE") {
      set_integ_type(conv_str2integ(out.at(1)));
      flag_int = true;
    }
    else if(out.at(0) == "DRAG_RAT") {
      flag_drag = true;
      drag_rat = std::stod(out.at(1));
    }
  }
  if(!flag_v) {
    set_field_calc(VEL_TYPE::MASON);
  }
  if(!flag_p) {
    set_vec_sp_calc(PROJ_TYPE::FIXED);
  }
  if(!flag_eg) {
    set_ext_gam_flag(EXT_GAM_TYPE::ZERO);
  }
  if(!flag_int) {
    set_integ_type(INTEG_TYPE::RK2);
  }
  if(flag_drag)
    f_calc.set_drag_rat(drag_rat);
}



// opts is a vector of strings in the form ...;d id;...
// where d and id are the dimension and id of the strata to be moved.
void vd_sim::set_c2move_opts(std::vector<std::string>& opts, 
                                                      bool f_v2move_in) {

  c_2move.clear();
  f_v2move = f_v2move_in;

  c_2move.reserve(opts.size());

  for(int i = 0; i < opts.size(); i++) {
    std::string temp("");
    std::string::size_type sz;
    temp = opts.at(i);

    int c_dim = std::stoi(temp, &sz);
    temp = temp.substr(sz);
    int c_id = std::stoi(temp, &sz);
    c_2move.push_back(std::make_pair(c_dim, c_id));
  }
}


void vd_sim::set_extract(bool on_off, std::string& opts) {
  extract_flag = on_off;
  ext_opts = opts;
}

void vd_sim::extract_data() {
/*
  ext_mesh e_mesh;
  e_mesh.set_mesh(m, c_base);

  for(int i = 0; i < m_opts.size(); i++) {
    e_mesh.set_opts(m_opts.at(i));

    if(m_opts.at(i).f_type == apf::SCALAR) {
      vd_series_ss data;
      //e_mesh.get_cell_bound_ext_angle(&data, 3, 1, &m_opts.c_bound.at(0));
      e_mesh.get_cell_meas(&data, &m_opts.at(i).c_top);
      std::stringstream ss;
      ss << time_curr*10000;

      data.title = "time";
      data.variable = ss.str();
      if(m_opts.at(i).first_ext) {
        write_whitespace_x_vert(&data, "./output/data.csv");
        m_opts.at(i).first_ext = false;
      }
      write_label_y_vert(&data, "./output/data.csv");
      //write_data(&data, tmp.c_str() );
    }
  }
*/
  vd_ext_inf v_ext_inf(m, c_base, &e_list, &f_calc, ext_opts, time_curr);
}
/*
void vd_sim::extract_MS() {
  ext_MS e_MS(m, c_base);
  e_MS.process_mesh();

  csvfile csvMS(outputMS);
  csvfile csvROC(outputROC);
  csvfile csvVOL(outputVOL);

  std::vector<double> roc(c_base->get_sz(3), 0);
  std::vector<apf::MeshEntity*> tet_bound(0);
  csvMS << time_curr;
  csvROC << time_curr;
  csvVOL << time_curr;

  std::cout << "t: " <<  time_curr << std::endl;
  for(int j = 0; j < c_base->get_sz(3); j++) {
    if(!c_base->is_free(3,j)) {
      tet_bound = vd_get_3c_tet_bound(m, e_list.e.at(3).at(j).at(3));
      roc.at(j) = vd_calc_roc(m, &tet_bound, &f_calc);
      std::cout << "3c" << j + 1 
                << " MS: " << e_MS.dVdt_c3[j+1]*f_calc.vdparam.v_mult 
                << " roc: " << roc.at(j) << std::endl; 
      csvMS << e_MS.dVdt_c3[j+1]*f_calc.vdparam.v_mult;
      csvROC << roc.at(j);
      csvVOL << vd_meas_set(m, &e_list.e.at(3).at(j).at(3));
    }
    else {
      csvMS << 0;
      csvROC << 0;
      csvVOL << 0;
    }
  }
  csvMS << endrow;
  csvROC << endrow;
  csvVOL << endrow;
}
*/
vd_param vd_sim::get_vd_param() {
  return f_calc.vdparam;
}

void vd_sim::set_vd_param(vd_param vd_par) {
  f_calc.set_vdparam(vd_par, m);
  //f_calc.vdparam = vd_par;
}

void vd_sim::set_ext_gam_flag(EXT_GAM_TYPE ext_in) {
  ext_gam_flag = ext_in;
}

field_calc* vd_sim::get_f_calc() {
  return &f_calc;
}

double vd_sim::get_sim_len() {
  return sim_len;
}

double vd_sim::get_len_topo() {
  return len_topo;
}

vd_entlist* vd_sim::get_elist() {
  return &e_list;
}

void vd_sim::refresh_elist() {
  e_list.refresh();
}

// Print the 2-3cell adjacencies around the 0cells. 
void vd_sim::print_0cellgraph() {
//  c_ins.print_graph();
}

void vd_sim::print_0cell_pc() {
//  c_ins.print_circ();
//  c_ins.print_path();
}

vd_sim::~vd_sim() {
//  tag_0cell_ins.clear();
  for (int i = 0; i < meas_below_map.size(); i++) {
    meas_below_map.at(i).clear();
  }
  for (int i = 0; i < t_below_map.size(); i++) {
    t_below_map.at(i).clear();
  }
  ins_rad_map.clear();
}

// Return true if shell membership disallows collapse.
bool vd_sim::chk_shell_issue(int d, int cell_id) {
  bool disallowed = false;
  if(tag_shell_chk.at(d-1)[cell_id+1]) {
    return rad_map_cc.at(d-1)[cell_id+1] or tag_shell_map.at(d-1)[cell_id+1];
  }
  if(proj_flag == (int) PROJ_TYPE::EXT_SHELL) {
    // Assuming no spurious cells exists, check if the exterior shell
    // membership prevents the collapse.
    // For each pair of same dim shells, check if the connecting shell
    // is also reprensented. Assume connectivity (No worm-like cells).

    // For lower dim cells, there should be a single lowest dim shell
    // that bounds all other shells. No need to check the current cell.
    std::vector< std::vector<int> > shell_list(3, 
                                            std::vector<int>(0));

    // Get the shells of each cell belonging to a shell.
    if(e_sh.chk_shell(d, cell_id)) {
      shell sh_curr = e_sh.get_shell(d, cell_id);
      shell_list.at(sh_curr.dim).push_back(sh_curr.id);
    }
    for(int dim_low = 0; dim_low < d; dim_low++) {
      ent_conn c_low;
      c_base->get_conn_dim(dim_low, d, cell_id, &c_low);
      for(int i = 0; i < c_low.conn.size(); i++) {
        if(e_sh.chk_shell(dim_low, c_low.conn.at(i))) {
          shell sh_curr = e_sh.get_shell(dim_low, c_low.conn.at(i));
          shell_list.at(sh_curr.dim).push_back(sh_curr.id);
        }
      }
    }
    // Get the unique shells.
    std::vector<int>::iterator it;
    for(int i = 0; i < 3; i++) {
      std::sort(shell_list.at(i).begin(), shell_list.at(i).end());
      it = std::unique(shell_list.at(i).begin(), shell_list.at(i).end());
      shell_list.at(i).resize(
                    std::distance(shell_list.at(i).begin(),it));
    }
    // For 2shells, check each pair for 1shell intersections.
    // If exists, it should be in the 1shell list. 
    // Otherwise, the collapse is not allowed.
    if(shell_list.at(2).size() > 1) {
      for(int i = 0; i < shell_list.at(2).size(); i++) {
        ent_conn sh_low1;
        e_sh.sh_base.get_conn_gmi(2, shell_list.at(2).at(i), &sh_low1);

        for(int j = i+1; j < shell_list.at(2).size(); j++) {
          ent_conn sh_low2;
          e_sh.sh_base.get_conn_gmi(2, shell_list.at(2).at(j), &sh_low2);
          for(int k = 0; k < sh_low2.conn.size(); k++) {
            it = std::find(sh_low1.conn.begin(), sh_low1.conn.end(), 
                                                    sh_low2.conn.at(k));
            if(it != sh_low1.conn.end()) {
              int sh_1 = *it;
              it = std::find(shell_list.at(1).begin(), 
                              shell_list.at(1).end(), sh_1);
              if(it == shell_list.at(1).end()) {
                std::cout << d << "c" << cell_id + 1 << std::endl;
                std::cout << "2shells " << shell_list.at(2).at(i) << " "
                          << shell_list.at(2).at(j)
                          << " are not connected and collapse invalid. " 
                          << std::endl;
                disallowed = true;
                i = shell_list.at(2).size();
                j = i;
                k = sh_low2.conn.size();
              }
            }
          }
        }
      }
    }
    if(shell_list.at(1).size() > 1) {
      for(int i = 0; i < shell_list.at(1).size(); i++) {
        ent_conn sh_low1;
        e_sh.sh_base.get_conn_gmi(1, shell_list.at(1).at(i), &sh_low1);

        for(int j = i+1; j < shell_list.at(1).size(); j++) {
          ent_conn sh_low2;
          e_sh.sh_base.get_conn_gmi(1, shell_list.at(1).at(j), &sh_low2);
          for(int k = 0; k < sh_low2.conn.size(); k++) {
            it = std::find(sh_low1.conn.begin(), sh_low1.conn.end(), 
                                                    sh_low2.conn.at(k));
            if(it != sh_low1.conn.end()) {
              int sh_0 = *it;
              it = std::find(shell_list.at(0).begin(), 
                              shell_list.at(0).end(), sh_0);
              if(it == shell_list.at(0).end()) {
                std::cout << d << "c" << cell_id + 1 << std::endl;
                std::cout << "1shells " << shell_list.at(1).at(i) << " "
                          << shell_list.at(1).at(j)
                          << " are not connected and collapse invalid. " 
                          << std::endl;
                disallowed = true;
                i = shell_list.at(1).size();
                j = i;
                k = sh_low2.conn.size();
              }
            }
          }
        }
      }
    }
  }
  tag_shell_chk.at(d-1)[cell_id+1] = true;
  tag_shell_map.at(d-1)[cell_id+1] = disallowed;

  return disallowed;
}

// Find cells below the containing cell radius threshold.
// TODO Some of these functions do the same sweep multiple times internally 
// to collect the bounding cell entities. They can be collected for efficiency.
std::pair<int, int> vd_sim::chk_small_cell_rad() {
  std::pair<int,int> c2c;

  double average_ent = ma::getAverageEdgeLength(m);

  std::cout << "Average, edge: " << average_ent << std::endl;
  // Attach the Runge Kutta fields.
  //f_calc.vd_att_rk4(m);

  for (int d = 3; d > 0; d--) {
    int d_sz = c_base->get_sz(d);
    for (int cell_id = 0; cell_id < d_sz; cell_id++) {
      //if((c_base->cell_sz(d, cell_id) > 0) ) {
      std::vector<ent_conn> adj_up(3-d);
      for (int up_dim = d+1; up_dim < 4; up_dim++) {
        c_base->get_conn_dim(up_dim, d, cell_id, &adj_up.at(up_dim-d-1));
      }

      if(!c_base->is_free(d, cell_id)) {
        rad_map_cc.at(d-1)[cell_id+1] = chk_shell_issue(d, cell_id);
        if(!rad_map.at(d-1)[cell_id+1] 
           and !rad_map_th.at(d-1)[cell_id+1] 
           and !rad_map_cc.at(d-1)[cell_id+1] 
           and chk_cell_rmv(d, cell_id+1)) {
          if(d < 3) {
            c2c = chk_cell_rmv_up(d, cell_id+1);
            if(c2c.first > -1) {
              if(!rad_map_th.at(c2c.first-1)[c2c.second]) {
                assert(rad_map_cc.at(c2c.first-1)[c2c.second]);
              }
              if(!rad_map_cc.at(c2c.first-1)[c2c.second]) {
                return c2c;
              }
              else {
                std::cout << "Upward adjacency cannot be collapsed." << std::endl;
              }
            }
          }

          c2c.first = d;
          c2c.second = cell_id+1;
          return c2c;
        }
      }
    }
  }
  c2c.first = 0;
  c2c.second = 0;

  return c2c;
}

std::pair<int, int> vd_sim::chk_cell_rmv_up(int d, int cell_id) {
  // The candidate cell should have an equivalent radius less than this multiple
  // of the collapsing cell.
  double up_rat = 1.4;

  ent_conn e_up;
  std::vector<apf::MeshEntity*> ev(0);
  std::vector<apf::MeshEntity*> ev_d(0);
  std::vector<apf::MeshEntity*> ev_e(0);
  std::vector<apf::MeshEntity*> ev_orig(0);

  std::map<apf::MeshEntity*, apf::Vector3> pos_old{};
  std::map<apf::MeshEntity*, apf::MeshEntity*> v_map{};

  std::vector<apf::MeshEntity*> ent_temp(0);
  std::vector<apf::MeshEntity*> vert_temp(0);

  vd_vert_set(m, &e_list.e.at(d).at(cell_id-1).at(d), &ev_orig);

  std::vector<apf::MeshEntity*>* es = &e_list.e.at(d).at(cell_id-1).at(0);
  // TODO This bit is going to be replaced with weighted PCA.
  //if(c_base->get_conn_sz(d, cell_id-1) > 1)
  //  vd_remove_set(&ev_orig, es);

  assert(ev_orig.size() > 0);
  int i_curr = 0;
  shell sh_min = shell(3,0);
  if(proj_flag == (int) PROJ_TYPE::EXT_SHELL) {
    for(int i = 0; i < ev_orig.size(); i++) {
      apf::ModelEntity* mdl = m->toModel(ev_orig.at(i));
      int em_type = m->getModelType(mdl);
      int em_tag = m->getModelTag(mdl);
      if(em_type < 3 and e_sh.chk_shell(em_type, em_tag - 1)) {
        shell sh_temp = e_sh.get_shell(em_type, em_tag - 1);
        if(sh_temp.dim < sh_min.dim) {
          sh_min = sh_temp;
          i_curr = i;
        }
      }
    }
  }
  else {
    apf::ModelEntity* mdl = m->toModel(ev_orig.at(0));
    int em_min = m->getModelType(mdl);
    for(int i = 0; i < ev_orig.size(); i++) {
      mdl = m->toModel(ev_orig.at(i));
      int em_type = m->getModelType(mdl);
      int em_tag = m->getModelTag(mdl);
      if(em_type < 3 and em_type < em_min) {
        i_curr = i;
        em_min = em_type;
      }
    }
  }
  for(int i = 0; i < ev_orig.size(); i++) {
    v_map[ev_orig.at(i)] = ev_orig.at(i_curr);
  }

  apf::Vector3 temp(0,0,0);

  for(int i = 0; i < ev_orig.size(); i++) {
    m->getPoint(ev_orig.at(i), 0, temp);
    pos_old[ev_orig.at(i)] = temp;
  }
  apf::Vector3 c_ctr = vd_get_center_e(m, &e_list.e.at(d).at(cell_id-1).at(d));
  if(proj_flag == (int) PROJ_TYPE::EXT_SHELL) {
    if(sh_min.dim < 3)
      c_ctr = e_sh.find_int_pos(sh_min, c_ctr);
  }

  for(int i = 0; i < ev_orig.size(); i++) {
    m->setPoint(ev_orig.at(i), 0, c_ctr);
  }

  std::cout << "Checking upward adjacencies of " << d << "c" << cell_id
            << " for collapse after collapsing the cell."
            << std::endl;
  bool found = false;
  int d_found = -1;
  int c_found = -1;
  for(int dim_up = 3; dim_up > d; dim_up--) {
    c_base->get_conn_dim_gmi(dim_up, d, cell_id, &e_up);
    for(int i = 0; i < e_up.conn.size(); i++) {
      int c_curr = e_up.conn.at(i);
      if(!rad_map.at(dim_up-1)[e_up.conn.at(i)]
         and !rad_map_cc.at(dim_up-1)[e_up.conn.at(i)]) {
        ent_temp = e_list.e.at(dim_up).at(c_curr-1).at(dim_up);
        vd_set_down(m, &ent_temp, &vert_temp, dim_up);
        double roc = f_calc.calc_roc_merg(m, &ent_temp, vert_temp, v_map);

        std::cout << "Max dimension less than threshold." << std::endl;
        std::cout << dim_up << "-cell" << c_curr << " rate of change " << roc 
                  << " r_equi " << sm_map.at(dim_up-1)[c_curr] << std::endl;
        std::cout << "dt: " << f_calc.vdparam.dt
                  << " %rate: " << roc/meas_map.at(dim_up-1)[c_curr]
                  << std::endl;
        //bool shrinking = (roc*f_calc.vdparam.dt/cell_meas_max < roc_tol);
        if(roc < -std::numeric_limits<double>::min() and
           sm_map.at(dim_up-1)[c_curr] < up_rat*sm_map.at(d-1)[cell_id]) {
          std::cout << dim_up << "-cell" << c_curr << " collapses with "
                    << d << "-cell" << cell_id << "."
                    << std::endl;

          found = true;
          d_found = dim_up;
          c_found = e_up.conn.at(i);
          i = e_up.conn.size();
          dim_up = d;
        }
      } 
    }

  }

  for(int i = 0; i < ev_orig.size(); i++) {
    m->setPoint(ev_orig.at(i), 0, pos_old[ev_orig.at(i)]);
  }
  if(found)
    return std::make_pair(d_found, c_found);
  else
    return std::make_pair(-1, -1);
}

apf::Vector3 vd_sim::get_cell_midpoint(int d, int cell_id) {
  // Get the midpoint
  std::vector<int> ext_c0(0);
  bool ext_cell = false;
  struct ent_conn e0;
  c_base->get_conn_dim_gmi(0, d, cell_id, &e0);

  ext_c0.clear();
  ext_c0.reserve(e0.conn.size());
  for(int i = 0; i < e0.conn.size(); i++) {
    if(c_base->get_cell_ext_gmi(0, e0.conn.at(i)) ) {
      ext_cell = true;
      ext_c0.push_back(e0.conn.at(i));
    }
  }
  std::vector<apf::Vector3> ext_c0_vec(0, apf::Vector3(0,0,0));
  ext_c0_vec.resize(ext_c0.size());

  for(int i = 0; i < ext_c0.size(); i++) {
    std::vector<apf::MeshEntity*>* e_pt = 
                                  &e_list.e.at(0).at(ext_c0.at(i)-1).at(0);
    assert(e_pt->size() == 1);
    m->getPoint(e_pt->at(0), 0, ext_c0_vec.at(i));
  }

  apf::Vector3 midpoint(0,0,0);
  if(ext_cell) {
    assert(ext_c0_vec.size() > 0);
    for(int i = 0; i < ext_c0_vec.size(); i++) {
      midpoint = midpoint + ext_c0_vec.at(i);
    }
    midpoint = midpoint/ext_c0_vec.size();
  }
  else
    midpoint = vd_get_center_e(m, &e_list.e.at(d).at(cell_id-1).at(d));
  return midpoint;
}


void vd_sim::clear_rad_map() {

  for (int i = 0; i < rad_map.size(); i++) {
    rad_map.at(i).clear();
  }
  for (int i = 0; i < rad_map_cc.size(); i++) {
    rad_map_cc.at(i).clear();
  }
  for (int i = 0; i < rad_map_th.size(); i++) {
    rad_map_th.at(i).clear();
  }
  for (int i = 0; i < sm_map.size(); i++) {
    sm_map.at(i).clear();
  }
  for (int i = 0; i < meas_map.size(); i++) {
    meas_map.at(i).clear();
  }
  //rad_map.clear();
}

// Given a cell, check if collapsing collapses a higher dim adjacency.
// If the collapse causes a disallowed topological transition, return <-1,-1>
// If the collapse doesn't cause any problems return the same id.
// If an upper dimensional cell can be collapsed, return it's id.

// It is allowed, if merging vertices collapse an edge belonging to an upper
// dim adjacency, if the collapse doesn't cause any additional topological 
// changes.

// This is in place, as with 2 edges per 1cell, meshadapt can create tets
// bounded by both edges. On a triangle bounded by both edges, the other edge  
// on the triangle will remain in the convex hull, and will also be removed.
std::pair<int, int> vd_sim::chk_col_cell_adj(int d, int cell_id) {

  std::vector<std::vector<apf::MeshEntity*> > ent(4, 
                                  std::vector<apf::MeshEntity*>(0));
  std::vector<std::map<int,bool> > col_map(3, std::map<int, bool> {});
 
  //vd_find_ent_geom(m, &ent.at(cell_col_dim), cell_col_id, cell_col_dim);
  ent.at(d) = e_list.e.at(d).at(cell_id-1).at(d);

  struct ent_conn* e_down = new ent_conn();
  for(int dim = d-1; dim > 0; dim--) {
    c_base->get_conn_dim_gmi(dim, d, cell_id, e_down);
    for(int i = 0; i < e_down->conn.size(); i++) {
      col_map.at(dim - 1)[e_down->conn.at(i)] = true;
    }
  }
  delete e_down;

  for(int dim = d; dim > 0; dim--) {
    vd_set_down(m, &ent.at(dim), &ent.at(dim-1));
  }

  std::vector<apf::MeshEntity*> ee(0);
  vd_set_up(m, &ent.at(0), &ee);

  apf::Downward d_v;
  for (int i = 0; i < ee.size(); i++) {
    int c_type = m->getModelType(m->toModel(ee.at(i)));
    int c_tag = m->getModelTag(m->toModel(ee.at(i)));
    if(!col_map.at(c_type - 1)[c_tag]) {
      m->getDownward(ee.at(i), 0, d_v);
      int i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[0]);
      if(i1 == -1) {
        i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[1]);
        assert(i1 > -1);
      }
      else {
        i1 = findIn(&ent.at(0), ent.at(0).size(), d_v[1]);
        if(i1 != -1) {
          std::cout << "Cannot collapse, entity outside cell bounded by merging vertices." << std::endl;
          return std::make_pair(-1, -1);
        }
      }
    }
  }

  for(int dim = 0; dim < 3; dim++) {
    vd_set_up(m, &ent.at(dim), &ent.at(dim+1));
  }
  return std::make_pair(d, cell_id);
}

// Check cell collapse based on change in energy associated with collapse.
// Collapse the cell in a trial mesh if a singular transition is caused by
// the collapse. 
// If negative, return true. Don't consider inverting elements, for now.

bool vd_sim::chk_cell_rmv(int d, int cell_id) {
  // If a cell has r_equi less than len_trans/rat_col, collapse it if it is 
  // shrinking.
  //double ins_len = len_trans/topo_rat;

  //if(c_base->get_conn_sz(d, cell_id-1) == 1 and d == 1)
  //  return false;
  double roc = f_calc.calc_roc(m, &e_list.e.at(d).at(cell_id-1).at(d));
  // Take the weighted average of ROC since reduced below threshold and current 
  // ROC multiplied by current time step divided by total time to yield weighted 
  // average ROC. 
  std::cout << d << "c" << cell_id 
            << " r_equi " << sm_map.at(d-1)[cell_id]
            << " roc before " << roc << std::endl;

  if(t_below_map.at(d-1)[cell_id] > std::numeric_limits<double>::min()) {
    roc = (meas_map.at(d-1)[cell_id] + roc * f_calc.vdparam.dt - meas_below_map.at(d-1)[cell_id])/(f_calc.vdparam.dt + time_curr - t_below_map.at(d-1)[cell_id]);
  }
  else {
    meas_below_map.at(d-1)[cell_id] = meas_map.at(d-1)[cell_id];
    t_below_map.at(d-1)[cell_id] = time_curr;
  }

  std::cout << "roc after " << roc << std::endl;
  std::cout << "dt: " << f_calc.vdparam.dt
            << " %rate: " << roc/meas_map.at(d-1)[cell_id]
            << std::endl;
  //bool shrinking = (roc*f_calc.vdparam.dt/cell_meas_max < roc_tol);
  bool shrinking = false;
  if(sm_map.at(d-1)[cell_id] < len_trans/ratio_col) {
    shrinking = (roc < -std::numeric_limits<double>::min());
  }
  else {
    shrinking = (roc < -std::numeric_limits<double>::min()
                           and chk_cell_shr(d, cell_id));
  }
  rad_map_th.at(d-1)[cell_id] = !shrinking;
  return (shrinking);
}

void vd_sim::set_mov_flag(bool mov_flag) {
  save_mov = mov_flag;
}

void vd_sim::set_sub_vtk_flag(bool vtk_flag) {
  sub_vtk = vtk_flag;
}

