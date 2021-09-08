#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include "topo_field.h"
#include "topo_manip.h"
#include "topo_disp.h"


VEL_TYPE conv_str2vel(const std::string& str) {
  if(str == "LAZAR") return VEL_TYPE::LAZAR;
  else if(str == "MASON") return VEL_TYPE::MASON;
  else if(str == "MASON_NBC") return VEL_TYPE::MASON_NBC;
  else if(str == "MASON_MIR") return VEL_TYPE::MASON_MIR;
  else if(str == "KUPRAT") return VEL_TYPE::KUPRAT;
  else if(str == "KUPRAT_NBC") return VEL_TYPE::KUPRAT_NBC;
  return VEL_TYPE::END;
}

INTEG_TYPE conv_str2integ(const std::string& str) {
  if(str == "EULER") return INTEG_TYPE::EULER;
  else if(str == "RK2") return INTEG_TYPE::RK2;
  else if(str == "RK4") return INTEG_TYPE::RK4;
  return INTEG_TYPE::END;
}

EXT_GAM_TYPE conv_str2extgam(const std::string& str) {
  if(str == "CONSTANT") return EXT_GAM_TYPE::CONSTANT;
  else if(str == "ZERO") return EXT_GAM_TYPE::ZERO;
  return EXT_GAM_TYPE::END;
}

PROJ_TYPE conv_str2proj(const std::string& str) {
  if(str == "FIXED") return PROJ_TYPE::FIXED;
  else if(str == "LOCAL_TRI") return PROJ_TYPE::LOCAL_TRI;
  else if(str == "LOCAL_SP") return PROJ_TYPE::LOCAL_SP;
  else if(str == "FLAT") return PROJ_TYPE::FLAT;
  else if(str == "EXT_SHELL") return PROJ_TYPE::EXT_SHELL;
  return PROJ_TYPE::END;
}

shell::shell() : dim(0), id(1) {
}

shell::shell(int dim_in, int id_in) : dim(dim_in), id(id_in) {
}

shell::~shell() {
}

shell::shell(const shell& that) : dim(that.dim), id(that.id) {
}

shell& shell::operator=(const shell& that) {
  dim = that.dim;
  id = that.id;
  return *this;
}

bool shell::operator==(const shell& that) const{
  return (dim == that.dim and id == that.id);
}

bool shell::operator<(const shell& that) const{
  if(dim != that.dim)
    return dim < that.dim;
  return (id < that.id);
}

// ext_shell used to preserve volume by a naive approach.
bool ext_shell::chk_shell(int dim, int tag) {
  if(dim == 0)
    return c0_bool[tag];

  else if(dim == 1)
    return c1_bool[tag];

  else if(dim == 2)
    return c2_bool[tag];
  return false;
}

shell ext_shell::get_shell(int dim, int tag) {
  assert(chk_shell(dim, tag));
  if(dim == 0)
    return c0_cor[tag];

  else if(dim == 1)
    return c1_side[tag];

  else if(dim == 2)
    return c2_face[tag];
}

void ext_shell::set_shell(int dim, int tag, shell sh, bool state) {
  //assert(tag > 0);
  //assert(dim < 3);
  assert(dim < 3 and dim > -1);
  //std::cout << dim << "c" << tag+1 << " setting shell " 
  //          << sh.dim << "sh" << sh.id << " " << state
  //          << std::endl;
  if(dim == 0) {
    c0_cor[tag] = sh;
    c0_bool[tag] = state;
  }

  else if(dim == 1) {
    c1_side[tag] = sh;
    c1_bool[tag] = state;
  }

  else if(dim == 2) {
    c2_face[tag] = sh;
    c2_bool[tag] = state;
  }
}

void ext_shell::set_pos(shell sh, apf::Vector3 pos) {
  if(sh.dim == 0)
    c0_cor_pos[sh] = pos;

  else if(sh.dim == 1)
    c1_side_pos[sh] = pos;

  else if(sh.dim == 2)
    c2_face_pos[sh] = pos;
  assert(sh.dim < 3 and sh.dim > -1);
}

void ext_shell::set_dir(shell sh, apf::Vector3 dir) {
  assert(sh.dim < 3 and sh.dim > 0);
  dir = norm_0(dir);
  if(sh.dim == 1)
    c1_side_dir[sh] = dir;

  else if(sh.dim == 2)
    c2_face_dir[sh] = dir;
}

apf::Vector3 ext_shell::get_shell_pos(shell sh) {
  assert(sh.dim < 3 and sh.dim > -1);
  if(sh.dim == 0)
    return c0_cor_pos[sh];

  else if(sh.dim == 1)
    return c1_side_pos[sh];

  else if(sh.dim == 2)
    return c2_face_pos[sh];
}

apf::Vector3 ext_shell::get_shell_dir(shell sh) {

  assert(sh.dim < 3 and sh.dim > 0);
  if(sh.dim == 1)
    return c1_side_dir[sh];

  else if(sh.dim == 2)
    return c2_face_dir[sh];
}

apf::Vector3 ext_shell::get_pos(int dim, int tag) {
  assert(chk_shell(dim, tag));
  return get_shell_pos(get_shell(dim, tag));
}

apf::Vector3 ext_shell::get_dir(int dim, int tag) {
  assert(chk_shell(dim, tag));
  return get_shell_dir(get_shell(dim, tag));
}

apf::Vector3 ext_shell::find_int_pos(int dim, int tag, apf::Vector3 pos) {
  assert(chk_shell(dim, tag));
  shell sh = get_shell(dim, tag);

  //std::cout << dim << "c" << tag+1 << " " 
  //          << sh.dim << "sh" << sh.id
  //          << get_shell_pos(sh) << " "
  //          << std::endl;
  //if(sh.dim > 0)
  //  std::cout << get_shell_dir(sh) << std::endl;

  if(sh.dim == 0) {
    return c0_cor_pos[sh];
  }
  // Move along the line.
  else if(sh.dim == 1) {
    apf::Vector3 n_curr = c1_side_dir[sh];
    return c1_side_pos[sh] + n_curr*
                           (n_curr*(pos - c1_side_pos[sh] ) );
  }
  // Move towards the plane.
  else {
    assert(sh.dim == 2);
    apf::Vector3 n_curr = c2_face_dir[sh];
    return pos + n_curr*(n_curr*(c2_face_pos[sh] - pos) );
  }
}

apf::Vector3 ext_shell::find_int_pos(shell sh, apf::Vector3 pos) {
  assert(sh.dim > -1 and sh.dim < 3);

  //std::cout << dim << "c" << tag+1 << " " 
  //          << sh.dim << "sh" << sh.id
  //          << get_shell_pos(sh) << " "
  //          << std::endl;
  //if(sh.dim > 0)
  //  std::cout << get_shell_dir(sh) << std::endl;

  if(sh.dim == 0) {
    return c0_cor_pos[sh];
  }
  // Move along the line.
  else if(sh.dim == 1) {
    apf::Vector3 n_curr = c1_side_dir[sh];
    return c1_side_pos[sh] + n_curr*
                           (n_curr*(pos - c1_side_pos[sh] ) );
  }
  // Move towards the plane.
  else {
    assert(sh.dim == 2);
    apf::Vector3 n_curr = c2_face_dir[sh];
    return pos + n_curr*(n_curr*(c2_face_pos[sh] - pos) );
  }
}

// Find the component along the boundary. If perpendicular, return zero vector.
apf::Vector3 ext_shell::find_para_dir(int dim, int tag, apf::Vector3 dir) {
  assert(chk_shell(dim, tag));
  shell sh = get_shell(dim, tag);

  apf::Vector3 vec(0,0,0);

  if(sh.dim == 0) {
    return vec;
  }
  else if(sh.dim == 1) {
    apf::Vector3 n_curr = c1_side_dir[sh];    
    return n_curr*(dir*n_curr);
  }
  // Move towards the plane.
  else {
    assert(sh.dim == 2);
    apf::Vector3 n_curr = c2_face_dir[sh];    
    return dir - n_curr*(dir*n_curr);
  }
}

apf::Vector3 ext_shell::find_int_traj(int dim, int tag, apf::Vector3 pos, 
                                                       apf::Vector3 dir) {
  assert(chk_shell(dim, tag));
  shell sh = get_shell(dim, tag);

  // Direction doesn't matter. To remain on the shell, project to the point.
  if(sh.dim == 0) {
    return c0_cor_pos[sh] - pos;
  }
  // If dir is parallel to the shell boundary, return the direction from the  
  // position, normal to the boundary.
  else if(sh.dim == 1) {
    apf::Vector3 temp(0,0,0);
    temp = c1_side_pos[sh] - pos;

    apf::Vector3 n_curr = c1_side_dir[sh];
    dir = norm_0(dir);
    if(std::fabs(n_curr*dir - 1) < std::numeric_limits<double>::epsilon() ) {
      return temp - n_curr*(temp*n_curr);
    }
    else {
      return dir*(temp*n_curr)/(dir*n_curr);
    }
  }
  else {
    assert(sh.dim == 2);
    apf::Vector3 n_curr = c2_face_dir[sh];
    dir = norm_0(dir);
    if(std::fabs(n_curr*dir) < std::numeric_limits<double>::epsilon() ) {
      dir = n_curr;
      return dir*(n_curr*(c2_face_pos[get_shell(dim, tag)] - pos));
    }
    else
      return dir*(n_curr*(c2_face_pos[get_shell(dim, tag)] - pos))/(dir*n_curr);
  }
}

// Given a position, return the normal displacement from the position to the 
// shell.
apf::Vector3 ext_shell::find_norm_disp(shell sh, apf::Vector3 pos) {
  // Direction doesn't matter. To remain on the shell, project to the point.
  if(sh.dim == 0) {
    return c0_cor_pos[sh] - pos;
  }
  // If dir is parallel to the shell boundary, return the direction from the  
  // position, normal to the boundary.
  else if(sh.dim == 1) {
    apf::Vector3 temp(0,0,0);
    temp = c1_side_pos[sh] - pos;

    apf::Vector3 n_curr = c1_side_dir[sh];
    return temp - n_curr*(temp*n_curr);
  }
  else {
    assert(sh.dim == 2);
    apf::Vector3 n_curr = c2_face_dir[sh];
    return n_curr*(n_curr*(c2_face_pos[sh] - pos));
  }
}

// Given a position, return the position mirrored by the shell.
apf::Vector3 ext_shell::find_pos_mir(int dim, int tag, apf::Vector3 pos) {
  assert(chk_shell(dim, tag));
  shell sh = get_shell(dim, tag);
  apf::Vector3 rel_pos(0,0,0);
  if(sh.dim == 0) {
    return c0_cor_pos[sh] - (c0_cor_pos[sh] - pos)*2;
  }
  // If dir is parallel to the shell boundary, return the direction from the  
  // position, normal to the boundary.
  else if(sh.dim == 1) {
    rel_pos = c1_side_pos[sh] - pos;
    rel_pos = rel_pos - c1_side_dir[sh]*(c1_side_dir[sh]*rel_pos);
    return c1_side_pos[sh] - rel_pos*2;
  }
  else {
    assert(sh.dim == 2);
    rel_pos = c2_face_pos[sh] - pos;

    rel_pos = c2_face_dir[sh]*(c2_face_dir[sh]*rel_pos);
    return c1_side_pos[sh] - rel_pos*2;
  }
}

// Given a direction, return the direction mirrored by the shell.
apf::Vector3 ext_shell::find_dir_mir(int dim, int tag, apf::Vector3 dir) {
  assert(chk_shell(dim, tag));
  shell sh = get_shell(dim, tag);
  apf::Vector3 n_dir(0,0,0);
  // Invert
  if(sh.dim == 0) {
    return dir*(-1);
  }
  // If dir is parallel to the shell boundary, return the direction from the  
  // position, normal to the boundary.
  else if(sh.dim == 1) {
    n_dir = c1_side_dir[sh]*(c1_side_dir[sh]*dir);
    return n_dir - (dir - n_dir)*2;
  }
  else {
    assert(sh.dim == 2);
    n_dir = c2_face_dir[sh]*(c2_face_dir[sh]*dir);
    return dir - n_dir*2;
  }
}

apf::Vector3 ext_shell::find_dir_mir(shell sh, apf::Vector3 dir) {
  apf::Vector3 n_dir(0,0,0);
  // Invert
  if(sh.dim == 0) {
    return dir*(-1);
  }
  // If dir is parallel to the shell boundary, return the direction from the  
  // position, normal to the boundary.
  else if(sh.dim == 1) {
    n_dir = c1_side_dir[sh]*(c1_side_dir[sh]*dir);
    return n_dir - (dir - n_dir)*2;
  }
  else {
    assert(sh.dim == 2);
    n_dir = c2_face_dir[sh]*(c2_face_dir[sh]*dir);
    return dir - n_dir*2;
  }
}


int ext_shell::find_new_cell_id(int dim, int tag, int tag_0cell) {
/*
//get the upward adjacencies of all the new 0-1 cells, removing the new cells 
//these will have the candidate shells

  // Downward adjacencies
  std::vector<ent_conn> e_d;
  std::vector<std::vector<ent_conn> > e_up;
  e_d.resize(dim);
  e_up.resize(dim);

  std::map<int, int> cor_0_map;


  //If a 0cell is exterior 

  //If a 0cell has multiple side 1cell upward adjacency, it belongs to the 
  //corner 0shell

  //If a 0cell has a single side 1cell upward adjacency, it belongs to the 
  //side 1shell

  for(int i = 0; i < dim; i++) {
    e_lens.c_base_curr->get_conn_dim_gmi(i, dim, tag, &e_d.at(i));

    e_up.at(i).resize(e_d.at(i).size());
    for(int j = 0; j < e_d.at(i).conn.size(); j++) {
      e_lens.c_base_curr->get_conn_dim_gmi(i+1, i, e_d.at(i).conn.at(j), 
                                                     &e_up.at(i).at(j));
    }
  }

  if(dim == 1) {
    for(int j = 0; j < e_up.at(0).size(); j++) {
      e_up.at(0).at(j).rem_ent(tag);
    }
  }
  else {
    assert(dim == 2);
    for(int i = 0; i < e_up.at(0).size(); i++) {
      for(int j = 0; j < e_d.at(1).size(); j++) {
        e_up.at(0).at(i).rem_ent(e_d.at(1).conn.at(j));
      }
    }

    for(int j = 0; j < e_up.at(1).size(); j++) {
      e_up.at(1).at(j).rem_ent(tag);
    }
  }

  for(int i = 0; i < e_d.at(0).conn.size(); i++) {
    if(e_lens.c_base_curr->get_cell_ext_gmi(e_d.at(0).conn.at(i)) ) {
      int cor_nbr = 0;
      int cor_id;
      for(int j = 0; j < e_up.at(0).at(i).size(); j++) {
        if(e_lens.c_base_curr->get_1c_corner(e_up.at(0).at(i).at(j) - 1) ) {
          cor_nbr = cor_nbr + 1;
          cor_id = e_up.at(0).at(i).at(j);
        }
      }

      if(cor_nbr == 0) {
        bool found = false;
        for(int j = 0; j < e_up.at(0).at(i).size(); j++) {
          if(e_sh->chk_shell(1, e_up.at(0).at(i).at(j) - 1) ) {
            shell sh = e_sh->get_shell(1, e_up.at(0).at(i).at(j) - 1);
            if(sh.dim == 2) {
              e_sh->set_shell(0, e_d.at(0).conn.at(i) - 1, sh);
              j = e_up.at(0).at(i).size();
              found = true;
            }
          }
        }
        assert(found);
      }
      else if(cor_nbr == 1) {
        shell sh = e_sh->get_shell(1, cor_id - 1);
        e_sh->set_shell(0, e_d.at(0).conn.at(i) - 1, sh);
      }
      else {
        shell sh = e_sh->get_shell(0, tag_0cell - 1);
        e_sh->set_shell(0, e_d.at(0).conn.at(i) - 1, sh);
      }

    }
    else {
      shell sh = e_sh->get_shell(0, tag_0cell - 1);
      e_sh->set_shell(0, e_d.at(0).conn.at(i) - 1, sh, false);
    }
  }

// If no side 1cell upward adjacency, check the 1cell adjecencies. It will belong
//to the shell of one of the exterior 1cells. They should all belong to the same
//shell. 

// If a 1cell is a side 1cell, check the bounding 0cells. If one is side 1shell, the 1cell belongs to the side 1shell

// Else, assert both do not belong to corner 0shell

//If a 1cell is a exterior 1cell, check the upper 2cell adjacencies. It will 
// belong to the shell of one of the exterior 2cells. They should all belong to 
// the same shell. 

  if(dim == 2) {
    for(int i = 0; i < e_d.at(1).conn.size(); i++) {
      ent_conn e0;
      e_lens.c_base_curr->get_down_gmi(1, e_d.at(1).conn.at(i), &e0);
      if(e_lens.c_base_curr->get_cell_ext_gmi(1, e_d.at(1).conn.at(i)) ) {
        shell sh1 = e_sh->get_shell(0, e0.conn.at(0) - 1);
        shell sh2 = e_sh->get_shell(0, e0.conn.at(1) - 1);
        if(sh1.dim > sh2.dim) {
          e_sh->set_shell(1, e_d.at(1).at(i).at(j) - 1, sh1);
        }
        else {
          e_sh->set_shell(1, e_d.at(1).at(i).at(j) - 1, sh2);
        }
      }
      else {
        shell sh = e_sh->get_shell(0, tag_0cell - 1);
        e_sh->set_shell(1, e_d.at(1).at(i).at(j) - 1, sh, false);
      }
    }
  }
*/
}


void ext_shell::clear() {

  c0_bool.clear();
  c1_bool.clear();
  c2_bool.clear();

  c0_cor.clear();
  c1_side.clear();
  c2_face.clear();

  c0_cor_pos.clear();
  c1_side_pos.clear();
  c2_face_pos.clear();

  c1_side_dir.clear();
  c2_face_dir.clear();
}

ext_shell::ext_shell() :
    c0_bool{}, c1_bool{}, c2_bool{},
    c0_cor{}, c1_side{}, c2_face{},
    c0_cor_pos{}, c1_side_pos{}, c2_face_pos{}, 
    c1_side_dir{}, c2_face_dir{},
    sh_base(0,0,0,0) {
  clear();
}

ext_shell::~ext_shell() {
//memory overflow either while constructing or destructing. the sh_base variable is not initialized correctly, it seems
  clear();
}

ext_shell::ext_shell(const ext_shell& that) :
    c0_bool{}, c1_bool{}, c2_bool{},
    c0_cor{}, c1_side{}, c2_face{},
    c0_cor_pos{}, c1_side_pos{}, c2_face_pos{}, 
    c1_side_dir{}, c2_face_dir{},
    sh_base(0,0,0,0) {
  clear();

  c0_bool = that.c0_bool;
  c1_bool = that.c1_bool;
  c2_bool = that.c2_bool;

  c0_cor = that.c0_cor;
  c1_side = that.c1_side;
  c2_face = that.c2_face;

    // The positions of the shell boundaries.
  c0_cor_pos = that.c0_cor_pos;
  c1_side_pos = that.c1_side_pos;
  c2_face_pos = that.c2_face_pos;

    // The directions of the shell sides and faces.
  c1_side_dir = that.c1_side_dir;
  c2_face_dir = that.c2_face_dir;

  sh_base = that.sh_base;
}

ext_shell& ext_shell::operator=(const ext_shell& that) {
  clear();
  c0_bool = that.c0_bool;
  c1_bool = that.c1_bool;
  c2_bool = that.c2_bool;

  c0_cor = that.c0_cor;
  c1_side = that.c1_side;
  c2_face = that.c2_face;

    // The positions of the shell boundaries.
  c0_cor_pos = that.c0_cor_pos;
  c1_side_pos = that.c1_side_pos;
  c2_face_pos = that.c2_face_pos;

    // The directions of the shell sides and faces.
  c1_side_dir = that.c1_side_dir;
  c2_face_dir = that.c2_face_dir;

  sh_base = that.sh_base;
  return *this;
}


shell_burner::shell_burner(ext_shell* e_sh_in, cell_base* cb_in, 
                           apf::Mesh2* m_in, vd_entlist* e_in) {
  e_sh = e_sh_in;
  c_base = cb_in;
  m = m_in;
  e_list = e_in;
}

// Burn the shell. Used to burn 1 and 2cells.
void shell_burner::burn_shell(int dim, int tag, shell sh) {

  assert(c_base->get_cell_ext(dim, tag));
  ent_conn* edown = new ent_conn();
  ent_conn* esame = new ent_conn();

  e_sh->set_shell(dim, tag, sh);
  c_base->get_conn(dim, tag, edown);

  for(int i = 0; i < edown->conn.size(); i++) {
    assert(c_base->get_cell_ext(dim - 1, edown->conn.at(i)) );
    if(!e_sh->chk_shell(dim - 1, edown->conn.at(i)) ) {
      e_sh->set_shell(dim - 1, edown->conn.at(i), sh);
      c_base->get_conn_dim(dim, dim - 1, edown->conn.at(i), esame);
      for(int j = 0; j < esame->conn.size(); j++) {
        // Along a certain side, there can be no bifurcations.
        bool c1_found = false;
        if(tag != esame->conn.at(j)) {
          if(dim == 1) {
            if(c_base->get_1c_corner(esame->conn.at(j)) and
               !e_sh->chk_shell(dim, esame->conn.at(j)) ) {
              assert(!c1_found);
              burn_shell(dim, esame->conn.at(j), sh);
              c1_found = true;
            }
          }
          else if(c_base->get_cell_ext(dim, esame->conn.at(j)) ) {
            if(!e_sh->chk_shell(dim, esame->conn.at(j)) ) {
              burn_shell(dim, esame->conn.at(j), sh);
            }
          }
        }
      }
    }
    else {
      shell sh_down = e_sh->get_shell(dim - 1, edown->conn.at(i));
      assert(sh_down.dim < sh.dim 
             or (sh_down.dim == sh.dim and sh_down.id == sh.id) );
      if(sh_down.dim == sh.dim - 1) {
        ent_conn sh_conn;
        e_sh->sh_base.get_conn_gmi(sh.dim, sh.id, &sh_conn);
        if(sh_conn.find_ent(sh_down.id) < 0) {
          sh_conn.add_ent(sh_down.id);
          e_sh->sh_base.set_conn_gmi(sh.dim, sh.id, &sh_conn);
        }
      }
    }
  }

  c_base->get_conn_dim(0, dim, tag, edown);
  for(int i = 0; i < edown->conn.size(); i++) {
    if(!e_sh->chk_shell(0, edown->conn.at(i)) ) {
      e_sh->set_shell(0, edown->conn.at(i), sh);
    }
  }

  delete edown;
  delete esame;

}

// Geometrically detect the side 1cells. Used to collect topology information 
// to preserve exterior geometry.
void shell_burner::collect_ext(double tol) {

  e_list->refresh();
  assert(e_list->e.size() > 0);

  c_base->clear_cor();
  e_sh->clear();

  std::map<int, apf::Vector3> norm_map{};
  std::map<int, int> cor_1c_count{};

  struct ent_conn e0;
  struct ent_conn e2;

  double z[3] = {0,0,0};
  apf::Vector3 n_vec;
  apf::Vector3 temp;
  n_vec.fromArray(z);

  for(int i = 0; i < c_base->get_sz(2); i++) {
    if(!c_base->is_free(2,i) and c_base->get_cell_ext(2, i)) {
      //std::cout << "Exterior 2c" << i+1 << " " 
      //          << e_list->e.at(2).at(i).size() << " "
      //          << std::endl;
      assert(e_list->e.at(2).at(i).size() > 0 and
         e_list->e.at(2).at(i).at(2).size() > 0);
      apf::MeshEntity* ent_curr = e_list->e.at(2).at(i).at(2).at(0);
      norm_map[i] = norm_0(vd_area_out(m, ent_curr));
    }
    //  else
    //    norm_map[i] = n_vec;
    //}
    e_sh->set_shell(2, i, shell(2, 1), false);
  }

  for(int i = 0; i < c_base->get_sz(1); i++) {
    c_base->set_cor_ext(i, false);
    e_sh->set_shell(1, i, shell(1, 1), false);
  }

  for(int i = 0; i < c_base->get_sz(0); i++) {
    cor_1c_count[i] = 0;
    e_sh->set_shell(0, i, shell(0, 1), false);
  }

  // If an exterior 1cell has less than 2 0cell adjacencies, do not consider 
  // it for exterior shell detection.
  std::map<int, bool> c1_skip{};

  for(int i = 0; i < c_base->get_sz(1); i++) {
    if(!c_base->is_free(1, i) and c_base->get_cell_ext(1, i)) {
      c_base->get_conn(1, i, &e0);

      if(e0.conn.size() < 2)
        c1_skip[i] = true;
      else {
        temp.fromArray(z);
        c_base->get_conn_dim(2, 1, i, &e2);
        if(e2.conn.size() > 0) {
          bool found = false;
          for(int j = 0; j < e2.conn.size(); j++) {
            if(c_base->get_cell_ext(2, e2.conn.at(j))) {
              if(!found) {
                temp = norm_map[e2.conn.at(j)];
                found = true;
              }
              else {
                temp = temp * (temp * norm_map[e2.conn.at(j)]);
              }
            }
          }
          assert(found);
          // If all exterior 2cells are in plane, face 1cell. Otherwise 
          // corner.
          //std::cout << "Exterior 1c" << i+1;
          if(!(std::fabs(temp.getLength()-1.) <= 
                      tol) ) {
//                      std::numeric_limits<double>::epsilon() ) ) {
            c_base->set_cor_ext(i, true);
            c_base->get_conn(1, i, &e0);
            for(int j = 0; j < e0.conn.size(); j++) {
              cor_1c_count[e0.conn.at(j)] = cor_1c_count[e0.conn.at(j)] + 1;
            }
            //std::cout << ", corner: temp " << temp << std::endl;
          }
          //else {
            //std::cout << ", plane: temp " << temp << std::endl;          
          //}
        }
      }
    }
  }

  e_sh->sh_base.clear();

  shell sh(0, 1);
  for(int i = 0; i < c_base->get_sz(0); i++) {
    //std::cout << "0c" << i+1 << " " << e_sh->chk_shell(0, i)
    //          << " cor_1c " << cor_1c_count[i];
    if(cor_1c_count[i] > 2) {
      e_sh->sh_base.add_free(0, 1);
      int sh_free = e_sh->sh_base.use_free_gmi(0);
      assert(sh_free == sh.id);

      //std::cout << " " << sh.dim << "sh" << sh.id;
      e_sh->set_shell(0, i, sh);
      sh.id = sh.id + 1;
    }
    //std::cout << std::endl;
  }

  sh.dim = 1;
  sh.id = 1;

  for(int i = 0; i < c_base->get_sz(1); i++) {
    if(c_base->get_1c_corner(i) ) {
      //std::cout << "Corner 1c" << i+1 << " ";
      if(!e_sh->chk_shell(1, i)) {
        e_sh->sh_base.add_free(1, 1);
        int sh_free = e_sh->sh_base.use_free_gmi(1);
        assert(sh_free == sh.id);

        burn_shell(1, i, sh);
        sh.id = sh.id + 1;
      }
      //else {
      //  std::cout << " burnt" << std::endl;
      //}
    }
  }

  sh.dim = 2;
  sh.id = 1;

  for(int i = 0; i < c_base->get_sz(2); i++) {
    if(c_base->get_cell_ext(2, i) and !e_sh->chk_shell(2, i)) {
      e_sh->sh_base.add_free(2, 1);
      int sh_free = e_sh->sh_base.use_free_gmi(2);
      assert(sh_free == sh.id);

      //std::cout << "Ext 2c" << i+1 << " " << sh.dim << "sh" << sh.id 
      //          << std::endl;
      burn_shell(2, i, sh);
      sh.id = sh.id + 1;
    }
  }

  e_sh->sh_base.add_free(3, 1);
  int sh_free = e_sh->sh_base.use_free_gmi(3);
  assert(sh_free == 1);

  ent_conn sh_conn;

  e_sh->sh_base.get_conn_gmi(3, 1, &sh_conn);
  for(int i = 0; i < e_sh->sh_base.get_sz(2); i++) {
    sh_conn.add_ent(i+1);
  }
  e_sh->sh_base.set_conn_gmi(3, 1, &sh_conn);

  collect_ext_dir();

}

void shell_burner::collect_ext_dir() {

  struct ent_conn e_adj;
  apf::Vector3 t_pos;
  apf::Vector3 t_dir;
  // If an exterior 1cell has less than 2 0cell adjacencies, do not consider 
  // it for exterior shell detection.
  std::map<int, bool> c1_skip{};

  // Assert no used exterior 1c with no connection, or a single bounding 0c.
  for(int i = 0; i < c_base->get_sz(1); i++) {
    if(!c_base->is_free(1, i) and c_base->get_cell_ext(1, i)) {
      c_base->get_conn(1, i, &e_adj);

      if(e_adj.conn.size() < 2)
        c1_skip[i] = true;
    }
  }

  // Burn the side 1c
  // Collect the line directions

  std::map<int, bool> shell_map{};
  for(int i = 0; i < c_base->get_sz(0); i++) {
    if(!c_base->is_free(0, i) and c_base->get_cell_ext(0, i) and e_sh->chk_shell(0, i)) {
      shell sh = e_sh->get_shell(0, i);
      if(sh.dim == 0 and !shell_map[sh.id]) {
        shell_map[sh.id] = true;
        assert(e_list->e.at(0).at(i).at(0).size() == 1);
        apf::MeshEntity* ent_curr = e_list->e.at(0).at(i).at(0).at(0);

        m->getPoint(ent_curr, 0, t_pos);
        std::cout << "0c" << i+1 << " 0sh" << sh.id
                  << " pos " << t_pos << std::endl;
        e_sh->set_pos(sh, t_pos);
      }
    }
  }
  shell_map.clear();

  for(int i = 0; i < c_base->get_sz(1); i++) {
    if(c_base->get_cell_ext(1, i) and e_sh->chk_shell(1, i)) {
      shell sh = e_sh->get_shell(1, i);
      if(sh.dim == 1 and !shell_map[sh.id]) {
        shell_map[sh.id] = true;
        assert(e_list->e.at(1).at(i).at(1).size() > 0);

        apf::MeshEntity* ent_curr = e_list->e.at(1).at(i).at(1).at(0);
        //apf::Downward d_v;

        //m->getDownward(ent_curr, 0, d_v);
        //m->getPoint(d_v[0], 0, t_pos);
        t_pos = vd_get_pos(m, ent_curr);
        e_sh->set_pos(sh, t_pos);

        t_dir = get_edge_dir(m, ent_curr);
        e_sh->set_dir(sh, t_dir);
        assert(t_dir.getLength() > std::numeric_limits<double>::epsilon());
        std::cout << "1c" << i+1 << " 1sh" << sh.id
                  << " pos " << t_pos 
                  << " dir " << t_dir 
                  << std::endl;

      }
    }
  }

  shell_map.clear();

  for(int i = 0; i < c_base->get_sz(2); i++) {
    if(c_base->get_cell_ext(2, i) and e_sh->chk_shell(2, i)) {
      shell sh = e_sh->get_shell(2, i);
      if(sh.dim == 2 and !shell_map[sh.id]) {
        shell_map[sh.id] = true;
        assert(e_list->e.at(2).at(i).at(2).size() > 0);

        apf::MeshEntity* ent_curr = e_list->e.at(2).at(i).at(2).at(0);
        //apf::Downward d_v;
        //m->getDownward(ent_curr, 0, d_v);
        //m->getPoint(d_v[0], 0, t_pos);
        t_pos = vd_get_pos(m, ent_curr);
        e_sh->set_pos(sh, t_pos);

        t_dir = norm_0(vd_area_out(m, e_list->e.at(2).at(i).at(2).at(0)));

        assert(t_dir.getLength() > std::numeric_limits<double>::epsilon());
        e_sh->set_dir(sh, t_dir);
        std::cout << "2c" << i+1 << " 2sh" << sh.id
                  << " pos " << t_pos 
                  << " dir " << t_dir 
                  << std::endl;
      }
    }
  }
  // Burn the exterior surfaces
  // Collect the surface normals
}

////////////////////////////////////////////////////

void boundary_energy::set_2cells(std::map<int, bool> c2s_in) {
  c2s_zeros = c2s_in;
}

boundary_energy::boundary_energy() : c2s_zeros{} {
}

boundary_energy::~boundary_energy() {
}

boundary_energy& boundary_energy::operator=(const boundary_energy& that) {
  c2s_zeros = that.c2s_zeros;
  return *this;
}

////////////////////////////////////////////////////

double field_calc::calc_roc(apf::Mesh2* m, 
                                    std::vector<apf::MeshEntity*>* ent) {
  assert(ent->size() > 0);
  int ent_type = m->getType(ent->at(0));
  int d = m->typeDimension[ent_type];

  assert(d > 0 and d < 4);
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  std::vector<apf::MeshEntity*> ents(0);
  ents = *ent;
  // Remove singular volume tets.
  if(d == 3) {
    int count = 0;
    for(int i = 0; i < ents.size() - count; i++) {
      if(vd_volume_tet(m, ents.at(i)) < std::numeric_limits<double>::min()) {
        apf::MeshEntity* temp = ents.at(ents.size()-count-1);
        ents.at(ents.size()-count-1) = ents.at(i);
        ents.at(i) = temp;
        count = count + 1;
        i = i - 1;
      }
    }
    ents.resize(ents.size()-count);
  }

  std::vector<apf::MeshEntity*> verts(0);
  vd_set_down(m, &ents, &verts, d);
  //vd_calc_rk4(m, vd_par, verts);

  //for(int i = 0; i < verts.size(); i++) {
  //  vd_upd_vel_field(m, verts.at(i));
  //}
  vd_upd_vel_field(m, &verts);

  return vd_calc_roc(m, &ents, this);
}

double field_calc::calc_roc_merg(apf::Mesh2* m, 
     std::vector<apf::MeshEntity*>* ent, std::vector<apf::MeshEntity*>& verts,
     std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map) {
  //assert(v_flag == VEL_TYPE::MASON);

  assert(ent->size() > 0);
  int ent_type = m->getType(ent->at(0));
  int d = m->typeDimension[ent_type];

  assert(d > 0 and d < 4);
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  vd_upd_vel_field_merg(m, verts, v_map);

//don't consider ents which have 0 area/volume, i.e. ents that are bounded by more than 1 merging vertices, so remove those from ent
  //std::map<apf::MeshEntity*, bool> col_map{};
  apf::Downward d_v;

  if(d > 1) {
    int count = 0;
    for(int i = 0; i < ent->size() - count; i++) {
      m->getDownward(ent->at(i), 0, d_v);
      bool found = false;
      for(int j = 0; j < d+1; j++) {
        if(v_map[d_v[j]] != NULL) {
          if(found) {
            //col_map[ent->at(i)] = true;
            apf::MeshEntity* temp = ent->at(ent->size()-count-1);
            ent->at(ent->size()-count-1) = ent->at(i);
            ent->at(i) = temp;
            count = count + 1;
            //if(i + count < ent->size() - 1)
            i = i - 1;
            //else
            //  i = ent->size();
            j = d+1;
          }
          else
            found = true;
        }
      }
    }
    ent->resize(ent->size()-count);
  }
  if(d == 3) {
    int count = 0;
    for(int i = 0; i < ent->size() - count; i++) {
      if(vd_volume_tet(m, ent->at(i)) < std::numeric_limits<double>::min()) {
        apf::MeshEntity* temp = ent->at(ent->size()-count-1);
        ent->at(ent->size()-count-1) = ent->at(i);
        ent->at(i) = temp;
        count = count + 1;
        i = i - 1;
      }
    }
    ent->resize(ent->size()-count);
  }

  // If no non-collapsing entity remains, return negative rate of change.
  if(ent->size() == 0)
    return -1;
  return vd_calc_roc(m, ent, this);
}

double field_calc::calc_roc(apf::Mesh2* m, apf::MeshEntity* ent) {
  int ent_type = m->getType(ent);
  int d = m->typeDimension[ent_type];

  assert(d > 0 and d < 4);

  apf::Downward d_v;

  double roc = 0;

  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  m->getDownward(ent, 0, d_v);

  for(int j = 0; j < d+1; j++) {
    vd_upd_vel_field(m, d_v[j], true);
  }
  return vd_calc_roc(m, ent, this);
}

apf::Vector3 field_calc::vd_calc_force(apf::Mesh2* m, apf::MeshEntity* v, bool fix) {
  apf::Vector3 force(0,0,0);
  force = EoM->vd_calc_force(v);
  if (fix and chk_vert_special(m, v)) {
    force = get_vec_special(m, v, force);
  }
  return force;
}

apf::Vector3 field_calc::vd_calc_force_tri(apf::Mesh2* m, apf::MeshEntity* v, std::vector<apf::MeshEntity*>* tris, bool fix) {
  apf::Vector3 force(0,0,0);
  force = EoM->vd_calc_force_tri(v, tris);
  if (fix and chk_vert_special(m, v)) {
    force = get_vec_special(m, v, force);
  }
  return force;
}

apf::Vector3 field_calc::vd_calc_force_edge(apf::Mesh2* m, apf::MeshEntity* v, std::vector<apf::MeshEntity*>* tris, bool fix) {
  apf::Vector3 force(0,0,0);
  force = EoM->vd_calc_force_edge(v, tris);
  if (fix and chk_vert_special(m, v)) {
    force = get_vec_special(m, v, force);
  }
  return force;
}

void field_calc::vd_upd_vel_field(apf::Mesh2* m, apf::MeshEntity* v, 
                                  bool drag_local) {
  if(t_int_flag == INTEG_TYPE::EULER)
    EoM->vd_upd_vel_field(v);
  else if(t_int_flag == INTEG_TYPE::RK2) {
    vd_calc_rk2(m, v, drag_local);
  }
  else {
    assert(t_int_flag == INTEG_TYPE::RK4);
    vd_calc_rk4(m, v, drag_local);
  }
  if (chk_vert_special(m, v)) {
    fix_vel_special(m, v);
  }
}

void field_calc::vd_upd_vel_field(apf::Mesh2* m, 
                                  std::vector<apf::MeshEntity*>* v, 
                                  bool drag_local) {
  if(t_int_flag == INTEG_TYPE::EULER)
    for(int i = 0; i < v->size(); i ++)
      EoM->vd_upd_vel_field(v->at(i));
  else if(t_int_flag == INTEG_TYPE::RK2) {
    vd_calc_rk2(m, v, drag_local);
  }
  else {
    assert(t_int_flag == INTEG_TYPE::RK4);
    vd_calc_rk4(m, v, drag_local);
  }
  for(int i = 0; i < v->size(); i ++) {
    if (chk_vert_special(m, v->at(i)))
      fix_vel_special(m, v->at(i));
  }
}

void field_calc::vd_upd_vel_field_merg(apf::Mesh2* m, 
                      std::vector<apf::MeshEntity*> &v, 
                      std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map, 
                      bool drag_local) {
  vd_upd_vel_field_mason_merg(m, v, v_map, this, drag_local);
  fix_vel_special_merg(m, v, v_map);
}

void field_calc::vd_upd_vel_field(apf::Mesh2* m, apf::MeshEntity* v, 
                  std::vector<apf::MeshEntity*>* tri_set, bool fix) {
  EoM->vd_upd_vel_field_tri(v, tri_set);
  if (fix and chk_vert_special(m, v)) {
    fix_vel_special(m, v);
  }
}

void field_calc::vd_upd_vel_field_edge(apf::Mesh2* m, apf::MeshEntity* v, 
                  std::vector<apf::MeshEntity*>* edge_set, bool fix) {
  EoM->vd_upd_vel_field_edge(v, edge_set);
  if (fix and chk_vert_special(m, v)) {
    fix_vel_special(m, v);
  }
}

void field_calc::vd_calc_vel(apf::Mesh2* m) {
  if(t_int_flag == INTEG_TYPE::EULER)
    EoM->calc_vel();
  else if(t_int_flag == INTEG_TYPE::RK2) {
    vd_calc_rk2(m);
  }
  else {
    assert(t_int_flag == INTEG_TYPE::RK4);
    vd_calc_rk4(m);
  }
/*
  apf::MeshEntity* v;
  apf::MeshIterator* it_e = m->begin(0);
  while (v = m->iterate(it_e)) {  
    if (chk_vert_special(m, v))
      fix_vel_special(m, v);
  }
  m->end(it_e);
*/
}

void field_calc::vd_calc_vel(apf::Mesh2* m, std::vector<apf::MeshEntity*>* verts) {
  if(t_int_flag == INTEG_TYPE::EULER)
    EoM->calc_vel(verts, false);
  else if(t_int_flag == INTEG_TYPE::RK2) {
    vd_calc_rk2(m, verts, false);
  }
  else {
    assert(t_int_flag == INTEG_TYPE::RK4);
    vd_calc_rk4(m, verts, false);
  }
}

// Using Runge Kutta2 method, calculate the position and velocity of the 
// vertices. 
void field_calc::vd_calc_rk2(apf::Mesh2* m) {

  apf::Field* rk1 = vd_att_vv_field(m, "rk1");
  assert(rk1);
  apf::Field* rk2 = vd_att_vv_field(m, "rk2");
  assert(rk2);
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  double h = vdparam.dt;
  // y0
  apf::Vector3 temp(0,0,0);
  apf::Vector3 point(0,0,0);

  apf::MeshEntity* v;
  apf::MeshIterator* it_e = m->begin(0);
  while (v = m->iterate(it_e)) {
    m->getPoint(v, 0, temp);
    apf::setVector(rk2, v, 0, temp);
  }
  m->end(it_e);

  EoM->calc_vel();
  it_e= m->begin(0);
  while (v = m->iterate(it_e)) {
    apf::getVector(vel_field, v, 0, temp);
    apf::Vector3 k1 = temp*h;                       // k1
    apf::setVector(rk1, v, 0, k1);
    apf::getVector(rk2, v, 0, point);
    m->setPoint(v, 0, point + k1/2);          // y1
  }
  m->end(it_e);

  EoM->calc_vel();

  it_e = m->begin(0);
  while (v = m->iterate(it_e)) {
    apf::getVector(vel_field, v, 0, temp);
    apf::Vector3 k2 = temp*h;                       // k4

    apf::getVector(rk2, v, 0, point);
    m->setPoint(v, 0, point);            // y_n
    point = k2/h;
    assert(!std::isnan(point.getLength()));
    apf::setVector(vel_field, v, 0, point);
  }
  m->end(it_e);
  //k1 = h*v0(y0)
  //y1 = y_0 + k1/2
  //y1 = y_0 + v0*h/2
  //k2 = h*v1(y1)
  //y2 = y_0 + k2/2
  //y2 = y_0 + v1*h/2
  //k3 = h*v2(y2)
  //y3 = y_0 + k3
  //y3 = y_0 + v2*h
  //k4 = h*v3(y3)
  //y_n = y0 + (k1+2*k2+2*k3+k4)/6
}

// Using Runge Kutta4 method, calculate the position and velocity of the 
// vertices. 
void field_calc::vd_calc_rk4(apf::Mesh2* m) {

  apf::Field* rk1 = vd_att_vv_field(m, "rk1");
  assert(rk1);
  apf::Field* rk2 = vd_att_vv_field(m, "rk2");
  assert(rk2);
  apf::Field* rk3 = vd_att_vv_field(m, "rk3");
  assert(rk3);
  apf::Field* rk4 = vd_att_vv_field(m, "rk4");
  assert(rk4);
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  double h = vdparam.dt;
  // y0
  apf::Vector3 temp(0,0,0);
  apf::Vector3 point(0,0,0);

  apf::MeshEntity* v;
  apf::MeshIterator* it_e = m->begin(0);
  while (v = m->iterate(it_e)) {
    m->getPoint(v, 0, temp);
    apf::setVector(rk4, v, 0, temp);
  }
  m->end(it_e);

  EoM->calc_vel();
  it_e= m->begin(0);
  while (v = m->iterate(it_e)) {
    apf::getVector(vel_field, v, 0, temp);
    apf::Vector3 k1 = temp*h;                       // k1
    apf::setVector(rk1, v, 0, k1);
    apf::getVector(rk4, v, 0, point);
    m->setPoint(v, 0, point + k1/2);          // y1
  }
  m->end(it_e);

  EoM->calc_vel();

  it_e= m->begin(0);
  while (v = m->iterate(it_e)) {
    apf::getVector(vel_field, v, 0, temp);
    apf::Vector3 k2 = temp*h;                       // k2
    apf::setVector(rk2, v, 0, k2);
    apf::getVector(rk4, v, 0, point);
    m->setPoint(v, 0, point + k2/2);          // y2
  }
  m->end(it_e);

  EoM->calc_vel();

  it_e= m->begin(0);
  while (v = m->iterate(it_e)) {
    apf::getVector(vel_field, v, 0, temp);
    apf::Vector3 k3 = temp*h;                       // k3
    apf::setVector(rk3, v, 0, k3);
    apf::getVector(rk4, v, 0, point);
    m->setPoint(v, 0, point + k3/2);          // y3
  }
  m->end(it_e);

  EoM->calc_vel();

  it_e = m->begin(0);
  while (v = m->iterate(it_e)) {
    apf::getVector(vel_field, v, 0, temp);
    apf::Vector3 k4 = temp*h;                       // k4

    //temp = (k1+2*k2+2*k3+k4)/6;
    apf::getVector(rk1, v, 0, temp);
    apf::getVector(rk2, v, 0, point);
    temp = temp + point*2;
    apf::getVector(rk3, v, 0, point);
    temp = temp + point*2 + k4;
    temp =  temp/6;

    apf::getVector(rk4, v, 0, point);
    m->setPoint(v, 0, point);            // y_n
    point = temp/h;
    assert(!std::isnan(point.getLength()));
    apf::setVector(vel_field, v, 0, point);
  }
  m->end(it_e);
  //k1 = h*v0(y0)
  //y1 = y_0 + k1/2
  //y1 = y_0 + v0*h/2
  //k2 = h*v1(y1)
  //y2 = y_0 + k2/2
  //y2 = y_0 + v1*h/2
  //k3 = h*v2(y2)
  //y3 = y_0 + k3
  //y3 = y_0 + v2*h
  //k4 = h*v3(y3)
  //y_n = y0 + (k1+2*k2+2*k3+k4)/6
}

// Using Runge Kutta2 method, calculate the position and velocity of the 
// vertices. 
void field_calc::vd_calc_rk2(apf::Mesh2* m, std::vector<apf::MeshEntity*>* v, 
                                      bool drag_local) {

  apf::Field* rk1 = m->findField("rk1");
  assert(rk1);
  apf::Field* rk2 = m->findField("rk2");
  assert(rk2);
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  double h = vdparam.dt;
  // y0
  apf::Vector3 temp(0,0,0);
  apf::Vector3 point(0,0,0);

  for(int i = 0; i < v->size(); i++) {
    m->getPoint(v->at(i), 0, temp);
    apf::setVector(rk2, v->at(i), 0, temp);
  }
  EoM->calc_vel(v, drag_local);

  for(int i = 0; i < v->size(); i++) {
    apf::getVector(vel_field, v->at(i), 0, temp);
    apf::Vector3 k1 = temp*h;                       // k1
    apf::setVector(rk1, v->at(i), 0, k1);
    apf::getVector(rk2, v->at(i), 0, point);
    m->setPoint(v->at(i), 0, point + k1/2);          // y1
  }

  EoM->calc_vel(v, drag_local);

  for(int i = 0; i < v->size(); i++) {
    apf::getVector(vel_field, v->at(i), 0, temp);
    apf::Vector3 k2 = temp*h;                       // k4

    apf::getVector(rk2, v->at(i), 0, point);
    m->setPoint(v->at(i), 0, point);            // y_n
    apf::setVector(vel_field, v->at(i), 0, k2/h);
  }
  for(int i = 0; i < v->size(); i++) {
    if(chk_vert_special(m, v->at(i))) {
      fix_vel_special(m, v->at(i));
    }
  }
}

// Using Runge Kutta4 method, calculate the position and velocity of the 
// vertices. 
void field_calc::vd_calc_rk4(apf::Mesh2* m, std::vector<apf::MeshEntity*> * v, 
                                        bool drag_local) {

  apf::Field* rk1 = m->findField("rk1");
  assert(rk1);
  apf::Field* rk2 = m->findField("rk2");
  assert(rk2);
  apf::Field* rk3 = m->findField("rk3");
  assert(rk3);
  apf::Field* rk4 = m->findField("rk4");
  assert(rk4);
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  double h = vdparam.dt;
  // y0
  apf::Vector3 temp(0,0,0);
  apf::Vector3 point(0,0,0);

  for(int i = 0; i < v->size(); i++) {
    m->getPoint(v->at(i), 0, temp);
    apf::setVector(rk4, v->at(i), 0, temp);
  }
  EoM->calc_vel(v, drag_local);

  for(int i = 0; i < v->size(); i++) {
    apf::getVector(vel_field, v->at(i), 0, temp);
    apf::Vector3 k1 = temp*h;                       // k1
    apf::setVector(rk1, v->at(i), 0, k1);
    apf::getVector(rk4, v->at(i), 0, point);
    m->setPoint(v->at(i), 0, point + k1/2);          // y1
  }

  EoM->calc_vel(v, drag_local);

  for(int i = 0; i < v->size(); i++) {
    apf::getVector(vel_field, v->at(i), 0, temp);
    apf::Vector3 k2 = temp*h;                       // k2
    apf::setVector(rk2, v->at(i), 0, k2);
    apf::getVector(rk4, v->at(i), 0, point);
    m->setPoint(v->at(i), 0, point + k2/2);          // y2
  }

  EoM->calc_vel(v, drag_local);

  for(int i = 0; i < v->size(); i++) {
    apf::getVector(vel_field, v->at(i), 0, temp);
    apf::Vector3 k3 = temp*h;                       // k3
    apf::setVector(rk3, v->at(i), 0, k3);
    apf::getVector(rk4, v->at(i), 0, point);
    m->setPoint(v->at(i), 0, point + k3/2);          // y3
  }

  EoM->calc_vel(v, drag_local);

  for(int i = 0; i < v->size(); i++) {
    apf::getVector(vel_field, v->at(i), 0, temp);
    apf::Vector3 k4 = temp*h;                       // k4

    //temp = (k1+2*k2+2*k3+k4)/6;
    apf::getVector(rk1, v->at(i), 0, temp);
    apf::getVector(rk2, v->at(i), 0, point);
    temp = temp + point*2;
    apf::getVector(rk3, v->at(i), 0, point);
    temp = temp + point*2 + k4;
    temp =  temp/6;

    apf::getVector(rk4, v->at(i), 0, point);
    m->setPoint(v->at(i), 0, point);            // y_n
    apf::setVector(vel_field, v->at(i), 0, temp/h);
  }
  for(int i = 0; i < v->size(); i++) {
    if(chk_vert_special(m, v->at(i))) {
      fix_vel_special(m, v->at(i));
    }
  }
}

void field_calc::vd_calc_rk2(apf::Mesh2* m, apf::MeshEntity* v, 
                                            bool drag_local) {

  apf::Field* rk1 = m->findField("rk1");
  assert(rk1);
  apf::Field* rk2 = m->findField("rk2");
  assert(rk2);
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  double h = vdparam.dt;
  // y0
  apf::Vector3 temp(0,0,0);
  apf::Vector3 point(0,0,0);

  m->getPoint(v, 0, temp);
  apf::setVector(rk2, v, 0, temp);
  EoM->vd_upd_vel_field(v, drag_local);

  apf::getVector(vel_field, v, 0, temp);
  apf::Vector3 k1 = temp*h;                       // k1
  apf::setVector(rk1, v, 0, k1);
  apf::getVector(rk2, v, 0, point);
  m->setPoint(v, 0, point + k1/2);          // y1

  EoM->vd_upd_vel_field(v, drag_local);

  apf::getVector(vel_field, v, 0, temp);
  apf::Vector3 k2 = temp*h;                       // k4

  apf::getVector(rk2, v, 0, point);
  m->setPoint(v, 0, point);            // y_n
  apf::setVector(vel_field, v, 0, k2/h);
}

// Using Runge Kutta4 method, calculate the position and velocity of the 
// vertices. 
void field_calc::vd_calc_rk4(apf::Mesh2* m, apf::MeshEntity* v, 
                                            bool drag_local) {

  apf::Field* rk1 = m->findField("rk1");
  assert(rk1);
  apf::Field* rk2 = m->findField("rk2");
  assert(rk2);
  apf::Field* rk3 = m->findField("rk3");
  assert(rk3);
  apf::Field* rk4 = m->findField("rk4");
  assert(rk4);
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  double h = vdparam.dt;
  // y0
  apf::Vector3 temp(0,0,0);
  apf::Vector3 point(0,0,0);

  m->getPoint(v, 0, temp);
  apf::setVector(rk4, v, 0, temp);
  EoM->vd_upd_vel_field(v, drag_local);

  apf::getVector(vel_field, v, 0, temp);
  apf::Vector3 k1 = temp*h;                       // k1
  apf::setVector(rk1, v, 0, k1);
  apf::getVector(rk4, v, 0, point);
  m->setPoint(v, 0, point + k1/2);          // y1

  EoM->vd_upd_vel_field(v, drag_local);

  apf::getVector(vel_field, v, 0, temp);
  apf::Vector3 k2 = temp*h;                       // k2
  apf::setVector(rk2, v, 0, k2);
  apf::getVector(rk4, v, 0, point);
  m->setPoint(v, 0, point + k2/2);          // y2

  EoM->vd_upd_vel_field(v, drag_local);

  apf::getVector(vel_field, v, 0, temp);
  apf::Vector3 k3 = temp*h;                       // k3
  apf::setVector(rk3, v, 0, k3);
  apf::getVector(rk4, v, 0, point);
  m->setPoint(v, 0, point + k3/2);          // y3

  EoM->vd_upd_vel_field(v, drag_local);
  apf::getVector(vel_field, v, 0, temp);
  apf::Vector3 k4 = temp*h;                       // k4

  //temp = (k1+2*k2+2*k3+k4)/6;
  apf::getVector(rk1, v, 0, temp);
  apf::getVector(rk2, v, 0, point);
  temp = temp + point*2;
  apf::getVector(rk3, v, 0, point);
  temp = temp + point*2 + k4;
  temp =  temp/6;

  apf::getVector(rk4, v, 0, point);
  m->setPoint(v, 0, point);            // y_n
  apf::setVector(vel_field, v, 0, temp/h);
}

// Given the meshes and the map from the m2 entities to m1 entities, 
// check and transfer the fields.
void field_calc::vd_trans_fields(apf::Mesh2* m1, apf::Mesh2* m2,
                    std::vector<std::vector<apf::MeshEntity*> > & ents,
                    std::map<apf::MeshEntity*, apf::MeshEntity*> &e_map) {
  apf::Field* vel_field1 = m1->findField("velocity_field");
  apf::Field* vel_field2 = m2->findField("velocity_field");
  assert(vel_field1 and vel_field2);

  apf::Vector3 temp(0,0,0);
  if(t_int_flag == INTEG_TYPE::EULER) {
  }
  else if(t_int_flag == INTEG_TYPE::RK2) {
    apf::Field* rk11 = m1->findField("rk1");
    apf::Field* rk12 = m1->findField("rk2");
    apf::Field* rk21 = m2->findField("rk1");
    apf::Field* rk22 = m2->findField("rk2");
    assert(rk11 and rk12 and rk21 and rk22);
    for(int i = 0; i < ents.at(0).size(); i++) {
      apf::MeshEntity* v2 = ents.at(0).at(i);
      apf::MeshEntity* v1 = e_map[v2];
      assert(v1 != NULL);

      apf::getVector(rk11, v1, 0, temp);
      apf::setVector(rk21, v2, 0, temp);
      apf::getVector(rk12, v1, 0, temp);
      apf::setVector(rk22, v2, 0, temp);
    }

  }
  else {
    apf::Field* rk11 = m1->findField("rk1");
    apf::Field* rk12 = m1->findField("rk2");
    apf::Field* rk13 = m1->findField("rk3");
    apf::Field* rk14 = m1->findField("rk4");
    apf::Field* rk21 = m2->findField("rk1");
    apf::Field* rk22 = m2->findField("rk2");
    apf::Field* rk23 = m2->findField("rk3");
    apf::Field* rk24 = m2->findField("rk4");
    assert(rk11 and rk12 and rk13 and rk14 and 
           rk21 and rk22 and rk23 and rk24);
    for(int i = 0; i < ents.at(0).size(); i++) {
      apf::MeshEntity* v2 = ents.at(0).at(i);
      apf::MeshEntity* v1 = e_map[v2];
      assert(v1 != NULL);

      apf::getVector(rk11, v1, 0, temp);
      apf::setVector(rk21, v2, 0, temp);
      apf::getVector(rk12, v1, 0, temp);
      apf::setVector(rk22, v2, 0, temp);
      apf::getVector(rk13, v1, 0, temp);
      apf::setVector(rk23, v2, 0, temp);
      apf::getVector(rk14, v1, 0, temp);
      apf::setVector(rk24, v2, 0, temp);
    }
  }
  vd_trans_tri_fields(m1, m2, ents, e_map);
}

// Given the meshes and the map from the m2 entities to m1 entities, 
// check and transfer the fields.
void field_calc::vd_trans_tri_fields(apf::Mesh2* m1, apf::Mesh2* m2,
                    std::vector<std::vector<apf::MeshEntity*> > & ents,
                    std::map<apf::MeshEntity*, apf::MeshEntity*> &e_map) {

  apf::Field* d21 = m1->findField("d2");
  apf::Field* gam21 = m1->findField("gam2");
  apf::Field* d22 = m1->findField("d2");
  apf::Field* gam22 = m1->findField("gam2");
  assert(d21 and d22 and gam21 and gam22);
  for(int i = 0; i < ents.at(2).size(); i++) {
    apf::MeshEntity* t2 = ents.at(2).at(i);
    apf::MeshEntity* t1 = e_map[t2];

    double temp = apf::getScalar(d21, t1, 0);
    apf::setScalar(d22, t2, 0, temp);
    temp = apf::getScalar(gam21, t1, 0);
    apf::setScalar(gam22, t2, 0, temp);
  }

}

void field_calc::vd_att_fields(apf::Mesh2* m) {
  apf::Field* vel_field = vd_att_vv_field(m, "velocity_field");
  if(t_int_flag == INTEG_TYPE::EULER) {
  }
  else if(t_int_flag == INTEG_TYPE::RK2) {
    apf::Field* rk1 = vd_att_vv_field(m, "rk1");
    apf::Field* rk2 = vd_att_vv_field(m, "rk2");
  }
  else {
    apf::Field* rk1 = vd_att_vv_field(m, "rk1");
    apf::Field* rk2 = vd_att_vv_field(m, "rk2");
    apf::Field* rk3 = vd_att_vv_field(m, "rk3");
    apf::Field* rk4 = vd_att_vv_field(m, "rk4");
  }
  std::cout << "RK fields attached." << std::endl;

  vd_att_tri_fields(m);
}

void field_calc::vd_del_fields(apf::Mesh2* m) {
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field != NULL);
  apf::destroyField(vel_field);

  if(t_int_flag == INTEG_TYPE::EULER) {
  }
  else if(t_int_flag == INTEG_TYPE::RK2) {
    apf::Field* rk1 = m->findField( "rk1");
    apf::Field* rk2 = m->findField( "rk2");
    assert(rk1 != NULL and rk2 != NULL);
    apf::destroyField(rk1);
    apf::destroyField(rk2);
  }
  else {
    apf::Field* rk1 = m->findField( "rk1");
    apf::Field* rk2 = m->findField( "rk2");
    apf::Field* rk3 = m->findField( "rk3");
    apf::Field* rk4 = m->findField( "rk4");
    assert(rk1 != NULL and rk2 != NULL);
    assert(rk3 != NULL and rk4 != NULL);
    apf::destroyField(rk1);
    apf::destroyField(rk2);
    apf::destroyField(rk3);
    apf::destroyField(rk4);
  }
  std::cout << "RK fields attached." << std::endl;

  del_tri_fields(m);
}

void field_calc::vd_att_fields(apf::Mesh2* m, apf::MeshEntity* v) {
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);
  apf::setVector(vel_field, v, 0, apf::Vector3(0,0,0));
  if(t_int_flag == INTEG_TYPE::EULER) {
  }
  else if(t_int_flag == INTEG_TYPE::RK2) {
    apf::Field* rk1 = m->findField("rk1");
    apf::Field* rk2 = m->findField("rk2");
    assert((rk1 and rk2));
    apf::setVector(rk1, v, 0, apf::Vector3(0,0,0));
    apf::setVector(rk2, v, 0, apf::Vector3(0,0,0));
  }
  else {
    apf::Field* rk1 = m->findField("rk1");
    apf::Field* rk2 = m->findField("rk2");
    apf::Field* rk3 = m->findField("rk3");
    apf::Field* rk4 = m->findField("rk4");
    assert((rk1 and rk2 and rk3 and rk4));
    apf::setVector(rk1, v, 0, apf::Vector3(0,0,0));
    apf::setVector(rk2, v, 0, apf::Vector3(0,0,0));
    apf::setVector(rk3, v, 0, apf::Vector3(0,0,0));
    apf::setVector(rk4, v, 0, apf::Vector3(0,0,0));
  }
  std::cout << "RK fields attached to " << v << "." << std::endl;

  apf::Field* d2_field = m->findField("d2");
  apf::Field* gam2_field = m->findField("gam2");
  if(gam2_field) {
    assert(d2_field);
    std::vector<apf::MeshEntity*> edge(0);
    std::vector<apf::MeshEntity*> tri(0);

    vd_set_up(m, v, &edge);
    vd_set_up(m, &edge, &tri);
    for(int i = 0; i < tri.size(); i++)
      vd_att_tri_fields(m, tri.at(i));
  }
}

void field_calc::vd_att_tri_fields(apf::Mesh2* m, apf::MeshEntity* tri) {
  apf::Field* d2_field = m->findField("d2");
  apf::Field* gam2_field = m->findField("gam2");

  if(m->getModelType(m->toModel(tri)) == 2) {
    apf::setScalar(gam2_field, tri, 0, vdparam.vd_gam);
    apf::setScalar(d2_field, tri, 0, 1/vdparam.mob);
  }
  else {
    apf::setScalar(gam2_field, tri, 0, 0);
    apf::setScalar(d2_field, tri, 0, 0);
  }
}

void field_calc::vd_att_tri_fields(apf::Mesh2* m) {
  apf::Field* d2_field = vd_att_ts_field(m, "d2");
  apf::Field* gam2_field = vd_att_ts_field(m, "gam2");

  upd_d2(m);
  upd_gam2(m);
}

void field_calc::del_tri_fields(apf::Mesh2* m) {
  apf::Field* d2_field = m->findField("d2");
  apf::Field* gam2_field = m->findField("gam2");

  assert(d2_field != NULL and gam2_field != NULL);
  apf::destroyField(d2_field);
  apf::destroyField(gam2_field);
}

bool field_calc::vd_chk_fields(apf::Mesh2* m) {
  apf::Field* vel_field = m->findField("velocity_field");
  if(!vel_field)
    return false;
  if(t_int_flag == INTEG_TYPE::EULER) {
  }
  else if(t_int_flag == INTEG_TYPE::RK2) {
    apf::Field* rk1 = m->findField("rk1");
    apf::Field* rk2 = m->findField("rk2");
    if(!rk1 or !rk2)
      return false;
  }
  else {
    apf::Field* rk1 = m->findField("rk1");
    apf::Field* rk2 = m->findField("rk2");
    apf::Field* rk3 = m->findField("rk3");
    apf::Field* rk4 = m->findField("rk4");
    if(!rk1 or !rk2 or !rk3 or !rk4)
      return false;
  }
  return true;
}

// Using Runge Kutta4 method, calculate the velocity of a set of vertices. Used
// before calc_roc to estimate the rate of change less affected by oscillations
// in the solution. Sets the vertices to their old positions, unlike the other
// rk4 implemented above.
void field_calc::vd_calc_rk4(apf::Mesh2* m, std::vector<apf::MeshEntity*> &verts) {

  apf::Field* rk1 = m->findField("rk1");
  if(!rk1)
    rk1 = vd_att_vv_field(m, "rk1");
  assert(rk1);
  apf::Field* rk2 = m->findField("rk2");
  if(!rk2)
    rk2 = vd_att_vv_field(m, "rk2");
  apf::Field* rk3 = m->findField("rk3");
  if(!rk3)
    rk3 = vd_att_vv_field(m, "rk3");
  apf::Field* rk4 = m->findField("rk4");
  if(!rk4)
    rk4 = vd_att_vv_field(m, "rk4");
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  double t_sub = find_min_t2(m, &verts);
  if(t_sub > std::numeric_limits<double>::min()) {
    vdparam.adj_dt(t_sub/8);
  }
  std::cout << "t_min " << t_sub << " " << vdparam.dt << std::endl;

  double h = vdparam.dt;
  // y0
  apf::Vector3 temp(0,0,0);
  apf::Vector3 point(0,0,0);

  for(int i = 0; i < verts.size(); i++) {
    m->getPoint(verts.at(i), 0, temp);
    apf::setVector(rk4, verts.at(i), 0, temp);
    //vd_upd_vel_field(m, verts.at(i), false);
  }
  EoM->calc_vel(&verts,false);

  for(int i = 0; i < verts.size(); i++) {
    apf::getVector(vel_field, verts.at(i), 0, temp);
    apf::Vector3 k1 = temp*h;                       // k1
    apf::setVector(rk1, verts.at(i), 0, k1);
    apf::getVector(rk4, verts.at(i), 0, point);
    m->setPoint(verts.at(i), 0, point + k1/2);          // y1
  }

  EoM->calc_vel(&verts,false);
//  for(int i = 0; i < verts.size(); i++) {
//    vd_upd_vel_field(m, verts.at(i), false);
//  }

  for(int i = 0; i < verts.size(); i++) {
    apf::getVector(vel_field, verts.at(i), 0, temp);
    apf::Vector3 k2 = temp*h;                       // k2
    apf::setVector(rk2, verts.at(i), 0, k2);
    apf::getVector(rk4, verts.at(i), 0, point);
    m->setPoint(verts.at(i), 0, point + k2/2);          // y2
  }

  EoM->calc_vel(&verts,false);

  //for(int i = 0; i < verts.size(); i++) {
  //  vd_upd_vel_field(m, verts.at(i), false);
  //}

  for(int i = 0; i < verts.size(); i++) {
    apf::getVector(vel_field, verts.at(i), 0, temp);
    apf::Vector3 k3 = temp*h;                       // k3
    apf::setVector(rk3, verts.at(i), 0, k3);
    apf::getVector(rk4, verts.at(i), 0, point);
    m->setPoint(verts.at(i), 0, point + k3/2);          // y3
  }

  //for(int i = 0; i < verts.size(); i++) {
  //  vd_upd_vel_field(m, verts.at(i), false);
  //}
  EoM->calc_vel(&verts,false);

  for(int i = 0; i < verts.size(); i++) {
    apf::getVector(vel_field, verts.at(i), 0, temp);
    apf::Vector3 k4 = temp*h;                       // k4

    //temp = (k1+2*k2+2*k3+k4)/6;
    apf::getVector(rk1, verts.at(i), 0, temp);
    std::cout << "v " << verts.at(i)
              << " rk1 " << temp;

    apf::getVector(rk2, verts.at(i), 0, point);
    std::cout << " rk2 " << point;
    temp = temp + point*2;
    apf::getVector(rk3, verts.at(i), 0, point);
    std::cout << " rk3 " << point << " rk4 " << k4 << std::endl;
    temp = temp + point*2 + k4;
    temp =  temp/6;

    apf::getVector(rk4, verts.at(i), 0, point);
    m->setPoint(verts.at(i), 0, point);            // y_n
    std::cout << "pos " << point << " vel " << temp/h << std::endl;

    apf::setVector(vel_field, verts.at(i), 0, temp/h);
  }

  //k1 = h*v0(y0)
  //y1 = y_0 + k1/2
  //y1 = y_0 + v0*h/2
  //k2 = h*v1(y1)
  //y2 = y_0 + k2/2
  //y2 = y_0 + v1*h/2
  //k3 = h*v2(y2)
  //y3 = y_0 + k3
  //y3 = y_0 + v2*h
  //k4 = h*v3(y3)
  //y_n = y0 + (k1+2*k2+2*k3+k4)/6
}

double field_calc::find_min_t2(apf::Mesh2* m, double t_set) {
  return vd_find_min_t2(m, calc_ext, t_set);
}

double field_calc::find_min_t_tet(apf::Mesh2* m) {
  return vd_find_min_t_tet(m, calc_ext);
}

double field_calc::find_min_t2(apf::Mesh2* m, std::vector<apf::MeshEntity*>* vert, double t_set) {
  return vd_find_min_t2(m, vert, calc_ext, t_set);
}

double field_calc::find_min_t2(apf::Mesh2* m, apf::MeshEntity* vert, double t_set) {
  return vd_find_min_t2(m, vert, calc_ext, t_set);
}

double field_calc::vd_get_max_vel_norm_int(apf::Mesh2* m) {
  vd_entlist e_list(m, c_base);

  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  apf::Vector3 dn(0,0,0);
  double norm_max = 0.;
//  for(int dim = 0; dim < 3; dim++) {
  int dim = 0;
  int sz = c_base->get_sz(dim);
  for(int i = 0; i < sz; i++) {
    if(!c_base->is_free(dim,i) and !c_base->get_cell_ext(dim, i)) {

      for(int j = 0; j < e_list.e.at(dim).at(i).at(0).size(); j++) {
        apf::MeshEntity* v_curr = e_list.e.at(dim).at(i).at(0).at(j);
        apf::getVector(vel_field, v_curr, 0, dn);
        double norm_curr = dn.getLength();
        if(norm_curr > norm_max)
          norm_max = norm_curr;
      }
    }
  }
  return norm_max;
}

apf::MeshEntity* field_calc::create_v(apf::Mesh2* m, apf::ModelEntity* mdl) {

  apf::MeshEntity* v_new = m->createVert(mdl);
  apf::Vector3 temp(0,0,0);
  apf::Field* vel_field = m->findField("velocity_field");
  if(!vel_field) {
    vd_att_fields(m);
    vel_field = m->findField("velocity_field");
  }

  assert(vel_field);
  apf::setVector(vel_field, v_new, 0, temp);

  if(t_int_flag == INTEG_TYPE::EULER) {
  }
  else if(t_int_flag == INTEG_TYPE::RK2) {
    apf::Field* rk1 = m->findField("rk1");
    assert(rk1);
    apf::Field* rk2 = m->findField("rk2");
    assert(rk2);
    apf::setVector(rk1, v_new, 0, temp);
    apf::setVector(rk2, v_new, 0, temp);
  }
  else {
    assert(t_int_flag == INTEG_TYPE::RK4);
    apf::Field* rk1 = m->findField("rk1");
    assert(rk1);
    apf::Field* rk2 = m->findField("rk2");
    assert(rk2);
    apf::Field* rk3 = m->findField("rk3");
    assert(rk3);
    apf::Field* rk4 = m->findField("rk4");
    assert(rk4);
    apf::setVector(rk1, v_new, 0, temp);
    apf::setVector(rk2, v_new, 0, temp);
    apf::setVector(rk3, v_new, 0, temp);
    apf::setVector(rk4, v_new, 0, temp);
  }
  return v_new;
}

////////////////////////////
// Energy and drag:
////////////////////////////
double field_calc::upd_gam2(apf::Mesh2* m) {
  apf::Field* gam2_field = m->findField("gam2");
  assert(gam2_field);

  int c_type, c_tag;
  apf::ModelEntity* mdl;
  apf::MeshEntity* tri;
  apf::MeshIterator* it_e = m->begin(2);
  while (tri = m->iterate(it_e)) {
    mdl = m->toModel(tri);
    c_type = m->getModelType(mdl);
    c_tag = m->getModelTag(mdl);
    if(c_type == 2) {
      if(!b_en.c2s_zeros[c_tag])
        apf::setScalar(gam2_field, tri, 0, vdparam.vd_gam);
      else
        apf::setScalar(gam2_field, tri, 0, 0);
    }
    else
      apf::setScalar(gam2_field, tri, 0, 0);

  }
  m->end(it_e);
}

void field_calc::set_vdparam(vd_param vdpar_in, apf::Mesh2* m) {
  vdparam = vdpar_in;
  upd_d2(m);
  upd_gam2(m);
}

double field_calc::upd_d2(apf::Mesh2* m) {
  apf::Field* d2_field = m->findField("d2");
  assert(d2_field);

  d2_glob = 1/vdparam.mob;

  int c_type, c_tag;
  apf::ModelEntity* mdl;
  apf::MeshEntity* tri;
  apf::MeshIterator* it_e = m->begin(2);
  while (tri = m->iterate(it_e)) {
    mdl = m->toModel(tri);
    c_type = m->getModelType(mdl);
    c_tag = m->getModelTag(mdl);
    if(c_type == 2) {
      apf::setScalar(d2_field, tri, 0, d2_glob);
    }
    else
      apf::setScalar(d2_field, tri, 0, 0);
  }
  m->end(it_e);
}

double field_calc::upd_gam2(apf::Mesh2* m, std::vector<apf::MeshEntity*>& tri, double val) {
  apf::Field* gam2_field = m->findField("gam2");
  assert(gam2_field);

  for(int i = 0; i < tri.size(); i++) {
    apf::setScalar(gam2_field, tri.at(i), 0, val);
  }
}

double field_calc::upd_d2(apf::Mesh2* m, std::vector<apf::MeshEntity*>& tri, double val) {
  apf::Field* d2_field = m->findField("d2");
  assert(d2_field);

  for(int i = 0; i < tri.size(); i++) {
    apf::setScalar(d2_field, tri.at(i), 0, val);
  }
}

////////////////////////////////////////////////////

double field_calc::d2(apf::Mesh2* m, apf::MeshEntity* tri) {
  apf::Field* d2_field = m->findField("d2");
  return apf::getScalar(d2_field, tri, 0);
}

double field_calc::gam2(apf::Mesh2* m, apf::MeshEntity* tri) {
  apf::Field* gam2_field = m->findField("gam2");
  return apf::getScalar(gam2_field, tri, 0);
}

double field_calc::en_tri(apf::Mesh2* m, apf::MeshEntity* tri) {
  apf::Field* gam2_field = m->findField("gam2");
  double area = vd_meas_ent(m, tri);
  return apf::getScalar(gam2_field, tri, 0)*area;
}


// Refresh the mesh and vd_elist lists of the EoM.
void field_calc::refresh_mesh() {
  EoM->refresh_elist();
}

field_calc::field_calc() :
  EoM(NULL), c_base(NULL),
  calc_ext(false), calc_corner(false), cb_flag(false),
  proj_flag((int)PROJ_TYPE::FIXED), e_sh(NULL),
  v_flag(VEL_TYPE::MASON), t_int_flag(INTEG_TYPE::EULER),
  e_sh_flag(false), drag_glob(-1), drag_rat(1000),
  d2_glob(-1),
  b_en(),
  vdparam() {
    std::cout << "Velocity calculation type: Mason." << std::endl;
    EoM = new vd_eqn_mason();
}

void field_calc::set_vel_type(VEL_TYPE DT) {
  if (DT == VEL_TYPE::LAZAR) {
    std::cout << "Velocity calculation type: Lazar." << std::endl;
    EoM = new vd_eqn_lazar();
  }
  else if (DT == VEL_TYPE::KUPRAT) {
    std::cout << "Velocity calculation type: Kuprat." << std::endl;
    EoM = new vd_eqn_kuprat();
  }
  else if (DT == VEL_TYPE::KUPRAT_NBC) {
    std::cout << "Velocity calculation type: Kuprat_NBC." << std::endl;
    EoM = new vd_eqn_kuprat_NBC();
  }
  else if (DT == VEL_TYPE::MASON) {
    std::cout << "Velocity calculation type: Mason." << std::endl;
    EoM = new vd_eqn_mason();
  }
  else if (DT == VEL_TYPE::MASON_MIR) {
    std::cout << "Velocity calculation type: Mason_mirror." << std::endl;
    EoM = new vd_eqn_mason_mir();
  }
  else if (DT == VEL_TYPE::MASON_NBC) {
    std::cout << "Velocity calculation type: Mason_NBC." << std::endl;
    EoM = new vd_eqn_mason_NBC();
  }
  else {
    std::cout << "Velocity calculation type: Mason." << std::endl;
    EoM = new vd_eqn_mason();
  }
}

field_calc::field_calc(VEL_TYPE DT) :
  EoM(NULL),
  calc_ext(false), calc_corner(false), cb_flag(false),
  proj_flag((int)PROJ_TYPE::FIXED), e_sh(NULL),
  v_flag(DT), t_int_flag(INTEG_TYPE::EULER),
  e_sh_flag(false), drag_glob(-1), drag_rat(1000),
  d2_glob(-1),
  b_en(),
  vdparam() {

  set_vel_type(DT);

  calc_ext = false;
  calc_corner = false;
  cb_flag = false;
  proj_flag == (int)PROJ_TYPE::FIXED;
  //c3_out = -1;
  e_sh = NULL;
  e_sh_flag = false;
}

field_calc::~field_calc() {
  delete EoM;
}

void field_calc::reload_cb(cell_base* c_in) {
  //m = NULL;
  c_base = c_in;
  cb_flag = true;
}
void field_calc::off_cb() {
  //m = NULL;
  c_base = NULL;
  cb_flag = false;
}

void field_calc::set_calc_ext(bool calc_ext_in) {
//int c3_out_in) {
  calc_ext = calc_ext_in;
  //c3_out = c3_out_in;
}

//int field_calc::get_c3_out() {
//  return c3_out;
//}

bool field_calc::get_ext() {
  return calc_ext;
}

bool field_calc::get_corner() {
  return calc_corner;
}

int field_calc::get_proj() {
  return proj_flag;
}

ext_shell field_calc::copy_e_sh() {
  return *e_sh;
}

ext_shell* field_calc::get_e_sh() {
  return e_sh;
}

void field_calc::set_e_sh(ext_shell* e_sh_in) {
  *e_sh = *e_sh_in;
}

bool field_calc::chk_skip(apf::Mesh2* m, apf::MeshEntity* vert) {
  assert(cb_flag);
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  if(em_type == 3)
    return true;

  else if (!get_ext() and chk_vert_special(m, vert))
    return true;
}

bool field_calc::chk_vert_special(apf::Mesh2* m, apf::MeshEntity* vert) {
  assert(cb_flag);
  //if(vd_chk_vert_dom(m, vert))
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  int em_tag = m->getModelTag(em);

  if(em_type < 3 and c_base->get_cell_ext_gmi(em_type, em_tag))
    return true;
}
bool field_calc::chk_0c_cor_gmi(int tag_0cell) {
  if(e_sh->chk_shell(0, tag_0cell - 1)) {
    shell sh = e_sh->get_shell(0, tag_0cell - 1);
    if(sh.dim == 0)
      return true;
  }
  return false;
}

void field_calc::set_vec_sp_calc(PROJ_TYPE PROJ, ext_shell* e_sh_in) {
  e_sh = NULL;
  e_sh_flag = false;

  proj_flag = (int)PROJ;
  assert(proj_flag < (int)PROJ_TYPE::END);

  if(proj_flag > (int)PROJ_TYPE::FIXED)
    calc_ext = true;

  if(proj_flag > (int)PROJ_TYPE::FIXED and 
     proj_flag < (int)PROJ_TYPE::END) {
    calc_corner = true;
  }

  if(proj_flag == (int)PROJ_TYPE::EXT_SHELL) {
    assert(e_sh_in != NULL);
    e_sh = e_sh_in;
    e_sh_flag = true;
  }

}

void field_calc::set_integ_type(INTEG_TYPE IT) {
  t_int_flag = IT;
}

// Set the drag coefficient for drag matrix regularization before inversion.
void field_calc::set_drag_glob(double drag_in) {
  drag_glob = drag_in;
}

double field_calc::get_drag_glob() {
  return drag_glob;
}

// Set the triangle drag coefficient.
void field_calc::set_d2_glob(double drag_in) {
  d2_glob = drag_in;
}

double field_calc::get_d2_glob() {
  return d2_glob;
}

double field_calc::get_drag_rat() {
  return drag_rat;
}

void field_calc::set_drag_rat(double rat_in) {
  drag_rat = rat_in;
}

void field_calc::dummy_func_stop() {
}

apf::Vector3 field_calc::get_vec_special(apf::Mesh2* m, apf::MeshEntity* vert,
                                         apf::Vector3 v_in) {

  apf::Vector3 dir(0,0,0);
  if(proj_flag == (int)PROJ_TYPE::FIXED)
    return dir;

  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  int em_tag = m->getModelTag(em);

  if(proj_flag < (int)PROJ_TYPE::EXT_SHELL) {

    //assert(c_base->get_cell_ext_gmi(em_type, em_tag));

    apf::Up up;
    apf::Downward down;
    std::vector<apf::MeshEntity*> es_edge(0);
    std::vector<apf::MeshEntity*> es_surf(0);
    std::vector<apf::MeshEntity*> es_elem(0);

    m->getUp(vert, up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);
    vd_set_up(m, &es_surf, &es_elem);

    //if(c3_out > -1)
    //  vd_rem_cell(m, &es_elem, c3_out, 3);
    apf::Vector3 vec_vp(0,0,0);
    vec_vp = get_volm_pres(m, vert, &es_elem, v_in);
    std::cout << "Volm projected " << vec_vp << std::endl;

    if(proj_flag == (int)PROJ_TYPE::LOCAL_TRI or 
       proj_flag == (int)PROJ_TYPE::LOCAL_SP) {
      std::cout << "Volm preserved " << std::endl;
      return vec_vp;
    }

    if((em_type == 0 and c_base->get_0c_corner_gmi(em_tag)) or
       (em_type == 1 and c_base->get_1c_corner_gmi(em_tag))) {

      // Project the velocity along all direction corner edge directions. 

      // Find one corner edge. Use its direction. 
      apf::Vector3 p1(0,0,0);
      dir = apf::Vector3(0,0,0);

      int sz = 0;
      for(int i = 0; i < es_edge.size(); i++) {

        apf::ModelEntity* em1 = m->toModel(es_edge.at(i));
        int em1_type = m->getModelType(em1);
        int em1_tag = m->getModelTag(em1);
        if(em1_type == 1 and c_base->get_1c_corner_gmi(em1_tag)) {
          m->getDownward(es_edge.at(i), 0, down);

          m->getPoint(down[0], 0, p1);
          m->getPoint(down[1], 0, dir);
          dir = dir - p1;

          if(dir.getLength() > std::numeric_limits<double>::epsilon()) {
            dir = norm_0(dir - p1);
            if(vec_vp*dir < 0)
              dir = dir*(-1);
            vec_vp = dir*(vec_vp*dir);
            //dummy_func_stop();
          }
        }
      }
      std::cout << "Side projected " << vec_vp << std::endl;
    }

    assert(proj_flag == (int)PROJ_TYPE::FLAT);
    return vec_vp;
  }
  else {
    assert(proj_flag == (int)PROJ_TYPE::EXT_SHELL);
    // Unlike the others that just fixes the vector, this projects the position
    // on the shell and projects the vector on the shell, at the same time.

    // The position should already be on the shell, but this step guarantees it.
    if(e_sh->chk_shell(em_type, em_tag-1)) {
      // If the shell is of the same dimensions, the stratum vertex will not be
      // moving. TODO this is to circumnavigate the "sliding" of vertices due 
      // to lack of drag for the moment. Playing with global drag is more 
      // finnicky.
      apf::Vector3 pos(0,0,0);
      shell sh = e_sh->get_shell(em_type, em_tag-1);
      if(sh.dim == em_type)
        return pos;
      m->getPoint(vert, 0, pos);
      //std::cout << "Pos_old " << pos
      //          << " pos ";
      pos = e_sh->find_int_pos(em_type, em_tag-1, pos);
      //std::cout << pos << std::endl;
      //m->setPoint(vert, 0, e_sh->find_int_pos(em_type, em_tag-1, pos));
      m->setPoint(vert, 0, pos);

      return e_sh->find_para_dir(em_type, em_tag-1, v_in);
    }
    return v_in;
  }
}

// Given a vertex, if the vertex is outside the bounds of the shell, push it 
// interior by distance to the bounding shell.
apf::Vector3 field_calc::fix_vert_bound(apf::Mesh2* m, apf::MeshEntity* vert, apf::Vector3 pos) {
  shell sh_curr;
  shell sh_b;
  ent_conn sh_d;
  int c_type = m->getModelType(m->toModel(vert));
  int c_tag = m->getModelTag(m->toModel(vert));

  apf::Vector3 temp1(0,0,0);
  apf::Vector3 temp2(0,0,0);
  apf::Vector3 pos_sh(0,0,0);
  apf::Vector3 pos_e(0,0,0);
  apf::Vector3 dir_e(0,0,0);

  // If member of a shell, check if the vertex passes through a bounding shell
  // by the displacement. For each bounding shell, check if the component of the
  // normal displacement has a larger magnitude than the normal distance of the
  // position to the shell. Scale the displacement such that the highest ratio 
  // of the normal displacement to the shell to the distance to the shell is 
  // equal to scale.
  if(e_sh->chk_shell(c_type, c_tag-1)) {
    sh_curr = e_sh->get_shell(c_type, c_tag-1);
    e_sh->sh_base.get_conn_gmi(sh_curr.dim, sh_curr.id, &sh_d);
    pos_sh = e_sh->get_shell_pos(sh_curr);

    if(sh_curr.dim == 2) {
      for(int i = 0; i < sh_d.conn.size(); i++) {
        sh_b = shell(sh_curr.dim-1, sh_d.conn.at(i));
        pos_e = e_sh->get_shell_pos(sh_b);
        dir_e = e_sh->get_shell_dir(sh_b);
        temp1 = (pos - pos_e);
        temp1 = temp1 - dir_e*(temp1*dir_e);
        temp2 = (pos_sh - pos_e);
        temp2 = temp2 - dir_e*(temp2*dir_e);
        if( temp1*temp2 < - std::numeric_limits<double>::min()) {
          pos = pos - temp1;
        }
      }
    }
    else {
      assert(sh_curr.dim == 1);
      for(int i = 0; i < sh_d.conn.size(); i++) {
        sh_b = shell(sh_curr.dim-1, sh_d.conn.at(i));
        pos_e = e_sh->get_shell_pos(sh_b);
        temp1 = (pos - pos_e);
        temp2 = (pos_sh - pos_e);
        if( temp1*temp2 < - std::numeric_limits<double>::min()) {
          pos = pos - temp1;
        }
      }
    }
  }
  // Else, check all face shells:
  else {
    for(int i = 0; i < e_sh->sh_base.get_sz(2); i++) {
      if(!e_sh->sh_base.is_free(2, i)) {
        for(int i = 0; i < sh_d.conn.size(); i++) {
          sh_b = shell(sh_curr.dim-1, sh_d.conn.at(i));
          pos_e = e_sh->get_shell_pos(sh_b);
          dir_e = e_sh->get_shell_dir(sh_b);
          temp1 = (pos - pos_e);
          temp1 = dir_e*(temp1*dir_e);
          temp2 = (pos_sh - pos_e);
          temp2 = dir_e*(temp2*dir_e);
          if( temp1*temp2 < - std::numeric_limits<double>::min()) {
            pos = pos - temp1;
          }
        }
      }
    }
  }
  return pos;
}

// Given a shell, a position and a displacement, fix the displacement vector
// such that the new position remains within bounds of the shell.
apf::Vector3 field_calc::fix_disp_bound(apf::Mesh2* m, apf::MeshEntity* vert,
                       apf::Vector3 pos, apf::Vector3 disp, double scale) {

  if(disp.getLength() < std::numeric_limits<double>::min())
    return disp;

  shell sh_curr;
  shell sh_b;
  ent_conn sh_d;
  int c_type = m->getModelType(m->toModel(vert));
  int c_tag = m->getModelTag(m->toModel(vert));

  apf::Vector3 disp_temp(0,0,0);
  apf::Vector3 disp_temp_dir(0,0,0);
  double disp_mag = disp.getLength();

  // If member of a shell, check if the vertex passes through a bounding shell
  // by the displacement. For each bounding shell, check if the component of the
  // normal displacement has a larger magnitude than the normal distance of the
  // position to the shell. Scale the displacement such that the highest ratio 
  // of the normal displacement to the shell to the distance to the shell is 
  // equal to scale.
  double sc_high = 0;
  if(e_sh->chk_shell(c_type, c_tag-1)) {
    sh_curr = e_sh->get_shell(c_type, c_tag-1);
    e_sh->sh_base.get_conn_gmi(sh_curr.dim, sh_curr.id, &sh_d);

    for(int i = 0; i < sh_d.conn.size(); i++) {
      sh_b = shell(sh_curr.dim-1, sh_d.conn.at(i));
      disp_temp = e_sh->find_norm_disp(sh_b, pos);
      double mag_curr = disp_temp.getLength();
      disp_temp_dir = norm_0(disp_temp);
      double disp_curr = disp_temp_dir*disp;
      if(disp_curr > std::numeric_limits<double>::min()) {
        // Remove the component of the displacement in the normal direction if the
        // distance is close to 0.
        if(mag_curr < std::numeric_limits<double>::epsilon()) {
          if(sh_b.dim == 0) {
            disp = apf::Vector3(0,0,0);
          }
          else {
            assert(sh_b.dim == 1);
            disp_temp = e_sh->get_shell_dir(sh_b);
            disp_temp_dir = norm_0(disp_temp);
            disp = disp_temp_dir*(disp_temp_dir*disp);
          }
        }
        else if(mag_curr < disp_curr) {
          disp = disp - disp_temp_dir*(disp_curr-mag_curr*scale);
        }

      }
      // If displacement component along the normal is zero or is away from the 
      // shell, no change.
      //else {
      //}
    }
  }
  // Else, check all face shells:
  else {
    for(int i = 0; i < e_sh->sh_base.get_sz(2); i++) {
      if(!e_sh->sh_base.is_free(2, i)) {
        sh_b = shell(2, i+1);
        disp_temp = e_sh->find_norm_disp(sh_b, pos);
        double mag_curr = disp_temp.getLength();
        disp_temp_dir = norm_0(disp_temp);
        double disp_curr = disp_temp_dir*disp;
        if(disp_curr > std::numeric_limits<double>::min()) {
          // Remove the component of the displacement in the normal direction if 
          // the distance is close to 0.
          if(mag_curr < std::numeric_limits<double>::epsilon()) {
            assert(sh_b.dim == 2);
            disp_temp_dir = e_sh->get_shell_dir(sh_b);
            disp = disp - disp_temp_dir*(disp_temp_dir*disp);
          }
          else if(mag_curr < disp_curr) {
            disp = disp - disp_temp_dir*(disp_curr-mag_curr*scale);
          }

        }
        // If displacement component along the normal is zero or is away from the 
        // shell, no change.
        //else {
        //}
      }
    }
  }

  return disp;

}

void field_calc::fix_vel_special(apf::Mesh2* m, apf::MeshEntity* vert) {
  apf::Vector3 temp(0,0,0);

  apf::Field* vel_field = m->findField("velocity_field");
  apf::getVector(vel_field, vert, 0, temp);
  if (chk_vert_special(m, vert)) {
    apf::setVector(vel_field, vert, 0, get_vec_special(m, vert, temp));
  }
}

void field_calc::fix_vel_special_merg(apf::Mesh2* m, 
                  std::vector<apf::MeshEntity*> &vert, 
                  std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map) {

  std::map<apf::MeshEntity*, bool> merg{};
  apf::MeshEntity* v_merg = NULL;
  for(int i = 0; i < vert.size(); i++) {
    if(v_map[vert.at(i)] != NULL) {
      merg[vert.at(i)] = true;
      if(v_merg == NULL)
        v_merg = v_map[vert.at(i)];
    }
  }
  assert(v_merg != NULL);

  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find(vert.begin(), vert.end(), v_merg);
  int i_m = std::distance(vert.begin(), it);

  apf::Field* vel_field = m->findField("velocity_field");
  apf::Vector3 temp(0,0,0);

  for(int i = 0; i < vert.size(); i++) {
    if(!merg[vert.at(i)]) {
      apf::getVector(vel_field, vert.at(i), 0, temp);
      if (chk_vert_special(m, vert.at(i))) {
        apf::setVector(vel_field, vert.at(i), 0, get_vec_special(m, vert.at(i), 
                                                                      temp));
      }
    }
  }
  apf::getVector(vel_field, vert.at(i_m), 0, temp);
  if (chk_vert_special(m, vert.at(i_m))) {
    temp = get_vec_special(m, vert.at(i_m), temp);
  }
  apf::setVector(vel_field, vert.at(i_m), 0, temp);
  for(int i = 0; i < vert.size(); i++) {
    if(merg[vert.at(i)]) {
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
}

void field_calc::corr_pos(apf::Mesh2* m) {
  //dummy_func_stop();
  assert(proj_flag == (int)PROJ_TYPE::EXT_SHELL);

  apf::MeshEntity* e_v;
  apf::ModelEntity* mdl;
  int c_type;
  int c_tag;
  apf::Vector3 pos;

  apf::MeshIterator* it_e = m->begin(0);
  while (e_v = m->iterate(it_e)) {
    mdl = m->toModel(e_v);
    c_type = m->getModelType(mdl);
    c_tag = m->getModelTag(mdl);
    if(e_sh->chk_shell(c_type, c_tag-1)) {
      assert(c_base->get_cell_ext_gmi(c_type, c_tag));
      //shell sh = e_sh->get_shell(c_type, c_tag-1);
      //std::cout << c_type << "c" << c_tag << " " 
      //          << sh.dim << "sh" << sh.id
      //          << e_sh->get_shell_pos(sh) << " "
      //          << std::endl;
      //if(sh.dim > 0)
      //  std::cout << e_sh->get_shell_dir(sh) << std::endl;

      m->getPoint(e_v, 0, pos);
      //std::cout << "Pos_old " << pos
      //          << " pos ";
      pos = e_sh->find_int_pos(c_type, c_tag-1, pos);
      //std::cout << pos << std::endl;
      m->setPoint(e_v, 0, pos);
    }
  }
  m->end(it_e);
}

void field_calc::corr_pos(apf::Mesh2* m, apf::MeshEntity* e_v) {
  //dummy_func_stop();
  assert(proj_flag == (int)PROJ_TYPE::EXT_SHELL);

  apf::ModelEntity* mdl = m->toModel(e_v);
  apf::Vector3 pos;

  int c_type = m->getModelType(mdl);
  int c_tag = m->getModelTag(mdl);

  if(e_sh->chk_shell(c_type, c_tag-1)) {
    assert(c_base->get_cell_ext_gmi(c_type, c_tag));

    m->getPoint(e_v, 0, pos);
    //std::cout << "Pos_old " << pos
    //          << " pos ";
    pos = e_sh->find_int_pos(c_type, c_tag-1, pos);
    //std::cout << pos << std::endl;
    m->setPoint(e_v, 0, pos);
  }
  //shell sh = e_sh->get_shell(c_type, c_tag-1);
  //std::cout << c_type << "c" << c_tag << " " 
  //          << sh.dim << "sh" << sh.id
  //          << e_sh->get_shell_pos(sh) << " "
  //          << std::endl;
  //if(sh.dim > 0)
  //  std::cout << e_sh->get_shell_dir(sh) << std::endl;

}

void field_calc::refresh(apf::Mesh2* m_in, cell_base* c_base_in, 
                                            vd_entlist* e_list_in) {
  EoM->refresh(m_in, c_base_in, this, e_list_in);
}

field_calc::field_calc(const field_calc& that) {
  delete EoM;
  set_vel_type(that.v_flag);

  calc_ext = that.calc_ext;
  calc_corner = that.calc_corner;

  proj_flag = that.proj_flag;

  cb_flag = that.cb_flag;
  c_base = that.c_base;

  e_sh = that.e_sh;
  e_sh_flag = that.e_sh_flag;
  drag_glob = that.drag_glob;
  drag_rat = that.drag_rat;

  v_flag = that.v_flag; 
  t_int_flag = that.t_int_flag;
}

field_calc& field_calc::operator=(const field_calc& that) {
  delete EoM;
  set_vel_type(that.v_flag);

  b_en = that.b_en;

  calc_ext = that.calc_ext;
  calc_corner = that.calc_corner;

  proj_flag = that.proj_flag;

  cb_flag = that.cb_flag;
  c_base = that.c_base;

  e_sh = that.e_sh;
  e_sh_flag = that.e_sh_flag;
  drag_glob = that.drag_glob;
  drag_rat = that.drag_rat;

  v_flag = that.v_flag; 
  t_int_flag = that.t_int_flag;

  vdparam = that.vdparam;

  return *this;
}

