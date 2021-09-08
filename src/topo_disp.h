#ifndef TOPO_DISP_H
#define TOPO_DISP_H

#include <algorithm>    

#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>

#include "topo_rand.h"
#include "topo_extinfo.h"
#include "topo_geom.h"

#include "topo_manip.h"

class vd_eqn_of_motion;
class vd_eqn_mason;
class vd_eqn_lazar;

// Containers for functions used in updating the mesh. Different velocity
// calculators could be specified using this container.
// TODO A similar object can be used in disc extension calculation to 
// update the displacement using a different approach.
// This can be incorporated into initialization parameters. For now they are 
// part of the program.
class vd_param {
  public:
    double mob;
    double dt;
    double vd_gam;
    double mult;
    double v_mult;

    vd_param() : mob(1.0), dt(0.0005), vd_gam(1.0) {
      v_mult = mob*vd_gam; 
      mult = dt*v_mult; 
    }
    // Adjust the parameters for mobility, boundary energy and time step.
    void adj_dt(double dt_new) { 
      dt = dt_new; 
      v_mult = mob*vd_gam; 
      mult = dt*v_mult; 
      printf("dt adjusted to %f.\n",dt);
    }
    void adj_mob(double mob_new) { 
      mob = mob_new; 
      v_mult = mob*vd_gam; 
      mult = dt*v_mult; 
      printf("Mobility adjusted to %f.\n", mob);
    }
    void adj_vd_gam(double vd_gam_new) { 
      vd_gam = vd_gam_new; 
      v_mult = mob*vd_gam; 
      mult = dt*v_mult; 
      printf("Boundary energy adjusted to %f.\n", vd_gam);
    }

    // Copy constructor
    vd_param(const vd_param& that) :
      mob(that.mob), dt(that.dt), vd_gam(that.vd_gam), v_mult(that.v_mult),
      mult(that.mult) {
    }
    // Copy
    vd_param& operator=(const vd_param& that) {
      mob = that.mob;
      dt = that.dt;
      vd_gam = that.vd_gam;
      v_mult = that.v_mult;
      mult = that.mult;
      return *this;
    }

};


class field_calc;

enum class VEL_TYPE {
  LAZAR,
  MASON,
  MASON_MIR,
  MASON_NBC,
  KUPRAT,
  KUPRAT_NBC,
  END
};

enum class INTEG_TYPE {
  EULER,
  RK2,
  RK4,
  END
};

enum class EXT_GAM_TYPE {
  CONSTANT,
  ZERO,
  END
};

enum class PROJ_TYPE {
  FIXED,
  LOCAL_TRI,
  LOCAL_SP,
  FLAT,
  EXT_SHELL,
  END
};

VEL_TYPE conv_str2vel(const std::string& str);
INTEG_TYPE conv_str2integ(const std::string& str);
EXT_GAM_TYPE conv_str2extgam(const std::string& str); 
PROJ_TYPE conv_str2proj(const std::string& str);

// ext_shell used to preserve volume by a naive approach.
class shell {
  private:

  public:
    int dim;
    int id;

    shell();
    shell(int dim_in, int id_in);
    ~shell();

    // Copy constructor
    shell(const shell& that);
    // Copy
    shell& operator=(const shell& that);

    // Comparison
    bool operator==(const shell& that) const;
    bool operator<(const shell& that) const;
};

// ext_shell used to preserve volume by a naive approach.
class ext_shell {
  private:
    std::map<int, bool> c0_bool;
    std::map<int, bool> c1_bool;
    std::map<int, bool> c2_bool;

    std::map<int, shell> c0_cor;
    std::map<int, shell> c1_side;
    std::map<int, shell> c2_face;

    // The positions of the shell boundaries.
    std::map<shell, apf::Vector3> c0_cor_pos;
    std::map<shell, apf::Vector3> c1_side_pos;
    std::map<shell, apf::Vector3> c2_face_pos;

    // The directions of the shell sides and faces.
    std::map<shell, apf::Vector3> c1_side_dir;
    std::map<shell, apf::Vector3> c2_face_dir;

  public:
    cell_base sh_base;

    // Check if on the shell. 
    bool chk_shell(int dim, int tag);
    // Return the shell id. 
    shell get_shell(int dim, int tag);

    // Set shell with the boundary id. 
    void set_shell(int dim, int tag, shell sh, bool state = true);
    // Set shell position. 
    void set_pos(shell sh, apf::Vector3 pos);
    // Set shell direction. 
    void set_dir(shell sh, apf::Vector3 dir);

    // Get the shell position for the given shell boundary.
    apf::Vector3 get_shell_pos(shell sh);
    // Get the shell direction for the given shell boundary.
    apf::Vector3 get_shell_dir(shell sh);

    // Get the shell position for the given exterior cell.
    apf::Vector3 get_pos(int dim, int tag);
    // Get the shell direction for the given exterior cell.
    apf::Vector3 get_dir(int dim, int tag);

    // Given an exterior cell:
    // Find the projection of the given direction on the shell.
    apf::Vector3 find_para_dir(int dim, int tag, apf::Vector3 dir);
    // Find the projection of the given position on the shell.
    apf::Vector3 find_int_pos(int dim, int tag, apf::Vector3 pos);
    apf::Vector3 find_int_pos(shell sh, apf::Vector3 pos);

    // Find the displacement from the position on the given direction to inter-
    // sect the shell.
    apf::Vector3 find_int_traj(int dim, int tag, apf::Vector3 pos, 
                                                apf::Vector3 dir);

    // Given a position, return the normal displacement from the position to the 
    // shell.
    apf::Vector3 find_norm_disp(shell sh, apf::Vector3 pos);

    // Given a position, return the position mirrored by the shell.
    apf::Vector3 find_pos_mir(int dim, int tag, apf::Vector3 pos);
    // Given a direction, return the direction mirrored by the shell.
    apf::Vector3 find_dir_mir(int dim, int tag, apf::Vector3 dir);
    apf::Vector3 find_dir_mir(shell sh, apf::Vector3 dir);

    // Used after insertion, to detect the shell id of the newly inserted 
    // exterior cells.
    int find_new_cell_id(int dim, int tag, int tag_0cell);

    void clear();

    ext_shell();
    ~ext_shell();

    // Copy constructor
    ext_shell(const ext_shell& that);
    // Copy
    ext_shell& operator=(const ext_shell& that);
};
class shell_burner {
  private:
    ext_shell* e_sh;
    cell_base* c_base;
    apf::Mesh2* m;
    vd_entlist* e_list;
  public:
    shell_burner(ext_shell* e_sh_in, cell_base* cb_in, apf::Mesh2* m_in, 
                  vd_entlist* e_in);
    void collect_ext(double tol = 0.01);
    void collect_ext_dir();
    void burn_shell(int dim, int tag, shell sh);
};

// field_calc uses this class to update energy and drag coefficients.
// Currently the values are constant for all simulation, except for specified
// 2cells.
class boundary_energy {
  private:

  public:
    std::map<int, bool> c2s_zeros;

    void set_2cells(std::map<int, bool> c2s_in);
    boundary_energy();
    ~boundary_energy();
    boundary_energy& operator=(const boundary_energy& that);

};

class field_calc {
  private:
    // cell_base flag. Exterior velocity fixing requires cell_base object when
    // exterior vertex positions or velocities are projected in some manner.
    bool cb_flag;

    // If the exterior is projected, and the type of calculation.
    bool calc_ext;
    bool calc_corner;
    int proj_flag;

    // TIme integration, Euler method, RK2, RK4
    VEL_TYPE v_flag;
    INTEG_TYPE t_int_flag;
    //EXT_GAM_TYPE ext_gam_flag;

    bool e_sh_flag;
    ext_shell* e_sh;

    cell_base* c_base;

    double drag_glob;
    double drag_rat;
    double d2_glob;

    void set_vel_type(VEL_TYPE DT);
    //vd_eqn_of_motion* EoM;

  public:
    vd_param vdparam;
    boundary_energy b_en;

    vd_eqn_of_motion* EoM;

    double calc_roc(apf::Mesh2* m, apf::MeshEntity* ent);
    double calc_roc(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ent);

    // Calculate the rate of change for a group of entities, where the velocity
    // of a set of vertices are calculated as a single vertex.
    // Used when deciding whether a higher dimensional cell is collapsing after
    // the collapse of an adjacent, lower dimensional cell. 
    // The bounding vertices of the collapsing cell are mapped to a single 
    // vertex. 
    double calc_roc_merg(apf::Mesh2* m, std::vector<apf::MeshEntity*> *ent, 
       std::vector<apf::MeshEntity*>& verts,
       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map);

    apf::Vector3 vd_calc_force(apf::Mesh2* m, apf::MeshEntity* v, bool fix = true);
    apf::Vector3 vd_calc_force_tri(apf::Mesh2* m, apf::MeshEntity* v, 
                    std::vector<apf::MeshEntity*> * tris, bool fix = true);
    apf::Vector3 vd_calc_force_edge(apf::Mesh2* m, apf::MeshEntity* v, 
                    std::vector<apf::MeshEntity*> * tris, bool fix = true);

    void vd_upd_vel_field(apf::Mesh2* m, std::vector<apf::MeshEntity*>* v, 
                                      bool drag_local = false);

    void vd_upd_vel_field(apf::Mesh2* m, apf::MeshEntity* v, 
                                   bool drag_local = false);

    // Considering some of the vertices to be merging, update the velocities.
    void vd_upd_vel_field_merg(apf::Mesh2* m, 
                          std::vector<apf::MeshEntity*> &v, 
                          std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map, 
                          bool drag_local = false);

    void vd_upd_vel_field(apf::Mesh2* m, apf::MeshEntity* v, 
            std::vector<apf::MeshEntity*>* tri_set, bool fix = true);

    void vd_upd_vel_field_edge(apf::Mesh2* m, apf::MeshEntity* v, 
            std::vector<apf::MeshEntity*>* edge_set, bool fix = true);

    void vd_calc_vel(apf::Mesh2* m);
    void vd_calc_vel(apf::Mesh2* m, std::vector<apf::MeshEntity*>* verts);

    // Using Runge Kutta2/4 method, calculate the position and velocity of the 
    // vertices. 
    void vd_calc_rk2(apf::Mesh2* m);
    void vd_calc_rk4(apf::Mesh2* m);

    void vd_calc_rk2(apf::Mesh2* m, std::vector<apf::MeshEntity*>* v, 
                                      bool drag_local);
    void vd_calc_rk4(apf::Mesh2* m, std::vector<apf::MeshEntity*>* v, 
                                      bool drag_local);

    void vd_calc_rk2(apf::Mesh2* m, apf::MeshEntity* v, 
                                      bool drag_local);
    void vd_calc_rk4(apf::Mesh2* m, apf::MeshEntity* v, 
                                      bool drag_local);

    //TODO All these field names and attributes should be generically kept in a 
    // fields list, rather than repeating similar calls in every function.

    void vd_att_fields(apf::Mesh2* m);
    void vd_del_fields(apf::Mesh2* m);
    void vd_att_fields(apf::Mesh2* m, apf::MeshEntity* v);
    void vd_att_tri_fields(apf::Mesh2* m, apf::MeshEntity* tri);

    void vd_att_tri_fields(apf::Mesh2* m);
    void del_tri_fields(apf::Mesh2* m);

    bool vd_chk_fields(apf::Mesh2* m);

    // Using Runge Kutta4 method, calculate the velocity of a set of vertices. 
    // Used before calc_roc to estimate the rate of change less affected by 
    // oscillations in the solution.
    void vd_calc_rk4(apf::Mesh2* m, std::vector<apf::MeshEntity*> &verts);

    // Given two meshes and a map from m1 to m2 entities, transfer the 
    // fields.
    void vd_trans_fields(apf::Mesh2* m1, apf::Mesh2* m2,
                    std::vector<std::vector<apf::MeshEntity*> > & ents,
                    std::map<apf::MeshEntity*, apf::MeshEntity*> &e_map);
    void vd_trans_tri_fields(apf::Mesh2* m1, apf::Mesh2* m2,
                    std::vector<std::vector<apf::MeshEntity*> > & ents,
                    std::map<apf::MeshEntity*, apf::MeshEntity*> &e_map);

    // Time step for the vertex to invert one of the neighboring tets.
    double find_min_t2(apf::Mesh2* m, double t_set = -1);
    // Tetrahedron rate of change of volume based time calculation.
    double find_min_t_tet(apf::Mesh2* m);
    double find_min_t2(apf::Mesh2* m, std::vector<apf::MeshEntity*>* vert, double t_set = -1);
    double find_min_t2(apf::Mesh2* m, apf::MeshEntity* vert, double t_set = -1);

    // Get the maximum of the norms of velocities of the interior vertices.
    // TODO this is to circumnavigate the unphysical high velocity magnitude
    // of the exterior vertices that doesn't affect simulation. Primarly used in
    // VD_sim_bench_GAUSS_evol.cc.
    double vd_get_max_vel_norm_int(apf::Mesh2* m);
    // Used as a wrapper for create_Vert to make sure the fields are 
    // generated correctly.
    apf::MeshEntity* create_v(apf::Mesh2* m, apf::ModelEntity* mdl);

    ////////////////////////////
    // Energy and drag:
    ////////////////////////////
    void set_vdparam(vd_param vdpar_in, apf::Mesh2* m);

    double upd_d2(apf::Mesh2* m);
    double upd_gam2(apf::Mesh2* m);

    // Update the drag coefficient and energy per unit area of the set of 
    // triangles.
    double upd_d2(apf::Mesh2* m, std::vector<apf::MeshEntity*>& tri, double val);
    double upd_gam2(apf::Mesh2* m, std::vector<apf::MeshEntity*>& tri, 
                                                                     double val);

    // Get surface energy of a triangles.
    double gam2(apf::Mesh2* m, apf::MeshEntity* tri);
    // Get drag coefficient of a triangles.
    double d2(apf::Mesh2* m, apf::MeshEntity* tri);

    // Get total surface energy of a triangle.
    double en_tri(apf::Mesh2* m, apf::MeshEntity* tri);

    // Refresh the mesh and vd_elist lists of the EoM.
    void refresh_mesh();

    void set_calc_ext(bool calc_ext_in);
    //, int c3_out_in = -1);
    //int get_c3_out();
    bool get_ext();
    bool get_corner();
    int get_proj();

    ext_shell copy_e_sh();
    ext_shell* get_e_sh();
    void set_e_sh(ext_shell* e_sh_in);

    // Check if the vertex is skipped. 
    bool chk_skip(apf::Mesh2* m, apf::MeshEntity* vert);
    // Check if the vertex velocity requires special processing. 
    bool chk_vert_special(apf::Mesh2* m, apf::MeshEntity* vert);

    bool chk_0c_cor_gmi(int tag_0cell);

    void set_vec_sp_calc(PROJ_TYPE PROJ, ext_shell* e_sh_in = NULL);

    void set_integ_type(INTEG_TYPE IT);

    void set_drag_glob(double drag_in);
    double get_drag_glob();
    double get_drag_rat();
    void set_drag_rat(double rat_in);

    void set_d2_glob(double drag_in);
    double get_d2_glob();

    void dummy_func_stop();
    apf::Vector3 get_vec_special(apf::Mesh2* m, apf::MeshEntity* vert,
                                         apf::Vector3 v_in);

    apf::Vector3 fix_vert_bound(apf::Mesh2* m, apf::MeshEntity* vert, 
                                                              apf::Vector3 pos);
    apf::Vector3 fix_disp_bound(apf::Mesh2* m, apf::MeshEntity* vert,
                       apf::Vector3 pos, apf::Vector3 disp, double scale = 0.9);

    void fix_vel_special(apf::Mesh2* m, apf::MeshEntity* vert);
    void fix_vel_special_merg(apf::Mesh2* m, 
                  std::vector<apf::MeshEntity*> &vert, 
                  std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map);

    void corr_pos(apf::Mesh2* m);
    void corr_pos(apf::Mesh2* m, apf::MeshEntity* e_v);

    void refresh(apf::Mesh2* m_in, cell_base* c_base_in, vd_entlist* e_list_in);

    field_calc();
    field_calc(VEL_TYPE DT);
    ~field_calc();

    void reload_cb(cell_base* c_in);
    void off_cb();

    // Copy constructor
    field_calc(const field_calc& that);
    // Copy
    field_calc& operator=(const field_calc& that);
};

#endif
