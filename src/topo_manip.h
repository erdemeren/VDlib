#ifndef TOPO_MANIP_H
#define TOPO_MANIP_H

#include <algorithm>    

#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>

#include "topo_rand.h"
#include "topo_extinfo.h"
#include "topo_geom.h"

#include "topo_disp.h"

class field_calc;
class vd_param;

struct plane{
  double a;
  apf::Vector3 norm;
  plane(): a(0), norm(0,0,1) {}
  void from_point(apf::Vector3* point, apf::Vector3* norm_dir) {
    norm = *norm_dir;
    norm = norm_0(norm);
    a = norm* *point;
  }
};

struct sphere{
  double r;
  apf::Vector3 ctr;
  sphere(): r(0), ctr(0,0,1) {}
  void from_point(apf::Vector3* point, double r_in) {
    ctr = *point;
    r = r_in;
  }
};

double const pi = 3.14159265358979323846;

//////////////////////////////////////
// Operations modifying the mesh:
//////////////////////////////////////
// TODO Remove outdated

// Perturb the position of a point, without inverting elements.
// Poor practice, but works for good quality tetrahedrons.
void vd_pert_point(apf::Mesh2* m, apf::MeshEntity* v1);

// Perturb the positions of all points in a mesh, without inverting elements.
void vd_pert_mesh(apf::Mesh2* m);

// Perturb the position of a point, by the ratio of the maximum distance along a 
// random direction.
void vd_pert_point_rat(apf::Mesh2* m, apf::MeshEntity* v, double rat = 0.5);
// Perturb the position of all points, by the ratio of the maximum distance along 
// a random direction.
void vd_pert_mesh_rat(apf::Mesh2* m, double rat = 0.5);
// Perturb the position of a set of points, by the ratio of the maximum distance 
// along a random direction.
void vd_pert_set_rat(apf::Mesh2* m, std::vector<apf::MeshEntity*>& v, double rat = 0.5);

// Given a mesh and axis, shrink the mesh.
void vd_shr_axis(apf::Mesh2* m, apf::Vector3 ax, double sp_size);
void vd_shr_axis2(apf::Mesh2* m, apf::Vector3 ax, double sp_size);
apf::Vector3 vd_shr_vec(apf::Vector3 vec, apf::Vector3 ax, double sp_size);

// TODO: Given a mesh and surface geometry, shrink the surface.
void vd_shr_surf(apf::Mesh2* m, int geom, float sp_size);

// Given a mesh and cell, shrink the cell by making the vertices belonging to 
// the cell and its lower dimensional adjacencies closer to the cell center 
// defined as the average position of all belonging vertices.
void vd_shr_cell(apf::Mesh2* m, int dim, int cell, float sp_size);

// Given a cell, try to make the cell shape equiaxed. This is achieved by
// rescaling the positions such that the difference between the maximum and 
// minimum coordinates x,y,z are made the same. It can be also achieved by
// PCA type of analysis.
bool vd_pert_cell_equi(apf::Mesh2* m, int cell_dim, int cell_id);

// Given a cell, project the outer vertices onto a sphere of mean width of the
// cell.
void vd_pert_cell_sphere(apf::Mesh2* m, int cell_dim, int cell_id);

// Given a mesh, project the outer surface points onto a sphere of size sp_size
// and center, and scale other points considering the farthest domain surface 
// point.
void vd_pert_sphere(apf::Mesh2* m, float sp_size);

// Given a mesh, project the vertex onto a sphere of size rho
// and center.
void vd_proj_v_sphere(apf::Mesh2* m, apf::MeshEntity* vert, apf::Vector3 ctr, double rho);

// Destroy the mesh entities inside an entity set.
void vd_dest_set(apf::Mesh2* m, std::vector<apf::MeshEntity*>* set_dest);

// Remove the elements and nodes above the cutting plane.
void vd_cut_mesh(apf::Mesh2* m, plane* cut_plane);

// Tag the entities above the cutting plane.
void vd_tag_cut(apf::Mesh2* m, plane* cut_plane);

// Tag the elements belonging to cells above the cutting plane.
void vd_cell_cut(apf::Mesh2* m, plane* cut_plane);
void vd_cell_cut_sph(apf::Mesh2* m, sphere* cut_sph);

//////////////////////////////////////
// Calculate velocity field:
//////////////////////////////////////
// TODO use an object oriented approach. 
// Input apf::Mesh2*, field_calc*, vd_entlist*
// Functions for calculating velocity and force.
// Optional inputs for a set of vertices, local drag modification.
class vd_eqn_of_motion {
  protected:
    apf::Mesh2* m;
    cell_base* c_base;
    field_calc* f_calc;
    vd_entlist* e_list;
    apf::Field* vel_field;

  public:
    // Update the velocity field at:
    // Every boundary vertex.
    virtual void calc_vel();
    virtual void calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local);
    virtual void calc_vel_curr(apf::MeshEntity* vert);

    // A single vertex.
    virtual void vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local = false);
    virtual apf::Vector3 vd_upd_vel_field_tri(apf::MeshEntity* vert, 
            std::vector<apf::MeshEntity*>* tris, bool drag_local = false);
    virtual apf::Vector3 vd_upd_vel_field_edge(apf::MeshEntity* vert, 
               std::vector<apf::MeshEntity*>* edges, bool drag_local = false);

    // A set of vertices while skipping labeled vertices.
    virtual void vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                           std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                           bool drag_local = false);

    // Calculate the force acting on a vertex.
    virtual apf::Vector3 vd_calc_force(apf::MeshEntity* vert);
    virtual apf::Vector3 vd_calc_force_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris);

    virtual apf::Vector3 vd_calc_force_edge(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* edges);

    vd_eqn_of_motion();
    vd_eqn_of_motion(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in);
    virtual ~vd_eqn_of_motion();

    void refresh(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in);
    void refresh_elist();

    // Copy constructor
    vd_eqn_of_motion(const vd_eqn_of_motion& that);
    // Copy
    vd_eqn_of_motion& operator=(const vd_eqn_of_motion& that);
};

// vd_eqn_mason
class vd_eqn_mason : public vd_eqn_of_motion {
  protected:
    std::vector<apf::MeshEntity*> es_edge;
    std::vector<apf::MeshEntity*> es_surf;

    apf::Up up;
    apf::Downward d_v;
    apf::Downward d_e;

    apf::Vector3 p_i;
    apf::Vector3 p_j;
    apf::Vector3 n_ij;
    apf::Vector3 p_ctr;
    apf::Vector3 temp;
    apf::Vector3 rhs;

    double average;

    apf::Field* vel_field;
    double drag_rat;

  public:
    void get_average_tri();
    void calc_vel_curr(apf::MeshEntity* vert);
    apf::Vector3 calc_vel_curr_tri(apf::MeshEntity* vert, 
             std::vector<apf::MeshEntity*>* tris);

    // Update the velocity field at:
    // Every boundary vertex.
    void calc_vel();
    void calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local);

    // A single vertex.
    void vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_tri(apf::MeshEntity* vert, 
             std::vector<apf::MeshEntity*>* tris, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_edge(apf::MeshEntity* vert, 
               std::vector<apf::MeshEntity*>* edges, bool drag_local = false);

    // A set of vertices while skipping labeled vertices.
    void vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                           std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                           bool drag_local = false);

    // Calculate the force acting on a vertex.
    apf::Vector3 vd_calc_force(apf::MeshEntity* vert);
    apf::Vector3 vd_calc_force_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris);

    apf::Vector3 vd_calc_force_edge(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* edges);
    vd_eqn_mason();
    vd_eqn_mason(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in);
    ~vd_eqn_mason();
    // Copy constructor
    vd_eqn_mason(const vd_eqn_mason& that);
    // Copy
    vd_eqn_mason& operator=(const vd_eqn_mason& that);
};

// Neumann BC
class vd_eqn_mason_NBC : public vd_eqn_of_motion {
  protected:
    std::vector<apf::MeshEntity*> es_edge;
    std::vector<apf::MeshEntity*> es_surf;

    apf::Up up;
    apf::Downward d_v;
    apf::Downward d_e;

    apf::Vector3 p_i;
    apf::Vector3 p_j;
    apf::Vector3 n_ij;
    apf::Vector3 p_ctr;
    apf::Vector3 temp;
    apf::Vector3 rhs;

    double average;

    apf::Field* vel_field;
    double drag_rat;

  public:
    void get_average_tri();
    void calc_vel_curr(apf::MeshEntity* vert);
    apf::Vector3 calc_vel_curr_tri(apf::MeshEntity* vert, 
             std::vector<apf::MeshEntity*>* tris);

    // Update the velocity field at:
    // Every boundary vertex.
    void calc_vel();
    void calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local);

    // A single vertex.
    void vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_tri(apf::MeshEntity* vert, 
             std::vector<apf::MeshEntity*>* tris, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_edge(apf::MeshEntity* vert, 
               std::vector<apf::MeshEntity*>* edges, bool drag_local = false);

    // A set of vertices while skipping labeled vertices.
    void vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                           std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                           bool drag_local = false);

    // Calculate the force acting on a vertex.
    apf::Vector3 vd_calc_force(apf::MeshEntity* vert);
    apf::Vector3 vd_calc_force_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris);

    apf::Vector3 vd_calc_force_edge(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* edges);
    vd_eqn_mason_NBC();
    vd_eqn_mason_NBC(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in);
    ~vd_eqn_mason_NBC();
    // Copy constructor
    vd_eqn_mason_NBC(const vd_eqn_mason_NBC& that);
    // Copy
    vd_eqn_mason_NBC& operator=(const vd_eqn_mason_NBC& that);
};

class vd_eqn_lazar : public vd_eqn_of_motion {
  protected:
    std::vector<apf::MeshEntity*> es_edge;
    std::vector<apf::MeshEntity*> es_surf;

    std::vector<apf::MeshEntity*> b_verts;

    // The 3-stratum considered to be inside when calculating the normal.
    std::map<apf::MeshEntity*, int> tri_tet_in;
    std::map<apf::MeshEntity*, double> edge_len;
    std::map<apf::MeshEntity*, apf::Vector3> tri_norm;
    std::map<apf::MeshEntity*, apf::Vector3> edge_pos;
    std::map<apf::MeshEntity*, apf::Vector3> tri_pos;
    std::map<int, double> a_equi_map;

    apf::Downward d_v;
    apf::Downward d_e;

    apf::Vector3 p_i;
    apf::Vector3 p_j;
    apf::Vector3 n_ij;
    apf::Vector3 p_ctr;
    apf::Vector3 temp;
    apf::Vector3 rhs;

    double average;

    apf::Field* vel_field;
    double drag_rat;

    void calc_vel_curr(apf::MeshEntity* vert);
  public:
    void get_average_tri();

    void set_ang_equi_cell(int c_id);
    void set_ang_equi_cells();
    void mark_edges_global();
    void mark_edges();
    void mark_tris_global();
    void mark_tris();
    void collect_orientations();
    apf::Vector3 get_eij(std::vector<apf::MeshEntity*> &es_bsurf, int c3_curr);
    void chk_surf_grain_int(std::vector<apf::MeshEntity*> &es_tris, std::vector<apf::MeshEntity*> &es_bsurf, int c3_curr);
    // Calculate the L and M quantities used in local Srolovitz-MacPherson.
    double vd_dVdt(apf::MeshEntity* vert, int c3_curr, apf::Vector3 &eij);
    double vd_dVdt_int(apf::MeshEntity* vert, int c3_curr, apf::Vector3 &eij);
    apf::Vector3 calc_dn_quad(apf::MeshEntity* vert, 
           std::vector<apf::Vector3>& eij, std::vector<double>& dVdt);
    apf::Vector3 calc_dn_surf(apf::MeshEntity* vert, 
           std::vector<apf::Vector3>& eij, std::vector<double>& dVdt);
    apf::Vector3 calc_dn_trip(apf::MeshEntity* vert, 
           std::vector<apf::Vector3>& eij, std::vector<double>& dVdt);

    apf::Vector3 calc_vel_curr_tri(apf::MeshEntity* vert, 
             std::vector<apf::MeshEntity*>* tris);

    // Update the velocity field at:
    // Every boundary vertex.
    void calc_vel();
    void calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local);

    // A single vertex.
    void vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_tri(apf::MeshEntity* vert, 
               std::vector<apf::MeshEntity*>* tris, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_edge(apf::MeshEntity* vert, 
               std::vector<apf::MeshEntity*>* edges, bool drag_local = false);

    // A set of vertices while skipping labeled vertices.
    void vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                           std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                           bool drag_local = false);

    // Calculate the force acting on a vertex.
    apf::Vector3 vd_calc_force(apf::MeshEntity* vert);
    apf::Vector3 vd_calc_force_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris);

    apf::Vector3 vd_calc_force_edge(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* edges);
    vd_eqn_lazar();
    vd_eqn_lazar(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in);
    ~vd_eqn_lazar();
    // Copy constructor
    vd_eqn_lazar(const vd_eqn_lazar& that);
    // Copy
    vd_eqn_lazar& operator=(const vd_eqn_lazar& that);
};

class vd_eqn_lazar_NBC : public vd_eqn_of_motion {
  protected:
    std::vector<apf::MeshEntity*> es_edge;
    std::vector<apf::MeshEntity*> es_surf;

    std::vector<apf::MeshEntity*> b_verts;

    // The 3-stratum considered to be inside when calculating the normal.
    std::map<apf::MeshEntity*, int> tri_tet_in;
    std::map<apf::MeshEntity*, double> edge_len;
    std::map<apf::MeshEntity*, apf::Vector3> tri_norm;
    std::map<apf::MeshEntity*, apf::Vector3> edge_pos;
    std::map<apf::MeshEntity*, apf::Vector3> tri_pos;
    std::map<int, double> a_equi_map;

    apf::Downward d_v;
    apf::Downward d_e;

    apf::Vector3 p_i;
    apf::Vector3 p_j;
    apf::Vector3 n_ij;
    apf::Vector3 p_ctr;
    apf::Vector3 temp;
    apf::Vector3 rhs;

    double average;

    apf::Field* vel_field;
    double drag_rat;

    void calc_vel_curr(apf::MeshEntity* vert);
  public:
    void get_average_tri();

    void set_ang_equi_cell(int c_id);
    void set_ang_equi_cells();
    void mark_edges_global();
    void mark_edges();
    void mark_tris_global();
    void mark_tris();
    void collect_orientations();
    apf::Vector3 get_eij(std::vector<apf::MeshEntity*> &es_bsurf, int c3_curr);
    void chk_surf_grain_int(std::vector<apf::MeshEntity*> &es_tris, std::vector<apf::MeshEntity*> &es_bsurf, int c3_curr);
    // Calculate the L and M quantities used in local Srolovitz-MacPherson.
    double vd_dVdt(apf::MeshEntity* vert, int c3_curr, apf::Vector3 &eij);
    double vd_dVdt_int(apf::MeshEntity* vert, int c3_curr, apf::Vector3 &eij);
    apf::Vector3 calc_dn_quad(apf::MeshEntity* vert, 
           std::vector<apf::Vector3>& eij, std::vector<double>& dVdt);
    apf::Vector3 calc_dn_surf(apf::MeshEntity* vert, 
           std::vector<apf::Vector3>& eij, std::vector<double>& dVdt);
    apf::Vector3 calc_dn_trip(apf::MeshEntity* vert, 
           std::vector<apf::Vector3>& eij, std::vector<double>& dVdt);

    apf::Vector3 calc_vel_curr_tri(apf::MeshEntity* vert, 
             std::vector<apf::MeshEntity*>* tris);

    // Update the velocity field at:
    // Every boundary vertex.
    void calc_vel();
    void calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local);

    // A single vertex.
    void vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_tri(apf::MeshEntity* vert, 
               std::vector<apf::MeshEntity*>* tris, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_edge(apf::MeshEntity* vert, 
               std::vector<apf::MeshEntity*>* edges, bool drag_local = false);

    // A set of vertices while skipping labeled vertices.
    void vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                           std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                           bool drag_local = false);

    // Calculate the force acting on a vertex.
    apf::Vector3 vd_calc_force(apf::MeshEntity* vert);
    apf::Vector3 vd_calc_force_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris);

    apf::Vector3 vd_calc_force_edge(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* edges);
    vd_eqn_lazar_NBC();
    vd_eqn_lazar_NBC(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in);
    ~vd_eqn_lazar_NBC();
    // Copy constructor
    vd_eqn_lazar_NBC(const vd_eqn_lazar_NBC& that);
    // Copy
    vd_eqn_lazar_NBC& operator=(const vd_eqn_lazar_NBC& that);
};

class vd_eqn_kuprat : public vd_eqn_of_motion {
  protected:
    std::vector<apf::MeshEntity*> es_edge;
    std::vector<apf::MeshEntity*> es_surf;

    apf::Up up;
    apf::Downward d_v;
    apf::Downward d_e;

    apf::Vector3 p_i;
    apf::Vector3 p_j;
    apf::Vector3 n_ij;
    apf::Vector3 p_ctr;
    apf::Vector3 temp;
    apf::Vector3 rhs;

    double average;

    apf::Field* vel_field;
    double drag_rat;

  public:
    void get_average_tri();
    void calc_vel_curr(apf::MeshEntity* vert);
    apf::Vector3 calc_vel_curr_tri(apf::MeshEntity* vert, 
             std::vector<apf::MeshEntity*>* tris);

    // Update the velocity field at:
    // Every boundary vertex.
    void calc_vel();
    void calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local);
    // A single vertex.
    void vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_tri(apf::MeshEntity* vert, 
             std::vector<apf::MeshEntity*>* tris, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_edge(apf::MeshEntity* vert, 
               std::vector<apf::MeshEntity*>* edges, bool drag_local = false);

    // A set of vertices while skipping labeled vertices.
    void vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                           std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                           bool drag_local = false);

    // Calculate the force acting on a vertex.
    apf::Vector3 vd_calc_force(apf::MeshEntity* vert);
    apf::Vector3 vd_calc_force_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris);

    apf::Vector3 vd_calc_force_edge(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* edges);
    vd_eqn_kuprat();
    vd_eqn_kuprat(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in);
    ~vd_eqn_kuprat();
    // Copy constructor
    vd_eqn_kuprat(const vd_eqn_kuprat& that);
    // Copy
    vd_eqn_kuprat& operator=(const vd_eqn_kuprat& that);
};

class vd_eqn_kuprat_NBC : public vd_eqn_of_motion {
  protected:
    std::vector<apf::MeshEntity*> es_edge;
    std::vector<apf::MeshEntity*> es_surf;

    apf::Up up;
    apf::Downward d_v;
    apf::Downward d_e;

    apf::Vector3 p_i;
    apf::Vector3 p_j;
    apf::Vector3 n_ij;
    apf::Vector3 p_ctr;
    apf::Vector3 temp;
    apf::Vector3 rhs;

    double average;

    apf::Field* vel_field;
    double drag_rat;

  public:
    void get_average_tri();
    void calc_vel_curr(apf::MeshEntity* vert);
    apf::Vector3 calc_vel_curr_tri(apf::MeshEntity* vert, 
             std::vector<apf::MeshEntity*>* tris);

    // Update the velocity field at:
    // Every boundary vertex.
    void calc_vel();
    void calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local);
    // A single vertex.
    void vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_tri(apf::MeshEntity* vert, 
             std::vector<apf::MeshEntity*>* tris, bool drag_local = false);
    apf::Vector3 vd_upd_vel_field_edge(apf::MeshEntity* vert, 
               std::vector<apf::MeshEntity*>* edges, bool drag_local = false);

    // A set of vertices while skipping labeled vertices.
    void vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                           std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                           bool drag_local = false);

    // Calculate the force acting on a vertex.
    apf::Vector3 vd_calc_force(apf::MeshEntity* vert);
    apf::Vector3 vd_calc_force_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris);

    apf::Vector3 vd_calc_force_edge(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* edges);
    vd_eqn_kuprat_NBC();
    vd_eqn_kuprat_NBC(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in);
    ~vd_eqn_kuprat_NBC();
    // Copy constructor
    vd_eqn_kuprat_NBC(const vd_eqn_kuprat_NBC& that);
    // Copy
    vd_eqn_kuprat_NBC& operator=(const vd_eqn_kuprat_NBC& that);
};


// PLACEHOLDER object for eqn of motion calculation using mirror planes. 
// Inherit from eqn of motion object potentially...
// break down the equation of motion calculation for Mason into parts.
// Per triangle. 

// PLACEHOLDER evaluate the forces and drag matrices of a triangle by mirroring 
// by a plane or not. Store the current values in the object.
class vd_eqn_mason_mir : public vd_eqn_of_motion {
  protected:
    std::vector<apf::MeshEntity*> es_edge;
    std::vector<apf::MeshEntity*> es_surf;

    apf::Up up;
    apf::Downward d_v;
    apf::Downward d_e;

    apf::Vector3 p_i;
    apf::Vector3 p_j;
    apf::Vector3 n_ij;
    apf::Vector3 p_ctr;
    apf::Vector3 temp;
    apf::Vector3 rhs;

    double average;

    apf::Field* vel_field;
    double drag_rat;

  public:
    void get_average_tri();
    void calc_vel_curr(apf::MeshEntity* vert);

    // Update the velocity field at:
    // Every boundary vertex.
    void calc_vel();
    void calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local);

    // A single vertex.
    void vd_upd_vel_field(apf::MeshEntity* vert);

    // A set of vertices while skipping labeled vertices.
    void vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                           std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                           bool drag_local = false);

    // Calculate the force acting on a vertex.
    apf::Vector3 vd_calc_force(apf::MeshEntity* vert);
    apf::Vector3 vd_calc_force_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris);

    apf::Vector3 vd_calc_force_edge(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* edges);
    vd_eqn_mason_mir();
    vd_eqn_mason_mir(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in);
    ~vd_eqn_mason_mir();
    // Copy constructor
    vd_eqn_mason_mir(const vd_eqn_mason_mir& that);
    // Copy
    vd_eqn_mason_mir& operator=(const vd_eqn_mason_mir& that);
};

typedef std::vector<apf::Matrix3x3 > v_rot;
typedef std::vector<apf::Vector3> v_mir;

#define V_MIR_INIT std::vector<apf::Vector3> (0, apf::Vector3(0,0,0))
#define V_ROT_INIT std::vector<apf::Matrix3x3> (0, apf::Matrix3x3(0,0,0,0,0,0,0,0,0))
// Exterior strata are defined mirror and rotation operations. The force applying
// triangles are rotated and reflected to be included in the velocity and force
// calculations. Currently used to find the equilibrium metastable shape of the
// truncated octahedron by only moving the local environment of a single 0-stratum
class vd_eqn_mason_mirrot : public vd_eqn_of_motion {
  protected:
    // Rotations and mirrors around the vertex on the boundary.
    // The entities around a corresponding boundary vertex are mirrored/rotated 
    // per set of mirrors/rotations
    // dim - id - vector of vector rotation matrices
    std::vector<std::vector<std::vector< v_rot > > > rot;
    // dim - id - vector of vector mirror axes
    std::vector<std::vector<std::vector< v_mir > > > mir;

    std::vector<apf::MeshEntity*> es_edge;
    std::vector<apf::MeshEntity*> es_surf;

    apf::Up up;
    apf::Downward d_v;
    apf::Downward d_e;

    apf::Vector3 p_i;
    apf::Vector3 p_j;
    apf::Vector3 n_ij;
    apf::Vector3 p_ctr;
    apf::Vector3 temp;
    apf::Vector3 rhs;

    double average;

    apf::Field* vel_field;
    double drag_rat;

  public:
    void get_average_tri();
    void calc_vel_curr(apf::MeshEntity* vert);

    // Update the velocity field at:
    // Every boundary vertex.
    void calc_vel();
    void calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local);

    // A single vertex.
    void vd_upd_vel_field(apf::MeshEntity* vert);

    // A set of vertices while skipping labeled vertices.
    void vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                           std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                           bool drag_local = false);

    // Calculate the force acting on a vertex.
    apf::Vector3 vd_calc_force(apf::MeshEntity* vert);
    apf::Vector3 vd_calc_force_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris);

    apf::Vector3 vd_calc_force_edge(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* edges);

    // Set the mirror and rotation operators for the boundary strata.
    void set_rot(std::vector<std::vector<std::vector< std::vector<apf::Matrix3x3 > > > > &rot_in);
    void set_mir(std::vector<std::vector<std::vector< std::vector<apf::Vector3 > > > > &mir_in);

    vd_eqn_mason_mirrot();
    vd_eqn_mason_mirrot(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in);
    ~vd_eqn_mason_mirrot();
    // Copy constructor
    vd_eqn_mason_mirrot(const vd_eqn_mason_mirrot& that);
    // Copy
    vd_eqn_mason_mirrot& operator=(const vd_eqn_mason_mirrot& that);
};

// A very large quality value for ill defined tetrahedra.
#define Q_TET_UNDEF 10e8

// Kuprat tetrahedra quality increasing forces. 
class vd_eqn_kuprat_volm : public vd_eqn_of_motion {
  protected:
    // The volume, area, lengths of entities. For calculation of quality.
    std::vector<std::map<apf::MeshEntity*, double> > VAL;
    // Triangle normal direction and positions
    std::map<apf::MeshEntity*, apf::Vector3 > a_norm;
    std::map<apf::MeshEntity*, apf::Vector3 > a_pos;
    // Tet quality. Varies from 324 to infinity for equilateral to badly shaped
    // tetrahedra. 
    std::map<apf::MeshEntity*, double> q_tet;

    apf::Up up;
    apf::Downward d_s;
    apf::Downward d_e;
    apf::Downward d_v;

    std::vector<apf::MeshEntity*> es_edge;
    std::vector<apf::MeshEntity*> es_surf;
    std::vector<apf::MeshEntity*> es_tet;

    // The quality restoring force is 2*q_tet/h_i * t_i
    // q_tet: quality of th e current tet.
    // h_i: Altitude of the tet w.r.t. current vertex
    // t_i: Altitude direction
    apf::Vector3 p_i;

    apf::Field* vel_field;
    double drag_rat;
    double average;

  public:
    void collect_VAL();
    double calc_quality(apf::MeshEntity* tet);
    void collect_q();

    void get_average_tri();

    // Update the velocity field at:
    // Every boundary vertex.
    void calc_vel_curr(apf::MeshEntity* vert);
    void calc_vel();
    void calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local);
    // A single vertex.
    void vd_upd_vel_field(apf::MeshEntity* vert);

    // A set of vertices while skipping labeled vertices.
    void vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                           std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                           bool drag_local = false);

    // Calculate the force acting on a vertex.
    apf::Vector3 vd_calc_force(apf::MeshEntity* vert);
    apf::Vector3 vd_calc_force_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris);

    apf::Vector3 vd_calc_force_edge(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* edges);

    vd_eqn_kuprat_volm();
    vd_eqn_kuprat_volm(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in);
    ~vd_eqn_kuprat_volm();
    // Copy constructor
    vd_eqn_kuprat_volm(const vd_eqn_kuprat_volm& that);
    // Copy
    vd_eqn_kuprat_volm& operator=(const vd_eqn_kuprat_volm& that);
};


//////////////////////////////////////
// Lazar algorithm:
//////////////////////////////////////
// Given dVn(allowed volumetric change), eij (sum of surface area vectors),
// using Lazar algorithm,
// calculate velocities for the surface vertices:
apf::Vector3 vd_calc_dn_surf (double* dVn, const apf::Vector3* eij);

// calculate velocities for the triple line vertices:
apf::Vector3 vd_calc_dn_trip (double* dVn, const apf::Vector3* eij);

// calculate velocities for the quadruple vertices:
apf::Vector3 vd_calc_dn_quad (double* dVn, const apf::Vector3* eij);

// Calculate the L and M quantities used in local Srolovitz-MacPherson.
void vd_calc_lm(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es_edge, int geom, double* lm);

// Given an entity set of triangles forming a surface, find the average area.
void vd_get_eij(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es_surf, int geom, apf::Vector3* eij);

// Calculate the dVn and eij quantities used in Lazar algorithm.
// Don't call with internal elements, although it is also defined for that. 
// LOOK repeated code, simplify.
void vd_calc_dVn(apf::Mesh2* m, apf::MeshEntity* vert, double mult, 
                 double* dVn, apf::Vector3* eij);

// Calculate velocity for the boundary vertices, using Lazar algorithm.
void vd_calc_vel_lazar (apf::Mesh2* m, vd_param* vd_par, bool ext = false);
void vd_calc_vel_lazar (apf::Mesh2* m, field_calc* f_calc);

// Update the velocity field at a vertex.
void vd_upd_vel_field_lazar(apf::Mesh2* m, apf::MeshEntity* vert, 
                              vd_param* par, bool ext = false);
void vd_upd_vel_field_lazar(apf::Mesh2* m, apf::MeshEntity* vert, field_calc* f_calc, bool drag_local);

//////////////////////////////////////
// Mason algorithm:
//////////////////////////////////////
// Calculate velocity for the boundary vertices, using Lazar algorithm.
void vd_calc_vel_mason(apf::Mesh2* m, vd_param* vd_par, bool ext = false);
void vd_calc_vel_mason (apf::Mesh2* m, field_calc* f_calc);

// Update the velocity field at a vertex.
void vd_upd_vel_field_mason(apf::Mesh2* m, apf::MeshEntity* vert, 
                              vd_param* par, bool drag_local = false);
void vd_upd_vel_field_mason(apf::Mesh2* m, apf::MeshEntity* vert, field_calc* f_calc, bool drag_local = false);

// Considering a set of vertices to be merging, update the velocities. 
// Consider the adjacencies of the merging vertices to be joint.
void vd_upd_vel_field_mason_merg(apf::Mesh2* m, 
      std::vector<apf::MeshEntity*> &vert,
      std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
      field_calc* f_calc, bool drag_local = false);

// Calculate the force acting on a vertex.
apf::Vector3 vd_calc_force_mason(apf::Mesh2* m, apf::MeshEntity* vert, 
                                                      field_calc* f_calc);
apf::Vector3 vd_calc_force_tri_mason(apf::Mesh2* m, apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris, field_calc* f_calc);

apf::Vector3 vd_calc_force_edge_mason(apf::Mesh2* m, apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges, field_calc* f_calc);

// Update the velocity field at a vertex, using only a set of triangles.
void vd_upd_vel_field_mason_tri(apf::Mesh2* m, apf::MeshEntity* vert, std::vector<apf::MeshEntity*>* tri_set, field_calc* f_calc, bool drag_local = false);

// Update the velocity field at a vertex, using only a set of triangles.
void vd_upd_vel_field_mason_edge(apf::Mesh2* m, apf::MeshEntity* vert, std::vector<apf::MeshEntity*>* edge_set, field_calc* f_calc, bool drag_local = false);

//////////////////////////////////////
// Kuprat algorithm:
//////////////////////////////////////

//////////////////////////////////////
// Apply velocity field:
//////////////////////////////////////
// TODO there is an actual function for that in SCOREC lib:
// void displaceMesh(Mesh2* m, Field* d, double factor=1.0); 
// Apply the displacement defined by the velocity field and time multiplier.  
// Can be used to revert positions.
void vd_apply_vel_field(apf::Mesh2* m, double mult = 1);

// TODO create a copy of the velocity field with only the set vertices assigned
// non zero values and use displaceMesh
// void displaceMesh(Mesh2* m, Field* d, double factor=1.0); 
// Apply the displacements to a set of vertices.
void vd_apply_vel_field(apf::Mesh2* m, std::vector<apf::MeshEntity* >* es, double mult = 1);

//////////////////////////////////////
// Time/multiplier required to invert elements:
//////////////////////////////////////
// Find the magnitude of the maximum displacement in the displacement field.
double vd_find_max_vel(apf::Mesh2* m);
double vd_find_min_t(apf::Mesh2* m);


// Tetrahedron based time step estimation:
double vd_calc_tet_t(apf::Mesh2* m, apf::MeshEntity* tet, 
                                              double &V, double &C);
// Calculate the time step required to invert the tetrahedron, based on the 
// first order expansion of the derivative of volume w.r.t. velocities of the
// vertices.
double vd_find_min_t_tet(apf::Mesh2* m, bool ext);

double vd_find_min_t_all(apf::Mesh2* m, bool ext, double t_set);

double vd_find_min_t2(apf::Mesh2* m, bool ext = false, double t_set = -1);
double vd_find_min_t2(apf::Mesh2* m, std::vector<apf::MeshEntity*>* vert, bool ext = false, double t_set = -1);

// Calculate the highest multiplier for vector field that would invert any 
// tetrahedra.
double vd_find_min_mult(apf::Mesh2* m, std::vector<apf::MeshEntity*>* vert,  
                        std::map<apf::MeshEntity*, bool> &tet_skip,
                        char const* vec_name = "velocity_field");

// Calculate the highest multiplier for vector field that would invert any 
// tetrahedra, that is not marked as to be skipped.
double vd_find_min_mult(apf::Mesh2* m, std::vector<apf::MeshEntity*>* vert, 
                        std::map<apf::MeshEntity*, bool> &tet_skip,
                        std::vector<apf::Vector3> & vect);

double vd_find_min_t2(apf::Mesh2* m, apf::MeshEntity* vert, bool ext = false, double t_set = -1);

// Given a vertex and a set of velocities, find the minimum allowable time.
double vd_find_min_t_v(apf::Mesh2* m, apf::MeshEntity* v, 
                             std::vector<apf::Vector3>* a_pos, 
                             std::vector<apf::Vector3>* area, 
                             std::vector<apf::Vector3>* v3);

//////////////////////////////////////
// Calculate rate of change:
//////////////////////////////////////
// First order expansion of the rate of change of volume. 
double vd_calc_roc(apf::Mesh2* m, apf::MeshEntity* ent, field_calc* f_calc);
double vd_calc_roc(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ent, 
                                                        field_calc* f_calc);

// Higher order terms included, assuming steady state velocity. 
double vd_calc_roc_vol_high(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ent, 
                                                        field_calc* f_calc, double dt);

//////////////////////////////////////
// Add fields:
//////////////////////////////////////
// Add a zero displacement field at a vertex.
void vd_add_0_disp_field(apf::Mesh2* m, apf::MeshEntity* vert);

// Add a zero velocity field at a vertex.
void vd_add_0_vel_field(apf::Mesh2* m, apf::MeshEntity* vert);

// Get velocity field at a vertex.
apf::Vector3 vd_get_v_vel_field(apf::Mesh2* m, apf::MeshEntity* vert);

// Get max norm of velocity.
double vd_get_max_vel_norm(apf::Mesh2* m);

//////////////////////////////////////
// NAN checks:
//////////////////////////////////////
bool isvecnan(apf::Vector3 const& v_in);
void itisnan();

//////////////////////////////////////
// Verbose operations:
//////////////////////////////////////
// Going over the tetrahedra, print the volume, edge lengths and velocities of 
// the vertices.
void vd_tet_vel(apf::Mesh2* m);


#endif
