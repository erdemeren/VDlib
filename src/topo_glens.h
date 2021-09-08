#ifndef TOPO_GLENS_H
#define TOPO_GLENS_H
#include <math.h>
#include <vector>

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
#include "topo_topo.h"

#include "topo_energy.h"

#include "topo_lens.h"

#include "topo_disp.h"
#include "topo_field.h"



// Adaptation of Jeremy's implementation of conjugate gradient with quadratic
// approximation of the function in the direction of the gradient for faster
// convergence. 
#define MAX_CG_ITER 250
#define MAX_LS_ITER 100

#define E_TOL 1.4901e-08
#define G_TOL 1.4901e-08

#define MAX_BND_STEP 5.
#define BND_MAG 1.618033988749895
#define PHI_SQ_INV 0.381966011250105

#define SQRT3_DIV2 0.866025403784439 
#define SQRT_EPS 1.4901e-08
#define EPS 2.2204e-16

#define SHIFT2(a, b, c) do { (a) = (b); (b) = (c); } while (0)
#define SHIFT3(a, b, c, d) do { (a) = (b); (b) = (c); (c) = (d); } while (0)

//
//before loop, initial volumes:
//vol_min = min(abs(V_j))
//vol_neg = min(V_j)

//vol_neg >= 0
//  vol_neg = -vol_scale

//w = -80*vol_scale/vol_neg
//v_th = abs(v_th_tol*vol_min/vol_scale)

//So, if min(V_j) > 0, w = 80

//In any case, even for negative volumes, A_j <= exp(80)

//during loop, volumes after collapse:

//if(w < -8*vol_scale/vol_neg or w > -160*vol_scale/vol_neg)
//  w = -80*vol_scale/vol_neg
//  v_th = abs(v_th_tol*vol_min/vol_scale)

//  vol_min = min(abs(V_j))
//  vol_neg = min(V_j)


// Container for conjugate gradient relaxation related variables.
class vd_relax_cont{
  private:
  public:
    apf::Mesh2* mesh;
    cell_base* c_base;
    field_calc* f_calc;
    bool ext_cell;
    bool upd_shell;
    apf::Vector3 midpoint;

    int dim_col;
    int id_col;

    std::vector<std::vector<apf::MeshEntity*> > ent;
    std::vector<apf::MeshEntity*> v_sp_set;

    std::vector<apf::MeshEntity*> ee;
    std::vector<apf::MeshEntity*> es;
    std::vector<apf::MeshEntity*> et;

    // Collapsing: 1, merging: 2, surrounding: 3.
    std::map<apf::MeshEntity*, int> e_map;

    std::map<apf::MeshEntity*, bool> tet_surr;
    std::map<apf::MeshEntity*, bool> tet_skip;
    std::map<apf::MeshEntity*, bool> tet_inv;
    std::map<apf::MeshEntity*, int> tet_i1;

    std::map<apf::MeshEntity*, bool> v_skip;
    std::map<apf::MeshEntity*, int> t_nbr;

    // Starting volume
    std::vector<double> vol_t;
    std::vector<std::vector<apf::MeshEntity*> > tet_v;
    std::vector<std::vector<int> > tet_v_id;

    std::vector<std::vector<int> > v_tet;
    // Adjacent split vertices of split vertices on the split surface.
    std::vector<std::vector<int> > v_vert;

    std::map<apf::MeshEntity*, apf::Vector3> pos_ori;
    std::map<apf::MeshEntity*, apf::Vector3> v_pos;
    std::map<apf::MeshEntity*, apf::Vector3> v_disp;

    std::vector<double> A;
    std::vector<double> C;
    std::vector<std::vector<apf::Vector3> > B;
    std::vector<double> exp_pw_tot;

    apf::Numbering* tagnumbering;
    apf::Numbering* tetset;

    apf::Field* f_field;

    double vol_scale;
    double vol_max;
    double vol_min;
    // Most negative volume, used to scale w, st. w*vol_t/vol_scale = 100;
    double vol_neg;
    double vol_neg_0;
    double v_th;
    double v_th_tol;
    double w;
    double r_edge_min;
    double epsilon;

    // Conjugate gradient related:
    std::vector<apf::Vector3> x0;
    std::vector<apf::Vector3> x1;

    std::vector<apf::Vector3> g0;
    std::vector<apf::Vector3> g1;
    std::vector<apf::Vector3> dir;
    std::vector<apf::Vector3> ndir;

    bool inv_flag;

    int iter_limit;
    int iter;
    double g_norm;
    double xa, xb, xc, xd, xe, xu;
    double fa, fb, fc, fd, fe, fu;
    double p, q, r, s, t, m, tol, tol2;
    bool inv_quad_step;
    double Phi0, Phi1;

    // Destructor:
    vd_relax_cont(apf::Mesh2* m_in, cell_base* cb_in, field_calc* f_calc_in,
                  int dim_in, int id_in,
                  std::vector<apf::MeshEntity*> &v_sp_set_in,
                  std::vector<std::vector<apf::MeshEntity*> > &ent_in,
                  std::vector<apf::MeshEntity*> &et_in, 
                  std::map<apf::MeshEntity*, int> &e_map_in,
                  bool ext_cell_in, 
                  bool upd_shell_in,
                  apf::Vector3 midpoint_in, double r_in, 
                  double v_th_tol_in = 0.00001);

    // Calculate the terms in the softmax function. B is the 1/6 of the sum of 
    // the inward pointing normals of the triangles adjacent to vertex v.
    apf::Vector3 calc_B(std::vector<apf::MeshEntity*> &v_set, 
            apf::MeshEntity* v, std::map<apf::MeshEntity*, apf::Vector3> &pos);

    // The rate of change of volume wrt. alfa. If v is specified, only consider
    // the residual at that vertex. iith vertex is fixed.
    double calc_C(std::vector<apf::MeshEntity*> &v, 
                          std::map<apf::MeshEntity*, apf::Vector3> &pos, int ii, 
                                                apf::MeshEntity* v_r = NULL);

    bool update_vol_inv_grad();

    // Calculate the potential.
    double calc_energy();
    void calc_Cs();
    void upd_grad(std::vector<apf::Vector3> &g_curr);

    void shift_v_sp_pos(std::vector<apf::Vector3> &x_curr,
                        double x_in);
    // Revert the split vertices to their original positions.
    void revert_ori();
    // Try to find a hull for non-inverting tets by conjugate gradient method.
    bool relax();

    ~vd_relax_cont();
};

/*
- Shrink the cell.
- Split all edges bounded by merging vertices not belonging to the collapsing
  cells. Use the sphere center as the centroid position.
- Find the cell membership of the last remaining vertex. If the collapsing 
  cell is bounded by a 0cell, make the remaining vertex member of the first
  0cell. Otherwise, going over the surrounding entities, check for cell 
  membership. Check 1/2/3cell membership. If there is a n-cell n-ent, there 
  should be only that n-cell. Use the lowest dim cell as replacement.
  Create the remaining vertex. Position at the centroid position.

- Looping over collapsing edges, assign each entity a merging or collapsing 
  flag. If collapsing, don't assign a merging flag.
  Collapsing entities will simply be removed. 
  Merging vertices will have the remaining vertex as replacement.
  Merging entities will have a replacing cell membership. For the merging 
  entities, call the rebuildElement using the new vertex to replace the merging
  vertices. There should only be a single vertex replaced.
  Map the old entity to the new entity. 
  Surrounding entities will be rebuilt in the same manner with replacement map.
  Remove all collapsing entities. Remove all merging and surrounding entities 
  that are not the same as their mapped entity.
*/

// It can also be connected to edge collapse operations.
class vd_glens{
  private:
    apf::Mesh2* m;
    gmi_model* mdl;

    struct cell_base* c_base;

    vd_lens* lens_split;
    vd_bipy* bipy_split;
    vd_sp_tet* tet_split;

    vd_param vd_par;

    bool vrfy_msh;
    bool inv_flag;
    bool shell_flag;
    bool precond_flag;

    // Flag to write VTK files after collapse.
    bool save_vtk;
    bool save_vtk_sub;
    std::string vtk_name;
    std::string smb_name;

    int cell_col_dim;
    int cell_col_id;
    apf::ModelEntity* mdl_col;

    bool c_flag;

    // Modify c_base during collapse or not. Set to false by vd_glens_trial 
    // during trial collapse.
    bool mod_c_base;

    bool ext_cell;
    std::vector<int> ext_c0;
    std::vector<apf::Vector3> ext_c0_vec;

    bool c0_flag;

    bool upd_shell;
    field_calc* f_calc;
    shell sh_min;

    // Fields used in creating the convex hull.
    apf::Field* ch_f;
    apf::Field* sp_f;
    apf::Field* a_f;
    // Force field used in relaxation.
    apf::Field* f_field;

    vd_entlist* e_list;

    std::vector<std::vector<apf::MeshEntity*> > ent;
    std::vector<std::vector<apf::ModelEntity*> > m_list;
    std::vector<std::vector<std::vector<apf::MeshEntity* > > > v_list;
    //std::vector<std::vector<apf::MeshEntity*> > ent_col;
    //std::vector<std::vector<apf::MeshEntity*> > ent_surr;

    // External edges that are collapsing:
    std::vector<apf::MeshEntity*> e_col_ext;

    std::map<apf::MeshEntity*, apf::Vector3> pos_merg;

    apf::MeshEntity* vert_ctr;
    apf::ModelEntity* vert_ctr_em;

    int type_vertex;
    int tag_vertex;
    // The id of the final vertex, if provided any. Otherwise tag_0c = -1.
    int tag_0c;

    // The ratio of the minimum edge length and the distance from the merging
    // vertices to the split vertices when creating the convex hull.
    double cvx_rat;
    // The ratio of the longest edge and the distance from the merging
    // vertices to the split vertices when creating the convex hull. If the edge
    // to be split is smaller than double of ratio times longest edge length, 
    // use half of the edge length instead.
    double cvx_rat2;

    // Collapsing: 1, merging: 2, surrounding: 3.
    std::map<apf::MeshEntity*, int> e_map;

    // If vertex is merging, map true.
    std::map<apf::MeshEntity*, bool> v_map;
    std::map<apf::MeshEntity*, bool> k_map;
    // Merging entities and associated remaining cell.
    std::map<apf::MeshEntity*, apf::ModelEntity*> m_map;

    // To obtain a convex hull, control the edge split ordering.
    // Map of edges to be split in the first loop.
    //std::map<apf::MeshEntity*, bool> e_surr_map;

    // Middle point of all vertices belonging to the cell to collapse.
    apf::Vector3 midpoint;

/*
    // The flag to determine how the split positions are set. By relaxation or 
    // by transfering from the trial mesh.
    bool sp_trial_load;
    // Split edge to split vertex map. Used to transfer split vertex positions
    // from the trial mesh to the actual mesh:
    std::map<apf::MeshEntity*, apf::MeshEntity*> sp_e2v_map;
    // Map from the edges on the actual mesh to the trial mesh.
    std::map<MeshEntity*, apf::MeshEntity*> e2e_map;
    std::map<MeshEntity*, apf::Vector3> &e_sp_pos_trial;
    // The positions of the split vertices, assoiated with the to be split edges.
    std::map<apf::MeshEntity*, apf::Vector3> e_sp_pos;
    // Given a map from edges on the actual mesh to the trial mesh, the map from
    // split edges to the positions of the split vertices on the trial mesh, 
    // upload the split vertex positions associated with the actual edges.
    void get_v_sp_pos(std::map<MeshEntity*, apf::MeshEntity*> &e2e_map_in, 
                      std::map<MeshEntity*, apf::Vector3> &e_sp_pos_trial_in); 
    bool set_v_sp_pos();
*/
    // Get the merging entities and their replacement. 
    void set_merg_map(apf::MeshEntity* e1, apf::MeshEntity* e2);
    void set_merg();
    void set_merg2();

    // Load the entities surrounding the collapsing entity.
    bool load_cell(); // TODO should be called from load_edge called within set_ent. Should work on a per lens basis.
 
    apf::MeshEntity* vd_rebuildElement(apf::MeshEntity* original,
                                              apf::BuildCallback* cb = 0);

    void get_vert();
    void recreate_ent();

    // Destroy the old entities after the rebuild operations.
    void destroy_ent();

    // CHECK AND UPDATE FUNCTIONS USED DURING COLLAPSE
    void set_merg_ctr();
    void get_merg_ori();
    void set_merg_ori();

    // Check if the surrounding elements become inverted after the collapse.
    // Print the inverting elements.
    bool chk_inversion(); //TODO update all merging vertices, surrounding entities joined together


    bool edge_comp(apf::MeshEntity* e1, apf::MeshEntity* e2, apf::Mesh2* m);
    void sort_edge(std::vector<apf::MeshEntity*>* edge_list);
    // Splitting extra merging edges, reloading, splitting merging topologies. 
    void extract_conn_map(cell_base* conn_map, 
                         std::vector<apf::MeshEntity*> ee,
                         std::vector<apf::MeshEntity*> es_in,
                         std::vector<apf::MeshEntity*> et);
    // Given a triangle and a list of intersection points lying on the triangle,
    // refine the triangle and the subsequent triangles recursively.
    void refine_tri(apf::MeshEntity* tri, std::vector<apf::Vector3> pts);

    bool corr_lens_new();
    void get_int_pos(std::vector<apf::MeshEntity*> ee, 
                           std::map<apf::MeshEntity*, apf::Vector3> &mer_pt,
                           std::map<apf::MeshEntity*, apf::Vector3> &int_pt,
                           std::map<apf::MeshEntity*, apf::Vector3> &vert_pt,
                           double r_cvx = -1);
    void get_int_pos2(std::vector<apf::MeshEntity*> ee, 
                           std::map<apf::MeshEntity*, apf::Vector3> &mer_pt,
                           std::map<apf::MeshEntity*, apf::Vector3> &int_pt,
                           std::map<apf::MeshEntity*, apf::Vector3> &vert_pt,
                           double r_cvx = -1);

    void split_ctr(std::vector<apf::MeshEntity*> &ee,
                   std::vector<apf::MeshEntity*> &es_in,
                   std::vector<apf::MeshEntity*> &et,
              std::vector<apf::MeshEntity*> &v_sp_set, 
              std::map<apf::MeshEntity*, apf::Vector3> &int_pt, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_sp_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_mer_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_ori_dir,
              double r_edge_min = -1);

    void split_cvx(std::vector<apf::MeshEntity*> &ee,
                   std::vector<apf::MeshEntity*> &es_in,
                   std::vector<apf::MeshEntity*> &et,
              std::vector<apf::MeshEntity*> &v_sp_set, 
              std::map<apf::MeshEntity*, apf::Vector3> &int_pt, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_sp_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_mer_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_ori_dir);


    void split_noncvx(std::vector<apf::MeshEntity*> &ee,
                   std::vector<apf::MeshEntity*> &es_in,
                   std::vector<apf::MeshEntity*> &et,
              std::vector<apf::MeshEntity*> &v_sp_set, 
              std::map<apf::MeshEntity*, apf::Vector3> &int_pt, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_sp_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_mer_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_ori_dir, double r_cvx);
    void split_noncvx2(std::vector<apf::MeshEntity*> &ee,
                   std::vector<apf::MeshEntity*> &es_in,
                   std::vector<apf::MeshEntity*> &et,
              std::vector<apf::MeshEntity*> &v_sp_set, 
              std::map<apf::MeshEntity*, apf::Vector3> &int_pt, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_sp_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_mer_pos, 
              std::map<apf::MeshEntity*, apf::Vector3> &v_ori_dir, double r_cvx);

    // The rate of change of volume wrt. alfa. If v is specified, only consider
    // the residual at that vertex.
    double calc_C(std::vector<apf::MeshEntity*> &v, 
      std::map<apf::MeshEntity*, apf::Vector3> &pos, 
                        int ii, apf::MeshEntity* v_r = NULL);
    // Calculate the terms in the softmax function. B is the 1/6 of the sum of 
    // the inward pointing normals of the triangles adjacent to vertex v.
    apf::Vector3 calc_B(std::vector<apf::MeshEntity*> &v_set, 
        apf::MeshEntity* v, std::map<apf::MeshEntity*, apf::Vector3> &pos);
    bool relax_inv_grad(std::vector<apf::MeshEntity*> &v_sp_set,
                              double r_edge_min);
    // Adaptation of Jeremy's implementation of conjugate gradient with inverse 
    // quadratic approximation of the function in the direction of the gradient 
    // for less evaluation of the gradient.
    bool chk_aperture(std::vector<apf::MeshEntity*> &v_sp_set);
    bool relax_inv_conjgrad(std::vector<apf::MeshEntity*> &v_sp_set,
                              double r_edge_min);
    double calc_energy(std::vector<apf::MeshEntity*> &v_sp_set,
                        std::vector<apf::MeshEntity*> &et,
                        std::vector<double> &C,
                  std::vector<std::vector<apf::MeshEntity*> > &tet_v,
                  std::vector<double> &exp_pw_tot,
                  std::map<apf::MeshEntity*, apf::Vector3> &v_pos,
                  std::map<apf::MeshEntity*, bool> &tet_skip,
                  std::map<apf::MeshEntity*, bool> &v_skip,
                  std::map<apf::MeshEntity*, int> &tet_i1,
                  double vol_scale, double w);

    void upd_grad(std::vector<apf::MeshEntity*> &v_sp_set,
                        std::vector<apf::MeshEntity*> &et,
                        std::vector<double> &A,
                        std::vector<std::vector<apf::Vector3> > &B, 
                        std::vector<double> &C,
                        std::map<apf::MeshEntity*, apf::Vector3> &pos_ori,
                        std::map<apf::MeshEntity*, apf::Vector3> &v_pos,
                        std::map<apf::MeshEntity*, apf::Vector3> &v_disp,
                        std::map<apf::MeshEntity*, bool> &tet_skip,
                        std::map<apf::MeshEntity*, bool> &v_skip,
                  std::vector<std::vector<apf::MeshEntity*> > &v_tet,
                  std::vector<std::vector<apf::MeshEntity*> > &tet_v,
                  std::vector<std::vector<int> > &tet_v_id,
                  std::vector<double> &vol_t,
                  std::vector<apf::Vector3> &g0,
                  std::vector<double> &exp_pw_tot,
                  double &v_th, double &vol_scale, double &f_norm, double &w);

    void shift_v_sp_pos(std::vector<apf::MeshEntity*> &v_sp_set,
                    std::vector<apf::Vector3> &x0,
                    std::vector<apf::Vector3> &ndir,
                    double xb);


    void relax_inv_dir(std::vector<apf::MeshEntity*> &et,  
 std::vector<apf::MeshEntity*> &es_in, std::vector<apf::MeshEntity*> &ee,
 std::vector<apf::MeshEntity*> &v_sp_set,
 std::map<apf::MeshEntity*, apf::Vector3> &v_sp_pos,
 std::map<apf::MeshEntity*, apf::Vector3> &v_ori_dir,
 std::map<apf::MeshEntity*, apf::Vector3> &v_mer_pos, double r_edge_min);

    void relax_inv(std::vector<apf::MeshEntity*> &et, 
      std::vector<apf::MeshEntity*> &es_in, std::vector<apf::MeshEntity*> &ee, 
 std::vector<apf::MeshEntity*> &v_sp_set,
 std::map<apf::MeshEntity*, apf::Vector3> &v_sp_pos,
 std::map<apf::MeshEntity*, apf::Vector3> &v_ori_dir,
 std::map<apf::MeshEntity*, apf::Vector3> &v_mer_pos, double r_edge_min);


    bool update_vol_inv_grad(std::vector<apf::MeshEntity*> &et,
                                   std::vector<double> &vol_t,
                                   std::map<apf::MeshEntity*, bool> &tet_inv,
                                   std::map<apf::MeshEntity*, bool> &tet_surr,
                                   std::map<apf::MeshEntity*, bool> &tet_skip,
                                   std::map<apf::MeshEntity*, int> &tet_i1,
                           std::vector<std::vector<apf::MeshEntity*> > &tet_v,
                   double &vol_neg, double &vol_min, double min_vol);

    bool corr_lens_new2();

    // Copy constructor
    vd_glens& operator=( const vd_glens& other );
    vd_glens(const vd_glens& that);

  public:


    // Check if all surrounding tris are non-merging during collapse. If not, 
    // preconditioning is necessary.
    bool chk_merg();
    // Check if all surrounding tets are non-inverting during collapse. If so, no 
    // preconditioning necessary.
    bool chk_precond(int dim, int cell_id);

    // Default constructor
    vd_glens(apf::Mesh2* msh, struct cell_base* c, vd_entlist* el_in);
    void l_clear();

    // Managae VTK and smb file output
    void set_files(const char* vtkFile, const char* meshFile);
    // TODO, internal lens should have different save files

    // Setup the collapse object. Return the cell dim and id of the new vertex.
    std::pair<int, int> set_cell(int dim, int cell);

    void evolve_boundaries();

    // Collapse the cell. Return the cell dim and id of the new vertex.
    std::pair<int, int> col_cell(int dim, int cell, int tag_0c_in = -1);

    // Return the cell dim and id of the new vertex.
    std::pair<int, int> get_cell();

    //TODO modify this, shrink, call vd_lens, collapse last entity, update topo
    void save_glens_vtk(const char* vtk_filename);

    void turn_save_vtk(bool onoff);
    void turn_save_vtk_sub(bool onoff);
    void set_verify(bool onoff);
    void set_precond(bool onoff);

    // Calculate if the given vertex is within the lens.
    apf::MeshEntity* get_vert_ctr();

    void createfields();
    void destroyfields();

    void set_field_calc(field_calc* f_calc_in, bool on_off = true);

    // Whether to check for inversion before and after collapses. Default is 
    // on, but in the case of collapsing unsuccessful insertions, it is turned
    // off.
    void set_inv(bool on_off);

    // Set whether to calculate split vertex positions by conjugate gradient or
    // using input positions.
    //void calc_sp_pos(bool on_off);

    // Test whether all edges can be intersected by the convex hull generated.
    bool test_cvx(int dim, int cell_id);

    // Set the simulation parameters.
    void set_vdpar(vd_param par_in);

    // Set the c_base modification flag.
    void set_mod_c_base(bool on_off);

    //bool tag_flag;
    //void vd_tag_glens();
    //void vd_start_tag();
    //bool vd_tag_ent_glens(apf::MeshEntity* ent, int tag);

    // Destructor:
    ~vd_glens();
};


// Given a cell to collapse, create a trial mesh with entitities adjacent to the
// bounding vertices of the cell. Check if preconditioning yields a collapseable
// configuration. If so, precondition the actual mesh and copy the positions.
// Assert it can be collapsed. 
class vd_glens_trial{
  private:
    apf::Mesh2* m;
    apf::Mesh2* m_trial;

    struct cell_base* c_base;
    field_calc* f_calc;

    vd_entlist* e_list;
    bool trial_load;
    bool precond_flag;

    bool vrfy_msh;

    int c_id;
    int c_dim;

    std::vector<apf::MeshEntity*> es;

    std::vector<std::vector<apf::MeshEntity*> > es_ent;

    std::vector<apf::MeshEntity*> vert_t;

    apf::MeshEntity* copyElement_pre(apf::MeshEntity* original,
                                            apf::BuildCallback* cb);
    void copy_ent();
    void set_cell(int cell_dim, int cell_id);

    void clear();
    void clear_ent();

  public:
    vd_glens_trial(apf::Mesh2* msh, struct cell_base* c, vd_entlist* e_list, field_calc* f_calc_in);

    //void set_field_calc(field_calc f_calc_in);

    std::pair<int, int> coll_cell(int cell_dim, int cell_id, int tag_0c_in = -1);
    void set_precond(bool on_off);

    void set_verify(bool onoff);
    // Destructor:
    ~vd_glens_trial();
};

// Given a cell to collapse, create a trial mesh with entitities adjacent to the
// bounding vertices of the cell. Check if preconditioning yields a collapseable
// configuration. Next check a successful insertion can be done. Check if the 
// inserted stratum has a different configuration.
class vd_glens_trial_ins{
  private:
    apf::Mesh2* m;
    apf::Mesh2* m_trial;

    struct cell_base* c_base;
    field_calc* f_calc;

    vd_entlist* e_list;
    bool trial_load;
    bool precond_flag;

    int c_id;
    int c_dim;

    std::vector<apf::MeshEntity*> es;

    std::vector<std::vector<apf::MeshEntity*> > es_ent;

    std::vector<apf::MeshEntity*> vert_t;

    apf::MeshEntity* copyElement_pre(apf::MeshEntity* original,
                                            apf::BuildCallback* cb);
    void copy_ent();
    void set_cell(int cell_dim, int cell_id);

    void clear();
    void clear_ent();

  public:
    vd_glens_trial_ins(apf::Mesh2* msh, struct cell_base* c, vd_entlist* e_list, field_calc* f_calc_in);

    //void set_field_calc(field_calc f_calc_in);

    std::pair<int, int> coll_cell(int cell_dim, int cell_id, int tag_0c_in = -1);
    void set_precond(bool on_off);

    // Destructor:
    ~vd_glens_trial_ins();
};


#endif
