#ifndef TOPO_EDISC_H
#define TOPO_EDISC_H

#include <vector>
#include <deque>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_extinfo.h"
#include "topo_topo.h"

#include "topo_graph.h"

#include "topo_lens.h"

#include "topo_entlist.h"

#include "topo_disp.h"

#include "topo_tmesh.h"

#include "topo_ma.h"

enum class WG_TYPE {
  TRI,
  EDGE,
  SPHERE,
  EVOLVE,
  END
};

// Expand disc object is used to make changes on the mesh compliant to the 
// addition of 1cell and 2cell topologies. 
// This is achieved:
// Selecting a 0cell.
// - Extraction of the 2-3cell adjacency map around the 0cell. [topo_graph]
// - Finding the circuits that split the graph into two disjoint graphs, each
//   containing at least a 2cell(3cell). [topo_graph] 
// - Selection of topological insertion (1-, 2cell) [could be implemented here]
// - Generation of a trial mesh, containing the entities surrounding the 
//   vertex of the 0cell [implement here]
// - Calculation of the energy change rate ((A_i-A_0)/sum(l_i)) 
//   [implement here]

// A disc is an object containing an ordered set of edges and another of 
// triangles, that can be used to construct the disc in a sequential manner.
// It also contains information whether it is ringlike, i.e. whether the last
// and the first triangle are connected.

// When expanding the disc, the cell memberships are designated as following:
// The disc expands, adding a new disc and a set of tets to fill between the 
// discs. 
// The new tets are connected by a set of triangles that will have their
// Cell memberships obtained by the cell memberships of the edges in the disc.
// The way this is obtained is related to how the disc is obtained:
//   - 3cell bounding triangles(1cell generation, preconditioning),
//   - 2cell bounding triangles(1cell generation),
//   - Path bounding triangles.

// The new disc triangles will belong to the cells of the new tets.
// The new disc edges will belong to the cells of the triangles bounding the 
// tets.

// ----------------------------------------------------------
// Operations used in finding triangles of discs.
// ----------------------------------------------------------

// Given a mesh, a center vertex and two edges that are bounded by this vertex,
// find the shortest path of triangles passing through a set of d-cells. Assumes
// the edges to be connected over entities belonging to that those d-cells.
// Append the triangles to the es_tri list.
bool vd_collect_tri_2cells(apf::Mesh2* m, apf::MeshEntity* v_ctr, 
                                   apf::MeshEntity* e_1, apf::MeshEntity* e_2,
                                   int cell_d, ent_conn* e_cells,
                                   std::vector<apf::MeshEntity*>* es_tri);

struct vd_disc {
  // The mesh of the disc, which can be used to find out about cell ids.
  apf::Mesh2* m;

  // The center vertex of the disc.
  apf::MeshEntity* v_ctr;
  // The motion directions of the new vertices when the disc expands.
  //apf::Vector3 vt;
  // The disc triangles.
  std::vector<apf::MeshEntity*> tri;
  // The disc edges.
  std::vector<apf::MeshEntity*> edge;

  std::vector<apf::Vector3> v_pos;

  // This stores the elements, neigboring the disc triangles. The top elements
  // are the elements that the disc will expand into. These can be used in the
  // burning algorithm.
  //std::vector<apf::MeshEntity*> elem_top;
  //std::vector<apf::MeshEntity*> elem_bottom;

  // Cell memberships of the disc edges.
  bool broken;
  apf::MeshEntity* e_end1;
  apf::MeshEntity* e_end2;

  // Cell memberships of the new disc triangles and tets.
  std::vector<apf::ModelEntity*> t_em;

  // Cell memberships of the edges and the triangles bounding the new tets.
  std::vector<apf::ModelEntity*> e_em;

  std::vector<apf::ModelEntity*> t2_em;

  std::map<apf::MeshEntity*, int> e_t_map1;
  std::map<apf::MeshEntity*, int> e_t_map2;

  std::map<apf::MeshEntity*, int> t_e_map1;
  std::map<apf::MeshEntity*, int> t_e_map2;

  vd_disc& operator=( const vd_disc& obj );

  vd_disc ();
  // Print the content of the disc.
  void print_cont();

  // Reset the disc.
  void clear();

  // Destructor.
  ~vd_disc();

};
void cp_tri(const std::vector<apf::MeshEntity*>* es_tri, vd_disc* disc);

// This is an extension of the disc, used to obtain the fins necessary for 
// 2cell insertion. Obtaining of the fin entities require a crawling in any 
// case, so this object readily implements it.
struct vd_fin : vd_disc {

  std::vector<int>* path;
  apf::MeshEntity* tri_first;
  apf::MeshEntity* edge_first;

  vd_fin (apf::Mesh2* m_in, std::vector<int>* path_in, 
            apf::MeshEntity* tri_in, apf::MeshEntity* edge_in);
};

// vd_elens is a container for discs and mesh entities used in lens expansion.
struct vd_elens {
  apf::Mesh2* m;
  std::vector<vd_disc> discs;

  bool init_flag;
  // Check if the map is nullified.
  bool null_flag;

  std::vector<apf::MeshEntity*> es_vert;
  std::vector<apf::MeshEntity*> es_edge;
  std::vector<apf::MeshEntity*> es_surf;
  std::vector<apf::MeshEntity*> es_elem;

  // Link entities in the entity sets with the slices associated with each disc.
  // The disc triangle and edge entities are already stored there.
  // -i-1 denotes disc triangles and edges. The slice entities are denoted by 
  // slice id+1, i.e. 1,2,...,n, where n is the number of slices.
  apf::MeshEntity* vert_ctr;
  apf::Vector3 pos;

  // Container for the new vertices, obtained by splitting the 0cell.
  std::vector<apf::MeshEntity*> vert_ctr_new;
  // Container for the new edges, obtained by splitting the 0cell. 
  // The inside of the new discs is filled with two lenses. The center edges
  // are seperated into the top and bottom lists.
  // The top and bottom center edges are connected by the new split vertices.
  std::vector<apf::MeshEntity*> edge_ctr_new_t;
  std::vector<apf::MeshEntity*> edge_ctr_new_b;
  std::vector<apf::MeshEntity*> vert_sp_new;

  // Used in 2cell insertion. 
  std::vector<std::vector<int > > c1;
  // Is a c1 repeated in multiple slices. Possible if c1 is bounded by a single
  // c0, in which case after the insertion there will be two c0 adjacencies of c1.
  std::vector<int> c1_mult;
  // Keep the center edges of each lens.
  std::vector<std::vector<apf::MeshEntity*> > edge_sp_new;
  // Keep the inside triangles of each lens. 
  std::vector<std::vector<apf::MeshEntity*> > tri_disc_new;

  // Stores the start and end edges of the 3-cell couple. If the start and end
  // 3-cells are both interior, the associated flags are true.
  int c3_f;
  int c3_e;

  bool e_3c_t_flag;
  bool e_3c_b_flag;
  apf::MeshEntity* e_3c_t;
  apf::MeshEntity* e_3c_b;

  apf::MeshEntity* v_2c;

  std::vector<int> new_cell0_id;
  std::vector<int> new_cell1_id;
  int new_cell2_id;

  cell_base* c_base_curr;

  // The motion directions of the new vertices when the disc expands.
  std::vector<apf::Vector3> vt;
  // The motion directions of the new split vertices, used in 2cell insertion.
  std::vector<apf::Vector3> split_vt;

  //std::vector<apf::Vector3> prev_mt_ctr;
  //std::vector<apf::Vector3> prev_mt_sp;

  // Each disc is associated with two slices and each slice is at most associated
  // with two paths.
  std::vector<int > slices;
  std::map<apf::MeshEntity*, int > slice_map;

  // Create map is a map between the old and the new entities of each disc.
  // The maps come in pairs, so create_map has size twice the nbr_of_disc.
  std::vector<std::map<apf::MeshEntity*, apf::MeshEntity* > > create_map;
  std::vector<std::map<apf::MeshEntity*, apf::MeshEntity*> > tri_tet_map;

  // Tri map is a map between the old edges of the disc and the new internal
  // triangles of the lens.
  std::vector<std::map<apf::MeshEntity*, apf::MeshEntity* > > tri_map;
  std::vector<std::map<apf::MeshEntity*, apf::MeshEntity* > > split_map;

  // Any entity, that is not associated with a disc, is stored in this map.
  std::map<apf::MeshEntity*, apf::MeshEntity* > create_map_rest;

  vd_elens();
  vd_elens(apf::Mesh2* m_in);

  // Change the mesh.
  void set_m(apf::Mesh2* m_in);

  // ----------------------------------------------------------
  // Operations used in finding triangles of discs.
  // ----------------------------------------------------------

  // Given a mesh, a center vertex and two edges that are bounded by this 
  // vertex, find the shortest path of triangles passing through a single 
  // 3cell. Assumes the edges to be connected over entities belonging to that 
  // 3cell. Do not cross edges on any disc.
  // Append the triangles to the es_tri list.
  bool vd_mesh_find_short(apf::MeshEntity* e_1, apf::MeshEntity* e_2,
                          std::vector<apf::MeshEntity*>* es_tri, int c3_id = 0);

  // Check the slice tetrahedra for negative volume entities after the 
  // insertion.
  long int vd_chk_neg_slice();

  // Print the content of the lens.
  void print_cont();

  void clear_old_ent();
  void clear_discs();
  void clear_new_vert();
  void clear_create_maps();

  // Reset the lens.
  void clear();
  void null_slice_map();

  ~vd_elens();
};

// Used in slicing the new discs.
class vd_disc_cut {
  private:
  public:
    // Intersection tolerance:
    double int_tol;
    circuit* circ_tup;
    std::pair<std::vector<int>, std::vector<int> >* ng_pt;
    std::vector<apf::MeshEntity* > c2_edge_new;
    std::vector<apf::ModelEntity* > c3_mdl;
    apf::ModelEntity* mdl_curr;
    std::vector<bool> c2_edge_flag;
    //std::vector<apf::MeshEntity* > c3_edge_new;

    std::vector<int> c2_cuts;
    std::vector<apf::MeshEntity*> c2_edges;
    std::vector<apf::Vector3> c2_cut_pos;
    std::vector<apf::MeshEntity*> ents;

    std::map<apf::MeshEntity*, bool> ent_burn;

    void clear();
    vd_disc_cut();
    ~vd_disc_cut();
};

class vd_edisc {
  private:

    // ---------------------------------------------------
    // Pointers to the mesh and cellbase objects:
    //struct cell_ins_chk* c_ins;

    struct cell_base* c_base;
    struct cell_base* c_base_init;
    struct cell_base* c_base_act;

    apf::Mesh2* m_main;

    field_calc f_calc;
    vd_entlist* e_list;
    vd_entlist* e_list_act;

    vd_cell_det* vd_cd;
    vd_3c_det* vd_3c;

    // Insert without comparing energies.
    bool skip_en;
    // TODO Consider spurious entities. Obsolete and kept as a placeholder.
    bool spur;

    // Exterior bounds the candidate 0cell.
    bool calc_ext;
    bool calc_corner;

    //ext_shell e_sh_save;
    shell sh_old;
    bool sh_ext;

    bool ext_0cell;
    bool cell_flag;

    // Used in exterior 2cell insertions.
    // Contains the shell membership of exterior 2cells. 
    std::vector<int> shell2_list;
    // Shells already connecting the exterior and the interior 3cell. A new 
    // insertion is only possible if there is an available shell. If more than
    // one exists, pick the one 
    std::map<int, bool> shell2_map;
    std::map<apf::MeshEntity*, bool > wing_map;

    // Exterior vertex.
    std::map<apf::MeshEntity*, bool> ext_new;
    // Slice contains the exterior.
    std::map<int, bool> ext_slice;
    // Path or circuit contains the exterior.
    std::map<int, bool> ext_pc;
    std::map<int, bool> ext_cor;
    // Used to preserve volume during insertions.
    std::map<int, int> ext_proj_type;
    std::vector<std::vector<apf::Vector3> > ext_proj_dir;

    // Fixes the velocities of the newly created vertices. Used in the first 
    // expansion after new entity generation, where some tets are zero volume. 
    // Currently assumes the boundaries to be flat.
    void clear_vel_map();
    void coll_proj_map();

    apf::Vector3 fix_vec_ctr_new(int i, apf::Vector3 vec_in);
    apf::Vector3 fix_vec_sp_new(int i, apf::Vector3 vec_in);

    apf::Vector3 calc_v_avg();
    apf::Vector3 calc_v_avg(std::map<apf::MeshEntity*, bool> sk_map);
    void fix_vel_new();

    // Record each step of the relaxation iteration.
    bool mov_flag;
    // Verbosity flag.
    bool verb_flag;

    bool isotropic;
    bool drag_flag;

    // Flags checking the step of cell insertion algorithm.
    // Mesh exists:
    bool m_ex;
    bool precond_ex;
    bool trial_ex;

    // Mesh loaded, or has entities:
    bool m_load;
    bool precond_load;
    bool trial_load;

    // Mainly used for tracking cell insertions and naming the saved vtk files.
    int trial_type;
    int trial_curr;
    int ngon_curr;

    // Path and circuit discs are loaded.
    bool disc_load;
    // Precondition successful.
    bool pre_flag;
    // Cell insertion successful.
    bool ins_flag;

    // Store the 1cell circuit number associated with minimum energy 
    // configuration. Also store the index of the 2cell path set associated 
    // with the minimum energy 2cell insertion. Compare these to apply an 
    // operation.
    int circ_min;
    int path_min;
    int ngon_min;

    int circ_sz;
    int pt_sz;

    apf::Vector3 pos_old;
    apf::Vector3 v_avg;

    double vel_max;

    // The cos(radians) threshold for expansion direction for the lenses. Used in
    // expansion inversion calculations.
    double exp_dir_cos_th;

    double en_min_1c;
    double en_min_2c;
    // A workaround for insertions around the 0strata where the energy minimizing
    // insertion is a spurious insertion and a similar insertion that is 
    // physically not likely can get inserted, which collapses in coming time 
    // iterations in ways that can slow down the simulation.
    // In this case, an insertion should only be allowed if the energy dissipation
    // rate is lower than the highest dissipation rate spurious insertion.
    // energy_th tracks the rate of that spurious insertion.
    // If no insertion can be found allow the simulation to relax.
    double energy_th;
    bool sh_used;

    double dt_curr;

    // Time of insertion
    double t_total;

    // Threshold for time step.
    double dt_max;

    // The numbers in circuits and paths refer to the indices of the 3cells and 
    // 2cells in these lists. The 2cells follow the 3cells, so the first 2cell
    // has the index cells3.n_fill in the path and circuit lists.
    struct ent_conn cells3;
    struct ent_conn cells2;

    // Container for lens extension.
    vd_elens e_lens;

    // Container for ngons in 2cell insertion.
    ngon_gmi ng;

    // At least a preconditioning is usually necessary before the 2cell 
    // insertion. A preconditioned copy of the modified entities can be kept for
    // easier resetting of the mesh between the trials.

    apf::Mesh2* m_precond;
    apf::Mesh2* m_trial;
    // Pointer to the active mesh(Actual, preconditioned or trial).
    apf::Mesh2* m_act;

    // Mesh center vertex, cell membership and id, copied center vertex.
    apf::MeshEntity* vert_ctr;
    //apf::MeshEntity* vert_ctr_pre;
    apf::ModelEntity* vert_ctr_em;
    int cell_id;

    // Length of the shortest edge. TODO Used to limit the displacements to 
    // apply during cell insertions.
    double len_sh;
    double len_sp;
    double l_min;
    // Ratios controlling insertion size compared to len_sp. To be controlled
    // at interface level to be consistent with collapse size.
    double rho_rat;
    double fudge_factor;
    double rat_init;

    double len_edge;

    // Controls whether sub vtk files are generated or not.
    bool sub_vtk;

    // The trial and preconditioned mesh center vertices.
    apf::MeshEntity* vert_ctr_cp;
    apf::MeshEntity* vert_ctr_precond;

    apf::MeshEntity* vert_ctr_act;

    // PATH and CIRCUIT information.
    // Containers for edges passing through disjoint set of triangles of 
    // 2cells.
    std::vector<apf::MeshEntity*> c2_edge;
    std::vector<apf::MeshEntity*> c3_edge;

    // Containers to store the elements that will be adjacent to the new 
    // vertices. Populated by a burning algorithm, by the elements of the trial
    // mesh. 
    // Tri burn is a reference for used triangles, including the fin triangles
    // bounding the paths.
    std::vector <std::vector < apf::MeshEntity* > > elem_slice;
    std::map<apf::MeshEntity*, bool> elem_burn;
    std::map<apf::MeshEntity*, bool> tri_burn;

    // The currently active elements.
    std::deque<apf::MeshEntity*> elem_ofire;

    // Called before expanding the lenses.
    void crt_new_vert_list();
    void crt_new_vert();
    void crt_new_cell();
    void asgn_new_vert_ext();

    // Called inside cell insertion routines. 

    // Collect the elements around the center vertex.
    void get_ent_set_new();
    // Subrouting for burning current slice of elements.
    bool burn_slice(int slice);

    // The locations of the discs and fins are arbitrary the resulting slices
    // can contain arbitrarily shaped triangles belonging to path and circuit 
    // 2cells. For consistency, do not consider these triangles in velocity 
    // calculation. 
    void burn_slice_wings(int disc_curr);

    // If it is a 1cell insertion, collect the two sets of elements, found on
    // different sides of the disc.
    // If it is a 2cell insertion, collect the sets of elements, belonging to 
    // the slices seperated by the triangles fins, bounding the used paths.
    void burn_trial();

    // Collect triangles of a the slice, belonging to a 2cell.
    void collect_slice_tris(
                        std::vector<std::vector<apf::MeshEntity*> > * es_tri);
    void collect_slice_tri(std::vector<apf::MeshEntity*>* es_tri, int slice);
    void collect_path_tri(std::vector<apf::MeshEntity*>* es_tri, int path_id);
    void collect_disc_path(int ngon_id);


    // The 2cell triangles belonging to each slice and adjacent paths.
    std::vector< std::vector<apf::MeshEntity* > > slice_tris;
    std::vector< std::vector<apf::MeshEntity* > > slice_edges;

    // In order to keep the path entities consistent, the graph node ids of the
    // disjoint entity sets should be fixed. This is done for the 
    // preconditioned mesh, which is translated to the trial and main meshes
    // when needed.
    std::map<apf::MeshEntity*, int> preid_map;
    std::map<apf::MeshEntity*, int> actid_map;
    std::map<apf::MeshEntity*, apf::MeshEntity*> main2pre_map;
    std::map<apf::MeshEntity*, apf::MeshEntity*> trial2pre_map;

    std::map<apf::MeshEntity*, apf::MeshEntity*> pre2main_map;
    std::map<apf::MeshEntity*, apf::MeshEntity*> pre2trial_map;

    void precond_map();
    void trans_map_pre2act();

    // Expand the discs into lenses and insert a set of triangles for the 2cell 
    // and tetrahedra for the 3-cells to be connected. 
    void expand_disc_path();
    //void collect_slice_path(std::vector<apf::MeshEntity*>* es_tri, int slice);

    void color_disc(int i);
    // Color the slice map of the ith slice.
    void color_slice(int i);

    // Entity destruction lists.

    // The length of the inserted edges. Used in rate of change calculation.
    // double ins_length;

    // ---------------------------------------------------
    // 

    // Get all surrounding entities:
    std::vector<apf::MeshEntity*> es_vert;
    std::vector<apf::MeshEntity*> es_edge;
    std::vector<apf::MeshEntity*> es_surf;
    std::vector<apf::MeshEntity*> es_elem;

    std::vector<apf::MeshEntity*> es_vert_act;
    std::vector<apf::MeshEntity*> es_edge_act;
    std::vector<apf::MeshEntity*> es_surf_act;
    std::vector<apf::MeshEntity*> es_elem_act;

    std::vector<apf::MeshEntity* > vert;
    std::vector<apf::MeshEntity* > vert_trial;

    // ---------------------------------------------------
    // The relevant mesh entities, 2cell insertion.
    // The central edges in the preconditioned mesh in the start and end 3cells.
    apf::MeshEntity* e_s_pre;
    apf::MeshEntity* e_e_pre;

    apf::MeshEntity* e_s_tri;
    apf::MeshEntity* e_e_tri;

    // The triangles around the e_s_tri. 
    //std::vector<vd_disc> disc_curr;

    // The triangles around the e_s_tri. 
    //apf::Up tri_start;

    // ---------------------------------------------------
    // Functions and variables used in energy change calculation. 
    // Energy before.
    double energy;
    double energy_bef;

    // Energy change associated with the 1cell and 2cell insertions.
    std::vector<double> w1;
    std::vector<std::vector<double> > w2;
    std::vector<double> e1;
    std::vector<std::vector<double> > e2;
    std::vector<std::vector<double> > w1_exp;
    std::vector<std::vector<std::vector<double> > > w2_exp;

    // Different physics can be implemented here to interact with fields 
    // attached to the mesh.
    //double (*energy_func_tri)(apf::Mesh2*, apf::MeshEntity*);
    //double (*energy_func_tet)(apf::Mesh2*, apf::MeshEntity*);
    double energy_func_tri(apf::Mesh2*, apf::MeshEntity*);
    double energy_func_tet(apf::Mesh2*, apf::MeshEntity*);

    // Copy entities from the original mesh.
    apf::MeshEntity* copyElement_pre(apf::MeshEntity* original,
                              apf::BuildCallback* cb);
    apf::MeshEntity* copyElement(apf::MeshEntity* original,
                              apf::BuildCallback* cb);

    // The trial mesh copier for dissipation rate calculation.
    apf::MeshEntity* copyElement_rpl(
          apf::MeshEntity* original,
          apf::MeshEntity* v_orig,
          apf::MeshEntity* v_repl,
          std::map<apf::MeshEntity*,apf::MeshEntity*> &p2t,
          std::map<apf::MeshEntity*,apf::MeshEntity*> &t2p,
          apf::BuildCallback* cb);

    // ---------------------------------------------------
    // Operations used in setting up the trials.

    // Set the vertex belonging to the selected 0cell. Calls set_trial, 
    // which preconditions.
    void set_vertex(apf::MeshEntity* vertex);

    // Set trial mesh, initial configuration.
    void set_trial();

    // ---------------------------------------------------
    // Operations used in insertion.

    bool edge_comp(apf::MeshEntity* e1, apf::MeshEntity* e2, apf::Mesh2* m);
    void sort_edge(std::vector<apf::MeshEntity*>* edge_list);
    //void sort_edge(std::vector<apf::MeshEntity*>* edge_list, 
    //                                             int left=-1, int right=-1);

    int partition_edge_list(std::vector<apf::MeshEntity*>* edge_list, int left, int right);

    // Check both disjoint graphs, if both of them contain a 3cell, return 0,  
    // if one of them contains a 3cell, return the index+1, else return -1.
    int cont_3cell_min (s_graph* circ_graph);

    // Initialize elens with the current active mesh and vertex.
    void init_elens();
    void reset_elens_ent();

    // Collect the disc triangles for the given split operation.
    void collect_disc();


    // Check the given set of triangles for members that 
    // belong to a 3cell, 
    // bound the center vertex and
    // have edges adjacent to central vertex with lower dimensional cell 
    // membership. 
    bool chk_span_surf(std::vector<apf::MeshEntity*> &es_surf);

    // Given a set of tetrahedra belonging to a disjoint 3stratum, collect the 
    // bounding 2stratum triangles.
    std::vector<apf::MeshEntity*> coll_bound_tri(
                              std::vector<apf::MeshEntity*> &es_tet);
    // Collect the disc triangles and the slice entities.
    void color_discs();

    // After calculating the new vertex motions, go over the discs to find 
    // inverting elements.
    void assert_inv_pre();
    bool detect_inv();
    bool detect_inv(int disc_id, apf::Vector3 v_pos, apf::Vector3 mot_dir);

    // Expand the disc:
    void expand_lens();
    // Fill inside the lens:
    void fill_lens();
    // Recreate the element if it has less then negative double point accuracy
    // volume.
    apf::MeshEntity* recreate_inv_elem(apf::MeshEntity* tet, 
                                                    apf::ModelEntity* mdl);

    void fill_lens_elem();
    void fill_3c_lens_elem();

    // Check the new interior lens elements for left handed curls.
    bool chk_ma_swap_lens();

    // Create the 2cell triangles and fill the openings inside the 3cells 
    // getting connected.
    void fill_void();

    int get_circ_type(int circ_in);
    circuit* get_circ(int circ_in);
    void get_circ_topo_gmi(std::vector<int>* c2_list, 
                                std::vector<int>* c3_list, circuit* circ);
    // Entity, path and circuit collection operations.
    void get_path_ngon(int c3_cp, ngon_gmi* ng);

    void conv_path_ngon_gmi(ngon_gmi* ng);
    int conv_path_2c_gmi(int n_id);
    int conv_path_3c_gmi(int n_id);

    // Collect the affected 1cells.
    void collect_cell();

    // Collect the affected 1cells on slices, without relying on mesh information.
    void collect_cell_wg(std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells);

    void upd_cell();
    // Without modfying the mesh.
    void upd_cell_wg();
    void upd_cell2();

    // Recreate the slice entities:
    void recreate_slices();
    void destroy_ent();

    // Going over all slices of all possible insertions, find the dt such that
    // the maximum of vertex motion is less than the half of the edge lengths.
    double calc_dt();

    // Calculate the velocities of the new 0cell vertices using the Mason 
    // algorithm, using the slice entities only.
    void calc_ctr_vel(bool fix);
    void calc_sp_vel(bool fix);

    // Move the new vertices by using the velocities in the e_lens.
    void find_max_disp();
    void find_ctr_vel();
    void find_sp_vel();

    void evolve_ins();

    void move_ctr();
    void move_sp();


    // Collect the triangles associated with each center and split vertex, without
    // generating new starum entities.
    void collect_tris_wg(
                std::vector<std::pair<std::vector<int >, std::vector<int > > >*
                                                                  path_cells,
                std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells);
    void collect_edges_wg(
                std::vector<std::pair<std::vector<int >, std::vector<int > > >*
                                                                  path_cells,
                std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells);

    void update_circpath_energies(
           std::vector<std::pair<std::vector<int >, std::vector<int > > >*
                                                                  path_cells);

    bool chk_inv_wg();

    // Calculate the velocities without generating new vertices.
    void calc_ctr_vel_wg();
    void calc_sp_vel_wg();

    // Calculate the positions of the new 0cell vertices using the Mason 
    // algorithm, by iterating until the direction of motions of vertices 
    // converge.
    void fix_skip(std::vector<apf::MeshEntity*>* ent_col, 
                              std::map<apf::MeshEntity*, bool> skip_map);
    std::map<apf::MeshEntity*, bool> calc_skip();


    // Calculate the highest multiplier for velocity field that would invert any 
    // tetrahedra. Position the new vertices such that nothing inverts. Used 
    // before the lenses are fully expanded.
    void find_pos_non_inv_init();

    // Calculate the highest multiplier for vector field that would invert any 
    // tetrahedra. Position the new vertices such that nothing inverts.
    void find_pos_non_inv();

    double calc_dt_ext_dir(std::vector<apf::Vector3> &v_ctr, 
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
                         double &ang_th, double rat_t);
    void move_vert_ext(std::vector<apf::Vector3> &v_ctr, 
                               std::vector<apf::Vector3> &v_ctr_old, 
                               std::vector<apf::Vector3> &v_sp,
                               std::vector<apf::Vector3> &v_sp_old, 
                               double &d_max);
    void find_ext_dir();
    void find_ext_dir3();

    // Find the minimum and maximum distances from the 1stratum centers for 
    // vertices bounding 1strata and 2stratum center for vertices bounding the 
    // 2stratum.
    std::pair<double, double> vd_find_sph_bound();
    double get_r_max(std::vector<double> &r_ctr, std::vector<double> &r_sp);
    void proj_around_ctr(double rho);

    // Insert an edge inside the 3cell. Used in preconditioning step to treat
    // disconnected 3cell couples.
    void insert_3celledge(int cell3_curr);

    // Default constructor removed.
    vd_edisc();

    void clear_map();
    void clear_ext();
    void clear();
    // Clear the trial mesh by removing the attached fields and removing its 
    // entities starting from highest dimension. 
    void clear_trial();

  public:
    // TODO Used to determine which type of dissipation rate calculation is 
    // better. Will be discarded later on.
    WG_TYPE wg_tag;
    // Calculate the dissipation rates generating new vertices.
    void calc_max_diss_wg();
    void calc_max_diss_trial_wg();

    // Change the active mesh in a centralized manner.
    void act_precond();
    void act_trial();
    void act_main();

    // Constructor.
    vd_edisc(apf::Mesh2* m_in, struct cell_base* c, vd_entlist* el_in);
    // Destructor.
    ~vd_edisc();

    // ---------------------------------------------------
    // Try 1cell and 2cell insertions. If one minimum energy configuration is 
    // found apply to the actual mesh.
    std::pair<int, int> try_insert();

    // ---------------------------------------------------
    // 1cell insertion
    // Going over the available circuits, try the insertion and store the change
    // in energy.
    void try_1cell();

    // ---------------------------------------------------
    // 2cell insertion
    void try_2cell();

    // ---------------------------------------------------
    // 2cell insertion interface.
    void insert_2cell(int c3_id, int ngon_id);

    // ---------------------------------------------------
    // 1cell insertion with the current circuit.
    void insert_1cell(int circ);

    // Check the folding of the loaded lens. If there are creases on the disc
    // (the area normals of subsequent triangles change sign) such that the 
    // disc expansion leads to inverted elements, fix the disc and return true.
    //bool fix_disc_fold();

    // By collecting the edges in the disjoint graph, find the entities to  
    // seperate with the top vertex.
    void get_disc_edges(std::vector<ent_conn>* cs_in, ent_conn* c2_circ, 
      ent_conn* c3_circ, vd_disc* disc_in);
    void get_disc_edges(circuit* circ_in, vd_disc* disc_in);

    // Select the 0cell to be split. 
    bool set_0cell(int cell, bool spur_in = false);

    // A wrapper for ma curl check implementation to check mesh validity. 
    // Check edges around the newly created vertices.
    void chk_ma_new();

    // Precondition the mesh for 1cell and 2cell insertions.
    // Precondition the input mesh for 2cell insertion. Preconditioning  
    // involves making sure that there is an internal edge for each 2cell and  
    // 3cell. This is to be called once on the trial mesh to prepare for the 
    // trials and once for the actual mesh for splitting.
    void precond_mesh();
    void collect_split_tri(std::vector<apf::MeshEntity* >* tri_2_sp, 
                                         vd_entlist_v* ent_list);
    void collect_split(std::vector<apf::MeshEntity* >* e_2_sp, 
                                         vd_entlist_v* ent_list);


    // The trial mesh generation for dissipation rate calculation.
    void copy_trial_wg(apf::MeshEntity* v_new, 
                                 std::vector<apf::MeshEntity*> &e_set, 
                              std::map<apf::MeshEntity*,apf::MeshEntity* > &p2t, 
                              std::map<apf::MeshEntity*,apf::MeshEntity* > &t2p);

    void create_trial_vert_wg(std::vector<apf::MeshEntity*> &e_set, 
                              std::map<apf::MeshEntity*,apf::MeshEntity* > &p2t, 
                              std::map<apf::MeshEntity*,apf::MeshEntity* > &t2p);

    void copy_slice_ents_wg(int slice, 
            std::vector<apf::MeshEntity*> & v_ctr_new,
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_ctr,
                std::map<apf::MeshEntity*, int> &s_map,
                std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells);
    void create_path_int_p_wg(int path, int slice,
                    std::vector<apf::MeshEntity*> & v_ctr_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_ctr,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map);
    void create_path_int_s_wg(int path, int slice,
                    std::vector<apf::MeshEntity*> & v_ctr_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_ctr,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map);

    // Create the outer triangles and the tetrahedra for path tetrahedra.
    void create_path_int_s_tet_wg(int path, apf::MeshEntity* e_curr,
                    apf::MeshEntity* tet, 
                    std::vector<apf::MeshEntity*> & v_ctr_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_ctr,
            std::map<apf::MeshEntity*, int> &s_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map,
                    std::vector<std::vector<apf::MeshEntity*> > &s_tets);
 
    apf::MeshEntity* create_path_int_p_tet_wg(int path, apf::MeshEntity* e_curr,
                    apf::MeshEntity* tet, 
                    std::vector<apf::MeshEntity*> & v_sp_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_sp,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_sp,
            std::map<apf::MeshEntity*, int> &s_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &p2sv_map);

    // Create the 3-cell tetrahedra for all spheres.
    // Using the vertex ordering of the tet on slices, replace the following:
    // the center vertex by the 2-cell vertex or it's copy
    // the vertex across triangle by the 0-cell vertex or it's copy
    // the 3-cell vertex by the copy of the 3-cell vertex
    // the last vertex by the 1-cell vertex 
    void create_3c_tet_wg(int path, apf::MeshEntity* e_curr,
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
                        std::vector<std::vector<apf::MeshEntity*> > &s_tets);

    // Creating of path related entities on path and slice spheres:
    void create_path_p_wg(int path, int tri_id, apf::MeshEntity* e_curr,
                    std::vector<apf::MeshEntity*> & v_sp_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_sp,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_sp,
            std::map<apf::MeshEntity*, int> &s_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &p2sv_map,
                    std::vector<std::vector<apf::MeshEntity*> > &p_tets_s1,
                    std::vector<std::vector<apf::MeshEntity*> > &p_tets_s2);


    void create_path_s_wg(int path, apf::MeshEntity* e_curr,
                    std::vector<apf::MeshEntity*> & v_ctr_new, 
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &p2t_ctr,  
            std::vector<std::map<apf::MeshEntity*,apf::MeshEntity* > > &t2p_ctr,
            std::map<apf::MeshEntity*, int> &s_map,
                    std::vector<std::map<int, apf::MeshEntity*> > &s2pv_map,
                    std::vector<std::vector<apf::MeshEntity*> > &s_tets);

    void copy_path_disc_wg(int path, 
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
                    std::vector<std::vector<apf::MeshEntity*> > &s_tets);

    void evolve_equi_wg(
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
                    apf::Vector3 &force_2c);

    //void evolve_equi_wg(apf::MeshEntity* v_ctr);

    // Update the positions on the spheres by updating the forces and velocities.
    void upd_pos_sph_wg(std::vector< std::pair< std::pair<int,int>, 
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
                    std::vector<std::vector<apf::MeshEntity*> > &s_tets);

    void expand_2c_void_wg(
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
                    std::vector<std::vector<apf::MeshEntity*> > &s_tets);

    void reload_trial_wg(
      std::vector<std::pair<std::vector<int >, std::vector<int > > >* path_cells,
                    std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells);

    // Destroy the trial mesh and reload.
    void reload_trial();

    void collect_pc();
    void overwrite_pc(vd_disc_cut* vd_cut);

    void refresh_e();

    bool burn_the_cut_c2(int k, vd_disc_cut* vd_cut, vd_plane* pl);
    apf::MeshEntity* burn_c3_edge(apf::MeshEntity* e_curr,
                                     vd_disc_cut* vd_cut, vd_plane* pl);
    apf::MeshEntity* burn_c3_tri(apf::MeshEntity* e_curr,
                                     vd_disc_cut* vd_cut, vd_plane* pl);

    std::pair<apf::MeshEntity*, apf::MeshEntity*> get_bound
                                          (apf::MeshEntity* edge_c2);
    bool burn_the_cut_c3(int k, vd_disc_cut* vd_cut);
    bool burn_the_cut_c3_old(int k, vd_disc_cut* vd_cut, vd_plane* pl);

    bool cut_sing_crease(vd_disc_cut* vd_cut);
    bool do_the_cuts(vd_disc_cut* vd_cut);
    void set_the_edges(vd_disc_cut* vd_cut);

    void reload_edges();
    void clear_edges();

    // Given a set of triangles, order the triangles and connecting edges
    // and store into a disc object.
    //bool ord_tri(const std::vector<apf::MeshEntity*>* es_tri, vd_disc* disc);
    //void cp_tri(const std::vector<apf::MeshEntity*>* es_tri, vd_disc* disc);

    // ---------------------------------------------------
    // Common insertion methods

    // Calculate the energy change. 
    // Having generated a trial mesh, calculate the rate of energy change.
    // Trial mesh may require addition of fields for energy related 
    // calculations. This should be governed by the vd_sim object.
    double calc_energy();

    double calc_energy_lens();
    double calc_perim_avg_per_side();


    // Given a 2stratum on one of the discs, find the vertex joining the edges on 
    // the lens belonging to the stratum. Move the vertex in a range without 
    // inverting the tets and calculate the power dissipation rate of insertion. 
    // This is to determine the effect of mesh anisotropy on the choice of 
    // insertion in energetically close insertions.
    void vary_vert_pos(int id);
    double calc_energy_diss_rate_sing();
    void calc_energy_diss_rate();

    double calc_force_lens(double rho);
    double calc_energy_vertices(apf::Mesh2* m_in,
                                     std::vector<apf::MeshEntity*>* v_in);


    // Given a set of tetrahedra belonging to a disjoint 3stratum, collect the 
    // bounding 2stratum triangles touching the central vertex. Calculate the 
    // anglecosine weighted sum of area normals. Also, return false if the 
    // anglecosine between the direction and one of plane normals is lower than a 
    // threshold.
    bool coll_bound_tri_dir(std::vector<apf::MeshEntity*> &es_tet,
                                      apf::Vector3 &dir, double ang_th = 0.999);

    // Set/get the allowed size of the insertions.
    void set_len_sh(double len_in);
    void set_rho_rat(double rat_in, double fudge_in);
    double get_len_sp();

    void set_field_calc(const field_calc& FC);

    // Apply the change on the actual mesh, by creating an entity with 
    // average entity width. (or slightly larger?).
    // This is called by cell_graph after loading the wanted circuit and 
    // operation. 
    void apply_change();

    // Set the flag for skipping energy dissipation rate calculation. Used in
    // applying a predetermined topological change.
    void set_en_skip(bool skip_in);

    // Set the simulation parameters.
    void set_vdpar(vd_param par_in);

    void set_isotropic(bool iso_in);

    // Set exterior options.
    void set_calc_ext(bool fix);
    void set_calc_corner(bool fix);

    void set_proj(PROJ_TYPE PROJ);

    // Visual output of the trial mesh. 
    void vtk_trial();
    void vtk_precond();
    void vtk_mesh();

    void set_sub_vtk_flag(bool vtk_flag);

    // Set whether to print terminal output related to substeps.
    void set_verbose_flag(bool flag_in);

    // Given a list of 2shells, find the joint shell.
    shell get_2shell_joint(std::vector<int>* shell_2_ids);

    std::vector<int> get_shell_id_list(int dim, std::vector<int>* tag);

    void assign_2sh(int sh_2id);
    void assign_01sh();
    // Based on the shell memberships, update the exterior labels of cells.
    void update_ext();
    void reset_ext();

    // Returns true if for the current 2-cell trial, there is a single spurious 
    // bounding 1-stratum. Also, if more than one bounding 1-cell belongs to the 
    // same 1shell, return true. Used in exterior insertions to skip migration 
    // only type insertions.
    bool chk_bound_spur();

    // The updated shell assignment. Considers unused shells and dissipation rates
    // using different shells when determining shell assignment.
    void determine_shell_gmi();
    // Without actually modifying the mesh:
    void determine_shell_gmi_wg();

    // Export the current dissipation rates to csvfile.
    void write_diss_csv(const char* outfile);
    void write_diss_exp_csv(const char* outfile);
    void write_ei_csv(const char* outfile);


    void update_shell_gmi();
    void update_shell_pos();

}; // class vd_edisc

// !!FUTURE!! ENERGY RELATED

// TODO implementation of euler characteristic. Implement isotropic case for 
// now.
/*
class vd_calc_euler {
  private:
  public:

};
*/

// TODO the following part could be customizable (virtual).
// This is a container for the energetically relevant entities and the fields
// to be considered in the calculation. 
/*
class vd_energetic_ent {
  private:
    std::vector<apf::MeshEntity*> entities;
    std::vector<apf::Field*> fields;

  public:

};
*/
#endif
