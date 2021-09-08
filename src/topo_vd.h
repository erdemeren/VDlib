#ifndef TOPO_VD_H
#define TOPO_VD_H

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include "ma.h"
#include <apfNumbering.h>
#include <apfShape.h>

#include <PCU.h>

#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_extinfo.h"

#include "topo_topo.h"

#include "topo_glens.h"

#include "topo_lens.h"

/* Topology manipulator:  */
#include "topo_manip.h"
#include "topo_disp.h"
/* Topology write:  */
#include "topo_write.h"
/* Topology geometry:  */
#include "topo_geom.h"
/* Topology field:  */
#include "topo_field.h"
#include "topo_ma.h"

#include "topo_edisc.h"

#include "topo_tess.h"

#include "topo_energy.h"

#include "topo_feat.h"

typedef std::pair<int, int> cell_2_col;

enum class T_COST_TYPE {
  T_COST_FLAG_ITER,
  T_COST_FLAG_ADAPT,
  T_COST_FLAG_INS,
  T_COST_FLAG_COL,
  T_COST_FLAG_AUX,
  END
};

// The flags used in defining the options in time input file passed to the 
// set_time_opts function.
// <TIME_OPTS_TYPE>::<option specifiers>;...
// TIME_FIXED; 0/1;
enum class TIME_OPTS_TYPE {
  TIME_START,      // Start time
  TIME_STEP,       // Time step
  TIME_FIXED,      // Fixed time step
  TIME_END,        // End time
  TIME_SUB_LIMIT,  // The minimum time step allowed before termination
  CONDITION_G_NBR, // End condition if specified based on remaining nbr of grains
  CONDITION_S_NBR, // End condition if specified based on remaining nbr of surfs
  CONDITION_L_NBR, // End condition if specified based on remaining nbr of lines
  END
};


// This contains the vertex dynamics simulation related parameters and 
// controls.
class vd_sim{
  private:
    field_calc f_calc;
    //vd_param vdparam;

    // Mesh and model objects and files.
    apf::Mesh2* m;

    // Mesh loaded flag.
    bool mesh_flag;
    bool load_flag;
    bool save_mov;

    // Loading a geometry and mesh for the first time. Used to distinguish 
    // between reloading just the mesh and loading a new model and mesh.
    bool first_time;

    bool ins_flag;
    bool col_flag;

    bool extract_flag;
    //std::vector<vd_ext_opts> m_opts;
    //std::string outputMS;
    //std::string outputROC;
    //std::string outputVOL;
    std::string ext_opts;

    // Flag for isotropic energy. In case of isotropic energy, insertion only
    // happens if modified euler characteristic is not satisfied.
    bool isotropic;

    // Flag to calculate exterior vertex motions.
    bool calc_ext;
    bool calc_corner;

    int proj_flag;
    INTEG_TYPE integ_flag;
    // The exterior triangle energy, CONSTANT or ZERO
    EXT_GAM_TYPE ext_gam_flag;

    // Flag that determines whether to move all boundary vertices vs. a select
    // set of vertices belonging to select boundary strata.
    bool f_v2move;
    std::vector<std::pair<int, int> > c_2move;

    // Store the average cell width for all dimensions for meshadapt.
    double sim_sz;
    // The largest width of the simulation cell. Finds minima and maxima in 
    // xyz coords as vectors and takes the 2-norm of the difference vector.
    double sim_len;

    // Percentage of change of corresponding length/area/volume of a 
    // cell vs. it's length/area/volume for the given time step for the cell  
    // to be collapsing. This should be large enough so that it works
    // even if the cell lenght is too small, and the small time step due to 
    // collapsing cell is the limiting factor in the multiplication.
    double roc_tol;
    // The ratio of the minimum time step vs the time for a vertex to invert 
    // any tet by its motion.
    double t_rat;

    //double cell_meas_max;
    //double cell_meas_min;
    std::vector<double> avg_cell;
    // The minimum equivalent radii of strata of all dimensions.
    std::vector<double> min_cell;
    // The minimum radius of candidate 0strata for insertions. If this number 
    // gets smaller than a fraction of the topological length scale use this 
    // length for calculation of good quality to prevent extreme refinement.
    double min_ins_rad;

    // To reduce the cell radius computations, map checked cells.
    std::vector<std::map<int, bool> > rad_map;
    // Map cells smaller than threshold, but not shrinking.
    std::vector<std::map<int, bool> > rad_map_th;
    // Map cells that are below the threshold and shrinking or not 
    // couldn't be collapsed for some reason.
    std::vector<std::map<int, bool> > rad_map_cc;

    // Used in 1cell based tapered mesh adaptation scheme.
    // Map the minimum distance to the center of weight for each boundary cell. 
    // The vertices connected to cell vertices by edges will be weighted such
    // that there is at least a vertex within the cell.
    // If the multiplication weight*edge_length > 1.5 reference length, 
    // edge is split. If it is weight*edge_length < 0.5 reference length
    // it is coarsened.
    std::vector<std::map<int, double> > sm_map;
    std::vector<std::map<int, double> > meas_map;
    // The size and the time the stratum got below threshold. If it grows above 
    // this value, the value and time are reset. Above threshold dt_below_map is 
    // set to -1.
    std::vector<std::map<int, double> > meas_below_map;
    std::vector<std::map<int, double> > t_below_map;

    // Track the preconditioning sphere radii for 0strata. If the sphere is 
    // shrinking between iterations, and is below half of the threshold, do the
    // insertion.
    std::map<int, double> ins_rad_map;
    // Mesh related objects.
    //vd_glens* g_lens;
    // Cell related objects.
    struct cell_base* c_base;
    //struct cell_ins_chk c_ins;
    //vd_edisc* ed;
    //vd_energy vd_e;

    vd_entlist e_list;

    // The exterior shell for open surface simulations.
    ext_shell e_sh;

    // Parameters related to file operations.
    std::string smb_name;
    std::string vtk_name;

    // Adaptation parameter.
    double adapt_len;
    double ad_param;

    // meshadapt thresholds for splitting and coarsening edges. Currently only 
    // works with ModelEdgeRefinerDist. 
    double coarse_th;
    double split_th;

    // Controls whether adaptation is used or not. 
    bool ad_flag;
    // Controls adaptations functions that fix shape to actually fix shape or 
    // not. By default this should be true.
    bool ad_fix_shape;
    // Flag to check whether to collect widths used in adaptation or not.
    bool len_f;
    // Flag Are stratum widths collected.
    bool len_col_f;
    ADAPT_TYPE ad_type;

    // Used in adaptation scheme ADAPT_TYPE::ADAPT_BOUND
    std::vector<std::pair<int, int> > cells_ref;
    double bound_len;

    // Length scale related to the median cell length of dimension specified by 
    // ad_type. 
    double len_topo;
    // Length scale that specifies insertion and collapse lengths for the 
    // current iteration. 
    double len_trans;
    // Flags and parameters related to topological operations. 
    bool correct;
    double edge_min;
    double surf_min;
    double vol_min;

    // List for splitting 0cells.
    std::vector <int > tag_0cell_ins;
    // Keep track of unstable 0cells.
    std::map<int, bool> tag_0cell_stable;
    // Maps to track whether a collapse can be applied based on shell membership.
    // tag_shell_chk: Stratum has been checked for collapse.
    // tag_shell_map: Stratum can be collapsed.
    std::vector<std::map<int, bool> > tag_shell_chk;
    std::vector<std::map<int, bool> > tag_shell_map;

    // Number of subiterations:
    int iter_sub;

    double dt_adapt;
    // Cell collapse threshold ratio.
    double ratio_ent;

    // The thresholds for insertion and collapse sizes. 
    // Minimum preconditioning sphere radius = len_trans/ratio_ins
    // Minimum r_equi after shrink for actual collapse: len_trans/ratio_col_sh
    // Minimum r_equi for collapse check using only roc: len_trans/ratio_col
    double ratio_ins;
    double ratio_col_sh; 
    double ratio_col;

    // adapt_mark_2c_dist adapt_len to layer thickness ratio.
    double lay_rat;

    // MeshAdapt related parameters.
    // Tet quality threshold for quality operation:
    double q_th;
    // Fixing q_th can make the simulation become stuck when the smallest line
    // length is multiple orders of magnitude smaller than the target edge length.
    // q_th should be determined by the smallest edge and target edge lengths.
    double ln_min;
    // Tet volume threshold for triangle splitting operation:
    double vol_th;
    // Internal functions:

  public:
    // Flags and parameters related to the simulation.
    // Time related
    bool fixed_time;
    // Name of the csv file for storing time cost of each subprocess.
    std::string time_cost_str;
    // Timer object to measure the time elapsed for each subprocess.
    // Slot 0: Equations of motion
    // Slot 1: Adaptations
    // Slot 2: Insertions
    // Slot 3: Collapses
    // Slot 4: Auxilliary processes
    vd_timer vd_tim;
    // True: Cost time is accumulated and to be written to file at the end of
    // simulation.
    // False: Write to file at the end of each iteration.
    bool cost_t_accum;

    double time_curr;
    double time_end;
    double dt_set;
    double t_sub;
    
    // Minimum time step limit
    bool limit_flag;
    double t_sub_limit;

    bool flag_g_nbr;
    int cond_g_nbr;
    bool flag_s_nbr;
    int cond_s_nbr;
    bool flag_l_nbr;
    int cond_l_nbr;

    // The VTK save period. If < 0., it will save every time vtk_save() is called.
    // 
    double dt_vtk_period;
    double t_save_last;

    void iter_step();
    void iter_cells();
    void shr_axis(apf::Vector3 ax, double shr);

    // Given a same dim list of 0-,1-,2-cells and an axis, take the weighted 
    // average of the velocity. Evolve the vertices of the boundary cells until
    // the relative velocity is below a threshold of the weighted average.
    void evolve_bound_ax_conv(std::vector<std::pair<int, int> > * cells, 
                              apf::Vector3 ax, double th = 0.1);

    // Given a 1- or 2-cell, and a position, calculate the weighted average 
    // radial velocity. Evolve the vertices of the boundary cells until
    // the relative radial velocity is below a threshold of the weighted average.
    void evolve_bound_rad_conv(std::pair<int, int> * cell, 
                              apf::Vector3 pos, double th = 0.1);

    void evolve_bound_ax(std::vector<std::pair<int, int> > * cells, 
                              apf::Vector3 ax, double th = 0.1);
    void evolve_bound_rad(std::pair<int, int> * cell, 
                              apf::Vector3 pos, double th = 0.1);


    bool chk_vert_val(apf::MeshEntity* vert);

    // TODO make it private
    // Topology organizers; functions that collect all the objects of 
    // different modules that must be handled together:
    std::pair<int, int> col_cell(int cell_dim,int cell_id);

    // Given a set of tetrahedra belonging to a disjoint 3stratum, collect the 
    // bounding 2stratum triangles touching the central vertex. Calculate the 
    // anglecosine weighted sum of area normals. Also, return false if the 
    // anglecosine between the direction and one of plane normals is lower than a 
    // threshold.
    bool coll_bound_tri_dir(apf::MeshEntity* vert,
                                  std::vector<apf::MeshEntity*> &es_tet,
                                  apf::Vector3 &dir, apf::Vector3 pos,
                                                    double ang_th = 0.999);
    // Given a 0cell vertex, find the minimum of the distances towards the 
    // disjoint set of entities belonging to a 3cell around that vertex. 
    double dist_3c_min(apf::MeshEntity* vert);

    // Try inserting cells around the given 0cell. Currently only inserts 1cell.
    bool ins_cell(int tag_0cell);

    // ma::adapt() sometimes crashes due to non-existing cavities. Check if this
    // results from ma processes or insertion or collapses on our end.
    //bool chk_ma_swap();

    // Check if 0cell is insertible. 
    // In case of isotropic grain boundary energy, disallow if modified euler 
    // characteristic is satisfied. 
    bool insertible(int tag_0cell);

    // Check the volumes of the mesh.
    bool check_ma();
    bool check_ma(std::vector<apf::MeshEntity*>* verts);

    bool check_ma_sgn();
    bool check_ma_sgn(std::vector<apf::MeshEntity*>* verts);

    // If c_ins and c_base are not synchronized, reload c_ins. Assuming each 
    // change in the cell complex is through vd_sim, and updates both c_ins
    // and c_base, only compare the total number of 0cells. Return true if 
    // c_ins is not synchronized and is updated.
    //bool sync_cins();

    // Default constructor:
    vd_sim();
    void set_mesh(const char* modelFile);
    void set_mesh(const char* modelFile, const char* meshFile);
    void set_mesh_smb(const char* modelFile, const char* meshFile);

    cell_base* get_c_base();
    //cell_ins_chk* get_c_ins();
    bool ret_ins_gmi(int c0_id);
    void ret_ins_gmi();


    void corr_0cell();
    void corr_1cell();

    void find_0cell();
    void start_sim();

    // Mesh save and load functions:
    void save_mesh(char const* FileName ="./tempmesh/beforecol");
    void reload_len();
    void check_model();
    void reload_mesh();

    bool comp_mesh_topo();

    // Set the vtk save interval parameters and check if enough time has passed
    // for automatic vtk saving.
    void set_save_vtk_interval(double dt_period_in);
    bool chk_save_vtk_interval();

    // Generates a vtk file at the simulation step.
    void save_vtk();
    // Generates a vtk per grain. Used in generating movies in a hackish way.
    void save_vtk_mov();
    void save_vtk_name(const char* vtk_out, bool smb_out = false);

    apf::Mesh2* get_mesh();
    void verify();

    // Get the length scale associated with insertions. TODO median 
    // approximation, median, average, set value should be options.
    // Entity based, cell width based should be options, as well.
    std::pair<double, double> get_length_scale();

    // Update the normal distance field for the vertices of the given cell, 
    // from the one lower dimensional entities on the bounding strata.
    void upd_dist_field(int c_dim, int c_id, 
                            std::map<apf::MeshEntity*, double> &v_dist);


    bool chk_cell_th(int d, int cell_id, double &r_equi);
    // In order to reduce the computational load in cell collapse check,
    // a cell is skipped if it has a size 5*average edge size. This is not an
    // exact comparison (edge length vs. 3cell volume (2cell area), but it 
    // should catch all positives, which is what is desired.
    std::vector<double> upd_cell_rad_sim();
    void adapt_mark_step();
    void adapt_mark_1c_min();

    void adapt_mark_0c_col(int d_col, int t_col);
    void adapt_mark_0c_min();
    void adapt_mark_bound_min();

    // Based on the normal distance from the 2strata, tag each vertex.
    // Use the distance field to adapt the mesh.
    // Can potentially be also used to detect transitions such as intersecting
    // boundaries (e.g. in the case of recrystallizing grains) 
    // Starting from each 2 stratum tri, mark the adjacent vertices as burnt
    // and set the weight 0.
    // Set the starting vertices in the first tri list.
    // Get the tet adjacencies of the tris, put them in a list. 
    // If tet unburnt, check:
    // If only 3 vertices are burnt, calc distance of the last vert. Assign tet
    // burnt. 
    // If the distance to the vert is less than threshold:
    //   Put the other tris in the next tri list. 
    // Put the vertex into the toburn list.
    // If 4 vertices are burnt, burn the tet.
    // Otherwise, keep the tet in the list by shifting
    // Burn the vertices in the to burn list.
    // Using the next tri list, repeat the process. 
    void adapt_mark_2c_dist(std::map<apf::MeshEntity*, double> &v_dist,
std::vector<double> &dist_c3);

    void adapt_iso();

    double get_lay_thk();
    // Return the median width of the 1cells/3cells divided by the adapt_len and
    // ad_param.
    double get_adapt_ln();

    void fix_low_q(); 
    // Given a cell, coarsen the edges emanating from bounding vertices within a
    // sphere with a radius twice the equivalent radius of the given cell. 
    void adapt_coarsen_cell(int cell_dim, int cell_id);

    // Given the low quality tets, set the edge belonging to the lowest  
    // dimensional cell to be coarsened. Only coarsen.
    void adapt_prob_edge(std::vector<apf::MeshEntity*> &tet); 
    // If adapt_prob_edge doesn't reduce the number of low quality tets, try 
    // marking all edges adjacent to the vertices to be coarsened.
    void adapt_prob_edge_all(std::vector<apf::MeshEntity*> &tet);
    void adapt_prob_edge_low(std::vector<apf::MeshEntity*> &tet);

    void adapt_prob_edge_sp_bound(std::vector<apf::MeshEntity*> &tet);

    // How many times to hit bound_manual before running the meshadapt fix for
    // low quality elements at the boundaries. It is not always run, but only
    // when the normal improvements for pancake-like elements do not work.
    int bound_manual_count_limit;
    int bound_manual_count;

    void adapt_prob_edge_col_bound_manual(std::vector<apf::MeshEntity*> &tet);
    void adapt_prob_edge_sp_bound_manual(std::vector<apf::MeshEntity*> &tet);
    // Given the low quality tets, set the longest edge of the tets to be split.
    // Only split.
    void adapt_prob_edge_sp(std::vector<apf::MeshEntity*> &tet); 
    void adapt_prob_edge_sp2(std::vector<apf::MeshEntity*> &tet); 
    void adapt_prob_edge_sp3(std::vector<apf::MeshEntity*> &tet); 

    double get_good_q();

    void adapt(bool avg_flag = true);
    void adapt_col(int dim_col, int tag_col);

    // Clean mesh related objects and flags.
    void clean_mesh();
    // Clean all objects in the heap.
    void clean_up();

    bool set_free_cells();

    void get_v2move(std::vector<apf::MeshEntity*>& v_2move);
    void set_c2move(std::vector<std::pair<int, int> >& c_2move_in, 
                                                  bool f_v2move_in);
    // Refine mesh.
    void ref_mesh_iso(double len);

    // Check whether a collapsing stratum is prevented from collapsing due to 
    // shell membership.
    bool chk_shell_issue(int d, int cell_id);

    // Find cells below the length, area or volume threshold.
    // For now the half of the average entity size is used as the threshold.
    std::pair<int, int> chk_small_cell_rad();
    std::pair<int, int> chk_cell_rmv_up(int d, int cell_id);
/*
    bool chk_col_energy(std::vector<apf::MeshEntity* >* ev, 
                                                  int d, int cell_id);
*/

    apf::Vector3 get_cell_midpoint(int d, int cell_id);

    // Check if the triangles surrounding the collapsing entity are farther 
    // than the collapsing cell vertices.
    // TODO recycle this or glens portion that does the same check.
    bool chk_cell_thresh(int d, int cell_id);
    void clear_rad_map();

    std::pair<int, int> chk_col_cell_adj(int d, int cell_id);
    bool chk_cell_rmv(int d, int cell_id);
    //bool chk_cell_shr(int d, int cell_id);
    //bool chk_cell_shr_rad(std::vector<apf::MeshEntity* >* ev, 
    //                                        int d, int cell_id);
    //bool chk_cell_shr_rad(int d, int cell_id);

    // Used with adaptation scheme ADAPT::ADAPT_BOUND. These cells are refined 
    bool set_adapt_bound_cells(double ratio, 
                               std::vector<std::pair<int, int> > &cells_ref_in);

    // Change the simulation parameters.
    void set_adapt_param(double param_in);
    void set_adapt_type(ADAPT_TYPE ad_in);
    // Change the vertex motion function.
    void set_field_calc(const field_calc& FC);
    void set_field_calc(VEL_TYPE vel_in);
    void set_vec_sp_calc(PROJ_TYPE PROJ);

    // Set whether the exterior vertices are moving or not.
    void set_calc_ext(bool calc_ext_in);

    // Set the integration type.
    void set_integ_type(INTEG_TYPE IT);

    // Mesh level exterior collection. The naive collection on c_base requires
    // side 1cells to be connected to a single 3cell.
    void update_shell_gmi(int dim, int tag, int tag_0cell);
    void burn_shell(int dim, int tag, shell sh);
    void collect_ext();
    void collect_ext_dir();

    // Set the time and whether it can be increased as mesh coarsens.
    void set_time(double time_in);
    void set_time(std::vector<std::string>& opts);

    double get_time();
    void set_end(float time_e);
    void set_dt(double dt_in, bool fixed_in = false);

    void set_iter_sub(int sub_nbr = 10);

    void write_time_cost();
    void set_time_cost_str(const std::string str_in, bool accumulate = true);

    // Set the cell split criteria:
    // 1- Based on modified euler characteristic,
    // 2- Based on decreasing energy.
    void set_euler_char(bool on_off);
    void set_ins_flag(bool on_off);
    void set_col_flag(bool on_off);

    void set_col_rat(double ratio_col_sh_in, double ratio_col_in);
    void set_ins_rat(double ratio_ins_in);

    void set_adapt(double param_in);
    void set_ad_flag(bool on_off);
    void set_ad_fix_flag(bool on_off);
    void set_ad_th(double coarse_in, double split_in);

    //void set_extract(bool on_off, std::vector<vd_ext_opts>* opts);
    //void set_extract_cells(int i, const std::vector<c_p>* c_top,
    //                    const std::vector<std::vector<c_p> >* c_bound);

    //void set_extract_MS(std::string MS_i, std::string ROC_i, std::string VOL_i, 
    //                    bool on_off = true);

    void set_equations(std::vector<std::string>& opts);
    void set_c2move_opts(std::vector<std::string>& opts, bool f_v2move_in);

    void set_extract(bool on_off, std::string& opts);
    // Extract data as specified in the opts file.
    void extract_data();
    // Extract data as related to MacPherson Srolovitz relation.
    //void extract_MS();

    // Test operations:
    // Adaptively shrink a cell to a size relative to its average entity size.
    double shr_cell(int cell_dim, int cell_id, double sp_size = 1);
    bool chk_cell_shr(int d, int cell_id);

    // Input functions:
    void set_options(const char* file_name);

    // Output functions:

    // Return the simulation parameters.
    vd_param get_vd_param();

    void set_ext_gam_flag(EXT_GAM_TYPE ext_in);
    void set_vd_param(vd_param vd_par);

    field_calc* get_f_calc();
    double get_sim_len();
    double get_len_topo();

    vd_entlist* get_elist();

    // Refresh the elist.
    void refresh_elist();

    // Print the 2-3cell adjacencies around the 0cells.
    void print_0cellgraph();
    void print_0cell_pc();

    ~vd_sim();

    // Set vtk output type. Right now, in order generate movies, the mesh of
    // each 3cell should be saved seperately. This flags controls whether they
    // are seperately saved or not.
    void set_mov_flag(bool mov_flag);

};

#endif
