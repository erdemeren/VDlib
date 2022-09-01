#ifndef TOPO_FEAT_H
#define TOPO_FEAT_H
#include <math.h>
#include <time.h>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_topo.h"
#include "topo_field.h"
#include "topo_entlist.h"

#include "topo_tmesh.h"

typedef std::pair<int, int> c_p;

int find(std::vector<c_p>* cell_list, int id);

class vd_ext_opts;
class vd_series_ss;
class vd_series_sv;
class vd_series_vs;
class vd_series_vv;
class ext_mesh;
class ext_topo;
class vd_processor;

// TODO vd_sim::extract_data() function should be restructured to accept custom
// vd_extractor objects.

// Boundary condition type. Mainly used in constant of mobility calculation.
// (2pi - 3acos(-1/3)) per adjacent 3-stratum for interior vertices.
// NEUMANN_EXT adds the exterior 0-strata to the integral of the Gaussian curv.
enum class VD_BC_TYPE {
  PERIODIC,
  PERIODIC_CONST,
  NEUMANN_EXT,
  NEUMANN_CONST,
  END
};

enum class VD_SEARCH {
  AVERAGE,
  TOTAL,
  MINIMUM,
  MAXIMUM,
  MAX_DIR,
  MIN_DIR,
  END
};

// Evaluator functions. 
apf::Vector3 vd_ext_pos(apf::Mesh2* m, apf::MeshEntity* e);
double vd_ext_field_db(apf::Mesh2* m, apf::MeshEntity* e, vd_ext_opts* opts);
long int vd_ext_field_li(apf::Mesh2* m, apf::MeshEntity* e, vd_ext_opts* opts);
apf::Vector3 vd_ext_field_v(apf::Mesh2* m, apf::MeshEntity* e, vd_ext_opts* opts);

void vd_calc_field_dummy(apf::Mesh2* m, apf::MeshEntity* e, vd_ext_opts* opts);

double vd_calc_curv(apf::Mesh2* m, apf::MeshEntity* e);

void write_data(vd_series_ss* data, const char* fname = "./output/data.csv");
void write_whitespace_x_vert(vd_series_ss* data, 
                                      const char* fname = "./output/data.csv");
void write_label_y_vert(vd_series_ss* data, 
                                      const char* fname = "./output/data.csv");

// TODO use templates, less ugly and probably easier to maintain// Used to keep options for what to extract from the mesh. 
class vd_series_ss {
  public:
    std::string x_name;
    std::string y_name;
    std::string title;
    std::string variable;
    std::vector<int> x;
    std::vector<double> y;

    void clear();

    // Copy constructor
    vd_series_ss(const vd_series_ss& that);
    // Copy
    vd_series_ss& operator=(const vd_series_ss& that);

    vd_series_ss();
    ~vd_series_ss();
};

// Used to keep options for what to extract from the mesh. 
class vd_series_sv {
  public:
    std::string x_name;
    std::string y_name;
    std::string title;
    std::string variable;
    std::vector<int> x;
    std::vector<apf::Vector3> y;

    void clear();

    // Copy constructor
    vd_series_sv(const vd_series_sv& that);
    // Copy
    vd_series_sv& operator=(const vd_series_sv& that);

    vd_series_sv();
    ~vd_series_sv();
};

// Used to keep options for what to extract from the mesh. 
class vd_series_vs {
  public:
    std::string x_name;
    std::string y_name;
    std::string title;
    std::string variable;
    std::vector<apf::Vector3> x;
    std::vector<int> y;

    void clear();

    // Copy constructor
    vd_series_vs(const vd_series_vs& that);
    // Copy
    vd_series_vs& operator=(const vd_series_vs& that);

    vd_series_vs();
    ~vd_series_vs();
};

// Used to keep options for what to extract from the mesh. 
class vd_series_vv {
  public:
    std::string x_name;
    std::string y_name;
    std::string title;
    std::string variable;
    std::vector<apf::Vector3> x;
    std::vector<apf::Vector3> y;

    void clear();

    // Copy constructor
    vd_series_vv(const vd_series_vv& that);
    // Copy
    vd_series_vv& operator=(const vd_series_vv& that);

    vd_series_vv();
    ~vd_series_vv();
};

// Used to measure process time. 
class vd_timer {
  private:
  public:
    std::vector<double> begin;
    std::vector<double> end; 
    std::vector<double> elapsed;

    void set_time(int slot);
    double get_time(int slot);
    double add_time(int slot);
    double get_elapsed(int slot);
    double get_elapsed_sum();
    void reset_elapsed();

    void resize(int slots);
    //vd_timer(int TIME_TYPE = CLOCK_PROCESS_CPUTIME_ID);
    vd_timer();
    vd_timer(int slots);
    ~vd_timer();
};

// Used to keep options for what to extract from the mesh. 
class vd_ext_opts {
  public:
    // Parameters
    bool n_flag;
    std::string n_name;
    std::string file_name;

    bool crt_flag;
    bool f_flag;
    int f_type;

    bool first_ext;

    // The cell list of interest.
    std::vector<c_p> c_top;
    std::vector<std::vector<c_p> > c_bound;

    // Defines the value selection for cell based extractions. Default is 
    // average.
    int search_type;
    // Searching direction.
    apf::Vector3 search_dir;

    // Functions
    void (*vd_calc_field_pt)(apf::Mesh2*, apf::MeshEntity*, vd_ext_opts*);
    void vd_calc_field(apf::Mesh2* m, apf::MeshEntity* e);

    void clear();

    void set_func(void (*func_pt_in)(apf::Mesh2*, apf::MeshEntity*, 
            vd_ext_opts*), char* n_in, int f_type_in, bool crt_in = true);
    void set_search(VD_SEARCH s_type, apf::Vector3 dir_in);

    void set_extract_cells(const std::vector<c_p>* c_top_in,
                        const std::vector<std::vector<c_p> >* c_bound_in);

    // Copy constructor
    vd_ext_opts();
    vd_ext_opts(const vd_ext_opts& that);
    // Copy
    vd_ext_opts& operator=(const vd_ext_opts& that);

    vd_ext_opts(char* field_name, bool crt_in = true);
    ~vd_ext_opts();
};

// Used to extract information from the mesh. 
// Mainly field information(including position).
class ext_mesh {

  public:
    vd_ext_opts m_opts;
    apf::Mesh2* m;
    cell_base* cb;

    void set_opts(vd_ext_opts opt_in);
    bool chk_field();
    void create_field();

    void set_mesh(apf::Mesh2* m_in, cell_base* cb_in);
    double get_avg_db(vd_entlist* e_list, int dim, int c_id);

    void calc_total_curv(vd_series_ss* data, int dim);

    void get_per_cell_dim(std::vector<vd_series_ss>* data, int dim);
    void get_per_cell_dim(std::vector<std::vector<vd_series_ss> >* data,
                                                                   int dim);
    void get_per_cell_all(std::vector<vd_series_ss>* data);
    void get_meas_cell_all(std::vector<vd_series_ss>* data);
    void get_cell_meas(vd_series_ss* data, std::vector<c_p>* c_top);

    // Get per bounding cell of dim_b of the given cell.
    void get_cell_bound(vd_series_ss* data, int dim, int c_id, 
                                          std::vector<c_p>* c_bound);
    void get_cell_bound(vd_series_sv* data, int dim, int c_id,
                                          std::vector<c_p>* c_bound);

    void get_cell_bound_ext_angle(vd_series_ss* data, int dim, int c_id,
                                          std::vector<c_p>* c_bound);


    void calc_field();
    void clear();

    ext_mesh();
    ~ext_mesh();
};

// Used to extract information from the topology. 
class ext_topo {

  public:
    cell_base* cb;
    void get_conn_dist(std::vector<vd_series_ss>* data);
    void get_per_cell_dim(std::vector<vd_series_ss>* data, int dim);

    void set_cbase(cell_base* cb_in);

    void clear();

    ext_topo();
    ~ext_topo();
};

// Extract MacPherson Srolovitz related quantities for each
// TODO include the exterior shell correctly.
class ext_MS {

  public:
    cell_base* cb;
    apf::Mesh2* m;
    vd_entlist e_list;

    //PROJ_TYPE proj_flag;
    //ext_shell e_sh;

    // Rate of change of volume of c3
    std::map<int, double> dVdt_c3;
    // Lengths of c1
    std::map<int, double> len_c1;

    // TODO Not used at the moment. c1c3 adjacencies and relative changes of 
    // surrounding c3 might yield some information. Using modified version of MS
    // to account for valencies higher than 3 for lines makes this somewhat  
    // obsolete...
    std::vector<std::vector<int> > c1c3;
    std::vector<std::vector<double> > dVdt_c1c3;

    std::vector<std::vector<int> > c2c3;
    std::vector<std::vector<double> > dVdt_c2c3;

    double check_ext_angle(apf::MeshEntity* tet, apf::MeshEntity* tri1, 
                           apf::MeshEntity* tri2, apf::MeshEntity* edge);

    void process_bound(int dim, int id);
    void process_mesh();

    void clear();

    ext_MS();
    ext_MS(apf::Mesh2* m_in, cell_base* cb_in);
    ~ext_MS();
};

#define ANG_EQUI_CONST (2*PI_L - 3*std::acos(-1./3.))
#define ANG_EQUI_CONST_EXT_0C (PI_L/2)
#define ANG_EQUI_CONST_EXT_1C (PI_L)
#define ANG_EQUI_CONST_EXT_2C (PI_L/3)

// Extract Gaussian curvature based constant of mobility.
// Iterate over boundary vertices, extract the disjoint 3-strata in the 
// neighborhood.
// TODO FUTURE Constant of mobility article reference...
// Both the actual contribution of the Gaussian curvature at the 0-strata and
// the constant of motion based on the equilibrium angular defect angle 
// assumption are calculated. 
// In the case of VD_BC_TYPEs ending with _CONST, the 0-cells under transitions 
// are not handled differently. In the normal case, as the transitioning 0-strata
// have unstable configurations, they are directly included in the integral of 
// the Gaussian curvature term.
class ext_GAUSS {

  public:
    cell_base* cb;
    apf::Mesh2* m;
    vd_entlist* e_list;
    field_calc* f_calc;

    VD_BC_TYPE ext_BC_type;

    apf::Field* ang_def_field;

    int dim;
    int c_id;
    int c3_nbr;

    apf::MeshEntity* v;
    apf::Vector3 temp_pos;
    double ang_def_equi;

    // The current vertex is exterior or not. Exterior vertices have their 
    // equilibrium angle added to constant of motion at the moment.
    bool ext_curr;

    // Validation:
    // CoM = gauss - (ang_defect + ang_defect_trans)
    // CoM = gauss - (ang_defect_const + ang_defect_dev + ang_defect_trans)
    // Constant of motion, integral of the Gaussian curvature over 2-strata
    // The contributions of the internal 2-strata added twice.
    double CoM;
    // Sum of the angular defects at 0-strata not under transition.
    double ang_defect;
    // Angular defect due to transition. This is used to distinguish the error
    // due to deviation from the assumed angular defect configuration and the
    // error due to unstable configuration during transition.
    double ang_defect_trans;
    // Angular defect as calculated by isotropic GBE constant angular defect 
    // assumption: f_0 * (2*pi - 3*acos(-1/3)), f_0 number of 0-strata
    double ang_defect_const;

    // The sum of the deviations of the angular defects from the constant GBE 
    // case.
    double ang_defect_dev;

    // The sum of the deviations of the angular defects around the interior 1cells
    double ang_defect_1c_int;

    // Updated during calc_gauss. Used in updating the ang_defect_const or 
    // ang_defect_const_ext and calculating the deviation term.
    double ang_equi;

    // Sum of the angular defects for the 3-strata. 
    // Should be 2pi * Euler characteristic
    // For validation.
    std::vector<double > c3_euler_mod;
    // Sum of the angular defects for the 0-strata. 
    // For testing.
    std::vector<double > c0_ang_def;
    // The deviations of the angular defects from the ideal configuration. 
    // Transitioning 0-strata are assigned -1 for filtering.
    std::vector<double > c0_ang_dev;
    //std::vector<std::vector<double > > c0_c3_ang;
    // The magnitude of velocity.
    std::vector<double > c0_vel;

    // The lengths and volumes of the 1- and 3-strata.
    std::vector<double > c1_len;
    std::vector<double > c3_vol;

    // The 1- and 3-strata connectivity lists of 0-strata.
    std::vector<std::vector<int> > c0_c1;
    std::vector<std::vector<int> > c0_c3;

    // Total of the Gaussian curvature. Should be f_3 * 4pi, f_3 number of 3strata
    double gauss;

    // The local neighborhood of a vertex.
    vd_3c_det* vd_3c;

    // The orientations of the boundary edges. For faster interior angle 
    // calculation.
    std::map<apf::MeshEntity*, apf::Vector3 > ori_map;
    std::map<apf::MeshEntity*, apf::Vector3 > pos_map;

    // The interior angles of triangles, at bounding vertices. Same t_id as 
    // vd_edisc ids.
    // ang_map.at(c2_id).at(t_id)[vertex]
    std::map<apf::MeshEntity*, std::map<apf::MeshEntity*, double > > ang_map;
    std::map<apf::MeshEntity*, std::map<apf::MeshEntity*, int > > ang_map_count;

    // For a given vertex, store the angles of the disjoint 2-strata.
    std::map<int, double> c2_map;
    std::map<int, double> c2_ang;

    // Return the interior angle associated with a triangle at the vertex.
    double get_ang(apf::MeshEntity* t);
    // Collect the angles associated with the adjacent 2-strata of the current
    // (0- or 1-stratum) vertex.
    void collect_ang();

    // Collect the local neighborhood of the current (0- or 1-stratum) vertex.
    void load_v();

    // Calculate the angular defects associated with a vertex. 
    // Add the results to the associated terms.
    // Doing the angular defect calculation per triangle would require 1/3 less 
    // dot products but doing it over vertices is less complex:
    // The 3-stratum valencies of the 0-strata are necessary to check anyways.
    // Add the contributions to adjacent 3-strata for validation.
    void process_0c_vertex();
    void process_ext_1c_vertex();
    void process_int_1c_vertex();
    void process_2c_vertex();

    double get_ang_val(apf::MeshEntity* tri, apf::MeshEntity* vert);
    int get_ang_count(apf::MeshEntity* tri, apf::MeshEntity* vert);


    void process_ext_edge(apf::MeshEntity* e);
    void process_tri(apf::MeshEntity* t);

    // Iterate over the vertices to calculate the invariant related quantities.
    void calc_gauss();

    // Collect the orientation of the edges.
    void collect_ori();

    // Get the 1- and 3-strata connections of 0-strata.
    void get_conn();
    // Get the lengths and volumes of the 1- and 3-strata.
    void get_len();

    void process_mesh();

    void set_bc_type(VD_BC_TYPE ext_BC_in);

    void clear();

    ext_GAUSS();
    ext_GAUSS(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* e_list_in, 
                                                    field_calc* f_calc_in);
    ~ext_GAUSS();
};

// Extract Gaussian curvature based constant of mobility.
// Extracts information per grain, i.e. the integral of Gaussian curvature for bounding 2-, 1-, and 0-strata. 
// Iterate over boundary vertices, extract the disjoint 3-strata in the 
// neighborhood.
// TODO FUTURE Constant of mobility article reference...
// Both the actual contribution of the Gaussian curvature at the 0-strata and
// the constant of motion based on the equilibrium angular defect angle 
// assumption are calculated. 
// In the case of VD_BC_TYPEs ending with _CONST, the 0-cells under transitions 
// are not handled differently. In the normal case, as the transitioning 0-strata
// have unstable configurations, they are directly included in the integral of 
// the Gaussian curvature term.
class ext_GAUSS_g {

  public:
    cell_base* cb;
    apf::Mesh2* m;
    vd_entlist* e_list;
    field_calc* f_calc;

    VD_BC_TYPE ext_BC_type;

    apf::Field* ang_def_field;

    int dim;
    int c_id;
    int c3_nbr;

    std::vector<int> g_id;

    apf::MeshEntity* v;
    apf::Vector3 temp_pos;
    double ang_def_equi;

    // The current vertex is exterior or not. Exterior vertices have their 
    // equilibrium angle added to constant of motion at the moment.
    bool ext_curr;

    // Validation:
    // CoM = Integral over 2-strata
    // ang_1c = Integral over 1-strata
    // ang_0c = Integral over 0-strata
    // Constant of motion, integral of the Gaussian curvature over 2-strata
    // of the selected grains.
    double CoM;
    double ang_0c;
    double ang_1c;

    // Updated during calc_gauss. Used in updating the ang_defect_const or 
    // ang_defect_const_ext and calculating the deviation term.
    double ang_equi;

    // Sum of the angular defects for the 3-strata. 
    // Should be 2pi * Euler characteristic
    // For validation.
    std::vector<double > c3_euler_mod;
    // Sum of the angular defects for the 0-strata. 
    // For testing.
    std::vector<double > c0_ang_def;
    // The deviations of the angular defects from the ideal configuration. 
    // Transitioning 0-strata are assigned -1 for filtering.
    std::vector<double > c0_ang_dev;
    //std::vector<std::vector<double > > c0_c3_ang;
    // The magnitude of velocity.
    std::vector<double > c0_vel;

    // The lengths and volumes of the 1- and 3-strata.
    std::vector<double > c1_len;
    std::vector<double > c3_vol;

    // The 1- and 3-strata connectivity lists of 0-strata.
    std::vector<std::vector<int> > c0_c1;
    std::vector<std::vector<int> > c0_c3;

    // Total of the Gaussian curvature. Should be f_3 * 4pi, f_3 number of 3strata
    double gauss;

    // The local neighborhood of a vertex.
    vd_3c_det* vd_3c;

    // The orientations of the boundary edges. For faster interior angle 
    // calculation.
    std::map<apf::MeshEntity*, apf::Vector3 > ori_map;
    std::map<apf::MeshEntity*, apf::Vector3 > pos_map;

    // The interior angles of triangles, at bounding vertices. Same t_id as 
    // vd_edisc ids.
    // ang_map.at(c2_id).at(t_id)[vertex]
    std::map<apf::MeshEntity*, std::map<apf::MeshEntity*, double > > ang_map;
    std::map<apf::MeshEntity*, std::map<apf::MeshEntity*, int > > ang_map_count;

    // For a given vertex, store the angles of the disjoint 2-strata.
    std::map<int, double> c2_map;
    std::map<int, double> c2_ang;

    // Return the interior angle associated with a triangle at the vertex.
    double get_ang(apf::MeshEntity* t);
    // Collect the angles associated with the adjacent 2-strata of the current
    // (0- or 1-stratum) vertex.
    void collect_ang();

    // Collect the local neighborhood of the current (0- or 1-stratum) vertex.
    void load_v();

    // Calculate the angular defects associated with a vertex. 
    // Add the results to the associated terms.
    // Doing the angular defect calculation per triangle would require 1/3 less 
    // dot products but doing it over vertices is less complex:
    // The 3-stratum valencies of the 0-strata are necessary to check anyways.
    // Add the contributions to adjacent 3-strata for validation.
    void process_0c_vertex();
    void process_ext_1c_vertex();
    void process_int_1c_vertex();
    void process_2c_vertex();

    double get_ang_val(apf::MeshEntity* tri, apf::MeshEntity* vert);
    int get_ang_count(apf::MeshEntity* tri, apf::MeshEntity* vert);


    void process_ext_edge(apf::MeshEntity* e);
    void process_tri(apf::MeshEntity* t);

    // Iterate over the vertices to calculate the invariant related quantities.
    void calc_gauss();

    // Collect the orientation of the edges.
    void collect_ori();

    // Get the 1- and 3-strata connections of 0-strata.
    void get_conn();
    // Get the lengths and volumes of the 1- and 3-strata.
    void get_len();

    void process_mesh();

    void set_bc_type(VD_BC_TYPE ext_BC_in);

    void clear();

    ext_GAUSS_g();
    ext_GAUSS_g(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* e_list_in, 
                                                    field_calc* f_calc_in);
    ~ext_GAUSS_g();
};

enum class EXT_TYPE {
  MS,        // Extract the volumes and the rates of changes of volume by MS, ROC
  VEL,       // Given a stratum, extract the weigted position
  POS,       // Given a stratum, extract the weigted position
  VEL_DIR,   // Given a stratum, extract the component of weigted velocity in dir
  POS_DIR,   // Given a stratum, extract the component of weigted position in dir
  CRV,   // Given an origin and 2stratum, extract mean curvature, radius and 
         // radial component of velocity (Area averaged)
  MEAS,  // Given a list of cells, extract the length/area/volume of the 
         // individual cells.
  AREA_T,// Extract the total area.
  GAUSS, // Given a mesh, extract the constant of motion.
  GAUSS_G, // Given a mesh, extract the constant of motion.
  ROC, // Given a list of strata, extract the rates of change of length, area, vol
  COUNT,     // Count of the d-dim entities of strata with c-dim c-tag.
             // COUNT; d-ent d-cell id-cell; ...;
  ANGLE1, // Given triplets of strata (a joint 0stratum and two 1strata joining 
          // at the stratum) extract the interior angle. Assumes single 
          // connection for each joining strata at the mesh level.
  ANGLE2, // Given quadruplet of strata (S1, S3, S20, S21) a joint 1stratum, a 
          // 3stratum and two 2strata joining at the 1stratum with the interior  
          // as the 3stratum) extract the interior angle. Assumes single 
          // connection for each joining strata at the mesh level.
          // ANGLE1; <ANGLE CSV FILE>; S1; S3; S20; S21
  ANGLE_CS, // Given quadruplet of strata (S3, S0/1, S1/2_0, S1/2_1), a 3stratum 
          // a joint 0/1stratum, and two bounding 1/2-strata (joining at the 
          // 1stratum with the interior as the 3stratum) and a plane intersecting 
          // the 0-stratum, extract the interior angle. Assumes single connection 
          // for each joining strata at the mesh level. Plane defined by point 
          // (P1 P2 P3), normal (N1 N2 N3). 
          // ANGLE_CS; <ANGLE CSV FILE>; S3; D_j S_j; D_b1 S_b1; D_b2 S_b2; P1 P2 P3; N1 N2 N3;
  HAUSDORFF, // Given a boundary 2-stratum, and a reference triangulation of the 
             // same 2-stratum as an ".off" file, extract the triangulation of 
             // the current mesh and run an exterior CGAL program to write the 
             // Hausdorff distance between the two triangulations to given .csv 
             // file.
  END
};

EXT_TYPE conv_str2ext(const std::string& str);

// Derived classes shouldn't require rewriting vd_sim. There should be an input 
// options file that will parsed and derived extraction classes called 
// accordingly.
class vd_extractor {

  protected:
    cell_base* cb;
    apf::Mesh2* m;
    vd_entlist* e_list;
    field_calc* f_calc;
    double t_curr;

    virtual void parse_opt(std::vector<std::string> &opts);
  public:
    vd_extractor();
    vd_extractor(apf::Mesh2* m_in, cell_base* cb_in, 
                             vd_entlist* e_list_in, field_calc* f_calc, 
                             std::vector<std::string> & opts_in, double t_in);
    virtual void write_csv();
    ~vd_extractor();
};

// A wrapper for ext_MS in style of vd_extractor. 
// TODO incorporate ext_MS into this.
class vd_ext_MS : protected vd_extractor {
  private:
    ext_MS e_MS;

    std::string OUTPUTMS;
    std::string OUTPUTROC;
    std::string OUTPUTVOL;
    double dt;
    void parse_opt(std::vector<std::string> &opts);
  public:
    vd_ext_MS(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc, 
                               std::vector<std::string> & opts_in,
                               double t_in);
    void write_csv();
};

// A wrapper for ext_GAUSS in style of vd_extractor. 
// TODO incorporate ext_GAUSS into this.
class vd_ext_GAUSS : protected vd_extractor {
  private:
    ext_GAUSS e_GAUSS;
    VD_BC_TYPE bc_type;

    std::string OUTPUTGAUSS;
    // The sum of the CoM and the angular defects at 0-strata for each grain. 
    // Should yield 2pi*Euler characteristic.
    std::string OUTPUTEULERMOD;
    // The deviations of the modified euler characteristic from the ideal 
    // expression for CoM + ang_def_equi*f_0.
    std::string OUTPUTERROR;
    // The deviations from the ideal configuration for the 0-strata.
    std::string OUTPUT0C_dev;
    // The lengths of the 1-strata.
    std::string OUTPUT1C;
    // The volumes of the 3-strata.
    std::string OUTPUT3C;
    // The connectivity of the 0-strata to the 1-strata.
    std::string OUTPUT1C_conn;
    // The connectivity of the 0-strata to the 3-strata.
    std::string OUTPUT3C_conn;

    // The velocities of the 0-strata vertices.
    std::string OUTPUT0C_vel;

    void parse_opt(std::vector<std::string> &opts);
  public:
    vd_ext_GAUSS(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc, 
                               std::vector<std::string> & opts_in,
                               double t_in);
    void write_csv();
};

// A wrapper for ext_GAUSS_g in style of vd_extractor. 
class vd_ext_GAUSS_g : protected vd_extractor {
  private:
    //ext_GAUSS_g e_GAUSS;
    VD_BC_TYPE bc_type;
    std::vector<int> g_id;
    std::vector<double> CoM;
    std::vector<double> ang_def_dev_0c;
    std::vector<double> ang_def_dev_1c;
    double CoM_avg;

    std::map<apf::MeshEntity*, apf::Vector3 > ori_map;
    std::map<apf::MeshEntity*, apf::Vector3 > pos_map;
    std::map<int, bool> c_pass;
    // Unique value per 2-strata, 4 values per 0-strata, 3 per 1-strata
    // ang_def_1c and ang_def_0c are to check for consistency of the Gaussian 
    // integral.
    std::map<int, double> ang_def_2c;
    std::map<int, double> ang_def_1c;
    std::map<int, double> ang_def_0c;

    // The sum of the CoM and the angular defects at 0-strata for each grain. 
    // Should yield 2pi*Euler characteristic.
    std::string OUTPUTGAUSS;

    void calc_ang_2c(int c_in);
    void collect_ori();
    void collect_ang();

    void parse_opt(std::vector<std::string> &opts);
  public:
    vd_ext_GAUSS_g(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc, 
                               std::vector<std::string> & opts_in,
                               double t_in);
    void write_csv();
};


// 
// ROC; <ROC CSV FILE>; S_DIM_1 S_ID_1; ...; S_DIM_N S_ID_N;
class vd_ext_ROC {

  public:
    cell_base* cb;
    apf::Mesh2* m;
    vd_entlist* e_list;
    field_calc* f_calc;

    std::vector<std::pair<int, int> > c_list;

    double t_curr;

    // The connectivity of the 0-strata to the 3-strata.
    std::string OUTPUTROC;

    void parse_opt(std::vector<std::string> &opts);
    void write_csv();

    vd_ext_ROC();
    vd_ext_ROC(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in);
};


// 
// COUNT; <COUNT CSV FILE>; d-ent S_DIM_1 S_ID_1; ...; S_DIM_N S_ID_N;
class vd_ext_COUNT {

  public:
    cell_base* cb;
    apf::Mesh2* m;
    vd_entlist* e_list;
    field_calc* f_calc;

    std::vector<std::vector<int> > c_list;

    double t_curr;

    // The connectivity of the 0-strata to the 3-strata.
    std::string OUTPUTCOUNT;

    void parse_opt(std::vector<std::string> &opts);
    void write_csv();

    vd_ext_COUNT();
    vd_ext_COUNT(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in);
};

// TODO implement this a a list of triplets:
// Given a triplet of strata (S0, S10, S11) a joint 0stratum and two 1strata  
// joining at the stratum) extract the interior angle. Assumes single 
// connection for each joining strata at the mesh level.
// ANGLE1; <ANGLE CSV FILE>; S0; S10; S11
class vd_ext_ANGLE1 {

  public:
    cell_base* cb;
    apf::Mesh2* m;
    vd_entlist* e_list;
    field_calc* f_calc;

    int s0;
    int s10;
    int s11;

    double angle;
    double t_curr;

    // The connectivity of the 0-strata to the 3-strata.
    std::string OUTPUTANGLE;

    void parse_opt(std::vector<std::string> &opts);
    void write_csv();

    vd_ext_ANGLE1();
    vd_ext_ANGLE1(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in);
};

// TODO implement this a a list of quadruplets:
// Given quadruplet of strata (S1, S3, S20, S21) a joint 1stratum, a 3stratum
// and two 2strata joining at the 1stratum with the interior as the 3stratum) 
// extract the interior angle. Assumes single 
// connection for each joining strata at the mesh level.
// ANGLE1; <ANGLE CSV FILE>; S1; S3; S20; S21
class vd_ext_ANGLE2 {

  public:
    cell_base* cb;
    apf::Mesh2* m;
    vd_entlist* e_list;
    field_calc* f_calc;

    int s1;
    int s3;
    int s20;
    int s21;

    double angle;
    double t_curr;

    // The connectivity of the 0-strata to the 3-strata.
    std::string OUTPUTANGLE;

    void parse_opt(std::vector<std::string> &opts);
    void write_csv();

    vd_ext_ANGLE2();
    vd_ext_ANGLE2(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in);
};

////////////////////////////////////////////////
// EXT_ANGLE_CS
// Calculate the exterior angle on the given cross-section plane, outside given 3-stratum between two bounding 1-/2-strata joined at a given 0-/1-stratum entity. 
// If the joint stratum is a 0-stratum, shift the plane to lie on the vertex.
// If one of the bounding strata is a 1-stratum and the 1-stratum doesn't exactly 
// lie in plane, find the 2-stratum entity that intersects the plane and calculate
// the angle using the 2-stratum.
class vd_ext_ANGLE_CS {

  public:
    cell_base* cb;
    apf::Mesh2* m;
    vd_entlist* e_list;
    field_calc* f_calc;

    apf::MeshEntity* e_b0;
    apf::MeshEntity* e_b1;
    apf::MeshEntity* e_j;

    int s3;
    // Stratum ids (apf style) of the bounding and joining strata.
    int sj;
    int sb0;
    int sb1;

    // Dimensions of the bounding and joining strata.
    int db0;
    int db1;
    int dj;

    apf::Vector3 P;
    apf::Vector3 N;

    double angle;
    double t_curr;
    double ang_tol;
    double dist_tol;

    // The connectivity of the 0-strata to the 3-strata.
    std::string OUTPUTANGLE;

    apf::MeshEntity* find_ent(std::vector<apf::MeshEntity*> &ents, apf::Vector3 &l_dir, apf::Vector3 &l_pos);
    bool pl_int_tri_bound_line(apf::Mesh2* m, apf::MeshEntity* tet, 
                                    apf::Vector3 p_pos, apf::Vector3 p_norm, 
                                    apf::Vector3 &l_dir, apf::Vector3 &l_pos);

    void parse_opt(std::vector<std::string> &opts);
    void write_csv();

    vd_ext_ANGLE_CS();
    vd_ext_ANGLE_CS(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in);
};


// Calculate the vertices per unit area.
double ext_ppa(apf::Mesh2* m, vd_entlist* e_list, std::vector<int>& s2_id);

void write_off(const char* fname, apf::Mesh2* m, vd_entlist* e_list, std::vector<int>& s2_id);

class vd_ext_HAUSDORFF : protected vd_extractor {
  private:
    // File name of the csv:
    std::string OUTPUTHAUSDORFF;
    // off file names of the reference and current mesh:
    std::string REF_OFF;
    std::string CURR_OFF;
    std::vector<int> s2_id;
    // Number of points per unit area for CGAL program:
    double nppa;

    void parse_opt(std::vector<std::string> &opts);
  public:
    vd_ext_HAUSDORFF(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc, 
                               std::vector<std::string> & opts_in,
                               double t_in);
    void write_csv();
};

class vd_ext_CRV : protected vd_extractor {
  private:
    std::string OUTPUTCRV;
    apf::Field* rad_field;
    apf::Field* curv_field;
    apf::Field* vel_r_field;
    apf::Field* vel_field;
    apf::Field* t_nbr_field;
    apf::Field* area_field;
    int cell_2;
    int cell_3;
    apf::Vector3 o;

    void parse_opt(std::vector<std::string> &opts);
  public:
    vd_ext_CRV(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc, 
                               std::vector<std::string> & opts_in,
                               double t_in);
    void write_csv();
};

class vd_ext_POS_D : protected vd_extractor {
  private:
    std::string OUTPUTPOS;
    apf::Field* pos_field;
    apf::Field* vel_field;
    apf::Field* e_nbr_field;
    int cell_dim;
    int cell_id;
    apf::Vector3 o;
    apf::Vector3 dir;

    void parse_opt(std::vector<std::string> &opts);
  public:
    vd_ext_POS_D(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc, 
                               std::vector<std::string> & opts_in,
                               double t_in);
    void write_csv();
};

class vd_ext_VEL_D : protected vd_extractor {
  private:
    std::string OUTPUTVEL;
    apf::Field* vel_s_field;
    apf::Field* vel_field;
    apf::Field* e_nbr_field;
    int cell_dim;
    int cell_id;
    apf::Vector3 dir;

    void parse_opt(std::vector<std::string> &opts);
  public:
    vd_ext_VEL_D(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc, 
                               std::vector<std::string> & opts_in,
                               double t_in);
    void write_csv();
};

// Extract the length/area/volume of cells.
class vd_ext_MEAS : protected vd_extractor {
  private:
    std::string OUTPUTMEAS;
    std::vector<std::pair<int, int> > cells;

    void parse_opt(std::vector<std::string> &opts);
  public:
    vd_ext_MEAS(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc, 
                               std::vector<std::string> & opts_in,
                               double t_in);
    void write_csv();
};

// Extract the total area of the boundaries.
class vd_ext_AREA_T : protected vd_extractor {
  private:
    std::string OUTPUTMEAS;

    void parse_opt(std::vector<std::string> &opts);
  public:
    vd_ext_AREA_T(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc, 
                               std::vector<std::string> & opts_in,
                               double t_in);
    void write_csv();
};

class vd_ext_POS : protected vd_extractor {
  private:
    std::string OUTPUTPOS;
    apf::Field* e_nbr_field;
    int cell_dim;
    int cell_id;
    apf::Vector3 o;

    void parse_opt(std::vector<std::string> &opts);
  public:
    vd_ext_POS(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc, 
                               std::vector<std::string> & opts_in,
                               double t_in);
    void write_csv();
};

class vd_ext_VEL : protected vd_extractor {
  private:
    std::string OUTPUTVEL;
    apf::Field* vel_field;
    apf::Field* e_nbr_field;
    int cell_dim;
    int cell_id;

    void parse_opt(std::vector<std::string> &opts);
  public:
    vd_ext_VEL(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc, 
                               std::vector<std::string> & opts_in,
                               double t_in);
    void write_csv();
};

// Given a vector of strings as extraction options specified as a string
// "EXT_TYPE; OUTFILE; OPT1; OPT2; OPT3;..."
// Parse the input to remove the EXT_TYPE and pass the rest of the input to an
// appropriate vd_extractor. Extract the information for all entries in vector. 
class vd_ext_inf {
  private:
    apf::Mesh2* m;
    cell_base* cb;
    vd_entlist* e_list;
    field_calc* f_calc;
    double t_curr;

    std::vector<std::vector<std::string> > opts;

    // Split opt_curr by separator ";". Call the appropriate vd_extractor with
    // splitted options.
    void write_csv();
  public:
    vd_ext_inf(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* e_list_in, 
                              field_calc* f_calc_in,    
                              std::string &opts_file, double t_in);
    ~vd_ext_inf();
};

/*
// Process the information extracted from the mesh.
// Uses extractors. Can be used to derive information such as the velocity at
// the the point satisfying a condition (e.g. at the maximum absolute mean curvature, center of the cell, or the mean velocity of the cell).
// Integral of the gaussian curvature, integral of total cell length, area, volume... So given a field (position, surface energy density) integrate using finite element basis functions (apf interface)
class vd_processor {
  public:
    void proc_field() {
      ext_mesh e_msh;
      e_msh.calc_field();

      // Average can have an integral description, as well
      if(Integrate) {
      }
    }

    vd_processor();
    ~vd_processor();
};

class vd_field_int : public Integrator
{
  public:
    double field_measure(MeshElement* e, char* field_name);
    apf::Vector3 field_measure(MeshElement* e, char* field_name);
apf::getValueType(field_curr);

    vd_field_int(int order):Integrator(order),m(0) {}
    void atPoint(Vector3 const&, double w, double f_val) {
      m += w*f_val;
    }
    void atPoint(Vector3 const&, double w, apf::Vector3 f_val) {
      v += w*f_val;
    }
    double m;
    apf::Vector3 v;
};

// Adapted from Integrator.process and measurer
double vd_field_int::field_measure_db(MeshElement* e, char* field_name) {
  apf::Field* field_curr = e_lens.m->findField(field_name);
  assert(field_curr);
  assert(apf::getValueType(field_curr) == apf::SCALAR);
  w = 0;

  this->inElement(e);
  int np = countIntPoints(e,this->order);

  for (int p=0; p < np; ++p)  {
    ipnode = p;
    apf::Vector3 point;
    getIntPoint(e,this->order,p,point);
    double w = getIntWeight(e,this->order,p);
    double f_val = apf::getScalar(field_curr, e, 0);

    this->atPoint(point,w,f_val);
  }
  this->outElement();
  return m;
}

// Adapted from Integrator.process and measurer
apf::Vector3 vd_field_int::field_measure_v(MeshElement* e, char* field_name) {
  apf::Field* field_curr = e_lens.m->findField(field_name);
  assert(field_curr);
  assert(apf::getValueType(field_curr) == apf::VECTOR);

  this->inElement(e);
  int np = countIntPoints(e,this->order);

  double z[3] = {0,0,0};
  v.fromArray(z);

  for (int p=0; p < np; ++p)  {
    ipnode = p;
    apf::Vector3 point;
    getIntPoint(e,this->order,p,point);
    double w = getIntWeight(e,this->order,p);
    apf::Vector3 f_val;
    apf::getVector(field_curr, e, 0, f_val);

    this->atPoint(point,w,f_val);
  }
  this->outElement();
  return v;
}
*/
#endif
