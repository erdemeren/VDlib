#ifndef TOPO_TMESH_H
#define TOPO_TMESH_H

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

// The cell detector object that works with expansions around 0cells with
// disjoint entity sets belonging to the same cell. This extracts the disjoint
// cells and tags them uniquely and extracts a 2/3-cell adjacency graph. 
class vd_cell_det {
  private:
    int c0_tag;

    int c_dim; 
    int d_ent;
    int id;

    bool calc_ext;
    bool calc_corner;

    bool ext_0cell;
    bool ext_corner;
    std::map<int, bool> ext_2cell;
    int ext_sz;
    int c3_nbr;

    apf::Mesh2* m;
    struct cell_base* cb;
    cell_graph cg;
    vd_entlist_v ent_list;

    ngon_gmi ng_local;

    // The cellbase structure for the disjoint entities. It is used to 
    // consistently tag entity sets between trials.
    struct cell_base* cb_dj;
    // The flag that checks to write to or read from the cb_dj when tagging 
    // the entity sets. 
    // True: fill the cb_dj and id the entity sets and their
    // False: use the ids and connectivity in the cb_dj. Assert the number of
    // cell entities is equal or larger than the number of cells bounded by the
    // actual 0cell.
    bool cb_write_flag;
    // The vertex disjoint 1cell id map. It is assigned from the vd_edisc 
    // object, for each trial and main meshes.
    std::map<apf::MeshEntity*, int> vert_c1_map;
    // The real c3 adjacencies of local c2s. Used in the case a number of c2s 
    // are bounded by the same c1s.
    std::vector<int> c2_c3real1;
    std::vector<int> c2_c3real2;

    apf::MeshEntity* v_0cell;

    apf::ModelEntity* mdl_curr;

    std::vector<std::vector<int > > cells;
    std::vector<std::vector<std::vector<int > > > adj;

    std::vector<apf::MeshEntity* > e_2c;
    std::map<int, apf::MeshEntity*> c3_edge;

    std::map<int, apf::Vector3> c2_dir;
    std::map<int, apf::Vector3> c3_dir;

    std::vector< circuit > circ_tup;
    std::vector< cell_path > path;

    std::vector<int> ent_sz;

    // [dim]<mdl, id> Disjoint segments of dim-cells
    std::vector<std::map<int, apf::ModelEntity*> > burn_mdl;
    std::map<apf::MeshEntity*, bool> burn_ent;
    // Same dimensional cell membership of the mesh entity. Checked and defined
    // for mesh entities belonging to the same dimensional cells.
    std::map<apf::MeshEntity*, int> ent_cell;

    std::vector<std::vector<apf::MeshEntity* > > ent_ofire;

    void burn_cell(int c_dim, int id, apf::MeshEntity* ent);

    // Copy constructor and assignment private:
    vd_cell_det(const vd_cell_det &obj);
    vd_cell_det& operator=( const vd_cell_det& obj );

  public:
    vd_cell_det();
    vd_cell_det(apf::Mesh2* m_in, struct cell_base* cb_in,
                                              apf::MeshEntity* v_0cell_in);
    void dummy_func_stop();
    void load_graph();

    // Compare the cell adjacencies with the connectivity of the disjoint sets of
    // entities. Used as a sanity check. 
    void compare_cells();

    // Using the ids of the bounding disjoint cells, return the id of dim-cell
    int find_cell_dj(int dim, apf::MeshEntity* ent, std::vector<int> bound);

    // Going over disjoint sets of dim-dimensional entities, assign their lower 
    // dimensional adjacencies bounded by the central vertex of the same disjoint 
    // set of entities with the same id. 
    void burn_lower();

    // TODO cg seems to assume the cells list to be ordered. This is a bug, but 
    // sort it out for now by reordering the cells before feeding into the 
    // cell_graph object.
    void reorder_cells();

    // Return the disjoint set id of the entity.
    int get_id(apf::MeshEntity* ent);

    // Burn the disjoint cells, while populating the disjoint cell_base.
    void burn_new();
    // Burn the disjoint cells, using the supplied 1cell vertex map and the
    // existing disjoint cell_base.
    void burn_from_cb();
    // Burn the disjoint cell to find the bounding disjoint cells.
    std::vector<int> burn_bound(int c_dim, apf::MeshEntity* ent);

    ~vd_cell_det();

    void get_circ_topo_gmi(std::vector<int>* c2_list, 
                                  std::vector<int>* c3_list, circuit* circ);
    void get_circ_topo_dis_gmi(std::vector<ent_conn>* cs, 
                                                circuit* circ, int g_id);
    void clear_path();
    void clear_mesh();
    void clear_topo();
    void clear();

    void reload(apf::Mesh2* m_in, struct cell_base* cb_in,
                                              apf::MeshEntity* v_0cell_in);

    std::vector<apf::MeshEntity*> get_dim_ent(int dim, int id);
    int get_3c_nbr();
    int get_3c_nbr_ext();
    int get_2c_nbr();
    // Return the edges associated with the disjoint 2cell triangles.
    void coll_2c_edges();
    void coll_3c_edges();
    std::vector<apf::MeshEntity*> get_c2_edges();
    std::vector<apf::MeshEntity*> get_c3_edges();
    void set_c2_edge(int ent_cell_curr, apf::MeshEntity* e_curr);
    void set_c3_edge(int ent_cell_curr, apf::MeshEntity* e_curr);

    void get_path_ngon(int c3_cp, ngon_gmi* ng);
    void conv_path_ngon_gmi(ngon_gmi* ng);
    int conv_path_2c_gmi(int n_id);
    int conv_path_3c_gmi(int n_id);
    void shift_ngon_gmi(ngon_gmi* ng);

    int get_circ_type(int c_in);
    bool get_ext();

    void set_calc_ext(bool fix);
    void set_calc_corner(bool fix);
    void set_proj(PROJ_TYPE PROJ);

    int get_circ_sz();
    int get_path_sz();
    circuit * get_circ(int t_curr);

    // Cell list treatment for corner 0cells
    void treat_cells();

    void update_path();

    // Return bounding disjoint strata of given stratum.
    std::vector<int> get_bound_ids(int dim, int id);
    // Extract the slices and paths from the current cell_graph object, for a 
    // given ngon or circuit.
    void set_ng_couple(int p_id);
    void find_slice_ng(int ng_id,
            std::vector<std::pair<std::vector<int >, std::vector<int > > >*
                                                              path_cells,
            std::vector< std::pair< std::pair<int,int>, 
                         std::vector<std::vector<int > > > >* slice_cells);

    void find_slice_circ(int circ_id, 
                std::vector<std::pair<std::vector<int >, std::vector<int > > >*
                                                                  path_cells,
                std::vector< std::pair< std::pair<int,int>, 
                             std::vector<std::vector<int > > > >* slice_cells);
    // ******************  CONSISTENCY OPS ************************** //
    apf::ModelEntity* get_mdl(int dim_in, int id_in);

    // Obtain the relative directions of the triangles/tets of a c2/c3.
    // Not used at the moment.
    void coll_c2_dir();
    void coll_c3_dir();

    // Get the stored relative directions of the triangles/tets of a c2/c3.
    // Not used at the moment.
    apf::Vector3 get_c2_dir(int c_id);
    apf::Vector3 get_c3_dir(int c_id);

    // To make the local indices of the disjoint cells consistent across all 
    // insertions, get the local id maps of c1 from the active vd_edisc object.
    void set_vert_c1_map(std::map<apf::MeshEntity*, int> map_in);
    std::map<apf::MeshEntity*, int> get_vert_c1_map();

    // ********************  FLAGS  **************************** //
    // Set if the existing local indices are to be used or a new local  
    // cell_base object is to be initialized.
    void set_cb_flag(bool on_off);

};

// Check if the cell is spurious using mesh information. This is useful when
// the regularity assumption of the cell complex does not hold, when upper dim.
// adjacencies don't exist.
// A cell is spurious if it is not bounding any upper dim cells or bounding 
// two upper dim cells bounding the same cells.
// This is only used when a newly created cell does not have any immediate 
// upper dim adjacencies after an insertion(1cell) or collapse(0cell), so that
// cell complex cannot be used.
// Assert that the lowest high dim. bounded cell is the only cell in its 
// dimension (So a single 2cell or if non a single 3cell).

class vd_spur_mesh {
  private:

    apf::MeshEntity* ent_c;
    int d_ent;
    apf::ModelEntity* mdl;

    int cell_repl;
    int c_nbr;

    apf::Mesh2* m;
    struct cell_base* cb;
    vd_entlist_v ent_list;


    void get_vert();
  public:
    vd_spur_mesh();
    vd_spur_mesh(apf::Mesh2* m_in, cell_base* cb_in,
                                              apf::MeshEntity* ent_c_in);

    void reload(apf::Mesh2* m_in, cell_base* cb_in,
                                              apf::MeshEntity* ent_c_in);

    // Return the lowest high dim. bounded cell. first = -1 if not spurious. 
    std::pair<int,int> get_cell_repl();

    //void fix_spur();
};

// The 3cell detector. Detects disjoint groups of entities belonging to 3cells
// TODO also check for the homology of the disjoint sets. Currently handles 
// only ball-like 3cells, doesn't support torus-like 3cells.
// Can be incorporated into vd_cell_det but that assumes a single edge per 
// 3cell at the moment.
class vd_3c_det {
  private:
    int c0_tag;

    int c_dim; 
    int d_ent;
    int id;

    int c3_nbr;

    apf::Mesh2* m;
    struct cell_base* cb;
    vd_entlist_v ent_list;

    // The cellbase structure for the disjoint entities. It is used to 
    // consistently tag entity sets between trials.
    struct cell_base* cb_dj;

    apf::MeshEntity* v_0cell;

    apf::ModelEntity* mdl_curr;

    std::vector<std::vector<int > > cells;
    std::vector<std::vector<std::vector<int > > > adj;

    std::vector<int> ent_sz;

    // The tetrahedra belonging to disjoint segments of 3strata.
    std::vector<std::vector<apf::MeshEntity* > > c3_ents;

    // [dim]<mdl, id> Disjoint segments of dim-cells
    std::vector<std::map<int, apf::ModelEntity*> > burn_mdl;
    std::map<apf::MeshEntity*, bool> burn_ent;
    // Same dimensional cell membership of the mesh entity. Checked and defined
    // for mesh entities belonging to the same dimensional cells.
    std::map<apf::MeshEntity*, int> ent_cell;

    std::vector<std::vector<apf::MeshEntity* > > ent_ofire;

    void burn_cell(int c_dim, int id, apf::MeshEntity* ent);

    // Copy constructor and assignment private:
    vd_3c_det(const vd_3c_det &obj);
    vd_3c_det& operator=( const vd_3c_det& obj );

  public:
    vd_3c_det();
    vd_3c_det(apf::Mesh2* m_in, struct cell_base* cb_in,
                                              apf::MeshEntity* v_0cell_in);
    void load_graph();

    // Using the ids of the bounding disjoint cells, return the id of dim-cell
    int find_cell_dj(int dim, apf::MeshEntity* ent, std::vector<int> bound);
    // Burn the disjoint cells, while populating the disjoint cell_base.
    void burn_new();

    ~vd_3c_det();
    void clear();

    void reload(apf::Mesh2* m_in, struct cell_base* cb_in,
                                              apf::MeshEntity* v_0cell_in);

    std::vector<apf::MeshEntity*> get_dim_ent(int dim, int id);
    // Return bounding disjoint strata of given stratum.
    std::vector<int> get_bound_ids(int dim, int id);

    std::vector<apf::MeshEntity*> get_tet_id(int id);
    int get_3c_nbr();

    // ******************  CONSISTENCY OPS ************************** //
    apf::ModelEntity* get_mdl(int dim_in, int id_in);
    int get_id(apf::MeshEntity* ent);

};

// The cell replacer replaces the cell membership of the entities with the
// new cell. 
class vd_cell_repl {
  private:
    apf::Mesh2* m;

    std::map<apf::ModelEntity*, apf::ModelEntity*> c_map;
    std::map<apf::MeshEntity*, apf::MeshEntity*> e_map;
    std::map<apf::ModelEntity*, bool> repl;
    std::vector<std::vector<apf::MeshEntity*> > dest_ent;

    std::vector<std::vector<std::vector<apf::MeshEntity*> > > down_list;
    std::vector<std::vector<apf::ModelEntity*> > mdl_list;

  public:
    vd_cell_repl(apf::Mesh2* m_in);
    // Update the mesh in the object.
    void update_m(apf::Mesh2* m_in);

    // Sets the cells.
    void set_cell(std::vector<c_elem>* in_cells);
    // Replaces the cells.
    void repl_cell();

    void clear();

};

class outer_cell_gen {
  private:
    apf::Mesh2* m;
    struct cell_base* cb;

    int c3_out;

    // Bounding 2 cell entities
    std::map<apf::ModelEntity*, bool> out;
    std::map<apf::ModelEntity*, apf::MeshEntity*> c0_map;
    std::map<apf::ModelEntity*, apf::MeshEntity*> c1_map;
    std::map<apf::ModelEntity*, apf::MeshEntity*> c2_map;

    std::map<apf::ModelEntity*, bool> repl;
    std::vector<std::vector<apf::MeshEntity*> > dest_ent;

    std::vector<std::vector<std::vector<apf::MeshEntity*> > > down_list;
    std::vector<std::vector<apf::ModelEntity*> > mdl_list;

    void collect_cell();

  public:
    outer_cell_gen(apf::Mesh2* m_in, struct cell_base* cb_in);

    // Find the outer shell by burning algorithm. It could be possible to have
    // holes. In that case, the algorithm will yield multiple shells. Each will
    // have an associated 3cell.
    // We potentially don't need the outer shell to be continuous.
    // Still, given two surfaces to be bounding the same outer cell, can one 
    // check if the surfaces are concave or convex?
    void create_cell();

    // Find the outer shell by burning algorithm. Assume no holes.
    // Create a 2cell vertex on each outer surface and create tets using the 
    // outer cell triangles and the new 2cell verts. 
    // Create outer cell corners, and new 1cell vertices. Call the new 1cells
    // beams.
    // Create fins connecting the 1cell edges to the new cell 1cell centers. 
    // Create triangles between the corner 0cell vertices and new edges 
    // connecting the new 2cell vertices and old 0cell corners, call them 
    // pillars.
    // Create tris on each new 2cell surface for each new 0cell corner.
    // Three at each new 0cell.
    // Pillar vertices, 0cell corner vertex.
    // For all tris between two new pillars, create tets with beam vertices.
    // Create three tets at each corner(one for each side of the pillar.
    // New 0cell corner vertex, pillar vertices, beam vertex.


    // It is not possible to diffrentiate the 0cell corners of the simulating 
    // box by topology alone. Geometric information is needed.
    // One can find the boundary center and outer pointing normals. If all
    // outer pointing normals point away from the shell center, it is convex 
    // and the new 3cell is an outer 3cell.
    // Else if all point towards the cell center, it is convex and the new 
    // 3cell is an inner 3cell, or hole.
    // Otherwise, it is not convex and current algorithm will not work.
    // If it is a hole, skip for now.

    // Create a bounding 2cell vertex for every 2cell on the boundary shell.
    // Create the new 3cell edges, tris, and tets using the new 2cell vertices 
    // and the old 2cell vertices and their bounding vertices, edges and 
    // old 2cell tris.
    // Call the edges connected to old 0cell, and 1cell subpoles,  
    // and the 2cell vertices pillars.
    // In addition, create a 2cell vertex and a 3cell edge for every 0cell on 
    // the boundary. Call these edges connected to old 0cell vertices poles. 
    // Create a vertex and two edges between pole end vertices. Call the vertex
    // the beam vertex.
    // Connect to the previous algorithm.

    // Position the new vertices a constant distance from each connected 2cell
    // vertex. This will make sure the outer boundary takes the shape of the 
    // inner boundary.


    // Actually, this:
    // It is not possible to diffrentiate the 0cell corners of the simulating 
    // box by topology alone. Geometric information is needed.
    // One can find the boundary center and outer pointing normals. If all
    // outer pointing normals point away from the shell center, it is convex 
    // and the new 3cell is an outer 3cell.
    // Else if all point towards the cell center, it is convex and the new 
    // 3cell is an inner 3cell, or hole.
    // Otherwise, it is not convex and current algorithm will not work.
    // If it is a hole, skip for now.

    // Copy the boundary vertices by scaling by a minor amount. Create tets
    // between the old and the new vertices.

    // Extend edges between the old and the new vertices. 
    // Each old edge will have two associated vertices as corner(one bottom 
    // one top). The new tets will be formed by these two vertices, and 
    // Top tet: two vertices, the other top and the new vertex
    // Bottom tet: two vertices, the other bottom and the old vertex

    void create_cell_out();

    // Gets the cells. If no 2cell exists (filled hole) get_2c_out returns -1.
    std::vector<int> get_3c_out();
    std::vector<int> get_2c_out();

    void clear();

};

// A truncated version of the cell_det object to check the validity of the mesh
// in the neighborhood of a vertex. vd_cell_det could also be used, but it 
// assumes in many places that the center vertex is a 0cell vertex. Use this
// for now. 
// TODO It should be possible to have the entities 
// of a stratum stratum to be connected multiple times to the given vertex,
// especially if it is a 0cell.

class vd_cell_val {
  private:
    int cell_tag;
    int cell_dim;

    int c_dim; 
    int d_ent;
    int id;

    bool valid;

    apf::Mesh2* m;
    struct cell_base* cb;
    vd_entlist_v ent_list;

    // The cellbase structure for the disjoint entities. It is used to 
    // consistently tag entity sets between trials.
    struct cell_base* cb_dj;

    apf::MeshEntity* v_ctr;
    apf::ModelEntity* mdl_curr;

    std::vector<std::vector<int > > cells;
    std::vector<std::vector<std::vector<int > > > adj;

    std::vector<int> ent_sz;

    // [dim]<mdl, id> Disjoint segments of dim-cells
    std::vector<std::map<int, apf::ModelEntity*> > burn_mdl;
    std::map<apf::MeshEntity*, bool> burn_ent;
    // Same dimensional cell membership of the mesh entity. Checked and defined
    // for mesh entities belonging to the same dimensional cells.
    std::map<apf::MeshEntity*, int> ent_cell;

    std::vector<std::vector<apf::MeshEntity* > > ent_ofire;

    void burn_cell(int c_dim, int id, apf::MeshEntity* ent);

    // Copy constructor and assignment private:
    vd_cell_val(const vd_cell_val &obj);
    vd_cell_val& operator=( const vd_cell_val& obj );

  public:
    vd_cell_val();
    vd_cell_val(apf::Mesh2* m_in, struct cell_base* cb_in,
                                              apf::MeshEntity* v_0cell_in);
    void dummy_func_stop();
    void load_graph();

    // Using the ids of the bounding disjoint cells, return the id of dim-cell
    int find_cell_dj(int dim, apf::MeshEntity* ent, std::vector<int> bound);
    // Burn the disjoint cells, while populating the disjoint cell_base.
    void burn_new();
    // Burn the disjoint cell to find the bounding disjoint cells.
    std::vector<int> burn_bound(int c_dim, apf::MeshEntity* ent);

    ~vd_cell_val();

    void clear();

    void reload(apf::Mesh2* m_in, struct cell_base* cb_in,
                                              apf::MeshEntity* v_0cell_in);

    std::vector<apf::MeshEntity*> get_dim_ent(int dim, int id);

    // ******************  CONSISTENCY OPS ************************** //
    apf::ModelEntity* get_mdl(int dim_in, int id_in);

    // ********************  CHECKS  **************************** //
    // Check if the edges and triangles around the vertex are allowed by the 
    // topology.
    bool vd_vert_valid();
};


#endif
