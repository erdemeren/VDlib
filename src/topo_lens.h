#ifndef TOPO_LENS_H
#define TOPO_LENS_H

#include <vector>

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

#include "topo_disp.h"

// TODO implementation of solution transfer for attached fields. Either, 
// implement a new solution transfer object by using the method implemented in
// SCOREC or try to use the implemented functions.

// vd_lens
// Lens object contains information about entities around a given edge. It is 
// designed to be used with a collapsing edge, member of a 1cell topology. 
// This edge is connected to vertices, at least one of which is a member of a
// 0cell. 
// It keeps entity and topological membership lists for merged edges and 
// removed and merged surfaces, collapsed elements.
// It also keeps modified cells and their level of modification:
//   - Complete removal, 
//   - Complete merge with another cell, 
//   - Pinch merge with another cell.

// It can also be connected to edge collapse operations.
class vd_lens{
  private:
    apf::Mesh2* m;
    gmi_model* mdl;

    struct cell_base* c_base;
    field_calc* f_calc;

    // Flag to write VTK files after collapse.
    bool save_vtk;
    std::string vtk_name;
    std::string smb_name;

    // Collapsing edge;
    apf::MeshEntity* edge_col;
    apf::MeshEntity* vert_ctr;
    apf::ModelEntity* vert_ctr_em;

    // New edges:
    std::pair<apf::MeshEntity*,apf::MeshEntity*> edge_new;

    int type_vertex;
    int tag_vertex;

    std::vector<apf::MeshEntity*> vert_mer;
    std::vector<apf::MeshEntity*> surf_col;
    std::vector<apf::MeshEntity*> elem_col;

    // Middle point of all vertices belonging to the cell to collapse.
    apf::Vector3 midpoint;

    // Rebuild the merging entities for the central disc. 
    apf::MeshEntity* build_repl(apf::MeshEntity* ent_in, apf::ModelEntity* mdl_in); 
    // Rebuild the splitting entities for the central disc. 
    apf::MeshEntity* build_repl(apf::MeshEntity* ent_in, apf::MeshEntity* v_repl);
    void l_clear();

    // Copy constructor
    vd_lens& operator=( const vd_lens& other );
    vd_lens(const vd_lens& that);

    // Go over the collapsing edges. Get the entities around them. For merging 
    // entities if one is collapsing, the other should as well.
    // If they are not collapsing, assign them merging and find the lowest 
    // cell dim cell one belongs to. If they are same dim, the cell should be 
    // the same.
    void split_ent();
  public:

    void split_lens(apf::MeshEntity* v_cls = NULL);
    void load_edge(apf::MeshEntity* edge_in);

    vd_lens(apf::Mesh2* msh, struct cell_base* c, field_calc* f_in = NULL);

    std::pair<int, int> get_cell();

    void set_field_calc(field_calc* f_in);

    std::pair<apf::MeshEntity*, apf::MeshEntity*> get_edge_new();

    void set_files(const char* vtkFile, const char* meshFile);
    void turn_save_vtk(bool onoff);
    apf::MeshEntity* get_vert_ctr();
    ~vd_lens();
};

// Given an edge, split it sequentially to yield one edge of the target length and another edge. Split the other edge until no edge longer than the target remains.
void split_edge_target_len(apf::Mesh2* m, cell_base* cb, field_calc* f_calc, apf::MeshEntity* edge, double len, int max_split = 4);

void split_edge_target_nbr(apf::Mesh2* m, cell_base* cb, field_calc* f_calc, apf::MeshEntity* edge, int nbr, int max_split = 4);

void vd_split_tet(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* tet,  field_calc* f_calc);

void vd_split_tet_vert(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* vert,  field_calc* f_calc);

// Splitter for a triangle (operations within a triangular bipyramid).
class vd_bipy{
  private:
    apf::Mesh2* m;
    gmi_model* mdl;

    struct cell_base* c_base;
    field_calc* f_calc;

    // Flag to write VTK files after collapse.
    bool save_vtk;
    std::string vtk_name;
    std::string smb_name;

    // Remaining entities
    std::vector<apf::MeshEntity*> vert_surr;
    apf::MeshEntity* vert_t;
    apf::MeshEntity* vert_b;

    // Splitting entities;
    apf::MeshEntity* tri_col;
    std::vector<apf::MeshEntity* > tets;

    // New entities;
    apf::MeshEntity* vert_ctr;
    apf::ModelEntity* vert_ctr_em;

    int type_vertex;
    int tag_vertex;

    // Edges and triangles from triangle:
    std::vector<apf::MeshEntity*> edge_tri;
    std::vector<apf::MeshEntity*> tri_tri;
    // Edges from tets:
    std::vector<apf::MeshEntity*> edges;

    // Middle point of all vertices belonging to the cell to collapse.
    apf::Vector3 midpoint;

    apf::MeshEntity* build_repl(); 
    // Rebuild the splitting entities for the central disc. 
    void l_clear();

    // Copy constructor
    vd_bipy& operator=( const vd_bipy& other );
    vd_bipy(const vd_bipy& that);

    // Go over the collapsing edges. Get the entities around them. For merging 
    // entities if one is collapsing, the other should as well.
    // If they are not collapsing, assign them merging and find the lowest 
    // cell dim cell one belongs to. If they are same dim, the cell should be 
    // the same.
    void split_ent();
  public:


    void split_bipy(apf::MeshEntity* v_cls = NULL);
    void load_tri(apf::MeshEntity* tri_in);

    vd_bipy(apf::Mesh2* msh, struct cell_base* c, field_calc* f_in = NULL);

    std::pair<int, int> get_cell();

    void set_field_calc(field_calc* f_in);

    void set_files(const char* vtkFile, const char* meshFile);
    void turn_save_vtk(bool onoff);
    apf::MeshEntity* get_vert_ctr();
    ~vd_bipy();
};
/*
// Splitter for a triangle and collapse of tet bounded by the triangle.
// (leads to collapse of the tet, recreation of the triangle by three triangles 
// and split of a triangular pyramid).
class vd_bimo{
  private:
    apf::Mesh2* m;
    gmi_model* mdl;

    struct cell_base* c_base;
    field_calc* f_calc;

    // Flag to write VTK files after collapse.
    bool save_vtk;
    std::string vtk_name;
    std::string smb_name;

    // Remaining entities
    std::vector<apf::MeshEntity*> vert_surr;
    apf::MeshEntity* vert_t;
    apf::MeshEntity* vert_b;

    // Splitting entities;
    apf::MeshEntity* tri_col;
    // Collapsing tet:
    apf::MeshEntity* tet_col;
    // All tets:
    std::vector<apf::MeshEntity* > tets;

    // New entities;
    apf::MeshEntity* vert_ctr;
    apf::ModelEntity* vert_ctr_em;

    int type_vertex;
    int tag_vertex;

    // Edges and triangles from triangle:
    std::vector<apf::MeshEntity*> edge_tri;
    std::vector<apf::MeshEntity*> tri_tri;
    // Edges from tets:
    std::vector<apf::MeshEntity*> edges;

    // Middle point of all vertices belonging to the cell to collapse.
    apf::Vector3 midpoint;

    apf::MeshEntity* build_repl(); 
    // Rebuild the splitting entities for the central disc. 
    void l_clear();

    // Copy constructor
    vd_bipy& operator=( const vd_bipy& other );
    vd_bipy(const vd_bipy& that);

    // Go over the collapsing edges. Get the entities around them. For merging 
    // entities if one is collapsing, the other should as well.
    // If they are not collapsing, assign them merging and find the lowest 
    // cell dim cell one belongs to. If they are same dim, the cell should be 
    // the same.
    void split_ent();
  public:


    void split_bipy(apf::MeshEntity* v_cls = NULL);
    void load_tri(apf::MeshEntity* tri_in, apf::MeshEntity* tet);

    vd_bipy(apf::Mesh2* msh, struct cell_base* c, field_calc* f_in = NULL);

    std::pair<int, int> get_cell();

    void set_field_calc(field_calc* f_in);

    void set_files(const char* vtkFile, const char* meshFile);
    void turn_save_vtk(bool onoff);
    apf::MeshEntity* get_vert_ctr();
    ~vd_bipy();
};
*/

// Splitter for a tet.
class vd_sp_tet{
  private:
    apf::Mesh2* m;
    gmi_model* mdl;

    struct cell_base* c_base;
    field_calc* f_calc;

    // Flag to write VTK files after collapse.
    bool save_vtk;
    std::string vtk_name;
    std::string smb_name;

    // Remaining entities
    std::vector<apf::MeshEntity*> vert_surr;

    // Splitting entities;
    apf::MeshEntity* tet_col;

    // New entities;
    apf::MeshEntity* vert_ctr;
    apf::ModelEntity* vert_ctr_em;

    int type_vertex;
    int tag_vertex;

    // Middle point of all vertices belonging to the cell to collapse.
    apf::Vector3 midpoint;

    // Rebuild the splitting entities for the central disc. 
    void l_clear();

    // Copy constructor
    vd_sp_tet& operator=( const vd_sp_tet& other );
    vd_sp_tet(const vd_sp_tet& that);

    // Go over the collapsing edges. Get the entities around them. For merging 
    // entities if one is collapsing, the other should as well.
    // If they are not collapsing, assign them merging and find the lowest 
    // cell dim cell one belongs to. If they are same dim, the cell should be 
    // the same.
    void split_ent();
  public:

    void split_tet(apf::MeshEntity* v_cls = NULL);
    void load_tet(apf::MeshEntity* tet_in);

    vd_sp_tet(apf::Mesh2* msh, struct cell_base* c, field_calc* f_in = NULL);

    std::pair<int, int> get_cell();
    void set_files(const char* vtkFile, const char* meshFile);

    void set_field_calc(field_calc* f_in);

    void turn_save_vtk(bool onoff);
    apf::MeshEntity* get_vert_ctr();
    ~vd_sp_tet();
};

void vd_split_edge(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity edge,  field_calc* f_calc = NULL);

// Given a triangle and a vertex belonging to the same 2-stratum, find the longest
// edge belonging to the same 2-stratum and not bounding the triangle. Split that
// edge if it exists.
void vd_split_edge_2c(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* tri, apf::MeshEntity* vert,  field_calc* f_calc = NULL);


// Collapsers:
// Edge collapser which assumes no topological change.
class vd_lens_col{
  private:
    apf::Mesh2* m;
    gmi_model* mdl;

    struct cell_base* c_base;
    field_calc* f_calc;

    // Flag to write VTK files after collapse.
    bool save_vtk;
    std::string vtk_name;
    std::string smb_name;

    // Collapsing edge;
    apf::MeshEntity* edge_col;

    // The remaining vertex;
    apf::MeshEntity* vert_ctr;
    apf::MeshEntity* vert_other;
    apf::ModelEntity* vert_ctr_em;

    int type_vertex;
    int tag_vertex;

    std::pair<bool, bool> valid;
    // If prs_bound is true, position the new vertex at the old position of the 
    // vertex with the smallest topological membership.
    bool prs_bound;

    std::vector<std::vector<apf::MeshEntity*> > edge_sur_vert;
    std::vector<std::vector<apf::MeshEntity*> > surf_sur_vert;
    std::vector<std::vector<apf::MeshEntity*> > elem_sur_vert;
    std::vector<apf::ModelEntity*> elem_sur_mdl;
    std::vector<apf::ModelEntity*> surf_sur_mdl;
    std::vector<apf::ModelEntity*> edge_sur_mdl;

    std::vector<std::vector<apf::MeshEntity*> > edge_mer_vert;
    std::vector<std::vector<apf::MeshEntity*> > surf_mer_vert;
    std::vector<apf::ModelEntity*> edge_mer_mdl;
    std::vector<apf::ModelEntity*> surf_mer_mdl;


    // The entities to be replaced adjacent to the vertex to be deleted.
    std::vector<apf::MeshEntity*> edge_mer;
    std::vector<apf::MeshEntity*> surf_mer;

    std::vector<apf::MeshEntity*> edge_sur;
    std::vector<apf::MeshEntity*> surf_sur;
    std::vector<apf::MeshEntity*> elem_sur;
    std::vector<apf::MeshEntity*> elem_sur_new;

    std::vector<apf::MeshEntity*> vert_mer;

    std::vector<apf::MeshEntity*> surf_col;
    std::vector<apf::MeshEntity*> elem_col;

    // Middle point of all vertices belonging to the cell to collapse.
    apf::Vector3 midpoint;

    // Rebuild the surrounding entities. 
    apf::MeshEntity* build_repl(std::vector<apf::MeshEntity*> &vert_in, apf::ModelEntity* mdl_in);

    apf::MeshEntity* build_repl(apf::MeshEntity* ent_in, apf::ModelEntity* mdl_in, apf::MeshEntity* v_repl);
    void l_clear();

    // Copy constructor
    vd_lens_col& operator=( const vd_lens_col& other );
    vd_lens_col(const vd_lens_col& that);

    // Check the stratum and shell memberships of the merging entities.
    // Get the stratum membership and position of the final vertex.
    // Check and fix unintended merging entities by splitting one of the 
    // associated edges.
    // Remove the collapsing entities, recreate the surrounding entities 
    // adjacent to the vertex to be deleted by replacing the vertex with the 
    // remaining vertex. 
    void col_ent();

    bool chk_ents(apf::MeshEntity* ec, apf::MeshEntity* e1, apf::MeshEntity* e2);
    bool chk_cell();
    bool chk_inv();
    void chk_adj();
    void fix_adj();

    void clear();

    void add_merge(apf::MeshEntity* e1, apf::MeshEntity* e2);
    void add_surr(apf::MeshEntity* e);
    void collect_merge();
  public:

    std::pair<bool, bool> col_lens(apf::MeshEntity* v_cls = NULL);
    std::pair<bool, bool> load_edge(apf::MeshEntity* edge_in);

    void assign_surr_map(std::map<apf::MeshEntity*, apf::MeshEntity*>& map_in);
    void assign_col_map(std::map<apf::MeshEntity*, bool>& map_in);

    void preserve_bound(bool flag_in);
    vd_lens_col(apf::Mesh2* msh, struct cell_base* c, field_calc* f_in = NULL);

    std::pair<int, int> get_cell();

    void set_field_calc(field_calc* f_in);

    void set_files(const char* vtkFile, const char* meshFile);
    void turn_save_vtk(bool onoff);
    apf::MeshEntity* get_vert_ctr();
    ~vd_lens_col();
};

std::pair<bool, bool> vd_col_edge(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* edge, std::map<apf::MeshEntity*, bool>& tet_surr_map, std::map<apf::MeshEntity*, bool>& tet_col_map, field_calc* f_calc = NULL, bool preserv_bound = false);

std::pair<bool, bool> vd_col_edge(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* edge,  field_calc* f_calc = NULL, bool preserv_bound = false);

std::pair<bool, bool> vd_col_edge_tet(apf::Mesh2* m, cell_base* c_base, apf::MeshEntity* tet, std::map<apf::MeshEntity*, apf::MeshEntity*>& tet_surr_map, std::map<apf::MeshEntity*, bool>& tet_col_map, field_calc* f_calc = NULL, bool preserv_bound = false);

#endif

