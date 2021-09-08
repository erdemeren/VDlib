#ifndef TOPO_EXT_H
#define TOPO_EXT_H

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_rand.h"

// Entity set prototype to keep twice upper or lower adjacencies, triple edges, 
// quadruple points, boundary triangles, etc. Similar to Up.

struct Entity_set
{
public:
  int n;
  apf::MeshEntity* e[6800];
  Entity_set() : n(0) {}
};

/*
struct Entity_set
{
public:
  int n;
  std::vector<apf::MeshEntity*> e;
  Entity_set() : n(0) {e.resize(68000)}
};
*/
void dummy_clear_stop();
void dummy_chk_stop();

// Vector implementation of findin, which is spurious but will help replacing 
// the Entity_set implementation.
int findIn(const std::vector<apf::MeshEntity* >* es, int n, apf::MeshEntity* e);
void copy_ent_set(std::vector<apf::MeshEntity*>* es, apf::Up up);
void copy_ent_set(std::vector<apf::MeshEntity*>* es, 
                                    apf::Downward down, int size);
void copy_ent_set(std::vector<apf::MeshEntity*>* es_copy, 
                  std::vector<apf::MeshEntity*>* es_in);
void init_ent_set(std::vector<apf::MeshEntity*>* es, apf::MeshEntity* e);
void rev_ent_set(std::vector<apf::MeshEntity*>* es);

// Given a set of entities, collect the bounding vertices of the set.
bool vd_vert_set(apf::Mesh2* m, std::vector<apf::MeshEntity* >* es_in, std::vector<apf::MeshEntity* >* es_out);

void print_ent_set(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es);
void find_elem_set_uniq(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es, std::vector<int>* g_list);

void vd_tag_set(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es, const char* set_name);

void vd_find_ent_geom(apf::Mesh2* m, std::vector<apf::MeshEntity* >* ep, 
                                    int cell_id, int dim,  int dim_e = -1);

void vd_find_cell_v(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev, int cell_id, int dim);
void vd_rem_cell(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev, int cell_id, int dim);
void vd_find_bound_v(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev, int cell_id, int dim);

void vd_find_edge_topo(apf::Mesh2* m, 
                std::vector<apf::MeshEntity* >* ep, int cell_id, int dim);

// Given a vertex, return the set of dim dimensional upward entity adjacencies 
// belonging to dim-cell-cell_id.
void vd_find_ent_topo(apf::Mesh2* m, apf::MeshEntity* v, 
                std::vector<apf::MeshEntity* >* ep, int cell_id, int dim);

// Given an entity set, find the entities member of a certain geometry.
// Store in ep*. If no geometry specified, return all, member of any geometry.
void vd_ext_ent_geom(apf::Mesh2* m, 
                    const std::vector<apf::MeshEntity* >* ep_in, 
                    std::vector<apf::MeshEntity* >* ep, int geom = 0);

// Given a mesh object, extract triple edges.
void vd_ext_trip(apf::Mesh2* m, std::vector<apf::MeshEntity* >* es_edge);

// Given a mesh object, extract boundary and domain surfaces.
void vd_ext_bsurf(apf::Mesh2* m, 
                    std::vector<apf::MeshEntity* >* es_surf, 
                    std::vector<apf::MeshEntity* >* es_dom);

// Given a mesh object, and domain surface list extract quadruple points, which
// lie on the interior, not on the domain boundary. For now, assume quadruple
// points indeed are members of at least four different grains.
// PLACEHOLDER for now...
void vd_ext_quad(apf::Mesh2* m, std::vector<apf::MeshEntity* >* es_quad);

void vd_remove_set(std::vector<apf::MeshEntity* >* ep1, 
                     const std::vector<apf::MeshEntity* >* ep2);
void vd_inter_set(std::vector<apf::MeshEntity* >* ep1, 
                      std::vector<apf::MeshEntity* >* ep2);

// Given two up adjacency entity arrays, obtain a list of unique entities.
void vd_merge_set(std::vector<apf::MeshEntity* >* ep1, 
                    std::vector<apf::MeshEntity* >* ep2);

// Given an entity set ep1, find one level higher adjacencies if possible.
// Store in ep*.
void vd_set_up2(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep1, 
                      std::vector<apf::MeshEntity*>* ep);

void vd_set_up(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep1, 
                      std::vector<apf::MeshEntity*>* ep);

// Given a entity, find one level higher adjacencies if possible.
// Store in ep*.
void vd_set_up(apf::Mesh2* m, apf::MeshEntity* e, 
                      std::vector<apf::MeshEntity*>* ep);

// Get the tetrahedra adjacent to the vertices bounding a given set of entities.
void vd_get_tet(apf::Mesh2* m, const std::vector<apf::MeshEntity*> & ep1, 
                      std::vector<apf::MeshEntity*> & ep);

void vd_set_down(apf::Mesh2* m, apf::MeshEntity* e, 
                              std::vector<apf::MeshEntity*>* ep, int level = 1);

// Given an entity set ep_in, find one level lower adjacencies if possible.
// Store in ep*.
void vd_set_down(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in, 
                      std::vector<apf::MeshEntity*>* ep, int level = 1);

void vd_set_down2(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in, 
                      std::vector<apf::MeshEntity*>* ep, int level = 1);

void vd_chk_edge_grain(apf::Mesh2* m, 
                        const std::vector<apf::MeshEntity* >* ep_in, 
                        std::vector<apf::MeshEntity* >* ep, int geom);

// Given an surface entity set ep_in, extract those, which rest on the grain 
// surface.
// Store in ep*.
// LOOK write functions that do something simple, and only that!
void vd_chk_surf_grain(apf::Mesh2* m, 
                        const std::vector<apf::MeshEntity* >* ep_in, 
                        std::vector<apf::MeshEntity* >* ep, int geom);

// Given an entity set ep_in, extract those, one of d dim upper adjacency of 
// which belongs to a certain geometry of that adjacencies own dimension.
// (i.e. a surface entity set, whose adjacent elements belong to geometry n.)
// Store in ep*.
// LOOK write functions that do something simple, and only that!
// So write a function, which extracts edges on a grain surface. 
void vd_chk_set_geom(apf::Mesh2* m, 
                      const std::vector<apf::MeshEntity* >* ep_in, 
                      std::vector<apf::MeshEntity* >* ep
                                        , int geom, int step = 1);

// Given an entity and a set of vertices, check if at least two vertices bound
// the entity.
bool vd_chk_ent_v(apf::Mesh2* m, apf::MeshEntity* ent, 
                                std::vector<apf::MeshEntity* >* es_v);

// Given an element entity set ep1, print lower adjacencies.
void vd_print_down(apf::Mesh2* m, const std::vector<apf::MeshEntity* >* ep_in);

void vd_print_vert_tri_x(apf::Mesh2* m, apf::MeshEntity* vert);
//----------------------------------------------------

void copy_ent_set(Entity_set* es, apf::Up up);
void copy_ent_set(Entity_set* es, apf::Downward down, int size);
void copy_ent_set(Entity_set* es_copy, Entity_set* es_in);
void copy_ent_set(Entity_set* es_copy, std::vector<apf::MeshEntity*>* es_in);
void copy_ent_set2(Entity_set* es, apf::MeshEntity** e, int size);
void init_ent_set(Entity_set* es, apf::MeshEntity* e);
void init_ent_set(Entity_set* es, std::vector<apf::MeshEntity*>* e);


// Revert entity set:
void rev_ent_set(Entity_set* es);

void print_ent_set(apf::Mesh2* m, Entity_set* es);

// TODO memory problem initializing this at each time step in boundary motion 
// calculations.
class m_lookup {
  public:
    // Given the vertex, return the adjacent edges within the triangle.
    int lookup_tri_vd [3][2];
    int lookup_tri_ed [3][2];
    // Given the vertex, return the adjacent triangles.
    int lookup_tet_surf [4][3];
    // Given the vertex, return the edges across within the tetrahedron.
    int lookup_v_tet_e_x [4][3];
    // Given the vertex, return the edges across within the triangle.
    int lookup_v_tri_e_x [3];
    // Given the edge, return the vertex across within the triangle.
    int lookup_v_tri_x_e [3];
    // Given the vertex, and an edge return the other edge adjacent to the 
    // vertex within the triangle.
    int lookup_v_tri_n_e [3][3];
    // Given two vertex indices, return the index of the other one. Not  
    // necessary but cleaner.
    int lookup_triverts2 [3][3];
    // Given two triangle indices, find the edge not contained within one of the
    // triangles.
    int lookup_tri_edges [4][4];

    m_lookup () {
      lookup_tri_vd[0][0] = 1; lookup_tri_vd[0][1] = 2;
      lookup_tri_vd[1][0] = 0; lookup_tri_vd[1][1] = 2;
      lookup_tri_vd[2][0] = 1; lookup_tri_vd[2][1] = 0;

      lookup_tri_ed[0][0] = 0; lookup_tri_ed[0][1] = 2;
      lookup_tri_ed[1][0] = 0; lookup_tri_ed[1][1] = 1;
      lookup_tri_ed[2][0] = 1; lookup_tri_ed[2][1] = 2;

      lookup_tet_surf[0][0] = 0; lookup_tet_surf[0][1] = 1; lookup_tet_surf[0][2] = 3;
      lookup_tet_surf[1][0] = 0; lookup_tet_surf[1][1] = 1; lookup_tet_surf[1][2] = 2;
      lookup_tet_surf[2][0] = 0; lookup_tet_surf[2][1] = 2; lookup_tet_surf[2][2] = 3;
      lookup_tet_surf[3][0] = 1; lookup_tet_surf[3][1] = 2; lookup_tet_surf[3][2] = 3;

      lookup_v_tet_e_x[0][0] = 1; lookup_v_tet_e_x[0][1] = 4; lookup_v_tet_e_x[0][2] = 5;
      lookup_v_tet_e_x[1][0] = 2; lookup_v_tet_e_x[1][1] = 5; lookup_v_tet_e_x[1][2] = 3;
      lookup_v_tet_e_x[2][0] = 0; lookup_v_tet_e_x[2][1] = 4; lookup_v_tet_e_x[2][2] = 3;
      lookup_v_tet_e_x[3][0] = 0; lookup_v_tet_e_x[3][1] = 1; lookup_v_tet_e_x[3][2] = 2;

      lookup_v_tri_e_x[0] = 1; lookup_v_tri_e_x[1] = 2; lookup_v_tri_e_x[2] = 0;
      lookup_v_tri_x_e[0] = 2; lookup_v_tri_x_e[1] = 0; lookup_v_tri_x_e[2] = 1;

      lookup_v_tri_n_e[0][0] = 2; lookup_v_tri_n_e[0][1] = -1; lookup_v_tri_n_e[0][2] = 0;
      lookup_v_tri_n_e[1][0] = 1; lookup_v_tri_n_e[1][1] = 0; lookup_v_tri_n_e[1][2] = -1;
      lookup_v_tri_n_e[2][0] = -1; lookup_v_tri_n_e[2][1] = 2; lookup_v_tri_n_e[2][2] = 1;

      lookup_triverts2[0][0] = -1; lookup_triverts2[0][1] = 2; lookup_triverts2[0][2] = 1;
      lookup_triverts2[1][0] = 2; lookup_triverts2[1][1] = -1; lookup_triverts2[1][2] = 0;
      lookup_triverts2[2][0] = 1; lookup_triverts2[2][1] = 0; lookup_triverts2[2][2] = -1;

      lookup_tri_edges[0][0] = -1; lookup_tri_edges[0][1] = 5; lookup_tri_edges[0][2] = 3; lookup_tri_edges[0][3] = 4;
      lookup_tri_edges[1][0] = 5; lookup_tri_edges[1][1] = -1; lookup_tri_edges[1][2] = 2; lookup_tri_edges[1][3] = 1;
      lookup_tri_edges[2][0] = 3; lookup_tri_edges[2][1] = 2; lookup_tri_edges[2][2] = -1; lookup_tri_edges[2][3] = 0;
      lookup_tri_edges[3][0] = 4; lookup_tri_edges[3][1] = 1; lookup_tri_edges[3][2] = 0; lookup_tri_edges[3][3] = -1;

    }
};


struct geom_list
{
public:
  int n;
  int geom[1600];
  geom_list() : n(0) {}
};

// To be defined in topo_topo.h.
struct ent_conn;

// Given an element entity set find the number of unique geometric memberships,
// of the same dimension.
void find_elem_set_uniq(apf::Mesh2* m, Entity_set* es, geom_list* g_list);

// Create tags for element IDs and geometric memberships.
void vd_tag_mesh(apf::Mesh2* m);

// Remove numberings attached to a mesh.
void vd_rem_tag(apf::Mesh2* m);

// Remove the given numbering attached to a mesh if exists.
void vd_rem_tag(apf::Mesh2* m, const char* set_name);

void vd_print_numberings(apf::Mesh2* m);

// To check if deformation causes geometric dissassociation. Given a mesh, 
// tag vertices based on their geometric membership, 
// i.e 0 for quad, 1 for trip, 2 for surf, 3 for internal.
void vd_tag_vert(apf::Mesh2* m);

// Given an entity set, create set tags.
// findNumbering is problematic, so destroy and rewrite.
//void vd_tag_vert_set(apf::Mesh2* m, Entity_set* es, const char* set_name);

// Given an entity set, create set tags.
// findNumbering is problematic, so destroy and rewrite.
void vd_tag_set(apf::Mesh2* m, Entity_set* es, const char* set_name);



// Prints triple edges.
void vd_print_trip(apf::Mesh2* m);

// Find the same dimensional entities member of a certain cell and dimension.
// Store in ep*.
void vd_find_ent_geom(apf::Mesh2* m, Entity_set* ep, int cell_id, int dim, 
                                                            int dim_e = -1);

// Check if the edge is bounded by two internal vertices, one internal and one 
// external or two external vertices. Return the vertex index of the internal 
// vertex if one is external (return 2, 0 or 1, -1).
// Store in ep*.
int vd_chk_edge_int(apf::Mesh2* m, apf::MeshEntity* edge
                                  , apf::ModelEntity* cell);

// Print the entities with their cell memberships.
void vd_print_ent(apf::Mesh2* m);

// Find the edges member of a certain cell and dimension.
// Store in ep*.
void vd_find_edge_topo(apf::Mesh2* m, Entity_set* ep, int cell_id, int dim);

// Find the edges bounded by two vertices. If does not exist return false.
bool vd_find_edge(apf::Mesh2* m, apf::MeshEntity** vert, 
                                 apf::MeshEntity** edge);

bool vd_find_edge(apf::Mesh2* m, apf::MeshEntity* v1, apf::MeshEntity* v2, 
                                 apf::MeshEntity** edge);

// Given an edge, a triangle and a vertex, find the edge bounding triangle, 
// joining the edge in the vertex. Return NULL if not possible.
apf::MeshEntity* vd_find_esjv(apf::Mesh2* m, apf::MeshEntity* edge, apf::MeshEntity* tri, apf::MeshEntity* v);

// Given two edges, return the joint vertex. If no joint vertex exists, return
// null
apf::MeshEntity* vd_find_vert(apf::Mesh2* m, apf::MeshEntity* edge1,
                                             apf::MeshEntity* edge2);

// Given an entity set, find the entities member of a certain geometry.
// Store in ep*. If no geometry specified, return all, member of any geometry.
void vd_ext_ent_geom(apf::Mesh2* m, const Entity_set* ep_in, Entity_set* ep, int geom = 0);

// LOOK: This could be hardwired to vd_calc_disp by a continue statement.
// Check if given surface, member of a surface geometry is on domain geometry
// surface.
bool vd_chk_surf_dom2(apf::Mesh2* m, apf::MeshEntity* b_surf);

// Check if given surface, member of a surface geometry is on domain geometry
// surface. Works mesh level by checking the number of element adjacencies.
// A single adjacency signief boundary triangle.
bool vd_chk_surf_dom(apf::Mesh2* m, apf::MeshEntity* b_surf);
bool vd_chk_surf_dom(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es_surf);

// Check if given vertex is on domain boundary.
bool vd_chk_vert_dom(apf::Mesh2* m, apf::MeshEntity* vert);

// Check if given entity(vertex, edge, surface) is on domain boundary.
bool vd_chk_ent_dom(apf::Mesh2* m, apf::MeshEntity* e);

// Given a mesh object, extract triple edges.
void vd_ext_trip(apf::Mesh2* m, Entity_set* es_edge);

// Given a mesh object, extract boundary and domain surfaces.
void vd_ext_bsurf(apf::Mesh2* m, Entity_set* es_surf, Entity_set* es_dom);

// Given a mesh object, and domain surface list extract quadruple points, which
// lie on the interior, not on the domain boundary. For now, assume quadruple
// points indeed are members of at least four different grains.
void vd_ext_quad(apf::Mesh2* m, Entity_set* es_quad);

// Given a vertex MeshEntity, find the surface couples that bound geometries 
// around that vertex, joined by edges adjacent to that vertex.
void vd_ext_surf_coup(apf::Mesh2* m, apf::MeshEntity* vert);

// Print the surface geometric adjacencies of given dimension adj_dim.
void vd_surf_adj(gmi_model* m, int ent_dim, int adj_dim);

// Print the tag of the model attached to the mesh entity.
int vd_mesh_mdl_adj(apf::Mesh2* m, apf::MeshEntity* e);

// Print the model entities attached to mesh vertices.
void vd_ext_do_mesh(apf::Mesh2* m);

// Print the adjacencies of an entity, of a distance n.
// This simply lists upper and lower adjacencies of an entity, and calls their
// adjacencies of distance n-1.
void vd_ent_deg_adj(apf::Mesh2* m, apf::MeshEntity* ent, int n);

// Given a mesh, going over the entities, print the lower and upper adjacencies
// and their type, whenever applicable.
void vd_mesh_up_adj(apf::Mesh2* m);

// Given a mesh, going over the entities, print the lower and upper adjacencies
// if they are not valid.
void vd_mesh_bad_adj(apf::Mesh2* m);

/////////////////////////////////////////////////////////////
// Finding lower and upper dimensional adjacencies:
/////////////////////////////////////////////////////////////

// Given a mesh, and an entity, print the lower and upper adjacencies
// if they are not valid. Return true if not valid.
bool vd_mesh_bad_ent(apf::Mesh2* m, apf::MeshEntity* e);

// Given a mesh entity, print the number of upper adjacencies and their type.
void vd_count_ent_up(apf::Mesh2* m, apf::MeshEntity* e);

// Given two entity sets, remove the elements of the second list from the 
// first or write the common elements into the first.
// swtch = 0, intersect; swtch = 1, remove.
// LOOK intersect would not work!
void vd_remove_set(Entity_set* ep1, Entity_set* ep2);
void vd_inter_set(Entity_set* ep1, Entity_set* ep2);

// Given two up adjacency entity arrays, obtain a list of unique entities.
void vd_merge_set(Entity_set* ep1, Entity_set* ep2);

// Given an entity set ep1, find one level higher adjacencies if possible.
// Store in ep*.
void vd_set_up(apf::Mesh2* m, const Entity_set* ep1, Entity_set* ep);

// Given an entity set ep_in, find one level lower adjacencies if possible.
// Store in ep*.
void vd_set_down(apf::Mesh2* m, const Entity_set* ep_in, Entity_set* ep, int level = 1);

//void vd_set_down(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in, std::vector<apf::MeshEntity*>* ep, int level = 1);

// For two d-dim entities, starting from d-1 dimensional adjacencies, 
// find the highest dimensional common lower adjacency. If none could be found
// return NULL.
apf::MeshEntity* vd_find_common_adj(apf::Mesh2* m, apf::MeshEntity* e1, 
                                             apf::MeshEntity* e2);

// For two edges, find the joining edge. assert one exists.
apf::MeshEntity* vd_find_tris_edge(apf::Mesh2* m, apf::MeshEntity* tri1, 
                                             apf::MeshEntity* tri2);


std::vector<apf::MeshEntity*> vd_get_3c_tet_bound(apf::Mesh2* m, std::vector<apf::MeshEntity*> & tet_in);

// Given an edge entity set ep_in, extract those, which rest on the grain 
// surface.
// (i.e. a surface entity set, whose adjacent elements belong to geometry n.)
// Store in ep*.
// LOOK write functions that do something simple, and only that!
void vd_chk_edge_grain(apf::Mesh2* m, const Entity_set* ep_in, Entity_set* ep
                                        , int geom);

// Given an surface entity set ep_in, extract those, which rest on the grain 
// surface.
// Store in ep*.
// LOOK write functions that do something simple, and only that!
void vd_chk_surf_grain(apf::Mesh2* m, const Entity_set* ep_in, Entity_set* ep
                                        , int geom);

// Given an entity set ep_in, extract those, one of d dim upper adjacency of 
// which belongs to a certain geometry of that adjacencies own dimension.
// (i.e. a surface entity set, whose adjacent elements belong to geometry n.)
// Store in ep*.
// LOOK write functions that do something simple, and only that!
// So write a function, which extracts edges on a grain surface. 
void vd_chk_set_geom(apf::Mesh2* m, const Entity_set* ep_in, Entity_set* ep
                                        , int geom, int step = 1);

// Given an element entity set ep1, print lower adjacencies.
void vd_print_down(apf::Mesh2* m, const Entity_set* ep_in);

void vd_print_down(apf::Mesh2* m, apf::MeshEntity* e_in);

std::pair<int,int> vd_cell_pair(apf::Mesh2* m, apf::MeshEntity* ent);

// Given a mesh, go over vertices, print number of entities in given upper 
// adjacency sets. 
void print_vert_adj(apf::Mesh2* m);

// Given a mesh, go over vertices, print their positions. 
void vd_print_pos_vert(apf::Mesh2* m);

// Given a set of positions, print them in MATLAB compatible matrix format. 
void vd_print_pos(std::vector<apf::Vector3> &v);
// Given a position, print them in MATLAB compatible matrix format. 
void vd_print_pos(apf::Vector3 &v);

// Tag the entities around a vertex and save the mesh.
void vd_save_vtk_vert(apf::Mesh2* m, std::vector<apf::MeshEntity*> * vert, 
                  char const* FileName ="./tempmesh/vert_set");

// Print the tets adjacent to the vertices bounding the given entity.
void vd_save_vtk_ent(apf::Mesh2* m, apf::MeshEntity* ent, char const* FileName ="./tempmesh/set");

// Print the tets adjacent to the vertices bounding the given set of entities.
void vd_save_vtk_ent(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ent, char const* FileName ="./tempmesh/set");

// Tag the entities around a vertex and save the mesh.
void vd_save_vtk_vert(apf::Mesh2* m, apf::MeshEntity* vert, 
                  char const* FileName ="./tempmesh/beforecol");

// Copy a set of entities from a mesh to a new mesh, and save into a VTK file.
void vd_save_vtk_set(apf::Mesh2* m, std::vector<Entity_set>* e_set, 
                  char const* FileName ="./tempmesh/set");

void vd_save_vtk_set(apf::Mesh2* m, std::vector< std::vector<apf::MeshEntity* > >* e_set, char const* FileName ="./tempmesh/set");

void safe_mkdir(const char* path);

void save_vtk_mov(apf::Mesh2* m, int t_in, const char* FileName ="./outputtemp/series");

// TODO remove. GDB doesn't print content of map containers. This is a quick 
// solution to print the map value associated with key.
double print_map_i2d(std::map<int, double>& map, int key);
#endif


