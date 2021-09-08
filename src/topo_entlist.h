#ifndef TOPO_ENTLIST_H
#define TOPO_ENTLIST_H

#include <vector>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <PCU.h>

#include "topo_extinfo.h"
#include "topo_topo.h"


// This contains the cell forward list of entities associated with the mesh.
// Main operations are:
// List the mesh entities by dimension and cell membership.
// Refresh the list.
// Find entities by cell.
// Find the number of entities in each cell.

// TODO In the future, the counting should be done in parallel.

// This does two loops over all entities. One for counting and one for registering.
class vd_entlist{
  private:
    apf::Mesh2* m;
    cell_base* c_base;

    // This list contains the cells of 
    //std::vector<std::vector<apf::ModelEntity* > > cell_lookup;
    // Temporary list to contain number of entities of each dimension for each 
    // cell.
    std::vector<std::vector<int > > ent_count;

    // Copy constructor and assignment private:
    vd_entlist(const vd_entlist &obj);
    vd_entlist& operator=( const vd_entlist& obj );

  public:
    // m dimensional cells containing n = 0, 1, 2, and 3-dimensional entities, 
    // with 3 >= m >= n. Cells are ordered in ascending dimension and index and 
    // entities are grouped in same cell lists of same dimensional entity lists.
    // m-dim cell:
    //   0-dim list,
    //   1-dim list,
    //   2-dim list,
    //   3-dim list.
    std::vector<std::vector<std::vector<std::vector<apf::MeshEntity* > > > > e;

    vd_entlist(apf::Mesh2* m_in = NULL, cell_base* c_in = NULL);

    void change_mesh(apf::Mesh2* m_in, cell_base* c_in);
    void refresh();
    // Given the entity dimension, and cell dimension and id, get the entities.
    void get_entities_gmi(int ent_dim, int cell_dim, int cell_id, 
          std::vector<apf::MeshEntity* >* es_ent);

    void clear();
    void clear_count();
    void clear_ent();

    ~vd_entlist();
};

// Used for registering the mesh entities around a vertex.
class vd_entlist_v{
  private:
    apf::Mesh2* m;
    std::vector<apf::MeshEntity*> v_ctr;
    cell_base* c_base;

    // This list contains the cells of 
    //std::vector<std::vector<apf::ModelEntity* > > cell_lookup;
    // Temporary list to contain number of entities of each dimension for each 
    // cell.
    std::vector<std::vector<int > > ent_count;
    void get_vert();

    // Copy constructor and assignment private:
    vd_entlist_v(const vd_entlist_v &obj);
    vd_entlist_v& operator=( const vd_entlist_v& obj );

  public:
    // m dimensional cells containing n = 0, 1, 2, and 3-dimensional entities, 
    // with 3 >= m >= n. Cells are ordered in ascending dimension and index and 
    // entities are grouped in same cell lists of same dimensional entity lists.
    // m-dim cell:
    //   0-dim list,
    //   1-dim list,
    //   2-dim list,
    //   3-dim list.
    // e.at(m).at(cell_id).at(n).at(ent_id)
    std::vector<std::vector<std::vector<std::vector<apf::MeshEntity* > > > > e;
    std::vector<std::vector<apf::MeshEntity* > > es;

    vd_entlist_v();
    vd_entlist_v(apf::Mesh2* m_in, apf::MeshEntity* v_in, cell_base* c_in);
    vd_entlist_v(apf::Mesh2* m_in, 
              std::vector<apf::MeshEntity*>* v_in, cell_base* c_in);

    void change_mesh(apf::Mesh2* m_in, apf::MeshEntity* v_in, cell_base* c_in);
    void change_mesh(apf::Mesh2* m_in, 
              std::vector<apf::MeshEntity*>* v_in, cell_base* c_in);
    // Given the entity dimension, and cell dimension and id, get the entities.
    void get_entities_gmi(int ent_dim, int cell_dim, int cell_id, 
          std::vector<apf::MeshEntity* >* es_ent);

    void print();
    void print_v_pos();
    void clear();
    void clear_count();
    void clear_ent();

    ~vd_entlist_v();
};

#endif
