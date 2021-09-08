#include <map>
#include <vector>
#include <iterator>
#include <algorithm>

// Obtained from apfVTK for safe_mkdir
// === includes for safe_mkdir ===
#include <reel.h>
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */
// 

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
#include "topo_rand.h"

#include "topo_geom.h"

// Given the vertex, return the adjacent edges within the triangle.
int vd_lookup_tri_ed [3][2] = {{0,2},{0,1},{1,2}};

void dummy_clear_stop() {
}

void dummy_chk_stop() {
}

// 
int findIn(const std::vector<apf::MeshEntity* >* es, int n, apf::MeshEntity* e) {
  std::vector<apf::MeshEntity*>::const_iterator it;
  it = find(es->begin(), es->end(), e);
  if(it == es->end())
    return -1;
  return std::distance(es->begin(), it);
}
void copy_ent_set(std::vector<apf::MeshEntity*>* es, apf::Up up) {
  es->resize(0);
  es->resize(up.n);
  for (int i = 0; i < up.n; i ++) {
    es->at(i) = up.e[i];
  }
}

void copy_ent_set(std::vector<apf::MeshEntity*>* es, 
                                    apf::Downward down, int size) {
  es->resize(size);
  for (int i = 0; i < size; i++) {
    es->at(i) = down[i];
  }
}

void init_ent_set(std::vector<apf::MeshEntity*>* es, apf::MeshEntity* e) {
  es->resize(1);
  es->at(0) = e;
}

void rev_ent_set(std::vector<apf::MeshEntity*>* es) {
  std::reverse(es->begin(), es->end());
}

bool vd_vert_set(apf::Mesh2* m, std::vector<apf::MeshEntity* >* es_in, std::vector<apf::MeshEntity* >* es_out) {

  if(es_in->size() > 0) {
    int e_type = m->getType(es_in->at(0));
    int d = m->typeDimension[e_type];
    if(d == 0)
      return false;

    //double t0 = MPI_Wtime();
    vd_set_down(m, es_in, es_out, d);
    //double t1 = MPI_Wtime();
    //std::cout << "vd_set_down " << t1 - t0 << std::endl;

    //for(int i = 0; i < es_out->size(); i++) {
    //  std::cout << es_out->at(i) << ", ";
    //}
    //std::cout << std::endl;

    //t0 = MPI_Wtime();
    //vd_set_down2(m, es_in, es_out, d);
    //t1 = MPI_Wtime();
    //std::cout << "vd_set_down2 " << t1 - t0 << std::endl;
    //for(int i = 0; i < es_out->size(); i++) {
    //  std::cout << es_out->at(i) << ", ";
    //}
    //std::cout << std::endl;

    return true;
  }
  else {
    return false;
  }
}

void print_ent_set(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es) {
  for (int i = 0; i < es->size(); i++) {
    int e_type = m->getType(es->at(i));
    int d = m->typeDimension[e_type];
    int type = m->getModelType(m->toModel(es->at(i)));
    int tag = m->getModelTag(m->toModel(es->at(i)));
    std::cout << d << "-ent " << es->at(i)
              << " " << type << "c" << tag << std::endl;
  }

}

// Given an element entity set find the number of unique geometric memberships,
// of the same dimension.
void find_elem_set_uniq(apf::Mesh2* m, 
                    std::vector<apf::MeshEntity*>* es,
                    std::vector<int>* g_list) {
  // e_type is from apf::Mesh:Type, 4 is tetrahedron.
  int e_type = m->getType(es->at(0));
  int d = m->typeDimension[e_type];

  std::map<int, bool> g_chk{};
  int count = 0;

  for (int i = 0; i < es->size(); i ++) {
    apf::ModelEntity* em = m->toModel(es->at(i));
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    // Check if geometric dimension is the same as the mesh entity dimension.
    if (em_type == d) {
      // Check if it already exists.
      if(!g_chk[em_tag]) {
        g_chk[em_tag] = true;
        count = count + 1;
      }
    }
  }

  g_list->resize(count);
  for (int i = 0; i < g_list->size(); i ++) {
    g_list->at(i) = 0;
  }

  g_chk.clear();
  count = 0;

  for (int i = 0; i < es->size(); i ++) {
    apf::ModelEntity* em = m->toModel(es->at(i));
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    // Check if geometric dimension is the same as the mesh entity dimension.
    if (em_type == d) {
      if(!g_chk[em_tag]) {
        g_chk[em_tag] = true;
        g_list->at(count) = em_tag;
        count = count + 1;
      }
    }
  }

}

void vd_tag_set(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es, const char* set_name) {

  if(es->size() == 0)
    return;

  int ent_type = m->getType(es->at(0));
  int d = m->typeDimension[ent_type];

  apf::Numbering* setnumbering;
  if (m->findNumbering(set_name)) {
    std::cout<< set_name << " tag already exists" << std::endl;
    setnumbering = m->findNumbering(set_name);
    apf::destroyNumbering(setnumbering);
    setnumbering = apf::createNumbering(m, set_name, apf::getConstant(d),1);
  }
  else {
    setnumbering = apf::createNumbering(m, set_name, apf::getConstant(d),1);
  }

  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(d);

  std::cout<< "Writing tag " << set_name << "...";
  while ((e = m->iterate(it))) 
  {
    if ((es->size() == 0) || (findIn(es, es->size(), e) == -1)) {
      apf::number(setnumbering, e, 0, 0, 0);
    }
    else {
      apf::number(setnumbering, e, 0, 0, 1);
    }
  }
  m->end(it);

  std::cout<< " Written." << std::endl;

}

// Find the dim_e dimensional entities member of a certain cell and dimension.
// Store in ep*. Default value of dim_e is the same as dim.
void vd_find_ent_geom(apf::Mesh2* m, 
                      std::vector<apf::MeshEntity*>* ep, int cell_id, int dim, 
                                            int dim_e) {

  ep->resize(0);

  if (dim_e == -1)
    dim_e = dim;
  //printf("Checking for %dcell%d", dim, cell_id);

  int cell_size = 0;
  apf::ModelEntity* cell = m->findModelEntity(dim, cell_id);

  apf::MeshEntity* e;

  int count = 0;

  apf::MeshIterator* it_e = m->begin(dim_e);
  while (e = m->iterate(it_e)) {

    if (cell == m->toModel(e)) {
      //printf("Vertex %p, cell %p, %dcell%d\n", e, m->toModel(e), 
      //        m->getModelType(m->toModel(e)), m->getModelTag(m->toModel(e)));
      count = count + 1;
    }
  }
  m->end(it_e);

  ep->reserve(count);

  it_e = m->begin(dim_e);
  while (e = m->iterate(it_e)) {

    if (cell == m->toModel(e)) {
      //printf("Vertex %p, cell %p, %dcell%d\n", e, m->toModel(e), 
      //        m->getModelType(m->toModel(e)), m->getModelTag(m->toModel(e)));
      ep->push_back(e);
    }
  }
  m->end(it_e);

}
/*
// Find the same dimensional entities member of a certain cell and dimension.
// Store in ep*.
void vd_find_ent_geom(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ep, int cell_id, int dim) {
  //printf("Checking for %dcell%d", dim, cell_id);
  ep->resize(0);

  int cell_size = 0;
  apf::ModelEntity* cell = m->findModelEntity(dim, cell_id);


  apf::MeshEntity* e;
  int count = 0;

  apf::MeshIterator* it_e = m->begin(dim);
  while (e = m->iterate(it_e)) {

    if (cell == m->toModel(e)) {
      //printf("Vertex %p, cell %p, %dcell%d\n", e, m->toModel(e), 
      //        m->getModelType(m->toModel(e)), m->getModelTag(m->toModel(e)));
      count = count + 1;
    }
  }
  m->end(it_e);

  ep->reserve(count);

  it_e = m->begin(dim);
  while (e = m->iterate(it_e)) {

    if (cell == m->toModel(e)) {
      //printf("Vertex %p, cell %p, %dcell%d\n", e, m->toModel(e), 
      //        m->getModelType(m->toModel(e)), m->getModelTag(m->toModel(e)));
      ep->push_back(e);
    }
  }
  m->end(it_e);
}
*/
// Find the edges member of a certain cell and dimension.
// Store in ep*.
void vd_find_edge_topo(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ep, 
                       int cell_id, int dim) {

  ep->resize(0);
  apf::ModelEntity* cell = m->findModelEntity(dim, cell_id);
  apf::MeshEntity* e;

  int count = 0;

  apf::MeshIterator* it_e = m->begin(1);
  while (e = m->iterate(it_e)) {

    if (cell == m->toModel(e)) {
      //int ent_type = m->getType(e);
      //int d = m->typeDimension[ent_type];
      //printf("%d-Ent %p, cell %p, %dcell%d\n", d, e, m->toModel(e), 
      //        m->getModelType(m->toModel(e)), m->getModelTag(m->toModel(e)));
      count = count + 1;
    }
  }
  m->end(it_e);

  ep->reserve(count);

  it_e = m->begin(1);
  while (e = m->iterate(it_e)) {

    if (cell == m->toModel(e)) {
      //printf("Vertex %p, cell %p, %dcell%d\n", e, m->toModel(e), 
      //        m->getModelType(m->toModel(e)), m->getModelTag(m->toModel(e)));
      ep->push_back(e);
    }
  }
  m->end(it_e);

}

void vd_find_ent_topo(apf::Mesh2* m, apf::MeshEntity* v, 
                std::vector<apf::MeshEntity* >* ep, int cell_id, int dim) {

  assert(dim < 4 and dim > 0);
  ep->resize(0);

  std::vector<std::vector<apf::MeshEntity* > > ent
          (0, std::vector<apf::MeshEntity* >(0));

  ent.resize(dim);
  vd_set_up(m, v, &ent.at(0));

  for(int i = 0; i < dim-1; i++)
    vd_set_up(m, &ent.at(i), &ent.at(i+1));

  ep->reserve(ent.at(dim-1).size());

  for(int i = 0; i < ent.at(dim-1).size(); i++) {
    apf::ModelEntity* em = m->toModel(ent.at(dim-1).at(i));
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);

    if(em_type == dim and em_tag == cell_id) {
      ep->push_back(ent.at(dim-1).at(i));
    }
  }
  dummy_clear_stop();

  for(int i = 0; i < ent.size(); i++)
    ent.at(i).clear();

}

// Given an entity set, find the entities member of a certain dimension and 
// topology.
// Store in ep*. If no geometry specified, return all, member of any topology 
// of the given dimension.
void vd_ext_ent_geom(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in, 
                        std::vector<apf::MeshEntity*>* ep, int geom) {

  ep->resize(0);
  ep->reserve(ep_in->size());

  int ent_type = m->getType(ep_in->at(0));
  int d = m->typeDimension[ent_type];

  if ((geom == 0)) {
    for (int i = 0; i < ep_in->size(); i++) {
      apf::ModelEntity* em = m->toModel(ep_in->at(i));
      int em_type = m->getModelType(em);
      int em_tag = m->getModelTag(em);

      if (d == em_type)
        ep->push_back(ep_in->at(i));
    }
  }
  else {

    for (int i = 0; i < ep_in->size(); i++) {
      apf::ModelEntity* em = m->toModel(ep_in->at(i));
      int em_type = m->getModelType(em);
      int em_tag = m->getModelTag(em);

      if ((geom == em_tag) && (d == em_type)) {
        ep->push_back(ep_in->at(i));
      }
    }
  }
}

// Given a mesh object, extract triple edges.
void vd_ext_trip(apf::Mesh2* m, std::vector<apf::MeshEntity* >* es_edge) {

  int dim = 1;
  apf::MeshEntity* edge;
  int trip_index = 0;
  es_edge->resize(0);

  apf::MeshIterator* it_e = m->begin(dim);
  while ((edge = m->iterate(it_e))) {
    apf::ModelEntity* em = m->toModel(edge);
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 1) {
      trip_index++;
    }
  }
  m->end(it_e);

  es_edge->reserve(trip_index);

  it_e = m->begin(dim);
  while ((edge = m->iterate(it_e))) {
    apf::ModelEntity* em = m->toModel(edge);
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 1) {
      es_edge->push_back(edge);
    }
  }
  m->end(it_e);

  printf("There are %d triple edges.\n", trip_index);
}

// Given a mesh object, extract boundary and domain surfaces.
void vd_ext_bsurf(apf::Mesh2* m, std::vector<apf::MeshEntity* >* es_surf, 
                                 std::vector<apf::MeshEntity* >* es_dom) {

  es_surf->resize(0);
  es_dom->resize(0);

  int dim = 2;

  apf::Up up;
  apf::MeshEntity* surf;

  int surf_index = 0;
  int dom_index = 0;

  apf::MeshIterator* it_s = m->begin(dim);
  while ((surf = m->iterate(it_s))) {
    apf::ModelEntity* em = m->toModel(surf);
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 2) {
      m->getUp(surf, up);
      if (up.n == 1) {
        dom_index ++;
      }
      else if (up.n == 2) {
        surf_index ++;
      }
      else {
        printf("Surface %p has inappropariate upper adjacency.\n",(void*)surf);
      }
    }
  }
  m->end(it_s);

  es_surf->reserve(surf_index);
  es_dom->reserve(dom_index);

  it_s = m->begin(dim);
  while ((surf = m->iterate(it_s))) {
    apf::ModelEntity* em = m->toModel(surf);
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 2) {
      m->getUp(surf, up);
      if (up.n == 1) {
        es_dom->push_back(surf);
      }
      else if (up.n == 2) {
        es_surf->push_back(surf);
      }
    }
  }
  m->end(it_s);
  printf("There are %d domain boundary triangles.\n",(int)es_dom->size());
  printf("There are %d grain boundary triangles.\n",(int)es_surf->size());
}

// LOOK: This is not necessary to find quadruple points. Currently a PLACEHOLDER
// Given a mesh object, and domain surface list extract quadruple points, which
// lie on the interior, not on the domain boundary. For now, assume quadruple
// points indeed are members of at least four different grains.
void vd_ext_quad(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es_quad) {
}


// Given two entity sets, remove the elements of the second list from the 
// first or write the common elements into the first.
// swtch = 0, intersect; swtch = 1, remove.
// TODO LOOK somehow, the first entry of the set of the removed entities is not 
// matched. 
void vd_remove_set(std::vector<apf::MeshEntity*>* ep1, 
                   const std::vector<apf::MeshEntity*>* ep2) {

  int rem_size = 0;

  for (int i = 0; i < ep1->size(); i++) {
    int i1 = findIn(ep2, ep2->size(), ep1->at(i));

    if ( i1 >= 0) {
      rem_size++;
    }
    else {
      ep1->at(i-rem_size) = ep1->at(i);
    }
  }
  ep1->resize(ep1->size() - rem_size);
}

void vd_inter_set(std::vector<apf::MeshEntity*>* ep1, 
                   std::vector<apf::MeshEntity*>* ep2) {

  int rem_size = 0;

  for (int i = 0; i < ep1->size(); i++) {
    int i1 = findIn(ep2, ep2->size(), ep1->at(i));
    if ( i1 < 0) {
      rem_size++;
    }
    else {
      ep1->at(i-rem_size) = ep1->at(i);
    }
  }

  printf("Set size:%d, remove size:%d, ", (int)ep1->size(), rem_size);
  ep1->resize(ep1->size() - rem_size);
  printf("New size%d.\n", (int)ep1->size());
}

// Given two up adjacency entity arrays, obtain a list of unique entities.
void vd_merge_set(std::vector<apf::MeshEntity*>* ep1, 
                  std::vector<apf::MeshEntity*>* ep2) {

  std::vector<apf::MeshEntity*> ep_list(0);
  ep_list.reserve(ep2->size());

  int i2 = 0;
  for (int i = 0; i < ep2->size(); i ++) {
    // printf("ep1->n is %d.\n", ep1->n);
    i2 = findIn(ep1, ep1->size(), ep2->at(i));
    // printf("i2 is %d Pointers: %p %p.\n", i2, ep2.e[i], ep1->e[i2]);
    if ( i2 == -1) {
      ep_list.push_back(ep2->at(i));
    }
  }

  ep1->reserve(ep1->size()+ep_list.size());
  for (int i = 0; i < ep_list.size(); i ++) {
    ep1->push_back(ep_list.at(i));
  }

}

// Given an entity set ep1, find one level higher adjacencies if possible.
// Store in ep*.
void vd_set_up2(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in, 
                              std::vector<apf::MeshEntity*>* ep) {

  apf::Up up;
  std::vector<apf::MeshEntity*> es(0);

  ep->resize(0);

  int up_tot = 0;

  for (int i = 0; i < ep_in->size(); i++) {
    m->getUp(ep_in->at(i), up);

    if (up.n > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      up_tot = up_tot + up.n;
    }
  }
  ep->reserve(up_tot);
  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->size(); i++) {
    m->getUp(ep_in->at(i), up);

    if (up.n > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      copy_ent_set(&es, up);
      vd_merge_set(ep, &es);
    }
  }
}

// Given an entity set ep1, find one level higher adjacencies if possible.
// Store in ep*.
void vd_set_up(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in, 
                              std::vector<apf::MeshEntity*>* ep) {

  apf::Up up;
  std::vector<apf::MeshEntity*> es(0);

  ep->resize(0);

  int up_tot = 0;

  for (int i = 0; i < ep_in->size(); i++) {
    m->getUp(ep_in->at(i), up);

    if (up.n > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      up_tot = up_tot + up.n;
    }
  }
  ep->reserve(up_tot);
  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->size(); i++) {
    m->getUp(ep_in->at(i), up);

    if (up.n > 0) {
      for (int i = 0; i < up.n; i++) {
        ep->push_back(up.e[i]);
      }
    }
  }

  std::sort(ep->begin(), ep->end());
  std::vector<apf::MeshEntity*>::iterator it;
  it = std::unique(ep->begin(), ep->end());
  ep->resize(std::distance(ep->begin(),it));
}

// Given an entity set ep1, find one level higher adjacencies if possible.
// Store in ep*.
void vd_set_up(apf::Mesh2* m, apf::MeshEntity* e, 
                              std::vector<apf::MeshEntity*>* ep) {

  apf::Up up;

  ep->resize(0);

  m->getUp(e, up);
  copy_ent_set(ep, up);
}

// Get the tetrahedra adjacent to the vertices bounding a given set of entities.
void vd_get_tet(apf::Mesh2* m, const std::vector<apf::MeshEntity*> & ep1, 
                      std::vector<apf::MeshEntity*> & ep) {
  if(ep1.size() > 0) {
    std::vector<apf::MeshEntity*> vert(0);
    std::vector<std::vector<apf::MeshEntity*> > ents(2, 
                                      std::vector<apf::MeshEntity*> (0));

    int e_type = m->getType(ep1.at(0));
    int d = m->typeDimension[e_type];

    vd_set_down(m, &ep1, &vert, d);
    vd_set_up(m, &vert, &ents.at(0));
    vd_set_up(m, &ents.at(0), &ents.at(1));
    vd_set_up(m, &ents.at(1), &ep);
  }
}

// Given an entity set ep1, find one level higher adjacencies if possible.
// Store in ep*.
void vd_set_down(apf::Mesh2* m, apf::MeshEntity* e, 
                              std::vector<apf::MeshEntity*>* ep, int level) {

  // Get the dimension of the entity type in the set:
  int ent_type = m->getType(e);
  int d = m->typeDimension[ent_type];

  apf::Downward down;

  ep->resize(0);

  int e_sz = m->getDownward(e, d-level, down);
  copy_ent_set(ep, down, e_sz);
}

// Given an entity set ep1, find one level lower adjacencies if possible.
// Store in ep*.
void vd_set_down2(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in,
                                std::vector<apf::MeshEntity*>* ep, int level) {

  ep->resize(0);
  if (ep_in->size() < 1)
    return;

  // Get the dimension of the entity type in the set:
  int ent_type = m->getType(ep_in->at(0));
  int d = m->typeDimension[ent_type];

  if (d < 1)
    return;

  apf::Downward down;
  m->getDownward(ep_in->at(0), d-level, down);
  int ent_type2 = m->getType(down[0]);

  int size = apf::Mesh::adjacentCount[ent_type][ent_type2];

  ep->reserve(ep_in->size()*size);

  // printf("Type is %d, dimension is %d\n",d, ent_type);
  std::vector<apf::MeshEntity*> es(0);

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->size(); i++) {
    int size = m->getDownward(ep_in->at(i), d-level, down);

    if (size > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      copy_ent_set(&es, down, size);
      vd_merge_set(ep, &es);
    }
  }
}

// Given an entity set ep1, find one level lower adjacencies if possible.
// Store in ep*.
void vd_set_down(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in,
                                std::vector<apf::MeshEntity*>* ep, int level) {

  ep->resize(0);
  if (ep_in->size() < 1)
    return;

  // Get the dimension of the entity type in the set:
  int ent_type = m->getType(ep_in->at(0));
  int d = m->typeDimension[ent_type];

  if (d < 1)
    return;

  apf::Downward down;
  m->getDownward(ep_in->at(0), d-level, down);
  int ent_type2 = m->getType(down[0]);

  int size = apf::Mesh::adjacentCount[ent_type][ent_type2];
  ep->resize(ep_in->size()*size);

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->size(); i++) {
    m->getDownward(ep_in->at(i), d-level, down);
    for (int j = 0; j < size; j++) {
      ep->at(i*size+j) = down[j];
    }
  }
  std::sort(ep->begin(), ep->end());
  std::vector<apf::MeshEntity*>::iterator it;
  it = std::unique(ep->begin(), ep->end());
  ep->resize(std::distance(ep->begin(),it));
}

// Given an entity set, the entity adjacencies of given dimension.
void vd_set_adj(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in,
                            std::vector<apf::MeshEntity*>* ep_adj, int dim) {

  ep_adj->resize(0);
  if (ep_in->size() < 1)
    return;

  // Get the dimension of the entity type in the set:
  int ent_type = m->getType(ep_in->at(0));
  int d = m->typeDimension[ent_type];

  if(d == dim) {
    std::vector<apf::MeshEntity*> temp(0);
    vd_set_down(m, ep_in, &temp);
    vd_set_up(m, &temp, ep_adj);
    vd_remove_set(ep_adj, ep_in);
  }
  else {
    std::vector<std::vector<apf::MeshEntity*> > temp(0, 
                                  std::vector<apf::MeshEntity*>(0));
    temp.resize(std::abs(d-dim));
    if(d < dim) {
      vd_set_up(m, ep_in, &temp.at(0));
      for(int i = 0; i < dim-d-1; i++) {
        vd_set_up(m, &temp.at(i), &temp.at(i+1));
      }
      *ep_adj = temp.at(dim-d-1);
    }
    else if(d > dim) {
      vd_set_down(m, ep_in, &temp.at(d-dim-1));
      for(int i = d-dim-1; i > 0; i--) {
        vd_set_down(m, &temp.at(i), &temp.at(i-1));
      }
      *ep_adj = temp.at(0);
    }
  }

}

// Given an edge entity set ep_in, extract those, which rest on the grain 
// surface.
// (i.e. a surface entity set, whose adjacent elements belong to geometry n.)
// Store in ep*.
// OUTDATED Should be using topological information. Extinfo functions don't, 
// so perhaps transfer to another(new) header. 
// LOOK write functions that do something simple, and only that!
void vd_chk_edge_grain(apf::Mesh2* m, 
                        const std::vector<apf::MeshEntity*>* ep_in, 
                        std::vector<apf::MeshEntity*>* ep, int grain) {
  apf::Up up;
  ep->resize(0);
  ep->reserve(ep_in->size());

  assert(m->getType(ep_in->at(0)) == 1);
  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->size(); i++) {
    apf::ModelEntity* em = m->toModel(ep_in->at(i)); 
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    // If triple edge or surface edge:
    if ((em_type == 1) || (em_type == 2)) {
      m->getUp(ep_in->at(i), up);
      std::vector<apf::MeshEntity*> es_surf(0);
      std::vector<apf::MeshEntity*> es_elem(0);
      copy_ent_set(&es_surf, up);
      vd_set_up(m, &es_surf, &es_elem);

      int j = 0;
      while (j < es_elem.size()) {
        em = m->toModel(es_elem.at(j)); 
        em_type = m->getModelType(em);
        em_tag = m->getModelTag(em);
        if ( em_tag == grain ) {
          ep->push_back(ep_in->at(i));
          j = es_elem.size();
        }
        j++;
      }
    }
  }

}

// Given an surface entity set ep_in, extract those, which rest on the grain 
// surface.
// Store in ep*.
// LOOK write functions that do something simple, and only that!
void vd_chk_surf_grain(apf::Mesh2* m, 
             const std::vector<apf::MeshEntity*>* ep_in,
             std::vector<apf::MeshEntity*>* ep, int geom) {
  apf::Up up;
  ep->resize(0);
  ep->reserve(ep_in->size());

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->size(); i++) {
    apf::ModelEntity* em = m->toModel(ep_in->at(i)); 
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 2) {

      m->getUp(ep_in->at(i), up);

      int j = 0;
      for(int j = 0; j < up.n; j++) {
        em = m->toModel(up.e[j]); 
        em_type = m->getModelType(em);
        em_tag = m->getModelTag(em);
        if ( em_tag == geom ) {
          ep->push_back(ep_in->at(i));
          j = up.n;
        }
      }
    }
  }
}

// Given an entity set ep_in, extract those, one of d dim upper adjacency of 
// which belongs to a certain geometry of that adjacencies own dimension.
// (i.e. a surface entity set, whose adjacent elements belong to geometry n.)
// Store in ep*.
// LOOK write functions that do something simple, and only that!
// So write a function, which extracts edges on a grain surface. 
void vd_chk_set_geom(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in, 
                    std::vector<apf::MeshEntity*>* ep, int geom, int step) {

  // Get the dimension of the entity set elements:
  int ent_type = m->getType(ep_in->at(0));
  int d = m->typeDimension[ent_type];

  apf::Up up;

  ep->resize(0);
  ep->reserve(ep_in->size());

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->size(); i++) {
    m->getUp(ep_in->at(i), up);
    std::vector<apf::MeshEntity*> es_up(0);
    if (step == 1) {
      copy_ent_set(&es_up, up);
    }
    else if (step == 2) {
      std::vector<apf::MeshEntity*> es_1(0);
      copy_ent_set(&es_1, up);
      vd_set_up(m, &es_1, &es_up);
    }
    int j = 0;
    while (j < es_up.size()) {
      apf::ModelEntity* em = m->toModel(es_up.at(j)); 
      int em_type = m->getModelType(em);
      int em_tag = m->getModelTag(em);
      // Check if the geometry is one dimension higher, tag matches and not
      // in list.
      if ((em_type == d+step) && (em_tag == geom)
           && (findIn(ep, ep->size(), es_up.at(j)) == -1) ) {
        ep->push_back(ep_in->at(i));
        j = up.n;
      }
      j ++;
    }
  }
}

bool vd_chk_ent_v(apf::Mesh2* m, apf::MeshEntity* ent, 
                                std::vector<apf::MeshEntity* >* es_v) {

  assert(es_v->size() > 0);
  int ent_type = m->getType(ent);
  int d = m->typeDimension[ent_type];

  apf::Downward dv;

  m->getDownward(ent, 0, dv);

  std::cout << d << "ent " << ent;
  bool first = false;

  for(int i = 0; i < d+1; i++) {
    int i1 = findIn(es_v, es_v->size(), dv[i]);
    std::cout << i1 << " ";
    if(i1 > -1) {
      if(first) {
        std::cout << "Found second. ";
        std::cout << std::endl;
        return true;
      }
      else {
        first = true;
        std::cout << "Found first. ";
      }
    }
  }
  std::cout << std::endl;

  return false;
}

// Given an element entity set ep1, print lower adjacencies.
void vd_print_down(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in) {

  // Get the dimension of the entity type in the set:
  int ent_type = m->getType(ep_in->at(0));
  int d = m->typeDimension[ent_type];
  printf("Type is %d, dimension is %d\n", ent_type, d);
  apf::Downward down;  

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->size(); i++) {
    int size = m->getDownward(ep_in->at(i), d-1, down);

    if (size > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      printf("The %dth element has surf adjacencies:", i);
      for (int j = 0; j < size; j++) {
        printf(" %p", (void*)down[j]);        
      }
      printf(".\n");
    }
    size = m->getDownward(ep_in->at(i), d-2, down);

    if (size > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      printf("The %dth element has edge adjacencies:", i);
      for (int j = 0; j < size; j++) {
        printf(" %p", (void*)down[j]);        
      }
      printf(".\n");
    }

    size = m->getDownward(ep_in->at(i), d-3, down);

    if (size > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      printf("The %dth element has node adjacencies:", i);
      for (int j = 0; j < size; j++) {
        printf(" %p", (void*)down[j]);        
      }
      printf(".\n");
    }
  }
}


// Given a vertex, print the triangles, their positions and cell memberships.
void vd_print_vert_tri_x(apf::Mesh2* m, apf::MeshEntity* vert) {
  int lookup_v_tri [4] = {2,3,1,0};

  // printf("Field initialized.\n");
  std::vector<apf::MeshEntity*> e_v(0);
  std::vector<apf::MeshEntity*> e_e(0);
  std::vector<apf::MeshEntity*> e_s(0);
  std::vector<apf::MeshEntity*> e_t(0);

  apf::Vector3 t_pos;
  apf::Vector3 v_pos;
  apf::Vector3 dn;
  apf::Vector3 sn;

  m->getPoint(vert, 0, v_pos);

  init_ent_set(&e_v, vert);
  vd_set_up(m, &e_v, &e_e);
  vd_set_up(m, &e_e, &e_s);
  vd_set_up(m, &e_s, &e_t);

  for(int i = 0; i < e_t.size(); i++) {
    apf::Downward d_t;
    apf::Downward d_v;
    m->getDownward(e_t.at(i), 2, d_t);
    m->getDownward(e_t.at(i), 0, d_v);
    int v_id = findIn(d_v, 4, vert);
    assert(v_id > -1);

    apf::MeshEntity* t_curr = d_t[lookup_v_tri[v_id]];
    sn = vd_area_out(m, t_curr, 0);

    t_pos = getLinearCentroid(m, t_curr);
    t_pos = t_pos-v_pos;

      double dist = t_pos*norm_0(t_pos);

      apf::MeshElement* me = createMeshElement(m, e_t.at(i));
      double vol = measure(me);
      destroyMeshElement(me);

      std::cout << "dist " << dist 
                << " v_pos " << v_pos
                << " t_pos " << t_pos
                << std::endl;

      apf::ModelEntity* vert_mdl = m->toModel (vert);
      int vert_type = m->getModelType (vert_mdl);
      int vert_tag = m->getModelTag (vert_mdl);
      std::cout << "vert: " << vert_type << "c" << vert_tag 
                << std::endl;
      vert_mdl = m->toModel (t_curr);
      vert_type = m->getModelType (vert_mdl);
      vert_tag = m->getModelTag (vert_mdl);
      std::cout << "tri: " << vert_type << "c" << vert_tag 
                << std::endl;

  }
}

//----------------------------------------------------

/*
Things to do:
1-)
It is possible to define multiple entity sets. Each have default of 1600 
elements. But it should possible to require more size, and most of the time 
size is a lot less. 
Write a class that will act as a global holder for 'Entity_set's. Delete all
'Entity_set's on scope end.
*/

/*
apf::Copy couldn't be used. Currently, Entity_set, similar to apf:Up is used.
Functions to initialize Entity_set structure from apf::Up and 
apf::MeshEntity. 
I couldn't create proper class definitions. Currently, these are used.
Later, merge into class or structure definitions. 
*/

void copy_ent_set(Entity_set* es, apf::Up up) {
  es->n = up.n;
  for (int i = 0; i < up.n; i ++) {
    es->e[i] = up.e[i];
  }
}

void copy_ent_set(Entity_set* es, apf::Downward down, int size) {
  es->n = size;
  for (int i = 0; i < size; i ++) {
    es->e[i] = down[i];
  }
}

// Copy the content of es_in to es_copy. 
void copy_ent_set(Entity_set* es_copy, Entity_set* es_in) {
  es_copy->n = es_in->n;
  for (int i = 0; i < es_copy->n; i ++) {
    es_copy->e[i] = es_in->e[i];
  }
}

// Copy the content of es_in to es_copy. 
void copy_ent_set(Entity_set* es_copy, std::vector<apf::MeshEntity*>* es_in) {
  es_copy->n = es_in->size();
  for (int i = 0; i < es_copy->n; i ++) {
    es_copy->e[i] = es_in->at(i);
  }
}

void copy_ent_set2(Entity_set* es, apf::MeshEntity** down, int size) {
  es->n = size;
  for (int i = 0; i < size; i ++) {
    es->e[i] = down[i];
  }
}

void init_ent_set(Entity_set* es, apf::MeshEntity* e) {
  es->n = 1;
  es->e[0] = e;
}

void init_ent_set(Entity_set* es, std::vector<apf::MeshEntity*>* e) {
  es->n = e->size();
  for(int i = 0; i < e->size(); i++)
    es->e[i] = e->at(i);
}

// Revert entity set:
void rev_ent_set(Entity_set* es) {
  //printf("Elem to revert: ");
  //for(int i = 0; i < es->n; i++)
  //  printf("%p, ", es->e[i]);

  apf::MeshEntity* e;
  int n = es->n;
  for (int i = 0; i < n/2; i++) {
    e = es->e[n-i-1];
    es->e[n-i-1] = es->e[i];
    es->e[i] = e;
  }
  //printf("Reverted: ");
  //for(int i = 0; i < es->n; i++)
  //  printf("%p, ", es->e[i]);

}

void print_ent_set(apf::Mesh2* m, Entity_set* es) {
  for (int i = 0; i < es->n; i++) {
    int e_type = m->getType(es->e[i]);
    int d = m->typeDimension[e_type];
    int type = m->getModelType(m->toModel(es->e[i]));
    int tag = m->getModelTag(m->toModel(es->e[i]));
    std::cout << d << "-ent " << es->e[i] 
              << " " << type << "c" << tag << std::endl;
  }

}

// Given an element entity set find the number of unique geometric memberships,
// of the same dimension.
void find_elem_set_uniq(apf::Mesh2* m, Entity_set* es, geom_list* g_list) {
  // e_type is from apf::Mesh:Type, 4 is tetrahedron.
  int e_type = m->getType(es->e[0]);
  int d = m->typeDimension[e_type];

  for (int i = 0; i < es->n; i ++) {
    apf::ModelEntity* em = m->toModel(es->e[i]);
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    // Check if geometric dimension is the same as the mesh entity dimension.
    if (em_type == d) {
      // Check if it already exists.
      int geom_reg_flag = 0;
      for (int j = 0; j < g_list->n; j++) {
        if (g_list->geom[j] == em_tag) {
          geom_reg_flag = 1;
        }
      }
      if (geom_reg_flag == 0) {
        g_list->geom[g_list->n] = em_tag;
        g_list->n++;
      }
    }
  }
}


// Given a mesh object, create element and element volume membership tags.
// findNumbering is problematic, so destroy and rewrite.
void vd_tag_mesh(apf::Mesh2* m) {
  int meshDimension = m->getDimension();
  apf::MeshIterator* it = m->begin(meshDimension);

  apf::Numbering* tagnumbering;
  if (m->findNumbering("grain_memb")) {
    //printf("grain_memb tag already exists.\n");
    tagnumbering = m->findNumbering("grain_memb");
    apf::destroyNumbering(tagnumbering);
    tagnumbering = apf::createNumbering(m, "grain_memb", apf::getConstant(meshDimension),1);
  }
  else {
    tagnumbering = apf::createNumbering(m, "grain_memb", apf::getConstant(meshDimension),1);
  }
  apf::Numbering* tagnumbering_ID;
  if (m->findNumbering("elem_ID")) {
    tagnumbering_ID = m->findNumbering("elem_ID");
    apf::destroyNumbering(tagnumbering_ID);
    tagnumbering_ID = apf::createNumbering(m, "elem_ID", apf::getConstant(meshDimension),1);
  }
  else {
    tagnumbering_ID = apf::createNumbering(m, "elem_ID", apf::getConstant(meshDimension),1);
  }

  apf::MeshEntity* e;

  int ID = 1;
  while ((e = m->iterate(it))) 
  {
    apf::number(tagnumbering, e, 0, 0, m->getModelTag(m->toModel(e)));
    apf::number(tagnumbering_ID, e, 0, 0, ID);
    ID++;
  }
  m->end(it);
}

// Remove numberings attached to a mesh.
void vd_rem_tag(apf::Mesh2* m, const char* set_name) {
  apf::Numbering* tagnumbering = m->findNumbering(set_name);
  if(tagnumbering)
    m->removeNumbering(tagnumbering);
}

// Remove numberings attached to a mesh.
void vd_rem_tag(apf::Mesh2* m) {

  while(m->countNumberings() > 0) {
    apf::Numbering* f = m->getNumbering(0);
    std::string n(getName(f));
    //std::cout<<"Name of the numbering being removed is "<< n<< std::endl;
    m->removeNumbering(f);
  }
}

void vd_print_numberings(apf::Mesh2* m) {
  for(int i = 0; i < m->countNumberings(); i++)
    std::cout << apf::getName(m->getNumbering(i)) << std::endl;
}

// To check if deformation causes geometric dissassociation. Given a mesh, 
// tag vertices based on their geometric membership, 
// i.e 0 for quad, 1 for trip, 2 for surf, 3 for internal.
void vd_tag_vert(apf::Mesh2* m) {
  apf::MeshIterator* it = m->begin(0);

  apf::Numbering* tagnumbering;
  if (m->findNumbering("ageom")) {
    printf("ageom tag already exists.\n");
    tagnumbering = m->findNumbering("ageom");
    apf::destroyNumbering(tagnumbering);
    tagnumbering = apf::createNumbering(m, "ageom", m->getShape(),1);
  }
  else {
    tagnumbering = apf::createNumbering(m, "ageom", m->getShape(),1);
  }

  apf::MeshEntity* e;

  while ((e = m->iterate(it))) 
  {
    apf::number(tagnumbering, e, 0, 0, m->getModelType(m->toModel(e)));
  }
  m->end(it);
}

// Given an entity set, create set tags.
// findNumbering is problematic, so destroy and rewrite.
void vd_tag_set(apf::Mesh2* m, Entity_set* es, const char* set_name) {

  int ent_type = m->getType(es->e[0]);
  int d = m->typeDimension[ent_type];

  apf::Numbering* setnumbering;
  if (m->findNumbering(set_name)) {
    std::cout<< set_name << " tag already exists" << std::endl;
    setnumbering = m->findNumbering(set_name);
    apf::destroyNumbering(setnumbering);
    setnumbering = apf::createNumbering(m, set_name, apf::getConstant(d),1);
  }
  else {
    setnumbering = apf::createNumbering(m, set_name, apf::getConstant(d),1);
  }

  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(d);

  std::cout<< "Writing tag " << set_name << "...";
  while ((e = m->iterate(it))) 
  {
    if ((es->n == 0) || (findIn(es->e, es->n, e) == -1)) {
      apf::number(setnumbering, e, 0, 0, 0);
    }
    else {
      apf::number(setnumbering, e, 0, 0, 1);
    }
  }
  m->end(it);

  std::cout<< " Written." << std::endl;

}


// Given a mesh object, print triple edges.
void vd_print_trip(apf::Mesh2* m) {
  int dim = 1;
  apf::MeshIterator* it_e = m->begin(dim);

  apf::MeshEntity* edge;
  while ((edge = m->iterate(it_e))) {
    apf::ModelEntity* em = m->toModel(edge);
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 1) {
      printf("Edge with pointer %p is a triple edge of geometric edge "
             "%d.\n", (void*)edge, em_tag);
    }
  }
  m->end(it_e);
}



// Check if the edge is bounded by two internal vertices, one internal and one 
// external or two external vertices. Return the vertex index of the internal 
// vertex if one is external (return 2, 0 or 1, -1).
// Store in ep*.
int vd_chk_edge_int(apf::Mesh2* m, apf::MeshEntity* edge
                                  , apf::ModelEntity* cell) {
  apf::Downward down;
  int downwardCount = m->getDownward(edge, 0, down);

  int c_type1 = m->getModelType(m->toModel(down[0]));
  int c_tag1 = m->getModelTag(m->toModel(down[0]));
  int c_type2 = m->getModelType(m->toModel(down[1]));
  int c_tag2 = m->getModelTag(m->toModel(down[1]));

  std::cout << c_type1 << "cell" << c_tag1 << " - "
            << c_type2 << "cell" << c_tag2 << std::endl; 
  if (cell == m->toModel(down[0])) {
    if (cell == m->toModel(down[1])) {
      return 2;
    }
    else {
      return 0;
    }
  }
  else {
    if (cell == m->toModel(down[1])) {
      return 1;
    }
    else {
      return -1;
    }
  }
}

// Print the entities with their cell memberships.
void vd_print_ent(apf::Mesh2* m) {

  //std::vector<apf::MeshEntity*> elem;
  apf::MeshIterator* it_e = m->begin(3);

  apf::MeshEntity* e;
  int res_sz = 0;
  while (e = m->iterate(it_e)) {
    res_sz++;
  }
  m->end(it_e);
  printf("elem:%d, ", res_sz);
  res_sz = 0;
  it_e = m->begin(2);
  while (e = m->iterate(it_e)) {
    res_sz++;
  }
  m->end(it_e);
  printf("surf:%d, ", res_sz);
  res_sz = 0;
  it_e = m->begin(1);
  while (e = m->iterate(it_e)) {
    res_sz++;
  }
  m->end(it_e);
  printf("edge:%d, ", res_sz);
  res_sz = 0;
  it_e = m->begin(0);
  while (e = m->iterate(it_e)) {
    res_sz++;
  }
  m->end(it_e);
  printf("vert:%d, ", res_sz);

  apf::Downward down;

  it_e = m->begin(3);
  while (e = m->iterate(it_e)) {
    apf::ModelEntity* cell = m->toModel(e);
    int cell_type = m->getModelType(cell);
    int cell_tag = m->getModelTag(cell);
    printf("Element %p, %dcell%d, ", (void*)e, cell_type, cell_tag);
    m->getDownward(e, 0, down);
    printf("Vertices: %p %p %p %p.\n", (void*)down[0], 
                              (void*)down[1], (void*)down[2], (void*)down[3]);
  }
  m->end(it_e);

  it_e = m->begin(2);
  while (e = m->iterate(it_e)) {
    apf::ModelEntity* cell = m->toModel(e);
    int cell_type = m->getModelType(cell);
    int cell_tag = m->getModelTag(cell);
    printf("Triangle %p, %dcell%d, ", (void*)e, cell_type, cell_tag);
    m->getDownward(e, 0, down);
    printf("Vertices: %p %p %p.\n", 
                              (void*)down[0], (void*)down[1], (void*)down[2]);
  }
  m->end(it_e);

  it_e = m->begin(1);
  while (e = m->iterate(it_e)) {
    apf::ModelEntity* cell = m->toModel(e);
    int cell_type = m->getModelType(cell);
    int cell_tag = m->getModelTag(cell);
    printf("Edge %p, %dcell%d, ", (void*)e, cell_type, cell_tag);
    m->getDownward(e, 0, down);
    printf("Vertices: %p %p.\n", (void*)down[0], (void*)down[1]);
  }
  m->end(it_e);

  apf::Vector3 vert_pos;
  it_e = m->begin(0);
  while (e = m->iterate(it_e)) {
    apf::ModelEntity* cell = m->toModel(e);
    int cell_type = m->getModelType(cell);
    int cell_tag = m->getModelTag(cell);
    printf("Vert %p, %dcell%d", (void*)e, cell_type, cell_tag);
    m->getPoint(e, 0, vert_pos);
    std::cout << " " << vert_pos << std::endl;
  }
  m->end(it_e);

/*
  int cell_type = m->getModelType(em);
  int cell_tag = m->getModelTag(em);

  if (cell == m->toModel(e)) {
    ep->e[cell_size] = e;
    cell_size++;
  }
*/

}

// Find the vertices bounding the given cell.
void vd_find_cell_v(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev, int cell_id, int dim) {
  std::vector<apf::MeshEntity*> es(0);

  vd_find_ent_geom(m, &es, cell_id, dim);
  //vd_find_ent_geom(m, &ev_d, cell_id, dim, 0);
  vd_vert_set(m, &es, ev);
}

// Given a set of simplicial entities, remove the ones belonging to the 
// cell.
void vd_rem_cell(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev, int cell_id, int dim) {
  std::vector<apf::MeshEntity*> ev_d(0);
  ev_d.reserve(ev->size());

  apf::ModelEntity* cell = m->findModelEntity(dim, cell_id);

  for(int i = 0; i < ev->size(); i++) {
    if (!(cell == m->toModel(ev->at(i)))) {
      //printf("Vertex %p, cell %p, %dcell%d\n", e, m->toModel(e), 
      //        m->getModelType(m->toModel(e)), m->getModelTag(m->toModel(e)));
      ev_d.push_back(ev->at(i));
    }
  }
  *ev = ev_d;
}

// Find the vertices bounding the given cell.
void vd_find_bound_v(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev, int cell_id, int dim) {
/*
  std::vector<apf::MeshEntity*> es;
  std::vector<apf::MeshEntity*> ev_d;

  vd_find_ent_geom(m, &es, cell_id, dim);
  vd_find_ent_geom(m, &ev_d, cell_id, dim, 0);
  vd_vert_set(m, &es, ev);
  vd_remove_set(ev, &ev_d);
*/
  vd_find_cell_v(m, ev, cell_id, dim);
  vd_rem_cell(m, ev, cell_id, dim);
}

// Find the dim_e dimensional entities member of a certain cell and dimension.
// Store in ep*. Default value of dim_e is the same as dim.
void vd_find_ent_geom(apf::Mesh2* m, Entity_set* ep, int cell_id, int dim, 
                                            int dim_e) {

  if (dim_e == -1)
    dim_e = dim;
  //printf("Checking for %dcell%d", dim, cell_id);

  int cell_size = 0;
  apf::ModelEntity* cell = m->findModelEntity(dim, cell_id);

  apf::MeshIterator* it_e = m->begin(dim_e);

  apf::MeshEntity* e;
  while (e = m->iterate(it_e)) {

    if (cell == m->toModel(e)) {
      //printf("Vertex %p, cell %p, %dcell%d\n", e, m->toModel(e), 
      //        m->getModelType(m->toModel(e)), m->getModelTag(m->toModel(e)));
      ep->e[cell_size] = e;
      cell_size++;
    }
  }
  m->end(it_e);
  ep->n = cell_size;

}

// Find the edges member of a certain cell and dimension.
// Store in ep*.
void vd_find_edge_topo(apf::Mesh2* m, Entity_set* ep, int cell_id, int dim) {

  int cell_size = 0;
  apf::ModelEntity* cell = m->findModelEntity(dim, cell_id);

  apf::MeshIterator* it_e = m->begin(1);

  apf::MeshEntity* e;
  while (e = m->iterate(it_e)) {

    if (cell == m->toModel(e)) {
      ep->e[cell_size] = e;
      cell_size++;
    }
  }
  m->end(it_e);
  ep->n = cell_size;

}

// Find the edges bounded by two vertices. If does not exist return false.
bool vd_find_edge(apf::Mesh2* m, apf::MeshEntity** vert, 
                                 apf::MeshEntity** edge) {

  apf::Up up_1;
  apf::Up up_2;

  m->getUp(vert[0], up_1);
  m->getUp(vert[1], up_2);

  for (int i = 0; i < up_1.n; i++) {
    for (int j = 0; j < up_2.n; j++) {
      if(up_1.e[i] == up_2.e[j]) {
        *edge = up_1.e[i];
        return true;
      }
    }
  }
  return false;
}

bool vd_find_edge(apf::Mesh2* m, apf::MeshEntity* v1, apf::MeshEntity* v2, 
                                 apf::MeshEntity** edge) {

  apf::Up up_1;
  apf::Up up_2;

  m->getUp(v1, up_1);
  m->getUp(v2, up_2);

  for (int i = 0; i < up_1.n; i++) {
    for (int j = 0; j < up_2.n; j++) {
      if(up_1.e[i] == up_2.e[j]) {
        *edge = up_1.e[i];
        return true;
      }
    }
  }
  return false;
}

apf::MeshEntity* vd_find_esjv(apf::Mesh2* m, apf::MeshEntity* edge, apf::MeshEntity* tri, apf::MeshEntity* v) {
  apf::Downward d_v;
  apf::Downward d_e;

  m->getDownward(tri, 0, d_v);
  m->getDownward(tri, 1, d_e);
  int v1 = findIn(d_v, 3, v);
  if(v1 == -1)
    return NULL;

  int e1 = vd_lookup_tri_ed[v1][0];
  int e2 = vd_lookup_tri_ed[v1][1];
  if(d_e[e1] == edge)
    return d_e[e2];
  else if(d_e[e2] == edge)
    return d_e[e1];
  return NULL;
}


apf::MeshEntity* vd_find_vert(apf::Mesh2* m, apf::MeshEntity* edge1,
                                             apf::MeshEntity* edge2) {
  apf::Downward d1;
  apf::Downward d2;
  assert(edge1 != edge2);
  m->getDownward(edge1, 0, d1);
  m->getDownward(edge2, 0, d2);

  if(d1[0] == d2[0]) {
    return d1[0];
  }
  else if(d1[1] == d2[0]) {
    return d1[1];
  }
  else if(d1[0] == d2[1]) {
    return d1[0];
  }
  else if(d1[1] == d2[1]) {
    return d1[1];
  }
  else
    return NULL;
}

// LOOK
// Given an entity set, find the entities member of a certain geometry.
// Store in ep*.
void vd_ext_ent_geom(apf::Mesh2* m, const Entity_set* ep_in, Entity_set* ep, int geom) {

  int ent_type = m->getType(ep_in->e[0]);
  int d = m->typeDimension[ent_type];

  int cell_size = 0;
  for (int i = 0; i < ep_in->n; i++) {
    apf::ModelEntity* em = m->toModel(ep_in->e[i]);
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);

    // If geometry is not specified, just check dimension. 
    if ((geom == 0) && (d == em_type)) {
      ep->e[cell_size] = ep_in->e[i];
      cell_size++;
    }
    else if ((geom == em_tag) && (d == em_type)) {
      ep->e[cell_size] = ep_in->e[i];
      cell_size++;
    }
  }
  ep->n = cell_size;
}

// LOOK: This could be hardwired to vd_calc_disp by a continue statement.
// Check if given surface, member of a surface geometry is on domain geometry
// surface.
bool vd_chk_surf_dom2(apf::Mesh2* m, apf::MeshEntity* b_surf) {
  apf::ModelEntity* em = m->toModel(b_surf);
  int em_type = m->getModelType(em);
  int em_tag = m->getModelTag(em);

  // printf("Checking if the surface (type %d) is on domain...\n", em_type);
  gmi_model* mdl = m->getModel();
  gmi_set* ent_adj;
  gmi_ent* m_surf = gmi_find(mdl, 2, em_tag);
  ent_adj = gmi_adjacent(mdl, m_surf, 3);
  // printf("%d number of adjacencies obtained.\n",ent_adj->n);

  if (ent_adj->n == 1) {
    // printf("On domain.\n");
    return true;
  }
  else {
    // printf("Not on domain.\n");
    return false;
  }

}

// Check if given surface, member of a surface geometry is on domain geometry
// surface. Works mesh level by checking the number of element adjacencies.
// A single adjacency signief boundary triangle.
bool vd_chk_surf_dom(apf::Mesh2* m, apf::MeshEntity* b_surf) {
  apf::Up up;
  m->getUp(b_surf, up);

  if (up.n == 1) {
    //printf("%p on domain.\n", b_surf);
    return true;
  }
  return false;

}

bool vd_chk_surf_dom(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es_surf) {
  for(int i = 0; i < es_surf->size(); i++) {
    // printf("%d: ",i);
    if(vd_chk_surf_dom(m, es_surf->at(i)))
      return true;
  }
}

// Check if given vertex is on domain boundary.
bool vd_chk_vert_dom(apf::Mesh2* m, apf::MeshEntity* vert) {

  apf::Up up;
  std::vector<apf::MeshEntity*> es_edge(0);
  std::vector<apf::MeshEntity*> es_surf(0);

  vd_set_up(m, vert, &es_edge);
  vd_set_up(m, &es_edge, &es_surf);

  // printf("Sets copied... Number of surface adjacencies is %d.\n", es_surf.n);

  // If any surface is on domain, then on domain.
  for(int i = 0; i < es_surf.size(); i++) {
    // printf("%d: ",i);
    if(vd_chk_surf_dom(m, es_surf.at(i)))
      return true;
  }
  return false;
}

// Check if given entity(vertex, edge, surface) is on domain boundary.
// LOOK: Doesn't work correctly for vertices.
bool vd_chk_ent_dom(apf::Mesh2* m, apf::MeshEntity* e) {
  apf::ModelEntity* em = m->toModel(e);
  int em_type = m->getModelType(em);
  int em_tag = m->getModelTag(em);
  // printf("Checking geometry type %d...\n", em_type);
  // Check geometric membership.
  // If surface, check volumetric geometric adjacencies. If there is one, 
  // on domain; if two, not on domain. If lower than surface, obtain entity
  // sets of surface adjacencies. Check all surfaces for the same.

  if (em_type == 3) {
    return false;
  }

  else if(em_type == 2) {
    return vd_chk_surf_dom(m, e);
  }

  else if(em_type == 1) {
    apf::Up up;
    Entity_set es_surf = Entity_set();

    m->getUp(e, up);
    copy_ent_set(&es_surf, up);
    // printf("Sets copied... Number of surface adjacencies is %d.\n", es_surf.n);
    // If any surface is on domain, then on domain.
    for(int i = 0; i < es_surf.n; i++) {
      em = m->toModel(es_surf.e[i]);
      em_type = m->getModelType(em);
      if (em_type == 3) {
      }
      else if(vd_chk_surf_dom(m, es_surf.e[i]))
        return true;
    }
    return false;
  }

  else {
    // em_type == 0
    apf::Up up;
    Entity_set es_edge = Entity_set();
    Entity_set es_surf = Entity_set();

    m->getUp(e, up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    // printf("Sets copied... Number of surface adjacencies is %d.\n", es_surf.n);

    // If any surface is on domain, then on domain.
    for(int i = 0; i < es_surf.n; i++) {
      // printf("%d: ",i);
      em = m->toModel(es_surf.e[i]);
      em_type = m->getModelType(em);
      if (em_type == 3) {
      }
      else if(vd_chk_surf_dom(m, es_surf.e[i]))
        return true;
    }
    return false;
  }

}


// Given a mesh object, extract triple edges.
void vd_ext_trip(apf::Mesh2* m, Entity_set* es_edge) {
  int dim = 1;

  apf::MeshEntity* edge;
  int trip_index = 0;
  apf::MeshIterator* it_e = m->begin(dim);
  while ((edge = m->iterate(it_e))) {
    apf::ModelEntity* em = m->toModel(edge);
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 1) {
      es_edge->e[trip_index] = edge;
      trip_index++;
    }
  }
  m->end(it_e);

  es_edge->n = trip_index;
  printf("There are %d triple edges.\n", trip_index);
  //int ent_type1 = m->getType(edge);
  //int ent_type2 = m->getType(es_edge->e[0]);
  //printf("The entity type is %d and the entity type in the set is %d.\n", ent_type1, ent_type2);
}

// Given a mesh object, extract boundary and domain surfaces.
void vd_ext_bsurf(apf::Mesh2* m, Entity_set* es_surf, Entity_set* es_dom) {
  int dim = 2;
  apf::MeshIterator* it_s = m->begin(dim);

  apf::Up up;
  apf::MeshEntity* surf;
  int surf_index = 0;
  int dom_index = 0;
  while ((surf = m->iterate(it_s))) {
    apf::ModelEntity* em = m->toModel(surf);
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 2) {
      m->getUp(surf, up);
      if (up.n == 1) {
        es_dom->e[dom_index] = surf;
        dom_index ++;
      }
      else if (up.n == 2) {
        es_surf->e[surf_index] = surf;
        surf_index ++;
      }
      else {
        printf("Surface %p has inappropariate upper adjacency.\n",(void*)surf);
      }
    }
  }
  m->end(it_s);
  es_surf->n = surf_index;
  es_dom->n = dom_index;
  printf("There are %d domain boundary triangles.\n",es_dom->n);
  printf("There are %d grain boundary triangles.\n",es_surf->n);
}

// LOOK: This is not necessary to find quadruple points. Still, this can be used
// as a template.
// Given a mesh object, and domain surface list extract quadruple points, which
// lie on the interior, not on the domain boundary. For now, assume quadruple
// points indeed are members of at least four different grains.
void vd_ext_quad(apf::Mesh2* m, Entity_set* es_quad) {
  int dim = 0;
  apf::MeshIterator* it = m->begin(dim);

  apf::Up up;
  apf::MeshEntity* vert;
  int quad_index = 0;
  Entity_set es_edge;
  Entity_set es_surf;
  Entity_set es_elem;

  while ((vert = m->iterate(it))) {
    // Going over vertices, find those at geometric vertices.
    apf::ModelEntity* em = m->toModel(vert);
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 0) {
      // es_edge = Entity_set();
      // es_surf = Entity_set();
      // es_elem = Entity_set();
      m->getUp(vert, up);

      // printf("Create up, size is %d, add to edge.\n", up.n);
      copy_ent_set(&es_edge, up);
      // printf("Edge of size %d. Setting surf up:\n", es_edge.n);
      vd_set_up(m, &es_edge, &es_surf);
      // printf("Set up to surf of size %d.\n", es_surf.n);
      // printf("Setting elem up:\n");
      vd_set_up(m, &es_surf, &es_elem);
      geom_list g_list = geom_list();
      find_elem_set_uniq(m, &es_elem, &g_list);
      // printf("Elem set, there are %d unique geometries.\n", uniq_geom);

      if (g_list.n > 3) {
        es_quad->e[quad_index] = vert;
        quad_index ++;
      }
    }
  }
  m->end(it);
  es_quad->n = quad_index;
  printf("Number of quadruple points is %d.\n",es_quad->n);
}

// Given a vertex MeshEntity, find the surface couples that bound geometries 
// around that vertex, joined by edges adjacent to that vertex. PLACEHOLDER
void vd_ext_surf_coup(apf::Mesh2* m, apf::MeshEntity* vert) {
}

// Print the surface geometric adjacencies of given dimension adj_dim.
void vd_surf_adj(gmi_model* m, int ent_dim, int adj_dim) {
  // For the surfaces in the model object, print the names of the adjacent
  // volumes.
  gmi_iter* it_s;
  gmi_ent* m_surf;
  it_s = gmi_begin(m, ent_dim);

  gmi_set* ent_adj;
  int i = 0;
  while (m_surf = gmi_next(m, it_s)) {
    ent_adj = gmi_adjacent(m, m_surf, adj_dim);
    // printf("The %d dimensional entity %d has the %d dimensional adjacencies", ent_dim, gmi_tag(m, m_surf), adj_dim);
    for (i = 0; i < ent_adj->n; i ++) {
      printf(", %d", gmi_tag(m, ent_adj->e[i]));
    }
    printf(".\n");
  }
  gmi_end(m, it_s);
}


// Templates for operations that relate and operate on mesh and model entities:

int vd_mesh_mdl_adj(apf::Mesh2* m, apf::MeshEntity* e) {
  // Given a mesh entity, find the geometry attached, print type and tag and
  // return the geometry tag.
  apf::ModelEntity* em = m->toModel(e);
  int em_type = m->getModelType(em);
  int em_tag = m->getModelTag(em);
  // printf("The %d dimensional entity has the tag %d.\n", em_type, em_tag);
  return em_tag;
}

void vd_ext_do_mesh(apf::Mesh2* m) {
  // Given a mesh, move over the vertex entities, find geom memb, 
  // print types and tags. Strangely, it goes over all entities. 
  // find virtual toModel decleration.
  apf::MeshIterator* it_s;
  int mesh_dim = 0;

  apf::MeshEntity* e;
  it_s = m->begin(mesh_dim);

  while (e = m->iterate(it_s)) {
    vd_mesh_mdl_adj(m, e);
  }
  m->end(it_s);
}

// Print the adjacencies of an entity, of a distance n.
// This simply lists upper and lower adjacencies of an entity, and calls their
// adjacencies of distance n-1, recursively.
void vd_ent_deg_adj(apf::Mesh2* m, apf::MeshEntity* ent, int n) {
  int ent_type = m->getType(ent);
  int d = m->typeDimension[ent_type];
  apf::ModelEntity* mdl_nv = m->toModel(ent);
  int type2 = m->getModelType(mdl_nv);
  int tag2 = m->getModelTag(mdl_nv);
  printf("%d ent %p (%dcell%d). \n", d, (void*)ent, type2, tag2);

  if (n > 0) {
    apf::Up up;
    apf::Downward down;
    if (d < 3) {
      m->getUp(ent, up);
      int upwardCount = up.n;
      printf("%d upward adj. \n", upwardCount);
      for (int i = 0; i< up.n ; i++)
        vd_ent_deg_adj(m, up.e[i], n-1);
    }
    if (d > 0) {
      int downwardCount = m->getDownward(ent, d-1, down);
      printf("%d downward adj. \n", downwardCount);
      for (int i = 0; i< downwardCount ; i++)
        vd_ent_deg_adj(m, down[i], n-1);
    }
  }
}

// Given a mesh, going over the entities, print the lower and upper adjacencies
// and their type, whenever applicable.
void vd_mesh_up_adj(apf::Mesh2* m) {

  apf::MeshEntity* e;
  apf::MeshIterator* it;
  int meshDimension = m->getDimension();
  for (int d = 0; d <= meshDimension; ++d) {
    int ID = 1;
    it = m->begin(d);
    while ((e = m->iterate(it))) {
      apf::Up up;
      apf::Downward down;
      if (d < 3) {
        m->getUp(e, up);
        int upwardCount = up.n;
        printf("Mesh element of dimension: %d ID: %d has %d upward adjacencies. \n", d, ID, upwardCount);      
      }
      if (d > 0) {
        int downwardCount = m->getDownward(e, d-1, down);
        printf("Mesh element of dimension: %d ID: %d has %d downward adjacencies. \n", d, ID, downwardCount);      
      }
      ID++;
    }
    m->end(it);
  }
}

// Given a mesh, going over the entities, print the lower and upper adjacencies
// if they are not valid.
void vd_mesh_bad_adj(apf::Mesh2* m) {

  apf::MeshEntity* e;
  apf::MeshIterator* it;
  it = m->begin(1);
  while ((e = m->iterate(it))) {
    apf::Downward down;
    int downwardCount = m->getDownward(e, 0, down);
    if (downwardCount != 2) {
      printf("Edge %p: %d down adj!\n", (void*)e, downwardCount);
      std::cout << vd_get_center(m,e) << std::endl;
    }
  }
  m->end(it);

  it = m->begin(2);
  while ((e = m->iterate(it))) {
    apf::Up up;
    apf::Downward down;
    m->getUp(e, up);
    int upwardCount = up.n;
    if ((upwardCount > 2) or (upwardCount < 1)) {
      printf("Triangle %p: %d upward adj!\n", (void*)e, upwardCount);
      for (int i = 0; i < upwardCount; i++)
        printf("%p, ", (void*)up.e[i]);
      std::cout << vd_get_center(m,e) << std::endl;
    }
    int downwardCount = m->getDownward(e, 1, down);
    if (downwardCount != 3) {
      printf("Triangle %p: %d down adj!\n", (void*)e, downwardCount);
      std::cout << vd_get_center(m,e) << std::endl;
    }
  }
  m->end(it);
  it = m->begin(3);
  while ((e = m->iterate(it))) {
    apf::Downward down;
    int downwardCount = m->getDownward(e, 2, down);
    if (downwardCount != 4) {
      printf("Tetrahedron %p: %d down adj!\n", (void*)e, downwardCount);
      std::cout << vd_get_center(m,e) << std::endl;
    }
  }
  m->end(it);
}

// Given a mesh, and an entity, print the lower and upper adjacencies
// if they are not valid. Return true if not valid.
bool vd_mesh_bad_ent(apf::Mesh2* m, apf::MeshEntity* e) {
  int ent_type = m->getType(e);
  int d = m->typeDimension[ent_type];

  if (d == 1) {
    apf::Downward down;
    int downwardCount = m->getDownward(e, 0, down);
    if (downwardCount != 2) {
      printf("Edge %p: %d down adj!\n", (void*)e, downwardCount);
      return true;
    }
  }
  else if (d == 2) {
    apf::Up up;
    m->getUp(e, up);
    int upwardCount = up.n;
    if ((upwardCount > 2) or (upwardCount < 1)) {
      printf("Triangle %p: %d upward adj!\n", (void*)e, upwardCount);
      return true;
    }

    apf::Downward down;
    int downwardCount = m->getDownward(e, 1, down);
    if (downwardCount != 3) {
      printf("Triangle %p: %d down adj!\n", (void*)e, downwardCount);
      return true;
    }
  }
  else if (d == 3) {
    apf::Downward down;
    int downwardCount = m->getDownward(e, 2, down);
    if (downwardCount != 4) {
      printf("Tetrahedron %p: %d down adj!\n", (void*)e, downwardCount);
      return true;
    }

  }
  return false;

}


// Given a mesh entity, print the number of upper adjacencies and their type.
void vd_count_ent_up(apf::Mesh2* m, apf::MeshEntity* e) {
  apf::MeshEntity* eu;
  int eu_n = m->countUpward(e);
  eu = m->getUpward(e, 0);
  // printf("The number of upper adjacencies is %d and the first member has type %d.\n", eu_n, m->getType(eu));
}


// Given two entity sets, remove the elements of the second list from the 
// first or write the common elements into the first.
// swtch = 0, intersect; swtch = 1, remove.
// TODO LOOK somehow, the first entry of the set of the removed entities is not 
// matched. 
void vd_remove_set(Entity_set* ep1, Entity_set* ep2) {

  int rem_size = 0;
  //for (int i = 0; i < ep2->n; i++) {
    // printf("%p, ", ep2->e[i]);
  //}
  // printf("\n");

  for (int i = 0; i < ep1->n; i++) {
    int i1 = findIn(ep2->e, ep2->n, ep1->e[i]);
    // printf("%p, ", ep1->e[i]);

    if ( i1 >= 0) {
      rem_size++;
      // printf("removed, ");
    }
    else {
      ep1->e[i-rem_size] = ep1->e[i];
      // printf("remains, ");
    }
    // printf("\n");
  }

  // printf("Set size:%d, remove size:%d, ", ep1->n, rem_size);
  ep1->n = ep1->n - rem_size;
  // printf("New size%d.\n", ep1->n);
}

void vd_inter_set(Entity_set* ep1, Entity_set* ep2) {

  int rem_size = 0;

  for (int i = 0; i < ep1->n; i++) {
    int i1 = findIn(ep2->e, ep2->n, ep1->e[i]);
    if ( i1 < 0) {
      rem_size++;
    }
    else {
      ep1->e[i-rem_size] = ep1->e[i];
    }
  }

  printf("Set size:%d, remove size:%d, ", ep1->n, rem_size);
  ep1->n = ep1->n - rem_size;
  printf("New size%d.\n", ep1->n);
}

// Given two up adjacency entity arrays, obtain a list of unique entities.
void vd_merge_set(Entity_set* ep1, Entity_set* ep2) {
  std::vector<apf::MeshEntity*> ep_list(ep2->n);

  int unique = 0;
  int i2 = 0;
  for (int i = 0; i < ep2->n; i ++) {
    // printf("ep1->n is %d.\n", ep1->n);
    i2 = findIn(ep1->e, ep1->n, ep2->e[i]);
    // printf("i2 is %d Pointers: %p %p.\n", i2, ep2.e[i], ep1->e[i2]);
    if ( i2 == -1) {
      ep_list.at(unique) = ep2->e[i];
      unique++;
    }
  }
  // printf("Unique is %d.\n",unique);

  for (int i = 0; i < unique; i ++) {
    ep1->e[i+ep1->n] = ep_list.at(i);
  }
  ep1->n = ep1->n+unique;

}

// Given an entity set ep1, find one level higher adjacencies if possible.
// Store in ep*.
void vd_set_up(apf::Mesh2* m, const Entity_set* ep_in, Entity_set* ep) {

  apf::Up up;
  Entity_set es;
  es.n = 0;

  ep->n = 0;

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->n; i++) {
    m->getUp(ep_in->e[i], up);

    if (up.n > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      copy_ent_set(&es, up);

      vd_merge_set(ep, &es);
    }
  }
}

// Given an entity set ep1, find one level lower adjacencies if possible.
// Store in ep*.
void vd_set_down(apf::Mesh2* m, const Entity_set* ep_in, Entity_set* ep, int level) {

  // Get the dimension of the entity type in the set:
  int ent_type = m->getType(ep_in->e[0]);
  int d = m->typeDimension[ent_type];
  // printf("Type is %d, dimension is %d\n",d, ent_type);
  apf::Downward down;  
  Entity_set es = Entity_set();

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->n; i++) {
    int size = m->getDownward(ep_in->e[i], d-level, down);

    if (size > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      copy_ent_set(&es, down, size);
      vd_merge_set(ep, &es);
    }
  }
}

// For two d-dim entities, starting from d-1 dimensional adjacencies, 
// find the highest dimensional common lower adjacency. If none could be found
// return NULL.
apf::MeshEntity* vd_find_common_adj(apf::Mesh2* m, apf::MeshEntity* e1, 
                                             apf::MeshEntity* e2) {
  int e_t1 = m->getType(e1);
  int e_t2 = m->getType(e2);
  int d1 = m->typeDimension[e_t1];
  int d2 = m->typeDimension[e_t2];

  apf::Downward dv1;
  apf::Downward dv2;
  apf::MeshEntity* common = NULL;
  for(int dim_low = std::min(d1, d2) - 1; dim_low > -1; dim_low--) {
    int sz1 = m->getDownward(e1, dim_low, dv1);
    int sz2 = m->getDownward(e2, dim_low, dv2);

    for(int i = 0; i < sz1; i++) {
      for(int j = 0; j < sz2; j++) {
        if(dv1[i] == dv2[j]) {
          common = dv1[i];
          i = sz1;
          j = sz2;
          dim_low = -1;
        }
      }
    }
  }
  return common;
}

// For two edges, find the joining edge. assert one exists.
apf::MeshEntity* vd_find_tris_edge(apf::Mesh2* m, apf::MeshEntity* tri1, 
                                             apf::MeshEntity* tri2) {
  int e_t1 = m->getType(tri1);
  int e_t2 = m->getType(tri2);
  int d1 = m->typeDimension[e_t1];
  int d2 = m->typeDimension[e_t2];
  assert(d1 == d2 and d1 == 2);

  apf::Downward dv1;
  apf::Downward dv2;
  apf::MeshEntity* common = NULL;
  m->getDownward(tri1, 1, dv1);
  m->getDownward(tri2, 1, dv2);

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      if(dv1[i] == dv2[j]) {
        common = dv1[i];
        i = 3;
        j = 3;
      }
    }
  }
  return common;
}

std::vector<apf::MeshEntity*> vd_get_3c_tet_bound(apf::Mesh2* m, 
                                std::vector<apf::MeshEntity*> & tet_in) {
  std::vector<apf::MeshEntity*> tet_bound(0);
  tet_bound.reserve(tet_in.size());
  apf::Downward down;
  for(int i = 0; i < tet_in.size(); i++) {
    m->getDownward(tet_in.at(i), 0, down);
    for(int j = 0; j < 4; j++) {
      if(m->getModelType(m->toModel(down[j])) < 3) {
        j = 4;
        tet_bound.push_back(tet_in.at(i));
      }
    }
  }
  return tet_bound;
}

/*
void vd_set_down(apf::Mesh2* m, const std::vector<apf::MeshEntity*>* ep_in, std::vector<apf::MeshEntity*>* ep, int level) {

  // Get the dimension of the entity type in the set:
  assert(ep_in->size() > 0);
  int ent_type = m->getType(ep_in->at(0));
  int d = m->typeDimension[ent_type];
  // printf("Type is %d, dimension is %d\n",d, ent_type);
  apf::Downward down;
  Entity_set es = Entity_set();

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->n; i++) {
    int size = m->getDownward(ep_in->e[i], d-level, down);

    if (size > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      copy_ent_set(&es, down, size);
      vd_merge_set(ep, &es);
    }
  }
}
*/

// Given an edge entity set ep_in, extract those, which rest on the grain 
// surface.
// (i.e. a surface entity set, whose adjacent elements belong to geometry n.)
// Store in ep*.
// LOOK write functions that do something simple, and only that!
void vd_chk_edge_grain(apf::Mesh2* m, const Entity_set* ep_in, Entity_set* ep
                                        , int geom) {
  apf::Up up;
  *ep = Entity_set();

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->n; i++) {
    apf::ModelEntity* em = m->toModel(ep_in->e[i]); 
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    // If triple edge or surface edge:
    if ((em_type == 1) || (em_type == 2)) {
      m->getUp(ep_in->e[i], up);
      Entity_set es_surf = Entity_set();
      Entity_set es_elem = Entity_set();
      copy_ent_set(&es_surf, up);
      vd_set_up(m, &es_surf, &es_elem);

      int j = 0;
      while (j < es_elem.n) {
        em = m->toModel(es_elem.e[j]); 
        em_type = m->getModelType(em);
        em_tag = m->getModelTag(em);
        if ( em_tag == geom ) {
          ep->e[ep->n] = ep_in->e[i];
          ep->n++;
          j = es_elem.n;
        }
        j ++;
      }
    }
  }

}


// Given an surface entity set ep_in, extract those, which rest on the grain 
// surface.
// Store in ep*.
// LOOK write functions that do something simple, and only that!
void vd_chk_surf_grain(apf::Mesh2* m, const Entity_set* ep_in, Entity_set* ep
                                        , int geom) {
  apf::Up up;
  *ep = Entity_set();

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->n; i++) {
    apf::ModelEntity* em = m->toModel(ep_in->e[i]); 
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 2) {

      m->getUp(ep_in->e[i], up);

      int j = 0;
      while (j < up.n) {
        em = m->toModel(up.e[j]); 
        em_type = m->getModelType(em);
        em_tag = m->getModelTag(em);
        if ( em_tag == geom ) {
          ep->e[ep->n] = ep_in->e[i];
          ep->n++;
          j = up.n;
        }
        j ++;
      }
    }
  }
}


// Given an entity set ep_in, extract those, one of d dim upper adjacency of 
// which belongs to a certain geometry of that adjacencies own dimension.
// (i.e. a surface entity set, whose adjacent elements belong to geometry n.)
// Store in ep*.
// LOOK write functions that do something simple, and only that!
// So write a function, which extracts edges on a grain surface. 
void vd_chk_set_geom(apf::Mesh2* m, const Entity_set* ep_in, Entity_set* ep
                                        , int geom, int step) {

  // Get the dimension of the entity set elements:
  int ent_type = m->getType(ep_in->e[0]);
  int d = m->typeDimension[ent_type];

  apf::Up up;
  *ep = Entity_set();

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->n; i++) {
    m->getUp(ep_in->e[i], up);
    Entity_set es_up = Entity_set();
    if (step == 1) {
      copy_ent_set(&es_up, up);
    }
    else if (step == 2) {
      Entity_set es_1 = Entity_set();
      copy_ent_set(&es_1, up);
      vd_set_up(m, &es_1, &es_up);
    }
    int j = 0;
    while (j < es_up.n) {
      apf::ModelEntity* em = m->toModel(es_up.e[j]); 
      int em_type = m->getModelType(em);
      int em_tag = m->getModelTag(em);
      // Check if the geometry is one dimension higher, tag matches and not
      // in list.
      if ((em_type == d+step) && (em_tag == geom)
           && (findIn(ep->e, ep->n, es_up.e[j]) == -1) ) {
        ep->e[ep->n] = ep_in->e[i];
        ep->n++;
        j = up.n;
      }
      j ++;
    }
  }

}

// Given an element entity set ep1, print lower adjacencies.
void vd_print_down(apf::Mesh2* m, const Entity_set* ep_in) {

  // Get the dimension of the entity type in the set:
  int ent_type = m->getType(ep_in->e[0]);
  int d = m->typeDimension[ent_type];
  printf("Type is %d, dimension is %d\n", ent_type, d);
  apf::Downward down;  

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < ep_in->n; i++) {
    int size = m->getDownward(ep_in->e[i], d-1, down);

    if (size > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      printf("The %dth element has surf adjacencies:", i);
      for (int j = 0; j < size; j++) {
        printf(" %p", (void*)down[j]);        
      }
      printf(".\n");
    }
    size = m->getDownward(ep_in->e[i], d-2, down);

    if (size > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      printf("The %dth element has edge adjacencies:", i);
      for (int j = 0; j < size; j++) {
        printf(" %p", (void*)down[j]);        
      }
      printf(".\n");
    }

    size = m->getDownward(ep_in->e[i], d-3, down);

    if (size > 0) {
      // Convert apf::Up to Entity_set. USE initialize class.
      printf("The %dth element has node adjacencies:", i);
      for (int j = 0; j < size; j++) {
        printf(" %p", (void*)down[j]);        
      }
      printf(".\n");
    }

  }
}

// Given a mesh entity, print the one dim lower adjacencies.
void vd_print_down(apf::Mesh2* m, apf::MeshEntity* e_in) {

  // Get the dimension of the entity type in the set:
  int ent_type = m->getType(e_in);
  int d = m->typeDimension[ent_type];
  printf("%p, type is %d, dimension is %d\n", (void*)e_in, ent_type, d);
  apf::Downward down;  

  std::cout << e_in << " "
            << m->getModelType(m->toModel(e_in)) << "c"
            << m->getModelTag(m->toModel(e_in))
            << std::endl;

  for(int dim = 0; dim < d; dim++) {

    int size = m->getDownward(e_in, dim, down);

    std::cout << "\t" << dim << "-ents" << std::endl;
    for(int i = 0; i < size; i++) {
      std::cout << "\t\t" << down[i] << " "
                << m->getModelType(m->toModel(down[i])) << "c"
                << m->getModelTag(m->toModel(down[i]))
                << std::endl;
    }
  }
}

std::pair<int,int> vd_cell_pair(apf::Mesh2* m, apf::MeshEntity* ent) {
  apf::ModelEntity* mdl = m->toModel(ent);
  int type = m->getModelType(mdl);
  int tag = m->getModelTag(mdl);
  return std::make_pair(type,tag);
}

// Given a mesh, go over vertices, print number of entities in given upper 
// adjacency sets. 
void print_vert_adj(apf::Mesh2* m) {

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;

  Entity_set es_vert = Entity_set();

  int ID = 0;
  while( (vert = m->iterate(it))) {
    init_ent_set(&es_vert, vert);
    printf("The current vertex entity list has number of entities %d.\n",es_vert.n);
    Entity_set es_edge = Entity_set();
    Entity_set es_surf = Entity_set();
    Entity_set es_elem = Entity_set();
    vd_set_up(m, &es_vert, &es_edge);
    printf("The current edge entity list has number of entities %d.\n",es_edge.n);
    // Surfaces
    vd_set_up(m, &es_edge, &es_surf);
    printf("The current surface entity list has number of entities %d.\n",es_surf.n);
    vd_set_up(m, &es_surf, &es_elem);
    printf("The current element entity list has number of entities %d.\n",es_elem.n);

    for (int i = 0; i < es_surf.n; i++) {
      int geom_memb = vd_mesh_mdl_adj(m, es_surf.e[i]);
      printf("Surface %d has geometric membership %d.\n",i, geom_memb);
    }
    for (int i = 0; i < es_elem.n; i++) {
      int geom_memb = vd_mesh_mdl_adj(m, es_elem.e[i]);
      printf("Surface %d has geometric membership %d.\n",i, geom_memb);
    }
    /*
    printf("Vert ID: %d pointer %p has %d number of surface adjacencies.\n", ID, vert, ep.n);
    printf("The first surface belongs to the volume %d.\n", m->getModelTag(m->toModel(ep.e[0])));
    ID++;
    */
  }
  m->end(it);
}

// Given a mesh, go over vertices, print number of entities in given upper 
// adjacency sets. 
void vd_print_pos_vert(apf::Mesh2* m) {

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;

  Entity_set es_vert = Entity_set();

  int ID = 0;
  apf::Vector3 pnt(0,0,0);
  while( (vert = m->iterate(it))) {
    m->getPoint(vert, 0, pnt);
    std::cout << vert << " pos " 
              << pnt
              << std::endl;
  }
  m->end(it);
}

// Given a set of positions, print them in MATLAB compatible matrix format. 
void vd_print_pos(std::vector<apf::Vector3> &v) {
  std::cout << "v = [";
  for(int i = 0; i < v.size(); i++) {
    std::cout << "   "
              << v.at(i)[0] << "," << v.at(i)[1] << "," << v.at(i)[2]
              << ";..." << std::endl;
  }
  std::cout << "]" << std::endl;
}

// Given a position, print them in MATLAB compatible matrix format. 
void vd_print_pos(apf::Vector3 &v) {
  std::cout << "v = ["
            << v[0] << "," << v[1] << "," << v[2]
            << "]" << std::endl;
}

// Tag the entities around a vertex and save the mesh.
void vd_save_vtk_vert(apf::Mesh2* m, apf::MeshEntity* vert, 
            char const* FileName) {

  // Get all surrounding entities:
  std::vector<apf::MeshEntity*> es_vert(0);
  std::vector<apf::MeshEntity*> es_edge(0);
  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_elem(0);

  init_ent_set(&es_vert, vert);

  // Get the tetrahedra:
  vd_set_up(m, &es_vert, &es_edge);
  vd_set_up(m, &es_edge, &es_surf);
  vd_set_up(m, &es_surf, &es_elem);

  // Save the entities before the transfer.
  vd_rem_tag(m);
  vd_tag_mesh(m);
  vd_tag_set(m, &es_elem, "Vertex_ent");

  apf::writeVtkFiles(FileName, m);
  vd_rem_tag(m);
  vd_tag_mesh(m);
}

// Tag the entities around a vertex and save the mesh.
void vd_save_vtk_vert(apf::Mesh2* m, std::vector<apf::MeshEntity*>* vert, 
            char const* FileName) {

  // Get all surrounding entities:
  std::vector<apf::MeshEntity*> es_vert(0);
  std::vector<apf::MeshEntity*> es_edge(0);
  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_elem(0);

  es_vert = *vert;

  std::cout << "Saving around vertices ";
  for(int i = 0; i < es_vert.size(); i++)
    std::cout << es_vert.at(i) << " ";
  std::cout << std::endl;

  // Get the tetrahedra:
  vd_set_up(m, &es_vert, &es_edge);
  vd_set_up(m, &es_edge, &es_surf);
  vd_set_up(m, &es_surf, &es_elem);

  // Save the entities before the transfer.
  vd_rem_tag(m);
  vd_tag_mesh(m);
  vd_tag_set(m, &es_elem, "Vertex_ent");

  apf::writeVtkFiles(FileName, m);
  vd_rem_tag(m);
  vd_tag_mesh(m);
}

void vd_save_vtk_ent(apf::Mesh2* m, apf::MeshEntity* ent, 
                                          char const* FileName) {

  std::cout << "Saving around entity " << ent << std::endl;

  int e_type = m->getType(ent);
  int d = m->typeDimension[e_type];
  // Get the verts:
  std::vector<apf::MeshEntity*> es_vert(0);
  std::vector<apf::MeshEntity*> es_edge(0);
  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_elem(0);

  es_edge.resize(1);
  es_edge.at(0) = ent;

  vd_set_down(m, &es_edge, &es_vert, d);
  //apf::Downward down;
  //int size = m->getDownward(ent, 0, down);
  //copy_ent_set(&vert, down, size);
  //vd_save_vtk_vert(m, &vert, FileName);
  // Get all surrounding entities:

  std::cout << "Saving around vertices ";
  for(int i = 0; i < es_vert.size(); i++)
    std::cout << es_vert.at(i) << " ";
  std::cout << std::endl;

  // Get the tetrahedra:
  vd_set_up(m, &es_vert, &es_edge);
  vd_set_up(m, &es_edge, &es_surf);
  vd_set_up(m, &es_surf, &es_elem);

  // Save the entities before the transfer.
  vd_rem_tag(m);
  vd_tag_mesh(m);

  vd_tag_set(m, &es_elem, "Vertex_ent");
  if(d == 0) {
    vd_tag_set(m, &es_elem, "ent_tet");
  }
  else if(d == 1) {
    vd_set_up(m, ent, &es_surf);
    vd_set_up(m, &es_surf, &es_elem);
    vd_tag_set(m, &es_elem, "ent_tet");
  }
  else if(d == 2) {
    vd_set_up(m, ent, &es_elem);
    vd_tag_set(m, &es_elem, "ent_tet");
  }
  else {
    es_elem.resize(1);
    es_elem.at(0) = ent;
    vd_tag_set(m, &es_elem, "ent_tet");
  }
  apf::writeVtkFiles(FileName, m);
}

void vd_save_vtk_ent(apf::Mesh2* m, std::vector<apf::MeshEntity* >* ent, 
                                          char const* FileName) {
  if(ent->size() == 0)
    return;
  std::cout << "Saving around entities " << std::endl;

  int e_type = m->getType(ent->at(0));
  int d = m->typeDimension[e_type];
  // Get the verts:
  std::vector<apf::MeshEntity*> es_vert(0);
  std::vector<apf::MeshEntity*> es_edge(0);
  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_elem(0);

  es_edge = *ent;

  vd_set_down(m, &es_edge, &es_vert, d);
  //apf::Downward down;
  //int size = m->getDownward(ent, 0, down);
  //copy_ent_set(&vert, down, size);
  //vd_save_vtk_vert(m, &vert, FileName);
  // Get all surrounding entities:

  std::cout << "Saving around vertices ";
  for(int i = 0; i < es_vert.size(); i++)
    std::cout << es_vert.at(i) << " ";
  std::cout << std::endl;

  // Get the tetrahedra:
  vd_set_up(m, &es_vert, &es_edge);
  vd_set_up(m, &es_edge, &es_surf);
  vd_set_up(m, &es_surf, &es_elem);

  // Save the entities before the transfer.
  vd_rem_tag(m);
  vd_tag_mesh(m);

  vd_tag_set(m, &es_elem, "Vertex_ent");
  if(d == 0) {
    vd_tag_set(m, &es_elem, "ent_tet");
  }
  else if(d == 1) {
    vd_set_up(m, ent, &es_surf);
    vd_set_up(m, &es_surf, &es_elem);
    vd_tag_set(m, &es_elem, "ent_tet");
  }
  else if(d == 2) {
    vd_set_up(m, ent, &es_elem);
    vd_tag_set(m, &es_elem, "ent_tet");
  }
  else {
    es_elem = *ent;
    vd_tag_set(m, &es_elem, "ent_tet");
  }
  apf::writeVtkFiles(FileName, m);
}

// Copy a set of entities from a mesh to a new mesh, and save into a VTK file.
// The entity sets are ordered by dimension, e_set->at(0) contains the vertices,
// e_set->at(1) contains the edges, etc.
void vd_save_vtk_set(apf::Mesh2* m, std::vector<Entity_set>* e_set, 
                  char const* FileName) {
  // Entity map, to keep entity mapping between the meshes.
  std::map<apf::MeshEntity*, apf::MeshEntity*> e_map{};
  apf::MeshEntity* e_new;
  apf::ModelEntity* mdl;

  gmi_register_null();

	apf::Mesh2* m_new = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);

  apf::Vector3 vert_pos;
  apf::Downward down;

  //std::cout<< "Number of vert " << e_set->at(0).n
  //         << "Edge " << e_set->at(1).n
  //         << "Tri " << e_set->at(2).n
  //         << "Tet " << e_set->at(3).n
  //         << std::endl;

  for(int i = 0; i < e_set->at(0).n; i++) {
    mdl = m->toModel(e_set->at(0).e[i]);
    int cell_dim = m->getModelType(mdl);
    int cell_id = m->getModelTag(mdl);
    mdl = m_new->findModelEntity(cell_dim, cell_id);

    e_new = m_new->createVert(mdl);

    e_map[e_set->at(0).e[i]] = e_new;

    m->getPoint(e_set->at(0).e[i], 0, vert_pos);
    //std::cout<< e_set->at(0).e[i] << " " << vert_pos << " " << e_new << std::endl;
    m_new->setPoint(e_new, 0, vert_pos);
  }

  for(int i = 0; i < e_set->at(1).n; i++) {
    mdl = m->toModel(e_set->at(1).e[i]);
    int cell_dim = m->getModelType(mdl);
    int cell_id = m->getModelTag(mdl);
    mdl = m_new->findModelEntity(cell_dim, cell_id);

    m->getDownward(e_set->at(1).e[i],0, down);

    //std::cout<< "Edge vertices" << std::endl;

    for (int j = 0; j < 2; j++) {
      //std::cout<< down[j];
      down[j] = e_map[down[j]];
      //std::cout<< " " << down[j] << std::endl;
    }

    e_new = makeOrFind(m_new, mdl, 1, down, 0);

    //std::cout<< e_set->at(1).e[i] << " " << e_new << std::endl;

    e_map[e_set->at(1).e[i]] = e_new;
  }

  for(int i = 0; i < e_set->at(2).n; i++) {
    mdl = m->toModel(e_set->at(2).e[i]);
    int cell_dim = m->getModelType(mdl);
    int cell_id = m->getModelTag(mdl);
    mdl = m_new->findModelEntity(cell_dim, cell_id);

    m->getDownward(e_set->at(2).e[i], 1, down);

    //std::cout<< "Tri vertices" << std::endl;
    for (int j = 0; j < 3; j++) {
      //std::cout << down[j];
      down[j] = e_map[down[j]];
      //std::cout<< " " << down[j] << std::endl;
    }

    e_new = makeOrFind(m_new, mdl, 2, down, 0);
    e_map[e_set->at(2).e[i]] = e_new;

    //std::cout<< e_set->at(2).e[i] << " " << e_new << std::endl;
  }

  for(int i = 0; i < e_set->at(3).n; i++) {
    mdl = m->toModel(e_set->at(3).e[i]);
    int cell_dim = m->getModelType(mdl);
    int cell_id = m->getModelTag(mdl);
    mdl = m_new->findModelEntity(cell_dim, cell_id);

    m->getDownward(e_set->at(3).e[i], 2, down);
    //std::cout<< "Tet vertices" << std::endl;
    for (int j = 0; j < 4; j++) {
      //std::cout<< down[j];
      down[j] = e_map[down[j]];
      //std::cout<< " " << down[j] << std::endl;
    }

    e_new = makeOrFind(m_new, mdl, 4, down, 0);
    //e_new = buildElement(m_new, mdl, 4, down, 0);
    e_map[e_set->at(3).e[i]] = e_new;

    //std::cout<< e_set->at(3).e[i] << " " << e_new << std::endl;
  }

  Entity_set e_set_new;
  for(int i = 0; i < e_set->at(3).n; i++) {
    e_set_new.e[i] = e_map[e_set->at(3).e[i]];
  }
  e_set_new.n = e_set->at(3).n;

  vd_rem_tag(m_new);
  vd_tag_mesh(m_new);

  vd_tag_set(m_new, &e_set_new, "Vertex_ent");

  apf::writeVtkFiles(FileName, m_new);

  apf::destroyMesh(m_new);

}

// Copy a set of entities from a mesh to a new mesh, and save into a VTK file.
// The entity sets are ordered by dimension, e_set->at(0) contains the vertices,
// e_set->at(1) contains the edges, etc.
void vd_save_vtk_set(apf::Mesh2* m, std::vector< std::vector<apf::MeshEntity* > >* e_set, char const* FileName) {
  // Entity map, to keep entity mapping between the meshes.
  std::map<apf::MeshEntity*, apf::MeshEntity*> e_map{};
  apf::MeshEntity* e_new;
  apf::ModelEntity* mdl;

  gmi_register_null();

	apf::Mesh2* m_new = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);

  apf::Vector3 vert_pos;
  apf::Downward down;

  //std::cout<< "Number of vert " << e_set->at(0).n
  //         << "Edge " << e_set->at(1).n
  //         << "Tri " << e_set->at(2).n
  //         << "Tet " << e_set->at(3).n
  //         << std::endl;

  for(int i = 0; i < e_set->at(0).size(); i++) {
    mdl = m->toModel(e_set->at(0).at(i));
    int cell_dim = m->getModelType(mdl);
    int cell_id = m->getModelTag(mdl);
    mdl = m_new->findModelEntity(cell_dim, cell_id);

    e_new = m_new->createVert(mdl);

    e_map[e_set->at(0).at(i)] = e_new;

    m->getPoint(e_set->at(0).at(i), 0, vert_pos);
    //std::cout<< e_set->at(0).e[i] << " " << vert_pos << " " << e_new << std::endl;
    m_new->setPoint(e_new, 0, vert_pos);
  }

  for(int i = 0; i < e_set->at(1).size(); i++) {
    mdl = m->toModel(e_set->at(1).at(i));
    int cell_dim = m->getModelType(mdl);
    int cell_id = m->getModelTag(mdl);
    mdl = m_new->findModelEntity(cell_dim, cell_id);

    m->getDownward(e_set->at(1).at(i),0, down);

    //std::cout<< "Edge vertices" << std::endl;

    for (int j = 0; j < 2; j++) {
      //std::cout<< down[j];
      down[j] = e_map[down[j]];
      //std::cout<< " " << down[j] << std::endl;
    }

    e_new = makeOrFind(m_new, mdl, 1, down, 0);

    //std::cout<< e_set->at(1).e[i] << " " << e_new << std::endl;

    e_map[e_set->at(1).at(i)] = e_new;
  }

  for(int i = 0; i < e_set->at(2).size(); i++) {
    mdl = m->toModel(e_set->at(2).at(i));
    int cell_dim = m->getModelType(mdl);
    int cell_id = m->getModelTag(mdl);
    mdl = m_new->findModelEntity(cell_dim, cell_id);

    m->getDownward(e_set->at(2).at(i), 1, down);

    //std::cout<< "Tri vertices" << std::endl;
    for (int j = 0; j < 3; j++) {
      //std::cout << down[j];
      down[j] = e_map[down[j]];
      //std::cout<< " " << down[j] << std::endl;
    }

    e_new = makeOrFind(m_new, mdl, 2, down, 0);
    e_map[e_set->at(2).at(i)] = e_new;

    //std::cout<< e_set->at(2).e[i] << " " << e_new << std::endl;
  }

  for(int i = 0; i < e_set->at(3).size(); i++) {
    mdl = m->toModel(e_set->at(3).at(i));
    int cell_dim = m->getModelType(mdl);
    int cell_id = m->getModelTag(mdl);
    mdl = m_new->findModelEntity(cell_dim, cell_id);

    m->getDownward(e_set->at(3).at(i), 2, down);
    //std::cout<< "Tet vertices" << std::endl;
    for (int j = 0; j < 4; j++) {
      //std::cout<< down[j];
      down[j] = e_map[down[j]];
      //std::cout<< " " << down[j] << std::endl;
    }

    e_new = makeOrFind(m_new, mdl, 4, down, 0);
    //e_new = buildElement(m_new, mdl, 4, down, 0);
    e_map[e_set->at(3).at(i)] = e_new;

    //std::cout<< e_set->at(3).e[i] << " " << e_new << std::endl;
  }

  std::vector<apf::MeshEntity* > e_set_new(0);
  e_set_new.resize(e_set->at(3).size());
  for(int i = 0; i < e_set->at(3).size(); i++) {
    e_set_new.at(i) = e_map[e_set->at(3).at(i)];
  }

  vd_rem_tag(m_new);
  vd_tag_mesh(m_new);

  vd_tag_set(m_new, &e_set_new, "Vertex_ent");

  apf::writeVtkFiles(FileName, m_new);

  apf::destroyMesh(m_new);

}

void safe_mkdir(const char* path) {
  mode_t const mode = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST)
  {
    reel_fail("MDS: could not create directory \"%s\"\n", path);
  }
}

void save_vtk_mov(apf::Mesh2* m, int t_in, const char* FileName) {
  vd_rem_tag(m);
  vd_tag_mesh(m);

  // Adding this part to generate a video or gif of the grain growth 
  // simulation in a hackish way.
  std::vector<std::vector<apf::MeshEntity*> > es_set_save
          (4, std::vector<apf::MeshEntity* >(0));

  struct gmi_iter* it;
  struct gmi_ent* e;

  gmi_model* mdl = m->getModel();
  it = gmi_begin(mdl, 3);

  while ((e = gmi_next(mdl, it))) {
    int tag = gmi_tag(mdl, e);
    vd_find_ent_geom(m, &es_set_save.at(3), tag, 3);
    for (int j = 3; j > 0; j--) {
      vd_set_down(m, &es_set_save.at(j), &es_set_save.at(j-1));
    }

    std::stringstream ss;
    ss.precision(2);
    ss << FileName << "cell" << tag << "t"
       << std::fixed << (double)(t_in);
    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    vd_save_vtk_set(m, &es_set_save, cstr);
  }
  gmi_end(mdl, it);

  vd_rem_tag(m);

}

// TODO remove. GDB doesn't print content of map containers. This is a quick 
// solution to print the map value associated with key.
double print_map_i2d(std::map<int, double>& mapd, int key) {
  return mapd[key];
}

