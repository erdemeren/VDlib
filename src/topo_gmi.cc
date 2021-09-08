#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_extinfo.h"
#include "topo_geom.h"

// The vertex order for edges and faces, as defined by MDS.
int vd_edge_vert [6][2] = {{0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}};
int vd_face_vert [4][3] = {{0,2,1}, {0,1,3}, {1,2,3}, {0,3,2}};

// Information retrival from and modifications to geometric model object.
// It should also modify the mesh adjacencies accordingly.

// Create a model object. 


// Create a mesh based on the model object (center point of every surface,
// center point of every grain...)


// Used in gmi_remap: Given an element, extract the geometric membership of 
// entities. See vd_gmi_remap.

void vd_ext_down_geo(apf::Mesh2* m, apf::MeshEntity** down, int size
                                                        , int* geom_type
                                                        , int* geom_tag ) {
  for (int i = 0; i < size; i++) {
    apf::ModelEntity* em = m->toModel(down[i]);
    geom_type[i] = m->getModelType(em);
    geom_tag[i] = m->getModelTag(em);
  }

}

void vd_ext_elem_geo(apf::Mesh2* m, apf::MeshEntity* elem, int* geom_type
                                                         , int* geom_tag) {

  apf::Downward surf;
  apf::Downward edge;
  apf::Downward vert;

  int size = m->getDownward(elem, 2, surf);
  size = m->getDownward(elem, 1, edge);
  size = m->getDownward(elem, 0, vert);

  vd_ext_down_geo(m, vert, 4, geom_type, geom_tag);
  vd_ext_down_geo(m, edge, 6, geom_type + 4, geom_tag + 4);
  vd_ext_down_geo(m, surf, 4, geom_type + 10, geom_tag + 10);
  vd_ext_down_geo(m, &elem, 1, geom_type + 14, geom_tag + 14);

}


// Try gmi functionalities: TODO should use buildElement if ever to be used.

// Destroy an element and recreate its subentities with a certain geometric
// membership. 
// geom_tag and geom_type lists are geometric membership information for 
// new entities, which are integer arrays of size 4 + 6 + 4 + 1(vert, edge, 
// face, elem). These are ordered such that they correspond to the MDS 
// structure.

void vd_gmi_remap(apf::Mesh2* m, apf::MeshEntity* elem, int* geom_type
                                                      , int* geom_tag) {

  // LOOK: Actually, this whole thing can be done with down.
  // But this might be more general in certain cases. ?Which are?
  printf("Element %p is to be destroyed and recreated.\n", (void *)elem);

  apf::Downward surf;
  apf::Downward edge;
  apf::Downward vert;

  int size = m->getDownward(elem, 2, surf);
  size = m->getDownward(elem, 1, edge);
  size = m->getDownward(elem, 0, vert);

  std::vector<apf::Vector3> vert_pos(4, apf::Vector3(0, 0, 0));

  // Find the lower adjacencies, store their geometry type and ID.
  // Store the vertex points.
  // Starting from below (is it necessary?) destroy the mesh entities.

  // Starting from below, create new entities.
  // apf::Mesh:TYPE vertex: 0, edge: 1, triangle: 2, tetrahedron: 4

  printf("Copying the vertices... ");
  for (int i = 0; i < 4; i++) {
    apf::ModelEntity* em = m->toModel(vert[i]);
    geom_type[i] = m->getModelType(em);
    geom_tag[i] = m->getModelTag(em);
    m->getPoint(vert[i], 0, vert_pos.at(i));
    printf("Vertex %p, geom type %d tag %d, ", (void *)vert[i], geom_type[i], 
                                                        geom_tag[i]);

  }
  printf(". Vertices copied. Destroying vertices...\n");
  printf("New vertices: ");
  for (int i = 0; i < 4; i++) {
    m->destroy(vert[i]);
    apf::ModelEntity* em = m->findModelEntity(geom_type[i], geom_tag[i]);
    vert[i] = m->createVert(em);
    // apf::reorderMdsMesh(m);
    m->setPoint(vert[i], 0, vert_pos.at(i));
    printf("Vertex %p, geom type %d tag %d, ", (void *)vert[i], geom_type[i], 
                                                        geom_tag[i]);
  }
  printf(". Destroyed vertices and recreated new ones.\n");

  printf("Copying the edges... ");
  for (int i = 0; i < 6; i++) {
    apf::ModelEntity* em = m->toModel(edge[i]);
    geom_type[i+4] = m->getModelType(em);
    geom_tag[i+4] = m->getModelTag(em);
    printf("Edge %p, geom type %d tag %d, ", (void *)edge[i], geom_type[i+4], 
                                                      geom_tag[i+4]);
  }
  printf(". Edges copied. Destroying edges...\n");

  // The downward adjacencies are ordered in MDS fashion. After destroying 
  // old edges, newly created edges are to be associated with correct vertices.
  // The order is given in MDS documentation, and applied in topo_geom.h, as
  // vd_edge_vert and vd_face_vert.
  apf::MeshEntity* edge_vert[2];

  printf("New edges: ");
  for (int i = 0; i < 6; i++) {
    edge_vert[0] = (apf::MeshEntity*) vert[vd_edge_vert[i][0]];
    edge_vert[1] = (apf::MeshEntity*) vert[vd_edge_vert[i][1]];

    m->destroy(edge[i]);
    apf::ModelEntity* em = m->findModelEntity(geom_type[i+4], geom_tag[i+4]);
    edge[i] = makeOrFind(m, em, 1, edge_vert);
    printf("Edge %p, geom type %d tag %d, ", (void *)edge[i], geom_type[i+4], 
                                                      geom_tag[i+4]);
  }

  printf(". Destroyed edges and recreated new ones.\n");

  printf("Copying the faces... ");
  for (int i = 0; i < 4; i++) {
    apf::ModelEntity* em = m->toModel(surf[i]);
    geom_type[i+10] = m->getModelType(em);
    geom_tag[i+10] = m->getModelTag(em);
    printf("Face %p, geom type %d tag %d, ", (void *)surf[i], geom_type[i+10], 
                                                      geom_tag[i+10]);
  }
  printf(". Faces copied. Destroying faces...\n");

  apf::MeshEntity* surf_vert[3];

  printf("New faces: ");
  for (int i = 0; i < 4; i++) {
    surf_vert[0] = (apf::MeshEntity*) vert[vd_face_vert[i][0]];
    surf_vert[1] = (apf::MeshEntity*) vert[vd_face_vert[i][1]];
    surf_vert[2] = (apf::MeshEntity*) vert[vd_face_vert[i][2]];

    m->destroy(surf[i]);
    apf::ModelEntity* em = m->findModelEntity(geom_type[i+10], geom_tag[i+10]);
    surf[i] = makeOrFind(m, em, 2, surf_vert);
    printf("Face %p, geom type %d tag %d, ", (void *)surf[i], geom_type[i+10], 
                                                      geom_tag[i+10]);
  }

  printf("Destroyed faces and recreated new ones.\n");

  printf("Copying the element... ");
  apf::ModelEntity* em = m->toModel(elem);

  printf("Element copied. Destroying element...\n");

  printf("Element %p, geom type %d tag %d, ", (void *)elem, geom_type[14], 
                                                    geom_tag[14]);
  m->destroy(elem);
  em = m->findModelEntity(geom_type[14], geom_tag[14]);
  elem = makeOrFind(m, em, 3, vert);
  printf("New element: ");
  printf("Element %p, geom type %d tag %d, ", (void *)elem, geom_type[14], 
                                                    geom_tag[14]);
  printf("Destroyed element and recreated the new one.\n");

}

// Destroy and recreate an element and its sub entities:

void vd_gmi_recreate(apf::Mesh2* m, apf::MeshEntity* elem) {
  // LOOK: Actually, this whole thing can be done with down.
  // But this might be more general in certain cases. ?Which are?
  printf("Element %p is to be destroyed and recreated.\n", (void *)elem);

  apf::Downward surf;
  apf::Downward edge;
  apf::Downward vert;

  int size = m->getDownward(elem, 2, surf);
  size = m->getDownward(elem, 1, edge);
  size = m->getDownward(elem, 0, vert);

  std::vector<apf::Vector3> vert_pos(4, apf::Vector3(0, 0, 0));
  std::vector<int> vert_em_type(4);
  std::vector<int> vert_em_tag(4);

  // Find the lower adjacencies, store their geometry type and ID.
  // Store the vertex points.
  // Starting from below (is it necessary?) destroy the mesh entities.

  // Starting from below, create new entities.
  // apf::Mesh:TYPE vertex: 0, edge: 1, triangle: 2, tetrahedron: 4

  printf("Copying the vertices... ");
  for (int i = 0; i < 4; i++) {
    apf::ModelEntity* em = m->toModel(vert[i]);
    vert_em_type.at(i) = m->getModelType(em);
    vert_em_tag.at(i) = m->getModelTag(em);
    m->getPoint(vert[i], 0, vert_pos.at(i));
    printf("Vertex %p, ", (void *)vert[i]);
  }
  printf(". Vertices copied. Destroying vertices...\n");
  printf("New vertices: ");
  for (int i = 0; i < 4; i++) {
    m->destroy(vert[i]);
    apf::ModelEntity* em = m->findModelEntity(vert_em_type[i], vert_em_tag[i]);
    vert[i] = m->createVert(em);
    m->setPoint(vert[i], 0, vert_pos.at(i));
    printf("Vertex %p, ", (void *)vert[i]);
  }
  printf(". Destroyed vertices and recreated new ones.\n");

  int edge_em_type[6];
  int edge_em_tag[6];
  printf("Copying the edges... ");
  for (int i = 0; i < 6; i++) {
    apf::ModelEntity* em = m->toModel(edge[i]);
    edge_em_type[i] = m->getModelType(em);
    edge_em_tag[i] = m->getModelTag(em);
    printf("Edge %p, ", (void *)edge[i]);
  }
  printf(". Edges copied. Destroying edges...\n");

  // The downward adjacencies are ordered in MDS fashion. After destroying 
  // old edges, newly created edges are to be associated with correct vertices.
  // The order is given in MDS documentation, and applied in topo_geom.h, as
  // vd_edge_vert and vd_face_vert.
  apf::MeshEntity* edge_vert[2];

  printf("New edges: ");
  for (int i = 0; i < 6; i++) {
    edge_vert[0] = (apf::MeshEntity*) vert[vd_edge_vert[i][0]];
    edge_vert[1] = (apf::MeshEntity*) vert[vd_edge_vert[i][1]];

    m->destroy(edge[i]);
    apf::ModelEntity* em = m->findModelEntity(edge_em_type[i], edge_em_tag[i]);
    edge[i] = makeOrFind(m, em, 1, edge_vert);
    printf("Edge %p, ", (void *)edge[i]);
  }

  printf(". Destroyed edges and recreated new ones.\n");

  int surf_em_type[4];
  int surf_em_tag[4];
  printf("Copying the faces... ");
  for (int i = 0; i < 4; i++) {
    apf::ModelEntity* em = m->toModel(surf[i]);
    surf_em_type[i] = m->getModelType(em);
    surf_em_tag[i] = m->getModelTag(em);
  }
  printf("Faces copied. Destroying faces...\n");

  apf::MeshEntity* surf_vert[3];

  for (int i = 0; i < 4; i++) {
    surf_vert[0] = (apf::MeshEntity*) vert[vd_face_vert[i][0]];
    surf_vert[1] = (apf::MeshEntity*) vert[vd_face_vert[i][1]];
    surf_vert[2] = (apf::MeshEntity*) vert[vd_face_vert[i][2]];

    m->destroy(surf[i]);
    apf::ModelEntity* em = m->findModelEntity(surf_em_type[i], surf_em_tag[i]);
    surf[i] = makeOrFind(m, em, 1, surf_vert);
  }

  printf("Destroyed faces and recreated new ones.\n");

  printf("Copying the element... ");
  apf::ModelEntity* em = m->toModel(elem);
  int elem_em_type = m->getModelType(em);
  int elem_em_tag = m->getModelTag(em);

  printf("Element copied. Destroying element...\n");
  m->destroy(elem);
  em = m->findModelEntity(elem_em_type, elem_em_tag);
  elem = makeOrFind(m, em, 1, vert);
  printf("Destroyed element and recreated the new one.\n");

  // m->destroy(MeshEntity* e)
  /* MeshEntity* buildElement(
    Mesh2* m,
    ModelEntity* c,
    int type,
    MeshEntity** verts,
    BuildCallback* cb = 0);
  */
  /* ModelEntity* me = findModelEntity(vert_em_type[i], vert_em_tag[i]);
     MeshEntity* vert_temp = m->createVert(me);
     m->setPoint(vert_temp, 0, vert_pos[i]);
     MeshEntity* createVertex(ModelEntity* c, Vector3 const& point,
        Vector3 const& param);
  */

}

// Destroy and recreate an element and its sub entities:

void vd_gmi_recreate2(apf::Mesh2* m, apf::MeshEntity* elem) {

  // LOOK: Actually, this whole thing can be done with down.
  // But this might be more general in certain cases. ?Which are?

  apf::Downward surf;
  Entity_set es_surf = Entity_set();
  Entity_set es_edge = Entity_set();
  Entity_set es_vert = Entity_set();

  int size = m->getDownward(elem, 2, surf);
  copy_ent_set(&es_surf, surf, size);
  vd_set_down(m, &es_surf, &es_edge, 1);
  vd_set_down(m, &es_edge, &es_vert, 1);

  std::vector<apf::Vector3> vert_pos(es_vert.n, apf::Vector3(0, 0, 0));
  std::vector<int> vert_em_type(es_vert.n);
  std::vector<int> vert_em_tag(es_vert.n);

  // Find the lower adjacencies, store their geometry type and ID.
  // Store the vertex points.
  // Starting from above (is it necessary?) destroy the mesh entities.

  // Starting from below, create new entities.
  // apf::Mesh:TYPE vertex: 0, edge: 1, triangle: 2, tetrahedron: 4
  printf("Copying the vertices... ");
  for (int i = 0; i < es_vert.n; i++) {
    apf::ModelEntity* em = m->toModel(es_vert.e[i]);
    vert_em_type.at(i) = m->getModelType(em);
    vert_em_tag.at(i) = m->getModelTag(em);
    m->getPoint(es_vert.e[i], 0, vert_pos.at(i));
  }
  printf("Vertices copied. Destroying vertices...\n");

  for (int i = 0; i < es_vert.n; i++) {
    m->destroy(es_vert.e[i]);
    apf::ModelEntity* em = m->findModelEntity(vert_em_type[i], vert_em_tag[i]);
    es_vert.e[i] = m->createVert(em);
    m->setPoint(es_vert.e[i], 0, vert_pos.at(i));
  }
  printf("Destroyed vertices and recreated new ones.\n");

  std::vector<int> edge_em_type(es_edge.n);
  std::vector<int> edge_em_tag(es_edge.n);
  printf("Copying the edges... ");
  for (int i = 0; i < es_edge.n; i++) {
    apf::ModelEntity* em = m->toModel(es_edge.e[i]);
    edge_em_type.at(i) = m->getModelType(em);
    edge_em_tag.at(i) = m->getModelTag(em);
  }
  printf("Edges copied. Destroying edges...\n");

  for (int i = 0; i < es_edge.n; i++) {
    m->destroy(es_edge.e[i]);
    apf::ModelEntity* em = m->findModelEntity(edge_em_type[i], edge_em_tag[i]);

  }
  printf("Destroyed vertices and recreated new ones.\n");

  // m->destroy(MeshEntity* e)
  /* MeshEntity* buildElement(
    Mesh2* m,
    ModelEntity* c,
    int type,
    MeshEntity** verts,
    BuildCallback* cb = 0);
  */
  /* ModelEntity* me = findModelEntity(vert_em_type[i], vert_em_tag[i]);
     MeshEntity* vert_temp = m->createVert(me);
     m->setPoint(vert_temp, 0, vert_pos[i]);
     MeshEntity* createVertex(ModelEntity* c, Vector3 const& point,
        Vector3 const& param);
  */

}


