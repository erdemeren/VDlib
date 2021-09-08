#include <cstdio>
#include <iostream>
#include <cstring>
#include <cassert>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include "topo_energy.h"

// The trial energy functions.
double energy_func_tri_trial(apf::Mesh2* m, apf::MeshEntity* tri) {

  int ent_type = m->getType(tri);
  int d = m->typeDimension[ent_type];
  int m_type = m->getModelType(m->toModel(tri));
  //std::cout << "Entity being measured: " << d << ", " << ent_type << ", "
  //          << m->getModelType(m->toModel(tri)) << ", "
  //          << m->getModelTag(m->toModel(tri)) << std::endl;

  apf::MeshElement* ee = createMeshElement(m, tri);
  double meas = measure(ee);
  destroyMeshElement(ee);

  //std::cout << "Area is " << meas << std::endl;

  if (d == 2 and m_type == 2) {
    return meas;
  }
  //std::cout << "Returning 0" << std::endl;

  return 0;
}

double energy_func_tet_trial(apf::Mesh2* m, apf::MeshEntity* tet) {
  //int ent_type = m->getType(tet);
  //int d = m->typeDimension[ent_type];

  //if (d == 3 and m->getModelType(m->toModel(tet)) == 3) {
  //  apf::MeshElement* ee = createMeshElement(m, tet);
  //  return measure(ee);
  //}
  return 0;
}

vd_energy::vd_energy() {
  //m = m_in;
  energy_func_tri = energy_func_tri_trial;
  energy_func_tet = energy_func_tet_trial;
}
/*
void vd_energy::update_mesh(apf::Mesh2* m_in) {
  m = m_in;
}
*/
double vd_energy::en_tet(apf::Mesh2* m, apf::MeshEntity* tet) {
  return (*energy_func_tet)(m, tet);
}

double vd_energy::en_tri(apf::Mesh2* m, apf::MeshEntity* tri) {
  return (*energy_func_tri)(m, tri);
}

