#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_extinfo.h"
#include "topo_field.h"

// This is a test function to attach a diagonally increasing scalar field. 
apf::Field* vd_att_test_field(apf::Mesh2* m) {

  apf::Field* temp_field = apf::createField(m, "temperature_field", 
                                              apf::SCALAR, m->getShape());

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) 
  {
    apf::Vector3 point;
    m->getPoint(e, 0, point);
    apf::setScalar(temp_field, e, 0, 0.1*point.getLength());
  }
  m->end(it);
  return temp_field;

}

// Attach the displacement field. 
apf::Field* vd_att_disp_field(apf::Mesh2* m) {

  apf::Field* disp_field;
  // Attach a vector field
  if (m->findField("displacement_field")) {
    disp_field = m->findField("displacement_field");
  }
  else {
    disp_field = apf::createField(m, "displacement_field", apf::VECTOR, 
                                                          m->getShape());
  }

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;

  apf::Vector3 disp;
  double abc[3] = {0,0,0}; disp.fromArray(abc);

  while ((e = m->iterate(it))) 
  {
    // This works with other mesh entities, where one can specify current 
    // node to be modified:
    apf::setVector(disp_field, e, 0, disp);
  }  
  m->end(it);
  return disp_field;
}

// Attach the velocity field. 
apf::Field* vd_att_vel_field(apf::Mesh2* m) {

  apf::Field* vel_field;
  // Attach a vector field
  if (m->findField("velocity_field")) {
    vel_field = m->findField("velocity_field");
  }
  else {
    vel_field = apf::createField(m, "velocity_field", apf::VECTOR, 
                                                      m->getShape());
  }

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;

  apf::Vector3 disp;
  double abc[3] = {0,0,0}; disp.fromArray(abc);

  while ((e = m->iterate(it))) 
  {
    // This works with other mesh entities, where one can specify current 
    // node to be modified:
    apf::setVector(vel_field, e, 0, disp);
  }  
  m->end(it);
  return vel_field;
}

// Attach the a scalar field over vertices. 
apf::Field* vd_att_vs_field(apf::Mesh2* m, const char* f_name) {

  apf::Field* field_curr;
  // Attach a vector field
  if (m->findField(f_name)) {
    field_curr = m->findField(f_name);
  }
  else {
    field_curr = apf::createField(m, f_name, apf::SCALAR, m->getShape());
  }

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;

  while ((e = m->iterate(it)))
  {
    // This works with other mesh entities, where one can specify current 
    // node to be modified:
    apf::setScalar(field_curr, e, 0, 0);
  }  
  m->end(it);
  return field_curr;
}

// Attach the displacement field. 
apf::Field* vd_att_vv_field(apf::Mesh2* m, const char* f_name) {

  apf::Field* field_curr;
  // Attach a vector field
  if (m->findField(f_name)) {
    field_curr = m->findField(f_name);
  }
  else {
    field_curr = apf::createField(m, f_name, apf::VECTOR, m->getShape());
  }

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;

  apf::Vector3 disp(0,0,0);

  while ((e = m->iterate(it)))
  {
    // This works with other mesh entities, where one can specify current 
    // node to be modified:
    apf::setVector(field_curr, e, 0, disp);
  }  
  m->end(it);
  return field_curr;
}

// Attach the displacement field. 
apf::Field* vd_att_vm_field(apf::Mesh2* m, const char* f_name) {

  apf::Field* field_curr;
  // Attach a vector field
  if (m->findField(f_name)) {
    field_curr = m->findField(f_name);
  }
  else {
    field_curr = apf::createField(m, f_name, apf::MATRIX, m->getShape());
  }

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;

  apf::Matrix3x3 zero(0,0,0,0,0,0,0,0,0);

  while ((e = m->iterate(it)))
  {
    // This works with other mesh entities, where one can specify current 
    // node to be modified:
    apf::setMatrix(field_curr, e, 0, zero);
  }  
  m->end(it);
  return field_curr;
}

// Attach a scalar field over triangles.
apf::Field* vd_att_es_field(apf::Mesh2* m, const char* f_name) {

  apf::Field* field_curr;
  // Attach a vector field
  if (m->findField(f_name)) {
    field_curr = m->findField(f_name);
  }
  else {
    field_curr = apf::createField(m, f_name, apf::SCALAR, apf::getConstant(1));
  }

  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;

  while ((e = m->iterate(it)))
  {
    // This works with other mesh entities, where one can specify current 
    // node to be modified:
    apf::setScalar(field_curr, e, 0, 0);
  }  
  m->end(it);
  return field_curr;
}

// Attach a scalar field over triangles.
apf::Field* vd_att_ts_field(apf::Mesh2* m, const char* f_name) {

  apf::Field* field_curr;
  // Attach a vector field
  if (m->findField(f_name)) {
    field_curr = m->findField(f_name);
  }
  else {
    field_curr = apf::createField(m, f_name, apf::SCALAR, apf::getConstant(2));
  }

  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;

  while ((e = m->iterate(it)))
  {
    // This works with other mesh entities, where one can specify current 
    // node to be modified:
    apf::setScalar(field_curr, e, 0, 0);
  }  
  m->end(it);
  return field_curr;
}

apf::Vector3 vd_get_vec(apf::Mesh2* m, apf::MeshEntity* v, const char* f_name) {
  apf::Field* vel_field;
  // Attach a vector field
  if (vel_field = m->findField(f_name)) {
    vel_field = m->findField(f_name);
    apf::Vector3 out(0,0,0);
    apf::getVector(vel_field, v, 0, out);
    return out;
  }
  else {
    return apf::Vector3(0,0,0);
  }
}

/*
// Attach the velocity field. 
apf::Field* vd_att_vel_field(apf::Mesh2* m) {

  apf::Field* imm_field;
  // Attach a vector field
  if (m->findField("immobile_field")) {
    imm_field = m->findField("immobile_field");
  }
  else {
    imm_field = apf::createField(m, "immobile_field", 1, apf::SCALAR);
  }

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;

  apf::Vector3 disp;
  double abc[3] = {0,0,0}; disp.fromArray(abc);

  while ((e = m->iterate(it))) 
  {
    // This works with other mesh entities, where one can specify current 
    // node to be modified:
    apf::setScalar(imm_field, e, 0, 1);
  }  
  m->end(it);
  return vel_field;
}
*/
