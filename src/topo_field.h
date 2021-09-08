#ifndef TOPO_FIELD_H
#define TOPO_FIELD_H

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

// This is a test function to attach a diagonally increasing scalar field. 
apf::Field* vd_att_test_field(apf::Mesh2* m);

// Attach the displacement field. 
apf::Field* vd_att_disp_field(apf::Mesh2* m);

// Attach the velocity field. 
apf::Field* vd_att_vel_field(apf::Mesh2* m);

// Attach the scalar, vectorial and tensor vertex field. 
apf::Field* vd_att_vs_field(apf::Mesh2* m, const char* f_name);
apf::Field* vd_att_vv_field(apf::Mesh2* m, const char* f_name);
apf::Field* vd_att_vm_field(apf::Mesh2* m, const char* f_name);

// Attach a scalar field over triangles.
apf::Field* vd_att_ts_field(apf::Mesh2* m, const char* f_name);

// Given a vector field and a vertex, return the field value at the vertex.
apf::Vector3 vd_get_vec(apf::Mesh2* m, apf::MeshEntity* v, const char* f_name);

// Attach immobile field. Used mainly with the artificial exterior grain  
// boundary, to preserve volume of inner grains.
//void vd_att_imm_field(apf::Mesh2* m, int c3_ext);

// Update immobile field.
//void vd_upd_imm_field(apf::Mesh2* m, int c3_ext);

#endif
