#ifndef TOPO_TESS_H
#define TOPO_TESS_H

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

void vd_create_tess(apf::Mesh2* m, const char* modelFile);

#endif
