#ifndef TOPO_GMI_H
#define TOPO_GMI_H

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_rand.h"


// Used in gmi_remap: Given an element, extract the geometric membership of 
// entities. See vd_gmi_remap.

void vd_ext_down_geo(apf::Mesh2* m, apf::MeshEntity** down, int size
                                                        , int* geom_type
                                                        , int* geom_tag );
void vd_ext_elem_geo(apf::Mesh2* m, apf::MeshEntity* elem, int* geom_type
                                                         , int* geom_tag);

// Try gmi functionalities: Destroy and recreate an element and its 
// sub entities:
// LOOK 
void vd_gmi_remap(apf::Mesh2* m, apf::MeshEntity* elem, int* geom_type
                                                      , int* geom_tag);

void vd_gmi_recreate(apf::Mesh2* m, apf::MeshEntity* elem);

void vd_gmi_recreate2(apf::Mesh2* m, apf::MeshEntity* elem);


#endif
