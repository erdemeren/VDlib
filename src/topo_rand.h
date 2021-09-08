#ifndef TOPO_RAND_H
#define TOPO_RAND_H

#include <apfMesh2.h>
#include <PCU.h>

#include <cstdlib>
#include <ctime>

void vd_init_rand();

void vd_rand_vect3(apf::Vector3& r);

float vd_norm_vect3(const apf::Vector3& r);

void vd_print_vect3(apf::Vector3& r);

#endif
