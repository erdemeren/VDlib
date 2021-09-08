#include <apfMesh2.h>
#include <PCU.h>

#include <cstdlib>
#include <ctime>

void vd_init_rand() {
  srand (static_cast <unsigned> (time(0)));
}

float vd_norm_vect3(const apf::Vector3& r) {
  float r_size = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
  return r_size;
}

void vd_rand_vect3(apf::Vector3& r) {
  r[0] = 2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)-0.5);
  r[1] = 2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)-0.5);
  r[2] = 2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)-0.5);
  float r_size = vd_norm_vect3(r);
  r[0] = r[0]/r_size;
  r[1] = r[1]/r_size;
  r[2] = r[2]/r_size;
}

void vd_print_vect3(apf::Vector3& r) {
  printf("Some random vector in the interval (-1.0, 1.0)^3: (%1.2f, %1.2f, %1.2f).\n", r[0], r[1], r[2]);  
}

