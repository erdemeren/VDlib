#ifndef TOPO_PCA_H
#define TOPO_PCA_H

#include <assert.h>

#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include "topo_extinfo.h"
// Container for the ellipsoid associated with a set of entities.
class ellipsoid {
  private:
  public:
    apf::Vector3 ctr;
    std::vector<apf::Vector3> e_v;
    std::vector<double> e_i;
    int dim;

    void clear();

    ellipsoid();
    ellipsoid(const ellipsoid &obj);
    ellipsoid& operator=( const ellipsoid& obj );
    friend std::ostream& operator<<(std::ostream& s, const ellipsoid& e);

};

void pca(const gsl_matrix* m_gsl, unsigned int L, std::vector<std::vector<double> >* e_v, std::vector<double>* e_i);

void pca(const std::vector<std::vector<double> >* m_in, unsigned int L, std::vector<std::vector<double> >* e_v, std::vector<double>* e_i);

ellipsoid pca(const std::vector<apf::Vector3>* pos_in, unsigned int L);

#endif
