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

#include "topo_pca.h"

void ellipsoid::clear() {
  double z[3] = {0,0,0};
  ctr.fromArray(z);
  dummy_clear_stop();

  e_v.clear();
  e_v.resize(3);
  for(int i = 0; i < 3; i++) {
    e_v.at(i).fromArray(z);
  }
  e_i.clear();
  e_i.resize(3);
  for(int i = 0; i < 3; i++) {
    e_i.at(i) = 0;
  }
}

ellipsoid::ellipsoid () : ctr(0, 0, 0), e_v(0, apf::Vector3(0, 0, 0)), e_i(0), 
                                                                      dim(0) {

  clear();
}

ellipsoid::ellipsoid(const ellipsoid &obj) :
  ctr(0, 0, 0), e_v(0, apf::Vector3(0, 0, 0)), e_i(0), dim(0) {

  assert(obj.e_v.size() == 3 and obj.e_i.size() == 3);
  clear();
  ctr = obj.ctr;
  for(int i = 0; i < 3; i++) {
    e_v.at(i) = obj.e_v.at(i);
    e_i.at(i) = obj.e_i.at(i);
  }
}

ellipsoid& ellipsoid::operator=( const ellipsoid& obj ) {
  assert(obj.e_v.size() == 3 and obj.e_i.size() == 3);
  clear();
  ctr = obj.ctr;
  for(int i = 0; i < 3; i++) {
    e_v.at(i) = obj.e_v.at(i);
    e_i.at(i) = obj.e_i.at(i);
  }
  dim = obj.dim;
  return *this;
}


std::ostream& operator<<(std::ostream& s, 
                     const ellipsoid& e) {
  s << "Ellipsoid ctr: " << e.ctr << std::endl;
  s << "\tMajor axes: " << std::endl;
  for(int i = 0; i < 3; i++) {
    s << "\t\t" << e.e_v.at(i) << " " << e.e_i.at(i) << std::endl;
  }

  return s;
}

// https://gist.github.com/microo8/4065693
// Assume the mean is subtracted.
void pca(const gsl_matrix* m_gsl, unsigned int L, std::vector<std::vector<double> >* e_v, std::vector<double>* e_i) {

    assert(m_gsl != NULL);
    assert(L > 0 && L < m_gsl->size2+1);
    /*
    @param data - matrix of data vectors, MxN matrix, each column is a data vector, M - dimension, N - data vector count
    @param L - dimension reduction
    */
    unsigned int i;
    unsigned int rows = m_gsl->size1;
    unsigned int cols = m_gsl->size2;
/*
    gsl_vector* mean = gsl_vector_alloc(rows);

    for(i = 0; i < rows; i++) {
        gsl_vector_set(mean, i, gsl_stats_mean(m_gsl->data + i * cols, 1, cols));
    }

    // Get mean-substracted data into matrix mean_substracted_data.
    gsl_matrix* mean_substracted_data = gsl_matrix_alloc(rows, cols);
    gsl_matrix_memcpy(mean_substracted_data, m_gsl);
    for(i = 0; i < cols; i++) {
        gsl_vector_view mean_substracted_point_view = gsl_matrix_column(mean_substracted_data, i);
        gsl_vector_sub(&mean_substracted_point_view.vector, mean);
    }
    gsl_vector_free(mean);

    // Compute Covariance matrix
    gsl_matrix* covariance_matrix = gsl_matrix_alloc(rows, rows);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / (double)(cols - 1), mean_substracted_data, mean_substracted_data, 0.0, covariance_matrix);
    gsl_matrix_free(mean_substracted_data);
*/
    gsl_matrix* covariance_matrix = gsl_matrix_alloc(rows, rows);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / (double)(cols - 1), m_gsl, m_gsl, 0.0, covariance_matrix);

    // Get eigenvectors, sort by eigenvalue.
    gsl_vector* eigenvalues = gsl_vector_alloc(rows);
    gsl_matrix* eigenvectors = gsl_matrix_alloc(rows, rows);
    gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(rows);
    gsl_eigen_symmv(covariance_matrix, eigenvalues, eigenvectors, workspace);
    gsl_eigen_symmv_free(workspace);
    gsl_matrix_free(covariance_matrix);

    // Sort the eigenvectors
    gsl_eigen_symmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);

    // Project the original dataset
    gsl_matrix* result = gsl_matrix_alloc(L, cols);
    gsl_matrix_view L_eigenvectors = gsl_matrix_submatrix(eigenvectors, 0, 0, rows, L);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &L_eigenvectors.matrix, m_gsl, 0.0, result);

    for(unsigned int ii = 0; ii < e_v->size(); ii++)
      e_v->at(ii).clear();

    e_v->clear();
    e_i->clear();

    e_v->resize(eigenvectors->size2);
    for(unsigned int ii = 0; ii < eigenvectors->size2; ii++) {
      e_v->at(ii).resize(eigenvectors->size1);
      for(unsigned int j = 0; j < eigenvectors->size1; j++)
        e_v->at(ii).at(j) = gsl_matrix_get (eigenvectors, j, ii);
    }

    e_i->resize(3);
    for(unsigned int ii = 0; ii < e_i->size(); ii++) {
      e_i->at(ii) = gsl_vector_get (eigenvalues, ii);
    }

    gsl_matrix_free(eigenvectors);
    gsl_vector_free(eigenvalues);
    gsl_matrix_free(result);

    // Result is n LxN matrix, each column is the original data vector with reduced dimension from M to L
}

// Assume the mean is subtracted.
// https://gist.github.com/microo8/4065693
void pca(const std::vector<std::vector<double> >* m_in, unsigned int L, std::vector<std::vector<double> >* e_v, std::vector<double>* e_i) {

  assert(m_in->size() > 0);
  assert(m_in->at(0).size() > 0);

  unsigned int M = m_in->size();
  unsigned int N = m_in->at(0).size();

  gsl_matrix* m_gsl = gsl_matrix_alloc (N, M);

  for (int i = 0; i < (int)(M); i++) {
    for (int j = 0; j < (int)(N); j++) {
      gsl_matrix_set (m_gsl, j, i, m_in->at(i).at(j));
    }
  }

  pca(m_gsl, L, e_v, e_i);
  gsl_matrix_free(m_gsl);
}

// https://gist.github.com/microo8/4065693
ellipsoid pca(const std::vector<apf::Vector3>* pos_in, unsigned int L) {

  assert(pos_in->size() > 0);
  assert(L > 0 and L < 4);

  ellipsoid ellip;

  unsigned int M = pos_in->size();
  unsigned int N = 3;

  gsl_matrix* m_gsl = gsl_matrix_alloc (N, M);

  for (int i = 0; i < (int)(M); i++) {
    for (int j = 0; j < (int)(N); j++) {
      gsl_matrix_set (m_gsl, j, i, pos_in->at(i)[j]);
    }
  }

  std::vector<std::vector<double> > e_vv(0, std::vector<double>(0));
  std::vector<double> e_i(0);

  pca(m_gsl, L, &e_vv, &e_i);
  gsl_matrix_free(m_gsl);

  for (int i = 0; i < e_vv.size(); i++) {
    for (int j = 0; j < e_vv.at(i).size(); j++) {
      ellip.e_v.at(i)[j] = e_vv.at(i).at(j);
    }
    ellip.e_i.at(i) = e_i.at(i);
  }

  for (int i = 0; i < e_vv.size(); i++) {
    e_vv.at(i).clear();
  }
  e_vv.clear();

  ellip.dim = L;

  return ellip;
}
