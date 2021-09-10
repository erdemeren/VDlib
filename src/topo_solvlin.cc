#include <assert.h>

#include <vector>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include "topo_solvlin.h"

// solv_kuprat

void solv_kuprat::numb_vert() {
  PCU_Barrier();

  int max_gn = 0;
  apf::MeshIterator *vertIter = m->begin(0);
  apf::MeshEntity *vert = nullptr;
  while ((vert = m->iterate(vertIter))) {
    //if (m->isOwned(vert) && m->getOwner(vert) == p_id && m->getModelType(m->toModel(vert)) != 3) {
    if (m->getModelType(m->toModel(vert)) != 3) {
      v2id[vert] = max_gn;
      v_map[vert] = true;
      id2v[max_gn] = vert;
      max_gn = max_gn + 1;
    }
  }
  m->end(vertIter);
  n = PCU_Max_Int(max_gn);
}

void solv_kuprat::numb_verts_in() {
  assert(verts != NULL);

  PCU_Barrier();

  int max_gn = 0;

  apf::MeshEntity *vert = nullptr;
  for (int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    if (m->getModelType(m->toModel(vert)) != 3) {
      v2id[vert] = max_gn;
      v_map[vert] = true;
      id2v[max_gn] = vert;
      max_gn = max_gn + 1;
    }
  }
  n = PCU_Max_Int(max_gn);
}

void solv_kuprat::dest_numb() {
}

void solv_kuprat::add_forces(apf::MeshEntity* tri) {
  apf::ModelEntity* em = m->toModel(tri);
  int em_type = m->getModelType(em);
  if(em_type == 2) {
    apf::Downward d_v;
    apf::Vector3 p_i(0,0,0);
    apf::Vector3 p_j(0,0,0);
    apf::Vector3 n_ij(0,0,0);
    apf::Vector3 force(0,0,0);

    std::vector<apf::Vector3> p(3, apf::Vector3(0,0,0));
    m->getDownward(tri, 0, d_v);

    m->getPoint(d_v[0], 0, p.at(0));
    m->getPoint(d_v[1], 0, p.at(1));
    m->getPoint(d_v[2], 0, p.at(2));

    assert(f_calc->d2(m, tri) > std::numeric_limits<double>::min());

    double f_coef = f_calc->gam2(m, tri)/(f_calc->d2(m, tri)*2);
    for(int i = 0; i < 3; i++) {
      if(v_map[d_v[i]]) {
        int n_curr = v2id[d_v[i]];
        p_i = p.at((i+1) % 3) - p.at(i);
        p_j = p.at((i+2) % 3) - p.at(i);
        n_ij = vd_cross(p_i, p_j);
        n_ij = norm_0(n_ij);
        // Not the actual force, but the quantity corresponding to the reduced 
        // mobility.
        force = vd_cross(n_ij, p_i-p_j)*f_coef;
        for(int ii = 0; ii < 3; ii++) {
          double temp = gsl_vector_get(f, 3*n_curr + ii);
          temp = temp + force[ii];
          assert(!std::isnan(temp));
          gsl_vector_set(f, 3*n_curr + ii, temp);
        }
      }
    }
  }
}

void solv_kuprat::collect_matrix_verts() {
  A = gsl_spmatrix_alloc(3*n, 3*n); /* triplet format */

  f = gsl_vector_alloc(3*n);        /* right hand side vector */
  u = gsl_vector_alloc(3*n);        /* solution vector */
  size_t i;
  gsl_vector_set_zero(f);

  /* construct the sparse matrix for the finite difference equation */

  apf::Vector3 p(0,0,0);
  apf::Vector3 a_norm(0,0,0);
  apf::Vector3 a_norm_proj(0,0,0);
  p = apf::Vector3(0,0,0);
  apf::Matrix3x3 J(0,0,0,0,0,0,0,0,0);
  apf::Matrix3x3 m_ij(0,0,0,0,0,0,0,0,0);

  double detJ = 0;
  apf::MeshIterator* it_e = m->begin(2);
  apf::MeshEntity* tri;

  std::vector<apf::MeshEntity*> edges(0);
  std::vector<apf::MeshEntity*> tris(0);
  vd_set_up(m, verts, &edges);
  vd_set_up(m, &edges, &tris);

  for (int i = 0; i < tris.size(); i++) {
    tri = tris.at(i);

    apf::ModelEntity* mdl = m->toModel(tri);
    bool isValid = true;
    if(m->isOwned(tri) && m->getType(tri) && apf::Mesh::TRIANGLE && m->getModelType(mdl) == 2) {

      apf::MeshElement* me = apf::createMeshElement(m,tri);
      apf::getJacobian(me,p,J);

      a_norm = norm_0(vd_area_out_n(m, tri));
      m_ij = tensorProduct(a_norm, a_norm);

      if(m->getDimension() == 3){
        detJ = apf::getJacobianDeterminant(J,2);
      }
      apf::Downward d_v;
      m->getDownward(tri, 0, d_v);

      int n[3];
      for(int i = 0; i < 3; i++) {
        if(v_map[d_v[i]])
          n[i] = v2id[d_v[i]];
        else
          n[i] = -1;
      }

      // Add the elements of the area normal outer product multiplied by 
      // the determinant of the Jacobian to the stiffness matrix.
      // A is 3N x 3N.
      for(int i = 0; i < 3; i++) {
        if(n[i] != -1) {
          for(int ii = 0; ii < 3; ii++) {
            for(int jj = 0; jj < 3; jj++) {
              double temp = gsl_spmatrix_get(A, 3*n[i]+ii, 3*n[i]+jj);
              //
              temp = temp + m_ij[ii][jj]*detJ/12;
              assert(!std::isnan(temp));
              gsl_spmatrix_set(A, 3*n[i]+ii, 3*n[i]+jj, temp);
            }
          }
        }
      }

      for(int i = 0; i < 2; i++) {
        for(int j = i + 1; j < 3; j++) {
          if(n[i] != -1 and n[j] != -1) {
            for(int ii = 0; ii < 3; ii++) {
              for(int jj = 0; jj < 3; jj++) {
                double temp = gsl_spmatrix_get(A, 3*n[i]+ii, 3*n[j]+jj);
                temp = temp + m_ij[ii][jj]*detJ/24;
                assert(!std::isnan(temp));
                gsl_spmatrix_set(A, 3*n[i]+ii, 3*n[j]+jj, temp);
                gsl_spmatrix_set(A, 3*n[j]+ii, 3*n[i]+jj, temp);
              }
            }
          }
        }
      }
      add_forces(tri);
    }
  }
  m->end(it_e);

  /* convert to compressed column format */
  C = gsl_spmatrix_ccs(A);
}

void solv_kuprat::collect_matrix() {
  A = gsl_spmatrix_alloc(3*n, 3*n); /* triplet format */

  f = gsl_vector_alloc(3*n);        /* right hand side vector */
  u = gsl_vector_alloc(3*n);        /* solution vector */

  gsl_vector_set_zero(f);
  size_t i;

  /* construct the sparse matrix for the finite difference equation */

  apf::Vector3 p(0,0,0);
  apf::Vector3 a_norm(0,0,0);
  p = apf::Vector3(0,0,0);
  apf::Matrix3x3 J(0,0,0,0,0,0,0,0,0);
  apf::Matrix3x3 m_ij(0,0,0,0,0,0,0,0,0);

  double detJ = 0;
  apf::MeshIterator* it_e = m->begin(2);
  apf::MeshEntity* tri;

  while(tri = m->iterate(it_e)) {
    apf::ModelEntity* mdl = m->toModel(tri);
    bool isValid = true;
    if(m->isOwned(tri) && m->getType(tri) && apf::Mesh::TRIANGLE && m->getModelType(mdl) == 2) {
      apf::MeshElement* me = apf::createMeshElement(m,tri);
      apf::getJacobian(me,p,J);

      a_norm = norm_0(vd_area_out_n(m, tri));
      m_ij = tensorProduct(a_norm, a_norm);

      if(m->getDimension() == 3){
        detJ = apf::getJacobianDeterminant(J,2);
      }
      apf::Downward d_v;
      m->getDownward(tri, 0, d_v);

      int n[3];
      for(int i = 0; i < 3; i++)
        n[i] = v2id[d_v[i]];

      // Add the elements of the area normal outer product multiplied by 
      // the determinant of the Jacobian to the stiffness matrix.
      // A is 3N x 3N.
      for(int i = 0; i < 3; i++) {
        for(int ii = 0; ii < 3; ii++) {
          for(int jj = 0; jj < 3; jj++) {
            double temp = gsl_spmatrix_get(A, 3*n[i]+ii, 3*n[i]+jj);
            temp = temp + m_ij[ii][jj]*detJ/12;
            assert(!std::isnan(temp));

            gsl_spmatrix_set(A, 3*n[i]+ii, 3*n[i]+jj, temp);
          }
        }
      }

      for(int i = 0; i < 2; i++) {
        for(int j = i + 1; j < 3; j++) {
          for(int ii = 0; ii < 3; ii++) {
            for(int jj = 0; jj < 3; jj++) {
              double temp = gsl_spmatrix_get(A, 3*n[i]+ii, 3*n[j]+jj);
              temp = temp + m_ij[ii][jj]*detJ/24;
              assert(!std::isnan(temp));
              gsl_spmatrix_set(A, 3*n[i]+ii, 3*n[j]+jj, temp);
              gsl_spmatrix_set(A, 3*n[j]+ii, 3*n[i]+jj, temp);
            }
          }
        }
      }
      add_forces(tri);
    }
  }
  m->end(it_e);

  /* convert to compressed column format */
  C = gsl_spmatrix_ccs(A);
}

void solv_kuprat::init(apf::Mesh2* m_in) {
  m = m_in;
  numb_vert();

  collect_matrix();
}

void solv_kuprat::init_verts(apf::Mesh2* m_in) {
  m = m_in;
  numb_verts_in();
  collect_matrix_verts();
}


void solv_kuprat::clear() {
  dest_numb();
  id2v.clear();
  v2id.clear();
  v_map.clear();

  gsl_spmatrix_free(A);
  gsl_spmatrix_free(C);
  gsl_vector_free(f);
  gsl_vector_free(u);
}

double solv_kuprat::solve(const double tol = 1.0e-6) {

  const size_t max_iter = 10; /* maximum iterations */
  const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
  gsl_splinalg_itersolve *work =
    gsl_splinalg_itersolve_alloc(T, 3*n, 0);
  size_t iter = 0;
  double residual;
  int status = GSL_CONTINUE;

  /* initial guess u = 0 */
  gsl_vector_set_zero(u);

  /* solve the system A u = f */
  while (status == GSL_CONTINUE && ++iter < max_iter) {
    status = gsl_splinalg_itersolve_iterate(C, f, tol, u, work);

    /* print out residual norm ||A*u - f|| */
    residual = gsl_splinalg_itersolve_normr(work);
    //fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);

    if (status == GSL_SUCCESS)
      fprintf(stderr, "Converged\n");
  }
  std::cout << "Nbr of iterations: " << iter << std::endl;

  size_t i;
  /* output solution */
  apf::Field* vel_field = m->findField("velocity_field");
  apf::Field* for_field = m->findField("force_field");

  apf::Vector3 vel(0,0,0);
  for (int i = 0; i < n; ++i) {
    apf::MeshEntity* vert = id2v[i];
    for (int ii = 0; ii < 3; ++ii) {
      vel[ii] = gsl_vector_get(u, 3*i+ii);
    }
    apf::setVector(vel_field, vert, 0, vel);
  }
  gsl_splinalg_itersolve_free(work);
  return residual;
}

solv_kuprat::solv_kuprat(apf::Mesh2* m_in, field_calc* f_calc_in) : verts(NULL), v_map{}, v2id{}, id2v{}  {
  f_calc = f_calc_in;
  init(m_in);
}

solv_kuprat::solv_kuprat(apf::Mesh2* m_in, field_calc* f_calc_in, std::vector<apf::MeshEntity*>* verts_in) : verts(NULL), v_map{}, v2id{}, id2v{}  {
  f_calc = f_calc_in;
  verts = verts_in;
  init_verts(m_in);
}

solv_kuprat::~solv_kuprat() {
  clear();
}

std::ostream& operator<<(std::ostream& s, 
                     const solv_kuprat& sl) {
  s << "Linear solver trial: " << std::endl;
  return s;
}

// solv_kuprat_NBC

void solv_kuprat_NBC::numb_vert() {
  PCU_Barrier();

  int max_gn = 0;
  apf::MeshIterator *vertIter = m->begin(0);
  apf::MeshEntity *vert = nullptr;
  while ((vert = m->iterate(vertIter))) {
    //if (m->isOwned(vert) && m->getOwner(vert) == p_id && m->getModelType(m->toModel(vert)) != 3) {
    if (m->getModelType(m->toModel(vert)) != 3) {
      v2id[vert] = max_gn;
      v_map[vert] = true;
      id2v[max_gn] = vert;
      max_gn = max_gn + 1;
    }
    if(f_calc->chk_vert_special(m, vert)) {
      v_sp[vert] = true;
    }
  }
  m->end(vertIter);
  n = PCU_Max_Int(max_gn);
}

void solv_kuprat_NBC::numb_verts_in() {
  assert(verts != NULL);

  PCU_Barrier();

  int max_gn = 0;

  apf::MeshEntity *vert = nullptr;
  for (int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    if (m->getModelType(m->toModel(vert)) != 3) {
      v2id[vert] = max_gn;
      v_map[vert] = true;
      id2v[max_gn] = vert;
      max_gn = max_gn + 1;
    }
    if(f_calc->chk_vert_special(m, vert)) {
      v_sp[vert] = true;
    }
  }
  n = PCU_Max_Int(max_gn);
}

void solv_kuprat_NBC::dest_numb() {
}

void solv_kuprat_NBC::add_forces(apf::MeshEntity* tri) {
  apf::ModelEntity* em = m->toModel(tri);
  int em_type = m->getModelType(em);
  if(em_type == 2) {
    apf::Downward d_v;
    apf::Vector3 p_i(0,0,0);
    apf::Vector3 p_j(0,0,0);
    apf::Vector3 n_ij(0,0,0);
    apf::Vector3 force(0,0,0);

    std::vector<apf::Vector3> p(3, apf::Vector3(0,0,0));
    m->getDownward(tri, 0, d_v);

    m->getPoint(d_v[0], 0, p.at(0));
    m->getPoint(d_v[1], 0, p.at(1));
    m->getPoint(d_v[2], 0, p.at(2));

    assert(f_calc->d2(m, tri) > std::numeric_limits<double>::min());
    double f_coef = f_calc->gam2(m, tri)/(f_calc->d2(m, tri)*2);
    for(int i = 0; i < 3; i++) {
      if(v_map[d_v[i]]) {
        int n_curr = v2id[d_v[i]];
        p_i = p.at((i+1) % 3) - p.at(i);
        p_j = p.at((i+2) % 3) - p.at(i);
        n_ij = vd_cross(p_i, p_j);
        n_ij = norm_0(n_ij);
        // Not the actual force, but the quantity corresponding to the reduced 
        // mobility.
        force = vd_cross(n_ij, p_i-p_j)*f_coef;
        if(f_calc->chk_vert_special(m, d_v[i])) {
          force = f_calc->get_vec_special(m, d_v[i], force);
        }
        for(int ii = 0; ii < 3; ii++) {
          double temp = gsl_vector_get(f, 3*n_curr + ii);
          temp = temp + force[ii];
          assert(!std::isnan(temp));
          gsl_vector_set(f, 3*n_curr + ii, temp);
        }
      }
    }
  }
}

void solv_kuprat_NBC::collect_matrix_verts() {
  A = gsl_spmatrix_alloc(3*n, 3*n); /* triplet format */

  f = gsl_vector_alloc(3*n);        /* right hand side vector */
  u = gsl_vector_alloc(3*n);        /* solution vector */
  size_t i;
  gsl_vector_set_zero(f);

  /* construct the sparse matrix for the finite difference equation */

  apf::Vector3 p(0,0,0);
  apf::Vector3 a_norm(0,0,0);
  apf::Vector3 a_norm_proj(0,0,0);
  p = apf::Vector3(0,0,0);
  apf::Matrix3x3 J(0,0,0,0,0,0,0,0,0);
  apf::Matrix3x3 m_ij(0,0,0,0,0,0,0,0,0);
  apf::Matrix3x3 m_ij_proj(0,0,0,0,0,0,0,0,0);

  double detJ = 0;
  apf::MeshIterator* it_e = m->begin(2);
  apf::MeshEntity* tri;

  std::vector<apf::MeshEntity*> edges(0);
  std::vector<apf::MeshEntity*> tris(0);
  vd_set_up(m, verts, &edges);
  vd_set_up(m, &edges, &tris);

  for (int i = 0; i < tris.size(); i++) {
    tri = tris.at(i);

    apf::ModelEntity* mdl = m->toModel(tri);
    bool isValid = true;
    if(m->isOwned(tri) && m->getType(tri) && apf::Mesh::TRIANGLE && m->getModelType(mdl) == 2) {

      apf::MeshElement* me = apf::createMeshElement(m,tri);
      apf::getJacobian(me,p,J);

      a_norm = norm_0(vd_area_out_n(m, tri));
      m_ij = tensorProduct(a_norm, a_norm);

      if(m->getDimension() == 3){
        detJ = apf::getJacobianDeterminant(J,2);
      }
      apf::Downward d_v;
      m->getDownward(tri, 0, d_v);

      int n[3];
      for(int i = 0; i < 3; i++) {
        if(v_map[d_v[i]])
          n[i] = v2id[d_v[i]];
        else
          n[i] = -1;
      }

      // Add the elements of the area normal outer product multiplied by 
      // the determinant of the Jacobian to the stiffness matrix.
      // A is 3N x 3N.
      for(int i = 0; i < 3; i++) {
        if(n[i] != -1) {
          m_ij_proj = m_ij;
          if(v_sp[d_v[i]]) {
            a_norm_proj = f_calc->get_vec_special(m, d_v[i], a_norm);
            m_ij_proj = tensorProduct(a_norm_proj, a_norm_proj);
          }

          for(int ii = 0; ii < 3; ii++) {
            for(int jj = 0; jj < 3; jj++) {
              double temp = gsl_spmatrix_get(A, 3*n[i]+ii, 3*n[i]+jj);
              //
              //temp = temp + m_ij[ii][jj]*detJ/12;
              temp = temp + m_ij_proj[ii][jj]*detJ/12;
              assert(!std::isnan(temp));
              gsl_spmatrix_set(A, 3*n[i]+ii, 3*n[i]+jj, temp);
            }
          }
        }
      }

      for(int i = 0; i < 2; i++) {
        for(int j = i + 1; j < 3; j++) {
          if(n[i] != -1 and n[j] != -1) {
            for(int ii = 0; ii < 3; ii++) {
              for(int jj = 0; jj < 3; jj++) {
                double temp = gsl_spmatrix_get(A, 3*n[i]+ii, 3*n[j]+jj);
                //temp = temp + m_ij[ii][jj]*detJ/24;
                temp = temp + m_ij_proj[ii][jj]*detJ/24;
                assert(!std::isnan(temp));
                gsl_spmatrix_set(A, 3*n[i]+ii, 3*n[j]+jj, temp);
                gsl_spmatrix_set(A, 3*n[j]+ii, 3*n[i]+jj, temp);
              }
            }
          }
        }
      }
      add_forces(tri);
    }
  }
  m->end(it_e);

  /* convert to compressed column format */
  C = gsl_spmatrix_ccs(A);
}

void solv_kuprat_NBC::collect_matrix() {
  A = gsl_spmatrix_alloc(3*n, 3*n); /* triplet format */

  f = gsl_vector_alloc(3*n);        /* right hand side vector */
  u = gsl_vector_alloc(3*n);        /* solution vector */

  gsl_vector_set_zero(f);
  size_t i;

  /* construct the sparse matrix for the finite difference equation */

  apf::Vector3 p(0,0,0);
  apf::Vector3 a_norm(0,0,0);
  p = apf::Vector3(0,0,0);
  apf::Matrix3x3 J(0,0,0,0,0,0,0,0,0);
  apf::Matrix3x3 m_ij(0,0,0,0,0,0,0,0,0);
  apf::Vector3 a_norm_proj(0,0,0);
  apf::Matrix3x3 m_ij_proj(0,0,0,0,0,0,0,0,0);

  double detJ = 0;
  apf::MeshIterator* it_e = m->begin(2);
  apf::MeshEntity* tri;

  while(tri = m->iterate(it_e)) {
    apf::ModelEntity* mdl = m->toModel(tri);
    bool isValid = true;
    if(m->isOwned(tri) && m->getType(tri) && apf::Mesh::TRIANGLE && m->getModelType(mdl) == 2) {
      apf::MeshElement* me = apf::createMeshElement(m,tri);
      apf::getJacobian(me,p,J);

      a_norm = norm_0(vd_area_out_n(m, tri));
      m_ij = tensorProduct(a_norm, a_norm);

      if(m->getDimension() == 3){
        detJ = apf::getJacobianDeterminant(J,2);
      }
      apf::Downward d_v;
      m->getDownward(tri, 0, d_v);

      int n[3];
      for(int i = 0; i < 3; i++)
        n[i] = v2id[d_v[i]];

      // Add the elements of the area normal outer product multiplied by 
      // the determinant of the Jacobian to the stiffness matrix.
      // A is 3N x 3N.
      for(int i = 0; i < 3; i++) {

        m_ij_proj = m_ij;
        if(v_sp[d_v[i]]) {
          a_norm_proj = f_calc->get_vec_special(m, d_v[i], a_norm);
          m_ij_proj = tensorProduct(a_norm_proj, a_norm_proj);
        }
        for(int ii = 0; ii < 3; ii++) {
          for(int jj = 0; jj < 3; jj++) {
            double temp = gsl_spmatrix_get(A, 3*n[i]+ii, 3*n[i]+jj);
            // 
            temp = temp + m_ij_proj[ii][jj]*detJ/12;
            assert(!std::isnan(temp));

            gsl_spmatrix_set(A, 3*n[i]+ii, 3*n[i]+jj, temp);
          }
        }
      }

      for(int i = 0; i < 2; i++) {
        for(int j = i + 1; j < 3; j++) {
          for(int ii = 0; ii < 3; ii++) {
            for(int jj = 0; jj < 3; jj++) {
              double temp = gsl_spmatrix_get(A, 3*n[i]+ii, 3*n[j]+jj);
              temp = temp + m_ij_proj[ii][jj]*detJ/24;
              //temp = temp + m_ij[ii][jj]*detJ/24;
              assert(!std::isnan(temp));
              gsl_spmatrix_set(A, 3*n[i]+ii, 3*n[j]+jj, temp);
              gsl_spmatrix_set(A, 3*n[j]+ii, 3*n[i]+jj, temp);
            }
          }
        }
      }
      add_forces(tri);
    }
  }
  m->end(it_e);

  /* convert to compressed column format */
  C = gsl_spmatrix_ccs(A);
}

// TODO there is a discrepency. Try representing the magnitudes of the drag and
// force terms for the simulation cell vertices according to the Neumann BC.
void solv_kuprat_NBC::process_ext_shell() {
  ext_shell* e_sh = f_calc->get_e_sh();

  std::vector<apf::MeshEntity*> v_adj(0);
  std::vector<apf::MeshEntity*> es_edge(0);

  apf::MeshIterator *vertIter = m->begin(0);
  apf::MeshEntity *vert = nullptr;
  while ((vert = m->iterate(vertIter))) {
    //if (m->isOwned(vert) && m->getOwner(vert) == p_id && m->getModelType(m->toModel(vert)) != 3) {
    int c_dim = m->getModelType(m->toModel(vert));
    int c_tag = m->getModelTag(m->toModel(vert))-1;
    if (e_sh->chk_shell(c_dim, c_tag)) {
      shell sh_curr = e_sh->get_shell(c_dim, c_tag);
      int n_i = v2id[vert];
      if(sh_curr.dim == 1) {
        // Multiply the integrals over the domains around the vertices at the 
        // boundaries.
        for(int ii = 0; ii < 3; ii++) {
          for(int jj = 0; jj < 3; jj++) {
            double temp = 4*gsl_spmatrix_get(A, 3*n_i+ii, 3*n_i+jj);
            gsl_spmatrix_set(A, 3*n_i+ii, 3*n_i+jj, temp);
          }
        }

        vd_set_up(m, vert, &es_edge);
        vd_set_down(m, &es_edge, &v_adj);
        for(int i = 0; i < v_adj.size(); i++) {
          if(v_adj.at(i) != vert and v_map[v_adj.at(i)]) {
            int n_adj = v2id[v_adj.at(i)];
            for(int ii = 0; ii < 3; ii++) {
              for(int jj = 0; jj < 3; jj++) {
                double temp = gsl_spmatrix_get(A, 3*n_i+ii, 3*n_adj+jj);
                temp = temp*4;
                gsl_spmatrix_set(A, 3*n_i+ii, 3*n_adj+jj, temp);
              }
            }
          }
        }
        for(int ii = 0; ii < 3; ii++) {
          double temp = 4*gsl_vector_get(f, 3*n_i + ii);
          assert(!std::isnan(temp));
          gsl_vector_set(f, 3*n_i + ii, temp);
        }
      }
      else if(sh_curr.dim == 2) {
        for(int ii = 0; ii < 3; ii++) {
          for(int jj = 0; jj < 3; jj++) {
            int n_i = v2id[vert];
            double temp = 2*gsl_spmatrix_get(A, 3*n_i+ii, 3*n_i+jj);
            gsl_spmatrix_set(A, 3*n_i+ii, 3*n_i+jj, temp);
          }
        }
        for(int i = 0; i < v_adj.size(); i++) {
          if(v_adj.at(i) != vert and v_map[v_adj.at(i)]) {
            int n_adj = v2id[v_adj.at(i)];
            for(int ii = 0; ii < 3; ii++) {
              for(int jj = 0; jj < 3; jj++) {
                double temp = gsl_spmatrix_get(A, 3*n_i+ii, 3*n_adj+jj);
                temp = temp*2;
                gsl_spmatrix_set(A, 3*n_i+ii, 3*n_adj+jj, temp);
              }
            }
          }
        }
        for(int ii = 0; ii < 3; ii++) {
          double temp = 2*gsl_vector_get(f, 3*n_i + ii);
          assert(!std::isnan(temp));
          gsl_vector_set(f, 3*n_i + ii, temp);
        }
      }
    }
  }
  m->end(vertIter);
}

void solv_kuprat_NBC::process_ext_shell_verts() {
    //if (m->isOwned(vert) && m->getOwner(vert) == p_id && m->getModelType(m->toModel(vert)) != 3) {
  std::vector<apf::MeshEntity*> v_adj(0);
  std::vector<apf::MeshEntity*> es_edge(0);

  ext_shell* e_sh = f_calc->get_e_sh();
  for (int i = 0; i < verts->size(); i++) {
    apf::MeshEntity* vert = verts->at(i);

    int c_dim = m->getModelType(m->toModel(vert));
    int c_tag = m->getModelTag(m->toModel(vert))-1;
    if (e_sh->chk_shell(c_dim, c_tag)) {
      int n_i = v2id[vert];
      shell sh_curr = e_sh->get_shell(c_dim, c_tag);

      vd_set_up(m, vert, &es_edge);
      vd_set_down(m, &es_edge, &v_adj);

      if(sh_curr.dim == 1) {
        int n_i = v2id[vert];
        for(int ii = 0; ii < 3; ii++) {
          for(int jj = 0; jj < 3; jj++) {
            double temp = 4*gsl_spmatrix_get(A, 3*n_i+ii, 3*n_i+jj);
            gsl_spmatrix_set(A, 3*n_i+ii, 3*n_i+jj, temp);
          }
        }
        for(int j = 0; j < v_adj.size(); j++) {
          if(v_adj.at(j) != vert and v_map[v_adj.at(j)]) {
            int n_adj = v2id[v_adj.at(j)];
            for(int ii = 0; ii < 3; ii++) {
              for(int jj = 0; jj < 3; jj++) {
                double temp = gsl_spmatrix_get(A, 3*n_i+ii, 3*n_adj+jj);
                temp = temp*4;
                gsl_spmatrix_set(A, 3*n_i+ii, 3*n_adj+jj, temp);
              }
            }
          }
        }

        for(int ii = 0; ii < 3; ii++) {
          double temp = 4*gsl_vector_get(f, 3*n_i + ii);
          assert(!std::isnan(temp));
          gsl_vector_set(f, 3*n_i + ii, temp);
        }
      }
      else if(sh_curr.dim == 2) {
        for(int ii = 0; ii < 3; ii++) {
          for(int jj = 0; jj < 3; jj++) {
            int n_i = v2id[vert];
            double temp = 2*gsl_spmatrix_get(A, 3*n_i+ii, 3*n_i+jj);
            gsl_spmatrix_set(A, 3*n_i+ii, 3*n_i+jj, temp);
          }
        }
        for(int j = 0; j < v_adj.size(); j++) {
          if(v_adj.at(j) != vert and v_map[v_adj.at(j)]) {
            int n_adj = v2id[v_adj.at(j)];
            for(int ii = 0; ii < 3; ii++) {
              for(int jj = 0; jj < 3; jj++) {
                double temp = gsl_spmatrix_get(A, 3*n_i+ii, 3*n_adj+jj);
                temp = temp*2;
                gsl_spmatrix_set(A, 3*n_i+ii, 3*n_adj+jj, temp);
              }
            }
          }
        }
        for(int ii = 0; ii < 3; ii++) {
          double temp = 2*gsl_vector_get(f, 3*n_i + ii);
          assert(!std::isnan(temp));
          gsl_vector_set(f, 3*n_i + ii, temp);
        }
      }
    }
  }
}

void solv_kuprat_NBC::init(apf::Mesh2* m_in) {
  m = m_in;
  numb_vert();

  collect_matrix();

  if(f_calc->get_proj() == (int) PROJ_TYPE::EXT_SHELL)
    process_ext_shell();
}

void solv_kuprat_NBC::init_verts(apf::Mesh2* m_in) {
  m = m_in;
  numb_verts_in();
  collect_matrix_verts();
  if(f_calc->get_proj() == (int) PROJ_TYPE::EXT_SHELL)
    process_ext_shell_verts();
}


void solv_kuprat_NBC::clear() {
  dest_numb();
  id2v.clear();
  v2id.clear();
  v_map.clear();
  v_sp.clear();

  gsl_spmatrix_free(A);
  gsl_spmatrix_free(C);
  gsl_vector_free(f);
  gsl_vector_free(u);
}

double solv_kuprat_NBC::solve(const double tol = 1.0e-6) {

  const size_t max_iter = 10; /* maximum iterations */
  const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
  gsl_splinalg_itersolve *work =
    gsl_splinalg_itersolve_alloc(T, 3*n, 0);
  size_t iter = 0;
  double residual;
  int status = GSL_CONTINUE;

  /* initial guess u = 0 */
  gsl_vector_set_zero(u);

  /* solve the system A u = f */
  while (status == GSL_CONTINUE && ++iter < max_iter) {
    status = gsl_splinalg_itersolve_iterate(C, f, tol, u, work);

    /* print out residual norm ||A*u - f|| */
    residual = gsl_splinalg_itersolve_normr(work);
    //fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);

    if (status == GSL_SUCCESS)
      fprintf(stderr, "Converged\n");
  }
  std::cout << "Nbr of iterations: " << iter << std::endl;
  size_t i;
  /* output solution */
  apf::Field* vel_field = m->findField("velocity_field");
  apf::Field* for_field = m->findField("force_field");

  apf::Vector3 vel(0,0,0);
  for (int i = 0; i < n; ++i) {
    apf::MeshEntity* vert = id2v[i];
    for (int ii = 0; ii < 3; ++ii) {
      vel[ii] = gsl_vector_get(u, 3*i+ii);
    }
    apf::setVector(vel_field, vert, 0, vel);
  }
  gsl_splinalg_itersolve_free(work);
  return residual;
}

solv_kuprat_NBC::solv_kuprat_NBC(apf::Mesh2* m_in, field_calc* f_calc_in) : verts(NULL), v_map{}, v2id{}, id2v{}, v_sp{}  {
  f_calc = f_calc_in;
  init(m_in);
}

solv_kuprat_NBC::solv_kuprat_NBC(apf::Mesh2* m_in, field_calc* f_calc_in, std::vector<apf::MeshEntity*>* verts_in) : verts(NULL), v_map{}, v2id{}, id2v{}, v_sp{}  {
  f_calc = f_calc_in;
  verts = verts_in;
  init_verts(m_in);
}

solv_kuprat_NBC::~solv_kuprat_NBC() {
  clear();
}

std::ostream& operator<<(std::ostream& s, 
                     const solv_kuprat_NBC& sl) {
  s << "Linear solver trial: " << std::endl;
  return s;
}
