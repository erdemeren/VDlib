#ifndef TOPO_SL_H
#define TOPO_SL_H

#include <assert.h>

#include <vector>
#include <stdio.h>
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

#include "topo_pca.h"
#include "topo_disp.h"

// https://www.gnu.org/software/gsl/doc/html/splinalg.html

// Container for the gsl based linear solver:
class solv_lin {
  private:
    solv_lin(const solv_lin &obj);
    solv_lin& operator=( const solv_lin& obj );
  public:
    int n;
    apf::Mesh2* m;

    gsl_spmatrix* A;
    gsl_spmatrix* C;
    gsl_vector *f;
    gsl_vector *u;

    apf::GlobalNumbering* gn;

    // Functions:
    void numb_vert();
    void dest_numb();
    void collect_matrix(apf::Mesh2* m_in);
    void init(apf::Mesh2* m_in);
    void clear();
    double solve(const double tol);

    solv_lin(apf::Mesh2* m);
    ~solv_lin();
    friend std::ostream& operator<<(std::ostream& s, const solv_lin& sl);

};

// Container for the gsl based linear solver:
class solv_kuprat {
  private:
    solv_kuprat(const solv_kuprat &obj);
    solv_kuprat& operator=( const solv_kuprat& obj );
  public:
    apf::Mesh2* m;
    field_calc* f_calc;

    int n;
    std::map<apf::MeshEntity*, int> v2id;
    std::map<apf::MeshEntity*, bool> v_map;
    std::map<int, apf::MeshEntity*> id2v;

    std::vector<apf::MeshEntity*>* verts;

    gsl_spmatrix* A;
    gsl_spmatrix* C;
    gsl_vector *f;
    gsl_vector *u;

    apf::GlobalNumbering* gn;

    // Functions:
    void numb_vert();
    void numb_verts_in();
    void dest_numb();
    void add_forces(apf::MeshEntity* tri);

    void collect_matrix_verts();
    void collect_matrix();
    void init(apf::Mesh2* m_in);
    void init_verts(apf::Mesh2* m_in);
    void clear();
    double solve(const double tol);

    solv_kuprat(apf::Mesh2* m, field_calc* f_calc_in);
    solv_kuprat(apf::Mesh2* m, field_calc* f_calc_in, std::vector<apf::MeshEntity*>* verts);

    ~solv_kuprat();

    friend std::ostream& operator<<(std::ostream& s, const solv_kuprat& sl);

};

// Container for the gsl based linear solver:
class solv_kuprat_NBC {
  private:
    solv_kuprat_NBC(const solv_kuprat_NBC &obj);
    solv_kuprat_NBC& operator=( const solv_kuprat_NBC& obj );
  public:
    apf::Mesh2* m;
    field_calc* f_calc;

    int n;
    std::map<apf::MeshEntity*, int> v2id;
    std::map<apf::MeshEntity*, bool> v_map;
    std::map<int, apf::MeshEntity*> id2v;

    std::map<apf::MeshEntity*, bool> v_sp;
    std::vector<apf::MeshEntity*>* verts;

    gsl_spmatrix* A;
    gsl_spmatrix* C;
    gsl_vector *f;
    gsl_vector *u;

    apf::GlobalNumbering* gn;

    // Functions:
    void numb_vert();
    void numb_verts_in();
    void dest_numb();
    void add_forces(apf::MeshEntity* tri);

    void collect_matrix_verts();
    void collect_matrix();

		void process_ext_shell();
		void process_ext_shell_verts();

    void init(apf::Mesh2* m_in);
    void init_verts(apf::Mesh2* m_in);
    void clear();
    double solve(const double tol);

    solv_kuprat_NBC(apf::Mesh2* m, field_calc* f_calc_in);
    solv_kuprat_NBC(apf::Mesh2* m, field_calc* f_calc_in, std::vector<apf::MeshEntity*>* verts);

    ~solv_kuprat_NBC();

    friend std::ostream& operator<<(std::ostream& s, const solv_kuprat_NBC& sl);

};
#endif
