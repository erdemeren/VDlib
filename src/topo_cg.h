#ifndef TOPO_CG_H
#define TOPO_CG_H

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_topo.h"
#include "topo_entlist.h"
#include "topo_disp.h"

/*

// The trial energy functions.
double cg_grad_trial();
double cg_energy_trial();
// The conjugate gradient implementation of Jeremy with customizable energy and  
// gradient functions and the definition of position.
// Energy and gradient should evaluate positions and parameters and 
class vd_cg{
  private:
    //apf::Mesh2* m;
    virtual void calc_grad() = 0;
    virtual void calc_energy() = 0;
  public:
    vd_cg();
    ~vd_cg() {};

};
*/

// Adaptation of Jeremy's implementation of conjugate gradient with quadratic
// approximation of the function in the direction of the gradient for faster
// convergence. 
#define cg_MAX_CG_ITER 250
#define cg_MAX_LS_ITER 100

#define cg_E_TOL_DEF 1.4901e-08
#define cg_G_TOL_DEF 1.4901e-08

#define cg_MAX_BND_STEP 5
#define cg_BND_MAG 1.618033988749895
#define cg_PHI_SQ_INV 0.381966011250105

#define cg_SQRT_EPS 1.4901e-08
#define cg_EPS 2.2204e-16

#define cg_SHIFT2(a, b, c) do { (a) = (b); (b) = (c); } while (0)
#define cg_SHIFT3(a, b, c, d) do { (a) = (b); (b) = (c); (c) = (d); } while (0)

// Container for conjugate gradient relaxation related variables.
class vd_relax_vel{
  private:
  public:
    apf::Mesh2* mesh;
    cell_base* c_base;
    field_calc* f_calc;

    vd_entlist e_list;

    apf::Field* f_field;

    // The set of moving vertices:
    std::vector<apf::MeshEntity*> v_set;
    std::vector<apf::MeshEntity*> t_set;
    std::vector<apf::MeshEntity*> tets;
    // Conjugate gradient related:
    std::vector<apf::Vector3> x0;
    std::vector<apf::Vector3> x1;

    std::vector<apf::Vector3> g0;
    std::vector<apf::Vector3> g1;
    std::vector<apf::Vector3> dir;
    std::vector<apf::Vector3> ndir;

    std::map<apf::MeshEntity*, apf::Vector3> v_pos;

    bool inv_flag;

    int iter_limit;
    int iter;
    double g_norm;
    double xa, xb, xc, xd, xe, xu;
    double fa, fb, fc, fd, fe, fu;
    double p, q, r, s, t, m, tol, tol2;
    bool inv_quad_step;
    double Phi0, Phi1;

    double e_th, g_th;

    // Constructor:
    vd_relax_vel(apf::Mesh2* m_in, cell_base* cb_in, field_calc* f_calc_in, double g_th_in = cg_G_TOL_DEF, double e_th_in = cg_E_TOL_DEF);

    bool chk_inv();
    // Save and revert positions
    void save_pos();
    void rev_pos();

    // Calculate the potential.
    double calc_energy();
    void upd_grad(std::vector<apf::Vector3> &g_curr);

    void shift_v_sp_pos(std::vector<apf::Vector3> &x_curr, double x_in);
    // Try to find a hull for non-inverting tets by conjugate gradient method.
    void relax();

    ~vd_relax_vel();
};


#endif
