#include <cstdio>
#include <iostream>
#include <cstring>
#include <cassert>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include "topo_cg.h"

vd_relax_vel::vd_relax_vel(apf::Mesh2* m_in, cell_base* cb_in, 
                    field_calc* f_calc_in, double g_th_in, double e_th_in) : 
                  // Initialization of mesh related information 
                  mesh(m_in), c_base(cb_in),
                  f_calc(f_calc_in), 
                  e_list(m_in, cb_in),
                  g_th(g_th_in), e_th(e_th_in), 
                  f_field(NULL),
                  // Conjugate gradient initialization
                  v_set(0),
                  t_set(0),
                  tets(0),
                  x0(0, apf::Vector3(0,0,0)),
                  x1(0, apf::Vector3(0,0,0)),
                  g0(0, apf::Vector3(0,0,0)),
                  g1(0, apf::Vector3(0,0,0)),
                  dir(0, apf::Vector3(0,0,0)),
                  ndir(0, apf::Vector3(0,0,0)),
                  v_pos{},

                  iter_limit(400),
                  iter(0), g_norm(0), 
                  xa(0), xb(0), xc(0), xd(0), xe(0), xu(0),
                  fa(0), fb(0), fc(0), fd(0), fe(0), fu(0),
                  p(0), q(0), r(0), s(0), t(0), m(0), tol(0), tol2(0),
                  inv_quad_step(false),
                  Phi0(0), Phi1(0) {

  f_field = mesh->findField("velocity_field");
  bool skip_ext = (f_calc->get_proj() == (int)PROJ_TYPE::FIXED);
  /////////////////////////////////////////////
  // Collect boundary vertices:              //
  /////////////////////////////////////////////

  int e_nbr = 0;
  int t_nbr = 0;
  int dim = 0;
  for(dim = 0; dim < 3; dim++) {
    for(int i = 0; i < e_list.e.at(dim).size(); i++) {
      if(!c_base->is_free(dim, i) and (!skip_ext or 
                                      !c_base->get_cell_ext(dim, i)))
        e_nbr = e_nbr + e_list.e.at(dim).at(i).at(0).size();
    }
  }

  dim = 2;
  for(int i = 0; i < e_list.e.at(dim).size(); i++) {
    if(!c_base->is_free(dim, i) and (!skip_ext or 
                                    !c_base->get_cell_ext(dim, i)))
      t_nbr = t_nbr + e_list.e.at(dim).at(i).at(2).size();
  }

  x0.reserve(e_nbr);
  x1.reserve(e_nbr);
  g0.reserve(e_nbr);
  g1.reserve(e_nbr);
  dir.reserve(e_nbr);
  ndir.reserve(e_nbr);
  v_set.reserve(e_nbr);
  t_set.reserve(t_nbr);

  apf::Vector3 zero(0,0,0);
  apf::Vector3 temp(0,0,0);

  for(dim = 0; dim < 3; dim++) {
    for(int i = 0; i < e_list.e.at(dim).size(); i++) {
      if(!c_base->is_free(dim, i) and (!skip_ext or 
                                      !c_base->get_cell_ext(dim, i))) {
        for(int j = 0; j < e_list.e.at(dim).at(i).at(0).size(); j++) {
          apf::MeshEntity* v_curr = e_list.e.at(dim).at(i).at(0).at(j);
          v_set.push_back(v_curr);
          mesh->getPoint(v_curr, 0, temp);
          x0.push_back(temp);
          x1.push_back(zero);
          g0.push_back(zero);
          g1.push_back(zero);
          dir.push_back(zero);
          ndir.push_back(zero);
        }
      }
    }
  }

  dim = 2;
  for(int i = 0; i < e_list.e.at(dim).size(); i++) {
    if(!c_base->is_free(dim, i) and (!skip_ext or 
                                    !c_base->get_cell_ext(dim, i))) {
      for(int j = 0; j < e_list.e.at(dim).at(i).at(dim).size(); j++) {
        apf::MeshEntity* t_curr = e_list.e.at(dim).at(i).at(dim).at(j);
        t_set.push_back(t_curr);
      }
    }
  }

  std::vector<apf::MeshEntity*> ee(0);
  std::vector<apf::MeshEntity*> et(0);
  vd_set_up(mesh, &v_set, &ee);
  vd_set_up(mesh, &ee, &et);
  vd_set_up(mesh, &et, &tets);

  apf::writeVtkFiles("./output/Relax_before", mesh);
}

// Check inversion.
bool vd_relax_vel::chk_inv() {
  for (int i = 0; i < tets.size(); i++) {
    if(!vd_volume_tet_sign(mesh, tets.at(i)))
      return true;
  }
  return false;
}

// Save vertex positions.
void vd_relax_vel::save_pos() {
  for (int i = 0; i < v_set.size(); i++) {
    v_pos[v_set.at(i)] = x0.at(i);
  }
}

// Revert vertex positions.
void vd_relax_vel::rev_pos() {
  for (int i = 0; i < v_set.size(); i++) {
    mesh->setPoint(v_set.at(i), 0, v_pos[v_set.at(i)]);
  }
}

// Calculate the potential.
double vd_relax_vel::calc_energy() {
  double Phi = 0;
  for (int i = 0; i < t_set.size(); i++) {
    Phi = f_calc->en_tri(mesh, t_set.at(i));
  }
  assert(!std::isnan(Phi));
  return Phi;
}

// TODO this works well enough, but it doesn't consider the cross terms
// grad_i phi_l and only considers grad_i phi_i. See SM of the method paper.
void vd_relax_vel::upd_grad(std::vector<apf::Vector3> &g_curr) {
  //f_calc->vd_upd_vel_field(mesh, &v_set);
  for (int i = 0; i < v_set.size(); i++) {
    dir.at(i) = f_calc->vd_calc_force(mesh, v_set.at(i));
    //apf::getVector(f_field, v_set.at(i), 0, dir.at(i));
    g_curr.at(i) = dir.at(i)*(-1);
  }
}

void vd_relax_vel::shift_v_sp_pos(std::vector<apf::Vector3> &x_curr,
                    double x_in) {
  apf::Vector3 pos(0,0,0);
  //apf::Vector3 temp(0,0,0);
  for(int i = 0; i < v_set.size(); i++) {
    if (f_calc->chk_vert_special(mesh, v_set.at(i)))
      pos = f_calc->get_vec_special(mesh, v_set.at(i), ndir.at(i)*x_in);
    else
      pos = ndir.at(i)*x_in;

    assert(!std::isnan(pos.getLength()));
    pos = x_curr.at(i) + pos;
    assert(!std::isnan(pos.getLength()));
    mesh->setPoint(v_set.at(i), 0, pos);
  }
}

void vd_relax_vel::relax() {
  save_pos();
  Phi0 = calc_energy();
  upd_grad(g0);

  for (int i = 0; i < v_set.size(); i++) {
    dir.at(i) = g0.at(i)*(-1);
  }
  g_norm = 0;
  for (int i = 0; i < v_set.size(); i++) {
    g_norm = g_norm + dir.at(i) * dir.at(i);
  }
  g_norm = std::sqrt(g_norm);

  for (int i = 0; i < v_set.size(); i++) {
    ndir.at(i) = dir.at(i) / g_norm;
  }

  for (int b = 0; b < cg_MAX_CG_ITER; ++b) {
    ////////////////////////////////////////////////////////////
    // bounds the line search, with ax < bx < cx and fa > fb < fc
    ////////////////////////////////////////////////////////////
    xa = 0.;
    fa = Phi0;

    xb = cg_SQRT_EPS;
    shift_v_sp_pos(x0, xb);
    fb = calc_energy();

    while (fb > fa && xb > cg_EPS) {
      // decrease step until energy is decreasing along ndir
      xb /= 10.;
      shift_v_sp_pos(x0, xb);
      fb = calc_energy();
    }

    // try inverse quadratic interpolation
    p = 0;
    for (int i = 0; i < v_set.size(); i++) {
      p = p - g0.at(i)*ndir.at(i);
    }
    p = p * xb;
    q = (fb - fa) + p;
    inv_quad_step = false;
    if (q > cg_EPS) {
      // parabola is concave up, find the minimum
      xc = (p * xb) / (2. * q);
      if (xc > (cg_MAX_BND_STEP + 1.) * xb) {
        // maximum step length
        inv_quad_step = true;
        // cg_MAX_BND_STEP could be too large when xb is not reduced.
        //if (xb/STEP_FIRST < 0.01)
        //  xc = (100 + 1.) * xb;
        //else
        //  xc = (cg_MAX_BND_STEP + 1.) * xb;
        xc = (cg_MAX_BND_STEP + 1.) * xb;

        shift_v_sp_pos(x0, xc);
        fc = calc_energy();
      } else if (xc > (cg_BND_MAG + 1.) * xb) {
        // normal step
        inv_quad_step = true;
        shift_v_sp_pos(x0, xc);
        fc = calc_energy();
        if (fc < fb) {
          // try to step past minimum
          cg_SHIFT2(xa, xb, xc);
          xc = xb + cg_BND_MAG * (xb - xa);
          cg_SHIFT2(fa, fb, fc);
          shift_v_sp_pos(x0, xc);
          fc = calc_energy();
        }
      } else if (xc > xa + cg_SQRT_EPS && xc < xb - cg_SQRT_EPS) {
        // minimum falls in (ax, bx)
        shift_v_sp_pos(x0, xc);
        fc = calc_energy();
        if (fc < fb) {
          // found bracket, all done
          inv_quad_step = true; 
          std::swap(xb, xc);
          std::swap(fb, fc);
        }
      }
    }
    if (!inv_quad_step) {
      // quadratic interpolation failed, conservative step
      xc = (cg_BND_MAG + 1.) * xb;
      shift_v_sp_pos(x0, xc);
      fc = calc_energy();
    }
    while (fc < fb) {
      // try inverse quadratic interpolation
      p = xc - xb;
      q = xa - xb;
      r = (fa - fb) * p;
      s = (fc - fb) * q;
      t = r - s;
      inv_quad_step = false;
      if (t > cg_EPS) {
        // parabola is concave up, find the minimum
        xd = xb + (r * p - s * q) / (2. * t);
        if (xd > xc + cg_MAX_BND_STEP * p) {
          // maximum step length
          inv_quad_step = true;
          xd = xc + cg_MAX_BND_STEP * p;

          shift_v_sp_pos(x0, xd);
          fd = calc_energy();
        } else if (xd > xc + cg_BND_MAG * p) {
          // normal step
          inv_quad_step = true;
          shift_v_sp_pos(x0, xd);
          fd = calc_energy();
          if (fd < fc) {
            // try to step past minimum
            cg_SHIFT3(xa, xb, xc, xd);
            xd = xc + cg_BND_MAG * (xc - xb);
            cg_SHIFT3(fa, fb, fc, fd);
            shift_v_sp_pos(x0, xd);
            fd = calc_energy();
          }
        } else if (xd > xb + cg_SQRT_EPS && xd < xc - cg_SQRT_EPS) {
          // minimum falls in (ax, bx)
          shift_v_sp_pos(x0, xd);
          fd = calc_energy();
          if (fd < fc) {
            // found bracket, all done
            inv_quad_step = true;
            cg_SHIFT2(xa, xb, xd);
            cg_SHIFT2(fa, fb, fd);
            break;
          }
        }
      }
      if (!inv_quad_step) {
        // quadratic interpolation failed, conservative step
        xd = xc + cg_BND_MAG * p;
        shift_v_sp_pos(x0, xd);
        fd = calc_energy();
      }
      // bookkeeping for next iteration
      cg_SHIFT3(xa, xb, xc, xd);
      cg_SHIFT3(fa, fb, fc, fd);
    }

    ////////////////////////////////////////////////////////////
    // Brent's method to find minimum along search direction, a translation
    // of the ALGOL 60 algorithm on page 79 of R. P. Brent, Algorithms for
    // Minimization Without Derivatives, 1973 with minor modifications. The
    // author gave permission to use this algorithm in private communication.
    ////////////////////////////////////////////////////////////
    // use values from bounding the line search
    if (fc < fa) {
      xd = xc;
      xe = xa;
      fd = fc;
      fe = fa;
    } else {
      xd = xa;
      xe = xc;
      fd = fa;
      fe = fc;
    }
    t = (xb < 0.5 * (xa + xc) ? xc : xa) - xb;
    s = cg_PHI_SQ_INV * t;
    int c = 0;
    for (c = 0; c < cg_MAX_LS_ITER; ++c) {
      m = 0.5 * (xa + xc);
      tol  = cg_SQRT_EPS * (fabs(xb) + 1.);
      tol2 = 2. * tol;
      // check stopping criterion
      if (fabs(xb - m) > tol2 - 0.5 * (xc - xa)) {
        inv_quad_step = false;
        if (fabs(t) > tol) {
          // inverse quadratic interpolation
          p = (xb - xd) * (fb - fe);
          q = (xb - xe) * (fb - fd);
          r = (xb - xe) * q - (xb - xd) * p;
          q = 2. * (q - p);
          if (q > 0.) { r = -r; }
          q = fabs(q);
          cg_SHIFT2(p, t, s);
          // mistake in ALGOL 60 routine, second condition is inverted
          if (fabs(r) < fabs(0.5 * q * p) && r > q * (xa - xb) && r < q * (xc - xb)) {
            // take inverse quadratic interpolation step
            inv_quad_step = true;
            s = r / q;
            xu = xb + s;
            // f should not be evaluated too close to xa or xc
            if (xu - xa < tol2 || xc - xu < tol2) {
              s = (xb < m ? tol : -tol);
            }
          }
        }
        if (!inv_quad_step) {
          // interpolation failed, take golden section step
          t = (xb < m ? xc : xa) - xb;
          s = cg_PHI_SQ_INV * t;
        }

        // f should not be evaluated too close to xb
        xu = xb + (fabs(s) >= tol ? s : (s > 0. ? tol : -tol));

        shift_v_sp_pos(x0, xu);
        fu = calc_energy();
        // bookkeeping for next iteration
        if (fu <= fb) {
          if (xu < xb) { xc = xb; } else { xa = xb; }
          cg_SHIFT3(xe, xd, xb, xu);
          cg_SHIFT3(fe, fd, fb, fu);
        } else {
          if (xu < xb) { xa = xu; } else { xc = xu; }
          if (fu <= fd || xd == xb) {
            cg_SHIFT2(xe, xd, xu);
            cg_SHIFT2(fe, fd, fu);
          } else if (fu <= fe || xe == xb || xe == xd) {
            xe = xu;
            fe = fu;
          }
        }
      }
      else {
        for (int i = 0; i < v_set.size(); i++) {
          x1.at(i) = x0.at(i) + ndir.at(i) * xb;
        }

        // found minimum, apply change and update energy
        Phi1 = fb;
        break;
      }
    }

    if (c == cg_MAX_LS_ITER)
      std::cerr << "RunSimulation: max LS iteration reached" << std::endl;

    ////////////////////////////////////////////////////////////
    // Conjugate gradient
    ////////////////////////////////////////////////////////////
    // check energy convergence
    if (fabs(Phi1 - Phi0) < e_th) {
      x0 = x1;
      break;
    }

    g_norm = 0;
    for (int i = 0; i < v_set.size(); i++) {
      g_norm = g_norm + dir.at(i) * dir.at(i);
    }
    g_norm = std::sqrt(g_norm);

    if(std::fabs(g_norm) < std::numeric_limits<double>::min())
      g_norm = 1;

    // check gradient convergence
    if (g_norm < g_th) {
      x0 = x1;
      break;
    }
    ////////////////////////////////////////////////////////////
    // Conjugate gradient
    ////////////////////////////////////////////////////////////
    // check energy convergence
    //if (fabs(Phi1 - Phi0) < E_TOL) {
    //  x0 = x1;
    //}

    // check gradient convergence
    shift_v_sp_pos(x1, 0);

    upd_grad(g1);

    // Direction update given by Y.H Dai, C.X. Kou, SIAM Journal of
    // Optimization, v 23, p 296-320, 2013
    for (int i = 0; i < v_set.size(); i++) {
      g0.at(i) = g0.at(i) - g1.at(i);
    }
    for (int i = 0; i < v_set.size(); i++) {
      x0.at(i) = x0.at(i) - x1.at(i);
    }

    double B0 = 0;
    double B1 = 0;
    double B0_a = 0;
    double B0_b = 0;
    double B0_c = 0;

    double B1_a = 0;
    double B1_b = 0;
    for (int i = 0; i < v_set.size(); i++) {
      B0_a = B0_a + g0.at(i)*g0.at(i);
      B0_b = B0_a + g1.at(i)*x0.at(i);
      B0_c = B0_a + g0.at(i)*x0.at(i);
    }
    if(std::fabs(B0_c) < std::numeric_limits<double>::min())
      B0_c = 1;
    B0 = B0_a*B0_b/B0_c;

    for (int i = 0; i < v_set.size(); i++) {
      B0_a = B0_a + g0.at(i)*g1.at(i);
      B0_b = B0_a + dir.at(i)*g0.at(i);

      B1_a = B0_a + g1.at(i)*dir.at(i);
      B1_b = B0_a + dir.at(i)*dir.at(i);
    }
    if(std::fabs(B0_b) < std::numeric_limits<double>::min())
      B0_b = 1;

    B0 = (B0_a - B0)/B0_b;

    if(std::fabs(B1_b) < std::numeric_limits<double>::min())
      B1_b = 1;

    B1 = B1_a / (2. * B1_b);

    double mult = fmax(B0, B1);
    for (int i = 0; i < v_set.size(); i++) {
      dir.at(i) = dir.at(i) * mult - g1.at(i);
    }

    g_norm = 0;
    for (int i = 0; i < v_set.size(); i++) {
      g_norm = g_norm + dir.at(i) * dir.at(i);
    }
    g_norm = std::sqrt(g_norm);

    if(std::fabs(g_norm) < std::numeric_limits<double>::min())
      g_norm = 1;

    for (int i = 0; i < v_set.size(); i++) {
      ndir.at(i) = dir.at(i) / g_norm;
    }

    // prepare for next iteration
    for (int i = 0; i < v_set.size(); i++) {
      x0.at(i) = x1.at(i);
    }
    shift_v_sp_pos(x0, 0);

    Phi0 = Phi1;
    for (int i = 0; i < v_set.size(); i++) {
      g0.at(i) = g1.at(i);
    }
  }
}

vd_relax_vel::~vd_relax_vel() {
}

