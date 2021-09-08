#ifndef TOPO_ENERGY_H
#define TOPO_ENERGY_H

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

// The trial energy functions.
double energy_func_tri_trial(apf::Mesh2* m, apf::MeshEntity* tri);
double energy_func_tet_trial(apf::Mesh2* m, apf::MeshEntity* tet);

class vd_energy{
  private:
    //apf::Mesh2* m;
    double (*energy_func_tri)(apf::Mesh2*, apf::MeshEntity*);
    double (*energy_func_tet)(apf::Mesh2*, apf::MeshEntity*);

  public:
    vd_energy();
    ~vd_energy() {};
    //void update_mesh(apf::Mesh2* m_in);

    double en_tet(apf::Mesh2* m, apf::MeshEntity* tet);
    double en_tri(apf::Mesh2* m, apf::MeshEntity* tri);

};

#endif
