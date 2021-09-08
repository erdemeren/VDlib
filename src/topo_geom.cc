#include <utility>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_extinfo.h"
#include "topo_topo.h"
#include "topo_pca.h"

#include "topo_entlist.h"
#include "topo_geom.h"

#include <math.h>  //acos()

// Used in printing the set total measure.
enum GEOM_TYPE {
  VERTEX,
  LENGTH,
  AREA,
  VOLUME,
  END
};
const char* geom_type[] = {"vertex", "length", "area", "volume"}; 

//#define PI_L std::acos(-1)  /* pi */

// LOOK This probably exists somewhere in apf. The order of vertices, used for
// normal calculations, from region_faces.jpg in mds:
// Face0 outward normal is 10 x 02
// Face1 outward normal is 01 x 13
// Face2 outward normal is 23 x 31
// Face3 outward normal is 32 x 20

// Normal Calculation: Order of vertices for used for constructing the 
// direction vectors.
int vd_vert_order [4][3] = {{1,0,2}, {0,1,3}, {2,3,1}, {3,2,0}};
// Given triangle id t1 and vertex id v1, order of remaining vertices {v2, v3}  
// such that the cross product (v2-v1) x (v3-v1) points into the tetrahedron.
int vd_vert_norm [4][4][2] = {{{1, 2}, {2, 0}, {0, 1}, {-1, -1}}, 
                              {{3, 1}, {0, 3}, {-1, -1}, {1, 0}}, 
                              {{-1, -1}, {3, 2}, {1, 3}, {2, 1}}, 
                              {{2, 3}, {-1, -1}, {3, 0}, {0, 2}}};

// Given two indices of vertices on a tet, return the other two.
int vd_tet_vv [4][4][2] = {{{-1,-1}, {2,3}, {1,3}, {1,2}}, 
                           {{2,3}, {-1,-1}, {0,3}, {0,2}},
                           {{1,3}, {0,3}, {-1,-1}, {0,1}},
                           {{1,2}, {0,2}, {0,1}, {-1,-1}}};

// The bounding vertices of the edges
int vd_tet_edge_vert [6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
// Given an edge on a tet, return the ids of triangles not bounded by the edge
int vd_tet_edge_c_tri [6][2] = {{2,3},{1,3},{1,2},{0,2},{0,3},{0,1}};

// Given the vertex index, return the index of the other two. Not necessary but
// cleaner.
int lookup_triverts [3][2] = {{1,2},{0,2},{0,1}};
int lookup_edgeverts [2] = {1, 0};
int lookup_triedges [3][2] = {{2,0},{0,1},{1,2}};

// Given the vertex, return the triangle across.
int lookup_tet_surfx [4] = {2, 3, 1, 0};
// Given the triangle, return the vertex across.
int lookup_tet_xsurf [4] = {3, 2, 0, 1};


// Given the vertex, return the edge across.
int lookup_tri_edgex [3] = {1, 2, 0};

// Given a vertex on a tri, return the indices of the other two vertices.
int lookup_tri_x_x [3][2] = {{1, 2}, {0, 2}, {0, 1}};

// Functions that derive geometric information from objects, such as length, 
// position, area, volume, cosine angle, normal direction.

// apfVector.h defines possible vector operations that will be of use.

// Vector operations from apf:
// Given a edge mesh entity, calculate its euclidian norm.
// The apf way of doing this is the following:
// apf::MeshElement* edge_elem = createMeshElement(m, edge);
// measure(edge_elem);

apf::Vector3 norm_0(apf::Vector3 v_in) {
  double norm = v_in.getLength();
  if (norm < std::numeric_limits<double>::min() ) {
    return apf::Vector3(0,0,0);
  }
  return v_in.normalize();
}

//LEFT HERE (vd_get_center is well used, so before replacing update all uses)
apf::Vector3 vd_get_center(apf::Mesh2* m, std::vector<apf::MeshEntity*>* set_in) {
  apf::Vector3 center(0,0,0);

  assert(set_in->size() != 0);
  int ent_type = m->getType(set_in->at(0));
  int d = m->typeDimension[ent_type];
  assert(d == 0);

  for (int i = 0; i < set_in->size(); i++) {
    apf::Vector3 point(0,0,0);
    m->getPoint(set_in->at(i), 0, point);
    center = center + point;
  }

  return center/set_in->size();
}

// Given a set of entities return the average of their positions, weighted by
// the length/area/volume of the entities.
apf::Vector3 vd_get_center_e(apf::Mesh2* m, std::vector<apf::MeshEntity*>* set_in) {
  apf::Vector3 center(0,0,0);
  apf::Vector3 p_center(0,0,0);

  assert(set_in->size() != 0);

  double meas = 0;
  double meas_tot = 0;
  for (int i = 0; i < set_in->size(); i++) {
    apf::Vector3 point(0,0,0);
    point = getLinearCentroid(m, set_in->at(i));

    //apf::MeshElement* ee = createMeshElement(m, set_in->at(i));
    //meas = measure(ee);
    //destroyMeshElement(ee);
    meas = vd_meas_ent(m, set_in->at(i));

    p_center = p_center + point;
    center = center + point*meas;
    meas_tot = meas_tot + meas;
  }
  if(meas_tot < std::numeric_limits<double>::min())
    return p_center/set_in->size();

  return center/meas_tot;
}

// Given two edges joined at v_ctr, find the plane defined by the edges.
apf::Vector3 vd_get_e_plane(apf::Mesh2* m, apf::MeshEntity* v_ctr, 
                              apf::MeshEntity* e1, apf::MeshEntity* e2) {
  apf::Vector3 temp1(0,0,0);
  apf::Vector3 temp2(0,0,0);
  apf::Vector3 ctr(0,0,0);

  m->getPoint(v_ctr, 0, ctr);

  apf::Downward d_v;
  m->getDownward(e1, 0, d_v);
  int v_id = findIn(d_v, 2, v_ctr);
  assert(v_id > -1);
  v_id = (v_id + 1) % 2;
  m->getPoint(d_v[v_id], 0, temp1);
  temp1 = temp1 - ctr;
  assert(temp1.getLength() > std::numeric_limits<double>::min());

  m->getDownward(e2, 0, d_v);
  v_id = findIn(d_v, 2, v_ctr);
  assert(v_id > -1);
  v_id = (v_id + 1) % 2;
  m->getPoint(d_v[v_id], 0, temp2);
  temp2 = temp2 - ctr;

  assert(temp2.getLength() > std::numeric_limits<double>::min());

  return norm_0(vd_cross(temp1, temp2));
}

// Given a set of entities, return the total length, area or volume.
double vd_meas_ent(apf::Mesh2* m, apf::MeshEntity* e) {

  int ent_type = m->getType(e);
  int d = m->typeDimension[ent_type];
  double meas = 0;

  if (d == 0) {
    //printf("The %s is not measurable.\n", geom_type[d]);
  }

  else if (d == 2) {
    //printf("The %s is not measurable.\n", geom_type[d]);
    meas = std::fabs(vd_area_out_n(m, e).getLength());
  }

  else if (d == 3) {
    //printf("The %s is not measurable.\n", geom_type[d]);
    meas = vd_volume_tet(m, e);
  }

  else {
    apf::MeshElement* ee = createMeshElement(m, e);
    meas = measure(ee);
    destroyMeshElement(ee);

    //printf("The %s is %f.\n", geom_type[d], meas);
  }
  return meas;
}

double vd_meas_ent_dim(apf::Mesh2* m, int dim) {

  double total = 0.;
  if(dim == 0)
    return total;
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  apf::ModelEntity* mdl;
  while ((e = m->iterate(it))) {
    mdl = m->toModel(e);
    if(dim == m->getModelType(mdl)) {
      total += vd_meas_ent(m, e);
    }
  }
  m->end(it);

  return total;
}

// Given a set of vertices, return the maximum radius from the 
// center of the cloud.
double vd_meas_rad_set(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev) {
  apf::Vector3 rad(0,0,0);
  apf::Vector3 ctr(0,0,0);
  apf::Vector3 temp(0,0,0);

  for(int i = 0; i < ev->size(); i++) {
    m->getPoint(ev->at(i), 0, temp);
    ctr = ctr + temp;
  }
  ctr = ctr/ev->size();

  for(int i = 0; i < ev->size(); i++) {
    m->getPoint(ev->at(i), 0, temp);
    temp = ctr - temp;
    if(rad.getLength() < temp.getLength() )
      rad = temp;
  }

  return rad.getLength();
}

// Given a set of vertices, return the maximum radius from the 
// center of the cloud.
double vd_meas_rad_set_min(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev) {
  double rad = -1;
  apf::Vector3 ctr(0,0,0);
  apf::Vector3 temp(0,0,0);

  for(int i = 0; i < ev->size(); i++) {
    m->getPoint(ev->at(i), 0, temp);
    ctr = ctr + temp;
  }
  ctr = ctr/ev->size();

  for(int i = 0; i < ev->size(); i++) {
    m->getPoint(ev->at(i), 0, temp);
    temp = ctr - temp;
    if(rad < 0 or rad > temp.getLength() )
      rad = temp.getLength();
  }

  return rad;
}
// Given a set of vertices, return the minimum and maximum radii from the 
// center of the cloud, barring points too close to the closest point.
std::pair<double, double> vd_meas_rad_set_minmax(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev) {
  if(ev->size() > 0) {
    apf::Vector3 rad_max(0,0,0);
    apf::Vector3 rad_min(0,0,0);
    apf::Vector3 rad_min2(0,0,0);
    apf::Vector3 ctr(0,0,0);
    apf::Vector3 temp(0,0,0);
    apf::Vector3 temp_min(0,0,0);

    double vect[3] = {0,0,0};
    ctr.fromArray(vect);
    rad_max.fromArray(vect);
    rad_min.fromArray(vect);
    rad_min2.fromArray(vect);
    temp.fromArray(vect);
    temp_min.fromArray(vect);

    double r_max = -1;
    double r_min = -1;

    for(int i = 0; i < ev->size(); i++) {
      m->getPoint(ev->at(i), 0, temp);
      ctr = ctr + temp;
    }
    ctr = ctr/ev->size();
    //std::cout << "ctr " << ctr << std::endl;

    for(int i = 0; i < ev->size(); i++) {
      m->getPoint(ev->at(i), 0, temp);
      temp = temp - ctr;
      //std::cout << temp << " r " << temp.getLength() << std::endl;
      if(r_max < temp.getLength() ) {
        rad_max = temp;
        r_max = rad_max.getLength();
        //std::cout << "\trad_max " << rad_max << " r " << r_max << std::endl;
      }
      if(r_min < 0 or r_min > temp.getLength() ) {
        rad_min = temp;
        r_min = rad_min.getLength();
        //std::cout << "\trad_min " << rad_min << " r " << r_min << std::endl;
      }
    }

    assert(r_max > -1);
    assert(r_min > -1);

    double r_min2 = -1;
    for(int i = 0; i < ev->size(); i++) {
      m->getPoint(ev->at(i), 0, temp);
      temp = temp - ctr;
      temp_min = temp - rad_min;

      //std::cout << "temp_min " << temp_min 
      //          << " r " << temp_min.getLength() << std::endl;
      if(r_min2 <= temp_min.getLength() ) {
        rad_min2 = temp_min;
        r_min2 = rad_min2.getLength();
        //std::cout << "\trad_min2 " << rad_min2 
        //          << " r " << r_min2 << std::endl;
      }
    }

    double min_th = 2*r_min;
    if(min_th > r_min2)
      min_th = 0.99*r_min2;
    //std::cout << "min_th " << min_th << std::endl;

    r_min2 = -1;
    for(int i = 0; i < ev->size(); i++) {
      m->getPoint(ev->at(i), 0, temp);
      temp = temp - ctr;
      temp_min = temp - rad_min;

      //std::cout << "temp " << temp 
      //          << " r " << temp.getLength() << std::endl;
      //std::cout << "temp_min " << temp_min 
      //          << " r " << temp_min.getLength() << std::endl;
      if(min_th <= temp_min.getLength() ) {
        if(r_min2 < 0 or r_min2 > temp.getLength() ) {
          rad_min2 = temp;
          r_min2 = rad_min2.getLength();
          //std::cout << "\trad_min2 " << rad_min2 
          //          << " r " << r_min2 << std::endl;
        }
      }
    }

    assert(r_min2 > -1);
    rad_min2 = rad_min2 - rad_min;
    r_min2 = rad_min2.getLength();

    return std::make_pair(r_max, r_min2/2);
  }
  return std::make_pair(0, 0);
  
}

/*
apf::Vector3 ellip_ax_val_pos(ellipsoid* ellip, apf::Vector3 pos_rel) {

  apf::Vector3 proj;
  double z[3] = {0,0,0};
  proj.fromArray(z);
  for(int i = 0; i < ellip.dim; i++) {
    proj.at(i) = pos_rel*ellip.e_v.at(i);
  }

  sc = proj.getLength();

  for(int i = 0; i < 3; i++) {
    proj.at(i) = pos_rel*ellip.e_v.at(i);
  }

  for(int i = ellip.dim; i < 3; i++) {
  }

}
*/
ellipsoid vd_meas_rad_set_minmax(apf::Mesh2* m, int c_dim, int c_id) {
// Pick a point inside the cell. Using that point, find the weighted center of
// the cell interior. Get the weighted relative positions of the bounding d-1 
// dim entities. Do a pca using the weighted positions.
// Using the bounding vertices, find the maximum width along the principle 
// axes.
  assert(c_dim > 0 and c_dim < 4);

  ellipsoid ellip;

  apf::Vector3 pos_in(0,0,0);
  std::vector<apf::Vector3> weight(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> pos(0, apf::Vector3(0,0,0));
  double meas;
  double meas_tot;

  std::vector<apf::MeshEntity*> es(0);
  std::vector<apf::MeshEntity*> es_down(0);
  std::vector<apf::MeshEntity*> es_in(0);

  vd_find_ent_geom(m, &es, c_id, c_dim, c_dim);
  vd_find_ent_geom(m, &es_in, c_id, c_dim, c_dim-1);

  assert(es.size() > 0);

  pos_in = vd_get_center_e(m, &es);

  vd_set_down(m, &es, &es_down);
  vd_rem_cell(m, &es_down, c_id, c_dim);
  //double z[3] = {0,0,0};
  if(c_dim == 1) {
    if(es_down.size() > 0) {
      assert(es_down.size() == 2);
      pos.resize(es_down.size());
      for(int i = 0; i < es_down.size(); i++) {
        m->getPoint(es_down.at(i), 0, pos.at(i));
        pos.at(i) = pos.at(i) - pos_in;
      }
    }
    else {
      pos.resize(es_in.size());
      for(int i = 0; i < es_in.size(); i++) {
        m->getPoint(es_in.at(i), 0, pos.at(i));
        pos.at(i) = pos.at(i) - pos_in;
      }
    }
  }
  else if(c_dim == 2) {
    meas_tot = 0;
    if(es_down.size() > 0) {
      weight.resize(es_down.size());
      pos.resize(es_down.size());
      for(int i = 0; i < es_down.size(); i++) {
        apf::MeshElement* ee = createMeshElement(m, es_down.at(i));
        meas = measure(ee);
        destroyMeshElement(ee);

        pos.at(i) = (getLinearCentroid(m, es_down.at(i)) - pos_in)*meas;
        weight.at(i) = norm_0(pos.at(i))*meas;
        meas_tot = meas_tot + meas;
      }
    }
    else {
      vd_get_apos(m, pos_in, &es_in, &weight, &pos);
      for(int i = 0; i < es_down.size(); i++) {
        meas = norm_0(pos.at(i))*weight.at(i);
        pos.at(i) = pos.at(i)*meas;
        meas_tot = meas_tot + meas;
      }
    }

    for(int i = 0; i < pos.size(); i++) {
      pos.at(i) = pos.at(i)/meas_tot;
    }
  }
  else {
    assert(es_down.size() > 0);
    vd_get_apos(m, pos_in, &es_down, &weight, &pos);

    meas_tot = 0;
    for(int i = 0; i < es_down.size(); i++) {
      meas = norm_0(pos.at(i))*weight.at(i);
      pos.at(i) = pos.at(i)*meas;
      meas_tot = meas_tot + meas;
    }

    for(int i = 0; i < pos.size(); i++) {
      pos.at(i) = pos.at(i)/meas_tot;
    }
  }

  std::vector<apf::Vector3> e_v(0, apf::Vector3(0,0,0));
  std::vector<double> e_i(0);

  ellip = pca(&pos, (unsigned int) (c_dim));

  std::vector<apf::MeshEntity*> es_v(0);
  vd_find_ent_geom(m, &es_v, c_id, c_dim, 0);

  if(es_down.size() > 0) {
    std::vector<std::vector<apf::MeshEntity*> > es_temp
          (0, std::vector<apf::MeshEntity*>(0));

    es_temp.resize(c_dim);
    vd_set_down(m, &es, &es_temp.at(c_dim-1));
    for(int i = c_dim-1; i > 0; i--) {
      vd_set_down(m, &es_temp.at(i), &es_temp.at(i-1));
    }
    // Remove the internal vertices(es_v) from all vertices.
    vd_remove_set(&es_temp.at(0), &es_v);

    pos.clear();
    pos.resize(es_temp.at(0).size());

    for(int i = 0; i < es_temp.at(0).size(); i++) {
      m->getPoint(es_temp.at(0).at(i), 0, pos.at(i));
      pos.at(i) = pos.at(i) - pos_in;
    }

  }
  else {
    pos.clear();
    pos.resize(es_v.size());

    for(int i = 0; i < es_v.size(); i++) {
      m->getPoint(es_v.at(i), 0, pos.at(i));
      pos.at(i) = pos.at(i) - pos_in;
    }
  }

  for(int i = 0; i < ellip.e_v.size(); i++) {
    double r_max = 0;

    for(int j = 0; j < pos.size(); j++) {
      //apf::Vector3 eig_val = ellip_ax_val_pos(&ellip, pos.at(j));
      double r_curr = std::fabs(pos.at(j)*ellip.e_v.at(i));
      if(r_max < r_curr)
        r_max = r_curr;
    }
    ellip.e_i.at(i) = r_max;
  }

  ellip.ctr = pos_in;
  return ellip;
}

ellipsoid vd_meas_rad_set_minmax(apf::Mesh2* m, vd_entlist* e_list, int c_dim, 
                                                  int c_id) {
// Pick a point inside the cell. Using that point, find the weighted center of
// the cell interior. Get the weighted relative positions of the bounding d-1 
// dim entities. Do a pca using the weighted positions.
// Using the bounding vertices, find the maximum width along the principle 
// axes.
  assert(c_dim > 0 and c_dim < 4);

  ellipsoid ellip;

  apf::Vector3 pos_in(0,0,0);
  std::vector<apf::Vector3> weight(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> pos(0, apf::Vector3(0,0,0));
  double meas;
  double meas_tot;

  std::vector<apf::MeshEntity*> es_down(0);
  std::vector<apf::MeshEntity*>* es_in = 
                                &e_list->e.at(c_dim).at(c_id-1).at(c_dim-1);


  assert(e_list->e.at(c_dim).at(c_id-1).at(c_dim).size() > 0);

  pos_in = vd_get_center_e(m, &e_list->e.at(c_dim).at(c_id-1).at(c_dim));

  vd_set_down(m, &e_list->e.at(c_dim).at(c_id-1).at(c_dim), &es_down);
  vd_rem_cell(m, &es_down, c_id, c_dim);
  //double z[3] = {0,0,0};
  if(c_dim == 1) {
    if(es_down.size() > 0) {
      assert(es_down.size() == 2);
      pos.resize(es_down.size());
      for(int i = 0; i < es_down.size(); i++) {
        m->getPoint(es_down.at(i), 0, pos.at(i));
        pos.at(i) = pos.at(i) - pos_in;
      }
    }
    else {
      pos.resize(es_in->size());
      for(int i = 0; i < es_in->size(); i++) {
        m->getPoint(es_in->at(i), 0, pos.at(i));
        pos.at(i) = pos.at(i) - pos_in;
      }
    }
  }
  else if(c_dim == 2) {
    meas_tot = 0;
    if(es_down.size() > 0) {
      weight.resize(es_down.size());
      pos.resize(es_down.size());
      for(int i = 0; i < es_down.size(); i++) {
        apf::MeshElement* ee = createMeshElement(m, es_down.at(i));
        meas = measure(ee);
        destroyMeshElement(ee);

        pos.at(i) = (getLinearCentroid(m, es_down.at(i)) - pos_in)*meas;
        weight.at(i) = norm_0(pos.at(i))*meas;
        meas_tot = meas_tot + meas;
      }
    }
    else {
      vd_get_apos(m, pos_in, es_in, &weight, &pos);
      for(int i = 0; i < es_down.size(); i++) {
        meas = norm_0(pos.at(i))*weight.at(i);
        pos.at(i) = pos.at(i)*meas;
        meas_tot = meas_tot + meas;
      }
    }

    for(int i = 0; i < pos.size(); i++) {
      pos.at(i) = pos.at(i)/meas_tot;
    }
  }
  else {
    assert(es_down.size() > 0);
    vd_get_apos(m, pos_in, &es_down, &weight, &pos);

    meas_tot = 0;
    for(int i = 0; i < es_down.size(); i++) {
      meas = norm_0(pos.at(i))*weight.at(i);
      pos.at(i) = pos.at(i)*meas;
      meas_tot = meas_tot + meas;
    }

    for(int i = 0; i < pos.size(); i++) {
      pos.at(i) = pos.at(i)/meas_tot;
    }
  }

  std::vector<apf::Vector3> e_v(0, apf::Vector3(0,0,0));
  std::vector<double> e_i(0);

  ellip = pca(&pos, (unsigned int) (c_dim));

  std::vector<apf::MeshEntity*>* es_v =
                                &e_list->e.at(c_dim).at(c_id-1).at(0);

  if(es_down.size() > 0) {
    std::vector<std::vector<apf::MeshEntity*> > es_temp
          (0, std::vector<apf::MeshEntity*>(0));

    es_temp.resize(c_dim);
    vd_set_down(m, &e_list->e.at(c_dim).at(c_id-1).at(c_dim), 
                                                  &es_temp.at(c_dim-1));
    for(int i = c_dim-1; i > 0; i--) {
      vd_set_down(m, &es_temp.at(i), &es_temp.at(i-1));
    }
    // Remove the internal vertices(es_v) from all vertices.
    vd_remove_set(&es_temp.at(0), es_v);

    pos.clear();
    pos.resize(es_temp.at(0).size());

    for(int i = 0; i < es_temp.at(0).size(); i++) {
      m->getPoint(es_temp.at(0).at(i), 0, pos.at(i));
      pos.at(i) = pos.at(i) - pos_in;
    }

  }
  else {
    pos.clear();
    pos.resize(es_v->size());

    for(int i = 0; i < es_v->size(); i++) {
      m->getPoint(es_v->at(i), 0, pos.at(i));
      pos.at(i) = pos.at(i) - pos_in;
    }
  }

  for(int i = 0; i < ellip.e_v.size(); i++) {
    double r_max = 0;

    for(int j = 0; j < pos.size(); j++) {
      //apf::Vector3 eig_val = ellip_ax_val_pos(&ellip, pos.at(j));
      double r_curr = std::fabs(pos.at(j)*ellip.e_v.at(i));
      if(r_max < r_curr)
        r_max = r_curr;
    }
    ellip.e_i.at(i) = r_max;
  }

  ellip.ctr = pos_in;
  return ellip;
}

// Given a set of entities, return the total length, area or volume.
double vd_meas_set(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es_in) {

  int ent_type = m->getType(es_in->at(0));
  int d = m->typeDimension[ent_type];
  double meas = 0;

  if (d == 0) {
    //printf("The %s is not measurable.\n", geom_type[d]);
    return meas;  
  }

  else if (d == 3) {
    for( int i = 0; i < es_in->size(); i++) {
      meas += vd_volume_tet(m,es_in->at(i));
    }

    //printf("The %s is %f.\n", geom_type[d], meas);
    return meas;
  }
  else {
    for( int i = 0; i < es_in->size(); i++) {
      apf::MeshElement* ee = createMeshElement(m, es_in->at(i));
      meas += measure(ee);
      destroyMeshElement(ee);
    }

    //printf("The %s is %f.\n", geom_type[d], meas);
    return meas;
  }
}


double vd_avg_meas(apf::Mesh2* m, std::vector<apf::MeshEntity*>* e_set) {

  double avg = 0;

  for(int i = 0; i < e_set->size(); i++) {
    apf::MeshElement* ee = createMeshElement(m, e_set->at(i));
    avg = avg + measure(ee);
    destroyMeshElement(ee);
  }

  return avg/e_set->size();

}

double vd_minmax_meas(apf::Mesh2* m, std::vector<apf::MeshEntity*>* e_set, bool min) {

  int c = 1;
  if (!min)
    c = -1;

  if(e_set->size() == 0)
    return 0;

  apf::MeshElement* ee = createMeshElement(m, e_set->at(0));
  double len_cross = c*measure(ee);
  destroyMeshElement(ee);
  for (int i = 1; i < e_set->size(); i++) {
    ee = createMeshElement(m, e_set->at(i));
    double len1 = c*measure(ee);
    destroyMeshElement(ee);
    if (len1 < len_cross)
      len_cross = len1;
  }
  len_cross = len_cross*c;

  return len_cross;

}

double vd_minmax_meas(apf::Mesh2* m, std::vector<apf::MeshEntity*>* e_tet, apf::MeshEntity* v, bool min) {

  apf::Vector3 v_pos(0,0,0);
  apf::Vector3 a_pos(0,0,0);
  apf::Vector3 area(0,0,0);

  m->getPoint(v, 0, v_pos);

  double dist = -1;
  for(int i = 0; i < e_tet->size(); i++) {
    apf::Downward d_v;
    apf::Downward d_s;
    m->getDownward(e_tet->at(i), 0, d_v);
    m->getDownward(e_tet->at(i), 2, d_s);

    int v_id = findIn(d_v, 4, v);

    a_pos = getLinearCentroid(m, d_s[lookup_tet_surfx[v_id]]);

    area = vd_area_out(m, d_s[lookup_tet_surfx[v_id]]);
    a_pos = a_pos - v_pos;

    double dist_curr = std::fabs(norm_0(area)*a_pos);

    if(min) {
      if(dist_curr < dist or dist < 0)
        dist = dist_curr;
    }
    else {
      if(dist_curr > dist)
        dist = dist_curr;
    }
  }

  return dist;

}

// Given a set of triangles joined at a vertex, find the negative of the 
// gradient of the total area, wrt. the vertex position.
apf::Vector3 vd_neg_grad(apf::Mesh2* m, std::vector<apf::MeshEntity*>* tri, apf::MeshEntity* vert) {
  apf::Vector3 neg_grad(0,0,0);
  apf::Downward dv;

  apf::Vector3 w1(0,0,0);
  apf::Vector3 w2(0,0,0);
  apf::Vector3 v(0,0,0);
  apf::Vector3 area_out(0,0,0);

  m->getPoint(vert, 0, v);

  for (int i = 0; i < tri->size(); i++) {

    m->getDownward(tri->at(i), 0, dv);
    int e1 = findIn(dv, 3, vert);
    m->getPoint(dv[lookup_triverts[e1][0]], 0, w1);
    m->getPoint(dv[lookup_triverts[e1][1]], 0, w2);
    area_out = norm_0(cross(w1 - v, w2 - v));
    neg_grad = neg_grad + cross(area_out, w1-w2);
    //apf::Vector3 dist = (w1+w2)/2;
    //neg_grad = neg_grad + cross(area_out, (w1-w2)/dist.getLength());
  }
  return neg_grad;
}

//-----------------------------------


// print the adjacencies of model entities of a given dimension. Can be 
// used as a template.
void vd_print_adj(apf::Mesh2* m, int dim) {
  gmi_model* mdl = m->getModel();
  gmi_iter* mdl_it = gmi_begin(mdl, dim);

  while(gmi_ent* em = gmi_next(mdl, mdl_it)) {
    printf("The model entity has tag %d.\n", gmi_tag(mdl, em));
    gmi_set* em_adj = gmi_adjacent(mdl, em, dim-1);
    printf("It has %d %d-dimensional adjacencies.", em_adj->n, dim-1);
    for (int i = 0; i < em_adj->n; i++) {
      printf("The %dth one has tag %d.", i, gmi_tag(mdl, em_adj->e[i]));
    }
    printf("\n");
    gmi_free_set(em_adj);
  }
  gmi_end(mdl, mdl_it);
}

// Given an entity set of cube corner vertices, return the center point.
apf::Vector3 vd_get_center(apf::Mesh2* m, apf::MeshEntity* ent) {
  apf::Downward down;
  int downwardCount = m->getDownward(ent, 0, down);
  Entity_set es;
  copy_ent_set(&es, down, downwardCount);
  return vd_get_center(m, &es);
}

// Given an entity set of cube corner vertices, return the center point.
apf::Vector3 vd_get_center(apf::Mesh2* m, Entity_set* set_in) {
  apf::Vector3 center(0,0,0);
  double vect[3] = {0,0,0};
  center.fromArray(vect);

  for (int i = 0; i < set_in->n; i++) {
    apf::Vector3 point(0,0,0);
    m->getPoint(set_in->e[i], 0, point);
    center = center + point;
  }
  return center/set_in->n;
}

// Adapted from getAverageEdgeLength in maSize.
// Get the average entity length, area or volume.
double getAverageEntSize(apf::Mesh2* m, int dim) {
  PCU_Barrier();

  if (dim == 0 or dim > 3)
    return 0;

  double sums[2];
  double& width_sum = sums[0];
  double& ent_count = sums[1];
  width_sum = 0;
  ent_count = 0;
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  while ((e = m->iterate(it)))
  {
    apf::MeshElement* ee = createMeshElement(m, e);
    width_sum += measure(ee);
    ent_count += 1.0;
    destroyMeshElement(ee);
  }
  m->end(it);
  PCU_Barrier();
  PCU_Add_Doubles(sums,2);
  if(ent_count == 0)
    return 0;
  else
    return width_sum / ent_count;
}

// Using a histogram of N bins, estimate the median of the entity 
// length / area / volume of given dimension. First find the maximum, and
// scale the histogram by using the maximum to find bin widths. Count the
// number of entities in each bin. 
// Assumes single mesh. TODO ask scorec maintainers about other best practices 
// for collecting info from multiple meshes. 
double getMedianEntSize_hist(apf::Mesh2* m, int dim, int N) {
  if (dim == 0 or dim > 3)
    return 0;

  std::vector<double> count(N, 0);
  double max_val = 0;
  double temp = 0;
  long int ent_count = 0;

  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  while ((e = m->iterate(it)))
  {
    apf::MeshElement* ee = createMeshElement(m, e);
    temp = measure(ee);
    if(temp > max_val)
      max_val = temp;
    ent_count += 1;
    destroyMeshElement(ee);
  }
  m->end(it);

  if(ent_count == 0)
    return 0;
  else {
    it = m->begin(dim);
    while ((e = m->iterate(it)))
    {
      apf::MeshElement* ee = createMeshElement(m, e);
      temp = measure(ee);
      int bin_id;
      if(std::fabs(temp-max_val) < std::numeric_limits<double>::epsilon())
        bin_id = N-1;
      else
        bin_id = floor(temp/max_val*N);
      count.at(bin_id) = count.at(bin_id) + 1;
      destroyMeshElement(ee);
    }
    m->end(it);

    long int nbr = 0;
    for(int i = 0; i < count.size(); i++) {
      nbr = nbr + count.at(i);
      if(nbr > ent_count/2) {
        if(ent_count%2 == 0)
          return max_val/N*(i+1)/2;
        else
          return max_val/N*i;
      }      
    }
    assert(nbr > ent_count/2);
  }
}

// Assumes single mesh. TODO ask scorec maintainers about other best practices 
// for collecting info from multiple meshes. 
double getMedianEntSize_hist(std::vector<double> &sz, int N) {
  std::vector<double> count(N, 0);
  double temp = 0;
  long int ent_count = 0;

  double max_val = 0;
  double min_val = 0;
  if(sz.size() > 0) {
    min_val = sz.at(0);
    max_val = sz.at(0);
  }

  for(int i = 0; i < sz.size(); i++) {
    if(sz.at(i) > max_val)
      max_val = sz.at(i);
    if(sz.at(i) < min_val)
      min_val = sz.at(i);
    ent_count += 1;
  }
  max_val = max_val - min_val;

  if(ent_count == 0)
    return 0;
  else {
    for(int i = 0; i < sz.size(); i++) {
      double temp = sz.at(i) - min_val;
      int bin_id;
      if(std::fabs(temp-max_val) < std::numeric_limits<double>::epsilon())
        bin_id = N-1;
      else
        bin_id = floor(temp/max_val*N);
      count.at(bin_id) = count.at(bin_id) + 1;
    }

    long int nbr = 0;
    for(int i = 0; i < count.size(); i++) {
      nbr = nbr + count.at(i);
      if(nbr > ent_count/2) {
        if(ent_count%2 == 0)
          return min_val + max_val/N*(i+1)/2;
        else
          return min_val + max_val/N*i;
      }      
    }
    assert(nbr > ent_count/2);
  }
}

/*
double getMedianEntSize_hist(apf::Mesh2* m, int dim, int N) {
  if (dim == 0 or dim > 3)
    return 0;

  double count[N] = {0};
  double max_val = 0;
  double temp = 0;
  long int ent_count = 0;

  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  while ((e = m->iterate(it)))
  {
    apf::MeshElement* ee = createMeshElement(m, e);
    temp = measure(ee);
    if(temp > max_val)
      max_val = temp;
    ent_count += 1;
    destroyMeshElement(ee);
  }
  m->end(it);
  PCU_Max_Double(max_val);
  PCU_Add_Long(ent_count);

  if(ent_count == 0)
    return 0;
  else {
    it = m->begin(dim);
    while ((e = m->iterate(it)))
    {
      apf::MeshElement* ee = createMeshElement(m, e);
      temp = measure(ee);
      int bin_id = floor(temp/max_val*N);
      count[bin_id] = count[bin_id] + 1;
      destroyMeshElement(ee);
    }
    m->end(it);

    return width_sum / ent_count;
  }
}
*/

// Adapted from getAverageEdgeLength in maSize.
// Get the average entity length, area or volume.
double getAverageEntSize(apf::Mesh2* m, apf::MeshEntity* v, int dim) {
  if (dim == 0 or dim > 3)
    return 0;

  double sums[2];
  double& width_sum = sums[0];
  double& ent_count = sums[1];
  width_sum = 0;
  ent_count = 0;

  std::vector< std::vector<apf::MeshEntity*> > ent
          (0, std::vector<apf::MeshEntity*>(0));

  ent.resize(dim);

  vd_set_up(m, v, &ent.at(0));
  for(int i = 1; i < dim; i++)
    vd_set_up(m, &ent.at(i-1), &ent.at(i));

  for(int i = 0; i < ent.at(dim-1).size(); i++) {
    apf::MeshElement* ee = createMeshElement(m, ent.at(dim-1).at(i));
    width_sum += measure(ee);
    ent_count += 1.0;
    destroyMeshElement(ee);
  }

  PCU_Add_Doubles(sums,2);
  if(ent_count == 0)
    return 0;
  else
    return width_sum / ent_count;
}

// Given the geometry dimension and tag, return the size of the geometry:
// total length, area or volume.
double vd_meas_geom(apf::Mesh2* m, int em_dim, int em_tag) {

  double meas = 0;

  // Check if dimensions are proper.
  if ((em_dim == 0) || (em_dim > 3)) {
    printf("The dimension %s is not appropriate.\n", geom_type[em_dim]);
    return meas;  
  }
  else {

    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(em_dim);
    
    while( (e = m->iterate(it))) {

      int ent_type = m->getType(e);
      apf::ModelEntity* em = m->toModel(e);
      int d = m->getModelType(em);
      int geom = m->getModelTag(em);

      // Check if the matched geometry is of the same dimension. 
      if ((d == em_dim) && (geom == em_tag)) {
        apf::MeshElement* ee = createMeshElement(m, e);
        meas = meas + measure(ee);
        destroyMeshElement(ee);
      }
    }
    m->end(it);
    //printf("The %s is %f.\n", geom_type[em_dim], meas);
    return meas;
  }

}

// Given the geometry dimension and tag, return the largest distance of the 
// vertices on bounding cells from the center.
double vd_rad_geom(apf::Mesh2* m, int em_dim, int em_tag) {

  std::vector<apf::MeshEntity*> es(0);
  std::vector<apf::MeshEntity*> ev(0);
  std::vector<apf::MeshEntity*> ev_d(0);

  vd_find_ent_geom(m, &es, em_tag, em_dim, em_dim);
  vd_find_ent_geom(m, &ev_d, em_tag, em_dim, 0);
  vd_vert_set(m, &es, &ev);

  vd_remove_set(&ev, &ev_d);
  return vd_meas_rad_set(m, &ev);
}

// TODO To count over the entities twice is unnecessary.
double vd_meas_geom_avg(apf::Mesh2* m, int em_dim, int em_tag) {

  std::vector<apf::MeshEntity*> es_edge(0);
  vd_find_ent_geom(m, &es_edge, em_tag, em_dim, em_dim);
  return vd_meas_geom(m, em_dim, em_tag)/es_edge.size();

}

// Given a set of entities, return the total length, area or volume.
double vd_meas_set(apf::Mesh2* m, Entity_set* es_in) {

  int ent_type = m->getType(es_in->e[0]);
  int d = m->typeDimension[ent_type];
  double meas = 0;

  if (d == 0) {
    //printf("The %s is not measurable.\n", geom_type[d]);
    return meas;  
  }

  else {
    for( int i = 0; i<es_in->n; i++) {
      apf::MeshElement* ee = createMeshElement(m, es_in->e[i]);
      meas += measure(ee);
      destroyMeshElement(ee);
    }

    //printf("The %s is %f.\n", geom_type[d], meas);
    return meas;
  }
}

// Given a set of entities, return the number of negative length, area, volume 
// entities.
long int vd_chk_neg_vert(apf::Mesh2* m, apf::MeshEntity* vert) {


  long int n = 0;
  std::vector<std::vector<apf::MeshEntity* > > e_inv(4, 
                                    std::vector<apf::MeshEntity*>(1) );
  e_inv.at(0).at(0) = vert;
  for(int k = 0; k < 3; k++) {
    vd_set_up(m, &e_inv.at(k), &e_inv.at(k+1));
  }

  bool inv = false;
  for (int k = 0; k < e_inv.at(3).size(); k++) {
    double meas_ee = vd_volume_tet(m, e_inv.at(3).at(k));
    if(meas_ee < std::numeric_limits<double>::min()) {

      int type1 = m->getModelType(m->toModel(e_inv.at(3).at(k)));
      int tag1 = m->getModelTag(m->toModel(e_inv.at(3).at(k)));
      printf("Element %p %dcell%d inverts, vol: %.5f.\n", 
                          (void*)e_inv.at(3).at(k), type1,tag1,meas_ee);
      std::cout << getLinearCentroid(m, e_inv.at(3).at(k)) << std::endl;
      inv = true;
      n++;
    }
  }

  return n;
}

// Given a set of entities, return the number of negative length, area, volume 
// entities.
long int vd_chk_neg(apf::Mesh2* m, std::vector<apf::MeshEntity*>* elem) {

  int ent_type = m->getType(elem->at(0));
  int d = m->typeDimension[ent_type];

  double meas = 0;
  long int n = 0;
  for( int i = 0; i< elem->size(); i++) {
    meas = vd_volume_tet(m, elem->at(i));
    if (meas < std::numeric_limits<double>::min()) {

      apf::ModelEntity* mdl = m->toModel(elem->at(i));
      std::cout << "Elem " << elem->at(i) << " " << m->getModelType(mdl) << "c"
            << m->getModelTag(mdl) << " is inverted with vol "
            << meas << std::endl;
      n++;
      print_tet(m, elem->at(i));
    }
  }

  return n;
}


// Given a set of entities, return the number of negative length, area, volume 
// entities.
long int vd_chk_neg(apf::Mesh2* m, apf::MeshEntity* elem) {

  double meas = vd_volume_tet(m, elem);

  if (meas < std::numeric_limits<double>::min()) {
    apf::ModelEntity* mdl = m->toModel(elem);
    std::cout << "Elem " << elem << " " << m->getModelType(mdl) << "c"
          << m->getModelTag(mdl) << " is inverted with vol "
          << meas << std::endl;
    print_tet(m, elem);
    return 1;
  }

  return 0;
}
// Given a set of entities, return the number of negative length, area, volume 
// entities.
long int vd_chk_neg(apf::Mesh2* m, Entity_set* elem) {

  int ent_type = m->getType(elem->e[0]);
  int d = m->typeDimension[ent_type];

  double meas = 0;
  long int n = 0;
  for( int i = 0; i< elem->n; i++) {
    meas = vd_volume_tet(m, elem->e[i]);

    if (meas < std::numeric_limits<double>::min()) {
      apf::ModelEntity* mdl = m->toModel(elem->e[i]);
      std::cout << "Elem " << elem->e[i] << " " << m->getModelType(mdl) << "c"
            << m->getModelTag(mdl) << " is inverted with vol "
            << meas << std::endl;
      print_tet(m, elem->e[i]);
      n++;
    }
  }

  return n;
}

// Given a mesh, return the number of negative volume entities.

long int vd_chk_neg(apf::Mesh2* m) {

  long int n = 0;
  double meas = 0;
  apf::MeshEntity* elem;

  apf::MeshIterator* it = m->begin(3);
  while(elem = m->iterate(it)) {

    meas = vd_volume_tet(m, elem);

    if (meas < std::numeric_limits<double>::min()) {
      apf::ModelEntity* mdl = m->toModel(elem);
      std::cout << "Elem " << elem << " " << m->getModelType(mdl) << "c"
            << m->getModelTag(mdl) << " is inverted with vol "
            << meas << std::endl;

      apf::Downward down;
      m->getDownward(elem, 0, down);
      bool sgn = vd_volume_tet_sign(m, down);
      std::cout << "Sign is " << sgn
                << " " <<  vd_volume_tet(m, down)
                << std::endl;
      //assert(!sgn);
      printf("Elem %p is inverted with vol %.5f.\n", (void*)elem, meas);
      n++;
      print_tet(m, elem);
    }
  }

  return n;
}

long int vd_chk_neg_sgn(apf::Mesh2* m) {

  long int n = 0;
  double meas = 0;
  apf::MeshEntity* elem;

  apf::MeshIterator* it = m->begin(3);
  while(elem = m->iterate(it)) {

    apf::Downward down;
    m->getDownward(elem, 0, down);
    bool sgn = vd_volume_tet_sign(m, down);
    if (!sgn) {
      double meas = vd_volume_tet(m, down);
      apf::ModelEntity* mdl = m->toModel(elem);
      std::cout << "Elem " << elem << " " << m->getModelType(mdl) << "c"
            << m->getModelTag(mdl) << " is inverted with vol "
            << meas << std::endl;

      //assert(!sgn);
      printf("Elem %p is inverted with vol %.5f.\n", (void*)elem, meas);
      n++;
      print_tet(m, elem);
    }
  }

  return n;
}

// Given and edge and a direction, check if splitting the edge and moving the  
// vertex to the position inverts any tetrahedra bounded by the edge.
bool vd_edge_inv(apf::Mesh2* m, apf::MeshEntity* edge, apf::Vector3 pos) {

  int ent_type = m->getType(edge);
  int d = m->typeDimension[ent_type];
  assert (d == 1);

  std::vector<apf::MeshEntity*> tris(0);
  std::vector<apf::MeshEntity*> tets(0);

  vd_set_up(m, edge, &tris);
  vd_set_up(m, &tris, &tets);

  apf::Downward d_e;
  apf::Downward d_t;

  apf::Vector3 i_ctr(0,0,0);
  apf::Vector3 area(0,0,0);
  apf::Vector3 a_pos(0,0,0);
  i_ctr = vd_get_pos(m, edge);

  for(int i = 0; i < tets.size(); i++) {
    m->getDownward(tets.at(i), 2, d_t);
    m->getDownward(tets.at(i), 1, d_e);
    int e1 = findIn(d_e, 6, edge);
    assert(e1 > -1);

    apf::MeshEntity* tri_curr = d_t[vd_tet_edge_c_tri[e1][0]];
    a_pos = vd_get_pos(m, tri_curr);
    area = vd_area_out_n(m, tri_curr);

    if((area*(pos-a_pos))*(area*(i_ctr-a_pos)) < 
        std::numeric_limits<double>::min())
      return true;

    tri_curr = d_t[vd_tet_edge_c_tri[e1][1]];
    a_pos = vd_get_pos(m, tri_curr);
    if((area*(pos-a_pos))*(area*(i_ctr-a_pos)) < 
        std::numeric_limits<double>::min())
      return true;

  }
  return false;
}

// Print the normals for elements inside the mesh.
void vd_print_norm(apf::Mesh2* m) {

  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* elem;
  while(elem = m->iterate(it)) {
    printf("%p\n",(void*)elem);
    int ent_type = m->getType(elem);
    int d = m->typeDimension[ent_type];
    printf("Entity dimension is %d.\n", d);
    if (d == 3) {
      apf::Downward s_down;
      apf::Downward v_down;

      int size = m->getDownward(elem, d-1, s_down);
      m->getDownward(elem, d-3, v_down);
      // Face0 outward normal is 10 x 02
      // Face1 outward normal is 01 x 13
      // Face2 outward normal is 23 x 31
      // Face3 outward normal is 32 x 20

      apf::Vector3 point0(0,0,0);
      apf::Vector3 point1(0,0,0);
      apf::Vector3 point2(0,0,0);
      apf::Vector3 point3(0,0,0);
      m->getPoint(v_down[0], 0, point0);
      m->getPoint(v_down[1], 0, point1);
      m->getPoint(v_down[2], 0, point2);
      m->getPoint(v_down[3], 0, point3);

      printf("For the element %p the outward pointing normals are:\n", 
                                                            (void*)elem);
      std::cout << s_down[0]<< " has the outward normal "
      << norm_0(cross(point0 - point1, point2 - point0)) <<std::endl;

      std::cout << s_down[1]<< " has the outward normal "
      << norm_0(cross(point1 - point0,point3 - point1)) <<std::endl;

      std::cout << s_down[2]<< " has the outward normal "
      << norm_0(cross(point3 - point2, point1 - point3)) <<std::endl;

      std::cout << s_down[3]<< " has the outward normal "
      << norm_0(cross(point2 - point3,point0 - point2)) <<std::endl;
    }
    else {
      printf("Dimension %d inappropriate.\n", d);
    }
  }
/*
  Entity_set es_elem = Entity_set();
  init_ent_set(&es_elem, elem);
  vd_print_down(m, &es_elem);
*/

  m->end(it);

}

// Given a point and an edge, find the projection of the point on the edge.
apf::Vector3 vd_projpt(apf::Mesh2* m, apf::MeshEntity* edge, apf::Vector3 pt) {
  apf::Vector3 e1(0,0,0);
  apf::Vector3 e2(0,0,0);

  apf::Downward d_v;
  m->getDownward(edge, 0, d_v);

  m->getPoint(d_v[0], 0, e1);
  m->getPoint(d_v[1], 0, e2);

  e2 = e2 - e1;
  e2 = norm_0(e2);

  return e1 + e2*(e2*(pt - e1 ) );
}

// Return the projection of pt on line e2-e1.
apf::Vector3 vd_projpt(apf::Vector3 e1, apf::Vector3 e2, apf::Vector3 pt) {
  e2 = e2 - e1;
  e2 = norm_0(e2);

  return e1 + e2*(e2*(pt - e1 ) );
}

// Given a volume geometry tag, and a boundary surface triangle, 
// return cross product of the edge vectors, pointing outwards from that 
// volume. It assumes the element entities belong to different geometries. 
apf::Vector3 vd_area_out(apf::Mesh2* m, apf::MeshEntity* surf, int geom) {
  apf::Up elem_up;

  m->getUp(surf, elem_up);

  apf::Vector3 dA(0,0,0);

  assert((elem_up.n == 1) || (elem_up.n == 2));

  apf::Downward s_down;
  apf::Downward v_down;

  m->getDownward(elem_up.e[0], 2, s_down);
  m->getDownward(elem_up.e[0], 0, v_down);

  int i1 = findIn(s_down, 4, surf);

  apf::Vector3 point0(0,0,0);
  apf::Vector3 point1(0,0,0);
  apf::Vector3 point2(0,0,0);
  apf::Vector3 proj(0,0,0);

  m->getPoint(v_down[vd_vert_order[i1][0]], 0, point0);
  m->getPoint(v_down[vd_vert_order[i1][1]], 0, point1);
  m->getPoint(v_down[vd_vert_order[i1][2]], 0, point2);

  proj = vd_projpt(point0, point1, point2);
  proj = point2 - proj;

  point1 = point1 - point0;
  double l1 = point1.getLength();
  double l2 = proj.getLength();

  proj = norm_0(proj);
  point1 = norm_0(point1);

  //dA = cross(proj, point1)*0.5; !Incorrect ordering
  dA = cross(point1, proj)*(0.5*l1*l2);

  apf::ModelEntity* em = m->toModel(elem_up.e[0]);
  int em_type = m->getModelType(em);
  int em_tag = m->getModelTag(em);

  if (elem_up.n == 1) {
    // printf("Domain surface triangle.\n");
    return dA;
  }

  else if (elem_up.n == 2) {
    // printf("Internal surface triangle.\n");
    if (em_tag == geom)
      return dA;
    else
      return dA*(-1.0);
  }

}

// Given a volume geometry tag, and a boundary surface triangle, 
// return cross product of area normal with the difference of edge vectors to 
// obtain the gradient of the area w.r.t. the vertex position. 
// It assumes the element entities belong to different geometries. 
// TODO The computational complexity can be reduced by potentially simplifying
// vector calculations.
apf::Vector3 vd_grad(apf::Mesh2* m, apf::MeshEntity* vert, apf::MeshEntity* surf, int geom) {
  apf::Up elem_up;

  m->getUp(surf, elem_up);

  apf::Vector3 dA(0,0,0);

  assert((elem_up.n == 1) || (elem_up.n == 2));

  apf::Downward s_down;
  apf::Downward v_down;

  m->getDownward(elem_up.e[0], 2, s_down);
  m->getDownward(elem_up.e[0], 0, v_down);

  int i1 = findIn(s_down, 4, surf);
  int v1 = findIn(v_down, 4, vert);

  apf::Vector3 point0(0,0,0);
  apf::Vector3 point1(0,0,0);
  apf::Vector3 point2(0,0,0);
  apf::Vector3 proj(0,0,0);
  apf::Vector3 temp(0,0,0);

  int v2 = vd_vert_norm[i1][v1][0];
  int v3 = vd_vert_norm[i1][v1][1];
  assert(v2 != -1);

  m->getPoint(v_down[v1], 0, point0);
  m->getPoint(v_down[v2], 0, point1);
  m->getPoint(v_down[v3], 0, point2);

  proj = vd_projpt(point0, point1, point2);
  proj = point2 - proj;

  temp = point1 - point0;

  proj = norm_0(proj);
  temp = norm_0(temp);

  dA = norm_0(vd_cross(proj, temp))*0.5;

  apf::ModelEntity* em = m->toModel(elem_up.e[0]);
  int em_type = m->getModelType(em);
  int em_tag = m->getModelTag(em);

  dA = vd_cross(dA, point1 - point2);

  if (elem_up.n == 1) {
    // printf("Domain surface triangle.\n");
    return dA;
  }

  else if (elem_up.n == 2) {
    // printf("Internal surface triangle.\n");
    if (em_tag == geom)
      return dA;
    else
      return dA*(-1.0);
  }
}

apf::Vector3 vd_area_out(apf::Mesh2* m, apf::MeshEntity* surf) {
  apf::Downward dv;
  apf::Vector3 v(0,0,0);
  apf::Vector3 w1(0,0,0);
  apf::Vector3 w2(0,0,0);

  m->getDownward(surf, 0, dv);
  m->getPoint(dv[0], 0, v);
  m->getPoint(dv[1], 0, w1);
  m->getPoint(dv[2], 0, w2);
  return cross(w1 - v, w2 - v)/2;
}

apf::Vector3 vd_area_out_n(apf::Mesh2* m, apf::MeshEntity* surf) {
  apf::Downward dv;
  apf::Vector3 v(0,0,0);
  apf::Vector3 w1(0,0,0);
  apf::Vector3 w2(0,0,0);

  m->getDownward(surf, 0, dv);
  m->getPoint(dv[0], 0, v);
  m->getPoint(dv[1], 0, w1);
  m->getPoint(dv[2], 0, w2);

  apf::Vector3 proj(0,0,0);

  proj = vd_projpt(v, w1, w2);
  proj = w2 - proj;

  w1 = w1 - v;
  double l1 = proj.getLength();
  double l2 = w1.getLength();

  proj = norm_0(proj);
  w1 = norm_0(w1);

  return cross(w1, proj)*l1*l2/2;
}

// Given a length 3 vector of point positions, calculate the normal direction 
// of the plane defined by the points.
apf::Vector3 vd_area_out_n(apf::Mesh2* m, std::vector<apf::Vector3>* pos) {
  apf::Vector3 norm(0, 0, 0);
  if(pos->size() != 3)
    return norm;

  apf::Vector3 v(0, 0, 0);
  apf::Vector3 w1(0, 0, 0);
  apf::Vector3 w2(0, 0, 0);

  w1 = pos->at(1) - pos->at(0);
  w2 = pos->at(2) - vd_projpt(pos->at(0), pos->at(1), pos->at(2));

  double l1 = w1.getLength();
  double l2 = w2.getLength();

  w1 = norm_0(w1);
  w2 = norm_0(w2);

  return cross(w1, w2)*l1*l2/2;
}

// Given a length 3 vector of point positions, calculate the normal direction 
// of the plane defined by the points.
apf::Vector3 vd_area_out_n(apf::Vector3 p0, apf::Vector3 p1, apf::Vector3 p2) {
  apf::Vector3 norm(0, 0, 0);

  apf::Vector3 v(0, 0, 0);
  apf::Vector3 w1(0, 0, 0);
  apf::Vector3 w2(0, 0, 0);

  w1 = p1 - p0;
  w2 = p2 - vd_projpt(p0, p1, p2);

  double l1 = w1.getLength();
  double l2 = w2.getLength();

  w1 = norm_0(w1);
  w2 = norm_0(w2);

  return cross(w1, w2)*l1*l2/2;
}

// Given a set of vertices consituting a triangle, calculate the area normal. 
apf::Vector3 vd_area_tri(apf::Mesh2* m, apf::MeshEntity** v) {
  apf::Vector3 point0(0, 0, 0);
  apf::Vector3 point1(0, 0, 0);
  apf::Vector3 point2(0, 0, 0);
  m->getPoint(v[0], 0, point0);
  m->getPoint(v[1], 0, point1);
  m->getPoint(v[2], 0, point2);
  return cross(point1 - point0, point2 - point0)*0.5; 
}

// Find the normal on plane p0p1p2, perpendicular to p0p1, towards p2. 
apf::Vector3 vd_dir_in_pl(apf::Vector3 p0, apf::Vector3 p1, apf::Vector3 p2) {
  apf::Vector3 n(0,0,0);
  apf::Vector3 dir(0,0,0);
  n = norm_0(p1 - p0);
  dir = p2 - p1;
  dir = norm_0(dir - n*(n*dir));
  return dir;
}

// Given a mesh, calculate the total volume. 
double vd_tot_volm(apf::Mesh2* m) {
  double volm = 0;
  double meas;
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* elem;
  while(elem = m->iterate(it)) {
    meas = vd_volume_tet(m, elem);
    volm = volm + meas;
  }
  m->end(it);
  return volm;
}

apf::Vector3 vd_cross(apf::Vector3 v1, apf::Vector3 v2) {
  apf::Vector3 temp1(0,0,0);
  apf::Vector3 temp2(0,0,0);

  temp2 = v2 - vd_projpt(apf::Vector3(0,0,0), v1, v2);
  double mag1 = v1.getLength();
  double mag2 = temp2.getLength();
  temp1 = norm_0(v1);
  temp2 = norm_0(temp2);

  temp1 = cross(temp1, temp2);
  return temp1*mag1*mag2;
  //return (point3-point0)*cross(point1 - point0, point2 - point0)/6*mag_min;
}

// Given a tetrahedron, calculate the volume. 
double vd_volume_tet(apf::Mesh2* m, apf::MeshEntity* tet) {
  int ent_type = m->getType(tet);
  assert(ent_type == apf::Mesh::TET);

  apf::Downward down;
  m->getDownward(tet, 0, down);
  return vd_volume_tet(m, down);
}

// Given a set of vertices consituting a tetrahedron, calculate the volume. 
double vd_volume_tet(apf::Mesh2* m, apf::MeshEntity** v) {
  std::vector<apf::Vector3> points(4, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> dir(3, apf::Vector3(0,0,0));
  std::vector<double> norms(3, 0);

  m->getPoint(v[0], 0, points.at(0));
  m->getPoint(v[1], 0, points.at(1));
  m->getPoint(v[2], 0, points.at(2));
  m->getPoint(v[3], 0, points.at(3));

  dir.at(0) = points.at(3) - points.at(0);
  //dir.at(1) = points.at(1) - points.at(0);
  //dir.at(2) = points.at(2) - points.at(0);

  dir.at(1) = points.at(1) - points.at(0);
  dir.at(2) = points.at(2) - vd_projpt(points.at(0), 
                                points.at(1), points.at(2));

  double mag_min = 1;
  double mag_curr = 0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      mag_curr = std::fabs(dir.at(i)[j]);
      if(mag_curr > std::numeric_limits<double>::min() and mag_curr < mag_min)
        mag_min = std::fabs(dir.at(i)[j]);
    }
  }


  for(int i = 0; i < 3; i++) {
    dir.at(i) = dir.at(i)/mag_min;
    norms.at(i) = dir.at(i).getLength();
    dir.at(i) = norm_0(dir.at(i));
  }

  apf::Vector3 cr(0, 0, 0);
  cr = cross(dir.at(1), dir.at(2));
  double cr_norm = cr.getLength();
  cr = norm_0(cr);

  double vol = dir.at(0)*cr*cr_norm;
  for(int i = 0; i < 3; i++) {
    vol = vol*norms.at(i);
  }
  //std::cout << std::endl;
  //std::cout << v[0] << " " << point0 << " "
  //          << v[1] << " " << point1 << " "
  //          << v[2] << " " << point2 << " "
  //          << v[3] << " " << point3 << std::endl;
  mag_min = mag_min*mag_min*mag_min;
  return vol/6*mag_min;
  //return (point3-point0)*cross(point1 - point0, point2 - point0)/6*mag_min;
}

// Given three vectors v03, v01, v02 calculate the triple product. 
double vd_trip_prod(apf::Vector3 v03, apf::Vector3 v01, apf::Vector3 v02) {

  std::vector<apf::Vector3> dir(3, apf::Vector3(0,0,0));
  std::vector<double> norms(3, 0);

  dir.at(0) = v03;
  dir.at(1) = v01;
  dir.at(2) = v02 - vd_projpt(apf::Vector3(0,0,0),  v01, v02);

  double mag_min = 1;
  double mag_curr = 0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      mag_curr = std::fabs(dir.at(i)[j]);
      if(mag_curr > std::numeric_limits<double>::min() and mag_curr < mag_min)
        mag_min = std::fabs(dir.at(i)[j]);
    }
  }

  for(int i = 0; i < 3; i++) {
    dir.at(i) = dir.at(i)/mag_min;
    norms.at(i) = dir.at(i).getLength();
    dir.at(i) = norm_0(dir.at(i));
  }

  apf::Vector3 cr(0, 0, 0);
  cr = cross(dir.at(1), dir.at(2));
  double cr_norm = cr.getLength();
  cr = norm_0(cr);

  double vol = dir.at(0)*cr*cr_norm;
  for(int i = 0; i < 3; i++) {
    vol = vol*norms.at(i);
  }
  //std::cout << std::endl;
  //std::cout << v[0] << " " << point0 << " "
  //          << v[1] << " " << point1 << " "
  //          << v[2] << " " << point2 << " "
  //          << v[3] << " " << point3 << std::endl;
  mag_min = mag_min*mag_min*mag_min;
  return vol/6*mag_min;
  //return (point3-point0)*cross(point1 - point0, point2 - point0)/6*mag_min;
}

// Given a set of vertices consituting a tetrahedron, calculate the volume. 
double vd_volume_tet(std::vector<apf::Vector3>* points) {
  assert(points->size() == 4);
  std::vector<apf::Vector3> dir(3, apf::Vector3(0,0,0));
  std::vector<double> norms(3, 0);

  dir.at(0) = points->at(3) - points->at(0);
  //dir.at(1) = points->at(1) - points->at(0);
  //dir.at(2) = points->at(2) - points->at(0);

  dir.at(1) = points->at(1) - points->at(0);
  dir.at(2) = points->at(2) - vd_projpt(points->at(0), 
                                points->at(1), points->at(2));
/*
  double mag_min = 1;
  double mag_curr = 0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      mag_curr = std::fabs(dir.at(i)[j]);
      if(mag_curr > std::numeric_limits<double>::min() and mag_curr < mag_min)
        mag_min = mag_curr;
    }
  }
*/
  double mag_min = -1;
  double mag_curr = 0;
  for(int i = 0; i < 3; i++) {
    mag_curr = dir.at(i).getLength();
    if(mag_min < - std::numeric_limits<double>::min() or
       (mag_curr > std::numeric_limits<double>::min() and mag_curr < mag_min))
      mag_min = mag_curr;
  }
  if(mag_min < std::numeric_limits<double>::min())
    return 0;

  for(int i = 0; i < 3; i++) {
    dir.at(i) = dir.at(i)/mag_min;
    norms.at(i) = dir.at(i).getLength();
    dir.at(i) = norm_0(dir.at(i));
  }

  apf::Vector3 cr(0, 0, 0);
  cr = cross(dir.at(1), dir.at(2));
  double cr_norm = cr.getLength();
  cr = norm_0(cr);

  double vol = dir.at(0)*cr*cr_norm;
  for(int i = 0; i < 3; i++) {
    vol = vol*norms.at(i);
  }
  //std::cout << std::endl;
  //std::cout << v[0] << " " << point0 << " "
  //          << v[1] << " " << point1 << " "
  //          << v[2] << " " << point2 << " "
  //          << v[3] << " " << point3 << std::endl;
  mag_min = mag_min*mag_min*mag_min;
  assert(!std::isnan(vol) or !std::isnan(mag_min));
  return vol/6*mag_min;
  //return (point3-point0)*cross(point1 - point0, point2 - point0)/6*mag_min;
}

// Given a tetrahedron, return true if the volume is positive. 
bool vd_volume_tet_sign(apf::Mesh2* m, apf::MeshEntity* tet) {
  int ent_type = m->getType(tet);
  assert(ent_type == apf::Mesh::TET);

  apf::Downward down;
  m->getDownward(tet, 0, down);
  return vd_volume_tet_sign(m, down);
}

// Given a set of vertices consituting a tetrahedron, return true if the vertex 
// ordering yields a positive volume.
bool vd_volume_tet_sign(apf::Mesh2* m, apf::MeshEntity** v) {
  std::vector<apf::Vector3> points(4, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> dir(3, apf::Vector3(0,0,0));
  std::vector<double> norms(3, 0);

  m->getPoint(v[0], 0, points.at(0));
  m->getPoint(v[1], 0, points.at(1));
  m->getPoint(v[2], 0, points.at(2));
  m->getPoint(v[3], 0, points.at(3));

  //std::cout << "p0 = [" << points.at(0)[0] << "," << points.at(0)[1] << "," 
  //          << points.at(0)[2] << "];" << std::endl;
  //std::cout << "p1 = [" << points.at(1)[0] << "," << points.at(1)[1] << "," 
  //          << points.at(1)[2] << "];" << std::endl;
  //std::cout << "p2 = [" << points.at(2)[0] << "," << points.at(2)[1] << "," 
  //          << points.at(2)[2] << "];" << std::endl;
  //std::cout << "p3 = [" << points.at(3)[0] << "," << points.at(3)[1] << "," 
  //          << points.at(3)[2] << "];" << std::endl;

  dir.at(0) = points.at(3) - points.at(0);
  //dir.at(1) = points.at(1) - points.at(0);
  //dir.at(2) = points.at(2) - points.at(0);

  dir.at(1) = points.at(1) - points.at(0);
  dir.at(2) = points.at(2) - vd_projpt(points.at(0), 
                                points.at(1), points.at(2));

  double mag_min = 1;
  double mag_curr = 0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      mag_curr = std::fabs(dir.at(i)[j]);
      if(mag_curr > std::numeric_limits<double>::min() and mag_curr < mag_min)
        mag_min = std::fabs(dir.at(i)[j]);
    }
  }

  for(int i = 0; i < 3; i++) {
    dir.at(i) = dir.at(i)/mag_min;
    norms.at(i) = dir.at(i).getLength();
    dir.at(i) = norm_0(dir.at(i));
  }

  apf::Vector3 cr(0, 0, 0);
  cr = cross(dir.at(1), dir.at(2));
  double cr_norm = cr.getLength();
  cr = norm_0(cr);

  double vol = dir.at(0)*cr*cr_norm;
  for(int i = 0; i < 3; i++) {
    vol = vol*norms.at(i);
  }
  //std::cout << std::endl;
  //std::cout << v[0] << " " << point0 << " "
  //          << v[1] << " " << point1 << " "
  //          << v[2] << " " << point2 << " "
  //          << v[3] << " " << point3 << std::endl;
  mag_min = mag_min*mag_min*mag_min;

  //std::cout << "dir0 = [" << dir.at(0)[0] << "," << dir.at(0)[1] << "," 
  //          << dir.at(0)[2] << "];" << std::endl;
  //std::cout << "dir1 = [" << dir.at(1)[0] << "," << dir.at(1)[1] << "," 
  //          << dir.at(1)[2] << "];" << std::endl;
  //std::cout << "dir2 = [" << dir.at(2)[0] << "," << dir.at(2)[1] << "," 
  //          << dir.at(2)[2] << "];" << std::endl;

  return vol/6*mag_min > std::numeric_limits<double>::min();

  //return norm_0(point3-point0)*
  //            cross(norm_0(point1 - point0), norm_0(point2 - point0))/6 
  //            > std::numeric_limits<double>::min();
}


void print_tet(apf::Mesh2* m, apf::MeshEntity* tet) {
  apf::Downward down;
  apf::Vector3 pos(0, 0, 0);

  apf::ModelEntity* em = m->toModel(tet);
  int em_type = m->getModelType(em);
  int em_tag = m->getModelTag(em);

  std::cout << "Tet " << tet << " "
            << em_type << "c" << em_tag << std::endl;

  m->getDownward(tet, 0, down);
  std::cout << "Pos " << getLinearCentroid(m, tet) << " volm "
            << vd_volume_tet(m, down)
            << std::endl;

  for(int i = 0; i < 4; i++) {
    em = m->toModel(down[i]);
    em_type = m->getModelType(em);
    em_tag = m->getModelTag(em);

    std::cout << "v " << down[i] << " "
              << em_type << "c" << em_tag;
    m->getPoint(down[i], 0, pos);
    std::cout << " pos(" << i+1 << ",:) = [" 
              << pos[0] << ", " << pos[1] << ", " << pos[2]
              << "];"
              << std::endl;
  }
}

// Print the topological and geometric information of the tets around an ent.
void print_tets_ent(apf::Mesh2* m, apf::MeshEntity* ent) {

  int ent_type = m->getType(ent);
  int d = m->typeDimension[ent_type];

  if(d == 3) {
    print_tet(m, ent);
    return;
  }

  std::vector< std::vector<apf::MeshEntity*> > ents
                                (3-d, std::vector<apf::MeshEntity*>(0));

  vd_set_up(m, ent, &ents.at(0));

  for(int i = 0; i < 3-d-1; i++) {
    vd_set_up(m, &ents.at(i), &ents.at(i+1));
  }
  for(int i = 0; i < ents.at(2-d).size(); i++) {
    print_tet(m, ents.at(2-d).at(i));
  }
}

// Print the linear centroid of the given element.
apf::Vector3 vd_print_pos(apf::Mesh2* m, apf::MeshEntity* ent) {
  apf::Vector3 e_pos(0,0,0);
  int ent_type = m->getType(ent);
  if(ent_type == 0) {
    m->getPoint(ent, 0, e_pos);
  }
  else
    e_pos = getLinearCentroid(m, ent);

  std::cout << m->typeDimension[ent_type] << "-ent " << ent 
            << " " << e_pos << std::endl;
  return e_pos;
}

// Print the linear centroid of the given element.
apf::Vector3 vd_get_pos(apf::Mesh2* m, apf::MeshEntity* ent) {
  apf::Vector3 e_pos(0,0,0);
  int ent_type = m->getType(ent);
  if(ent_type == 0) {
    m->getPoint(ent, 0, e_pos);
  }
  else
    e_pos = getLinearCentroid(m, ent);

  return e_pos;
}

// Print the position of the other vertex of the edge.
apf::Vector3 vd_edge_pos_other(apf::Mesh2* m, apf::MeshEntity* edge, 
                               apf::MeshEntity* vert) {
  apf::Downward down;
  apf::Vector3 pos(0,0,0);
  m->getDownward(edge, 0, down);
  int i1 = findIn(down, 2, vert);
  assert(i1 > -1);
  i1 = (i1 + 1)%2;
  m->getPoint(down[i1], 0, pos);
  return pos;
}

apf::Vector3 get_pos_v_x_t_in_list(apf::Mesh2* m, apf::MeshEntity* tri_curr, std::vector<apf::MeshEntity*> & es_ent) {
  apf::Up up;
  apf::Vector3 pos_v(0,0,0);
  apf::Downward d_t;
  apf::Downward d_v;

  m->getUp(tri_curr, up);
  m->getDownward(up.e[0], 2, d_t);
  int t1 = findIn(d_t, 4, tri_curr);
  int v1 = lookup_tet_xsurf[t1];
  m->getDownward(up.e[0], 0, d_v);
  int e1 = findIn(&es_ent, es_ent.size(), d_v[v1]);
  if(e1 > -1) {
    pos_v = vd_get_pos(m, d_v[v1]);
  }
  else {
    m->getDownward(up.e[1], 2, d_t);
    t1 = findIn(d_t, 4, tri_curr);
    v1 = lookup_tet_xsurf[t1];
    m->getDownward(up.e[1], 0, d_v);
    e1 = findIn(&es_ent, es_ent.size(), d_v[v1]);
    assert(e1 > -1);
    pos_v = vd_get_pos(m, d_v[v1]);
  }
  return pos_v;
}

double vd_avg_len(apf::Mesh2* m, apf::MeshEntity* vertex) {

  apf::Up edges;

  m->getUp(vertex, edges);

  double avg = 0;

  for(int i = 0; i < edges.n; i++) {
    apf::MeshElement* ee = createMeshElement(m, edges.e[i]);
    avg = avg + measure(ee);
    destroyMeshElement(ee);
  }
  return avg/edges.n;

}

apf::Vector3 get_edge_dir(apf::Mesh2* m, apf::MeshEntity* edge) {
  apf::Downward down;
  apf::Vector3 p1(0,0,0);
  apf::Vector3 dir(0,0,0);

  m->getDownward(edge, 0, down);

  m->getPoint(down[0], 0, p1);
  m->getPoint(down[1], 0, dir);

  dir = dir - p1;

  dir = norm_0(dir);
  return dir;
}

double vd_avg_meas(apf::Mesh2* m, Entity_set* e_set) {

  double avg = 0;

  for(int i = 0; i < e_set->n; i++) {
    apf::MeshElement* ee = createMeshElement(m, e_set->e[i]);
    avg = avg + measure(ee);
    destroyMeshElement(ee);
  }

  return avg/e_set->n;

}

double vd_minmax_meas(apf::Mesh2* m, Entity_set* e_set, bool min) {

  int c = 1;
  if (!min)
    c = -1;

  if(e_set->n == 0)
    return 0;

  apf::MeshElement* ee = createMeshElement(m, e_set->e[0]);
  double len_cross = c*measure(ee);
  destroyMeshElement(ee);
  for (int i = 1; i < e_set->n; i++) {
    ee = createMeshElement(m, e_set->e[i]);
    double len1 = c*measure(ee);
    destroyMeshElement(ee);
    if (len1 < len_cross)
      len_cross = len1;
  }
  len_cross = len_cross*c;

  return len_cross;

}



// Print the normals for elements inside the mesh, using area out.
void vd_print_norm2(apf::Mesh2* m) {

  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* elem;
  while(elem = m->iterate(it)) {
    printf("%p\n",(void*)elem);
    apf::ModelEntity* em = m->toModel(elem);
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);

    apf::Downward s_down;
    m->getDownward(elem, 2, s_down);
    apf::Vector3 dA(0, 0, 0);
    printf("For the element %p the outward pointing normals are:\n", 
                                                                (void*)elem);
    for (int i = 0; i < 4; i++) {
      dA = vd_area_out(m, s_down[i], em_tag);
      std::cout << s_down[i]<< " has the outward normal "
      << norm_0(dA) <<std::endl;
    }
  }
  m->end(it);

}

// Given two vectors, calculate the inner angle.
double vd_inner_angle(apf::Vector3 v1, apf::Vector3 v2) {
  v1 = norm_0(v1);
  v2 = norm_0(v2);

  return acos(std::min(std::max(v1*v2,-1.0),1.0));
  //return std::fmod(acos(std::min(std::max(v1*v2,-1.0),1.0)) + PI_L, 2*PI_L);
}

// Given two adjacent boundary surfaces and the geometry of the inside volume, 
// find the exterior angle defined between those surfaces.
double vd_ext_angle(apf::Mesh2* m, apf::MeshEntity* surf0, 
                                  apf::MeshEntity* surf1, int geom) {

  double angle = 0;
  apf::MeshEntity* edge = vd_find_tris_edge(m, surf0, surf1);
  apf::Vector3 l_dir(0,0,0);
  apf::Vector3 l_pos(0,0,0);
  apf::Vector3 pos0(0,0,0);
  apf::Vector3 pos1(0,0,0);
  l_dir = get_edge_dir(m, edge);
  l_pos = vd_get_pos(m, edge);
  pos0 = vd_get_pos(m, surf0);
  pos1 = vd_get_pos(m, surf1);

  apf::Up elem_up;
  m->getUp(surf0, elem_up);
  assert((elem_up.n == 1) || (elem_up.n == 2));

  apf::ModelEntity* em = m->toModel(elem_up.e[0]);
  int em_tag = m->getModelTag(em);
  if(em_tag == geom)
    angle = vd_ext_angle(pos0, pos1, l_dir, l_pos, vd_get_pos(m, elem_up.e[0]));
  else {
    assert(elem_up.n == 2);
    em = m->toModel(elem_up.e[1]);
    em_tag = m->getModelTag(em);
    angle = vd_ext_angle(pos0, pos1, l_dir, l_pos, vd_get_pos(m, elem_up.e[1]));
  }
  return angle;
}


// Given the outward normal and a point on each planes bounding the same volume, 
// find the exterior angle. Assume all points are coplanar and 
// the planes intersect at zero point.
double vd_ext_angle_n(apf::Vector3 pos1, apf::Vector3 pos2,
                    apf::Vector3 norm1, apf::Vector3 norm2) {
  // Range of acos 0 to PI
  double ang_inner = std::acos(std::min(std::max(norm1*norm2,-1.0),1.0));
  // The triangles are coplanar:
  if(std::fabs(ang_inner - PI_L) < std::numeric_limits<double>::epsilon())
    return 0.;
  // Check sides of points with respect to the other plane:
  double dir1 = norm1*pos2;

  // The interior is within the smaller arc:
  if(dir1 < - std::numeric_limits<double>::min())
    return ang_inner;
  else
    return -ang_inner;
}

// Given two planes and points on each plane, the direction of the line that
// joins the planes, a center point on one side of the planes, find the exterior
// angle to the angle that lies on the side of the planes with the center point.
double vd_ext_angle(apf::Vector3 pos1, apf::Vector3 pos2, apf::Vector3 l_dir, 
                                apf::Vector3 l_pos, apf::Vector3 ctr) {

  apf::Vector3 norm0(0, 0, 0);
  norm0 = norm_0(norm_0(l_dir));
  // TODO apf::Vector3::normalize causes float point accuracy issues and gives
  // non-zero exterior angle for some coplanar triangles. Double normalization
  // seems to mitigate this issue, but it shouldn't be necessary...
  pos1 = pos1 - l_pos;
  pos2 = pos2 - l_pos;
  ctr = ctr - l_pos;

  pos1 = pos1 - norm0*(pos1*norm0);
  pos2 = pos2 - norm0*(pos2*norm0);
  ctr = ctr - norm0*(ctr*norm0);

  pos1 = norm_0(norm_0(pos1));
  pos2 = norm_0(norm_0(pos2));
  ctr = norm_0(norm_0(ctr));

  double ang_inner = std::acos(std::min(std::max(pos1*pos2,-1.0),1.0));
  double ang_inner2 = std::acos(std::min(std::max(pos2*pos1,-1.0),1.0));

  std::cout << ang_inner << " " << ang_inner2 << std::endl;
  // The triangles are coplanar:
  if(std::fabs(ang_inner - PI_L) < std::numeric_limits<double>::epsilon())
    return 0.;

  // If the angle between the inner position and one of the ends is larger than,  
  // the inner angle, the exterior angle is negative:
  double ang_ctr = std::acos(std::min(std::max(pos1*ctr,-1.0),1.0));
  if((ang_inner - ang_ctr) < std::numeric_limits<double>::epsilon()
     or std::fabs(ang_ctr - PI_L) < std::numeric_limits<double>::epsilon() )
    return ang_inner - PI_L;

  // Check the order of end points:
  double dir1 = vd_cross(pos1, pos2)*l_dir;
  double dir2 = vd_cross(pos1, ctr)*l_dir;

  // The ctr vertex is within the inner arc.
  if(dir1*dir2 > std::numeric_limits<double>::min())
    return ang_inner;
  // The ctr vertex is within the outer arc.
  else
    return ang_inner - PI_L;
}

// TODO this doesn't look like a good approach. 
// Given a vertex belonging to a 1cell, find two triangles bounding the 3-cell.
// Return the exterior angle. The triangles found in this fashion do not 
// necessarily share an edge. Also if 3cell has disjoint sets around v, it can
// also find triangles that do not bound the same disjoint set. 
double vd_ext_angle(apf::Mesh2* m, apf::MeshEntity* v, int cell) {
/*
  apf::ModelEntity* mdl = m->toModel(v);
  apf::ModelEntity* mdl_3c = m->findModelEntity(3, cell);
  int type = m->getModelType(mdl);
  int tag = m->getModelTag(mdl);

  assert(type == 1);

  std::vector<apf::MeshEntity*> es_e(0);
  std::vector<apf::MeshEntity*> es_e_cp(0);

  vd_set_up(m, v, &es_e);

  es_e_cp.reserve(es_e.size());
  for (int i = 0; i < es_e.size(); i++) {
    if(m->toModel(es_e.at(i)) == mdl)
      es_e_cp.push_back(es_e.at(i));
  }

  std::vector<apf::MeshEntity*> es_s(0);

  vd_set_up(m, &es_e_cp, &es_s);

  bool found1 = false;
  bool found2 = false;
  apf::MeshEntity* s1;
  apf::MeshEntity* s2;
  apf::Up up;
  for (int i = 0; i < es_s.size(); i++) {
    int type = m->getModelType(m->toModel(es_s.at(i)) );
    if(type == 2) {
      m->getUp(es_s.at(i), up);
      for (int j = 0; j < up.n; j++) {
        if(m->toModel(up.e[j]) == mdl_3c) {
          if(found1) {
            found2 = true;
            s2 = es_s.at(i);
            i = es_s.size();
          }
          else {
            found1 = true;
            s1 = es_s.at(i);
          }
          j = up.n;
        }
      }
    }
  }

  assert(found1 and found2);

  return vd_ext_angle(m, s1, s2, cell);
*/
  return 0;
}


// Given a vertex on a 2cell, and an interior 3cell, find the mean curvature
// of the surface associated with the normal outward from the 3cell.
double vd_mean_curv(apf::Mesh2* m, apf::MeshEntity* v, int cell) {

  apf::ModelEntity* mdl = m->toModel(v);
  int type = m->getModelType(mdl);
  int tag = m->getModelTag(mdl);

  if(type != 2)
    return 0;

  std::vector<apf::MeshEntity*> es_e(0);
  std::vector<apf::MeshEntity*> es_s(0);
  std::vector<apf::MeshEntity*> es_s_cp(0);

  vd_set_up(m, v, &es_e);
  vd_set_up(m, &es_e, &es_s);

  for (int i = 0; i < es_s.size(); i++) {
    if(m->toModel(es_s.at(i)) == mdl)
      es_s_cp.push_back(es_s.at(i));
  }

  apf::Downward d_v;

  apf::Vector3 neg_grad(0, 0, 0);
  for (int i = 0; i < es_s_cp.size(); i++) {
/*
    m->getDownward(es_s_cp.at(i), 0, d_v);
    int v1 = findIn(d_v, 3, v);
    m->getPoint(d_v[v1], 0, temp);
    m->getPoint(d_v[lookup_triverts[v1][0]], 0, w1);
    m->getPoint(d_v[lookup_triverts[v1][1]], 0, w2);

    area_out = vd_area_out(m, es_s_cp.at(i), cell);
    // The vertex ordering for the positively oriented surface might be different
    // than the current ordering.
    if(area_out*vd_cross(w1-temp, w2-temp) < -std::numeric_limits<double>::min())
      area_out = area_out*(-1);
    neg_grad = neg_grad + cross(area_out, w1-w2);
    //apf::Vector3 dist = (w1+w2)/2;
    //neg_grad = neg_grad + cross(area_out, (w1-w2)/dist.getLength());
*/
    neg_grad = neg_grad + vd_grad(m, v, es_s_cp.at(i), cell);
  }

  //return -neg_grad.getLength()/2;
  return -neg_grad.getLength();
}

// Given a vertex and an interior 3stratum, find the Gaussian curvature
// of the surfaces of the 3cell joining at the vertex by angular defect formula.
// If there are disjoint components of the 3stratum, return the sum of the angular
// defects due to each disjoint component. 
// This function lets calculation of angular defect at a vertex, but considering
// that information of repeated neighborhood is needed for the invariant 
// calculation, this is not the main function used to extract Gaussian curvature.
// See topo_feat.h. 
double vd_gauss_curv(apf::Mesh2* m, apf::MeshEntity* v, int cell) {
}

// Given a set of triangles joined at a vertex, find the negative of the 
// gradient of the total area, wrt. the vertex position.
apf::Vector3 vd_neg_grad(apf::Mesh2* m, Entity_set* tri, apf::MeshEntity* vert) {
  apf::Vector3 neg_grad(0, 0, 0);
  double vect[3] = {0,0,0};
  neg_grad.fromArray(vect);
  apf::Downward dv;

  apf::Vector3 w1(0, 0, 0);
  apf::Vector3 w2(0, 0, 0);
  apf::Vector3 v(0, 0, 0);
  apf::Vector3 area_out(0, 0, 0);

  m->getPoint(vert, 0, v);

  for (int i = 0; i < tri->n; i++) {

    m->getDownward(tri->e[i], 0, dv);
    int e1 = findIn(dv, 3, vert);
    m->getPoint(dv[lookup_triverts[e1][0]], 0, w1);
    m->getPoint(dv[lookup_triverts[e1][1]], 0, w2);
    area_out = norm_0(cross(w1 - v, w2 - v));
    neg_grad = neg_grad + cross(area_out, w1-w2);
    //apf::Vector3 dist = (w1+w2)/2;
    //neg_grad = neg_grad + cross(area_out, (w1-w2)/dist.getLength());
  }
  return neg_grad;
}

// Given a distance field and a tet on which the field is defined, 
// approximate the distance defined on one the bounding vertices, approximating
// the gradient of the distance field by using the distances on the other 
// vertices and assuming a linear variation of the distance across the tet. 
double approx_dist(apf::Mesh2* m, std::map<apf::MeshEntity*, double> & dist, apf::MeshEntity* tet, apf::MeshEntity* vert) {
  apf::Downward d_v;
  m->getDownward(tet, 0, d_v);
  return approx_dist(m, dist, d_v, vert);
}

// Given a distance field and a tet on which the field is defined, 
// approximate the distance defined on one the bounding vertices, approximating
// the gradient of the distance field by using the distances on the other 
// vertices and assuming a linear variation of the distance across the tet. 
// TODO for this to work, all distances should be from the same boundary!
// 1c and 2c adjacencies should be considered
double approx_dist(apf::Mesh2* m, std::map<apf::MeshEntity*, double> & dist, apf::Downward d_v, apf::MeshEntity* vert) {
  std::vector<apf::Vector3> pos(4, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> n(3, apf::Vector3(0,0,0));
  apf::Vector3 p(0,0,0);
  apf::Vector3 n_g(0,0,0);

  int i0 = findIn(d_v, 4, vert);
  for(int i = 0; i < 4; i++) {
    m->getPoint(d_v[i], 0, pos.at(i));
  }
  int ids[3] = {1, 2, 3};

  int i1 = (i0+1) %4;
  int i2 = (i0+2) %4;
  int i3 = (i0+3) %4;

  // Replace pos.at(i1) or pos.at(i1) with p to obtain an right angled triangle
  p = vd_projpt(pos.at(i1), pos.at(i2), pos.at(i3));
  n.at(0) = pos.at(i3) - p;
  n.at(1) = pos.at(i1) - p;
  n.at(2) = pos.at(i0) - p;

  double a = n.at(0).getLength();
  double b = n.at(1).getLength();
  double c = n.at(2).getLength();

  n.at(0) = norm_0(n.at(0));
  n.at(1) = norm_0(n.at(1));
  n.at(2) = norm_0(n.at(2));

  double dist_p = dist[d_v[i1]] + (p - pos.at(i1)).getLength()
     /(pos.at(i2) - pos.at(i1)).getLength()*(dist[d_v[i2]] - dist[d_v[i1]]);

  n_g[0] = (dist[d_v[i3]] - dist_p)/a;
  n_g[1] = (dist[d_v[i1]] - dist_p)/b;
  assert(std::fabs(n_g[0]) > -std::numeric_limits<double>::min());
  assert(std::fabs(n_g[0]) < 1);
  assert(std::fabs(n_g[1]) > -std::numeric_limits<double>::min());
  assert(std::fabs(n_g[1]) < 1);

  n_g[2] = std::sqrt(1 - n_g[0]*n_g[0] - n_g[1]*n_g[1]);
  n_g = n.at(0)*n_g[0] + n.at(1)*n_g[1] + n.at(2)*n_g[2];

  return dist_p + (pos.at(i0) - p)*n_g;
}


// Given two vertices A and B that form a vector V from the second one to the 
// first one, and two other vertices C and D, check if the area vector obtained  
// by taking the cross products of vectors V_C and V_D from the center point of 
// A and B to the other two vertices is positively aligned with the vector V.
bool vd_vec_align(apf::Mesh2* m, apf::MeshEntity** vert) {

  //std::cout << "Checking alignment" << std::endl;
  apf::Vector3 VA(0, 0, 0);
  apf::Vector3 VB(0, 0, 0);
  apf::Vector3 VC(0, 0, 0);
  apf::Vector3 VD(0, 0, 0);

  apf::Vector3 ctr(0, 0, 0);

  m->getPoint(vert[0], 0, VA);
  m->getPoint(vert[1], 0, VB);
  m->getPoint(vert[2], 0, VC);
  m->getPoint(vert[3], 0, VD);

  ctr = (VA + VB)/2;

  if (cross(VC - ctr, VD - ctr)*(VA-VB) > 0) {
    //std::cout << "Aligned" << std::endl;
    return true;
  }
  else {
    //std::cout << "Not aligned" << std::endl;
    return false;
  }
}

// Given four triangles that can form a tetrahedron, reorder them so that the 
// tetrahedron will have a positive volume.
// Simply check if the volume with the given triangle order is positive. If not,
// switch the positions of two of the triangles.
//bool vd_vec_align(apf::Mesh2* m, apf::MeshEntity** tri) {
//}

// Given a triangle and a vertex, calculate the angle at the vertex corner in 
// radians.
double vd_int_angle(apf::Mesh2* m, apf::MeshEntity* tri, apf::MeshEntity* v) {

  apf::Downward down;
  apf::Vector3 vec0(0, 0, 0);
  apf::Vector3 vec1(0, 0, 0);
  apf::Vector3 vec2(0, 0, 0);

  m->getDownward(tri, 0, down);
  int v1 = findIn(down, 3, v);

  m->getPoint(v, 0, vec0);
  m->getPoint(down[lookup_triverts[v1][0]], 0, vec1);
  m->getPoint(down[lookup_triverts[v1][1]], 0, vec2);

  vec1 = vec1 - vec0;
  vec2 = vec2 - vec0;
  vec1 = vec1/sqrt(vec1*vec1);
  vec2 = vec2/sqrt(vec2*vec2);

  return std::acos(vec1*vec2);

}

double vd_int_angle(apf::Mesh2* m, apf::MeshEntity* e_1, apf::MeshEntity* e_2,
                                                           apf::MeshEntity* v) {
  apf::Downward down;
  apf::Vector3 vec0(0, 0, 0);
  apf::Vector3 vec1(0, 0, 0);
  apf::Vector3 vec2(0, 0, 0);

  m->getDownward(e_1, 0, down);
  int v1 = findIn(down, 2, v);

  m->getPoint(v, 0, vec0);
  m->getPoint(down[lookup_edgeverts[v1]], 0, vec1);

  m->getDownward(e_2, 0, down);
  v1 = findIn(down, 2, v);
  m->getPoint(down[lookup_edgeverts[v1]], 0, vec2);

  vec1 = vec1 - vec0;
  vec2 = vec2 - vec0;
  vec1 = vec1/sqrt(vec1*vec1);
  vec2 = vec2/sqrt(vec2*vec2);

  return std::acos(vec1*vec2);

}

// Given the vertex and tetrahedra, find the distance to the triangle across.
double vd_find_dist(apf::Mesh2* m, apf::MeshEntity* tet, apf::MeshEntity* v) {

  apf::MeshElement* ee = createMeshElement(m, tet);
  double V_curr = measure(ee);
  destroyMeshElement(ee);

  apf::Downward d_v;
  apf::Downward d_s;
  m->getDownward(tet, 0, d_v);
  m->getDownward(tet, 2, d_s);

  // Calculate the maximum displacement possible before inverting.
  int v_id = findIn(d_v, 4, v);
  ee = createMeshElement(m,d_s[lookup_tet_surfx[v_id]]);

  double A_curr = measure(ee);
  destroyMeshElement(ee);

  return V_curr/A_curr*3;
}

// Skip given tetrahedra. Used mainly in collapse preconditioning, where there
// are inverting elements around vertex. In that case, do not consider the 
// distances to triangles within inverted tetrahedra. If none found, return
// -1.
double vd_dist_v_x(apf::Mesh2* m, apf::MeshEntity* v, apf::Vector3 dir, 
            std::map<apf::MeshEntity*, bool> &skip, apf::Vector3* int_pt) {
  if(dir.getLength() < std::numeric_limits<double>::min())
    return -1;

  apf::Vector3 o(0,0,0);
  apf::Vector3 proj(0,0,0);
  std::vector<apf::MeshEntity*> ee(0);
  std::vector<apf::MeshEntity*> es(0);
  std::vector<apf::MeshEntity*> et(0);

  m->getPoint(v, 0, o);

  vd_set_up(m, v, &ee);
  vd_set_up(m, &ee, &es);
  vd_set_up(m, &es, &et);

  apf::Downward d_v;
  apf::Downward d_t;

  bool hit = false;
  for(int i = 0; i < et.size(); i++) {
    if(!skip[et.at(i)] ) {
      m->getDownward(et.at(i), 0, d_v);
      m->getDownward(et.at(i), 2, d_t);
      int i1 = findIn(d_v, 4, v);
      assert(i1 > -1);
      apf::MeshEntity* tri_curr = d_t[lookup_tet_surfx[i1]];
      if(tri_int_line(m, tri_curr, o, dir, &proj)) {
        hit = true;
        if(int_pt != NULL)
          *int_pt = proj;
        return (proj - o).getLength();
      }
    }
  }
  return -1;
}

// Given a vertex v and a direction, find the distance to the nearest triangle 
// along the direction. Assume direction vector is non-zero and all tetrahedra 
// are positive volume. If direction is zero, return -1.
double vd_dist_v_x(apf::Mesh2* m, apf::MeshEntity* v, apf::Vector3 dir, 
                                                    apf::Vector3* int_pt) {
  if(dir.getLength() < std::numeric_limits<double>::min())
    return -1;

  apf::Vector3 o(0,0,0);
  apf::Vector3 proj(0,0,0);
  std::vector<apf::MeshEntity*> ee(0);
  std::vector<apf::MeshEntity*> es(0);
  std::vector<apf::MeshEntity*> es_x(0);
  std::vector<apf::MeshEntity*> et(0);

  m->getPoint(v, 0, o);

  vd_set_up(m, v, &ee);
  vd_set_up(m, &ee, &es);
  vd_set_up(m, &es, &et);
  vd_set_down(m, &et, &es_x);
  vd_remove_set(&es_x, &es);

  bool hit = false;
  for(int i = 0; i < es_x.size(); i++) {
    if(tri_int_line(m, es_x.at(i), o, dir, &proj)) {
      hit = true;
      if(int_pt != NULL)
        *int_pt = proj;
      return (proj - o).getLength();
    }
  }
  assert(hit);
}

// Given a vertex v find the displacement to the nearest edge 
// along the plane perpendicular to the edge.
apf::Vector3 vd_dist_v_ex(apf::Mesh2* m, apf::MeshEntity* v) {
  apf::Vector3 o(0,0,0);
  apf::Vector3 e_pos(0,0,0);
  apf::Vector3 dir(0, 0, 0);
  apf::Vector3 disp(0, 0, 0);
  double dist = -1;
  apf::Downward d_e;
  apf::Downward d_v;

  std::vector<apf::MeshEntity*> ee(0);
  std::vector<apf::MeshEntity*> es(0);

  m->getPoint(v, 0, o);

  vd_set_up(m, v, &ee);
  vd_set_up(m, &ee, &es);

  for(int i = 0; i < es.size(); i++) {
    m->getDownward(es.at(i), 0, d_v);
    m->getDownward(es.at(i), 1, d_e);
    int v_id = findIn(d_v, 3, v);
    assert(v_id > -1);

    apf::MeshEntity* e_x = d_e[lookup_tri_edgex[v_id]];
    m->getPoint(d_v[(v_id+1)%3], 0, e_pos);
    m->getPoint(d_v[(v_id+2)%3], 0, dir);
    dir = norm_0(dir - e_pos);
    e_pos = getLinearCentroid(m, e_x);
    e_pos = e_pos - o;
    e_pos = e_pos - dir*(dir*e_pos);
    double dist_curr = e_pos.getLength();
    if(dist < std::numeric_limits<double>::min()) {
      dist = dist_curr;
      disp = e_pos;
    }
  }
  assert(dist > - std::numeric_limits<double>::min());
  return disp;
}

// Given a vertex v find the displacement to the nearest triangle 
// along the surface normal. Assume tetrahedra are positive volume.
apf::Vector3 vd_dist_v_x(apf::Mesh2* m, apf::MeshEntity* v) {
  apf::Vector3 o(0,0,0);
  apf::Vector3 a_pos(0,0,0);
  apf::Vector3 area(0, 0, 0);
  apf::Vector3 disp(0, 0, 0);
  double dist = -1;

  std::vector<apf::MeshEntity*> ee(0);
  std::vector<apf::MeshEntity*> es(0);
  std::vector<apf::MeshEntity*> es_x(0);
  std::vector<apf::MeshEntity*> et(0);

  m->getPoint(v, 0, o);

  vd_set_up(m, v, &ee);
  vd_set_up(m, &ee, &es);
  vd_set_up(m, &es, &et);
  vd_set_down(m, &et, &es_x);
  vd_remove_set(&es_x, &es);

  for(int i = 0; i < es_x.size(); i++) {
    area = vd_area_out(m, es_x.at(i));
    a_pos = getLinearCentroid(m, es_x.at(i));
    a_pos = a_pos - o;
    area = norm_0(area);
    double dist_curr = area*a_pos;
    // If the area normal intersects outside the triangle, worst case scenario
    // there is a triangle closer than this. If farther, this is conservative.
    // TODO still not the best implementation.
    if(dist_curr < - std::numeric_limits<double>::min()) {
      area = area*(-1);
      dist_curr = -dist_curr;
    }
    if(dist < std::numeric_limits<double>::min()) {
      dist = dist_curr;
      disp = area*dist;
    }
  }
  assert(dist > - std::numeric_limits<double>::min());
  return disp;
}

// Given a vertex v, a set of triangles across and a direction, find the 
// displacement to the nearest triangle plane along the given direction.
double vd_dist_v_x_pl_dir(apf::Mesh2* m, apf::MeshEntity* v, 
                        std::vector<apf::MeshEntity*> &es_x, apf::Vector3 dir) {

  apf::Vector3 o(0,0,0);
  apf::Vector3 a_pos(0,0,0);
  apf::Vector3 area(0, 0, 0);
  apf::Vector3 disp(0, 0, 0);
  double dist = -1;

  m->getPoint(v, 0, o);

  for(int i = 0; i < es_x.size(); i++) {
    area = vd_area_out(m, es_x.at(i));
    a_pos = getLinearCentroid(m, es_x.at(i));
    a_pos = a_pos - o;
    area = norm_0(area);
    double dist_curr = area*a_pos;
    // If the area normal intersects outside the triangle, worst case scenario
    // there is a triangle closer than this. If farther, this is conservative.
    // TODO still not the best implementation.
    if(dist_curr < - std::numeric_limits<double>::min()) {
      area = area*(-1);
      dist_curr = -dist_curr;
    }
    double a_cos = area*dir;
    if(a_cos > std::numeric_limits<double>::min()) {
      dist_curr = dist_curr/a_cos;

      if(dist < std::numeric_limits<double>::min() or dist > dist_curr) {
        dist = dist_curr;
      }
    }
  }
  assert(dist > - std::numeric_limits<double>::min());
  return dist;
}

// Given the vertex and tetrahedron, find the distance to the triangle across, 
// weighted by the surface across. So, a vector with magnitude of volume, in 
// the area normal direction.
apf::Vector3 vd_find_dist_w(apf::Mesh2* m, apf::MeshEntity* tet, 
                                                        apf::MeshEntity* v) {

  apf::MeshElement* ee = createMeshElement(m, tet);
  double V_curr = measure(ee);
  destroyMeshElement(ee);

  apf::Downward d_v;
  apf::Downward d_s;
  m->getDownward(tet, 0, d_v);
  m->getDownward(tet, 2, d_s);

  // Calculate the maximum displacement possible before inverting.
  int v_id = findIn(d_v, 4, v);

  apf::Vector3 v_pos(0, 0, 0);
  m->getPoint(v, 0, v_pos);
  apf::Vector3 a_pos(0, 0, 0);
  a_pos = getLinearCentroid(m, d_s[lookup_tet_surfx[v_id]]);

  apf::Vector3 area(0, 0, 0);
  area = vd_area_out(m, d_s[lookup_tet_surfx[v_id]]);
  v_pos = a_pos - v_pos;
  area = norm_0(area);

  v_pos = norm_0(v_pos);

  if(area*v_pos < 0)
    return area*V_curr*(-1);
  else
    return area*V_curr;

}

void vd_get_apos(apf::Mesh2* m, apf::MeshEntity* v, 
                         std::vector<apf::MeshEntity*>* tet, 
                         std::vector<apf::Vector3>* area, 
                         std::vector<apf::Vector3>* a_pos) {
  a_pos->resize(0);
  area->resize(0);

  a_pos->resize(tet->size());
  area->resize(tet->size());

  double z[3] = {0,0,0};

  apf::Vector3 v_pos(0, 0, 0);
  m->getPoint(v, 0, v_pos);
  for(int i = 0; i < tet->size(); i++) {

    apf::Downward d_v;
    apf::Downward d_s;
    m->getDownward(tet->at(i), 0, d_v);
    m->getDownward(tet->at(i), 2, d_s);

    // Calculate the maximum displacement possible before inverting.
    int v_id = findIn(d_v, 4, v);

    a_pos->at(i).fromArray(z);
    a_pos->at(i) = getLinearCentroid(m, d_s[lookup_tet_surfx[v_id]]);

    area->at(i).fromArray(z);
    area->at(i) = vd_area_out(m, d_s[lookup_tet_surfx[v_id]]);
    a_pos->at(i) = a_pos->at(i) - v_pos;

    double dist_curr = norm_0(area->at(i))*a_pos->at(i);

    if(dist_curr < 0) {
      area->at(i) = area->at(i)*(-1);
    }
    //std::cout << "relative pos " << a_pos->at(i) << " area " << area->at(i)
    //          << std::endl;
  }
}

// Find whether a point lies inside a triangle or not. Assumes the point and 
// triangle are more or less coplanar.
bool pt_on_tri(apf::Mesh2* m, apf::Vector3 v_pos, apf::MeshEntity* tri) {
  apf::Downward d_v;
  std::vector<apf::Vector3> pos(3, apf::Vector3(0,0,0));
  apf::Vector3 e_pos(0,0,0);
  apf::Vector3 e_dir(0,0,0);

  apf::Vector3 temp1(0,0,0);
  apf::Vector3 temp2(0,0,0);

  m->getDownward(tri, 0, d_v);

  for(int i = 0; i < 3; i++) {
    m->getPoint(d_v[i], 0, pos.at(i));
  }

  for(int i1 = 0; i1 < 3; i1++) {
    int i2 = (i1+1) % 3;
    int i3 = (i1+2) % 3;
    e_pos = (pos.at(i1) + pos.at(i2)) /2;
    e_dir = norm_0((pos.at(i1) - pos.at(i2)));

    temp1 = (v_pos - e_pos);
    temp1 = temp1 - e_dir*(e_dir*temp1);

    temp2 = (pos.at(i3) - e_pos);
    temp2 = temp2 - e_dir*(e_dir*temp2);

    // Does the point and other  vertex lie on different sides of the edge:
    if( temp1*temp2 < - std::numeric_limits<double>::min())
      return false;
  }

  return true;
}

// Given a position enclosed by a set of triangles, collect the triangle areas
// positions relative to the position.
// Only consider the distances to the triangles, if the intersection with the 
// plane lies inside the triangle.
std::pair<double, double> vd_get_apos(apf::Mesh2* m, apf::Vector3 v_pos, 
                         std::vector<apf::MeshEntity*>* tri, 
                         std::vector<apf::Vector3>* area, 
                         std::vector<apf::Vector3>* a_pos) {
  a_pos->resize(0);
  area->resize(0);

  a_pos->resize(tri->size());
  area->resize(tri->size());

  double dist_min = -1;
  double dist_max = -1;

  for(int i = 0; i < tri->size(); i++) {

    a_pos->at(i) = getLinearCentroid(m, tri->at(i));

    area->at(i) = vd_area_out(m, tri->at(i));
    a_pos->at(i) = a_pos->at(i) - v_pos;

    double dist_curr = norm_0(area->at(i))*a_pos->at(i);

    if(dist_curr < 0) {
      area->at(i) = area->at(i)*(-1);
    }
    apf::Vector3 n_curr(0,0,0);
    apf::Vector3 p_int(0,0,0);
    n_curr = norm_0(area->at(i));
    p_int = v_pos + n_curr*(dist_curr);

    if( pt_on_tri(m, v_pos, tri->at(i)) ) {

      dist_curr = std::fabs(dist_curr);
      if(dist_min < 0 or dist_curr < dist_min)
        dist_min = dist_curr;
      if(dist_max < 0 or dist_curr > dist_max)
        dist_max = dist_curr;
      //std::cout << "relative pos " << a_pos->at(i) << " area " << area->at(i)
      //          << std::endl;
    }
  }
  return std::make_pair(dist_min, dist_max);
}

/*

// Check whether the given triangles intersect/lie within the sphere. If not, 
// return the normal distance to the closest triangle.
std::pair<bool, double> tri_int_sphere(apf::Mesh2* m, apf::Vector3 v_pos, 
                         double r, std::vector<apf::MeshEntity*>* tri) {

  double dist_min = -1;
  double dist_max = -1;

  apf::Vector3 a_pos(0,0,0);
  apf::Vector3 area(0,0,0);

  for(int i = 0; i < tri->size(); i++) {

    a_pos = getLinearCentroid(m, tri->at(i));

    area = vd_area_out(m, tri->at(i));
    a_pos = a_pos - v_pos;

    double dist_curr = norm_0(area)*a_pos;

    if(dist_curr < 0) {
      area = area*(-1);
    }
    apf::Vector3 n_curr(0,0,0);
    apf::Vector3 p_int(0,0,0);
    n_curr = norm_0(area);
    p_int = v_pos + n_curr*(dist_curr);

    std::vector<bool> out(3,false);

    apf::Downward d_v;
    std::vector<apf::Vector3> pos(3, apf::Vector3(0,0,0));
    apf::Vector3 e_pos(0,0,0);

    m->getDownward(tri->at(i), 0, d_v);

    for(int i = 0; i < 3; i++) {
      m->getPoint(d_v[i], 0, pos.at(i));
    }

    bool found = false;
    for(int i1 = 0; i1 < 3; i1++) {
      int i2 = (i1+1) % 3;
      int i3 = (i1+2) % 3;
      e_pos = (pos.at(i1) + pos.at(i2)) /2;
      // Does the point and other  vertex lie on different sides of the edge:
      if( (p_int - e_pos) * (pos.at(i3) - e_pos) 
                          < - std::numeric_limits<double>::min()) {
        out.at(i1) = true;
        found = true;
      }
    }

    if(!found) {
    }
    
    if( pt_on_tri(m, v_pos, tri->at(i)) ) {

      dist_curr = std::fabs(dist_curr);
      if(dist_min < 0 or dist_curr < dist_min)
        dist_min = dist_curr;
      if(dist_max < 0 or dist_curr > dist_max)
        dist_max = dist_curr;
      //std::cout << "relative pos " << a_pos->at(i) << " area " << area->at(i)
      //          << std::endl;
    }
    else {
    }
  }
  return std::make_pair(dist_min, dist_max);
}
*/


// Get the radius of the largest sphere that doesn't intersect any triangles
// across the vertex and the radius of the largest sphere without intersecting
// a single triangle.
std::pair<double, double> vert_dist_sphere(apf::Mesh2* m, 
                                             apf::MeshEntity* vert) {

  apf::Vector3 v_pos(0,0,0);
  std::vector<apf::MeshEntity*> es_e_ctr(0);
  std::vector<apf::MeshEntity*> es_s_ctr(0);
  std::vector<apf::MeshEntity*> es_ee_ctr(0);
  m->getPoint(vert, 0, v_pos);

  vd_set_up(m, vert, &es_e_ctr);
  vd_set_up(m, &es_e_ctr, &es_s_ctr);
  vd_set_up(m, &es_s_ctr, &es_ee_ctr);

  std::vector<apf::MeshEntity*> tri(0);
  tri.reserve(es_ee_ctr.size());
  apf::Downward d_t;
  apf::Downward d_v;
  int lookup_ts [4] = {2, 3, 1, 0};

  for(int i = 0; i < es_ee_ctr.size(); i++) {
    m->getDownward(es_ee_ctr.at(i), 0, d_v);
    m->getDownward(es_ee_ctr.at(i), 2, d_t);

    int i1 = findIn(d_v, 4, vert);
    assert(i1 > -1);

    tri.push_back(d_t[lookup_ts[i1]]);
  }
  return tri_dist_sphere(m, v_pos, &tri);
}

// So consider the projection of the normal on triangle plane and consider the 
// edges. If the projection falls outside of the triangle there can be at most 
// 2 edges that the point will lie outside of.
// if 1, take the normal distance of the projection on the edge
// if 2, take the distance of the projection to the corner point
// the allowed radius is the hypothenuse of the normal distance to the plane
// and the distance from the above condition.
std::pair<double, double> tri_dist_sphere(apf::Mesh2* m, apf::Vector3 v_pos, 
                         std::vector<apf::MeshEntity*>* tri) {

  double dist_min = -1;
  double dist_max = -1;

  apf::Vector3 a_pos(0,0,0);
  apf::Vector3 area(0,0,0);

  //std::cout << "v_pos = [" 
  //              << v_pos[0] << ", " << v_pos[1] << ", " 
  //              << v_pos[2] << "];" << std::endl;

  for(int i = 0; i < tri->size(); i++) {

    a_pos = getLinearCentroid(m, tri->at(i));

    area = vd_area_out(m, tri->at(i));
    a_pos = a_pos - v_pos;


    //std::cout << "Tri " << tri->at(i) << std::endl;
    //std::cout << "a_pos " << a_pos << " area " << area << std::endl;


    double dist_curr = norm_0(area)*a_pos;

    if(dist_curr < 0) {
      area = area*(-1);
      dist_curr = - dist_curr;
    }
    apf::Vector3 n_curr(0,0,0);
    apf::Vector3 p_int(0,0,0);
    n_curr = norm_0(area);
    p_int = v_pos + n_curr*(dist_curr);


    //std::cout << "p_int = [" 
    //          << p_int[0] << ", " << p_int[1] << ", " 
    //          << p_int[2] << "];" << std::endl;

    apf::Downward d_v;
    apf::Downward d_e;
    std::vector<apf::Vector3> pos(3, apf::Vector3(0,0,0));
    apf::Vector3 e_pos(0,0,0);
    apf::Vector3 e_dir(0,0,0);
    apf::Vector3 int_pos(0,0,0);

    apf::Vector3 temp1(0,0,0);
    apf::Vector3 temp2(0,0,0);

    m->getDownward(tri->at(i), 0, d_v);
    m->getDownward(tri->at(i), 1, d_e);

    for(int i = 0; i < 3; i++) {
      m->getPoint(d_v[i], 0, pos.at(i));

      //std::cout << "v(" << i << ",:) = [" << pos.at(i)[0] << ", " 
      //          << pos.at(i)[1] << ", " << pos.at(i)[2] << "];"
      //          << std::endl;

    }

    int out_id = -1;
    bool out1 = false;
    bool out2 = false;
    // Loop over vertices. 
    for(int i1 = 0; i1 < 3; i1++) {
      int i2 = (i1+1) % 3;
      int i3 = (i1+2) % 3;
      e_pos = (pos.at(i2) + pos.at(i3)) /2;
      e_dir = norm_0(pos.at(i2) - pos.at(i3));
      temp1 = (p_int - e_pos);
      temp1 = temp1 - e_dir*(temp1*e_dir);
      temp2 = (pos.at(i1) - e_pos);
      temp2 = temp2 - e_dir*(temp2*e_dir);
      // Does the point and other  vertex lie on different sides of the edge:
      if(temp1*temp2 < - std::numeric_limits<double>::min()) {
        // Intersects a vertex:
        if(out1) {
          assert(!out2);
          apf::MeshEntity* v_joint = vd_find_vert(m, d_e[(i1 + 1) % 3], 
                                                     d_e[out_id]);

          assert(v_joint != NULL);
          out2 = true;
          m->getPoint(v_joint, 0, int_pos);
        }
        // Intersects the edge:
        else {
          int_pos = e_pos + e_dir*((p_int - e_pos)*e_dir);
          out_id = (i1 + 1) % 3;
          out1 = true;
        }
      }
    }

    if(out1 or out2)
      dist_curr = (int_pos - v_pos).getLength();

    if(dist_min < 0 or dist_curr < dist_min)
      dist_min = dist_curr;
    if(dist_max < 0 or dist_curr > dist_max)
      dist_max = dist_curr;


    //std::cout << "r_curr = " << dist_curr << std::endl;
  }

  return std::make_pair(dist_min, dist_max);

}

std::pair<double, double> vd_find_dist_v(apf::Mesh2* m, 
                                     std::vector<apf::MeshEntity*>* edge, 
                                     apf::MeshEntity* v) {

  assert(edge->size() > 0);

  double dist_min = -1;
  double dist_max = -1;

  apf::Vector3 temp1(0, 0, 0);
  apf::Vector3 temp2(0, 0, 0);
  apf::Downward d_v;

  double l_curr;
  for(int i = 0; i < edge->size(); i++) {
    m->getDownward(edge->at(i), 0, d_v);
    assert(findIn(d_v, 2, v) > -1);

    m->getPoint(d_v[0], 0, temp1);
    m->getPoint(d_v[1], 0, temp2);

    l_curr = (temp1-temp2).getLength();
    //std::cout << d_v[0] << " " << temp1 << " "
    //          << d_v[1] << " " << temp2 << std::endl;
    //std::cout << "l_curr " << l_curr << std::endl;

    if(dist_min < 0 or dist_min > l_curr) {
      dist_min = l_curr;
    }

    if(dist_max < 0 or dist_max < l_curr) {
      dist_max = l_curr;
    }
  }
  return std::make_pair(dist_min, dist_max);

}

// Used in exterior insertions.
std::pair<double, double> vd_find_dist_v(apf::Mesh2* m, 
                                     std::vector<apf::MeshEntity*>* edge, 
                             apf::MeshEntity* v, apf::Vector3 v_pos) {

  assert(edge->size() > 0);

  double dist_min = -1;
  double dist_max = -1;

  apf::Vector3 temp1(0, 0, 0);
  apf::Downward d_v;

  double l_curr;
  for(int i = 0; i < edge->size(); i++) {
    m->getDownward(edge->at(i), 0, d_v);
    int i1 = findIn(d_v, 2, v);
    assert(i1 > -1);
    int i2 = (i1 + 1) % 2;
    m->getPoint(d_v[i2], 0, temp1);

    l_curr = (temp1-v_pos).getLength();
    //std::cout << d_v[0] << " " << temp1 << " "
    //          << d_v[1] << " " << temp2 << std::endl;
    //std::cout << "l_curr " << l_curr << std::endl;

    if(dist_min < 0 or dist_min > l_curr) {
      dist_min = l_curr;
    }

    if(dist_max < 0 or dist_max < l_curr) {
      dist_max = l_curr;
    }
  }
  return std::make_pair(dist_min, dist_max);

}


// Given the vertex and triangles, find the displacement to the edges 
// across.
std::pair<double, double> vd_find_dist_e(apf::Mesh2* m, 
                                  std::vector<apf::MeshEntity*>* tri, 
                                      apf::MeshEntity* v) {

  assert(tri->size() > 0);

  double dist_min = -1;
  double dist_max = -1;

  apf::Vector3 e_pos(0, 0, 0);

  apf::Vector3 v_pos(0, 0, 0);
  m->getPoint(v, 0, v_pos);
  for(int i = 0; i < tri->size(); i++) {

    apf::Downward d_ve;

    apf::Downward d_v;
    apf::Downward d_e;

    m->getDownward(tri->at(i), 0, d_v);
    m->getDownward(tri->at(i), 1, d_e);

    int v_id = findIn(d_v, 3, v);
    assert(v_id > -1);

    e_pos = getLinearCentroid(m, d_e[lookup_tri_edgex[v_id]]);
    e_pos = e_pos - v_pos;

    m->getDownward(d_e[lookup_tri_edgex[v_id]], 0, d_ve);
    apf::Vector3 v_1(0, 0, 0);
    apf::Vector3 v_2(0, 0, 0);
    m->getPoint(d_ve[0], 0, v_1);
    m->getPoint(d_ve[1], 0, v_2);

    v_1 = v_2 - v_1;
    v_1 = norm_0(v_1);

    e_pos = e_pos - v_1*(e_pos*v_1);
    double dist_curr = e_pos.getLength();

    if(dist_min < 0 or dist_min > dist_curr) {
      dist_min = dist_curr;
    }

    if(dist_max < 0 or dist_max < dist_curr) {
      dist_max = dist_curr;
    }
  }
  
  return std::make_pair(dist_min, dist_max);
}


// Given the vertex and triangles, find the displacement to the edges 
// across, from v_pos. Used in exterior insertions.
std::pair<double, double> vd_find_dist_e(apf::Mesh2* m, 
                                  std::vector<apf::MeshEntity*>* tri, 
                          apf::MeshEntity* v, apf::Vector3 v_pos) {

  assert(tri->size() > 0);

  double dist_min = -1;
  double dist_max = -1;

  apf::Vector3 e_pos(0, 0, 0);

  for(int i = 0; i < tri->size(); i++) {

    apf::Downward d_ve;

    apf::Downward d_v;
    apf::Downward d_e;

    m->getDownward(tri->at(i), 0, d_v);
    m->getDownward(tri->at(i), 1, d_e);

    int v_id = findIn(d_v, 3, v);
    assert(v_id > -1);

    e_pos = getLinearCentroid(m, d_e[lookup_tri_edgex[v_id]]);
    e_pos = e_pos - v_pos;

    m->getDownward(d_e[lookup_tri_edgex[v_id]], 0, d_ve);
    apf::Vector3 v_1(0, 0, 0);
    apf::Vector3 v_2(0, 0, 0);
    m->getPoint(d_ve[0], 0, v_1);
    m->getPoint(d_ve[1], 0, v_2);

    v_1 = v_2 - v_1;
    v_1 = norm_0(v_1);

    e_pos = e_pos - v_1*(e_pos*v_1);
    double dist_curr = e_pos.getLength();

    if(dist_min < 0 or dist_min > dist_curr) {
      dist_min = dist_curr;
    }

    if(dist_max < 0 or dist_max < dist_curr) {
      dist_max = dist_curr;
    }
  }
  
  return std::make_pair(dist_min, dist_max);
}

// Given the vertex position and edges, find the displacement to the vertices 
// across.
std::pair<double, double> vd_find_dist_e(apf::Mesh2* m, 
                                  std::vector<apf::MeshEntity*>* edge, 
                                      apf::Vector3 v_pos) {

  assert(edge->size() > 0);

  double dist_min = -1;
  double dist_max = -1;

  apf::Vector3 e_pos(0, 0, 0);
  apf::Downward d_ve;

  for(int i = 0; i < edge->size(); i++) {

    e_pos = getLinearCentroid(m, edge->at(i));
    e_pos = e_pos - v_pos;

    m->getDownward(edge->at(i), 0, d_ve);
    apf::Vector3 v_1(0, 0, 0);
    apf::Vector3 v_2(0, 0, 0);
    m->getPoint(d_ve[0], 0, v_1);
    m->getPoint(d_ve[1], 0, v_2);

    v_1 = v_2 - v_1;
    v_1 = norm_0(v_1);

    e_pos = e_pos - v_1*(e_pos*v_1);
    double dist_curr = e_pos.getLength();

    if(dist_min < 0 or dist_min > dist_curr) {
      dist_min = dist_curr;
    }

    if(dist_max < 0 or dist_max < dist_curr) {
      dist_max = dist_curr;
    }
  }
  
  return std::make_pair(dist_min, dist_max);
}

// Given the vertex and tetrahedra, find the total of the area in the 
// direction of triangle centroid position, divided by the distance of the 
// triangle. This approximates all triangles to be on a unit sphere and finds
// a position in the center of the volume enclosed by the area. This doesn't 
// work if the volume of the tets is pancake like or has holes in it.
apf::Vector3 vd_find_dist_w(apf::Mesh2* m, std::vector<apf::MeshEntity*>* tet, 
                                                        apf::MeshEntity* v) {

  std::vector<apf::Vector3> a_pos(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> area(0, apf::Vector3(0,0,0));

  a_pos.resize(tet->size());
  area.resize(tet->size());

  apf::Vector3 v_tot(0, 0, 0);
  double vect[3] = {0,0,0};
  v_tot.fromArray(vect);

  double a_tot = 0;

  apf::Vector3 v_pos(0, 0, 0);
  m->getPoint(v, 0, v_pos);
  for(int i = 0; i < tet->size(); i++) {

    apf::Downward d_v;
    apf::Downward d_s;
    m->getDownward(tet->at(i), 0, d_v);
    m->getDownward(tet->at(i), 2, d_s);

    // Calculate the maximum displacement possible before inverting.
    int v_id = findIn(d_v, 4, v);

    a_pos.at(i) = getLinearCentroid(m, d_s[lookup_tet_surfx[v_id]]);

    //apf::MeshElement* ee = createMeshElement(m, d_s[lookup_tet_surfx[v_id]]);
    //double a_curr = measure(ee);
    //destroyMeshElement(ee);

    area.at(i) = vd_area_out(m, d_s[lookup_tet_surfx[v_id]]);
    //std::cout << "Area " << area.at(i).getLength() << std::endl;
    a_pos.at(i) = a_pos.at(i) - v_pos;

    double dist_curr = norm_0(area.at(i))*a_pos.at(i);

    if(dist_curr < 0) {
      area.at(i) = area.at(i)*(-1);
      dist_curr = norm_0(area.at(i))*a_pos.at(i);
    }

    double norm_apos = a_pos.at(i).getLength();
    //v_tot = v_tot + a_pos.at(i)*area.at(i).getLength();
    if (norm_apos < std::numeric_limits<double>::epsilon() )
      norm_apos = 1;

    v_tot = v_tot + norm_0(a_pos.at(i))
                *(norm_0(a_pos.at(i))*area.at(i)/norm_apos);
    a_tot = a_tot + area.at(i).getLength();
  }
  //return v_tot/a_tot;

  v_tot = norm_0(v_tot);

  double dist_min = -1;
  for(int i = 0; i < tet->size(); i++) {

    m->getPoint(v, 0, v_pos);
    double cos_c = norm_0(area.at(i))*v_tot;
    double a_curr = norm_0(area.at(i))*a_pos.at(i);

    if(cos_c > 0) {
      double l_curr = a_curr/cos_c;
      if(l_curr < dist_min or dist_min < 0)
        dist_min = l_curr;
    }
  }

  return v_tot*dist_min;
}

// Given the vertex and triangles, find the total of the length in the 
// direction of triangle centroid position, divided by the distance of the 
// edge across. This approximates all edges to be on a unit sphere and finds
// a position in the center of the area enclosed by the perimeter. This doesn't 
// work if the area of the tris is disc like or has holes in it.
apf::Vector3 vd_find_dist_w_tri(apf::Mesh2* m, 
                    std::vector<apf::MeshEntity*>* tri, apf::MeshEntity* v) {

  std::vector<apf::Vector3> e_pos(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> len(0, apf::Vector3(0,0,0));
  apf::Vector3 temp(0, 0, 0);

  e_pos.resize(tri->size());
  len.resize(tri->size());

  apf::Vector3 a_tot(0, 0, 0);
  double vect[3] = {0,0,0};
  a_tot.fromArray(vect);

  double e_tot = 0;

  apf::Vector3 v_pos(0, 0, 0);
  m->getPoint(v, 0, v_pos);
  for(int i = 0; i < tri->size(); i++) {

    apf::Downward d_v;
    apf::Downward d_e;
    m->getDownward(tri->at(i), 0, d_v);
    m->getDownward(tri->at(i), 1, d_e);

    // Calculate the maximum displacement possible before inverting.
    int v_id = findIn(d_v, 3, v);

    e_pos.at(i) = getLinearCentroid(m, d_e[lookup_tri_edgex[v_id]]);

    //apf::MeshElement* ee = createMeshElement(m, d_e[lookup_tri_edgex[v_id]]);
    //double a_curr = measure(ee);
    //destroyMeshElement(ee);

    m->getPoint(d_v[lookup_tri_x_x[v_id][0]], 0, len.at(i));
    m->getPoint(d_v[lookup_tri_x_x[v_id][1]], 0, temp);
    len.at(i) = len.at(i) - temp;

    len.at(i) = cross(vd_area_out(m, d_e[lookup_tri_edgex[v_id]]), len.at(i));
    //std::cout << "Len " << len.at(i).getLength() << std::endl;
    e_pos.at(i) = e_pos.at(i) - v_pos;

    double dist_curr = norm_0(len.at(i))*e_pos.at(i);

    if(dist_curr < 0) {
      len.at(i) = len.at(i)*(-1);
      dist_curr = norm_0(len.at(i))*e_pos.at(i);
    }

    double norm_epos = e_pos.at(i).getLength();
    //a_tot = a_tot + e_pos.at(i)*len.at(i).getLength();
    if (norm_epos < std::numeric_limits<double>::epsilon() )
      norm_epos = 1;

    a_tot = a_tot + norm_0(e_pos.at(i))
                *(norm_0(e_pos.at(i))*len.at(i)/norm_epos);
    e_tot = e_tot + len.at(i).getLength();
  }
  //return a_tot/e_tot;

  a_tot = norm_0(a_tot);

  double dist_min = -1;
  for(int i = 0; i < tri->size(); i++) {

    m->getPoint(v, 0, v_pos);
    double cos_c = norm_0(len.at(i))*a_tot;
    double a_curr = norm_0(len.at(i))*e_pos.at(i);

    if(cos_c > 0) {
      double l_curr = a_curr/cos_c;
      if(l_curr < dist_min or dist_min < 0)
        dist_min = l_curr;
    }
  }

  return a_tot*dist_min;
}

// Find the distance to the centroid.
double vd_find_dist_c(apf::Mesh2* m, apf::MeshEntity* tet, apf::MeshEntity* v) {

  apf::Vector3 v_pos(0, 0, 0);

  m->getPoint(v, 0, v_pos);
  apf::Vector3 c_pos(0, 0, 0);
  c_pos = getLinearCentroid(m, tet);

  v_pos = v_pos-c_pos;
  return sqrt(v_pos*v_pos);
}

// Given a vertex, set of bound tetrahedra, and a displacement or velocity, 
// remove the component of the displacement or velocity that is not volume 
// preserving for the given collection of tetrahedra.
apf::Vector3 get_volm_pres(apf::Mesh2* m, apf::MeshEntity* v, 
                      std::vector<apf::MeshEntity*>* tet, apf::Vector3 v_in) {

  apf::Vector3 area(0,0,0);
  apf::Vector3 a_pos(0,0,0);
  apf::Vector3 area_eff(0,0,0);

  apf::Vector3 v_pos(0,0,0);
  m->getPoint(v, 0, v_pos);

  double z[3] = {0,0,0};

  for(int i = 0; i < tet->size(); i++) {

    apf::Downward d_v;
    apf::Downward d_s;
    m->getDownward(tet->at(i), 0, d_v);
    m->getDownward(tet->at(i), 2, d_s);

    // Calculate the maximum displacement possible before inverting.
    int v_id = findIn(d_v, 4, v);

    a_pos.fromArray(z);
    a_pos = getLinearCentroid(m, d_s[lookup_tet_surfx[v_id]]);

    area.fromArray(z);
    area = vd_area_out(m, d_s[lookup_tet_surfx[v_id]]);
    a_pos = a_pos - v_pos;

    double dist_curr = norm_0(area)*a_pos;

    if(dist_curr < 0) {
      area = area*(-1);
    }
    //std::cout << "relative pos " << a_pos->at(i) << " area " << area->at(i)
    //          << std::endl;
    area_eff = area_eff + area;
  }
  //std::cout << " a_eff " <<  area_eff << std::endl;

  area_eff = norm_0(area_eff);

  //std::cout << "v_in " << v_in 
  //          << " perp " << area_eff*(area_eff*v_in)
  //          << std::endl;

  v_in = v_in - area_eff*(area_eff*v_in);
  return v_in;
}

apf::Vector3 get_volm_pres2(apf::Mesh2* m, apf::MeshEntity* v, 
                      std::vector<apf::MeshEntity*>* tet, apf::Vector3 v_in) {
  std::vector<apf::Vector3> area(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> a_pos(0, apf::Vector3(0,0,0));

  vd_get_apos(m, v, tet, &area, &a_pos);

  apf::Vector3 area_eff(0,0,0);

  for(int i = 0; i < area.size(); i++) {
    area_eff = area_eff + area.at(i);
  }
  std::cout << " a_eff " <<  area_eff << std::endl;

  area_eff = norm_0(area_eff);

  std::cout << "v_in " << v_in 
            << " perp " << area_eff*(area_eff*v_in)
            << std::endl;

  v_in = v_in - area_eff*(area_eff*v_in);
  return v_in;
}

apf::Vector3 rem_norm_comp(apf::Mesh2* m, std::vector<apf::MeshEntity*>* tri, 
                                  apf::Vector3 v_in) {
  for(int i = 0; i < tri->size(); i++) {
    apf::Vector3 norm(0, 0, 0);
    norm = norm_0(vd_area_out(m, tri->at(i)));
    v_in = v_in - norm*(norm*v_in);
  }
  return v_in;
}

// Given an entity, print out the vertices and their positions. 
void vd_print_vert(apf::Mesh2* m, apf::MeshEntity* ent) {

  int ent_type = m->getType(ent);
  if( m->typeDimension[ent_type] > 0) {
    apf::Downward d_v;

    apf::Vector3 point(0,0,0);
    int dc = m->getDownward(ent, 0, d_v);
    std::cout << ent << ", vertices: " << std::endl;
    for(int i = 0; i < dc; i++) {
      m->getPoint(d_v[i], 0, point);
      std::cout << d_v[i] << " v(" << i+1 << ",:) = [" 
                << point[0] << ", " << point[1] << "," << point[2] << "];"
                << std::endl;
    }
    std::cout << std::endl;
  }
}

std::vector<double> upd_cell_rad(apf::Mesh2* m) {

  double average_ent = ma::getAverageEdgeLength(m);

  std::cout << "Average, edge: " << average_ent << std::endl;
  double cell_meas_max;
  double cell_meas_min;

  std::vector<double> avg_cell(0);
  avg_cell.resize(3);

  std::vector<apf::MeshEntity*> es(0);
  std::vector<apf::MeshEntity*> ev(0);
  std::vector<apf::MeshEntity*> ev_e(0);
  std::vector<apf::MeshEntity*> ev_d(0);

  gmi_model* mdl = m->getModel();

  struct gmi_iter* it;
  struct gmi_ent* e;

  for (int d = 3; d > 0; d--) {
    int counter = 0;

    it = gmi_begin(mdl, d);
    while ((e = gmi_next(mdl, it))) {
      int cell_id = gmi_tag(mdl, e);

      vd_find_ent_geom(m, &es, cell_id, d);
      if(es.size() > 0) {
        std::cout << d << "c" << cell_id << std::endl;
        //vd_find_ent_geom(m, &ev_d, cell_id, dim, 0);
        vd_vert_set(m, &es, &ev);
        ev_e = ev;

        // Store the average cell width for all dimensions for meshadapt.
        double cell_meas_max;
        double cell_meas_min;

        // TODO This bit is going to be replaced with weighted PCA.
        vd_find_ent_geom(m, &ev_d, cell_id, d, 0);

        gmi_set* em_adj = gmi_adjacent(mdl, gmi_find(mdl, d, cell_id), d-1);
        if(em_adj->n > 1)
          vd_remove_set(&ev_e, &ev_d);
        gmi_free_set(em_adj);

        std::tie(cell_meas_max, cell_meas_min) = vd_meas_rad_set_minmax(m, 
                                                                      &ev_e);
        counter++;
        std::cout << cell_meas_min << " counter " << counter << std::endl;

        avg_cell.at(d-1) = avg_cell.at(d-1) + cell_meas_min;
      }
    }
    gmi_end(mdl, it);

    avg_cell.at(d-1) = avg_cell.at(d-1)/counter;
    std::cout << "Average " << d << "c " << avg_cell.at(d-1) << std::endl;
  }
  return avg_cell;
}

std::vector<double> upd_cell_rad(apf::Mesh2* m, vd_entlist* e_list) {
  double average_ent = ma::getAverageEdgeLength(m);

  std::cout << "Average, edge: " << average_ent << std::endl;
  double cell_meas_max;
  double cell_meas_min;

  std::vector<double> avg_cell(3);

  std::vector<apf::MeshEntity*> ev(0);
  std::vector<apf::MeshEntity*> ev_e(0);
  std::vector<apf::MeshEntity*> ev_d(0);

  gmi_model* mdl = m->getModel();

  struct gmi_iter* it;
  struct gmi_ent* e;

  for (int d = 3; d > 0; d--) {

    it = gmi_begin(mdl, d);
    int counter = 0;
    while ((e = gmi_next(mdl, it))) {
      int cell_id = gmi_tag(mdl, e);

      if(e_list->e.at(d).at(cell_id-1).at(d).size() > 0) {
        std::cout << d << "c" << cell_id << std::endl;
        //vd_find_ent_geom(m, &ev_d, cell_id, dim, 0);
        vd_vert_set(m, &e_list->e.at(d).at(cell_id-1).at(d), &ev);
        ev_e = ev;

        // Store the average cell width for all dimensions for meshadapt.
        double cell_meas_max;
        double cell_meas_min;

        gmi_set* em_adj = gmi_adjacent(mdl, gmi_find(mdl, d, cell_id), d-1);
        if(em_adj->n > 1)
          vd_remove_set(&ev_e, &e_list->e.at(d).at(cell_id-1).at(0));
        gmi_free_set(em_adj);

        std::tie(cell_meas_max, cell_meas_min) = vd_meas_rad_set_minmax(m, 
                                                                      &ev_e);
        counter++;
        std::cout << cell_meas_min << " counter " << counter << std::endl;

        avg_cell.at(d-1) = avg_cell.at(d-1) + cell_meas_min;
      }
    }
    gmi_end(mdl, it);

    avg_cell.at(d-1) = avg_cell.at(d-1)/counter;
    std::cout << "Average " << d << "c " << avg_cell.at(d-1) << std::endl;
  }
  return avg_cell;

}

// Given a cloud of vertices, check if the plane intersects the cloud.
int pl_int_verts(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ent, 
                       apf::Vector3 p_pos, apf::Vector3 p_norm, double tol) {
  // TODO
  // Tolerance is not implemented as it would require change of logic ordering.
  tol = std::fabs(tol);
  assert(tol < 1);
  assert(ent->size() > 1);
  int ent_type = m->getType(ent->at(0));
  assert(ent_type == apf::Mesh::VERTEX);

  apf::Vector3 p1(0,0,0);
  apf::Vector3 p2(0,0,0);
  m->getPoint(ent->at(0), 0, p1);

  p1 = p1 - p_pos;
  double dir = p1*p_norm;
  if(std::fabs(dir) < std::numeric_limits<double>::epsilon()) {
    std::cout << "Intersecting at a point" << std::endl;
    return -1;
  }
  else {
    for(int i = 1; i < ent->size(); i++) {
      m->getPoint(ent->at(i), 0, p2);
      p2 = p2 - p_pos;
      double dir2 = p2*p_norm;
      if(std::fabs(dir2) < std::numeric_limits<double>::epsilon()) {
        std::cout << "Intersecting at a point" << std::endl;
        return -i-1;
      }
      else if (dir*dir2 < 0) {
        std::cout << "Intersecting" << std::endl;
        return 1;
      }
    }
    return 0;
  }

}

// Given an edge, check if the plane intersects the edge.
int pl_int_edge(apf::Mesh2* m, apf::MeshEntity* ent, apf::Vector3 p_pos,
                                  apf::Vector3 p_norm, double tol) {
  int ent_type = m->getType(ent);
  assert(ent_type == apf::Mesh::EDGE);
  apf::Downward down;

  apf::Vector3 p1(0,0,0);
  apf::Vector3 p2(0,0,0);
  apf::Vector3 mid(0,0,0);
  m->getDownward(ent, 0, down);
  m->getPoint(down[0], 0, p1);
  m->getPoint(down[1], 0, p2);

  mid = (p1 + p2)/2;
  p1 = (p1 - mid)*(1-tol) + mid;
  p2 = (p2 - mid)*(1-tol) + mid;

  p1 = p1 - p_pos;
  p2 = p2 - p_pos;

  if(std::fabs(p1*p_norm) < std::numeric_limits<double>::epsilon()) {
    std::cout << "Intersecting at a point" << std::endl;
    return -1;
  }
  else if (std::fabs(p2*p_norm) < std::numeric_limits<double>::epsilon()) {
    std::cout << "Intersecting at a point" << std::endl;
    return -2;
  }
  else if(p1*p_norm * (p2*p_norm) < 0) {
    std::cout << "Intersecting" << std::endl;
    return 1;
  }
  else {
    std::cout << "Not intersecting" << std::endl;
    return 0;
  }
}

// Given an edge, check if the plane intersects the edge.
void pl_int_edge(apf::Mesh2* m, apf::MeshEntity* ent, vd_plane* pl, 
                               vd_inter* v_int, double tol) {

  int ent_type = m->getType(ent);
  assert(ent_type == apf::Mesh::EDGE);
  apf::Downward down;

  apf::Vector3 p1(0,0,0);
  apf::Vector3 p2(0,0,0);
  m->getDownward(ent, 0, down);
  m->getPoint(down[0], 0, p1);
  m->getPoint(down[1], 0, p2);

  // Assuming same length (roughly as we are doing additional divisions by 
  // cutting), this could be used(roughly) to calculate how small the angle 
  // between the vector lying in the intersection of the plane and the 
  // triangle and the direction vectors along the edges bounding the triangle
  // and bounded by the central vertex. 
  double w = std::fabs((p1 - p2)*pl->norm);

  if(std::fabs((p1 - pl->pos)*pl->norm) 
                            < tol*w) {
    //                        < std::numeric_limits<double>::epsilon()) {
    std::cout << "Intersecting at a point" << std::endl;
    v_int->cut = -1;
    v_int->pos = p1;
    std::cout << v_int->pos << std::endl;
  }
  else if(std::fabs((p2 - pl->pos)*pl->norm) 
                            < tol*w) {
    //                        < std::numeric_limits<double>::epsilon()) {
    std::cout << "Intersecting at a point" << std::endl;
    v_int->cut = -2;
    v_int->pos = p2;
    std::cout << v_int->pos << std::endl;
  }
  else if(((p1 - pl->pos)*pl->norm) * ((p2 - pl->pos)*pl->norm) < 0) {
    std::cout << "Intersecting" << std::endl;
    v_int->cut = 1;
    apf::Vector3 norm_dir(0, 0, 0);
    norm_dir = norm_0(p1-p2);
    double ang_cos = norm_dir*pl->norm;
    v_int->pos = p2+norm_dir*(pl->norm*(pl->pos-p2))/ang_cos;
    std::cout << v_int->pos << std::endl;
  }
  else {
    std::cout << "Not intersecting" << std::endl;
    v_int->cut = 0;
    double z[3] = {0,0,0};
    v_int->pos.fromArray(z);
  }
}

// Given two points, check if the plane lies between them.
void pl_int_2pts(apf::Mesh2* m, apf::Vector3 p1, apf::Vector3 p2, 
                                              vd_plane* pl, vd_inter* v_int) {
  if(std::fabs((p1 - pl->pos)*pl->norm) 
                            < std::numeric_limits<double>::min()) {
    //std::cout << "Intersecting at a point" << std::endl;
    v_int->cut = -1;
    v_int->pos = p1;
    std::cout << v_int->pos << std::endl;
  }
  else if(std::fabs((p2 - pl->pos)*pl->norm) 
                            < std::numeric_limits<double>::min()) {
    //std::cout << "Intersecting at a point" << std::endl;
    v_int->cut = -2;
    v_int->pos = p2;
    std::cout << v_int->pos << std::endl;
  }
  else if(((p1 - pl->pos)*pl->norm) * ((p2 - pl->pos)*pl->norm) <  
                                      std::numeric_limits<double>::min()) {
    //std::cout << "Intersecting" << std::endl;
    v_int->cut = 1;
    apf::Vector3 norm_dir(0,0,0);
    norm_dir = norm_0(p1-p2);
    double ang_cos = norm_dir*pl->norm;
    v_int->pos = p2+norm_dir*(pl->norm*(pl->pos-p2))/ang_cos;
    std::cout << v_int->pos << std::endl;
  }
  else {
    //std::cout << "Not intersecting" << std::endl;
    v_int->cut = 0;
    double z[3] = {0,0,0};
    v_int->pos.fromArray(z);
  }
}

// Given a triangle, check if the plane intersects the triangle.
int pl_int_tri(apf::Mesh2* m, apf::MeshEntity* ent, apf::Vector3 p_pos, 
                                                      apf::Vector3 p_norm) {
  int ent_type = m->getType(ent);
  assert(ent_type == apf::Mesh::TRIANGLE);
  apf::Downward down;

  apf::Vector3 p1(0,0,0);
  apf::Vector3 p2(0,0,0);
  m->getDownward(ent, 0, down);
  m->getPoint(down[0], 0, p1);

  p1 = p1 - p_pos;
  double dir = p1*p_norm;
  if(std::fabs(dir) < std::numeric_limits<double>::epsilon()) {
    std::cout << "Intersecting at a point" << std::endl;
    return -1;
  }
  else {
    for(int i = 1; i < 3; i++) {
      m->getPoint(down[i], 0, p2);
      p2 = p2 - p_pos;
      double dir2 = p2*p_norm;
      if(std::fabs(dir2) < std::numeric_limits<double>::epsilon()) {
        std::cout << "Intersecting at a point" << std::endl;
        return -i-1;
      }
      else if (dir*dir2 < 0) {
        std::cout << "Intersecting" << std::endl;
        return 1;
      }
    }
    return 0;
  }
}

// Find the intersection of a line and a plane. Assume they intersect.
apf::Vector3 pl_int_line(apf::Vector3 l_pos, apf::Vector3 l, 
                          apf::Vector3 pl_pos, apf::Vector3 pl_norm) {
  apf::Vector3 p_rel(0,0,0);
  apf::Vector3 l_dir(0,0,0);
  apf::Vector3 n(0,0,0);

  p_rel = pl_pos - l_pos;
  n = norm_0(pl_norm);
  l_dir = norm_0(l);

  double dist = p_rel*n;
  double a_cos = n*l_dir;
  return l_pos + l_dir*dist/a_cos;
}

// Does the line passing from o to l intersect the triangles v1v2v3. 
bool tri_int_line(apf::Vector3 v1, apf::Vector3 v2, apf::Vector3 v3, 
                  apf::Vector3 o, apf::Vector3 l, apf::Vector3* int_pt) {

  apf::Vector3 temp1(0,0,0);
  apf::Vector3 temp2(0,0,0);
  temp2 = (v1+v2+v3)/3;

  // The area normal is outwards.
  temp1 = vd_area_out_n(v1, v2, v3);
  temp2 = temp2 - o;

  double dist_curr = norm_0(temp1)*temp2;

  if(dist_curr < 0) {
    temp1 = temp1*(-1);
  }
  apf::Vector3 n_curr(0,0,0);
  apf::Vector3 p_int(0,0,0);

  n_curr = norm_0(temp1);
  l = norm_0(l-o);
  temp1 = n_curr*(n_curr*l);
  temp2 = n_curr*(n_curr*temp2);
  std::cout << "v(1,:) = [" << v1[0] << "," << v1[1] << "," << v1[2] << "];"
            << "v(2,:) = [" << v2[0] << "," << v2[1] << "," << v2[2] << "];"
            << "v(3,:) = [" << v3[0] << "," << v3[1] << "," << v3[2] << "];"
            << "n = [" << n_curr[0] << "," << n_curr[1] << "," << n_curr[2] << "];"
            << "l = [" << l[0] << "," << l[1] << "," << l[2] << "];" 
            <<std::endl;
  // Is the line directed towards the triangle. If not return false.
  if(temp1*temp2 < std::numeric_limits<double>::min()) {
    //std::cout << "temp1*temp2 " << std::endl;
    return false;
  }
  // Is the line direction inplane.
  if(std::fabs(n_curr*l) < std::numeric_limits<double>::min()) {
    //std::cout << "n_curr*l " << std::endl;
    return false;
  }

  p_int = pl_int_line(o, l, v1, n_curr);

  apf::Vector3 e_pos(0,0,0);
  apf::Vector3 e_dir(0,0,0);

  e_pos = (v1 + v2)/2;
  e_dir = norm_0(v2 - v1);
  temp1 = (p_int - e_pos);
  temp1 = temp1 - e_dir*(temp1*e_dir);
  temp2 = (v3 - e_pos);
  temp2 = temp2 - e_dir*(temp2*e_dir);
  if( temp1*temp2 < - std::numeric_limits<double>::min()) {
    //std::cout << "v1v2 - v3" << temp1*temp2 << std::endl;
    return false;
  }

  e_pos = (v2 + v3)/2;
  e_dir = norm_0(v3 - v2);
  temp1 = (p_int - e_pos);
  temp1 = temp1 - e_dir*(temp1*e_dir);
  temp2 = (v1 - e_pos);
  temp2 = temp2 - e_dir*(temp2*e_dir);
  if( temp1*temp2 < - std::numeric_limits<double>::min()) {
    //std::cout << "v2v3 - v1" << temp1*temp2 << std::endl;
    return false;
  }

  e_pos = (v1 + v3)/2;
  e_dir = norm_0(v3 - v1);
  temp1 = (p_int - e_pos);
  temp1 = temp1 - e_dir*(temp1*e_dir);
  temp2 = (v2 - e_pos);
  temp2 = temp2 - e_dir*(temp2*e_dir);
  if( temp1*temp2 < - std::numeric_limits<double>::min()) {
    //std::cout << "v1v3 - v2" << temp1*temp2 << std::endl;
    return false;
  }

  if(int_pt != NULL)
    *int_pt = p_int;
  return true;
}


// Does the line passing from o to l intersect the triangles v1v2v3. 
bool tri_int_line(apf::Mesh2* m, apf::MeshEntity* tri,
                  apf::Vector3 o, apf::Vector3 l, apf::Vector3* int_pt) {
  apf::Vector3 v1(0,0,0);
  apf::Vector3 v2(0,0,0);
  apf::Vector3 v3(0,0,0);

  int ent_type = m->getType(tri);
  int d = m->typeDimension[ent_type];
  assert(d == 2);
  apf::Downward d_v;
  m->getDownward(tri, 0, d_v);
  m->getPoint(d_v[0], 0, v1);
  m->getPoint(d_v[1], 0, v2);
  m->getPoint(d_v[2], 0, v3);

  return tri_int_line(v1, v2, v3, o, l, int_pt);
}


// Does a point lie on a triangle. Assume coplanar.
bool pt_int_tri(apf::Mesh2* m, apf::MeshEntity* tri, apf::Vector3 pt, double tol) {
  int ent_type = m->getType(tri);
  int d = m->typeDimension[ent_type];
  assert(d == 2);
  apf::Vector3 v1(0,0,0);
  apf::Vector3 v2(0,0,0); 
  apf::Vector3 v3(0,0,0);
  apf::Downward d_v;
  m->getDownward(tri, 0, d_v);
  m->getPoint(d_v[0], 0, v1);
  m->getPoint(d_v[1], 0, v2);
  m->getPoint(d_v[2], 0, v3);

  return pt_int_tri(v1,v2,v3,pt, tol);
}
// Does a point lie on a triangle. Assume coplanar.
bool pt_int_tri(apf::Vector3 v1, apf::Vector3 v2, apf::Vector3 v3, 
                  apf::Vector3 pt, double tol) {
  tol = std::fabs(tol);

  apf::Vector3 e_pos(0,0,0);
  apf::Vector3 e_dir(0,0,0);
  apf::Vector3 temp1(0,0,0);
  apf::Vector3 temp2(0,0,0);

  e_pos = (v1 + v2)*(1-tol)/2 + v3*tol;
  e_dir = norm_0(v2 - v1);
  temp1 = (pt - e_pos);
  temp1 = temp1 - e_dir*(temp1*e_dir);
  temp2 = (v3 - e_pos);
  temp2 = temp2 - e_dir*(temp2*e_dir);
  if( temp1*temp2 < - std::numeric_limits<double>::min())
    return false;
  e_pos = (v2 + v3)*(1-tol)/2 + v1*tol;
  e_dir = norm_0(v3 - v2);
  temp1 = (pt - e_pos);
  temp1 = temp1 - e_dir*(temp1*e_dir);
  temp2 = (v1 - e_pos);
  temp2 = temp2 - e_dir*(temp2*e_dir);
  if( temp1*temp2 < - std::numeric_limits<double>::min())
    return false;

  e_pos = (v1 + v3)*(1-tol)/2 + v2*tol;
  e_dir = norm_0(v3 - v1);
  temp1 = (pt - e_pos);
  temp1 = temp1 - e_dir*(temp1*e_dir);
  temp2 = (v2 - e_pos);
  temp2 = temp2 - e_dir*(temp2*e_dir);
  if( temp1*temp2 < - std::numeric_limits<double>::min())
    return false;

  return true;
}

// TODO maShape getSliverCode seems to be able to do triangle projection

// Does the line passing from o to l intersect the triangles v1v2v3. 
// tol is the percentage shift to be applied to the edge center position
// when calculating whether the intersection is within the triangle to make the
// triangle effectively smaller.
bool tri_int_edge(apf::Vector3 v1, apf::Vector3 v2, apf::Vector3 v3, 
           apf::Vector3 o, apf::Vector3 l, double tol, apf::Vector3* int_pt) {
  tol = std::fabs(tol);
  apf::Vector3 l_dir(0,0,0);

  l_dir = l-o;
  double e_len = l_dir.getLength();
  l_dir = norm_0(l_dir);

  apf::Vector3 temp1(0,0,0);
  apf::Vector3 temp2(0,0,0);
  temp2 = (v1+v2+v3)/3;

  // The area normal is outwards.
  temp1 = vd_area_out_n(v1, v2, v3);
  temp2 = temp2 - o;

  apf::Vector3 n_curr(0,0,0);
  n_curr = norm_0(temp1);
  double dist_curr = n_curr*temp2;

  if(dist_curr < -std::numeric_limits<double>::min()) {
    n_curr = n_curr*(-1);
  }
  apf::Vector3 p_int(0,0,0);

  temp1 = n_curr*(n_curr*l_dir);
  temp2 = n_curr*(n_curr*temp2);
  // Is the line directed towards the triangle. If not return false.
  if(temp1*temp2 < std::numeric_limits<double>::min())
    return false;
  // Is the line direction inplane.
  if(std::fabs(n_curr*l_dir) < std::numeric_limits<double>::min())
    return false;
  // The intersection of the plane with the line.
  p_int = pl_int_line(o, l_dir, v1, n_curr);
  double i_len = (p_int - o)*l_dir;
  if(i_len < std::numeric_limits<double>::min() or 
      i_len > e_len)
    return false;

  std::vector<bool> out(3,false);

  apf::Vector3 e_pos(0,0,0);
  apf::Vector3 e_dir(0,0,0);

  e_pos = (v1 + v2)*(1-tol)/2 + v3*tol;
  e_dir = norm_0(v2 - v1);
  temp1 = (p_int - e_pos);
  temp1 = temp1 - e_dir*(temp1*e_dir);
  temp2 = (v3 - e_pos);
  temp2 = temp2 - e_dir*(temp2*e_dir);
  if( temp1*temp2 < - std::numeric_limits<double>::min())
    return false;
  e_pos = (v2 + v3)*(1-tol)/2 + v1*tol;
  e_dir = norm_0(v3 - v2);
  temp1 = (p_int - e_pos);
  temp1 = temp1 - e_dir*(temp1*e_dir);
  temp2 = (v1 - e_pos);
  temp2 = temp2 - e_dir*(temp2*e_dir);
  if( temp1*temp2 < - std::numeric_limits<double>::min())
    return false;

  e_pos = (v1 + v3)*(1-tol)/2 + v2*tol;
  e_dir = norm_0(v3 - v1);
  temp1 = (p_int - e_pos);
  temp1 = temp1 - e_dir*(temp1*e_dir);
  temp2 = (v2 - e_pos);
  temp2 = temp2 - e_dir*(temp2*e_dir);
  if( temp1*temp2 < - std::numeric_limits<double>::min())
    return false;

  if(int_pt != NULL)
    *int_pt = p_int;
  return true;
}

// Find the intersection of the line segment and the sphere, centered at c_pos
// with radius r. Assume the line segment starts at o, inside the sphere, and 
// ends at l outside the sphere.
apf::Vector3 sph_int_line(apf::Vector3 c_pos, double r, 
                          apf::Vector3 l, apf::Vector3 o) {
  apf::Vector3 oc(0,0,0);

  l = l - o;
  l = norm_0(l);

  oc = o-c_pos;
/*
  std::cout << " o " << o << " "
            << " l " << l
            << " oc " << oc
            << std::endl;
*/
  double a = l*oc;
  double b = oc.getLength();
  b = a*a - b*b + r*r;

  //std::cout << "a " << a << " b " << b << std::endl;

  assert(b > 0);

  double d = - a + sqrt(b);
  assert(d > 0);
  return o+l*d;
}

// Calculate the incenter of an entity.
apf::Vector3 inctr_ent(apf::Mesh2* m, apf::MeshEntity* ent) {
  apf::Vector3 i_ctr(0,0,0);
  apf::Downward d_v;

  if(m->getType(ent) == apf::Mesh::VERTEX) {
    m->getPoint(ent, 0, i_ctr);
  }
  else if(m->getType(ent) == apf::Mesh::EDGE) {
    i_ctr = vd_get_pos(m, ent);
  }
  else if(m->getType(ent) == apf::Mesh::TRIANGLE) {
    std::vector<apf::Vector3> v(3, apf::Vector3(0,0,0));
    std::vector<apf::Vector3> dir(3, apf::Vector3(0,0,0));
    std::vector<apf::Vector3> n(3, apf::Vector3(0,0,0));
    std::vector<double> Ax(3);
    std::vector<double> Ay(3);
    std::vector<double> a(3);

    apf::Vector3 o1(0,0,0);
    apf::Vector3 o2(0,0,0);

    m->getDownward(ent, 0, d_v);
    for(int i = 0; i < 3; i++) {
      m->getPoint(d_v[i], 0, v.at(i));
    }

    dir.at(0) = v.at(1) - v.at(2);
    if(dir.at(0).getLength() < std::numeric_limits<double>::min()) {
      i_ctr = v.at(1)/2 + v.at(2)/2;
      return i_ctr;
    }
    dir.at(1) = v.at(0) - v.at(2);
    if(dir.at(1).getLength() < std::numeric_limits<double>::min()) {
      i_ctr = v.at(0)/2 + v.at(2)/2;
      return i_ctr;
    }
    dir.at(2) = v.at(0) - v.at(1);
    if(dir.at(2).getLength() < std::numeric_limits<double>::min()) {
      i_ctr = v.at(0)/2 + v.at(1)/2;
      return i_ctr;
    }

    double p = 0;
    for(int i = 0; i < 3; i++) {
      a.at(i) = dir.at(i).getLength();
      p = p + a.at(i);
    }

    n.at(0) = norm_0(dir.at(0));
    n.at(1) = norm_0(dir.at(1));
    n.at(2) = norm_0(dir.at(2));

    o1 = dir.at(1) - n.at(2)*(dir.at(1)*n.at(2));
    o1 = norm_0(o1);
    o2 = norm_0(dir.at(2));

    for(int i = 0; i < 3; i++) {
      Ax.at(i) = v.at(i)*o1;
      Ay.at(i) = v.at(i)*o2;
    }

    double i_localx = 0;
    double i_localy = 0;
    for(int i = 0; i < 3; i++) {
      i_localx = i_localx + a.at(i)/p*Ax.at(i);
      i_localy = i_localy + a.at(i)/p*Ay.at(i);
    }

    i_ctr = v.at(2) + o1*i_localx + o2*i_localy;

  }
  // Find the bisecting planes at two disconnected edges and a third edge.
  // Find their intersection.
  else {
    assert(m->getType(ent) == apf::Mesh::TET);

    // The indices of the vertices and edges to find bisecting planes. 
    // 0: origin vertex, 1: Id of the joint edge, 2-3: id of the
    // other edges bounded by the origin vertex 
    int bisect_trip [3][4] = {{0,3,0,2}, {0,0,2,3}, {1,1,0,4}};

    std::vector<apf::Vector3> v(4, apf::Vector3(0,0,0));
    std::vector<apf::Vector3> dir(6, apf::Vector3(0,0,0));
    std::vector<apf::Vector3> l(3, apf::Vector3(0,0,0));
    std::vector<apf::Vector3> b(3, apf::Vector3(0,0,0));
    std::vector<apf::Vector3> n(3, apf::Vector3(0,0,0));
    //apf::Vector3 l(0,0,0);
    apf::Vector3 e1(0,0,0);
    apf::Vector3 e2(0,0,0);

    m->getDownward(ent, 0, d_v);

    apf::Downward d_t;
    m->getDownward(ent, 0, d_t);
    for(int i = 0; i < 4; i++) {
      m->getPoint(d_v[i], 0, v.at(i));
    }
    for(int i = 0; i < 6; i++) {
      int i0 = vd_tet_edge_vert[i][0];
      int i1 = vd_tet_edge_vert[i][1];
      dir.at(i) = v.at(i1) - v.at(i0);
      dir.at(i) = norm_0(dir.at(i));
      // If tet is singular, return the incenter of one of the triangles
      // not bounded by the current edge.
      if(dir.at(i).getLength() < std::numeric_limits<double>::min())
        return inctr_ent(m, d_t[vd_tet_edge_c_tri[i][0]]);
    }

    for(int i = 0; i < 3; i++) {
      int i0 = bisect_trip[i][0];
      int e_id1 = bisect_trip[i][1];
      int e_id2 = bisect_trip[i][2];
      int e_id3 = bisect_trip[i][3];
      l.at(i) = dir.at(e_id1);

      // The edge direction vectors should point out of the corner vertex.
      if(vd_tet_edge_vert[e_id2][0] < i0)
        e1 = dir.at(e_id2)*(-1);
      else
        e1 = dir.at(e_id2);
      if(vd_tet_edge_vert[e_id3][0] < i0)
        e2 = dir.at(e_id3)*(-1);
      else
        e2 = dir.at(e_id3);

      e1 = e1 - l.at(i)*(e1*l.at(i));
      e2 = e2 - l.at(i)*(e1*l.at(i));
      e1 = norm_0(e1);
      e2 = norm_0(e2);

      b.at(i) = norm_0(e1+e2);
      n.at(i) = norm_0(cross(l.at(i), b.at(i)));
      assert(n.at(i).getLength() > std::numeric_limits<double>::min());
      std::cout << "Plane " << i
                << " pt " << v.at(i0)
                << " norm " << n.at(i)
                << std::endl;
    }
    e1 = norm_0(cross(n.at(0), n.at(1)));
    assert(e1.getLength() > std::numeric_limits<double>::min());

    // Intersect line starting at v0 along e1 with the plane 3.
    i_ctr = pl_int_line(v.at(bisect_trip[0][0]), e1, 
                        v.at(bisect_trip[2][0]), n.at(2));
  }
  return i_ctr;
}


// Calculate the incenter of a set of vertices. If v.size() > 4, return the 
// center of mass.
apf::Vector3 inctr_ent(std::vector<apf::Vector3>* v) {
  apf::Vector3 i_ctr(0,0,0);
  assert(v->size() > 0);
  if(v->size() == 1) {
    return v->at(0);
  }
  else if(v->size() == 2) {
    return v->at(0)/2 + v->at(1)/2;
  }
  else if(v->size() == 3) {
    std::vector<apf::Vector3> dir(3, apf::Vector3(0,0,0));
    std::vector<apf::Vector3> n(3, apf::Vector3(0,0,0));
    std::vector<double> Ax(3);
    std::vector<double> Ay(3);
    std::vector<double> a(3);

    apf::Vector3 o1(0,0,0);
    apf::Vector3 o2(0,0,0);

    dir.at(0) = v->at(1) - v->at(2);
    if(dir.at(0).getLength() < std::numeric_limits<double>::min()) {
      i_ctr = v->at(1)/2 + v->at(2)/2;
      return i_ctr;
    }
    dir.at(1) = v->at(0) - v->at(2);
    if(dir.at(1).getLength() < std::numeric_limits<double>::min()) {
      i_ctr = v->at(0)/2 + v->at(2)/2;
      return i_ctr;
    }
    dir.at(2) = v->at(0) - v->at(1);
    if(dir.at(2).getLength() < std::numeric_limits<double>::min()) {
      i_ctr = v->at(0)/2 + v->at(1)/2;
      return i_ctr;
    }

    double p = 0;
    for(int i = 0; i < 3; i++) {
      a.at(i) = dir.at(i).getLength();
      p = p + a.at(i);
    }

    n.at(0) = norm_0(dir.at(0));
    n.at(1) = norm_0(dir.at(1));
    n.at(2) = norm_0(dir.at(2));

    o1 = dir.at(1) - n.at(2)*(dir.at(1)*n.at(2));
    o1 = norm_0(o1);
    o2 = norm_0(dir.at(2));

    for(int i = 0; i < 3; i++) {
      Ax.at(i) = v->at(i)*o1;
      Ay.at(i) = v->at(i)*o2;
    }

    double i_localx = 0;
    double i_localy = 0;
    for(int i = 0; i < 3; i++) {
      i_localx = i_localx + a.at(i)/p*Ax.at(i);
      i_localy = i_localy + a.at(i)/p*Ay.at(i);
    }

    i_ctr = v->at(2) + o1*i_localx + o2*i_localy;

  }
  // Find the bisecting planes at two disconnected edges and a third edge.
  // Find their intersection.
  else if(v->size() == 4) {

    // The indices of the vertices and edges to find bisecting planes. 
    // 0: origin vertex, 1: Id of the joint edge, 2-3: id of the
    // other edges bounded by the origin vertex 
    int bisect_trip [3][4] = {{0,3,0,2}, {0,0,2,3}, {1,1,0,4}};

    std::vector<apf::Vector3> dir(6, apf::Vector3(0,0,0));
    std::vector<apf::Vector3> l(3, apf::Vector3(0,0,0));
    std::vector<apf::Vector3> b(3, apf::Vector3(0,0,0));
    std::vector<apf::Vector3> n(3, apf::Vector3(0,0,0));
    //apf::Vector3 l(0,0,0);
    apf::Vector3 e1(0,0,0);
    apf::Vector3 e2(0,0,0);

    for(int i = 0; i < 6; i++) {
      int i0 = vd_tet_edge_vert[i][0];
      int i1 = vd_tet_edge_vert[i][1];
      dir.at(i) = v->at(i1) - v->at(i0);
      dir.at(i) = norm_0(dir.at(i));
      // If tet is singular, return the incenter of one of the triangles
      // not bounded by the current edge.
      if(dir.at(i).getLength() < std::numeric_limits<double>::min()) {
        l.at(0) = v->at(i0);
        l.at(1) = v->at(vd_tet_vv[i0][i1][0]);
        l.at(2) = v->at(vd_tet_vv[i0][i1][1]);
        return inctr_ent(&l);
      }
    }

    for(int i = 0; i < 3; i++) {
      int i0 = bisect_trip[i][0];
      int e_id1 = bisect_trip[i][1];
      int e_id2 = bisect_trip[i][2];
      int e_id3 = bisect_trip[i][3];
      l.at(i) = dir.at(e_id1);

      // The edge direction vectors should point out of the corner vertex.
      if(vd_tet_edge_vert[e_id2][0] < i0)
        e1 = dir.at(e_id2)*(-1);
      else
        e1 = dir.at(e_id2);
      if(vd_tet_edge_vert[e_id3][0] < i0)
        e2 = dir.at(e_id3)*(-1);
      else
        e2 = dir.at(e_id3);

      e1 = e1 - l.at(i)*(e1*l.at(i));
      e2 = e2 - l.at(i)*(e1*l.at(i));
      e1 = norm_0(e1);
      e2 = norm_0(e2);

      b.at(i) = norm_0(e1+e2);
      n.at(i) = norm_0(cross(l.at(i), b.at(i)));
      assert(n.at(i).getLength() > std::numeric_limits<double>::min());
    }
    e1 = norm_0(cross(n.at(0), n.at(1)));
    assert(e1.getLength() > std::numeric_limits<double>::min());

    // Intersect line starting at v0 along e1 with the plane 3.
    i_ctr = pl_int_line(v->at(bisect_trip[0][0]), e1, 
                        v->at(bisect_trip[2][0]), n.at(3));
  }
  return i_ctr;
}


// Given pos, return pos mirrored by plane defined by ori and dir.
apf::Vector3 mir_pos_pl(apf::Vector3 ori, apf::Vector3 dir_m, apf::Vector3 pos) {
  return ori - dir_m*(dir_m*(ori - pos))*2;
}

// Given pos, return pos mirrored by line defined by ori and dir.
apf::Vector3 mir_pos_ln(apf::Vector3 ori, apf::Vector3 dir_m, apf::Vector3 pos) {
  apf::Vector3 rel_pos(0,0,0);
  rel_pos = ori - pos;
  return ori - (rel_pos - dir_m*(dir_m*rel_pos))*2;
}

// Given dir, return dir mirrored by plane defined dir_m.
apf::Vector3 mir_dir_pl(apf::Vector3 dir_m, apf::Vector3 dir) {
  return dir - (dir_m*(dir_m*dir))*2;
}

// Given dir, return dir mirrored by line defined by dir_m.
apf::Vector3 mir_dir_ln(apf::Vector3 dir_m, apf::Vector3 dir) {
  apf::Vector3 n_dir(0,0,0);
  n_dir = dir_m*(dir_m*dir);
  return n_dir - (dir - n_dir)*2;
}

// Given pos, return pos rotated around ori by rot.
apf::Vector3 rot_pos(apf::Vector3 ori, apf::Matrix3x3* rot, apf::Vector3 pos) {
  return ori + (*rot) * (pos-ori);
}

// Given dir, return dir rotated by rot.
apf::Vector3 rot_dir(apf::Matrix3x3* rot, apf::Vector3 dir) {
  return (*rot) * dir;
}

// Given dir and angle, determine the rotation matrix.
apf::Matrix3x3 rot_matrix_dir(apf::Vector3 dir, double angle) {
  double cos_ang = std::cos(angle);
  double sin_ang = std::sin(angle);

  double r11 = cos_ang + dir[0]*dir[0]*(1-cos_ang);
  double r12 = dir[0]*dir[1]*(1-cos_ang) - dir[2]*sin_ang;
  double r13 = dir[0]*dir[2]*(1-cos_ang) + dir[1]*sin_ang;

  double r21 = dir[0]*dir[1]*(1-cos_ang) + dir[2]*sin_ang;
  double r22 = cos_ang + dir[1]*dir[1]*(1-cos_ang);
  double r23 = dir[1]*dir[2]*(1-cos_ang) - dir[0]*sin_ang;

  double r31 = dir[0]*dir[2]*(1-cos_ang) - dir[1]*sin_ang;
  double r32 = dir[1]*dir[2]*(1-cos_ang) + dir[0]*sin_ang;
  double r33 = cos_ang + dir[2]*dir[2]*(1-cos_ang);

  return apf::Matrix3x3(r11, r12, r13, r21, r22, r23, r31, r32, r33);
}

// Find the center and radius of circle passing through three points.
// https://stackoverflow.com/questions/13977354/build-circle-from-3-points-in-3d-space-implementation-in-c-or-c
// TODO check if correct
/*
std::pair<apf::Vector3, double> three_pt_circ(apf::Vector3& p0, apf::Vector3& p1,
                                              apf::Vector3& p2) {
  apf::Vector3 t = p1 - p0;
  apf::Vector3 u = p2 - p0;
  apf::Vector3 v = p2 - p1;

  // triangle normal
  apf::Vector3 w = vd_cross(t, u);
  double wsl = w.getLength();

  // helpers ??
  double iwsl2 = 1.0 / (2.0*wsl);
  double tt = t*t;
  double uu = u*u;

  // result circle
  apf::Vector3 circCenter = p1 + (u*tt*(u*v) - t*uu*(t*v)) * iwsl2;
  double circRadius = sqrt(tt * uu * (v*v) * iwsl2*0.5);
  apf::Vector3 circAxis = w / sqrt(wsl);
  return std::make_pair(circCenter, circRadius);
}
*/
