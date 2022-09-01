#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <ma.h>

#include "topo_manip.h"
#include "topo_field.h"
#include "topo_rand.h"
#include "topo_extinfo.h"
#include "topo_geom.h"
#include "topo_solvlin.h"

#include <algorithm>
#include <math.h>       /* acos, isnan */

#include <limits>       /* numeric_limits */

int lookup_tri_e [3][2] = {{0,2},{0,1},{1,2}};
//int lookup_tri_vd [3][2] = {{1,2},{0,2},{1,0}};
int lookup_tri_vd [3][2] = {{2,1},{0,2},{1,0}};
int lookup_v_tri [4] = {2,3,1,0};

apf::Matrix3x3 const M_EYE(1,0,0,0,1,0,0,0,1);

// Perturb the position of a point, without inverting elements.
// Poor practice, but works for good quality tetrahedrons.
void vd_pert_point(apf::Mesh2* m, apf::MeshEntity* v1) {

  apf::Up up;
  apf::MeshEntity* edge;
  m->getUp(v1, up);
  int upwardCount = up.n;
  float edge_min = 10.0;
  int edge_index;

  apf::MeshEntity* v2;
  apf::Vector3 point1(0,0,0);
  apf::Vector3 point2(0,0,0);

  m->getPoint(v1, 0, point1);
  //printf("The position of the vertex1: (%1.2f, %1.2f, %1.2f).\n", point1[0], point1[1], point1[2]);
  for (int i = 0; i < upwardCount; i++) {
    edge = up.e[i];
    v2 = apf::getEdgeVertOppositeVert(m, edge, v1);
    m->getPoint(v2, 0, point2);
    //printf("The position of the vertex2: (%1.2f, %1.2f, %1.2f).\n", point2[0], point2[1], point2[2]);
    point1 = point2 - point1;
    if ( point1.getLength() < edge_min) {
      edge_min = point1.getLength();
      edge_index = i;
    }
    // m->setPoint(v1, 1, point1);
  }
  apf::Vector3 r(0,0,0);
  vd_rand_vect3(r);
  printf("Edgemin is %1.2f.\n", edge_min);
  r = r * edge_min/10.1;
  printf("The displacement applied to vertex1: (%1.2f, %1.2f, %1.2f).\n", r[0], r[1], r[2]);
  m->getPoint(v1, 0, point1);
  m->setPoint(v1, 0, point1 + r);

}

// Perturb the positions of all points in a mesh, without inverting elements.
void vd_pert_mesh(apf::Mesh2* m) {

  int d = 0;
  int ID = 1;
  apf::MeshEntity* v1;
  apf::MeshIterator* it = m->begin(d);

  while ((v1 = m->iterate(it))) {
    vd_pert_point(m, v1);
    ID++;
  }

  m->end(it);
}

// Perturb the position of a point, by the ratio of the maximum distance along a 
// random direction.
void vd_pert_point_rat(apf::Mesh2* m, apf::MeshEntity* v, double rat) {

  apf::Vector3 r(0,0,0);
  apf::Vector3 point(0,0,0);
  vd_rand_vect3(r);
  r = norm_0(r);
  double dist = vd_dist_v_x(m, v, r);

  m->getPoint(v, 0, point);
  m->setPoint(v, 0, point + r*rat*dist);

}

// Perturb the position of all points, by the ratio of the maximum distance along 
// a random direction.
void vd_pert_mesh_rat(apf::Mesh2* m, double rat) {
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    vd_pert_point_rat(m, v, rat);
  }
  m->end(it);
}

// Perturb the position of a set of points, by the ratio of the maximum distance 
// along a random direction.
void vd_pert_set_rat(apf::Mesh2* m, std::vector<apf::MeshEntity*>& v, double rat) {
  for(int i = 0; i < v.size(); i++) {
    vd_pert_point_rat(m, v.at(i), rat);
  }
}


// Given a mesh and axis, shrink the mesh. Assumes shrinking around 0,0,0
void vd_shr_axis(apf::Mesh2* m, apf::Vector3 ax, double sp_size) {
  apf::Vector3 point(0,0,0);
  double proj;

  ax = norm_0(ax);

  apf::MeshEntity* v;
  apf::MeshIterator* v_it = m->begin(0);

  while ((v = m->iterate(v_it))) {
    m->getPoint(v, 0, point);
    proj = ax*point;
    point = point+ax*proj*(sp_size-1);
    m->setPoint(v, 0, point);
  }

  m->end(v_it);
}

// Given a mesh and axis, shrink the mesh. Assumes shrinking around 0,0,0
void vd_shr_axis2(apf::Mesh2* m, apf::Vector3 ax, double sp_size) {
  apf::Vector3 point(0,0,0);
  apf::Vector3 disp(0,0,0);
  apf::Vector3 p_zero(0,0,0);

  double proj;

  ax = norm_0(ax);

  apf::MeshEntity* v;
  apf::MeshIterator* v_it = m->begin(0);

  while ((v = m->iterate(v_it))) {
    m->getPoint(v, 0, point);
    disp = (point-p_zero);

    proj = ax*disp;
    point = (disp-ax*proj)+ax*proj*sp_size+p_zero;
    m->setPoint(v, 0, point);
  }

  m->end(v_it);
}

// Given a mesh and axis, shrink the mesh.
apf::Vector3 vd_shr_vec(apf::Vector3 vec, apf::Vector3 ax, double sp_size) {
  ax = ax*(vec*ax);
  return (vec-ax)+ax*sp_size;
}

// TODO: Given a mesh and surface geometry, shrink the surface.
void vd_shr_surf(apf::Mesh2* m, int geom, float sp_size) {
  std::vector<apf::MeshEntity*> surf(0);
  std::vector<apf::MeshEntity*> vert(0);
  std::vector<apf::MeshEntity*> quad(0);

  vd_find_ent_geom(m, &surf, geom, 2);
  vd_set_down(m, &surf, &vert, 2);
  vd_ext_ent_geom(m, &vert, &quad, -1);

  apf::Vector3 surf_center(0,0,0);
  surf_center = vd_get_center(m, &quad);

  // Up to some point, this part is same with pert_sphere, so create another
  // function.
  apf::Vector3 point(0,0,0);
  apf::Vector3 disp(0,0,0);

  float distance = 0;
  for (int i = 0; i < quad.size(); i++) {
    m->getPoint(quad.at(i), 0, point);
    disp = (point-surf_center);
    if (disp.getLength() > distance)
      distance = disp.getLength();
  }

  // Find the scale to displace the central vertices:
  float scale = sp_size/distance;

  // Displace the domain vertices:
  for (int i = 0; i < vert.size(); i ++) {
    m->getPoint(vert.at(i), 0, point);    
    disp = (point-surf_center);
    distance = disp.getLength();

    m->setPoint(vert.at(i), 0, surf_center + disp*sp_size/distance);
  }

}

// Given a mesh and cell, shrink the cell by making the vertices belonging to 
// the cell and its lower dimensional adjacencies closer to the cell center 
// defined as the average position of all belonging vertices.
void vd_shr_cell(apf::Mesh2* m, int dim, int cell, float sp_size) {
  std::vector<apf::MeshEntity*> vert(0);
  std::vector<apf::MeshEntity*> es_cell(0);

  vd_find_ent_geom(m, &es_cell, cell, dim);

  double cell_meas = vd_meas_geom(m, dim, cell);
  printf("Size before shrinking %.5f... ", cell_meas);
  vd_set_down(m, &es_cell, &vert, dim);

  apf::Vector3 cell_center(0,0,0);
  cell_center = vd_get_center(m, &vert);

  // Up to some point, this part is same with pert_sphere, so create another
  // function.
  apf::Vector3 point(0,0,0);
  apf::Vector3 disp(0,0,0);

  float distance = 0;
  for (int i = 0; i < vert.size(); i++) {
    m->getPoint(vert.at(i), 0, point);
    disp = (point-cell_center);
    if (disp.getLength() > distance)
      distance = disp.getLength();
  }

  // Find the scale to displace the central vertices:
  float scale = sp_size/distance;

  // Displace the domain vertices:
  for (int i = 0; i < vert.size(); i ++) {
    m->getPoint(vert.at(i), 0, point);    
    disp = (point-cell_center);

    m->setPoint(vert.at(i), 0, cell_center + disp*sp_size);
  }

  cell_meas = vd_meas_geom(m, dim, cell);
  printf("Afterwards: %.5f. \n", cell_meas);

}

// Given a cell, try to make the cell shape equiaxed. This is achieved by
// rescaling the positions such that the difference between the maximum and 
// minimum coordinates x,y,z are made the same. It can be also achieved by
// PCA type of analysis.
bool vd_pert_cell_equi(apf::Mesh2* m, int cell_dim, int cell_id) {
  if (cell_dim > 0) {
    std::vector<apf::MeshEntity*> entit(0);
    std::vector<apf::MeshEntity*> es_cell_ext(0);
    std::vector<apf::MeshEntity*> es_cell_int(0);

    vd_find_ent_geom(m, &entit, cell_id, cell_dim);

    if(entit.size() > 0) {
      if (cell_dim == 3) {
        vd_set_down(m, &entit, &es_cell_ext, 3);
      }

      if (cell_dim == 2) {
        vd_set_down(m, &entit, &es_cell_ext, 2);
      }

      if (cell_dim == 1) {
        vd_set_down(m, &entit, &es_cell_ext, 1);
      }

      es_cell_int.reserve(es_cell_ext.size());

      for (int i = 0; i < es_cell_ext.size(); i++) {
        apf::ModelEntity* vert_mdl = m->toModel (es_cell_ext.at(i));
        int vert_type = m->getModelType (vert_mdl);
        int vert_tag = m->getModelTag (vert_mdl);
        if (vert_type == cell_dim and vert_tag == cell_id) {
          es_cell_int.push_back(es_cell_ext.at(i));
        }
      }

      vd_remove_set(&es_cell_ext, &es_cell_int);

      apf::Vector3 position(0,0,0);
      m->getPoint(es_cell_ext.at(0), 0, position);

      double coords[6]; // x_min x_max y_min y_max z_min z_max
      coords[0] = position[0]; coords[1] = position[0];
      coords[2] = position[1]; coords[3] = position[1];
      coords[4] = position[2]; coords[5] = position[2];

      std::cout << "Position is " << position << std::endl;

      for (int i = 1; i < es_cell_ext.size(); i++) {
        m->getPoint(es_cell_ext.at(i), 0, position);
        std::cout << "Position is " << position << std::endl;
        if (position[0] < coords[0])
          coords[0] = position[0];
        if (position[0] > coords[1])
          coords[1] = position[0];
        if (position[1] < coords[2])
          coords[2] = position[1];
        if (position[1] > coords[3])
          coords[3] = position[1];
        if (position[2] < coords[4])
          coords[4] = position[2];
        if (position[2] > coords[5])
          coords[5] = position[2];
      }

      printf("The min max coord are (%.3f-%.3f,%.3f-%.3f,%.3f-%.3f). \n",
                 coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);

      double scale[3]; // x_scale y_scale z_scale
      double diff[3]; // x_scale y_scale z_scale

      diff[0] = coords[1] - coords[0];
      diff[1] = coords[3] - coords[2];
      diff[2] = coords[5] - coords[4];

      if (diff[0] <= 0.) {
        diff[0] = 0;
      }
      if (diff[1] <= 0.) {
        diff[1] = 0;
      }
      if (diff[2] <= 0.) {
        diff[2] = 0;
      }

      double diff_mean = (diff[0]+diff[1]+diff[2])/3;

      scale[0] = diff[0]/diff_mean;
      scale[1] = diff[1]/diff_mean;
      scale[2] = diff[2]/diff_mean;

      if (diff[0] <= 0.) {
        scale[0] = 1;
      }
      if (diff[1] <= 0.) {
        scale[1] = 1;
      }
      if (diff[2] <= 0.) {
        scale[2] = 1;
      }

      //std::cout << "Scale are before " << scale[0] << ", " <<
      //                                   scale[1] << ", " <<
      //                                   scale[2] << std::endl;

      // Further reduce in order to prevent inversions.
      double scale_max = 0;
      for (int i = 0; i < 3; i++) {
        if (scale[i] > scale_max)
          scale_max = scale[i];
      }

      scale [0] = scale[0]/scale_max;
      scale [1] = scale[1]/scale_max;
      scale [2] = scale[2]/scale_max;

      //std::cout << "Scale are " << scale[0] << ", " <<
      //                             scale[1] << ", " <<
      //                             scale[2] << std::endl;

      apf::Vector3 cell_center(0,0,0);
      cell_center = vd_get_center(m, &es_cell_ext);
      apf::Vector3 disp(0,0,0);

      // Before changing the positions of the vertices, store them. If the changes 
      // introduce element inversion, revert the changes and return false.
      // If no inversion occurs, return true.

      std::vector<apf::Vector3> v_pos_e(0, apf::Vector3(0,0,0));
      std::vector<apf::Vector3> v_pos_i(0, apf::Vector3(0,0,0));
      v_pos_e.resize(es_cell_ext.size());
      v_pos_i.resize(es_cell_int.size());

      for (int i = 0; i < es_cell_ext.size(); i++) {
        m->getPoint(es_cell_ext.at(i), 0, position);
        v_pos_e.at(i) = position;
        //std::cout << "Position before " << position << std::endl;
        disp = position - cell_center;
        double vect[3];
        disp.toArray(vect);
        vect[0] = vect[0]/scale[0];
        vect[1] = vect[1]/scale[1];
        vect[2] = vect[2]/scale[2];
        disp.fromArray(vect);
        position = cell_center+disp;
        m->setPoint(es_cell_ext.at(i), 0, position);
        //std::cout << "Position after " << position << std::endl;
      }

      for (int i = 0; i < es_cell_int.size(); i++) {
        m->getPoint(es_cell_int.at(i), 0, position);
        v_pos_i.at(i) = position;
        //std::cout << "Position before " << position << std::endl;
        disp = position - cell_center;
        double vect[3];
        disp.toArray(vect);
        vect[0] = vect[0]/scale[0];
        vect[1] = vect[1]/scale[1];
        vect[2] = vect[2]/scale[2];
        disp.fromArray(vect);
        position = cell_center+disp;
        m->setPoint(es_cell_int.at(i), 0, position);
        //std::cout << "Position after " << position << std::endl;
      }
      if (vd_chk_neg(m) > 0) {
        for (int i = 0; i < es_cell_ext.size(); i++) {
          m->setPoint(es_cell_ext.at(i), 0, v_pos_e.at(i));
        }
        for (int i = 0; i < es_cell_int.size(); i++) {
          m->setPoint(es_cell_int.at(i), 0, v_pos_i.at(i));
        }

        return false;
      }

    }
    return true;
  }
  return false;
}

// Given a cell, project the outer vertices onto a sphere of mean width of the
// cell.
void vd_pert_cell_sphere(apf::Mesh2* m, int cell_dim, int cell_id) {
  std::vector<apf::MeshEntity*> entit(0);
  std::vector<apf::MeshEntity*> es_cell_ext(0);
  std::vector<apf::MeshEntity*> es_cell_int(0);

  vd_find_ent_geom(m, &entit, cell_id, cell_dim);

  if (cell_dim == 3) {
    vd_set_down(m, &entit, &es_cell_ext, 3);
  }

  if (cell_dim == 2) {
    vd_set_down(m, &entit, &es_cell_ext, 2);
  }

  if (cell_dim == 1) {
    vd_set_down(m, &entit, &es_cell_ext, 1);
  }

  es_cell_int.reserve(es_cell_ext.size());
  for (int i = 0; i < es_cell_ext.size(); i++) {
    apf::ModelEntity* vert_mdl = m->toModel (es_cell_ext.at(i));
    int vert_type = m->getModelType (vert_mdl);
    int vert_tag = m->getModelTag (vert_mdl);
    if (vert_type == cell_dim and vert_tag == cell_id) {
      es_cell_int.push_back(es_cell_ext.at(i));
    }
  }

  vd_remove_set(&es_cell_ext, &es_cell_int);

  printf("%dcell%d external vert:%d, internal vert:%d\n", cell_dim, cell_id, 
                                              (int)(es_cell_ext.size()), 
                                              (int)(es_cell_int.size()) );

  apf::Vector3 cell_center(0,0,0);
  cell_center = vd_get_center(m, &es_cell_ext);
  float radius = 0;

  apf::Vector3 position(0,0,0);
  apf::Vector3 disp(0,0,0);

  for (int i = 0; i < es_cell_ext.size(); i++) {
    m->getPoint(es_cell_ext.at(i), 0, position);
    disp = position - cell_center;
    radius = radius + sqrt(disp*disp);
  }
  radius = radius/es_cell_ext.size();

  for (int i = 0; i < es_cell_ext.size(); i++) {
    m->getPoint(es_cell_ext.at(i), 0, position);
    disp = position - cell_center;
    float coef = radius/sqrt(disp*disp);
    position = cell_center + disp*coef;
    m->setPoint(es_cell_ext.at(i), 0, position);
  }

  vd_rem_tag(m);
  vd_tag_mesh(m);
  apf::writeVtkFiles("./output/sphere", m);

  float disp_max = 0;
  for (int i = 0; i < es_cell_int.size(); i++) {
    float disp_curr = 0;
    m->getPoint(es_cell_int.at(i), 0, position);
    disp = position - cell_center;
    disp_curr = sqrt(disp*disp);
    if (disp_max < disp_curr) {
      disp_max = disp_curr;
    }
  }

  for (int i = 0; i < es_cell_int.size(); i++) {
    m->getPoint(es_cell_int.at(i), 0, position);
    disp = position - cell_center;
    float coef = radius/(disp_max*2);
    position = cell_center + disp*coef;
    m->setPoint(es_cell_ext.at(i), 0, position);
  }

  vd_rem_tag(m);
  vd_tag_mesh(m);
  apf::writeVtkFiles("./output/sphere2", m);
}

// Given a mesh, project the outer surface points onto a sphere of size sp_size
// and center, and scale other points considering the farthest domain surface 
// point.
void vd_pert_sphere(apf::Mesh2* m, float sp_size) {

  // Extract the domain geometric vertices:
  // Find the surface and domain surface triangle elements:
  std::vector<apf::MeshEntity*> es_gsurf(0);
  std::vector<apf::MeshEntity*> es_dom(0);

  vd_ext_bsurf(m, &es_gsurf, &es_dom);  

  // Find the vertices belonging to domain surface triangles:
  std::vector<apf::MeshEntity*> vert_dom(0);
  std::vector<apf::MeshEntity*> vert_dom_quad(0);
  vd_set_down(m, &es_dom, &vert_dom, 2);
  vd_ext_ent_geom(m, &vert_dom, &vert_dom_quad);

  // Find the center, and for the cube corners, maximum displacement from 
  // the center.
  apf::Vector3 center(0,0,0);
  center = vd_get_center(m, &vert_dom_quad);

  apf::Vector3 point(0,0,0);
  apf::Vector3 disp(0,0,0);

  float distance = 0;
  for (int i = 0; i < vert_dom_quad.size(); i++) {
    m->getPoint(vert_dom_quad.at(i), 0, point);
    disp = (point-center);
    if (disp.getLength() > distance)
      distance = disp.getLength();
  }

  // Find the scale to displace the central vertices:
  float scale = sp_size/distance;

  // Displace the domain vertices:
  for (int i = 0; i < vert_dom.size(); i ++) {
    m->getPoint(vert_dom.at(i), 0, point);    
    disp = (point-center);
    distance = disp.getLength();

    m->setPoint(vert_dom.at(i), 0, center + disp*sp_size/distance);
  }
  // Displace the internal vertices:
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;

  while ((vert = m->iterate(it))) {
    if (findIn(&vert_dom, vert_dom.size(), vert) == -1) {
      m->getPoint(vert, 0, point);    
      disp = (point-center);
      m->setPoint(vert, 0, center + disp*scale);
    }
  }

  m->end(it);

}

void vd_proj_v_sphere(apf::Mesh2* m, apf::MeshEntity* vert, 
                          apf::Vector3 ctr, double rho){
  apf::Vector3 point(0,0,0);
  m->getPoint(vert, 0, point);
  point = (point-ctr);
  if(point.getLength() <= std::numeric_limits<double>::min())
    point = point*0;
  else
    point = norm_0(point)*rho;

  m->setPoint(vert, 0, ctr + point);
}

// Destroy the mesh entities inside an entity set.
void vd_dest_set(apf::Mesh2* m, std::vector<apf::MeshEntity*>* set_dest) {
  for (int i = 0; i < set_dest->size(); i++) {
    m->destroy(set_dest->at(i));
  }
}

// Tag the entities above the cutting plane.
void vd_cut_mesh(apf::Mesh2* m, plane* cut_plane) {
  //double eps = 0.0001;

  // Find vertices above cutting plane, find upper adjacencies.
  // Destroy every entity.
  std::vector<apf::MeshEntity*> vert_cut(0);

  apf::MeshEntity* vert;

  int count = 0;
  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    apf::Vector3 point(0,0,0);
    m->getPoint(vert, 0, point);    
    //std::cout<< "Current position" << point << "plane parameter: " << point*cut_plane->norm <<std::endl;

    if ( cut_plane->a < point*cut_plane->norm ) {
      count = count + 1;
    }
  }
  m->end(it);

  vert_cut.reserve(count);

  it = m->begin(0);
  while ((vert = m->iterate(it))) {
    apf::Vector3 point(0,0,0);
    m->getPoint(vert, 0, point);    
    //std::cout<< "Current position" << point << "plane parameter: " << point*cut_plane->norm <<std::endl;

    if ( cut_plane->a < point*cut_plane->norm ) {
      vert_cut.push_back(vert);
    }
  }

  m->end(it);

  std::vector<apf::MeshEntity*> edge_cut(0);
  std::vector<apf::MeshEntity*> surf_cut(0);
  std::vector<apf::MeshEntity*> elem_cut(0);

  vd_set_up(m, &vert_cut, &edge_cut);
  vd_set_up(m, &edge_cut, &surf_cut);
  vd_set_up(m, &surf_cut, &elem_cut);

  printf("There are %d vertices, \n", (int)(vert_cut.size()));
  printf("There are %d edges, \n", (int)(edge_cut.size()));
  printf("There are %d triangles, \n", (int)(surf_cut.size()));
  printf("There are %d elements, above the cut surface. \n", 
                                               (int)(elem_cut.size()));

  vd_dest_set(m, &elem_cut);
  vd_dest_set(m, &surf_cut);
  vd_dest_set(m, &edge_cut);
  vd_dest_set(m, &vert_cut);
}


// Tag the entities above the cutting plane.
void vd_tag_cut(apf::Mesh2* m, plane* cut_plane) {

  // Find vertices above cutting plane, find upper adjacencies.
  // Destroy every entity.
  std::vector<apf::MeshEntity*> vert_cut(0);


  apf::MeshEntity* vert;
  int count = 0;

  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    apf::Vector3 point(0,0,0);
    m->getPoint(vert, 0, point);    
    //std::cout<< "Current position" << point << "plane parameter: " << point*cut_plane->norm <<std::endl;

    if ( cut_plane->a < point*cut_plane->norm ) {
      count = count + 1;
    }
  }
  m->end(it);

  vert_cut.reserve(count);

  it = m->begin(0);
  while ((vert = m->iterate(it))) {
    apf::Vector3 point(0,0,0);
    m->getPoint(vert, 0, point);    
    //std::cout<< "Current position" << point << "plane parameter: " << point*cut_plane->norm <<std::endl;

    if ( cut_plane->a < point*cut_plane->norm ) {
      vert_cut.push_back(vert);
    }
  }

  m->end(it);

  printf("%d vertex entities are cut.\n", (int)(vert_cut.size()));

  std::vector<apf::MeshEntity*> edge_cut(0);
  std::vector<apf::MeshEntity*> surf_cut(0);
  std::vector<apf::MeshEntity*> elem_cut(0);

  vd_set_up(m, &vert_cut, &edge_cut);
  printf("%d edge entities are cut.\n", (int)(edge_cut.size()));
  vd_set_up(m, &edge_cut, &surf_cut);
  printf("%d surf entities are cut.\n", (int)(surf_cut.size()));
  vd_set_up(m, &surf_cut, &elem_cut);
  printf("%d elem entities are cut.\n", (int)(elem_cut.size()));

  vd_tag_set(m, &elem_cut, "Cutelem");
  //vd_tag_set(m, &surf_cut, "Cutsurf");
  //vd_tag_set(m, &edge_cut, "Cutedge");
  vd_tag_set(m, &vert_cut, "Cutvert");

}

// Tag the entities above the cutting plane.
void vd_cell_cut(apf::Mesh2* m, plane* cut_plane) {

  struct gmi_model* mdl = m->getModel();
  struct gmi_iter* g_it = gmi_begin(mdl, 3);
  struct gmi_ent* g_ent;

  int nbr_cell = 0;
  while ((g_ent = gmi_next(mdl, g_it)))
    nbr_cell++;
  gmi_end(mdl, g_it);

  //std::cout << "cut plane: " << cut_plane->a << " = "
  //          << cut_plane->norm[0] << "x + "
  //          << cut_plane->norm[1] << "y + "
  //          << cut_plane->norm[2] << "z"
  //          << std::endl;
  // If any cell has any element below the plane, it is not drawn.
  std::vector<int> cell_draw(nbr_cell);
  for(int i = 0; i < nbr_cell; i++)
    cell_draw.at(i) = 0;
  printf("Number of cells: %d", nbr_cell);

  //apf::Numbering* tetnumbering;
  //if (m->findNumbering("Cuttet")) {
  //  std::cout<< "Cuttet" << " tag already exists" << std::endl;
  //  tetnumbering = m->findNumbering("Cuttet");
  //}
  //else {
  //  tetnumbering = apf::createNumbering(m, "Cuttet", apf::getConstant(3),1);
  //}
  // Find elements with their centroids above the cutting plane.

  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* elem;

  while ((elem = m->iterate(it))) {
    apf::Vector3 point(0,0,0);
    point = getLinearCentroid(m, elem);
    //std::cout << "p: " << point 
    //          << " distance " << point*cut_plane->norm
    //          << std::endl;
    if ( cut_plane->a < point*cut_plane->norm ) {
      int tag = m->getModelTag(m->toModel(elem));
      cell_draw.at(tag - 1) = 1;
      //apf::number(tetnumbering, elem, 0, 0, 0);
    }
    //else
      //apf::number(tetnumbering, elem, 0, 0, 1);
  }
  m->end(it);

  PCU_Barrier();
  PCU_Max_Ints(&cell_draw.front(), nbr_cell);

  //for(int i = 0; i < nbr_cell; i++)
  //  printf("Cell%d: %d\n", i, cell_draw[i]);

/*
  Entity_set elem_cut;

  it = m->begin(3);
  while ((elem = m->iterate(it))) {
    int tag = m->getModelTag(m->toModel(elem));
    //printf("Elem %p, tag: %d\n", elem, tag);
    if (cell_draw[tag] == 1) {
      elem_cut.at(elem_cut.n) = elem;
      elem_cut.n++;
    }
  }
  m->end(it);

  vd_tag_set(m, &elem_cut, "Cutcell");
*/
  apf::Numbering* setnumbering;
  if (m->findNumbering("Cutcell")) {
    std::cout<< "Cutcell" << " tag already exists" << std::endl;
    setnumbering = m->findNumbering("Cutcell");
  }
  else {
    setnumbering = apf::createNumbering(m, "Cutcell", apf::getConstant(3),1);
  }

  it = m->begin(3);

  std::cout<< "Writing tag " << "Cutcell" << "...";
  while ((elem = m->iterate(it))) 
  {
    int tag = m->getModelTag(m->toModel(elem));
    if (! isNumbered(setnumbering, elem, 0, 0)) {
      if (cell_draw.at(tag - 1) == 1) {
        apf::number(setnumbering, elem, 0, 0, 0);
      }
      else {
        apf::number(setnumbering, elem, 0, 0, 1);
      }
    }
  }
  m->end(it);
  apf::writeVtkFiles("./output/cut_plane", m);

}

// Tag the entities above the cutting plane.
void vd_cell_cut_sph(apf::Mesh2* m, sphere* cut_sph) {

  struct gmi_model* mdl = m->getModel();
  struct gmi_iter* g_it = gmi_begin(mdl, 3);
  struct gmi_ent* g_ent;

  int nbr_cell = 0;
  while ((g_ent = gmi_next(mdl, g_it)))
    nbr_cell++;
  gmi_end(mdl, g_it);

  //std::cout << "cut plane: " << cut_plane->a << " = "
  //          << cut_plane->norm[0] << "x + "
  //          << cut_plane->norm[1] << "y + "
  //          << cut_plane->norm[2] << "z"
  //          << std::endl;
  // If any cell has any element below the plane, it is not drawn.
  std::vector<int> cell_draw(nbr_cell);
  for(int i = 0; i < nbr_cell; i++)
    cell_draw.at(i) = 0;
  printf("Number of cells: %d", nbr_cell);

  //apf::Numbering* tetnumbering;
  //if (m->findNumbering("Cuttet")) {
  //  std::cout<< "Cuttet" << " tag already exists" << std::endl;
  //  tetnumbering = m->findNumbering("Cuttet");
  //}
  //else {
  //  tetnumbering = apf::createNumbering(m, "Cuttet", apf::getConstant(3),1);
  //}
  // Find elements with their centroids above the cutting plane.

  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* elem;

  while ((elem = m->iterate(it))) {
    apf::Vector3 point(0,0,0);
    point = getLinearCentroid(m, elem);
    //std::cout << "p: " << point 
    //          << " distance " << point*cut_plane->norm
    //          << std::endl;
    point = point - cut_sph->ctr;
    if ( cut_sph->r > point.getLength() ) {
      int tag = m->getModelTag(m->toModel(elem));
      cell_draw.at(tag - 1) = 1;
      //apf::number(tetnumbering, elem, 0, 0, 0);
    }
    //else
      //apf::number(tetnumbering, elem, 0, 0, 1);
  }
  m->end(it);

  PCU_Barrier();
  PCU_Max_Ints(&cell_draw.front(), nbr_cell);

  //for(int i = 0; i < nbr_cell; i++)
  //  printf("Cell%d: %d\n", i, cell_draw[i]);

/*
  Entity_set elem_cut;

  it = m->begin(3);
  while ((elem = m->iterate(it))) {
    int tag = m->getModelTag(m->toModel(elem));
    //printf("Elem %p, tag: %d\n", elem, tag);
    if (cell_draw[tag] == 1) {
      elem_cut.at(elem_cut.n) = elem;
      elem_cut.n++;
    }
  }
  m->end(it);

  vd_tag_set(m, &elem_cut, "Cutcell");
*/
  apf::Numbering* setnumbering;
  if (m->findNumbering("Cutcell")) {
    std::cout<< "Cutcell" << " tag already exists" << std::endl;
    setnumbering = m->findNumbering("Cutcell");
  }
  else {
    setnumbering = apf::createNumbering(m, "Cutcell", apf::getConstant(3),1);
  }

  it = m->begin(3);

  std::cout<< "Writing tag " << "Cutcell" << "...";
  while ((elem = m->iterate(it))) 
  {
    int tag = m->getModelTag(m->toModel(elem));
    if (! isNumbered(setnumbering, elem, 0, 0)) {
      if (cell_draw.at(tag - 1) == 1) {
        apf::number(setnumbering, elem, 0, 0, 0);
      }
      else {
        apf::number(setnumbering, elem, 0, 0, 1);
      }
    }
  }
  m->end(it);
  apf::writeVtkFiles("./output/cut_sph", m);

}

//////////////////////////////////////
// Calculate velocity field:
//////////////////////////////////////
//////////////////////////////////////
// vd_eqn_of_motion
//////////////////////////////////////

void vd_eqn_of_motion::calc_vel() {}
void vd_eqn_of_motion::calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local) {}
void vd_eqn_of_motion::calc_vel_curr(apf::MeshEntity* vert) {}
// A single vertex.
void vd_eqn_of_motion::vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local) {}
apf::Vector3 vd_eqn_of_motion::vd_upd_vel_field_tri(apf::MeshEntity* vert, 
        std::vector<apf::MeshEntity*>* tris, bool drag_local) {
  return apf::Vector3(0,0,0);
}
apf::Vector3 vd_eqn_of_motion::vd_upd_vel_field_edge(apf::MeshEntity* vert, 
           std::vector<apf::MeshEntity*>* edges, bool drag_local) {
  return apf::Vector3(0,0,0);
}

// A set of vertices while skipping labeled vertices.
void vd_eqn_of_motion::vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                       bool drag_local) {
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_eqn_of_motion::vd_calc_force(apf::MeshEntity* vert) {
  return apf::Vector3(0,0,0);
}

apf::Vector3 vd_eqn_of_motion::vd_calc_force_tri(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris) {
  return apf::Vector3(0,0,0);
}

apf::Vector3 vd_eqn_of_motion::vd_calc_force_edge(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges) {
  return apf::Vector3(0,0,0);
}


vd_eqn_of_motion::vd_eqn_of_motion() : m(NULL), c_base(NULL), f_calc(NULL), 
                                       e_list(NULL),
                                       vel_field(NULL), drag_rat(1000) {
}

vd_eqn_of_motion::vd_eqn_of_motion(apf::Mesh2* m_in, cell_base* c_base_in, 
                            field_calc* f_calc_in, 
                            vd_entlist* e_list_in) : 
                    m(m_in), f_calc(f_calc_in), e_list(e_list_in),
                                       vel_field(NULL), drag_rat(1000) {
}

vd_eqn_of_motion::~vd_eqn_of_motion() {
}

void vd_eqn_of_motion::refresh(apf::Mesh2* m_in, cell_base* c_base_in, field_calc* f_calc_in, vd_entlist* e_list_in) {
  m = m_in;
  f_calc = f_calc_in;
  e_list = e_list_in;
  c_base = c_base_in;
}

void vd_eqn_of_motion::refresh_elist() {
  e_list->refresh();
}


vd_eqn_of_motion::vd_eqn_of_motion(const vd_eqn_of_motion& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  vel_field = that.vel_field;

  c_base = that.c_base;
}

vd_eqn_of_motion& vd_eqn_of_motion::operator=(const vd_eqn_of_motion& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  vel_field = that.vel_field;

  c_base = that.c_base;
  return *this;
}

void vd_eqn_of_motion::set_drag_rat(double rat_in) {
  drag_rat = rat_in;
}

//////////////////////////////////////
// vd_eqn_mason
//////////////////////////////////////
// Update the velocity field at:
// Every boundary vertex.


vd_eqn_mason::vd_eqn_mason() : vd_eqn_of_motion(),
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
}

vd_eqn_mason::vd_eqn_mason(apf::Mesh2* m_in, cell_base* c_base_in,
                            field_calc* f_calc_in, 
                            vd_entlist* e_list_in) :
                         vd_eqn_of_motion(m_in, c_base_in, f_calc_in, e_list_in), 
                           es_edge(0), es_surf(0), 
                           p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                           p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                           average(-1.), drag_rat(1000) {
  m = m_in;
  f_calc = f_calc_in;
  e_list = e_list_in;
  c_base = c_base_in;
}

vd_eqn_mason::vd_eqn_mason(const vd_eqn_mason& that) :
                               vd_eqn_of_motion(), 
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;
}

vd_eqn_mason& vd_eqn_mason::operator=(const vd_eqn_mason& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;

  es_edge.clear();
  es_surf.clear();
  apf::Vector3 zero(0,0,0);
  p_i = zero;
  p_j = zero;
  n_ij = zero;
  p_ctr = zero;
  temp = zero;
  rhs = zero;
  average = that.average;
  vel_field = that.vel_field;
  drag_rat = that.drag_rat; 
  return *this;
}

void vd_eqn_mason::get_average_tri() {
  average = f_calc->get_drag_glob();
  if(average < - std::numeric_limits<double>::min())
    average = getAverageEntSize(m, 2);
}

void vd_eqn_mason::calc_vel_curr(apf::MeshEntity* vert) {

  apf::Vector3 zero(0,0,0);
  apf::Vector3 area_proj(0,0,0);
  apf::Vector3 f_temp(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);
  //double norm_sum = 0;
  //double len_avg = 0;
  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* em = m->toModel(es_surf.at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      //area_proj = area_proj + n_ij;

      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      //norm_sum = norm_sum + norm;
      //len_avg = len_avg + (p_i.getLength() + p_j.getLength())/2;
      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      f_temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, es_surf.at(i));
      rhs = rhs + f_temp;
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, es_surf.at(i));
    }
  }
  //len_avg = len_avg/es_surf.size();
  //area_proj = area_proj/6;

  m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;
  // 2017ActeMateMason:
  // rhs = F*2
  // m_eij = D*6
  // apf::invert(m_eij) : Mobility tensor/6.
  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
}

apf::Vector3 vd_eqn_mason::calc_vel_curr_tri(apf::MeshEntity* vert,
                           std::vector<apf::MeshEntity*>* tris) {

  apf::Vector3 zero(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < tris->size(); i++) {
    apf::ModelEntity* em = m->toModel(tris->at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                            f_calc->gam2(m, tris->at(i));
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, tris->at(i));
    }
  }
  m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;

  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
  return temp;
}

void vd_eqn_mason::calc_vel() {
  vel_field = m->findField("velocity_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      calc_vel_curr(vert);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  m->end(it);
  printf("Velocity field finished.\n");
}

void vd_eqn_mason::calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local) {
  vel_field = m->findField("velocity_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  for(int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      vd_upd_vel_field(vert, drag_local);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  printf("Velocity field finished.\n");
}

// A single vertex.
void vd_eqn_mason::vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  if(drag_local) {
    average = getAverageEntSize(m, vert, 2);
  }
  else
    get_average_tri();
  calc_vel_curr(vert);
}

// A single vertex.
apf::Vector3 vd_eqn_mason::vd_upd_vel_field_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  if(drag_local) {
    average = getAverageEntSize(m, vert, 2);
  }
  else
    get_average_tri();
  return calc_vel_curr_tri(vert, tris);
}

apf::Vector3 vd_eqn_mason::vd_upd_vel_field_edge(apf::MeshEntity* vert, 
           std::vector<apf::MeshEntity*>* edges, bool drag_local) {
  std::cout << "upd_vel_field_edge not defined!" << std::endl;
  return apf::Vector3(0,0,0);
}

// A set of vertices while skipping labeled vertices. Merging vertices contribute
// to the target vertex and they all have the same resultant force and velocity.
void vd_eqn_mason::vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                       bool drag_local) {
  vel_field = m->findField("velocity_field");

  std::map<apf::MeshEntity*, bool> merg{};
  apf::MeshEntity* v_merg = NULL;
  for(int i = 0; i < vert.size(); i++) {
    if(v_map[vert.at(i)] != NULL) {
      merg[vert.at(i)] = true;
      if(v_merg == NULL)
        v_merg = v_map[vert.at(i)];
    }
  }
  assert(v_merg != NULL);

  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find(vert.begin(), vert.end(), v_merg);
  assert(it != vert.end());
  int i_m = std::distance(vert.begin(), it);

  std::vector<apf::Vector3> rhs(vert.size(), apf::Vector3(0,0,0));
  std::vector<apf::Matrix3x3> m_eij(vert.size(), 
                              apf::Matrix3x3(0,0,0,0,0,0,0,0,0));
  apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

  std::vector<double> average_v(vert.size());
  if(drag_local) {
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = getAverageEntSize(m, vert.at(i), 2);
    }
  }
  else {
    double avg_temp = 0;
    get_average_tri();
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = average;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    m->getUp(vert.at(i), up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    m->getPoint(vert.at(i), 0, p_ctr);

    for(int j = 0; j < es_surf.size(); j++) {
      apf::ModelEntity* em = m->toModel(es_surf.at(j));
      int em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(es_surf.at(j), 0, d_v);

        int v1 = findIn(d_v, 3, vert.at(i));
        m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
        m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;

        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);

        if(!merg[vert.at(i)]) {
          rhs.at(i) = rhs.at(i) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i) = m_eij.at(i) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
        else {
          rhs.at(i_m) = rhs.at(i_m) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i_m) = m_eij.at(i_m) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
      }
    }
    if(!merg[vert.at(i)]) {
      m_eij.at(i) = m_eij.at(i) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
    else {
      m_eij.at(i_m) = m_eij.at(i_m) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    if(!merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i))*rhs.at(i)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
  for(int i = 0; i < vert.size(); i++) {
    if(merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i_m))*rhs.at(i_m)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_eqn_mason::vd_calc_force(apf::MeshEntity* vert) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < es_surf.size(); i++) {
    em = m->toModel(es_surf.at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);

      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_mason::vd_calc_force_tri(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  for(int i = 0; i < tris->size(); i++) {
    em = m->toModel(tris->at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {

      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, tris->at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_mason::vd_calc_force_edge(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  std::vector<apf::MeshEntity*> t_set(0);

  for(int i = 0; i < edges->size(); i++) {
    vd_set_up(m, edges->at(i), &t_set);
    for(int j = 0; j < t_set.size(); j++) {
      em = m->toModel(t_set.at(j));
      em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(t_set.at(j), 0, d_v);
        m->getDownward(t_set.at(j), 1, d_e);

        int v1 = findIn(d_v, 3, vert);
        int ii = lookup_tri_e[v1][0];
        int jj = lookup_tri_e[v1][1];
        int iv = lookup_tri_vd[v1][0];
        int jv = lookup_tri_vd[v1][1];
        if(d_e[ii] != edges->at(i)) {
          ii = lookup_tri_e[v1][1];
          jj = lookup_tri_e[v1][0];
          iv = lookup_tri_vd[v1][1];
          jv = lookup_tri_vd[v1][0];
        }
        m->getPoint(d_v[iv], 0, p_i);
        m->getPoint(d_v[jv], 0, p_j);

        rhs = rhs + vd_dir_in_pl(p_ctr, p_i, p_j)* 
                                  f_calc->gam2(m, t_set.at(j))*
                                  (p_i - p_ctr).getLength();
      }
    }
  }
  return rhs/2;
}

vd_eqn_mason::~vd_eqn_mason() {
}


//////////////////////////////////////
// vd_eqn_mason_NBC
//////////////////////////////////////
// Update the velocity field at:
// Every boundary vertex.


vd_eqn_mason_NBC::vd_eqn_mason_NBC() : vd_eqn_of_motion(),
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
}

vd_eqn_mason_NBC::vd_eqn_mason_NBC(apf::Mesh2* m_in, cell_base* c_base_in,
                            field_calc* f_calc_in, 
                            vd_entlist* e_list_in) :
                         vd_eqn_of_motion(m_in, c_base_in, f_calc_in, e_list_in), 
                           es_edge(0), es_surf(0), 
                           p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                           p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                           average(-1.), drag_rat(1000) {
  m = m_in;
  f_calc = f_calc_in;
  e_list = e_list_in;
  c_base = c_base_in;
}

vd_eqn_mason_NBC::vd_eqn_mason_NBC(const vd_eqn_mason_NBC& that) :
                               vd_eqn_of_motion(), 
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;
}

vd_eqn_mason_NBC& vd_eqn_mason_NBC::operator=(const vd_eqn_mason_NBC& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;

  es_edge.clear();
  es_surf.clear();
  apf::Vector3 zero(0,0,0);
  p_i = zero;
  p_j = zero;
  n_ij = zero;
  p_ctr = zero;
  temp = zero;
  rhs = zero;
  average = that.average;
  vel_field = that.vel_field;
  drag_rat = that.drag_rat; 
  return *this;
}

void vd_eqn_mason_NBC::get_average_tri() {
  average = f_calc->get_drag_glob();
  if(average < - std::numeric_limits<double>::min())
    average = getAverageEntSize(m, 2);
}

void vd_eqn_mason_NBC::calc_vel_curr(apf::MeshEntity* vert) {

  apf::Vector3 zero(0,0,0);
  apf::Vector3 area_proj(0,0,0);
  apf::Vector3 f_temp(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);
  //double norm_sum = 0;
  //double len_avg = 0;
  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* em = m->toModel(es_surf.at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      //area_proj = area_proj + n_ij;

      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      //norm_sum = norm_sum + norm;
      //len_avg = len_avg + (p_i.getLength() + p_j.getLength())/2;
      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      f_temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, es_surf.at(i));
      if(f_calc->chk_vert_special(m, vert)) {
        f_temp = f_calc->get_vec_special(m, vert, f_temp);
        n_ij = f_calc->get_vec_special(m, vert, n_ij);
      }
      rhs = rhs + f_temp;
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, es_surf.at(i));
    }
  }
  //len_avg = len_avg/es_surf.size();
  //area_proj = area_proj/6;

  m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;
  // 2017ActeMateMason:
  // rhs = F*2
  // m_eij = D*6
  // apf::invert(m_eij) : Mobility tensor/6.
  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
}

apf::Vector3 vd_eqn_mason_NBC::calc_vel_curr_tri(apf::MeshEntity* vert,
                           std::vector<apf::MeshEntity*>* tris) {

  apf::Vector3 zero(0,0,0);
  apf::Vector3 f_temp(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < tris->size(); i++) {
    apf::ModelEntity* em = m->toModel(tris->at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      f_temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, tris->at(i));
      if(f_calc->chk_vert_special(m, vert)) {
        f_temp = f_calc->get_vec_special(m, vert, f_temp);
        n_ij = f_calc->get_vec_special(m, vert, n_ij);
      }
      rhs = rhs + f_temp;
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, tris->at(i));
    }
  }
  m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;

  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
  return temp;
}

void vd_eqn_mason_NBC::calc_vel() {
  vel_field = m->findField("velocity_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      calc_vel_curr(vert);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  m->end(it);
  printf("Velocity field finished.\n");
}

void vd_eqn_mason_NBC::calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::Field* drag_field = vd_att_vm_field(m, "drag_field");
  apf::Field* force_field = vd_att_vv_field(m, "force_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);
  apf::Vector3 rhs(0,0,0);

  vd_set_up(m, verts, &es_edge);
  vd_set_up(m, &es_edge, &es_surf);

  std::vector<apf::Vector3> p(3, apf::Vector3(0,0,0));
  apf::Vector3 p_i(0,0,0);
  apf::Vector3 p_j(0,0,0);

  apf::Vector3 f_curr(0,0,0);
  apf::Vector3 f_temp(0,0,0);

  apf::Vector3 n_ij(0,0,0);
  apf::Vector3 n_ij_temp(0,0,0);

  apf::Matrix3x3 m_temp(0,0,0,0,0,0,0,0,0);
  // TODO don't need to calculate the edge directions twice for each edge...
  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* em = m->toModel(es_surf.at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      double gam_curr = f_calc->gam2(m, es_surf.at(i));
      double drag_curr = f_calc->d2(m, es_surf.at(i));

      m->getDownward(es_surf.at(i), 0, d_v);

      for(int j = 0; j < 3; j++) {
        m->getPoint(d_v[j], 0, p.at(j));
      }

      // vertex 1
      int v1 = 0;
      p_i = p.at(lookup_tri_vd[v1][0]) - p.at(v1);
      p_j = p.at(lookup_tri_vd[v1][1]) - p.at(v1);

      n_ij = vd_cross(p_i, p_j);
      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);
      n_ij_temp = n_ij;
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      f_curr = vd_cross(n_ij, p_i-p_j)*gam_curr;
      f_temp = f_curr;

      vert = d_v[v1];
      if(f_calc->chk_vert_special(m, vert)) {
        f_temp = f_calc->get_vec_special(m, vert, f_curr);
        n_ij_temp = f_calc->get_vec_special(m, vert, n_ij);
      }
      apf::getVector(force_field, vert, 0, rhs);
      rhs = rhs + f_temp;
      apf::setVector(force_field, vert, 0, rhs);

      apf::getMatrix(drag_field, vert, 0, m_temp);
      m_temp = m_temp + tensorProduct(n_ij_temp, n_ij_temp)*norm*drag_curr;
      apf::setMatrix(drag_field, vert, 0, m_temp);

      //m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*drag_curr;

      // vertex 2
      v1 = 1;
      p_i = p.at(lookup_tri_vd[v1][0]) - p.at(v1);
      p_j = p.at(lookup_tri_vd[v1][1]) - p.at(v1);

      f_curr = vd_cross(n_ij, p_i-p_j)*gam_curr;
      f_temp = f_curr;

      vert = d_v[v1];
      if(f_calc->chk_vert_special(m, vert)) {
        f_temp = f_calc->get_vec_special(m, vert, f_curr);
        n_ij_temp = f_calc->get_vec_special(m, vert, n_ij);
      }
      apf::getVector(force_field, vert, 0, rhs);
      rhs = rhs + f_temp;
      apf::setVector(force_field, vert, 0, rhs);

      apf::getMatrix(drag_field, vert, 0, m_temp);
      m_temp = m_temp + tensorProduct(n_ij_temp, n_ij_temp)*norm*drag_curr;
      apf::setMatrix(drag_field, vert, 0, m_temp);

      // vertex 3
      v1 = 2;
      vert = d_v[v1];

      p_i = p.at(lookup_tri_vd[v1][0]) - p.at(v1);
      p_j = p.at(lookup_tri_vd[v1][1]) - p.at(v1);
      f_curr = vd_cross(n_ij, p_i-p_j)*gam_curr;
      f_temp = f_curr;

      if(f_calc->chk_vert_special(m, vert)) {
        f_temp = f_calc->get_vec_special(m, vert, f_curr);
        n_ij_temp = f_calc->get_vec_special(m, vert, n_ij);
      }
      apf::getVector(force_field, vert, 0, rhs);
      rhs = rhs + f_temp;
      apf::setVector(force_field, vert, 0, rhs);

      apf::getMatrix(drag_field, vert, 0, m_temp);
      m_temp = m_temp + tensorProduct(n_ij_temp, n_ij_temp)*norm*drag_curr;
      apf::setMatrix(drag_field, vert, 0, m_temp);
    }
  }

  apf::Matrix3x3 m_sing(0,0,0,0,0,0,0,0,0);
  m_sing = M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;
  for(int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    apf::getVector(force_field, verts->at(i), 0, rhs);
    apf::getMatrix(drag_field, verts->at(i), 0, m_temp);
    m_temp = m_temp + m_sing;
    temp = apf::invert(m_temp)*rhs*3;
    assert(!std::isnan(temp.getLength()));
    apf::setVector(vel_field, verts->at(i), 0, temp);

    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  apf::destroyField(drag_field);
  apf::destroyField(force_field);

/*
  for(int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      vd_upd_vel_field(vert, drag_local);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
*/
  printf("Velocity field finished.\n");
}

// A single vertex.
void vd_eqn_mason_NBC::vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  if(drag_local) {
    average = getAverageEntSize(m, vert, 2);
  }
  else
    get_average_tri();
  calc_vel_curr(vert);
}

// A single vertex.
apf::Vector3 vd_eqn_mason_NBC::vd_upd_vel_field_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  if(drag_local) {
    average = getAverageEntSize(m, vert, 2);
  }
  else
    get_average_tri();
  return calc_vel_curr_tri(vert, tris);
}

apf::Vector3 vd_eqn_mason_NBC::vd_upd_vel_field_edge(apf::MeshEntity* vert, 
           std::vector<apf::MeshEntity*>* edges, bool drag_local) {
  std::cout << "upd_vel_field_edge not defined!" << std::endl;
  return apf::Vector3(0,0,0);
}

// A set of vertices while skipping labeled vertices. Merging vertices contribute
// to the target vertex and they all have the same resultant force and velocity.
void vd_eqn_mason_NBC::vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                       bool drag_local) {
  vel_field = m->findField("velocity_field");

  std::map<apf::MeshEntity*, bool> merg{};
  apf::MeshEntity* v_merg = NULL;
  for(int i = 0; i < vert.size(); i++) {
    if(v_map[vert.at(i)] != NULL) {
      merg[vert.at(i)] = true;
      if(v_merg == NULL)
        v_merg = v_map[vert.at(i)];
    }
  }
  assert(v_merg != NULL);

  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find(vert.begin(), vert.end(), v_merg);
  assert(it != vert.end());
  int i_m = std::distance(vert.begin(), it);

  apf::Vector3 f_temp(0,0,0);
  std::vector<apf::Vector3> rhs(vert.size(), apf::Vector3(0,0,0));
  std::vector<apf::Matrix3x3> m_eij(vert.size(), 
                              apf::Matrix3x3(0,0,0,0,0,0,0,0,0));
  apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

  std::vector<double> average_v(vert.size());
  if(drag_local) {
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = getAverageEntSize(m, vert.at(i), 2);
    }
  }
  else {
    double avg_temp = 0;
    get_average_tri();
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = average;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    m->getUp(vert.at(i), up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    m->getPoint(vert.at(i), 0, p_ctr);

    for(int j = 0; j < es_surf.size(); j++) {
      apf::ModelEntity* em = m->toModel(es_surf.at(j));
      int em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(es_surf.at(j), 0, d_v);

        int v1 = findIn(d_v, 3, vert.at(i));
        m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
        m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;

        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);

        if(!merg[vert.at(i)]) {
          f_temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, es_surf.at(j));
          if(f_calc->chk_vert_special(m, vert.at(i))) {
            f_temp = f_calc->get_vec_special(m, vert.at(i), f_temp);
            n_ij = f_calc->get_vec_special(m, vert.at(i), n_ij);
          }
          rhs.at(i) = rhs.at(i) + f_temp;
          m_eij.at(i) = m_eij.at(i) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
        else {
          f_temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, es_surf.at(j));
          if(f_calc->chk_vert_special(m, vert.at(i_m))) {
            f_temp = f_calc->get_vec_special(m, vert.at(i_m), f_temp);
            n_ij = f_calc->get_vec_special(m, vert.at(i_m), n_ij);
          }
          rhs.at(i_m) = rhs.at(i_m) + f_temp;
          m_eij.at(i_m) = m_eij.at(i_m) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
      }
    }
    if(!merg[vert.at(i)]) {
      m_eij.at(i) = m_eij.at(i) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
    else {
      m_eij.at(i_m) = m_eij.at(i_m) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    if(!merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i))*rhs.at(i)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
  for(int i = 0; i < vert.size(); i++) {
    if(merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i_m))*rhs.at(i_m)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_eqn_mason_NBC::vd_calc_force(apf::MeshEntity* vert) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < es_surf.size(); i++) {
    em = m->toModel(es_surf.at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);

      temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, es_surf.at(i));
      if(f_calc->chk_vert_special(m, vert)) {
        temp = f_calc->get_vec_special(m, vert, temp);
      }

      rhs = rhs + temp;
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_mason_NBC::vd_calc_force_tri(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  for(int i = 0; i < tris->size(); i++) {
    em = m->toModel(tris->at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {

      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);
      temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, tris->at(i));
      if(f_calc->chk_vert_special(m, vert)) {
        temp = f_calc->get_vec_special(m, vert, temp);
      }
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, tris->at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_mason_NBC::vd_calc_force_edge(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  std::vector<apf::MeshEntity*> t_set(0);

  for(int i = 0; i < edges->size(); i++) {
    vd_set_up(m, edges->at(i), &t_set);
    for(int j = 0; j < t_set.size(); j++) {
      em = m->toModel(t_set.at(j));
      em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(t_set.at(j), 0, d_v);
        m->getDownward(t_set.at(j), 1, d_e);

        int v1 = findIn(d_v, 3, vert);
        int ii = lookup_tri_e[v1][0];
        int jj = lookup_tri_e[v1][1];
        int iv = lookup_tri_vd[v1][0];
        int jv = lookup_tri_vd[v1][1];
        if(d_e[ii] != edges->at(i)) {
          ii = lookup_tri_e[v1][1];
          jj = lookup_tri_e[v1][0];
          iv = lookup_tri_vd[v1][1];
          jv = lookup_tri_vd[v1][0];
        }
        m->getPoint(d_v[iv], 0, p_i);
        m->getPoint(d_v[jv], 0, p_j);

        temp = vd_dir_in_pl(p_ctr, p_i, p_j)* 
                                  f_calc->gam2(m, t_set.at(j))*
                                  (p_i - p_ctr).getLength();
        if(f_calc->chk_vert_special(m, vert)) {
          temp = f_calc->get_vec_special(m, vert, temp);
        }

        rhs = rhs + temp;
      }
    }
  }
  return rhs/2;
}

vd_eqn_mason_NBC::~vd_eqn_mason_NBC() {
}


//////////////////////////////////////
// vd_eqn_mason_NBC_drag
//////////////////////////////////////
// Modified version of vd_eqn_mason_NBC that replaces identity matrix multiplied by the drag term and the average triangle area/100 by an anisotropic drag tensor which is a sum of the outer products of the 2tratum edge directions across the vertex. The drag tensor is normalized by its largest eigenvalue and multiplied by the drag term and the average triangle area.
// This heuristic effectively reduces the sliding.
// TODO Include a way of specifying the drag term calculation in each equation of motion implementation rather than in a separate class.
// Update the velocity field at:
// Every boundary vertex.


vd_eqn_mason_NBC_drag::vd_eqn_mason_NBC_drag() : vd_eqn_of_motion(),
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
}

vd_eqn_mason_NBC_drag::vd_eqn_mason_NBC_drag(apf::Mesh2* m_in, cell_base* c_base_in,
                            field_calc* f_calc_in, 
                            vd_entlist* e_list_in) :
                         vd_eqn_of_motion(m_in, c_base_in, f_calc_in, e_list_in), 
                           es_edge(0), es_surf(0), 
                           p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                           p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                           average(-1.), drag_rat(1000) {
  m = m_in;
  f_calc = f_calc_in;
  e_list = e_list_in;
  c_base = c_base_in;
}

vd_eqn_mason_NBC_drag::vd_eqn_mason_NBC_drag(const vd_eqn_mason_NBC_drag& that) :
                               vd_eqn_of_motion(), 
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;
}

vd_eqn_mason_NBC_drag& vd_eqn_mason_NBC_drag::operator=(const vd_eqn_mason_NBC_drag& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;

  es_edge.clear();
  es_surf.clear();
  apf::Vector3 zero(0,0,0);
  p_i = zero;
  p_j = zero;
  n_ij = zero;
  p_ctr = zero;
  temp = zero;
  rhs = zero;
  average = that.average;
  vel_field = that.vel_field;
  drag_rat = that.drag_rat; 
  return *this;
}

void vd_eqn_mason_NBC_drag::get_average_tri() {
  average = f_calc->get_drag_glob();
  if(average < - std::numeric_limits<double>::min())
    average = getAverageEntSize(m, 2);
}

void vd_eqn_mason_NBC_drag::calc_vel_curr(apf::MeshEntity* vert) {

  apf::Vector3 zero(0,0,0);
  apf::Vector3 area_proj(0,0,0);
  apf::Vector3 f_temp(0,0,0);
  apf::Vector3 lat_temp(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);
  apf::Matrix3x3 m_lat(0,0,0,0,0,0,0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);
  //double norm_sum = 0;
  //double len_avg = 0;
  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* em = m->toModel(es_surf.at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      //area_proj = area_proj + n_ij;

      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      //norm_sum = norm_sum + norm;
      //len_avg = len_avg + (p_i.getLength() + p_j.getLength())/2;
      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      f_temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, es_surf.at(i));
      if(f_calc->chk_vert_special(m, vert)) {
        f_temp = f_calc->get_vec_special(m, vert, f_temp);
        n_ij = f_calc->get_vec_special(m, vert, n_ij);
      }
      rhs = rhs + f_temp;
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, es_surf.at(i));
      lat_temp = p_i - p_j;
      m_lat = m_lat + tensorProduct(lat_temp, lat_temp);
    }
  }
  //len_avg = len_avg/es_surf.size();
  //area_proj = area_proj/6;
  apf::Vector3 eigenVectors[3];
  double eigenValues[3];
  eigen(m_lat, eigenVectors, eigenValues);
  int max_id = std::max_element(eigenValues, eigenValues + 3) - eigenValues;
  int min_id = std::min_element(eigenValues, eigenValues + 3) - eigenValues;
  if(eigenValues[min_id] < std::numeric_limits<double>::min())
    m_lat = m_lat + tensorProduct(eigenVectors[min_id], eigenVectors[min_id])/100.;

  if(eigenValues[max_id] > std::numeric_limits<double>::min())
    m_lat = m_lat/eigenValues[max_id];
  else
    m_lat = apf::Matrix3x3(1,0,0,0,1,0,0,0,1);

  m_eij = m_eij + m_lat*average*average*f_calc->get_d2_glob()/drag_rat;
  //m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;
  // 2017ActeMateMason:
  // rhs = F*2
  // m_eij = D*6
  // apf::invert(m_eij) : Mobility tensor/6.
  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
}

apf::Vector3 vd_eqn_mason_NBC_drag::calc_vel_curr_tri(apf::MeshEntity* vert,
                           std::vector<apf::MeshEntity*>* tris) {

  apf::Vector3 zero(0,0,0);
  apf::Vector3 f_temp(0,0,0);
  apf::Vector3 lat_temp(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);
  apf::Matrix3x3 m_lat(0,0,0,0,0,0,0,0,0);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < tris->size(); i++) {
    apf::ModelEntity* em = m->toModel(tris->at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      f_temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, tris->at(i));
      if(f_calc->chk_vert_special(m, vert)) {
        f_temp = f_calc->get_vec_special(m, vert, f_temp);
        n_ij = f_calc->get_vec_special(m, vert, n_ij);
      }
      rhs = rhs + f_temp;
      lat_temp = p_i - p_j;
      m_lat = m_lat + tensorProduct(lat_temp, lat_temp);
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, tris->at(i));
    }
  }

  apf::Vector3 eigenVectors[3];
  double eigenValues[3];
  eigen(m_lat, eigenVectors, eigenValues);
  int max_id = std::max_element(eigenValues, eigenValues + 3) - eigenValues;
  int min_id = std::min_element(eigenValues, eigenValues + 3) - eigenValues;
  if(eigenValues[min_id] < std::numeric_limits<double>::min())
    m_lat = m_lat + tensorProduct(eigenVectors[min_id], eigenVectors[min_id])/100.;

  if(eigenValues[max_id] > std::numeric_limits<double>::min())
    m_lat = m_lat/eigenValues[max_id];
  else
    m_lat = apf::Matrix3x3(1,0,0,0,1,0,0,0,1);

  m_eij = m_eij + m_lat*average*average*f_calc->get_d2_glob()/drag_rat;
  //m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;

  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
  return temp;
}

void vd_eqn_mason_NBC_drag::calc_vel() {
  vel_field = m->findField("velocity_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      calc_vel_curr(vert);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  m->end(it);
  printf("Velocity field finished.\n");
}

void vd_eqn_mason_NBC_drag::calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::Field* drag_field = vd_att_vm_field(m, "drag_field");
  apf::Field* drag_lat_field = vd_att_vm_field(m, "drag_lat_field");
  apf::Field* force_field = vd_att_vv_field(m, "force_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);
  apf::Vector3 rhs(0,0,0);

  vd_set_up(m, verts, &es_edge);
  vd_set_up(m, &es_edge, &es_surf);

  std::vector<apf::Vector3> p(3, apf::Vector3(0,0,0));
  apf::Vector3 p_i(0,0,0);
  apf::Vector3 p_j(0,0,0);

  apf::Vector3 f_curr(0,0,0);
  apf::Vector3 f_temp(0,0,0);
  apf::Vector3 lat_temp(0,0,0);

  apf::Vector3 n_ij(0,0,0);
  apf::Vector3 n_ij_temp(0,0,0);

  apf::Matrix3x3 m_temp(0,0,0,0,0,0,0,0,0);
  apf::Matrix3x3 m_lat(0,0,0,0,0,0,0,0,0);

  // TODO don't need to calculate the edge directions twice for each edge...
  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* em = m->toModel(es_surf.at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      double gam_curr = f_calc->gam2(m, es_surf.at(i));
      double drag_curr = f_calc->d2(m, es_surf.at(i));

      m->getDownward(es_surf.at(i), 0, d_v);

      for(int j = 0; j < 3; j++) {
        m->getPoint(d_v[j], 0, p.at(j));
      }

      // vertex 1
      int v1 = 0;
      p_i = p.at(lookup_tri_vd[v1][0]) - p.at(v1);
      p_j = p.at(lookup_tri_vd[v1][1]) - p.at(v1);

      n_ij = vd_cross(p_i, p_j);
      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);
      n_ij_temp = n_ij;
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      f_curr = vd_cross(n_ij, p_i-p_j)*gam_curr;
      f_temp = f_curr;

      vert = d_v[v1];
      if(f_calc->chk_vert_special(m, vert)) {
        f_temp = f_calc->get_vec_special(m, vert, f_curr);
        n_ij_temp = f_calc->get_vec_special(m, vert, n_ij);
      }
      apf::getVector(force_field, vert, 0, rhs);
      rhs = rhs + f_temp;
      apf::setVector(force_field, vert, 0, rhs);

      apf::getMatrix(drag_field, vert, 0, m_temp);
      m_temp = m_temp + tensorProduct(n_ij_temp, n_ij_temp)*norm*drag_curr;
      apf::setMatrix(drag_field, vert, 0, m_temp);
      apf::getMatrix(drag_field, vert, 0, m_lat);
      lat_temp = p_i - p_j;
      m_lat = m_lat + tensorProduct(lat_temp, lat_temp);
      apf::setMatrix(drag_lat_field, vert, 0, m_lat);

      //m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*drag_curr;

      // vertex 2
      v1 = 1;
      p_i = p.at(lookup_tri_vd[v1][0]) - p.at(v1);
      p_j = p.at(lookup_tri_vd[v1][1]) - p.at(v1);

      f_curr = vd_cross(n_ij, p_i-p_j)*gam_curr;
      f_temp = f_curr;

      vert = d_v[v1];
      if(f_calc->chk_vert_special(m, vert)) {
        f_temp = f_calc->get_vec_special(m, vert, f_curr);
        n_ij_temp = f_calc->get_vec_special(m, vert, n_ij);
      }
      apf::getVector(force_field, vert, 0, rhs);
      rhs = rhs + f_temp;
      apf::setVector(force_field, vert, 0, rhs);

      apf::getMatrix(drag_field, vert, 0, m_temp);
      m_temp = m_temp + tensorProduct(n_ij_temp, n_ij_temp)*norm*drag_curr;
      apf::setMatrix(drag_field, vert, 0, m_temp);
      apf::getMatrix(drag_field, vert, 0, m_lat);
      lat_temp = p_i - p_j;
      m_lat = m_lat + tensorProduct(lat_temp, lat_temp);
      apf::setMatrix(drag_lat_field, vert, 0, m_lat);

      // vertex 3
      v1 = 2;
      vert = d_v[v1];

      p_i = p.at(lookup_tri_vd[v1][0]) - p.at(v1);
      p_j = p.at(lookup_tri_vd[v1][1]) - p.at(v1);
      f_curr = vd_cross(n_ij, p_i-p_j)*gam_curr;
      f_temp = f_curr;

      if(f_calc->chk_vert_special(m, vert)) {
        f_temp = f_calc->get_vec_special(m, vert, f_curr);
        n_ij_temp = f_calc->get_vec_special(m, vert, n_ij);
      }
      apf::getVector(force_field, vert, 0, rhs);
      rhs = rhs + f_temp;
      apf::setVector(force_field, vert, 0, rhs);

      apf::getMatrix(drag_field, vert, 0, m_temp);
      m_temp = m_temp + tensorProduct(n_ij_temp, n_ij_temp)*norm*drag_curr;
      apf::setMatrix(drag_field, vert, 0, m_temp);
      apf::getMatrix(drag_field, vert, 0, m_lat);
      lat_temp = p_i - p_j;
      m_lat = m_lat + tensorProduct(lat_temp, lat_temp);
      apf::setMatrix(drag_lat_field, vert, 0, m_lat);
    }
  }

  apf::Matrix3x3 m_sing(0,0,0,0,0,0,0,0,0);
  m_sing = M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;
  for(int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    apf::getVector(force_field, verts->at(i), 0, rhs);
    apf::getMatrix(drag_field, verts->at(i), 0, m_temp);
    apf::getMatrix(drag_lat_field, verts->at(i), 0, m_lat);
    m_temp = m_temp + m_sing;

    apf::Vector3 eigenVectors[3];
    double eigenValues[3];
    eigen(m_lat, eigenVectors, eigenValues);
    int max_id = std::max_element(eigenValues, eigenValues + 3) - eigenValues;
    int min_id = std::min_element(eigenValues, eigenValues + 3) - eigenValues;
    if(eigenValues[min_id] < std::numeric_limits<double>::min())
      m_lat = m_lat + tensorProduct(eigenVectors[min_id], eigenVectors[min_id])/100.;

    if(eigenValues[max_id] > std::numeric_limits<double>::min())
      m_lat = m_lat/eigenValues[max_id];
    else
      m_lat = apf::Matrix3x3(1,0,0,0,1,0,0,0,1);

    m_temp = m_temp + m_lat*average*average*f_calc->get_d2_glob()/drag_rat;

    temp = apf::invert(m_temp)*rhs*3;
    assert(!std::isnan(temp.getLength()));
    apf::setVector(vel_field, verts->at(i), 0, temp);

    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  apf::destroyField(drag_field);
  apf::destroyField(drag_lat_field);
  apf::destroyField(force_field);

/*
  for(int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      vd_upd_vel_field(vert, drag_local);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
*/
  printf("Velocity field finished.\n");
}

// A single vertex.
void vd_eqn_mason_NBC_drag::vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  if(drag_local) {
    average = getAverageEntSize(m, vert, 2);
  }
  else
    get_average_tri();
  calc_vel_curr(vert);
}

// A single vertex.
apf::Vector3 vd_eqn_mason_NBC_drag::vd_upd_vel_field_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  if(drag_local) {
    average = getAverageEntSize(m, vert, 2);
  }
  else
    get_average_tri();
  return calc_vel_curr_tri(vert, tris);
}

apf::Vector3 vd_eqn_mason_NBC_drag::vd_upd_vel_field_edge(apf::MeshEntity* vert, 
           std::vector<apf::MeshEntity*>* edges, bool drag_local) {
  std::cout << "upd_vel_field_edge not defined!" << std::endl;
  return apf::Vector3(0,0,0);
}

// A set of vertices while skipping labeled vertices. Merging vertices contribute
// to the target vertex and they all have the same resultant force and velocity.
void vd_eqn_mason_NBC_drag::vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                       bool drag_local) {
  vel_field = m->findField("velocity_field");

  std::map<apf::MeshEntity*, bool> merg{};
  apf::MeshEntity* v_merg = NULL;
  for(int i = 0; i < vert.size(); i++) {
    if(v_map[vert.at(i)] != NULL) {
      merg[vert.at(i)] = true;
      if(v_merg == NULL)
        v_merg = v_map[vert.at(i)];
    }
  }
  assert(v_merg != NULL);

  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find(vert.begin(), vert.end(), v_merg);
  assert(it != vert.end());
  int i_m = std::distance(vert.begin(), it);

  apf::Vector3 f_temp(0,0,0);
  std::vector<apf::Vector3> rhs(vert.size(), apf::Vector3(0,0,0));
  std::vector<apf::Matrix3x3> m_eij(vert.size(), 
                              apf::Matrix3x3(0,0,0,0,0,0,0,0,0));
  apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

  std::vector<double> average_v(vert.size());
  if(drag_local) {
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = getAverageEntSize(m, vert.at(i), 2);
    }
  }
  else {
    double avg_temp = 0;
    get_average_tri();
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = average;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    m->getUp(vert.at(i), up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    m->getPoint(vert.at(i), 0, p_ctr);

    for(int j = 0; j < es_surf.size(); j++) {
      apf::ModelEntity* em = m->toModel(es_surf.at(j));
      int em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(es_surf.at(j), 0, d_v);

        int v1 = findIn(d_v, 3, vert.at(i));
        m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
        m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;

        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);

        if(!merg[vert.at(i)]) {
          f_temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, es_surf.at(j));
          if(f_calc->chk_vert_special(m, vert.at(i))) {
            f_temp = f_calc->get_vec_special(m, vert.at(i), f_temp);
            n_ij = f_calc->get_vec_special(m, vert.at(i), n_ij);
          }
          rhs.at(i) = rhs.at(i) + f_temp;
          m_eij.at(i) = m_eij.at(i) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
        else {
          f_temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, es_surf.at(j));
          if(f_calc->chk_vert_special(m, vert.at(i_m))) {
            f_temp = f_calc->get_vec_special(m, vert.at(i_m), f_temp);
            n_ij = f_calc->get_vec_special(m, vert.at(i_m), n_ij);
          }
          rhs.at(i_m) = rhs.at(i_m) + f_temp;
          m_eij.at(i_m) = m_eij.at(i_m) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
      }
    }
    if(!merg[vert.at(i)]) {
      m_eij.at(i) = m_eij.at(i) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
    else {
      m_eij.at(i_m) = m_eij.at(i_m) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    if(!merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i))*rhs.at(i)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
  for(int i = 0; i < vert.size(); i++) {
    if(merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i_m))*rhs.at(i_m)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_eqn_mason_NBC_drag::vd_calc_force(apf::MeshEntity* vert) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < es_surf.size(); i++) {
    em = m->toModel(es_surf.at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);

      temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, es_surf.at(i));
      if(f_calc->chk_vert_special(m, vert)) {
        temp = f_calc->get_vec_special(m, vert, temp);
      }

      rhs = rhs + temp;
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_mason_NBC_drag::vd_calc_force_tri(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  for(int i = 0; i < tris->size(); i++) {
    em = m->toModel(tris->at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {

      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);
      temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, tris->at(i));
      if(f_calc->chk_vert_special(m, vert)) {
        temp = f_calc->get_vec_special(m, vert, temp);
      }
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, tris->at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_mason_NBC_drag::vd_calc_force_edge(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  std::vector<apf::MeshEntity*> t_set(0);

  for(int i = 0; i < edges->size(); i++) {
    vd_set_up(m, edges->at(i), &t_set);
    for(int j = 0; j < t_set.size(); j++) {
      em = m->toModel(t_set.at(j));
      em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(t_set.at(j), 0, d_v);
        m->getDownward(t_set.at(j), 1, d_e);

        int v1 = findIn(d_v, 3, vert);
        int ii = lookup_tri_e[v1][0];
        int jj = lookup_tri_e[v1][1];
        int iv = lookup_tri_vd[v1][0];
        int jv = lookup_tri_vd[v1][1];
        if(d_e[ii] != edges->at(i)) {
          ii = lookup_tri_e[v1][1];
          jj = lookup_tri_e[v1][0];
          iv = lookup_tri_vd[v1][1];
          jv = lookup_tri_vd[v1][0];
        }
        m->getPoint(d_v[iv], 0, p_i);
        m->getPoint(d_v[jv], 0, p_j);

        temp = vd_dir_in_pl(p_ctr, p_i, p_j)* 
                                  f_calc->gam2(m, t_set.at(j))*
                                  (p_i - p_ctr).getLength();
        if(f_calc->chk_vert_special(m, vert)) {
          temp = f_calc->get_vec_special(m, vert, temp);
        }

        rhs = rhs + temp;
      }
    }
  }
  return rhs/2;
}

vd_eqn_mason_NBC_drag::~vd_eqn_mason_NBC_drag() {
}

//////////////////////////////////////
// vd_eqn_kuprat
//////////////////////////////////////
// Update the velocity field at:
// Every boundary vertex.
// TODO calc_vel is global, the rest of the functions are local and the same as 
// vd_eqn_mason.

vd_eqn_kuprat::vd_eqn_kuprat() : vd_eqn_of_motion(),
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
}

vd_eqn_kuprat::vd_eqn_kuprat(apf::Mesh2* m_in, cell_base* c_base_in,
                            field_calc* f_calc_in, 
                            vd_entlist* e_list_in) :
                         vd_eqn_of_motion(m_in, c_base_in, f_calc_in, e_list_in), 
                           es_edge(0), es_surf(0), 
                           p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                           p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                           average(-1.), vel_field(NULL), drag_rat(1000) {
  m = m_in;
  f_calc = f_calc_in;
  e_list = e_list_in;
  c_base = c_base_in;
}

vd_eqn_kuprat::vd_eqn_kuprat(const vd_eqn_kuprat& that) :
                               vd_eqn_of_motion(), 
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;
}

vd_eqn_kuprat& vd_eqn_kuprat::operator=(const vd_eqn_kuprat& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;

  es_edge.clear();
  es_surf.clear();
  apf::Vector3 zero(0,0,0);
  p_i = zero;
  p_j = zero;
  n_ij = zero;
  p_ctr = zero;
  temp = zero;
  rhs = zero;
  average = that.average;
  vel_field = that.vel_field;
  drag_rat = that.drag_rat; 
  return *this;
}

void vd_eqn_kuprat::get_average_tri() {
  average = f_calc->get_drag_glob();
  if(average < - std::numeric_limits<double>::min())
    average = getAverageEntSize(m, 2);
}

void vd_eqn_kuprat::calc_vel_curr(apf::MeshEntity* vert) {
  vel_field = m->findField("velocity_field");

  apf::Vector3 zero(0,0,0);
  apf::Vector3 area_proj(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);
  double norm_sum = 0;
  double len_avg = 0;
  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* em = m->toModel(es_surf.at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      area_proj = area_proj + n_ij;

      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      norm_sum = norm_sum + norm;
      len_avg = len_avg + (p_i.getLength() + p_j.getLength())/2;
      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                            f_calc->gam2(m, es_surf.at(i));
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, es_surf.at(i));
    }
  }
  len_avg = len_avg/es_surf.size();
  area_proj = area_proj/6;

  m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;
  // 2017ActeMateMason:
  // rhs = F*2
  // m_eij = D*6
  // apf::invert(m_eij) : Mobility tensor/6.
  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
}

apf::Vector3 vd_eqn_kuprat::calc_vel_curr_tri(apf::MeshEntity* vert,
                           std::vector<apf::MeshEntity*>* tris) {

  apf::Vector3 zero(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < tris->size(); i++) {
    apf::ModelEntity* em = m->toModel(tris->at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                            f_calc->gam2(m, tris->at(i));
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, tris->at(i));
    }
  }
  m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;

  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
  return temp;
}

void vd_eqn_kuprat::calc_vel() {
  vel_field = m->findField("velocity_field");
  solv_kuprat sl(m, f_calc);

  apf::MeshEntity* vert;
  apf::MeshIterator* it = m->begin(0);
  int count = 0;
  while ((vert = m->iterate(it))) {
    count = count + 1;
  }
  m->end(it);

  // The tolerance is defined as |Ax - b| <= tol * |b|
  // So it should actually decrease with the number of nodes, rather than 
  // increasing!
  //double residual_target = 1.0e-3*std::sqrt(count);
  double residual_target = 1.0e-6/std::sqrt(count);
  double residual = sl.solve(residual_target);

  apf::Vector3 zero(0,0,0);

  it = m->begin(0);
  while ((vert = m->iterate(it))) {
    apf::ModelEntity* em = m->toModel(vert);
    int em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  m->end(it);

  printf("Velocity field finished target %.4E tol %.4E.\n", residual_target, residual);
}

void vd_eqn_kuprat::calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  solv_kuprat sl(m, f_calc, verts);
  // The tolerance is defined as |Ax - b| <= tol * |b|
  // So it should actually decrease with the number of nodes, rather than 
  // increasing!
  //double residual_target = 1.0e-3*std::sqrt(count);
  double residual_target = 1.0e-6/std::sqrt(verts->size());
  double residual = sl.solve(residual_target);

  for(int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    apf::ModelEntity* em = m->toModel(vert);
    int em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  printf("Velocity field finished target %.4E tol %.4E.\n", residual_target, residual);
}

// A single vertex.
void vd_eqn_kuprat::vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  if(drag_local) {
    average = getAverageEntSize(m, vert, 2);
  }
  else
    get_average_tri();
  calc_vel_curr(vert);
}

// A single vertex.
apf::Vector3 vd_eqn_kuprat::vd_upd_vel_field_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  if(drag_local) {
    average = getAverageEntSize(m, vert, 2);
  }
  else
    get_average_tri();
  return calc_vel_curr_tri(vert, tris);
}

apf::Vector3 vd_eqn_kuprat::vd_upd_vel_field_edge(apf::MeshEntity* vert, 
           std::vector<apf::MeshEntity*>* edges, bool drag_local) {
  std::cout << "upd_vel_field_edge not defined!" << std::endl;
  return apf::Vector3(0,0,0);
}


// A set of vertices while skipping labeled vertices. Merging vertices contribute
// to the target vertex and they all have the same resultant force and velocity.
void vd_eqn_kuprat::vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                       bool drag_local) {
  vel_field = m->findField("velocity_field");

  std::map<apf::MeshEntity*, bool> merg{};
  apf::MeshEntity* v_merg = NULL;
  for(int i = 0; i < vert.size(); i++) {
    if(v_map[vert.at(i)] != NULL) {
      merg[vert.at(i)] = true;
      if(v_merg == NULL)
        v_merg = v_map[vert.at(i)];
    }
  }
  assert(v_merg != NULL);

  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find(vert.begin(), vert.end(), v_merg);
  assert(it != vert.end());
  int i_m = std::distance(vert.begin(), it);

  std::vector<apf::Vector3> rhs(vert.size(), apf::Vector3(0,0,0));
  std::vector<apf::Matrix3x3> m_eij(vert.size(), 
                              apf::Matrix3x3(0,0,0,0,0,0,0,0,0));
  apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

  std::vector<double> average_v(vert.size());
  if(drag_local) {
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = getAverageEntSize(m, vert.at(i), 2);
    }
  }
  else {
    double avg_temp = 0;
    get_average_tri();
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = average;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    m->getUp(vert.at(i), up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    m->getPoint(vert.at(i), 0, p_ctr);

    for(int j = 0; j < es_surf.size(); j++) {
      apf::ModelEntity* em = m->toModel(es_surf.at(j));
      int em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(es_surf.at(j), 0, d_v);

        int v1 = findIn(d_v, 3, vert.at(i));
        m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
        m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;

        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);

        if(!merg[vert.at(i)]) {
          rhs.at(i) = rhs.at(i) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i) = m_eij.at(i) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
        else {
          rhs.at(i_m) = rhs.at(i_m) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i_m) = m_eij.at(i_m) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
      }
    }
    if(!merg[vert.at(i)]) {
      m_eij.at(i) = m_eij.at(i) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
    else {
      m_eij.at(i_m) = m_eij.at(i_m) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    if(!merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i))*rhs.at(i)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
  for(int i = 0; i < vert.size(); i++) {
    if(merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i_m))*rhs.at(i_m)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_eqn_kuprat::vd_calc_force(apf::MeshEntity* vert) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < es_surf.size(); i++) {
    em = m->toModel(es_surf.at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);

      // This is the actual force, rather than the quantity corresponding to the  
      // reduced mobility.
      rhs = rhs + vd_cross(n_ij, p_i-p_j)*
                                f_calc->gam2(m, es_surf.at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_kuprat::vd_calc_force_tri(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  for(int i = 0; i < tris->size(); i++) {
    em = m->toModel(tris->at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {

      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, tris->at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_kuprat::vd_calc_force_edge(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  std::vector<apf::MeshEntity*> t_set(0);

  for(int i = 0; i < edges->size(); i++) {
    vd_set_up(m, edges->at(i), &t_set);
    for(int j = 0; j < t_set.size(); j++) {
      em = m->toModel(t_set.at(j));
      em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(t_set.at(j), 0, d_v);
        m->getDownward(t_set.at(j), 1, d_e);

        int v1 = findIn(d_v, 3, vert);
        int ii = lookup_tri_e[v1][0];
        int jj = lookup_tri_e[v1][1];
        int iv = lookup_tri_vd[v1][0];
        int jv = lookup_tri_vd[v1][1];
        if(d_e[ii] != edges->at(i)) {
          ii = lookup_tri_e[v1][1];
          jj = lookup_tri_e[v1][0];
          iv = lookup_tri_vd[v1][1];
          jv = lookup_tri_vd[v1][0];
        }
        m->getPoint(d_v[iv], 0, p_i);
        m->getPoint(d_v[jv], 0, p_j);

        rhs = rhs + vd_dir_in_pl(p_ctr, p_i, p_j)* 
                                  f_calc->gam2(m, t_set.at(j))*
                                  (p_i - p_ctr).getLength();
      }
    }
  }
  return rhs/2;
}

vd_eqn_kuprat::~vd_eqn_kuprat() {
}


//////////////////////////////////////
// vd_eqn_kuprat_NBC
//////////////////////////////////////
// Update the velocity field at:
// Every boundary vertex.
// TODO calc_vel is global, the rest of the functions are local and the same as 
// vd_eqn_mason.

vd_eqn_kuprat_NBC::vd_eqn_kuprat_NBC() : vd_eqn_of_motion(),
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
}

vd_eqn_kuprat_NBC::vd_eqn_kuprat_NBC(apf::Mesh2* m_in, cell_base* c_base_in,
                            field_calc* f_calc_in, 
                            vd_entlist* e_list_in) :
                         vd_eqn_of_motion(m_in, c_base_in, f_calc_in, e_list_in), 
                           es_edge(0), es_surf(0), 
                           p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                           p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                           average(-1.), drag_rat(1000) {
  m = m_in;
  f_calc = f_calc_in;
  e_list = e_list_in;
  c_base = c_base_in;
}

vd_eqn_kuprat_NBC::vd_eqn_kuprat_NBC(const vd_eqn_kuprat_NBC& that) :
                               vd_eqn_of_motion(), 
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;
}

vd_eqn_kuprat_NBC& vd_eqn_kuprat_NBC::operator=(const vd_eqn_kuprat_NBC& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;

  es_edge.clear();
  es_surf.clear();
  apf::Vector3 zero(0,0,0);
  p_i = zero;
  p_j = zero;
  n_ij = zero;
  p_ctr = zero;
  temp = zero;
  rhs = zero;
  average = that.average;
  vel_field = that.vel_field;
  drag_rat = that.drag_rat; 
  return *this;
}

void vd_eqn_kuprat_NBC::get_average_tri() {
  average = f_calc->get_drag_glob();
  if(average < - std::numeric_limits<double>::min())
    average = getAverageEntSize(m, 2);
}

void vd_eqn_kuprat_NBC::calc_vel_curr(apf::MeshEntity* vert) {

  vel_field = m->findField("velocity_field");
  apf::Vector3 zero(0,0,0);
  apf::Vector3 area_proj(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);
  double norm_sum = 0;
  double len_avg = 0;
  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* em = m->toModel(es_surf.at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      area_proj = area_proj + n_ij;

      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      norm_sum = norm_sum + norm;
      len_avg = len_avg + (p_i.getLength() + p_j.getLength())/2;
      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                            f_calc->gam2(m, es_surf.at(i));
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, es_surf.at(i));
    }
  }
  len_avg = len_avg/es_surf.size();
  area_proj = area_proj/6;

  m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;
  // 2017ActeMateMason:
  // rhs = F*2
  // m_eij = D*6
  // apf::invert(m_eij) : Mobility tensor/6.
  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
}

apf::Vector3 vd_eqn_kuprat_NBC::calc_vel_curr_tri(apf::MeshEntity* vert,
                           std::vector<apf::MeshEntity*>* tris) {

  vel_field = m->findField("velocity_field");
  apf::Vector3 zero(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < tris->size(); i++) {
    apf::ModelEntity* em = m->toModel(tris->at(i));
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                            f_calc->gam2(m, tris->at(i));
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, tris->at(i));
    }
  }
  m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;

  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
  return temp;
}

void vd_eqn_kuprat_NBC::calc_vel() {
  vel_field = m->findField("velocity_field");
  solv_kuprat_NBC sl(m, f_calc);

  apf::MeshEntity* vert;
  apf::MeshIterator* it = m->begin(0);
  int count = 0;
  while ((vert = m->iterate(it))) {
    count = count + 1;
  }
  m->end(it);

  // The tolerance is defined as |Ax - b| <= tol * |b|
  // So it should actually decrease with the number of nodes, rather than 
  // increasing!
  //double residual_target = 1.0e-3*std::sqrt(count);
  double residual_target = 1.0e-6/std::sqrt(count);
  double residual = sl.solve(residual_target);

  apf::Vector3 zero(0,0,0);

  it = m->begin(0);
  while ((vert = m->iterate(it))) {
    apf::ModelEntity* em = m->toModel(vert);
    int em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  m->end(it);

  printf("Velocity field finished target %.4E tol %.4E.\n", residual_target, residual);
}

void vd_eqn_kuprat_NBC::calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  solv_kuprat_NBC sl(m, f_calc, verts);
  // The tolerance is defined as |Ax - b| <= tol * |b|
  // So it should actually decrease with the number of nodes, rather than 
  // increasing!
  //double residual_target = 1.0e-3*std::sqrt(count);
  double residual_target = 1.0e-6/std::sqrt(verts->size());
  double residual = sl.solve(residual_target);

  for(int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    apf::ModelEntity* em = m->toModel(vert);
    int em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  printf("Velocity field finished target %.4E tol %.4E.\n", residual_target, residual);
}

// A single vertex.
void vd_eqn_kuprat_NBC::vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  if(drag_local) {
    average = getAverageEntSize(m, vert, 2);
  }
  else
    get_average_tri();
  calc_vel_curr(vert);
}

// A single vertex.
apf::Vector3 vd_eqn_kuprat_NBC::vd_upd_vel_field_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  if(drag_local) {
    average = getAverageEntSize(m, vert, 2);
  }
  else
    get_average_tri();
  return calc_vel_curr_tri(vert, tris);
}

apf::Vector3 vd_eqn_kuprat_NBC::vd_upd_vel_field_edge(apf::MeshEntity* vert, 
           std::vector<apf::MeshEntity*>* edges, bool drag_local) {
  std::cout << "upd_vel_field_edge not defined!" << std::endl;
  return apf::Vector3(0,0,0);
}


// A set of vertices while skipping labeled vertices. Merging vertices contribute
// to the target vertex and they all have the same resultant force and velocity.
void vd_eqn_kuprat_NBC::vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                       bool drag_local) {
  vel_field = m->findField("velocity_field");

  std::map<apf::MeshEntity*, bool> merg{};
  apf::MeshEntity* v_merg = NULL;
  for(int i = 0; i < vert.size(); i++) {
    if(v_map[vert.at(i)] != NULL) {
      merg[vert.at(i)] = true;
      if(v_merg == NULL)
        v_merg = v_map[vert.at(i)];
    }
  }
  assert(v_merg != NULL);

  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find(vert.begin(), vert.end(), v_merg);
  assert(it != vert.end());
  int i_m = std::distance(vert.begin(), it);

  std::vector<apf::Vector3> rhs(vert.size(), apf::Vector3(0,0,0));
  std::vector<apf::Matrix3x3> m_eij(vert.size(), 
                              apf::Matrix3x3(0,0,0,0,0,0,0,0,0));
  apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

  std::vector<double> average_v(vert.size());
  if(drag_local) {
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = getAverageEntSize(m, vert.at(i), 2);
    }
  }
  else {
    double avg_temp = 0;
    get_average_tri();
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = average;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    m->getUp(vert.at(i), up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    m->getPoint(vert.at(i), 0, p_ctr);

    for(int j = 0; j < es_surf.size(); j++) {
      apf::ModelEntity* em = m->toModel(es_surf.at(j));
      int em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(es_surf.at(j), 0, d_v);

        int v1 = findIn(d_v, 3, vert.at(i));
        m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
        m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;

        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);

        if(!merg[vert.at(i)]) {
          rhs.at(i) = rhs.at(i) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i) = m_eij.at(i) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
        else {
          rhs.at(i_m) = rhs.at(i_m) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i_m) = m_eij.at(i_m) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
      }
    }
    if(!merg[vert.at(i)]) {
      m_eij.at(i) = m_eij.at(i) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
    else {
      m_eij.at(i_m) = m_eij.at(i_m) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    if(!merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i))*rhs.at(i)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
  for(int i = 0; i < vert.size(); i++) {
    if(merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i_m))*rhs.at(i_m)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_eqn_kuprat_NBC::vd_calc_force(apf::MeshEntity* vert) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < es_surf.size(); i++) {
    em = m->toModel(es_surf.at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);

      // This is the actual force, rather than the quantity corresponding to the  
      // reduced mobility.
      rhs = rhs + vd_cross(n_ij, p_i-p_j)*
                                f_calc->gam2(m, es_surf.at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_kuprat_NBC::vd_calc_force_tri(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  for(int i = 0; i < tris->size(); i++) {
    em = m->toModel(tris->at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {

      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, tris->at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_kuprat_NBC::vd_calc_force_edge(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  std::vector<apf::MeshEntity*> t_set(0);

  for(int i = 0; i < edges->size(); i++) {
    vd_set_up(m, edges->at(i), &t_set);
    for(int j = 0; j < t_set.size(); j++) {
      em = m->toModel(t_set.at(j));
      em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(t_set.at(j), 0, d_v);
        m->getDownward(t_set.at(j), 1, d_e);

        int v1 = findIn(d_v, 3, vert);
        int ii = lookup_tri_e[v1][0];
        int jj = lookup_tri_e[v1][1];
        int iv = lookup_tri_vd[v1][0];
        int jv = lookup_tri_vd[v1][1];
        if(d_e[ii] != edges->at(i)) {
          ii = lookup_tri_e[v1][1];
          jj = lookup_tri_e[v1][0];
          iv = lookup_tri_vd[v1][1];
          jv = lookup_tri_vd[v1][0];
        }
        m->getPoint(d_v[iv], 0, p_i);
        m->getPoint(d_v[jv], 0, p_j);

        rhs = rhs + vd_dir_in_pl(p_ctr, p_i, p_j)* 
                                  f_calc->gam2(m, t_set.at(j))*
                                  (p_i - p_ctr).getLength();
      }
    }
  }
  return rhs/2;
}

vd_eqn_kuprat_NBC::~vd_eqn_kuprat_NBC() {
}

//////////////////////////////////////
// vd_eqn_lazar
//////////////////////////////////////
// Update the velocity field at:
// Every boundary vertex.
vd_eqn_lazar::vd_eqn_lazar() : vd_eqn_of_motion(),
                               es_edge(0), es_surf(0), 
                               b_verts(0),
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               tri_tet_in{}, edge_len{}, tri_norm{}, 
                               edge_pos{}, tri_pos{}, a_equi_map{},
                               average(-1.), vel_field(NULL), drag_rat(1000) {
}

vd_eqn_lazar::vd_eqn_lazar(apf::Mesh2* m_in, cell_base* c_base_in, 
                            field_calc* f_calc_in, 
                            vd_entlist* e_list_in) :
                         vd_eqn_of_motion(m_in, c_base_in, f_calc_in, e_list_in), 
                           es_edge(0), es_surf(0), 
                           b_verts(0),
                           p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                           p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                           tri_tet_in{}, edge_len{}, tri_norm{}, 
                           edge_pos{}, tri_pos{}, a_equi_map{},
                           average(-1.), drag_rat(1000) {
  m = m_in;
  f_calc = f_calc_in;
  e_list = e_list_in;
  c_base = c_base_in;
}

vd_eqn_lazar::vd_eqn_lazar(const vd_eqn_lazar& that) :
                               vd_eqn_of_motion(), 
                               es_edge(0), es_surf(0), 
                               b_verts(0),
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               tri_tet_in{}, edge_len{}, tri_norm{}, 
                               edge_pos{}, tri_pos{}, a_equi_map{},
                               average(-1.), vel_field(NULL), drag_rat(1000) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;
}

vd_eqn_lazar& vd_eqn_lazar::operator=(const vd_eqn_lazar& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;

  es_edge.clear();
  es_surf.clear();

  tri_tet_in.clear();
  edge_len.clear();
  tri_norm.clear(); 
  edge_pos.clear();
  tri_pos.clear();
  a_equi_map.clear();

  b_verts = that.b_verts;

  apf::Vector3 zero(0,0,0);
  p_i = zero;
  p_j = zero;
  n_ij = zero;
  p_ctr = zero;
  temp = zero;
  rhs = zero;
  average = that.average;
  vel_field = that.vel_field;
  drag_rat = that.drag_rat; 
  return *this;
}

void vd_eqn_lazar::get_average_tri() {
  average = f_calc->get_drag_glob();
  if(average < - std::numeric_limits<double>::min())
    average = getAverageEntSize(m, 2);
}

// Calculate the line equilibrium angle. Assume 2pi/3 for interior (would be 2pi/n for higher valency edges where n is the number of grains), pi/2 for edges on cubic sim. cell edges and faces. 
// c_id is cellbase indices rather than gmi indices.
void vd_eqn_lazar::set_ang_equi_cell(int c_id) {

  ext_shell* e_sh = f_calc->get_e_sh();

  shell sh(3,-1); 
  if(e_sh->chk_shell(1, c_id)) {
    sh = e_sh->get_shell(1, c_id);
    if(sh.dim == 3)
      a_equi_map[c_id + 1] = pi/3;
    else if(sh.dim == 2 or sh.dim == 1)
      a_equi_map[c_id + 1] = pi/2;
  }
  else
    a_equi_map[c_id + 1] = pi/3;
}

void vd_eqn_lazar::set_ang_equi_cells() {
  if(f_calc->get_proj() == (int) PROJ_TYPE::EXT_SHELL) {
    for (int c_id = 0; c_id < c_base->get_sz(1); c_id++) {
      if(!c_base->is_free(1, c_id)) {
        set_ang_equi_cell(c_id);
      }
    }
  }
  else {
    for (int c_id = 0; c_id < c_base->get_sz(1); c_id++) {
      if(!c_base->is_free(1, c_id)) {
        a_equi_map[c_id + 1] = pi/3;
      }
    }
  }
}

// Calculate the lenghts and positions of all boundary edges.
void vd_eqn_lazar::mark_edges_global() {
  apf::MeshEntity* e_curr;
  apf::MeshIterator* it = m->begin(1);
  while ((e_curr = m->iterate(it))) {
    apf::ModelEntity* em = m->toModel(e_curr);
    int em_type = m->getModelType(em);
    if(!(em_type == 3)) {
      edge_len[e_curr] = vd_meas_ent(m, e_curr);
      edge_pos[e_curr] = vd_get_pos(m, e_curr);
    }
  }
  m->end(it);
}

// Calculate the lenghts and positions of boundary edges associated with the selected vertices.
void vd_eqn_lazar::mark_edges() {
  vd_set_up(m, &b_verts, &es_edge);

  for(int i = 0; i < es_edge.size(); i++) {
    apf::ModelEntity* mdl = m->toModel(es_edge.at(i));
    int em_type = m->getModelType(mdl);

    if(!(em_type == 3)) {
      edge_len[es_edge.at(i)] = vd_meas_ent(m, es_edge.at(i));
      edge_pos[es_edge.at(i)] = vd_get_pos(m, es_edge.at(i));
    }
  }
}

// Calculate the orientations and positions of all boundary triangles.
void vd_eqn_lazar::mark_tris_global() {
  apf::Up up_tet;
  apf::MeshEntity* tri;
  apf::MeshIterator* it = m->begin(2);
  while ((tri = m->iterate(it))) {
    apf::ModelEntity* em = m->toModel(tri);
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getUp(tri, up_tet);
      em = m->toModel(up_tet.e[0]);
      int c3_tag = m->getModelTag(em);
      tri_tet_in[tri] = c3_tag;
      tri_pos[tri] = vd_get_pos(m, tri);
      tri_norm[tri] = norm_0(vd_area_out(m, tri, c3_tag));
    }
  }
  m->end(it);
}

// Calculate the orientations and positions of boundary tris associated with the selected vertices.
void vd_eqn_lazar::mark_tris() {
  apf::Up up_tet;
  vd_set_up(m, &es_edge, &es_surf);

  for(int i = 0; i < es_surf.size(); i++) {
    apf::MeshEntity* tri = es_surf.at(i);
    apf::ModelEntity* em = m->toModel(tri);
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getUp(tri, up_tet);
      em = m->toModel(up_tet.e[0]);
      int c3_tag = m->getModelTag(em);
      tri_tet_in[tri] = c3_tag;
      tri_pos[tri] = vd_get_pos(m, tri);
      tri_norm[tri] = norm_0(vd_area_out(m, tri, c3_tag));
    }
  }
}

// Calculate the L and M quantities used in local Srolovitz-MacPherson.
void vd_eqn_lazar::collect_orientations() {
  if(b_verts.size() == 0) {
    mark_edges();
    mark_tris();
  }
  else {
    mark_edges_global();
    mark_tris_global();
  }
  set_ang_equi_cells();
}

// Given an entity set of triangles forming a surface, find the approximate area.
apf::Vector3 vd_eqn_lazar::get_eij(std::vector<apf::MeshEntity*> &es_bsurf, int c3_curr) {
  //apf::Up up_tet;
  apf::Vector3 eij(0,0,0);
  for (int i = 0; i < es_bsurf.size(); i++) {
    //double coef = 1.;
    //m->getUp(es_bsurf.at(i), up_tet);
    //if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr)
    //  coef = 1.;
    //else
    //  coef = -1.;
    //eij = eij + vd_area_out(m, es_bsurf.at(i), c3_curr)*coef;
    eij = eij + vd_area_out(m, es_bsurf.at(i), c3_curr);
    // std::cout << i << "The surface is directed " << dA << std::endl;
  }
  return eij;
}

// Calculate the L and M quantities used in local Srolovitz-MacPherson.
double vd_eqn_lazar::vd_dVdt(apf::MeshEntity* vert, int c3_curr, apf::Vector3 &eij) {
  // First, obtain the edges adjacent to the current grain.
  double ln = 0.;

  std::vector<apf::MeshEntity*> es_bedge(0);
  std::vector<apf::MeshEntity*> es_tris(0);
  std::vector<apf::MeshEntity*> es_bsurf(0);
  apf::Up up_tet;

  apf::Vector3 pos_line(0,0,0);
  apf::Vector3 pos_tri1(0,0,0);
  apf::Vector3 pos_tri2(0,0,0);
  apf::Vector3 norm1(0,0,0);
  apf::Vector3 norm2(0,0,0);
  double angle_curr = 0.;

  // Get edges on grain boundary.
  vd_set_up(m, vert, &es_edge);
  vd_chk_edge_grain(m, &es_edge, &es_bedge, c3_curr);

  // Going over the edges, find the surface adjacencies on the grain specified.
  // There are two of them. Calculate the exterior angle, use in L. If triple
  // edge, add the length to M.
  for (int i = 0; i < es_bedge.size(); i++) {
    double len = vd_meas_ent(m, es_bedge.at(i));

    vd_set_up(m, es_bedge.at(i), &es_tris);
    vd_chk_surf_grain(m, &es_tris, &es_bsurf, c3_curr);

    pos_line = edge_pos[es_bedge.at(i)];
    pos_tri1 = tri_pos[es_bsurf.at(0)] - pos_line;
    pos_tri2 = tri_pos[es_bsurf.at(1)] - pos_line; 

    m->getUp(es_bsurf.at(0), up_tet);
    if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr)
      norm1 = tri_norm[es_bsurf.at(0)];
    else
      norm1 = tri_norm[es_bsurf.at(0)]*(-1);
    m->getUp(es_bsurf.at(1), up_tet);
    if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr)
      norm2 = tri_norm[es_bsurf.at(1)];
    else
      norm2 = tri_norm[es_bsurf.at(1)]*(-1);
    angle_curr = vd_ext_angle_n(pos_tri1, pos_tri2, norm1, norm2);

    // L(D)*2pi = sum(e_i*alpha_i)
    ln = ln + len*angle_curr;
    // printf("the exterior angle is %2.2f L_n(D) = %f.\n", angle, lm[0]);
    apf::ModelEntity* em = m->toModel(es_bedge.at(i));
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    // N(D) = sum(e_i*alpha_equi)
    if(em_type == 1)
      ln = ln - len*a_equi_map[em_tag];
  }
  vd_set_up(m, &es_bedge, &es_tris);
  vd_chk_surf_grain(m, &es_tris, &es_bsurf, c3_curr);
  eij = get_eij(es_bsurf, c3_curr);

  return f_calc->vdparam.v_mult*ln/(-2);
}

// Only collect interior triangles. Specially used in 
void vd_eqn_lazar::chk_surf_grain_int(std::vector<apf::MeshEntity*> &es_tris, std::vector<apf::MeshEntity*> &es_bsurf, int c3_curr) {
  es_bsurf.clear();
  es_bsurf.reserve(es_tris.size());
  apf::Up up_tet;

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < es_tris.size(); i++) {
    apf::MeshEntity* tri_curr = es_tris.at(i);
    apf::ModelEntity* em = m->toModel(tri_curr); 
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 2) {
      if(!c_base->get_cell_ext_gmi(em_type, em_tag)) {
        m->getUp(tri_curr, up_tet);
        assert(up_tet.n == 2);
        if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr or m->getModelTag(m->toModel(up_tet.e[1])) == c3_curr)
          es_bsurf.push_back(tri_curr);
      }
    }
  }
}

// Calculate the L and M quantities used in local Srolovitz-MacPherson, only using interior boundaries. Used in calculations involving exterior vertices.
double vd_eqn_lazar::vd_dVdt_int(apf::MeshEntity* vert, int c3_curr, apf::Vector3 &eij) {
  apf::Up up_tet;
  // First, obtain the edges adjacent to the current grain.
  double ln = 0.;

  std::vector<apf::MeshEntity*> es_bedge(0);
  std::vector<apf::MeshEntity*> es_tris(0);
  std::vector<apf::MeshEntity*> es_bsurf(0);

  apf::Vector3 pos_line(0,0,0);
  apf::Vector3 pos_tri1(0,0,0);
  apf::Vector3 pos_tri2(0,0,0);
  apf::Vector3 norm1(0,0,0);
  apf::Vector3 norm2(0,0,0);
  double angle_curr = 0.;

  // Get edges on grain boundary.
  vd_set_up(m, vert, &es_edge);
  vd_chk_edge_grain(m, &es_edge, &es_bedge, c3_curr);

  bool sp_vert = f_calc->chk_vert_special(m, vert);
  // Going over the edges, find the surface adjacencies on the grain specified.
  // There are two of them. Calculate the exterior angle, use in L. If triple
  // edge, add the length to M.
  for (int i = 0; i < es_bedge.size(); i++) {
    double len = vd_meas_ent(m, es_bedge.at(i));

    vd_set_up(m, es_bedge.at(i), &es_tris);
    vd_chk_surf_grain(m, &es_tris, &es_bsurf, c3_curr);

    pos_line = edge_pos[es_bedge.at(i)];
    pos_tri1 = tri_pos[es_bsurf.at(0)] - pos_line;
    pos_tri2 = tri_pos[es_bsurf.at(1)] - pos_line; 

    m->getUp(es_bsurf.at(0), up_tet);
    if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr)
      norm1 = tri_norm[es_bsurf.at(0)];
    else
      norm1 = tri_norm[es_bsurf.at(0)]*(-1);
    m->getUp(es_bsurf.at(1), up_tet);
    if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr)
      norm2 = tri_norm[es_bsurf.at(1)];
    else
      norm2 = tri_norm[es_bsurf.at(1)]*(-1);
    angle_curr = vd_ext_angle_n(pos_tri1, pos_tri2, norm1, norm2);

    // L(D)*2pi = sum(e_i*alpha_i)
    ln = ln + len*angle_curr;
    // printf("the exterior angle is %2.2f L_n(D) = %f.\n", angle, lm[0]);
    apf::ModelEntity* em = m->toModel(es_bedge.at(i));
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);

    //shell sh(3,-1); 
    //assert(e_sh->chk_shell(em_type, em_tag-1));
    //sh = e_sh->get_shell(em_type, em_tag-1);

    if(em_type == 1)
      // N(D) = sum(e_i*alpha_equi)
      ln = ln - len*a_equi_map[em_tag];
      // printf("M_n(D) = %f.\n", lm[1]);
  }
  vd_set_up(m, &es_bedge, &es_tris);
  chk_surf_grain_int(es_tris, es_bsurf, c3_curr);
  eij = get_eij(es_bsurf, c3_curr);

  if(f_calc->get_proj() == (int) PROJ_TYPE::EXT_SHELL and sp_vert) {
    eij = f_calc->get_vec_special(m, vert, eij);
  }

  return f_calc->vdparam.v_mult*ln/(-2);
  // std::cout << "Grain: " << geom << std::endl;
  // printf("L_n(D) = %f M_n(D) = %f.\n",lm[0], lm[1]);
}

apf::Vector3 vd_eqn_lazar::calc_dn_quad(apf::MeshEntity* vert, std::vector<apf::Vector3>& eij, std::vector<double>& dVdt) {
  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  apf::Vector3 dVn_v(0,0,0);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      m_eij[i][j] = eij.at(i)[j];
    }
    dVn_v[i] = dVdt.at(i);
  }

  apf::Vector3 disp(0,0,0);
  disp = apf::invert(m_eij)*dVn_v*3;
  if(std::isnan(disp.getLength()))
    return apf::Vector3(0,0,0);
  else
    return disp;
}

apf::Vector3 vd_eqn_lazar::calc_dn_surf(apf::MeshEntity* vert, std::vector<apf::Vector3>& eij, std::vector<double>& dVdt) {
  double k = dVdt.at(0)/(eij.at(0)*eij.at(0));
  return eij.at(0)*k*3;
  //std::cout << eij[0] << eij[1] << eij[2] << 
  //", dVn: " << dVn[0] << ", " << dVn[1] << ", " << dVn[2] << 
  //std::endl;
}

apf::Vector3 vd_eqn_lazar::calc_dn_trip(apf::MeshEntity* vert, std::vector<apf::Vector3>& eij, std::vector<double>& dVdt) {

  eij.at(2) = vd_cross(eij.at(0),eij.at(1));
  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  apf::Vector3 dVn_v(dVdt.at(0),dVdt.at(1), 0);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      m_eij[i][j] = eij.at(i)[j];
    }
  }

  apf::Vector3 disp(0,0,0);
  disp = apf::invert(m_eij)*dVn_v*3;
  if(std::isnan(disp.getLength()))
    return apf::Vector3(0,0,0);
  else
    return disp;
}

void vd_eqn_lazar::calc_vel_curr(apf::MeshEntity* vert) {
  assert(vel_field);
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  int em_tag = m->getModelTag(em);

  apf::Vector3 temp(0,0,0);
  if(em_type == 3) {
    apf::setVector(vel_field, vert, 0, temp);
    return;
  }
  //std::vector<apf::Vector3> eij(3, apf::Vector3(0,0,0));
  //std::vector<double> dVdt(3, 0.);
  std::vector<apf::Vector3> eijs(3, apf::Vector3(0,0,0));
  std::vector<double> dVdts(3, 0.);

  struct ent_conn e_3c;
  c_base->get_conn_dim_gmi(3, em_type, em_tag, &e_3c);

  if(!c_base->get_cell_ext_gmi(em_type, em_tag)) {
    if(e_3c.conn.size() > 3 - em_type + 1)
      temp = apf::Vector3(0,0,0);
    assert(e_3c.conn.size()  == 3 - em_type + 1);
    for(int i = 0; i < e_3c.conn.size()-1; i++) {
      dVdts.at(i) = vd_dVdt(vert, e_3c.conn.at(i), eijs.at(i));
    }
    if(em_type == 2)
      temp = calc_dn_surf(vert, eijs, dVdts);
    else if(em_type == 1)
      temp = calc_dn_trip(vert, eijs, dVdts);
    else
      temp = calc_dn_quad(vert, eijs, dVdts);
  }
  else {
    if(f_calc->get_proj() == (int) PROJ_TYPE::EXT_SHELL) {

      ext_shell* e_sh = f_calc->get_e_sh();
      assert(e_sh);

      shell sh(3,-1); 
      assert(e_sh->chk_shell(em_type, em_tag-1));
      sh = e_sh->get_shell(em_type, em_tag-1);

      if(sh.dim > em_type) {
        for(int i = 0; i < e_3c.conn.size()-1; i++) {
          dVdts.at(i) = vd_dVdt_int(vert, e_3c.conn.at(i), eijs.at(i));
        }
        if(sh.dim - em_type == 2)
          temp = calc_dn_trip(vert, eijs, dVdts);
        else if(sh.dim - em_type == 1)
          temp = calc_dn_surf(vert, eijs, dVdts);
        else
          temp = apf::Vector3(0,0,0);
      }
      else
        temp = apf::Vector3(0,0,0);

    }
    else {
      temp = apf::Vector3(0,0,0);
    }
  }
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
}

void vd_eqn_lazar::calc_vel() {
  vel_field = m->findField("velocity_field");
  b_verts.clear();
  collect_orientations();

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      calc_vel_curr(vert);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  m->end(it);
  printf("Velocity field finished.\n");
}

void vd_eqn_lazar::calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local) {
  vel_field = m->findField("velocity_field");
  b_verts = *verts;
  collect_orientations();

  apf::MeshEntity* vert;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  for(int i = 0; i < b_verts.size(); i++) {
    vert = b_verts.at(i);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      calc_vel_curr(vert);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  printf("Velocity field finished.\n");
  b_verts.clear();
}


// A single vertex.
void vd_eqn_lazar::vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::Vector3 zero(0,0,0);
  apf::setVector(vel_field, vert, 0, zero);
}

// A single vertex.
apf::Vector3 vd_eqn_lazar::vd_upd_vel_field_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris, bool drag_local) {
  return apf::Vector3(0,0,0);
}

apf::Vector3 vd_eqn_lazar::vd_upd_vel_field_edge(apf::MeshEntity* vert, 
           std::vector<apf::MeshEntity*>* edges, bool drag_local) {
  return apf::Vector3(0,0,0);
}

// A set of vertices while skipping labeled vertices.
void vd_eqn_lazar::vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                       bool drag_local) {
  apf::Up up;

  vel_field = m->findField("velocity_field");

  std::map<apf::MeshEntity*, bool> merg{};
  apf::MeshEntity* v_merg = NULL;
  for(int i = 0; i < vert.size(); i++) {
    if(v_map[vert.at(i)] != NULL) {
      merg[vert.at(i)] = true;
      if(v_merg == NULL)
        v_merg = v_map[vert.at(i)];
    }
  }
  assert(v_merg != NULL);

  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find(vert.begin(), vert.end(), v_merg);
  assert(it != vert.end());
  int i_m = std::distance(vert.begin(), it);

  std::vector<apf::Vector3> rhs(vert.size(), apf::Vector3(0,0,0));
  std::vector<apf::Matrix3x3> m_eij(vert.size(), 
                              apf::Matrix3x3(0,0,0,0,0,0,0,0,0));
  apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

  std::vector<double> average_v(vert.size());
  if(drag_local) {
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = getAverageEntSize(m, vert.at(i), 2);
    }
  }
  else {
    double avg_temp = 0;
    get_average_tri();
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = average;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    m->getUp(vert.at(i), up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    m->getPoint(vert.at(i), 0, p_ctr);

    for(int j = 0; j < es_surf.size(); j++) {
      apf::ModelEntity* em = m->toModel(es_surf.at(j));
      int em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(es_surf.at(j), 0, d_v);

        int v1 = findIn(d_v, 3, vert.at(i));
        m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
        m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;

        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);

        if(!merg[vert.at(i)]) {
          rhs.at(i) = rhs.at(i) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i) = m_eij.at(i) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
        else {
          rhs.at(i_m) = rhs.at(i_m) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i_m) = m_eij.at(i_m) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
      }
    }
    if(!merg[vert.at(i)]) {
      m_eij.at(i) = m_eij.at(i) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
    else {
      m_eij.at(i_m) = m_eij.at(i_m) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    if(!merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i))*rhs.at(i)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
  for(int i = 0; i < vert.size(); i++) {
    if(merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i_m))*rhs.at(i_m)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_eqn_lazar::vd_calc_force(apf::MeshEntity* vert) {
  return apf::Vector3(0,0,0);
}

apf::Vector3 vd_eqn_lazar::vd_calc_force_edge(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges) {
  return apf::Vector3(0,0,0);
}

apf::Vector3 vd_eqn_lazar::vd_calc_force_tri(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris) {
  return apf::Vector3(0,0,0);
}

vd_eqn_lazar::~vd_eqn_lazar() {
}
/*
//////////////////////////////////////
// vd_eqn_lazar_NBC
//////////////////////////////////////
// Update the velocity field at:
// Every boundary vertex.
vd_eqn_lazar_NBC::vd_eqn_lazar_NBC() : vd_eqn_of_motion(),
                               es_edge(0), es_surf(0), 
                               b_verts(0),
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               tri_tet_in{}, edge_len{}, tri_norm{}, 
                               edge_pos{}, tri_pos{}, a_equi_map{},
                               average(-1.), vel_field(NULL), drag_rat(1000) {
}

vd_eqn_lazar_NBC::vd_eqn_lazar_NBC(apf::Mesh2* m_in, cell_base* c_base_in, 
                            field_calc* f_calc_in, 
                            vd_entlist* e_list_in) :
                         vd_eqn_of_motion(m_in, c_base_in, f_calc_in, e_list_in), 
                           es_edge(0), es_surf(0), 
                           b_verts(0),
                           p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                           p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                           tri_tet_in{}, edge_len{}, tri_norm{}, 
                           edge_pos{}, tri_pos{}, a_equi_map{},
                           average(-1.), drag_rat(1000) {
  m = m_in;
  f_calc = f_calc_in;
  e_list = e_list_in;
  c_base = c_base_in;
}

vd_eqn_lazar_NBC::vd_eqn_lazar_NBC(const vd_eqn_lazar_NBC& that) :
                               vd_eqn_of_motion(), 
                               es_edge(0), es_surf(0), 
                               b_verts(0),
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               tri_tet_in{}, edge_len{}, tri_norm{}, 
                               edge_pos{}, tri_pos{}, a_equi_map{},
                               average(-1.), vel_field(NULL), drag_rat(1000) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;
}

vd_eqn_lazar_NBC& vd_eqn_lazar_NBC::operator=(const vd_eqn_lazar_NBC& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;
  c_base = that.c_base;

  es_edge.clear();
  es_surf.clear();

  tri_tet_in.clear();
  edge_len.clear();
  tri_norm.clear(); 
  edge_pos.clear();
  tri_pos.clear();
  a_equi_map.clear();

  b_verts = that.b_verts;

  apf::Vector3 zero(0,0,0);
  p_i = zero;
  p_j = zero;
  n_ij = zero;
  p_ctr = zero;
  temp = zero;
  rhs = zero;
  average = that.average;
  vel_field = that.vel_field;
  drag_rat = that.drag_rat; 
  return *this;
}

void vd_eqn_lazar_NBC::get_average_tri() {
  average = f_calc->get_drag_glob();
  if(average < - std::numeric_limits<double>::min())
    average = getAverageEntSize(m, 2);
}

// Calculate the line equilibrium angle. Assume 2pi/3 for interior (would be 2pi/n for higher valency edges where n is the number of grains), pi/2 for edges on cubic sim. cell edges and faces. 
// c_id is cellbase indices rather than gmi indices.
void vd_eqn_lazar_NBC::set_ang_equi_cell(int c_id) {

  ext_shell* e_sh = f_calc->get_e_sh();

  shell sh(3,-1); 
  if(e_sh->chk_shell(1, c_id)) {
    sh = e_sh->get_shell(1, c_id);
    if(sh.dim == 3)
      a_equi_map[c_id + 1] = pi/3;
    else if(sh.dim == 2 or sh.dim == 1)
      a_equi_map[c_id + 1] = pi/2;
  }
  else
    a_equi_map[c_id + 1] = pi/3;
}

void vd_eqn_lazar_NBC::set_ang_equi_cells() {
  if(f_calc->get_proj() == (int) PROJ_TYPE::EXT_SHELL) {
    for (int c_id = 0; c_id < c_base->get_sz(1); c_id++) {
      if(!c_base->is_free(1, c_id)) {
        set_ang_equi_cell(c_id);
      }
    }
  }
  else {
    for (int c_id = 0; c_id < c_base->get_sz(1); c_id++) {
      if(!c_base->is_free(1, c_id)) {
        a_equi_map[c_id + 1] = pi/3;
      }
    }
  }
}

// Calculate the lenghts and positions of all boundary edges.
void vd_eqn_lazar_NBC::mark_edges_global() {
  apf::MeshEntity* e_curr;
  apf::MeshIterator* it = m->begin(1);
  while ((e_curr = m->iterate(it))) {
    apf::ModelEntity* em = m->toModel(e_curr);
    int em_type = m->getModelType(em);
    if(!(em_type == 3)) {
      edge_len[e_curr] = vd_meas_ent(m, e_curr);
      edge_pos[e_curr] = vd_get_pos(m, e_curr);
    }
  }
  m->end(it);
}

// Calculate the lenghts and positions of boundary edges associated with the selected vertices.
void vd_eqn_lazar_NBC::mark_edges() {
  vd_set_up(m, &b_verts, &es_edge);

  for(int i = 0; i < es_edge.size(); i++) {
    apf::ModelEntity* mdl = m->toModel(es_edge.at(i));
    int em_type = m->getModelType(mdl);

    if(!(em_type == 3)) {
      edge_len[es_edge.at(i)] = vd_meas_ent(m, es_edge.at(i));
      edge_pos[es_edge.at(i)] = vd_get_pos(m, es_edge.at(i));
    }
  }
}

// Calculate the orientations and positions of all boundary triangles.
void vd_eqn_lazar_NBC::mark_tris_global() {
  apf::Up up_tet;
  apf::MeshEntity* tri;
  apf::MeshIterator* it = m->begin(2);
  while ((tri = m->iterate(it))) {
    apf::ModelEntity* em = m->toModel(tri);
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getUp(tri, up_tet);
      em = m->toModel(up_tet.e[0]);
      int c3_tag = m->getModelTag(em);
      tri_tet_in[tri] = c3_tag;
      tri_pos[tri] = vd_get_pos(m, tri);
      tri_norm[tri] = norm_0(vd_area_out(m, tri, c3_tag));
    }
  }
  m->end(it);
}

// Calculate the orientations and positions of boundary tris associated with the selected vertices.
void vd_eqn_lazar_NBC::mark_tris() {
  apf::Up up_tet;
  vd_set_up(m, &es_edge, &es_surf);

  for(int i = 0; i < es_surf.size(); i++) {
    apf::MeshEntity* tri = es_surf.at(i);
    apf::ModelEntity* em = m->toModel(tri);
    int em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getUp(tri, up_tet);
      em = m->toModel(up_tet.e[0]);
      int c3_tag = m->getModelTag(em);
      tri_tet_in[tri] = c3_tag;
      tri_pos[tri] = vd_get_pos(m, tri);
      tri_norm[tri] = norm_0(vd_area_out(m, tri, c3_tag));
    }
  }
}

// Calculate the L and M quantities used in local Srolovitz-MacPherson.
void vd_eqn_lazar_NBC::collect_orientations() {
  if(b_verts.size() == 0) {
    mark_edges();
    mark_tris();
  }
  else {
    mark_edges_global();
    mark_tris_global();
  }
  set_ang_equi_cells();
}

// Given an entity set of triangles forming a surface, find the approximate area.
apf::Vector3 vd_eqn_lazar_NBC::get_eij(std::vector<apf::MeshEntity*> &es_bsurf, int c3_curr) {
  apf::Up up_tet;
  apf::Vector3 eij(0,0,0);
  for (int i = 0; i < es_bsurf.size(); i++) {
    //double coef = 1.;
    //if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr)
    //  coef = 1.;
    //else
    //  coef = -1.;
    //eij = eij + vd_area_out(m, es_bsurf.at(i), c3_curr)*coef;
    eij = eij + vd_area_out(m, es_bsurf.at(i), c3_curr);
    // std::cout << i << "The surface is directed " << dA << std::endl;
  }
  return eij;
}

// Calculate the L and M quantities used in local Srolovitz-MacPherson.
double vd_eqn_lazar_NBC::vd_dVdt(apf::MeshEntity* vert, int c3_curr, apf::Vector3 &eij) {
  // First, obtain the edges adjacent to the current grain.
  double ln = 0.;

  std::vector<apf::MeshEntity*> es_bedge(0);
  std::vector<apf::MeshEntity*> es_tris(0);
  std::vector<apf::MeshEntity*> es_bsurf(0);
  apf::Up up_tet;

  apf::Vector3 pos_line(0,0,0);
  apf::Vector3 pos_tri1(0,0,0);
  apf::Vector3 pos_tri2(0,0,0);
  apf::Vector3 norm1(0,0,0);
  apf::Vector3 norm2(0,0,0);
  double angle_curr = 0.;

  // Get edges on grain boundary.
  vd_set_up(m, vert, &es_edge);
  vd_chk_edge_grain(m, &es_edge, &es_bedge, c3_curr);

  // Going over the edges, find the surface adjacencies on the grain specified.
  // There are two of them. Calculate the exterior angle, use in L. If triple
  // edge, add the length to M.
  for (int i = 0; i < es_bedge.size(); i++) {
    double len = vd_meas_ent(m, es_bedge.at(i));

    vd_set_up(m, es_bedge.at(i), &es_tris);
    vd_chk_surf_grain(m, &es_tris, &es_bsurf, c3_curr);

    pos_line = edge_pos[es_bedge.at(i)];
    pos_tri1 = tri_pos[es_bsurf.at(0)] - pos_line;
    pos_tri2 = tri_pos[es_bsurf.at(1)] - pos_line; 

    m->getUp(es_bsurf.at(0), up_tet);
    if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr)
      norm1 = tri_norm[es_bsurf.at(0)];
    else
      norm1 = tri_norm[es_bsurf.at(0)]*(-1);
    m->getUp(es_bsurf.at(1), up_tet);
    if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr)
      norm2 = tri_norm[es_bsurf.at(1)];
    else
      norm2 = tri_norm[es_bsurf.at(1)]*(-1);
    angle_curr = vd_ext_angle_n(pos_tri1, pos_tri2, norm1, norm2);

    // L(D)*2pi = sum(e_i*alpha_i)
    ln = ln + len*angle_curr;
    // printf("the exterior angle is %2.2f L_n(D) = %f.\n", angle, lm[0]);
    apf::ModelEntity* em = m->toModel(es_bedge.at(i));
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    // N(D) = sum(e_i*alpha_equi)
    if(em_type == 1)
      ln = ln - len*a_equi_map[em_tag];
  }
  vd_set_up(m, &es_bedge, &es_tris);
  vd_chk_surf_grain(m, &es_tris, &es_bsurf, c3_curr);
  eij = get_eij(es_bsurf, c3_curr);

  return f_calc->vdparam.v_mult*ln/(-2);
}

// Only collect interior triangles. Specially used in 
void vd_eqn_lazar_NBC::chk_surf_grain_int(std::vector<apf::MeshEntity*> &es_tris, std::vector<apf::MeshEntity*> &es_bsurf, int c3_curr) {
  es_bsurf.clear();
  es_bsurf.reserve(es_tris.size());
  apf::Up up_tet;

  // Going over the entities of the set, create a list of their upper 
  // adjacencies.
  for (int i = 0; i < es_tris.size(); i++) {
    apf::MeshEntity* tri_curr = es_tris.at(i);
    apf::ModelEntity* em = m->toModel(tri_curr); 
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if (em_type == 2) {
      if(!c_base->get_cell_ext_gmi(em_type, em_tag)) {
        m->getUp(tri_curr, up_tet);
        assert(up_tet.n == 2);
        if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr or m->getModelTag(m->toModel(up_tet.e[1])) == c3_curr)
          es_bsurf.push_back(tri_curr);
      }
    }
  }
}

// Calculate the L and M quantities used in local Srolovitz-MacPherson, only using interior boundaries. Used in calculations involving exterior vertices.
double vd_eqn_lazar_NBC::vd_dVdt_int(apf::MeshEntity* vert, int c3_curr, apf::Vector3 &eij) {
  apf::Up up_tet;
  // First, obtain the edges adjacent to the current grain.
  double ln = 0.;

  std::vector<apf::MeshEntity*> es_bedge(0);
  std::vector<apf::MeshEntity*> es_tris(0);
  std::vector<apf::MeshEntity*> es_bsurf(0);

  apf::Vector3 pos_line(0,0,0);
  apf::Vector3 pos_tri1(0,0,0);
  apf::Vector3 pos_tri2(0,0,0);
  apf::Vector3 norm1(0,0,0);
  apf::Vector3 norm2(0,0,0);
  double angle_curr = 0.;

  // Get edges on grain boundary.
  vd_set_up(m, vert, &es_edge);
  vd_chk_edge_grain(m, &es_edge, &es_bedge, c3_curr);

  // Going over the edges, find the surface adjacencies on the grain specified.
  // There are two of them. Calculate the exterior angle, use in L. If triple
  // edge, add the length to M.
  for (int i = 0; i < es_bedge.size(); i++) {
    double len = vd_meas_ent(m, es_bedge.at(i));

    vd_set_up(m, es_bedge.at(i), &es_tris);
    vd_chk_surf_grain(m, &es_tris, &es_bsurf, c3_curr);

    pos_line = edge_pos[es_bedge.at(i)];
    pos_tri1 = tri_pos[es_bsurf.at(0)] - pos_line;
    pos_tri2 = tri_pos[es_bsurf.at(1)] - pos_line; 

    m->getUp(es_bsurf.at(0), up_tet);
    if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr)
      norm1 = tri_norm[es_bsurf.at(0)];
    else
      norm1 = tri_norm[es_bsurf.at(0)]*(-1);
    m->getUp(es_bsurf.at(1), up_tet);
    if(m->getModelTag(m->toModel(up_tet.e[0])) == c3_curr)
      norm2 = tri_norm[es_bsurf.at(1)];
    else
      norm2 = tri_norm[es_bsurf.at(1)]*(-1);
    angle_curr = vd_ext_angle_n(pos_tri1, pos_tri2, norm1, norm2);

    // L(D)*2pi = sum(e_i*alpha_i)
    ln = ln + len*angle_curr;
    // printf("the exterior angle is %2.2f L_n(D) = %f.\n", angle, lm[0]);
    apf::ModelEntity* em = m->toModel(es_bedge.at(i));
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if(em_type == 1)
      // N(D) = sum(e_i*alpha_equi)
      ln = ln - len*a_equi_map[em_tag];
      // printf("M_n(D) = %f.\n", lm[1]);
  }
  vd_set_up(m, &es_bedge, &es_tris);
  chk_surf_grain_int(es_tris, es_bsurf, c3_curr);
  eij = get_eij(es_bsurf, c3_curr);

  return f_calc->vdparam.v_mult*ln/(-2);
  // std::cout << "Grain: " << geom << std::endl;
  // printf("L_n(D) = %f M_n(D) = %f.\n",lm[0], lm[1]);
}

apf::Vector3 vd_eqn_lazar_NBC::calc_dn_quad(apf::MeshEntity* vert, std::vector<apf::Vector3>& eij, std::vector<double>& dVdt) {
  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  apf::Vector3 dVn_v(0,0,0);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      m_eij[i][j] = eij.at(i)[j];
    }
    dVn_v[i] = dVdt.at(i);
  }

  apf::Vector3 disp(0,0,0);
  disp = apf::invert(m_eij)*dVn_v*3;
  if(std::isnan(disp.getLength()))
    return apf::Vector3(0,0,0);
  else
    return disp;
}

apf::Vector3 vd_eqn_lazar_NBC::calc_dn_surf(apf::MeshEntity* vert, std::vector<apf::Vector3>& eij, std::vector<double>& dVdt) {
  double k = dVdt.at(0)/(eij.at(0)*eij.at(0));
  return eij.at(0)*k*3;
  //std::cout << eij[0] << eij[1] << eij[2] << 
  //", dVn: " << dVn[0] << ", " << dVn[1] << ", " << dVn[2] << 
  //std::endl;
}

apf::Vector3 vd_eqn_lazar_NBC::calc_dn_trip(apf::MeshEntity* vert, std::vector<apf::Vector3>& eij, std::vector<double>& dVdt) {

  eij.at(2) = vd_cross(eij.at(0),eij.at(1));
  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  apf::Vector3 dVn_v(dVdt.at(0),dVdt.at(1), 0);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      m_eij[i][j] = eij.at(i)[j];
    }
  }

  apf::Vector3 disp(0,0,0);
  disp = apf::invert(m_eij)*dVn_v*3;
  if(std::isnan(disp.getLength()))
    return apf::Vector3(0,0,0);
  else
    return disp;
}

void vd_eqn_lazar_NBC::calc_vel_curr(apf::MeshEntity* vert) {
  assert(vel_field);
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  int em_tag = m->getModelTag(em);

  apf::Vector3 temp(0,0,0);
  if(em_type == 3) {
    apf::setVector(vel_field, vert, 0, temp);
    return;
  }
  //std::vector<apf::Vector3> eij(3, apf::Vector3(0,0,0));
  //std::vector<double> dVdt(3, 0.);
  std::vector<apf::Vector3> eijs(3, apf::Vector3(0,0,0));
  std::vector<double> dVdts(3, 0.);

  struct ent_conn e_3c;
  c_base->get_conn_dim_gmi(3, em_type, em_tag, &e_3c);

  if(!c_base->get_cell_ext_gmi(em_type, em_tag)) {
    if(e_3c.conn.size() > 3 - em_type + 1)
      temp = apf::Vector3(0,0,0);
    assert(e_3c.conn.size()  == 3 - em_type + 1);
    for(int i = 0; i < e_3c.conn.size()-1; i++) {
      dVdts.at(i) = vd_dVdt(vert, e_3c.conn.at(i), eijs.at(i));
    }
    if(em_type == 2)
      temp = calc_dn_surf(vert, eijs, dVdts);
    else if(em_type == 1)
      temp = calc_dn_trip(vert, eijs, dVdts);
    else
      temp = calc_dn_quad(vert, eijs, dVdts);
  }
  else {
    if(f_calc->get_proj() == (int) PROJ_TYPE::EXT_SHELL) {

      ext_shell* e_sh = f_calc->get_e_sh();
      assert(e_sh);

      shell sh(3,-1); 
      assert(e_sh->chk_shell(em_type, em_tag-1));
      sh = e_sh->get_shell(em_type, em_tag-1);

      if(sh.dim > em_type) {
        for(int i = 0; i < e_3c.conn.size()-1; i++) {
          dVdts.at(i) = vd_dVdt_int(vert, e_3c.conn.at(i), eijs.at(i));
        }
        if(sh.dim - em_type == 2)
          temp = calc_dn_trip(vert, eijs, dVdts);
        else if(sh.dim - em_type == 1)
          temp = calc_dn_surf(vert, eijs, dVdts);
        else
          temp = apf::Vector3(0,0,0);
      }
      else
        temp = apf::Vector3(0,0,0);

    }
    else {
      temp = apf::Vector3(0,0,0);
    }
  }
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
}

void vd_eqn_lazar_NBC::calc_vel() {
  vel_field = m->findField("velocity_field");
  b_verts.clear();
  collect_orientations();

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      calc_vel_curr(vert);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  m->end(it);
  printf("Velocity field finished.\n");
}

void vd_eqn_lazar_NBC::calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local) {
  vel_field = m->findField("velocity_field");
  b_verts = *verts;
  collect_orientations();

  apf::MeshEntity* vert;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  for(int i = 0; i < b_verts.size(); i++) {
    vert = b_verts.at(i);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      calc_vel_curr(vert);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  printf("Velocity field finished.\n");
  b_verts.clear();
}


// A single vertex.
void vd_eqn_lazar_NBC::vd_upd_vel_field(apf::MeshEntity* vert, bool drag_local) {
  vel_field = m->findField("velocity_field");
  apf::Vector3 zero(0,0,0);
  apf::setVector(vel_field, vert, 0, zero);
}

// A single vertex.
apf::Vector3 vd_eqn_lazar_NBC::vd_upd_vel_field_tri(apf::MeshEntity* vert, 
                           std::vector<apf::MeshEntity*>* tris, bool drag_local) {
  return apf::Vector3(0,0,0);
}

apf::Vector3 vd_eqn_lazar_NBC::vd_upd_vel_field_edge(apf::MeshEntity* vert, 
           std::vector<apf::MeshEntity*>* edges, bool drag_local) {
  return apf::Vector3(0,0,0);
}

// A set of vertices while skipping labeled vertices.
void vd_eqn_lazar_NBC::vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                       bool drag_local) {
  apf::Up up;

  vel_field = m->findField("velocity_field");

  std::map<apf::MeshEntity*, bool> merg{};
  apf::MeshEntity* v_merg = NULL;
  for(int i = 0; i < vert.size(); i++) {
    if(v_map[vert.at(i)] != NULL) {
      merg[vert.at(i)] = true;
      if(v_merg == NULL)
        v_merg = v_map[vert.at(i)];
    }
  }
  assert(v_merg != NULL);

  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find(vert.begin(), vert.end(), v_merg);
  assert(it != vert.end());
  int i_m = std::distance(vert.begin(), it);

  std::vector<apf::Vector3> rhs(vert.size(), apf::Vector3(0,0,0));
  std::vector<apf::Matrix3x3> m_eij(vert.size(), 
                              apf::Matrix3x3(0,0,0,0,0,0,0,0,0));
  apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

  std::vector<double> average_v(vert.size());
  if(drag_local) {
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = getAverageEntSize(m, vert.at(i), 2);
    }
  }
  else {
    double avg_temp = 0;
    get_average_tri();
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = average;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    m->getUp(vert.at(i), up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    m->getPoint(vert.at(i), 0, p_ctr);

    for(int j = 0; j < es_surf.size(); j++) {
      apf::ModelEntity* em = m->toModel(es_surf.at(j));
      int em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(es_surf.at(j), 0, d_v);

        int v1 = findIn(d_v, 3, vert.at(i));
        m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
        m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;

        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);

        if(!merg[vert.at(i)]) {
          rhs.at(i) = rhs.at(i) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i) = m_eij.at(i) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
        else {
          rhs.at(i_m) = rhs.at(i_m) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i_m) = m_eij.at(i_m) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
      }
    }
    if(!merg[vert.at(i)]) {
      m_eij.at(i) = m_eij.at(i) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
    else {
      m_eij.at(i_m) = m_eij.at(i_m) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    if(!merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i))*rhs.at(i)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
  for(int i = 0; i < vert.size(); i++) {
    if(merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i_m))*rhs.at(i_m)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_eqn_lazar_NBC::vd_calc_force(apf::MeshEntity* vert) {
  return apf::Vector3(0,0,0);
}

apf::Vector3 vd_eqn_lazar_NBC::vd_calc_force_edge(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges) {
}

apf::Vector3 vd_eqn_lazar_NBC::vd_calc_force_tri(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris) {
}

vd_eqn_lazar_NBC::~vd_eqn_lazar_NBC() {
}
*/
//////////////////////////////////////
// vd_eqn_lazar_mod: Implementation of eqn. of motion by Lazar et. al. with 
// support for high valency junctions and exterior shell implementation
//////////////////////////////////////

///////////////////////////////////////////////////
//  vd_eqn_mason_mir : Mirror boundary conditions applied on the exterior shell
///////////////////////////////////////////////////

vd_eqn_mason_mir::vd_eqn_mason_mir() : vd_eqn_of_motion(),
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
}

vd_eqn_mason_mir::vd_eqn_mason_mir(apf::Mesh2* m_in, cell_base* c_base_in,
                            field_calc* f_calc_in, 
                            vd_entlist* e_list_in) :
                         vd_eqn_of_motion(m_in, c_base_in, f_calc_in, e_list_in), 
                           es_edge(0), es_surf(0), 
                           p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                           p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                           average(-1.), drag_rat(1000) {
  m = m_in;
  f_calc = f_calc_in;
  e_list = e_list_in;
  c_base = c_base_in;
}

vd_eqn_mason_mir::vd_eqn_mason_mir(const vd_eqn_mason_mir& that) :
                               vd_eqn_of_motion(), 
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;

  c_base = that.c_base;
}

vd_eqn_mason_mir& vd_eqn_mason_mir::operator=(const vd_eqn_mason_mir& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;

  c_base = that.c_base;

  es_edge.clear();
  es_surf.clear();
  apf::Vector3 zero(0,0,0);
  p_i = zero;
  p_j = zero;
  n_ij = zero;
  p_ctr = zero;
  temp = zero;
  rhs = zero;
  average = that.average;
  vel_field = that.vel_field;
  drag_rat = that.drag_rat; 
  return *this;
}

void vd_eqn_mason_mir::get_average_tri() {
  average = f_calc->get_drag_glob();
  if(average < - std::numeric_limits<double>::min())
    average = getAverageEntSize(m, 2);
}

void vd_eqn_mason_mir::calc_vel_curr(apf::MeshEntity* vert) {

  vel_field = m->findField("velocity_field");
  apf::Vector3 zero(0,0,0);
  apf::Vector3 rhs_temp(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);

  ext_shell* e_sh = f_calc->get_e_sh();

  shell sh(3,-1); 
  apf::ModelEntity* mdl = m->toModel(vert);
  int dim_v = m->getModelType(mdl);
  int id_v = m->getModelTag(mdl);

  std::vector<shell> sh_adj(0, shell(3,-1));
  if(e_sh->chk_shell(dim_v, id_v-1)) {
    sh = e_sh->get_shell(dim_v, id_v-1);

    if(sh.dim == 0) {
      apf::setVector(vel_field, vert, 0, zero);
      return;
    }
    else if(sh.dim == 1) {
      struct ent_conn e_conn;
      e_sh->sh_base.get_conn_dim_gmi(2, 1, sh.id, &e_conn);
      sh_adj.resize(e_conn.conn.size());
      for(int i = 0; i < e_conn.conn.size(); i++) {
        sh_adj.at(i) = shell(2, e_conn.conn.at(i));
      }
    }
    else if(sh.dim == 2) {
      sh_adj.resize(1);
      sh_adj.at(0) = sh;
    }
    else {
      sh_adj.resize(0);
    }
  }

  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* em = m->toModel(es_surf.at(i));
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      rhs_temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, es_surf.at(i));
      rhs = rhs + rhs_temp;
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, es_surf.at(i));
      if(!e_sh->chk_shell(em_type, em_tag - 1)) {
        for(int j = 0; j < sh_adj.size(); j++) {
          apf::Vector3 n_dir(0,0,0);
          //n_dir = e_sh->find_dir_mir(sh_adj.at(j).dim, sh_adj.at(j).id, n_ij);
          //rhs = rhs + e_sh->find_dir_mir(sh_adj.at(j).dim, sh_adj.at(j).id, 
          //                                                           rhs_temp);
          n_dir = e_sh->find_dir_mir(sh_adj.at(j), n_ij);
          rhs = rhs + e_sh->find_dir_mir(sh_adj.at(j), rhs_temp);
          m_eij = m_eij + tensorProduct(n_dir, n_dir)*norm*
                                          f_calc->d2(m, es_surf.at(i));
        }
      }
    }
  }
  m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;

  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
}

void vd_eqn_mason_mir::calc_vel() {
  vel_field = m->findField("velocity_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      calc_vel_curr(vert);
    }
  }
  m->end(it);
  printf("Velocity field finished.\n");
}

void vd_eqn_mason_mir::calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local) {
  vel_field = m->findField("velocity_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  for(int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      vd_upd_vel_field(vert);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  printf("Velocity field finished.\n");
}

// A single vertex.
void vd_eqn_mason_mir::vd_upd_vel_field(apf::MeshEntity* vert) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  get_average_tri();
  calc_vel_curr(vert);
}

// A set of vertices while skipping labeled vertices.
void vd_eqn_mason_mir::vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                       bool drag_local) {
  vel_field = m->findField("velocity_field");

  std::map<apf::MeshEntity*, bool> merg{};
  apf::MeshEntity* v_merg = NULL;
  for(int i = 0; i < vert.size(); i++) {
    if(v_map[vert.at(i)] != NULL) {
      merg[vert.at(i)] = true;
      if(v_merg == NULL)
        v_merg = v_map[vert.at(i)];
    }
  }
  assert(v_merg != NULL);

  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find(vert.begin(), vert.end(), v_merg);
  assert(it != vert.end());
  int i_m = std::distance(vert.begin(), it);

  std::vector<apf::Vector3> rhs(vert.size(), apf::Vector3(0,0,0));
  std::vector<apf::Matrix3x3> m_eij(vert.size(), 
                              apf::Matrix3x3(0,0,0,0,0,0,0,0,0));
  apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

  std::vector<double> average_v(vert.size());
  if(drag_local) {
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = getAverageEntSize(m, vert.at(i), 2);
    }
  }
  else {
    double avg_temp = 0;
    get_average_tri();
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = average;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    m->getUp(vert.at(i), up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    m->getPoint(vert.at(i), 0, p_ctr);

    for(int j = 0; j < es_surf.size(); j++) {
      apf::ModelEntity* em = m->toModel(es_surf.at(j));
      int em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(es_surf.at(j), 0, d_v);

        int v1 = findIn(d_v, 3, vert.at(i));
        m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
        m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;

        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);

        if(!merg[vert.at(i)]) {
          rhs.at(i) = rhs.at(i) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i) = m_eij.at(i) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
        else {
          rhs.at(i_m) = rhs.at(i_m) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i_m) = m_eij.at(i_m) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
      }
    }
    if(!merg[vert.at(i)]) {
      m_eij.at(i) = m_eij.at(i) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
    else {
      m_eij.at(i_m) = m_eij.at(i_m) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    if(!merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i))*rhs.at(i)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
  for(int i = 0; i < vert.size(); i++) {
    if(merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i_m))*rhs.at(i_m)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_eqn_mason_mir::vd_calc_force(apf::MeshEntity* vert) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < es_surf.size(); i++) {
    em = m->toModel(es_surf.at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);

      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_mason_mir::vd_calc_force_tri(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  for(int i = 0; i < tris->size(); i++) {
    em = m->toModel(tris->at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {

      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, tris->at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_mason_mir::vd_calc_force_edge(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  std::vector<apf::MeshEntity*> t_set(0);

  for(int i = 0; i < edges->size(); i++) {
    vd_set_up(m, edges->at(i), &t_set);
    for(int j = 0; j < t_set.size(); j++) {
      em = m->toModel(t_set.at(j));
      em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(t_set.at(j), 0, d_v);
        m->getDownward(t_set.at(j), 1, d_e);

        int v1 = findIn(d_v, 3, vert);
        int ii = lookup_tri_e[v1][0];
        int jj = lookup_tri_e[v1][1];
        int iv = lookup_tri_vd[v1][0];
        int jv = lookup_tri_vd[v1][1];
        if(d_e[ii] != edges->at(i)) {
          ii = lookup_tri_e[v1][1];
          jj = lookup_tri_e[v1][0];
          iv = lookup_tri_vd[v1][1];
          jv = lookup_tri_vd[v1][0];
        }
        m->getPoint(d_v[iv], 0, p_i);
        m->getPoint(d_v[jv], 0, p_j);

        rhs = rhs + vd_dir_in_pl(p_ctr, p_i, p_j)* 
                                  f_calc->gam2(m, t_set.at(j))*
                                  (p_i - p_ctr).getLength();
      }
    }
  }
  return rhs/2;
}

vd_eqn_mason_mir::~vd_eqn_mason_mir() {
}

///////////////////////////////////////////////////
//  vd_eqn_mason_mirrot : Mirror and rotation BC applied on the exterior strata
///////////////////////////////////////////////////

void vd_eqn_mason_mirrot::set_rot(std::vector<std::vector<std::vector< 
                                  std::vector<apf::Matrix3x3 > > > >  &rot_in) {
  int sz_dim = rot_in.size();
  assert(sz_dim == 3);
  rot.resize(sz_dim);
  
  for(int dim = 0; dim < sz_dim; dim++) {
    int sz = rot_in.at(dim).size();
    rot.at(dim).resize(sz);
    for(int i = 0; i < sz; i++) {
      int sz_curr = rot_in.at(dim).at(i).size();
      rot.at(dim).at(i).resize(sz_curr);
      for(int j = 0; j < sz_curr; j++) {
        rot.at(dim).at(i) = rot_in.at(dim).at(i);
      }
    }
  }
}

void vd_eqn_mason_mirrot::set_mir(std::vector<std::vector<std::vector< 
                                  std::vector<apf::Vector3 > > > >  &mir_in) {
  int sz_dim = mir_in.size();
  assert(sz_dim == 3);
  mir.resize(sz_dim);
  
  for(int dim = 0; dim < sz_dim; dim++) {
    int sz = mir_in.at(dim).size();
    mir.at(dim).resize(sz);
    for(int i = 0; i < sz; i++) {
      int sz_curr = mir_in.at(dim).at(i).size();
      mir.at(dim).at(i).resize(sz_curr);
      for(int j = 0; j < sz_curr; j++) {
        mir.at(dim).at(i).at(j) = mir_in.at(dim).at(i).at(j);
      }
    }
  }
}

vd_eqn_mason_mirrot::vd_eqn_mason_mirrot() : vd_eqn_of_motion(),
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               rot(0, std::vector< std::vector< v_rot > > (0,
                                      std::vector< v_rot > (0, V_ROT_INIT) ) ),
                               mir(0, std::vector< std::vector< v_mir > > (0,
                                      std::vector< v_mir > (0, V_MIR_INIT) ) ),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
}

vd_eqn_mason_mirrot::vd_eqn_mason_mirrot(apf::Mesh2* m_in, cell_base* c_base_in, 
                            field_calc* f_calc_in, 
                            vd_entlist* e_list_in) :
                         vd_eqn_of_motion(m_in, c_base_in, f_calc_in, e_list_in), 
                           es_edge(0), es_surf(0), 
                           p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                           p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                           rot(0, std::vector< std::vector< v_rot > > (0,
                                  std::vector< v_rot > (0, V_ROT_INIT) ) ),
                           mir(0, std::vector< std::vector< v_mir > > (0,
                                  std::vector< v_mir > (0, V_MIR_INIT) ) ),
                           average(-1.), drag_rat(1000) {
  m = m_in;
  f_calc = f_calc_in;
  e_list = e_list_in;

  c_base = c_base_in;
}

vd_eqn_mason_mirrot::vd_eqn_mason_mirrot(const vd_eqn_mason_mirrot& that) :
                               vd_eqn_of_motion(), 
                               es_edge(0), es_surf(0), 
                               p_i(0,0,0), p_j(0,0,0), n_ij(0,0,0), 
                               p_ctr(0,0,0), temp(0,0,0), rhs(0,0,0),
                               rot(0, std::vector< std::vector< v_rot > > (0,
                                      std::vector< v_rot > (0, V_ROT_INIT) ) ),
                               mir(0, std::vector< std::vector< v_mir > > (0,
                                      std::vector< v_mir > (0, V_MIR_INIT) ) ),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;

  c_base = that.c_base;
}

vd_eqn_mason_mirrot& vd_eqn_mason_mirrot::operator=(const vd_eqn_mason_mirrot& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;

  c_base = that.c_base;

  es_edge.clear();
  es_surf.clear();
  apf::Vector3 zero(0,0,0);
  p_i = zero;
  p_j = zero;
  n_ij = zero;
  p_ctr = zero;
  temp = zero;
  rhs = zero;
  average = that.average;
  vel_field = that.vel_field;
  drag_rat = that.drag_rat; 
  return *this;
}

void vd_eqn_mason_mirrot::get_average_tri() {
  average = f_calc->get_drag_glob();
  if(average < - std::numeric_limits<double>::min())
    average = getAverageEntSize(m, 2);
}

void vd_eqn_mason_mirrot::calc_vel_curr(apf::MeshEntity* vert) {
  apf::Vector3 zero(0,0,0);
  apf::Vector3 rhs_temp(0,0,0);
  rhs = zero;
  temp = zero;

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);

  ext_shell* e_sh = f_calc->get_e_sh();

  shell sh(3,-1); 
  apf::ModelEntity* mdl = m->toModel(vert);
  int dim_v = m->getModelType(mdl);
  int id_v = m->getModelTag(mdl);



  for(int i = 0; i < es_surf.size(); i++) {
    apf::ModelEntity* em = m->toModel(es_surf.at(i));
    int em_type = m->getModelType(em);
    int em_tag = m->getModelTag(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      double norm = n_ij.getLength();
      n_ij = norm_0(n_ij);

      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
      // |t_i| t_i_hat terms are represented by p_i - p_j
      rhs_temp = vd_cross(n_ij, p_i-p_j)*f_calc->gam2(m, es_surf.at(i));
      rhs = rhs + rhs_temp;
      m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, es_surf.at(i));
      if(c_base->get_cell_ext_gmi(em_type, em_tag)) {
        int sz = rot.at(em_type-1).at(em_tag-1).size();

        for(int j = 0; j < sz; j++) {
          int sz_curr = rot.at(em_type-1).at(em_tag-1).at(j).size();
          apf::Vector3 n_dir(0,0,0);
          apf::Vector3 rhs_curr(0,0,0);
          n_dir = n_ij;
          rhs_curr = rhs_temp;
          for(int k = 0; k < sz_curr; k++) {
            apf::Matrix3x3* rot_curr = &rot.at(em_type-1).at(em_tag-1).at(j).at(k);
            n_dir = rot_dir(rot_curr, n_dir);
            rhs_curr = rot_dir(rot_curr, rhs_curr);
          }
          rhs = rhs + rhs_curr;
          m_eij = m_eij + tensorProduct(n_dir, n_dir)*norm*
                                          f_calc->d2(m, es_surf.at(i));
        }

        sz = mir.at(em_type-1).at(em_tag-1).size();
        for(int j = 0; j < sz; j++) {
          int sz_curr = mir.at(em_type-1).at(em_tag-1).at(j).size();
          apf::Vector3 n_dir(0,0,0);
          apf::Vector3 rhs_curr(0,0,0);
          n_dir = n_ij;
          rhs_curr = rhs_temp;
          for(int k = 0; k < sz_curr; k++) {
            apf::Vector3 mir_curr = mir.at(em_type-1).at(em_tag-1).at(j).at(k);
            n_dir = mir_dir_pl(mir_curr, n_dir);
            rhs_curr = mir_dir_pl(mir_curr, rhs_curr);
          }
          rhs = rhs + rhs_curr;
          m_eij = m_eij + tensorProduct(n_dir, n_dir)*norm*
                                          f_calc->d2(m, es_surf.at(i));
        }
      }
    }
  }
  m_eij = m_eij + M_EYE*average*average*f_calc->get_d2_glob()/drag_rat;

  temp = apf::invert(m_eij)*rhs*3;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
}

void vd_eqn_mason_mirrot::calc_vel() {
  vel_field = m->findField("velocity_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    calc_vel_curr(vert);
  }
  m->end(it);
  printf("Velocity field finished.\n");
}

// A single vertex.
void vd_eqn_mason_mirrot::vd_upd_vel_field(apf::MeshEntity* vert) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  get_average_tri();
  calc_vel_curr(vert);
}

void vd_eqn_mason_mirrot::calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local) {
  vel_field = m->findField("velocity_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  for(int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      vd_upd_vel_field(vert);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  printf("Velocity field finished.\n");
}

// A set of vertices while skipping labeled vertices.
void vd_eqn_mason_mirrot::vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                       bool drag_local) {
  vel_field = m->findField("velocity_field");

  std::map<apf::MeshEntity*, bool> merg{};
  apf::MeshEntity* v_merg = NULL;
  for(int i = 0; i < vert.size(); i++) {
    if(v_map[vert.at(i)] != NULL) {
      merg[vert.at(i)] = true;
      if(v_merg == NULL)
        v_merg = v_map[vert.at(i)];
    }
  }
  assert(v_merg != NULL);

  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find(vert.begin(), vert.end(), v_merg);
  assert(it != vert.end());
  int i_m = std::distance(vert.begin(), it);

  std::vector<apf::Vector3> rhs(vert.size(), apf::Vector3(0,0,0));
  std::vector<apf::Matrix3x3> m_eij(vert.size(), 
                              apf::Matrix3x3(0,0,0,0,0,0,0,0,0));
  apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

  std::vector<double> average_v(vert.size());
  if(drag_local) {
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = getAverageEntSize(m, vert.at(i), 2);
    }
  }
  else {
    double avg_temp = 0;
    get_average_tri();
    for(int i = 0; i < vert.size(); i++) {
      average_v.at(i) = average;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    m->getUp(vert.at(i), up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    m->getPoint(vert.at(i), 0, p_ctr);

    for(int j = 0; j < es_surf.size(); j++) {
      apf::ModelEntity* em = m->toModel(es_surf.at(j));
      int em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(es_surf.at(j), 0, d_v);

        int v1 = findIn(d_v, 3, vert.at(i));
        m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
        m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;

        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);

        if(!merg[vert.at(i)]) {
          rhs.at(i) = rhs.at(i) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i) = m_eij.at(i) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
        else {
          rhs.at(i_m) = rhs.at(i_m) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i_m) = m_eij.at(i_m) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
      }
    }
    if(!merg[vert.at(i)]) {
      m_eij.at(i) = m_eij.at(i) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
    else {
      m_eij.at(i_m) = m_eij.at(i_m) + m_eye*average_v.at(i)*average_v.at(i)*f_calc->get_d2_glob()/drag_rat;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    if(!merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i))*rhs.at(i)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
  for(int i = 0; i < vert.size(); i++) {
    if(merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i_m))*rhs.at(i_m)*3;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_eqn_mason_mirrot::vd_calc_force(apf::MeshEntity* vert) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < es_surf.size(); i++) {
    em = m->toModel(es_surf.at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {
      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);

      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_mason_mirrot::vd_calc_force_tri(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  for(int i = 0; i < tris->size(); i++) {
    em = m->toModel(tris->at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {

      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, tris->at(i));
    }
  }
  return rhs/2;
}

apf::Vector3 vd_eqn_mason_mirrot::vd_calc_force_edge(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges) {
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  temp = apf::Vector3(0,0,0);
  rhs = apf::Vector3(0,0,0);

  m->getPoint(vert, 0, p_ctr);
  std::vector<apf::MeshEntity*> t_set(0);

  for(int i = 0; i < edges->size(); i++) {
    vd_set_up(m, edges->at(i), &t_set);
    for(int j = 0; j < t_set.size(); j++) {
      em = m->toModel(t_set.at(j));
      em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(t_set.at(j), 0, d_v);
        m->getDownward(t_set.at(j), 1, d_e);

        int v1 = findIn(d_v, 3, vert);
        int ii = lookup_tri_e[v1][0];
        int jj = lookup_tri_e[v1][1];
        int iv = lookup_tri_vd[v1][0];
        int jv = lookup_tri_vd[v1][1];
        if(d_e[ii] != edges->at(i)) {
          ii = lookup_tri_e[v1][1];
          jj = lookup_tri_e[v1][0];
          iv = lookup_tri_vd[v1][1];
          jv = lookup_tri_vd[v1][0];
        }
        m->getPoint(d_v[iv], 0, p_i);
        m->getPoint(d_v[jv], 0, p_j);

        rhs = rhs + vd_dir_in_pl(p_ctr, p_i, p_j)* 
                                  f_calc->gam2(m, t_set.at(j))*
                                  (p_i - p_ctr).getLength();
      }
    }
  }
  return rhs/2;
}

vd_eqn_mason_mirrot::~vd_eqn_mason_mirrot() {
}


///////////////////////////////////////////////////
//  vd_eqn_kuprat_volm : Mirror and rotation BC applied on the exterior strata
///////////////////////////////////////////////////

void vd_eqn_kuprat_volm::collect_VAL() {
  
  for(int dim = 1; dim < 4; dim++) {
    for(int i = 0; i < e_list->e.at(dim).size(); i++) {
      for(int dim_low = dim; dim_low > 0; dim_low--) {
        for(int j = 0; j < e_list->e.at(dim).at(i).at(dim_low).size(); j++) {
          apf::MeshEntity* e_curr = e_list->e.at(dim).at(i).at(dim_low).at(j);
          VAL.at(dim_low - 1)[e_curr] = vd_meas_ent(m, e_curr);
        }
      }
    }
  }
  for(int dim = 2; dim < 4; dim++) {
    for(int i = 0; i < e_list->e.at(dim).size(); i++) {
      for(int j = 0; j < e_list->e.at(dim).at(i).at(2).size(); j++) {
        apf::MeshEntity* e_curr = e_list->e.at(dim).at(i).at(2).at(j);
        a_norm[e_curr] = norm_0(vd_area_out_n(m,  e_curr));
        a_pos[e_curr] = vd_get_pos(m,  e_curr);
      }
    }
  }
}

double vd_eqn_kuprat_volm::calc_quality(apf::MeshEntity* tet) {
  m->getDownward(tet, 1, d_e);
  m->getDownward(tet, 2, d_s);
  // Sum of squares of lengths and areas of bounding edges and triangles.
  double l_sq = 0;
  double a_sq = 0;
  for(int i = 0; i < 6; i++) {
    double l_curr = VAL.at(0)[d_e[i]];
    l_sq = l_sq + l_curr*l_curr;
  }

  for(int i = 0; i < 4; i++) {
    double a_curr = VAL.at(1)[d_s[i]];
    a_sq = a_sq + a_curr*a_curr;
  }
  double v_sq = VAL.at(2)[tet]*VAL.at(2)[tet];
  if(v_sq < std::numeric_limits<double>::min())
    return Q_TET_UNDEF;
  return l_sq*a_sq/v_sq;
}

void vd_eqn_kuprat_volm::collect_q() {
  int sz = e_list->e.at(3).size();
  for(int i = 0; i < sz; i++) {
    int sz_tet = e_list->e.at(3).at(i).at(3).size();
    for(int j = 0; j < sz_tet; j++) {
      apf::MeshEntity* tet_curr = e_list->e.at(3).at(i).at(3).at(j);
      q_tet[tet_curr] = calc_quality(tet_curr);
    }
  }
}


vd_eqn_kuprat_volm::vd_eqn_kuprat_volm() : vd_eqn_of_motion(),
                              es_edge(0), es_surf(0), es_tet(0), 
                              VAL(3, std::map<apf::MeshEntity*, double>{}),
                              a_norm{}, a_pos{}, q_tet{},
                              p_i(0,0,0),
                              average(-1.), vel_field(NULL), drag_rat(1000) {
}

vd_eqn_kuprat_volm::vd_eqn_kuprat_volm(apf::Mesh2* m_in, cell_base* c_base_in, 
                            field_calc* f_calc_in, 
                            vd_entlist* e_list_in) :
                         vd_eqn_of_motion(m_in, c_base_in, f_calc_in, e_list_in), 
                           es_edge(0), es_surf(0), es_tet(0), 
                           VAL(3, std::map<apf::MeshEntity*, double>{}),
                           a_norm{}, a_pos{}, q_tet{},
                           p_i(0,0,0),
                           average(-1.), drag_rat(1000) {
  m = m_in;
  f_calc = f_calc_in;
  e_list = e_list_in;

  c_base = c_base_in;
}

vd_eqn_kuprat_volm::vd_eqn_kuprat_volm(const vd_eqn_kuprat_volm& that) :
                               vd_eqn_of_motion(), 
                               es_edge(0), es_surf(0), es_tet(0), 
                               VAL(3, std::map<apf::MeshEntity*, double>{}),
                               a_norm{}, a_pos{}, q_tet{},
                               p_i(0,0,0),
                               average(-1.), vel_field(NULL), drag_rat(1000) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;

  c_base = that.c_base;
}

vd_eqn_kuprat_volm& vd_eqn_kuprat_volm::operator=(const vd_eqn_kuprat_volm& that) {
  m = that.m;
  f_calc = that.f_calc;
  e_list = that.e_list;

  c_base = that.c_base;

  es_edge.clear();
  es_surf.clear();
  es_tet.clear();

  apf::Vector3 zero(0,0,0);
  p_i = zero;

  average = that.average;
  vel_field = that.vel_field;
  drag_rat = that.drag_rat; 
  return *this;
}

void vd_eqn_kuprat_volm::get_average_tri() {
  average = f_calc->get_drag_glob();
  if(average < - std::numeric_limits<double>::min())
    average = getAverageEntSize(m, 2);
}

void vd_eqn_kuprat_volm::calc_vel_curr(apf::MeshEntity* vert) {

  apf::Vector3 temp(0,0,0);

  vd_set_up(m, vert, &es_edge);
  vd_set_up(m, &es_edge, &es_surf);
  vd_set_up(m, &es_surf, &es_tet);

  m->getPoint(vert, 0, p_i);

  for(int i = 0; i < es_tet.size(); i++) {
    m->getDownward(es_tet.at(i), 0, d_v);
    m->getDownward(es_tet.at(i), 2, d_s);

    int v_id = findIn(d_v, 4, vert);
    assert(v_id > -1);
    apf::MeshEntity* t_curr = d_s[lookup_v_tri[v_id]];

    //apf::Vector3 a_nc = a_norm[t_curr];
    //apf::Vector3 a_pc = a_pos[t_curr];
    //double q_curr = q_tet[es_tet.at(i)];
    //a_nc*2*q_curr/h_i

    double h_i = a_norm[t_curr] * (p_i - a_pos[t_curr]);
    // Force away from the tetrahedron as the q_tet increases.
    if(std::fabs(h_i) < std::numeric_limits<double>::min())
      temp = temp + (p_i - a_pos[t_curr]) * q_tet[es_tet.at(i)]*200;
    else
      temp = temp + a_norm[t_curr]*2*q_tet[es_tet.at(i)]/h_i;

  }
  apf::setVector(vel_field, vert, 0, temp);
}

void vd_eqn_kuprat_volm::calc_vel() {
  collect_VAL();
  collect_q();

  vel_field = m->findField("velocity_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();
  es_tet.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(em_type == 3) {
      apf::setVector(vel_field, vert, 0, zero);
      calc_vel_curr(vert);
    }
    else
      apf::setVector(vel_field, vert, 0, zero);
  }
  m->end(it);
  printf("Velocity field finished.\n");
}


void vd_eqn_kuprat_volm::calc_vel(std::vector<apf::MeshEntity*>* verts, bool drag_local) {
  vel_field = m->findField("velocity_field");

  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  es_edge.clear();
  es_surf.clear();

  get_average_tri();
  apf::Vector3 zero(0,0,0);

  for(int i = 0; i < verts->size(); i++) {
    vert = verts->at(i);
    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {
      apf::setVector(vel_field, vert, 0, zero);
    }
    else {
      vd_upd_vel_field(vert);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }
  }
  printf("Velocity field finished.\n");
}


// A single vertex.
void vd_eqn_kuprat_volm::vd_upd_vel_field(apf::MeshEntity* vert) {
  vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  get_average_tri();
  calc_vel_curr(vert);
}

// A set of vertices while skipping labeled vertices.
void vd_eqn_kuprat_volm::vd_upd_vel_field_merg(std::vector<apf::MeshEntity*> &vert,
                       std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
                       bool drag_local) {
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_eqn_kuprat_volm::vd_calc_force(apf::MeshEntity* vert) {
  return apf::Vector3(0,0,0);
}

apf::Vector3 vd_eqn_kuprat_volm::vd_calc_force_tri(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* tris) {
  return apf::Vector3(0,0,0);
}

apf::Vector3 vd_eqn_kuprat_volm::vd_calc_force_edge(apf::MeshEntity* vert, 
                       std::vector<apf::MeshEntity*>* edges) {
  return apf::Vector3(0,0,0);
}

vd_eqn_kuprat_volm::~vd_eqn_kuprat_volm() {
}

///////////////////////////////////////////////////
//  OTHER FUNCTIONS
///////////////////////////////////////////////////

// TODO for constants 1/6: OK
// Given dVn(allowed volumetric change), eij (sum of surface area vectors),
// using Lazar algorithm, (1/2 from surface area vectors, so multiply with 3)
// calculate displacements for the surface vertices:
// TODO correct this for surface correction!
apf::Vector3 vd_calc_dn_surf (double* dVn, const apf::Vector3* eij) {
  double k = dVn[0]/(eij[0]*eij[0]);
  return eij[0]*k*3;
  //std::cout << eij[0] << eij[1] << eij[2] << 
  //", dVn: " << dVn[0] << ", " << dVn[1] << ", " << dVn[2] << 
  //std::endl;
}

// calculate displacements for the triple line vertices:
apf::Vector3 vd_calc_dn_trip (double* dVn, const apf::Vector3* eij) {

  // The direction vector:
  apf::Vector3 line_v = vd_cross(eij[0],eij[1]);
  // The normal should be perpendicular to this line.
  double eij0[3];
  double eij1[3];
  double eij2[3];
  eij[0].toArray(eij0);
  eij[1].toArray(eij1);

  apf::Matrix3x3 m_eij(eij0[0],eij0[1],eij0[2],eij1[0],eij1[1],eij1[2],
                       line_v[0],line_v[1],line_v[2]);

  double dVn_mod[3] = {dVn[0],dVn[1],0};
  apf::Vector3 dVn_v(dVn_mod);

  //std::cout << eij[0] << eij[1] << line_v << 
  //", dVn: " << dVn_mod[0] << ", " << dVn_mod[1] << ", " << dVn_mod[2] << 
  //std::endl;

  apf::Vector3 disp(0,0,0);
  disp = apf::invert(m_eij)*dVn_v*3;
  if(std::isnan(disp.getLength()))
    return apf::Vector3(0,0,0);
  else
    return disp;
}

// calculate displacements for the quadruple vertices:
apf::Vector3 vd_calc_dn_quad (double* dVn, const apf::Vector3* eij) {
  // Finding the displacement requires solving three equations for three
  // grains. Use apfMatrix and apfVector:
  // LOOK: weird initialization, but no alternative was visible.
  double eij0[3];
  double eij1[3];
  double eij2[3];
  eij[0].toArray(eij0);
  eij[1].toArray(eij1);
  eij[2].toArray(eij2);
  // printf("eij array to be turned into a matrix...\n");
  apf::Matrix3x3 m_eij(eij0[0],eij0[1],eij0[2],eij1[0],eij1[1],eij1[2],
                       eij2[0],eij2[1],eij2[2]);

  // printf("Matrix obtained. Obtain the dVn vector...\n");
  apf::Vector3 dVn_v(0,0,0);

  //std::cout << eij[0] << eij[1] << eij[2] << 
  //", dVn: " << dVn[0] << ", " << dVn[1] << ", " << dVn[2] << 
  //std::endl;

  // printf("Calculate dn... ");
  // dVn_v = apf::invert(m_eij)*dVn_v*3;
  // printf("Calculated.\n");
  // return dVn_v;
  apf::Vector3 disp(0,0,0);
  disp = apf::invert(m_eij)*dVn_v*3;
  if(std::isnan(disp.getLength()))
    return apf::Vector3(0,0,0);
  else
    return disp;
}

// Calculate the L and M quantities used in local Srolovitz-MacPherson.
void vd_calc_lm(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es_edge, int geom, double* lm) {
  // First, obtain the edges adjacent to the current grain.
  lm[0] = 0;
  lm[1] = 0;

  std::vector<apf::MeshEntity*> es_bedge(0);
  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_bsurf(0);
  apf::Up up;

  // Get edges on grain boundary.
  vd_chk_edge_grain(m, es_edge, &es_bedge, geom);

  // Going over the edges, find the surface adjacencies on the grain specified.
  // There are two of them. Calculate the exterior angle, use in L. If triple
  // edge, add the length to M.
  for (int i = 0; i < es_bedge.size(); i++) {
    apf::MeshElement* ee = createMeshElement(m, es_bedge.at(i));
    double len = measure(ee);
    destroyMeshElement(ee);

    m->getUp(es_bedge.at(i),up);
    copy_ent_set(&es_surf, up);
    vd_chk_surf_grain(m, &es_surf, &es_bsurf, geom);
    //apf::Vector3 norm0 = vd_area_out(m, es_bsurf.at(0), geom);
    //norm0 = norm_0(norm0);
    //apf::Vector3 norm1 = vd_area_out(m, es_bsurf.at(1), geom);
    //norm1 = norm_0(norm1);

    double angle = vd_ext_angle(m, es_bsurf.at(0), es_bsurf.at(1), geom);
    lm[0] = lm[0] + len*angle;
    // printf("the exterior angle is %2.2f L_n(D) = %f.\n", angle, lm[0]);
    apf::ModelEntity* em = m->toModel(es_bedge.at(i));
    int em_type = m->getModelType(em);
    if(em_type == 1)
      lm[1] = lm[1] + len;
      // printf("M_n(D) = %f.\n", lm[1]);

  }
  lm[0] = lm[0]/(4*pi);
  lm[1] = lm[1]/2;

  // std::cout << "Grain: " << geom << std::endl;
  // printf("L_n(D) = %f M_n(D) = %f.\n",lm[0], lm[1]);
}

// Given an entity set of triangles forming a surface, find the approximate area.
void vd_get_eij(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es_surf, int geom, apf::Vector3* eij) {

  //printf("Number of surfaces is %d.\n", es_surf->n);

  double v[3] = {0,0,0};
  eij->fromArray(v);

  apf::Vector3 dA(0,0,0);
  for (int i = 0; i < es_surf->size(); i++) {
    dA = vd_area_out(m, es_surf->at(i), geom);
    // std::cout << i << "The surface is directed " << dA << std::endl;
    eij[0] = eij[0] + dA;
  }
}

// Calculate the dVn and eij quantities used in Lazar algorithm.
// Don't call with internal elements, although it is also defined for that. 
// TODO repeated code, simplify.
void vd_calc_dVn(apf::Mesh2* m, apf::MeshEntity* vert, double mult, 
                 double* dVn, apf::Vector3* eij) {

  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
  int em_tag = m->getModelTag(em);

  double lm[2] = {0,0};

  // Check if the vertex is quad, trip, surf.
  // If quad, we need to evaluate dVn and eij for 3 grains. 
  // If trip, 2, if surf 1.

  // So first obtain the geometry list of elements. 
  // (TODO this can become faster.)

  apf::Up up;

  std::vector<apf::MeshEntity*> es_edge(0);
  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_elem(0);
  // To keep surfaces adjacent to current geometry.
  std::vector<apf::MeshEntity*> es_bsurf(0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);
  vd_set_up(m, &es_surf, &es_elem);

  std::vector<int> geom_l(0);
  find_elem_set_uniq(m, &es_elem, &geom_l);

  //apf::Vector3 point(0,0,0);
  //m->getPoint(vert, 0, point);
  //std::cout << "Position: " << point << std::endl;

  // For each grain, check the following.
  // Calculate L, m, dVn.
  // Calculate eij.
  if (em_type == 0) {
    //printf("A quadruple vertex.\n");
    vd_calc_lm(m, &es_edge, geom_l.at(0), lm);
    dVn[0] = -2*pi*mult*(lm[0]-lm[1]/6);
    vd_calc_lm(m, &es_edge, geom_l.at(1), lm);
    dVn[1] = -2*pi*mult*(lm[0]-lm[1]/6);
    vd_calc_lm(m, &es_edge, geom_l.at(2), lm);
    dVn[2] = -2*pi*mult*(lm[0]-lm[1]/6);
    //std::cout << "Geometries are " << geom_l.geom[0] << geom_l.geom[1] 
    //<< geom_l.geom[2] << std::endl;
    // Once the edge set is finished, calculate dVn.
    // Going over the surf set, check geometric membership type.
    // To obtain the eij, first find surfaces adjacent to selected volumes.
    vd_chk_surf_grain(m, &es_surf, &es_bsurf, geom_l.at(0));
    vd_get_eij(m, &es_bsurf, geom_l.at(0), eij);

    vd_chk_surf_grain(m, &es_surf, &es_bsurf, geom_l.at(1));
    vd_get_eij(m, &es_bsurf, geom_l.at(1), eij+1);

    vd_chk_surf_grain(m, &es_surf, &es_bsurf, geom_l.at(2));
    vd_get_eij(m, &es_bsurf, geom_l.at(2), eij+2);
  }
  else if (em_type == 1) {
    //printf("A triple line vertex.\n");
    vd_calc_lm(m, &es_edge, geom_l.at(0), lm);
    dVn[0] = -2*pi*mult*(lm[0]-lm[1]/6);
    vd_calc_lm(m, &es_edge, geom_l.at(1), lm);
    dVn[1] = -2*pi*mult*(lm[0]-lm[1]/6);

    //std::cout << "Geometries are " << geom_l.geom[0] << geom_l.geom[1] 
    //<< std::endl;
    // To obtain the eij, first find surfaces adjacent to selected volumes.
    vd_chk_surf_grain(m, &es_surf, &es_bsurf, geom_l.at(0));
    vd_get_eij(m, &es_bsurf, geom_l.at(0), eij);

    vd_chk_surf_grain(m, &es_surf, &es_bsurf, geom_l.at(1));
    vd_get_eij(m, &es_bsurf, geom_l.at(1), eij+1);
  }
  else if (em_type == 2) {
    //printf("A surface vertex.\n");
    vd_calc_lm(m, &es_edge, geom_l.at(0), lm);
    dVn[0] = -2*pi*mult*(lm[0]-lm[1]/6);

    //std::cout << "Geometry is " << geom_l.geom[0] << std::endl;
    // To obtain the eij, first find surfaces adjacent to selected volumes.
    vd_chk_surf_grain(m, &es_surf, &es_bsurf, geom_l.at(0));
    vd_get_eij(m, &es_bsurf, geom_l.at(0), eij);
  }
  else {
    // Shouldn't come here if properly executed, but this is here for 
    // completeness.
    //printf("An internal vertex.\n");
  }

}

// Calculate displacements for the vertices, using Lazar algorithm.
void vd_calc_vel_lazar (apf::Mesh2* m, vd_param* vd_par, bool ext) {
  apf::Field* vel_field = m->findField("velocity_field");

  printf("Field initialized.\n");
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;

  while ((vert = m->iterate(it))) {
    apf::ModelEntity* em = m->toModel(vert);
    int em_type = m->getModelType(em);

    double dVn[3] = {0,0,0};
    apf::Vector3 eij[3];
    apf::Vector3 dn(0,0,0);
    eij[0].fromArray(dVn); eij[1].fromArray(dVn); eij[2].fromArray(dVn);
    // printf("dVn and eij initialized.\n");
    // std::cout << eij[0] << eij[1] << eij[2] << ", dVn: " << *dVn << std::endl;
    if((em_type == 3) || (vd_chk_vert_dom(m, vert))) {
      // If on domain boundary or internal, don't do anything.
      // TODO inefficiency in the code due to passing internal vertices.
      // Alternative?
      // printf("On domain or internal: %d.\n", em_type);
    }
    else {
      // em_up = gmi_adjacent(mdl, em_vert, 1);
      // Calculate physically possible volume changes and directions.
      // printf("Geometric membership is %d.\n", em_type);
      vd_calc_dVn(m, vert, vd_par->v_mult, dVn, eij);
      // printf("dVn and eij calculated:\n");
      // Calculate the motion for the vertex.
      if(em_type == 0)
        dn = vd_calc_dn_quad (dVn, eij);
      else if(em_type == 1)
        dn = vd_calc_dn_trip (dVn, eij);
      else if(em_type == 2)
        dn = vd_calc_dn_surf (dVn, eij);
      //else
        //std::cout << "Internal vertex." << std::endl;
      // printf("Displacement calculated...");

    // Set the field.
    //std::cout << "dn is " << dn << std::endl;
    assert(!std::isnan(dn.getLength()));
    apf::setVector(vel_field, vert, 0, dn);
    // printf("Velocity set.\n");
    }
  }
  m->end(it);
  printf("Velocity field finished.\n");
}

// Update the displacement field at a vertex.
void vd_upd_vel_field_lazar(apf::Mesh2* m, apf::MeshEntity* vert, vd_param* vd_par, bool ext) {

  apf::Field* vel_field = m->findField("velocity_field");

  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  double dVn[3] = {0,0,0};
  apf::Vector3 eij[3];
  apf::Vector3 dn(0,0,0);
  dn.fromArray(dVn);
  eij[0].fromArray(dVn); eij[1].fromArray(dVn); eij[2].fromArray(dVn);
  // printf("dVn and eij initialized.\n");
  // std::cout << eij[0] << eij[1] << eij[2] << ", dVn: " << *dVn << std::endl;
  if((em_type == 3) || (vd_chk_vert_dom(m, vert))) {
    // If on domain boundary or internal, don't do anything.
    // TODO inefficiency in the code due to passing internal vertices.
    // Alternative?
    // printf("On domain or internal: %d.\n", em_type);
  }
  else {
    // em_up = gmi_adjacent(mdl, em_vert, 1);
    // Calculate physically possible volume changes and directions.
    //printf("Cell membership is %d.\n", em_type);
    vd_calc_dVn(m, vert, vd_par->v_mult, dVn, eij);
    //printf("dVn and eij calculated:\n");
    // Calculate the motion for the vertex.
    if(em_type == 0)
      dn = vd_calc_dn_quad (dVn, eij);
    else if(em_type == 1)
      dn = vd_calc_dn_trip (dVn, eij);
    else if(em_type == 2)
      dn = vd_calc_dn_surf (dVn, eij);
    //else
      //std::cout << "Internal vertex." << std::endl;
     //printf("Velocity calculated...");
  }
  // Set the field.
  //std::cout << "dn is " << dn << std::endl;
  assert(!std::isnan(dn.getLength()));
  apf::setVector(vel_field, vert, 0, dn);
}

// Calculate the force acting on a vertex. TODO considers energy to be unity.
apf::Vector3 vd_calc_force_mason(apf::Mesh2* m, apf::MeshEntity* vert, 
                                                  field_calc* f_calc) {

  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  std::vector<apf::MeshEntity*> es_edge(0);
  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_2cell(0);
  apf::Up up;

  apf::Downward d_v;
  apf::Downward d_e;

  apf::Vector3 p_i(0,0,0);
  apf::Vector3 p_j(0,0,0);
  apf::Vector3 n_ij(0,0,0);

  apf::Vector3 p_ctr(0,0,0);

  apf::Vector3 temp(0,0,0);

  apf::Vector3 rhs(0,0,0);

  m->getUp(vert, up);
  copy_ent_set(&es_edge, up);
  vd_set_up(m, &es_edge, &es_surf);

  m->getPoint(vert, 0, p_ctr);

  for(int i = 0; i < es_surf.size(); i++) {
    em = m->toModel(es_surf.at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {
      //std::cout << es_surf.at(i) << " 2c" << m->getModelTag(em) 
      //          << std::endl;

      m->getDownward(es_surf.at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][0]], 1, p_i);
      //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][1]], 1, p_j);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      //std::cout << "p_ctr " << p_ctr << " " 
      //          << "p_i " << p_i << " "
      //          << "p_j " << p_j << " " << std::endl;

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      //std::cout << "v_i " << p_i << " "
      //          << "v_j " << p_j << " " << std::endl;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);

      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(i));
    }
  }
  return rhs/2;
}

// Calculate the force acting on a vertex. 
apf::Vector3 vd_calc_force_tri_mason(apf::Mesh2* m, apf::MeshEntity* vert, 
                     std::vector<apf::MeshEntity*>* tris, field_calc* f_calc) {

  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  apf::Vector3 p_i(0,0,0);
  apf::Vector3 p_j(0,0,0);
  apf::Vector3 n_ij(0,0,0);

  apf::Vector3 p_ctr(0,0,0);

  apf::Vector3 temp(0,0,0);

  apf::Vector3 rhs(0,0,0);

  m->getPoint(vert, 0, p_ctr);

  apf::Downward d_v;
  for(int i = 0; i < tris->size(); i++) {
    em = m->toModel(tris->at(i));
    em_type = m->getModelType(em);
    if(em_type == 2) {
      //std::cout << tris->at(i) << " 2c" << m->getModelTag(em) 
      //          << std::endl;

      m->getDownward(tris->at(i), 0, d_v);

      int v1 = findIn(d_v, 3, vert);
      //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][0]], 1, p_i);
      //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][1]], 1, p_j);
      m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
      m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

      //std::cout << "p_ctr " << p_ctr << " " 
      //          << "p_i " << p_i << " "
      //          << "p_j " << p_j << " " << std::endl;

      p_i = p_i - p_ctr;
      p_j = p_j - p_ctr;

      //std::cout << "v_i " << p_i << " "
      //          << "v_j " << p_j << " " << std::endl;

      n_ij = vd_cross(p_i, p_j);
      n_ij = norm_0(n_ij);

      //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
      rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, tris->at(i));
    }
  }
  return rhs/2;
}

// Calculate the force acting on a vertex.
apf::Vector3 vd_calc_force_edge_mason(apf::Mesh2* m, apf::MeshEntity* vert, 
                     std::vector<apf::MeshEntity*>* edges, field_calc* f_calc) {

  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  apf::Vector3 p_i(0,0,0);
  apf::Vector3 p_j(0,0,0);
  apf::Vector3 n_ij(0,0,0);

  apf::Vector3 p_ctr(0,0,0);

  apf::Vector3 temp(0,0,0);

  apf::Vector3 rhs(0,0,0);

  m->getPoint(vert, 0, p_ctr);

  apf::Downward d_v;
  apf::Downward d_e;
  std::vector<apf::MeshEntity*> t_set(0);

  for(int i = 0; i < edges->size(); i++) {
    vd_set_up(m, edges->at(i), &t_set);
    for(int j = 0; j < t_set.size(); j++) {
      em = m->toModel(t_set.at(j));
      em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(t_set.at(j), 0, d_v);
        m->getDownward(t_set.at(j), 1, d_e);

        int v1 = findIn(d_v, 3, vert);
        int ii = lookup_tri_e[v1][0];
        int jj = lookup_tri_e[v1][1];
        int iv = lookup_tri_vd[v1][0];
        int jv = lookup_tri_vd[v1][1];
        if(d_e[ii] != edges->at(i)) {
          ii = lookup_tri_e[v1][1];
          jj = lookup_tri_e[v1][0];
          iv = lookup_tri_vd[v1][1];
          jv = lookup_tri_vd[v1][0];
        }
        m->getPoint(d_v[iv], 0, p_i);
        m->getPoint(d_v[jv], 0, p_j);

        rhs = rhs + vd_dir_in_pl(p_ctr, p_i, p_j)* 
                                  f_calc->gam2(m, t_set.at(j))*
                                  (p_i - p_ctr).getLength();
      }
    }
  }
  return rhs/2;
}

// Calculate displacements for the vertices, using Lazar algorithm.
void vd_calc_vel_lazar (apf::Mesh2* m, field_calc* f_calc) {
  apf::Field* vel_field = m->findField("velocity_field");

  printf("Field initialized.\n");
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;

  while ((vert = m->iterate(it))) {
    apf::ModelEntity* em = m->toModel(vert);
    int em_type = m->getModelType(em);

    double dVn[3] = {0,0,0};
    apf::Vector3 eij[3];
    apf::Vector3 dn(0,0,0);

    eij[0].fromArray(dVn); eij[1].fromArray(dVn); eij[2].fromArray(dVn);
    // printf("dVn and eij initialized.\n");
    // std::cout << eij[0] << eij[1] << eij[2] << ", dVn: " << *dVn << std::endl;
    if((em_type == 3) || (vd_chk_vert_dom(m, vert))) {
      // If on domain boundary or internal, don't do anything.
      // TODO inefficiency in the code due to passing internal vertices.
      // Alternative?
      // printf("On domain or internal: %d.\n", em_type);
    }
    else {
      // em_up = gmi_adjacent(mdl, em_vert, 1);
      // Calculate physically possible volume changes and directions.
      // printf("Geometric membership is %d.\n", em_type);
      vd_calc_dVn(m, vert, f_calc->vdparam.v_mult, dVn, eij);
      // printf("dVn and eij calculated:\n");
      // Calculate the motion for the vertex.
      if(em_type == 0)
        dn = vd_calc_dn_quad (dVn, eij);
      else if(em_type == 1)
        dn = vd_calc_dn_trip (dVn, eij);
      else if(em_type == 2)
        dn = vd_calc_dn_surf (dVn, eij);
      //else
        //std::cout << "Internal vertex." << std::endl;
      // printf("Displacement calculated...");

    // Set the field.
    //std::cout << "dn is " << dn << std::endl;
    assert(!std::isnan(dn.getLength()));
    apf::setVector(vel_field, vert, 0, dn);
    // printf("Velocity set.\n");
    }
  }
  m->end(it);
  printf("Velocity field finished.\n");
}

// Update the displacement field at a vertex.
void vd_upd_vel_field_lazar(apf::Mesh2* m, apf::MeshEntity* vert, field_calc* f_calc, bool drag_local) {

  apf::Field* vel_field = m->findField("velocity_field");

  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  double dVn[3] = {0,0,0};
  apf::Vector3 eij[3];
  apf::Vector3 dn(0,0,0);

  eij[0].fromArray(dVn); eij[1].fromArray(dVn); eij[2].fromArray(dVn);
  // printf("dVn and eij initialized.\n");
  // std::cout << eij[0] << eij[1] << eij[2] << ", dVn: " << *dVn << std::endl;
  if((em_type == 3) || (vd_chk_vert_dom(m, vert))) {
    // If on domain boundary or internal, don't do anything.
    // TODO inefficiency in the code due to passing internal vertices.
    // Alternative?
    // printf("On domain or internal: %d.\n", em_type);
  }
  else {
    // em_up = gmi_adjacent(mdl, em_vert, 1);
    // Calculate physically possible volume changes and directions.
    //printf("Cell membership is %d.\n", em_type);
    vd_calc_dVn(m, vert, f_calc->vdparam.v_mult, dVn, eij);
    //printf("dVn and eij calculated:\n");
    // Calculate the motion for the vertex.
    if(em_type == 0)
      dn = vd_calc_dn_quad (dVn, eij);
    else if(em_type == 1)
      dn = vd_calc_dn_trip (dVn, eij);
    else if(em_type == 2)
      dn = vd_calc_dn_surf (dVn, eij);
    //else
      //std::cout << "Internal vertex." << std::endl;
     //printf("Velocity calculated...");
  }
  // Set the field.
  //std::cout << "dn is " << dn << std::endl;
  assert(!std::isnan(dn.getLength()));
  apf::setVector(vel_field, vert, 0, dn);
}

/*
// Calculate displacements for the surface vertices, using Mason algorithm.
void vd_calc_disp_mason (apf::Mesh2* m, vd_param* vd_par) {
  apf::Field* disp_field = m->findField("displacement_field");

  printf("Field initialized.\n");
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;

  while ((vert = m->iterate(it))) {
    vd_upd_disp_field_mason(m, vert, vd_par);
  }
  m->end(it);
  printf("Displacement field finished.\n");
}
*/

// Calculate velocities for the boundary vertices, using Mason algorithm.
void vd_calc_vel_mason (apf::Mesh2* m, field_calc* f_calc) {
  apf::Field* vel_field = m->findField("velocity_field");

  printf("Field initialized.\n");
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;

  apf::ModelEntity* em;
  int em_type;

  std::vector<apf::MeshEntity*> es_edge(0);
  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_2cell(0);
  apf::Up up;

  apf::Downward d_v;
  apf::Downward d_e;

  apf::Vector3 p_i(0,0,0);
  apf::Vector3 p_j(0,0,0);
  apf::Vector3 n_ij(0,0,0);

  apf::Vector3 p_ctr(0,0,0);
  apf::Vector3 temp(0,0,0);
  apf::Vector3 rhs(0,0,0);

  double average;

  average = f_calc->get_drag_glob();
  if(average < - std::numeric_limits<double>::min())
    average = getAverageEntSize(m, 2);

  while ((vert = m->iterate(it))) {

    em = m->toModel(vert);
    em_type = m->getModelType(em);
    if(f_calc->chk_skip(m, vert)) {

      double z[3] = {0,0,0};
      temp.fromArray(z);
      apf::setVector(vel_field, vert, 0, temp);
    }
    else {
      double z[3] = {0,0,0};
      rhs.fromArray(z);
      temp.fromArray(z);

      apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);
      apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

      m->getUp(vert, up);
      copy_ent_set(&es_edge, up);
      vd_set_up(m, &es_edge, &es_surf);

      m->getPoint(vert, 0, p_ctr);

      for(int i = 0; i < es_surf.size(); i++) {
        em = m->toModel(es_surf.at(i));
        em_type = m->getModelType(em);
        if(em_type == 2) {
          m->getDownward(es_surf.at(i), 0, d_v);

          int v1 = findIn(d_v, 3, vert);
          //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][0]], 1, p_i);
          //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][1]], 1, p_j);
          m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
          m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
          p_i = p_i - p_ctr;
          p_j = p_j - p_ctr;

          n_ij = vd_cross(p_i, p_j);
          double norm = n_ij.getLength();
          n_ij = norm_0(n_ij);

          //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
          // 1/2 |t_i| SUM_ij (n_ij_hat x t_i_hat)gamma(n_ij_hat))
          // |t_i| t_i_hat terms are represented by p_i - p_j
          rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(i));
          m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(i));
        }
      }
      m_eij = m_eij + m_eye*average*average*f_calc->get_d2_glob()/f_calc->get_drag_rat();

      temp = apf::invert(m_eij)*rhs*3;
      //std::cout << "dn is " << temp << std::endl;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert, 0, temp);
    }
    if(f_calc->chk_vert_special(m, vert)) {
      f_calc->fix_vel_special(m, vert);
    }

  }
  m->end(it);
  printf("Velocity field finished.\n");
}

// Update the velocity field at a vertex.
void vd_upd_vel_field_mason(apf::Mesh2* m, apf::MeshEntity* vert, 
                                      field_calc* f_calc, bool drag_local) {

  //m_lookup m_lu;
  //std::cout << "Update " << vert << " using Mason" << std::endl;
  apf::Field* vel_field = m->findField("velocity_field");
  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  std::vector<apf::MeshEntity*> es_edge(0);
  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_2cell(0);
  apf::Up up;

  apf::Downward d_v;
  apf::Downward d_e;

  apf::Vector3 p_i(0,0,0);
  apf::Vector3 p_j(0,0,0);
  apf::Vector3 n_ij(0,0,0);

  apf::Vector3 p_ctr(0,0,0);

  apf::Vector3 temp(0,0,0);

  apf::Vector3 rhs(0,0,0);
/*
  if(em_type == 3) {
    printf("Internal: %d.\n", em_type);
  }
  else if (vd_chk_vert_dom(m, vert)) {
    // If on domain boundary or internal, the volume should be preserved.
    // One way to achieve that is to calculate the motion and project it 
    // on the outer 2-cell planes that the vertex belongs to.
    printf("On domain or internal: %d.\n", em_type);
    if(ext) {
      // If on domain boundary or internal, the volume should be preserved.
      // One way to achieve that is to calculate the motion and project it 
      // on the outer 2-cell planes that the vertex belongs to.
      // printf("On domain or internal: %d.\n", em_type);
      std::vector<apf::MeshEntity*> es_elem;

      vd_upd_vel_field_mason(m, vert, vd_par);
      m->getUp(vert, up);
      copy_ent_set(&es_edge, up);
      vd_set_up(m, &es_edge, &es_surf);
      vd_set_up(m, &es_surf, &es_elem);
      vd_rem_cell(m, &es_elem, f_calc.c3_out, 3);
      apf::getVector(vel_field, vert, 0, temp);
      apf::setVector(vel_field, vert, 0, get_volm_pres(m, vert, 
                                                         &es_elem, temp));
    }
  }
  else {
*/
  double z[3] = {0,0,0};
    rhs.fromArray(z);
    temp.fromArray(z);

    apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);
    apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

    m->getUp(vert, up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    m->getPoint(vert, 0, p_ctr);

    for(int i = 0; i < es_surf.size(); i++) {
      em = m->toModel(es_surf.at(i));
      em_type = m->getModelType(em);
      if(em_type == 2) {
        //std::cout << es_surf.at(i) << " 2c" << m->getModelTag(em) 
        //          << std::endl;

        m->getDownward(es_surf.at(i), 0, d_v);

        int v1 = findIn(d_v, 3, vert);
        //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][0]], 1, p_i);
        //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][1]], 1, p_j);
        m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
        m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

        //std::cout << "p_ctr " << p_ctr << " " 
        //          << "p_i " << p_i << " "
        //          << "p_j " << p_j << " " << std::endl;

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;

        //std::cout << "v_i " << p_i << " "
        //          << "v_j " << p_j << " " << std::endl;

        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);

        //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
        rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(i));
        //assert(!std::isnan(n_ij.getLength());
        //assert(!std::isnan(rhs.getLength());
        //std::cout << "Tri " << es_surf.at(i) 
        //          << " nij: " << n_ij
        //          << " vd_cross: " << vd_cross(n_ij, p_i-p_j)
        //          << " norm: " << norm
        //          << " rhs: " << rhs;
        m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                    f_calc->d2(m, es_surf.at(i));
      }
    }
    double average;
    if(drag_local)
      average = getAverageEntSize(m, vert, 2);
    else {
      average = f_calc->get_drag_glob();
      if(average < - std::numeric_limits<double>::min())
        average = getAverageEntSize(m, 2);
    }
    //std::cout << std::endl;
    //std::cout << "m_eij " << m_eij << std::endl;
    //std::cout << "added m_eye " << m_eij << std::endl;

    m_eij = m_eij + m_eye*average*average*f_calc->get_d2_glob()/f_calc->get_drag_rat();
    //m_eij = apf::invert(m_eij);
    //std::cout << "inverted " << m_eij << std::endl;
    //std::cout << "rhs " << rhs << std::endl;
    //temp = m_eij*rhs;
    //std::cout << "temp " << temp << std::endl;
    //temp = temp*3*f_calc->vdparam.v_mult;
    //std::cout << "scaled " << temp << std::endl;
    temp = apf::invert(m_eij)*rhs*3;
    std::cout << "dn is " << temp << std::endl;
    assert(!std::isnan(temp.getLength()));
    apf::setVector(vel_field, vert, 0, temp);
  //}
}

// Update the velocity field at a vertex.
void vd_upd_vel_field_mason_merg(apf::Mesh2* m, 
      std::vector<apf::MeshEntity*> &vert,
      std::map<apf::MeshEntity*, apf::MeshEntity*> &v_map,
      field_calc* f_calc, bool drag_local) {
  //m_lookup m_lu;
  //std::cout << "Update " << vert << " using Mason" << std::endl;
  apf::Field* vel_field = m->findField("velocity_field");

  std::map<apf::MeshEntity*, bool> merg{};
  apf::MeshEntity* v_merg = NULL;
  for(int i = 0; i < vert.size(); i++) {
    if(v_map[vert.at(i)] != NULL) {
      merg[vert.at(i)] = true;
      if(v_merg == NULL)
        v_merg = v_map[vert.at(i)];
    }
  }
  assert(v_merg != NULL);

  std::vector<apf::MeshEntity*>::iterator it;
  it = std::find(vert.begin(), vert.end(), v_merg);
  assert(it != vert.end());
  int i_m = std::distance(vert.begin(), it);

  std::vector<apf::MeshEntity*> es_edge(0);
  std::vector<apf::MeshEntity*> es_surf(0);
  std::vector<apf::MeshEntity*> es_2cell(0);
  apf::Up up;

  apf::Downward d_v;
  apf::Downward d_e;

  apf::Vector3 p_i(0,0,0);
  apf::Vector3 p_j(0,0,0);
  apf::Vector3 n_ij(0,0,0);

  apf::Vector3 p_ctr(0,0,0);

  apf::Vector3 temp(0,0,0);

  std::vector<apf::Vector3> rhs(vert.size(), apf::Vector3(0,0,0));
  std::vector<apf::Matrix3x3> m_eij(vert.size(), 
                              apf::Matrix3x3(0,0,0,0,0,0,0,0,0));
  apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

  std::vector<double> average(vert.size());
  if(drag_local) {
    for(int i = 0; i < vert.size(); i++) {
      average.at(i) = getAverageEntSize(m, vert.at(i), 2);
    }
  }
  else {
    double avg_temp = 0;
    if(f_calc->get_drag_glob() < - std::numeric_limits<double>::min())
      avg_temp = getAverageEntSize(m, 2);
    else
      avg_temp = f_calc->get_drag_glob();
    for(int i = 0; i < vert.size(); i++) {
      average.at(i) = avg_temp;
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    m->getUp(vert.at(i), up);
    copy_ent_set(&es_edge, up);
    vd_set_up(m, &es_edge, &es_surf);

    m->getPoint(vert.at(i), 0, p_ctr);

    for(int j = 0; j < es_surf.size(); j++) {
      apf::ModelEntity* em = m->toModel(es_surf.at(j));
      int em_type = m->getModelType(em);
      if(em_type == 2) {
        //std::cout << es_surf.at(i) << " 2c" << m->getModelTag(em) 
        //          << std::endl;

        m->getDownward(es_surf.at(j), 0, d_v);

        int v1 = findIn(d_v, 3, vert.at(i));
        //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][0]], 1, p_i);
        //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][1]], 1, p_j);
        m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
        m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);

        //std::cout << "p_ctr " << p_ctr << " " 
        //          << "p_i " << p_i << " "
        //          << "p_j " << p_j << " " << std::endl;

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;

        //std::cout << "v_i " << p_i << " "
        //          << "v_j " << p_j << " " << std::endl;

        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);

        //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
        if(!merg[vert.at(i)]) {
          rhs.at(i) = rhs.at(i) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i) = m_eij.at(i) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
        else {
          rhs.at(i_m) = rhs.at(i_m) + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, es_surf.at(j));
          m_eij.at(i_m) = m_eij.at(i_m) + tensorProduct(n_ij, n_ij)*norm*
                                              f_calc->d2(m, es_surf.at(j));
        }
      }
    }
    if(!merg[vert.at(i)]) {
      m_eij.at(i) = m_eij.at(i) + m_eye*average.at(i)*average.at(i)*f_calc->get_d2_glob()/f_calc->get_drag_rat();
    }
    else {
      m_eij.at(i_m) = m_eij.at(i_m) + m_eye*average.at(i)*average.at(i)*f_calc->get_d2_glob()/f_calc->get_drag_rat();
    }
  }

  for(int i = 0; i < vert.size(); i++) {
    if(!merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i))*rhs.at(i)*3;
      std::cout << "dn is " << temp << std::endl;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
  for(int i = 0; i < vert.size(); i++) {
    if(merg[vert.at(i)]) {
      temp = apf::invert(m_eij.at(i_m))*rhs.at(i_m)*3;
      std::cout << "dn is " << temp << std::endl;
      assert(!std::isnan(temp.getLength()));
      apf::setVector(vel_field, vert.at(i), 0, temp);
    }
  }
}

// Update the velocity field at a vertex, using only a set of triangles.
void vd_upd_vel_field_mason_tri(apf::Mesh2* m, apf::MeshEntity* vert, std::vector<apf::MeshEntity*>* tri_set, field_calc* f_calc, bool drag_local) {
  //m_lookup m_lu;
  //std::cout << "Update " << vert << " using Mason" << std::endl;
  apf::Field* vel_field = m->findField("velocity_field");

  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);
/*
  if(em_type == 3) {
    printf("Internal: %d.\n", em_type);
  }
  else if (vd_chk_vert_dom(m, vert)) {
    // If on domain boundary or internal, the volume should be preserved.
    // One way to achieve that is to calculate the motion and project it 
    // on the outer 2-cell planes that the vertex belongs to.
    printf("On domain or internal: %d.\n", em_type);
    if(ext) {
      // If on domain boundary or internal, the volume should be preserved.
      // One way to achieve that is to calculate the motion and project it 
      // on the outer 2-cell planes that the vertex belongs to.
      // printf("On domain or internal: %d.\n", em_type);
      std::vector<apf::MeshEntity*> es_elem;

      vd_upd_vel_field_mason(m, vert, vd_par);
      m->getUp(vert, up);
      copy_ent_set(&es_edge, up);
      vd_set_up(m, &es_edge, &es_surf);
      vd_set_up(m, &es_surf, &es_elem);

      vd_rem_cell(m, &es_elem, f_calc.c3_out, 3);
      apf::getVector(vel_field, vert, 0, temp);
      apf::setVector(vel_field, vert, 0, get_volm_pres(m, vert, 
                                                         &es_elem, temp));
    }
  }
  else {
*/
    apf::Downward d_v;
    apf::Downward d_e;

    apf::Vector3 p_i(0,0,0);
    apf::Vector3 p_j(0,0,0);
    apf::Vector3 n_ij(0,0,0);

    apf::Vector3 p_ctr(0,0,0);

    apf::Vector3 temp(0,0,0);

    apf::Vector3 rhs(0,0,0);
    double z[3] = {0,0,0};
    rhs.fromArray(z);
    temp.fromArray(z);

    apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);
    apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

    m->getPoint(vert, 0, p_ctr);

    //std::cout << "num of tri " << tri_set->size() << std::endl;

    for(int i = 0; i < tri_set->size(); i++) {
      em = m->toModel(tri_set->at(i));
      em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(tri_set->at(i), 0, d_v);

        int v1 = findIn(d_v, 3, vert);
        //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][0]], 1, p_i);
        //m->getPoint(d_v[m_lu.lookup_tri_vd[v1][1]], 1, p_j);
        if(v1 > -1) {

          //std::cout << "i " << d_v[lookup_tri_vd[v1][0]] << std::endl; 
          //std::cout << "j " << d_v[lookup_tri_vd[v1][1]] << std::endl; 

          m->getPoint(d_v[lookup_tri_vd[v1][0]], 0, p_i);
          m->getPoint(d_v[lookup_tri_vd[v1][1]], 0, p_j);
          p_i = p_i - p_ctr;
          p_j = p_j - p_ctr;

          //std::cout << "p_ctr " << p_ctr << std::endl; 
          //std::cout << "p_i " << p_i << " p_j " << p_j << std::endl; 

          n_ij = vd_cross(p_i, p_j);
          double norm = n_ij.getLength();
          n_ij = norm_0(n_ij);

          rhs = rhs + vd_cross(n_ij, p_i-p_j)* 
                                f_calc->gam2(m, tri_set->at(i));
          //rhs = rhs + vd_cross(n_ij, p_i) + vd_cross(n_ij, p_j);
          //std::cout << "Tri " << tri_set->at(i) 
          //          << " nij: " << n_ij
          //          << " rhs: " << vd_cross(n_ij, p_i-p_j);
          m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                          f_calc->d2(m, tri_set->at(i));
        }
      }
    }
    double average;
    if(drag_local)
      average = getAverageEntSize(m, vert, 2);
    else {
      average = f_calc->get_drag_glob();
      if(average < - std::numeric_limits<double>::min())
        average = getAverageEntSize(m, 2);
    }

    m_eij = m_eij + m_eye*average*average*f_calc->get_d2_glob()/f_calc->get_drag_rat();

    temp = apf::invert(m_eij)*rhs*3;
    std::cout << "dn is " << temp << std::endl;
    assert(!std::isnan(temp.getLength()));
    apf::setVector(vel_field, vert, 0, temp);
  //}

/*
    if(ext_0cell) {
      e_lens.vt.at(i) = rem_norm_comp(e_lens.m, &tri_slice.at(i), 
                                  e_lens.vt.at(i));
    }
*/
}

// Update the velocity field at a vertex, using only a set of triangles.
void vd_upd_vel_field_mason_edge(apf::Mesh2* m, apf::MeshEntity* vert, std::vector<apf::MeshEntity*>* edge_set, field_calc* f_calc, bool drag_local) {
  //m_lookup m_lu;
  //std::cout << "Update " << vert << " using Mason" << std::endl;
  apf::Field* vel_field = m->findField("velocity_field");

  apf::ModelEntity* em = m->toModel(vert);
  int em_type = m->getModelType(em);

  apf::Downward d_v;
  apf::Downward d_e;

  apf::Vector3 p_i(0,0,0);
  apf::Vector3 p_j(0,0,0);
  apf::Vector3 n_ij(0,0,0);

  apf::Vector3 p_ctr(0,0,0);

  apf::Vector3 temp(0,0,0);

  apf::Vector3 rhs(0,0,0);
  double z[3] = {0,0,0};
  rhs.fromArray(z);
  temp.fromArray(z);

  apf::Matrix3x3 m_eij(0,0,0,0,0,0,0,0,0);
  apf::Matrix3x3 m_eye(1,0,0,0,1,0,0,0,1);

  m->getPoint(vert, 0, p_ctr);

  std::vector<apf::MeshEntity*> t_set(0);

  for(int i = 0; i < edge_set->size(); i++) {
    vd_set_up(m, edge_set->at(i), &t_set);
    for(int j = 0; j < t_set.size(); j++) {
      em = m->toModel(t_set.at(j));
      em_type = m->getModelType(em);
      if(em_type == 2) {
        m->getDownward(t_set.at(j), 0, d_v);
        m->getDownward(t_set.at(j), 1, d_e);

        int v1 = findIn(d_v, 3, vert);
        int ii = lookup_tri_e[v1][0];
        int jj = lookup_tri_e[v1][1];
        int iv = lookup_tri_vd[v1][0];
        int jv = lookup_tri_vd[v1][1];
        if(d_e[ii] != edge_set->at(i)) {
          ii = lookup_tri_e[v1][1];
          jj = lookup_tri_e[v1][0];
          iv = lookup_tri_vd[v1][1];
          jv = lookup_tri_vd[v1][0];
          assert(d_e[ii] == edge_set->at(i));
        }
        m->getPoint(d_v[iv], 0, p_i);
        m->getPoint(d_v[jv], 0, p_j);

        rhs = rhs + vd_dir_in_pl(p_ctr, p_i, p_j)* 
                                  f_calc->gam2(m, t_set.at(j))*
                                  (p_i - p_ctr).getLength();

        p_i = p_i - p_ctr;
        p_j = p_j - p_ctr;
        n_ij = vd_cross(p_i, p_j);
        double norm = n_ij.getLength();
        n_ij = norm_0(n_ij);
        m_eij = m_eij + tensorProduct(n_ij, n_ij)*norm*
                                        f_calc->d2(m, t_set.at(j))/6;
      }
    }
  }
  double average;
  if(drag_local)
    average = getAverageEntSize(m, vert, 2);
  else {
    average = f_calc->get_drag_glob();
    if(average < - std::numeric_limits<double>::min())
      average = getAverageEntSize(m, 2);
  }

  m_eij = m_eij + m_eye*average*average*f_calc->get_d2_glob()/f_calc->get_drag_rat();

  temp = apf::invert(m_eij)*rhs;
  std::cout << "dn is " << temp << std::endl;
  assert(!std::isnan(temp.getLength()));
  apf::setVector(vel_field, vert, 0, temp);
}

// TODO Reverting positions can be done by apf::field. So some of the 
// functionality is obsolete.
// Apply the displacements determined by the velocity field, multiplied by time
// step. Can be used to revert positions.
void vd_apply_vel_field(apf::Mesh2* m, double mult) {
  apf::Field* vel_field = m->findField("velocity_field");

  // printf("Field initialized.\n");
  apf::MeshEntity* vert;
  apf::MeshIterator* it = m->begin(0);

  while ((vert = m->iterate(it))) {

    apf::Vector3 point(0,0,0);
    apf::Vector3 dn(0,0,0);

    apf::getVector(vel_field, vert, 0, dn);
    m->getPoint(vert, 0, point);
    m->setPoint(vert, 0, point + dn*mult);
    // apf::setVector(disp_field, vert, 0, dn*0);
  }
  m->end(it);
  printf("Positions updated.\n");
}

void vd_apply_vel_field(apf::Mesh2* m, std::vector<apf::MeshEntity* >* es, double mult) {
  apf::Field* vel_field = m->findField("velocity_field");

  for (int i = 0; i < es->size(); i++) {

    apf::Vector3 point(0,0,0);
    apf::Vector3 dn(0,0,0);

    apf::getVector(vel_field, es->at(i), 0, dn);
    m->getPoint(es->at(i), 0, point);
    m->setPoint(es->at(i), 0, point + dn*mult);
    // apf::setVector(disp_field, vert, 0, dn*0);
  }
  printf("Positions updated.\n");
}
/*
double vd_find_min_t(apf::Mesh2* m) {
  apf::Field* vel_field = m->findField("velocity_field");

  // printf("Field initialized.\n");
  apf::MeshEntity* vert;
  std::vector<apf::MeshEntity*> e_v(0);
  std::vector<apf::MeshEntity*> e_e(0);
  std::vector<apf::MeshEntity*> e_s(0);
  std::vector<apf::MeshEntity*> e_t(0);

  apf::Vector3 t_pos;
  apf::Vector3 v_pos;
  apf::Vector3 dn;

  double t_min = -1;
  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    if(!(vd_chk_vert_dom(m, vert))) {
      apf::getVector(vel_field, vert, 0, dn);

      m->getPoint(vert, 0, v_pos);

      init_ent_set(&e_v, vert);
      vd_set_up(m, &e_v, &e_e);
      vd_set_up(m, &e_e, &e_s);
      vd_set_up(m, &e_s, &e_t);

      // The tetrahedra vertex is moving into.
      int tet_id = 0;
      double comp = 0.;

      //std::cout << "v_pos " << v_pos << " dn " << dn << std::endl;
      for(int i = 0; i < e_t.size(); i++) {

        t_pos = getLinearCentroid(m, e_t.at(i));
        t_pos = t_pos-v_pos;
        t_pos = norm_0(t_pos);
        if(t_pos*dn > comp)
          tet_id = i;
        //std::cout << "\tt_pos " << t_pos << std::endl;

      }

      t_pos = getLinearCentroid(m, e_t.at(tet_id));
      t_pos = t_pos-v_pos;

      double dist = t_pos.getLength();
      //m->getPoint(vert, 0, point);
      double t = dist/(norm_0(t_pos)*dn);
      t = fabs(t);
      //std::cout << "t_curr " << t << std::endl;

      //std::cout << vert << ", vel: " << vel << std::endl;
      if (t_min < 0 or t < t_min) {
        t_min = t;
      }
    }
    // apf::setVector(disp_field, vert, 0, dn*0);
  }
  m->end(it);

  printf("Minimum time is %.3f.\n", t_min);
  return t_min;
}

*/
// Tetrahedron based time step estimation:
double vd_calc_tet_t(apf::Mesh2* m, apf::MeshEntity* tet, 
                                              double &V, double &C) {
  int v_sense[4] = {-1, 1, -1, 1};

  int ent_type = m->getType(tet);
  int d = m->typeDimension[ent_type];
  assert(d == 3);

  apf::Vector3 r(0,0,0);
  std::vector<apf::Vector3> pos(4, apf::Vector3(0,0,0));
  C = -1;
  V = 0;

  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);
  apf::Downward d_v;
  m->getDownward(tet, 0, d_v);
          
  for(int i = 0; i < 4; i++) {
    m->getPoint(d_v[i], 0, pos.at(i));
  }

  V = vd_trip_prod(pos.at(3) - pos.at(0), 
                   pos.at(1) - pos.at(0), 
                    pos.at(2) - pos.at(0));

  double C_curr = 0;
  double C_tot = 0;

  for(int i0 = 0; i0 < 4; i0++) {
    int i1 = (i0+1) %4;
    int i2 = (i0+2) %4;
    int i3 = (i0+3) %4;
    apf::getVector(vel_field, d_v[i0], 0, r);
    // vd_trip_prod(03, 01, 02): (01 x 02)*03
    // If volume is decreasing, consider the rate of decrease:
    C_curr = vd_trip_prod(r, pos.at(i2) - pos.at(i1), 
                    pos.at(i3) - pos.at(i1))*v_sense[i0]*(-1);
    if(C_curr > std::numeric_limits<double>::min()) {
      if(C < -std::numeric_limits<double>::min() or C_curr < C) {
        C = C_curr;
      }
    }
    C_tot = C_tot + C_curr;
  }
  //std::cout << "V: " << V << " C_min: " << C 
  //          << " C_tot " << C_tot << std::endl;
  if(C > std::numeric_limits<double>::min())
    return V/C;
  else
    return -1;

}

// Calculate the time step required to invert the tetrahedron, based on the 
// first order expansion of the derivative of volume w.r.t. velocities of the
// vertices.
double vd_find_min_t_tet(apf::Mesh2* m, bool ext) {
  // printf("Field initialized.\n");
  apf::MeshEntity* tet;

  double t_min = -1;
  double t_curr = -1;
  double V = 0;
  double C = 0;

  apf::MeshIterator* it = m->begin(3);
  while ((tet = m->iterate(it))) {
    t_curr = vd_calc_tet_t(m, tet, V, C);
    if (t_curr > std::numeric_limits<double>::min() and 
        (t_curr < t_min or t_min < -std::numeric_limits<double>::min())) {
      t_min = t_curr;
      //std::cout << tet << " V " << V
      //          << " C " << C
      //          << " t " << t_curr << "\n"
      //          << " t_pos " << getLinearCentroid(m, tet)
      //          << std::endl;

      apf::ModelEntity* mdl = m->toModel (tet);
      int c_type = m->getModelType (mdl);
      int c_tag = m->getModelTag (mdl);
      //std::cout << c_type << "c" << c_tag 
      //          << std::endl;
    }
  }
  m->end(it);

  assert(t_min > std::numeric_limits<double>::min());
  printf("Minimum time is %.9f.\n", t_min);
  return t_min;
}

// Vertex based time step estimation:
// This is actually the correct way of calculating the min time, by using all
// triangles across. So any triangle the vertex moves towards is checked for
// inversion.
double vd_find_min_t_all(apf::Mesh2* m, bool ext, double t_set) {
  apf::Field* vel_field = m->findField("velocity_field");

  // printf("Field initialized.\n");
  apf::MeshEntity* vert;
  std::vector<apf::MeshEntity*> e_v(0);
  std::vector<apf::MeshEntity*> e_e(0);
  std::vector<apf::MeshEntity*> e_s(0);
  std::vector<apf::MeshEntity*> e_t(0);

  apf::Vector3 t_pos(0,0,0);
  apf::Vector3 v_pos(0,0,0);
  apf::Vector3 dn(0,0,0);
  apf::Vector3 sn(0,0,0);

  double t_min = t_set;
  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {

    if(ext or !(vd_chk_vert_dom(m, vert))) {
      apf::getVector(vel_field, vert, 0, dn);

      m->getPoint(vert, 0, v_pos);

      init_ent_set(&e_v, vert);
      vd_set_up(m, &e_v, &e_e);
      vd_set_up(m, &e_e, &e_s);
      vd_set_up(m, &e_s, &e_t);

      // The tetrahedra vertex is moving into.
      int tet_id = 0;
      double comp = 0.;

      //std::cout << "v_pos " << v_pos << " dn " << dn << std::endl;
      for(int i = 0; i < e_t.size(); i++) {
        apf::Downward d_t;
        apf::Downward d_v;
        m->getDownward(e_t.at(i), 2, d_t);
        m->getDownward(e_t.at(i), 0, d_v);
        int v_id = findIn(d_v, 4, vert);
        assert(v_id > -1);

        apf::MeshEntity* t_curr = d_t[lookup_v_tri[v_id]];

        apf::Downward d_vv;
        m->getDownward(t_curr, 0, d_vv);
        int v_id1 = findIn(d_vv, 4, vert);
        assert(v_id1 == -1);

        sn = vd_area_out_n(m, t_curr);

        t_pos = getLinearCentroid(m, t_curr);
        t_pos = t_pos-v_pos;

        assert(t_pos.getLength() > std::numeric_limits<double>::min());
        sn = norm_0(sn);

        if(sn*t_pos < -std::numeric_limits<double>::min())
          sn = sn*(-1);

        double dist = t_pos*sn;
        double v = dn*sn;
        if(v > std::numeric_limits<double>::min()) {
          double t = dist/v;
          t = fabs(t);
          //std::cout << "dist " << dist
          //          << " v " << v
          //          << " t_curr " << t << std::endl;

          //apf::MeshElement* me = createMeshElement(m, e_t.at(i));
          //double vol = measure(me);
          //destroyMeshElement(me);

          //std::cout << "tet vol " << vol;
          //me = createMeshElement(m, t_curr);
          //vol = measure(me);
          //destroyMeshElement(me);

          //std::cout << " tri surf " << vol << std::endl;
          if (t_min < 0 or t < t_min) {
            t_min = t;
            std::cout << "dist " << dist
                      << " v " << v 
                      << " t " << t << "\n"
                      << " v_pos " << v_pos
                      << " t_pos " << t_pos
                      << " dn " << dn
                      << " sn " << sn
                      << std::endl;

            apf::ModelEntity* vert_mdl = m->toModel (vert);
            int vert_type = m->getModelType (vert_mdl);
            int vert_tag = m->getModelTag (vert_mdl);
            std::cout << "vert: " << vert << " "
                      << vert_type << "c" << vert_tag 
                      << std::endl;
            vert_mdl = m->toModel (t_curr);
            vert_type = m->getModelType (vert_mdl);
            vert_tag = m->getModelTag (vert_mdl);
            std::cout << "tri: " << vert_type << "c" << vert_tag 
                      << std::endl;
          }
        }
      }
      // apf::setVector(disp_field, vert, 0, dn*0);
    }
  }
  m->end(it);

  if(t_min < -std::numeric_limits<double>::min())
    t_min = 0;

  printf("Minimum time is %.9f.\n", t_min);
  return t_min;
}

// Vertex based time step estimation:
// This is actually the correct way of calculating the min time, by using all
// triangles across. So any triangle the vertex moves towards is checked for
// inversion.
double vd_find_min_t2(apf::Mesh2* m, bool ext, double t_set) {
  apf::Field* vel_field = m->findField("velocity_field");

  // printf("Field initialized.\n");
  apf::MeshEntity* vert;
  std::vector<apf::MeshEntity*> e_v(0);
  std::vector<apf::MeshEntity*> e_e(0);
  std::vector<apf::MeshEntity*> e_s(0);
  std::vector<apf::MeshEntity*> e_t(0);

  apf::Vector3 t_pos(0,0,0);
  apf::Vector3 v_pos(0,0,0);
  apf::Vector3 dn(0,0,0);
  apf::Vector3 sn(0,0,0);

  double t_min = t_set;
  apf::MeshIterator* it = m->begin(0);
  while ((vert = m->iterate(it))) {
    if(m->getModelType(m->toModel(vert)) != 3) {

      if(ext or !(vd_chk_vert_dom(m, vert))) {
        apf::getVector(vel_field, vert, 0, dn);

        m->getPoint(vert, 0, v_pos);

        init_ent_set(&e_v, vert);
        vd_set_up(m, &e_v, &e_e);
        vd_set_up(m, &e_e, &e_s);
        vd_set_up(m, &e_s, &e_t);

        // The tetrahedra vertex is moving into.
        int tet_id = 0;
        double comp = 0.;

        //std::cout << "v_pos " << v_pos << " dn " << dn << std::endl;
        for(int i = 0; i < e_t.size(); i++) {
          apf::Downward d_t;
          apf::Downward d_v;
          m->getDownward(e_t.at(i), 2, d_t);
          m->getDownward(e_t.at(i), 0, d_v);
          int v_id = findIn(d_v, 4, vert);
          assert(v_id > -1);

          apf::MeshEntity* t_curr = d_t[lookup_v_tri[v_id]];

          apf::Downward d_vv;
          m->getDownward(t_curr, 0, d_vv);
          int v_id1 = findIn(d_vv, 4, vert);
          assert(v_id1 == -1);

          sn = vd_area_out_n(m, t_curr);

          t_pos = getLinearCentroid(m, t_curr);
          t_pos = t_pos-v_pos;

          assert(t_pos.getLength() > std::numeric_limits<double>::min());
          sn = norm_0(sn);

          if(sn*t_pos < -std::numeric_limits<double>::min())
            sn = sn*(-1);

          double dist = t_pos*sn;
          double v = dn*sn;
          if(v > std::numeric_limits<double>::min()) {
            double t = dist/v;
            t = fabs(t);
            //std::cout << "dist " << dist
            //          << " v " << v
            //          << " t_curr " << t << std::endl;

            //apf::MeshElement* me = createMeshElement(m, e_t.at(i));
            //double vol = measure(me);
            //destroyMeshElement(me);

            //std::cout << "tet vol " << vol;
            //me = createMeshElement(m, t_curr);
            //vol = measure(me);
            //destroyMeshElement(me);

            //std::cout << " tri surf " << vol << std::endl;
            if (t_min < 0 or t < t_min) {
              t_min = t;
              std::cout << "dist " << dist
                        << " v " << v 
                        << " t " << t << "\n"
                        << " v_pos " << v_pos
                        << " t_pos " << t_pos
                        << " dn " << dn
                        << " sn " << sn
                        << std::endl;

              apf::ModelEntity* vert_mdl = m->toModel (vert);
              int vert_type = m->getModelType (vert_mdl);
              int vert_tag = m->getModelTag (vert_mdl);
              std::cout << "vert: " << vert << " " 
                        << vert_type << "c" << vert_tag 
                        << std::endl;
              vert_mdl = m->toModel (t_curr);
              vert_type = m->getModelType (vert_mdl);
              vert_tag = m->getModelTag (vert_mdl);
              std::cout << "tri: " << vert_type << "c" << vert_tag 
                        << std::endl;
            }
          }
        }
        //std::cout << "\tt_min " << t_min << std::endl;

      }
      // apf::setVector(disp_field, vert, 0, dn*0);
    }
  }
  m->end(it);

  if(t_min < -std::numeric_limits<double>::min())
    t_min = 0;

  printf("Minimum time is %.9f.\n", t_min);
  return t_min;
}

// This is actually the correct way of calculating the min time, by using all
// triangles across. So any triangle the vertex moves towards is checked for
// inversion.
double vd_find_min_t2(apf::Mesh2* m, std::vector<apf::MeshEntity*>* vert, bool ext, double t_set) {
  apf::Field* vel_field = m->findField("velocity_field");

  // printf("Field initialized.\n");
  //std::vector<apf::MeshEntity*> e_v;
  std::vector<apf::MeshEntity*> e_e(0);
  std::vector<apf::MeshEntity*> e_s(0);
  std::vector<apf::MeshEntity*> e_t(0);

  apf::Vector3 t_pos(0,0,0);
  apf::Vector3 v_pos(0,0,0);
  apf::Vector3 dn(0,0,0);
  apf::Vector3 sn(0,0,0);

  double t_min = t_set;
  for(int j = 0; j < vert->size(); j++) {
    if(ext or !(vd_chk_vert_dom(m, vert->at(j)))) {
      apf::getVector(vel_field, vert->at(j), 0, dn);

      m->getPoint(vert->at(j), 0, v_pos);

      //init_ent_set(&e_v, vert->at(j));
      vd_set_up(m, vert->at(j), &e_e);
      vd_set_up(m, &e_e, &e_s);
      vd_set_up(m, &e_s, &e_t);

      // The tetrahedra vertex is moving into.
      int tet_id = 0;
      double comp = 0.;

      //std::cout << "v_pos " << v_pos << " dn " << dn << std::endl;
      for(int i = 0; i < e_t.size(); i++) {
        apf::Downward d_t;
        apf::Downward d_v;
        m->getDownward(e_t.at(i), 2, d_t);
        m->getDownward(e_t.at(i), 0, d_v);
        int v_id = findIn(d_v, 4, vert->at(j));
        assert(v_id > -1);

        apf::MeshEntity* t_curr = d_t[lookup_v_tri[v_id]];
        sn = vd_area_out_n(m, t_curr);

        t_pos = getLinearCentroid(m, t_curr);
        t_pos = t_pos-v_pos;

        assert(t_pos.getLength() > std::numeric_limits<double>::min());
        sn = norm_0(sn);

        if(sn*t_pos < -std::numeric_limits<double>::min())
          sn = sn*(-1);

        double dist = t_pos*sn;
        double v = dn*sn;
        if(v > std::numeric_limits<double>::min()) {
          double t = dist/v;
          t = fabs(t);
          //std::cout << "dist " << dist
          //          << " v " << v
          //          << " t_curr " << t << std::endl;

          //apf::MeshElement* me = createMeshElement(m, e_t.at(i));
          //double vol = measure(me);
          //destroyMeshElement(me);

          //std::cout << "tet vol " << vol;
          //me = createMeshElement(m, t_curr);
          //vol = measure(me);
          //destroyMeshElement(me);

          //std::cout << " tri surf " << vol << std::endl;
          if (t_min < 0 or t < t_min) {
            t_min = t;
            std::cout << "dist " << dist << " v " << v 
                      << " dn " << dn
                      << " t " << t << "\n"
                      << " v_pos " << v_pos
                      << " t_pos " << t_pos
                      << std::endl;

            apf::ModelEntity* vert_mdl = m->toModel (vert->at(j));
            int vert_type = m->getModelType (vert_mdl);
            int vert_tag = m->getModelTag (vert_mdl);
            std::cout << "vert: " << vert << " "
                      << vert_type << "c" << vert_tag 
                      << std::endl;
            vert_mdl = m->toModel (t_curr);
            vert_type = m->getModelType (vert_mdl);
            vert_tag = m->getModelTag (vert_mdl);
            std::cout << "tri: " << vert_type << "c" << vert_tag 
                      << std::endl;
          }
        }
      }
      //std::cout << "\tt_min " << t_min << std::endl;

    }
    // apf::setVector(disp_field, vert, 0, dn*0);
  }

  if(t_min < -std::numeric_limits<double>::min())
    t_min = 0;

  printf("Minimum time is %.9f.\n", t_min);
  return t_min;
}

// This is actually the correct way of calculating the min time, by using all
// triangles across. So any triangle the vertex moves towards is checked for
// inversion.
double vd_find_min_t2(apf::Mesh2* m, apf::MeshEntity* vert, bool ext, double t_set) {
  apf::Field* vel_field = m->findField("velocity_field");

  // printf("Field initialized.\n");
  std::vector<apf::MeshEntity*> e_v(0);
  std::vector<apf::MeshEntity*> e_e(0);
  std::vector<apf::MeshEntity*> e_s(0);
  std::vector<apf::MeshEntity*> e_t(0);

  apf::Vector3 t_pos(0,0,0);
  apf::Vector3 v_pos(0,0,0);
  apf::Vector3 dn(0,0,0);
  apf::Vector3 sn(0,0,0);

  double t_min = t_set;
  if(ext or !(vd_chk_vert_dom(m, vert))) {
    apf::getVector(vel_field, vert, 0, dn);

    m->getPoint(vert, 0, v_pos);

    init_ent_set(&e_v, vert);
    vd_set_up(m, &e_v, &e_e);
    vd_set_up(m, &e_e, &e_s);
    vd_set_up(m, &e_s, &e_t);

    // The tetrahedra vertex is moving into.
    int tet_id = 0;
    double comp = 0.;

    //std::cout << "v_pos " << v_pos << " dn " << dn << std::endl;
    for(int i = 0; i < e_t.size(); i++) {
      apf::Downward d_t;
      apf::Downward d_v;
      m->getDownward(e_t.at(i), 2, d_t);
      m->getDownward(e_t.at(i), 0, d_v);
      int v_id = findIn(d_v, 4, vert);
      assert(v_id > -1);

      apf::MeshEntity* t_curr = d_t[lookup_v_tri[v_id]];
      sn = vd_area_out_n(m, t_curr);

      t_pos = getLinearCentroid(m, t_curr);
      t_pos = t_pos-v_pos;

      sn = norm_0(sn);
      if(sn*t_pos < -std::numeric_limits<double>::min())
        sn = sn*(-1);
      double dist = t_pos*sn;
      double v = dn*sn;
      if(v > std::numeric_limits<double>::min()) {
        double t = dist/v;
        t = fabs(t);
        //std::cout << "dist " << dist
        //          << " v " << v
        //          << " t_curr " << t << std::endl;

        //apf::MeshElement* me = createMeshElement(m, e_t.at(i));
        //double vol = measure(me);
        //destroyMeshElement(me);

        //std::cout << "tet vol " << vol;
        //me = createMeshElement(m, t_curr);
        //vol = measure(me);
        //destroyMeshElement(me);

        //std::cout << " tri surf " << vol << std::endl;
        t_min = t;
        std::cout << "dist " << dist << " v " << v 
                  << " t " << t << "\n"
                  << " v_pos " << v_pos
                  << " t_pos " << t_pos
                  << " dn " << dn
                  << " sn " << sn
                  << std::endl;

        apf::ModelEntity* vert_mdl = m->toModel (vert);
        int vert_type = m->getModelType (vert_mdl);
        int vert_tag = m->getModelTag (vert_mdl);
        std::cout << "vert: " << vert << " "
                  << vert_type << "c" << vert_tag 
                  << std::endl;
        vert_mdl = m->toModel (t_curr);
        vert_type = m->getModelType (vert_mdl);
        vert_tag = m->getModelTag (vert_mdl);
        std::cout << "tri: " << vert_type << "c" << vert_tag 
                  << std::endl;

      }
    }
    //std::cout << "\tt_min " << t_min << std::endl;

  }
  // apf::setVector(disp_field, vert, 0, dn*0);

  if(t_min < -std::numeric_limits<double>::min())
    t_min = 0;

  printf("Minimum time is %.9f.\n", t_min);
  return t_min;
}

// Calculate the highest multiplier for vector field that would invert any 
// tetrahedra, that is not marked as to be skipped.
double vd_find_min_mult(apf::Mesh2* m, std::vector<apf::MeshEntity*>* vert, 
                        std::map<apf::MeshEntity*, bool> &tet_skip,
                        char const* vec_name) {
  apf::Field* vec_field = m->findField(vec_name);
  assert(vec_field != NULL);

  // printf("Field initialized.\n");
  //std::vector<apf::MeshEntity*> e_v;
  std::vector<apf::MeshEntity*> e_e(0);
  std::vector<apf::MeshEntity*> e_s(0);
  std::vector<apf::MeshEntity*> e_t(0);

  apf::Vector3 t_pos(0,0,0);
  apf::Vector3 v_pos(0,0,0);
  apf::Vector3 dn(0,0,0);
  apf::Vector3 sn(0,0,0);

  double t_min = -1;
  for(int j = 0; j < vert->size(); j++) {
    //if(!(vd_chk_vert_dom(m, vert->at(j)))) {
      apf::getVector(vec_field, vert->at(j), 0, dn);

      m->getPoint(vert->at(j), 0, v_pos);

      //init_ent_set(&e_v, vert->at(j));
      vd_set_up(m, vert->at(j), &e_e);
      vd_set_up(m, &e_e, &e_s);
      vd_set_up(m, &e_s, &e_t);

      // The tetrahedra vertex is moving into.
      int tet_id = 0;
      double comp = 0.;

      //std::cout << "v_pos " << v_pos << " dn " << dn << std::endl;
      for(int i = 0; i < e_t.size(); i++) {
        if(!tet_skip[e_t.at(i)]) {
          apf::Downward d_t;
          apf::Downward d_v;
          m->getDownward(e_t.at(i), 2, d_t);
          m->getDownward(e_t.at(i), 0, d_v);
          int v_id = findIn(d_v, 4, vert->at(j));
          assert(v_id > -1);

          apf::MeshEntity* t_curr = d_t[lookup_v_tri[v_id]];
          sn = vd_area_out(m, t_curr, 0);

          t_pos = getLinearCentroid(m, t_curr);
          t_pos = t_pos-v_pos;

          assert(t_pos.getLength() > std::numeric_limits<double>::min());
          sn = norm_0(sn);

          if(sn*t_pos < - std::numeric_limits<double>::min()
              and std::fabs(sn*t_pos) > std::numeric_limits<double>::min())
            sn = sn*(-1);

          double dist = t_pos*sn;
          double v = dn*sn;
          if(v > std::numeric_limits<double>::min() 
              and std::fabs(v) > std::numeric_limits<double>::min()) {
            double t = dist/v;
            t = fabs(t);
            //std::cout << "dist " << dist
            //          << " v " << v
            //          << " t_curr " << t << std::endl;

            //apf::MeshElement* me = createMeshElement(m, e_t.at(i));
            //double vol = measure(me);
            //destroyMeshElement(me);

            //std::cout << "tet vol " << vol;
            //me = createMeshElement(m, t_curr);
            //vol = measure(me);
            //destroyMeshElement(me);

            //std::cout << " tri surf " << vol << std::endl;
            if (t_min < 0 or t < t_min) {
              t_min = t;
              std::cout << "dist " << dist << " v " << v 
                        << " t " << t << "\n"
                        << " v_pos " << v_pos
                        << " t_pos " << t_pos
                        << std::endl;

              apf::ModelEntity* vert_mdl = m->toModel (vert->at(j));
              int vert_type = m->getModelType (vert_mdl);
              int vert_tag = m->getModelTag (vert_mdl);
              std::cout << "vert: " << vert_type << "c" << vert_tag 
                        << std::endl;
              vert_mdl = m->toModel (t_curr);
              vert_type = m->getModelType (vert_mdl);
              vert_tag = m->getModelTag (vert_mdl);
              std::cout << "tri: " << vert_type << "c" << vert_tag 
                        << std::endl;
            }
          }
        }
      }
      //std::cout << "\tt_min " << t_min << std::endl;

    //}
    // apf::setVector(disp_field, vert, 0, dn*0);
  }
  if(t_min < -std::numeric_limits<double>::min())
    t_min = 0;

  printf("Minimum multiplier is %.9f.\n", t_min);
  return t_min;
}

// Calculate the highest multiplier for vector field that would invert any 
// tetrahedra, that is not marked as to be skipped.
double vd_find_min_mult(apf::Mesh2* m, std::vector<apf::MeshEntity*>* vert, 
                        std::map<apf::MeshEntity*, bool> &tet_skip,
                        std::vector<apf::Vector3> & vect) {

  // printf("Field initialized.\n");
  //std::vector<apf::MeshEntity*> e_v;
  std::vector<apf::MeshEntity*> e_e(0);
  std::vector<apf::MeshEntity*> e_s(0);
  std::vector<apf::MeshEntity*> e_t(0);

  apf::Vector3 t_pos(0,0,0);
  apf::Vector3 v_pos(0,0,0);
  apf::Vector3 sn(0,0,0);

  double t_min = -1;
  for(int j = 0; j < vert->size(); j++) {
    m->getPoint(vert->at(j), 0, v_pos);

    //init_ent_set(&e_v, vert->at(j));
    vd_set_up(m, vert->at(j), &e_e);
    vd_set_up(m, &e_e, &e_s);
    vd_set_up(m, &e_s, &e_t);

    // The tetrahedra vertex is moving into.
    int tet_id = 0;
    double comp = 0.;

    for(int i = 0; i < e_t.size(); i++) {
      if(!tet_skip[e_t.at(i)]) {
        apf::Downward d_t;
        apf::Downward d_v;
        m->getDownward(e_t.at(i), 2, d_t);
        m->getDownward(e_t.at(i), 0, d_v);
        int v_id = findIn(d_v, 4, vert->at(j));
        assert(v_id > -1);

        apf::MeshEntity* t_curr = d_t[lookup_v_tri[v_id]];
        sn = vd_area_out(m, t_curr, 0);

        t_pos = getLinearCentroid(m, t_curr);
        t_pos = t_pos-v_pos;

        assert(t_pos.getLength() > std::numeric_limits<double>::min());
        sn = norm_0(sn);

        if(sn*t_pos < - std::numeric_limits<double>::min()
            and std::fabs(sn*t_pos) > std::numeric_limits<double>::min())
          sn = sn*(-1);

        double dist = t_pos*sn;
        double v = vect.at(j)*sn;
        if(v > std::numeric_limits<double>::min() 
            and std::fabs(v) > std::numeric_limits<double>::min()) {
          double t = dist/v;
          t = fabs(t);
          //std::cout << "dist " << dist
          //          << " v " << v
          //          << " t_curr " << t << std::endl;

          //apf::MeshElement* me = createMeshElement(m, e_t.at(i));
          //double vol = measure(me);
          //destroyMeshElement(me);

          //std::cout << "tet vol " << vol;
          //me = createMeshElement(m, t_curr);
          //vol = measure(me);
          //destroyMeshElement(me);

          //std::cout << " tri surf " << vol << std::endl;
          if (t_min < 0 or t < t_min) {
            t_min = t;
            std::cout << "dist " << dist << " v " << v 
                      << " t " << t << "\n"
                      << " v_pos " << v_pos
                      << " t_pos " << t_pos
                      << std::endl;

            apf::ModelEntity* vert_mdl = m->toModel (vert->at(j));
            int vert_type = m->getModelType (vert_mdl);
            int vert_tag = m->getModelTag (vert_mdl);
            std::cout << "vert: " << vert_type << "c" << vert_tag 
                      << std::endl;
            vert_mdl = m->toModel (t_curr);
            vert_type = m->getModelType (vert_mdl);
            vert_tag = m->getModelTag (vert_mdl);
            std::cout << "tri: " << vert_type << "c" << vert_tag 
                      << std::endl;
          }
        }
      }
    }
  
  }
  if(t_min < -std::numeric_limits<double>::min())
    t_min = 0;

  printf("Minimum multiplier is %.9f.\n", t_min);
  return t_min;
}

double vd_find_max_vel(apf::Mesh2* m) {
  apf::Field* vel_field = m->findField("velocity_field");

  // printf("Field initialized.\n");
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;

  double vel_max = 0.;
  while ((vert = m->iterate(it))) {

    apf::Vector3 point(0,0,0);
    apf::Vector3 dn(0,0,0);

    apf::getVector(vel_field, vert, 0, dn);
    //m->getPoint(vert, 0, point);
    double vel = sqrt(dn*dn);
    //std::cout << vert << ", vel: " << vel << std::endl;
    if (vel > vel_max) {
      vel_max = vel;
    }
    // apf::setVector(disp_field, vert, 0, dn*0);
  }
  m->end(it);
  printf("Maximum velocity is %.3f.\n", vel_max);
  return vel_max;
}


// Given a vertex and a set of velocities, find the minimum allowable time.
double vd_find_min_t_v(apf::Mesh2* m, apf::MeshEntity* v, 
                             std::vector<apf::Vector3>* a_pos, 
                             std::vector<apf::Vector3>* area, 
                             std::vector<apf::Vector3>* v3) {

  double t_min = -1;
  for(int k = 0; k < v3->size(); k++) {
    //std::cout << "velocity(" << k << ") = " << v3->at(k) << std::endl;

    for(int i = 0; i < a_pos->size(); i++) {

      apf::Vector3 norm = norm_0(area->at(i));
      double cos_c = norm*v3->at(k);
      double a_curr = norm*a_pos->at(i);
/*
      std::cout << "\trelative pos " << a_pos->at(i) 
                << " area " << norm
                << " v_proj " << cos_c << " dist " << a_curr
                << std::endl;
*/
      if(cos_c > std::numeric_limits<double>::min()) {
        double t_curr = a_curr/cos_c;
        if(t_curr < t_min or t_min < 0)
          t_min = t_curr;
      }
    }
  }

  return t_min;
}

// Given an entity and velocity field attached to its vertices, calculate the
// rate of change in length/area/volume (first order approximation).
double vd_calc_roc(apf::Mesh2* m, apf::MeshEntity* ent, field_calc* f_calc) {

  int ent_type = m->getType(ent);
  int d = m->typeDimension[ent_type];
  assert(d > 0 and d < 4);

  apf::Downward d_v;
  m->getDownward(ent, 0, d_v);

  std::vector<apf::Vector3> pos(0, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> vel(0, apf::Vector3(0,0,0));
  vel.resize(d+1);
  pos.resize(d+1);

  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  for(int i = 0; i < d+1; i++) {
    m->getPoint(d_v[i], 0, pos.at(i));
    apf::getVector(vel_field, d_v[i], 0, vel.at(i));
  }

  double roc = 0;

  if(d == 1) {
    apf::Vector3 norm1 = norm_0(pos.at(1) - pos.at(0));
    apf::Vector3 norm2 = vel.at(1) - vel.at(0);
    double len = norm2.getLength();

    norm2 = norm_0(norm2);
    roc = norm1*norm2*len;
  }
  else if(d == 2) {
    apf::Vector3 p01 = (pos.at(1) - pos.at(0));
    apf::Vector3 p02 = (pos.at(2) - pos.at(0));

    double len1 = p01.getLength();
    double len2 = p02.getLength();

    p01 = norm_0(p01);
    p02 = norm_0(p02);

    apf::Vector3 a_norm = vd_cross(p01, p02);
    double area = a_norm.getLength();
    a_norm = norm_0(a_norm);

    apf::Vector3 v01 = (vel.at(1) - vel.at(0));
    apf::Vector3 v02 = (vel.at(2) - vel.at(0));

    double s1 = v01.getLength();
    double s2 = v02.getLength();

    v01 = norm_0(v01);
    v02 = norm_0(v02);

    roc = ((vd_cross(v01, p02)*len2*s1 + 
                 vd_cross(p01, v02)*len1*s2)*a_norm)/2;
  }
  else {
    apf::Vector3 p01 = (pos.at(1) - pos.at(0));
    apf::Vector3 p02 = (pos.at(2) - pos.at(0));
    apf::Vector3 p03 = (pos.at(3) - pos.at(0));

    double len1 = p01.getLength();
    double len2 = p02.getLength();
    double len3 = p03.getLength();

    p01 = norm_0(p01);
    p02 = norm_0(p02);
    p03 = norm_0(p03);

    //apf::Vector3 a_norm = cross(p01, p02);
    //double area = a_norm.getLength();
    //a_norm = norm_0(a_norm);

    apf::Vector3 v01 = (vel.at(1) - vel.at(0));
    apf::Vector3 v02 = (vel.at(2) - vel.at(0));
    apf::Vector3 v03 = (vel.at(3) - vel.at(0));

    double s1 = v01.getLength();
    double s2 = v02.getLength();
    double s3 = v03.getLength();

    v01 = norm_0(v01);
    v02 = norm_0(v02);
    v03 = norm_0(v03);

    assert(vd_cross(p01, p02)*p03 > 0);
    roc = ((vd_cross(v01, p02)*len2*s1 
                  + vd_cross(p01, v02)*len1*s2)*p03*len3
                  + vd_cross(p01, p02)*v03*len1*len2*s3)/6;
  }
  assert(!std::isnan(roc));
  return roc;
}

// Given an entity set and velocity field attached to its vertices, calculate 
// the total rate of change in length/area/volume of the set of entities.
double vd_calc_roc(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ent, 
                                                        field_calc* f_calc) {

  assert(ent->size() > 0);
  int ent_type = m->getType(ent->at(0));
  int d = m->typeDimension[ent_type];

  assert(d > 0 and d < 4);

  apf::Downward d_v;

  std::vector<apf::Vector3> pos(d+1, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> vel(d+1, apf::Vector3(0,0,0));

  double roc = 0;

  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);
/*
  // Find the length/area/volume averaged velocity.
  apf::Vector3 mean_v(0,0,0);
  double vol_tot = 0;
  for(int i = 0; i < ent->size(); i++) {
    m->getDownward(ent->at(i), 0, d_v);
    double vol = 0;
    if(d == 3)
      vol = vd_volume_tet(m, d_v);
    else {
      apf::MeshElement* ee = createMeshElement(m, ent->at(i));
      vol = measure(ee);
      destroyMeshElement(ee);
    }

    if(vol > std::numeric_limits<double>::min()) {
      for(int j = 0; j < d+1; j++) {
        apf::getVector(vel_field, d_v[j], 0, vel.at(j));
        mean_v = mean_v + vel.at(j)*vol/(d+1);
      }
      vol_tot = vol_tot + vol;
    }
  }
  if(vol_tot > std::numeric_limits<double>::min()) {
    mean_v = mean_v / vol_tot;
  }
  else
    return -1;
*/
  for(int i = 0; i < ent->size(); i++) {
    m->getDownward(ent->at(i), 0, d_v);

    for(int j = 0; j < d+1; j++) {
      m->getPoint(d_v[j], 0, pos.at(j));
      apf::getVector(vel_field, d_v[j], 0, vel.at(j));
      //vel.at(j) = vel.at(j) - mean_v;
    }

    if(d == 1) {
      apf::Vector3 norm1 = norm_0(pos.at(1) - pos.at(0));
      apf::Vector3 norm2 = vel.at(1) - vel.at(0);
      double len = norm2.getLength();

      norm2 = norm_0(norm2);
      roc = roc + norm1*norm2*len;
    }
    else if(d == 2) {
      apf::Vector3 p01 = (pos.at(1) - pos.at(0));
      apf::Vector3 p02 = (pos.at(2) - pos.at(0));

      double len1 = p01.getLength();
      double len2 = p02.getLength();

      p01 = norm_0(p01);
      p02 = norm_0(p02);

      apf::Vector3 a_norm = vd_cross(p01, p02);
      double area = a_norm.getLength();
      a_norm = norm_0(a_norm);

      apf::Vector3 v01 = (vel.at(1) - vel.at(0));
      apf::Vector3 v02 = (vel.at(2) - vel.at(0));

      double s1 = v01.getLength();
      double s2 = v02.getLength();

      v01 = norm_0(v01);
      v02 = norm_0(v02);

      //roc = roc + (cross(v01, p02)*len2*s1 + cross(p01, v02)*len1*s2)*a_norm;
      roc = roc + ((vd_cross(v01, p02)*len2*s1 + 
                 vd_cross(p01, v02)*len1*s2)*a_norm)/2;
    }
    else {
      apf::Vector3 p01 = (pos.at(1) - pos.at(0));
      apf::Vector3 p02 = (pos.at(2) - pos.at(0));
      apf::Vector3 p03 = (pos.at(3) - pos.at(0));

      double len1 = p01.getLength();
      double len2 = p02.getLength();
      double len3 = p03.getLength();

      p01 = norm_0(p01);
      p02 = norm_0(p02);
      p03 = norm_0(p03);

      //apf::Vector3 anorm = cross(p01, p02);
      //double area = anorm.getLength();
      //a_norm = norm_0(a_norm);

      apf::Vector3 v01 = (vel.at(1) - vel.at(0));
      apf::Vector3 v02 = (vel.at(2) - vel.at(0));
      apf::Vector3 v03 = (vel.at(3) - vel.at(0));

      double s1 = v01.getLength();
      double s2 = v02.getLength();
      double s3 = v03.getLength();

      v01 = norm_0(v01);
      v02 = norm_0(v02);
      v03 = norm_0(v03);

      assert(vd_cross(p01, p02)*p03 > 0);
      //roc =  roc + (cross(v01, p02)*len2*s1 +cross(p01,v02)*len1*s2)*p03*len3
      //                        + cross(p01, p02)*v03*len1*len2*s3;
      roc =  roc + ((vd_cross(v01, p02)*len2*s1 
                    + vd_cross(p01,v02)*len1*s2)*p03*len3
                    + vd_cross(p01, p02)*v03*len1*len2*s3)/6;
    }
  }
  assert(!std::isnan(roc));
  return roc;
}

// Given an entity set and velocity field attached to its vertices, calculate 
// the total rate of change in length/area/volume of the set of entities.
double vd_calc_roc_vol_high(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ent, 
                                                        field_calc* f_calc, 
                                                        double dt) {

  assert(ent->size() > 0);
  int ent_type = m->getType(ent->at(0));
  int d = m->typeDimension[ent_type];

  if(d < 3)
    return 0;

  apf::Downward d_v;

  std::vector<apf::Vector3> pos(d+1, apf::Vector3(0,0,0));
  std::vector<apf::Vector3> vel(d+1, apf::Vector3(0,0,0));

  double roc = 0;

  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  double dt2 = dt*dt;
  for(int i = 0; i < ent->size(); i++) {
    m->getDownward(ent->at(i), 0, d_v);

    for(int j = 0; j < d+1; j++) {
      m->getPoint(d_v[j], 0, pos.at(j));
      apf::getVector(vel_field, d_v[j], 0, vel.at(j));
      //vel.at(j) = vel.at(j) - mean_v;
    }

    apf::Vector3 p01 = (pos.at(1) - pos.at(0));
    apf::Vector3 p02 = (pos.at(2) - pos.at(0));
    apf::Vector3 p03 = (pos.at(3) - pos.at(0));

    double len1 = p01.getLength();
    double len2 = p02.getLength();
    double len3 = p03.getLength();

    p01 = norm_0(p01);
    p02 = norm_0(p02);
    p03 = norm_0(p03);

    //apf::Vector3 anorm = cross(p01, p02);
    //double area = anorm.getLength();
    //a_norm = norm_0(a_norm);

    apf::Vector3 v01 = (vel.at(1) - vel.at(0));
    apf::Vector3 v02 = (vel.at(2) - vel.at(0));
    apf::Vector3 v03 = (vel.at(3) - vel.at(0));

    double s1 = v01.getLength();
    double s2 = v02.getLength();
    double s3 = v03.getLength();

    v01 = norm_0(v01);
    v02 = norm_0(v02);
    v03 = norm_0(v03);

    assert(vd_cross(p01, p02)*p03 > 0);
    //roc =  roc + (cross(v01, p02)*len2*s1 +cross(p01,v02)*len1*s2)*p03*len3
    //                        + cross(p01, p02)*v03*len1*len2*s3;
    double dVdt = ((vd_cross(v01, p02)*len2*s1 
                  + vd_cross(p01,v02)*len1*s2)*p03*len3
                  + vd_cross(p01, p02)*v03*len1*len2*s3);
    double dVdt2 = ((vd_cross(v01, p02)*len2*s1 
                  + vd_cross(p01,v02)*len1*s2)*v03*s3
                  + vd_cross(v01, v02)*p03*s1*s2*len3)*2;
    double dVdt3 = vd_cross(v01,v02)*s1*s2*v03*s3*6;
    roc =  roc + (dVdt + dVdt2*dt + dVdt3*dt2)/6;
  }
  assert(!std::isnan(roc));
  return roc;
}


// Add a zero displacement field at a vertex.
void vd_add_0_disp_field(apf::Mesh2* m, apf::MeshEntity* vert) {
  apf::Field* disp_field = m->findField("displacement_field");
  assert(disp_field);
  apf::Vector3 dn(0,0,0);
  apf::setVector(disp_field, vert, 0, dn);
}

// Add a zero velocity field at a vertex.
void vd_add_0_vel_field(apf::Mesh2* m, apf::MeshEntity* vert) {
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);
  apf::Vector3 dn(0,0,0);
  apf::setVector(vel_field, vert, 0, dn);
}

apf::Vector3 vd_get_v_vel_field(apf::Mesh2* m, apf::MeshEntity* vert) {
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);
  apf::Vector3 dn(0,0,0);
  apf::getVector(vel_field, vert, 0, dn);
  return dn;
}

double vd_get_max_vel_norm(apf::Mesh2* m) {
  apf::Field* vel_field = m->findField("velocity_field");
  assert(vel_field);

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;
  apf::Vector3 dn(0,0,0);
  double norm_max = 0.;
  while ((vert = m->iterate(it))) {
    apf::getVector(vel_field, vert, 0, dn);
    double norm_curr = dn.getLength();
    if(norm_curr > norm_max)
      norm_max = norm_curr;
  }
  m->end(it);

  return norm_max;
}

bool isvecnan(apf::Vector3 const& v_in) {
  if(std::isnan(v_in[0]) or std::isnan(v_in[1]) or std::isnan(v_in[2])) {
    itisnan();
    return true;
  }
  return false;
}

void itisnan() {
}

void vd_tet_vel(apf::Mesh2* m) {

  apf::MeshEntity* elem;
  apf::MeshElement* ee;

  apf::MeshIterator* it = m->begin(3);

  apf::Field* vel_field = m->findField("velocity_field");

  double dt_min = 10;
  while ((elem = m->iterate(it))) {

    ee = createMeshElement(m, elem);
    double vol = measure(ee);
    std::cout << elem << ", vol: " << vol << std::endl;
    destroyMeshElement(ee);

    apf::Downward d_v;
    apf::Downward d_s;

    m->getDownward(elem, 0, d_v);
    m->getDownward(elem, 2, d_s);
    for(int i = 0; i < 4; i++) {
      ee = createMeshElement(m, d_s[i]);
      double area = measure(ee);
      std::cout << d_s[i] << ", area: " << area << std::endl;
      destroyMeshElement(ee);

      apf::Vector3 point(0,0,0);
      apf::Vector3 dn(0,0,0);
      apf::getVector(vel_field, d_v[lookup_v_tri[i]], 0, dn);
      m->getPoint(d_v[lookup_v_tri[i]], 0, point);
      double vel = sqrt(dn*dn);
      std::cout << "\t" << d_v[lookup_v_tri[i]] << "pos: " << point 
                << ", vel: " << dn << ", mag: " << vel << std::endl;

      double dt;
      if(vel > std::numeric_limits<double>::min() and
         area > std::numeric_limits<double>::min())
      dt = 3*vol/area/vel;

      if(dt < dt_min)
        dt_min = dt;

      std::cout << "\ttime " << dt << std::endl;
    }

    // apf::setVector(disp_field, vert, 0, dn*0);
  }
  m->end(it);
  std::cout << "dt_min " << dt_min << std::endl;

}

