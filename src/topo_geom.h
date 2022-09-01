#ifndef TOPO_GEOM_H
#define TOPO_GEOM_H

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include "ma.h"

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_extinfo.h"
#include "topo_topo.h"
#include "topo_entlist.h"

#define PI_L std::acos(-1.)  /* pi */

struct vd_plane {
  apf::Vector3 pos;
  apf::Vector3 norm;
};

struct vd_inter {
  int cut;
  apf::Vector3 pos;
};

apf::Vector3 norm_0(apf::Vector3 v_in);

// Given an entity, get its vertices and return the center point.
//apf::Vector3 vd_get_center(apf::Mesh2* m, apf::MeshEntity* ent);

// Given an entity set of cube corner vertices, return the center point.
apf::Vector3 vd_get_center(apf::Mesh2* m, std::vector<apf::MeshEntity*>* set_in);
apf::Vector3 vd_get_center_e(apf::Mesh2* m, std::vector<apf::MeshEntity*>* set_in);

// Given two edges joined at v_ctr, find the plane defined by the edges.
apf::Vector3 vd_get_e_plane(apf::Mesh2* m, apf::MeshEntity* v_ctr, apf::MeshEntity* e1, apf::MeshEntity* e2);

// Given an entity, return the length, area or volume.
double vd_meas_ent(apf::Mesh2* m, apf::MeshEntity* e);

// Given the dimension, return the total length/area/volume of the entities 
// belonging to strata of the same dimension.
double vd_meas_ent_dim(apf::Mesh2* m, int dim);

// Given a set of vertices, return the radius from the center of mass.
double vd_meas_rad_set(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev);
// Given a set of vertices, return the radius from the center of mass.
double vd_meas_rad_set_min(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev);

// Given a set of entities, return the minimum radius of the containing sphere 
// positioned at the center. Get the position of minimum distance 
// from the center. Get the second position, that is twice that distance from 
// the first position and is minimum among these.
std::pair<double, double> vd_meas_rad_set_minmax(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ev);

// Given a set of entities, return the total length, area or volume.
double vd_meas_set(apf::Mesh2* m, std::vector<apf::MeshEntity*>* es_in);


double vd_avg_meas(apf::Mesh2* m, std::vector<apf::MeshEntity*>* e_set);

// Given an tet set, return the min(max) displacement to triangles across the
// center vertex.
double vd_minmax_meas(apf::Mesh2* m, std::vector<apf::MeshEntity*>* e_tet, apf::MeshEntity* v, bool min = true);

// Given an entity set, return the min(max) length, area, volume of the 
// entities.
double vd_minmax_meas(apf::Mesh2* m, std::vector<apf::MeshEntity*>* e_set, bool min = true);

// Given a set of triangles joined at a vertex, find the 
// gradient of the total area, wrt. the vertex position.
apf::Vector3 vd_grad(apf::Mesh2* m, std::vector<apf::MeshEntity*>* tri, apf::MeshEntity* vert);

// Given a distance field and a tet on which the field is defined, 
// approximate the distance defined on one the bounding vertices, approximating
// the gradient of the distance field by using the distances on the other 
// vertices and assuming a linear variation of the distance across the tet. 
double approx_dist(apf::Mesh2* m, std::map<apf::MeshEntity*, double> & dist, apf::MeshEntity* tet, apf::MeshEntity* vert);
double approx_dist(apf::Mesh2* m, std::map<apf::MeshEntity*, double> & dist, apf::Downward d_v, apf::MeshEntity* vert);

//-----------------------------------

// For information retrival, geometric modifications, etc. it might be wise to:
// Define an object, Grain, that keeps geometry tag, average orientation, 
// functions to get elements, do calculations, etc. Assuming geometry will be 
// modified properly, it is not necessary to keep element list. Surfaces can 
// be obtained from gmi, which can be used to obtain grain adjacencies.

// print the adjacencies of model entities of a given dimension. Can be 
// used as a template.
void vd_print_adj(apf::Mesh2* m, int dim);

// Given an entity, get its vertices and return the center point.
apf::Vector3 vd_get_center(apf::Mesh2* m, apf::MeshEntity* ent);

// Given an entity set of cube corner vertices, return the center point.
apf::Vector3 vd_get_center(apf::Mesh2* m, Entity_set* set_in);

// Adapted from getAverageEdgeLength in maSize.
// Get the average entity length, area or volume.
double getMedianEntSize_hist(apf::Mesh2* m, int dim, int N = 100);
double getMedianEntSize_hist(std::vector<double> &sz, int N = 100);

double getAverageEntSize(apf::Mesh2* m, int dim);
double getAverageEntSize(apf::Mesh2* m, apf::MeshEntity* v, int dim);

// Given the geometry dimension and tag, return the size of the geometry:
// total length, area or volume.
double vd_meas_geom(apf::Mesh2* m, int em_dim, int em_tag);
// Average entity size in geometry.
double vd_meas_geom_avg(apf::Mesh2* m, int em_dim, int em_tag);

// Given a set of entities, return the total length, area or volume.
double vd_meas_set(apf::Mesh2* m, Entity_set* es_in);

// Given a vertex, return the number of negative volume tet adjacencies and 
// print out the body center and cell membership of the tets.
long int vd_chk_neg_vert(apf::Mesh2* m, apf::MeshEntity* vert);

// TODO silly implementation. Given a tet, return 1 if negative volume.
long int vd_chk_neg(apf::Mesh2* m, apf::MeshEntity* elem);

// Given a vector of entities, return the number of negative volume elements.
long int vd_chk_neg(apf::Mesh2* m, std::vector<apf::MeshEntity*>* elem);

// Given a vector of entities, return the number of negative volume elements.
long int vd_chk_neg(apf::Mesh2* m, Entity_set* elem);

// Going over the mesh elements, return the number of negative volume elements,
// using custom function that performs better with small edge lengths.
long int vd_chk_neg_sgn(apf::Mesh2* m);

// Going over the mesh elements, return the number of negative volume elements.
long int vd_chk_neg(apf::Mesh2* m);

// Given and edge and a direction, check if splitting the edge and moving the  
// vertex to the position inverts any tetrahedra bounded by the edge.
bool vd_edge_inv(apf::Mesh2* m, apf::MeshEntity* edge, apf::Vector3 pos);

// Print the normals for elements inside the mesh.
void vd_print_norm(apf::Mesh2* m);

// Given a point and an edge, find the projection of the point on the edge.
apf::Vector3 vd_projpt(apf::Mesh2* m, apf::MeshEntity* edge, apf::Vector3 pt);
apf::Vector3 vd_projpt(apf::Vector3 e1, apf::Vector3 e2, apf::Vector3 pt);

// Given a volume geometry tag, and a boundary surface triangle, 
// return cross product of area normal with the difference of edge vectors to 
// obtain the gradient of the area w.r.t. the vertex position. 
// It assumes the element entities belong to different geometries. 
// TODO The computational complexity can be reduced by potentially simplifying
// vector calculations.
apf::Vector3 vd_neg_grad(apf::Mesh2* m, apf::MeshEntity* vert, apf::MeshEntity* surf, int geom);

// Given a volume geometry tag, and a boundary surface triangle, 
// return cross product of the edge vectors, pointing outwards from that 
// volume.
apf::Vector3 vd_area_out(apf::Mesh2* m, apf::MeshEntity* surf, int geom);

apf::Vector3 vd_area_out(apf::Mesh2* m, apf::MeshEntity* surf);

apf::Vector3 vd_area_out_n(apf::Mesh2* m, apf::MeshEntity* surf);

apf::Vector3 vd_area_out_n(apf::Mesh2* m, std::vector<apf::Vector3>* pos);
apf::Vector3 vd_area_out_n(apf::Vector3 p0, apf::Vector3 p1, apf::Vector3 p2);

// Given a set of vertices consituting a triangle, calculate the area normal. 
apf::Vector3 vd_area_tri(apf::Mesh2* m, apf::MeshEntity** v);

// Find the normal on plane p0p1p2, perpendicular to p0p1, towards p2. 
apf::Vector3 vd_dir_in_pl(apf::Vector3 p0, apf::Vector3 p1, apf::Vector3 p2);

// Print the topological and geometric information of the tets around an ent.
void print_tets_ent(apf::Mesh2* m, apf::MeshEntity* ent);

// Given a tetrahedron, print the geometric and topological membership 
// information.
void print_tet(apf::Mesh2* m, apf::MeshEntity* tet);

// Given a mesh, calculate the total volume. 
double vd_tot_volm(apf::Mesh2* m);

// Cross product using orthonormalized vectors to mitigate numerical errors.
apf::Vector3 vd_cross(apf::Vector3 v1, apf::Vector3 v2);

// Given three vectors v03, v01, v02 calculate the triple product. 
double vd_trip_prod(apf::Vector3 v03, apf::Vector3 v01, apf::Vector3 v02);

// Given a set of vertices consituting a tetrahedron, calculate the volume. 
double vd_volume_tet(apf::Mesh2* m, apf::MeshEntity* tet);
double vd_volume_tet(apf::Mesh2* m, apf::MeshEntity** v);
double vd_volume_tet(std::vector<apf::Vector3>* points);

// Given a tetrahedron, return true if the volume is positive. 
bool vd_volume_tet_sign(apf::Mesh2* m, apf::MeshEntity* tet);
bool vd_volume_tet_sign(apf::Mesh2* m, apf::MeshEntity** v);

// Print the linear centroid of the given element.
apf::Vector3 vd_print_pos(apf::Mesh2* m, apf::MeshEntity* ent);
apf::Vector3 vd_get_pos(apf::Mesh2* m, apf::MeshEntity* ent);

// Return area normal of the triangle towards interior of the tet. 
apf::Vector3 vd_get_pos(apf::Mesh2* m, apf::MeshEntity* tet, 
                                       apf::MeshEntity* tri);

// Print the position of the other vertex of the edge.
apf::Vector3 vd_edge_pos_other(apf::Mesh2* m, apf::MeshEntity* edge, 
                               apf::MeshEntity* vert);

// Get the position of the vertex across the triangle on the adjacent tetrahedra
// in the es_ent. Assume there is strictly one such vertex.
apf::Vector3 get_pos_v_x_t_in_list(apf::Mesh2* m, apf::MeshEntity* tri_curr, std::vector<apf::MeshEntity*> & es_ent);

// Given a vertex, return the average edge length of the edges bounded by the
// vertex. 
double vd_avg_len(apf::Mesh2* m, apf::MeshEntity* vertex);

// Given an entity set, return the average length, area, volume of the entities.
double vd_avg_meas(apf::Mesh2* m, Entity_set* e_set);

apf::Vector3 get_edge_dir(apf::Mesh2* m, apf::MeshEntity* edge);

// Given an entity set, return the min(max) length, area, volume of the 
// entities.
double vd_minmax_meas(apf::Mesh2* m, Entity_set* e_set, bool min = true);
//double vd_minmax_meas(apf::Mesh2* m, std::vector<apf::MeshEntity*>* e_set, bool min = true);

// Print the normals for elements inside the mesh, using area out.
void vd_print_norm2(apf::Mesh2* m);

// Given two vectors, calculate the inner angle.
double vd_inner_angle(apf::Vector3 v1, apf::Vector3 v2);

// Given an interior direction and two directions, all lying on the same plane,
// find the angle of the arc that the interior direction dissects.
double vd_int_angle_n(apf::Vector3 norm_in, apf::Vector3 norm1, apf::Vector3 norm2, double tol = 10e-1);

// Given two adjacent boundary surfaces and the geometry of the inside volume, 
// find the exterior angle defined between those surfaces.
double vd_ext_angle(apf::Mesh2* m, apf::MeshEntity* surf1, 
                                  apf::MeshEntity* surf2, int geom);

// TODO this doesn't look like a good approach. 
// Given a vertex belonging to a 1cell, find two triangles bounding the 3-cell.
// Return the exterior angle. The triangles found in this fashion do not 
// necessarily share an edge. Also if 3cell has disjoint sets around v, it can
// also find triangles that do not bound the same disjoint set. 
double vd_ext_angle(apf::Mesh2* m, apf::MeshEntity* v, int cell);

// Given the outward normal and a point on each planes bounding the same volume, 
// find the exterior angle that joins. Assume all points are coplanar and 
// the planes intersect at zero point.
double vd_ext_angle_n(apf::Vector3 pos1, apf::Vector3 pos2,
                    apf::Vector3 norm1, apf::Vector3 norm2);

// Given two planes and points on each plane, the direction of the line that
// joins the planes, a center point on one side of the planes, find the exterior
// angle to the angle that lies on the side of the planes with the center point.
double vd_ext_angle(apf::Vector3 pos1, apf::Vector3 pos2, apf::Vector3 l_dir, 
                                apf::Vector3 l_pos, apf::Vector3 ctr);

// Given a vertex on a 2cell, and an interior 3cell, find the mean curvature
// of the surface associated with the normal outward from the 3cell.
double vd_mean_curv(apf::Mesh2* m, apf::MeshEntity* v, int cell);

// Given a vertex and an interior 3cell, find the Gaussian curvature
// of the surfaces of the 3cell joining at the vertex by angular defect formula.
double vd_gauss_curv(apf::Mesh2* m, apf::MeshEntity* v, int cell);

// Given a set of triangles joined at a vertex, find the negative of the 
// gradient of the total area, wrt. the vertex position.
apf::Vector3 vd_neg_grad(apf::Mesh2* m, Entity_set* tri, apf::MeshEntity* vert);

// Given two vertices A and B that form a vector V from the second one to the 
// first one, and two other vertices C and D, check if the area vector obtained  
// by taking the cross products of vectors V_C and V_D from the center point of 
// A and B to the other two vertices is positively aligned with the vector V.
bool vd_vec_align(apf::Mesh2* m, apf::MeshEntity** vert);

// Given a triangle and a vertex, calculate the angle at the vertex corner in 
// radians.
double vd_int_angle(apf::Mesh2* m, apf::MeshEntity* tri, apf::MeshEntity* v);

// Given a two edges bounded by a vertex v, calculate the angle at the vertex  
// corner in radians.
double vd_int_angle(apf::Mesh2* m, apf::MeshEntity* e_1, apf::MeshEntity* e_2, 
                                                            apf::MeshEntity* v);

// Given a set of triangles bounded by a vertex, find the smallest and largest
// distance to the edges across.
std::pair<double, double> vd_find_dist_e(apf::Mesh2* m, 
                                     std::vector<apf::MeshEntity*>* tri, 
                                     apf::MeshEntity* v);
// Used in exterior insertions.
std::pair<double, double> vd_find_dist_e(apf::Mesh2* m, 
                                     std::vector<apf::MeshEntity*>* tri, 
                              apf::MeshEntity* v, apf::Vector3 v_pos);

std::pair<double, double> vd_find_dist_e(apf::Mesh2* m, 
                                  std::vector<apf::MeshEntity*>* edge, 
                                      apf::Vector3 v_pos);

// Given a set of edges bounded by a vertex, find the smallest and largest
// distance to the vertices across.
std::pair<double, double> vd_find_dist_v(apf::Mesh2* m, 
                                     std::vector<apf::MeshEntity*>* edge, 
                                     apf::MeshEntity* v);

std::pair<double, double> vd_find_dist_v(apf::Mesh2* m, 
                                     std::vector<apf::MeshEntity*>* edge, 
                              apf::MeshEntity* v, apf::Vector3 v_pos);

// Skip given tetrahedra. Used mainly in collapse preconditioning, where there
// are inverting elements around vertex. In that case, do not consider the 
// distances to triangles within inverted tetrahedra. If none found, return
// -1.
double vd_dist_v_x(apf::Mesh2* m, apf::MeshEntity* v, apf::Vector3 dir, 
        std::map<apf::MeshEntity*, bool> &skip, apf::Vector3* int_pt = NULL);

// Given a vertex v and a direction, find the distance to the nearest triangle 
// along the direction. Assume direction vector is non-zero and all tetrahedra 
// are positive volume.
double vd_dist_v_x(apf::Mesh2* m, apf::MeshEntity* v, apf::Vector3 dir, 
                                                  apf::Vector3* int_pt= NULL);

// Given a vertex v find the displacement to the nearest edge 
// along the plane perpendicular to the edge.
apf::Vector3 vd_dist_v_ex(apf::Mesh2* m, apf::MeshEntity* v);
// Given a vertex v find the displacement to the nearest triangle 
// along the surface normal. Assume tetrahedra are positive volume.
apf::Vector3 vd_dist_v_x(apf::Mesh2* m, apf::MeshEntity* v);

// Given a vertex v, a set of triangles across and a direction, find the 
// displacement to the nearest triangle plane along the given direction.
double vd_dist_v_x_pl_dir(apf::Mesh2* m, apf::MeshEntity* v, 
                        std::vector<apf::MeshEntity*> &es_x, apf::Vector3 dir);

apf::Vector3 vd_find_dist_w(apf::Mesh2* m, 
                    std::vector<apf::MeshEntity* >* tet, apf::MeshEntity* v);
apf::Vector3 vd_find_dist_w_tri(apf::Mesh2* m, 
                    std::vector<apf::MeshEntity*>* tri, apf::MeshEntity* v);

apf::Vector3 vd_find_dist_w(apf::Mesh2* m, apf::MeshEntity* tet, 
                                                         apf::MeshEntity* v);
double vd_find_dist(apf::Mesh2* m, apf::MeshEntity* tet, apf::MeshEntity* v);
double vd_find_dist_c(apf::Mesh2* m, apf::MeshEntity* tet, apf::MeshEntity* v);

apf::Vector3 get_volm_pres(apf::Mesh2* m, apf::MeshEntity* v, 
                      std::vector<apf::MeshEntity*>* tet, apf::Vector3 v_in);

apf::Vector3 get_volm_pres2(apf::Mesh2* m, apf::MeshEntity* v, 
                      std::vector<apf::MeshEntity*>* tet, apf::Vector3 v_in);

apf::Vector3 rem_norm_comp(apf::Mesh2* m, std::vector<apf::MeshEntity*>* tri, 
                                  apf::Vector3 v_in);

void vd_get_apos(apf::Mesh2* m, apf::MeshEntity* v, 
                         std::vector<apf::MeshEntity*>* tet, 
                         std::vector<apf::Vector3>* area, 
                         std::vector<apf::Vector3>* a_pos);

bool pt_on_tri(apf::Mesh2* m, apf::Vector3 v_pos, apf::MeshEntity* tri);

std::pair<double, double> vd_get_apos(apf::Mesh2* m, apf::Vector3 v_pos, 
                         std::vector<apf::MeshEntity*>* tet, 
                         std::vector<apf::Vector3>* area, 
                         std::vector<apf::Vector3>* a_pos);

// Check whether the given triangles intersect/lie within the sphere. If not, 
// return the normal distance to the closest triangle. OBSOLETE
//std::pair<bool, double> tri_int_sphere(apf::Mesh2* m, apf::Vector3 v_pos, 
//                         double r, std::vector<apf::MeshEntity*>* tri);

// Get the radius of the largest sphere that doesn't intersect any triangles
// across the vertex and the radius of the largest sphere without intersecting
// a single triangle.
std::pair<double, double> vert_dist_sphere(apf::Mesh2* m, apf::MeshEntity* v);

// Given a point and a set of triangles, find the radius of the largest sphere
// centered at the point that at most is tangent to one of the triangle edges 
// or intersects at a corner point.
std::pair<double, double> tri_dist_sphere(apf::Mesh2* m, apf::Vector3 v_pos, 
                      std::vector<apf::MeshEntity*>* tri);

void vd_print_vert(apf::Mesh2* m, apf::MeshEntity* ent);

std::vector<double> upd_cell_rad(apf::Mesh2* m);
std::vector<double> upd_cell_rad(apf::Mesh2* m, vd_entlist* e_list);

// Return 0 if doesn't intersect, 1 if intersects, -n-1, if it intersects the
// nth vertex. If intersects multiple points, return the -n_low-1, where n_low
// is the lowest index.
// Given a cloud of vertices, check if the plane intersects the cloud.
int pl_int_verts(apf::Mesh2* m, std::vector<apf::MeshEntity*>* ent, apf::Vector3 p_pos, apf::Vector3 p_norm, double tol = 0);

// Given two points, check if the plane lies between them.
void pl_int_2pts(apf::Mesh2* m, apf::Vector3 p1, apf::Vector3 p2, 
                                              vd_plane* pl, vd_inter* v_int);
// Given an edge, check if the plane intersects the edge.
int pl_int_edge(apf::Mesh2* m, apf::MeshEntity* ent, apf::Vector3 p_pos, apf::Vector3 p_norm, double tol = 0);

void pl_int_edge(apf::Mesh2* m, apf::MeshEntity* ent, vd_plane* pl, vd_inter* v_int, double tol = 0);

// Given a triangle, check if the plane intersects the triangle.
int pl_int_tri(apf::Mesh2* m, apf::MeshEntity* ent, apf::Vector3 p_pos, apf::Vector3 p_norm);

// Calculate if a triangle intersects a plane. If so, set the direction of the 
// intersection line and return whether the intersection is valid.
// If all lie on one side of the plane no intersection, return false.
// If all lie on the plane, return false.
// If two lie on the plane and one on one side, return true.
// If two lie on one side and one on the line, return false.
// If two lie on one side and one on the other, true.
// If return value is true, set l_dir as cross product of the p_norm and area 
// norm.
bool pl_int_tri_line(apf::Mesh2* m, apf::MeshEntity* ent, 
                                    apf::Vector3 p_pos, apf::Vector3 p_norm, 
                                    apf::Vector3 &l_dir, apf::Vector3 &l_pos, double tol = -1);

// Given a tetrahedron and a plane, find a bounding triangle that intersects the
// plane. Return true if found and set the line of intersection.
bool pl_int_tet_line(apf::Mesh2* m, apf::MeshEntity* tet, 
                                    apf::Vector3 p_pos, apf::Vector3 p_norm, 
                                    apf::Vector3 &l_dir, apf::Vector3 &l_pos);

// Find the intersection of a line and a plane. Assume they intersect.
apf::Vector3 pl_int_line(apf::Vector3 l_pos, apf::Vector3 l, 
                          apf::Vector3 pl_pos, apf::Vector3 pl_norm);

// Return true if line intersects the triangle.
bool tri_int_line(apf::Vector3 v1, apf::Vector3 v2, apf::Vector3 v3, 
                  apf::Vector3 o, apf::Vector3 l, apf::Vector3* int_pt = NULL);

bool tri_int_line(apf::Mesh2* m, apf::MeshEntity* tri,
                  apf::Vector3 o, apf::Vector3 l, apf::Vector3* int_pt = NULL);

// Does a point lie on a triangle. Assume coplanar.
bool pt_int_tri(apf::Mesh2* m, apf::MeshEntity* tri, apf::Vector3 pt, double tol = 0.01);
bool pt_int_tri(apf::Vector3 v1, apf::Vector3 v2, apf::Vector3 v3, 
                  apf::Vector3 pt, double tol = 0.01);

// Return true if line intersects the triangle.
bool tri_int_edge(apf::Vector3 v1, apf::Vector3 v2, apf::Vector3 v3, 
                  apf::Vector3 o, apf::Vector3 l, 
                  double tol = 0.01, apf::Vector3* int_pt = NULL);

// Find the intersection of the line segment and the sphere, centered at c_pos
// with radius r. Assume the line segment starts at o, inside the sphere, and 
// l points from o to a position outside the sphere.
apf::Vector3 sph_int_line(apf::Vector3 c_pos, double r, 
                          apf::Vector3 l, apf::Vector3 o);

// Calculate the incenter of an entity.
// TODO resolve potential bug, some elements get inverted when splitting using 
// the incenter
apf::Vector3 inctr_ent(apf::Mesh2* m, apf::MeshEntity* ent);
apf::Vector3 inctr_ent(std::vector<apf::Vector3>* v);

// Given pos, return pos mirrored by plane defined by ori and dir.
apf::Vector3 mir_pos_pl(apf::Vector3 ori, apf::Vector3 dir_m, apf::Vector3 pos);
// Given pos, return pos mirrored by line defined by ori and dir.
apf::Vector3 mir_pos_ln(apf::Vector3 ori, apf::Vector3 dir_m, apf::Vector3 pos);
// Given dir, return dir mirrored by plane defined dir_m.
apf::Vector3 mir_dir_pl(apf::Vector3 dir_m, apf::Vector3 dir);
// Given dir, return dir mirrored by line defined by dir_m.
apf::Vector3 mir_dir_ln(apf::Vector3 dir_m, apf::Vector3 dir);

// Given pos, return pos rotated around ori by rot.
apf::Vector3 rot_pos(apf::Vector3 ori, apf::Matrix3x3* rot, apf::Vector3 pos);
// Given dir, return dir rotated by rot.
apf::Vector3 rot_dir(apf::Matrix3x3* rot, apf::Vector3 dir);

// Given dir and angle, determine the rotation matrix.
apf::Matrix3x3 rot_matrix_dir(apf::Vector3 dir, double angle);

// Invert a 3x3 matrix. Uses vd_cross instead of cross which is prone to numerical round-off errors.
apf::Matrix3x3 vd_invert(apf::Matrix3x3& mat_in);

// Find the center and radius of circle passing through three points.
std::pair<apf::Vector3, double> three_pt_circ(std::vector<apf::Vector3>& v);

// Find the center and radius of sphere passing through four points.
std::pair<apf::Vector3, double> four_pt_sph(std::vector<apf::Vector3>& v);

#endif
