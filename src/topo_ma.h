#ifndef TOPO_MA_H
#define TOPO_MA_H

#include "apf.h"
#include "ma.h"
#include "maShape.h"
#include "maSize.h"

#include "topo_topo.h"
#include "topo_entlist.h"

std::pair<apf::MeshEntity*, double> calc_valid_q(apf::Mesh2* m, ma::SizeField* s_f);

double get_low_q(apf::Mesh2* m, ma::SizeField* s_f, std::vector<apf::MeshEntity*> &tets, double q_th = 10e-5);
bool vd_chk_vol_valid(apf::Mesh2* m, apf::MeshEntity* tet, apf::MeshEntity* ent);

double calc_good_vol(apf::Mesh2* m, double ref_len, double m_len);
double calc_good_q(apf::Mesh2* m, double ref_len, double m_len);

double calc_q(apf::Mesh2* m, std::vector<apf::Vector3> &v);

// ADAPT_1CELL: Uniform refinement based on average 1cell length. 
// ADAPT_3CELL: Uniform refinement based on average 3cell radius. 
// ADAPT_STEP: Stepwise refinement based on average 3cell radius, boundaries
// refined more than interior.
// ADAPT_STEP_1CELL: Stepwise refinement based on average 3cell radius, 
// boundaries refined more than interior and 1cell has two edge at the minimum.

// ADAPT_SPLIT: Adapt using only edge split. 
// ADAPT_COL: Adapt using only edge coarsening. 
// ADAPT_Q: Adapt using only fix shape. 
// ADAPT_BALANCE: Adapt using PARMA and ZOLTAN parallel load balancers. 
// ADAPT_ALL: Adapt using all flags open.
enum class INPUT_TYPE {
  ADAPT_SPLIT,
  ADAPT_COL,
  ADAPT_Q,
  END
};

// Model based sizefields:
//  EDGE_COL: ModelEdgeCollapse
//  EDGE_SPLIT: ModelEdgeSplit
//  EDGE_REFINER: ModelEdgeRefiner
//  EDGE_REF_DIST: ModelEdgeRefinerDist
enum class MA_SIZE_TYPE {
  EDGE_COL,
  EDGE_SPLIT,
  EDGE_REFINER,
  EDGE_REF_DIST,
  END
};


// MA tester, to check if the current state of the mesh leads to cavityexists 
// crashes. Adapted from maEdgeSwap.cc
#define MAX_VERTS 7
/* this class represents a loop of vertices around
   an edge to be swapped. Most of its code is dedicated
   to finding and storing this loop using adjacency searches
   from the edge and one key face attached to it */
class MaLoop
{
  public:
    void init(apf::Mesh2* m)
    {
      mesh = m;
    }
    void setEdge(apf::MeshEntity* e)
    {
      edge = e;
      mesh->getDownward(e,0,edge_verts);
    }
    /* returns true iff going from vert to the other
       vertex in the tet around the edge curls in the ev[0]-ev[1]
       direction. This establishes an orientation for
       the loop and its triangles with respect to the edge
       so that the generated tets are not left-handed */
    bool isCurlOk(apf::MeshEntity* vert, apf::MeshEntity* tet)
    {
      apf::MeshEntity* tv0[4];
      mesh->getDownward(tet,0,tv0);
      int n = findIn(tv0,4,vert);
      apf::MeshEntity* tv[4];
      ma::rotateTet(tv0,n*3,tv);
      assert(tv[0]==vert);
      int a = findIn(tv+1,3,edge_verts[0]);
      int b = findIn(tv+1,3,edge_verts[1]);
      if (b == ((a+1)%3))
        return true;
      assert(a == ((b+1)%3));
      return false;
    }
    apf::MeshEntity* findInitialTet(apf::MeshEntity* vert, apf::MeshEntity* face)
    {
      apf::Up tets;
      mesh->getUp(face,tets);
      for (int i=0; i < tets.n; ++i)
      {
        apf::MeshEntity* tet = tets.e[i];
        if (isCurlOk(vert,tet))
          return tet;
      }
      return 0;
    }
    void findFromFace(apf::MeshEntity* startFace)
    {
      apf::MeshEntity* face = startFace;
      apf::MeshEntity* tet = 0;
      model = 0;
      size = 0;
      while (1)
      {
        apf::MeshEntity* v = ma::getTriVertOppositeEdge(mesh,face,edge);
/* there may be more than MAX_VERTS: just walk over them.
   overflow will be checked by users of SwapLoop
   (i.e. SwapCavity::findGoodTriangulation) */
        if (size < MAX_VERTS)
          verts[size] = v;
        ++size;
        if ( ! tet)
          tet = findInitialTet(v,face);
        else
          tet = getOtherTet(tet,face);
        if ( ! tet)
          break;
/* all the tets here should be in the same model region,
   record it for classifying new ones */
        if ( ! model)
          model = mesh->toModel(tet);
        else if (model != mesh->toModel(tet))
          break; /* this is analogous to (!tet) when
                    we cross a non-manifold face */
        face = getOtherFace(face,tet);
        assert(face);
        if (face == startFace)
          break;
      }
    }
    apf::MeshEntity* getOtherTet(apf::MeshEntity* oldTet, apf::MeshEntity* face)
    {
      apf::Up tets;
      mesh->getUp(face,tets);
      for (int i=0; i < tets.n; ++i)
      {
        apf::MeshEntity* newTet = tets.e[i];
        if (newTet != oldTet)
          return newTet;
      }
      return 0;
    }
    apf::MeshEntity* getOtherFace(apf::MeshEntity* oldFace, apf::MeshEntity* tet)
    {
      apf::MeshEntity* tf[4];
      mesh->getDownward(tet,2,tf);
      for (int i=0; i < 4; ++i)
      {
        if (tf[i] == oldFace) continue;
        apf::MeshEntity* fe[3];
        mesh->getDownward(tf[i],1,fe);
        if (findIn(fe,3,edge) != -1)
          return tf[i];
      }
      return 0;
    }
    int getSize() {return size;}
    apf::MeshEntity* getVert(int i) {return verts[i];}
    apf::MeshEntity* getEdgeVert(int i) {return edge_verts[i];}
    apf::ModelEntity* getModel() {return model;}
  private:
    apf::Mesh2* mesh;
    apf::MeshEntity* edge;
    apf::MeshEntity* edge_verts[2];
    int size;
    apf::MeshEntity* verts[7];
    apf::ModelEntity* model;
};

/* this class represents the full cavity around the
   loop of vertices. It is responsible for trying
   new triangulations as efficiently as possible
   and creating the first one that works */
class MaCavity
{
  public:
    void init(apf::Mesh2* m)
    {
      mesh = m;
      loop.init(m);
    }
    bool setFromEdge(apf::MeshEntity* edge)
    {
      return setFromEdgeAndFace(edge,mesh->getUpward(edge,0));
    }
/* returns true iff there are tets in the cavity */
    bool setFromEdgeAndFace(apf::MeshEntity* edge, apf::MeshEntity* face)
    {
      loop.setEdge(edge);
      loop.findFromFace(face);
      return loop.getSize() > 1;
    }
  private:
    apf::Mesh2* mesh;
    MaLoop loop;
};

class MaSwap3Dcheck {
  public:
    MaSwap3Dcheck(apf::Mesh2* m);
    ~MaSwap3Dcheck();
    bool run(apf::MeshEntity* e);
    bool run_all_surf(apf::MeshEntity* e);
  private:
    apf::Mesh2* mesh;
    apf::MeshEntity* edge;
    MaCavity halves[2];
    bool cavityExists[2];
};

// ma::adapt() sometimes crashes due to non-existing cavities. Check if this
// results from ma processes or insertion or collapses on our end.
// Uses MaSwap3Dcheck
bool chk_ma_swap(apf::Mesh2* m);

bool chk_ma_swap(apf::Mesh2* m, std::vector<apf::MeshEntity*> &vert_in);

bool chk_ma_swap_all(apf::Mesh2* m);

bool chk_ma_swap_all(apf::Mesh2* m, std::vector<apf::MeshEntity*> &vert_in);

// A wrapper for ma::Input that uses INPUT_TYPE to specify the fields of ma::Input
class vd_ma_input {
  public:
    int maximumIterations;
    bool shouldCoarsen;
    bool shouldFixShape;
    bool shouldForceAdaptation;
    bool shouldPrintQuality;
    bool shouldCheckQualityForDoubleSplits;
    double validQuality;
    double goodQuality;
    double maximumImbalance;
    double maximumEdgeRatio;
    bool shouldRunPreZoltan;
    bool shouldRunPreZoltanRib;
    bool shouldRunPreParma;
    bool shouldRunMidZoltan;
    bool shouldRunMidParma;
    bool shouldRunPostZoltan;
    bool shouldRunPostZoltanRib;
    bool shouldRunPostParma;
    bool shouldTurnLayerToTets;
    bool shouldCleanupLayer;
    bool shouldRefineLayer;
    bool shouldCoarsenLayer;
    bool splitAllLayerEdges;

    vd_ma_input();
    ~vd_ma_input();

    // Set the flags according to the input type.
    void set_def();
    void set_template(INPUT_TYPE IT_IN);
    void flip_template(INPUT_TYPE IT_IN);

    // Get the input fields from ma::Input.
    void get_input(ma::Input* in);
    // Set the input fields to ma::Input.
    void set_input(ma::Input* in);

    // Copy constructor
    vd_ma_input(const vd_ma_input& that);
    // Copy
    vd_ma_input& operator=(const vd_ma_input& that);

};

// Given an ma::Input, replace the old size field. 
void repl_sz_field(ma::Input* in, apf::Mesh2* m, MA_SIZE_TYPE MA_T);

// A wrapper for ma::Input object with common predefined use cases.
// Using INPUT_TYPE flags as input, the input flags can be set for a certain
// option.
// The order of INPUT_TYPEs is preserved in adaptation and each subadaptation is
// done individually.
class vd_input_iso {
  private:

  public:
    std::vector<INPUT_TYPE> IT;
    //vd_input_iso(apf::Mesh2* m, ma::IsotropicFunction* sf, ma::IdentitySizeField* ref = NULL);
    vd_input_iso();
    ~vd_input_iso();

    //void set_valid(double v_in);
    //void set_good(double g_in);

    //void set_sizefield(ma::IdentitySizeField* ref);
    // Set the input flags such that only the flags associated with the 
    // INPUT_TYPE are on. The order of adaptation is preserved.
    void set_input(std::vector<INPUT_TYPE> &IT_in);
    // Get the pointer to the input: 
    //ma::Input* get_input();
};

enum class ADAPT_TYPE {
  ADAPT_1CELL,
  ADAPT_3CELL,
  ADAPT_BOUND,
  ADAPT_CURVE,
  ADAPT_STEP,
  ADAPT_STEP_1CELL,
  END
};

/*
// Used to keep options for how to adapt the mesh. 
class ad_opts {
  public:
    // Parameters
    bool ad_flag;
    double ad_param;
    ADAPT_TYPE ad_type;
    // Constructor
    ad_opts();
    ad_opts(bool flag_in, double param_in, ADAPT_TYPE type_in);
    // Copy constructor
    ad_opts(const ad_opts& that);
    // Copy
    ad_opts& operator=(const ad_opts& that);

    ~ad_opts();
};

class vd_adapter {
  public:
    // Parameters
    ad_opts opts;
    // Constructor
    vd_adapter();
    vd_adapter(ad_opts opts_in);
    // Copy constructor
    vd_adapter(const vd_adapter& that);
    // Copy
    vd_adapter& operator=(const vd_adapter& that);

    ~vd_adapter();
};
*/

// Container for keeping calculations of distances of a given stratum vertices 
// from its bounding strata. Calculates for bounding 2s for 3c, 1s for 2s and
// 0s for 1s. 
// Assumes strata are not ring-like in all dimensions 
// (ring for 1s, sphere for 2s, projection of the 3sphere for 3s).
class DistBurner {
  public:
    apf::Mesh2* m;
    cell_base* c_base;

    double sim_len;
    double dist_max;

    // A point on the boundary entities. Position of one of the bounding 
    // vertices.
    std::map<apf::MeshEntity*, apf::Vector3> v_pos;
    std::map<apf::MeshEntity*, apf::Vector3> e_pos;
    std::map<apf::MeshEntity*, apf::Vector3> t_pos;

    // e_dir points towards the edge interior from the e_pos.
    std::map<apf::MeshEntity*, apf::Vector3> e_dir;
    // t_dir points towards the 3stratum tetrahedra from the t_pos.
    std::map<apf::MeshEntity*, apf::Vector3> t_dir;

    // Distance field over vertices.
    std::map<apf::MeshEntity*, double> v_dist;

    // The list of burning entities, bounded by c_dim burned and one unburned 
    // vertex
    std::vector<apf::MeshEntity*> front;

    // The 1-dim lower adjacencies of the bounding stratum entities.
    std::map<apf::MeshEntity*, std::vector<apf::MeshEntity*> > e_d;
    // The vertex adjacencies of the bounding stratum entities. For bounding 1s
    // they are the same so no need to allocate.
    std::map<apf::MeshEntity*, std::vector<apf::MeshEntity*> > e_v;

    // The vertex associated with the c-dim dimensional stratum entity. 
    std::map<apf::MeshEntity*, apf::MeshEntity* > v_map;
    std::map<apf::MeshEntity*, int > v_id;

    // For the current front entity and vertex, keep the map of the considered
    // boundary entities.
    std::map<apf::MeshEntity*, bool > burn_curr;

    // Adjacency maps for 1-dim lower adjacencies of the c_dim-1 dimensional 
    // bounding entities. e.g. burning 3stratum vertices, bounding entities are
    // triangles, which will be bounded by boundary edges. v_map is from edges
    // to triangles on the boundary. 
    // Unless self tangential at an edge or vertex, 3strata will have exactly 
    // two per edge. 2s might not be homological to a disc, in which case
    // there will be some single adjacency boundary entities. 
    std::map<apf::MeshEntity*, apf::MeshEntity* > e_map1;
    std::map<apf::MeshEntity*, apf::MeshEntity* > e_map2;

    // The burned entities. Vertices are burnt when the equidistance surface
    // passes through. The others are burnts once the associated vertex is 
    // burnt or the entity is only bound by burnt vertices.
    std::map<apf::MeshEntity*, bool> b_map;

    // Set the cell for the distances to be calculated:
    void set_cell(int d_in, int id_in, apf::Mesh2* m_in, vd_entlist & e_list, 
                                      cell_base* c_base_in, double len);
    // Collect adjacencies:
    void collect_adj(vd_entlist & e_list);

    // Starting from the current entity on the front and moving over the 
    // adjacent entities, find the smallest distance boundary element 
    // reachable from the starting element for the vertex associated 
    // with the starting entity.
    void calc_v_3c(apf::MeshEntity* tet);
    void calc_v_2c(apf::MeshEntity* tri);
    void calc_v_1c(apf::MeshEntity* edge);

    apf::MeshEntity* tri_proj_v(apf::MeshEntity* v, 
                      apf::MeshEntity* tri, apf::Vector3& int_pt);

    double calc_d_tri(apf::MeshEntity* v, apf::MeshEntity* tri);
    double calc_d_edge(apf::MeshEntity* v, apf::MeshEntity* edge);
    double calc_d_vert(apf::MeshEntity* v, apf::MeshEntity* vert);

    // Burn the vertices below (1+dist_tol)*dist_min distance from the boundary
    // remove the burnt entities in the front and add new ones.
    void step_front();

    // Run step_front until the front is empty.
    void burn_dist();

    int c_dim;
    int c_id;

    // The percent difference of the maximum distance being burned to the  
    // minimum distance being burned in the current loop.  
    double dist_tol;

    void clear();

    DistBurner();
    ~DistBurner();

};

// Used to set the condition to split edges.
class ModelEdgeRefiner : public ma::IdentitySizeField {
  public:
    // Coarse map, tagging edges of the cells not to be split.
    std::vector<std::map<int, bool> > coarse_map;
    bool split_all;

  ModelEdgeRefiner(ma::Mesh* m);
  bool shouldCollapse(ma::Entity* edge);
  bool shouldSplit(ma::Entity* edge);
  ~ModelEdgeRefiner();
    // TODO: if entity is a mesh edge spanning a model edge return true
};

// Used to set the condition to split edges.
class ModelEdgeRefinerDist : public ma::IdentitySizeField {
  public:
    // Coarse map, tagging edges of the cells not to be split.
    std::vector<std::map<int, bool> > coarse_map;
    apf::Field* field_step;
    bool split_all;
    double split_target;
    double coarse_target;

  ModelEdgeRefinerDist(ma::Mesh* m);
  double measure_dist(ma::Entity* edge);
  bool shouldCollapse(ma::Entity* edge);
  bool shouldSplit(ma::Entity* edge);
  void set_target_th(double coarse_in, double split_in);
  ~ModelEdgeRefinerDist();
    // TODO: if entity is a mesh edge spanning a model edge return true
};

// Used to only collapse edges.
class ModelEdgeCollapse : public ma::IdentitySizeField {
  public:
    // Coarse map, tagging edges of the cells not to be split.
    std::map<apf::MeshEntity*, bool > coarse_map;

  ModelEdgeCollapse(ma::Mesh* m);
  bool shouldCollapse(ma::Entity* edge);
  bool shouldSplit(ma::Entity* edge);
  ~ModelEdgeCollapse();
    // TODO: if entity is a mesh edge spanning a model edge return true
};

// Used to only collapse edges.
class ModelEdgeSplit : public ma::IdentitySizeField {
  public:
    // Coarse map, tagging edges of the cells not to be split.
    std::map<apf::MeshEntity*, bool > split_map;

  ModelEdgeSplit(ma::Mesh* m);
  bool shouldCollapse(ma::Entity* edge);
  bool shouldSplit(ma::Entity* edge);
  ~ModelEdgeSplit();
    // TODO: if entity is a mesh edge spanning a model edge return true
};

class ModelEdgeRefinerGrain : public ma::IdentitySizeField {
  public:
    // Coarse map, tagging edges of the cells not to be split.
    std::vector<std::map<int, bool> > coarse_map;
    apf::Field* field_step;
    bool split_all;

  ModelEdgeRefinerGrain(ma::Mesh* m);
  double measure_dist(ma::Entity* edge);
  bool shouldCollapse(ma::Entity* edge);
  bool shouldSplit(ma::Entity* edge);
  ~ModelEdgeRefinerGrain();
    // TODO: if entity is a mesh edge spanning a model edge return true
};

// Var
class ModelEdgeRefinerVarying : public ma::IdentitySizeField {
  public:
    // Equivalent length map, per cell of dim > 0. 
    std::vector<std::map<int, double> > len_map;
    // len*rat_map*len_map gives the effective length of the edge. 
    // rat_map is 1 for all existing cells, by default.
    std::vector<std::map<int, double> > rat_map;

    double split_target;
    double coarse_target;

  ModelEdgeRefinerVarying(ma::Mesh* m);
  void set_all_cell(double len);
  void set_bound_cell(double len, std::vector<std::pair<int,int> >* cells);
  double measure_dist(ma::Entity* edge);
  bool shouldCollapse(ma::Entity* edge);
  bool shouldSplit(ma::Entity* edge);
  void set_target_th(double coarse_in, double split_in);

  ~ModelEdgeRefinerVarying();
    // TODO: if entity is a mesh edge spanning a model edge return true
};

// Here is a template for meshadapt:
class Linear : public ma::IsotropicFunction
{
  public:
    Linear(ma::Mesh* m, double sz)
    {
      sizing = sz;

      mesh = m;
      average = ma::getAverageEdgeLength(m);
    }
    virtual double getValue(ma::Entity* v)
    {
      //return average*sizing;
      return 1/sizing;
    }
  private:
    ma::Mesh* mesh;
    double average;
    double sizing;
};

// Stepwise graded refinement. Each vertex has a different weight that is 
// divided by to the reference length, sizing.
// Default: The vertices on the boundaries have weight 1,
// adjacent vertices have 2, and the farthest vertices have 3.
// Based on this, the interior is refined less.
class Step : public ma::IsotropicFunction
{
  public:
    Step(ma::Mesh* m, double sz);
    double getValue(ma::Entity* v);
  private:
    ma::Mesh* mesh;
    double average;
    double sizing;
    apf::Field* field_step;
};

// Stepwise graded refinement. .
class Step_ns : public ma::IsotropicFunction
{
  public:
    Step_ns(ma::Mesh* m, double sz);
    double getValue(ma::Entity* v);
  private:
    ma::Mesh* mesh;
    double average;
    double sizing;
    apf::Field* field_step;
};

// Stepwise graded refinement using anisotropic function.
class Step_aniso : public ma::AnisotropicFunction
{
  public:
    Step_aniso(ma::Mesh* m, double sz);
    double getValue(ma::Entity* v);
  private:
    ma::Mesh* mesh;
    double average;
    double sizing;
    apf::Field* field_step;
};


void vd_adapt_0cell(apf::Mesh2* m, struct cell_base* c);
void vd_adapt(apf::Mesh2* m);

// Collect 0cells, collect their positions. Keep the list of 0cells of each 3cell. Keep the list of bounded 3cells for each 1cell and 2cell. Check cell of every edge. Check the distance to bounding 0cells of the 3cells the edge bounds. Depending on the closest distance, scale the length. The shortest edge should have the length 2*dist. If beyond a threshold(avg length/2), just keep it average length.

// TODO There are many operations across all entities such as this or minimum 
// time calculation. There might be ways of combining them without making 
// things too complicated.

class Linear_0cell : public ma::IsotropicFunction
{
  public:
    Linear_0cell(ma::Mesh* m, struct cell_base* c,
                                          std::vector<double>* avg_cell_in);
    void reload(ma::Mesh* m, struct cell_base* c, 
                                          std::vector<double>* avg_cell_in);

    virtual double getValue(ma::Entity* v);
    double getDistance(ma::Entity* v);
    ~Linear_0cell();

  private:
    struct cell_base* c_base;
    ma::Mesh* mesh;
    double average;
    std::vector<double> avg_cell;

    // 0cell positions.
    std::vector<apf::MeshEntity*> cell0;
    std::vector<apf::Vector3> cell0_pos;
    // 0cell adjacencies of 1- and 2-cells and 3-cells.
    std::vector<std::vector<std::vector<int > > > cnc0;

    void refresh();
    void load_0cell_ent();
    void load_0cell_pos();
    void load_0cell();
    void clear();

};


#endif
