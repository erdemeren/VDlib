#include <tgmath.h> 
#include<algorithm>
#include "apf.h"
#include "ma.h"
#include "maMesh.h"

#include "topo_extinfo.h"
#include "topo_geom.h"
#include "topo_ma.h"


// Get the lowest quality element and it's quality.
std::pair<apf::MeshEntity*, double> calc_valid_q(apf::Mesh2* m, 
                                            ma::SizeField* s_f) {
  apf::MeshEntity* e_valid = NULL;
  double valid = 1;
  apf::MeshEntity* elem;
  apf::MeshIterator* it = m->begin(3);
  while(elem = m->iterate(it)) {
    double q_temp = measureTetQuality(m, s_f, elem);
    if(q_temp < valid) {
      e_valid = elem;
      valid = q_temp;
    }
  }
  m->end(it);

  if(valid < std::numeric_limits<double>::min()) {
    valid = std::fabs(valid);
  }

  return std::make_pair(e_valid,valid/5);
}

// Collect elements with quality less than threshold.
double get_low_q(apf::Mesh2* m, ma::SizeField* s_f, std::vector<apf::MeshEntity*> &tets, double q_th) {
  std::map<apf::MeshEntity*, bool> low{};
  std::map<apf::MeshEntity*, double> qual{};
  int nbr = 0;

  double valid = 1;
  tets.clear();
  apf::MeshEntity* elem;
  apf::MeshIterator* it = m->begin(3);
  while(elem = m->iterate(it)) {
    double q_temp = measureTetQuality(m, s_f, elem);
    if(q_temp < q_th) {
      nbr = nbr + 1;
      low[elem] = true;
      qual[elem] = q_temp;
    }
    if(q_temp < valid)
      valid = q_temp;
  }
  m->end(it);

  tets.reserve(nbr);
  it = m->begin(3);
  while(elem = m->iterate(it)) {
    if(low[elem]) {
      tets.push_back(elem);
    }
  }
  m->end(it);

  if(valid < std::numeric_limits<double>::min()) {
    valid = std::fabs(valid);
  }

  return valid;
}

double calc_good_vol(apf::Mesh2* m, double ref_len, double m_len) {
  std::vector<apf::Vector3> v(4, apf::Vector3(0,0,0));
  apf::Downward d_v;
  apf::Downward d_e;
  apf::Downward d_s;
  v.at(0) = apf::Vector3(0,0,0);
  v.at(1) = apf::Vector3(m_len,0,0);
  v.at(2) = apf::Vector3(0,ref_len,0);
  v.at(3) = apf::Vector3(0,0,ref_len);

  // Create a new entity to calculate it's quailty.
  apf::MeshEntity* elem;
  apf::MeshIterator* it = m->begin(3);
  elem = m->iterate(it);
  m->end(it);
  apf::ModelEntity* mdl = m->toModel(elem);

  for(int i = 0; i < 4; i++) {
    d_v[i] = m->createVert(mdl);
    m->setPoint(d_v[i], 0, v.at(i));
  }

  elem = buildElement(m, mdl, apf::Mesh::TET, d_v, 0);

  m->getDownward(elem, 1, d_e);
  m->getDownward(elem, 2, d_s);

  double vol = vd_volume_tet(m, elem);

  // Destroy the entities:
  m->destroy(elem);
  for(int i = 0; i < 4; i++)
    m->destroy(d_s[i]);
  for(int i = 0; i < 6; i++)
    m->destroy(d_e[i]);
  for(int i = 0; i < 4; i++)
    m->destroy(d_v[i]);
  return vol;
}

bool vd_chk_vol_valid(apf::Mesh2* m, apf::MeshEntity* tet, apf::MeshEntity* ent) {
  int ent_type = m->getType(ent);
  int d = m->typeDimension[ent_type];
  assert(d > -1 and d < 3);
  bool valid = true;
  if(d > 0) {
    std::vector<apf::Vector3> v(4, apf::Vector3(0,0,0));
    apf::Downward d_v;
    apf::Downward d_v2;

    m->getDownward(tet, 0, d_v);
    std::map<apf::MeshEntity*, apf::Vector3> old_pos{};
    for(int i = 0; i < 4; i++) {
      m->getPoint(d_v[i], 0, v.at(i));
      old_pos[d_v[i]] = v.at(i);
    }

    int dc = m->getDownward(tet, d, d_v2);
    int e1 = findIn(d_v2, dc, ent);
    assert(e1 > -1);
    dc = m->getDownward(ent, 0, d_v2);

    apf::Vector3 midpoint(0,0,0);
    midpoint = vd_get_pos(m, ent);

    for(int j = 0; j < dc; j++) {
      m->setPoint(d_v2[j], 0, midpoint);
      bool vol = vd_volume_tet_sign(m, tet);
      valid = (valid and vol);
      m->setPoint(d_v2[j], 0, old_pos[d_v2[j]]);
    }
  }
  return valid;
}

double calc_good_q(apf::Mesh2* m, double ref_len, double m_len) {
  std::vector<apf::Vector3> v(4, apf::Vector3(0,0,0));
  apf::Downward d_v;
  apf::Downward d_e;
  apf::Downward d_s;
  v.at(0) = apf::Vector3(0,0,0);
  v.at(1) = apf::Vector3(m_len,0,0);
  v.at(2) = apf::Vector3(0,ref_len,0);
  v.at(3) = apf::Vector3(0,0,ref_len);

  // Create a new entity to calculate it's quailty.
  apf::MeshEntity* elem;
  apf::MeshIterator* it = m->begin(3);
  elem = m->iterate(it);
  m->end(it);
  apf::ModelEntity* mdl = m->toModel(elem);

  for(int i = 0; i < 4; i++) {
    d_v[i] = m->createVert(mdl);
    m->setPoint(d_v[i], 0, v.at(i));
  }

  elem = buildElement(m, mdl, apf::Mesh::TET, d_v, 0);

  m->getDownward(elem, 1, d_e);
  m->getDownward(elem, 2, d_s);

  ma::IdentitySizeField* ref_0c = new ma::IdentitySizeField(m);
  double qual = measureTetQuality(m, ref_0c, elem);
  delete ref_0c;

  // Destroy the entities:
  m->destroy(elem);
  for(int i = 0; i < 4; i++)
    m->destroy(d_s[i]);
  for(int i = 0; i < 6; i++)
    m->destroy(d_e[i]);
  for(int i = 0; i < 4; i++)
    m->destroy(d_v[i]);
  return qual;
}

double calc_q(apf::Mesh2* m, std::vector<apf::Vector3> &v) {
  apf::Downward d_v;
  apf::Downward d_s;
  apf::Downward d_e;

  // Create a new entity to calculate it's quailty.
  apf::MeshEntity* elem;
  apf::MeshIterator* it = m->begin(3);
  elem = m->iterate(it);
  m->end(it);
  apf::ModelEntity* mdl = m->toModel(elem);

  for(int i = 0; i < 4; i++) {
    d_v[i] = m->createVert(mdl);
    m->setPoint(d_v[i], 0, v.at(i));
  }

  elem = buildElement(m, mdl, apf::Mesh::TET, d_v, 0);

  m->getDownward(elem, 1, d_e);
  m->getDownward(elem, 2, d_s);

  ma::IdentitySizeField* ref_0c = new ma::IdentitySizeField(m);
  double qual = measureTetQuality(m, ref_0c, elem);
  delete ref_0c;

  // Destroy the entities:
  m->destroy(elem);
  for(int i = 0; i < 4; i++)
    m->destroy(d_s[i]);
  for(int i = 0; i < 6; i++)
    m->destroy(d_e[i]);
  for(int i = 0; i < 4; i++)
    m->destroy(d_v[i]);
  return qual;
}

MaSwap3Dcheck::MaSwap3Dcheck(apf::Mesh2* m): mesh(m) {
  halves[0].init(m);
  halves[1].init(m);
  edge = 0;
}

MaSwap3Dcheck::~MaSwap3Dcheck() {}

bool MaSwap3Dcheck::run(apf::MeshEntity* e) {
  if (ma::isOnModelEdge(mesh,e))
    return false;
  edge = e;
  if (ma::isOnModelFace(mesh,edge)) {
  }
  else {
    cavityExists[0] = halves[0].setFromEdge(edge);
    if(!cavityExists[0]) {
      MaLoop loop;
      loop.init(mesh);
      loop.setEdge(edge);

      std::vector<apf::MeshEntity*> tets(0);
      std::vector<apf::MeshEntity*> surf(0);
      vd_set_up(mesh, e, &surf);
      vd_set_up(mesh, &surf, &tets);

      apf::ModelEntity* mdl = mesh->toModel(edge);
      std::cout << "Curl is not OK! " << edge << " "
                << mesh->getModelType(mdl) << "c"
                << mesh->getModelTag(mdl) 
                << std::endl;
      apf::Downward dv;
      apf::Up up;
      mesh->getDownward(edge, 0, dv);
      mdl = mesh->toModel(dv[0]);
      std::cout << "Verts: " << dv[0] << " "
                << mesh->getModelType(mdl) << "c"
                << mesh->getModelTag(mdl) 
                << std::endl;

      mdl = mesh->toModel(dv[1]);
      std::cout << dv[1] << " "
                << mesh->getModelType(mdl) << "c"
                << mesh->getModelTag(mdl) 
                << std::endl;

      for(int i = 0; i < tets.size(); i++) {
        mdl = mesh->toModel(tets.at(i));
        std::cout << "Tet " << tets.at(i) << " "
                  << mesh->getModelType(mdl) << "c" << mesh->getModelTag(mdl)
                  << " " << vd_get_pos(mesh, tets.at(i))
                  << std::endl;
        vd_print_down(mesh, tets.at(i));
      }
      for(int i = 0; i < surf.size(); i++) {
        mdl = mesh->toModel(surf.at(i));
        std::cout << "Surf " << surf.at(i)  << " "
                  << mesh->getModelType(mdl) << "c" << mesh->getModelTag(mdl)
                  << " " << vd_get_pos(mesh, surf.at(i))
                  << std::endl;
        vd_print_down(mesh, surf.at(i));

        vd_set_up(mesh, surf.at(i), &tets);
        apf::MeshEntity* v = ma::getTriVertOppositeEdge(mesh,surf.at(i),edge);
        for(int j = 0; j < tets.size(); j++) {
          mdl = mesh->toModel(tets.at(j));
          std::cout << "Tet " << tets.at(j) << " "
                    << mesh->getModelType(mdl) << "c" << mesh->getModelTag(mdl)
                    << std::endl;
          std::cout << "Curl:" << loop.isCurlOk(v,tets.at(j)) << std::endl;
        }
      }
    }
    assert(cavityExists[0]);
  }
  return true;
}

bool MaSwap3Dcheck::run_all_surf(apf::MeshEntity* e) {
  bool curl_ok = true;
  if (ma::isOnModelEdge(mesh,e))
    return curl_ok;
  edge = e;
  if (ma::isOnModelFace(mesh,edge)) {
  }
  else {
    MaLoop loop;
    loop.init(mesh);
    loop.setEdge(edge);

    std::vector<apf::MeshEntity*> tets(0);
    std::vector<apf::MeshEntity*> surf(0);
    vd_set_up(mesh, e, &surf);
    vd_set_up(mesh, &surf, &tets);

    apf::ModelEntity* mdl = mesh->toModel(edge);
    std::cout << edge << " "
              << mesh->getModelType(mdl) << "c"
              << mesh->getModelTag(mdl) 
              << std::endl;
    apf::Downward dv;
    apf::Up up;
    mesh->getDownward(edge, 0, dv);
    mdl = mesh->toModel(dv[0]);
    std::cout << "Verts: " << dv[0] << " "
              << mesh->getModelType(mdl) << "c"
              << mesh->getModelTag(mdl) 
              << std::endl;

    mdl = mesh->toModel(dv[1]);
    std::cout << dv[1] << " "
              << mesh->getModelType(mdl) << "c"
              << mesh->getModelTag(mdl) 
              << std::endl;

    for(int i = 0; i < tets.size(); i++) {
      mdl = mesh->toModel(tets.at(i));
      std::cout << "Tet " << tets.at(i) << " "
                << mesh->getModelType(mdl) << "c" << mesh->getModelTag(mdl)
                << " " << vd_get_pos(mesh, tets.at(i))
                << std::endl;
      vd_print_down(mesh, tets.at(i));
    }
    for(int i = 0; i < surf.size(); i++) {
      mdl = mesh->toModel(surf.at(i));
      std::cout << "Surf " << surf.at(i)  << " "
                << mesh->getModelType(mdl) << "c" << mesh->getModelTag(mdl)
                << " " << vd_get_pos(mesh, surf.at(i))
                << std::endl;
      vd_print_down(mesh, surf.at(i));

      vd_set_up(mesh, surf.at(i), &tets);
      bool curls_tot = false;
      apf::MeshEntity* v = ma::getTriVertOppositeEdge(mesh,surf.at(i),edge);
      for(int j = 0; j < tets.size(); j++) {
        mdl = mesh->toModel(tets.at(j));
        std::cout << "Tet " << tets.at(j) << " "
                  << mesh->getModelType(mdl) << "c" << mesh->getModelTag(mdl)
                  << std::endl;
        bool curl_curr = loop.isCurlOk(v,tets.at(j));
        std::cout << "Curl:" << curl_curr << std::endl;
        curls_tot = curls_tot or curl_curr;
      }
      curl_ok = curl_ok and curls_tot;
    }
  }
  return curl_ok;
}

bool chk_ma_swap(apf::Mesh2* m) {
  MaSwap3Dcheck swap_chk(m);

  apf::MeshEntity* e;
  apf::ModelEntity* mdl;
  std::cout << "Doing swap check for cavityexists crashes" << std::endl;

  apf::MeshIterator* it = m->begin(1);
  while ((e = m->iterate(it))) {
    mdl = m->toModel(e);
    if(m->getModelType(mdl) == 3) {
      if(!swap_chk.run(e)) {
        return false;
      }
    }
  }
  m->end(it);
  return true;
}

bool chk_ma_swap(apf::Mesh2* m, std::vector<apf::MeshEntity*> &vert_in) {
  MaSwap3Dcheck swap_chk(m);

  apf::ModelEntity* mdl;
  std::cout << "Doing swap check for cavityexists crashes" << std::endl;

  std::vector<apf::MeshEntity*> edge(0);
  vd_set_up(m, &vert_in, &edge);
  apf::MeshEntity* e;
  for(int i = 0; i < edge.size(); i++) {
    e = edge.at(i);
    mdl = m->toModel(e);
    if(m->getModelType(mdl) == 3) {
      if(!swap_chk.run(e)) {
        return false;
      }
    }
  }
  return true;
}

bool chk_ma_swap_all(apf::Mesh2* m) {
  MaSwap3Dcheck swap_chk(m);

  apf::MeshEntity* e;
  apf::ModelEntity* mdl;
  std::cout << "Doing swap check for cavityexists crashes" << std::endl;

  bool curl_ok = true;

  apf::MeshIterator* it = m->begin(1);
  while ((e = m->iterate(it))) {
    mdl = m->toModel(e);
    if(m->getModelType(mdl) == 3) {
      bool curl_curr = swap_chk.run_all_surf(e);
      curl_ok = curl_ok and curl_curr;
    }
  }
  m->end(it);
  return curl_ok;
}

bool chk_ma_swap_all(apf::Mesh2* m, std::vector<apf::MeshEntity*> &vert_in) {
  MaSwap3Dcheck swap_chk(m);

  apf::ModelEntity* mdl;
  std::cout << "Doing swap check for cavityexists crashes" << std::endl;

  bool curl_ok = true;

  std::vector<apf::MeshEntity*> edge(0);
  vd_set_up(m, &vert_in, &edge);
  apf::MeshEntity* e;
  for(int i = 0; i < edge.size(); i++) {
    e = edge.at(i);
    mdl = m->toModel(e);
    if(m->getModelType(mdl) == 3) {
      bool curl_curr = swap_chk.run_all_surf(e);
      curl_ok = curl_ok and curl_curr;
    }
  }
  return curl_ok;
}

// Given an ma::Input, replace the old size field. 
void repl_sz_field(ma::Input* in, apf::Mesh2* m, MA_SIZE_TYPE MA_T) {
  if(MA_T == MA_SIZE_TYPE::EDGE_COL) {
    delete in->sizeField;
    ModelEdgeCollapse* ref = new ModelEdgeCollapse(m);
    in->sizeField = ref;
  }
  else if(MA_T == MA_SIZE_TYPE::EDGE_SPLIT) {
    delete in->sizeField;
    ModelEdgeSplit* ref = new ModelEdgeSplit(m);
    in->sizeField = ref;
  }
  else if(MA_T == MA_SIZE_TYPE::EDGE_REFINER) {
    delete in->sizeField;
    ModelEdgeRefiner* ref = new ModelEdgeRefiner(m);
    in->sizeField = ref;
  }
  else {
    assert(MA_T == MA_SIZE_TYPE::EDGE_REF_DIST);
    delete in->sizeField;
    ModelEdgeRefinerDist* ref = new ModelEdgeRefinerDist(m);
    in->sizeField = ref;
  }
}

void vd_ma_input::set_def() {
  maximumIterations = 3;
  shouldCoarsen = true;
  shouldFixShape = true;
  shouldForceAdaptation = false;
  shouldPrintQuality = true;
  goodQuality = 0.027;
  maximumEdgeRatio = 2.0;
  shouldCheckQualityForDoubleSplits = false;
  validQuality = 1e-10;
  maximumImbalance = 1.10;
  shouldRunPreZoltan = false;
  shouldRunPreZoltanRib = false;
  shouldRunPreParma = false;
  shouldRunMidZoltan = false;
  shouldRunMidParma = false;
  shouldRunPostZoltan = false;
  shouldRunPostZoltanRib = false;
  shouldRunPostParma = false;
  shouldTurnLayerToTets = false;
  shouldCleanupLayer = false;
  shouldRefineLayer = false;
  shouldCoarsenLayer = false;
  splitAllLayerEdges = false;
}

void vd_ma_input::flip_template(INPUT_TYPE IT_IN) {
  assert(IT_IN < INPUT_TYPE::END);
  if(IT_IN == INPUT_TYPE::ADAPT_SPLIT) {
    shouldRefineLayer = true;
  }
  else if(IT_IN == INPUT_TYPE::ADAPT_COL) {
    shouldRunPreZoltan = false;
    shouldRunMidParma = false;
    shouldRunPostParma = false;
    shouldRefineLayer = false;

    shouldCoarsen = true;
    shouldCoarsenLayer = true;

    shouldFixShape = false;
  }
  else {
    assert(IT_IN == INPUT_TYPE::ADAPT_Q);
    shouldRunPreZoltan = true;
    shouldRunMidParma = true;
    shouldRunPostParma = true;
    shouldRefineLayer = true;

    shouldFixShape = true;
  }

}

void vd_ma_input::set_template(INPUT_TYPE IT_IN) {
  set_def();
  flip_template(IT_IN);
}

// Get the input fields from ma::Input.
void vd_ma_input::get_input(ma::Input* in) {
  maximumIterations = in->maximumIterations;
  shouldCoarsen = in->shouldCoarsen;
  shouldFixShape = in->shouldFixShape;
  shouldForceAdaptation = in->shouldForceAdaptation;
  shouldPrintQuality = in->shouldPrintQuality;
  goodQuality = in->goodQuality;
  maximumEdgeRatio = in->maximumEdgeRatio;
  shouldCheckQualityForDoubleSplits = in->shouldCheckQualityForDoubleSplits;
  validQuality = in->validQuality;
  maximumImbalance = in->maximumImbalance;
  shouldRunPreZoltan = in->shouldRunPreZoltan;
  shouldRunPreZoltanRib = in->shouldRunPreZoltanRib;
  shouldRunPreParma = in->shouldRunPreParma;
  shouldRunMidZoltan = in->shouldRunMidZoltan;
  shouldRunMidParma = in->shouldRunMidParma;
  shouldRunPostZoltan = in->shouldRunPostZoltan;
  shouldRunPostZoltanRib = in->shouldRunPostZoltanRib;
  shouldRunPostParma = in->shouldRunPostParma;
  shouldTurnLayerToTets = in->shouldTurnLayerToTets;
  shouldCleanupLayer = in->shouldCleanupLayer;
  shouldRefineLayer = in->shouldRefineLayer;
  shouldCoarsenLayer = in->shouldCoarsenLayer;
  splitAllLayerEdges = in->splitAllLayerEdges;
}

// Transfer the input fields to ma::Input.
void vd_ma_input::set_input(ma::Input* in) {
  in->maximumIterations = maximumIterations;
  in->shouldCoarsen = shouldCoarsen;
  in->shouldFixShape = shouldFixShape;
  in->shouldForceAdaptation = shouldForceAdaptation;
  in->shouldPrintQuality = shouldPrintQuality;
  in->goodQuality = goodQuality;
  in->maximumEdgeRatio = maximumEdgeRatio;
  in->shouldCheckQualityForDoubleSplits = shouldCheckQualityForDoubleSplits;
  in->validQuality = validQuality;
  in->maximumImbalance = maximumImbalance;
  in->shouldRunPreZoltan = shouldRunPreZoltan;
  in->shouldRunPreZoltanRib = shouldRunPreZoltanRib;
  in->shouldRunPreParma = shouldRunPreParma;
  in->shouldRunMidZoltan = shouldRunMidZoltan;
  in->shouldRunMidParma = shouldRunMidParma;
  in->shouldRunPostZoltan = shouldRunPostZoltan;
  in->shouldRunPostZoltanRib = shouldRunPostZoltanRib;
  in->shouldRunPostParma = shouldRunPostParma;
  in->shouldTurnLayerToTets = shouldTurnLayerToTets;
  in->shouldCleanupLayer = shouldCleanupLayer;
  in->shouldRefineLayer = shouldRefineLayer;
  in->shouldCoarsenLayer = shouldCoarsenLayer;
  in->splitAllLayerEdges = splitAllLayerEdges;
}

// Copy constructor
vd_ma_input::vd_ma_input(const vd_ma_input& that) {
  maximumIterations = that.maximumIterations;
  shouldCoarsen = that.shouldCoarsen;
  shouldFixShape = that.shouldFixShape;
  shouldForceAdaptation = that.shouldForceAdaptation;
  shouldPrintQuality = that.shouldPrintQuality;
  goodQuality = that.goodQuality;
  maximumEdgeRatio = that.maximumEdgeRatio;
  shouldCheckQualityForDoubleSplits = that.shouldCheckQualityForDoubleSplits;
  validQuality = that.validQuality;
  maximumImbalance = that.maximumImbalance;
  shouldRunPreZoltan = that.shouldRunPreZoltan;
  shouldRunPreZoltanRib = that.shouldRunPreZoltanRib;
  shouldRunPreParma = that.shouldRunPreParma;
  shouldRunMidZoltan = that.shouldRunMidZoltan;
  shouldRunMidParma = that.shouldRunMidParma;
  shouldRunPostZoltan = that.shouldRunPostZoltan;
  shouldRunPostZoltanRib = that.shouldRunPostZoltanRib;
  shouldRunPostParma = that.shouldRunPostParma;
  shouldTurnLayerToTets = that.shouldTurnLayerToTets;
  shouldCleanupLayer = that.shouldCleanupLayer;
  shouldRefineLayer = that.shouldRefineLayer;
  shouldCoarsenLayer = that.shouldCoarsenLayer;
  splitAllLayerEdges = that.splitAllLayerEdges;
}
// Copy
vd_ma_input& vd_ma_input::operator=(const vd_ma_input& that) {
  maximumIterations = that.maximumIterations;
  shouldCoarsen = that.shouldCoarsen;
  shouldFixShape = that.shouldFixShape;
  shouldForceAdaptation = that.shouldForceAdaptation;
  shouldPrintQuality = that.shouldPrintQuality;
  goodQuality = that.goodQuality;
  maximumEdgeRatio = that.maximumEdgeRatio;
  shouldCheckQualityForDoubleSplits = that.shouldCheckQualityForDoubleSplits;
  validQuality = that.validQuality;
  maximumImbalance = that.maximumImbalance;
  shouldRunPreZoltan = that.shouldRunPreZoltan;
  shouldRunPreZoltanRib = that.shouldRunPreZoltanRib;
  shouldRunPreParma = that.shouldRunPreParma;
  shouldRunMidZoltan = that.shouldRunMidZoltan;
  shouldRunMidParma = that.shouldRunMidParma;
  shouldRunPostZoltan = that.shouldRunPostZoltan;
  shouldRunPostZoltanRib = that.shouldRunPostZoltanRib;
  shouldRunPostParma = that.shouldRunPostParma;
  shouldTurnLayerToTets = that.shouldTurnLayerToTets;
  shouldCleanupLayer = that.shouldCleanupLayer;
  shouldRefineLayer = that.shouldRefineLayer;
  shouldCoarsenLayer = that.shouldCoarsenLayer;
  splitAllLayerEdges = that.splitAllLayerEdges;
  return *this;
}

vd_ma_input::vd_ma_input() : 
  maximumIterations(3),
  shouldCoarsen(true),
  shouldFixShape(true),
  shouldForceAdaptation(false),
  shouldPrintQuality(true),
  goodQuality(0.027),
  maximumEdgeRatio(2.0),
  shouldCheckQualityForDoubleSplits(false),
  validQuality(1e-10),
  maximumImbalance(1.10),
  shouldRunPreZoltan(false),
  shouldRunPreZoltanRib(false),
  shouldRunPreParma(false),
  shouldRunMidZoltan(false),
  shouldRunMidParma(false),
  shouldRunPostZoltan(false),
  shouldRunPostZoltanRib(false),
  shouldRunPostParma(false),
  shouldTurnLayerToTets(false),
  shouldCleanupLayer(false),
  shouldRefineLayer(false),
  shouldCoarsenLayer(false),
  splitAllLayerEdges(false) {
}
vd_ma_input::~vd_ma_input() {}

vd_input_iso::vd_input_iso() {
}
/*
// A wrapper for ma::Input object with common predefined use cases.
// Using INPUT_TYPE flags as input, the input flags can be set for a certain
// option, or the flags associated with a INPUT_TYPE can be flipped.
vd_input_iso::vd_input_iso(apf::Mesh2* m, ma::IsotropicFunction* sf, 
                                  ma::IdentitySizeField* ref) : 
                                                                in(NULL) {
  in = ma::configure(m, sf);

  if(ref != NULL) {
    delete in->sizeField;
    in->sizeField = ref;
  }

  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;

  in->shouldCoarsenLayer = true;

  in->shouldFixShape = true;
  in->goodQuality = 0.2;
  in->validQuality = 10e-5;

  std::cout << "Minimum quality is " << in->validQuality 
            << ", good quality is " << in->goodQuality  
            << std::endl;

  in->maximumEdgeRatio = 12;
  in->maximumIterations = 3;
}
*/

//void vd_input_iso::set_valid(double v_in) {
//  in->validQuality = v_in;
//  std::cout << "Set valid quality to " << in->validQuality 
//            << std::endl;
//}
//void vd_input_iso::set_good(double g_in) {
//  in->goodQuality = g_in;
//  std::cout << "Set good quality to " << in->goodQuality 
//            << std::endl;
//}

vd_input_iso::~vd_input_iso() {
}

//void vd_input_iso::set_sizefield(ma::IdentitySizeField* ref) {
//  assert(ref != NULL);
//  delete in->sizeField;
//  in->sizeField = ref;
//}

void vd_input_iso::set_input(std::vector<INPUT_TYPE> &IT_IN) {
  IT = IT_IN;
}

//ma::Input* vd_input_iso::get_input() {
//  return in;
//}

// Container for keeping calculations of distances of a given stratum vertices 
// from its bounding strata. Calculates for bounding 2s for 3c, 1s for 2s and
// 0s for 1s. 

// Set the cell for the distances to be calculated:
void DistBurner::set_cell(int d_in, int id_in, apf::Mesh2* m_in, 
                                                      vd_entlist & e_list, 
                                  cell_base* c_base_in, double len) {
  clear();
  c_dim = d_in;
  c_id = id_in;

  sim_len = len;

  m = m_in;
  c_base = c_base_in;

  collect_adj(e_list);

  burn_dist();
//consider the distance to entity centers, in addition to vertices. 
}

// Collect adjacencies:
void DistBurner::collect_adj(vd_entlist & e_list) {
  assert(e_list.e.at(c_dim).at(c_id).at(c_dim).size() > 0);
  front.reserve(e_list.e.at(c_dim).at(c_id).at(c_dim).size());

  // Link the downward adjacencies of the same dimensional entities of the
  // calculated stratum to it's entities of one dimension lower. 
  int d_lower = c_dim - 1;

  apf::Downward d_e;
  for(int i = 0; i < e_list.e.at(c_dim).at(c_id).at(c_dim).size(); i++) {
    apf::MeshEntity* e_curr = e_list.e.at(c_dim).at(c_id).at(c_dim).at(i);
    m->getDownward(e_curr, d_lower, d_e);
    for(int j = 0; j < c_dim + 1; j++) {
      if(e_map1[d_e[j]] == NULL) {
        e_map1[d_e[j]] = e_curr;
      }
      else {
        assert(e_map2[d_e[j]] == NULL);
        e_map2[d_e[j]] = e_curr;
      }
    }
  }

  ent_conn* edown = new ent_conn();
  c_base->get_conn(c_dim, c_id, edown);

  // TODO for now, assume no ring-like strata. In future, consider the distance
  // from the weighted center.
  assert(edown->conn.size() > 0);

  // Collect the geometric information about boundary entities.
  int lookup_tet_x_surf [4] = {3, 2, 0, 1};
  int lookup_v_tri_x_e [3] = {2,0,1};

  apf::Vector3 temp(0,0,0);
  apf::Downward d_v;
  apf::Up up;

  for(int c_down = 0; c_down < e_list.e.at(d_lower).size(); c_down++) {
    if(e_list.e.at(d_lower).at(c_down).at(d_lower).size() > 0) {
      // Walk over downward adjancency lists:
      for(int i = 0; i < e_list.e.at(d_lower).at(c_down).at(d_lower).size(); 
                                                                        i++) {
        apf::MeshEntity* e_curr = e_list.e.at(d_lower).at(c_down)
                                                       .at(d_lower).at(i);
        m->getUp(e_curr, up);
        // The geometric information is assigned. Used in lower dim boundary
        // entities to prevent spurious assignment. 
        std::map<apf::MeshEntity*, bool> asgn{};
        bool found = false;
        for(int j = 0; j < up.n; j++) {
          apf::ModelEntity* mdl = m->toModel(up.e[j]);
          int d_up = m->getModelType(mdl);
          int id_up = m->getModelTag(mdl) - 1;
          if(d_up == c_dim and id_up == c_id) {
            assert(!found);
            found = true;
            front.push_back(up.e[j]);

            m->getDownward(up.e[j], d_lower, d_e);
            int e1 = findIn(d_e, c_dim+1, e_curr);
            // Collecting entities bounding a 1cell:
            if(d_lower == 0) {
              // The information about the bounding vertex is collected:
              m->getPoint(e_curr, 0, temp);
              v_pos[e_curr] = temp;
              asgn[e_curr] = true;

              // The vertex of the front entity is collected:
              int v1 = (e1 + 1) % 2;
              v_map[up.e[j]] = d_e[v1];
              v_id[up.e[j]] = v1;
              m->getPoint(d_e[v1], 0, temp);
              v_pos[d_e[v1]] = temp;
            }
            // Collecting entities bounding a 2cell:
            else if(d_lower == 1) {
              // The vertex of the front entity is collected:
              m->getDownward(up.e[j], 0, d_v);
              int v1 = lookup_v_tri_x_e[e1];
              v_map[up.e[j]] = d_v[v1];
              v_id[up.e[j]] = v1;

              // Position of the interior vertex.
              m->getPoint(d_v[v1], 0, temp);
              v_pos[d_v[v1]] = temp;

              // The information about the bounding entities are collected:
              apf::Vector3 temp2(0,0,0);
              int v2 = (v1+1)%3;
              m->getPoint(d_v[v2], 0, temp);
              v_pos[d_v[v2]] = temp;
              asgn[d_v[v2]] = true;

              v2 = (v1+2)%3;
              m->getPoint(d_v[v2], 0, temp2);
              v_pos[d_v[v2]] = temp2;
              asgn[d_v[v2]] = true;

              e_pos[e_curr] = temp;
              e_dir[e_curr] = norm_0(temp2 - temp);
              asgn[e_curr] = true;
            }
            else {
              // The vertex of the front entity is collected:
              m->getDownward(up.e[j], 0, d_v);
              int v1 = lookup_tet_x_surf[e1];
              v_map[up.e[j]] = d_v[v1];
              v_id[up.e[j]] = v1;

              // Position of the interior vertex.
              m->getPoint(d_v[v1], 0, temp);
              v_pos[d_v[v1]] = temp;

              // The information about the bounding entities are collected:
              t_pos[e_curr] = apf::Vector3(0,0,0);
              for(int k = 0; k < 3; k++) {
                int v2 = (v1+k+1)%4;
                m->getPoint(d_v[v2], 0, temp);
                v_pos[d_v[v2]] = temp;
                asgn[d_v[v2]] = true;
                t_pos[e_curr] = t_pos[e_curr] + temp;
              }

              m->getPoint(d_v[v1], 0, temp);
              t_dir[e_curr] = vd_area_out_n(m, e_curr);
              if((temp - t_pos[e_curr])*t_dir[e_curr] < 
                              - std::numeric_limits<double>::min())
                t_dir[e_curr] = t_dir[e_curr]*(-1);

              m->getDownward(e_curr, 1, d_e);
              for(int k = 0; k < 3; k++) {
                m->getDownward(d_e[k], 0, d_v);
                assert(asgn[d_v[0]]);
                assert(asgn[d_v[1]]);
                m->getPoint(d_v[0], 0, temp);
                e_pos[d_e[k]] = temp;

                m->getPoint(d_v[1], 0, temp);
                e_dir[d_e[k]] = norm_0(temp - e_pos[d_e[k]]);
                asgn[d_e[k]] = true;
              }
            }
          }
        }
        assert(found);
      }
    }
  }
  delete edown;
}

// Starting from the current entity on the front and moving over the 
// adjacent entities, find the smallest distance boundary element 
// reachable from the starting element for the vertex associated 
// with the starting entity.
void DistBurner::calc_v_3c(apf::MeshEntity* tet) {
  apf::Downward d_v;
/*
apf::Vector3 int_pt(0,0,0);
apf::MeshEntity* last;
apf::MeshEntity* next;
  next = tri_proj_v(v_map[tet], last, int_pt);
  if(last == next)
*/

}

void DistBurner::calc_v_2c(apf::MeshEntity* tri) {
}

void DistBurner::calc_v_1c(apf::MeshEntity* edge) {
}


// Starting from given triangle, return the entity the closest projection of 
// the vertex onto the boundary lies on. If the projection is on the triangle,
// return the intersection point.
apf::MeshEntity* DistBurner::tri_proj_v(apf::MeshEntity* v, apf::MeshEntity* tri, apf::Vector3& int_pt) {
  apf::MeshEntity* ent = tri;

  apf::Vector3 temp1(0,0,0);
  apf::Vector3 temp2(0,0,0);

  apf::Vector3 edir(0,0,0);
  apf::Vector3 epos(0,0,0);

  int_pt = v_pos[v] - t_pos[tri];
  int_pt = int_pt - t_dir[tri]*(int_pt*t_dir[tri]) + t_pos[tri];

  for(int i = 0; i < 3; i++) {
    
    edir = e_dir[e_d[tri].at(i)];
    epos = e_pos[e_d[tri].at(i)];

    temp1 = (int_pt - epos);
    temp1 = temp1 - edir*(temp1*edir);
    temp2 = (t_pos[tri] - epos);
    temp2 = temp2 - edir*(temp2*edir);
    // If outside get the first edge that the projection lies outside of.
    if( temp1*temp2 < - std::numeric_limits<double>::min()) {
      //std::cout << "v1v2 - t_pos[tri]" << temp1*temp2 << std::endl;

      if(e_map1[e_d[tri].at(i)] == tri)
        ent = e_map2[e_d[tri].at(i)];
      else {
        assert(e_map2[e_d[tri].at(i)] == tri);
        ent = e_map1[e_d[tri].at(i)];
      }
      i = 3;
    }
  }
  return ent;
}


// Calculate the distance from the vertex to a given entity.
double DistBurner::calc_d_tri(apf::MeshEntity* v, apf::MeshEntity* tri) {
  return 0.;
}

double DistBurner::calc_d_edge(apf::MeshEntity* v, apf::MeshEntity* edge) {
  return 0.;
}
double DistBurner::calc_d_vert(apf::MeshEntity* v, apf::MeshEntity* vert) {
  apf::Downward d_v;
  return 0.;
}


// Burn the vertices below (1+dist_tol)*dist_min distance from the boundary
// remove the burnt entities in the front and add new ones.
void DistBurner::step_front() {
}

// Run step_front until the front is empty.
void DistBurner::burn_dist() {
/*
  dist_max = sim_len;
  while(front.size() > 0) {
    // The end burnt number of tets are burnt and not to be considered.
    int burnt = 0;
    double dist_min = sim_len;
    for(int i = 0; i < front.size() - burnt; i++) {
      apf::MeshEntity* t_curr = front.at(i);
      assert(!t_burn[t_curr]);
      m->getDownward(t_curr, 0, d_v);
      bool found = false;
      for(int j = 0; j < 4; j++) {
        if(!v_burn[d_v[j]]) {
          assert(!found);
          found = true;
          v_map[t_curr] = d_v[j];
          v_id[t_curr] = j;
        }
      }
      assert(found);
      t_dist[t_curr] = approx_dist(m, v_dist, d_v, v_map[t_curr]);
      if(t_dist[t_curr] < v_dist[v_map[t_curr]]) {
        v_dist[v_map[t_curr]] = t_dist[t_curr];
        if(v_dist[v_map[t_curr]] < dist_min)
          dist_min = v_dist[v_map[t_curr]];
      }
    }
    int added = 0;
    double dist_th = dist_min*(1+dist_tol);
    for(int i = 0; i < front.size() - burnt - added; i++) {
      apf::MeshEntity* t_curr = front.at(i);
      if(v_dist[v_map[t_curr]] < dist_th) {
        v_burn[v_map[t_curr]] = true;
      }
    }

    for(int i = 0; i < front.size() - burnt - added; i++) {
      apf::MeshEntity* t_curr = front.at(i);
      assert(tv_id[t_curr] > -1);
      assert(v_map[t_curr] != NULL);

      if(v_dist[v_map[t_curr]] < dist_th) {
        m->getDownward(t_curr, 2, d_t);
        int v_id = tv_id[t_curr];
        for(int j = 0; j < 3; j++) {
          int t_id = lookup_tet_surf[v_id][j];
          apf::MeshEntity* tri_curr = d_t[t_id];

          m->getUp(tri_curr, up);
          if(up.n == 2) {
            apf::MeshEntity* t_next;
            if(up.e[0] == t_curr)
              t_next = up.e[1];
            else
              t_next = up.e[0];
            if(!t_burn[t_next]) {
              int t1 = findIn(&front, front.size(), t_next);
              if(t1 == -1) {
                m->getDownward(t_next, 0, d_v);
                bool found1 = false;
                bool found2 = false;
                for(int j = 0; j < 4; j++) {
                  if(!v_burn[d_v[j]]) {
                    if(found1)
                      found2 = true;
                    else
                      found1 = true;
                  }
                }
                assert(!found2);
                if(found1) {
                  if(burnt > 0) {
                    int i_b = front.size() - burnt;
                    assert(!t_burn[front.at(i_b-1)]);
                    assert(t_burn[front.at(i_b)]);
                    front.push_back(front.at(i_b));
                    front.at(i_b) = t_next;
                    //m->getDownward(t_next, 2, d_t2);
                    //int t2 = findIn(d_t2, 4, tri_curr);
                    //assert(t2 > -1);
                    //v_map[t_next] = d_v[lookup_tet_surf_x[t2]];
                    //t_dist[t_next] = approx_dist(m, v_dist, d_v, 
                    //                              v_map[t_curr]);
                  }
                  else
                    front.push_back(t_next);

                  tv_id[t_next] = -1;
                  v_map[t_next] = NULL;

                  added = added + 1;
                }
                else {
                  t_burn[t_next] = true;
                }
              }
            }
          }
        }
        t_burn[t_curr] = true;
        int i_ub = front.size() - burnt - added - 1;
        int i_ad = front.size() - burnt - 1;
        if(i_ub != i)
          assert(!t_burn[front.at(i_ub)]);
        if(added > 0)
          assert(!t_burn[front.at(i_ub+1)]);
        if(burnt > 0)
          assert(t_burn[front.at(i_ub+added+1)]);

        front.at(i) = front.at(i_ub);
        front.at(i_ub) = front.at(i_ad);
        front.at(i_ad) = t_curr;
        burnt = burnt + 1;
        i = i-1;
      }
    }
    front.resize(front.size() - burnt);
  }
*/

}


void DistBurner::clear() {
  v_pos.clear(); 
  e_pos.clear();
  t_pos.clear();
  e_dir.clear();
  t_dir.clear();
  v_dist.clear();
  front.clear();

  burn_curr.clear();
  e_d.clear();
  e_v.clear();

  v_map.clear();
  v_id.clear();

  e_map1.clear();
  e_map2.clear();
  b_map.clear();
}

DistBurner::DistBurner() : 
  m(NULL), v_pos{}, e_pos{}, t_pos{}, e_dir{}, t_dir{}, v_dist{},
  burn_curr{}, 
  front(0), e_d{}, e_v{}, v_map{}, v_id{}, e_map1{}, e_map2{}, b_map{} {
}

DistBurner::~DistBurner() {
  clear();
}


ModelEdgeRefiner::ModelEdgeRefiner(ma::Mesh* m) : coarse_map(3, std::map<int, bool> {}), IdentitySizeField(m), split_all(false) {
}

bool ModelEdgeRefiner::shouldCollapse(ma::Entity* edge) {
  int dim_e = mesh->getModelType(mesh->toModel(edge));

  if(coarse_map.at(dim_e - 1)[mesh->getModelTag(mesh->toModel(edge))])
    return true;
  return false;
}

bool ModelEdgeRefiner::shouldSplit(ma::Entity* edge) {
  apf::Downward d_v;
  mesh->getDownward(edge, 0, d_v);
  int dim_e = mesh->getModelType(mesh->toModel(edge));
  int dim_v0 = mesh->getModelType(mesh->toModel(d_v[0]));
  int dim_v1 = mesh->getModelType(mesh->toModel(d_v[1]));

  if(coarse_map.at(dim_e - 1)[mesh->getModelTag(mesh->toModel(edge))])
    return false;
  else if(dim_e > dim_v0 and dim_e > dim_v1) {
    if(dim_e == 1)
      return true;
    else if(split_all)
      return true;
  }
  return false;
}

ModelEdgeRefiner::~ModelEdgeRefiner() {
  for(int i = 0; i < 3; i++) {
    coarse_map.at(i).clear();
  }
  coarse_map.clear();
}

ModelEdgeRefinerDist::ModelEdgeRefinerDist(ma::Mesh* m) : coarse_map(3, std::map<int, bool> {}),  IdentitySizeField(m), split_all(false),
split_target(1.5), coarse_target(0.5) {
  field_step = m->findField("adapt_step");
  assert(field_step != NULL);
}

double ModelEdgeRefinerDist::measure_dist(ma::Entity* edge) {

  // TODO Use integrators of apf to integrate the field directly. For linear
  // elements and for isotropic metrics, it is equivalent. 
  apf::MeshElement* me = apf::createMeshElement(mesh, edge);
  double len = measure(edge);
  apf::destroyMeshElement(me);

  apf::Downward d_v;
  mesh->getDownward(edge, 0, d_v);
  double dist = apf::getScalar(field_step, d_v[0], 0);
  dist = (dist + apf::getScalar(field_step, d_v[1], 0))/2;
  assert(dist > std::numeric_limits<double>::min());
  //return average*dist*sizing;
  return len*dist;
}

bool ModelEdgeRefinerDist::shouldCollapse(ma::Entity* edge) {
  int dim_e = mesh->getModelType(mesh->toModel(edge));

  if(coarse_map.at(dim_e - 1)[mesh->getModelTag(mesh->toModel(edge))])
    return true;
  return this->measure_dist(edge) < coarse_target;
}

bool ModelEdgeRefinerDist::shouldSplit(ma::Entity* edge) {
  apf::Downward d_v;
  mesh->getDownward(edge, 0, d_v);
  int dim_e = mesh->getModelType(mesh->toModel(edge));
  int dim_v0 = mesh->getModelType(mesh->toModel(d_v[0]));
  int dim_v1 = mesh->getModelType(mesh->toModel(d_v[1]));
  //double v0 = apf::getScalar(field_step, d_v[0], 0);
  //double v1 = apf::getScalar(field_step, d_v[1], 0);

  //double len = vd_meas_ent(mesh, edge);
  //return (1 > v0) or (1 > v1);
  //return (2*len > v0) or (2*len > v1);
  //return (len > v0/4+v1/4);
  if(coarse_map.at(dim_e - 1)[mesh->getModelTag(mesh->toModel(edge))])
    return false;
  //else if(dim_e == 1 and dim_e > dim_v0 and dim_e > dim_v1)
  else if(dim_e > dim_v0 and dim_e > dim_v1) {
    if(dim_e == 1)
      return true;
    else if(split_all)
      return true;
  }
  return this->measure_dist(edge) > split_target;
  //else if
}

void ModelEdgeRefinerDist::set_target_th(double coarse_in, double split_in) {
  split_target = split_in;
  coarse_target = coarse_in;
}

ModelEdgeRefinerDist::~ModelEdgeRefinerDist() {
  for(int i = 0; i < 3; i++) {
    coarse_map.at(i).clear();
  }
  coarse_map.clear();
}

ModelEdgeCollapse::ModelEdgeCollapse(ma::Mesh* m) : coarse_map{}, 
                                                    IdentitySizeField(m) {
}

bool ModelEdgeCollapse::shouldCollapse(ma::Entity* edge) {
  int dim_e = mesh->getModelType(mesh->toModel(edge));

  if(coarse_map[edge])
    return true;
  return false;
}

bool ModelEdgeCollapse::shouldSplit(ma::Entity* edge) {
  return false;
}
ModelEdgeCollapse::~ModelEdgeCollapse() {

  coarse_map.clear();
}

ModelEdgeSplit::ModelEdgeSplit(ma::Mesh* m) : split_map{}, 
                                                    IdentitySizeField(m) {
}

bool ModelEdgeSplit::shouldCollapse(ma::Entity* edge) {
  return false;
}

bool ModelEdgeSplit::shouldSplit(ma::Entity* edge) {
  if(split_map[edge])
    return true;
  return false;
}

ModelEdgeSplit::~ModelEdgeSplit() {
  split_map.clear();
}

/////////////////////////////////////
// ModelEdgeRefinerGrain
/////////////////////////////////////
ModelEdgeRefinerGrain::ModelEdgeRefinerGrain(ma::Mesh* m) : coarse_map(3, std::map<int, bool> {}), IdentitySizeField(m), split_all(false) {
  field_step = m->findField("adapt_step");
  assert(field_step != NULL);
}

double ModelEdgeRefinerGrain::measure_dist(ma::Entity* edge) {

  // TODO Use integrators of apf to integrate the field directly. For linear
  // elements and for isotropic metrics, it is equivalent. 
  apf::MeshElement* me = apf::createMeshElement(mesh, edge);
  double len = measure(edge);
  apf::destroyMeshElement(me);

  apf::Downward d_v;
  mesh->getDownward(edge, 0, d_v);
  double dist = apf::getScalar(field_step, d_v[0], 0);
  dist = (dist + apf::getScalar(field_step, d_v[1], 0))/2;
  assert(dist > std::numeric_limits<double>::min());
  //return average*dist*sizing;
  return len*dist;
}

bool ModelEdgeRefinerGrain::shouldCollapse(ma::Entity* edge) {
  int dim_e = mesh->getModelType(mesh->toModel(edge));

  if(coarse_map.at(dim_e - 1)[mesh->getModelTag(mesh->toModel(edge))])
    return true;
  return this->measure_dist(edge) < 0.5;
}

bool ModelEdgeRefinerGrain::shouldSplit(ma::Entity* edge) {
  apf::Downward d_v;
  mesh->getDownward(edge, 0, d_v);
  int dim_e = mesh->getModelType(mesh->toModel(edge));
  int dim_v0 = mesh->getModelType(mesh->toModel(d_v[0]));
  int dim_v1 = mesh->getModelType(mesh->toModel(d_v[1]));
  //double v0 = apf::getScalar(field_step, d_v[0], 0);
  //double v1 = apf::getScalar(field_step, d_v[1], 0);

  //double len = vd_meas_ent(mesh, edge);
  //return (1 > v0) or (1 > v1);
  //return (2*len > v0) or (2*len > v1);
  //return (len > v0/4+v1/4);
  if(coarse_map.at(dim_e - 1)[mesh->getModelTag(mesh->toModel(edge))])
    return false;
  else if(dim_e > dim_v0 and dim_e > dim_v1) {
    if(dim_e == 1)
      return true;
    else if(split_all)
      return true;
  }
  return this->measure_dist(edge) > 1.5;
  //else if
}
ModelEdgeRefinerGrain::~ModelEdgeRefinerGrain() {
  for(int i = 0; i < 3; i++) {
    coarse_map.at(i).clear();
  }
  coarse_map.clear();
}

/////////////////////////////////////
// ModelEdgeRefinerVarying
/////////////////////////////////////
ModelEdgeRefinerVarying::ModelEdgeRefinerVarying(ma::Mesh* m) : len_map(3, std::map<int, double> {}), rat_map(3, std::map<int, double> {}), IdentitySizeField(m), split_target(1.5), coarse_target(0.5) {
  gmi_model* mdl = m->getModel();
  struct gmi_iter* it_gmi;
  struct gmi_ent* e_gmi;
  struct gmi_set* s_gmi;

  apf::MeshEntity* e;

  for(int dim = 1; dim < 4; dim++) {
    apf::MeshIterator* it = m->begin(dim);
    while ((e = m->iterate(it))) {
      apf::ModelEntity* mdl_curr = m->toModel(e);
      int c_type = m->getModelType(mdl_curr);
      int c_tag = m->getModelTag(mdl_curr);
      double meas_curr = vd_meas_ent(m, e);
      len_map.at(c_type - 1)[c_tag] = len_map.at(c_type - 1)[c_tag] + meas_curr;
    }
    m->end(it);
  }

  double pow[3] = {1., 0.5, 1.0/3};
  double coef[3] = {1./2, 1./8, 1./16};
  for(int dim = 1; dim < 4; dim++) {
    it_gmi = gmi_begin(mdl, dim);
    while ((e_gmi = gmi_next(mdl, it_gmi))) {
      int c_tag = gmi_tag(mdl,e_gmi);
      len_map.at(dim - 1)[c_tag] = std::pow(len_map.at(dim - 1)[c_tag]*coef[dim-1], pow[dim-1]);
      rat_map.at(dim - 1)[c_tag] = 1;
    }
    gmi_end(mdl, it_gmi);
  }
}

// Set the rat_map for boundary cells such that they are refined to yield edge 
// length len.
void ModelEdgeRefinerVarying::set_all_cell(double len) {

  gmi_model* mdl = mesh->getModel();
  struct gmi_iter* it_gmi;
  struct gmi_ent* e_gmi;
  struct gmi_set* s_gmi;

  for(int dim = 1; dim < 4; dim++) {
    it_gmi = gmi_begin(mdl, dim);
    while ((e_gmi = gmi_next(mdl, it_gmi))) {
      int c_tag = gmi_tag(mdl,e_gmi);
      rat_map.at(dim - 1)[c_tag] = len_map.at(dim - 1)[c_tag]/len;
    }
    gmi_end(mdl, it_gmi);
  }
}

// Set the rat_map for boundary cells such that they are refined to yield edge 
// length len.
void ModelEdgeRefinerVarying::set_bound_cell(double len,
                                std::vector<std::pair<int,int> >* cells) {
  if(cells != NULL) {
    for(int i = 0; i < cells->size(); i++) {
      int dim = cells->at(i).first;
      int c_tag = cells->at(i).second;
      rat_map.at(dim - 1)[c_tag] = len_map.at(dim - 1)[c_tag]/len;
    }
  }
  else {
    gmi_model* mdl = mesh->getModel();
    struct gmi_iter* it_gmi;
    struct gmi_ent* e_gmi;
    struct gmi_set* s_gmi;

    for(int dim = 1; dim < 3; dim++) {
      it_gmi = gmi_begin(mdl, dim);
      while ((e_gmi = gmi_next(mdl, it_gmi))) {
        int c_tag = gmi_tag(mdl,e_gmi);
        rat_map.at(dim - 1)[c_tag] = len_map.at(dim - 1)[c_tag]/len;
      }
      gmi_end(mdl, it_gmi);
    }
  }
}

double ModelEdgeRefinerVarying::measure_dist(ma::Entity* edge) {

  // TODO Use integrators of apf to integrate the field directly. For linear
  // elements and for isotropic metrics, it is equivalent. 
  apf::MeshElement* me = apf::createMeshElement(mesh, edge);
  double len = measure(edge);
  apf::destroyMeshElement(me);

  apf::ModelEntity* mdl = mesh->toModel(edge);
  int dim_e = mesh->getModelType(mdl);
  int tag_e = mesh->getModelTag(mdl);

  if(len_map.at(dim_e - 1)[tag_e] < std::numeric_limits<double>::min())
    return len;
  return len/len_map.at(dim_e - 1)[tag_e]*rat_map.at(dim_e - 1)[tag_e];
}

bool ModelEdgeRefinerVarying::shouldCollapse(ma::Entity* edge) {
  return this->measure_dist(edge) < coarse_target;
}

bool ModelEdgeRefinerVarying::shouldSplit(ma::Entity* edge) {
  apf::Downward d_v;
  mesh->getDownward(edge, 0, d_v);
  int dim_e = mesh->getModelType(mesh->toModel(edge));
  int dim_v0 = mesh->getModelType(mesh->toModel(d_v[0]));
  int dim_v1 = mesh->getModelType(mesh->toModel(d_v[1]));
  //double v0 = apf::getScalar(field_step, d_v[0], 0);
  //double v1 = apf::getScalar(field_step, d_v[1], 0);

  //double len = vd_meas_ent(mesh, edge);
  //return (1 > v0) or (1 > v1);
  //return (2*len > v0) or (2*len > v1);
  //return (len > v0/4+v1/4);
  if(dim_e > dim_v0 and dim_e > dim_v1) {
    if(dim_e == 1)
      return true;
  }
  return this->measure_dist(edge) > split_target;
}

void ModelEdgeRefinerVarying::set_target_th(double coarse_in, double split_in) {
  split_target = split_in;
  coarse_target = coarse_in;
}

ModelEdgeRefinerVarying::~ModelEdgeRefinerVarying() {
  for(int i = 0; i < 3; i++) {
    len_map.at(i).clear();
  }
  len_map.clear();
}

////////////////////////////
// STEP
////////////////////////////

Step::Step(ma::Mesh* m, double sz) {
  if(sz > std::numeric_limits<double>::min())
    sizing = sz;
  else 
    sizing = 1;
  mesh = m;
  average = ma::getAverageEdgeLength(m);
  field_step = m->findField("adapt_step");
  assert(field_step);
}
double Step::getValue(ma::Entity* v) {
  double dist = apf::getScalar(field_step, v, 0);
  assert(dist > std::numeric_limits<double>::min());
  //return average*dist*sizing;
  return dist/sizing;
}

Step_ns::Step_ns(ma::Mesh* m, double sz) {
  if(sz > std::numeric_limits<double>::min())
    sizing = sz;
  else 
    sizing = 1;
  mesh = m;
  average = ma::getAverageEdgeLength(m);
  field_step = m->findField("adapt_step");
  assert(field_step);
}

double Step_ns::getValue(ma::Entity* v) {
  double dist = apf::getScalar(field_step, v, 0);
  assert(dist > std::numeric_limits<double>::min());
  //return average*dist*sizing;
  return dist;
}

void vd_adapt_0cell(apf::Mesh2* m, struct cell_base* c) {
  std::vector<double> avg_cell(0);
  avg_cell = upd_cell_rad(m);
  Linear_0cell sf(m, c, &avg_cell);
  ma::Input* in = ma::configure(m, &sf);
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;
  ma::adapt(in);
}

void vd_adapt(apf::Mesh2* m) {
  Linear sf(m, 1.1);
  ma::Input* in = ma::configure(m, &sf);
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;
  ma::adapt(in);
}

Linear_0cell::Linear_0cell(ma::Mesh* m, struct cell_base* c,
                                          std::vector<double>* avg_cell_in) :
  mesh(NULL), average(0), avg_cell(0), cell0(0), 
  cell0_pos(0, apf::Vector3(0, 0, 0)),
  cnc0(0, std::vector<std::vector<int > > (0, std::vector<int >(0) ) ) {
  avg_cell = *avg_cell_in;

  mesh = m;
  c_base = c;

  refresh();
}
void Linear_0cell::reload(ma::Mesh* m, struct cell_base* c, 
                                          std::vector<double>* avg_cell_in) {
  avg_cell = *avg_cell_in;

  mesh = m;
  c_base = c;
  refresh();
}

double Linear_0cell::getValue(ma::Entity* v) {
  return getDistance(v);
  //return getDistance(v)*sizing;
}

double Linear_0cell::getDistance(ma::Entity* v) {
  double distance = average+1;
  double temp_len;

  apf::Vector3 e_pos;
  apf::Vector3 temp;

  mesh->getPoint(v, 0, e_pos);

  int ent_type = mesh->getType(v);
  int d = mesh->typeDimension[ent_type];

  apf::ModelEntity* mdl = mesh->toModel(v);
  int c_type = mesh->getModelType(mdl);
  int c_tag = mesh->getModelTag(mdl);

  //std::cout << c_type << "c" << c_tag << " "
  //          << d << "-ent" << v << std::endl;
  for(int i = 0; i < cnc0.at(c_type).at(c_tag-1).size(); i++) {
    int c_id = cnc0.at(c_type).at(c_tag-1).at(i);
    if(cell0.at(c_id) != NULL) {
      //std::cout << "v " << cell0.at(c_id)
      //          << " pos " << cell0_pos.at(c_id) << std::endl;

      temp = e_pos - cell0_pos.at(c_id);
      temp_len = temp.getLength();
      if(distance > temp_len);
        distance = temp_len;
    }
  }
  if(distance < average) {

    //if(c_type == 0)
      return distance*0.7;
    //else
    //  return distance;
    //return distance*3;
    //return average*3;
    //return average*exp(-average/distance);
  }
  else
    //return avg_cell.at(c_type-1);
    return avg_cell.at(0)*0.7;
    //return average;
}

Linear_0cell::~Linear_0cell() {
  clear();
}


void Linear_0cell::refresh() {
  clear();
  load_0cell();
  average = ma::getAverageEdgeLength(mesh);
}

void Linear_0cell::load_0cell_ent() {
  cell0_pos.resize(c_base->get_sz(0));
  cell0.resize(c_base->get_sz(0));
  double z[3] = {0,0,0};

  for(int i = 0; i < cell0_pos.size(); i++) {
    cell0_pos.at(i).fromArray(z);
  }

  for(int i = 0; i < cell0.size(); i++) {
    cell0.at(i) = NULL;
  }

  apf::MeshIterator* it_e = mesh->begin(0);
  //std::cout << "Iterator." << std::endl;
  apf::MeshEntity* v;

  while (v = mesh->iterate(it_e)) {
    apf::ModelEntity* mdl = mesh->toModel(v);
    int type = mesh->getModelType(mdl);
    int tag = mesh->getModelTag(mdl);

    if(type == 0) {
      assert(cell0.at(tag-1) == NULL);
      cell0.at(tag-1) = v;
    }

  }
  mesh->end(it_e);

}

void Linear_0cell::load_0cell_pos() {
  apf::Vector3 vert_pos;
  for(int i = 0; i < cell0.size(); i++) {
    if(cell0.at(i) != NULL) {
      mesh->getPoint(cell0.at(i), 0, vert_pos);
      cell0_pos.at(i) = vert_pos;
    }
  }

}

// TODO very ugly notation.
void Linear_0cell::load_0cell() {
  load_0cell_ent();
  load_0cell_pos();

  cnc0.at(3).resize(c_base->get_sz(3));
  cnc0.at(2).resize(c_base->get_sz(2));
  cnc0.at(1).resize(c_base->get_sz(1));
  cnc0.at(0).resize(c_base->get_sz(0));

  struct ent_conn e_cover;
  int max_sz = 1;
  for(int i = 0; i < cnc0.at(3).size(); i++) {

    c_base->get_conn_lower(0, 3, i, &e_cover);
    cnc0.at(3).at(i).reserve(e_cover.conn.size());
    for(int j = 0; j < e_cover.conn.size(); j++) {
      cnc0.at(3).at(i).push_back(e_cover.conn.at(j));
    }
    if(max_sz < e_cover.conn.size())
      max_sz = e_cover.conn.size();
  }

  for(int i = 0; i < cnc0.at(2).size(); i++) {
    cnc0.at(2).at(i).reserve(2*max_sz);
  }

  for(int i = 0; i < cnc0.at(1).size(); i++) {
    cnc0.at(1).at(i).reserve(7*max_sz);
  }

  for(int i = 0; i < cnc0.at(0).size(); i++) {
    c_base->get_conn_dim(0, 0, i, &e_cover);
    cnc0.at(0).at(i).reserve(e_cover.conn.size());
    for(int j = 0; j < e_cover.conn.size(); j++) {
      if(i != e_cover.conn.at(j))
        cnc0.at(0).at(i).push_back(e_cover.conn.at(j));
    }
  }

  // Go over 2cell adj. of 3cell, check for 0cell adj of the 3cell in the 0cell
  // list of the 2cell.
  for(int i = 0; i < cnc0.at(3).size(); i++) {
    //std::cout << "3c-" << i + 1 << std::endl;
    c_base->get_conn_lower(2, 3, i, &e_cover);
    for(int j = 0; j < e_cover.conn.size(); j++) {
      //std::cout << "2c-" << e_cover.conn.at(j) + 1 << std::endl;
      for(int k = 0; k < cnc0.at(3).at(i).size(); k++) {
        std::vector<int>::iterator c_st;
        std::vector<int>::iterator c_end;
        c_st = cnc0.at(2).at(e_cover.conn.at(j)).begin();
        c_end = cnc0.at(2).at(e_cover.conn.at(j)).end();
        int c0 = cnc0.at(3).at(i).at(k);
        //std::cout << "0c-" << c0 + 1 << ", " 
        //          << (std::find(c_st, c_end, c0) != c_end) << ", ";
        if(std::find(c_st, c_end, c0) == c_end) {
          cnc0.at(2).at(e_cover.conn.at(j)).push_back(c0);
        }
      }
      //std::cout << std::endl;
    }
    //std::cout << std::endl;

    c_base->get_conn_lower(1, 3, i, &e_cover);
    for(int j = 0; j < e_cover.conn.size(); j++) {
      for(int k = 0; k < cnc0.at(3).at(i).size(); k++) {
        std::vector<int>::iterator c_st;
        std::vector<int>::iterator c_end;
        c_st = cnc0.at(1).at(e_cover.conn.at(j)).begin();
        c_end = cnc0.at(1).at(e_cover.conn.at(j)).end();
        int c0 = cnc0.at(3).at(i).at(k);
        if(std::find(c_st, c_end, c0) == c_end) {
          cnc0.at(1).at(e_cover.conn.at(j)).push_back(c0);
        }
      }
    }
  }

  for(int i = 0; i < cnc0.size(); i++) {
    for(int j = 0; j < cnc0.at(i).size(); j++) {
      //std::cout << i << "c" << j+1 << std::endl;
      for(int k = 0; k < cnc0.at(i).at(j).size(); k++) {
        //std::cout << "0" << "c" << cnc0.at(i).at(j).at(k) + 1 << ", ";
      }
      //std::cout << std::endl;
    }
  }
}

void Linear_0cell::clear() {

  cell0_pos.clear();
  for(int i = 0; i < cnc0.size(); i++) {
    for(int j = 0; j < cnc0.at(i).size(); j++) {
      cnc0.at(i).at(j).clear();
    }
    cnc0.at(i).clear();
  }
  cnc0.clear();
  cnc0.resize(4);
}
