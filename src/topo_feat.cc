#include <cstdlib>
#include <algorithm>
#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_entlist.h"

#include "topo_geom.h"

#include "topo_write.h"
#include "topo_disp.h"

#include "topo_feat.h"

int find(std::vector<c_p>* cell_list, int id) {
  for(int i = 0; i < cell_list->size(); i++) {
    if(cell_list->at(i).second == id)
      return i;
  }
  return -1;
}

void vd_calc_field_dummy(apf::Mesh2* m, apf::MeshEntity* e,
                                                 vd_ext_opts* opts) {
}

void write_data(vd_series_ss* data, const char* fname) {
  csvfile csv(fname);
  csv << data->x_name << data->y_name
      << endrow;

  if(data->x.size() > 0) {
    for(int i = 0; i < data->x.size(); i++) {
      csv << data->x.at(i);
      csv << data->y.at(i);
      csv << endrow;
    }
  }

}

// Used in conjunction with write_label_y_vert
// Print together
// title\xlabel | xvalues
// step1        | yvalues1
// step2        | yvalues2
void write_whitespace_x_vert(vd_series_ss* data, const char* fname) {
  csvfile csv(fname);

  std::stringstream ss;
  ss << data->title << "\\" << data->x_name;

  csv << ss.str();
  if(data->x.size() > 0) {
    for(int i = 0; i < data->x.size(); i++) {
      csv << data->x.at(i);
    }
  }
  csv << endrow;
}

void write_label_y_vert(vd_series_ss* data, const char* fname) {
  csvfile csv(fname);

  if(data->y.size() > 0) {
    csv << data->variable;
    for(int i = 0; i < data->y.size(); i++) {
      csv << data->y.at(i);
    }
    csv << endrow;
  }
}

// TODO use templates, less ugly and probably easier to maintain
// vd_series
void vd_series_ss::clear() {
  dummy_clear_stop();
  x_name.clear();
  y_name.clear();
  title.clear();
  variable.clear();
  x.clear();
  y.clear();
}

vd_series_ss::vd_series_ss(const vd_series_ss& that) {
  dummy_clear_stop();

  x_name = that.x_name;
  y_name = that.y_name;
  title.clear();
  variable.clear();
  x = that.x;
  y = that.y;
}
vd_series_ss& vd_series_ss::operator=(const vd_series_ss& that) {
  dummy_clear_stop();
  x_name = that.x_name;
  y_name = that.y_name;
  title.clear();
  variable.clear();
  x = that.x;
  y = that.y;
  return *this;
}

vd_series_ss::vd_series_ss() {
  x_name.reserve(0);
  y_name.reserve(0);
  title.reserve(0);
  variable.reserve(0);
  x.reserve(0);
  y.reserve(0);
}

vd_series_ss::~vd_series_ss() {
  clear();
}

void vd_series_sv::clear() {
  dummy_clear_stop();
  x_name.clear();
  y_name.clear();
  title.clear();
  variable.clear();
  x.clear();
  y.clear();
}

vd_series_sv::vd_series_sv(const vd_series_sv& that) {
  dummy_clear_stop();
  x_name = that.x_name;
  y_name = that.y_name;
  title.clear();
  variable.clear();
  x = that.x;
  y = that.y;
}

vd_series_sv& vd_series_sv::operator=(const vd_series_sv& that) {
  dummy_clear_stop();
  x_name = that.x_name;
  y_name = that.y_name;
  title = that.title;
  variable = that.variable;
  x = that.x;
  y = that.y;
  return *this;
}

vd_series_sv::vd_series_sv() {
  x_name.reserve(0);
  y_name.reserve(0);
  title.reserve(0);
  variable.reserve(0);
  x.reserve(0);
  y.reserve(0);
}

vd_series_sv::~vd_series_sv() {
  clear();
}

void vd_series_vs::clear() {
  dummy_clear_stop();
  x_name.clear();
  y_name.clear();
  title.clear();
  variable.clear();
  x.clear();
  y.clear();
}

vd_series_vs::vd_series_vs(const vd_series_vs& that) {
  x_name = that.x_name;
  y_name = that.y_name;
  title = that.title;
  variable = that.variable;
  x = that.x;
  y = that.y;
}
vd_series_vs& vd_series_vs::operator=(const vd_series_vs& that) {
  x_name = that.x_name;
  y_name = that.y_name;
  title = that.title;
  variable = that.variable;
  x = that.x;
  y = that.y;
  return *this;
}

vd_series_vs::vd_series_vs() {
  x_name.reserve(0);
  y_name.reserve(0);
  title.reserve(0);
  variable.reserve(0);
  x.reserve(0);
  y.reserve(0);
}

vd_series_vs::~vd_series_vs() {
  clear();
}

void vd_series_vv::clear() {
  dummy_clear_stop();
  x_name.clear();
  y_name.clear();
  title.clear();
  variable.clear();
  x.clear();
  y.clear();
}

vd_series_vv::vd_series_vv(const vd_series_vv& that) {
  x_name = that.x_name;
  y_name = that.y_name;
  title = that.title;
  variable = that.variable;
  x = that.x;
  y = that.y;
}
vd_series_vv& vd_series_vv::operator=(const vd_series_vv& that) {
  x_name = that.x_name;
  y_name = that.y_name;
  title = that.title;
  variable = that.variable;
  x = that.x;
  y = that.y;
  return *this;
}

vd_series_vv::vd_series_vv() {
  x_name.reserve(0);
  y_name.reserve(0);
  title.reserve(0);
  variable.reserve(0);
  x.reserve(0);
  y.reserve(0);
}

vd_series_vv::~vd_series_vv() {
  clear();
}

// vd_timer
void vd_timer::set_time(int slot) {
  begin.at(slot) = PCU_Time();
  //clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
}

// Get time elapsed for timer slot.
double vd_timer::get_time(int slot) {
  //clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
  //long sec = end.tv_sec - begin.tv_sec;
  //long nsec = end.tv_nsec - begin.tv_nsec;
  //elapsed = sec + nsec*1e-9;
  end.at(slot) = PCU_Time();
  elapsed.at(slot) = end.at(slot) - begin.at(slot);
  return elapsed.at(slot);
}

// Get time elapsed for timer slot.
double vd_timer::add_time(int slot) {
  //clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
  //long sec = end.tv_sec - begin.tv_sec;
  //long nsec = end.tv_nsec - begin.tv_nsec;
  //elapsed = sec + nsec*1e-9;
  end.at(slot) = PCU_Time();
  elapsed.at(slot) = elapsed.at(slot) + end.at(slot) - begin.at(slot);
  return elapsed.at(slot);
}

// Return the last recorded elapsed for timer slot.
double vd_timer::get_elapsed(int slot) {
  return elapsed.at(slot);
}

double vd_timer::get_elapsed_sum() {
  double sum = 0.;
  for(int i = 0; i < elapsed.size(); i++)
    sum = sum + elapsed.at(i);
  return sum;
}

void vd_timer::resize(int slots) {
  elapsed.resize(slots);
  end.resize(slots);
  begin.resize(slots);
}

void vd_timer::reset_elapsed() {
  for(int i = 0; i < elapsed.size(); i++)
    elapsed.at(i) = 0;
}

/*
void vd_timer::write2csv(csvfile>& csv_in) {
  for(int i = 0; i < elapsed.size(); i++)
    csv_in << elapsed.at(i);
  csv_in << endrow;
}
*/

vd_timer::vd_timer() : elapsed(1, 0.), begin(1, 0.), end(1, 0.) {
  set_time(0);
}

//vd_timer::vd_timer(int TIME_TYPE_IN) : elapsed(0.) {
vd_timer::vd_timer(int slots) : elapsed(1, 0.), begin(1, 0.), end(1, 0.) {
/*
  TIME_TYPE = TIME_TYPE_IN;

  if(TIME_TYPE == CLOCK_REALTIME or
     TIME_TYPE == CLOCK_MONOTONIC or 
     TIME_TYPE == CLOCK_PROCESS_CPUTIME_ID or 
     TIME_TYPE == CLOCK_THREAD_CPUTIME_ID or 
     TIME_TYPE == CLOCK_MONOTONIC_RAW or 
     TIME_TYPE == CLOCK_REALTIME_COARSE or 
     TIME_TYPE == CLOCK_MONOTONIC_COARSE or 
     TIME_TYPE == CLOCK_BOOTTIME or 
     TIME_TYPE == CLOCK_REALTIME_ALARM or 
     TIME_TYPE == CLOCK_BOOTTIME_ALARM or 
     TIME_TYPE == CLOCK_TAI) {
  }
  else {
    TIME_TYPE = CLOCK_PROCESS_CPUTIME_ID;
  }
*/
  elapsed.resize(slots);
  end.resize(slots);
  begin.resize(slots);
  for(int i = 0; i < elapsed.size(); i++)
    set_time(i);
}
vd_timer::~vd_timer() {
}

// vd_ext_opts
void vd_ext_opts::vd_calc_field(apf::Mesh2* m, apf::MeshEntity* e) {
  if(!vd_calc_field_pt == 0)
    vd_calc_field_pt(m, e, this);
}

void vd_ext_opts::clear() {
  dummy_clear_stop();

  n_name.clear();
  n_flag = false;
  crt_flag = false;
  f_flag = false;
  vd_calc_field_pt = NULL;

  c_top.clear();
  for(int i = 0; i < c_bound.size(); i++)
    c_bound.at(i).clear();
  c_bound.clear();
}

void vd_ext_opts::set_func(void (*func_pt_in)(apf::Mesh2*, apf::MeshEntity*, vd_ext_opts*), char* n_in, int f_type_in, bool crt_in) {
  if(func_pt_in != NULL and n_in != NULL) {
    vd_calc_field_pt = func_pt_in;
    n_name = n_in;
    n_flag = true;
    f_type = f_type_in;
    f_flag = true;
    crt_flag = crt_in;
  }
}

void vd_ext_opts::set_search(VD_SEARCH s_type, apf::Vector3 dir_in) {
  search_type = (int) s_type;
  search_dir = dir_in;
}

void vd_ext_opts::set_extract_cells(const std::vector<c_p>* c_top_in,
                      const std::vector<std::vector<c_p> >* c_bound_in) {

  c_top.clear();
  for(int i = 0; i < c_bound.size(); i++)
    c_bound.at(i).clear();
  c_bound.clear();

  c_top = *c_top_in;
  c_bound.resize(c_bound_in->size());
  for(int i = 0; i < c_bound_in->size(); i++)
    c_bound.at(i) = c_bound_in->at(i);
}

vd_ext_opts::vd_ext_opts() {
  n_name.reserve(0);
  file_name.reserve(0);

  c_top.reserve(0);
  c_bound.reserve(0);

  vd_calc_field_pt = NULL;
  n_name = "";
  file_name = "./output/canvas.csv";
  n_flag = false;
  f_type = -1;
  f_flag = false;
  crt_flag = false;

  search_type = (int) VD_SEARCH::AVERAGE;
  double z[3] = {0,0,0};   
  search_dir.fromArray(z);
}

vd_ext_opts::vd_ext_opts(const vd_ext_opts& that) {
  n_name.reserve(0);
  file_name.reserve(0);

  c_top.reserve(0);
  c_bound.reserve(0);

  n_name = that.n_name;
  file_name = that.file_name;
  n_flag = that.n_flag;
  crt_flag = that.crt_flag;
  f_flag = that.f_flag;
  f_type = that.f_type;

  vd_calc_field_pt = that.vd_calc_field_pt;

  search_type = that.search_type;
  search_dir = that.search_dir;

  set_extract_cells(&that.c_top, &that.c_bound); 
}

vd_ext_opts& vd_ext_opts::operator=(const vd_ext_opts& that) {
  n_name = that.n_name;
  file_name = that.file_name;
  n_flag = that.n_flag;
  crt_flag = that.crt_flag;
  f_flag = that.f_flag;
  f_type = that.f_type;

  vd_calc_field_pt = that.vd_calc_field_pt;

  search_type = that.search_type;
  search_dir = that.search_dir;

  set_extract_cells(&that.c_top, &that.c_bound);
  return *this;
}

vd_ext_opts::vd_ext_opts(char* field_name, bool crt_in) {
  n_name.reserve(0);
  file_name.reserve(0);

  c_top.reserve(0);
  c_bound.reserve(0);

  n_name = field_name;
  file_name = "./output/canvas.csv";
  n_flag = true;
  crt_flag = crt_in;
  f_type = -1;

  first_ext = true;

  f_flag = false;
  vd_calc_field_pt = NULL;

  search_type = (int) VD_SEARCH::AVERAGE;
  double z[3] = {0,0,0};   
  search_dir.fromArray(z);
}

vd_ext_opts::~vd_ext_opts() {
  clear();
}

void ext_mesh::set_opts(vd_ext_opts opt_in) {
  m_opts = opt_in;
}

bool ext_mesh::chk_field() {
  assert(m_opts.n_flag);
  assert(m_opts.f_type == apf::SCALAR or m_opts.f_type == apf::VECTOR);
  return m->findField(m_opts.n_name.c_str());
}

bool ext_mesh::create_field() {
  if(!chk_field()) {
    if(m_opts.f_type == apf::SCALAR) {
      vd_att_vs_field(m, m_opts.n_name.c_str());
    }
    else if(m_opts.f_type == apf::VECTOR) {
      vd_att_vv_field(m, m_opts.n_name.c_str());
    }
  }
}

void ext_mesh::set_mesh(apf::Mesh2* m_in, cell_base* cb_in) {
  m = m_in;
  cb = cb_in;
}

// Based on the options, process the field information related to the ents.
// Length, area, volume weighted average.
double ext_mesh::get_avg_db(vd_entlist* e_list, int dim, int c_id) {
  double meas_tot = 0;
  double average = 0;

  std::map<apf::MeshEntity*, double> weight;
  std::vector<apf::MeshEntity*>* e_pt = &e_list->e.at(dim).at(c_id-1).at(0);
  for(int i = 0; i < e_pt->size(); i++) {
    weight[e_pt->at(i)] = 0;
  }

  e_pt = &e_list->e.at(dim).at(c_id-1).at(dim);

  apf::Downward d_v;
  apf::ModelEntity* mdl = m->findModelEntity(dim, c_id);

  if(e_pt->size() == 1) {
    return 0;
  }
  else {
    for(int i = 0; i < e_pt->size(); i++) {
      double meas_curr = vd_meas_ent(m, e_pt->at(i));
      m->getDownward(e_pt->at(i), 0, d_v);
      int in_cell = 0;
      for(int j = 0; j < dim+1; j++) {
        if(mdl == m->toModel(d_v[j]) ) {
          in_cell = in_cell + 1;
        }
      }
      assert(in_cell > 0);
      for(int j = 0; j < dim+1; j++) {
        if(mdl == m->toModel(d_v[j]) ) {
          weight[d_v[j]] = weight[d_v[j]] + meas_curr/in_cell;
        }
      }
      meas_tot = meas_tot + meas_curr;
    }

    apf::Field* field_curr = m->findField(m_opts.n_name.c_str());
    e_pt = &e_list->e.at(dim).at(c_id-1).at(0);
    for(int i = 0; i < e_pt->size(); i++) {
      average = average + weight[e_pt->at(i)]*
                                apf::getScalar(field_curr, e_pt->at(i), 0);
    }
    return average/meas_tot;
  }
}
/*
apf::Vector3 ext_mesh::get_avg_v(vd_entlist* e_list, int dim, int c_id) {

}

double ext_mesh::get_int_db(std::vector<apf::MeshEntity*>* ents) {

}
apf::Vector3 ext_mesh::get_int_v(std::vector<apf::MeshEntity*>* ents) {

}


void ext_mesh::get_per_cell_dim(std::vector<vd_series_ss>* data, int dim) {
}
void ext_mesh::calc_total_curv(vd_series_ss* data, int dim) {
}

void ext_mesh::get_per_cell_all(std::vector<vd_series_ss>* data) {
}
*/

void ext_mesh::get_meas_cell_all(std::vector<vd_series_ss>* data) {
  dummy_clear_stop();

  for(int i = 0; i < data->size(); i++) {
    data->at(i).clear();
  }
  data->clear();

  data->resize(3);
  for(int i = 0; i < 3; i++) {
    data->at(i).x.resize(cb->get_sz(i+1));
    data->at(i).y.resize(cb->get_sz(i+1));
    std::stringstream ss;
    ss << i+1 << "cell";
    data->at(i).x_name = ss.str();
    for(int j = 0; j < data->at(i).x.size(); j++) {
      data->at(i).x.at(j) = j + 1;
    }
    if(i == 0) {
      data->at(i).y_name = "length";
    }
    else if(i == 1) {
      data->at(i).y_name = "area";
    }
    else {
      data->at(i).y_name = "volume";
    }
  }

  for(int dim = 1; dim < 4; dim++) {
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(dim);

    //std::cout << dim << "-ents:" << std::endl; 
    while ((e = m->iterate(it))) {
      apf::ModelEntity* mdl = m->toModel(e);
      int type = m->getModelType(mdl);
      int tag = m->getModelTag(mdl);
      //std::cout << type << "c" << tag << " " << e; 
      if(type == dim) {
        double meas = vd_meas_ent(m, e);
        //std::cout << " " << meas;
        data->at(dim-1).y.at(tag-1) = data->at(dim-1).y.at(tag-1) + meas;
      }
      //std::cout<< std::endl;
    }
    m->end(it);
  }
}

void ext_mesh::get_cell_meas(vd_series_ss* data, std::vector<c_p>* c_top) {
  create_field();
  data->clear();
  if (c_top->size() > 0) {
    std::stringstream ss;
    data->x_name = ss.str();
    data->y_name = m_opts.n_name;

    data->x.resize(c_top->size());
    data->y.resize(c_top->size());

    int dim = c_top->at(0).first; 
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(dim);

    std::cout << dim << "-ents:" << std::endl; 
    while ((e = m->iterate(it))) {
      apf::ModelEntity* mdl = m->toModel(e);
      int type = m->getModelType(mdl);
      int tag = m->getModelTag(mdl);
      std::cout << type << "c" << tag << " " << e; 

      int i1 = find(c_top, tag);
      if(type == dim and i1 > -1) {
        double meas = vd_meas_ent(m, e);
        std::cout << " " << meas;
        data->y.at(i1) = data->y.at(i1) + meas;
      }
      std::cout<< std::endl;
    }
    m->end(it);
  }
}

//void ext_mesh::get_cell_curv(vd_series_ss* data, std::vector<c_p>* c_top, std::vector<std::vector<c_p> > * c_bound) {
//vd_mean_curv(apf::Mesh2* m, apf::MeshEntity* v, int cell)
//}

void ext_mesh::get_cell_bound(vd_series_sv* data, int dim, int c_id,
                                          std::vector<c_p>* c_bound) {
  create_field();
  data->clear();
  if (c_bound->size() > 0) {
    std::stringstream ss;
    ss << 3 << "cell" << c_id;
    data->x_name = ss.str();
    data->y_name = m_opts.n_name;

    int d_b = c_bound->at(0).first;
    ent_conn e_dim;
    cb->get_conn_dim_gmi(d_b, dim, c_id, &e_dim);
    for(int i = 0; i < c_bound->size(); i++) {
      assert(e_dim.find_ent(c_bound->at(i).second) > -1);
    }

    vd_entlist e_list(m, cb);
    for(int i = 0; i < c_bound->size(); i++) {
      std::vector<apf::MeshEntity*>* e_pt = 
              &e_list.e.at(d_b).at(c_bound->at(i).second-1).at(0);

      for(int j = 0; j < e_pt->size(); j++) {
        m_opts.vd_calc_field(m, e_pt->at(j));
      }
    }
    apf::Field* field_curr = m->findField(m_opts.n_name.c_str());

    for(int i = 0; i < c_bound->size(); i++) {
      if(m_opts.search_type == (int) VD_SEARCH::AVERAGE) {
      }
      else if(m_opts.search_type == (int) VD_SEARCH::TOTAL) {
      }
      else if(m_opts.search_type == (int) VD_SEARCH::MINIMUM) {
      }
      else if(m_opts.search_type == (int) VD_SEARCH::MAXIMUM) {
      }
      else if(m_opts.search_type == (int) VD_SEARCH::MAX_DIR) {
      }
      else {
        assert(m_opts.search_type == (int) VD_SEARCH::MIN_DIR);
      }
    }
  }
}

void ext_mesh::get_cell_bound_ext_angle(vd_series_ss* data, int dim, int c_id,
                                          std::vector<c_p>* c_bound) {
/*
  create_field();
  data->clear();
  if (c_bound->size() > 0) {
    std::stringstream ss;
    ss << 3 << "cell" << c_id;
    data->x_name = ss.str();
    data->y_name = m_opts.n_name;

    data->x.resize(c_bound->size());
    data->y.resize(c_bound->size());

    apf::Field* field_curr = m->findField(m_opts.n_name.c_str());

    int d_b = c_bound->at(0).first;
    ent_conn e_dim;
    cb->get_conn_dim_gmi(d_b, dim, c_id, &e_dim);
    for(int i = 0; i < c_bound->size(); i++) {
      assert(e_dim.find_ent(c_bound->at(i).second) > -1);
    }

    vd_entlist e_list(m, cb);
    for(int i = 0; i < c_bound->size(); i++) {
      std::vector<apf::MeshEntity*>* e_pt = 
              &e_list.e.at(d_b).at(c_bound->at(i).second-1).at(0);

      for(int j = 0; j < e_pt->size(); j++) {
        apf::setScalar(field_curr, e_pt->at(j), 0, 
                            vd_ext_angle(m, e_pt->at(j), c_id));
      }
    }

    for(int i = 0; i < c_bound->size(); i++) {
      if(m_opts.search_type == (int) VD_SEARCH::AVERAGE) {
        data->x.at(i) = c_bound->at(i).second;
        data->y.at(i) = get_avg_db(&e_list, 1, c_bound->at(i).second);
      }
      else if(m_opts.search_type == (int) VD_SEARCH::TOTAL) {
      }
      else if(m_opts.search_type == (int) VD_SEARCH::MINIMUM) {
      }
      else if(m_opts.search_type == (int) VD_SEARCH::MAXIMUM) {
      }
      else if(m_opts.search_type == (int) VD_SEARCH::MAX_DIR) {
      }
      else {
        assert(m_opts.search_type == (int) VD_SEARCH::MIN_DIR);
      }
    }
  }
*/
}

void ext_mesh::calc_field() {
  if(!m_opts.vd_calc_field_pt == 0) {
    assert(m_opts.f_type > -1);
    create_field();

    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(0);

    while ((e = m->iterate(it))) {
      m_opts.vd_calc_field(m, e);
    }
    m->end(it);
  }
}

void ext_mesh::clear() {
}
ext_mesh::ext_mesh() {
}
ext_mesh::~ext_mesh() {
}

double ext_MS::check_ext_angle(apf::MeshEntity* tet, apf::MeshEntity* tri1, apf::MeshEntity* tri2, apf::MeshEntity* edge) {
  apf::Vector3 pos_tet(0,0,0);
  apf::Vector3 pos_tri1(0,0,0);
  apf::Vector3 pos_tri2(0,0,0);
  apf::Vector3 dir_line(0,0,0);
  apf::Vector3 pos_line(0,0,0);

  pos_tet = vd_get_pos(m, tet);
  pos_tri1 = vd_get_pos(m, tri1);
  pos_tri2 = vd_get_pos(m, tri2);
  pos_line = vd_get_pos(m, edge);
  dir_line = get_edge_dir(m, edge);

  return vd_ext_angle(pos_tri1, pos_tri2, dir_line, pos_line, pos_tet);
}

void ext_MS::process_bound(int dim, int id) {
  std::vector<apf::MeshEntity*> tris(0);
  std::vector<apf::MeshEntity*> tets(0);

  int e_nbr = e_list.e.at(dim).at(id).at(1).size();
  std::map<apf::MeshEntity*, int> tri_tet1{};
  std::map<apf::MeshEntity*, int> tri_tet2{};
  std::map<apf::MeshEntity*, int> tet_tri1{};
  std::map<apf::MeshEntity*, int> tet_tri2{};
  apf::Up up;

  // This could be faster if done in a top down approach, by looping over 3cells
  // and counting for each 1 and 2 cell number of 3cell adjacencies.
  // That number will be used to define the equilibrium angles in MS relation
  // and can be used to 
  ent_conn e_3c;

  cb->get_conn_dim(3, dim, id, &e_3c);
  bool disjoint = false;
  apf::Vector3 pos_tet(0,0,0);
  apf::Vector3 pos_tri1(0,0,0);
  apf::Vector3 pos_tri2(0,0,0);
  apf::Vector3 norm1(0,0,0);
  apf::Vector3 norm2(0,0,0);
  apf::Vector3 pos_line(0,0,0);
  apf::Vector3 dir_line(0,0,0);
  double len = 0;

  for(int i = 0; i < e_nbr; i++) {
    apf::MeshEntity* e_curr = e_list.e.at(dim).at(id).at(1).at(i);

    dir_line = get_edge_dir(m, e_curr);
    pos_line = vd_get_pos(m, e_curr);
    len = vd_meas_ent(m, e_curr);

    vd_set_up(m, e_curr, &tris);
    vd_set_up(m, &tris, &tets);

    tri_tet1.clear();
    tri_tet2.clear();
    tet_tri1.clear();
    tet_tri2.clear();

    bool broken = false;
    int end1 = -1;
    int end2 = -1;
    int start = -1;

    for(int j = 0; j < tris.size(); j++) {

      apf::MeshEntity* tri_curr = tris.at(j);

      apf::ModelEntity* mdl_tri = m->toModel(tris.at(j));
      if(start == -1 and m->getModelType(mdl_tri) == 2) {
        start = j;
      }
      m->getUp(tri_curr, up);
      bool found_1 = false;
      bool found_2 = false;
      for (int k = 0; k < up.n; k++) {
        std::vector<apf::MeshEntity*>::iterator it_e;
        it_e = std::find(tets.begin(), 
                         tets.end(), up.e[k]);
        // If on the disc, replace it with the new edge for the top disc.
        if(it_e != tets.end()) {
          if(!found_1) {
            found_1 = true;
            tri_tet1[tri_curr] = std::distance(tets.begin(), it_e)+1;
            if(tet_tri1[*it_e] == 0) {
              tet_tri1[*it_e] = j+1;
            }
            else {
              assert(tet_tri2[*it_e] == 0);
              tet_tri2[*it_e] = j+1;
            }
          }
          else {
            assert(!found_2);
            found_2 = true;
            tri_tet2[tri_curr] = std::distance(tets.begin(), it_e)+1;
            if(tet_tri1[*it_e] == 0) {
              tet_tri1[*it_e] = j+1;
            }
            else {
              assert(tet_tri2[*it_e] == 0);
              tet_tri2[*it_e] = j+1;
            }
            k = up.n;
          }
        }
      }
      if(up.n == 1) {
        broken = true;
        if(end1 > -1) {
          assert(end2 == -1);
          end2 = j;
        }
        else {
          end1 = j;
        }
      }
    }
    assert(start > -1);

    if(broken) {
      start = end1;
    }

    // Starting from a 2stratum triangle, go over tets until a second 2stratum
    // triangle is reached. Add the bounded 3stratum and the bounding triangle
    // ids to the lists. Do this until all tets are scanned.
    std::vector<apf::MeshEntity*> c3s(0);
    std::vector<std::pair<int, int> > bounds(0, std::make_pair(-1,-1));
    bounds.reserve(tets.size());
    c3s.reserve(tets.size());

    int curr = start;
    apf::MeshEntity* tri_curr = tris.at(curr);
    int tet_id = tri_tet1[tri_curr] - 1;
    apf::MeshEntity* tet_curr = tets.at(tet_id);
    int next = tet_tri1[tet_curr] - 1;
    if(next == curr) {
      next = tet_tri2[tet_curr] - 1;
    }
    apf::MeshEntity* tet_next;

    apf::MeshEntity* tri_next = tris.at(next);
    apf::ModelEntity* mdl_tri2 = m->toModel(tri_next);

    int bound1 = curr;
    for(int i = 0; i < tets.size() - 1; i++) {
      if(m->getModelType(mdl_tri2) == 2) {
        c3s.push_back(tet_curr);
        bounds.push_back(std::make_pair(bound1, next));
        bound1 = next;
      }
      tet_id = tri_tet1[tri_next] - 1;
      tet_next = tets.at(tet_id);
      if(tet_next == tet_curr) {
        assert(tri_tet2[tri_next] != 0);
        tet_id = tri_tet2[tri_next] - 1;
        tet_next = tets.at(tet_id);
      }
      curr = next;
      next = tet_tri1[tet_next] - 1;
      if(next == curr) {
        next = tet_tri2[tet_next] - 1;
      }
      tri_next = tris.at(next);
      mdl_tri2 = m->toModel(tri_next);
      tet_curr = tet_next;
    }
    assert(m->getModelType(mdl_tri2) == 2);
    c3s.push_back(tet_curr);
    bounds.push_back(std::make_pair(bound1, next));
    bound1 = next;


    if(e_3c.conn.size() != c3s.size())
      disjoint = true;

    double ang_equi_ext;
    if(dim == 1) {
      // B: interior angle, a: exterior angle, n: Number of 3strata around 1strata
      // (180 - a)*n = B
      if(cb->get_1c_corner(id)) {
        ang_equi_ext = PI_L - PI_L/2/c3s.size();
      }
      // Exterior:
      else if(cb->get_cell_ext(dim, id)) {
        ang_equi_ext = PI_L - PI_L/c3s.size();
      }
      // Interior: 2PI/3 with 3 grains
      else {
        ang_equi_ext = PI_L - 2*PI_L/c3s.size();
      }
    }
    double ang_sum = 0;
    double dVdt_sum = 0;
    double angle_curr = 0;
    apf::MeshEntity* tri1;
    apf::MeshEntity* tri2;
    int c3_curr = 0;
    for(int j = 0; j < c3s.size(); j++) {
      c3_curr = m->getModelTag(m->toModel(c3s.at(j)));
      //pos_tet = vd_get_pos(m, c3s.at(j));
      tri1 = tris.at(bounds.at(j).first);
      tri2 = tris.at(bounds.at(j).second);
      pos_tri1 = vd_get_pos(m, tri1) - pos_line;
      pos_tri2 = vd_get_pos(m, tri2) - pos_line; 
      norm1 = norm_0(vd_area_out(m, tri1, c3_curr));
      norm2 = norm_0(vd_area_out(m, tri2, c3_curr));
      angle_curr = vd_ext_angle_n(pos_tri1, pos_tri2, norm1, norm2);

      ang_sum = ang_sum + PI_L - angle_curr;
      // Ln term * 2PI = len*angle_curr, per edge
      double Ln = len*angle_curr;
      dVdt_c3[c3_curr] = dVdt_c3[c3_curr] - Ln;
      dVdt_sum = dVdt_sum - Ln;
      if(dim == 1) {
        // Mn term * 2PI/6 = len*ang_equi_ext, per edge
        double Mn = len*ang_equi_ext;
        dVdt_c3[c3_curr] = dVdt_c3[c3_curr] + Mn;
        dVdt_sum = dVdt_sum + Mn;
      }
    }

    //if(std::fabs(dVdt_sum) > std::numeric_limits<double>::min())
    //  std::cout << "The total change of volume is " << dVdt_sum << "!" 
    //            << " angle deviation is " << ang_sum - 
    //                                         (PI_L - ang_equi_ext)*c3s.size()
    //            << std::endl;
  }
  //if(disjoint) {
  //  std::cout << "A 3stratum has disjoint components around " << dim << "c" 
  //            << id + 1 << std::endl;
  //}
}

void ext_MS::process_mesh() {
  e_list.change_mesh(m, cb);

  for(int i = 0; i < cb->get_sz(3); i++) {
    dVdt_c3[i+1] = 0;
  }
  for(int i = 0; i < cb->get_sz(2); i++) {
    if(!cb->is_free(2, i)) {
      process_bound(2, i);
    }
  }
  for(int i = 0; i < cb->get_sz(1); i++) {
    if(!cb->is_free(1, i)) {
      process_bound(1, i);
    }
  }
}

void ext_MS::clear() {
  len_c1.clear();
  dVdt_c3.clear();

  for(int i = 0; i < c1c3.size(); i++) {
    c1c3.at(i).clear();
  }
  c1c3.clear();
  for(int i = 0; i < dVdt_c1c3.size(); i++) {
    dVdt_c1c3.at(i).clear();
  }
  dVdt_c1c3.clear();

  for(int i = 0; i < c2c3.size(); i++) {
    c2c3.at(i).clear();
  }
  c2c3.clear();
  for(int i = 0; i < dVdt_c2c3.size(); i++) {
    dVdt_c2c3.at(i).clear();
  }
  dVdt_c2c3.clear();
}

ext_MS::ext_MS() : cb(NULL), m(NULL), e_list(), 
                   //proj_flag(PROJ_TYPE::FIXED), 
                   dVdt_c3{}, len_c1{},
                   c1c3(0, std::vector<int>(0)),
                   dVdt_c1c3(0, std::vector<double>(0)),
                   c2c3(0, std::vector<int>(0)),
                   dVdt_c2c3(0, std::vector<double>(0)) {
}

ext_MS::ext_MS(apf::Mesh2* m_in, cell_base* cb_in) : cb(cb_in), m(m_in), 
                   e_list(), 
                   //proj_flag(PROJ_TYPE::FIXED),
                   dVdt_c3{}, len_c1{},
                   c1c3(0, std::vector<int>(0)),
                   dVdt_c1c3(0, std::vector<double>(0)),
                   c2c3(0, std::vector<int>(0)),
                   dVdt_c2c3(0, std::vector<double>(0)) {
}

ext_MS::~ext_MS() {
  clear();
}

/////////////////////////////////
// ext_GAUSS
/////////////////////////////////
double ext_GAUSS::get_ang(apf::MeshEntity* t) {
  apf::Downward d_v;
  apf::Downward d_e;

  m->getDownward(t, 0, d_v);
  m->getDownward(t, 1, d_e);

  int lookup_tri_ed [3][2] = {{2,0},{0,1},{1,2}};

  int v1 = findIn(d_v, 3, v);
  int e1 = lookup_tri_ed[v1][0];
  int e2 = lookup_tri_ed[v1][1];

  double coef = 1.;
  apf::Vector3 temp(0,0,0);
  temp = norm_0(pos_map[d_e[e1]] - temp_pos);
  if(temp*ori_map[d_e[e1]] < -std::numeric_limits<double>::min())
    coef = -coef;
  temp = norm_0(pos_map[d_e[e2]] - temp_pos);
  if(temp*ori_map[d_e[e2]] < -std::numeric_limits<double>::min())
    coef = -coef;
  double ang_cos = ori_map[d_e[e2]]*ori_map[d_e[e1]]*coef;
  return std::acos(std::min(std::max(ang_cos,-1.0),1.0));
  //assert(ang_cos >= -1. and ang_cos <= 1.);
  //return std::acos(ang_cos);
}

void ext_GAUSS::collect_ang() {
  std::vector<apf::MeshEntity*> t_list(0);
  for(int c3_curr = 0; c3_curr < c3_nbr; c3_curr++) {
    std::vector<int> c2_bound(0);
    c2_bound = vd_3c->get_bound_ids(3, c3_curr);
    for(int i = 0; i < c2_bound.size(); i++) {
      if(!c2_map[c2_bound.at(i)]) {
        double ang_curr = 0.;
        t_list = vd_3c->get_dim_ent(2, c2_bound.at(i));
        for(int j = 0; j < t_list.size(); j++) {
          //ang_curr = ang_curr + get_ang(t_list.at(j));
          ang_curr = ang_curr + ang_map[t_list.at(j)][v];
          ang_map_count[t_list.at(j)][v] = ang_map_count[t_list.at(j)][v] + 1;
        }
        c2_map[c2_bound.at(i)] = true;
        c2_ang[c2_bound.at(i)] = ang_curr;
      }
    }
  }
}

void ext_GAUSS::load_v() {
  vd_3c->reload(m, cb, v);
  c3_nbr = vd_3c->get_3c_nbr();
  m->getPoint(v, 0, temp_pos);
  c2_map.clear();
  c2_ang.clear();
  collect_ang();
}

void ext_GAUSS::process_0c_vertex() {
  load_v();

  for(int c3_curr = 0; c3_curr < c3_nbr; c3_curr++) {
    apf::ModelEntity* mdl_3c = vd_3c->get_mdl(3, c3_curr);
    int c3_id = m->getModelTag(mdl_3c);
    std::vector<int> c2_bound(0);
    c2_bound = vd_3c->get_bound_ids(3, c3_curr);

    double ang_curr = 0.;
    for(int i = 0; i < c2_bound.size(); i++) {
      ang_curr = ang_curr + c2_ang[c2_bound.at(i)];
    }
    double ang_def_curr = 2*PI_L - ang_curr;
    c3_euler_mod.at(c3_id - 1) = c3_euler_mod.at(c3_id - 1) + ang_def_curr;
    c0_ang_def.at(c_id) = c0_ang_def.at(c_id) + ang_def_curr;

    if(ext_curr) {
      CoM = CoM + ang_def_curr;
//    if(ext_curr and 
//       ext_BC_type == VD_BC_TYPE::NEUMANN_EXT) {
//      CoM = CoM + ang_def_curr;
//    }
//    else if(ext_curr and 
//      (f_calc->get_proj() == (int)PROJ_TYPE::FIXED)) {
//     CoM = CoM + ang_def_curr;
      c0_ang_dev.at(c_id) = -1;
    }
    else {
      // TODO for exterior vertices, use the insertible 0-stratum detection
      // from topo_vd.cc.
      if(c3_nbr == 4) {
        ang_defect = ang_defect + ang_def_curr;
        ang_defect_const = ang_defect_const + ang_def_equi;
        ang_defect_dev = ang_defect_dev + (ang_def_curr - ang_def_equi);

        c0_ang_dev.at(c_id) = ang_def_curr - ang_def_equi;
      }
      else {
        ang_defect_trans = ang_defect_trans + ang_def_curr;
        c0_ang_dev.at(c_id) = -1;
      }
    }
  }
  apf::Field* vel_field = m->findField("velocity_field");
  apf::Vector3 v_temp(0,0,0);
  apf::getVector(vel_field, v, 0, v_temp);

  c0_vel.at(c_id) = v_temp.getLength();
}

void ext_GAUSS::process_ext_1c_vertex() {
  load_v();

  for(int c3_curr = 0; c3_curr < c3_nbr; c3_curr++) {
    apf::ModelEntity* mdl_3c = vd_3c->get_mdl(3, c3_curr);
    int c3_id = m->getModelTag(mdl_3c);
    std::vector<int> c2_bound(0);
    c2_bound = vd_3c->get_bound_ids(3, c3_curr);

    double ang_curr = 0.;
    for(int i = 0; i < c2_bound.size(); i++) {
      ang_curr = ang_curr + c2_ang[c2_bound.at(i)];
    }
    double ang_def_curr = 2*PI_L - ang_curr;
    c3_euler_mod.at(c3_id - 1) = c3_euler_mod.at(c3_id - 1) + ang_def_curr;
    CoM = CoM + ang_def_curr;
  }
}

void ext_GAUSS::process_int_1c_vertex() {
  load_v();

  for(int c3_curr = 0; c3_curr < c3_nbr; c3_curr++) {
    apf::ModelEntity* mdl_3c = vd_3c->get_mdl(3, c3_curr);
    int c3_id = m->getModelTag(mdl_3c);
    std::vector<int> c2_bound(0);
    c2_bound = vd_3c->get_bound_ids(3, c3_curr);

    double ang_curr = 0.;
    for(int i = 0; i < c2_bound.size(); i++) {
      ang_curr = ang_curr + c2_ang[c2_bound.at(i)];
    }
    double ang_def_curr = 2*PI_L - ang_curr;
    c3_euler_mod.at(c3_id - 1) = c3_euler_mod.at(c3_id - 1) + ang_def_curr;
    ang_defect_1c_int = ang_defect_1c_int + ang_def_curr;
  }
}

void ext_GAUSS::process_2c_vertex() {
  m->getPoint(v, 0, temp_pos);
  std::vector<apf::MeshEntity*> edges(0);
  std::vector<apf::MeshEntity*> tris(0);
  //std::map<apf::MeshEntity*, double> ang_map{};

  vd_set_up(m, v, &edges);
  vd_set_up(m, &edges, &tris);

  ent_conn e_conn;
  cb->get_conn_dim(3, 2, c_id, &e_conn);

  std::vector<apf::MeshEntity*> t_list(0);
  for(int i = 0; i  < e_conn.conn.size(); i++) {
    int c3_curr = e_conn.conn.at(i);

    double ang_curr = 0.;
    for(int j = 0; j < tris.size(); j++) {
      apf::ModelEntity* mdl_tri = m->toModel(tris.at(j));
      if(m->getModelType(mdl_tri) == 2) {
        //ang_curr = ang_curr + get_ang(t_list.at(j));
        ang_curr = ang_curr + ang_map[tris.at(j)][v];
        ang_map_count[tris.at(j)][v] = ang_map_count[tris.at(j)][v] + 1;
      }
    }
    double ang_def_curr = 2*PI_L - ang_curr;
    CoM = CoM + ang_def_curr;
    c3_euler_mod.at(c3_curr) = c3_euler_mod.at(c3_curr) + ang_def_curr;

    apf::setScalar(ang_def_field, v, 0, ang_def_curr);
  }
}

double ext_GAUSS::get_ang_val(apf::MeshEntity* tri, apf::MeshEntity* vert) {
  return ang_map[tri][vert];
}

int ext_GAUSS::get_ang_count(apf::MeshEntity* tri, apf::MeshEntity* vert) {
  return ang_map_count[tri][vert];
}

void ext_GAUSS::calc_gauss() {
  // Check the total of the gaussian curvatures by summing up the interior angles 
  // of the triangles bounding a 3-stratum and subtracting from 2*pi*n_v, where
  // n_v is the number of vertices bounding the 3-stratum.

  ang_def_field = vd_att_vs_field(m, "a_def");

  dim = 3;
  struct ent_conn* e_c = new ent_conn();
  for(c_id = 0; c_id < cb->get_sz(dim); c_id++) {
    if(!cb->is_free(dim, c_id)) {
      int v_nbr = 0;
      cb->get_conn_dim(0, dim, c_id, e_c);
      v_nbr = v_nbr + e_c->conn.size();

      cb->get_conn_dim(1, dim, c_id, e_c);
      for(int i = 0; i < e_c->conn.size(); i++) {
        v_nbr = v_nbr + e_list->e.at(1).at(e_c->conn.at(i)).at(0).size();
      }

      cb->get_conn_dim(2, dim, c_id, e_c);
      for(int i = 0; i < e_c->conn.size(); i++) {
        v_nbr = v_nbr + e_list->e.at(2).at(e_c->conn.at(i)).at(0).size();
      }

      apf::Downward d_v;
      double ang_glob = 0.;
      for(int i = 0; i < e_c->conn.size(); i++) {
        int c2_c = e_c->conn.at(i);
        std::vector<apf::MeshEntity*>* ents = &e_list->e.at(2).at(c2_c).at(2);
        for(int j = 0; j < ents->size(); j++) {
          m->getDownward(ents->at(j), 0, d_v);
          //std::cout << "Tri " << ents->at(i) << std::endl;
          for(int k = 0; k < 3; k++) {
            ang_glob = ang_glob + ang_map[ents->at(j)][d_v[k]];
          }
        }
      }
      std::cout << "3c" << c_id + 1 << " " << 2*v_nbr*PI_L - ang_glob 
                << std::endl;
    }
  }
  delete e_c;

  dim = 0;
  int sz = cb->get_sz(dim);
  bool has_ext = (f_calc->get_proj() == (int)PROJ_TYPE::FIXED);
  bool has_shell = (f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL);
  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      std::vector<apf::MeshEntity*>* verts = &e_list->e.at(dim).at(c_id).at(0);
      if(has_ext and cb->get_cell_ext(dim, c_id)) {
        ext_curr = true;
      }
      else if(has_shell and f_calc->get_e_sh()->chk_shell(dim, c_id) ) {
        ext_curr = true;
        shell sh = f_calc->get_e_sh()->get_shell(dim, c_id);
        if(sh.dim == 0) {
          ang_def_equi = ANG_EQUI_CONST_EXT_0C;
        }
        else if(sh.dim == 1) {
          ang_def_equi = ANG_EQUI_CONST_EXT_1C;
        }
        else if(sh.dim == 2)  {
          ang_def_equi = ANG_EQUI_CONST_EXT_2C;
        }
      }
      else {
        ext_curr = false;
        ang_def_equi = ANG_EQUI_CONST;
      }
      assert(verts->size() == 1);
      v = verts->at(0);
      process_0c_vertex();
    }
  }

  dim = 1;
  sz = cb->get_sz(dim);
  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      std::vector<apf::MeshEntity*>* verts = &e_list->e.at(dim).at(c_id).at(0);
      if(cb->get_cell_ext(dim, c_id)) {
        for(int i = 0; i < verts->size(); i++) {
          v = verts->at(i);
          process_ext_1c_vertex();
        }
      }
      else {
        for(int i = 0; i < verts->size(); i++) {
          v = verts->at(i);
          process_int_1c_vertex();
        }
      }
    }
  }

  dim = 2;
  sz = cb->get_sz(dim);
  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      std::vector<apf::MeshEntity*>* verts = &e_list->e.at(dim).at(c_id).at(0);
      //if(cb->get_cell_ext(dim, c_id)) {
        // TODO This is currently a placeholder. Would be used for free surface
        // BC that allows grooving. 
        // For Neumann BC, this is the Gaussian curv. is 0.
        //for(int j = 0; j < verts->size(); j++) {
        //  process_ext_2c_vertex(verts->at(j));
        //}
      //}
      //else {
        for(int i = 0; i < verts->size(); i++) {
          v = verts->at(i);
          process_2c_vertex();
        }
      //}
    }
  }

  apf::writeVtkFiles("./output/gauss", m);
  // Check the interior angle usage counts:

/*
  dim = 2;
  sz = cb->get_sz(dim);
  apf::Downward d_v;
  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      std::vector<apf::MeshEntity*>* tris = &e_list->e.at(dim).at(c_id).at(2);
      for(int i = 0; i < tris->size(); i++) {
        m->getDownward(tris->at(i), 0, d_v);
        for(int j = 0; j < 3; j++) {
          if(ang_map_count[tris->at(i)][d_v[j]] != 2) {
            apf::ModelEntity* mdl = m->toModel(tris->at(i));
            int cell_dim = m->getModelType(mdl);
            int cell_id = m->getModelTag(mdl);
            assert(cell_dim == 2 and cell_id == c_id + 1);
            std::cout << cell_dim << "c" << cell_id 
                      << " Tri " << tris->at(i) << " ";
            mdl = m->toModel(d_v[j]);
            cell_dim = m->getModelType(mdl);
            cell_id = m->getModelTag(mdl);
            std::cout << cell_dim << "c" << cell_id
                      << " v " << d_v[j] << " used "
                      << ang_map_count[tris->at(i)][d_v[j]] << std::endl;
          }
        }
      }
    }
  }
*/
}

void ext_GAUSS::collect_ori() {
  for(dim = 1; dim < 3; dim++) {
    int sz = cb->get_sz(dim);
    for(c_id = 0; c_id < sz; c_id++) {
      if(!cb->is_free(dim, c_id)) {
        std::vector<apf::MeshEntity*>* ents = &e_list->e.at(dim).at(c_id).at(1);
        for(int i = 0; i < ents->size(); i++) {
          ori_map[ents->at(i)] = get_edge_dir(m, ents->at(i));
          pos_map[ents->at(i)] = vd_get_pos(m, ents->at(i));
        }
      }
    }
  }

  apf::Downward d_v;
  dim = 2;
  int sz = cb->get_sz(dim);
  //std::cout << "Total interior angle deviation calculation." << std::endl;
  double ang_glob = 0.;
  int t_nbr = 0;
  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      std::vector<apf::MeshEntity*>* ents = &e_list->e.at(dim).at(c_id).at(2);
      for(int i = 0; i < ents->size(); i++) {
        m->getDownward(ents->at(i), 0, d_v);
        //std::cout << "Tri " << ents->at(i) << std::endl;
        double ang_tot = 0;
        for(int j = 0; j < 3; j++) {
          v = d_v[j];
          m->getPoint(v, 0, temp_pos);
          double ang_curr = get_ang(ents->at(i));
          //std::cout << "\t" << ang_curr << std::endl;
          ang_map[ents->at(i)][v] = ang_curr;
          ang_tot = ang_tot + ang_map[ents->at(i)][v];
        }
        //std::cout << "Tri " << ents->at(i) 
        //          << " " << PI_L - ang_tot << std::endl;
        ang_glob = ang_glob + ang_tot;
      }
      t_nbr = t_nbr + ents->size();
    }
  }
  //std::cout << "Total of interior angles " << ang_glob
  //          << " " << " n_tris*pi - ang_glob " << t_nbr*PI_L - ang_glob
  //          << std::endl;

}

void ext_GAUSS::get_conn() {
  struct ent_conn* e_c = new ent_conn();

  int dim = 0;
  int sz = cb->get_sz(dim);

  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      cb->get_conn_dim_gmi(3, dim, c_id+1, e_c);
      c0_c3.at(c_id) = e_c->conn;
      cb->get_conn_dim_gmi(1, dim, c_id+1, e_c);
      c0_c1.at(c_id) = e_c->conn;
    }
  }
  delete e_c;
}

void ext_GAUSS::get_len() {
  int dim = 1;
  int sz = cb->get_sz(dim);
  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      c1_len.at(c_id) = vd_meas_set(m, &e_list->e.at(dim).at(c_id).at(dim));
    }
  }
  dim = 3;
  sz = cb->get_sz(dim);
  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      c3_vol.at(c_id) = vd_meas_set(m, &e_list->e.at(dim).at(c_id).at(dim));
    }
  }
}

void ext_GAUSS::process_mesh() {
  clear();
  collect_ori();
  c3_euler_mod.resize(cb->get_sz(3));

  c0_ang_def.resize(cb->get_sz(0));
  c0_ang_dev.resize(cb->get_sz(0));
  c0_vel.resize(cb->get_sz(0));

  c0_c1.resize(cb->get_sz(0));
  c0_c3.resize(cb->get_sz(0));

  c1_len.resize(cb->get_sz(1));
  c3_vol.resize(cb->get_sz(3));

  calc_gauss();
  // Get the 1- and 3-strata connections of 0-strata.
  get_conn();
  // Collect the lengths and volumes of 1- and 3-strata.
  get_len();

  gauss = CoM + ang_defect + ang_defect_trans;
}

void ext_GAUSS::set_bc_type(VD_BC_TYPE ext_BC_in) {
  ext_BC_type = ext_BC_in;
}

// Going over each vertex, 
//double ext_GAUSS::validation() {
//}

void ext_GAUSS::clear() {
  CoM = 0.;
  ang_defect = 0.;
  ang_defect_trans = 0.;
  ang_defect_const = 0.;
  ang_defect_dev = 0.;
  gauss = 0.;

  ang_defect_1c_int = 0;

  ori_map.clear();
  ang_map.clear();
  ang_map_count.clear();

  pos_map.clear();
  c3_euler_mod.clear();

  c0_ang_def.clear();
  c0_ang_dev.clear();
  c0_vel.clear();
  c2_map.clear();
  c2_ang.clear();

  for(int i = 0; i < c0_c1.size(); i++)
    c0_c1.at(i).clear();
  c0_c1.clear();
  for(int i = 0; i < c0_c3.size(); i++)
    c0_c3.at(i).clear();
  c0_c3.clear();
}

ext_GAUSS::ext_GAUSS() : cb(NULL), m(NULL), f_calc(NULL), e_list(NULL),
                   vd_3c(NULL),
                   c2_map{}, c2_ang{},
                   c3_euler_mod(0), 
                   c0_ang_def(0), c0_ang_dev(0), c0_vel(0),
                   c0_c1(0, std::vector<int>(0)), 
                   c0_c3(0, std::vector<int>(0)), 
                   c3_vol(0),
                   c1_len(0),
                   c3_nbr(0), c_id(-1),
                   ori_map{}, pos_map{},
                   ang_map{},
                   ang_map_count{},
                   CoM(0.), ang_defect(0.), ang_defect_trans(0.),
                   ang_defect_const(0.), ang_defect_dev(0.), gauss(0.),
                   ang_def_equi(0.),
                   temp_pos(0,0,0),
                   ext_BC_type(VD_BC_TYPE::NEUMANN_EXT) {
}

ext_GAUSS::ext_GAUSS(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* e_list_in, 
                     field_calc* f_calc_in) : 
                   cb(cb_in), m(m_in), f_calc(f_calc_in),
                   e_list(e_list_in), 
                   vd_3c(NULL),
                   c2_map{}, c2_ang{},
                   c3_euler_mod(0), 
                   c0_ang_def(0), c0_ang_dev(0), c0_vel(0),
                   c0_c1(0, std::vector<int>(0)), 
                   c0_c3(0, std::vector<int>(0)), 
                   c3_vol(0),
                   c1_len(0),
                   c3_nbr(0), c_id(-1),
                   ori_map{}, pos_map{},
                   ang_map{},
                   ang_map_count{},
                   CoM(0.), ang_defect(0.), ang_defect_trans(0.),
                   ang_defect_const(0.), ang_defect_dev(0.), gauss(0.),
                   ang_def_equi(0.),
                   temp_pos(0,0,0),
                   ext_BC_type(VD_BC_TYPE::NEUMANN_EXT) {

  vd_3c = new vd_3c_det();
}

ext_GAUSS::~ext_GAUSS() {
  clear();
  delete vd_3c;
}
/*
/////////////////////////////////
// ext_GAUSS_g
/////////////////////////////////
double ext_GAUSS_g::get_ang(apf::MeshEntity* t) {
  apf::Downward d_v;
  apf::Downward d_e;

  m->getDownward(t, 0, d_v);
  m->getDownward(t, 1, d_e);

  int lookup_tri_ed [3][2] = {{2,0},{0,1},{1,2}};

  int v1 = findIn(d_v, 3, v);
  int e1 = lookup_tri_ed[v1][0];
  int e2 = lookup_tri_ed[v1][1];

  double coef = 1.;
  apf::Vector3 temp(0,0,0);
  temp = norm_0(pos_map[d_e[e1]] - temp_pos);
  if(temp*ori_map[d_e[e1]] < -std::numeric_limits<double>::min())
    coef = -coef;
  temp = norm_0(pos_map[d_e[e2]] - temp_pos);
  if(temp*ori_map[d_e[e2]] < -std::numeric_limits<double>::min())
    coef = -coef;
  double ang_cos = ori_map[d_e[e2]]*ori_map[d_e[e1]]*coef;
  return std::acos(std::min(std::max(ang_cos,-1.0),1.0));
  //assert(ang_cos >= -1. and ang_cos <= 1.);
  //return std::acos(ang_cos);
}

void ext_GAUSS_g::collect_ang() {
  std::vector<apf::MeshEntity*> t_list(0);
  for(int i = 0; i < g_id.size(); i++) {
    int c3_curr = g_id.at(i);
    std::vector<int> c2_bound(0);
    c2_bound = vd_3c->get_bound_ids(3, c3_curr);
    for(int i = 0; i < c2_bound.size(); i++) {
      if(!c2_map[c2_bound.at(i)]) {
        double ang_curr = 0.;
        t_list = vd_3c->get_dim_ent(2, c2_bound.at(i));
        for(int j = 0; j < t_list.size(); j++) {
          ang_curr = ang_curr + ang_map[t_list.at(j)][v];
          ang_map_count[t_list.at(j)][v] = ang_map_count[t_list.at(j)][v] + 1;
        }
        c2_map[c2_bound.at(i)] = true;
        c2_ang[c2_bound.at(i)] = ang_curr;
      }
    }
  }
}

void ext_GAUSS_g::load_v() {
  vd_3c->reload(m, cb, v);
  c3_nbr = vd_3c->get_3c_nbr();
  m->getPoint(v, 0, temp_pos);
  c2_map.clear();
  c2_ang.clear();
  collect_ang();
}

void ext_GAUSS_g::process_0c_vertex() {
  load_v();

  for(int c3_curr = 0; c3_curr < c3_nbr; c3_curr++) {
    apf::ModelEntity* mdl_3c = vd_3c->get_mdl(3, c3_curr);
    int c3_id = m->getModelTag(mdl_3c);
    std::vector<int> c2_bound(0);
    c2_bound = vd_3c->get_bound_ids(3, c3_curr);

    double ang_curr = 0.;
    for(int i = 0; i < c2_bound.size(); i++) {
      ang_curr = ang_curr + c2_ang[c2_bound.at(i)];
    }
    double ang_def_curr = 2*PI_L - ang_curr;
    c3_euler_mod.at(c3_id - 1) = c3_euler_mod.at(c3_id - 1) + ang_def_curr;
    c0_ang_def.at(c_id) = c0_ang_def.at(c_id) + ang_def_curr;

    if(ext_curr) {
      CoM = CoM + ang_def_curr;
//    if(ext_curr and 
//       ext_BC_type == VD_BC_TYPE::NEUMANN_EXT) {
//      CoM = CoM + ang_def_curr;
//    }
//    else if(ext_curr and 
//      (f_calc->get_proj() == (int)PROJ_TYPE::FIXED)) {
//     CoM = CoM + ang_def_curr;
      c0_ang_dev.at(c_id) = -1;
    }
    else {
      // TODO for exterior vertices, use the insertible 0-stratum detection
      // from topo_vd.cc.
      if(c3_nbr == 4) {
        ang_defect = ang_defect + ang_def_curr;
        ang_defect_const = ang_defect_const + ang_def_equi;
        ang_defect_dev = ang_defect_dev + (ang_def_curr - ang_def_equi);

        c0_ang_dev.at(c_id) = ang_def_curr - ang_def_equi;
      }
      else {
        ang_defect_trans = ang_defect_trans + ang_def_curr;
        c0_ang_dev.at(c_id) = -1;
      }
    }
  }
  apf::Field* vel_field = m->findField("velocity_field");
  apf::Vector3 v_temp(0,0,0);
  apf::getVector(vel_field, v, 0, v_temp);

  c0_vel.at(c_id) = v_temp.getLength();
}

void ext_GAUSS_g::process_ext_1c_vertex() {
  load_v();

  for(int c3_curr = 0; c3_curr < c3_nbr; c3_curr++) {
    apf::ModelEntity* mdl_3c = vd_3c->get_mdl(3, c3_curr);
    int c3_id = m->getModelTag(mdl_3c);
    std::vector<int> c2_bound(0);
    c2_bound = vd_3c->get_bound_ids(3, c3_curr);

    double ang_curr = 0.;
    for(int i = 0; i < c2_bound.size(); i++) {
      ang_curr = ang_curr + c2_ang[c2_bound.at(i)];
    }
    double ang_def_curr = 2*PI_L - ang_curr;
    c3_euler_mod.at(c3_id - 1) = c3_euler_mod.at(c3_id - 1) + ang_def_curr;
    CoM = CoM + ang_def_curr;
  }
}

void ext_GAUSS_g::process_int_1c_vertex() {
  load_v();

  for(int c3_curr = 0; c3_curr < c3_nbr; c3_curr++) {
    apf::ModelEntity* mdl_3c = vd_3c->get_mdl(3, c3_curr);
    int c3_id = m->getModelTag(mdl_3c);
    std::vector<int> c2_bound(0);
    c2_bound = vd_3c->get_bound_ids(3, c3_curr);

    double ang_curr = 0.;
    for(int i = 0; i < c2_bound.size(); i++) {
      ang_curr = ang_curr + c2_ang[c2_bound.at(i)];
    }
    double ang_def_curr = 2*PI_L - ang_curr;
    c3_euler_mod.at(c3_id - 1) = c3_euler_mod.at(c3_id - 1) + ang_def_curr;
    ang_defect_1c_int = ang_defect_1c_int + ang_def_curr;
  }
}

void ext_GAUSS_g::process_2c_vertex() {
  m->getPoint(v, 0, temp_pos);
  std::vector<apf::MeshEntity*> edges(0);
  std::vector<apf::MeshEntity*> tris(0);
  //std::map<apf::MeshEntity*, double> ang_map{};

  vd_set_up(m, v, &edges);
  vd_set_up(m, &edges, &tris);

  ent_conn e_conn;
  cb->get_conn_dim(3, 2, c_id, &e_conn);

  std::vector<apf::MeshEntity*> t_list(0);
  for(int i = 0; i  < e_conn.conn.size(); i++) {
    int c3_curr = e_conn.conn.at(i);

    double ang_curr = 0.;
    for(int j = 0; j < tris.size(); j++) {
      apf::ModelEntity* mdl_tri = m->toModel(tris.at(j));
      if(m->getModelType(mdl_tri) == 2) {
        ang_curr = ang_curr + ang_map[tris.at(j)][v];
        ang_map_count[tris.at(j)][v] = ang_map_count[tris.at(j)][v] + 1;
      }
    }
    double ang_def_curr = 2*PI_L - ang_curr;
    CoM = CoM + ang_def_curr;
    c3_euler_mod.at(c3_curr) = c3_euler_mod.at(c3_curr) + ang_def_curr;

    apf::setScalar(ang_def_field, v, 0, ang_def_curr);
  }
}

double ext_GAUSS_g::get_ang_val(apf::MeshEntity* tri, apf::MeshEntity* vert) {
  return ang_map[tri][vert];
}

int ext_GAUSS_g::get_ang_count(apf::MeshEntity* tri, apf::MeshEntity* vert) {
  return ang_map_count[tri][vert];
}

void ext_GAUSS_g::calc_gauss() {
  // Check the total of the gaussian curvatures by summing up the interior angles 
  // of the triangles bounding a 3-stratum and subtracting from 2*pi*n_v, where
  // n_v is the number of vertices bounding the 3-stratum.

  ang_def_field = vd_att_vs_field(m, "a_def");

  dim = 3;
  struct ent_conn* e_c = new ent_conn();
  for(int i = 0; i < g_id.size(); i++) {
    c_id = g_id.at(i);
    assert(!cb->is_free(dim, c_id));
    int v_nbr = 0;
    cb->get_conn_dim(0, dim, c_id, e_c);
    v_nbr = v_nbr + e_c->conn.size();

    cb->get_conn_dim(1, dim, c_id, e_c);
    for(int i = 0; i < e_c->conn.size(); i++) {
      v_nbr = v_nbr + e_list->e.at(1).at(e_c->conn.at(i)).at(0).size();
    }

    cb->get_conn_dim(2, dim, c_id, e_c);
    for(int i = 0; i < e_c->conn.size(); i++) {
      v_nbr = v_nbr + e_list->e.at(2).at(e_c->conn.at(i)).at(0).size();
    }

    apf::Downward d_v;
    double ang_glob = 0.;
    for(int i = 0; i < e_c->conn.size(); i++) {
      int c2_c = e_c->conn.at(i);
      std::vector<apf::MeshEntity*>* ents = &e_list->e.at(2).at(c2_c).at(2);
      for(int j = 0; j < ents->size(); j++) {
        m->getDownward(ents->at(j), 0, d_v);
        //std::cout << "Tri " << ents->at(i) << std::endl;
        for(int k = 0; k < 3; k++) {
          ang_glob = ang_glob + ang_map[ents->at(j)][d_v[k]];
        }
      }
    }
    std::cout << "3c" << c_id + 1 << " " << 2*v_nbr*PI_L - ang_glob 
              << std::endl;
  }

  bool has_ext = (f_calc->get_proj() == (int)PROJ_TYPE::FIXED);
  bool has_shell = (f_calc->get_proj() == (int)PROJ_TYPE::EXT_SHELL);

  std::map<int, bool> c_pass{};

  dim = 0;

  for(int i = 0; i < g_id.size(); i++) {
    int c3_curr = g_id.at(i);
    cb->get_conn_dim(0, 3, c3_curr, e_c);
    for(int j = 0; j < e_c->conn.size(); j++) {
      c_id = e_c->conn.at(j);

      if(!c_pass[c_id]) {
        std::vector<apf::MeshEntity*>* verts = &e_list->e.at(dim).at(c_id).at(0);
        if(has_ext and cb->get_cell_ext(dim, c_id)) {
          ext_curr = true;
        }
        else if(has_shell and f_calc->get_e_sh()->chk_shell(dim, c_id) ) {
          ext_curr = true;
          shell sh = f_calc->get_e_sh()->get_shell(dim, c_id);
          if(sh.dim == 0) {
            ang_def_equi = ANG_EQUI_CONST_EXT_0C;
          }
          else if(sh.dim == 1) {
            ang_def_equi = ANG_EQUI_CONST_EXT_1C;
          }
          else if(sh.dim == 2)  {
            ang_def_equi = ANG_EQUI_CONST_EXT_2C;
          }
        }
        else {
          ext_curr = false;
          ang_def_equi = ANG_EQUI_CONST;
        }
        assert(verts->size() == 1);
        v = verts->at(0);
        process_0c_vertex();
        c_pass[c_id] = true;
      }
    }
  }
  c_pass.clear();

  dim = 1;
  for(int i = 0; i < g_id.size(); i++) {
    int c3_curr = g_id.at(i);
    cb->get_conn_dim(1, 3, c3_curr, e_c);
    for(int j = 0; j < e_c->conn.size(); j++) {
      c_id = e_c->conn.at(j);

      if(!c_pass[c_id]) {
        std::vector<apf::MeshEntity*>* verts = &e_list->e.at(dim).at(c_id).at(0);
        if(cb->get_cell_ext(dim, c_id)) {
          for(int i = 0; i < verts->size(); i++) {
            v = verts->at(i);
            process_ext_1c_vertex();
          }
        }
        else {
          for(int i = 0; i < verts->size(); i++) {
            v = verts->at(i);
            process_int_1c_vertex();
          }
        }
      }
      c_pass[c_id] = true;
    }
  }
  c_pass.clear();
  dim = 2;
  for(int i = 0; i < g_id.size(); i++) {
    int c3_curr = g_id.at(i);
    cb->get_conn_dim(1, 3, c3_curr, e_c);
    for(int j = 0; j < e_c->conn.size(); j++) {
      c_id = e_c->conn.at(j);
      if(!c_pass[c_id]) {

        std::vector<apf::MeshEntity*>* verts = &e_list->e.at(dim).at(c_id).at(0);
        //if(cb->get_cell_ext(dim, c_id)) {
          // TODO This is currently a placeholder. Would be used for free surface
          // BC that allows grooving. 
          // For Neumann BC, this is the Gaussian curv. is 0.
          //for(int j = 0; j < verts->size(); j++) {
          //  process_ext_2c_vertex(verts->at(j));
          //}
        //}
        //else {
          for(int i = 0; i < verts->size(); i++) {
            v = verts->at(i);
            process_2c_vertex();
          }
        //}
      }
      c_pass[c_id] = true;
    }
  }
  delete e_c;

  apf::writeVtkFiles("./output/gauss", m);
}

void ext_GAUSS_g::collect_ori() {
  for(dim = 1; dim < 3; dim++) {
    int sz = cb->get_sz(dim);
    for(c_id = 0; c_id < sz; c_id++) {
      if(!cb->is_free(dim, c_id)) {
        std::vector<apf::MeshEntity*>* ents = &e_list->e.at(dim).at(c_id).at(1);
        for(int i = 0; i < ents->size(); i++) {
          ori_map[ents->at(i)] = get_edge_dir(m, ents->at(i));
          pos_map[ents->at(i)] = vd_get_pos(m, ents->at(i));
        }
      }
    }
  }

  apf::Downward d_v;
  dim = 2;
  int sz = cb->get_sz(dim);
  //std::cout << "Total interior angle deviation calculation." << std::endl;
  double ang_glob = 0.;
  int t_nbr = 0;
  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      std::vector<apf::MeshEntity*>* ents = &e_list->e.at(dim).at(c_id).at(2);
      for(int i = 0; i < ents->size(); i++) {
        m->getDownward(ents->at(i), 0, d_v);
        //std::cout << "Tri " << ents->at(i) << std::endl;
        double ang_tot = 0;
        for(int j = 0; j < 3; j++) {
          v = d_v[j];
          m->getPoint(v, 0, temp_pos);
          double ang_curr = get_ang(ents->at(i));
          //std::cout << "\t" << ang_curr << std::endl;
          ang_map[ents->at(i)][v] = ang_curr;
          ang_tot = ang_tot + ang_map[ents->at(i)][v];
        }
        //std::cout << "Tri " << ents->at(i) 
        //          << " " << PI_L - ang_tot << std::endl;
        ang_glob = ang_glob + ang_tot;
      }
      t_nbr = t_nbr + ents->size();
    }
  }
  //std::cout << "Total of interior angles " << ang_glob
  //          << " " << " n_tris*pi - ang_glob " << t_nbr*PI_L - ang_glob
  //          << std::endl;

}

void ext_GAUSS_g::get_conn() {
  struct ent_conn* e_c = new ent_conn();

  int dim = 0;
  int sz = cb->get_sz(dim);

  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      cb->get_conn_dim_gmi(3, dim, c_id+1, e_c);
      c0_c3.at(c_id) = e_c->conn;
      cb->get_conn_dim_gmi(1, dim, c_id+1, e_c);
      c0_c1.at(c_id) = e_c->conn;
    }
  }
  delete e_c;
}

void ext_GAUSS_g::get_len() {
  int dim = 1;
  int sz = cb->get_sz(dim);
  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      c1_len.at(c_id) = vd_meas_set(m, &e_list->e.at(dim).at(c_id).at(dim));
    }
  }
  dim = 3;
  sz = cb->get_sz(dim);
  for(c_id = 0; c_id < sz; c_id++) {
    if(!cb->is_free(dim, c_id)) {
      c3_vol.at(c_id) = vd_meas_set(m, &e_list->e.at(dim).at(c_id).at(dim));
    }
  }
}

void ext_GAUSS_g::process_mesh() {
  clear();
  collect_ori();
  c3_euler_mod.resize(cb->get_sz(3));

  c0_ang_def.resize(cb->get_sz(0));
  c0_ang_dev.resize(cb->get_sz(0));
  c0_vel.resize(cb->get_sz(0));

  c0_c1.resize(cb->get_sz(0));
  c0_c3.resize(cb->get_sz(0));

  c1_len.resize(cb->get_sz(1));
  c3_vol.resize(cb->get_sz(3));

  calc_gauss();
  // Get the 1- and 3-strata connections of 0-strata.
  get_conn();
  // Collect the lengths and volumes of 1- and 3-strata.
  get_len();

  gauss = CoM + ang_defect + ang_defect_trans;
}

void ext_GAUSS_g::set_bc_type(VD_BC_TYPE ext_BC_in) {
  ext_BC_type = ext_BC_in;
}

// Going over each vertex, 
//double ext_GAUSS_g::validation() {
//}

void ext_GAUSS_g::clear() {
  CoM = 0.;
  ang_1c = 0;
  ang_0c = 0;

  ori_map.clear();
  ang_map.clear();
  ang_map_count.clear();

  pos_map.clear();
  c3_euler_mod.clear();

  c0_ang_def.clear();
  c0_ang_dev.clear();
  c0_vel.clear();
  c2_map.clear();
  c2_ang.clear();

  for(int i = 0; i < c0_c1.size(); i++)
    c0_c1.at(i).clear();
  c0_c1.clear();
  for(int i = 0; i < c0_c3.size(); i++)
    c0_c3.at(i).clear();
  c0_c3.clear();
}

ext_GAUSS_g::ext_GAUSS_g() : cb(NULL), m(NULL), f_calc(NULL), e_list(NULL),
                   vd_3c(NULL),
                   c2_map{}, c2_ang{},
                   c3_euler_mod(0), 
                   c0_ang_def(0), c0_ang_dev(0), c0_vel(0),
                   c0_c1(0, std::vector<int>(0)), 
                   c0_c3(0, std::vector<int>(0)), 
                   c3_vol(0),
                   c1_len(0),
                   c3_nbr(0), c_id(-1),
                   ori_map{}, pos_map{},
                   ang_map{},
                   ang_map_count{},
                   CoM(0.), ang_defect(0.), ang_defect_trans(0.),
                   ang_defect_const(0.), ang_defect_dev(0.), gauss(0.),
                   ang_def_equi(0.),
                   temp_pos(0,0,0),
                   g_id(0),
                   ext_BC_type(VD_BC_TYPE::NEUMANN_EXT) {
}

ext_GAUSS_g::ext_GAUSS_g(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* e_list_in, 
                     field_calc* f_calc_in) : 
                   cb(cb_in), m(m_in), f_calc(f_calc_in),
                   e_list(e_list_in), 
                   vd_3c(NULL),
                   c2_map{}, c2_ang{},
                   c3_euler_mod(0), 
                   c0_ang_def(0), c0_ang_dev(0), c0_vel(0),
                   c0_c1(0, std::vector<int>(0)), 
                   c0_c3(0, std::vector<int>(0)), 
                   c3_vol(0),
                   c1_len(0),
                   c3_nbr(0), c_id(-1),
                   ori_map{}, pos_map{},
                   ang_map{},
                   ang_map_count{},
                   CoM(0.), ang_defect(0.), ang_defect_trans(0.),
                   ang_defect_const(0.), ang_defect_dev(0.), gauss(0.),
                   ang_def_equi(0.),
                   temp_pos(0,0,0),
                   g_id(0),
                   ext_BC_type(VD_BC_TYPE::NEUMANN_EXT) {

  vd_3c = new vd_3c_det();
}

ext_GAUSS_g::~ext_GAUSS_g() {
  clear();
  delete vd_3c;
}
*/

/////////////////////////////////
// vd_extractor
/////////////////////////////////

EXT_TYPE conv_str2ext(const std::string& str) {
  if(str == "MS") return EXT_TYPE::MS;
  else if(str == "VEL") return EXT_TYPE::VEL;
  else if(str == "POS") return EXT_TYPE::POS;
  else if(str == "VEL_DIR") return EXT_TYPE::VEL_DIR;
  else if(str == "POS_DIR") return EXT_TYPE::POS_DIR;
  else if(str == "CRV") return EXT_TYPE::CRV;
  else if(str == "MEAS") return EXT_TYPE::MEAS;
  else if(str == "AREA_T") return EXT_TYPE::AREA_T;
  else if(str == "GAUSS") return EXT_TYPE::GAUSS;
  else if(str == "GAUSS_G") return EXT_TYPE::GAUSS_G;
  else if(str == "ROC") return EXT_TYPE::ROC;
  else if(str == "COUNT") return EXT_TYPE::COUNT;
  else if(str == "ANGLE1") return EXT_TYPE::ANGLE1;
  else if(str == "ANGLE2") return EXT_TYPE::ANGLE2;
  return EXT_TYPE::END;
}

VD_BC_TYPE conv_str2bc(const std::string& str) {
  if(str == "PERIODIC") return VD_BC_TYPE::PERIODIC;
  else if(str == "PERIODIC_CONST") return VD_BC_TYPE::PERIODIC_CONST;
  else if(str == "NEUMANN_EXT") return VD_BC_TYPE::NEUMANN_EXT;
  else if(str == "NEUMANN_CONST") return VD_BC_TYPE::NEUMANN_CONST;
  return VD_BC_TYPE::END;
}

void vd_extractor::parse_opt(std::vector<std::string> &opts) {
}

vd_extractor::vd_extractor() : 
                            m(NULL), cb(NULL), e_list(NULL), 
                            t_curr(0) {
}

vd_extractor::vd_extractor(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                            m(NULL), cb(NULL), e_list(NULL), f_calc(NULL),
                            t_curr(t_in) {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  parse_opt(opts_in);
  write_csv();
}

void vd_extractor::write_csv() {
}

vd_extractor::~vd_extractor() {
//
}

vd_ext_MS::vd_ext_MS(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                                e_MS(m_in, cb_in), 
                                OUTPUTMS(""),
                                OUTPUTROC(""),
                                OUTPUTVOL("")  {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  parse_opt(opts_in);
  write_csv();
}

// A wrapper for ext_MS in style of vd_extractor. 
// TODO incorporate ext_MS into this.
void vd_ext_MS::parse_opt(std::vector<std::string> &opts) {
  e_MS.process_mesh();
  OUTPUTMS = opts.at(1);
  OUTPUTROC = opts.at(2);
  OUTPUTVOL = opts.at(3);
  if(opts.size() > 4)
    dt = std::stod(opts.at(4));
  else
    dt = -1;
}

void vd_ext_MS::write_csv() {
  csvfile csvMS(OUTPUTMS);
  csvfile csvROC(OUTPUTROC);
  csvfile csvVOL(OUTPUTVOL);

  std::vector<double> roc(cb->get_sz(3), 0);
  std::vector<apf::MeshEntity*> tet_bound(0);
  csvMS << t_curr;
  csvROC << t_curr;
  csvVOL << t_curr;

  std::cout << "t: " <<  t_curr << std::endl;
  for(int j = 0; j < cb->get_sz(3); j++) {
    if(!cb->is_free(3,j)) {
      tet_bound = vd_get_3c_tet_bound(m, e_list->e.at(3).at(j).at(3));
      if(dt < 0.) {
        roc.at(j) = vd_calc_roc(m, &tet_bound, f_calc);
      }
      else {
        roc.at(j) = vd_calc_roc_vol_high(m, &tet_bound, f_calc, dt);
      }
      std::cout << "3c" << j + 1 
                << " MS: " << e_MS.dVdt_c3[j+1]*f_calc->vdparam.v_mult 
                << " roc: " << roc.at(j) << std::endl; 
      csvMS << e_MS.dVdt_c3[j+1]*f_calc->vdparam.v_mult;
      csvROC << roc.at(j);
      csvVOL << vd_meas_set(m, &e_list->e.at(3).at(j).at(3));
    }
    else {
      csvMS << 0;
      csvROC << 0;
      csvVOL << 0;
    }
  }
  csvMS << endrow;
  csvROC << endrow;
  csvVOL << endrow;
}

vd_ext_GAUSS::vd_ext_GAUSS(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                                e_GAUSS(m_in, cb_in, e_list_in, f_calc_in), 
                                OUTPUTGAUSS(""), 
                                OUTPUTEULERMOD(""), OUTPUTERROR(""),
                                OUTPUT0C_dev(""), 
                                OUTPUT1C(""), OUTPUT3C(""), 
                                OUTPUT1C_conn(""), OUTPUT3C_conn(""),
                                OUTPUT0C_vel("") {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  parse_opt(opts_in);
  write_csv();
}

// A wrapper for ext_MS in style of vd_extractor. 
// TODO incorporate ext_MS into this.
void vd_ext_GAUSS::parse_opt(std::vector<std::string> &opts) {
  OUTPUTGAUSS = opts.at(1);
  OUTPUTEULERMOD = opts.at(2);
  OUTPUTERROR = opts.at(3);

  OUTPUT0C_dev = opts.at(4);
  OUTPUT1C = opts.at(5);
  OUTPUT3C = opts.at(6);
  OUTPUT1C_conn = opts.at(7);
  OUTPUT3C_conn = opts.at(8);
  OUTPUT0C_vel = opts.at(9);

  bc_type = conv_str2bc(opts.at(10));
  e_GAUSS.set_bc_type(bc_type);
  e_GAUSS.process_mesh();
}

void vd_ext_GAUSS::write_csv() {
  csvfile csvGAUSS(OUTPUTGAUSS);

  csvGAUSS << t_curr;
  // Integral of the GC over the 2-stratum vertices. Also includes exterior if
  // bc_type == VD_BC_TYPE::NEUMANN_EXT
  csvGAUSS << e_GAUSS.CoM;
  // Sum of the angular defects at 0-strata not under transition.
  csvGAUSS << e_GAUSS.ang_defect;
  // Angular defect due to transition. This is used to distinguish the error
  // due to deviation from the assumed angular defect configuration and the
  // error due to unstable configuration during transition.
  csvGAUSS << e_GAUSS.ang_defect_trans;
  // Angular defect as calculated by isotropic GBE constant angular defect 
  // assumption: f_0 * (2*pi - 3*acos(-1/3)), f_0 number of 0-strata
  csvGAUSS << e_GAUSS.ang_defect_const;
  // The sum of the deviations of the angular defects from the constant GBE 
  // case.
  csvGAUSS << e_GAUSS.ang_defect_dev;
  // Gaussian curvature.
  csvGAUSS << e_GAUSS.gauss;

  // Integral of the Gaussian curvature corresponding to the interior 1strata.
  csvGAUSS << e_GAUSS.ang_defect_1c_int;

  csvGAUSS << endrow;

  csvfile csvEM(OUTPUTEULERMOD);
  csvEM << t_curr;
  for(int i = 0; i < e_GAUSS.c3_euler_mod.size(); i++)
    csvEM << e_GAUSS.c3_euler_mod.at(i);
  csvEM << endrow;

  csvfile csv0Cdev(OUTPUT0C_dev);
  csv0Cdev << t_curr;
  for(int i = 0; i < e_GAUSS.c0_ang_dev.size(); i++)
    csv0Cdev << e_GAUSS.c0_ang_dev.at(i);
  csv0Cdev << endrow;

  csvfile csv0Cvel(OUTPUT0C_vel);
  csv0Cvel << t_curr;
  for(int i = 0; i < e_GAUSS.c0_vel.size(); i++)
    csv0Cvel << e_GAUSS.c0_vel.at(i);
  csv0Cvel << endrow;

  // The 1-/3-strata connectivities of the 0-strata.
  csvfile csv1c(OUTPUT1C_conn);
  csv1c << t_curr;
  for(int i = 0; i < e_GAUSS.c0_c1.size(); i++) {
    std::string conn("");
    for(int j = 0; j < e_GAUSS.c0_c1.at(i).size(); j++) {
      conn = conn + std::to_string(e_GAUSS.c0_c1.at(i).at(j)) + " ";
    }
    csv1c << conn;
  }
  csv1c << endrow;

  csvfile csv3c(OUTPUT3C_conn);
  csv3c << t_curr;
  for(int i = 0; i < e_GAUSS.c0_c3.size(); i++) {
    std::string conn("");
    for(int j = 0; j < e_GAUSS.c0_c3.at(i).size(); j++) {
      conn = conn + std::to_string(e_GAUSS.c0_c3.at(i).at(j)) + " ";
    }
    csv3c << conn;
  }
  csv3c << endrow;

  // The 1-/3-strata lengths/volumes.
  csvfile csv1clen(OUTPUT1C);
  csv1clen << t_curr;
  for(int i = 0; i < e_GAUSS.c1_len.size(); i++) {
    csv1clen << e_GAUSS.c1_len.at(i);
  }
  csv1clen << endrow;

  csvfile csv3cvol(OUTPUT3C);
  csv3cvol << t_curr;
  for(int i = 0; i < e_GAUSS.c3_vol.size(); i++) {
    csv3cvol << e_GAUSS.c3_vol.at(i);
  }
  csv3cvol << endrow;
}


vd_ext_GAUSS_g::vd_ext_GAUSS_g(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                                g_id(0), 
                                CoM(0), 
                                CoM_avg(0.),
                                ang_def_dev_0c(0.),
                                ang_def_dev_1c(0.),
                                ori_map{}, pos_map{}, c_pass{}, 
                                ang_def_2c{}, ang_def_1c{}, ang_def_0c{},
                                OUTPUTGAUSS("") {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  parse_opt(opts_in);
  write_csv();
}


void vd_ext_GAUSS_g::collect_ori() {
  for(int dim = 1; dim < 3; dim++) {
    int sz = cb->get_sz(dim);
    for(int c_id = 0; c_id < sz; c_id++) {
      if(!cb->is_free(dim, c_id)) {
        std::vector<apf::MeshEntity*>* ents = &e_list->e.at(dim).at(c_id).at(1);
        for(int i = 0; i < ents->size(); i++) {
          ori_map[ents->at(i)] = get_edge_dir(m, ents->at(i));
          pos_map[ents->at(i)] = vd_get_pos(m, ents->at(i));
        }
      }
    }
  }
}

// To calculate the integral of the gaussian curvature corresponding to the 2-stratum, calculate the total of the angles around the 2-stratum vertices and subtract that from 2*pi*n_v, where n_v is the number of 2-stratum vertices.
void vd_ext_GAUSS_g::calc_ang_2c(int c_in) {
  if(c_pass[c_in])
    return;

  apf::ModelEntity* mdl_2c = m->findModelEntity(2, c_in+1);

  double ang_tot = 0.;
  double ang_dev_1c = 0.;
  double ang_dev_0c = 0.;

  apf::Downward d_v;
  apf::Downward d_e;
  int lookup_tri_ed [3][2] = {{2,0},{0,1},{1,2}};

  std::vector<apf::MeshEntity*> * t_list = &e_list->e.at(2).at(c_in).at(2);
  apf::MeshEntity* t_curr;

  apf::Vector3 temp_pos(0,0,0);
  apf::Vector3 temp(0,0,0);

  for(int i = 0; i < t_list->size(); i++) {
    t_curr = t_list->at(i);
    m->getDownward(t_curr, 0, d_v);
    m->getDownward(t_curr, 1, d_e);

    for(int v1 = 0; v1 < 3; v1++) {
      int e1 = lookup_tri_ed[v1][0];
      int e2 = lookup_tri_ed[v1][1];
      apf::ModelEntity* mdl_v = m->toModel(d_v[v1]);

      m->getPoint(d_v[v1], 0, temp_pos);

      double coef = 1.;
      temp = norm_0(pos_map[d_e[e1]] - temp_pos);
      if(temp*ori_map[d_e[e1]] < -std::numeric_limits<double>::min())
        coef = -coef;
      temp = norm_0(pos_map[d_e[e2]] - temp_pos);
      if(temp*ori_map[d_e[e2]] < -std::numeric_limits<double>::min())
        coef = -coef;
      double ang_cos = ori_map[d_e[e2]]*ori_map[d_e[e1]]*coef;

      double ang_curr = std::acos(std::min(std::max(ang_cos,-1.0),1.0));

      if(mdl_v == mdl_2c) {
        ang_tot = ang_tot + ang_curr;
      }
      else if(m->getModelType(mdl_v) == 1) {
        int c_curr = m->getModelTag(mdl_v) - 1;
        ang_def_1c[c_curr] = ang_def_1c[c_curr] + ang_curr;
      }
      else {
        int c_curr = m->getModelTag(mdl_v) - 1;
        ang_def_0c[c_curr] = ang_def_0c[c_curr] + ang_curr;
      }
    }
  }

  int v_sz = e_list->e.at(2).at(c_in).at(0).size();
  ang_def_2c[c_in] = 2*v_sz*PI_L - ang_tot;
 
  c_pass[c_in] = true;
}

void vd_ext_GAUSS_g::collect_ang() {
  struct ent_conn* e_c = new ent_conn();

  for(int i = 0; i < g_id.size(); i++) {
    int c3_curr = g_id.at(i);
    cb->get_conn_dim(1, 3, c3_curr, e_c);

    for(int j = 0; j < e_c->conn.size(); j++) {
      int c_id = e_c->conn.at(j);
      std::vector<apf::MeshEntity*> * ent_list = &e_list->e.at(1).at(c_id).at(1);
      for(int k = 0; k < ent_list->size(); k++) {
        ori_map[ent_list->at(k)] = get_edge_dir(m, ent_list->at(k));
        pos_map[ent_list->at(k)] = vd_get_pos(m, ent_list->at(k));
      }
    }
  }

  apf::Downward d_v;
  apf::Downward d_e;

  for(int i = 0; i < g_id.size(); i++) {
    int c3_curr = g_id.at(i);
    cb->get_conn_dim(2, 3, c3_curr, e_c);
    ang_def_1c.clear();
    ang_def_0c.clear();

    for(int j = 0; j < e_c->conn.size(); j++) {
      int c_id = e_c->conn.at(j);
      calc_ang_2c(c_id);
      CoM.at(i) = CoM.at(i) + ang_def_2c[c_id];
    }
    CoM_avg = CoM_avg + CoM.at(i);

    double ang_def_0c_dev = 0.;
    double ang_def_1c_dev = 0.;
    cb->get_conn_dim(1, 3, c3_curr, e_c);
    for(int j = 0; j < e_c->conn.size(); j++) {
      int c_id = e_c->conn.at(j);
      int v_sz = e_list->e.at(1).at(c_id).at(0).size();
      ang_def_1c_dev = ang_def_1c_dev + 2*v_sz*PI_L - ang_def_1c[c_id];
    }

    cb->get_conn_dim(0, 3, c3_curr, e_c);
    for(int j = 0; j < e_c->conn.size(); j++) {
      int c_id = e_c->conn.at(j);
      int v_sz = e_list->e.at(0).at(c_id).at(0).size();
      ang_def_0c_dev = ang_def_0c_dev + 2*v_sz*PI_L - ang_def_0c[c_id];
    }
    std::cout << "Grain " << g_id.at(i) + 1 << " "
              << CoM.at(i) + ang_def_1c_dev + ang_def_0c_dev << std::endl;
    std::cout << "0c_dev " << ang_def_0c_dev - e_c->conn.size()*ANG_EQUI_CONST
              << std::endl;
    std::cout << "1c_dev " << ang_def_1c_dev << std::endl;
    ang_def_dev_0c.at(i) = ang_def_0c_dev;
    ang_def_dev_1c.at(i) = ang_def_1c_dev;
  }
  CoM_avg = CoM_avg/g_id.size();
  delete e_c;
}

// A wrapper for ext_MS in style of vd_extractor. 
// TODO incorporate ext_MS into this.
void vd_ext_GAUSS_g::parse_opt(std::vector<std::string> &opts) {
  OUTPUTGAUSS = opts.at(1);

  bc_type = conv_str2bc(opts.at(2));

  std::string temp("");
  std::string::size_type sz;
  temp = opts.at(3);
  while(!temp.empty()) {
    int s3_curr = std::stoi(temp, &sz) - 1;
    g_id.push_back(s3_curr);
    temp = temp.substr(sz);
  }
  CoM.resize(g_id.size());
  ang_def_dev_0c.resize(g_id.size());
  ang_def_dev_1c.resize(g_id.size());
  for(int i = 0; i < CoM.size(); i++) {
    CoM.at(i) = 0.;
    ang_def_dev_0c.at(i) = 0.;
    ang_def_dev_1c.at(i) = 0.;
  }
  collect_ori();
  collect_ang();
}

void vd_ext_GAUSS_g::write_csv() {
  csvfile csvGAUSS(OUTPUTGAUSS);

  csvGAUSS << t_curr;
  // Integral of the GC over the 2-stratum vertices. Also includes exterior if
  // bc_type == VD_BC_TYPE::NEUMANN_EXT
  csvGAUSS << CoM_avg;
  for(int i = 0; i < CoM.size(); i++) {
    csvGAUSS << CoM.at(i);
    csvGAUSS << ang_def_dev_0c.at(i);
    csvGAUSS << ang_def_dev_1c.at(i);
  }

  csvGAUSS << endrow;
}

////////////////////////////////////////////////
// EXT_ROC: Rate of change of length, area, volume of given cells.

vd_ext_ROC::vd_ext_ROC(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                                c_list(0, std::make_pair(1,1)),
                                OUTPUTROC("") {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  parse_opt(opts_in);
  write_csv();
}

void vd_ext_ROC::parse_opt(std::vector<std::string> &opts) {
  OUTPUTROC = opts.at(1);

  std::string temp("");
  std::string::size_type sz;
  c_list.resize(opts.size() - 2);
  for(int i = 2; i < opts.size(); i++) {
    temp = opts.at(i);
    c_list.at(i-2).first = std::stod(temp, &sz);
    temp = temp.substr(sz);
    c_list.at(i-2).second = std::stod(temp, &sz);
  }
}

void vd_ext_ROC::write_csv() {
  csvfile csvROC(OUTPUTROC);

  csvROC << t_curr;
  for(int i = 0; i < c_list.size(); i++) {
    int c_dim = c_list.at(i).first;
    int c_id = c_list.at(i).second;
    if(c_dim > 0 and e_list->e.at(c_dim).at(c_id-1).at(c_dim).size() > 0)
      csvROC << f_calc->calc_roc(m, &e_list->e.at(c_dim).at(c_id-1).at(c_dim));
    else
      csvROC << 0.;
  }
  csvROC << endrow;
}

////////////////////////////////////////////////
// EXT_COUNT: Rate of change of length, area, volume of given cells.

vd_ext_COUNT::vd_ext_COUNT(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                                c_list(0, std::vector<int>(0)),
                                OUTPUTCOUNT("") {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  parse_opt(opts_in);
  write_csv();
}

void vd_ext_COUNT::parse_opt(std::vector<std::string> &opts) {
  OUTPUTCOUNT = opts.at(1);

  std::string temp("");
  std::string::size_type sz;
  c_list.resize(opts.size() - 2);
  for(int i = 2; i < opts.size(); i++) {
    c_list.at(i-2).resize(3);
    temp = opts.at(i);
    c_list.at(i-2).at(0) = std::stod(temp, &sz);
    temp = temp.substr(sz);
    c_list.at(i-2).at(1) = std::stod(temp, &sz);
    temp = temp.substr(sz);
    c_list.at(i-2).at(2) = std::stod(temp, &sz);
  }
}

void vd_ext_COUNT::write_csv() {
  csvfile csvCOUNT(OUTPUTCOUNT);

  csvCOUNT << t_curr;
  for(int i = 0; i < c_list.size(); i++) {
    int e_dim = c_list.at(i).at(0);
    int c_dim = c_list.at(i).at(1);
    int c_id = c_list.at(i).at(2)-1;
    if(c_id == -2) {
      int count = 0;
      for(int i = 0; i < e_list->e.at(c_dim).size(); i++) {
        if(!cb->is_free(c_dim, i)) {
          count = count + e_list->e.at(c_dim).at(i).at(e_dim).size();
        }
      }
      csvCOUNT << count;
    }
    else {
      if(!cb->is_free(c_dim, c_id))
        csvCOUNT << e_list->e.at(c_dim).at(c_id).at(e_dim).size();
      else
        csvCOUNT << 0.;
    }
  }
  csvCOUNT << endrow;
}

////////////////////////////////////////////////
// EXT_ANGLE1

vd_ext_ANGLE1::vd_ext_ANGLE1(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                                angle(0.), 
                                OUTPUTANGLE("") {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  parse_opt(opts_in);
  write_csv();
}

// A wrapper for ext_MS in style of vd_extractor. 
// TODO incorporate ext_MS into this.
void vd_ext_ANGLE1::parse_opt(std::vector<std::string> &opts) {
  OUTPUTANGLE = opts.at(1);

  s0 = std::stoi(opts.at(2));
  s10 = std::stoi(opts.at(3));
  s11 = std::stoi(opts.at(4));

  assert(e_list->e.at(0).at(s0-1).at(0).size() == 1);
  ent_conn e_s;
  cb->get_conn_gmi(1, s10, &e_s);
  assert(e_s.chk_ent(s0));
  cb->get_conn_gmi(1, s11, &e_s);
  assert(e_s.chk_ent(s0));

  apf::MeshEntity* v = e_list->e.at(0).at(s0-1).at(0).at(0);
  apf::MeshEntity* e_10;
  apf::MeshEntity* e_11;
  std::vector<apf::MeshEntity*> edges(0);
  vd_set_up(m, v, &edges);

  bool s10_found = false;
  bool s11_found = false;
  for(int i = 0; i < edges.size(); i++) {
    apf::ModelEntity* mdl = m->toModel(edges.at(i));
    if(m->getModelType(mdl) == 1) {
      if(m->getModelTag(mdl) == s10) {
        assert(!s10_found);
        s10_found = true;
        e_10 = edges.at(i);
      }
      else if(m->getModelTag(mdl) == s11) {
        assert(!s11_found);
        s11_found = true;
        e_11 = edges.at(i);
      }
    }
  }
  assert(s10_found and s11_found);

  apf::Vector3 pos_s0(0,0,0);
  apf::Vector3 dir_s10(0,0,0);
  apf::Vector3 dir_s11(0,0,0);

  pos_s0 = vd_get_pos(m, v);
  dir_s10 = vd_get_pos(m, e_10);
  dir_s11 = vd_get_pos(m, e_11);
  dir_s10 = norm_0(dir_s10 - pos_s0);
  dir_s11 = norm_0(dir_s11 - pos_s0);

  angle = vd_inner_angle(dir_s10, dir_s11);
}

void vd_ext_ANGLE1::write_csv() {
  csvfile csvANGLE(OUTPUTANGLE);

  csvANGLE << t_curr;
  csvANGLE << angle;
  csvANGLE << endrow;
}

////////////////////////////////////////////////
// EXT_ANGLE2

vd_ext_ANGLE2::vd_ext_ANGLE2(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                                angle(0.), 
                                OUTPUTANGLE("") {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  parse_opt(opts_in);
  write_csv();
}

// A wrapper for ext_MS in style of vd_extractor. 
// TODO incorporate ext_MS into this.
void vd_ext_ANGLE2::parse_opt(std::vector<std::string> &opts) {
  OUTPUTANGLE = opts.at(1);

  s1 = std::stoi(opts.at(2));
  s3 = std::stoi(opts.at(3));
  s20 = std::stoi(opts.at(4));
  s21 = std::stoi(opts.at(5));

  assert(e_list->e.at(1).at(s1-1).at(1).size() > 0);
  ent_conn e_s;
  cb->get_conn_dim_gmi(1, 3, s3, &e_s);
  assert(e_s.chk_ent(s1));
  cb->get_conn_dim_gmi(2, 3, s3, &e_s);
  assert(e_s.chk_ent(s20));
  assert(e_s.chk_ent(s21));

  apf::Vector3 pos_line(0,0,0);
  apf::Vector3 pos_t_20(0,0,0);
  apf::Vector3 pos_t_21(0,0,0);
  apf::Vector3 norm0(0,0,0);
  apf::Vector3 norm1(0,0,0);

  std::vector<apf::MeshEntity*> * edges = &e_list->e.at(1).at(s1-1).at(1);
  std::vector<apf::MeshEntity*> tris(0);

  angle = 0.;
  double w_tot = 0.;

  for(int i = 0; i < edges->size(); i++) {
    apf::MeshEntity* e = edges->at(i);
    apf::MeshEntity* t_20;
    apf::MeshEntity* t_21;

    double w_curr = vd_meas_ent(m, e);
    w_tot = w_tot + w_curr;

    vd_set_up(m, e, &tris);
    bool s20_found = false;
    bool s21_found = false;
    for(int j = 0; j < tris.size(); j++) {
      apf::ModelEntity* mdl = m->toModel(tris.at(j));
      if(m->getModelType(mdl) == 2) {
        if(m->getModelTag(mdl) == s20) {
          assert(!s20_found);
          s20_found = true;
          t_20 = tris.at(j);
        }
        else if(m->getModelTag(mdl) == s21) {
          assert(!s21_found);
          s21_found = true;
          t_21 = tris.at(j);
        }
      }
    }
    assert(s20_found and s21_found);
    pos_line =  vd_get_pos(m, e);
    pos_t_20 = vd_get_pos(m, t_20) - pos_line;
    pos_t_21 = vd_get_pos(m, t_21) - pos_line; 
    norm0 = norm_0(vd_area_out(m, t_20, s3));
    norm1 = norm_0(vd_area_out(m, t_21, s3));
    // TODO For some reason vd_ext_angle_n works for MS calculation but it seems 
    // to be wrong. Adding here instead of subtracting...
    // TODO it seems the vd_area_out was actually pointing inward. Reverting to
    // addition again.
    double angle_curr = PI_L - vd_ext_angle_n(pos_t_20, pos_t_21, norm0, norm1);
    //double angle_curr = PI_L + vd_ext_angle_n(pos_t_20, pos_t_21, norm0, norm1);
    angle = angle + w_curr*angle_curr;
  }

  angle = angle/w_tot;
}

void vd_ext_ANGLE2::write_csv() {
  csvfile csvANGLE(OUTPUTANGLE);

  csvANGLE << t_curr;
  csvANGLE << angle;
  csvANGLE << endrow;
}


////////////////////////////////////////////////

double ext_ppa(apf::Mesh2* m, vd_entlist* e_list, std::vector<int>& s2_id) {
  double A_tot = 0.;
  for(int i = 0; i < s2_id.size(); i++) {
    for(int j = 0; j < e_list->e.at(2).at(s2_id.at(i)).at(2).size(); j++) {
      A_tot = A_tot + vd_meas_ent(m, e_list->e.at(2).at(s2_id.at(i)).at(2).at(j));
    }
  }
  int v_count = 0;
  for(int i = 0; i < s2_id.size(); i++) {
    v_count = v_count + e_list->e.at(2).at(s2_id.at(i)).at(0).size();
  }
  double ppa = v_count/A_tot;
  if(std::isnan(ppa))
    return 0;
  return ppa;
}

////////////////////////////////////////////////
void write_off(const char* fname, apf::Mesh2* m, vd_entlist* e_list, std::vector<int>& s2_id) {

  std::vector<apf::MeshEntity*> t_2c(0);
  for(int i = 0; i < s2_id.size(); i++)
    t_2c.reserve(t_2c.size() + e_list->e.at(2).at(s2_id.at(i)).at(2).size());
  for(int i = 0; i < s2_id.size(); i++) {
    for(int j = 0; j < e_list->e.at(2).at(s2_id.at(i)).at(2).size(); j++) {
      t_2c.push_back(e_list->e.at(2).at(s2_id.at(i)).at(2).at(j));
    }
  }

  std::vector<apf::MeshEntity*> v_2c(0);
  std::vector<apf::MeshEntity*> temp(0);
  vd_set_down(m, &t_2c, &temp);
  vd_set_down(m, &temp, &v_2c);

  std::ofstream fs_;
  fs_.exceptions(std::ios::failbit | std::ios::badbit);
  //fs_.open(fname, std::ios_base::app);
  fs_.open(fname, std::ofstream::out | std::ofstream::trunc);

  apf::Vector3 pt(0,0,0);

  fs_ << "OFF\r\n";
  fs_ << v_2c.size() << " " << t_2c.size() << " 0\r\n";
  fs_ << std::setprecision(6);
  for(int i = 0; i < v_2c.size(); i++) {
    m->getPoint(v_2c.at(i), 0, pt);
    fs_ << std::fixed << pt[0] << " " << std::fixed << pt[1] << " " 
        << std::fixed << pt[2] << "\r\n";
  }

  apf::Downward d_v;
  apf::Downward d_e;
  std::map<apf::MeshEntity*, int> e_t_nbr{};
  std::map<apf::MeshEntity*, apf::MeshEntity*> e_t_map1{};
  std::map<apf::MeshEntity*, apf::MeshEntity*> e_t_map2{};

  bool disc_like = true;
  for(int i = 0; i < t_2c.size(); i++) {
    m->getDownward(t_2c.at(i), 1, d_e);
    for(int j = 0; j < 3; j++) {
      e_t_nbr[d_e[j]] = e_t_nbr[d_e[j]] + 1;
      if(!e_t_map1[d_e[j]])
        e_t_map1[d_e[j]] = t_2c.at(i);
      else if(!e_t_map2[d_e[j]])
        e_t_map2[d_e[j]] = t_2c.at(i);
      else {
        disc_like = false;
        i = t_2c.size();
        j = 3;
      }
      if(e_t_nbr[d_e[j]] > 2) {
        disc_like = false;
        i = t_2c.size();
        j = 3;
      }
    }
  }

  if(!disc_like)
    return;
  // Store the current ordering of vertices. 
  // If an orientation of any vertex pair is used, the ordering of the vertices
  // should be flipped. For a surface composed of edge sharing triangles, 
  // this should yield a consistent arrangement. 
  std::map<std::pair<int,int>, bool> e_used{};
  std::vector<int> v_id(3);
  std::map<apf::MeshEntity*, bool> t_burned{};
  std::vector<std::vector<apf::MeshEntity*> > burning(2, 
                                      std::vector<apf::MeshEntity*> (0));
  burning.at(0).reserve(t_2c.size());
  burning.at(1).reserve(t_2c.size());
  int id_curr = 0;
  int id_next = 1;
  for(int i = 0; i < t_2c.size(); i++) {
    if(!t_burned[t_2c.at(i)]) {
      burning.at(id_curr).clear();
      burning.at(id_curr).push_back(t_2c.at(i));
      t_burned[burning.at(id_curr).back()] = true;

      while(!burning.at(id_curr).empty()) {
        apf::MeshEntity* t_curr = burning.at(id_curr).back();
        m->getDownward(t_curr, 1, d_e);

        // Add current triangle:
        m->getDownward(t_curr, 0, d_v);
        fs_ << 3 << " ";
        for(int j = 0; j < 3; j++) {
          int v1 = findIn(&v_2c, 3, d_v[j]);
          v_id.at(j) = v1;
        }
        bool found = false;
        for(int j = 0; j < 3; j++) {
          int v2 = v_id.at((j + 1) % 3);
          std::pair<int, int> e_id = std::make_pair(v_id.at(j), v2);
          if(e_used[e_id]) {
            found = true;
            j = 3;
          }
        }
        if(found) {
          int temp = v_id.at(0);
          v_id.at(0) = v_id.at(2);
          v_id.at(2) = temp;
        }
        for(int j = 0; j < 3; j++) {
          int v2 = v_id.at((j + 1) % 3);
          std::pair<int, int> e_id = std::make_pair(v_id.at(j), v2);
          e_used[e_id] = true;
          fs_ << v_id.at(j) << " ";
        }
        fs_ << "\r\n";

        // Collect the next triangles:
        for(int j = 0; j < 3; j++) {
          apf::MeshEntity* t_next = NULL;
          if(e_t_map1[d_e[j]] != t_curr) {
            assert(e_t_map2[d_e[j]] == t_curr);
            t_next = e_t_map1[d_e[j]];
          }
          else {
            assert(e_t_map2[d_e[j]] != NULL);
            t_next = e_t_map2[d_e[j]];
          }
          if(t_next != NULL and !t_burned[t_next]) {
            burning.at(id_next).push_back(t_next);
            t_burned[t_next] = true;
          }
        }
        burning.at(id_curr).pop_back();
        if(burning.at(id_curr).empty()) {
          id_next = id_curr;
          id_curr = (id_next + 1) % 2;
        }
      }
    }
  }

  fs_.flush();
  fs_.close();
}

vd_ext_HAUSDORFF::vd_ext_HAUSDORFF(apf::Mesh2* m_in, cell_base* cb_in, 
                                   vd_entlist* e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                                s2_id(0), 
                                REF_OFF(""),
                                CURR_OFF(""),
                                OUTPUTHAUSDORFF("")  {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  parse_opt(opts_in);
  write_csv();
}

// A wrapper for ext_MS in style of vd_extractor. 
// TODO incorporate ext_MS into this.
void vd_ext_HAUSDORFF::parse_opt(std::vector<std::string> &opts) {
  OUTPUTHAUSDORFF = opts.at(1);
  REF_OFF = opts.at(2);
  CURR_OFF = opts.at(3);
  s2_id.reserve(opts.at(4).size());

  std::string temp("");
  std::string::size_type sz;
  temp = opts.at(4);
  while(!temp.empty()) {
    int s2_curr = std::stoi(temp, &sz);
    s2_id.push_back(s2_curr);
    temp = temp.substr(sz);
  }

  nppa = std::stod(opts.at(5));

  write_off(CURR_OFF.c_str(), m, e_list, s2_id);
}

void vd_ext_HAUSDORFF::write_csv() {
  std::stringstream ss;
  ss.precision(2);
  ss << "./hausdorff " << REF_OFF << " " << CURR_OFF << " " << nppa << " " 
     << std::setprecision(6) << std::fixed << t_curr
     << " " << OUTPUTHAUSDORFF;

  std::string tmp = ss.str();
  std::system(tmp.c_str());
}

////////////////////////////////////////////////

vd_ext_CRV::vd_ext_CRV(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                               rad_field(NULL), curv_field(NULL), 
                               vel_r_field(NULL), vel_field(NULL), 
                               area_field(NULL),
                               t_nbr_field(NULL),
                               OUTPUTCRV(""),
                               cell_2(-1), cell_3(-1), o(0,0,0)
                                                {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  rad_field = vd_att_vs_field(m, "rad");
  curv_field = vd_att_vs_field(m, "curv");
  vel_r_field = vd_att_vs_field(m, "vel");
  vel_field = m->findField("velocity_field");
  area_field = vd_att_vs_field(m, "area");
  t_nbr_field = vd_att_vs_field(m, "t_nbr");

  if(!vel_field) {
    f_calc->vd_att_fields(m);
    f_calc->vd_calc_vel(m);
    vel_field = m->findField("velocity_field");
    assert(vel_field);
  }
  parse_opt(opts_in);
  write_csv();

  apf::destroyField(rad_field);
  apf::destroyField(curv_field);
  apf::destroyField(vel_r_field);
  apf::destroyField(t_nbr_field);
  apf::destroyField(area_field);
}

// Parse the options for extracting curvature, radial position and velocity 
// associated with a 2stratum and an interior 3stratum.
// opts: EXT_TYPE::CRV; cell_2; cell_3; o[0] 0[1] o[2]; 
void vd_ext_CRV::parse_opt(std::vector<std::string> &opts) {
  assert(opts.size() == 5);
  OUTPUTCRV = opts.at(1);
  cell_2 = std::stoi(opts.at(2)) - 1;
  cell_3 = std::stoi(opts.at(3));
  std::string temp("");
  std::string::size_type sz;
  temp = opts.at(4);
  o[0] = std::stod(temp, &sz);
  temp = temp.substr(sz);
  o[1] = std::stod(temp, &sz);
  temp = temp.substr(sz);
  o[2] = std::stod(temp, &sz);

  assert(rad_field);
  assert(curv_field);
  assert(vel_r_field);
  assert(vel_field);
  assert(area_field);

  std::vector<apf::MeshEntity*>* v_2c = &e_list->e.at(2).at(cell_2).at(0);
  std::vector<apf::MeshEntity*>* t_2c = &e_list->e.at(2).at(cell_2).at(2);

  // In order to calculate area weighted quantities defined at vertices, 
  // count number of adjacent triangles for the boundary vertices.
  apf::Downward d_v;
  for(int i = 0; i < t_2c->size(); i++) {
    m->getDownward(t_2c->at(i), 0, d_v);
    for(int j = 0; j < 3; j++) {
      int nbr_curr = apf::getScalar(t_nbr_field, d_v[j], 0);
      apf::setScalar(t_nbr_field, d_v[j], 0, nbr_curr + 1);
    }
  }
  apf::Vector3 rad(0,0,0);
  apf::Vector3 vel(0,0,0);

  for(int i = 0; i < v_2c->size(); i++) {
    double curv_curr = vd_mean_curv(m, v_2c->at(i), cell_3);
    apf::setScalar(curv_field, v_2c->at(i), 0, curv_curr);
    m->getPoint(v_2c->at(i), 0, rad);
    rad = rad - o;
    double r_curr = rad.getLength();
    apf::setScalar(rad_field, v_2c->at(i), 0, r_curr);
    apf::getVector(vel_field, v_2c->at(i), 0, vel);
    double v_curr = rad*vel;
    apf::setScalar(vel_r_field, v_2c->at(i), 0, v_curr);
    apf::setScalar(area_field, v_2c->at(i), 0, 0);
  }
  struct ent_conn e_1c;
  cb->get_conn(2, cell_2, &e_1c);
  for(int i = 0; i < e_1c.conn.size(); i++) {
    int cell_1 = e_1c.conn.at(i);
    std::vector<apf::MeshEntity*>* v_1c = &e_list->e.at(1).at(cell_1).at(0); 
    for(int j = 0; j < v_1c->size(); j++) {
      double curv_curr = vd_mean_curv(m, v_1c->at(j), cell_3);
      apf::setScalar(curv_field, v_1c->at(j), 0, curv_curr);
      m->getPoint(v_1c->at(j), 0, rad);
      rad = rad - o;
      double r_curr = rad.getLength();
      apf::setScalar(rad_field, v_1c->at(j), 0, r_curr);
      apf::getVector(vel_field, v_1c->at(j), 0, vel);
      double v_curr = rad*vel;
      apf::setScalar(vel_r_field, v_1c->at(j), 0, v_curr);
      apf::setScalar(area_field, v_1c->at(j), 0, 0);
    }
  }
}

void vd_ext_CRV::write_csv() {

  std::vector<apf::MeshEntity*>* t_2c = &e_list->e.at(2).at(cell_2).at(2);
  std::vector<apf::MeshEntity*>* v_2c = &e_list->e.at(2).at(cell_2).at(0);
  if(t_2c->size() == 0)
    return;
  apf::Downward d_v;
  double area_tot = 0;
  double curv_w = 0;
  double r_w = 0;
  double vel_w = 0;
  for(int i = 0; i < t_2c->size(); i++) {
    double a_curr = vd_area_out_n(m, t_2c->at(i)).getLength();
    area_tot = area_tot + a_curr;
    m->getDownward(t_2c->at(i), 0, d_v);
    for(int j = 0; j < 3; j++) {
      //int nbr_curr = apf::getScalar(t_nbr_field, d_v[j], 0);
      double area_vert = apf::getScalar(area_field, d_v[j], 0);
      apf::setScalar(area_field, d_v[j], 0, area_vert + a_curr);
      double v_curr = apf::getScalar(vel_r_field, d_v[j], 0);
      //double curv_curr = apf::getScalar(curv_field, d_v[j], 0);
      double rad_curr = apf::getScalar(rad_field, d_v[j], 0);
      vel_w = vel_w + a_curr*v_curr;
      r_w = r_w + a_curr*rad_curr;
    }
  }
  r_w = r_w/area_tot/3;
  vel_w = vel_w/area_tot/3;

  double area_total_2c = 0.;
  for(int i = 0; i < v_2c->size(); i++) {
    double a_curr = apf::getScalar(area_field, v_2c->at(i), 0);
    area_total_2c = area_total_2c + a_curr;
    double curv_curr = apf::getScalar(curv_field, v_2c->at(i), 0);
    curv_w = curv_w + a_curr*curv_curr;
  }
  if(v_2c->size() == 0)
    curv_w = 0.;
  else
    curv_w = curv_w/area_total_2c/3;

  csvfile csvCRV(OUTPUTCRV);

  double ROC = f_calc->calc_roc(m, t_2c);

  csvCRV << t_curr;
  csvCRV << curv_w << r_w << vel_w << area_tot << ROC;
  csvCRV << endrow;
}

vd_ext_POS_D::vd_ext_POS_D(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                               pos_field(NULL),
                               e_nbr_field(NULL),
                               OUTPUTPOS(""),
                               cell_dim(-1), cell_id(-1), o(0,0,0), dir(0,0,0)
                                                {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  pos_field = vd_att_vs_field(m, "pos");
  e_nbr_field = vd_att_vs_field(m, "e_nbr");

  parse_opt(opts_in);
  write_csv();

  apf::destroyField(pos_field);
  apf::destroyField(e_nbr_field);
}

// Parse the options for extracting curvature, radial position and velocity 
// associated with a 2stratum and an interior 3stratum.
// opts: EXT_TYPE::CRV; cell_2; cell_3; o[0] 0[1] o[2]; 
void vd_ext_POS_D::parse_opt(std::vector<std::string> &opts) {
  assert(opts.size() == 6);
  OUTPUTPOS = opts.at(1);
  cell_dim = std::stoi(opts.at(2));
  cell_id = std::stoi(opts.at(3)) - 1;
  std::string temp("");
  std::string::size_type sz;
  temp = opts.at(4);
  o[0] = std::stod(temp, &sz);
  temp = temp.substr(sz);
  o[1] = std::stod(temp, &sz);
  temp = temp.substr(sz);
  o[2] = std::stod(temp, &sz);
  temp = opts.at(5);
  dir[0] = std::stod(temp, &sz);
  temp = temp.substr(sz);
  dir[1] = std::stod(temp, &sz);
  temp = temp.substr(sz);
  dir[2] = std::stod(temp, &sz);
  dir = norm_0(dir);

  assert(pos_field);

  if(cell_dim > 0) {
    std::vector<apf::MeshEntity*>* v_c = &e_list->e.at(cell_dim).at(cell_id).at(0);
    std::vector<apf::MeshEntity*>* e_c = &e_list->e.at(cell_dim).at(cell_id).at(cell_dim);

    apf::Vector3 pos(0,0,0);
    // In order to calculate area weighted quantities defined at vertices, 
    // count number of adjacent triangles for the boundary vertices.
    apf::Downward d_v;
    for(int i = 0; i < e_c->size(); i++) {
      m->getDownward(e_c->at(i), 0, d_v);
      for(int j = 0; j < cell_dim + 1; j++) {
        m->getPoint(d_v[j], 0, pos);
        pos = pos - o;
        double r_curr = pos*dir;
        apf::setScalar(pos_field, d_v[j], 0, r_curr);
      }
    }
  }
}

void vd_ext_POS_D::write_csv() {

  std::vector<apf::MeshEntity*>* e_c = &e_list->e.at(cell_dim).at(cell_id).at(cell_dim);
  if(cell_dim > 0) {
    apf::Downward d_v;
    double w_tot = 0;
    double r_w = 0;
    int sz = cell_dim + 1;
    for(int i = 0; i < e_c->size(); i++) {
      double w_curr = vd_meas_ent(m, e_c->at(i));
      w_tot = w_tot + w_curr;
      m->getDownward(e_c->at(i), 0, d_v);

      double r_sum = 0;
      for(int j = 0; j < sz; j++) {
        double rad_curr = apf::getScalar(pos_field, d_v[j], 0);
        r_sum = r_sum + rad_curr;
      }
      r_w = r_w + w_curr*r_sum/sz;
    }
    r_w = r_w/w_tot;

    csvfile csvPOS(OUTPUTPOS);

    csvPOS << t_curr;
    csvPOS << r_w;
    csvPOS << endrow;
  }
  else {
    assert(e_c->size() == 1);
    apf::Vector3 pos(0,0,0);
    m->getPoint(e_c->at(0), 0, pos);
    double r_curr = pos*dir;
    csvfile csvPOS(OUTPUTPOS);

    csvPOS << t_curr;
    csvPOS << r_curr;
    csvPOS << endrow;
  }
}

vd_ext_VEL_D::vd_ext_VEL_D(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                               vel_s_field(NULL), vel_field(NULL), 
                               e_nbr_field(NULL),
                               OUTPUTVEL(""),
                               cell_dim(-1), cell_id(-1), dir(0,0,0)
                                                {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;
  e_list->refresh();

  t_curr = t_in;

  vel_s_field = vd_att_vs_field(m, "vel");
  vel_field = m->findField("velocity_field");
  e_nbr_field = vd_att_vs_field(m, "e_nbr");

  if(!vel_field) {
    f_calc->vd_att_fields(m);
    f_calc->vd_calc_vel(m);
    vel_field = m->findField("velocity_field");
    assert(vel_field);
  }
  parse_opt(opts_in);
  write_csv();

  apf::destroyField(vel_s_field);
  apf::destroyField(e_nbr_field);
}

// Parse the options for extracting curvature, radial position and velocity 
// associated with a 2stratum and an interior 3stratum.
// opts: EXT_TYPE::CRV; cell_2; cell_3; o[0] 0[1] o[2]; 
void vd_ext_VEL_D::parse_opt(std::vector<std::string> &opts) {
  assert(opts.size() == 5);
  OUTPUTVEL = opts.at(1);
  cell_dim = std::stoi(opts.at(2));
  cell_id = std::stoi(opts.at(3)) - 1;
  std::string temp("");
  std::string::size_type sz;
  temp = opts.at(4);
  dir[0] = std::stod(temp, &sz);
  temp = temp.substr(sz);
  dir[1] = std::stod(temp, &sz);
  temp = temp.substr(sz);
  dir[2] = std::stod(temp, &sz);
  dir = norm_0(dir);

  assert(vel_s_field);
  assert(vel_field);
  if(cell_dim > 0) {
    std::vector<apf::MeshEntity*>* v_c = &e_list->e.at(cell_dim).at(cell_id).at(0);
    std::vector<apf::MeshEntity*>* e_c = &e_list->e.at(cell_dim).at(cell_id).at(cell_dim);

    // In order to calculate area weighted quantities defined at vertices, 
    // count number of adjacent triangles for the boundary vertices.
    apf::Downward d_v;

    apf::Vector3 vel(0,0,0);

    for(int i = 0; i < v_c->size(); i++) {
      apf::getVector(vel_field, v_c->at(i), 0, vel);
      double v_curr = vel*dir;
      apf::setScalar(vel_s_field, v_c->at(i), 0, v_curr);
    }
  }
}

void vd_ext_VEL_D::write_csv() {

  std::vector<apf::MeshEntity*>* e_c = &e_list->e.at(cell_dim).at(cell_id).at(cell_dim);

  if(cell_dim > 0) {

    apf::Downward d_v;
    double w_tot = 0;
    double v_w = 0;
    int sz = cell_dim + 1;
    for(int i = 0; i < e_c->size(); i++) {
      double w_curr = vd_meas_ent(m, e_c->at(i));
      w_tot = w_tot + w_curr;
      m->getDownward(e_c->at(i), 0, d_v);
      double v_sum = 0;
      for(int j = 0; j < sz; j++) {
        double vel_curr = apf::getScalar(vel_s_field, d_v[j], 0);
        v_sum = v_sum + vel_curr;
      }
      v_w = v_w + w_curr*v_sum/sz;
    }
    v_w = v_w/w_tot;

    csvfile csvVEL(OUTPUTVEL);

    csvVEL << t_curr;
    csvVEL << v_w;
    csvVEL << endrow;
  }
  else {
    assert(e_c->size() == 1);
    apf::Vector3 vel(0,0,0);
    apf::getVector(vel_field, e_c->at(0), 0, vel);
    double v_curr = vel*dir;
    csvfile csvVEL(OUTPUTVEL);

    csvVEL << t_curr;
    csvVEL << v_curr;
    csvVEL << endrow;
  }
}

///////////////////
// vd_ext_meas
///////////////////

vd_ext_MEAS::vd_ext_MEAS(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                               OUTPUTMEAS(""),
                               cells(0, std::make_pair(-1,-1)) {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  parse_opt(opts_in);
  write_csv();

}

// Parse the options for extracting curvature, radial position and velocity 
// associated with a 2stratum and an interior 3stratum.
// opts: EXT_TYPE::CRV; cell_2; cell_3; o[0] 0[1] o[2]; 
void vd_ext_MEAS::parse_opt(std::vector<std::string> &opts) {
  assert(opts.size() > 2);
  OUTPUTMEAS = opts.at(1);
  cells.reserve(opts.size() - 2);
  for(int i = 2; i < opts.size(); i++) {
    std::string temp("");
    std::string::size_type sz;
    temp = opts.at(i);
    int c_dim = std::stoi(temp, &sz);
    temp = temp.substr(sz);
    int c_id = std::stoi(temp, &sz);

    cells.push_back(std::make_pair(c_dim, c_id));
  }

}

void vd_ext_MEAS::write_csv() {
  csvfile csvMEAS(OUTPUTMEAS);
  csvMEAS << t_curr;

  for(int i = 0; i < cells.size(); i++) {
    int cell_dim = cells.at(i).first;
    int cell_id = cells.at(i).second - 1;
    std::vector<apf::MeshEntity*>* e_c = &e_list->e.at(cell_dim).at(cell_id).at(cell_dim);
    csvMEAS << vd_meas_set(m, e_c);
  }
  csvMEAS << endrow;
}
///////////////////
// vd_ext_AREA_T
///////////////////

vd_ext_AREA_T::vd_ext_AREA_T(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                               OUTPUTMEAS("") {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  parse_opt(opts_in);
  write_csv();

}

// Parse the options for extracting curvature, radial position and velocity 
// associated with a 2stratum and an interior 3stratum.
// opts: EXT_TYPE::CRV; cell_2; cell_3; o[0] 0[1] o[2]; 
void vd_ext_AREA_T::parse_opt(std::vector<std::string> &opts) {
  OUTPUTMEAS = opts.at(1);
}

void vd_ext_AREA_T::write_csv() {
  csvfile csvMEAS(OUTPUTMEAS);
  csvMEAS << t_curr;

  apf::MeshIterator* it = m->begin(2);
  double area_tot = 0.;
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    apf::ModelEntity* mdl = m->toModel(e);
    int type = m->getModelType(mdl);
    int tag = m->getModelTag(mdl);
    if(type == 2) {
      area_tot = area_tot + vd_meas_ent(m, e);
    }
    //std::cout<< std::endl;
  }
  m->end(it);
  csvMEAS << area_tot;

  csvMEAS << endrow;
}
///////////////////
// vd_ext_pos
///////////////////

vd_ext_POS::vd_ext_POS(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                               e_nbr_field(NULL),
                               OUTPUTPOS(""),
                               cell_dim(-1), cell_id(-1), o(0,0,0)
                                                {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  e_nbr_field = vd_att_vs_field(m, "e_nbr");

  parse_opt(opts_in);
  write_csv();

  apf::destroyField(e_nbr_field);
}

// Parse the options for extracting curvature, radial position and velocity 
// associated with a 2stratum and an interior 3stratum.
// opts: EXT_TYPE::CRV; cell_2; cell_3; o[0] 0[1] o[2]; 
void vd_ext_POS::parse_opt(std::vector<std::string> &opts) {
  assert(opts.size() == 6);
  OUTPUTPOS = opts.at(1);
  cell_dim = std::stoi(opts.at(2));
  cell_id = std::stoi(opts.at(3)) - 1;
  std::string temp("");
  std::string::size_type sz;
  temp = opts.at(4);
  o[0] = std::stod(temp, &sz);
  temp = temp.substr(sz);
  o[1] = std::stod(temp, &sz);
  temp = temp.substr(sz);
  o[2] = std::stod(temp, &sz);

}

void vd_ext_POS::write_csv() {

  std::vector<apf::MeshEntity*>* e_c = &e_list->e.at(cell_dim).at(cell_id).at(cell_dim);
  if(cell_dim > 0) {
    apf::Downward d_v;
    double w_tot = 0;
    double curv_w = 0;
    apf::Vector3 pos_w(0,0,0);
    double vel_w = 0;
    for(int i = 0; i < e_c->size(); i++) {
      double w_curr = vd_meas_ent(m, e_c->at(i));
      w_tot = w_tot + w_curr;
      m->getDownward(e_c->at(i), 0, d_v);
      for(int j = 0; j < cell_dim+1; j++) {
        apf::Vector3 pos_curr(0,0,0);
        m->getPoint(d_v[j], 0, pos_curr);
        pos_curr = pos_curr - o;
        pos_w = pos_w + pos_curr*w_curr/(cell_dim+1);
      }
    }
    pos_w = pos_w/w_tot;

    csvfile csvPOS(OUTPUTPOS);

    csvPOS << t_curr;
    csvPOS << pos_w[0] << pos_w[1] << pos_w[2];
    csvPOS << endrow;
  }
  else {
    assert(e_c->size() == 1);
    apf::Vector3 pos(0,0,0);
    m->getPoint(e_c->at(0), 0, pos);
    pos = pos - o;

    csvfile csvPOS(OUTPUTPOS);

    csvPOS << t_curr;
    csvPOS << pos[0] << pos[1] << pos[2];
    csvPOS << endrow;
  }
}

////////////////////////////
//  vd_ext_vel
////////////////////////////

vd_ext_VEL::vd_ext_VEL(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                               e_list_in, field_calc* f_calc_in, 
                               std::vector<std::string> & opts_in,
                               double t_in) : 
                               vel_field(NULL), 
                               e_nbr_field(NULL),
                               OUTPUTVEL(""),
                               cell_dim(-1), cell_id(-1)
                                                {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  vel_field = m->findField("velocity_field");
  e_nbr_field = vd_att_vs_field(m, "e_nbr");

  if(!vel_field) {
    f_calc->vd_att_fields(m);
    f_calc->vd_calc_vel(m);
    vel_field = m->findField("velocity_field");
    assert(vel_field);
  }
  parse_opt(opts_in);
  write_csv();

  apf::destroyField(e_nbr_field);
}

// Parse the options for extracting curvature, radial position and velocity 
// associated with a 2stratum and an interior 3stratum.
// opts: EXT_TYPE::CRV; cell_2; cell_3; o[0] 0[1] o[2]; 
void vd_ext_VEL::parse_opt(std::vector<std::string> &opts) {
  assert(opts.size() == 5);
  OUTPUTVEL = opts.at(1);
  cell_dim = std::stoi(opts.at(2));
  cell_id = std::stoi(opts.at(3)) - 1;

  assert(vel_field);
}

void vd_ext_VEL::write_csv() {

  std::vector<apf::MeshEntity*>* e_c = &e_list->e.at(cell_dim).at(cell_id).at(cell_dim);

  apf::Vector3 vel(0,0,0);
  if(cell_dim > 0) {
    apf::Vector3 vel_curr(0,0,0);

    apf::Downward d_v;
    double w_tot = 0;
    double v_w = 0;
    for(int i = 0; i < e_c->size(); i++) {
      double w_curr = vd_meas_ent(m, e_c->at(i));
      w_tot = w_tot + w_curr;
      m->getDownward(e_c->at(i), 0, d_v);
      for(int j = 0; j < cell_dim+1; j++) {
        apf::getVector(vel_field, e_c->at(0), 0, vel_curr);
        vel = vel + vel_curr*w_curr/(cell_dim+1);
      }
    }
    v_w = v_w/w_tot;

    csvfile csvVEL(OUTPUTVEL);

    csvVEL << t_curr;
    csvVEL << v_w;
    csvVEL << endrow;
  }
  else {
    assert(e_c->size() == 1);
    apf::getVector(vel_field, e_c->at(0), 0, vel);
    csvfile csvVEL(OUTPUTVEL);

    csvVEL << t_curr;
    csvVEL << vel[0] << vel[1] << vel[2];
    csvVEL << endrow;
  }
}

////////////////////////////
//  vd_ext_inf
////////////////////////////

void vd_ext_inf::write_csv() {
  for(int i = 0; i < opts.size(); i++) {
    assert(opts.at(i).size() > 1);
    EXT_TYPE et_curr = conv_str2ext(opts.at(i).at(0));
    assert(et_curr < EXT_TYPE::END);
    if(et_curr == EXT_TYPE::MS) {
      vd_ext_MS ext_MS(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::VEL) {
    }
    else if(et_curr == EXT_TYPE::POS) {
    }
    else if(et_curr == EXT_TYPE::VEL_DIR) {
      vd_ext_VEL_D ext_VELDIR(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::POS_DIR) {
      vd_ext_POS_D ext_POSDIR(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::CRV) {
      vd_ext_CRV ext_CRV(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::MEAS) {
      vd_ext_MEAS ext_MEAS(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::AREA_T) {
      vd_ext_AREA_T ext_AREA(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::GAUSS) {
      vd_ext_GAUSS ext_GAUSS(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::GAUSS_G) {
      vd_ext_GAUSS_g ext_GAUSS_g(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::ROC) {
      vd_ext_ROC ext_ROC(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::COUNT) {
      vd_ext_COUNT ext_COUNT(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::ANGLE1) {
      vd_ext_ANGLE1 ext_ANGLE(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::ANGLE2) {
      vd_ext_ANGLE2 ext_ANGLE(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
    else if(et_curr == EXT_TYPE::HAUSDORFF) {
      vd_ext_HAUSDORFF ext_HAUSDORFF(m, cb, e_list, f_calc, opts.at(i), t_curr);
    }
  }
}

vd_ext_inf::vd_ext_inf(apf::Mesh2* m_in, cell_base* cb_in, vd_entlist* 
                                          e_list_in, field_calc* f_calc_in, 
                                          std::string &opts_file,
                                          double t_in) : 
                                          m(NULL), cb(NULL),
                                          e_list(NULL), f_calc(NULL), t_curr(0),
                                     opts(0, std::vector<std::string> (0,"")) {
  m = m_in;
  cb = cb_in;
  e_list = e_list_in;
  f_calc = f_calc_in;

  t_curr = t_in;

  ReadNames(opts_file, ";", opts);
  write_csv();

}

vd_ext_inf::~vd_ext_inf() {
}

/*
// Process the information extracted from the mesh.
// Uses extractors. Can be used to derive information such as the velocity at
// the the point satisfying a condition (e.g. at the maximum absolute mean curvature, center of the cell, or the mean velocity of the cell).
// Integral of the gaussian curvature, integral of total cell length, area, volume... So given a field (position, surface energy density) integrate using finite element basis functions (apf interface)
class vd_processor {
  public:
    void proc_field() {
      ext_mesh e_msh;
      e_msh.calc_field();

      // Average can have an integral description, as well
      if(Integrate) {
      }
    }

    vd_processor();
    ~vd_processor();
};

class vd_field_int : public Integrator
{
  public:
    double field_measure(MeshElement* e, char* field_name);
    apf::Vector3 field_measure(MeshElement* e, char* field_name);
apf::getValueType(field_curr);

    vd_field_int(int order):Integrator(order),m(0) {}
    void atPoint(Vector3 const&, double w, double f_val) {
      m += w*f_val;
    }
    void atPoint(Vector3 const&, double w, apf::Vector3 f_val) {
      v += w*f_val;
    }
    double m;
    apf::Vector3 v;
};

// Adapted from Integrator.process and measurer
double vd_field_int::field_measure_db(MeshElement* e, char* field_name) {
  apf::Field* field_curr = e_lens.m->findField(field_name);
  assert(field_curr);
  assert(apf::getValueType(field_curr) == apf::SCALAR);
  w = 0;

  this->inElement(e);
  int np = countIntPoints(e,this->order);

  for (int p=0; p < np; ++p)  {
    ipnode = p;
    apf::Vector3 point;
    getIntPoint(e,this->order,p,point);
    double w = getIntWeight(e,this->order,p);
    double f_val = apf::getScalar(field_curr, e, 0);

    this->atPoint(point,w,f_val);
  }
  this->outElement();
  return m;
}

// Adapted from Integrator.process and measurer
apf::Vector3 vd_field_int::field_measure_v(MeshElement* e, char* field_name) {
  apf::Field* field_curr = e_lens.m->findField(field_name);
  assert(field_curr);
  assert(apf::getValueType(field_curr) == apf::VECTOR);

  this->inElement(e);
  int np = countIntPoints(e,this->order);

  double z[3] = {0,0,0};
  v.fromArray(z);

  for (int p=0; p < np; ++p)  {
    ipnode = p;
    apf::Vector3 point;
    getIntPoint(e,this->order,p,point);
    double w = getIntWeight(e,this->order,p);
    apf::Vector3 f_val;
    apf::getVector(field_curr, e, 0, f_val);

    this->atPoint(point,w,f_val);
  }
  this->outElement();
  return v;
}
*/
