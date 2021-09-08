#include <iostream>
#include <algorithm>    // std::find
#include <assert.h>

#include "gmi_base.h"
#include <gmi_mesh.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>

#include "topo_topo.h"
#include "topo_extinfo.h"
// ---------------------------------------------------------
// ent_conn is a container for storing cell information.
// ---------------------------------------------------------

// ---------------------------------------------------------
// Connection list operations. Check, add to the end, remove tag, remove first
// ---------------------------------------------------------

// Remove list of entities from the conn list. 
void rem_ents(std::vector<int>* conn, std::vector<int>* tags) {
  if(conn->size() == 0 or tags->size() == 0)
    return;
  int rmv_nbr = 0;
  for(int i = 0; i < conn->size() - rmv_nbr; i++) {
    std::vector<int>::iterator it = std::find(tags->begin(), 
                                              tags->end(), conn->at(i));
    if(it != tags->end()) {
      int last_id = conn->size() - 1 - rmv_nbr;
      int temp = conn->at(last_id);
      conn->at(last_id) = conn->at(i);
      conn->at(i) = temp;
      rmv_nbr = rmv_nbr + 1;
      i = i - 1;
    }
  }
  conn->resize(conn->size() - rmv_nbr);

}

bool compc_elem(const c_elem &a, const c_elem &b) { 
  if (a.first.first == b.first.first)
    return a.first.second < b.first.second;
  return a.first.first < b.first.first; 
}

void sort_celem(std::vector<c_elem>* c_set, int left, int right) {
  if(right == -1 and left == -1) {
    left = 0;
    right = c_set->size()-1;
  }
  //std::cout << "Sorting edge set " << e_set << std::endl;
  //std::cout << left << " " << right << std::endl;

  if (left < right) {
    int part = partition_celem(c_set, left, right);
    sort_celem(c_set, left, part - 1);
    sort_celem(c_set, part + 1, right);
  }
}
//https://en.wikipedia.org/wiki/Quicksort
//Function to determine the partitions
// partitions the array and returns the middle subscript
int partition_celem(std::vector<c_elem>* c_set, int left, int right) {
  c_elem pivot = c_set->at(right);
  // move the mid point value to the front.
  int i = left-1;
  int j = left;
  for (; j < right; j++) {
    if(compc_elem(c_set->at(j), pivot)) {
      i++;
      std::swap(c_set->at(i), c_set->at(j));
    }
  }
  std::swap(c_set->at(i+1),c_set->at(right));
  return i + 1;
}

// Check if entity exists, return true if does.
bool ent_conn::chk_ent(int tag) {
  for (int i = 0; i < conn.size(); i++) {
    if (tag == conn.at(i))
      return true;
  }
  return false;
}

int ent_conn::find_ent(int tag) {
  for (int i = 0; i < conn.size(); i++) {
    if (tag == conn.at(i))
      return i;
  }
  return -1;
}

// Add entity to the end of the removal list.
void ent_conn::add_ent(int tag) {
  if(!chk_ent(tag)) {
    conn.push_back(tag);
  }
}

// Remove entity with given tag from the removal list. 
// Shift later entities.
void ent_conn::rem_ent(int tag) {
  std::vector<int>::iterator it = std::find(conn.begin(), conn.end(), tag);
  if(it != conn.end())
    conn.erase(it);
}

// Entity indices higher than given tag are reduced by one.
void ent_conn::shift_ind(int tag) {
  for (int i = 0; i < conn.size(); i++) {
    if (conn.at(i) > tag) {
      conn.at(i) = conn.at(i) - 1;
    }
  }
}

// Try replacing tag2 with tag_new. Assumes at most single occurance. 
// Used for 1 cells. 
bool ent_conn::repl_tag(int tag_new, int tag_old) {
  for (int i = 0; i < conn.size(); i++) {
    if (conn.at(i) == tag_old) {
      conn.at(i) = tag_new;
      return true;
    }
  }
  return false;
}

// Quick Sort Functions for Ascending Order. TODO This is a buggy implementation!
// (2 Functions) 
//http://codereview.stackexchange.com/questions/77782/quick-sort-implementation


void ent_conn::quicksort(int left, int right) {
  if(right == -1 and left == -1) {
    left = 0;
    right = conn.size()-1;
  }

  // left > right
  if (left < right) {
    int part = partition(left, right);
    quicksort(left, part - 1);
    quicksort(part + 1, right);
  }
}

int ent_conn::partition(int left, int right) {
  int pivot = conn.at(right);
  // move the mid point value to the front.
  int i = left-1;
  int j = left;
  int temp;
  for (; j < right; j++) {
    if(conn.at(j) < pivot) {
      i++;
      temp = conn.at(i);
      conn.at(i) = conn.at(j);
      conn.at(j) = temp;
    }
  }
  temp = conn.at(right);
  conn.at(right) = conn.at(i+1);
  conn.at(i+1) = temp;
  return i + 1;
}

// Merge another connection list into this list. Assumes the second list 
// contains only unique elements. 
//int ent_conn::merge(struct ent_conn* e_conn2)
//{
//  if (conn.size() + ent_conn2->conn.size() < n) {
//    quicksort(0,conn.size());
//    e_conn2->quicksort(0,e_conn2->conn.size());
//    for (int i = 0; i < ent_conn2->conn.size(); i++) {
//     if chk_ent(ent_conn2->conn.at(i)) {
//        rem_ent(ent_conn2->conn.at(i);
//      }
//    }
//
//    for (int i = 0; i < ent_conn2->conn.size(); i++) {
//      if chk_ent(ent_conn2->conn.at(i)) {
//        rem_ent(ent_conn2->conn.at(i);
//      }
//    }
//    
//  }
//  else {
//    printf("Container size too small.\n");
//  }
//}

int ent_conn::chk_intsct(ent_conn* e2, ent_conn* intsct) {
  std::vector<int>::iterator it;
  std::sort (conn.begin(),conn.end());
  std::sort (e2->conn.begin(),e2->conn.end());

  intsct->conn.clear();
  intsct->conn.reserve(conn.size());

  std::set_intersection(conn.begin(), conn.end(), e2->conn.begin(), 
                             e2->conn.end(), 
                         std::inserter(intsct->conn, intsct->conn.begin()));
  return intsct->conn.size();
}

void ent_conn::unique() {
  quicksort();
  assert(conn.size() > 0);
  std::vector<int>::iterator it;

  it = std::unique (conn.begin(), conn.end()); 
  conn.resize(std::distance(conn.begin(),it));
}

void ent_conn::clear() {
  conn.clear();
}

void ent_conn::reserve(int n_new) {

  conn.reserve(n_new);
}

void ent_conn::resize(int n_new) {

  conn.resize(n_new);
}

void ent_conn::print() {
  for(int i = 0; i < conn.size(); i++) {
    std::cout << conn.at(i) + 1 << ", ";
  }
  std::cout << std::endl;
}

void ent_conn::append(ent_conn* e2) {
  conn.insert(conn.end(), e2->conn.begin(), e2->conn.end());
}


// ---------------------------------------------------------
// cell_adder can be used to add new cells.
// ---------------------------------------------------------
cell_adder::cell_adder() {
  e_new_id.resize(4);
  e_new_con.resize(3);
  e_new_add.resize(4);

  min_sz = 4;
}

void cell_adder::add_cell(int dim, int tag, std::vector<int>* new_con, 
                                            std::vector<int>* new_add) {
  //std::cout << "Adding " << dim << "c" << tag + 1 << std::endl;

  e_new_id.at(dim).push_back(tag);
  std::vector<int> conn(0);
  conn = *new_con;
  if(dim > 0)
    e_new_con.at(dim-1).push_back(conn);

  conn = *new_add;
  e_new_add.at(dim).push_back(conn);
}

void cell_adder::clear() {
  dummy_clear_stop();
  for(int i = 0; i < e_new_id.size(); i++) {
    e_new_id.at(i).clear();
  }
  e_new_id.clear();

  for(int i = 0; i < e_new_con.size(); i++) {
    for(int j = 0; j < e_new_con.at(i).size(); j++) {
      e_new_con.at(i).at(j).clear();
    }
    e_new_con.at(i).clear();
  }
  e_new_con.clear();

  for(int i = 0; i < e_new_add.size(); i++) {
    for(int j = 0; j < e_new_add.at(i).size(); j++) {
      e_new_add.at(i).at(j).clear();
    }
    e_new_add.at(i).clear();
  }
  e_new_add.clear();

  e_new_id.resize(4);
  e_new_con.resize(3);
  e_new_add.resize(4);
}

// ---------------------------------------------------------
// cell_base is the container for topological information.
// ---------------------------------------------------------

// ---------------------------------------------------------
// Regular constructor.
// File constructor.
// agm constructor.
// ---------------------------------------------------------

cell_base::cell_base(int n0, int n1, int n2, int n3) : 
  max_sz(100),
  cell_1_sz(0), cell_2_sz(0), cell_3_sz(0),
  con_1(0), con_2(0), con_3(0), 
  cell_free(4, std::vector<int>(0)),
  c_fmap(4, std::map<int, bool>{}), 
  ext_spur(false),
  ext_map(0, std::map<int, bool>{} ), cor_1c_map{}, cor_0c_map{}, 
  fix_list(0)
 {

  n[0] = n0;
  n[1] = n1;
  n[2] = n2;
  n[3] = n3;

  max_sz = 100;
  //cell_1_sz.clear();
  //cell_2_sz.clear();
  //cell_3_sz.clear();

  cell_1_sz.resize(n[1]);
  cell_2_sz.resize(n[2]);
  cell_3_sz.resize(n[3]);

  //con_1.clear();
  //con_2.clear();
  //con_3.clear();

  con_1.resize(2*n[1]);
  con_2.resize(max_sz*n[2]);
  con_3.resize(max_sz*n[3]);

  for(int dim = 0; dim < 4; dim++) {
    cell_free.at(dim).reserve(n[dim]);
    cell_free.at(dim).resize(0);
  }

  collect_ext();
  fix_free();

  ext_spur = false;
}

void cell_base::load_file(FILE* f) {
  load_flag = 1;
  file = f;

  int i,k;
  int tag;

  std::cout<<"Reading entity counts..."<<std::endl;
  // read entity counts
  gmi_fscanf(f, 4, "%d %d %d %d", &n[3], &n[2], &n[1], &n[0]);
  printf("N are... V:%d E:%d S:%d R:%d.\n", n[0], n[1], n[2], n[3]);

  // Set the memory:
  max_sz = 100;

  //cell_1_sz.clear();
  //cell_2_sz.clear();
  //cell_3_sz.clear();

  cell_1_sz.resize(n[1]);
  cell_2_sz.resize(n[2]);
  cell_3_sz.resize(n[3]);

  //con_1.clear();
  //con_2.clear();
  //con_3.clear();

  con_1.resize(2*n[1]);
  con_2.resize(max_sz*n[2]);
  con_3.resize(max_sz*n[3]);

  for(int dim = 0; dim < 4; dim++) {
    cell_free.at(dim).reserve(n[dim]);
    cell_free.at(dim).resize(0);
  }
  // End of memory setup.

  ent_conn e_con;

  std::cout<<"Entity counts read..."<<std::endl;
  // bounding box
  gmi_fscanf(f, 0, "%*f %*f %*f");
  gmi_fscanf(f, 0, "%*f %*f %*f");

  // 0cells //
  for (i = 0; i < n[0]; ++i) {
    int v1;
    gmi_fscanf(f, 1, "%d %*f %*f %*f", &v1);
    // std::cout<<"Vertex " << v1 <<std::endl;
  }

  std::cout<<"0-Cells read..."<<std::endl;
  // 1cells //
  for (i = 0; i < n[1]; ++i) {
    int tag1, tag2;
    gmi_fscanf(f, 3, "%d %d %d", &tag, &tag1, &tag2);
    printf("%d %d %d. ", tag-1, tag1, tag2);
    if (tag1 == -42) {
      if (tag2 == -42)
        e_con.clear();
      else {
        e_con.resize(1);
        e_con.conn.at(0) = tag2-1;
      }
    }
    else if (tag2 == -42) {
      e_con.resize(1);
      e_con.conn.at(0) = tag1-1;
    }

    set_conn(1, tag-1, &e_con);
    print_conn(1, tag-1);
  }
  std::cout<<"1-Cells read..."<<std::endl;

  int c_sz;
  // 2cells //
  for (i = 0; i < n[2]; ++i) {
    gmi_fscanf(f, 1, "%d %*d", &tag);
    gmi_fscanf(f, 1, "%d", &c_sz);
    e_con.conn.resize(c_sz);

    for (k = 0; k < e_con.conn.size(); ++k) {
      // tag, direction //
      gmi_fscanf(f, 1, "%d %*d", &e_con.conn.at(k));
      e_con.conn.at(k)--;
    }

    set_conn(2, tag-1, &e_con);
    print_conn(2, tag-1);
  }

  std::cout<<"2-Cells read..."<<std::endl;
  // 3cells //
  for (i = 0; i < n[3]; ++i) {
    gmi_fscanf(f, 1, "%d %*d", &tag);
    gmi_fscanf(f, 1, "%d", &c_sz);
    e_con.conn.resize(c_sz);

    for (k = 0; k < e_con.conn.size(); ++k) {
      // tag, direction //
      gmi_fscanf(f, 1, "%d %*d", &e_con.conn.at(k));
      e_con.conn.at(k)--;
    }

    if (e_con.conn.size() == 0) {
      cell_free.at(3).push_back(tag);
    }

    set_conn(3, tag-1, &e_con);
    print_conn(3, tag-1);
  }
  std::cout<<"3-Cells read..."<<std::endl;

  // Find the free 0cells:
  //struct ent_conn e_upper;
  //for (i = 0; i < n[0]; i++) {
  //  get_conn_dim(1, 0, i, &e_upper);
  //  if (e_upper.conn.size() == 0)
  //    cell_0_free.push_back(i);
  //}

  collect_ext();
  fix_free();

  ext_spur = false;
}

cell_base::cell_base(FILE* f) :
  max_sz(100),
  cell_1_sz(0), cell_2_sz(0), cell_3_sz(0),
  con_1(0), con_2(0), con_3(0), 
  cell_free(4, std::vector<int>(0)),
  c_fmap(4, std::map<int, bool>{}), 
  ext_spur(false),
  ext_map(0, std::map<int, bool>{} ), cor_1c_map{}, cor_0c_map{}, 
  fix_list(0)
 {
  load_file(f);
}

cell_base::cell_base(const char* modelFile) :
  max_sz(100),
  cell_1_sz(0), cell_2_sz(0), cell_3_sz(0),
  con_1(0), con_2(0), con_3(0), 
  cell_free(4, std::vector<int>(0)),
  c_fmap(4, std::map<int, bool>{}),  
  ext_spur(false),
  ext_map(0, std::map<int, bool>{} ), cor_1c_map{}, cor_0c_map{}, 
  fix_list(0)
 {

  FILE* f = fopen(modelFile, "r");
  load_file(f);
}

cell_base::cell_base(struct gmi_model* m) : 
  max_sz(100),
  cell_1_sz(0), cell_2_sz(0), cell_3_sz(0),
  con_1(0), con_2(0), con_3(0), 
  cell_free(4, std::vector<int>(0)),
  c_fmap(4, std::map<int, bool>{}), 
  ext_spur(false),
  ext_map(0, std::map<int, bool>{} ), cor_1c_map{}, cor_0c_map{}, 
  fix_list(0)
{

  load_flag = 0;
  mdl = m;

  struct gmi_iter* it;
  struct gmi_ent* e;
  struct gmi_set* s;

  int i;
  /* entity counts */
  n[0] = m->n[0];
  n[1] = m->n[1];
  n[2] = m->n[2];
  n[3] = m->n[3];

  // Set the memory:
  max_sz = 100;

  //cell_1_sz.clear();
  //cell_2_sz.clear();
  //cell_3_sz.clear();

  cell_1_sz.resize(n[1]);
  cell_2_sz.resize(n[2]);
  cell_3_sz.resize(n[3]);

  //con_1.clear();
  //con_2.clear();
  //con_3.clear();

  con_1.resize(2*n[1]);
  con_2.resize(max_sz*n[2]);
  con_3.resize(max_sz*n[3]);

  for(int dim = 0; dim < 4; dim++) {
    cell_free.at(dim).reserve(n[dim]);
    cell_free.at(dim).resize(0);
  }

  ent_conn e_con;
  int tag;
  // 0cells

  it = gmi_begin(m, 0);

  while ((e = gmi_next(m, it))) {
    s = gmi_adjacent(m, e, 1);
    tag = gmi_tag(m, e) - 1;
    if (s->n == 0) {
      cell_free.at(0).push_back(tag);
      c_fmap.at(0)[tag] = true;
    }
    gmi_free_set(s);
  }
  gmi_end(m, it);


  // 1cells
  it = gmi_begin(m, 1);

  while ((e = gmi_next(m, it))) {
    s = gmi_adjacent(m, e, 0);
    tag = gmi_tag(m, e) - 1;
    if (s->n == 0) {
      cell_free.at(1).push_back(tag);
      c_fmap.at(1)[tag] = true;
      e_con.resize(0);
    }
    else if (s->n == 1) {
      e_con.resize(1);
      e_con.conn.at(0) = gmi_tag(m, s->e[0]) - 1;
      //e_con.conn.at(1) = e_con.conn.at(0);
    }
    else {
      e_con.resize(2);
      e_con.conn.at(0) = gmi_tag(m, s->e[0]) - 1;
      e_con.conn.at(1) = gmi_tag(m, s->e[1]) - 1;
    }
    gmi_free_set(s);
    set_conn(1, tag, &e_con);
  }
  gmi_end(m, it);

  // 2cells
  it = gmi_begin(m, 2);
  while ((e = gmi_next(m, it))) {
    tag = gmi_tag(m, e) - 1;
    s = gmi_adjacent(m, e, 1);
    e_con.conn.resize(s->n);
    for (i = 0; i < e_con.conn.size(); ++i)
      e_con.conn.at(i) = gmi_tag(m, s->e[i]) - 1;

    if (s->n == 0) {
      cell_free.at(2).push_back(tag);
      c_fmap.at(2)[tag] = true;
    }

    gmi_free_set(s);
    set_conn(2, tag, &e_con);
  }
  gmi_end(m, it);

  // 3cells
  it = gmi_begin(m, 3);
  while ((e = gmi_next(m, it))) {
    tag = gmi_tag(m, e) - 1;
    s = gmi_adjacent(m, e, 2);
    e_con.conn.resize(s->n);
    for (i = 0; i < e_con.conn.size(); ++i)
      e_con.conn.at(i) = gmi_tag(m, s->e[i]) - 1;

    if (s->n == 0) {
      cell_free.at(3).push_back(tag);
      c_fmap.at(3)[tag] = true;
    }

    gmi_free_set(s);
    set_conn(3, tag, &e_con);
  }
  gmi_end(m, it);

  collect_ext();
  fix_free();

  ext_spur = false;
}
/*
// This is to collect the topological entities used in the mesh.
cell_base::cell_base(apf::Mesh2* m) : 
  max_sz(100),
  cell_1_sz(0), cell_2_sz(0), cell_3_sz(0),
  con_1(0), con_2(0), con_3(0), 
  cell_0_free(0), cell_1_free(0), cell_2_free(0), cell_3_free(0),
  ext_spur(false),
  ext_map(0, std::map<int, bool>{} ), cor_1c_map{}, cor_0c_map{}, 
  fix_list(0) {

  load_flag = 0;

  apf::MeshIterator* it;
  apf::MeshEntity* ent;

  for(int dim = 0; dim < 4; dim++) {
    int count = 0;
    it = m->begin(dim);

    std::map<apf::ModelEntity*, bool> mdl_map{};
    while (ent = m->iterate(it)) {
      apf::ModelEntity* mdl_ent = m->toModel(ent);
      if(m->getModelType(mdl_ent) == dim) {
        int id = m->getModelTag(mdl_ent);
        assert(!mdl_map[id]);
        mdl_map[id] = true;
        if(count < id - 1)
          count = id - 1;
      }
    }
    m->end(it);

    n[dim] = count;
  }

  // Downward adjacency lists.
  std::vector<std::list<int> > > conn_1c(n[1], std::list<int> (0));
  std::vector<std::list<int> > > conn_2c(n[2], std::list<int> (0));
  std::vector<std::list<int> > > conn_3c(n[3], std::list<int> (0));

  cell_1_sz.resize(n[1]);
  cell_2_sz.resize(n[2]);
  cell_3_sz.resize(n[3]);

  con_1.resize(2*n[1]);
  con_2.resize(max_sz*n[2]);
  con_3.resize(max_sz*n[3]);

  cell_0_free.reserve(n[0]);
  cell_1_free.reserve(n[1]);
  cell_2_free.reserve(n[2]);
  cell_3_free.reserve(n[3]);

  ent_conn e_con;
  int tag;
  // 0cells
  int dim = 1;
  it = m->begin(dim);
  std::vector<std::map<int, bool> > conn_map(n[1], std::map<int, bool> {});
  apf::Downward d_v;
  while (ent = m->iterate(it)) {
    apf::ModelEntity* mdl_ent = m->toModel(ent);
    if(m->getModelType(mdl_ent) == dim) {
      int id = m->getModelTag(mdl_ent);
      m->getDownward(e, 0, d_v);

      for(int i = 0; i < dim + 1; i++) {
        mdl_ent = m->toModel(d_v[i]);
        if(m->getModelType(mdl_ent) == dim - 1) {
          int id_d = m->getModelTag(mdl_ent);
          if(!conn_map.at(id - 1)[id_d]) {
            conn_map.at(id - 1)[id_d] = true;
          }
        }
      }

      if(count < id - 1)
        count = id - 1;
    }
  }
  m->end(it);


  collect_ext();
  fix_free();

  ext_spur = false;
}
*/
// Reload the topology.
void cell_base::reload_topo() {
  printf("Load flag is:%d\n", load_flag);

  if (load_flag == 0) 
    reload_topo(mdl);
  else if (load_flag == 1)
    reload_topo(file);
  else
    printf("No reload option set.\n");
}

void cell_base::reload_topo(struct gmi_model* m) {
  clear();
  struct gmi_iter* it;
  struct gmi_ent* e;
  struct gmi_set* s;

  int i;
  //entity counts
  n[0] = m->n[0];
  n[1] = m->n[1];
  n[2] = m->n[2];
  n[3] = m->n[3];

  // Set the memory:
  max_sz = 100;

  dummy_clear_stop();

  cell_1_sz.clear();
  cell_2_sz.clear();
  cell_3_sz.clear();

  cell_1_sz.resize(n[1]);
  cell_2_sz.resize(n[2]);
  cell_3_sz.resize(n[3]);

  con_1.clear();
  con_2.clear();
  con_3.clear();

  con_1.resize(2*n[1]);
  con_2.resize(max_sz*n[2]);
  con_3.resize(max_sz*n[3]);

  for(int dim = 0; dim < 4; dim++) {
    cell_free.at(dim).reserve(n[dim]);
    cell_free.at(dim).resize(0);
  }

  ent_conn e_con;

  int tag;
  // 0cells

  it = gmi_begin(m, 0);

  while ((e = gmi_next(m, it))) {
    s = gmi_adjacent(m, e, 1);
    tag = gmi_tag(m, e) - 1;
    if (s->n == 0) {
      cell_free.at(0).push_back(tag);
      c_fmap.at(0)[tag] = true;
    }
    gmi_free_set(s);
  }
  gmi_end(m, it);

  // 1cells
  it = gmi_begin(m, 1);

  while ((e = gmi_next(m, it))) {
    s = gmi_adjacent(m, e, 0);
    tag = gmi_tag(m, e) - 1;
    if (s->n == 0) {
      e_con.conn.clear();
      cell_free.at(1).push_back(tag);
      c_fmap.at(1)[tag] = true;
    }
    else if (s->n == 1) {
      e_con.resize(2);
      e_con.conn.at(0) = gmi_tag(m, s->e[0]) - 1;
      e_con.conn.at(1) = e_con.conn.at(0);
    }
    else {
      e_con.resize(2);
      e_con.conn.at(0) = gmi_tag(m, s->e[0]) - 1;
      e_con.conn.at(1) = gmi_tag(m, s->e[1]) - 1;
    }
    gmi_free_set(s);
    set_conn(1, tag, &e_con);
  }
  gmi_end(m, it);

  // 2cells
  it = gmi_begin(m, 2);
  while ((e = gmi_next(m, it))) {
    tag = gmi_tag(m, e) - 1;
    s = gmi_adjacent(m, e, 1);
    e_con.resize(s->n);
    for (i = 0; i < e_con.conn.size(); ++i)
      e_con.conn.at(i) = gmi_tag(m, s->e[i]) - 1;

    if (s->n == 0) {
      cell_free.at(2).push_back(tag);
      c_fmap.at(2)[tag] = true;
    }

    gmi_free_set(s);
    set_conn(2, tag, &e_con);
  }
  gmi_end(m, it);

  // 3cells
  it = gmi_begin(m, 3);
  while ((e = gmi_next(m, it))) {
    tag = gmi_tag(m, e) - 1;
    s = gmi_adjacent(m, e, 2);
    e_con.resize(s->n);
    for (i = 0; i < e_con.conn.size(); ++i)
      e_con.conn.at(i) = gmi_tag(m, s->e[i]) - 1;

    if (s->n == 0) {
      cell_free.at(3).push_back(tag);
      c_fmap.at(3)[tag] = true;
    }

    gmi_free_set(s);
    set_conn(3, tag, &e_con);
  }
  gmi_end(m, it);

  std::cout << "Cell base reloaded." << std::endl;

  collect_ext();
  fix_free();

}

void cell_base::reload_topo(FILE* f) {
  load_file(f);
}

// ---------------------------------------------------------
// Get sizes for connection and cell entity lists.
// ---------------------------------------------------------

int cell_base::get_sz(int dim) {
  if (dim < 4 and dim >= 0) {
    return n[dim];
  }
  else {
    return 0;
  }
}

int cell_base::cell_sz(int dim, int tag) {

  if (dim == 1) {
    return cell_1_sz.at(tag);
  }
  else if (dim == 2) {
    return cell_2_sz.at(tag);
  }
  else if (dim == 3) {
    return cell_3_sz.at(tag);
  }
  else {
    printf("Cell_sz dim %d is not appropriate.\n", dim);
  }
}

int cell_base::get_conn_sz(int dim, int tag) {
  if (dim < 4 and dim > 0) {
    if (n[dim] > tag)
      return cell_sz(dim,tag);
    else
      return 0;
  }
  else
    return 0;
}



// ---------------------------------------------------------
// Connection list operations, get and set using ent_conn.
// ---------------------------------------------------------
// LOOK merge all dims
// Set the connection list of a given entity. 
bool cell_base::set_conn(int dim, int tag, struct ent_conn* e_con) {

  if (dim > 3 or dim < 1 or tag > n[dim])
    return false;

  //std::cout << "Setting " << dim << "c" << tag + 1 << " ";
  //for (int i = 0; i < e_con->conn.size(); i++)
  //  std::cout << dim - 1 << "c" << e_con->conn.at(i) + 1 << " ";
  //std::cout << std::endl;

  //else if (dim == 1) {
  if (dim == 1) {
    cell_1_sz.at(tag) = e_con->conn.size();
    for (int i = 0; i < e_con->conn.size(); i++) {
      con_1.at(2*tag + i) = e_con->conn.at(i);
    }
    return true;
  }

  else if (dim == 2) {

    assert(e_con->conn.size() < max_sz);
    cell_2_sz.at(tag) = e_con->conn.size();
    for (int i = 0; i < e_con->conn.size(); i++) {
      con_2.at(max_sz*tag + i) = e_con->conn.at(i);
    }
    return true;
  }

  else if (dim == 3) {

    assert(e_con->conn.size() < max_sz);
    cell_3_sz.at(tag) = e_con->conn.size();
    for (int i = 0; i < e_con->conn.size(); i++) {
      con_3.at(max_sz*tag + i) = e_con->conn.at(i);
    }
    return true;
  }

  else {
    std::cout<<"Dimension of entity not correct."<<std::endl;
    return false;
  }

}

// Get the connection list of a given entity. 
bool cell_base::get_conn(int dim, int tag, struct ent_conn* e_con) {

  if (dim > 3 or dim < 1 or tag > n[dim]) {
    e_con->resize(0);
    return false;
  }

  else if (dim == 1) {
    //std::cout << "Tag " << tag << " n " << n[dim]
    //          << " sz " << cell_1_sz.at(tag) 
    //          << std::endl;
    e_con->resize(cell_1_sz.at(tag));
    for (int i = 0; i < cell_1_sz.at(tag); i++)
      e_con->conn.at(i) = con_1.at(2*tag + i);
    return true;
  }

  else if (dim == 2) {
    // printf("%d.\n",conn.at(0));
    e_con->resize(cell_2_sz.at(tag));
    for (int i = 0; i < e_con->conn.size(); i++) {
      e_con->conn.at(i) = con_2.at(max_sz*tag + i);
    }
    return true;
  }

  else if (dim == 3) {
    // printf("%d.\n",conn.at(0));
    e_con->resize(cell_3_sz.at(tag));
    for (int i = 0; i < e_con->conn.size(); i++) {
      e_con->conn.at(i) = con_3.at(max_sz*tag + i);
    }
    return true;
  }

  else {
    std::cout<<"Dimension of entity not correct."<<std::endl;
    return false;
  }

}

// GMI versions of the get_conn and set_conn. Indices are shifted by 1.

bool cell_base::set_conn_gmi(int dim, int tag, struct ent_conn* e_con) {
  for (int i = 0; i < e_con->conn.size(); i++)
    e_con->conn.at(i) = e_con->conn.at(i) - 1;
  set_conn(dim, tag-1, e_con);
  for (int i = 0; i < e_con->conn.size(); i++)
    e_con->conn.at(i) = e_con->conn.at(i) + 1;
}

bool cell_base::get_conn_gmi(int dim, int tag, struct ent_conn* e_con) {

  get_conn(dim, tag-1, e_con);

  for (int i = 0; i < e_con->conn.size(); i++)
    e_con->conn.at(i) = e_con->conn.at(i) + 1;

}

// Get the given lower dimensional connection list of a given entity. 
bool cell_base::get_conn_lower(int dim_adj, int dim, int tag, 
                                                struct ent_conn* e_cover) {

  struct ent_conn e_con;
  struct ent_conn e_next;

  e_cover->conn.clear();

  if(dim_adj == dim-1) {
    get_conn(dim, tag, e_cover);
  }
  else {
    get_conn(dim, tag, &e_con);

    for (int i = 0; i < e_con.conn.size(); i++) {
      get_conn_lower(dim_adj, dim-1, e_con.conn.at(i), &e_next); 
      for (int j = 0; j < e_next.conn.size(); j++) {
        if(!e_cover->chk_ent(e_next.conn.at(j))) {
          e_cover->add_ent(e_next.conn.at(j));
        }
      }
    }
  }
  return true;
}

// Get the given lower dimensional connection list of a given entity. 
bool cell_base::get_conn_lower(int dim_adj, int dim, struct ent_conn* e_set,
                                                struct ent_conn* e_cover) {

  struct ent_conn e_con;
  struct ent_conn e_next;

  if(dim_adj == dim-1) {
    int n_set = 0;
    for (int i = 0; i < e_set->conn.size(); i++) {
      n_set = n_set + cell_sz(dim, e_set->conn.at(i));
    }
    e_cover->resize(n_set);
    for (int i = 0; i < e_set->conn.size(); i++) {
      get_conn(dim, e_set->conn.at(i), &e_con);
      e_cover->append(&e_con);
    }
    e_cover->unique();
  }
  else {
    e_cover->resize(get_sz(dim_adj));

    for (int i = 0; i < e_set->conn.size(); i++) {
      get_conn(dim, e_set->conn.at(i), &e_con);
      for (int j = 0; j < e_con.conn.size(); j++) {
        get_conn_lower(dim_adj, dim-1, e_con.conn.at(i), &e_next); 
        e_cover->append(&e_next);
      }
    }
    e_cover->unique();
  }
  return true;
}

// TODO
// Get the given dimensional connection list of a given entity. 
bool cell_base::get_conn_dim(int dim_adj, int dim, int tag, 
                                                struct ent_conn* e_cover) {

  if (dim_adj > 3 or dim_adj < 0 or dim > 3 or dim < 0)
    return false;
  // Check the dim order. If dim > dim2, simply merge lower dim until dim2 is 
  // reached.

  // Else, put tag in a list. Collect higher dim containing list elements, 
  // until dim2 is reached. 

  // It feels this could be simpler in the dim < dim2 case, as the lower dims 
  // to upper dim connections are not as numerous as the upper dim to lower dim
  // connections. Still, this might be a limitation of the hierarchical cell
  // complex.

  e_cover->conn.clear();

  struct ent_conn e_next;
  struct ent_conn e_con;
  if (dim_adj == dim) {
    // Get the cell connections. Get the cell connections of same dim. cells.
    // Compare to find common members. If found, add to e_cover. 
    if (dim == 3) {
      // Get the 2cell adjacencies. TODO fix the indices
      get_conn(dim, tag, &e_con);
      for (int i = 0; i < get_sz(3); i++) {
        get_conn(dim, i, &e_next);
        for(int cell_bound = 0; cell_bound < e_con.conn.size(); cell_bound++) {
          if (e_next.chk_ent(e_con.conn.at(cell_bound))) {
            cell_bound = e_con.conn.size();
            e_cover->add_ent(i);
          }
        }
      }
      // Remove the tag from the list.
      e_cover->rem_ent(tag);
    }
    // If dim < 3, go one lever higher, chk_ent(tag). If found, add all cells
    // to the cover list, not found on the cover list. Remove tag.
    else {
      for (int i = 0; i < get_sz(dim+1); i++) {
        get_conn(dim+1, i, &e_con);
        if (e_con.chk_ent(tag)) {
          for (int j = 0; j < e_con.conn.size(); j++) {
            if (e_con.conn.at(j) != tag and 
                !e_cover->chk_ent(e_con.conn.at(j))) {
              e_cover->add_ent(e_con.conn.at(j));
            }
          }
        }
      }
    }
  }

  else if (dim_adj < dim) {
    get_conn_lower(dim_adj, dim, tag, e_cover);
  }
  // dim_adj > dim.
  else {
    struct ent_conn e_con;
    for (int i = 0; i < get_sz(dim_adj); i++) {
      // Get the dim dimensional adjacencies of ith dim_adj dimensional cell.
      // If it contains the actual cell, add the ith cell to the cover.
      get_conn_lower(dim, dim_adj, i, &e_con);
      if(e_con.chk_ent(tag)) {
        e_cover->add_ent(i);
      }
    }
  }

  return true;
}

bool cell_base::get_conn_dim_gmi(int dim_adj, int dim, int tag, 
                                                struct ent_conn* e_cover) {

  bool flag = get_conn_dim(dim_adj, dim, tag-1, e_cover);

  for (int i = 0; i < e_cover->conn.size(); i++)
    e_cover->conn.at(i)++;
  return flag;
}

// Get the list of 1cell neighbors of a 0cell, connected to a specific 3cell. 
// This can be used during 0/3cell detachment, to obtain the 1cells connected 
// to the new 0cell and those connected to the old 0cell.

bool cell_base::get_conn_13(int cell0, int cell3, 
                                                struct ent_conn* e_adj,
                                                struct ent_conn* e_not) {

  e_adj->clear();
  e_not->clear();
  // Create two empty 1cell list, e_adj, e_not.
  // Get the 2cell conn list of the 3cell.
  int sz = get_sz(1);
  struct ent_conn e_1cell;
  for (int i = 0; i < sz; i++) {
    get_conn(1, i, &e_1cell);
    if(e_1cell.chk_ent(cell0))
      e_not->add_ent(i);
  }

  // Going over the 1cells, check 0cell adjacency. Add them to e_not.
  // Going over the 2cell list, check for each 1cell adjacency of the current 
  // 2cell in the e_not list. If exists, move from e_not to e_adj.
  struct ent_conn e_3cell;
  get_conn(3, cell3, &e_3cell);

  struct ent_conn e_2cell;
  for (int i = 0; i < e_3cell.conn.size(); i++) {
    get_conn(2, e_3cell.conn.at(i), &e_2cell);
    for (int j = 0; j < e_2cell.conn.size(); j++) {
      get_conn(1, e_2cell.conn.at(j), &e_1cell);
      for (int k = 0; k < e_not->conn.size(); k++) {
        if(e_1cell.chk_ent(e_not->conn.at(k))) {
          e_adj->add_ent(e_not->conn.at(k));
          e_not->rem_ent(e_not->conn.at(k));
        }
      }
    }
  }

}

bool cell_base::get_conn_12(int cell0, int cell2, struct ent_conn* e_adj) {

  e_adj->clear(); 

  struct ent_conn e_2cell;
  struct ent_conn e_1cell;
  get_conn(2, cell2, &e_2cell);

  for (int i = 0; i < e_2cell.conn.size(); i++) {
    get_conn(1, e_2cell.conn.at(i), &e_1cell);
    if(e_1cell.chk_ent(cell0)) {
      e_adj->add_ent(e_2cell.conn.at(i));
    }
  }

  assert(e_adj->conn.size() == 0 or e_adj->conn.size() == 2);
}

// Given a d-cell-cc check if the cell is bounded by 0cell0.
bool cell_base::get_conn_0(int cell0, int c_dim, int c_id) {
  assert(c_dim > 0 and c_dim < 4);
  struct ent_conn e_cc;

  if(c_dim == 1) {
    get_conn(1, c_id, &e_cc);
    if(e_cc.chk_ent(cell0)) {
      return true;
    }
    return false;
  }
  else {
    //struct ent_conn e_sub;
    get_conn(c_dim, c_id, &e_cc);
    for (int j = 0; j < e_cc.conn.size(); j++) {
      if(get_conn_0(cell0, c_dim-1, e_cc.conn.at(j)) )
        return true;
    }
    return false;
  }
}

// Given two cells, check if one is bounded by the other.
bool cell_base::chk_conn_d(int c_dim1, int c_id1, int c_dim2, int c_id2) {
  assert(c_dim1 > -1 and c_dim1 < 4);
  assert(c_dim2 > -1 and c_dim2 < 4);
  if(c_dim1 == c_dim2)
    return c_id1 == c_id2;

  if(c_dim1 > c_dim2) {
    int temp = c_dim1;
    c_dim1 = c_dim2;
    c_dim2 = temp;

    temp = c_id1;
    c_id1 = c_id2;
    c_id2 = temp;
  }

  struct ent_conn e_cc;

  if(c_dim2 - c_dim1 == 1) {
    get_conn(c_dim2, c_id2, &e_cc);
    return e_cc.chk_ent(c_id1);
  }

  else {
    struct ent_conn e_sub;
    get_conn(c_dim2, c_id2, &e_cc);
    for (int j = 0; j < e_cc.conn.size(); j++) {
      if(chk_conn_d(c_dim1, c_id1, c_dim2-1, e_cc.conn.at(j)) )
        return true;
    }
    return false;
  }
}

// Given two cells, check if one is bounded by the other.
bool cell_base::chk_conn_d_gmi(int c_dim1, int c_id1, int c_dim2, int c_id2) {
  return chk_conn_d(c_dim1, c_id1-1, c_dim2, c_id2-1);
}

// Given two cells of the same dim, check if they have a joint boundary.
bool cell_base::chk_conn_joint(int c_dim, int c_id1, int c_id2, int dim_other) {
  assert(c_dim > -1 and c_dim < 4);
  if(dim_other == -1)
    dim_other = c_dim - 1;
  assert(dim_other > -1);

  struct ent_conn e1;
  struct ent_conn e2;

  get_conn_dim(dim_other, c_dim, c_id1, &e1);
  get_conn_dim(dim_other, c_dim, c_id2, &e2);

  for (int i = 0; i < e1.conn.size(); i++) {
    if(e2.chk_ent(e1.conn.at(i))) {
      return true;
    }
  }
  return false;
}

// Given two cells of the same dim, check if they have a joint boundary.
bool cell_base::chk_conn_joint_gmi(int c_dim, int c_id1, int c_id2, 
                                                      int dim_other) {
  return chk_conn_joint_gmi(c_dim, c_id1 - 1, c_id2 - 1, dim_other);
}

// Given a d-cell-c mark the bounding cells  as adj_c and cells bounded by 
// bounding 0cells not belonging to the bounding set as end_c.
void cell_base::get_end_gmi(int c_dim, int c_id, 
                        std::vector<std::map<int, bool> >& end_c, 
                        std::vector<std::map<int, bool> >& adj_c) {
  assert(adj_c.size() == 3 and end_c.size() == 3);
  for(int i = 0; i < 3; i++) {
    adj_c.at(i).clear();
    end_c.at(i).clear();
  }

  std::map<int, bool> burn_3c{};

  ent_conn* e_0c = new ent_conn();
  ent_conn* e_3c = new ent_conn();
  ent_conn* e_down = new ent_conn();

  for(int dim_low = c_dim - 1; dim_low > - 1; dim_low--) {
    get_conn_dim_gmi(dim_low, c_dim, c_id, e_0c);
    for(int i = 0; i < e_0c->conn.size(); i++) {
      adj_c.at(dim_low)[e_0c->conn.at(i)] = true;
    }
  }

  for(int i = 0; i < e_0c->conn.size(); i++) {
    get_conn_dim_gmi(3, 0, e_0c->conn.at(i), e_3c);
    for(int j = 0; j < e_3c->conn.size(); j++) {
      if(!burn_3c[e_3c->conn.at(j)]) {
        burn_3c[e_3c->conn.at(j)] = true;
        for(int dim_low = 0; dim_low < 3; dim_low++) {
          get_conn_dim_gmi(dim_low, 3, e_3c->conn.at(j), e_down);
          for(int k = 0; k < e_down->conn.size(); k++) {
            if(!adj_c.at(dim_low)[e_down->conn.at(k)]
              or !end_c.at(dim_low)[e_down->conn.at(k)])
              end_c.at(dim_low)[e_down->conn.at(k)] = true;
          }
        }
      }
    }
  }
  delete e_0c;
  delete e_3c;
  delete e_down;
}


// Given two cells of the same dim, check if they have a joint boundary.
ent_conn cell_base::get_conn_joint(int c_dim, int c_id1, int c_id2, int dim_other) {
  if(dim_other == -1)
    dim_other = c_dim - 1;
  assert(dim_other > -1);

  struct ent_conn e1;
  struct ent_conn e2;
  struct ent_conn ej;

  get_conn_dim(dim_other, c_dim, c_id1, &e1);
  get_conn_dim(dim_other, c_dim, c_id2, &e2);

  if(e1.conn.size() < e1.conn.size()) {
    for (int i = 0; i < e1.conn.size(); i++) {
      if(e2.chk_ent(e1.conn.at(i)))
        ej.add_ent(e1.conn.at(i));
    }
  }
  else {
    for (int i = 0; i < e1.conn.size(); i++) {
      if(e2.chk_ent(e1.conn.at(i)))
        ej.add_ent(e1.conn.at(i));
    }
  }
  ej.quicksort();
  return ej;
}

// Given two cells of the same dim, check if they have a joint boundary.
ent_conn cell_base::get_conn_joint_gmi(int c_dim, int c_id1, int c_id2, int dim_other) {
  struct ent_conn ej;
  ej = get_conn_joint(c_dim, c_id1-1, c_id2-1);
  for (int i = 0; i < ej.conn.size(); i++) {
    ej.conn.at(i) = ej.conn.at(i) + 1;
  }
  return ej;
}

// Given a 1cell, check if it is on the boundary corner.
bool cell_base::chk_cell_ext_gmi(int c_dim, int c_id) {
  if(c_dim == 0)
    return chk_0cell_ext_gmi(c_id);
  else if(c_dim == 1)
    return chk_1cell_ext_gmi(c_id);
  else if(c_dim == 2)
    return chk_2cell_ext_gmi(c_id);
}

bool cell_base::get_cell_ext(int c_dim, int c_id) {
  return ext_map.at(c_dim)[c_id];
}

bool cell_base::get_cell_ext_gmi(int c_dim, int c_id) {
  return ext_map.at(c_dim)[c_id-1];
}

bool cell_base::get_cell_bound_ext(int c_dim, int c_id) {
  if(c_dim < 3 and ext_map.at(c_dim)[c_id])
    return true;
  else {
    ent_conn e_0c;
    get_conn_dim(0, c_dim, c_id, &e_0c);
    for(int i = 0; i < e_0c.conn.size(); i++) {
      if(ext_map.at(0)[e_0c.conn.at(i)])
        return true;
    }
  }
  return false;
}

bool cell_base::get_cell_bound_ext_gmi(int c_dim, int c_id) {
  return get_cell_bound_ext(c_dim, c_id - 1);
}

bool cell_base::get_1c_corner_gmi(int c_id) {
  return cor_1c_map[c_id - 1];
}

bool cell_base::get_1c_corner(int c_id) {
  return cor_1c_map[c_id];
}

bool cell_base::get_0c_corner_gmi(int c_id) {
  return cor_0c_map[c_id - 1];
}

bool cell_base::get_0c_corner(int c_id) {
  return cor_0c_map[c_id];
}

// Given a 0cell, check if it is on the boundary.
bool cell_base::chk_0cell_ext_gmi(int c_id) {
  struct ent_conn e_2cells;
  struct ent_conn e_3cells;
  get_conn_dim_gmi(2, 0, c_id, &e_2cells);
  for (int i = 0; i < e_2cells.conn.size(); i++) {
    get_conn_dim_gmi(3, 2, e_2cells.conn.at(i), &e_3cells);
    if(e_3cells.conn.size() != 2)
      return true;
  }
  return false;
}

// Given a 0cell, check if it is on the boundary.
bool cell_base::chk_1cell_ext_gmi(int c_id) {
  return chk_1cell_ext_c(c_id-1);
}

// Given a 1cell, check if it is on the boundary corner.
bool cell_base::chk_1cell_ext_c(int c_id) {
  struct ent_conn e_2cells;
  struct ent_conn e_3cells;
  get_conn_dim(2, 1, c_id, &e_2cells);
  for (int i = 0; i < e_2cells.conn.size(); i++) {
    get_conn_dim_gmi(3, 2, e_2cells.conn.at(i), &e_3cells);
    if(e_3cells.conn.size() != 2)
      return true;
  }
  return false;
}

// Given a 1cell, check if it is on the boundary corner.
bool cell_base::chk_2cell_ext_gmi(int c_id) {
  struct ent_conn e_3cells;
  get_conn_dim_gmi(3, 2, c_id, &e_3cells);
  if(e_3cells.conn.size() != 2)
    return true;
  return false;
}

// Given a 1cell, check if it is on the boundary corner.
bool cell_base::chk_2cell_ext(int c_id) {
  struct ent_conn e_3cells;
  get_conn_dim(3, 2, c_id, &e_3cells);
  if(e_3cells.conn.size() == 1)
    return true;
  return false;
}

// Get the list of d-cells neighbors of a 0cell. d-cell > 0 and d-cell < 4.
void cell_base::get_conn_0d(int cell0, int c_dim, struct ent_conn* e_adj) {
  assert(c_dim > 0 and c_dim < 4);

  if(c_dim == 1) {
    struct ent_conn e_cc;
    for (int i = 0; i < get_sz(1); i++) {
      get_conn(1, i, &e_cc);
      if(e_cc.chk_ent(cell0)) {
        e_adj->add_ent(i);
      }
    }
  }

  else {

    for (int i = 0; i < get_sz(c_dim); i++) {
      if(get_conn_0(cell0, c_dim, i) )
          e_adj->add_ent(i);
    }
  }
}

// Given a d_list-cell list, remove the cells in the list that do not bound 
// or are not bounded by c_dim-cell-c_id.
void cell_base::rem_conn(int c_dim, int c_id, int d_list, 
                                    struct ent_conn* e_adj) {

  //std::cout << "Checking " << c_dim << "cell" << c_id + 1
  //          << "vs " << d_list<< "cells" << std::endl;

  int shift = 0;
  int i = 0;
  while(i < e_adj->conn.size()-shift) {
    if (chk_conn_d(c_dim, c_id, d_list, e_adj->conn.at(i+shift))) {
      e_adj->conn.at(i) = e_adj->conn.at(i+shift);
      i++;
    }
    else {
      //std::cout << "Skipping " << e_adj->conn.at(i+shift)+1 << std::endl;
      shift++;
    }
  }
  e_adj->resize(e_adj->conn.size() - shift);
}

void cell_base::rem_conn_gmi(int c_dim, int c_id, int d_list, 
                                    struct ent_conn* e_adj) {
  for(int i = 0; i < e_adj->conn.size(); i++)
    e_adj->conn.at(i) = e_adj->conn.at(i) - 1;

  rem_conn(c_dim, c_id-1, d_list, e_adj);

  for(int i = 0; i < e_adj->conn.size(); i++)
    e_adj->conn.at(i) = e_adj->conn.at(i) + 1;
}


// Get the list of d-cells neighbors of a 0cell. d-cell > 0 and d-cell < 4.
void cell_base::get_conn_0d_gmi(int cell0, int c_dim, struct ent_conn* e_adj) {
  get_conn_0d(cell0-1, c_dim, e_adj);
  for(int i = 0; i < e_adj->conn.size(); i++)
    e_adj->conn.at(i) = e_adj->conn.at(i) + 1;
}

/*
// Given a 1cell and a list of 2cells, find the 2cells accessible by 
// starting from the 1cell and going over 1 and 2cells, touching the 0cell.
// Remove the disconnected 2cells.
void cell_base::get_conn_12d_gmi(int cell0, int cell1,
                                                 struct ent_conn* e_2cells) {
  for(int i = 0; i < e_2cells->conn.size(); i++) {
    if(chk_conn_d_gmi(1, cell1, 2, e_2cells->conn.at(i));
  }
get_conn_0(int cell0, int c_dim, int c_id)
}
*/

// Given a 2,3-cell, find the set of collapsing 1cells in order to reduce
// all all bounding 2-cells to be bounded by at most three 1-cells. 
// Fills them into the ent_conn list. 
void cell_base::get_col_1cell(int cell_dim, int cell_id, int tag_0cell, 
                                            struct ent_conn* e_conn) {

  ent_conn cell_cand; // Max 100 1cells
  e_conn->clear();
  e_conn->reserve(get_sz(1));

  if(cell_dim == 3) {
    // Get the bounding 2 cells.
    ent_conn cell_2; // Max 100 2cells
    ent_conn cell_1; // Max 100 1cells
    get_conn(cell_dim, cell_id, &cell_2);
    // Add the bounding 1 cells of these 2 cells to the e_conn list.
    for (int i = 0; i < cell_2.conn.size(); i++) {
      get_conn(2, cell_2.conn.at(i), &cell_1);
      // Add the available 1cells to the list.
      for (int j = 0; j < cell_1.conn.size(); j++) {
        if (!cell_cand.chk_ent(cell_1.conn.at(j)))
          cell_cand.add_ent(cell_1.conn.at(j));
      }
    }
  }
  else if(cell_dim == 2) {
    ent_conn cell_1; // Max 100 1cells
    // Add the bounding 1 cells of the 2 cell to the e_conn list.
    get_conn(cell_dim, cell_id, &cell_cand);
  }
  else if(cell_dim == 1) {
    cell_cand.add_ent(cell_id);
  }
  else {
    printf("Get_col_1cell: Dimension not appropriate.\n");
  }

  printf("Candidates are:\n");
  for (int i = 0; i < cell_cand.conn.size(); i++) {
    printf("1cell%d, ", cell_cand.conn.at(i)+1);
  }
  printf("\n");

  // Check whether the candidate 1cell is collapsable. If so collapse.
  for (int i = 0; i < cell_cand.conn.size(); i++) {
    ent_conn cell_1; // Max 100 2cells
    bool collapse = true;
    int tag = cell_cand.conn.at(i);
    int sz_2 = get_sz(2);
    printf("1cell%d\n",tag+1);
    for (int j = 0; j < sz_2; j++) {
      get_conn(2, j, &cell_1);
      // If 1cell is contained within this 2cell which has 3 or less 1cells 
      // bounding, this 1cell cannot be collapsed. Not all 4 1cell 2cells are 
      // catched.
      if ((cell_1.chk_ent(tag)) and cell_1.conn.size() < 4) {
        print_conn(2, j);
        j = sz_2;
        collapse = false;
      }
    }
    if (collapse) {
      e_conn->add_ent(tag);
      //coll_cell(1, tag, tag_0cell);
    }
  }

  //reload_topo();

}

// Calls the get_col_1cell by using the dmg notation.
void cell_base::get_col_1cell_gmi(int cell_dim, int cell_id, int tag_0cell, 
                                            struct ent_conn* e_conn) {
  get_col_1cell(cell_dim, cell_id-1, tag_0cell-1, e_conn);

  for(int i = 0; i < e_conn->conn.size(); i++) {
    e_conn->conn.at(i) = e_conn->conn.at(i) + 1;
  }

}

// ---------------------------------------------------------
// Removal and merge operations, eject and collapse (merges two lower dim).
// These update the necessary cell lists.
// ---------------------------------------------------------
/*
// Eject the entity, naive removal. Updates the same dim cell list.
bool cell_base::ej_cell(int dim, int tag) {
  // printf("Ejecting %dcell %d:\n", dim, tag);
  ent_conn e_con;
  int sz = get_sz(dim);

  if ((not sz > 0) or tag > n[dim]) {
    printf("Not valid.\n");
    return false;
  }
  else {
    // printf("Shifting %dcell entities.\n", dim);
    for (int i = tag; i < sz-1; i++) {
      // printf("Entity %d, ", i);
      get_conn(dim, i+1, &e_con);
      set_conn(dim, i, &e_con);
    }
    // printf("Shifted %dcell entities. ", dim);
    n[dim]--;
    // printf("Number of entites reduced.\n");
    return true;
  }
}
*/
// Remove the entity. Go over the upper cell connection list, check if tag 
// exists. If so, remove tag. Change the indices of the entities with 
// higher index.
bool cell_base::rem_cell(int dim, int tag) {

  if (free_cell(dim, tag)) {
    if (get_sz(dim+1) > 0) {
      ent_conn e_con;
      ent_conn e_up;
      get_conn_dim(dim+1, dim, tag, &e_up);

      for (int i = 0; i < e_up.conn.size(); i++) {
        get_conn(dim+1, e_up.conn.at(i), &e_con);
        while (e_con.chk_ent(tag)) {
          e_con.rem_ent(tag);
        }
        // In any case, shift indices. Not anymore.
        // e_con.shift_ind(tag);
        set_conn(dim+1, e_up.conn.at(i), &e_con);
      }
    }
    struct ent_conn e_empty; 
    set_conn(dim, tag, &e_empty);

    if(dim < 3)
      ext_map.at(dim)[tag] = false;
    if(dim == 1)
      cor_1c_map[tag] = false;
    if(dim == 0)
      cor_0c_map[tag] = false;
    return true;
  }
  return false;
}
/*
bool cell_base::rem_cell(int dim, int tag) {

  if (free_cell(dim, tag)) {
    if (get_sz(dim+1) > 0) {
      ent_conn e_con;
      for (int i = 0; i < n[dim+1]; i++) {
        get_conn(dim+1, i, &e_con);
        while (e_con.chk_ent(tag)) {
          e_con.rem_ent(tag);
        }
        // In any case, shift indices. Not anymore.
        // e_con.shift_ind(tag);
        set_conn(dim+1, i, &e_con);
      }
    }
    struct ent_conn e_empty; 
    set_conn(dim, tag, &e_empty);

    if(dim < 3)
      ext_map.at(dim)[tag] = false;
    if(dim == 1)
      cor_1c_map[tag] = false;
    if(dim == 0)
      cor_0c_map[tag] = false;
    return true;
  }
  return false;
}
*/
// Replace cell1 with cell2. Go over the upper cell connection list, check 
// if tag exists. If so, check if the replacement also exists. If not, 
// replace tag. 
bool cell_base::repl_cell(int dim, int tag_new, int tag_old) {
  if (get_sz(dim+1) > 0) {
    ent_conn e_con;
    for (int i = 0; i < n[dim+1]; i++) {
      get_conn(dim+1, i, &e_con);
      if (e_con.chk_ent(tag_old)) { // tag1 does exist.
        if (e_con.chk_ent(tag_new)) { // tag2 does exist.
          e_con.rem_ent(tag_old);
        }
        else {
          e_con.repl_tag(tag_new, tag_old);
        }
      }
      // In any case, shift indices. Not anymore.
      // e_con.shift_ind(tag);
      set_conn(dim+1, i, &e_con);
    }
  }
}

// Extend the cell list of given dim. Add the new cells to the free list. 
// This should be called when a given free cell list is empty.
bool cell_base::add_free(int dim, int add_sz) {
  //std::cout << "Adding free " << dim << " " << add_sz << std::endl;
/*
  if(get_free(dim) > -1) {
    printf("add_free called for non-empty %ddim free list.\n", dim);
    return false;
  }
  else {
*/
    std::vector<int>* cell_sz;
    int cont_sz = 0;
    if (dim == 0) {
    }
    else if (dim == 1) {
      cont_sz = 2;
      cell_sz = &cell_1_sz;
      con_1.resize(con_1.size()+add_sz*cont_sz);
      cell_sz->reserve(n[dim] + add_sz);
    }
    else if (dim == 2) {
      cont_sz = max_sz;
      cell_sz = &cell_2_sz;
      con_2.resize(con_2.size()+add_sz*max_sz);
      cell_sz->reserve(n[dim] + add_sz);
    }
    else if (dim == 3) {
      cont_sz = max_sz;
      cell_sz = &cell_3_sz;
      con_3.resize(con_3.size()+add_sz*max_sz);
      cell_sz->reserve(n[dim] + add_sz);
    }

    cell_free.at(dim).reserve(n[dim] + add_sz);

    // Adding the new cells to the free cell list. The new ones are added such
    // that the lower index ones are added to the end.
    for (int i = 1; i < add_sz + 1; i++) {
      cell_free.at(dim).push_back(n[dim] + add_sz - i);
      c_fmap.at(dim)[n[dim] + add_sz - i] = true;
    }
    // Adding the new cell also requires update of the connection list size.
    // The unused ones should have 0.
    if(dim != 0) {
      for (int i = 0; i < add_sz; i++) {
        cell_sz->push_back(0);
      }
    }

    if(add_sz > 0)
      n[dim] = n[dim] + add_sz;
    //std::cout << "Added free " << dim << " " << get_free_sz(dim) << std::endl;
    return true;
  //}
}

// Return whether the cell is in free list.
bool cell_base::is_free(int dim, int cell_id) {
/*
  if (dim == 0)
    return (std::find(cell_0_free.begin(),
                      cell_0_free.end(), cell_id) != cell_0_free.end());
  else if (dim == 1)
    return (std::find(cell_1_free.begin(),
                      cell_1_free.end(), cell_id) != cell_1_free.end());
  else if (dim == 2)
    return (std::find(cell_2_free.begin(),
                      cell_2_free.end(), cell_id) != cell_2_free.end());
  else if (dim == 3)
    return (std::find(cell_3_free.begin(),
                      cell_3_free.end(), cell_id) != cell_3_free.end());
*/
  return c_fmap.at(dim)[cell_id];
}

// Return the index of the free cell of given dim. If empty, return -1.
int cell_base::get_free(int dim) {
  int free_size = get_free_sz(dim);
  if (free_size > 0) {
    return cell_free.at(dim).back();
  }
  else
    return -1;
}

int cell_base::get_free_sz(int dim) {
  int free_size = 0.;
  if (dim > -1 and dim < 4)
    free_size = cell_free.at(dim).size();

  return free_size;
}

int cell_base::use_free_gmi(int dim) {
  return use_free(dim)+1;
}

// Return the index of the free cell of given dim. If empty, add new cells,
// and return the index of the new free cell.
int cell_base::use_free(int dim) {
  int fr_cell = get_free(dim);
  if (fr_cell > -1) {
    cell_free.at(dim).pop_back();
    c_fmap.at(dim)[fr_cell] = false;
    //std::cout << "Using free " << dim << "c" << fr_cell+1 << std::endl;
    //std::cout << "Left " << get_free_sz(dim) << std::endl;
    return fr_cell;
  }
  else {
    //add_free(dim, 10);
    //return use_free(dim);
    return -1;
  }
}

// TODO Obsolete
/*
void cell_base::remove_free(int dim, int tag) {
  std::vector<int>::iterator i_it;
  if (dim == 0) {
    i_it = std::find(cell_0_free.begin(), cell_0_free.end(), tag);
    assert(i_it != cell_0_free.end());
    cell_0_free.erase(i_it);
  }
  else if (dim == 1) {
    i_it = std::find(cell_1_free.begin(), cell_1_free.end(), tag);
    assert(i_it != cell_1_free.end());
    cell_1_free.erase(i_it);
  }
  else if (dim == 2) {
    i_it = std::find(cell_2_free.begin(), cell_2_free.end(), tag);
    assert(i_it != cell_2_free.end());
    cell_2_free.erase(i_it);
  }
  else if (dim == 3) {
    i_it = std::find(cell_3_free.begin(), cell_3_free.end(), tag);
    assert(i_it != cell_3_free.end());
    cell_3_free.erase(i_it);
  }
}
*/

// After removal from cell list, add the freed cell to the end of the free
// list.
bool cell_base::free_cell(int dim, int tag) {
  ent_conn e_con;
  int sz = get_sz(dim);

  if ((not sz > 0) or (tag < 0 or tag > n[dim])) {
    printf("Not valid, sz:%d, tag:%d, dim:%d, n:%d.\n", sz, tag, dim, n[dim]);
    return false;
  }
  if (is_free(dim,tag)) {
    return false;
  }

  if (dim > -1 and dim < 4) {
    cell_free.at(dim).push_back(tag);
    c_fmap.at(dim)[tag] = true;
  }
  return true;
}


bool cell_base::chk_free(int dim, int cell_id) {
  struct ent_conn e_cover;
  bool up = true;
  bool down = true;

  if(dim < 3) {
    get_conn_dim(dim+1, dim, cell_id, &e_cover);
    up = (e_cover.conn.size() > 0);
    if (up)
      return false;
    else
      up = true;
  }
  else {
    up = true;
  }

  if(dim > 0) {
    get_conn_dim(dim-1, dim, cell_id, &e_cover);
    down = (e_cover.conn.size() > 0);
    if (down)
      return false;
    else
      down = true;
  }
  else {
    down = true;
  }

  return (down and up);
}

void cell_base::fix_free() {

  for(int dim = 0; dim < 4; dim++) {
    cell_free.at(dim).clear();
    c_fmap.at(dim).clear();
  }

  struct ent_conn e_down;

  for (int dim = 3; dim > -1; dim--) {
    int sz = get_sz(dim);
    for (int i = 0; i < sz; i++) {
      c_fmap.at(dim)[i] = true;
    }
  }

  int dim = 3;
  int sz = get_sz(dim);
  for (int i = 0; i < sz; i++) {
    get_conn(dim, i, &e_down);
    if(e_down.conn.size() == 0) {
      cell_free.at(dim).push_back(i);
      c_fmap.at(dim)[i] = true;
    }
    else {
      c_fmap.at(dim)[i] = false;
      for (int j = 0; j < e_down.conn.size(); j++) {
        c_fmap.at(dim-1)[e_down.conn.at(j)] = false;
      }
    }
  }

  for (int dim = 2; dim > 0; dim--) {
    sz = get_sz(dim);
    for (int i = 0; i < sz; i++) {
      // If a stratum has no upper adjacency, it is assumed to be free. 
      // Else, mark the one lower dimensional adjacencies as not free.
      if(!c_fmap.at(dim)[i]) {
        get_conn(dim, i, &e_down);
        for (int j = 0; j < e_down.conn.size(); j++) {
          c_fmap.at(dim-1)[e_down.conn.at(j)] = false;
        }
      }
      else {
        cell_free.at(dim).push_back(i);
        c_fmap.at(dim)[i] = true;
      }
    }
  }

  dim = 0;
  sz = get_sz(dim);
  for (int i = 0; i < sz; i++) {
    if(c_fmap.at(dim)[i]) {
      cell_free.at(dim).push_back(i);
    }
  }
}



// Collapse a given cell. Connected 0_cell are merged. All other lower 
// dimensional cells are removed. The remaining 0cell is the input 0cell.
void cell_base::coll_cell(int dim, int tag, int tag_0cell) {
  if (dim == 1) {
    //printf("Collapsing 1cell %d.\n", tag+1);
    ent_conn e_con;
    ent_conn e_con_loop;
    ent_conn e_con_1c;
    if (get_conn(dim, tag, &e_con)) {
      bool ext_curr = false;
      for(int i = 0; i < e_con.conn.size(); i++) {
        if(get_cell_ext(0, e_con.conn.at(i)))
          ext_curr = true;
      }
      if (e_con.conn.size() > 0) {
        if (tag_0cell < 0) {
          tag_0cell = e_con.conn.at(0);
        }
        if(tag_0cell != e_con.conn.at(0)) {
          get_conn_dim(1, 0, e_con.conn.at(0), &e_con_1c);
          for (int i = 0; i < e_con_1c.conn.size(); i++) {
            if (tag != e_con_1c.conn.at(i)) {
              get_conn(dim, e_con_1c.conn.at(i), &e_con_loop);
              if(e_con_loop.chk_ent(tag_0cell))
                e_con_loop.rem_ent(e_con.conn.at(0));
              else
                e_con_loop.repl_tag(tag_0cell, e_con.conn.at(0));
              set_conn(dim, e_con_1c.conn.at(i), &e_con_loop);
            }
          }
          rem_cell(0, e_con.conn.at(0));
        }
        if (e_con.conn.size() == 2 and (e_con.conn.at(1) != tag_0cell)) {
          get_conn_dim(1, 0, e_con.conn.at(1), &e_con_1c);
          for (int i = 0; i < e_con_1c.conn.size(); i++) {
            if (tag != e_con_1c.conn.at(i)) {
              get_conn(dim, e_con_1c.conn.at(i), &e_con_loop);
              if(e_con_loop.chk_ent(tag_0cell))
                e_con_loop.rem_ent(e_con.conn.at(1));
              else
                e_con_loop.repl_tag(tag_0cell, e_con.conn.at(1));
              set_conn(dim, e_con_1c.conn.at(i), &e_con_loop);
            }
          }
          rem_cell(0, e_con.conn.at(1));
        }
      }
      set_ext(0, tag_0cell, ext_curr);
    }
    rem_cell(dim, tag);
  }

  // When collapsing higher dim cells, first collapse lower dim cells.
  else if (dim == 2 or dim == 3) {
    //printf("Collapsing %dcell %d.\n", dim, tag+1);
    ent_conn e_con;
    ent_conn e_con_loop;
    if (get_conn(dim, tag, &e_con)) {
      for (int i = 0; i < e_con.conn.size(); i++) {
        coll_cell(dim-1, e_con.conn.at(i), tag_0cell);
        // e_con.shift_ind(e_con.conn.at(i));
      }
      rem_cell(dim, tag);
    }
  }
}
/*
// Collapse a given cell. Connected 0_cell are merged. All other lower 
// dimensional cells are removed. The remaining 0cell is the input 0cell.
void cell_base::coll_cell(int dim, int tag, int tag_0cell) {
  if (dim == 1) {
    //printf("Collapsing 1cell %d.\n", tag+1);
    ent_conn e_con;
    ent_conn e_con_loop;
    if (get_conn(dim, tag, &e_con)) {
      bool ext_curr = false;
      for(int i = 0; i < e_con.conn.size(); i++) {
        if(get_cell_ext(0, e_con.conn.at(i)))
          ext_curr = true;
      }
      if (e_con.conn.size() > 0) {
        if (tag_0cell < 0) {
          tag_0cell = e_con.conn.at(0);
        }
        else if(tag_0cell != e_con.conn.at(0)) {
          for (int i = 0; i < n[1]; i++) {
            if (get_conn(dim, i, &e_con_loop)) {
              if(e_con_loop.chk_ent(tag_0cell))
                e_con_loop.rem_ent(e_con.conn.at(0));
              else
                e_con_loop.repl_tag(tag_0cell, e_con.conn.at(0));
              set_conn(dim, i, &e_con_loop);
            }
          }
          rem_cell(0, e_con.conn.at(0));
        }
        if (e_con.conn.size() == 2 and e_con.conn.at(1) != tag_0cell) {
          for (int i = 0; i < n[1]; i++) {
            if (get_conn(dim, i, &e_con_loop)) {
              if(e_con_loop.chk_ent(tag_0cell))
                e_con_loop.rem_ent(e_con.conn.at(1));
              else
                e_con_loop.repl_tag(tag_0cell, e_con.conn.at(1));
              set_conn(dim, i, &e_con_loop);
            }
          }
          rem_cell(0, e_con.conn.at(1));
        }
      }
      set_ext(0, tag_0cell, ext_curr);
    }
    rem_cell(dim, tag);
  }

  // When collapsing higher dim cells, first collapse lower dim cells.
  else if (dim == 2 or dim == 3) {
    //printf("Collapsing %dcell %d.\n", dim, tag+1);
    ent_conn e_con;
    ent_conn e_con_loop;
    if (get_conn(dim, tag, &e_con)) {
      for (int i = 0; i < e_con.conn.size(); i++) {
        coll_cell(dim-1, e_con.conn.at(i), tag_0cell);
        // e_con.shift_ind(e_con.conn.at(i));
      }
      rem_cell(dim, tag);
    }
  }
}
*/
void cell_base::coll_cell_gmi(int dim, int tag, int tag_0cell) {
  coll_cell(dim, tag-1, tag_0cell-1);
}

// After collapse, empty cells should be checked and removed.

void cell_base::clean_cells() {
  clean_dim(1);
  clean_dim(2);
  clean_dim(3);
}

void cell_base::clean_dim(int dim) {
  int i = 0;

  while (i < n[dim]) {
    int sz_curr = get_conn_sz(dim,i);
    if (sz_curr == 0) {
      printf("%dcell %d is empty. ", dim, i+1);
/*
      if (dim == 1)
        printf("Size:%d.\n", cell_1_sz.at(i));
      else if (dim == 2)
        printf("Size:%d.\n", cell_2_sz.at(i));
      else if (dim == 3)
        printf("Size:%d.\n", cell_3_sz.at(i));
*/
      printf("Size:%d.\n", sz_curr);

      rem_cell(dim,i);
    }
    else {
      i++;
    }
  }
}

// ---------------------------------------------------------
// Insert operations.
// ---------------------------------------------------------
void cell_base::ins_cell(cell_adder* ca) {
  assert(ca->dim > 0 and ca->dim < 3);
  struct ent_conn e_con;

  std::cout << "Inserting ";
  for(int dim = 1; dim <= ca->dim; dim++) {
    std::cout << dim << "c " << std::endl;
    for(int i = 0; i < ca->e_new_con.at(dim-1).size(); i++) {
      std::cout << "\t" << ca->e_new_id.at(dim).at(i) + 1 << std::endl;
      e_con = ca->e_new_con.at(dim-1).at(i);
      set_conn(dim, ca->e_new_id.at(dim).at(i), &e_con);
      std::cout << "New conn:" << std::endl;
      e_con.print();

      std::cout << "Adding to " << dim + 1 << "cells" << std::endl;
      for(int j = 0; j < ca->e_new_add.at(dim).at(i).size(); j++) {
        std::cout << "\t\t" << ca->e_new_add.at(dim).at(i).at(j) + 1 
                  << std::endl;
        get_conn(dim+1, ca->e_new_add.at(dim).at(i).at(j), &e_con);
        e_con.add_ent(ca->e_new_id.at(dim).at(i));
        set_conn(dim+1, ca->e_new_add.at(dim).at(i).at(j), &e_con);
        std::cout << "adj:" << std::endl;
        e_con.print();
        if(dim < 2) {
          if(ext_map.at(dim+1)[ca->e_new_add.at(dim).at(i).at(j)])
            ext_map.at(dim)[ca->e_new_id.at(dim).at(i)] = true;
        }
      }
    }
  }

  std::cout << "0c: " << ca->tag_0 + 1 << std::endl;
  for(int i = 0; i < ca->e_new_id.at(0).size(); i++) {
    std::cout << "0c" << ca->e_new_id.at(0).at(i) + 1 << std::endl;
    for(int j = 0; j < ca->e_new_add.at(0).at(i).size(); j++) {
      get_conn(1, ca->e_new_add.at(0).at(i).at(j), &e_con);
      e_con.repl_tag(ca->e_new_id.at(0).at(i), ca->tag_0);
      set_conn(1, ca->e_new_add.at(0).at(i).at(j), &e_con);
      e_con.print();

      if(ext_map.at(1)[ca->e_new_add.at(0).at(i).at(j)])
        ext_map.at(0)[ca->e_new_id.at(0).at(i)] = true;
    }
  }
}

// ---------------------------------------------------------
// Determination of the exterior boundaries. OBSOLETE as this needs mesh 
// information. TODO
// ---------------------------------------------------------
void cell_base::collect_ext() {

  for(int i = 0; i < ext_map.size(); i++) {
    ext_map.at(i).clear();
  }
  ext_map.clear();
  ext_map.resize(3);

  struct ent_conn e1;
  struct ent_conn e0;

  for(int i = 0; i < n[2]; i++) {
    if(chk_2cell_ext(i)) {
      get_conn(2, i, &e1);
      for(int j = 0; j < e1.conn.size(); j++) {
        get_conn(1, e1.conn.at(j), &e0);
        for(int k = 0; k < e0.conn.size(); k++) {
          ext_map.at(0)[e0.conn.at(k)] = true;
        }
        ext_map.at(1)[e1.conn.at(j)] = true;
      }
      ext_map.at(2)[i] = true;
    }
  }

  std::map<int, int> c3_adj{};
  for(int i = 0; i < n[3]; i++) {
    get_conn_dim(1, 3, i, &e1);
    for(int j = 0; j < e1.conn.size(); j++) {
      c3_adj[e1.conn.at(j)] = c3_adj[e1.conn.at(j)] + 1;
    }
  }
  for(int i = 0; i < n[1]; i++) {
    //std::cout << "1c" << i + 1
    //          << " c3_adj " << c3_adj[i] << std::endl;
    if(c3_adj[i] == 1) {
      cor_1c_map[i] = true;
      get_conn(1, i, &e0);
      for(int j = 0; j < e0.conn.size(); j++) {
        //std::cout << " 0c" << e0.conn.at(j) + 1;
        cor_0c_map[e0.conn.at(j)] = true;
      }
    }
    //std::cout << std::endl;
  }

}

void cell_base::set_ext_spur(bool state) {
  ext_spur = state;
}

void cell_base::set_ext(int c_dim, int c_tag, bool ext) {
  if(!is_free(c_dim, c_tag)) {
    ext_map.at(c_dim)[c_tag] = ext;
  }
}

void cell_base::set_ext_gmi(int c_dim, int c_tag, bool ext) {
  if(!is_free(c_dim, c_tag-1)) {
    ext_map.at(c_dim)[c_tag-1] = ext;
  }
}

// 1cell exterior corner flag.
void cell_base::set_cor_ext(int c_tag, bool ext) {
  if(!is_free(1, c_tag)) {
    cor_1c_map[c_tag] = ext;
    if(ext) {
      struct ent_conn e0;
      get_conn(1, c_tag, &e0);
      for(int i = 0; i < e0.conn.size(); i++) {
        cor_0c_map[e0.conn.at(i)] = ext;
      }
    }
  }
}

void cell_base::set_cor_ext_0c(int c_tag, bool ext) {
  if(!is_free(0, c_tag))
    cor_0c_map[c_tag] = ext;
}

void cell_base::set_cor_ext_0c_gmi(int c_tag, bool ext) {
  set_cor_ext_0c(c_tag - 1, ext);
}

void cell_base::set_cor_ext_gmi(int c_tag, bool ext) {
  set_cor_ext(c_tag-1, ext);
}

void cell_base::clear_cor() {
  cor_0c_map.clear();
  cor_1c_map.clear();
}

void cell_base::set_cor_ext(std::vector<int>* c1_in) {
  for(int i = 0; i < c1_in->size(); i++) {
    set_cor_ext(c1_in->at(i), true);
  }
}

void cell_base::set_cor_ext_gmi(std::vector<int>* c1_in) {
  for(int i = 0; i < c1_in->size(); i++) {
    set_cor_ext(c1_in->at(i)-1, true);
  }
}

// ---------------------------------------------------------
// Detection of spurious cells.
// ---------------------------------------------------------
// Merge two given cells. The second cell is replaced and removed. 
void cell_base::merg_cell(int dim, int tag1, int tag2) {
  struct ent_conn e1;
  struct ent_conn e2;

  get_conn(dim, tag1, &e1);
  get_conn(dim, tag2, &e2);

  e1.append(&e2);
  e1.unique();

  set_conn(dim, tag1, &e1);
  rem_cell(dim, tag2);
}

// Check the bounded cells, if none exists return true.
bool cell_base::chk_spur_up(int dim, int tag) {
  struct ent_conn e_up;
  get_conn_dim(dim+1, dim, tag, &e_up);
  if(e_up.conn.size() == 0)
    return true;
  return false;
}

bool cell_base::chk_spur(int dim, int tag) {
  struct ent_conn e_up;
  struct ent_conn e_down;
  assert(dim < 3);

  get_conn_dim(dim+1, dim, tag, &e_up);
  get_conn(dim, tag, &e_down);

  // 0cell: 3 > 1cell adj
  // 1cell: 3 > 2cell adj
  // 2cell: Both 3cell adj are same.
  if(dim == 2) {
    if(e_up.conn.size() == 2)
      return (e_up.conn.at(0) == e_up.conn.at(1));
  }
  else {
    //if(!chk_1cell_ext_c(tag) and e_up.conn.size() < 3)
    if(e_up.conn.size() < 3) {
      if(dim == 0)
        return true;
      return !get_1c_corner(tag);
    }
  }
  return false;
}

bool cell_base::fix_spur(int dim, int tag, bool clean) {
  if(is_free(dim,tag))
    return false;

  struct ent_conn e_up;
  struct ent_conn e_adj;
  std::vector<ent_conn> e_down(dim, ent_conn());
  assert(dim < 3);

  get_conn_dim(dim+1, dim, tag, &e_up);
  for(int dim_low = 0; dim_low < dim; ++dim_low)
    get_conn_dim(dim_low, dim, tag, &e_down.at(dim_low));

  if(clean) {
    clear_fix();
  }

  // 0cell: 3 > 1cell adj
  // 1cell: 3 > 2cell adj
  // 2cell: Both 3cell adj are same.
  bool fixed = false;

  if(!get_cell_ext(dim, tag)) {
    if(dim == 2) {
      if(e_up.conn.size() == 2 and e_up.conn.at(0) == e_up.conn.at(1)) {
        //get_conn(dim+1, e_up.conn.at(0), &e_adj);
        //e_adj.rem_ent(tag);
        //set_conn(dim+1, e_up.conn.at(0), &e_adj);
        c_elem celem_curr = std::make_pair(
                              std::make_pair(dim, tag + 1),
                              std::make_pair(dim+1, e_up.conn.at(0) + 1) );
        fix_list.push_back(celem_curr);

        rem_cell(dim, tag);
        fixed = true;
      }
    }
    else {
      if(e_up.conn.size() < 3) {
        assert(e_up.conn.size() == 2);
        c_elem celem_curr = std::make_pair(
                              std::make_pair(dim, tag + 1),
                              std::make_pair(dim+1, e_up.conn.at(0) + 1) );
        fix_list.push_back(celem_curr);

        rem_cell(dim, tag);

        //if(e_up.conn.size() == 2) {
          merg_cell(dim+1, e_up.conn.at(0), e_up.conn.at(1));
          celem_curr = std::make_pair(std::make_pair(dim+1, e_up.conn.at(1)+1),
                                  std::make_pair(dim+1, e_up.conn.at(0) + 1) );
          fix_list.push_back(celem_curr);
        //}

        fixed = true;
      }
    }
  }
  else {
    if(ext_spur) {
      if(dim < 2 and e_up.conn.size() == 2) {
        if((dim == 1 and !get_1c_corner(tag)) 
            or (dim == 0)) {
          //  or (dim == 0 and !get_0c_corner(tag))) {
          c_elem celem_curr = std::make_pair(
                                std::make_pair(dim, tag + 1),
                                std::make_pair(dim+1, e_up.conn.at(0) + 1) );
          fix_list.push_back(celem_curr);

          rem_cell(dim, tag);

          merg_cell(dim+1, e_up.conn.at(0), e_up.conn.at(1));
          celem_curr = std::make_pair(std::make_pair(dim+1, e_up.conn.at(1)+1),
                                  std::make_pair(dim+1, e_up.conn.at(0) + 1) );
          fix_list.push_back(celem_curr);

          fixed = true;
        }
      }
    }
  }
  if(dim > 0) {
    for(int dim_low = dim - 1; dim_low > -1; dim_low--) {
      for(int i = 0; i < e_down.at(dim_low).conn.size(); i++) {
        bool fix_curr = fix_spur(dim_low, e_down.at(dim_low).conn.at(i), false);
        fixed = (fixed or fix_curr);
      }
    }
  }
  return fixed;
}

// The fixes to spurious cells occur recursively, so the fix_list may contain
// replacing cells that are replaced by others. Going over the fix list, replace
// the intermediate replacements with the final ones.
void cell_base::process_fix_list() {
  std::map<std::pair<int, int>, bool> replaced{};
  std::map<std::pair<int, int>, int> rp_id{};

  for(int i = 0; i < fix_list.size(); i++) {
    replaced[fix_list.at(i).first] = true;
    rp_id[fix_list.at(i).first] = i;
  }

  for(int i = 0; i < fix_list.size(); i++) {
    if(replaced[fix_list.at(i).second]) {
      int i_rpl = rp_id[fix_list.at(i).second];
      fix_list.at(i).second = fix_list.at(i_rpl).second;
      i = i - 1;
    }
  }
}

/*

bool cell_base::fix_spur2(int dim, int tag, bool clean) {
  struct ent_conn e_up;
  struct ent_conn e_down;
  struct ent_conn e_adj;
  assert(dim < 3);

  get_conn_dim(dim+1, dim, tag, &e_up);
  get_conn(dim, tag, &e_down);

  if(clean) {
    clear_fix();
  }

  // 0cell: 3 > 1cell adj
  // 1cell: 3 > 2cell adj
  // 2cell: Both 3cell adj are same.

  if(!get_cell_ext(dim, tag)) {
    if(dim == 2) {
      if(e_up.conn.size() == 2 and e_up.conn.at(0) == e_up.conn.at(1)) {
        //get_conn(dim+1, e_up.conn.at(0), &e_adj);
        //e_adj.rem_ent(tag);
        //set_conn(dim+1, e_up.conn.at(0), &e_adj);
        c_elem celem_curr = std::make_pair(
                              std::make_pair(dim, tag + 1),
                              std::make_pair(dim+1, e_up.conn.at(0) + 1) );
        fix_list.push_back(celem_curr);

        rem_cell(dim, tag);

        for(int i = 0; i < e_down.conn.size(); i++)
          fix_spur(dim-1, e_down.conn.at(i), false);
        return true;
      }
    }
    else {
      if(e_up.conn.size() < 3) {
        c_elem celem_curr = std::make_pair(
                              std::make_pair(dim, tag + 1),
                              std::make_pair(dim+1, e_up.conn.at(0) + 1) );
        fix_list.push_back(celem_curr);

        rem_cell(dim, tag);

        if(e_up.conn.size() == 2) {
          merg_cell(dim+1, e_up.conn.at(0), e_up.conn.at(1));
          celem_curr = std::make_pair(std::make_pair(dim+1, e_up.conn.at(1)+1),
                                  std::make_pair(dim+1, e_up.conn.at(0) + 1) );
          fix_list.push_back(celem_curr);
        }

        for(int i = 0; i < e_down.conn.size(); i++)
          fix_spur(dim-1, e_down.conn.at(i), false);
        return true;
      }
    }
  }
  else {
    if(ext_spur) {
      if(dim < 2 and e_up.conn.size() == 2) {
        if((dim == 1 and !get_1c_corner(tag)) 
            or (dim == 0)) {
          //  or (dim == 0 and !get_0c_corner(tag))) {
          c_elem celem_curr = std::make_pair(
                                std::make_pair(dim, tag + 1),
                                std::make_pair(dim+1, e_up.conn.at(0) + 1) );
          fix_list.push_back(celem_curr);

          rem_cell(dim, tag);

          merg_cell(dim+1, e_up.conn.at(0), e_up.conn.at(1));
          celem_curr = std::make_pair(std::make_pair(dim+1, e_up.conn.at(1)+1),
                                  std::make_pair(dim+1, e_up.conn.at(0) + 1) );
          fix_list.push_back(celem_curr);

          for(int i = 0; i < e_down.conn.size(); i++)
            fix_spur(dim-1, e_down.conn.at(i), false);
          return true;
        }
      }
    }
  }
  return false;
}
*/
bool cell_base::chk_2c(int tag) {
  struct ent_conn e_adj;
  struct ent_conn e_1c;
  get_conn_0d(tag, 1, &e_adj);

  for(int i = 0; i < e_adj.conn.size(); i++) {
    //std::cout << "1c" << e_adj.conn.at(i) + 1 << std::endl;
    get_conn(1, e_adj.conn.at(i), &e_1c);
    //std::cout << "\t0c" << e_1c.conn.at(0) + 1 << " "
    //          << e_1c.conn.at(1) + 1 << std::endl;
    //if(e_1c.conn.at(0) == e_1c.conn.at(1) ) {
    if(get_conn_sz(1, e_adj.conn.at(i)) == 1) {
      assert(e_1c.conn.at(0) == tag);
      return true;
    }
  }
  return false;
}

bool cell_base::chk_2c_gmi(int tag) {
  return chk_2c(tag-1);
}

// ---------------------------------------------------------
// Cell generation operations.
// ---------------------------------------------------------
// Given a collection of 2cells, add a new 1cell to all lists, by splitting
// the 0cell.
bool cell_base::add_1cell(int tag_0cell, struct ent_conn* two_cell) {
  // First, check if all 2cells actually contain the 0cell:
  for (int i = 0; i < two_cell->conn.size(); i++) {
    struct ent_conn e_cover;
    get_conn_lower(0, 2, two_cell->conn.at(i), &e_cover);
    if(!e_cover.chk_ent(tag_0cell) ) {
      return false;
    }
  }
  int tag_1cell;
  int tag_0cell_new;
  if(get_free(1) > -1)
    tag_1cell = use_free(1);
  else {
    add_free(1, 10);
    tag_1cell = use_free(1);
  }
  if(get_free(0) > -1)
    tag_0cell_new = use_free(0);
  else {
    add_free(0, 10);
    tag_0cell_new = use_free(0);
  }

  struct ent_conn e_con;
  get_conn(1, tag_1cell, &e_con);

  e_con.conn.clear();
  e_con.conn.resize(2);
  e_con.conn.at(0) = tag_0cell;
  e_con.conn.at(1) = tag_0cell_new;
  set_conn(1, tag_1cell, &e_con);

  for (int i = 0; i < two_cell->conn.size(); i++) {
    get_conn(2, two_cell->conn.at(i), &e_con);
    e_con.add_ent(tag_1cell);
    set_conn(2, two_cell->conn.at(i), &e_con);
  }
}


// ---------------------------------------------------------
// Modification list operations. Mainly used in lens collapse.
// ---------------------------------------------------------
/*
// Set a topological modification. 
void cell_base::add_mod(int dim, int tag, int flag) {
  if (dim == 0) {
    cell_0_mod.push_back(tag);
    cell_0_flag.push_back(flag);
  }
  if (dim == 1) {
    cell_1_mod.push_back(tag);
    cell_1_flag.push_back(flag);
  }
  if (dim == 2) {
    cell_2_mod.push_back(tag);
    cell_2_flag.push_back(flag);
  }
  if (dim == 3) {
    cell_3_mod.push_back(tag);
    cell_3_flag.push_back(flag);
  }
}

// Apply list of modifications, starting from above.
void cell_base::apply_mod() {
  while(cell_3_mod.size() > 0) {
    if(cell_3_flag.back() == -1) {
      rem_cell(3, cell_3_mod.back());
    }
    else {
      repl_cell(3, cell_3_mod.back(), cell_3_flag.back());
    }
    cell_3_free.push_back(cell_3_mod.back());
    cell_3_mod.pop_back();
    cell_3_flag.pop_back();
  }

  while(cell_2_mod.size() > 0) {
    if(cell_2_flag.back() == -1) {
      rem_cell(2, cell_2_mod.back());
    }
    else {
      repl_cell(2, cell_2_mod.back(), cell_2_flag.back());
    }
    cell_2_free.push_back(cell_2_mod.back());
    cell_2_mod.pop_back();
    cell_2_flag.pop_back();
  }

  while(cell_1_mod.size() > 0) {
    if(cell_1_flag.back() == -1) {
      rem_cell(1, cell_1_mod.back());
    }
    else {
      repl_cell(1, cell_1_mod.back(), cell_1_flag.back());
    }
    cell_1_free.push_back(cell_1_mod.back());
    cell_1_mod.pop_back();
    cell_1_flag.pop_back();
  }

  while(cell_0_mod.size() > 0) {
    if(cell_0_flag.back() == -1) {
      rem_cell(0, cell_0_mod.back());
    }
    else {
      repl_cell(0, cell_0_mod.back(), cell_0_flag.back());
    }
    cell_0_free.push_back(cell_0_mod.back());
    cell_0_mod.pop_back();
    cell_0_flag.pop_back();
  }

}

// Clean list of modifications.
void cell_base::clean_mod() {
  cell_0_mod.clear();
  cell_0_flag.clear();
  cell_1_mod.clear();
  cell_1_flag.clear();
  cell_2_mod.clear();
  cell_2_flag.clear();
  cell_3_mod.clear();
  cell_3_flag.clear();
}
*/
// ---------------------------------------------------------
// Print information.
// ---------------------------------------------------------

void cell_base::print_conn(int dim, int tag) {
  ent_conn e_con;
  if (dim < 1 or dim > 3 or tag > n[dim]) {
    printf("Tag or dim not valid.\n");
  }

  else if (get_conn(dim, tag, &e_con)) {
    printf("%dcell %d has %d %dcells", dim, tag+1, (int)e_con.conn.size(), dim-1);
    for (int i = 0; i < e_con.conn.size(); i++) {
      printf(" %d", e_con.conn.at(i)+1);
    }
    printf(".\n");
  }
}

ent_conn cell_base::print_conn_e(int dim, int tag) {
  ent_conn e_con;
  if (dim < 1 or dim > 3 or tag > n[dim]) {
    printf("Tag or dim not valid.\n");
  }

  get_conn(dim, tag, &e_con);
  return e_con;
}

ent_conn cell_base::print_conn_dim(int dim, int tag, int dim_adj) {
  ent_conn e_con;
  if (dim < 1 or dim > 3 or tag > n[dim]) {
    printf("Tag or dim not valid.\n");
  }

  get_conn_dim(dim_adj, dim, tag, &e_con);
  return e_con;
}

// Print entities:
void cell_base::print_ent() {
  ent_conn e_con;
  printf("There are %d 0, %d 1, %d 2, %d 3 cells.\n", get_sz(0), 
                                    get_sz(1), get_sz(2), get_sz(3));
  for (int i = 0; i < get_sz(1); i++) {
    if (get_conn(1, i, &e_con))
      printf("1cell %d has %d 0cells", i+1, (int)e_con.conn.size());
      for (int j = 0; j < e_con.conn.size(); j++) {
        printf(" %d", e_con.conn.at(j)+1);
      }
      printf(".\n");
  }

  for (int i = 0; i < get_sz(2); i++) {
    if (get_conn(2, i, &e_con)) {
      printf("2cell %d has %d 1cells", i+1, (int)e_con.conn.size());
      for (int j = 0; j < e_con.conn.size(); j++) {
        printf(" %d", e_con.conn.at(j)+1);
      }
      printf(".\n");
    }
  }

  for (int i = 0; i < get_sz(3); i++) {
    if (get_conn(3, i, &e_con)) {
      printf("3cell %d has %d surfaces", i+1, (int)e_con.conn.size());
      for (int j = 0; j < e_con.conn.size(); j++) {
        printf(" %d", e_con.conn.at(j)+1);
      }
      printf(".\n");
    }
  }
}

/*
// Print cells to be modified:
void cell_base::print_mod() {

  for (int i = 0; i < cell_0_mod.size(); i++) {
    printf("0cell%d, to be ", cell_0_mod.at(i)+1);
    if (cell_0_flag.at(i) == -1)
      printf("destroyed.\n");
    else
      printf("replaced by %d.\n", cell_0_flag.at(i));
  }
  for (int i = 0; i < cell_1_mod.size(); i++) {
    printf("1cell%d, to be ", cell_1_mod.at(i)+1);
    if (cell_1_flag.at(i) == -1)
      printf("destroyed.\n");
    else
      printf("replaced by %d.\n", cell_1_flag.at(i));
  }
  for (int i = 0; i < cell_2_mod.size(); i++) {
    printf("0cell%d, to be ", cell_2_mod.at(i)+1);
    if (cell_2_flag.at(i) == -1)
      printf("destroyed.\n");
    else
      printf("replaced by %d.\n", cell_2_flag.at(i));
  }
  for (int i = 0; i < cell_3_mod.size(); i++) {
    printf("3cell%d, to be ", cell_3_mod.at(i)+1);
    if (cell_3_flag.at(i) == -1)
      printf("destroyed.\n");
    else
      printf("replaced by %d.\n", cell_3_flag.at(i)+1);
  }

}
*/

// ---------------------------------------------------------
// Output operations. write_dmg, extract agm.
// ---------------------------------------------------------

void cell_base::vd_write_dmg(const char* filename) {
  FILE* f = fopen(filename, "w");
  ent_conn e_con;
  int i;
  /* entity counts */
  fprintf(f, "%d %d %d %d\n", n[3], n[2], n[1], n[0]);
  /* bounding box */
  fprintf(f, "0 0 0\n");
  fprintf(f, "0 0 0\n");
  /* vertices */
  for (int j = 0; j < n[0]; j++) {
    fprintf(f, "%d 0 0 0\n", j+1);
  }
  /* edges */
  for (int j = 0; j < n[1]; j++) {
    fprintf(f, "%d ", j+1);
    get_conn(1, j, &e_con);
    // In our algorithm, there cannot be zero number of elements. Still keep it.
    if (e_con.conn.size() == 0)
      fprintf(f,"-42 -42\n");
    // LOOK. -42 may require further handling. Check if 1cell handling is 
    // appropriate.
    else if (e_con.conn.size() == 1)
      fprintf(f,"%d -42\n", e_con.conn.at(0)+1);
    else
      fprintf(f, "%d %d\n",
          e_con.conn.at(0)+1, e_con.conn.at(1)+1);    
  }
  /* faces */
  for (int j = 0; j < n[2]; j++) {
    /* from gmi_write_dmg: we're going to cheat a bit here and
       treat all edges as one giant loop */
    fprintf(f, "%d 1\n", j+1);
    get_conn(2, j, &e_con);

    fprintf(f, "%d\n", (int)e_con.conn.size());
    for (i = 0; i < e_con.conn.size(); ++i)
      fprintf(f, " %d 0\n", e_con.conn.at(i)+1);
  }

  /* regions */
  for (int j = 0; j < n[3]; j++) {
    /* from gmi_write_dmg: same sort of cheat, all faces are one shell */
    fprintf(f, "%d 1\n", j+1);
    get_conn(3, j, &e_con);

    fprintf(f, "%d\n", (int)e_con.conn.size());
    for (i = 0; i < e_con.conn.size(); ++i)
      fprintf(f, " %d 0\n", e_con.conn.at(i)+1);
  }
  fclose(f);
}

// TODO Destroy the previous mesh if there is any.
// Reload the mesh with the model. LOOK Might be better to remove this from 
// here in the future. Used in 1cell collapse ATM.
// TODO Not used, may be removed.
/*
void cell_base::vd_set_mdl(apf::Mesh2* m, const char* meshFile) {
  char modelFile[50] = "temp.dmg";
  vd_write_dmg(modelFile);
  m = apf::loadMdsFromGmsh(gmi_load(modelFile), meshFile);
}

void cell_base::vd_set_mdl_smb(apf::Mesh2* m, const char* meshFile) {
  char modelFile[50] = "temp.dmg";
  vd_write_dmg(modelFile);

  m = apf::loadMdsMesh(gmi_load(modelFile), meshFile);
}
*/
// ---------------------------------------------------------
// Constructor and destructor.
// ---------------------------------------------------------

// Copy constructor
cell_base& cell_base::operator=( const cell_base& other ) {

  load_flag = other.load_flag;
  mdl = other.mdl;
  file = other.file;

  printf("Copying cell complex.\n");
  n[0] = other.n[0];
  n[1] = other.n[1];
  n[2] = other.n[2];
  n[3] = other.n[3];

  max_sz = other.max_sz;

  //cell_1_sz.clear();
  //cell_2_sz.clear();
  //cell_3_sz.clear();

  cell_1_sz.resize(n[1]);
  cell_2_sz.resize(n[2]);
  cell_3_sz.resize(n[3]);

  //con_1.clear();
  //con_2.clear();
  //con_3.clear();

  con_1.resize(2*n[1]);
  con_2.resize(max_sz*n[2]);
  con_3.resize(max_sz*n[3]);

  for(int dim = 0; dim < 4; dim++) {
    cell_free.at(dim).reserve(n[dim]);
    cell_free.at(dim).resize(0);
  }

  printf("Created pointers, fill the lists.\n");

  cell_1_sz = other.cell_1_sz;
  con_1 = other.con_1;

  cell_2_sz = other.cell_2_sz;
  con_2 = other.con_2;

  cell_3_sz = other.cell_3_sz;
  con_3 = other.con_3;

  for(int dim = 0; dim < 4; dim++) {
    cell_free.at(dim) = other.cell_free.at(dim);
    c_fmap.at(dim) = other.c_fmap.at(dim);
  }

/*
  for (int i = 0; i < n[1]; i++) {
    cell_1_sz.at(i) = other.cell_1_sz.at(i);
    for (int j = 0; j < cell_1_sz.at(i); j++) {
      con_1.at(2*i+j) = other.con_1.at(2*i+j);
    }
  }
  for (int i = 0; i < n[2]; i++) {
    cell_2_sz.at(i) = other.cell_2_sz.at(i);
    for (int j = 0; j < cell_2_sz.at(i); j++) {
      con_2.at(max_sz*i+j) = other.con_2.at(max_sz*i+j);
    }
  }
  for (int i = 0; i < n[3]; i++) {
    cell_3_sz.at(i) = other.cell_3_sz.at(i);
    for (int j = 0; j < cell_3_sz.at(i); j++) {
      con_3.at(max_sz*i+j) = other.con_3.at(max_sz*i+j);
    }
  }

  // Copy free cells.
  cell_0_free.resize(other.cell_0_free.size());
  for (int i = 0; i < cell_0_free.size(); i++) {
    cell_0_free.at(i) = other.cell_0_free.at(i);
  }
  cell_1_free.resize(other.cell_1_free.size());
  for (int i = 0; i < cell_1_free.size(); i++) {
    cell_1_free.at(i) = other.cell_1_free.at(i);
  }
  cell_2_free.resize(other.cell_2_free.size());
  for (int i = 0; i < cell_2_free.size(); i++) {
    cell_2_free.at(i) = other.cell_2_free.at(i);
  }
  cell_3_free.resize(other.cell_3_free.size());
  for (int i = 0; i < cell_3_free.size(); i++) {
    cell_3_free.at(i) = other.cell_3_free.at(i);
  }
*/
  fix_free();

  for(int i = 0; i < ext_map.size(); i++) {
    ext_map.at(i).clear();
  }
  ext_map.clear();
  ext_map.resize(3);

  for(int i = 0; i < ext_map.size(); i++) {
    ext_map.at(i) = other.ext_map.at(i);
  }
  cor_1c_map = other.cor_1c_map;
  cor_0c_map = other.cor_0c_map;

  ext_spur = other.ext_spur;

  return *this;
}


// Copy constructor
cell_base& cell_base::copy_local(cell_base& other, int c0_id ) {
  clear();

  load_flag = other.load_flag;
  mdl = other.mdl;
  file = other.file;

  ent_conn c_3;
  ent_conn c_2;
  ent_conn c_1;
  ent_conn c_0;
  other.get_conn_dim(3, 0, c0_id, &c_3);
  other.get_conn_lower(2, 3, &c_3, &c_2);
  other.get_conn_lower(1, 2, &c_2, &c_1);
  other.get_conn_lower(0, 1, &c_1, &c_0);

  printf("Copying cell complex.\n");
  std::vector<int>::iterator it;

  it = std::max_element(c_0.conn.begin(), c_0.conn.end());
  n[0] = *it+1;
  it = std::max_element(c_1.conn.begin(), c_1.conn.end());
  n[1] = *it+1;
  it = std::max_element(c_2.conn.begin(), c_2.conn.end());
  n[2] = *it+1;
  it = std::max_element(c_3.conn.begin(), c_3.conn.end());
  n[3] = *it+1;

  max_sz = other.max_sz;

  //cell_1_sz.clear();
  //cell_2_sz.clear();
  //cell_3_sz.clear();

  cell_1_sz.resize(n[1]);
  cell_2_sz.resize(n[2]);
  cell_3_sz.resize(n[3]);

  //con_1.clear();
  //con_2.clear();
  //con_3.clear();

  con_1.resize(2*n[1]);
  con_2.resize(max_sz*n[2]);
  con_3.resize(max_sz*n[3]);

  for(int dim = 0; dim < 4; dim++) {
    cell_free.at(dim).reserve(n[dim]);
    cell_free.at(dim).resize(0);
  }

  printf("Created pointers, fill the lists.\n");

  int c_curr = -1;
  for (int i = 0; i < c_3.conn.size(); i++) {
    c_curr = c_3.conn.at(i);
    cell_3_sz.at(c_curr) = other.cell_3_sz.at(c_curr);
    for (int j = 0; j < cell_3_sz.at(c_curr); j++) {
      con_3.at(max_sz*c_curr+j) = other.con_3.at(max_sz*c_curr+j);
    }
  }

  for (int i = 0; i < c_2.conn.size(); i++) {
    c_curr = c_2.conn.at(i);
    cell_2_sz.at(c_curr) = other.cell_2_sz.at(c_curr);
    for (int j = 0; j < cell_2_sz.at(c_curr); j++) {
      con_2.at(max_sz*c_curr+j) = other.con_2.at(max_sz*c_curr+j);
    }
  }

  for (int i = 0; i < c_1.conn.size(); i++) {
    c_curr = c_1.conn.at(i);
    cell_1_sz.at(c_curr) = other.cell_1_sz.at(c_curr);
    for (int j = 0; j < cell_1_sz.at(c_curr); j++) {
      con_1.at(2*c_curr+j) = other.con_1.at(2*c_curr+j);
    }
  }

  fix_free();

  if (get_free_sz(0) < 10)
    add_free(0, 10 - get_free_sz(0));
  if (get_free_sz(1) < 10)
    add_free(1, 10 - get_free_sz(1));
  if (get_free_sz(2) < 5)
    add_free(2, 5 - get_free_sz(2));

  for(int i = 0; i < ext_map.size(); i++) {
    ext_map.at(i).clear();
  }
  ext_map.clear();
  ext_map.resize(3);
  cor_0c_map.clear();
  cor_1c_map.clear();

  for(int i = 0; i < c_0.conn.size(); i++) {
    c_curr = c_0.conn.at(i);
    ext_map.at(0)[c_curr] = other.get_cell_ext(0, c_curr);
    cor_0c_map[c_curr] = other.get_0c_corner(c_curr);
  }

  for(int i = 0; i < c_1.conn.size(); i++) {
    c_curr = c_1.conn.at(i);
    ext_map.at(1)[c_curr] = other.get_cell_ext(1, c_curr);
    cor_1c_map[c_curr] = other.get_1c_corner(c_curr);
  }
  for(int i = 0; i < c_2.conn.size(); i++) {
    c_curr = c_2.conn.at(i);
    ext_map.at(2)[c_curr] = other.get_cell_ext(2, c_curr);
  }

  ext_spur = other.ext_spur;

  return *this;
}

void cell_base::clear_fix() {
/*
  for(int i = 0; i < fix_list.size(); i++) {
    fix_list.at(i).first.clear();
    fix_list.at(i).second.clear();
    fix_list.at(i).clear();
  }
*/
  fix_list.clear();
}

void cell_base::clear() {
  n[0] = 0;
  n[1] = 0;
  n[2] = 0;
  n[3] = 0;
  dummy_clear_stop();

  cell_1_sz.clear();
  cell_2_sz.clear();
  cell_3_sz.clear();

  con_1.clear();
  con_2.clear();
  con_3.clear();

  for(int dim = 0; dim < 4; dim++) {
    cell_free.at(dim).clear();
    c_fmap.at(dim).clear();
  }
/*
  cell_0_mod.clear();
  cell_1_mod.clear();
  cell_2_mod.clear();
  cell_3_mod.clear();

  cell_0_flag.clear();
  cell_1_flag.clear();
  cell_2_flag.clear();
  cell_3_flag.clear();
*/
  clear_fix();
  cor_1c_map.clear();
  cor_0c_map.clear();
  for(int i = 0; i < ext_map.size(); i++) {
    ext_map.at(i).clear();
  }
  ext_map.clear();
  ext_map.resize(3);
}

// Destructor
cell_base::~cell_base() {

  clear();
}

/*
// Adapted from gmi_base_read_dmg
// void vd_read_dmg(FILE* f)
struct cell_base cell_base::vd_read_dmg(FILE* f)
{
  int n[4];
  int i,k;
  int tag;

  std::cout<<"Reading entity counts..."<<std::endl;
  // read entity counts
  gmi_fscanf(f, 4, "%d %d %d %d", &n[3], &n[2], &n[1], &n[0]);
  printf("N are... V:%d E:%d S:%d R:%d.\n", n[0], n[1], n[2], n[3]);
  cell_base cell_comp(n[0],n[1],n[2],n[3]);
  ent_conn e_con;

  std::cout<<"Entity counts read..."<<std::endl;
  // bounding box
  gmi_fscanf(f, 0, "%*f %*f %*f");
  gmi_fscanf(f, 0, "%*f %*f %*f");

  // vertices 
  for (i = 0; i < n[0]; ++i) {
    int v1;
    gmi_fscanf(f, 1, "%d %*f %*f %*f", &v1);
    // std::cout<<"Vertex " << v1 <<std::endl;
  }

  std::cout<<"0-Cells read..."<<std::endl;
  // edges //
  for (i = 0; i < n[1]; ++i) {
    e_con.conn.size() = 2;
    gmi_fscanf(f, 3, "%d %d %d", &tag, &e_con.conn.at(0), &e_con.conn.at(1));
    printf("%d %d %d. ", tag-1, e_con.conn.at(0), e_con.conn.at(1));
    e_con.conn.at(0)--;
    e_con.conn.at(1)--;
    cell_comp.set_conn(1, tag-1, &e_con);
    cell_comp.print_conn(1, tag-1);
  }
  std::cout<<"1-Cells read..."<<std::endl;

  // faces //
  for (i = 0; i < n[2]; ++i) {
    gmi_fscanf(f, 1, "%d %*d", &tag);
    gmi_fscanf(f, 1, "%d", &e_con.conn.size());
    for (k = 0; k < e_con.conn.size(); ++k) {
      // tag, direction //
      gmi_fscanf(f, 1, "%d %*d", &e_con.conn.at(k));
      e_con.conn.at(k)--;
    }
    cell_comp.set_conn(2, tag-1, &e_con);
    cell_comp.print_conn(2, tag-1);
  }

  std::cout<<"2-Cells read..."<<std::endl;
  // regions //
  for (i = 0; i < n[3]; ++i) {
    gmi_fscanf(f, 1, "%d %*d", &tag);
    gmi_fscanf(f, 1, "%d", &e_con.conn.size());
    for (k = 0; k < e_con.conn.size(); ++k) {
      // tag, direction //
      gmi_fscanf(f, 1, "%d %*d", &e_con.conn.at(k));
      e_con.conn.at(k)--;
    }
    cell_comp.set_conn(3, tag-1, &e_con);
    cell_comp.print_conn(3, tag-1);
  }
  std::cout<<"3-Cells read..."<<std::endl;
  return cell_comp;

}
*/
