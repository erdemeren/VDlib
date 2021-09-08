
#include <vector>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <PCU.h>

#include "topo_extinfo.h"
#include "topo_topo.h"

#include "topo_entlist.h"

vd_entlist::vd_entlist(apf::Mesh2* m_in, cell_base* c_in) :
  m(NULL), c_base(NULL), ent_count(0, std::vector<int >(0) ),
  e(0, std::vector<std::vector<std::vector<apf::MeshEntity* > > > (0,
       std::vector<std::vector<apf::MeshEntity* > > (0,
       std::vector<apf::MeshEntity* > (0) ) ) ) {

  m = m_in;
  c_base = c_in;
  refresh();
}

void vd_entlist::change_mesh(apf::Mesh2* m_in, cell_base* c_in) {
  m = m_in;
  c_base = c_in;

  refresh();
}

void vd_entlist::refresh() {

  if(m == NULL)
    return;

  clear_ent();

  e.resize(4);

  for (int d_cell = 0; d_cell < 4; d_cell++)
    e.at(d_cell).resize(c_base->get_sz(d_cell));

  // Entities cannot belong to a cell of lower dimension. So, each cell 
  // contains d_cell+1 entity lists; for each mesh entity of dimensions it can 
  // contain. 
  for (int d_cell = 0; d_cell < 4; d_cell++) {
    for (int tag_cell = 0; tag_cell < e.at(d_cell).size(); tag_cell++) {
      e.at(d_cell).at(tag_cell).resize(d_cell+1);
    }
  }
 
  // Count all entities across all dimensions.
  for (int d_ent = 0; d_ent < 4; d_ent++) {
    clear_count();
    ent_count.resize(4);
    for (int d_cell = d_ent; d_cell < 4; d_cell++) {
      ent_count.at(d_cell).resize(c_base->get_sz(d_cell));
      for (int i = 0; i < ent_count.at(d_cell).size(); i++) {
        ent_count.at(d_cell).at(i) = 0;
      }
    }
    apf::MeshIterator* it_e = m->begin(d_ent);
    apf::MeshEntity* ent;

    while ((ent = m->iterate(it_e))) {
      int type = m->getType(ent);
      int d = apf::Mesh::typeDimension[type];
      int c_dim = m->getModelType(m->toModel(ent));
      int c_tag = m->getModelTag(m->toModel(ent));
      //std::cout << d << "-dim" << ent << " " 
      //          << c_dim << "-cell" << c_tag << std::endl;
      ent_count.at(c_dim).at(c_tag-1) = ent_count.at(c_dim).at(c_tag-1) + 1;
    }

    m->end(it_e);

    // Iterate over all of the cells, reserve the ith dimension entity list.
    for (int d_cell = d_ent; d_cell < 4; d_cell++) {
      for (int tag_cell = 0; tag_cell < ent_count[d_cell].size(); tag_cell++) {
        e.at(d_cell).at(tag_cell).at(d_ent).
                          reserve(ent_count.at(d_cell).at(tag_cell));
        //std::cout << d_cell << "cell" << tag_cell+1 << "has " 
        //  << ent_count[d_cell][tag_cell] << " " << d_ent << "-dim entities" 
        //  << std::endl;
      }
    }
  }


  // Register all entities across all dimensions.
  for (int d_ent = 0; d_ent < 4; d_ent++) {
    apf::MeshIterator* it_e = m->begin(d_ent);
    apf::MeshEntity* ent;

    while ((ent = m->iterate(it_e))) {
      int type = m->getType(ent);
      int d = apf::Mesh::typeDimension[type];
      int c_dim = m->getModelType(m->toModel(ent));
      int c_tag = m->getModelTag(m->toModel(ent));

      //assert(d==d_ent);
      e.at(c_dim).at(c_tag-1).at(d_ent).push_back(ent);
    }

    m->end(it_e);
  }

}

// Given the entity dimension, and cell dimension and id, get the entities.
void vd_entlist::get_entities_gmi(int ent_dim, int cell_dim, int cell_id, 
      std::vector<apf::MeshEntity*>* es_ent) {
  es_ent->resize(e.at(cell_dim).at(cell_id-1).at(ent_dim).size());

  for (int i = 0; i < es_ent->size(); i++) {
    es_ent->at(i) = e.at(cell_dim).at(cell_id-1).at(ent_dim).at(i);
  }

}

void vd_entlist::clear() {
  clear_ent();
}
void vd_entlist::clear_count() {

  for (int i = 0; i < ent_count.size(); i++) {
    ent_count.at(i).clear();
  }
  ent_count.clear();
}

void vd_entlist::clear_ent() {

  std::cout << "Clearing " << ent_count.size() << " " << e.size() << std::endl;
  clear_count();

  for (int i = 0; i < e.size(); i++) {
    for (int j = 0; j < e.at(i).size(); j++) {
      for (int k = 0; k < e.at(i).at(j).size(); k++) {
          e.at(i).at(j).at(k).clear();
          assert((int)e.at(i).at(j).at(k).capacity() > -1);
          assert((int)e.at(i).at(j).at(k).size() == 0);
      }
      e.at(i).at(j).clear();
      assert((int)e.at(i).at(j).capacity() > -1);
      assert((int)e.at(i).at(j).size() == 0);
    }
    e.at(i).clear();
    assert((int)e.at(i).capacity() > -1);
    assert((int)e.at(i).size() == 0);
  }
  e.clear();
  assert((int)e.capacity() > -1);
  assert((int)e.size() == 0);
}

vd_entlist::~vd_entlist() {
  clear();
}

vd_entlist_v::vd_entlist_v() :
  m(NULL), c_base(NULL), ent_count(0, std::vector<int > (0) ),
  v_ctr(0),
  e(0, std::vector<std::vector<std::vector<apf::MeshEntity* > > > (0,
       std::vector<std::vector<apf::MeshEntity* > > (0,
       std::vector<apf::MeshEntity* > (0) ) ) ),
  es(0, std::vector<apf::MeshEntity* > (0) ) {
}

vd_entlist_v::vd_entlist_v(apf::Mesh2* m_in, apf::MeshEntity* v_in, cell_base* c_in) :
  m(NULL), c_base(NULL), ent_count(0, std::vector<int > (0) ),
  v_ctr(0),
  e(0, std::vector<std::vector<std::vector<apf::MeshEntity* > > > (0,
       std::vector<std::vector<apf::MeshEntity* > > (0,
       std::vector<apf::MeshEntity* > (0) ) ) ),
  es(0, std::vector<apf::MeshEntity* > (0) ) {

  m = m_in;

  v_ctr.clear();
  v_ctr.resize(1);
  v_ctr[0] = v_in;
  c_base = c_in;

  get_vert();
}

vd_entlist_v::vd_entlist_v(apf::Mesh2* m_in, 
                    std::vector<apf::MeshEntity*>* v_in, cell_base* c_in) :
  m(NULL), c_base(NULL), ent_count(0, std::vector<int >(0) ),
  v_ctr(0),
  e(0, std::vector<std::vector<std::vector<apf::MeshEntity* > > > (0,
       std::vector<std::vector<apf::MeshEntity* > > (0,
       std::vector<apf::MeshEntity* > (0) ) ) ),
  es(0, std::vector<apf::MeshEntity* > (0) ) {

  m = m_in;

  v_ctr.clear();
  v_ctr.reserve(v_in->size());
  std::copy(v_in->begin(), v_in->end(), 
             std::inserter(v_ctr, v_ctr.end() ));

  c_base = c_in;

  get_vert();
}

void vd_entlist_v::change_mesh(apf::Mesh2* m_in, apf::MeshEntity* v_in, cell_base* c_in) {
  m = m_in;
  c_base = c_in;

  v_ctr.clear();
  v_ctr.resize(1);
  v_ctr[0] = v_in;

  get_vert();
}

void vd_entlist_v::change_mesh(apf::Mesh2* m_in, 
                    std::vector<apf::MeshEntity*>* v_in, cell_base* c_in) {
  m = m_in;
  c_base = c_in;

  v_ctr.clear();
  v_ctr.reserve(v_in->size());
  std::copy(v_in->begin(), v_in->end(), 
             std::inserter(v_ctr, v_ctr.end() ));

  get_vert();
}

void vd_entlist_v::get_vert() {
  clear_ent();

  e.resize(4);

  for (int d_cell = 0; d_cell < 4; d_cell++)
    e.at(d_cell).resize(c_base->get_sz(d_cell));

  // Entities cannot belong to a cell of lower dimension. So, each cell 
  // contains d_cell+1 entity lists; for each mesh entity of dimensions it can 
  // contain. 
  for (int d_cell = 0; d_cell < 4; d_cell++) {
    for (int tag_cell = 0; tag_cell < e.at(d_cell).size(); tag_cell++) {
      e.at(d_cell).at(tag_cell).resize(d_cell+1);
    }
  }

  es.clear();
  es.resize(4);

  //std::cout << "Getting entities around the vertices ";
  //for (int i = 0; i < v_ctr.size(); i++) {
  //  std::cout << v_ctr.at(i) << " [" 
  //            << m->getModelType(m->toModel(v_ctr.at(i))) 
  //            << "c" << m->getModelTag(m->toModel(v_ctr.at(i))) << "] ";
  //}
  //std::cout << std::endl;

  es[0] = v_ctr;

  for (int i = 0; i < 3; i++) {
    vd_set_up(m, &es.at(i), &es.at(i+1));
  }
  // This shouldn't set down, as the edisc object assumes all vd_entlist_v 
  // entities to be bounded by center vertex.
  //for (int i = 2; i > -1; i--) {
  //  vd_set_down(m, &es[i+1], &es[i]);
  //}

  // Count all entities across all dimensions.
  for (int d_ent = 0; d_ent < 4; d_ent++) {
    clear_count();
    ent_count.resize(4);
    for (int d_cell = d_ent; d_cell < 4; d_cell++) {
      //std::cout << d_cell << "-dim cells: " 
      //          << c_base->get_sz(d_cell) << std::endl;
      ent_count.at(d_cell).resize(c_base->get_sz(d_cell));
      for (int i = 0; i < ent_count.at(d_cell).size(); i++) {
        ent_count.at(d_cell).at(i) = 0;
      }
    }

    for (int i = 0; i < es.at(d_ent).size(); i++) {
      int type = m->getType(es.at(d_ent).at(i));
      int d = apf::Mesh::typeDimension[type];
      //std::cout << d << "-dim" << es.at(d_ent).at(i) << " ";

      int c_dim = m->getModelType(m->toModel(es.at(d_ent).at(i)));
      int c_tag = m->getModelTag(m->toModel(es.at(d_ent).at(i)));
      //std::cout << c_dim << "-cell" << c_tag << std::endl;
      ent_count.at(c_dim).at(c_tag-1) = ent_count.at(c_dim).at(c_tag-1) + 1;
    }

    // Iterate over all of the cells, reserve the ith dimension entity list.
    for (int d_cell = d_ent; d_cell < 4; d_cell++) {
      for (int tag_cell = 0; tag_cell < ent_count.at(d_cell).size(); tag_cell++) {
        e.at(d_cell).at(tag_cell).at(d_ent).reserve(
                                        ent_count.at(d_cell).at(tag_cell));
        //std::cout << d_cell << "cell" << tag_cell+1 << "has " 
        //  << ent_count[d_cell][tag_cell] << " " 
        //  << d_ent << "-dim entities" << std::endl;
      }
    }
  }


  // Register all entities across all dimensions.
  for (int d_ent = 0; d_ent < 4; d_ent++) {
    for (int i = 0; i < es.at(d_ent).size(); i++) {
      int type = m->getType(es.at(d_ent).at(i));
      int d = apf::Mesh::typeDimension[type];
      int c_dim = m->getModelType(m->toModel(es.at(d_ent).at(i)));
      int c_tag = m->getModelTag(m->toModel(es.at(d_ent).at(i)));

      //std::cout << c_dim << "cell" << c_tag << " "
      //    << d_ent << " type " 
      //    << es.at(d_ent).at(i) << std::endl;

      e.at(c_dim).at(c_tag-1).at(d_ent).push_back(es.at(d_ent).at(i));
    }
  }
}

// Given the entity dimension, and cell dimension and id, get the entities.
void vd_entlist_v::get_entities_gmi(int ent_dim, int cell_dim, int cell_id, 
      std::vector<apf::MeshEntity* >* es_ent) {
  es_ent->resize(e.at(cell_dim).at(cell_id-1).at(ent_dim).size());
  for (int i = 0; i < es_ent->size(); i++) {
    es_ent->at(i) = e.at(cell_dim).at(cell_id-1).at(ent_dim).at(i);
  }

}


void vd_entlist_v::print() {
  std::cout << "Printing the cell entity lists around vertices ";
  for(int i = 0; i < v_ctr.size(); i++) {
    std::cout << v_ctr[i] 
              << " " << m->getModelType(m->toModel(v_ctr[i])) << "c"
              << m->getModelTag(m->toModel(v_ctr[i]))
              << std::endl;
  }

  for (int c_dim = 3; c_dim > -1; c_dim--) {
    for (int c_tag = 0; c_tag < c_base->get_sz(c_dim); c_tag++) {
      std::cout << c_dim << "cell"
                << c_tag+1 << std::endl;
      for (int d_ent = 0; d_ent < c_dim+1; d_ent++) {
        for (int i = 0; i < e.at(c_dim).at(c_tag).at(d_ent).size(); i++) {
          std::cout << "\t" << d_ent << "-dim " 
                  << e.at(c_dim).at(c_tag).at(d_ent).at(i) << ", ";
        }
      }
      std::cout << std::endl;

    }
  }
}

void vd_entlist_v::print_v_pos() {
  std::cout << "Printing the cell entity lists vertex positions ";
  for(int i = 0; i < v_ctr.size(); i++) {
    std::cout << v_ctr[i] 
              << " " << m->getModelType(m->toModel(v_ctr[i])) << "c"
              << m->getModelTag(m->toModel(v_ctr[i]))
              << std::endl;
  }

  for (int c_dim = 3; c_dim > -1; c_dim--) {
    for (int c_tag = 0; c_tag < c_base->get_sz(c_dim); c_tag++) {
      std::cout << c_dim << "cell"
                << c_tag+1 << std::endl;
      int d_ent = 0;
      for (int i = 0; i < e.at(c_dim).at(c_tag).at(d_ent).size(); i++) {
        apf::Vector3 pos;
        m->getPoint(e.at(c_dim).at(c_tag).at(d_ent).at(i), 0, pos);
        std::cout << "\t" << d_ent << "-dim " 
                << e.at(c_dim).at(c_tag).at(d_ent).at(i) << ", ";
        std::cout << pos << std::endl;
      }
      std::cout << std::endl;

    }
  }
}

void vd_entlist_v::clear() {

  v_ctr.clear();
  clear_ent();
}

void vd_entlist_v::clear_count() {

  for (int i = 0; i < ent_count.size(); i++) {
    ent_count.at(i).clear();
  }
  ent_count.clear();
}

void vd_entlist_v::clear_ent() {

  for (int i = 0; i < ent_count.size(); i++) {
    ent_count.at(i).clear();
  }
  ent_count.clear();

  for (int i = 0; i < e.size(); i++) {
    for (int j = 0; j < e.at(i).size(); j++) {
      for (int k = 0; k < e.at(i).at(j).size(); k++) {
        e.at(i).at(j).at(k).clear();
        assert((int)e.at(i).at(j).at(k).capacity() > -1);
        assert((int)e.at(i).at(j).at(k).size() == 0);
      }
      e.at(i).at(j).clear();
      assert((int)e.at(i).at(j).capacity() > -1);
      assert((int)e.at(i).at(j).size() == 0);
    }
    e.at(i).clear();
    assert((int)e.at(i).capacity() > -1);
    assert((int)e.at(i).size() == 0);
  }
  e.clear();
  es.clear();
}

vd_entlist_v::~vd_entlist_v() {
  clear();
}


