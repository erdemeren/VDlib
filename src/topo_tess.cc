
#include <cstdio>
#include <cstring>
#include <cassert>
#include <iostream>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>

#include <apfNumbering.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_mesh.h>

#include "topo_extinfo.h"
#include "topo_tess.h"

// Adapted from gmi_file.c
static int starts_with(char const* with, char const* s) {
  int lw;
  int ls;
  lw = strlen(with);
  ls = strlen(s);
  if (ls < lw)
    return 0;
  return strncmp(with, s, lw) == 0;
}

static void seek_marker(FILE* f, char const* marker) {
  char* line = 0;
  size_t linecap = 0;
  while (-1 != gmi_getline(&line, &linecap, f))
    if (starts_with(marker, line))
      return;
}

void read_0c_pos(std::vector<apf::Vector3>* pos, FILE* f) {
  int n;
  int tag;
  int i;
  int j;

  seek_marker(f, " **vertex");
  gmi_fscanf(f, 1, "%d", &n);

  pos->clear();
  pos->resize(n);

  float fz[3] = {0,0,0};
  double z[3] = {0,0,0};
  //std::cout << "0c positions:" << std::endl;
  for (i = 0; i < n; ++i) {
    apf::Vector3 v;
    gmi_fscanf(f, 4, "%d %f %f %f %*d", &tag, &fz[0], &fz[1], &fz[2]);
    // It doesn't make sense to cast float to double, but I couldn't find
    // the correct format for reading double using vfscanf...
    z[0] = (double)fz[0];
    z[1] = (double)fz[1];
    z[2] = (double)fz[2];
    pos->at(tag-1).fromArray(z);
    //std::cout << "0c" << i+1 << " " << pos->at(tag-1) << std::endl;
  }
}
// gmi_file.c part ends.

void vd_create_tess(apf::Mesh2* m, const char* modelFile) {

  struct gmi_model* mdl = m->getModel();

  std::vector<apf::Vector3> pos_0c(0);
  std::map<int, apf::Vector3> pos_1c{};
  std::map<int, apf::Vector3> pos_2c{};

  FILE* fi = fopen(modelFile, "r");
  read_0c_pos(&pos_0c, fi);

  std::map<int, apf::MeshEntity*> c0_m{};
  //std::map<int, apf::MeshEntity*> c1_m;
  std::map<int, apf::MeshEntity*> c2_m{};
  apf::Downward dv;

  struct gmi_iter* it;
  struct gmi_ent* e;

  struct gmi_set* s;
  struct gmi_set* s_up;
  struct gmi_set* s_v;
  struct gmi_set* s_1c;

  double z[3] = {0,0,0};
  apf::Vector3 v_temp;

  it = gmi_begin(mdl, 0);

  int tag;

  while ((e = gmi_next(mdl, it))) {
    s = gmi_adjacent(mdl, e, 1);
    tag = gmi_tag(mdl, e);

    if (s->n != 0) {
      //std::cout << "0c" << tag << std::endl;
      apf::ModelEntity* vert_em = m->findModelEntity(0, tag);
      apf::MeshEntity* vert = m->createVert(vert_em);
      m->setPoint(vert, 0, pos_0c.at(tag-1));
      c0_m[tag] = vert;

      //std::cout << vert << " " << pos_0c.at(tag-1) 
      //          << " c0_m " << c0_m[tag]
      //          << std::endl;
    }
    gmi_free_set(s);
  }
  gmi_end(mdl, it);

  m->acceptChanges();

  // 1cells
  it = gmi_begin(mdl, 1);
  while ((e = gmi_next(mdl, it))) {
    tag = gmi_tag(mdl, e);
    s = gmi_adjacent(mdl, e, 0);
    s_up = gmi_adjacent(mdl, e, 2);
    if (s->n != 0) {
      assert(s->n == 2 and s_up->n != 0);

      v_temp.fromArray(z);
      //std::cout << "1c" << tag << " vertices: ";
      for(int i = 0; i < s->n; i++) {
        int c0_tag = gmi_tag(mdl, s->e[i]);
        dv[i] = c0_m[c0_tag];
        v_temp = v_temp + pos_0c.at(c0_tag-1);

        //std::cout << "0c" << c0_tag << "["<< dv[i] << "] ";
      }
      //std::cout << std::endl;

      v_temp = v_temp/s->n;
      pos_1c[tag] = v_temp;

      apf::ModelEntity* edge_em = m->findModelEntity(1, tag);
      apf::MeshEntity* edge = buildElement(m, edge_em, apf::Mesh::EDGE, dv);

      //std::cout << "1c" << tag << " " << edge << std::endl;
    }
    else {
      // Don't support unbounded 1-cells.
      assert(s_up->n == 0);
    }

    gmi_free_set(s);
    gmi_free_set(s_up);
  }
  gmi_end(mdl, it);
  m->acceptChanges();

  // 2cells

  it = gmi_begin(mdl, 2);

  while ((e = gmi_next(mdl, it))) {
    tag = gmi_tag(mdl, e);
    s = gmi_adjacent(mdl, e, 1);
    s_up = gmi_adjacent(mdl, e, 3);

    if (s->n != 0) {
      assert(s_up->n != 0);

      //std::cout << "2c" << tag << std::endl;

      v_temp.fromArray(z);
      for(int i = 0; i < s->n; i++) {
        v_temp = v_temp + pos_1c[gmi_tag(mdl, s->e[i])];
      }
      v_temp = v_temp/s->n;

      apf::ModelEntity* c2_em = m->findModelEntity(2, tag);
      apf::MeshEntity* vert = m->createVert(c2_em);
      m->setPoint(vert, 0, v_temp);

      c2_m[tag] = vert;
      pos_2c[tag] = v_temp;

      //std::cout << vert << " " << v_temp << std::endl;

      for(int i = 0; i < s->n; i++) {
        s_v = gmi_adjacent(mdl, s->e[i], 0);

        //std::cout << "\t1c " << gmi_tag(mdl, s->e[i]) << std::endl;
        //std::cout << "\tVertices ";
        for(int j = 0; j < s_v->n; j++) {
          int c0_tag = gmi_tag(mdl, s_v->e[j]);
          dv[j] = c0_m[c0_tag];
          //std::cout << "0c" << c0_tag << "["<< dv[j] << "] ";
        }
        dv[2] = vert;

        // Internal edges
        apf::Downward dv2;
        dv2[0] = vert;
        dv2[1] = dv[0];
        buildElement(m, c2_em, apf::Mesh::EDGE, dv2);
        dv2[1] = dv[1];
        buildElement(m, c2_em, apf::Mesh::EDGE, dv2);

        m->acceptChanges();
        //std::cout << "2c" << tag << "["<< dv[2] << "] " << std::endl;

        apf::MeshEntity* tri = buildElement(m, c2_em, 
                                    apf::Mesh::TRIANGLE, dv);
        //std::cout << "2c" << tag << " " << tri << std::endl;
        m->acceptChanges();

        gmi_free_set(s_v);
      }
    }
    else {
      // Don't support unbounded 1-cells.
      assert(s_up->n == 0);
    }

    gmi_free_set(s);
    gmi_free_set(s_up);
  }
  gmi_end(mdl, it);

  m->acceptChanges();

  // 3cells

  it = gmi_begin(mdl, 3);

  while ((e = gmi_next(mdl, it))) {
    tag = gmi_tag(mdl, e);
    s = gmi_adjacent(mdl, e, 2);

    if (s->n != 0) {
      //std::cout << "\n3c" << tag << std::endl;

      v_temp.fromArray(z);
      for(int i = 0; i < s->n; i++) {
        v_temp = v_temp + pos_2c[gmi_tag(mdl, s->e[i])];
      }
      v_temp = v_temp/s->n;

      apf::ModelEntity* c3_em = m->findModelEntity(3, tag);
      apf::MeshEntity* vert = m->createVert(c3_em);
      m->setPoint(vert, 0, v_temp);

      //std::cout << vert << " " << v_temp << std::endl;

      for(int i = 0; i < s->n; i++) {
        int c2_tag = gmi_tag(mdl, s->e[i]);
        //std::cout << "\t2c" << c2_tag << std::endl;
        s_1c = gmi_adjacent(mdl, s->e[i], 1);
        for(int j = 0; j < s_1c->n; j++) {
          int c1_tag = gmi_tag(mdl, s_1c->e[j]);
          //std::cout << "\t\t1c " << c1_tag << std::endl;

          s_v = gmi_adjacent(mdl, s_1c->e[j], 0);

          //std::cout << "\t\tVertices ";
          for(int k = 0; k < s_v->n; k++) {
            int c0_tag = gmi_tag(mdl, s_v->e[k]);
            dv[k] = c0_m[c0_tag];
            //std::cout << "0c" << c0_tag << "["<< dv[k] << "] ";
          }
          dv[2] = c2_m[c2_tag];
          //std::cout << "2c" << c2_tag << "["<< dv[2] << "] ";
          dv[3] = vert;
          //std::cout << "3c" << tag << "["<< dv[3] << "] ";
          //std::cout << std::endl;

          apf::Downward dv2;
          dv2[0] = vert;
          dv2[1] = dv[0];
          buildElement(m, c3_em, apf::Mesh::EDGE, dv2);
          dv2[1] = dv[1];
          buildElement(m, c3_em, apf::Mesh::EDGE, dv2);
          dv2[1] = dv[2];
          buildElement(m, c3_em, apf::Mesh::EDGE, dv2);

          m->acceptChanges();
          dv2[1] = dv[0];
          dv2[2] = dv[1];
          buildElement(m, c3_em, apf::Mesh::TRIANGLE, dv2);
          dv2[1] = dv[1];
          dv2[2] = dv[2];
          buildElement(m, c3_em, apf::Mesh::TRIANGLE, dv2);
          dv2[1] = dv[0];
          dv2[2] = dv[2];
          buildElement(m, c3_em, apf::Mesh::TRIANGLE, dv2);

          m->acceptChanges();

          dv2[0] = vert;
          dv2[1] = dv[0];
          dv2[2] = dv[1];
          dv2[3] = dv[2];
          //bool found3 = false;
          //for(int k = 0; k < 4; k++) {
          //  int c_dim = m->getModelType(m->toModel(dv2[k]));
          //  int c_id = m->getModelTag(m->toModel(dv2[k]));
          //  if(c_dim == 3) {
          //    assert(!found3);
          //    found3 = true;
          //  }
          //  std::cout << dv2[k] << " " << c_dim << "c" << c_id << ", ";
          //}
          //std::cout << std::endl;

          apf::MeshEntity* tet = buildElement(m, c3_em, 
                                    apf::Mesh::TET, dv2);
          //m->getDownward(tet, 0, dv);
          //std::cout << "tet " << tet << std::endl;
          //found3 = false;
          //for(int k = 0; k < 4; k++) {
          //  int c_dim = m->getModelType(m->toModel(dv[k]));
          //  int c_id = m->getModelTag(m->toModel(dv[k]));
          //  if(c_dim == 3) {
          //    assert(!found3);
          //    found3 = true;
          //  }
          //  std::cout << dv[k] << " " << c_dim << "c" << c_id << ", ";
          //}
          //std::cout << std::endl;
/*
          apf::Downward down;
          m->getDownward(tet, 0, down);
          for(int k = 0; k < 4; k++) {
            m->getPoint(down[k], 0, v_temp);
            //std::cout << "\t\t\t" << down[k] << " " << v_temp << std::endl;
          }
*/
          apf::MeshElement* ee = createMeshElement(m, tet);
          double meas = measure(ee);
          destroyMeshElement(ee);

          //std::cout << tet << " meas " << meas << std::endl;
          if(meas < 0.) {
            m->destroy(tet);
            apf::MeshEntity* e_temp = dv2[1];
            dv2[1] = dv2[2];
            dv2[2] = e_temp;
            //std::cout << "replaced tet " << tet;
            tet = buildElement(m, c3_em, apf::Mesh::TET, dv2);
            //m->getDownward(tet, 0, dv);

            //std::cout << "with " << tet << std::endl;
            //found3 = false;
            //for(int k = 0; k < 4; k++) {
            //  int c_dim = m->getModelType(m->toModel(dv[k]));
            //  int c_id = m->getModelTag(m->toModel(dv[k]));
            //  if(c_dim == 3) {
            //    assert(!found3);
            //    found3 = true;
            //  }
            //  std::cout << dv[k] << " " << c_dim << "c" << c_id << ", ";
            //}
            //std::cout << std::endl;
/*
            m->getDownward(tet, 0, down);
            for(int k = 0; k < 4; k++) {
              m->getPoint(down[k], 0, v_temp);
              std::cout << "\t\t\t" << down[k] << " " << v_temp << std::endl;
            }
*/
            ee = createMeshElement(m, tet);
            meas = measure(ee);
            destroyMeshElement(ee);

            //std::cout << tet << " meas " << meas << std::endl;
            assert(meas > 0);
          }

          gmi_free_set(s_v);
        }
        gmi_free_set(s_1c);
      }
    }
    else {
      // Don't support unbounded 1-cells for now. 
      assert(s_up->n == 0);
    }

    gmi_free_set(s);
  }
  gmi_end(mdl, it);

  m->acceptChanges();

  fclose(fi);
}

