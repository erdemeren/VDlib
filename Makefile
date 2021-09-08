CC=mpicxx
COMPFLAGS=-g -lmpi -std=gnu++11 -gdwarf-2 -g3

SCORECDIR = <SCOREC DIRECTORY>
ZOL_PARMETISDIR = <ZOLTAN/PARMETIS DIRECTORY>
GSL_DIR = <GSL DIRECTORY>

LIBDIR = -I$(SCORECDIR)/include -I$(ZOL_PARMETISDIR)/include -I$(GSL_DIR)/include -L$(SCORECDIR)/lib -L$(ZOL_PARMETISDIR)/lib -L$(GSL_DIR)/lib
LIBS = -lma -lparma -lapf_zoltan -lmds -lapf -lgmi -llion -lcrv -lmth -lph -lsam -lspr -lpumi -lpcu -lzoltan -lparmetis -lmetis -lgsl -lgslcblas

SRC_DIR = src
TEST_DIR = src/test

#------------------------
#  Source code for current project
#------------------------ topo_vd.cc
SOURCES = topo_solvlin.cc topo_rand.cc topo_pca.cc topo_ma.cc topo_chull.cc topo_disp.cc topo_manip.cc topo_extinfo.cc topo_write.cc topo_geom.cc topo_feat.cc topo_cg.cc topo_field.cc topo_gmi.cc topo_topo.cc topo_entlist.cc topo_tmesh.cc topo_lens.cc topo_glens.cc topo_edisc.cc topo_graph.cc topo_energy.cc topo_tess.cc topo_vd.cc    
# topo_edisc.cc 
TEST_SOURCES = topo_rand.cc topo_disp.cc topo_manip.cc topo_extinfo.cc topo_write.cc topo_geom.cc topo_field.cc topo_gmi.cc topo_topo.cc topo_tmesh.cc topo_lens.cc topo_glens.cc topo_vd.cc

# This calls .cc.o:
OBJS = $(patsubst %.cc, $(SRC_DIR)/%.o, $(SOURCES))
TEST_OBJS = $(patsubst %.cc, $(TEST_DIR)/%.o, $(TEST_SOURCES))
OBJS2 = $(SOURCES:.cc=.o)

#---------------------------
# simpler
#---------------------------
.cc.o:
	$(CC) $(COMPFLAGS) $(LIBDIR) $(LIBS) -c $< -o $@

all: $(OBJS)
	mkdir -p lib
	ar rc ./lib/libvdtopo.a $(OBJS)	

test: $(OBJS) $(TEST_OBJS)
	ar rc ./lib/libvdtopo.a $(OBJS)	$(TEST_OBJS)	

clean:
	rm -fr lib/libvdtopo.a src/*.o src/test/*.o src/*~  src/test/*~

