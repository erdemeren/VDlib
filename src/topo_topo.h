#ifndef TOPO_TOPO_H
#define TOPO_TOPO_H

#include <iostream>
#include <assert.h>
#include <fstream>

#include "gmi_base.h"
#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>

// TODO cell insertion operations, should work with topo_graph
// TODO ent_conn insert element should check maximum size. If above, resize, 
// copy, etc.
// ---------------------------------------------------------
// Functions and objects to read(extract), manipulate, write(set) dmg(agm) 
// objects. Currently only works with dmg.
// ---------------------------------------------------------

// ---------------------------------------------------------
// ent_conn is a container for storing cell information.
// ---------------------------------------------------------

class c_elem {
  private:
  public:
    std::pair<int, int> first;
    std::pair<int, int> second;

    // Default constructor:
    c_elem() : first(-1,-1), second(-1,-1) {}

    // Copy constructor:
    c_elem(const c_elem &obj) : first(obj.first), second(obj.second) {}
    c_elem(const std::pair<std::pair<int, int>, std::pair<int, int> > &obj) : 
                                      first(obj.first), second(obj.second) {}
    c_elem& operator=( const c_elem& obj ) {
      first = obj.first;
      second = obj.second;
      return *this;
    }
    c_elem& operator=( const std::pair<std::pair<int, int>, 
                                              std::pair<int, int> > & obj ) {
      first = obj.first;
      second = obj.second;
      return *this;
    }
    bool operator==(const c_elem& obj)const {
      return (first == obj.first and second == obj.second);
    }
    bool operator==(const std::pair<std::pair<int, int>, 
                                       std::pair<int, int> > & obj)const {
      return (first == obj.first and second == obj.second);
    }
    ~c_elem() {}
};

// Remove list of entities from the ent_conn list. 
void rem_ents(std::vector<int>* conn, std::vector<int>* tags);

bool compc_elem(const c_elem &a, const c_elem &b);
void sort_celem(std::vector<c_elem >* c_set, int left = -1, int right = -1);
int partition_celem(std::vector<c_elem >* c_set, int left, int right);

class ent_conn {
  public:
    std::vector<int> conn;

    // Default constructor
    ent_conn() : conn(0) {
      conn.reserve(100);
      conn.clear();
    }

    ent_conn(int n_cap) : conn(0) {
      conn.reserve(n_cap);
      conn.clear();
    }
    // Copy constructor
    ent_conn(const ent_conn& that) : conn(that.conn) {
      //conn = that.conn;
    }
    // Copy
    ent_conn& operator=(const ent_conn& that) {
      conn = that.conn;
    }

    ent_conn& operator=(const std::vector<int>& that) {
      conn = that;
    }

    // ---------------------------------------------------------
    // Connection list operations. Check, add to the end, remove, remove first
    // ---------------------------------------------------------

    // Check if entity exists, return true if it does.
    bool chk_ent(int tag);
    int find_ent(int tag);

    // Add entity to the end of the removal list.
    void add_ent(int tag);
    // Remove entity with given tag from the removal list. 
    // Shift later entities.
    void rem_ent(int tag);

    // Entity indices higher than given tag are reduced by one.
    void shift_ind(int tag);

    // Try replacing tag1 with tag2. Assumes at most single occurance. 
    // Used for 1 cells. 
    bool repl_tag(int tag_new, int tag_old);

    // Quick Sort Functions for Ascending Order 
    // (2 Functions) 
    //http://codereview.stackexchange.com/questions/77782/quick-sort-implementation
    void quicksort(int left = -1, int right = -1);
    //Function to determine the partitions
    // partitions the array and returns the middle subscript
    int partition(int left, int right);

    // Check if two entity connection lists intersect. Return the index. If not 
    // return 0.
    int chk_intsct(ent_conn* e2, ent_conn* intsct);

    // Sort and keep unique entities.
    void unique();

    // Append the non existing elements of a second ent_conn list.
    void append(ent_conn* e2);

    void clear();
    void reserve(int n_new);
    void resize(int n_new);

    void print();

    // Destructor
    ~ent_conn() {
      //delete[] conn;
    }
};

// Cell modification related objects.
// As the mesh is modified, keep track of the cell modifications to make for 
// each mesh.
// They are modified externally and fed into the cell_base split_0cell object.

// This object is used to add a 2cell by splitting a 0cell, to connect two 
// 3cells that joined at the 0cell and were not connected before.
class cell_adder {
  private:
  public:
    cell_adder();
    void add_cell(int dim, int tag, std::vector<int>* new_con, 
                                            std::vector<int>* new_add);
    int dim;
    int tag;
    int tag_0;

    int min_sz;

    std::vector<std::vector<int> > e_new_id;
    // New 0-, 1-, 2-cells will have connection lists specified here.
    // Assuming 0-cell split, the connections will be new cells, as well.
    // [dim][index][conn]
    std::vector<std::vector<std::vector<int> > > e_new_con;

    // New 0-, 1-, 2-cells will be added to 1-, 2-, 3-cells. These are assumed
    // to include new cells.
    // [dim][index][conn] 
    std::vector<std::vector<std::vector<int> > > e_new_add;

    void clear();
};

class cell_base {
  private:

    // Store the source information. Used in reloading. 0 is gmi model, 1 is
    // file, -1 for none and cannot be reloaded.
    int load_flag;
    struct gmi_model* mdl;
    FILE* file;

    int n[4];
    // In order to simplify cell connection access structure, the cell sizes are
    // stored in an array. 1_cell structures are always of size 1 or 2.
    // A cell may also be empty or unused.
    //int max_sz;
    //std::vector<int> cell_1_sz;
    //std::vector<int> cell_2_sz;
    //std::vector<int> cell_3_sz;

    std::vector<std::vector<std::vector<int> > > con_up;
    std::vector<std::vector<std::vector<int> > > con;

    //std::vector<int> con_1; // Size is twice the number of 1cells.
    //std::vector<int> con_2; // Size is not bounded: A 2cell can have any number 
                            // of 1cells. Let's fix it at max_sz. The number of 
                            // 1cells a 2cell is connected to is stored in 
                            // 2_cell_sz. 
    //std::vector<int> con_3; // Similarly, a 3cell can have any number of 2cells.

    std::vector<std::map<int, bool> > ext_map;
    bool ext_spur;
    // Corner 1cells, bounding a single 3cell. 
    std::map<int, bool> cor_1c_map;
    // 0cells adjacent to corner 1cells.
    std::map<int, bool> cor_0c_map;

    // Free cell lists. These are LIFO vectors.
    // Maps for faster checking rather than vector search. The vectors store the
    // the useable ids.
    // TODO Better make these vectors of maps and vectors of vectors.
    std::vector<std::map<int, bool> > c_fmap;

    std::vector<std::vector<int> > cell_free;

    cell_base(const cell_base& that);

  public:
    // List of cells to fix and the cells they turn into.
    std::vector<c_elem > fix_list;

    // ---------------------------------------------------------
    // Empty constructor. Takes as input the sizes of each cell dimension.
    // File constructor.
    // agm constructor.
    // ---------------------------------------------------------

    cell_base(int n0, int n1, int n2, int n3);

    void load_file(FILE* f);

    cell_base(FILE* f);
    cell_base(const char* modelFile);

    cell_base(struct gmi_model* m);
    cell_base(apf::Mesh2* m);

    void reload_topo();
    void reload_topo(struct gmi_model* m);
    void reload_topo(FILE* f);

    // ---------------------------------------------------------
    // Get sizes for connection and cell entity lists.
    // ---------------------------------------------------------

    int get_sz(int dim);

    int cell_sz(int dim, int tag);

    int get_conn_sz(int dim, int tag);

    // ---------------------------------------------------------
    // Connection list operations, get and set using ent_conn.
    // ---------------------------------------------------------

    // Set the connection list of a given entity. 
    bool set_conn(int dim, int tag, struct ent_conn* e_con);
    bool set_conn_up(int dim, int tag, struct ent_conn* e_con);

    // Set the connection list of a given entity. 
    bool set_conn_repl(int dim, int tag_new, int tag_old, struct ent_conn* e_con);

    // Get the connection list of a given entity. 
    bool get_conn(int dim, int tag, struct ent_conn* e_con);

    // Get the connection list of a given entity. 
    bool get_conn_up(int dim, int tag, struct ent_conn* e_con);

    // GMI versions of get_conn and set_conn. Indices are shifted by one.
    bool set_conn_gmi(int dim, int tag, struct ent_conn* e_con);
    bool set_conn_up_gmi(int dim, int tag, struct ent_conn* e_con);

    // Get the connection list of a given entity. 
    bool get_conn_gmi(int dim, int tag, struct ent_conn* e_con);

    // Get the given lower dimensional connection list of a given entity. 
    bool get_conn_lower(int dim_adj, int dim, int tag, 
                                                struct ent_conn* e_cover);
    bool get_conn_upper(int dim_adj, int dim, int tag, 
                                                struct ent_conn* e_cover);

    // Get the given lower dimensional connection list of a given set of same
    // dimensional cells. 
    bool get_conn_lower(int dim_adj, int dim, struct ent_conn* e_set, 
                                                struct ent_conn* e_cover);
    bool get_conn_upper(int dim_adj, int dim, struct ent_conn* e_set, 
                                                struct ent_conn* e_cover);

    bool get_conn_lower_gmi(int dim_adj, int dim, struct ent_conn* e_set, 
                                                struct ent_conn* e_cover);

    // Get the given dimensional connection list of a given entity. 
    bool get_conn_dim(int dim_adj, int dim, int tag, struct ent_conn* e_cover);

    bool get_conn_dim_gmi(int dim_adj, int dim, int tag, 
                                                  struct ent_conn* e_cover);

    // Get the list of 1cell neighbors of a 0cell, connected to a specific 
    // 3cell. 
    bool get_conn_13(int cell0, int cell3, struct ent_conn* e_adj,
                                           struct ent_conn* e_not);

    bool get_conn_12(int cell0, int cell2, struct ent_conn* e_adj);

    // Given a d-cell-cc check if the cell is bounded by 0cell0.
    bool get_conn_0(int cell0, int c_dim, int c_id);

    // Given a d_list-cell list, remove the cells in the list that do not 
    // bound or are not bounded by c_dim-cell-c_id.
    void rem_conn(int c_dim, int c_id, int d_list, struct ent_conn* e_adj);
    void rem_conn_gmi(int c_dim, int c_id, int d_list, struct ent_conn* e_adj);

    // Given two cells, check if one is bounded by the other.
    bool chk_conn_d(int c_dim1, int c_id1, int c_dim2, int c_id2);

    // Given two cells, check if one is bounded by the other.
    bool chk_conn_d_gmi(int c_dim1, int c_id1, int c_dim2, int c_id2);

    // Given two cells of the same dim, check if they have a joint boundary.
    bool chk_conn_joint(int c_dim, int c_id1, int c_id2, int dim_other = -1);

    // Given two cells of the same dim, check if they have a joint boundary.
    bool chk_conn_joint_gmi(int c_dim, int c_id1, int c_id2, int dim_other = -1);

    // Given two cells of the same dim, get the joint boundary strata of
    // given dim.
    ent_conn get_conn_joint(int c_dim, int c_id1, int c_id2, int dim_other = -1);
    ent_conn get_conn_joint_gmi(int c_dim, int c_id1, int c_id2, int dim_other = -1);

    // Given two cells of the same dim, check if they have a joint boundary.
    bool chk_conn_joint_gmi(int c_dim, int c_id1, int c_id2);

    // Given a d-cell-c mark the bounding cells  as adj_c and cells bounded by 
    // bounding 0cells not belonging to the bounding set as end_c.
    void get_end_gmi(int c_dim, int c_id, 
                            std::vector<std::map<int, bool> >& end_c, 
                            std::vector<std::map<int, bool> >& adj_c);


    // Given a 0cell, check if it is on the boundary.
    bool chk_cell_ext_gmi(int c_dim, int c_id);
    // Given a cell, check if it is on the boundary using the map.
    bool get_cell_ext(int c_dim, int c_id);
    bool get_cell_ext_gmi(int c_dim, int c_id);

    // Given a cell, check if it is on the boundary or is bounded by a cell on 
    // the boundary using the map.
    bool get_cell_bound_ext(int c_dim, int c_id);
    bool get_cell_bound_ext_gmi(int c_dim, int c_id);

    // Using the map, if 1cell is a corner 1cell, return true. 
    bool get_1c_corner_gmi(int c_id);
    bool get_1c_corner(int c_id);

    bool get_0c_corner_gmi(int c_id);
    bool get_0c_corner(int c_id);

    // Looking over the adjacencies, check if the cell is on the exterior. 
    bool chk_0cell_ext_gmi(int c_id);
    bool chk_1cell_ext_gmi(int c_id);
    bool chk_2cell_ext_gmi(int c_id);
    bool chk_2cell_ext(int c_id);


    // Given a 1cell, check if it is on the boundary corner.
    bool chk_1cell_ext_c(int c_id);

    // Get the list of d-cell neighbors of a 0cell. d-cell > 0 and d-cell < 4.
    void get_conn_0d(int cell0, int c_dim, struct ent_conn* e_adj);

    // Get the list of d-cell neighbors of a 0cell. d-cell > 0 and d-cell < 4.
    void get_conn_0d_gmi(int cell0, int c_dim, struct ent_conn* e_adj);

    // Given a 1cell and a list of 2cells, find the 2cells accessible by 
    // starting from the 1cell and going over 1 and 2cells, touching the 0cell.
    // Remove the disconnected 2cells.
    // void get_conn_12d_gmi(int cell0, int cell1, struct ent_conn* e_2cells);

    // Given a 2,3-cell, find the set of collapsing 1cells in order to reduce
    // all all bounding 2-cells to be bounded by at most three 1-cells. 
    // Fills them into the ent_conn list.
    // Used in cell collapse algorithm. 
    void get_col_1cell(int cell_dim, int cell_id, int tag_0cell, struct ent_conn* e_conn);

    // Calls the get_col_1cell by using the dmg notation.
    void get_col_1cell_gmi(int cell_dim, int cell_id, int tag_0cell, struct ent_conn* e_conn);

    // ---------------------------------------------------------
    // Removal and eject, maintainance of the free cell lists.
    // ---------------------------------------------------------

    // Eject the cell, naive removal. Updates the same dim cell list.
    bool ej_cell(int dim, int tag);

    // Remove the cell. Go over the upper cell connection list, check if tag 
    // exists. If so, remove tag. 
    bool rem_cell(int dim, int tag);

    // Replace cell1 with cell2. Go over the upper cell connection list, check 
    // if tag exists. If so, check if the replacement also exists. If not, 
    // replace tag. 
    bool repl_cell(int dim, int tag_new, int tag_old);

    // Extend the cell list of given dim. Add the new cells to the free list. 
    // This should be called when a given free cell list is empty.
    bool add_free(int dim, int add_sz);

    // Return true if the cell is in the free list.
    // Can be used to skip unused cells.
    bool is_free(int dim, int cell_id);

    // Return the index of the free cell of given dim. If empty, return -1.
    int get_free(int dim);

    // Return the index of the free cell of given dim. If empty, return -1.
    // If not, remove the index of the cell from free list.
    int use_free_gmi(int dim);
    int use_free(int dim);

    int get_free_sz(int dim);

    // Remove the desired free cell from the free cell list.
    void remove_free(int dim, int tag);

    // After removal from cell list, add the freed cell to the end of the free
    // list.
    bool free_cell(int dim, int tag);

    // If the cell contains no lower dimensional or upper dimensional 
    // adjacencies, return true.
    bool chk_free(int dim, int cell_id);
    // Going over the cell list, find the cells that are free. Move them to the
    // free list. This is to be called when initializing the mesh to 
    // synchronize the cell_base structure with the mesh geometryz.
    void fix_free();

    // ---------------------------------------------------------
    // Merge operations and collapse (merges two lower dim).
    // These update the necessary cell lists.
    // ---------------------------------------------------------

    // Collapse a given cell. Connected 0_cell are merged. All other lower 
    // dimensional cells are removed.
    void coll_cell(int dim, int tag, int tag_0cell);
    void coll_cell_gmi(int dim, int tag, int tag_0cell);

    // Collapse a given triangular 2cell across a given 1cell.
    void coll_2cell(int tag2cell, int tag1cell);

    // After collapse check and remove empty cells.

    void clean_cells();

    void clean_dim(int dim);

    // ---------------------------------------------------------
    // Insert operations.
    // ---------------------------------------------------------
    // Insert a 2-cell or 1-cell using the cell_adder object.
    void ins_cell(cell_adder* ca);

    // ---------------------------------------------------------
    // Determination of the exterior boundaries.
    // ---------------------------------------------------------
    // Find the external and internal boundaries.
    void collect_ext();
    void set_ext_spur(bool state);
    void set_ext(int c_dim, int c_tag, bool ext);
    void set_ext_gmi(int c_dim, int c_tag, bool ext);

    void set_cor_ext(int c_tag, bool ext);
    void set_cor_ext_gmi(int c_tag, bool ext);

    void set_cor_ext_0c(int c_tag, bool ext);
    void set_cor_ext_0c_gmi(int c_tag, bool ext);

    void clear_cor();
    void set_cor_ext(std::vector<int>* c1_in);
    void set_cor_ext_gmi(std::vector<int>* c1_in);

    void keep_vect_vect(std::vector<int> &v_mod, std::vector<int> const &v_keep);
    // ---------------------------------------------------------
    // Detection of spurious cells.
    // ---------------------------------------------------------
    // Merge two given cells. This is mainly used when a cell is spurious. 
    void merg_cell(int dim, int tag1, int tag2);

    // Check if the cell is spurious. A spurious cell has lower than allowed
    // adjacencies(upper or lower) or bounds two cells that bound same cells
    // or bounds the same cell twice.
    bool chk_spur_up(int dim, int tag);
    bool chk_spur(int dim, int tag);
    // A spurious cell will have accompanying lower and higher dimensional 
    // spurious cells. A spurious cell is removed by merging the two higher
    // dimensional adjacencies, if they exist. 
    // 0cell: 3 > 1cell adj
    // 1cell: 3 > 2cell adj
    // 2cell: Both 3cell adj are same.
    // Save the lower dimensional cells of the spurious cell.
    // After merging, check whether the lower dimensional adjacencies are 
    // spurious and fix.
    bool fix_spur(int dim, int tag, bool clean = true);

    // The fixes to spurious cells occur recursively, so the fix_list may contain
    // replacing cells that are replaced by others. Going over the fix list, 
    // replace the intermediate replacements with the final ones.
    void process_fix_list();

    // Check if a spurious cell can be inserted. If a 1cell adjacent to 0celltag
    // has a single 0cell adjacency (0celltag), a spurious insertion may take
    // place. The cell_graph needs to be constructed from the mesh, rather than
    // the cell_base object.
    bool chk_2c(int tag);
    bool chk_2c_gmi(int tag);

    // ---------------------------------------------------------
    // Cell generation operations.
    // ---------------------------------------------------------
    // Given a collection of 2cells, add a new 1cell to all lists, by splitting
    // the 0cell.
    bool add_1cell(int tag_0cell, struct ent_conn* two_cell);

    // ---------------------------------------------------------
    // Modification list operations. Mainly used in lens collapse.
    // ---------------------------------------------------------
    // These operations mainly do removal and replacement. We could have used 
    // higher level functions provided in merge and collapse section, but some
    // decisions (e.g. whether a 3cell is to be collapsed) are taken at lens 
    // level. So this is more straightforward.
    // ---------------------------------------------------------
    // Set a topological modification. 
    void add_mod(int dim, int tag, int flag);

    // Apply list of modifications, starting from above.
    void apply_mod();

    // Clean list of modifications.
    void clean_mod();

    // ---------------------------------------------------------
    // Print information.
    // ---------------------------------------------------------

    void print_conn(int dim, int tag);
    ent_conn print_conn_e(int dim, int tag);
    ent_conn print_conn_dim(int dim, int tag, int dim_adj);

    // Print entities:
    void print_ent();

    // Print cells to be modified:
    void print_mod();

    // ---------------------------------------------------------
    // Output operations. write_dmg, extract agm, set mesh model.
    // ---------------------------------------------------------
    void vd_write_dmg(const char* filename);
    void write_dot(const char* dotname);

    // void vd_ext_agm(struct gmi_base* m);
    void vd_set_mdl(apf::Mesh2* m, const char* meshFile);
    void vd_set_mdl_smb(apf::Mesh2* m, const char* meshFile);

    // ---------------------------------------------------------
    // Constructor and destructor.
    // ---------------------------------------------------------

    // Copy constructor
    cell_base& operator=( const cell_base& other );
    cell_base& copy_local( cell_base& other, int c0_id );

    void clear();
    void clear_fix();
    // Destructor
    ~cell_base();
};


/*
// Adapted from gmi_base_read_dmg
// void vd_read_dmg(FILE* f)
struct cell_base cell_base::vd_read_dmg(FILE* f);
*/

#endif
