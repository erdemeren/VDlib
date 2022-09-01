#          MOST USED FUNCTIONS AND OBJECTS           #

## MPI ##
From https://www.open-mpi.org/doc/v1.8/
Mpi libraries are already included in Scorec headers.

#### Starting function calls: ####
int MPI_Init(&argc,&argv)
Initializes the MPI execution environment.

#### Closing function calls: ####
int MPI_Finalize()
Terminates MPI execution environment. 

## SCOREC ##

The SCOREC functions are explained in http://scorec.rpi.edu/~dibanez/core/. The following are the ones that are most used in VDlib to interface with the SCOREC library.


### PCU ###

pcu.h

#### Starting function calls: ####

int PCU_Comm_Init()
Initializes the PCU library. 
  This function must be called by all MPI processes before calling any other 
  PCU functions. MPI_Init or MPI_Init_thread should be called before this 
  function. 

int PCU_Comm_Self()

Returns the communication rank of the calling thread.

When called from a non-threaded MPI process, this function is equivalent to MPI_Comm_rank(MPI_COMM_WORLD,rank).

Ranks are consecutive from 0 to $pt-1$ for a program with $p$ processes and $t$ threads per process. Ranks are contiguous within a process, so that the $t$ threads in process $i$ are numbered from $ti$ to $ti+t-1$.

#### Closing function calls: ####

int PCU_Comm_Free()

Frees all PCU library structures.

This function must be called by all MPI processes after all other calls to PCU, and before calling MPI_Finalize.
  
#### Thread synchronization function calls: ####
void PCU_Barrier()
Blocking barrier over all threads. 

void PCU_Thrd_Barrier()
Blocks all threads of a process until all have hit the barrier.

Allreduce operations that synchronizes an array over all threads.
This function must be called by all ranks at the same time. p must point to an array of n variables. After this call, p[i] will contain the total/maximum/average/... of all p[i]'s given by each rank.

void PCU_Max_Ints(int* p, size_t n)
Performs an Allreduce maximum of int arrays.

### APF ###
From the headers apf.h, apfMesh.h, apfMesh2.h, apfNumbering.h, apfVector.h, apfMatrix.h.

#### Mesh related ####

apf::Mesh2* m: 
Interface to a mesh part. All mesh related operations are interfaced through
this object.

##### Input and output operations: #####

Mesh2 * apf::loadMdsMesh (const char* modelFile, const char* meshfile)
Load an MDS mesh from files.

Mesh2 * apf::loadMdsFromGmsh(gmi_load(const char* modelFile), const char* meshFile)
Load an MDS mesh from Neper and gmsh compatible files.

void writeNative (const char *fileName)
Write the underlying mesh into a set of smb files.

apf::writeVtkFiles(const char *prefix, apf::Mesh* m);
Write a set of parallel VTK Unstructured Mesh files from an apf::Mesh with ASCII encoding. Used to transfer the mesh into a visual format.

##### Entity related: #####
apf::MeshEntity* e;
Pointer to an entity in the mesh hierarchy.

int type = apf::MeshEntity::Mesh::getType(apf::MeshEntity* e)
Get the type of the mesh entity. 0, 1, 2, and 4 for vertices, edges, triangles
and tetrahedra, respectively.

int apf::Mesh::typeDimension[type]
Get the dimension of the mesh entity. 0, 1, 2, and 3 for type 0, 1, 2 and 4.

apf::MeshEntity* makeOrFind (apf::Mesh2 *m, apf::ModelEntity *c, int type, apf::MeshEntity **down, BuildCallback *cb=0)
Returns already existing entities, or creates a new entity with the downward adjacencies down and attached to the topology c.

apf::Mesh2::destroy(apf::MeshEntity* e);
Destroy a mesh entity. 

apf::Vector3 apf::getLinearCentroid(apf::Mesh* m, apf::MeshEntity* e)
Get the average of the entity's vertex coordinates. This also works if given just a vertex, so its a convenient way to get the centroid of any entity.

##### Finalizing the mesh: #####
apf::destroyMesh(m);
Destroys an apf::Mesh. 

##### Reposition vertices: #####
apf::Mesh::getPoint(apf::MeshEntity* e, int node, apf::Vector3& point);
Returns the coordinates of a node on a mesh entity.

apf::Mesh::setPoint(apf::MeshEntity* e, int node, apf::Vector3& point);
Set the coordinates of a mesh vertex.

##### Iterators: #####
apf::MeshIterator* it
Iterator used in iterating over mesh entities of a given dimension.

apf::MeshIterator* apf::Mesh::begin(int dimension)
Begin iteration over elements of one dimension.

apf::MeshEntity* apf::Mesh::iterate(apf::MeshIterator* it)
Iterate and return the next entity.

void apf::Mesh::end(apf::MeshIterator* it)
End and destroy the iterator.

#### Containers, functions for finding related mesh entities ####
apf::Up up;
Statically sized container for upward adjacency queries. Capped at 400 entities for now. 

apf::Downward down
A mesh entity pointer pointing to an array for keeping downward adjacencies.

int m->getDownward(MeshEntity* e, int dimension, apf::MeshEntity** adjacent)
Get the downward adjacencies of dimension and store them in adjacent. Return the number of downward adjacencies.

m->getUp(apf::MeshEntity* e, apf::Up up)
Get the upward adjacencies of dimension d+1.

int i1 = findIn(apf::MeshEntity** a, int ent_nbr, apf::MeshEntity* b)
Find the entity b in array pointed by a of size ent_nbr. Returns -1 if there is no match.

apf::MeshEntity* apf::getEdgeVertOppositeVert(apf::Mesh2* m, apf::MeshEntity* edge, apf::MeshEntity* v1)
Returns the vertex opposite of v1 bounding edge.

#### Fields attached to the mesh entities ####

apf::Field* apf::createField(apf::Mesh* m, const char* name, int valueType, apf::FieldShape* shape)
Create a field from any builtin or user defined FieldShape. Can be used to store scalar or tensorial fields over entities.

apf::Field* apf::Mesh::findField(const char* name)
Lookup a field by its unique name.

apf::setScalar(apf::Field *f, apf::MeshEntity *e, int node, double value)
Set a nodal value of a scalar field. 

double getScalar(Field* f, MeshEntity* e, int node);
Get a nodal value of a scalar field. 
	
apf::setVector(apf::Field *f, apf::MeshEntity *e, int node, apf::Vector3 const &value)
Set the nodal value of a vector field. 

apf::getVector(apf::Field *f, apf::MeshEntity *e, int node, apf::Vector3 &value)
Get the nodal value of a vector field. 

void getMatrix(Field* f, MeshEntity* e, int node, Matrix3x3& value);
Get the nodal value of a matrix field. 

void setMatrix(Field* f, MeshEntity* e, int node, Matrix3x3 const& value);
Set the nodal value of a matrix field. 

#### Element related ####
##### MeshElement, and measurement of entities. #####
apf::MeshElement* createMeshElement(apf::Mesh* m, apf::MeshEntity* e);
Creates a mesh element over an entity. Must be followed by destroyMeshElement(MeshElement* e) before leaving the scope the element was created in.

apf::MeshEntity * apf::getMeshEntity (apf::MeshElement* me)
Retrieve the mesh entity associated with an apf::MeshElement. 

apf::measure(MeshElement* ee);
Measures the volume, area, or length of a Mesh Element. 

void destroyMeshElement(MeshElement* e);
Destroy the mesh element over an entity. 

#### Vectors and matrices ####

apf::Matrix3x3::Matrix3x3 ( double  a11, double  a12, double  a13, double  a21,double  a22, doublea23, double  a31, double  a32, double  a33)

3x3 matrix object, that is initialized component-wise.

apf::Matrix< 3, 3 > invert (apf::Matrix< 3, 3 > const &m)
Inverts a 3x3 matrix and returns the inverted form. 

apf::Vector3::Vector3 ( double  a, double  b, double  c ) 
3 vector object, that is initialized component-wise.

#### Model related ####
apf::ModelEntity* apf::Mesh::toModel(apf::MeshEntity* e)
Get geometric classification of the mesh entity e.

int apf::Mesh::getModelType(apf::ModelEntity* me)
Return the model entity dimension.

int apf::Mesh::getModelTag(apf::ModelEntity* me)
Get the dimension-unique model entity identifier.

apf::ModelEntity* apf::Mesh::findModelEntity(int type, int tag)
Get the model entity by dimension and identifier.

#### Numberings ####

apf::Numbering* tagnumbering;
Used in storing information attached to the entities.

apf::Numbering* apf::createNumbering(apf::Mesh *mesh, const char *name, apf::FieldShape *shape, int components)
Create a generally-defined Numbering. 

apf::destroyNumbering (apf::Numbering *n)
Destroy a Numbering. 

apf::number(apf::Numbering *n, apf::MeshEntity *e, int node, int component, int number)
Number a degree of freedom of numbering n.

apf::Mesh::removeNumbering(apf::Numbering *n)
Remove numbering n attached to the mesh.

### GMI ###
From the header gmi.h.

#### Model related ####

gmi::agm::agm_ent_type 
0 for Vertex, 1, 2 and 3 for region

gmi_model* mdl
The basic structure for all GMI models.

gmi_model* apf::Mesh::getModel()
Returns the pointer to the model atttached to the mesh.

#### Entities and iterators ####

gmi_ent *e
Pointer to a model entity

gmi_iter* mdl_it
Iterator over model entities

gmi_iter* gmi_begin(struct gmi_model* mdl, int dim);
Begin an iterator over model entities of one dimension.

gmi_ent* gmi_next(struct gmi_model* mdl, gmi_iter* mdl_it)
Dereference and then increment an interator. 

void gmi_end(struct gmi_model* mdl, gmi_iter* mdl_it);
Free an iterator.

void gmi_tag(struct gmi_model *m, struct gmi_ent *e)
Get the tag of a model entity.

##### Entity set #####

struct gmi_set* s
a set of model entities. Should be created only with gmi_make_set and always call gmi_free_set after receiving one.

void gmi_free_set (struct gmi_set *s)
Free a gmi_set

gmi_set* gmi_adjacent(struct gmi_model *m, struct gmi_ent *e, int dim)
Query model entity adjacencies 

### MA ###
ma.h

ma::Entity *ent
Same thing as a apf::MeshEntity*.

#### Metric for edge length: ####
maSize.h

ma::IsotropicFunction
User-defined Isotropic size function. 

ma::AnisotropicFunction
User-defined Anisotropic size function. 

Inheriting from this class and overloading getValue, one can create custom metrics for allowable edge length.

double ma::IsotropicFunction::getValue (ma::Entity *vert)
Get the desired element size at this vertex 

void ma::AnisotropicFunction::getValue ( apf::MeshEntity*  vert, apf::Matrix &  r, apf::Vector &  h) 
Get the size field at this vertex.

#### User Configuration for MeshAdapt: ####

maSolutionTransfer.h
------------------------------------------------------

ma::SolutionTransfer* s
User-defined solution transfer base. Use to transfer attached fields to adapted
mesh entities.

maInput.h
------------------------------------------------------

ma::Input* in 
User configuration for MeshAdapt. Look into ma::Input class for options.

ma::Input* ma::configure(ma::configure (ma::Mesh*  m, ma::IsotropicFunction * 
f, ma::SolutionTransfer *  s = 0)
Generate a configuration based on an Isotropic function. 

ma::Input* ma::configure(ma::configure (ma::Mesh *  m, ma::AnisotropicFunction * f, ma::SolutionTransfer*  s = 0)
Generate a configuration based on an anisotropic function. 

ma::Input* ma::configure (ma::Mesh*  m, apf::Field *  size, ma::SolutionTransfer *  s = 0)
Generate a configuration based on an isotropic field.

#### Running MeshAdapt: ####
maAdapt.h
------------------------------------------------------

void ma::adapt (ma::Mesh *m, ma::IsotropicFunction *f, ma::SolutionTransfer *s=0)
Adapt based on an isotropic function.

void ma::adapt (ma::Mesh *m, ma::AnisotropicFunction *f, ma::SolutionTransfer *s=0)
Adapt based on an anisotropic function.

void ma::adapt (ma::Input *in)
Adapt with custom configuration.

