# VDLIB

VDLIB is an open source library for simulating grain growth. It is based on  [SCOREC](https://github.com/SCOREC/core), a massively parallelizable open source finite element library developed by Rensselaer Polytechnic Institute. The long term goal of VDLIB is modeling general microstructure evolution by inclusion of additional physics such as crystal plasticity. The user communicates with the objects of the vd_sim class through input files for specifying equations of motion, time integration scheme, and information to extract from the simulation. Example programs are provided in [VDtrials](https://github.com/erdemeren/VDtrials).

## Headers ##

topo_cg:
Relax vertices using a conjugate gradient descent in the energy lanscape. Used to find equilibrium mesh profiles of stable microstructures.

topo_disp:
Wrapper for equations of motion.

topo_edisc:
Mesh level insertion operations.

topo_energy:
Currently not used. Class to define energy functions.

topo_entlist:
Reverse topology - entity lookup class.

topo_extinfo:
Collect mesh related information.

topo_feat:
Extract information to csv files.

topo_field:
Functions that interface with apf::Field to add fields to entities.

topo_geom:
Functions to extract or calculate geometric information.

topo_glens:
Generalized lens collapse operations.

topo_graph:
Graph based operations used to detect insertions.

topo_lens:
Edge, triangle, and tetrahedron splitters.

topo_ma:
Meshadapt related classes and functions.

topo_manip:
Equations of motion implementations and functions that change vertex positions.

topo_rand:
Generate randomized vectors.

topo_solvlin:
GSL implementation to solve Kuprats' equations of motion.

topo_tess:
Generate a mesh from input tess files. Assume convex elements. The SCOREC equivalent inverts elements.

topo_tmesh:
Classes used to extract topological information from the mesh.

topo_topo:
Describes the topological representation cell_base which is a misnomer as the stratification is not limited to cell complexes.

topo_vd:
Mainly contains vd_sim, the simulation interface.

topo_write:
Some of the file input and output functions.

## Reference ##
More detailed information about the method is available on the arXiv:
[Topological transitions during grain growth on a finite element mesh](https://arxiv.org/abs/2101.12321).

## Dependencies ##
### MPI ###
This version is built with [OpenMPI v. 3.1.2](https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.2.tar.gz).

### GSL ###
This version is built with [GSL v. 2.5](https://ftp.sotirov-bg.net/pub/mirrors/gnu/gsl/gsl-2.5.tar.gz).

### SCOREC ###
This version is built with a specific commit of SCOREC and has not been updated and is not guaranteed to work with the later versions. The commit hash is provided below:

```
git clone -n https://github.com/SCOREC/core.git
git checkout c341ecd5bf908369d21cca095b067ce6ef4f4e19
```

## License ##
VDlib is licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or https://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or https://opensource.org/licenses/MIT)

at your option.

## Contact ##

* This code base is developed by Erdem Eren, ereren@ucdavis.edu 
* Other community or team contact: 
* Coordinator: Jeremy Mason jkmason@ucdavis.edu
