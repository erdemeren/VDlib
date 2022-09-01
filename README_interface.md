## Headers ##

topo_cg:
Relax vertices using a conjugate gradient descent in the energy lanscape. Used to find equilibrium mesh profiles of stable microstructures.

topo_disp:
Mainly defines f_calc class, the wrapper for equations of motion.

topo_edisc:
Mesh level insertion operations.

topo_energy:
Currently not used. Class to define energy functions to be interfaced through f_calc. 

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

