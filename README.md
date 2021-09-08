# VDLIB

VDLIB is an open source library for simulating grain growth. It is based on  [SCOREC](https://github.com/SCOREC/core), a massively parallelizable open source finite element library developed by Rensselaer Polytechnic Institute. The long term goal of VDLIB is modeling general microstructure evolution by inclusion of additional physics such as crystal plasticity. The user communicates with the objects of the vd_sim class through input files for specifying equations of motion, time integration scheme, and information to extract from the simulation. Example programs are provided in [VDtrials](https://github.com/erdemeren/VDtrials).

## Reference ##
More detailed information about the method is available on the arXiv:
[Topological transitions during grain growth on a finite element mesh](https://arxiv.org/abs/2101.12321)

## Dependencies ##
### MPI ###

### GSL ###

### SCOREC ###
This version is built with a specific commit of SCOREC and has not been updated to work with the later versions. The commit hash is provided below:

git clone -n https://github.com/SCOREC/core.git
git checkout c341ecd5bf908369d21cca095b067ce6ef4f4e19

## Contact ##

* This code base is developed by Erdem Eren, ereren@ucdavis.edu 
* Other community or team contact: 
* Coordinator: Jeremy Mason jkmason@ucdavis.edu
