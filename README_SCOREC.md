## Installation:
First follow the installation of SCOREC test meshes which are required for running SCOREC tests and documents. 
```
export installDIR=/desired/directory/
cd $installDIR
```

### SCORECMESH: Contains trial meshes by SCOREC. Some are missing during the tests.
```
git clone https://github.com/SCOREC/pumi-meshes.git scorecmesh
```

### SCORECDOCS:
```
git clone https://github.com/SCOREC/docs.git "scorecdocs"
cd scorecdocs
mkdir build
cd build
cmake ..
```
These documents may not be up to date. To compile the old PUMI user guide follow the steps (Interestingly, the *fig* folder is defined to be in the parent directory, so move the *fig* folder to the parent directory):
```
cd ..
mv ./userGuide/fig ./fig     
latexmk -r ./latexmkrc.linux
```
*PUMI.pdf* in the *out* subdirectory explains the basic mesh structure.

# Follow README_ParmetisZoltan.md
For installation instructions, please follow [README_ParmetisZoltan](README_ParmetisZoltan.md).

### Setting the environmental variables:

```
export MPIHOME=/openmpi/install/directory
% Prefix is the install directory for the binaries
export PREFIX=$installDIR/SCOREC
export PARMETIS_HOME=$installDIR/install # (Assuming ParMETIS and Zoltan are installed in $installDIR/install)
export ZOLTAN_HOME=$installDIR/install
```
# SCORECCORE
The particular commit can be downloaded as provided below:
```
export scorecsrc=$installDIR"/scoreccore"
sudo rm -rf $scorecsrc
wget https://github.com/SCOREC/core/archive/c341ecd5bf908369d21cca095b067ce6ef4f4e19.zip -O scoreccore.zip
unzip -d ./ scoreccore.zip
mv core-c341ecd5bf908369d21cca095b067ce6ef4f4e19 scoreccore
```

```
cd $scorecsrc
rm -rf build
mkdir build
cd build
```

#### Remove -Werror flag from "\$scorecsrc/build/cmake/bob.cmake"
#### If the number of entities exceeds 10000000, -DMDS_ID_TYPE=long
```
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_C_COMPILER=$MPIHOME/bin/mpicc -DCMAKE_CXX_COMPILER=$MPIHOME/bin/mpicxx -DMPIRUN=$MPIHOME/bin/mpirun -DCMAKE_C_FLAGS="-O2 -g -Wall" -DCMAKE_CXX_FLAGS="-O2 -g -Wall" -DMDS_ID_TYPE=long -DIS_TESTING=True  -DENABLE_ZOLTAN=ON  -DZOLTAN_LIBRARY=$ZOLTAN_HOME/lib/libzoltan.a  -DZOLTAN_INCLUDE_DIR=$ZOLTAN_HOME/include  -DPARMETIS_LIBRARY=$PARMETIS_HOME/lib/libparmetis.a -DMETIS_LIBRARY=$PARMETIS_HOME/lib/libmetis.a  -DPARMETIS_LIBRARY=$PARMETIS_HOME/lib/libparmetis.a  -DPARMETIS_INCLUDE_DIR=$PARMETIS_HOME/include -DMETIS_INCLUDE_DIR=$PARMETIS_HOME/include  -DCMAKE_BUILD_TYPE=Debug -DIS_TESTING=True -DMESHES=$installDIR/scorecmesh

make 

cd test
ctest -VV >> $scorecsrc/testlog.txt

cd ..
sudo make install
```
### Useful functions:
For some of the most commonly used SCOREC functions in VDlib, please visit [README_SCOREC_mostused](README_SCOREC_mostused.md).
