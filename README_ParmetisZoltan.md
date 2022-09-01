# ParMETIS and Zoltan installation

```
export SRC=$installDIR
export TGT=$installDIR"/install"
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$TGT
sudo rm -rf $TGT
mkdir $TGT
mkdir $SRC
cd $SRC
```
# ParMETIS

```
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz

sudo rm -rf parmetis-4.0.3
tar xzf parmetis-4.0.3.tar.gz && cd parmetis-4.0.3
```
# IF 64 bit
```
sed -i -e 's:IDXTYPEWIDTH 32:IDXTYPEWIDTH 64:' metis/include/metis.h
```

# Building ParMETIS

```
cd build
cmake ../ \
-DCMAKE_INSTALL_PREFIX=$TGT \
-DMETIS_PATH=$PWD/../metis \
-DGKLIB_PATH=$PWD/../metis/GKlib \
-DCMAKE_C_COMPILER=mpicc \
-DCMAKE_CXX_COMPILER=mpicxx 

make 
sudo make install
sudo cp ../metis/include/metis.h $TGT/include/.
sudo cp libmetis/libmetis.a $TGT/lib/.
```

# Zoltan 
(seems to check --enable-threads )
config/tac_arg_enable_option.m4:dnl will test for --enable-threads when configure is run.  If it is defined (and not set to "no")
But --enable-threads is unrecognized by the configure.

```
cd $SRC
wget https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.83.tar.gz

sudo rm -r Zoltan-3.83
tar xzf v3.83.tar.gz
cd Zoltan-3.83
```

'tar: Ignoring unknown extended header keyword `SCHILY.dev'' It is said these are related to MAC OS, and can be safely ignored: http://lifeonubuntu.com/tar-errors-ignoring-unknown-extended-header-keyword/

```
mkdir build 
cd build

../configure \
--prefix=$TGT \
--with-parmetis \
--with-parmetis-libdir=$TGT/lib \
--with-parmetis-incdir=$TGT/include \
--enable-mpi \
--disable-examples \
--with-gnumake \
--with-id-type=ulong \
CXXCPP='mpicxx -E ' \
CC=mpicc \
CXX=mpicxx

make everything 
sudo make install
```

