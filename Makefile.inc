# PYTHIA configuration file.
# Generated on Fri Aug 13 09:21:24 CST 2021 with the user supplied options:
# --with-root=/cvmfs/juno.ihep.ac.cn/sl6_amd64_gcc494/Pre-Release/J18v2r1-branch/ExternalLibs/ROOT/6.12.06
# --with-fastjet3=/hep/home/yili/Software/fastjet340

# Install directory prefixes.
PREFIX_BIN=/hep/home/yili/Software/pythia8306/bin
PREFIX_INCLUDE=/hep/home/yili/Software/pythia8306/include
PREFIX_LIB=/hep/home/yili/Software/pythia8306/lib
PREFIX_SHARE=/hep/home/yili/Software/pythia8306/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=/cvmfs/lhaaso.ihep.ac.cn/sl6_amd64_gcc494/gcc-install/bin/g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=/sdufs/home/yili/Software/fastjet340/bin/
FASTJET3_INCLUDE=-I/sdufs/home/yili/Software/fastjet340/include
FASTJET3_LIB=-L/sdufs/home/yili/Software/fastjet340/lib -Wl,-rpath,/sdufs/home/yili/Software/fastjet340/lib -lfastjet

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=false
LHAPDF6_CONFIG=
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/cvmfs/juno.ihep.ac.cn/sl6_amd64_gcc494/Pre-Release/J18v2r1-branch/ExternalLibs/ROOT/6.12.06/bin/
ROOT_INCLUDE=-I/cvmfs/juno.ihep.ac.cn/sl6_amd64_gcc494/Pre-Release/J18v2r1-branch/ExternalLibs/ROOT/6.12.06/include
ROOT_LIB=-L/cvmfs/juno.ihep.ac.cn/sl6_amd64_gcc494/Pre-Release/J18v2r1-branch/ExternalLibs/ROOT/6.12.06/lib -Wl,-rpath,/cvmfs/juno.ihep.ac.cn/sl6_amd64_gcc494/Pre-Release/J18v2r1-branch/ExternalLibs/ROOT/6.12.06/lib -lCore

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
