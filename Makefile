# Locations of external dependencies
ROOT_DIR := /home/a1648400

PLIK_HI_L_FILE_DIR := $(ROOT_DIR)/staging/plc_2.0/hi_l/plik
# PLIK_HI_L_FILE_DIR := $(ROOT_DIR)/staging/plc_2.0/hi_l/plik_lite
PLIK_LOW_L_FILE_DIR := $(ROOT_DIR)/staging/plc_2.0/low_l/bflike
CLASS_PBH_FILE_DIR := $(ROOT_DIR)/staging/pbh_bsplines
HYREC_FILE_DIR := $(ROOT_DIR)/staging/hyrec_data

PLIK_DIR := $(ROOT_DIR)/plc-2.0
CLASS_DIR := $(ROOT_DIR)/class
MULTINEST_DIR := $(ROOT_DIR)/MultiNest_v3.10
CFITSIO_DIR := $(ROOT_DIR)/cfitsio
DIVER_DIR := $(ROOT_DIR)/diver

# Location of pc_multinest
PC_MULTINEST_DIR := $(PWD)
WORK_DIR := $(PC_MULTINEST_DIR)/build

MACHINE = intel

# Flags for the C++ compiler
ifeq ($(MACHINE),gnu)
CPPC := g++
OPT_FLAGS := -fopenmp -Wall -fPIC -O2 #-ffast-math -march=native
else
CPPC := icpc
OPT_FLAGS := -qopenmp -Wall -fPIC -O3 -march=native
endif
# turn off warnings using -w
INC_FLAGS := -I$(PC_MULTINEST_DIR)/include -I$(PC_MULTINEST_DIR)/params -I$(CLASS_DIR)/include -I$(PLIK_DIR)/include -I$(CFITSIO_DIR)/include -I$(DIVER_DIR)/include -I$(CLASS_DIR)/hyrec/include

BATCH_PLC_FLAGS = -DHAVE_PYEMBED=1 -DHAVE_PYTHON_H=1 -DHAS_LAPACK -DLAPACK_CLIK -DNOHEALPIX -DCLIK_LENSING -D'CLIKSVNVERSION="6dc2a8cf3965 MAKEFILE"' -DCAMSPEC_V1
ifeq ($(MACHINE),gnu)
BATCH_LIB_FLAGS = -L$(CFITSIO_DIR)/lib -L$(PLIK_DIR)/lib -L$(CLASS_DIR) -ldl -lpthread -lcfitsio -lclik -lclass
else
BATCH_LIB_FLAGS = -L$(CFITSIO_DIR)/lib -L$(PLIK_DIR)/lib -L$(CLASS_DIR) -ldl -lintlc -limf -lsvml -liomp5 -lifport -lifcoremt -lpthread -lcfitsio -lclik -lclass -lirc #-lclik_mkl
endif

# Flags for pc_propagate
CLASS_INC_FLAGS := -I$(CLASS_DIR)/cpp -I$(CLASS_DIR)/include
BATCH_CLASS_FLAGS = -Wl,-rpath,$(CLASS_DIR)

PC_MULTINEST_DEFS = #-g -DDBUG #-DLITE_HI_L #-DBAO_LIKE

# Flags for the Fortran compiler which compiles the .o files into the final binary when adding MultiNest
FC_LIBS = -L$(MULTINEST_DIR) -lnest3 -lstdc++ -L$(DIVER_DIR)/lib -ldiver
FC_FLAGS :=

ifeq ($(MACHINE),gnu)
FC := gfortran
else
FC := ifort
FC_LIBS += #-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lmkl_scalapack_ilp64 -lmkl_blacs_intelmpi_ilp64
endif

# Flags for the MPI compilers
ifeq ($(MACHINE),gnu)
FC_MPI := mpifort
FC_MPI_FLAGS := -ffree-line-length-none -DMPI
else
FC_MPI := mpifort
FC_MPI_FLAGS := -nofor-main -DMPI
endif

# CPPC_MPI := mpic++
CPPC_MPI := g++
CPPC_MPI_FLAGS := #-DMPI

# Own object files to link against
# PC_OBJS = PLCPack.o ClikObject.o ClikPar.o Param.o
PC_OBJS = pbh_io.o hyrec_io.o init_plc.o loglike.o multinest_loglike.o ClassEngine.o Engine.o
PC_DIVER_INC = pc_diver.h diver_loglike.h loglike.h

# CLASS object files to link against
CLASS_OBJS = ClassEngine.o Engine.o

ALL_OBJS = $(PC_OBJS) $(CLASS_OBJS)


########################################################################

PC_MULTINEST_FILES = -D'PLIK_HI_L_FILE_DIR="$(PLIK_HI_L_FILE_DIR)"' -D'PLIK_LOW_L_FILE_DIR="$(PLIK_LOW_L_FILE_DIR)"' -D'CLASS_PBH_FILE_DIR="$(CLASS_PBH_FILE_DIR)"' -D'HYREC_FILE_DIR="$(HYREC_FILE_DIR)"' -DHYREC

PC_BUILD_OBJS = $(addprefix $(WORK_DIR)/, $(ALL_OBJS))
PC_DIVER_LOC = $(addprefix include/, $(PC_DIVER_INC))

# Work directory building
.base:
	if ! [ -e $(WORK_DIR) ]; then mkdir $(WORK_DIR) ; fi;
	touch $(WORK_DIR)/.base

vpath %.cc $(PC_MULTINEST_DIR)/src
vpath %.o $(WORK_DIR)
vpath .base $(WORK_DIR)

# Various dependencies for targets

all: pc_multinest pc_multinest_mpi pc_speedtest

Engine.o:

ClassEngine.o: Engine.o

hyrec_io.o:

pbh_io.o:

init_plc.o: ClassEngine.o pbh_io.o hyrec_io.o

loglike.o: init_plc.o

multinest_loglike.o: loglike.o

test_diver.o: init_plc.o

test_multinest.o: multinest_loglike.o

pc_diver.o: init_plc.o

pc_multinest_mpi.o: multinest_loglike.o


# Compilation commands

%.o: %.cc .base
	cd $(WORK_DIR); $(CPPC) -c -o $*.o $< $(OPT_FLAGS) $(PC_MULTINEST_DEFS) $(PC_MULTINEST_FILES) $(INC_FLAGS)

# MultiNest-compatible program compiled with MPI
pc_multinest_mpi: pc_multinest_mpi.o $(PC_BUILD_OBJS)
	$(FC_MPI) $(FC_MPI_FLAGS) -o $@ $(addprefix $(WORK_DIR)/,$(notdir $^)) $(OPT_FLAGS) $(PC_MULTINEST_DEFS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS) $(FC_LIBS)
	rm -rf output/*

pc_multinest_mpi.o: ../pc_multinest.cc $(PLIK_DIR)/src/clik.c $(PC_INC) .base
	cd $(WORK_DIR); $(CPPC) -c -o $@ $< $(OPT_FLAGS) $(PC_MULTINEST_DEFS) $(PC_MULTINEST_FILES) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS)

pc_diver: pc_diver.o $(PC_BUILD_OBJS)
	$(CPPC) -o $@ $(addprefix $(WORK_DIR)/,$(notdir $^)) $(OPT_FLAGS) $(PC_MULTINEST_DEFS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS) $(FC_LIBS)
	rm -rf output/*

pc_diver.o: ../pc_diver.cc $(PLIK_DIR)/src/clik.c $(PC_DIVER_LOC) .base
	cd $(WORK_DIR); $(CPPC) -c -o $@ $< $(OPT_FLAGS) $(PC_MULTINEST_DEFS) $(PC_MULTINEST_FILES) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS)

test_diver: test_diver.o $(PC_BUILD_OBJS)
	$(CPPC) -o $@ $(addprefix $(WORK_DIR)/,$(notdir $^)) $(OPT_FLAGS) $(PC_MULTINEST_DEFS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS)

test_diver.o: ../test_diver.cc $(PLIK_DIR)/src/clik.c $(PC_DIVER_LOC) .base
	cd $(WORK_DIR); $(CPPC) -c -o $@ $< $(OPT_FLAGS) $(PC_MULTINEST_DEFS) $(PC_MULTINEST_FILES) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS)

test_multinest: test_multinest.o $(PC_BUILD_OBJS)
	$(FC_MPI) $(FC_MPI_FLAGS) -o $@ $(addprefix $(WORK_DIR)/,$(notdir $^)) $(OPT_FLAGS) $(PC_MULTINEST_DEFS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS) $(FC_LIBS)
	rm -rf output/*

test_multinest.o: ../test_multinest.cc $(PLIK_DIR)/src/clik.c $(PC_INC) .base
	cd $(WORK_DIR); $(CPPC) -c -o $@ $< $(OPT_FLAGS) $(PC_MULTINEST_DEFS) $(PC_MULTINEST_FILES) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS)

output_chain_files: output_chain_files.cc
	$(CPPC) -o $@ $< $(OPT_FLAGS) $(PC_MULTINEST_DEFS) $(INC_FLAGS)

clean:
	rm -rf build/* output/*
