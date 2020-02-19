# Locations of external dependencies
ROOT_DIR := /home/a1648400

# PLIK_HI_L_FILE_DIR := $(ROOT_DIR)/staging/plc_2.0/hi_l/plik
PLIK_HI_L_FILE_DIR := $(ROOT_DIR)/staging/plc_2.0/hi_l/plik_lite
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
OUTPUT_DIR := $(PC_MULTINEST_DIR)/output

MACHINE = intel

# Flags for the C++ compiler
ifeq ($(MACHINE),gnu)
CPPC := g++
C_OPT_FLAGS := -std=c++11 -fopenmp -Wall -fPIC -O3# -ffast-math -march=native
else
CPPC := icpc
C_OPT_FLAGS := -std=c++11 -qopenmp -Wall -fPIC -O3 -march=native
endif

INC_FLAGS := -Iinclude -Iparams -I$(CLASS_DIR)/include -I$(PLIK_DIR)/include -I$(CFITSIO_DIR)/include -I$(DIVER_DIR)/include -I$(CLASS_DIR)/hyrec/include

BATCH_PLC_FLAGS = -DHAVE_PYEMBED=1 -DHAVE_PYTHON_H=1 -DHAS_LAPACK -DLAPACK_CLIK -DNOHEALPIX -DCLIK_LENSING -D'CLIKSVNVERSION="6dc2a8cf3965 MAKEFILE"' -DCAMSPEC_V1
ifeq ($(MACHINE),gnu)
BATCH_LIB_FLAGS = -L$(CFITSIO_DIR)/lib -L$(PLIK_DIR)/lib -L$(CLASS_DIR) -ldl -lpthread -lcfitsio -lclik -lclass_pbh
else
BATCH_LIB_FLAGS = -L$(CFITSIO_DIR)/lib -L$(PLIK_DIR)/lib -L$(CLASS_DIR) -ldl -lintlc -limf -lsvml -liomp5 -lifport -lifcoremt -lpthread -lcfitsio -lclik -lclass_pbh -lirc# -lclik_mkl
endif

PC_MULTINEST_DEFS = -DLITE_HI_L -g# -DDBUG# -DBAO_LIKE

# Flags for the Fortran compiler which compiles the .o files into the final binary when adding MultiNest
FC_LIBS = -L$(MULTINEST_DIR) -lnest3 -lstdc++ -L$(DIVER_DIR)/lib -ldiver

ifeq ($(MACHINE),gnu)
FC := gfortran
F_OPT_FLAGS := -fopenmp -fPIC -O3
else
FC := ifort
F_OPT_FLAGS := -qopenmp -fPIC -O3 -march=native
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

# Own object files to link against
PC_OBJS = init_plc.o pc_loglike.o multinest_loglike.o pbh_io.o hyrec_io.o

# CLASS object files to link against
CLASS_OBJS = ClassEngine.o Engine.o

ALL_OBJS = $(PC_OBJS) $(CLASS_OBJS)


########################################################################

PC_MULTINEST_FILES = -D'PLIK_HI_L_FILE_DIR="$(PLIK_HI_L_FILE_DIR)"' -D'PLIK_LOW_L_FILE_DIR="$(PLIK_LOW_L_FILE_DIR)"' -D'CLASS_PBH_FILE_DIR="$(CLASS_PBH_FILE_DIR)"' -D'HYREC_FILE_DIR="$(HYREC_FILE_DIR)"' -DHYREC

PC_BUILD_OBJS = $(ALL_OBJS)

all: .base pc_multinest_mpi

# Work directory building
.base:
	if ! [ -e $(WORK_DIR) ]; then mkdir $(WORK_DIR) ; fi;
	if ! [ -e $(OUTPUT_DIR) ]; then mkdir $(OUTPUT_DIR) ; fi;
	touch $(WORK_DIR)/.base

vpath %.cpp src:tests
vpath %.o $(WORK_DIR)
vpath .base $(WORK_DIR)

# Compilation commands

%.o: %.cpp .base
	$(CPPC) $(C_OPT_FLAGS) -c $< -o $(addprefix $(WORK_DIR)/, $(notdir $*.o)) $(PC_MULTINEST_DEFS) $(PC_MULTINEST_FILES) $(INC_FLAGS)

# MultiNest-compatible program compiled with MPI
pc_multinest_mpi: pc_multinest_mpi.o $(PC_BUILD_OBJS)
	$(FC_MPI) $(FC_MPI_FLAGS) -o $@ $(addprefix $(WORK_DIR)/,$(notdir $^)) $(F_OPT_FLAGS) $(PC_MULTINEST_DEFS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS) $(FC_LIBS)

pc_multinest_mpi.o: pc_multinest.cpp .base
	$(CPPC) -c -o $(addprefix $(WORK_DIR)/, $@) $< $(C_OPT_FLAGS) $(PC_MULTINEST_DEFS) $(PC_MULTINEST_FILES) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS)

pc_diver: pc_diver.o $(PC_BUILD_OBJS)
	$(CPPC) -o $@ $(addprefix $(WORK_DIR)/,$(notdir $^)) $(C_OPT_FLAGS) $(PC_MULTINEST_DEFS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS) $(FC_LIBS)

test_diver: test_diver.o $(PC_BUILD_OBJS)
	$(CPPC) -o $@ $(addprefix $(WORK_DIR)/,$(notdir $^)) $(C_OPT_FLAGS) $(PC_MULTINEST_DEFS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS)

test_multinest: test_multinest.o $(PC_BUILD_OBJS)
	$(FC_MPI) $(FC_MPI_FLAGS) -o $@ $(addprefix $(WORK_DIR)/,$(notdir $^)) $(C_OPT_FLAGS) $(PC_MULTINEST_DEFS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS) $(FC_LIBS)

point_multinest: point_multinest.o $(PC_BUILD_OBJS)
	$(FC_MPI) $(FC_MPI_FLAGS) -o $@ $(addprefix $(WORK_DIR)/,$(notdir $^)) $(C_OPT_FLAGS) $(PC_MULTINEST_DEFS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS) $(FC_LIBS)

output_chain_files: output_chain_files.o
	$(CPPC) -o $@ $(addprefix $(WORK_DIR)/,$(notdir $<)) $(C_OPT_FLAGS) $(PC_MULTINEST_DEFS) $(INC_FLAGS)

clean:
	rm -rf $(WORK_DIR)

clean-output:
	rm -rf $(OUTPUT_DIR)
