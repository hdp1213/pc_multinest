HOME_DIR := /coepp/cephfs/adl/harryp
BATCH_DIR := /data/harryp

PLIK_DIR := $(HOME_DIR)/plc-2.0
CLASS_DIR := $(BATCH_DIR)/class
MULTINEST_DIR := $(HOME_DIR)/MultiNest_v3.10

# Flags for the C++ compiler
CPPC := g++
# turn off warnings using -w
OPT_FLAGS := -O4 -ffast-math -fopenmp -w
INC_FLAGS := -I$(CLASS_DIR)/cpp -I$(CLASS_DIR)/include -I$(PLIK_DIR)/include
CLASS_INC_FLAGS := -I$(CLASS_DIR)/cpp -I$(CLASS_DIR)/include

BATCH_PLC_FLAGS = -DHAVE_PYEMBED=1 -DHAVE_PYTHON_H=1 -DHAS_LAPACK -DLAPACK_CLIK -m64 -Wl,-rpath,$(BATCH_DIR)/lib -Wl,-rpath,$(PLIK_DIR) -Wl,-rpath,$(CLASS_DIR)
BATCH_LIB_FLAGS = -llapack -lblas -ldl -lgfortran -lgomp -lclik -lcfitsio -L$(BATCH_DIR)/lib
BATCH_CLASS_FLAGS = -Wl,-rpath,$(CLASS_DIR)
BATCH_MPI_FLAGS = -lmpi

# Flags for the Fortran compiler which compiles the .o files into the final binary when adding MultiNest
FC := gfortran
FC_FLAGS := -ffree-line-length-none
FC_LIBS := -L$(MULTINEST_DIR) -lnest3 -lstdc++

# Flags for the MPI compilers
FC_MPI := mpifort
FC_MPI_FLAGS := -ffree-line-length-none -DMPI
CPPC_MPI := g++
CPPC_MPI_FLAGS := #-DMPI

# Own object files to link against
PC_OBJS = PLCPack.o ClikObject.o
PC_INC = loglike.h ClikPar.h

# CLASS object files to link against
CLASS_SOURCE = input.o background.o thermodynamics.o perturbations.o primordial.o nonlinear.o transfer.o spectra.o lensing.o
CLASS_TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o  evolver_ndf15.o arrays.o parser.o quadrature.o hyperspherical.o common.o
HYREC = hyrectools.o helium.o hydrogen.o history.o
CLASS_LIB = -L$(CLASS_DIR)
OUTPUT = output.o
CPP = ClassEngine.o Engine.o

# Ultra Giga Mega Apex Alpha Omega Pack of CLASS objects to link
CLASS_BUILD_OBJS = $(addprefix $(CLASS_DIR)/build/, $(CLASS_SOURCE) $(CLASS_TOOLS) $(HYREC) $(OUTPUT))
CLASS_CPP_OBJS = $(addprefix $(CLASS_DIR)/cpp/, $(CPP))

# Same thing, but with CLASS sources to compile myself (lame)
CPP_SRC = $(addprefix $(CLASS_DIR)/cpp/, $(addsuffix .cc, $(basename $(CPP))))

vpath %.cc $(CLASS_DIR)/cpp
vpath %.o $(CLASS_DIR)/cpp


# Various dependencies for targets

all: pc_multinest pc_multinest_mpi pc_speedtest

$(CLASS_DIR)/cpp/ClassEngine.o: $(CLASS_DIR)/cpp/Engine.o

PLCPack.o: ClikObject.o


# Compilation commands

%.o: %.cc
	$(CPPC) -c -o $*.o $< $(OPT_FLAGS) $(INC_FLAGS) $(CLASS_LIB)

# MultiNest-compatible program 
# Must link with gfortran!!!
pc_multinest: pc_multinest.o $(CLASS_CPP_OBJS) $(PC_OBJS)
	$(FC) $(FC_FLAGS) -o pc_multinest pc_multinest.o $(CLASS_BUILD_OBJS) $(CLASS_CPP_OBJS) $(PC_OBJS) $(OPT_FLAGS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS) $(FC_LIBS)
	rm -f output/*

pc_multinest.o: pc_multinest.cc $(PLIK_DIR)/src/clik.c $(CPP_SRC) $(PC_INC)
	$(CPPC) -c -o pc_multinest.o pc_multinest.cc $(OPT_FLAGS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS)

# MultiNest-compatible program compiled with MPI
pc_multinest_mpi: pc_multinest_mpi.o $(CLASS_CPP_OBJS) $(PC_OBJS)
	$(FC_MPI) $(FC_MPI_FLAGS) -o pc_multinest_mpi pc_multinest_mpi.o $(CLASS_BUILD_OBJS) $(CLASS_CPP_OBJS) $(PC_OBJS) $(OPT_FLAGS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS) $(FC_LIBS)
	rm -f output/*

pc_multinest_mpi.o: pc_multinest.cc $(PLIK_DIR)/src/clik.c $(CPP_SRC) $(PC_INC)
	$(CPPC) -c -o pc_multinest_mpi.o pc_multinest.cc $(OPT_FLAGS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS)

# Speedtest program
pc_speedtest: pc_speedtest.o $(CLASS_CPP_OBJS) $(PC_OBJS)
	$(FC) $(FC_FLAGS) -o pc_speedtest pc_speedtest.o $(CLASS_BUILD_OBJS) $(CLASS_CPP_OBJS) $(PC_OBJS) $(OPT_FLAGS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS) $(FC_LIBS)
	rm -f output/*

pc_speedtest.o: pc_speedtest.cc $(PLIK_DIR)/src/clik.c $(CPP_SRC) $(PC_INC)
	$(CPPC) -c -o pc_speedtest.o pc_speedtest.cc $(OPT_FLAGS) $(INC_FLAGS) $(BATCH_PLC_FLAGS) $(BATCH_LIB_FLAGS)

# Run CLASS quickly to calculate any derived parameters you want
pc_propagate: pc_propagate.o $(CLASS_CPP_OBJS)
	$(CPPC_MPI) $(CPPC_MPI_FLAGS) -o pc_propagate pc_propagate.o $(CLASS_BUILD_OBJS) $(CLASS_CPP_OBJS) $(OPT_FLAGS) $(CLASS_INC_FLAGS) $(BATCH_CLASS_FLAGS)

pc_propagate.o: pc_propagate.cc $(CPP_SRC)
	$(CPPC_MPI) $(CPPC_MPI_FLAGS) -c -o pc_propagate.o pc_propagate.cc $(OPT_FLAGS) $(CLASS_INC_FLAGS) $(BATCH_CLASS_FLAGS)

clean:
	rm -f *.o output/*
	rm -f $(CLASS_DIR)/cpp/*.o
