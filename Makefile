# Data directories
ROOT_DIR := /home/a1648400

HI_L_LIKE = full

ifeq ($(HI_L_LIKE),lite)
PLIK_HI_L_FILE_DIR := $(ROOT_DIR)/staging/plc_2.0/hi_l/plik_lite
else
PLIK_HI_L_FILE_DIR := $(ROOT_DIR)/staging/plc_2.0/hi_l/plik
endif
PLIK_LOW_L_FILE_DIR := $(ROOT_DIR)/staging/plc_2.0/low_l/bflike
CLASS_PBH_FILE_DIR := $(ROOT_DIR)/staging/pbh_bsplines
HYREC_FILE_DIR := $(ROOT_DIR)/staging/hyrec_data

PLIK_DIR := $(ROOT_DIR)/plc-2.0
CLASS_DIR := $(ROOT_DIR)/class
MULTINEST_DIR := $(ROOT_DIR)/MultiNest_v3.10
CFITSIO_DIR := $(ROOT_DIR)/cfitsio

# Project directories
APP_DIR := apps
SRC_DIR := src
BIN_DIR := bin
TEST_DIR := test
BUILD_DIR := build
INC_DIR := include

HEADEXT := hpp
SRCEXT := cpp

# Macros
MACROS := -D'PLIK_HI_L_FILE_DIR="$(PLIK_HI_L_FILE_DIR)"' -D'PLIK_LOW_L_FILE_DIR="$(PLIK_LOW_L_FILE_DIR)"' -D'CLASS_PBH_FILE_DIR="$(CLASS_PBH_FILE_DIR)"' -D'HYREC_FILE_DIR="$(HYREC_FILE_DIR)"' -DHYREC

ifeq ($(HI_L_LIKE),lite)
MACROS += -DLITE_HI_L
endif

all: .base $(BIN_DIR)/base_planck

.base:
	@if ! [ -e $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi;
	@if ! [ -e $(BUILD_DIR) ]; then mkdir $(BUILD_DIR); fi;
	@touch $(BUILD_DIR)/.base

vpath %.$(SRCEXT) $(APP_DIR):$(SRC_DIR):$(TEST_DIR)
vpath %.o $(BUILD_DIR)
vpath .base $(BUILD_DIR)

# Compilers
CC = icpc
FC = mpiifort

# Flags
OPTFLAGS = -O3
CCDBUG = -g -Wall -Wextra# -Wcast-align
FCDBUG = -g

CCFLAGS = -fPIC -std=c++11 -qopenmp $(OPTFLAGS) $(CCDBUG)
FCFLAGS = -fPIC -nofor-main -DMPI $(OPTFLAGS) $(FCDBUG)
LDFLAGS = -fPIC -std=c++11 $(OPTFLAGS) $(CCDBUG)

# Includes and libraries
INCLUDES = -I$(INC_DIR) -I$(CLASS_DIR)/include -I$(PLIK_DIR)/include -I$(CLASS_DIR)/hyrec/include
LIBS = -L$(PLIK_DIR)/lib -L$(CLASS_DIR) -liomp5 -lclik -lclass_pbh
FCLIBS := -L$(MULTINEST_DIR)/lib -lnest3 -lstdc++

# Object files to compile
SOURCES := $(notdir $(shell find $(SRC_DIR) -type f -name *.$(SRCEXT)))
HEADERS := $(notdir $(shell find $(INC_DIR) -type f -name *.$(HEADEXT)))
TESTS := $(notdir $(shell find $(TEST_DIR) -type f -name *.$(SRCEXT)))

OBJECTS := $(SOURCES:.$(SRCEXT)=.o)
TESTOBJECTS := $(TESTS:.$(SRCEXT)=.o)
MAIN = main.o

# Rules
%.o: %.$(SRCEXT) .base
	@echo "Compiling" $@ "..."
	@$(CC) $(CCFLAGS) $(INCLUDES) $(MACROS) -c -o $(BUILD_DIR)/$@ $<

$(BIN_DIR)/base_planck: $(MAIN) $(OBJECTS)
	@echo "Linking" $@ "..."
	@$(FC) $(FCFLAGS) -o $@ $(addprefix $(BUILD_DIR)/, $(notdir $^)) $(LIBS) $(FCLIBS)

$(BIN_DIR)/profile: profile.o $(OBJECTS)
	@echo "Linking" $@ "..."
	@$(CC) $(LDFLAGS) -o $@ $(addprefix $(BUILD_DIR)/, $(notdir $^)) $(LIBS)

$(BIN_DIR)/test_clikobject: test_ClikObject.o $(OBJECTS)
	@echo "Linking" $@ "..."
	@$(CC) $(LDFLAGS) -o $@ $(addprefix $(BUILD_DIR)/, $(notdir $^)) $(LIBS)

$(BIN_DIR)/test_classengine: test_ClassEngine.o $(OBJECTS)
	@echo "Linking" $@ "..."
	@$(CC) $(LDFLAGS) -o $@ $(addprefix $(BUILD_DIR)/, $(notdir $^)) $(LIBS)

$(BIN_DIR)/test_likelihood: test_likelihood.o $(OBJECTS)
	@echo "Linking" $@ "..."
	@$(CC) $(LDFLAGS) -o $@ $(addprefix $(BUILD_DIR)/, $(notdir $^)) $(LIBS)

$(BIN_DIR)/test_point: test_point.o $(OBJECTS)
	@echo "Linking" $@ "..."
	@$(CC) $(LDFLAGS) -o $@ $(addprefix $(BUILD_DIR)/, $(notdir $^)) $(LIBS)

test: $(BIN_DIR)/test_clikobject $(BIN_DIR)/test_classengine $(BIN_DIR)/test_likelihood $(BIN_DIR)/test_point

profile: $(BIN_DIR)/profile

clean:
	@echo "Cleaning project ..."
	@rm -rf $(BUILD_DIR)/ $(BIN_DIR)/*

.PHONY: all clean profile
