#
#        FILE: Makefile
#      AUTHOR: Paul Bartholomew <ptb08@ic.ac.uk>
# DESCRIPTION: Makefile for Py4Incompact3D: runs tests and builds documentation.
#

DEFS = -DDOUBLE_PREC

FC = mpif90
FFLAGS = -cpp -O3 -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -fbacktrace -ffree-line-length-none
F2PY = python3 -m numpy.f2py --f90exec=$(FC)

OBJ_DIR = obj

SRC = ./py4incompact3d.f90

DEP = $(D2DIR)/decomp_2d.f90
OBJ := $(DEP:%.f90=%.o)
OBJ := $(patsubst $(D2DIR)/%.o,%.o,$(OBJ))
OBJ := $(addprefix $(OBJ_DIR)/,$(OBJ))


decomp2d: $(OBJ_DIR) $(SRC) $(OBJ)
	$(info $(OBJ))
	$(F2PY) --f90flags="$(FFLAGS) $(DEFS)" -I$(OBJ_DIR) -m $@ -c $(filter-out $(OBJ_DIR),$^)

$(OBJ_DIR)/%.o: %.f90
	$(FC) $(FFLAGS) -J$(OBJ_DIR) -c $^ -o $@

$(OBJ_DIR)/%.o: $(D2DIR)/%.f90
	$(FC) $(FFLAGS) -J$(OBJ_DIR) -c $^ -o $@ -I$(D2DIR)

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

all: test

test:
	python3 -m unittest discover

.PHONY: doc
doc:
	make -C doc/ latexpdf


.PHONY: clean
clean:
	rm -f *.mod *.o *.so
	rm $(OBJ)
