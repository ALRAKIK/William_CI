# --------------------------- #
#      Basic directories      #
# --------------------------- #
BDIR = bin
ODIR = obj
SDIR = src
MODULE_DIR = $(SDIR)/module

# --------------------------- #
#     Create directories      #
# --------------------------- #
$(shell mkdir -p $(BDIR) $(ODIR))

# --------------------------- #
#      Compiler and flags     #
# --------------------------- #

# ========================= #
#    BUILD CONFIGURATION
# ========================= #
BUILD ?= local

# ======================================== #
# Compiler and flags based on build type   #
# ======================================== #
ifeq ($(BUILD),cluster)
    FC = gfortran
    FFLAGS = -Wall -Wno-unused -Wno-unused-dummy-argument -O3 -lblas -g -fbacktrace -fcheck=all -fimplicit-none -lm -ffree-line-length-none -fopenmp -ftree-vectorize -ffast-math -frecursive
    MODDIR = -J$(ODIR) -I$(ODIR) -I/usr/local/include
    LIBS = -fopenmp src/lib/libquadpack.a src/lib/liblapack.a src/lib/librefblas.a -ltrexio
    EXECUTABLE = $(BDIR)/CI_pop
    $(info Building for CLUSTER environment)
else
    FC = gfortran
    FFLAGS = -Wall -Wno-unused -Wno-unused-dummy-argument -O3 -g -fbacktrace -fcheck=all -fimplicit-none -fopenmp -funroll-loops -ftree-vectorize -ffast-math -frecursive -Wno-conversion
    MODDIR = -J$(ODIR) -I$(ODIR) -I/usr/local/include
    LIBS = src/lib/libquadpack.a -lblas -llapack \
           -L/usr/local/lib -ltrexio \
           -lm -fopenmp
    EXECUTABLE = $(BDIR)/CI_pop
    $(info Building for LOCAL environment)
endif

# =================================================== #
MODULES_IN_ORDER = constants_module \
                   u_polynomials \
                   boys \
                   trexio_f \
                   Torus_init \
                   bessel_iv_log \
                   files \
                   precompute_t \
                   Tools \
                   Lookup_table \
                   1D_bessel_Nmax \
                   atom_basis \
                   keywords \
                   ERI_mod \
                   bessel_derivative \
                   filter \
                   heaviside \
                   unitcell \
                   functions \
                   K_blocks \
                   MCDE \
# =================================================== #

# Find all module files
MODULE_SRC = $(wildcard $(MODULE_DIR)/*.f90)
MODULE_OBJ = $(addprefix $(ODIR)/, $(addsuffix .o, $(MODULES_IN_ORDER)))

# Find all other source files excluding modules
OTHER_SRC = $(shell find $(SDIR) -name "*.f90" -not -path "$(MODULE_DIR)/*")
OTHER_OBJ = $(patsubst $(SDIR)/%.f90,$(ODIR)/%.o,$(OTHER_SRC))

# Main program file
MAIN_SRC = main.f90
MAIN_OBJ = $(patsubst %.f90,$(ODIR)/%.o,$(MAIN_SRC))

# All object files
ALL_OBJ = $(MODULE_OBJ) $(OTHER_OBJ) $(MAIN_OBJ)

# Default target
all: local

# Target for modules only
modules: $(MODULE_OBJ)

# Compile the executable
$(EXECUTABLE): $(OTHER_OBJ) $(MAIN_OBJ)
	$(FC) -o $@ $(ALL_OBJ) $(LIBS)

# Rule to compile module files
$(ODIR)/%.o: $(MODULE_DIR)/%.f90
	@mkdir -p $(ODIR)
	$(FC) $(MODDIR) -c $< -o $@ $(FFLAGS)

# Rule to compile other source files
$(ODIR)/%.o: $(SDIR)/%.f90
	@mkdir -p $(dir $@)
	$(FC) $(MODDIR) -c $< -o $@ $(FFLAGS)

# Rule to compile main program file
$(ODIR)/%.o: %.f90
	$(FC) $(MODDIR) -c $< -o $@ $(FFLAGS)

# Clean everything
clean:
	rm -rf $(ODIR) $(EXECUTABLE)
	mkdir -p $(ODIR) $(BDIR)

# Build targets
local: clean
	@$(MAKE) BUILD=local modules -j1
	@$(MAKE) BUILD=local $(EXECUTABLE) -j

cluster: clean
	@$(MAKE) BUILD=cluster modules -j1
	@$(MAKE) BUILD=cluster $(EXECUTABLE) -j

fast:
	@$(MAKE) BUILD=local modules -j1
	@$(MAKE) BUILD=local $(EXECUTABLE) -j

.PHONY: all modules clean local cluster