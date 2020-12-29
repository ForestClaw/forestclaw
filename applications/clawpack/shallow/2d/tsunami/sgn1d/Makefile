# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = classic                  # Clawpack package to use
EXE = xclaw                         # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first

# Environment variable FC should be set to fortran compiler, e.g. gfortran

LAPACK_DIR=/Users/calhoun/software/LAPACK/local

# Compiler flags can be specified here or set as an environment variable
# -ffpe-trap=invalid,zero,overflow,underflow,denormal
# -ffpe-trap=invalid,zero,overflow      # conventional usage
# FFLAGS ?=  -O2 -g  -fbacktrace -ffpe-trap=invalid,zero,overflow -L${LAPACK_DIR}/lib -llapack -lblas
FFLAGS ?=  -O0 -g -ffpe-trap=invalid,zero,overflow  -L${LAPACK_DIR}/lib -llapack -lblas

# ---------------------------------
# package sources for this program:
# ---------------------------------

include $(CLAW)/classic/src/1d/Makefile.classic_1d


# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = setaux.f inlinelimiter.f bc1.f 

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------
MODULES = \

SOURCES = \
  qinit.f90 \
  setprob.f90 \
  bc1.f \
  ../swe/bathy.f90 \
  ../swe/geoclaw_riemann_utils.f \
  philim_user.f90 \
  inlinelimiter.f \
  src1.f90 \
  setaux.f90 \
  sgn.f90 \
  tri_solve.f90 \
  rp1_geoclaw.f90

# $(CLAW)/riemann/src/geoclaw_riemann_utils.f


#  ../all/rp1_shallow_bathymetry_fwave.f90
 



#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

