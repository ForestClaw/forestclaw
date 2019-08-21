
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = amrclaw                  # Clawpack package to use
EXE = xamr                          # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?=  -O2

# ---------------------------------
# package sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
include $(AMRLIB)/Makefile.amr_2d

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = setaux.f90 flag2refine2.f90 qad.f

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

MODULES = \

SOURCES = \
  setprob.f \
  qinit.f \
  q0.f \
  fdisc.f \
  torus_setaux.f90 \
  fclaw_utils.f \
  flag2refine2.f90 \
  fclaw2d_metric_default_fort.f \
  qad_choice.f90 \
  qad_original.f \
  qad_modified.f \
  qad_new.f90 \
  rpn2qad_flux.f90 \
  rpn2cons_fw_manifold.f \
  rpt2cons_manifold.f \
  rpn2.f \
  rpt2.f \
  ../all/torus_velocity.f \
  ../all/torus_basis.f \
  ../all/mapping_terms.f \
  ../all/mapc2m_torus.f \
  ../all/cellave2.f

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

