# -*- Makefile -*-

include_HEADERS += -I$(AMR_WAVEPROP) -I$(AMR_WAVEPROP_FORTRAN)

solver_SOURCES += \
	src/solvers/amr_waveprop/amr_waveprop.cpp \
	src/solvers/amr_waveprop/fortran_source/limiter.f \
	src/solvers/amr_waveprop/fortran_source/src2.f \
	src/solvers/amr_waveprop/fortran_source/bc2.f \
	src/solvers/amr_waveprop/fortran_source/flux2.f \
	src/solvers/amr_waveprop/fortran_source/step2.f \
	src/solvers/amr_waveprop/fortran_source/ClawPatch2.f
