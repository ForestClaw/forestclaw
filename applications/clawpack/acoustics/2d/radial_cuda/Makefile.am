## Makefile.am in applications/clawpack/acoustics/2d/radial_cuda

bin_PROGRAMS += applications/clawpack/acoustics/2d/radial_cuda/radial

applications_clawpack_acoustics_2d_radial_cuda_radial_SOURCES = \
	applications/clawpack/acoustics/2d/radial_cuda/radial.cpp \
	applications/clawpack/acoustics/2d/radial_cuda/radial_user.cpp \
	applications/clawpack/acoustics/2d/radial_cuda/radial_options.c \
	applications/clawpack/acoustics/2d/radial_cuda/radial_user.h \
	applications/clawpack/acoustics/2d/radial_cuda/fclaw2d_map_pillowdisk5.c \
	applications/clawpack/acoustics/2d/radial_cuda/setprob.f \
	applications/clawpack/acoustics/2d/radial_cuda/user_4.6/qinit.f \
	applications/clawpack/acoustics/2d/radial_cuda/user_5.0/qinit.f90 \
	applications/clawpack/acoustics/2d/radial_cuda/user_5.0/tag4refinement.f90 \
	applications/clawpack/acoustics/2d/radial_cuda/user_5.0/tag4coarsening.f90 \
	applications/clawpack/acoustics/2d/radial_cuda/user_5.0/setaux.f90 \
	applications/clawpack/acoustics/2d/rp/rpn2ac.f \
	applications/clawpack/acoustics/2d/rp/acoustics_user_fort.h \
	applications/clawpack/acoustics/2d/rp/rpt2ac.f \
	applications/clawpack/acoustics/2d/rp/rpn2_acoustics.f90 \
	applications/clawpack/acoustics/2d/rp/rpt2_acoustics.f90 \
	applications/clawpack/acoustics/2d/rp/rpn2_acoustics_mapped.f90 \
	applications/clawpack/acoustics/2d/rp/rpt2_acoustics_mapped.f90 \
	applications/clawpack/acoustics/2d/radial_cuda/user_cuda/qinit.f90 \
	applications/clawpack/acoustics/2d/radial_cuda/user_cuda/rpn2acoustics.cu

EXTRA_DIST += \
	applications/clawpack/acoustics/2d/radial_cuda/fclaw_options.ini \
	applications/clawpack/acoustics/2d/radial_cuda/regression.ini \
	applications/clawpack/acoustics/2d/radial_cuda/afterframe.m \
	applications/clawpack/acoustics/2d/radial_cuda/setplot2.m \
	applications/clawpack/acoustics/2d/radial_cuda/mapc2m.m \
	applications/clawpack/acoustics/2d/radial_cuda/map1d.m \
	applications/clawpack/acoustics/2d/radial_cuda/make_data.py \
	applications/clawpack/acoustics/2d/radial_cuda/make_plots.py \
	applications/clawpack/acoustics/2d/radial_cuda/1drad/setplot.py \
	applications/clawpack/acoustics/2d/radial_cuda/1drad/setrun.py


applications_clawpack_acoustics_2d_radial_cuda_radial_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

ACOUSTICS_RADIAL_DIR = applications/clawpack/acoustics/2d/radial_cuda

applications_clawpack_acoustics_2d_radial_cuda_radial_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CUDACLAW_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD) \
		$(ACOUSTICS_RADIAL_DIR)/user_cuda/devicelink.o


## nvcc -dlink has to be called after everything for the executable has been compiled.
## this generates a single object file with all of the linked  device code.
$(ACOUSTICS_RADIAL_DIR)/user_cuda/devicelink.o : \
	$(ACOUSTICS_RADIAL_DIR)/user_cuda/rpn2acoustics.o \
	src/solvers/fc2d_cudaclaw/libcudaclaw.la
		$(NVCC) $(CUDA_CFLAGS) -dlink  \
		$(ACOUSTICS_RADIAL_DIR)/user_cuda/rpn2acoustics.o \
		src/solvers/fc2d_cudaclaw/.libs/libcudaclaw.a \
		-o $(ACOUSTICS_RADIAL_DIR)/user_cuda/devicelink.o 
