## Makefile.am in applications/cudaclaw/acoustics/2d/radial

bin_PROGRAMS += applications/cudaclaw/acoustics/2d/radial/radial

applications_cudaclaw_acoustics_2d_radial_radial_SOURCES = \
	applications/cudaclaw/acoustics/2d/radial/radial.cpp \
	applications/cudaclaw/acoustics/2d/radial/radial_user.cpp \
	applications/cudaclaw/acoustics/2d/radial/radial_options.c \
	applications/cudaclaw/acoustics/2d/radial/radial_user.h \
	applications/cudaclaw/acoustics/2d/radial/user_cuda/qinit.f90 \
	applications/cudaclaw/acoustics/2d/radial/user_cuda/rpn2acoustics.cu

EXTRA_DIST += \
	applications/cudaclaw/acoustics/2d/radial/fclaw_options.ini \
	applications/cudaclaw/acoustics/2d/radial/regression.ini \
	applications/cudaclaw/acoustics/2d/radial/afterframe.m \
	applications/cudaclaw/acoustics/2d/radial/setplot2.m \
	applications/cudaclaw/acoustics/2d/radial/map1d.m \
	applications/cudaclaw/acoustics/2d/radial/make_data.py \
	applications/cudaclaw/acoustics/2d/radial/make_plots.py \
	applications/cudaclaw/acoustics/2d/radial/1drad/setplot.py \
	applications/cudaclaw/acoustics/2d/radial/1drad/setrun.py


applications_cudaclaw_acoustics_2d_radial_radial_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CUDACLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

ACOUSTICS_RADIAL_DIR = applications/cudaclaw/acoustics/2d/radial

applications_cudaclaw_acoustics_2d_radial_radial_LDADD = \
        $(LDADD) \
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
