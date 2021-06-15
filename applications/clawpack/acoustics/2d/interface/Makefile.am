
## Makefile.am in applications/clawpack/acoustics_2d_interface

bin_PROGRAMS += applications/clawpack/acoustics/2d/interface/interface

applications_clawpack_acoustics_2d_interface_interface_SOURCES = \
	applications/clawpack/acoustics/2d/interface/interface.cpp \
	applications/clawpack/acoustics/2d/interface/interface_options.c \
	applications/clawpack/acoustics/2d/interface/interface_user.cpp \
	applications/clawpack/acoustics/2d/interface/interface_user.h \
	applications/clawpack/acoustics/2d/interface/p0.f \
	applications/clawpack/acoustics/2d/interface/setprob.f \
	applications/clawpack/acoustics/2d/interface/value_exceeds_th.f90 \
	applications/clawpack/acoustics/2d/interface/user_4.6/setaux.f \
	applications/clawpack/acoustics/2d/interface/user_4.6/qinit.f \
	applications/clawpack/acoustics/2d/interface/user_4.6/tag4refinement.f \
	applications/clawpack/acoustics/2d/interface/user_4.6/tag4coarsening.f \
	applications/clawpack/acoustics/2d/interface/user_5.0/setaux.f90 \
	applications/clawpack/acoustics/2d/interface/user_5.0/qinit.f90 \
	applications/clawpack/acoustics/2d/interface/user_5.0/tag4refinement.f90 \
	applications/clawpack/acoustics/2d/interface/user_5.0/tag4coarsening.f90 \
	applications/clawpack/acoustics/2d/rp/acoustics_user_fort.h \
	applications/clawpack/acoustics/2d/rp/rpn2acv.f \
	applications/clawpack/acoustics/2d/rp/rpt2acv.f \
	applications/clawpack/acoustics/2d/rp/rpn2_vc_acoustics.f90 \
	applications/clawpack/acoustics/2d/rp/rpt2_vc_acoustics.f90

EXTRA_DIST += \
	applications/clawpack/acoustics/2d/interface/fclaw_options.ini \
	applications/clawpack/acoustics/2d/interface/regression.ini \
	applications/clawpack/acoustics/2d/interface/afterframe.m \
	applications/clawpack/acoustics/2d/interface/setplot2.m

applications_clawpack_acoustics_2d_interface_interface_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_acoustics_2d_interface_interface_LDADD = \
        $(LDADD)  \
        $(FCLAW_CLAWPACK46_LDADD)  \
        $(FCLAW_CLAWPACK5_LDADD)  \
        $(FCLAW_CLAWPATCH_LDADD)  \
        $(FCLAW_LDADD)
