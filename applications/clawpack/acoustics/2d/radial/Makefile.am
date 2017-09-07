## Makefile.am in applications/clawpack/acoustics_2d_radial

bin_PROGRAMS += applications/clawpack/acoustics/2d/radial/radial

applications_clawpack_acoustics_2d_radial_radial_SOURCES = \
	applications/clawpack/acoustics/2d/radial/radial_user.cpp \
	applications/clawpack/acoustics/2d/radial/radial_user.h \
	applications/clawpack/acoustics/2d/radial/radial.cpp \
	applications/clawpack/acoustics/2d/radial/fclaw2d_map_pillowdisk5.c \
	applications/clawpack/acoustics/2d/radial/setprob.f \
	applications/clawpack/acoustics/2d/radial/user_4.6/qinit.f \
	applications/clawpack/acoustics/2d/radial/user_5.0/qinit.f90 \
	applications/clawpack/acoustics/2d/radial/user_5.0/tag4refinement.f90 \
	applications/clawpack/acoustics/2d/radial/user_5.0/tag4coarsening.f90 \
	applications/clawpack/acoustics/2d/radial/user_5.0/setaux.f90 \
	applications/clawpack/acoustics/2d/rp/rpn2ac.f \
	applications/clawpack/acoustics/2d/rp/acoustics_user_fort.h \
	applications/clawpack/acoustics/2d/rp/rpt2ac.f \
	applications/clawpack/acoustics/2d/rp/rpn2_acoustics.f90 \
	applications/clawpack/acoustics/2d/rp/rpt2_acoustics.f90 \
	applications/clawpack/acoustics/2d/rp/rpn2_acoustics_mapped.f90 \
	applications/clawpack/acoustics/2d/rp/rpt2_acoustics_mapped.f90


applications_clawpack_acoustics_2d_radial_radial_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_clawpack_acoustics_2d_radial_radial_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
