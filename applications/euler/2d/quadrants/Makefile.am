## Makefile.am in applications/euler/2d/quadrants

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/euler/2d/quadrants/quadrants

applications_euler_2d_quadrants_quadrants_SOURCES = \
	applications/euler/2d/quadrants/quadrants.cpp \
	applications/euler/2d/quadrants/quadrants_user.h \
	applications/euler/2d/quadrants/quadrants_user.cpp \
	applications/euler/2d/quadrants/fclaw2d_map_nomap.c \
	applications/euler/2d/quadrants/fclaw2d_map_fivepatch.c \
	applications/euler/2d/quadrants/setprob.f \
	applications/euler/2d/quadrants/user_4.6/qinit.f \
	applications/euler/2d/quadrants/user_5.0/qinit.f90 \
	applications/euler/2d/all/clawpack_user.h \
	applications/euler/2d/all/clawpack46_rpn2eu4.f \
	applications/euler/2d/all/clawpack46_rpt2eu.f \
	applications/euler/2d/all/clawpack5_rpn2eu4.f90 \
	applications/euler/2d/all/clawpack5_rpt2eu.f90

applications_euler_2d_quadrants_quadrants_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_euler_2d_quadrants_quadrants_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD)  $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
