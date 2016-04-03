## Makefile.am in applications/advection/2d/swirl50

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/advection/2d/swirl50/swirl50

applications_advection_2d_swirl50_swirl50_SOURCES = \
	applications/advection/2d/swirl50/swirl50_user.cpp \
	applications/advection/2d/swirl50/swirl50_user.h \
	applications/advection/2d/swirl50/swirl50.cpp \
	applications/advection/2d/swirl50/fclaw2d_map_nomap.c \
	applications/advection/2d/swirl50/psi.f \
	applications/advection/2d/swirl50/user_4.6/tag4refinement.f \
	applications/advection/2d/swirl50/user_4.6/bc2.f \
	applications/advection/2d/swirl50/user_4.6/rpn2ad1.f \
	applications/advection/2d/swirl50/user_4.6/rpt2ad1.f \
	applications/advection/2d/swirl50/user_4.6/qinit.f \
	applications/advection/2d/swirl50/user_4.6/b4step2.f \
	applications/advection/2d/swirl50/user_4.6/setaux.f \
	applications/advection/2d/swirl50/user_4.6/setprob.f \
	applications/advection/2d/swirl50/user_4.6/swirl_output.f

applications_advection_2d_swirl50_swirl50_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK50_CPPFLAGS) $(FCLAW_DUMMY_CPPFLAGS)

applications_advection_2d_swirl50_swirl50_LDADD = \
        $(FCLAW_CLAWPACK50_LDADD) $(FCLAW_DUMMY_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fclaw2d_clawpack4.6/libclawpack5.0.la
