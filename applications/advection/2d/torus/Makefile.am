## Makefile.am in applications/advection/2d/sphere

bin_PROGRAMS += applications/advection/2d/torus/torus \
		applications/advection/2d/torus/torthem

applications_advection_2d_torus_torus_SOURCES = \
	applications/advection/2d/torus/torus.cpp \
	applications/advection/2d/torus/torus_common.c \
	applications/advection/2d/torus/torus_user.cpp \
	applications/advection/2d/torus/fclaw2d_map_nomap.c \
	applications/advection/2d/torus/fclaw2d_map_cart.c \
	applications/advection/2d/torus/fclaw2d_map_annulus.c \
	applications/advection/2d/torus/fclaw2d_map_torus.c \
	applications/advection/2d/torus/fclaw2d_map_latlong.c \
	applications/advection/2d/torus/setprob.f \
	applications/advection/2d/torus/psi.f \
	applications/advection/2d/torus/user_4.6/rpn2noncons.f \
	applications/advection/2d/torus/user_4.6/rpt2noncons.f \
	applications/advection/2d/torus/user_4.6/qinit.f \
	applications/advection/2d/torus/user_4.6/setaux.f \
	applications/advection/2d/torus/user_4.6/setaux_manifold.f

applications_advection_2d_torus_torus_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_advection_2d_torus_torus_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)

applications_advection_2d_torus_torthem_SOURCES = \
	applications/advection/2d/torus/torthem.c \
	applications/advection/2d/torus/torus_common.c

## Due to C++ libraries we need to force C++ linking for a C main program
nodist_EXTRA_applications_advection_2d_torus_torthem_SOURCES = dummy.cxx

applications_advection_2d_torus_torthem_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_advection_2d_torus_torthem_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
