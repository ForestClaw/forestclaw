## Makefile.am in applications/advection/2d/sphere

bin_PROGRAMS += applications/advection/2d/torus/torus

applications_advection_2d_torus_torus_SOURCES = \
	applications/advection/2d/torus/torus_user.cpp \
	applications/advection/2d/torus/torus.cpp \
	applications/advection/2d/torus/fclaw2d_map_cart.c \
	applications/advection/2d/torus/fclaw2d_map_annulus.c \
	applications/advection/2d/torus/fclaw2d_map_torus.c \
	applications/advection/2d/torus/fclaw2d_map_latlong.c \
	applications/advection/2d/torus/setprob.f \
	applications/advection/2d/torus/psi.f \
	applications/advection/2d/torus/user_4.6/tag4refinement.f \
	applications/advection/2d/torus/user_4.6/rpn2noncons.f \
	applications/advection/2d/torus/user_4.6/rpt2noncons.f \
	applications/advection/2d/torus/user_4.6/qinit_manifold.f \
	applications/advection/2d/torus/user_4.6/setaux_manifold.f \
	applications/advection/2d/torus/user_4.6/torus_output.f

applications_advection_2d_torus_torus_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_advection_2d_torus_torus_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
