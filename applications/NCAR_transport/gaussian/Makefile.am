## Makefile.am in applications/NCAR_transport/gaussian

bin_PROGRAMS += applications/NCAR_transport/gaussian/gaussian

include_HEADERS += applications/NCAR_transport/gaussian/gaussian_user.h

applications_NCAR_transport_gaussian_gaussian_SOURCES = \
	applications/NCAR_transport/gaussian/gaussian.cpp \
	applications/NCAR_transport/gaussian/gaussian_user.cpp \
	applications/NCAR_transport/gaussian/setprob.f \
	applications/NCAR_transport/gaussian/user_4.6/qinit.f \
	applications/NCAR_transport/gaussian/user_5.0/qinit.f90 \
	applications/NCAR_transport/gaussian/user_4.6/tag4refinement.f \
	applications/NCAR_transport/gaussian/user_5.0/tag4refinement.f \
	applications/NCAR_transport/gaussian/user_4.6/tag4coarsening.f \
	applications/NCAR_transport/gaussian/user_5.0/tag4coarsening.f \
	applications/NCAR_transport/all/velocity_wind.f \
	applications/NCAR_transport/all/gaussian.f \
	applications/NCAR_transport/all/transport_user.h \
	applications/NCAR_transport/all/transport_user.cpp \
	applications/NCAR_transport/all/transport_options.h \
	applications/NCAR_transport/all/transport_domain.c \
	applications/NCAR_transport/all/transport_options.c \
	applications/NCAR_transport/all/user46_b4step2_manifold.f \
	applications/NCAR_transport/all/user5_b4step2_manifold.f90 \
	applications/NCAR_transport/all/fclaw2d_map_cubedsphere.c \
	applications/NCAR_transport/all/fclaw2d_map_pillowsphere.c \
	applications/advection/2d/all/user46_setaux_manifold.f \
	applications/advection/2d/all/user5_setaux_manifold.f90 \
	applications/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/advection/2d/all/clawpack5_rpt2adv_manifold.f90

applications_NCAR_transport_gaussian_gaussian_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_NCAR_transport_gaussian_gaussian_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD) $(FCLAW_CLAWPACK5_LDADD)
