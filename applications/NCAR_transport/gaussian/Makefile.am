## Makefile.am in applications/NCAR_transport/gaussian

bin_PROGRAMS += applications/NCAR_transport/gaussian/gaussian

include_HEADERS += applications/NCAR_transport/gaussian/gaussian_user.h

applications_NCAR_transport_gaussian_gaussian_SOURCES = \
	applications/NCAR_transport/gaussian/gaussian_user.cpp \
	applications/NCAR_transport/gaussian/gaussian.cpp \
	applications/NCAR_transport/gaussian/fclaw2d_map_cubedsphere.c \
	applications/NCAR_transport/gaussian/fclaw2d_map_pillowsphere.c \
	applications/NCAR_transport/gaussian/user_4.6/sphere_compute_error_fort.f \
	applications/NCAR_transport/all/setprob_transport.f \
	applications/NCAR_transport/all/qinit_transport.f \
	applications/NCAR_transport/all/setaux_transport.f \
	applications/NCAR_transport/all/b4step2_transport.f \
	applications/NCAR_transport/all/slotted_disk.f \
	applications/NCAR_transport/all/gaussian.f \
	applications/NCAR_transport/all/cosine_bell.f \
	applications/NCAR_transport/all/velocity_wind.f \
	applications/NCAR_transport/all/rpn2noncons.f \
	applications/NCAR_transport/all/rpt2noncons.f \
	applications/NCAR_transport/all/exact_solution.f

applications_NCAR_transport_gaussian_gaussian_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_NCAR_transport_gaussian_gaussian_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
