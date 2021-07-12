bin_PROGRAMS += applications/clawpack/transport/2d/torus/torus

TESTS += applications/clawpack/transport/2d/torus/regressions.sh

applications_clawpack_transport_2d_torus_torus_SOURCES = \
	applications/clawpack/transport/2d/torus/torus.cpp \
	applications/clawpack/transport/2d/torus/torus_options.c \
	applications/clawpack/transport/2d/torus/torus_user.h \
	applications/clawpack/transport/2d/torus/torus_user.cpp \
	applications/clawpack/transport/2d/torus/setprob.f \
	applications/clawpack/transport/2d/torus/torus_velocity.f \
	applications/clawpack/transport/2d/torus/torus_psi.f \
	applications/clawpack/transport/2d/torus/fdisc.f \
	applications/clawpack/transport/2d/torus/fclaw2d_map_torus.c \
	applications/clawpack/transport/2d/torus/user_4.6/qinit.f \
	applications/clawpack/transport/2d/torus/user_5.0/qinit.f90 \
	applications/clawpack/transport/2d/all/transport_patch_setup_manifold.c \
	applications/clawpack/transport/2d/all/transport_b4step2_manifold.c \
	applications/clawpack/transport/2d/all/transport46_setaux_manifold.f \
	applications/clawpack/transport/2d/all/transport46_set_velocity.f \
	applications/clawpack/transport/2d/all/clawpack46_rpn2cons_fw_manifold.f \
    applications/clawpack/transport/2d/all/clawpack46_rpt2cons_manifold.f \
	applications/clawpack/transport/2d/all/transport5_setaux_manifold.f \
	applications/clawpack/transport/2d/all/transport5_set_velocity.f \
	applications/clawpack/transport/2d/all/clawpack5_rpn2cons_fw_manifold.f \
    applications/clawpack/transport/2d/all/clawpack5_rpt2cons_manifold.f \
	applications/clawpack/transport/2d/all/rpn2cons_update_manifold.f


EXTRA_DIST += \
	applications/clawpack/transport/2d/torus/fclaw_options.ini \
	applications/clawpack/transport/2d/torus/regression.ini \
	applications/clawpack/transport/2d/torus/afterframe.m \
	applications/clawpack/transport/2d/torus/setplot2.m \
	applications/clawpack/transport/2d/torus/mapc2m.m 

applications_clawpack_transport_2d_torus_torus_CPPFLAGS = \
	    -I@top_srcdir@/applications/clawpack/transport/2d/all  \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_transport_2d_torus_torus_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
