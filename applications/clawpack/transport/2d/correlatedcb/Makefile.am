## Makefile.am in applications/clawpack/transport/correlatedcb

bin_PROGRAMS += applications/clawpack/transport/2d/correlatedcb/correlatedcb

TESTS += applications/clawpack/transport/2d/correlatedcb/regressions.sh

applications_clawpack_transport_2d_correlatedcb_correlatedcb_SOURCES = \
	applications/clawpack/transport/2d/correlatedcb/correlatedcb.cpp \
	applications/clawpack/transport/2d/correlatedcb/cosine_bell.f \
	applications/clawpack/transport/2d/correlatedcb/user_4.6/qinit.f \
	applications/clawpack/transport/2d/correlatedcb/user_5.0/qinit.f90 \
	applications/clawpack/transport/2d/all/fclaw2d_map_cubedsphere.c \
	applications/clawpack/transport/2d/all/fclaw2d_map_pillowsphere.c \
	applications/clawpack/transport/2d/all/transport_setprob.f \
	applications/clawpack/transport/2d/all/velocity_wind.f \
	applications/clawpack/transport/2d/all/transport_util.f \
	applications/clawpack/transport/2d/all/transport_user.cpp \
	applications/clawpack/transport/2d/all/transport_domain.c \
	applications/clawpack/transport/2d/all/transport_options.c \
	applications/clawpack/transport/2d/all/user46_b4step2_manifold.f \
	applications/clawpack/transport/2d/all/user5_b4step2_manifold.f90 \
	applications/clawpack/advection/2d/all/advection_user_fort.h \
	applications/clawpack/advection/2d/all/user46_setaux_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90

EXTRA_DIST += \
	applications/clawpack/transport/2d/correlatedcb/fclaw_options.ini \
	applications/clawpack/transport/2d/correlatedcb/regression.ini \
	applications/clawpack/transport/2d/correlatedcb/afterframe.m \
	applications/clawpack/transport/2d/correlatedcb/setplot2.m \
	applications/clawpack/transport/2d/correlatedcb/mapc2m.m \
	applications/clawpack/transport/2d/correlatedcb/plot_mixing.m

applications_clawpack_transport_2d_correlatedcb_correlatedcb_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_transport_2d_correlatedcb_correlatedcb_LDADD = \
	    $(LDADD) \
	    $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD)  \
        $(FCLAW_LDADD)
