## Makefile.am in applications/clawpack/advection/correlatedcb

bin_PROGRAMS += applications/clawpack/advection/2d/correlatedcb/correlatedcb

TESTS += applications/clawpack/advection/2d/correlatedcb/regressions.txt

applications_clawpack_advection_2d_correlatedcb_correlatedcb_SOURCES = \
	applications/clawpack/advection/2d/correlatedcb/correlatedcb.cpp \
	applications/clawpack/advection/2d/correlatedcb/correlatedcb_user.cpp \
	applications/clawpack/advection/2d/correlatedcb/correlatedcb_options.c \
	applications/clawpack/advection/2d/correlatedcb/cosine_bell.f \
	applications/clawpack/advection/2d/correlatedcb/user_4.6/qinit.f \
	applications/clawpack/advection/2d/correlatedcb/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/all/fclaw2d_map_cubedsphere.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_pillowsphere.c \
	applications/clawpack/advection/2d/all/advection_patch_setup_manifold.c \
	applications/clawpack/advection/2d/all/advection_b4step2_manifold.c \
	applications/clawpack/advection/2d/all/NCAR_tracer_utils.f \
	applications/clawpack/advection/2d/all/user46_b4step2_manifold.f \
	applications/clawpack/advection/2d/all/user5_b4step2_manifold.f90 \
	applications/clawpack/advection/2d/all/advection_user_fort.h \
	applications/clawpack/advection/2d/all/user46_setaux_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90

EXTRA_DIST += \
	applications/clawpack/advection/2d/correlatedcb/fclaw_options.ini \
	applications/clawpack/advection/2d/correlatedcb/regression.ini \
	applications/clawpack/advection/2d/correlatedcb/afterframe.m \
	applications/clawpack/advection/2d/correlatedcb/setplot2.m \
	applications/clawpack/advection/2d/correlatedcb/mapc2m.m \
	applications/clawpack/advection/2d/correlatedcb/plot_mixing.m

applications_clawpack_advection_2d_correlatedcb_correlatedcb_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_advection_2d_correlatedcb_correlatedcb_LDADD = \
	    $(LDADD) \
	    $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD)  \
        $(FCLAW_LDADD)
