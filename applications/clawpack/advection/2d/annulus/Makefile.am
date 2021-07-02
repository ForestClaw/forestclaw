bin_PROGRAMS += applications/clawpack/advection/2d/annulus/annulus


applications_clawpack_advection_2d_annulus_annulus_SOURCES = \
	applications/clawpack/advection/2d/annulus/annulus.cpp \
	applications/clawpack/advection/2d/annulus/annulus_options.c \
	applications/clawpack/advection/2d/annulus/annulus_user.h \
	applications/clawpack/advection/2d/annulus/annulus_user.cpp \
	applications/clawpack/advection/2d/annulus/fclaw2d_map_annulus.c \
	applications/clawpack/advection/2d/annulus/setprob.f \
	applications/clawpack/advection/2d/annulus/fdisc.f \
	applications/clawpack/advection/2d/annulus/psi.f \
	applications/clawpack/advection/2d/annulus/user_4.6/qinit.f \
	applications/clawpack/advection/2d/annulus/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/all/clawpatch46_tag4refinement.f \
	applications/clawpack/advection/2d/all/clawpatch46_tag4coarsening.f \
	applications/clawpack/advection/2d/all/clawpatch5_tag4refinement.f \
	applications/clawpack/advection/2d/all/clawpatch5_tag4coarsening.f \
	applications/clawpack/advection/2d/all/user46_setaux_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90

EXTRA_DIST += \
	applications/clawpack/advection/2d/annulus/fclaw_options.ini \
	applications/clawpack/advection/2d/annulus/regression.ini \
	applications/clawpack/advection/2d/annulus/afterframe.m \
	applications/clawpack/advection/2d/annulus/setplot2.m \
	applications/clawpack/advection/2d/annulus/mapc2m.m

applications_clawpack_advection_2d_annulus_annulus_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_advection_2d_annulus_annulus_LDADD = \
		$(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
	    $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)

TESTS += applications/clawpack/advection/2d/annulus/regressions.sh