## Makefile.am in applications/clawpack/advection/2d/hemisphere

bin_PROGRAMS += applications/clawpack/advection/2d/hemisphere/hemisphere


applications_clawpack_advection_2d_hemisphere_hemisphere_SOURCES = \
	applications/clawpack/advection/2d/hemisphere/hemisphere.cpp \
	applications/clawpack/advection/2d/hemisphere/hemisphere_options.c \
	applications/clawpack/advection/2d/hemisphere/hemisphere_user.cpp \
	applications/clawpack/advection/2d/hemisphere/hemisphere_user.h \
	applications/clawpack/advection/2d/hemisphere/fclaw2d_map_pillowsphere.c \
	applications/clawpack/advection/2d/hemisphere/fclaw2d_map_pillowsphere5.c \
	applications/clawpack/advection/2d/hemisphere/setprob.f \
	applications/clawpack/advection/2d/hemisphere/fdisc.f \
	applications/clawpack/advection/2d/hemisphere/psi.f \
	applications/clawpack/advection/2d/hemisphere/user_4.6/qinit.f \
	applications/clawpack/advection/2d/hemisphere/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/all/user46_setaux_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpatch46_tag4refinement.f \
	applications/clawpack/advection/2d/all/clawpatch46_tag4coarsening.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpatch5_tag4refinement.f \
	applications/clawpack/advection/2d/all/clawpatch5_tag4coarsening.f

EXTRA_DIST += \
	applications/clawpack/advection/2d/hemisphere/fclaw_options.ini \
	applications/clawpack/advection/2d/hemisphere/regression.ini \
	applications/clawpack/advection/2d/hemisphere/afterframe.m \
	applications/clawpack/advection/2d/hemisphere/setplot2.m \
	applications/clawpack/advection/2d/hemisphere/mapc2m.m

applications_clawpack_advection_2d_hemisphere_hemisphere_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_advection_2d_hemisphere_hemisphere_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)

TESTS += applications/clawpack/advection/2d/hemisphere/regressions.sh