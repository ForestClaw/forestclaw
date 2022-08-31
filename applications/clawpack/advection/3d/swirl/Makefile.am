## Makefile.am in applications/clawpack/advection/3d/swirl

bin_PROGRAMS += applications/clawpack/advection/3d/swirl/swirl

TESTS += applications/clawpack/advection/3d/swirl/regressions.sh

applications_clawpack_advection_3d_swirl_swirl_SOURCES = \
	applications/clawpack/advection/3d/swirl/swirl.cpp \
	applications/clawpack/advection/3d/swirl/swirl_user.cpp \
	applications/clawpack/advection/3d/swirl/swirl_user.h \
	applications/clawpack/advection/3d/swirl/swirl_map_extrude.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_cart.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_fivepatch.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_bilinear.c \
	applications/clawpack/advection/3d/swirl/swirl_options.c \
	applications/clawpack/advection/3d/swirl/setprob.f90 \
	applications/clawpack/advection/3d/swirl/user_exceeds_th.f90 \
	applications/clawpack/advection/3d/all/rpn3_vc_advection.f90 \
	applications/clawpack/advection/3d/all/rpt3_vc_advection.f90 \
	applications/clawpack/advection/3d/all/rptt3_vc_advection.f90 \
	applications/clawpack/advection/3d/swirl/user_4.6/qinit.f90 \
	applications/clawpack/advection/3d/swirl/user_4.6/swirl_setaux.f90 \
	applications/clawpack/advection/3d/swirl/user_4.6/swirl_setaux_manifold.f90

EXTRA_DIST += \
	applications/clawpack/advection/3d/swirl/fclaw_options.ini \
	applications/clawpack/advection/3d/swirl/mappings.ini \
	applications/clawpack/advection/3d/swirl/regression.ini \
	applications/clawpack/advection/3d/swirl/regression_map.ini \
	applications/clawpack/advection/3d/swirl/afterframe.m \
	applications/clawpack/advection/3d/swirl/read_vars.m \
	applications/clawpack/advection/3d/swirl/mapc2p.m \
	applications/clawpack/advection/3d/swirl/setplot3.m

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness
applications_clawpack_advection_3d_swirl_swirl_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK3_46_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


applications_clawpack_advection_3d_swirl_swirl_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK3_46_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
