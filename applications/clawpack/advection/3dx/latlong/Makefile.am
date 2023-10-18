bin_PROGRAMS += applications/clawpack/advection/3dx/latlong/latlong

TESTS += applications/clawpack/advection/3dx/latlong/regressions.sh

applications_clawpack_advection_3d_latlong_latlong_SOURCES = \
	applications/clawpack/advection/3dx/latlong/latlong.cpp \
	applications/clawpack/advection/3dx/latlong/latlong_options.c \
	applications/clawpack/advection/3dx/latlong/latlong_user.h \
	applications/clawpack/advection/3dx/latlong/latlong_user.cpp \
	applications/clawpack/advection/3dx/latlong/latlong_map_extrude.c \
	applications/clawpack/advection/3dx/latlong/fdisc.f \
	applications/clawpack/advection/3dx/latlong/setprob.f90 \
	applications/clawpack/advection/3dx/latlong/psi.f90 \
	applications/clawpack/advection/3dx/latlong/user_4.6/qinit.f90 \
	applications/clawpack/advection/3dx/all/advection_patch_setup_manifold.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_latlong.c \
	applications/clawpack/advection/3dx/all/claw3_setaux_manifold.f90 \
	applications/clawpack/advection/3dx/all/rpn3_vc_advection.f90 \
	applications/clawpack/advection/3dx/all/rpt3_vc_advection.f90 \
	applications/clawpack/advection/3dx/all/rptt3_vc_advection.f90

EXTRA_DIST += \
	applications/clawpack/advection/3dx/latlong/fclaw_options.ini \
	applications/clawpack/advection/3dx/latlong/regression.ini \
	applications/clawpack/advection/3dx/latlong/afterframe.m \
	applications/clawpack/advection/3dx/latlong/setplot2.m \
	applications/clawpack/advection/3dx/latlong/mapc2m.m

applications_clawpack_advection_3d_latlong_latlong_CPPFLAGS = \
	$(AM_CPPFLAGS) \
	$(FCLAW_CLAWPACK46_CPPFLAGS) \
	$(FCLAW_CLAWPATCH_CPPFLAGS)


applications_clawpack_advection_3d_latlong_latlong_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD)  \
        $(FCLAW_LDADD)
