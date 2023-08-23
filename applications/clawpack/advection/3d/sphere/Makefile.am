bin_PROGRAMS += applications/clawpack/advection/3d/sphere/sphere

TESTS += applications/clawpack/advection/3d/sphere/regressions.txt

applications_clawpack_advection_3d_sphere_sphere_SOURCES = \
	applications/clawpack/advection/3d/sphere/sphere.cpp \
	applications/clawpack/advection/3d/sphere/sphere_user.cpp \
	applications/clawpack/advection/3d/sphere/sphere_user.h \
	applications/clawpack/advection/3d/sphere/sphere_options.c \
	applications/clawpack/advection/3d/sphere/sphere_map_extrude.c \
	applications/clawpack/advection/3d/sphere/setprob.f \
	applications/clawpack/advection/3d/sphere/fdisc.f \
	applications/clawpack/advection/3d/sphere/psi.f90 \
	applications/clawpack/advection/3d/sphere/user_4.6/qinit.f90 \
	applications/clawpack/advection/2d/all/fclaw2d_map_cubedsphere.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_pillowsphere.c \
	applications/clawpack/advection/3d/all/claw3_setaux_manifold.f90 \
	applications/clawpack/advection/3d/all/advection_patch_setup_manifold.c  \
	applications/clawpack/advection/3d/all/rpn3_vc_advection.f90 \
	applications/clawpack/advection/3d/all/rpt3_vc_advection.f90 \
	applications/clawpack/advection/3d/all/rptt3_vc_advection.f90

EXTRA_DIST += \
	applications/clawpack/advection/3d/sphere/fclaw_options.ini \
	applications/clawpack/advection/3d/sphere/regression.ini \
	applications/clawpack/advection/3d/sphere/afterframe.m \
	applications/clawpack/advection/3d/sphere/read_vars.m \
	applications/clawpack/advection/3d/sphere/parallelpartitions.m \
	applications/clawpack/advection/3d/sphere/setplot3.m \
	applications/clawpack/advection/3d/sphere/mapc2p.m

applications_clawpack_advection_3d_sphere_sphere_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_advection_3d_sphere_sphere_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
