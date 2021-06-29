## Makefile.am in applications/clawpack/transport/2d/sphere

bin_PROGRAMS += applications/clawpack/transport/2d/sphere/sphere

applications_clawpack_transport_2d_sphere_sphere_SOURCES = \
	applications/clawpack/transport/2d/sphere/sphere.cpp \
	applications/clawpack/transport/2d/sphere/sphere_user.cpp \
	applications/clawpack/transport/2d/sphere/sphere_user.h \
	applications/clawpack/transport/2d/sphere/sphere_options.c \
	applications/clawpack/transport/2d/sphere/fclaw2d_map_cubedsphere.c \
	applications/clawpack/transport/2d/sphere/fclaw2d_map_pillowsphere.c \
	applications/clawpack/transport/2d/sphere/setprob.f \
	applications/clawpack/transport/2d/sphere/q0.f \
	applications/clawpack/transport/2d/sphere/fdisc.f \
	applications/clawpack/transport/2d/sphere/sphere_basis.f \
	applications/clawpack/transport/2d/sphere/mapping_terms.f \
	applications/clawpack/transport/2d/sphere/sphere_velocity.f \
	applications/clawpack/transport/2d/sphere/qexact.f \
	applications/clawpack/transport/2d/sphere/dopri5.f \
	applications/clawpack/transport/2d/sphere/user_4.6/sphere_output.f \
	applications/clawpack/transport/2d/sphere/user_4.6/compute_error.f \
	applications/clawpack/transport/2d/sphere/user_4.6/qinit.f \
	applications/clawpack/transport/2d/sphere/user_4.6/sphere_setaux.f \
	applications/clawpack/transport/2d/sphere/user_4.6/rpn2cons_fw_manifold.f \
	applications/clawpack/transport/2d/sphere/user_4.6/rpt2cons_manifold.f \
	applications/clawpack/transport/2d/sphere/user_4.6/rpn2_cons_update_manifold.f \
	applications/clawpack/transport/2d/sphere/user_4.6/tag4coarsening.f \
	applications/clawpack/transport/2d/sphere/user_4.6/tag4refinement.f

EXTRA_DIST += \
	applications/clawpack/transport/2d/sphere/fclaw_options.ini \
	applications/clawpack/transport/2d/sphere/afterframe.m \
	applications/clawpack/transport/2d/sphere/setplot2.m

applications_clawpack_transport_2d_sphere_sphere_CPPFLAGS = \
	    -I@top_srcdir@/applications/clawpack/advection/2d/all  \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_transport_2d_sphere_sphere_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)

## @top_builddir@/src/solvers/fclaw2d_clawpack46/libclawpack5.la
