## Makefile.am in applications/paper/sphere

bin_PROGRAMS += applications/paper/sphere/sphere

applications_paper_sphere_sphere_SOURCES = \
	applications/paper/sphere/sphere.cpp \
	applications/paper/sphere/sphere_user.cpp \
	applications/paper/sphere/sphere_options.c \
	applications/paper/sphere/fclaw2d_map_cubedsphere.c \
	applications/paper/sphere/fclaw2d_map_pillowsphere.c \
	applications/paper/sphere/setprob.f \
	applications/paper/sphere/q0.f \
	applications/paper/sphere/fdisc.f \
	applications/paper/sphere/sphere_basis.f \
	applications/paper/sphere/mapping_terms.f \
	applications/paper/sphere/sphere_velocity.f \
	applications/paper/sphere/qexact.f \
	applications/paper/sphere/dopri5.f \
	applications/paper/sphere/user_4.6/sphere_output.f \
	applications/paper/sphere/user_4.6/compute_error.f \
	applications/paper/sphere/user_4.6/qinit.f \
	applications/paper/sphere/user_4.6/sphere_setaux.f \
	applications/paper/sphere/user_4.6/rpn2cons_fw_manifold.f \
	applications/paper/sphere/user_4.6/rpt2cons_manifold.f \
	applications/paper/sphere/user_4.6/rpn2_cons_update_manifold.f \
	applications/paper/sphere/user_4.6/tag4coarsening.f \
	applications/paper/sphere/user_4.6/tag4refinement.f

EXTRA_DIST += \
	applications/paper/sphere/fclaw_options.ini \
	applications/paper/sphere/afterframe.m \
	applications/paper/sphere/setplot2.m

applications_paper_sphere_sphere_CPPFLAGS = \
	    -I@top_srcdir@/applications/clawpack/advection/2d/all  \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_paper_sphere_sphere_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)

## @top_builddir@/src/solvers/fclaw2d_clawpack46/libclawpack5.la
