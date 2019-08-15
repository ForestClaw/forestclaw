## Makefile.am in applications/paper/square

bin_PROGRAMS += applications/paper/square/square

applications_paper_square_square_SOURCES = \
	applications/paper/square/square.cpp \
	applications/paper/square/square_user.cpp \
	applications/paper/square/square_user.h \
	applications/paper/square/square_options.c \
	applications/paper/square/fclaw2d_map_identity.c \
	applications/paper/square/fclaw2d_map_cart.c \
	applications/paper/square/fclaw2d_map_fivepatch.c \
	applications/paper/square/fclaw2d_map_bilinear.c \
	applications/paper/square/setprob.f \
	applications/paper/square/q0.f \
	applications/paper/square/fdisc.f \
	applications/paper/square/square_basis.f \
	applications/paper/square/mapping_terms.f \
	applications/paper/square/square_velocity.f \
	applications/paper/square/qexact.f \
	applications/paper/square/dopri5.f \
	applications/paper/square/user_4.6/square_output.f \
	applications/paper/square/user_4.6/compute_error.f \
	applications/paper/square/user_4.6/qinit.f \
	applications/paper/square/user_4.6/square_setaux.f \
	applications/paper/square/user_4.6/rpn2cons_fw_manifold.f \
	applications/paper/square/user_4.6/rpt2cons_manifold.f \
	applications/paper/square/user_4.6/rpn2cons_update_manifold.f \
	applications/paper/square/user_4.6/tag4coarsening.f \
	applications/paper/square/user_4.6/tag4refinement.f

EXTRA_DIST += \
	applications/paper/square/fclaw_options.ini \
	applications/paper/square/afterframe.m \
	applications/paper/square/setplot2.m

applications_paper_square_square_CPPFLAGS = \
	    -I@top_srcdir@/applications/clawpack/advection/2d/all  \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_paper_square_square_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)

## @top_builddir@/src/solvers/fclaw2d_clawpack46/libclawpack5.la
