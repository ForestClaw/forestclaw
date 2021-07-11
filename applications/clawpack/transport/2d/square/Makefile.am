## Makefile.am in applications/clawpack/transport/2d/square

bin_PROGRAMS += applications/clawpack/transport/2d/square/square

applications_clawpack_transport_2d_square_square_SOURCES = \
	applications/clawpack/transport/2d/square/square.cpp \
	applications/clawpack/transport/2d/square/square_user.cpp \
	applications/clawpack/transport/2d/square/square_user.h \
	applications/clawpack/transport/2d/square/square_options.c \
	applications/clawpack/transport/2d/all/fclaw2d_map_identity.c \
	applications/clawpack/transport/2d/all/fclaw2d_map_cart.c \
	applications/clawpack/transport/2d/all/fclaw2d_map_fivepatch.c \
	applications/clawpack/transport/2d/all/fclaw2d_map_bilinear.c \
	applications/clawpack/transport/2d/square/setprob.f \
	applications/clawpack/transport/2d/square/q0.f \
	applications/clawpack/transport/2d/square/fdisc.f \
	applications/clawpack/transport/2d/square/square_basis.f \
	applications/clawpack/transport/2d/square/mapping_terms.f \
	applications/clawpack/transport/2d/square/square_velocity.f \
	applications/clawpack/transport/2d/square/rpn2cons_update_manifold.f \
	applications/clawpack/transport/2d/square/qexact.f \
	applications/clawpack/transport/2d/square/dopri5.f \
	applications/clawpack/transport/2d/square/user_4.6/qinit.f \
	applications/clawpack/transport/2d/square/user_5.0/qinit.f \
	applications/clawpack/transport/2d/all/transport_patch_setup_manifold.c \
	applications/clawpack/transport/2d/all/transport_b4step2_manifold.c \
	applications/clawpack/transport/2d/all/clawpatch46_setaux_manifold.f \
	applications/clawpack/transport/2d/all/transport46_set_velocity.f \
	applications/clawpack/transport/2d/square/user_4.6/square_output.f \
	applications/clawpack/transport/2d/square/user_4.6/compute_error.f \
	applications/clawpack/transport/2d/square/user_4.6/rpn2cons_fw_manifold.f \
	applications/clawpack/transport/2d/square/user_4.6/rpt2cons_manifold.f \
	applications/clawpack/transport/2d/square/user_5.0/square_output.f \
	applications/clawpack/transport/2d/square/user_5.0/compute_error.f \
	applications/clawpack/transport/2d/square/user_5.0/square_setaux.f \
	applications/clawpack/transport/2d/square/user_5.0/rpn2cons_fw_manifold.f \
	applications/clawpack/transport/2d/square/user_5.0/rpt2cons_manifold.f

EXTRA_DIST += \
	applications/clawpack/transport/2d/square/fclaw_options.ini \
	applications/clawpack/transport/2d/square/afterframe.m \
	applications/clawpack/transport/2d/square/setplot2.m

applications_clawpack_transport_2d_square_square_CPPFLAGS = \
	    -I@top_srcdir@/applications/clawpack/advection/2d/all  \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_transport_2d_square_square_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)

## @top_builddir@/src/solvers/fclaw2d_clawpack46/libclawpack5.la
