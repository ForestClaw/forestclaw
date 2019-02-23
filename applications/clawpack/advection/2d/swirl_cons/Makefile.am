## Makefile.am in applications/clawpack/advection/2d/swirl_cons

bin_PROGRAMS += applications/clawpack/advection/2d/swirl_cons/swirlcons

applications_clawpack_advection_2d_swirl_cons_swirlcons_SOURCES = \
	applications/clawpack/advection/2d/swirl_cons/swirlcons.cpp \
	applications/clawpack/advection/2d/swirl_cons/swirlcons_user.cpp \
	applications/clawpack/advection/2d/swirl_cons/swirlcons_user.h \
	applications/clawpack/advection/2d/swirl_cons/swirlcons_options.c \
	applications/clawpack/advection/2d/swirl_cons/fclaw2d_map_cart.c \
	applications/clawpack/advection/2d/swirl_cons/fclaw2d_map_fivepatch.c \
	applications/clawpack/advection/2d/swirl_cons/fclaw2d_map_bilinear.c \
	applications/clawpack/advection/2d/swirl_cons/setprob.f \
	applications/clawpack/advection/2d/swirl_cons/fdisc.f \
	applications/clawpack/advection/2d/swirl_cons/psi.f \
	applications/clawpack/advection/2d/swirl_cons/velocity.f \
	applications/clawpack/advection/2d/swirl_cons/qexact.f \
	applications/clawpack/advection/2d/swirl_cons/dopri5.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/compute_error.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/qinit.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/swirl46_setaux.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/rpn2cons_wd.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/rpn2cons_fw.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/rpn2cons_ec.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/rpt2cons.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/rpn2cons_fw_manifold.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/rpn2cons_ec_manifold.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/rpt2cons_manifold.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/rpn2cons_update_manifold.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/tag4coarsening.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/tag4refinement.f \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f 

EXTRA_DIST += \
	applications/clawpack/advection/2d/swirl_cons/fclaw_options.ini \
	applications/clawpack/advection/2d/swirl_cons/afterframe.m \
	applications/clawpack/advection/2d/swirl_cons/setplot2.m

applications_clawpack_advection_2d_swirl_cons_swirlcons_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_advection_2d_swirl_cons_swirlcons_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)

## @top_builddir@/src/solvers/fclaw2d_clawpack46/libclawpack5.la
