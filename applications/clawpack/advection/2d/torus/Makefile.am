bin_PROGRAMS += applications/clawpack/advection/2d/torus/torus

applications_clawpack_advection_2d_torus_torus_SOURCES = \
	applications/clawpack/advection/2d/torus/torus.cpp \
	applications/clawpack/advection/2d/torus/torus_options.c \
	applications/clawpack/advection/2d/torus/torus_user.h \
	applications/clawpack/advection/2d/torus/torus_user.cpp \
	applications/clawpack/advection/2d/torus/fclaw2d_map_torus.c \
	applications/clawpack/advection/2d/torus/fclaw2d_map_cart.c \
	applications/clawpack/advection/2d/torus/setprob.f \
	applications/clawpack/advection/2d/torus/velocity.f \
	applications/clawpack/advection/2d/torus/fdisc.f \
	applications/clawpack/advection/2d/torus/qexact.f \
	applications/clawpack/advection/2d/torus/user_4.6/qinit.f \
	applications/clawpack/advection/2d/torus/user_4.6/compute_error.f \
	applications/clawpack/advection/2d/torus/user_4.6/torus46_setaux.f \
	applications/clawpack/advection/2d/torus/user_4.6/torus46_output.f \
	applications/clawpack/advection/2d/torus/user_4.6/tag4refinement.f \
	applications/clawpack/advection/2d/torus/user_4.6/tag4coarsening.f \
	applications/clawpack/advection/2d/torus/user_4.6/rpn2cons_fw_manifold.f \
    applications/clawpack/advection/2d/torus/user_4.6/rpt2cons_manifold.f \
	applications/clawpack/advection/2d/torus/user_4.6/rpn2cons_update_manifold.f \
	applications/clawpack/advection/2d/torus/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/torus/user_5.0/torus5_output.f90 \
	applications/clawpack/advection/2d/torus/user_5.0/setaux_manifold.f90 \
	applications/clawpack/advection/2d/torus/user_5.0/tag4refinement.f90 \
	applications/clawpack/advection/2d/torus/user_5.0/tag4coarsening.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90




EXTRA_DIST += \
	applications/clawpack/advection/2d/torus/fclaw_options.ini \
	applications/clawpack/advection/2d/torus/regression.ini \
	applications/clawpack/advection/2d/torus/afterframe.m \
	applications/clawpack/advection/2d/torus/setplot2.m \
	applications/clawpack/advection/2d/torus/mapc2m.m \
	applications/clawpack/advection/2d/torus/torus_soln.m \
	applications/clawpack/advection/2d/torus/torus_iso_init.m \
	applications/clawpack/advection/2d/torus/hsmooth.m 

applications_clawpack_advection_2d_torus_torus_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_advection_2d_torus_torus_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
