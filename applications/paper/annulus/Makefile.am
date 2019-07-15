bin_PROGRAMS += applications/paper/annulus/annulus


applications_paper_annulus_annulus_SOURCES = \
	applications/paper/annulus/annulus.cpp \
	applications/paper/annulus/annulus_options.c \
	applications/paper/annulus/annulus_user.h \
	applications/paper/annulus/annulus_user.cpp \
	applications/paper/annulus/fclaw2d_map_annulus.c \
	applications/paper/annulus/setprob.f \
	applications/paper/annulus/dopri5.f \
	applications/paper/annulus/qexact.f \
	applications/paper/annulus/fdisc.f \
	applications/paper/annulus/psi.f \
	applications/paper/annulus/annulus_psi.f \
	applications/paper/annulus/velocity.f \
	applications/paper/annulus/user_4.6/qinit.f \
	applications/paper/annulus/user_4.6/rpn2cons_fw_manifold.f \
	applications/paper/annulus/user_4.6/rpt2cons_manifold.f \
	applications/paper/annulus/user_4.6/rpn2cons_update_manifold.f \
	applications/paper/annulus/user_4.6/annulus46_setaux.f \
	applications/paper/annulus/user_4.6/annulus46_output.f \
	applications/paper/annulus/user_4.6/tag4refinement.f \
	applications/paper/annulus/user_4.6/tag4coarsening.f \
	applications/paper/annulus/user_4.6/compute_error.f \
	applications/paper/annulus/user_4.6/annulus46_rpt2adv_manifold.f \
	applications/paper/annulus/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90

EXTRA_DIST += \
	applications/paper/annulus/regression.ini \
	applications/paper/annulus/afterframe.m \
	applications/paper/annulus/setplot2.m \
	applications/paper/annulus/mapc2m.m

applications_paper_annulus_annulus_CPPFLAGS = \
	    -I@top_srcdir@/applications/clawpack/advection/2d/all  \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_paper_annulus_annulus_LDADD = \
		$(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
	    $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
