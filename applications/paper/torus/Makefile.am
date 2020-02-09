bin_PROGRAMS += applications/paper/torus/torus

applications_paper_torus_torus_SOURCES = \
	applications/paper/torus/torus.cpp \
	applications/paper/torus/torus_options.c \
	applications/paper/torus/torus_user.h \
	applications/paper/torus/torus_user.cpp \
	applications/paper/torus/setprob.f \
	applications/paper/torus/compute_exact.f90 \
	applications/paper/torus/fclaw2d_map_torus.c \
	applications/paper/all/mapc2m_torus.f \
	applications/paper/all/torus_velocity.f \
	applications/paper/all/torus_basis.f \
	applications/paper/all/mapping_terms.f \
	applications/paper/torus/dopri5.f \
	applications/paper/torus/qexact.f \
	applications/paper/torus/q0.f \
	applications/paper/torus/fdisc.f \
	applications/paper/torus/user_4.6/qinit.f \
	applications/paper/torus/user_4.6/compute_error.f \
	applications/paper/torus/user_4.6/torus_setaux.f \
	applications/paper/torus/user_4.6/torus46_output.f \
	applications/paper/torus/user_4.6/torus_tag4coarsening.f \
	applications/paper/torus/user_4.6/torus_tag4refinement.f \
	applications/paper/torus/user_4.6/rpn2cons_fw_manifold.f \
    applications/paper/torus/user_4.6/rpt2cons_manifold.f \
	applications/paper/torus/user_4.6/rpn2cons_update_manifold.f




EXTRA_DIST += \
	applications/paper/torus/fclaw_options.ini \
	applications/paper/torus/afterframe.m \
	applications/paper/torus/setplot2.m \
	applications/paper/torus/mapc2m.m \
	applications/paper/torus/torus_soln.m \
	applications/paper/torus/torus_iso_init.m \
	applications/paper/torus/hsmooth.m 

applications_paper_torus_torus_CPPFLAGS = \
	    -I@top_srcdir@/applications/clawpack/advection/2d/all  \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_paper_torus_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
