bin_PROGRAMS += applications/paper/cylinder/cylinder

applications_paper_cylinder_cylinder_SOURCES = \
	applications/paper/cylinder/cylinder.cpp \
	applications/paper/cylinder/cylinder_options.c \
	applications/paper/cylinder/cylinder_user.h \
	applications/paper/cylinder/cylinder_user.cpp \
	applications/paper/cylinder/setprob.f \
	applications/paper/cylinder/compute_exact.f90 \
	applications/paper/cylinder/fclaw2d_map_cylinder.c \
	applications/paper/all/mapc2m_cylinder.f \
	applications/paper/all/cylinder_velocity.f \
	applications/paper/all/cylinder_basis.f \
	applications/paper/all/mapping_terms.f \
	applications/paper/cylinder/dopri5.f \
	applications/paper/cylinder/qexact.f \
	applications/paper/cylinder/q0.f \
	applications/paper/cylinder/fdisc.f \
	applications/paper/cylinder/user_4.6/qinit.f \
	applications/paper/cylinder/user_4.6/compute_error.f \
	applications/paper/cylinder/user_4.6/cylinder_setaux.f \
	applications/paper/cylinder/user_4.6/cylinder46_output.f \
	applications/paper/cylinder/user_4.6/cylinder_tag4coarsening.f \
	applications/paper/cylinder/user_4.6/cylinder_tag4refinement.f \
	applications/paper/cylinder/user_4.6/rpn2cons_fw_manifold.f \
    applications/paper/cylinder/user_4.6/rpt2cons_manifold.f \
	applications/paper/cylinder/user_4.6/rpn2cons_update_manifold.f




EXTRA_DIST += \
	applications/paper/cylinder/fclaw_options.ini \
	applications/paper/cylinder/afterframe.m \
	applications/paper/cylinder/setplot2.m \
	applications/paper/cylinder/mapc2m.m \
	applications/paper/cylinder/cylinder_soln.m \
	applications/paper/cylinder/cylinder_iso_init.m \
	applications/paper/cylinder/hsmooth.m 

applications_paper_cylinder_cylinder_CPPFLAGS = \
	    -I@top_srcdir@/applications/clawpack/advection/2d/all  \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_paper_cylinder_cylinder_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
