bin_PROGRAMS += applications/paper/torus_conservation/torus

applications_paper_torus_conservation_torus_SOURCES = \
	applications/paper/torus_conservation/torus.cpp \
	applications/paper/torus_conservation/torus_options.c \
	applications/paper/torus_conservation/torus_user.h \
	applications/paper/torus_conservation/torus_user.cpp \
	applications/paper/torus_conservation/fclaw2d_map_torus.c \
	applications/paper/torus_conservation/setprob.f \
	applications/paper/torus_conservation/torus_velocity.f \
v	applications/paper/torus_conservation/torus_psi.f \
	applications/paper/torus_conservation/fdisc.f \
	applications/paper/torus_conservation/user_4.6/qinit.f \
	applications/paper/torus_conservation/user_4.6/torus46_setaux.f \
	applications/paper/torus_conservation/user_4.6/tag4refinement.f \
	applications/paper/torus_conservation/user_4.6/tag4coarsening.f \
	applications/paper/torus_conservation/user_4.6/rpn2cons_fw_manifold.f \
    applications/paper/torus_conservation/user_4.6/rpt2cons_manifold.f \
	applications/paper/torus_conservation/user_4.6/rpn2cons_update_manifold.f




EXTRA_DIST += \
	applications/paper/torus_conservation/fclaw_options.ini \
	applications/paper/torus_conservation/afterframe.m \
	applications/paper/torus_conservation/setplot2.m \
	applications/paper/torus_conservation/mapc2m.m

applications_paper_torus_conservation_torus_CPPFLAGS = \
	    -I@top_srcdir@/applications/clawpack/advection/2d/all  \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_paper_torus_conservation_torus_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
