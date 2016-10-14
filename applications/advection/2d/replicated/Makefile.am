bin_PROGRAMS += applications/advection/2d/replicated/torus

applications_advection_2d_replicated_torus_SOURCES = \
	applications/advection/2d/replicated/torus.cpp \
	applications/advection/2d/replicated/torus_user.h \
	applications/advection/2d/replicated/torus_user.cpp \
	applications/advection/2d/replicated/fclaw2d_map_nomap.c \
	applications/advection/2d/replicated/setprob.f \
	applications/advection/2d/replicated/psi.f \
	applications/advection/2d/replicated/fdisc.f \
	applications/advection/2d/replicated/user_4.6/qinit.f \
	applications/advection/2d/replicated/user_5.0/qinit.f90 \
	applications/advection/2d/all/clawpack46_rpn2adv.f \
	applications/advection/2d/all/clawpack46_rpt2adv.f \
	applications/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/advection/2d/all/clawpack46_setaux.f \
	applications/advection/2d/all/clawpack5_setaux.f90

applications_advection_2d_replicated_torus_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_advection_2d_replicated_torus_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
