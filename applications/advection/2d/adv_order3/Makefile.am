## Makefile.am in applications/advection/2d/adv_order3

bin_PROGRAMS += applications/advection/2d/adv_order3/adv

applications_advection_2d_adv_order3_adv_SOURCES = \
	applications/advection/2d/adv_order3/adv_user.cpp \
	applications/advection/2d/adv_order3/adv_user.h \
	applications/advection/2d/adv_order3/adv.cpp \
	applications/advection/2d/adv_order3/fclaw2d_map_nomap.c \
	applications/advection/2d/adv_order3/psi.f \
	applications/advection/2d/adv_order3/setprob.f \
	applications/advection/2d/all/clawpack46_setaux.f \
	applications/advection/2d/all/clawpack_user.h \
	applications/advection/2d/all/clawpack46_rpn2adv.f \
	applications/advection/2d/all/clawpack46_rpt2adv.f \
	applications/advection/2d/adv_order3/user_4.6/qinit.f \
	applications/advection/2d/adv_order3/user_4.6/b4step2.f \
	applications/advection/2d/adv_order3/user_4.6/tag4refinement.f \
	applications/advection/2d/adv_order3/user_4.6/tag4coarsening.f

applications_advection_2d_adv_order3_adv_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS) \
	$(FCLAW_DUMMY_CPPFLAGS)

applications_advection_2d_adv_order3_adv_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_DUMMY_LDADD)  $(FCLAW_CLAWPACK5_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fclaw2d_clawpack46/libclawpack5.la
