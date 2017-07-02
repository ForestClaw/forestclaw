## Makefile.am in applications/clawpack/advection/2d/swirl_cons

bin_PROGRAMS += applications/clawpack/advection/2d/swirl_cons/swirlcons

applications_clawpack_advection_2d_swirl_cons_swirlcons_SOURCES = \
	applications/clawpack/advection/2d/swirl_cons/swirlcons_user.cpp \
	applications/clawpack/advection/2d/swirl_cons/swirlcons_user.h \
	applications/clawpack/advection/2d/swirl_cons/swirlcons.cpp \
	applications/clawpack/advection/2d/swirl_cons/psi.f \
	applications/clawpack/advection/2d/swirl_cons/setprob.f \
	applications/clawpack/advection/2d/all/clawpack46_setaux.f \
	applications/clawpack/advection/2d/all/clawpack_user.h \
	applications/clawpack/advection/2d/all/clawpack5_setaux.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/qinit.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/b4step2.f \
	applications/clawpack/advection/2d/swirl_cons/user_4.6/swirlcons_b4step2.f \
	applications/clawpack/advection/2d/swirl_cons/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/swirl_cons/user_5.0/b4step2.f90 \
	applications/clawpack/advection/2d/swirl_cons/Riemann_solvers/rpn2cons_cc.f \
	applications/clawpack/advection/2d/swirl_cons/Riemann_solvers/rpt2cons_cc.f \
	applications/clawpack/advection/2d/swirl_cons/swirlcons_bc2.f

applications_clawpack_advection_2d_swirl_cons_swirlcons_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_clawpack_advection_2d_swirl_cons_swirlcons_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fclaw2d_clawpack46/libclawpack5.la
