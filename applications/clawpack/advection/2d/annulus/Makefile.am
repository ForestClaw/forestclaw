bin_PROGRAMS += applications/clawpack/advection/2d/annulus/annulus


applications_clawpack_advection_2d_annulus_annulus_SOURCES = \
	applications/clawpack/advection/2d/annulus/annulus.cpp \
	applications/clawpack/advection/2d/annulus/annulus_user.h \
	applications/clawpack/advection/2d/annulus/annulus_user.cpp \
	applications/clawpack/advection/2d/annulus/fclaw2d_map_annulus.c \
	applications/clawpack/advection/2d/annulus/setprob.f \
	applications/clawpack/advection/2d/annulus/fdisc.f \
	applications/clawpack/advection/2d/annulus/psi.f \
	applications/clawpack/advection/2d/annulus/user_4.6/qinit.f \
	applications/clawpack/advection/2d/annulus/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/all/user46_setaux_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90

applications_clawpack_advection_2d_annulus_annulus_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_clawpack_advection_2d_annulus_annulus_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
