bin_PROGRAMS += applications/advection/2d/latlong/latlong


applications_advection_2d_latlong_latlong_SOURCES = \
	applications/advection/2d/latlong/latlong.cpp \
	applications/advection/2d/latlong/latlong_user.h \
	applications/advection/2d/latlong/latlong_user.cpp \
	applications/advection/2d/latlong/fclaw2d_map_latlong.c \
	applications/advection/2d/latlong/fdisc.f \
	applications/advection/2d/latlong/setprob.f \
	applications/advection/2d/latlong/psi.f \
	applications/advection/2d/latlong/user_4.6/qinit.f \
	applications/advection/2d/latlong/user_5.0/qinit.f90 \
	applications/advection/2d/all/user46_setaux_manifold.f \
	applications/advection/2d/all/user5_setaux_manifold.f90 \
	applications/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/advection/2d/all/clawpack5_rpt2adv_manifold.f90

applications_advection_2d_latlong_latlong_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_advection_2d_latlong_latlong_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
