## Makefile.am in applications/advection/2.5d/swirl

bin_PROGRAMS += applications/advection/2.5d/swirl/swirl

applications_advection_2_5d_swirl_swirl_SOURCES = \
	applications/advection/2.5d/swirl/swirl_user.cpp \
	applications/advection/2.5d/swirl/swirl_user.h \
	applications/advection/2.5d/swirl/swirl.cpp \
	applications/advection/2.5d/swirl/fclaw2d_map_nomap.c \
	applications/advection/2.5d/swirl/psi.f \
	applications/advection/2.5d/swirl/setprob.f \
	applications/advection/2.5d/swirl/user_5.0/qinit.f90

applications_advection_2_5d_swirl_swirl_CPPFLAGS = \
    $(AM_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS) $(FCLAW_CLAWPACK5_3D_CPPFLAGS) \
	$(FCLAW_DUMMY_CPPFLAGS)

applications_advection_2_5d_swirl_swirl_LDADD = \
    $(FCLAW_CLAWPACK5_LDADD) $(FCLAW_DUMMY_LDADD) $(FCLAW_CLAWPACK5_3D_LDADD) $(LDADD)
