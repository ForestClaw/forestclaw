## Makefile.am in applications/advection/2d/hemisphere

bin_PROGRAMS += applications/advection/2d/hemisphere/hemisphere


applications_advection_2d_hemisphere_hemisphere_SOURCES = \
	applications/advection/2d/hemisphere/hemisphere_user.cpp \
	applications/advection/2d/hemisphere/hemisphere_user.h \
	applications/advection/2d/hemisphere/hemisphere.cpp \
	applications/advection/2d/hemisphere/velocity_wind.f \
	applications/advection/2d/hemisphere/fclaw2d_map_pillowsphere.c \
	applications/advection/2d/hemisphere/fclaw2d_map_pillowsphere5.c \
	applications/advection/2d/hemisphere/setprob.f \
	applications/advection/2d/hemisphere/fdisc.f \
	applications/advection/2d/hemisphere/user_4.6/rpn2.f \
	applications/advection/2d/hemisphere/user_4.6/rpt2.f \
	applications/advection/2d/hemisphere/user_4.6/qinit.f \
	applications/advection/2d/hemisphere/user_4.6/setaux_manifold.f \
	applications/advection/2d/hemisphere/user_5.0/rpn2.f90 \
	applications/advection/2d/hemisphere/user_5.0/rpt2.f90 \
	applications/advection/2d/hemisphere/user_5.0/qinit.f90 \
	applications/advection/2d/hemisphere/user_5.0/setaux_manifold.f90

applications_advection_2d_hemisphere_hemisphere_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_advection_2d_hemisphere_hemisphere_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
