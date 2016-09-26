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
	applications/advection/2d/hemisphere/user_4.6/clawpack46_rpn2.f \
	applications/advection/2d/hemisphere/user_4.6/clawpack46_rpt2.f \
	applications/advection/2d/hemisphere/user_4.6/clawpack46_qinit.f \
	applications/advection/2d/hemisphere/user_4.6/hemisphere46_setaux_manifold.f

applications_advection_2d_hemisphere_hemisphere_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_advection_2d_hemisphere_hemisphere_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
