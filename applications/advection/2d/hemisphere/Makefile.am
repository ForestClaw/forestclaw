## Makefile.am in applications/advection/2d/hemisphere

bin_PROGRAMS += applications/advection/2d/hemisphere/hemisphere

applications_advection_2d_hemisphere_hemisphere_SOURCES = \
	applications/advection/2d/hemisphere/hemisphere_user.cpp \
	applications/advection/2d/hemisphere/hemisphere.cpp \
	applications/advection/2d/hemisphere/velocity_wind.f \
	applications/advection/2d/hemisphere/fclaw2d_map_pillowsphere.f \
	applications/advection/2d/hemisphere/fclaw2d_map_pillowsphere5.f \
	applications/advection/2d/hemisphere/setprob.f \
	applications/advection/2d/hemisphere/user_4.6/tag4refinement.f \
	applications/advection/2d/hemisphere/user_4.6/hemisphere_output.f \
	applications/advection/2d/hemisphere/user_4.6/rpn2noncons.f \
	applications/advection/2d/hemisphere/user_4.6/rpt2noncons.f \
	applications/advection/2d/hemisphere/user_4.6/qinit_manifold.f \
	applications/advection/2d/hemisphere/user_4.6/b4step2_manifold.f \
	applications/advection/2d/hemisphere/user_4.6/setaux_manifold.f

applications_advection_2d_hemisphere_hemisphere_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_advection_2d_hemisphere_hemisphere = -D TRAPFPE

# applications_advection_2d_hemisphere_hemisphere_FFLAGS = \
#	$(AM_FFLAGS)

applications_advection_2d_hemisphere_hemisphere_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
