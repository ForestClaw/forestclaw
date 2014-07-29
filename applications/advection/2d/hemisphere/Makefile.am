## Makefile.am in applications/advection/2d/hemisphere



bin_PROGRAMS += applications/advection/2d/hemisphere/hemisphere



applications_advection_2d_hemisphere_hemisphere_SOURCES = \
	applications/advection/2d/hemisphere/hemisphere_user.cpp \
	applications/advection/2d/hemisphere/hemisphere.cpp \
	applications/advection/2d/hemisphere/ghost_twopatch.f \
	applications/advection/2d/hemisphere/tag4refinement.f \
	applications/advection/2d/hemisphere/mapc2m.f \
	applications/advection/2d/hemisphere/velocity_wind.f \
	applications/advection/2d/hemisphere/rpn2noncons.f \
	applications/advection/2d/hemisphere/rpt2noncons.f \
	applications/advection/2d/hemisphere/qinit_manifold.f \
	applications/advection/2d/hemisphere/b4step2_manifold.f \
	applications/advection/2d/hemisphere/setaux_manifold.f \
	applications/advection/2d/hemisphere/setprob.f

applications_advection_2d_hemisphere_hemisphere_CPPFLAGS = \
	$(FCLAW_CPPFLAGS) $(AM_CPPFLAGS)

applications_advection_2d_hemisphere_hemisphere_FFLAGS = $(FCLAW_FFLAGS)

applications_advection_2d_hemisphere_hemisphere_LDADD = \
	@top_builddir@/src/solvers/fclaw2d_clawpack4.6/libclawpack4.6.la \
	$(LDADD)
