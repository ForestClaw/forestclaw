## Makefile.am in applications/advection/2d/sphere

bin_PROGRAMS += applications/advection/2d/sphere/sphere

applications_advection_2d_sphere_sphere_SOURCES = \
	applications/advection/2d/sphere/sphere_user.cpp \
	applications/advection/2d/sphere/sphere.cpp \
	applications/advection/2d/sphere/ghost_twopatch.f \
	applications/advection/2d/sphere/tag4refinement.f \
	applications/advection/2d/sphere/mapc2m.f \
	applications/advection/2d/sphere/velocity_wind.f \
	applications/advection/2d/sphere/rpn2noncons.f \
	applications/advection/2d/sphere/rpt2noncons.f \
	applications/advection/2d/sphere/qinit_manifold.f \
	applications/advection/2d/sphere/b4step2_manifold.f \
	applications/advection/2d/sphere/setaux_manifold.f \
	applications/advection/2d/sphere/sphere_output.f \
	applications/advection/2d/sphere/setprob.f

applications_advection_2d_sphere_sphere_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

# applications_advection_2d_sphere_sphere_FFLAGS = \
#	$(AM_FFLAGS)

applications_advection_2d_sphere_sphere_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
