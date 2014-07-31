## Makefile.am in applications/advection/2d/cubed_sphere

bin_PROGRAMS += applications/advection/2d/cubed_sphere/cubed_sphere

applications_advection_2d_cubed_sphere_cubed_sphere_SOURCES = \
	applications/advection/2d/cubed_sphere/cubed_sphere_user.cpp \
	applications/advection/2d/cubed_sphere/cubed_sphere.cpp \
	applications/advection/2d/cubed_sphere/ghost_twopatch.f \
	applications/advection/2d/cubed_sphere/tag4refinement.f \
	applications/advection/2d/cubed_sphere/mapc2m.f \
	applications/advection/2d/cubed_sphere/velocity_wind.f \
	applications/advection/2d/cubed_sphere/rpn2noncons.f \
	applications/advection/2d/cubed_sphere/rpt2noncons.f \
	applications/advection/2d/cubed_sphere/qinit_manifold.f \
	applications/advection/2d/cubed_sphere/b4step2_manifold.f \
	applications/advection/2d/cubed_sphere/setaux_manifold.f \
	applications/advection/2d/cubed_sphere/sphere_output.f \
	applications/advection/2d/cubed_sphere/setprob.f

applications_advection_2d_cubed_sphere_cubed_sphere_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

# applications_advection_2d_cubed_sphere_cubed_sphere_FFLAGS = \
#	$(AM_FFLAGS)

applications_advection_2d_cubed_sphere_cubed_sphere_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
