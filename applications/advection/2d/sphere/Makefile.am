## Makefile.am in applications/advection/2d/cubed_sphere

bin_PROGRAMS += applications/advection/2d/cubed_sphere/cubed_sphere

applications_advection_2d_cubed_sphere_cubed_sphere_SOURCES = \
	applications/advection/2d/cubed_sphere/cubed_sphere_user.cpp \
	applications/advection/2d/cubed_sphere/cubed_sphere.cpp \
	applications/advection/2d/cubed_sphere/fclaw2d_map_cubedsphere.c \
	applications/advection/2d/cubed_sphere/fclaw2d_map_pillowsphere.c \
	applications/advection/2d/cubed_sphere/velocity_wind.f \
	applications/advection/2d/cubed_sphere/setprob.f \
	applications/advection/2d/cubed_sphere/user_4.6/tag4refinement.f \
	applications/advection/2d/cubed_sphere/user_4.6/rpn2noncons.f \
	applications/advection/2d/cubed_sphere/user_4.6/rpt2noncons.f \
	applications/advection/2d/cubed_sphere/user_4.6/qinit_manifold.f \
	applications/advection/2d/cubed_sphere/user_4.6/b4step2_manifold.f \
	applications/advection/2d/cubed_sphere/user_4.6/setaux_manifold.f \
	applications/advection/2d/cubed_sphere/user_4.6/sphere_output.f

applications_advection_2d_cubed_sphere_cubed_sphere_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_advection_2d_cubed_sphere_cubed_sphere = -D TRAPFPE

# applications_advection_2d_cubed_sphere_cubed_sphere_FFLAGS = \
#	$(AM_FFLAGS)

applications_advection_2d_cubed_sphere_cubed_sphere_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
