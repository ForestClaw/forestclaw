bin_PROGRAMS += applications/advection/2d/sphere/sphere


applications_advection_2d_sphere_sphere_SOURCES = \
	applications/advection/2d/sphere/sphere_user.cpp \
	applications/advection/2d/sphere/sphere_user.h \
	applications/advection/2d/sphere/sphere.cpp \
	applications/advection/2d/sphere/fclaw2d_map_cubedsphere.c \
	applications/advection/2d/sphere/fclaw2d_map_pillowsphere.c \
	applications/advection/2d/sphere/psi.f \
	applications/advection/2d/sphere/setprob.f \
	applications/advection/2d/sphere/fdisc.f \
	applications/advection/2d/sphere/user_4.6/qinit.f \
	applications/advection/2d/sphere/user_5.0/qinit.f90 \
	applications/advection/2d/sphere/user_4.6/tag4refinement.f \
	applications/advection/2d/sphere/user_4.6/tag4coarsening.f \
	applications/advection/2d/sphere/user_5.0/tag4refinement.f \
	applications/advection/2d/sphere/user_5.0/tag4coarsening.f \
	applications/advection/2d/all/user46_setaux_manifold.f \
	applications/advection/2d/all/user5_setaux_manifold.f90 \
	applications/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/advection/2d/all/clawpack5_rpt2adv_manifold.f90

applications_advection_2d_sphere_sphere_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_advection_2d_sphere_sphere_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
