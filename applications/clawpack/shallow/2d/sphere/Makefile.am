## Makefile.am in applications/clawpack/shallow/2d/sphere

bin_PROGRAMS += applications/clawpack/shallow/2d/sphere/sphere


applications_clawpack_shallow_2d_sphere_sphere_SOURCES = \
	applications/clawpack/shallow/2d/sphere/sphere.cpp \
	applications/clawpack/shallow/2d/sphere/sphere_user.h \
	applications/clawpack/shallow/2d/sphere/sphere_user.cpp \
	applications/clawpack/shallow/2d/sphere/fclaw2d_map_latlong.c \
	applications/clawpack/shallow/2d/sphere/fclaw2d_map_cubedsphere.c \
	applications/clawpack/shallow/2d/sphere/setprob.f \
	applications/clawpack/shallow/2d/sphere/write_latlong_data.f \
	applications/clawpack/shallow/2d/sphere/user_4.6/qinit.f \
	applications/clawpack/shallow/2d/sphere/user_4.6/gaussian.f \
	applications/clawpack/shallow/2d/sphere/user_4.6/setaux_sphere.f \
	applications/clawpack/shallow/2d/sphere/user_4.6/b4step2.f \
	applications/clawpack/shallow/2d/sphere/user_4.6/rpn2swsphere.f \
	applications/clawpack/shallow/2d/sphere/user_4.6/rpt2swq.f

applications_clawpack_shallow_2d_sphere_sphere_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_clawpack_shallow_2d_sphere_sphere_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
