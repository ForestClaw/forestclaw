## Makefile.am in applications/shallow/2d/sphere

bin_PROGRAMS += applications/shallow/2d/sphere/sphere

applications_shallow_2d_sphere_sphere_SOURCES = \
	applications/shallow/2d/sphere/sphere.cpp \
	applications/shallow/2d/sphere/sphere_user.cpp \
	applications/shallow/2d/sphere/fclaw2d_map_latlong.c \
	applications/shallow/2d/sphere/fclaw2d_map_cubedsphere.c \
	applications/shallow/2d/sphere/setprob.f \
	applications/shallow/2d/sphere/write_latlong_data.f \
	applications/shallow/2d/sphere/user_4.6/qinit.f \
	applications/shallow/2d/sphere/user_4.6/gaussian.f \
	applications/shallow/2d/sphere/user_4.6/setaux_sphere.f \
	applications/shallow/2d/sphere/user_4.6/b4step2.f \
	applications/shallow/2d/sphere/user_4.6/rpn2swsphere.f \
	applications/shallow/2d/sphere/user_4.6/rpt2swq.f

applications_shallow_2d_sphere_sphere_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_shallow_2d_sphere_sphere_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
