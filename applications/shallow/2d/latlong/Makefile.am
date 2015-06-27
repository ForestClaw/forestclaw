## Makefile.am in applications/shallow/2d/latlong

bin_PROGRAMS += applications/shallow/2d/latlong/latlong

applications_shallow_2d_latlong_latlong_SOURCES = \
	applications/shallow/2d/latlong/latlong.cpp \
	applications/shallow/2d/latlong/latlong_user.cpp \
	applications/shallow/2d/latlong/fclaw2d_map_latlong.c \
	applications/shallow/2d/latlong/fclaw2d_map_cubedsphere.c \
	applications/shallow/2d/latlong/setprob.f \
	applications/shallow/2d/latlong/write_latlong_data.f \
	applications/shallow/2d/latlong/user_4.6/qinit.f \
	applications/shallow/2d/latlong/user_4.6/gaussian.f \
	applications/shallow/2d/latlong/user_4.6/setaux_sphere.f \
	applications/shallow/2d/latlong/user_4.6/b4step2.f \
	applications/shallow/2d/latlong/user_4.6/rpn2swsphere.f \
	applications/shallow/2d/latlong/user_4.6/rpt2swq.f

applications_shallow_2d_latlong_latlong_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_shallow_2d_latlong_latlong_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
