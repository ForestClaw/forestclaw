## Makefile.am in applications/metric/2d/sphere

bin_PROGRAMS += applications/metric/2d/sphere/metric


applications_metric_2d_sphere_metric_SOURCES = \
	applications/metric/2d/sphere/metric.cpp \
	applications/metric/2d/sphere/metric_user.cpp \
	applications/metric/2d/sphere/metric_diagnostics.cpp \
	applications/metric/2d/sphere/tag4refinement.f \
	applications/metric/2d/sphere/metric_output.f \
	applications/metric/2d/sphere/initialize.f \
	applications/metric/2d/sphere/compute_error.f \
	applications/metric/2d/sphere/setprob.f \
	applications/metric/2d/sphere/diagnostics.f \
	applications/metric/2d/sphere/fclaw2d_map_nomap.c \
	applications/metric/2d/sphere/fclaw2d_map_cart.c \
	applications/metric/2d/sphere/fclaw2d_map_fivepatch.c \
	applications/metric/2d/sphere/fclaw2d_map_pillowdisk.c \
	applications/metric/2d/sphere/fclaw2d_map_pillowdisk5.c \
	applications/metric/2d/sphere/fclaw2d_map_squareddisk.c \
	applications/metric/2d/sphere/fclaw2d_map_pillowsphere.c \
	applications/metric/2d/sphere/fclaw2d_map_cubedsphere.c \
	applications/metric/2d/sphere/fclaw2d_map_torus.c

applications_advection_2d_metric_metric_CPPFLAGS = \
        $(AM_CPPFLAGS)


# applications_advection_2d_metric_metric_CPPFLAGS += -D TRAPFPE

# applications_advection_2d_metric_metric_CPPFLAGS += -D MPI_DEBUG


applications_sphere_metric_LDADD = $(LDADD) $(FLIBS)
