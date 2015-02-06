## Makefile.am in applications/metric/2d/all_mappings

bin_PROGRAMS += applications/metric/2d/all_mappings/metric


applications_metric_2d_all_mappings_metric_SOURCES = \
	applications/metric/2d/all_mappings/metric.cpp \
	applications/metric/2d/all_mappings/metric_user.cpp \
	applications/metric/2d/all_mappings/metric_diagnostics.cpp \
	applications/metric/2d/all_mappings/tag4refinement.f \
	applications/metric/2d/all_mappings/metric_output.f \
	applications/metric/2d/all_mappings/initialize.f \
	applications/metric/2d/all_mappings/compute_error.f \
	applications/metric/2d/all_mappings/setprob.f \
	applications/metric/2d/all_mappings/diagnostics.f \
	applications/metric/2d/all_mappings/fclaw2d_map_nomap.c \
	applications/metric/2d/all_mappings/fclaw2d_map_cart.c \
	applications/metric/2d/all_mappings/fclaw2d_map_fivepatch.c \
	applications/metric/2d/all_mappings/fclaw2d_map_pillowdisk.c \
	applications/metric/2d/all_mappings/fclaw2d_map_pillowdisk5.c \
	applications/metric/2d/all_mappings/fclaw2d_map_squareddisk.c \
	applications/metric/2d/all_mappings/fclaw2d_map_pillowsphere.c \
	applications/metric/2d/all_mappings/fclaw2d_map_cubedsphere.c \
	applications/metric/2d/all_mappings/fclaw2d_map_torus.c

applications_metric_2d_all_mappings_metric_CPPFLAGS = \
	        $(AM_CPPFLAGS)

# This isn't used;  it is just included here to make sure that we are not
# relying on any clawpack features that are not available.
applications_advection_2d_swirl_swirl_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
