## Makefile.am in applications/metric/2d/sphere

bin_PROGRAMS += applications/metric/2d/sphere/sphere


applications_metric_2d_sphere_sphere_SOURCES = \
	applications/metric/2d/sphere/sphere_user.cpp \
	applications/metric/2d/sphere/sphere.cpp \
	applications/metric/2d/sphere/ghost_twopatch.f \
	applications/metric/2d/sphere/tag4refinement.f \
	applications/metric/2d/sphere/output_metric.f \
	applications/metric/2d/sphere/mapc2m.f \
	applications/metric/2d/sphere/setprob.f

applications_metric_2d_sphere_sphere_FFLAGS = \
	$(AM_FFLAGS) -I@top_srcdir@/src/fortran_source2d
