## Makefile.am in applications/metric/2d/mesh

bin_PROGRAMS += applications/metric/2d/mesh/mesh

applications_metric_2d_mesh_mesh_SOURCES = \
	applications/metric/2d/mesh/mesh_user.cpp \
	applications/metric/2d/mesh/mesh_user.h \
	applications/metric/2d/mesh/mesh.cpp \
	applications/metric/2d/mesh/fclaw2d_map_nomap.c \
	applications/metric/2d/mesh/psi.f \
	applications/metric/2d/mesh/setprob.f \
	applications/metric/2d/all/clawpack46_setaux.f \
	applications/metric/2d/all/clawpack_user.h \
	applications/metric/2d/all/clawpack5_setaux.f90 \
	applications/metric/2d/all/clawpack46_rpn2adv.f \
	applications/metric/2d/all/clawpack46_rpt2adv.f \
	applications/metric/2d/all/clawpack5_rpn2adv.f90 \
	applications/metric/2d/all/clawpack5_rpt2adv.f90 \
	applications/metric/2d/mesh/user_4.6/qinit.f \
	applications/metric/2d/mesh/user_4.6/b4step2.f \
	applications/metric/2d/mesh/user_5.0/qinit.f90 \
	applications/metric/2d/mesh/user_5.0/b4step2.f90

applications_metric_2d_mesh_mesh_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS) \
	$(FCLAW_DUMMY_CPPFLAGS)

applications_metric_2d_mesh_mesh_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_DUMMY_LDADD)  $(FCLAW_CLAWPACK5_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fclaw2d_clawpack46/libclawpack5.la
