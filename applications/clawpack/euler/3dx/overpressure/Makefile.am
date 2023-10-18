## Makefile.am in applications/clawpack/euler/3dx/overpressure

bin_PROGRAMS += applications/clawpack/euler/3dx/overpressure/overpressure

TESTS += applications/clawpack/euler/3dx/overpressure/regressions.sh

applications_clawpack_euler_3d_overpressure_overpressure_SOURCES = \
	applications/clawpack/euler/3dx/overpressure/overpressure.cpp \
	applications/clawpack/euler/3dx/overpressure/overpressure_options.c \
	applications/clawpack/euler/3dx/overpressure/overpressure_user.cpp \
	applications/clawpack/euler/3dx/overpressure/overpressure_user.h \
	applications/clawpack/euler/3dx/overpressure/overpressure_map_extrude.c \
	applications/clawpack/euler/3dx/overpressure/fdisc.f90 \
	applications/clawpack/euler/3dx/overpressure/setprob_module.f90 \
	applications/clawpack/euler/3dx/overpressure/user_4.6/qinit.f90 \
	applications/clawpack/euler/3dx/overpressure/user_4.6/bc3_mapped.f90 \
	applications/clawpack/euler/3dx/all/euler3d_pressure_exceeds_th.f90 \
	applications/clawpack/euler/3dx/all/euler3d_setaux_manifold.f90 \
	applications/clawpack/euler/3dx/all/euler3d_user_fort.h \
	applications/clawpack/euler/3dx/rp/gethexinfo.f90 \
	applications/clawpack/euler/3dx/rp/rpn3_euler.f90 \
	applications/clawpack/euler/3dx/rp/rpn3_euler_mapped.f90 \
	applications/clawpack/euler/3dx/rp/rpt3_euler.f90 \
	applications/clawpack/euler/3dx/rp/rpt3_euler_mapped.f90 \
	applications/clawpack/euler/3dx/rp/rptt3_euler.f90 \
	applications/clawpack/euler/3dx/rp/rptt3_euler_mapped.f90 \
	applications/clawpack/advection/2d/all/fclaw2d_map_cart.c  \
	applications/clawpack/advection/2d/all/fclaw2d_map_latlong.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_cubedsphere.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_pillowsphere.c 

EXTRA_DIST += \
	applications/clawpack/euler/3dx/overpressure/fclaw_options.ini \
	applications/clawpack/euler/3dx/overpressure/regression.ini \
	applications/clawpack/euler/3dx/overpressure/regressions.sh \
	applications/clawpack/euler/3dx/overpressure/afterframe.m \
	applications/clawpack/euler/3dx/overpressure/mapc2p.m \
	applications/clawpack/euler/3dx/overpressure/read_vars.m \
	applications/clawpack/euler/3dx/overpressure/setplot3.m

applications_clawpack_euler_3d_overpressure_overpressure_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_euler_3d_overpressure_overpressure_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)


