## Makefile.am in applications/clawpack/advection/3dx/filament

bin_PROGRAMS += applications/clawpack/advection/3dx/filament/filament

TESTS += applications/clawpack/advection/3dx/filament/regressions.sh

applications_clawpack_advection_3d_filament_filament_SOURCES = \
	applications/clawpack/advection/3dx/filament/filament.cpp \
	applications/clawpack/advection/3dx/filament/filament_options.c \
	applications/clawpack/advection/3dx/filament/filament_user.cpp \
	applications/clawpack/advection/3dx/filament/filament_map_extrude.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_cart.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_fivepatch.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_bilinear.c \
	applications/clawpack/advection/3dx/filament/user_exceeds_th.f90 \
	applications/clawpack/advection/3dx/filament/psi.f \
	applications/clawpack/advection/3dx/filament/fdisc.f \
	applications/clawpack/advection/3dx/filament/setprob.f \
	applications/clawpack/advection/3dx/filament/user_4.6/qinit.f90 \
	applications/clawpack/advection/3dx/filament/user_4.6/setaux.f90 \
	applications/clawpack/advection/3dx/all/advection_patch_setup_manifold.c \
	applications/clawpack/advection/3dx/all/claw3_setaux_manifold.f90 \
	applications/clawpack/advection/3dx/all/rpn3_vc_advection.f90 \
	applications/clawpack/advection/3dx/all/rpt3_vc_advection.f90 \
	applications/clawpack/advection/3dx/all/rptt3_vc_advection.f90

EXTRA_DIST += \
	applications/clawpack/advection/3dx/filament/fclaw_options.ini \
	applications/clawpack/advection/3dx/filament/mapping.ini \
	applications/clawpack/advection/3dx/filament/regression.ini \
	applications/clawpack/advection/3dx/filament/regression_map.ini \
	applications/clawpack/advection/3dx/filament/afterframe.m \
	applications/clawpack/advection/3dx/filament/setplot3.m  \
	applications/clawpack/advection/3dx/filament/mapc2p.m \
	applications/clawpack/advection/3dx/filament/read_vars.m

applications_clawpack_advection_3d_filament_filament_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK3_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)
	    

# applications_clawpack_advection_3d_filament_filament_CPPFLAGS += -D MPI_DEBUG


applications_clawpack_advection_3d_filament_filament_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK3_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)