## Makefile.am in applications/clawpack/advection/2d/filament

bin_PROGRAMS += applications/clawpack/advection/2d/filament/filament

TESTS += applications/clawpack/advection/2d/filament/regressions.sh

applications_clawpack_advection_2d_filament_filament_SOURCES = \
	applications/clawpack/advection/2d/filament/filament.cpp \
	applications/clawpack/advection/2d/filament/filament_options.c \
	applications/clawpack/advection/2d/filament/filament_user.cpp \
	applications/clawpack/advection/2d/filament/filament_user.h \
	applications/clawpack/advection/2d/filament/psi.f \
	applications/clawpack/advection/2d/filament/fdisc.f \
	applications/clawpack/advection/2d/filament/setprob.f \
	applications/clawpack/advection/2d/filament/user_4.6/qinit.f \
	applications/clawpack/advection/2d/filament/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/all/advection_user.cpp \
	applications/clawpack/advection/2d/all/fclaw2d_map_cart.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_fivepatch.c \
	applications/clawpack/advection/2d/all/user46_setaux_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack46_setaux.f \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack5_setaux.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90

EXTRA_DIST += \
	applications/clawpack/advection/2d/filament/fclaw_options.ini \
	applications/clawpack/advection/2d/filament/regression.ini \
	applications/clawpack/advection/2d/filament/afterframe.m \
	applications/clawpack/advection/2d/filament/setplot2.m \
	applications/clawpack/advection/2d/filament/mapc2m.m \
	applications/clawpack/advection/2d/filament/filament_soln.m \
	applications/clawpack/advection/2d/filament/create_filament_soln_tikz.m \
	applications/clawpack/advection/2d/filament/filament.ipynb \
	applications/clawpack/advection/2d/filament/make_plots.py

applications_clawpack_advection_2d_filament_filament_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)
	    

# applications_clawpack_advection_2d_filament_filament_CPPFLAGS += -D MPI_DEBUG


applications_clawpack_advection_2d_filament_filament_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)