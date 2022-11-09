## Makefile.am in applications/clawpack/advection/2d/filament_swirl

bin_PROGRAMS += applications/clawpack/advection/2d/filament_swirl/filament_swirl

TESTS += applications/clawpack/advection/2d/filament_swirl/regressions.sh

applications_clawpack_advection_2d_filament_swirl_filament_swirl_SOURCES = \
	applications/clawpack/advection/2d/filament_swirl/filament_swirl.cpp \
	applications/clawpack/advection/2d/filament_swirl/psi.f \
	applications/clawpack/advection/2d/filament_swirl/filament/filament.cpp \
	applications/clawpack/advection/2d/filament_swirl/filament/filament_options.c \
	applications/clawpack/advection/2d/filament_swirl/filament/filament_user.cpp \
	applications/clawpack/advection/2d/filament_swirl/filament/filament_user.h \
	applications/clawpack/advection/2d/filament_swirl/filament/psi.f \
	applications/clawpack/advection/2d/filament_swirl/filament/fdisc.f \
	applications/clawpack/advection/2d/filament_swirl/filament/setprob.f \
	applications/clawpack/advection/2d/filament_swirl/filament/user_4.6/qinit.f \
	applications/clawpack/advection/2d/filament_swirl/filament/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/filament_swirl/swirl/swirl_user.cpp \
	applications/clawpack/advection/2d/filament_swirl/swirl/swirl_user.h \
	applications/clawpack/advection/2d/filament_swirl/swirl/swirl_options.c \
	applications/clawpack/advection/2d/filament_swirl/swirl/swirl.cpp \
	applications/clawpack/advection/2d/filament_swirl/swirl/psi.f \
	applications/clawpack/advection/2d/filament_swirl/swirl/setprob.f \
	applications/clawpack/advection/2d/filament_swirl/swirl/user_4.6/qinit.f \
	applications/clawpack/advection/2d/filament_swirl/swirl/user_4.6/b4step2.f \
	applications/clawpack/advection/2d/filament_swirl/swirl/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/filament_swirl/swirl/user_5.0/b4step2.f90 \
	applications/clawpack/advection/2d/all/advection_patch_setup_manifold.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_cart.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_fivepatch.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_bilinear.c \
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
	applications/clawpack/advection/2d/filament_swirl/fclaw_options.ini \
	applications/clawpack/advection/2d/filament_swirl/regression.ini \
	applications/clawpack/advection/2d/filament_swirl/filament/afterframe.m \
	applications/clawpack/advection/2d/filament_swirl/filament/setplot2.m \
	applications/clawpack/advection/2d/filament_swirl/filament/mapc2m.m \
	applications/clawpack/advection/2d/filament_swirl/filament/filament_soln.m \
	applications/clawpack/advection/2d/filament_swirl/filament/create_filament_soln_tikz.m \
	applications/clawpack/advection/2d/filament_swirl/filament/filament.ipynb \
	applications/clawpack/advection/2d/filament_swirl/filament/make_plots.py \
	applications/clawpack/advection/2d/filament_swirl/swirl/afterframe.m \
	applications/clawpack/advection/2d/filament_swirl/swirl/setplot2.m \
	applications/clawpack/advection/2d/filament_swirl/swirl/swirl.ipynb \
	applications/clawpack/advection/2d/filament_swirl/swirl/make_plots.py

applications_clawpack_advection_2d_filament_swirl_filament_swirl_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)
	    

# applications_clawpack_advection_2d_filament_swirl_filament_swirl_CPPFLAGS += -D MPI_DEBUG


applications_clawpack_advection_2d_filament_swirl_filament_swirl_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)