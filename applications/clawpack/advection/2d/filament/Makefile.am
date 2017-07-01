## Makefile.am in applications/clawpack/advection/2d/filament

bin_PROGRAMS += applications/clawpack/advection/2d/filament/filament


applications_clawpack_advection_2d_filament_filament_SOURCES = \
	applications/clawpack/advection/2d/filament/filament_user.cpp \
	applications/clawpack/advection/2d/filament/filament_user.h \
	applications/clawpack/advection/2d/filament/filament.cpp \
	applications/clawpack/advection/2d/filament/fclaw2d_map_cart.c \
	applications/clawpack/advection/2d/filament/fclaw2d_map_fivepatch.c \
	applications/clawpack/advection/2d/filament/psi.f \
	applications/clawpack/advection/2d/filament/fdisc.f \
	applications/clawpack/advection/2d/filament/setprob.f \
	applications/clawpack/advection/2d/filament/user_4.6/qinit.f \
	applications/clawpack/advection/2d/filament/user_4.6/tag4refinement.f \
	applications/clawpack/advection/2d/filament/user_4.6/tag4coarsening.f \
	applications/clawpack/advection/2d/filament/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/filament/user_5.0/tag4refinement.f \
	applications/clawpack/advection/2d/filament/user_5.0/tag4coarsening.f \
	applications/clawpack/advection/2d/all/user46_setaux_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack46_setaux.f \
	applications/clawpack/advection/2d/all/clawpack5_setaux.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv.f \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90


applications_clawpack_advection_2d_filament_filament_CPPFLAGS = -I../all/ \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

# applications_clawpack_advection_2d_filament_filament_CPPFLAGS += -D MPI_DEBUG


applications_clawpack_advection_2d_filament_filament_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD)  $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
