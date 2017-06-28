## Makefile.am in applications/advection/2d/filament

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/advection/2d/filament/filament


applications_advection_2d_filament_filament_SOURCES = \
	applications/advection/2d/filament/filament_user.cpp \
	applications/advection/2d/filament/filament_user.h \
	applications/advection/2d/filament/filament.cpp \
	applications/advection/2d/filament/fclaw2d_map_cart.c \
	applications/advection/2d/filament/fclaw2d_map_fivepatch.c \
	applications/advection/2d/filament/psi.f \
	applications/advection/2d/filament/fdisc.f \
	applications/advection/2d/filament/setprob.f \
	applications/advection/2d/filament/user_4.6/qinit.f \
	applications/advection/2d/filament/user_4.6/tag4refinement.f \
	applications/advection/2d/filament/user_4.6/tag4coarsening.f \
	applications/advection/2d/filament/user_5.0/qinit.f90 \
	applications/advection/2d/filament/user_5.0/tag4refinement.f \
	applications/advection/2d/filament/user_5.0/tag4coarsening.f \
	applications/advection/2d/all/user46_setaux_manifold.f \
	applications/advection/2d/all/user5_setaux_manifold.f90 \
	applications/advection/2d/all/clawpack46_setaux.f \
	applications/advection/2d/all/clawpack5_setaux.f90 \
	applications/advection/2d/all/clawpack46_rpn2adv.f \
	applications/advection/2d/all/clawpack46_rpt2adv.f \
	applications/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/advection/2d/all/clawpack5_rpn2adv.f90 \
	applications/advection/2d/all/clawpack5_rpt2adv.f90 \
	applications/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/advection/2d/all/clawpack5_rpt2adv_manifold.f90


applications_advection_2d_filament_filament_CPPFLAGS = -I../all/ \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

# applications_advection_2d_filament_filament_CPPFLAGS += -D MPI_DEBUG


applications_advection_2d_filament_filament_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD)  $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
