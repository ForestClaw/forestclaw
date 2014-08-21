## Makefile.am in applications/advection/2d/filament

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/advection/2d/filament/filament

applications_advection_2d_filament_filament_SOURCES = \
	applications/advection/2d/filament/filament_user.cpp \
	applications/advection/2d/filament/filament.cpp \
	applications/advection/2d/filament/fclaw2d_map_cart.c \
	applications/advection/2d/filament/fclaw2d_map_pillowdisk.c \
	applications/advection/2d/filament/fclaw2d_map_squareddisk.c \
	applications/advection/2d/filament/psi.f \
	applications/advection/2d/filament/setprob.f \
	applications/advection/2d/filament/user_4.6/tag4refinement.f \
	applications/advection/2d/filament/user_4.6/rpn2noncons.f \
	applications/advection/2d/filament/user_4.6/rpt2noncons.f \
	applications/advection/2d/filament/user_4.6/qinit.f \
	applications/advection/2d/filament/user_4.6/setaux_manifold.f \
	applications/advection/2d/filament/user_4.6/filament_output.f


applications_advection_2d_filament_filament_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

# applications_advection_2d_filament_filament_CPPFLAGS += -D TRAPFPE
# applications_advection_2d_filament_filament_CPPFLAGS += -D TRAPFPE_OSX

# applications_advection_2d_filament_filament_CPPFLAGS += -D MPI_DEBUG


applications_advection_2d_filament_filament_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
## @top_builddir@/src/solvers/fclaw2d_clawpack4.6/libclawpack4.6.la
