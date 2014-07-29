## Makefile.am in applications/advection/2d/swirl

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/advection/2d/no_solver/no_solver

applications_advection_2d_no_solver_no_solver_SOURCES = \
	applications/advection/2d/no_solver/no_solver.cpp \
	applications/advection/2d/no_solver/no_solver_user.cpp \
	applications/advection/2d/no_solver/user_4.6/ghost_twopatch.f \
	applications/advection/2d/no_solver/user_4.6/tag4refinement.f \
	applications/advection/2d/no_solver/user_4.6/initialize.f

applications_advection_2d_no_solver_no_solver_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)

applications_advection_2d_no_solver_no_solver_LDADD = $(LDADD) $(FLIBS)
