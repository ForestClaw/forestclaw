## Makefile.am in applications/no_solver

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/no_solver/no_solver

applications_no_solver_no_solver_SOURCES = \
	applications/no_solver/no_solver.cpp \
	applications/no_solver/no_solver_user.cpp \
	applications/no_solver/user_4.6/tag4refinement.f \
	applications/no_solver/user_4.6/initialize.f

# applications_no_solver_no_solver_CPPFLAGS = \
#       $(AM_CPPFLAGS)

applications_no_solver_no_solver_LDADD = $(LDADD) $(FLIBS)
