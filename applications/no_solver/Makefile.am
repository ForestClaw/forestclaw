## Makefile.am in applications/no_solver

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/no_solver/no_solver

applications_no_solver_no_solver_SOURCES = \
	applications/no_solver/no_solver.cpp \
	applications/no_solver/no_solver_user.cpp \
	applications/no_solver/fclaw2d_map_nomap.c \
	applications/no_solver/fclaw2d_map_cart.c \
	applications/no_solver/fclaw2d_map_fivepatch.c \
	applications/no_solver/fclaw2d_map_pillowdisk.c \
	applications/no_solver/fclaw2d_map_pillowdisk5.c \
	applications/no_solver/tag4refinement.f \
	applications/no_solver/initialize.f

applications_no_solver_no_solver_CPPFLAGS = \
	       $(AM_CPPFLAGS)

# applications_no_solver_no_solver_LDADD = $(LDADD) $(FLIBS)

# Not needed, but included to check compile
applications_no_solver_no_solver_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
