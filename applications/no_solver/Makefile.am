## Makefile.am in applications/no_solver

bin_PROGRAMS += applications/no_solver/no_solver

applications_no_solver_no_solver_SOURCES = \
	applications/no_solver/no_solver.cpp \
    applications/no_solver/no_solver_user.cpp \
    applications/no_solver/no_solver_options.c \
    applications/no_solver/gem3d_output_mesh.c \
	applications/no_solver/user_4.6/initialize.f

applications_no_solver_no_solver_CPPFLAGS = \
	       $(AM_CPPFLAGS) \
           $(FCLAW_CLAWPATCH_CPPFLAGS)

# Not needed, but included to check compile
applications_no_solver_no_solver_LDADD = \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
