## Makefile.am in applications/euler/2d/quadrants

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/euler/2d/quadrants/quadrants

applications_euler_2d_quadrants_quadrants_SOURCES = \
	applications/euler/2d/quadrants/quadrants.cpp \
	applications/euler/2d/quadrants/quadrants_user.cpp \
	applications/euler/2d/quadrants/tag4refinement.f \
	applications/euler/2d/quadrants/quadrants_output.f \
	applications/euler/2d/quadrants/fclaw2d_map_nomap.c \
	applications/euler/2d/quadrants/fclaw2d_map_cart.c \
	applications/euler/2d/quadrants/fclaw2d_map_fivepatch.c \
	applications/euler/2d/quadrants/rpn2eu4.f \
	applications/euler/2d/quadrants/rpt2eu.f \
	applications/euler/2d/quadrants/qinit.f \
	applications/euler/2d/quadrants/setprob.f

applications_euler_2d_quadrants_quadrants_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

# applications_euler_2d_quadrants_quadrants_CPPFLAGS += -D TRAPFPE

# applications_euler_2d_quadrants_quadrants_CPPFLAGS += -D MPI_DEBUG

applications_euler_2d_quadrants_quadrants_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
