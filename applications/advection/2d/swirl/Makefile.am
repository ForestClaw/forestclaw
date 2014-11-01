## Makefile.am in applications/advection/2d/swirl

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/advection/2d/swirl/swirl

applications_advection_2d_swirl_swirl_SOURCES = \
	applications/advection/2d/swirl/swirl_user.cpp \
	applications/advection/2d/swirl/swirl.cpp \
	applications/advection/2d/swirl/fclaw2d_map_nomap.c \
	applications/advection/2d/swirl/psi.f \
	applications/advection/2d/swirl/user_4.6/tag4refinement.f \
	applications/advection/2d/swirl/user_4.6/rpn2ad1.f \
	applications/advection/2d/swirl/user_4.6/rpt2ad1.f \
	applications/advection/2d/swirl/user_4.6/qinit.f \
	applications/advection/2d/swirl/user_4.6/b4step2.f \
	applications/advection/2d/swirl/user_4.6/setaux.f \
	applications/advection/2d/swirl/user_4.6/setprob.f \
	applications/advection/2d/swirl/user_4.6/swirl_output.f

applications_advection_2d_swirl_swirl_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_advection_2d_swirl_swirl_CPPFLAGS += -D TRAPFPE

# applications_advection_2d_swirl_swirl_CPPFLAGS += -D MPI_DEBUG

applications_advection_2d_swirl_swirl_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fclaw2d_clawpack4.6/libclawpack4.6.la
