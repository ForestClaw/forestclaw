## Makefile.am in applications/advection/2d/swirl

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

## We should probably include .am files only from the toplevel Makefile.am,
## and not recursively as is done here.
## include src/solvers/amr_waveprop/Makefile.am

## We need to find a solution about how to compile the waveprop files without
## creating a name conflict with other solvers.

bin_PROGRAMS += applications/advection/2d/swirl/swirl

applications_advection_2d_swirl_swirl_SOURCES = \
	applications/advection/2d/swirl/swirl_user.cpp \
	applications/advection/2d/swirl/swirl.cpp \
	applications/advection/2d/swirl/ghost_twopatch.f \
	applications/advection/2d/swirl/tag4refinement.f \
	applications/advection/2d/swirl/psi.f \
	applications/advection/2d/swirl/rpn2ad1.f \
	applications/advection/2d/swirl/rpt2ad1.f \
	applications/advection/2d/swirl/qinit.f \
	applications/advection/2d/swirl/b4step2.f \
	applications/advection/2d/swirl/setaux.f \
	applications/advection/2d/swirl/setprob.f \
	applications/advection/2d/swirl/swirl_output.f

applications_advection_2d_swirl_swirl_CPPFLAGS = \
	$(AM_CPPFLAGS) -I@top_srcdir@/src/solvers/amr_waveprop

applications_advection_2d_swirl_swirl_LDADD = @top_builddir@/local/lib/libwaveprop.la $(LDADD)
