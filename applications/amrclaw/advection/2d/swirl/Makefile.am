## Makefile.am in applications/amrclaw/advection/2d/swirl

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/amrclaw/advection/2d/swirl/swirl

applications_amrclaw_advection_2d_swirl_swirl_SOURCES = \
	applications/amrclaw/advection/2d/swirl/swirl_user.cpp \
	applications/amrclaw/advection/2d/swirl/swirl.cpp \
	applications/amrclaw/advection/2d/swirl/fclaw2d_map_nomap.c \
	applications/amrclaw/advection/2d/swirl/psi.f \
	applications/amrclaw/advection/2d/swirl/setprob.f \
	applications/amrclaw/advection/2d/swirl/user_4.x/tag4refinement.f \
	applications/amrclaw/advection/2d/swirl/user_4.x/qinit.f \
	applications/amrclaw/advection/2d/swirl/user_4.x/b4step2.f \
	applications/amrclaw/advection/2d/swirl/user_4.x/setaux.f \
	applications/amrclaw/advection/2d/swirl/user_4.x/swirl_output.f \
	applications/amrclaw/advection/2d/rp/rpn2ad1.f \
	applications/amrclaw/advection/2d/rp/rpt2ad1.f

applications_amrclaw_advection_2d_swirl_swirl_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_amrclaw_advection_2d_swirl_swirl_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
