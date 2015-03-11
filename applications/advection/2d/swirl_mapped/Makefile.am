## Makefile.am in applications/advection/2d/swirl_mapped

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/advection/2d/swirl_mapped/swirl_mapped

applications_advection_2d_swirl_mapped_swirl_mapped_SOURCES = \
	applications/advection/2d/swirl_mapped/swirl_user.cpp \
	applications/advection/2d/swirl_mapped/swirl.cpp \
	applications/advection/2d/swirl_mapped/fclaw2d_map_nomap.c \
	applications/advection/2d/swirl_mapped/fclaw2d_map_cart.c \
	applications/advection/2d/swirl_mapped/fclaw2d_map_identity.c \
	applications/advection/2d/swirl_mapped/fclaw2d_map_fivepatch.c \
	applications/advection/2d/swirl_mapped/setprob.f \
	applications/advection/2d/swirl_mapped/psi.f \
	applications/advection/2d/swirl_mapped/user_4.6/setaux_manifold.f \
	applications/advection/2d/swirl_mapped/user_4.6/b4step2_manifold.f \
	applications/advection/2d/swirl_mapped/user_4.6/qinit.f \
	applications/advection/2d/swirl_mapped/user_4.6/setaux.f \
	applications/advection/2d/swirl_mapped/user_4.6/rpn2ad1.f \
	applications/advection/2d/swirl_mapped/user_4.6/rpt2ad1.f

applications_advection_2d_swirl_mapped_swirl_mapped_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_advection_2d_swirl_mapped_swirl_mapped_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
