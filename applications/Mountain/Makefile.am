## Makefile.am in applications/Mountain

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/Mountain/mountain

applications_Mountain_mountain_SOURCES = \
	applications/Mountain/mountain_user.cpp \
	applications/Mountain/mountain.cpp \
	applications/Mountain/mountain_fort.f \
	applications/Mountain/setaux_manifold.f \
	applications/Mountain/fclaw2d_map_mountain.c \
	applications/Mountain/fclaw2d_map_identity.c \
	applications/Mountain/setprob.f \
	applications/Mountain/psi.f \
	applications/Mountain/user_4.6/b4step2_manifold.f \
	applications/Mountain/user_4.6/tag4refinement.f \
	applications/Mountain/user_4.6/rpn2ad1.f \
	applications/Mountain/user_4.6/rpt2ad1.f \
	applications/Mountain/user_4.6/qinit.f \
	applications/Mountain/user_4.6/mountain_output.f

applications_Mountain_mountain_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

# applications_advection_2d_mountain_mapped_swirl_CPPFLAGS += -D MPI_DEBUG

applications_Mountain_mountain_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
