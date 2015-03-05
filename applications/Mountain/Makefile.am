## Makefile.am in applications/Mountain

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/Mountain/mountain

applications_Mountain_mountain_SOURCES = \
	applications/Mountain/mountain_user.cpp \
	applications/Mountain/mountain.cpp \
	applications/Mountain/mountain_fort.f \
	applications/Mountain/setprob.f \
	applications/Mountain/psi.f \
	applications/Mountain/fclaw2d_map_mountain.c \
	applications/Mountain/fclaw2d_map_identity.c \
	applications/Mountain/user_4.x/qinit.f \
	applications/Mountain/user_4.x/setaux.f \
	applications/Mountain/user_4.x/rpn2ad1.f \
	applications/Mountain/user_4.x/rpt2ad1.f \
	applications/Mountain/user_4.x/tag4refinement.f \
	applications/Mountain/user_4.x/tag4coarsening.f

applications_Mountain_mountain_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_Mountain_mountain_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
