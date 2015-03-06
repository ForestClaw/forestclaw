## Makefile.am in applications/amrclaw/acoustics_2d_interface

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/amrclaw/acoustics/2d/interface/interface

applications_amrclaw_acoustics_2d_interface_interface_SOURCES = \
	applications/amrclaw/acoustics/2d/interface/interface_user.cpp \
	applications/amrclaw/acoustics/2d/interface/interface.cpp \
	applications/amrclaw/acoustics/2d/interface/p0.f \
	applications/amrclaw/acoustics/2d/interface/fclaw2d_map_nomap.c \
	applications/amrclaw/acoustics/2d/interface/setprob.f \
	applications/amrclaw/acoustics/2d/interface/user_4.x/setaux.f \
	applications/amrclaw/acoustics/2d/interface/user_4.x/qinit.f \
	applications/amrclaw/acoustics/2d/interface/user_4.x/tag4refinement.f \
	applications/amrclaw/acoustics/2d/rp/rpn2acv.f \
	applications/amrclaw/acoustics/2d/rp/rpt2acv.f


applications_amrclaw_acoustics_2d_interface_interface_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_amrclaw_acoustics_2d_interface_interface_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
