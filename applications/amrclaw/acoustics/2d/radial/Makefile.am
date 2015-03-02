## Makefile.am in applications/amrclaw/acoustics_2d_radial

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/amrclaw/acoustics/2d/radial/radial

applications_amrclaw_acoustics_2d_radial_radial_SOURCES = \
	applications/amrclaw/acoustics/2d/radial/radial_user.cpp \
	applications/amrclaw/acoustics/2d/radial/radial.cpp \
	applications/amrclaw/acoustics/2d/radial/fclaw2d_map_nomap.c \
	applications/amrclaw/acoustics/2d/radial/fclaw2d_map_squareddisk.c \
	applications/amrclaw/acoustics/2d/radial/setprob.f \
	applications/amrclaw/acoustics/2d/radial/user_4.x/tag4refinement.f \
	applications/amrclaw/acoustics/2d/radial/user_4.x/qinit.f \
	applications/amrclaw/acoustics/2d/rp/rpn2ac.f \
	applications/amrclaw/acoustics/2d/rp/rpt2ac.f


applications_amrclaw_acoustics_2d_radial_radial_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_amrclaw_acoustics_2d_radial_radial_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
