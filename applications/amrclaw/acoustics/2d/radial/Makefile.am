## Makefile.am in applications/amrclaw/acoustics_2d_radial

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/amrclaw/acoustics/2d/radial/radial

applications_amrclaw_acoustics_2d_radial_radial_SOURCES = \
	applications/amrclaw/acoustics/2d/radial/radial_user.cpp \
	applications/amrclaw/acoustics/2d/radial/radial.cpp \
	applications/amrclaw/acoustics/2d/radial/fclaw2d_map_nomap.c \
	applications/amrclaw/acoustics/2d/radial/fclaw2d_map_pillowdisk5.c \
	applications/amrclaw/acoustics/2d/radial/setprob.f \
	applications/amrclaw/acoustics/2d/radial/user_4.6/qinit.f \
	applications/amrclaw/acoustics/2d/radial/user_5.0/qinit.f90 \
	applications/amrclaw/acoustics/2d/radial/user_5.0/tag4refinement.f90 \
	applications/amrclaw/acoustics/2d/radial/user_5.0/tag4coarsening.f90 \
	applications/amrclaw/acoustics/2d/radial/user_5.0/setaux.f90 \
	applications/amrclaw/acoustics/2d/rp/rpn2ac.f \
	applications/amrclaw/acoustics/2d/rp/rpt2ac.f \
	applications/amrclaw/acoustics/2d/rp/rpn2_acoustics.f90 \
	applications/amrclaw/acoustics/2d/rp/rpt2_acoustics.f90 \
	applications/amrclaw/acoustics/2d/rp/rpn2_acoustics_mapped.f90 \
	applications/amrclaw/acoustics/2d/rp/rpt2_acoustics_mapped.f90


applications_amrclaw_acoustics_2d_radial_radial_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_amrclaw_acoustics_2d_radial_radial_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
