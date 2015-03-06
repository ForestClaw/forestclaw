## Makefile.am in applications/amrclaw/euler/2d/quadrants

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/amrclaw/euler/2d/quadrants/quadrants

applications_amrclaw_euler_2d_quadrants_quadrants_SOURCES = \
	applications/amrclaw/euler/2d/quadrants/quadrants.cpp \
	applications/amrclaw/euler/2d/quadrants/quadrants_user.cpp \
	applications/amrclaw/euler/2d/quadrants/fclaw2d_map_nomap.c \
	applications/amrclaw/euler/2d/quadrants/setprob.f \
	applications/amrclaw/euler/2d/quadrants/user_4.x/qinit.f \
	applications/amrclaw/euler/2d/rp/rpn2eu3.f \
	applications/amrclaw/euler/2d/rp/rpn2eu4.f \
	applications/amrclaw/euler/2d/rp/rpn2eu5.f \
	applications/amrclaw/euler/2d/rp/rpt2eu.f

applications_amrclaw_euler_2d_quadrants_quadrants_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_amrclaw_euler_2d_quadrants_quadrants_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
