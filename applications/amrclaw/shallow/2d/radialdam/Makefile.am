## Makefile.am in applications/amrclaw/shallow/2d/radialdam

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/amrclaw/shallow/2d/radialdam/radialdam

applications_amrclaw_shallow_2d_radialdam_radialdam_SOURCES = \
	applications/amrclaw/shallow/2d/radialdam/radialdam.cpp \
	applications/amrclaw/shallow/2d/radialdam/radialdam_user.cpp \
	applications/amrclaw/shallow/2d/radialdam/fclaw2d_map_nomap.c \
	applications/amrclaw/shallow/2d/radialdam/setprob.f \
	applications/amrclaw/shallow/2d/radialdam/fdisc.f \
	applications/amrclaw/shallow/2d/radialdam/user_4.x/tag4refinement.f \
	applications/amrclaw/shallow/2d/radialdam/user_4.x/radialdam_output.f \
	applications/amrclaw/shallow/2d/radialdam/user_4.x/qinit.f \
	applications/amrclaw/shallow/2d/rp/rpn2sw.f \
	applications/amrclaw/shallow/2d/rp/rpt2sw.f

applications_amrclaw_shallow_2d_radialdam_radialdam_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_amrclaw_shallow_2d_radialdam_radialdam_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
