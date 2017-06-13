## Makefile.am in applications/amrclaw/shallow/2d/radialdam

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/amrclaw/shallow/2d/radialdam/radialdam

applications_amrclaw_shallow_2d_radialdam_radialdam_SOURCES = \
	applications/amrclaw/shallow/2d/radialdam/radialdam.cpp \
	applications/amrclaw/shallow/2d/radialdam/radialdam_user.cpp \
	applications/amrclaw/shallow/2d/radialdam/radialdam_user.h \
	applications/amrclaw/shallow/2d/radialdam/fclaw2d_map_nomap.c \
	applications/amrclaw/shallow/2d/radialdam/fclaw2d_map_pillowdisk.c \
	applications/amrclaw/shallow/2d/radialdam/setprob.f \
	applications/amrclaw/shallow/2d/radialdam/fdisc.f \
	applications/amrclaw/shallow/2d/radialdam/user_4.6/qinit.f \
	applications/amrclaw/shallow/2d/radialdam/user_5.0/qinit.f90 \
	applications/amrclaw/shallow/2d/radialdam/user_5.0/setaux.f90 \
	applications/amrclaw/shallow/2d/radialdam/user_5.0/tag4refinement.f90 \
	applications/amrclaw/shallow/2d/radialdam/user_5.0/tag4coarsening.f90 \
	applications/amrclaw/shallow/2d/rp/clawpack_user.h \
	applications/amrclaw/shallow/2d/rp/clawpack46_rpn2.f \
	applications/amrclaw/shallow/2d/rp/clawpack46_rpt2.f \
	applications/amrclaw/shallow/2d/rp/clawpack5_rpn2.f90 \
	applications/amrclaw/shallow/2d/rp/clawpack5_rpt2.f90 \
	applications/amrclaw/shallow/2d/rp/clawpack5_rpn2_manifold.f90 \
	applications/amrclaw/shallow/2d/rp/clawpack5_rpt2_manifold.f90

applications_amrclaw_shallow_2d_radialdam_radialdam_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_amrclaw_shallow_2d_radialdam_radialdam_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
