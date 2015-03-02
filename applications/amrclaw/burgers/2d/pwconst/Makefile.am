## Makefile.am in applications/amrclaw/burgers_2d_pwconst

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/amrclaw/burgers/2d/pwconst/pwconst

applications_amrclaw_burgers_2d_pwconst_pwconst_SOURCES = \
	applications/amrclaw/burgers/2d/pwconst/pwconst_user.cpp \
	applications/amrclaw/burgers/2d/pwconst/pwconst.cpp \
	applications/amrclaw/burgers/2d/pwconst/fclaw2d_map_nomap.c \
	applications/amrclaw/burgers/2d/pwconst/user_4.x/qinit.f \
	applications/amrclaw/burgers/2d/rp/rpn2bu.f \
	applications/amrclaw/burgers/2d/rp/rpt2bu.f


applications_amrclaw_burgers_2d_pwconst_pwconst_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_amrclaw_burgers_2d_pwconst_pwconst_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
