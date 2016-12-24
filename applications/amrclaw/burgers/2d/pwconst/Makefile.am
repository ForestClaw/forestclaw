bin_PROGRAMS += applications/amrclaw/burgers/2d/pwconst/pwconst

applications_amrclaw_burgers_2d_pwconst_pwconst_SOURCES = \
	applications/amrclaw/burgers/2d/pwconst/pwconst_user.cpp \
	applications/amrclaw/burgers/2d/pwconst/pwconst_user.h \
	applications/amrclaw/burgers/2d/pwconst/pwconst.cpp \
	applications/amrclaw/burgers/2d/pwconst/fclaw2d_map_nomap.c \
	applications/amrclaw/burgers/2d/pwconst/user_4.6/qinit.f \
	applications/amrclaw/burgers/2d/pwconst/user_5.0/qinit.f90 \
	applications/amrclaw/burgers/2d/rp/clawpack46_rpn2_burgers.f  \
	applications/amrclaw/burgers/2d/rp/clawpack46_rpt2_burgers.f  \
	applications/amrclaw/burgers/2d/rp/clawpack5_rpn2_burgers.f90 \
	applications/amrclaw/burgers/2d/rp/clawpack5_rpt2_burgers.f90


applications_amrclaw_burgers_2d_pwconst_pwconst_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_amrclaw_burgers_2d_pwconst_pwconst_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
