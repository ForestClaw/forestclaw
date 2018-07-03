bin_PROGRAMS += applications/clawpack/burgers/2d/pwconst/pwconst

applications_clawpack_burgers_2d_pwconst_pwconst_SOURCES = \
	applications/clawpack/burgers/2d/pwconst/pwconst.cpp \
	applications/clawpack/burgers/2d/pwconst/pwconst_options.c \
	applications/clawpack/burgers/2d/pwconst/pwconst_user.cpp \
	applications/clawpack/burgers/2d/pwconst/pwconst_user.h \
	applications/clawpack/burgers/2d/pwconst/user_4.6/qinit.f \
	applications/clawpack/burgers/2d/pwconst/user_5.0/qinit.f90 \
	applications/clawpack/burgers/2d/rp/clawpack46_rpn2_burgers.f  \
	applications/clawpack/burgers/2d/rp/clawpack46_rpt2_burgers.f  \
	applications/clawpack/burgers/2d/rp/clawpack5_rpn2_burgers.f90 \
	applications/clawpack/burgers/2d/rp/clawpack5_rpt2_burgers.f90

EXTRA_DIST += \
	applications/clawpack/burgers/2d/pwconst/fclaw_options.ini \
	applications/clawpack/burgers/2d/pwconst/regression.ini \
	applications/clawpack/burgers/2d/pwconst/afterframe.m \
	applications/clawpack/burgers/2d/pwconst/setplot2.m

applications_clawpack_burgers_2d_pwconst_pwconst_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_burgers_2d_pwconst_pwconst_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
