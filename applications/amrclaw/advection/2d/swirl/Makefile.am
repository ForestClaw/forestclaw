# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/amrclaw/advection/2d/swirl/swirl

applications_amrclaw_advection_2d_swirl_swirl_SOURCES = \
	applications/amrclaw/advection/2d/swirl/swirl_user.cpp \
	applications/amrclaw/advection/2d/swirl/swirl_user.h \
	applications/amrclaw/advection/2d/swirl/swirl.cpp \
	applications/amrclaw/advection/2d/swirl/psi.f \
	applications/amrclaw/advection/2d/swirl/setprob.f \
	applications/amrclaw/advection/2d/swirl/user_4.6/qinit.f \
	applications/amrclaw/advection/2d/swirl/user_4.6/b4step2.f \
	applications/amrclaw/advection/2d/swirl/user_4.6/setaux.f \
	applications/amrclaw/advection/2d/swirl/user_5.0/qinit.f90 \
	applications/amrclaw/advection/2d/swirl/user_5.0/b4step2.f90 \
	applications/amrclaw/advection/2d/swirl/user_5.0/setaux.f90 \
	applications/amrclaw/advection/2d/rp/clawpack46_rpn2.f \
	applications/amrclaw/advection/2d/rp/clawpack46_rpt2.f \
	applications/amrclaw/advection/2d/rp/clawpack5_rpn2.f90 \
	applications/amrclaw/advection/2d/rp/clawpack5_rpt2.f90

applications_amrclaw_advection_2d_swirl_swirl_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_amrclaw_advection_2d_swirl_swirl_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
