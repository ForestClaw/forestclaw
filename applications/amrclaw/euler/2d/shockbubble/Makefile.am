## Makefile.am in applications/amrclaw/euler/2d/shockbubble

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/amrclaw/euler/2d/shockbubble/shockbubble


applications_amrclaw_euler_2d_shockbubble_shockbubble_SOURCES = \
	applications/amrclaw/euler/2d/shockbubble/shockbubble.cpp \
	applications/amrclaw/euler/2d/shockbubble/shockbubble_user.h \
	applications/amrclaw/euler/2d/shockbubble/shockbubble_user.cpp \
	applications/amrclaw/euler/2d/shockbubble/fdisc.f \
	applications/amrclaw/euler/2d/shockbubble/setprob.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.6/qinit.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.6/setaux.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.6/bc2.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.6/src2.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.6/tag4coarsening.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.6/tag4refinement.f \
	applications/amrclaw/euler/2d/shockbubble/user_5.0/qinit.f90 \
	applications/amrclaw/euler/2d/shockbubble/user_5.0/setaux.f90 \
	applications/amrclaw/euler/2d/shockbubble/user_5.0/bc2.f90 \
	applications/amrclaw/euler/2d/shockbubble/user_5.0/src2.f90 \
	applications/amrclaw/euler/2d/shockbubble/user_5.0/tag4coarsening.f90 \
	applications/amrclaw/euler/2d/shockbubble/user_5.0/tag4refinement.f90 \
	applications/amrclaw/euler/2d/rp/clawpack46_rpn2_euler3.f \
	applications/amrclaw/euler/2d/rp/clawpack46_rpn2_euler4.f \
	applications/amrclaw/euler/2d/rp/clawpack46_rpn2_euler5.f \
	applications/amrclaw/euler/2d/rp/clawpack46_rpt2.f \
	applications/amrclaw/euler/2d/rp/clawpack5_rpn2_euler4.f90 \
	applications/amrclaw/euler/2d/rp/clawpack5_rpn2_euler5.f90 \
	applications/amrclaw/euler/2d/rp/clawpack5_rpt2.f90

applications_amrclaw_euler_2d_shockbubble_shockbubble_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_amrclaw_euler_2d_shockbubble_shockbubble_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
