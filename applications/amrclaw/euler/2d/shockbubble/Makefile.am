## Makefile.am in applications/amrclaw/euler/2d/shockbubble

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/amrclaw/euler/2d/shockbubble/shockbubble

applications_amrclaw_euler_2d_shockbubble_shockbubble_SOURCES = \
	applications/amrclaw/euler/2d/shockbubble/shockbubble.cpp \
	applications/amrclaw/euler/2d/shockbubble/shockbubble_user.cpp \
	applications/amrclaw/euler/2d/shockbubble/fdisc.f \
	applications/amrclaw/euler/2d/shockbubble/setprob.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.x/qinit.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.x/setaux.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.x/bc2.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.x/src2.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.x/tag4coarsening.f \
	applications/amrclaw/euler/2d/shockbubble/user_4.x/tag4refinement.f \
	applications/amrclaw/euler/2d/rp/rpn2eu3.f \
	applications/amrclaw/euler/2d/rp/rpn2eu4.f \
	applications/amrclaw/euler/2d/rp/rpn2eu5.f \
	applications/amrclaw/euler/2d/rp/rpt2eu.f \
	applications/amrclaw/euler/2d/rp/rpt2eu5.f

applications_amrclaw_euler_2d_shockbubble_shockbubble_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_amrclaw_euler_2d_shockbubble_shockbubble_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
