bin_PROGRAMS += applications/clawpack/shallow/2d/wavetank/wavetank

applications_clawpack_shallow_2d_wavetank_wavetank_SOURCES = \
	applications/clawpack/shallow/2d/wavetank/wavetank.cpp \
	applications/clawpack/shallow/2d/wavetank/wavetank_options.c \
	applications/clawpack/shallow/2d/wavetank/wavetank_user.cpp \
	applications/clawpack/shallow/2d/wavetank/wavetank_user.h \
	applications/clawpack/shallow/2d/wavetank/setprob.f \
	applications/clawpack/shallow/2d/wavetank/bathy.f90 \
	applications/clawpack/shallow/2d/wavetank/user_4.6/qinit.f90 \
	applications/clawpack/shallow/2d/wavetank/user_4.6/setaux.f90 \
	applications/clawpack/shallow/2d/wavetank/user_4.6/tag4refinement.f90 \
	applications/clawpack/shallow/2d/wavetank/user_4.6/tag4coarsening.f90 \
	applications/clawpack/shallow/2d/sgn/sgn_operator.cpp \
	applications/clawpack/shallow/2d/sgn/sgn_patch_operator.cpp \
	applications/clawpack/shallow/2d/sgn/sgn_options.c \
	applications/clawpack/shallow/2d/sgn/sgn_run.c \
	applications/clawpack/shallow/2d/sgn/sgn.c \
	applications/clawpack/shallow/2d/sgn/fort_4.6/sgn_apply_bc.f90 \
	applications/clawpack/shallow/2d/sgn/fort_4.6/sgn_eval_bc.f90 \
	applications/clawpack/shallow/2d/sgn/fort_4.6/sgn_fort_rhs.f90 \
	applications/clawpack/shallow/2d/sgn/fort_4.6/sgn_update_q.f90 \
	applications/clawpack/shallow/2d/rp/swe_simple_riemann.f90 \
	applications/clawpack/shallow/2d/rp/rpn2_tsunami.f90 \
	applications/clawpack/shallow/2d/rp/rpt2_tsunami.f90 \
	applications/clawpack/shallow/2d/rp/rpn2_geoclaw.f90 \
	applications/clawpack/shallow/2d/rp/rpt2_geoclaw.f90 \
	applications/clawpack/shallow/2d/rp/geoclaw_riemann_utils.f90

EXTRA_DIST += \
	applications/clawpack/shallow/2d/wavetank/fclaw_options.ini \
	applications/clawpack/shallow/2d/wavetank/regression.ini \
	applications/clawpack/shallow/2d/wavetank/afterframe.m \
	applications/clawpack/shallow/2d/wavetank/setplot2.m

applications_clawpack_shallow_2d_wavetank_wavetank_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_shallow_2d_wavetank_wavetank_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
