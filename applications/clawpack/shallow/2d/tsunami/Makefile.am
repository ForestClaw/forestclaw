bin_PROGRAMS += applications/clawpack/shallow/2d/tsunami/tsunami

applications_clawpack_shallow_2d_tsunami_tsunami_SOURCES = \
	applications/clawpack/shallow/2d/tsunami/tsunami.cpp \
	applications/clawpack/shallow/2d/tsunami/tsunami_options.c \
	applications/clawpack/shallow/2d/tsunami/tsunami_user.cpp \
	applications/clawpack/shallow/2d/tsunami/tsunami_user.h \
	applications/clawpack/shallow/2d/tsunami/tsunami_run.c \
	applications/clawpack/shallow/2d/tsunami/sgn_operator.cpp \
	applications/clawpack/shallow/2d/tsunami/sgn_patch_operator.cpp \
	applications/clawpack/shallow/2d/tsunami/setprob.f \
	applications/clawpack/shallow/2d/tsunami/sgn1d/sgn.f90 \
	applications/clawpack/shallow/2d/tsunami/sgn1d/tri_solve.f90 \
	applications/clawpack/shallow/2d/tsunami/user_4.6/qinit.f90 \
	applications/clawpack/shallow/2d/tsunami/user_4.6/setaux.f90 \
	applications/clawpack/shallow/2d/tsunami/user_4.6/tag4refinement.f90 \
	applications/clawpack/shallow/2d/tsunami/user_4.6/tag4coarsening.f90 \
	applications/clawpack/shallow/2d/tsunami/user_4.6/sgn_apply_bc.f90 \
	applications/clawpack/shallow/2d/tsunami/user_4.6/sgn_eval_bc.f90 \
	applications/clawpack/shallow/2d/tsunami/user_4.6/sgn_fort_rhs.f90 \
	applications/clawpack/shallow/2d/tsunami/user_4.6/sgn_update_q.f90 \
	applications/clawpack/shallow/2d/tsunami/bathy.f90 \
	applications/clawpack/shallow/2d/rp/swe_simple_riemann.f90 \
	applications/clawpack/shallow/2d/rp/rpn2_tsunami.f90 \
	applications/clawpack/shallow/2d/rp/rpt2_tsunami.f90 \
	applications/clawpack/shallow/2d/rp/rpn2_geoclaw.f90 \
	applications/clawpack/shallow/2d/rp/rpt2_geoclaw.f90 \
	applications/clawpack/shallow/2d/rp/geoclaw_riemann_utils.f90

EXTRA_DIST += \
	applications/clawpack/shallow/2d/tsunami/fclaw_options.ini \
	applications/clawpack/shallow/2d/tsunami/regression.ini \
	applications/clawpack/shallow/2d/tsunami/afterframe.m \
	applications/clawpack/shallow/2d/tsunami/setplot2.m

applications_clawpack_shallow_2d_tsunami_tsunami_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_shallow_2d_tsunami_tsunami_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
