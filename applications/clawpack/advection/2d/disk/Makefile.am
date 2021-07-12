## Makefile.am in applications/clawpack/advection/2d/disk

bin_PROGRAMS += applications/clawpack/advection/2d/disk/disk


applications_clawpack_advection_2d_disk_disk_SOURCES = \
	applications/clawpack/advection/2d/disk/disk.cpp \
	applications/clawpack/advection/2d/disk/disk_options.c \
	applications/clawpack/advection/2d/disk/disk_user.cpp \
	applications/clawpack/advection/2d/disk/disk_user.h \
	applications/clawpack/advection/2d/disk/psi.f \
	applications/clawpack/advection/2d/disk/fdisc.f \
	applications/clawpack/advection/2d/disk/setprob.f \
	applications/clawpack/advection/2d/disk/user_4.6/qinit.f \
	applications/clawpack/advection/2d/disk/user_5.0/qinit.f90 \
	applications/clawpack/advection/2d/all/advection_patch_setup_manifold.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_pillowdisk.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_pillowdisk5.c \
	applications/clawpack/advection/2d/all/user46_setaux_manifold.f \
	applications/clawpack/advection/2d/all/user5_setaux_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/clawpack/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/clawpack/advection/2d/all/clawpack5_rpt2adv_manifold.f90

EXTRA_DIST += \
	applications/clawpack/advection/2d/disk/fclaw_options.ini \
	applications/clawpack/advection/2d/disk/regression.ini \
	applications/clawpack/advection/2d/disk/afterframe.m \
	applications/clawpack/advection/2d/disk/setplot2.m \
	applications/clawpack/advection/2d/disk/mapc2m.m \
	applications/clawpack/advection/2d/disk/filament_soln.m \
	applications/clawpack/advection/2d/disk/create_filament_soln_tikz.m

applications_clawpack_advection_2d_disk_disk_CPPFLAGS = -I../all/ \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK46_CPPFLAGS) \
        $(FCLAW_CLAWPACK5_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_advection_2d_disk_disk_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK46_LDADD)  \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)

TESTS += applications/clawpack/advection/2d/disk/regressions.sh