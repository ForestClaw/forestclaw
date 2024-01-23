## Makefile.am in applications/clawpack/advection/3dx/disk

bin_PROGRAMS += applications/clawpack/advection/3dx/disk/disk


applications_clawpack_advection_3dx_disk_disk_SOURCES = \
	applications/clawpack/advection/3dx/disk/disk.cpp \
	applications/clawpack/advection/3dx/disk/disk_options.c \
	applications/clawpack/advection/3dx/disk/disk_user.cpp \
	applications/clawpack/advection/3dx/disk/disk_map_extrude.c \
	applications/clawpack/advection/3dx/disk/psi.f \
	applications/clawpack/advection/3dx/disk/fdisc.f90 \
	applications/clawpack/advection/3dx/disk/setprob.f \
	applications/clawpack/advection/3dx/disk/user_4.6/qinit.f90 \
	applications/clawpack/advection/2d/all/fclaw2d_map_pillowdisk.c \
	applications/clawpack/advection/2d/all/fclaw2d_map_pillowdisk5.c \
	applications/clawpack/advection/3dx/all/claw3_setaux_manifold.f90 \
	applications/clawpack/advection/3dx/all/advection_patch_setup_manifold.c  \
	applications/clawpack/advection/3dx/all/rpn3_vc_advection.f90 \
	applications/clawpack/advection/3dx/all/rpt3_vc_advection.f90 \
	applications/clawpack/advection/3dx/all/rptt3_vc_advection.f90

EXTRA_DIST += \
	applications/clawpack/advection/3dx/disk/fclaw_options.ini \
	applications/clawpack/advection/3dx/disk/regression.ini \
	applications/clawpack/advection/3dx/disk/afterframe.m \
	applications/clawpack/advection/3dx/disk/read_vars.m \
	applications/clawpack/advection/3dx/disk/setplot3.m \
	applications/clawpack/advection/3dx/disk/mapc2p.m \
	applications/clawpack/advection/3dx/disk/filament_soln.m

applications_clawpack_advection_3dx_disk_disk_CPPFLAGS = -I../all/ \
        $(AM_CPPFLAGS) \
        $(FCLAW_CLAWPACK3_46_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

applications_clawpack_advection_3dx_disk_disk_LDADD = \
        $(LDADD) \
        $(FCLAW_CLAWPACK3_46_LDADD)  \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)

TESTS += applications/clawpack/advection/3dx/disk/regressions.sh