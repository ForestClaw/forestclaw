## Makefile.am in applications/advection/2d/disk

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/advection/2d/disk/disk


applications_advection_2d_disk_disk_SOURCES = \
	applications/advection/2d/disk/disk_user.cpp \
	applications/advection/2d/disk/disk_user.h \
	applications/advection/2d/disk/disk.cpp \
	applications/advection/2d/disk/fclaw2d_map_pillowdisk.c \
	applications/advection/2d/disk/fclaw2d_map_pillowdisk5.c \
	applications/advection/2d/disk/psi.f \
	applications/advection/2d/disk/fdisc.f \
	applications/advection/2d/disk/setprob.f \
	applications/advection/2d/disk/user_4.6/qinit.f \
	applications/advection/2d/disk/user_4.6/tag4refinement.f \
	applications/advection/2d/disk/user_4.6/tag4coarsening.f \
	applications/advection/2d/disk/user_5.0/qinit.f90 \
	applications/advection/2d/disk/user_5.0/tag4refinement.f \
	applications/advection/2d/disk/user_5.0/tag4coarsening.f \
	applications/advection/2d/all/user46_setaux_manifold.f \
	applications/advection/2d/all/clawpack46_rpn2adv_manifold.f \
	applications/advection/2d/all/clawpack46_rpt2adv_manifold.f \
	applications/advection/2d/all/user5_setaux_manifold.f90 \
	applications/advection/2d/all/clawpack5_rpn2adv_manifold.f90 \
	applications/advection/2d/all/clawpack5_rpt2adv_manifold.f90


applications_advection_2d_disk_disk_CPPFLAGS = -I../all/ \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

# applications_advection_2d_disk_disk_CPPFLAGS += -D MPI_DEBUG


applications_advection_2d_disk_disk_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD)  $(FCLAW_CLAWPACK5_LDADD) $(LDADD)
