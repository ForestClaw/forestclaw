## Makefile.am in applications/NCAR_transport/correlatedcb

bin_PROGRAMS += applications/NCAR_transport/correlatedcb/correlatedcb

applications_NCAR_transport_correlatedcb_correlatedcb_SOURCES = \
	applications/NCAR_transport/correlatedcb/correlatedcb_user.cpp \
	applications/NCAR_transport/correlatedcb/correlatedcb.cpp \
	applications/NCAR_transport/correlatedcb/correlatedcb_output.f \
	applications/NCAR_transport/correlatedcb/tag4refinement.f \
	applications/NCAR_transport/correlatedcb/fclaw2d_map_cubedsphere.c \
	applications/NCAR_transport/correlatedcb/fclaw2d_map_pillowsphere.c \
	applications/NCAR_transport/all/setprob_transport.f \
	applications/NCAR_transport/all/qinit_transport.f \
	applications/NCAR_transport/all/slotted_disk.f \
	applications/NCAR_transport/all/gaussian.f \
	applications/NCAR_transport/all/cosine_bell.f \
	applications/NCAR_transport/all/velocity_wind.f \
	applications/NCAR_transport/all/setaux_transport.f \
	applications/NCAR_transport/all/rpn2noncons.f \
	applications/NCAR_transport/all/rpt2noncons.f \
	applications/NCAR_transport/all/exact_solution.f \
	applications/NCAR_transport/all/b4step2_transport.f

applications_NCAR_transport_correlatedcb_correlatedcb_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_NCAR_transport_correlatedcb_correlatedcb_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
