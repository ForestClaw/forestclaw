## Makefile.am in applications/advection/2d/sphere

bin_PROGRAMS += applications/advection/2d/torthem/torthem \
		applications/advection/2d/torthem/torthem

applications_advection_2d_torthem_torthem_SOURCES = \
	applications/advection/2d/torthem/torthem.c \
	applications/advection/2d/torthem/torus_common.h \
	applications/advection/2d/torthem/torus_common.c

applications_advection_2d_torthem_torthem_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_advection_2d_torthem_torthem_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) $(FCLAW_CLAWPACK5_LDADD) $(LDADD)

## Due to C++ libraries we need to force C++ linking for a C main program
nodist_EXTRA_applications_advection_2d_torthem_torthem_SOURCES = dummy.cxx
