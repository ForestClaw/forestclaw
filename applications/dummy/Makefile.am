## Makefile.am in applications/dummy

bin_PROGRAMS += applications/dummy/fclaw2d_dummy

applications_dummy_fclaw2d_dummy_SOURCES = \
	applications/dummy/fclaw2d_dummy.c

# TODO: we don't need this, but for testing activating this line should link
# DAC : Still need to replace "link by file name" to an all pointer approach.
# Until I fix that, we can't link clawpack solvers, since they will look for a
# bunch of files that need to be present.
# applications_dummy_fclaw2d_dummy_LDADD = $(FCLAW_CLAWPACK46_LDADD) $(LDADD)
