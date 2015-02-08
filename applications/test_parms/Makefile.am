## Makefile.am in applications/test_parms

# A single hash comment like this one WILL show up in the Makefile.
## A double hash comment like this one will NOT be copied to the Makefile.

bin_PROGRAMS += applications/test_parms/test_parms

applications_test_parms_test_parms_SOURCES = \
	applications/test_parms/test_parms.cpp

applications_test_parms_test_parms_CPPFLAGS = \
       $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS)

applications_test_parms_test_parms_LDADD = \
	$(LDADD) $(FCLAW_CLAWPACK46_LDADD)
