## Makefile.am in applications/elliptic/phasefield

check_PROGRAMS += test/test

TESTS = test/run_test.sh

test_test_SOURCES = \
	test/test.cpp \
	$(TEST_SOURCES)

if FCLAW_ENABLE_THUNDEREGG
test_test_SOURCES += \
	test/utils/DomainReader.h \
	test/utils/DomainReader.cpp 
endif

test_test_CPPFLAGS = -DDATADIR=\"$(srcdir)/\" \
        $(AM_CPPFLAGS) \
        $(FCLAW_THUNDEREGG_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS) \
        $(FCLAW_CPPFLAGS)
test_test_LDADD = \
        $(LDADD) \
        $(FCLAW_THUNDEREGG_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)