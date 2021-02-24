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

test_test_CPPFLAGS = $(TEST_CPPFLAGS) $(AM_CPPFLAGS) -DDATADIR=\"$(srcdir)/\"
test_test_LDADD = $(TEST_LDADD) $(AM_LDADD)
