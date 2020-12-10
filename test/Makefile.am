## Makefile.am in applications/elliptic/phasefield

check_PROGRAMS += test/test

TESTS = test/run_test.sh

test_test_SOURCES = \
	test/test.cpp \
	test/utils/DomainReader.h \
	test/utils/DomainReader.cpp \
	$(TEST_SOURCES)

test_test_CPPFLAGS = $(TEST_CPPFLAGS) -DDATADIR=\"$(srcdir)/\"
test_test_LDADD = $(TEST_LDADD)
