## Makefile.am in applications/elliptic/phasefield

check_LTLIBRARIES = test/libtestutils.la

TEST_EXTENTIONS = .MPI1 .TEST
MPI1_LOG_COMPILER = @top_srcdir@/test/run_test.sh

TEST_SRC_DIR = @top_srcdir@

test_libtestutils_la_SOURCES = \
	test/test.cpp

if FCLAW_ENABLE_THUNDEREGG
test_libtestutils_la_SOURCES += \
	test/utils/DomainReader.h \
	test/utils/DomainReader.cpp 
endif

test_libtestutils_la_CPPFLAGS = \
	-I@top_srcdir@/test \
	$(AM_CPPFLAGS) \
	$(FCLAW_CPPFLAGS)

if FCLAW_ENABLE_THUNDEREGG
test_libtestutils_la_CPPFLAGS += \
	$(FCLAW_THUNDEREGG_CPPFLAGS)
endif