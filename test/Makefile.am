## Makefile.am in applications/elliptic/phasefield

check_LTLIBRARIES = test/libtestutils.la

#TEST_EXTENTIONS += .MPI1 .TEST

TEST_EXTENSIONS   = .sh .TEST .MPI1 .txt
SH_LOG_COMPILER   = $(SHELL)
AM_SH_LOG_FLAGS   = 
MPI1_LOG_DRIVER   = @top_srcdir@/test/doc-test-driver
AM_MPI1_LOG_FLAGS =
TEST_LOG_DRIVER   = @top_srcdir@/test/doc-test-driver
TXT_LOG_DRIVER    = @top_srcdir@/test/regression-driver
AM_TEST_LOG_FLAGS =
LOG_COMPILER = $(SHELL)
AM_LOG_FLAGS = 

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