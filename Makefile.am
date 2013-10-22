## Toplevel Makefile.am for ForestClaw

ACLOCAL_AMFLAGS = -I config -I $(top_srcdir)@FCLAW_TOPLEVEL_MOD@/sc/config
AM_CPPFLAGS = -I@top_srcdir@/src \
 -I$(top_builddir)@FCLAW_TOPLEVEL_MOD@/sc/src \
 -I$(top_srcdir)@FCLAW_TOPLEVEL_MOD@/sc/src
DISTCLEANFILES = _configs.sed src/fclaw_config.h
EXTRA_DIST = bootstrap fclawindent config/fclaw_include.m4 config/git-version-gen
LDADD =
TESTS =
bin_PROGRAMS =
check_PROGRAMS =
include_HEADERS =
lib_LTLIBRARIES =
nodist_include_HEADERS =

LDADD += $(FCLAW_LDADD) $(top_builddir)@FCLAW_TOPLEVEL_MOD@/sc/src/libsc.la
SUBDIRS = $(FCLAW_SUBDIRS)

## install m4 macros in the correct directory
fclawaclocaldir = $(datadir)/aclocal
dist_fclawaclocal_DATA = config/fclaw_include.m4

## setup test environment
if FCLAW_MPIRUN
LOG_COMPILER = @FCLAW_MPIRUN@
AM_LOG_FLAGS = @FCLAW_MPI_TEST_FLAGS@
endif

## have make dist enforce a clean source tree
dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
	test "x$(VERSION)" = "x`@top_srcdir@/config/git-version-gen \
		@top_srcdir@/.tarball-version`" || \
	((echo "Stale version;" ; echo "Please run:" ; \
	  echo "     cd @top_srcdir@ && ./bootstrap" ; \
	  echo "before make dist") 1>&2 ; exit 1)

## using nonrecursive make
include src/Makefile.am
include applications/dummy/Makefile.am
