## Toplevel Makefile.am for forestclaw

ACLOCAL_AMFLAGS = -I config -I $(top_srcdir)@FCLAW_TOPLEVEL_MOD@/sc/config
AM_CPPFLAGS = -I@top_srcdir@/src \
 -I$(top_builddir)@FCLAW_TOPLEVEL_MOD@/sc/src \
 -I$(top_srcdir)@FCLAW_TOPLEVEL_MOD@/sc/src
DISTCLEANFILES = _configs.sed src/fclaw_config.h
EXTRA_DIST = bootstrap fclawindent config/fclaw_include.m4 config/git-version-gen
LDADD =
bin_PROGRAMS =
include_HEADERS =
lib_LTLIBRARIES =
nodist_include_HEADERS =

LDADD += $(FCLAW_LDADD) $(top_builddir)@FCLAW_TOPLEVEL_MOD@/sc/src/libsc.a
SUBDIRS = $(FCLAW_SUBDIRS)

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
