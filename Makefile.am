## Toplevel Makefile.am for forestclaw

ACLOCAL_AMFLAGS = -I config -I sc/config
AM_CPPFLAGS =
DISTCLEANFILES =
EXTRA_DIST =
lib_LTLIBRARIES =

DISTCLEANFILES += _configs.sed src/fclaw_config.h

EXTRA_DIST += bootstrap config/fclaw_include.m4 config/git-version-gen

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
	test "x$(VERSION)" = "x`@top_srcdir@/config/git-version-gen \
		@top_srcdir@/.tarball-version`" || \
	((echo "Stale version;" ; echo "Please run:" ; \
	  echo "     cd @top_srcdir@ && ./bootstrap" ; \
	  echo "before make dist") 1>&2 ; exit 1)

## using nonrecursive make
include src/Makefile.am
