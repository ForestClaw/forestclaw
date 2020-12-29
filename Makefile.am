
# This file is part of ForestClaw.
# Makefile.am in toplevel directory

# include makefiles from installed subpackages
ACLOCAL_AMFLAGS = -I config @FCLAW_SC_AMFLAGS@ @FCLAW_P4EST_AMFLAGS@
if FCLAW_SC_MK_USE
@FCLAW_SC_MK_INCLUDE@
endif
if FCLAW_P4EST_MK_USE
@FCLAW_P4EST_MK_INCLUDE@
endif

# initialize empty variables
AM_CPPFLAGS =
AM_FFLAGS =
CLEANFILES =
DISTCLEANFILES =
EXTRA_DIST =
LDADD = $(FLIBS)
LINT_CSOURCES =
TESTS =
TEST_SOURCES =
TEST_CPPFLAGS =
TEST_LDADD =
bin_PROGRAMS =
check_PROGRAMS =
include_HEADERS =
lib_LTLIBRARIES =
nodist_include_HEADERS =
noinst_HEADERS =
noinst_PROGRAMS =
noinst_LIBRARIES =
sysconf_DATA =

# use this if you want to link in ForestClaw without autotools
sysconf_DATA += Makefile.fclaw.mk
CLEANFILES += Makefile.fclaw.mk
Makefile.fclaw.mk : Makefile.fclaw.pre
	cat $< | \
        sed -e 's,{\(\(.*prefix\|sysconfdir\)\)},{fclaw_\1},g' \
            -e 's,^\(\(.*prefix\|sysconfdir\) *=\),fclaw_\1,g' > $@

# install ForestClaw m4 macros in the correct directory
fclawaclocaldir = $(datadir)/aclocal
dist_fclawaclocal_DATA = config/fclaw_include.m4

# install example .ini files in a dedicated directory
fclawinidir = $(datadir)/ini
dist_fclawini_DATA =
AM_CPPFLAGS += -DFCLAW_INIDIR=\"$(fclawinidir)\"

# setup test environment
if FCLAW_MPIRUN
LOG_COMPILER = @FCLAW_MPIRUN@
AM_LOG_FLAGS = @FCLAW_MPI_TEST_FLAGS@
endif

# recurse only into subpackages
SUBDIRS = @FCLAW_SC_SUBDIR@ @FCLAW_P4EST_SUBDIR@
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += \
  bootstrap fclawindent build-aux/git-version-gen build-aux/git2cl
DISTCLEANFILES += \
        _configs.sed src/fclaw_config.h @FCLAW_DISTCLEAN@
.PHONY: lint ChangeLog

# --------------------- Solvers and examples ---------------------

# --------------------------------------------
# ForestClaw 
# --------------------------------------------
include src/Makefile.am

# Examples
include Makefile.apps

# Tests
include test/Makefile.am

# ------------------------------------------------------------------

# revision control and ChangeLog
ChangeLog:
	(GIT_DIR=@top_srcdir@/.git git log > .ChangeLog.tmp && \
         cat .ChangeLog.tmp | @top_srcdir@/build-aux/git2cl > ChangeLog) ; \
        rm -f .ChangeLog.tmp

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
	test "x$(VERSION)" = "x`@top_srcdir@/build-aux/git-version-gen \
          @top_srcdir@/.tarball-version`" || \
        ((echo "Stale version;" ; echo "Please run:" ; \
          echo "     cd @top_srcdir@ && ./bootstrap" ; \
          echo "before make dist") 1>&2 ; rm -r $(distdir) ; exit 1)
if FCLAW_DIST_DENY
	@echo "-----------------------------------------------------"
	@echo "make dist does not work with external subpackages"
	@echo "-----------------------------------------------------"
	rm -rf $(distdir)
	@exit 1
endif

maintainer-clean-local:
	rm -rf doxygen

clean-local:
	rm -f ChangeLog *vtu *.visit

mostlyclean-local:
	rm -f *.mod



doxygen: Doxyfile
	doxygen

.PHONY: doxygen
