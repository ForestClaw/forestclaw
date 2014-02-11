
# This file is part of ForestClaw.
# Makefile.am in src
# included non-recursively from toplevel directory

libforestclaw_generated_headers = src/fclaw_config.h
libforestclaw_installed_headers = src/fclaw2d_base.h
libforestclaw_internal_headres =
libforestclaw_compiled_sources = src/fclaw2d_base.c

# this variable is used for headers that are not publicly installed
FCLAW_CPPFLAGS =

lib_LTLIBRARIES += src/libforestclaw.la
src_libforestclaw_la_SOURCES = \
 $(libforestclaw_internal_headers) \
 $(libforestclaw_compiled_sources)
src_libforestclaw_la_CPPFLAGS = $(AM_CPPFLAGS) $(FCLAW_CPPFLAGS)
LDADD += @top_builddir@/src/libforestclaw.la @FCLAW_P4EST_LDADD@ @FCLAW_SC_LDADD@
EXTRA_src_libforestclaw_la_DEPENDENCIES = @FCLAW_P4EST_LDADD@ @FCLAW_SC_LDADD@

nodist_include_HEADERS += $(libforestclaw_generated_headers)
include_HEADERS += $(libforestclaw_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src @FCLAW_SC_CPPFLAGS@ @FCLAW_P4EST_CPPFLAGS@
