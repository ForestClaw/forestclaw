## Makefile.am in src

libforestclaw_generated_headers = src/fclaw_config.h
libforestclaw_installed_headers = src/fclaw2d_base.h
libforestclaw_compiled_sources = src/fclaw2d_base.c

include_HEADERS += $(libforestclaw_installed_headers)
nodist_include_HEADERS += $(libforestclaw_generated_headers)

lib_LTLIBRARIES += src/libforestclaw.la
src_libforestclaw_la_SOURCES = \
 $(libforestclaw_installed_headers) $(libforestclaw_compiled_sources)

FCLAW_LDADD = @top_builddir@/src/libforestclaw.la
