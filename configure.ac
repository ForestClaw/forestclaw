dnl
dnl This file is part of ForestClaw.
dnl

dnl We need a forwarding email for both Donna and Carsten.
AC_INIT([forestclaw],
        [m4_esyscmd([config/git-version-gen .tarball-version])],
        [burstedde@ins.uni-bonn.de])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([src/pre_config.h])
AC_CONFIG_SRCDIR([src/fclaw2d_base.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([config])
AC_PREFIX_DEFAULT([$PWD/local])
AX_PREFIX_CONFIG_H([src/fclaw_config.h], [FCLAW])
AM_INIT_AUTOMAKE([subdir-objects])
AM_SILENT_RULES
SC_VERSION([FCLAW])

# Check for subpackage sc
echo "o---------------------------------------"
echo "| Checking for subpackages"
echo "o---------------------------------------"

AC_MSG_NOTICE([Current srcdir: $srcdir])
AC_MSG_CHECKING([whether we are a toplevel package])
if test -e "$srcdir/sc/src/sc.h" ; then
AC_MSG_RESULT([yes])
FCLAW_IS_TOPLEVEL=true
FCLAW_TOPLEVEL_MOD=
FCLAW_SUBDIRS=sc
AC_DEFINE([IS_TOPLEVEL], [1], [We are a toplevel package])
AC_CONFIG_SUBDIRS([sc])
else
if test -e "$srcdir/../sc/src/sc.h" ; then
AC_MSG_RESULT([no])
FCLAW_IS_TOPLEVEL=
FCLAW_TOPLEVEL_MOD=/..
FCLAW_SUBDIRS=
else
AC_MSG_ERROR([Could not find package sc])
fi
fi
AM_CONDITIONAL([IS_TOPLEVEL], [test x$FCLAW_IS_TOPLEVEL = xtrue])
dnl In Makefile.am files we can use
dnl $(top_srcdir)@FCLAW_TOPLEVEL_MOD@/path and it should work
AC_SUBST([FCLAW_TOPLEVEL_MOD])
AC_SUBST([FCLAW_SUBDIRS])

# Checks for MPI C compiler and other programs (no F77 at this point!).
echo "o---------------------------------------"
echo "| Checking MPI and other programs"
echo "o---------------------------------------"

SC_MPI_CONFIG([FCLAW])
SC_MPI_ENGAGE([FCLAW])
# This is needed for compatibility with automake >= 1.12
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
AM_PROG_CC_C_O
LT_INIT

# Checks for libraries.
echo "o---------------------------------------"
echo "| Checking libraries"
echo "o---------------------------------------"

SC_REQUIRE_LIB([m], [fabs])

dnl dnl dnl # Checks for header files.
dnl dnl dnl echo "o---------------------------------------"
dnl dnl dnl echo "| Checking headers"
dnl dnl dnl echo "o---------------------------------------"
dnl dnl dnl 
dnl dnl dnl AC_CHECK_HEADERS([execinfo.h signal.h sys/time.h sys/types.h time.h])

# Checks for functions.
echo "o---------------------------------------"
echo "| Checking functions"
echo "o---------------------------------------"

dnl dnl dnl AC_CHECK_FUNCS([backtrace backtrace_symbols])

dnl dnl dnl AC_PROG_F77([gfortran g77 fl32 f77 fort77 xlf f90 xlf90])
dnl dnl dnl AC_F77_WRAPPERS
dnl dnl dnl AC_F77_LIBRARY_LDFLAGS
dnl dnl dnl dgemm=;AC_F77_FUNC(dgemm)

# Checks for typedefs, structures, and compiler characteristics.
echo "o---------------------------------------"
echo "| Checking keywords and types"
echo "o---------------------------------------"

AC_C_BIGENDIAN
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([unsigned long])
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

AC_DEFINE_UNQUOTED(CC,          ["${CC}"],          [C compiler])
AC_DEFINE_UNQUOTED(CFLAGS,      ["${CFLAGS}"],      [C compiler flags])
AC_DEFINE_UNQUOTED(CPP,         ["${CPP}"],         [C preprocessor])
AC_DEFINE_UNQUOTED(CPPFLAGS,    ["${CPPFLAGS}"],    [C preprocessor flags])
AC_DEFINE_UNQUOTED(F77,         ["${F77}"],         [F77 compiler])
AC_DEFINE_UNQUOTED(FFLAGS,      ["${FFLAGS}"],      [F77 compiler flags])
AC_DEFINE_UNQUOTED(LDFLAGS,     ["${LDFLAGS}"],     [Linker flags])
AC_DEFINE_UNQUOTED(BLAS_LIBS,   ["${BLAS_LIBS}"],   [BLAS libraries])
AC_DEFINE_UNQUOTED(LAPACK_LIBS, ["${LAPACK_LIBS}"], [LAPACK libraries])
AC_DEFINE_UNQUOTED(LIBS,        ["${LIBS}"],        [Libraries])
AC_DEFINE_UNQUOTED(FLIBS,       ["${FLIBS}"],       [Fortran libraries])

echo "o----------------------------------"
echo "| Results for libsc are"
echo "o----------------------------------"
echo "| CC:          $CC"
echo "| CFLAGS:      $CFLAGS"
echo "| CPP:         $CPP"
echo "| CPPFLAGS:    $CPPFLAGS"
echo "| F77:         $F77"
echo "| FFLAGS:      $FFLAGS"
echo "| LDFLAGS:     $LDFLAGS"
echo "| BLAS_LIBS:   $BLAS_LIBS"
echo "| LAPACK_LIBS: $LAPACK_LIBS"
echo "| LIBS:        $LIBS"
echo "| FLIBS:       $FLIBS"
echo "o----------------------------------"

# Create output files.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
