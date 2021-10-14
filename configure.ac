dnl
dnl This file is part of ForestClaw.
dnl

dnl TODO: We need a forwarding email for both Donna and Carsten.
AC_INIT([ForestClaw],
        [m4_esyscmd([build-aux/git-version-gen .tarball-version])],
        [burstedde@ins.uni-bonn.de])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([src/pre_config.h])
AC_CONFIG_SRCDIR([src/fclaw_base.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([config])
AC_PREFIX_DEFAULT([$PWD/local])
AX_PREFIX_CONFIG_H([src/fclaw_config.h], [FCLAW])
AC_CONFIG_FILES([test/test_config.h])
AM_INIT_AUTOMAKE([parallel-tests subdir-objects])
AM_SILENT_RULES
SC_VERSION([FCLAW])

m4_include([fdep/fortran_dependencies.m4])
FDEP_F90_GNU_MAKE_DEPS

echo "o---------------------------------------"
echo "| Checking options"
echo "o---------------------------------------"

FCLAW_ARG_ENABLE([debug], [Enable debug mode (assertions and extra checks)],
                 [DEBUG])
                 
FCLAW_ARG_ENABLE([geoclaw], [Enable geoclaw library and applications],
                 [GEOCLAW])

FCLAW_ARG_ENABLE([ash3d], [enable ash3d library and applications],
                 [ASH3D])

FCLAW_ARG_ENABLE([clawpack], [Enable clawpack libraries 4.6 and 5.0 and applications],
                 [CLAWPACK])

FCLAW_ARG_ENABLE([clawpack3], [Enable clawpack 3d libraries and applications],
                 [CLAWPACK3])

FCLAW_ARG_ENABLE([cudaclaw], [Enable cudaclaw library and application],
                 [CUDACLAW])

FCLAW_ARG_ENABLE([development], [Enable general applications under development],
                 [DEVELOPMENT])

FCLAW_ARG_ENABLE([paper], [Enable examples applications from paper],
                 [PAPER])

FCLAW_ARG_ENABLE([activeflux], [Enable activeflux project],
                 [ACTIVEFLUX])

FCLAW_ARG_ENABLE([thunderegg], [Enable thunderegg solver],
                 [THUNDEREGG])

FCLAW_ARG_ENABLE([hps], [Enable hps solver],
                 [HPS])

echo "o---------------------------------------"
echo "| Checking MPI and related programs"
echo "o---------------------------------------"

SC_MPI_CONFIG([FCLAW], [yes], [yes])
SC_MPI_ENGAGE([FCLAW])
# This is needed for compatibility with automake >= 1.12
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
SC_PROG_LINT
LT_INIT

echo "o---------------------------------------"
echo "| Checking libraries"
echo "o---------------------------------------"

SC_CHECK_LIBRARIES([FCLAW])
P4EST_CHECK_LIBRARIES([FCLAW])
FCLAW_CHECK_LIBRARIES([FCLAW])

echo "o---------------------------------------"
echo "| Checking headers"
echo "o---------------------------------------"s

AC_CHECK_HEADERS([fenv.h signal.h])

echo "o---------------------------------------"
echo "| Checking functions"
echo "o---------------------------------------"

AC_CHECK_FUNCS([feenableexcept])

echo "o---------------------------------------"
echo "| Checking subpackages"
echo "o---------------------------------------"

FCLAW_DIST_DENY=
FCLAW_DISTCLEAN=
SC_AS_SUBPACKAGE([FCLAW])
export P4EST_SC_SOURCE=../sc
P4EST_AS_SUBPACKAGE([FCLAW])
AM_CONDITIONAL([FCLAW_DIST_DENY], [test "x$FCLAW_DIST_DENY" != x])
AC_SUBST([FCLAW_DISTCLEAN])

AX_CHECK_CUDA()

# Print summary.
AC_DEFINE_UNQUOTED(CPP,         ["${CPP}"],         [C preprocessor])
AC_DEFINE_UNQUOTED(CPPFLAGS,    ["${CPPFLAGS}"],    [C preprocessor flags])
AC_DEFINE_UNQUOTED(F77,         ["${F77}"],         [F77 compiler])
AC_DEFINE_UNQUOTED(FFLAGS,      ["${FFLAGS}"],      [F77 compiler flags])
AC_DEFINE_UNQUOTED(FCFLAGS,     ["${FCFLAGS}"],     [FC compiler flags])
AC_DEFINE_UNQUOTED(CC,          ["${CC}"],          [C compiler])
AC_DEFINE_UNQUOTED(CFLAGS,      ["${CFLAGS}"],      [C compiler flags])
AC_DEFINE_UNQUOTED(CXX,         ["${CXX}"],         [C++ compiler])
AC_DEFINE_UNQUOTED(CXXFLAGS,    ["${CXXFLAGS}"],    [C++ compiler flags])
AC_DEFINE_UNQUOTED(NVCC,        ["${NVCC}"],        [Nvidia C/C++ compiler])
AC_DEFINE_UNQUOTED(CUDA_CFLAGS, ["${CUDA_CFLAGS}"], [Nvidia C/C++ flags])
AC_DEFINE_UNQUOTED(CUDA_LDFLAGS,["${CUDA_LDFLAGS}"],[Nvidia LD flags])
AC_DEFINE_UNQUOTED(LDFLAGS,     ["${LDFLAGS}"],     [Linker flags])
AC_DEFINE_UNQUOTED(LIBS,        ["${LIBS}"],        [Libraries])
AC_DEFINE_UNQUOTED(FLIBS,       ["${FLIBS}"],       [Fortran libraries])

echo "o----------------------------------"
echo "| Results for ForestClaw are"
echo "o----------------------------------"
echo "| CPP:             $CPP"
echo "| CPPFLAGS:        $CPPFLAGS"
echo "| F77:             $F77"
echo "| FFLAGS:          $FFLAGS"
echo "| FCFLAGS:         $FCFLAGS"
echo "| CC:              $CC"
echo "| CFLAGS:          $CFLAGS"
echo "| CXX:             $CXX"
echo "| CXXFLAGS:        $CXXFLAGS"
echo "| NVCC:            $NVCC"
echo "| CUDA_CFLAGS:     $CUDA_CFLAGS"
echo "| CUDA_LDFLAGS:    $CUDA_LDFLAGS"
echo "| LDFLAGS:         $LDFLAGS"
echo "| LIBS:            $LIBS"
echo "| FLIBS:           $FLIBS"
echo "o----------------------------------"

# Create output files.
AC_CONFIG_FILES([Makefile Makefile.fclaw.pre Doxyfile])
AC_OUTPUT

# Final messages.
SC_FINAL_MESSAGES([FCLAW])
P4EST_FINAL_MESSAGES([FCLAW])
FCLAW_FINAL_MESSAGES([FCLAW])
