dnl Copyright 2015 Lorenz Hüdepohl
dnl
dnl This file is part of fdep and licensed under the MIT license
dnl see the file LICENSE for more information
dnl
AC_DEFUN([FDEP_F90_GNU_MAKE_DEPS],[
AC_MSG_CHECKING([for GNU make])
for a in "$MAKE" make gmake gnumake ; do
        if test -z "$a" ; then continue ; fi ;
        if  ( sh -c "$a --version" 2> /dev/null | grep GNU  2>&1 > /dev/null ) ;  then
                _fdep_gnu_make_command=$a ;
                break;
        fi
done ;
AC_MSG_RESULT([$_fdep_gnu_make_command])
if test x$_fdep_gnu_make_command = x ; then
	AC_MSG_ERROR([Need GNU Make])
fi
AC_SUBST([FORTRAN_MODULE_DEPS], ["
CLEANFILES +=
include ${srcdir}/fdep/fortran_dependencies.mk
"])
AM_SUBST_NOTMAKE([FORTRAN_MODULE_DEPS])
])
