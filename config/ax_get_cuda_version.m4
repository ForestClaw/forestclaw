# ===========================================================================
#     https://www.gnu.org/software/autoconf-archive/ax_split_version.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_GET_CUDA_VERSION
#
# DESCRIPTION
#
#   Gets the CUDA version.
#   Sets the variables.
#
# LICENSE
#
#   Copyright (c) 2008 Tom Howard <tomhoward@users.sf.net> AX_SPLIT_VERSION
#   Copyright (c) 2022 Scott Aiton - modified to check for cuda version
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 10

AC_DEFUN([AX_GET_CUDA_VERSION],[
    AC_REQUIRE([AC_PROG_SED])
    AC_REQUIRE([AC_PROG_GREP])
    CUDA_VERSION=`$NVCC --version | $GREP release | $SED 's/^.*V//'`
    CUDA_MAJOR_VERSION=`echo "$CUDA_VERSION" | $SED 's/\([[^.]][[^.]]*\).*/\1/'`
    CUDA_MINOR_VERSION=`echo "$CUDA_VERSION" | $SED 's/[[^.]][[^.]]*.\([[^.]][[^.]]*\).*/\1/'`
    AC_MSG_CHECKING([CUDA version])
    AC_MSG_RESULT([$CUDA_VERSION])
    AC_MSG_CHECKING([CUDA Major version])
    AC_MSG_RESULT([$CUDA_MAJOR_VERSION])
    AC_MSG_CHECKING([CUDA Minor version])
    AC_MSG_RESULT([$CUDA_MINOR_VERSION])
])
