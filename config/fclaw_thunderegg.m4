dnl FCLAW_CHECK_THUNDEREGG(PREFIX)
dnl Check for the Thunderegg library
dnl
AC_DEFUN([FCLAW_CHECK_THUNDEREGG], [

SC_ARG_WITH_PREFIX([thunderegg], [specify thunderegg library], [TE], [$1])
if test "x$$1_WITH_TE" != xno ; then
  export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$$1_WITH_TE/lib/pkgconfig
  PKG_CHECK_MODULES([THUNDEREGG],[ThunderEgg])
  AC_MSG_RESULT([successful])
fi
])

