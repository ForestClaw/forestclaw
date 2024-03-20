dnl FCLAW_CHECK_HDF5(PREFIX)
dnl Check for the Hdf5 library
dnl
AC_DEFUN([FCLAW_CHECK_HDF5], [

SC_ARG_WITH_PREFIX([hdf5], [specify hdf5 library], [HDF5], [$1])
if test "x$$1_WITH_HDF5" != xno ; then
  export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$$1_WITH_HDF5/lib/pkgconfig
  PKG_CHECK_MODULES([HDF5],[hdf5])
  AC_MSG_RESULT([successful])
fi
])

