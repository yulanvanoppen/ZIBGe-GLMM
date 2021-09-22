dnl Process this file with Autoconf to produce a configure script/

AC_PREREQ([2.68])

AC_INIT([JAGS-ZIBG], [1.0], [yulanvanoppen@gmail.com], [JAGS-ZIBG-MODULE])
JAGS_MAJOR=4
JAGS_MINOR=3
AC_SUBST(JAGS_MAJOR)
AC_SUBST(JAGS_MINOR)

AC_CANONICAL_HOST
dnl The following lines check if the required files exist
dnl (more AC_CONFIG_SRCDIR directives may be needed - for all relevant files)
AC_CONFIG_SRCDIR([src/distributions/DZIBG.cc])
AC_CONFIG_SRCDIR([src/distributions/DBEEG.cc])
AC_CONFIG_MACRO_DIR([m4])
dnl The configure process creates a header file called config.h
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

dnl libtool and ltdl configureation
LT_PREREQ(2.2.6)
LT_CONFIG_LTDL_DIR([libltdl])
LT_INIT([dlopen disable-static win32-dll])
LTDL_INIT([recursive])

dnl Indicate C++
AC_PROG_CXX

dnl Optionally, reference the Rmath library
dnl AC_DEFINE(MATHLIB_STANDALONE, 1, [Define if you have standalone R math library])

case "${host_os}" in 
    mingw*)
        win=true ;;
    *)
        win=false ;;
esac
AM_CONDITIONAL(WINDOWS, test x$win = xtrue)

jagsmoddir=/usr/lib/x86_64-linux-gnu/JAGS/modules-${JAGS_MAJOR}
AC_SUBST(jagsmoddir)

AC_CONFIG_FILES([
    Makefile
    libltdl/Makefile
    src/Makefile
    src/distributions/Makefile
])
AC_OUTPUT