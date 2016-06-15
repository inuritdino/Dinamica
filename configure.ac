#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Dinamica], [1.0], [elias.potapov@gmail.com],[dinamica],[dl.sv.nongnu.org/releases/din/])
AC_COPYRIGHT([Copyright 2013, 2014 Elias Potapov])
AC_REVISION([$Revision: 1.6 $])
AC_CONFIG_SRCDIR([src/din_main.c])
AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

#We build a library so we need the following
AC_PROG_RANLIB
# Checks for programs.
AC_PROG_CC
#AC_PROG_INSTALL
AC_PROG_LN_S
AM_PROG_AR
dnl We have to include the Gnuplot support.
dnl First, check if the --enable-gnuplot argument was supplied
AC_ARG_ENABLE(gnuplot,AC_HELP_STRING([--enable-gnuplot=<val>],
			[enables the Gnuplot support, yes is default]),[],[])
dnl If the argument was there, start actual check for the program
if test "$enable_gnuplot" == "" ; then
   AC_CHECK_PROG(GNUPLOT, gnuplot, yes)
   if test x"$GNUPLOT" != x"yes"; then
      AC_MSG_ERROR([[Gnuplot is NOT found. Gnuplot is included into the default
      installation. If you do not need it with Dinamica put --enable-gnuplot=no,
      otherwise consider changing the PATH variable to include a path to the
      gnuplot executable.]])
      dnl AC_DEFINE(HAVE_GNUPLOT,0,[Define to 1 if we have Gnuplot. 0 otherwise.])
   else
      AC_DEFINE(HAVE_GNUPLOT,1,[Define to 1 if we have Gnuplot. 0 otherwise.])
   fi
else
   if test x"$enable_gnuplot" != x"no"; then
      AC_CHECK_PROG(GNUPLOT, gnuplot, yes)
      if test x"$GNUPLOT" != x"yes"; then
      	 AC_MSG_ERROR([[Gnuplot is NOT found, although you required it.
      	 If there is no Gnuplot in the system put --enable-gnuplot=no.]])
      else
      	 AC_DEFINE(HAVE_GNUPLOT,1,[Define to 1 if we have Gnuplot. 0 otherwise.])
      fi
   else
         AC_DEFINE(HAVE_GNUPLOT,0,[Define to 1 if we have Gnuplot. 0 otherwise.])
   fi
fi

dnl Set the paths that Dinamica will use for the compilation
AC_DEFINE_UNQUOTED(DIN_CPPFLAGS,"$CPPFLAGS",[Path to the headers Dinamica will use])
AC_DEFINE_UNQUOTED(DIN_LDFLAGS,"$LDFLAGS",[Path to the libraries Dinamica will use])
dnl LIBDIN_PATH is appended to the -L flag of the dinamical compiler. Whether the
dnl user specified the ${prefix} is checked in the din_main.c.
AC_DEFINE_UNQUOTED(LIBDIN_PATH,"${prefix}/lib",[Path to the libdin.a])

# Checks for libraries.
AC_CHECK_LIB([m], [pow],[],[AC_MSG_ERROR([Error: libm not found])])
AC_CHECK_LIB([gslcblas], [cblas_sdsdot],[],[AC_MSG_ERROR([Error: gslcblas not found])],[-lgsl])
AC_CHECK_LIB([gsl], [gsl_log1p],[],[AC_MSG_ERROR([Error: gsl not found])],[-lgslcblas])
AC_CHECK_LIB([gsl], [gsl_odeiv2_step_alloc],[],[AC_MSG_ERROR([Error: perhaps wrong version of GSL])],[-lgslcblas])


# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h time.h gsl/gsl_multiroots.h],[],
			   [AC_MSG_ERROR([Error: headers cannot be found])])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([bzero memset pow putenv sqrt strcasecmp strchr strtol clock])

AC_CONFIG_FILES([
	Makefile
	src/Makefile
	doc/Makefile
])
AC_OUTPUT
