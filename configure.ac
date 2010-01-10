# Configure.ac to produce configure file by autoconf
#
# Copyright (C) 2010 Chunhua Wu
# Lincese: GPL(>=2.0)
#

AC_PREREQ([2.64])

AC_INIT([DESCRIPTION])

# Checks for GCC 
AC_PROG_CC

AC_CONFIG_SRCDIR([src])

AC_CHECK_HEADER([gsl/gsl_version.h], ,
             [AC_MSG_ERROR([GNU Scientific Library Not Found.])])

# Checks for libraries.
AC_CHECK_LIB(gsl, gsl_linalg_cholesky_invert,
             [LIBS="${LIBS} -lgsl -lgslcblas"], 
             [AC_MSG_ERROR([GNU Scientific Library Version Not New enough, Suggest V >= 1.12.])],
             [-lgslcblas])

# Checks for header files.
 AC_CHECK_HEADERS([math.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.


AC_OUTPUT()
