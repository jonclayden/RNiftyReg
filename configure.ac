#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([RNiftyReg], [2.6.1], [code@clayden.org])

AC_LANG([C++])

# Get compiler flags from R
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

# Checks for programs.
AC_PROG_CXX
# AC_PROG_CC

# Checks for OpenMP.
m4_include([m4/ax_openmp.m4])
AX_OPENMP

AC_SUBST(OPENMP_CXXFLAGS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
