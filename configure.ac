# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.62])
DUNE_AC_INIT # gets module version from dune.module file
AC_CONFIG_SRCDIR([dune/functions/common/differentiablefunction.hh])
AC_CONFIG_HEADERS([config.h])

# we need no more than the standard DE-stuff
# this module depends on dune-localfunctions dune-grid dune-typetree
# this implies checking for [dune-common], [dune-geometry], [dune-localfunctions], [dune-grid], [dune-typetree]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  cmake/Makefile
  cmake/modules/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/functions/Makefile
  dune/functions/common/Makefile
  dune/functions/common/test/Makefile
  dune/functions/functionspacebases/Makefile
  dune/functions/functionspacebases/test/Makefile
  dune/functions/gridfunctions/Makefile
  dune/functions/gridfunctions/test/Makefile
  dune/functions/test/Makefile
  m4/Makefile
  dune-functions.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
