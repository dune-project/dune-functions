AC_PREREQ(2.50)
DUNE_AC_INIT
AC_CONFIG_SRCDIR([dune/functions/common/differentiablefunction.hh])
AC_CONFIG_HEADERS([config.h])

DUNE_CHECK_ALL

AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  cmake/Makefile
  cmake/pkg/Makefile
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

DUNE_SUMMARY_ALL
