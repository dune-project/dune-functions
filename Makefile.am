EXTRA_DIST = dune.module
EXTRA_DIST += CMakeLists.txt dune-functions-config.cmake.in dune-functions-version.cmake.in

SUBDIRS = cmake doc dune m4

DISTCHECK_CONFIGURE_FLAGS = --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT)  CXX="$(CXX)" CC="$(CC)"

# don't follow the full GNU-standard
# we need automake 1.9 or newer
AUTOMAKE_OPTIONS = foreign 1.9

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules
