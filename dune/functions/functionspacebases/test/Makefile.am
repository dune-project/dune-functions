TESTS =  gridviewfunctionspacebasistest

check_PROGRAMS = $(TESTS)

gridviewfunctionspacebasistest_SOURCES  = gridviewfunctionspacebasistest.cc
gridviewfunctionspacebasistest_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS)
gridviewfunctionspacebasistest_LDFLAGS  = $(AM_LDFLAGS) $(DUNE_LDFLAGS) $(DUNEMPILDFLAGS)
gridviewfunctionspacebasistest_LDADD    = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(DUNEMPILIBS) $(LDADD)

EXTRA_DIST = CMakeLists.txt

CLEANFILES =

include $(top_srcdir)/am/global-rules
