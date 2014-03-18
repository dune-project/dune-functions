TESTS =  discretescalarglobalbasisfunctiontest

check_PROGRAMS = $(TESTS)

discretescalarglobalbasisfunctiontest_SOURCES  = discretescalarglobalbasisfunctiontest.cc
discretescalarglobalbasisfunctiontest_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS)
discretescalarglobalbasisfunctiontest_LDFLAGS  = $(AM_LDFLAGS) $(DUNE_LDFLAGS) $(DUNEMPILDFLAGS)
discretescalarglobalbasisfunctiontest_LDADD    = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(DUNEMPILIBS) $(LDADD)

EXTRA_DIST = CMakeLists.txt

CLEANFILES =

include $(top_srcdir)/am/global-rules
