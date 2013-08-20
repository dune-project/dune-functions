TESTS = differentiablefunctiontest

check_PROGRAMS = $(TESTS)

differentiablefunctiontest_SOURCES = differentiablefunctiontest.cc
differentiablefunctiontest_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS)
differentiablefunctiontest_LDFLAGS = $(AM_LDFLAGS) $(DUNE_LDFLAGS) $(DUNEMPILDFLAGS)
differentiablefunctiontest_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(DUNEMPILIBS) $(LDADD)

EXTRA_DIST = CMakeLists.txt

CLEANFILES =

include $(top_srcdir)/am/global-rules
