TESTS = hello

check_PROGRAMS = $(TESTS)

hello_SOURCES = hello.cc
hello_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS)
hello_LDFLAGS = $(AM_LDFLAGS) $(DUNE_LDFLAGS) $(DUNEMPILDFLAGS)
hello_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(DUNEMPILIBS) $(LDADD)

EXTRA_DIST = CMakeLists.txt

CLEANFILES =

include $(top_srcdir)/am/global-rules
