noinst_PROGRAMS = testSp3Flags testSp3Read testInterpolate

MCXXFLAGS = \
  -std=c++17 \
  -g \
  -pg \
  -Wno-unused-variable \
  -Wall \
  -Wextra \
  -Werror \
  -pedantic \
  -W \
  -Wshadow \
  -O2 \
  -march=native \
  -DDEBUG

AM_LDFLAGS = -lggeodesy -lggdatetime

testSp3Flags_SOURCES   = test_sp3_flags.cc
testSp3Flags_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testSp3Flags_LDADD     = $(top_srcdir)/src/libcppsp3.la

testSp3Read_SOURCES   = test_sp3_read.cc
testSp3Read_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testSp3Read_LDADD     = $(top_srcdir)/src/libcppsp3.la

testInterpolate_SOURCES   = test_sv_interpolation.cc
testInterpolate_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testInterpolate_LDADD     = $(top_srcdir)/src/libcppsp3.la
