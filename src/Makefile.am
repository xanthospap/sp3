##
##  The library will be called: (lib)"cppdoris".
## ----------------------------------------------
##
lib_LTLIBRARIES = libcppsp3.la

##
##  Only use the following flags for debuging purposes.
## ----------------------------------------------------
##
libcppsp3_la_CXXFLAGS = \
-std=c++17 \
-g \
-pg \
-Wall \
-Wextra \
-Werror \
-pedantic \
-W \
-Wshadow \
-Winline \
-Wdisabled-optimization \
-march=native \
-DDEBUG

##
##  Header files (distributed) installed in /$(includedir)/$(package name).
## ------------------------------------------------------------------------
##
pkginclude_HEADERS = \
  sp3_flag_details.hpp \
  sp3flag.hpp \
  satellite.hpp \
  sp3.hpp

##
##  Source files (distributed).
## ----------------------------
##
dist_libcppsp3_la_SOURCES = \
  sp3flag.cpp \
  sp3_read_header.cpp \
  sp3.cpp
