# ===========================================================================
#        http://autoconf-archive.cryp.to/ac_cxx_cxxflags_std_lang.html
# ===========================================================================
#
# SYNOPSIS
#
#   AC_CXX_CXXFLAGS_STD_LANG(CXX-FLAGS)
#
# DESCRIPTION
#
#   Append to variable CXX-FLAGS the set of compile-time flags that should
#   be passed to the C++ compiler in order to enable use of C++ features as
#   defined in the ANSI C++ standard (eg. use of standard iostream classes
#   in the `std' namespace, etc.). Note that if you use GNU Libtool you may
#   need to prefix each of those flags with `-Xcompiler' so that Libtool
#   doesn't discard them (see Libtool's manual and `AC_LIBTOOLIZE_CFLAGS').
#
# LICENSE
#
#   Copyright (c) 2008 Ludovic Court�s <ludo@chbouib.org>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AC_CXX_CXXFLAGS_STD_LANG],
  [AC_REQUIRE([AC_CXX_COMPILER_VENDOR])
   case "$ac_cv_cxx_compiler_vendor" in
     sgi)    $1="$$1 -LANG:std -exceptions";;
     hp)     $1="$$1 -AA";;
   esac])
