# ===========================================================================
#   http://autoconf-archive.cryp.to/ac_cxx_full_specialization_syntax.html
# ===========================================================================
#
# SYNOPSIS
#
#   AC_CXX_FULL_SPECIALIZATION_SYNTAX
#
# DESCRIPTION
#
#   If the compiler recognizes the full specialization syntax, define
#   HAVE_FULL_SPECIALIZATION_SYNTAX.
#
# LICENSE
#
#   Copyright (c) 2008 Todd Veldhuizen
#   Copyright (c) 2008 Luc Maisonobe <luc@spaceroots.org>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AC_CXX_FULL_SPECIALIZATION_SYNTAX],
[AC_CACHE_CHECK(whether the compiler recognizes the full specialization syntax,
ac_cv_cxx_full_specialization_syntax,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
template<class T> class A        { public : int f () const { return 1; } };
template<>        class A<float> { public:  int f () const { return 0; } };],[
A<float> a; return a.f();],
 ac_cv_cxx_full_specialization_syntax=yes, ac_cv_cxx_full_specialization_syntax=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_full_specialization_syntax" = yes; then
  AC_DEFINE(HAVE_FULL_SPECIALIZATION_SYNTAX,,
            [define if the compiler recognizes the full specialization syntax])
fi
])
