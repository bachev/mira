dnl @synopsis AC_CHECK_CC_OPT(flag, cachevarname)
dnl
dnl @obsoleted Use CFLAGS/CXXFLAGS related macros as soon as possible.
dnl
dnl AC_CHECK_CC_OPT(-fvomit-frame,vomitframe) would show a message as
dnl like "checking wether gcc accepts -fvomit-frame ... no" and sets
dnl the shell-variable $vomitframe to either "-fvomit-frame" or (in
dnl this case) just a simple "". In many cases you would then call
dnl AC_SUBST(_fvomit_frame_,$vomitframe) to create a substitution that
dnl could be fed as "CFLAGS = @_funsigned_char_@ @_fvomit_frame_@.
dnl
dnl In consequence this function is much more general than their
dnl specific counterparts like ac_cxx_rtti.m4 that will test for
dnl -fno-rtti -fno-exceptions.
dnl
dnl This macro will be obsolete in the very near future - there should
dnl be two macros AX_CFLAGS_OPTION and AX_CXXFLAGS_OPTION that will add
dnl directly to the CFLAGS/CXXFLAGS unless a ac_subst-variable was
dnl given. Remind me of doing so if I forget about writing it all up -
dnl to use $CC-cc directly in here instead of ac_compile macro is bad
dnl anyway.
dnl
dnl @category Obsolete
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([AC_CHECK_CC_OPT],
[AC_CACHE_CHECK(whether ${CC-cc} accepts [$1], [$2],
[AC_SUBST($2)
echo 'void f(){}' > conftest.c
if test -z "`${CC-cc} -c $1 conftest.c 2>&1`"; then
  $2="$1"
else
  $2=""
fi
rm -f conftest*
])])
