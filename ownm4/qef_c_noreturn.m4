dnl @synopsis QEF_C_NORETURN
dnl
dnl @obsoleted Test "__GNUC__+0 >= 2 && __GNUC_MINOR__+0 >= 5" instead.
dnl
dnl Check if we can use GCC's __noreturn__ attribute in prototypes to
dnl indicate that functions never return. This can be used by the
dnl compiler to do some extra optimizations.
dnl
dnl FUNCATTR_NORETURN is defined as what we should put at the end of
dnl function prototypes to achieve this. If the compiler doesn't
dnl support it then it is defined as empty.
dnl
dnl An example of a a function's prototype and implementation using
dnl this macro:
dnl
dnl   void this_function_never_returns (void) FUNCATTR_NORETURN;
dnl
dnl   void this_function_never_returns (void) {
dnl      exit (0);
dnl   }
dnl
dnl @category Obsolete
dnl @author Geoff Richards <ctzgpr@scs.leeds.ac.uk>
dnl @version 2005-01-23
dnl @license AllPermissive

AC_DEFUN([QEF_C_NORETURN],
[AC_REQUIRE([AC_PROG_CC])
AC_MSG_CHECKING(whether the C compiler (${CC-cc}) accepts noreturn attribute)
AC_CACHE_VAL(qef_cv_c_noreturn,
[qef_cv_c_noreturn=no
AC_TRY_COMPILE(
[#include <stdio.h>
void f (void) __attribute__ ((noreturn));
void f (void)
{
   exit (1);
}
], [
   f ();
],
[qef_cv_c_noreturn="yes";  FUNCATTR_NORETURN_VAL="__attribute__ ((noreturn))"],
[qef_cv_c_noreturn="no";   FUNCATTR_NORETURN_VAL="/* will not return */"])
])

AC_MSG_RESULT($qef_cv_c_noreturn)
AC_DEFINE_UNQUOTED(FUNCATTR_NORETURN, $FUNCATTR_NORETURN_VAL)
])dnl
