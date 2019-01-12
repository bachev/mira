dnl @synopsis PETI_WITH_DMALLOC
dnl
dnl @obsoleted Renamed to AX_WITH_DMALLOC.
dnl
dnl Let the user enable/disable dmalloc library support. See
dnl <http://www.dmalloc.org/> for further details.
dnl
dnl Calling this macro defines a user command-line flag
dnl "--with-dmalloc". Furthermore, "-IPREFIX/include" will be added to
dnl "$CPPFLAGS", "-LPREFIX/lib" to "$LDFLAGS", and "-DDEBUG_DMALLOC"
dnl and "-DDMALLOC_FUNC_CHECK" to "$CPPFLAGS".
dnl
dnl To enable dmalloc support in your code, add this snippet to tho
dnl appropriate header files:
dnl
dnl   #ifdef DEBUG_DMALLOC
dnl   #  include <dmalloc.h>
dnl   #endif
dnl
dnl @category Obsolete
dnl @author Peter Simons <simons@cryp.to>
dnl @version 2006-06-04
dnl @license AllPermissive

AC_DEFUN([PETI_WITH_DMALLOC], [
AC_MSG_CHECKING(whether to use the dmalloc library)
AC_ARG_WITH(dmalloc,
[  --with-dmalloc[=PREFIX]  Compile with dmalloc library],
if test "$withval" = "" -o "$withval" = "yes"; then
    ac_cv_dmalloc="/usr/local"
else
    ac_cv_dmalloc="$withval"
fi
AC_MSG_RESULT(yes)
CPPFLAGS="$CPPFLAGS -DDEBUG_DMALLOC -DDMALLOC_FUNC_CHECK -I$ac_cv_dmalloc/include"
LDFLAGS="$LDFLAGS -L$ac_cv_dmalloc/lib"
LIBS="$LIBS -ldmalloc"
,AC_MSG_RESULT(no))
])dnl
