dnl @synopsis PETI_ENABLE_DYNAMIC_LINK
dnl
dnl @obsoleted Deprecated with the advent of GNU Libtool.
dnl
dnl Add a command-line flag to enable/disable dynamic linking.
dnl
dnl Calling this macro adds the flag "--enable-dynamic-link" to
dnl command-line. When disabled, the compiler/linker flag "-static" is
dnl added to "$LDFLAGS". The default is dynamic linkage.
dnl
dnl @category Obsolete
dnl @author Peter Simons <simons@cryp.to>
dnl @version 2006-06-04
dnl @license AllPermissive

AC_DEFUN([PETI_ENABLE_DYNAMIC_LINK], [
AC_MSG_CHECKING(what kind of binaries we shall create)
AC_ARG_ENABLE(dynamic-link,
[  --enable-dynamic-link   Create dynamically-linked binaries (default)],
if test "$enableval" = "yes"; then
    AC_MSG_RESULT(dynamically linked)
else
    LDFLAGS="$LDFLAGS -static"
    AC_MSG_RESULT(statically linked)
fi,
AC_MSG_RESULT(dynamically linked))
])dnl
