dnl @synopsis AC_ECHO_N(PATH)
dnl
dnl this is somewhat like the macro _AC_ECHO_N from autoconf 2.4x
dnl defined here for use in autoconf 2.1x, add the _ when you use 2.4x
dnl
dnl @category Misc
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

# _AC_ECHO_N(STRING, [FD = AS_MESSAGE_FD])
# ------------------------------------
# Same as _AS_ECHO, but echo doesn't return to a new line.
AC_DEFUN([AC_ECHO_N],[dnl
[echo $ac_n "$1""$ac_c" >&AC_FD_MSG])
