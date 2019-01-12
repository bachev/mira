dnl @synopsis AC_DEFINE_PATH_STYLE ([defvar-name])
dnl
dnl _AC_DEFINE(PATH_STYLE) describing the filesys interface. The value
dnl is numeric, where the basetype is encoded as 16 = dos/win, 32 =
dnl unix, 64 = url/www, 0 = other
dnl
dnl note that there could be a combination of the values that should
dnl lead you to accept multiple forms of PATH_SEP and DIR_SEP
dnl
dnl @category Misc
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([AC_DEFINE_PATH_STYLE], [
AC_CACHE_CHECK([for path style], ac_cv_path_style,
[
  if test -z "$ac_cv_path_style"; then
    case "$target_os" in
      *djgpp | *mingw32* | *emx*) ac_cv_path_style="dos" ;;
      *) ac_cv_path_style="unix" ;;    # it is just the default ;-)
    esac
    if test "$ac_cv_path_style" = "unix" ; then
      _exec_prefix=`eval "echo $exec_prefix"`
      _exec_prefix=`eval "echo $_exec_prefix"`
      case "$_exec_prefix" in
        *:*) ac_cv_path_style="url" ;;
        *\\) ac_cv_path_style="dos" ;;
      esac
    fi
  fi
])
case "$ac_cv_path_style" in
  *dos*)  ac_define_path_style=16 ;;
  *unix*) ac_define_path_style=32 ;;
  *url*)  ac_define_path_style=64 ;;
  *mac*)  ac_define_path_style=128 ;;
  *) ac_define_path_style=1
esac
AC_DEFINE_UNQUOTED(PATH_STYLE,$ac_define_path_style,dnl
   [path style 16=dos 32=unix 64=url 128=mac])dnl
])
