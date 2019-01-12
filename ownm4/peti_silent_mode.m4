dnl @synopsis PETI_SILENT_MODE(on|off)
dnl
dnl @obsoleted Renamed to AX_SILENT_MODE.
dnl
dnl Temporarily disable console output. For example:
dnl
dnl   PETI_SILENT_MODE(on)    dnl be silent
dnl   AC_PROG_CXX
dnl   PETI_SILENT_MODE(off)   dnl talk to me again
dnl   AC_PROG_RANLIB
dnl
dnl Many thanks to Paolo Bonzini for proposing this macro.
dnl
dnl @category Obsolete
dnl @author Peter Simons <simons@cryp.to>
dnl @version 2006-06-04
dnl @license AllPermissive

AC_DEFUN([PETI_SILENT_MODE],
  [
  case "$1" in
    on)
      exec 6>/dev/null
      ;;
    off)
      exec 6>&1
      ;;
    *)
      AC_MSG_ERROR(Silent mode can only be switched "on" or "off".)
      ;;
  esac
  ])dnl
