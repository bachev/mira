dnl @synopsis MERK_PROG_TCL([min-version])
dnl
dnl Searches for tcl (tclsh and wish) in PATH and checks which version
dnl is installed. The macro bails out if either tcl is not found or the
dnl minimum version is not satisfied, unless minimum version is "0".
dnl
dnl Example:
dnl
dnl  MERK_PROG_TCL
dnl
dnl This checks for tcl and if not found, exits with an error. If
dnl found, it prints tcl path and version number.
dnl
dnl  MERK_PROG_TCL([8.0])
dnl
dnl Checks for tcl and exits with an error if its not found or the
dnl version is below 8.0.
dnl
dnl @category InstalledPackages
dnl @author Uwe Mayer <merkosh@hadiko.de>
dnl @version 2005-01-21
dnl @license GPLWithACException

AC_DEFUN([MERK_PROG_TCL], [
#-- check for tclsh in PATH
AC_PATH_PROG([TCLSH], [tclsh], [no])
if [(test x"$TCLSH" == x"no") && (test x"$1" != x"0")]; then
        AC_MSG_ERROR([tclsh not found])
fi

#-- check for wish in PATH
AC_PATH_PROG([WISH], [wish], [no])

#-- check vor tcl version
AC_MSG_CHECKING([tcl version])
tmp=$(tempfile)
echo ["puts [set tcl_version]"] >"$tmp"
version=$(tclsh "$tmp")
rm "$tmp"
AC_MSG_RESULT([$version])

#-- compare tcl version with min-version
required=$1
if [(test x"$1" != x"") && (test "${required/./}" -gt "${version/./}")]; then
        AC_MSG_ERROR([tcl version $1 required])
fi
])dnl
