dnl @synopsis AC_SUBST_DIR(VARNAME, [DIR])
dnl
dnl @obsoleted The macro AC_DEFINE_DIR provides this functionality directly now.
dnl
dnl This macro substitutes (with AC_SUBST) VARNAME with the expansion
dnl of itself or the DIR variable if specified, taking care of fixing
dnl up ${prefix} and such.
dnl
dnl Side effect: VARNAME is replaced with the expansion.
dnl
dnl AC_SUBST_DIR bases on Alexandre Oliva's AC_DEFINE_DIR macro.
dnl
dnl Examples:
dnl
dnl    AC_SUBST_DIR(DATADIR)
dnl
dnl @category Obsolete
dnl @author Stepan Kasal <kasal@ucw.cz>
dnl @author Mathias Hasselmann <mathias.hasselmann@gmx.de>
dnl @version 2005-07-29
dnl @license AllPermissive

AC_DEFUN([AC_SUBST_DIR], [
        ifelse($2,,,$1="[$]$2")
        $1=`(
            test "x$prefix" = xNONE && prefix="$ac_default_prefix"
            test "x$exec_prefix" = xNONE && exec_prefix="${prefix}"
            eval $1=\""[$]$1"\"
            eval echo \""[$]$1"\"
        )`
        AC_SUBST($1)
])
