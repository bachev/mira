dnl @synopsis AC_PROG_JPEG2PS
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if jpeg2ps is installed. If jpeg2ps is installed,
dnl it set $jpeg2ps to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_JPEG2PS],[
AC_CHECK_PROGS(jpeg2ps,[jpeg2ps],no)
export jpeg2ps;
if test $jpeg2ps = "no" ;
then
	AC_MSG_ERROR([Unable to find a jpeg2ps application]);
fi
AC_SUBST(jpeg2ps)
])
