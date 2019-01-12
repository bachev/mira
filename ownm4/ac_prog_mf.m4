dnl @synopsis AC_PROG_MF
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if mf is installed. If mf is installed, it set $mf
dnl to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_MF],[
AC_CHECK_PROGS(mf,[mf mfw mf-nowin],no)
export mf;
if test $mf = "no" ;
then
	AC_MSG_ERROR([Unable to find a mf application]);
fi
AC_SUBST(mf)
])
