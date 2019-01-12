dnl @synopsis AC_PROG_FIG2DEV
dnl
dnl @obsoleted Replaced by VL_PROG_FIG2DEV.
dnl
dnl This macro test if fig2dev is installed. If fig2dev is installed,
dnl it set $fig2dev to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2005-01-25
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_FIG2DEV],[
AC_CHECK_PROGS(fig2dev,[fig2dev],no)
export fig2dev;
if test $fig2dev = "no" ;
then
	AC_MSG_ERROR([Unable to find a fig2dev application]);
fi
AC_SUBST(fig2dev)
])
