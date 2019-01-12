dnl @synopsis AC_PROG_BIBTEX
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if bibtex is installed. If bibtex is installed, it
dnl set $bibtex to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_BIBTEX],[
AC_CHECK_PROGS(bibtex,bibtex,no)
export bibtex;
if test $bibtex = "no" ;
then
	AC_MSG_ERROR([Unable to find a BiTex application]);
fi
AC_SUBST(bibtex)
])
