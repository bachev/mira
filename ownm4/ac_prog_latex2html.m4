dnl @synopsis AC_PROG_LATEX2HTML
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if latex2html is installed. If latex2html is
dnl installed, it set $latex2html to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_LATEX2HTML],[
AC_CHECK_PROGS(latex2html,[latex2html],no)
export latex2html;
if test $latex2html = "no" ;
then
	AC_MSG_ERROR([Unable to find a LaTeX2html application]);
fi
AC_SUBST(latex2html)
])
