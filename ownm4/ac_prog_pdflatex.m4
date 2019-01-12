dnl @synopsis AC_PROG_PDFLATEX
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if pdflatex is installed. If pdflatex is installed,
dnl it set $pdflatex to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_PDFLATEX],[
AC_CHECK_PROGS(pdflatex,[pdflatex],no)
export pdflatex;
if test $pdflatex = "no" ;
then
	AC_MSG_ERROR([Unable to find a PDFLaTeX application]);
fi
AC_SUBST(pdflatex)
])
