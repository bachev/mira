dnl @synopsis AC_PROG_GNUPLOT
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if gnuplot is installed. If gnuplot is installed,
dnl it set $gnuplot to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_GNUPLOT],[
AC_CHECK_PROGS(gnuplot,[gnuplot],no)
export gnuplot;
if test $gnuplot = "no" ;
then
	AC_MSG_ERROR([Unable to find a gnuplot application]);
fi
AC_SUBST(gnuplot)
])
