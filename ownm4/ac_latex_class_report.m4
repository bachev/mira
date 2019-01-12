dnl @synopsis AC_LATEX_CLASS_REPORT
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl same as AC_LATEX_CLASS(report,report)
dnl
dnl @category Obsolete
dnl @author Boretti Mathieu <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_LATEX_CLASS_REPORT],[
AC_LATEX_CLASS(report,report)
if test $report = "no";
then
    AC_MSG_ERROR([Unable to find the report class])
fi
])
