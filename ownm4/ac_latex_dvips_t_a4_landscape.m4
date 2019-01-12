dnl @synopsis AC_LATEX_DVIPS_T_A4_LANDSCAPE
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl same as AC_LATEX_DVIPS_T(a4,dvips_t_a4_landscape,on)
dnl
dnl @category Obsolete
dnl @author Boretti Mathieu <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_LATEX_DVIPS_T_A4_LANDSCAPE],[
AC_REQUIRE([AC_LATEX_DVIPS_T_A4])
AC_LATEX_DVIPS_T(a4,dvips_t_a4_landscape,on)
if test $dvips_t_a4_landscape = "no";
then
    AC_MSG_ERROR([Unable to find the -t a4 -t landscape option in dvips])
fi
])
