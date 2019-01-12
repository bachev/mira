dnl @synopsis AC_LATEX_DVIPS_T_LETTER_LANDSCAPE
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl same as AC_LATEX_DVIPS_T(letter,dvips_t_letter_landscape,on)
dnl
dnl @category Obsolete
dnl @author Boretti Mathieu <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_LATEX_DVIPS_T_LETTER_LANDSCAPE],[
AC_REQUIRE([AC_LATEX_DVIPS_T_LETTER])
AC_LATEX_DVIPS_T(letter,dvips_t_letter_landscape,on)
if test $dvips_t_letter_landscape = "no";
then
    AC_MSG_ERROR([Unable to find the -t letter -t landscape option in dvips])
fi
])
