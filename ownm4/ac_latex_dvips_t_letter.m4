dnl @synopsis AC_LATEX_DVIPS_T_LETTER
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl same as AC_LATEX_DVIPS_T(letter,dvips_t_letter)
dnl
dnl @category Obsolete
dnl @author Boretti Mathieu <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_LATEX_DVIPS_T_LETTER],[
AC_LATEX_DVIPS_T(letter,dvips_t_letter)
if test $dvips_t_letter = "no";
then
    AC_MSG_ERROR([Unable to find the -t letter option in dvips])
fi
])
