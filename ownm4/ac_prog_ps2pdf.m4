dnl @synopsis AC_PROG_PS2PDF
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if ps2pdf is installed. If ps2pdf is installed, it
dnl set $ps2pdf to the right value
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_PS2PDF],[
AC_CHECK_PROGS(ps2pdf,[ps2pdf14 ps2pdf13 ps2pdf12 ps2pdf],no)
export ps2pdf;
if test $ps2pdf = "no" ;
then
	AC_MSG_ERROR([Unable to find a ps2pdf application]);
fi
AC_SUBST(ps2pdf)
])
