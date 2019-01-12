dnl @synopsis AC_LATEX_PACKAGE_AMSMATH
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if \usepackage{amsmath,amsfonts} works. If yes, it
dnl set $amsmath="\usepackage{amsmath,amsfonts}" Else if
dnl \usepackage{amstex} works, set $amsmath="\usepackage{amstex}" else
dnl ERROR
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_LATEX_PACKAGE_AMSMATH],[
AC_LATEX_CLASS_BOOK
AC_CACHE_CHECK([for amsmath],[ac_cv_latex_package_f_amsmath],[
_AC_LATEX_TEST([
\documentclass{book}
\usepackage{amsmath,amsfonts}
\begin{document}
\end{document}
],[ac_cv_latex_package_f_amsmath])
if test $ac_cv_latex_package_f_amsmath = "yes" ;
then
    [ac_cv_latex_package_f_amsmath]="\\usepackage{amsmath,amsfonts}" ; export [ac_cv_latex_package_f_amsmath] ;
else
    _AC_LATEX_TEST([
    \documentclass{book}
    \usepackage{amstex}
    \begin{document}
    \end{document}
    ],[ac_cv_latex_package_f_amsmath])
    if test $ac_cv_latex_package_f_amsmath = "yes" ;
    then
        [ac_cv_latex_package_f_amsmath]="\\usepackage{amstex}" ; export [ac_cv_latex_package_f_amsmath] ;
    else
        AC_MSG_ERROR([Unable to find amsmath])
    fi
fi
])
amsmath=$[ac_cv_latex_package_f_amsmath]; export amsmath;
AC_SUBST(amsmath)
])
