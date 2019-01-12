dnl @synopsis AC_LATEX_PACKAGE_FONTENC
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if \usepackage[T1]{fontenc} works. If yes it set
dnl $fontenc="T1" else if \usepackage[OT1]{fontenc} works, set
dnl $fontenc="OT1" else ERROR
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

define(_AC_LATEX_PACKAGE_FONTENC_INTERNE,[
changequote(*, !)dnl
\documentclass{book}
\usepackage[$1]{fontenc}
\begin{document}
\end{document}
changequote([, ])dnl

])

AC_DEFUN([AC_LATEX_PACKAGE_FONTENC],[
    AC_LATEX_CLASS_BOOK
    AC_CACHE_CHECK([for fontenc],[ac_cv_latex_package_fontenc_opt],[
        _AC_LATEX_TEST([_AC_LATEX_PACKAGE_FONTENC_INTERNE(T1)],[ac_cv_latex_package_fontenc_opt])
        if test $ac_cv_latex_package_fontenc_opt = "yes" ;
        then
            ac_cv_latex_package_fontenc_opt="T1"; export ac_cv_latex_package_fontenc_opt;
        else
            _AC_LATEX_TEST([_AC_LATEX_PACKAGE_FONTENC_INTERNE(OT1)],[ac_cv_latex_package_fontenc_opt])
            if test $ac_cv_latex_package_fontenc_opt = "yes" ;
            then
                ac_cv_latex_package_fontenc_opt="OT1"; export ac_cv_latex_package_fontenc_opt;
            fi
        fi

    ])
    if test $ac_cv_latex_package_fontenc_opt = "no" ;
    then
        AC_MSG_ERROR([Unable to use fontenc with T1 nor OT1])
    fi
    fontenc=$ac_cv_latex_package_fontenc_opt ; export fontenc ;
    AC_SUBST(fontenc)
])
