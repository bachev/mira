dnl @synopsis AC_LATEX_INPUT(<package>,<class>,<variable>)
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if package in <class> exists and set <variable> to
dnl the right value (yes or no) Use \input instance of \usepackage
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_LATEX_PACKAGE_INPUT],[
if test "$[ac_cv_latex_class_]translit($2,[-],[_])" = "" ;
then
	AC_LATEX_CLASS($2,boretti_classesansparametre)
	export boretti_classesansparametre;
else
	boretti_classesansparametre=$[ac_cv_latex_class_]translit($2,[-],[_]) ;
	export boretti_classesansparemetre;
fi;
if test $boretti_classesansparametre = "no" ;
then
    AC_MSG_ERROR([Unable to find $1 class])
fi
AC_CACHE_CHECK([for $1 in class $2, using input insteance of usepackage],[ac_cv_latex_i_]translit($1,[-],[_])[_]translit($2,[-],[_]),[
AC_LATEX_TEST([
\documentclass{$2}
\input $1
\begin{document}
\end{document}
],[ac_cv_latex_i_]translit($1,[-],[_])[_]translit($2,[-],[_]))
]
$3=$[ac_cv_latex_i_]translit($1,[-],[_])[_]translit($2,[-],[_]); export $3;
AC_SUBST($3)
])
