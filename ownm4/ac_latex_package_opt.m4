dnl @synopsis AC_LATEX_PACKAGE_OPT(<package>,<class>,<variable>,<option>)
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if package in <class> with option <option> exists
dnl and set <variable> to the right value (yes or no)
dnl
dnl @category Obsolete
dnl @author Mathieu Boretti <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

define(_AC_LATEX_PACKAGE_OPT_INTERNE,[
changequote(*, !)dnl
\documentclass{$2}
\usepackage[$3]{$1}
\begin{document}
\end{document}
changequote([, ])dnl

])

AC_DEFUN([AC_LATEX_PACKAGE_OPT],[
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
AC_CACHE_CHECK([for $1 in class $2 with $4 as option],[ac_cv_latex_]translit($1,[-],[_])[_]translit($2,[-],[_])[_]translit($4,[-],[_]),[
_AC_LATEX_TEST([
_AC_LATEX_PACKAGE_OPT_INTERNE($1,$2,$4)
],[ac_cv_latex_]translit($1,[-],[_])[_]translit($2,[-],[_])[_]translit($4,[-],[_]))
])
$3=$[ac_cv_latex_]translit($1,[-],[_])[_]translit($2,[-],[_])[_]translit($4,[-],[_]); export $3;
AC_SUBST($3)
])
