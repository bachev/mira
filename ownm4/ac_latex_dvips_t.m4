dnl @synopsis AC_LATEX_DVIPS_T(<paper>,<var>) or AC_LATEX_DVIPS_T(<paper>,<var>,on|off)
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl This macro test if dvips -o ... -t <paper> works. When using the on
dnl option, test if dvips -o ... -t <paper> -t landscape works. if it
dnl works, set $var to yes, else $var="no"
dnl
dnl @category Obsolete
dnl @author Boretti Mathieu <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

AC_DEFUN([AC_LATEX_DVIPS_T],[
AC_REQUIRE([AC_LATEX_CLASS_BOOK])
if test "$3" = "on" ;
then
_ac_latex_dvips_local=" -t landscape" ; export _ac_latex_dvips_local ;
else
_ac_latex_dvips_local=" " ; export _ac_latex_dvips_local ;
fi
AC_CACHE_CHECK([for option -t $1 $_ac_latex_dvips_local with dvips],[ac_cv_dvips_t_]translit($1,[-],[_])[_]translit($3,[-],[_]),[
rm -rf .dvips
mkdir .dvips
cd .dvips
cat > test.tex << EOF
\documentclass{book}
\begin{document}
Test
\end{document}
EOF
$latex test.tex 1>/dev/null 2>&1
[ac_cv_dvips_t_]translit($1,[-],[_])[_]translit($3,[-],[_])="yes"; export [ac_cv_dvips_t_]translit($1,[-],[_])[_]translit($3,[-],[_]);
$dvips -o test.ps test.dvi -t $1 $_ac_latex_dvips_local 2>&1 1>/dev/null | (grep "dvips: no match for papersize" 1>/dev/null 2>&1 && [ac_cv_dvips_t_]translit($1,[-],[_])[_]translit($3,[-],[_])="no"; export [ac_cv_dvips_t_]translit($1,[-],[_])[_]translit($3,[-],[_]))
cd ..
])
$2=$[ac_cv_dvips_t_]translit($1,[-],[_])[_]translit($3,[-],[_]); export $2;
AC_SUBST($2)
])
