dnl @synopsis AC_LATEX_CLASSES([<class1>,<class2>,...],<var>)
dnl
dnl @obsoleted Replaced by the newer ACLTX_XXX set of macros.
dnl
dnl Test if class1 exists and if not class2 and so and set $var to the
dnl right value
dnl
dnl  AC_LATEX_CLASSES([allo,book,bnjour],book)
dnl  should set $book="book"
dnl
dnl @category Obsolete
dnl @author Boretti Mathieu <boretti@eig.unige.ch>
dnl @version 2006-07-16
dnl @license GPLWithACException

define(_AC_LATEX_CLASSES_INTERNE,[
	ifelse($#,1,[],$#,2,[
		AC_LATEX_CLASS($2,$1)
	],[
		AC_LATEX_CLASS($2,$1)
		if test "$$1" = "yes";
		then
			$1=$2 ; export $1 ;
		else
			_AC_LATEX_CLASSES_INTERNE($1,m4_shift(m4_shift($@)))
		fi;
	])
])

AC_DEFUN([AC_LATEX_CLASSES],[
	_AC_LATEX_CLASSES_INTERNE($2,$1)
	AC_SUBST($2)
])
