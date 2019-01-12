dnl @synopsis AC_CXX_NCEG_RESTRICT
dnl
dnl If the compiler supports the Numerical C Extensions Group restrict
dnl keyword, define HAVE_NCEG_RESTRICT.
dnl
dnl @category Cxx
dnl @author Todd Veldhuizen
dnl @author Luc Maisonobe <luc@spaceroots.org>
dnl @version 2004-02-04
dnl @license AllPermissive

AC_DEFUN([AC_CXX_NCEG_RESTRICT],
[AC_CACHE_CHECK(whether the compiler supports the Numerical C Extensions Group restrict keyword,
ac_cv_cxx_nceg_restrict,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
void add(int length, double * restrict a,
         const double * restrict b, const double * restrict c)
{ for (int i=0; i < length; ++i) a[i] = b[i] + c[i];}
],[double a[10], b[10], c[10];
for (int i=0; i < 10; ++i) { a[i] = 0.0;  b[i] = 0.0; c[i] = 0.0;}
add(10,a,b,c);
return 0;],
 ac_cv_cxx_nceg_restrict=yes, ac_cv_cxx_nceg_restrict=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_nceg_restrict" = yes; then
  AC_DEFINE(HAVE_NCEG_RESTRICT,,
            [define if  the compiler supports the Numerical C Extensions Group restrict keyword])
fi
])
