# ===========================================================================
#           http://autoconf-archive.cryp.to/acx_blas_f77_func.html
# ===========================================================================
#
# SYNOPSIS
#
#   ACX_BLAS_F77_FUNC([ACTION-IF-PASS[, ACTION-IF-FAIL[, ACTION-IF-CROSS-COMPILING]])
#   ACX_BLAS_WITH_F77_FUNC([ACTION-IF-FOUND-AND-PASS[, ACTION-IF-NOT-FOUND-OR-FAIL]])
#
# DESCRIPTION
#
#   These macros are intended as a supplement to the ACX_BLAS macro, to
#   verify that BLAS functions are properly callable from Fortran. This is
#   necessary, for example, if you want to build the LAPACK library on top
#   of the BLAS.
#
#   ACX_BLAS_F77_FUNC uses the defined BLAS_LIBS and Fortran environment to
#   check for compatibility, and takes a specific action in case of success,
#   resp. failure, resp. cross-compilation.
#
#   ACX_BLAS_WITH_F77_FUNC is a drop-in replacement wrapper for ACX_BLAS
#   that calls ACX_BLAS_F77_FUNC after detecting a BLAS library and rejects
#   it on failure (i.e. pretends that no library was found).
#
# LICENSE
#
#   Copyright (c) 2008 Jaroslav Hajek <highegg@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([ACX_BLAS_F77_FUNC], [
AC_PREREQ(2.50)
AC_REQUIRE([ACX_BLAS])

# F77 call-compatibility checks
if test "$cross_compiling" = yes ; then
	ifelse($3, ,$1,$3)
elif test x"$acx_blas_ok" = xyes; then
	LIBS="$BLAS_LIBS $LIBS"
	AC_LANG_PUSH(Fortran 77)
# LSAME check (LOGICAL return values)
	AC_MSG_CHECKING([whether LSAME is called correctly from Fortran])
	AC_RUN_IFELSE(AC_LANG_PROGRAM(,[[
      logical lsame,w
      external lsame
      character c1,c2
      c1 = 'A'
      c2 = 'B'
      w = lsame(c1,c2)
      if (w) stop 1
      w = lsame(c1,c1)
      if (.not. w) stop 1
      ]]),[acx_blas_lsame_fcall_ok=yes],
	[acx_blas_lsame_fcall_ok=no])
	AC_MSG_RESULT([$acx_blas_lsame_fcall_ok])
# ISAMAX check (INTEGER return values)
	AC_MSG_CHECKING([whether ISAMAX is called correctly from Fortran])
	AC_RUN_IFELSE(AC_LANG_PROGRAM(,[[
      integer isamax,i
      external isamax
      real a(2)
      a(1) = 1e0
      a(2) = -2e0
      i = isamax(2,a,1)
      if (i.ne.2) stop 1
      ]]),[acx_blas_isamax_fcall_ok=yes],
	[acx_blas_isamax_fcall_ok=no])
	AC_MSG_RESULT([$acx_blas_isamax_fcall_ok])
# SDOT check (REAL return values)
	AC_MSG_CHECKING([whether DDOT is called correctly from Fortran])
	AC_RUN_IFELSE(AC_LANG_PROGRAM(,[[
      real sdot,a(1),b(1),w
      external sdot
      a(1) = 1e0
      b(1) = 2e0
      w = sdot(1,a,1,b,1)
      if (w .ne. a(1)*b(1)) stop 1
      ]]),[acx_blas_sdot_fcall_ok=yes],
	[acx_blas_sdot_fcall_ok=no])
	AC_MSG_RESULT([$acx_blas_sdot_fcall_ok])
# DDOT check (DOUBLE return values)
	AC_MSG_CHECKING([whether DDOT is called correctly from Fortran])
	AC_RUN_IFELSE(AC_LANG_PROGRAM(,[[
      double precision ddot,a(1),b(1),w
      external ddot
      a(1) = 1d0
      b(1) = 2d0
      w = ddot(1,a,1,b,1)
      if (w .ne. a(1)*b(1)) stop 1
      ]]),[acx_blas_ddot_fcall_ok=yes],
	[acx_blas_ddot_fcall_ok=no])
	AC_MSG_RESULT([$acx_blas_ddot_fcall_ok])
# CDOTU check (COMPLEX return values)
	AC_MSG_CHECKING([whether CDOTU is called correctly from Fortran])
	AC_RUN_IFELSE(AC_LANG_PROGRAM(,[[
      complex cdotu,a(1),b(1),w
      external cdotu
      a(1) = cmplx(1e0,1e0)
      b(1) = cmplx(1e0,2e0)
      w = cdotu(1,a,1,b,1)
      if (w .ne. a(1)*b(1)) stop 1
      ]]),[acx_blas_cdotu_fcall_ok=yes],
	[acx_blas_cdotu_fcall_ok=no])
	AC_MSG_RESULT([$acx_blas_cdotu_fcall_ok])
# ZDOTU check (DOUBLE COMPLEX return values)
	AC_MSG_CHECKING([whether ZDOTU is called correctly from Fortran])
	AC_RUN_IFELSE(AC_LANG_PROGRAM(,[[
      double complex zdotu,a(1),b(1),w
      external zdotu
      a(1) = dcmplx(1d0,1d0)
      b(1) = dcmplx(1d0,2d0)
      w = zdotu(1,a,1,b,1)
      if (w .ne. a(1)*b(1)) stop 1
      ]]),[acx_blas_zdotu_fcall_ok=yes],
	[acx_blas_zdotu_fcall_ok=no])
	AC_MSG_RESULT([$acx_blas_zdotu_fcall_ok])

	AC_LANG_POP(Fortran 77)

# if any of the tests failed, reject the BLAS library
	if test $acx_blas_lsame_fcall_ok = yes \
		-a $acx_blas_sdot_fcall_ok = yes \
		-a $acx_blas_ddot_fcall_ok = yes \
		-a $acx_blas_cdotu_fcall_ok = yes \
		-a $acx_blas_zdotu_fcall_ok = yes ; then
		acx_blas_f77_func_ok=yes;
		$1
	else
		acx_blas_f77_func_ok=no;
		$2
	fi
fi

])dnl ACX_BLAS_F77_FUNC

AC_DEFUN([ACX_BLAS_WITH_F77_FUNC], [
AC_PREREQ(2.50)
ACX_BLAS([# disable special action], [])
if test x$acx_blas_ok = xyes ; then
	ACX_BLAS_F77_FUNC(
	[ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])],
	[acx_blas_ok=no; BLAS_LIBS=])
fi
if test x$acx_blas_ok = xno ; then
	$2
fi
])dnl ACX_BLAS_WITH_F77_FUNC
