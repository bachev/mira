dnl @synopsis AC_PROG_CC_WARNINGS([ANSI])
dnl
dnl @obsoleted Renamed to VL_PROG_CC_WARNINGS, or use AX_CFLAGS_WARN_ALL.
dnl
dnl Enables a reasonable set of warnings for the C compiler.
dnl Optionally, if the first argument is nonempty, turns on flags which
dnl enforce and/or enable proper ANSI C if such are known with the
dnl compiler used.
dnl
dnl Currently this macro knows about GCC, Solaris C compiler, Digital
dnl Unix C compiler, C for AIX Compiler, HP-UX C compiler, IRIX C
dnl compiler, NEC SX-5 (Super-UX 10) C compiler, and Cray J90 (Unicos
dnl 10.0.0.8) C compiler.
dnl
dnl @category Obsolete
dnl @author Ville Laurikari <vl@iki.fi>
dnl @version 2005-01-25
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_CC_WARNINGS], [
  ansi=$1
  if test -z "$ansi"; then
    msg="for C compiler warning flags"
  else
    msg="for C compiler warning and ANSI conformance flags"
  fi
  AC_CACHE_CHECK($msg, ac_cv_prog_cc_warnings, [
    if test -n "$CC"; then
      cat > conftest.c <<EOF
int main(int argc, char **argv) { return 0; }
EOF

      dnl GCC
      if test "$GCC" = "yes"; then
        if test -z "$ansi"; then
          ac_cv_prog_cc_warnings="-Wall"
        else
          ac_cv_prog_cc_warnings="-Wall -ansi -pedantic"
        fi

      dnl Solaris C compiler
      elif $CC -flags 2>&1 | grep "Xc.*strict ANSI C" > /dev/null 2>&1 &&
           $CC -c -v -Xc conftest.c > /dev/null 2>&1 &&
           test -f conftest.o; then
        if test -z "$ansi"; then
          ac_cv_prog_cc_warnings="-v"
        else
          ac_cv_prog_cc_warnings="-v -Xc"
        fi

      dnl HP-UX C compiler
      elif $CC > /dev/null 2>&1 &&
           $CC -c -Aa +w1 conftest.c > /dev/null 2>&1 &&
           test -f conftest.o; then
        if test -z "$ansi"; then
          ac_cv_prog_cc_warnings="+w1"
        else
          ac_cv_prog_cc_warnings="+w1 -Aa"
        fi

      dnl Digital Unix C compiler
      elif $CC -c -verbose -w0 -warnprotos -std1 conftest.c > /dev/null 2>&1 &&
           test -f conftest.o; then
        if test -z "$ansi"; then
          ac_cv_prog_cc_warnings="-verbose -w0 -warnprotos"
        else
          ac_cv_prog_cc_warnings="-verbose -w0 -warnprotos -std1"
        fi

      dnl C for AIX Compiler
      elif $CC 2>&1 | grep AIX > /dev/null 2>&1 &&
           $CC -c -qlanglvl=ansi -qinfo=all conftest.c > /dev/null 2>&1 &&
           test -f conftest.o; then
        if test -z "$ansi"; then
          ac_cv_prog_cc_warnings="-qsrcmsg -qinfo=all:noppt:noppc:noobs:nocnd"
        else
          ac_cv_prog_cc_warnings="-qsrcmsg -qinfo=all:noppt:noppc:noobs:nocnd -qlanglvl=ansi"
        fi

      dnl IRIX C compiler
      elif $CC -c -fullwarn -ansi -ansiE conftest.c > /dev/null 2>&1 &&
           test -f conftest.o; then
        if test -z "$ansi"; then
          ac_cv_prog_cc_warnings="-fullwarn"
        else
          ac_cv_prog_cc_warnings="-fullwarn -ansi -ansiE"
        fi

      dnl The NEC SX-5 (Super-UX 10) C compiler
      elif $CC -c -pvctl[,]fullmsg -Xc conftest.c > /dev/null 2>&1 &&
           test -f conftest.o; then
        if test -z "$ansi"; then
          ac_cv_prog_cc_warnings="-pvctl[,]fullmsg"
        else
          ac_cv_prog_cc_warnings="-pvctl[,]fullmsg -Xc"
        fi

      dnl The Cray J90 (Unicos 10.0.0.8) C compiler
      elif $CC -c -h msglevel 2 conftest.c > /dev/null 2>&1 &&
           test -f conftest.o; then
        if test -z "$ansi"; then
          ac_cv_prog_cc_warnings="-h msglevel 2"
        else
          ac_cv_prog_cc_warnings="-h msglevel 2 -h conform"
        fi

      fi
      rm -f conftest.*
    fi
    if test -n "$ac_cv_prog_cc_warnings"; then
      CFLAGS="$CFLAGS $ac_cv_prog_cc_warnings"
    else
      ac_cv_prog_cc_warnings="unknown"
    fi
  ])
])
