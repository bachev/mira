dnl @synopsis AC_PROG_CC_NO_WRITEABLE_STRINGS(substvar [,hard])
dnl
dnl @obsoleted Use AX_CFLAGS_NO_WRITABLE_STRINGS.
dnl
dnl Try to find a compiler option that warns when a stringliteral is
dnl used in a place that could potentially modify the address. This
dnl should warn on giving an stringliteral to a function that asks of a
dnl non-const-modified char-pointer.
dnl
dnl The sanity check is done by looking at string.h which has a set of
dnl strcpy definitions that should be defined with const-modifiers to
dnl not emit a warning in all so many places.
dnl
dnl Currently this macro knows about GCC. hopefully will evolve to use:
dnl Solaris C compiler, Digital Unix C compiler, C for AIX Compiler,
dnl HP-UX C compiler, and IRIX C compiler.
dnl
dnl @category Obsolete
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_CC_NO_WRITEABLE_STRINGS], [
  pushdef([CV],ac_cv_prog_cc_no_writeable_strings)dnl
  hard=$2
  if test -z "$hard"; then
    msg="C to warn about writing to stringliterals"
  else
    msg="C to prohibit any write to stringliterals"
  fi
  AC_CACHE_CHECK($msg, CV, [
  cat > conftest.c <<EOF
#include <string.h>
int main (void)
{
   char test[[16]];
   if (strcpy (test, "test")) return 0;
   return 1;
}
EOF
  dnl GCC
  if test "$GCC" = "yes";
  then
        if test -z "$hard"; then
            CV="-Wwrite-strings"
        else
            CV="-fno-writable-strings -Wwrite-strings"
        fi

        if test -n "`${CC-cc} -c $CV conftest.c 2>&1`" ; then
            CV="suppressed: string.h"
        fi

  dnl Solaris C compiler
  elif  $CC -flags 2>&1 | grep "Xc.*strict ANSI C" > /dev/null 2>&1 &&
        $CC -c -xstrconst conftest.c > /dev/null 2>&1 &&
        test -f conftest.o
  then
        # strings go into readonly segment
        CV="-xstrconst"

        rm conftest.o
        if test -n "`${CC-cc} -c $CV conftest.c 2>&1`" ; then
             CV="suppressed: string.h"
        fi

  dnl HP-UX C compiler
  elif  $CC > /dev/null 2>&1 &&
        $CC -c +ESlit conftest.c > /dev/null 2>&1 &&
        test -f conftest.o
  then
       # strings go into readonly segment
        CV="+ESlit"

        rm conftest.o
        if test -n "`${CC-cc} -c $CV conftest.c 2>&1`" ; then
             CV="suppressed: string.h"
        fi

  dnl Digital Unix C compiler
  elif ! $CC > /dev/null 2>&1 &&
        $CC -c -readonly_strings conftest.c > /dev/null 2>&1 &&
        test -f conftest.o
  then
       # strings go into readonly segment
        CV="-readonly_strings"

        rm conftest.o
        if test -n "`${CC-cc} -c $CV conftest.c 2>&1`" ; then
             CV="suppressed: string.h"
        fi

  dnl C for AIX Compiler

  dnl IRIX C compiler
        # -use_readonly_const is the default for IRIX C,
        # puts them into .rodata, but they are copied later.
        # need to be "-G0 -rdatashared" for strictmode but
        # I am not sure what effect that has really.

  fi
  rm -f conftest.*
  ])
  if test -z "[$]$1" ; then
    if test -n "$CV" ; then
      case "$CV" in
        suppressed*) $1="" ;; # known but suppressed
        *)  $1="$CV" ;;
      esac
    fi
  fi
  AC_SUBST($1)
  popdef([CV])dnl
])
