dnl @synopsis AC_PROG_CC_STRICT_PROTOTYPES(substvar [,hard])
dnl
dnl @obsoleted Use AX_CFLAGS_STRICT_PROTOTYPES.
dnl
dnl Try to find a compiler option that warns when a function prototype
dnl is not fully defined. Enable it only if the the system headers are
dnl reasonably clean with respect to compiling with strict-prototypes.
dnl
dnl The sanity check is done by looking at sys/signal.h which has a set
dnl of macro-definitions SIG_DFL and SIG_IGN that are cast to the local
dnl signal-handler type. If that signal-handler type is not fully
dnl qualified then the system headers are not seen as strictly
dnl prototype clean.
dnl
dnl Currently this macro knows about GCC. hopefully will evolve to use:
dnl Solaris C compiler, Digital Unix C compiler, C for AIX Compiler,
dnl HP-UX C compiler, and IRIX C compiler.
dnl
dnl @category Obsolete
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_CC_STRICT_PROTOTYPES], [
  pushdef([CV], ac_cv_prog_cc_strict_prototypes)dnl
  hard=$2
  if test -z "$hard"; then
    msg="C to warn at nonstrict prototypes"
  else
    msg="C to require strict prototypes"
  fi
  AC_CACHE_CHECK($msg, CV, [
  cat > conftest.c <<EOF
#include <sys/signal.h>
int main (void)
{
   if (signal (SIGINT, SIG_IGN) == SIG_DFL) return 0;
   return 1;
}
EOF

  dnl GCC
  if test "$GCC" = "yes"; then
    if test -z "$hard"; then
      CV="-Wstrict-prototypes"
    else
      CV="-fstrict-prototypes -Wstrict-prototypes"
    fi

    if test -n "`${CC-cc} -c $CV conftest.c 2>&1`" ; then
      CV="suppressed...sys/stat.h"
    fi

  dnl Solaris C compiler

  dnl HP-UX C compiler

  dnl Digital Unix C compiler

  dnl C for AIX Compiler

  dnl IRIX C compiler

  fi
  rm -f conftest.*
  ])
  if test -z "[$]$1" ; then
    if test -n "$CV" ; then
      case "$CV" in
        *...*) $1="" ;; # known but suppressed
        *)  $1="$CV" ;;
      esac
    fi
  fi
  AC_SUBST($1)
  popdef([CV])dnl
])
