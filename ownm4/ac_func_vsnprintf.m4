dnl @synopsis AC_FUNC_VSNPRINTF
dnl
dnl @obsoleted Use AC_FUNC_SNPRINTF.
dnl
dnl Check whether there is a reasonably sane vsnprintf() function
dnl installed. "Reasonably sane" in this context means never clobbering
dnl memory beyond the buffer supplied, and having a sensible return
dnl value. It is explicitly allowed not to NUL-terminate the return
dnl value, however.
dnl
dnl @category Obsolete
dnl @author Gaute Strokkenes <gs234@cam.ac.uk>
dnl @version 2005-01-25
dnl @license GPLWithACException

AC_DEFUN([AC_FUNC_VSNPRINTF],
[AC_CACHE_CHECK(for working vsnprintf,
  ac_cv_func_vsnprintf,
[AC_TRY_RUN(
[#include <stdio.h>
#include <stdarg.h>

int
doit(char * s, ...)
{
  char buffer[32];
  va_list args;
  int r;

  buffer[5] = 'X';

  va_start(args, s);
  r = vsnprintf(buffer, 5, s, args);
  va_end(args);

  /* -1 is pre-C99, 7 is C99. */

  if (r != -1 && r != 7)
    exit(1);

  /* We deliberately do not care if the result is NUL-terminated or
     not, since this is easy to work around like this.  */

  buffer[4] = 0;

  /* Simple sanity check.  */

  if (strcmp(buffer, "1234"))
    exit(1);

  if (buffer[5] != 'X')
    exit(1);

  exit(0);
}

int
main(void)
{
  doit("1234567");
  exit(1);
}], ac_cv_func_vsnprintf=yes, ac_cv_func_vsnprintf=no, ac_cv_func_vsnprintf=no)])
dnl Note that the default is to be pessimistic in the case of cross compilation.
dnl If you know that the target has a sensible vsnprintf(), you can get around this
dnl by setting ac_func_vsnprintf to yes, as described in the Autoconf manual.
if test $ac_cv_func_vsnprintf = yes; then
  AC_DEFINE(HAVE_WORKING_VSNPRINTF, 1,
            [Define if you have a version of the `vsnprintf' function
             that honours the size argument and has a proper return value.])
fi
])# AC_FUNC_VSNPRINTF
