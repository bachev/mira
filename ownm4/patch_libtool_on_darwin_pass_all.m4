dnl @synopsis PATCH_LIBTOOL_ON_DARWIN_PASS_ALL
dnl
dnl libtool 1.4.x on darwin uses a lib_check with a file_magic that
dnl tests for "Mach-O dynamically linked shared library". However, this
dnl is the file_magic for ".dylib" sharedlibraries but not for ".so"
dnl sharedlibraries. They have another "file -L" result of "Mach-O
dnl bundle ppc", which has an annoying result: when a a module (a .so)
dnl is dependent on another module (another .so) then libtool will
dnl error out and say that the import-module was not found where in
dnl fact it is available. It does not even try to call the real linker.
dnl
dnl Later libtool generations have changed the processing, the import
dnl file_check has been changed from "file_magic" to "pass_all". This
dnl ac-macro does a similar thing: it checks for the darwin host, it
dnl checks for the check_method, and when it was not "pass_all" then we
dnl set it to "deplibs_check_method=pass_all"
dnl
dnl @category Misc
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([PATCH_LIBTOOL_ON_DARWIN_PASS_ALL],
[# libtool-1.4 specific, on darwin set deplibs_check_method=pass_all
case "$host_os" in
  darwin*)
    if grep "^deplibs_check_method=.*file_magic" libtool >/dev/null ; then
AC_MSG_RESULT(patching libtool to set deplibs_check_method=pass_all)
      test -f libtool.old || (mv libtool libtool.old && cp libtool.old libtool)
      sed -e '/^deplibs_check_method=/s/=.*/="pass_all"/' libtool >libtool.new
      (test -s libtool.new || rm libtool.new) 2>/dev/null
      test -f libtool.new && mv libtool.new libtool # not 2>/dev/null !!
      test -f libtool     || mv libtool.old libtool
    fi
  ;;
esac
])
