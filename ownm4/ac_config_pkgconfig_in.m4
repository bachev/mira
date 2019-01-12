dnl @synopsis AC_CONFIG_PKGCONFIG_IN [(LIBRARY [, DESCRIPTION [, DESTINATION]])]
dnl
dnl @obsoleted Deprecated with the advent of pkg-config.
dnl
dnl Creates a pkg-config meta-data file for a library with the
dnl following fields: Name, Description, Requires, Conflicts, Version,
dnl Libs, and Cflags.
dnl
dnl This macro automates the construction of a pkg-config .pc meta-data
dnl file; you don't even need to distribute one along. Place this macro
dnl in your configure.ac, et voila, you got one that you want to
dnl install.
dnl
dnl The options:
dnl
dnl  $1 = LIBRARY       e.g. gtk, ncurses
dnl  $2 = DESCRIPTION   one line description of library
dnl  $3 = DESTINATION   directory path
dnl
dnl It is suggested that the following CFLAGS and LIBS variables are
dnl used in your configure.ac. library_libs is *essential*.
dnl library_cflags is important, but not always needed. If they do not
dnl exist, defaults will be taken from LIBRARY_CFLAGS, LIBRARY_LIBS
dnl (should be -llibrary *only*) and LIBRARY_LIBDEPS (-l options for
dnl libraries your library depends upon. LIBLIBRARY_LIBS is simply
dnl $LIBRARY_LIBS $LIBRARY_LIBDEPS. NB. LIBRARY and library are the
dnl name of your library, in upper and lower case repectively e.g. GTK,
dnl gtk.
dnl
dnl  LIBRARY_CFLAGS:    cflags for compiling libraries and example progs
dnl  LIBRARY_LIBS:      libraries for linking programs
dnl  LIBRARY_LIBDEPS*:  libraries for linking libraries against (needed
dnl                     to link -static
dnl  LIBRARY_REQUIRES:  packages required by your library
dnl  LIBRARY_CONFLICTS: packages to conflict with your library
dnl  library_cflags*:   cflags to store in library.pc
dnl  library_libs*:     libs to store in library.pc
dnl  LIBLIBRARY_LIBS:   libs to link programs IN THIS PACKAGE ONLY against
dnl  LIBRARY_VERSION*:  the version of your library (x.y.z recommended)
dnl    *=required if you want sensible output, otherwise they will be
dnl      *guessed* (DWIM, but usually correct)
dnl
dnl There is also an AC_SUBST(LIBRARY_PKGCONFIG) that will be set to
dnl the name of the meta-data file. Use as:
dnl
dnl  install-data-local: install-pkgconfig
dnl  install-pkgconfig:
dnl     $(mkinstalldirs) $(DESTDIR)$(prefix)/lib/pkgconfig
dnl     $(INSTALL_DATA) @LIBRARY_PKGCONFIG@ $(DESTDIR)$(prefix)/lib/pkgconfig
dnl
dnl Or, if using automake:
dnl
dnl  pkgconfigdatadir = $(prefix)/lib/pkgconfig
dnl  pkgconfigdata_DATA = @LIBRARY_PKGCONFIG@
dnl
dnl Example usage:
dnl
dnl  GIMPPPRINT_LIBS="-lgimpprint"
dnl  AC_CHECK_LIB(m,pow,
dnl               GIMPPRINT_DEPLIBS="${GIMPPRINT_DEPLIBS} -lm")
dnl  AC_CONFIG_PKGCONFIG_IN([gimpprint], [GIMP Print Top Quality Printer Drivers], [src/main])
dnl  AC_CONFIG_FILES(gimpprint.pc)
dnl
dnl  If using automake with the "--add-missing" option, the call to
dnl  AC_CONFIG_FILES will cause it to abort because it is looking for the
dnl  .pc.in file which has not yet been created (because configure itself
dnl  hasn't been run). To get around this, don't use AC_CONFIG_FILES to
dnl  generate it.  Instead, use this in Makefile.am:
dnl
dnl  @LIBRARY_PKGCONFIG@ : @LIBRARY_PKGCONFIG@.in
dnl  		./config-status --file=@LIBRARY_PKGCONFIG@
dnl
dnl @category Obsolete
dnl @author Roger Leigh <roger@whinlatter.uklinux.net>
dnl @author Diab Jerius <djerius@cfa.harvard.edu>
dnl @version 2006-05-27
dnl @license AllPermissive

## AC_CONFIG_PKGCONFIG_IN(LIBRARY, DESCRIPTION, DESTINATION)
## ---------------------------------------------------------
## Create a pkg-config meta-data file for LIBRARY.  Include a one-line
## DESCRIPTION.  The meta-data file will be created in a DESTINATION
## directory.
AC_DEFUN([AC_CONFIG_PKGCONFIG_IN],
[# create a pkg-config meta-data file ($1.pc.in)
m4_pushdef([PKGCONFIG_DIR], [m4_if([$3], , , [$3/])])
PKGCONFIG_FILE_PC="PKGCONFIG_DIR[]$1.pc"
PKGCONFIG_FILE="$PKGCONFIG_FILE_PC[].in"
AC_SUBST(target)dnl
AC_SUBST(host)dnl
AC_SUBST(build)dnl
dnl create directory if it does not preexist
m4_if([$3], , , [AS_MKDIR_P([$3])])
AC_MSG_NOTICE([creating $PKGCONFIG_FILE])
dnl we're going to need uppercase, lowercase and user-friendly versions of the
dnl string `MODULE'
m4_pushdef([MODULE_UP], m4_translit([$1], [a-z], [A-Z]))dnl
m4_pushdef([MODULE_DOWN], m4_translit([$1], [A-Z], [a-z]))dnl
if test -z "$MODULE_DOWN[]_cflags" ; then
  if test -n "$MODULE_UP[]_CFLAGS" ; then
      MODULE_DOWN[]_cflags="$MODULE_UP[]_CFLAGS"
  else
dnl    AC_MSG_WARN([variable `MODULE_DOWN[]_cflags' undefined])
    MODULE_DOWN[]_cflags=''
  fi
fi
AC_SUBST(MODULE_DOWN[]_cflags)dnl
if test -z "$MODULE_DOWN[]_libs" ; then
  if test -n "$MODULE_UP[]_LIBS" ; then
    MODULE_DOWN[]_libs="$MODULE_UP[]_LIBS"
  else
    AC_MSG_WARN([variable `MODULE_DOWN[]_libs' and `MODULE_UP[]_LIBS' undefined])
    MODULE_DOWN[]_libs='-l$1'
  fi
  if test -n "$MODULE_UP[]_LIBDEPS" ; then
    MODULE_DOWN[]_libs="$MODULE_DOWN[]_libs $MODULE_UP[]_LIBDEPS"
  fi
fi
AC_SUBST(MODULE_DOWN[]_libs)dnl
AC_SUBST(MODULE_UP[]_REQUIRES)
AC_SUBST(MODULE_UP[]_CONFLICTS)
if test -z "$MODULE_UP[]_VERSION" ; then
  AC_MSG_WARN([variable `MODULE_UP[]_VERSION' undefined])
  MODULE_UP[]_VERSION="$VERSION"
fi
AC_SUBST(MODULE_UP[]_VERSION)dnl
echo 'prefix=@prefix@' >$PKGCONFIG_FILE
echo 'exec_prefix=@exec_prefix@' >>$PKGCONFIG_FILE
echo 'libdir=@libdir@' >>$PKGCONFIG_FILE
echo 'includedir=@includedir@' >>$PKGCONFIG_FILE
echo '' >>$PKGCONFIG_FILE
echo 'Name: @PACKAGE_NAME@' >>$PKGCONFIG_FILE
echo 'Description: $2' >>$PKGCONFIG_FILE
if test -n "$MODULE_UP[]_REQUIRES" ; then
  echo 'Requires: @MODULE_UP[]_REQUIRES@' >>$PKGCONFIG_FILE
fi
if test -n "$MODULE_UP[]_CONFLICTS" ; then
  echo 'Conflicts: @MODULE_UP[]_CONFLICTS@' >>$PKGCONFIG_FILE
fi
echo 'Version: @MODULE_UP[]_VERSION@' >>$PKGCONFIG_FILE
echo 'Libs: -L${libdir} @MODULE_DOWN[]_libs@' >>$PKGCONFIG_FILE
echo 'Cflags: -I${includedir} @MODULE_DOWN[]_cflags@' >>$PKGCONFIG_FILE
m4_pushdef([PKGCONFIG_UP], [m4_translit([$1], [a-z], [A-Z])])dnl
PKGCONFIG_UP[]_PKGCONFIG="$PKGCONFIG_FILE_PC"
AC_SUBST(PKGCONFIG_UP[]_PKGCONFIG)
m4_popdef([PKGCONFIG_UP])
m4_popdef([MODULE_DOWN])dnl
m4_popdef([MODULE_UP])dnl
m4_popdef([PKGCONFIG_DIR])dnl
])
