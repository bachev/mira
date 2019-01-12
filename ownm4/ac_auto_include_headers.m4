dnl @synopsis AC_AUTO_INCLUDE_HEADERS (HEADER-FILE, INCLUDE-FILE ...)
dnl
dnl Given a HEADER-FILE and a space-separated list of INCLUDE-FILEs,
dnl AC_AUTO_INCLUDE_HEADERS will append to HEADER-FILE a conditional
dnl #include for each INCLUDE-FILE. For instance, the following macro
dnl call:
dnl
dnl    AC_AUTO_INCLUDE_HEADERS([config-inc.h], [sys/foobar.h])
dnl
dnl will append the following text to config-inc.h:
dnl
dnl    #ifdef HAVE_SYS_FOOBAR_H
dnl    # include <sys/foobar.h>
dnl    #endif
dnl
dnl AC_AUTO_INCLUDE_HEADERS makes it easy to auto-generate a single
dnl header file that can then be #include'd by multiple files in a
dnl project. Because the #ifdef's are appended to HEADER-FILE, it's
dnl also convenient to include additional text in that file. For
dnl instance:
dnl
dnl    cat <<\CIH_EOF > config-inc.h
dnl    /* This file was generated automatically by configure. */
dnl
dnl    #ifndef _CONFIG_INC_H_
dnl    #define _CONFIG_INC_H_
dnl
dnl    #include <stdio.h>
dnl
dnl    CIH_EOF
dnl    AC_AUTO_INCLUDE_HEADERS([config-inc.h], [arpa/inet.h dlfcn.h errno.h])
dnl    echo "#endif" >> config-inc.h
dnl
dnl Here's an easy way to get a complete list of header files from
dnl config.h:
dnl
dnl    cat config.h | perl -ane '/ HAVE_\S+_H / && do {$_=$F[$#F-1]; s/^HAVE_//; s/_H/.h/; s|_|/|g; tr/A-Z/a-z/; print "$_ "}'
dnl
dnl You can then manually edit the resulting list.
dnl
dnl @category InstalledPackages
dnl @author Scott Pakin <pakin@uiuc.edu>
dnl @version 2002-03-04
dnl @license AllPermissive

AC_DEFUN([AC_AUTO_INCLUDE_HEADERS],
[touch $1
for ac_auto_include_file in $2; do
  ac_auto_include_have=`echo $ac_auto_include_file | sed 'y%./*abcdefghijklmnopqrstuvwxyz%__PABCDEFGHIJKLMNOPQRSTUVWXYZ%;s%[^_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789]%_%g'`
  echo "#ifdef HAVE_$ac_auto_include_have" >> $1
  echo "# include <$ac_auto_include_file>" >> $1
  echo "#endif" >> $1
  echo "" >> $1
done])
