dnl @synopsis AC_SET_DEFAULT_PATHS_DLLSYSTEM
dnl
dnl @obsoleted AC_SET_DEFAULT_PATHS_SYSTEM is even more intelligent.
dnl
dnl This macro diverts all subpaths to either /bin/.. or /share/..
dnl
dnl @category Obsolete
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([AC_SET_DEFAULT_PATHS_DLLSYSTEM],
[AC_REQUIRE([AC_CANONICAL_HOST])
case ${host_os} in
  *cygwin* | *mingw* | *uwin* | *djgpp | *emx*)
     AC_MSG_RESULT(changing default paths for win/dos target...yes)
     test "$ac_default_prefix" = "/usr/local" && ac_default_prefix="/programs"
     # on win/dos, .exe .dll and .cfg live in the same directory
     bindir=`echo $bindir |sed -e '/^..exec_prefix/s:/bin$:/${PACKAGE}:'`
     libdir=`echo $libdir |sed -e 's:^..exec_prefix./lib$:${bindir}:'`
     sbindir=`echo $sbindir |sed -e 's:^..exec_prefix./sbin$:${bindir}:'`
     sysconfdir=`echo $sysconfdir |sed -e 's:^..prefix./etc$:${bindir}:'`
     libexecdir=`echo $libexecdir |sed -e 's:^..exec_prefix./libexec$:${bindir}/system:'`
     # help-files shall be set with --infodir
     # leave datadir as /share
     infodir=`echo $infodir |sed -e 's:^..prefix./info$:${datadir}/help:'`
     mandir=`echo $mandir |sed -e 's:^..prefix./man$:${datadir}/help:'`
     includedir=`echo $includedir |sed -e 's:..prefix./include$:${datadir}/include:'`
     sharedstatedir=`echo $sharedstatedir |sed -e 's:..prefix./com$:${datadir}/common:'`
     localstatedir=`echo $localstatedir |sed -e 's:..prefix./var$:${datadir}/local:'`
  ;;
esac
])
