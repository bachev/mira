dnl @synopsis AC_SET_RELEASEINFO_VERSIONINFO [(VERSION)]
dnl
dnl   default $1 = $VERSION
dnl
dnl check the $VERSION number and cut the two last digit-sequences off
dnl which will form a -version-info in a @VERSIONINFO@ ac_subst while
dnl the rest is going to the -release name in a @RELEASEINFO@ ac_subst.
dnl
dnl you should keep these two seperate - the release-name may contain
dnl alpha-characters and can be modified later with extra release-hints
dnl e.g. RELEASEINFO="$RELEASINFO-debug" for a debug version of your
dnl lib.
dnl
dnl example: a VERSION="2.4.18" will be transformed into "-release 2
dnl -version-info 4:18" and for a linux-target this will tell libtool
dnl to install the lib as "libmy.so libmy.la libmy.a libmy-2.so.4
dnl libmy-2.so.4.0.18" and executables will get link-resolve-infos for
dnl libmy-2.so.4 - therefore the patch-level is ignored during ldso
dnl linking, and ldso will use the one with the highest patchlevel.
dnl Using just "-release $(VERSION)" during libtool-linking would not
dnl do that - omitting the -version-info will libtool install libmy.so
dnl libmy.la libmy.a libmy-2.4.18.so and executables would get
dnl hardlinked with the 2.4.18 version of your lib.
dnl
dnl This background does also explain the default dll name for a win32
dnl target : libtool will choose to make up libmy-2-4.dll for this
dnl version spec.
dnl
dnl this macro does set the three parts
dnl VERSION_REL.VERSION_REQ.VERSION_REL from the VERSION-spec but does
dnl not ac_subst them like the two INFOs. If you prefer a two-part
dnl VERSION-spec, the VERSION_REL will still be set, either to the
dnl host_cpu or just a simple "00". You may add sublevel parts like
dnl "1.4.2-ac5" where the sublevel is just killed from these
dnl versioninfo/releasinfo substvars.
dnl
dnl @category Misc
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([AC_SET_RELEASEINFO_VERSIONINFO],
[# ------ AC SET RELEASEINFO VERSIONINFO --------------------------------
AC_MSG_CHECKING(version info)
  VERSION_REQ=`echo ifelse( $1, , $VERSION, $1 )` # VERSION_TMP really...
  VERSION_REL=`echo $VERSION_REQ | sed -e 's/[[.]][[^.]]*$//'`  # delete micro
  VERSION_REV=`echo $VERSION_REQ | sed -e "s/^$VERSION_REL.//"` # the rest
  VERSION_REQ=`echo $VERSION_REL | sed -e 's/.*[[.]]//'`  # delete prefix now
  VERSION_REV=`echo $VERSION_REV | sed -e 's/[[^0-9]].*//'` # 5-p4 -> 5
  if test "$VERSION_REQ" != "$VERSION_REL" ; then # three-part version...
  VERSION_REL=`echo $VERSION_REL | sed -e "s/.$VERSION_REQ\$//"`
  else # or has been two-part version - try using host_cpu if available
  VERSION_REL="00" ; test "_$host_cpu" != "_" && VERSION_REL="$host_cpu"
  fi
  RELEASEINFO="-release $VERSION_REL"
  VERSIONINFO="-version-info $VERSION_REQ:$VERSION_REV"
AC_MSG_RESULT([$RELEASEINFO $VERSIONINFO])
AC_SUBST([RELEASEINFO])
AC_SUBST([VERSIONINFO])
])
