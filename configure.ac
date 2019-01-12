# Process this file with autoconf to produce a configure script.

#LT_PREREQ([2.4.2])

AC_INIT(mira, m4_esyscmd([build-aux/git-version-gen isshort .tarball-version]))
AC_CONFIG_MACRO_DIR([m4])
#AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/stdinc/defines.H)
AM_INIT_AUTOMAKE([dist-bzip2 tar-ustar])

#AM_CONFIG_HEADER(config.h)

AM_CONDITIONAL([ISDEBIAN], [test x$DEBIAN_MAINTAINER != x])

systemuname=`uname`

# The internal echo for /bin/sh on darwin doesn't support -n
# therefore we have to use /bin/echo.
if test "${systemuname}" = "Darwin"; then
    ECHO_N="/bin/echo -n";
else
    ECHO_N="echo -n";
fi

# we'll need that quite soon
AC_PROG_SED

AC_ARG_ENABLE(mirastatic,
	AC_HELP_STRING([--enable-mirastatic],
	[enables building of static binaries, default is no]))

# let's make life simple if user has everything in standard homebrew
AC_ARG_WITH(brew,
	AC_HELP_STRING([--with-brew],
	[if you have build system installed via homebrew]))

AC_ARG_WITH(compiler,
	AC_HELP_STRING([--with-compiler],
	[either gcc or clang]))


BREWBASE="brewbase-not-set"
if test "x${with_brew}" = "xyes" ; then
  AC_CHECK_PROG([HAVE_BREW], [brew],[yes],[no])
  if test "${HAVE_XXD}" = "no"; then
    echo
    echo "Could not find the program 'brew', but you wanted brew. That's odd."
    echo "Aborting."
  fi

  BREWBASE=`brew --prefix`
fi

extradyldpath=""
if test "x${enable_mirastatic}" = "xyes" ; then
  AC_SUBST(MIRASTATIC,"static")
  AC_SUBST(MIRASTATICLD,"-all-static")
  extradyldpath=${BREWBASE}/lib:
else
  enable_mirastatic="no"
  AC_SUBST(MIRASTATIC,"dynamic")
  AC_SUBST(MIRASTATICLD,"")
fi

AC_SUBST(EXTRADYLDPATH,${extradyldpath})

# ------------------------------ C / C++ compiler and flags
# get the C/C++ compiler flags
CFLAGS="${CFLAGS} "
CXXFLAGS="${CXXFLAGS} "

AMCXXFLAGS=""

bachCCOMPILER=""
bachCCsource=""

case "${CC}" in
*gcc*)
  bachCCOMPILER=gcc
  bachCCsource="CC env variable"
;;
*clang*)
  bachCCOMPILER=clang
  bachCCsource="CC env variable"
;;
esac

if test "x${with_compiler}" != "x" ; then
  case "${with_compiler}" in
  *gcc*)
    bachCCOMPILER=gcc
    bachCCsource="--with-compiler flag"
  ;;
  *clang*)
    bachCCOMPILER=clang
    bachCCsource="--with-compiler flag"
  ;;
  *)
    echo "Unknown compiler given via --with-compiler"
    exit 1
  ;;
  esac
fi

if test "x${bachCCOMPILER}" = "x" ; then
  bachCCOMPILER="gcc"
  if test "${systemuname}" = "Darwin"; then
    bachCCOMPILER="clang"
  fi
  bachCCsource="default setting"
fi

AC_MSG_CHECKING(compiler family)
AC_MSG_RESULT(${bachCCOMPILER})
AC_MSG_CHECKING(compiler determined by)
AC_MSG_RESULT(${bachCCsource})

case "${bachCCOMPILER}" in
gcc)
  ccname="gcc"
  cxxname="g++"
  ccminversion=6.1
;;
clang)
  ccname="clang"
  cxxname="clang++"
  ccminversion=3.5
;;
*)
  echo "Could not find out what compiler we're using?"
  exit 1
;;
esac

# if using brew, try to find out correct compiler binary name
#  and set CC / CXX if they were not set
if test "x${with_brew}" = "xyes" ; then
  case "${bachCCOMPILER}" in
  gcc)
    bDIR=${BREWBASE}/opt/gcc/bin
  ;;
  clang)
    bDIR=${BREWBASE}/opt/llvm/bin
  ;;
  *)
    echo "The compiler is not GCC nor clang, no idea how to set that up with brew."
    exit 100;
  ;;
  esac
  ename="${bDIR}/${ccname}"
  AC_CHECK_FILE([$ename],[CCexe=yes],[CCexe=no])
  # if not found, check for versions with "dash versionnumber"
  if test "${CCexe}" = "no"; then
    AC_MSG_CHECKING(checking for compiler with version number as executable)
    tname=`ls -1 ${ename}-[[0-9]]* | head -1`
    if test "x$tname" = "x"; then
      echo "Did not find executable ${ccname}, nor variations with -number."
      echo "Aborting."
      exit 1
    else
      # will drop begin of string up to last occur of `-`
      #  via https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash
      #  (and is not restricted to bash!)
      vnum=${tname##*-}
      ccname="${ccname}-$vnum"
      cxxname="${cxxname}-$vnum"
    fi
    AC_MSG_RESULT(${vnum})
  fi
fi

AC_PROG_CC([$ccname])
AC_PROG_CXX([$cxxname])

# circumvent a bug (don't remember which and where it strikes)
if test "$bachCCOMPILER" != "gcc" ; then
  # No, autoconf, no! Bad autoconf!
  # We do not have GCC! I'm telling you!
  # (though, to be fair, it's not entirely your fault)
  GCC="no"
fi

# get compiler version
if test $bachCCOMPILER = gcc; then
  ccversion=`$CC -dumpversion`
  cxxversion=`$CXX -dumpversion`
elif test $bachCCOMPILER = clang; then
  ccversion=`$CC --version | head -1 | cut -f 3 -d ' '`
  cxxversion=`$CXX --version | head -1 | cut -f 3 -d ' '`
else
  echo "Need to find out how to get version string for ${bachCCOMPILER}, contact author."
  exit 1
fi

if test "x${ccversion}" != "x${cxxversion}"; then
  echo "${CC} is version ${ccversion} while ${CXX} is version ${cxxversion}"
  echo "This is not a sane environment, both versions should be the same (>=${ccminversion}). Stopping here."
  exit 10
fi

ccgoodversion="no"
AX_COMPARE_VERSION([$ccversion],[ge],[$ccminversion],[ccgoodversion="yes"])

if test "${ccgoodversion}" = "no"; then
  echo
  echo "When using ${bachCCOMPILER}, you need to have a version >= ${ccminversion}. You have only ${ccversion}, aborting."
  exit 10
fi

# OK, MIRA wants C++14, makes life a whole lot easier
# Also, do this very early in the process as some compiler tests need those
#  flags already set correctly

case "${bachCCOMPILER}" in
gcc)
  CXXFLAGS="${CXXFLAGS} -std=c++14"
;;
clang)
  CXXFLAGS="${CXXFLAGS} -std=c++14 -stdlib=libc++"
;;
*)
  echo "The compiler is not GCC nor clang and I need a flag to enable at least C++14 standard, please contact the author!"
  exit 100;
;;
esac

echo
echo "+++++++++++++++++++++++++++++++++++++"
echo "systemuname: ${systemuname}"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "+++++++++++++++++++++++++++++++++++++"

if test "x${with_brew}" = "xyes" ; then
  case "${bachCCOMPILER}" in
  gcc)
    tDIR=${BREWBASE}/opt/gcc
  ;;
  clang)
    tDIR=${BREWBASE}/opt/llvm
  ;;
  *)
    echo "The compiler is not GCC nor clang, no idea how to set that up with brew."
    exit 100;
  ;;
  esac
  CFLAGS="$CFLAGS -I${tDIR}/include"
  CXXFLAGS="$CXXFLAGS -I${tDIR}/include"
  LDFLAGS="$LDFLAGS -L${tDIR}/lib -Wl,-rpath,${tDIR}/lib"

  tDIR=${BREWBASE}/opt/flex
  CFLAGS="$CFLAGS -I${tDIR}/include"
  CXXFLAGS="$CXXFLAGS -I${tDIR}/include"
  LDFLAGS="$LDFLAGS -L${tDIR}/lib"

  if test "x${with_boost}" = "x" ; then
    with_boost=${BREWBASE}/opt/boost
  fi

  echo "After brew flag setting:"
  echo "+++++++++++++++++++++++++++++++++++++"
  echo "CFLAGS now: ${CFLAGS}"
  echo "CXXFLAGS now: ${CXXFLAGS}"
  echo "LDFLAGS now: ${LDFLAGS}"
  echo "+++++++++++++++++++++++++++++++++++++"
fi

#
# Seen on systems where third party programs/libs were installed in a
#  separate non-standard directory like "/home/foo/public":
#
# The programs may be found if they're in the path, but the configure
#  script will have no idea it needs to use /home/foo/public/include as
#  additional include for the C/C++ compiler and the lib directory for
#  linking.
#
# Therefore, do that manually

AC_MSG_CHECKING(whether to add --prefix as include path)
haveuserprefix="no"
if test "x${prefix}" != "xNONE" ; then
  CFLAGS="$CFLAGS -I${prefix}/include"
  CXXFLAGS="$CXXFLAGS -I${prefix}/include"
  LDFLAGS="$LDFLAGS -L${prefix}/lib"
  haveuserprefix="yes"
fi
AC_MSG_RESULT(${haveuserprefix})

echo
echo "+++++++++++++++++++++++++++++++++++++"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "+++++++++++++++++++++++++++++++++++++"


# we want OpenMP
AX_OPENMP
have_openmp=yes
if test "x${OPENMP_CFLAGS}" = "x"; then
  have_openmp=no
fi
# *sigh* Some autoconf versions apparently set OPENMP c flags, but not c++
if test "x${OPENMP_CFLAGS}" != "x${OPENMP_CXXFLAGS}"; then
  if test "x${OPENMP_CXXFLAGS}" = "x"; then
    echo "Rescuing C++ flags for OpenMP"
    OPENMP_CXXFLAGS=$OPENMP_CFLAGS
  fi
fi
echo "+++++++++++++++++++++++++++++++++++++"
echo "OpenMP CFLAGS: ${OPENMP_CFLAGS}"
echo "OpenMP CXXFLAGS: ${OPENMP_CXXFLAGS}"
echo "+++++++++++++++++++++++++++++++++++++"


# checks for BOOST and libraries behind the checks for 64 bit as some
#  rescue values and compiler flags depend on the 64 bit enabled/not
#  enabled result

# ---------------------------------------------------------------------------
# adapted from EMBOSS configure.in
# ---------------------------------------------------------------------------

echo
echo "+++++++++++++++++++++++++++++++++++++"
echo "systemuname: ${systemuname}"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "+++++++++++++++++++++++++++++++++++++"

AC_ARG_ENABLE(64,
	AC_HELP_STRING([--enable-64],
	[enables 64 bit compiling, use 'no' for 32 bit on 64 bit machines, default is dependend of your compiler defaults]))

# check whether compiler handles 64 bit types
AC_CHECK_TYPE(int64_t,has_int64_t="yes",has_int64_t="no")

# silent check for 64 bit data type
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
     #include <sys/types.h>
     #include <limits.h>
]], [[
     #if LONG_MAX!=0x7fffffffffffffff
     #error The long definition is other than 8 byte. Try setting right compiler options.
     #endif
]])],[
  silenthas_int64_t="yes"
],[
  silenthas_int64_t="no"
])

# if 64 bit was not forbidden or silent check for 64 succeeded
if test "${enable_64}" != "no" || test "${silenthas_int64_t}" = "yes" ; then

  # Test for Linux 64 bit
  if test "${systemuname}" = "Linux"; then
    CXXFLAGS="-DAJ_Linux64 $CXXFLAGS"
  fi

  # Test for FreeBSD 64 bit
  if test "${systemuname}" = "FreeBSD"; then
    CXXFLAGS="-DAJ_FreeBSD64 $CXXFLAGS"
  fi

  # Test for OSF1 64 bit (its always this anyway)
  if test "${systemuname}" = "OSF1"; then
    CXXFLAGS="-DAJ_OSF164 $CXXFLAGS"
  fi

  # Test for SGI 64 bit
  if test "${systemuname}" = "IRIX" || test "${systemuname}" = "IRIX64"; then
    CXXFLAGS="-DAJ_IRIX64 $CXXFLAGS"
  fi

  # Test for AIX 64 bit
  if test "${systemuname}" = "AIX"; then
    CXXFLAGS="-DAJ_AIX64 $CXXFLAGS"
    if test $bachCCOMPILER != gcc; then
      CC="$CC -q64"
    fi
    NM="nm -B -X 64"
    AR="ar -X 64"
  fi

  # Test for HP-UX 64 bit
  if test "${systemuname}" = "HP-UX"; then
    if test $bachCCOMPILER != gcc; then
      CC="$CC +DD64"
    fi
    CXXFLAGS="-DAJ_HPUX64 $CXXFLAGS"
  fi

  # Test for Solaris 64 bit
  if test "${systemuname}" = "SunOS"; then
    if test $bachCCOMPILER = gcc; then
      CFLAGS="$CFLAGS -m64"
      CXXFLAGS="$CXXFLAGS -m64"
      LDFLAGS="$LDFLAGS -m64"
    fi
  fi

fi

echo
echo "+++++++++++++++++++++++++++++++++++++"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "+++++++++++++++++++++++++++++++++++++"

# 'official' check for 64 bit data type
AC_MSG_CHECKING(whether long has 64bit)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
     #include <sys/types.h>
     #include <limits.h>
]], [[
     #if LONG_MAX!=0x7fffffffffffffffL
     #error The long definition is other than 8 byte. Try setting right compiler options.
     #endif
]])],[
  has_int64_t="yes"
  AC_MSG_RESULT(yes)
],[
  has_int64_t="no"
  AC_MSG_RESULT(no)
])

## Test if --enable-64 is given
## If not, set enable_64 to "yes" for when the CPU handles 64 bit types,
##  else set to "no"
AC_MSG_CHECKING(whether enabling 64bit compilation support)
if test "${enable_64+set}" != "set"; then
   enable_64=${has_int64_t}
fi
AC_MSG_RESULT(${enable_64})

AC_SUBST(BITPACKAGE,${build_cpu})
AC_DEFINE_UNQUOTED(ENABLE64)

# circumvent the botched linker with explicit template instantiation on OSX gcc
if test $bachCCOMPILER = gcc; then
  if test "${systemuname}" = "Darwin"; then
    AC_DEFINE_UNQUOTED(BOTCHEDEXPLICITINSTANTIATIONLINKER)
  fi
fi


# having this before AC_CANONICAL_BUILD sets the correct CPU info in ${build_os}
#
# and we need 10.6 instead of 10.5 as for 10.5, the linker creates a binary like
# this
#  /usr/lib/libgcc_s.1.dylib (compatibility version 1.0.0, current version 1669.0.0)
#  /opt/localwgcc48/lib/libgcc_s.1.dylib (compatibility version 1.0.0, current version 1.0.0)
#  /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 169.3.0)
#
# which does not run on a couple of other OSX installations
#
# with clang as compiler, we even need 10.7

if test "${systemuname}" = "Darwin"; then
  case "${CC}" in
    *clang*)
      CFLAGS="-mmacosx-version-min=10.7 ${CFLAGS}"
      CXXFLAGS="-mmacosx-version-min=10.7 ${CXXFLAGS}"
    ;;
    *)
      CFLAGS="-mmacosx-version-min=10.6 ${CFLAGS}"
      CXXFLAGS="-mmacosx-version-min=10.6 ${CXXFLAGS}"
    ;;
  esac
fi

AC_CANONICAL_BUILD
#AC_CANONICAL_TARGET
#AC_CANONICAL_HOST

# AC_CANONICAL_BUILD on Darwin needed the -arch flag set above to correctly
#  set build_cpu. So only now we can set the BITPACKAGE

AC_SUBST(BITPACKAGE,${build_cpu})


# ---------------------------------------------------------------------------
# Set optimisation flags
# ---------------------------------------------------------------------------
if test "${enable_optimisations+set}" != "set"; then
  enable_optimisations="yes"
fi
AC_ARG_ENABLE(optimisations,
	AC_HELP_STRING([--enable-optimisations],
	[enables compiler optimisations, default is yes]))

AC_ARG_ENABLE(native,
	AC_HELP_STRING([--enable-native],
	[compile only for your specific processor, default is no]))

if test "${enable_native+set}" != "set"; then
  enable_native="no"
fi


if test "${enable_optimisations}" = "yes" ; then
  AC_MSG_CHECKING(for compiler optimisation settings)

  if test $bachCCOMPILER = gcc; then
    # For GCC, be conservative for all plattforms except the ones where
    #  -O3 is known to work

    OPTFLAGS="-O2"
    if test "${systemuname}" = "Linux"; then
      OPTFLAGS="-O3 -funroll-loops"
    fi

    if test "${systemuname}" = "Darwin"; then
      OPTFLAGS="-O3 -funroll-loops"
    fi

    # native?
    if test "${enable_native}" = "yes" ; then
      OPTFLAGS+=" -march=native -mtune=native"
#    else
#      # tuning for the compiler we currently build on is not a too bad idea
#      # or is it?
#      OPTFLAGS+=" -mtune=native"
    fi

  else

    # default
    OPTFLAGS="-O2"

    # Test for AIX optimisation
    if test "${systemuname}" = "AIX"; then
      OPTFLAGS="-O3 -qstrict -qarch=auto -qtune=auto"
    fi

    # test for SGI optimisation
    if test "${systemuname}" = "IRIX" || test "${systemuname}" = "IRIX64" ; then
      LD="/usr/bin/ld -IPA"
      OPTFLAGS="-O3"
    fi

    # Test for HP-UX optimisation
    if test "${systemuname}" = "HP-UX"; then
      OPTFLAGS="-fast"
    fi

    # Test for OSF1 optimisation
    if test "${systemuname}" = "OSF1"; then
    OPTFLAGS="-fast -U_FASTMATH"
    fi

    # test for Solaris optimisation
    if test "${systemuname}" = "SunOS"; then
      OPTFLAGS="-O"
    # test for 64 bit ptr here (see Solaris 64bit above)
      if test "${enable_64}" = "yes" ; then
        OPTFLAGS="-xtarget=ultra -xarch=v9"
      fi
    fi

    # Test for Linux optimisation
    if test "${systemuname}" = "Linux"; then
      OPTFLAGS="-O2"
    fi

    # Test for FreeBSD optimisation
    if test "${systemuname}" = "FreeBSD"; then
      OPTFLAGS="-O2"
    fi

    # Intel MacOSX requires reduced optimisation for PCRE code
    if test "${systemuname}" = "Darwin"; then
      if test "`uname -a | grep i386`"; then
        OPTFLAGS="-O1"
      else
        OPTFLAGS="-O2"
      fi
    fi
  fi

  CFLAGS="${CFLAGS} ${OPTFLAGS}"
  CXXFLAGS="${CXXFLAGS} ${OPTFLAGS}"

  AC_MSG_RESULT(${OPTFLAGS})
fi
# End of optimisation tests
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Set special system/compiler combo flags
# ---------------------------------------------------------------------------

if test $bachCCOMPILER != gcc; then
   COMBOFLAGS=""

  # Test for AIX optimisation
  if test "${systemuname}" = "AIX"; then
   COMBOFLAGS=""
  fi

  # test for SGI optimisation
  if test "${systemuname}" = "IRIX" || test "${systemuname}" = "IRIX64" ; then
   COMBOFLAGS=""
  fi

  # Test for HP-UX optimisation
  if test "${systemuname}" = "HP-UX"; then
    COMBOFLAGS="-Aa +nostl"
  fi

  # Test for OSF1 optimisation
  if test "${systemuname}" = "OSF1"; then
   COMBOFLAGS=""
  fi

  # test for Solaris optimisation
  if test "${systemuname}" = "SunOS"; then
   COMBOFLAGS=""
  fi

  # Test for Linux optimisation
  if test "${systemuname}" = "Linux"; then
   COMBOFLAGS=""
  fi

  # Test for FreeBSD optimisation
  if test "${systemuname}" = "FreeBSD"; then
   COMBOFLAGS=""
  fi

  # Intel MacOSX requires reduced optimisation for PCRE code
  if test "${systemuname}" = "Darwin"; then
    COMBOFLAGS=""
  fi

  CFLAGS="${CFLAGS} ${COMBOFLAGS}"
  CXXFLAGS="${CXXFLAGS} ${COMBOFLAGS}"
fi

# End of special combo flags
# ---------------------------------------------------------------------------



# Add extensions to Solaris for some reentrant functions
if test "${systemuname}" = "SunOS"; then
  CFLAGS="$CFLAGS -D__EXTENSIONS__"
fi

# Building on OSX is ... difficult.
BUNDLETARGET=nobundle
if test "${systemuname}" = "Darwin"; then
  # this is a solution from http://hackage.haskell.org/trac/ghc/ticket/5019
  # to the exception handling problem of gcc and the changed OSX functionality
  if test $bachCCOMPILER = gcc; then
    LDFLAGS="-Wl,-no_compact_unwind ${LDFLAGS}"
  fi

  BUNDLETARGET=darwin
fi
AC_SUBST(BUNDLETARGET,${BUNDLETARGET})


echo
echo "+++++++++++++++++++++++++++++++++++++"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "+++++++++++++++++++++++++++++++++++++"


# ---------------------------------------------------------------------------
# Generic things
# ---------------------------------------------------------------------------

# Libtool
LT_INIT

# Checks for other programs.

AC_PROG_MAKE_SET
AC_PROG_LN_S
# we want flex, and not only lex
AC_CHECK_PROG([FLEX], [flex],[yes],[no])
# but this seems to be needed as well to set some include file correctly
# (seen on SuSE SLES 9)
AC_PROG_LEX
# and we need flex++
AC_CHECK_PROG([FLEXPP],[flex++],[yes],[no])
# and xxd
AC_CHECK_PROG([HAVE_XXD], [xxd],[yes],[no])

# when using homebrew on Mac, take the xsltproc installed
# additionally, have it use the XML_CATALOG_FILES env variable as it does
#  not automagically find the catalog files like on Unix systems.
XSLTPROC_BIN=xsltproc
XML_CATALOG=
if test "x${with_brew}" = "xyes" ; then
  XSLTPROC_BIN=${BREWBASE}/opt/libxslt/bin/xsltproc
  XML_CATALOG=XML_CATALOG_FILES=${BREWBASE}/etc/xml/catalog
fi
AC_SUBST(XSLTPROC_BIN,${XSLTPROC_BIN})
AC_SUBST(XML_CATALOG,${XML_CATALOG})

AC_CHECK_PROG([mHAVE_DBLATEX], [dblatex],[1],[0])
AM_CONDITIONAL([HAVE_DBLATEX], [test $mHAVE_DBLATEX = 1])

# Checks for header files.
AC_HEADER_STDC
#AC_CHECK_HEADERS([limits.h stdlib.h string.h strings.h unistd.h expat.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([floor memset pow sqrt])

# Enable large file support
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

# -----------------------------------------------------------------------------
# check for isblank(3)
# copied from aclocal.m4 file for sudo
#
# Copyright (c) 1994-1996, 1998-2008
#         Todd C. Miller <Todd.Miller@courtesan.com>
#
#    Permission to use, copy, modify, and distribute this software for any
#    purpose with or without fee is hereby granted, provided that the above
#    copyright notice and this permission notice appear in all copies.

AC_DEFUN([SUDO_FUNC_ISBLANK],
  [AC_CACHE_CHECK([for isblank], sudo_cv_func_isblank,
    [AC_TRY_LINK([#include <ctype.h>], [(void)isblank(1);],
    sudo_cv_func_isblank=yes, sudo_cv_func_isblank=no)])
] [
  if test "$sudo_cv_func_isblank" = "yes"; then
    AC_DEFINE(HAVE_ISBLANK, 1, [Define if you have isblank(3).])
  fi
])

# end of copied block
# -----------------------------------------------------------------------------

SUDO_FUNC_ISBLANK

echo
echo "+++++++++++++++++++++++++++++++++++++"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "+++++++++++++++++++++++++++++++++++++"


# check BOOST before anything else, it's a prime candidate
#  for being installed by the system (/usr/...) but as well
#  in another (newer?) version somewhere else.
# so make sure the include flags from boost come first

# BOOST >= ??? and <= 1.60 have an error which leads to Mac binaries
#   compiled on new OSX versions (El Capitan) not being able to run
#   on old OSX versions
if test "${systemuname}" = "Darwin"; then
  AX_BOOST_BASE([1.61.0])
else
  AX_BOOST_BASE([1.48.0])
fi

echo
echo "+++++++++++++++++++++++++++++++++++++"
echo " BOOST_CPPFLAGS:    ${BOOST_CPPFLAGS}"
echo " BOOST LDFLAGS:     ${BOOST_LDFLAGS}"
echo " ac_boost_path:     ${ac_boost_path}"
echo " ac_boost_lib_path: ${ac_boost_lib_path}"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "+++++++++++++++++++++++++++++++++++++"

# if the user gave --with-boost but no --with-boost-libdir, set it automatically
# because if not, the checks which need to link further down might fail
if test "${ac_boost_path}" != "" ; then
  if test "${ac_boost_lib_path}" = "" ; then
    ac_boost_lib_path=${ac_boost_path}/lib
    echo "Setting rescue value for ac_boost_lib_path: ${ac_boost_lib_path}"
    echo "If the rescue value is wrong, use --with-boost-libdir=/path/to/dir/with/boostlibs"
    # replace or append? atm replace
    BOOST_LDFLAGS="-L${ac_boost_lib_path}"
  fi
fi

# I have no idea why I need this ... but looks like
#  the autoconf macros have problems in automatically finding
#  the best paths on some systems
#
# So ... if we want boost (we do) and something went wrong with
#  automatic configuration, try setting rescue values which might work

possibleboostproblem="no"

if test "${BOOST_CPPFLAGS}" = "" ; then

  possibleboostproblem="yes"

  echo "Automatic recognition of Boost CPP flags failed? I have"
  echo "////// BOOST_CPPFLAGS:   ${BOOST_CPPFLAGS}"
  echo "////// ac_boost_path:    ${ac_boost_path}"
  echo "Setting rescue values!"

  if test "${systemuname}" = "SunOS" ; then
    BOOST_CPPFLAGS="-pthreads"
  elif test "x$build_os" = "xming32" ; then
    BOOST_CPPFLAGS="-mthreads"
  else
    BOOST_CPPFLAGS="-pthread"
  fi

  if test "${ac_boost_path}" = "" ; then
    BOOST_CPPFLAGS="${BOOST_CPPFLAGS} -I $/usr/include"
    echo "If the rescue below is wrong, use --with-boost=/path/to/boost"
  else
    BOOST_CPPFLAGS="${BOOST_CPPFLAGS} -I ${ac_boost_path}/include"
  fi

  echo "Choosen: ${BOOST_CPPFLAGS}"
  #LDFLAGS="${LDFLAGS} ${BOOST_THREAD_LIB}"
fi

if test "${BOOST_LDFLAGS}" = "" ; then
  possibleboostproblem="yes"
  echo "Automatic recognition of Boost LD failed? I have"
  echo "////// BOOST_CPPFLAGS:    ${BOOST_CPPFLAGS}"
  echo "////// BOOST LDFLAGS:     ${BOOST_LDFLAGS}"
  echo "////// ac_boost_path:     ${ac_boost_path}"
  echo "////// ac_boost_lib_path: ${ac_boost_lib_path}"

  if test "${ac_boost_lib_path}" = "" ; then
    if test "${enable_64}" = "yes"; then
      BOOST_LDFLAGS="-L/usr/lib64"
    else
      BOOST_LDFLAGS="-L/usr/lib"
    fi
    echo "If the rescue value below is wrong, use --with-boost-libdir=/path/to/dir/with/boostlibs"
  else
    BOOST_LDFLAGS="-L${ac_boost_lib_path}"
  fi

  echo "Setting rescue value: ${BOOST_LDFLAGS}"
fi


# some of the BOOST checks below break on older systems with newer BOOSTs
# e.g., boost_thread implicitly needs boost_clock and boost_system,
#  but the linker on KUbuntu 9.10 somehow does not find them automatically
# Way out: for the boost checks, set the linker rpath to include boost lib path
#  (only valid for gcc)
oldCFLAGS=${CFLAGS}
if test $bachCCOMPILER = gcc; then
  CFLAGS="${CFLAGS} -Wl,-rpath,${ac_boost_lib_path}"
fi

echo
echo "+++++++++++++++++++++++++++++++++++++"
echo " BOOST_CPPFLAGS:    ${BOOST_CPPFLAGS}"
echo " BOOST LDFLAGS:     ${BOOST_LDFLAGS}"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "+++++++++++++++++++++++++++++++++++++"

AX_BOOST_THREAD

echo
echo "+++++++++++++++++++++++++++++++++++++"
echo " BOOST_CPPFLAGS:    ${BOOST_CPPFLAGS}"
echo " BOOST LDFLAGS:     ${BOOST_LDFLAGS}"
echo " ac_boost_path:     ${ac_boost_path}"
echo " ac_boost_lib_path: ${ac_boost_lib_path}"
echo " BOOST_THREAD_LIB:  ${BOOST_THREAD_LIB}"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "+++++++++++++++++++++++++++++++++++++"

if test "${BOOST_THREAD_LIB}" = "" ; then
  echo "Automatic recognition of Boost thread lib failed? I have"
  echo "////// BOOST_CPPFLAGS:    ${BOOST_CPPFLAGS}"
  echo "////// BOOST LDFLAGS:     ${BOOST_LDFLAGS}"
  echo "////// ac_boost_path:     ${ac_boost_path}"
  echo "////// ac_boost_lib_path: ${ac_boost_lib_path}"
  echo "////// BOOST_THREAD_LIB:  ${BOOST_THREAD_LIB}"
  echo "////// BOOST_ROOT:        ${BOOST_ROOT}"

  case "${BOOST_LDFLAGS}" in
    *stage* )
    	 isstaged="yes"
    break;
    ;;
  esac

  if test "${isstaged}" = "yes" ; then
      echo
      echo
      echo "You seem to have a staged BOOST."
      echo
      if test "${BOOST_ROOT}" = "" ; then
	  echo "Please set \$BOOST_ROOT in your environment (via export or setenv, depending"
	  echo "on your shell) to the directory with the staged boost."
	  echo "E.g.: export BOOST_ROOT=\"/path/to/boost_1_39_0\""
      else
	  echo "You set \$BOOST_ROOT to: ${BOOST_ROOT}"
	  echo "Maybe this is not correct?"
	  echo
	  echo "If it is correct, please do not use a staged BOOST but a fully"
  	  echo "installed one. Apparently, using a staged BOOST does not work"
	  echo "on your system (sorry)."
      fi
  else
      echo "Please give the configure script the location of your BOOST"
      echo " libraries. Or, if needed, the name of your thread library (without"
      echo " 'lib' in front):"
      echo "   --with-boost[=DIR]"
      echo "or --with-boost-libdir=LIB_DIR"
      echo "or --with-boost-thread[=special-lib]"
  fi
  echo
  AC_MSG_ERROR([[Stopping here as BOOST thread lib must be available at this stage. Sorry.]])
fi

CXXFLAGS="${CXXFLAGS} ${BOOST_CPPFLAGS}"
LDFLAGS="${LDFLAGS} ${BOOST_LDFLAGS}"
LIBS="${LIBS} ${BOOST_THREAD_LIB}"

echo "+++++++++++++++++++++++++++++++++++++"
echo " BOOST_CPPFLAGS:    ${BOOST_CPPFLAGS}"
echo " BOOST LDFLAGS:     ${BOOST_LDFLAGS}"
echo " ac_boost_path:     ${ac_boost_path}"
echo " ac_boost_lib_path: ${ac_boost_lib_path}"
echo " BOOST_THREAD_LIB:  ${BOOST_THREAD_LIB}"
echo
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "LIBS now: ${LIBS}"
echo "+++++++++++++++++++++++++++++++++++++"

# BOOST_FILESYSTEM will need BOOST_SYSTEM on Macs
AX_BOOST_SYSTEM

echo "+++++++++++++++++++++++++++++++++++++"
echo " BOOST_CPPFLAGS:    ${BOOST_CPPFLAGS}"
echo " BOOST LDFLAGS:     ${BOOST_LDFLAGS}"
echo " ac_boost_path:     ${ac_boost_path}"
echo " ac_boost_lib_path: ${ac_boost_lib_path}"
echo " BOOST_SYSTEM_LIB:   ${BOOST_SYSTEM_LIB}"
echo

LIBS="${LIBS} ${BOOST_SYSTEM_LIB}"

echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LIBS now: ${LIBS}"
echo "+++++++++++++++++++++++++++++++++++++"

AX_BOOST_FILESYSTEM

echo "+++++++++++++++++++++++++++++++++++++"
echo " BOOST_CPPFLAGS:    ${BOOST_CPPFLAGS}"
echo " BOOST LDFLAGS:     ${BOOST_LDFLAGS}"
echo " ac_boost_path:     ${ac_boost_path}"
echo " ac_boost_lib_path: ${ac_boost_lib_path}"
echo " BOOST_FILESYSTEM_LIB:   ${BOOST_FILESYSTEM_LIB}"
echo

LIBS="${LIBS} ${BOOST_FILESYSTEM_LIB}"

echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LIBS now: ${LIBS}"
echo "+++++++++++++++++++++++++++++++++++++"


AX_BOOST_IOSTREAMS

echo "+++++++++++++++++++++++++++++++++++++"
echo " BOOST_CPPFLAGS:    ${BOOST_CPPFLAGS}"
echo " BOOST LDFLAGS:     ${BOOST_LDFLAGS}"
echo " ac_boost_path:     ${ac_boost_path}"
echo " ac_boost_lib_path: ${ac_boost_lib_path}"
echo " BOOST_IOSTREAMS_LIB:   ${BOOST_IOSTREAMS_LIB}"
echo

# The "-lz" here is no error: else it does not link on Trinity, no matter where I put -lz
#  in the Makefile-am of src/progs :-(

LIBS="${LIBS} ${BOOST_IOSTREAMS_LIB}"

echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LIBS now: ${LIBS}"
echo "+++++++++++++++++++++++++++++++++++++"



# check whether rt library is needed. BOOST thread needs
#  clock_gettime, on some systems it's in rt lib which is
#  NOT automatically used
# problem: AC_CHECK_LIB prepends -lrt to LIBS, we need
# to have it *after* the boost libs
oldLIBS=$LIBS
AC_CHECK_LIB(rt, clock_gettime)
if test "${LIBS}" != "${oldLIBS}"; then
  LIBS="$oldLIBS -lrt"
  echo "+++++++++++++++++++++++++++++++++++++"
  echo "Need to specifically use rt library"
  echo "LIBS now: ${LIBS}"
  echo "+++++++++++++++++++++++++++++++++++++"
fi


CFLAGS=${oldCFLAGS}

# check for expat and zlib

echo "BOOST is through, yay."
echo "+++++++++++++++++++++++++++++++++++++"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "LIBS now: ${LIBS}"
echo "+++++++++++++++++++++++++++++++++++++"

# BaCh 2013-03-08: AX_LIB_EXPAT serial 9 overwrites $LIBS! Make sure to use at least serial 10

AX_LIB_EXPAT
if test "${HAVE_EXPAT}" = "yes"; then
  CFLAGS="$CFLAGS ${EXPAT_CFLAGS}"
  CXXFLAGS="$CXXFLAGS ${EXPAT_CFLAGS}"
  LDFLAGS="${LDFLAGS} ${EXPAT_LDFLAGS}"
  LIBS="${LIBS} ${EXPAT_LIBS}"
fi


echo
echo "+++++++++++++++++++++++++++++++++++++"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "LIBS now: ${LIBS}"
echo "+++++++++++++++++++++++++++++++++++++"

# check for zlib
AX_CHECK_ZLIB()
LIBS="${LIBS} -lm -lz"

echo
echo "+++++++++++++++++++++++++++++++++++++"
echo "CFLAGS now: ${CFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "LIBS now: ${LIBS}"
echo "+++++++++++++++++++++++++++++++++++++"

AC_LANG_SAVE
AC_LANG_C
AC_CHECK_LIB([z], [gzoffset], [zlib_havegzoffset=yes], [zlib_havegzoffset=no])
AC_LANG_RESTORE
if test "$zlib_havegzoffset" = "yes"; then
  AC_DEFINE_UNQUOTED(HAVE_GZOFFSET)
fi

have_bz2="no"
AC_CHECK_LIB([bz2], [BZ2_bzRead], [have_bz2="yes"], [])

if test "${have_bz2}" = "yes"; then
  LIBS="${LIBS} -lbz2"
fi

#ZLIB_VERSION=`$EGREP "define ZLIB_VERSION" $ZLIB_HOME/include/zlib.h | $SED -e 's/[[^0-9\.]]//g'`
#AC_MSG_RESULT([$ZLIB_VERSION])

# ---------------------------------------------------------------------------
# check for lex/flex version
# ---------------------------------------------------------------------------
flexgoodversion="unknown"
if test "${LEX}" = "flex"; then
  AC_MSG_CHECKING(for flex version compatible with MIRA)
  # following line shamelessly stolen from the PHP build package
  flexversion=`${LEX} -V -v --version | ${SED} -e 's/^.* //'`

  flex_major=`echo ${flexversion} | cut -d '.' -f 1`
  flex_minor=`echo ${flexversion} | cut -d '.' -f 2`
  flex_revision=`echo ${flexversion} | cut -d '.' -f 3`

  # everything <= 2.5.35 is bad (compile error on clang 7)
  # tested with 2.6.4, that works, so let's take that.
  if test ${flex_major} -lt 2; then
    flexgoodversion="no"
  else
    if test ${flex_minor} -lt 6; then
      # this is the <2.6 branch
      flexgoodversion="no"
    else
      if test ${flex_revision} -lt 4; then
        flexgoodversion="unknown"
      else
        flexgoodversion="yes"
      fi
    fi
  fi

  AC_MSG_RESULT(${flexgoodversion}, ${flexversion} does not work for MIRA)
else
  flexgoodversion='no'
fi


# ---------------------------------------------------------------------------
# Set flag whether building environment is known to work
# Additionally, if special things are needed, set them here.
# ---------------------------------------------------------------------------
buildenvworks="unknown"
if test "${systemuname}" = "Darwin"; then
  buildenvworks="yes"
fi
if test "${systemuname}" = "Linux"; then
  buildenvworks="yes"
fi
if test "${systemuname}" = "SunOS"; then
  buildenvworks="yes"
fi
if test "${systemuname}" = "FreeBSD"; then
  buildenvworks="yes"
fi




# Test if --enable-debug given
AC_ARG_ENABLE(debug,
	      AC_HELP_STRING([--enable-debug=no/yes],
              [use 'yes' to enable debugging information in the executable)]))
if test "${enable_debug}" = "yes" ; then
  CFLAGS="$CFLAGS -g"
  CXXFLAGS="$CXXFLAGS -g"
fi

# Test if --enable-cxxdebug given
AC_ARG_ENABLE(cxxdebug,
	      AC_HELP_STRING([--enable-cxxdebug=no/yes],
              [use 'yes' to enable debugging version of STL objects)]))
if test "${enable_cxxdebug}" = "yes" ; then
  CXXFLAGS="$CXXFLAGS -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC"
  # CXXFLAGS="$CXXFLAGS -D_GLIBCXX_DEBUG"
fi



# Test if --enable-profile given
AC_ARG_ENABLE(profile,
	      AC_HELP_STRING([--enable-profile=no/yes],
              [use 'yes' to have profile information in the executable)]))

AC_MSG_CHECKING(for profiling settings)
if test "${enable_profile}" = "yes" ; then
  CFLAGS="$CFLAGS -pg"
  CXXFLAGS="$CXXFLAGS -pg"
  LDFLAGS="${LDFLAGS} -pg"
else
  enable_profile="no"
fi
AC_MSG_RESULT(${enable_profile})

## -flto (link-time-optimiser) for gcc
## not always supported yet, does not work well will -pg
## therefore enable it only if -pg not used
#if test $bachCCOMPILER = gcc; then
#  if test "${enable_profile}" != "yes" ; then
#    AC_MSG_CHECKING(for using gcc -flto)
#    canuseflto="no";
#    oldCFLAGS="$CFLAGS"
#    CFLAGS="$CFLAGS -flto"
#    AC_TRY_COMPILE(, [return 0;], [canuseflto="yes"],)
#    CFLAGS="$oldCFLAGS"
#    AC_MSG_RESULT(${canuseflto})
#    if test "${canuseflto}" = "yes" ; then
#      CFLAGS="$CFLAGS -flto"
#      CXXFLAGS="$CXXFLAGS -flto"
#    fi
#  fi
#fi



# Turn off irritating linker warnings in IRIX
if test "${systemuname}" = "IRIX" || test "${systemuname}" = "IRIX64"; then
  CFLAGS="-Wl,-LD_MSG:off=85:off=84:off=16:off=134 $CFLAGS"
  CXXFLAGS="-Wl,-LD_MSG:off=85:off=84:off=16:off=134 $CXXFLAGS"
fi


# Test if --enable-warnings given
AC_ARG_ENABLE(warnings,
	      AC_HELP_STRING([--enable-warnings=no/yes],
              [use 'yes' to enable warnings (-W and other options on gcc compiler)]))
if test "${enable_warnings}" = "yes" ; then
  if test $bachCCOMPILER = gcc; then

# that one pukes way too often with boost 1.50
# -Wold-style-cast

    CXXFLAGS="${CXXFLAGS} -W -Wstrict-prototypes -fdiagnostics-show-option -Wextra -Wshadow -fno-strict-aliasing -ansi -Wcast-align -Wcast-qual -Wchar-subscripts -Wpointer-arith -Wredundant-decls -Wwrite-strings -Woverloaded-virtual -Wsign-promo -Wformat -Wimplicit -Wmissing-braces -Wparentheses -Wswitch -Wswitch-default -Wswitch-enum -Wunused -Wfloat-equal -Wendif-labels -Wconversion -Wdisabled-optimization -Wsequence-point"
    CFLAGS="${CFLAGS} -W -Wstrict-prototypes -fno-strict-aliasing -Wextra -ansi -Wbad-function-cast -Wcast-align -Wcast-qual -Wchar-subscripts -Wmissing-prototypes -Wnested-externs -Wpointer-arith -Wredundant-decls -Wshadow -Wstrict-prototypes -Wwrite-strings -Wformat -Wimplicit -Wmissing-braces -Wparentheses -Wswitch -Wswitch-default -Wswitch-enum -Wunused -Wfloat-equal -Wendif-labels -Wconversion -Wdisabled-optimization -Wsequence-point"

    # unused values may be a problem (e.g. using == in an assignment),
    #  but ...
    # this one pukes on OSX High Sierra with:
    #   /usr/include/signal.h:125:26: error: right operand of comma operator has no effect [-Werror=unused-value]
    #   #define sigemptyset(set) (*(set) = 0, 0)
    if test "${systemuname}" != "Darwin"; then
      CXXFLAGS="${CXXFLAGS} -Werror=unused-value"
      CFLAGS="${CXXFLAGS} -Werror=unused-value"
    fi

  fi
fi

if test $bachCCOMPILER = gcc; then
  CXXFLAGS="${CXXFLAGS} -Werror=uninitialized -Werror=return-type -Werror=parentheses"
  CFLAGS="${CFLAGS} -Werror=uninitialized -Werror=return-type "
fi

# hack for CygWin/gcc where mktmpnam and a couple of others are not impolemented in strict ansi
# and -std=c++11 implies strict ansi
if test $bachCCOMPILER = gcc; then
  if test $build_os = cygwin ; then
    CFLAGS="${CFLAGS} -U__STRICT_ANSI__"
    CXXFLAGS="${CXXFLAGS} -U__STRICT_ANSI__"
  fi
fi



# ---------------------------------------------------------------------------
# MIRA specific flags and settings
# ---------------------------------------------------------------------------

# Check if we have to compile mira in quiet mode
# Non-quiet mode has some extra output on stdout during assembly

AC_MSG_CHECKING(whether mira should compile in public quiet mode)
def_publicquietmira="yes"
AC_ARG_ENABLE([publicquietmira],
	AC_HELP_STRING([--enable-publicquietmira=yes/no],
	[use 'no' for a chatty version of MIRA, default is yes]),
	[publicquietmira=$enableval],
	[publicquietmira=$def_publicquietmira])
if test "$publicquietmira" = "no"; then
   CXXFLAGS="-DBE_AWARE_THAT_THIS_VERSION_OF_MIRA_MIGHT_BE_A_BIT_CHATTY ${CXXFLAGS}"
else
   CXXFLAGS="-DPUBLICQUIET ${CXXFLAGS}"
fi
AC_MSG_RESULT($publicquietmira)



# compiles in some basic checks into mira that look for sanity within
#  some functions
# leaving this on "yes" is encouraged, impact on run time is minimal
AC_MSG_CHECKING(whether mira should compile with bound tracking)
def_boundtracking="yes"
AC_ARG_ENABLE([boundtracking],
	AC_HELP_STRING([--enable-boundtracking=yes/no],
	[use 'yes' for some additional value bound tracking in MIRA, default is yes]),
	[boundtracking=$enableval],
	[boundtracking=$def_boundtracking])
if test "$boundtracking" = "yes"; then
   AC_DEFINE_UNQUOTED(BOUNDTRACKFLAG)
fi
AC_MSG_RESULT($boundtracking)


# compiles in some basic checks into mira that look for sanity within
#  some functions
# as with boundchecking, leaving this on "yes" is encouraged,
#  impact on run time is minimal
AC_MSG_CHECKING(whether mira should compile with bug tracking)
def_bugtracking="yes"
AC_ARG_ENABLE([bugtracking],
	AC_HELP_STRING([--enable-bugtracking=yes/no],
	[use 'yes' for some additional bug tracking in MIRA, default is yes]),
	[bugtracking=$enableval],
	[bugtracking=$def_bugtracking])
if test "$bugtracking" = "yes"; then
   AC_DEFINE_UNQUOTED(BUGTRACKFLAG)
fi
AC_MSG_RESULT($bugtracking)


# compiles in extra checks into mira that look for sanity within
#  some functions
# contrary to the other checks above, this extended bugtracking *is*
#  costly and impacts runtime noticeably, use only for bughunting on
#  core dumps or similar
AC_MSG_CHECKING(whether mira should compile with extended bug tracking)
def_extendedbugtracking="no"
AC_ARG_ENABLE([extendedbugtracking],
	AC_HELP_STRING([--enable-extendedbugtracking=no/yes],
	[use 'yes' for some additional bug tracking in MIRA, default is no]),
	[extendedbugtracking=$enableval],
	[extendedbugtracking=$def_extendedbugtracking])
if test "$extendedbugtracking" = "yes"; then
   AC_DEFINE_UNQUOTED(PARANOIABUGTRACKFLAG)
fi
AC_MSG_RESULT($extendedbugtracking)


# traces every function call and exit on stdout
# setting function trace to "yes" is *highly* discouraged except for
#  toy projects
AC_MSG_CHECKING(whether mira should compile with function trace)
def_functrace="no"
AC_ARG_ENABLE([functrace],
	AC_HELP_STRING([--enable-functrace=no/yes],
	[use 'yes' for a function trace in MIRA, default is no]),
	[functrace=$enableval],
	[functrace=$def_functrace])
if test "$functrace" = "yes"; then
   AC_DEFINE_UNQUOTED(TRACEFLAG)
fi
AC_MSG_RESULT($functrace)


AC_ARG_ENABLE(miramemorc,
	AC_HELP_STRING([--enable-miramemorc],
	[enables building with memory overrun checks, default is no]))
AM_CONDITIONAL([BUILDMIRAMEMORC], [test x$enable_miramemorc = xyes])

if test "$enable_miramemorc" = "yes"; then
   CXXFLAGS="-DMIRAMEMORC ${CXXFLAGS}"
fi


# -Winline -Wmissing-noreturn???
# -Weffc++ (C++ only, from Scott Myers' book)

#CXXFLAGS="${CXXFLAGS} -march=native"
#CFLAGS="${CFLAGS} -march=native"

#CXXFLAGS="${CXXFLAGS} -g -O0 -Wall"
#CFLAGS="${CFLAGS} -g -O0 -Wall"

# -Werror ???
# CXXFLAGS="-O3 -funroll-loops -march=i686 -fprofile-arcs"
# CFLAGS="-O3 -funroll-loops -march=i686 -fprofile-arcs"
# CXXFLAGS="-O3 -funroll-loops -march=i686 -fbranch-probabilities"
# CFLAGS="-O3 -funroll-loops -march=i686 -fbranch-probabilities"
# from KDE3.3 konstruct
# -falign-functions=4 -fomit-frame-pointer -mfancy-math-387 -mcpu=athlon-xp



# Various settings
# AC_COMPILE_WARNINGS

#AC_CONFIG_FILES([Makefile])

# OK, get those OpenMP flags into the compile flags
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"


AC_OUTPUT(Makefile \
		    src/Makefile \
		    src/other/Makefile \
		    src/other/sls/Makefile \
		    src/other/sls/mira-createsls/Makefile \
		    src/other/sls/mira-createsls/bin/Makefile \
		    src/stdinc/Makefile \
		    src/debuggersupport/Makefile \
		    src/memorc/Makefile \
		    src/errorhandling/Makefile \
		    src/util/Makefile \
		    src/io/Makefile \
		    src/mira/Makefile \
		    src/caf/Makefile \
		    src/support/Makefile \
		    src/scripts/Makefile \
		    src/modules/Makefile \
		    src/version/Makefile \
		    src/progs/Makefile \
		    man/Makefile \
		    doc/Makefile \
		    doc/docbook/Makefile \
		    distribution/Makefile
         )

# write docbook version

cat <<EOF >doc/docbook/versionfile
<?xml version="1.0" ?>
<!DOCTYPE book PUBLIC "-//OASIS//DTD DocBook XML V4.5//EN" "http://www.docbook.org/xml/4.5/docbookx.dtd">
 <releaseinfo>MIRA Version ${VERSION}</releaseinfo>
EOF


# give summary of the most important settings
MIRAVERSION=`build-aux/git-version-gen .tarball-version`_${build_os}_${BITPACKAGE}_${MIRASTATIC}

echo ${VERSION} >.version

echo
echo "--------------------------------------------------------------------------------"
echo "                       Summary of MIRA configuration"
echo "--------------------------------------------------------------------------------"
echo

possibleproblems="no"
if test "${HAVE_XXD}" = "no"; then
  echo
  echo "Could not find the program 'xxd'. It is usually part of 'vim', please"
  echo "install it."
  possibleproblems="yes"
fi
if test "${HAVE_EXPAT}" = "no"; then
  echo
  echo "Could not find expat. Please install the Expat library."
  echo "Have a look at your distribution which should already have"
  echo " an already prepackaged version for you to install."
  echo "If not, have a look at http://sourceforge.net/projects/expat/"
  echo
  echo "If you have expat installed at a non-standard location, please tell"
  echo " this to the configure script with the '--with-expat*' options."
  possibleproblems="yes"
fi
if test "${LEX}" = "no"; then
  echo
  echo "lex/flex does not seem to exist. Please install GNU flex (at least 2.6.4)."
  possibleproblems="yes"
fi
if test "${flexgoodversion}" = "no"; then
  if test "${LEX}" = "flex"; then
    echo
    echo "Your 'flex' is version ${flexversion} which is known NOT to work with MIRA"
    echo "Please install newer GNU flex (2.6.4 is known to work)."
    echo "See http://flex.sourceforge.net/"
  else
    echo
    echo "Your lex is not flex, but is seen as: ${LEX}"
    echo "Please install GNU flex (at least 2.6.4)"
  fi
  possibleproblems="yes"
fi

if test "${FLEXPP}" = "no"; then
  echo
  echo "You have flex installed, but no executable 'flex++'."
  echo "This is easy to fix: create a link from 'flex++' to flex. This can be made"
  echo "in the directory where flex resides or, if you do not have sufficient rights"
  echo "for that exact directory, in any directory which is in your PATH variable."
  possibleproblems="yes"
fi

if test "${HAVE_BOOST_THREAD}" = "yes"; then
  echo
  echo "I couldn't find the C++ Boost libraries. They're really needed."
  echo "Install at least 1.44 from the package repository of your"
  echo " distribution/vendor or download and install manually from"
  echo " http://www.boost.org/"
  possibleproblems="yes"
fi

if test "${possibleproblems}" = "yes"; then
  echo
  echo "Packages/programs needed for the compilation were not found."
  echo "Please install these first and then try again."
  echo
  echo "If you installed them, perhaps the programs are not in your PATH?"
  echo "Or you forgot to give the configure script options to tell where"
  echo "the packages are?"
  echo
  exit 10
fi

if test "${possibleboostproblem}" = "" ; then
  echo
  echo "Hmmm ... looks like there might be a problem with"
  echo "the C++ Boost libraries. Please use one of the Boost"
  echo "configuration options, giving them the correct location."
  echo "   --with-boost[=DIR]"
  echo "or --with-boost-libdir=LIB_DIR"
  echo "or --with-boost-thread[=special-lib]"
  possibleproblems="yes"
fi

if test "${LEX}" = "lex"; then
  echo "You have 'lex' on your system, not 'flex'. This may be a problem."
  echo
  possibleproblems="yes"
fi

if test "${LEX}" = "flex"; then
  if test "${flexgoodversion}" = "unknown"; then
    echo "Your 'flex' is version ${flexversion}, the earliest version known to work is 2.6.4"
    echo "This might be a problem, but does not need to be."
    echo
    possibleproblems="yes"
  fi

  tst_flex_version=`flex -V -v --version | ${SED} -e 's/^.* //'`
  tst_flexpp_version=`flex++ -V -v --version | ${SED} -e 's/^.* //'`
  if test "$tst_flex_version" != "$tst_flexpp_version"; then
    echo "flex installed has version  : $tst_flex_version"
    echo "flex++ installed has version: $tst_flexpp_version"
    echo "These versions are not the same and you might run into trouble."
    echo
    possibleproblems="yes"
  fi
fi

if test "${zlib_havegzoffset}" = "no"; then
  possibleproblems="yes"
  echo "While zlib has been found, the version which is used is apparently"
  echo "lower than 1.2.4. Progress bars for FASTQ loading will not work."
  echo "Please make sure you have a good version of zlib."
  echo
fi

if test "${have_openmp}" = "no"; then
  possibleproblems="yes"
  echo "Did not find a way to compile with OpenMP. This is not critical,"
  echo "but some parts of mirabait will run slower."
  echo "Please make sure you have OpenMP installed."
  echo
fi

echo

echo "Building version ................................. ${MIRAVERSION}"
echo

echo "We are building on ............................... ${systemuname}"

echo "CPU supports 64 bit? ............................. ${has_int64_t}"
echo "Compiler family .................................. ${bachCCOMPILER}"

if test "${systemuname}" = "Darwin"; then
  if test $bachCCOMPILER != clang; then
    possibleproblems="yes"
    echo
    echo "System is Apple, but compiler is not clang. I was able to compile,"
    echo "but the program crashed upon invocation."
    echo
  fi
fi

$ECHO_N "Compiler optimisation flags used ................. "
if test "${enable_optimisations}" = "yes"; then
  echo ${OPTFLAGS}
else
  echo "flags disabled on request"
fi

if test "${enable_debug}" = "yes"; then
  $ECHO_N "Debug information in executable?.................. "
  echo "${enable_debug}"
fi

if test "${enable_miramemorc}" = "yes"; then
  $ECHO_N "Memory overrun checks in executable? ............. "
  echo "${enable_debug}"
fi

if test "${enable_cxxdebug}" = "yes"; then
  $ECHO_N "Debug version of C++ STL library? ................ "
  echo "${enable_debug}"
fi

if test "${enable_profile}" = "yes"; then
  $ECHO_N "Profile information in executable? ............... "
  echo "${enable_profile}"
fi
if test "${has_int64_t}" = "yes"; then
  echo "Do we compile in 64 bit? ......................... ${enable_64}"
  if test "${enable_64}" = "no"; then
    echo "Known and tested 32 bit compilation support? ..... ${tested32bitsupport}"
    if test "${tested32bitsupport}" = "no"; then
      possibleproblems="yes"
      echo "THE LINE ABOVE MEANS WE MIGHT HAVE A PROBLEM FOR COMPILATION!"
    fi
  fi
fi

echo "Building static? ................................. ${enable_mirastatic}"


$ECHO_N "BOOST expected to compile and link? .............. "
if test "${possibleboostproblem}" = "yes"; then
  echo "maybe"
  possibleproblems="yes"
else
  echo "yes"
fi
echo "Have OpenMP? ..................................... ${have_openmp}"
echo "Have complete zlib support? ...................... ${zlib_havegzoffset}"

echo
echo "Can MIRA be built on this system? ................ ${buildenvworks}"
if test "${buildenvworks}" != "yes"; then
  possibleproblems="yes"
  if test "${buildenvworks}" = "unknown"; then
    echo
    echo "No information regarding builds in this environment available. Please report"
    echo "success or failure."
    echo
  else
    echo
    echo "This system has been reported as having problems."
    echo "Please help fixing if you can."
    echo
  fi
fi

echo
echo "Using CFLAGS   : ${CFLAGS}"
echo "Using CXXFLAGS : ${CXXFLAGS}"
echo "Using LDFLAGS  : ${LDFLAGS}"
echo "Using LIBS     : ${LIBS}"
echo

if test "${possibleproblems}" = "no"; then
  echo "Looks good, we do not expect compilation problems."
else
  echo "System: `uname -a`"
  echo "Compiler: '${CC} -v' output:"
  ${CC} -v
  echo
  echo

  if test "${possibleboostproblem}" = "yes"; then
    echo "There was a problem in configuring Boost completely automatically."
    echo "The configure script tried to guess some directories, but in case"
    echo " there is a problem in compiling or linking with BOOST, you"
    echo " might want to test a few other ./configure options for Boost."
  fi
fi

echo "--------------------------------------------------------------------------------"

if test "${possibleproblems}" = "no"; then
  echo
  echo "Type 'make' to make the binaries"
  echo "and then 'make install' to install binaries and man pages."
  echo "Type 'make docinstall' to install extensive documentation as HTML and PDF."
else
  echo "Ooops, some functionality might be missing or we might encounter some"
  echo "problems in compilation/linking. Please give it a try and report success"
  echo "or failure."
  echo
  echo "Type 'make' to make the binaries and if no error occurred,"
  echo "then 'make install' to install them."
  echo
  echo "In case of failure, please help to fix if you can."
fi
