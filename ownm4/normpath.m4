dnl @synopsis adl_NORMALIZE_PATH(VARNAME, [REFERENCE_STRING])
dnl
dnl Perform some cleanups on the value of $VARNAME (interpreted as a
dnl path):
dnl
dnl   - empty paths are changed to '.'
dnl   - trailing slashes are removed
dnl   - repeated slashes are squeezed except a leading doubled slash '//'
dnl     (which might indicate a networked disk on some OS).
dnl
dnl REFERENCE_STRING is used to turn '/' into '\' and vice-versa: if
dnl REFERENCE_STRING contains some backslashes, all slashes and
dnl backslashes are turned into backslashes, otherwise they are all
dnl turned into slashes.
dnl
dnl This makes processing of DOS filenames quite easier, because you
dnl can turn a filename to the Unix notation, make your processing, and
dnl turn it back to original notation.
dnl
dnl   filename='A:\FOO\\BAR\'
dnl   old_filename="$filename"
dnl   # Switch to the unix notation
dnl   adl_NORMALIZE_PATH([filename], ["/"])
dnl   # now we have $filename = 'A:/FOO/BAR' and we can process it as if
dnl   # it was a Unix path.  For instance let's say that you want
dnl   # to append '/subpath':
dnl   filename="$filename/subpath"
dnl   # finally switch back to the original notation
dnl   adl_NORMALIZE_PATH([filename], ["$old_filename"])
dnl   # now $filename equals to 'A:\FOO\BAR\subpath'
dnl
dnl One good reason to make all path processing with the unix
dnl convention is that backslashes have a special meaning in many
dnl cases. For instance
dnl
dnl   expr 'A:\FOO' : 'A:\Foo'
dnl
dnl will return 0 because the second argument is a regex in which
dnl backslashes have to be backslashed. In other words, to have the two
dnl strings to match you should write this instead:
dnl
dnl   expr 'A:\Foo' : 'A:\\Foo'
dnl
dnl Such behavior makes DOS filenames extremely unpleasant to work
dnl with. So temporary turn your paths to the Unix notation, and revert
dnl them to the original notation after the processing. See the macro
dnl adl_COMPUTE_RELATIVE_PATHS for a concrete example of this.
dnl
dnl REFERENCE_STRING defaults to $VARIABLE, this means that slashes
dnl will be converted to backslashes if $VARIABLE already contains some
dnl backslashes (see $thirddir below).
dnl
dnl   firstdir='/usr/local//share'
dnl   seconddir='C:\Program Files\\'
dnl   thirddir='C:\home/usr/'
dnl   adl_NORMALIZE_PATH([firstdir])
dnl   adl_NORMALIZE_PATH([seconddir])
dnl   adl_NORMALIZE_PATH([thirddir])
dnl   # $firstdir = '/usr/local/share'
dnl   # $seconddir = 'C:\Program Files'
dnl   # $thirddir = 'C:\home\usr'
dnl
dnl @category Misc
dnl @author Alexandre Duret-Lutz <duret_g@epita.fr>
dnl @version 2001-05-25
dnl @license GPLWithACException

AC_DEFUN([adl_NORMALIZE_PATH],
[case ":[$]$1:" in
# change empty paths to '.'
  ::) $1='.' ;;
# strip trailing slashes
  :*[[\\/]]:) $1=`echo "[$]$1" | sed 's,[[\\/]]*[$],,'` ;;
  :*:) ;;
esac
# squeze repeated slashes
case ifelse($2,,"[$]$1",$2) in
# if the path contains any backslashes, turn slashes into backslashes
 *\\*) $1=`echo "[$]$1" | sed 's,\(.\)[[\\/]][[\\/]]*,\1\\\\,g'` ;;
# if the path contains slashes, also turn backslashes into slashes
 *) $1=`echo "[$]$1" | sed 's,\(.\)[[\\/]][[\\/]]*,\1/,g'` ;;
esac])
