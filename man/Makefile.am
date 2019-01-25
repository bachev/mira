# mira, mirabait and miraconvert are the three "official" MIRA programs
#  useful to users. All others should not be documented atm, they're
#  helper programs in various stages of imperfectness.

# manpages are created dynamically via help2man

dist_pkgdata_DATA = mira.1 mirabait.1 miraconvert.1
pkgdatadir=$(prefix)/share/man/man1

# the detour via the 'ttt' file and sed serves two purposes
# 1. if 'make manpages' is invoked but no binary was constructed, an empty
#    .1 file would be created by help2man which make would see as good enough
#    target in the next call and not re-create a valid manpage. Redirecting
#    to 'ttt' file negates that
# 2. When creating a source distribution archive via 'make dist', removes
#    the _0+g... string from the version number.
#    Rationale: have valid manpages included by default in the source package
#     when building from a git checkout directly from a tagged version
%.1:
	help2man --no-info $(top_srcdir)/src/progs/$* >ttt
	sed -e 's/\(MIRA[^ ]* V[0-9.]\{1,\}\)_0+g[0-9a-f]\{1,\}/\1/g' ttt >$@
	-rm ttt

mclean:
	-rm -f ttt $(dist_pkgdata_DATA)

distclean-local: mclean

manpages: mclean $(dist_pkgdata_DATA)
