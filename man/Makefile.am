# mira, mirabait and miraconvert are the three "official" MIRA programs
#  useful to users. All others should not be documented atm, they're
#  helper programs in various stages of imperfectness.

dist_pkgdata_DATA = mira.1 mirabait.1 miraconvert.1
pkgdatadir=$(prefix)/share/man/man1

%.1:
	help2man --no-info $(top_srcdir)/src/progs/$* >$@

mclean:
	-rm -f $(dist_pkgdata_DATA)

distclean-local: mclean

manpages: mclean $(dist_pkgdata_DATA)
