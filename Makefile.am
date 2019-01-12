ACLOCAL_AMFLAGS = -I m4

BUILT_SOURCES = $(top_srcdir)/.version
$(top_srcdir)/.version:
	echo $(VERSION) > $@-t && mv $@-t $@


AM_CXXFLAGS = $(BOOST_CPPFLAGS)

SUBDIRS = src man
EXTRA_DIST = README_MAINTAINERS GETTING_STARTED THANKS build-aux


binaries: src
	make -C src all

manpages:
	make -C man manpages

docs: manpages doc doc/docbook
	make -C doc/docbook docs

docinstall:
	make -C doc/docbook install

distrib: binaries docs distribution
	make -C distribution distrib

distribin: binaries
	make -C distribution distribin

bclean:
	-dot_clean -n .
	-find . -name "*~" -exec rm -rf {} \;
	-find . -name ".pure" -exec rm -rf {} \;
	-find . -name "._*" -exec rm -rf {} \;
	-find . -name "*.o" -exec rm {} \;
	-find . -name "*.a" -exec rm {} \;
	-find . -name "bla_*" -exec rm {} \;
	-find . -name "bin.*" -exec rm {} \;
	-find . -name "log.*" -exec rm {} \;
	-find . -name "*.bak" -exec rm {} \;
	-find . -name "ttt*" -exec rm {} \;
	-find . -name "TTT*" -exec rm {} \;
	-find . -name "core" -exec rm {} \;
	-find . -name "*.flex.C" -exec rm {} \;
	-find . -name "*flexer.cc" -exec rm {} \;
	-find . -type d -name "*.dSYM" -exec rm -rf {} \;

distclean-local: bclean
	-find . -name ".deps" -exec rm -rf {} \;
	-rm -rf m4
	-rm -rf Makefile Makefile.in aclocal.m4 config.status config.h config.h.in config.log configure autom4te.cache config.guess config.sub missing libtool
	-rm -rf compile depcomp install-sh ltmain.sh
	-rm -rf config

clean-local: bclean
	make -C distribution clean
	make -C doc clean

cvsclean: bclean
	find . -name "CVS" -exec rm -rf {} \;

#       -find . -path './arch' -prune -o -perm -a+x -type f -exec rm {} \;
#	 -rm -rf .deps *.bak *.out *.o *.d *.flex.C core ttt* log.* tca.*
#	 -rm -rf main main2
#	 -rm -rf *.aux *.dvi *.log ii_files .ix* *_1000*





# tar to pick up the other directories
# then remove any GIT subdirectories and files
# DON'T do that in an one-liner (at least it breaks for me)

dist-hook:
	$(ECHO_N) `build-aux/git-version-gen isshort doesnotexist` > $(distdir)/.tarball-version
	cp $(distdir)/.tarball-version $(distdir)/.version
	( cd doc/docbook; make clean; cd ../.. )
	tar cBf - doc | ( cd $(distdir); tar xBf - )
	rm -rf $(distdir)/doc/docbook/db_guides
	rm -rf $(distdir)/doc/docbook/util
	rm -rf $(distdir)/doc/other
	( cd doc/docbook; make versionfile; cd ../.. )
	${XML_CATALOG} ${XSLTPROC_BIN} --xinclude --output $(distdir)/README_build.html $(top_srcdir)/doc/docbook/mira.xsl $(top_srcdir)/doc/docbook/chap_installation_part.xml
	cp -r $(top_srcdir)/doc/docbook/doccss $(distdir)
	tar cBf - minidemo | ( cd $(distdir); tar xBf - )
	( cd distribution; make clean; cd .. )
	tar cBf - distribution | ( cd $(distdir); tar xBf - )
	find $(distdir) -name ".git*" -exec rm -rf {} \;
