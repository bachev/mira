GITSHORT=`$(top_srcdir)/build-aux/git-version-gen isshort $(top_srcdir)/.tarball-version`
MIRAVERSION=${GITSHORT}_${build_os}_${BITPACKAGE}_${MIRASTATIC}
MIRADIR=${PACKAGE_TARNAME}_${MIRAVERSION}

BUNDLETARGET=darwin

DYPATH=DYLD_LIBRARY_PATH=$(EXTRADYLDPATH)$(DYLD_LIBRARY_PATH)

3rdparty:
	rm -rf 3rdparty
	cp -R -L ../src/3rdparty .
	tar cvjf mira_3rdparty_`date +%d-%m-%Y`.tar.bz2 3rdparty

libbundle-nobundle:
	@echo "Do not need to bundle on this platform"

libbundle-darwin:
	@echo "Using Macdylibbundler to create OSX compatible dylib references"
	(cd ${MIRADIR}; $(DYPATH) dylibbundler -b -od -d lib -p @executable_path/../lib/ -x bin/mira; cd ..)

dbdata:
	mkdir -p ${MIRADIR}/dbdata
	cp ../src/other/sls/README.txt ${MIRADIR}/dbdata
	cp ../src/other/sls/mira-install-sls-rrna.sh ${MIRADIR}/dbdata
	cp ../src/other/sls/rfam_rrna-21-1.sls.gz ${MIRADIR}/dbdata

copybin:
	@echo "Making distribution for ${BITPACKAGE}"
	-mkdir -p ${MIRADIR}/bin
	cp ../src/progs/mira ${MIRADIR}/bin
	cp ../src/other/sls/mira-install-sls-rrna.sh ${MIRADIR}/bin
	ln -s -f mira ${MIRADIR}/bin/mirabait
	ln -s -f mira ${MIRADIR}/bin/miraconvert
	ln -s -f mira ${MIRADIR}/bin/miramem
	ln -s -f mira ${MIRADIR}/bin/miramer

copydocs:
	-mkdir -p ${MIRADIR}/docs
	-cp -R -L ../doc/docbook/Defi*.html ${MIRADIR}/docs
	-(cd ${MIRADIR}/docs; ln -s DefinitiveGuideToMIRA.html index.html)
	-cp -R -L ../doc/docbook/*.pdf ${MIRADIR}/docs
	cp -R -L ../doc/docbook/bookfigures ${MIRADIR}/docs
	cp -R -L ../doc/docbook/images ${MIRADIR}/docs
	cp -R -L ../doc/docbook/doccss ${MIRADIR}/docs
	-mkdir -p ${MIRADIR}/man/man1
	cp -R -L ../man/*.1 ${MIRADIR}/man/man1

distribin: clean copybin libbundle-${BUNDLETARGET}
	strip ${MIRADIR}/bin/mira[a-z]*
	tar cvjf ${PACKAGE_TARNAME}_binonly_${MIRAVERSION}.tar.bz2 ${MIRADIR}

distridbg: clean copybin libbundle-${BUNDLETARGET}
	tar cvjf ${PACKAGE_TARNAME}_dbgonly_${MIRAVERSION}.tar.bz2 ${MIRADIR}

distrib: clean dbdata copybin libbundle-${BUNDLETARGET} copydocs
	cp ../THANKS ${MIRADIR}
	cp README ${MIRADIR}
	cp LICENCE ${MIRADIR}
	cp -R -L ../src/support ${MIRADIR}
	rm -rf ${MIRADIR}/support/Makefile*
	cp ../src/mira/CHANGES.txt ${MIRADIR}
	cp ../src/mira/CHANGES_old.txt ${MIRADIR}

	tar cvjf ${PACKAGE_TARNAME}_${MIRAVERSION}_fullsymbols.tar.bz2 ${MIRADIR}
	strip ${MIRADIR}/bin/mira[a-z]*
	tar cvjf ${PACKAGE_TARNAME}_${MIRAVERSION}.tar.bz2 ${MIRADIR}


clean:
	-rm -rf ${MIRADIR} *.tar *.gz
	rm -rf mira_*
	rm -rf mira-*
	rm -rf 3rdparty*
