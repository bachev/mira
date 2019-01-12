
GITBUILD=`$(top_srcdir)/build-aux/git-version-gen isshort $(top_srcdir)/.tarball-version`

versionfile: $(top_srcdir)/.version
	@echo "Making versionfile"
	@cp version.stub versionfile
	@echo \<releaseinfo\>MIRA Version ${GITBUILD}\</releaseinfo\> >>versionfile


CHAPTERS = chap_specialparams_part.xml chap_denovo_part.xml chap_mapping_part.xml chap_seqtechdesc_part.xml chap_bitsandpieces_part chap_mirautils_part chap_est_part chap_faq_part chap_hard_part chap_logfiles_part chap_maf_part chap_reference_part  chap_installation_part chap_intro_part chap_results_part chap_sra_part chap_seqadvice_part chap_bonus_part chap_preface_part chap_3p_ghbambus_part.xml

CHAP_XML = $(CHAPTERS:%_part=%_part.xml)
CHAP_HTML = $(CHAPTERS:%_part=%_part.html)

EXTRA_DIST = ${CHAP_XML} book_definitiveguide.xml bookfigures book.xsl doccss images mira.xsl versionfile

if HAVE_DBLATEX
FINALDOCFILES = 3rdparty_MIRA.html 3rdparty_MIRA.pdf DefinitiveGuideToMIRA.html DefinitiveGuideToMIRA.pdf
else
FINALDOCFILES = 3rdparty_MIRA.html DefinitiveGuideToMIRA.html
endif


clean-local:
	rm -rf *html *pdf versionfile

%_part.html: %_part.xml
	@echo "Making .html out of _part.xml" $@
	${XML_CATALOG} ${XSLTPROC_BIN} --nonet --xinclude --output $@  $(top_srcdir)/doc/docbook/mira.xsl $<


chapters: ${CHAP_HTML} versionfile

book_3rdparty_html: versionfile ${CHAP_XML} book_3rdparty.xml
	${XML_CATALOG} ${XSLTPROC_BIN} --nonet --xinclude --output 3rdparty_MIRA.html $(top_srcdir)/doc/docbook/book.xsl $(top_srcdir)/doc/docbook/book_3rdparty.xml

book_3rdparty_pdf: versionfile ${CHAP_XML} book_3rdparty.xml
	-${XML_CATALOG} dblatex --xslt-opts=--nonet --output=3rdparty_MIRA.pdf $(top_srcdir)/doc/docbook/book_3rdparty.xml

book_definitive_html: versionfile ${CHAP_XML} book_definitiveguide.xml
	${XML_CATALOG} ${XSLTPROC_BIN} --nonet --xinclude --output DefinitiveGuideToMIRA.html $(top_srcdir)/doc/docbook/book.xsl $(top_srcdir)/doc/docbook/book_definitiveguide.xml
book_definitive_pdf: versionfile ${CHAP_XML} book_definitiveguide.xml
	-${XML_CATALOG} dblatex --xslt-opts=--nonet --output=DefinitiveGuideToMIRA.pdf $(top_srcdir)/doc/docbook/book_definitiveguide.xml

book_definitive: versionfile ${CHAP_XML} book_definitive_html book_definitive_pdf


books_html: versionfile book_definitive_html book_3rdparty_html

books_pdf: versionfile book_definitive_pdf book_3rdparty_pdf

if HAVE_DBLATEX
books: versionfile books_html books_pdf
else
books: versionfile books_html
endif

all: versionfile books chapters

docs: versionfile all

install-exec-hook: docs
	mkdir -p $(prefix)/share/$(PACKAGE)/docs && \
        cp $(FINALDOCFILES) $(prefix)/share/$(PACKAGE)/docs && \
        cp -r bookfigures $(prefix)/share/$(PACKAGE)/docs && \
        cp -r doccss $(prefix)/share/$(PACKAGE)/docs && \
        cp -r images $(prefix)/share/$(PACKAGE)/docs
