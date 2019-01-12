
EXTRA_DIST = mira-install-sls-rrna.sh README.txt rfam_rrna-21-1.sls.gz


SUBDIRS = \
	mira-createsls

MHSDIR=$(prefix)/share/mira/mhs

install-exec-hook:
	cp mira-install-sls-rrna.sh $(DESTDIR)$(bindir) && \
	mkdir -p $(MHSDIR) && \
	$(DESTDIR)$(bindir)/mira-install-sls-rrna.sh $(MHSDIR) rfam_rrna-21-1.sls.gz && \
	cd $(MHSDIR) && \
	$(LN_S) -f rfam_rrna-21-1.mhs.gz filter_default_rrna.mhs.gz

uninstall-hook:
	rm $(DESTDIR)$(bindir)/mira-install-sls-rrna.sh && \
	rm -f $(prefix)/share/mira
