# Help the people at Debian for a simpler and reproducible building
if ISDEBIAN
CITARGET=debian
else
CITARGET=std
endif

AM_CPPFLAGS = -I$(top_srcdir)/src $(all_includes)
AM_CXXFLAGS = $(BOOST_CPPFLAGS)

bin_PROGRAMS= mira
noinst_PROGRAMS= miratest

ciclean:
	-@rm -f compileinfo.itxt

compileinfo-std:
	@/bin/echo -n "Compiled by: " >compileinfo.itxt.xxd
	@whoami >>compileinfo.itxt.xxd
	@date >>compileinfo.itxt.xxd
	@/bin/echo -n "On: " >>compileinfo.itxt.xxd
	@uname -a >>compileinfo.itxt.xxd

compileinfo-debian:
	@/bin/echo -n "Compiled by: " >compileinfo.itxt.xxd
	@$$DEBIAN_MAINTAINER >>compileinfo.itxt.xxd
	@echo $$SOURCE_DATE_EPOCH >>compileinfo.itxt.xxd

compileinfo.itxt: ciclean compileinfo-$(CITARGET)
	@echo "Making .C out of .xxd" $@
	@echo "const char compileinfo2[] = {" >compileinfo.itxt.xxd.C
	@$(SHELL) -ec 'grep -v "^#" compileinfo.itxt.xxd | xxd -i >>compileinfo.itxt.xxd.C'
	@echo ",0 };" >>compileinfo.itxt.xxd.C
	@echo "const char * compileinfo=compileinfo2;" >>compileinfo.itxt.xxd.C
	@rm compileinfo.itxt.xxd

compileinfo.itxt.xxd.C: ciclean compileinfo.itxt

# This tell make that this file is generated, helping it to set dependencies right
# The net effect in this case: binaries always get linked in every make call,
#  therefore including current compile infor text
BUILT_SOURCES=compileinfo.itxt.xxd.C

noinst_HEADERS= quirks.H

AM_LDFLAGS= -L../version -L../modules -L../debuggersupport -L../memorc -L../io/ -L../util/ -L../errorhandling/ -L../mira -L../caf $(BOOST_LDFLAGS) $(MIRASTATICLD)

if BUILDMIRAMEMORC
MEMORC= -lmemorc
endif

MIRALIBS= -lmiraversion -lmiramodules -lmiradebuggersupport $(MEMORC) -lmira -lmiraestass -lmiracaf -lmirautil -lmiradptools -lmirafio -lmiraerrorhandling -lmiradebuggersupport -lmirafmttext



links:
	-rm -f mirabait$(EXEEXT)
	-$(LN_S) mira$(EXEEXT) mirabait$(EXEEXT)
	-rm -f miraconvert$(EXEEXT)
	-$(LN_S) mira$(EXEEXT) miraconvert$(EXEEXT)

install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	rm -f miramem$(EXEEXT) && \
	$(LN_S) mira$(EXEEXT) miramem$(EXEEXT) && \
	rm -f mirabait$(EXEEXT) && \
	$(LN_S) mira$(EXEEXT) mirabait$(EXEEXT) &&\
	rm -f miraconvert$(EXEEXT) && \
	$(LN_S) mira$(EXEEXT) miraconvert$(EXEEXT) && \
	rm -f miramer$(EXEEXT) && \
	$(LN_S) mira$(EXEEXT) miramer$(EXEEXT)

uninstall-hook:
	cd $(DESTDIR)$(bindir) && \
	rm -f miramem$(EXEEXT) && \
	rm -f mirabait$(EXEEXT) && \
	rm -f miraconvert$(EXEEXT) && \
	rm -f miramer$(EXEEXT)

clean-local:
	-rm -f mirabait$(EXEEXT)
	-rm -f miraconvert$(EXEEXT)


strip:
	du -sk $(bin_PROGRAMS)
	strip $(bin_PROGRAMS)
	du -sk $(bin_PROGRAMS)

privinst:
	cp $(bin_PROGRAMS) /home/bach/pub/bin

binclean:
	rm -rf $(bin_PROGRAMS)
	rm -rf $(noinst_PROGRAMS)

mira.C: links
mira_SOURCES= mira.C quirks.C  compileinfo.itxt.xxd.C
mira_LDADD= $(MIRALIBS)

miratest.C:
miratest_SOURCES= miratest.C  compileinfo.itxt.xxd.C
miratest_LDADD= $(MIRALIBS)
