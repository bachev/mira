
GITBUILD=`$(top_srcdir)/build-aux/git-version-gen $(top_srcdir)/.tarball-version`

BUILT_SOURCES = mversion.H
EXTRA_DIST = mversion.stub

noinst_LIBRARIES = libmiraversion.a


libmiraversion_a_SOURCES= \
	mversion.C

noinst_HEADERS = mversion.H

target_not_existing_always_run:
	@echo >/dev/null

#make sure mversion.H gets made before mversion.C or else things break
mversion.C: mversion.H

mversion.H: mversion.stub target_not_existing_always_run
	@cp mversion.stub mversion.H
	@echo \#define GITBUILD \"${GITBUILD}\" >>mversion.H
	@echo \#define MIRALIBMVERSION GITBUILD >>mversion.H
	@echo \#define MIRAMVERSION MIRALIBMVERSION \" \" MLV_QUIET MLV_TIME\"\" >>mversion.H


clean-local:
	-@rm mversion.H
