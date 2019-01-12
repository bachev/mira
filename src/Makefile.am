GITBUILD=`$(top_srcdir)/build-aux/git-version-gen $(top_srcdir)/.tarball-version`

EXTRA_DIST = version.H

SUBDIRS = \
	debuggersupport\
	io\
	stdinc\
	errorhandling\
	util\
	mira\
	caf\
	memorc\
	modules\
	support\
	scripts\
	version\
	progs\
	other

# some botched test install might create a "share" here
# get rid of it
clean-local:
	-rm -rf share
