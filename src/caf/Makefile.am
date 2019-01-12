AM_CPPFLAGS = -I$(top_srcdir)/src $(all_includes)

clean-local:
	rm -rf *flexer.cc

noinst_LIBRARIES = libmiracaf.a
libmiracaf_a_SOURCES= caf.C caf_flexer.ll
noinst_HEADERS= caf.H caf_tokens.h
