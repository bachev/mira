INCLUDES = -I$(top_srcdir)/src $(all_includes)

clean-local:
	rm -rf *flexer.cc

noinst_LIBRARIES = libmiraEdIt.a
libmiraEdIt_a_FLEXSRC= 
nodist_libmiraEdIt_a_SOURCES= $(libmiraEdIt_a_FLEXSRC:.flex=.flex.C)
libmiraEdIt_a_SOURCES= hypothesen.C parameters.C faultRegion.C parameters_flexer.ll
noinst_HEADERS= hypothesen.H parameters.H parameters_tokens.h



