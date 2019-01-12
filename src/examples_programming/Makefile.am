
AM_CXXFLAGS = $(BOOST_CPPFLAGS)

INCLUDES = -I$(top_srcdir)/src $(all_includes)


noinst_PROGRAMS= mira_101 read_101 readpool_101

AM_LDFLAGS= -L../io/ -L../util/ -L../errorhandling/ -L../mira -L../examine/ -L../EdIt/ -L../caf/ -L../knn_abi373 -L../knn_alf $(BOOST_LDFLAGS)


mira_101.C: 
mira_101_SOURCES= mira_101.C
mira_101_LDADD= -lmira  -lmiraEdIt -lmiraExamine -lmirasupport -lmiraestass -lmiraerrorhandling  -lmirautil -lmiradptools -lmirafio -lmiracaf -lmiraKNN_abi373 -lm -lexpat -lz $(BOOST_THREAD_LIB)  $(BOOST_REGEX_LIB) $(BLOODYDARWINHACK)

read_101.C: 
read_101_SOURCES= read_101.C
read_101_LDADD= -lmirasupport -lmiraerrorhandling  -lmirautil -lmiradptools -lmirafio -lm -lz

readpool_101.C: 
readpool_101_SOURCES= readpool_101.C
readpool_101_LDADD= -lmirasupport -lmiraerrorhandling  -lmirautil -lmiradptools -lmirafio -lm -lz -lexpat



