INCLUDES = -I$(top_srcdir)/src $(all_includes)

noinst_PROGRAMS= tstRead
# EdIt_abi373 EdIt_alf

AM_LDFLAGS= -L../../../io/ -L../../../util/ -L../../../errorhandling/ -L../../../mira -L../../../examine/ -L../../../EdIt/ -L../../../caf/ -L../../../knn_abi373 -L../../../knn_alf 

tstRead_SOURCES= tstRead.C
tstRead_LDADD= -lmiraerrorhandling  -lmirautil -lmira -lmirasupport -lmiradptools -lmirafio

