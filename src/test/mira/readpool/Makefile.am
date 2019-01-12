INCLUDES = -I$(top_srcdir)/src $(all_includes)

noinst_PROGRAMS= tstReadPool
# EdIt_abi373 EdIt_alf

AM_LDFLAGS= -L../../../io/ -L../../../util/ -L../../../errorhandling/ -L../../../mira -L../../../examine/ -L../../../EdIt/ -L../../../caf/ -L../../../knn_abi373 -L../../../knn_alf 

tstReadPool_SOURCES= tstReadPool.C
tstReadPool_LDADD= -lmiraerrorhandling  -lmira -lmirasupport -lmiradptools -lmirautil -lmirafio -lz -lexpat

