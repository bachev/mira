AM_CPPFLAGS = -I$(top_srcdir)/src $(all_includes)

AM_LDFLAGS= -L../io/ -L../util/ -L../errorhandling/ -L../mira -L../caf/ $(BOOST_LDFLAGS) $(MIRASTATICLD)

if BUILDMIRAMEMORC
MEMORC= -L. -lmemorc
endif

noinst_LIBRARIES = libmemorc.a
libmemorc_a_SOURCES= memorc.C
noinst_HEADERS= memorc.H

noinst_PROGRAMS= main

main.C: memorc.o
main_SOURCES= main.C
main_LDADD= $(MEMORC) -lmiraerrorhandling

#main_LDADD= -L../memorc -lmemorc -L../util/ -lmirautil -L../errorhandling/ -lmiraerrorhandling -L../mira/ -lmira -lmirasupport -lpreventinitfiasco
