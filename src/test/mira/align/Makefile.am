INCLUDES = -I$(top_srcdir)/src $(all_includes)

noinst_PROGRAMS= tstAlign
# EdIt_abi373 EdIt_alf

MIRALIBS= -lmiramodules -lmiradebuggersupport $(MEMORC) -lmira -lmiraestass -lmiracaf -lmirautil -lmiradptools -lmirafio -lmiraerrorhandling -lmiradebuggersupport -lmirafmttext $(LINKEDITLIB)

AM_LDFLAGS=  $(BLOODYDARWINSTATLIBDIRHACK) -L../../../modules -L../../../debuggersupport -L../../../memorc -L../../../io/ -L../../../util/ -L../../../errorhandling/ -L../../../mira -L../../../caf $(LINKEDITDIR) $(BOOST_LDFLAGS) -all-static


tstAlign_SOURCES= tstAlign.C
tstAlign_LDADD= $(MIRALIBS)




#convert_project: libinstall convert_project.C
##       $(CPP) $(CCFLAGS) $(INC_DIRS) makeFASTA.C -o makeFASTA  $(LIB) $(LIB_DIR
#S)  -lscflookclips -lClip -lerrorhandling  -lutil -lmyio -lm
#	  $(CPP) $(CCFLAGS) convert_project.C -o convert_project $(INC_DIRS) $(LI
#B) $(LIB_DIRS)  -lerrorhandling  -lutil -lmyio -lcaf -lEdIt -lExamine -lKNN_jena
# -lm
#	  cp convert_project $(TARGET_BINDIR)
