INCLUDES = -I$(top_srcdir)/src $(all_includes)

noinst_PROGRAMS= tstgeneralio

AM_LDFLAGS= -L../../../io/ -L../../../util/ -L../../../errorhandling/ -L../../../mira -L../../../examine/ -L../../../EdIt/ -L../../../caf/ -L../../../knn_abi373 -L../../../knn_alf 

tstgeneralio_SOURCES= tstgeneralio.C
tstgeneralio_LDADD= -lmiraerrorhandling  -lmirautil -lmirafio



#convert_project: libinstall convert_project.C
##       $(CPP) $(CCFLAGS) $(INC_DIRS) makeFASTA.C -o makeFASTA  $(LIB) $(LIB_DIR
#S)  -lscflookclips -lClip -lerrorhandling  -lutil -lmyio -lm
#	  $(CPP) $(CCFLAGS) convert_project.C -o convert_project $(INC_DIRS) $(LI
#B) $(LIB_DIRS)  -lerrorhandling  -lutil -lmyio -lcaf -lEdIt -lExamine -lKNN_jena
# -lm
#	  cp convert_project $(TARGET_BINDIR)
