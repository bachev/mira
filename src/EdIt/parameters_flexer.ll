%option noyywrap
%option c++
%option outfile="lex.yy.c"
%option prefix="EP"

%{
/*
 * Written by Thomas Pfisterer
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Thomas Pfisterer
 *
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 */

#include <fstream>
#include "parameters_tokens.h"

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1110
#endif

  //  int myretcode=0;

%}

ANID [A-Za-z][A-Za-z0-9]*
NAME [A-Za-z0-9_#.]+
INT [+\-]?[0-9]+
FLOAT [0-9]+"."[0-9]+

%x GET_UNTIL_NEWLINE
%x TAG_MODE
%x FILE_MODE
%x VERB_MODE
%x CI_FILE_NAME
%x CO_FILE_NAME
%x HO_FILE_NAME
%x TO_FILE_NAME
%x LO_FILE_NAME
%x STRAND_MODE
%x COMMAND_MODE
%%



<*>"-HELP" |
<*>"-h"                           { return EP_HELP; }

<*>"-TAG" |
<*>"-T"                           { BEGIN(TAG_MODE); return yylex();}
<TAG_MODE>"alter_tag" |
<TAG_MODE>"at"                    { return EP_ALTER_TAG;}
<TAG_MODE>"insert_tag" |
<TAG_MODE>"it"                    { return EP_INSERT_TAG;}
<TAG_MODE>"delete_tag" |
<TAG_MODE>"dt"                    { return EP_DELETE_TAG;}
<TAG_MODE>"all"                   { return EP_ALL_TAGS;}
<TAG_MODE>"no"                    { return EP_NO_TAGS;}
<TAG_MODE>"delete_asterisk" |
<TAG_MODE>"ast"		          { return EP_DELETE_ASTERISK; }
<TAG_MODE>"consensus" |
<TAG_MODE>"co"                    { return EP_CONSENSUS_TAG; }

<*>"-FILE" |
<*>"-F"                           { BEGIN(FILE_MODE); return yylex();} 
<FILE_MODE>"caf" | 
<FILE_MODE>"c"			  { BEGIN(CI_FILE_NAME); }
<FILE_MODE>"out" |
<FILE_MODE>"o"			  { BEGIN(CO_FILE_NAME); }
<FILE_MODE>"log" |
<FILE_MODE>"l"                    { BEGIN(LO_FILE_NAME); }
<FILE_MODE>"html" |
<FILE_MODE>"h"                    { BEGIN(HO_FILE_NAME); }
<FILE_MODE>"text" |
<FILE_MODE>"t"                    { BEGIN(TO_FILE_NAME); }
<FILE_MODE>"void" |
<FILE_MODE>"v"                    { return EP_outfile_void; }


<CI_FILE_NAME>{NAME}    { BEGIN(FILE_MODE); return EP_caf_in_file;}
<CO_FILE_NAME>{NAME}	{ BEGIN(FILE_MODE); return EP_caf_out_file; }
<HO_FILE_NAME>{NAME}    { BEGIN(FILE_MODE); return EP_html_out_file;}
<TO_FILE_NAME>{NAME}    { BEGIN(FILE_MODE); return EP_text_out_file;}
<LO_FILE_NAME>{NAME}	{ BEGIN(FILE_MODE); return EP_log_out_file;}	


<*>"-VERB" |
<*>"-V"                      { BEGIN(VERB_MODE); return yylex(); }
<VERB_MODE>"v" |
<VERB_MODE>"verbose"	     { return EP_VERBOSE; }
<VERB_MODE>"d" |
<VERB_MODE>"debug"           { return EP_DEBUG; }


<*>"-STRAND" |
<*>"-S"                      { BEGIN(STRAND_MODE); return yylex(); }
<STRAND_MODE>"double" |
<STRAND_MODE>"d"             { return EP_make_double; }
<STRAND_MODE>"score" |  
<STRAND_MODE>"s"             { return EP_doublestrand_minscore; }   
<STRAND_MODE>"cover" |
<STRAND_MODE>"c"             { return EP_mincoverage_single; }


<*>"-DO" |
<*>"-D"                      { BEGIN(COMMAND_MODE); return yylex(); }
<COMMAND_MODE>"all"          { return EP_DO_ALL;       }
<COMMAND_MODE>"nop"          { return EP_DO_NOP;       }
<COMMAND_MODE>"hypo"         { return EP_DO_HYPO;      }
<COMMAND_MODE>"eval"         { return EP_DO_EVAL;      }
<COMMAND_MODE>"contig"       { return EP_DO_CONTIG;    }
<COMMAND_MODE>"range"        { return EP_DO_RANGE;     }
<COMMAND_MODE>"strict_on"    { return EP_DO_STRICT;    }
<COMMAND_MODE>"strict_off"   { return EP_DO_NOSTRICT;  }		     
<COMMAND_MODE>"low"          { return EP_DO_LOWQUAL;   }
<COMMAND_MODE>"threshold" |
<COMMAND_MODE>"t"            { return EP_DO_THRESHOLD; }
<COMMAND_MODE>"col" |
<COMMAND_MODE>"column"       { return EP_DO_COLUMNS;   }
<COMMAND_MODE>"reg" |
<COMMAND_MODE>"region"       { return EP_DO_REGIONS;   }
<COMMAND_MODE>"small"        { return EP_DO_SMALL_REGIONS;    }

<*>"-=BEGIN0=-" {BEGIN(0);}     /* Hack for parsing command lines
                                   where each blank resets the state */
<*>{ANID}  {return EP_ANID;}
<*>{FLOAT} {return EP_FLOAT;}
<*>{INT}   {return EP_INT;}

<*>":"
<*>"="

<*>[\t\n ]   {/* munch these, if not recognised earlier */ }

%%

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1110
#endif


