%option noyywrap
%option c++
%option outfile="lex.yy.c"
%option prefix="CAF"

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
 *
 */

#include <cstdio>

#include "caf_tokens.h"
#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1110
#endif

int intcount=0;
int opa=0;
%}

%option yylineno


ws         [ \t]
blank      [" "]
letter     [A-Za-z]
seqletter  [A-Za-z\-\*]
digit      [0-9]
alphanum   {letter}|{digit}
identifier ({letter}|{digit}|[\.\-\,\_\+\#\|])+
identifierextended ({letter}|{digit}|[\%\&\/\.\-\,\_\+\#\|":"])+
number     ([\-]{digit})|{digit}({digit})*
comment    \/\/
sequence   ^Sequence{ws}*[":"]{ws}*
dna        ^DNA{ws}*[":"]{ws}*
quality    ^BaseQuality{blank}*[":"]{blank}*


%x NUMBER_2_TEXT
%x NUMBER_4
%x IDENT_NUMBER_4
%x TAG_MODE
%x QUOTED_TEXT
%x DNA_MODE
%x DNA_SEQ
%x QUALITY_MODE
%x IDENTIFIER_MODE
%x TEXT
%%

<*>^"Align_to_SCF "    { BEGIN(NUMBER_4);
		         intcount = 0;
			 return(token_align_SCF); }

<*>^"Assembled_from "  { BEGIN(IDENT_NUMBER_4);
		         intcount = 0;
		         return(token_assembled); }

<*>^"Clone_vec "       { BEGIN(NUMBER_2_TEXT);
		         intcount = 0;
		         return(token_clone_vector); }

<*>^"Clipping "        { BEGIN(NUMBER_2_TEXT);
		         intcount = 0;
		         return(token_clipping); }

<*>^"Clone "	       { BEGIN(QUOTED_TEXT);
		         return(token_clone);}

<*>^"Insert_size "     { BEGIN(NUMBER_4);
		         intcount = 2;
		         return(token_insert_size); }

<*>^"SCF_File "        { BEGIN(QUOTED_TEXT);
		         return(token_scf_file);   }

<*>^"Seq_vec "         { BEGIN(NUMBER_2_TEXT);
		         intcount = 0;
		         return(token_seq_vector);   }

<*>^"Sequencing_vector " { BEGIN(QUOTED_TEXT);
		           return(token_sequencing_vector); }

<*>^"Staden_id "       { BEGIN(QUOTED_TEXT);
		         return(token_staden_id);  }

<*>^"Tag "	       { BEGIN(TAG_MODE);
		         intcount = 0;
			 return(token_tag); }

<*>{dna}               { BEGIN(DNA_MODE);
                         return(token_sequence); }

<*>{sequence}          { BEGIN(IDENTIFIER_MODE);
		         return(token_sequencename); }

<*>{quality}           { BEGIN(QUALITY_MODE);
		         return(token_quality);   }

<*>^"Asped "	       { BEGIN(QUOTED_TEXT);
		         return(token_asped);      }

<*>^"Base_caller "     { BEGIN(QUOTED_TEXT);
			 return(token_base_caller);  }

<*>^"Template "        { BEGIN(QUOTED_TEXT);
		         return(token_template);   }

<*>^"Ligation_no "     { BEGIN(QUOTED_TEXT);
		         return(token_ligation);   }

<*>^"ProcessStatus "   { BEGIN(QUOTED_TEXT);
			 return(token_pstatus);    }

<*>^"Primer "          { BEGIN QUOTED_TEXT;
	                 return(token_primer);     }

<*>^"Stolen "	       { BEGIN(QUOTED_TEXT);
			 return(token_stolen);    }



<DNA_MODE>{identifierextended} { BEGIN(DNA_SEQ); return(token_identifier); }

<DNA_SEQ>\n{ws}*\n       { BEGIN(0); return(token_ende); }
<DNA_SEQ>\n              { }
<DNA_SEQ>{seqletter}*    { return(token_identifier); }


<TAG_MODE>{number}           { intcount++;
		               BEGIN(NUMBER_2_TEXT);
		               return(token_number); }
<TAG_MODE>{identifier}       { BEGIN(NUMBER_2_TEXT);
		               return(token_identifier); }

<NUMBER_4>{ws}*
<NUMBER_4>{number}     { if (++intcount == 4) BEGIN(0);
  		         return(token_number);  }


<NUMBER_2_TEXT>{ws}*
<NUMBER_2_TEXT>{number} { if (++intcount == 2) BEGIN(QUOTED_TEXT);
  		          return(token_number);  }
<NUMBER_2_TEXT>QUAL|SVEC|CONT|CVEC



<IDENT_NUMBER_4>{identifierextended}  { BEGIN(NUMBER_4);
			        return(token_identifier); }

<QUOTED_TEXT>[\"]   { BEGIN(TEXT); opa = 1; }
<QUOTED_TEXT>{ws}
<QUOTED_TEXT>\n     { yyless(yyleng-1); BEGIN(0); return(token_quoted_text);}
<QUOTED_TEXT>.      { yymore(); BEGIN(TEXT); opa = 0; }


<TEXT>[\"]   { yymore(); }
<TEXT>\n     { yyless(yyleng-1-opa);
  	       BEGIN(0); return(token_quoted_text); }
<TEXT>.      { yymore(); }


<IDENTIFIER_MODE>{identifierextended} { BEGIN(0); return(token_identifier); }


<QUALITY_MODE>\n{ws}*\n         { BEGIN(0); return(token_ende); }
<QUALITY_MODE>{number}          { return(token_number);     }
<QUALITY_MODE>{identifierextended}      { return(token_identifier); }





<*>^Is_read          { return(token_type_read);   }
<*>^Is_contig        { return(token_type_contig); }
<*>^Is_group         { return(token_type_group);  }
<*>^Is_assembly      { return(token_type_assembly); }
<*>^Padded           { return(token_padded);     }
<*>^Unpadded         { return(token_unpadded);   }

<*>^Dye              { return(token_dye);        }
<*>Dye_terminator   { return(token_dye_terminator); }
<*>Dye_primer       { return(token_dye_primer); }


<*>^Strand{blank}Forward   { return(token_forward); }
<*>^Strand{blank}Reverse   { return(token_reverse); }
<*>{comment}        { while( yyinput() != '\n'); }
<*>\n
<*>\0               { return(token_error); }
<*>.
<<EOF>>             { yyterminate();  }
%%

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1110
#endif
