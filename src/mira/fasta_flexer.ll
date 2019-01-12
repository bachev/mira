%option noyywrap
%option prefix="FA"

%{
/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2000 and later by Bastien Chevreux
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

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1110
#endif

  //  int myretcode=0;

%}

ANID [A-Za-z][A-Za-z0-9_\.\-]*
DNASEQ [ACGTNacgtn-]*

%x SCANFORNAME_MODE
%x SCANDNA_MODE
%x SEARCHDNA_MODE
%x GGG
%%

^> {BEGIN(SCANFORNAME_MODE);}
<SCANFORNAME_MODE>{ANID} { BEGIN(SEARCHDNA_MODE);
                           cout << "seenq: ";ECHO; cout << endl;
			   return 1;
			 }
<SEARCHDNA_MODE>.*\n {cout << "scanning now "; ECHO; BEGIN(GGG);}
<GGG>[AaCcGgTtNnXRYMSKWHBVDxrymskwhbvd\-\n ]* { yymore();}
<GGG>^">" {BEGIN(SCANFORNAME_MODE); yyleng-=10; return 0;}
<GGG>.     {cout << "unknown character in fasta sequence:" << yytext<<endl;return 666;}

<*>[\t\n ]   {}

%%

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1110
#endif
