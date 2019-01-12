%option noyywrap
%option c++
%option outfile="lex.yy.c"
%option prefix="EXP"

%{
//#include <sstream>
#include <fstream>
#include "tokens.h"

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1110
#endif

  //  int myretcode=0;
  int tag_intcounter=0;

//<GET_UNTIL_EOSQ>[AaCcGgTtNn\-\n ] { yymore();}
//<GET_UNTIL_EOSQ>[XRYMSKWHBVDxrymskwhbvd*]  { yymore();}
//<GET_UNTIL_EOSQ>^"//" {BEGIN(0); return EXPT_EOSQ;}
//<GET_UNTIL_EOSQ>.     {return EXPT_ILLBASEINSQ;}


%}

ANID [A-Za-z\_][A-Za-z0-9\_\.]*
INT [0-9]+
FLOAT [0-9]+"."[0-9]+

%x GET_UNTIL_NEWLINE
%x GET_UNTIL_EOSQ
%x AV_MODE
%x TAG_MODE
%x TAG_MODE_EOL
%x ON_MODE1
%x ON_MODE2
%%

<<EOF>>  { BEGIN(0); yyterminate();}

^"SF   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_SF; }

^"TG   "/[^" "]  { BEGIN(TAG_MODE); tag_intcounter=0; return EXPT_TG; }
<TAG_MODE>[\+\-=] { return EXPT_PME;}
<TAG_MODE>".."    {return EXPT_PP;}
<TAG_MODE>{ANID}  {return EXPT_ANID;}
<TAG_MODE>{INT}   {
                    if(++tag_intcounter==2) BEGIN(TAG_MODE_EOL);
                    return EXPT_INT;
		  }
<TAG_MODE>" "
<TAG_MODE>\n      {BEGIN(0); return EXPT_EOL;}
<TAG_MODE_EOL>.*/\n {BEGIN(0); return EXPT_TAGONELINER;}
^"TG   "" "{5} { BEGIN(GET_UNTIL_NEWLINE); return EXPT_TGC; }
^"SQ   " { BEGIN(GET_UNTIL_EOSQ); return EXPT_SQ;}
^"SQ"" "* { BEGIN(GET_UNTIL_EOSQ); return EXPT_SQ;}
<GET_UNTIL_EOSQ>[AaCcGgTtNn\- XRYMSKWHBVDxrymskwhbvd*]* {return EXPT_SQseq;}
<GET_UNTIL_EOSQ>\n     { }
<GET_UNTIL_EOSQ>"//"    {BEGIN(0); return EXPT_EOSQ;}
<GET_UNTIL_EOSQ>.     {return EXPT_ILLBASEINSQ;}

^"AV   "        { BEGIN(AV_MODE); return EXPT_AV;}
<AV_MODE>" "
<AV_MODE>{INT}  { return EXPT_INT;}
<AV_MODE>\n     { BEGIN(0); return EXPT_EOL;}

^"ON   "        { BEGIN(ON_MODE1); return EXPT_ON;}
<ON_MODE1>" "
<ON_MODE1>{INT}".."  { BEGIN(ON_MODE2); return EXPT_RANGELOW;}
<ON_MODE1>{INT}  { return EXPT_INT;}
<ON_MODE1>\n     { BEGIN(0); return EXPT_EOL;}
<ON_MODE2>{INT}  { BEGIN(ON_MODE1); return EXPT_RANGEHIGH;}
<ON_MODE2>\n     { BEGIN(0); return EXPT_EOL;}

^"CF   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_CF;}
^"CN   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_CN;}
^"CV   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_CV;}
^"DT   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_DT;}
^"EN   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_EN;}
^"ID   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_ID;}
^"LE   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_LE;}
^"LI   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_LI;}
^"LN   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_LN;}
^"LT   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_LT;}
^"MA   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_MT;}
^"MC   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_MT;}
^"MT   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_MT;}
^"MN   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_MN;}
^"PS   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_PS;}
^"SV   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_SV;}
^"SS   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_SS;}
^"TN   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_TN;}
^"OP   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_OP;}
^"BC   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_BC;}

^"CS   " { return EXPT_CS;}
^"SI   " { return EXPT_SI;}

^"AQ   " { return EXPT_AQ;}
^"SL   " { return EXPT_SL;}
^"SR   " { return EXPT_SR;}
^"QL   " { return EXPT_QL;}
^"QR   " { return EXPT_QR;}
^"CL   " { return EXPT_CL;}
^"CR   " { return EXPT_CR;}
^"CH   " { return EXPT_CH;}

^"PR   " { return EXPT_PR;}
^"SC   " { return EXPT_SC;}
^"SP   " { return EXPT_SP;}

^"Blast database" {
                     // error in U13a02a09.p1
                     BEGIN(GET_UNTIL_NEWLINE); yylex(); return EXPT_JENAQUIRK;
                  }

^[A-Z][A-Z]"   " { BEGIN(GET_UNTIL_NEWLINE); return EXPT_UNKNOWN;}

<GET_UNTIL_NEWLINE>.*/\n { BEGIN(0); return EXPT_EOL;}


".."    {return EXPT_PP;}
"//"    {return EXPT_SS;}
{ANID}  {return EXPT_ANID;}
{FLOAT} {return EXPT_FLOAT;}
{INT}   {return EXPT_INT;}
"\-"    {return EXPT_MINUS;}
"\+"    {return EXPT_PLUS;}

[\n ]   {/* munch these, if not recognised earlier */ }

%%

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1110
#endif

//<GET_UNTIL_EOSQ>[AaCcGgTtNn\-\n XRYMSKWHBVDxrymskwhbvd*]*  { yymore();}
