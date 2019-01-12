/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2005 and later by Bastien Chevreux
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


#include "util/dptools.H"

#include <cstdlib>     // exit() under Cygwin

using std::cout;
using std::cerr;
using std::cin;
using std::endl;

#define CEBUG(bla)


/*
 * nsTranslationTables is needed to translate a codon triplet into an AA code
 *  and back
 *
 * from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1
 *
 * description of fields:
 * for each translation table: 5*64 characters (+ 0 at end)
 * 1st line:  AA code
 * 2nd line:  whether start codon (M) or not (-)
 * 3rd to 5th line: base 1 to 3 of the codon
 */

char dptools::nsTranslationTables[23][321]={
// transl_table=1
  "\
FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\
---M---------------M---------------M----------------------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Vertebrate Mitochondrial Code (transl_table=2)
  "\
FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG\
--------------------------------MMMM---------------M------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// The Yeast Mitochondrial Code (transl_table=3)
  "\
FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG\
----------------------------------MM----------------------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
  "\
FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\
--MM---------------M------------MMMM---------------M------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Invertebrate Mitochondrial Code (transl_table=5)
  "\
FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG\
---M----------------------------MMMM---------------M------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
  "\
FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\
-----------------------------------M----------------------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Table 7 phased out and deleted, historical
  "\
FFLLSSSSYY**CC*WLLLLPPPPHHAARRRRIIIMTTTTNNKKSSIIIVVVHHHHCCAABBBB\
---M---------------M------------------------------MM------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Table 8 phased out and deleted, historical
  "\
FFLLSSSSYY**AA*RIIIIPPPPRRFFTTTTHHHGRRRRYYPPCCRRVVVVAAAADDEEGGGG\
---M----------------MMMM---------------------------M------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Echinoderm and Flatworm Mitochondrial Code (transl_table=9)
  "\
FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG\
-----------------------------------M---------------M------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Euplotid Nuclear Code (transl_table=10)
  "\
FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\
-----------------------------------M----------------------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Bacterial and Plant Plastid Code (transl_table=11)
  "\
FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\
---M---------------M------------MMMM---------------M------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Alternative Yeast Nuclear Code (transl_table=12)
  "\
FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\
-------------------M---------------M----------------------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Ascidian Mitochondrial Code (transl_table=13)
  "\
FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG\
---M------------------------------MM---------------M------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Alternative Flatworm Mitochondrial Code (transl_table=14)
  "\
FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG\
-----------------------------------M----------------------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Blepharisma Nuclear Code (transl_table=15)
  "\
FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\
-----------------------------------M----------------------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Chlorophycean Mitochondrial Code (transl_table=16)
  "\
FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\
-----------------------------------M----------------------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Table 17 not defined
  "\
",
// Table 18 not defined
  "\
",
// Table 19 not defined
  "\
",
// Table 20 not defined
  "\
",
// Trematode Mitochondrial Code (transl_table=21)
  "\
FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG\
-----------------------------------M---------------M------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Scenedesmus obliquus mitochondrial Code (transl_table=22)
  "\
FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\
-----------------------------------M----------------------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
",
// Thraustochytrium Mitochondrial Code (transl_table=23)
  "\
FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\
--------------------------------M--M---------------M------------\
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG\
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG\
"
};


// MIRA always uses "ACGT*" when building groups etc.
// From index to base is esy (e.g. "ACGT*"[2]), but
//  the reverse is easiest with this table and
//  the e.g. dptools::getIndexOfBase('a') inline function.
uint8 dptools::nsbase2index_translationtable[256]=
{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,  // 0x41 = A
  0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,   // 0x61 = a
  0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

char dptools::nsvalid_ACGTbases[256]=
{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 'A', 0, 'C', 0, 0, 0, 'G', 0, 0, 0, 0, 0, 0, 0, 0,  // 0x41 = A
  0, 0, 0, 0, 'T', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 'a', 0, 'c', 0, 0, 0, 'g', 0, 0, 0, 0, 0, 0, 0, 0,   // 0x61 = a
  0, 0, 0, 0, 't', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};


char dptools::nsvalid_bases[256]=
{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 'A', 0, 'C', 0, 0, 0, 'G', 0, 0, 0, 0, 0, 0, 'N', 0,  // 0x41 = A
  0, 0, 0, 0, 'T', 0, 0, 0, 'X', 0, 0, 0, 0, 0, 0, 0,
  0, 'a', 0, 'c', 0, 0, 0, 'g', 0, 0, 0, 0, 0, 0, 'n', 0,   // 0x61 = a
  0, 0, 0, 0, 't', 0, 0, 0, 'x', 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};



// '@' is a base with no coverage in MIRA, not strictly a valid
//  base, but the complement of it is ... surprise, '@'
char dptools::nscomplement_bases[256]=
{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '*', 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  '@', 'T', 0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0,  // 0x41 = A
  0, 0, 0, 0, 'A', 0, 0, 0, 'X', 0, 0, 0, 0, 0, 0, 0,
  0, 't', 0, 'g', 0, 0, 0, 'c', 0, 0, 0, 0, 0, 0, 'n', 0,   // 0x61 = a
  0, 0, 0, 0, 'a', 0, 0, 0, 'x', 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

char dptools::nsvalidIUPAC_bases[256]=
{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 'A', 'B', 'C', 'D', 0, 0, 'G', 'H', 0, 0, 'K', 0, 'M', 'N', 0,  // 0x41 = A
  0, 0, 'R', 'S', 'T', 0, 'V', 'W', 'X', 'Y', 0, 0, 0, 0, 0, 0,
  0, 'a', 'b', 'c', 'd', 0, 0, 'g', 'h', 0, 0, 'k', 0, 'm', 'n', 0,   // 0x61 = a
  0, 0, 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};


//Table 2. Definition of complementary DNA symbols
//  (*) some complements are equal to the original
//
//   Symbol     A B C D G H K M R S  T/U V W  Y     N  X
//   Complement T V G H C D M K Y S*  A  B W* R     N* X*

char dptools::nscomplementIUPAC_bases[256]=
{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '*', 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  '@', 'T', 'V', 'G', 'H', 0, 0, 'C', 'D', 0, 0, 'M', 0, 'K', 'N', 0,  // 0x41 = A
  0, 0, 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 0, 0, 0, 0, 0, 0,
  0, 't', 'v', 'g', 'h', 0, 0, 'c', 'd', 0, 0, 'm', 0, 'k', 'n', 0,   // 0x61 = a
  0, 0, 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};


// A C G T
// 3 2 1 0
// 8 4 2 1

//        0
//  T/U   1   T
//  G     2   G
//  K     3   GT
//  C     4   C
//  Y     5   CT
//  S     6   CG
//  B     7   CGT
//  A     8   A
//  W     9   AT
//  R     a   AG
//  D     b   AGT
//  M     c   AC
//  H     d   ACT
//  V     e   ACG
//  N     f   ACGT

// X=N


// N and X not present in this table, they'd kill all
uint8 dptools::nsIUPAC_basebitmask[256]=
{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 8, 7, 4, 0xb, 0, 0, 2, 0xd, 0, 0, 3, 0, 0xc, 0, 0,  // 0x41 = A
  0, 0, 0xa, 6, 1, 1, 0xe, 9, 0, 5, 0, 0, 0, 0, 0, 0,
  0, 8, 7, 4, 0xb, 0, 0, 2, 0xd, 0, 0, 3, 0, 0xc, 0, 0,   // 0x61 = a
  0, 0, 0xa, 6, 1, 1, 0xe, 9, 0, 5, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// like above, but N is present
//  (X not, as one realy does not know what it could be)
uint8 dptools::nsIUPAC_basebitmaskN[256]=
{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 8, 7, 4, 0xb, 0, 0, 2, 0xd, 0, 0, 3, 0, 0xc, 0xf, 0,  // 0x41 = A
  0, 0, 0xa, 6, 1, 1, 0xe, 9, 0, 5, 0, 0, 0, 0, 0, 0,
  0, 8, 7, 4, 0xb, 0, 0, 2, 0xd, 0, 0, 3, 0, 0xc, 0xf, 0,   // 0x61 = a
  0, 0, 0xa, 6, 1, 1, 0xe, 9, 0, 5, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};


// TODO: check if [0] shouldn't be a '*'
char dptools::nsIUPACbitmask[16]=
{
  'N',
  'T',
  'G',
  'K',
  'C',
  'Y',
  'S',
  'B',
  'A',
  'W',
  'R',
  'D',
  'M',
  'H',
  'V',
  'N'
};


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

uint8 dptools::getIndexOfCodon(const char b1, const char b2, const char b3)
{
  uint8 retvalue=64; // 64 == invalid characters encountered (non-ACGT)
  if(isValidACGTBase(b1) && isValidACGTBase(b2) && isValidACGTBase(b3)){
    retvalue=getIndexOfBase(b1) << 4;
    retvalue+=getIndexOfBase(b2) << 2;
    retvalue+=getIndexOfBase(b3);
  }
  return retvalue;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void dptools::codon2AminoAcids_wrapped(const uint8 transl_table, char base1, char base2, char base3, std::vector<std::string> & codonvariants, std::vector<char> & aaresult, std::vector<bool> & isstartresult)
{
  bool hadIUPAC=false;
  char * aacode=&nsTranslationTables[0][0]+transl_table*321;

  if(base1=='X' || base1=='@') base1='N';
  if(base2=='X' || base2=='@') base2='N';
  if(base3=='X' || base3=='@') base3='N';

  if(base1=='N' && base2=='N' && base3=='N'){
    // this is a shortcut for the case 'nnn' or '@@@' or 'n@n' etc.
    aaresult.push_back('X');
    isstartresult.push_back(false);
    return;
  }

  // the switches are all fall throughs till T
  switch(toupper(base1)) {
  case 'G': aacode+=16;
  case 'A': aacode+=16;
  case 'C': aacode+=16;
  case 'T': break;
  default: {
    std::vector<char> v;
    if(isValidIUPACBase(base1)){
      v=getNucleicAcidFromIUPAC(base1);
    }else{
      v.push_back('A');
      v.push_back('C');
      v.push_back('G');
      v.push_back('T');
    }
    if(v.size()==0) return;
    for(uint8 i=0; i<v.size(); i++){
      codon2AminoAcids_wrapped(transl_table,
			       v[i],
			       base2,
			       base3,
			       codonvariants,
			       aaresult,
			       isstartresult);
    }
    hadIUPAC=true;
  }
  }

  switch(toupper(base2)) {
  case 'G': aacode+=4;
  case 'A': aacode+=4;
  case 'C': aacode+=4;
  case 'T': break;
  default: {
    std::vector<char> v;
    if(isValidIUPACBase(base2)){
      v=getNucleicAcidFromIUPAC(base2);
    }else{
      v.push_back('A');
      v.push_back('C');
      v.push_back('G');
      v.push_back('T');
    }
    if(v.size()==0) return;
    for(uint8 i=0; i<v.size(); i++){
      codon2AminoAcids_wrapped(transl_table,
			       base1,
			       v[i],
			       base3,
			       codonvariants,
			       aaresult,
			       isstartresult);
    }
    hadIUPAC=true;
  }
  }

  switch(toupper(base3)) {
  case 'G': aacode+=1;
  case 'A': aacode+=1;
  case 'C': aacode+=1;
  case 'T': break;
  default: {
    std::vector<char> v;
    if(isValidIUPACBase(base3)){
      v=getNucleicAcidFromIUPAC(base3);
    }else{
      v.push_back('A');
      v.push_back('C');
      v.push_back('G');
      v.push_back('T');
    }
    if(v.size()==0) return;
    for(uint8 i=0; i<v.size(); i++){
      codon2AminoAcids_wrapped(transl_table,
			       base1,
			       base2,
			       v[i],
			       codonvariants,
			       aaresult,
			       isstartresult);
    }
    hadIUPAC=true;
  }
  }

  // if there was a IUPAC in one of the bases, it has been resolved
  //  and added in the recursive loops
  if(hadIUPAC) return;

  // add the codon as string to the variants
  codonvariants.resize(codonvariants.size()+1);
  codonvariants.back().push_back(base1);
  codonvariants.back().push_back(base2);
  codonvariants.back().push_back(base3);

  // the following code is to make sure that the AA and the status whether
  //  it's a start or not is only added to the array if it is not already in
  const char aafound=*aacode;
  aacode+=64;
  bool aaisstart=false;
  if(*aacode=='M') {
    aaisstart=true;
  }

  uint8 i=0;
  bool alreadyinarray=false;
  for(;i<aaresult.size();i++){
    if(aaresult[i]==aafound
       && aaisstart==isstartresult[i]){
      alreadyinarray=true;
      break;
    }
  }

  if(!alreadyinarray){
    aaresult.push_back(aafound);
    isstartresult.push_back(aaisstart);
  }

  return;
}



/*************************************************************************
 *
 * given a translation table and a DNA codon (base123), give back in a vector
 *  the amino acid(s) this codon translates to and whether the|each
 *  AA is a start AA or not
 *
 * Each base of the codon can be a IUPAC base (contain ambiguities)
 *
 *************************************************************************/

void dptools::codon2AminoAcids(const uint8 transl_table, const char base1, const char base2, const char base3, std::vector<std::string> & codonvariants, std::vector<char> & aaresult, std::vector<bool> & isstartresult)
{
  codonvariants.clear();
  aaresult.clear();
  isstartresult.clear();

  bool tttableok=false;
  if(transl_table<=23){
    switch(transl_table){
    case 0:
    case 7:
    case 8:
    case 17:
    case 18:
    case 19:
    case 20: break;
    default: {
      tttableok=true;
      // decrease transl_table by one to adapt to the c array
      codon2AminoAcids_wrapped(transl_table-1,
			       static_cast<char>(toupper(base1)),
			       static_cast<char>(toupper(base2)),
			       static_cast<char>(toupper(base3)),
			       codonvariants,
			       aaresult,
			       isstartresult);
      // not the slightest idea, why gcc 4.3.2 says:
      //  warning: conversion to ‘uint8’ from ‘int’ may alter its value [-Wconversion]
    }
    }

    // if aaresult.size() is >1: check if all results are the same and if yes, collapse
    // e.g.: ATH -> III = I
    if(aaresult.size()>1){
      auto checkbase=aaresult.front();
      bool allsame=true;
      for(auto x : aaresult){
	if(x!=checkbase){
	  allsame=false;
	  break;
	}
      }
      if(allsame){
	aaresult.resize(1);
	isstartresult.resize(1);
      }
    }
  }
  if(!tttableok){
    cout << "\n\nInternal error, unknown translation table " << static_cast<uint16>(transl_table) << endl;
    exit(1);
  }

  return;
}

/*************************************************************************
 *
 * dnaToProtein has become a misnomer: it also extracts the DNA from
 *  a larger piece of DNA. Rework.
 *
 * from a piece of DNA (which might contain gap characters), extract
 *  a given piece (store it with no gaps, in 5' -> 3' direction),
 *  translate it into protein (also stored).
 *
 * if musttranslate && orfend is == orfstart, then extract according to protein
 *  translation until a terminator codon is found (translated to '*'),
 *  else extract strictly between orfstart and orfend
 *
 *************************************************************************/

//#define CEBUG(bla)  {cout << bla; cout.flush();}
void dptools::dnaToProtein(const std::string & dna, std::string & protresult, std::string & dnaresult, uint32 orfstart, uint32 orfend, const int8 orfdir, const uint8 transl_table, const uint8 codonstart, bool musttranslate)
{
  //CEBUG("dna: " << dna << endl);
  CEBUG("orfstart   : " << orfstart << endl);
  CEBUG("orfdir     : " << static_cast<int16>(orfdir) << endl);
  CEBUG("codonstart : " << static_cast<uint16>(codonstart) << endl);
  CEBUG("transl_table : " << static_cast<uint16>(transl_table) << endl);

  if(codonstart<1 || codonstart>3) {
    cerr << "Internal error dnaToProtein(): codonstart " << static_cast<uint16>(codonstart) << " != {1,2,3}\n";
    exit(10);
  }
  if(transl_table>=23) {
    cerr << "Internal error dnaToProtein(): transl_table " << static_cast<uint16>(transl_table) << " >= 23\n";
    exit(10);
  }

  protresult.clear();
  dnaresult.clear();

  char codon[3];
  std::vector<char> aa;
  std::vector<bool> startres;
  std::vector<std::string> codonvariants; // unused, but codon2AminoAcids() needs it

  bool stopatstopcodon=(orfstart==orfend);

  if(orfdir>=0){
    if(orfend<orfstart) std::swap(orfstart, orfend);
    if(musttranslate && orfstart==orfend){
      orfend=static_cast<uint32>(dna.size());
    }

    for(uint32 pos=orfstart+codonstart-1; pos<orfend && pos<dna.size();){
      CEBUG("\npos1: " << pos);

      aa.clear();
      startres.clear();

      uint8 nuc=0;
      for(; nuc<3 && pos < orfend; ++pos){
	CEBUG("  pos2: " << pos);
	if(dna[pos]!='*'){
	  CEBUG(dna[pos]);
	  codon[nuc]=dna[pos];
	  dnaresult.push_back(dna[pos]);
	  ++nuc;
	}
      }
      CEBUG(" // ");
      if(musttranslate && nuc == 3){
	CEBUG(" c2aa ");
	CEBUG("codon: " << static_cast<char>(codon[0]) << static_cast<char>(codon[1]) << static_cast<char>(codon[2]));
	codon2AminoAcids(transl_table,
			 codon[0],
			 codon[1],
			 codon[2],
			 codonvariants,
			 aa,
			 startres);
	// aa vector with 0 or more than 1 entry means: oooops
	if(aa.empty()){
	  cout << "\naa.empty() 1 ???" << endl;
	  exit(1);
	}
	if(aa.size()!=1) {
	  protresult+='X';
	  CEBUG("\taa: X\n");
	}else{
	  CEBUG("\taa: " << aa[0] << endl);
	  protresult+=aa[0];

	  // stop translation on stop
	  if(aa[0]=='*' && stopatstopcodon) orfend=0;;
	}
      }
    }
  }else{
    if(orfstart<orfend) std::swap(orfstart, orfend);

    for(int32 pos=orfstart-codonstart+1; pos>=0 && pos>=static_cast<int32>(orfend);){
      if(musttranslate&& orfstart==orfend){
	orfend=0;
      }
      CEBUG("\npos1: " << pos);

      aa.clear();
      startres.clear();

      uint8 nuc=0;
      for(; nuc<3 && pos>=0; --pos){
	CEBUG("pos2: " << pos << "   ");
	if(dna[pos]!='*'){
	  codon[nuc]=getComplementIUPACBase(dna[pos]);
	  dnaresult.push_back(codon[nuc]);
	  ++nuc;
	}
      }
      CEBUG(" \\ ");
      if(musttranslate && nuc == 3) {
	CEBUG("codon: " << static_cast<char>(codon[0]) << static_cast<char>(codon[1]) << static_cast<char>(codon[2]));

	codon2AminoAcids(transl_table,
			 codon[0],
			 codon[1],
			 codon[2],
			 codonvariants,
			 aa,
			 startres);
	if(aa.empty()){
	  cout << "\naa.empty() 2 ???" << endl;
	  exit(1);
	}
	if(aa.size()!=1) {
	  protresult+='X';
	  CEBUG("\taa: X\n");
	}else{
	  CEBUG("   aa: " << aa[0] << endl);
	  protresult+=aa[0];
	  if(aa[0]=='*' && stopatstopcodon) orfend=static_cast<uint32>(dna.size());
	}
      }
    }
  }
}
//#define CEBUG(bla)




/*************************************************************************
 *
 * works only, when whichdnapos does not hit a '*'
 *
 *
 *************************************************************************/


//#define CEBUG(bla)  {cout << bla; cout.flush();}
void dptools::infoOnAAatDNAPos(const std::string & dna, const uint32 whichdnapos, const uint32 orfstart, const int8 orfdir, const uint8 transl_table, const uint8 codonstart, std::string & codon, std::vector<std::string> & codonvariants, std::vector<char> & aminoacid, std::vector<bool> & isstart, int32 & aanumber, int8 & posinaa)
{
  if(codonstart<1 || codonstart>3) {
    cerr << "Internal error infoOnAAatDNAPos(): codonstart " << static_cast<uint16>(codonstart) << " != {1,2,3}\n";
    exit(10);
  }
  if(transl_table>=23) {
    cerr << "Internal error infoOnAAatDNAPos(): transl_table " << static_cast<uint16>(transl_table) << " >= 23\n";
    exit(10);
  }

  CEBUG("dna: " << dna << endl);
  CEBUG("whichdnapos: " << whichdnapos << endl);
  CEBUG("orfstart   : " << orfstart << endl);
  CEBUG("orfdir     : " << static_cast<int16>(orfdir) << endl);
  CEBUG("codonstart : " << static_cast<uint16>(codonstart) << endl);

  codon.clear();
  codonvariants.clear();

  aminoacid.clear();
  isstart.clear();

  aanumber=0;
  posinaa=-1;

  char actcodon[3];
  if(orfdir>0){
    for(uint32 pos=orfstart+codonstart-1; pos<dna.size(); ){

      CEBUG("pos1a: " << pos << "   ");

      uint8 nuc=0;
      for(nuc=0; nuc<3 && pos<dna.size(); pos++){
	CEBUG("pos2: " << pos << "   ");
	if(dna[pos]!='*'){
	  actcodon[nuc]=dna[pos];
	  nuc++;
	}
	if(pos==whichdnapos) posinaa=nuc;
      }
      CEBUG(" // ");
      if(nuc == 3){
	aanumber++;
	if(posinaa >= 0) {
	  CEBUG("codon: " << static_cast<char>(actcodon[0]) << static_cast<char>(actcodon[1]) << static_cast<char>(actcodon[2]));
	  codon2AminoAcids(transl_table,
			   actcodon[0],
			   actcodon[1],
			   actcodon[2],
			   codonvariants,
			   aminoacid,
			   isstart);
	  codon+=actcodon[0];
	  codon+=actcodon[1];
	  codon+=actcodon[2];

	  CEBUG("Returning.\n");

	  return;
	}
	CEBUG("\n");
      }
    }
  }else{
    for(int32 pos=orfstart-codonstart+1; pos>=0; ){
      CEBUG("pos1b: " << pos << "   ");

      uint8 nuc=0;
      for(nuc=0; nuc<3 && pos>=0; pos--){
	CEBUG("pos2: " << pos << "   ");
	if(dna[pos]!='*'){
	  actcodon[nuc]=getComplementIUPACBase(dna[pos]);
	  nuc++;
	}
	if(pos==static_cast<int32>(whichdnapos)) posinaa=nuc;
      }
      CEBUG(" \\ ");
      if(nuc == 3) {
	aanumber++;
	if(posinaa >= 0){
	  CEBUG("codon: " << static_cast<char>(actcodon[0]) << static_cast<char>(actcodon[1]) << static_cast<char>(actcodon[2]));
	  codon2AminoAcids(transl_table,
			   actcodon[0],
			   actcodon[1],
			   actcodon[2],
			   codonvariants,
			   aminoacid,
			   isstart);
	  codon+=actcodon[0];
	  codon+=actcodon[1];
	  codon+=actcodon[2];

	  CEBUG("Returning.\n");

	  return;
	}
      }
      CEBUG("\n");
    }
  }
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * given two proteins that start at the same place, calculates a protein
 *  identity number
 *
 *************************************************************************/

double dptools::calcProteinIdentity(const std::string & p1, const std::string & p2)
{
  bool p1running=true;
  bool p2running=true;
  uint32 identicals=0;
  uint32 notidenticals=0;
  uint32 totalaa=static_cast<uint32>(std::max(p1.size(), p2.size()));

  for(uint32 i=0;i<totalaa;i++){
    if(p1running && p2running) {
      if(i<p1.size() && i<p2.size()) {
	if(toupper(p1[i]) == toupper(p2[i])){
	  identicals++;
	}
      }
    }else if(p1running||p2running){
      notidenticals++;
    }else{
      // oh, both proteins have ended (have a * for whatever reason in
      //  the sequence)
      // this stops the loop
      totalaa=i;
    }
    if(i>=p1.size()
       || (i<p1.size() && p1[i]=='*')) p1running=false;
    if(i>=p2.size()
       || (i<p2.size() && p2[i]=='*')) p2running=false;
  }

  if(totalaa==0) return 1.0;

  return (100.0*identicals/totalaa);
}





/*************************************************************************
 *
 *
 *
 *************************************************************************/

bool dptools::isCodonStart(const uint8 transl_table, const char base1, const char base2, const char base3)
{
  std::vector<char> aaresult;  // unused, but codon2AminoAcids() needs it
  std::vector<std::string> codonvariants;  // unused, but codon2AminoAcids() needs it

  std::vector<bool> isstartresult;

  if(transl_table==0 || transl_table>23) return false;
  switch(transl_table){
  case 7:
  case 8:
  case 17:
  case 18:
  case 19:
  case 20: return false;
  default: {
    // decrease transl_table by one to adapt to the c array
    codon2AminoAcids_wrapped(transl_table-1,
			     static_cast<char>(toupper(base1)),
			     static_cast<char>(toupper(base2)),
			     static_cast<char>(toupper(base3)),
			     codonvariants,
			     aaresult,
			     isstartresult);
    // not the slightest idea, why gcc 4.3.2 says:
    //  warning: conversion to ‘uint8’ from ‘int’ may alter its value [-Wconversion]

  }
  }

  if(aaresult.empty()) return false;

  return isstartresult[0];
}
