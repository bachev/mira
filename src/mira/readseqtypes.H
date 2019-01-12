
// These enums live in ReadGroupLib
// should be included ONLY from ReadGroupLib.H

enum read_sequencing_type : unsigned int {
  SEQTYPE_SANGER=0,
    SEQTYPE_454GS20,
    SEQTYPE_IONTORRENT,
    SEQTYPE_PACBIOHQ,
    SEQTYPE_PACBIOLQ,
    SEQTYPE_TEXT,
    SEQTYPE_SOLEXA,
    SEQTYPE_ABISOLID,
    SEQTYPE_END
};

// cludge to get the C(++) preprocessor know the number of sequencing types, adapt
//  according to the enum array above

#ifndef CPP_READ_SEQTYPE_END
#define CPP_READ_SEQTYPE_END 8
#endif


enum read_naming_scheme : unsigned int {
  SCHEME_UNKNOWN=0,
    SCHEME_SANGER,
    SCHEME_TIGR,
    SCHEME_FR,          // .f, .r
    SCHEME_SOLEXA,      // /1, /2
    SCHEME_STLOUIS,
    SCHEME_SRARAW,
    SCHEME_EMPTY,      // nothing will be set
    SCHEME_NONE
};


// have unknown = 0, FF > 0 and FR/RF <0
enum splace_code {
  SPLACE_UNUSED_LOW=-3,
  SPLACE_RF,
  SPLACE_FR,
  SPLACE_UNKNOWN,
  SPLACE_SF,
  SPLACE_SB,
  SPLACE_SU,
  SPLACE_UNUSED_HIGH
};
