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

#include "scf_look.H"
#include "globals.H"


/* --------------------- */
/*       afh_Info        */
/* --------------------- */

/*
afh_Info::afh_Info(SCF_look *r)
{
  read   = r;
  status = fhs_UNHANDLED;
  next   = nullptr;

  original_base   = 'X';
  hypothesis_base = 'X';
}
*/


afh_Info::afh_Info(SCF_look *r, const int32 dbPos, const int32 scfPos,
		   const char hypoBase, const char oldBase, const int32 fClass)
{
  theTrace = nullptr;
  contig_read = nullptr;

  status   = fhs_UNHANDLED;
  confirmed= -1.0;
  next     = nullptr;
  question = fClass;

  dbPosition       = dbPos;
  scfPosition      = scfPos;
  hypothesis_base  = hypoBase;
  original_base    = oldBase;
  original_quality = 99;
}


afh_Info::afh_Info(const Contig::contigread_t *r, const int32 dbPos,
		   const int32 scfPos,  const char hypoBase,
		   const char oldBase, const int32 fClass)
{
  theTrace    = nullptr;
  contig_read = r;

  status   = fhs_UNHANDLED;
  confirmed= -1.0;
  next     = nullptr;
  question = fClass;

  dbPosition      = dbPos;
  scfPosition     = scfPos;
  hypothesis_base = hypoBase;
  original_base   = oldBase;
  original_quality = 99;
}



afh_Info::~afh_Info()
{
  if (theTrace != nullptr) {
    //    SCFBuffer::bufferDelete(theTrace);
    delete theTrace;
  }
  theTrace = nullptr;
}


/*
afh_Info::afh_Info(afh_Info const &other)
{
  clear();
  *this=other;
}
*/


bool afh_Info::compareWith(afh_Info &a)
{
  return (contig_read == a.getContigRead());
}



SCF_look* afh_Info::performEditOnRead(afh_Info &a)
{
  SCF_look *unedited;

  if (theTrace != nullptr) {
    delete theTrace;
  }
  theTrace = nullptr;

  unedited = a.getTraceData();

  if (unedited == nullptr) {
    // we can not perform any edits because there is no trace.
    // e.g. if the scf-file is missing.
    return nullptr;
  }

  theTrace  = new SCF_look;
  *theTrace = *unedited;

  int32 operation = a.getQuestion();
  char  old_base  = a.getOriginalBase();
  char  new_base  = a.getHypothesisBase();


#ifdef VERBOSE
  cout << "VORHER" << endl;
  theTrace->showDBArea(a.dbPos(), 5, cout);
  cout << "\n" << old_base << " => " << new_base << endl;
#endif



  if (operation  == fhc_SYNTACTIC) {
    if (old_base  == '.' || old_base == ' ') {

#ifdef VERBOSE
      cout << "fhc_SYNTACTIC -> fhc_MISSING" << endl;
#endif
      operation = fhc_MISSING;
    }
    if (new_base == '.' || new_base == ' ') {
#ifdef VERBOSE
      cout << "fhc_SYNTACTIC -> fhc_ADDITIONAL" << endl;
#endif

      //      assert(old_base == '*');
      //theTrace->deleteBaseDB(a.dbPos());
      operation = fhc_ADDITIONAL;
    }
  }

  switch (operation) {
  case fhc_OVERCALL:
  case fhc_ADDITIONAL:

#ifdef VERBOSE
    cout << "Delete Base" << endl;
    cout << *this;
#endif

    if (new_base != '*') {
      theTrace->deleteBaseDB(a.dbPos());
    } else {
      theTrace->alterBaseDB(a.dbPos(), new_base);
    }
    break;

  case fhc_UNDERCALL:
  case fhc_MISSING:
#ifdef VERBOSE
    cout << "Insert Base" << endl;
    cout << *this;
#endif

    if (old_base != '*') {
      theTrace->insertBaseDB(a.dbPos(), a.scfPos(), new_base);
    } else {
      if (theTrace->getDBBase(a.dbPos()) == '*') {
	theTrace->alterBaseDB(a.dbPos(), new_base);
      } else {
	if (a.dbPos() > 0 && theTrace->getDBBase(a.dbPos())-1 == '*') {
	  theTrace->alterBaseDB(a.dbPos()-1, new_base);
	} else {
	  theTrace->alterBaseDB(a.dbPos()+1, new_base);
	}
      }
    }
    break;

  case fhc_WRONG:
#ifdef VERBOSE
    cout << "Alter Base" << endl;
    cout << *this;
#endif
    theTrace->alterBaseDB(a.dbPos(), new_base);
    break;

#ifdef VERBOSE
  default:
    cout << "Syntactic" << endl;
#endif

  }

#ifdef VERBOSE
  cout << "\nNACHHER:" << endl;
  theTrace->showDBArea(a.dbPos(), 5, cout);
  cout << endl << endl;
#endif

  return theTrace;
}



const char *afh_Info::getReadName() const
{
  if (theTrace != nullptr) {
    return theTrace->getFileName();
  }
  if (contig_read != nullptr) {
    return (contig_read->read).getName().c_str();
  }
  return nullptr;
}



const char* afh_Info::getFaultClassName() const
{
 static const char *fhClassName[] = {
  "Undefined",
  "Overcall",
  "Additional Call",
  "Wrong Call",
  "Undercall",
  "Missing Call",
  "Syntactic",
  "Correct Call",
  "Correct Gap"
  };

 if (question < 0 || question > 8) {
   return fhClassName[0];
 } else {
   return fhClassName[question];
 }
}



bool afh_Info::getIsExtendedAlterOperation()
{
  switch (question) {
  case fhc_WRONG: return true;
  case fhc_OVERCALL:
  case fhc_ADDITIONAL:
    return hypothesis_base == '*';
  case fhc_UNDERCALL:
  case fhc_MISSING:
    return original_base == '*';
  }
  return false;
}




bool afh_Info::loadRead()
{
  if (confirmed < 0) {
    if (theTrace == nullptr) {
      try {
	//	theTrace = new SCF_look;
	//*theTrace = *(ScfBuffer::bufferRead(contig_read->read,
	//				   contig_read->direction));

	theTrace = &(ScfBuffer::bufferReadCopy(contig_read->read,
					       contig_read->direction));
      }
      catch (Notify n) {
	theTrace  = nullptr;
	reject(0);
      }
    }
  }
  return (theTrace != nullptr);
}



float afh_Info::eval(const bool verbose)
{
  if (theTrace == nullptr) {
    loadRead();
  }

  if (theTrace != nullptr) {

    if (verbose) {
      confirmed =  evaluate(scfPosition, dbPosition, question,
			    hypothesis_base, theTrace, cout);
    } else {
      confirmed =  evaluate(scfPosition, dbPosition, question,
			    hypothesis_base, theTrace);
    }

  }
  return confirmed;
}



const char* afh_Info::getFaultStatusName() const
{
  static const char *fhStatusName[] = {
    "UNDEFINED",
    "UNHANDLED",
    "CONFIRMED",
    "REJECTED",
    "UNDECIDABLE"
  };

  if (status < 0 || status > 4) {
    return fhStatusName[0];
  } else {
    return fhStatusName[status];
  }
}


bool afh_Info::isConfirmed() const
{
  return (status == fhs_CONFIRMED);
}


void afh_Info::confirm(const float confirmValue) {
  status  = fhs_CONFIRMED;
  confirmed = confirmValue;
}


void afh_Info::reject(const float confirmValue) {
  status  = fhs_REJECTED;
  confirmed = confirmValue;
}



afh_Info* afh_Info::nextAfh() {
  return next;
}



int16 afh_Info::isSupportingHypotheses()
{
  return (question == fhc_GAP || question == fhc_CORRECT);
}



bool afh_Info::isNPlusHypotheses()
{
  return (isUndefinedBase(original_base) && hypothesis_base != '*');
}


bool afh_Info::isNToBaseHypotheses()
{
  return (isUndefinedBase(original_base) && isRealBase(hypothesis_base));
}



int16 afh_Info::isReversed() {
  if (theTrace != nullptr) {
    return theTrace->isReversed();
  }
  if (contig_read != nullptr) {
    return (contig_read->direction < 0);
  }
  return 0;
}



int16 afh_Info::appendAfh(afh_Info *a) {
  afh_Info *xx;

  if (a == nullptr) { return 1; }

  xx = this;
  while (xx->next != nullptr) { xx = xx->next; }

  xx->next = a;
  a->next = nullptr;

  return 0;
}


void afh_Info::setDBPos(int32 newPos)
{
  dbPosition = newPos;
}

void afh_Info::setSCFPos(int32 newPos)
{
  scfPosition = newPos;
}


// verschiebe alle editpositionen hinter hypothese wenn sie den selben
// read betreffen und die position des edits nach pos liegt

void afh_Info::shiftLeftAfterPos(int32 startPos)
{
  afh_Info *l = this;

  cout << "shift left" << endl;
  while (l != nullptr) {
    if (l->compareWith(*this))  {
      if (l->dbPos() >= startPos) {
	l->setDBPos(l->dbPos() - 1);

	/*
	if (l->isReversed()) {
	  l->setSCFPos(l->scfPos() + 1);
	} else {
	  l->setSCFPos(l->scfPos() - 1);
	}
	*/
      }
    }
    l = l->next;
  }
}



void afh_Info::shiftRightAfterPos(int32 startPos)
{
  afh_Info *l = this;

  cout << "shift right" << endl;
  while (l != nullptr) {
    if (l->compareWith(*this))  {
      if (l->dbPos() >= startPos) {
	l->setDBPos(l->dbPos() + 1);
	/*
	if (l->isReversed()) {
	  l->setSCFPos(l->scfPos() - 1);
	} else {
	  l->setSCFPos(l->scfPos() + 1);
	}
	*/
      }
    }
    l = l->next;
  }
}



char afh_Info::dbBase()
{
  if (theTrace == nullptr) {
    throw Notify(Notify::FATAL, "dbBase", "SCF-File not loaded!!");
  }
  return theTrace->getDBBase(dbPosition);
}



// -------------------------------------------------------


ostream &operator<<(ostream &ostr, afh_Info const &i)
{
  ostr << "DB: " << i.dbPosition;
  ostr << "\tSCF: " << i.scfPosition;
  ostr << "\t" << i.getReadName();
  ostr << "\t" << i.getFaultClassName();
  ostr << "\t" << i.getFaultStatusName();
  ostr << "\t" << i.original_base << "=>" << i.hypothesis_base;
  ostr << "\tconf " << i.confirmed
       << " BQ " << ((int)(i.original_quality))
       << " ID" << i.scfInfoId << "\t" << i.theTrace << endl;

  if (i.next != nullptr) { ostr << *(i.next); }

  return ostr;
}



// Is there something in the chemistry that supports that the
// correct bases are before-base-after?
// For reversed reads take the complement of the DB-bases.
int16 afh_Info::chemistryEvaluation(char before, char base, char after)
{
  int result = 100;

  if (isReversed()) {
    if ((toupper(after) == 'A') && (toupper(base) == 'G')) {
      result = 200;
      DEBUG_EDIT("Supported by chemistry (reversed)\n");
    }
  } else {
    if ((toupper(before) == 'A') && (toupper(base) == 'G')) {
      DEBUG_EDIT("Supported by chemistry\n");
      result = 200;
    }
  }

  return result;
}
