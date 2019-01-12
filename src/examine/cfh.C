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
 *
 * SCF examination methods
 *
 * Written by Thomas Pfisterer
 *
 *
 */


#include "examine/scf_look.H"


/* --------------------- */
/*       cfh_Info        */
/* --------------------- */

void cfh_Info::foolCompiler(void)
{
#include "stdinc/foolcompiler.C"
}



cfh_Info::cfh_Info()
{
  first = nullptr;
  last  = nullptr;

  rejected_faults = 0;
  rejected_n = 0;
  confirmed_positive = 0;
  confirmed_reverse = 0;
  positive = 0;
  reversed = 0;
}


cfh_Info::cfh_Info(afh_Info *a)
{
  first = a;
  last  = a;

  rejected_faults = 0;
  rejected_n = 0;
  confirmed_positive = 0;
  confirmed_reverse = 0;
  positive = 0;
  reversed = 0;
}


cfh_Info::~cfh_Info()
{
  deleteAllBuffered();
}



void cfh_Info::appendAfh(afh_Info *a)
{
  if (last == nullptr) {
    first = a;
    last  = a;
  } else {
    last->appendAfh(a);
    last = a;
  }
}


void cfh_Info::deleteAllBuffered()
{
  afh_Info *weg;

  while (first != nullptr) {
    weg = first;
    first = first->nextAfh();

    delete weg;
  }

  first = nullptr;
  last  = nullptr;
}



void cfh_Info::countInformation()
{
  afh_Info *loop = first;

  confirmed_reverse = 0;
  confirmed_positive = 0;
  reversed = 0;
  positive = 0;
  rejected_faults = 0;
  rejected_n = 0;

  while (loop != nullptr) {
    // DEBUG_EDIT("Status: " << loop->getStatus()) ;
    // DEBUG_EDIT(" Direction: " << loop->isReversed() << " " << endl);

    if (loop->isSupportingHypotheses()){
      if (loop->getStatus() == fhs_CONFIRMED) {
	if (loop->isReversed()) {
	  confirmed_reverse++;
	} else {
	  confirmed_positive++;
	}
      }
      if (loop->isReversed()) {
	reversed++;
      } else {
	positive++;
      }
    } else {
	// good, fault-hypotheses
	if (loop->getStatus() != fhs_CONFIRMED) {
	  if (false == loop->isNToBaseHypotheses()) {
	    rejected_faults++;
	  } else {
	    //rejected_faults++;
	    rejected_n++;
	  }
	} else {
	  if (loop->isReversed()) {
	    confirmed_reverse++;
	  } else {
	    confirmed_positive++;
	  }
	}
    }
    loop = loop->nextAfh();
  }
}



bool cfh_Info::needsInformationOfType(
               const int16 question, const int16 direction)
{
  countInformation();
#ifdef HYPOTHESES_VERBOSE
  cout << *this;
#endif

  // Cannot become true;

  if (rejected_faults != 0) return false;

  if (question == fhc_GAP || question == fhc_CORRECT) {
    // Supporting Hypotheses....
    if (confirmed_positive > 0 && confirmed_reverse > 0) {
      return false;
    }
    if (confirmed_positive > 0) return (direction < 0);
    if (confirmed_reverse  > 0) return (direction > 0);
    return true;
  } else {
    // Fault hypotheses: always
    return true;
  }
}



// ****************************************************
// examineCfhPartial:
// examines if a hypotheses is confirmed and if
// the hypotheses can become true. Returns false
// if there is no possibiliy to confirm the Cfh.
// ****************************************************

bool cfh_Info::examineCfhPartial(bool &isConfirmed)
{
  countInformation();

#ifdef HYPOTHESES_VERBOSE
  cout << *this;
#endif


  //  isConfirmed = (rejected_faults == 0 && confirmed_positive > 0 &&
  //		 confirmed_reverse > 0);
  isConfirmed = rejected_faults == 0;

  return (rejected_faults == 0);

}



ostream & operator<<(ostream &ostr, cfh_Info const &i)
{
  ostr << "===== cfh_info ===== " << endl;
  ostr << "rejected_faults   : " << i.rejected_faults;
  ostr << "   rejected_n   " << i.rejected_n;
  ostr << "   conf+   " << i.confirmed_positive
       << "  (" << i.positive << ")";
  ostr << "   conf- : " << i.confirmed_reverse
       << "  (" << i.reversed << ")";
  ostr << "\tScore: " << i.rel_score << endl;

  if (i.first != nullptr) {
    ostr << *(i.first);
  }

  ostr << "==================== " << endl;
  return ostr;
}



// Evaluate the afh's and confirm if they are above the given threshold
// if complete == false evaluate until a hypotheses is rejected
// otherwise evaluate all of them.

bool cfh_Info::eval(const bool complete, const bool verbose,
		    const float threshold, bool strict_N)
{

  afh_Info* l = first;
  bool  result = true;
  float f;
  float score = 0;
  int32 count = 0;


#ifndef RUNONLY
  if (verbose==true) {
    cout << "Evaluate: " << endl;
    cout << *this << endl;
  }
#endif

  while (l != nullptr && (complete || result)) {
    if (loadAfhTraces(l)) {

      f = l->eval(verbose);

      if (f > threshold) {
	l->confirm();

	if (l->getFaultClass() != fhc_SYNTACTIC) {
	  score = score + 1.0;
	  count++;
	}
      } else {
	if (l->isNPlusHypotheses() && strict_N == false) {
          l->confirm();
        }

	if (l->isNToBaseHypotheses() == false || strict_N == true) {
          score = score + f;
          count++;
          result = false;
        }
      }
    } else {
      result = false;
      l->reject();
    }
    l = l->nextAfh();
  }


  if (l == nullptr && count > 0) {
    rel_score = score / count;
  } else {
    rel_score = -1.0;
  }

  return result;
}


//
// confirm only a single class of hypothesis
// above a given threshold.

bool cfh_Info::singleEval(const float threshold, int32 fhc_class)
{
  afh_Info* l = first;
  bool hasSinglePositive = false;
  float f;

  while (l != nullptr) {
    f = l->eval(false);
    if (l->getQuestion() == fhc_class && f > threshold) {
      l->confirm();
      hasSinglePositive = true;
    } else {
      l->reject();
    }
    l = l->nextAfh();
  }

  return hasSinglePositive;
}


//
// confirm only a single class of hypothesis
// above a given threshold.

bool cfh_Info::evalExtendedAlterOperations(const float threshold)
{
  afh_Info* l = first;
  bool hasSinglePositive = false;
  float f;

  while (l != nullptr) {
    f = l->eval(false);

    if (l->getIsExtendedAlterOperation() && f > threshold) {
      l->confirm();
      hasSinglePositive = true;
    } else {
      l->reject();
    }
    l = l->nextAfh();
  }

  return hasSinglePositive;
}


bool cfh_Info::confirmQuality(const int32 threshold)
{
  afh_Info* l = first;

  while (l != nullptr) {
    if (l->getOriginalBaseQuality() != 0 &&
        l->getOriginalBaseQuality() < threshold) {
      l->confirm(1.0);
    }
    l = l->nextAfh();
  }

  return true;
}


//   If more of the afh in a cfh edit the same base, there is a problem:
//   edit operations should be performed on the SCF_look before the next
//   operation is examined - but we would like to avoid to copy the object
//   if not absolutely necessary.
//   Thus we test for an afh if another operation affected the same read.
//
//


bool cfh_Info::loadAfhTraces(afh_Info *a)
{
  afh_Info* lastEdit = readEditedBefore(a);

  if (lastEdit == nullptr) {
    // cout << "New Trace" << endl;
    return a->loadRead();
  } else {
    // cout << "Copy Trace" << endl;
    a->performEditOnRead(*lastEdit);
  }

  return true;
}



afh_Info* cfh_Info::readEditedBefore(afh_Info *a)
{
  afh_Info *l = first;
  afh_Info *found = nullptr;


  while (l != a) {
    if (l->compareWith(*a))  {
      found = l;
    }
    l = l->nextAfh();
  }

  return found;

}
