/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* Punti di Gauss */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "gauss.h"


/* Dati dei punti */
const doublereal dGaussPnt[] = {      
   /*  1 point,  2nd order */  0.000000000000000,
   /*  2 point,  4th order */ -0.577350269189626,  0.577350269189626,
   /*  3 point,  6th order */ -0.774596669241483,  0.000000000000000,  0.774596669241483,
   /*  4 point,  8th order */ -0.861136311594053, -0.339981043584856,  0.339981043584856,  0.861136311594053,
   /*  5 point, 10th order */ -0.906179845938664, -0.538469310105683,  0.000000000000000,  0.538469310105683,  0.906179845938664,
   /*  6 point, 12th order */ -0.932469514203152, -0.661209386466264, -0.238619186083197,  0.238619186083197,  0.661209386466264,  0.932469514203152,
   /*  7 point, 14th order */ -0.949107912342758, -0.741531185599394, -0.405845151377397,  0.000000000000000,  0.405845151377397,  0.741531185599394,  0.949107912342758, 
   /*  8 point, 16th order */ -0.960289856497536, -0.796666477413627, -0.525532409916329, -0.183434642495650,  0.183434642495650,  0.525532409916329,  0.796666477413627,  0.960289856497536,
   /*  9 point, 18th order */ -0.968160239507626, -0.836031107326636, -0.613371432700590, -0.324253423403809,  0.000000000000000,  0.324253423403809,  0.613371432700591,  0.836031107326636,  0.968160239507626,
   /* 10 point, 20th order */ -0.973906528517171, -0.865063366688985, -0.679409568299025, -0.433395394129247, -0.148874338981631,  0.148874338981631,  0.433395394129247,  0.679409568299024,  0.865063366688984,  0.973906528517172
   /* 11 point, 22nd order */ 
};

const doublereal dGaussWght[] = {   
   /*  1 point,  2nd order */  2.000000000000000,
   /*  2 point,  4nd order */  1.000000000000000,  1.000000000000000,
   /*  3 point,  6nd order */  0.555555555555556,  0.888888888888889,  0.555555555555556,
   /*  4 point,  8nd order */  0.347854845137454,  0.652145154862546,  0.652145154862546,  0.347854845137454,
   /*  5 point, 10nd order */  0.236926885056189,  0.478628670499367,  0.568888888888889,  0.478628670499367,  0.236926885056189,
   /*  6 point, 12th order */  0.171324492379171,  0.360761573048139,  0.467913934572690,  0.360761573048139,  0.171324492379171, 
   /*  7 point, 14th order */  0.129484966168870,  0.279705391489277,  0.381830050505119,  0.417959183673469,  0.381830050505119,  0.279705391489277,  0.129484966168870,
   /*  8 point, 16th order */  0.101228536290376,  0.222381034453375,  0.313706645877887,  0.362683783378362,  0.362683783378362,  0.313706645877887,  0.222381034453375,  0.101228536290376, 
   /*  9 point, 18th order */  0.081274388361575,  0.180648160694858,  0.260610696402935,  0.312347077040002,  0.330239355001259,  0.312347077040002,  0.260610696402936,  0.180648160694857,  0.081274388361575,
   /* 10 point, 20th order */  0.066671344308688,  0.149451349150581,  0.219086362515981,  0.269266719309996,  0.295524224714752,  0.295524224714752,  0.269266719309996,  0.219086362515982,  0.149451349150581,  0.066671344308688
   /* 11 point, 22nd order */ 
};

const doublereal dTrapezoidPnt[] = {
   0.,
   -1., 1.,
   -1., 0., 1.,
   -1., -.3333333333333333, .3333333333333333, 1.,
   -1., -.5, 0., .5, 1.
};

const doublereal dTrapezoidWght[] = {
   2.,
   1., 1.,
   .3333333333333333, 1.333333333333333, .3333333333333333,
   .25, .75, .75, .25,
   .1555555555555556, .7111111111111111, .2666666666666667,
       .7111111111111111, .1555555555555556
};

/* GaussData - begin */

GaussData::GaussData(integer iN)
: NumIntData((iN < 1 ? 1 : (iN > 10 ? 10 : iN))), pdPnt(NULL), pdWght(NULL) {
   ASSERT(iNum >= 1 && iNum <= 10);
   integer i = iNum*(iNum-1)/2-1;
   (doublereal*&)pdPnt = (doublereal*)dGaussPnt+i;
   (doublereal*&)pdWght = (doublereal*)dGaussWght+i;
}


doublereal GaussData::dGetPnt(integer i) const
{
   ASSERT(i > 0 && i <= iNum);
   return pdPnt[i];
}

   
doublereal GaussData::dGetWght(integer i) const
{   
   ASSERT(i > 0 && i <= iNum);
   return pdWght[i];
}


PntWght GaussData::Get(integer i) const
{
   ASSERT(i > 0 && i <= iNum);
   return PntWght(pdPnt[i], pdWght[i]);
}


const doublereal* GaussData::pdGetPnt(void) const
{
   return pdPnt+1;
}


const doublereal* GaussData::pdGetWght(void) const
{
   return pdWght+1;
}

/* GaussData - end */


/* TrapezoidData - begin */

TrapezoidData::TrapezoidData(integer iN)
: NumIntData((iN < 1 ? 1 : (iN > 5 ? 5 : iN))), pdPnt(NULL), pdWght(NULL) {
   integer i = iNum*(iNum-1)/2-1;
   (doublereal*&)pdPnt = (doublereal*)dTrapezoidPnt+i;
   (doublereal*&)pdWght = (doublereal*)dTrapezoidWght+i;
}


doublereal TrapezoidData::dGetPnt(integer i) const
{
   ASSERT(i > 0 && i <= iNum);
   return pdPnt[i];
}

   
doublereal TrapezoidData::dGetWght(integer i) const
{   
   ASSERT(i > 0 && i <= iNum);
   return pdWght[i];
}


PntWght TrapezoidData::Get(integer i) const
{
   ASSERT(i > 0 && i <= iNum);
   return PntWght(pdPnt[i], pdWght[i]);
}


const doublereal* TrapezoidData::pdGetPnt(void) const
{
   return pdPnt+1;
}


const doublereal* TrapezoidData::pdGetWght(void) const
{
   return pdWght+1;
}

/* TrapezoidData - end */


/* GaussDataIterator - begin */

GaussDataIterator::GaussDataIterator(integer iN)
: GaussData(iN), iCurr(1)
{
   NO_OP;
}


doublereal GaussDataIterator::dGetFirst(integer i) const
{
   ASSERT(i == 0 || i == 1);
   
   (integer&)iCurr = 1;
   if(i == 0) {
      return pdPnt[1];
   }
   /* if(i == 1) */
   return pdWght[1];   
}


PntWght GaussDataIterator::GetFirst(void) const
{
   (integer&)iCurr = 1;
   return PntWght(pdPnt[1], pdWght[1]);
}


flag GaussDataIterator::fGetNext(doublereal& d, integer i) const
{
   ASSERT(i == 0 || i == 1);
   
   (integer&)iCurr += 1;
   if(iCurr > iNum) {
      return flag(0);
   }
   
   if(i == 0) {
      d = pdPnt[iCurr];
   } else /* if(i == 1) */ {
      d = pdWght[iCurr];
   }   
   
   return flag(1);
}


doublereal GaussDataIterator::dGetCurrPnt(void) const
{
   ASSERT(iCurr >= 1 && iCurr <= iNum); 

   return pdPnt[iCurr];
}


doublereal GaussDataIterator::dGetCurrWght(void) const
{   
   ASSERT(iCurr >= 1 && iCurr <= iNum); 
   
   return pdWght[iCurr];
}


flag GaussDataIterator::fGetNext(PntWght& PW) const
{
   (integer&)iCurr += 1;
   if(iCurr > iNum) {
      return flag(0);
   }
   
   PW = PntWght(pdPnt[iCurr], pdWght[iCurr]);
   
   return flag(1);
}

/* GaussdataIterator - end */


/* NumIntIterator - begin */

NumIntIterator::NumIntIterator(NumIntData& d)
: iCurr(1), data(d)
{
   NO_OP;
}


NumIntIterator::~NumIntIterator(void)
{
	NO_OP;
}


doublereal NumIntIterator::dGetFirst(integer i) const
{
   ASSERT(i == 0 || i == 1);
   
   (integer&)iCurr = 1;
   if(i == 0) {
      return data.dGetPnt(1);
   }
   /* else if(i == 1) */
   return data.dGetWght(1);   
}


PntWght NumIntIterator::GetFirst(void) const
{
   (integer&)iCurr = 1;
   return PntWght(data.Get(1));
}


flag NumIntIterator::fGetNext(doublereal& d, integer i) const
{
   ASSERT(i == 0 || i == 1);
   
   (integer&)iCurr += 1;
   if(iCurr > data.iGetNum()) {
      return flag(0);
   }
   
   if(i == 0) {
      d = data.dGetPnt(iCurr);
   } else /* if(i == 1) */ {
      d = data.dGetWght(iCurr);
   }   
   
   return flag(1);
}


doublereal NumIntIterator::dGetCurrPnt(void) const
{
   ASSERT(iCurr >= 1 && iCurr <= data.iGetNum()); 

   return data.dGetPnt(iCurr);
}


doublereal NumIntIterator::dGetCurrWght(void) const
{   
   ASSERT(iCurr >= 1 && iCurr <= data.iGetNum()); 
   
   return data.dGetWght(iCurr);
}


flag NumIntIterator::fGetNext(PntWght& PW) const
{
   (integer&)iCurr += 1;
   if(iCurr > data.iGetNum()) {
      return flag(0);
   }
   
   PW = PntWght(data.Get(iCurr));
   
   return flag(1);
}

/* NumIntIterator - end */


