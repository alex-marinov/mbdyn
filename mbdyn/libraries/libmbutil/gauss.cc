/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <gauss.h>

/* Dati dei punti */
const doublereal dGaussPnt[] = {      
   0.,
   -.5773502691896, .5773502691896,
   -.7745966692415, 0., .7745966692415,
   -.8611363115941, -.3399810435849, .3399810435849, .8611363115941,
   -.90617985, -.53846931, 0., .53846931, .90617985
};

const doublereal dGaussWght[] = {   
   2.,
   1.,1.,
   .5555555555556,.8888888888889,.5555555555556,
   .3478548451375,.6521451548625,.6521451548625,.3478548451375,
   .23692689,.47862867,.568888888889,.47862867,.23692689
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
: NumIntData((iN < 1 ? 1 : (iN > 5 ? 5 : iN))), pdPnt(NULL), pdWght(NULL) {
   ASSERT(iNum >= 1 && iNum <= 5);
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


