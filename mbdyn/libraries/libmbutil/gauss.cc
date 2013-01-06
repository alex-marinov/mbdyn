/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* FIXME: better use double precision numbers ... */

/* Dati dei punti */
const doublereal dGaussPnt[] = {      
   /*  1 point,  2nd order */  0.0000000000000,
   /*  2 point,  4th order */ -0.5773502691896,  0.5773502691896,
   /*  3 point,  6th order */ -0.7745966692415,  0.0000000000000,  0.7745966692415,
   /*  4 point,  8th order */ -0.8611363115941, -0.3399810435849,  0.3399810435849,  0.8611363115941,
   /*  5 point, 10th order */ -0.90617985, -0.53846931,  0.00000000,  0.53846931,  0.90617985,
   /*  6 point, 12th order */ -0.93246951, -0.66120939, -0.23861919,  0.23861919,  0.66120939,  0.93246951,
   /*  7 point, 14th order */ -0.94910791, -0.74153119, -0.40584515,  0.00000000,  0.40584515,  0.74153119,  0.94910791,
   /*  8 point, 16th order */ -0.96028986, -0.79666648, -0.52553241, -0.18343464,  0.18343464,  0.52553241,  0.79666648,  0.96028986,
   /*  9 point, 18th order */ -0.96816024, -0.83603111, -0.61337143, -0.32425342,  0.00000000,  0.32425342,  0.61337143,  0.83603111,  0.96816024,
   /* 10 point, 20th order */ -0.97390653, -0.86506337, -0.67940957, -0.43339539, -0.14887434,  0.14887434,  0.43339539,  0.67940957,  0.86506337,  0.97390653, 
   /* 11 point, 22nd order */ 
};

const doublereal dGaussWght[] = {   
   /*  1 point,  2nd order */  2.0000000000000,
   /*  2 point,  4nd order */  1.0000000000000,  1.0000000000000,
   /*  3 point,  6nd order */  0.5555555555556,  0.8888888888889,  0.5555555555556,
   /*  4 point,  8nd order */  0.3478548451375,  0.6521451548625,  0.6521451548625,  0.3478548451375,
   /*  5 point, 10nd order */  0.23692689,  0.47862867,  0.56888889,  0.47862867,  0.23692689,
   /*  6 point, 12th order */  0.17132449,  0.36076157,  0.46791393,  0.46791393,  0.36076157,  0.17132449,
   /*  7 point, 14th order */  0.12948497,  0.27970539,  0.38183005,  0.41795918,  0.38183005,  0.27970539,  0.12948497,
   /*  8 point, 16th order */  0.10122854,  0.22238103,  0.31370665,  0.36268378,  0.36268378,  0.31370665,  0.22238103,  0.10122854,
   /*  9 point, 18th order */  0.08127439,  0.18064816,  0.26061070,  0.31234708,  0.33023936,  0.31234708,  0.26061070,  0.18064816,  0.08127439,
   /* 10 point, 20th order */  0.06667134,  0.14945135,  0.21908636,  0.26926672,  0.29552422,  0.29552422,  0.26926672,  0.21908636,  0.14945135,  0.06667134,
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


