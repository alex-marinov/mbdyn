/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* Dati e metodi relativi ai punti di integrazione di Gauss */

#ifndef GAUSS_H
#define GAUSS_H

#include <ac/f2c.h>

#include <myassert.h>

/* Struttura punto di Gauss.
 * Ha due campi: la coordinata ed il peso.
 * Possiede inoltre due metodi per l'accesso in lettura ed un costruttore.
 * Se si desidera una maggiore robustezza, la si puo' dichiarare come
 * classe. In tale caso l'accesso in scrittura e' assicurato grazie al
 * costruttore (e quindi e' riservato a chi genera il punto), mentre la 
 * lettura avviene esclusivamente attraverso i metodi.
 */

struct PntWght {
   doublereal dPnt;
   doublereal dWght;

   PntWght(doublereal dP, doublereal dW) : dPnt(dP), dWght(dW) { NO_OP; };
   doublereal dGetPnt(void) const { return dPnt; };
   doublereal dGetWght(void) const { return dWght; };
};


/* Versione classe
class PntWght {   
 private:
   doublereal dPnt;
   doublereal dWght;
   
 public:
   PntWght(doublereal dP, doublereal dW) : dPnt(dP), dWght(dW) { NO_OP; };
   doublereal dGetPnt(void) const { return dPnt; };
   doublereal dGetWght(void) const { return dWght; };
};
 */


/* Classe set di punti di Gauss.
 * Viene costruita a partire dal numero di punti che si desidera.
 * Contiene metodi per la lettura dei singoli punti e pesi, come reali
 * o come strutture di punti.
 * Inoltre consente l'accesso ai vettori dei punti (inutile e sconsigliato)
 * 
 * I dati sono limitati. Sono disponibili per set da 1 a 5 punti.
 * Se il costruttore viene chiamato con un numero minore di 1, 
 * vengono restituiti i dati relativi ad 1 punto; analogamente, 
 * se il costruttore viene chiamato con un numero maggiore di 5,
 * vengono restituiti i dati relativi a 5 punti.
 * Se si utilizza l'iteratore descritto piu' avanti, a parte la perdita
 * di precisione che ne consegue, l'integrazione non risente quindi
 * di un valore al di fuori del range previsto.
 */

class NumIntData {
 protected:
   integer iNum;
 public:
   NumIntData(const integer& i) : iNum(i) {};
   virtual ~NumIntData(void) {};
   virtual integer iGetNum(void) const { return iNum; };
   virtual doublereal dGetPnt(integer i) const = 0;
   virtual doublereal dGetWght(integer i) const = 0;
   virtual PntWght Get(integer i) const = 0;
};

class GaussData : public NumIntData {
 protected:
   const doublereal* const pdPnt;
   const doublereal* const pdWght;
   
 public:   
   GaussData(integer iN);
   doublereal dGetPnt(integer i) const;
   doublereal dGetWght(integer i) const;
   PntWght Get(integer i) const;   
   const doublereal* pdGetPnt(void) const;   
   const doublereal* pdGetWght(void) const;
};

class TrapezoidData : public NumIntData {
 protected:
   const doublereal* const pdPnt;
   const doublereal* const pdWght;
   
 public:   
   TrapezoidData(integer iN);
   doublereal dGetPnt(integer i) const;
   doublereal dGetWght(integer i) const;
   PntWght Get(integer i) const;   
   const doublereal* pdGetPnt(void) const;   
   const doublereal* pdGetWght(void) const;
};


/* Classe iteratore sui punti di Gauss.
 * Dal momento che in ogni caso la classe set di punti di Gauss in ogni caso
 * e' un'interfaccia per l'accesso diretto ai dati, l'iteratore e' derivato
 * dalla classe GaussData in contrasto con la filosofia degli iteratori.
 * In questo modo lo stesso oggetto consente l'accesso diretto ed iterativo
 * ai dati.
 * L'utilizzo con la struttura di dati del singolo punto e' consigliato
 * perche' consente costrutti del tipo:
 * 
 *     GaussDataIterator gdi(n);
 *     PntWght p = gdi.GetFirst();
 *     do {
 *        //... use p.dGetPnt(), p.dGetWght() ...
 *     } while(gdi.fGetNext(p));
 * 
 * che risultano molto essenziali ed eleganti. E' possibile un uso analogo,
 * ma meno elegante, dei singoli punti o pesi, qualora fossero richiesti
 * solo gli uni o gli altri:
 * 
 *     GaussDataIterator gdi(n);
 *     doublereal dp = gdi.dGetFirst();
 *     do {
 *        doublereal dw = gdi.dGetCurrWght();
 *        //...
 *     } while(gdi.fGetNext(dp));
 * 
 * ATTENZIONE: se si usa dGetCurrPnt() o dGetCurrWght() dopo che fGetNext()
 * ha ritornato FALSE, il risultato e' unpredictable in quanto il contatore
 * e' volutamente out of range.
 */



class GaussDataIterator : public GaussData {
 protected:
   integer iCurr;
 public:
   GaussDataIterator(integer iN);
   doublereal dGetFirst(integer i = 0) const;
   PntWght GetFirst(void) const;
   flag fGetNext(doublereal& d, integer i = 0) const;
   doublereal dGetCurrPnt(void) const;
   doublereal dGetCurrWght(void) const;
   flag fGetNext(PntWght& PW) const;
};


class NumIntIterator {
 protected:
   integer iCurr;
   NumIntData& data;
 public:
   NumIntIterator(NumIntData& d);
   virtual ~NumIntIterator(void);
   virtual doublereal dGetFirst(integer i = 0) const;
   virtual PntWght GetFirst(void) const;
   virtual flag fGetNext(doublereal& d, integer i = 0) const;
   virtual doublereal dGetCurrPnt(void) const;
   virtual doublereal dGetCurrWght(void) const;
   virtual flag fGetNext(PntWght& PW) const;
};

#endif /* GAUSS_H */
