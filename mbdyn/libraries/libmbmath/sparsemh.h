/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

/*
 * sparse matrix managing
 */

#ifndef SPARSEMH_H
#define SPARSEMH_H

#define SPARSE_MATRIX_NO_RESET

#include <myassert.h>
#include <mynewmem.h>
#include <except.h>
#include <spdata.h>
#include <solman.h>
#include <submat.h>

/* classi dichiarate in questo file */
class SparseData;           /* gestore di sparsita' */
class SparseMatrixHandler;  /* gestore matrice sparsa (assemblaggio) */

/*
 * Union usata per impaccare il vettore di indici di riga e colonna
 * nel gestore di sparsita' (i vettori sono integer = long int, mentre
 * in questo modo gli indici massimi sono 65535)
 */

/* SparseData - begin */

/*
 * Gestore di entita' sparse con hash + linked list.
 * usa spazio di lavoro messo a disposizione da altri.
 */

class SparseData {
public:
   	class ErrNoRoom {};

	union uPacVec {
		integer iInt;
		struct {
			unsigned short int ir;
			unsigned short int ic;
		} sRC;
	};

private:

protected:
   	integer iMaxSize;   	/* Dimensione della memoria allocata */
   	integer iCurSize;   	/* Dimensione corrente */
   	integer iNumItem;   	/* Numero di entries correnti */
   	integer** ppiTable; 	/* Vettore della linked list */
   	integer** ppiKeys;  	/* Vettore delle keys */

	integer iv[3];		/* internal data for fortran routines :( */

public:

   	/*
    	 * Costruttore.
	 * iSize -   dimensione degli spazi di lavoro
	 * piTable - puntatore a punt. ad un array di interi di dim. iSize,
	 *           che viene usato come tabella di stato del campo,
	 * piKeys  - puntatore a punt. ad un array delle stesse dimensioni del
	 *           precedente, che viene usato per contenere le keys
	 */
   	SparseData(integer iSize, integer** ppiTmpTable, integer** ppiTmpKeys);

   	/* Distruttore (banalissimo: non c'e' nulla da distruggere) */
   	~SparseData(void);

	/* Inizializza i vettori con il valore -(iSize+1) ecc. */
	void ResetVec(void);

	/* Risetta la dimensione corrente */
	bool SetCurSize(integer i);

   	/*
	 * Trova la posizione del termine dato da iKey
	 * Ritorna un termine che rappresenta la posizione
	 * della key di ingresso nel vettore piKey;
	 * - se il termine e' positivo, la key esiste ed e' valida,
	 *   oppure l'inserzione e' riuscita correttamente;
	 * - se il termine e' negativo, la key e' stata definita
	 * e poi cancellata (attualmente non supportato)
	 * - se il termine e' zero l'inserzione non e' riuscita
	 *   (ad esempio perche' si e' esaurito lo spazio)
	 *
	 * FIXME: the call to kd01b creates room for the iKey index,
	 * 	  thus the matrix gets filled up if it is completely
	 *	  traversed.
	 */
   	 inline integer iGetIndex(integer iKey) {
      		ASSERT(ppiTable != NULL);
      		ASSERT(ppiKeys != NULL);
      		ASSERT(*ppiTable != NULL);
      		ASSERT(*ppiKeys != NULL);
      		ASSERT(iKey > 0);
#ifdef DEBUG_MEMMANAGER
      		ASSERT(defaultMemoryManager.fIsPointerToBlock(*ppiTable));
      		ASSERT(defaultMemoryManager.fIsPointerToBlock(*ppiKeys));
#endif /* DEBUG_MEMMANAGER */

      		integer iFree = 0;
      		__FC_DECL__(kd01b)(iv, *ppiTable, *ppiKeys, &iKey, &iFree);
      		if (iFree == 0) {
			std::cerr << "SparseData: there's no room left "
				"in matrix" << std::endl;
	 		THROW(SparseData::ErrNoRoom());
      		}

      		return iFree;
   	};

   	inline integer iGetIndexIfAvailable(integer iKey) {
      		ASSERT(ppiTable != NULL);
      		ASSERT(ppiKeys != NULL);
      		ASSERT(*ppiTable != NULL);
      		ASSERT(*ppiKeys != NULL);
      		ASSERT(iKey > 0);
#ifdef DEBUG_MEMMANAGER
      		ASSERT(defaultMemoryManager.fIsPointerToBlock(*ppiTable));
      		ASSERT(defaultMemoryManager.fIsPointerToBlock(*ppiKeys));
#endif /* DEBUG_MEMMANAGER */

      		integer iFree = 0;
      		__FC_DECL__(kd01c)(iv, *ppiTable, *ppiKeys, &iKey, &iFree);
      		if (iFree == 0) {
			std::cerr << "SparseData: there's no room left "
				"in matrix" << std::endl;
	 		THROW(SparseData::ErrNoRoom());
      		}

      		return iFree;
   	};
};

/* SparseData - end */


/* SparseMatrixHandler - begin */

/*
 * Gestore di matrici sparse; usa spazio messo a disposizione da altri;
 * usa il gestore di sparsita' SparseData
 */

class SparseMatrixHandler : public MatrixHandler {
private:
   	integer iWorkSpaceSize; /* dimensione del workspace */
   	integer iCurSize; /* dimensione corrente del workspace */
   	integer iNumItem; /* numero di elementi effettivamente contenuti */

   	integer iMatSize; /* Ordine della matrice (supposta quadrata) */

   	SparseData* pSD;  /* Puntatore all'oggetto SparseData da utilizzare */

protected:
   	integer** ppiRow; /* punt. a punt. ad array di int. di dim. iMaxSize */
   	integer** ppiCol; /* '' '' */
   	doublereal** ppdMat; /* p. a p. ad array di reali di dim. iMaxSize */
   	doublereal* pdMatm1; /* p. a p. ad array di reali di dim. iMaxSize */

public:

   	/*
	 * Costruttore
	 * iMSize   - ordine della matrice;
	 * ppiRow   - puntatore ad array di interi di dimensione iSize,
	 *            contiene la tabella usata da SparseData;
	 * ppiCol   - puntatore ad array di interi di dimensione iSize,
	 *            contiene le keys usate da SparseData
	 * ppdMat   - puntatore ad array di reali, contiene la matrice;
	 * iWSSize  - dimensione del workspace
	 */
   	SparseMatrixHandler(integer iMSize,
		       	    integer** ppiTmpRow,
			    integer** ppiTmpCol,
			    doublereal** ppdTmpMat,
			    integer iWSSize);

   	/* Distruttore banale - non c'e' nulla da distruggere */
   	~SparseMatrixHandler(void);

   	/* Usato per il debug */
   	virtual void IsValid(void) const;

   	/*
	 * Inizializzazione della matrice: la riempe con il valore desiderato
	 * ed inizializza il gestore di sparsita'
	 */
   	void Init(const doublereal& dResetVal = 0.);

	/* Risetta la dimensione corrente */
	integer iGetCurSize(void) const;

	/* Risetta la dimensione corrente */
	bool SetCurSize(integer i);

	/* Compatta la matrice */
	integer iPacVec(void);

   	/* Restituisce un puntatore all'array di reali della matrice */
   	inline doublereal* pdGetMat(void) const {
		return *ppdMat;
	};

   	/* Restituisce un puntatore al vettore delle righe */
   	inline integer* piGetRows(void) const {
		return *ppiRow;
	};

   	/* Restituisce un puntatore al vettore delle colonne */
   	inline integer* piGetCols(void) const {
		return *ppiCol;
	};

   	/*
	 * Impacchetta la matrice; restituisce il numero di elementi
	 * diversi da zero
	 */
   	integer PacMat(void);

   	/* Inserisce un coefficiente */
   	inline flag fPutCoef(integer iRow, integer iCol,
			     const doublereal& dCoef);

   	/* Incrementa un coefficiente - se non esiste lo crea */
   	inline flag fIncCoef(integer iRow, integer iCol,
			     const doublereal& dCoef);

   	/* Incrementa un coefficiente - se non esiste lo crea */
   	inline flag fDecCoef(integer iRow, integer iCol,
			     const doublereal& dCoef);

   	/* Restituisce un coefficiente - zero se non e' definito */
   	inline const doublereal& dGetCoef(integer iRow, integer iCol) const;

	virtual const doublereal&
	operator () (integer iRow, integer iCol) const = 0;

	virtual doublereal&
	operator () (integer iRow, integer iCol) = 0;

   	/* dimensioni */
   	virtual integer iGetNumRows(void) const {
      		return iMatSize;
   	};

   	virtual integer iGetNumCols(void) const {
      		return iMatSize;
   	};

	std::ostream& SparseOutput(std::ostream& out,
			int w = 20) const;
};



inline flag
SparseMatrixHandler::fPutCoef(integer iRow,
		      	      integer iCol,
			      const doublereal& dCoef)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

   	ASSERT((iRow > 0) && (iRow <= iMatSize));
   	ASSERT((iCol > 0) && (iCol <= iMatSize));

   	if (dCoef != doublereal(0.)) {
      		SparseData::uPacVec uPV;
      		uPV.sRC.ir = (unsigned short int)(iRow);
      		uPV.sRC.ic = (unsigned short int)(iCol);

      		integer iField = uPV.iInt;
      		integer iReturnFlag = abs(pSD->iGetIndex(iField));
      		pdMatm1[iReturnFlag] = dCoef;

      		return flag(0);
   	}

   	return flag(1);
}

inline flag
SparseMatrixHandler::fIncCoef(integer iRow,
			      integer iCol,
			      const doublereal& dCoef)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

   	ASSERT((iRow > 0) && (iRow <= iMatSize));
   	ASSERT((iCol > 0) && (iCol <= iMatSize));

   	if (dCoef != doublereal(0.)) {
      		SparseData::uPacVec uPV;
      		uPV.sRC.ir = (unsigned short int)(iRow);
      		uPV.sRC.ic = (unsigned short int)(iCol);

      		integer iField = uPV.iInt;
      		integer iReturnFlag = pSD->iGetIndex(iField);
#ifndef SPARSE_MATRIX_NO_RESET
      		iReturnFlag = abs(iReturnFlag);
      		pdMatm1[iReturnFlag] += dCoef;
#else /* SPARSE_MATRIX_NO_RESET */
		if (iReturnFlag < 0) {
			/*
			 * already in: add
			 */
			pdMatm1[-iReturnFlag] += dCoef;

		} else {
			/*
			 * first insert: set (so there's no need
			 * to reset the matrix, another 1% speedup)
			 */
			pdMatm1[iReturnFlag] = dCoef;
		}
#endif /* SPARSE_MATRIX_NO_RESET */

      		return flag(0);
   	}

   	return flag(1);
}

inline flag
SparseMatrixHandler::fDecCoef(integer iRow,
		      	      integer iCol,
			      const doublereal& dCoef)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

   	ASSERT((iRow > 0) && (iRow <= iMatSize));
   	ASSERT((iCol > 0) && (iCol <= iMatSize));

   	if (dCoef != doublereal(0.)) {
      		SparseData::uPacVec uPV;
      		uPV.sRC.ir = (unsigned short int)(iRow);
      		uPV.sRC.ic = (unsigned short int)(iCol);

      		integer iField = uPV.iInt;
      		integer iReturnFlag = pSD->iGetIndex(iField);
#ifndef SPARSE_MATRIX_NO_RESET
      		iReturnFlag = abs(iReturnFlag);
      		pdMatm1[iReturnFlag] -= dCoef;
#else /* SPARSE_MATRIX_NO_RESET */
		if (iReturnFlag < 0) {
			/*
			 * already in: sub
			 */
			pdMatm1[-iReturnFlag] -= dCoef;

		} else {
			/*
			 * first insert: set (so there's no need
			 * to reset the matrix, another 1% speedup)
			 */
			pdMatm1[iReturnFlag] = dCoef;
		}
#endif /* SPARSE_MATRIX_NO_RESET */

      		return flag(0);
   	}

   	return flag(1);
}

inline const doublereal&
SparseMatrixHandler::dGetCoef(integer iRow, integer iCol) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

   	ASSERT((iRow > 0) && (iRow <= iMatSize));
   	ASSERT((iCol > 0) && (iCol <= iMatSize));

   	SparseData::uPacVec uPV;
   	uPV.sRC.ir = short(iRow);
   	uPV.sRC.ic = short(iCol);

	/*
	 * warning: the key is allocated if it does not exist,
	 * so the matrix gets filled if traversed
	 */
   	integer iField = uPV.iInt;
   	integer iReturnFlag = pSD->iGetIndexIfAvailable(iField);
	/* pSD->iGetIndex throws an exception if no room left in matrix */
	if (iReturnFlag < 0) {
      		return pdMatm1[-iReturnFlag];
	}

   	return ::dZero; /* zero! */
}

inline doublereal&
SparseMatrixHandler::operator () (integer iRow, integer iCol)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

   	ASSERT((iRow > 0) && (iRow <= iMatSize));
   	ASSERT((iCol > 0) && (iCol <= iMatSize));

   	SparseData::uPacVec uPV;
   	uPV.sRC.ir = short(iRow);
   	uPV.sRC.ic = short(iCol);

	/*
	 * warning: the key is allocated if it does not exist,
	 * so the matrix gets filled if traversed
	 */
   	integer iField = uPV.iInt;
   	integer iReturnFlag = abs(pSD->iGetIndex(iField));
      	return pdMatm1[iReturnFlag];
}

inline const doublereal&
SparseMatrixHandler::operator () (integer iRow, integer iCol) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

   	ASSERT((iRow > 0) && (iRow <= iMatSize));
   	ASSERT((iCol > 0) && (iCol <= iMatSize));

   	SparseData::uPacVec uPV;
   	uPV.sRC.ir = short(iRow);
   	uPV.sRC.ic = short(iCol);

	/*
	 * warning: the key is allocated if it does not exist,
	 * so the matrix gets filled if traversed
	 */
   	integer iField = uPV.iInt;
   	integer iReturnFlag = pSD->iGetIndex(iField);

   	if (iReturnFlag < 0) {
      		return pdMatm1[-iReturnFlag];

	}

   	return ::dZero; /* zero! */
}

/* SparseMatrixHandler - end */

#endif /* SPARSEMH_H */

