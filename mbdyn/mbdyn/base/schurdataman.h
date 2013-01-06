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

/* Schur Data Manager */

/* 
 * Copyright 1999-2013 Giuseppe Quaranta <quaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */
 
 /* E' una classe derivata da Data Manager che suddivide
  * i dati in blocchi e li assegna ai vari processi.
  */
 
#ifndef SCHURDATAMAN_H
#define SCHURDATAMAN_H

#include <assert.h>
#include <mynewmem.h>


#include <dataman.h>
#include <mbpar.h>

#ifdef USE_MPI
#include "ac/mpi.h"
#endif /* USE_MPI */
class Solver;

class SchurDataManager : public DataManager {
public:
	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	}; 

private:
	int iTotVertices;	/* # totale dei vertici: nodi = elementi */ 
	Elem** ppMyElems;	/* Lista elementi assegnati a questo processo */
	mutable VecIter<Elem *> MyElemIter;

	int iNumLocElems;	/* # di elementi assegnati a questo processo */
	Elem** ppMyIntElems;
	int iNumIntElems;

	Node** ppMyNodes;	/* Lista dei nodi assegnati a questo processo */
	int iNumLocNodes;	/* # di nodi assegnati a questo processo */
	integer* pLocalDofs;
	int iNumLocDofs;
	integer* pLocalIntDofs;
	int iNumIntDofs;
	Node** ppIntNodes;
	int iNumIntNodes;
	int iNumMyInt;
	integer* pMyIntDofs;
	unsigned int* pLabelsList;	/* str. di serv. per CreatePartition */
	enum {
		WEIGHT_NONE = 0x00U,
		WEIGHT_VERTICES = 0x01U,
		WEIGHT_COMM = 0x02U,
		WEIGHT_EDGES = 0x04U
	};
	unsigned wgtflag;
	int* pParAmgProcs; 

	enum PartitionLibrary {
		PARTITIONER_UNKNOWN,
		PARTITIONER_MANUAL,
		PARTITIONER_METIS,
		PARTITIONER_CHACO,

#if defined(USE_CHACO)
		PARTITIONER_DEFAULT = PARTITIONER_CHACO,
#elif defined(USE_METIS)
		PARTITIONER_DEFAULT = PARTITIONER_METIS,
#else /* ! USE_CHACO && ! USE_METIS */
		PARTITIONER_DEFAULT = PARTITIONER_MANUAL,
#endif /* ! USE_CHACO && ! USE_METIS */

		PARTITIONER_LAST
	} Partitioner;

protected:
#ifdef USE_MPI
	MPI::Intracomm DataComm; 
	MPI::Intracomm* pIndVelComm;
#endif /* USE_MPI */
	int MyRank, DataCommSize;

	Node** ppExpCntNodes;
	Elem** ppExpCntElems;
	int iTotalExpConnections;

private:

	/* compatta il vettore delle adiacenze */
	void Pack(int* pList, int dim);

	/* Inizializza le varie liste che vengono utilizzate 
	 * in CreatePartition */
	void InitList(int* list, int dim, int value);
	void InitList(float* list, int dim, int value);

	Node** SearchNode(Node** ppFirst, int dim, unsigned int& label);

	void OutputPartition(void);

public:
	/* Costruttore Inizializza i dati */
	SchurDataManager(MBDynParser& HP,
			unsigned OF,
			Solver* pS,
			doublereal dInitialTime,
			const char* sOutputFileName,
			const char* sInputFileName,
			bool bAbortAfterInput);

	/* Distruttore libera la memoria */
	~SchurDataManager(void);

	/* Crea la partizione o la legge dal file fornito */
	void CreatePartition(void);

	/* Assembla il residuo */
	void AssRes(VectorHandler& ResHdl, doublereal dCoef) 
		throw(ChangedEquationStructure);

	/* Assembla lo jacobiano */
	void AssJac(MatrixHandler& JacHdl, doublereal dCoef);

	/* Funzioni di aggiornamento dati durante la simulazione */
	void DerivativesUpdate(void) const;
	void BeforePredict(VectorHandler& X, VectorHandler& XP, 
			VectorHandler& XPrev, VectorHandler& XPPrev) const;
	void AfterPredict(void) const;
	void Update(void) const;
	void AfterConvergence(void) const;

	/* stampa i risultati */
	void Output(bool force = false) const;

	enum DofType { TOTAL = 1, LOCAL = 2, INTERNAL = 3, MYINTERNAL = 4 };

	integer HowManyDofs(DofType who) const;
  	integer* GetDofsList(DofType who) const;

	Dof* pGetDofsList(void) const;
};

#endif /* SCHURDATAMAN_H */

