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

/* Schur Data Manager */

/* 
 * Copyright 1999-2003 Giuseppe Quaranta <giuquaranta@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */
 
 /* E' una classe derivata da Data Manager che suddivide i dati in blocchi e li
   assegna ai vari processi.
  */
 
#ifndef SCHURDATAMAN_H
#define SCHURDATAMAN_H

#include <assert.h>
#include <mynewmem.h>


#include <dataman.h>
#include <mbpar.h>

#ifdef USE_MPI
#include <mpi++.h>
#endif /* USE_MPI */

class SchurDataManager : public DataManager {

 public:
  class ErrGeneric {}; 

 private:
  int iTotVertices;               /* # totale dei vertici: nodi = elementi */ 
  Elem** ppMyElems;                  /* Lista degli elementi assegnati a questo processo */
  
  int iNumLocElems;               /* # di elementi assegnati a questo processo */
  Elem** ppMyIntElems;
  int iNumIntElems;

  Node** ppMyNodes;                 /* Lista dei nodi assegnati a questo processo */
  int iNumLocNodes;               /* # di nodi assegnati a questo processo */
  integer* LocalDofs;
  int iNumLocDofs;
  integer* LocalIntDofs;
  int iNumIntDofs;
  Node** ppIntNodes;
  int iNumIntNodes;
  int iNumMyInt;
  integer* pMyIntDofs;
  unsigned int* pLabelsList;     /* struttura di servizio per CreatePartition*/
  int wgtflag;  
  int* pParAmgProcs; 
	
 protected:
#ifdef USE_MPI
  MPI::Intracomm DataComm; 
  MPI::Intracomm* pRotorComm;
#endif /* USE_MPI */
  int MyRank, DataCommSize;

  Node** ppExpCntNodes;
  Elem** ppExpCntElems;
  int iTotalExpConnections;
  
 private: 

  /* compatta il vettore delle adiacenze */
  void Pack(int* pList, int dim);

  /* Inizializza le varie liste che vengono utilizzate in CreatePartition */
  void InitList(int* list, int dim, int value);
  
  Node** SearchNode(Node** ppFirst, int dim, unsigned int& label);

  void OutputPartition(void);

 public:
  /* Costruttore Inizializza i dati */
  SchurDataManager(MBDynParser& HP, 
		  unsigned OF,
		  doublereal dInitialTime,
		  const char* sInputFileName,
		  const char* sOutputFileName,
		  flag fAbortAfterInput);

  /* Distruttore libera la memoria */
  ~SchurDataManager(void);


  /* Crea la partizione o la legge dal file fornito */
  void CreatePartition(void);


  /* Assembla il residuo */
  void AssRes(VectorHandler& ResHdl, doublereal dCoef);

  /* Assembla lo jacobiano */
  void AssJac(MatrixHandler& JacHdl, doublereal dCoef);


  /* Funzioni di aggiornamento dati durante la simulazione */
  void DerivativesUpdate(void) const;
  void BeforePredict(VectorHandler& X, VectorHandler& XP, VectorHandler& XPrev, VectorHandler& XPPrev) const;
  void AfterPredict(void) const;
  void Update(void)const;
  void AfterConvergence(void) const;

 /* stampa i risultati */
  void Output(void) const;

  enum DofType { TOTAL = 1, LOCAL = 2, INTERNAL = 3, MYINTERNAL = 4 };
  integer HowManyDofs(DofType who) {
    switch(who) {
    case TOTAL:
      return iTotDofs;
      
    case LOCAL:
      return iNumLocDofs;
      
    case INTERNAL:
      return iNumIntDofs;
      
    case MYINTERNAL:
      return iNumMyInt;
      
    default:
      std::cerr << "SchurDataManager::HowManyDofs: illegal request (" 
	      << unsigned(who) << ")" << std::endl;
      THROW(ErrGeneric());
    } 
  };

  integer* GetDofsList(DofType who) {
    switch(who) {
    case LOCAL:
      return LocalDofs;
      
    case INTERNAL:
      return LocalIntDofs;
      
    case MYINTERNAL: 
      return pMyIntDofs;

    default:
      std::cerr << "SchurDataManager::GetDofsList: illegal request ("
	      << unsigned(who) << ")" << std::endl;
      THROW(ErrGeneric());
    } 
  };
  
  Dof* pGetDofsList(void) {
    return pDofs;
  };
};
#endif /* SCHURDATAMAN_H */

