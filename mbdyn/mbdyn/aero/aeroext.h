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
 
 /* Aerodynmic External Element by Giuseppe Quaranta(C) -  December 2002
  * <quaranta@aero.polimi.it>	
  */

#ifndef AEROEXT_H
#define AEROEXT_H

#ifdef USE_EXTERNAL
#include "aerodyn.h"
#include "ac/mpi.h"

#include "myassert.h"
#include "except.h"

#include "modal.h"
/* AerodynamicExternal  */


class AerodynamicExternal :
	virtual public Elem,
	public AerodynamicElem,
	public DriveOwner
{
private:
	MyVectorHandler*   pdBuffer;    /* buffer per lo scambio dei dati riguardanti 
				           le posizioni e le forze */
	MyVectorHandler*   pdBufferVel; /* buffer per lo scambio dei dati riguardanti le velocita' */
	integer            NodeN;
	const StructNode** ppNode;      /* vettore di puntatori a i nodi che compongono il blocco */
	const doublereal*  pRefLength;  /* lunghezza di riferimento per scalare i punti definiti nella matrice offset 
					   associati a ciascun nodo*/
					 
	int                OffN;	/* numero di punto rigidamente collegati con coascun nodo */
	Mat3xN*            pOffsetVectors; /* posizioni dei punto collegati con ciascun nodo */
	MPI::Intercomm*    pInterfComm;    /* Intercomunicatore con il codice di interfacia */
	MPI::Request* 	   pRecReq;
	bool               VelFlag;
	bool		   MomFlag;
	bool               SentFlag;
	Mat3xN*            pForce;    	/* output forze  */
	Mat3xN*            pMoms;    	/* output momenti  */
	unsigned int*      pLabList;    /* lista label nodi */

public:
	
	AerodynamicExternal(unsigned int uLabel,
		const DofOwner *pDO,
		const DriveCaller* pDC,
		int NN,
		const StructNode** ppN,
		const doublereal* RefL,
		MPI::Intercomm* IC,
		flag fOut,
		bool VF,
		bool MF);

	AerodynamicExternal(unsigned int uLabel,
		const DofOwner *pDO,
		const DriveCaller* pDC, 
		int NN,
		const StructNode** ppN,
		const doublereal* RefL,
		MPI::Intercomm* IC,
		int	ON,
		Mat3xN* OV,
		flag fOut,
		bool VF,
		bool MF);
		
	virtual ~AerodynamicExternal(void);
	
	bool NeedsAirProperties(void) const
	{ return true; };
	
	/* Tipo dell'elemento (usato per debug ecc.) */
   	Elem::Type GetElemType(void) const 
	{ return Elem::AERODYNAMIC;};  
	
	   /* Contributo al file di restart */
   	virtual std::ostream& Restart(std::ostream& out) const 
	{ return out << std::endl; };
	
	virtual void BeforePredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ ,
		VectorHandler& /* XPrev */ ,
		VectorHandler& /* XPPrev */ ) const { NO_OP; };
		
		
	virtual void AfterPredict(VectorHandler& X  , 
		VectorHandler&  XP  );
		
	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 0;
		*piNumCols = 0;
	};
	
	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
	       doublereal /* dCoef */ ,
	       const VectorHandler& /* XCurr */ ,
	       const VectorHandler& /* XPrimeCurr */ ) {
	       	DEBUGCOUTFNAME("AerodynamicExternal::AssJac");
		WorkMat.SetNullMatrix();
		return WorkMat;
	};
	
	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
	
	virtual void Update(const VectorHandler&  XCurr  , 
		const VectorHandler& XPrimeCurr );
	
	/*
	 * Elaborazione stato interno dopo la convergenza
	 */
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP) { NO_OP; };

	/*
	 * output; si assume che ogni tipo di elemento sappia, attraverso
	 * l'OutputHandler, dove scrivere il proprio output
	 */
	virtual void Output(OutputHandler& OH) const;
 
 	/* Tipo di elemento aerodinamico */
	virtual AerodynamicElem::Type GetAerodynamicElemType(void) const {
		return AerodynamicElem::AERODYNAMICEXTERNAL;
	};
	
	/* *******PER IL SOLUTORE PARALLELO******** */        
	/*
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
	virtual int GetNumConnectedNodes(void) const {
		return NodeN;
	};
	
	
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(NodeN);
		for (int i = 0; i < NodeN; i++) {
			connectedNodes[i] = ppNode[i];
		}
	};

private:
	void ConstructAndInitialize(void);
	
	void Send(const VectorHandler& X  , 
		const VectorHandler&  XP  );

};




/* AerodynamicExternalModal */
class AerodynamicExternalModal :
	virtual public Elem,
	public AerodynamicElem,
	public DriveOwner
{
protected:
	
	MyVectorHandler*   pdBuffer;    /* buffer per lo scambio dei dati riguardanti 
				           le posizioni e le forze */
	MyVectorHandler*   pdBufferVel; /* buffer per lo scambio dei dati riguardanti le velocita' */
	
	Modal* 	   	   pModal;
	int 		   ModalNodes;
	MPI::Intercomm*    pInterfComm;  /* Intercomunicatore con il codice di interfaccia */
	MPI::Prequest*     pSenReq;
	MPI::Prequest* 	   pRecReq;
	bool               VelFlag;
	bool		   MomFlag;
	bool               SentFlag;
	MyVectorHandler*   pForce;    /* buffer per l'output delle forze */

public:
	
	AerodynamicExternalModal(unsigned int uLabel,
		const DofOwner* pDO,
		const DriveCaller* pDC,
		Modal* pM,
		MPI::Intercomm* IC,
		flag fOut,
		bool VelFlag,
		bool MomFlag);

	virtual ~AerodynamicExternalModal(void);
	
	bool NeedsAirProperties(void) const
	{ return false; };
	
	/* Tipo dell'elemento (usato per debug ecc.) */
   	Elem::Type GetElemType(void) const 
	{ return Elem::AERODYNAMIC;};  
	
	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 0;
		*piNumCols = 0;
	};
	
	   /* Contributo al file di restart */
   	virtual std::ostream& Restart(std::ostream& out) const 
	{ return out << std::endl; };

	virtual void BeforePredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ ,
		VectorHandler& /* XPrev */ ,
		VectorHandler& /* XPPrev */ ) const { NO_OP; };
		
		
	virtual void AfterPredict(VectorHandler& X  , 
		VectorHandler&  XP  );
		
		
	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
	       doublereal /* dCoef */ ,
	       const VectorHandler& /* XCurr */ ,
	       const VectorHandler& /* XPrimeCurr */ ) {
	       	DEBUGCOUTFNAME("AerodynamicExternalModal::AssJac");
		WorkMat.SetNullMatrix();
		return WorkMat;
	};
	
	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
	       
	virtual void Update(const VectorHandler&  XCurr  , 
		const VectorHandler& XPrimeCurr );
	/*
	 * Elaborazione stato interno dopo la convergenza
	 */ 
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP)
			{ NO_OP; };

	/*
	 * output; si assume che ogni tipo di elemento sappia, attraverso
	 * l'OutputHandler, dove scrivere il proprio output
	 */
	virtual void Output(OutputHandler& OH) const;
 
 	/* Tipo di elemento aerodinamico */
	virtual AerodynamicElem::Type GetAerodynamicElemType(void) const {
		return AerodynamicElem::AERODYNAMICEXTERNAL;
	};
	
	/* *******PER IL SOLUTORE PARALLELO******** */        
	/*
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
	virtual int GetNumConnectedNodes(void) const {
		return 1;
	};
	
	
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(1);
		connectedNodes[0] = pModal->pGetModalNode();
	};
	
private:
	void Send(const VectorHandler& X, const VectorHandler& XP);
};

class DataManager;
class MBDynParser;

extern Elem *
ReadAerodynamicExternal(DataManager* pDM, MBDynParser& HP,
	const DofOwner *pDO, unsigned int uLabel);

extern Elem *
ReadAerodynamicExternalModal(DataManager* pDM, MBDynParser& HP,
	const DofOwner *pDO, unsigned int uLabel);

#endif /* USE_EXTERNAL */
#endif /* AEROEXT_H */
