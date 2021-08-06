/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

/* GENEL filtro scalare analogico */

#ifndef GENFILT_H
#define GENFILT_H

#include "genel.h"
#include "scalarvalue.h"

/* GenelStateSpaceSISO - begin */

class GenelStateSpaceSISO : public Genel {
protected:
	ScalarDof SD_y;		/* uscita */
	ScalarValue *SV_u;	/* ingresso */

	unsigned int iNumDofs;

	doublereal* pdE;
	doublereal* pdA;
	doublereal* pdB;
	doublereal* pdC;
	doublereal dD;

	doublereal* pdX;
	doublereal* pdXP;

public:
	GenelStateSpaceSISO(unsigned int uLabel, const DofOwner* pDO,
		const ScalarDof& y, ScalarValue* u,
		unsigned int Order,
		doublereal *pE,
		doublereal* pA, doublereal* pB,
		doublereal* pC, doublereal D,
		bool bBalance,
		doublereal *pdX0,
		doublereal *pdXP0,
		flag fOutput);

	virtual ~GenelStateSpaceSISO(void);

	virtual unsigned int iGetNumDof(void) const;

	/* esegue operazioni sui dof di proprieta' dell'elemento */
	virtual DofOrder::Order GetDofType(unsigned int i) const;

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo di Genel */
	virtual Genel::Type GetGenelType(void) const {
		return Genel::STATESPACESISO;
	};

	/* Dimensioni del workspace */
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ );

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	/* output; si assume che ogni tipo di elemento sappia, attraverso
	 * l'OutputHandler, dove scrivere il proprio output */
	virtual void Output(OutputHandler& OH) const;

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* GenelStateSpaceSISO - end */


/* GenelStateSpaceMIMO - begin */

class GenelStateSpaceMIMO : public Genel {
protected:
	unsigned int iNumOutputs;
	unsigned int iNumInputs;
	ScalarDof* pvSD_y; /* uscite */
	std::vector<ScalarValue *> SV_u; /* ingressi */

	unsigned int iNumDofs;

	doublereal* pdE;
	doublereal* pdA;
	doublereal* pdB;
	doublereal* pdC;
	doublereal* pdD;

	doublereal* pdX;
	doublereal* pdXP;

public:
	GenelStateSpaceMIMO(unsigned int uLabel, const DofOwner* pDO,
		unsigned int iNumOut, const ScalarDof* y,
		std::vector<ScalarValue *>& u,
		unsigned int Order,
		doublereal* pE,
		doublereal* pA, doublereal* pB,
		doublereal* pC, doublereal* pD,
		bool bBalance,
		doublereal *pdX0,
		doublereal *pdXP0,
		flag fOutput);

	virtual ~GenelStateSpaceMIMO(void);

	virtual unsigned int iGetNumDof(void) const;

	/* esegue operazioni sui dof di proprieta' dell'elemento */
	virtual DofOrder::Order GetDofType(unsigned int i) const;

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo di Genel */
	virtual Genel::Type GetGenelType(void) const {
		return Genel::STATESPACEMIMO;
	};

	/* Dimensioni del workspace */
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ );

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	/* output; si assume che ogni tipo di elemento sappia, attraverso
	 * l'OutputHandler, dove scrivere il proprio output */
	virtual void Output(OutputHandler& OH) const;

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* GenelStateSpaceMIMO - end */

#endif // GENFILT_H
