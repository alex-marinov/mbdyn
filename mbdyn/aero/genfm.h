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

#ifndef GENFM_H
#define GENFM_H

/* Generic force/moment element based on table lookup */
/* originated from GARTEUR HC AG-16 activity to model generic aerodynamic
 * forces related to rigid body motion as a function of dynamic pressure
 * angle of attack and sideslip angle, with data provided by DLR for Bo105
 * fuselage and tail empennages */

#include "aerodyn.h"

struct GenericAerodynamicData {
	std::string name;

	bool bAlphaFirst;

	/* angle of attack and sideslip angle datapoints number */
	int nAlpha;
	int nBeta;

	std::vector<doublereal> Alpha;
	std::vector<doublereal> Beta;

	struct GenericAerodynamicCoef {
		doublereal	dCoef[6];

		GenericAerodynamicCoef(void);
		GenericAerodynamicCoef(const GenericAerodynamicCoef& c);
		GenericAerodynamicCoef operator + (const GenericAerodynamicCoef& c) const;
		GenericAerodynamicCoef operator - (const GenericAerodynamicCoef& c) const;
		GenericAerodynamicCoef operator * (const doublereal& d) const;
		GenericAerodynamicCoef operator / (const doublereal& d) const;
	};

	std::vector<std::vector<GenericAerodynamicCoef> > Data;
};

class GenericAerodynamicForce :
	virtual public Elem,
	public AerodynamicElem,
	public InitialAssemblyElem
{
protected:
	/* Node the forces are applied to */
	const StructNode* pNode;
	/* Reference surface and length used to dimensionalize
	 * non-dimensional coefficients */
	const doublereal dRefSurface;
	const doublereal dRefLength;
	const bool bAlphaFirst;
	/* Offset of aerodynamic center with respect to node position */
	const Vec3 tilde_f;
	/* orientation of aerodynamic reference frame with respect to node */
	const Mat3x3 tilde_Ra;	/* Rotaz. del sistema aerodinamico al nodo */

	/* force and moment */
	Vec3 tilde_F;
	Vec3 tilde_M;
	Vec3 F;
	Vec3 M;

	// persistent
	doublereal dAlpha, dBeta;

	/* aerodynamic data */
	GenericAerodynamicData *pData;

	/* Assemblaggio residuo */
	void AssVec(SubVectorHandler& WorkVec);

public:
	GenericAerodynamicForce(unsigned int uLabel,
		const DofOwner *pDO,
		const StructNode* pN,
		const Vec3& fTmp, const Mat3x3& RaTmp,
		doublereal dS, doublereal dL, bool bAlphaFirst,
		GenericAerodynamicData *pData,
		flag fOut);
	virtual ~GenericAerodynamicForce(void);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo dell'elemento (usato per debug ecc.) */
	virtual Elem::Type GetElemType(void) const;

	/* funzioni proprie */

	/* Dimensioni del workspace */
	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);

	/*
	 * output; si assume che ogni tipo di elemento sappia, attraverso
	 * l'OutputHandler, dove scrivere il proprio output
	 */
	virtual void Output(OutputHandler& OH) const;

	/* Dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;

	/* Numero di GDL iniziali */
	virtual unsigned int iGetInitialNumDof(void) const;

	/* Dimensioni del workspace */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);

	/* Tipo di elemento aerodinamico */
	virtual AerodynamicElem::Type GetAerodynamicElemType(void) const;

	/*
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;

	/* describes the dimension of components of equation */
    virtual std::ostream& DescribeEq(std::ostream& out,
		  const char *prefix = "",
		  bool bInitial = false) const;
};

extern Elem *
ReadGenericAerodynamicForce(DataManager* pDM, MBDynParser& HP,
	const DofOwner *pDO, unsigned int uLabel);

/* GenericAerodynamicForce - end */

#endif /* GENFM_H */

