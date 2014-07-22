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

/* inertia element */

#ifndef INERTIA_H
#define INERTIA_H

#include "elem.h"
#include "strnode.h"
#include "gravity.h"

/* CenterOfMass - begin */

class CenterOfMass {
protected:
	std::set<const ElemGravityOwner *> elements;

	mutable doublereal dMass;
	mutable Vec3 S;
	mutable Mat3x3 J;

	/* Center of mass kinematics and inertia */
	Vec3 X_cm;
	Vec3 V_cm;
	Vec3 Omega_cm;
	Mat3x3 J_cm;

	virtual void Collect_int(void);

	virtual std::ostream& Output_int(std::ostream& out) const;

public:
	/* Costruttore definitivo (da mettere a punto) */
	CenterOfMass(std::set<const ElemGravityOwner *>& elements);
	virtual ~CenterOfMass(void);
};

/* CenterOfMass - end */

/* Inertia - begin */

class Inertia :
virtual public Elem, public ElemGravityOwner, public InitialAssemblyElem, public CenterOfMass {
protected:
	/* reference point */
	unsigned flags;
	Vec3 X0;
	Mat3x3 R0;
	Mat3x3 J0;

	Mat3x3 R_princ;
	Vec3 J_princ;

	/* momento statico */
	Vec3 GetS_int(void) const;

	/* momento d'inerzia */
	Mat3x3 GetJ_int(void) const;

	virtual void Collect_int(void);

	virtual std::ostream& Output_int(std::ostream& out) const;

public:
	/* Costruttore definitivo (da mettere a punto) */
	Inertia(unsigned int uL, const std::string& sN, std::set<const ElemGravityOwner *>& elements,
		const Vec3& x0, const Mat3x3& r0, std::ostream& log, flag fOut);

	virtual ~Inertia(void);

	/* massa totale */
	doublereal dGetM(void) const;

	/* Tipo dell'elemento (usato solo per debug ecc.) */
	virtual Elem::Type GetElemType(void) const;

	/* Numero gdl durante l'assemblaggio iniziale */
	virtual unsigned int iGetInitialNumDof(void) const;

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void Output(OutputHandler& OH) const;   

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Dimensione del workspace durante l'assemblaggio iniziale.
	 * Occorre tener conto del numero di dof che l'elemento definisce
	 * in questa fase e dei dof dei nodi che vengono utilizzati.
	 * Sono considerati dof indipendenti la posizione e la velocita'
	 * dei nodi */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

	/* Usata per inizializzare la quantita' di moto */
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

    virtual unsigned int iGetNumPrivData(void) const;
    virtual unsigned int iGetPrivDataIdx(const char *s) const;
    virtual doublereal dGetPrivData(unsigned int i) const;
};

/* Inertia - end */

#endif // INERTIA_H

