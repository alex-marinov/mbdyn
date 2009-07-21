/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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

/* Cyclocopter inflow models */

#ifndef CYCLOCOPTER_H
#define CYCLOCOPTER_H

#include "indvel.h"

/* CyclocopterKARI - begin */

/*

From:

"A New VTOL UAV Cyclocopter with Cycloidal Blades System",

Chul Yong Yun, Illkyung Park,
Ho Yong Lee, Jai Sang Jung, In Seong Hwang,
Seung Jo Kim,
Sung Nam Jung

Presented at the American Helicopter Society 60th Annual Forum,
Baltimore, MD, June 7-10, 2004.
*/

class CyclocopterKARI
: virtual public Elem, public InducedVelocity {
protected:
	const StructNode* pRotor;
	Mat3x3 RRot;

public:
	CyclocopterKARI(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, ResForceSet **ppres, flag fOut);
	virtual ~CyclocopterKARI(void);

	virtual InducedVelocity::Type GetInducedVelocityType(void) const;

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	// Relativo ai ...WithDofs
	virtual void SetInitialValue(VectorHandler& X);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

#if 0
	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(unsigned int uL, const Vec3& F, const Vec3& M, const Vec3& X);
#endif

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(const Vec3& X) const;

	// *******PER IL SOLUTORE PARALLELO********
	// Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	// utile per l'assemblaggio della matrice di connessione fra i dofs
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	// ************************************************
};

/* CyclocopterKARI - end */

class DataManager;
class MBDynParser;

extern Elem*
ReadCyclocopter(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO, 
	unsigned int uLabel,
	const StructNode* pC,
	const Mat3x3& rrot,
	const StructNode* pR);

#endif /* CYCLOCOPTER_H */

