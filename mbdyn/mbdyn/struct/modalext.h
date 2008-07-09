/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2008
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

/* Forza */

#ifndef MODALEXT_H
#define MODALEXT_H

#include <vector>
#include <string>
#include "extforce.h"
#include "modal.h"

/* ExtModalForceBase - begin */

class ExtModalForceBase {
public:
	enum {
		EMF_NONE		= 0x0U,
		EMF_RIGID		= 0x01U,
		EMF_MODAL		= 0x02U,

		EMF_RIGID_DETECT	= 0x10U,
		EMF_MODAL_DETECT	= 0x20U,

		EMF_RIGID_DETECT_MASK	= (EMF_RIGID_DETECT|EMF_RIGID),
		EMF_MODAL_DETECT_MASK	= (EMF_MODAL_DETECT|EMF_MODAL),

		EMF_DETECT_MASK		= (EMF_RIGID_DETECT|EMF_MODAL_DETECT),

		EMF_ERR			= 0x10000000U
	};

	virtual ~ExtModalForceBase(void);

	virtual unsigned
	Recv(std::istream& fin, unsigned uFlags, unsigned& uLabel,
		Vec3& f, Vec3& m, std::vector<doublereal>& f) = 0;

	virtual void
	Send(std::ostream& fout, unsigned uFlags, unsigned uLabel,
		const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
		const std::vector<doublereal>& q,
		const std::vector<doublereal>& qP) = 0;
};

/* ExtModalForceBase - end */

/* ModalExt - begin */

class ModalExt : virtual public Elem, public ExtForce {
protected:
	Modal *pModal;
	bool bOutputAccelerations;
	ExtModalForceBase *pEMF;
	unsigned uFlags;

	Vec3 F, M;
	std::vector<doublereal> f;

	// Temporary?
	std::vector<doublereal> q;
	std::vector<doublereal> qP;

	void Send(std::ostream& out, bool bAfterConvergence = false);
	void Recv(std::istream& in);
   
public:
	/* Costruttore */
	ModalExt(unsigned int uL,
		Modal *pmodal,
		bool bOutputAccelerations,
		ExtFileHandlerBase *pEFH,
		ExtModalForceBase *pEMF,
		int iCoupling,
		flag fOut);

	virtual ~ModalExt(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const { 
		return Force::EXTERNALMODAL; 
	};
 
	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
		*piNumRows = (pModal->pGetModalNode() ? 6 : 0)
			+ pModal->uGetNModes();
		*piNumCols = 1;
	};

	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);     

	virtual void Output(OutputHandler& OH) const;

	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes)
	{
		if (pModal->pGetModalNode()) {
			connectedNodes.resize(1);
			connectedNodes[0] = pModal->pGetModalNode();
		} else {
			connectedNodes.resize(0);
		}
	};
	/* ************************************************ */
};

/* ModalExt - end */

class DataManager;
class MBDynParser;

extern Elem*
ReadModalExtForce(DataManager* pDM, 
       MBDynParser& HP, 
       unsigned int uLabel);

#endif /* MODALEXT_H */

