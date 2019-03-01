/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2017
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
	enum BitMask {
		EMF_NONE		= 0x0U,
		EMF_RIGID		= 0x1U,
		EMF_MODAL		= 0x2U,
		EMF_ALL			= (EMF_RIGID | EMF_MODAL),

		EMF_RIGID_DETECT	= 0x10U,
		EMF_MODAL_DETECT	= 0x20U,

		EMF_RIGID_DETECT_MASK	= (EMF_RIGID_DETECT | EMF_RIGID),
		EMF_MODAL_DETECT_MASK	= (EMF_MODAL_DETECT | EMF_MODAL),

		EMF_DETECT_MASK		= (EMF_RIGID_DETECT | EMF_MODAL_DETECT),

		EMF_ERR			= 0x10000000U
	};

	virtual ~ExtModalForceBase(void);

	virtual bool
	Prepare(ExtFileHandlerBase *pEFH, unsigned uLabel, bool bRigid, unsigned uModes) = 0;

	/*
	 * Interface:
	 *	f	force oriented wrt/ the global reference frame
	 *	m	moment oriented wrt/ the global reference frame,
	 *		the pole is the node
	 *
	 *	fv	modal forces
	 */
	virtual unsigned
	Recv(ExtFileHandlerBase *pEFH, unsigned uFlags, unsigned& uLabel,
		Vec3& f, Vec3& m, std::vector<doublereal>& fv) = 0;

	/*
	 * Interface:
	 *	x	position of node wrt/ the global reference frame
	 *	R	orientation of node wrt/ the global reference frame
	 *	v	velocity of node wrt/ the global reference frame
	 *	w	angular velocity of node wrt/ the global reference frame
	 *
	 *	q	modal coordinates
	 *	qP	derivatives of modal coordinates
	 */
	virtual void
	Send(ExtFileHandlerBase *pEFH, unsigned uFlags, unsigned uLabel,
		const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
		const std::vector<doublereal>& q,
		const std::vector<doublereal>& qP) = 0;
};

/* ExtModalForceBase - end */

/* ExtModalForce - begin */

class ExtModalForce : public ExtModalForceBase {
public:
	ExtModalForce(void);
	virtual ~ExtModalForce(void);

	virtual bool
	Prepare(ExtFileHandlerBase *pEFH, unsigned uLabel, bool bRigid, unsigned uModes);

	virtual unsigned
	Recv(ExtFileHandlerBase *pEFH, unsigned uFlags, unsigned& uLabel,
		Vec3& f, Vec3& m, std::vector<doublereal>& fv);

	virtual void
	Send(ExtFileHandlerBase *pEFH, unsigned uFlags, unsigned uLabel,
		const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
		const std::vector<doublereal>& q,
		const std::vector<doublereal>& qP);

protected:
	virtual unsigned
	RecvFromStream(std::istream& inf, unsigned uFlags, unsigned& uLabel,
		Vec3& f, Vec3& m, std::vector<doublereal>& fv);
	virtual unsigned
	RecvFromFileDes(int infd, int recv_flags, unsigned uFlags, unsigned& uLabel,
		Vec3& f, Vec3& m, std::vector<doublereal>& fv);

	virtual void
	SendToStream(std::ostream& outf, unsigned uFlags, unsigned uLabel,
		const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
		const std::vector<doublereal>& q,
		const std::vector<doublereal>& qP);
	virtual void
	SendToFileDes(int outfd, int send_flags, unsigned uFlags, unsigned uLabel,
		const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
		const std::vector<doublereal>& q,
		const std::vector<doublereal>& qP);
};

/* ExtModalForce - end */

/* ModalExt - begin */

class ModalExt : virtual public Elem, public ExtForce {
protected:
	const Modal *pModal;
	// if pModal != 0 then pNode == pModal->pGetModalNode()
	const StructNode *pNode;
	bool bOutputAccelerations;
	ExtModalForceBase *pEMF;
	unsigned uFlags;

	Vec3 F, M;
	std::vector<doublereal> f;

	// Temporary?
	std::vector<doublereal> q;
	std::vector<doublereal> qP;

	bool Prepare(ExtFileHandlerBase *pEFH);
	void Send(ExtFileHandlerBase *pEFH, ExtFileHandlerBase::SendWhen when);
	void Recv(ExtFileHandlerBase *pEFH);

public:
	/* Costruttore */
	ModalExt(unsigned int uL,
		DataManager *pDM,
		const Modal *pmodal,
		const StructNode *pnode,
		bool bOutputAccelerations,
		ExtFileHandlerBase *pEFH,
		ExtModalForceBase *pEMF,
		bool bSendAfterPredict,
		int iCoupling,
		ExtModalForceBase::BitMask bm,
		flag fOut);

	virtual ~ModalExt(void);

	/* Tipo di forza */
	virtual Force::Type GetForceType(void) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual void Output(OutputHandler& OH) const;

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
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

