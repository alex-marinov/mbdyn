/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

/* Elementi elettrici */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "strnode.h"
#include "elecnode.h"
#include "displacement.h"

/* DispMeasure - begin */

/* Costruttore */
DispMeasure::DispMeasure(unsigned int uL, const DofOwner* pDO,
	const StructNode* pS1, const StructNode* pS2,
	const ScalarDifferentialNode* pA,
	const Vec3& Tmpf1, const Vec3& Tmpf2,
	flag fOut)
: Elem(uL, fOut),
Electric(uL, pDO, fOut),
pStrNode1(pS1), pStrNode2(pS2), pAbsNode(pA),
f1(Tmpf1), f2(Tmpf2)
{
	ASSERT(pStrNode1 != NULL);
	ASSERT(pStrNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pStrNode2 != NULL);
	ASSERT(pStrNode2->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pAbsNode != NULL);
	ASSERT(pAbsNode->GetNodeType() == Node::ABSTRACT);
}

/* Distruttore banale */
DispMeasure::~DispMeasure(void)
{
	NO_OP;
}

/* Contributo al file di restart */
std::ostream&
DispMeasure::Restart(std::ostream& out) const
{
	Electric::Restart(out) << ", displacement, "
		<< pStrNode1->GetLabel() << ", "
		"reference, node, ", f1.Write(out, ", ") << ", "
		<< pStrNode2->GetLabel() << ", "
		"reference, node, ", f2.Write(out, ", ") << ", "
		<< pAbsNode->GetLabel() << ';'
		<< std::endl;
	return out;
}

/* Costruisce il contributo allo jacobiano */
VariableSubMatrixHandler&
DispMeasure::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DispMeasure::AssJac()" << std::endl);

	/* Casting di WorkMat */
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.ResizeReset(1, 0);

	/* Indici delle equazioni */
	integer iAbstractIndex = pAbsNode->iGetFirstIndex()+1;

	WM.PutItem(1, iAbstractIndex, iAbstractIndex, dCoef);

	return WorkMat;
}

/* Costruisce il contributo al residuo */
SubVectorHandler&
DispMeasure::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DispMeasure::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	WorkVec.Resize(1);

	integer iAbstractIndex = pAbsNode->iGetFirstIndex() + 1;

	WorkVec.PutRowIndex(1, iAbstractIndex);

	Vec3 x1 = pStrNode1->GetXCurr() + pStrNode1->GetRCurr()*f1;
	Vec3 x2 = pStrNode2->GetXCurr() + pStrNode2->GetRCurr()*f2;

	doublereal a = pAbsNode->dGetX();

	doublereal d = (x2 - x1).Norm();

	WorkVec.PutCoef(1, d - a);

	return WorkVec;
}

void
DispMeasure::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 1;
	*piNumCols = 1;
}

/* Setta i valori iniziali delle variabili (e fa altre cose)
 * prima di iniziare l'integrazione */
void
DispMeasure::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	integer iIndex = pAbsNode->iGetFirstIndex() + 1;

	/* distanza */
	Vec3 ff1 = pStrNode1->GetRCurr()*f1;
	Vec3 ff2 = pStrNode2->GetRCurr()*f2;

	Vec3 x1 = pStrNode1->GetXCurr()+ff1;
	Vec3 x2 = pStrNode2->GetXCurr()+ff2;

	doublereal d = (x2 - x1).Norm();

	const_cast<ScalarDifferentialNode *>(pAbsNode)->SetX(d);
	X.PutCoef(iIndex, d);

	/* velocita' */
	Vec3 v1 = pStrNode1->GetVCurr() + (pStrNode1->GetWCurr()).Cross(ff1);
	Vec3 v2 = pStrNode2->GetVCurr() + (pStrNode2->GetWCurr()).Cross(ff2);

	doublereal v = (v2 - v1).Norm();

	const_cast<ScalarDifferentialNode *>(pAbsNode)->SetXPrime(v);

	XP.PutCoef(iIndex, v);
}

/* DispMeasure - end */
