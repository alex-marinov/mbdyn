/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "stroutput.h"

/* StructOutputManip - begin */

StructOutputManip::StructOutputManip(void)
{
	NO_OP;
}

StructOutputManip::~StructOutputManip(void)
{
	NO_OP;
}

/* StructOutputManip - end */

/* StructOutputEnd - begin */

StructOutputEnd::StructOutputEnd(unsigned uLabel, flag fOut)
: Elem(uLabel, fOut)
{
	NO_OP;
}

StructOutputEnd::~StructOutputEnd(void)
{
	NO_OP;
}

Elem::Type
StructOutputEnd::GetElemType(void) const
{
	return Elem::SOCKETSTREAM_OUTPUT;
}

void
StructOutputEnd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

SubVectorHandler&
StructOutputEnd::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(0);

	return WorkVec;
}

VariableSubMatrixHandler& 
StructOutputEnd::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

/* StructOutputEnd - end */

/* StructOutputStart - begin */

void
StructOutputStart::ManipulateInit(void)
{
	dynamic_cast<StructOutputManip *>(pElem)->ManipulateInit(data);
}

void
StructOutputStart::Manipulate(void)
{
	dynamic_cast<StructOutputManip *>(pElem)->Manipulate(data);
}

StructOutputStart::StructOutputStart(const Elem *pE)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
NestedElem(pE)
{
	NO_OP;
}

StructOutputStart::~StructOutputStart(void)
{
	NO_OP;
}

void
StructOutputStart::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{	
	ManipulateInit();
	NestedElem::SetValue(pDM, X, XP, ph);
}

void
StructOutputStart::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	Manipulate();
}

/* StructOutputStart - end */

/* StructOutput - begin */

void
StructOutput::ManipulateInit(const GeometryData& data)
{
	dynamic_cast<StructOutputManip *>(pElem)->ManipulateInit(data);
}

void
StructOutput::Manipulate(const GeometryData& data)
{
	dynamic_cast<StructOutputManip *>(pElem)->Manipulate(data);
}

StructOutput::StructOutput(const Elem *pE)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
NestedElem(pE)
{
	NO_OP;
}

StructOutput::~StructOutput(void)
{
	NO_OP;
}

/* StructOutput - end */

/* StructOutputCollect - begin */

StructOutputCollectBase::StructOutputCollectBase(const Elem *pE,
	unsigned uFlags,
	const std::vector<const StructNode *>& nodes)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutputStart(pE),
Nodes(nodes)
{
	// Initialize communication structure
	data.uFlags = uFlags;
	data.data.resize(Nodes.size());

	// Initialize labels
	for (unsigned i = 0; i < Nodes.size(); i++) {
		data.data[i].uLabel = Nodes[i]->GetLabel();
	}
}

StructOutputCollectBase::~StructOutputCollectBase(void)
{
	NO_OP;
}

void
StructOutputCollectBase::ManipulateInit(void)
{
	Manipulate_int();
	StructOutputStart::ManipulateInit();
}

void
StructOutputCollectBase::Manipulate(void)
{
	Manipulate_int();
	StructOutputStart::Manipulate();
}

void
StructOutputCollect::Manipulate_int(void)
{
	ASSERT(data.data.size() == Nodes.size());
		
	if (data.uFlags & GeometryData::X) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].X = Nodes[i]->GetXCurr();
		}
	}

	if (data.uFlags & GeometryData::R) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].R = Nodes[i]->GetRCurr();
		}
	}

	if (data.uFlags & GeometryData::V) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].V = Nodes[i]->GetVCurr();
		}
	}

	if (data.uFlags & GeometryData::W) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].W = Nodes[i]->GetWCurr();
		}
	}

	if (data.uFlags & GeometryData::XPP) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].XPP = Nodes[i]->GetXPPCurr();
		}
	}

	if (data.uFlags & GeometryData::WP) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].WP = Nodes[i]->GetWPCurr();
		}
	}
	
}

StructOutputCollect::StructOutputCollect(const Elem *pE,
	unsigned uFlags,
	const std::vector<const StructNode *>& nodes)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutputCollectBase(pE, uFlags, nodes)
{
	NO_OP;
}

StructOutputCollect::~StructOutputCollect(void)
{
	NO_OP;
}

std::ostream&
StructOutputCollect::Restart(std::ostream& out) const
{
	return out << "# StructOutputCollect(" << GetLabel() << "): "
		"not implemented yet!" << std::endl;
}

void
StructOutputCollectRelative::Manipulate_int(void)
{
	ASSERT(data.data.size() == Nodes.size());

	// Note: "Relative" in the sense that everything is oriented
	// according to the orientation of the frame of the reference node,
	// and positions are also relative to the position of the reference
	// node in that frame
	//
	//	~p   = RRef^T * (x - xRef)
	//	~R   = RRef^T * R
	//	~v   = RRef^T * v
	//	~w   = RRef^T * w
	//	~xpp = RRef^T * xpp
	//	~wp  = RRef^T * wp

	// If the really relative kinematics are needed, implement this:
	//
	//	~p   = RRef^T * (x - xRef)
	//	~R   = RRef^T * R
	//	~v   = RRef^T * (v - vRef)
	//	~w   = RRef^T * (w - wRef)
	//	~xpp = RRef^T * (xpp - xppRef)
	//	~wp  = RRef^T * (wp - wpRef)

	const Vec3& XRef = pRefNode->GetXCurr();
	Mat3x3 RRefT = pRefNode->GetRCurr().Transpose();

	if (data.uFlags & GeometryData::X) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].X = RRefT*(Nodes[i]->GetXCurr() - XRef);
		}
	}

	if (data.uFlags & GeometryData::R) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].R = RRefT*Nodes[i]->GetRCurr();
		}
	}

	if (data.uFlags & GeometryData::V) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].V = RRefT*Nodes[i]->GetVCurr();
		}
	}

	if (data.uFlags & GeometryData::W) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].W = RRefT*Nodes[i]->GetWCurr();
		}
	}

	if (data.uFlags & GeometryData::XPP) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].XPP = RRefT*Nodes[i]->GetXPPCurr();
		}
	}

	if (data.uFlags & GeometryData::WP) {
		for (unsigned i = 0; i < Nodes.size(); i++) {
			data.data[i].WP = RRefT*Nodes[i]->GetWPCurr();
		}
	}
}

StructOutputCollectRelative::StructOutputCollectRelative(const Elem *pE,
	unsigned uFlags,
	const StructNode *pRefNode,
	const std::vector<const StructNode *>& nodes)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutputCollectBase(pE, uFlags, nodes),
pRefNode(pRefNode)
{
	NO_OP;
}

StructOutputCollectRelative::~StructOutputCollectRelative(void)
{
	NO_OP;
}

std::ostream&
StructOutputCollectRelative::Restart(std::ostream& out) const
{
	return out << "# StructOutputCollectRelative(" << GetLabel() << "): "
		"not implemented yet!" << std::endl;
}

static Elem *
ReadStructOutputCollect(DataManager *pDM, MBDynParser& HP, const Elem *pNestedElem)
{
	unsigned uFlags = 0;

	uFlags = (GeometryData::X)|(GeometryData::R);

	std::cout << uFlags << std::endl;

	const StructNode *pRefNode = 0;
	if (HP.IsKeyWord("reference" "node")) {
		pRefNode = dynamic_cast<const StructNode*>(pDM->ReadNode(HP, Node::STRUCTURAL));
	}

	int n = HP.GetInt();
	if (n <= 0) {
		silent_cerr("StructOutputCollect(" << pNestedElem->GetLabel() << "): "
			"illegal node number " << n
			<< " at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	std::vector<const StructNode *> Nodes(n);
	for (int i = 0; i < n; i++ ) {
		Nodes[i] = dynamic_cast<StructNode*>(pDM->ReadNode(HP, Node::STRUCTURAL));
	}

	Elem *pEl = 0;
	if (pRefNode == 0) {
		SAFENEWWITHCONSTRUCTOR(pEl, StructOutputCollect,
			StructOutputCollect(pNestedElem, uFlags, Nodes));

	} else {
		SAFENEWWITHCONSTRUCTOR(pEl, StructOutputCollectRelative,
			StructOutputCollectRelative(pNestedElem, uFlags, pRefNode, Nodes));
	}

	return pEl;
}

/* StructOutputCollect - end */

/* StructOutputInterpBase - begin */

void
StructOutputInterpBase::ManipulateInit(const GeometryData& mb_data)
{
	
	// gather interpolation data
	ManipulateInit_int(mb_data);

	// propago i dati interpolati
	StructOutput::ManipulateInit(fem_data);

}

void
StructOutputInterpBase::Manipulate(const GeometryData& mb_data)
{
       
	// manipulate data
	Manipulate_int(mb_data);

	// propago i dati interpolati
	StructOutput::Manipulate(fem_data);

}

StructOutputInterpBase::StructOutputInterpBase(const Elem *pE)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutput(pE)
{
	NO_OP;
}

StructOutputInterpBase::~StructOutputInterpBase(void)
{
	NO_OP;
}

/* StructOutputInterpBase - end */

/* StructOutputInterp - begin */

void
StructOutputInterp::ManipulateInit_int(const GeometryData& mb_data)
{

	// use orig_data and fem_data to cook interpolation matrix
	std::cout << "Interpolation matrix computation" << std::endl;
	std::cout << "mb_data.size=" << mb_data.data.size()*(Adj.size()+1)
		<< " fem_data.size=" << fem_data.data.size() << std::endl;

	// Allocating matrix H
	SAFENEWWITHCONSTRUCTOR(pH, SpMapMatrixHandler ,
		SpMapMatrixHandler(fem_data.data.size() , mb_data.data.size()*(Adj.size()+1)));

	// Computing matrix H
	std::cout << "Interpolation matrix: " << pH->iGetNumRows() << " x " << pH->iGetNumCols() << std::endl;

#ifdef USE_X_ANN
	if (Adj.size() == 0){
		pInt->Interpolate(fem_data, mb_data, pH);
	} else {
		pInt->Interpolate_Adj(fem_data, mb_data, pH, Adj);
	}
#endif // USE_X_ANN
}

void
StructOutputInterp::Manipulate_int(const GeometryData& mb_data)
{

	fem_data.uFlags = mb_data.uFlags;

	if (mb_data.uFlags & GeometryData::X) {
		// Unpacking data	
		FullMatrixHandler* pX_mb;
		FullMatrixHandler* pX_fem;
		SAFENEWWITHCONSTRUCTOR( pX_mb, FullMatrixHandler, FullMatrixHandler(mb_data.data.size()*(Adj.size()+1),3) );
		SAFENEWWITHCONSTRUCTOR( pX_fem, FullMatrixHandler, FullMatrixHandler(fem_data.data.size(),3) );
		for (unsigned i = 0 ; i < mb_data.data.size() ; i++){
			for (unsigned j = 1; j <= 3 ; j++) {
				pX_mb->PutCoef(i*(Adj.size()+1)+1, j, mb_data.data[i].X(j) );
			}
			// X on adjoint
			// x = x0 + R*delta_x
			for (unsigned i_a = 0; i_a < Adj.size() ; i_a++){
				Vec3 tmp = mb_data.data[i].X + mb_data.data[i].R*Adj[i_a];
				for (unsigned j = 1; j <= 3 ; j++) {
					pX_mb->PutCoef(i*(Adj.size()+1)+1+i_a+1, j, tmp(j) );
				}
			}
		}	
		// Matrix Product
		pH->MatMatMul(pX_fem,*pX_mb);
		// Packing data
		for (unsigned i = 0 ; i < fem_data.data.size() ; i++){
			for (int j = 1; j <= 3 ; j++) {
				fem_data.data[i].X(j) = pX_fem->dGetCoef(i+1,j);
			}
		}
	}

	if (mb_data.uFlags & GeometryData::R) {
		// Unpacking data	
		FullMatrixHandler* pR_mb;
		FullMatrixHandler* pR_fem;
		SAFENEWWITHCONSTRUCTOR( pR_mb, FullMatrixHandler, FullMatrixHandler(mb_data.data.size()*(Adj.size()+1),9) );
		SAFENEWWITHCONSTRUCTOR( pR_fem, FullMatrixHandler, FullMatrixHandler(fem_data.data.size(),9) );
		for (unsigned i = 0 ; i < mb_data.data.size() ; i++){
			for (unsigned j = 1; j <= 3 ; j++) {
				for (unsigned k = 1; k <= 3; k++) {
					pR_mb->PutCoef(i*(Adj.size()+1)+1, 3*(j-1)+k, mb_data.data[i].R(j,k) );
					// R on adjoint
					// R = R0
					for (unsigned i_a = 0; i_a < Adj.size() ; i_a++){
						pR_mb->PutCoef(i*(Adj.size()+1)+1+i_a+1,3*(j-1)+k, mb_data.data[i].R(j,k) );	
					}
				}
			}
		}	
		// Matrix Product
		pH->MatMatMul(pR_fem,*pR_mb);
		// Packing data
		for (unsigned i = 0 ; i < fem_data.data.size() ; i++){
			for (unsigned j = 1; j <= 3 ; j++) {
				for (unsigned k = 1 ; k <= 3; k++) {
					fem_data.data[i].R(j,k) = pR_fem->dGetCoef(i+1,(j-1)*3+k);
				}
			}
		}
	}

	if (mb_data.uFlags & GeometryData::V) {
		// Unpacking data	
		FullMatrixHandler* pV_mb;
		FullMatrixHandler* pV_fem;
		SAFENEWWITHCONSTRUCTOR( pV_mb, FullMatrixHandler, FullMatrixHandler(mb_data.data.size()*(Adj.size()+1),3) );
		SAFENEWWITHCONSTRUCTOR( pV_fem, FullMatrixHandler, FullMatrixHandler(fem_data.data.size(),3) );
		for (unsigned i = 0 ; i < mb_data.data.size() ; i++){
			for (int j = 1; j <= 3 ; j++) {
				pV_mb->PutCoef(i+1, j, mb_data.data[i].V(j) );
				// V on adjoint
				// V = V0 + W*delta_x
				for (unsigned i_a = 0; i_a < Adj.size() ; i_a++){
					Vec3 tmp = mb_data.data[i].V + mb_data.data[i].W.Cross(Adj[i_a]);
					for (unsigned j = 1; j <= 3 ; j++) {
						pV_mb->PutCoef(i*(Adj.size()+1)+1+i_a+1, j, tmp(j) );
					}	
				}
			}
		}	
		// Matrix Product
		pH->MatMatMul(pV_fem,*pV_mb);
		// Packing data
		for (unsigned i = 0 ; i < fem_data.data.size() ; i++){
			for (int j = 1; j <= 3 ; j++) {
				fem_data.data[i].V(j) = pV_fem->dGetCoef(i+1,j);
			}
		}
	}

	if (mb_data.uFlags & GeometryData::W) {
		// Unpacking data	
		FullMatrixHandler* pW_mb;
		FullMatrixHandler* pW_fem;
		SAFENEWWITHCONSTRUCTOR( pW_mb, FullMatrixHandler, FullMatrixHandler(mb_data.data.size()*(Adj.size()+1),3) );
		SAFENEWWITHCONSTRUCTOR( pW_fem, FullMatrixHandler, FullMatrixHandler(fem_data.data.size(),3) );
		for (unsigned i = 0 ; i < mb_data.data.size() ; i++){
			for (unsigned j = 1; j <= 3 ; j++) {
				pW_mb->PutCoef(i+1, j, mb_data.data[i].W(j) );
				// W on adjoint
				// W = W0 
				for (unsigned i_a = 0; i_a < Adj.size() ; i_a++){
					pW_mb->PutCoef(i*(Adj.size()+1)+1+i_a+1, j,  mb_data.data[i].W(j)  );				
				}
			}
		}	
		// Matrix Product
		pH->MatMatMul(pW_fem,*pW_mb);
		// Packing data
		for (unsigned i = 0 ; i < fem_data.data.size() ; i++){
			for (int j = 1; j <= 3 ; j++) {
				fem_data.data[i].W(j) = pW_fem->dGetCoef(i+1,j);
			}
		}
	}

	if (mb_data.uFlags & GeometryData::XPP) {
		// Unpacking data	
		FullMatrixHandler* pXPP_mb;
		FullMatrixHandler* pXPP_fem;
		SAFENEWWITHCONSTRUCTOR( pXPP_mb, FullMatrixHandler, FullMatrixHandler(mb_data.data.size()*(Adj.size()+1),3) );
		SAFENEWWITHCONSTRUCTOR( pXPP_fem, FullMatrixHandler, FullMatrixHandler(fem_data.data.size(),3) );
		for (unsigned i = 0 ; i < mb_data.data.size() ; i++){
			for (int j = 1; j <= 3 ; j++) {
				pXPP_mb->PutCoef(i+1, j, mb_data.data[i].XPP(j) );
			}
			// XPP on adjoint
			// XPP = XPP0 + WP*delta_x + W*W*delta_x
			for (unsigned i_a = 0; i_a < Adj.size() ; i_a++){
				Vec3 tmp = mb_data.data[i].XPP + mb_data.data[i].WP.Cross(Adj[i_a])
				           + mb_data.data[i].W.Cross(mb_data.data[i].W.Cross(Adj[i_a]));
				for (unsigned j = 1; j <= 3 ; j++) {
					pXPP_mb->PutCoef(i*(Adj.size()+1)+1+i_a+1, j, tmp(j) );
				}	
			}
		}	
		// Matrix Product
		pH->MatMatMul(pXPP_fem,*pXPP_mb);
		// Packing data
		for (unsigned i = 0 ; i < fem_data.data.size() ; i++){
			for (int j = 1; j <= 3 ; j++) {
				fem_data.data[i].XPP(j) = pXPP_fem->dGetCoef(i+1,j);
			}
		}
	}

	if (mb_data.uFlags & GeometryData::WP) {
		// Unpacking data	
		FullMatrixHandler* pWP_mb;
		FullMatrixHandler* pWP_fem;
		SAFENEWWITHCONSTRUCTOR( pWP_mb, FullMatrixHandler, FullMatrixHandler(mb_data.data.size()*(Adj.size()+1),3) );
		SAFENEWWITHCONSTRUCTOR( pWP_fem, FullMatrixHandler, FullMatrixHandler(fem_data.data.size(),3) );
		for (unsigned i = 0 ; i < mb_data.data.size() ; i++){
			for (unsigned j = 1; j <= 3 ; j++) {
				pWP_mb->PutCoef(i+1, j, mb_data.data[i].WP(j) );
				// W on adjoint
				// W = W0 
				for (unsigned i_a = 0; i_a < Adj.size() ; i_a++){
					pWP_mb->PutCoef(i*(Adj.size()+1)+1+i_a+1, j,  mb_data.data[i].WP(j) );				
				}
			}
		}	
		// Matrix Product
		pH->MatMatMul(pWP_fem,*pWP_mb);
		// Packing data
		for (unsigned i = 0 ; i < fem_data.data.size() ; i++){
			for (int j = 1; j <= 3 ; j++) {
				fem_data.data[i].WP(j) = pWP_fem->dGetCoef(i+1,j);
			}
		}
	}
}

StructOutputInterp::StructOutputInterp(const Elem *pE,
	bool bQuad,
	int RBForder,
	int NearNodes,
	int Nadj,
	double dL)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutputInterpBase(pE)
{
#ifdef USE_X_ANN
	SAFENEWWITHCONSTRUCTOR(pInt, MLSP ,
		MLSP(NearNodes, RBForder, bQuad));
#endif // USE_X_ANN
	// Adjoint Nodes Vector
	if (Nadj != 0){
		if (Nadj==3){
			//Adj[0] = Vec3(dL/5);
			//Adj[1] = Vec3(0.,dL/5);
			//Adj[2] = Vec3(0.,0.,dL/5);
			Adj.push_back(Vec3(dL/5));
			Adj.push_back(Vec3(0.,dL/5));
			Adj.push_back(Vec3(0.,0.,dL/5));
		} else {
			Adj[0] = Vec3(dL/5);
			Adj[1] = Vec3(-dL/5);
			Adj[2] = Vec3(0.,dL/5);
			Adj[3] = Vec3(0.,-dL/5);
			Adj[4] = Vec3(0.,0.,dL/5);
			Adj[5] = Vec3(0.,0.,-dL/5);
		}
	}
}

StructOutputInterp::~StructOutputInterp(void)
{
	NO_OP;
}

std::ostream&
StructOutputInterp::Restart(std::ostream& out) const
{
	return out << "# StructOutputInterp(" << GetLabel() << "): "
		"not implemented yet!" << std::endl;
}

/* StructOutputInterp - end */

/* StructOutputInterpOP2 - begin */

void
StructOutputInterpOP2::ManipulateInit_int(const GeometryData& mb_data)
{
	// read op2 file
	std::ifstream inf(infilename.c_str());
	if (!inf) {
		silent_cerr("StructOutputWrite(" << GetLabel() << "): "
			"unable to open input file \"" << infilename << "\""
			<< std::endl);
		throw ErrGeneric();
	}
	std::cout << "letto il file..." << std::endl;
	StructOutputInterp::ManipulateInit_int(mb_data);
}

StructOutputInterpOP2::StructOutputInterpOP2(const Elem *pE,
					const std::string& infilename,
					bool bQuad,
					int RBForder,
					int NearNodes,
					int Nadj,
					double dL)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutputInterp(pE, bQuad, RBForder, NearNodes, Nadj, dL),
infilename(infilename)
{
	NO_OP;
}

StructOutputInterpOP2::~StructOutputInterpOP2(void)
{
	NO_OP;
}

std::ostream&
StructOutputInterpOP2::Restart(std::ostream& out) const
{
	return out << "# StructOutputInterpOP2(" << GetLabel() << "): "
		"not implemented yet!" << std::endl;
}

/* StrucitOutputInterpOP2 - end */

/* StructOutputInterpNative - begin */

void
StructOutputInterpNative::ManipulateInit_int(const GeometryData& mb_data)
{

	int Nnod;

	//fem_data.uFlags = GeometryData::X;
	
	// read plain ASCII file
	std::ifstream inf(infilename.c_str());
	if (!inf) {
		silent_cerr("StructOutputWrite(" << GetLabel() << "): "
			"unable to open input file \"" << infilename << "\""
			<< std::endl);
		throw ErrGeneric();
	}

	// start reading
	inf >> Nnod;

	silent_cout("StructOutputWrite(" << GetLabel() << "): "
		"Nnod=" << Nnod << std::endl);

	fem_data.data.resize(Nnod);
	for (int i = 0; i < Nnod ; i++) {
		inf >> fem_data.data[i].uLabel
			>> fem_data.data[i].X(1)
			>> fem_data.data[i].X(2)
			>> fem_data.data[i].X(3);
	}
	StructOutputInterp::ManipulateInit_int(mb_data);
}

StructOutputInterpNative::StructOutputInterpNative(const Elem *pE,
					const std::string& infilename,
					bool bQuad,
					int RBForder,
					int NearNodes,
					int Nadj,
					double dL)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutputInterp(pE, bQuad, RBForder, NearNodes, Nadj, dL),
infilename(infilename)
{
	NO_OP;
}

StructOutputInterpNative::~StructOutputInterpNative(void)
{
	NO_OP;
}

std::ostream&
StructOutputInterpNative::Restart(std::ostream& out) const
{
	return out << "# StructOutputInterpOP2(" << GetLabel() << "): "
		"not implemented yet!" << std::endl;
}

/* StructOutputInterpOP2 - end */

static Elem *
ReadStructOutputInterp(DataManager *pDM, MBDynParser& HP, const Elem *pNestedElem)
{
#ifndef USE_X_ANN
	silent_cerr("StructOutputInterp(" << pNestedElem->GetLabel() << "): "
		"need libANN..." << std::endl);
#endif // ! USE_X_ANN

	bool bOp2;
	if (HP.IsKeyWord("OP2")) {
		bOp2 = true;

	} else if (HP.IsKeyWord("native")) {
		bOp2 = false;

	} else {
		silent_cerr("StructOutputInterp(" << pNestedElem->GetLabel() << "): "
			"illegal input file format"
			" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	// collect parameters and pass to constructor
	const char *s = HP.GetFileName();
	
	if (s == 0) {
		silent_cerr("StructOutputInterp(" << pNestedElem->GetLabel() << "): "
			"unable to read input file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	std::string infilename(s);

	bool bQuad;
	if (HP.IsKeyWord("quadratic")) {
		bQuad = true;

	} else if (HP.IsKeyWord("linear")) {
		bQuad = false;

	} else {
		silent_cerr("StructOutputInterp(" << pNestedElem->GetLabel() << "): "
			"illegal interpolation order"
			" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	int RBForder = HP.GetInt();
	switch (RBForder) {
		case 0:
		case 2:
		case 4:
			break;

		default:
			silent_cerr("StructOutputInterp(" << pNestedElem->GetLabel() << "): "
				"illegal RBF order"
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric();
	}

	int NearNodes = HP.GetInt();
	if (NearNodes <= 0) {
		silent_cerr("StructOutputInterp(" << pNestedElem->GetLabel() << "): "
			"illegal interpolation order"
			" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();

	} else if (bQuad && NearNodes < 10) {
		silent_cerr("Insufficient near nodes for quadratic interpolation: "
			"increasing to 10" << std::endl);
		NearNodes = 10;

	} else if (!bQuad && NearNodes < 4) {
		silent_cerr("Insufficient near nodes for linear interpolation: "
		"increasing to 4" << std::endl);
		NearNodes = 4;
	}

	int Nadj = 0;
	double dL = 0.;
	if (HP.IsKeyWord("adjoint")) {
		Nadj = HP.GetInt();
		switch (Nadj) {
			case 0:
			case 3:
			case 6:
				break;
			default:
				silent_cerr("StructOutputInterp(" << pNestedElem->GetLabel() << "): "
				"illegal number of adjoint nodes"
				" at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric();
		}
		dL = HP.GetReal();
		if (dL <= 0.){
			silent_cerr("StructOutputInterp(" << pNestedElem->GetLabel() << "): "
			"reference length must be a positive value"
			" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}
	}

	Elem *pEl = 0;
	if (bOp2) {
		SAFENEWWITHCONSTRUCTOR(pEl, StructOutputInterpOP2,
			StructOutputInterpOP2(pNestedElem, infilename,
				bQuad, RBForder, NearNodes, Nadj, dL));

	} else {
		SAFENEWWITHCONSTRUCTOR(pEl, StructOutputInterpNative,
			StructOutputInterpNative(pNestedElem, infilename,
				bQuad, RBForder, NearNodes, Nadj, dL));
	}

	return pEl;
}

/* StructOutputWrite - begin */

StructOutputWriteBase::StructOutputWriteBase(unsigned uLabel,
	const std::string& outfilename,
	bool bNoClobberOut,
	flag fOut)
: Elem(uLabel, fOut),
StructOutputEnd(uLabel, fOut),
outfilename(outfilename),
bNoClobberOut(bNoClobberOut)
{
	NO_OP;
}

StructOutputWriteBase::~StructOutputWriteBase(void)
{
	NO_OP;
}

/* StructOutputWriteBase - end */

/* StructOutputWrite - begin */

void
StructOutputWrite::ManipulateInit(const GeometryData& data)
{
	NO_OP;
}

void
StructOutputWrite::Manipulate(const GeometryData& data)
{
	// open file
	std::ofstream outf(outfilename.c_str(),std::ios::app);
	if (!outf) {
		silent_cerr("StructOutputWrite(" << GetLabel() << "): "
			"unable to open output file \"" << outfilename << "\""
			<< std::endl);
		throw ErrGeneric();
	}

	if (iPrecision) {
		outf.precision(iPrecision);
	}
	outf.setf(std::ios::scientific);

	// header
	outf
		<< "# StructOutputWrite(" << GetLabel() << ")" << std::endl
		<< "# Label";

	if (data.uFlags & GeometryData::X) {
		outf << " X(3)";
	}
	
	if (data.uFlags & GeometryData::R) {
		outf << " R(3,3)";
	}
	
	if (data.uFlags & GeometryData::V) {
		outf << " V(3)";
	}
	
	if (data.uFlags & GeometryData::W) {
		outf << " W(3)";
	}
	
	if (data.uFlags & GeometryData::XPP) {
		outf << " XPP(3)";
	}
	
	if (data.uFlags & GeometryData::WP) {
		outf << " WP(3)";
	}
	
	if (data.uFlags & GeometryData::F) {
		outf << " F(3)";
	}
	
	if (data.uFlags & GeometryData::M) {
		outf << " M(3)";
	}
	outf << std::endl;

	// data
	for (std::vector<Geometry>::const_iterator i = data.data.begin();
		i != data.data.end(); i++)
	{
		outf << i->uLabel;

		if (data.uFlags & GeometryData::X) {
			outf << " " << i->X;
		}

		if (data.uFlags & GeometryData::R) {
			outf << " " << i->R;
		}

		if (data.uFlags & GeometryData::V) {
			outf << " " << i->V;
		}

		if (data.uFlags & GeometryData::W) {
			outf << " " << i->W;
		}

		if (data.uFlags & GeometryData::XPP) {
			outf << " " << i->XPP;
		}

		if (data.uFlags & GeometryData::WP) {
			outf << " " << i->WP;
		}

		if (data.uFlags & GeometryData::F) {
			outf << " " << i->F;
		}

		if (data.uFlags & GeometryData::M) {
			outf << " " << i->M;
		}

		outf << std::endl;
	}
}

StructOutputWrite::StructOutputWrite(unsigned uLabel,
	const std::string& outfilename,
	bool bNoClobberOut,
	int iPrecision,
	flag fOut)
: Elem(uLabel, fOut),
StructOutputWriteBase(uLabel, outfilename, bNoClobberOut, fOut),
iPrecision(iPrecision)
{
	NO_OP;
}

StructOutputWrite::~StructOutputWrite(void)
{
	NO_OP;
}

std::ostream&
StructOutputWrite::Restart(std::ostream& out) const
{
	return out << "# StructOutputWrite(" << GetLabel() << "): "
		"not implemented yet!" << std::endl;
}

static Elem *
ReadStructOutputWrite(DataManager *pDM, MBDynParser& HP, unsigned int uLabel)
{
	const char *s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("StructOutputWrite(" << uLabel << "): "
			"unable to read output file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	std::string outfilename(s);

	bool bNoClobberOut = false;
	if (HP.IsKeyWord("no" "clobber")) {
		bNoClobberOut = true;
	}

	int iPrecision = 0;
	if (HP.IsKeyWord("precision")) {
		if (!HP.IsKeyWord("default")) {
			iPrecision = HP.GetInt();
		}

		if (iPrecision < 0) {
			silent_cerr("ExtForce(" << uLabel << "): "
				"invalid precision value "
				"\"" << iPrecision << "\""
				" at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric();
		}
	}

#if 0
	flag fOut = pDM->fReadOutput(HP, Elem::SOCKETSTREAM_OUTPUT);
#endif
	flag fOut(1);

	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, StructOutputWrite,
		StructOutputWrite(uLabel, outfilename, bNoClobberOut, iPrecision, fOut));

	return pEl;
}

/* StructOutputWrite - end */

/* StructOutputWriteNASTRAN - begin */

void
StructOutputWriteNASTRAN::ManipulateInit(const GeometryData& data)
{
	NO_OP;
}

void
StructOutputWriteNASTRAN::Manipulate(const GeometryData& data)
{
	NO_OP;
}

StructOutputWriteNASTRAN::StructOutputWriteNASTRAN(unsigned uLabel,
	const std::string& outfilename,
	bool bNoClobberOut,
	flag fOut)
: Elem(uLabel, fOut),
StructOutputWriteBase(uLabel, outfilename, bNoClobberOut, fOut)
{
	NO_OP;
}

StructOutputWriteNASTRAN::~StructOutputWriteNASTRAN(void)
{
	NO_OP;
}

std::ostream&
StructOutputWriteNASTRAN::Restart(std::ostream& out) const
{
	return out << "# StructOutputWriteNASTRAN(" << GetLabel() << "): "
		"not implemented yet!" << std::endl;
}

static Elem *
ReadStructOutputWriteNASTRAN(DataManager *pDM, MBDynParser& HP, unsigned int uLabel)
{
	return 0;
}

/* StructOutputWriteNASTRAN - end */

Elem *
ReadStructOutput(DataManager *pDM, MBDynParser& HP, unsigned int uLabel)
{
	Elem *pEl = 0;

	// Put end types here
	if (HP.IsKeyWord("write")) {
		pEl = ReadStructOutputWrite(pDM, HP, uLabel);

	} else if (HP.IsKeyWord("write" "NASTRAN")) {
		pEl = ReadStructOutputWriteNASTRAN(pDM, HP, uLabel);
	}

	if (dynamic_cast<StructOutputEnd *>(pEl) == 0) {
		silent_cerr("StructOutput(" << uLabel << "): "
			"first element must be \"end\"-typed "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric();
	}

	// Add other end types here...

	while (HP.IsArg()) {
		// Put manip types here
		if (HP.IsKeyWord("interpolate")) {
			pEl = ReadStructOutputInterp(pDM, HP, pEl);
		 
		// Add other manip types here...

		// Put start types here
		} else if (HP.IsKeyWord("collect")) {
			pEl = ReadStructOutputCollect(pDM, HP, pEl);

		// Add other start types here...

		// Default
		} else  {
			silent_cerr("StructOutput(" << uLabel << "): "
				"unknown type at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric();
		}

		if (dynamic_cast<StructOutputStart *>(pEl) != 0) {
			break;
		}
	}

	return pEl;
}

