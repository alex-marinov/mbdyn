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

/*
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2020(-2020) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#include <unordered_set>
#include <vector>

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */
#include <dataman.h>
#include <userelem.h>

#include <gradient.h>
#include <matvec.h>
#include <matvecass.h>
#include "module-triangular_contact.h"

#ifdef USE_AUTODIFF

class LugreFriction
{
	public:
	LugreFriction(const DataManager* pDM);
	LugreFriction(const LugreFriction& oFrictModel)=default;
	~LugreFriction();
	void ParseInput(const Elem* pParent,
			MBDynParser& HP);
	void GetFrictionForce(const grad::Vector<doublereal, 2>& U,
			      doublereal p,
			      grad::Vector<doublereal, 2>& tau);
	template <grad::index_type N>
	void GetFrictionForce(const grad::Vector<grad::Gradient<N>, 2>& U,
			      const grad::Gradient<N>& p,
			      grad::Vector<grad::Gradient<N>, 2>& tau);

	void AfterConvergence();
	private:
	template <typename T>
	void GetFrictionForceTpl(const grad::Vector<T, 2>& U,
				 const T& p,
				 grad::Vector<T, 2>& tau);

	void SaveStictionState(const grad::Vector<doublereal, 2>& z,
			       const grad::Vector<doublereal, 2>& zP);
	template <grad::index_type N>
	void SaveStictionState(const grad::Vector<grad::Gradient<N>, 2>& z,
			       const grad::Vector<grad::Gradient<N>, 2>& zP);

	const DataManager* const pDM;
	grad::Matrix<doublereal, 2, 2> Mk, Mk2, invMk2_sigma0, Ms, Ms2, sigma0, sigma1;
	doublereal beta, vs, gamma;
	grad::Vector<doublereal, 2> zPrev, zCurr, zPPrev, zPCurr;
	doublereal tPrev, tCurr;
};

LugreFriction::LugreFriction(const DataManager* pDM)
	:pDM(pDM),
	 beta(1.),
	 vs(0.),
	 gamma(1.)
{
	tCurr = tPrev = pDM->dGetTime();
}

LugreFriction::~LugreFriction()
{

}

void LugreFriction::ParseInput(const Elem* pParent, MBDynParser& HP)
{
	using namespace grad;

	if (HP.IsKeyWord("method")) {
		if (HP.IsKeyWord("explicit" "euler")) {
			beta = 0.;
		} else if (HP.IsKeyWord("implicit" "euler")) {
			beta = 1.;
		} else if (HP.IsKeyWord("trapezoidal" "rule")) {
			beta = 0.5;
		} else if (HP.IsKeyWord("custom")) {
			beta = HP.GetReal();
		} else {
			silent_cerr("triangular surface contact("
				    << pParent->GetLabel()
				    << "): keyword \"explicit euler\", "
				    "\"implicit euler\", "
				    "\"trapezoidal rule\" or \"custom\" "
				    "expected at line "
				    << HP.GetLineData() << std::endl);

			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!(HP.IsKeyWord("coulomb" "friction" "coefficient")
	      || HP.IsKeyWord("coulomb" "friction" "coefficient" "x"))) {
		silent_cerr("triangular surface contact("
			    << pParent->GetLabel()
			    << "): keyword \"coulomb friction coefficient\""
			    " or \"coulomb friction coefficient x\" expected at line "
			    << HP.GetLineData() << std::endl);

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const doublereal mukx = HP.GetReal();

	doublereal muky;

	if (HP.IsKeyWord("coulomb" "friction" "coefficient" "y")) {
		muky = HP.GetReal();
	} else {
		muky = mukx;
	}

	doublereal musx, musy;

	if (HP.IsKeyWord("static" "friction" "coefficient")
	    || HP.IsKeyWord("static" "friction" "coefficient" "x")) {
		musx = HP.GetReal();

		if (HP.IsKeyWord("static" "friction" "coefficient" "y")) {
			musy = HP.GetReal();
		} else {
			musy = musx;
		}
	} else {
		musx = mukx;
		musy = muky;
	}

	if (HP.IsKeyWord("sliding" "velocity" "coefficient")) {
		vs = HP.GetReal();
	} else {
		vs = 1.;
	}

	if (HP.IsKeyWord("sliding" "velocity" "exponent")) {
		gamma = HP.GetReal();
	} else {
		gamma = 1.;
	}

	if (!(HP.IsKeyWord("micro" "slip" "stiffness") || HP.IsKeyWord("micro" "slip" "stiffness" "x"))) {
		silent_cerr("triangular surface contact("
			    << pParent->GetLabel()
			    << "): keyword \"micro slip stiffness\" or \"micro slip stiffness x\" expected at line "
			    << HP.GetLineData() << std::endl);

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const doublereal sigma0x = HP.GetReal();

	doublereal sigma0y;

	if (HP.IsKeyWord("micro" "slip" "stiffness" "y")) {
		sigma0y = HP.GetReal();
	} else {
		sigma0y = sigma0x;
	}

	doublereal sigma1x, sigma1y;

	if (HP.IsKeyWord("micro" "slip" "damping") || HP.IsKeyWord("micro" "slip" "damping" "x")) {
		sigma1x = HP.GetReal();

		if (HP.IsKeyWord("micro" "slip" "damping" "y")) {
			sigma1y = HP.GetReal();
		} else {
			sigma1y = sigma1x;
		}
	} else {
		sigma1x = 0.;
		sigma1y = 0.;
	}

	Mk(1, 1) = mukx;
	Mk(2, 2) = muky;

	Mk2 = Mk * Mk;

	Ms(1, 1) = musx;
	Ms(2, 2) = musy;

	Ms2 = Ms * Ms;

	sigma0(1, 1) = sigma0x;
	sigma0(2, 2) = sigma0y;

	sigma1(1, 1) = sigma1x;
	sigma1(2, 2) = sigma1y;

	invMk2_sigma0 = Inv(Mk2) * sigma0;
}

void LugreFriction::GetFrictionForce(const grad::Vector<doublereal, 2>& U,
				     doublereal p,
				     grad::Vector<doublereal, 2>& tau)
{
	GetFrictionForceTpl(U, p, tau);
}

template <grad::index_type N>
void LugreFriction::GetFrictionForce(const grad::Vector<grad::Gradient<N>, 2>& U,
				     const grad::Gradient<N>& p,
				     grad::Vector<grad::Gradient<N>, 2>& tau)
{
	GetFrictionForceTpl(U, p, tau);
}

void LugreFriction::AfterConvergence()
{
	tPrev = tCurr;
	zPrev = zCurr;
	zPPrev = zPCurr;
}

template <typename T>
void LugreFriction::GetFrictionForceTpl(const grad::Vector<T, 2>& U,
					const T& p,
					grad::Vector<T, 2>& tau)
{
	using namespace grad;

	typedef Matrix<doublereal, 2, 2> CMat2x2;
	typedef Vector<doublereal, 2> CVec2;
	typedef Matrix<T, 2, 2> VMat2x2;
	typedef Vector<T, 2> VVec2;

	const VVec2 Ueff = U * doublereal(p > 0.);

	const T norm_Ueff = Dot(Ueff, Ueff);

	T kappa;

	if (norm_Ueff == 0.) {
		kappa = 0.;
	} else {
		const VVec2 Mk_U = Mk * Ueff;
		const VVec2 Ms_U = Ms * Ueff;
		const VVec2 Mk2_U = Mk2 * Ueff;
		const VVec2 Ms2_U = Ms2 * Ueff;
		const T norm_Mk2_U = sqrt(Dot(Mk2_U, Mk2_U));
		const T a0 = norm_Mk2_U / sqrt(Dot(Mk_U, Mk_U));
		const T a1 = sqrt(Dot(Ms2_U, Ms2_U)) / sqrt(Dot(Ms_U, Ms_U));
		const T g = a0 + (a1 - a0) * exp(-pow(sqrt(norm_Ueff) / vs, gamma));

		kappa = norm_Mk2_U / g;
	}

	tCurr = pDM->dGetTime();

	const doublereal dt = tCurr - tPrev;
	const VMat2x2 A = invMk2_sigma0 * kappa;
	const VMat2x2 B = A * (beta * dt) + CMat2x2(1., 0., 0., 1);

	const Vector<T, 2> zP = Inv(B) * VVec2(Ueff - A * CVec2(zPrev + zPPrev * ((1 - beta) * dt)));
	const Vector<T, 2> z = zPrev + (zP * beta + zPPrev * (1 - beta)) * dt;

	SaveStictionState(z, zP);

	tau = (sigma0 * z + sigma1 * zP) * p;
}

void LugreFriction::SaveStictionState(const grad::Vector<doublereal, 2>& z,
				      const grad::Vector<doublereal, 2>& zP)
{
	zCurr = z;
	zPCurr = zP;
}

template <grad::index_type N>
void LugreFriction::SaveStictionState(const grad::Vector<grad::Gradient<N>, 2>&,
				      const grad::Vector<grad::Gradient<N>, 2>&)
{
	// Do Nothing
}

class TriangSurfContact: virtual public Elem, public UserDefinedElem
{
	public:
	TriangSurfContact(unsigned uLabel, const DofOwner *pDO,
			  DataManager* pDM, MBDynParser& HP);
	virtual ~TriangSurfContact(void);
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
	template <typename T>
	inline void
	AssRes(grad::GradientAssVec<T>& WorkVec,
	       doublereal dCoef,
	       const grad::GradientVectorHandler<T>& XCurr,
	       const grad::GradientVectorHandler<T>& XPrimeCurr,
	       enum grad::FunctionCall func);
	template <typename T>
	inline void
	InitialAssRes(grad::GradientAssVec<T>& WorkVec,
		      const grad::GradientVectorHandler<T>& XCurr,
		      enum grad::FunctionCall func);
	virtual void AfterConvergence(const VectorHandler& X,
				      const VectorHandler& XP);
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	std::ostream& Restart(std::ostream& out) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		      const VectorHandler& XCurr);
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
	virtual void SetValue(DataManager*, VectorHandler&, VectorHandler&, SimulationEntity::Hints*);
	virtual unsigned int iGetInitialNumDof() const;
	private:
	struct TargetVertex {
		Mat3x3 R = ::Zero3x3;
		Vec3 o = ::Zero3;
	};

	struct TargetFace {
		static constexpr integer iNumVertices = 3;

		TargetFace(const DataManager* pDM)
			:oc(::Zero3),
			 oFrictionModel(pDM) {
		}

		Vec3 oc;
		std::array<TargetVertex, iNumVertices> rgVert;
		LugreFriction oFrictionModel;
	};

	struct ContactNode {
		ContactNode(const StructNode* pNode,
			    const Vec3& vOffset,
			    doublereal dRadius,
			    integer iNumFaces)
			:pContNode(pNode),
			 vOffset(vOffset),
			 dRadius(dRadius)
			{
				rgTargetFaces.reserve(iNumFaces);
			}

		void Reset() {
			rgTargetFaces.clear();
		}

		const StructNode* const pContNode;
		std::vector<TargetFace*> rgTargetFaces;
		const Vec3 vOffset;
		const doublereal dRadius;
	};

	void ContactSearch();

	doublereal GetContactForce(doublereal dz, grad::LocalDofMap*) const {
		return (*pCL)(dz);
	}

	template <grad::index_type N>
	grad::Gradient<N> GetContactForce(const grad::Gradient<N>& dz, grad::LocalDofMap* pDofMap) const {
		using namespace grad;

		Gradient<N> F2n;

		F2n.SetValuePreserve((*pCL)(dz.dGetValue()));

		const doublereal dF_dz = pCL->ComputeDiff(dz.dGetValue());

		const index_type iStartIndex = dz.iGetStartIndexLocal();
		const index_type iEndIndex = dz.iGetEndIndexLocal();

		F2n.DerivativeResizeReset(pDofMap, iStartIndex, iEndIndex, MapVectorBase::LOCAL, 0.);

		for (index_type i = iStartIndex; i < iEndIndex; ++i) {
			F2n.SetDerivativeLocal(i, dF_dz * dz.dGetDerivativeLocal(i));
		}

		return F2n;
	}

	template <typename T>
	void
	UnivAssRes(grad::GradientAssVec<T>& WorkVec,
		   doublereal dCoef,
		   const grad::GradientVectorHandler<T>& XCurr,
		   enum grad::FunctionCall func);

	std::vector<TargetFace> rgTargetMesh;
	std::vector<ContactNode> rgContactMesh;
	const StructNode* pTargetNode;
	doublereal dSearchRadius;
	std::unordered_set<ContactNode*> oNodeSet;
	grad::LocalDofMap oDofMap;
	const DifferentiableScalarFunction* pCL;
	static constexpr grad::index_type iNumDofGradient = 13;

	enum class FrictionModel {
		None,
		Lugre
	} eFrictionModel;
};

TriangSurfContact::TriangSurfContact(unsigned uLabel, const DofOwner *pDO,
				     DataManager* pDM, MBDynParser& HP)
	:Elem(uLabel, flag(0)),
	 UserDefinedElem(uLabel, pDO),
	 pTargetNode(nullptr),
	 dSearchRadius(std::numeric_limits<doublereal>::max()),
	 pCL(nullptr),
	 eFrictionModel(FrictionModel::None)
{
	if (HP.IsKeyWord("help"))
	{
		silent_cout("Module: triangular surface contact\n"
			    "This element implements unilateral contact between an arbitrary rigid body and a set of nodes\n");

		if (!HP.IsArg()) {
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!HP.IsKeyWord("target" "node")) {
		silent_cerr("triangular surface contact(" << uLabel << "): keyword \"target node\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pTargetNode = pDM->ReadNode<StructNode, Node::STRUCTURAL>(HP);

	if (!HP.IsKeyWord("penalty" "function")) {
		silent_cerr("triangular surface contact(" << uLabel << "): keyword \"penalty function\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pCL = dynamic_cast<const DifferentiableScalarFunction*>(ParseScalarFunction(HP, pDM));

	if (!pCL) {
		silent_cerr("triangular surface contact(" << uLabel << "): scalar function is not differentiable at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("search" "radius")) {
		dSearchRadius = HP.GetReal();
	}

	if (dSearchRadius <= 0.) {
		silent_cerr("triangular surface contact(" << uLabel << "): invalid value for search radius at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	TargetFace oCurrFace(pDM);

	if (HP.IsKeyWord("friction" "model")) {
		if (HP.IsKeyWord("lugre")) {
			eFrictionModel = FrictionModel::Lugre;
			oCurrFace.oFrictionModel.ParseInput(this, HP);
		} else if (!HP.IsKeyWord("none")) {
			silent_cerr("triangular surface(" << uLabel << "): keyword \"lugre\" or \"none\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	ReferenceFrame oRefTarget(pTargetNode);

	if (!HP.IsKeyWord("number" "of" "target" "vertices")) {
		silent_cerr("triangular surface contact(" << uLabel << "): keyword \"number of target vertices\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const integer iNumVertices = HP.GetInt();

	if (iNumVertices <= 0) {
		silent_cerr("triangular surface contact(" << uLabel << "): at least one vertex is required at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::vector<Vec3> rgVertices;

	rgVertices.reserve(iNumVertices);

	for (integer i = 0; i < iNumVertices; ++i) {
		rgVertices.emplace_back(HP.GetPosRel(oRefTarget));
	}

	if (!HP.IsKeyWord("number" "of" "target" "faces")) {
		silent_cerr("triangular surface contact(" << uLabel << "): keyword \"number of target faces\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const integer iNumFaces = HP.GetInt();

	if (iNumFaces <= 0) {
		silent_cerr("triangular surface contact(" << uLabel << "): at least one face is required at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	rgTargetMesh.reserve(iNumFaces);

	for (integer i = 0; i < iNumFaces; ++i) {
		oCurrFace.oc = ::Zero3;

		for (integer j = 0; j < TargetFace::iNumVertices; ++j) {
			const integer iVertex = HP.GetInt() - 1;

			if (iVertex < 0 || iVertex >= iNumVertices) {
				silent_cerr("triangular surface contact(" << uLabel << "): vertex index out of range (1:" << iNumVertices << ") at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			oCurrFace.rgVert[j].o = rgVertices[iVertex];
			oCurrFace.oc += rgVertices[iVertex];
		}

		oCurrFace.oc /= TargetFace::iNumVertices;

		for (integer j = 0; j < TargetFace::iNumVertices; ++j) {
			static const integer e1_idx[3][2] = {{1, 0}, {2, 1}, {0, 2}};
			static const integer e2_idx[3][2] = {{2, 0}, {0, 1}, {1, 2}};
			Vec3 e1 = oCurrFace.rgVert[e1_idx[j][1]].o - oCurrFace.rgVert[e1_idx[j][0]].o;
			Vec3 e2 = oCurrFace.rgVert[e2_idx[j][1]].o - oCurrFace.rgVert[e2_idx[j][0]].o;
			Vec3 e3 = e1.Cross(e2);
			e2 = e3.Cross(e1);

			doublereal e1n = e1.Norm(), e2n = e2.Norm(), e3n = e3.Norm();

			if (!(e1n && e2n && e3n)) {
				silent_cerr("triangular surface contact(" << uLabel << "): singular face geometry detected at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			oCurrFace.rgVert[j].R = Mat3x3(e1 / e1n, e2 / e2n, e3 / e3n);
		}

		rgTargetMesh.emplace_back(oCurrFace);
	}

	if (!HP.IsKeyWord("number" "of" "contact" "nodes")) {
		silent_cerr("triangular surface contact(" << uLabel << "): keyword \"number of contact nodes\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const integer iNumContactNodes = HP.GetInt();

	if (iNumContactNodes <= 0) {
		silent_cerr("triangular surface contact(" << uLabel << "): invalid number of contact nodes at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	rgContactMesh.reserve(iNumContactNodes);

	for (integer i = 0; i < iNumContactNodes; ++i) {
		const StructNode* pContNode = pDM->ReadNode<StructNode, Node::STRUCTURAL>(HP);
		const Vec3 o1 = HP.GetPosRel(ReferenceFrame(pContNode));
		const doublereal dRadius = HP.GetReal();

		if (dRadius < 0.) {
			silent_cerr("triangular surface contact(" << uLabel << "): invalid value for radius at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		rgContactMesh.emplace_back(pContNode, o1, dRadius, iNumFaces);
	}




	oNodeSet.reserve(iNumContactNodes);
}

TriangSurfContact::~TriangSurfContact(void)
{

}

void TriangSurfContact::SetValue(DataManager*, VectorHandler&, VectorHandler&, SimulationEntity::Hints*)
{

}

unsigned int TriangSurfContact::iGetInitialNumDof() const
{
	return 0;
}

std::ostream& TriangSurfContact::Restart(std::ostream& out) const
{
	return out;
}

int TriangSurfContact::iGetNumConnectedNodes(void) const
{
	return rgContactMesh.size() + 1;
}

void TriangSurfContact::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.clear();
	connectedNodes.reserve(iGetNumConnectedNodes());

	connectedNodes.push_back(pTargetNode);

	for (const auto& rCont: rgContactMesh) {
		connectedNodes.push_back(rCont.pContNode);
	}
}

void TriangSurfContact::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	constexpr integer iNumRes = iNumDofGradient;

	*piNumRows = -(rgContactMesh.size() * iNumRes);
	*piNumCols = iNumRes;
}

void
TriangSurfContact::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	WorkSpaceDim(piNumRows, piNumCols);
}

void TriangSurfContact::ContactSearch()
{
	const Vec3& X2 = pTargetNode->GetXCurr();
	const Mat3x3& R2 = pTargetNode->GetRCurr();

	for (auto& rNode: rgContactMesh) {
		rNode.Reset();
	}

	for (auto& rFace: rgTargetMesh) {
		oNodeSet.clear();

		for (auto& rNode: rgContactMesh) {
			const Vec3& X1 = rNode.pContNode->GetXCurr();
			const Mat3x3& R1 = rNode.pContNode->GetRCurr();
			const Vec3& o1 = rNode.vOffset;
			doublereal dDist = (X2 + R2 * rFace.oc - X1 - R1 * o1).Norm() - rNode.dRadius;

			if (dDist <= dSearchRadius) {
				oNodeSet.insert(&rNode);
			}
		}

		for (auto pNode: oNodeSet) {
			const Vec3& X1 = pNode->pContNode->GetXCurr();
			const Mat3x3& R1 = pNode->pContNode->GetRCurr();
			const Vec3& o1 = pNode->vOffset;
			const Vec3 dX = R2.MulTV(X1 + R1 * o1 - X2);

			bool bInsert = true;

			for (auto& rVert: rFace.rgVert) {
				const Vec3& o2i = rVert.o;
				const Mat3x3& R2i = rVert.R;
				const doublereal dy = R2i.GetCol(2).Dot(dX - o2i);

				if (dy > 0) {
					bInsert = false;
					break;
				}
			}

			if (bInsert) {
				pNode->rgTargetFaces.push_back(&rFace);
			}
		}
	}
}

template <typename T>
inline void
TriangSurfContact::AssRes(grad::GradientAssVec<T>& WorkVec,
			  doublereal dCoef,
			  const grad::GradientVectorHandler<T>& XCurr,
			  const grad::GradientVectorHandler<T>& XPrimeCurr,
			  enum grad::FunctionCall func)
{
	UnivAssRes(WorkVec, dCoef, XCurr, func);
}

template <typename T>
inline void
TriangSurfContact::InitialAssRes(grad::GradientAssVec<T>& WorkVec,
				 const grad::GradientVectorHandler<T>& XCurr,
				 enum grad::FunctionCall func)
{
	UnivAssRes(WorkVec, 1., XCurr, func);
}

template <typename T>
inline void
TriangSurfContact::UnivAssRes(grad::GradientAssVec<T>& WorkVec,
			      doublereal dCoef,
			      const grad::GradientVectorHandler<T>& XCurr,
			      enum grad::FunctionCall func)
{
	using namespace grad;

	if (func & RESIDUAL_FLAG) {
		ContactSearch();
	}

	const integer iFirstIndex2 = pTargetNode->iGetFirstMomentumIndex();

	Vector<T, 3> X1, X2, F1, M1, F2, M2;
	Vector<T, 3> XP1, XP2, omega1, omega2, dV;
	Vector<T, 2> U, tau;
	Matrix<T, 3, 3> R1, R2;

	for (const auto& rNode: rgContactMesh) {
		F1 = M1 = F2 = M2 = ::Zero3;

		oDofMap.Reset();

		for (auto pTargetFace: rNode.rgTargetFaces) {
			rNode.pContNode->GetXCurr(X1, dCoef, func, &oDofMap);
			rNode.pContNode->GetRCurr(R1, dCoef, func, &oDofMap);

			pTargetNode->GetXCurr(X2, dCoef, func, &oDofMap);
			pTargetNode->GetRCurr(R2, dCoef, func, &oDofMap);

			const Vec3& o1 = rNode.vOffset;
			const Mat3x3& R2i = pTargetFace->rgVert[0].R;
			const Vec3& o2i = pTargetFace->rgVert[0].o;

			const Vector<T, 3> n3 = R2 * R2i.GetCol(3);
			const Vector<T, 3> l1 = R1 * o1;
			const Vector<T, 3> l2 = X1 + l1 - X2;
			const T dz = Dot(n3, Vector<T, 3>(l2 - R2 * o2i));
			const Vector<T, 3> l2c = l2 + n3 * dz;
			const Vector<T, 3> l1c = l1 - n3 * dz;
			const T pz = dz - rNode.dRadius;
			const T F2in = GetContactForce(pz, &oDofMap);
			Vector<T, 3> F2i = n3 * F2in;

			if (eFrictionModel == FrictionModel::Lugre) {
				pTargetNode->GetVCurr(XP2, dCoef, func, &oDofMap);
				pTargetNode->GetWCurr(omega2, dCoef, func, &oDofMap);
				rNode.pContNode->GetVCurr(XP1, dCoef, func, &oDofMap);
				rNode.pContNode->GetWCurr(omega1, dCoef, func, &oDofMap);

				dV = Transpose(R2) * Vector<T, 3>(XP1 + Cross(omega1, l1c) - XP2 - Cross(omega2, l2c));

				for (index_type i = 1; i <= 2; ++i) {
					U(i) = Dot(Direct(R2i.GetCol(i)), dV);
				}

				pTargetFace->oFrictionModel.GetFrictionForce(U, T(fabs(F2in)), tau);

				for (index_type i = 1; i <= 2; ++i) {
					F2i += R2 * (Direct(R2i.GetCol(i)) * tau(i));
				}
			}

			F2 += F2i;
			M2 += Cross(l2c, F2i);
			F1 -= F2i;
			M1 -= Cross(l1c, F2i);
		}

		const integer iFirstIndex1 = rNode.pContNode->iGetFirstMomentumIndex();

		WorkVec.AddItem(iFirstIndex1 + 1, F1);
		WorkVec.AddItem(iFirstIndex1 + 4, M1);
		WorkVec.AddItem(iFirstIndex2 + 1, F2);
		WorkVec.AddItem(iFirstIndex2 + 4, M2);
	}
}

VariableSubMatrixHandler&
TriangSurfContact::AssJac(VariableSubMatrixHandler& WorkMat,
			  doublereal dCoef,
			  const VectorHandler& XCurr,
			  const VectorHandler& XPrimeCurr)
{
	using namespace grad;

	GradientAssVec<Gradient<iNumDofGradient> >::AssJac(this, WorkMat.SetSparse(), dCoef, XCurr, XPrimeCurr, REGULAR_JAC, &oDofMap);

	return WorkMat;
}

SubVectorHandler&
TriangSurfContact::AssRes(SubVectorHandler& WorkVec,
			  doublereal dCoef,
			  const VectorHandler& XCurr,
			  const VectorHandler& XPrimeCurr)
{
	using namespace grad;

	GradientAssVec<doublereal>::AssRes(this, WorkVec, dCoef, XCurr, XPrimeCurr, REGULAR_RES);

	return WorkVec;
}

VariableSubMatrixHandler&
TriangSurfContact::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				 const VectorHandler& XCurr)
{
	using namespace grad;

	GradientAssVec<Gradient<iNumDofGradient> >::InitialAssJac(this, WorkMat.SetSparse(), XCurr, INITIAL_ASS_JAC, &oDofMap);

	return WorkMat;
}

SubVectorHandler&
TriangSurfContact::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
	using namespace grad;

	GradientAssVec<doublereal>::InitialAssRes(this, WorkVec, XCurr, INITIAL_ASS_RES);

	return WorkVec;
}

void TriangSurfContact::AfterConvergence(const VectorHandler& X,
					 const VectorHandler& XP)
{
	if (eFrictionModel != FrictionModel::None) {
		for (auto& oTargetFace: rgTargetMesh) {
			oTargetFace.oFrictionModel.AfterConvergence();
		}
	}
}
#endif

bool triangular_contact_set(void)
{
#ifdef USE_AUTODIFF
	UserDefinedElemRead *rf = new UDERead<TriangSurfContact>;

	if (!SetUDE("triangular" "surface" "contact", rf))
	{
		delete rf;
		return false;
	}
#endif

	return true;
}

#ifndef STATIC_MODULES

extern "C"
{

    int module_init(const char *module_name, void *pdm, void *php)
    {
	    if (!triangular_contact_set())
	    {
		    silent_cerr("contact: "
				"module_init(" << module_name << ") "
				"failed" << std::endl);

		    return -1;
	    }

	    return 0;
    }

}

#endif // ! STATIC_MODULE
