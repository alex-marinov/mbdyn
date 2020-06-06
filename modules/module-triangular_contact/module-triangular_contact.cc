/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2020
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

class LugreData
{
public:
     friend class LugreState;

     LugreData(const DataManager* pDM);
     LugreData(const LugreData& oData)=default;

     void ParseInput(const Elem* pParent,
		     MBDynParser& HP);

     doublereal dGetTime() const {
	  return pDM->dGetTime();
     }

private:
     const DataManager* const pDM;
     grad::Matrix<doublereal, 2, 2> Mk, Mk2, invMk2_sigma0, Ms, Ms2, sigma0, sigma1;
     doublereal beta, vs, gamma;
};

class LugreState
{
public:
     LugreState(const LugreData* pData)
	  :pData(pData) {}

     template <typename T>
     void GetFrictionForce(doublereal dt,
			   const grad::Vector<T, 2>& U,
			   const T& p,
			   grad::Vector<T, 2>& tau);

     void AfterConvergence();

     LugreState& operator=(const LugreState& oState) {
	  ASSERT(oState.pData == pData);

	  zPrev = oState.zPrev;
	  zCurr = oState.zCurr;
	  zPPrev = oState.zPPrev;
	  zPCurr = oState.zPCurr;

	  return *this;
     }

     void Project(const Mat3x3& R1, const Mat3x3& R2, const LugreState& oState2) {
	  for (integer i = 1; i <= 2; ++i) {
	       zCurr(i) = zPCurr(i) = 0.;
	  }

	  for (integer i = 1; i <= 2; ++i) {
	       for (integer j = 1; j <= 2; ++j) {
		    const doublereal R1TR2ij = R1.GetCol(i).Dot(R2.GetCol(j));
		    zCurr(i) += R1TR2ij * oState2.zCurr(j);
		    zPCurr(i) += R1TR2ij * oState2.zPCurr(j);
	       }
	  }
     }

     std::ostream& Print(std::ostream& os, const Mat3x3& R, const Vec3& o) const {
	  os << "o=" << o << std::endl;
	  os << "z=" << R.GetCol(1) * zCurr(1) + R.GetCol(2) * zCurr(2) << std::endl;
	  os << "zP=" << R.GetCol(1) * zPCurr(1) + R.GetCol(2) * zPCurr(2) << std::endl;
	  return os;
     }

private:
     void SaveStictionState(const grad::Vector<doublereal, 2>& z,
			    const grad::Vector<doublereal, 2>& zP);

     template <grad::index_type N>
     void SaveStictionState(const grad::Vector<grad::Gradient<N>, 2>& z,
			    const grad::Vector<grad::Gradient<N>, 2>& zP);

     const LugreData* const pData;
     grad::Vector<doublereal, 2> zPrev, zCurr, zPPrev, zPCurr;
};

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
     struct TargetFace;

     struct TargetVertex {
	  Mat3x3 R = ::Zero3x3;
	  Vec3 o = ::Zero3;
	  integer iVertex = -1;
	  integer iEdgeNeighbor = -1;
	  TargetFace* pNeighbor = nullptr;
     };

     struct TargetEdge {
	  TargetEdge(integer iNode1, integer iNode2)
	       :rgNodeIdx{iNode1 < iNode2 ? iNode1 : iNode2, iNode1 < iNode2 ? iNode2 : iNode1} {
	  }

	  TargetEdge(const TargetEdge&) = default;

	  bool operator==(const TargetEdge& oEdge) const {
	       return rgNodeIdx[0] == oEdge.rgNodeIdx[0] && rgNodeIdx[1] == oEdge.rgNodeIdx[1];
	  }

	  static constexpr integer iNumNodes = 2;
	  const std::array<integer, iNumNodes> rgNodeIdx;
     };

     struct TargetEdgeHash {
	  std::size_t operator()(const TargetEdge& oEdge) const {
	       return (std::hash<integer>{}(oEdge.rgNodeIdx[0]) << 1) ^ std::hash<integer>{}(oEdge.rgNodeIdx[1]);
	  }
     };

     struct TargetFace {
	  static constexpr integer iNumVertices = 3;

	  TargetFace(const DataManager* pDM);

	  Vec3 oc;
	  doublereal r;
	  std::array<TargetVertex, iNumVertices> rgVert;
	  integer iNumNeighbors;
	  LugreData oFrictData;
     };

     struct ContactPair {
	  ContactPair(const LugreData* pFrictData,
		      std::size_t iVertex,
		      const std::array<doublereal, TargetFace::iNumVertices>& vy)
	       :oFrictState(pFrictData), iVertex(iVertex), vy(vy) {
	  }
	  ContactPair(const ContactPair&)=default;

	  LugreState oFrictState;
	  const std::size_t iVertex;
	  std::array<doublereal, TargetFace::iNumVertices>  vy;
     };

     struct ContactVertex {
	  ContactVertex(const Vec3& o1, doublereal r1)
	       :o1(o1), r1(r1) {
	  }

	  ContactVertex(const ContactVertex&) = default;

	  Vec3 o1;
	  doublereal r1;
     };

     struct ContactNode {
	  ContactNode(const StructNode* pNode,
		      std::vector<ContactVertex>&& v,
		      integer iNumFaces);

	  const StructNode* const pContNode;
	  const std::vector<ContactVertex> rgVertices;
	  std::unordered_map<TargetFace*, ContactPair> rgContCurr, rgContPrev;
     };

     void ContactSearch();

     doublereal GetContactForce(doublereal dz, grad::LocalDofMap*) const;

     template <grad::index_type N>
     grad::Gradient<N> GetContactForce(const grad::Gradient<N>& dz, grad::LocalDofMap* pDofMap) const;

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
     grad::LocalDofMap oDofMap;
     const DifferentiableScalarFunction* pCL;
     const DataManager* const pDM;
     doublereal tCurr, tPrev;
     static constexpr grad::index_type iNumDofGradient = 13;

     enum class FrictionModel {
	  None,
	  Lugre
     } eFrictionModel;
};

LugreData::LugreData(const DataManager* pDM)
     :pDM(pDM),
      beta(1.),
      vs(0.),
      gamma(1.)
{
}

void LugreData::ParseInput(const Elem* pParent, MBDynParser& HP)
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


void LugreState::AfterConvergence()
{
     zPrev = zCurr;
     zPPrev = zPCurr;
}

template <typename T>
void LugreState::GetFrictionForce(const doublereal dt,
				  const grad::Vector<T, 2>& U,
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
	  const VVec2 Mk_U = pData->Mk * Ueff;
	  const VVec2 Ms_U = pData->Ms * Ueff;
	  const VVec2 Mk2_U = pData->Mk2 * Ueff;
	  const VVec2 Ms2_U = pData->Ms2 * Ueff;
	  const T norm_Mk2_U = sqrt(Dot(Mk2_U, Mk2_U));
	  const T a0 = norm_Mk2_U / sqrt(Dot(Mk_U, Mk_U));
	  const T a1 = sqrt(Dot(Ms2_U, Ms2_U)) / sqrt(Dot(Ms_U, Ms_U));
	  const T g = a0 + (a1 - a0) * exp(-pow(sqrt(norm_Ueff) / pData->vs, pData->gamma));

	  kappa = norm_Mk2_U / g;
     }

     const VMat2x2 A = pData->invMk2_sigma0 * kappa;
     const VMat2x2 B = A * (pData->beta * dt) + CMat2x2(1., 0., 0., 1);

     const Vector<T, 2> zP = Inv(B) * VVec2(Ueff - A * CVec2(zPrev + zPPrev * ((1 - pData->beta) * dt)));
     const Vector<T, 2> z = zPrev + (zP * pData->beta + zPPrev * (1 - pData->beta)) * dt;

     SaveStictionState(z, zP);

     tau = (pData->sigma0 * z + pData->sigma1 * zP) * p;
}

void LugreState::SaveStictionState(const grad::Vector<doublereal, 2>& z,
				   const grad::Vector<doublereal, 2>& zP)
{
     zCurr = z;
     zPCurr = zP;
}

template <grad::index_type N>
void LugreState::SaveStictionState(const grad::Vector<grad::Gradient<N>, 2>&,
				   const grad::Vector<grad::Gradient<N>, 2>&)
{
     // Do nothing
}

TriangSurfContact::TargetFace::TargetFace(const DataManager* pDM)
     :oc(::Zero3),
      r(0.),
      iNumNeighbors(0),
      oFrictData(pDM)
{
}

TriangSurfContact::ContactNode::ContactNode(const StructNode* pNode,
					    std::vector<ContactVertex>&& rgVert,
					    integer iNumFaces)
     :pContNode(pNode),
      rgVertices(std::move(rgVert))
{
     rgContCurr.reserve(iNumFaces);
     rgContPrev.reserve(iNumFaces);
}

TriangSurfContact::TriangSurfContact(unsigned uLabel, const DofOwner *pDO,
				     DataManager* pDM, MBDynParser& HP)
     :Elem(uLabel, flag(0)),
      UserDefinedElem(uLabel, pDO),
      pTargetNode(nullptr),
      dSearchRadius(std::numeric_limits<doublereal>::max()),
      pCL(nullptr),
      pDM(pDM),
      eFrictionModel(FrictionModel::None)
{
     tCurr = tPrev = pDM->dGetTime();

     if (HP.IsKeyWord("help"))
     {
	  silent_cout("Module: triangular surface contact\n"
		      "This element implements unilateral contact between an arbitrary rigid body, "
		      "represented by a triangular mesh, and a set of nodes\n");

	  if (!HP.IsArg()) {
	       throw NoErr(MBDYN_EXCEPT_ARGS);
	  }
     }

     if (!HP.IsKeyWord("target" "node")) {
	  silent_cerr("triangular surface contact(" << uLabel
		      << "): keyword \"target node\" expected at line "
		      << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     pTargetNode = pDM->ReadNode<StructNode, Node::STRUCTURAL>(HP);

     if (!HP.IsKeyWord("penalty" "function")) {
	  silent_cerr("triangular surface contact(" << uLabel
		      << "): keyword \"penalty function\" expected at line "
		      << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     pCL = dynamic_cast<const DifferentiableScalarFunction*>(ParseScalarFunction(HP, pDM));

     if (!pCL) {
	  silent_cerr("triangular surface contact(" << uLabel
		      << "): scalar function is not differentiable at line "
		      << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     if (HP.IsKeyWord("search" "radius")) {
	  dSearchRadius = HP.GetReal();
     }

     if (dSearchRadius <= 0.) {
	  silent_cerr("triangular surface contact(" << uLabel
		      << "): invalid value for search radius at line "
		      << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     TargetFace oCurrFace(pDM);

     if (HP.IsKeyWord("friction" "model")) {
	  if (HP.IsKeyWord("lugre")) {
	       eFrictionModel = FrictionModel::Lugre;
	       oCurrFace.oFrictData.ParseInput(this, HP);
	  } else if (!HP.IsKeyWord("none")) {
	       silent_cerr("triangular surface(" << uLabel
			   << "): keyword \"lugre\" or \"none\" expected at line "
			   << HP.GetLineData() << std::endl);
	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }
     }

     ReferenceFrame oRefTarget(pTargetNode);

     if (!HP.IsKeyWord("number" "of" "target" "vertices")) {
	  silent_cerr("triangular surface contact(" << uLabel
		      << "): keyword \"number of target vertices\" expected at line "
		      << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     const integer iNumVertices = HP.GetInt();

     if (iNumVertices <= 0) {
	  silent_cerr("triangular surface contact(" << uLabel
		      << "): at least one vertex is required at line "
		      << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     std::vector<Vec3> rgTargetVert;

     rgTargetVert.reserve(iNumVertices);

     for (integer i = 0; i < iNumVertices; ++i) {
	  rgTargetVert.emplace_back(HP.GetPosRel(oRefTarget));
     }

     if (!HP.IsKeyWord("number" "of" "target" "faces")) {
	  silent_cerr("triangular surface contact(" << uLabel
		      << "): keyword \"number of target faces\" expected at line "
		      << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     const integer iNumFaces = HP.GetInt();

     if (iNumFaces <= 0) {
	  silent_cerr("triangular surface contact(" << uLabel
		      << "): at least one face is required at line "
		      << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     rgTargetMesh.reserve(iNumFaces);

     std::unordered_multimap<TargetEdge, TargetFace*, TargetEdgeHash> rgEdgeMap;

     if (eFrictionModel != FrictionModel::None) {
	  rgEdgeMap.reserve(iNumFaces * TargetFace::iNumVertices);
     }

     static const integer rgEdges[3][2] = {{0, 1}, {1, 2}, {2, 0}};

     for (integer i = 0; i < iNumFaces; ++i) {
	  oCurrFace.oc = ::Zero3;

	  for (integer j = 0; j < TargetFace::iNumVertices; ++j) {
	       const integer iVertex = HP.GetInt() - 1;

	       if (iVertex < 0 || iVertex >= iNumVertices) {
		    silent_cerr("triangular surface contact(" << uLabel
				<< "): vertex index out of range (1:"
				<< iNumVertices << ") at line "
				<< HP.GetLineData() << std::endl);
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       oCurrFace.rgVert[j].iVertex = iVertex;
	       oCurrFace.rgVert[j].o = rgTargetVert[iVertex];
	       oCurrFace.oc += rgTargetVert[iVertex];
	  }

	  oCurrFace.oc /= TargetFace::iNumVertices;

	  oCurrFace.r = 0.;

	  for (integer j = 0; j < TargetFace::iNumVertices; ++j) {
	       doublereal lj = (oCurrFace.rgVert[j].o - oCurrFace.oc).Norm();
	       oCurrFace.r = std::max(oCurrFace.r, lj);
	  }

	  for (integer j = 0; j < TargetFace::iNumVertices; ++j) {
	       static const integer e1_idx[3][2] = {{1, 0}, {2, 1}, {0, 2}};
	       static const integer e2_idx[3][2] = {{2, 0}, {0, 1}, {1, 2}};
	       Vec3 e1 = oCurrFace.rgVert[e1_idx[j][0]].o - oCurrFace.rgVert[e1_idx[j][1]].o;
	       Vec3 e2 = oCurrFace.rgVert[e2_idx[j][0]].o - oCurrFace.rgVert[e2_idx[j][1]].o;
	       Vec3 e3 = e1.Cross(e2);
	       e2 = e3.Cross(e1);

	       doublereal e1n = e1.Norm(), e2n = e2.Norm(), e3n = e3.Norm();

	       if (!(e1n && e2n && e3n)) {
		    silent_cerr("triangular surface contact(" << uLabel
				<< "): singular face geometry detected at line "
				<< HP.GetLineData() << std::endl);
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       oCurrFace.rgVert[j].R = Mat3x3(e1 / e1n, e2 / e2n, e3 / e3n);
	  }

	  rgTargetMesh.emplace_back(oCurrFace);

	  if (eFrictionModel != FrictionModel::None) {
	       TargetFace* pFace = &rgTargetMesh.back();

	       for (integer j = 0; j < TargetFace::iNumVertices; ++j) {
		    TargetEdge oEdge{pFace->rgVert[rgEdges[j][0]].iVertex,
				     pFace->rgVert[rgEdges[j][1]].iVertex};
		    rgEdgeMap.emplace(oEdge, pFace);
	       }
	  }
     }

     if (eFrictionModel != FrictionModel::None) {
	  for (auto& rCurr: rgEdgeMap) {
	       auto oNeighbors = rgEdgeMap.equal_range(rCurr.first);

	       for (auto itNeighbor = oNeighbors.first; itNeighbor != oNeighbors.second; ++itNeighbor) {
		    if (itNeighbor->second == rCurr.second) {
			 continue;
		    }

		    if (rCurr.second->iNumNeighbors >= TargetFace::iNumVertices) {
			 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		    }

		    for (integer iEdge = 0; iEdge < TargetFace::iNumVertices; ++iEdge) {
			 TargetEdge oEdge{rCurr.second->rgVert[rgEdges[iEdge][0]].iVertex,
					  rCurr.second->rgVert[rgEdges[iEdge][1]].iVertex};

			 if (oEdge == itNeighbor->first) {
			      ASSERT(rCurr.second->rgVert[iEdge].pNeighbor == nullptr);

			      rCurr.second->rgVert[iEdge].pNeighbor = itNeighbor->second;

			      for (integer jEdge = 0; jEdge < TargetFace::iNumVertices; ++jEdge) {
				   TargetEdge oEdgeNeighbor{itNeighbor->second->rgVert[rgEdges[jEdge][0]].iVertex,
							    itNeighbor->second->rgVert[rgEdges[jEdge][1]].iVertex};

				   if (oEdgeNeighbor == oEdge) {
					rCurr.second->rgVert[iEdge].iEdgeNeighbor = jEdge;
					break;
				   }
			      }

			      rCurr.second->iNumNeighbors++;

			      ASSERT(rCurr.second->rgVert[iEdge].iEdgeNeighbor >= 0);
			      ASSERT(rCurr.second->rgVert[iEdge].iEdgeNeighbor < TargetFace::iNumVertices);

			      break;
			 }
		    }
	       }
	  }
     }

     if (!HP.IsKeyWord("number" "of" "contact" "nodes")) {
	  silent_cerr("triangular surface contact(" << uLabel
		      << "): keyword \"number of contact nodes\" expected at line "
		      << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     const integer iNumContactNodes = HP.GetInt();

     if (iNumContactNodes <= 0) {
	  silent_cerr("triangular surface contact(" << uLabel
		      << "): invalid number of contact nodes at line "
		      << HP.GetLineData() << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     rgContactMesh.reserve(iNumContactNodes);

     for (integer i = 0; i < iNumContactNodes; ++i) {
	  const StructNode* pContNode = pDM->ReadNode<StructNode, Node::STRUCTURAL>(HP);
	  const ReferenceFrame oRefFrame(pContNode);
	  const integer iNumVertices = HP.IsKeyWord("number" "of" "contact" "vertices") ? HP.GetInt() : 1;

	  std::vector<ContactVertex> rgVertices;

	  rgVertices.reserve(iNumVertices);

	  for (integer j = 1; j <= iNumVertices; ++j) {
	       const Vec3 o1 = HP.IsKeyWord("offset") ? HP.GetPosRel(oRefFrame) : Zero3;

	       const doublereal r1 = HP.IsKeyWord("radius") ? HP.GetReal() : 0.;

	       if (r1 < 0.) {
		    silent_cerr("triangular surface contact(" << uLabel
				<< "): invalid value for radius at line "
				<< HP.GetLineData() << std::endl);
		    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       }

	       rgVertices.emplace_back(o1, r1);
	  }

	  rgContactMesh.emplace_back(pContNode, std::move(rgVertices), iNumFaces);
     }
}

TriangSurfContact::~TriangSurfContact(void)
{

}

doublereal TriangSurfContact::GetContactForce(doublereal dz, grad::LocalDofMap*) const {
     return (*pCL)(dz);
}

template <grad::index_type N>
grad::Gradient<N> TriangSurfContact::GetContactForce(const grad::Gradient<N>& dz, grad::LocalDofMap* pDofMap) const
{
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
     *piNumRows = -(rgContactMesh.size() * iNumDofGradient);
     *piNumCols = iNumDofGradient;
}

void
TriangSurfContact::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = -(rgContactMesh.size() * iNumDofGradient);
     *piNumCols = 2 * iNumDofGradient;
}

void TriangSurfContact::ContactSearch()
{
     const Vec3& X2 = pTargetNode->GetXCurr();
     const Mat3x3& R2 = pTargetNode->GetRCurr();
     const Vec3& X2P = pTargetNode->GetVCurr();
     const Vec3& omega2 = pTargetNode->GetWCurr();

     for (auto& rNode: rgContactMesh) {
	  rNode.rgContCurr.clear();
     }

     for (auto& rFace: rgTargetMesh) {
	  for (auto& rNode: rgContactMesh) {
	       const Vec3& X1 = rNode.pContNode->GetXCurr();
	       const Mat3x3& R1 = rNode.pContNode->GetRCurr();
	       const Vec3& X1P = rNode.pContNode->GetVCurr();
	       const Vec3& omega1 = rNode.pContNode->GetWCurr();

	       for (std::size_t iVertex = 0; iVertex < rNode.rgVertices.size(); ++iVertex) {
		    const Vec3& o1 = rNode.rgVertices[iVertex].o1;
		    const doublereal r1 = rNode.rgVertices[iVertex].r1;
		    const Vec3 R1_o1 = R1 * o1;
		    const Vec3 l1 = X1 + R1_o1 - X2;
		    const doublereal dDist = (l1 - R2 * rFace.oc).Norm() - r1 - rFace.r;

		    if (dDist > dSearchRadius) {
			 continue;
		    }

		    const Vec3 dX = R2.MulTV(l1);

		    bool bInsert = true;

		    for (integer i = 0; i < TargetFace::iNumVertices; ++i) {
			 const Vec3& o2i = rFace.rgVert[i].o;
			 const Mat3x3& R2i = rFace.rgVert[i].R;
			 const doublereal dy = R2i.GetCol(2).Dot(dX - o2i);

			 if (dy < 0.) {
			      bInsert = false;
			      break;
			 }
		    }

		    if (bInsert) {
			 std::array<doublereal, TargetFace::iNumVertices> vy = {0};

			 if (eFrictionModel != FrictionModel::None) {
			      const Vec3 l1P = X1P + omega1.Cross(R1_o1) - X2P;
			      const Vec3 dXP = R2.MulTV(l1P - omega2.Cross(l1));

			      for (integer i = 0; i < TargetFace::iNumVertices; ++i) {
				   const Mat3x3& R2i = rFace.rgVert[i].R;
				   vy[i] = R2i.GetCol(2).Dot(dXP);
			      }
			 }

			 rNode.rgContCurr.emplace(&rFace, ContactPair{&rFace.oFrictData, iVertex, vy});
		    }
	       }
	  }
     }

     if (eFrictionModel != FrictionModel::None) {
	  for (auto& rNode: rgContactMesh) {
	       for (auto itCurr = rNode.rgContCurr.begin(); itCurr != rNode.rgContCurr.end(); ++itCurr) {
		    // Search for current contact face in previous set of contact faces
		    auto itPrev = rNode.rgContPrev.find(itCurr->first);

		    // Did the current contact face already exist in the previous contact set?
		    if (itPrev != rNode.rgContPrev.end()) {
			 // Reuse previous stiction states
			 itCurr->second.oFrictState = itPrev->second.oFrictState;
		    } else {
			 // Check if we can project the stiction states from neighbor faces
			 doublereal vymaxCurr = -std::numeric_limits<doublereal>::max();
			 integer iEdgeCurr = -1;

			 for (integer iEdge = 0; iEdge < TargetFace::iNumVertices; ++iEdge) {
			      const doublereal vyCurr = itCurr->second.vy[iEdge];

			      if (vyCurr >= vymaxCurr) {
				   vymaxCurr = vyCurr;
				   iEdgeCurr = iEdge;
			      }
			 }

			 if (iEdgeCurr < 0) {
			      continue;
			 }

			 itPrev = rNode.rgContPrev.find(itCurr->first->rgVert[iEdgeCurr].pNeighbor);

			 if (itPrev == rNode.rgContPrev.end()) {
			      continue;
			 }

			 // Project previous stiction state from the neighbor
			 const Mat3x3& R2iPrev = itPrev->first->rgVert[0].R;
			 const Mat3x3& R2iCurr = itCurr->first->rgVert[0].R;

			 itCurr->second.oFrictState.Project(R2iCurr, R2iPrev, itPrev->second.oFrictState);
		    }
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

     tCurr = pDM->dGetTime();

     const doublereal dt = tCurr - tPrev;
     const integer iFirstIndex2 = pTargetNode->iGetFirstMomentumIndex();

     Vector<T, 3> X1, X2, F1, M1, F2, M2;
     Vector<T, 3> XP1, XP2, omega1, omega2;
     Vector<T, 2> U, tau;
     Matrix<T, 3, 3> R1, R2;

     for (auto& rNode: rgContactMesh) {
	  F1 = M1 = F2 = M2 = ::Zero3;

	  oDofMap.Reset();

	  rNode.pContNode->GetXCurr(X1, dCoef, func, &oDofMap);
	  rNode.pContNode->GetRCurr(R1, dCoef, func, &oDofMap);

	  pTargetNode->GetXCurr(X2, dCoef, func, &oDofMap);
	  pTargetNode->GetRCurr(R2, dCoef, func, &oDofMap);

	  if (eFrictionModel != FrictionModel::None) {
	       rNode.pContNode->GetVCurr(XP1, dCoef, func, &oDofMap);
	       rNode.pContNode->GetWCurr(omega1, dCoef, func, &oDofMap);

	       pTargetNode->GetVCurr(XP2, dCoef, func, &oDofMap);
	       pTargetNode->GetWCurr(omega2, dCoef, func, &oDofMap);
	  }

	  for (auto& oContPair: rNode.rgContCurr) {
	       const TargetFace& oTargetFace = *oContPair.first;
	       const ContactVertex& oVertex = rNode.rgVertices[oContPair.second.iVertex];
	       const Vec3& o1 = oVertex.o1;
	       const doublereal r1 = oVertex.r1;
	       const Mat3x3& R2i = oTargetFace.rgVert[0].R;
	       const Vec3& o2i = oTargetFace.rgVert[0].o;

	       const Vector<T, 3> n3 = R2 * R2i.GetCol(3);
	       const Vector<T, 3> l1 = R1 * o1;
	       const Vector<T, 3> l2 = X1 + l1 - X2;
	       const T dz = Dot(n3, Vector<T, 3>(l2 - R2 * o2i));
	       const Vector<T, 3> l2c = l2 - n3 * dz;
	       const Vector<T, 3> l1c = l1 - n3 * dz;
	       const T pz = dz - r1;
	       const T F2in = GetContactForce(pz, &oDofMap);
	       Vector<T, 3> F2i = n3 * F2in;

	       if (eFrictionModel != FrictionModel::None) {
		    LugreState& oFrictState = oContPair.second.oFrictState;

		    const Vector<T, 3> dV = XP1 + Cross(omega1, l1c) - XP2 - Cross(omega2, l2c);
		    const Vector<T, 3> dV_R2 = Transpose(R2) * dV;

		    for (index_type i = 1; i <= 2; ++i) {
			 U(i) = Dot(Direct(R2i.GetCol(i)), dV_R2);
		    }

		    oFrictState.GetFrictionForce(dt, U, T(fabs(F2in)), tau);

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

     GradientAssVec<Gradient<2 * iNumDofGradient> >::InitialAssJac(this, WorkMat.SetSparse(), XCurr, INITIAL_ASS_JAC, &oDofMap);

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
	  tPrev = tCurr;

	  for (auto& oContNode: rgContactMesh) {
	       for (auto& oContPair: oContNode.rgContCurr) {
		    oContPair.second.oFrictState.AfterConvergence();
	       }

	       oContNode.rgContPrev = std::move(oContNode.rgContCurr);
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
