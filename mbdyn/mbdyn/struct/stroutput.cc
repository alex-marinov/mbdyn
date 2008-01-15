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
StructOutputStart::Manipulate(const GeometryData& data)
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

/* StructOutput - end */

/* StructOutput - begin */

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
StructOutputCollectBase::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	Manipulate_int();
}

void
StructOutputCollectBase::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	Manipulate_int();
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

	StructOutputStart::Manipulate(data);
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

	StructOutputStart::Manipulate(data);
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

	uFlags = GeometryData::X;

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
StructOutputInterpBase::InterpInit(void)
{
	// gather interpolation data
	InterpInit_int();

	// initialize interpolation data
}

void
StructOutputInterpBase::Manipulate(const GeometryData& orig_data)
{
	// do the interpolation at each call

	StructOutput::Manipulate(data);
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
StructOutputInterp::InterpInit_int(void)
{
	NO_OP;
}

StructOutputInterp::StructOutputInterp(const Elem *pE)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutputInterpBase(pE)
{
	NO_OP;
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
StructOutputInterpOP2::InterpInit_int(void)
{
	// read op2 file
	NO_OP;
}

StructOutputInterpOP2::StructOutputInterpOP2(const Elem *pE)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutputInterpBase(pE)
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

static Elem *
ReadStructOutputInterp(DataManager *pDM, MBDynParser& HP, const Elem *pNestedElem)
{
	// collect parameters and pass to constructor

	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, StructOutputInterpOP2,
		StructOutputInterpOP2(pNestedElem));

	return pEl;
}

/* StructOutputInterp - end */

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
StructOutputWrite::Manipulate(const GeometryData& data)
{
	// open file
	std::ofstream outf(outfilename.c_str());
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
	outf << "# Label";
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

