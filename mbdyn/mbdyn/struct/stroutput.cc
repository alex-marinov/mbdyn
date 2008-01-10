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

/* StructOutputEnd - end */

/* StructOutputStart - begin */

void
StructOutputStart::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP,
	GeometryData& data)
{
	dynamic_cast<StructOutputManip *>(pElem)->AfterConvergence(X, XP, data);
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
StructOutput::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP,
	GeometryData& data)
{
	dynamic_cast<StructOutputManip *>(pElem)->AfterConvergence(X, XP, data);
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

void
StructOutputCollect::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP,
	GeometryData& data)
{
#if 0
	if (pRefNode) {
		// TODO
		//
		//	~p = RRef^T * (x - xRef)
		//	~R = RRef^T * R
		//	~v = RRef^T * (v - vRef)	(RRef^T * v ?)
		//	~w = RRef^T * (w - wRef)	(RRef^T * w ?)

	} else {
		data.data.resize(Nodes.size());
		
		for (unsigned i = 0; i < Nodes.size(); i++) {
			if (data.uFlags & GeometryData::X) {
				data.data[i].X = Nodes[i]->GetXCurr();
			}
			if (data.uFlags & GeometryData::R) {
				data.data[i].R = Nodes[i]->GetRCurr();
			}
			if (data.uFlags & GeometryData::V) {
				data.data[i].V = Nodes[i]->GetVCurr();
			}
			if (data.uFlags & GeometryData::W) {
				data.data[i].W = Nodes[i]->GetWCurr();
			}

			if (data.uFlags & GeometryData::XPP) {
				data.data[i].XPP = Nodes[i]->GetXPPCurr();
			}
			if (data.uFlags & GeometryData::WP) {
				data.data[i].WP = Nodes[i]->GetWPCurr();
			}
		}
	}
#endif

	StructOutputStart::AfterConvergence(X, XP, data);
}

StructOutputCollect::StructOutputCollect(const Elem *pE)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutputStart(pE)
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
StructOutputCollect::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	// Finto!
	GeometryData data;
	AfterConvergence(X, XP, data);
}

static Elem *
ReadStructOutputCollect(DataManager *pDM, MBDynParser& HP, unsigned int uLabel)
{
	return 0;
}

/* StructOutputCollect - end */

/* StructOutputInterp - begin */

void
StructOutputInterp::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP,
	GeometryData& data)
{
	StructOutput::AfterConvergence(X, XP, data);
}

StructOutputInterp::StructOutputInterp(const Elem *pE)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutput(pE)
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

static Elem *
ReadStructOutputInterp(DataManager *pDM, MBDynParser& HP, unsigned int uLabel)
{
	return 0;
}

/* StructOutputInterp - end */

/* StructOutputWrite - begin */

void
StructOutputWrite::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP,
	GeometryData& data)
{
	NO_OP;
}

StructOutputWrite::StructOutputWrite(const Elem *pE)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutputEnd(pE->GetLabel(), pE->fToBeOutput())
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
	return 0;
}

/* StructOutputWrite - end */

/* StructOutputWriteNASTRAN - begin */

void
StructOutputWriteNASTRAN::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP,
	GeometryData& data)
{
	NO_OP;
}

StructOutputWriteNASTRAN::StructOutputWriteNASTRAN(const Elem *pE)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
StructOutputEnd(pE->GetLabel(), pE->fToBeOutput())
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
	if (HP.IsKeyWord("collect")) {
		return ReadStructOutputCollect(pDM, HP, uLabel);

	} else if (HP.IsKeyWord("interpolate")) {
		return ReadStructOutputInterp(pDM, HP, uLabel);

	} else if (HP.IsKeyWord("write")) {
		return ReadStructOutputWrite(pDM, HP, uLabel);

	} else if (HP.IsKeyWord("write" "NASTRAN")) {
		return ReadStructOutputWriteNASTRAN(pDM, HP, uLabel);
	}

	silent_cerr("StructOutput(" << uLabel << "): "
		"unknown type at line " << HP.GetLineData() << std::endl);

	return 0;
}

