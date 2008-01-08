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

/* StructOutput - begin */

StructOutput::StructOutput(const Elem *pE)
: Elem(pE->GetLabel(), 1),
NestedElem(pE)
{

}

StructOutput::~StructOutput(void)
{

}

/* StructOutput - end */

/* StructOutputCollect - begin */

void
StructOutputCollect::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP,
	GeometryData& data)
{
}

StructOutputCollect::StructOutputCollect(const Elem *pE)
: Elem(pE->GetLabel(), 1),
StructOutput(pE)
{
}

StructOutputCollect::~StructOutputCollect(void)
{
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
}

StructOutputInterp::StructOutputInterp(const Elem *pE)
: Elem(pE->GetLabel(), 1),
StructOutput(pE)
{
}

StructOutputInterp::~StructOutputInterp(void)
{
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
}

StructOutputWrite::StructOutputWrite(const Elem *pE)
: Elem(pE->GetLabel(), 1),
StructOutput(pE)
{
}

StructOutputWrite::~StructOutputWrite(void)
{
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
}

StructOutputWriteNASTRAN::StructOutputWriteNASTRAN(const Elem *pE)
: Elem(pE->GetLabel(), 1),
StructOutput(pE)
{
}

StructOutputWriteNASTRAN::~StructOutputWriteNASTRAN(void)
{
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

