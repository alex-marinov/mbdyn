/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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
 * Portions Copyright Michele Attolico <attolico@aero.polimi.it>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "streamoutelem.h"
#include "geomdata.h"
#include "socketstreammotionelem.h"

/* StreamOutElem - begin */

StreamOutElem::StreamOutElem(unsigned int uL,
	const std::string& name,
	unsigned int oe)
: Elem(uL, flag(0)),
name(name),
OutputEvery(oe), OutputCounter(0)
{
	ASSERT(OutputEvery > 0);
}

StreamOutElem::~StreamOutElem(void)
{
	NO_OP;
}

Elem::Type
StreamOutElem::GetElemType(void) const
{
	return Elem::SOCKETSTREAM_OUTPUT;
}

void
StreamOutElem::WorkSpaceDim(integer* piRows, integer* piCols) const
{
	*piRows = 0;
	*piCols = 0;
}

SubVectorHandler&
StreamOutElem::AssRes(SubVectorHandler& WorkVec, doublereal dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkVec.Resize(0);
	return WorkVec;
}

VariableSubMatrixHandler& 
StreamOutElem::AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef,
		const VectorHandler& X, const VectorHandler& XP)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

/* StreamOutElem - end */


/* StreamContent - begin */

StreamContent::StreamContent(void)
{
	NO_OP;
}

StreamContent::~StreamContent(void)
{
	NO_OP;
}

void *
StreamContent::GetBuf(void) const
{
	return (void *)&buf[0];
}

int
StreamContent::GetSize(void) const
{
	return buf.size();
}

/* StreamContent - end */


/* StreamContentValue - begin */

StreamContentValue::StreamContentValue(const std::vector<ScalarValue *>& v)
: Values(v)
{
	ASSERT(Values.size() > 0);

	buf.resize(sizeof(doublereal)*Values.size());
}

StreamContentValue::~StreamContentValue(void)
{
	for (std::vector<ScalarValue *>::iterator i = Values.begin();
		i != Values.end(); ++i)
	{
		delete *i;
	}
}

void
StreamContentValue::Prepare(void)
{
	char *curbuf = &buf[0];
	for (std::vector<ScalarValue *>::iterator i = Values.begin(); i != Values.end(); ++i) {
		/* assign value somewhere into mailbox buffer */
		doublereal v = (*i)->dGetValue();

		doublereal *dbuf = (doublereal *)curbuf;
		dbuf[0] = v;

		curbuf += sizeof(doublereal);
	}
}

unsigned
StreamContentValue::GetNumChannels(void) const
{
	return Values.size();
}

/* StreamContentValue - end */

/*
 * NOTE: the use of type to determine what type of stream contents
 * to read is only a hack; needs improvements
 */

StreamContent*
ReadStreamContent(DataManager *pDM, MBDynParser& HP, StreamContent::Type type)
{
	switch (type) {
	case StreamContent::UNKNOWN:
		type = StreamContent::VALUES;
		// fallthru

	case StreamContent::VALUES:
		if (HP.IsKeyWord("motion")) {
			type = StreamContent::MOTION;

		} else if (!HP.IsKeyWord("values")) {
			silent_cerr("ReadStreamContent: "
				"missing type, assuming \"values\" "
				"at line " << HP.GetLineData() << std::endl);
		}
		break;

	default:
		break;
	}

	StreamContent *pSC = 0;
	switch (type) {
	case StreamContent::VALUES: {
		int nch = HP.GetInt();
		if (nch <= 0) {
			silent_cerr("illegal number of channels for StreamContent "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		std::vector<ScalarValue *> Values(nch);
		ReadScalarValues(pDM, HP, Values);

		SAFENEWWITHCONSTRUCTOR(pSC, StreamContentValue, StreamContentValue(Values));
		} break;

	case StreamContent::MOTION: {
		unsigned uFlags = GeometryData::X;
		if (HP.IsKeyWord("output" "flags")) {
			uFlags = 0;
			while (true) {
				if (HP.IsKeyWord("position")) {
					if (uFlags & GeometryData::X) {
						silent_cerr("StreamContent: "
							"position flag already defined "
							"at line " << HP.GetLineData()
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					uFlags |= GeometryData::X;
	
				} else if (HP.IsKeyWord("orientation" "matrix")) {
					if (uFlags & GeometryData::ORIENTATION_MASK) {
						silent_cerr("StreamContent: "
							"orientation flag already defined "
							"at line " << HP.GetLineData()
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					uFlags |= GeometryData::R;
	
				} else if (HP.IsKeyWord("orientation" "matrix" "transpose")) {
					if (uFlags & GeometryData::ORIENTATION_MASK) {
						silent_cerr("StreamContent: "
							"orientation flag already defined "
							"at line " << HP.GetLineData()
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					uFlags |= GeometryData::RT;
	
				} else if (HP.IsKeyWord("velocity")) {
					if (uFlags & GeometryData::V) {
						silent_cerr("StreamContent: "
							"velocity flag already defined "
							"at line " << HP.GetLineData()
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					uFlags |= GeometryData::V;
	
				} else if (HP.IsKeyWord("angular" "velocity")) {
					if (uFlags & GeometryData::W) {
						silent_cerr("StreamContent: "
							"angular velocity flag already defined "
							"at line " << HP.GetLineData()
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					uFlags |= GeometryData::W;
	
				} else {
					break;
				}
			}
		}
	
		std::vector<const StructNode *> nodes;
		if (HP.IsKeyWord("all")) {
			/* FIXME: todo */
	
		} else {
			while (HP.IsArg()) {
				nodes.insert(nodes.end(), pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP));
			}
		}
	
		SAFENEWWITHCONSTRUCTOR(pSC, StreamContentMotion, StreamContentMotion(uFlags, nodes));
		} break;

	default:
		silent_cerr("ReadStreamContent: unknown type "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pSC;
}
