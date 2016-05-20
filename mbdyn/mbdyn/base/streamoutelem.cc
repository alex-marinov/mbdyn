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

/*
 * Portions Copyright Michele Attolico <attolico@aero.polimi.it>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "streamoutelem.h"
#include "geomdata.h"
#include "socketstream_out_elem.h"
#include "bufferstream_out_elem.h"
#include "socketstreammotionelem.h"
#include "bufmod.h"

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

StreamContent::StreamContent(size_t size, StreamContent::Modifier *pMod)
: buf(size), m_pMod(pMod)
{
	if (m_pMod == 0) {
		m_pMod = new StreamContent::Copy(size, &buf[0]);

	} else {
		m_pMod->Set(size, &buf[0]);
	}
}

StreamContent::~StreamContent(void)
{
	if (m_pMod != 0) {
		delete m_pMod;
		m_pMod = 0;
	}
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

const void *
StreamContent::GetOutBuf(void) const
{
	return m_pMod->GetOutBuf();
}

int
StreamContent::GetOutSize(void) const
{
	return m_pMod->GetOutSize();
}

StreamContent::Modifier::~Modifier(void)
{
	NO_OP;
}

StreamContent::Copy::Copy(size_t size, const char *outbuf)
: m_size(size), m_outbuf(outbuf)
{
	NO_OP;
}

void
StreamContent::Copy::Set(size_t size, const char *outbuf)
{
	m_size = size;
	m_outbuf = outbuf;
}

void
StreamContent::Copy::Modify(void)
{
	NO_OP;
}

const void *
StreamContent::Copy::GetOutBuf(void) const
{
	return m_outbuf;
}

int
StreamContent::Copy::GetOutSize(void) const
{
	return m_size;
}

class StreamContentCopyCast : public StreamContent::Modifier {
protected:
	size_t m_size;
	const char *m_buf;
	std::vector<char> m_outbuf;
	std::vector<BufCast *> m_data;

public:
	StreamContentCopyCast(size_t size, const char *buf, size_t outsize, const std::vector<BufCast *>& data);

	void Set(size_t size, const char *buf);
	void Modify(void);

	const void *GetOutBuf(void) const;
	int GetOutSize(void) const;
};

StreamContentCopyCast::StreamContentCopyCast(size_t size, const char *buf, size_t outsize, const std::vector<BufCast *>& data)
: m_size(size), m_buf(buf), m_outbuf(outsize), m_data(data)
{
#ifdef DEBUG
	std::vector<BufCast *>::const_iterator i = data.end();
	ASSERT(i > data.begin());
	--i;
	size_t minsize = (*i)->offset() + (*i)->size();
	ASSERT(outsize >= minsize);
#endif
}

void
StreamContentCopyCast::Set(size_t size, const char *buf)
{
	m_size = size;
	m_buf = buf;

	// FIXME: what about outbuf?
	// in principle, size is sizeof(doublereal)*nChannels
	ASSERT(m_size == sizeof(doublereal)*m_data.size());
}

void
StreamContentCopyCast::Modify(void)
{
	ASSERT(m_size == m_data.size()*sizeof(doublereal));

	doublereal *rbuf = (doublereal *)&m_buf[0];
	for (size_t i = 0; i < m_data.size(); i++) {
		m_data[i]->uncast(&m_outbuf[0], rbuf[i]);
	}
}

const void *
StreamContentCopyCast::GetOutBuf(void) const
{
	return &m_outbuf[0];
}

int
StreamContentCopyCast::GetOutSize(void) const
{
	return m_outbuf.size();
}

/* StreamContent - end */


/* StreamContentValue - begin */

StreamContentValue::StreamContentValue(const std::vector<ScalarValue *>& v,
	StreamContent::Modifier *pMod)
: StreamContent(sizeof(doublereal)*v.size(), pMod), Values(v)
{
	ASSERT(Values.size() > 0);
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

	m_pMod->Modify();
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

StreamContent::Modifier *
ReadStreamContentModifier(MBDynParser& HP, integer nCh)
{
	StreamContent::Modifier *pSCM(0);

	if (HP.IsKeyWord("copy" "cast")) {
		std::vector<BufCast *> data(nCh);
		ReadBufCast(HP, data);
		size_t minsize = data[data.size() - 1]->offset() + data[data.size() - 1]->size();
		size_t size = minsize;
		if (HP.IsKeyWord("size")) {
			integer i = HP.GetInt();
			if (i <= 0) {
				silent_cerr("ReadStreamContentModifier: invalid size " << i
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			size = size_t(i);
			if (size < minsize) {
				silent_cerr("ReadStreamContentModifier: size " << size
					<< " is less than min size " << minsize
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		pSCM = new StreamContentCopyCast(0, 0, size, data);

	} else if (!HP.IsKeyWord("copy")) {
		silent_cerr("ReadStreamContentModifier: unknown modifier type at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pSCM;
}

StreamContent *
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

		StreamContent::Modifier *pMod(0);
		if (HP.IsKeyWord("modifier")) {
			pMod = ReadStreamContentModifier(HP, nch);
		}

		SAFENEWWITHCONSTRUCTOR(pSC, StreamContentValue, StreamContentValue(Values, pMod));
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

		// FIXME: right now, we don't modify stream motion stuff	
		StreamContent::Modifier *pMod(0);

		SAFENEWWITHCONSTRUCTOR(pSC, StreamContentMotion, StreamContentMotion(uFlags, nodes, pMod));
		} break;

	default:
		silent_cerr("ReadStreamContent: unknown type "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pSC;
}

/* StreamOutEcho - begin */

StreamOutEcho::StreamOutEcho(std::string& sOutFileName, int iPrecision, doublereal dShift)
: sOutFileName(sOutFileName),
iPrecision(iPrecision),
dShift(dShift)
{
	NO_OP;
}

StreamOutEcho::~StreamOutEcho(void)
{
	NO_OP;
}

bool
StreamOutEcho::Init(const std::string& msg, unsigned uLabel, unsigned nChannels)
{
	outFile.open(sOutFileName.c_str());
	if (!outFile) {
		silent_cerr(msg << "(" << uLabel << "): "
			"unable to open echo file \"" << sOutFileName << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (iPrecision > 0) {
		outFile.precision(iPrecision);
	}
	outFile.setf(std::ios::scientific);

	outFile
		<< "# generated by " << msg << "(" << uLabel << ")"
		<< std::endl;
	if (nChannels == 1) {
		outFile
			<< "# Channel #1"
			<< std::endl;

	} else {
		outFile
			<< "# Channels #1-" << nChannels
			<< std::endl;
	}

	return true;
}

void
StreamOutEcho::Echo(const doublereal *pbuf, unsigned size)
{
	outFile << pbuf[0];
	for (unsigned i = 1; i < size; i++) {
		outFile << " " << pbuf[i];
	}
	outFile << std::endl;
}

StreamOutEcho *
ReadStreamOutEcho(MBDynParser& HP)
{
	StreamOutEcho *pSOE(0);

	std::string sOutFileName;
	int iPrecision = -1;
	doublereal dShift = 0.;

	if (HP.IsKeyWord("echo")) {
		const char *s = HP.GetFileName();
		if (s == NULL) {
			silent_cerr("ReadStreamOutEcho: "
				"unable to parse echo file name "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		sOutFileName = s;

		if (HP.IsKeyWord("precision")) {
			iPrecision = HP.GetInt();
			if (iPrecision <= 0) {
				silent_cerr("ReadStreamOutEcho: "
					"invalid echo precision " << iPrecision
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		if (HP.IsKeyWord("shift")) {
			dShift = HP.GetReal();
		}

		pSOE = new StreamOutEcho(sOutFileName, iPrecision, dShift);
	}

	return pSOE;
}

/* StreamOutEcho - end */

Elem *
ReadOutputElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel, StreamOutElem::Type eType, StreamContent::Type sType)
{
	Elem *pE(0);

	if (eType == StreamOutElem::UNDEFINED) {
		sType = StreamContent::UNKNOWN;

		if (HP.IsKeyWord("socket" "stream")) {
			eType = StreamOutElem::SOCKETSTREAM;

		} else if (HP.IsKeyWord("buffer" "stream")) {
			eType = StreamOutElem::BUFFERSTREAM;

		} else {
			silent_cerr("ReadOutputElem(" << uLabel << "): "
				"unknown \"type\" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	switch (eType) {
	case StreamOutElem::SOCKETSTREAM:
		pE = ReadSocketStreamElem(pDM, HP, uLabel, sType);
		break;

	case StreamOutElem::BUFFERSTREAM:
		pE = ReadBufferStreamElem(pDM, HP, uLabel, sType);
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pE;
}
