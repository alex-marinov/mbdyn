/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#ifndef STREAMOUTELEM_H
#define STREAMOUTELEM_H

#include "elem.h"
#include "scalarvalue.h"
#include "bufmod.h"

#define UNIX_PATH_MAX	108
#define DEFAULT_PORT	9011 /* intentionally unassigned by IANA */
#define DEFAULT_HOST	"127.0.0.1"

/* StreamOutElem - begin */

class StreamOutElem : virtual public Elem {
public:
	enum Type {
		UNDEFINED = -1,

		RTAI,
		SOCKETSTREAM,
		BUFFERSTREAM
	};

protected:
	std::string name;

	/* output with lower ratio */
	unsigned int OutputEvery;
	mutable unsigned int OutputCounter;

public:
   	StreamOutElem(unsigned int uL, const std::string& name, unsigned int oe);
			
   	virtual ~StreamOutElem(void);

	virtual Elem::Type GetElemType(void) const;
	virtual void WorkSpaceDim(integer* piRows, integer* piCols) const;
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec, doublereal dCoef,
			const VectorHandler& X, const VectorHandler& XP);
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef,
			const VectorHandler& X, const VectorHandler& XP);
};

/* StreamOutElem - end */

/* StreamContent - begin */

class StreamContent {
public:
	enum Type {
		UNKNOWN = -1,

		VALUES = 0,
		MOTION = 1,

		LASTTYPE
	};

	class Modifier {
	public:
		virtual ~Modifier(void);

		virtual void Set(size_t size, const char *buf) = 0;
		virtual void Modify(void) = 0;

		virtual const void *GetOutBuf(void) const = 0;
		virtual int GetOutSize(void) const = 0;
	};

	class Copy : public Modifier {
	protected:
		size_t m_size;
		const char *m_outbuf;

	public:
		Copy(size_t size, const char *buf);

		void Set(size_t size, const char *buf);
		void Modify(void);

		virtual const void *GetOutBuf(void) const;
		virtual int GetOutSize(void) const;
	};

protected:
	/* Stream buffer */
	std::vector<char> buf;

	Modifier *m_pMod;

public:
	StreamContent(size_t size, Modifier *pMod);
	virtual ~StreamContent(void);

	void *GetBuf(void) const;
	int GetSize(void) const;

	const void *GetOutBuf(void) const;
	int GetOutSize(void) const;

	virtual void Prepare(void) = 0;
	virtual unsigned GetNumChannels(void) const = 0;

};

extern StreamContent*
ReadStreamContent(DataManager *pDM, MBDynParser& HP, StreamContent::Type type);

StreamContent::Modifier *
ReadStreamContentModifier(MBDynParser& HP, integer nCh);

/* StreamContent - end */

/* StreamContentValue - begin */

class StreamContentValue : public StreamContent {
protected:
	std::vector<ScalarValue *> Values;

public:
	StreamContentValue(const std::vector<ScalarValue *>& v,
		StreamContent::Modifier *pMod);
	virtual ~StreamContentValue(void);

	void Prepare(void);
	unsigned GetNumChannels(void) const;
};

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


/* StreamContentValue - end */

/* StreamOutEcho - begin */

class StreamOutEcho {
private:
	std::string sOutFileName;
	std::ofstream outFile;
	int iPrecision;
	doublereal dShift;

public:
	StreamOutEcho(std::string& sOutFileName, int iPrecision, doublereal dShift);
	~StreamOutEcho(void);
	bool Init(const std::string& msg, unsigned uLabel, unsigned nChannels);
	void Echo(const doublereal *pbuf, unsigned nChannels);
};

extern StreamOutEcho *
ReadStreamOutEcho(MBDynParser& HP);

/* StreamOutEcho - end */

extern Elem *
ReadOutputElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel,
	StreamOutElem::Type eType, StreamContent::Type sType);

/*----------------------------------------------------------------------------
management of 'content type' for stream output element ('values','motion', etc)
------------------------------------------------------------------------------

Rearranged by Luca Conti (May 2017) on the basis of previous existing code
(fully working, just rearranged).

Edited in order to apply the same mechanism with 'readers' and 'maps' (std::map)
  already in use for constitutive laws and drives
*/

/* stream output content type reader: every content type must inherit
from this struct and implement its own Read method */
struct StreamOutputContentTypeReader {
	virtual StreamContent* Read(DataManager* pDM, MBDynParser& HP) = 0;
};

/* default content type options: reader for 'motion' */
struct MotionContentTypeReader : public StreamOutputContentTypeReader {
	virtual StreamContent* Read(DataManager* pDM, MBDynParser& HP);
};

/* default content type options: reader for 'values' */
struct ValuesContentTypeReader : public StreamOutputContentTypeReader {
	virtual StreamContent* Read(DataManager* pDM, MBDynParser& HP);
};

/* bag of content-type readers - every content type is registered inside
of it by using SetStreamOutputContentType(...) */
typedef std::map<std::string,StreamOutputContentTypeReader*> StreamOutputContentTypeMap;
extern StreamOutputContentTypeMap streamOutputContentTypeMap;

/* stream output content type parsing checker: allows the parser
to understand if the next keyword is a content type */
struct StreamOutputContentTypeWordSetType : public HighParser::WordSet {
	virtual bool IsWord(const std::string& s) const;
};

extern StreamOutputContentTypeWordSetType streamOutputContentTypeWordSet;

/* registration function: call it to register a new content type -
it is also used by InitStreamOutputContentTypes(void) in order to load default
options ('values','motion') */
bool SetStreamOutputContentType(const char *name, StreamOutputContentTypeReader *rf);

/* initialization counter: counts how many times InitStreamOutputContentTypes
is called */
extern unsigned streamOutputInitFunctionCalls;

/* initialization function: called by mbpar.cc in order to load default
options ('values','motion') */
void InitStreamOutputContentTypes(void);

/* deallocation of all content types in streamOutputContentTypeMap - called by mbpar.cc */
void DestroyStreamOutputContentTypes(void);

#endif /* STREAMOUTELEM_H */
