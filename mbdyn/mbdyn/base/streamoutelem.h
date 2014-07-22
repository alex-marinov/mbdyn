/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#define UNIX_PATH_MAX	108
#define DEFAULT_PORT	9011 /* intentionally unassigned by IANA */
#define DEFAULT_HOST	"127.0.0.1"

/* StreamOutElem - begin */

class StreamOutElem : virtual public Elem {
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

protected:
	/* Stream buffer */
	std::vector<char> buf;

public:
	StreamContent(void);
	virtual ~StreamContent(void);

	void *GetBuf(void) const;
	int GetSize(void) const;

	virtual void Prepare(void) = 0;
	virtual unsigned GetNumChannels(void) const = 0;
};

extern StreamContent*
ReadStreamContent(DataManager *pDM, MBDynParser& HP, StreamContent::Type type);

/* StreamContent - end */

/* StreamContentValue - begin */

class StreamContentValue : public StreamContent {
protected:
	std::vector<ScalarValue *> Values;

public:
	StreamContentValue(const std::vector<ScalarValue *>& v);
	virtual ~StreamContentValue(void);

	void Prepare(void);
	unsigned GetNumChannels(void) const;
};

/* StreamContentValue - end */

#endif /* STREAMOUTELEM_H */

