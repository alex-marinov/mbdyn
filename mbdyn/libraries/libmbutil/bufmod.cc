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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

//#include "unistd.h"
//#include "ac/f2c.h"

#include "parser.h"
#include "bufmod.h"

BufCast::BufCast(size_t offset)
: m_offset(offset)
{
	NO_OP;
}

BufCast::~BufCast(void)
{
	NO_OP;
}

template <class T>
class TBufCast : public BufCast {
public:
	TBufCast(size_t offset) : BufCast(offset) {};

	size_t size(void) const {
		return sizeof(T);
	};

	size_t offset(void) const {
		return m_offset;
	};

	doublereal cast(const void *pFrom) const {
		const char *p = &((const char *)pFrom)[m_offset];
		return doublereal(*((T *)p));
	};

	void uncast(void *pTo, doublereal d) const {
		char *p = &((char *)pTo)[m_offset];
		((T *)p)[0] = T(d);
	};

	BufCast *copy(size_t offset) const {
		return new TBufCast<T>(offset);
	};
};

// TODO: network-independent TBufCast 

static BufCast *
ReadOneBufCast(HighParser& HP, size_t& offset)
{
	BufCast *pBC(0);

	if (HP.IsKeyWord("int8_t")) {
		pBC = new TBufCast<int8_t>(offset);
		offset += sizeof(int8_t);

	} else if (HP.IsKeyWord("uint8_t")) {
		pBC = new TBufCast<uint8_t>(offset);
		offset += sizeof(uint8_t);

	} else if (HP.IsKeyWord("int16_t")) {
		pBC = new TBufCast<int16_t>(offset);
		offset += sizeof(int16_t);

	} else if (HP.IsKeyWord("uint16_t")) {
		pBC = new TBufCast<uint16_t>(offset);
		offset += sizeof(uint16_t);

	} else if (HP.IsKeyWord("int32_t")) {
		pBC = new TBufCast<int32_t>(offset);
		offset += sizeof(int32_t);

	} else if (HP.IsKeyWord("uint32_t")) {
		pBC = new TBufCast<uint32_t>(offset);
		offset += sizeof(uint32_t);

	} else if (HP.IsKeyWord("float")) {
		pBC = new TBufCast<float>(offset);
		offset += sizeof(float);

	} else if (HP.IsKeyWord("double")) {
		pBC = new TBufCast<double>(offset);
		offset += sizeof(double);

	} else if (HP.IsKeyWord("skip")) {
		integer skip = HP.GetInt();
		if (skip < 0) {
			silent_cerr("ReadBufCast: invalid number of bytes " << skip
				<< " to be skipped at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		offset += skip;

	} else {
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pBC;
}

void
ReadBufCast(HighParser& HP, std::vector<BufCast *>& data)
{
	if (HP.IsKeyWord("all")) {
		size_t size(0), offset(0);
		BufCast *pBC = ReadOneBufCast(HP, size);
		data[0] = pBC;
		for (size_t i = 1; i < data.size(); i++) {
			offset += size;
			data[i] = data[i - 1]->copy(offset);
		}

	} else {
		size_t offset(0);
		for (size_t i = 0; i < data.size(); i++) {
retry:;
			BufCast *pBC = ReadOneBufCast(HP, offset);
			if (pBC == 0) {
				// got skip
				goto retry;
			}
			data[i] = pBC;
		}
	}
}

